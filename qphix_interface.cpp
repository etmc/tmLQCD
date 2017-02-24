/***********************************************************************
 *
 * Copyright (C) 2015 Mario Schroeck
 *
 * This file is part of tmLQCD.
 *
 * tmLQCD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * tmLQCD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with tmLQCD.  If not, see <http://www.gnu.org/licenses/>.
 *
 ***********************************************************************/
/***********************************************************************
*
* File qphix_interface.c
*
* Author: Mario Schroeck <mario.schroeck@roma3.infn.it>
* 
* Last changes: 03/2015
*
*
* Integration of the QUDA inverter for multi-GPU usage
*
* The externally accessible functions are
*
*   void _initQphix( int verbose )
*     Initializes the QUDA library. Carries over the lattice size and the 
*     MPI process grid and thus must be called after initializing MPI (and 
*     after 'read_infile(argc,argv)').
*     Memory for the QUDA gaugefield on the host is allocated but not filled 
*     yet (the latter is done in _loadGaugeQphix(), see below).
*     Performance critical settings are done here and can be changed.
*     Input parameter: verbose (0=SILENT, 1=SUMMARIZE, 2=VERBOSE).
*
*   void _endQphix()
*     Finalizes the QUDA library. Call before MPI_Finalize().
*
*   void _loadGaugeQphix()
*     Copies and reorders the gaugefield on the host and copies it to the GPU.
*     Must be called between last changes on the gaugefield (smearing etc.)
*     and first call of the inverter. In particular, 'boundary(const double kappa)'
*     must be called before if nontrivial boundary conditions are to be used since
*     those will be applied directly to the gaugefield.
*
*   double tmcgne_qphix(int nmx,double res,int k,int l,int *status,int *ifail)
*     The same functionality as 'tmcgne' (see tmcg.c) but inversion is performed on 
*     the GPU using QUDA. Final residuum check is performed on the host (CPU)
*     with the function 'void tmQnohat_dble(int k,int l)' (see tmdirac.c).
*
*   void tmQnohat_qphix(int k, int l)
*     The implementation of the QUDA equivalent of 'tmQnohat_dble'. 
*
**************************************************************************/

#include "qphix_interface.h"

// #undef SEEK_SET
// #undef SEEK_CUR
// #undef SEEK_END

// include mpi.h first
#include <mpi.h>
#include "global.h"
extern "C" {
#include "boundary.h"
#include "linalg/convert_eo_to_lexic.h"
#include "solver/solver.h"
//#include "solver/solver_field.h"
#include "gettime.h"
#include "update_backward_gauge.h"
}

#include "timeDslashNoQDP.h"
#include <omp.h>
#include "qphix/wilson.h"
#if 1
#include "qphix/blas.h"
#include "qphix/invcg.h"
//#include "qphix/invbicgstab.h"
#include "qphix/print_utils.h"
#endif

#include <cstdlib>
#include <cstring>
#include <cfloat>

using namespace std;
using namespace QPhiX;

#ifndef QPHIX_SOALEN
#define QPHIX_SOALEN 8
#endif

#if defined(QPHIX_MIC_SOURCE)
#define VECLEN_SP 16
#define VECLEN_HP 16
#define VECLEN_DP 8
#endif


#if defined(QPHIX_AVX_SOURCE)
#define VECLEN_SP 8
#define VECLEN_DP 4
#endif

#if defined(QPHIX_SCALAR_SOURCE)
#define VECLEN_SP 1
#define VECLEN_DP 1
#endif

#if defined(QPHIX_QPX_SOURCE)
#define VECLEN_SP 4
#define VECLEN_DP 4
#endif

#ifdef QMP_COMMS
#include <qmp.h>
#endif

int By;
int Bz;
int NCores;
int Sy;
int Sz;
int PadXY;
int PadXYZ;
int MinCt;
int N_simt;
bool compress12;
QphixPrec precision;

int subLattSize[4];
int lattSize[4];

// Hardwire these for now.
int iters = 1;
int qmp_geom[4]={1,1,1,1};


template<typename T>
struct rsdTarget {
  static const double value;
};

template<>
const double rsdTarget<half>::value = (double)(1.0e-4);

template<>
const double rsdTarget<float>::value = (double)(1.0e-7);

template<>
const double rsdTarget<double>::value = (double)(1.0e-12);


// define order of the spatial indices
// default is LX-LY-LZ-T, see below def. of local lattice size, this is related to
// the gamma basis transformation from tmLQCD -> UKQCD
// for details see https://github.com/lattice/qphix/issues/157
#define USE_LZ_LY_LX_T 0

// TRIVIAL_BC are trivial (anti-)periodic boundary conditions,
// i.e. 1 or -1 on last timeslice
// tmLQCD uses twisted BC, i.e. phases on all timeslices.
// if using TRIVIAL_BC: can't compare inversion result to tmLQCD
// if not using TRIVIAL_BC: BC will be applied to gauge field,
// can't use 12 parameter reconstruction
#define TRIVIAL_BC 0

// final check of residual with DD functions on the CPU
#define FINAL_RESIDUAL_CHECK_CPU_DD 1


#define MAX(a,b) ((a)>(b)?(a):(b))


enum EvenOdd { Even = 0, Odd };
enum class spinor_type { half, full };

void print_ptr(spinor_type type, void* ptr, string text)
{
	cout << endl << text << endl;
	double* show_ptr = (double*) ptr;
  long unsigned vol = 24*lattSize[0]*lattSize[1]*lattSize[2]*lattSize[3];
  if(type == spinor_type::half) vol /= 2;

	for(long unsigned i=0; i<vol; ++i) {
		if(fabs(show_ptr[i]) > DBL_EPSILON) {
			cout << i << " : " << show_ptr[i] << endl;
		}
	}
	cout << endl;
}

void print_gauge_ptr(void* ptr, string text)
{
	cout << endl << text << endl;
	double* show_ptr = (double*) ptr;
	// for(int i=0; i<36; ++i) {
	for(int i=0; i<36*lattSize[0]*lattSize[1]*lattSize[2]*lattSize[3]; ++i) {
		if(fabs(show_ptr[i]) > DBL_EPSILON) {
			cout << i << " : " << show_ptr[i] << endl;
		}
	}
	cout << endl;
}

// function that maps coordinates in the communication grid to MPI ranks
int commsMap(const int *coords, void *fdata)
{
#if USE_LZ_LY_LX_T
  int n[4] = {coords[3], coords[2], coords[1], coords[0]};
#else
  int n[4] = {coords[3], coords[0], coords[1], coords[2]};
#endif

  int rank = 0;

#ifdef TM_USE_MPI
  MPI_Cart_rank( g_cart_grid, n, &rank );
#endif

  return rank;
}


void _initQphix(int argc, char **argv, int By_, int Bz_, int NCores_, int Sy_, int Sz_, int PadXY_, int PadXYZ_, int MinCt_, int c12, QphixPrec precision_)
{
	// Global Lattice Size
	lattSize[0] = LX*g_nproc_x;
	lattSize[1] = LY*g_nproc_y;
	lattSize[2] = LZ*g_nproc_z;
	lattSize[3] = T*g_nproc_t;

	// Local Lattice Size
	subLattSize[0] = LX;
	subLattSize[1] = LY;
	subLattSize[2] = LZ;
	subLattSize[3] = T;

	By = By_;
	Bz = Bz_;
	NCores = NCores_;
	Sy = Sy_;
	Sz = Sz_;
	PadXY = PadXY_;
	PadXYZ = PadXYZ_;
	MinCt = MinCt_;
	N_simt = Sy_*Sz_;
	compress12 = c12;
	precision = precision_;

	omp_set_num_threads(NCores*Sy*Sz);

#ifdef QPHIX_QMP_COMMS
	// Initialize QMP
	QMP_thread_level_t prv;
	if( QMP_init_msg_passing(&argc, &argv, QMP_THREAD_SINGLE, &prv) != QMP_SUCCESS ) {
		QMP_error("Failed to initialize QMP\n");
		abort();

	}
	if ( QMP_is_primary_node() ) {
		printf("QMP IS INITIALIZED\n");
	}

	// Declare the logical topology
	if ( QMP_declare_logical_topology(qmp_geom, 4)!= QMP_SUCCESS ) {
		QMP_error("Failed to declare QMP Logical Topology\n");
		abort();
	}
#endif

	// QPhiX::masterPrintf("# Values used at initialisation:\n");
	// QPhiX::masterPrintf("#  By = %d\n", By);
	// QPhiX::masterPrintf("#  Bz = %d\n", Bz);
	// QPhiX::masterPrintf("#  NCores = %d\n", NCores);
	// QPhiX::masterPrintf("#  Sy = %d\n", Sy);
	// QPhiX::masterPrintf("#  Sz = %d\n", Sz);
	// QPhiX::masterPrintf("#  PadXY = %d\n", PadXY);
	// QPhiX::masterPrintf("#  PadXYZ = %d\n", PadXYZ);
	// QPhiX::masterPrintf("#  MinCt = %d\n", MinCt);
	// QPhiX::masterPrintf("#  N_simt = %d\n", N_simt);
	// QPhiX::masterPrintf("#  compress12 = %d\n", compress12);
	// QPhiX::masterPrintf("#  precision = %d\n", precision);

	// QPhiX::masterPrintf("# Declared QMP Topology: %d %d %d %d\n\n",
	// qmp_geom[0], qmp_geom[1], qmp_geom[2], qmp_geom[3]);


#ifdef QPHIX_QPX_SOURCE
	if( thread_bind ) {
		QPhiX::setThreadAffinity(NCores_user, Sy_user*Sz_user);
	}
	QPhiX::reportAffinity();
#endif
}

// Finalize the QPhiX library
void _endQphix()
{
}

// Set QPhiX gauge field to unit gauge
void unit_gauge_QPhiX(double* qphix_gauge)
{
  double startTime = gettime();
  int Nz = 2;
  int Ns = 4;
  int Nc = 3;

  for( int x0=0; x0<T; x0++ )
    for( int x1=0; x1<LX/2; x1++ )
      for( int x2=0; x2<LY; x2++ )
        for( int x3=0; x3<LZ; x3++ )
        {
          int qphix_idx = x1 + LX/2*x2 + LX/2*LY*x3 + LX/2*LY*LZ*x0;

          for (int dim=0; dim<4; dim++) { // tmLQCD \mu
            for (int dir=0; dir<2; dir++) { // backward/forward
              for (int c1=0; c1<Nc; c1++) {
                for (int c2=0; c2<Nc; c2++) {
                  for (int z=0; z<Nz; z++) {
                    int q_mu = 2 * dim + dir;
                    int q_inner_idx = z + c2*Nz + c1*Nz*Nc + q_mu*Nz*Nc*Nc + qphix_idx*8*Nc*Nc*Nz;
                    if ( c1 == c2 && z == RE )
                    // if ( c1 == c2 ) // 1. + I
                      qphix_gauge[q_inner_idx] = 1.0;
                    else
                      qphix_gauge[q_inner_idx] = 0.0;
                  }
                }
              }
            }
          }

        } // volume

  double endTime = gettime();
  double diffTime = endTime - startTime;
  printf("time spent in unit_gauge_QPhiX: %f secs\n", diffTime);
}

// Reorder an odd tmLQCD spinor to an (odd) QPhiX spinor
void reorder_gauge_toQphix(double* qphix_gauge_cb0, double* qphix_gauge_cb1)
{
  double startTime = gettime();
  double *in;
  double *out;
  int change_dim[4] = {3, 0, 1, 2};
  int Nz = 2;
  int Ns = 4;
  int Nc = 3;
  int tm_idx = 0;

  in = (double*) &g_gauge_field[0][0].c00;

  // now copy and reorder from tempSpinor to spinor
  for( int x0=0; x0<T; x0++ )
    for( int x1=0; x1<LX; x1++ )
      for( int x2=0; x2<LY; x2++ )
        for( int x3=0; x3<LZ; x3++ )
        {
          int qphix_idx = x1/2 + LX/2*x2 + LX/2*LY*x3 + LX/2*LY*LZ*x0;

          if( (x0+x1+x2+x3)&1 ) { // cb1
            out = qphix_gauge_cb1;
          } else { // cb0
            out = qphix_gauge_cb0;
          }

          for (int dim=0; dim<4; dim++) { // tmLQCD \mu
            for (int dir=0; dir<2; dir++) { // backward/forward

              if(dir==0) { // this is the adjoint gauge field to be
                tm_idx = g_idn[ g_ipt[x0][x1][x2][x3] ][dim];
              } else {
                tm_idx = g_ipt[x0][x1][x2][x3];
              }

              for (int c1=0; c1<Nc; c1++) {
                for (int c2=0; c2<Nc; c2++) {
                  for (int z=0; z<Nz; z++) {
                    int q_mu = 2 * change_dim[dim] + dir;
                    int t_inner_idx = z + c2*Nz + c1*Nz*Nc +  dim*Nz*Nc*Nc + tm_idx*Nz*Nc*Nc*4;
                    // QPHIX gauge field transposed w.r.t. tmLQCD:
                    int q_inner_idx = z + c1*Nz + c2*Nz*Nc + q_mu*Nz*Nc*Nc + qphix_idx*8*Nc*Nc*Nz;
                    out[q_inner_idx] = in[t_inner_idx];
                  }
                }
              }

            } // dir
          } // dim

        } // volume

  double endTime = gettime();
  double diffTime = endTime - startTime;
  printf("time spent in reorder_spinor_toQphix: %f secs\n", diffTime);
}

// Reorder an odd tmLQCD spinor to an (odd) QPhiX spinor
void reorder_spinor_toQphix(double* tm_spinor, double* qphix_spinor_cb0, double* qphix_spinor_cb1)
{
  double startTime = gettime();
  double *in;
  double *out;
  int Nz = 2;
  int Ns = 4;
  int Nc = 3;

  double K1[4] = {1.0, -1.0, -1.0, 1.0};
  int change_spin[4] = {3, 2, 1, 0};

  // now copy and reorder from tempSpinor to spinor
  for( int x0=0; x0<T; x0++ )
    for( int x1=0; x1<LX; x1++ )
      for( int x2=0; x2<LY; x2++ )
        for( int x3=0; x3<LZ; x3++ )
        {
          int qphix_idx = x1/2 + LX/2*x2 + LX/2*LY*x3 + LX/2*LY*LZ*x0;
          int tm_idx    = g_ipt[x0][x1][x2][x3];

          in  = tm_spinor + Ns*Nc*Nz * tm_idx;

          if( (x0+x1+x2+x3)&1 ) { // cb1
            out = qphix_spinor_cb1 + Ns*Nc*Nz * qphix_idx;
          } else { // cb0
            out = qphix_spinor_cb0 + Ns*Nc*Nz * qphix_idx;
          }

          // gamma basis transformation
          for (int s=0; s<Ns; s++) { // QPhiX spin index
            for (int c=0; c<Nc; c++) {
              for (int z=0; z<Nz; z++) {
                int qphix_internal_idx  = z + s*Nz + c*Nz*Ns;
                int tm_internal_idx     = z + c*Nz + change_spin[s]*Nz*Nc;
                out[qphix_internal_idx] = K1[s] * in[tm_internal_idx];
              }
            }
          }

        } // volume

  double endTime = gettime();
  double diffTime = endTime - startTime;
  printf("time spent in reorder_spinor_toQphix: %f secs\n", diffTime);
}

// Reorder an (odd) QPhiX spinor to an odd tmLQCD spinor and multiply output spinor (tm) with normFac
void reorder_spinor_fromQphix(double* tm_spinor, double* qphix_spinor_cb0, double* qphix_spinor_cb1, double normFac = 1.0)
{
  double startTime = gettime();
  double *in;
  double *out;
  int Nz = 2;
  int Ns = 4;
  int Nc = 3;

  double K1[4] = {1.0, -1.0, -1.0, 1.0};
  int change_spin[4] = {3, 2, 1, 0};

  // Reorder from QPhiX to tmLQCD
  for( int x0=0; x0<T; x0++ )
    for( int x1=0; x1<LX; x1++ )
      for( int x2=0; x2<LY; x2++ )
        for( int x3=0; x3<LZ; x3++ )
        {
          int qphix_idx = x1/2 + LX/2*x2 + LX/2*LY*x3 + LX/2*LY*LZ*x0;
          int tm_idx    = g_ipt[x0][x1][x2][x3];

          if( (x0+x1+x2+x3)&1 ) { // cb1
            in = qphix_spinor_cb1 + Ns*Nc*Nz * qphix_idx;
          } else { // cb0
            in = qphix_spinor_cb0 + Ns*Nc*Nz * qphix_idx;
          }

          out = tm_spinor + Ns*Nc*Nz * tm_idx;

          // gamma basis transformation
          for (int s=0; s<Ns; s++) { // tmlQCD spin index
            for (int c=0; c<Nc; c++) {
              for (int z=0; z<Nz; z++) {
                int tm_internal_idx    = z + c*Nz + s*Nz*Nc;
                int qphix_internal_idx = z + change_spin[s]*Nz + c*Nz*Ns;
                out[tm_internal_idx]   = normFac * K1[s] * in[qphix_internal_idx];
              }
            }
          }

        } // volume

  double endTime = gettime();
  double diffTime = endTime - startTime;
  printf("time spent in reorder_spinor_fromQphix: %f secs\n", diffTime);
}


template<typename FT, int V, int S, bool compress>
void
invert(spinor * const tmlqcd_out, spinor * const tmlqcd_in, const int max_iter, double eps_sq, const int rel_prec)
{
  typedef typename Geometry<FT,V,S,compress>::SU3MatrixBlock QGauge;
  typedef typename Geometry<FT,V,S,compress>::FourSpinorBlock QSpinor;

  bool verbose = true;

  // Work out the size of checkerboarded X-dimension
  int X1h = lattSize[0]/2;
  int Nx = lattSize[0];
  int Ny = lattSize[1];
  int Nz = lattSize[2];
  int Nt = lattSize[3];

  int lX1h = subLattSize[0]/2;
  int lY = subLattSize[1];
  int lZ = subLattSize[2];
  int lT = subLattSize[3];

  // Diagnostic information:
  masterPrintf("VECLEN=%d SOALEN=%d\n", V, S);
  masterPrintf("Global Lattice Size = ");
  for(int mu=0; mu < 4; mu++){
     masterPrintf(" %d", lattSize[mu]);
  }
  masterPrintf("\n");

  masterPrintf("Local Lattice Size = ");
  for(int mu=0; mu < 4; mu++){
    masterPrintf(" %d", subLattSize[mu]);
  }
  masterPrintf("\n");

  masterPrintf("Block Sizes: By= %d Bz=%d\n", By, Bz);
  masterPrintf("Cores = %d\n", NCores);
  masterPrintf("SMT Grid: Sy=%d Sz=%d\n", Sy, Sz);
  masterPrintf("Pad Factors: PadXY=%d PadXYZ=%d\n", PadXY, PadXYZ);
  masterPrintf("Threads_per_core = %d\n", N_simt);
  masterPrintf("Initializing QPhiX Dslash\n");

  // Create Scalar Dslash Class
  double t_boundary = (FT)(1);
  double coeff_s    = (FT)(1);
  double coeff_t    = (FT)(1);
  Geometry<FT,V,S,compress> geom(subLattSize, By, Bz, NCores, Sy, Sz, PadXY, PadXYZ, MinCt);
  Dslash<FT,V,S,compress> DQPhiX(&geom, t_boundary, coeff_s, coeff_t);


	/************************
	 *                      *
	 *     GAUGE FIELDS     *
	 *                      *
	************************/

  // Allocate data for the gauges
  QGauge* u_packed[2];
  QGauge *packed_gauge_cb0 = (QGauge*) geom.allocCBGauge(); // Links emanating from ODD sites (cb=0)
  QGauge *packed_gauge_cb1 = (QGauge*) geom.allocCBGauge(); // Links emanating from EVEN sites (cb=1)
  u_packed[0] = packed_gauge_cb0;
  u_packed[1] = packed_gauge_cb1;

	// Initialise unit QGauge field
  masterPrintf("Initializing Unit QGauge Field: ");
  int nvecs = geom.nVecs();
  int nyg = geom.nGY();
  int Pxy = geom.getPxy();
  int Pxyz = geom.getPxyz();

  double start = omp_get_wtime();
	reorder_gauge_toQphix((double*) u_packed[0], (double*) u_packed[1]); // uses global tmlQCD gauge field as input
	double end = omp_get_wtime();
	masterPrintf(" QGauge reordering took: %g sec\n", end - start);


	/************************
	 *                      *
	 *     SPINOR FIELDS    *
	 *                      *
	************************/

	// Allocate data for the spinors
	QSpinor* qphix_in[2];
	QSpinor* qphix_out[2];
  QSpinor *packed_spinor_in_cb0 = (QSpinor*) geom.allocCBFourSpinor();
  QSpinor *packed_spinor_in_cb1 = (QSpinor*) geom.allocCBFourSpinor();
  QSpinor *packed_spinor_out_cb0 = (QSpinor*) geom.allocCBFourSpinor();
  QSpinor *packed_spinor_out_cb1 = (QSpinor*) geom.allocCBFourSpinor();
  qphix_in[0]  = packed_spinor_in_cb0;
  qphix_in[1]  = packed_spinor_in_cb1;
  qphix_out[0] = packed_spinor_out_cb0;
  qphix_out[1] = packed_spinor_out_cb1;

	// Reorder input spinor from tmLQCD to QPhiX
	reorder_spinor_toQphix( (double*) tmlqcd_in, (double*) qphix_in[0], (double*) qphix_in[1] );

	masterPrintf("Zeroing out output spinor: ");
	start = omp_get_wtime();
#pragma omp parallel for collapse(4)
	for(int cb=0; cb < 2; cb++) {
    for(int t=0; t < lT; t++) {
      for(int z=0; z < lZ; z++) {
        for(int y=0; y < lY; y++) {
          for(int s=0; s < nvecs; s++) {
            for(int spin=0; spin < 4; spin++) {
              for(int col=0; col < 3; col++)  {
                for(int x=0; x < S; x++) {
                  int ind = t*Pxyz+z*Pxy+y*nvecs+s; //((t*Nz+z)*Ny+y)*nvecs+s;
                  qphix_out[cb][ind][col][spin][0][x] = rep<FT,double> (0);
                  qphix_out[cb][ind][col][spin][1][x] = rep<FT,double> (0);
                }
              }
            }
          }
        }
      }
    }
	}
	end = omp_get_wtime();
	masterPrintf(" %g sec\n", end -start);

	// Apply QPhiX Dslash to qphix_in spinors
	DQPhiX.dslash(qphix_out[1], qphix_in[0], u_packed[1], /* isign == non-conjugate */ 1, /* cb == */ 1);
	DQPhiX.dslash(qphix_out[0], qphix_in[1], u_packed[0], /* isign == non-conjugate */ 1, /* cb == */ 0);

	// Reorder spinor fields back to tmLQCD
	reorder_spinor_fromQphix((double*) tmlqcd_out, (double*) qphix_out[0], (double*) qphix_out[1], (1.*g_kappa));
	reorder_spinor_fromQphix((double*) tmlqcd_in,  (double*) qphix_in[0],  (double*) qphix_in[1]);


	masterPrintf("Cleaning up\n");

	geom.free(packed_gauge_cb0);
	geom.free(packed_gauge_cb1);
	geom.free(packed_spinor_in_cb0);
	geom.free(packed_spinor_in_cb1);
	geom.free(packed_spinor_out_cb0);
	geom.free(packed_spinor_out_cb1);
}

int invert_qphix(spinor * const odd_out, spinor * const odd_in, const int max_iter, double eps_sq, const int rel_prec)
{
	if( precision == DOUBLE_PREC ) {
		if ( QPHIX_SOALEN > VECLEN_DP ) {
			masterPrintf("SOALEN=%d is greater than the double prec VECLEN=%d\n", QPHIX_SOALEN, VECLEN_DP);
			abort();
		}
		masterPrintf("TIMING IN DOUBLE PRECISION \n");
		if ( compress12 ) {
			invert<double,VECLEN_DP,QPHIX_SOALEN,true>(odd_out,odd_in,max_iter,eps_sq,rel_prec);
		}
		else {
			invert<double,VECLEN_DP,QPHIX_SOALEN,false>(odd_out,odd_in,max_iter,eps_sq,rel_prec);
		}
	}
	else if ( precision == FLOAT_PREC ) {
		if ( QPHIX_SOALEN > VECLEN_SP ) {
			masterPrintf("SOALEN=%d is greater than the single prec VECLEN=%d\n", QPHIX_SOALEN,VECLEN_SP);
			abort();
		}
		masterPrintf("TIMING IN SINGLE PRECISION \n");
		if ( compress12 ) {
			invert<float,VECLEN_SP,QPHIX_SOALEN,true>(odd_out,odd_in,max_iter,eps_sq,rel_prec);
		}
		else {
			invert<float,VECLEN_SP,QPHIX_SOALEN,false>(odd_out,odd_in,max_iter,eps_sq,rel_prec);
		}
	}
#if defined(QPHIX_MIC_SOURCE)
	else if ( precision == HALF_PREC ) {
		if ( QPHIX_SOALEN > VECLEN_HP ) {
			masterPrintf("SOALEN=%d is greater than the single prec VECLEN=%d\n", QPHIX_SOALEN,VECLEN_SP);
			abort();
		}
		masterPrintf("TIMING IN HALF PRECISION \n");
		if ( compress12 ) {
			invert<half,VECLEN_HP,QPHIX_SOALEN,true>(odd_out,odd_in,max_iter,eps_sq,rel_prec);
		}
		else {
			invert<half,VECLEN_HP,QPHIX_SOALEN,false>(odd_out,odd_in,max_iter,eps_sq,rel_prec);
		}
	}
#endif

	return 0;
}

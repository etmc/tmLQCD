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
*
* Notes:
*
* Minimum QUDA version is 0.7.0 (see https://github.com/lattice/qphix/issues/151
* and https://github.com/lattice/qphix/issues/157).
*
* To enable compilation of the same code for QUDA usage and standard non-QUDA usage, 
* all calls of these functions should be wrapped in precompiler switches of the form
*
*   #ifdef QUDA
*     ...
*   #endif  
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

int subLattSize[4] = {4,4,4,4};
int lattSize[4] = {4,4,4,4};

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

// gauge and invert paramameter structs; init. in _initQphix()
//QudaGaugeParam  gauge_param;
//QudaInvertParam inv_param;

// pointer to a temp. spinor, used for reordering etc.
double *tempSpinor;
// needed if even_odd_flag set
double *fullSpinor1;
double *fullSpinor2;

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

#if 0
template<typename FT, int V, int S, bool compress>
void
runTest(const int lattSize[], const int qmp_geom[])
{

  typedef typename Geometry<FT,V,S,compress>::SU3MatrixBlock Gauge;
  typedef typename Geometry<FT,V,S,compress>::FourSpinorBlock Spinor;

  bool verbose = false;

  // Work out local lattice size
  for(int mu=0; mu < 4; mu++){
    subLattSize[mu]=lattSize[mu]/qmp_geom[mu];
  }

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


  masterPrintf("Initializing Dslash\n");

  double t_boundary=(FT)(1);
  double coeff_s = (FT)(1);
  double coeff_t = (FT)(1);

  // Create Scalar Dslash Class
  Geometry<FT,V,S,compress> geom(subLattSize, By, Bz, NCores, Sy, Sz, PadXY, PadXYZ, MinCt);

//  Dslash<FT,V,S,compress> D32(&geom, t_boundary, coeff_s,coeff_t);


  // Allocate data for the gauges
  Gauge* packed_gauge_cb0 = (Gauge*)geom.allocCBGauge();
  Gauge* packed_gauge_cb1 = (Gauge*)geom.allocCBGauge();


  Gauge* u_packed[2];
  u_packed[0] = packed_gauge_cb0;
  u_packed[1] = packed_gauge_cb1;

  double factor=0.08;

  masterPrintf("Initializing Fake Gauge Field: ");
  int nvecs = geom.nVecs();
  int nyg = geom.nGY();
  int Pxy = geom.getPxy();
  int Pxyz = geom.getPxyz();

  double start = omp_get_wtime();

#pragma omp parallel for collapse(4)
  for(int t = 0; t < lT; t++) {
    for(int z = 0; z < lZ; z++) {
      for(int y = 0; y < lY; y++) {
	for(int s = 0; s < nvecs; s++) {
	  for(int mu = 0; mu < 8; mu++) {
	    for(int c = 0; c < (compress ? 2 : 3) ; c++) {
	      for(int c2 = 0; c2 < 3; c2++) {
		for(int x = 0; x < S; x++) {

		  int block = (t*Pxyz+z*Pxy)/nyg+(y/nyg)*nvecs+s;

		  // This will work out to be between 0 and veclen
		  int xx = (y%nyg)*S+x;

		  double d1=factor*(drand48()-0.5);
		  double d2=factor*(drand48()-0.5);
		  double d3=factor*(drand48()-0.5);
		  double d4=factor*(drand48()-0.5);

		  if( c == c2 ) {
		    u_packed[0][block][mu][c][c2][RE][xx]=rep<FT,double>((double)1 + d1);
		    u_packed[1][block][mu][c][c2][RE][xx]=rep<FT,double>((double)1 + d3);
		  }
		  else {
		     u_packed[0][block][mu][c][c2][RE][xx]=rep<FT,double>(d1);
		     u_packed[1][block][mu][c][c2][RE][xx]=rep<FT,double>(d3);
		  }

		  u_packed[0][block][mu][c][c2][IM][xx]=rep<FT,double>(d2);
		  u_packed[1][block][mu][c][c2][IM][xx]=rep<FT,double>(d4);

		}
	      }
	    } // row
	  }
	}
      }
    }
  }

  if ( !compress ) {
#pragma omp parallel for collapse(4)
    for(int t = 0; t < lT; t++) {
      for(int z = 0; z < lZ; z++) {
	for(int y = 0; y < lY; y++) {
	  for(int s = 0; s < nvecs; s++) {

	    int block = (t*Pxyz+z*Pxy)/nyg+(y/nyg)*nvecs+s;

	    for(int mu = 0; mu < 8; mu++) {
	      for(int row =0; row < 2; row++) {

		double norm_row_cb0[V];
		double norm_row_cb1[V];

		for(int x = 0; x < V; x++) {
		  norm_row_cb0[x]=0;
		  norm_row_cb1[x]=0;
		}

		// This will work out to be between 0 and veclen
		// Accumulate the norms
		for(int col=0; col < 3; col++) {
		  for(int x=0; x < S; x++){
		    int xx = (y%nyg)*S+x;

		    double u0_re = rep<double,FT>( u_packed[0][block][mu][row][col][RE][xx] );
		    double u1_re = rep<double,FT>( u_packed[1][block][mu][row][col][RE][xx] );

		    double u0_im = rep<double,FT>( u_packed[0][block][mu][row][col][IM][xx] );
		    double u1_im = rep<double,FT>( u_packed[1][block][mu][row][col][IM][xx] );

		    norm_row_cb0[xx] += ( u0_re
					 *u0_re)
		      +(u0_im
			*u0_im);

		    norm_row_cb1[xx]+= (u1_re*u1_re) + (u1_im*u1_im);

		  } // x
		} // col

		for(int x=0; x < V; x++) {
		  norm_row_cb0[x] = sqrt(norm_row_cb0[x]);
		  norm_row_cb1[x] = sqrt(norm_row_cb1[x]);
		}

		// Normalize each component.

		for(int col=0; col < 3; col++) {
		  for(int x=0; x < S; x++) {
		    int xx = (y%nyg)*S+x;

		    double u0_re = rep<double,FT>( u_packed[0][block][mu][row][col][RE][xx] )/norm_row_cb0[xx];
		    double u1_re = rep<double,FT>( u_packed[1][block][mu][row][col][RE][xx] )/norm_row_cb1[xx];

		    double u0_im = rep<double,FT>( u_packed[0][block][mu][row][col][IM][xx] )/norm_row_cb0[xx];
		    double u1_im = rep<double,FT>( u_packed[1][block][mu][row][col][IM][xx] )/norm_row_cb1[xx];

		    u_packed[0][block][mu][row][col][RE][xx]=rep<FT,double>(u0_re);
		    u_packed[0][block][mu][row][col][IM][xx]=rep<FT,double>(u0_im);

		    u_packed[1][block][mu][row][col][RE][xx]=rep<FT,double>(u1_re);
		    u_packed[1][block][mu][row][col][IM][xx]=rep<FT,double>(u1_im);
		  } // x
		} // col
	      } // row

	      {
		for(int x=0; x < S; x++) {
		  // 3rd row reconstruction.
		  int xx = (y%nyg)*S+x;
		  double ar=rep<double,FT>(u_packed[0][block][mu][0][0][RE][xx]);
		  double ai=rep<double,FT>(u_packed[0][block][mu][0][0][IM][xx]);

		  double br=rep<double,FT>(u_packed[0][block][mu][0][1][RE][xx]);
		  double bi=rep<double,FT>(u_packed[0][block][mu][0][1][IM][xx]);

		  double cr=rep<double,FT>(u_packed[0][block][mu][0][2][RE][xx]);
		  double ci=rep<double,FT>(u_packed[0][block][mu][0][2][IM][xx]);

		  double dr=rep<double,FT>(u_packed[0][block][mu][1][0][RE][xx]);
		  double di=rep<double,FT>(u_packed[0][block][mu][1][0][IM][xx]);

		  double er=rep<double,FT>(u_packed[0][block][mu][1][1][RE][xx]);
		  double ei=rep<double,FT>(u_packed[0][block][mu][1][1][IM][xx]);

		  double fr=rep<double,FT>(u_packed[0][block][mu][1][2][RE][xx]);
		  double fi=rep<double,FT>(u_packed[0][block][mu][1][2][IM][xx]);

		  u_packed[0][block][mu][2][0][RE][xx]=rep<FT,double>(br*fr-bi*fi-er*cr+ei*ci);
		  u_packed[0][block][mu][2][0][IM][xx]=rep<FT,double>(er*ci+ei*cr-br*fi-bi*fr);
		  u_packed[0][block][mu][2][1][RE][xx]=rep<FT,double>(dr*cr-di*ci-ar*fr+ai*fi);
		  u_packed[0][block][mu][2][1][IM][xx]=rep<FT,double>(ar*fi+ai*fr-dr*ci-di*cr);
		  u_packed[0][block][mu][2][2][RE][xx]=rep<FT,double>(ar*er-ai*ei-dr*br+di*bi);
		  u_packed[0][block][mu][2][2][IM][xx]=rep<FT,double>(dr*bi+di*br-ar*ei-ai*er);
		}
	      }

	      {
		for(int x=0; x < S; x++) {
		  int xx = (y%nyg)*S+x;
		  // 3rd row reconstruction.
		  double ar=rep<double,FT>(u_packed[1][block][mu][0][0][RE][xx]);
		  double ai=rep<double,FT>(u_packed[1][block][mu][0][0][IM][xx]);

		  double br=rep<double,FT>(u_packed[1][block][mu][0][1][RE][xx]);
		  double bi=rep<double,FT>(u_packed[1][block][mu][0][1][IM][xx]);

		  double cr=rep<double,FT>(u_packed[1][block][mu][0][2][RE][xx]);
		  double ci=rep<double,FT>(u_packed[1][block][mu][0][2][IM][xx]);

		  double dr=rep<double,FT>(u_packed[1][block][mu][1][0][RE][xx]);
		  double di=rep<double,FT>(u_packed[1][block][mu][1][0][IM][xx]);

		  double er=rep<double,FT>(u_packed[1][block][mu][1][1][RE][xx]);
		  double ei=rep<double,FT>(u_packed[1][block][mu][1][1][IM][xx]);

		  double fr=rep<double,FT>(u_packed[1][block][mu][1][2][RE][xx]);
		  double fi=rep<double,FT>(u_packed[1][block][mu][1][2][IM][xx]);

		  u_packed[1][block][mu][2][0][RE][xx]=rep<FT,double>(br*fr-bi*fi-er*cr+ei*ci);
		  u_packed[1][block][mu][2][0][IM][xx]=rep<FT,double>(er*ci+ei*cr-br*fi-bi*fr);
		  u_packed[1][block][mu][2][1][RE][xx]=rep<FT,double>(dr*cr-di*ci-ar*fr+ai*fi);
		  u_packed[1][block][mu][2][1][IM][xx]=rep<FT,double>(ar*fi+ai*fr-dr*ci-di*cr);
		  u_packed[1][block][mu][2][2][RE][xx]=rep<FT,double>(ar*er-ai*ei-dr*br+di*bi);
		  u_packed[1][block][mu][2][2][IM][xx]=rep<FT,double>(dr*bi+di*br-ar*ei-ai*er);
		} // x
	      }

	    } // mu
	  } // s
	} // y
      } // z
    } // t

  } // end if ! compress

  double end = omp_get_wtime();
  masterPrintf(" %g sec\n", end - start);
  // Allocate data for the spinors

  Spinor* p_even=(Spinor*)geom.allocCBFourSpinor();
  Spinor* p_odd=(Spinor*)geom.allocCBFourSpinor();
  Spinor* c_even=(Spinor*)geom.allocCBFourSpinor();
  Spinor* c_odd=(Spinor*)geom.allocCBFourSpinor();


  // Point to the second block of the array. Now there is padding on both ends.
  Spinor *psi_s[2] = { p_even, p_odd };
  Spinor *chi_s[2] = { c_even, c_odd };


  masterPrintf("Filling Input spinor: ");


  start=omp_get_wtime();
#pragma omp parallel for collapse(4)
  for(int t=0; t < lT; t++) {
    for(int z=0; z < lZ; z++) {
      for(int y=0; y < lY; y++) {
	for(int s=0; s < nvecs; s++) {
	  for(int spin=0; spin < 4; spin++) {
	    for(int col=0; col < 3; col++)  {
	      for(int x=0; x < S; x++) {
		double d1=drand48()-0.5;
		double d2=drand48()-0.5;
		double d3=drand48()-0.5;
		double d4=drand48()-0.5;

		int ind = t*Pxyz+z*Pxy+y*nvecs+s; //((t*Nz+z)*Ny+y)*nvecs+s;
		psi_s[0][ind][col][spin][0][x] = rep<FT,double>(d1);
		psi_s[0][ind][col][spin][1][x] = rep<FT,double>(d2);
		psi_s[1][ind][col][spin][0][x] = rep<FT,double>(d3);
		psi_s[1][ind][col][spin][1][x] = rep<FT,double>(d4);
	      }
	    }
	  }
	}
      }
    }
  }

  end = omp_get_wtime();
  masterPrintf(" %g sec\n", end - start);

  masterPrintf("Zeroing output spinor: ");
  start = omp_get_wtime();
#pragma omp parallel for collapse(4)
  for(int t=0; t < lT; t++) {
    for(int z=0; z < lZ; z++) {
      for(int y=0; y < lY; y++) {
	for(int s=0; s < nvecs; s++) {
	  for(int spin=0; spin < 4; spin++) {
	    for(int col=0; col < 3; col++)  {
	      for(int x=0; x < S; x++) {
		double d=0;
		int ind = t*Pxyz+z*Pxy+y*nvecs+s; //((t*Nz+z)*Ny+y)*nvecs+s;
		chi_s[0][ind][col][spin][0][x] = rep<FT,double>(d);
		chi_s[0][ind][col][spin][1][x] = rep<FT,double>(d);
		chi_s[1][ind][col][spin][0][x] = rep<FT,double>(d);
		chi_s[1][ind][col][spin][1][x] = rep<FT,double>(d);
	      }
	    }
	  }
	}
      }
    }
  }
  end = omp_get_wtime();
  masterPrintf(" %g sec\n", end -start);

#if 0
  // Go through the test cases -- apply SSE dslash versus, QDP Dslash
  for(int isign=1; isign >= -1; isign -=2) {
    for(int cb=0; cb < 2; cb++) {
      int source_cb = 1 - cb;
      int target_cb = cb;
      masterPrintf("Timing on cb=%d isign=%d\n", cb, isign);
      masterPrintf("=============================\n");

      for(int repeat=0; repeat < 3; repeat++) {
	double start = omp_get_wtime();

	for(int i=0; i < iters; i++) {
	  // Apply Optimized Dslash
	  D32.dslash(chi_s[target_cb],
		     psi_s[source_cb],
		     u_packed[target_cb],
		     isign,
		     target_cb);
	}

	double end = omp_get_wtime();
	double time = end - start;
	CommsUtils::sumDouble(&time);
	time /= (double)CommsUtils::numNodes();

	masterPrintf("\t timing %d of 3\n", repeat);
	masterPrintf("\t %d iterations in %e seconds\n", iters, time);
	masterPrintf("\t %e usec/iteration\n", 1.0e6*time/(double)iters);
	double Gflops = 1320.0f*(double)(iters)*(double)(X1h*Ny*Nz*Nt)/1.0e9;
	double perf = Gflops/time;
	masterPrintf("\t Performance: %g GFLOPS total\n", perf);
      }
    }
  }
#endif


#if 1
  masterPrintf("Creating Wilson Op\n");
  double Mass=0.1;
  EvenOddWilsonOperator<FT, V, S,compress> M(Mass, u_packed, &geom, t_boundary, coeff_s, coeff_t);

  // Go through the test cases -- apply SSE dslash versus, QDP Dslash
  for(int isign=1; isign >= -1; isign -=2) {

    masterPrintf("Timing M: isign=%d\n",  isign);
    masterPrintf("=============================\n");

    for(int repeat=0; repeat < 3; repeat++) {
      double start = omp_get_wtime();

      for(int i=0; i < iters; i++) {
	// Apply Optimized Dslash
	M(chi_s[0],
	  psi_s[0],
	  isign);
      }

      double end = omp_get_wtime();
      double time = end - start;

      CommsUtils::sumDouble(&time);
      time /= (double)CommsUtils::numNodes();

      masterPrintf("\t timing %d of 3\n", repeat);
      masterPrintf("\t %d iterations in %e seconds\n", iters, time);
      masterPrintf("\t %e usec/iteration\n", 1.0e6*time/(double)iters);
      double flops_per_iter = 1320.0f*2.0 + 24.0*3.0;
      double Gflops = flops_per_iter*(double)(iters)*(double)(X1h*Ny*Nz*Nt)/1.0e9;
      double perf = Gflops/time;
      masterPrintf("\t Performance: %g GFLOPS total\n", perf);
      masterPrintf("\t              %g GFLOPS / node\n", perf/(double)CommsUtils::numNodes());

    }

  }
#endif

#if 1
  double rsd_target=rsdTarget<FT>::value;
  int max_iters=5000;
  int niters;
  double rsd_final;
  int len = (geom.getPxyz()*geom.Nt()*sizeof(Spinor))/sizeof(FT);
  FT *c_s0 = (FT *)chi_s[0];

  {
    masterPrintf("Creating Solver\n");
    InvCG<FT,V,S, compress> solver(M, max_iters);

    masterPrintf("Tuning Solver\n");
    solver.tune();

    for(int solve = 0; solve < 5; solve++ ) {
      masterPrintf("Starting solver\n");
      unsigned long site_flops=0;
      unsigned long mv_apps=0;

      FT *psi_0 = (FT *)psi_s[0];
#if defined(__INTEL_COMPILER)
#pragma simd
#endif
#pragma omp parallel for
      for(int i=0; i < len; i++) {
	c_s0[i] = rep<FT,double>(0);
      psi_0[i]= rep<FT,double>(0);
      }

#pragma omp parallel for collapse(4)
      for(int t=0; t < lT; t++) {
	for(int z=0; z < lZ; z++) {
	  for(int y=0; y < lY; y++) {
	    for(int s=0; s < nvecs; s++) {
	      for(int spin=0; spin < 4; spin++) {
		for(int col=0; col < 3; col++)  {
		  for(int x=0; x < S; x++) {

		    int ind = t*Pxyz+z*Pxy+y*nvecs+s; //((t*Nz+z)*Ny+y)*nvecs+s;
		    int x_coord = s*S + x;
		    double d1 = drand48()-0.5;
		    double d2 = drand48()-0.5;
		    double d3 = drand48()-0.5;
		    double d4 = drand48()-0.5;

		    psi_s[0][ind][col][spin][0][x] = rep<FT,double>(d1);
		    psi_s[0][ind][col][spin][1][x] = rep<FT,double>(d2);
		    psi_s[1][ind][col][spin][0][x] = rep<FT,double>(d3);
		    psi_s[1][ind][col][spin][1][x] = rep<FT,double>(d4);
		  }
		}
	      }
	    }
	  }
	}
      }

      start = omp_get_wtime();
      solver(chi_s[0], psi_s[0], rsd_target, niters, rsd_final, site_flops, mv_apps,1,verbose);
      end = omp_get_wtime();


      unsigned long num_cb_sites=X1h*Ny*Nz*Nt;
      unsigned long total_flops = (site_flops + (72+2*1320)*mv_apps)*num_cb_sites;
      masterPrintf("Solver Time=%g(s)\n", (end-start));
      masterPrintf("CG GFLOPS=%g\n", 1.0e-9*(double)(total_flops)/(end -start));

    }
  } // Solver
#endif

#if 0
  {
    masterPrintf("Creating BiCGStab Solver\n");
    InvBiCGStab<FT,V,S,compress> solver2(M, max_iters);
    masterPrintf("Tuning BiCGStab Solver\n");
    solver2.tune();

    for(int solve =0; solve < 5; solve++) {
      unsigned long site_flops;
      unsigned long mv_apps;

      FT *psi_0 = (FT *)psi_s[0];
#if defined(__INTEL_COMPILER)
#pragma simd
#endif
#pragma omp parallel for
      for(int i=0; i < len; i++) {
	c_s0[i] = rep<FT,double>(0);
	psi_0[i]= rep<FT,double>(0);
      }

#pragma omp parallel for collapse(4)
      for(int t=0; t < lT; t++) {
	for(int z=0; z < lZ; z++) {
	  for(int y=0; y < lY; y++) {
	    for(int s=0; s < nvecs; s++) {
	      for(int spin=0; spin < 4; spin++) {
		for(int col=0; col < 3; col++)  {
		  for(int x=0; x < S; x++) {

		    int ind = t*Pxyz+z*Pxy+y*nvecs+s; //((t*Nz+z)*Ny+y)*nvecs+s;
		    int x_coord = s*S + x;

		    double d1 = drand48()-0.5;
		    double d2 = drand48()-0.5;
		    double d3 = drand48()-0.5;
		    double d4 = drand48()-0.5;

		    psi_s[0][ind][col][spin][0][x] = rep<FT,double>(d1);
		    psi_s[0][ind][col][spin][1][x] = rep<FT,double>(d2);
		    psi_s[1][ind][col][spin][0][x] = rep<FT,double>(d3);
		    psi_s[1][ind][col][spin][1][x] = rep<FT,double>(d4);

		  }
		}
	      }
	    }
	  }
	}
      }

      start = omp_get_wtime();
      solver2(chi_s[0], psi_s[0], rsd_target, niters, rsd_final, site_flops, mv_apps,1,verbose);
      end = omp_get_wtime();


      unsigned long num_cb_sites=X1h*Ny*Nz*Nt;
      unsigned long total_flops = (site_flops + (72+2*1320)*mv_apps)*num_cb_sites;

      masterPrintf("Solver Time=%g(s)\n", (end-start));
      masterPrintf("BICGSTAB GFLOPS=%g\n", 1.0e-9*(double)(total_flops)/(end -start));
    }
  }
#endif

  masterPrintf("Cleaning up\n");


  geom.free(packed_gauge_cb0);
  geom.free(packed_gauge_cb1);
  geom.free(p_even);
  geom.free(p_odd);
  geom.free(c_even);
  geom.free(c_odd);
}
#endif

#if 0
  template<typename FT, int veclen, int soalen, bool compress, typename QDPGauge>
    void qdp_pack_gauge(const QDPGauge& u,
			typename Geometry<FT,veclen,soalen,compress>::SU3MatrixBlock *u_cb0,
			typename Geometry<FT,veclen,soalen,compress>::SU3MatrixBlock *u_cb1,
			Geometry<FT,veclen,soalen,compress>& s)
  {
    // Get the subgrid latt size.
    int Nt = s.Nt();
    int Nz = s.Nz();
    int Ny = s.Ny();
    int nvecs = s.nVecs();
    int nyg = s.nGY();
    int Pxy = s.getPxy();
    int Pxyz = s.getPxyz();


    // Shift the lattice to get U(x-mu)
    QDPGauge u_minus(4);
    for(int mu=0; mu < 4; mu++) {
      u_minus[mu] = shift(u[mu], BACKWARD, mu);
    }


#pragma omp parallel for collapse(4)
    for(int t = 0; t < Nt; t++) {
      for(int z = 0; z < Nz; z++) {
	for(int y = 0; y < Ny; y++) {
	  for(int s = 0; s < nvecs; s++) {
	    for(int mu = 0; mu < 4; mu++) {
	      int outer_c = 3;
	      if ( compress ) {
		outer_c = 2;
	      }
	      for(int c = 0; c < outer_c; c++) {
		for(int c2 = 0; c2 < 3; c2++) {
		  for(int x = 0; x < soalen; x++) {

		    //#ifndef USE_PACKED_GAUGES
		    //int xx = x;
		    //int block = ((t*Nz+z)*Ny+y)*nvecs+s;

		    //#endif
		    //#else // USE_PACKED_GAUGES
		    int block = (t*Pxyz+z*Pxy)/nyg+(y/nyg)*nvecs+s;
		    int xx = (y%nyg)*soalen+x;
		    // #endif // USE_PACKED_GAUGES

		    int qdpsite = x + soalen*(s + nvecs*(y + Ny*(z + Nz*t)));
		    u_cb0[block][2*mu][c][c2][0][xx] = u_minus[mu].elem(rb[0].start() + qdpsite).elem().elem(c2,c).real();
		    u_cb0[block][2*mu][c][c2][1][xx] = u_minus[mu].elem(rb[0].start() + qdpsite).elem().elem(c2,c).imag();
		    u_cb0[block][2*mu+1][c][c2][0][xx] = u[mu].elem(rb[0].start() + qdpsite).elem().elem(c2,c).real();
		    u_cb0[block][2*mu+1][c][c2][1][xx] = u[mu].elem(rb[0].start() + qdpsite).elem().elem(c2,c).imag();


		    u_cb1[block][2*mu][c][c2][0][xx] = u_minus[mu].elem(rb[1].start() + qdpsite).elem().elem(c2,c).real();
		    u_cb1[block][2*mu][c][c2][1][xx] = u_minus[mu].elem(rb[1].start() + qdpsite).elem().elem(c2,c).imag();
		    u_cb1[block][2*mu+1][c][c2][0][xx] = u[mu].elem(rb[1].start() + qdpsite).elem().elem(c2,c).real();
		    u_cb1[block][2*mu+1][c][c2][1][xx] = u[mu].elem(rb[1].start() + qdpsite).elem().elem(c2,c).imag();
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
#endif

  template<typename FT, int veclen, int soalen, bool compress>
  void qdp_pack_cb_spinor(const double* psi_in,
			  typename Geometry<FT,veclen,soalen, compress>::FourSpinorBlock* psi,
			  Geometry<FT,veclen,soalen,compress>& s,
			  int cb)
  {
    int Nt = s.Nt();
    int Nz = s.Nz();
    int Ny = s.Ny();
    int Nxh = s.Nxh();
    int nvecs = s.nVecs();
    int Pxy = s.getPxy();
    int Pxyz = s.getPxyz();

//    double pionr[Nt];
//    for( int t = 0; t < Nt; t++ )
//		pionr[t] = 0.0;

#pragma omp parallel for collapse(4)
    for(int t=0; t < Nt; t++) {
      for(int z=0; z < Nz; z++) {
	for(int y=0; y < Ny; y++) {
	  for(int s=0; s < nvecs; s++) {
	    for(int spin=0; spin < 4; spin++) {
	      for(int col=0; col < 3; col++)  {
		for(int x=0; x < soalen; x++) {

		  int ind = t*Pxyz+z*Pxy+y*nvecs+s; //((t*Nz+z)*Ny+y)*nvecs+s;
		  int x_coord = s*soalen + x;
		  int qdp_ind = ((t*Nz + z)*Ny + y)*Nxh + x_coord;

		  int oddBit = (t+y+z) & 1;
		  int evenBit = (oddBit?0:1);

		  int tm_t = t;
		  int tm_z = z;
		  int tm_y = y;
		  int tm_x = x_coord*2+oddBit;
//
//		  if( tm_t%2 != cb ) tm_z++;
//		  if( tm_z%2 != cb ) tm_y++;
//		  if( tm_y%2 != cb ) tm_x++;

		  int tm_idx = g_ipt[tm_t][tm_x][tm_y][tm_z];
		  int tm_ieo = g_lexic2eosub[ tm_idx ];

		    psi[ind][col][spin][0][x] = psi_in[24*tm_ieo+6*spin+2*col+0];
		    		//psi_in.elem(rb[cb].start()+qdp_ind).elem(spin).elem(col).real();
		    psi[ind][col][spin][1][x] = psi_in[24*tm_ieo+6*spin+2*col+0];
		    		//psi_in.elem(rb[cb].start()+qdp_ind).elem(spin).elem(col).imag();

		  }
		}
	      }
	    }
	  }
	}
      }

  }

//  template<typename FT, int veclen, int soalen, bool compress>
//  void qdp_pack_spinor(const double** psi_in,
//		       typename Geometry<FT,veclen,soalen, compress>::FourSpinorBlock* psi_even,
//		       typename Geometry<FT,veclen,soalen, compress>::FourSpinorBlock* psi_odd,
//		       Geometry<FT,veclen,soalen,compress>& s)
//  {
//    qdp_pack_cb_spinor(psi_in,psi_even,s,0);
//    qdp_pack_cb_spinor(psi_in,psi_odd,s,1);
//  }

  template<typename FT, int veclen, int soalen, bool compress>
    void qdp_unpack_cb_spinor(typename Geometry<FT,veclen,soalen,compress>::FourSpinorBlock* chi_packed,
			      double* chi,
			      Geometry<FT,veclen,soalen,compress>& s,
			      int cb)
  {
    int Nt = s.Nt();
    int Nz = s.Nz();
    int Ny = s.Ny();
    int Nxh = s.Nxh();
    int nvecs = s.nVecs();
    int Pxy = s.getPxy();
    int Pxyz = s.getPxyz();

//    double pionr[Nt];
//    for( int t = 0; t < Nt; t++ )
//		pionr[t] = 0.0;

#pragma omp parallel for collapse(4)
    for(int t=0; t < Nt; t++) {
      for(int z=0; z < Nz; z++) {
	for(int y=0; y < Ny; y++) {
	  for(int s=0; s < nvecs; s++) {
	    for(int spin=0; spin < 4; spin++) {
	      for(int col=0; col < 3; col++)  {
		for(int x=0; x < soalen; x++) {

		  int ind = t*Pxyz+z*Pxy+y*nvecs+s; //((t*Nz+z)*Ny+y)*nvecs+s;
		  int x_coord = s*soalen + x;
		  int qdp_ind = ((t*Nz + z)*Ny + y)*Nxh + x_coord;

		  int oddBit = (t+y+z) & 1;
		  int evenBit = (oddBit?0:1);

		  int tm_t = t;
		  int tm_z = z;
		  int tm_y = y;
		  int tm_x = x_coord*2+oddBit;
//
//		  if( tm_t%2 != cb ) tm_z++;
//		  if( tm_z%2 != cb ) tm_y++;
//		  if( tm_y%2 != cb ) tm_x++;

		  int tm_idx = g_ipt[tm_t][tm_x][tm_y][tm_z];
		  int tm_ieo = g_lexic2eosub[ tm_idx ];

//		  if( tm_x==0 && tm_y==0 && tm_z==0 ) {
//			if(  col==0 && spin==0 )
//		  if( cb== 1)
//			  masterPrintf("%d\t%e\tcb=%d, t=%d,x=%d,y=%d,z=%d\n",t,chi_packed[ind][col][spin][0][x]*chi_packed[ind][col][spin][0][x]
//										   +chi_packed[ind][col][spin][1][x]*chi_packed[ind][col][spin][1][x],cb, t,tm_x,tm_y,tm_z);

//			  pionr[tm_t] = chi_packed[ind][col][spin][0][x]*chi_packed[ind][col][spin][0][x]
//				         +chi_packed[ind][col][spin][1][x]*chi_packed[ind][col][spin][1][x];
//		  }

		  //chi.elem(rb[cb].start()+qdp_ind).elem(spin).elem(col).real()
		  chi[24*tm_ieo+6*spin+2*col+0] = chi_packed[ind][col][spin][0][x];
		  //chi.elem(rb[cb].start()+qdp_ind).elem(spin).elem(col).imag()
		  chi[24*tm_ieo+6*spin+2*col+1] = chi_packed[ind][col][spin][1][x];

		}
	      }
	    }
	  }
	}
      }
    }
//    for( int t = 0; t < Nt; t++ )
//    	printf("%i\t%e\n", t, pionr[t]);

  }

//  template<typename FT, int veclen, int soalen, bool compress>
//    void qdp_unpack_spinor(typename Geometry<FT,veclen,soalen,compress>::FourSpinorBlock* chi_even,
//			   typename Geometry<FT,veclen,soalen,compress>::FourSpinorBlock* chi_odd,
//			   double* chi,
//			   Geometry<FT,veclen,soalen,compress>& s)
//  {
//    qdp_unpack_cb_spinor(chi_even,chi,s,0);
//    qdp_unpack_cb_spinor(chi_odd,chi,s,1);
//  }

template<typename FT, int V, int S, bool compress>
void
invert(spinor * const P, spinor * const Q, const int max_iter, double eps_sq, const int rel_prec )
{

  typedef typename Geometry<FT,V,S,compress>::SU3MatrixBlock Gauge;
  typedef typename Geometry<FT,V,S,compress>::FourSpinorBlock Spinor;

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


  masterPrintf("Initializing Dslash\n");

  double t_boundary=(FT)(1);
  double coeff_s = (FT)(1);
  double coeff_t = (FT)(1);

  // Create Scalar Dslash Class
  Geometry<FT,V,S,compress> geom(subLattSize, By, Bz, NCores, Sy, Sz, PadXY, PadXYZ, MinCt);

//  Dslash<FT,V,S,compress> D32(&geom, t_boundary, coeff_s,coeff_t);


  // Allocate data for the gauges
  Gauge* packed_gauge_cb0 = (Gauge*)geom.allocCBGauge();
  Gauge* packed_gauge_cb1 = (Gauge*)geom.allocCBGauge();


  Gauge* u_packed[2];
  u_packed[0] = packed_gauge_cb0;
  u_packed[1] = packed_gauge_cb1;

#if 0
  qdp_pack_gauge<>(g_gauge_field, packed_gauge_cb0,packed_gauge_cb1, geom);
#else

  double factor=0.08;

  masterPrintf("Initializing Fake Gauge Field: ");
  int nvecs = geom.nVecs();
  int nyg = geom.nGY();
  int Pxy = geom.getPxy();
  int Pxyz = geom.getPxyz();

  double start = omp_get_wtime();

#pragma omp parallel for collapse(4)
  for(int t = 0; t < lT; t++) {
    for(int z = 0; z < lZ; z++) {
      for(int y = 0; y < lY; y++) {
	for(int s = 0; s < nvecs; s++) {
	  for(int mu = 0; mu < 8; mu++) {
	    for(int c = 0; c < (compress ? 2 : 3) ; c++) {
	      for(int c2 = 0; c2 < 3; c2++) {
		for(int x = 0; x < S; x++) {

		  int block = (t*Pxyz+z*Pxy)/nyg+(y/nyg)*nvecs+s;

		  // This will work out to be between 0 and veclen
		  int xx = (y%nyg)*S+x;

		  double d1=factor*(drand48()-0.5);
		  double d2=factor*(drand48()-0.5);
		  double d3=factor*(drand48()-0.5);
		  double d4=factor*(drand48()-0.5);

		  if( c == c2 ) {
		    u_packed[0][block][mu][c][c2][RE][xx]=rep<FT,double>((double)1);//rep<FT,double>((double)1 + d1);
		    u_packed[1][block][mu][c][c2][RE][xx]=rep<FT,double>((double)1);//rep<FT,double>((double)1 + d3);
		  }
		  else {
		     u_packed[0][block][mu][c][c2][RE][xx]=rep<FT,double>((double)0);//rep<FT,double>(d1);
		     u_packed[1][block][mu][c][c2][RE][xx]=rep<FT,double>((double)0);//rep<FT,double>(d3);
		  }

		  u_packed[0][block][mu][c][c2][IM][xx]=rep<FT,double>((double)0);//rep<FT,double>(d2);
		  u_packed[1][block][mu][c][c2][IM][xx]=rep<FT,double>((double)0);//rep<FT,double>(d4);

		}
	      }
	    } // row
	  }
	}
      }
    }
  }

  if ( !compress ) {
#pragma omp parallel for collapse(4)
    for(int t = 0; t < lT; t++) {
      for(int z = 0; z < lZ; z++) {
	for(int y = 0; y < lY; y++) {
	  for(int s = 0; s < nvecs; s++) {

	    int block = (t*Pxyz+z*Pxy)/nyg+(y/nyg)*nvecs+s;

	    for(int mu = 0; mu < 8; mu++) {
	      for(int row =0; row < 2; row++) {

		double norm_row_cb0[V];
		double norm_row_cb1[V];

		for(int x = 0; x < V; x++) {
		  norm_row_cb0[x]=0;
		  norm_row_cb1[x]=0;
		}

		// This will work out to be between 0 and veclen
		// Accumulate the norms
		for(int col=0; col < 3; col++) {
		  for(int x=0; x < S; x++){
		    int xx = (y%nyg)*S+x;

		    double u0_re = rep<double,FT>( u_packed[0][block][mu][row][col][RE][xx] );
		    double u1_re = rep<double,FT>( u_packed[1][block][mu][row][col][RE][xx] );

		    double u0_im = rep<double,FT>( u_packed[0][block][mu][row][col][IM][xx] );
		    double u1_im = rep<double,FT>( u_packed[1][block][mu][row][col][IM][xx] );

		    norm_row_cb0[xx] += ( u0_re
					 *u0_re)
		      +(u0_im
			*u0_im);

		    norm_row_cb1[xx]+= (u1_re*u1_re) + (u1_im*u1_im);

		  } // x
		} // col

		for(int x=0; x < V; x++) {
		  norm_row_cb0[x] = sqrt(norm_row_cb0[x]);
		  norm_row_cb1[x] = sqrt(norm_row_cb1[x]);
		}

		// Normalize each component.

		for(int col=0; col < 3; col++) {
		  for(int x=0; x < S; x++) {
		    int xx = (y%nyg)*S+x;

		    double u0_re = rep<double,FT>( u_packed[0][block][mu][row][col][RE][xx] )/norm_row_cb0[xx];
		    double u1_re = rep<double,FT>( u_packed[1][block][mu][row][col][RE][xx] )/norm_row_cb1[xx];

		    double u0_im = rep<double,FT>( u_packed[0][block][mu][row][col][IM][xx] )/norm_row_cb0[xx];
		    double u1_im = rep<double,FT>( u_packed[1][block][mu][row][col][IM][xx] )/norm_row_cb1[xx];

		    u_packed[0][block][mu][row][col][RE][xx]=rep<FT,double>(u0_re);
		    u_packed[0][block][mu][row][col][IM][xx]=rep<FT,double>(u0_im);

		    u_packed[1][block][mu][row][col][RE][xx]=rep<FT,double>(u1_re);
		    u_packed[1][block][mu][row][col][IM][xx]=rep<FT,double>(u1_im);
		  } // x
		} // col
	      } // row

	      {
		for(int x=0; x < S; x++) {
		  // 3rd row reconstruction.
		  int xx = (y%nyg)*S+x;
		  double ar=rep<double,FT>(u_packed[0][block][mu][0][0][RE][xx]);
		  double ai=rep<double,FT>(u_packed[0][block][mu][0][0][IM][xx]);

		  double br=rep<double,FT>(u_packed[0][block][mu][0][1][RE][xx]);
		  double bi=rep<double,FT>(u_packed[0][block][mu][0][1][IM][xx]);

		  double cr=rep<double,FT>(u_packed[0][block][mu][0][2][RE][xx]);
		  double ci=rep<double,FT>(u_packed[0][block][mu][0][2][IM][xx]);

		  double dr=rep<double,FT>(u_packed[0][block][mu][1][0][RE][xx]);
		  double di=rep<double,FT>(u_packed[0][block][mu][1][0][IM][xx]);

		  double er=rep<double,FT>(u_packed[0][block][mu][1][1][RE][xx]);
		  double ei=rep<double,FT>(u_packed[0][block][mu][1][1][IM][xx]);

		  double fr=rep<double,FT>(u_packed[0][block][mu][1][2][RE][xx]);
		  double fi=rep<double,FT>(u_packed[0][block][mu][1][2][IM][xx]);

		  u_packed[0][block][mu][2][0][RE][xx]=rep<FT,double>(br*fr-bi*fi-er*cr+ei*ci);
		  u_packed[0][block][mu][2][0][IM][xx]=rep<FT,double>(er*ci+ei*cr-br*fi-bi*fr);
		  u_packed[0][block][mu][2][1][RE][xx]=rep<FT,double>(dr*cr-di*ci-ar*fr+ai*fi);
		  u_packed[0][block][mu][2][1][IM][xx]=rep<FT,double>(ar*fi+ai*fr-dr*ci-di*cr);
		  u_packed[0][block][mu][2][2][RE][xx]=rep<FT,double>(ar*er-ai*ei-dr*br+di*bi);
		  u_packed[0][block][mu][2][2][IM][xx]=rep<FT,double>(dr*bi+di*br-ar*ei-ai*er);
		}
	      }

	      {
		for(int x=0; x < S; x++) {
		  int xx = (y%nyg)*S+x;
		  // 3rd row reconstruction.
		  double ar=rep<double,FT>(u_packed[1][block][mu][0][0][RE][xx]);
		  double ai=rep<double,FT>(u_packed[1][block][mu][0][0][IM][xx]);

		  double br=rep<double,FT>(u_packed[1][block][mu][0][1][RE][xx]);
		  double bi=rep<double,FT>(u_packed[1][block][mu][0][1][IM][xx]);

		  double cr=rep<double,FT>(u_packed[1][block][mu][0][2][RE][xx]);
		  double ci=rep<double,FT>(u_packed[1][block][mu][0][2][IM][xx]);

		  double dr=rep<double,FT>(u_packed[1][block][mu][1][0][RE][xx]);
		  double di=rep<double,FT>(u_packed[1][block][mu][1][0][IM][xx]);

		  double er=rep<double,FT>(u_packed[1][block][mu][1][1][RE][xx]);
		  double ei=rep<double,FT>(u_packed[1][block][mu][1][1][IM][xx]);

		  double fr=rep<double,FT>(u_packed[1][block][mu][1][2][RE][xx]);
		  double fi=rep<double,FT>(u_packed[1][block][mu][1][2][IM][xx]);

		  u_packed[1][block][mu][2][0][RE][xx]=rep<FT,double>(br*fr-bi*fi-er*cr+ei*ci);
		  u_packed[1][block][mu][2][0][IM][xx]=rep<FT,double>(er*ci+ei*cr-br*fi-bi*fr);
		  u_packed[1][block][mu][2][1][RE][xx]=rep<FT,double>(dr*cr-di*ci-ar*fr+ai*fi);
		  u_packed[1][block][mu][2][1][IM][xx]=rep<FT,double>(ar*fi+ai*fr-dr*ci-di*cr);
		  u_packed[1][block][mu][2][2][RE][xx]=rep<FT,double>(ar*er-ai*ei-dr*br+di*bi);
		  u_packed[1][block][mu][2][2][IM][xx]=rep<FT,double>(dr*bi+di*br-ar*ei-ai*er);
		} // x
	      }

	    } // mu
	  } // s
	} // y
      } // z
    } // t

  } // end if ! compress

  double end = omp_get_wtime();
  masterPrintf(" %g sec\n", end - start);
#endif

  // Allocate data for the spinors
  Spinor* p_even=(Spinor*)geom.allocCBFourSpinor();
  Spinor* p_odd=(Spinor*)geom.allocCBFourSpinor();
  Spinor* c_even=(Spinor*)geom.allocCBFourSpinor();
  Spinor* c_odd=(Spinor*)geom.allocCBFourSpinor();
  Spinor* prep_p_even=(Spinor*)geom.allocCBFourSpinor();
  Spinor* prep_p_odd=(Spinor*)geom.allocCBFourSpinor();


  // Point to the second block of the array. Now there is padding on both ends.
  Spinor *psi_s[2] = { p_even, p_odd };
  Spinor *chi_s[2] = { c_even, c_odd };
  Spinor *prep_psi_s[2] = { prep_p_even, prep_p_odd };

//  spinor ** solver_field = NULL;
//  const int nr_sf = 2;
//  init_solver_field(&solver_field, VOLUMEPLUSRAND, nr_sf);

#if 1
// if(eo) convert to lexic
//  convert_lexic_to_eo(solver_field[0], solver_field[1], Q);

//  qdp_pack_spinor<>((const double**)solver_field, p_even, p_odd, geom);
  // need only odd sites
  qdp_pack_cb_spinor((const double*)Q,p_odd,geom,1);
#else
  masterPrintf("Filling Input spinor: ");


  start=omp_get_wtime();
#pragma omp parallel for collapse(4)
  for(int t=0; t < lT; t++) {
    for(int z=0; z < lZ; z++) {
      for(int y=0; y < lY; y++) {
	for(int s=0; s < nvecs; s++) {
	  for(int spin=0; spin < 4; spin++) {
	    for(int col=0; col < 3; col++)  {
	      for(int x=0; x < S; x++) {
		double d1=drand48()-0.5;
		double d2=drand48()-0.5;
		double d3=drand48()-0.5;
		double d4=drand48()-0.5;

		int ind = t*Pxyz+z*Pxy+y*nvecs+s; //((t*Nz+z)*Ny+y)*nvecs+s;
		psi_s[0][ind][col][spin][0][x] = rep<FT,double>(0/*d1*/);
		psi_s[0][ind][col][spin][1][x] = rep<FT,double>(0/*d2*/);
		psi_s[1][ind][col][spin][0][x] = rep<FT,double>(0/*d3*/);
		psi_s[1][ind][col][spin][1][x] = rep<FT,double>(0/*d4*/);

		if( ind==0 && col==0 && spin==0 && x==0 )
			psi_s[0][ind][col][spin][0][x] = rep<FT,double>(1);
	      }
	    }
	  }
	}
      }
    }
  }

  end = omp_get_wtime();
  masterPrintf(" %g sec\n", end - start);
#endif

  masterPrintf("Zeroing output spinor: ");
  start = omp_get_wtime();
#pragma omp parallel for collapse(4)
  for(int t=0; t < lT; t++) {
    for(int z=0; z < lZ; z++) {
      for(int y=0; y < lY; y++) {
	for(int s=0; s < nvecs; s++) {
	  for(int spin=0; spin < 4; spin++) {
	    for(int col=0; col < 3; col++)  {
	      for(int x=0; x < S; x++) {
		double d=0;
		int ind = t*Pxyz+z*Pxy+y*nvecs+s; //((t*Nz+z)*Ny+y)*nvecs+s;
		chi_s[0][ind][col][spin][0][x] = rep<FT,double>(d);
		chi_s[0][ind][col][spin][1][x] = rep<FT,double>(d);
		chi_s[1][ind][col][spin][0][x] = rep<FT,double>(d);
		chi_s[1][ind][col][spin][1][x] = rep<FT,double>(d);
	      }
	    }
	  }
	}
      }
    }
  }
  end = omp_get_wtime();
  masterPrintf(" %g sec\n", end -start);


#if 1
  masterPrintf("Creating Wilson Op\n");
  double Mass = 1.0/(2.0*g_kappa) - 4.0;
  masterPrintf("kappa=%f => Mass=%f\n",g_kappa,Mass);
  EvenOddWilsonOperator<FT, V, S,compress> M(Mass, u_packed, &geom, t_boundary, coeff_s, coeff_t);
#endif

#if 1
//  double rsd_target=rsdTarget<FT>::value;
//  int max_iters=5000;
  int niters;
  double rsd_final;
  int len = (geom.getPxyz()*geom.Nt()*sizeof(Spinor))/sizeof(FT);
  FT *c_s0 = (FT *)chi_s[0];

  {
    masterPrintf("Creating Solver\n");
    InvCG<FT,V,S, compress> solver(M, max_iter);

    masterPrintf("Tuning Solver\n");
    solver.tune();

    for(int solve = 0; solve < 1; solve++ ) {
      masterPrintf("Starting solver\n");
      unsigned long site_flops=0;
      unsigned long mv_apps=0;

//      FT *psi_0 = (FT *)psi_s[0];
//#if defined(__INTEL_COMPILER)
//#pragma simd
//#endif
//#pragma omp parallel for
//      for(int i=0; i < len; i++) {
//	c_s0[i] = rep<FT,double>(0);
//      psi_0[i]= rep<FT,double>(0);
//      }
//
//#pragma omp parallel for collapse(4)
//      for(int t=0; t < lT; t++) {
//	for(int z=0; z < lZ; z++) {
//	  for(int y=0; y < lY; y++) {
//	    for(int s=0; s < nvecs; s++) {
//	      for(int spin=0; spin < 4; spin++) {
//		for(int col=0; col < 3; col++)  {
//		  for(int x=0; x < S; x++) {
//
//		    int ind = t*Pxyz+z*Pxy+y*nvecs+s; //((t*Nz+z)*Ny+y)*nvecs+s;
//		    int x_coord = s*S + x;
//		    double d1 = drand48()-0.5;
//		    double d2 = drand48()-0.5;
//		    double d3 = drand48()-0.5;
//		    double d4 = drand48()-0.5;
//
//		    psi_s[0][ind][col][spin][0][x] = rep<FT,double>(d1);
//		    psi_s[0][ind][col][spin][1][x] = rep<FT,double>(d2);
//		    psi_s[1][ind][col][spin][0][x] = rep<FT,double>(d3);
//		    psi_s[1][ind][col][spin][1][x] = rep<FT,double>(d4);
//		  }
//		}
//	      }
//	    }
//	  }
//	}
//      }

#if 0
      // Prepare the source for preconditioning
      Spinor* tmp_even=(Spinor*)geom.allocCBFourSpinor();
      Spinor* tmp_odd=(Spinor*)geom.allocCBFourSpinor();
      Spinor *tmp_s[2] = { tmp_even, tmp_odd };

      // p_even -> M_{ee}^{-1} p_even
      p_even *= 1.0/(4.0+Mass); // TODO refactor

      // tmp_odd = M_{oe} p_even
      int target_cb = 1, source_cb = 0;
      M->dslash(tmp_s[target_cb], p_even, u_packed[target_cb], 1, target_cb);
      tmp_s[target_cb] *= -0.5;

      // this is now p_odd' = p_odd - M_{oe} M_{ee}^{-1} p_even
      psi_s[1] -= tmp_s[1];
      // i.e. we have psi' = L^{-1} psi

      // trivially solve for c_even'
      // c_even' = M_{ee}^{-1} p_even
      chi_s[0] = 1.0/(4.0+Mass) * psi_s[0]; // TODO refactor
#endif

      // prepare source for CG: psi' -> M^\dagger psi'
//      M(prep_psi_s[1], psi_s[1] , -1);

      // solve for c_odd'
      start = omp_get_wtime();
      solver(chi_s[1], psi_s[1], eps_sq, niters, rsd_final, site_flops, mv_apps,1,verbose);
      end = omp_get_wtime();


      unsigned long num_cb_sites=X1h*Ny*Nz*Nt;
      unsigned long total_flops = (site_flops + (72+2*1320)*mv_apps)*num_cb_sites;
      masterPrintf("Solver Time=%g(s)\n", (end-start));
      masterPrintf("CG GFLOPS=%g\n", 1.0e-9*(double)(total_flops)/(end -start));

#if 0
      // need to reconstruct solution (preconditioning)
      // tmp_even = M_{eo} c_odd
      target_cb = 0; source_cb = 1;
      M->dslash(tmp_s[target_cb], chi_s[source_cb], u_packed[target_cb], 1, target_cb);
      tmp_s[target_cb] *= -0.5;

      // tmp_even -> M_{ee}^{-1} tmp_even
      tmp_s[0] *= 1.0/(4.0+Mass); // TODO refactor

      chi_s[0] -= tmp_s[0];
#endif

      // unpack
//      qdp_unpack_spinor<>(chi_s[0], chi_s[1], (double*)P, geom);
      // only odd sites
      qdp_unpack_cb_spinor(chi_s[1],(double*)P,geom,1);

      //(const double**)solver_field, p_even, p_odd,
//      convert_eo_to_lexic(P, solver_field[0], solver_field[1]); //solver_field[0], solver_field[1], Q);
      //if(eo) convert from lexic to eo
    }
  } // Solver
#endif



  masterPrintf("Cleaning up\n");

//  finalize_solver(solver_field, nr_sf);

  geom.free(packed_gauge_cb0);
  geom.free(packed_gauge_cb1);
  geom.free(p_even);
  geom.free(p_odd);
  geom.free(c_even);
  geom.free(c_odd);
  geom.free(prep_p_even);
  geom.free(prep_p_odd);

}


void _initQphix(int argc, char **argv, int By_, int Bz_, int NCores_, int Sy_, int Sz_, int PadXY_, int PadXYZ_, int MinCt_, int c12, QphixPrec precision_)
{
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

  QPhiX::masterPrintf("Declared QMP Topology: %d %d %d %d\n",
		qmp_geom[0], qmp_geom[1], qmp_geom[2], qmp_geom[3]);


#ifdef QPHIX_QPX_SOURCE
  if( thread_bind ) {
    QPhiX::setThreadAffinity(NCores_user, Sy_user*Sz_user);
  }
  QPhiX::reportAffinity();
#endif

  QPhiX::masterPrintf("Launching TestCase\n");

  // Launch the test case.
//  timeDslashNoQDP test(By_user, Bz_user, NCores_user, Sy_user, Sz_user, PadXY_user, PadXYZ_user, MinCt_user,  iters, compress12, prec_user);

//  test.run(nrow_in, qmp_geometry);
//#ifdef QPHIX_QMP_COMMS
//  QMP_finalize_msg_passing();
//#endif



//  if( verbose > 0 )
//    printf("\nDetected QUDA version %d.%d.%d\n\n", QUDA_VERSION_MAJOR, QUDA_VERSION_MINOR, QUDA_VERSION_SUBMINOR);
//
////  error_root((QUDA_VERSION_MAJOR == 0 && QUDA_VERSION_MINOR < 7),1,"_initQphix [qphix_interface.c]","minimum QUDA version required is 0.7.0 (for support of chiral basis and removal of bug in mass normalization with preconditioning).");
//
//
//
//  gauge_param = newQudaGaugeParam();
//  inv_param = newQudaInvertParam();
//
//  // *** QUDA parameters begin here (may be modified)
//  QudaDslashType dslash_type = QUDA_CLOVER_WILSON_DSLASH;
//  QudaPrecision cpu_prec  = QUDA_DOUBLE_PRECISION;
//  QudaPrecision cuda_prec = QUDA_DOUBLE_PRECISION;
//  QudaPrecision cuda_prec_sloppy = QUDA_SINGLE_PRECISION;
//  QudaPrecision cuda_prec_precondition = QUDA_HALF_PRECISION;
//#if TRIVIAL_BC
//  gauge_param.t_boundary = ( abs(phase_0)>0.0 ? QUDA_ANTI_PERIODIC_T : QUDA_PERIODIC_T );
//  QudaReconstructType link_recon = 12;
//  QudaReconstructType link_recon_sloppy = 12;
//#else
//  gauge_param.t_boundary = QUDA_PERIODIC_T; // BC will be applied to gaugefield
//  QudaReconstructType link_recon = 18;
//  QudaReconstructType link_recon_sloppy = 18;
//#endif
//  QudaTune tune = QUDA_TUNE_YES;
//
//
//  // *** the remainder should not be changed for this application
//  // local lattice size
//#if USE_LZ_LY_LX_T
//  gauge_param.X[0] = LZ;
//  gauge_param.X[1] = LY;
//  gauge_param.X[2] = LX;
//  gauge_param.X[3] = T;
//#else
//  gauge_param.X[0] = LX;
//  gauge_param.X[1] = LY;
//  gauge_param.X[2] = LZ;
//  gauge_param.X[3] = T;
//#endif
//
//  inv_param.Ls = 1;
//
//  gauge_param.anisotropy = 1.0;
//  gauge_param.type = QUDA_WILSON_LINKS;
//  gauge_param.gauge_order = QUDA_QDP_GAUGE_ORDER;
//
//  gauge_param.cpu_prec = cpu_prec;
//  gauge_param.cuda_prec = cuda_prec;
//  gauge_param.reconstruct = link_recon;
//  gauge_param.cuda_prec_sloppy = cuda_prec_sloppy;
//  gauge_param.reconstruct_sloppy = link_recon_sloppy;
//  gauge_param.cuda_prec_precondition = cuda_prec_precondition;
//  gauge_param.reconstruct_precondition = link_recon_sloppy;
//  gauge_param.gauge_fix = QUDA_GAUGE_FIXED_NO;
//
//  inv_param.dslash_type = dslash_type;
//
//
//  // offsets used only by multi-shift solver
//  inv_param.num_offset = 4;
//  double offset[4] = {0.01, 0.02, 0.03, 0.04};
//  for (int i=0; i<inv_param.num_offset; i++) inv_param.offset[i] = offset[i];
//
//  inv_param.matpc_type = QUDA_MATPC_EVEN_EVEN;// QUDA_MATPC_EVEN_EVEN_ASYMMETRIC;
//  inv_param.solution_type = QUDA_MAT_SOLUTION;//QUDA_MAT_SOLUTION;
//
//  inv_param.dagger = QUDA_DAG_NO;
//  inv_param.mass_normalization = QUDA_KAPPA_NORMALIZATION;
//  inv_param.solver_normalization = QUDA_DEFAULT_NORMALIZATION;
//  inv_param.solve_type = QUDA_NORMOP_PC_SOLVE;
//  inv_param.inv_type = QUDA_CG_INVERTER;//QUDA_CG_INVERTER;
//
//  inv_param.pipeline = 0;
//  inv_param.gcrNkrylov = 10;
//
//  // require both L2 relative and heavy quark residual to determine convergence
//  inv_param.residual_type = (QudaResidualType)(QUDA_L2_RELATIVE_RESIDUAL | QUDA_HEAVY_QUARK_RESIDUAL);
//  inv_param.tol_hq = 1.0;//1e-3; // specify a tolerance for the residual for heavy quark residual
//  inv_param.reliable_delta = 1e-2; // ignored by multi-shift solver
//
//  // domain decomposition preconditioner parameters
//  inv_param.inv_type_precondition = QUDA_INVALID_INVERTER;
//  inv_param.schwarz_type = QUDA_ADDITIVE_SCHWARZ;
//  inv_param.precondition_cycle = 1;
//  inv_param.tol_precondition = 1e-1;
//  inv_param.maxiter_precondition = 10;
//  inv_param.verbosity_precondition = QUDA_SILENT;
//  inv_param.cuda_prec_precondition = cuda_prec_precondition;
//  inv_param.omega = 1.0;
//
//  inv_param.cpu_prec = cpu_prec;
//  inv_param.cuda_prec = cuda_prec;
//  inv_param.cuda_prec_sloppy = cuda_prec_sloppy;
//  inv_param.preserve_source = QUDA_PRESERVE_SOURCE_YES;
//  inv_param.gamma_basis = QUDA_CHIRAL_GAMMA_BASIS;
//  inv_param.dirac_order = QUDA_DIRAC_ORDER;
//
//  inv_param.input_location = QUDA_CPU_FIELD_LOCATION;
//  inv_param.output_location = QUDA_CPU_FIELD_LOCATION;
//
//  inv_param.tune = tune ? QUDA_TUNE_YES : QUDA_TUNE_NO;
//
//  gauge_param.ga_pad = 0; // 24*24*24/2;
//  inv_param.sp_pad = 0; // 24*24*24/2;
//  inv_param.cl_pad = 0; // 24*24*24/2;
//
//  // For multi-GPU, ga_pad must be large enough to store a time-slice
//  int x_face_size = gauge_param.X[1]*gauge_param.X[2]*gauge_param.X[3]/2;
//  int y_face_size = gauge_param.X[0]*gauge_param.X[2]*gauge_param.X[3]/2;
//  int z_face_size = gauge_param.X[0]*gauge_param.X[1]*gauge_param.X[3]/2;
//  int t_face_size = gauge_param.X[0]*gauge_param.X[1]*gauge_param.X[2]/2;
//  int pad_size =MAX(x_face_size, y_face_size);
//  pad_size = MAX(pad_size, z_face_size);
//  pad_size = MAX(pad_size, t_face_size);
//  gauge_param.ga_pad = pad_size;
//
//  if( verbose == 0 )
//    inv_param.verbosity = QUDA_SILENT;
//  else if( verbose == 1 )
//    inv_param.verbosity = QUDA_SUMMARIZE;
//  else
//    inv_param.verbosity = QUDA_VERBOSE;
//
//
//  // declare the grid mapping used for communications in a multi-GPU grid
//#if USE_LZ_LY_LX_T
//  int grid[4] = {NPROC3, NPROC2, NPROC1, NPROC0};
//#else
//  int grid[4] = {g_nproc_x, g_nproc_y, g_nproc_z, g_nproc_t};
//#endif
//
//  initCommsGridQuda(4, grid, commsMap, NULL);
//
//  // alloc gauge_qphix
//  size_t gSize = (gauge_param.cpu_prec == QUDA_DOUBLE_PRECISION) ? sizeof(double) : sizeof(float);
//
//  for (int dir = 0; dir < 4; dir++) {
//    gauge_qphix[dir] = (double*) malloc(VOLUME*18*gSize);
////    error_root((gauge_qphix[dir] == NULL),1,"_initQphix [qphix_interface.c]","malloc for gauge_qphix[dir] failed");
//  }
//
//  // alloc space for a temp. spinor, used throughout this module
//  tempSpinor  = (double*)malloc( VOLUME*24*sizeof(double) );
//
//  // only needed if even_odd_flag set
//  fullSpinor1 = (double*)malloc( VOLUME*24*sizeof(double) );
//  fullSpinor2 = (double*)malloc( VOLUME*24*sizeof(double) );
//
////  error_root((tempSpinor == NULL),1,"reorder_spinor_toQphix [qphix_interface.c]","malloc for tempSpinor failed");
//
//  // initialize the QUDA library
//  initQuda(-1);
}

// finalize the QUDA library
void _endQphix()
{
//  freeGaugeQuda();
//  free((void*)tempSpinor);
//  endQuda();
}


void _loadGaugeQphix()
{
//  if( inv_param.verbosity > QUDA_SILENT )
//    printf("\nCalled _loadGaugeQphix\n\n");
//
//  // update boundary if necessary
////   if( query_flags(UDBUF_UP2DATE)!=1 )
////    copy_bnd_ud();
//
//  _Complex double tmpcplx;
//
//  size_t gSize = (gauge_param.cpu_prec == QUDA_DOUBLE_PRECISION) ? sizeof(double) : sizeof(float);
//
//  // now copy and reorder
//  for( int x0=0; x0<T; x0++ )
//    for( int x1=0; x1<LX; x1++ )
//      for( int x2=0; x2<LY; x2++ )
//        for( int x3=0; x3<LZ; x3++ )
//        {
//          /* ipt[x3+LZ*x2+LY*LZ*x1+LX*LY*LZ*x0] is the index of the
//             point on the local lattice with cartesian coordinates
//             (x0,x1,x2,x3) */
//
//#if USE_LZ_LY_LX_T
//          int j = x3 + LZ*x2 + LY*LZ*x1 + LX*LY*LZ*x0;
//          int tm_idx   = g_ipt[x0][x1][x2][x3];
//#else
//          int j = x1 + LX*x2 + LY*LX*x3 + LZ*LY*LX*x0;
//          int tm_idx   = x3 + LZ*x2 + LY*LZ*x1 + LX*LY*LZ*x0;//g_ipt[x0][x3][x2][x1];
//#endif
//
//          int oddBit = (x0+x1+x2+x3) & 1;
//          int qphix_idx = 18*(oddBit*VOLUME/2+j/2);
//
//
//#if USE_LZ_LY_LX_T
//            memcpy( &(gauge_qphix[0][qphix_idx]), pud[tm_idx][3], 18*gSize);
//            memcpy( &(gauge_qphix[1][qphix_idx]), pud[tm_idx][2], 18*gSize);
//            memcpy( &(gauge_qphix[2][qphix_idx]), pud[tm_idx][1], 18*gSize);
//            memcpy( &(gauge_qphix[3][qphix_idx]), pud[tm_idx][0], 18*gSize);
//#else
//            memcpy( &(gauge_qphix[0][qphix_idx]), &(g_gauge_field[tm_idx][1]), 18*gSize);
//            memcpy( &(gauge_qphix[1][qphix_idx]), &(g_gauge_field[tm_idx][2]), 18*gSize);
//            memcpy( &(gauge_qphix[2][qphix_idx]), &(g_gauge_field[tm_idx][3]), 18*gSize);
//            memcpy( &(gauge_qphix[3][qphix_idx]), &(g_gauge_field[tm_idx][0]), 18*gSize);
//
//#if !(TRIVIAL_BC)
//            // apply boundary conditions
//            for( int i=0; i<9; i++ )
//            {
//            	tmpcplx = gauge_qphix[0][qphix_idx+2*i] + I*gauge_qphix[0][qphix_idx+2*i+1];
//            	tmpcplx *= -phase_1/g_kappa;
//            	gauge_qphix[0][qphix_idx+2*i]   = creal(tmpcplx);
//            	gauge_qphix[0][qphix_idx+2*i+1] = cimag(tmpcplx);
//
//            	tmpcplx = gauge_qphix[1][qphix_idx+2*i] + I*gauge_qphix[1][qphix_idx+2*i+1];
//            	tmpcplx *= -phase_2/g_kappa;
//            	gauge_qphix[1][qphix_idx+2*i]   = creal(tmpcplx);
//            	gauge_qphix[1][qphix_idx+2*i+1] = cimag(tmpcplx);
//
//            	tmpcplx = gauge_qphix[2][qphix_idx+2*i] + I*gauge_qphix[2][qphix_idx+2*i+1];
//            	tmpcplx *= -phase_3/g_kappa;
//            	gauge_qphix[2][qphix_idx+2*i]   = creal(tmpcplx);
//            	gauge_qphix[2][qphix_idx+2*i+1] = cimag(tmpcplx);
//
//            	tmpcplx = gauge_qphix[3][qphix_idx+2*i] + I*gauge_qphix[3][qphix_idx+2*i+1];
//            	tmpcplx *= -phase_0/g_kappa;
//            	gauge_qphix[3][qphix_idx+2*i]   = creal(tmpcplx);
//            	gauge_qphix[3][qphix_idx+2*i+1] = cimag(tmpcplx);
//            }
//#endif
//
//#endif
//        }
//
//  loadGaugeQuda((void*)gauge_qphix, &gauge_param);
}


// reorder spinor to QPhiX format (odd sites only)
void reorder_spinor_toQphix( double* sp )
{
  double startTime = gettime();
  // alloc space for a temp. spinor, used throughout this module TODO put somewhere central
  tempSpinor  = (double*)malloc( (VOLUME)*24*sizeof(double) );
  memcpy( tempSpinor, sp, (VOLUME/2)*24*sizeof(double) );
  double *in,*out;
  int Ns=4,Nc=3;
  double K1[4] = {-1.0,1.0,1.0,-1.0};
  int s1[4] = {3,2,1,0};

  // now copy and reorder from tempSpinor to spinor 
  for( int x0=0; x0<T; x0++ )
    for( int x1=0; x1<LX; x1++ )
      for( int x2=0; x2<LY; x2++ )
        for( int x3=0; x3<LZ; x3++ )
        {
#if USE_LZ_LY_LX_T
          int j = x3 + LZ*x2 + LY*LZ*x1 + LX*LY*LZ*x0;
          int tm_idx   = g_ipt[x0][x1][x2][x3];
#else
          int j = x3 + LZ*x2 + LY*LZ*x1 + LX*LY*LZ*x0;
          int tm_idx = g_lexic2eosub[ g_ipt[x0][x1][x2][x3] ];
#endif
          // is this an odd site?
          if((x0+x1+x2+x3+g_proc_coords[3]*LZ+g_proc_coords[2]*LY
          + g_proc_coords[0]*T+g_proc_coords[1]*LX)%2 != 0) {
           // gamma basis transformation
            in  = tempSpinor + 24*tm_idx;
            out = sp + 24*(j/2);

            for (int s=0; s<Ns; s++) {
              for (int c=0; c<Nc; c++) {
                for (int z=0; z<2; z++) {
                  out[(s*Nc+c)*2+z] = K1[s]*in[(s1[s]*Nc+c)*2+z];
                }
              }
            }
          }

//          int oddBit = (x0+x1+x2+x3) & 1;
//          int qphix_idx = 18*(oddBit*VOLUME/2+j/2);
//          memcpy( &(spinor[24*(oddBit*VOLUME/2+j/2)]), &(tempSpinor[24*tm_idx]), 24*sizeof(double));

        }
        
  double endTime = gettime();
  double diffTime = endTime - startTime;
  printf("time spent in reorder_spinor_toQphix: %f secs\n", diffTime);
}

// reorder spinor from QPhiX format (odd sites only)
void reorder_spinor_fromQphix( double* sp, double normFac=1.0 )
{
  double startTime = gettime();
  memcpy( tempSpinor, sp, (VOLUME/2)*24*sizeof(double) );
  double *in,*out;
  int Ns=4,Nc=3;
  double K1[4] = {-1.0,1.0,1.0,-1.0};
  int s1[4] = {3,2,1,0};

  // now copy and reorder from tempSpinor to spinor 
  for( int x0=0; x0<T; x0++ )
    for( int x1=0; x1<LX; x1++ )
      for( int x2=0; x2<LY; x2++ )
        for( int x3=0; x3<LZ; x3++ )
        {
#if USE_LZ_LY_LX_T
          int j = x3 + LZ*x2 + LY*LZ*x1 + LX*LY*LZ*x0;
          int tm_idx   = g_ipt[x0][x1][x2][x3];
#else
          int j = x3 + LZ*x2 + LY*LZ*x1 + LX*LY*LZ*x0;
          int tm_idx = g_lexic2eosub[ g_ipt[x0][x1][x2][x3] ];
#endif
          // is this an odd site?
          if((x0+x1+x2+x3+g_proc_coords[3]*LZ+g_proc_coords[2]*LY
          + g_proc_coords[0]*T+g_proc_coords[1]*LX)%2 != 0) {
           // gamma basis transformation
            in  = tempSpinor + 24*(j/2);
            out = sp + 24*tm_idx;

            for (int s=0; s<Ns; s++) {
              for (int c=0; c<Nc; c++) {
                for (int z=0; z<2; z++) {
                  out[(s*Nc+c)*2+z] = K1[s]*in[(s1[s]*Nc+c)*2+z] * normFac;
                }
              }
            }
          }

//          int oddBit = (x0+x1+x2+x3) & 1;
//          int qphix_idx = 18*(oddBit*VOLUME/2+j/2);
//          memcpy( &(spinor[24*(oddBit*VOLUME/2+j/2)]), &(tempSpinor[24*tm_idx]), 24*sizeof(double));

        }
  double endTime = gettime();
  double diffTime = endTime - startTime;
  printf("time spent in reorder_spinor_fromQphix: %f secs\n", diffTime);
}

// if even_odd_flag set
void M_full_qphix(spinor * const Even_new, spinor * const Odd_new,  spinor * const Even, spinor * const Odd)
{
//  inv_param.kappa = g_kappa;
//  inv_param.mu = fabs(g_mu);
//  inv_param.epsilon = 0.0;
//
//  // IMPORTANT: use opposite TM flavor since gamma5 -> -gamma5 (until LXLYLZT prob. resolved)
//  inv_param.twist_flavor = (g_mu < 0.0 ? QUDA_TWIST_PLUS : QUDA_TWIST_MINUS);
//  inv_param.Ls = (inv_param.twist_flavor == QUDA_TWIST_NONDEG_DOUBLET ||
//       inv_param.twist_flavor == QUDA_TWIST_DEG_DOUBLET ) ? 2 : 1;
//
//  void *spinorIn  = (void*)fullSpinor1;
//  void *spinorOut = (void*)fullSpinor2;
//
//  // reorder spinor
//  convert_eo_to_lexic( spinorIn, Even, Odd );
//  reorder_spinor_toQphix( (double*)spinorIn, inv_param.cpu_prec );
//
//  // multiply
////   inv_param.solution_type = QUDA_MAT_SOLUTION;
//  MatQuda( spinorOut, spinorIn, &inv_param);
//
//  // reorder spinor
////  reorder_spinor_fromQphix( (double*)spinorIn,  inv_param.cpu_prec );
////  convert_lexic_to_eo( Even, Odd, spinorIn );
//
//  reorder_spinor_fromQphix( (double*)spinorOut, inv_param.cpu_prec );
//  convert_lexic_to_eo( Even_new, Odd_new, spinorOut );
}


// no even-odd
void D_psi_qphix(spinor * const P, spinor * const Q)
{
//  inv_param.kappa = g_kappa;
//  inv_param.mu = fabs(g_mu);
//  inv_param.epsilon = 0.0;
//
//  // IMPORTANT: use opposite TM flavor since gamma5 -> -gamma5 (until LXLYLZT prob. resolved)
//  inv_param.twist_flavor = (g_mu < 0.0 ? QUDA_TWIST_PLUS : QUDA_TWIST_MINUS);
//  inv_param.Ls = (inv_param.twist_flavor == QUDA_TWIST_NONDEG_DOUBLET ||
//       inv_param.twist_flavor == QUDA_TWIST_DEG_DOUBLET ) ? 2 : 1;
//
//  void *spinorIn  = (void*)Q;
//  void *spinorOut = (void*)P;
//
//  // reorder spinor
//  reorder_spinor_toQphix( (double*)spinorIn, inv_param.cpu_prec );
//
//  // multiply
////   inv_param.solution_type = QUDA_MAT_SOLUTION;
//  MatQuda( spinorOut, spinorIn, &inv_param);
//
//  // reorder spinor
//  reorder_spinor_fromQphix( (double*)spinorIn,  inv_param.cpu_prec );
//  reorder_spinor_fromQphix( (double*)spinorOut, inv_param.cpu_prec );
}


int invert_qphix(spinor * const P, spinor * const Q, const int max_iter, double eps_sq, const int rel_prec )
{
  // lattice size
  subLattSize[0] = LX;
  subLattSize[1] = LY;
  subLattSize[2] = LZ;
  subLattSize[3] = T;

  // Work out global lattice size
  lattSize[0] = LX*g_nproc_x;
  lattSize[1] = LY*g_nproc_y;
  lattSize[2] = LZ*g_nproc_z;
  lattSize[3] = T*g_nproc_t;

  // reorder spinor
  reorder_spinor_toQphix( (double*)Q );

  if ( precision == FLOAT_PREC ) {
    if ( QPHIX_SOALEN > VECLEN_SP ) {
      masterPrintf("SOALEN=%d is greater than the single prec VECLEN=%d\n", QPHIX_SOALEN,VECLEN_SP);
      abort();
    }

    masterPrintf("TIMING IN SINGLE PRECISION \n");
    if ( compress12 ) {
      invert<float,VECLEN_SP,QPHIX_SOALEN,true>(P,Q,max_iter,eps_sq,rel_prec);
    }
    else {
      invert<float,VECLEN_SP,QPHIX_SOALEN,false>(P,Q,max_iter,eps_sq,rel_prec);
    }
  }

#if 1
#if defined(QPHIX_MIC_SOURCE)
  if ( precision == HALF_PREC ) {
    if ( QPHIX_SOALEN > VECLEN_HP ) {
      masterPrintf("SOALEN=%d is greater than the single prec VECLEN=%d\n", QPHIX_SOALEN,VECLEN_SP);
      abort();
    }


    masterPrintf("TIMING IN HALF PRECISION \n");
    if ( compress12 ) {
      invert<half,VECLEN_HP,QPHIX_SOALEN,true>(P,Q,max_iter,eps_sq,rel_prec);
    }
    else {
      invert<half,VECLEN_HP,QPHIX_SOALEN,false>(P,Q,max_iter,eps_sq,rel_prec);
    }
  }
#endif

  if( precision == DOUBLE_PREC ) {
    if ( QPHIX_SOALEN > VECLEN_DP ) {
      masterPrintf("SOALEN=%d is greater than the double prec VECLEN=%d\n", QPHIX_SOALEN, VECLEN_DP);
      abort();
    }

    masterPrintf("TIMING IN DOUBLE PRECISION \n");
    if ( compress12 ) {
      invert<double,VECLEN_DP,QPHIX_SOALEN,true>(P,Q,max_iter,eps_sq,rel_prec);
    }
    else {
      invert<double,VECLEN_DP,QPHIX_SOALEN,false>(P,Q,max_iter,eps_sq,rel_prec);
    }
  }
#endif



//  double startTime = gettime();
//
//  if( inv_param.verbosity > QUDA_SILENT )
//    printf("\nCalled invert_qphix\n\n");
//
//  inv_param.maxiter = max_iter;
//
//  void *spinorIn  = (void*)Q; // source
//  void *spinorOut = (void*)P; // solution
//
//  //ranlxd(spinorIn,VOLUME*24); // a random source for debugging
//
//  // get initial rel. residual (rn) and residual^2 (nrm2) from DD
//  int iteration;
////  rn = getResidualDD(k,l,&nrm2);
////  double tol = res*rn;
//  inv_param.tol = sqrt(rel_prec);
//
//  // these can be set individually
//  for (int i=0; i<inv_param.num_offset; i++) {
//    inv_param.tol_offset[i] = inv_param.tol;
//    inv_param.tol_hq_offset[i] = inv_param.tol_hq;
//  }
//  inv_param.kappa = g_kappa;
//  inv_param.mu = fabs(g_mu);
//  inv_param.epsilon = 0.0;
//  // IMPORTANT: use opposite TM flavor since gamma5 -> -gamma5 (until LXLYLZT prob. resolved)
////  inv_param.twist_flavor = (g_mu < 0.0 ? QUDA_TWIST_PLUS : QUDA_TWIST_MINUS);
//  inv_param.Ls = (inv_param.twist_flavor == QUDA_TWIST_NONDEG_DOUBLET ||
//       inv_param.twist_flavor == QUDA_TWIST_DEG_DOUBLET ) ? 2 : 1;
//
//  if (inv_param.dslash_type == QUDA_CLOVER_WILSON_DSLASH || inv_param.dslash_type == QUDA_TWISTED_CLOVER_DSLASH) {
//    inv_param.clover_cpu_prec = QUDA_DOUBLE_PRECISION;
//    inv_param.clover_cuda_prec = QUDA_DOUBLE_PRECISION;
//    inv_param.clover_cuda_prec_sloppy = QUDA_SINGLE_PRECISION;
//    inv_param.clover_cuda_prec_precondition = QUDA_HALF_PRECISION;
//    inv_param.clover_order = QUDA_PACKED_CLOVER_ORDER;
//    inv_param.clover_coeff = g_c_sw*g_kappa;
//
//    loadCloverQuda(NULL, NULL, &inv_param);
//    printf("kappa = %f, c_sw = %f\n", inv_param.kappa, inv_param.clover_coeff);
//  }
//
//  // reorder spinor
//  reorder_spinor_toQphix( (double*)spinorIn, inv_param.cpu_prec );
//
//  // perform the inversion
//  invertQuda(spinorOut, spinorIn, &inv_param);
//
//  if( inv_param.verbosity == QUDA_VERBOSE )
//    printf("Device memory used:\n   Spinor: %f GiB\n    Gauge: %f GiB\n",
//	 inv_param.spinorGiB, gauge_param.gaugeGiB);
//	if( inv_param.verbosity > QUDA_SILENT )
//    printf("Done: %i iter / %g secs = %g Gflops\n",
//	 inv_param.iter, inv_param.secs, inv_param.gflops/inv_param.secs);
//
//  // number of CG iterations
//  iteration = inv_param.iter;
//
  // reorder spinor
  reorder_spinor_fromQphix( (double*)P, 1.0/g_kappa );
  reorder_spinor_fromQphix( (double*)Q );
//
//  double endTime = gettime();
//  double diffTime = endTime - startTime;
//  printf("time spent in tmcgne_qphix: %f secs\n", diffTime);
//
//  if(iteration >= max_iter) return(-1);
//  	  return(iteration);


	return 0;
}


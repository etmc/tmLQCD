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
*     and first call of the inverter. In particular, 'boundary(const double
*     kappa)' must be called before if nontrivial boundary conditions are to 
*     be used since those will be applied directly to the gaugefield.
*
*   double tmcgne_qphix(int nmx,double res,int k,int l,int *status,int *ifail)
*     The same functionality as 'tmcgne' (see tmcg.c) but inversion is performed
*     on the GPU using QUDA. Final residuum check is performed on the host (CPU)
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
#include "config.h"
#include "global.h"
extern "C" {
#include "boundary.h"
#include "linalg/convert_eo_to_lexic.h"
#include "solver/solver.h"
//#include "solver/solver_field.h"
#include "gettime.h"
#include "update_backward_gauge.h"
}
#ifdef TM_USE_OMP
#include <omp.h>
#endif
#include <cfloat>
#include <cstdlib>
#include <cstring>
#include "qphix/invbicgstab.h"
#include "qphix/invcg.h"
#include "qphix/print_utils.h"
#include "qphix/qphix_config.h"
#include "qphix/wilson.h"
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

#if defined(QPHIX_AVX2_SOURCE)
#define VECLEN_SP 8
#define VECLEN_DP 4
#endif

#if defined(QPHIX_AVX_SOURCE)
#define VECLEN_SP 8
#define VECLEN_DP 4
#endif

#if defined(QPHIX_SSE_SOURCE)
#define VECLEN_SP 4
#define VECLEN_DP 2
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
int qmp_geom[4] = {1, 1, 1, 1};

template <typename T>
struct rsdTarget {
  static const double value;
};

template <>
const double rsdTarget<half>::value = (double)(1.0e-4);

template <>
const double rsdTarget<float>::value = (double)(1.0e-7);

template <>
const double rsdTarget<double>::value = (double)(1.0e-12);



void _initQphix(int argc, char **argv, int By_, int Bz_, int NCores_, int Sy_, int Sz_, int PadXY_,
                int PadXYZ_, int MinCt_, int c12, QphixPrec precision_) {
  // Global Lattice Size
  lattSize[0] = LX * g_nproc_x;
  lattSize[1] = LY * g_nproc_y;
  lattSize[2] = LZ * g_nproc_z;
  lattSize[3] = T * g_nproc_t;

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
  N_simt = Sy_ * Sz_;
  compress12 = c12;
  precision = precision_;

  omp_set_num_threads(NCores * Sy * Sz);

#ifdef QPHIX_QMP_COMMS
  // Declare the logical topology
  qmp_geom[0] = g_nproc_x;
  qmp_geom[1] = g_nproc_y;
  qmp_geom[2] = g_nproc_z;
  qmp_geom[3] = g_nproc_t;
  if (QMP_declare_logical_topology(qmp_geom, 4) != QMP_SUCCESS) {
    QMP_error("Failed to declare QMP Logical Topology\n");
    abort();
  }
#endif

#ifdef QPHIX_QPX_SOURCE
  if (thread_bind) {
    QPhiX::setThreadAffinity(NCores_user, Sy_user * Sz_user);
  }
  QPhiX::reportAffinity();
#endif
}


// Finalize the QPhiX library
void _endQphix() {}


// Reorder the tmLQCD gauge field to a cb0 and a cb1 QPhiX gauge field
template <typename FT, int VECLEN, int SOALEN, bool compress12>
void reorder_gauge_toQphix(Geometry<FT, VECLEN, SOALEN, compress12> &geom, double *qphix_gauge_cb0,
                           double *qphix_gauge_cb1) {
  double startTime = gettime();
  double *in;
  double *out;

  // Number of elements in spin, color & complex
  // Here c1 is QPhiX's outer color, and c2 the inner one
  int Ns = 4;
  int Nc1 = compress12 ? 2 : 3;
  int Nc2 = 3;
  int Nz = 2;

  // Geometric parameters for QPhiX data layout
  auto nyg = geom.nGY();
  auto nVecs = geom.nVecs();
  auto Pxy = geom.getPxy();
  auto Pxyz = geom.getPxyz();

  // This is needed to translate between the different
  // orderings of "\mu" in tmlQCD and QPhiX, respectively
  int change_dim[4] = {3, 0, 1, 2};

  // Get the base pointer for the (global) tmlQCD gauge field
  uint64_t tm_idx = 0;
  in = (double *)&g_gauge_field[0][0].c00;

  // This will loop over the entire lattice and calculate
  // the array and internal indices for both tmlQCD & QPhiX
  for (uint64_t t = 0; t < T; t++)
    for (uint64_t x = 0; x < LX; x++)
      for (uint64_t y = 0; y < LY; y++)
        for (uint64_t z = 0; z < LZ; z++) {
          // These are the QPhiX SIMD vector in checkerboarded x direction
          // (up to LX/2), the index inside one single Structure of Arrays (SOA)
          // and the internal position inside the ("packed") SIMD vector
          uint64_t SIMD_vector = (x / 2) / SOALEN;
          uint64_t x_one_SOA = (x / 2) % SOALEN;
          uint64_t x_internal = (y % nyg) * SOALEN + x_one_SOA;

          // Calculate the array index in QPhiX, given a global
          // lattice index (t,x,y,z). This is also called the
          // "block" in QPhiX and the QPhiX/QDP packers.
          uint64_t qphix_idx = (t * Pxyz + z * Pxy) / nyg + (y / nyg) * nVecs + SIMD_vector;

          if ((t + x + y + z) & 1)
            out = qphix_gauge_cb1;  // cb1
          else
            out = qphix_gauge_cb0;  // cb0

          for (int dim = 0; dim < 4; dim++)  // dimension == tmLQCD \mu
          {
            for (int dir = 0; dir < 2; dir++)  // direction == backward/forward
            {
              if (dir == 0) tm_idx = g_idn[g_ipt[t][x][y][z]][dim];  // this is the adjoint
                                                                     // gauge field to be
                                                                     // (backwards shift)
              else
                tm_idx = g_ipt[t][x][y][z];  // this is the normal gauge field
                                             // to be (same lattice site)

              for (int c1 = 0; c1 < Nc1; c1++)    // QPhiX convention color 1 (runs up to 2 or 3)
                for (int c2 = 0; c2 < Nc2; c2++)  // QPhiX convention color 2 (always runs up to 3)
                  for (int z = 0; z < Nz; z++) {
                    // Note:
                    // -----
                    // 1. \mu in QPhiX runs from 0..7 for all eight neighbouring
                    // links.
                    //    Here, the ordering of the direction (backward/forward)
                    //    is the same
                    //    for tmlQCD and QPhiX, but we have to change the
                    //    ordering of the dimensions.
                    int q_mu = 2 * change_dim[dim] + dir;

                    // 2. QPhiX gauge field matrices are transposed w.r.t.
                    // tmLQCD.
                    // 3. tmlQCD always uses 3x3 color matrices (Nc2*Nc2).
                    uint64_t t_inner_idx = z + c1 * Nz + c2 * Nz * Nc2 + dim * Nz * Nc2 * Nc2 +
                                           tm_idx * Nz * Nc2 * Nc2 * 4;
                    uint64_t q_inner_idx = x_internal + z * VECLEN + c2 * VECLEN * Nz +
                                           c1 * VECLEN * Nz * Nc2 + q_mu * VECLEN * Nz * Nc2 * Nc1 +
                                           qphix_idx * VECLEN * Nz * Nc2 * Nc1 * 8;

                    out[q_inner_idx] = in[t_inner_idx];
                  }

            }  // direction
          }    // dimension

        }  // volume

  double endTime = gettime();
  double diffTime = endTime - startTime;
  masterPrintf("  time spent in reorder_gauge_toQphix: %f secs\n", diffTime);
}

// Reorder tmLQCD spinor to a cb0 and cb1 QPhiX spinor
template <typename FT, int VECLEN, int SOALEN, bool compress12>
void reorder_spinor_toQphix(Geometry<FT, VECLEN, SOALEN, compress12> &geom, double *tm_spinor,
                            double *qphix_spinor_cb0, double *qphix_spinor_cb1) {
  double startTime = gettime();
  double *in;
  double *out;

  // Number of elements in spin, color & complex
  int Ns = 4;
  int Nc = 3;
  int Nz = 2;

  // Geometric parameters for QPhiX data layout
  auto nVecs = geom.nVecs();
  auto Pxy = geom.getPxy();
  auto Pxyz = geom.getPxyz();

  // This is needed to translate between the different
  // gamma bases tmlQCD and QPhiX are using
  int change_sign[4] = {1, -1, -1, 1};
  int change_spin[4] = {3, 2, 1, 0};

  // This will loop over the entire lattice and calculate
  // the array and internal indices for both tmlQCD & QPhiX
  for (uint64_t t = 0; t < T; t++)
    for (uint64_t x = 0; x < LX; x++)
      for (uint64_t y = 0; y < LY; y++)
        for (uint64_t z = 0; z < LZ; z++) {
          // These are the QPhiX SIMD vector in checkerboarded x direction
          // (up to LX/2) and the internal position inside the SIMD vector
          uint64_t SIMD_vector = (x / 2) / SOALEN;
          uint64_t x_internal = (x / 2) % SOALEN;

          // Calculate the array index in tmlQCD & QPhiX,
          // given a global lattice index (t,x,y,z)
          uint64_t qphix_idx = t * Pxyz + z * Pxy + y * nVecs + SIMD_vector;
          uint64_t tm_idx = g_ipt[t][x][y][z];

          // Calculate base point for every spinor field element (tmlQCD) or
          // for every SIMD vector of spinors, a.k.a FourSpinorBlock (QPhiX),
          // which will depend on the checkerboard (cb)
          in = tm_spinor + Ns * Nc * Nz * tm_idx;
          if ((t + x + y + z) & 1)
            out = qphix_spinor_cb1 + SOALEN * Nz * Nc * Ns * qphix_idx;  // cb1
          else
            out = qphix_spinor_cb0 + SOALEN * Nz * Nc * Ns * qphix_idx;  // cb0

          // Copy the internal elements, performing a gamma basis transformation
          for (int spin = 0; spin < Ns; spin++)  // QPhiX spin index
            for (int color = 0; color < Nc; color++)
              for (int z = 0; z < Nz; z++)  // RE or IM
              {
                uint64_t qId =
                    x_internal + z * SOALEN + spin * SOALEN * Nz + color * SOALEN * Nz * Ns;
                uint64_t tId = z + color * Nz + change_spin[spin] * Nz * Nc;

                out[qId] = change_sign[spin] * in[tId];
              }

        }  // volume

  double endTime = gettime();
  double diffTime = endTime - startTime;
  masterPrintf("  time spent in reorder_spinor_toQphix: %f secs\n", diffTime);
}

// Reorder a cb0 and cb1 QPhiX spinor to a tmLQCD spinor
template <typename FT, int VECLEN, int SOALEN, bool compress12>
void reorder_spinor_fromQphix(Geometry<FT, VECLEN, SOALEN, compress12> &geom, double *tm_spinor,
                              double *qphix_spinor_cb0, double *qphix_spinor_cb1,
                              double normFac = 1.0) {
  double startTime = gettime();
  double *in;
  double *out;

  // Number of elements in spin, color & complex
  int Ns = 4;
  int Nc = 3;
  int Nz = 2;

  // Geometric parameters for QPhiX data layout
  auto nVecs = geom.nVecs();
  auto Pxy = geom.getPxy();
  auto Pxyz = geom.getPxyz();

  // This is needed to translate between the different
  // gamma bases tmlQCD and QPhiX are using
  int change_sign[4] = {1, -1, -1, 1};
  int change_spin[4] = {3, 2, 1, 0};

  // This will loop over the entire lattice and calculate
  // the array and internal indices for both tmlQCD & QPhiX
  for (uint64_t t = 0; t < T; t++)
    for (uint64_t x = 0; x < LX; x++)
      for (uint64_t y = 0; y < LY; y++)
        for (uint64_t z = 0; z < LZ; z++) {
          // These are the QPhiX SIMD vector in checkerboarded x direction
          // (up to LX/2) and the internal position inside the SIMD vector
          uint64_t SIMD_vector = (x / 2) / SOALEN;
          uint64_t x_internal = (x / 2) % SOALEN;

          // Calculate the array index in tmlQCD & QPhiX,
          // given a global lattice index (t,x,y,z)
          uint64_t qphix_idx = t * Pxyz + z * Pxy + y * nVecs + SIMD_vector;
          uint64_t tm_idx = g_ipt[t][x][y][z];

          // Calculate base point for every spinor field element (tmlQCD) or
          // for every SIMD vector of spinors, a.k.a FourSpinorBlock (QPhiX),
          // which will depend on the checkerboard (cb)
          if ((t + x + y + z) & 1)
            in = qphix_spinor_cb1 + SOALEN * Nz * Nc * Ns * qphix_idx;  // cb1
          else
            in = qphix_spinor_cb0 + SOALEN * Nz * Nc * Ns * qphix_idx;  // cb0
          out = tm_spinor + Ns * Nc * Nz * tm_idx;

          // Copy the internal elements, performing a gamma basis transformation
          for (int spin = 0; spin < Ns; spin++)  // tmlQCD spin index
            for (int color = 0; color < Nc; color++)
              for (int z = 0; z < Nz; z++)  // RE or IM
              {
                uint64_t qId = x_internal + z * SOALEN + change_spin[spin] * SOALEN * Nz +
                               color * SOALEN * Nz * Ns;
                uint64_t tId = z + color * Nz + spin * Nz * Nc;

                out[tId] = normFac * change_sign[spin] * in[qId];
              }

        }  // volume

  double endTime = gettime();
  double diffTime = endTime - startTime;
  masterPrintf("  time spent in reorder_spinor_fromQphix: %f secs\n", diffTime);
}

template <typename FT, int V, int S, bool compress>
void D_psi(spinor* tmlqcd_out, const spinor* tmlqcd_in) {

  typedef typename Geometry<FT, V, S, compress>::SU3MatrixBlock QGauge;
  typedef typename Geometry<FT, V, S, compress>::FourSpinorBlock QSpinor;

  /************************
   *                      *
   *     Diagnostic       *
   *     Information      *
   *         &            *
   *    Creating Dslash   *
   *                      *
  ************************/

  masterPrintf("VECLEN=%d SOALEN=%d\n", V, S);
  masterPrintf("# Declared QMP Topology: %d %d %d %d\n", qmp_geom[0], qmp_geom[1], qmp_geom[2],
               qmp_geom[3]);
  masterPrintf("Global Lattice Size = ");
  for (int mu = 0; mu < 4; mu++) {
    masterPrintf(" %d", lattSize[mu]);
  }
  masterPrintf("\n");
  masterPrintf("Local Lattice Size = ");
  for (int mu = 0; mu < 4; mu++) {
    masterPrintf(" %d", subLattSize[mu]);
  }
  masterPrintf("\n");
  masterPrintf("Block Sizes: By= %d Bz=%d\n", By, Bz);
  masterPrintf("Cores = %d\n", NCores);
  masterPrintf("SMT Grid: Sy=%d Sz=%d\n", Sy, Sz);
  masterPrintf("Pad Factors: PadXY=%d PadXYZ=%d\n", PadXY, PadXYZ);
  masterPrintf("Threads_per_core = %d\n", N_simt);
  masterPrintf("MinCt = %d\n", MinCt);
  masterPrintf("Initializing QPhiX Dslash\n");

  // Create Dslash Class
  double t_boundary = (FT)(1);
  double coeff_s = (FT)(1);
  double coeff_t = (FT)(1);
  Geometry<FT, V, S, compress> geom(subLattSize, By, Bz, NCores, Sy, Sz, PadXY, PadXYZ, MinCt);
  Dslash<FT, V, S, compress> DQPhiX(&geom, t_boundary, coeff_s, coeff_t);

  /************************
   *                      *
   *     GAUGE FIELDS     *
   *                      *
  ************************/

  // Allocate data for the gauge fields
  QGauge *u_packed[2];
  QGauge *packed_gauge_cb0 =
      (QGauge *)geom.allocCBGauge();  // Links emanating from ODD sites (cb=0)
  QGauge *packed_gauge_cb1 =
      (QGauge *)geom.allocCBGauge();  // Links emanating from EVEN sites (cb=1)
  u_packed[0] = packed_gauge_cb0;
  u_packed[1] = packed_gauge_cb1;

  // Reorder (global) input gauge field from tmLQCD to QPhiX
  reorder_gauge_toQphix(geom, (double *)u_packed[0],
                        (double *)u_packed[1]);  // uses global tmlQCD gauge field as input

  /************************
   *                      *
   *     SPINOR FIELDS    *
   *                      *
  ************************/

  // Allocate data for the even/odd (checkerboarded) QPhiX spinors
  QSpinor *qphix_in[2];
  QSpinor *qphix_out[2];
  QSpinor *packed_spinor_in_cb0 = (QSpinor *)geom.allocCBFourSpinor();
  QSpinor *packed_spinor_in_cb1 = (QSpinor *)geom.allocCBFourSpinor();
  QSpinor *packed_spinor_out_cb0 = (QSpinor *)geom.allocCBFourSpinor();
  QSpinor *packed_spinor_out_cb1 = (QSpinor *)geom.allocCBFourSpinor();
  qphix_in[0] = packed_spinor_in_cb0;
  qphix_in[1] = packed_spinor_in_cb1;
  qphix_out[0] = packed_spinor_out_cb0;
  qphix_out[1] = packed_spinor_out_cb1;

  // Reorder input spinor from tmLQCD to QPhiX
  reorder_spinor_toQphix(geom, (double *)tmlqcd_in, (double *)qphix_in[0], (double *)qphix_in[1]);

  // Apply QPhiX Dslash to qphix_in spinors
  DQPhiX.dslash(qphix_out[1], qphix_in[0], u_packed[1],
                /* isign == non-conjugate */ 1, /* cb == */
                1);
  DQPhiX.dslash(qphix_out[0], qphix_in[1], u_packed[0],
                /* isign == non-conjugate */ 1, /* cb == */
                0);

  // Reorder spinor fields back to tmLQCD
  reorder_spinor_fromQphix(geom, (double *)tmlqcd_out, (double *)qphix_out[0],
                           (double *)qphix_out[1], (1. * g_kappa));
  reorder_spinor_fromQphix(geom, (double *)tmlqcd_in, (double *)qphix_in[0], (double *)qphix_in[1]);

  masterPrintf("Cleaning up\n");
  geom.free(packed_gauge_cb0);
  geom.free(packed_gauge_cb1);
  geom.free(packed_spinor_in_cb0);
  geom.free(packed_spinor_in_cb1);
  geom.free(packed_spinor_out_cb0);
  geom.free(packed_spinor_out_cb1);
}

void D_psi_qphix(spinor* tmlqcd_out, const spinor* tmlqcd_in) {
  if (precision == QPHIX_DOUBLE_PREC) {
    if (QPHIX_SOALEN > VECLEN_DP) {
      masterPrintf("SOALEN=%d is greater than the double prec VECLEN=%d\n", QPHIX_SOALEN,
                   VECLEN_DP);
      abort();
    }
    masterPrintf("TIMING IN DOUBLE PRECISION \n");
    if (compress12) {
      D_psi<double, VECLEN_DP, QPHIX_SOALEN, true>(tmlqcd_out, tmlqcd_in);
    } else {
      D_psi<double, VECLEN_DP, QPHIX_SOALEN, false>(tmlqcd_out, tmlqcd_in);
    }
  } else if (precision == QPHIX_FLOAT_PREC) {
    if (QPHIX_SOALEN > VECLEN_SP) {
      masterPrintf("SOALEN=%d is greater than the single prec VECLEN=%d\n", QPHIX_SOALEN,
                   VECLEN_SP);
      abort();
    }
    masterPrintf("TIMING IN SINGLE PRECISION \n");
    if (compress12) {
      D_psi<float, VECLEN_SP, QPHIX_SOALEN, true>(tmlqcd_out, tmlqcd_in);
    } else {
      D_psi<float, VECLEN_SP, QPHIX_SOALEN, false>(tmlqcd_out, tmlqcd_in);
    }
  }
#if defined(QPHIX_MIC_SOURCE)
  else if (precision == QPHIX_HALF_PREC) {
    if (QPHIX_SOALEN > VECLEN_HP) {
      masterPrintf("SOALEN=%d is greater than the single prec VECLEN=%d\n", QPHIX_SOALEN,
                   VECLEN_SP);
      abort();
    }
    masterPrintf("TIMING IN HALF PRECISION \n");
    if (compress12) {
      D_psi<half, VECLEN_HP, QPHIX_SOALEN, true>(tmlqcd_out, tmlqcd_in);
    } else {
      D_psi<half, VECLEN_HP, QPHIX_SOALEN, false>(tmlqcd_out, tmlqcd_in);
    }
  }
#endif
}

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
* Author: Mario Schroeck <mario.schroeck@roma3.infn.it>,
*         Peter Labus <Peter.Labus@sissa.it>
*
* Last changes: 03/2017
*
*
* Integration of the QPhiX library for Intel Xeon Phi usage
*
* The externally accessible functions are
*
*   void _initQphix(int argc, char **argv,
*                   int By_, int Bz_, int NCores_,
*                   int Sy_, int Sz_, int PadXY_,
*                   int PadXYZ_, int MinCt_, int c12, QphixPrec precision_)
*     Initializes the QPhiX library. Carries over the lattice size and the
*     MPI process grid and thus must be called after initializing MPI (and
*     after 'read_infile(argc,argv)').
*
*   void _endQphix()
*     Finalizes the QPhiX library. Call before MPI_Finalize().
*
**************************************************************************/

#include "qphix_interface.h"

// #undef SEEK_SET
// #undef SEEK_CUR
// #undef SEEK_END

#include "qphix_base_classes.hpp"

#ifdef TM_USE_MPI
// include mpi.h first
#include <mpi.h>
#endif

#include "config.h"
#include "global.h"
extern "C" {
#include "boundary.h"
#include "linalg/convert_eo_to_lexic.h"
#include "solver/solver.h"
#include "solver/solver_params.h"
#include "solver/solver_field.h"
#include "gettime.h"
#include "update_backward_gauge.h"
}
#ifdef TM_USE_OMP
#include <omp.h>
#endif
#include <cfloat>
#include <cstdlib>
#include <cstring>
#include <qphix/blas_new_c.h>
#include <qphix/invcg.h>
#include <qphix/invbicgstab.h>
#include <qphix/print_utils.h>
#include <qphix/qphix_config.h>
#include <qphix/wilson.h>
#include <qphix/twisted_mass.h>

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

#ifdef QPHIX_QMP_COMMS
#include <qmp.h>
#endif

int qphix_input_By = 0;
int qphix_input_Bz = 0;
int qphix_input_NCores = 0;
int qphix_input_Sy = 0;
int qphix_input_Sz = 0;
int qphix_input_PadXY = 0;
int qphix_input_PadXYZ = 0;
int qphix_input_MinCt= 1;

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
QphixPrec qphix_precision;

int subLattSize[4];
int lattSize[4];
int qmp_geom[4];

template <typename T>
struct rsdTarget {
  static const double value;
};

template <>
const double rsdTarget<half>::value = 1.0e-4;

template <>
const double rsdTarget<float>::value = 1.0e-7;

template <>
const double rsdTarget<double>::value = 1.0e-08;

void checkQphixInputParameters();

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
  qphix_precision = precision_;

  // this is called within tmLQCD already, it should not be called here
  //omp_set_num_threads(NCores * Sy * Sz);

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
void reorder_gauge_to_QPhiX(Geometry<FT, VECLEN, SOALEN, compress12> &geom, FT *qphix_gauge_cb0,
                            FT *qphix_gauge_cb1) {
  double startTime = gettime();
  double *in;
  FT *out;

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
  in = reinterpret_cast<double*>(&g_gauge_field[0][0].c00);

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
  masterPrintf("  time spent in reorder_gauge_to_QPhiX: %f secs\n", diffTime);
}

// Reorder tmLQCD spinor to a cb0 and cb1 QPhiX spinor
template <typename FT, int VECLEN, int SOALEN, bool compress12>
void reorder_spinor_to_QPhiX(Geometry<FT, VECLEN, SOALEN, compress12> &geom, double const *tm_spinor,
                             FT *qphix_spinor_cb0, FT *qphix_spinor_cb1) {
  double startTime = gettime();

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
          const double *in = tm_spinor + Ns * Nc * Nz * tm_idx;
          FT *out;
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
  masterPrintf("  time spent in reorder_spinor_to_QPhiX: %f secs\n", diffTime);
}

// Reorder a cb0 and cb1 QPhiX spinor to a tmLQCD spinor
template <typename FT, int VECLEN, int SOALEN, bool compress12>
void reorder_spinor_from_QPhiX(Geometry<FT, VECLEN, SOALEN, compress12> &geom, double *tm_spinor,
                               FT const *qphix_spinor_cb0, FT const *qphix_spinor_cb1,
                               double normFac = 1.0) {
  double startTime = gettime();

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
          const FT *in;
          if ((t + x + y + z) & 1)
            in = qphix_spinor_cb1 + SOALEN * Nz * Nc * Ns * qphix_idx;  // cb1
          else
            in = qphix_spinor_cb0 + SOALEN * Nz * Nc * Ns * qphix_idx;  // cb0
          double *out = tm_spinor + Ns * Nc * Nz * tm_idx;

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
  masterPrintf("  time spent in reorder_spinor_from_QPhiX: %f secs\n", diffTime);
}

// Apply the Dslash to a full tmlQCD spinor and return a full tmlQCD spinor
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

  // tmLQCD only stores kappa, QPhiX uses the mass. Convert here.
  double const mass = 1 / (2.0 * g_kappa) - 4;

#if 0 // Change the operator to use here.
  tmlqcd::WilsonDslash<FT, V, S, compress> concrete_dslash(&geom, t_boundary, coeff_s, coeff_t, mass);
#else
  tmlqcd::WilsonTMDslash<FT, V, S, compress> concrete_dslash(&geom, t_boundary, coeff_s, coeff_t, mass, 1.21);
#endif

  tmlqcd::Dslash<FT, V, S, compress> &polymorphic_dslash = concrete_dslash;


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

  // Reorder (global) input gauge field from tmLQCD to QPhiX,
  // which uses global tmlQCD gauge field as input
  reorder_gauge_to_QPhiX(geom,
                         reinterpret_cast<FT*>(u_packed[0]),
                         reinterpret_cast<FT*>(u_packed[1]));

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

  QSpinor *tmp_spinor = (QSpinor *)geom.allocCBFourSpinor();

  // Reorder input spinor from tmLQCD to QPhiX
  reorder_spinor_to_QPhiX(geom,
                          reinterpret_cast<double const*>(tmlqcd_in),
                          reinterpret_cast<FT*>(qphix_in[0]),
                          reinterpret_cast<FT*>(qphix_in[1]));

  // Apply QPhiX Dslash to qphix_in spinors
  polymorphic_dslash.dslash(qphix_out[1],
                            qphix_in[0],
                            u_packed[1],
                            /* isign == non-conjugate */ 1, /* cb == */
                            1);
  polymorphic_dslash.dslash(qphix_out[0],
                            qphix_in[1],
                            u_packed[0],
                            /* isign == non-conjugate */ 1, /* cb == */
                            0);

  if (std::is_same<decltype(concrete_dslash), tmlqcd::WilsonTMDslash<FT, V, S, compress>>::value) {
    for (int cb : {0, 1}) {
      polymorphic_dslash.A_chi(tmp_spinor, qphix_out[cb], 1);
      copySpinor(qphix_out[cb], tmp_spinor, geom, 1);
    }
  }

  // Reorder spinor fields back to tmLQCD
  reorder_spinor_from_QPhiX(geom,
                            reinterpret_cast<double*>(tmlqcd_out),
                            reinterpret_cast<FT*>(qphix_out[0]),
                            reinterpret_cast<FT*>(qphix_out[1]),
                            (1. * g_kappa));

  masterPrintf("Cleaning up\n");

  geom.free(packed_gauge_cb0);
  geom.free(packed_gauge_cb1);
  geom.free(packed_spinor_in_cb0);
  geom.free(packed_spinor_in_cb1);
  geom.free(packed_spinor_out_cb0);
  geom.free(packed_spinor_out_cb1);
  geom.free(tmp_spinor);
}

// Templatized even-odd preconditioned solver using QPhiX Library
template <typename FT, int V, int S, bool compress>
int invert_eo_qphix_helper(spinor * const tmlqcd_even_out,
    spinor * const tmlqcd_odd_out,
    spinor * const tmlqcd_even_in,
    spinor * const tmlqcd_odd_in,
    const double precision,
    const int max_iter,
    const int solver_flag,
    const int rel_prec,
    solver_params_t solver_params,
    const CompressionType compression) {

  typedef typename Geometry<FT, V, S, compress>::SU3MatrixBlock QGauge;
  typedef typename Geometry<FT, V, S, compress>::FourSpinorBlock QSpinor;

  /************************
   *                      *
   *    SETUP GEOMETRY    *
   *                      *
  ************************/

  // _initQphix should have been called at this point

  masterPrintf("# VECLEN = %d, SOALEN = %d\n", V, S);
  masterPrintf("# Declared QMP Topology: %d %d %d %d\n",
      qmp_geom[0], qmp_geom[1], qmp_geom[2], qmp_geom[3]);
  masterPrintf("# Global Lattice Size = %d %d %d %d\n",
      lattSize[0], lattSize[1], lattSize[2], lattSize[3]);
  masterPrintf("#  Local Lattice Size = %d %d %d %d\n",
      subLattSize[0], subLattSize[1], subLattSize[2], subLattSize[3]);
  masterPrintf("# Block Sizes: By = %d, Bz = %d\n", By, Bz);
  masterPrintf("# Cores = %d\n", NCores);
  masterPrintf("# SMT Grid: Sy = %d, Sz = %d\n", Sy, Sz);
  masterPrintf("# Pad Factors: PadXY = %d, PadXYZ = %d\n", PadXY, PadXYZ);
  masterPrintf("# Threads_per_core = %d\n", N_simt);
  masterPrintf("# MinCt = %d\n", MinCt);

  // Create a Geometry Class
  masterPrintf("# Initializing QPhiX Geometry...\n");
  Geometry<FT, V, S, compress> geom(subLattSize, By, Bz, NCores, Sy, Sz, PadXY, PadXYZ, MinCt);
  masterPrintf("# ...done.\n");

  /************************
   *                      *
   *     GAUGE FIELDS     *
   *                      *
  ************************/

  masterPrintf("# Allocating and preparing gauge fields...\n");

  // Allocate data for the gauge fields
  QGauge *u_packed[2];
  QGauge *packed_gauge_cb0 = (QGauge *)geom.allocCBGauge();
  QGauge *packed_gauge_cb1 = (QGauge *)geom.allocCBGauge();
  u_packed[0] = packed_gauge_cb0;
  u_packed[1] = packed_gauge_cb1;

  // Reorder (global) input gauge field from tmLQCD to QPhiX
  reorder_gauge_to_QPhiX(geom,
                         reinterpret_cast<FT*>(u_packed[0]),
                         reinterpret_cast<FT*>(u_packed[1]));

  masterPrintf("# ...done.\n");

  /************************
   *                      *
   *     SPINOR FIELDS    *
   *                      *
  ************************/

  masterPrintf("# Allocating fermion fields...\n");

  // Allocate data for the even/odd (checkerboarded) QPhiX in/out spinors
  QSpinor *packed_spinor_in_cb0 = (QSpinor *)geom.allocCBFourSpinor();
  QSpinor *packed_spinor_in_cb1 = (QSpinor *)geom.allocCBFourSpinor();
  QSpinor *packed_spinor_out_cb0 = (QSpinor *)geom.allocCBFourSpinor();
  QSpinor *packed_spinor_out_cb1 = (QSpinor *)geom.allocCBFourSpinor();
  QSpinor *qphix_in[2];
  QSpinor *qphix_out[2];
  qphix_in[0] = packed_spinor_in_cb0;
  qphix_in[1] = packed_spinor_in_cb1;
  qphix_out[0] = packed_spinor_out_cb0;
  qphix_out[1] = packed_spinor_out_cb1;

  // Allocate data for odd (cb0) QPhiX prepared in spinor
  // and a buffer for the CG solver (to do the M^dagger matrix
  // multiplication after the solve)
  QSpinor *qphix_in_prepared = (QSpinor *)geom.allocCBFourSpinor();
  QSpinor *qphix_buffer = (QSpinor *)geom.allocCBFourSpinor();

  // Allocate tmlQCD full spinor buffer to store lexic in- and output spinor
  spinor **solver_fields = nullptr;
  const int nr_solver_fields = 1;
  init_solver_field(&solver_fields, VOLUMEPLUSRAND, nr_solver_fields);
  spinor *tmlqcd_full_buffer = solver_fields[0];

  masterPrintf("# ...done.\n");

  /************************************************
   *                                              *
   *    SETUP DSLASH / FERMION MATRIX / SOLVER    *
   *                                              *
  ************************************************/

  // Time Boundary Conditions and Anisotropy Coefficents
  const double t_boundary = 1.0;
  const double coeff_s = 1.0;
  const double coeff_t = 1.0;

  // The Wilson mass re-express in terms of \kappa
  const double mass = 1.0 / (2.0 * g_kappa) - 4.0;

  // Create a Wilson Dslash for source preparation
  // and solution reconstruction
  masterPrintf("# Creating plain QPhiX Wilson Dslash (for preparation)...\n");
  QPhiX::Dslash<FT, V, S, compress>* WilsonDslash =
    new QPhiX::Dslash<FT, V, S, compress>(&geom, t_boundary, coeff_s, coeff_t);
  masterPrintf("# ...done.\n");

  // Create a Dslash & an even-odd preconditioned Fermion Matrix object,
  // depending on the chosen fermion action
  tmlqcd::Dslash<FT, V, S, compress>* DslashQPhiX;
  EvenOddLinearOperator<FT, V, S, compress>* FermionMatrixQPhiX;
  if( g_mu != 0.0 && g_c_sw > 0.0 ) { // TWISTED-MASS-CLOVER
    // TODO: Implement me!
    masterPrintf("# TWISTED-MASS-CLOVER CASE NOT YET IMPLEMENTED!\n");
    masterPrintf(" Aborting...\n");
    abort();
  } else if( g_mu != 0.0 ) { // TWISTED-MASS
    masterPrintf("# Creating QPhiX Twisted Mass Wilson Dslash...\n");
    const double TwistedMass = - g_mu / (2.0 * g_kappa);
    DslashQPhiX = new tmlqcd::WilsonTMDslash<FT, V, S, compress>
      (&geom, t_boundary, coeff_s, coeff_t, mass, TwistedMass);
    masterPrintf("# ...done.\n");

    masterPrintf("# Creating QPhiX Twisted Mass Wilson Fermion Matrix...\n");
    FermionMatrixQPhiX = new EvenOddTMWilsonOperator<FT, V, S, compress>
      (mass, TwistedMass, u_packed, &geom, t_boundary, coeff_s, coeff_t);
    masterPrintf("# ...done.\n");
  } else if( g_c_sw > 0.0 ) { // WILSON CLOVER
    // TODO: Implement me!
    masterPrintf("# WILSON CLOVER CASE NOT YET IMPLEMENTED!\n");
    masterPrintf(" Aborting...\n");
    abort();
  } else { // WILSON
    masterPrintf("# Creating QPhiX Wilson Dslash...\n");
    DslashQPhiX = new tmlqcd::WilsonDslash<FT, V, S, compress>(&geom, t_boundary, coeff_s, coeff_t, mass);
    masterPrintf("# ...done.\n");

    masterPrintf("# Creating QPhiX Wilson Fermion Matrix...\n");
    FermionMatrixQPhiX = new EvenOddWilsonOperator<FT, V, S, compress>
      (mass, u_packed, &geom, t_boundary, coeff_s, coeff_t);
    masterPrintf("# ...done.\n");
  }

  // Create a Linear Solver Object
  AbstractSolver<FT, V, S, compress>* SolverQPhiX;
  if(solver_flag == CG) {
    masterPrintf("# Creating CG Solver...\n");
    SolverQPhiX = new InvCG<FT, V, S, compress> (*FermionMatrixQPhiX, max_iter);
  } else if(solver_flag == BICGSTAB) {
    masterPrintf("# Creating BiCGStab Solver...\n");
    SolverQPhiX = new InvBiCGStab<FT, V, S, compress> (*FermionMatrixQPhiX, max_iter);
  } else {
    // TODO: Implement multi-shift CG, Richardson multi-precision
    masterPrintf(" Solver not yet supported by QPhiX!\n");
    masterPrintf(" Aborting...\n");
    abort();
  }
  masterPrintf("# ...done.\n");

  // Set number of BLAS threads by hand.
  // In case some implements the tune routines in QPhiX
  // this may be updated...
  masterPrintf("# Setting number of BLAS threads...\n");
  const int n_blas_simt = N_simt;
  masterPrintf("# ...done.\n");

  /************************
   *                      *
   *    PREPARE SOURCE    *
   *                      *
  ************************/

  masterPrintf("# Preparing odd source...\n");

  // 1. Reorder input spinor from tmLQCD to QPhiX:
  // a) Merge the even & odd tmlQCD input spinors to a full spinor
  // b) Convert full tmlQCD spinor to a cb0 & cb1 QPhiX spinor

  convert_eo_to_lexic(tmlqcd_full_buffer, // new full spinor
                      tmlqcd_even_in,     // even spinor
                      tmlqcd_odd_in);     // odd spinor

  reorder_spinor_to_QPhiX(geom,
                          (double*) tmlqcd_full_buffer,
                          reinterpret_cast<FT*>(qphix_in[0]),
                          reinterpret_cast<FT*>(qphix_in[1]));

  // 2. Prepare the odd (cb0) source
  //
  //      \tilde b_o = 1/2 Dslash^{Wilson}_oe A^{-1}_{ee} b_e + b_o
  //
  // in three steps:
  // a) Apply A^{-1} to b_e and save result in qphix_buffer
  // b) Apply Wilson Dslash to qphix_buffer and save result in qphix_in_prepared
  // c) Apply AYPX to rescale last result (=y) and add b_o (=x)

  DslashQPhiX->A_inv_chi(qphix_buffer, // out spinor
                         qphix_in[1],  // in spinor
                         1);           // non-conjugate
  WilsonDslash->dslash(qphix_in_prepared, // out spinor
                       qphix_buffer,      // in spinor
                       u_packed[0],       // gauge field on target cb
                       1,                 // non-conjugate
                       0);                // target cb == odd
  QPhiX::aypx(0.5, qphix_in[0], qphix_in_prepared, geom, n_blas_simt);

  masterPrintf("# ...done.\n");

  /************************
   *                      *
   *     SOLVE ON CB0     *
   *                      *
  ************************/

  masterPrintf("# Calling the solver...\n");

  // Set variables need for solve
  bool verbose = true;
  int niters = -1;
  double rsd_final = -1.0;
  uint64_t site_flops = -1;
  uint64_t mv_apps = -1;

  // Set the right precision for the QPhiX solver
  double rhs_norm2 = 1.0;
  QPhiX::norm2Spinor(rhs_norm2, qphix_in_prepared, geom, n_blas_simt);
  const double RsdTarget = sqrt(precision / rhs_norm2);

  // Calling the solver
  double start = omp_get_wtime();
  if(solver_flag == CG) {

    // USING CG:
    // We are solving
    //   M M^dagger qphix_buffer = qphix_in_prepared
    // here, that is, isign = -1 for the QPhiX CG solver.
    // After that multiply with M^dagger:
    //   qphix_out[0] = M^dagger M^dagger^-1 M^-1 qphix_in_prepared
    (*SolverQPhiX)(qphix_buffer, qphix_in_prepared, RsdTarget, niters, rsd_final,
                   site_flops, mv_apps, -1, verbose);
    (*FermionMatrixQPhiX)(qphix_out[0], qphix_buffer, /* conjugate */ -1);

  } else if(solver_flag == BICGSTAB) {

    // USING BiCGStab:
    // Solve M qphix_out[0] = qphix_in_prepared, directly.
    (*SolverQPhiX)(qphix_out[0], qphix_in_prepared, RsdTarget, niters, rsd_final,
                   site_flops, mv_apps, 1, verbose);

  }
  double end = omp_get_wtime();

  uint64_t num_cb_sites = lattSize[0]/2 * lattSize[1] * lattSize[2] * lattSize[3];
  uint64_t total_flops = (site_flops + (72 + 2*1320) * mv_apps) * num_cb_sites;
  masterPrintf("# Solver Time = %g sec\n", (end-start));
  masterPrintf("# Performance in GFLOPS = %g\n", 1.0e-9 * total_flops / (end-start));

  /**************************
   *                        *
   *  RECONSTRUCT SOLUTION  *
   *                        *
  **************************/

  masterPrintf("# Reconstruction even solution...\n");

  // 1. Reconstruct the even (cb1) solution
  //
  //      x_e = A^{-1}_{ee} (b_e + 1/2 Dslash^{Wilson}_eo x_o)
  //
  // in three steps:
  // b) Apply Wilson Dslash to x_o and save result in qphix_buffer
  // c) Apply AYPX to rescale last result (=y) and add b_e (=x)
  // c) Apply A^{-1} to qphix_buffer and save result in x_e

  WilsonDslash->dslash(qphix_buffer, // out spinor
                       qphix_out[0], // in spinor
                       u_packed[1],  // gauge field on target cb
                       1,            // non-conjugate
                       1);           // target cb == even
  QPhiX::aypx(0.5, qphix_in[1], qphix_buffer, geom, n_blas_simt);
  DslashQPhiX->A_inv_chi(qphix_out[1], // out spinor
                         qphix_buffer, // in spinor
                         1);           // non-conjugate

  // 2. Reorder spinor fields back to tmLQCD, rescaling by a factor 1/(2*\kappa)

  reorder_spinor_from_QPhiX(geom,
                            reinterpret_cast<double*>(tmlqcd_full_buffer),
                            reinterpret_cast<FT*>(qphix_out[0]),
                            reinterpret_cast<FT*>(qphix_out[1]),
                            1.0 / (2.0 * g_kappa));

  convert_lexic_to_eo(tmlqcd_even_out,     // new even spinor
                      tmlqcd_odd_out,      // new odd spinor
                      tmlqcd_full_buffer); // full spinor

  masterPrintf("# ...done.\n");

  /******************
   *                *
   *    CLEAN UP    *
   *                *
  ******************/

  masterPrintf("# Cleaning up\n");

  geom.free(packed_gauge_cb0);
  geom.free(packed_gauge_cb1);
  geom.free(packed_spinor_in_cb0);
  geom.free(packed_spinor_in_cb1);
  geom.free(packed_spinor_out_cb0);
  geom.free(packed_spinor_out_cb1);
  geom.free(qphix_in_prepared);
  geom.free(qphix_buffer);
  finalize_solver(solver_fields, nr_solver_fields);

  // FIXME: This should be called properly somewhere else
  _endQphix();

  masterPrintf("# ...done.\n\n");

  return niters;
}


// Template wrapper for the Dslash operator call-able from C code
void D_psi_qphix(spinor* tmlqcd_out, const spinor* tmlqcd_in) {
  if (qphix_precision == QPHIX_DOUBLE_PREC) {
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
  } else if (qphix_precision == QPHIX_FLOAT_PREC) {
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
  else if (qphix_precision == QPHIX_HALF_PREC) {
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

// Template wrapper for Full Solver call-able from C code, return number of iterations
int invert_eo_qphix(spinor * const Even_new,
    spinor * const Odd_new,
    spinor * const Even,
    spinor * const Odd,
    const double precision,
    const int max_iter,
    const int solver_flag,
    const int rel_prec,
    solver_params_t solver_params,
    const SloppyPrecision sloppy,
    const CompressionType compression) {


  checkQphixInputParameters();

  if (precision < rsdTarget<double>::value) {
    if (QPHIX_SOALEN > VECLEN_DP) {
      masterPrintf("SOALEN=%d is greater than the double prec VECLEN=%d\n", QPHIX_SOALEN,
                   VECLEN_DP);
      abort();
    }
    masterPrintf("# INITIALIZING QPHIX SOLVER\n");
    masterPrintf("# USING DOUBLE PRECISION\n");
    _initQphix(0, nullptr, qphix_input_By, qphix_input_Bz, qphix_input_NCores,
           qphix_input_Sy, qphix_input_Sz,
           qphix_input_PadXY, qphix_input_PadXYZ, qphix_input_MinCt,
           compression, QPHIX_DOUBLE_PREC);

    if (compress12) {
      return invert_eo_qphix_helper<double, VECLEN_DP, QPHIX_SOALEN, true>
        (Even_new,
         Odd_new,
         Even,
         Odd,
         precision,
         max_iter,
         solver_flag,
         rel_prec,
         solver_params,
         compression);
    } else {
      return invert_eo_qphix_helper<double, VECLEN_DP, QPHIX_SOALEN, false>
        (Even_new,
         Odd_new,
         Even,
         Odd,
         precision,
         max_iter,
         solver_flag,
         rel_prec,
         solver_params,
         compression);
    }
#ifdef QPHIX_MIC_SOURCE
  } else if (precision < rsdTarget<float>::value) {
#else
  } else {
#endif
    if (QPHIX_SOALEN > VECLEN_SP) {
      masterPrintf("SOALEN=%d is greater than the single prec VECLEN=%d\n", QPHIX_SOALEN,
                   VECLEN_SP);
      abort();
    }
    masterPrintf("# INITIALIZING QPHIX SOLVER\n");
    masterPrintf("# USING SINGLE PRECISION\n");
    _initQphix(0, nullptr, qphix_input_By, qphix_input_Bz, qphix_input_NCores,
           qphix_input_Sy, qphix_input_Sz,
           qphix_input_PadXY, qphix_input_PadXYZ, qphix_input_MinCt,
           compression, QPHIX_FLOAT_PREC);

    if (compress12) {
      return invert_eo_qphix_helper<float, VECLEN_SP, QPHIX_SOALEN, true>
        (Even_new,
         Odd_new,
         Even,
         Odd,
         precision,
         max_iter,
         solver_flag,
         rel_prec,
         solver_params,
         compression);
    } else {
      return invert_eo_qphix_helper<float, VECLEN_SP, QPHIX_SOALEN, false>
        (Even_new,
         Odd_new,
         Even,
         Odd,
         precision,
         max_iter,
         solver_flag,
         rel_prec,
         solver_params,
         compression);
    }
  }
#if defined(QPHIX_MIC_SOURCE)
  else if (precision < rsdTarget<half>::value) {
    if (QPHIX_SOALEN > VECLEN_HP) {
      masterPrintf("SOALEN=%d is greater than the single prec VECLEN=%d\n", QPHIX_SOALEN,
                   VECLEN_SP);
      abort();
    }
    masterPrintf("# INITIALIZING QPHIX SOLVER\n");
    masterPrintf("# USING HALF PRECISION\n");
    _initQphix(0, nullptr, qphix_input_By, qphix_input_Bz, qphix_input_NCores,
           qphix_input_Sy, qphix_input_Sz,
           qphix_input_PadXY, qphix_input_PadXYZ, qphix_input_MinCt,
           compression, QPHIX_HALF_PREC);
    
    if (compress12) {
      return invert_eo_qphix_helper<half, VECLEN_SP, QPHIX_SOALEN, true>
        (Even_new,
         Odd_new,
         Even,
         Odd,
         precision,
         max_iter,
         solver_flag,
         rel_prec,
         solver_params,
         compression);
    } else {
      return invert_eo_qphix_helper<half, VECLEN_SP, QPHIX_SOALEN, false>
        (Even_new,
         Odd_new,
         Even,
         Odd,
         precision,
         max_iter,
         solver_flag,
         rel_prec,
         solver_params,
         compression);
    }
  }
#endif

  return -1;
}

void checkQphixInputParameters() {
  if( qphix_input_By == 0 || qphix_input_Bz == 0){
    masterPrintf("QPHIX Error: By and Bz may not be 0 ! Aborting.\n");
    abort();
  }
  if( qphix_input_NCores * qphix_input_Sy * qphix_input_Sz != omp_num_threads ){
    masterPrintf("QPHIX Error: NCores * Sy * Sz != ompnumthreads ! Aborting.\n");
    abort();
  }
} 

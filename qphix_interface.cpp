/***********************************************************************
 *
 * Copyright (C) 2015 Mario Schroeck
 *               2016 Peter Labus
 *               2017 Peter Labus, Martin Ueding, Bartosz Kostrzewa
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

#include "qphix_interface.h"
#include "qphix_base_classes.hpp"
#include "qphix_interface_utils.hpp"
#include "qphix_types.h"
#include "qphix_veclen.h"

#ifdef TM_USE_MPI
#include <mpi.h>
#endif

#include "config.h"
#include "global.h"
extern "C" {
#include "boundary.h"
#include "geometry_eo.h"
#include "gettime.h"
#include "linalg/convert_eo_to_lexic.h"
#include "linalg/diff.h"
#include "linalg/square_norm.h"
#include "linalg/square_norm.h"
#include "operator/clovertm_operators.h"
#include "operator_types.h"
#include "solver/solver.h"
#include "solver/solver_field.h"
#include "solver/solver_params.h"
#include "xchange/xchange_gauge.h"
}
#ifdef TM_USE_OMP
#include <omp.h>
#endif
#include <qphix/blas_new_c.h>
#include <qphix/invbicgstab.h>
#include <qphix/invcg.h>
#include <qphix/print_utils.h>
#include <qphix/qphix_config.h>
#include <qphix/twisted_mass.h>
#include <qphix/wilson.h>
#include <cfloat>
#include <cstdlib>
#include <cstring>
#include <vector>

using namespace tmlqcd;

QphixParams_t qphix_input;

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
QphixPrec_t qphix_precision;

int subLattSize[4];
int lattSize[4];
int qmp_geom[4];
int qmp_tm_map[4];

template <typename T>
struct rsdTarget {
  static const double value;
};

template <>
const double rsdTarget<QPhiX::half>::value = 1.0e-3;

template <>
const double rsdTarget<float>::value = 1.0e-8;

void _initQphix(int argc, char **argv, QphixParams_t params, int c12, QphixPrec_t precision_) {
  static bool qmp_topo_initialised = false;

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

  By = params.By;
  Bz = params.Bz;
  NCores = params.NCores;
  Sy = params.Sy;
  Sz = params.Sz;
  PadXY = params.PadXY;
  PadXYZ = params.PadXYZ;
  MinCt = params.MinCt;
  N_simt = Sy * Sz;
  if (c12 == 8) {
    QPhiX::masterPrintf(
        "# INFO QphiX: 8-parameter gauge compression not supported, using two row compression "
        "instead!\n");
    c12 = 12;
  }
  compress12 = c12 == 12 ? true : false;
  qphix_precision = precision_;

#ifdef QPHIX_QMP_COMMS
  // Declare the logical topology
  if (!qmp_topo_initialised) {
    // the QMP topology is the one implied by the number of processes in each
    // dimension as required by QPHIX ( x fastest to t slowest running )
    qmp_geom[0] = g_nproc_x; 
    qmp_geom[1] = g_nproc_y;
    qmp_geom[2] = g_nproc_z;
    qmp_geom[3] = g_nproc_t;

    // in order for the topologies to agree between tmLQCD and QPhiX, the dimensions need to be permuted
    // since Z is fastest in tmLQCD and X is second-slowest
    qmp_tm_map[0] = 2;
    qmp_tm_map[1] = 1;
    qmp_tm_map[2] = 0;
    qmp_tm_map[3] = 3; 
    if (QMP_declare_logical_topology_map(qmp_geom, 4, qmp_tm_map, 4) != QMP_SUCCESS) {
      QMP_error("Failed to declare QMP Logical Topology\n");
      abort();
    }
    // longish test to check if the logical coordinates are correctly mapped
    if(g_debug_level >= 5){
      for(int proc = 0; proc < g_nproc; proc++){
        if( proc == g_proc_id ){
          const int coordinates[4] = { g_proc_coords[1], g_proc_coords[2], g_proc_coords[3], g_proc_coords[0] };
          int id = QMP_get_node_number_from( coordinates );
          int* qmp_coords = QMP_get_logical_coordinates_from(id);
          fflush(stdout);
          printf("QMP id: %3d x:%3d y:%3d z:%3d t:%3d\n", id, 
                 qmp_coords[0], qmp_coords[1], qmp_coords[2], qmp_coords[3]); 
          printf("MPI id: %3d x:%3d y:%3d z:%3d t:%3d\n\n", g_proc_id, 
                 g_proc_coords[1], g_proc_coords[2], g_proc_coords[3], g_proc_coords[0]);
          free(qmp_coords);
          fflush(stdout);
          MPI_Barrier(MPI_COMM_WORLD);
        } else {
          MPI_Barrier(MPI_COMM_WORLD);
        }
      }
    }
    qmp_topo_initialised = true;
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

// Reorder the tmLQCD clover field to a Wilson QPhiX clover field
template <typename FT, int VECLEN, int SOALEN, bool compress12>
void reorder_clover_to_QPhiX(QPhiX::Geometry<FT, VECLEN, SOALEN, compress12> &geom,
                             FT *qphix_clover) {
  const double startTime = gettime();

  std::vector<uint64_t> const lookup_table{
      0,  8,  16, 72, 80, 88,  2,  3,  4,  5,  36, 37, 38, 39, 40, 41, 10,  11,
      42, 43, 44, 45, 46, 47,  48, 49, 50, 51, 52, 53, 74, 75, 76, 77, 82,  83,
      18, 26, 34, 90, 98, 106, 20, 21, 22, 23, 54, 55, 56, 57, 58, 59, 28,  29,
      60, 61, 62, 63, 64, 65,  66, 67, 68, 69, 70, 71, 92, 93, 94, 95, 100, 101};

  assert(lookup_table.size() == 72 &&
         "Number of elements in the lookup table needs to exactly match the number of elements in "
         "the clover data structure.");

  // Need to rescale for QPhiX by alpha
  const double alpha = 1.0 / (2.0 * g_kappa);

  // Number of elements in spin, color & complex
  const int Ns = 4;
  const int Nc = 3;
  const int Nz = 2;

  // Geometric parameters for QPhiX data layout
  const auto nyg = geom.nGY();
  const auto nVecs = geom.nVecs();
  const auto Pxy = geom.getPxy();
  const auto Pxyz = geom.getPxyz();

  // Get the base pointer for the (global) tmlQCD clover field
  const double *in = reinterpret_cast<double *>(&sw[0][0][0].c00);

  FT *out = qphix_clover;

  // This will loop over the entire lattice and calculate
  // the array and internal indices for both tmlQCD & QPhiX
  for (uint64_t t = 0; t < T; t++)
    for (uint64_t x = 0; x < LX; x++)
      for (uint64_t y = 0; y < LY; y++)
        for (uint64_t z = 0; z < LZ; z++) {
          // Only copy the odd-odd elements of the clover term
          // XXX This looks like `0` is `even`. In other parts of the QPhiX code, it seems like `1`
          // is `even`. Either way, there is a magic number here which should be refactored into a
          // global constant, probably best an enum. Then also this should be a function like
          // `is_even` or `is_odd` to make the code more readable and less error prone.
          if ((t + x + y + z) & 1 == 0) continue;

          const uint64_t tm_idx = Index(t, x, y, z);

          // These are the QPhiX SIMD vector in checkerboarded x direction
          // (up to LX/2), the index inside one single Structure of Arrays (SOA)
          // and the internal position inside the ("packed") SIMD vector
          const uint64_t SIMD_vector = (x / 2) / SOALEN;
          const uint64_t x_one_SOA = (x / 2) % SOALEN;
          const uint64_t x_internal = (y % nyg) * SOALEN + x_one_SOA;

          // Calculate the array index in QPhiX, given a global
          // lattice index (t,x,y,z). This is also called the
          // "block" in QPhiX and the QPhiX/QDP packers.
          const uint64_t qphix_idx = (t * Pxyz + z * Pxy) / nyg + (y / nyg) * nVecs + SIMD_vector;

          // Calculate the index where the tile for a given site begins
          const uint64_t qphix_base = qphix_idx * VECLEN * 2 * (6 + 15 * Nz);
          const uint64_t tm_base = tm_idx * 3 * 2 * Nc * Nc * Nz;

          // FIXME Is there any better way to do that?
          for (uint64_t lt_idx = 0; lt_idx < lookup_table.size(); ++lt_idx) {
            out[qphix_base + lt_idx * VECLEN + x_internal] =
                alpha * in[tm_base + lookup_table[lt_idx]];
          }

        }  // volume

  const double endTime = gettime();
  const double diffTime = endTime - startTime;
  QPhiX::masterPrintf("  time spent in reorder_clover_to_QPhiX: %f secs\n", diffTime);
}

template <typename FT, int VECLEN, int SOALEN, bool compress12>
void reorder_gauge_to_QPhiX(QPhiX::Geometry<FT, VECLEN, SOALEN, compress12> &geom,
                            typename QPhiX::Geometry<FT, VECLEN, SOALEN, compress12>::SU3MatrixBlock *qphix_gauge_cb0,
                            typename QPhiX::Geometry<FT, VECLEN, SOALEN, compress12>::SU3MatrixBlock *qphix_gauge_cb1){
  const double startTime = gettime();

  // Number of elements in spin, color & complex
  // Here c1 is QPhiX's outer color, and c2 the inner one
  const int Ns = 4;
  const int Nc1 = compress12 ? 2 : 3;
  const int Nc2 = 3;
  const int Nz = 2;

  // Geometric parameters for QPhiX data layout
  const auto ngy = geom.nGY();
  const auto nVecs = geom.nVecs();
  const auto Pxy = geom.getPxy();
  const auto Pxyz = geom.getPxyz();

  // This is needed to translate between the different
  // orderings of the direction index "\mu" in tmlQCD
  // and QPhiX, respectively
  // in qphix, the Dirac operator is applied in the order
  //   -+x -> -+y -> -+z -> -+t
  // while tmlqcd does
  //   -+t -> -+x -> -+y -> -+z
  // same as the lattice ordering
  // The mappingn between the application dimensions is thus:
  //  tmlqcd_dim(t(0) -> x(1) -> y(2) -> z(3)) = qphix_dim( t(3) -> x(0) -> y(1) -> z(2) )
  const int change_dim[4] = {1, 2, 3, 0};

  // Get the base pointer for the (global) tmlQCD gauge field
  xchange_gauge(g_gauge_field);
  const double *in = reinterpret_cast<double *>(&g_gauge_field[0][0].c00);

#pragma omp parallel for collapse(4)
  for (int64_t t = 0; t < T; t++)
    for (int64_t z = 0; z < LZ; z++)
      for (int64_t y = 0; y < LY; y++)
        for (int64_t v = 0; v < nVecs; v++){

          int64_t block = (t * Pxyz + z * Pxy) / ngy + (y / ngy) * nVecs + v;

          for (int dim = 0; dim < 4; dim++)  // dimension == QPhiX \mu
            for (int c1 = 0; c1 < Nc1; c1++)    // QPhiX convention color 1 (runs up to 2 or 3)
              for (int c2 = 0; c2 < Nc2; c2++)  // QPhiX convention color 2 (always runs up to 3)
                for (int x_soa = 0; x_soa < SOALEN; x_soa++){
                  int64_t xx = ( y % ngy ) * SOALEN + x_soa;
                  int64_t q_cb_x_coord = x_soa + v*SOALEN;
                  int64_t tm_x_coord_cb0 = q_cb_x_coord * 2 + ( ((t + y + z) & 1) ^ 0 );
                  int64_t tm_x_coord_cb1 = q_cb_x_coord * 2 + ( ((t + y + z) & 1) ^ 1);
                  
                  int64_t tm_idx_cb0;
                  int64_t tm_idx_cb1;
                 
                  // backward / forward
                  for (int dir = 0; dir < 2; dir++) {
                    if (dir == 0){
                      tm_idx_cb0 = g_idn[g_ipt[t][tm_x_coord_cb0][y][z]][change_dim[dim]];
                      tm_idx_cb1 = g_idn[g_ipt[t][tm_x_coord_cb1][y][z]][change_dim[dim]];
                    } else {
                      tm_idx_cb0 = g_ipt[t][tm_x_coord_cb0][y][z];
                      tm_idx_cb1 = g_ipt[t][tm_x_coord_cb1][y][z];
                    }
                    for (int reim = 0; reim < Nz; reim++) {
                      // Note:
                      // -----
                      // 1. \mu in QPhiX runs from 0..7 for all eight neighbouring
                      // links.
                      //    Here, the ordering of the direction (backward/forward)
                      //    is the same
                      //    for tmlQCD and QPhiX, but we have to change the
                      //    ordering of the dimensions.
                      int q_mu = 2 * dim + dir;

                      // 2. QPhiX gauge field matrices are transposed w.r.t.
                      // tmLQCD.
                      // 3. tmlQCD always uses 3x3 color matrices (Nc2*Nc2).
                      int64_t t_inner_idx_cb0 = reim + c1 * Nz + c2 * Nz * Nc2 +
                                                change_dim[dim] * Nz * Nc2 * Nc2 + tm_idx_cb0 * Nz * Nc2 * Nc2 * 4;
                      int64_t t_inner_idx_cb1 = reim + c1 * Nz + c2 * Nz * Nc2 +
                                                change_dim[dim] * Nz * Nc2 * Nc2 + tm_idx_cb1 * Nz * Nc2 * Nc2 * 4;
                      qphix_gauge_cb0[block][q_mu][c1][c2][reim][xx] = in[t_inner_idx_cb0];
                      qphix_gauge_cb1[block][q_mu][c1][c2][reim][xx] = in[t_inner_idx_cb1];
                    }
                  }
                } // for(dim,c1,c2,x_soa)
        } // outer loop (t,z,y,v)
       
  const double diffTime = gettime() - startTime;
  QPhiX::masterPrintf("  time spent in reorder_gauge_to_QPhiX: %f secs\n", diffTime);
}

// Reorder tmLQCD eo-spinor to a FourSpinorBlock QPhiX spinor on the given checkerboard
template <typename FT, int VECLEN, int SOALEN, bool compress12>
void reorder_eo_spinor_to_QPhiX(
    QPhiX::Geometry<FT, VECLEN, SOALEN, compress12> &geom, double const *const tm_eo_spinor,
    typename QPhiX::Geometry<FT, VECLEN, SOALEN, compress12>::FourSpinorBlock *qphix_spinor,
    const int cb) {
  const double startTime = gettime();

  const int Ns = 4;
  const int Nc = 3;
  const int Nz = 2;

  const auto nVecs = geom.nVecs();
  const auto Pxy = geom.getPxy();
  const auto Pxyz = geom.getPxyz();
  const auto Nxh = geom.Nxh();

  // This is needed to translate between the different
  // gamma bases tmlQCD and QPhiX are using
  // (note, this is a 4x4 matrix with 4 non-zero elements)
  const int change_sign[4] = {1, -1, -1, 1};
  const int change_spin[4] = {3, 2, 1, 0};

#pragma omp parallel for collapse(4)
  for (int64_t t = 0; t < T; t++) {
    for (int64_t z = 0; z < LZ; z++) {
      for (int64_t y = 0; y < LY; y++) {
        for (int64_t v = 0; v < nVecs; v++) {
          for (int col = 0; col < Nc; col++) {
            for (int q_spin = 0; q_spin < Ns; q_spin++) {
              for (int x_soa = 0; x_soa < SOALEN; x_soa++) {
                int64_t q_ind = t * Pxyz + z * Pxy + y * nVecs + v;
                int64_t q_cb_x_coord = v * SOALEN + x_soa;
                // when t+y+z is odd and we're on an odd (1) checkerboard OR
                // when t+y+z is even and we're on an even (0) checkerboard
                // the full x coordinate is 2*x_cb
                // otherwise, it is 2*x_cb+1
                int64_t tm_x_coord = q_cb_x_coord * 2 + (((t + y + z) & 1) ^ cb);
                // exchange x and z dimensions
                int64_t tm_eo_ind = g_lexic2eosub[g_ipt[t][tm_x_coord][y][z]];
                int64_t tm_eo_offset =
                    tm_eo_ind * Nc * Ns * Nz + change_spin[q_spin] * Nc * Nz + Nz * col;
                for (int reim = 0; reim < 2; reim++) {
                  qphix_spinor[q_ind][col][q_spin][reim][x_soa] =
                      change_sign[q_spin] * tm_eo_spinor[tm_eo_offset + reim];
                }
              }
            }
          }
        }
      }
    }
  }
}

template <typename FT, int VECLEN, int SOALEN, bool compress12>
void reorder_eo_spinor_from_QPhiX(
    QPhiX::Geometry<FT, VECLEN, SOALEN, compress12> &geom, double *const tm_eo_spinor,
    typename QPhiX::Geometry<FT, VECLEN, SOALEN, compress12>::FourSpinorBlock *qphix_spinor,
    const int cb, double normFac = 1.0) {
  const double startTime = gettime();

  const int Ns = 4;
  const int Nc = 3;
  const int Nz = 2;

  const auto nVecs = geom.nVecs();
  const auto Pxy = geom.getPxy();
  const auto Pxyz = geom.getPxyz();
  const auto Nxh = geom.Nxh();

  // This is needed to translate between the different
  // gamma bases tmlQCD and QPhiX are using
  // (note, this is a 4x4 matrix with 4 non-zero elements)
  const int change_sign[4] = {1, -1, -1, 1};
  const int change_spin[4] = {3, 2, 1, 0};

#pragma omp parallel for collapse(4)
  for (int64_t t = 0; t < T; t++) {
    for (int64_t z = 0; z < LZ; z++) {
      for (int64_t y = 0; y < LY; y++) {
        for (int64_t v = 0; v < nVecs; v++) {
          for (int col = 0; col < Nc; col++) {
            for (int q_spin = 0; q_spin < Ns; q_spin++) {
              for (int x_soa = 0; x_soa < SOALEN; x_soa++) {
                int64_t q_ind = t * Pxyz + z * Pxy + y * nVecs + v;
                int64_t q_cb_x_coord = v * SOALEN + x_soa;
                // when t+y+z is odd and we're on an odd checkerboard (1) OR
                // when t+y+z is even and we're on an even (0) checkerboard
                // the full x coordinate is 2*x_cb
                // otherwise, it is 2*x_cb+1
                int64_t tm_x_coord = q_cb_x_coord * 2 + (((t + y + z) & 1) ^ cb);
                // exchange x and z dimensions
                int64_t tm_eo_ind = g_lexic2eosub[g_ipt[t][tm_x_coord][y][z]];
                int64_t tm_eo_offset =
                    tm_eo_ind * Nc * Ns * Nz + change_spin[q_spin] * Nc * Nz + Nz * col;
                for (int reim = 0; reim < 2; reim++) {
                  tm_eo_spinor[tm_eo_offset + reim] =
                      change_sign[q_spin] * normFac * qphix_spinor[q_ind][col][q_spin][reim][x_soa];
                }
              }
            }
          }
        }
      }
    }
  }
}

// Reorder a full tmLQCD spinor to a cb0 and cb1 QPhiX spinor
template <typename FT, int VECLEN, int SOALEN, bool compress12>
void reorder_spinor_to_QPhiX(QPhiX::Geometry<FT, VECLEN, SOALEN, compress12> &geom,
                             double const *tm_spinor, FT *qphix_spinor_cb0, FT *qphix_spinor_cb1) {
  const double startTime = gettime();

  // Number of elements in spin, color & complex
  const int Ns = 4;
  const int Nc = 3;
  const int Nz = 2;

  // Geometric parameters for QPhiX data layout
  const auto nVecs = geom.nVecs();
  const auto Pxy = geom.getPxy();
  const auto Pxyz = geom.getPxyz();

  // This is needed to translate between the different
  // gamma bases tmlQCD and QPhiX are using
  const int change_sign[4] = {1, -1, -1, 1};
  const int change_spin[4] = {3, 2, 1, 0};

// This will loop over the entire lattice and calculate
// the array and internal indices for both tmlQCD & QPhiX
#pragma omp parallel for collapse(4)
  for (uint64_t t = 0; t < T; t++)
    for (uint64_t x = 0; x < LX; x++)
      for (uint64_t y = 0; y < LY; y++)
        for (uint64_t z = 0; z < LZ; z++) {
          // These are the QPhiX SIMD vector in checkerboarded x direction
          // (up to LX/2) and the internal position inside the SIMD vector
          const uint64_t SIMD_vector = (x / 2) / SOALEN;
          const uint64_t x_internal = (x / 2) % SOALEN;

          // Calculate the array index in tmlQCD & QPhiX,
          // given a global lattice index (t,x,y,z)
          const uint64_t qphix_idx = t * Pxyz + z * Pxy + y * nVecs + SIMD_vector;
          const uint64_t tm_idx = g_ipt[t][x][y][z];

          // Calculate base point for every spinor field element (tmlQCD) or
          // for every SIMD vector of spinors, a.k.a FourSpinorBlock (QPhiX),
          // which will depend on the checkerboard (cb)
          const double *in = tm_spinor + Ns * Nc * Nz * tm_idx;
          FT *out;
          if ((t + x + y + z) & 1)
            out = qphix_spinor_cb1 + SOALEN * Nz * Nc * Ns * qphix_idx;  // odd -> cb1
          else
            out = qphix_spinor_cb0 + SOALEN * Nz * Nc * Ns * qphix_idx;  // even -> cb0

          // Copy the internal elements, performing a gamma basis transformation
          for (int spin = 0; spin < Ns; spin++)  // QPhiX spin index
            for (int color = 0; color < Nc; color++)
              for (int z = 0; z < Nz; z++)  // RE or IM
              {
                const uint64_t qId =
                    x_internal + z * SOALEN + spin * SOALEN * Nz + color * SOALEN * Nz * Ns;
                const uint64_t tId = z + color * Nz + change_spin[spin] * Nz * Nc;

                out[qId] = change_sign[spin] * in[tId];
              }

        }  // volume

  const double endTime = gettime();
  const double diffTime = endTime - startTime;
  if (g_debug_level > 1)
    QPhiX::masterPrintf("# QPHIX-interface:  time spent in reorder_spinor_to_QPhiX: %f secs\n",
                        diffTime);
}

// Reorder a cb0 and cb1 QPhiX spinor to a full tmLQCD spinor
template <typename FT, int VECLEN, int SOALEN, bool compress12>
void reorder_spinor_from_QPhiX(QPhiX::Geometry<FT, VECLEN, SOALEN, compress12> &geom,
                               double *tm_spinor, FT const *qphix_spinor_cb0,
                               FT const *qphix_spinor_cb1, double normFac = 1.0) {
  const double startTime = gettime();

  // Number of elements in spin, color & complex
  const int Ns = 4;
  const int Nc = 3;
  const int Nz = 2;

  // Geometric parameters for QPhiX data layout
  const auto nVecs = geom.nVecs();
  const auto Pxy = geom.getPxy();
  const auto Pxyz = geom.getPxyz();

  // This is needed to translate between the different
  // gamma bases tmlQCD and QPhiX are using
  const int change_sign[4] = {1, -1, -1, 1};
  const int change_spin[4] = {3, 2, 1, 0};

// This will loop over the entire lattice and calculate
// the array and internal indices for both tmlQCD & QPhiX
#pragma omp parallel for collapse(4)
  for (uint64_t t = 0; t < T; t++)
    for (uint64_t x = 0; x < LX; x++)
      for (uint64_t y = 0; y < LY; y++)
        for (uint64_t z = 0; z < LZ; z++) {
          // These are the QPhiX SIMD vector in checkerboarded x direction
          // (up to LX/2) and the internal position inside the SIMD vector
          const uint64_t SIMD_vector = (x / 2) / SOALEN;
          const uint64_t x_internal = (x / 2) % SOALEN;

          // Calculate the array index in tmlQCD & QPhiX,
          // given a global lattice index (t,x,y,z)
          const uint64_t qphix_idx = t * Pxyz + z * Pxy + y * nVecs + SIMD_vector;
          const uint64_t tm_idx = g_ipt[t][x][y][z];

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
                const uint64_t qId = x_internal + z * SOALEN + change_spin[spin] * SOALEN * Nz +
                                     color * SOALEN * Nz * Ns;
                const uint64_t tId = z + color * Nz + spin * Nz * Nc;

                out[tId] = normFac * change_sign[spin] * in[qId];
              }

        }  // volume

  const double endTime = gettime();
  const double diffTime = endTime - startTime;
  QPhiX::masterPrintf("  time spent in reorder_spinor_from_QPhiX: %f secs\n", diffTime);
}

// Apply the full QPhiX fermion matrix to checkerboarded tm spinors
template <typename FT, int V, int S, bool compress>
void Mfull_helper(spinor *Even_out, spinor *Odd_out, const spinor *Even_in, const spinor *Odd_in,
                  const op_type_t op_type) {
  // TODO: this should use handles for gauge and spinors because these are definitely temporary
  // objects
  typedef typename QPhiX::Geometry<FT, V, S, compress>::SU3MatrixBlock QGauge;
  typedef typename QPhiX::Geometry<FT, V, S, compress>::FourSpinorBlock QSpinor;

  if (g_debug_level > 1) tmlqcd::printQphixDiagnostics(V, S, compress);

  // Create Dslash Class
  double t_boundary = (FT)(1);
  double coeff_s = (FT)(1);
  double coeff_t = (FT)(1);

  QPhiX::Geometry<FT, V, S, compress> geom(subLattSize, By, Bz, NCores, Sy, Sz, PadXY, PadXYZ,
                                           MinCt);

  double const mass = 1 / (2.0 * g_kappa) - 4;

  tmlqcd::Dslash<FT, V, S, compress> *polymorphic_dslash;

  if (op_type == WILSON) {
    polymorphic_dslash =
        new tmlqcd::WilsonDslash<FT, V, S, compress>(&geom, t_boundary, coeff_s, coeff_t, mass);
  } else if (op_type == TMWILSON) {
    polymorphic_dslash = new tmlqcd::WilsonTMDslash<FT, V, S, compress>(
        &geom, t_boundary, coeff_s, coeff_t, mass, -g_mu / (2.0 * g_kappa));
  } else if (op_type == CLOVER && g_mu <= DBL_EPSILON) {
    QPhiX::masterPrintf(
        "tmlqcd::Mfull_helper; Wilson clover operator pass-through not implemented yet\n");
    abort();
  } else if (op_type == CLOVER && g_mu > DBL_EPSILON) {
    QPhiX::masterPrintf(
        "tmlqcd::Mfull_helper; Twisted clover operator pass-through not implemented yet\n");
    abort();
  } else {
    QPhiX::masterPrintf("tmlqcd::Mfull_helper; No such operator type: %d\n", op_type);
    abort();
  }

  /************************
   *                      *
   *     GAUGE FIELDS     *
   *                      *
  ************************/

  // Allocate data for the gauge fields
  QGauge *u_packed[2];
  QGauge *packed_gauge_cb0 =
      (QGauge *)geom.allocCBGauge();  // Links emanating from EVEN sites (cb=0)
  QGauge *packed_gauge_cb1 =
      (QGauge *)geom.allocCBGauge();  // Links emanating from ODD sites (cb=1)
  u_packed[cb_even] = packed_gauge_cb0;
  u_packed[cb_odd] = packed_gauge_cb1;

  reorder_gauge_to_QPhiX(geom, u_packed[cb_even], u_packed[cb_odd]);

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
  qphix_in[cb_even] = packed_spinor_in_cb0;
  qphix_in[cb_odd] = packed_spinor_in_cb1;
  qphix_out[cb_even] = packed_spinor_out_cb0;
  qphix_out[cb_odd] = packed_spinor_out_cb1;

  QSpinor *tmp_spinor = (QSpinor *)geom.allocCBFourSpinor();

  reorder_eo_spinor_to_QPhiX(geom, reinterpret_cast<double const *const>(Even_in),
                             qphix_in[cb_even], cb_even);
  reorder_eo_spinor_to_QPhiX(geom, reinterpret_cast<double const *const>(Odd_in), qphix_in[cb_odd],
                             cb_odd);

  // Apply QPhiX Mfull
  polymorphic_dslash->plain_dslash(qphix_out[cb_odd], qphix_in[cb_even], u_packed[cb_odd],
                                   /* isign == non-conjugate */ 1, cb_odd);
  polymorphic_dslash->plain_dslash(qphix_out[cb_even], qphix_in[cb_odd], u_packed[cb_even],
                                   /* isign == non-conjugate */ 1, cb_even);
  for (int cb : {0, 1}) {
    polymorphic_dslash->A_chi(tmp_spinor, qphix_in[cb], 1);
    QPhiX::aypx(-0.5, tmp_spinor, qphix_out[cb], geom, 1);
  }

  reorder_eo_spinor_from_QPhiX(geom, reinterpret_cast<double *>(Even_out), qphix_out[cb_even],
                               cb_even, 2.0 * g_kappa);
  reorder_eo_spinor_from_QPhiX(geom, reinterpret_cast<double *>(Odd_out), qphix_out[cb_odd], cb_odd,
                               2.0 * g_kappa);

  geom.free(packed_gauge_cb0);
  geom.free(packed_gauge_cb1);
  geom.free(packed_spinor_in_cb0);
  geom.free(packed_spinor_in_cb1);
  geom.free(packed_spinor_out_cb0);
  geom.free(packed_spinor_out_cb1);
  geom.free(tmp_spinor);
  delete (polymorphic_dslash);
}

// Templated even-odd preconditioned solver using QPhiX Library
template <typename FT, int V, int S, bool compress>
int invert_eo_qphix_helper(spinor *const tmlqcd_even_out, spinor *const tmlqcd_odd_out,
                           spinor *const tmlqcd_even_in, spinor *const tmlqcd_odd_in,
                           const double precision, const int max_iter, const int solver_flag,
                           const int rel_prec, solver_params_t solver_params,
                           const CompressionType compression) {
  // TODO: it would perhaps be beneficial to keep the fields resident
  typedef typename QPhiX::Geometry<FT, V, S, compress>::SU3MatrixBlock QGauge;
  typedef typename QPhiX::Geometry<FT, V, S, compress>::FourSpinorBlock QSpinor;

  /************************
   *                      *
   *    SETUP GEOMETRY    *
   *                      *
  ************************/

  if (g_debug_level > 1) {
    tmlqcd::printQphixDiagnostics(V, S, compress);
  }

  // Create a Geometry Class
  QPhiX::masterPrintf("# Initializing QPhiX Geometry...\n");
  QPhiX::Geometry<FT, V, S, compress> geom(subLattSize, By, Bz, NCores, Sy, Sz, PadXY, PadXYZ,
                                           MinCt);
  QPhiX::masterPrintf("# ...done.\n");

  /************************
   *                      *
   *     GAUGE FIELDS     *
   *                      *
  ************************/

  QPhiX::masterPrintf("# Allocating and preparing gauge fields...\n");

  // Allocate data for the gauge fields
  QGauge *u_packed[2];
  QGauge *packed_gauge_cb0 = (QGauge *)geom.allocCBGauge();
  QGauge *packed_gauge_cb1 = (QGauge *)geom.allocCBGauge();
  u_packed[cb_even] = packed_gauge_cb0;
  u_packed[cb_odd] = packed_gauge_cb1;

  // Reorder (global) input gauge field from tmLQCD to QPhiX
  reorder_gauge_to_QPhiX(geom, u_packed[cb_even], u_packed[cb_odd]);

  QPhiX::masterPrintf("# ...done.\n");

  /************************
   *                      *
   *     SPINOR FIELDS    *
   *                      *
  ************************/

  QPhiX::masterPrintf("# Allocating fermion fields...\n");

  // Allocate data for the even/odd (checkerboarded) QPhiX in/out spinors
  QSpinor *packed_spinor_in_cb0 = (QSpinor *)geom.allocCBFourSpinor();
  QSpinor *packed_spinor_in_cb1 = (QSpinor *)geom.allocCBFourSpinor();
  QSpinor *packed_spinor_out_cb0 = (QSpinor *)geom.allocCBFourSpinor();
  QSpinor *packed_spinor_out_cb1 = (QSpinor *)geom.allocCBFourSpinor();
  QSpinor *qphix_in[2];
  QSpinor *qphix_out[2];
  qphix_in[cb_even] = packed_spinor_in_cb0;
  qphix_in[cb_odd] = packed_spinor_in_cb1;
  qphix_out[cb_even] = packed_spinor_out_cb0;
  qphix_out[cb_odd] = packed_spinor_out_cb1;

  // Allocate data for odd (cb1) QPhiX prepared in spinor
  // and a buffer for the CG solver (to do the M^dagger matrix
  // multiplication after the solve)
  QSpinor *qphix_in_prepared = (QSpinor *)geom.allocCBFourSpinor();
  QSpinor *qphix_buffer = (QSpinor *)geom.allocCBFourSpinor();

  QPhiX::masterPrintf("# ...done.\n");

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
  QPhiX::masterPrintf("# Creating plain QPhiX Wilson Dslash (for preparation)...\n");
  QPhiX::Dslash<FT, V, S, compress> *WilsonDslash =
      new QPhiX::Dslash<FT, V, S, compress>(&geom, t_boundary, coeff_s, coeff_t);
  QPhiX::masterPrintf("# ...done.\n");

  // Create a Dslash & an even-odd preconditioned Fermion Matrix object,
  // depending on the chosen fermion action
  tmlqcd::Dslash<FT, V, S, compress> *DslashQPhiX;
  QPhiX::EvenOddLinearOperator<FT, V, S, compress> *FermionMatrixQPhiX;
  if (g_mu != 0.0 && g_c_sw > 0.0) {  // TWISTED-MASS-CLOVER
    // TODO: Implement me!
    QPhiX::masterPrintf("# TWISTED-MASS-CLOVER CASE NOT YET IMPLEMENTED!\n");
    QPhiX::masterPrintf(" Aborting...\n");
    abort();
  } else if (g_mu != 0.0) {  // TWISTED-MASS
    QPhiX::masterPrintf("# Creating QPhiX Twisted Mass Wilson Dslash...\n");
    const double TwistedMass = -g_mu / (2.0 * g_kappa);
    DslashQPhiX = new tmlqcd::WilsonTMDslash<FT, V, S, compress>(&geom, t_boundary, coeff_s,
                                                                 coeff_t, mass, TwistedMass);
    QPhiX::masterPrintf("# ...done.\n");

    QPhiX::masterPrintf("# Creating QPhiX Twisted Mass Wilson Fermion Matrix...\n");
    FermionMatrixQPhiX = new QPhiX::EvenOddTMWilsonOperator<FT, V, S, compress>(
        mass, TwistedMass, u_packed, &geom, t_boundary, coeff_s, coeff_t);
    QPhiX::masterPrintf("# ...done.\n");
  } else if (g_c_sw > 0.0) {  // WILSON CLOVER
    // TODO: Implement me!
    QPhiX::masterPrintf("# WILSON CLOVER CASE NOT YET IMPLEMENTED!\n");
    QPhiX::masterPrintf(" Aborting...\n");
    abort();
  } else {  // WILSON
    QPhiX::masterPrintf("# Creating QPhiX Wilson Dslash...\n");
    DslashQPhiX =
        new tmlqcd::WilsonDslash<FT, V, S, compress>(&geom, t_boundary, coeff_s, coeff_t, mass);
    QPhiX::masterPrintf("# ...done.\n");

    QPhiX::masterPrintf("# Creating QPhiX Wilson Fermion Matrix...\n");
    FermionMatrixQPhiX = new QPhiX::EvenOddWilsonOperator<FT, V, S, compress>(
        mass, u_packed, &geom, t_boundary, coeff_s, coeff_t);
    QPhiX::masterPrintf("# ...done.\n");
  }

  // Create a Linear Solver Object
  QPhiX::AbstractSolver<FT, V, S, compress> *SolverQPhiX;
  if (solver_flag == CG) {
    QPhiX::masterPrintf("# Creating CG Solver...\n");
    SolverQPhiX = new QPhiX::InvCG<FT, V, S, compress>(*FermionMatrixQPhiX, max_iter);
  } else if (solver_flag == BICGSTAB) {
    QPhiX::masterPrintf("# Creating BiCGStab Solver...\n");
    SolverQPhiX = new QPhiX::InvBiCGStab<FT, V, S, compress>(*FermionMatrixQPhiX, max_iter);
  } else {
    // TODO: Implement multi-shift CG, Richardson multi-precision
    QPhiX::masterPrintf(" Solver not yet supported by QPhiX!\n");
    QPhiX::masterPrintf(" Aborting...\n");
    abort();
  }
  QPhiX::masterPrintf("# ...done.\n");

  // Set number of BLAS threads by hand.
  // In case some implements the tune routines in QPhiX
  // this may be updated...
  QPhiX::masterPrintf("# Setting number of BLAS threads...\n");
  const int n_blas_simt = N_simt;
  QPhiX::masterPrintf("# ...done.\n");

  /************************
   *                      *
   *    PREPARE SOURCE    *
   *                      *
  ************************/

  QPhiX::masterPrintf("# Preparing odd source...\n");

  reorder_eo_spinor_to_QPhiX(geom, reinterpret_cast<double const *const>(tmlqcd_even_in),
                             qphix_in[cb_even], cb_even);
  reorder_eo_spinor_to_QPhiX(geom, reinterpret_cast<double const *const>(tmlqcd_odd_in),
                             qphix_in[cb_odd], cb_odd);

  // 2. Prepare the odd (cb1) source
  //
  //      \tilde b_o = 1/2 Dslash^{Wilson}_oe A^{-1}_{ee} b_e + b_o
  //
  // in three steps:
  // a) Apply A^{-1} to b_e and save result in qphix_buffer
  // b) Apply Wilson Dslash to qphix_buffer and save result in qphix_in_prepared
  // c) Apply AYPX to rescale last result (=y) and add b_o (=x)

  DslashQPhiX->A_inv_chi(qphix_buffer,       // out spinor
                         qphix_in[cb_even],  // in spinor
                         1);                 // non-conjugate
  WilsonDslash->dslash(qphix_in_prepared,    // out spinor
                       qphix_buffer,         // in spinor
                       u_packed[cb_odd],     // gauge field on target cb
                       1,                    // non-conjugate
                       cb_odd);              // target cb
  QPhiX::aypx(0.5, qphix_in[cb_odd], qphix_in_prepared, geom, n_blas_simt);

  QPhiX::masterPrintf("# ...done.\n");

  /************************
   *                      *
   *   SOLVE ON ODD CB    *
   *                      *
  ************************/

  QPhiX::masterPrintf("# Calling the solver...\n");

  // Set variables need for solve
  bool verbose = g_debug_level > 2 ? true : false;
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
  if (solver_flag == CG) {
    // USING CG:
    // We are solving
    //   M M^dagger qphix_buffer = qphix_in_prepared
    // here, that is, isign = -1 for the QPhiX CG solver.
    // After that multiply with M^dagger:
    //   qphix_out[1] = M^dagger M^dagger^-1 M^-1 qphix_in_prepared
    (*SolverQPhiX)(qphix_buffer, qphix_in_prepared, RsdTarget, niters, rsd_final, site_flops,
                   mv_apps, -1, verbose);
    (*FermionMatrixQPhiX)(qphix_out[cb_odd], qphix_buffer, /* conjugate */ -1);

  } else if (solver_flag == BICGSTAB) {
    // USING BiCGStab:
    // Solve M qphix_out[1] = qphix_in_prepared, directly.
    (*SolverQPhiX)(qphix_out[cb_odd], qphix_in_prepared, RsdTarget, niters, rsd_final, site_flops,
                   mv_apps, 1, verbose);
  }
  double end = omp_get_wtime();

  uint64_t num_cb_sites = lattSize[0] / 2 * lattSize[1] * lattSize[2] * lattSize[3];
  // FIXME: this needs to be adjusted depending on the operator used
  uint64_t total_flops = (site_flops + (72 + 2 * 1320) * mv_apps) * num_cb_sites;
  QPhiX::masterPrintf("# Solver Time = %g sec\n", (end - start));
  QPhiX::masterPrintf("# Performance in GFLOPS = %g\n", 1.0e-9 * total_flops / (end - start));

  /**************************
   *                        *
   *  RECONSTRUCT SOLUTION  *
   *                        *
  **************************/

  QPhiX::masterPrintf("# Reconstruction even solution...\n");

  // 1. Reconstruct the even (cb1) solution
  //
  //      x_e = A^{-1}_{ee} (b_e + 1/2 Dslash^{Wilson}_eo x_o)
  //
  // in three steps:
  // b) Apply Wilson Dslash to x_o and save result in qphix_buffer
  // c) Apply AYPX to rescale last result (=y) and add b_e (=x)
  // c) Apply A^{-1} to qphix_buffer and save result in x_e

  WilsonDslash->dslash(qphix_buffer,       // out spinor
                       qphix_out[cb_odd],  // in spinor (solution on odd cb)
                       u_packed[cb_even],  // gauge field on target cb
                       1,                  // non-conjugate
                       cb_even);           // target cb == even
  QPhiX::aypx(0.5, qphix_in[0], qphix_buffer, geom, n_blas_simt);
  DslashQPhiX->A_inv_chi(qphix_out[cb_even],  // out spinor
                         qphix_buffer,        // in spinor
                         1);                  // non-conjugate

  // 2. Reorder spinor fields back to tmLQCD, rescaling by a factor 1/(2*\kappa)
  //    to account for the operator normalisation in QPhiX

  reorder_eo_spinor_from_QPhiX(geom, reinterpret_cast<double *const>(tmlqcd_even_out),
                               qphix_out[cb_even], cb_even, 1.0 / (2.0 * g_kappa));
  reorder_eo_spinor_from_QPhiX(geom, reinterpret_cast<double *const>(tmlqcd_odd_out),
                               qphix_out[cb_odd], cb_odd, 1.0 / (2.0 * g_kappa));

  QPhiX::masterPrintf("# ...done.\n");

  /******************
   *                *
   *    CLEAN UP    *
   *                *
  ******************/

  QPhiX::masterPrintf("# Cleaning up\n");

  geom.free(packed_gauge_cb0);
  geom.free(packed_gauge_cb1);
  geom.free(packed_spinor_in_cb0);
  geom.free(packed_spinor_in_cb1);
  geom.free(packed_spinor_out_cb0);
  geom.free(packed_spinor_out_cb1);
  geom.free(qphix_in_prepared);
  geom.free(qphix_buffer);

  // FIXME: This should be called properly somewhere else
  _endQphix();

  QPhiX::masterPrintf("# ...done.\n\n");

  return niters;
}

// Template wrapper for the Dslash operator call-able from C code
void Mfull_qphix(spinor *Even_out, spinor *Odd_out, const spinor *Even_in, const spinor *Odd_in,
                 const op_type_t op_type) {
  tmlqcd::checkQphixInputParameters(qphix_input);
  // FIXME: two-row gauge compression and double precision hard-coded
  _initQphix(0, nullptr, qphix_input, 12, QPHIX_DOUBLE_PREC);

  if (qphix_precision == QPHIX_DOUBLE_PREC) {
    if (QPHIX_SOALEN > VECLEN_DP) {
      QPhiX::masterPrintf("SOALEN=%d is greater than the double prec VECLEN=%d\n", QPHIX_SOALEN,
                          VECLEN_DP);
      abort();
    }
    QPhiX::masterPrintf("TESTING IN DOUBLE PRECISION \n");
    if (compress12) {
      Mfull_helper<double, VECLEN_DP, QPHIX_SOALEN, true>(Even_out, Odd_out, Even_in, Odd_in,
                                                          op_type);
    } else {
      Mfull_helper<double, VECLEN_DP, QPHIX_SOALEN, false>(Even_out, Odd_out, Even_in, Odd_in,
                                                           op_type);
    }
  } else if (qphix_precision == QPHIX_FLOAT_PREC) {
    if (QPHIX_SOALEN > VECLEN_SP) {
      QPhiX::masterPrintf("SOALEN=%d is greater than the single prec VECLEN=%d\n", QPHIX_SOALEN,
                          VECLEN_SP);
      abort();
    }
    QPhiX::masterPrintf("TESTING IN SINGLE PRECISION \n");
    if (compress12) {
      Mfull_helper<float, VECLEN_SP, QPHIX_SOALEN, true>(Even_out, Odd_out, Even_in, Odd_in,
                                                         op_type);
    } else {
      Mfull_helper<float, VECLEN_SP, QPHIX_SOALEN, false>(Even_out, Odd_out, Even_in, Odd_in,
                                                          op_type);
    }
  }
#if (defined(QPHIX_MIC_SOURCE) || defined(QPHIX_AVX512_SOURCE))
  else if (qphix_precision == QPHIX_HALF_PREC) {
    if (QPHIX_SOALEN > VECLEN_HP) {
      QPhiX::masterPrintf("SOALEN=%d is greater than the half prec VECLEN=%d\n", QPHIX_SOALEN,
                          VECLEN_HP);
      abort();
    }
    QPhiX::masterPrintf("TESTING IN HALF PRECISION \n");
    if (compress12) {
      Mfull_helper<QPhiX::half, VECLEN_HP, QPHIX_SOALEN, true>(Even_out, Odd_out, Even_in, Odd_in,
                                                               op_type);
    } else {
      Mfull_helper<QPhiX::half, VECLEN_HP, QPHIX_SOALEN, false>(Even_out, Odd_out, Even_in, Odd_in,
                                                                op_type);
    }
  }
#endif
}

// Template wrapper for QPhiX solvers callable from C code, return number of iterations
int invert_eo_qphix(spinor *const Even_new, spinor *const Odd_new, spinor *const Even,
                    spinor *const Odd, const double precision, const int max_iter,
                    const int solver_flag, const int rel_prec, solver_params_t solver_params,
                    const SloppyPrecision sloppy, const CompressionType compression) {
  tmlqcd::checkQphixInputParameters(qphix_input);

  double target_precision = precision;
  double src_norm = (square_norm(Even, VOLUME / 2, 1) + square_norm(Odd, VOLUME / 2, 1));
  double precision_lambda = target_precision / src_norm;
  if (rel_prec == 1) {
    QPhiX::masterPrintf("# QPHIX: Using relative precision\n");
    target_precision = precision * src_norm;
    precision_lambda = precision;
  }
  QPhiX::masterPrintf("# QPHIX: precision_lambda: %g, target_precision: %g\n\n", precision_lambda,
                      target_precision);

#if (defined(QPHIX_MIC_SOURCE) || defined(QPHIX_AVX512_SOURCE))
  if (sloppy == SLOPPY_HALF || precision_lambda >= rsdTarget<QPhiX::half>::value) {
    if (QPHIX_SOALEN > VECLEN_HP) {
      QPhiX::masterPrintf("SOALEN=%d is greater than the half prec VECLEN=%d\n", QPHIX_SOALEN,
                          VECLEN_HP);
      abort();
    }
    QPhiX::masterPrintf("# INITIALIZING QPHIX SOLVER\n");
    QPhiX::masterPrintf("# USING HALF PRECISION\n");
    _initQphix(0, nullptr, qphix_input, compression, QPHIX_HALF_PREC);

    if (compress12) {
      return invert_eo_qphix_helper<QPhiX::half, VECLEN_HP, QPHIX_SOALEN, true>(
          Even_new, Odd_new, Even, Odd, target_precision, max_iter, solver_flag, rel_prec,
          solver_params, compression);
    } else {
      return invert_eo_qphix_helper<QPhiX::half, VECLEN_HP, QPHIX_SOALEN, false>(
          Even_new, Odd_new, Even, Odd, target_precision, max_iter, solver_flag, rel_prec,
          solver_params, compression);
    }
  }
#else
  if (sloppy == SLOPPY_HALF) {
    QPhiX::masterPrintf("QPHIX interface: half precision not supported on this architecture!\n");
    abort();
  } else
#endif
  if (sloppy == SLOPPY_SINGLE || precision_lambda >= rsdTarget<float>::value) {
    if (QPHIX_SOALEN > VECLEN_SP) {
      QPhiX::masterPrintf("SOALEN=%d is greater than the single prec VECLEN=%d\n", QPHIX_SOALEN,
                          VECLEN_SP);
      abort();
    }
    QPhiX::masterPrintf("# INITIALIZING QPHIX SOLVER\n");
    QPhiX::masterPrintf("# USING SINGLE PRECISION\n");
    _initQphix(0, nullptr, qphix_input, compression, QPHIX_FLOAT_PREC);

    if (compress12) {
      return invert_eo_qphix_helper<float, VECLEN_SP, QPHIX_SOALEN, true>(
          Even_new, Odd_new, Even, Odd, target_precision, max_iter, solver_flag, rel_prec,
          solver_params, compression);
    } else {
      return invert_eo_qphix_helper<float, VECLEN_SP, QPHIX_SOALEN, false>(
          Even_new, Odd_new, Even, Odd, target_precision, max_iter, solver_flag, rel_prec,
          solver_params, compression);
    }
  } else {
    if (QPHIX_SOALEN > VECLEN_DP) {
      QPhiX::masterPrintf("SOALEN=%d is greater than the double prec VECLEN=%d\n", QPHIX_SOALEN,
                          VECLEN_DP);
      abort();
    }
    QPhiX::masterPrintf("# INITIALIZING QPHIX SOLVER\n");
    QPhiX::masterPrintf("# USING DOUBLE PRECISION\n");
    _initQphix(0, nullptr, qphix_input, compression, QPHIX_DOUBLE_PREC);

    if (compress12) {
      return invert_eo_qphix_helper<double, VECLEN_DP, QPHIX_SOALEN, true>(
          Even_new, Odd_new, Even, Odd, target_precision, max_iter, solver_flag, rel_prec,
          solver_params, compression);
    } else {
      return invert_eo_qphix_helper<double, VECLEN_DP, QPHIX_SOALEN, false>(
          Even_new, Odd_new, Even, Odd, target_precision, max_iter, solver_flag, rel_prec,
          solver_params, compression);
    }
  }  // if( sloppy || target_precision )
  return -1;
}

void tmlqcd::checkQphixInputParameters(const QphixParams_t &params) {
  if (params.MinCt == 0) {
    QPhiX::masterPrintf("QPHIX Error: MinCt cannot be 0! Minimal value: 1. Aborting.\n");
    abort();
  }
  if (params.By == 0 || params.Bz == 0) {
    QPhiX::masterPrintf("QPHIX Error: By and Bz may not be 0! Minimal value: 1. Aborting.\n");
    abort();
  }
  if (params.NCores * params.Sy * params.Sz != omp_num_threads) {
    QPhiX::masterPrintf("QPHIX Error: NCores * Sy * Sz != ompnumthreads ! Aborting.\n");
    abort();
  }
}

void tmlqcd::printQphixDiagnostics(int VECLEN, int SOALEN, bool compress) {
  QPhiX::masterPrintf("# QphiX: VECLEN=%d SOALEN=%d\n", VECLEN, SOALEN);
  
  QPhiX::masterPrintf("# QphiX: Declared QMP Topology (xyzt):");
  for(int mu = 0; mu < 4; mu++) QPhiX::masterPrintf(" %d", qmp_geom[mu]);
  QPhiX::masterPrintf("\n");
  
  QPhiX::masterPrintf("# QphiX: Mapping of dimensions QMP -> tmLQCD (xyzt):");
  for(int mu = 0; mu < 4; mu++) QPhiX::masterPrintf(" %d->%d", mu, qmp_tm_map[mu] );
  QPhiX::masterPrintf("\n");

  QPhiX::masterPrintf("# QphiX: Global Lattice Size (xyzt) = ");
  for (int mu = 0; mu < 4; mu++) {
    QPhiX::masterPrintf(" %d", lattSize[mu]);
  }
  QPhiX::masterPrintf("\n");
  QPhiX::masterPrintf("# QphiX: Local Lattice Size (xyzt) = ");
  for (int mu = 0; mu < 4; mu++) {
    QPhiX::masterPrintf(" %d", subLattSize[mu]);
  }
  QPhiX::masterPrintf("\n");
  QPhiX::masterPrintf("# QphiX: Block Sizes: By= %d Bz=%d\n", By, Bz);
  QPhiX::masterPrintf("# QphiX: Cores = %d\n", NCores);
  QPhiX::masterPrintf("# QphiX: SMT Grid: Sy=%d Sz=%d\n", Sy, Sz);
  QPhiX::masterPrintf("# QphiX: Pad Factors: PadXY=%d PadXYZ=%d\n", PadXY, PadXYZ);
  QPhiX::masterPrintf("# QphiX: Threads_per_core = %d\n", N_simt);
  QPhiX::masterPrintf("# QphiX: MinCt = %d\n", MinCt);
  if (compress) {
    QPhiX::masterPrintf("# QphiX: Using two-row gauge compression (compress12)\n");
  }
}

void testSpinorPackers(spinor *Even_out, spinor *Odd_out, const spinor *const Even_in,
                       const spinor *const Odd_in) {
  tmlqcd::checkQphixInputParameters(qphix_input);
  // FIXME: two-row gauge compression and double precision hard-coded
  _initQphix(0, nullptr, qphix_input, 12, QPHIX_DOUBLE_PREC);

  QPhiX::Geometry<double, VECLEN_SP, QPHIX_SOALEN, true> geom(subLattSize, By, Bz, NCores, Sy, Sz,
                                                              PadXY, PadXYZ, MinCt);

  auto qphix_cb_even = QPhiX::makeFourSpinorHandle(geom);
  auto qphix_cb_odd = QPhiX::makeFourSpinorHandle(geom);

  spinor **tmp;
  init_solver_field(&tmp, VOLUME / 2, 2);

  reorder_eo_spinor_to_QPhiX(geom, reinterpret_cast<double const *const>(Even_in),
                             qphix_cb_even.get(), cb_even);
  reorder_eo_spinor_to_QPhiX(geom, reinterpret_cast<double const *const>(Odd_in),
                             qphix_cb_odd.get(), cb_odd);

  reorder_eo_spinor_from_QPhiX(geom, reinterpret_cast<double *>(Even_out), qphix_cb_even.get(),
                               cb_even, 1.0);
  reorder_eo_spinor_from_QPhiX(geom, reinterpret_cast<double *>(Odd_out), qphix_cb_odd.get(),
                               cb_odd, 1.0);

  diff(tmp[0], Even_out, Even_in, VOLUME / 2);
  diff(tmp[1], Odd_out, Odd_in, VOLUME / 2);
  double l2norm = square_norm(tmp[0], VOLUME / 2, 1) + square_norm(tmp[1], VOLUME / 2, 1);
  QPhiX::masterPrintf("QPHIX eo spinor packer back and forth difference L2 norm: %lf\n", l2norm);
  finalize_solver(tmp, 2);
}

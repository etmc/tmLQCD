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
#include "qphix_interface.hpp"
#include "qphix_interface_utils.hpp"
#include "qphix_types.h"
#include "qphix_veclen.h"

#ifdef TM_USE_MPI
#include <mpi.h>
#endif

extern "C" {
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "boundary.h"
#include "geometry_eo.h"
#include "gettime.h"
#include "global.h"
#include "struct_accessors.h"
#include "linalg/convert_eo_to_lexic.h"
#include "linalg/diff.h"
#include "linalg/square_norm.h"
#include "linalg/square_norm.h"
#include "misc_types.h"
// for the normalisation of the heavy doublet when running
// RHMC
#include "phmc.h"
#include "start.h"
#include "operator/clover_leaf.h"
#include "operator/clovertm_operators.h"
#include "operator_types.h"
#include "operator/Hopping_Matrix.h"
#include "solver/matrix_mult_typedef.h"
#include "solver/solver_types.h"
#include "solver/solver.h"
#include "solver/solver_field.h"
#include "solver/solver_params.h"
#include "xchange/xchange_gauge.h"
}
#ifdef TM_USE_OMP
#include <omp.h>
#endif
#include <qphix/blas_new_c.h>
#include <qphix/clover.h>
#include <qphix/invbicgstab.h>
#include <qphix/invcg.h>
#include <qphix/inv_richardson_multiprec.h>
#include <qphix/inv_dummy_hermtest.h>
#include <qphix/minvcg.h>
#include <qphix/ndtm_reuse_operator.h>
#include <qphix/ndtm_reuse_operator_clover.h>
#include <qphix/print_utils.h>
#include <qphix/qphix_config.h>
#include <qphix/twisted_mass.h>
#include <qphix/twisted_mass_clover.h>
#include <qphix/wilson.h>
#include <cfloat>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <complex>
#include <cmath>

using namespace tmlqcd;

tm_QPhiXParams_t qphix_input;

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
QphixPrec_t qphix_inner_precision;


int subLattSize[4];
int lattSize[4];
int qmp_geom[4];
int qmp_tm_map[4];

// angles for boundary phases, values come from read_input
extern double X0, X1, X2, X3;

bool use_tbc[4];
double tbc_phases[4][2];
// we always use twisted boundary conditions, which means that we are always
// periodic in time and any possible anti-periodicity is implemented via
// the phase
double constexpr t_boundary = 1.0;

template <typename T>
struct rsdTarget {
  static const double value;
};

template <>
const double rsdTarget<QPhiX::half>::value = 1.0e-3;

template <>
const double rsdTarget<float>::value = 1.0e-8;


void _initQphix(int argc, char **argv, tm_QPhiXParams_t params, int c12, QphixPrec_t precision_, QphixPrec_t inner_precision_) {
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

  // extract twisted boundary conditions
  for (int dim = 0; dim < 4; dim++) {
    bool dim_tbc = false;
    double dim_phase[2] = {1.0, 0.0};
    if (dim == 0) {
      dim_tbc = (fabs(X1) > DBL_EPSILON);
      dim_phase[0] = creal(-phase_1 / g_kappa);
      dim_phase[1] = cimag(-phase_1 / g_kappa);
    } else if (dim == 1) {
      dim_tbc = (fabs(X2) > DBL_EPSILON);
      dim_phase[0] = creal(-phase_2 / g_kappa);
      dim_phase[1] = cimag(-phase_2 / g_kappa);
    } else if (dim == 2) {
      dim_tbc = (fabs(X3) > DBL_EPSILON);
      dim_phase[0] = creal(-phase_3 / g_kappa);
      dim_phase[1] = cimag(-phase_3 / g_kappa);
    } else if (dim == 3) {
      dim_tbc = (fabs(X0) > DBL_EPSILON);
      dim_phase[0] = creal(-phase_0 / g_kappa);
      dim_phase[1] = cimag(-phase_0 / g_kappa);
    }
    use_tbc[dim] = dim_tbc;
    tbc_phases[dim][0] = dim_phase[0];
    tbc_phases[dim][1] = dim_phase[1];
  }

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
  qphix_inner_precision = inner_precision_;

#ifdef QPHIX_QMP_COMMS
  // Declare the logical topology
  if (!qmp_topo_initialised) {
    // the QMP topology is the one implied by the number of processes in each
    // dimension as required by QPHIX ( x fastest to t slowest running )
    qmp_geom[0] = g_nproc_x;
    qmp_geom[1] = g_nproc_y;
    qmp_geom[2] = g_nproc_z;
    qmp_geom[3] = g_nproc_t;

    // in order for the topologies to agree between tmLQCD and QPhiX, the dimensions need to be
    // permuted
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
    if (g_debug_level >= 5) {
      for (int proc = 0; proc < g_nproc; proc++) {
        if (proc == g_proc_id) {
          const int coordinates[4] = {g_proc_coords[1], g_proc_coords[2], g_proc_coords[3],
                                      g_proc_coords[0]};
          int id = QMP_get_node_number_from(coordinates);
          int *qmp_coords = QMP_get_logical_coordinates_from(id);
          fflush(stdout);
          printf("QMP id: %3d x:%3d y:%3d z:%3d t:%3d\n", id, qmp_coords[0], qmp_coords[1],
                 qmp_coords[2], qmp_coords[3]);
          printf("MPI id: %3d x:%3d y:%3d z:%3d t:%3d\n\n", g_proc_id, g_proc_coords[1],
                 g_proc_coords[2], g_proc_coords[3], g_proc_coords[0]);
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

void _initQphix(int argc, char **argv, tm_QPhiXParams_t params, int c12, QphixPrec_t precision_){
  _initQphix(argc, argv, params, c12, precision_, precision_);
}

// Finalize the QPhiX library
void _endQphix() {}

template <typename FT, int VECLEN, int SOALEN, bool compress12>
void reorder_clover_to_QPhiX(
    QPhiX::Geometry<FT, VECLEN, SOALEN, compress12> &geom,
    typename QPhiX::Geometry<FT, VECLEN, SOALEN, compress12>::CloverBlock *qphix_clover, int cb,
    bool inverse, bool fl_offdiag = false) {
  const double startTime = gettime();

  /* the spin-colour clover term in sw_term and the corresponding inverse
   * in sw_inv are stored in the tmLQCD gamma basis.
   * When we translate spinors to QPhiX, we apply a transformation V to the tmLQCD
   * spinor and then apply the same transformation to the output spinor
   * ( we have V^dagger = V and V*V = 1 )
   * Thus, in order to translate the clover field, we need to copy
   *   (1+T)' = V*(1+T)*V, where T is the spin-colour clover-term
   * This way, the clover term will be in the correct gamma basis.
   *
   * The tmLQCD clover term is stored in half-spinor blocks of colour matrices
   * for which we need to work out what (1+T)'=V*(1+T)*V implies.
   * Below, each sAB represents one 3x3 colour matrix
   *
   *                +s33 -s32    0    0
   *  T' = V*T*V =  -s23 +s22    0    0
   *                   0    0 +s11 -s10
   *                   0    0 -s01 +s00
   *
   * Such that the half-spinor blocks are inverted and within these, the ordering is
   * reversed. Note that the off-diagonal 3x3 colour blocks are hermitian conjugate to
   * each other and this is preserved by the transformation.
   *
   * The QPhiX (Wilson) clover term is stored as 12 reals on the diagonal
   * in two 6-element vectors, one for each half-spinor spin pair
   * and two sets of off-diagonal complex components.
   *
   * In addition, colour matrices are transposed in QPhiX.
   *
   * The tmLQCD clover term is stored as:
   *
   *      s00 s01
   *          s11
   * T =          s22 s23
   *                  s33
   *
   * with indexing
   *
   *     sw[0][0] sw[1][0]
   *              sw[2][0]
   *                       sw[0][1] sw[1][1]
   *                                sw[2][1]
   *
   * The inverse has four su3 blocks instead and is indexed
   *     sw_inv[0][0] sw_inv[1][0]
   *     sw_inv[3][0] sw_inv[2][0]
   *                               sw_inv[0][1] sw_inv[1][1]
   *                               sw_inv[3][1] sw_inv[2][1]
   *
   * where blocks sw_inv[3][0] and sw_inv[3][1] are relevant only when mu > 0
   * 
   * There is a special case for the non-degenerate twisted clover operator. The flavour-off-diagonal
   * components of the inverse clover term do not have an imaginary part on the spin-colour diagonal.
   * They can thus be stored as CloverBlock, which is done in the QPhiX implementation
   * of the ND tmclover operator.
   * 
   * As a hack, this inverse is prepared by sw_invert_epsbar and placed in to the last
   * VOLUME/2 sites of sw_inv. Reading from there is triggered by the boolean
   * fl_offdiag.
   */

  // rescale to get clover term (or its inverse) in the physical normalisation
  // rather than the kappa normalisation
  const double scale = inverse ? 2.0 * g_kappa : 1.0 / (2.0 * g_kappa);
  su3 ***tm_clover = inverse ? sw_inv : sw;

  // Number of elements in spin, color & complex
  const int Ns = 4;
  const int Nc = 3;
  const int Nz = 2;

  // Geometric parameters for QPhiX data layout
  const auto ngy = geom.nGY();
  const auto nVecs = geom.nVecs();
  const auto Pxy = geom.getPxy();
  const auto Pxyz = geom.getPxyz();

  // packer for Wilson clover (real diagonal + complex upper-triangular)
  /* for the index in the off_diagN arrays, we map to an index in the su3 struct
  * keeping in mind complex conjugation
  * The off-diagonal in QPhiX is stored as follows:
  *
  * 0 1 3 6 10
  *   2 4 7 11
  *     5 8 12
  *       9 13
  *         14
  *
  * which we are going to map to su3 in blocks
  *
  *     0* 1*
  *        2*
  *
  * 3   4  5
  * 6   7  8
  * 10 11 12
  *
  *   9* 13*
  *      14*
  *
  * where the asterisk indicates complex conjugation. As a linear array then,
  * these mappings are:
  *
  */
  const int od_su3_offsets[15] = {Nz,
                                  2 * Nz,            //     0 1
                                  Nc * Nz + 2 * Nz,  //       2

                                  0,
                                  Nz,
                                  2 * Nz,  // 3  4  5
                                  Nc * Nz,
                                  Nc * Nz + Nz,
                                  Nc * Nz + 2 * Nz,  // 6  7  8

                                  Nz,  //     9

                                  2 * Nc * Nz,
                                  2 * Nc * Nz + Nz,
                                  2 * Nc * Nz + 2 * Nz,  // 10 11 12

                                  2 * Nz,
                                  Nc * Nz + 2 * Nz};  // 13 14

#pragma omp parallel for collapse(4)
  for (int64_t t = 0; t < T; t++) {
    for (int64_t z = 0; z < LZ; z++) {
      for (int64_t y = 0; y < LY; y++) {
        for (int64_t v = 0; v < nVecs; v++) {
          int64_t block = (t * Pxyz + z * Pxy) / ngy + (y / ngy) * nVecs + v;

          for (int64_t x_soa = 0; x_soa < SOALEN; x_soa++) {
            int64_t xx = (y % ngy) * SOALEN + x_soa;
            int64_t q_cb_x_coord = x_soa + v * SOALEN;
            int64_t tm_x_coord = q_cb_x_coord * 2 + (((t + y + z) & 1) ^ cb);

            // the inverse of the clover term is in even-odd ordering
            // while the clover term itself is lexicographically ordered
            // for the special case of the nd tmclover operator, the inverse of the flavour off-diagonal
            // components is stored in the last VOLUME/2 elements of sw_inv
            int64_t tm_idx =
                (inverse ? g_lexic2eosub[g_ipt[t][tm_x_coord][y][z]] : g_ipt[t][tm_x_coord][y][z]) +
                ( (inverse && fl_offdiag) ? VOLUME/2 : 0 );

            int b_idx;

            //             we begin with the diagonal elements in CloverBlock
            for (int d = 0; d < 6; d++) {
              //               choose the block in sw which corresponds to the block in T'
              b_idx = d < 3 ? 2 : 0;
              //               get the right colour components
              qphix_clover[block].diag1[d][xx] =
                QPhiX::rep<FT, double>(
                  *(reinterpret_cast<double const *const>(&tm_clover[tm_idx][b_idx][1].c00) +
                    (Nc * Nz + Nz) * (d % 3)) *
                  scale
                );

              qphix_clover[block].diag2[d][xx] =
                QPhiX::rep<FT, double>(
                  *(reinterpret_cast<double const *const>(&tm_clover[tm_idx][b_idx][0].c00) +
                    (Nc * Nz + Nz) * (d % 3)) *
                  scale
                );
            }

            b_idx = 2;  // s33 and s11
            for (int od : {0, 1, 2}) {
              for (int reim : {0, 1}) {
                qphix_clover[block].off_diag1[od][reim][xx] =
                  QPhiX::rep<FT, double>(
                    (reim == 1 ? -1.0 : 1.0) *
                    *(reinterpret_cast<double const *const>(&tm_clover[tm_idx][b_idx][1].c00) +
                      od_su3_offsets[od] + reim) *
                    scale
                  );

                qphix_clover[block].off_diag2[od][reim][xx] =
                  QPhiX::rep<FT, double>( 
                    (reim == 1 ? -1.0 : 1.0) *
                    *(reinterpret_cast<double const *const>(&tm_clover[tm_idx][b_idx][0].c00) +
                      od_su3_offsets[od] + reim) *
                    scale 
                  ); 
              }
            }

            b_idx = 1;  // s32 and s10
            for (int od : {3, 4, 5, 6, 7, 8, 10, 11, 12}) {
              for (int reim : {0, 1}) {
                qphix_clover[block].off_diag1[od][reim][xx] = 
                  QPhiX::rep<FT, double>(
                    *(reinterpret_cast<double const *const>(&tm_clover[tm_idx][b_idx][1].c00) +
                      od_su3_offsets[od] + reim) *
                    (-scale) 
                  );

                qphix_clover[block].off_diag2[od][reim][xx] =
                  QPhiX::rep<FT, double>(
                    *(reinterpret_cast<double const *const>(&tm_clover[tm_idx][b_idx][0].c00) +
                      od_su3_offsets[od] + reim) *
                    (-scale)
                  );
              }
            }

            b_idx = 0;  // s22 and s00
            for (int od : {9, 13, 14}) {
              for (int reim : {0, 1}) {
                qphix_clover[block].off_diag1[od][reim][xx] =
                  QPhiX::rep<FT, double>(
                    (reim == 1 ? -1.0 : 1.0) *
                    *(reinterpret_cast<double const *const>(&tm_clover[tm_idx][b_idx][1].c00) +
                      od_su3_offsets[od] + reim) *
                    scale
                  );

                qphix_clover[block].off_diag2[od][reim][xx] =
                  QPhiX::rep<FT, double>(
                    (reim == 1 ? -1.0 : 1.0) *
                    *(reinterpret_cast<double const *const>(&tm_clover[tm_idx][b_idx][0].c00) +
                      od_su3_offsets[od] + reim) *
                    scale
                  );
              }
            }

          }  // x_soa
        }    // for(v)
      }      // for(y)
    }        // for(z)
  }          // for(t)

  const double diffTime = gettime() - startTime;
  if (g_debug_level > 1) {
    QPhiX::masterPrintf(
        "# QPHIX-interface: time spent in reorder_clover_to_QPhiX (CloverBlock): %f secs\n",
        diffTime);
  }
}

template <typename FT, int VECLEN, int SOALEN, bool compress12>
void reorder_clover_to_QPhiX(
    QPhiX::Geometry<FT, VECLEN, SOALEN, compress12> &geom,
    typename QPhiX::Geometry<FT, VECLEN, SOALEN, compress12>::FullCloverBlock *qphix_clover[2],
    int cb, bool inverse) {
  const double startTime = gettime();

  /* the spin-colour clover term in sw_term and the corresponding inverse
   * in sw_inv are stored in the tmLQCD gamma basis.
   * When we translate spinors to QPhiX, we apply a transformation V to the tmLQCD
   * spinor and then apply the same transformation to the output spinor
   * ( we have V^dagger = V and V*V = 1 )
   * Thus, in order to translate the clover field, we need to copy
   *   (1+T)' = V*(1+T)*V, where T is the spin-colour clover-term
   * This way, the clover term will be in the correct gamma basis.
   *
   * The tmLQCD clover term is stored in half-spinor blocks of colour matrices
   * for which we need to work out what (1+T)'=V*(1+T)*V implies.
   * Below, each sAB represents one 3x3 colour matrix
   *
   *                +s33 -s32    0    0
   *  T' = V*T*V =  -s23 +s22    0    0
   *                   0    0 +s11 -s10
   *                   0    0 -s01 +s00
   *
   * Such that the half-spinor blocks are inverted and within these, the ordering is
   * reversed. Note that the off-diagonal 3x3 colour blocks are hermitian conjugate to
   * each other and this is preserved by the transformation.
   *
   * The QPhiX (tmclover) clover term and its inverse are stored as a pair of full
   * 6x6 complex matrices which are multiplied with the spinor in exactly the same way
   * as in tmLQCD.
   *
   * The tmLQCD clover term is stored as:
   *
   *      s00 s01
   *          s11
   * T =          s22 s23
   *                  s33
   *
   * with indexing
   *
   *     sw[0][0] sw[1][0]
   *              sw[2][0]
   *                       sw[0][1] sw[1][1]
   *                                sw[2][1]
   *
   * The inverse has four su3 blocks instead and is indexed
   *     sw_inv[0][0] sw_inv[1][0]
   *     sw_inv[3][0] sw_inv[2][0]
   *                               sw_inv[0][1] sw_inv[1][1]
   *                               sw_inv[3][1] sw_inv[2][1]
   *
   * where blocks sw_inv[3][0] and sw_inv[3][1] are relevant only when mu > 0   *
   */

  // rescale to get clover term (or its inverse) in the physical normalisation
  // rather than the kappa normalisation
  const double scale = inverse ? 2.0 * g_kappa : 1.0 / (2.0 * g_kappa);
  su3 ***tm_clover = inverse ? sw_inv : sw;

  // Number of elements in spin, color & complex
  const int Ns = 4;
  const int Nc = 3;
  const int Nz = 2;

  const double amu = g_mu / (2.0 * g_kappa);

  // Geometric parameters for QPhiX data layout
  const auto ngy = geom.nGY();
  const auto nVecs = geom.nVecs();
  const auto Pxy = geom.getPxy();
  const auto Pxyz = geom.getPxyz();

#pragma omp parallel for collapse(4)
  for (int64_t t = 0; t < T; t++) {
    for (int64_t z = 0; z < LZ; z++) {
      for (int64_t y = 0; y < LY; y++) {
        for (int64_t v = 0; v < nVecs; v++) {
          int64_t block = (t * Pxyz + z * Pxy) / ngy + (y / ngy) * nVecs + v;

          for (int64_t x_soa = 0; x_soa < SOALEN; x_soa++) {
            int64_t xx = (y % ngy) * SOALEN + x_soa;
            int64_t q_cb_x_coord = x_soa + v * SOALEN;
            int64_t tm_x_coord = q_cb_x_coord * 2 + (((t + y + z) & 1) ^ cb);

            //             the inverse of the clover term is in even-odd ordering
            //             while the clover term itself is lexicographically ordered
            int64_t tm_idx =
                inverse ? g_lexic2eosub[g_ipt[t][tm_x_coord][y][z]] : g_ipt[t][tm_x_coord][y][z];

            for (int fl : {0, 1}) {
              if (inverse && fl == 1) {
                // the inverse clover term for the second flavour is stored at an offset
                tm_idx += VOLUME / 2;
              }
              for (int q_hs : {0, 1}) {
                auto &hs_block =
                    ((q_hs == 0) ? qphix_clover[fl][block].block1 : qphix_clover[fl][block].block2);
                for (int q_sc1 = 0; q_sc1 < 6; q_sc1++) {
                  for (int q_sc2 = 0; q_sc2 < 6; q_sc2++) {
                    const int q_s1 = q_sc1 / 3;
                    const int q_s2 = q_sc2 / 3;
                    const int q_c1 = q_sc1 % 3;
                    const int q_c2 = q_sc2 % 3;

                    // invert in spin as required by V*T*V
                    const int t_hs = 1 - q_hs;
                    // the indices inside the half-spinor are also inverted
                    // (which transposes them, of course)
                    const int t_s1 = 1 - q_s1;
                    const int t_s2 = 1 - q_s2;
                    // carry out the mapping from T' to T, keeping in mind that for the inverse
                    // there are four blocks also on the tmLQCD side, otherwise there are just three
                    const int t_b_idx = t_s1 + t_s2 + ((inverse && t_s1 == 1 && t_s2 == 0) ? 2 : 0);
                    for (int reim : {0, 1}) {
                      hs_block[q_sc1][q_sc2][reim][xx] =
                        QPhiX::rep<FT,double>(
                          scale *
                              // off-diagonal (odd-numbered) blocks change sign
                              (t_b_idx & 1 ? (-1.0) : 1.0) *
                              // if not doing the inverse and in the bottom-left block, need to
                              // complex conjugate
                              ((!inverse && (t_s1 == 1 && t_s2 == 0) && reim == 1) ? -1.0 : 1.0) *
                              *(reinterpret_cast<double const *const>(
                                    &(tm_clover[tm_idx][t_b_idx][t_hs].c00)) +
                                // if not doing the inverse and in the bottom-left block, transpose
                                // in colour
                                // because we're actually reading out of the top-right block
                                Nz * ((!inverse && (t_s1 == 1 && t_s2 == 0)) ? Nc * q_c2 + q_c1
                                                                             : Nc * q_c1 + q_c2) +
                                reim) +
                          // in the QPhiX gamma basis, the twisted quark mass enters with the
                          // opposite
                          // sign for consistency
                          ((!inverse && q_sc1 == q_sc2 && q_hs == 0 && reim == 1)
                               ? -amu * (1 - 2 * fl)
                               : 0) +
                          ((!inverse && q_sc1 == q_sc2 && q_hs == 1 && reim == 1)
                               ? amu * (1 - 2 * fl)
                               : 0)
                        );
                    }
                  }  // q_sc2
                }    // q_sc1
              }      // q_hs
            }        // fl

          }  // x_soa
        }    // for(v)
      }      // for(y)
    }        // for(z)
  }          // for(t)

  const double diffTime = gettime() - startTime;
  if (g_debug_level > 1) {
    QPhiX::masterPrintf(
        "# QPHIX-interface: time spent in reorder_clover_to_QPhiX (FullCloverBlock): %f secs\n",
        diffTime);
  }
}

template <typename FT, int VECLEN, int SOALEN, bool compress12>
void reorder_gauge_to_QPhiX(
    QPhiX::Geometry<FT, VECLEN, SOALEN, compress12> &geom,
    typename QPhiX::Geometry<FT, VECLEN, SOALEN, compress12>::SU3MatrixBlock *qphix_gauge_cb0,
    typename QPhiX::Geometry<FT, VECLEN, SOALEN, compress12>::SU3MatrixBlock *qphix_gauge_cb1) {
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
        for (int64_t v = 0; v < nVecs; v++) {
          int64_t block = (t * Pxyz + z * Pxy) / ngy + (y / ngy) * nVecs + v;

          for (int dim = 0; dim < 4; dim++)     // dimension == QPhiX \mu
            for (int c1 = 0; c1 < Nc1; c1++)    // QPhiX convention color 1 (runs up to 2 or 3)
              for (int c2 = 0; c2 < Nc2; c2++)  // QPhiX convention color 2 (always runs up to 3)
                for (int x_soa = 0; x_soa < SOALEN; x_soa++) {
                  int64_t xx = (y % ngy) * SOALEN + x_soa;
                  int64_t q_cb_x_coord = x_soa + v * SOALEN;
                  int64_t tm_x_coord_cb0 = q_cb_x_coord * 2 + (((t + y + z) & 1) ^ 0);
                  int64_t tm_x_coord_cb1 = q_cb_x_coord * 2 + (((t + y + z) & 1) ^ 1);

                  int64_t tm_idx_cb0;
                  int64_t tm_idx_cb1;

                  // backward / forward
                  for (int dir = 0; dir < 2; dir++) {
                    if (dir == 0) {
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

                      qphix_gauge_cb0[block][q_mu][c1][c2][reim][xx] = QPhiX::rep<FT, double>(
                        su3_get_elem(&(g_gauge_field[tm_idx_cb0][change_dim[dim]]), c2, c1, reim ) );
                      qphix_gauge_cb1[block][q_mu][c1][c2][reim][xx] = QPhiX::rep<FT, double>(
                        su3_get_elem(&(g_gauge_field[tm_idx_cb1][change_dim[dim]]), c2, c1, reim ) );
                    }
                  }
                }  // for(dim,c1,c2,x_soa)
        }          // outer loop (t,z,y,v)

  const double diffTime = gettime() - startTime;
  if (g_debug_level > 1) {
    QPhiX::masterPrintf("# QPHIX-interface: time spent in reorder_gauge_to_QPhiX: %f secs\n",
                        diffTime);
  }
}

// Reorder tmLQCD eo-spinor to a FourSpinorBlock QPhiX spinor on the given checkerboard
template <typename FT, int VECLEN, int SOALEN, bool compress12>
void reorder_eo_spinor_to_QPhiX(
    QPhiX::Geometry<FT, VECLEN, SOALEN, compress12> &geom, spinor const *const tm_eo_spinor,
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

                for (int reim = 0; reim < 2; reim++) {
                  qphix_spinor[q_ind][col][q_spin][reim][x_soa] =
                    QPhiX::rep<FT, double>(
                         change_sign[q_spin] * spinor_get_elem( &(tm_eo_spinor[tm_eo_ind]), 
                                                                change_spin[q_spin],
                                                                col,
                                                                reim
                                                              )
                    );
                }
              }
            }
          }
        }
      }
    }
  }
  const double diffTime = gettime() - startTime;
  if (g_debug_level > 1) {
    QPhiX::masterPrintf("# QPHIX-interface: time spent in reorder_eo_spinor_to_QPhiX: %f secs\n",
                        diffTime);
  }
}

template <typename FT, int VECLEN, int SOALEN, bool compress12>
void reorder_eo_spinor_from_QPhiX(
    QPhiX::Geometry<FT, VECLEN, SOALEN, compress12> &geom, spinor* tm_eo_spinor,
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

                spinor_set_elem( &(tm_eo_spinor[tm_eo_ind]),
                                 change_spin[q_spin],
                                 col,
                                 change_sign[q_spin] * normFac * QPhiX::rep<double, FT>(
                                  qphix_spinor[q_ind][col][q_spin][0][x_soa]
                                 ),
                                 change_sign[q_spin] * normFac * QPhiX::rep<double, FT>(
                                  qphix_spinor[q_ind][col][q_spin][1][x_soa]
                                 )
                               );
              }
            }
          }
        }
      }
    }
  }
  const double diffTime = gettime() - startTime;
  if (g_debug_level > 1) {
    QPhiX::masterPrintf("# QPHIX-interface: time spent in reorder_eo_spinor_from_QPhiX: %f secs\n",
                        diffTime);
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

                out[qId] = QPhiX::rep<FT, double>( change_sign[spin] * in[tId] ); 
              }

        }  // volume

  const double diffTime = gettime() - startTime;
  if (g_debug_level > 1) {
    QPhiX::masterPrintf("# QPHIX-interface: time spent in reorder_spinor_to_QPhiX: %f secs\n",
                        diffTime);
  }
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

                out[tId] = QPhiX::rep<double, FT>( normFac * change_sign[spin] * in[qId] );
              }

        }  // volume

  const double diffTime = gettime() - startTime;
  if (g_debug_level > 1) {
    QPhiX::masterPrintf("# QPHIX-interface: time spent in reorder_spinor_from_QPhiX: %f secs\n",
                        diffTime);
  }
}

template <typename FT, int V, int S, bool compress12,
          typename FT_inner, int V_inner, int S_inner, bool compress12_inner>
void pack_nd_clover(QPhiX::Geometry<FT, V, S, compress12> &geom,
                   QPhiX::Geometry<FT_inner, V_inner, S_inner, compress12_inner> &geom_inner,
                   typename QPhiX::Geometry<FT, V, S, compress12>::FullCloverBlock *full_invclov[2],
                   typename QPhiX::Geometry<FT, V, S, compress12>::CloverBlock *invclov_odiag,
                   typename QPhiX::Geometry<FT, V, S, compress12>::CloverBlock *clov,
                   typename QPhiX::Geometry<FT_inner, V_inner, S_inner, compress12_inner>::FullCloverBlock *full_invclov_inner[2],
                   typename QPhiX::Geometry<FT_inner, V_inner, S_inner, compress12_inner>::CloverBlock *invclov_odiag_inner,
                   typename QPhiX::Geometry<FT_inner, V_inner, S_inner, compress12_inner>::CloverBlock *clov_inner,
                   const int cb,
                   bool pack_inner){
  
  typedef typename QPhiX::Geometry<FT, V, S, compress12>::CloverBlock QClover;
  typedef typename QPhiX::Geometry<FT, V, S, compress12>::FullCloverBlock QFullClover;
  typedef typename QPhiX::Geometry<FT_inner, V_inner, S_inner, compress12_inner>::CloverBlock QClover_inner;
  typedef typename QPhiX::Geometry<FT_inner, V_inner, S_inner, compress12_inner>::FullCloverBlock QFullClover_inner;

  double start = gettime(); 
  reorder_clover_to_QPhiX(geom, clov, cb, false);
  if(pack_inner){
    reorder_clover_to_QPhiX(geom_inner, clov_inner, cb, false);
  }
  
  sw_invert_epsbar(g_epsbar);
  reorder_clover_to_QPhiX(geom, invclov_odiag, 1-cb, true, true);
  if(pack_inner){
    reorder_clover_to_QPhiX(geom_inner, invclov_odiag_inner, 1-cb, true, true);
  }

  // no minus sign here, the difference in the sign of gamma5 
  // is taken care of internally
  sw_invert_mubar(g_mubar);
  reorder_clover_to_QPhiX(geom, full_invclov, 1-cb, true);
  if(pack_inner){
    reorder_clover_to_QPhiX(geom_inner, full_invclov_inner, 1-cb, true);
  }

  sw_invert_nd(g_mubar*g_mubar-g_epsbar*g_epsbar);
  
  if(g_debug_level > 1){
    QPhiX::masterPrintf("# QPHIX-inteface: ND TMClover clover-field packing took %.4lf seconds\n", gettime()-start);
  }
}

// Due to github issue #404, the helper functions to apply the full QPhiX operator
// are currently disabled because they conflict with the new interfaces in QPhiX
// itself. If required, these should be rewritten to use these interfaces
// rather than the base classes in qphix_base_classes.hpp

// Apply the full QPhiX fermion matrix to checkerboarded tm spinors
//template <typename FT, int V, int S, bool compress>
//void Mfull_helper(spinor *Even_out, spinor *Odd_out, const spinor *Even_in, const spinor *Odd_in,
//                  const op_type_t op_type) {
//  // TODO: this should use handles for gauge and spinors because these are definitely temporary
//  // objects
//  typedef typename QPhiX::Geometry<FT, V, S, compress>::SU3MatrixBlock QGauge;
//  typedef typename QPhiX::Geometry<FT, V, S, compress>::FourSpinorBlock QSpinor;
//  typedef typename QPhiX::Geometry<FT, V, S, compress>::CloverBlock QClover;
//  typedef typename QPhiX::Geometry<FT, V, S, compress>::FullCloverBlock QFullClover;
//
//  if (g_debug_level > 1) tmlqcd::printQphixDiagnostics(V, S, compress, V, S, compress);
//
//  double coeff_s = (FT)(1);
//  double coeff_t = (FT)(1);
//
//  QPhiX::Geometry<FT, V, S, compress> geom(subLattSize, By, Bz, NCores, Sy, Sz, PadXY, PadXYZ,
//                                           MinCt);
//
//  // Wilson mass
//  double mass = 1 / (2.0 * g_kappa) - 4;
//
//  tmlqcd::Dslash<FT, V, S, compress> *polymorphic_dslash;
//
//  QGauge *u_packed[2];
//  QSpinor *qphix_in[2];
//  QSpinor *qphix_out[2];
//
//  QClover *clover[2];
//  QClover *inv_clover[2];
//
//  QFullClover *inv_fullclover[2][2];
//
//  QSpinor *tmp_spinor = (QSpinor *)geom.allocCBFourSpinor();
//  for (int cb : {0, 1}) {
//    u_packed[cb] = (QGauge *)geom.allocCBGauge();
//    qphix_in[cb] = (QSpinor *)geom.allocCBFourSpinor();
//    qphix_out[cb] = (QSpinor *)geom.allocCBFourSpinor();
//    clover[cb] = nullptr;
//    inv_clover[cb] = nullptr;
//    for (int fl : {0, 1}) {
//      inv_fullclover[cb][fl] = nullptr;
//    }
//  }
//  reorder_gauge_to_QPhiX(geom, u_packed[cb_even], u_packed[cb_odd]);
//
//  if (op_type == WILSON) {
//    polymorphic_dslash = new tmlqcd::WilsonDslash<FT, V, S, compress>(
//        &geom, t_boundary, coeff_s, coeff_t, mass, use_tbc, tbc_phases);
//  } else if (op_type == TMWILSON) {
//    polymorphic_dslash = new tmlqcd::WilsonTMDslash<FT, V, S, compress>(
//        &geom, t_boundary, coeff_s, coeff_t, mass, -g_mu / (2.0 * g_kappa), use_tbc, tbc_phases);
//  } else if (op_type == CLOVER && fabs(g_mu) <= DBL_EPSILON) {
//    for (int cb : {0, 1}) {
//      clover[cb] = (QClover *)geom.allocCBClov();
//      inv_clover[cb] = (QClover *)geom.allocCBClov();
//
//      reorder_clover_to_QPhiX(geom, clover[cb], cb, false);
//      sw_invert(cb, 0);
//      reorder_clover_to_QPhiX(geom, inv_clover[cb], cb, true);
//    }
//
//    polymorphic_dslash = new tmlqcd::WilsonClovDslash<FT, V, S, compress>(
//        &geom, t_boundary, coeff_s, coeff_t, mass, clover, inv_clover, use_tbc, tbc_phases);
//
//  } else if (op_type == CLOVER && fabs(g_mu) > DBL_EPSILON) {
//    for (int cb : {0, 1}) {
//      clover[cb] = (QClover *)geom.allocCBClov();
//      for (int fl : {0, 1}) {
//        inv_fullclover[cb][fl] = (QFullClover *)geom.allocCBFullClov();
//      }
//      reorder_clover_to_QPhiX(geom, clover[cb], cb, false);
//      sw_invert(cb, g_mu);
//      reorder_clover_to_QPhiX(geom, inv_fullclover[cb], cb, true);
//    }
//
//    polymorphic_dslash = new tmlqcd::WilsonClovTMDslash<FT, V, S, compress>(
//        &geom, t_boundary, coeff_s, coeff_t, mass, -g_mu / (2.0 * g_kappa), clover,
//        inv_fullclover, use_tbc, tbc_phases);
//
//  } else {
//    QPhiX::masterPrintf("tmlqcd::Mfull_helper; No such operator type: %d\n", op_type);
//    abort();
//  }
//
////   reorder_eo_spinor_to_QPhiX(geom, reinterpret_cast<double const *const>(Even_in),
////                              qphix_in[cb_even], cb_even);
////   reorder_eo_spinor_to_QPhiX(geom, reinterpret_cast<double const *const>(Odd_in), qphix_in[cb_odd],
////                              cb_odd);
//  reorder_eo_spinor_to_QPhiX(geom, Even_in,
//                             qphix_in[cb_even], cb_even);
//  reorder_eo_spinor_to_QPhiX(geom, Odd_in, qphix_in[cb_odd],
//                             cb_odd);
//  // Apply QPhiX Mfull
//  polymorphic_dslash->plain_dslash(qphix_out[cb_odd], qphix_in[cb_even], u_packed[cb_odd],
//                                   /* isign == non-conjugate */ 1, cb_odd);
//  polymorphic_dslash->plain_dslash(qphix_out[cb_even], qphix_in[cb_odd], u_packed[cb_even],
//                                   /* isign == non-conjugate */ 1, cb_even);
//  for (int cb : {0, 1}) {
//    polymorphic_dslash->A_chi(tmp_spinor, qphix_in[cb], 1, cb);
//    QPhiX::aypx(-0.5, tmp_spinor, qphix_out[cb], geom, 1);
//  }
//
//  reorder_eo_spinor_from_QPhiX(geom, Even_out, qphix_out[cb_even],
//                               cb_even, 2.0 * g_kappa);
//  reorder_eo_spinor_from_QPhiX(geom, Odd_out, qphix_out[cb_odd], cb_odd,
//                               2.0 * g_kappa);
//
//  geom.free(tmp_spinor);
//  for (int cb : {0, 1}) {
//    geom.free(u_packed[cb]);
//    geom.free(qphix_in[cb]);
//    geom.free(qphix_out[cb]);
//    geom.free(clover[cb]);
//    geom.free(inv_clover[cb]);
//    for (int fl : {0, 1}) {
//      geom.free(inv_fullclover[cb][fl]);
//    }
//  };
//  delete (polymorphic_dslash);
//}

// Templated even-odd preconditioned solver using QPhiX Library
template <typename FT, int V, int S, bool compress, 
          typename FT_inner = FT, int V_inner = V, int S_inner = S, bool compress_inner = compress>
int invert_eo_qphix_helper(std::vector< std::vector < spinor* > > &tmlqcd_odd_out, 
                           std::vector< std::vector < spinor* > > &tmlqcd_odd_in,
                           const double target_precision, const int max_iter, const int solver_flag,
                           solver_params_t solver_params, const int num_flavour) {
  // TODO: it would perhaps be beneficial to keep the fields resident
  typedef typename QPhiX::Geometry<FT, V, S, compress>::SU3MatrixBlock QGauge;
  typedef typename QPhiX::Geometry<FT, V, S, compress>::FourSpinorBlock QSpinor;
  typedef typename QPhiX::FourSpinorHandle<FT, V, S, compress> QSpinorHandle;
  typedef typename QPhiX::Geometry<FT, V, S, compress>::CloverBlock QClover;
  typedef typename QPhiX::Geometry<FT, V, S, compress>::FullCloverBlock QFullClover;

  typedef typename QPhiX::Geometry<FT_inner, V_inner, S_inner, compress_inner>::SU3MatrixBlock QGauge_inner;
  typedef typename QPhiX::Geometry<FT_inner, V_inner, S_inner, compress_inner>::FourSpinorBlock QSpinor_inner;
  typedef typename QPhiX::FourSpinorHandle<FT_inner, V_inner, S_inner, compress_inner> QSpinorHandle_inner;
  typedef typename QPhiX::Geometry<FT_inner, V_inner, S_inner, compress_inner>::CloverBlock QClover_inner;
  typedef typename QPhiX::Geometry<FT_inner, V_inner, S_inner, compress_inner>::FullCloverBlock QFullClover_inner;

  /************************
   *                      *
   *    SETUP GEOMETRY    *
   *                      *
  ************************/

  if (g_debug_level > 1) {
    tmlqcd::printQphixDiagnostics(V, S, compress, V_inner, S_inner, compress_inner);
  }

  QPhiX::Geometry<FT, V, S, compress> geom(subLattSize, By, Bz, NCores, Sy, Sz, PadXY, PadXYZ,
                                           MinCt);

  // we always create the inner geometry, the overhead should be small...
  QPhiX::Geometry<FT_inner, V_inner, S_inner, compress_inner> geom_inner(
      subLattSize, By, Bz, NCores, Sy, Sz, PadXY, PadXYZ, MinCt);

  // Set number of BLAS threads by hand.
  // In case some implements the tune routines in QPhiX
  // this may be updated...
  QPhiX::masterPrintf("# Setting number of BLAS threads...\n");
  const int n_blas_simt = N_simt;
  QPhiX::masterPrintf("# ...done.\n");

  // Anisotropy Coefficents
  const double coeff_s = 1.0;
  const double coeff_t = 1.0;

  // The Wilson mass
  const double mass = 1.0 / (2.0 * g_kappa) - 4.0;

  // Set variables need for solve
  bool verbose = g_debug_level > 2 ? true : false;
  int niters = -1;
  int niters2 = 0;
  double rsd_final = -1.0;
  uint64_t site_flops = 0;
  uint64_t site_flops2 = 0;
  uint64_t mv_apps = 0;
  uint64_t mv_apps2 = 0;
  
  double start_time;
  double end_time;
  
  // support for multi-shift solves via the length of the output vector,
  // which counts the shifts on the outer index and the flavour on the inner index
  const int num_shifts = tmlqcd_odd_out.size();
  std::vector < double > shifts; shifts.resize( num_shifts );
  std::vector <double> RsdTargetArr; RsdTargetArr.resize(num_shifts);
  std::vector <double> RsdFinalArr; RsdFinalArr.resize(num_shifts);

  double rescale = 0.5 / g_kappa;
  // the inverse of M M^dag, as required for the HMC, comes with a factor of alpha^2
  if (solver_params.solution_type == TM_SOLUTION_M_MDAG) {
    rescale *= rescale;
  }

  std::vector<QSpinorHandle> q_spinor_handles;

  QGauge *u_packed[2] = {nullptr, nullptr};
  QGauge_inner *u_packed_inner[2] = {nullptr, nullptr};
  for (int cb : {0, 1}) {
    u_packed[cb] = (QGauge *)geom.allocCBGauge();
  }
  // Reorder (global) input gauge field from tmLQCD to QPhiX
  reorder_gauge_to_QPhiX(geom, u_packed[cb_even], u_packed[cb_odd]);
  
  // for mixed solvers, we also need the gauge field in the inner precision
  if( solver_is_mixed(solver_flag) ){
    for(int cb : {0, 1}) {
      u_packed_inner[cb] = (QGauge_inner *)geom_inner.allocCBGauge();
    }
    reorder_gauge_to_QPhiX(geom_inner, u_packed_inner[cb_even], u_packed_inner[cb_odd]);
  }

  if (num_flavour == 1) {
    constexpr int nf = 1;
    std::vector < QSpinor* > qphix_in; qphix_in.resize( 1 );
    std::vector < QSpinor* > qphix_out; qphix_out.resize( num_shifts );
    QSpinor *qphix_buffer;

    QClover *qphix_clover = nullptr;
    QClover *qphix_inv_clover = nullptr;
    
    QClover_inner *qphix_clover_inner = nullptr;
    QClover_inner *qphix_inv_clover_inner = nullptr;

    QFullClover *qphix_inv_fullclover[2] = {nullptr, nullptr};
    
    QFullClover_inner *qphix_inv_fullclover_inner[2] = {nullptr, nullptr};

    q_spinor_handles.push_back(makeFourSpinorHandle(geom));
    qphix_in[0] = q_spinor_handles.back().get();

    for( int shift = 0; shift < num_shifts; shift++ ) {
      q_spinor_handles.push_back(makeFourSpinorHandle(geom));
      qphix_out[shift] = q_spinor_handles.back().get();
    }

    q_spinor_handles.push_back(makeFourSpinorHandle(geom));
    qphix_buffer = q_spinor_handles.back().get();

    QPhiX::EvenOddLinearOperator<FT, V, S, compress> *FermionMatrixQPhiX = nullptr;
    QPhiX::EvenOddLinearOperator<FT_inner, V_inner, S_inner, compress_inner> *InnerFermionMatrixQPhiX = nullptr;
    if ( ( fabs(g_mu) > DBL_EPSILON ) && g_c_sw > DBL_EPSILON) {  // TWISTED-MASS-CLOVER
      qphix_clover = (QClover *)geom.allocCBClov();
      for (int fl : {0, 1}) {
        qphix_inv_fullclover[fl] = (QFullClover *)geom.allocCBFullClov(); 
      }
      reorder_clover_to_QPhiX(geom, qphix_clover, cb_odd, false);
      reorder_clover_to_QPhiX(geom, qphix_inv_fullclover, cb_even, true);
      
      QPhiX::masterPrintf("# Creating QPhiX Twisted Clover Fermion Matrix...\n");
      FermionMatrixQPhiX = new QPhiX::EvenOddTMCloverOperator<FT, V, S, compress>(
          u_packed, qphix_clover, qphix_inv_fullclover, &geom, t_boundary, coeff_s, coeff_t,
          use_tbc, tbc_phases, -0.5*(g_mu3+g_mu)/g_kappa);
      if( solver_is_mixed(solver_flag) ){
        qphix_clover_inner = (QClover_inner *)geom_inner.allocCBClov();
        for( int fl : {0, 1} ){
          qphix_inv_fullclover_inner[fl] = (QFullClover_inner *)geom_inner.allocCBFullClov();
        }
        reorder_clover_to_QPhiX(geom_inner, qphix_clover_inner, cb_odd, false);        
        reorder_clover_to_QPhiX(geom_inner, qphix_inv_fullclover_inner, cb_even, true);
        InnerFermionMatrixQPhiX = new QPhiX::EvenOddTMCloverOperator<FT_inner, V_inner, S_inner, compress_inner>(
          u_packed_inner, qphix_clover_inner, qphix_inv_fullclover_inner, &geom_inner, t_boundary, coeff_s, coeff_t,
          use_tbc, tbc_phases, -0.5*(g_mu3+g_mu)/g_kappa);
      }
      QPhiX::masterPrintf("# ...done.\n");
    } else if ( fabs(g_mu) > DBL_EPSILON ) {  // TWISTED-MASS
      const double TwistedMass = -g_mu / (2.0 * g_kappa);
      QPhiX::masterPrintf("# Creating QPhiX Twisted Mass Wilson Fermion Matrix...\n");
      FermionMatrixQPhiX = new QPhiX::EvenOddTMWilsonOperator<FT, V, S, compress>(
          mass, TwistedMass, u_packed, &geom, t_boundary, coeff_s, coeff_t, use_tbc, tbc_phases);
      QPhiX::masterPrintf("# ...done.\n");
      if( solver_is_mixed(solver_flag) ){
        InnerFermionMatrixQPhiX = new QPhiX::EvenOddTMWilsonOperator<FT_inner, V_inner, S_inner, compress_inner>(
            mass, TwistedMass, u_packed_inner, &geom_inner, t_boundary, coeff_s, coeff_t, use_tbc, tbc_phases);
      }
    } else if (g_c_sw > DBL_EPSILON) {  // WILSON CLOVER
      qphix_clover = (QClover *)geom.allocCBClov();
      qphix_inv_clover = (QClover *)geom.allocCBClov();

      reorder_clover_to_QPhiX(geom, qphix_clover, cb_odd, false);
      reorder_clover_to_QPhiX(geom, qphix_inv_clover, cb_even, true);

      QPhiX::masterPrintf("# Creating QPhiX Wilson Clover Fermion Matrix...\n");
      FermionMatrixQPhiX = new QPhiX::EvenOddCloverOperator<FT, V, S, compress>(
          u_packed, qphix_clover, qphix_inv_clover, &geom, t_boundary, coeff_s, coeff_t, use_tbc,
          tbc_phases, -0.5*g_mu3/g_kappa);
      if( solver_is_mixed(solver_flag) ){
        qphix_clover_inner = (QClover_inner *)geom_inner.allocCBClov();
        qphix_inv_clover_inner = (QClover_inner *)geom_inner.allocCBClov();
        reorder_clover_to_QPhiX(geom_inner, qphix_clover_inner, cb_odd, false);        
        reorder_clover_to_QPhiX(geom_inner, qphix_inv_clover_inner, cb_even, true);
        InnerFermionMatrixQPhiX = new QPhiX::EvenOddCloverOperator<FT_inner, V_inner, S_inner, compress_inner>(
          u_packed_inner, qphix_clover_inner, qphix_inv_clover_inner, &geom_inner, t_boundary, coeff_s, coeff_t,
          use_tbc, tbc_phases, -0.5*g_mu3/g_kappa);
      }
      QPhiX::masterPrintf("# ...done.\n");

    } else {  // WILSON
      QPhiX::masterPrintf("# Creating QPhiX Wilson Fermion Matrix...\n");
      FermionMatrixQPhiX = new QPhiX::EvenOddWilsonOperator<FT, V, S, compress>(
          mass, u_packed, &geom, t_boundary, coeff_s, coeff_t, use_tbc, tbc_phases);
      if( solver_is_mixed(solver_flag) ){
        InnerFermionMatrixQPhiX = new QPhiX::EvenOddWilsonOperator<FT_inner, V_inner, S_inner, compress_inner>(
            mass, u_packed_inner, &geom_inner, t_boundary, coeff_s, coeff_t, use_tbc, tbc_phases);
      }
      QPhiX::masterPrintf("# ...done.\n");
    }

    // Create a Linear Solver Object
    QPhiX::AbstractSolver<FT, V, S, compress> *SolverQPhiX = nullptr;
    QPhiX::AbstractSolver<FT_inner, V_inner, S_inner, compress_inner> *InnerSolverQPhiX = nullptr;
    QPhiX::AbstractMultiSolver<FT, V, S, compress, nf> *MultiSolverQPhiX = nullptr;
    if( solver_flag == DUMMYHERMTEST ) {
      QPhiX::masterPrintf("# QPHIX: Creating dummy solver for hermiticity test...\n");
      SolverQPhiX =
        new QPhiX::InvDummyHermTest<FT, V, S, compress,
                                    typename QPhiX::EvenOddLinearOperator<FT, V, S, compress> >(
              *FermionMatrixQPhiX, max_iter);
    } else if (solver_flag == CG) {
      QPhiX::masterPrintf("# QPHIX: Creating CG solver...\n");
      SolverQPhiX = new QPhiX::InvCG<FT, V, S, compress>(*FermionMatrixQPhiX, max_iter);
    } else if (solver_flag == BICGSTAB) {
      QPhiX::masterPrintf("# QPHIX: Creating BiCGStab solver...\n");
      SolverQPhiX = new QPhiX::InvBiCGStab<FT, V, S, compress>(*FermionMatrixQPhiX, max_iter);
    } else if ( solver_flag == MIXEDCG ) {
      // TODO: probably need to adjust inner solver iterations here...
      QPhiX::masterPrintf("# QPHIX: Creating mixed-precision CG solver...\n");
      InnerSolverQPhiX = new QPhiX::InvCG<FT_inner, V_inner, S_inner, compress_inner>(*InnerFermionMatrixQPhiX, max_iter);
      const bool MMdag = true;
      SolverQPhiX = new QPhiX::InvRichardsonMultiPrec<FT, V, S, compress, FT_inner, V_inner, S_inner, compress_inner, MMdag>(
          *FermionMatrixQPhiX, *InnerSolverQPhiX, solver_params.mcg_delta, max_iter);
    } else if (solver_flag == MIXEDBICGSTAB ) {
      QPhiX::masterPrintf("# QPHIX: Creating mixed-precision BICGCGSTAB solver...\n");
      InnerSolverQPhiX = new QPhiX::InvBiCGStab<FT_inner, V_inner, S_inner, compress_inner>(*InnerFermionMatrixQPhiX, max_iter);
      const bool MMdag = false;
      SolverQPhiX = new QPhiX::InvRichardsonMultiPrec<FT, V, S, compress, FT_inner, V_inner, S_inner, compress_inner, MMdag>(
          *FermionMatrixQPhiX, *InnerSolverQPhiX, solver_params.mcg_delta, max_iter);
    } else if (solver_flag == CGMMS ) {
      QPhiX::masterPrintf("# QPHIX: Creating multi-shift CG solver ...\n");
      MultiSolverQPhiX = new QPhiX::MInvCG<FT, V, S, compress>( *FermionMatrixQPhiX, max_iter, num_shifts );
    } else {
      QPhiX::masterPrintf(" Solver not yet supported by QPhiX!\n");
      QPhiX::masterPrintf(" Aborting...\n");
      abort();
    }
    QPhiX::masterPrintf("# ...done.\n");

//     reorder_eo_spinor_to_QPhiX(geom, reinterpret_cast<double const *const>(tmlqcd_odd_in[0][0]),
//                                qphix_in[0], cb_odd);
    reorder_eo_spinor_to_QPhiX(geom, tmlqcd_odd_in[0][0],
                               qphix_in[0], cb_odd);
    QPhiX::masterPrintf("# Calling the solver...\n");

    // Set the right precision for the QPhiX solver
    // we get target_precision externally and and is given such, that it's either
    // already relative or absolute
    // Most QPhiX solvers allow setting absolute or relative residual
    // by passing an appropriate flag, but this is not true for the multi-shift solver.
    // As a result, we follow that solver and call ALL solvers with
    // QPhiX::RELATIVE, which gives results consistent with tmLQCD in all cases. 
    double rhs_norm2 = 1.0;
    QPhiX::norm2Spinor(rhs_norm2, qphix_in[0], geom, n_blas_simt);
    const double RsdTarget = sqrt(target_precision / rhs_norm2);

    // Calling the solver
    start_time = gettime();
    if ( solver_flag == DUMMYHERMTEST ){
      random_spinor_field_eo(tmlqcd_odd_out[0][0], 0, RN_GAUSS);
      reorder_eo_spinor_to_QPhiX(geom, tmlqcd_odd_out[0][0], qphix_buffer, cb_odd);      
      for(int isign : {-1, 1} ){
        (*SolverQPhiX)(qphix_buffer, qphix_in[0], RsdTarget, niters, rsd_final, site_flops, mv_apps, isign,
                       verbose, cb_odd, QPhiX::RELATIVE);
      }
      QPhiX::copySpinor(qphix_out[0], qphix_buffer, geom, n_blas_simt);
    } else if (solver_flag == CG || solver_flag == MIXEDCG || solver_flag == RGMIXEDCG) {
      // USING CG:
      // We are solving
      //   M M^dagger qphix_buffer = qphix_in_prepared
      // here, that is, isign = -1 for the QPhiX CG solver.
      (*SolverQPhiX)(qphix_buffer, qphix_in[0], RsdTarget, niters, rsd_final, site_flops, mv_apps, -1,
                     verbose, cb_odd, QPhiX::RELATIVE);
      // After that. if required by the solution type, multiply with M^dagger:
      //   qphix_out[1] = M^dagger ( M^dagger^-1 M^-1 ) qphix_in_prepared
      if (solver_params.solution_type == TM_SOLUTION_M) {
        (*FermionMatrixQPhiX)(qphix_out[0], qphix_buffer, /* conjugate */ -1);
        mv_apps++;
      } else {
        QPhiX::copySpinor(qphix_out[0], qphix_buffer, geom, n_blas_simt);
      }
    } else if (solver_flag == CGMMS ){
      // TODO: handle the residuals properly
      if(g_debug_level > 2 ) QPhiX::masterPrintf("# QPHIX CGMMS: shifts: \n");
      for( int shift = 0; shift < num_shifts; shift++ ){
        RsdTargetArr[shift] = RsdTarget;
        RsdFinalArr[shift] = -1.0;
        shifts[shift] = solver_params.shifts[shift]*solver_params.shifts[shift]/(4*g_kappa*g_kappa);
        if(g_debug_level > 2 ) QPhiX::masterPrintf("# QPHIX CGMMS: shift[%d] = %.6e\n", shift, shifts[shift]);
      }
      if(g_debug_level > 2 ) QPhiX::masterPrintf("\n");
      (*MultiSolverQPhiX)(qphix_out.data(), qphix_in[0], num_shifts, shifts.data(), 
                          RsdTargetArr.data(), niters, RsdFinalArr.data(), site_flops, mv_apps, -1, verbose );
      rsd_final = RsdFinalArr[0];
    } else if (solver_flag == BICGSTAB || solver_flag == MIXEDBICGSTAB) {
      (*SolverQPhiX)(qphix_buffer, qphix_in[0], RsdTarget, niters, rsd_final, site_flops, mv_apps, 1,
                     verbose, cb_odd, QPhiX::RELATIVE);
      // for M^dagger^-1 M^-1 solution type, need to call BiCGstab twice
      if (solver_params.solution_type == TM_SOLUTION_M_MDAG) {
        (*SolverQPhiX)(qphix_out[0], qphix_buffer, RsdTarget, niters2, rsd_final, site_flops, mv_apps2,
                       -1, verbose, cb_odd, QPhiX::RELATIVE);
      } else {
        QPhiX::copySpinor(qphix_out[0], qphix_buffer, geom, n_blas_simt);
      }
    }
    end_time = gettime();
    
    for(int shift = 0; shift < num_shifts; shift++ ){
      reorder_eo_spinor_from_QPhiX(geom, tmlqcd_odd_out[shift][0],
                                   qphix_out[shift], cb_odd, rescale);
    }

    QPhiX::masterPrintf("# QPHIX: ...done.\n");
    QPhiX::masterPrintf("# QPHIX: Cleaning up\n");
    delete (FermionMatrixQPhiX);
    delete (InnerFermionMatrixQPhiX);
    delete (SolverQPhiX);
    delete (InnerSolverQPhiX);
    delete (MultiSolverQPhiX);
    // on KNL, it seems that munmap is problematic, so we check for nullptr
    if(qphix_clover) geom.free(qphix_clover);
    if(qphix_inv_clover) geom.free(qphix_inv_clover);
    if(qphix_clover_inner) geom_inner.free(qphix_clover_inner);
    if(qphix_inv_clover_inner) geom_inner.free(qphix_inv_clover_inner);
    for (int fl : {0, 1}) {
      if(qphix_inv_fullclover[fl]) geom.free(qphix_inv_fullclover[fl]);
      if(qphix_inv_fullclover_inner[fl]) geom_inner.free(qphix_inv_fullclover_inner[fl]);
    }
    QPhiX::masterPrintf("# QPHIX: ...done.\n\n");

  } else if (num_flavour == 2) {
    // for explicit template arguments
    constexpr int nf = 2;

    QSpinor *qphix_in[2];
    std::vector < QSpinor** > qphix_out;
    qphix_out.resize( num_shifts );
    for( int shift = 0; shift < num_shifts; shift++ ){
      qphix_out[shift] = new QSpinor*[2];
      for (int fl : {0, 1}) {
        q_spinor_handles.push_back(makeFourSpinorHandle(geom));
        qphix_out[shift][fl] = q_spinor_handles.back().get();
      }
    }
    
    QSpinor *qphix_buffer[2];
    for (int fl : {0, 1}) {
      q_spinor_handles.push_back(makeFourSpinorHandle(geom));
      qphix_in[fl] = q_spinor_handles.back().get();
      q_spinor_handles.push_back(makeFourSpinorHandle(geom));
      qphix_buffer[fl] = q_spinor_handles.back().get();
    }
    
    QClover *qphix_clover = nullptr;
    QClover_inner *qphix_clover_inner = nullptr;
    
    QClover *qphix_invclov_odiag = nullptr;
    QClover_inner *qphix_invclov_odiag_inner = nullptr;

    QFullClover *qphix_inv_fullclover[2] = {nullptr, nullptr};
    QFullClover_inner *qphix_inv_fullclover_inner[2] = {nullptr, nullptr};

    QPhiX::TwoFlavEvenOddLinearOperator<FT, V, S, compress> *TwoFlavFermionMatrixQPhiX = nullptr;
    QPhiX::TwoFlavEvenOddLinearOperator<FT_inner, V_inner, S_inner, compress_inner> *InnerTwoFlavFermionMatrixQPhiX = nullptr;
    
    if (g_c_sw > DBL_EPSILON) {  // DBCLOVER
      qphix_clover = (QClover *)geom.allocCBClov();
      qphix_invclov_odiag = (QClover *)geom.allocCBClov();
      if( solver_is_mixed(solver_flag) ){
        qphix_clover_inner = (QClover_inner *)geom_inner.allocCBClov();
        qphix_invclov_odiag_inner = (QClover_inner *)geom_inner.allocCBClov();
      }

      for (int fl : {0, 1}) {
        qphix_inv_fullclover[fl] = (QFullClover *)geom.allocCBFullClov();
        if( solver_is_mixed(solver_flag) ){
          qphix_inv_fullclover_inner[fl] = (QFullClover_inner*)geom_inner.allocCBFullClov();
        }
      }

      pack_nd_clover(geom, geom_inner, 
                     qphix_inv_fullclover, qphix_invclov_odiag, qphix_clover,
                     qphix_inv_fullclover_inner, qphix_invclov_odiag_inner, qphix_clover_inner,
                     cb_odd,
                     solver_is_mixed(solver_flag));

      QPhiX::masterPrintf("# QPHIX: Creating two-flavour QPhiX Wilson Twisted Clover Fermion Matrix...\n");
      TwoFlavFermionMatrixQPhiX = new QPhiX::EvenOddNDTMCloverReuseOperator<FT, V, S, compress>(
          -0.5 * g_mubar / g_kappa, 0.5 * g_epsbar / g_kappa, 
          u_packed, qphix_clover, qphix_invclov_odiag, qphix_inv_fullclover,
          &geom, t_boundary,
          coeff_s, coeff_t, use_tbc, tbc_phases);
      if( solver_is_mixed(solver_flag) ){
        InnerTwoFlavFermionMatrixQPhiX = new QPhiX::EvenOddNDTMCloverReuseOperator<FT_inner, V_inner, S_inner, compress_inner>(
            -0.5 * g_mubar / g_kappa, 0.5 * g_epsbar / g_kappa, 
            u_packed_inner, qphix_clover_inner, qphix_invclov_odiag_inner, qphix_inv_fullclover_inner,
            &geom_inner, t_boundary,
            coeff_s, coeff_t, use_tbc, tbc_phases);
      }
    } else {  // DBTMWILSON
      QPhiX::masterPrintf("# QPHIX: Creating two-flavour QPhiX Wilson Twisted Mass Fermion Matrix...\n");
      TwoFlavFermionMatrixQPhiX = new QPhiX::EvenOddNDTMWilsonReuseOperator<FT, V, S, compress>(
          mass, -0.5 * g_mubar / g_kappa, 0.5 * g_epsbar / g_kappa, u_packed, &geom, t_boundary,
          coeff_s, coeff_t, use_tbc, tbc_phases);
      if( solver_is_mixed(solver_flag) ){
        InnerTwoFlavFermionMatrixQPhiX = new QPhiX::EvenOddNDTMWilsonReuseOperator<
                                               FT_inner, V_inner, S_inner, compress_inner>(
          mass, -0.5 * g_mubar / g_kappa, 0.5 * g_epsbar / g_kappa, u_packed_inner, &geom_inner, t_boundary,
          coeff_s, coeff_t, use_tbc, tbc_phases);
      }
    }

    //
    QPhiX::AbstractSolver<FT, V, S, compress, nf> *TwoFlavSolverQPhiX = nullptr;
    QPhiX::AbstractSolver<FT_inner, V_inner, S_inner, compress_inner, nf> *InnerTwoFlavSolverQPhiX = nullptr;
    QPhiX::AbstractMultiSolver<FT, V, S, compress, nf> *TwoFlavMultiSolverQPhiX = nullptr;
    if( solver_flag == DUMMYHERMTEST ) {
      QPhiX::masterPrintf("# QPHIX: Creating dummy solver for hermiticity test...\n");
      TwoFlavSolverQPhiX =
        new QPhiX::InvDummyHermTest<FT, V, S, compress,
                                    typename QPhiX::TwoFlavEvenOddLinearOperator<FT, V, S, compress> >(
              *TwoFlavFermionMatrixQPhiX, max_iter);
    } else if (solver_flag == CG) {
      QPhiX::masterPrintf("# QPHIX: Creating CG solver...\n");
      TwoFlavSolverQPhiX =
          new QPhiX::InvCG<FT, V, S, compress,
                           typename QPhiX::TwoFlavEvenOddLinearOperator<FT, V, S, compress> >(
              *TwoFlavFermionMatrixQPhiX, max_iter);
    } else if (solver_flag == BICGSTAB) {
      QPhiX::masterPrintf("# QPHIX: Creating BiCGstab solver...\n");
      TwoFlavSolverQPhiX =
          new QPhiX::InvBiCGStab<FT, V, S, compress,
                                 typename QPhiX::TwoFlavEvenOddLinearOperator<FT, V, S, compress> >(
              *TwoFlavFermionMatrixQPhiX, max_iter);
    } else if (solver_flag == MIXEDCG) {
      QPhiX::masterPrintf("# QPHIX: Creating mixed-precision CG solver...\n");
      InnerTwoFlavSolverQPhiX = new QPhiX::InvCG<
          FT_inner, V_inner, S_inner, compress_inner,
          typename QPhiX::TwoFlavEvenOddLinearOperator<FT_inner, V_inner, S_inner, compress_inner> 
        >(
            *InnerTwoFlavFermionMatrixQPhiX, max_iter
         );
      const bool MMdag = true;
      TwoFlavSolverQPhiX = new QPhiX::InvRichardsonMultiPrec
        <
          FT, V, S, compress, 
          FT_inner, V_inner, S_inner, compress_inner,
          MMdag, typename QPhiX::TwoFlavEvenOddLinearOperator<FT, V, S, compress> 
        >(
            *TwoFlavFermionMatrixQPhiX,
            *InnerTwoFlavSolverQPhiX,
            solver_params.mcg_delta,
            max_iter
        );
    } else if ( solver_flag == CGMMSND ) {
      QPhiX::masterPrintf("# QPHIX: Creating multi-shift CG solver...\n");
      TwoFlavMultiSolverQPhiX = new QPhiX::MInvCG<FT, V, S, compress, 
                                                  typename QPhiX::TwoFlavEvenOddLinearOperator<FT, V, S, compress>
                                                  >( *TwoFlavFermionMatrixQPhiX, max_iter, num_shifts );
    } else {
      QPhiX::masterPrintf(" Solver not yet supported by QPhiX!\n");
      QPhiX::masterPrintf(" Aborting...\n");
      abort();
    }
    QPhiX::masterPrintf("# QPHIX: ...done.\n");

    for (int fl : {0, 1}) {
//       reorder_eo_spinor_to_QPhiX(geom, reinterpret_cast<double const *const>(tmlqcd_odd_in[0][fl]),
//                                  qphix_in[fl], cb_odd);
      reorder_eo_spinor_to_QPhiX(geom, tmlqcd_odd_in[0][fl],
                                 qphix_in[fl], cb_odd);      
    }

    QPhiX::masterPrintf("# QPHIX: Calling the solver...\n");

    // Set the right precision for the QPhiX solver
    // we get target_precision externally and and is given such, that it's either
    // already relative or absolute
    // Most QPhiX solvers allow setting absolute or relative residual
    // by passing an appropriate flag, but this is not true for the multi-shift solver.
    // As a result, we follow that solver and call ALL solvers with
    // QPhiX::RELATIVE, which gives results consistent with tmLQCD in all cases. 
    double rhs_norm2 = 1.0;
    QPhiX::norm2Spinor<FT, V, S, compress, nf>(rhs_norm2, qphix_in, geom, n_blas_simt);
    const double RsdTarget = sqrt(target_precision / rhs_norm2);

    // Calling the solver
    start_time = gettime();
    if ( solver_flag == DUMMYHERMTEST ){
      for(int fl : {0, 1}){
        random_spinor_field_eo(tmlqcd_odd_out[0][fl], 0, RN_GAUSS);
        reorder_eo_spinor_to_QPhiX(geom, tmlqcd_odd_out[0][fl], qphix_buffer[fl], cb_odd);      
      }
      for( int isign : {-1, 1} ){
        (*TwoFlavSolverQPhiX)(qphix_buffer, qphix_in, RsdTarget, niters, rsd_final, site_flops, mv_apps, isign,
                       verbose, cb_odd, QPhiX::RELATIVE);
      }
      QPhiX::copySpinor<FT, V, S, compress, nf>(qphix_out[0], qphix_buffer, geom, n_blas_simt);
    } else if (solver_flag == CG || solver_flag == MIXEDCG) {
      // USING CG:
      // We are solving
      //   M M^dagger qphix_buffer = qphix_in_prepared
      // here, that is, isign = -1 for the QPhiX CG solver.
      (*TwoFlavSolverQPhiX)(qphix_buffer, qphix_in, RsdTarget, niters, rsd_final, site_flops, mv_apps, -1,
                     verbose, cb_odd, QPhiX::RELATIVE);
      // After that. if required by the solution type, multiply with M^dagger:
      //   qphix_out[1] = M^dagger M^dagger^-1 M^-1 qphix_in_prepared
      if (solver_params.solution_type == TM_SOLUTION_M) {
        (*TwoFlavFermionMatrixQPhiX)(qphix_out[0], qphix_buffer, /* conjugate */ -1);
        mv_apps++;
      } else {
        QPhiX::copySpinor<FT, V, S, compress, nf>(qphix_out[0], qphix_buffer, geom, n_blas_simt);
      }
    } else if (solver_flag == BICGSTAB || solver_flag == MIXEDBICGSTAB) {
      (*TwoFlavSolverQPhiX)(qphix_buffer, qphix_in, RsdTarget, niters, rsd_final, site_flops, mv_apps, 1,
                            verbose, cb_odd, QPhiX::RELATIVE);
      // for M^dagger^-1 M^-1 solution type, need to call BiCGstab twice
      if (solver_params.solution_type == TM_SOLUTION_M_MDAG) {
        (*TwoFlavSolverQPhiX)(qphix_out[0], qphix_buffer, RsdTarget, niters2, rsd_final, site_flops, mv_apps2,
                              -1, verbose, cb_odd, QPhiX::RELATIVE);
      } else {
        QPhiX::copySpinor<FT, V, S, compress, nf>(qphix_out[0], qphix_buffer, geom, n_blas_simt);
      }
    } else if (solver_flag == CGMMSND ){
      // TODO: handle the residuals properly
      if(g_debug_level > 2 ) QPhiX::masterPrintf("# QPHIX CGMMSND: shifts: \n");
      // tmLQCD weights the operator with 1/maxev in the RHMC relative to the shifts
      // we will do this externally on the inverse (in monomial_solve) and thus need to weight
      // the shifts by maxev^2
      const double maxev_sq = (1.0/phmc_invmaxev)*(1.0/phmc_invmaxev);
      for( int shift = 0; shift < num_shifts; shift++ ){
        RsdTargetArr[shift] = RsdTarget;
        RsdFinalArr[shift] = -1.0;
        shifts[shift] = maxev_sq *
                        solver_params.shifts[shift]*solver_params.shifts[shift] /
                        (4*g_kappa*g_kappa); 
        if(g_debug_level > 2 ) QPhiX::masterPrintf("# [%d] = %lf\n", shift, shifts[shift]);
      }
      if(g_debug_level > 2 ) QPhiX::masterPrintf("\n");
      (*TwoFlavMultiSolverQPhiX)(qphix_out.data(), qphix_in, num_shifts, shifts.data(), 
                          RsdTargetArr.data(), niters, RsdFinalArr.data(), site_flops, mv_apps, -1, verbose );
      rsd_final = RsdFinalArr[0];
    }
    end_time = gettime();

    for(int shift = 0; shift < num_shifts; shift++){
      for (int fl : {0, 1}) {
        reorder_eo_spinor_from_QPhiX(geom, tmlqcd_odd_out[shift][fl],
                                    qphix_out[shift][fl], cb_odd, rescale);
      }
    }

    delete TwoFlavFermionMatrixQPhiX;
    delete InnerTwoFlavFermionMatrixQPhiX;
    delete InnerTwoFlavSolverQPhiX;
    delete TwoFlavMultiSolverQPhiX;
    delete TwoFlavSolverQPhiX;
    for( int shift = 0; shift < num_shifts; shift++ ){
      delete[] qphix_out[shift];
    }
    
    if(qphix_clover) geom.free(qphix_clover);
    if(qphix_invclov_odiag) geom.free(qphix_invclov_odiag);
    if(qphix_clover_inner) geom_inner.free(qphix_clover_inner);
    if(qphix_invclov_odiag_inner) geom_inner.free(qphix_invclov_odiag_inner);
    for (int fl : {0, 1}) {
      if(qphix_inv_fullclover[fl]) geom.free(qphix_inv_fullclover[fl]);
      if(qphix_inv_fullclover_inner[fl]) geom_inner.free(qphix_inv_fullclover_inner[fl]);
    }    

  } else { // if(num_flavour)
    // complain, this number of flavours is not valid
  } // if(num_flavour)

  for (int cb : {0, 1}) {
    if(u_packed[cb]) geom.free(u_packed[cb]);
    if(u_packed_inner[cb]) geom_inner.free(u_packed_inner[cb]);
  }

  // FIXME: This should be called properly somewhere else
  _endQphix();

  QPhiX::masterPrintf("# ...done.\n\n");

  uint64_t num_cb_sites = lattSize[0] / 2 * lattSize[1] * lattSize[2] * lattSize[3];
  // FIXME: this needs to be adjusted depending on the operator used
  uint64_t op_flops_per_site = 1320;
  uint64_t total_flops = (site_flops + site_flops2 + 
                          (2 * num_flavour * op_flops_per_site) * (mv_apps+mv_apps2) ) * num_cb_sites;
  QPhiX::masterPrintf("# QPHIX: Solver Time = %g sec\n", (end_time - start_time));
  QPhiX::masterPrintf("# QPHIX: Performance in GFLOPS = %g\n\n", 1.0e-9 * total_flops / (end_time - start_time));  

  if( solver_is_mixed(solver_flag) ){
    // the mixed solver reports the outer iterations, we would like to get
    // some better total
    niters = mv_apps / 2;
    if( solver_flag == MIXEDBICGSTAB && solver_params.solution_type == TM_SOLUTION_M_MDAG ){
      niters2 = mv_apps2 / 2;
    }
  }
  // solver did not converge in maximum number of iterations
  // FIXME: non-convergence does not work correctly yet
  if( (niters+niters2) > max_iter ){
    niters = -1;
    niters2 = 0;
  }
  return (niters+niters2);
}

// Due to github issue #404, the helper functions to apply the full QPhiX operator
// are currently disabled because they conflict with the new interfaces in QPhiX
// itself. If required, these should be rewritten to use these interfaces
// rather than the base classes in qphix_base_classes.hpp

// Template wrapper for the Dslash operator call-able from C code
//void Mfull_qphix(spinor *Even_out, spinor *Odd_out, const spinor *Even_in, const spinor *Odd_in,
//                 const op_type_t op_type) {
//  tmlqcd::checkQphixInputParameters(qphix_input);
//  // FIXME: two-row gauge compression and double precision hard-coded
//  _initQphix(0, nullptr, qphix_input, 12, QPHIX_DOUBLE_PREC);
//
//  if (qphix_precision == QPHIX_DOUBLE_PREC) {
//    if (QPHIX_SOALEN > VECLEN_DP) {
//      QPhiX::masterPrintf("SOALEN=%d is greater than the double prec VECLEN=%d\n", QPHIX_SOALEN,
//                          VECLEN_DP);
//      abort();
//    }
//    QPhiX::masterPrintf("TESTING IN DOUBLE PRECISION \n");
//    if (compress12) {
//      Mfull_helper<double, VECLEN_DP, QPHIX_SOALEN, true>(Even_out, Odd_out, Even_in, Odd_in,
//                                                          op_type);
//    } else {
//      Mfull_helper<double, VECLEN_DP, QPHIX_SOALEN, false>(Even_out, Odd_out, Even_in, Odd_in,
//                                                           op_type);
//    }
//  } else if (qphix_precision == QPHIX_FLOAT_PREC) {
//    if (QPHIX_SOALEN > VECLEN_SP) {
//      QPhiX::masterPrintf("SOALEN=%d is greater than the single prec VECLEN=%d\n", QPHIX_SOALEN,
//                          VECLEN_SP);
//      abort();
//    }
//    QPhiX::masterPrintf("TESTING IN SINGLE PRECISION \n");
//    if (compress12) {
//      Mfull_helper<float, VECLEN_SP, QPHIX_SOALEN, true>(Even_out, Odd_out, Even_in, Odd_in,
//                                                         op_type);
//    } else {
//      Mfull_helper<float, VECLEN_SP, QPHIX_SOALEN, false>(Even_out, Odd_out, Even_in, Odd_in,
//                                                          op_type);
//    }
//  }
//#if (defined(QPHIX_MIC_SOURCE) || defined(QPHIX_AVX512_SOURCE))
//  else if (qphix_precision == QPHIX_HALF_PREC) {
//    if (QPHIX_SOALEN > VECLEN_HP) {
//      QPhiX::masterPrintf("SOALEN=%d is greater than the half prec VECLEN=%d\n", QPHIX_SOALEN,
//                          VECLEN_HP);
//      abort();
//    }
//    QPhiX::masterPrintf("TESTING IN HALF PRECISION \n");
//    if (compress12) {
//      Mfull_helper<QPhiX::half, VECLEN_HP, QPHIX_SOALEN, true>(Even_out, Odd_out, Even_in, Odd_in,
//                                                               op_type);
//    } else {
//      Mfull_helper<QPhiX::half, VECLEN_HP, QPHIX_SOALEN, false>(Even_out, Odd_out, Even_in, Odd_in,
//                                                                op_type);
//    }
//  }
//#endif
//}

// we have a unified interface for n-flavour inversions, but we need to provide wrappers
// which can be called by the tmLQCD solver drivers for one and two-flavour inversions
int invert_eo_qphix_oneflavour(spinor *Odd_out_1f, spinor *Odd_in_1f, const int max_iter,
                               const double precision, const int solver_flag, const int rel_prec,
                               const solver_params_t solver_params, const SloppyPrecision sloppy,
                               const CompressionType compression) {
  const int num_flavour = 1;
  const int num_shifts = 1;
  std::vector< std::vector < spinor* > > Odd_out;
  std::vector< std::vector < spinor* > > Odd_in;
  
  Odd_out.resize( num_shifts ); Odd_out[0].resize( num_flavour );
  Odd_in.resize( 1 ); Odd_in[0].resize( num_flavour );
  
  Odd_in[0][0] = Odd_in_1f;
  Odd_out[0][0] = Odd_out_1f;
  
  return invert_eo_qphix_nflavour_mshift(Odd_out, Odd_in, precision, max_iter,
                                         solver_flag, rel_prec,
                                         solver_params, sloppy, compression, num_flavour);
}

int invert_eo_qphix_oneflavour_mshift(spinor **Odd_out_1f, spinor *Odd_in_1f, const int max_iter,
                                      const double precision, const int solver_flag, const int rel_prec,
                                      const solver_params_t solver_params, const SloppyPrecision sloppy,
                                      const CompressionType compression) {
  // even though the default is set to 1, guard against zeroes
  const int num_shifts = solver_params.no_shifts == 0 ? 1 : solver_params.no_shifts;
  const int num_flavour = 1;
  std::vector< std::vector < spinor* > > Odd_out;
  std::vector< std::vector < spinor* > > Odd_in;
  
  Odd_out.resize( num_shifts );
  Odd_in.resize( 1 ); Odd_in[0].resize( num_flavour );
  
  Odd_in[0][0] = Odd_in_1f;
  for( int shift = 0; shift < num_shifts; shift++ ){
    Odd_out[shift].resize( num_flavour );
    Odd_out[shift][0] = Odd_out_1f[shift];
  }

  return invert_eo_qphix_nflavour_mshift(Odd_out, Odd_in, precision, max_iter,
                                        solver_flag, rel_prec,
                                        solver_params, sloppy, compression, num_flavour);
}

// Template wrapper for QPhiX solvers callable from C code, return number of iterations
int invert_eo_qphix_twoflavour(spinor *Odd_out_s, spinor *Odd_out_c, spinor *Odd_in_s,
                               spinor *Odd_in_c, const int max_iter, const double precision,
                               const int solver_flag, const int rel_prec,
                               const solver_params_t solver_params, const SloppyPrecision sloppy,
                               const CompressionType compression) {
  const int num_flavour = 2;
  const int num_shifts = 1;
  std::vector< std::vector < spinor* > > Odd_out;
  std::vector< std::vector < spinor* > > Odd_in;
  
  Odd_out.resize( num_shifts ); Odd_out[0].resize( num_flavour );
  Odd_in.resize( 1 ); Odd_in[0].resize( num_flavour );
  
  Odd_in[0][0] = Odd_in_s; 
  Odd_in[0][1] = Odd_in_c;
  
  Odd_out[0][0] = Odd_out_s;
  Odd_out[0][1] = Odd_out_c;
  
  return invert_eo_qphix_nflavour_mshift(Odd_out, Odd_in, precision, max_iter,
                                         solver_flag, rel_prec,
                                         solver_params, sloppy, compression, num_flavour);
}

int invert_eo_qphix_twoflavour_mshift(spinor **Odd_out_s, spinor **Odd_out_c, spinor *Odd_in_s,
                                      spinor *Odd_in_c, const int max_iter, const double precision,
                                      const int solver_flag, const int rel_prec,
                                      const solver_params_t solver_params, const SloppyPrecision sloppy,
                                      const CompressionType compression) {
  // even though the default is set to 1, guard against zeroes
  const int num_shifts = solver_params.no_shifts == 0 ? 1 : solver_params.no_shifts;  
  const int num_flavour = 2;
  std::vector< std::vector < spinor* > > Odd_out;
  std::vector< std::vector < spinor* > > Odd_in;
  
  Odd_out.resize( num_shifts );
  Odd_in.resize( 1 ); Odd_in[0].resize( num_flavour );
  
  Odd_in[0][0] = Odd_in_s; 
  Odd_in[0][1] = Odd_in_c;
  
  for( int shift = 0; shift < num_shifts; shift++ ){
    Odd_out[shift].resize( num_flavour );
    Odd_out[shift][0] = Odd_out_s[shift];
    Odd_out[shift][1] = Odd_out_c[shift];
  }
  
  return invert_eo_qphix_nflavour_mshift(Odd_out, Odd_in, precision, max_iter,
                                         solver_flag, rel_prec,
                                         solver_params, sloppy, compression, num_flavour);
}

// Template wrapper for QPhiX solvers callable from C code, return number of iterations
// the interface is prepared for multi-rhs solves, hence the double vector for the input
int invert_eo_qphix_nflavour_mshift(std::vector< std::vector< spinor* > > &Odd_out, 
                                    std::vector< std::vector< spinor* > > &Odd_in, 
                                    const double precision,
                                    const int max_iter,
                                    const int solver_flag, 
                                    const int rel_prec,
                                    solver_params_t solver_params,
                                    const SloppyPrecision sloppy, const CompressionType compression,
                                    const int num_flavour) {
  tmlqcd::checkQphixInputParameters(qphix_input);
  double target_precision = precision;
  double src_norm = 0.0;
  for (int f = 0; f < num_flavour; ++f) {
    src_norm += square_norm(Odd_in[0][f], VOLUME / 2, 1);
  }
  // we use "precision_lambda" to determine if a system can be solved in half or float
  // precision (when a fixed-precision solver is used)
  double precision_lambda = target_precision / src_norm;
  if (rel_prec == 1) {
    QPhiX::masterPrintf("# QPHIX: Using relative precision\n");
    target_precision = precision * src_norm;
    precision_lambda = precision;
  }
  QPhiX::masterPrintf("# QPHIX: precision_lambda: %g, target_precision: %g\n\n", precision_lambda,
                      target_precision);

  // mixed solvers require inner and outer precisions, which we specify explicitly here
  if( solver_is_mixed(solver_flag) ){
#if (defined(QPHIX_MIC_SOURCE) || defined(QPHIX_AVX512_SOURCE))
    if( sloppy == SLOPPY_HALF ){
      if( QPHIX_SOALEN > VECLEN_DP || QPHIX_SOALEN > VECLEN_HP ){
        QPhiX::masterPrintf("SOALEN=%d is greater than the half prec VECLEN=%d or the double prec VECLEN=%d\n", QPHIX_SOALEN,
                            VECLEN_HP, VECLEN_DP);
        abort();
      }
      QPhiX::masterPrintf("# INITIALIZING QPHIX MIXED SOLVER\n");
      QPhiX::masterPrintf("# USING DOUBLE-HALF PRECISION\n");
      _initQphix(0, nullptr, qphix_input, compression, QPHIX_DOUBLE_PREC, QPHIX_HALF_PREC);
      if(compress12) {
        return invert_eo_qphix_helper<
          double, VECLEN_DP, QPHIX_SOALEN, true, QPhiX::half, VECLEN_HP, QPHIX_SOALEN, true>(
            Odd_out, Odd_in, target_precision, max_iter, solver_flag, solver_params, num_flavour);
      } else {
        return invert_eo_qphix_helper<
          double, VECLEN_DP, QPHIX_SOALEN, false, QPhiX::half, VECLEN_HP, QPHIX_SOALEN, false>(
            Odd_out, Odd_in, target_precision, max_iter, solver_flag, solver_params, num_flavour);
      }
    } else
#else
    if (sloppy == SLOPPY_HALF) {
      QPhiX::masterPrintf("QPHIX interface: half precision not supported on this architecture!\n");
      abort();
    } else
#endif
    if( sloppy == SLOPPY_SINGLE ) {
      if( QPHIX_SOALEN > VECLEN_DP || QPHIX_SOALEN > VECLEN_SP ){
        QPhiX::masterPrintf("SOALEN=%d is greater than the single prec VECLEN=%d or the double prec VECLEN=%d\n", QPHIX_SOALEN,
                            VECLEN_SP, VECLEN_DP);
        abort();
      }
      QPhiX::masterPrintf("# INITIALIZING QPHIX MIXED SOLVER\n");
      QPhiX::masterPrintf("# USING DOUBLE-SINGLE PRECISION\n");
      _initQphix(0, nullptr, qphix_input, compression, QPHIX_DOUBLE_PREC, QPHIX_FLOAT_PREC);
      if(compress12) {
        return invert_eo_qphix_helper<
          double, VECLEN_DP, QPHIX_SOALEN, true, float, VECLEN_SP, QPHIX_SOALEN, true>(
            Odd_out, Odd_in, target_precision, max_iter, solver_flag, solver_params, num_flavour);
      } else {
        return invert_eo_qphix_helper<
          double, VECLEN_DP, QPHIX_SOALEN, false, float, VECLEN_SP, QPHIX_SOALEN, false>(
            Odd_out, Odd_in, target_precision, max_iter, solver_flag, solver_params, num_flavour);
      }
    } else { // if(sloppy)
      if( QPHIX_SOALEN > VECLEN_DP ){
        QPhiX::masterPrintf("SOALEN=%d is greater than the double prec VECLEN=%d\n", QPHIX_SOALEN,
                            VECLEN_DP);
        abort();
      }
      QPhiX::masterPrintf("# INITIALIZING QPHIX MIXED SOLVER\n");
      QPhiX::masterPrintf("# USING DOUBLE-DOUBLE PRECISION\n");
      _initQphix(0, nullptr, qphix_input, compression, QPHIX_DOUBLE_PREC, QPHIX_DOUBLE_PREC);
      if(compress12) {
        return invert_eo_qphix_helper<
          double, VECLEN_DP, QPHIX_SOALEN, true>(
            Odd_out, Odd_in, target_precision, max_iter, solver_flag, solver_params, num_flavour);
      } else {
        return invert_eo_qphix_helper<
          double, VECLEN_DP, QPHIX_SOALEN, false>(
            Odd_out, Odd_in, target_precision, max_iter, solver_flag, solver_params, num_flavour);
      }
    } // if( sloppy )
  } else { // if( solver_is_mixed )
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
            Odd_out, Odd_in, target_precision, max_iter, solver_flag, solver_params, num_flavour);
      } else {
        return invert_eo_qphix_helper<QPhiX::half, VECLEN_HP, QPHIX_SOALEN, false>(
            Odd_out, Odd_in, target_precision, max_iter, solver_flag, solver_params, num_flavour);
      }
    } else
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
            Odd_out, Odd_in, target_precision, max_iter, solver_flag, solver_params, num_flavour);
      } else {
        return invert_eo_qphix_helper<float, VECLEN_SP, QPHIX_SOALEN, false>(
            Odd_out, Odd_in, target_precision, max_iter, solver_flag, solver_params, num_flavour);
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
            Odd_out, Odd_in, target_precision, max_iter, solver_flag, solver_params, num_flavour);
      } else {
        return invert_eo_qphix_helper<double, VECLEN_DP, QPHIX_SOALEN, false>(
            Odd_out, Odd_in, target_precision, max_iter, solver_flag, solver_params, num_flavour);
      }
    }  // if( sloppy || target_precision )
  } // if ( solver_flag == *MIXEDCG )
  return -1;
}

void tmlqcd::checkQphixInputParameters(const tm_QPhiXParams_t &params) {
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

void tmlqcd::printQphixDiagnostics(int VECLEN, int SOALEN, bool compress, int VECLEN_inner, int SOALEN_inner, bool compress_inner) {
  QPhiX::masterPrintf("# QphiX: VECLEN=%d SOALEN=%d VECLEN_inner=%d, SOALEN_inner=%d\n", VECLEN, SOALEN, VECLEN_inner, SOALEN_inner);

  QPhiX::masterPrintf("# QphiX: Declared QMP Topology (xyzt):");
  for (int mu = 0; mu < 4; mu++) QPhiX::masterPrintf(" %d", qmp_geom[mu]);
  QPhiX::masterPrintf("\n");

  QPhiX::masterPrintf("# QphiX: Mapping of dimensions QMP -> tmLQCD (xyzt):");
  for (int mu = 0; mu < 4; mu++) QPhiX::masterPrintf(" %d->%d", mu, qmp_tm_map[mu]);
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
  if (compress_inner) {
    QPhiX::masterPrintf("# QphiX: Inner solver using two-row gauge compression (compress12)\n");
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

//   reorder_eo_spinor_to_QPhiX(geom, reinterpret_cast<double const *const>(Even_in),
//                              qphix_cb_even.get(), cb_even);
//   reorder_eo_spinor_to_QPhiX(geom, reinterpret_cast<double const *const>(Odd_in),
//                              qphix_cb_odd.get(), cb_odd);
  reorder_eo_spinor_to_QPhiX(geom, Even_in,
                             qphix_cb_even.get(), cb_even);
  reorder_eo_spinor_to_QPhiX(geom, Odd_in,
                             qphix_cb_odd.get(), cb_odd);
  
  reorder_eo_spinor_from_QPhiX(geom, Even_out, qphix_cb_even.get(),
                               cb_even, 1.0);
  reorder_eo_spinor_from_QPhiX(geom, Odd_out, qphix_cb_odd.get(),
                               cb_odd, 1.0);

  diff(tmp[0], Even_out, Even_in, VOLUME / 2);
  diff(tmp[1], Odd_out, Odd_in, VOLUME / 2);
  double l2norm = square_norm(tmp[0], VOLUME / 2, 1) + square_norm(tmp[1], VOLUME / 2, 1);
  QPhiX::masterPrintf("QPHIX eo spinor packer back and forth difference L2 norm: %lf\n", l2norm);
  finalize_solver(tmp, 2);
}

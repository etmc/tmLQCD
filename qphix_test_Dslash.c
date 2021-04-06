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

#ifdef HAVE_CONFIG_H
#include <tmlqcd_config.h>
#endif
#ifdef TM_USE_QPHIX
#include <qphix/qphix_config.h>
#endif
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#if (defined BGL && !defined BGP)
#include <rts.h>
#endif
#ifdef TM_USE_MPI
#include <mpi.h>
#ifdef HAVE_LIBLEMON
#include <io/gauge.h>
#include <io/params.h>
#endif
#endif
#ifdef TM_USE_OMP
#include <omp.h>
#include "init/init_openmp.h"
#endif
#ifdef QPHIX_QMP_COMMS
#include <qmp.h>
#endif
#include "boundary.h"
#include "gamma.h"
#include "geometry_eo.h"
#include "gettime.h"
#include "global.h"
#include "init/init.h"
#include "init/init.h"
#include "invert_clover_eo.h"
#include "invert_eo.h"
#include "linalg/assign_add_mul_r.h"
#include "linalg/convert_eo_to_lexic.h"
#include "linalg/diff_and_square_norm.h"
#include "linalg/square_norm.h"
#include "mpi_init.h"
#include "operator.h"
#include "operator/D_psi.h"
#include "operator/Hopping_Matrix.h"
#include "operator/Hopping_Matrix_nocom.h"
#include "operator/clover_leaf.h"
#include "operator/clovertm_operators.h"
#include "operator/clovertm_operators.h"
#include "operator/tm_operators.h"
#include "prepare_source.h"
#include "qphix_interface.h"
#include "ranlxd.h"
#include "read_input.h"
#include "solver/cg_her.h"
#include "solver/solver_field.h"
#include "start.h"
#include "su3.h"
#include "su3adj.h"
#include "test/check_geometry.h"
#include "update_backward_gauge.h"
#include "xchange/xchange.h"
#include "struct_accessors.h"

int check_xchange();
double compare_spinors(spinor* s1, spinor* s2);

int main(int argc, char* argv[]) {
  int j;
#ifdef HAVE_LIBLEMON
  paramsXlfInfo* xlfInfo;
#endif
  int status = 0;

  static double tm_t1, tm_t2, q_t1, q_t2;

  DUM_DERI = 8;
  DUM_MATRIX = DUM_DERI + 5;
  NO_OF_SPINORFIELDS = DUM_MATRIX + 4;

  /* Set the input file */
  char input_filename[500];
  snprintf(input_filename, 500, "test_Dslash.input");

  init_parallel_and_read_input(argc, argv, input_filename);
  tmlqcd_mpi_init(argc, argv);
  g_dbw2rand = 0;

#ifdef _GAUGE_COPY
  init_gauge_field(VOLUMEPLUSRAND, 1);
#else
  init_gauge_field(VOLUMEPLUSRAND, 0);
#endif

  init_geometry_indices(VOLUMEPLUSRAND);
  j = init_spinor_field(VOLUMEPLUSRAND, NO_OF_SPINORFIELDS);
  if (j != 0) {
    fprintf(stderr, "Not enough memory for spinor fields! Aborting...\n");
    exit(0);
  }

  if (g_proc_id == 0) {
    fprintf(stdout, "# The number of processes is %d \n", g_nproc);
    printf("# The lattice size is %d x %d x %d x %d\n", (int)(T * g_nproc_t), (int)(LX * g_nproc_x),
           (int)(LY * g_nproc_y), (int)(g_nproc_z * LZ));
    printf("# The local lattice size is %d x %d x %d x %d\n", (int)(T), (int)(LX), (int)(LY),
           (int)LZ);
    if (even_odd_flag) {
      printf("# testing the even/odd preconditioned Dirac operator\n");
    } else {
      printf("# testing the standard Dirac operator\n");
    }
    fflush(stdout);
  }

  /* define the geometry */
  geometry();

#ifdef _USE_HALFSPINOR
  j = init_dirac_halfspinor();
  if (j != 0) {
    fprintf(stderr, "Not enough memory for halfspinor fields! Aborting...\n");
    exit(0);
  }
  j = init_dirac_halfspinor32();
  if (j != 0) {
    fprintf(stderr, "Not enough memory for 32-Bit halfspinor fields! Aborting...\n");
    exit(0);
  }
#if (defined _PERSISTENT)
  init_xchange_halffield();
#endif
#endif

  status = check_geometry();
  if (status != 0) {
    fprintf(stderr, "Checking if geometry failed. Unable to proceed.\nAborting....\n");
    exit(1);
  }

  start_ranlux(1, 123456);
  if (startoption == 0) {
    unit_g_gauge_field();  // unit 3x3 colour matrices
  } else {
    random_gauge_field(1, g_gauge_field);
  }

// g_gauge_field[ g_ipt[0][0][0][1] ][0].c00 = 1.0;
// g_gauge_field[ g_ipt[0][0][0][1] ][0].c01 = 0.0;
// g_gauge_field[ g_ipt[0][0][0][1] ][0].c02 = 0.0;
// g_gauge_field[ g_ipt[0][0][0][1] ][0].c10 = 0.0;
// g_gauge_field[ g_ipt[0][0][0][1] ][0].c11 = 1.0;
// g_gauge_field[ g_ipt[0][0][0][1] ][0].c12 = 0.0;
// g_gauge_field[ g_ipt[0][0][0][1] ][0].c20 = 0.0;
// g_gauge_field[ g_ipt[0][0][0][1] ][0].c21 = 0.0;
// g_gauge_field[ g_ipt[0][0][0][1] ][0].c22 = 1.0;

#ifdef TM_USE_MPI
  /*For parallelization: exchange the gaugefield */
  xchange_gauge(g_gauge_field);
#endif

  g_update_gauge_copy = 1;
#ifdef _GAUGE_COPY
  update_backward_gauge(g_gauge_field);
#endif

  init_operators();

  spinor** qphix_out_cb_spinors;
  init_solver_field(&qphix_out_cb_spinors, VOLUME / 2, 2);

  spinor** tmp;
  init_solver_field(&tmp, VOLUME, 2);

  double* difference_l2norm = calloc(no_operators, sizeof(double));

  /* we will loop over the operators defined in the input file
   * and first apply the tmLQCD operator to the test spinor, then
   * the QPhiX operator and then compare */
  for (int op_id = 0; op_id < no_operators; ++op_id) {
    operator* op =& operator_list[op_id];
    op_set_globals(op_id);
    if (op->type == CLOVER || op->type == DBCLOVER) {
      sw_term((const su3**)g_gauge_field, op->kappa, op->c_sw);
      sw_invert(EE, op->mu);
    }
    boundary(g_kappa);
    // check BC
    if (g_proc_id == 0) {
      printf("\nphase_0 = %f + I*%f\n", creal(phase_0), cimag(phase_0));
      printf("phase_1 = %f + I*%f\n", creal(phase_1), cimag(phase_1));
      printf("phase_2 = %f + I*%f\n", creal(phase_2), cimag(phase_2));
      printf("phase_3 = %f + I*%f\n\n", creal(phase_3), cimag(phase_3));
    }
    /* depending on what has been set in the input file, this will create
     * 1) a point source at source_location, spin/colour corresponding to index_start
     * 2) a volume source
     * 3) a time-slice source
     * for the given operator */
    prepare_source(0 /*nstore*/, 0 /*isample*/, index_start, op_id, 0 /*read_source_flag*/,
                   source_location, 12345 /* seed */);

#ifdef TM_USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    tm_t1 = gettime();
    op->applyM(op->prop0, op->prop1, op->sr0, op->sr1);
    // Hopping_Matrix(OE, op->prop0, op->sr1);
    // Hopping_Matrix(EO, op->prop1, op->sr0);
    tm_t2 = gettime();

#ifdef TM_USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    q_t1 = gettime();
    Mfull_qphix(qphix_out_cb_spinors[0], qphix_out_cb_spinors[1], op->sr0, op->sr1, op->type);
    q_t2 = gettime();

    double squarenorm = square_norm(op->sr0, VOLUME / 2, 1) + square_norm(op->sr1, VOLUME / 2, 1);
    if (g_proc_id == 0) {
      printf("  ||source||^2 = %e\n\n", squarenorm);
      fflush(stdout);
    }

    // print L2-norm of result:
    squarenorm = square_norm(op->prop0, VOLUME / 2, 1) + square_norm(op->prop1, VOLUME / 2, 1);
    if (g_proc_id == 0) {
      printf("\n\n");
      printf("# -------------------------------------------- #\n\n");
      printf("# Dslash 1 (tmLQCD) op_type=%d:\n", op->type);
      printf("# ====================\n\n");
      printf("  ||result_1||^2 = %.16e\n", squarenorm);
      printf("  Time for MV mult: %e\n", tm_t2 - tm_t1);
      fflush(stdout);
    }

    // print L2-norm of result:
    squarenorm = square_norm(qphix_out_cb_spinors[0], VOLUME / 2, 1) +
                 square_norm(qphix_out_cb_spinors[1], VOLUME / 2, 1);
    if (g_proc_id == 0) {
      printf("\n\n");
      printf("# -------------------------------------------- #\n\n");
      printf("# Dslash 2 (QPhiX) op_type=%d:\n", op->type);
      printf("# ====================\n\n");
      printf("  ||result_2||^2 = %.16e\n", squarenorm);
      printf("  Time for MV mult: %e\n", q_t2 - q_t1);
      fflush(stdout);
    }

    convert_eo_to_lexic(tmp[0], op->prop0, op->prop1);
    convert_eo_to_lexic(tmp[1], qphix_out_cb_spinors[0], qphix_out_cb_spinors[1]);

    difference_l2norm[op_id] = compare_spinors(tmp[0], tmp[1]);

  }  // for(op_id)

  int failed = 0;
  for (int op_id = 0; op_id < no_operators; op_id++) {
    if (g_proc_id == 0) {
      printf("op_id: %d, |diff|^2 = %.16e\n", op_id, difference_l2norm[op_id]);
    }
    // check if the l2 norm of the difference is tolerable up to rounding
    if (difference_l2norm[op_id] > 2 * g_nproc * VOLUME * DBL_EPSILON) {
      failed = 1;
    }
  }

  free(difference_l2norm);
  finalize_solver(qphix_out_cb_spinors, 2);
  finalize_solver(tmp, 2);
#ifdef TM_USE_OMP
  free_omp_accumulators();
#endif
  free_gauge_field();
  free_geometry_indices();
  free_spinor_field();
  free_moment_field();
#ifdef TM_USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
#endif
  return (failed);
}

double compare_spinors(spinor* s1, spinor* s2) {
#ifdef TM_USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  int coords[4];
  int x, y, z, t, id = 0;
  // list non-zero elements in spinors, but only if the source type was a point source
  // otherwise the output is overwhelming
  if (SourceInfo.type == SRC_TYPE_POINT) {
    if (g_proc_id == 0) printf("\n OUTPUT TMLQCD vs QPHIX SPINOR (tmlQCD format):\n");
    if (g_proc_id == 0)
      printf("g_proc_id | T=%3d LX=%3d LY=%3d LZ=%3d %26s", g_nproc_t * T, g_nproc_x * LX,
             g_nproc_y * LY, g_nproc_z * LZ, " ");
    if (g_proc_id == 0)
      printf("T=%3d LX=%3d LY=%3d LZ=%3d \n", g_nproc_t * T, g_nproc_x * LX, g_nproc_y * LY,
             g_nproc_z * LZ);
    for (int t_global = 0; t_global < g_nproc_t * T; t_global++) {
      coords[0] = t_global / T;
      for (int x_global = 0; x_global < g_nproc_x * LX; x_global++) {
        coords[1] = x_global / LX;
        for (int y_global = 0; y_global < g_nproc_y * LY; y_global++) {
          coords[2] = y_global / LY;
          for (int z_global = 0; z_global < g_nproc_z * LZ; z_global++) {
            coords[3] = z_global / LZ;
#ifdef TM_USE_MPI
            MPI_Cart_rank(g_cart_grid, coords, &id);
#endif
            if (g_proc_id == id) {
              t = t_global - g_proc_coords[0] * T;
              x = x_global - g_proc_coords[1] * LX;
              y = y_global - g_proc_coords[2] * LY;
              z = z_global - g_proc_coords[3] * LZ;
              int idx = g_ipt[t][x][y][z];
              for (int sc = 0; sc < 24; sc++) {
                double e_tmlqcd = spinor_get_elem_linear(&s2[idx],sc/2,sc%2);
                double e_qphix = spinor_get_elem_linear(&s1[idx],sc/2,sc%2);
                
                if (fabs(e_tmlqcd) > 2 * DBL_EPSILON ||
                    fabs(e_qphix) > 2 * DBL_EPSILON) {
                  fflush(stdout);
                  printf("%9d | %5d %6d %6d %6d s%1d c%1d reim%1d : %+5lf %2s", g_proc_id, t_global,
                         x_global, y_global, z_global, sc / 6, (sc / 2) % 3, sc % 2, e_tmlqcd ,
                         " ");
                  printf("%5d %6d %6d %6d s%1d c%1d reim%1d : %+5lf", t_global, x_global, y_global,
                         z_global, sc / 6, (sc / 2) % 3, sc % 2, e_qphix);
                  if (fabs(e_tmlqcd - e_qphix) > 2 * DBL_EPSILON) printf(" !!! ");
                  printf("\n");
                }
              }
            }
#ifdef TM_USE_MPI
            MPI_Barrier(MPI_COMM_WORLD);
#endif
          }  // z
        }    // y
      }      // x
    }        // t
  }          // if( SourceInfo.type == SRC_TYPE_POINT )

#ifdef TM_USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  if (g_proc_id == 0) {
    printf("\n");
    printf("# Comparison tmLQCD vs QPhiX:\n");
    printf("# ===========================\n\n");
  }

  if (g_proc_id == 0) printf("\n OUTPUT TMLQCD vs QPHIX SPINOR (tmlQCD format):\n");
  if (g_proc_id == 0)
    printf("g_proc_id | T=%3d LX=%3d LY=%3d LZ=%3d \n", g_nproc_t * T, g_nproc_x * LX,
           g_nproc_y * LY, g_nproc_z * LZ);
  double squarenorm = diff_and_square_norm(s1, s2, VOLUME);

#ifdef TM_USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  id = 0;
  for (int t_global = 0; t_global < g_nproc_t * T; t_global++) {
    coords[0] = t_global / T;
    for (int x_global = 0; x_global < g_nproc_x * LX; x_global++) {
      coords[1] = x_global / LX;
      for (int y_global = 0; y_global < g_nproc_y * LY; y_global++) {
        coords[2] = y_global / LY;
        for (int z_global = 0; z_global < g_nproc_z * LZ; z_global++) {
          coords[3] = z_global / LZ;
#ifdef TM_USE_MPI
          MPI_Cart_rank(g_cart_grid, coords, &id);
#endif
          if (g_proc_id == id) {
            t = t_global - g_proc_coords[0] * T;
            x = x_global - g_proc_coords[1] * LX;
            y = y_global - g_proc_coords[2] * LY;
            z = z_global - g_proc_coords[3] * LZ;
            int idx = g_ipt[t][x][y][z];
            for (int sc = 0; sc < 24; sc++) {
              double e_diff = spinor_get_elem_linear(&s1[idx],sc/2,sc%2);
              // when a volume source is used, these will be zero up to significant rounding
              // we account for that by the scaling of DBL_EPSILON
              if (fabs(e_diff) > 8 * 24 * DBL_EPSILON) {
                fflush(stdout);
                printf("%9d | %5d %6d %6d %6d s%1d c%1d reim%1d : %+5lf\n", g_proc_id, t_global,
                       x_global, y_global, z_global, sc / 6, (sc / 2) % 3, sc % 2, e_diff);
              }
            }
          }
#ifdef TM_USE_MPI
          MPI_Barrier(MPI_COMM_WORLD);
#endif
        }  // z
      }    // y
    }      // x
  }        // t

  if (g_proc_id == 0) {
    printf("\n  ||result_1 - result_2||^2 = %e\n\n", squarenorm);
    fflush(stdout);
  }
  return squarenorm;
}

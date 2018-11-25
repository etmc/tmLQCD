/***********************************************************************
 *
 * Copyright (C) 2014 Carsten Urbach
 *               2018 Bartosz Kostrzewa
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
 * invert wrapper for using tmLQCD as a library
 *
 * Author: Carsten Urbach
 *         curbach@gmx.de
 *
 *******************************************************************************/

#include <math.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <float.h>
#ifdef TM_USE_MPI
#include <mpi.h>
#endif
#ifdef TM_USE_OMP
#include <omp.h>
#endif
#include "geometry_eo.h"
#include "getopt.h"
#include "git_hash.h"
#include "global.h"
#include "linalg_eo.h"
#ifdef TM_USE_MPI
#include "xchange/xchange.h"
#endif
#ifdef TM_USE_QUDA
#include "quda_interface.h"
#endif
#ifdef TM_USE_QPHIX
#include "qphix_interface.h"
#endif
#include <io/gauge.h>
#include <io/utils.h>
#include "boundary.h"
#include "fatal_error.h"
#include "include/tmLQCD.h"
#include "init/init.h"
#include "invert_eo.h"
#include "linalg/convert_eo_to_lexic.h"
#include "measure_gauge_action.h"
#include "mpi_init.h"
#include "operator.h"
#include "operator/clover_leaf.h"
#include "read_input.h"
#include "sighandler.h"
#include "start.h"
#include "include/tmLQCD.h"
#include "fatal_error.h"
#include "misc_types.h"

#ifdef HAVE_GPU
extern void init_mixedsolve_eo(su3** gf);
extern void init_mixedsolve(su3** gf);
extern void finalize_mixedsolve();
extern void init_gpu_fields(int need_momenta);
extern void finalize_gpu_fields();
#include "GPU/cudadefs.h"
#ifdef TEMPORALGAUGE
#include "temporalgauge.h"
#endif
#endif

#define CONF_FILENAME_LENGTH 500

// we want to be able to modify the boundary phase angles
extern double X0, X1, X2, X3;

static int tmLQCD_invert_initialised = 0;

int tmLQCD_init_parallel_and_read_input(int argc, char* argv[], const int _verbose,
                                        char const * const input_filename){
  verbose = _verbose;
  init_parallel_and_read_input(argc, argv, input_filename);
  return(0);
}

int tmLQCD_invert_init(int argc, char* argv[], const int _verbose, const int external_id) {
  init_critical_globals(TM_PROGRAM_EXTERNAL); 
  
  DUM_DERI = 8;
  DUM_MATRIX = DUM_DERI + 5;
  NO_OF_SPINORFIELDS = DUM_MATRIX + 3;
  // 4 extra fields (corresponding to DUM_MATRIX+0..5) for deg. and ND matrix mult.
  NO_OF_SPINORFIELDS_32 = 6;

  NO_OF_SPINORFIELDS_32 = 6;

  // in read_input.h
  verbose = _verbose;
  g_use_clover_flag = 0;

#ifdef TM_USE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &g_proc_id);
#else
  g_proc_id = 0;
#endif

  /* Read the input file */
  if ((read_input("invert.input")) != 0) {
    fprintf(stderr, "tmLQCD_init_invert: Could not find input file: invert.input\nAborting...");
    return (-1);
  }

#ifndef TM_USE_MPI
  if (subprocess_flag) g_external_id = external_id;
#endif

#ifdef TM_USE_OMP
  init_openmp();
#endif

  tmlqcd_mpi_init(argc, argv);
  g_dbw2rand = 0;
  for (int j = 0; j < no_operators; j++)
    if (!operator_list[j].even_odd_flag) even_odd_flag = 0;

#ifdef _GAUGE_COPY
  int j = init_gauge_field(VOLUMEPLUSRAND, 1);
  j += init_gauge_field_32(VOLUMEPLUSRAND, 1);
#else
  int j = init_gauge_field(VOLUMEPLUSRAND, 0);
  j += init_gauge_field_32(VOLUMEPLUSRAND, 0);
#endif
  if (j != 0) {
    fprintf(stderr, "tmLQCD_init_invert: Not enough memory for gauge_fields! Aborting...\n");
    return (-1);
  }
  j = init_geometry_indices(VOLUMEPLUSRAND);
  if (j != 0) {
    fprintf(stderr, "tmLQCD_init_invert: Not enough memory for geometry indices! Aborting...\n");
    return (-1);
  }
  if (!lowmem_flag) {
    if (even_odd_flag) {
      j = init_spinor_field(VOLUMEPLUSRAND / 2, NO_OF_SPINORFIELDS);
      j += init_spinor_field_32(VOLUMEPLUSRAND / 2, NO_OF_SPINORFIELDS_32);
    } else {
      j = init_spinor_field(VOLUMEPLUSRAND, NO_OF_SPINORFIELDS);
      j += init_spinor_field_32(VOLUMEPLUSRAND, NO_OF_SPINORFIELDS_32);
    }
    if (j != 0) {
      fprintf(stderr, "tmLQCD_init_invert: Not enough memory for spinor fields! Aborting...\n");
      return (-1);
    }
  }

  if (g_cart_id == 0) {
    FILE* parameterfile = parameterfile = fopen("tmLQCD-libwrapper.para", "w");
    write_first_messages(parameterfile, "tmLQCD lib-wrapper", git_hash);
    fclose(parameterfile);
  }

  // define the geometry
  geometry();

  // initialise the operators
  init_operators();

#ifdef HAVE_GPU
  if (usegpu_flag) {
    if (even_odd_flag) {
      init_mixedsolve_eo(g_gauge_field);
    } else {
      init_mixedsolve(g_gauge_field);
    }
#ifdef GPU_DOUBLE
    /*init double fields w/o momenta*/
    init_gpu_fields(0);
#endif
#ifdef TEMPORALGAUGE
    int retval;
    if ((retval = init_temporalgauge(VOLUME, g_gauge_field)) != 0) {
      if (g_proc_id == 0)
        printf("tmLQCD_init_invert: Error while initializing temporal gauge. Aborting...\n");
      exit(200);
    }
#endif
  }  // usegpu_flag
#endif

#ifdef _USE_HALFSPINOR
  if (!lowmem_flag) {
    j = init_dirac_halfspinor();
    if (j != 0) {
      fprintf(stderr, "tmLQCD_init_invert: Not enough memory for halffield! Aborting...\n");
      return (-1);
    }
    j = init_dirac_halfspinor32();
    if (j != 0) {
      fprintf(stderr, "tmLQCD_init_invert: Not enough memory for 32-bit halffield! Aborting...\n");
      return (-1);
    }
#if (defined _PERSISTENT)
    if (even_odd_flag) init_xchange_halffield();
#endif
  }
#endif
  tmLQCD_invert_initialised = 1;
  return (0);
}

int tmLQCD_read_gauge(const int nconfig) {
  char conf_filename[CONF_FILENAME_LENGTH];
  if (!tmLQCD_invert_initialised) {
    fprintf(stderr, "tmLQCD_read_gauge: tmLQCD_inver_init must be called first. Aborting...\n");
    return (-1);
  }

  int n_written =
      snprintf(conf_filename, CONF_FILENAME_LENGTH, "%s.%.4d", gauge_input_filename, nconfig);
  if (n_written < 0 || n_written >= CONF_FILENAME_LENGTH) {
    char error_message[500];
    snprintf(error_message, 500,
             "Encoding error or gauge configuration filename "
             "longer than %d characters! See wrapper/lib_wrapper.c CONF_FILENAME_LENGTH\n",
             CONF_FILENAME_LENGTH);
    fatal_error(error_message, "tmLQCD_read_gauge");
  }
  int j = 0;
  if (g_cart_id == 0) {
    printf("#\n# Trying to read gauge field from file %s.\n", conf_filename);
    fflush(stdout);
  }
  if ((j = read_gauge_field(conf_filename, g_gauge_field)) != 0) {
    fprintf(stderr, "tmLQCD_read_gauge: Error %d while reading gauge field from %s\n ...\n", j,
            conf_filename);
    return (-1);
  }
  if (g_cart_id == 0) {
    printf("# Finished reading gauge field.\n");
    fflush(stdout);
  }

  // set the global nstore parameter
  nstore = nconfig;

#ifdef TM_USE_MPI
  if (!lowmem_flag) {
    xchange_gauge(g_gauge_field);
  }
#endif
  if (!lowmem_flag) {
    convert_32_gauge_field(g_gauge_field_32, g_gauge_field, VOLUMEPLUSRAND);
  }

  double plaquette = measure_plaquette((const su3** const)g_gauge_field) / (6. * VOLUME * g_nproc);
  if (g_cart_id == 0) {
    printf("# The computed plaquette value is %.16e.\n", plaquette);
  }
  return (0);
}

int tmLQCD_invert(double* const propagator, double* const source, const int op_id,
                  const int write_prop) {
  unsigned int index_start = 0;
  g_mu = 0.;

  if (lowmem_flag && g_proc_id == 0) {
    printf(
        "!!! WARNING: you are calling tmLQCD_invert in \'lowmem\' mode.\n Did you make sure that "
        "all required fields are allocated and initialised??\n");
  }

  if (!tmLQCD_invert_initialised) {
    fprintf(stderr, "tmLQCD_invert: tmLQCD_inver_init must be called first. Aborting...\n");
    return (-1);
  }

  if (op_id < 0 || op_id >= no_operators) {
    fprintf(stderr, "tmLQCD_invert: op_id=%d not in valid range. Aborting...\n", op_id);
    return (-1);
  }

  operator_list[op_id].sr0 = g_spinor_field[0];
  operator_list[op_id].sr1 = g_spinor_field[1];
  operator_list[op_id].prop0 = g_spinor_field[2];
  operator_list[op_id].prop1 = g_spinor_field[3];

  zero_spinor_field(operator_list[op_id].prop0, VOLUME / 2);
  zero_spinor_field(operator_list[op_id].prop1, VOLUME / 2);

  // convert to even/odd order
  convert_lexic_to_eo(operator_list[op_id].sr0, operator_list[op_id].sr1, (spinor*)source);

  // invert
  operator_list[op_id].inverter(op_id, index_start, write_prop);

  // convert back to lexicographic order
  convert_eo_to_lexic((spinor*)propagator, operator_list[op_id].prop0, operator_list[op_id].prop1);

  return (0);
}

int tmLQCD_invert_eo(double* const Odd_out, double* const Odd_in, const int op_id){
  unsigned int index_start = 0;
  if (!tmLQCD_invert_initialised) {
    fprintf(stderr, "tmLQCD_invert_eo: tmLQCD_inver_init must be called first. Aborting...\n");
    return (-1);
  }
  if (op_id < 0 || op_id >= no_operators) {
    fprintf(stderr, "tmLQCD_invert_eo: op_id=%d not in valid range. Aborting...\n", op_id);
    return (-2);
  }

  // here we need to provide a way to call the selected solver directly on the
  // e/o precon prepared problem

  return (0);
}

int tmLQCD_finalise() {
#ifdef TM_USE_OMP
  free_omp_accumulators();
#endif

#ifdef HAVE_GPU
  if (usegpu_flag) {
    finalize_mixedsolve();
#ifdef GPU_DOUBLE
    finalize_gpu_fields();
#endif
#ifdef TEMPORALGAUGE
    finalize_temporalgauge();
#endif
  }
#endif

#ifdef TM_USE_QUDA
  _endQuda();
#endif

  free_gauge_field();
  free_geometry_indices();
  if (!lowmem_flag) {
    free_gauge_field_32();
    free_spinor_field();
    free_spinor_field_32();
    free_moment_field();
    free_chi_spinor_field();
  }
#ifdef TM_USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  return (0);
}

int tmLQCD_get_lat_params(tmLQCD_lat_params* params) {
  if (!tmLQCD_invert_initialised) {
    fprintf(stderr, "tmLQCD_get_lat_params: tmLQCD_inver_init must be called first. Aborting...\n");
    return (-1);
  }

  params->LX = LX;
  params->LY = LY;
  params->LZ = LZ;
  params->T = T;
  params->nstore = nstore;
  params->nsave = Nsave;
  params->no_operators = no_operators;
  return (0);
}

int tmLQCD_get_mpi_params(tmLQCD_mpi_params* params) {
  if (!tmLQCD_invert_initialised) {
    fprintf(stderr, "tmLQCD_get_mpi_params: tmLQCD_inver_init must be called first. Aborting...\n");
    return (-1);
  }

  params->nproc = g_nproc;
  params->nproc_t = g_nproc_t;
  params->nproc_x = g_nproc_x;
  params->nproc_y = g_nproc_y;
  params->nproc_z = g_nproc_z;
  params->cart_id = g_cart_id;
  params->proc_id = g_proc_id;
  params->time_rank = g_mpi_time_rank;
  params->omp_num_threads = omp_num_threads;
  params->proc_coords[0] = g_proc_coords[0];
  params->proc_coords[1] = g_proc_coords[1];
  params->proc_coords[2] = g_proc_coords[2];
  params->proc_coords[3] = g_proc_coords[3];
  // g_cart_grid only exists when MPI is in use
#ifdef TM_USE_MPI
  params->cart_grid = g_cart_grid;
#else
  params->cart_grid = -1;
#endif
  return(0);
}

int tmLQCD_get_gauge_field_pointer(double** gf) {
  if (!tmLQCD_invert_initialised) {
    fprintf(
        stderr,
        "tmLQCD_get_gauge_field_pointer: tmLQCD_invert_init must be called first. Aborting...\n");
    return (-1);
  }
#ifdef TM_USE_MPI
  xchange_gauge(g_gauge_field);
#endif
  if (!lowmem_flag) {
    convert_32_gauge_field(g_gauge_field_32, g_gauge_field, VOLUMEPLUSRAND);
  }

  *gf = (double*)g_gauge_field[0];

  return (0);
}

int tmLQCD_get_op_params(tmLQCD_op_params* params, const int op_id) {
  if (!tmLQCD_invert_initialised) {
    fprintf(stderr, "tmLQCD_get_op_params: tmLQCD_invert_init must be called first. Aborting...\n");
    return (-1);
  }
  if (op_id < 0 || op_id >= no_operators) {
    fprintf(stderr, "tmLQCD_get_op_params: op_id=%d not in valid range. Aborting...\n", op_id);
    return (-2);
  }

  params->kappa = operator_list[op_id].kappa;
  params->mu = 0.5 * operator_list[op_id].mu / operator_list[op_id].kappa;
  params->mubar = 0.5 * operator_list[op_id].mubar / operator_list[op_id].kappa;
  params->epsbar = 0.5 * operator_list[op_id].epsbar / operator_list[op_id].kappa;
  params->c_sw = operator_list[op_id].c_sw;

  params->theta_x = X1;
  params->theta_y = X2;
  params->theta_z = X3;
  params->theta_t = X0;
  return (0);
}

int tmLQCD_set_op_params(tmLQCD_op_params const* const params, const int op_id) {
  if (!tmLQCD_invert_initialised) {
    fprintf(stderr, "tmLQCD_set_op_params: tmLQCD_invert_init must be called first. Aborting...\n");
    return (-1);
  }
  if (op_id < 0 || op_id >= no_operators) {
    fprintf(stderr, "tmLQCD_set_op_params: op_id=%d not in valid range. Aborting...\n", op_id);
    return (-2);
  }

  operator_list[op_id].kappa = params->kappa;
  operator_list[op_id].mu = 2 * params->kappa * params->mu;
  operator_list[op_id].mubar = 2 * params->kappa * params->mubar;
  operator_list[op_id].epsbar = 2 * params->kappa * params->epsbar;
  operator_list[op_id].c_sw = params->c_sw;

  X1 = params->theta_x;
  X2 = params->theta_y;
  X3 = params->theta_z;
  X0 = params->theta_t;
  boundary(params->kappa);
  return (0);
}

#ifdef TM_USE_QPHIX
int tmLQCD_invert_qphix_direct(double * const Odd_out, double * const Odd_in, const int op_id){
  static double clover_term_c_sw = -1.0;
  static double clover_term_kappa = -1.0;
  static double inv_clover_term_mu = 0.0;
  
  if (!tmLQCD_invert_initialised) {
    fprintf(stderr, "tmLQCD_invert: tmLQCD_inver_init must be called first. Aborting...\n");
    return (-1);
  }

  if (op_id < 0 || op_id >= no_operators) {
    fprintf(stderr, "tmLQCD_invert: op_id=%d not in valid range. Aborting...\n", op_id);
    return (-2);
  }

  op_backup_restore_globals(TM_BACKUP_GLOBALS);
  op_set_globals(op_id);
  boundary(g_kappa);
  
  // if this is a clover operator, we check if the clover term and its inverse
  // which might be currently in memory are consistent with the requested parameters
  // if not, we recompute them
  if( operator_list[op_id].type == CLOVER ){
    if( (fabs(g_c_sw - clover_term_c_sw) > 2*DBL_EPSILON ) ||
        (fabs(g_kappa - clover_term_kappa) > 2*DBL_EPSILON ) ){
      sw_term((const su3**)g_gauge_field,
              g_kappa,
              g_c_sw);
      clover_term_kappa = g_kappa;
      clover_term_c_sw = g_c_sw;
    }
    if( fabs(g_mu - inv_clover_term_mu) > 2*DBL_EPSILON ){
      // in the hopping matrix and elsewhere "even_even" is called "EE"
      // but it's only defined in operator/HoppingMatrix.h, we don't want to
      // include that here
      const int even_even = 0;
      sw_invert(even_even, g_mu);
      inv_clover_term_mu = g_mu;
    }
  }
  
  int niter = invert_eo_qphix_oneflavour((spinor *const) Odd_out, (spinor* const) Odd_in,
                                         operator_list[op_id].maxiter,
                                         operator_list[op_id].eps_sq,
                                         operator_list[op_id].solver,
                                         operator_list[op_id].rel_prec,
                                         operator_list[op_id].solver_params,
                                         operator_list[op_id].sloppy_precision,
                                         operator_list[op_id].compression_type
                                         );
  op_backup_restore_globals(TM_RESTORE_GLOBALS);
  return(niter);
}
#endif

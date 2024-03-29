/***********************************************************************
 *
 * Copyright (C) 2012 Carsten Urbach, Albert Deuzeman, Bartosz Kostrzewa
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
 * naive pion correlator for twisted mass QCD
 *
 *******************************************************************************/

#include <lime.h>
#ifdef HAVE_CONFIG_H
#include "tmlqcd_config.h"
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <signal.h>
#ifdef TM_USE_MPI
#include <mpi.h>
#endif
#ifdef TM_USE_OMP
#include <omp.h>
#endif
#include "global.h"
#include "git_hash.h"
#include "read_input.h"
#include "getopt.h"
#include "geometry_eo.h"
#include "start.h"
#include "measure_gauge_action.h"
#ifdef TM_USE_MPI
#include "xchange/xchange.h"
#endif
#include "sighandler.h"
#include "mpi_init.h"
#include "sighandler.h"
#include "boundary.h"
#include "init/init.h"
#include "monomial/monomial.h"
#include "phmc.h"
#include "block.h"
#include "operator.h"
#include "solver/generate_dfl_subspace.h"
#include "io/gauge.h"
#include "meas/measurements.h"
#ifdef TM_USE_QUDA
#include "quda_interface.h"
#endif

#define CONF_FILENAME_LENGTH 500

extern int nstore;
int check_geometry();

static void usage(const tm_ExitCode_t exit_code);
static void process_args(int argc, char *argv[], char ** input_filename, char ** filename);
static void set_default_filenames(char ** input_filename, char ** filename);

int main(int argc, char *argv[])
{
  FILE *parameterfile = NULL;
  int j, i, ix = 0, isample = 0, op_id = 0;
  char datafilename[206];
  char parameterfilename[206];
  char conf_filename[CONF_FILENAME_LENGTH];
  char * input_filename = NULL;
  char * filename = NULL;
  double plaquette_energy;

  init_critical_globals(TM_PROGRAM_OFFLINE_MEASUREMENT);  

#ifdef _KOJAK_INST
#pragma pomp inst init
#pragma pomp inst begin(main)
#endif

#if (defined SSE || defined SSE2 || SSE3)
  signal(SIGILL, &catch_ill_inst);
#endif

  DUM_DERI = 8;
  DUM_MATRIX = DUM_DERI + 5;
#if ((defined BGL && defined XLC) || defined _USE_TSPLITPAR)
  NO_OF_SPINORFIELDS = DUM_MATRIX + 3;
#else
  NO_OF_SPINORFIELDS = DUM_MATRIX + 3;
#endif

  verbose = 0;
  g_use_clover_flag = 0;

  process_args(argc,argv,&input_filename,&filename);
  set_default_filenames(&input_filename, &filename);
  init_parallel_and_read_input(argc, argv, input_filename);

  /* this DBW2 stuff is not needed for the inversion ! */
  if (g_dflgcr_flag == 1) {
    even_odd_flag = 0;
  }
  if (Nsave == 0) {
    Nsave = 1;
  }

  if (g_running_phmc) {
    NO_OF_SPINORFIELDS = DUM_MATRIX + 8;
  }

  tmlqcd_mpi_init(argc, argv);

  /* starts the single and double precision random number */
  /* generator                                            */
  start_ranlux(rlxd_level, random_seed^nstore);
  
  /* we need to make sure that we don't have even_odd_flag = 1 */
  /* if any of the operators doesn't use it                    */
  /* in this way even/odd can still be used by other operators */
  for(j = 0; j < no_operators; j++) if(!operator_list[j].even_odd_flag) even_odd_flag = 0;

#ifndef TM_USE_MPI
  g_dbw2rand = 0;
#endif

#ifdef _GAUGE_COPY
  j = init_gauge_field(VOLUMEPLUSRAND + g_dbw2rand, 1);
#else
  j = init_gauge_field(VOLUMEPLUSRAND + g_dbw2rand, 0);
#endif
  if (j != 0) {
    fprintf(stderr, "Not enough memory for gauge_fields! Aborting...\n");
    exit(-1);
  }
  j = init_geometry_indices(VOLUMEPLUSRAND + g_dbw2rand);
  if (j != 0) {
    fprintf(stderr, "Not enough memory for geometry indices! Aborting...\n");
    exit(-1);
  }
  if (no_monomials > 0) {
    if (even_odd_flag) {
      j = init_monomials(VOLUMEPLUSRAND / 2, even_odd_flag);
    }
    else {
      j = init_monomials(VOLUMEPLUSRAND, even_odd_flag);
    }
    if (j != 0) {
      fprintf(stderr, "Not enough memory for monomial pseudo fermion fields! Aborting...\n");
      exit(-1);
    }
  }
  if (even_odd_flag) {
    j = init_spinor_field(VOLUMEPLUSRAND / 2, NO_OF_SPINORFIELDS);
  }
  else {
    j = init_spinor_field(VOLUMEPLUSRAND, NO_OF_SPINORFIELDS);
  }
  if (j != 0) {
    fprintf(stderr, "Not enough memory for spinor fields! Aborting...\n");
    exit(-1);
  }

  if (g_running_phmc) {
    j = init_chi_spinor_field(VOLUMEPLUSRAND / 2, 20);
    if (j != 0) {
      fprintf(stderr, "Not enough memory for PHMC Chi fields! Aborting...\n");
      exit(-1);
    }
  }

  g_mu = g_mu1;

  if (g_cart_id == 0) {
    /*construct the filenames for the observables and the parameters*/
    strncpy(datafilename, filename, 200);
    strcat(datafilename, ".data");
    strncpy(parameterfilename, filename, 200);
    strcat(parameterfilename, ".para");

    parameterfile = fopen(parameterfilename, "w");
    write_first_messages(parameterfile, "invert", git_hash);
    fclose(parameterfile);
  }

  /* define the geometry */
  geometry();
  int status = check_geometry();

  if (status != 0) {
    fprintf(stderr, "Checking of geometry failed. Unable to proceed.\nAborting....\n");
    exit(1);
  }

  /* define the boundary conditions for the fermion fields */
  boundary(g_kappa);

  phmc_invmaxev = 1.;

  init_operators();
  
  /* list and initialize measurements*/
  if(g_proc_id == 0) {
    printf("\n");
    for(int j = 0; j < no_measurements; j++) {
      printf("# measurement id %d, type = %d\n", j, measurement_list[j].type);
    }
  }
  init_measurements();  

  /* this could be maybe moved to init_operators */
#ifdef _USE_HALFSPINOR
  j = init_dirac_halfspinor();
  if (j != 0) {
    fprintf(stderr, "Not enough memory for halffield! Aborting...\n");
    exit(-1);
  }
  if (g_sloppy_precision_flag == 1) {
    j = init_dirac_halfspinor32();
    if (j != 0)
    {
      fprintf(stderr, "Not enough memory for 32-bit halffield! Aborting...\n");
      exit(-1);
    }
  }
#  if (defined _PERSISTENT)
  if (even_odd_flag)
    init_xchange_halffield();
#  endif
#endif

  for (j = 0; j < Nmeas; j++) {
    int n_written = snprintf(conf_filename, CONF_FILENAME_LENGTH, "%s.%.4d", gauge_input_filename, nstore);
    if( n_written < 0 || n_written > CONF_FILENAME_LENGTH ){
      char error_message[500];
      snprintf(error_message,
               500,
               "Encoding error or gauge configuration filename "
               "longer than %d characters! See offline_measurement.c CONF_FILENAME_LENGTH\n", 
               CONF_FILENAME_LENGTH);
      fatal_error(error_message, "offline_measurement.c");
    }
    if (g_cart_id == 0) {
      printf("#\n# Trying to read gauge field from file %s in %s precision.\n",
            conf_filename, (gauge_precision_read_flag == 32 ? "single" : "double"));
      fflush(stdout);
    }
    if( (i = read_gauge_field(conf_filename,g_gauge_field)) !=0) {
      fprintf(stderr, "Error %d while reading gauge field from %s\n Aborting...\n", i, conf_filename);
      exit(-2);
    }
  
    if (g_cart_id == 0) {
      printf("# Finished reading gauge field.\n");
      fflush(stdout);
    }

#ifdef TM_USE_MPI
    xchange_gauge(g_gauge_field);
    update_tm_gauge_exchange(&g_gauge_state);
#endif

    /*compute the energy of the gauge field*/
    plaquette_energy = measure_plaquette( (const su3** const) g_gauge_field);

    if (g_cart_id == 0) {
      printf("# The computed plaquette value is %e.\n", plaquette_energy / (6.*VOLUME*g_nproc));
      fflush(stdout);
    }

    if (g_cart_id == 0) {
      fprintf(stdout, "#\n"); /*Indicate starting of the operator part*/
    }

    
    /* offline measurements */
    measurement * meas;
    for(int imeas = 0; imeas < no_measurements; imeas++){
      meas = &measurement_list[imeas];
      if (g_proc_id == 0) {
        fprintf(stdout, "#\n# Beginning offline measurement.\n");
      }
      meas->measurefunc(nstore, imeas, even_odd_flag);
    }      
    nstore += Nsave;
  }

#ifdef TM_USE_OMP
  free_omp_accumulators();
#endif
#ifdef TM_USE_QUDA
  _endQuda();
#endif

  free_blocks();
  free_dfl_subspace();
  free_geometry_indices();
  free_spinor_field();

  free_chi_spinor_field();

  free(filename);
  free(input_filename);

#ifdef TM_USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
#endif
  return(0);
  
  
#ifdef _KOJAK_INST
#pragma pomp inst end(main)
#endif
}

static void usage(const tm_ExitCode_t exit_code)
{ 
  if( g_proc_id == 0 ){
    fprintf(stdout, "Offline version of the online measurements for twisted mass QCD\n");
    fprintf(stdout, "Version %s \n\n", TMLQCD_PACKAGE_VERSION);
    fprintf(stdout, "Please send bug reports to %s\n", TMLQCD_PACKAGE_BUGREPORT);
    fprintf(stdout, "Usage:   invert [options]\n");
    fprintf(stdout, "Options: [-f input-filename]\n");
    fprintf(stdout, "         [-v] more verbosity\n");
    fprintf(stdout, "         [-V] print version information and exit\n");
    fprintf(stdout, "         [-m level] request MPI thread level 'single' or 'multiple' (default: 'single')\n");
    fprintf(stdout, "         [-h|-? this help]\n");
  }
  exit(exit_code);
}

static void process_args(int argc, char *argv[], char ** input_filename, char ** filename) {
  int c;
  while ((c = getopt(argc, argv, "h?vVf:o:")) != -1) {
    switch (c) {
      case 'f':
        *input_filename = calloc(200, sizeof(char));
        strncpy(*input_filename, optarg, 200);
        break;
      case 'v':
        verbose = 1;
        break;
      case 'V':
        if(g_proc_id == 0) {
          fprintf(stdout,"%s %s\n",TMLQCD_PACKAGE_STRING,git_hash);
        }
        exit(TM_EXIT_SUCCESS);
        break;
      case 'm':
        if( !strcmp(optarg, "single") ){
          g_mpi_thread_level = TM_MPI_THREAD_SINGLE;
        } else if ( !strcmp(optarg, "multiple") ) {
          g_mpi_thread_level = TM_MPI_THREAD_MULTIPLE;
        } else {
          tm_debug_printf(0, 0, "[offline_measurement process_args]: invalid input for -m command line argument\n");
          usage(TM_EXIT_INVALID_CMDLINE_ARG);
        }
        break;
      case 'h':
      case '?':
      default:
        usage(TM_EXIT_SUCCESS);
        break;
    }
  }
}

static void set_default_filenames(char ** input_filename, char ** filename) {
  if( *input_filename == NULL ) {
    *input_filename = calloc(28, sizeof(char));
    strcpy(*input_filename,"offline_measurement.input");
  }
  
  if( *filename == NULL ) {
    *filename = calloc(7, sizeof(char));
    strcpy(*filename,"output");
  } 
}


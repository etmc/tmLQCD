/***********************************************************************
 *
 * Copyright (C) 2021 Ferenc Pittler
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
 *
 * Computing the gauge force only
 *
 * Author: Ferenc Pittler
 *         f.pittler@cyi.ac.cy
 *         
 *******************************************************************************/
#include "lime.h"
#if HAVE_CONFIG_H
#include<tmlqcd_config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include <signal.h>
#include <unistd.h>
#ifdef TM_USE_MPI
#include <mpi.h>
#endif
#ifdef TM_USE_OMP
#include <omp.h>
#endif
#include "global.h"
#include "git_hash.h"
#include "io/params.h"
#include "io/gauge.h"
#include "getopt.h"
#include "ranlxd.h"
#include "geometry_eo.h"
#include "start.h"
#include "measure_gauge_action.h"
#include "measure_rectangles.h"
#ifdef TM_USE_MPI
#include "xchange/xchange.h"
#endif
#include "read_input.h"
#include "mpi_init.h"
#include "sighandler.h"
#include "update_tm.h"
#include "init/init.h"
#include "test/check_geometry.h"
#include "boundary.h"
#include "phmc.h"
#include "solver/solver.h"
#include "monomial/monomial.h"
#include "integrator.h"
#include "sighandler.h"
#include "meas/measurements.h"
#ifdef DDalphaAMG
#include "DDalphaAMG_interface.h"
#endif
#ifdef TM_USE_QUDA
#  include "quda_interface.h"
#endif

extern int nstore;

int const rlxdsize = 105;

static void usage(const tm_ExitCode_t exit_code);
static void process_args(int argc, char *argv[], char ** input_filename, char ** filename);
static void set_default_filenames(char ** input_filename, char ** filename);

int main(int argc,char *argv[]) {

  FILE *parameterfile=NULL, *countfile=NULL;
  char *filename = NULL;
  char datafilename[206];
  char parameterfilename[206];
  char gauge_filename[50];
  char nstore_filename[50];
  char tmp_filename[50];
  char *input_filename = NULL;
  int status = 0, accept = 0;
  int j,ix,mu, trajectory_counter=0;
  unsigned int const io_max_attempts = 5; /* Make this configurable? */
  unsigned int const io_timeout = 5; /* Make this configurable? */
  struct timeval t1;

  /* Energy corresponding to the Gauge part */
  double plaquette_energy = 0., rectangle_energy = 0.;
  /* Acceptance rate */
  int Rate=0;
  /* Do we want to perform reversibility checks */
  /* See also return_check_flag in read_input.h */
  int return_check = 0;

  paramsXlfInfo *xlfInfo;

/* For online measurements */
  measurement * meas;
  int imeas;

  init_critical_globals(TM_PROGRAM_HMC_TM);  
  
#ifdef _KOJAK_INST
#pragma pomp inst init
#pragma pomp inst begin(main)
#endif

#if (defined SSE || defined SSE2 || SSE3)
  signal(SIGILL,&catch_ill_inst);
#endif

  strcpy(gauge_filename,"conf.save");
  strcpy(nstore_filename,".nstore_counter");
  strcpy(tmp_filename, ".conf.tmp");

  verbose = 1;
  g_use_clover_flag = 0;

  process_args(argc,argv,&input_filename,&filename);
  set_default_filenames(&input_filename,&filename);

  init_parallel_and_read_input(argc, argv, input_filename);

  DUM_DERI = 4;
  DUM_MATRIX = DUM_DERI+7;
  if(g_running_phmc) {
    NO_OF_SPINORFIELDS = DUM_MATRIX+8;
  }
  else {
    NO_OF_SPINORFIELDS = DUM_MATRIX+6;
  }
  DUM_BI_DERI = 6;
  DUM_BI_SOLVER = DUM_BI_DERI+7;

  DUM_BI_MATRIX = DUM_BI_SOLVER+6;
  NO_OF_BISPINORFIELDS = DUM_BI_MATRIX+6;
  
  //4 extra fields (corresponding to DUM_MATRIX+0..5) for deg. and ND matrix mult.
  NO_OF_SPINORFIELDS_32 = 6;
  
  tmlqcd_mpi_init(argc, argv);

#ifndef TM_USE_MPI
  g_dbw2rand = 0;
#endif
  
  
  g_mu = g_mu1;
  
#ifdef _GAUGE_COPY
  status = init_gauge_field(VOLUMEPLUSRAND + g_dbw2rand, 1);
  status += init_gauge_field_32(VOLUMEPLUSRAND + g_dbw2rand, 1);
#else
  status = init_gauge_field(VOLUMEPLUSRAND + g_dbw2rand, 0);
  status += init_gauge_field_32(VOLUMEPLUSRAND + g_dbw2rand, 0);   
#endif
  /* need temporary gauge field for gauge reread checks and in update_tm */
  status += init_gauge_tmp(VOLUME);

  status += init_gauge_fg(VOLUME);

  if (status != 0) {
    fprintf(stderr, "Not enough memory for gauge_fields! Aborting...\n");
    exit(0);
  }
  j = init_geometry_indices(VOLUMEPLUSRAND + g_dbw2rand);
  if (j != 0) {
    fprintf(stderr, "Not enough memory for geometry_indices! Aborting...\n");
    exit(0);
  }
  j = init_moment_field(VOLUME, VOLUMEPLUSRAND + g_dbw2rand);
  if (j != 0) {
    fprintf(stderr, "Not enough memory for moment fields! Aborting...\n");
    exit(0);
  }

  if(g_running_phmc) {
    j = init_bispinor_field(VOLUME/2, NO_OF_BISPINORFIELDS);
    if (j!= 0) {
      fprintf(stderr, "Not enough memory for bi-spinor fields! Aborting...\n");
      exit(0);
    }
  }

  /*construct the filenames for the observables and the parameters*/
  strncpy(datafilename,filename,200);  
  strcat(datafilename,".data");
  strncpy(parameterfilename,filename,200);  
  strcat(parameterfilename,".para");

  if(g_proc_id == 0){
    parameterfile = fopen(parameterfilename, "a");
    write_first_messages(parameterfile, "hmc", git_hash);
  }

  /* define the geometry */
  geometry();

  /* define the boundary conditions for the fermion fields */
  boundary(g_kappa);

  status = check_geometry();

  if (status != 0) {
    fprintf(stderr, "Checking of geometry failed. Unable to proceed.\nAborting....\n");
    exit(1);
  }


  /* Initialise random number generator */
  start_ranlux(rlxd_level, random_seed^trajectory_counter);

  /* hot */
  random_gauge_field(reproduce_randomnumber_flag, g_gauge_field);

  /*For parallelization: exchange the gaugefield */
#ifdef TM_USE_MPI
  xchange_gauge(g_gauge_field);
  update_tm_gauge_exchange(&g_gauge_state);
#endif
    
  
    
  if(even_odd_flag) {
    j = init_monomials(VOLUMEPLUSRAND/2, even_odd_flag);
  }
  else {
    j = init_monomials(VOLUMEPLUSRAND, even_odd_flag);
  }
  if (j != 0) {
    fprintf(stderr, "Not enough memory for monomial pseudo fermion fields! Aborting...\n");
    exit(0);
  }

  init_integrator();

  if(g_proc_id == 0) {
    for(j = 0; j < no_monomials; j++) {
      printf("# monomial id %d type = %d timescale %d\n", j, monomial_list[j].type, monomial_list[j].timescale);
    }
  }

  plaquette_energy = measure_plaquette( (const su3**) g_gauge_field);
  if(g_rgi_C1 > 0. || g_rgi_C1 < 0.) {
    rectangle_energy = measure_rectangles( (const su3**) g_gauge_field);
    if(g_proc_id == 0){
      fprintf(parameterfile,"# Computed rectangle value: %14.12f.\n",rectangle_energy/(12.*VOLUME*g_nproc));
    }
  }
  //eneg = g_rgi_C0 * plaquette_energy + g_rgi_C1 * rectangle_energy;

  if(g_proc_id == 0) {
    fprintf(parameterfile,"# Computed plaquette value: %14.12f.\n", plaquette_energy/(6.*VOLUME*g_nproc));
    printf("# Computed plaquette value: %14.12f.\n", plaquette_energy/(6.*VOLUME*g_nproc));
    fclose(parameterfile);
  }

  /* set df0 to zero */
  for(ix = 0; ix < VOLUMEPLUSRAND; ix++){
    for(mu=0; mu<4; mu++){
      df0[ix][mu].d1=0.;
      df0[ix][mu].d2=0.;
      df0[ix][mu].d3=0.;
      df0[ix][mu].d4=0.;
      df0[ix][mu].d5=0.;
      df0[ix][mu].d6=0.;
      df0[ix][mu].d7=0.;
      df0[ix][mu].d8=0.;
    }
  }
  double c1=0.;


    if(g_proc_id == 0) {
      printf("#\n# Starting trajectory no %d\n", trajectory_counter);
    }

    return_check = return_check_flag && (trajectory_counter%return_check_interval == 0);

    hamiltonian_field_t hf;

    hf.gaugefield = g_gauge_field;
    hf.momenta = moment;
    hf.derivative = df0;
    hf.update_gauge_copy = g_update_gauge_copy;
    hf.traj_counter = 0;
    integrator_set_fields(&hf);

    monomial_list[0].c0=1-8*c1;
    monomial_list[0].c1=c1;
    if (c1!=0){
     monomial_list[0].use_rectangles=1;
    }
    else{
      monomial_list[0].use_rectangles=0;
    }
    g_beta=1;
    
    monomial_list[0].derivativefunction = &gauge_derivative;
    monomial_list[0].derivativefunction(0, &hf);

    double sum2 =0.;
   for(int i = 0; i < VOLUME; i++) {
     for(int mu = 0; mu < 4; mu++) {
       sum2 = _su3adj_square_norm(hf.derivative[i][mu]);
       printf("i = %d  mu = %d sum %e \n", i, mu, sum2);
     }
   }



#ifdef TM_USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif

#ifdef TM_USE_OMP
  free_omp_accumulators();
#endif
  free_gauge_tmp();
  free_gauge_field();
  free_gauge_field_32();  
  free_geometry_indices();
  free_spinor_field();
  free_spinor_field_32();  
  free_moment_field();
  free_monomials();
  if(g_running_phmc) {
    free_bispinor_field();
    free_chi_spinor_field();
  }
  free(input_filename);
  free(filename);
  free(SourceInfo.basename);
  free(PropInfo.basename);
#ifdef TM_USE_QUDA
  _endQuda();
#endif
#ifdef TM_USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
#endif
  return(0);
#ifdef _KOJAK_INST
#pragma pomp inst end(main)
#endif
}

static void usage(const tm_ExitCode_t exit_code){
  if(g_proc_id == 0){
    fprintf(stdout, "HMC for Wilson twisted mass QCD\n");
    fprintf(stdout, "Version %s \n\n", TMLQCD_PACKAGE_VERSION);
    fprintf(stdout, "Please send bug reports to %s\n", TMLQCD_PACKAGE_BUGREPORT);
    fprintf(stdout, "Usage:   hmc_tm [options]\n");
    fprintf(stdout, "Options: [-f input-filename]  default: hmc.input\n");
    fprintf(stdout, "         [-o output-filename] default: output\n");
    fprintf(stdout, "         [-v] more verbosity\n");
    fprintf(stdout, "         [-V] print version information and exit\n");
    fprintf(stdout, "         [-m level] request MPI thread level 'single' or 'multiple' (default: 'single')\n");
    fprintf(stdout, "         [-h|-? this help]\n");
  }
  exit(exit_code);
}

static void process_args(int argc, char *argv[], char ** input_filename, char ** filename) {
  int c;
  while ((c = getopt(argc, argv, "h?vVf:o:m:")) != -1) {
    switch (c) {
      case 'f':
        *input_filename = calloc(200, sizeof(char));
        strncpy(*input_filename, optarg, 200);
        break;
      case 'o':
        *filename = calloc(200, sizeof(char));
        strncpy(*filename, optarg, 200);
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
          tm_debug_printf(0, 0, "[hmc_tm process_args]: invalid input for -m command line argument\n");
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
    *input_filename = calloc(13, sizeof(char));
    strcpy(*input_filename,"hmc.input");
  }
  
  if( *filename == NULL ) {
    *filename = calloc(7, sizeof(char));
    strcpy(*filename,"output");
  } 
}


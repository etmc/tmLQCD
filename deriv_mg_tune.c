/***********************************************************************
 *
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008 Carsten Urbach
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
 * Hybrid-Monte-Carlo for twisted mass QCD
 *
 * Author: Carsten Urbach
 *         urbach@physik.fu-berlin.de
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

  init_critical_globals(TM_PROGRAM_DERIV_MG_TUNE);  
  
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
  tm_stopwatch_push(&g_timers, "DERIV_MG_TUNE", "");

  if(nstore == -1) {
    countfile = fopen(nstore_filename, "r");
    if(countfile != NULL) {
      j = fscanf(countfile, "%d %d %s\n", &nstore, &trajectory_counter, gauge_input_filename);
      if(j < 1) nstore = 0;
      if(j < 2) trajectory_counter = 0;
      fclose(countfile);
    }
    else {
      nstore = 0;
      trajectory_counter = 0;
    }
  }
  
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
  if(even_odd_flag) {
    j = init_spinor_field(VOLUMEPLUSRAND/2, NO_OF_SPINORFIELDS);
    j += init_spinor_field_32(VOLUMEPLUSRAND/2, NO_OF_SPINORFIELDS_32);      
  }
  else {
    j = init_spinor_field(VOLUMEPLUSRAND, NO_OF_SPINORFIELDS);
    j += init_spinor_field_32(VOLUMEPLUSRAND, NO_OF_SPINORFIELDS_32);    
  }
  if (j != 0) {
    fprintf(stderr, "Not enough memory for spinor fields! Aborting...\n");
    exit(0);
  }
  if(even_odd_flag) {
    j = init_csg_field(VOLUMEPLUSRAND/2);
  }
  else {
    j = init_csg_field(VOLUMEPLUSRAND);
  }
  if (j != 0) {
    fprintf(stderr, "Not enough memory for csg fields! Aborting...\n");
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

  /* define the geometry */
  geometry();

  /* define the boundary conditions for the fermion fields */
  boundary(g_kappa);

  status = check_geometry();

  if (status != 0) {
    fprintf(stderr, "Checking of geometry failed. Unable to proceed.\nAborting....\n");
    exit(1);
  }


#ifdef _USE_HALFSPINOR
  j = init_dirac_halfspinor();
  if (j!= 0) {
    fprintf(stderr, "Not enough memory for halffield! Aborting...\n");
    exit(-1);
  }

  j = init_dirac_halfspinor32();
  if (j != 0)
  {
    fprintf(stderr, "Not enough memory for 32-bit halffield! Aborting...\n");
    exit(-1);
  } 
  
#  if (defined _PERSISTENT)
  init_xchange_halffield();
#  endif
#endif

  /* Initialise random number generator */
  start_ranlux(rlxd_level, random_seed^trajectory_counter);

  /* Set up the gauge field */
  /* continue and restart */
  if(startoption==3 || startoption == 2) {
    if(g_proc_id == 0) {
      printf("# Trying to read gauge field from file %s in %s precision.\n",
            gauge_input_filename, (gauge_precision_read_flag == 32 ? "single" : "double"));
      fflush(stdout);
    }
    if( (status = read_gauge_field(gauge_input_filename,g_gauge_field)) != 0) {
      fprintf(stderr, "Error %d while reading gauge field from %s\nAborting...\n", status, gauge_input_filename);
      exit(-2);
    }

    if (g_proc_id == 0){
      printf("# Finished reading gauge field.\n");
      fflush(stdout);
    }
  }
  else if (startoption == 1) {
    /* hot */
    random_gauge_field(reproduce_randomnumber_flag, g_gauge_field);
  }
  else if(startoption == 0) {
    /* cold */
    unit_g_gauge_field();
  }

  /*For parallelization: exchange the gaugefield */
#ifdef TM_USE_MPI
  xchange_gauge(g_gauge_field);
  update_tm_gauge_exchange(&g_gauge_state);
#endif
    
  /*Convert to a 32 bit gauge field, after xchange*/
  convert_32_gauge_field(g_gauge_field_32, g_gauge_field, VOLUMEPLUSRAND + g_dbw2rand);
#ifdef TM_USE_MPI
  update_tm_gauge_exchange(&g_gauge_state_32);
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

  if( no_monomials > 1 ){
    int have_cltrlog = 0;
    for ( int mnl = 0; mnl < no_monomials; mnl++ ){
      have_cltrlog = (monomial_list[mnl].type == CLOVERTRLOG) || have_cltrlog;
    }
    if( (have_cltrlog && no_monomials > 2) || (!have_cltrlog && no_monomials > 1) ){ 
      fprintf(stderr, "deriv_mg_tune: only one determinant or determinant ratio monomial may be defined! Aborting...\n");
      exit(129);
    }
  }

  const monomial * const mnl = &monomial_list[0];
  const int mnl_type = mnl->type;
  if( !(mnl_type == DET || mnl_type == DETRATIO || mnl_type == CLOVERDET || mnl_type == CLOVERDETRATIO) ){
    fprintf(stderr, "deriv_mg_tune: only DET, DETRATIO, CLOVERDET and CLOVERDETRATIO monomials are supported! Aborting...\n");
    exit(130);
  }

  if(g_proc_id == 0) {
    for(j = 0; j < no_monomials; j++) {
      printf("# monomial id %d type = %d timescale %d\n", j, monomial_list[j].type, monomial_list[j].timescale);
    }
  }

  plaquette_energy = measure_plaquette( (const su3**) g_gauge_field);
  if(g_proc_id == 0) {
    printf("# Computed plaquette value: %14.12f.\n", plaquette_energy/(6.*VOLUME*g_nproc));
  }

  /* set ddummy to zero */
  for(ix = 0; ix < VOLUMEPLUSRAND; ix++){
    for(mu=0; mu<4; mu++){
      ddummy[ix][mu].d1=0.;
      ddummy[ix][mu].d2=0.;
      ddummy[ix][mu].d3=0.;
      ddummy[ix][mu].d4=0.;
      ddummy[ix][mu].d5=0.;
      ddummy[ix][mu].d6=0.;
      ddummy[ix][mu].d7=0.;
      ddummy[ix][mu].d8=0.;
    }
  }

  hamiltonian_field_t hf;
  hf.gaugefield = g_gauge_field;
  hf.momenta = moment;
  hf.derivative = df0; 

  if( mnl_type == CLOVERDET || mnl_type == CLOVERDETRATIO ){
    mnl_backup_restore_globals(TM_BACKUP_GLOBALS);
    g_mu = mnl->mu;
    g_mu3 = mnl->rho;
    g_c_sw = mnl->c_sw;
    g_kappa = mnl->kappa;
    boundary(mnl->kappa);
    init_sw_fields();
    sw_term((const su3**)hf.gaugefield, mnl->kappa, mnl->c_sw);
    mnl_backup_restore_globals(TM_RESTORE_GLOBALS);
  }

  /* Loop for measurements */
  for(j = 0; j < Nmeas; j++) {
    random_spinor_field_eo(mnl->pf, mnl->rngrepro, RN_GAUSS);
    monomial_list[0].derivativefunction(0, &hf); 
  }

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

  tm_stopwatch_pop(&g_timers, 0, 1, "");

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
    fprintf(stdout, "MG tuner for derivative monomials in Wilson twisted mass QCD\n");
    fprintf(stdout, "Version %s \n\n", TMLQCD_PACKAGE_VERSION);
    fprintf(stdout, "Please send bug reports to %s\n", TMLQCD_PACKAGE_BUGREPORT);
    fprintf(stdout, "Usage:   deriv_mg_tune [options]\n");
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


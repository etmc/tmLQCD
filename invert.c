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
 * invert for twisted mass QCD
 *
 * Author: Carsten Urbach
 *         urbach@physik.fu-berlin.de
 *
 *******************************************************************************/

#include"lime.h"
#ifdef HAVE_CONFIG_H
# include<config.h>
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
# include <omp.h>
#endif
#include "global.h"
#include "git_hash.h"
#include "getopt.h"
#include "linalg_eo.h"
#include "geometry_eo.h"
#include "start.h"
/*#include "eigenvalues.h"*/
#include "measure_gauge_action.h"
#ifdef TM_USE_MPI
#include "xchange/xchange.h"
#endif
#include <io/utils.h>
#include "source_generation.h"
#include "read_input.h"
#include "mpi_init.h"
#include "sighandler.h"
#include "boundary.h"
#include "solver/solver.h"
#include "init/init.h"
#include "init/init_gauge_tmp.h"
#include "smearing/stout.h"
#include "invert_eo.h"
#include "monomial/monomial.h"
#include "ranlxd.h"
#include "phmc.h"
#include "operator/D_psi.h"
#include "little_D.h"
#include "reweighting_factor.h"
#include "linalg/convert_eo_to_lexic.h"
#include "block.h"
#include "operator.h"
#include "sighandler.h"
#include "solver/generate_dfl_subspace.h"
#include "prepare_source.h"
#include <io/params.h>
#include <io/gauge.h>
#include <io/spinor.h>
#include <io/utils.h>
#include "solver/dirac_operator_eigenvectors.h"
#include "source_generation.h"
#include "P_M_eta.h"
#include "operator/tm_operators.h"
#include "operator/Dov_psi.h"
#include "solver/spectral_proj.h"
#ifdef TM_USE_QUDA
#  include "quda_interface.h"
#endif
#ifdef TM_USE_QPHIX
#  include "qphix_interface.h"
#endif
#ifdef DDalphaAMG
#  include "DDalphaAMG_interface.h"
#endif
#include "meas/measurements.h"
#include "source_generation.h"
#include "expo.h"

#define CONF_FILENAME_LENGTH 500

extern int nstore;
int check_geometry();

static void usage();
static void process_args(int argc, char *argv[], char ** input_filename, char ** filename);
static void set_default_filenames(char ** input_filename, char ** filename);
static void invert_compute_modenumber();

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
  struct stout_parameters params_smear;

#ifdef _KOJAK_INST
#pragma pomp inst init
#pragma pomp inst begin(main)
#endif

#if (defined SSE || defined SSE2 || SSE3)
  signal(SIGILL, &catch_ill_inst);
#endif

  DUM_DERI = 8;
  DUM_MATRIX = DUM_DERI + 5;
  NO_OF_SPINORFIELDS = DUM_MATRIX + 4;

  //4 extra fields (corresponding to DUM_MATRIX+0..5) for deg. and ND matrix mult.  
  NO_OF_SPINORFIELDS_32 = 6;

  verbose = 0;
  g_use_clover_flag = 0;


  process_args(argc,argv,&input_filename,&filename);
  set_default_filenames(&input_filename, &filename);

  init_parallel_and_read_input(argc, argv, input_filename);

  /* this DBW2 stuff is not needed for the inversion ! */
  if (g_dflgcr_flag == 1) {
    even_odd_flag = 0;
  }
  g_rgi_C1 = 0;
  if (Nsave == 0) {
    Nsave = 1;
  }

  if (g_running_phmc) {
    NO_OF_SPINORFIELDS = DUM_MATRIX + 8;
  }

  tmlqcd_mpi_init(argc, argv);

  g_dbw2rand = 0;

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
  j = init_gauge_field(VOLUMEPLUSRAND, 1);
  j += init_gauge_field_32(VOLUMEPLUSRAND, 1);
#else
  j = init_gauge_field(VOLUMEPLUSRAND, 0);
  j += init_gauge_field_32(VOLUMEPLUSRAND, 0);  
#endif
  if(restoresu3_flag) {
    j += init_gauge_tmp(VOLUMEPLUSRAND);
  }
 
  if (j != 0) {
    fprintf(stderr, "Not enough memory for gauge_fields! Aborting...\n");
    exit(-1);
  }
  j = init_geometry_indices(VOLUMEPLUSRAND);
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
    j += init_spinor_field_32(VOLUMEPLUSRAND / 2, NO_OF_SPINORFIELDS_32);   
  }
  else {
    j = init_spinor_field(VOLUMEPLUSRAND, NO_OF_SPINORFIELDS);
    j += init_spinor_field_32(VOLUMEPLUSRAND, NO_OF_SPINORFIELDS_32);   
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
  /* for mixed precision solvers, the 32 bit halfspinor field must always be there */
  j = init_dirac_halfspinor32();
  if (j != 0)
  {
    fprintf(stderr, "Not enough memory for 32-bit halffield! Aborting...\n");
    exit(-1);
  }
#  if (defined _PERSISTENT)
  if (even_odd_flag)
    init_xchange_halffield();
#  endif
#endif

  for (j = 0; j < Nmeas; j++) {
    int n_written = snprintf(conf_filename, CONF_FILENAME_LENGTH, "%s.%.4d", gauge_input_filename, nstore);
    if( n_written < 0 || n_written >= CONF_FILENAME_LENGTH ){
      char error_message[500];
      snprintf(error_message,
               500,
               "Encoding error or gauge configuration filename "
               "longer than %d characters! See invert.c CONF_FILENAME_LENGTH\n", 
               CONF_FILENAME_LENGTH);
      fatal_error(error_message, "invert.c");
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
    if (restoresu3_flag) {
      if (g_cart_id == 0) 
        printf("# Restoring SU(3) matrices.\n");
      for(int ix=0;ix<VOLUME;ix++) {
        for(int mu=0;mu<4;mu++){
          su3 *v, *w;
          v=&(g_gauge_field[ix][mu]);
          w=&(gauge_tmp[ix][mu]);
          _su3_assign(*w,*v);
          restoresu3_in_place(v);
        }
      }
    }

    if (g_cart_id == 0) {
      printf("# Finished reading gauge field.\n");
      fflush(stdout);
    }
#ifdef TM_USE_MPI
    xchange_gauge(g_gauge_field);
    if (restoresu3_flag) {
      xchange_gauge(gauge_tmp);
    }
#endif
    /*Convert to a 32 bit gauge field, after xchange*/
    convert_32_gauge_field(g_gauge_field_32, g_gauge_field, VOLUMEPLUSRAND);
    /*compute the energy of the gauge field*/
    plaquette_energy = measure_plaquette( (const su3**) g_gauge_field);

    if (g_cart_id == 0) {
      printf("# The computed plaquette value is %e.\n", plaquette_energy / (6.*VOLUME*g_nproc));
      fflush(stdout);
    }

    if (restoresu3_flag) {
      double plaquette_old = measure_plaquette( (const su3**) gauge_tmp);
      if (g_cart_id == 0) {
        printf("# The computed plaquette value before restoring SU(3) is %e\n which differ from the new one of %e.\n",
               plaquette_old / (6.*VOLUME*g_nproc), (plaquette_energy-plaquette_old) / (6.*VOLUME*g_nproc));
        fflush(stdout);
      }
    }

    if (use_stout_flag == 1){
      params_smear.rho = stout_rho;
      params_smear.iterations = stout_no_iter;
/*       if (stout_smear((su3_tuple*)(g_gauge_field[0]), &params_smear, (su3_tuple*)(g_gauge_field[0])) != 0) */
/*         exit(1) ; */
      g_update_gauge_copy = 1;
      plaquette_energy = measure_plaquette( (const su3**) g_gauge_field);

      if (g_cart_id == 0) {
        printf("# The plaquette value after stouting is %e\n", plaquette_energy / (6.*VOLUME*g_nproc));
        fflush(stdout);
      }
    }

    /* if any measurements are defined in the input file, do them here */
    measurement * meas;
    for(int imeas = 0; imeas < no_measurements; imeas++){
      meas = &measurement_list[imeas];
      if (g_proc_id == 0) {
        fprintf(stdout, "#\n# Beginning online measurement.\n");
      }
      meas->measurefunc(nstore, imeas, even_odd_flag);
    }

    if (reweighting_flag == 1) {
      reweighting_factor(reweighting_samples, nstore);
    }

    /* Compute minimal eigenvalues, if wanted */
    if (compute_evs != 0) {
      eigenvalues(&no_eigenvalues, 5000, eigenvalue_precision,
                  0, compute_evs, nstore, even_odd_flag);
    }
    if (phmc_compute_evs != 0) {
#ifdef TM_USE_MPI
      MPI_Finalize();
#endif
      return(0);
    }

    /* Compute the mode number or topological susceptibility using spectral projectors, if wanted*/
    if(compute_modenumber != 0 || compute_topsus !=0){
      invert_compute_modenumber(); 
    }

    //  set up blocks if Deflation is used 
    if (g_dflgcr_flag) 
      init_blocks(nblocks_t, nblocks_x, nblocks_y, nblocks_z);
    
    if(SourceInfo.type == SRC_TYPE_VOL || SourceInfo.type == SRC_TYPE_PION_TS || SourceInfo.type == SRC_TYPE_GEN_PION_TS) {
      index_start = 0;
      index_end = 1;
    }

    g_precWS=NULL;
    if(use_preconditioning == 1){
      /* todo load fftw wisdom */
#if (defined HAVE_FFTW ) && !( defined TM_USE_MPI)
      loadFFTWWisdom(g_spinor_field[0],g_spinor_field[1],T,LX);
#else
      use_preconditioning=0;
#endif
    }

    if (g_cart_id == 0) {
      fprintf(stdout, "#\n"); /*Indicate starting of the operator part*/
    }
    for(op_id = 0; op_id < no_operators; op_id++) {
      boundary(operator_list[op_id].kappa);
      g_kappa = operator_list[op_id].kappa; 
      g_mu = operator_list[op_id].mu;
      g_c_sw = operator_list[op_id].c_sw;
      // DFLGCR and DFLFGMRES
      if(operator_list[op_id].solver == DFLGCR || operator_list[op_id].solver == DFLFGMRES) {
        generate_dfl_subspace(g_N_s, VOLUME, reproduce_randomnumber_flag);
      }

      if(use_preconditioning==1 && PRECWSOPERATORSELECT[operator_list[op_id].solver]!=PRECWS_NO ){
        printf("# Using preconditioning with treelevel preconditioning operator: %s \n",
              precWSOpToString(PRECWSOPERATORSELECT[operator_list[op_id].solver]));
        /* initial preconditioning workspace */
        operator_list[op_id].precWS=(spinorPrecWS*)malloc(sizeof(spinorPrecWS));
        spinorPrecWS_Init(operator_list[op_id].precWS,
                  operator_list[op_id].kappa,
                  operator_list[op_id].mu/2./operator_list[op_id].kappa,
                  -(0.5/operator_list[op_id].kappa-4.),
                  PRECWSOPERATORSELECT[operator_list[op_id].solver]);
        g_precWS = operator_list[op_id].precWS;

        if(PRECWSOPERATORSELECT[operator_list[op_id].solver] == PRECWS_D_DAGGER_D) {
          fitPrecParams(op_id);
        }
      }

      for(isample = 0; isample < no_samples; isample++) {
        for (ix = index_start; ix < index_end; ix++) {
          if (g_cart_id == 0) {
            fprintf(stdout, "#\n"); /*Indicate starting of new index*/
          }
          /* we use g_spinor_field[0-7] for sources and props for the moment */
          /* 0-3 in case of 1 flavour  */
          /* 0-7 in case of 2 flavours */
          prepare_source(nstore, isample, ix, op_id, read_source_flag, source_location, random_seed);
          //randmize initial guess for eigcg if needed-----experimental
          if( (operator_list[op_id].solver == INCREIGCG) && (operator_list[op_id].solver_params.eigcg_rand_guess_opt) ){ //randomize the initial guess
              gaussian_volume_source( operator_list[op_id].prop0, operator_list[op_id].prop1,isample,ix,0); //need to check this
          } 
          operator_list[op_id].inverter(op_id, index_start, operator_list[op_id].write_prop_flag);
        }
      }


      if(use_preconditioning==1 && operator_list[op_id].precWS!=NULL ){
        /* free preconditioning workspace */
        spinorPrecWS_Free(operator_list[op_id].precWS);
        free(operator_list[op_id].precWS);
      }

      if(operator_list[op_id].type == OVERLAP){
        free_Dov_WS();
      }

    }
    nstore += Nsave;
  }

#ifdef TM_USE_OMP
  free_omp_accumulators();
#endif
  free_blocks();
  free_dfl_subspace();
  free_gauge_tmp();
  free_gauge_field();
  free_gauge_field_32();
  free_geometry_indices();
  free_spinor_field();
  free_spinor_field_32();  
  free_moment_field();
  free_chi_spinor_field();
  free(filename);
  free(input_filename);
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

static void usage()
{
  fprintf(stdout, "Inversion for EO preconditioned Wilson twisted mass QCD\n");
  fprintf(stdout, "Version %s \n\n", PACKAGE_VERSION);
  fprintf(stdout, "Please send bug reports to %s\n", PACKAGE_BUGREPORT);
  fprintf(stdout, "Usage:   invert [options]\n");
  fprintf(stdout, "Options: [-f input-filename]\n");
  fprintf(stdout, "         [-o output-filename]\n");
  fprintf(stdout, "         [-v] more verbosity\n");
  fprintf(stdout, "         [-h|-? this help]\n");
  fprintf(stdout, "         [-V] print version information and exit\n");
  exit(0);
}

static void process_args(int argc, char *argv[], char ** input_filename, char ** filename) {
  int c;
  while ((c = getopt(argc, argv, "h?vVf:o:")) != -1) {
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
          fprintf(stdout,"%s %s\n",PACKAGE_STRING,git_hash);
        }
        exit(0);
        break;
      case 'h':
      case '?':
      default:
        if( g_proc_id == 0 ) {
          usage();
        }
        break;
    }
  }
}

static void set_default_filenames(char ** input_filename, char ** filename) {
  if( *input_filename == NULL ) {
    *input_filename = calloc(13, sizeof(char));
    strcpy(*input_filename,"invert.input");
  }
  
  if( *filename == NULL ) {
    *filename = calloc(7, sizeof(char));
    strcpy(*filename,"output");
  } 
}

static void invert_compute_modenumber() {
  spinor * s_ = calloc(no_sources_z2*VOLUMEPLUSRAND+1, sizeof(spinor));
  spinor ** s  = calloc(no_sources_z2, sizeof(spinor*));
  if(s_ == NULL) { 
    printf("Not enough memory in %s: %d",__FILE__,__LINE__); exit(42); 
  }
  if(s == NULL) { 
    printf("Not enough memory in %s: %d",__FILE__,__LINE__); exit(42); 
  }
  for(int i = 0; i < no_sources_z2; i++) {
    s[i] = (spinor*)(((unsigned long int)(s_)+ALIGN_BASE)&~ALIGN_BASE)+i*VOLUMEPLUSRAND;
    random_spinor_field_lexic(s[i], reproduce_randomnumber_flag,RN_Z2);
	
    if(g_proc_id == 0) {
      printf("source %d \n", i);
    }
	
    if(compute_modenumber != 0){
      mode_number(s[i], mstarsq);
    }
	  
    if(compute_topsus !=0) {
      top_sus(s[i], mstarsq);
    }
  }
  free(s);
  free(s_);
}


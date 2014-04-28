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

#define MAIN_PROGRAM

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
#ifdef _USE_MPI
#include <mpi.h>
#endif
#ifdef OMP
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
#ifdef _USE_MPI
#include "xchange/xchange.h"
#endif
#include <io/utils.h>
#include "read_input.h"
#include "mpi_init.h"
#include "sighandler.h"
#include "boundary.h"
#include "solver/solver.h"
#include "init/init.h"
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
#include "solver/dfl_projector.h"
#include "solver/generate_dfl_subspace.h"
#include "prepare_source.h"
#include <io/params.h>
#include <io/gauge.h>
#include <io/spinor.h>
#include <io/utils.h>
#include "solver/dirac_operator_eigenvectors.h"
#include "P_M_eta.h"
#include "operator/tm_operators.h"
#include "operator/Dov_psi.h"
#include "solver/spectral_proj.h"

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

extern int nstore;
int check_geometry();

static void usage();
static void process_args(int argc, char *argv[], char ** input_filename, char ** filename);
static void set_default_filenames(char ** input_filename, char ** filename);

int main(int argc, char *argv[])
{
  FILE *parameterfile = NULL;
  int j, i, ix = 0, isample = 0, op_id = 0;
  char datafilename[206];
  char parameterfilename[206];
  char conf_filename[50];
  char * input_filename = NULL;
  char * filename = NULL;
  double plaquette_energy;
  struct stout_parameters params_smear;
  spinor **s, *s_;

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

#ifdef _USE_MPI

#  ifdef OMP
  int mpi_thread_provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_SERIALIZED, &mpi_thread_provided);
#  else
  MPI_Init(&argc, &argv);
#  endif

  MPI_Comm_rank(MPI_COMM_WORLD, &g_proc_id);
#else
  g_proc_id = 0;
#endif

  process_args(argc,argv,&input_filename,&filename);
  set_default_filenames(&input_filename, &filename);

  /* Read the input file */
  if( (j = read_input(input_filename)) != 0) {
    fprintf(stderr, "Could not find input file: %s\nAborting...\n", input_filename);
    exit(-1);
  }

#ifdef OMP
  init_openmp();
#endif

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
  start_ranlux(rlxd_level, random_seed);

  /* we need to make sure that we don't have even_odd_flag = 1 */
  /* if any of the operators doesn't use it                    */
  /* in this way even/odd can still be used by other operators */
  for(j = 0; j < no_operators; j++) if(!operator_list[j].even_odd_flag) even_odd_flag = 0;

#ifndef _USE_MPI
  g_dbw2rand = 0;
#endif

#ifdef _GAUGE_COPY
  j = init_gauge_field(VOLUMEPLUSRAND, 1);
#else
  j = init_gauge_field(VOLUMEPLUSRAND, 0);
#endif
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

  /* define the boundary conditions for the fermion fields */
  boundary(g_kappa);

  phmc_invmaxev = 1.;

  init_operators();


if(usegpu_flag){
  if(even_odd_flag){
      init_mixedsolve_eo(g_gauge_field);
  }
  else{
    init_mixedsolve(g_gauge_field);
  }
  #ifdef GPU_DOUBLE
    /*init double fields w/o momenta*/
    init_gpu_fields(0);
  #endif
  #ifdef TEMPORALGAUGE
    int retval;
    if((retval=init_temporalgauge(VOLUME, g_gauge_field)) !=0){
      if(g_proc_id == 0) printf("Error while initializing temporal gauge. Aborting...\n");   
      exit(200);
    }
  #endif
}//usegpu_flag
  
  
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
    sprintf(conf_filename, "%s.%.4d", gauge_input_filename, nstore);
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
#ifdef _USE_MPI
    xchange_gauge(g_gauge_field);
#endif

    /*compute the energy of the gauge field*/
    plaquette_energy = measure_plaquette( (const su3**) g_gauge_field);

    if (g_cart_id == 0) {
      printf("# The computed plaquette value is %e.\n", plaquette_energy / (6.*VOLUME*g_nproc));
      fflush(stdout);
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

    if (reweighting_flag == 1) {
      reweighting_factor(reweighting_samples, nstore);
    }

    /* Compute minimal eigenvalues, if wanted */
    if (compute_evs != 0) {
      eigenvalues(&no_eigenvalues, 5000, eigenvalue_precision,
                  0, compute_evs, nstore, even_odd_flag);
    }
    if (phmc_compute_evs != 0) {
#ifdef _USE_MPI
      MPI_Finalize();
#endif
      return(0);
    }

    /* Compute the mode number or topological susceptibility using spectral projectors, if wanted*/

    if(compute_modenumber != 0 || compute_topsus !=0){
      
      s_ = calloc(no_sources_z2*VOLUMEPLUSRAND+1, sizeof(spinor));
      s  = calloc(no_sources_z2, sizeof(spinor*));
      if(s_ == NULL) { 
	printf("Not enough memory in %s: %d",__FILE__,__LINE__); exit(42); 
      }
      if(s == NULL) { 
	printf("Not enough memory in %s: %d",__FILE__,__LINE__); exit(42); 
      }
      
      
      for(i = 0; i < no_sources_z2; i++) {
#if (defined SSE3 || defined SSE2 || defined SSE)
        s[i] = (spinor*)(((unsigned long int)(s_)+ALIGN_BASE)&~ALIGN_BASE)+i*VOLUMEPLUSRAND;
#else
        s[i] = s_+i*VOLUMEPLUSRAND;
#endif
	
        random_spinor_field_lexic(s[i], reproduce_randomnumber_flag,RN_Z2);
	
/* 	what is this here needed for?? */
/*         spinor *aux_,*aux; */
/* #if ( defined SSE || defined SSE2 || defined SSE3 ) */
/*         aux_=calloc(VOLUMEPLUSRAND+1, sizeof(spinor)); */
/*         aux = (spinor *)(((unsigned long int)(aux_)+ALIGN_BASE)&~ALIGN_BASE); */
/* #else */
/*         aux_=calloc(VOLUMEPLUSRAND, sizeof(spinor)); */
/*         aux = aux_; */
/* #endif */
	
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


    /* move to operators as well */
    if (g_dflgcr_flag == 1) {
      /* set up deflation blocks */
      init_blocks(nblocks_t, nblocks_x, nblocks_y, nblocks_z);

      /* the can stay here for now, but later we probably need */
      /* something like init_dfl_solver called somewhere else  */
      /* create set of approximate lowest eigenvectors ("global deflation subspace") */

      /*       g_mu = 0.; */
      /*       boundary(0.125); */
      generate_dfl_subspace(g_N_s, VOLUME, reproduce_randomnumber_flag);
      /*       boundary(g_kappa); */
      /*       g_mu = g_mu1; */

      /* Compute little Dirac operators */
      /*       alt_block_compute_little_D(); */
      if (g_debug_level > 0) {
        check_projectors(reproduce_randomnumber_flag);
        check_local_D(reproduce_randomnumber_flag);
      }
      if (g_debug_level > 1) {
        check_little_D_inversion(reproduce_randomnumber_flag);
      }

    }
    if(SourceInfo.type == 1) {
      index_start = 0;
      index_end = 1;
    }

    g_precWS=NULL;
    if(use_preconditioning == 1){
      /* todo load fftw wisdom */
#if (defined HAVE_FFTW ) && !( defined _USE_MPI)
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
      g_mu = 0.;

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
          prepare_source(nstore, isample, ix, op_id, read_source_flag, source_location);
          operator_list[op_id].inverter(op_id, index_start, 1);
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

#ifdef OMP
  free_omp_accumulators();
#endif


if(usegpu_flag){
  finalize_mixedsolve();
  #ifdef GPU_DOUBLE
    finalize_gpu_fields();
  #endif
    #ifdef TEMPORALGAUGE
      finalize_temporalgauge();
    #endif
}



  free_blocks();
  free_dfl_subspace();
  free_gauge_field();
  free_geometry_indices();
  free_spinor_field();
  free_moment_field();
  free_chi_spinor_field();
  free(filename);
  free(input_filename);
#ifdef _USE_MPI
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


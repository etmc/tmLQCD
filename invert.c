/* $Id$ */
/*******************************************************************************
*
* Hybrid-Monte-Carlo for twisted mass QCD
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
#ifdef MPI
#include <mpi.h>
#endif
#include "global.h"
#include "getopt.h"
#include "linalg_eo.h"
#include "geometry_eo.h"
#include "start.h"
/*#include "eigenvalues.h"*/
#include "observables.h"
#ifdef MPI
#include "xchange.h"
#endif
#include "io.h"
#include "io_utils.h"
#include "propagator_io.h"
#include "gauge_io.h"
#include "read_input.h"
#include "mpi_init.h"
#include "sighandler.h"
#include "boundary.h"
#include "solver/solver.h"
#include "init_gauge_field.h"
#include "init_geometry_indices.h"
#include "init_spinor_field.h"
#include "init_moment_field.h"
#include "init_dirac_halfspinor.h"
#include "xchange_halffield.h"
#include "stout_smear.h"
#include "invert_eo.h"
#include "phmc.h"
#include "D_psi.h"
#include "little_D.h"
#include "linalg/convert_eo_to_lexic.h"
#include "block.h"
#include "solver/dfl_projector.h"
#include "solver/generate_dfl_subspace.h"

void usage(){
  fprintf(stdout, "Inversion for EO preconditioned Wilson twisted mass QCD\n");
  fprintf(stdout, "Version %s \n\n", PACKAGE_VERSION);
  fprintf(stdout, "Please send bug reports to %s\n", PACKAGE_BUGREPORT);
  fprintf(stdout, "Usage:   invert [options]\n");
  fprintf(stdout, "Options: [-f input-filename]\n");
  fprintf(stdout, "         [-o output-filename]\n");
  fprintf(stdout, "         [-h|-? this help]\n");
  exit(0);
}

extern int nstore;
int check_geometry();

int main(int argc,char *argv[]) {

  FILE *parameterfile=NULL, *ifs=NULL;
  int c, iter, j, ix=0, is=0, ic=0, err=0;
  char * filename = NULL;
  char datafilename[50];
  char parameterfilename[50];
  char conf_filename[50];
  char * input_filename = NULL;
  double plaquette_energy;
  double ratime, retime;
  
  int nsiter;
  complex * a1, * a2;

#ifdef _GAUGE_COPY
  int kb=0;
#endif
  double nrm1, nrm2;
  double atime=0., etime=0.;
#ifdef _KOJAK_INST
#pragma pomp inst init
#pragma pomp inst begin(main)
#endif

  DUM_DERI = 6;
  /* DUM_DERI + 2 is enough (not 7) */
  DUM_SOLVER = DUM_DERI+2;
  DUM_MATRIX = DUM_SOLVER+6;
  /* DUM_MATRIX + 2 is enough (not 6) */
  NO_OF_SPINORFIELDS = DUM_MATRIX+2;

  verbose = 0;
  g_use_clover_flag = 0;
  g_dflgcr_flag = 1;

#ifdef MPI
  MPI_Init(&argc, &argv);
#endif

  while ((c = getopt(argc, argv, "h?f:o:")) != -1) {
    switch (c) {
    case 'f': 
      input_filename = calloc(200, sizeof(char));
      strcpy(input_filename,optarg);
      break;
    case 'o':
      filename = calloc(200, sizeof(char));
      strcpy(filename,optarg);
      break;
    case 'h':
    case '?':
    default:
      usage();
      break;
    }
  }
  if(input_filename == NULL){
    input_filename = "hmc.input";
  }
  if(filename == NULL){
    filename = "output";
  } 

  /* Read the input file */
  read_input(input_filename);
  /* this DBW2 stuff is not needed for the inversion ! */
  if(g_dflgcr_flag == 1) {
    even_odd_flag = 0;
  }
  g_rgi_C1 = 0;
  if(Nsave == 0){
    Nsave = 1;
  }

  mpi_init(argc, argv);

  g_dbw2rand = 0;

#ifndef MPI
  g_dbw2rand = 0;
#endif

#ifdef _GAUGE_COPY
  j = init_gauge_field(VOLUMEPLUSRAND, 1);
#else
  j = init_gauge_field(VOLUMEPLUSRAND, 0);
#endif
  if ( j!= 0) {
    fprintf(stderr, "Not enough memory for gauge_fields! Aborting...\n");
    exit(-1);
  }
  j = init_geometry_indices(VOLUMEPLUSRAND);
  if ( j!= 0) {
    fprintf(stderr, "Not enough memory for geometry indices! Aborting...\n");
    exit(-1);
  }
  if(even_odd_flag) {
    j = init_spinor_field(VOLUMEPLUSRAND/2, NO_OF_SPINORFIELDS);
  }
  else {
    j = init_spinor_field(VOLUMEPLUSRAND, NO_OF_SPINORFIELDS);
  }
  if ( j!= 0) {
    fprintf(stderr, "Not enough memory for spinor fields! Aborting...\n");
    exit(-1);
  }

  g_mu = g_mu1; 
  if(g_proc_id == 0){    
    /*construct the filenames for the observables and the parameters*/
    strcpy(datafilename,filename);  strcat(datafilename,".data");
    strcpy(parameterfilename,filename);  strcat(parameterfilename,".para");
    
    parameterfile=fopen(parameterfilename, "w");
    write_first_messages(parameterfile, 0, 1);
  }

  /* define the geometry */
  geometry();

  /* define the boundary conditions for the fermion fields */
  boundary(g_kappa);

#ifdef _USE_HALFSPINOR
  j = init_dirac_halfspinor();
  if ( j!= 0) {
    fprintf(stderr, "Not enough memory for halffield! Aborting...\n");
    exit(-1);
  }
  if(g_sloppy_precision_flag == 1) {
    j = init_dirac_halfspinor32();
    if ( j!= 0) {
      fprintf(stderr, "Not enough memory for 32-Bit halffield! Aborting...\n");
      exit(-1);
    }
  }
#  if (defined _PERSISTENT)
  init_xchange_halffield();
#  endif
#endif

  for(j=0;j<Nmeas; j++) {
    sprintf(conf_filename,"%s.%.4d", gauge_input_filename, nstore);
    if (g_proc_id == 0){
      printf("Reading Gauge field from file %s\n", conf_filename); fflush(stdout);
    }
    read_lime_gauge_field(conf_filename);
    if (g_proc_id == 0){
      printf("done!\n"); fflush(stdout);
    }
#ifdef MPI
    xchange_gauge();
#endif

    /*compute the energy of the gauge field*/
    plaquette_energy = measure_gauge_action();

    if(g_proc_id == 0) {
      printf("The plaquette value is %e\n", plaquette_energy/(6.*VOLUME*g_nproc)); fflush(stdout);
    }

    if(use_stout_flag == 1) {
      if( stout_smear_gauge_field(stout_rho , stout_no_iter) != 0 ) exit(1) ;

      plaquette_energy = measure_gauge_action();

      if(g_proc_id == 0) {
	printf("The plaquette value after stouting is %e\n", plaquette_energy/(6.*VOLUME*g_nproc)); fflush(stdout);
      }
    }

    /* Compute minimal eigenvalues, if wanted */
    if(compute_evs != 0) {
      eigenvalues(&no_eigenvalues, max_solver_iterations, eigenvalue_precision, 
		  0, compute_evs, nstore, even_odd_flag);
    }
    if(phmc_compute_evs != 0) {
#ifdef MPI
      MPI_Finalize();
#endif
      return(0);
    }

    if(g_dflgcr_flag == 1) {
      /* set up deflation blocks */
      init_blocks();

      /* the can stay here for now, but later we probably need */
      /* something like init_dfl_solver called somewhere else  */
      /* create set of approximate lowest eigenvectors ("global deflation subspace") */

      init_dfl_subspace(g_N_s);
      generate_dfl_subspace(g_N_s, VOLUME);

      for (nsiter = 0; nsiter < g_N_s; nsiter++) { 
        /* add it to the basis */
        split_global_field(block_list[0].basis[nsiter], block_list[1].basis[nsiter], dfl_fields[nsiter]);
      }

      /* perform local orthonormalization */
      block_orthonormalize(block_list);
      block_orthonormalize(block_list+1);

      /* Compute little Dirac operators */
      block_compute_little_D_diagonal(block_list);
      block_compute_little_D_diagonal(block_list + 1);
      block_compute_little_D_offdiagonal();

      check_projectors();
      check_little_D_inversion();


      a1 = calloc(2*9*g_N_s, sizeof(complex));
      a2 = calloc(2*9*g_N_s, sizeof(complex));
      a1[0].re = 1.;
      lgcr(a2, a1, 10, 5, 1.e-15, 1, 2*g_N_s, 2*9*g_N_s, &little_D);
      free(a1);
      free(a2);
      /* TODO Generate projectors */
    }

    for(ix = index_start; ix < index_end; ix++) {
      is = (ix / 3);
      ic = (ix % 3);
      if(read_source_flag == 0) {
        if(source_location == 0)
          source_spinor_field(g_spinor_field[0], g_spinor_field[1], is, ic);
        else
          source_spinor_field_point_from_file(g_spinor_field[0], g_spinor_field[1], is, ic, source_location);
      }
      else {
#ifdef MPI
	ratime = MPI_Wtime();
#else
	ratime = (double)clock()/(double)(CLOCKS_PER_SEC);
#endif
	if(source_splitted) {
	  sprintf(conf_filename,"%s.%.2d", source_input_filename, ix);
	  if(g_proc_id == 0) {
	    printf("Reading source from %s\n", conf_filename);
	  }
	  read_source(g_spinor_field[0], g_spinor_field[1], conf_filename, source_format_flag, 0);
	}
	else {
	  sprintf(conf_filename,"%s", source_input_filename);
	  if(g_proc_id == 0) {
	    printf("Reading source no %d from %s\n", ix, conf_filename);
	  }
	  read_source(g_spinor_field[0], g_spinor_field[1], conf_filename, source_format_flag, ix);
	}
#ifdef MPI
	retime = MPI_Wtime();
#else
	retime = (double)clock()/(double)(CLOCKS_PER_SEC);
#endif
	if(g_proc_id == 0) {
	  printf("time for reading source was %e seconds\n", retime-ratime);
	}
      }
      if(g_proc_id == 0) {printf("mu = %e\n", g_mu);}

      if(propagator_splitted) {
	sprintf(conf_filename,"%s.%.2d.inverted", source_input_filename, ix);
      }
      else {
	sprintf(conf_filename,"%s.inverted", source_input_filename);
      }

      /* If the solver is _not_ CG we might read in */
      /* here some better guess                     */
      /* This also works for re-iteration           */
      if(solver_flag != CG && solver_flag != PCG) {
	ifs = fopen(conf_filename, "r");
	if(ifs != NULL) {
	  if(g_proc_id == g_stdio_proc){
	    printf("# Trying to read guess from file %s\n", conf_filename);
	    fflush(stdout);
	  }
	  fclose(ifs);
	  err = 0;
	  iter = get_propagator_type(conf_filename);
	  if(iter > -1 ) {
	    if(iter == 1) {
	      if(propagator_splitted){
		err = read_lime_spinor(g_spinor_field[2], g_spinor_field[3], conf_filename, 0);
	      }
	      else {
		err = read_lime_spinor(g_spinor_field[2], g_spinor_field[3], conf_filename, ix);
	      }
	    }
	    else if(iter == 0 ) {
	      if(propagator_splitted){
		err = read_lime_spinor(g_spinor_field[2], g_spinor_field[3], conf_filename, 0);
	      }
	      else {
		err = read_lime_spinor(g_spinor_field[2], g_spinor_field[3], conf_filename, 2*ix);
	      }
	    }
	    if(g_kappa != 0.) {
	      mul_r(g_spinor_field[3], 1./(2*g_kappa), g_spinor_field[3], VOLUME/2);
	      mul_r(g_spinor_field[2], 1./(2*g_kappa), g_spinor_field[2], VOLUME/2);
	    }
	  }
	  else {
	    /* trying cmi format */
	    if(g_proc_id == g_stdio_proc){
	      printf("# Trying cmi format instead of ETMC standard\n");
	      fflush(stdout);
	    }
	    err = read_source(g_spinor_field[2], g_spinor_field[3], conf_filename, 1, 0);
	  }
	  if(err != 0) {
	    zero_spinor_field(g_spinor_field[3],VOLUME/2);
	  }
	}
	else {
	  zero_spinor_field(g_spinor_field[3],VOLUME/2);
	}
      }
      else {
	zero_spinor_field(g_spinor_field[3],VOLUME/2);
      }

#ifdef MPI
      atime = MPI_Wtime();
#else
      atime = (double)clock()/(double)(CLOCKS_PER_SEC);
#endif
       iter = invert_eo(g_spinor_field[2], g_spinor_field[3], g_spinor_field[0], g_spinor_field[1],  
 		       solver_precision, max_solver_iterations, solver_flag,g_relative_precision_flag, 
 		       sub_evs_cg_flag, even_odd_flag); 
#ifdef MPI
      etime = MPI_Wtime();
#else
      etime = (double)clock()/(double)(CLOCKS_PER_SEC);
#endif

      /* To write in standard format */
      /* we have to mult. by 2*kappa */
      if(write_prop_format_flag != 11 && g_kappa != 0.) {
	mul_r(g_spinor_field[2], (2*g_kappa), g_spinor_field[2], VOLUME/2);
	mul_r(g_spinor_field[3], (2*g_kappa), g_spinor_field[3], VOLUME/2);
      }

      if(write_prop_format_flag < 10) {
	if(propagator_splitted) {
	  write_propagator_type(write_prop_format_flag, conf_filename);
	  write_xlf_info(plaquette_energy/(6.*VOLUME*g_nproc), nstore, conf_filename, 1);
	  /* write the source depending on format */
	  if(write_prop_format_flag == 1) {
	    write_source(g_spinor_field[0], g_spinor_field[1], conf_filename, 1, 32);
	  }
	}
	else {
	  if(ix == index_start) {
	    write_propagator_type(write_prop_format_flag, conf_filename);
	  }
	  write_xlf_info(plaquette_energy/(6.*VOLUME*g_nproc), nstore, conf_filename, 1);
	  /* write the source depending on format */
	  if(write_prop_format_flag == 1) {
	    write_source(g_spinor_field[0], g_spinor_field[1], conf_filename, 1, 32);
	  }
	}
      }
#ifdef MPI
      ratime = MPI_Wtime();
#else
      ratime = (double)clock()/(double)(CLOCKS_PER_SEC);
#endif
      write_propagator(g_spinor_field[2], g_spinor_field[3], conf_filename, 1, 
		       prop_precision_flag, write_prop_format_flag);
      
#ifdef MPI
      retime = MPI_Wtime();
#else
      retime = (double)clock()/(double)(CLOCKS_PER_SEC);
#endif
      if(g_proc_id == 0) {
	printf("time for writing prop was %e seconds\n", retime-ratime);
      }
      /* Check the result */
      M_full(g_spinor_field[4], g_spinor_field[5], g_spinor_field[2], g_spinor_field[3]); 

      if(write_prop_format_flag != 11) {
	if(g_kappa != 0.) {
	  mul_r(g_spinor_field[4], 1./(2*g_kappa), g_spinor_field[4], VOLUME/2);  
	  mul_r(g_spinor_field[5], 1./(2*g_kappa), g_spinor_field[5], VOLUME/2); 
	}
      }

      diff(g_spinor_field[4], g_spinor_field[4], g_spinor_field[0], VOLUME/2); 
      diff(g_spinor_field[5], g_spinor_field[5], g_spinor_field[1], VOLUME/2); 

      nrm1 = square_norm(g_spinor_field[4], VOLUME/2); 
      nrm2 = square_norm(g_spinor_field[5], VOLUME/2); 

      if(g_proc_id == 0) {
	printf("Inversion for source %d done in %d iterations, residue = %e!\n", ix, iter, nrm1+nrm2);
	printf("Inversion done in %1.2e sec. \n", etime-atime);
        write_inverter_info(nrm1+nrm2, iter, 0, 1, conf_filename);
      }
    }
    nstore+=Nsave;
  }
#ifdef MPI
  MPI_Finalize();
#endif
  free_blocks();
  free_dfl_subspace();
  free_gauge_field();
  free_geometry_indices();
  free_spinor_field();
  free_moment_field();
  return(0);
#ifdef _KOJAK_INST
#pragma pomp inst end(main)
#endif
}

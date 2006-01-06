/* $Id$ */
/*******************************************************************************
*
* File hybrid.c
*
* Hybrid-Monte-Carlo for twisted mass QCD
*
* Author: Carsten Urbach
*         urbach@physik.fu-berlin.de
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <signal.h>
#ifdef MPI
#include <mpi.h>
#endif
#include "getopt.h"
#include "global.h"
#include "linalg_eo.h"
#include "geometry_eo.h"
#include "start.h"
#include "observables.h"
#ifdef MPI
#include "xchange.h"
#endif
#include "io.h"
#include "read_input.h"
#include "mpi_init.h"
#include "sighandler.h"
#include "update_tm.h"
#include "boundary.h"
#include "solver/solver.h"
#include "init_gauge_field.h"
#include "init_geometry_indices.h"
#include "init_spinor_field.h"
#include "init_moment_field.h"
#include "invert_eo.h"


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

  FILE *parameterfile=NULL, *countfile=NULL, *ifs=NULL;
  int c, iter, j, ix=0, is=0, ic=0, write_counter=0;
  char * filename = NULL;
  char datafilename[50];
  char parameterfilename[50];
  char conf_filename[50];
  char * input_filename = NULL;
  char * nstore_filename = ".nstore_counter";
#ifdef _GAUGE_COPY
  int kb=0;
#endif
  double nrm1, nrm2;
#ifdef MPI
  double atime=0., etime=0.;
#endif

  verbose = 0;
  g_use_clover_flag = 0;
  g_nr_of_psf = 1;

#ifndef XLC 
  signal(SIGUSR1,&catch_del_sig);
  signal(SIGUSR2,&catch_del_sig);
  signal(SIGTERM,&catch_del_sig);
  signal(SIGXCPU,&catch_del_sig);
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
  g_rgi_C1 = 0;

  mpi_init(argc, argv);

  if(nstore == -1) {
    countfile = fopen(nstore_filename, "r");
    if(countfile != NULL) {
      fscanf(countfile, "%d\n", &nstore);
      fclose(countfile);
    }
    else {
      nstore = 0;
    }
    write_counter=1;
  }

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
    exit(0);
  }
  j = init_geometry_indices(VOLUMEPLUSRAND);
  if ( j!= 0) {
    fprintf(stderr, "Not enough memory for geometry indices! Aborting...\n");
    exit(0);
  }
  j = init_spinor_field(VOLUMEPLUSRAND/2, NO_OF_SPINORFIELDS);
  if ( j!= 0) {
    fprintf(stderr, "Not enough memory for spinor fields! Aborting...\n");
    exit(0);
  }

  g_mu = g_mu1; 
  if(g_proc_id == 0){    
    
/*     fscanf(fp6,"%s",filename); */
    /*construct the filenames for the observables and the parameters*/
    strcpy(datafilename,filename);  strcat(datafilename,".data");
    strcpy(parameterfilename,filename);  strcat(parameterfilename,".para");
    
    parameterfile=fopen(parameterfilename, "w");
    printf("# This is the hmc code for twisted Mass Wilson QCD\n\nVersion %s\n", PACKAGE_VERSION);
#ifdef SSE
    printf("# The code was compiled with SSE instructions\n");
#endif
#ifdef SSE2
    printf("# The code was compiled with SSE2 instructions\n");
#endif
#ifdef P4
    printf("# The code was compiled for Pentium4\n");
#endif
#ifdef OPTERON
    printf("# The code was compiled for AMD Opteron\n");
#endif
#ifdef _NEW_GEOMETRY
    printf("# The code was compiled with -D_NEW_GEOMETRY\n");
#endif
#ifdef _GAUGE_COPY
    printf("# The code was compiled with -D_GAUGE_COPY\n");
#endif
    printf("# The lattice size is %d x %d x %d x %d\n",
	   (int)(T*g_nproc_t), (int)(LX*g_nproc_x), (int)(LY), (int)(LZ));
    printf("# The local lattice size is %d x %d x %d x %d\n", 
	   (int)(T), (int)(LX), (int)(LY),(int) LZ);
    printf("# beta = %f , kappa= %f, mu= %f \n",g_beta,g_kappa,g_mu);
    printf("# mus = %f, %f, %f\n", g_mu1, g_mu2, g_mu3);
    
    fprintf(parameterfile, "The lattice size is %d x %d^3\n", (int)(g_nproc_t*T),(int)(L));
    fprintf(parameterfile, "The local lattice size is %d x %d^3\n", (int)(T),(int)(L));
    fprintf(parameterfile, "g_beta = %f , g_kappa= %f, g_kappa*csw/8= %f g_mu = %f \n",g_beta,g_kappa,g_ka_csw_8, g_mu);
    fprintf(parameterfile, "boundary of fermion fields (t,x,y,z): %f %f %f %f \n",X0,X1,X2,X3);
    fprintf(parameterfile, "ITER_MAX_BCG=%d, EPS_SQ0=%e, EPS_SQ1=%e EPS_SQ2=%e, EPS_SQ3=%e \n"
	    ,ITER_MAX_BCG,EPS_SQ0,EPS_SQ1,EPS_SQ2,EPS_SQ3);
    fprintf(parameterfile,"dtau=%f, Nsteps=%d, Nmeas=%d, Nskip=%d, integtyp=%d, nsmall=%d \n",
	    dtau,Nsteps,Nmeas,Nskip,integtyp,nsmall);
    fprintf(parameterfile,"mu = %f, mu2=%f\n ", g_mu, g_mu2);
    
  }

  /* define the geometry */
  geometry();

  /* define the boundary conditions for the fermion fields */
  boundary();

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
#ifdef _GAUGE_COPY
    /* set the backward gauge field */
    for(ix = 0; ix < VOLUME;ix++) {
      kb=g_idn[ix][0];
      _su3_assign(g_gauge_field_back[ix][0],g_gauge_field[kb][0]);
      kb=g_idn[ix][1];
      _su3_assign(g_gauge_field_back[ix][1],g_gauge_field[kb][1]);
      kb=g_idn[ix][2];
      _su3_assign(g_gauge_field_back[ix][2],g_gauge_field[kb][2]);
      kb=g_idn[ix][3];
      _su3_assign(g_gauge_field_back[ix][3],g_gauge_field[kb][3]);
    }
#endif  
    if(source_format_flag == 0) {
      sprintf(conf_filename,"%s%.2d.%.4d", "prop.mass", mass_number, nstore);
    }
    else if(source_format_flag == 1) {
      sprintf(conf_filename,"%s.inverted", source_input_filename);
    }
    
    for(ix = index_start; ix < index_end; ix++) {
      is = (ix / 3);
      ic = (ix % 3);
      if(read_source_flag == 0) {
	source_spinor_field(g_spinor_field[0], g_spinor_field[1], is, ic);
      }
      else {
	if(source_format_flag == 0) {
	  sprintf(conf_filename,"%s%.2d.is%.1dic%.1d.%.4d", source_input_filename, mass_number, is, ic, nstore); 
	  read_spinorfield_eo_time(g_spinor_field[0], g_spinor_field[1], conf_filename); 
	}
	else if(source_format_flag == 1) {
	  printf("Reading source from %s\n",source_input_filename);
	  read_spinorfield_cm_single(g_spinor_field[0], g_spinor_field[1], source_input_filename, source_time_slice, 1);
	}
      }

      if(g_proc_id == 0) {printf("mu = %e\n", g_mu);}

      if(source_format_flag == 0) {
	sprintf(conf_filename,"%s%.2d.is%.1dic%.1d.%.4d", "prop.mass", mass_number, is, ic, nstore);
      }
      else if(source_format_flag == 1) {
	sprintf(conf_filename, "%s.inverted", source_input_filename);
      }

      /* If the solver is _not_ CG we might read in */
      /* here some better guess                     */
      /* This also works for re-iteration           */
      if(solver_flag != CG) {
	ifs = fopen(conf_filename, "r");
	if(ifs != NULL) {
	  if(g_proc_id == g_stdio_proc){
	    printf("Reading in from file %s\n", filename);
	    fflush(stdout);
	  }
	  fclose(ifs);
	  
	  read_spinorfield_eo_time(g_spinor_field[2], g_spinor_field[3], conf_filename);
	  mul_r(g_spinor_field[3], 1./(2*g_kappa), g_spinor_field[3], VOLUME/2);
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
#endif
      iter = invert_eo(g_spinor_field[2], g_spinor_field[3], g_spinor_field[0], g_spinor_field[1], 
		       solver_precision, max_solver_iterations, solver_flag,g_relative_precision_flag);
#ifdef MPI
      etime = MPI_Wtime();
#endif

      if(source_format_flag == 0) {
	/* To write in standard format */
	/* we have to mult. by 2*kappa */
	mul_r(g_spinor_field[2], (2*g_kappa), g_spinor_field[2], VOLUME/2);  
	mul_r(g_spinor_field[3], (2*g_kappa), g_spinor_field[3], VOLUME/2);
	write_spinorfield_eo_time_p(g_spinor_field[2], g_spinor_field[3], conf_filename, 0);
      }
      else if(source_format_flag == 1) {
	write_spinorfield_cm_single(g_spinor_field[2], g_spinor_field[3], conf_filename);
      }

      /* Check the result */
      M_full(g_spinor_field[4], g_spinor_field[5], g_spinor_field[2], g_spinor_field[3]); 
      if(source_format_flag == 0) {
	mul_r(g_spinor_field[4], 1./(2*g_kappa), g_spinor_field[4], VOLUME/2);  
	mul_r(g_spinor_field[5], 1./(2*g_kappa), g_spinor_field[5], VOLUME/2); 
      }
      diff(g_spinor_field[4], g_spinor_field[4], g_spinor_field[0], VOLUME/2); 
      diff(g_spinor_field[5], g_spinor_field[5], g_spinor_field[1], VOLUME/2); 

      nrm1 = square_norm(g_spinor_field[4], VOLUME/2); 
      nrm2 = square_norm(g_spinor_field[5], VOLUME/2); 

      if(g_proc_id == 0) {
	printf("Inversion for is = %d, ic = %d done in %d iterations, residue = %e!\n", is, ic, iter, nrm1+nrm2);
#ifdef MPI
	printf("Inversion done in %e sec. (MPI_Wtime)\n", etime-atime);
#endif
      }
    }
    nstore++;
    if(g_proc_id == 0) {
      countfile = fopen(nstore_filename, "w");
      fprintf(countfile, "%d\n", nstore);
      fclose(countfile);
    }
  }

#ifdef MPI
  MPI_Finalize();
#endif
  free_gauge_field();
  free_geometry_indices();
  free_spinor_field();
  free_moment_field();
  return(0);
}

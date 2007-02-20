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
#include "observables.h"
#ifdef MPI
#include "xchange.h"
#endif
#include "io.h"
#include "read_input.h"
#include "mpi_init.h"
#include "sighandler.h"
#include "boundary.h"
#include "solver/solver.h"
#include "init_gauge_field.h"
#include "init_geometry_indices.h"
#include "init_spinor_field.h"
#include "init_moment_field.h"
#include "update_backward_gauge.h"
#include "tm_operators.h"
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

  FILE *parameterfile=NULL, *ifs=NULL;
  int c, iter, j, ix=0, is=0, ic=0;
  char * filename = NULL;
  char datafilename[50];
  char parameterfilename[50];
  char conf_filename[50];
  char * input_filename = NULL;
  double plaquette_energy;
#ifdef _GAUGE_COPY
  int kb=0;
#endif
  double nrm1, nrm2;
#ifdef MPI
  double atime=0., etime=0.;
#endif

  DUM_DERI = 6;
  /* DUM_DERI + 2 is enough (not 7) */
  DUM_SOLVER = DUM_DERI+2;
  DUM_MATRIX = DUM_SOLVER+0;
  /* DUM_MATRIX + 2 is enough (not 6) */
  NO_OF_SPINORFIELDS = DUM_MATRIX+2;

  verbose = 0;
  g_use_clover_flag = 0;
  g_nr_of_psf = 1;

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
  g_rgi_C1 = 0;
  if(Nskip == 0){
    Nskip = 1;
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
    write_first_messages(parameterfile, 0, 1);
  }

  /* define the geometry */
  geometry();

  /* define the boundary conditions for the fermion fields */
  boundary();


  sprintf(conf_filename,"%s.%.4d", gauge_input_filename, nstore);
  if (g_proc_id == 0) {
    printf("Reading Gauge field from file %s\n", conf_filename); fflush(stdout);
  }
  if(gauge_precision_read_flag == 64) {
    read_lime_gauge_field(conf_filename);
  }
  else if(gauge_precision_read_flag == 32){
      read_lime_gauge_field_singleprec(conf_filename);
  }
  if (g_proc_id == 0){
      printf("done!\n"); fflush(stdout);
  }
#ifdef MPI
  xchange_gauge();
#endif
#ifdef _GAUGE_COPY
  update_backward_gauge();
#endif
    
    /*compute the energy of the gauge field*/
  plaquette_energy = measure_gauge_action();

  if(g_proc_id == 0) {
      printf("The plaquette value is %e\n", plaquette_energy/(6.*VOLUME*g_nproc)); fflush(stdout);
  }

  if(source_format_flag == 0) {
    sprintf(conf_filename,"%s", source_input_filename); 
    if(g_proc_id == 0) {
      printf("Reading source from %s\n", conf_filename);
    }
    read_spinorfield_eo_time(g_spinor_field[0], g_spinor_field[1], conf_filename); 
  }
  else if(source_format_flag == 1) {
    sprintf(conf_filename,"%s", source_input_filename);
    if(g_proc_id == 0) {
      printf("Reading source from %s\n", conf_filename);
    }
    read_spinorfield_cm_single(g_spinor_field[0], g_spinor_field[1], conf_filename, -1, 1);
  }
  if(g_proc_id == 0) {printf("mu = %e\n", g_mu);}
  
  if(source_format_flag == 0) {
      sprintf(conf_filename,"%s.applied", source_input_filename);
  }
  else if(source_format_flag == 1) {
      sprintf(conf_filename, "%s.applied", source_input_filename);
  }
  
#ifdef MPI
  atime = MPI_Wtime();
#endif
  /* Now apply (1+i mu gamma_5)^-1 (M-1) Nmeas times */
  for(j=0;j<Nmeas/2; j++) {
      M_minus_1_timesC(g_spinor_field[4], g_spinor_field[5],
		       g_spinor_field[0], g_spinor_field[1]);
      M_minus_1_timesC(g_spinor_field[0], g_spinor_field[1],
		       g_spinor_field[4], g_spinor_field[5]);
  }
  
#ifdef MPI
  etime = MPI_Wtime();
#endif
  
  if(source_format_flag == 0) {
      /* To write in standard format */
      /* we have to mult. by 2*kappa */
      mul_r(g_spinor_field[2], (2*g_kappa), g_spinor_field[0], VOLUME/2);  
      mul_r(g_spinor_field[3], (2*g_kappa), g_spinor_field[1], VOLUME/2);
      write_spinorfield_eo_time_p(g_spinor_field[2], g_spinor_field[3], conf_filename, 0);
  }
  else if(source_format_flag == 1) {
      write_spinorfield_cm_single(g_spinor_field[0], g_spinor_field[1], conf_filename);
  }
  
  if(g_proc_id == 0) {
#ifdef MPI
      printf("Done in %e sec. (MPI_Wtime)\n", etime-atime);
#endif
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

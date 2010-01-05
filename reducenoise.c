/***********************************************************************
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
 ***********************************************************************/
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

#include <lime.h>
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
#include <io/utils.h>
#include <io/gauge.h>
#include <io/spinor.h>
#include "read_input.h"
#include "mpi_init.h"
#include "sighandler.h"
#include "boundary.h"
#include "solver/solver.h"
#include "init_gauge_field.h"
#include "init_geometry_indices.h"
#include "init_spinor_field.h"
#include "init_dirac_halfspinor.h"
#include "init_moment_field.h"
#include "update_backward_gauge.h"
#include "tm_operators.h"
#include "stout_smear.h"
#include "invert_eo.h"
#include "gamma.h"

void usage(){
  fprintf(stdout, "Application of four times the Hopping Matrix to a source\n");
  fprintf(stdout, "for noise reduction\n");
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
  WRITER *writer = NULL;
  paramsSourceFormat *sourceFormat = NULL;
  char * xlfmessage = NULL;
  char * gaugelfn = NULL;
  DML_Checksum gaugecksum;


  DUM_DERI = 6;
  /* DUM_DERI + 2 is enough (not 7) */
  DUM_SOLVER = DUM_DERI+2;
  DUM_MATRIX = DUM_SOLVER+0;
  /* DUM_MATRIX + 2 is enough (not 6) */
  NO_OF_SPINORFIELDS = DUM_MATRIX+2;

  verbose = 0;
  g_use_clover_flag = 0;

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
  if(Nsave == 0){
    Nsave = 1;
  }
  tmlqcd_mpi_init(argc, argv);

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
  j = init_spinor_field(VOLUMEPLUSRAND/2, NO_OF_SPINORFIELDS);
  if ( j!= 0) {
    fprintf(stderr, "Not enough memory for spinor fields! Aborting...\n");
    exit(-1);
  }

  g_mu = g_mu1;
  if(g_proc_id == 0){    
    
/*     fscanf(fp6,"%s",filename); */
    /*construct the filenames for the observables and the parameters*/
    strcpy(datafilename,filename);  strcat(datafilename,".data");
    strcpy(parameterfilename,filename);  strcat(parameterfilename,".para");
    
    parameterfile=fopen(parameterfilename, "w");
    write_first_messages(parameterfile, 1);
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


  sprintf(conf_filename,"%s.%.4d", gauge_input_filename, nstore);
  if (g_proc_id == 0) {
    printf("Reading Gauge field from file %s\n", conf_filename); fflush(stdout);
  }
  read_gauge_field(conf_filename, &gaugecksum, &xlfmessage, &gaugelfn);
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

  if(use_stout_flag == 1) {
    if( stout_smear_gauge_field(stout_rho , stout_no_iter) != 0 ) {
      if(g_proc_id == 0) {
	fprintf(stderr, "not enough memory in stou_smear_gauge_field (stout_smear.c)\n");
      }
      exit(-1) ;
    }
    
    plaquette_energy = measure_gauge_action();
    
    if(g_proc_id == 0) {
      printf("The plaquette value after stouting is %e\n", plaquette_energy/(6.*VOLUME*g_nproc)); 
      fflush(stdout);
    }
  }

  if(g_proc_id == 0) {
    printf("Reading source from %s\n", conf_filename);
  }
  read_spinor(g_spinor_field[0], g_spinor_field[1], conf_filename, 0);
    if(g_proc_id == 0) printf("...done\n", conf_filename);

  if(g_proc_id == 0) {printf("mu = %e\n", g_mu/2./g_kappa);}
  
  sprintf(conf_filename, "%s.applied", source_input_filename);
  
  gamma5(g_spinor_field[0], g_spinor_field[0], VOLUME/2);
  gamma5(g_spinor_field[1], g_spinor_field[1], VOLUME/2);
  
#ifdef MPI
  atime = MPI_Wtime();
#endif
  /* Now apply (1+i mu gamma_5)^-1 (M-1) Nmeas times */
  for(j = 0; j < Nmeas/2; j++) {
    M_minus_1_timesC(g_spinor_field[4], g_spinor_field[5],
		     g_spinor_field[0], g_spinor_field[1]);
    M_minus_1_timesC(g_spinor_field[0], g_spinor_field[1],
		     g_spinor_field[4], g_spinor_field[5]);
  }

  gamma5(g_spinor_field[0], g_spinor_field[0], VOLUME/2);
  gamma5(g_spinor_field[1], g_spinor_field[1], VOLUME/2);

  
#ifdef MPI
  etime = MPI_Wtime();
#endif
  
  if(g_proc_id == 0) {
    printf("Wrinting to file %s\n", conf_filename);
  }

  mul_r(g_spinor_field[2], (2*g_kappa), g_spinor_field[0], VOLUME/2);  
  mul_r(g_spinor_field[3], (2*g_kappa), g_spinor_field[1], VOLUME/2);
  construct_writer(&writer, conf_filename);
  sourceFormat = construct_paramsSourceFormat(32, 1, 4, 3);
  write_source_format(writer, sourceFormat);
  write_spinor(writer, &g_spinor_field[2], &g_spinor_field[3], 1, 32);
  free(sourceFormat);

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

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
#include "measure_gauge_action.h"
#ifdef MPI
#include "xchange/xchange.h"
#endif
#include "read_input.h"
#include "mpi_init.h"
#include "sighandler.h"
#include "boundary.h"
#include "solver/solver.h"
#include "init/init.h"
#include "smearing/stout.h"
#include "su3spinor.h"
#include "invert_eo.h"
#include "operator/D_psi.h"
#include "linalg/convert_eo_to_lexic.h"



void usage(){
  fprintf(stdout, "Code for locallity check of the Dirac operator\n");
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

  FILE *parameterfile=NULL;
  int c, j, is=0, ic=0;
  int x, X, y, Y, z, Z, t, tt, i, sum;
  char * filename = NULL;
  char datafilename[50];
  char parameterfilename[50];
  char conf_filename[50];
  char * input_filename = NULL;
  double plaquette_energy, nrm;
  double * norm;
  struct stout_parameters params_smear;
  
#ifdef _GAUGE_COPY
  int kb=0;
#endif
#ifdef MPI
  double atime=0., etime=0.;
#endif
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
  /* here we want no even/odd preconditioning */
  even_odd_flag = 0;

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
    write_first_messages(parameterfile, "check_locality", "NA");
  }

  /* define the geometry */
  geometry();

  /* define the boundary conditions for the fermion fields */
  boundary();

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
  norm = (double*)calloc(3.*LX/2.+T/2., sizeof(double));

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
    xchange_gauge(g_gauge_field);
#endif

    /* Compute minimal eigenvalues, if wanted */
    if(compute_evs != 0) {
      eigenvalues(&no_eigenvalues, 1000, eigenvalue_precision, 0, compute_evs, nstore, even_odd_flag);
    }
    /*compute the energy of the gauge field*/
    plaquette_energy = measure_plaquette(g_gauge_field);

    if(g_proc_id == 0) {
      printf("The plaquette value is %e\n", plaquette_energy/(6.*VOLUME*g_nproc)); fflush(stdout);
    }
    if (use_stout_flag == 1){
      params_smear.rho = stout_rho;
      params_smear.iterations = stout_no_iter;
      if (stout_smear((su3_tuple*)(g_gauge_field[0]), &params_smear, (su3_tuple*)(g_gauge_field[0])) != 0)
        exit(1) ;
      g_update_gauge_copy = 1;
      plaquette_energy = measure_plaquette(g_gauge_field);

      if (g_proc_id == 0) {
        printf("# The plaquette value after stouting is %e\n", plaquette_energy / (6.*VOLUME*g_nproc));
        fflush(stdout);
      }
    }

    source_spinor_field(g_spinor_field[0], g_spinor_field[1], 0, 0);
    convert_eo_to_lexic(g_spinor_field[DUM_DERI], g_spinor_field[0], g_spinor_field[1]);
    D_psi(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI]);
    if(even_odd_flag) {
      i = invert_eo(g_spinor_field[2], g_spinor_field[3], g_spinor_field[0], g_spinor_field[1], 
		    solver_precision, max_solver_iterations, solver_flag, g_relative_precision_flag,
		    sub_evs_cg_flag, even_odd_flag, 0, NULL, -1);
      convert_eo_to_lexic(g_spinor_field[DUM_DERI+1], g_spinor_field[2], g_spinor_field[3]);
    }


    for(i = 0; i < 3*LX/2+T/2; i++){
      norm[i] = 0.;
    }
    
    for(x = 0; x < LX; x++){
      if(x > LX/2) X = LX-x;
      else X = x;
      for(y = 0; y < LY; y++){
	if(y > LY/2) Y = LY-y;
	else Y = y;
	for(z = 0; z < LZ; z++){
	  if(z > LZ/2) Z = LZ-z;
	  else Z = z;
	  for(t = 0; t < T; t++){
	    if(t > T/2) tt = T - t;
	    else tt = t;
	    sum = X + Y + Z + tt;
	    _spinor_norm_sq(nrm, g_spinor_field[DUM_DERI+1][ g_ipt[t][x][y][z] ]);
/* 	    _spinor_norm_sq(nrm, qprop[0][0][1][ g_ipt[t][x][y][z] ]); */
 	    printf("%e %e\n", creal(g_spinor_field[DUM_DERI+1][ g_ipt[t][x][y][z] ].s0.c0), cimag(g_spinor_field[DUM_DERI+1][ g_ipt[t][x][y][z] ].s0.c0));
	    nrm = sqrt( nrm );
	    printf("%1.12e\n", nrm);
	    if(nrm > norm[sum]) norm[sum] = nrm;
	  }
	}
      }
    }
    
    for(i = 0; i < 3*L/2+T/2; i++){
      printf("%d %1.12e\n", i, norm[i]);
    }
    printf("\n");
    
    nstore+=Nsave;
  }

#ifdef MPI
  MPI_Finalize();
#endif
  free_gauge_field();
  free_geometry_indices();
  free_spinor_field();
  free_moment_field();
  return(0);
#ifdef _KOJAK_INST
#pragma pomp inst end(main)
#endif
}

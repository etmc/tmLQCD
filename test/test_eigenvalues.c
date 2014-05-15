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
 *
 *
 * Main for testing the Eigenvalues computation using bispinors 
 *
 * Author: Thomas Chiarappa
 *         Thomas.Chiarappa@mib.infn.it
 *
 *******************************************************************************/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include <signal.h>
#ifdef MPI
#include <mpi.h>
#endif
#include "global.h"
#include "getopt.h"
#include "ranlxd.h"
#include "geometry_eo.h"
#include "start.h"
/*
#include "clover_eo.h"
*/
#include "observables.h"
#include "measure_rectangles.h"
#ifdef MPI
#include "xchange.h"
#endif
#include "io.h"
#include "read_input.h"
#include "mpi_init.h"
#include "sighandler.h"
#include "monomial/moment_energy.h"
#include "update_tm.h"
#include "init/init.h"
#include "test/check_geometry.h"
#include "boundary.h"
#include "polyakov_loop.h"

#include "solver/eigenvalues_bi.h"

char * Version = "2.3.5";


void usage(){
  fprintf(stderr, "hmc for Wilson twisted mass QCD\n\n");
  fprintf(stderr, "Usage: [-f input-filename]\n");
  fprintf(stderr, "Usage: [-o output-filename]\n");
  exit(1);
}

extern int nstore;

int main(int argc,char *argv[]) {
 
  FILE *parameterfile=NULL,*rlxdfile=NULL, *countfile=NULL;
  char * filename = NULL;
  char datafilename[50];
  char parameterfilename[50];
  char gauge_filename[50];
  char * nstore_filename = ".nstore_counter";
  char * input_filename = NULL;
  int rlxd_state[105];
  int j,ix,mu;
  int k;
  struct timeval t1;

  int g_nev, max_iter_ev;
  double stop_prec_ev;


  /* Energy corresponding to the Gauge part */
  double eneg = 0., plaquette_energy = 0., rectangle_energy = 0.;
  /* Acceptance rate */
  int Rate=0;
  /* Do we want to perform reversibility checks */
  /* See also return_check_flag in read_input.h */
  int return_check = 0;
  /* For getopt */
  int c;

  /* For the Polyakov loop: */
  int dir = 2;
  _Complex double pl, pl4;

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

  mpi_init(argc, argv);

  if(Nsave == 0){
    Nsave = 1;
  }
  if(nstore == -1) {
    countfile = fopen(nstore_filename, "r");
    if(countfile != NULL) {
      fscanf(countfile, "%d\n", &nstore);
      fclose(countfile);
    }
    else {
      nstore = 0;
    }
  }
  
  if(g_rgi_C1 == 0.) {
    g_dbw2rand = 0;
  }
#ifndef MPI
  g_dbw2rand = 0;
#endif

  /* Reorder the mu parameter and the number of iterations */
  if(g_mu3 > 0.) {
    g_mu = g_mu1;
    g_mu1 = g_mu3;
    g_mu3 = g_mu;

    j = int_n[1];
    int_n[1] = int_n[3];
    int_n[3] = j;

    j = g_csg_N[0];
    g_csg_N[0] = g_csg_N[4];
    g_csg_N[4] = j;
    g_csg_N[6] = j;
    if(fabs(g_mu3) > 0) {
      g_csg_N[6] = 0;
    }

    g_nr_of_psf = 3;
  }
  else if(g_mu2 > 0.) {
    g_mu = g_mu1;
    g_mu1 = g_mu2;
    g_mu2 = g_mu;

    int_n[3] = int_n[1];
    int_n[1] = int_n[2];
    int_n[2] = int_n[3];

    /* For chronological inverter */
    g_csg_N[4] = g_csg_N[0];
    g_csg_N[0] = g_csg_N[2];
    g_csg_N[2] = g_csg_N[4];
    if(fabs(g_mu2) > 0) {
      g_csg_N[4] = 0;
    }
    g_csg_N[6] = 0;

    g_nr_of_psf = 2;
  }
  else {
    g_csg_N[2] = g_csg_N[0];
    if(fabs(g_mu2) > 0) {
      g_csg_N[2] = 0;
    }
    g_csg_N[4] = 0;
    g_csg_N[6] = 0;
  }

  for(j = 0; j < g_nr_of_psf+1; j++) {
    if(int_n[j] == 0) int_n[j] = 1;
  }
  if(g_nr_of_psf == 3) {
    g_eps_sq_force = g_eps_sq_force1;
    g_eps_sq_force1 = g_eps_sq_force3;
    g_eps_sq_force3 = g_eps_sq_force;
    g_eps_sq_acc = g_eps_sq_acc1;
    g_eps_sq_acc1 = g_eps_sq_acc3;
    g_eps_sq_acc3 = g_eps_sq_acc;
  }
  if(g_nr_of_psf == 2) {
    g_eps_sq_force = g_eps_sq_force1;
    g_eps_sq_force1 = g_eps_sq_force2;
    g_eps_sq_force2 = g_eps_sq_force;
    g_eps_sq_acc = g_eps_sq_acc1;
    g_eps_sq_acc1 = g_eps_sq_acc2;
    g_eps_sq_acc2 = g_eps_sq_acc;
  }
  g_mu = g_mu1;
  g_eps_sq_acc = g_eps_sq_acc1;
  g_eps_sq_force = g_eps_sq_force1;


#ifdef _GAUGE_COPY
  j = init_gauge_field(VOLUMEPLUSRAND + g_dbw2rand, 1);
#else
  j = init_gauge_field(VOLUMEPLUSRAND + g_dbw2rand, 0);
#endif
  if ( j!= 0) {
    fprintf(stderr, "Not enough memory for gauge_fields! Aborting...\n");
    exit(0);
  }
  j = init_geometry_indices(VOLUMEPLUSRAND + g_dbw2rand);
  if ( j!= 0) {
    fprintf(stderr, "Not enough memory for geometry_indices! Aborting...\n");
    exit(0);
  }
  j = init_spinor_field(VOLUMEPLUSRAND/2, NO_OF_SPINORFIELDS);
  if ( j!= 0) {
    fprintf(stderr, "Not enough memory for spinor fields! Aborting...\n");
    exit(0);
  }

  j = init_bispinor_field(VOLUME/2, NO_OF_SPINORFIELDS);


  j = init_csg_field(VOLUMEPLUSRAND/2, g_csg_N);
  if ( j!= 0) {
    fprintf(stderr, "Not enough memory for csg fields! Aborting...\n");
    exit(0);
  }
  j = init_moment_field(VOLUME, VOLUMEPLUSRAND);
  if ( j!= 0) {
    fprintf(stderr, "Not enough memory for moment fields! Aborting...\n");
    exit(0);
  }

  zero_spinor_field(g_spinor_field[DUM_DERI+4],VOLUME/2);
  zero_spinor_field(g_spinor_field[DUM_DERI+5],VOLUME/2);
  zero_spinor_field(g_spinor_field[DUM_DERI+6],VOLUME/2);
 

  if(g_proc_id == 0){
    
/*     fscanf(fp6,"%s",filename); */
    /*construct the filenames for the observables and the parameters*/
    strcpy(datafilename,filename);  strcat(datafilename,".data");
    strcpy(parameterfilename,filename);  strcat(parameterfilename,".para");
    
    parameterfile=fopen(parameterfilename, "w");
    printf("# This is the hmc code for twisted Mass Wilson QCD\n\nVersion %s\n", Version);
#ifdef SSE
    printf("# The code was compiled with SSE instructions\n");
#endif
#ifdef SSE2
    printf("# The code was compiled with SSE2 instructions\n");
#endif
#ifdef SSE3
    printf("# The code was compiled with SSE3 instructions\n");
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
    printf("# beta = %f , kappa= %f\n", g_beta, g_kappa);
    printf("# mus = %f, %f, %f\n", g_mu1, g_mu2, g_mu3);
    printf("# int_n_gauge = %d, int_n_ferm1 = %d, int_n_ferm2 = %d, int_n_ferm3 = %d\n", 
	    int_n[0], int_n[1], int_n[2], int_n[3]);
    printf("# g_rgi_C0 = %f, g_rgi_C1 = %f\n", g_rgi_C0, g_rgi_C1);
    printf("# Number of pseudo-fermion fields: %d\n", g_nr_of_psf);
    printf("# g_eps_sq_force = %e, g_eps_sq_acc = %e\n", g_eps_sq_force, g_eps_sq_acc);
    printf("# Integration scheme: ");
    if(integtyp == 1) printf("leap-frog (single time scale)\n");
    if(integtyp == 2) printf("Sexton-Weingarten (single time scale)\n");
    if(integtyp == 3) printf("leap-frog (multiple time scales)\n");
    if(integtyp == 4) printf("Sexton-Weingarten (multiple time scales)\n");
    if(integtyp == 5) printf("higher order and leap-frog (multiple time scales)\n");
    printf("# Using %s precision for the inversions!\n", 
	   g_relative_precision_flag ? "relative" : "absolute");
    printf("# Using in chronological inverter for spinor_field 1,2,3 a history of %d, %d, %d, respectively\n", 
	   g_csg_N[0], g_csg_N[2], g_csg_N[4]);


    fprintf(parameterfile, "The lattice size is %d x %d x %d x %d\n", (int)(g_nproc_t*T), (int)(g_nproc_x*LX), (int)(LY), (int)(LZ));
    fprintf(parameterfile, "The local lattice size is %d x %d x %d x %d\n", (int)(T), (int)(LX), (int)(LY), (int)(LZ));
    fprintf(parameterfile, "g_beta = %f , g_kappa= %f, g_kappa*csw/8= %f \n",g_beta,g_kappa,g_ka_csw_8);
    fprintf(parameterfile, "boundary of fermion fields (t,x,y,z): %f %f %f %f \n",X0,X1,X2,X3);
    fprintf(parameterfile, "EPS_SQ0=%e, EPS_SQ1=%e EPS_SQ2=%e, EPS_SQ3=%e \n"
	    ,EPS_SQ0,EPS_SQ1,EPS_SQ2,EPS_SQ3);
    fprintf(parameterfile, "g_eps_sq_force = %e, g_eps_sq_acc = %e\n", g_eps_sq_force, g_eps_sq_acc);
    fprintf(parameterfile, "dtau=%f, Nsteps=%d, Nmeas=%d, Nsave=%d, integtyp=%d, nsmall=%d \n",
	    dtau,Nsteps,Nmeas,Nsave,integtyp,nsmall);
    fprintf(parameterfile, "mu = %f, mu2=%f, mu3=%f\n ", g_mu, g_mu2, g_mu3);
    fprintf(parameterfile, "int_n_gauge = %d, int_n_ferm1 = %d, int_n_ferm2 = %d, int_n_ferm3 = %d\n ", 
	    int_n[0], int_n[1], int_n[2], int_n[3]);
    fprintf(parameterfile, "g_rgi_C0 = %f, g_rgi_C1 = %f\n", g_rgi_C0, g_rgi_C1);
    fprintf(parameterfile, "# Number of pseudo-fermion fields: %d\n", g_nr_of_psf);
    fprintf(parameterfile, "# Integration scheme: ");
    if(integtyp == 1) fprintf(parameterfile, "leap-frog (single time scale)\n");
    if(integtyp == 2) fprintf(parameterfile, "Sexton-Weingarten (single time scale)\n");
    if(integtyp == 3) fprintf(parameterfile, "leap-frog (multiple time scales)\n");
    if(integtyp == 4) fprintf(parameterfile, "Sexton-Weingarten (multiple time scales)\n");
    if(integtyp == 5) fprintf(parameterfile, "higher order and leap-frog (multiple time scales)\n");
    fprintf(parameterfile, "Using %s precision for the inversions!\n", 
	   g_relative_precision_flag ? "relative" : "absolute");
    fprintf(parameterfile, "Using in chronological inverter for spinor_field 1,2,3 a history of %d, %d, %d, respectively\n", 
	   g_csg_N[0], g_csg_N[2], g_csg_N[4]);
    fflush(stdout); fflush(parameterfile);
  }

  /* define the geometry */
  geometry();

  /* define the boundary conditions for the fermion fields */
  boundary();

  check_geometry();

  if(g_proc_id == 0) {
#if defined GEOMETRIC
    if(g_proc_id==0) fprintf(parameterfile,"The geometric series is used as solver \n\n");
#else
    if(g_proc_id==0) fprintf(parameterfile,"The BICG_stab is used as solver \n\n");
#endif
    fflush(parameterfile);
  }
  
  /* Continue */
  if(startoption == 3){
    rlxdfile = fopen(rlxd_input_filename,"r");
    if(rlxdfile != NULL) {
      if(g_proc_id == 0) {
	fread(rlxd_state,sizeof(rlxd_state),1,rlxdfile);
      }
    }
    else {
      if(g_proc_id == 0) {
	printf("%s does not exist, switching to restart...\n", rlxd_input_filename);
      }
      startoption = 2;
    }
    fclose(rlxdfile);
    if(startoption != 2) {
      if(g_proc_id == 0) {
	rlxd_reset(rlxd_state);
	printf("Reading Gauge field from file %s\n", gauge_input_filename); fflush(stdout);
      }
      
      read_gauge_field_time_p(gauge_input_filename,g_gauge_field);
    }
  }
  if(startoption != 3){
    /* Initialize random number generator */
    if(g_proc_id == 0) {
      rlxd_init(1, random_seed);
      /* hot */
      if(startoption == 1) {
	random_gauge_field();
      }
      rlxd_get(rlxd_state);
#ifdef MPI
      MPI_Send(&rlxd_state[0], 105, MPI_INT, 1, 99, MPI_COMM_WORLD);
      MPI_Recv(&rlxd_state[0], 105, MPI_INT, g_nproc-1, 99, MPI_COMM_WORLD, &status);
      rlxd_reset(rlxd_state);
#endif
    }
#ifdef MPI
    else {
      MPI_Recv(&rlxd_state[0], 105, MPI_INT, g_proc_id-1, 99, MPI_COMM_WORLD, &status);
      rlxd_reset(rlxd_state);
      /* hot */
      if(startoption == 1) {
	random_gauge_field();
      }
      k=g_proc_id+1; 
      if(k==g_nproc){
	k=0;
      }
      rlxd_get(rlxd_state);
      MPI_Send(&rlxd_state[0], 105, MPI_INT, k, 99, MPI_COMM_WORLD);
    }
#endif

    /* Cold */
    if(startoption == 0) {
      unit_g_gauge_field();
    }
    /* Restart */
    else if(startoption == 2) {
      if (g_proc_id == 0){
	printf("Reading Gauge field from file %s\n", gauge_input_filename); fflush(stdout);
      }
      read_gauge_field_time_p(gauge_input_filename,g_gauge_field);
    }

  }

  /*For parallelization: exchange the gaugefield */
#ifdef MPI
  xchange_gauge(g_gauge_field);
#endif
#ifdef _GAUGE_COPY
  update_backward_gauge();
#endif

  /*compute the energy of the gauge field*/
  plaquette_energy=measure_gauge_action();
  if(g_rgi_C1 > 0. || g_rgi_C1 < 0.) {
    rectangle_energy = measure_rectangles();
    if(g_proc_id==0){
      fprintf(parameterfile,"#First rectangle value: %14.12f \n",rectangle_energy/(12.*VOLUME*g_nproc));
    }
  }
  eneg = g_rgi_C0 * plaquette_energy + g_rgi_C1 * rectangle_energy;
  
  /* Measure and print the Polyakov loop: */
  polyakov_loop(&pl, dir);

  if(g_proc_id==0){
    fprintf(parameterfile,"#First plaquette value: %14.12f \n", plaquette_energy/(6.*VOLUME*g_nproc));
    fprintf(parameterfile,"#First Polyakov loop value in %d-direction |L(%d)|= %14.12f \n",
	    dir, dir, cabs(pl));
  }

  dir=3;
  polyakov_loop(&pl, dir);
  if(g_proc_id==0){
    fprintf(parameterfile,"#First Polyakov loop value in %d-direction |L(%d)|= %14.12f \n",
	    dir, dir, cabs(pl));
    fclose(parameterfile);
  }

  /* set ddummy to zero */
  for(ix = 0; ix < VOLUME+RAND; ix++){
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

  if(g_proc_id == 0) {
    gettimeofday(&t1,NULL);
    countfile = fopen("history_hmc_tm", "a");
    fprintf(countfile, "!!! Timestamp %ld, Nsave = %d, g_mu = %e, g_mu1 = %e, g_mu_2 = %e, g_mu3 = %e, beta = %f, kappa = %f, C1 = %f, int0 = %d, int1 = %d, int2 = %d, int3 = %d, g_eps_sq_force = %e, g_eps_sq_acc = %e, ", 
	    t1.tv_sec, Nsave, g_mu, g_mu1, g_mu2, g_mu3, g_beta, g_kappa, g_rgi_C1, 
	    int_n[0], int_n[1], int_n[2], int_n[3], g_eps_sq_force, g_eps_sq_acc); 
    fprintf(countfile, "Nsteps = %d, dtau = %e, tau = %e, integtyp = %d, rel. prec. = %d\n", 
	    Nsteps, dtau, tau, integtyp, g_relative_precision_flag);
    fclose(countfile);
  }



     /* HERE THE CALLS FOR SOME EIGENVALUES */

  /* for lowest
  g_nev = 10;
  */

  /* for largest
  */
  g_nev = 10;

  max_iter_ev = 1000;
  stop_prec_ev = 1.e-10;

  if(g_proc_id==0) {

  printf(" Values of   mu = %e     mubar = %e     eps = %e     precision = %e  \n \n", g_mu, g_mubar, g_epsbar, stop_prec_ev);

  }

  eigenvalues(&g_nev, operator_flag, max_iter_ev, stop_prec_ev);

  g_nev = 4;

  max_iter_ev = 200;
  stop_prec_ev = 1.e-03;

  max_eigenvalues(&g_nev, operator_flag, max_iter_ev, stop_prec_ev);

  if(g_proc_id==0) {

  printf(" Values of   mu = %e     mubar = %e     eps = %e     precision = %e  \n \n", g_mu, g_mubar, g_epsbar, stop_prec_ev);

  /*
  printf(" Values of   mu = %e     precision = %e  \n \n", g_mu, stop_prec_ev);
  */

  }

   /* END OF EIGENVALUES CALLS */


  if(g_proc_id==0) {
    rlxd_get(rlxd_state);
    rlxdfile=fopen("last_state","w");
    fwrite(rlxd_state,sizeof(rlxd_state),1,rlxdfile);
    fclose(rlxdfile);

    printf("Acceptance Rate was: %e Prozent\n", 100.*(double)Rate/(double)Nmeas);
    fflush(stdout);
    parameterfile = fopen(parameterfilename, "a");
    fprintf(parameterfile, "Acceptance Rate was: %e Prozent\n", 100.*(double)Rate/(double)Nmeas);
    fclose(parameterfile);
  }
#ifdef MPI
  MPI_Finalize();
#endif
  free_gauge_tmp();
  free_gauge_field();
  free_geometry_indices();
  free_spinor_field();
  free_bispinor_field();  
  free_moment_field();
  return(0);
}


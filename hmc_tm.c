/* $Id$ */
/*******************************************************************************
*
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
#include <sys/time.h>
#include <string.h>
#include <signal.h>
#ifdef MPI
# include <mpi.h>
#endif
#include "global.h"
#include "getopt.h"
#include "ranlxd.h"
#include "geometry_eo.h"
#include "start.h"
#include "observables.h"
#include "measure_rectangles.h"
#ifdef MPI
# include "xchange.h"
#endif
#include "io.h"
#include "read_input.h"
#include "mpi_init.h"
#include "sighandler.h"
/* #include "hybrid_update.h" */
#include "update_tm.h"
#include "update_backward_gauge.h"
#include "init_gauge_field.h"
#include "init_geometry_indices.h"
#include "init_spinor_field.h"
#include "init_moment_field.h"
#include "init_gauge_tmp.h"
#include "init_dirac_halfspinor.h"
#include "xchange_halffield.h"
#include "test/check_geometry.h"
#include "boundary.h"
#include "phmc.h"
#include "solver/solver.h"
#include "polyakov_loop.h"

void usage(){
  fprintf(stdout, "HMC for Wilson twisted mass QCD\n");
  fprintf(stdout, "Version %s \n\n", PACKAGE_VERSION);
  fprintf(stdout, "Please send bug reports to %s\n", PACKAGE_BUGREPORT);
  fprintf(stdout, "Usage:   hmc_tm [options]\n");
  fprintf(stdout, "Options: [-f input-filename]  default: hmc.input\n");
  fprintf(stdout, "         [-o output-filename] default: output\n");
  fprintf(stdout, "         [-h|-? this help]\n");
  exit(0);
}


extern int nstore;
const int rlxdsize = 105;

int main(int argc,char *argv[]) {
 
  FILE *parameterfile=NULL,*rlxdfile=NULL, *countfile=NULL;
  char * filename = NULL;
  char datafilename[50];
  char parameterfilename[50];
  char * gauge_filename = "conf.save";
  char * nstore_filename = ".nstore_counter";
  char * tmp_filename = ".conf.tmp";
  char * input_filename = NULL;
  char command_string[300];
  int rlxd_state[105];
  int j,ix,mu, trajectory_counter=1;
  int k;
  struct timeval t1;
  double x;

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
  complex pl, pl4;
#ifdef _KOJAK_INST
#pragma pomp inst init
#pragma pomp inst begin(main)
#endif

  g_running_phmc = 0;
  DUM_DERI = 6;
  DUM_SOLVER = DUM_DERI+7;
  DUM_MATRIX = DUM_SOLVER+6;
  NO_OF_SPINORFIELDS = DUM_MATRIX+6;

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

  mpi_init(argc, argv);

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

    x = lambda[1];
    lambda[1] = lambda[3];
    lambda[3] = x;

    j = g_csg_N[0];
    g_csg_N[0] = g_csg_N[4];
    g_csg_N[4] = j;
    g_csg_N[6] = j;
    if(ITER_MAX_BCG == 0 || fabs(g_mu3) > 0) {
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

    lambda[3] = lambda[1];
    lambda[1] = lambda[2];
    lambda[2] = lambda[3];

    /* For chronological inverter */
    g_csg_N[4] = g_csg_N[0];
    g_csg_N[0] = g_csg_N[2];
    g_csg_N[2] = g_csg_N[4];
    if(ITER_MAX_BCG == 0 || fabs(g_mu2) > 0) {
      g_csg_N[4] = 0;
    }
    g_csg_N[6] = 0;

    g_nr_of_psf = 2;
  }
  else {
    g_csg_N[2] = g_csg_N[0];
    if(ITER_MAX_BCG == 0 || fabs(g_mu2) > 0) {
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
  if(even_odd_flag)
  {
    j = init_spinor_field(VOLUMEPLUSRAND/2, NO_OF_SPINORFIELDS);
  }
  else
  {
    j = init_spinor_field(VOLUMEPLUSRAND, NO_OF_SPINORFIELDS);
  }
  if ( j!= 0) {
    fprintf(stderr, "Not enough memory for spinor fields! Aborting...\n");
    exit(0);
  }
  if(even_odd_flag)
  {
    j = init_csg_field(VOLUMEPLUSRAND/2, g_csg_N);
  }
  else
  {
    j = init_csg_field(VOLUMEPLUSRAND, g_csg_N);
  }
  if ( j!= 0) {
    fprintf(stderr, "Not enough memory for csg fields! Aborting...\n");
    exit(0);
  }
  j = init_moment_field(VOLUME, VOLUMEPLUSRAND);
  if ( j!= 0) {
    fprintf(stderr, "Not enough memory for moment fields! Aborting...\n");
    exit(0);
  }

  zero_spinor_field(g_spinor_field[DUM_DERI+4],VOLUME);
  zero_spinor_field(g_spinor_field[DUM_DERI+5],VOLUME);
  zero_spinor_field(g_spinor_field[DUM_DERI+6],VOLUME);

  if(use_stout_flag == 1)
    init_stout_smear_vars(VOLUMEPLUSRAND, stout_no_iter);

  /*construct the filenames for the observables and the parameters*/
  strcpy(datafilename,filename);  strcat(datafilename,".data");
  strcpy(parameterfilename,filename);  strcat(parameterfilename,".para");
  
  if(g_proc_id == 0){    
    parameterfile = fopen(parameterfilename, "a");
    write_first_messages(parameterfile, integtyp, 0);
  }
  /* define the geometry */
  geometry();
  
  /* define the boundary conditions for the fermion fields */
  boundary();
  
  check_geometry();

#ifdef _USE_HALFSPINOR
  j = init_dirac_halfspinor();
  if ( j!= 0) {
    fprintf(stderr, "Not enough memory for halffield! Aborting...\n");
    exit(0);
  }
  if(g_sloppy_precision_flag == 1) {
    init_dirac_halfspinor32();
  }
#  if (defined _PERSISTENT)
  init_xchange_halffield();
#  endif
#endif

  /* Initialise random number generator */
  /* Continue */
  if(startoption == 3) {
    if( (j = read_rlxd_state(gauge_input_filename, rlxd_state, rlxdsize)) == -1) {
      if(g_proc_id == 0) {
	printf("no rlxd_state found in %s, switching to restart...\n", gauge_input_filename);
	fflush(stdout);
      }
      startoption = 2;
    }
    else {
      rlxd_reset(rlxd_state);
    }
  }
  /* restart, hot and cold */
  if(startoption != 3) {
    rlxd_init(1, random_seed + g_proc_id*97);
  }
  
  /* Set up the gauge field */
  /* continue and restart */
  if(startoption==3 || startoption == 2) {
    if(g_proc_id == 0) {
      printf("# Reading Gauge field from file %s in %d Bit\n", 
	     gauge_input_filename, gauge_precision_read_flag); 
      fflush(stdout);
    }
    if(gauge_precision_read_flag == 64) {
      read_lime_gauge_field(gauge_input_filename);
    }
    else if(gauge_precision_read_flag == 32){
      read_lime_gauge_field_singleprec(gauge_input_filename);
    }
    if (g_proc_id == 0){
      printf("# done!\n"); fflush(stdout);
    }
  }
  else if (startoption == 1) {
    /* hot */
    random_gauge_field();
  }
  else if(startoption == 0) {
    /* cold */
    unit_g_gauge_field();    
  }
  

  /*For parallelization: exchange the gaugefield */
#ifdef MPI
  xchange_gauge();
#endif
#ifdef _GAUGE_COPY
  update_backward_gauge();
#endif

  /*compute the energy of the gauge field*/
  plaquette_energy=measure_gauge_action();
  if(g_rgi_C1 > 0. || g_rgi_C1 < 0.) {
    rectangle_energy = measure_rectangles();
    if(g_proc_id==0){
      fprintf(parameterfile,"# First rectangle value: %14.12f \n",rectangle_energy/(12.*VOLUME*g_nproc));
    }
  }
  eneg = g_rgi_C0 * plaquette_energy + g_rgi_C1 * rectangle_energy;
  
  /* Measure and print the Polyakov loop: */
  polyakov_loop(&pl, dir);

  if(g_proc_id==0){
    fprintf(parameterfile,"# First plaquette value: %14.12f \n", plaquette_energy/(6.*VOLUME*g_nproc));
    printf("# First plaquette value: %14.12f \n", plaquette_energy/(6.*VOLUME*g_nproc));
    fprintf(parameterfile,"# First Polyakov loop value in %d-direction |L(%d)|= %14.12f \n",
	    dir, dir, sqrt(pl.re*pl.re+pl.im*pl.im));
  }

  dir=3;
  polyakov_loop(&pl, dir);
  if(g_proc_id==0){
    fprintf(parameterfile,"# First Polyakov loop value in %d-direction |L(%d)|= %14.12f \n",
	    dir, dir, sqrt(pl.re*pl.re+pl.im*pl.im));
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
    fprintf(countfile, "!!! Timestamp %ld, Nskip = %d, g_mu = %e, g_mu1 = %e, g_mu_2 = %e, g_mu3 = %e, beta = %f, kappa = %f, C1 = %f, int0 = %d, int1 = %d, int2 = %d, int3 = %d, g_eps_sq_force = %e, g_eps_sq_acc = %e, ", 
	    t1.tv_sec, Nskip, g_mu, g_mu1, g_mu2, g_mu3, g_beta, g_kappa, g_rgi_C1, 
	    int_n[0], int_n[1], int_n[2], int_n[3], g_eps_sq_force, g_eps_sq_acc); 
    fprintf(countfile, "Nsteps = %d, dtau = %e, tau = %e, integtyp = %d, rel. prec. = %d\n", 
	    Nsteps, dtau, tau, integtyp, g_relative_precision_flag);
    fclose(countfile);
  }

  /* Loop for measurements */
  for(j=0;j<Nmeas;j++) {
    if(g_proc_id == 0) {
      printf("# Starting trajectory no %d\n", trajectory_counter);
    }

    if(return_check_flag == 1 && (j+1)%return_check_interval == 0) return_check = 1;
    else return_check = 0;

    Rate += update_tm(integtyp, &plaquette_energy, &rectangle_energy, datafilename, 
		      dtau, Nsteps, nsmall, tau, int_n, return_check, lambda, reproduce_randomnumber_flag);


    /* Measure the Polyakov loop in direction 2 and 3:*/
    polyakov_loop(&pl, 2); 
    polyakov_loop(&pl4, 3);  
    
    /* Save gauge configuration all Nskip times */
    if((Nskip !=0) && (trajectory_counter%Nskip == 0) && (trajectory_counter!=0)) {
      sprintf(gauge_filename,"conf.%.4d", nstore);
      if(g_proc_id == 0) {
        countfile = fopen("history_hmc_tm", "a");
	fprintf(countfile, "%.4d, measurement %d of %d, Nskip = %d, Plaquette = %e, |L(%d)| = %e, |L(%d)| = %e trajectory nr = %d\n", 
		nstore, j, Nmeas, Nskip, plaquette_energy/(6.*VOLUME*g_nproc),
		2, sqrt(pl.re*pl.re+pl.im*pl.im),
		dir, sqrt(pl4.re*pl4.re+pl4.im*pl4.im), trajectory_counter);
	fclose(countfile);
      }
      nstore ++;
    }
    else if(write_cp_flag == 1 || (j >= (Nmeas -1))) {
      sprintf(gauge_filename,"conf.save");
    }
    
    if(((Nskip !=0) && (trajectory_counter%Nskip == 0) && (trajectory_counter!=0)) || (write_cp_flag == 1) || (j >= (Nmeas - 1))) {
      /* Write the gauge configuration first to a temporary file */
      if(gauge_precision_write_flag == 64) {
  	write_lime_gauge_field( tmp_filename , plaquette_energy/(6.*VOLUME*g_nproc), trajectory_counter);
      }
      else if(gauge_precision_write_flag == 32) {
	write_lime_gauge_field_singleprec( tmp_filename , plaquette_energy/(6.*VOLUME*g_nproc), trajectory_counter);
      }
      /*  write the status of the random number generator to a file */
      if(g_proc_id==0) {
	rlxd_get(rlxd_state);
 	write_rlxd_state(tmp_filename, rlxd_state, rlxdsize);
      }
      
      /* Now move it! */
      if(g_proc_id == 0) {
  	rename(tmp_filename, gauge_filename);
        countfile = fopen(nstore_filename, "w");
        fprintf(countfile, "%d %d %s\n", nstore, trajectory_counter+1, gauge_filename);
        fclose(countfile);
      }
    }
    if(g_proc_id == 0) {
      verbose = 1;
    }
    ix = reread_input("hmc.reread");
    if(g_proc_id == 0) {
      verbose = 0;
    }

#ifdef MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    if(ix == 0 && g_proc_id == 0) {
      countfile = fopen("history_hmc_tm", "a");
      fprintf(countfile, "# Changed parameter according to hmc.reread: measurment %d of %d\n", j, Nmeas); 
      fclose(countfile);
      printf("# Changed parameter according to hmc.reread (see stdout): measurment %d of %d\n", j, Nmeas); 
      remove("hmc.reread");
    }
    trajectory_counter++;
  }

  if(g_proc_id==0) {
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
  free_moment_field();
  if(use_stout_flag == 1)
    free_stout_smear_vars();

  return(0);
#ifdef _KOJAK_INST
#pragma pomp inst end(main)
#endif
}

static char const rcsid[] = "$Id$";

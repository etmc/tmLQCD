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
#include "hybrid_update.h"
#include "update_tm.h"
#include "init_gauge_field.h"
#include "init_geometry_indices.h"
#include "init_spinor_field.h"
#include "init_moment_field.h"
#include "init_gauge_tmp.h"
#include "test/check_geometry.h"
#include "boundary.h"
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

int main(int argc,char *argv[]) {
 
  FILE *parameterfile=NULL,*rlxdfile=NULL, *countfile=NULL;
  char * filename = NULL;
  char datafilename[50];
  char parameterfilename[50];
  char gauge_filename[50];
  char * nstore_filename = ".nstore_counter";
  char * input_filename = NULL;
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

  if(Nskip == 0){
    Nskip = 1;
  }
  if(nstore == -1) {
    countfile = fopen(nstore_filename, "r");
    if(countfile != NULL) {
      fscanf(countfile, "%d %d\n", &nstore, &trajectory_counter);
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
  j = init_spinor_field(VOLUMEPLUSRAND/2, NO_OF_SPINORFIELDS);
  if ( j!= 0) {
    fprintf(stderr, "Not enough memory for spinor fields! Aborting...\n");
    exit(0);
  }
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
    printf("# This is the hmc code for twisted Mass Wilson QCD\n\nVersion %s\n", PACKAGE_VERSION);
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
	   (int)(T*g_nproc_t), (int)(LX*g_nproc_x), (int)(LY*g_nproc_y), (int)(LZ*g_nproc_z));
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
    if(integtyp == 6) printf("second order minimal norm (velocity version, multiple time scales)\n");
    if(integtyp == 7) printf("second order minimal norm (position version, multiple time scales)\n");
    printf("# Using %s precision for the inversions!\n", 
	   g_relative_precision_flag ? "relative" : "absolute");
    printf("# Using in chronological inverter for spinor_field 1,2,3 a history of %d, %d, %d, respectively\n", 
	   g_csg_N[0], g_csg_N[2], g_csg_N[4]);


    fprintf(parameterfile, "The lattice size is %d x %d x %d x %d\n", (int)(g_nproc_t*T), (int)(g_nproc_x*LX), 
	    (int)(g_nproc_y*LY), (int)(g_nproc_z*LZ));
    fprintf(parameterfile, "The local lattice size is %d x %d x %d x %d\n", (int)(T), (int)(LX), (int)(LY), (int)(LZ));
    fprintf(parameterfile, "g_beta = %f , g_kappa= %f, g_kappa*csw/8= %f \n",g_beta,g_kappa,g_ka_csw_8);
    fprintf(parameterfile, "boundary of fermion fields (t,x,y,z): %f %f %f %f \n",X0,X1,X2,X3);
    fprintf(parameterfile, "ITER_MAX_BCG=%d, EPS_SQ0=%e, EPS_SQ1=%e EPS_SQ2=%e, EPS_SQ3=%e \n"
	    ,ITER_MAX_BCG,EPS_SQ0,EPS_SQ1,EPS_SQ2,EPS_SQ3);
    fprintf(parameterfile, "g_eps_sq_force = %e, g_eps_sq_acc = %e\n", g_eps_sq_force, g_eps_sq_acc);
    fprintf(parameterfile, "dtau=%f, Nsteps=%d, Nmeas=%d, Nskip=%d, integtyp=%d, nsmall=%d \n",
	    dtau,Nsteps,Nmeas,Nskip,integtyp,nsmall);
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
      read_lime_gauge_field(gauge_input_filename);
/*       read_gauge_field_time_p(gauge_input_filename); */
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
      read_lime_gauge_field(gauge_input_filename);
/*       read_gauge_field_time_p(gauge_input_filename); */
    }

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
      fprintf(parameterfile,"#First rectangle value: %14.12f \n",rectangle_energy/(12.*VOLUME*g_nproc));
    }
  }
  eneg = g_rgi_C0 * plaquette_energy + g_rgi_C1 * rectangle_energy;
  
  /* Measure and print the Polyakov loop: */
  polyakov_loop(&pl, dir);

  if(g_proc_id==0){
    fprintf(parameterfile,"#First plaquette value: %14.12f \n", plaquette_energy/(6.*VOLUME*g_nproc));
    fprintf(parameterfile,"#First Polyakov loop value in %d-direction |L(%d)|= %14.12f \n",
	    dir, dir, sqrt(pl.re*pl.re+pl.im*pl.im));
  }

  dir=3;
  polyakov_loop(&pl, dir);
  if(g_proc_id==0){
    fprintf(parameterfile,"#First Polyakov loop value in %d-direction |L(%d)|= %14.12f \n",
	    dir, dir, sqrt(pl.re*pl.re+pl.im*pl.im));
    fclose(parameterfile);
  }

  /* compute the energy of the determinant term */
  /* needed for exact continuation of the run, since evamax and eva use
     random numbers */ 
  if(startoption == 2 && g_proc_id == 0){
    rlxd_reset(rlxd_state);
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
    if(return_check_flag == 1 && (j+1)%return_check_interval == 0) return_check = 1;
    else return_check = 0;

    Rate += update_tm(integtyp, &plaquette_energy, &rectangle_energy, datafilename, 
		      dtau, Nsteps, nsmall, tau, int_n, return_check, lambda);

    /* Measure the Polyakov loop in direction 2 and 3:*/
    polyakov_loop(&pl, 2); 
    polyakov_loop(&pl4, 3);  
    
    /* Save gauge configuration all Nskip times */
    if((trajectory_counter%Nskip == 0) && (trajectory_counter!=0)) {
      sprintf(gauge_filename,"%s.%.4d", "conf", nstore);
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
    else {
      sprintf(gauge_filename,"%s", "conf.save");
    }
    countfile = fopen(nstore_filename, "w");
    fprintf(countfile, "%d %d\n", nstore, trajectory_counter+1);
    fclose(countfile);
    verbose = 1;
    ix = reread_input("hmc.reread");
    verbose = 0;

    write_lime_gauge_field( gauge_filename , plaquette_energy/(6.*VOLUME*g_nproc), trajectory_counter);
    /*  write the status of the random number generator on a file */
    if(g_proc_id==0) {
      rlxd_get(rlxd_state);
      rlxdfile=fopen("rlxd_state","w");
      fwrite(rlxd_state,sizeof(rlxd_state),1,rlxdfile);
      fclose(rlxdfile);
    }

#ifdef MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    if(ix == 0 && g_proc_id == 0) {
      countfile = fopen("history_hmc_tm", "a");
      fprintf(countfile, "# Changed parameter according to hmc.reread: measurment %d of %d\n", j, Nmeas); 
      fclose(countfile);
      printf("# Changed parameter according to hmc.reread (see stdout): measurment %d of %d\n", j, Nmeas); 
      system("rm hmc.reread");
    }
    trajectory_counter++;
  }
  /* write the gauge configuration to the file last_configuration */
  write_lime_gauge_field( "last_configuration" , plaquette_energy/(6.*VOLUME*g_nproc), trajectory_counter);

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
  free_moment_field();
  return(0);
}

static char const rcsid[] = "$Id$";

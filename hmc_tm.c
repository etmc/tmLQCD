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
#include "ranlxd.h"
#include "geometry_eo.h"
#include "start.h"
#include "clover_eo.h"
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

char * Version = "1.2";

int check_geometry();

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

  /* Energy corresponding to the Gauge part */
  double eneg=0.;
  /* Acceptance rate */
  int Rate=0;

  int c;

  verbose = 0;
  g_use_clover_flag = 0;
  g_nr_of_psf = 1;
 
  signal(SIGUSR1,&catch_del_sig);
  signal(SIGUSR2,&catch_del_sig);
  signal(SIGTERM,&catch_del_sig);
  signal(SIGXCPU,&catch_del_sig);
  
  mpi_init(argc, argv);

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

  q_off = 0.;
  q_off2 = 0.;
 
  /* Reorder the mu parameter */
  if(g_mu3 > 0.) {
    g_mu = g_mu1;
    g_mu1 = g_mu3;
    g_mu3 = g_mu;
    g_nr_of_psf = 3;
  }
  else if(g_mu2 > 0.) {
    g_mu = g_mu1;
    g_mu1 = g_mu2;
    g_mu2 = g_mu;
    g_nr_of_psf = 2;
  }
  g_mu = g_mu1;

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
#ifdef P4
    printf("# The code was compiled for pentium4\n");
#endif
    printf("# The lattice size is %d x %d^3\n",(int)(T)*g_nproc,(int)(L));
    printf("# The local lattice size is %d x %d^3\n",(int)(T),(int)(L));
    printf("# beta = %f , kappa= %f, mu= %f \n",g_beta,g_kappa,g_mu);
    printf("# mus = %f, %f, %f\n", g_mu1, g_mu2, g_mu3);
    
    fprintf(parameterfile, "The lattice size is %d x %d^3\n", (int)(g_nproc*T),(int)(L));
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

  ix = check_geometry();

  if(g_proc_id == 0) {
#if defined GEOMETRIC
    if(g_proc_id==0) fprintf(parameterfile,"The geometric series is used as solver \n\n");
#else
    if(g_proc_id==0) fprintf(parameterfile,"The BICG_stab is used as solver \n\n");
#endif
  }
  
  /* Continue */
  if(startoption == 3){
    if (g_proc_id == 0){
      rlxdfile=fopen(rlxd_input_filename,"r");
      fread(rlxd_state,sizeof(rlxd_state),1,rlxdfile);
      fclose(rlxdfile);
      rlxd_reset(rlxd_state);
      printf("Reading Gauge field from file %s\n", gauge_input_filename); fflush(stdout);
    }

    read_gauge_field_time_p(gauge_input_filename);
  }
  else {
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
      read_gauge_field_time_p(gauge_input_filename);
    }

  }

  /*For parallelization: exchange the gaugefield */
#ifdef MPI
  xchange_gauge();
#endif

  /*compute the energy of the gauge field*/
  eneg=measure_gauge_action();

  if(g_proc_id==0){
    fprintf(parameterfile,"#First plaquette value: %14.12f \n",eneg/(6.*VOLUME*g_nproc));
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

  /* Loop for measurements */
  for(j=0;j<Nmeas;j++) {

    Rate += update_tm(integtyp, &eneg, datafilename);

    /* Save gauge configuration all Nskip times */
    if((j+1)%Nskip == 0) {
      sprintf(gauge_filename,"%s.%.4d", "conf", nstore);
      nstore ++;
      countfile = fopen(nstore_filename, "w");
      fprintf(countfile, "%d\n", nstore);
      fclose(countfile);
      write_gauge_field_time_p( gauge_filename );

      /*  write the status of the random number generator on a file */
      if(g_proc_id==0) {
	rlxd_get(rlxd_state);
	rlxdfile=fopen("rlxd_state","w");
	fwrite(rlxd_state,sizeof(rlxd_state),1,rlxdfile);
	fclose(rlxdfile);
      }
    }
  }
  /* write the gauge configuration to the file last_configuration */
  write_gauge_field_time_p( "last_configuration" );
  
  if(g_proc_id==0) {
    rlxd_get(rlxd_state);
    rlxdfile=fopen("last_state","w");
    fwrite(rlxd_state,sizeof(rlxd_state),1,rlxdfile);
    fclose(rlxdfile);
  }
  if(g_proc_id == 0) {
    printf("Acceptance Rate was: %e Prozent\n", 100.*(double)Rate/(double)Nmeas);
    fflush(stdout);
    parameterfile = fopen(parameterfilename, "a");
    fprintf(parameterfile, "Acceptance Rate was: %e Prozent\n", 100.*(double)Rate/(double)Nmeas);
    fclose(parameterfile);
  }
#ifdef MPI
  MPI_Finalize();
#endif
  return(0);
}

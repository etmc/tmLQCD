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
#ifdef MPI
#include <mpi.h>
#endif
#include "getopt.h"
#include "global.h"
#include "su3.h"
#include "su3adj.h"
#include "expo.h"
#include "ranlxd.h"
#include "geometry_eo.h"
#include "start.h"
#include "linalg_eo.h"
#include "linsolve.h"
#include "clover_eo.h"
#include "observables.h"
#include "hybrid_update.h"
#ifdef MPI
#include "xchange.h"
#endif
#include "sw.h"
#include "io.h"
#include "read_input.h"
#include "tm_operators.h"
#include "bicgstabell.h"
#include "boundary.h"

char * Version = "0.9";

void usage(){
  fprintf(stderr, "hmc for Wilson twisted mass QCD\n\n");
  fprintf(stderr, "Usage: [-f input-filename]\n");
  fprintf(stderr, "Usage: [-o output-filename]\n");
  exit(1);
}

extern int nstore;

int main(int argc,char *argv[]) {
 
  FILE *datafile=NULL,*parameterfile=NULL,*rlxdfile=NULL;
  char * filename = NULL;
  char filename1[50];
  char filename2[50];
  char filename3[50];
  int rlxd_state[105];
  int idis0;
  int j,ix,mu;
  int i,k;

  double yy[1];
/*   static double step; */
/*   static double q_off,q_off2; */
  double eneg, enegx, enep, enepx, dh;
  double enerphi0,enerphi0x;
  static su3 gauge_tmp[VOLUME][4];
  su3 *v,*w;

  int Rate=0;
#ifdef MPI
  int  namelen;
  char processor_name[MPI_MAX_PROCESSOR_NAME];
#endif
  char * input_filename = NULL;
  int c;

  verbose = 0;
  g_use_clover_flag = 0;
 
  
#ifdef MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &g_nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &g_proc_id);
  MPI_Get_processor_name(processor_name, &namelen);
  g_cart_grid = MPI_COMM_WORLD;

  fprintf(stdout,"Process %d of %d on %s\n",
	  g_proc_id, g_nproc, processor_name);
  fflush(stdout);
#else
  g_nproc = 1;
  g_proc_id = 0;
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
 
  q_off = 0.;
  q_off2 = 0.;
 
  if(g_proc_id == 0){
    
/*     fscanf(fp6,"%s",filename); */
    /*construct the filenames for the observables and the parameters*/
    strcpy(filename1,filename);  strcat(filename1,".data");
    strcpy(filename2,filename);  strcat(filename2,".para");
    
    datafile=fopen(filename1, "w");
    parameterfile=fopen(filename2, "w");
    
    if(g_proc_id == 0){
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
    }
    printf("The lattice size is %d x %d^3\n",(int)(T)*g_nproc,(int)(L));
    printf("The local lattice size is %d x %d^3\n",(int)(T),(int)(L));
    printf("g_beta = %f , g_kappa= %f, g_kappa*csw/8= %f \n\n",g_beta,g_kappa,g_ka_csw_8);
    
    fprintf(parameterfile,"The lattice size is %d x %d^3\n\n",(int)(g_nproc*T),(int)(L));
    fprintf(parameterfile,"g_beta = %f , g_kappa= %f, g_kappa*csw/8= %f g_mu = %f \n\n",g_beta,g_kappa,g_ka_csw_8, g_mu);
    fprintf(parameterfile,"boundary %f %f %f %f \n \n",X0,X1,X2,X3);
    fprintf(parameterfile,"ITER_MAX_BCG=%d, EPS_SQ0=%e, EPS_SQ1=%e EPS_SQ2=%e, EPS_SQ3=%e \n\n"
	    ,ITER_MAX_BCG,EPS_SQ0,EPS_SQ1,EPS_SQ2,EPS_SQ3);
    fprintf(parameterfile,"q_off=%f, q_off2=%f, dtau=%f, Nsteps=%d, Nmeas=%d, Nskip=%d, integtyp=%d, nsmall=%d \n\n",
	    q_off,q_off2,dtau,Nsteps,Nmeas,Nskip,integtyp,nsmall);
    fprintf(parameterfile,"mu = %f", g_mu);
    
  }

  /* define the geometry */
  geometry();

  /* define the boundary conditions for the fermion fields */
  boundary();

#if defined GEOMETRIC
  if(g_proc_id==0) fprintf(parameterfile,"The geometric series is used as solver \n\n");
#else
  if(g_proc_id==0) fprintf(parameterfile,"The BICG_stab is used as solver \n\n");
#endif
  
  /* Continue */
  if(startoption == 3){
    if (g_proc_id == 0){
      rlxdfile=fopen("rlxd_state","r");
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

  for(j=0;j<Nmeas;j++) {
    /* copy the gauge field to gauge_tmp */
    for(ix=0;ix<VOLUME;ix++) { 
      for(mu=0;mu<4;mu++) {
	v=&g_gauge_field[ix][mu];
	w=&gauge_tmp[ix][mu];
	_su3_assign(*w,*v);
      }
    }

    /* initialize the pseudo-fermion fields    */
    /* depending on q_off and q_off2 we use    */
    /* one, two or three pseudo-fermion fields */
    random_spinor_field(2);
    /* compute the square of the norm */
    enerphi0 = square_norm(2, VOLUME/2);

    /* apply the fermion matrix to the spinor */
    Qtm_plus_psi(0, 2);

    /* initialize the momenta */
    enep=ini_momenta();

    /*run the trajectory*/
    if(integtyp == 1) {
      /* Leap-frog integration scheme */
      leap_frog(q_off, q_off2, dtau, Nsteps, nsmall); 
    }
    else if(integtyp == 2) {
      /* Sexton Weingarten integration scheme */
      sexton(q_off, q_off2, dtau, Nsteps, nsmall);
    }

    /*perform the accept-reject-step*/
    enepx=moment_energy();
    enegx=measure_gauge_action();
    /*compute the energy contributions from the pseudo-fermions */

    zero_spinor_field(2);
    idis0=bicg(2, 0, q_off, EPS_SQ0);

    enerphi0x=square_norm(2, VOLUME/2);

    /* Compute the energy difference */
    dh=+enepx - g_beta*enegx - enep + g_beta*eneg
      + enerphi0x - enerphi0; 
      
    /* the random number is only taken at node zero and then distributed to 
       the other sites */
    if(g_proc_id==0) {
      ranlxd(yy,1);
#ifdef MPI
      for(i = 1; i < g_nproc; i++) {
	MPI_Send(&yy[0], 1, MPI_DOUBLE, i, 31, MPI_COMM_WORLD);
      }
#endif
    }
#ifdef MPI
    else{
      MPI_Recv(&yy[0], 1, MPI_DOUBLE, 0, 31, MPI_COMM_WORLD, &status);
    }
#endif
    if(exp(-dh) > yy[0]) {
      /* accept */
      Rate += 1;
      eneg=enegx;
      /* put the links back to SU(3) group */
      for(ix=0;ix<VOLUME;ix++) { 
	for(mu=0;mu<4;mu++) { 
	  /* this is MIST */
	  v=&g_gauge_field[ix][mu];
	  *v=restoresu3(*v); 
	}
      }
    }
    else {
      /* reject: copy gauge_tmp to g_gauge_field */
      for(ix=0;ix<VOLUME;ix++) {
	for(mu=0;mu<4;mu++){
	  /* Auch MIST */
	  v=&g_gauge_field[ix][mu];
	  w=&gauge_tmp[ix][mu];
	  _su3_assign(*v,*w);
	}
      }
    }
#ifdef MPI
    xchange_gauge();
#endif

    if(g_proc_id==0){
      fprintf(datafile,"%14.12f %14.12f %e %d %d %d \n",
	      eneg/(6.*VOLUME*g_nproc),dh,exp(-dh),
	      idis0,count00,count01);
      fflush(datafile);
    }
    /* Save gauge configuration all Nskip times */
    if((j+1)%Nskip == 0) {
      sprintf(filename3,"%s.%.4d", "conf", nstore);
      nstore ++;
      write_gauge_field_time_p( filename3 );

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
    fclose(parameterfile);
    fclose(datafile);
  }
#ifdef MPI
  MPI_Finalize();
#endif
  return(0);
}

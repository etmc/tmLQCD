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
#include "xchange.h"
#include "sw.h"
#include "io.h"
#include "read_input.h"
#include "tm_operators.h"
#include "bicgstabell.h"
#include "Hopping_Matrix.h"
#include "gamma.h"
#include "mpi_init.h"
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
 
  FILE *fp1=NULL,*fp2=NULL,*fp4=NULL;
  char * filename = NULL;
  char filename1[50];
  char filename2[50];
  char filename3[50];
  int rlxd_state[105];
  int idis0;
  int j,ix,mu;
  int i,k;
#ifdef _GAUGE_COPY
  int kb=0;
#endif

  double yy[1];
/*   static double step; */
/*   static double q_off,q_off2; */
  double eneg,enegx,enep,enepx,dh;
  double enerphi0,enerphi0x;
  static su3 gauge_tmp[VOLUME][4];
  su3 *v,*w;

  int Rate=0;
  int  namelen;
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  char * input_filename = NULL;
  int c;

  /* Test zeugs */
  int x0, x1, x2, x3, f, g, h;
  complex Trace[VOLUME][3];


  verbose = 0;
  g_use_clover_flag = 0;
 
  
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
 
  q_off = 0.;
  q_off2 = 0.;
 
  if(g_proc_id == 0){
    
/*     fscanf(fp6,"%s",filename); */
    /*construct the filenames for the observables and the parameters*/
    strcpy(filename1,filename);  strcat(filename1,".data");
    strcpy(filename2,filename);  strcat(filename2,".para");
    
    fp1=fopen(filename1, "w");
    fp2=fopen(filename2, "w");
    
    if(g_proc_id == 0){
      printf("# This is the hmc code for twisted Mass Wilson QCD\n\n Version %s", Version);
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
      printf("# The code was compiled for pentium4\n");
#endif
    }
    printf("The lattice size is %d x %d^3\n",(int)(T)*g_nproc,(int)(L));
    printf("The local lattice size is %d x %d^3\n",(int)(T),(int)(L));
    printf("g_beta = %f , g_kappa= %f, g_kappa*csw/8= %f \n\n",g_beta,g_kappa,g_ka_csw_8);
    
    fprintf(fp2,"The lattice size is %d x %d^3\n\n",(int)(g_nproc*T),(int)(L));
    fprintf(fp2,"g_beta = %f , g_kappa= %f, g_kappa*csw/8= %f g_mu = %f \n\n",g_beta,g_kappa,g_ka_csw_8, g_mu);
    fprintf(fp2,"boundary %f %f %f %f \n \n",X0,X1,X2,X3);
    fprintf(fp2,"ITER_MAX_BCG=%d, EPS_SQ0=%e, EPS_SQ1=%e EPS_SQ2=%e, EPS_SQ3=%e \n\n"
	    ,ITER_MAX_BCG,EPS_SQ0,EPS_SQ1,EPS_SQ2,EPS_SQ3);
    fprintf(fp2,"q_off=%f, q_off2=%f, dtau=%f, Nsteps=%d, Nmeas=%d, Nskip=%d, integtyp=%d, nsmall=%d \n\n",
	    q_off,q_off2,dtau,Nsteps,Nmeas,Nskip,integtyp,nsmall);
    fprintf(fp2,"mu = %f", g_mu);
    
  }

  /* define the geometry */
  geometry();

  /* define the boundary conditions for the fermion fields */
  boundary();

#if defined GEOMETRIC
  if(g_proc_id==0) fprintf(fp2,"The geometric series is used as solver \n\n");
#else
  if(g_proc_id==0) fprintf(fp2,"The BICG_stab is used as solver \n\n");
#endif
  
  /* Continue */
  if(startoption == 3){
    if (g_proc_id == 0){
      fp4=fopen("rlxd_state","r");
      fread(rlxd_state,sizeof(rlxd_state),1,fp4);
      fclose(fp4);
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
      MPI_Send(&rlxd_state[0], 105, MPI_INT, 1, 99, MPI_COMM_WORLD);
      MPI_Recv(&rlxd_state[0], 105, MPI_INT, g_nproc-1, 99, MPI_COMM_WORLD, &status);
      rlxd_reset(rlxd_state);
    }
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

  unit_g_gauge_field();   
  read_gauge_field_time_p("last_configuration");
#ifdef MPI
  xchange_gauge();
#endif
#ifdef _GAUGE_COPY
  /* set the backward gauge field */
  for(ix = 0; ix < VOLUME+RAND;ix++) {
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
/*   g_kappa*=5; */
  boundary();
  if(g_proc_id == 0){
    for(x0 = 0; x0 < 4; x0++){
      printf("%e %e\n%e %e\n%e %e\n%e %e\n%e %e\n%e %e\n%e %e\n%e %e\n%e %e\n\n", 
	     g_gauge_field[g_ipt[1][1][1][2]][x0].c00.re, 
	     g_gauge_field[g_ipt[1][1][1][2]][x0].c00.im,
	     g_gauge_field[g_ipt[1][1][1][2]][x0].c01.re, 
	     g_gauge_field[g_ipt[1][1][1][2]][x0].c01.im,
	     g_gauge_field[g_ipt[1][1][1][2]][x0].c02.re, 
	     g_gauge_field[g_ipt[1][1][1][2]][x0].c02.im,
	     g_gauge_field[g_ipt[1][1][1][2]][x0].c10.re, 
	     g_gauge_field[g_ipt[1][1][1][2]][x0].c10.im,
	     g_gauge_field[g_ipt[1][1][1][2]][x0].c11.re, 
	     g_gauge_field[g_ipt[1][1][1][2]][x0].c11.im,
	     g_gauge_field[g_ipt[1][1][1][2]][x0].c12.re, 
	     g_gauge_field[g_ipt[1][1][1][2]][x0].c12.im,
	     g_gauge_field[g_ipt[1][1][1][2]][x0].c20.re, 
	     g_gauge_field[g_ipt[1][1][1][2]][x0].c20.im,
	     g_gauge_field[g_ipt[1][1][1][2]][x0].c21.re, 
	     g_gauge_field[g_ipt[1][1][1][2]][x0].c21.im,
	     g_gauge_field[g_ipt[1][1][1][2]][x0].c22.re, 
	     g_gauge_field[g_ipt[1][1][1][2]][x0].c22.im
	     );
    }
  }
  eneg=measure_gauge_action();
  if(g_proc_id==0){
    printf("#First plaquette value: %14.12f \n",eneg/(6.*VOLUME*g_nproc));
  }
  for(x0 = 0; x0 < VOLUME; x0++){
    for(x1 = 0; x1 < 3; x1 ++){
      Trace[x0][x1].re = 0.;  
      Trace[x0][x1].im = 0.;  
    }
  }


  zero_spinor_field(2);
  zero_spinor_field(1);
  if(g_proc_id == 0){
    spinor_field[1][g_lexic2eo[g_ipt[1][2][1][2]]].s0.c0.re = 1.;
  }

  Hopping_Matrix(EO, 0, 2);
  mul_one_pm_imu(1, 1.);
  assign_add_mul_r(spinor_field[1], spinor_field[0], -1., VOLUME/2);
/*   gamma5(1, 1);     */
  gamma3(4, 1, VOLUME/2);  

  zero_spinor_field(3);
  zero_spinor_field(2);
  if(g_proc_id == 0) {
    spinor_field[2][g_lexic2eo[g_ipt[1][2][1][2]]].s0.c0.re = 1.; 
  }
  Hopping_Matrix(OE, 0, 2);
  mul_one_pm_imu(3, 1.);
  assign_add_mul_r(spinor_field[3], spinor_field[0], -1., VOLUME/2);
/*   gamma5(3, 3);   */
  gamma3(5, 3, VOLUME/2);  
  if(g_proc_id == 0 || g_proc_id == 1){
    for (x0=0;x0<T;x0++){
      for (x1=0;x1<L;x1++){
	for (x2=0;x2<L;x2++){
	  for (x3=0;x3<L;x3++){
	    if((x0+x1+x2+x3+g_proc_id*T)%2==0) {
	      f = 4;
	      g = g_lexic2eo[ g_ipt[x0][x1][x2][x3] ];
	      h = g;
	    } 
	    else {
	      f = 5;
	      h = g_lexic2eo[ g_ipt[x0][x1][x2][x3] ] - RAND/2 ;
	      g = h - VOLUME/2;
	    }
	    
	    Trace[h][0].re += spinor_field[f][g].s0.c0.re;
	    Trace[h][1].re += spinor_field[f][g].s0.c1.re;
	    Trace[h][2].re += spinor_field[f][g].s0.c2.re;
	    Trace[h][0].im += spinor_field[f][g].s0.c0.im;
	    Trace[h][1].im += spinor_field[f][g].s0.c1.im;
	    Trace[h][2].im += spinor_field[f][g].s0.c2.im;
	  }
	}
      }
    }
    
    zero_spinor_field(2);
    zero_spinor_field(1);
    if(g_proc_id == 0){
      spinor_field[1][g_lexic2eo[g_ipt[1][2][1][2]]].s1.c0.re = 1.;
    }
    
    Hopping_Matrix(EO, 0, 2);
    mul_one_pm_imu(1, 1.);
    assign_add_mul_r(spinor_field[1], spinor_field[0], -1., VOLUME/2);
    /*   gamma5(1, 1);   */
    gamma3(4, 1, VOLUME/2);
    
    zero_spinor_field(3);
    zero_spinor_field(2);
    if(g_proc_id == 0){
      spinor_field[2][g_lexic2eo[g_ipt[1][2][1][2]]].s1.c0.re = 1.; 
    }
    
    Hopping_Matrix(OE, 0, 2);
    mul_one_pm_imu(3, 1.);
    assign_add_mul_r(spinor_field[3], spinor_field[0], -1., VOLUME/2);
/*     gamma5(3, 3);   */
    gamma3(5, 3, VOLUME/2);
    
    for (x0=0;x0<T;x0++){
      for (x1=0;x1<L;x1++){
	for (x2=0;x2<L;x2++){
	  for (x3=0;x3<L;x3++){
	    if((x0+x1+x2+x3+g_proc_id*T)%2==0) {
	      f = 4;
	      g = g_lexic2eo[ g_ipt[x0][x1][x2][x3] ];
	      h = g;
	    } 
	    else {
	      f = 5;
	      h = g_lexic2eo[ g_ipt[x0][x1][x2][x3] ] - RAND/2 ;
	      g = h - VOLUME/2;
	    }
	    
	    Trace[h][0].re += spinor_field[f][g].s1.c0.re;
	    Trace[h][1].re += spinor_field[f][g].s1.c1.re;
	    Trace[h][2].re += spinor_field[f][g].s1.c2.re;
	    Trace[h][0].im += spinor_field[f][g].s1.c0.im;
	    Trace[h][1].im += spinor_field[f][g].s1.c1.im;
	    Trace[h][2].im += spinor_field[f][g].s1.c2.im;
	  }
	}
      }
    }
    
    zero_spinor_field(2);
    zero_spinor_field(1);
    if(g_proc_id == 0){
      spinor_field[1][g_lexic2eo[g_ipt[1][2][1][2]]].s2.c0.re = 1.;
    }
    
    Hopping_Matrix(EO, 0, 2);
    mul_one_pm_imu(1, 1.);
    assign_add_mul_r(spinor_field[1], spinor_field[0], -1., VOLUME/2);
    /*   gamma5(1, 1);   */
    gamma3(4, 1, VOLUME/2);
    
    zero_spinor_field(3);
    zero_spinor_field(2);
    if(g_proc_id == 0){
      spinor_field[2][g_lexic2eo[g_ipt[1][2][1][2]]].s2.c0.re = 1.; 
    }
    
    Hopping_Matrix(OE, 0, 2);
    mul_one_pm_imu(3, 1.);
    assign_add_mul_r(spinor_field[3], spinor_field[0], -1., VOLUME/2);
/*     gamma5(3, 3);   */
    gamma3(5, 3, VOLUME/2);
    
    for (x0=0;x0<T;x0++){
      for (x1=0;x1<L;x1++){
	for (x2=0;x2<L;x2++){
	  for (x3=0;x3<L;x3++){
	    if((x0+x1+x2+x3+g_proc_id*T)%2==0) {
	      f = 4;
	      g = g_lexic2eo[ g_ipt[x0][x1][x2][x3] ];
	      h = g;
	    } 
	    else {
	      f = 5;
	      h = g_lexic2eo[ g_ipt[x0][x1][x2][x3] ] - RAND/2 ;
	      g = h - VOLUME/2;
	    }
	    
	    Trace[h][0].re += spinor_field[f][g].s2.c0.re;
	    Trace[h][1].re += spinor_field[f][g].s2.c1.re;
	    Trace[h][2].re += spinor_field[f][g].s2.c2.re;
	    Trace[h][0].im += spinor_field[f][g].s2.c0.im;
	    Trace[h][1].im += spinor_field[f][g].s2.c1.im;
	    Trace[h][2].im += spinor_field[f][g].s2.c2.im;
	  }
	}
      }
    }
    
    zero_spinor_field(2);
    zero_spinor_field(1);
    if(g_proc_id == 0){
      spinor_field[1][g_lexic2eo[g_ipt[1][2][1][2]]].s3.c0.re = 1.;
    }
    
    Hopping_Matrix(EO, 0, 2);
    mul_one_pm_imu(1, 1.);
    assign_add_mul_r(spinor_field[1], spinor_field[0], -1., VOLUME/2);
    /*   gamma5(1, 1);   */
    gamma3(4, 1, VOLUME/2);
    
    zero_spinor_field(3);
    zero_spinor_field(2);
    if(g_proc_id == 0){
      spinor_field[2][g_lexic2eo[g_ipt[1][2][1][2]]].s3.c0.re = 1.; 
    }
    
    Hopping_Matrix(OE, 0, 2);
    mul_one_pm_imu(3, 1.);
    assign_add_mul_r(spinor_field[3], spinor_field[0], -1., VOLUME/2);
/*     gamma5(3, 3);   */
    gamma3(5, 3, VOLUME/2);
    
    for (x0=0;x0<T;x0++){
      for (x1=0;x1<L;x1++){
	for (x2=0;x2<L;x2++){
	  for (x3=0;x3<L;x3++){
	    if((x0+x1+x2+x3+g_proc_id*T)%2==0) {
	      f = 4;
	      g = g_lexic2eo[ g_ipt[x0][x1][x2][x3] ];
	      h = g;
	    } 
	    else {
	      f = 5;
	      h = g_lexic2eo[ g_ipt[x0][x1][x2][x3] ] - RAND/2 ;
	      g = h - VOLUME/2;
	    }
	    
	    Trace[h][0].re += spinor_field[f][g].s3.c0.re;
	    Trace[h][1].re += spinor_field[f][g].s3.c1.re;
	    Trace[h][2].re += spinor_field[f][g].s3.c2.re;
	    Trace[h][0].im += spinor_field[f][g].s3.c0.im;
	    Trace[h][1].im += spinor_field[f][g].s3.c1.im;
	    Trace[h][2].im += spinor_field[f][g].s3.c2.im;
	  }
	}
      }
    }


    
    if(g_proc_id != 5){
      for (x0=0;x0<T;x0++){
	for (x1=0;x1<L;x1++){
	  for (x2=0;x2<L;x2++){
	    for (x3=0;x3<L;x3++){
	      if((x0+x1+x2+x3+g_proc_id*T)%2==0) {
		g = g_lexic2eo[ g_ipt[x0][x1][x2][x3] ];
		h = g;
	      } 
	      else {
		h = g_lexic2eo[ g_ipt[x0][x1][x2][x3] ] - RAND/2 ;
		g = h - VOLUME/2;
	      }
	      printf("(x,y,z,t,c) (%d, %d, %d, %d, 0) %e %e\n", x1, x2, x3, x0+g_proc_id*T, Trace[h][0].re, Trace[h][0].im);
	      printf("(x,y,z,t,c) (%d, %d, %d, %d, 1) %e %e\n", x1, x2, x3, x0+g_proc_id*T, Trace[h][1].re, Trace[h][1].im);
	      printf("(x,y,z,t,c) (%d, %d, %d, %d, 2) %e %e\n", x1, x2, x3, x0+g_proc_id*T, Trace[h][2].re, Trace[h][2].im);
	    }
	  }
	}
      }
    }
    
    MPI_Finalize();
    exit(0);
  }
  
  if(g_proc_id == 0) {
    for (x0=0;x0<T;x0++){
      for (x1=0;x1<L;x1++){
	for (x2=0;x2<L;x2++){
	  for (x3=0;x3<L;x3++){
	    if((x0+x1+x2+x3+g_proc_id*T)%2==0) {
	      f = 1;
	      g = g_lexic2eo[ g_ipt[x0][x1][x2][x3] ];
	    } 
	    else {
	      f = 3;
	      g = g_lexic2eo[ g_ipt[x0][x1][x2][x3] ] - (VOLUME+RAND)/2;
	    }
	    printf("(t,x,y,z,s,c) (%d, %d, %d, %d, 0 0) %e %e\n", x0+g_proc_id*2, x1, x2, x3, spinor_field[f][g].s0.c0.re, spinor_field[f][g].s0.c0.im);
	    printf("(t,x,y,z,s,c) (%d, %d, %d, %d, 0 1) %e %e\n", x0+g_proc_id*2, x1, x2, x3, spinor_field[f][g].s0.c1.re, spinor_field[f][g].s0.c1.im);
	    printf("(t,x,y,z,s,c) (%d, %d, %d, %d, 0 2) %e %e\n", x0+g_proc_id*2, x1, x2, x3, spinor_field[f][g].s0.c2.re, spinor_field[f][g].s0.c2.im);
	    printf("(t,x,y,z,s,c) (%d, %d, %d, %d, 1 0) %e %e\n", x0+g_proc_id*2, x1, x2, x3, spinor_field[f][g].s1.c0.re, spinor_field[f][g].s1.c0.im);
	    printf("(t,x,y,z,s,c) (%d, %d, %d, %d, 1 1) %e %e\n", x0+g_proc_id*2, x1, x2, x3, spinor_field[f][g].s1.c1.re, spinor_field[f][g].s1.c1.im);
	    printf("(t,x,y,z,s,c) (%d, %d, %d, %d, 1 2) %e %e\n", x0+g_proc_id*2, x1, x2, x3, spinor_field[f][g].s1.c2.re, spinor_field[f][g].s1.c2.im);
	    printf("(t,x,y,z,s,c) (%d, %d, %d, %d, 2 0) %e %e\n", x0+g_proc_id*2, x1, x2, x3, spinor_field[f][g].s2.c0.re, spinor_field[f][g].s2.c0.im);
	    printf("(t,x,y,z,s,c) (%d, %d, %d, %d, 2 1) %e %e\n", x0+g_proc_id*2, x1, x2, x3, spinor_field[f][g].s2.c1.re, spinor_field[f][g].s2.c1.im);
	    printf("(t,x,y,z,s,c) (%d, %d, %d, %d, 2 2) %e %e\n", x0+g_proc_id*2, x1, x2, x3, spinor_field[f][g].s2.c2.re, spinor_field[f][g].s2.c2.im);
	    printf("(t,x,y,z,s,c) (%d, %d, %d, %d, 3 0) %e %e\n", x0+g_proc_id*2, x1, x2, x3, spinor_field[f][g].s3.c0.re, spinor_field[f][g].s3.c0.im);
	    printf("(t,x,y,z,s,c) (%d, %d, %d, %d, 3 1) %e %e\n", x0+g_proc_id*2, x1, x2, x3, spinor_field[f][g].s3.c1.re, spinor_field[f][g].s3.c1.im);
	    printf("(t,x,y,z,s,c) (%d, %d, %d, %d, 3 2) %e %e\n", x0+g_proc_id*2, x1, x2, x3, spinor_field[1][g].s3.c2.re, spinor_field[f][g].s3.c2.im);
	    printf("\n");
	    fflush(stdout);
	  }
	}
      }
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);

  if(g_proc_id == 3) {
    for (x0=0;x0<T;x0++){
      for (x1=0;x1<L;x1++){
	for (x2=0;x2<L;x2++){
	  for (x3=0;x3<L;x3++){
	    if((x0+x1+x2+x3+g_proc_id*T)%2==0){
	      f = 1;
	      g = g_lexic2eo[ g_ipt[x0][x1][x2][x3] ];
	    } 
	    else{
	      f = 3;
	      g = g_lexic2eo[ g_ipt[x0][x1][x2][x3] ] - (VOLUME+RAND)/2;
	    }
	    printf("(t,x,y,z,s,c) (%d, %d, %d, %d, 0 0) %e %e\n", x0+2, x1, x2, x3, spinor_field[f][g].s0.c0.re, spinor_field[f][g].s0.c0.im);
	    printf("(t,x,y,z,s,c) (%d, %d, %d, %d, 0 1) %e %e\n", x0+2, x1, x2, x3, spinor_field[f][g].s0.c1.re, spinor_field[f][g].s0.c1.im);
	    printf("(t,x,y,z,s,c) (%d, %d, %d, %d, 0 2) %e %e\n", x0+2, x1, x2, x3, spinor_field[f][g].s0.c2.re, spinor_field[f][g].s0.c2.im);
	    printf("(t,x,y,z,s,c) (%d, %d, %d, %d, 1 0) %e %e\n", x0+2, x1, x2, x3, spinor_field[f][g].s1.c0.re, spinor_field[f][g].s1.c0.im);
	    printf("(t,x,y,z,s,c) (%d, %d, %d, %d, 1 1) %e %e\n", x0+2, x1, x2, x3, spinor_field[f][g].s1.c1.re, spinor_field[f][g].s1.c1.im);
	    printf("(t,x,y,z,s,c) (%d, %d, %d, %d, 1 2) %e %e\n", x0+2, x1, x2, x3, spinor_field[f][g].s1.c2.re, spinor_field[f][g].s1.c2.im);
	    printf("(t,x,y,z,s,c) (%d, %d, %d, %d, 2 0) %e %e\n", x0+2, x1, x2, x3, spinor_field[f][g].s2.c0.re, spinor_field[f][g].s2.c0.im);
	    printf("(t,x,y,z,s,c) (%d, %d, %d, %d, 2 1) %e %e\n", x0+2, x1, x2, x3, spinor_field[f][g].s2.c1.re, spinor_field[f][g].s2.c1.im);
	    printf("(t,x,y,z,s,c) (%d, %d, %d, %d, 2 2) %e %e\n", x0+2, x1, x2, x3, spinor_field[f][g].s2.c2.re, spinor_field[f][g].s2.c2.im);
	    printf("(t,x,y,z,s,c) (%d, %d, %d, %d, 3 0) %e %e\n", x0+2, x1, x2, x3, spinor_field[f][g].s3.c0.re, spinor_field[f][g].s3.c0.im);
	    printf("(t,x,y,z,s,c) (%d, %d, %d, %d, 3 1) %e %e\n", x0+2, x1, x2, x3, spinor_field[f][g].s3.c1.re, spinor_field[f][g].s3.c1.im);
	    printf("(t,x,y,z,s,c) (%d, %d, %d, %d, 3 2) %e %e\n", x0+2, x1, x2, x3, spinor_field[1][g].s3.c2.re, spinor_field[f][g].s3.c2.im);
	    printf("\n");
	  }
	}
      }
    }
  }


  MPI_Finalize();
  return(0);
}

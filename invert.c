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
#include "solver/solver.h"
#include "invert_eo.h"

char * Version = "2.1";


void usage(){
  fprintf(stderr, "Inversion for EO preconditioned Wilson twisted mass QCD\n\n");
  fprintf(stderr, "Usage: [-f input-filename]\n");
  fprintf(stderr, "Usage: [-o output-filename]\n");
  exit(1);
}

extern int nstore;

int check_geometry();

int main(int argc,char *argv[]) {

  FILE *parameterfile=NULL;
  int c, iter, j, ix=0, is=0, ic=0;
  char * filename = NULL;
  char datafilename[50];
  char parameterfilename[50];
  char gauge_filename[50];
  char conf_filename[50];
  char * input_filename = NULL;
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
    printf("# The lattice size is %d x %d^3\n",(int)(T)*g_nproc_t,(int)(L));
    printf("# The local lattice size is %d x %d^3\n",(int)(T),(int)(L));
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
    read_gauge_field_time_p(conf_filename);
if (g_proc_id == 0){
      printf("done!\n"); fflush(stdout);
    }
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

    sprintf(conf_filename,"%s.%.4d", "prop.mass00", nstore);
    
    for(ix = index_start; ix < index_end; ix++) {
      is = (ix / 3);
      ic = (ix % 3);
      source_spinor_field(spinor_field[0], spinor_field[1], is, ic);

      printf("mu = %e\n", g_mu);

#ifdef MPI
      atime = MPI_Wtime();
#endif
      iter = invert_eo(spinor_field[2], spinor_field[3], spinor_field[0], spinor_field[1], 
		       1.e-15, ITER_MAX_CG, solver_flag);
#ifdef MPI
      etime = MPI_Wtime();
#endif
      mul_r(spinor_field[2], (2*g_kappa), spinor_field[2], VOLUME/2);  
      mul_r(spinor_field[3], (2*g_kappa), spinor_field[3], VOLUME/2);
      if(ix == 0) {
	write_spinorfield_eo_time_p(spinor_field[2], spinor_field[3], conf_filename, 0);
      }
      else {
	write_spinorfield_eo_time_p(spinor_field[2], spinor_field[3], conf_filename, 1);
      }
    
      M_full(spinor_field[4], spinor_field[5], spinor_field[2], spinor_field[3]); 
      mul_r(spinor_field[4], 1./(2*g_kappa), spinor_field[4], VOLUME/2);  
      mul_r(spinor_field[5], 1./(2*g_kappa), spinor_field[5], VOLUME/2); 
      diff(spinor_field[4], spinor_field[4], spinor_field[0], VOLUME/2); 
      diff(spinor_field[5], spinor_field[5], spinor_field[1], VOLUME/2); 

      nrm1 = square_norm(spinor_field[4], VOLUME/2); 
      nrm2 = square_norm(spinor_field[5], VOLUME/2); 

      printf("Inversion for is = %d, ic = %d done in %d iterations, residue = %e!\n", is, ic, iter, nrm1+nrm2);
      printf("Inversion done in %e sec.\n", etime-atime);
    }
    nstore++;
  }

  return(0);
}

/* $Id$ */
/*******************************************************************************
*
*
* source generation main file 
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
#ifdef MPI
# include <mpi.h>
#endif
#include "global.h"
#include "getopt.h"
#include "geometry_eo.h"
#include "start.h"
#include "io.h"
#include "read_input.h"
#include "mpi_init.h"
#include "source_generation.h"
#include "init_geometry_indices.h"
#include "linalg_eo.h"
#include "init_spinor_field.h"


void usage() {
  fprintf(stdout, "Code to generate stochastic sources\n");
  fprintf(stdout, "Version %s \n\n", PACKAGE_VERSION);
  fprintf(stdout, "Please send bug reports to %s\n", PACKAGE_BUGREPORT);
  fprintf(stdout, "Usage:   gen_sources [options]\n");
  fprintf(stdout, "Options: -L spatial lattice size\n");
  fprintf(stdout, "         -T temporal lattice size\n");
  fprintf(stdout, "         -o output-filename basename [optional, default source]\n");
  fprintf(stdout, "         -n configuration number [optional, default 0]\n");
  fprintf(stdout, "         -s sample number [optional, default 0]\n");
  fprintf(stdout, "         -t start timslice [optional, default 0]\n");
  fprintf(stdout, "         -S spatial spacing [optional, default 1]\n");
  fprintf(stdout, "         -P temporal spacing [optional, default T]\n");
  fprintf(stdout, "         -N produce nucleon sources [optional, default meson]\n");
  fprintf(stdout, "         -p use plain output filename [default, complex]\n");
  fprintf(stdout, "         -O pion only -> wallsource at start timeslice \n");
  fprintf(stdout, "         -E extendent source for pion only \n");
  fprintf(stdout, "         -h|-? this help \n\n");
  fprintf(stdout, "plain output file (-p) corresponds to basename.00 - basename.11\n");
  fprintf(stdout, "complex ones (no -p) to basename.samplenr.gaugenr.tsnr.00 - 11\n");
  exit(0);
}


extern int nstore;
const int rlxdsize = 105;

int main(int argc,char *argv[]) {
 
  char spinorfilename[100];
  char * filename = NULL;
  int sample=0, ts=0, ss=1, typeflag = 1, t0=0, formatflag = 0, piononly, ext_sourceflag = 0;
  int is, ic, j, tt=0, filenameflag = 0;
  complex co;
  int c;

  verbose = 0;
  g_use_clover_flag = 0;
  g_nr_of_psf = 1;
  nstore = 0;
  L=0;
  T=0;
  
#ifdef MPI
  MPI_Init(&argc, &argv);
#endif


  while ((c = getopt(argc, argv, "h?NCpOEo:L:T:n:t:s:S:P:")) != -1) {
    switch (c) {
    case 'L':
      L = atoi(optarg);
      LX = L;
      LY = L;
      LZ = L;
      break;
    case 'T':
      T = atoi(optarg);
      T_global = T;
      break;
    case 'N':
      typeflag = 0;
      break;
    case 'O':
      piononly = 1;
      break;
    case 'n':
      nstore = atoi(optarg);
      break;
    case 's':
      sample = atoi(optarg);
      break;
    case 't':
      t0 = atoi(optarg);
      break;
    case 'S':
      ss = atoi(optarg);
      break;
    case 'P':
      ts = atoi(optarg);
      break;
    case 'C':
      formatflag = 1;
      break;
    case 'o':
      filename = calloc(200, sizeof(char));
      strcpy(filename,optarg);
      break;
    case 'E':
      ext_sourceflag = 1;
      break;
    case 'p':
      filenameflag = 1;
      break;
    case 'h':
    case '?':
    default:
      usage();
      break;
    }
  }
  if(ts == 0) {
    ts = T;
  }
  if(filename == NULL){
    filename = "source";
  } 
  if(L==0 || T==0) {
    if(g_proc_id == 0) {
      fprintf(stderr, "L and T must be specified! Aborting...\n");
      fflush( stderr );
    }
    exit(1);
  }

  mpi_init(argc, argv);

  j = init_geometry_indices(VOLUMEPLUSRAND);
  if ( j!= 0) {
    fprintf(stderr, "Not enough memory for geometry_indices! Aborting...\n");
    exit(0);
  }
  if(!ext_sourceflag) {
    j = init_spinor_field(VOLUMEPLUSRAND/2, 2);
  }
  else {
    j = init_spinor_field(VOLUMEPLUSRAND/2, 4);
  }
  if ( j!= 0) {
    fprintf(stderr, "Not enough memory for spinor fields! Aborting...\n");
    exit(0);
  }

  /* define the geometry */
  geometry();
  
  if(!piononly) {
    for(is = 0; is < 4; is ++) {
      for(ic = 0; ic < 3; ic++) {
	if(!filenameflag) {
	  sprintf(spinorfilename, "%s.%.4d.%.4d.%.2d.%.2d", filename, nstore, sample, tt, 3*is+ic); 
	}
	else {
	  sprintf(spinorfilename, "%s.%.2d", filename, 3*is+ic); 
	}
	printf("Generating source %s!\n", spinorfilename);
	fflush(stdout);
	
	source_generation_nucleon(g_spinor_field[0], g_spinor_field[1], 
				  is, ic, t0, ts, ss, sample, nstore, typeflag);
	
	co = scalar_prod(g_spinor_field[1], g_spinor_field[1], VOLUME/2);
	if(formatflag == 1) {
	  write_spinorfield_cm_single(g_spinor_field[0], g_spinor_field[1], spinorfilename);
	}
	else {
	  write_spinorfield_eo_time_p(g_spinor_field[0], g_spinor_field[1], spinorfilename, 0);
	}
      }
    }
  }
  else {
    if(!ext_sourceflag) {
      if(!filenameflag) {
	sprintf(spinorfilename, "%s.%.4d.%.4d.%.2d", filename, nstore, sample, tt); 
      }
      else {
	sprintf(spinorfilename, "%s", filename); 
      }
      printf("Generating source %s!\n", spinorfilename);
      fflush(stdout);
      source_generation_pion_only(g_spinor_field[0], g_spinor_field[1], 
				  t0, sample, nstore);
      
      co = scalar_prod(g_spinor_field[1], g_spinor_field[1], VOLUME/2);
      if(formatflag == 1) {
	write_spinorfield_cm_single(g_spinor_field[0], g_spinor_field[1], spinorfilename);
      }
      else {
	write_spinorfield_eo_time_p(g_spinor_field[0], g_spinor_field[1], spinorfilename, 0);
      }
    }
    else {
      if(!filenameflag) {
        sprintf(spinorfilename, "%s.%.4d.%.4d.%.2d.inverted", filename, nstore, sample, tt);
      }
      else {
        sprintf(spinorfilename, "%s.inverted", filename);
      }
      if(formatflag == 1) {
        read_spinorfield_cm_single(g_spinor_field[0], g_spinor_field[1], spinorfilename, -1, 1);
      }
      else {
        read_spinorfield_eo_time(g_spinor_field[0], g_spinor_field[1], spinorfilename);
      }
      extended_pion_source(g_spinor_field[2], g_spinor_field[3],
			   g_spinor_field[0], g_spinor_field[1],
			   (t0+T/2)%T, 0., 0., 0.);
      if(!filenameflag) {
	sprintf(spinorfilename, "g%s.%.4d.%.4d.%.2d", filename, nstore, sample, tt); 
      }
      else {
	sprintf(spinorfilename, "g%s", filename); 
      }
      if(formatflag == 1) {
        write_spinorfield_cm_single(g_spinor_field[2], g_spinor_field[3], spinorfilename);
      }
      else {
        write_spinorfield_eo_time_p(g_spinor_field[2], g_spinor_field[3], spinorfilename, 0);
      }

    }
  }

#ifdef MPI
  MPI_Finalize();
#endif
  free_geometry_indices();
  free_spinor_field();
  return(0);
}

static char const rcsid[] = "$Id$";

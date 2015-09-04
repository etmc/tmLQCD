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
/*******************************************************************************
*
*
* source generation main file 
*
* Author: Carsten Urbach
*         urbach@physik.fu-berlin.de
*
*******************************************************************************/

#include "lime.h"
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
#ifdef OMP
# include <omp.h>
#endif
#include "global.h"
#include "getopt.h"
#include "geometry_eo.h"
#include "start.h"
#include <io/utils.h>
#include "read_input.h"
#include "mpi_init.h"
#include "source_generation.h"
#include "init/init.h"
#include "linalg_eo.h"
#include "phmc.h"

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
  fprintf(stdout, "         -E extended source for pion only \n");
  fprintf(stdout, "         -d double precision \n");
  fprintf(stdout, "         -a store all sources in one file\n");
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
  int sample=0, ts=0, ss=1, typeflag = 1, t0=0, piononly = 0, ext_sourceflag = 0;
  int is, ic, j, filenameflag = 0, appendflag = 0;
  complex co;
  int c;
  int prec=32;

  verbose = 0;
  g_use_clover_flag = 0;
  nstore = 0;
  L=0;
  T=0;
  
#ifdef MPI
  MPI_Init(&argc, &argv);
#endif

#ifdef OMP
  /* FIXME: in principle this should not be set like this as it could result
    in thread oversubscription when more than one process is run locally
    unfortunately, there does not seem to be a standard way to determine
    the number of "local" MPI processes  */
  omp_num_threads = omp_get_max_threads();
  init_openmp();
#endif

  while ((c = getopt(argc, argv, "h?NCpOEdao:L:T:n:t:s:S:P:")) != -1) {
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
    case 'd':
      prec = 64;
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
    case 'a':
      appendflag = 1;
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

  tmlqcd_mpi_init(argc, argv);

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
	if(!filenameflag && !appendflag) {
	  sprintf(spinorfilename, "%s.%.4d.%.4d.%.2d.%.2d", filename, nstore, sample, t0, 3*is+ic); 
	}
	else if(!filenameflag && appendflag) {
	  sprintf(spinorfilename, "%s.%.4d.%.4d.%.2d", filename, nstore, sample, t0); 
	}
	else{
	  sprintf(spinorfilename, "%s.%.2d", filename, 3*is+ic); 
	}
	if(!appendflag || (is == 0 && ic ==0)) {
	  printf("Generating source %s!\n", spinorfilename);
	  fflush(stdout);
	}
	
	source_generation_nucleon(g_spinor_field[0], g_spinor_field[1], 
				  is, ic, t0, ts, ss, sample, nstore, typeflag);
	
	co = scalar_prod(g_spinor_field[1], g_spinor_field[1], VOLUME/2, 1);
	if((is == 0 && ic == 0) || appendflag == 0) {
	  write_source_type(0, spinorfilename);
	}
	write_source(g_spinor_field[0], g_spinor_field[1], spinorfilename, 1, prec);
      }
    }
  }
  else {
    if(!ext_sourceflag) {
      if(!filenameflag) {
	sprintf(spinorfilename, "%s.%.4d.%.4d.%.2d", filename, nstore, sample, t0); 
      }
      else {
	sprintf(spinorfilename, "%s", filename); 
      }
      printf("Generating source %s!\n", spinorfilename);
      fflush(stdout);
      source_generation_pion_only(g_spinor_field[0], g_spinor_field[1], 
				  t0, sample, nstore);
      
      co = scalar_prod(g_spinor_field[1], g_spinor_field[1], VOLUME/2, 1);
      write_source_type(0, spinorfilename);
      write_source(g_spinor_field[0], g_spinor_field[1], spinorfilename, 1, prec);
    }
    else {
      if(!filenameflag) {
        sprintf(spinorfilename, "%s.%.4d.%.4d.%.2d.inverted", filename, nstore, sample, t0);
      }
      else {
        sprintf(spinorfilename, "%s.inverted", filename);
      }
      read_lime_spinor(g_spinor_field[0], g_spinor_field[1], spinorfilename, 0);

      printf("Generating ext. pion source %s!\n", spinorfilename);
      extended_pion_source(g_spinor_field[2], g_spinor_field[3],
			   g_spinor_field[0], g_spinor_field[1],
			   t0, 0., 0., 0.);
      if(!filenameflag) {
	sprintf(spinorfilename, "g%s.%.4d.%.4d.%.2d", filename, nstore, sample, t0); 
      }
      else {
	sprintf(spinorfilename, "g%s", filename); 
      }
      write_source_type(0, spinorfilename);
      write_source(g_spinor_field[2], g_spinor_field[3], spinorfilename, 1, prec);
    }
  }

#ifdef MPI
  MPI_Finalize();
#endif
  free_geometry_indices();
  free_spinor_field();
  return(0);
}


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
/* 
 *  Program for computing the eigensystem of the Laplacian operator
 * Authors Luigi Scorzato, Marco Cristoforetti
 *
 *
 *******************************************************************************/

#ifdef HAVE_CONFIG_H
# include "config.h"
#else
#error "no config.h"
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#if (defined BGL && !defined BGP)
#  include <rts.h>
#endif
#ifdef MPI
# include <mpi.h>
#endif
#include "global.h"
#include <io/params.h>
#include <io/gauge.h>
#include "su3.h"
#include "ranlxd.h"
#include "geometry_eo.h"
#include "read_input.h"
#include "start.h"
#include "xchange/xchange.h"
#include "init/init.h"
#include "mpi_init.h"
#include "solver/eigenvalues_Jacobi.h"

int main(int argc,char *argv[])
{
  int tslice,j,k;
  char conf_filename[50];
  
#ifdef MPI
  MPI_Init(&argc, &argv);
#endif
  
  /* Read the input file */
  read_input("LapH.input");
  
  tmlqcd_mpi_init(argc, argv);
  
  if(g_proc_id==0) {
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
#ifdef _GAUGE_COPY
    printf("# The code was compiled with -D_GAUGE_COPY\n");
#endif
#ifdef BGL
    printf("# The code was compiled for Blue Gene/L\n");
#endif
#ifdef BGP
    printf("# The code was compiled for Blue Gene/P\n");
#endif
#ifdef _USE_HALFSPINOR
    printf("# The code was compiled with -D_USE_HALFSPINOR\n");
#endif    
#ifdef _USE_SHMEM
    printf("# the code was compiled with -D_USE_SHMEM\n");
#  ifdef _PERSISTENT
    printf("# the code was compiled for persistent MPI calls (halfspinor only)\n");
#  endif
#endif
#ifdef MPI
#  ifdef _NON_BLOCKING
    printf("# the code was compiled for non-blocking MPI calls (spinor and gauge)\n");
#  endif
#endif
    printf("\n");
    fflush(stdout);
  }
  

#ifndef WITHLAPH
  printf(" Error: WITHLAPH not defined");
  exit(0);
#endif
#ifdef MPI
#ifndef _INDEX_INDEP_GEOM
  printf(" Error: _INDEX_INDEP_GEOM not defined");
  exit(0);
#endif
#ifndef _USE_TSPLITPAR
  printf(" Error: _USE_TSPLITPAR not defined");
  exit(0);
#endif
#endif
#ifdef FIXEDVOLUME
  printf(" Error: FIXEDVOLUME not allowed");
  exit(0);
#endif

  
  init_gauge_field(VOLUMEPLUSRAND + g_dbw2rand, 0);
  init_geometry_indices(VOLUMEPLUSRAND + g_dbw2rand);

  if(g_proc_id == 0) {
    fprintf(stdout,"The number of processes is %d \n",g_nproc);
    printf("# The lattice size is %d x %d x %d x %d\n",
	   (int)(T*g_nproc_t), (int)(LX*g_nproc_x), (int)(LY*g_nproc_y), (int)(g_nproc_z*LZ));
    printf("# The local lattice size is %d x %d x %d x %d\n", 
	   (int)(T), (int)(LX), (int)(LY),(int) LZ);
    printf("# Computing LapH eigensystem \n");

    fflush(stdout);
  }
  
  /* define the geometry */
  geometry();

  start_ranlux(1, 123456);

  /* Read Gauge field */
  sprintf(conf_filename, "%s.%.4d", gauge_input_filename, nstore);
  if (g_cart_id == 0) {
    printf("#\n# Trying to read gauge field from file %s in %s precision.\n",
	   conf_filename, (gauge_precision_read_flag == 32 ? "single" : "double"));
    fflush(stdout);
  }
  if( (j = read_gauge_field(conf_filename,g_gauge_field)) !=0) {
    fprintf(stderr, "Error %d while reading gauge field from %s\n Aborting...\n", j, conf_filename);
    exit(-2);
  }

  
  if (g_cart_id == 0) {
    printf("# Finished reading gauge field.\n");
    fflush(stdout);
  }
  
#ifdef MPI
  /*For parallelization: exchange the gaugefield */
  xchange_gauge(g_gauge_field);
#endif
  
  /* Init Jacobi field */
  init_jacobi_field(SPACEVOLUME+SPACERAND,3);

#ifdef MPI
  {
     /* for debugging in parallel set i_gdb = 0 */
    volatile int i_gdb = 8;
    char hostname[256];
    gethostname(hostname, sizeof(hostname));
    printf("PID %d on %s ready for attach\n", getpid(), hostname);
    fflush(stdout);
    if(g_cart_id == 0){
      while (0 == i_gdb){
	sleep(5);
      }
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
#endif

  for (k=0 ; k<3 ; k++)
    random_jacobi_field(g_jacobi_field[k],SPACEVOLUME);


  /* Compute LapH Eigensystem */
  
  for(tslice=0; tslice<T; tslice++){ 
    eigenvalues_Jacobi(&no_eigenvalues,5000, eigenvalue_precision,0,tslice,nstore);
  }
  
#ifdef MPI
  MPI_Finalize();
#endif
  return(0);
}

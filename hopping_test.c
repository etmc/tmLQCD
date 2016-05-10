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
* Test program for the even-odd preconditioned Wilson-Dirac operator
*
*
*******************************************************************************/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <complex.h>
#if (defined BGL && !defined BGP)
#  include <rts.h>
#endif
#ifdef MPI
# include <mpi.h>
# ifdef HAVE_LIBLEMON
#  include <io/params.h>
#  include <io/gauge.h>
# endif
#endif
#include "su3.h"
#include "su3adj.h"
#include "su3spinor.h"
#include "ranlxd.h"
#include "geometry_eo.h"
#include "read_input.h"
#include "start.h"
#include "boundary.h"
#include "operator/Hopping_Matrix.h"
#include "operator/Hopping_Matrix_nocom.h"
#include "operator/tm_operators.h"
#include "global.h"
#include "xchange/xchange.h"
#include "init/init.h"
#include "test/check_geometry.h"
#include "operator/D_psi.h"
#include "phmc.h"
#include "mpi_init.h"
#include "io/io_cm.h"

#ifdef PARALLELT
#  define SLICE (LX*LY*LZ/2)
#elif defined PARALLELXT
#  define SLICE ((LX*LY*LZ/2)+(T*LY*LZ/2))
#elif defined PARALLELXYT
#  define SLICE ((LX*LY*LZ/2)+(T*LY*LZ/2) + (T*LX*LZ/2))
#elif defined PARALLELXYZT
#  define SLICE ((LX*LY*LZ/2)+(T*LY*LZ/2) + (T*LX*LZ/2) + (T*LX*LY/2))
#elif defined PARALLELX
#  define SLICE ((LY*LZ*T/2))
#elif defined PARALLELXY
#  define SLICE ((LY*LZ*T/2) + (LX*LZ*T/2))
#elif defined PARALLELXYZ
#  define SLICE ((LY*LZ*T/2) + (LX*LZ*T/2) + (LX*LY*T/2))
#endif


#define  MAX(A, B)  ((A) > (B) ? (A) : (B))
#define  MIN(A, B)  ((A) < (B) ? (A) : (B))

#if (defined BGL && !defined BGP)
static double clockspeed=1.0e-6/700.0;

double bgl_wtime() {
  return ( rts_get_timebase() * clockspeed );
}
#else
# ifdef MPI
double bgl_wtime() { return(MPI_Wtime()); }
# else
double bgl_wtime() { return(0); }
# endif
#endif

int check_xchange();

int main(int argc,char *argv[])
{
  int j,j_max,k,k_max = 2;
  paramsXlfInfo *xlfInfo;
  int ix, n, *nn,*mm,i;
  double delta, deltamax;
  spinor rsp;
  int status = 0;
#ifdef MPI
  DUM_DERI = 6;
  DUM_SOLVER = DUM_DERI+2;
  DUM_MATRIX = DUM_SOLVER+6;
  NO_OF_SPINORFIELDS = DUM_MATRIX+2;

  MPI_Init(&argc, &argv);
#endif
  g_rgi_C1 = 1.; 
  
  /* Read the input file */
  read_input("hopping_test.input");
  
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
#ifdef _INDEX_INDEP_GEOM
    printf("# the code was compiled with index independent geometry\n");
#endif
#ifdef MPI
#  ifdef _NON_BLOCKING
    printf("# the code was compiled for non-blocking MPI calls (spinor and gauge)\n");
#  endif
#  ifdef _USE_TSPLITPAR
    printf("# the code was compiled with tsplit parallelization\n");
#  endif
#endif
    printf("\n");
    fflush(stdout);
  }

  
#ifdef _GAUGE_COPY
  init_gauge_field(VOLUMEPLUSRAND + g_dbw2rand, 1);
#else
  init_gauge_field(VOLUMEPLUSRAND + g_dbw2rand, 0);
#endif
  init_geometry_indices(VOLUMEPLUSRAND + g_dbw2rand);

  if(even_odd_flag) {
    j = init_spinor_field(VOLUMEPLUSRAND/2, 2*k_max+1);
  }
  else {
    j = init_spinor_field(VOLUMEPLUSRAND, 2*k_max);
  }

  if ( j!= 0) {
    fprintf(stderr, "Not enough memory for spinor fields! Aborting...\n");
    exit(0);
  }
  j = init_moment_field(VOLUME, VOLUMEPLUSRAND);
  if ( j!= 0) {
    fprintf(stderr, "Not enough memory for moment fields! Aborting...\n");
    exit(0);
  }
  
  if(g_proc_id == 0) {
    fprintf(stdout,"The number of processes is %d \n",g_nproc);
    printf("# The lattice size is %d x %d x %d x %d\n",
	   (int)(T*g_nproc_t), (int)(LX*g_nproc_x), (int)(LY*g_nproc_y), (int)(g_nproc_z*LZ));
    printf("# The local lattice size is %d x %d x %d x %d\n", 
	   (int)(T), (int)(LX), (int)(LY),(int) LZ);
    if(even_odd_flag) {
      printf("# testinging the even/odd preconditioned Dirac operator\n");
    }
    else {
      printf("# testinging the standard Dirac operator\n");
    }
    fflush(stdout);
  }
  
  /* define the geometry */
  geometry();
  /* define the boundary conditions for the fermion fields */
  boundary(g_kappa);

#ifdef _USE_HALFSPINOR
  j = init_dirac_halfspinor();
  if ( j!= 0) {
    fprintf(stderr, "Not enough memory for halfspinor fields! Aborting...\n");
    exit(0);
  }
  if(g_sloppy_precision_flag == 1) {
    g_sloppy_precision = 1;
    j = init_dirac_halfspinor32();
    if ( j!= 0) {
      fprintf(stderr, "Not enough memory for 32-Bit halfspinor fields! Aborting...\n");
      exit(0);
    }
  }
#  if (defined _PERSISTENT)
  init_xchange_halffield();
#  endif
#endif  

  status = check_geometry();
  if (status != 0) {
    fprintf(stderr, "Checking of geometry failed. Unable to proceed.\nAborting....\n");
    exit(1);
  }

#if (defined MPI && !(defined _USE_SHMEM))
  check_xchange(); 
#endif

  start_ranlux(1, 123456);

  xlfInfo = construct_paramsXlfInfo(0.5, 0);

  random_gauge_field(reproduce_randomnumber_flag, g_gauge_field);
  if ( startoption == 2 ) {  /* restart */ 
    write_gauge_field(gauge_input_filename,gauge_precision_write_flag,xlfInfo);
  } else if ( startoption == 0 ) { /* cold */
    unit_g_gauge_field();
  } else if (startoption == 3 ) { /* continue */
    read_gauge_field(gauge_input_filename,g_gauge_field);
  } else if ( startoption == 1 ) { /* hot */
  }


#ifdef MPI
  /*For parallelization: exchange the gaugefield */
  xchange_gauge(g_gauge_field);
#endif

  if(even_odd_flag) {
    /*initialize the pseudo-fermion fields*/
    j_max=1;
    for (k = 0; k < k_max; k++) {
      random_spinor_field_eo(g_spinor_field[k], reproduce_randomnumber_flag, RN_GAUSS);
    }

    if (read_source_flag == 2) { /* save */
      /* even first, odd second */
      write_spinorfield_cm_single(g_spinor_field[0],g_spinor_field[1],SourceInfo.basename);  
    }	else if (read_source_flag == 1) { /* yes */
      /* even first, odd second */
      read_spinorfield_cm_single(g_spinor_field[0],g_spinor_field[1],SourceInfo.basename,-1,0); 
# if (!defined MPI)
      if (write_cp_flag == 1) {
	strcat(SourceInfo.basename,".2");
	read_spinorfield_cm_single(g_spinor_field[2],g_spinor_field[3],SourceInfo.basename,-1,0); 

	nn=(int*)calloc(VOLUME,sizeof(int));
	if((void*)nn == NULL) return(100);
	mm=(int*)calloc(VOLUME,sizeof(int));
	if((void*)mm == NULL) return(100);

	n=0;
	deltamax=0.0;
	for(ix=0;ix<VOLUME/2;ix++){
	  (rsp.s0).c0 = (g_spinor_field[2][ix].s0).c0 - (g_spinor_field[0][ix].s0).c0;
	  (rsp.s0).c1 = (g_spinor_field[2][ix].s0).c1 - (g_spinor_field[0][ix].s0).c1;
	  (rsp.s0).c2 = (g_spinor_field[2][ix].s0).c2 - (g_spinor_field[0][ix].s0).c2;
	  (rsp.s1).c0 = (g_spinor_field[2][ix].s1).c0 - (g_spinor_field[0][ix].s1).c0;
	  (rsp.s1).c1 = (g_spinor_field[2][ix].s1).c1 - (g_spinor_field[0][ix].s1).c1;
	  (rsp.s1).c2 = (g_spinor_field[2][ix].s1).c2 - (g_spinor_field[0][ix].s1).c2;
	  (rsp.s2).c0 = (g_spinor_field[2][ix].s2).c0 - (g_spinor_field[0][ix].s2).c0;
	  (rsp.s2).c1 = (g_spinor_field[2][ix].s2).c1 - (g_spinor_field[0][ix].s2).c1;
	  (rsp.s2).c2 = (g_spinor_field[2][ix].s2).c2 - (g_spinor_field[0][ix].s2).c2;
	  (rsp.s3).c0 = (g_spinor_field[2][ix].s3).c0 - (g_spinor_field[0][ix].s3).c0;
	  (rsp.s3).c1 = (g_spinor_field[2][ix].s3).c1 - (g_spinor_field[0][ix].s3).c1;
	  (rsp.s3).c2 = (g_spinor_field[2][ix].s3).c2 - (g_spinor_field[0][ix].s3).c2;
	  _spinor_norm_sq(delta,rsp);
	  if (delta > 1.0e-12) {
	    nn[n] = g_eo2lexic[ix];
	    mm[n]=ix;
	    n++;	    
	  }
	  if(delta>deltamax) deltamax=delta;
	}
	if (n>0){
	  printf("mismatch in even spincolorfield in %d points:\n",n);
	  for(i=0; i< MIN(n,1000); i++){
	    printf("%d,(%d,%d,%d,%d):%f vs. %f\n",nn[i],g_coord[nn[i]][0],g_coord[nn[i]][1],g_coord[nn[i]][2],g_coord[nn[i]][3],creal((g_spinor_field[2][mm[i]].s0).c0), creal((g_spinor_field[0][mm[i]].s0).c0));fflush(stdout);
	  }
	}
	n = 0;
	for(ix=0;ix<VOLUME/2;ix++){
	  (rsp.s0).c0 = (g_spinor_field[3][ix].s0).c0 - (g_spinor_field[1][ix].s0).c0;
	  (rsp.s0).c1 = (g_spinor_field[3][ix].s0).c1 - (g_spinor_field[1][ix].s0).c1;
	  (rsp.s0).c2 = (g_spinor_field[3][ix].s0).c2 - (g_spinor_field[1][ix].s0).c2;
	  (rsp.s1).c0 = (g_spinor_field[3][ix].s1).c0 - (g_spinor_field[1][ix].s1).c0;
	  (rsp.s1).c1 = (g_spinor_field[3][ix].s1).c1 - (g_spinor_field[1][ix].s1).c1;
	  (rsp.s1).c2 = (g_spinor_field[3][ix].s1).c2 - (g_spinor_field[1][ix].s1).c2;
	  (rsp.s2).c0 = (g_spinor_field[3][ix].s2).c0 - (g_spinor_field[1][ix].s2).c0;
	  (rsp.s2).c1 = (g_spinor_field[3][ix].s2).c1 - (g_spinor_field[1][ix].s2).c1;
	  (rsp.s2).c2 = (g_spinor_field[3][ix].s2).c2 - (g_spinor_field[1][ix].s2).c2;
	  (rsp.s3).c0 = (g_spinor_field[3][ix].s3).c0 - (g_spinor_field[1][ix].s3).c0;
	  (rsp.s3).c1 = (g_spinor_field[3][ix].s3).c1 - (g_spinor_field[1][ix].s3).c1;
	  (rsp.s3).c2 = (g_spinor_field[3][ix].s3).c2 - (g_spinor_field[1][ix].s3).c2;
	  _spinor_norm_sq(delta,rsp);
	  if (delta > 1.0e-12) {
	    nn[n]=g_eo2lexic[ix+(VOLUME+RAND)/2];
	    mm[n]=ix;
	    n++;	    
	  }
	  if(delta>deltamax) deltamax=delta;
	}
	if (n>0){
	  printf("mismatch in odd spincolorfield in %d points:\n",n);
	  for(i=0; i< MIN(n,1000); i++){
	    printf("%d,(%d,%d,%d,%d):%f vs. %f\n",nn[i],g_coord[nn[i]][0],g_coord[nn[i]][1],g_coord[nn[i]][2],g_coord[nn[i]][3],creal(g_spinor_field[3][mm[i]].s0.c0), creal(g_spinor_field[1][mm[i]].s0.c0));fflush(stdout);
	  }
	}
	printf("max delta=%e",deltamax);fflush(stdout);
      }
# endif
    }
    
    if (read_source_flag > 0 && write_cp_flag == 0) { /* read-source yes or nobutsave; checkpoint no */
      /* first spinorial arg is output, the second is input */
      Hopping_Matrix(1, g_spinor_field[1], g_spinor_field[0]);      /*ieo=1 M_{eo}*/
      Hopping_Matrix(0, g_spinor_field[0], g_spinor_field[1]);      /*ieo=0 M_{oe}*/
      strcat(SourceInfo.basename,".out");
      write_spinorfield_cm_single(g_spinor_field[0],g_spinor_field[1],SourceInfo.basename);
      printf("Check-field printed. Exiting...\n");
      fflush(stdout);
    }

#ifdef MPI
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
#endif
  }

  free_gauge_field();
  free_geometry_indices();
  free_spinor_field();
  free_moment_field();
  return(0);
}

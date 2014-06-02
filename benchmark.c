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
* Benchmark program for the even-odd preconditioned Wilson-Dirac operator
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
#ifdef OMP
# include <omp.h>
# include "init/init_openmp.h"
#endif
#include "gettime.h"
#include "su3.h"
#include "su3adj.h"
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

int check_xchange();

int main(int argc,char *argv[])
{
  int j,j_max,k,k_max = 1;
#ifdef HAVE_LIBLEMON
  paramsXlfInfo *xlfInfo;
#endif
  int status = 0;
  
  static double t1,t2,dt,sdt,dts,qdt,sqdt;
  double antioptaway=0.0;

#ifdef MPI
  static double dt2;
  
  DUM_DERI = 6;
  DUM_SOLVER = DUM_DERI+2;
  DUM_MATRIX = DUM_SOLVER+6;
  NO_OF_SPINORFIELDS = DUM_MATRIX+2;

#  ifdef OMP
  int mpi_thread_provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_SERIALIZED, &mpi_thread_provided);
#  else
  MPI_Init(&argc, &argv);
#  endif
  MPI_Comm_rank(MPI_COMM_WORLD, &g_proc_id);

#else
  g_proc_id = 0;
#endif

  g_rgi_C1 = 1.; 

    /* Read the input file */
  if((status = read_input("benchmark.input")) != 0) {
    fprintf(stderr, "Could not find input file: benchmark.input\nAborting...\n");
    exit(-1);
  }

#ifdef OMP
  init_openmp();
#endif

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
    printf("# The code was compiled with -D_USE_SHMEM\n");
#  ifdef _PERSISTENT
    printf("# The code was compiled for persistent MPI calls (halfspinor only)\n");
#  endif
#endif
#ifdef MPI
#  ifdef _NON_BLOCKING
    printf("# The code was compiled for non-blocking MPI calls (spinor and gauge)\n");
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
  j = init_moment_field(VOLUME, VOLUMEPLUSRAND + g_dbw2rand);
  if ( j!= 0) {
    fprintf(stderr, "Not enough memory for moment fields! Aborting...\n");
    exit(0);
  }
  
  if(g_proc_id == 0) {
    fprintf(stdout,"# The number of processes is %d \n",g_nproc);
    printf("# The lattice size is %d x %d x %d x %d\n",
	   (int)(T*g_nproc_t), (int)(LX*g_nproc_x), (int)(LY*g_nproc_y), (int)(g_nproc_z*LZ));
    printf("# The local lattice size is %d x %d x %d x %d\n", 
	   (int)(T), (int)(LX), (int)(LY),(int) LZ);
    if(even_odd_flag) {
      printf("# benchmarking the even/odd preconditioned Dirac operator\n");
    }
    else {
      printf("# benchmarking the standard Dirac operator\n");
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
  random_gauge_field(reproduce_randomnumber_flag, g_gauge_field);

#ifdef MPI
  /*For parallelization: exchange the gaugefield */
  xchange_gauge(g_gauge_field);
#endif

  if(even_odd_flag) {
    /*initialize the pseudo-fermion fields*/
    j_max=2048;
    sdt=0.;
    for (k = 0; k < k_max; k++) {
      random_spinor_field_eo(g_spinor_field[k], reproduce_randomnumber_flag, RN_GAUSS);
    }
    
    while(sdt < 30.) {
#ifdef MPI
      MPI_Barrier(MPI_COMM_WORLD);
#endif
      t1 = gettime();
      antioptaway=0.0;
      for (j=0;j<j_max;j++) {
        for (k=0;k<k_max;k++) {
          Hopping_Matrix(0, g_spinor_field[k+k_max], g_spinor_field[k]);
          Hopping_Matrix(1, g_spinor_field[2*k_max], g_spinor_field[k+k_max]);
          antioptaway+=creal(g_spinor_field[2*k_max][0].s0.c0);
        }
      }
      t2 = gettime();
      dt = t2-t1;
#ifdef MPI
      MPI_Allreduce (&dt, &sdt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
      sdt = dt;
#endif
      qdt=dt*dt;
#ifdef MPI
      MPI_Allreduce (&qdt, &sqdt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
      sqdt = qdt;
#endif
      sdt=sdt/((double)g_nproc);
      sqdt=sqrt(sqdt/g_nproc-sdt*sdt);
      j_max*=2;
    }
    j_max=j_max/2;
    dts=dt;
    sdt=1.0e6f*sdt/((double)(k_max*j_max*(VOLUME)));
    sqdt=1.0e6f*sqdt/((double)(k_max*j_max*(VOLUME)));
    
    if(g_proc_id==0) {
      printf("# The following result is just to make sure that the calculation is not optimized away: %e\n", antioptaway);
      printf("# Total compute time %e sec, variance of the time %e sec. (%d iterations).\n", sdt, sqdt, j_max);
      printf("# Communication switched on:\n# (%d Mflops [%d bit arithmetic])\n", (int)(1608.0f/sdt),(int)sizeof(spinor)/3);
#ifdef OMP
      printf("# Mflops per OpenMP thread ~ %d\n",(int)(1608.0f/(omp_num_threads*sdt)));
#endif
      printf("\n");
      fflush(stdout);
    }
    
#ifdef MPI
    /* isolated computation */
    t1 = gettime();
    antioptaway=0.0;
    for (j=0;j<j_max;j++) {
      for (k=0;k<k_max;k++) {
        Hopping_Matrix_nocom(0, g_spinor_field[k+k_max], g_spinor_field[k]);
        Hopping_Matrix_nocom(1, g_spinor_field[2*k_max], g_spinor_field[k+k_max]);
        antioptaway += creal(g_spinor_field[2*k_max][0].s0.c0);
      }
    }
    t2 = gettime();
    dt2 = t2-t1;
    /* compute the bandwidth */
    dt=dts-dt2;
    MPI_Allreduce (&dt, &sdt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    sdt=sdt/((double)g_nproc);
    MPI_Allreduce (&dt2, &dt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    dt=dt/((double)g_nproc);
    dt=1.0e6f*dt/((double)(k_max*j_max*(VOLUME)));
    if(g_proc_id==0) {
      printf("# The following result is printed just to make sure that the calculation is not optimized away: %e\n",antioptaway);
      printf("# Communication switched off: \n# (%d Mflops [%d bit arithmetic])\n", (int)(1608.0f/dt),(int)sizeof(spinor)/3);
#ifdef OMP
      printf("# Mflops per OpenMP thread ~ %d\n",(int)(1608.0f/(omp_num_threads*dt)));
#endif
      printf("\n"); 
      fflush(stdout);
    }
    sdt=sdt/((double)k_max);
    sdt=sdt/((double)j_max);
    sdt=sdt/((double)(2*SLICE));
    if(g_proc_id==0) {
      printf("# The size of the package is %d bytes.\n",(SLICE)*192);
#ifdef _USE_HALFSPINOR
      printf("# The bandwidth is %5.2f + %5.2f MB/sec\n", 192./sdt/1024/1024, 192./sdt/1024./1024);
#else
      printf("# The bandwidth is %5.2f + %5.2f MB/sec\n", 2.*192./sdt/1024/1024, 2.*192./sdt/1024./1024);
#endif
    }
#endif
    fflush(stdout);
  }
  else {
    /* the non even/odd case now */
    /*initialize the pseudo-fermion fields*/
    j_max=1;
    sdt=0.;
    for (k=0;k<k_max;k++) {
      random_spinor_field_lexic(g_spinor_field[k], reproduce_randomnumber_flag, RN_GAUSS);
    }
    
    while(sdt < 3.) {
#ifdef MPI
      MPI_Barrier(MPI_COMM_WORLD);
#endif
      t1 = gettime();
      for (j=0;j<j_max;j++) {
        for (k=0;k<k_max;k++) {
          D_psi(g_spinor_field[k+k_max], g_spinor_field[k]);
          antioptaway+=creal(g_spinor_field[k+k_max][0].s0.c0);
        }
      }
      t2 = gettime();
      dt=t2-t1;
#ifdef MPI
      MPI_Allreduce (&dt, &sdt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
      sdt = dt;
#endif
      qdt=dt*dt;
#ifdef MPI
      MPI_Allreduce (&qdt, &sqdt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
      sqdt = qdt;
#endif
      sdt=sdt/((double)g_nproc);
      sqdt=sqrt(sqdt/g_nproc-sdt*sdt);
      j_max*=2;
    }
    j_max=j_max/2;
    dts=dt;
    sdt=1.0e6f*sdt/((double)(k_max*j_max*(VOLUME)));
    sqdt=1.0e6f*sqdt/((double)(k_max*j_max*(VOLUME)));

    if(g_proc_id==0) {
      printf("# The following result is just to make sure that the calculation is not optimized away: %e\n", antioptaway);
      printf("# Total compute time %e sec, variance of the time %e sec. (%d iterations).\n", sdt, sqdt, j_max);
      printf("\n# (%d Mflops [%d bit arithmetic])\n", (int)(1680.0f/sdt),(int)sizeof(spinor)/3);
#ifdef OMP
      printf("# Mflops per OpenMP thread ~ %d\n",(int)(1680.0f/(omp_num_threads*sdt)));
#endif
      printf("\n"); 
      fflush(stdout);
    }
  }
#ifdef HAVE_LIBLEMON
  if(g_proc_id==0) {
    printf("# Performing parallel IO test ...\n");
  }
  xlfInfo = construct_paramsXlfInfo(0.5, 0);
  write_gauge_field( "conf.test", 64, xlfInfo);
  free(xlfInfo);
  if(g_proc_id==0) {
    printf("# done ...\n");
  }
#endif


#ifdef OMP
  free_omp_accumulators();
#endif
  free_gauge_field();
  free_geometry_indices();
  free_spinor_field();
  free_moment_field();
#ifdef MPI
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
#endif
  return(0);
}

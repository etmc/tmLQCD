/* $Id$ */
/*******************************************************************************
*
* Benchmark program for the even-odd preconditioned Wilson-Dirac operator
*
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
#include <string.h>
#if (defined BGL)
#  include <rts.h>
#endif
#include "su3.h"
#include "su3adj.h"
#include "ranlxd.h"
#include "geometry_eo.h"
#include "read_input.h"
#include "start.h"
#include "boundary.h"
#include "Hopping_Matrix.h"
#include "Hopping_Matrix_nocom.h"
#include "global.h"
#include "xchange.h"
#include "init_gauge_field.h"
#include "init_geometry_indices.h"
#include "init_spinor_field.h"
#include "init_moment_field.h"
#include "init_dirac_halfspinor.h"
#include "update_backward_gauge.h"
#include "test/check_geometry.h"
#include "mpi_init.h"

#ifdef PARALLELT
#  define SLICE (LX*LY*LZ/2)
#elif defined PARALLELXT
#  define SLICE ((LX*LY*LZ/2)+(T*LY*LZ/2))
#elif defined PARALLELXYT
#  define SLICE ((LX*LY*LZ/2)+(T*LY*LZ/2) + (T*LX*LZ))
#elif defined PARALLELXYZT
#  define SLICE ((LX*LY*LZ/2)+(T*LY*LZ/2) + (T*LX*LZ) + (T*LX*LY))
#endif

#ifdef BGL
static double clockspeed=1.0e-6/700.0;

double bgl_wtime() {
  return ( rts_get_timebase() * clockspeed );
}
#endif

int check_xchange();

int main(int argc,char *argv[])
{
  int j,j_max,k,k_max = 5;
#ifdef _GAUGE_COPY
  int kb=0;
#endif
  
  
  static double t1,t2,dt,sdt,dts,qdt,sqdt;
#ifdef MPI
  static double dt2;
  int rlxd_state[105];
  
  MPI_Init(&argc, &argv);
#endif
  g_rgi_C1 = 1.; 
  
  /* Read the input file */
  read_input("benchmark.input");
  
  mpi_init(argc, argv);
  
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
    
    printf("\n");
    fflush(stdout);
  }
  
  
  
#ifdef _GAUGE_COPY
  init_gauge_field(VOLUMEPLUSRAND + g_dbw2rand, 1);
#else
  init_gauge_field(VOLUMEPLUSRAND + g_dbw2rand, 0);
#endif
  init_geometry_indices(VOLUMEPLUSRAND + g_dbw2rand);
  j = init_spinor_field(VOLUMEPLUSRAND/2, 3*k_max);
  init_dirac_halfspinor();

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
    fflush(stdout);
  }
  
  /* define the geometry */
  geometry();
  /* define the boundary conditions for the fermion fields */
  boundary();
  
  check_geometry();
#ifdef MPI
  check_xchange(); 
#endif
  
  if(reproduce_randomnumber_flag == 1) {
    /* here we generate exactly the same configuration as for the 
       single node simulation */
    if(g_proc_id==0) {
      rlxd_init(1, 123456);   
      random_gauge_field();
      /*  send the state of the random-number generator to 1 */
#ifdef MPI
      rlxd_get(rlxd_state);
      MPI_Send((void*)rlxd_state, 105, MPI_INT, 1, 99, MPI_COMM_WORLD);
#endif
    }
#ifdef MPI
    else {
      /* recieve the random number state form g_proc_id-1 */
      MPI_Recv((void*)rlxd_state, 105, MPI_INT, g_proc_id-1, 99, MPI_COMM_WORLD, &status);
      rlxd_reset(rlxd_state);
      random_gauge_field();
      /* send the random number state to g_proc_id+1 */
      k=g_proc_id+1; 
      if(k==g_nproc) k=0;
      rlxd_get(rlxd_state);
      MPI_Send((void*)rlxd_state, 105, MPI_INT, k, 99, MPI_COMM_WORLD);
    }
    if(g_proc_id==0) {
      MPI_Recv((void*)rlxd_state, 105, MPI_INT,g_nproc-1,99, MPI_COMM_WORLD, &status);
      rlxd_reset(rlxd_state);
    }
#endif
  }
  else {
    rlxd_init(1, 123456 + g_proc_id*97);
    random_gauge_field();
  }

#ifdef MPI
  /*For parallelization: exchange the gaugefield */
  xchange_gauge();
#endif


#ifdef _GAUGE_COPY
  update_backward_gauge();
#endif

  /*initialize the pseudo-fermion fields*/
  j_max=1;
  sdt=0.;
  for (k=0;k<k_max;k++) {
    random_spinor_field(g_spinor_field[k], VOLUME/2, reproduce_randomnumber_flag);
  }

  while(sdt < 30.) {
#ifdef MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
#if defined BGL
    t1 = bgl_wtime();
#else
    t1=(double)clock();
#endif
    for (j=0;j<j_max;j++) {
      for (k=0;k<k_max;k++) {
	Hopping_Matrix(0, g_spinor_field[k+k_max], g_spinor_field[k]);
	Hopping_Matrix(1, g_spinor_field[k+2*k_max], g_spinor_field[k+k_max]);
      }
    }
#if defined BGL
    t2 = bgl_wtime();
    dt = t2 - t1;
#else
    t2=(double)clock();
    dt=(t2-t1)/((double)(CLOCKS_PER_SEC));
#endif
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
    printf("total time %e sec, Variance of the time %e sec \n",sdt,sqdt);
    printf("\n");
    printf(" (%d Mflops [%d bit arithmetic])\n",
	   (int)(1320.0f/sdt),(int)sizeof(spinor)/3);
    printf("\n");
    fflush(stdout);
  }

#ifdef MPI
  /* isolated computation */
#if defined BGL
    t1 = bgl_wtime();
#else
    t1=(double)clock();
#endif
  for (j=0;j<j_max;j++) {
    for (k=0;k<k_max;k++) {
      Hopping_Matrix_nocom(0, g_spinor_field[k+k_max], g_spinor_field[k]);
      Hopping_Matrix_nocom(1, g_spinor_field[k+2*k_max], g_spinor_field[k+k_max]);
    }
  }
#if defined BGL
    t2 = bgl_wtime();
    dt2 = t2 - t1;
#else
    t2=(double)clock();
    dt2=(t2-t1)/((double)(CLOCKS_PER_SEC));
#endif
  /* compute the bandwidth */
  dt=dts-dt2;
  MPI_Allreduce (&dt, &sdt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  sdt=sdt/((double)g_nproc);
  MPI_Allreduce (&dt2, &dt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  dt=dt/((double)g_nproc);
  dt=1.0e6f*dt/((double)(k_max*j_max*(VOLUME)));
  if(g_proc_id==0) {
    printf("communication switched off \n");
    printf(" (%d Mflops [%d bit arithmetic])\n",
	   (int)(1320.0f/dt),(int)sizeof(spinor)/3);
    printf("\n");
    fflush(stdout);
  }
  sdt=sdt/((double)k_max);
  sdt=sdt/((double)j_max);
  sdt=sdt/((double)(2*SLICE));
  if(g_proc_id==0) {
    printf("The size of the package is %d Byte \n",(SLICE)*192);
    printf("The bandwidth is %5.2f + %5.2f   MB/sec\n",
	   2.*192./sdt/1024/1024, 2.*192./sdt/1024./1024);
    fflush(stdout);
  }
#endif
  fflush(stdout);
#ifdef MPI
  MPI_Finalize();
#endif
  free_gauge_field();
  free_geometry_indices();
  free_spinor_field();
  free_moment_field();
  return(0);
}

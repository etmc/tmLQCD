/*******************************************************************************
*
* File hybrid.c
*
* Benchmark program for the even-odd preconditioned Wilson-Dirac operator
*
* Author: Martin Hasenbusch
* Date: Wed, Aug 29, 2001 02:06:26 PM
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
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
#include "test/check_geometry.h"
#include "mpi_init.h"

#ifndef PARALLELXT
#define SLICE (LX*LY*LZ/2)
#else
#define SLICE ((LX*LY*LZ/2)+(T*LY*LZ/2))
#endif

int check_xchange();

int main(int argc,char *argv[])
{
  int j,j_max,k,k_max;
#ifdef _GAUGE_COPY
  int kb=0;
#endif


  static double t1,t2,dt,sdt,dts,qdt,sqdt;
#ifdef MPI
  static double dt2;
  int rlxd_state[105];
#endif

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
#ifdef _NEW_GEOMETRY
    printf("# The code was compiled with -D_NEW_GEOMETRY\n");
#endif
#ifdef _GAUGE_COPY
    printf("# The code was compiled with -D_GAUGE_COPY\n");
#endif

    printf("\n");
    fflush(stdout);
  }

  if(g_rgi_C1 == 0.) {
    g_dbw2rand = 0;
  }
#ifndef MPI
  g_dbw2rand = 0;
#endif

#ifdef _GAUGE_COPY
  init_gauge_field(VOLUMEPLUSRAND, 1);
#else
  init_gauge_field(VOLUMEPLUSRAND, 0);
#endif
  init_geometry_indices(VOLUMEPLUSRAND);
  j = init_spinor_field(VOLUMEPLUSRAND/2, NO_OF_SPINORFIELDS);
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
    fprintf(stdout,"The local lattice size is %d x %d x %d^2 \n\n",T,LX,L);
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

  /*For parallelization: exchange the gaugefield */
  xchange_gauge();
#endif


#ifdef _GAUGE_COPY
  /* set the backward gauge field */
  for(k=0;k<VOLUME;k++) {
    kb=g_idn[k][0];
    _su3_assign(g_gauge_field_back[k][0],g_gauge_field[kb][0]);
    kb=g_idn[k][1];
    _su3_assign(g_gauge_field_back[k][1],g_gauge_field[kb][1]);
    kb=g_idn[k][2];
    _su3_assign(g_gauge_field_back[k][2],g_gauge_field[kb][2]);
    kb=g_idn[k][3];
    _su3_assign(g_gauge_field_back[k][3],g_gauge_field[kb][3]);
  }
#endif

  /*initialize the pseudo-fermion fields*/
  k_max=(NO_OF_SPINORFIELDS)/3;
  j_max=1;
  sdt=0.;
  for (k=0;k<k_max;k++) {
    random_spinor_field(k);
  }

  while(sdt<30.) {
#ifdef MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    t1=(double)clock();
    for (j=0;j<j_max;j++) {
      for (k=0;k<k_max;k++) {
	Hopping_Matrix(0, spinor_field[k+k_max], spinor_field[k]);
	Hopping_Matrix(1, spinor_field[k+2*k_max], spinor_field[k+k_max]);
      }
    }
    t2=(double)clock();

    dt=(t2-t1)/((double)(CLOCKS_PER_SEC));
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
    sdt=sdt/g_nproc;
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
  t1=(double)clock();
  for (j=0;j<j_max;j++) {
    for (k=0;k<k_max;k++) {
      Hopping_Matrix_nocom(0, spinor_field[k+k_max], spinor_field[k]);
      Hopping_Matrix_nocom(1, spinor_field[k+2*k_max], spinor_field[k+k_max]);
    }
  }
  t2=(double)clock();

  dt2=(t2-t1)/((double)(CLOCKS_PER_SEC));
  /* compute the bandwidth */
  dt=dts-dt2;
  MPI_Allreduce (&dt, &sdt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  sdt=sdt/g_nproc;
  MPI_Allreduce (&dt2, &dt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  dt=dt/g_nproc;
  dt=1.0e6f*dt/((double)(k_max*j_max*(VOLUME)));
  if(g_proc_id==0) {
    printf("communication switched off \n");
    printf(" (%d Mflops [%d bit arithmetic])\n",
	   (int)(1320.0f/dt),(int)sizeof(spinor)/3);
    printf("\n");
    fflush(stdout);
  }
  sdt=sdt/k_max;
  sdt=sdt/j_max;
  sdt=sdt/(2*SLICE);
  if(g_proc_id==0) {
    printf("The size of the package is %d Byte \n",(SLICE)*192);
    printf("The bandwidth is %5.2f + %5.2f   MB/sec\n",
	   0.000001*2.*192./sdt, 0.000001*2.*192./sdt);
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

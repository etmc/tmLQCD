/*******************************************************************************
 * $Id$
 *
 * File time_Hopping_Matrix.c
 *
 * Timing of the program Hopping_Matrix
 *
 * Author: Carsten Urbach
 *         urbach@physik.fu-berlin.de
 *
 *******************************************************************************/

#define MAIN_PROGRAM
#ifdef XLC
#define CACHE_SIZE 1500000
#else
#define CACHE_SIZE 256000
#endif

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "global.h"
#include "su3.h"
#include "ranlxd.h"
#include "geometry_eo.h"
#include "start.h"
#include "Hopping_Matrix.h"
#include "linalg_eo.h"
#include "boundary.h"
#if defined MPI
#include "xchange.h"
#endif

int main(int argc,char *argv[]){
  int i,k,kmax,n,count;
  float t1,t2,dt;
  double m0;
  int a,b;

#if defined MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &g_nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &g_proc_id);

  if(g_proc_id == 0) {
    printf("g_proc_id=%d g_nproc=%d \n",
	   g_proc_id,g_nproc);
  }
#else
  g_proc_id = 0;
  g_nproc = 1;
#endif
  if(g_proc_id == 0) {
    printf("\n");
    printf("Timing of Hopping_Matrix (random spinor and gauge fields)\n");
    printf("------------------------------------------------\n");
    
    printf("\n");
    printf("The lattice size is %d x %d^3\n\n",T*g_nproc, L);
  }
  kmax=((int)(CACHE_SIZE))/(96*(int)(VOLUME));
  if (kmax==0){
    kmax=1;
  }
  if(kmax>10){
    kmax=10;
  }


  n=30000000/(int)(VOLUME/2);
  if (n<2){
    n=2;
  }
  g_kappa = 0.125;
  g_mu = 0.1;

  rlxd_init(1,123456);   
  geometry();
  boundary();
  random_gauge_field();
#if defined MPI
  xchange_gauge();
#endif

  for (k=0;k<2*kmax;k++){
    random_spinor_field(k);
  }

  k=0;
  t1=(float)clock();
  for (count=0;count<n;count++){      
    Hopping_Matrix(EO, k+kmax, k);
    k++;
    if (k==kmax){
      k=0;
    }
  }      
  t2=(float)clock();

  dt=(t2-t1)/((float)(CLOCKS_PER_SEC));
  dt=1.0e6f*dt/((float)(n*(VOLUME/2)));

  printf("Time per lattice point: %4.3f micro sec",dt);
  printf(" (%d Mflops [%d bit arithmetic])\n",
	 (int)(1392.0f/dt),(int)sizeof(spinor)/3);   
  printf("\n");

  a = (int)(1392.0f/dt);
  MPI_Reduce(&a, &b, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  if(g_proc_id == 0) {
    printf("Sum: %d, Average: %d Mflops\n", b, b/g_nproc);
    printf("Number of applications of Hopping_Matrix was %d!\n", n);
  }
#if defined MPI
   MPI_Finalize();
#endif
  return(0);
}

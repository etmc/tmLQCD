/*******************************************************************************
 * $Id$
 *
 * File check_geometry.c
 *
 * Consistency of the index arrays ipt, iup and idn
 *
 * Author: Martin Luescher <luscher@mail.desy.de>
 * Date: 24.10.2000
 *
 *******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#ifdef MPI
#include <mpi.h>
#endif
#include "global.h"
#include "geometry_eo.h"
#ifdef MPI
#include "mpi_init.h"
#endif

int main(int argc,char *argv[])
{
  int ix;
  int * itest;
  int x0,x1,x2,x3;
  int iy0,iy1,iy2,iy3;
  int iz0,iz1,iz2,iz3;   
#ifdef MPI
  int r1=0, r2=0;
#endif

  mpi_init(argc, argv);

  itest = calloc(VOLUME+RAND, sizeof(int));
  printf("\n");
  printf("Test of the geometry programs \n");
  printf("------------------------------\n");

  printf("\n");
  printf("The lattice size is %d x %d^3 \n\n",(int)(T),(int)(L));
   
  geometry();

  for (ix=0;ix<VOLUMEPLUSRAND;ix++){
    itest[ix]=0;
  }
  

  for (x0 = 0; x0 < T; x0++){
    for (x1 = 0; x1 < LX; x1++){
      for (x2 = 0; x2 < LY; x2++){
	for (x3 = 0; x3 < LZ; x3++){
	  ix=g_ipt[x0][x1][x2][x3];

	  if ((ix<0)||(ix>=VOLUME)){
	    printf("The index ipt is out of range\n");
	    printf("Program aborted\n");
#ifdef MPI
	    MPI_Finalize();
#endif
	    exit(0);
	  }
               
	  itest[ix]+=1;
	}
      }
    }
  }

  for (ix = 0; ix < VOLUME; ix++){
    if (itest[ix]!=1){
      printf("The index ipt is not one-to-one\n");
      printf("Program aborted\n");
#ifdef MPI
      MPI_Finalize();
#endif
      exit(0);
    }
  }

  for (x0 = 0; x0 < T; x0++){
    for (x1 = 0; x1 < LX; x1++){
      for (x2 = 0; x2 < LY; x2++){
	for (x3 = 0; x3 < LZ; x3++){
	  ix=g_ipt[x0][x1][x2][x3];

	  iy0=g_iup[ix][0];
#if (defined PARALLELT || defined PARALLELXT )
	  if(x0!=T-1) {
	    iz0=g_ipt[(x0+1)%T][x1][x2][x3];
	  }
	  else {
	    iz0 = iy0;
	    itest[iy0]++;
	  }
#else
	  iz0=g_ipt[(x0+1)%T][x1][x2][x3];
#endif     
                  
	  iy1=g_iup[ix][1];
#if (defined PARALLELXT)
	  if(x1 !=0) {
	    iz1=g_ipt[x0][(x1+1)%LX][x2][x3];
	  }
	  else {
	    iz1 = iy1;
	    itest[iy1]++;
	  }
#else
	  iz1=g_ipt[x0][(x1+1)%LX][x2][x3];
#endif

	  iy2=g_iup[ix][2];
	  iz2=g_ipt[x0][x1][(x2+1)%LY][x3];

	  iy3=g_iup[ix][3];
	  iz3=g_ipt[x0][x1][x2][(x3+1)%LZ];               
               
	  if ((iy0!=iz0)||(iy1!=iz1)||(iy2!=iz2)||(iy3!=iz3)){
	    printf("The index iup is incorrect\n");
	    printf("%d %d %d %d\n", (iy0!=iz0), (iy1!=iz1), (iy2!=iz2), (iy3!=iz3));
	    printf("Program aborted\n");
#ifdef MPI
	    MPI_Finalize();
#endif
	    exit(0);
	  }

	  iy0=g_idn[ix][0];
#if (defined PARALLELT || defined PARALLELXT )
	  if(x0 !=0) {
	    iz0=g_ipt[(x0+T-1)%T][x1][x2][x3];
	  }
	  else {
	    iz0 = iy0;
	    itest[iy0]++;
	  }
#else
	  iz0=g_ipt[(x0+T-1)%T][x1][x2][x3];
#endif
	  iy1=g_idn[ix][1];
#if (defined PARALLELXT)
	  if(x1 !=0) {
	    iz1=g_ipt[x0][(x1+LX-1)%LX][x2][x3];
	  }
	  else {
	    iz1 = iy1;
	    itest[iy1]++;
	  }
#else
	  iz1=g_ipt[x0][(x1+LX-1)%LX][x2][x3];
#endif
	  iy2=g_idn[ix][2];
	  iz2=g_ipt[x0][x1][(x2+LY-1)%LY][x3];

	  iy3=g_idn[ix][3];
	  iz3=g_ipt[x0][x1][x2][(x3+LZ-1)%LZ];               
               
	  if ((iy0!=iz0)||(iy1!=iz1)||(iy2!=iz2)||(iy3!=iz3)){
	    printf("The index idn is incorrect\n");
	    printf("%d %d %d %d\n", (iy0!=iz0), (iy1!=iz1), (iy2!=iz2), (iy3!=iz3));
	    printf("Program aborted\n");
#ifdef MPI
	    MPI_Finalize();
#endif
	    exit(0);
	  }
	}
      }
    }
  }

  for (ix = VOLUME; ix < VOLUMEPLUSRAND; ix++){
    if (itest[ix]!=1){
      printf("The boundary is not correctly used\n");
      printf("Program aborted\n");
#ifdef MPI
      MPI_Finalize();
#endif
      exit(0);
    }
  }

  printf("The lattice is correctly mapped by the index arrays\n");
  printf("\n");

#ifdef MPI
  MPI_Finalize();
#endif
  return(0);
}

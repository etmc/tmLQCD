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

/* int main(int argc,char *argv[]) */
int check_geometry()
{
  int ix, j;
  int * stest;
  int * itest;
  int x0,x1,x2,x3;
  int iy0,iy1,iy2,iy3;
  int iz0,iz1,iz2,iz3;   
#ifdef MPI
  int r1=0, r2=0;
#endif
  
/*   mpi_init(argc, argv); */
  
  itest = calloc(VOLUME+RAND, sizeof(int));
  stest = calloc((VOLUME+RAND)/2, sizeof(int));
  printf("\n");
  printf("Test of the geometry programs \n");
  printf("------------------------------\n");
  
  printf("\n");
  printf("The lattice size is %d x %d^3 \n\n",(int)(T),(int)(L));
  
/*   geometry(); */
  
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
    if (itest[ix]!=1) {
      printf("The boundary is not correctly used\n");
      printf("Program aborted\n");
#ifdef MPI
      MPI_Finalize();
#endif
      exit(0);
    }
  }
  

  /* check of EO geometry */

  for (ix=0;ix<VOLUMEPLUSRAND/2;ix++){
    stest[ix]=0;
  }
  
  for(j = 0; j < VOLUME/2; j++) {
    ix = trans2[j];
    
    iy0 = g_idn[ix][0];
    iz0 = trans1[iy0] - VOLUMEPLUSRAND/2;
    if(iz0 > VOLUMEPLUSRAND/2 || iz0 < 0) {
      printf("There is a problem with EO geometry in direction 0-\n");
      printf("%d\n", iz0);
#ifdef MPI
      MPI_Finalize();
#endif
      exit(0);
    }
    stest[iz0] += 1;
    
    iy0 = g_iup[ix][0];
    iz0 = trans1[iy0] - VOLUMEPLUSRAND/2;
    if(iz0 > VOLUMEPLUSRAND/2 || iz0 < 0) {
      printf("There is a problem with EO geometry in direction 0+\n");
#ifdef MPI
      MPI_Finalize();
#endif
      exit(0);
    }
    stest[iz0] += 1;
    
    iy1 = g_idn[ix][1];
    iz1 = trans1[iy1] - VOLUMEPLUSRAND/2;
    if(iz1 > VOLUMEPLUSRAND/2  || iz1 < 0) {
      printf("There is a problem with EO geometry in direction 1-\n");
#ifdef MPI
      MPI_Finalize();
#endif
      exit(0);
    }
    stest[iz1] += 1;
    
    iy1 = g_iup[ix][1];
    iz1 = trans1[iy1] - VOLUMEPLUSRAND/2;
    if(iz1 > VOLUMEPLUSRAND/2  || iz1 < 0) {
      printf("There is a problem with EO geometry in direction 1+\n");
#ifdef MPI
      MPI_Finalize();
#endif
      exit(0);
    }
    stest[iz1] += 1;
    
    iy2 = g_idn[ix][2];
    iz2 = trans1[iy2] - VOLUMEPLUSRAND/2;
    if(iz2 > VOLUMEPLUSRAND/2 || iz2 < 0) {
      printf("There is a problem with EO geometry in direction 2-\n");
#ifdef MPI
      MPI_Finalize();
#endif
      exit(0);
    }
    stest[iz2] += 1;
    
    iy2 = g_iup[ix][2];
    iz2 = trans1[iy2] - VOLUMEPLUSRAND/2;
    if(iz2 > VOLUMEPLUSRAND/2 || iz2 < 0) {
      printf("There is a problem with EO geometry in direction 2+\n");
#ifdef MPI
      MPI_Finalize();
#endif
      exit(0);
    }
    stest[iz2] += 1;
    
    
    iy3 = g_idn[ix][3];
    iz3 = trans1[iy3] - VOLUMEPLUSRAND/2;
    if(iz3 > VOLUMEPLUSRAND/2 || iz3 < 0) {
      printf("There is a problem with EO geometry in direction 3-\n");
#ifdef MPI
      MPI_Finalize();
#endif
      exit(0);
    }
    stest[iz3] += 1;
    
    iy3 = g_iup[ix][3];
    iz3 = trans1[iy3] - VOLUMEPLUSRAND/2;
    if(ix == 0 && iy3 == 1) {
      printf("%d\n", iz3);
    }
    if(iz3 > VOLUMEPLUSRAND/2 || iz3 < 0) {
      printf("There is a problem with EO geometry in direction 3+\n");
#ifdef MPI
      MPI_Finalize();
#endif
      exit(0);
    }
    stest[iz3] += 1;
  }
  iz0 = 0;
  for(j = 0; j < (VOLUME)/2; j++) {
    iz0 += stest[j];
  }
  if(iz0 != 8*(VOLUME)/2-RAND/2) {
    printf("There is a problem in the the even odd geometry\n");
    printf("%d is not equal to 8*(VOLUME)/2-RAND/2=%d\n", iz0, 8*(VOLUME)/2-RAND/2);
#ifdef MPI
    MPI_Finalize();
#endif
    exit(0); 
  }

  for(j = VOLUME; j < VOLUMEPLUSRAND/2; j++) {
    if(stest[j] != 1) {
      printf("There is a problem in the boundary of the even odd geometry\n");
#ifdef MPI
      MPI_Finalize();
#endif
      exit(0);
    }
  }
  

  for (ix=0;ix<VOLUMEPLUSRAND/2;ix++){
    stest[ix]=0;
  }

  for(j = VOLUMEPLUSRAND/2; j < (VOLUME+VOLUMEPLUSRAND)/2; j++) {
    ix = trans2[j];
    
    iy0 = g_idn[ix][0];
    iz0 = trans1[iy0];
    if(iz0 > VOLUMEPLUSRAND/2 || iz0 < 0) {
      printf("There is a problem with EO geometry in direction 0-\n");
      printf("%d\n", iz0);
#ifdef MPI
      MPI_Finalize();
#endif
      exit(0);
    }
    stest[iz0] += 1;
    
    iy0 = g_iup[ix][0];
    iz0 = trans1[iy0];
    if(iz0 > VOLUMEPLUSRAND/2 || iz0 < 0) {
      printf("There is a problem with EO geometry in direction 0+\n");
#ifdef MPI
      MPI_Finalize();
#endif
      exit(0);
    }
    stest[iz0] += 1;
    
    iy1 = g_idn[ix][1];
    iz1 = trans1[iy1];
    if(iz1 > VOLUMEPLUSRAND/2  || iz1 < 0) {
      printf("There is a problem with EO geometry in direction 1-\n");
#ifdef MPI
      MPI_Finalize();
#endif
      exit(0);
    }
    stest[iz1] += 1;
    
    iy1 = g_iup[ix][1];
    iz1 = trans1[iy1];
    if(iz1 > VOLUMEPLUSRAND/2  || iz1 < 0) {
      printf("There is a problem with EO geometry in direction 1+\n");
#ifdef MPI
      MPI_Finalize();
#endif
      exit(0);
    }
    stest[iz1] += 1;
    
    iy2 = g_idn[ix][2];
    iz2 = trans1[iy2];
    if(iz2 > VOLUMEPLUSRAND/2 || iz2 < 0) {
      printf("There is a problem with EO geometry in direction 2-\n");
#ifdef MPI
      MPI_Finalize();
#endif
      exit(0);
    }
    stest[iz2] += 1;
    
    iy2 = g_iup[ix][2];
    iz2 = trans1[iy2];
    if(iz2 > VOLUMEPLUSRAND/2 || iz2 < 0) {
      printf("There is a problem with EO geometry in direction 2+\n");
#ifdef MPI
      MPI_Finalize();
#endif
      exit(0);
    }
    stest[iz2] += 1;
    
    
    iy3 = g_idn[ix][3];
    iz3 = trans1[iy3];
    if(iz3 > VOLUMEPLUSRAND/2 || iz3 < 0) {
      printf("There is a problem with EO geometry in direction 3-\n");
#ifdef MPI
      MPI_Finalize();
#endif
      exit(0);
    }
    stest[iz3] += 1;
    
    iy3 = g_iup[ix][3];
    iz3 = trans1[iy3];
 if(ix == 0 && iy3 == 1) {
      printf("%d\n", iz3);
    }
    if(iz3 > VOLUMEPLUSRAND/2 || iz3 < 0) {
      printf("There is a problem with EO geometry in direction 3+\n");
#ifdef MPI
      MPI_Finalize();
#endif
      exit(0);
    }
    stest[iz3] += 1;
  }
  iz0 = 0;
  for(j = 0; j < (VOLUME)/2; j++) {
    iz0 += stest[j];
  }
  if(iz0 != 8*(VOLUME)/2-RAND/2) {
    printf("There is a problem in the the even odd geometry\n");
    printf("%d is not equal to 8*(VOLUME)/2-RAND/2=%d\n", iz0, 8*(VOLUME)/2-RAND/2);
#ifdef MPI
    MPI_Finalize();
#endif
    exit(0); 
  }

  for(j = VOLUME; j < VOLUMEPLUSRAND/2; j++) {
    if(stest[j] != 1) {
      printf("There is a problem in the boundary of the even odd geometry\n");
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

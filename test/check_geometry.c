/*******************************************************************************
 * $Id$
 *
 * File check_geometry.c
 *
 * Consistency of the index arrays ipt, iup and idn
 *
 * Author: Carsten Urbach <urbach@physik.fu-berlin.de>
 *         using a file of Martin Luescher as template
 *
 *******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#ifdef MPI
#include <mpi.h>
#endif
#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include "global.h"
#include "geometry_eo.h"
#include "test/check_geometry.h"
#ifdef MPI
#include "mpi_init.h"
#endif

int check_geometry()
{
  int ix, j;
  int * stest;
  int * itest;
  int x0,x1,x2,x3;
  int iy0,iy1,iy2,iy3;
  int iz0,iz1,iz2,iz3;   
  
  itest = calloc(VOLUMEPLUSRAND + g_dbw2rand, sizeof(int));
  stest = calloc((VOLUMEPLUSRAND)/2, sizeof(int));

  for (ix=0;ix<VOLUMEPLUSRAND + g_dbw2rand;ix++){
    itest[ix]=0;
  }


  
  for (x0 = 0; x0 < T; x0++){
    for (x1 = 0; x1 < LX; x1++){
      for (x2 = 0; x2 < LY; x2++){
	for (x3 = 0; x3 < LZ; x3++){
	  ix=g_ipt[x0][x1][x2][x3];

	  if ((ix < 0) || (ix >= VOLUME)) {
	    printf("The index ipt is out of range (%d, %d, %d, %d) ix = %d\n", x0, x1, x2, x3, ix);
	    printf("Program aborted\n");
#ifdef MPI
	    MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
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
      printf("The index ipt is not one-to-one %d\n", itest[ix]);
      printf("Program aborted\n");
#ifdef MPI
      MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
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
	    iz0 = g_ipt[T][x1][x2][x3];
	    itest[iy0]++;
	  }
#else
	  iz0=g_ipt[(x0+1)%T][x1][x2][x3];
#endif     
                  
	  iy1=g_iup[ix][1];
#if (defined PARALLELXT)
	  if(x1 !=LX-1) {
	    iz1=g_ipt[x0][(x1+1)%LX][x2][x3];
	  }
	  else {
	    iz1=g_ipt[x0][LX][x2][x3];
	    itest[iy1]++;
	    if(iy1 < VOLUME + 2*LX*LY*LZ || iy1 >= VOLUME + 2*LX*LY*LZ + T*LY*LZ) {
	      printf("Boundary for x direction up is wrong %d %d %d\n", iy1 < VOLUME, iy1 > VOLUMEPLUSRAND, iy1);
	      printf("Program aborted\n");
#ifdef MPI
	      MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
#endif
	      exit(0);
	    }
	  }
#else
	  iz1=g_ipt[x0][(x1+1)%LX][x2][x3];
#endif

	  iy2=g_iup[ix][2];
	  iz2=g_ipt[x0][x1][(x2+1)%LY][x3];

	  iy3=g_iup[ix][3];
	  iz3=g_ipt[x0][x1][x2][(x3+1)%LZ];               
               
	  if ((iy0!=iz0)||(iy1!=iz1)||(iy2!=iz2)||(iy3!=iz3)||
	      (g_idn[iy0][0]!=ix)||(g_idn[iy1][1]!=ix)||(g_idn[iy2][2]!=ix)||(g_idn[iy3][3]!=ix)) {
	    printf("The index iup is incorrect\n");
	    printf("%d %d %d %d\n", (iy0!=iz0), (iy1!=iz1), (iy2!=iz2), (iy3!=iz3));
	    printf("Program aborted\n");
#ifdef MPI
	    MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
#endif
	    exit(0);
	  }

	  iy0=g_idn[ix][0];
#if (defined PARALLELT || defined PARALLELXT )
	  if(x0 !=0) {
	    iz0=g_ipt[(x0+T-1)%T][x1][x2][x3];
	  }
	  else {
	    iz0 = g_ipt[T+1][x1][x2][x3];;
	    itest[iy0]++;
	  }
#else
	  iz0=g_ipt[(x0+T-1)%T][x1][x2][x3];
#endif
	  iy1=g_idn[ix][1];
	  if(g_iup[iy1][1]!=ix) {
	    printf("Hallo\n");
	  }
#if (defined PARALLELXT)
	  if(x1 !=0) {
	    iz1=g_ipt[x0][(x1+LX-1)%LX][x2][x3];
	  }
	  else {
	    iz1 = g_ipt[x0][LX+1][x2][x3];;
	    itest[iy1]++;
	    if(iy1 < VOLUME + 2*LX*LY*LZ + T*LY*LZ || iy1 > VOLUME + RAND) {
	      printf("Boundary for x direction is wrong %d %d %d\n", iy1 < VOLUME, iy1 > VOLUMEPLUSRAND, iy1);
	      printf("Program aborted\n");
#ifdef MPI
	      MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
#endif
	      exit(0);
	    }
	  }
#else
	  iz1=g_ipt[x0][(x1+LX-1)%LX][x2][x3];
#endif
	  iy2=g_idn[ix][2];
	  iz2=g_ipt[x0][x1][(x2+LY-1)%LY][x3];

	  iy3=g_idn[ix][3];
	  iz3=g_ipt[x0][x1][x2][(x3+LZ-1)%LZ];               
               
	  if ((iy0!=iz0)||(iy1!=iz1)||(iy2!=iz2)||(iy3!=iz3)||
	      (g_iup[iy0][0]!=ix)||(g_iup[iy1][1]!=ix)||(g_iup[iy2][2]!=ix)||(g_iup[iy3][3]!=ix)) {
	    printf("The index idn is incorrect\n");
	    printf("%d %d %d %d\n", (iy0!=iz0), (iy1!=iz1), (iy2!=iz2), (iy3!=iz3));
	    printf("Program aborted\n");
#ifdef MPI
	    MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
#endif
	    exit(0);
	  }
	}
      }
    }
  }
  
  for (ix = VOLUME; ix < (VOLUME+RAND); ix++){
    if (itest[ix]!=1) {
      printf("The boundary is not correctly used itest = %d ix = %d \n", itest[ix], ix);
      printf("Program aborted\n");
#ifdef MPI
      MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
#endif
      exit(0); 
    }
  }

#ifdef MPI
  if(g_dbw2rand > 0) {
    for (x1 = 0; x1 < LX; x1++) {
      for (x2 = 0; x2 < LY; x2++) {
	for (x3 = 0; x3 < LZ; x3++) {
	  x0 = T;
	  ix=g_ipt[x0][x1][x2][x3];

	  iy0=g_iup[ix][0];
	  if(iy0 < VOLUMEPLUSRAND || iy0 > VOLUMEPLUSRAND+g_dbw2rand) {
	    printf("The DBW2 boundary is not correctly mapped in up t-direction  %d %d %d %d %d %d\n", x0, x1, x2, x3, ix, iy0);
	    printf("Program aborted\n");
	    MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	    exit(0); 
	  }
	  itest[iy0]++;
	  if (itest[iy0]>1) {
	    printf("The DBW2 boundary is not correctly used itest = %d (%d %d %d %d) iy0 = %d ix = %d\n", itest[iy0], x0, x1, x2, x3, iy0, ix);
	    printf("Program aborted\n");
	    MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	    exit(0); 
	  }
	  
	  x0 = T+1;
	  ix=g_ipt[x0][x1][x2][x3];

	  iy0=g_idn[ix][0];
	  if(iy0 < VOLUMEPLUSRAND || iy0 > VOLUMEPLUSRAND+g_dbw2rand) {
	    printf("The DBW2 boundary is not correctly mapped in down t-direction %d %d %d %d %d %d\n", x0, x1, x2, x3, ix, iy0);
	    printf("Program aborted\n");
	    MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	    exit(0); 
	  }
	  itest[iy0]++;
	  if (itest[iy0]>1) {
	    printf("The DBW2 boundary is not correctly used itest = %d (%d %d %d %d) iy0 = %d ix = %d\n", itest[iy0], x0, x1, x2, x3, iy0, ix);
	    printf("Program aborted\n");
	    MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	    exit(0); 
	  }
	}
      }
    }
#ifdef PARALLELXT
    for (x0 = 0; x0 < T+2; x0++) {
      for (x2 = 0; x2 < LY; x2++) {
	for (x3 = 0; x3 < LZ; x3++) {
	  x1 = LX;
	  ix = g_ipt[x0][x1][x2][x3];

	  iy1=g_iup[ix][1];
	  if(iy1 < VOLUMEPLUSRAND || iy1 > VOLUMEPLUSRAND+g_dbw2rand) {
	    printf("The DBW2 boundary is not correctly mapped in up x-direction %d %d %d %d %d %d\n", x0, x1, x2, x3, ix, iy1);
	    printf("Program aborted\n");
	    MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	    exit(0); 
	  }
	  itest[iy1]++;
	  if (itest[iy1]>1) {
	    printf("The DBW2 boundary is not correctly used itest = %d (%d %d %d %d) iy1 = %d ix = %d \n", itest[iy1], x0, x1, x2, x3, iy1, ix);
	    printf("Program aborted\n");
	    MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	    exit(0); 
	  }

	  if(x0 == T) {
	    iy0 = g_iup[ix][0];
	    if(iy0 < VOLUMEPLUSRAND || iy0 > VOLUMEPLUSRAND+g_dbw2rand) {
	      printf("The DBW2 boundary is not correctly mapped in down t-direction up x %d %d %d %d %d %d\n", x0, x1, x2, x3, ix, iy0);
	      printf("Program aborted\n");
	      MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	      exit(0); 
	    }
	    itest[iy0]++;
	    if (itest[iy0]>1) {
	    printf("The DBW2 boundary is not correctly used itest = %d (%d %d %d %d) iy0 = %d ix = %d\n", itest[iy0], x0, x1, x2, x3, iy0, ix);
	      printf("Program aborted\n");
	      MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	      exit(0); 
	    }
	  }
	  if(x0 == T+1) {
	    iy0 = g_idn[ix][0];
	    if(iy0 < VOLUMEPLUSRAND || iy0 > VOLUMEPLUSRAND+g_dbw2rand) {
	      printf("The DBW2 boundary is not correctly mapped in down t-direction up x %d %d %d %d %d %d\n", x0, x1, x2, x3, ix, iy0);
	      printf("Program aborted\n");
	      MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	      exit(0); 
	    }
	    itest[iy0]++;
	    if (itest[iy0]>1) {
	    printf("The DBW2 boundary is not correctly used itest = %d (%d %d %d %d) iy0 = %d ix = %d\n", itest[iy0], x0, x1, x2, x3, iy0, ix);
	      printf("Program aborted\n");
	      MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	      exit(0); 
	    }
	  }


	  x1 = LX+1;
	  ix = g_ipt[x0][x1][x2][x3];

	  iy1=g_idn[ix][1];
	  if(iy1 < VOLUMEPLUSRAND || iy1 > VOLUMEPLUSRAND+g_dbw2rand) {
	    printf("The DBW2 boundary is not correctly mapped in down x-direction %d %d %d %d %d %d\n", x0, x1, x2, x3, ix, iy1);
	    printf("Program aborted\n");
	    MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	    exit(0); 
	  }
	  itest[iy1]++;
	  if (itest[iy1]>1) {
	    printf("The DBW2 boundary is not correctly used itest = %d (%d %d %d %d) iy1 = %d ix = %d \n", itest[iy1], x0, x1, x2, x3, iy1, ix);
	    printf("Program aborted\n");
	    MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	    exit(0); 
	  }

	  if(x0 == T) {
	    iy0 = g_iup[ix][0];
	    if(iy0 < VOLUMEPLUSRAND || iy0 > VOLUMEPLUSRAND+g_dbw2rand) {
	      printf("The DBW2 boundary is not correctly mapped in down t-direction down x %d %d %d %d %d %d\n", x0, x1, x2, x3, ix, iy0);
	      printf("Program aborted\n");
	      MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	      exit(0); 
	    }
	    itest[iy0]++;
	    if (itest[iy0]>1) {
	      printf("The DBW2 boundary is not correctly used itest = %d (%d %d %d %d) iy0 = %d ix = %d \n", itest[iy0], x0, x1, x2, x3, iy0, ix);
	      printf("Program aborted\n");
	      MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	      exit(0); 
	    }
	  }
	  if(x0 == T+1) {
	    iy0 = g_idn[ix][0];
	    if(iy0 < VOLUMEPLUSRAND || iy0 > VOLUMEPLUSRAND+g_dbw2rand) {
	      printf("The DBW2 boundary is not correctly mapped in down t-direction down x %d %d %d %d %d %d\n", x0, x1, x2, x3, ix, iy0);
	      printf("Program aborted\n");
	      MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	      exit(0); 
	    }
	    itest[iy0]++;
	    if (itest[iy0]>1) {
	      printf("The DBW2 boundary is not correctly used itest = %d (%d %d %d %d) iy0 = %d ix = %d \n", itest[iy0], x0, x1, x2, x3, iy0, ix);
	      printf("Program aborted\n");
	      MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	      exit(0); 
	    }
	  }
	}
      }
    }
#endif
  }

  for (ix = VOLUMEPLUSRAND; ix < (VOLUMEPLUSRAND)+g_dbw2rand; ix++){
    if (itest[ix]!=1) {
      printf("The DBW2 boundary is not correctly used itest = %d ix = %d \n", itest[ix], ix);
      printf("Program aborted\n");
#ifdef MPI
      MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
#endif
      exit(0); 
    }
  }
#endif

  /* check of EO geometry */

  for (ix=0;ix<VOLUMEPLUSRAND/2;ix++){
    stest[ix]=0;
  }
  
  for(j = 0; j < VOLUME/2; j++) {
    ix = g_eo2lexic[j];
    
    iy0 = g_idn[ix][0];
    iz0 = g_lexic2eosub[iy0];
    if(iz0 > VOLUMEPLUSRAND/2 || iz0 < 0) {
      printf("There is a problem with EO geometry in direction 0-\n");
      printf("%d\n", iz0);
#ifdef MPI
      MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
#endif
      exit(0);
    }
    stest[iz0] += 1;
    
    iy0 = g_iup[ix][0];
    iz0 = g_lexic2eosub[iy0];
    if(iz0 > VOLUMEPLUSRAND/2 || iz0 < 0) {
      printf("There is a problem with EO geometry in direction 0+\n");
#ifdef MPI
      MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
#endif
      exit(0);
    }
    stest[iz0] += 1;
    
    iy1 = g_idn[ix][1];
    iz1 = g_lexic2eosub[iy1];
    if(iz1 > VOLUMEPLUSRAND/2  || iz1 < 0) {
      printf("There is a problem with EO geometry in direction 1-\n");
#ifdef MPI
      MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
#endif
      exit(0);
    }
    stest[iz1] += 1;
    
    iy1 = g_iup[ix][1];
    iz1 = g_lexic2eosub[iy1];
    if(iz1 >= VOLUMEPLUSRAND/2  || iz1 < 0) {
      printf("There is a problem with EO geometry in direction 1+\n");
#ifdef MPI
      MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
#endif
      exit(0);
    }
    stest[iz1] += 1;
    
    iy2 = g_idn[ix][2];
    iz2 = g_lexic2eosub[iy2];
    if(iz2 > VOLUMEPLUSRAND/2 || iz2 < 0) {
      printf("There is a problem with EO geometry in direction 2-\n");
#ifdef MPI
      MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
#endif
      exit(0);
    }
    stest[iz2] += 1;
    
    iy2 = g_iup[ix][2];
    iz2 = g_lexic2eosub[iy2];
    if(iz2 > VOLUMEPLUSRAND/2 || iz2 < 0) {
      printf("There is a problem with EO geometry in direction 2+\n");
#ifdef MPI
      MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
#endif
      exit(0);
    }
    stest[iz2] += 1;
    
    
    iy3 = g_idn[ix][3];
    iz3 = g_lexic2eosub[iy3];
    if(iz3 > VOLUMEPLUSRAND/2 || iz3 < 0) {
      printf("There is a problem with EO geometry in direction 3-\n");
#ifdef MPI
      MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
#endif
      exit(0);
    }
    stest[iz3] += 1;
    
    iy3 = g_iup[ix][3];
    iz3 = g_lexic2eosub[iy3];
    if(iz3 > VOLUMEPLUSRAND/2 || iz3 < 0) {
      printf("There is a problem with EO geometry in direction 3+\n");
#ifdef MPI
      MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
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
    printf("There is a problem in the first part of the even odd geometry\n");
    printf("%d is not equal to 8*(VOLUME)/2-RAND/2=%d\n", iz0, 8*(VOLUME)/2-RAND/2);
#ifdef MPI
    MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
#endif
    exit(0); 
  }

  for(j = VOLUME/2; j < (VOLUME+RAND)/2; j++) {
    if(stest[j] != 1) {
      printf("There is a problem in the first boundary of the even odd geometry\n");
#ifdef MPI
      MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
#endif
      exit(0);
    }
  }
  

  for (ix=0;ix<VOLUMEPLUSRAND/2;ix++){
    stest[ix]=0;
  }

  for(j = (VOLUME+RAND)/2; j < VOLUME+RAND/2; j++) {
    ix = g_eo2lexic[j];
    
    iy0 = g_idn[ix][0];
    iz0 = g_lexic2eosub[iy0];
    if(iz0 > VOLUMEPLUSRAND/2 || iz0 < 0) {
      printf("There is a problem with EO geometry in direction 0-\n");
#ifdef MPI
      MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
#endif
      exit(0);
    }
    stest[iz0] += 1;
    
    iy0 = g_iup[ix][0];
    iz0 = g_lexic2eosub[iy0];
    if(iz0 > VOLUMEPLUSRAND/2 || iz0 < 0) {
      printf("There is a problem with EO geometry in direction 0+\n");
#ifdef MPI
      MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
#endif
      exit(0);
    }
    stest[iz0] += 1;
    
    iy1 = g_idn[ix][1];
    iz1 = g_lexic2eosub[iy1];
    if(iz1 > VOLUMEPLUSRAND/2  || iz1 < 0) {
      printf("There is a problem with EO geometry in direction 1-\n");
#ifdef MPI
      MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
#endif
      exit(0);
    }
    stest[iz1] += 1;
    
    iy1 = g_iup[ix][1];
    iz1 = g_lexic2eosub[iy1];
    if(iz1 > VOLUMEPLUSRAND/2  || iz1 < 0) {
      printf("There is a problem with EO geometry in direction 1+\n"); 
#ifdef MPI
      MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize(); 
#endif
      exit(0); 
    }
    stest[iz1] += 1; 
    
    iy2 = g_idn[ix][2];
    iz2 = g_lexic2eosub[iy2];
    if(iz2 > VOLUMEPLUSRAND/2 || iz2 < 0) {
      printf("There is a problem with EO geometry in direction 2-\n");
#ifdef MPI
      MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
#endif
      exit(0);
    }
    stest[iz2] += 1;
    
    iy2 = g_iup[ix][2];
    iz2 = g_lexic2eosub[iy2];
    if(iz2 > VOLUMEPLUSRAND/2 || iz2 < 0) {
      printf("There is a problem with EO geometry in direction 2+\n");
#ifdef MPI
      MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
#endif
      exit(0);
    }
    stest[iz2] += 1;
    
    
    iy3 = g_idn[ix][3];
    iz3 = g_lexic2eosub[iy3];
    if(iz3 > VOLUMEPLUSRAND/2 || iz3 < 0) {
      printf("There is a problem with EO geometry in direction 3-\n");
#ifdef MPI
      MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
#endif
      exit(0);
    }
    stest[iz3] += 1;
    
    iy3 = g_iup[ix][3];
    iz3 = g_lexic2eosub[iy3];
    if(iz3 > VOLUMEPLUSRAND/2 || iz3 < 0) {
      printf("There is a problem with EO geometry in direction 3+\n");
#ifdef MPI
      MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
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
    printf("There is a problem in the second part of the even odd geometry\n");
    printf("%d is not equal to 8*(VOLUME)/2-RAND/2=%d\n", iz0, 8*(VOLUME)/2-RAND/2);
#ifdef MPI
    MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
#endif
    exit(0); 
  }

  for(j = VOLUME/2; j < (VOLUME+RAND)/2; j++) {
    if(stest[j] != 1) {
      printf("There is a problem in the second boundary of the even odd geometry\n");
#ifdef MPI
      MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
#endif
      exit(0);
    }
  }

  if(g_proc_id == 0 ) {
    printf("The lattice is correctly mapped by the index arrays\n");
    printf("\n"); 
  }

  free(stest);
  free(itest);

  return(0);
}

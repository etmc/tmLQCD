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
#include "global.h"
#include "geometry_eo.h"

int main(void)
{
  int ix;
  int * itest;
  int x0,x1,x2,x3;
  int iy0,iy1,iy2,iy3;
  int iz0,iz1,iz2,iz3;   

  itest = calloc(VOLUME, sizeof(int));
  printf("\n");
  printf("Test of the geometry programs \n");
  printf("So far only for serial version\n");
  printf("------------------------------\n");

  printf("\n");
  printf("The lattice size is %d x %d^3 \n\n",(int)(T),(int)(L));
   
  geometry();

  for (ix=0;ix<VOLUME;ix++){
    itest[ix]=0;
  }
   
  for (x0=0;x0<T;x0++){
    for (x1=0;x1<L;x1++){
      for (x2=0;x2<L;x2++){
	for (x3=0;x3<L;x3++){
	  ix=g_ipt[x0][x1][x2][x3];

	  if ((ix<0)||(ix>=VOLUME)){
	    printf("The index ipt is out of range\n");
	    printf("Program aborted\n");
	    exit(0);
	  }
               
	  itest[ix]+=1;
	}
      }
    }
  }

  for (ix=0;ix<VOLUME;ix++){
    if (itest[ix]!=1){
      printf("The index ipt is not one-to-one\n");
      printf("Program aborted\n");
      exit(0);
    }
  }

  for (x0=0;x0<T;x0++){
    for (x1=0;x1<L;x1++){
      for (x2=0;x2<L;x2++){
	for (x3=0;x3<L;x3++){
	  ix=g_ipt[x0][x1][x2][x3];

	  iy0=g_iup[ix][0];
	  iz0=g_ipt[(x0+1)%T][x1][x2][x3];
                       
	  iy1=g_iup[ix][1];
	  iz1=g_ipt[x0][(x1+1)%L][x2][x3];

	  iy2=g_iup[ix][2];
	  iz2=g_ipt[x0][x1][(x2+1)%L][x3];

	  iy3=g_iup[ix][3];
	  iz3=g_ipt[x0][x1][x2][(x3+1)%L];               
               
	  if ((iy0!=iz0)||(iy1!=iz1)||(iy2!=iz2)||(iy3!=iz3)){
	    printf("The index iup is incorrect\n");
	    printf("Program aborted\n");
	    exit(0);
	  }

	  iy0=g_idn[ix][0];
	  iz0=g_ipt[(x0+T-1)%T][x1][x2][x3];
                       
	  iy1=g_idn[ix][1];
	  iz1=g_ipt[x0][(x1+L-1)%L][x2][x3];

	  iy2=g_idn[ix][2];
	  iz2=g_ipt[x0][x1][(x2+L-1)%L][x3];

	  iy3=g_idn[ix][3];
	  iz3=g_ipt[x0][x1][x2][(x3+L-1)%L];               
               
	  if ((iy0!=iz0)||(iy1!=iz1)||(iy2!=iz2)||(iy3!=iz3)){
	    printf("The index idn is incorrect\n");
	    printf("Program aborted\n");
	    exit(0);
	  }
	}
      }
    }
  }

  printf("The lattice is correctly mapped by the index arrays\n");
  printf("\n");
  return(0);
}

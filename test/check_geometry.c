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
 * File check_geometry.c
 *
 * Consistency of the index arrays ipt, iup and idn
 *
 * Author: Carsten Urbach <urbach@physik.fu-berlin.de>
 *         using a file of Martin Luescher as template
 *
 *******************************************************************************/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#ifdef MPI
#include <mpi.h>
#endif
#include "global.h"
#include "geometry_eo.h"
#include "test/check_geometry.h"
#ifdef MPI
#include "mpi_init.h"
#endif

#if defined _INDEX_INDEP_GEOM

int check_geometry()
{
#ifdef XLC
#pragma execution_frequency(very_low)
#endif
  int ix, j;
  int * stest;
  int * itest;
  int x0,x1,x2,x3;
  int iy0,iy1,iy2,iy3;
  int iz0,iz1,iz2,iz3;   
  int bndcnt = 0;
  int ext_t=0, ext_x=0, ext_y=0, ext_z=0;

#if ( defined PARALLELT || defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT )
  ext_t=2;
#endif
#if ( defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT || defined PARALLELX || defined PARALLELXY || defined PARALLELXYZ )
  ext_x=2;
#endif
#if ( defined PARALLELXYT || defined PARALLELXYZT || defined PARALLELXY || defined PARALLELXYZ )
  ext_y=2;
#endif
#if ( defined PARALLELXYZT || defined PARALLELXYZ )
  ext_z=2;
#endif




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
            return(-1);
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
      return(-1);
    }
  }

  for (x0 = 0; x0 < T; x0++){
    for (x1 = 0; x1 < LX; x1++){
      for (x2 = 0; x2 < LY; x2++){
	for (x3 = 0; x3 < LZ; x3++){
	  ix=g_ipt[x0][x1][x2][x3];

	  iy0=g_iup[ix][0];
#if (defined PARALLELT || defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT)
	  if(x0!=T-1) {
	    iz0=g_ipt[(x0+1)%T][x1][x2][x3];
	  }
	  else {
	    iz0 = g_ipt[T][x1][x2][x3];
	    itest[iy0]++;
	    if(iy0 <  gI_L_0_0_0  || iy0 >= gI_L_0_0_0 + LX*LY*LZ) {
	      printf("Boundary for time direction up is wrong %d %d %d\n", 
		     iy0 < gI_L_0_0_0, iy0 >= gI_L_0_0_0 + LX*LY*LZ, iy0);
          return(-1);
	    }
	  }
#else
	  iz0=g_ipt[(x0+1)%T][x1][x2][x3];
#endif     
                  
	  iy1=g_iup[ix][1];
#if (defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT || defined PARALLELX || defined PARALLELXY || defined PARALLELXYZ)
	  if(x1 !=LX-1) {
	    iz1=g_ipt[x0][(x1+1)%LX][x2][x3];
	  }
	  else {
	    iz1=g_ipt[x0][LX][x2][x3];
	    itest[iy1]++;
	    if(iy1 < gI_0_L_0_0 || iy1 >= gI_0_L_0_0 + T*LY*LZ) {
	      printf("Boundary for x direction up is wrong %d %d %d\n", 
		     iy1 < gI_0_L_0_0, iy1 >= gI_0_L_0_0 + T*LY*LZ, iy1);
            return(-1);
	    }
	  }
#else
	  iz1=g_ipt[x0][(x1+1)%LX][x2][x3];
#endif

	  iy2=g_iup[ix][2];
#if (defined PARALLELXYT || defined PARALLELXYZT || defined PARALLELXY || defined PARALLELXYZ )
	  if(x2 !=LY-1) {
	    iz2=g_ipt[x0][x1][(x2+1)%LY][x3];
	  }
	  else {
	    iz2=g_ipt[x0][x1][LY][x3];
	    itest[iy2]++;
	    if(iy2 < gI_0_0_L_0 || iy2 >= gI_0_0_L_0 + T*LX*LZ) {
	      printf("Boundary for y direction up is wrong %d %d %d\n", 
		     iy2 < gI_0_0_L_0, iy2 > gI_0_0_L_0 + T*LX*LZ, iy2);
            return(-1);
	    }
	  }
#else
	  iz2=g_ipt[x0][x1][(x2+1)%LY][x3];
#endif

	  iy3=g_iup[ix][3];
#if ( defined PARALLELXYZT || defined PARALLELXYZ )
	  if(x3 !=LZ-1) {
	    iz3=g_ipt[x0][x1][x2][(x3+1)%LZ];
	  }
	  else {
	    iz3=g_ipt[x0][x1][x2][LZ];
	    itest[iy3]++;
	    if(iy3 < gI_0_0_0_L || iy3 >= gI_0_0_0_L + T*LX*LY) {
	      printf("Boundary for z direction up is wrong %d %d %d\n", 
		     iy3 < gI_0_0_0_L, iy3 > gI_0_0_0_L+ T*LX*LY, iy3);
            return(-1);
	    }
	  }
#else
	  iz3=g_ipt[x0][x1][x2][(x3+1)%LZ];
#endif
               
	  if ((iy0!=iz0)||(iy1!=iz1)||(iy2!=iz2)||(iy3!=iz3)||
	      (g_idn[iy0][0]!=ix)||(g_idn[iy1][1]!=ix)||(g_idn[iy2][2]!=ix)||(g_idn[iy3][3]!=ix)) {
	    printf("The index iup is incorrect\n");
	    printf("%d %d %d %d\n", (iy0!=iz0), (iy1!=iz1), (iy2!=iz2), (iy3!=iz3));
	    printf("%d %d %d %d %d %d\n", x0, x1, x2, x3, iy1, iz1);
            return(-1);
	  }

	  iy0=g_idn[ix][0];
#if (defined PARALLELT || defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT)
	  if(x0 !=0) {
	    iz0=g_ipt[(x0+T-1)%T][x1][x2][x3];
	  }
	  else {
	    iz0 = g_ipt[T+1][x1][x2][x3];;
	    itest[iy0]++;
	    if(iy0 < gI_m1_0_0_0  || iy0 >= gI_m1_0_0_0  + LX*LY*LZ) {
	      printf("Boundary for time direction is wrong %d %d %d\n", 
		     iy0 < gI_m1_0_0_0, iy0 >= gI_m1_0_0_0 + LX*LY*LZ, iy0);
            return(-1);
	    }
	  }
#else
	  iz0=g_ipt[(x0+T-1)%T][x1][x2][x3];
#endif

	  iy1=g_idn[ix][1];
#if (defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT || defined PARALLELX || defined PARALLELXY || defined PARALLELXYZ)
	  if(x1 !=0) {
	    iz1=g_ipt[x0][(x1+LX-1)%LX][x2][x3];
	  }
	  else {
	    iz1 = g_ipt[x0][LX+1][x2][x3];
	    itest[iy1]++;
	    if(iy1 < gI_0_m1_0_0 || iy1 >= gI_0_m1_0_0 + T*LY*LZ) {
	      printf("Boundary for x direction is wrong %d %d %d\n", 
		     iy1 < gI_0_m1_0_0, iy1 >= gI_0_m1_0_0 + T*LY*LZ, iy1);
            return(-1);
	    }
	  }
#else
	  iz1=g_ipt[x0][(x1+LX-1)%LX][x2][x3];
#endif
	  iy2=g_idn[ix][2];
#if (defined PARALLELXYT || defined PARALLELXYZT || defined PARALLELXY || defined PARALLELXYZ )
	  if(x2 !=0) {
	    iz2=g_ipt[x0][x1][(x2+LY-1)%LY][x3];
	  }
	  else {
	    iz2 = g_ipt[x0][x1][LY+1][x3];
	    itest[iy2]++;
	    if(iy2 < gI_0_0_m1_0 || iy2 >= gI_0_0_m1_0 + T*LX*LZ) {
	      printf("Boundary for y direction is wrong %d %d %d\n", 
		     iy2 < gI_0_0_m1_0, iy2 >= gI_0_0_m1_0 + T*LX*LZ, iy2);
            return(-1);
	    }
	  }
#else
	  iz2=g_ipt[x0][x1][(x2+LY-1)%LY][x3];
#endif

	  iy3=g_idn[ix][3];
#if ( defined PARALLELXYZT || defined PARALLELXYZ )
	  if(x3 !=0) {
	    iz3=g_ipt[x0][x1][x2][(x3+LZ-1)%LZ];
	  }
	  else {
	    iz3 = g_ipt[x0][x1][x2][LZ+1];
	    itest[iy3]++;
	    if(iy3 < gI_0_0_0_m1 || iy3 >= gI_0_0_0_m1 + T*LX*LY) {
	      printf("Boundary for z direction is wrong %d %d %d\n", 
		     iy3 < gI_0_0_0_m1, iy3 >= gI_0_0_0_m1 + T*LX*LY, iy3);
	      printf("%d %d %d %d %d\n", x0, x1, x2, x3, ix);
            return(-1);
	    }
	  }
#else
	  iz3=g_ipt[x0][x1][x2][(x3+LZ-1)%LZ];
#endif
               
	  if ((iy0!=iz0)||(iy1!=iz1)||(iy2!=iz2)||(iy3!=iz3)||
	      (g_iup[iy0][0]!=ix)||(g_iup[iy1][1]!=ix)||(g_iup[iy2][2]!=ix)||(g_iup[iy3][3]!=ix)) {
	    printf("The index idn is incorrect\n");
	    printf("%d %d %d %d\n", (iy0!=iz0), (iy1!=iz1), (iy2!=iz2), (iy3!=iz3));
	    printf("%d %d %d %d\n", iy1, iz1, iy2, iz2);
	    printf("%d %d %d %d %d\n", x0, x1, x2, x3, ix);
            return(-1);
	  }

  /* The edges */
  /* In case of PARALLELT or PARALLELX there is actually no edge to take care of */
#if ((defined PARALLELXT) || (defined PARALLELXYT) || (defined PARALLELXYZT))
	  if(x0 == 0) {
	    iy0 = g_idn[ g_idn[ix][1] ][0];
	    if(x1 != 0) {
	      iz0 = g_ipt[T+1][(x1+LX-1)%LX][x2][x3];
	    }
	    else {
	      iz0 = g_ipt[T+1][LX+1][x2][x3];
	      itest[iy0]++;
	    }
	    if(iz0 != iy0) {
	      printf("Edge -t -x has an error\n");
            return(-1);
	    }

	    iy0 = g_idn[ g_iup[ix][1] ][0]; 
	    if(x1 != LX-1) {
	      iz0 = g_ipt[T+1][(x1+1)%LX][x2][x3];
	    }
	    else {
	      iz0 = g_ipt[T+1][LX][x2][x3];
	      itest[iy0]++;
	    }
	    if(iz0 != iy0) {
	      printf("Edge -t +x has an error\n");
            return(-1);
	    }
	  }

	  if(x0 == T-1) {
	    iy0 = g_iup[ g_idn[ix][1] ][0];
	    if(x1 != 0) {
	      iz0 = g_ipt[T][(x1+LX-1)%LX][x2][x3];
	    }
	    else {
	      iz0 = g_ipt[T][LX+1][x2][x3];
	      itest[iy0]++;
	    }
	    if(iz0 != iy0) {
	      printf("Edge +t -x has an error\n");
            return(-1);
	    }

	    iy0 = g_iup[ g_iup[ix][1] ][0]; 
	    if(x1 != LX-1) {
	      iz0 = g_ipt[T][(x1+1)%LX][x2][x3];
	    }
	    else {
	      iz0 = g_ipt[T][LX][x2][x3];
	      itest[iy0]++;
	    }
	    if(iz0 != iy0) {
	      printf("Edge +t +x has an error\n");
            return(-1);
	    }
	  }

#endif

#if (defined PARALLELXYT || defined PARALLELXYZT)
	  if(x0 == 0) {
	    iy0 = g_idn[ g_idn[ix][2] ][0]; 
	    if(x2 != 0) {
	      iz0 = g_ipt[T+1][x1][(x2+LY-1)%LY][x3];
	    }
	    else {
	      iz0 = g_ipt[T+1][x1][LY+1][x3];
	      itest[iy0]++;
	    }
	    if(iz0 != iy0) {
	      printf("Edge -t -y has an error\n");
            return(-1);
	    }

	    iy0 = g_idn[ g_iup[ix][2] ][0]; 
	    if(x2 != LY-1) {
	      iz0 = g_ipt[T+1][x1][(x2+1)%LY][x3];
	    }
	    else {
	      iz0 = g_ipt[T+1][x1][LY][x3];
	      itest[iy0]++;
	    }
	    if(iz0 != iy0) {
	      printf("Edge -t +y has an error\n");
            return(-1);
	    }
	  }

	  if(x0 == T-1) {
	    iy0 = g_iup[ g_idn[ix][2] ][0];
	    if(x2 != 0) {
	      iz0 = g_ipt[T][x1][(x2+LY-1)%LY][x3];
	    }
	    else {
	      iz0 = g_ipt[T][x1][LY+1][x3];
	      itest[iy0]++;
	    }
	    if(iz0 != iy0) {
	      printf("Edge +t -y has an error\n");
            return(-1);
	    }

	    iy0 = g_iup[ g_iup[ix][2] ][0]; 
	    if(x2 != LY-1) {
	      iz0 = g_ipt[T][x1][(x2+1)%LY][x3];
	    }
	    else {
	      iz0 = g_ipt[T][x1][LY][x3];
	      itest[iy0]++;
	    }
	    if(iz0 != iy0) {
	      printf("Edge +t +y has an error\n");
            return(-1);
	    }
	  }

#endif
#if (defined PARALLELXYT || defined PARALLELXYZT || defined PARALLELXY || defined PARALLELXYZ)

	  if(x1 == 0) {
	    iy0 = g_idn[ g_idn[ix][2] ][1]; 
	    if(x2 != 0) {
	      iz0 = g_ipt[x0][LX+1][(x2+LY-1)%LY][x3];
	    }
	    else {
	      iz0 = g_ipt[x0][LX+1][LY+1][x3];
	      itest[iy0]++;
	    }
	    if(iz0 != iy0) {
	      printf("Edge -x -y has an error\n");
            return(-1);
	    }
	    iy0 = g_idn[ g_iup[ix][2] ][1]; 
	    if(x2 != LY-1) {
	      iz0 = g_ipt[x0][LX+1][(x2+1)%LY][x3];
	    }
	    else {
	      iz0 = g_ipt[x0][LX+1][LY][x3];
	      itest[iy0]++;
	    }
	    if(iz0 != iy0) {
	      printf("Edge -x +y has an error\n");
            return(-1);
	    }	    
	  }
	  if(x1 == LX-1) {
	    iy0 = g_iup[ g_idn[ix][2] ][1];
	    if(x2 != 0) {
	      iz0 = g_ipt[x0][LX][(x2+LY-1)%LY][x3];
	    }
	    else {
	      iz0 = g_ipt[x0][LX][LY+1][x3];
	      itest[iy0]++;
	    }
	    if(iz0 != iy0) {
	      printf("Edge +x -y has an error\n");
            return(-1);
	    }

	    iy0 = g_iup[ g_iup[ix][2] ][1]; 
	    if(x2 != LY-1) {
	      iz0 = g_ipt[x0][LX][(x2+1)%LY][x3];
	    }
	    else {
	      iz0 = g_ipt[x0][LX][LY][x3];
	      itest[iy0]++;
	    }
	    if(iz0 != iy0) {
	      printf("Edge +x +y has an error\n");
            return(-1);
	    }
	  }
#endif
#if defined PARALLELXYZT
	  if(x0 == 0) {
	    iy0 = g_idn[ g_idn[ix][3] ][0]; 
	    if(x3 != 0) {
	      iz0 = g_ipt[T+1][x1][x2][(x3+LZ-1)%LZ];
	    }
	    else {
	      iz0 = g_ipt[T+1][x1][x2][LZ+1];
	      itest[iy0]++;
	      if(itest[iy0]>1) {
		printf("Edge -t -z has itest = %d at %d, %d, %d, %d, iy0 = %d\n", itest[iy0], x0, x1, x2, x3, iy0);
            return(-1);
	      }
	    }
	    if(iz0 != iy0) {
	      printf("Edge -t -z has an error\n");
	      printf("Program aborted\n");
#  ifdef MPI
	      MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
#  endif
	      exit(0);
	    }

	    iy0 = g_idn[ g_iup[ix][3] ][0]; 
	    if(x3 != LZ-1) {
	      iz0 = g_ipt[T+1][x1][x2][(x3+1)%LZ];
	    }
	    else {
	      iz0 = g_ipt[T+1][x1][x2][LZ];
	      itest[iy0]++;
	      if(itest[iy0]>1) {
		printf("Edge -t +z has itest = %d at %d, %d, %d, %d, iy0 = %d\n", itest[iy0], x0, x1, x2, x3, iy0);
		printf("ix = %d, iz0 = %d %d\n", ix, iz0, g_iup[ix][3]); 
            return(-1);
	      }
	    }
	    if(iz0 != iy0) {
	      printf("Edge -t +z has an error\n");
	      printf("Program aborted\n");
#  ifdef MPI
	      MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
#  endif
	      exit(0);
	    }
	  }

	  if(x0 == T-1) {
	    iy0 = g_iup[ g_idn[ix][3] ][0];
	    if(x3 != 0) {
	      iz0 = g_ipt[T][x1][x2][(x3+LZ-1)%LZ];
	    }
	    else {
	      iz0 = g_ipt[T][x1][x2][LZ+1];
	      itest[iy0]++;
	      if(itest[iy0]>1) {
		printf("Edge +t -z has itest = %d at %d, %d, %d, %d, iy0 = %d\n", itest[iy0], x0, x1, x2, x3, iy0);
            return(-1);
	      }
	    }
	    if(iz0 != iy0) {
	      printf("Edge +t -z has an error\n");
            return(-1);
	    }

	    iy0 = g_iup[ g_iup[ix][0] ][3]; 
	    if(x3 != LZ-1) {
	      iz0 = g_ipt[T][x1][x2][(x3+1)%LZ];
	    }
	    else {
	      iz0 = g_ipt[T][x1][x2][LZ];
	      itest[iy0]++;
	      if(itest[iy0]>1) {
		printf("Edge +t +z has itest = %d at %d, %d, %d, %d, iy0 = %d\n", itest[iy0], x0, x1, x2, x3, iy0);
            return(-1);
	      }
	    }
	    if(iz0 != iy0) {
	      printf("Edge +t +z has an error\n");
            return(-1);
	    }
	  }

#endif
#if ( defined PARALLELXYZT || defined PARALLELXYZ )

	  if(x1 == 0) {
	    iy0 = g_idn[ g_idn[ix][3] ][1]; 
	    if(x3 != 0) {
	      iz0 = g_ipt[x0][LX+1][x2][(x3+LZ-1)%LZ];
	    }
	    else {
	      iz0 = g_ipt[x0][LX+1][x2][LZ+1];
	      itest[iy0]++;
	      if(itest[iy0]>1) {
		printf("Edge -x -z has itest = %d at %d, %d, %d, %d, iy0 = %d\n", itest[iy0], x0, x1, x2, x3, iy0);
            return(-1);
	      }
	    }
	    if(iz0 != iy0) {
	      printf("Edge -x -z has an error\n");
            return(-1);
	    }
	    iy0 = g_idn[ g_iup[ix][3] ][1]; 
	    if(x3 != LZ-1) {
	      iz0 = g_ipt[x0][LX+1][x2][(x3+1)%LZ];
	    }
	    else {
	      iz0 = g_ipt[x0][LX+1][x2][LZ];
	      itest[iy0]++;
	      if(itest[iy0]>1) {
		printf("Edge -x +z has itest = %d at %d, %d, %d, %d, iy0 = %d\n", itest[iy0], x0, x1, x2, x3, iy0);
            return(-1);
	      }
	    }
	    if(iz0 != iy0) {
	      printf("Edge -x +z has an error\n");
            return(-1);
	    }	    
	  }
	  if(x1 == LX-1) {
	    iy0 = g_iup[ g_idn[ix][3] ][1];
	    if(x3 != 0) {
	      iz0 = g_ipt[x0][LX][x2][(x3+LZ-1)%LZ];
	    }
	    else {
	      iz0 = g_ipt[x0][LX][x2][LZ+1];
	      itest[iy0]++;
	      if(itest[iy0]>1) {
		printf("Edge +x -z has itest = %d at %d, %d, %d, %d, iy0 = %d\n", itest[iy0], x0, x1, x2, x3, iy0);
            return(-1);
	      }
	    }
	    if(iz0 != iy0) {
	      printf("Edge +x -z has an error\n");
            return(-1);
	    }

	    iy0 = g_iup[ g_iup[ix][3] ][1]; 
	    if(x3 != LZ-1) {
	      iz0 = g_ipt[x0][LX][x2][(x3+1)%LZ];
	    }
	    else {
	      iz0 = g_ipt[x0][LX][x2][LZ];
	      itest[iy0]++;
	      if(itest[iy0]>1) {
		printf("Edge +x +z has itest = %d at %d, %d, %d, %d, iy0 = %d iz0 = %d giup[%d][3] = %d ix = %d\n",
		       itest[iy0], x0, x1, x2, x3, iy0, iz0, ix, g_iup[ix][3], ix);
            return(-1);
	      }
	    }
	    if(iz0 != iy0) {
	      printf("Edge +x +z has an error\n");
            return(-1);
	    }
	  }

	  if(x2 == 0) {
	    iy0 = g_idn[ g_idn[ix][3] ][2]; 
	    if(x3 != 0) {
	      iz0 = g_ipt[x0][x1][LY+1][(x3+LZ-1)%LZ];
	    }
	    else {
	      iz0 = g_ipt[x0][x1][LY+1][LZ+1];
	      itest[iy0]++;
	      if(itest[iy0]>1) {
		printf("Edge -y -z has itest = %d at %d, %d, %d, %d, iy0 = %d\n", itest[iy0], x0, x1, x2, x3, iy0);
            return(-1);
	      }
	    }
	    if(iz0 != iy0) {
	      printf("Edge -y -z has an error\n");
            return(-1);
	    }
	    iy0 = g_idn[ g_iup[ix][3] ][2]; 
	    if(x3 != LZ-1) {
	      iz0 = g_ipt[x0][x1][LY+1][(x3+1)%LZ];
	    }
	    else {
	      iz0 = g_ipt[x0][x1][LY+1][LZ];
	      itest[iy0]++;
	      if(itest[iy0]>1) {
		printf("Edge -y +z has itest = %d at %d, %d, %d, %d, iy0 = %d\n", itest[iy0], x0, x1, x2, x3, iy0);
            return(-1);
	      }
	    }
	    if(iz0 != iy0) {
	      printf("Edge -y +z has an error\n");
            return(-1);
	    }	    
	  }
	  if(x2 == LY-1) {
	    iy0 = g_iup[ g_idn[ix][3] ][2];
	    if(x3 != 0) {
	      iz0 = g_ipt[x0][x1][LY][(x3+LZ-1)%LZ];
	    }
	    else {
	      iz0 = g_ipt[x0][x1][LY][LZ+1];
	      itest[iy0]++;
	      if(itest[iy0]>1) {
		printf("Edge +y -z has itest = %d at %d, %d, %d, %d, iy0 = %d\n", itest[iy0], x0, x1, x2, x3, iy0);
            return(-1);
	      }
	    }
	    if(iz0 != iy0) {
	      printf("Edge +y -z has an error\n");
            return(-1);
	    }

	    iy0 = g_iup[ g_iup[ix][3] ][2]; 
	    if(x3 != LZ-1) {
	      iz0 = g_ipt[x0][x1][LY][(x3+1)%LZ];
	    }
	    else {
	      iz0 = g_ipt[x0][x1][LY][LZ];
	      itest[iy0]++;
	      if(itest[iy0]>1) {
		printf("Edge +y +z has itest = %d at %d, %d, %d, %d, iy0 = %d\n", itest[iy0], x0, x1, x2, x3, iy0);
            return(-1);
	      }
	    }
	    if(iz0 != iy0) {
	      printf("Edge +y +z has an error\n");
            return(-1);
	    }
	  }
#endif
	}
      }
    }
  }
  
  for (ix = VOLUME; ix < (VOLUME+RAND+EDGES); ix++){
    if (itest[ix]!=1) {
      printf("The boundary is not correctly used itest = %d ix = %d %d %d %d\n", itest[ix], ix, VOLUME, RAND, EDGES);
            return(-1);
    }
  }

  for (ix=0;ix<VOLUMEPLUSRAND + g_dbw2rand;ix++){
    itest[ix]=0;
  }

#ifdef MPI
  if(g_dbw2rand > 0) {

#if (defined PARALLELT || defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT)
    for (x1 = 0; x1 < LX; x1++) {
      for (x2 = 0; x2 < LY; x2++) {
	for (x3 = 0; x3 < LZ; x3++) {
	  x0 = T;
	  ix=g_ipt[x0][x1][x2][x3];

	  iy0=g_iup[ix][0];
	  if(iy0 < gI_Lp1_0_0_0 || iy0 >= gI_Lp1_0_0_0 + LX*LY*LZ) {
	    printf("The DBW2 boundary is not correctly mapped in up t-direction  %d %d %d %d %d %d\n", x0, x1, x2, x3, ix, iy0);
            return(-1);
	  }
	  itest[iy0]++;
	  if (itest[iy0]>1) {
	    printf("The DBW2 boundary is not correctly used itest = %d (%d %d %d %d) iy0 = %d ix = %d\n", itest[iy0], x0, x1, x2, x3, iy0, ix);
            return(-1);
	  }
	  
	  x0 = T+1;
	  ix=g_ipt[x0][x1][x2][x3];

	  iy0=g_idn[ix][0];
	  if(iy0 < gI_m2_0_0_0 || iy0 >= gI_m2_0_0_0 + LX*LY*LZ) {
	    printf("The DBW2 boundary is not correctly mapped in down t-direction %d %d %d %d %d %d\n", x0, x1, x2, x3, ix, iy0);
            return(-1);
	  }
	  itest[iy0]++;
	  if (itest[iy0]>1) {
	    printf("The DBW2 boundary is not correctly used itest = %d (%d %d %d %d) iy0 = %d ix = %d\n", itest[iy0], x0, x1, x2, x3, iy0, ix);
            return(-1);
	  }
	}
      }
    }
#endif

#if (defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT || defined PARALLELX || defined PARALLELXY || defined PARALLELXYZ )
    for (x0 = 0; x0 < T+ext_t; x0++) {
      for (x2 = 0; x2 < LY; x2++) {
	for (x3 = 0; x3 < LZ; x3++) {
	  x1 = LX;
	  ix = g_ipt[x0][x1][x2][x3];

	  iy1=g_iup[ix][1];
	  if((iy1 < gI_0_Lp1_0_0 || iy1 >= gI_0_Lp1_0_0 + T*LY*LZ) && x0 < T) {
	    printf("The DBW2 boundary is not correctly mapped in up x-direction %d %d %d %d %d %d\n", x0, x1, x2, x3, ix, iy1);
            return(-1);
	  }
	  itest[iy1]++;
	  if (itest[iy1]>1) {
	    printf("The DBW2 boundary is not correctly used up x itest = %d (%d %d %d %d) iy1 = %d ix = %d \n", itest[iy1], x0, x1, x2, x3, iy1, ix);
            return(-1);
	  }

	  if(x0 == T) {
	    iy0 = g_iup[ix][0];
	    if(iy0 < gI_Lp1_L_0_0 || iy0 >= gI_Lp1_L_0_0 + LY*LZ) {
	      printf("The DBW2 boundary is not correctly mapped in up t-direction up x %d %d %d %d %d %d\n", x0, x1, x2, x3, ix, iy0);
            return(-1);
	    }
	    itest[iy0]++;
	    if (itest[iy0]>1) {
	    printf("The DBW2 boundary is not correctly used up t up x itest = %d (%d %d %d %d) iy0 = %d ix = %d\n", itest[iy0], x0, x1, x2, x3, iy0, ix);
            return(-1);
	    }
	  }
	  if(x0 == T+1) {
	    iy0 = g_idn[ix][0];
	    if(iy0 < gI_m2_L_0_0 || iy0 >= gI_m2_L_0_0 + LY*LZ) {
	      printf("The DBW2 boundary is not correctly mapped in down t-direction up x %d %d %d %d %d %d\n", x0, x1, x2, x3, ix, iy0);
            return(-1);
	    }
	    itest[iy0]++;
	    if (itest[iy0]>1) {
	    printf("The DBW2 boundary is not correctly used down t up x itest = %d (%d %d %d %d) iy0 = %d ix = %d\n", itest[iy0], x0, x1, x2, x3, iy0, ix);
            return(-1);
	    }
	  }


	  x1 = LX+1;
	  ix = g_ipt[x0][x1][x2][x3];

	  iy1=g_idn[ix][1];
	  if((iy1 < gI_0_m2_0_0 || iy1 >= gI_0_m2_0_0 + T*LY*LZ) && x0 < T) {
	    printf("The DBW2 boundary is not correctly mapped in down x-direction %d %d %d %d %d %d\n", x0, x1, x2, x3, ix, iy1);
            return(-1);
	  }
	  itest[iy1]++;
	  if (itest[iy1]>1) {
	    printf("The DBW2 boundary is not correctly used down x itest = %d (%d %d %d %d) iy1 = %d ix = %d \n", itest[iy1], x0, x1, x2, x3, iy1, ix);
            return(-1);
	  }

	  if(x0 == T) {
	    iy0 = g_iup[ix][0];
	    if(iy0 < gI_Lp1_m1_0_0 || iy0 >= gI_Lp1_m1_0_0 + LY*LZ) {
	      printf("The DBW2 boundary is not correctly mapped in up t-direction down x %d %d %d %d %d %d\n", x0, x1, x2, x3, ix, iy0);
            return(-1);
	    }
	    itest[iy0]++;
	    if (itest[iy0]>1) {
	      printf("The DBW2 boundary is not correctly used up t down x itest = %d (%d %d %d %d) iy0 = %d ix = %d \n", itest[iy0], x0, x1, x2, x3, iy0, ix);
return(-1);	    }
	  }
	  if(x0 == T+1) {
	    iy0 = g_idn[ix][0];
	    if(iy0 < gI_m2_m1_0_0 || iy0 >= gI_m2_m1_0_0 + LY*LZ) {
	      printf("The DBW2 boundary is not correctly mapped in down t-direction down x %d %d %d %d %d %d\n", x0, x1, x2, x3, ix, iy0);
return(-1);
		}
	    itest[iy0]++;
	    if (itest[iy0]>1) {
	      printf("The DBW2 boundary is not correctly used down t down xitest = %d (%d %d %d %d) iy0 = %d ix = %d \n", itest[iy0], x0, x1, x2, x3, iy0, ix);
return(-1);
	    }
	  }
	}
      }
    }
#endif

#if ( defined PARALLELXYT || defined PARALLELXYZT || defined PARALLELXY || defined PARALLELXYZ )

    for (x0 = 0; x0 < T+ext_t; x0++) {
      for (x1 = 0; x1 < LX+ext_x; x1++) {
	for (x3 = 0; x3 < LZ; x3++) {
	  if(x0 < T || x1 < LX) {
	    x2 = LY;
	    ix = g_ipt[x0][x1][x2][x3];
	    
	    iy2=g_iup[ix][2];
	    if((iy2 < gI_0_0_Lp1_0 || iy2 >= gI_0_0_Lp1_0 + T*LX*LZ) 
	       && x0 < T && x1 < LX) {
	      printf("The DBW2 boundary is not correctly mapped in up y-direction %d %d %d %d %d %d\n", 
		     x0, x1, x2, x3, ix, iy2);
return(-1);
	    }
	    itest[iy2]++;
	    if (itest[iy2]>1) {
	      printf("The DBW2 boundary is not correctly used up y itest = %d (%d %d %d %d) iy2 = %d ix = %d \n", 
		     itest[iy2], x0, x1, x2, x3, iy2, ix);
return(-1);
	    }
	    
	    if(x0 == T && x1 < LX) {
	      iy0 = g_iup[ix][0];
	      if(iy0 < gI_Lp1_0_L_0 || iy0 >= gI_Lp1_0_L_0 + LX*LZ) {
		printf("The DBW2 boundary is not correctly mapped in up t-direction up y %d %d %d %d %d %d\n", 
		       x0, x1, x2, x3, ix, iy0);
return(-1);
	      }
	      itest[iy0]++;
	      if (itest[iy0]>1) {
		printf("The DBW2 boundary is not correctly used up t up y itest = %d (%d %d %d %d) iy0 = %d ix = %d\n", 
		       itest[iy0], x0, x1, x2, x3, iy0, ix);
return(-1);
	      }
	    }
	    if(x0 == T+1 && x1 < LX) {
	      iy0 = g_idn[ix][0];
	      if(iy0 < gI_m2_0_L_0 || iy0 >= gI_m2_0_L_0 + LX*LZ) {
		printf("The DBW2 boundary is not correctly mapped in down t-direction up y %d %d %d %d %d %d\n", 
		       x0, x1, x2, x3, ix, iy0);
return(-1);
	      }
	      itest[iy0]++;
	      if (itest[iy0]>1) {
		printf("The DBW2 boundary is not correctly used down t up y itest = %d (%d %d %d %d) iy0 = %d ix = %d\n", 
		       itest[iy0], x0, x1, x2, x3, iy0, ix);
return(-1);
	      }
	    }
	    
	    if(x1 == LX && x0 < T) {
	      iy1 = g_iup[ix][1];
	      if(iy1 < gI_0_Lp1_L_0 || iy1 >= gI_0_Lp1_L_0 + T*LZ) {
		printf("The DBW2 boundary is not correctly mapped in up x-direction up y %d %d %d %d %d %d\n", 
		       x0, x1, x2, x3, ix, iy1);
return(-1);
	      }
	      itest[iy1]++;
	      if (itest[iy1]>1) {
		printf("The DBW2 boundary is not correctly used x up y up itest = %d (%d %d %d %d) iy1 = %d ix = %d\n", 
		       itest[iy1], x0, x1, x2, x3, iy1, ix);
return(-1);
	      }
	    }
	    if(x1 == LX+1 && x0 < T) {
	      iy1 = g_idn[ix][1];
	      if(iy1 < gI_0_m2_L_0 || iy1 >= gI_0_m2_L_0 + T*LZ) {
		printf("The DBW2 boundary is not correctly mapped in down x-direction up y %d %d %d %d %d %d\n", 
		       x0, x1, x2, x3, ix, iy1);
return(-1);
	      }
	      itest[iy1]++;
	      if (itest[iy1]>1) {
		printf("The DBW2 boundary is not correctly used x down y up itest = %d (%d %d %d %d) iy1 = %d ix = %d\n", 
		       itest[iy1], x0, x1, x2, x3, iy1, ix);
return(-1);
	      }
	    }
	    
	    x2 = LY+1;
	    ix = g_ipt[x0][x1][x2][x3];
	    iy2=g_idn[ix][2];
	    if((iy2 < gI_0_0_m2_0 || iy2 >= gI_0_0_m2_0 + T*LX*LZ )
	       && x0 < T && x1 < LX) {
	      printf("The DBW2 boundary is not correctly mapped in down y-direction %d %d %d %d %d %d\n", 
		     x0, x1, x2, x3, ix, iy2);
return(-1);
	    }
	    itest[iy2]++;
	    if (itest[iy2]>1) {
	      printf("The DBW2 boundary is not correctly used down y itest = %d (%d %d %d %d) iy2 = %d ix = %d \n", 
		     itest[iy2], x0, x1, x2, x3, iy2, ix);
return(-1);
	    }
	    
	    if(x0 == T && x1 < LX) {
	      iy0 = g_iup[ix][0];
	      if(iy0 < gI_Lp1_0_m1_0 || iy0 >= gI_Lp1_0_m1_0 + LX*LZ) {
		printf("The DBW2 boundary is not correctly mapped in up t-direction down y %d %d %d %d %d %d\n", 
		       x0, x1, x2, x3, ix, iy0);
return(-1);
	      }
	      itest[iy0]++;
	      if (itest[iy0]>1) {
		printf("The DBW2 boundary is not correctly used up t down y itest = %d (%d %d %d %d) iy0 = %d ix = %d \n", 
		       itest[iy0], x0, x1, x2, x3, iy0, ix);
return(-1);
	      }
	    }
	    if(x0 == T+1 && x1 < LX) {
	      iy0 = g_idn[ix][0];
	      if(iy0 < gI_m2_0_m1_0 || iy0 >= gI_m2_0_m1_0 + LX*LZ ) {
		printf("The DBW2 boundary is not correctly mapped in down t-direction down y %d %d %d %d %d %d\n", 
		       x0, x1, x2, x3, ix, iy0);
return(-1);
	      }
	      itest[iy0]++;
	      if (itest[iy0]>1) {
		printf("The DBW2 boundary is not correctly used down t down y itest = %d (%d %d %d %d) iy0 = %d ix = %d \n", 
		       itest[iy0], x0, x1, x2, x3, iy0, ix);
return(-1);
	      }
	    }
	    if(x1 == LX && x0 < T) {
	      iy1 = g_iup[ix][1];
	      if(iy1 < gI_0_Lp1_m1_0 || iy1 >= gI_0_Lp1_m1_0 + T*LZ) {
		printf("The DBW2 boundary is not correctly mapped in up x-direction down y %d %d %d %d %d %d\n", 
		       x0, x1, x2, x3, ix, iy1);
return(-1);
	      }
	      itest[iy1]++;
	      if (itest[iy1]>1) {
		printf("The DBW2 boundary is not correctly used up x down y itest = %d (%d %d %d %d) iy1 = %d ix = %d\n", 
		       itest[iy1], x0, x1, x2, x3, iy1, ix);
return(-1);
	      }
	    }
	    if(x1 == LX+1 && x0 < T) {
	      iy1 = g_idn[ix][1];
	      if(iy1 < gI_0_m2_m1_0 || iy1 >= gI_0_m2_m1_0 + T*LZ) {
		printf("The DBW2 boundary is not correctly mapped in down x-direction down y %d %d %d %d %d %d\n", 
		       x0, x1, x2, x3, ix, iy1);
return(-1);
	      }
	      itest[iy1]++;
	      if (itest[iy1]>1) {
		printf("The DBW2 boundary is not correctly used down x down y itest = %d (%d %d %d %d) iy1 = %d ix = %d\n", 
		       itest[iy1], x0, x1, x2, x3, iy1, ix);
return(-1);
	      }
	    }
	  }
	}
      }
    }
#endif
#if ( defined PARALLELXYZT || defined PARALLELXYZ )
    for (x0 = 0; x0 < T+ext_t; x0++) {
      for (x1 = 0; x1 < LX+ext_x; x1++) {
	for (x2 = 0; x2 < LY+ext_y; x2++) {
	  bndcnt = 0;
	  if(x0 >= T) bndcnt++;
	  if(x1 >= LX) bndcnt++;
	  if(x2 >= LY) bndcnt++;
	  if(bndcnt < 2) {
	    x3 = LZ;
	    ix = g_ipt[x0][x1][x2][x3];
	    
	    iy3=g_iup[ix][3];
	    if(((iy3 < gI_0_0_0_Lp1 || iy3 >= gI_0_0_0_Lp1  + T*LX*LY) && bndcnt == 0) ||
	       (x0 == T && (iy3 < gI_L_0_0_Lp1 || iy3 >= gI_L_0_0_Lp1 + LX*LY )) 
	      ){ 
	      printf("The DBW2 boundary is not correctly mapped in up z-direction %d %d %d %d %d %d\n", 
		     x0, x1, x2, x3, ix, iy3);
return(-1);
	    }
	    itest[iy3]++;
	    if (itest[iy3]>1) {
	      printf("The DBW2 boundary is not correctly used up z itest = %d (%d %d %d %d) iy3 = %d ix = %d \n", 
		     itest[iy3], x0, x1, x2, x3, iy3, ix);
return(-1);
	    }
	    if(x0 == T) {
	      iy0 = g_iup[ix][0];
	      if(iy0 < gI_Lp1_0_0_L || iy0 >= gI_Lp1_0_0_L + LX*LY) {
		printf("The DBW2 boundary is not correctly mapped in up t-direction up z %d %d %d %d %d %d\n", 
		       x0, x1, x2, x3, ix, iy0);
return(-1);
	      }
	      itest[iy0]++;
	      if (itest[iy0]>1) {
		printf("The DBW2 boundary is not correctly used up t up z itest = %d (%d %d %d %d) iy0 = %d ix = %d\n", 
		       itest[iy0], x0, x1, x2, x3, iy0, ix);
return(-1);
	      }
	    }
	    if(x0 == T+1) {
	      iy0 = g_idn[ix][0];
	      if(iy0 < gI_m2_0_0_L || iy0 >= gI_m2_0_0_L + LX*LY) {
		printf("The DBW2 boundary is not correctly mapped in down t-direction up z %d %d %d %d %d %d\n", 
		       x0, x1, x2, x3, ix, iy0);
return(-1);
	      }
	      itest[iy0]++;
	      if (itest[iy0]>1) {
		printf("The DBW2 boundary is not correctly used down t up z itest = %d (%d %d %d %d) iy0 = %d ix = %d\n", 
		       itest[iy0], x0, x1, x2, x3, iy0, ix);
return(-1);
	      }
	    }
	    
	    if(x1 == LX) {
	      iy1 = g_iup[ix][1];
	      if(iy1 < gI_0_Lp1_0_L || iy1 >= gI_0_Lp1_0_L + T*LY) {
		printf("The DBW2 boundary is not correctly mapped in up x-direction up z %d %d %d %d %d %d\n", 
		       x0, x1, x2, x3, ix, iy1);
return(-1);
	      }
	      itest[iy1]++;
	      if (itest[iy1]>1) {
		printf("The DBW2 boundary is not correctly used x up z up itest = %d (%d %d %d %d) iy1 = %d ix = %d\n", 
		       itest[iy1], x0, x1, x2, x3, iy1, ix);
return(-1);
	      }
	    }
	    if(x1 == LX+1) {
	      iy1 = g_idn[ix][1];
	      if(iy1 < gI_0_m2_0_L || iy1 >= gI_0_m2_0_L + T*LY ) {
		printf("The DBW2 boundary is not correctly mapped in down x-direction up z %d %d %d %d %d %d\n", 
		       x0, x1, x2, x3, ix, iy1);
return(-1);
	      }
	      itest[iy1]++;
	      if (itest[iy1]>1) {
		printf("The DBW2 boundary is not correctly used x down z up itest = %d (%d %d %d %d) iy1 = %d ix = %d\n", 
		       itest[iy1], x0, x1, x2, x3, iy1, ix);
return(-1);
	      }
	    }

	    if(x2 == LY) {
	      iy2 = g_iup[ix][2];
	      if(iy2 < gI_0_0_Lp1_L || iy2 >= gI_0_0_Lp1_L + T*LX ) {
		printf("The DBW2 boundary is not correctly mapped in up y-direction up z %d %d %d %d %d %d\n", 
		       x0, x1, x2, x3, ix, iy2);
return(-1);
	      }
	      itest[iy2]++;
	      if (itest[iy2]>1) {
		printf("The DBW2 boundary is not correctly used y up z up itest = %d (%d %d %d %d) iy2 = %d ix = %d\n", 
		       itest[iy2], x0, x1, x2, x3, iy2, ix);
return(-1);
	      }
	    }
	    if(x2 == LY+1) {
	      iy2 = g_idn[ix][2];
	      if(iy2 < gI_0_0_m2_L || iy2 >= gI_0_0_m2_L + T*LX ) {
		printf("The DBW2 boundary is not correctly mapped in down y-direction up z %d %d %d %d %d %d\n", 
		       x0, x1, x2, x3, ix, iy2);
return(-1);
	      }
	      itest[iy2]++;
	      if (itest[iy2]>1) {
		printf("The DBW2 boundary is not correctly used y down z up itest = %d (%d %d %d %d) iy2 = %d ix = %d\n", 
		       itest[iy2], x0, x1, x2, x3, iy2, ix);
return(-1);
	      }
	    }
	    
	    
	    x3 = LZ+1;
	    ix = g_ipt[x0][x1][x2][x3];
	    iy3=g_idn[ix][3];
	    if(((iy3 < gI_0_0_0_m2 || iy3 >= gI_0_0_0_m2 + T*LX*LY) && bndcnt == 0) ||
	       (x0 == T && (iy3 < gI_L_0_0_m2 || iy3 >= gI_L_0_0_m2 + LX*LY)) ||
	       (x0 == T+1 && (iy3 < gI_m1_0_0_m2 || iy3 >= gI_m1_0_0_m2 + LX*LY)) 
		) {
	      printf("The DBW2 boundary is not correctly mapped in down z-direction %d %d %d %d %d %d\n", 
		     x0, x1, x2, x3, ix, iy3);
return(-1);
	    }
	    itest[iy3]++;
	    if (itest[iy3]>1) {
	      printf("The DBW2 boundary is not correctly used down z itest = %d (%d %d %d %d) iy3 = %d ix = %d \n", 
		     itest[iy3], x0, x1, x2, x3, iy3, ix);
return(-1);
	    }

	    if(x0 == T) {
	      iy0 = g_iup[ix][0];
	      if(iy0 < gI_Lp1_0_0_m1 || iy0 >= gI_Lp1_0_0_m1 + LX*LY) {
		printf("The DBW2 boundary is not correctly mapped in up t-direction down z %d %d %d %d %d %d\n", 
		       x0, x1, x2, x3, ix, iy0);
return(-1);
	      }
	      itest[iy0]++;
	      if (itest[iy0]>1) {
		printf("The DBW2 boundary is not correctly used up t down z itest = %d (%d %d %d %d) iy0 = %d ix = %d \n", 
		       itest[iy0], x0, x1, x2, x3, iy0, ix);
return(-1);
	      }
	    }
	    if(x0 == T+1) {
	      iy0 = g_idn[ix][0];
	      if(iy0 < gI_m2_0_0_m1 || iy0 >= gI_m2_0_0_m1 + LX*LY) {
		printf("The DBW2 boundary is not correctly mapped in down t-direction down z %d %d %d %d %d %d\n", 
		       x0, x1, x2, x3, ix, iy0);
return(-1);
	      }
	      itest[iy0]++;
	      if (itest[iy0]>1) {
		printf("The DBW2 boundary is not correctly used down t down z itest = %d (%d %d %d %d) iy0 = %d ix = %d \n", 
		       itest[iy0], x0, x1, x2, x3, iy0, ix);
return(-1);
	      }
	    }
	    if(x1 == LX) {
	      iy1 = g_iup[ix][1];
	      if(iy1 < gI_0_Lp1_0_m1 || iy1 >= gI_0_Lp1_0_m1 + T*LY ) {
		printf("The DBW2 boundary is not correctly mapped in up x-direction down z %d %d %d %d %d %d\n", 
		       x0, x1, x2, x3, ix, iy1);
return(-1);
	      }
	      itest[iy1]++;
	      if (itest[iy1]>1) {
		printf("The DBW2 boundary is not correctly used up x down z itest = %d (%d %d %d %d) iy1 = %d ix = %d\n", 
		       itest[iy1], x0, x1, x2, x3, iy1, ix);
return(-1);
	      }
	    }
	    if(x1 == LX+1) {
	      iy1 = g_idn[ix][1];
	      if(iy1 < gI_0_m2_0_m1 || iy1 >= gI_0_m2_0_m1 + T*LY ) {
		printf("The DBW2 boundary is not correctly mapped in down x-direction down z %d %d %d %d %d %d\n", 
		       x0, x1, x2, x3, ix, iy1);
return(-1);
	      }
	      itest[iy1]++;
	      if (itest[iy1]>1) {
		printf("The DBW2 boundary is not correctly used down x down z itest = %d (%d %d %d %d) iy1 = %d ix = %d\n", 
		       itest[iy1], x0, x1, x2, x3, iy1, ix);
return(-1);
	      }
	    }

	    if(x2 == LY) {
	      iy2 = g_iup[ix][2];
	      if(iy2 < gI_0_0_Lp1_m1 || iy2 >= gI_0_0_Lp1_m1 + T*LX ) {
		printf("The DBW2 boundary is not correctly mapped in up y-direction down z %d %d %d %d %d %d\n", 
		       x0, x1, x2, x3, ix, iy2);
return(-1);
	      }
	      itest[iy2]++;
	      if (itest[iy2]>1) {
		printf("The DBW2 boundary is not correctly used y up z down itest = %d (%d %d %d %d) iy2 = %d ix = %d\n", 
		       itest[iy2], x0, x1, x2, x3, iy2, ix);
return(-1);
	      }
	    }
	    if(x2 == LY+1) {
	      iy2 = g_idn[ix][2];
	      if(iy2 < gI_0_0_m2_m1 || iy2 >= gI_0_0_m2_m1 + T*LX ) {
		printf("The DBW2 boundary is not correctly mapped in down y-direction down z %d %d %d %d %d %d\n", 
		       x0, x1, x2, x3, ix, iy2);
return(-1);
	      }
	      itest[iy2]++;
	      if (itest[iy2]>1) {
		printf("The DBW2 boundary is not correctly used y down z down itest = %d (%d %d %d %d) iy2 = %d ix = %d\n", 
		       itest[iy2], x0, x1, x2, x3, iy2, ix);
return(-1);
	      }
	    }
	  }
	}
      }
    }
#endif
  } /* end of if dbw2>0  */
  for (ix = VOLUMEPLUSRAND; ix < (VOLUMEPLUSRAND) + g_dbw2rand; ix++){ 
    if (itest[ix]!=1) {
      printf("The DBW2 boundary is not correctly used itest = %d ix = %d \n", itest[ix], ix);
return(-1);
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
return(-1);
    }
    stest[iz0] += 1;
    
    iy0 = g_iup[ix][0];
    iz0 = g_lexic2eosub[iy0];
    if(iz0 > VOLUMEPLUSRAND/2 || iz0 < 0) {
      printf("There is a problem with EO geometry in direction 0+\n");
return(-1);
    }
    stest[iz0] += 1;
    
    iy1 = g_idn[ix][1];
    iz1 = g_lexic2eosub[iy1];
    if(iz1 > VOLUMEPLUSRAND/2  || iz1 < 0) {
      printf("There is a problem with EO geometry in direction 1-\n");
return(-1);
    }
    stest[iz1] += 1;
    
    iy1 = g_iup[ix][1];
    iz1 = g_lexic2eosub[iy1];
    if(iz1 >= VOLUMEPLUSRAND/2  || iz1 < 0) {
      printf("There is a problem with EO geometry in direction 1+\n");
return(-1);
    }
    stest[iz1] += 1;
    
    iy2 = g_idn[ix][2];
    iz2 = g_lexic2eosub[iy2];
    if(iz2 > VOLUMEPLUSRAND/2 || iz2 < 0) {
      printf("There is a problem with EO geometry in direction 2-\n");
return(-1);
    }
    stest[iz2] += 1;
    
    iy2 = g_iup[ix][2];
    iz2 = g_lexic2eosub[iy2];
    if(iz2 > VOLUMEPLUSRAND/2 || iz2 < 0) {
      printf("There is a problem with EO geometry in direction 2+\n");
return(-1);
    }
    stest[iz2] += 1;
    
    
    iy3 = g_idn[ix][3];
    iz3 = g_lexic2eosub[iy3];
    if(iz3 > VOLUMEPLUSRAND/2 || iz3 < 0) {
      printf("There is a problem with EO geometry in direction 3-\n");
return(-1);
    }
    stest[iz3] += 1;
    
    iy3 = g_iup[ix][3];
    iz3 = g_lexic2eosub[iy3];
    if(iz3 > VOLUMEPLUSRAND/2 || iz3 < 0) {
      printf("There is a problem with EO geometry in direction 3+\n");
return(-1);
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
return(-1);
  }

  for(j = VOLUME/2; j < (VOLUME+RAND)/2; j++) {
    if(stest[j] != 1) {
      printf("There is a problem in the first boundary of the even odd geometry\n");
return(-1);
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
return(-1);
    }
    stest[iz0] += 1;
    
    iy0 = g_iup[ix][0];
    iz0 = g_lexic2eosub[iy0];
    if(iz0 > VOLUMEPLUSRAND/2 || iz0 < 0) {
      printf("There is a problem with EO geometry in direction 0+\n");
return(-1);
    }
    stest[iz0] += 1;
    
    iy1 = g_idn[ix][1];
    iz1 = g_lexic2eosub[iy1];
    if(iz1 > VOLUMEPLUSRAND/2  || iz1 < 0) {
      printf("There is a problem with EO geometry in direction 1-\n");
return(-1);
    }
    stest[iz1] += 1;
    
    iy1 = g_iup[ix][1];
    iz1 = g_lexic2eosub[iy1];
    if(iz1 > VOLUMEPLUSRAND/2  || iz1 < 0) {
      printf("There is a problem with EO geometry in direction 1+\n"); 
return(-1);
    }
    stest[iz1] += 1; 
    
    iy2 = g_idn[ix][2];
    iz2 = g_lexic2eosub[iy2];
    if(iz2 > VOLUMEPLUSRAND/2 || iz2 < 0) {
      printf("There is a problem with EO geometry in direction 2-\n");
return(-1);
    }
    stest[iz2] += 1;
    
    iy2 = g_iup[ix][2];
    iz2 = g_lexic2eosub[iy2];
    if(iz2 > VOLUMEPLUSRAND/2 || iz2 < 0) {
      printf("There is a problem with EO geometry in direction 2+\n");
return(-1);
    }
    stest[iz2] += 1;
    
    
    iy3 = g_idn[ix][3];
    iz3 = g_lexic2eosub[iy3];
    if(iz3 > VOLUMEPLUSRAND/2 || iz3 < 0) {
      printf("There is a problem with EO geometry in direction 3-\n");
return(-1);
    }
    stest[iz3] += 1;
    
    iy3 = g_iup[ix][3];
    iz3 = g_lexic2eosub[iy3];
    if(iz3 > VOLUMEPLUSRAND/2 || iz3 < 0) {
      printf("There is a problem with EO geometry in direction 3+\n");
return(-1);
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
return(-1);
  }

  for(j = VOLUME/2; j < (VOLUME+RAND)/2; j++) {
    if(stest[j] != 1) {
      printf("There is a problem in the second boundary of the even odd geometry\n");
return(-1);
    }
  }

  if(g_proc_id == 0 ) {
    printf("# The lattice is correctly mapped by the index arrays\n\n");
  }
  fflush(stdout);

  free(stest);
  free(itest);

  return(0);
}

#else /* _INDEX_INDEP_GEOM */

int check_geometry()
{
#ifdef XLC
#pragma execution_frequency(very_low)
#endif
  int ix, j;
  int * stest;
  int * itest;
  int x0,x1,x2,x3;
  int iy0,iy1,iy2,iy3;
  int iz0,iz1,iz2,iz3;   
  int bndcnt = 0;
  
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
return(-1);
	  }
               
	  itest[ix]+=1;
	}
      }
    }
  }

  for (ix = 0; ix < VOLUME; ix++){
    if (itest[ix]!=1){
      printf("The index ipt is not one-to-one %d\n", itest[ix]);
return(-1);
    }
  }

  for (x0 = 0; x0 < T; x0++){
    for (x1 = 0; x1 < LX; x1++){
      for (x2 = 0; x2 < LY; x2++){
	for (x3 = 0; x3 < LZ; x3++){
	  ix=g_ipt[x0][x1][x2][x3];

	  iy0=g_iup[ix][0];
#if (defined PARALLELT || defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT)
	  if(x0!=T-1) {
	    iz0=g_ipt[(x0+1)%T][x1][x2][x3];
	  }
	  else {
	    iz0 = g_ipt[T][x1][x2][x3];
	    itest[iy0]++;
	    if(iy0 < VOLUME  || iy0 >= VOLUME + LX*LY*LZ) {
	      printf("Boundary for time direction up is wrong %d %d %d\n", 
		     iy0 < VOLUME, iy0 >= VOLUME+LX*LY*LZ, iy0);
return(-1);
	    }
	  }
#else
	  iz0=g_ipt[(x0+1)%T][x1][x2][x3];
#endif     
                  
	  iy1=g_iup[ix][1];
#if (defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT)
	  if(x1 !=LX-1) {
	    iz1=g_ipt[x0][(x1+1)%LX][x2][x3];
	  }
	  else {
	    iz1=g_ipt[x0][LX][x2][x3];
	    itest[iy1]++;
	    if(iy1 < VOLUME + 2*LX*LY*LZ || iy1 >= VOLUME + 2*LX*LY*LZ + T*LY*LZ) {
	      printf("Boundary for x direction up is wrong %d %d %d\n", 
		     iy1 < VOLUME + 2*LX*LY*LZ, iy1 >= VOLUME + 2*LX*LY*LZ + T*LY*LZ, iy1);
return(-1);
	    }
	  }
#else
	  iz1=g_ipt[x0][(x1+1)%LX][x2][x3];
#endif

	  iy2=g_iup[ix][2];
#if (defined PARALLELXYT || defined PARALLELXYZT)
	  if(x2 !=LY-1) {
	    iz2=g_ipt[x0][x1][(x2+1)%LY][x3];
	  }
	  else {
	    iz2=g_ipt[x0][x1][LY][x3];
	    itest[iy2]++;
	    if(iy2 < VOLUME + 2*LX*LY*LZ + 2*T*LY*LZ|| iy2 >= VOLUME + 2*LX*LY*LZ + 2*T*LY*LZ + T*LX*LZ) {
	      printf("Boundary for y direction up is wrong %d %d %d\n", 
		     iy2 < VOLUME  + 2*LX*LY*LZ + 2*T*LY*LZ, iy2 > VOLUME + 2*LX*LY*LZ + 2*T*LY*LZ + T*LX*LZ, iy2);
return(-1);
	    }
	  }
#else
	  iz2=g_ipt[x0][x1][(x2+1)%LY][x3];
#endif

	  iy3=g_iup[ix][3];
#if defined PARALLELXYZT
	  if(x3 !=LZ-1) {
	    iz3=g_ipt[x0][x1][x2][(x3+1)%LZ];
	  }
	  else {
	    iz3=g_ipt[x0][x1][x2][LZ];
	    itest[iy3]++;
	    if(iy3 < VOLUME + 2*LX*LY*LZ + 2*T*LY*LZ + 2*T*LX*LZ || iy3 >= VOLUME + 2*LX*LY*LZ + 2*T*LY*LZ + 2*T*LX*LZ + T*LX*LY) {
	      printf("Boundary for z direction up is wrong %d %d %d\n", 
		     iy3 < VOLUME  + 2*LX*LY*LZ + 2*T*LY*LZ + 2*T*LX*LZ, iy3 > VOLUME + 2*LX*LY*LZ + 2*T*LY*LZ + 2*T*LX*LZ + T*LX*LY, iy3);
return(-1);
	    }
	  }
#else
	  iz3=g_ipt[x0][x1][x2][(x3+1)%LZ];
#endif
               
	  if ((iy0!=iz0)||(iy1!=iz1)||(iy2!=iz2)||(iy3!=iz3)||
	      (g_idn[iy0][0]!=ix)||(g_idn[iy1][1]!=ix)||(g_idn[iy2][2]!=ix)||(g_idn[iy3][3]!=ix)) {
	    printf("The index iup is incorrect\n");
	    printf("%d %d %d %d\n", (iy0!=iz0), (iy1!=iz1), (iy2!=iz2), (iy3!=iz3));
	    printf("%d %d %d %d %d %d\n", x0, x1, x2, x3, iy1, iz1);
return(-1);
	  }

	  iy0=g_idn[ix][0];
#if (defined PARALLELT || defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT)
	  if(x0 !=0) {
	    iz0=g_ipt[(x0+T-1)%T][x1][x2][x3];
	  }
	  else {
	    iz0 = g_ipt[T+1][x1][x2][x3];;
	    itest[iy0]++;
	    if(iy0 < VOLUME + LX*LY*LZ  || iy0 >= VOLUME + 2*LX*LY*LZ) {
	      printf("Boundary for time direction is wrong %d %d %d\n", 
		     iy0 < VOLUME + LX*LY*LZ, iy0 >= VOLUME + 2*LX*LY*LZ, iy0);
return(-1);
	    }
	  }
#else
	  iz0=g_ipt[(x0+T-1)%T][x1][x2][x3];
#endif

	  iy1=g_idn[ix][1];
#if (defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT)
	  if(x1 !=0) {
	    iz1=g_ipt[x0][(x1+LX-1)%LX][x2][x3];
	  }
	  else {
	    iz1 = g_ipt[x0][LX+1][x2][x3];
	    itest[iy1]++;
	    if(iy1 < VOLUME + 2*LX*LY*LZ + T*LY*LZ || iy1 >= VOLUME + 2*LX*LY*LZ + 2*T*LY*LZ) {
	      printf("Boundary for x direction is wrong %d %d %d\n", 
		     iy1 < VOLUME + 2*LX*LY*LZ + T*LY*LZ, iy1 >= VOLUME + 2*LX*LY*LZ + 2*T*LY*LZ, iy1);
return(-1);
	    }
	  }
#else
	  iz1=g_ipt[x0][(x1+LX-1)%LX][x2][x3];
#endif
	  iy2=g_idn[ix][2];
#if (defined PARALLELXYT || defined PARALLELXYZT)
	  if(x2 !=0) {
	    iz2=g_ipt[x0][x1][(x2+LY-1)%LY][x3];
	  }
	  else {
	    iz2 = g_ipt[x0][x1][LY+1][x3];
	    itest[iy2]++;
	    if(iy2 < VOLUME + 2*LX*LY*LZ + 2*T*LY*LZ + T*LX*LZ|| iy2 >= VOLUME  + 2*LX*LY*LZ + 2*T*LY*LZ + 2*T*LX*LZ) {
	      printf("Boundary for y direction is wrong %d %d %d\n", 
		     iy2 < VOLUME + 2*LX*LY*LZ + 2*T*LY*LZ + T*LX*LZ, iy2 >= VOLUME  + 2*LX*LY*LZ + 2*T*LY*LZ + 2*T*LX*LZ, iy2);
return(-1);
	    }
	  }
#else
	  iz2=g_ipt[x0][x1][(x2+LY-1)%LY][x3];
#endif

	  iy3=g_idn[ix][3];
#if defined PARALLELXYZT
	  if(x3 !=0) {
	    iz3=g_ipt[x0][x1][x2][(x3+LZ-1)%LZ];
	  }
	  else {
	    iz3 = g_ipt[x0][x1][x2][LZ+1];
	    itest[iy3]++;
	    if(iy3 < VOLUME + 2*LX*LY*LZ + 2*T*LY*LZ + 2*T*LX*LZ + T*LX*LY|| iy3 >= VOLUME + RAND) {
	      printf("Boundary for z direction is wrong %d %d %d\n", 
		     iy3 < VOLUME + 2*LX*LY*LZ + 2*T*LY*LZ + 2*T*LX*LZ + T*LX*LY, iy3 >= VOLUME + RAND, iy3);
	      printf("%d %d %d %d %d\n", x0, x1, x2, x3, ix);
return(-1);
	    }
	  }
#else
	  iz3=g_ipt[x0][x1][x2][(x3+LZ-1)%LZ];
#endif
               
	  if ((iy0!=iz0)||(iy1!=iz1)||(iy2!=iz2)||(iy3!=iz3)||
	      (g_iup[iy0][0]!=ix)||(g_iup[iy1][1]!=ix)||(g_iup[iy2][2]!=ix)||(g_iup[iy3][3]!=ix)) {
	    printf("The index idn is incorrect\n");
	    printf("%d %d %d %d\n", (iy0!=iz0), (iy1!=iz1), (iy2!=iz2), (iy3!=iz3));
	    printf("%d %d %d %d\n", iy1, iz1, iy2, iz2);
	    printf("%d %d %d %d %d\n", x0, x1, x2, x3, ix);
return(-1);
	  }

  /* The edges */
  /* In case of PARALLELT there is actually no edge to take care of */
#if ((defined PARALLELXT) || (defined PARALLELXYT) || (defined PARALLELXYZT))
	  if(x0 == 0) {
	    iy0 = g_idn[ g_idn[ix][1] ][0];
	    if(x1 != 0) {
	      iz0 = g_ipt[T+1][(x1+LX-1)%LX][x2][x3];
	    }
	    else {
	      iz0 = g_ipt[T+1][LX+1][x2][x3];
	      itest[iy0]++;
	    }
	    if(iz0 != iy0) {
	      printf("Edge -t -x has an error\n");
return(-1);
	    }

	    iy0 = g_idn[ g_iup[ix][1] ][0]; 
	    if(x1 != LX-1) {
	      iz0 = g_ipt[T+1][(x1+1)%LX][x2][x3];
	    }
	    else {
	      iz0 = g_ipt[T+1][LX][x2][x3];
	      itest[iy0]++;
	    }
	    if(iz0 != iy0) {
	      printf("Edge -t +x has an error\n");
return(-1);
	    }
	  }

	  if(x0 == T-1) {
	    iy0 = g_iup[ g_idn[ix][1] ][0];
	    if(x1 != 0) {
	      iz0 = g_ipt[T][(x1+LX-1)%LX][x2][x3];
	    }
	    else {
	      iz0 = g_ipt[T][LX+1][x2][x3];
	      itest[iy0]++;
	    }
	    if(iz0 != iy0) {
	      printf("Edge +t -x has an error\n");
return(-1);
	    }

	    iy0 = g_iup[ g_iup[ix][1] ][0]; 
	    if(x1 != LX-1) {
	      iz0 = g_ipt[T][(x1+1)%LX][x2][x3];
	    }
	    else {
	      iz0 = g_ipt[T][LX][x2][x3];
	      itest[iy0]++;
	    }
	    if(iz0 != iy0) {
	      printf("Edge +t +x has an error\n");
return(-1);
	    }
	  }

#endif

#if (defined PARALLELXYT || defined PARALLELXYZT)
	  if(x0 == 0) {
	    iy0 = g_idn[ g_idn[ix][2] ][0]; 
	    if(x2 != 0) {
	      iz0 = g_ipt[T+1][x1][(x2+LY-1)%LY][x3];
	    }
	    else {
	      iz0 = g_ipt[T+1][x1][LY+1][x3];
	      itest[iy0]++;
	    }
	    if(iz0 != iy0) {
	      printf("Edge -t -y has an error\n");
return(-1);
	    }

	    iy0 = g_idn[ g_iup[ix][2] ][0]; 
	    if(x2 != LY-1) {
	      iz0 = g_ipt[T+1][x1][(x2+1)%LY][x3];
	    }
	    else {
	      iz0 = g_ipt[T+1][x1][LY][x3];
	      itest[iy0]++;
	    }
	    if(iz0 != iy0) {
	      printf("Edge -t +y has an error\n");
return(-1);
	    }
	  }

	  if(x0 == T-1) {
	    iy0 = g_iup[ g_idn[ix][2] ][0];
	    if(x2 != 0) {
	      iz0 = g_ipt[T][x1][(x2+LY-1)%LY][x3];
	    }
	    else {
	      iz0 = g_ipt[T][x1][LY+1][x3];
	      itest[iy0]++;
	    }
	    if(iz0 != iy0) {
	      printf("Edge +t -y has an error\n");
return(-1);
	    }

	    iy0 = g_iup[ g_iup[ix][2] ][0]; 
	    if(x2 != LY-1) {
	      iz0 = g_ipt[T][x1][(x2+1)%LY][x3];
	    }
	    else {
	      iz0 = g_ipt[T][x1][LY][x3];
	      itest[iy0]++;
	    }
	    if(iz0 != iy0) {
	      printf("Edge +t +y has an error\n");
return(-1);
	    }
	  }

	  if(x1 == 0) {
	    iy0 = g_idn[ g_idn[ix][2] ][1]; 
	    if(x2 != 0) {
	      iz0 = g_ipt[x0][LX+1][(x2+LY-1)%LY][x3];
	    }
	    else {
	      iz0 = g_ipt[x0][LX+1][LY+1][x3];
	      itest[iy0]++;
	    }
	    if(iz0 != iy0) {
	      printf("Edge -x -y has an error\n");
return(-1);
	    }
	    iy0 = g_idn[ g_iup[ix][2] ][1]; 
	    if(x2 != LY-1) {
	      iz0 = g_ipt[x0][LX+1][(x2+1)%LY][x3];
	    }
	    else {
	      iz0 = g_ipt[x0][LX+1][LY][x3];
	      itest[iy0]++;
	    }
	    if(iz0 != iy0) {
	      printf("Edge -x +y has an error\n");
return(-1);
	    }	    
	  }
	  if(x1 == LX-1) {
	    iy0 = g_iup[ g_idn[ix][2] ][1];
	    if(x2 != 0) {
	      iz0 = g_ipt[x0][LX][(x2+LY-1)%LY][x3];
	    }
	    else {
	      iz0 = g_ipt[x0][LX][LY+1][x3];
	      itest[iy0]++;
	    }
	    if(iz0 != iy0) {
	      printf("Edge +x -y has an error\n");
return(-1);
	    }

	    iy0 = g_iup[ g_iup[ix][2] ][1]; 
	    if(x2 != LY-1) {
	      iz0 = g_ipt[x0][LX][(x2+1)%LY][x3];
	    }
	    else {
	      iz0 = g_ipt[x0][LX][LY][x3];
	      itest[iy0]++;
	    }
	    if(iz0 != iy0) {
	      printf("Edge +x +y has an error\n");
return(-1);
	    }
	  }
#endif
#if defined PARALLELXYZT
	  if(x0 == 0) {
	    iy0 = g_idn[ g_idn[ix][3] ][0]; 
	    if(x3 != 0) {
	      iz0 = g_ipt[T+1][x1][x2][(x3+LZ-1)%LZ];
	    }
	    else {
	      iz0 = g_ipt[T+1][x1][x2][LZ+1];
	      itest[iy0]++;
	      if(itest[iy0]>1) {
		printf("Edge -t -z has itest = %d at %d, %d, %d, %d, iy0 = %d\n", itest[iy0], x0, x1, x2, x3, iy0);
return(-1);
	      }
	    }
	    if(iz0 != iy0) {
	      printf("Edge -t -z has an error\n");
return(-1);
	    }

	    iy0 = g_idn[ g_iup[ix][3] ][0]; 
	    if(x3 != LZ-1) {
	      iz0 = g_ipt[T+1][x1][x2][(x3+1)%LZ];
	    }
	    else {
	      iz0 = g_ipt[T+1][x1][x2][LZ];
	      itest[iy0]++;
	      if(itest[iy0]>1) {
		printf("Edge -t +z has itest = %d at %d, %d, %d, %d, iy0 = %d\n", itest[iy0], x0, x1, x2, x3, iy0);
		printf("ix = %d, iz0 = %d %d\n", ix, iz0, g_iup[ix][3]); 
return(-1);
	      }
	    }
	    if(iz0 != iy0) {
	      printf("Edge -t +z has an error\n");
return(-1);
	    }
	  }

	  if(x0 == T-1) {
	    iy0 = g_iup[ g_idn[ix][3] ][0];
	    if(x3 != 0) {
	      iz0 = g_ipt[T][x1][x2][(x3+LZ-1)%LZ];
	    }
	    else {
	      iz0 = g_ipt[T][x1][x2][LZ+1];
	      itest[iy0]++;
	      if(itest[iy0]>1) {
		printf("Edge +t -z has itest = %d at %d, %d, %d, %d, iy0 = %d\n", itest[iy0], x0, x1, x2, x3, iy0);
return(-1);
	      }
	    }
	    if(iz0 != iy0) {
	      printf("Edge +t -z has an error\n");
return(-1);
	    }

	    iy0 = g_iup[ g_iup[ix][0] ][3]; 
	    if(x3 != LZ-1) {
	      iz0 = g_ipt[T][x1][x2][(x3+1)%LZ];
	    }
	    else {
	      iz0 = g_ipt[T][x1][x2][LZ];
	      itest[iy0]++;
	      if(itest[iy0]>1) {
		printf("Edge +t +z has itest = %d at %d, %d, %d, %d, iy0 = %d\n", itest[iy0], x0, x1, x2, x3, iy0);
return(-1);
	      }
	    }
	    if(iz0 != iy0) {
	      printf("Edge +t +z has an error\n");
return(-1);
	    }
	  }

	  if(x1 == 0) {
	    iy0 = g_idn[ g_idn[ix][3] ][1]; 
	    if(x3 != 0) {
	      iz0 = g_ipt[x0][LX+1][x2][(x3+LZ-1)%LZ];
	    }
	    else {
	      iz0 = g_ipt[x0][LX+1][x2][LZ+1];
	      itest[iy0]++;
	      if(itest[iy0]>1) {
		printf("Edge -x -z has itest = %d at %d, %d, %d, %d, iy0 = %d\n", itest[iy0], x0, x1, x2, x3, iy0);
return(-1);
	      }
	    }
	    if(iz0 != iy0) {
	      printf("Edge -x -z has an error\n");
return(-1);
	    }
	    iy0 = g_idn[ g_iup[ix][3] ][1]; 
	    if(x3 != LZ-1) {
	      iz0 = g_ipt[x0][LX+1][x2][(x3+1)%LZ];
	    }
	    else {
	      iz0 = g_ipt[x0][LX+1][x2][LZ];
	      itest[iy0]++;
	      if(itest[iy0]>1) {
		printf("Edge -x +z has itest = %d at %d, %d, %d, %d, iy0 = %d\n", itest[iy0], x0, x1, x2, x3, iy0);
return(-1);
	      }
	    }
	    if(iz0 != iy0) {
	      printf("Edge -x +z has an error\n");
return(-1);
	    }	    
	  }
	  if(x1 == LX-1) {
	    iy0 = g_iup[ g_idn[ix][3] ][1];
	    if(x3 != 0) {
	      iz0 = g_ipt[x0][LX][x2][(x3+LZ-1)%LZ];
	    }
	    else {
	      iz0 = g_ipt[x0][LX][x2][LZ+1];
	      itest[iy0]++;
	      if(itest[iy0]>1) {
		printf("Edge +x -z has itest = %d at %d, %d, %d, %d, iy0 = %d\n", itest[iy0], x0, x1, x2, x3, iy0);
return(-1);
	      }
	    }
	    if(iz0 != iy0) {
	      printf("Edge +x -z has an error\n");
return(-1);
	    }

	    iy0 = g_iup[ g_iup[ix][3] ][1]; 
	    if(x3 != LZ-1) {
	      iz0 = g_ipt[x0][LX][x2][(x3+1)%LZ];
	    }
	    else {
	      iz0 = g_ipt[x0][LX][x2][LZ];
	      itest[iy0]++;
	      if(itest[iy0]>1) {
		printf("Edge +x +z has itest = %d at %d, %d, %d, %d, iy0 = %d iz0 = %d giup[%d][3] = %d ix = %d\n",
		       itest[iy0], x0, x1, x2, x3, iy0, iz0, ix, g_iup[ix][3], ix);
return(-1);
	      }
	    }
	    if(iz0 != iy0) {
	      printf("Edge +x +z has an error\n");
return(-1);
	    }
	  }

	  if(x2 == 0) {
	    iy0 = g_idn[ g_idn[ix][3] ][2]; 
	    if(x3 != 0) {
	      iz0 = g_ipt[x0][x1][LY+1][(x3+LZ-1)%LZ];
	    }
	    else {
	      iz0 = g_ipt[x0][x1][LY+1][LZ+1];
	      itest[iy0]++;
	      if(itest[iy0]>1) {
		printf("Edge -y -z has itest = %d at %d, %d, %d, %d, iy0 = %d\n", itest[iy0], x0, x1, x2, x3, iy0);
return(-1);
	      }
	    }
	    if(iz0 != iy0) {
	      printf("Edge -y -z has an error\n");
return(-1);
	    }
	    iy0 = g_idn[ g_iup[ix][3] ][2]; 
	    if(x3 != LZ-1) {
	      iz0 = g_ipt[x0][x1][LY+1][(x3+1)%LZ];
	    }
	    else {
	      iz0 = g_ipt[x0][x1][LY+1][LZ];
	      itest[iy0]++;
	      if(itest[iy0]>1) {
		printf("Edge -y +z has itest = %d at %d, %d, %d, %d, iy0 = %d\n", itest[iy0], x0, x1, x2, x3, iy0);
return(-1);
	      }
	    }
	    if(iz0 != iy0) {
	      printf("Edge -y +z has an error\n");
return(-1);
	    }	    
	  }
	  if(x2 == LY-1) {
	    iy0 = g_iup[ g_idn[ix][3] ][2];
	    if(x3 != 0) {
	      iz0 = g_ipt[x0][x1][LY][(x3+LZ-1)%LZ];
	    }
	    else {
	      iz0 = g_ipt[x0][x1][LY][LZ+1];
	      itest[iy0]++;
	      if(itest[iy0]>1) {
		printf("Edge +y -z has itest = %d at %d, %d, %d, %d, iy0 = %d\n", itest[iy0], x0, x1, x2, x3, iy0);
return(-1);
	      }
	    }
	    if(iz0 != iy0) {
	      printf("Edge +y -z has an error\n");
return(-1);
	    }

	    iy0 = g_iup[ g_iup[ix][3] ][2]; 
	    if(x3 != LZ-1) {
	      iz0 = g_ipt[x0][x1][LY][(x3+1)%LZ];
	    }
	    else {
	      iz0 = g_ipt[x0][x1][LY][LZ];
	      itest[iy0]++;
	      if(itest[iy0]>1) {
		printf("Edge +y +z has itest = %d at %d, %d, %d, %d, iy0 = %d\n", itest[iy0], x0, x1, x2, x3, iy0);
return(-1);
	      }
	    }
	    if(iz0 != iy0) {
	      printf("Edge +y +z has an error\n");
return(-1);
	    }
	  }
#endif
	}
      }
    }
  }
  
  for (ix = VOLUME; ix < (VOLUME+RAND+EDGES); ix++){
    if (itest[ix]!=1) {
      printf("The boundary is not correctly used itest = %d ix = %d %d %d %d\n", itest[ix], ix, VOLUME, RAND, EDGES);
return(-1);
    }
  }

  for (ix=0;ix<VOLUMEPLUSRAND + g_dbw2rand;ix++){
    itest[ix]=0;
  }

#ifdef MPI
  if(g_dbw2rand > 0) {
    for (x1 = 0; x1 < LX; x1++) {
      for (x2 = 0; x2 < LY; x2++) {
	for (x3 = 0; x3 < LZ; x3++) {
	  x0 = T;
	  ix=g_ipt[x0][x1][x2][x3];

	  iy0=g_iup[ix][0];
	  if(iy0 < VOLUMEPLUSRAND || iy0 >= VOLUMEPLUSRAND + LX*LY*LZ) {
	    printf("The DBW2 boundary is not correctly mapped in up t-direction  %d %d %d %d %d %d\n", x0, x1, x2, x3, ix, iy0);
return(-1);
	  }
	  itest[iy0]++;
	  if (itest[iy0]>1) {
	    printf("The DBW2 boundary is not correctly used itest = %d (%d %d %d %d) iy0 = %d ix = %d\n", itest[iy0], x0, x1, x2, x3, iy0, ix);
return(-1);
	  }
	  
	  x0 = T+1;
	  ix=g_ipt[x0][x1][x2][x3];

	  iy0=g_idn[ix][0];
	  if(iy0 < VOLUMEPLUSRAND + LX*LZ*LY|| iy0 >= VOLUMEPLUSRAND + 2*LX*LY*LZ) {
	    printf("The DBW2 boundary is not correctly mapped in down t-direction %d %d %d %d %d %d\n", x0, x1, x2, x3, ix, iy0);
return(-1);
	  }
	  itest[iy0]++;
	  if (itest[iy0]>1) {
	    printf("The DBW2 boundary is not correctly used itest = %d (%d %d %d %d) iy0 = %d ix = %d\n", itest[iy0], x0, x1, x2, x3, iy0, ix);
return(-1);
	  }
	}
      }
    }

#if (defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT)
    for (x0 = 0; x0 < T+2; x0++) {
      for (x2 = 0; x2 < LY; x2++) {
	for (x3 = 0; x3 < LZ; x3++) {
	  x1 = LX;
	  ix = g_ipt[x0][x1][x2][x3];

	  iy1=g_iup[ix][1];
	  if((iy1 < VOLUMEPLUSRAND + 2*LX*LY*LZ || iy1 >= VOLUMEPLUSRAND + 2*LX*LY*LZ + T*LY*LZ) && x0 < T) {
	    printf("The DBW2 boundary is not correctly mapped in up x-direction %d %d %d %d %d %d\n", x0, x1, x2, x3, ix, iy1);
return(-1);
	  }
	  itest[iy1]++;
	  if (itest[iy1]>1) {
	    printf("The DBW2 boundary is not correctly used up x itest = %d (%d %d %d %d) iy1 = %d ix = %d \n", itest[iy1], x0, x1, x2, x3, iy1, ix);
return(-1);
	  }

	  if(x0 == T) {
	    iy0 = g_iup[ix][0];
	    if(iy0 < VOLUMEPLUSRAND + RAND || iy0 >= VOLUMEPLUSRAND + RAND + LY*LZ) {
	      printf("The DBW2 boundary is not correctly mapped in up t-direction up x %d %d %d %d %d %d\n", x0, x1, x2, x3, ix, iy0);
return(-1);
	    }
	    itest[iy0]++;
	    if (itest[iy0]>1) {
	    printf("The DBW2 boundary is not correctly used up t up x itest = %d (%d %d %d %d) iy0 = %d ix = %d\n", itest[iy0], x0, x1, x2, x3, iy0, ix);
return(-1);
	    }
	  }
	  if(x0 == T+1) {
	    iy0 = g_idn[ix][0];
	    if(iy0 < VOLUMEPLUSRAND + RAND + 2*LY*LZ|| iy0 >= VOLUMEPLUSRAND + RAND + 3*LY*LZ) {
	      printf("The DBW2 boundary is not correctly mapped in down t-direction up x %d %d %d %d %d %d\n", x0, x1, x2, x3, ix, iy0);
return(-1);
	    }
	    itest[iy0]++;
	    if (itest[iy0]>1) {
	    printf("The DBW2 boundary is not correctly used down t up x itest = %d (%d %d %d %d) iy0 = %d ix = %d\n", itest[iy0], x0, x1, x2, x3, iy0, ix);
return(-1);
	    }
	  }


	  x1 = LX+1;
	  ix = g_ipt[x0][x1][x2][x3];

	  iy1=g_idn[ix][1];
	  if((iy1 < VOLUMEPLUSRAND +  2*LX*LY*LZ + T*LY*LZ|| iy1 >= VOLUMEPLUSRAND + 2*LX*LY*LZ + 2*T*LY*LZ) && x0 < T) {
	    printf("The DBW2 boundary is not correctly mapped in down x-direction %d %d %d %d %d %d\n", x0, x1, x2, x3, ix, iy1);
return(-1);
	  }
	  itest[iy1]++;
	  if (itest[iy1]>1) {
	    printf("The DBW2 boundary is not correctly used down x itest = %d (%d %d %d %d) iy1 = %d ix = %d \n", itest[iy1], x0, x1, x2, x3, iy1, ix);
return(-1);
	  }

	  if(x0 == T) {
	    iy0 = g_iup[ix][0];
	    if(iy0 < VOLUMEPLUSRAND + RAND + LY*LZ || iy0 >= VOLUMEPLUSRAND + RAND + 2*LY*LZ) {
	      printf("The DBW2 boundary is not correctly mapped in up t-direction down x %d %d %d %d %d %d\n", x0, x1, x2, x3, ix, iy0);
return(-1);
	    }
	    itest[iy0]++;
	    if (itest[iy0]>1) {
	      printf("The DBW2 boundary is not correctly used up t down x itest = %d (%d %d %d %d) iy0 = %d ix = %d \n", itest[iy0], x0, x1, x2, x3, iy0, ix);
return(-1);
	    }
	  }
	  if(x0 == T+1) {
	    iy0 = g_idn[ix][0];
	    if(iy0 < VOLUMEPLUSRAND + RAND + 3*LY*LZ || iy0 >= VOLUMEPLUSRAND + RAND + 4*LY*LZ) {
	      printf("The DBW2 boundary is not correctly mapped in down t-direction down x %d %d %d %d %d %d\n", x0, x1, x2, x3, ix, iy0);
return(-1);
	    }
	    itest[iy0]++;
	    if (itest[iy0]>1) {
	      printf("The DBW2 boundary is not correctly used down t down xitest = %d (%d %d %d %d) iy0 = %d ix = %d \n", itest[iy0], x0, x1, x2, x3, iy0, ix);
return(-1);
	    }
	  }
	}
      }
    }
#endif

#if (defined PARALLELXYT || defined PARALLELXYZT)

    for (x0 = 0; x0 < T+2; x0++) {
      for (x1 = 0; x1 < LX+2; x1++) {
	for (x3 = 0; x3 < LZ; x3++) {
	  if(x0 < T || x1 < LX) {
	    x2 = LY;
	    ix = g_ipt[x0][x1][x2][x3];
	    
	    iy2=g_iup[ix][2];
	    if((iy2 < VOLUMEPLUSRAND + 2*LX*LY*LZ + 2*T*LY*LZ || iy2 >= VOLUMEPLUSRAND + 2*LX*LY*LZ + 2*T*LY*LZ + T*LX*LZ) 
	       && x0 < T && x1 < LX) {
	      printf("The DBW2 boundary is not correctly mapped in up y-direction %d %d %d %d %d %d\n", 
		     x0, x1, x2, x3, ix, iy2);
return(-1);
	    }
	    itest[iy2]++;
	    if (itest[iy2]>1) {
	      printf("The DBW2 boundary is not correctly used up y itest = %d (%d %d %d %d) iy2 = %d ix = %d \n", 
		     itest[iy2], x0, x1, x2, x3, iy2, ix);
return(-1);
	    }
	    
	    if(x0 == T && x1 < LX) {
	      iy0 = g_iup[ix][0];
	      if(iy0 < VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ || iy0 >= VOLUMEPLUSRAND+g_dbw2rand) {
		printf("The DBW2 boundary is not correctly mapped in up t-direction up y %d %d %d %d %d %d\n", 
		       x0, x1, x2, x3, ix, iy0);
return(-1);
	      }
	      itest[iy0]++;
	      if (itest[iy0]>1) {
		printf("The DBW2 boundary is not correctly used up t up y itest = %d (%d %d %d %d) iy0 = %d ix = %d\n", 
		       itest[iy0], x0, x1, x2, x3, iy0, ix);
return(-1);
	      }
	    }
	    if(x0 == T+1 && x1 < LX) {
	      iy0 = g_idn[ix][0];
	      if(iy0 < VOLUMEPLUSRAND || iy0 >= VOLUMEPLUSRAND+g_dbw2rand) {
		printf("The DBW2 boundary is not correctly mapped in down t-direction up y %d %d %d %d %d %d\n", 
		       x0, x1, x2, x3, ix, iy0);
return(-1);
	      }
	      itest[iy0]++;
	      if (itest[iy0]>1) {
		printf("The DBW2 boundary is not correctly used down t up y itest = %d (%d %d %d %d) iy0 = %d ix = %d\n", 
		       itest[iy0], x0, x1, x2, x3, iy0, ix);
return(-1);
	      }
	    }
	    
	    if(x1 == LX && x0 < T) {
	      iy1 = g_iup[ix][1];
	      if(iy1 < VOLUMEPLUSRAND || iy1 >= VOLUMEPLUSRAND + g_dbw2rand) {
		printf("The DBW2 boundary is not correctly mapped in up x-direction up y %d %d %d %d %d %d\n", 
		       x0, x1, x2, x3, ix, iy1);
return(-1);
	      }
	      itest[iy1]++;
	      if (itest[iy1]>1) {
		printf("The DBW2 boundary is not correctly used x up y up itest = %d (%d %d %d %d) iy1 = %d ix = %d\n", 
		       itest[iy1], x0, x1, x2, x3, iy1, ix);
return(-1);
	      }
	    }
	    if(x1 == LX+1 && x0 < T) {
	      iy1 = g_idn[ix][1];
	      if(iy1 < VOLUMEPLUSRAND || iy1 >= VOLUMEPLUSRAND+g_dbw2rand) {
		printf("The DBW2 boundary is not correctly mapped in down x-direction up y %d %d %d %d %d %d\n", 
		       x0, x1, x2, x3, ix, iy1);
return(-1);
	      }
	      itest[iy1]++;
	      if (itest[iy1]>1) {
		printf("The DBW2 boundary is not correctly used x down y up itest = %d (%d %d %d %d) iy1 = %d ix = %d\n", 
		       itest[iy1], x0, x1, x2, x3, iy1, ix);
return(-1);
	      }
	    }
	    
	    
	    x2 = LY+1;
	    ix = g_ipt[x0][x1][x2][x3];
	    iy2=g_idn[ix][2];
	    if((iy2 < VOLUMEPLUSRAND + 2*LX*LY*LZ + 2*T*LY*LZ + T*LX*LZ|| iy2 >= VOLUMEPLUSRAND + 2*LX*LY*LZ + 2*T*LY*LZ + 2*T*LX*LZ )
	       && x0 < T && x1 < LX) {
	      printf("The DBW2 boundary is not correctly mapped in down y-direction %d %d %d %d %d %d\n", 
		     x0, x1, x2, x3, ix, iy2);
return(-1);
	    }
	    itest[iy2]++;
	    if (itest[iy2]>1) {
	      printf("The DBW2 boundary is not correctly used down y itest = %d (%d %d %d %d) iy2 = %d ix = %d \n", 
		     itest[iy2], x0, x1, x2, x3, iy2, ix);
return(-1);
	    }
	    
	    if(x0 == T && x1 < LX) {
	      iy0 = g_iup[ix][0];
	      if(iy0 < VOLUMEPLUSRAND || iy0 >= VOLUMEPLUSRAND+g_dbw2rand) {
		printf("The DBW2 boundary is not correctly mapped in up t-direction down y %d %d %d %d %d %d\n", 
		       x0, x1, x2, x3, ix, iy0);
return(-1);
	      }
	      itest[iy0]++;
	      if (itest[iy0]>1) {
		printf("The DBW2 boundary is not correctly used up t down y itest = %d (%d %d %d %d) iy0 = %d ix = %d \n", 
		       itest[iy0], x0, x1, x2, x3, iy0, ix);
return(-1);
	      }
	    }
	    if(x0 == T+1 && x1 < LX) {
	      iy0 = g_idn[ix][0];
	      if(iy0 < VOLUMEPLUSRAND || iy0 >= VOLUMEPLUSRAND+g_dbw2rand) {
		printf("The DBW2 boundary is not correctly mapped in down t-direction down y %d %d %d %d %d %d\n", 
		       x0, x1, x2, x3, ix, iy0);
return(-1);
	      }
	      itest[iy0]++;
	      if (itest[iy0]>1) {
		printf("The DBW2 boundary is not correctly used down t down y itest = %d (%d %d %d %d) iy0 = %d ix = %d \n", 
		       itest[iy0], x0, x1, x2, x3, iy0, ix);
return(-1);
	      }
	    }
	    if(x1 == LX && x0 < T) {
	      iy1 = g_iup[ix][1];
	      if(iy1 < VOLUMEPLUSRAND || iy1 >= VOLUMEPLUSRAND+g_dbw2rand) {
		printf("The DBW2 boundary is not correctly mapped in up x-direction down y %d %d %d %d %d %d\n", 
		       x0, x1, x2, x3, ix, iy1);
return(-1);
	      }
	      itest[iy1]++;
	      if (itest[iy1]>1) {
		printf("The DBW2 boundary is not correctly used up x down y itest = %d (%d %d %d %d) iy1 = %d ix = %d\n", 
		       itest[iy1], x0, x1, x2, x3, iy1, ix);
return(-1);
	      }
	    }
	    if(x1 == LX+1 && x0 < T) {
	      iy1 = g_idn[ix][1];
	      if(iy1 < VOLUMEPLUSRAND || iy1 >= VOLUMEPLUSRAND+g_dbw2rand) {
		printf("The DBW2 boundary is not correctly mapped in down x-direction down y %d %d %d %d %d %d\n", 
		       x0, x1, x2, x3, ix, iy1);
return(-1);
	      }
	      itest[iy1]++;
	      if (itest[iy1]>1) {
		printf("The DBW2 boundary is not correctly used down x down y itest = %d (%d %d %d %d) iy1 = %d ix = %d\n", 
		       itest[iy1], x0, x1, x2, x3, iy1, ix);
return(-1);
	      }
	    }
	  }
	}
      }
    }
#endif
#ifdef PARALLELXYZT
    for (x0 = 0; x0 < T+2; x0++) {
      for (x1 = 0; x1 < LX+2; x1++) {
	for (x2 = 0; x2 < LY+2; x2++) {
	  bndcnt = 0;
	  if(x0 >= T) bndcnt++;
	  if(x1 >= LX) bndcnt++;
	  if(x2 >= LY) bndcnt++;
	  if(bndcnt < 2) {
	    x3 = LZ;
	    ix = g_ipt[x0][x1][x2][x3];
	    
	    iy3=g_iup[ix][3];
	    if(((iy3 < VOLUMEPLUSRAND + 2*LX*LY*LZ + 2*T*LY*LZ + 2*T*LX*LZ || 
	       iy3 >= VOLUMEPLUSRAND + 2*LX*LY*LZ + 2*T*LY*LZ + 2*T*LX*LZ + T*LX*LY) &&
	       bndcnt == 0) ||
	       (x0 == T && (iy3 < VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 4*LX*LY ||
		iy3 >= VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 5*LX*LY )) 
/* 	       ||(x0 == T+1 && (iy3 < VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 5*LX*LY || */
/* 			      iy3 >= VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 6*LX*LY)) */
	      ){ 
	      printf("The DBW2 boundary is not correctly mapped in up z-direction %d %d %d %d %d %d\n", 
		     x0, x1, x2, x3, ix, iy3);
return(-1);
	    }
	    itest[iy3]++;
	    if (itest[iy3]>1) {
	      printf("The DBW2 boundary is not correctly used up z itest = %d (%d %d %d %d) iy3 = %d ix = %d \n", 
		     itest[iy3], x0, x1, x2, x3, iy3, ix);
return(-1);
	    }
	    
	    if(x0 == T) {
	      iy0 = g_iup[ix][0];
	      if(iy0 < VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ || 
	         iy0 >= VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + LX*LY) {
		printf("The DBW2 boundary is not correctly mapped in up t-direction up z %d %d %d %d %d %d\n", 
		       x0, x1, x2, x3, ix, iy0);
return(-1);
	      }
	      itest[iy0]++;
	      if (itest[iy0]>1) {
		printf("The DBW2 boundary is not correctly used up t up z itest = %d (%d %d %d %d) iy0 = %d ix = %d\n", 
		       itest[iy0], x0, x1, x2, x3, iy0, ix);
return(-1);
	      }
	    }
	    if(x0 == T+1) {
	      iy0 = g_idn[ix][0];
	      if(iy0 < VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + LX*LY || 
                 iy0 >= VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 2*LX*LY) {
		printf("The DBW2 boundary is not correctly mapped in down t-direction up z %d %d %d %d %d %d\n", 
		       x0, x1, x2, x3, ix, iy0);
return(-1);
	      }
	      itest[iy0]++;
	      if (itest[iy0]>1) {
		printf("The DBW2 boundary is not correctly used down t up z itest = %d (%d %d %d %d) iy0 = %d ix = %d\n", 
		       itest[iy0], x0, x1, x2, x3, iy0, ix);
return(-1);
	      }
	    }
	    
	    if(x1 == LX) {
	      iy1 = g_iup[ix][1];
	      if(iy1 < VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 4*T*LY|| 
	         iy1 >= VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 5*T*LY) {
		printf("The DBW2 boundary is not correctly mapped in up x-direction up z %d %d %d %d %d %d\n", 
		       x0, x1, x2, x3, ix, iy1);
return(-1);
	      }
	      itest[iy1]++;
	      if (itest[iy1]>1) {
		printf("The DBW2 boundary is not correctly used x up z up itest = %d (%d %d %d %d) iy1 = %d ix = %d\n", 
		       itest[iy1], x0, x1, x2, x3, iy1, ix);
return(-1);
	      }
	    }
	    if(x1 == LX+1) {
	      iy1 = g_idn[ix][1];
	      if(iy1 < VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 6*T*LY || 
	         iy1 >= VOLUMEPLUSRAND + + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 7*T*LY ) {
		printf("The DBW2 boundary is not correctly mapped in down x-direction up z %d %d %d %d %d %d\n", 
		       x0, x1, x2, x3, ix, iy1);
return(-1);
	      }
	      itest[iy1]++;
	      if (itest[iy1]>1) {
		printf("The DBW2 boundary is not correctly used x down z up itest = %d (%d %d %d %d) iy1 = %d ix = %d\n", 
		       itest[iy1], x0, x1, x2, x3, iy1, ix);
return(-1);
	      }
	    }

	    if(x2 == LY) {
	      iy2 = g_iup[ix][2];
	      if(iy2 < VOLUMEPLUSRAND  + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 8*T*LY + 4*T*LX || 
	         iy2 >= VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 8*T*LY + 5*T*LX ) {
		printf("The DBW2 boundary is not correctly mapped in up y-direction up z %d %d %d %d %d %d\n", 
		       x0, x1, x2, x3, ix, iy2);
return(-1);
	      }
	      itest[iy2]++;
	      if (itest[iy2]>1) {
		printf("The DBW2 boundary is not correctly used y up z up itest = %d (%d %d %d %d) iy2 = %d ix = %d\n", 
		       itest[iy2], x0, x1, x2, x3, iy2, ix);
return(-1);
	      }
	    }
	    if(x2 == LY+1) {
	      iy2 = g_idn[ix][2];
	      if(iy2 < VOLUMEPLUSRAND  + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 8*T*LY + 6*T*LX || 
	         iy2 >= VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 8*T*LY + 7*T*LX ) {
		printf("The DBW2 boundary is not correctly mapped in down y-direction up z %d %d %d %d %d %d\n", 
		       x0, x1, x2, x3, ix, iy2);
return(-1);
	      }
	      itest[iy2]++;
	      if (itest[iy2]>1) {
		printf("The DBW2 boundary is not correctly used y down z up itest = %d (%d %d %d %d) iy2 = %d ix = %d\n", 
		       itest[iy2], x0, x1, x2, x3, iy2, ix);
return(-1);
	      }
	    }
	    
	    
	    x3 = LZ+1;
	    ix = g_ipt[x0][x1][x2][x3];
	    iy3=g_idn[ix][3];
	    if(((iy3 < VOLUMEPLUSRAND + 2*LX*LY*LZ + 2*T*LY*LZ + 2*T*LX*LZ + T*LX*LY|| 
	       iy3 >= VOLUMEPLUSRAND + 2*LX*LY*LZ + 2*T*LY*LZ + 2*T*LX*LZ + 2*T*LX*LY) &&
	       bndcnt == 0) ||
	       (x0 == T && (iy3 < VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 6*LX*LY||
			    iy3 >= VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 7*LX*LY)) ||
	       (x0 == T+1 && (iy3 < VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 7*LX*LY||
			    iy3 >= VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY)) 
		) {
	      printf("The DBW2 boundary is not correctly mapped in down z-direction %d %d %d %d %d %d\n", 
		     x0, x1, x2, x3, ix, iy3);
return(-1);
	    }
	    itest[iy3]++;
	    if (itest[iy3]>1) {
	      printf("The DBW2 boundary is not correctly used down z itest = %d (%d %d %d %d) iy3 = %d ix = %d \n", 
		     itest[iy3], x0, x1, x2, x3, iy3, ix);
return(-1);
	    }
	    
	    if(x0 == T) {
	      iy0 = g_iup[ix][0];
	      if(iy0 < VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 2*LX*LY ||
	         iy0 >= VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 3*LX*LY) {
		printf("The DBW2 boundary is not correctly mapped in up t-direction down z %d %d %d %d %d %d\n", 
		       x0, x1, x2, x3, ix, iy0);
return(-1);
	      }
	      itest[iy0]++;
	      if (itest[iy0]>1) {
		printf("The DBW2 boundary is not correctly used up t down z itest = %d (%d %d %d %d) iy0 = %d ix = %d \n", 
		       itest[iy0], x0, x1, x2, x3, iy0, ix);
return(-1);
	      }
	    }
	    if(x0 == T+1) {
	      iy0 = g_idn[ix][0];
	      if(iy0 < VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 3*LX*LY || 
	         iy0 >= VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 4*LX*LY) {
		printf("The DBW2 boundary is not correctly mapped in down t-direction down z %d %d %d %d %d %d\n", 
		       x0, x1, x2, x3, ix, iy0);
return(-1);
	      }
	      itest[iy0]++;
	      if (itest[iy0]>1) {
		printf("The DBW2 boundary is not correctly used down t down z itest = %d (%d %d %d %d) iy0 = %d ix = %d \n", 
		       itest[iy0], x0, x1, x2, x3, iy0, ix);
return(-1);
	      }
	    }
	    if(x1 == LX) {
	      iy1 = g_iup[ix][1];
	      if(iy1 < VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 5*T*LY || 
	         iy1 >= VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 6*T*LY ) {
		printf("The DBW2 boundary is not correctly mapped in up x-direction down z %d %d %d %d %d %d\n", 
		       x0, x1, x2, x3, ix, iy1);
return(-1);
	      }
	      itest[iy1]++;
	      if (itest[iy1]>1) {
		printf("The DBW2 boundary is not correctly used up x down z itest = %d (%d %d %d %d) iy1 = %d ix = %d\n", 
		       itest[iy1], x0, x1, x2, x3, iy1, ix);
return(-1);
	      }
	    }
	    if(x1 == LX+1) {
	      iy1 = g_idn[ix][1];
	      if(iy1 < VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 7*T*LY || 
	         iy1 >= VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 8*T*LY ) {
		printf("The DBW2 boundary is not correctly mapped in down x-direction down z %d %d %d %d %d %d\n", 
		       x0, x1, x2, x3, ix, iy1);
return(-1);
	      }
	      itest[iy1]++;
	      if (itest[iy1]>1) {
		printf("The DBW2 boundary is not correctly used down x down z itest = %d (%d %d %d %d) iy1 = %d ix = %d\n", 
		       itest[iy1], x0, x1, x2, x3, iy1, ix);
return(-1);
	      }
	    }

	    if(x2 == LY) {
	      iy2 = g_iup[ix][2];
	      if(iy2 < VOLUMEPLUSRAND  + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 8*T*LY + 5*T*LX || 
	         iy2 >= VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 8*T*LY + 6*T*LX ) {
		printf("The DBW2 boundary is not correctly mapped in up y-direction down z %d %d %d %d %d %d\n", 
		       x0, x1, x2, x3, ix, iy2);
return(-1);
	      }
	      itest[iy2]++;
	      if (itest[iy2]>1) {
		printf("The DBW2 boundary is not correctly used y up z down itest = %d (%d %d %d %d) iy2 = %d ix = %d\n", 
		       itest[iy2], x0, x1, x2, x3, iy2, ix);
return(-1);
	      }
	    }
	    if(x2 == LY+1) {
	      iy2 = g_idn[ix][2];
	      if(iy2 < VOLUMEPLUSRAND  + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 8*T*LY + 7*T*LX || 
	         iy2 >= VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 8*T*LY + 8*T*LX ) {
		printf("The DBW2 boundary is not correctly mapped in down y-direction down z %d %d %d %d %d %d\n", 
		       x0, x1, x2, x3, ix, iy2);
return(-1);
	      }
	      itest[iy2]++;
	      if (itest[iy2]>1) {
		printf("The DBW2 boundary is not correctly used y down z down itest = %d (%d %d %d %d) iy2 = %d ix = %d\n", 
		       itest[iy2], x0, x1, x2, x3, iy2, ix);
return(-1);
	      }
	    }
	  }
	}
      }
    }
#endif
  }
  for (ix = VOLUMEPLUSRAND; ix < (VOLUMEPLUSRAND) + g_dbw2rand; ix++){ 
    if (itest[ix]!=1) {
      printf("The DBW2 boundary is not correctly used itest = %d ix = %d \n", itest[ix], ix);
return(-1);
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
return(-1);
    }
    stest[iz0] += 1;
    
    iy0 = g_iup[ix][0];
    iz0 = g_lexic2eosub[iy0];
    if(iz0 > VOLUMEPLUSRAND/2 || iz0 < 0) {
      printf("There is a problem with EO geometry in direction 0+\n");
return(-1);
    }
    stest[iz0] += 1;
    
    iy1 = g_idn[ix][1];
    iz1 = g_lexic2eosub[iy1];
    if(iz1 > VOLUMEPLUSRAND/2  || iz1 < 0) {
      printf("There is a problem with EO geometry in direction 1-\n");
return(-1);
    }
    stest[iz1] += 1;
    
    iy1 = g_iup[ix][1];
    iz1 = g_lexic2eosub[iy1];
    if(iz1 >= VOLUMEPLUSRAND/2  || iz1 < 0) {
      printf("There is a problem with EO geometry in direction 1+\n");
return(-1);
    }
    stest[iz1] += 1;
    
    iy2 = g_idn[ix][2];
    iz2 = g_lexic2eosub[iy2];
    if(iz2 > VOLUMEPLUSRAND/2 || iz2 < 0) {
      printf("There is a problem with EO geometry in direction 2-\n");
return(-1);
    }
    stest[iz2] += 1;
    
    iy2 = g_iup[ix][2];
    iz2 = g_lexic2eosub[iy2];
    if(iz2 > VOLUMEPLUSRAND/2 || iz2 < 0) {
      printf("There is a problem with EO geometry in direction 2+\n");
return(-1);
    }
    stest[iz2] += 1;
    
    
    iy3 = g_idn[ix][3];
    iz3 = g_lexic2eosub[iy3];
    if(iz3 > VOLUMEPLUSRAND/2 || iz3 < 0) {
      printf("There is a problem with EO geometry in direction 3-\n");
return(-1);
    }
    stest[iz3] += 1;
    
    iy3 = g_iup[ix][3];
    iz3 = g_lexic2eosub[iy3];
    if(iz3 > VOLUMEPLUSRAND/2 || iz3 < 0) {
      printf("There is a problem with EO geometry in direction 3+\n");
return(-1);
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
return(-1);
  }

  for(j = VOLUME/2; j < (VOLUME+RAND)/2; j++) {
    if(stest[j] != 1) {
      printf("There is a problem in the first boundary of the even odd geometry\n");
return(-1);
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
return(-1);
    }
    stest[iz0] += 1;
    
    iy0 = g_iup[ix][0];
    iz0 = g_lexic2eosub[iy0];
    if(iz0 > VOLUMEPLUSRAND/2 || iz0 < 0) {
      printf("There is a problem with EO geometry in direction 0+\n");
return(-1);
    }
    stest[iz0] += 1;
    
    iy1 = g_idn[ix][1];
    iz1 = g_lexic2eosub[iy1];
    if(iz1 > VOLUMEPLUSRAND/2  || iz1 < 0) {
      printf("There is a problem with EO geometry in direction 1-\n");
return(-1);
    }
    stest[iz1] += 1;
    
    iy1 = g_iup[ix][1];
    iz1 = g_lexic2eosub[iy1];
    if(iz1 > VOLUMEPLUSRAND/2  || iz1 < 0) {
      printf("There is a problem with EO geometry in direction 1+\n"); 
return(-1);
    }
    stest[iz1] += 1; 
    
    iy2 = g_idn[ix][2];
    iz2 = g_lexic2eosub[iy2];
    if(iz2 > VOLUMEPLUSRAND/2 || iz2 < 0) {
      printf("There is a problem with EO geometry in direction 2-\n");
return(-1);
    }
    stest[iz2] += 1;
    
    iy2 = g_iup[ix][2];
    iz2 = g_lexic2eosub[iy2];
    if(iz2 > VOLUMEPLUSRAND/2 || iz2 < 0) {
      printf("There is a problem with EO geometry in direction 2+\n");
return(-1);
    }
    stest[iz2] += 1;
    
    
    iy3 = g_idn[ix][3];
    iz3 = g_lexic2eosub[iy3];
    if(iz3 > VOLUMEPLUSRAND/2 || iz3 < 0) {
      printf("There is a problem with EO geometry in direction 3-\n");
return(-1);
    }
    stest[iz3] += 1;
    
    iy3 = g_iup[ix][3];
    iz3 = g_lexic2eosub[iy3];
    if(iz3 > VOLUMEPLUSRAND/2 || iz3 < 0) {
      printf("There is a problem with EO geometry in direction 3+\n");
return(-1);
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
return(-1);
  }

  for(j = VOLUME/2; j < (VOLUME+RAND)/2; j++) {
    if(stest[j] != 1) {
      printf("There is a problem in the second boundary of the even odd geometry\n");
return(-1);
    }
  }

  if(g_proc_id == 0 ) {
    printf("# The lattice is correctly mapped by the index arrays\n\n");
  }
  fflush(stdout);

  free(stest);
  free(itest);

  return(0);
}

#endif /* _INDEX_INDEP_GEOM */

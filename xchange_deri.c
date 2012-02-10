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
 * 
 * exchange routines for gauge fields
 *
 * Author: Carsten Urbach 
 *
 **********************************************************/


#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#ifdef MPI
# include <mpi.h>
#endif
#include "global.h"
#include "mpi_init.h"
#include "su3.h"
#include "su3adj.h"
#include "xchange_deri.h"

/* this if statement will be removed in future and _INDEX_INDEP_GEOM will be the default */
#ifdef _INDEX_INDEP_GEOM

void xchange_deri(su3adj ** const df)
{
#  ifdef MPI
  int ix,mu, t, y, z, x;
  MPI_Status status;

#    if (defined PARALLELT || defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT )
  /* send the data to the neighbour on the left in time direction */
  /* recieve the data from the neighbour on the right in time direction */
  MPI_Sendrecv(&df[gI_m1_0_0_0][0].d1,    1, deri_time_slice_cont, g_nb_t_dn, 43,
	       &ddummy[gI_Lm1_0_0_0][0].d1, 1, deri_time_slice_cont, g_nb_t_up, 43,
	       g_cart_grid, &status);

  /* add ddummy to df */
  for(x = 0; x < LX; x++) {
    for(y = 0; y < LY; y++) {
      for(z = 0; z < LZ; z++) {
	ix = g_ipt[T-1][x][y][z];
	for(mu=0;mu<4;mu++){ 
	  df[ix][mu].d1 += ddummy[ix][mu].d1;
	  df[ix][mu].d2 += ddummy[ix][mu].d2;
	  df[ix][mu].d3 += ddummy[ix][mu].d3;
	  df[ix][mu].d4 += ddummy[ix][mu].d4;
	  df[ix][mu].d5 += ddummy[ix][mu].d5;
	  df[ix][mu].d6 += ddummy[ix][mu].d6;
	  df[ix][mu].d7 += ddummy[ix][mu].d7;
	  df[ix][mu].d8 += ddummy[ix][mu].d8;
	}
      }
    }
  }

  /* send the data to the neighbour on the right is not needed*/
#    endif
#    if (defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT || defined PARALLELX || defined PARALLELXY || defined PARALLELXYZ )

  /* send the data to the neighbour on the left in x direction */
  /* recieve the data from the neighbour on the right in x direction */
  MPI_Sendrecv(&df[gI_0_m1_0_0][0],    1, deri_x_slice_cont, g_nb_x_dn, 44,
	       &ddummy[gI_0_Lm1_0_0][0],             1, deri_x_slice_gath, g_nb_x_up, 44,
	       g_cart_grid, &status);
  /* add ddummy to df */
  for(t = 0; t < T; t++) {
    for(y = 0; y < LY; y++) {
      for(z = 0; z < LZ; z++) {
	ix = g_ipt[t][LX-1][y][z];
	for(mu=0;mu<4;mu++){
	  df[ix][mu].d1 += ddummy[ix][mu].d1;
	  df[ix][mu].d2 += ddummy[ix][mu].d2;
	  df[ix][mu].d3 += ddummy[ix][mu].d3;
	  df[ix][mu].d4 += ddummy[ix][mu].d4;
	  df[ix][mu].d5 += ddummy[ix][mu].d5;
	  df[ix][mu].d6 += ddummy[ix][mu].d6;
	  df[ix][mu].d7 += ddummy[ix][mu].d7;
	  df[ix][mu].d8 += ddummy[ix][mu].d8;
	}
      }
    }
  }
  /* send the data to the neighbour on the right is not needed*/  
#    endif

#    if (defined PARALLELXYT || defined PARALLELXYZT || defined PARALLELXY || defined PARALLELXYZ )
  /* send the data to the neighbour on the left in y direction */
  /* recieve the data from the neighbour on the right in y direction */
  MPI_Sendrecv((void*)df[gI_0_0_m1_0], 
	       1, deri_y_slice_cont, g_nb_y_dn, 45,
	       (void*)ddummy[gI_0_0_Lm1_0],
	       1, deri_y_slice_gath, g_nb_y_up, 45,
	       g_cart_grid, &status);
  /* add ddummy to df */
  for(t = 0; t < T; t++) {
    for(x = 0; x < LX; x++) {
      for(z = 0; z < LZ; z++) {
	ix = g_ipt[t][x][LY-1][z];
	for(mu=0;mu<4;mu++){
	  df[ix][mu].d1 += ddummy[ix][mu].d1;
	  df[ix][mu].d2 += ddummy[ix][mu].d2;
	  df[ix][mu].d3 += ddummy[ix][mu].d3;
	  df[ix][mu].d4 += ddummy[ix][mu].d4;
	  df[ix][mu].d5 += ddummy[ix][mu].d5;
	  df[ix][mu].d6 += ddummy[ix][mu].d6;
	  df[ix][mu].d7 += ddummy[ix][mu].d7;
	  df[ix][mu].d8 += ddummy[ix][mu].d8;
	}
      }
    }
  }
  /* send the data to the neighbour on the right is not needed*/  
#    endif

#    if (defined PARALLELXYZT || defined PARALLELXYZ )
  /* send the data to the neighbour on the left in z direction */
  /* recieve the data from the neighbour on the right in z direction */
  MPI_Sendrecv((void*)df[gI_0_0_0_m1], 
	       1, deri_z_slice_cont, g_nb_z_dn, 46,
	       (void*)ddummy[gI_0_0_0_Lm1],
	       1, deri_z_slice_gath, g_nb_z_up, 46,
	       g_cart_grid, &status);
  /* add ddummy to df */
  for(t = 0; t < T; t++) {
    for(x = 0; x < LX; x++) {
      for(y = 0; y < LY; y++) {
	ix = g_ipt[t][x][y][LZ-1];
	for(mu=0;mu<4;mu++){
	  df[ix][mu].d1 += ddummy[ix][mu].d1;
	  df[ix][mu].d2 += ddummy[ix][mu].d2;
	  df[ix][mu].d3 += ddummy[ix][mu].d3;
	  df[ix][mu].d4 += ddummy[ix][mu].d4;
	  df[ix][mu].d5 += ddummy[ix][mu].d5;
	  df[ix][mu].d6 += ddummy[ix][mu].d6;
	  df[ix][mu].d7 += ddummy[ix][mu].d7;
	  df[ix][mu].d8 += ddummy[ix][mu].d8;
	}
      }
    }
  }
  /* send the data to the neighbour on the right is not needed*/  
#    endif
  return;
#  endif /* MPI */
}

#else /* _INDEX_INDEP_GEOM */

void xchange_deri(su3adj ** const df)
{
#  ifdef MPI
  int ix,mu, t, y, z, x;
  MPI_Status status;
  /* send the data to the neighbour on the left in time direction */
  /* recieve the data from the neighbour on the right in time direction */
  MPI_Sendrecv((void*)df[(T+1)*LX*LY*LZ],     1, deri_time_slice_cont, g_nb_t_dn, 40,
	       (void*)ddummy[(T-1)*LX*LY*LZ], 1, deri_time_slice_cont, g_nb_t_up, 40,
	       g_cart_grid, &status);

  /* add ddummy to df */
  for(x = 0; x < LX; x++) {
    for(y = 0; y < LY; y++) {
      for(z = 0; z < LZ; z++) {
	ix = g_ipt[T-1][x][y][z];
	for(mu=0;mu<4;mu++){ 
	  df[ix][mu].d1 += ddummy[ix][mu].d1;
	  df[ix][mu].d2 += ddummy[ix][mu].d2;
	  df[ix][mu].d3 += ddummy[ix][mu].d3;
	  df[ix][mu].d4 += ddummy[ix][mu].d4;
	  df[ix][mu].d5 += ddummy[ix][mu].d5;
	  df[ix][mu].d6 += ddummy[ix][mu].d6;
	  df[ix][mu].d7 += ddummy[ix][mu].d7;
	  df[ix][mu].d8 += ddummy[ix][mu].d8;
	}
      }
    }
  }

  /* send the data to the neighbour on the right in time direction needed for clover */

  MPI_Sendrecv((void*)df[T*L*L*L], 1, deri_time_slice_cont, g_nb_t_up, 41,
	       (void*)ddummy[0],   1, deri_time_slice_cont, g_nb_t_dn, 41,
	       g_cart_grid, &status);

  /* add ddummy to df */
  for(x = 0; x < LX; x++) {
    for(y = 0; y < LY; y++) {
      for(z = 0; z < LZ; z++) {
	ix = g_ipt[0][x][y][z];
	for(mu = 0; mu < 4; mu++) {
	  df[ix][mu].d1+=ddummy[ix][mu].d1;
	  df[ix][mu].d2+=ddummy[ix][mu].d2;
	  df[ix][mu].d3+=ddummy[ix][mu].d3;
	  df[ix][mu].d4+=ddummy[ix][mu].d4;
	  df[ix][mu].d5+=ddummy[ix][mu].d5;
	  df[ix][mu].d6+=ddummy[ix][mu].d6;
	  df[ix][mu].d7+=ddummy[ix][mu].d7;
	  df[ix][mu].d8+=ddummy[ix][mu].d8;
	}
      }
    }
  }

#    if (defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT)

  /* send the data to the neighbour on the left in x direction */
  /* recieve the data from the neighbour on the right in x direction */
  MPI_Sendrecv((void*)df[(T+2)*LX*LY*LZ + T*LY*LZ], 1, deri_x_slice_cont, g_nb_x_dn, 42,
	       (void*)ddummy[(LX-1)*LY*LZ],         1, deri_x_slice_gath, g_nb_x_up, 42,
	       g_cart_grid, &status);
  /* add ddummy to df */
  for(t = 0; t < T; t++) {
    for(y = 0; y < LY; y++) {
      for(z = 0; z < LZ; z++) {
	ix = g_ipt[t][LX-1][y][z];
	for(mu = 0; mu < 4; mu++) {
	  df[ix][mu].d1 += ddummy[ix][mu].d1;
	  df[ix][mu].d2 += ddummy[ix][mu].d2;
	  df[ix][mu].d3 += ddummy[ix][mu].d3;
	  df[ix][mu].d4 += ddummy[ix][mu].d4;
	  df[ix][mu].d5 += ddummy[ix][mu].d5;
	  df[ix][mu].d6 += ddummy[ix][mu].d6;
	  df[ix][mu].d7 += ddummy[ix][mu].d7;
	  df[ix][mu].d8 += ddummy[ix][mu].d8;
	}
      }
    }
  }
  /* send the data to the neighbour on the right needed for clover */  

  MPI_Sendrecv((void*)df[(T+2)*LX*LY*LZ], 1, deri_x_slice_cont, g_nb_x_up, 43,
	       (void*)ddummy[0],          1, deri_x_slice_gath, g_nb_x_dn, 43,
	       g_cart_grid, &status);
  /* add ddummy to df */
  for(t = 0; t < T; t++) {
    for(y = 0; y < LY; y++) {
      for(z = 0; z < LZ; z++) {
	ix = g_ipt[t][0][y][z];
	for(mu = 0; mu < 4; mu++) {
	  df[ix][mu].d1 += ddummy[ix][mu].d1;
	  df[ix][mu].d2 += ddummy[ix][mu].d2;
	  df[ix][mu].d3 += ddummy[ix][mu].d3;
	  df[ix][mu].d4 += ddummy[ix][mu].d4;
	  df[ix][mu].d5 += ddummy[ix][mu].d5;
	  df[ix][mu].d6 += ddummy[ix][mu].d6;
	  df[ix][mu].d7 += ddummy[ix][mu].d7;
	  df[ix][mu].d8 += ddummy[ix][mu].d8;
	}
      }
    }
  }


  /* end of ifdef PARALLELXT || PARALLELXYT */
#    endif

#    if (defined PARALLELXYT || defined PARALLELXYZT)
  /* send the data to the neighbour on the left in y direction */
  /* recieve the data from the neighbour on the right in y direction */
  MPI_Sendrecv((void*)df[VOLUME + 2*LZ*(LX*LY + T*LY) + T*LX*LZ], 
	       1, deri_y_slice_cont, g_nb_y_dn, 44,
	       (void*)ddummy[(LY-1)*LZ],
	       1, deri_y_slice_gath, g_nb_y_up, 44,
	       g_cart_grid, &status);
  /* add ddummy to df */
  for(t = 0; t < T; t++) {
    for(x = 0; x < LX; x++) {
      for(z = 0; z < LZ; z++) {
	ix = g_ipt[t][x][LY-1][z];
	for(mu = 0; mu < 4; mu++) {
	  df[ix][mu].d1 += ddummy[ix][mu].d1;
	  df[ix][mu].d2 += ddummy[ix][mu].d2;
	  df[ix][mu].d3 += ddummy[ix][mu].d3;
	  df[ix][mu].d4 += ddummy[ix][mu].d4;
	  df[ix][mu].d5 += ddummy[ix][mu].d5;
	  df[ix][mu].d6 += ddummy[ix][mu].d6;
	  df[ix][mu].d7 += ddummy[ix][mu].d7;
	  df[ix][mu].d8 += ddummy[ix][mu].d8;
	}
      }
    }
  }
  /* send the data to the neighbour on the right needed for clover*/  

  MPI_Sendrecv((void*)df[VOLUME + 2*LZ*(LX*LY + T*LY)], 
	       1, deri_y_slice_cont, g_nb_y_up, 45,
	       (void*)ddummy[0],
	       1, deri_y_slice_gath, g_nb_y_dn, 45,
	       g_cart_grid, &status);
  /* add ddummy to df */
  for(t = 0; t < T; t++) {
    for(x = 0; x < LX; x++) {
      for(z = 0; z < LZ; z++) {
	ix = g_ipt[t][x][0][z];
	for(mu = 0; mu < 4; mu++) {
	  df[ix][mu].d1 += ddummy[ix][mu].d1;
	  df[ix][mu].d2 += ddummy[ix][mu].d2;
	  df[ix][mu].d3 += ddummy[ix][mu].d3;
	  df[ix][mu].d4 += ddummy[ix][mu].d4;
	  df[ix][mu].d5 += ddummy[ix][mu].d5;
	  df[ix][mu].d6 += ddummy[ix][mu].d6;
	  df[ix][mu].d7 += ddummy[ix][mu].d7;
	  df[ix][mu].d8 += ddummy[ix][mu].d8;
	}
      }
    }
  }


  /* end of ifdef PARALLELXYT */
#    endif

#    ifdef PARALLELXYZT
  /* send the data to the neighbour on the left in y direction */
  /* recieve the data from the neighbour on the right in y direction */
  MPI_Sendrecv((void*)df[VOLUME + 2*LX*LY*LZ + 2*T*LY*LZ + 2*T*LX*LZ + T*LX*LY], 
	       1, deri_z_slice_cont, g_nb_z_dn, 46,
	       (void*)ddummy[LZ-1],
	       1, deri_z_slice_gath, g_nb_z_up, 46,
	       g_cart_grid, &status);
  /* add ddummy to df */
  for(t = 0; t < T; t++) {
    for(x = 0; x < LX; x++) {
      for(y = 0; y < LY; y++) {
	ix = g_ipt[t][x][y][LZ-1];
	for(mu=0;mu<4;mu++){
	  df[ix][mu].d1 += ddummy[ix][mu].d1;
	  df[ix][mu].d2 += ddummy[ix][mu].d2;
	  df[ix][mu].d3 += ddummy[ix][mu].d3;
	  df[ix][mu].d4 += ddummy[ix][mu].d4;
	  df[ix][mu].d5 += ddummy[ix][mu].d5;
	  df[ix][mu].d6 += ddummy[ix][mu].d6;
	  df[ix][mu].d7 += ddummy[ix][mu].d7;
	  df[ix][mu].d8 += ddummy[ix][mu].d8;
	}
      }
    }
  }
  /* send the data to the neighbour on the right is not needed*/  

  MPI_Sendrecv((void*)df[VOLUME + 2*LX*LY*LZ + 2*T*LY*LZ + 2*T*LX*LZ], 
	       1, deri_z_slice_cont, g_nb_z_up, 47,
	       (void*)ddummy[0],
	       1, deri_z_slice_gath, g_nb_z_dn, 47,
	       g_cart_grid, &status);
  /* add ddummy to df */
  for(t = 0; t < T; t++) {
    for(x = 0; x < LX; x++) {
      for(y = 0; y < LY; y++) {
	ix = g_ipt[t][x][y][0];
	for(mu=0;mu<4;mu++){
	  df[ix][mu].d1 += ddummy[ix][mu].d1;
	  df[ix][mu].d2 += ddummy[ix][mu].d2;
	  df[ix][mu].d3 += ddummy[ix][mu].d3;
	  df[ix][mu].d4 += ddummy[ix][mu].d4;
	  df[ix][mu].d5 += ddummy[ix][mu].d5;
	  df[ix][mu].d6 += ddummy[ix][mu].d6;
	  df[ix][mu].d7 += ddummy[ix][mu].d7;
	  df[ix][mu].d8 += ddummy[ix][mu].d8;
	}
      }
    }
  }

  /* end of ifdef PARALLELXYT */
#    endif
  return;
#  endif
}

#endif /* _INDEX_INDEP_GEOM */


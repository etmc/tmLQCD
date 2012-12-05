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

inline void addup_ddummy(su3adj** const df, const int ix, const int iy) {
  for(int mu = 0; mu < 4; mu++) {
    df[ix][mu].d1 += ddummy[iy][mu].d1;
    df[ix][mu].d2 += ddummy[iy][mu].d2;
    df[ix][mu].d3 += ddummy[iy][mu].d3;
    df[ix][mu].d4 += ddummy[iy][mu].d4;
    df[ix][mu].d5 += ddummy[iy][mu].d5;
    df[ix][mu].d6 += ddummy[iy][mu].d6;
    df[ix][mu].d7 += ddummy[iy][mu].d7;
    df[ix][mu].d8 += ddummy[iy][mu].d8;
  }
  return;
}

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

  /* send the data to the neighbour on the right is needed for the */
  /* clover case, so this needs fixing here! */
#    endif /* (defined PARALLELT || defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT ) */
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
  /* send the data to the neighbour on the right is needed for the */
  /* clover case, so this needs fixing here! */
#    endif /* (defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT || defined PARALLELX || defined PARALLELXY || defined PARALLELXYZ ) */

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
  /* send the data to the neighbour on the right is needed for the */
  /* clover case, so this needs fixing here! */
#    endif /* (defined PARALLELXYT || defined PARALLELXYZT || defined PARALLELXY || defined PARALLELXYZ ) */

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
  /* send the data to the neighbour on the right is needed for the */
  /* clover case, so this needs fixing here! */
#    endif /* (defined PARALLELXYZT || defined PARALLELXYZ ) */
#  endif /* MPI */
  return;
}

#else /* _INDEX_INDEP_GEOM */

void xchange_deri(su3adj ** const df)
{
#  ifdef MPI
  int ix,iy, t, y, z, x;
  MPI_Status status;
#    if (defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT)
  /* The edges need to come first */

  /* send the data to the neighbour on the left in t direction */
  /* recieve the data from the neighbour on the right in t direction */
  /* is on the x-boundary: xt-edge */
  MPI_Sendrecv((void*)df[VOLUME + RAND + 2*LY*LZ], 1, deri_xt_edge_cont, g_nb_t_dn, 492,
	       (void*)ddummy[0],                   1, deri_xt_edge_cont, g_nb_t_up, 492, 
	       g_cart_grid, &status);

  /* add ddummy to df */
  for(y = 0; y < LY; y++) {
    for(z = 0; z < LZ; z++) {
      ix = g_iup[ g_ipt[T-1][LX-1][y][z] ][1];
      iy = y*LZ + z;
      addup_ddummy(df, ix, iy);

      ix = g_idn[ g_ipt[T-1][0][y][z] ][1];
      iy = LY*LZ + y*LZ + z;
      addup_ddummy(df, ix, iy);
    }
  }

  /* send the data to the neighbour on the right in t direction */
  /* recieve the data from the neighbour on the left in t direction */
  /* xt-edge */
  MPI_Sendrecv((void*)df[VOLUME + RAND], 1, deri_xt_edge_cont, g_nb_t_up, 493,
	       (void*)ddummy[0],         1, deri_xt_edge_cont, g_nb_t_dn, 493,
	       g_cart_grid, &status);

  /* add ddummy to df */
  for(y = 0; y < LY; y++) {
    for(z = 0; z < LZ; z++) {
      ix = g_iup[ g_ipt[0][LX-1][y][z] ][1];
      iy = y*LZ + z;
      addup_ddummy(df, ix, iy);

      ix = g_idn[ g_ipt[0][0][y][z] ][1];
      iy = LY*LZ + y*LZ + z;
      addup_ddummy(df, ix, iy);
    }
  }

#    endif /* (defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT) */

#    if (defined PARALLELXYT || defined PARALLELXYZT)
  /* edges */

  /* send the data to the neighbour on the left in x direction */
  /* recieve the data from the neighbour on the right in x direction */
  /* is on the y-Rand -> yx-edge*/
  MPI_Sendrecv((void*)df[VOLUME + RAND + 4*LZ*LY + 2*T*LZ], 1, deri_yx_edge_cont, g_nb_x_dn, 494,
	       (void*)ddummy[0],                            1, deri_yx_edge_cont, g_nb_x_up, 494, 
	       g_cart_grid, &status);

  /* add ddummy to df */
  for(t = 0; t < T; t++) {
    for(z = 0; z < LZ; z++) {
      ix = g_iup[ g_ipt[t][LX-1][LY-1][z] ][2];
      iy = t*LZ + z;
      addup_ddummy(df, ix, iy);

      ix = g_idn[ g_ipt[t][LX-1][0][z] ][2];
      iy = T*LZ + t*LZ + z;
      addup_ddummy(df, ix, iy);
    }
  }


  /* send the data to the neighbour on the right in x direction */
  /* recieve the data from the neighbour on the left in x direction */
  /* yx-edge */
  MPI_Sendrecv((void*)df[VOLUME + RAND + 4*LZ*LY], 1, deri_yx_edge_cont, g_nb_x_up, 495,
	       (void*)ddummy[0],                   1, deri_yx_edge_cont, g_nb_x_dn, 495,
	       g_cart_grid, &status);

  /* add ddummy to df */
  for(t = 0; t < T; t++) {
    for(z = 0; z < LZ; z++) {
      ix = g_iup[ g_ipt[t][0][LY-1][z] ][2];
      iy = t*LZ + z;
      addup_ddummy(df, ix, iy);

      ix = g_idn[ g_ipt[t][0][0][z] ][2];
      iy = T*LZ + t*LZ + z;
      addup_ddummy(df, ix, iy);
    }
  }

  /* send the data to the neighbour on the left in y direction */
  /* recieve the data from the neighbour on the right in y direction */
  /* is on the t-Rand -> ty-edge*/
  MPI_Sendrecv((void*)df[VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 2*LX*LZ], 1, deri_ty_edge_cont, g_nb_y_dn, 496,
	       (void*)ddummy[0],                                      1, deri_ty_edge_cont, g_nb_y_up, 496, 
	       g_cart_grid, &status);

  /* add ddummy to df */
  for(x = 0; x < LX; x++) {
    for(z = 0; z < LZ; z++) {
      ix = g_iup[ g_ipt[T-1][x][LY-1][z] ][0];
      iy = x*LZ + z;
      addup_ddummy(df, ix, iy);

      ix = g_idn[ g_ipt[0][x][LY-1][z] ][0];
      iy = LX*LZ + x*LZ + z;
      addup_ddummy(df, ix, iy);
    }
  }


  /* send the data to the neighbour on the right in y direction */
  /* recieve the data from the neighbour on the left in y direction */
  /* ty-edge */
  MPI_Sendrecv((void*)df[VOLUME + RAND + 4*LY*LZ + 4*T*LZ], 1, deri_ty_edge_cont, g_nb_y_up, 497,
	       (void*)ddummy[0],                            1, deri_ty_edge_cont, g_nb_y_dn, 497,
	       g_cart_grid, &status);

  /* add ddummy to df */
  for(x = 0; x < LX; x++) {
    for(z = 0; z < LZ; z++) {
      ix = g_iup[ g_ipt[T-1][x][0][z] ][0];
      iy = x*LZ + z;
      addup_ddummy(df, ix, iy);

      ix = g_idn[ g_ipt[0][x][0][z] ][0];
      iy = LX*LZ + x*LZ + z;
      addup_ddummy(df, ix, iy);
    }
  }

#    endif /* (defined PARALLELXYT || defined PARALLELXYZT) */

#    ifdef PARALLELXYZT

  /* send the data to the neighbour on the left in x direction */
  /* recieve the data from the neighbour on the right in x direction */
  /* xz-edge */
  MPI_Sendrecv((void*)df[VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 4*LX*LZ + 2*T*LY], 
	       1, deri_zx_edge_cont, g_nb_x_dn, 498,
	       (void*)ddummy[0],                                      
	       1, deri_zx_edge_cont, g_nb_x_up, 498, 
	       g_cart_grid, &status);

  /* add ddummy to df */
  for(t = 0; t < T; t++) {
    for(y = 0; y < LY; y++) {
      ix = g_iup[ g_ipt[t][LX-1][y][LZ-1] ][3];
      iy = t*LY + y;
      addup_ddummy(df, ix, iy);

      ix = g_idn[ g_ipt[t][LX-1][y][0] ][3];
      iy = T*LY + t*LY + y;
      addup_ddummy(df, ix, iy);
    }
  }


  /* send the data to the neighbour on the right in x direction */
  /* recieve the data from the neighbour on the left in x direction */
  /* xz-edge */
  MPI_Sendrecv((void*)df[VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 4*LX*LZ], 
	       1, deri_zx_edge_cont, g_nb_x_up, 499,
	       (void*)ddummy[0],                            
	       1, deri_zx_edge_cont, g_nb_x_dn, 499,
	       g_cart_grid, &status);

  /* add ddummy to df */
  for(t = 0; t < T; t++) {
    for(y = 0; y < LY; y++) {
      ix = g_iup[ g_ipt[t][0][y][LZ-1] ][3];
      iy = t*LY + y;
      addup_ddummy(df, ix, iy);

      ix = g_idn[ g_ipt[t][0][y][0] ][3];
      iy = T*LY + t*LY + y;
      addup_ddummy(df, ix, iy);
    }
  }

  /* send the data to the neighbour on the left in z direction */
  /* recieve the data from the neighbour on the right in z direction */
  /* tz-edge */
  MPI_Sendrecv((void*)df[VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 4*LX*LZ + 4*T*LY + 2*LX*LY], 
	       1, deri_tz_edge_cont, g_nb_z_dn, 500,
	       (void*)ddummy[0],
	       1, deri_tz_edge_cont, g_nb_z_up, 500, 
	       g_cart_grid, &status);

  /* add ddummy to df */
  for(x = 0; x < LX; x++) {
    for(y = 0; y < LY; y++) {
      ix = g_iup[ g_ipt[T-1][x][y][LZ-1] ][0];
      iy = x*LY + y;
      addup_ddummy(df, ix, iy);

      ix = g_idn[ g_ipt[0][x][y][LZ-1] ][0];
      iy = LX*LY + x*LY + y;
      addup_ddummy(df, ix, iy);
    }
  }


  /* send the data to the neighbour on the right in z direction */
  /* recieve the data from the neighbour on the left in z direction */
  /* tz-edge */
  MPI_Sendrecv((void*)df[VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 4*LX*LZ + 4*T*LY], 
	       1, deri_tz_edge_cont, g_nb_z_up, 501,
	       (void*)ddummy[0],                            
	       1, deri_tz_edge_cont, g_nb_z_dn, 501,
	       g_cart_grid, &status);

  /* add ddummy to df */
  for(x = 0; x < LX; x++) {
    for(y = 0; y < LY; y++) {
      ix = g_iup[ g_ipt[T-1][x][y][0] ][0];
      iy = x*LY + y;
      addup_ddummy(df, ix, iy);

      ix = g_idn[ g_ipt[0][x][y][0] ][0];
      iy = LX*LY + x*LY + y;
      addup_ddummy(df, ix, iy);
    }
  }

  /* send the data to the neighbour on the left in y direction */
  /* recieve the data from the neighbour on the right in y direction */
  /* zy-edge */
  MPI_Sendrecv((void*)df[VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 4*LX*LZ + 4*T*LY + 4*LX*LY + 2*T*LX], 
	       1, deri_zy_edge_cont, g_nb_y_dn, 502,
	       (void*)ddummy[0],
	       1, deri_zy_edge_cont, g_nb_y_up, 502, 
	       g_cart_grid, &status);

  /* add ddummy to df */
  for(t = 0; t < T; t++) {
    for(x = 0; x < LX; x++) {
      ix = g_iup[ g_ipt[t][x][LY-1][LZ-1] ][3];
      iy = t*LX + x;
      addup_ddummy(df, ix, iy);

      ix = g_idn[ g_ipt[t][x][LY-1][0] ][3];
      iy = T*LX + t*LX + x;
      addup_ddummy(df, ix, iy);
    }
  }


  /* send the data to the neighbour on the right in y direction */
  /* recieve the data from the neighbour on the left in y direction */
  /* zy-edge */
  MPI_Sendrecv((void*)df[VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 4*LX*LZ + 4*T*LY + 4*LX*LY], 
	       1, deri_zy_edge_cont, g_nb_y_up, 503,
	       (void*)ddummy[0],
	       1, deri_zy_edge_cont, g_nb_y_dn, 503,
	       g_cart_grid, &status);

  /* add ddummy to df */
  for(t = 0; t < T; t++) {
    for(x = 0; x < LX; x++) {
      ix = g_iup[ g_ipt[t][x][0][LZ-1] ][3];
      iy = t*LX + x;
      addup_ddummy(df, ix, iy);

      ix = g_idn[ g_ipt[t][x][0][0] ][3];
      iy = T*LX + t*LX + x;
      addup_ddummy(df, ix, iy);
    }
  }

#    endif /* PARALLELXYZT */

  // now the normal boundaries

  /* send the data to the neighbour on the left in time direction */
  /* recieve the data from the neighbour on the right in time direction */
  MPI_Sendrecv((void*)df[(T+1)*LX*LY*LZ],     1, deri_time_slice_cont, g_nb_t_dn, 40,
	       (void*)ddummy[0],              1, deri_time_slice_cont, g_nb_t_up, 40,
	       g_cart_grid, &status);

  /* add ddummy to df */
  for(x = 0; x < LX; x++) {
    for(y = 0; y < LY; y++) {
      for(z = 0; z < LZ; z++) {
	ix = g_ipt[T-1][x][y][z];
	iy = x*LY*LZ + y*LZ + z;
	addup_ddummy(df, ix, iy);
      }
    }
  }

  /* send the data to the neighbour on the right in time direction needed for clover */

  MPI_Sendrecv((void*)df[T*LX*LY*LZ], 1, deri_time_slice_cont, g_nb_t_up, 41,
	       (void*)ddummy[0],      1, deri_time_slice_cont, g_nb_t_dn, 41,
	       g_cart_grid, &status);

  /* add ddummy to df */
  for(x = 0; x < LX; x++) {
    for(y = 0; y < LY; y++) {
      for(z = 0; z < LZ; z++) {
	ix = g_ipt[0][x][y][z];
	iy = x*LY*LZ + y*LZ + z;
	addup_ddummy(df, ix, iy);
      }
    }
  }



#    if (defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT)
  /* send the data to the neighbour on the left in x direction */
  /* recieve the data from the neighbour on the right in x direction */
  MPI_Sendrecv((void*)df[(T+2)*LX*LY*LZ + T*LY*LZ], 1, deri_x_slice_cont, g_nb_x_dn, 42,
	       (void*)ddummy[0],         1, deri_x_slice_cont, g_nb_x_up, 42,
	       g_cart_grid, &status);
  /* add ddummy to df */
  for(t = 0; t < T; t++) {
    for(y = 0; y < LY; y++) {
      for(z = 0; z < LZ; z++) {
	ix = g_ipt[t][LX-1][y][z];
	iy = t*LY*LZ + y*LZ + z;
	addup_ddummy(df, ix, iy);
      }
    }
  }

  /* send the data to the neighbour on the right needed for clover */  
  /* and receive from the one on the left                          */
  MPI_Sendrecv((void*)df[(T+2)*LX*LY*LZ], 1, deri_x_slice_cont, g_nb_x_up, 43,
	       (void*)ddummy[0],          1, deri_x_slice_cont, g_nb_x_dn, 43,
	       g_cart_grid, &status);
  /* add ddummy to df */
  for(t = 0; t < T; t++) {
    for(y = 0; y < LY; y++) {
      for(z = 0; z < LZ; z++) {
	ix = g_ipt[t][0][y][z];
	iy = t*LY*LZ + y*LZ + z;
	addup_ddummy(df, ix, iy);
      }
    }
  }

#    endif /* (defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT) */


#    if (defined PARALLELXYT || defined PARALLELXYZT)

  /* send the data to the neighbour on the left in y direction */
  /* recieve the data from the neighbour on the right in y direction */
  MPI_Sendrecv((void*)df[VOLUME + 2*LZ*(LX*LY + T*LY) + T*LX*LZ], 
	       1, deri_y_slice_cont, g_nb_y_dn, 44,
	       (void*)ddummy[0],
	       1, deri_y_slice_cont, g_nb_y_up, 44,
	       g_cart_grid, &status);
  /* add ddummy to df */
  for(t = 0; t < T; t++) {
    for(x = 0; x < LX; x++) {
      for(z = 0; z < LZ; z++) {
	ix = g_ipt[t][x][LY-1][z];
	iy = t*LX*LZ + x*LZ + z;
	addup_ddummy(df, ix, iy);
      }
    }
  }
  /* send the data to the neighbour on the right needed for clover*/  

  MPI_Sendrecv((void*)df[VOLUME + 2*LZ*(LX*LY + T*LY)], 
	       1, deri_y_slice_cont, g_nb_y_up, 45,
	       (void*)ddummy[0],
	       1, deri_y_slice_cont, g_nb_y_dn, 45,
	       g_cart_grid, &status);
  /* add ddummy to df */
  for(t = 0; t < T; t++) {
    for(x = 0; x < LX; x++) {
      for(z = 0; z < LZ; z++) {
	ix = g_ipt[t][x][0][z];
	iy = t*LX*LZ + x*LZ + z;
	addup_ddummy(df, ix, iy);
      }
    }
  }


#    endif /* (defined PARALLELXYT || defined PARALLELXYZT) */

#    ifdef PARALLELXYZT
  /* send the data to the neighbour on the left in y direction */
  /* recieve the data from the neighbour on the right in y direction */
  MPI_Sendrecv((void*)df[VOLUME + 2*LX*LY*LZ + 2*T*LY*LZ + 2*T*LX*LZ + T*LX*LY], 
	       1, deri_z_slice_cont, g_nb_z_dn, 46,
	       (void*)ddummy[0],
	       1, deri_z_slice_cont, g_nb_z_up, 46,
	       g_cart_grid, &status);
  /* add ddummy to df */
  for(t = 0; t < T; t++) {
    for(x = 0; x < LX; x++) {
      for(y = 0; y < LY; y++) {
	ix = g_ipt[t][x][y][LZ-1];
	iy = t*LX*LY + x*LY + y;
	addup_ddummy(df, ix, iy);
      }
    }
  }
  /* send the data to the neighbour on the right needed for clover */  

  MPI_Sendrecv((void*)df[VOLUME + 2*LX*LY*LZ + 2*T*LY*LZ + 2*T*LX*LZ], 
	       1, deri_z_slice_cont, g_nb_z_up, 47,
	       (void*)ddummy[0],
	       1, deri_z_slice_cont, g_nb_z_dn, 47,
	       g_cart_grid, &status);
  /* add ddummy to df */
  for(t = 0; t < T; t++) {
    for(x = 0; x < LX; x++) {
      for(y = 0; y < LY; y++) {
	ix = g_ipt[t][x][y][0];
	iy = t*LX*LY + x*LY + y;
	addup_ddummy(df, ix, iy);
      }
    }
  }

#    endif /* PARALLELXYZT */
#  endif /* MPI */
  return;
}

#endif /* _INDEX_INDEP_GEOM */


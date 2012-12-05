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

/**********************************************************
 * 
 * exchange routines for half spinor fields
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
#include <string.h>
#ifdef MPI
# include <mpi.h>
#endif
#include "global.h"
#if (defined XLC && defined BGL)
#  include "bgl.h"
#endif
#include "mpi_init.h"
#include "su3.h"
#include "init/init_dirac_halfspinor.h"
#include "xchange_halffield.h"

#if (defined _USE_HALFSPINOR)

#if (defined _PERSISTENT)

MPI_Request prequests[16];

/* 2. */
void init_xchange_halffield() {

#  ifdef MPI

#  ifdef PARALLELT
  int reqcount = 4;
#  elif defined PARALLELXT
  int reqcount = 8;
#  elif defined PARALLELXYT
  int reqcount = 12;
#  elif defined PARALLELXYZT
  int x0=0, x1=0, x2=0, ix=0;
  int reqcount = 16;
#  endif
#  if (defined XLC && defined BGL)
  __alignx(16, HalfSpinor);
#  endif

  /* send the data to the neighbour on the right in t direction */
  /* recieve the data from the neighbour on the left in t direction */
  MPI_Send_init((void*)(sendBuffer), LX*LY*LZ*12/2, MPI_DOUBLE, 
		g_nb_t_up, 81, g_cart_grid, &prequests[0]);
  
  MPI_Recv_init((void*)(recvBuffer + LX*LY*LZ/2), LX*LY*LZ*12/2, MPI_DOUBLE, 
		g_nb_t_dn, 81, g_cart_grid, &prequests[1]);

  /* send the data to the neighbour on the left in t direction */
  /* recieve the data from the neighbour on the right in t direction */
  MPI_Send_init((void*)(sendBuffer + LX*LY*LZ/2), LX*LY*LZ*12/2, MPI_DOUBLE, 
		g_nb_t_dn, 82, g_cart_grid, &prequests[2]);

  MPI_Recv_init((void*)(recvBuffer), LX*LY*LZ*12/2, MPI_DOUBLE, 
		g_nb_t_up, 82, g_cart_grid, &prequests[3]);

#    if (defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT)

  /* send the data to the neighbour on the right in x direction */
  /* recieve the data from the neighbour on the left in x direction */
  MPI_Send_init((void*)(sendBuffer + LX*LY*LZ), T*LY*LZ*12/2, MPI_DOUBLE, 
		g_nb_x_up, 91, g_cart_grid, &prequests[4]);

  MPI_Recv_init((void*)(recvBuffer + LX*LY*LZ + T*LY*LZ/2), T*LY*LZ*12/2, MPI_DOUBLE,
		g_nb_x_dn, 91, g_cart_grid, &prequests[5]);

  /* send the data to the neighbour on the left in x direction */
  /* recieve the data from the neighbour on the right in x direction */  
  MPI_Send_init((void*)(sendBuffer + LX*LY*LZ + T*LY*LZ/2), T*LY*LZ*12/2, MPI_DOUBLE,
		g_nb_x_dn, 92, g_cart_grid, &prequests[6]);

  MPI_Recv_init((void*)(recvBuffer + LX*LY*LZ), T*LY*LZ*12/2, MPI_DOUBLE,
		g_nb_x_up, 92, g_cart_grid, &prequests[7]);
#    endif
    
#    if (defined PARALLELXYT || defined PARALLELXYZT)
    /* send the data to the neighbour on the right in y direction */
    /* recieve the data from the neighbour on the left in y direction */
  MPI_Send_init((void*)(sendBuffer + LX*LY*LZ + T*LY*LZ), T*LX*LZ*12/2, MPI_DOUBLE, 
		g_nb_y_up, 101, g_cart_grid, &prequests[8]);

  MPI_Recv_init((void*)(recvBuffer + LX*LY*LZ + T*LY*LZ + T*LX*LZ/2), T*LX*LZ*12/2, MPI_DOUBLE, 
		g_nb_y_dn, 101, g_cart_grid, &prequests[9]);

    /* send the data to the neighbour on the leftt in y direction */
    /* recieve the data from the neighbour on the right in y direction */
  MPI_Send_init((void*)(sendBuffer + LX*LY*LZ + T*LY*LZ + T*LX*LZ/2), T*LX*LZ*12/2, MPI_DOUBLE, 
		g_nb_y_dn, 102, g_cart_grid, &prequests[10]);

  MPI_Recv_init((void*)(recvBuffer + LX*LY*LZ + T*LY*LZ), T*LX*LZ*12/2, MPI_DOUBLE, 
		g_nb_y_up, 102, g_cart_grid, &prequests[11]);
#    endif
    
#    if (defined PARALLELXYZT)
  /* send the data to the neighbour on the right in z direction */
  /* recieve the data from the neighbour on the left in z direction */
  MPI_Send_init((void*)(sendBuffer + LX*LY*LZ + T*LY*LZ + T*LX*LZ), 
		T*LX*LY*12/2, MPI_DOUBLE, g_nb_z_up, 503, g_cart_grid, &prequests[12]); 

  MPI_Recv_init((void*)(recvBuffer + LX*LY*LZ + T*LY*LZ + T*LX*LZ + T*LX*LY/2), 
		T*LX*LY*12/2, MPI_DOUBLE, g_nb_z_dn, 503, g_cart_grid, &prequests[13]); 

  /* send the data to the neighbour on the left in z direction */
  /* recieve the data from the neighbour on the right in z direction */
  MPI_Send_init((void*)(sendBuffer + LX*LY*LZ + T*LY*LZ + T*LX*LZ + T*LX*LY/2), 
		12*T*LX*LY/2, MPI_DOUBLE, g_nb_z_dn, 504, g_cart_grid, &prequests[14]); 

  MPI_Recv_init((void*)(recvBuffer + LX*LY*LZ + T*LY*LZ + T*LX*LZ), 
		T*LX*LY*12/2, MPI_DOUBLE, g_nb_z_up, 504, g_cart_grid, &prequests[15]); 
#  endif
#  endif /* MPI */
  return;
}

/* 3. */
void xchange_halffield() {
#  ifdef MPI

  MPI_Status status[16];
#    ifdef PARALLELT
  int reqcount = 4;
#    elif defined PARALLELXT
  int reqcount = 8;
#    elif defined PARALLELXYT
  int reqcount = 12;
#    elif defined PARALLELXYZT
  int x0=0, x1=0, x2=0, ix=0;
  int reqcount = 16;
#    endif
#    if (defined XLC && defined BGL)
  __alignx(16, HalfSpinor);
#    endif
  MPI_Startall(reqcount, prequests);

  MPI_Waitall(reqcount, prequests, status); 
#  endif /* MPI */
  return;
}

#else /* def (_USE_SHMEM || _PERSISTENT) */ 

# if defined _INDEX_INDEP_GEOM

/* 4. -IIG */
void xchange_halffield() {

#  ifdef MPI

  MPI_Request requests[16];
  MPI_Status status[16];
#  if ((defined PARALLELT) || (defined PARALLELX))
  int reqcount = 4;
#  elif ((defined PARALLELXT) || (defined PARALLELXY))
  int reqcount = 8;
#  elif ((defined PARALLELXYT) || (defined PARALLELXYZ))
  int reqcount = 12;
#  elif defined PARALLELXYZT
  int reqcount = 16;
#  endif
#  if (defined XLC && defined BGL)
  __alignx(16, HalfSpinor);
#  endif

#ifdef _KOJAK_INST
#pragma pomp inst begin(xchangehalf)
#endif

#    if (defined PARALLELT || defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT )
  /* send the data to the neighbour on the right in t direction */
  /* recieve the data from the neighbour on the left in t direction */
  MPI_Isend((void*)(sendBuffer + g_HS_shift_t), LX*LY*LZ*12/2, MPI_DOUBLE, 
	    g_nb_t_up, 81, g_cart_grid, &requests[0]);
  MPI_Irecv((void*)(recvBuffer + g_HS_shift_t + LX*LY*LZ/2), LX*LY*LZ*12/2, MPI_DOUBLE, 
	    g_nb_t_dn, 81, g_cart_grid, &requests[1]);
  /* send the data to the neighbour on the left in t direction */
  /* recieve the data from the neighbour on the right in t direction */
  MPI_Isend((void*)(sendBuffer + g_HS_shift_t + LX*LY*LZ/2), LX*LY*LZ*12/2, MPI_DOUBLE, 
	    g_nb_t_dn, 82, g_cart_grid, &requests[2]);
  MPI_Irecv((void*)(recvBuffer + g_HS_shift_t), LX*LY*LZ*12/2, MPI_DOUBLE, 
	    g_nb_t_up, 82, g_cart_grid, &requests[3]);
#    endif
#    if (defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT || defined PARALLELX || defined PARALLELXY || defined PARALLELXYZ )
  /* send the data to the neighbour on the right in x direction */
  /* recieve the data from the neighbour on the left in x direction */
  MPI_Isend((void*)(sendBuffer + g_HS_shift_x), T*LY*LZ*12/2, MPI_DOUBLE, 
	    g_nb_x_up, 91, g_cart_grid, &requests[4]);
  MPI_Irecv((void*)(recvBuffer + g_HS_shift_x + T*LY*LZ/2), T*LY*LZ*12/2, MPI_DOUBLE,
	    g_nb_x_dn, 91, g_cart_grid, &requests[5]);
  /* send the data to the neighbour on the left in x direction */
  /* recieve the data from the neighbour on the right in x direction */  
  MPI_Isend((void*)(sendBuffer + g_HS_shift_x + T*LY*LZ/2), T*LY*LZ*12/2, MPI_DOUBLE,
 	    g_nb_x_dn, 92, g_cart_grid, &requests[6]);
  MPI_Irecv((void*)(recvBuffer + g_HS_shift_x), T*LY*LZ*12/2, MPI_DOUBLE,
 	    g_nb_x_up, 92, g_cart_grid, &requests[7]);
#    endif
#    if (defined PARALLELXYT || defined PARALLELXYZT || defined PARALLELXY || defined PARALLELXYZ )
    /* send the data to the neighbour on the right in y direction */
    /* recieve the data from the neighbour on the left in y direction */
  MPI_Isend((void*)(sendBuffer + g_HS_shift_y), T*LX*LZ*12/2, MPI_DOUBLE, 
	    g_nb_y_up, 101, g_cart_grid, &requests[8]);
  MPI_Irecv((void*)(recvBuffer + g_HS_shift_y + T*LX*LZ/2), T*LX*LZ*12/2, MPI_DOUBLE, 
	    g_nb_y_dn, 101, g_cart_grid, &requests[9]);
    /* send the data to the neighbour on the leftt in y direction */
    /* recieve the data from the neighbour on the right in y direction */
  MPI_Isend((void*)(sendBuffer + g_HS_shift_y + T*LX*LZ/2), T*LX*LZ*12/2, MPI_DOUBLE, 
	    g_nb_y_dn, 102, g_cart_grid, &requests[10]);
  MPI_Irecv((void*)(recvBuffer + g_HS_shift_y), T*LX*LZ*12/2, MPI_DOUBLE, 
	    g_nb_y_up, 102, g_cart_grid, &requests[11]);
#    endif
#    if (defined PARALLELXYZT || defined PARALLELXYZ )
  /* send the data to the neighbour on the right in z direction */
  /* recieve the data from the neighbour on the left in z direction */
  MPI_Isend((void*)(sendBuffer + g_HS_shift_z), T*LX*LY*12/2, MPI_DOUBLE, 
		g_nb_z_up, 503, g_cart_grid, &requests[12]);
  MPI_Irecv((void*)(recvBuffer + g_HS_shift_z + T*LX*LY/2), T*LX*LY*12/2, MPI_DOUBLE, 
		g_nb_z_dn, 503, g_cart_grid, &requests[13]); 
  /* send the data to the neighbour on the left in z direction */
  /* recieve the data from the neighbour on the right in z direction */
  MPI_Isend((void*)(sendBuffer + g_HS_shift_z + T*LX*LY/2), 12*T*LX*LY/2, MPI_DOUBLE, 
		g_nb_z_dn, 504, g_cart_grid, &requests[14]);
  MPI_Irecv((void*)(recvBuffer + g_HS_shift_z), T*LX*LY*12/2, MPI_DOUBLE, 
		g_nb_z_up, 504, g_cart_grid, &requests[15]); 
#    endif

  MPI_Waitall(reqcount, requests, status); 
#  endif /* MPI */
  return;

#ifdef _KOJAK_INST
#pragma pomp inst end(xchangehalf)
#endif
}

# else /* _INDEX_INDEP_GEOM */

/* 4. */
void xchange_halffield() {

#  ifdef MPI

  MPI_Request requests[16];
  MPI_Status status[16];
#  ifdef PARALLELT
  int reqcount = 4;
#  elif defined PARALLELXT
  int reqcount = 8;
#  elif defined PARALLELXYT
  int reqcount = 12;
#  elif defined PARALLELXYZT
  int reqcount = 16;
#  endif
#  if (defined XLC && defined BGL)
  __alignx(16, HalfSpinor);
#  endif

#ifdef _KOJAK_INST
#pragma pomp inst begin(xchangehalf)
#endif
  /* send the data to the neighbour on the right in t direction */
  /* recieve the data from the neighbour on the left in t direction */
  MPI_Isend((void*)(sendBuffer), LX*LY*LZ*12/2, MPI_DOUBLE, 
	    g_nb_t_up, 81, g_cart_grid, &requests[0]);
  MPI_Irecv((void*)(recvBuffer + LX*LY*LZ/2), LX*LY*LZ*12/2, MPI_DOUBLE, 
	    g_nb_t_dn, 81, g_cart_grid, &requests[1]);

  /* send the data to the neighbour on the left in t direction */
  /* recieve the data from the neighbour on the right in t direction */
  MPI_Isend((void*)(sendBuffer+ LX*LY*LZ/2), LX*LY*LZ*12/2, MPI_DOUBLE, 
	    g_nb_t_dn, 82, g_cart_grid, &requests[2]);
  MPI_Irecv((void*)(recvBuffer), LX*LY*LZ*12/2, MPI_DOUBLE, 
	    g_nb_t_up, 82, g_cart_grid, &requests[3]);

#    if (defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT)

  /* send the data to the neighbour on the right in x direction */
  /* recieve the data from the neighbour on the left in x direction */
  MPI_Isend((void*)(sendBuffer + LX*LY*LZ), T*LY*LZ*12/2, MPI_DOUBLE, 
	    g_nb_x_up, 91, g_cart_grid, &requests[4]);
  MPI_Irecv((void*)(recvBuffer+ LX*LY*LZ + T*LY*LZ/2), T*LY*LZ*12/2, MPI_DOUBLE,
	    g_nb_x_dn, 91, g_cart_grid, &requests[5]);

  /* send the data to the neighbour on the left in x direction */
  /* recieve the data from the neighbour on the right in x direction */  
  MPI_Isend((void*)(sendBuffer + LX*LY*LZ + T*LY*LZ/2), T*LY*LZ*12/2, MPI_DOUBLE,
 	    g_nb_x_dn, 92, g_cart_grid, &requests[6]);
  MPI_Irecv((void*)(recvBuffer + LX*LY*LZ), T*LY*LZ*12/2, MPI_DOUBLE,
 	    g_nb_x_up, 92, g_cart_grid, &requests[7]);
#    endif
    
#    if (defined PARALLELXYT || defined PARALLELXYZT)
  /* send the data to the neighbour on the right in y direction */
  /* recieve the data from the neighbour on the left in y direction */
  MPI_Isend((void*)(sendBuffer + LX*LY*LZ + T*LY*LZ), T*LX*LZ*12/2, MPI_DOUBLE, 
	    g_nb_y_up, 101, g_cart_grid, &requests[8]);
  MPI_Irecv((void*)(recvBuffer + LX*LY*LZ + T*LY*LZ + T*LX*LZ/2), T*LX*LZ*12/2, MPI_DOUBLE, 
	    g_nb_y_dn, 101, g_cart_grid, &requests[9]);
  
  /* send the data to the neighbour on the leftt in y direction */
  /* recieve the data from the neighbour on the right in y direction */
  MPI_Isend((void*)(sendBuffer + LX*LY*LZ + T*LY*LZ + T*LX*LZ/2), T*LX*LZ*12/2, MPI_DOUBLE, 
	    g_nb_y_dn, 102, g_cart_grid, &requests[10]);
  MPI_Irecv((void*)(recvBuffer + LX*LY*LZ + T*LY*LZ), T*LX*LZ*12/2, MPI_DOUBLE, 
	    g_nb_y_up, 102, g_cart_grid, &requests[11]);
#    endif
  
#    if (defined PARALLELXYZT)
  /* send the data to the neighbour on the right in z direction */
  /* recieve the data from the neighbour on the left in z direction */
  MPI_Isend((void*)(sendBuffer + LX*LY*LZ + T*LY*LZ + T*LX*LZ), 
	    T*LX*LY*12/2, MPI_DOUBLE, g_nb_z_up, 503, g_cart_grid, &requests[12]);
  MPI_Irecv((void*)(recvBuffer + LX*LY*LZ + T*LY*LZ + T*LX*LZ + T*LX*LY/2), 
	    T*LX*LY*12/2, MPI_DOUBLE, g_nb_z_dn, 503, g_cart_grid, &requests[13]); 
  
  /* send the data to the neighbour on the left in z direction */
  /* recieve the data from the neighbour on the right in z direction */
  MPI_Isend((void*)(sendBuffer + LX*LY*LZ + T*LY*LZ + T*LX*LZ + T*LX*LY/2), 
	    12*T*LX*LY/2, MPI_DOUBLE, g_nb_z_dn, 504, g_cart_grid, &requests[14]);
  MPI_Irecv((void*)(recvBuffer + LX*LY*LZ + T*LY*LZ + T*LX*LZ), 
	    T*LX*LY*12/2, MPI_DOUBLE, g_nb_z_up, 504, g_cart_grid, &requests[15]); 
#    endif
  
  MPI_Waitall(reqcount, requests, status); 
#  endif /* MPI */
  return;
  
#ifdef _KOJAK_INST
#pragma pomp inst end(xchangehalf)
#endif
}

# endif /* _INDEX_INDEP_GEOM */

#endif /* def (_USE_SHMEM || _PERSISTENT) */ 


# if defined _INDEX_INDEP_GEOM
// IIG xchange_halffield32 still Missing
# else // defined _INDEX_INDEP_GEOM
/* 32-2. */
void xchange_halffield32() {

#  ifdef MPI

  MPI_Request requests[16];
  MPI_Status status[16];
#  ifdef PARALLELT
  int reqcount = 4;
#  elif defined PARALLELXT
  int reqcount = 8;
#  elif defined PARALLELXYT
  int reqcount = 12;
#  elif defined PARALLELXYZT
  int reqcount = 16;
#  endif
#ifdef _KOJAK_INST
#pragma pomp inst begin(xchangehalf32)
#endif
#  if (defined XLC && defined BGL)
  __alignx(16, HalfSpinor32);
#  endif

  /* send the data to the neighbour on the right in t direction */
  /* recieve the data from the neighbour on the left in t direction */
  MPI_Isend((void*)(sendBuffer32), LX*LY*LZ*12/2, MPI_FLOAT, 
	    g_nb_t_up, 81, g_cart_grid, &requests[0]);
  MPI_Irecv((void*)(recvBuffer32 + LX*LY*LZ/2), LX*LY*LZ*12/2, MPI_FLOAT, 
	    g_nb_t_dn, 81, g_cart_grid, &requests[1]);

  /* send the data to the neighbour on the left in t direction */
  /* recieve the data from the neighbour on the right in t direction */
  MPI_Isend((void*)(sendBuffer32 + LX*LY*LZ/2), LX*LY*LZ*12/2, MPI_FLOAT, 
	    g_nb_t_dn, 82, g_cart_grid, &requests[2]);
  MPI_Irecv((void*)(recvBuffer32), LX*LY*LZ*12/2, MPI_FLOAT, 
	    g_nb_t_up, 82, g_cart_grid, &requests[3]);

#    if (defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT)

  /* send the data to the neighbour on the right in x direction */
  /* recieve the data from the neighbour on the left in x direction */
  MPI_Isend((void*)(sendBuffer32 + LX*LY*LZ), T*LY*LZ*12/2, MPI_FLOAT, 
	    g_nb_x_up, 91, g_cart_grid, &requests[4]);
  MPI_Irecv((void*)(recvBuffer32 + LX*LY*LZ + T*LY*LZ/2), T*LY*LZ*12/2, MPI_FLOAT,
	    g_nb_x_dn, 91, g_cart_grid, &requests[5]);

  /* send the data to the neighbour on the left in x direction */
  /* recieve the data from the neighbour on the right in x direction */  
  MPI_Isend((void*)(sendBuffer32 + LX*LY*LZ + T*LY*LZ/2), T*LY*LZ*12/2, MPI_FLOAT,
 	    g_nb_x_dn, 92, g_cart_grid, &requests[6]);
  MPI_Irecv((void*)(recvBuffer32 + LX*LY*LZ), T*LY*LZ*12/2, MPI_FLOAT,
 	    g_nb_x_up, 92, g_cart_grid, &requests[7]);
#    endif
  
#    if (defined PARALLELXYT || defined PARALLELXYZT)
  /* send the data to the neighbour on the right in y direction */
  /* recieve the data from the neighbour on the left in y direction */
  MPI_Isend((void*)(sendBuffer32 + LX*LY*LZ + T*LY*LZ), T*LX*LZ*12/2, MPI_FLOAT, 
	    g_nb_y_up, 101, g_cart_grid, &requests[8]);
  MPI_Irecv((void*)(recvBuffer32 + LX*LY*LZ + T*LY*LZ + T*LX*LZ/2), T*LX*LZ*12/2, MPI_FLOAT, 
	    g_nb_y_dn, 101, g_cart_grid, &requests[9]);

  /* send the data to the neighbour on the leftt in y direction */
  /* recieve the data from the neighbour on the right in y direction */
  MPI_Isend((void*)(sendBuffer32 + LX*LY*LZ + T*LY*LZ + T*LX*LZ/2), T*LX*LZ*12/2, MPI_FLOAT, 
	    g_nb_y_dn, 102, g_cart_grid, &requests[10]);
  MPI_Irecv((void*)(recvBuffer32 + LX*LY*LZ + T*LY*LZ), T*LX*LZ*12/2, MPI_FLOAT, 
	    g_nb_y_up, 102, g_cart_grid, &requests[11]);
#    endif
    
#    if (defined PARALLELXYZT)
  /* send the data to the neighbour on the right in z direction */
  /* recieve the data from the neighbour on the left in z direction */
  MPI_Isend((void*)(sendBuffer32 + LX*LY*LZ + T*LY*LZ + T*LX*LZ), 
	    T*LX*LY*12/2, MPI_FLOAT, g_nb_z_up, 503, g_cart_grid, &requests[12]);
  MPI_Irecv((void*)(recvBuffer32 + LX*LY*LZ + T*LY*LZ + T*LX*LZ + T*LX*LY/2), 
	    T*LX*LY*12/2, MPI_FLOAT, g_nb_z_dn, 503, g_cart_grid, &requests[13]); 

  /* send the data to the neighbour on the left in z direction */
  /* recieve the data from the neighbour on the right in z direction */
  MPI_Isend((void*)(sendBuffer32 + LX*LY*LZ + T*LY*LZ + T*LX*LZ + T*LX*LY/2), 
	    12*T*LX*LY/2, MPI_FLOAT, g_nb_z_dn, 504, g_cart_grid, &requests[14]);
  MPI_Irecv((void*)(recvBuffer32 + LX*LY*LZ + T*LY*LZ + T*LX*LZ), 
	    T*LX*LY*12/2, MPI_FLOAT, g_nb_z_up, 504, g_cart_grid, &requests[15]); 
#    endif

  MPI_Waitall(reqcount, requests, status); 
#  endif /* MPI */
  return;
#ifdef _KOJAK_INST
#pragma pomp inst end(xchangehalf32)
#endif
}
# endif /* defined _INDEX_INDEP_GEOM */
#endif /* defined _USE_HALFSPINOR */













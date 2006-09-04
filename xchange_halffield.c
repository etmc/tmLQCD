/* $Id$ */

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
#ifdef _USE_SHMEM
# include <mpp/shmem.h>
#endif
#include "global.h"
#if (defined XLC && defined BGL)
#  include "bgl.h"
#endif
#include "mpi_init.h"
#include "su3.h"
#include "init_dirac_halfspinor.h"
#include "xchange_halffield.h"


#if (defined _USE_SHMEM && defined _USE_HALFSPINOR)
# include <mpp/shmem.h>
void xchange_halffield() {

#ifdef _KOJAK_INST
#pragma pomp inst begin(xchangehalf)
#endif
#  ifdef MPI

  shmem_barrier_all();

  shmem_double_put((double*)(HalfSpinor + 4*VOLUME + RAND/2 + LX*LY*LZ/2), 
		   (double*)(HalfSpinor + 4*VOLUME),
                   (LX*LY*LZ*6), g_nb_t_up);
  shmem_double_put((double*)(HalfSpinor + 4*VOLUME + RAND/2), 
		   (double*)(HalfSpinor + 4*VOLUME + LX*LY*LZ/2),
                   (LX*LY*LZ*6), g_nb_t_dn);

#    if (defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT)
  shmem_double_put((double*)(HalfSpinor + 4*VOLUME + RAND/2 + LX*LY*LZ + T*LY*LZ/2), 
		   (double*)(HalfSpinor + 4*VOLUME + LX*LY*LZ),
                   (T*LY*LZ*6), g_nb_x_up);
  shmem_double_put((double*)(HalfSpinor + 4*VOLUME + RAND/2 + LX*LY*LZ),
		   (double*)(HalfSpinor + 4*VOLUME + LX*LY*LZ + T*LY*LZ/2),
                   (T*LY*LZ*6), g_nb_x_dn);

#    endif

#    if (defined PARALLELXYT || defined PARALLELXYZT)
  shmem_double_put((double*)(HalfSpinor + 4*VOLUME + RAND/2 + LX*LY*LZ + T*LY*LZ + T*LX*LZ/2),
		   (double*)(HalfSpinor + 4*VOLUME + LX*LY*LZ + T*LY*LZ),
                   (T*LX*LZ*6), g_nb_y_up);
  shmem_double_put((double*)(HalfSpinor + 4*VOLUME + RAND/2 + LX*LY*LZ + T*LY*LZ),
		   (double*)(HalfSpinor + 4*VOLUME + LX*LY*LZ + T*LY*LZ + T*LX*LZ/2),
                   (T*LX*LZ*6), g_nb_y_dn);

#    endif

#    if (defined PARALLELXYZT)
  shmem_double_put((double*)(HalfSpinor + 4*VOLUME + RAND/2 + LX*LY*LZ + T*LY*LZ + T*LX*LZ + T*LX*LY/2),
		   (double*)(HalfSpinor + 4*VOLUME + LX*LY*LZ + T*LY*LZ + T*LX*LZ),
                   (T*LX*LY*6), g_nb_z_up);
  shmem_double_put((double*)(HalfSpinor + 4*VOLUME + RAND/2 + LX*LY*LZ + T*LY*LZ + T*LX*LZ),
		   (double*)(HalfSpinor + 4*VOLUME + LX*LY*LZ + T*LY*LZ + T*LX*LZ + T*LX*LY/2),
                   (T*LX*LY*6), g_nb_z_dn);

#    endif

  shmem_barrier_all();
#  endif
  return;
#ifdef _KOJAK_INST
#pragma pomp inst end(xchangehalf)
#endif
}
#else
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
  int x0=0, x1=0, x2=0, ix=0;
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
  MPI_Isend((void*)(HalfSpinor + 4*VOLUME), LX*LY*LZ*12/2, MPI_DOUBLE, 
	    g_nb_t_up, 81, g_cart_grid, &requests[0]);

  MPI_Irecv((void*)(HalfSpinor + 4*VOLUME + RAND/2 + LX*LY*LZ/2), LX*LY*LZ*12/2, MPI_DOUBLE, 
	    g_nb_t_dn, 81, g_cart_grid, &requests[1]);

  /* send the data to the neighbour on the left in t direction */
  /* recieve the data from the neighbour on the right in t direction */
  MPI_Isend((void*)(HalfSpinor + 4*VOLUME + LX*LY*LZ/2), LX*LY*LZ*12/2, MPI_DOUBLE, 
	    g_nb_t_dn, 82, g_cart_grid, &requests[2]);

  MPI_Irecv((void*)(HalfSpinor + 4*VOLUME + RAND/2), LX*LY*LZ*12/2, MPI_DOUBLE, 
	    g_nb_t_up, 82, g_cart_grid, &requests[3]);

#    if (defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT)

  /* send the data to the neighbour on the right in x direction */
  /* recieve the data from the neighbour on the left in x direction */
  MPI_Isend((void*)(HalfSpinor + 4*VOLUME + LX*LY*LZ), T*LY*LZ*12/2, MPI_DOUBLE, 
	    g_nb_x_up, 91, g_cart_grid, &requests[4]);

  MPI_Irecv((void*)(HalfSpinor + 4*VOLUME + RAND/2 + LX*LY*LZ + T*LY*LZ/2), T*LY*LZ*12/2, MPI_DOUBLE,
	    g_nb_x_dn, 91, g_cart_grid, &requests[5]);

  /* send the data to the neighbour on the left in x direction */
  /* recieve the data from the neighbour on the right in x direction */  
  MPI_Isend((void*)(HalfSpinor + 4*VOLUME + LX*LY*LZ + T*LY*LZ/2), T*LY*LZ*12/2, MPI_DOUBLE,
 	    g_nb_x_dn, 92, g_cart_grid, &requests[6]);

  MPI_Irecv((void*)(HalfSpinor + 4*VOLUME + RAND/2 + LX*LY*LZ), T*LY*LZ*12/2, MPI_DOUBLE,
 	    g_nb_x_up, 92, g_cart_grid, &requests[7]);
#    endif
    
#    if (defined PARALLELXYT || defined PARALLELXYZT)
    /* send the data to the neighbour on the right in y direction */
    /* recieve the data from the neighbour on the left in y direction */
  MPI_Isend((void*)(HalfSpinor + 4*VOLUME + LX*LY*LZ + T*LY*LZ), T*LX*LZ*12/2, MPI_DOUBLE, 
	    g_nb_y_up, 101, g_cart_grid, &requests[8]);

  MPI_Irecv((void*)(HalfSpinor + 4*VOLUME + RAND/2 + LX*LY*LZ + T*LY*LZ + T*LX*LZ/2), T*LX*LZ*12/2, MPI_DOUBLE, 
	    g_nb_y_dn, 101, g_cart_grid, &requests[9]);

    /* send the data to the neighbour on the leftt in y direction */
    /* recieve the data from the neighbour on the right in y direction */
  MPI_Isend((void*)(HalfSpinor + 4*VOLUME + LX*LY*LZ + T*LY*LZ + T*LX*LZ/2), T*LX*LZ*12/2, MPI_DOUBLE, 
	    g_nb_y_dn, 102, g_cart_grid, &requests[10]);

  MPI_Irecv((void*)(HalfSpinor + 4*VOLUME + RAND/2 + LX*LY*LZ + T*LY*LZ), T*LX*LZ*12/2, MPI_DOUBLE, 
	    g_nb_y_up, 102, g_cart_grid, &requests[11]);
#    endif
    
#    if (defined PARALLELXYZT)
  /* send the data to the neighbour on the right in z direction */
  /* recieve the data from the neighbour on the left in z direction */
  MPI_Isend((void*)(HalfSpinor + 4*VOLUME + LX*LY*LZ + T*LY*LZ + T*LX*LZ), 
	    T*LX*LY*12/2, MPI_DOUBLE, g_nb_z_up, 503, g_cart_grid, &requests[12]); 

  MPI_Irecv((void*)(HalfSpinor + 4*VOLUME + RAND/2 + LX*LY*LZ + T*LY*LZ + T*LX*LZ + T*LX*LY/2), 
	    T*LX*LY*12/2, MPI_DOUBLE, g_nb_z_dn, 503, g_cart_grid, &requests[13]); 

  /* send the data to the neighbour on the left in z direction */
  /* recieve the data from the neighbour on the right in z direction */
  MPI_Isend((void*)(HalfSpinor + 4*VOLUME + LX*LY*LZ + T*LY*LZ + T*LX*LZ + T*LX*LY/2), 
	    12*T*LX*LY/2, MPI_DOUBLE, g_nb_z_dn, 504, g_cart_grid, &requests[14]); 

  MPI_Irecv((void*)(HalfSpinor + 4*VOLUME + RAND/2 + LX*LY*LZ + T*LY*LZ + T*LX*LZ), 
	    T*LX*LY*12/2, MPI_DOUBLE, g_nb_z_up, 504, g_cart_grid, &requests[15]); 
#    endif

  MPI_Waitall(reqcount, requests, status); 
#  endif
  return;

#ifdef _KOJAK_INST
#pragma pomp inst end(xchangehalf)
#endif
}

#endif

#if (defined _USE_SHMEM && defined _USE_HALFSPINOR)
# include <mpp/shmem.h>
void xchange_halffield32() {

#ifdef _KOJAK_INST
#pragma pomp inst begin(xchangehalf32)
#endif
#  ifdef MPI

  shmem_barrier_all();

  shmem_float_put((float*)(HalfSpinor32 + 4*VOLUME + RAND/2 + LX*LY*LZ/2), 
		   (float*)(HalfSpinor32 + 4*VOLUME),
                   (LX*LY*LZ*6), g_nb_t_up);
  shmem_float_put((float*)(HalfSpinor32 + 4*VOLUME + RAND/2), 
		   (float*)(HalfSpinor32 + 4*VOLUME + LX*LY*LZ/2),
                   (LX*LY*LZ*6), g_nb_t_dn);

#    if (defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT)
  shmem_float_put((float*)(HalfSpinor32 + 4*VOLUME + RAND/2 + LX*LY*LZ + T*LY*LZ/2), 
		   (float*)(HalfSpinor32 + 4*VOLUME + LX*LY*LZ),
                   (T*LY*LZ*6), g_nb_x_up);
  shmem_float_put((float*)(HalfSpinor32 + 4*VOLUME + RAND/2 + LX*LY*LZ),
		   (float*)(HalfSpinor32 + 4*VOLUME + LX*LY*LZ + T*LY*LZ/2),
                   (T*LY*LZ*6), g_nb_x_dn);

#    endif

#    if (defined PARALLELXYT || defined PARALLELXYZT)
  shmem_float_put((float*)(HalfSpinor32 + 4*VOLUME + RAND/2 + LX*LY*LZ + T*LY*LZ + T*LX*LZ/2),
		   (float*)(HalfSpinor32 + 4*VOLUME + LX*LY*LZ + T*LY*LZ),
                   (T*LX*LZ*6), g_nb_y_up);
  shmem_float_put((float*)(HalfSpinor32 + 4*VOLUME + RAND/2 + LX*LY*LZ + T*LY*LZ),
		   (float*)(HalfSpinor32 + 4*VOLUME + LX*LY*LZ + T*LY*LZ + T*LX*LZ/2),
                   (T*LX*LZ*6), g_nb_y_dn);

#    endif

#    if (defined PARALLELXYZT)
  shmem_float_put((float*)(HalfSpinor32 + 4*VOLUME + RAND/2 + LX*LY*LZ + T*LY*LZ + T*LX*LZ + T*LX*LY/2),
		   (float*)(HalfSpinor32 + 4*VOLUME + LX*LY*LZ + T*LY*LZ + T*LX*LZ),
                   (T*LX*LY*6), g_nb_z_up);
  shmem_float_put((float*)(HalfSpinor32 + 4*VOLUME + RAND/2 + LX*LY*LZ + T*LY*LZ + T*LX*LZ),
		   (float*)(HalfSpinor32 + 4*VOLUME + LX*LY*LZ + T*LY*LZ + T*LX*LZ + T*LX*LY/2),
                   (T*LX*LY*6), g_nb_z_dn);

#    endif

  shmem_barrier_all();
#  endif
  return;
#ifdef _KOJAK_INST
#pragma pomp inst end(xchangehalf32)
#endif
}

#else

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
  int x0=0, x1=0, x2=0, ix=0;
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
  MPI_Isend((void*)(HalfSpinor32 + 4*VOLUME), LX*LY*LZ*12/2, MPI_FLOAT, 
	    g_nb_t_up, 81, g_cart_grid, &requests[0]);

  MPI_Irecv((void*)(HalfSpinor32 + 4*VOLUME + RAND/2 + LX*LY*LZ/2), LX*LY*LZ*12/2, MPI_FLOAT, 
	    g_nb_t_dn, 81, g_cart_grid, &requests[1]);

  /* send the data to the neighbour on the left in t direction */
  /* recieve the data from the neighbour on the right in t direction */
  MPI_Isend((void*)(HalfSpinor32 + 4*VOLUME + LX*LY*LZ/2), LX*LY*LZ*12/2, MPI_FLOAT, 
	    g_nb_t_dn, 82, g_cart_grid, &requests[2]);

  MPI_Irecv((void*)(HalfSpinor32 + 4*VOLUME + RAND/2), LX*LY*LZ*12/2, MPI_FLOAT, 
	    g_nb_t_up, 82, g_cart_grid, &requests[3]);

#    if (defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT)

  /* send the data to the neighbour on the right in x direction */
  /* recieve the data from the neighbour on the left in x direction */
  MPI_Isend((void*)(HalfSpinor32 + 4*VOLUME + LX*LY*LZ), T*LY*LZ*12/2, MPI_FLOAT, 
	    g_nb_x_up, 91, g_cart_grid, &requests[4]);

  MPI_Irecv((void*)(HalfSpinor32 + 4*VOLUME + RAND/2 + LX*LY*LZ + T*LY*LZ/2), T*LY*LZ*12/2, MPI_FLOAT,
	    g_nb_x_dn, 91, g_cart_grid, &requests[5]);

  /* send the data to the neighbour on the left in x direction */
  /* recieve the data from the neighbour on the right in x direction */  
  MPI_Isend((void*)(HalfSpinor32 + 4*VOLUME + LX*LY*LZ + T*LY*LZ/2), T*LY*LZ*12/2, MPI_FLOAT,
 	    g_nb_x_dn, 92, g_cart_grid, &requests[6]);

  MPI_Irecv((void*)(HalfSpinor32 + 4*VOLUME + RAND/2 + LX*LY*LZ), T*LY*LZ*12/2, MPI_FLOAT,
 	    g_nb_x_up, 92, g_cart_grid, &requests[7]);
#    endif
    
#    if (defined PARALLELXYT || defined PARALLELXYZT)
    /* send the data to the neighbour on the right in y direction */
    /* recieve the data from the neighbour on the left in y direction */
  MPI_Isend((void*)(HalfSpinor32 + 4*VOLUME + LX*LY*LZ + T*LY*LZ), T*LX*LZ*12/2, MPI_FLOAT, 
	    g_nb_y_up, 101, g_cart_grid, &requests[8]);

  MPI_Irecv((void*)(HalfSpinor32 + 4*VOLUME + RAND/2 + LX*LY*LZ + T*LY*LZ + T*LX*LZ/2), T*LX*LZ*12/2, MPI_FLOAT, 
	    g_nb_y_dn, 101, g_cart_grid, &requests[9]);

    /* send the data to the neighbour on the leftt in y direction */
    /* recieve the data from the neighbour on the right in y direction */
  MPI_Isend((void*)(HalfSpinor32 + 4*VOLUME + LX*LY*LZ + T*LY*LZ + T*LX*LZ/2), T*LX*LZ*12/2, MPI_FLOAT, 
	    g_nb_y_dn, 102, g_cart_grid, &requests[10]);

  MPI_Irecv((void*)(HalfSpinor32 + 4*VOLUME + RAND/2 + LX*LY*LZ + T*LY*LZ), T*LX*LZ*12/2, MPI_FLOAT, 
	    g_nb_y_up, 102, g_cart_grid, &requests[11]);
#    endif
    
#    if (defined PARALLELXYZT)
  /* send the data to the neighbour on the right in z direction */
  /* recieve the data from the neighbour on the left in z direction */
  MPI_Isend((void*)(HalfSpinor32 + 4*VOLUME + LX*LY*LZ + T*LY*LZ + T*LX*LZ), 
	    T*LX*LY*12/2, MPI_FLOAT, g_nb_z_up, 503, g_cart_grid, &requests[12]); 

  MPI_Irecv((void*)(HalfSpinor32 + 4*VOLUME + RAND/2 + LX*LY*LZ + T*LY*LZ + T*LX*LZ + T*LX*LY/2), 
	    T*LX*LY*12/2, MPI_FLOAT, g_nb_z_dn, 503, g_cart_grid, &requests[13]); 

  /* send the data to the neighbour on the left in z direction */
  /* recieve the data from the neighbour on the right in z direction */
  MPI_Isend((void*)(HalfSpinor32 + 4*VOLUME + LX*LY*LZ + T*LY*LZ + T*LX*LZ + T*LX*LY/2), 
	    12*T*LX*LY/2, MPI_FLOAT, g_nb_z_dn, 504, g_cart_grid, &requests[14]); 

  MPI_Irecv((void*)(HalfSpinor32 + 4*VOLUME + RAND/2 + LX*LY*LZ + T*LY*LZ + T*LX*LZ), 
	    T*LX*LY*12/2, MPI_FLOAT, g_nb_z_up, 504, g_cart_grid, &requests[15]); 
#    endif

  MPI_Waitall(reqcount, requests, status); 
#  endif
  return;
#ifdef _KOJAK_INST
#pragma pomp inst end(xchangehalf32)
#endif
}
#endif
static char const rcsid[] = "$Id$";

/* $Id$ */

/**********************************************************
 * 
 * exchange routines for spinor fields
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
#include "xchange_field.h"


/* exchanges the field  l */
void xchange_field(spinor * const l) {

#ifdef MPI

  MPI_Status status;
  /* send the data to the neighbour on the left */
  /* recieve the data from the neighbour on the right */
  MPI_Sendrecv((void*)l,                1, field_time_slice_cont, g_nb_t_dn, 81,
	       (void*)(l+T*LX*LY*LZ/2), 1, field_time_slice_cont, g_nb_t_up, 81,
	       g_cart_grid, &status);

  /* send the data to the neighbour on the right */
  /* recieve the data from the neighbour on the left */
  MPI_Sendrecv((void*)(l+(T-1)*LX*LY*LZ/2), 1, field_time_slice_cont, g_nb_t_up, 82,
	       (void*)(l+(T+1)*LX*LY*LZ/2), 1, field_time_slice_cont, g_nb_t_dn, 82,
	       g_cart_grid, &status);

# if (defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT)
  /* send the data to the neighbour on the left in x direction */
  /* recieve the data from the neighbour on the right in x direction */
  MPI_Sendrecv((void*)l,                    1, field_x_slice_gath, g_nb_x_dn, 91, 
 	       (void*)(l+(T+2)*LX*LY*LZ/2), 1, field_x_slice_cont, g_nb_x_up, 91,
	       g_cart_grid, &status);

  /* send the data to the neighbour on the right in x direction */
  /* recieve the data from the neighbour on the left in x direction */  
  MPI_Sendrecv((void*)(l+(LX-1)*LY*LZ/2),               1, field_x_slice_gath, g_nb_x_up, 92, 
	       (void*)(l+((T+2)*LX*LY*LZ + T*LY*LZ)/2), 1, field_x_slice_cont, g_nb_x_dn, 92,
	       g_cart_grid, &status);

# endif

# if (defined PARALLELXYT || defined PARALLELXYZT)
  /* send the data to the neighbour on the left in y direction */
  /* recieve the data from the neighbour on the right in y direction */
  MPI_Sendrecv((void*)l,                                  1, field_y_slice_gath, g_nb_y_dn, 101, 
 	       (void*)(l+((T+2)*LX*LY*LZ + 2*T*LY*LZ)/2), 1, field_y_slice_cont, g_nb_y_up, 101,
	       g_cart_grid, &status);

  /* send the data to the neighbour on the right in y direction */
  /* recieve the data from the neighbour on the left in y direction */  
  MPI_Sendrecv((void*)(l+(LY-1)*LZ/2),                              1, field_y_slice_gath, g_nb_y_up, 102, 
	       (void*)(l+((T+2)*LX*LY*LZ + 2*T*LY*LZ + T*LX*LZ)/2), 1, field_y_slice_cont, g_nb_y_dn, 102,
	       g_cart_grid, &status);

# endif

# if (defined PARALLELXYZT)
  /* send the data to the neighbour on the left in z direction */
  /* recieve the data from the neighbour on the right in z direction */
  MPI_Sendrecv((void*)l,
	       1, field_z_slice_gath, g_nb_z_dn, 103, 
 	       (void*)(l+((T+2)*LX*LY*LZ + 2*T*LY*LZ +2*T*LX*LZ)/2), 
	       1, field_z_slice_cont, g_nb_z_up, 103,
	       g_cart_grid, &status);

  /* send the data to the neighbour on the right in y direction */
  /* recieve the data from the neighbour on the left in y direction */  
  MPI_Sendrecv((void*)(l+(LZ-1)/2), 
	       1, field_z_slice_gath, g_nb_z_up, 104, 
	       (void*)(l+((T+2)*LX*LY*LZ + 2*T*LY*LZ + 2*T*LX*LZ + T*LX*LY)/2), 
	       1, field_z_slice_cont, g_nb_z_dn, 104,
	       g_cart_grid, &status);

# endif

#endif
  return;
}

static char const rcsid[] = "$Id$";

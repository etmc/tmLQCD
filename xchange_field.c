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
#include <string.h>
#ifdef MPI
# include <mpi.h>
#endif
#include "global.h"
#include "mpi_init.h"
#include "su3.h"
#include "xchange_field.h"


/* exchanges the field  l */
void xchange_field(spinor * const l, const int ieo) {

#ifdef PARALLELXYZT
  int x0=0, x1=0, x2=0, ix=0;
#endif
#ifdef MPI

#  if (defined BGL && defined XLC)
  __alignx(16, l);
#  endif

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
  /* fill buffer ! */
  /* This is now depending on whether the field is */
  /* even or odd */
  ix = 0;
  for(x0 = 0; x0 < T; x0++) {
    for(x1 = 0; x1 < LX; x1++) {
      for(x2 = 0; x2 < LY; x2++) {
 	if((x0 + x1 + x2 +  
 	    g_proc_coords[0]*T + g_proc_coords[1]*LX +  
 	    g_proc_coords[2]*LY + g_proc_coords[3]*LZ)%2==(ieo+1)%2) { 
	  memcpy((void*)(field_buffer_z+ix), (void*)&l[ g_lexic2eosub[ g_ipt[x0][x1][x2][0]] ], sizeof(spinor));
	  ix++;
	}
      }
    }
  }

  /* send the data to the neighbour on the left in z direction */
  /* recieve the data from the neighbour on the right in z direction */
  MPI_Sendrecv((void*)field_buffer_z,
	       1, field_z_slice_cont, g_nb_z_dn, 503, 
 	       (void*)(l+(VOLUME/2 + LX*LY*LZ + T*LY*LZ +T*LX*LZ)), 
	       1, field_z_slice_cont, g_nb_z_up, 503,
	       g_cart_grid, &status);

  ix = 0;
  for(x0 = 0; x0 < T; x0++) {
    for(x1 = 0; x1 < LX; x1++) {
      for(x2 = 0; x2 < LY; x2++) {
	if((x0 + x1 + x2 + (LZ-1) + 
	    g_proc_coords[0]*T + g_proc_coords[1]*LX +  
	    g_proc_coords[2]*LY + g_proc_coords[3]*LZ)%2==(ieo+1)%2) { 
	  memcpy((void*)(field_buffer_z+ix), (void*)&l[ g_lexic2eosub[ g_ipt[x0][x1][x2][LZ-1]] ], sizeof(spinor));
	  ix++;
	}
      }
    }
  }
  /* send the data to the neighbour on the right in y direction */
  /* recieve the data from the neighbour on the left in y direction */  
  MPI_Sendrecv((void*)field_buffer_z, 
	       1, field_z_slice_cont, g_nb_z_up, 504, 
	       (void*)(l+(VOLUME + 2*LX*LY*LZ + 2*T*LY*LZ + 2*T*LX*LZ + T*LX*LY)/2), 
	       1, field_z_slice_cont, g_nb_z_dn, 504,
	       g_cart_grid, &status);

# endif

#endif
  return;
}

static char const rcsid[] = "$Id$";

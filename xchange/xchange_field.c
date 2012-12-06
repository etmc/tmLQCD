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
#ifdef _USE_SHMEM
# include <mpp/shmem.h>
#endif

#include "global.h"
#if (defined XLC && defined BGL)
#  include "bgl.h"
#endif
#include "mpi_init.h"
#include "su3.h"
#include "xchange_field.h"

#if (defined XLC && defined PARALLELXYZT)
#pragma disjoint(*field_buffer_z2, *field_buffer_z)
#endif

/* this version uses non-blocking MPI calls */
#if (defined _NON_BLOCKING)

/* this is the version independent of the content of the function Index  */
/* this if statement will be removed in future and _INDEX_INDEP_GEOM will be the default */
# ifdef  _INDEX_INDEP_GEOM

void xchange_field(spinor * const l, const int ieo) {

#ifdef MPI
  MPI_Request requests[16];
  MPI_Status status[16];
#endif
  int ireq;
#  if ( defined PARALLELT || defined PARALLELX )
  int reqcount = 4;
#  elif ( defined PARALLELXT || defined PARALLELXY )
  int reqcount = 8;
#  elif ( defined PARALLELXYT || defined PARALLELXYZ )
  int reqcount = 12;
#  elif defined PARALLELXYZT
  int ix=0;
  int reqcount = 16;
#  endif

#ifdef _KOJAK_INST
#pragma pomp inst begin(xchangefield)
#endif
#  if (defined BGL && defined XLC)
  __alignx(16, l);
#  endif

#  ifdef MPI


  /* In 4 dimensions there are two processors sharing the   */
  /* communication bandwidth. So the first should start     */
  /* in forward direction, the second in backward direction */
  /* This might only work if the third direction is         */
  /* parallelised only on the node                          */
  if(g_proc_coords[3]%2 == 0) {

    ireq=0;

#    if (defined PARALLELT || defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT )
    /* send the data to the neighbour on the left */
    /* recieve the data from the neighbour on the right */
    MPI_Isend((void*)(l+g_1st_t_int_dn), 1, field_time_slice_cont, g_nb_t_dn, 81, g_cart_grid, &requests[ireq]);
    MPI_Irecv( (void*)(l+g_1st_t_ext_up), 1, field_time_slice_cont, g_nb_t_up, 81, g_cart_grid, &requests[ireq+1]);
    ireq=ireq+4;
#    endif
    
#    if (defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT || defined PARALLELX || defined PARALLELXY || defined PARALLELXYZ )
    /* send the data to the neighbour on the left in x direction */
    /* recieve the data from the neighbour on the right in x direction */
    MPI_Isend((void*)(l+g_1st_x_int_dn), 1, field_x_slice_gath, g_nb_x_dn, 91, g_cart_grid,  &requests[ireq]);
    MPI_Irecv((void*)(l+g_1st_x_ext_up), 1, field_x_slice_cont, g_nb_x_up, 91, g_cart_grid, &requests[ireq+1]);
    ireq=ireq+4;
#    endif
    
#    if (defined PARALLELXYT || defined PARALLELXYZT || defined PARALLELXY || defined PARALLELXYZ )
    /* send the data to the neighbour on the left in y direction */
    /* recieve the data from the neighbour on the right in y direction */
    MPI_Isend((void*)(l+g_1st_y_int_dn), 1, field_y_slice_gath, g_nb_y_dn, 101, g_cart_grid, &requests[ireq]);
    MPI_Irecv((void*)(l+g_1st_y_ext_up), 1, field_y_slice_cont, g_nb_y_up, 101, g_cart_grid, &requests[ireq+1]);
    ireq=ireq+4;
#    endif

#    if (defined PARALLELXYZT || defined PARALLELXYZ )
    /* This is now depending on whether the field is even or odd */
    /* send the data to the neighbour on the left in z direction */
    /* recieve the data from the neighbour on the right in z direction */


    if(ieo == 1) { 
      MPI_Isend((void*)(l+g_1st_z_int_dn),1,field_z_slice_even_dn,g_nb_z_dn,503,g_cart_grid,&requests[ireq]);
      MPI_Irecv((void*)(l+g_1st_z_ext_up),1,field_z_slice_cont,g_nb_z_up,503,g_cart_grid,&requests[ireq+1]); 
    } else { 
      MPI_Isend((void*)(l+g_1st_z_int_dn),1,field_z_slice_odd_dn,g_nb_z_dn,503,g_cart_grid,&requests[ireq]);
      MPI_Irecv((void*)(l+g_1st_z_ext_up),1,field_z_slice_cont,g_nb_z_up,503,g_cart_grid,&requests[ireq+1]); 
    } 

#    endif


    ireq=2;

#    if (defined PARALLELT || defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT )
    /* send the data to the neighbour on the right */
    /* recieve the data from the neighbour on the left */
    MPI_Isend((void*)(l+g_1st_t_int_up), 1, field_time_slice_cont, g_nb_t_up, 82, g_cart_grid, &requests[ireq]);
    MPI_Irecv((void*)(l+g_1st_t_ext_dn), 1, field_time_slice_cont, g_nb_t_dn, 82, g_cart_grid, &requests[ireq+1]);
    ireq=ireq+4;
#    endif
    
#    if (defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT || defined PARALLELX || defined PARALLELXY || defined PARALLELXYZ )
    /* send the data to the neighbour on the right in x direction */
    /* recieve the data from the neighbour on the left in x direction */  
    MPI_Isend((void*)(l+g_1st_x_int_up), 1, field_x_slice_gath, g_nb_x_up, 92, g_cart_grid, &requests[ireq]);
    MPI_Irecv((void*)(l+g_1st_x_ext_dn), 1, field_x_slice_cont, g_nb_x_dn, 92, g_cart_grid, &requests[ireq+1]);
    ireq=ireq+4;
#    endif
    
#    if (defined PARALLELXYT || defined PARALLELXYZT || defined PARALLELXY || defined PARALLELXYZ )
    /* send the data to the neighbour on the right in y direction */
    /* recieve the data from the neighbour on the left in y direction */  
    MPI_Isend((void*)(l+g_1st_y_int_up), 1, field_y_slice_gath, g_nb_y_up, 102, g_cart_grid, &requests[ireq]);
    MPI_Irecv((void*)(l+g_1st_y_ext_dn), 1, field_y_slice_cont, g_nb_y_dn, 102, g_cart_grid, &requests[ireq+1]);
    ireq=ireq+4;
#    endif
    
#    if ( defined PARALLELXYZT || defined PARALLELXYZ )
  /* send the data to the neighbour on the right in z direction */
  /* recieve the data from the neighbour on the left in z direction */  
    if(ieo == 1) { 
      MPI_Isend((void*)(l+g_1st_z_int_up),1,field_z_slice_even_up,g_nb_z_up,504,g_cart_grid,&requests[ireq]);
      MPI_Irecv((void*)(l+g_1st_z_ext_dn),1,field_z_slice_cont,g_nb_z_dn,504,g_cart_grid,&requests[ireq+1]); 
    } else {  
      MPI_Isend((void*)(l+g_1st_z_int_up),1,field_z_slice_odd_up,g_nb_z_up,504,g_cart_grid,&requests[ireq]);
      MPI_Irecv((void*)(l+g_1st_z_ext_dn),1,field_z_slice_cont,g_nb_z_dn,504,g_cart_grid,&requests[ireq+1]); 
    } 
#    endif

  } else {
    ireq=0;
    
#    if (defined PARALLELT || defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT )
    /* send the data to the neighbour on the right */
    /* recieve the data from the neighbour on the left */
    MPI_Isend((void*)(l+g_1st_t_int_up), 1, field_time_slice_cont, g_nb_t_up, 82, g_cart_grid, &requests[ireq]);
    MPI_Irecv((void*)(l+g_1st_t_ext_dn), 1, field_time_slice_cont, g_nb_t_dn, 82, g_cart_grid, &requests[ireq+1]);
    ireq=ireq+4;
#    endif

#    if (defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT || defined PARALLELX || defined PARALLELXY || defined PARALLELXYZ )
    /* send the data to the neighbour on the right in x direction */
    /* recieve the data from the neighbour on the left in x direction */  
    MPI_Isend((void*)(l+g_1st_x_int_up), 1, field_x_slice_gath, g_nb_x_up, 92, g_cart_grid, &requests[ireq]);
    MPI_Irecv((void*)(l+g_1st_x_ext_dn), 1, field_x_slice_cont, g_nb_x_dn, 92, g_cart_grid, &requests[ireq+1]);
    ireq=ireq+4;
#    endif
    
#    if (defined PARALLELXYT || defined PARALLELXYZT || defined PARALLELXY || defined PARALLELXYZ )
    /* send the data to the neighbour on the right in y direction */
    /* recieve the data from the neighbour on the left in y direction */  
    MPI_Isend((void*)(l+g_1st_y_int_up), 1, field_y_slice_gath, g_nb_y_up, 102, g_cart_grid, &requests[ireq]);
    MPI_Irecv((void*)(l+g_1st_y_ext_dn), 1, field_y_slice_cont, g_nb_y_dn, 102, g_cart_grid, &requests[ireq+1]);
    ireq=ireq+4;
#    endif

#    if (defined PARALLELXYZT || defined PARALLELXYZ )
    /* This is now depending on whether the field is even or odd */
    /* send the data to the neighbour on the left in z direction */
    /* recieve the data from the neighbour on the right in z direction */
    if(ieo == 1) { 
      MPI_Isend((void*)(l+g_1st_z_int_dn),1,field_z_slice_even_dn,g_nb_z_dn,503,g_cart_grid,&requests[ireq]);
      MPI_Irecv((void*)(l+g_1st_z_ext_up),1,field_z_slice_cont,g_nb_z_up,503,g_cart_grid,&requests[ireq+1]);
    } else { 
      MPI_Isend((void*)(l+g_1st_z_int_dn),1,field_z_slice_odd_dn,g_nb_z_dn,503,g_cart_grid,&requests[ireq]);
      MPI_Irecv((void*)(l+g_1st_z_ext_up),1,field_z_slice_cont,g_nb_z_up,503,g_cart_grid,&requests[ireq+1]);
    } 
#    endif    

    ireq=2;

#    if (defined PARALLELT || defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT )
    /* send the data to the neighbour on the left */
    /* recieve the data from the neighbour on the right */
    MPI_Isend((void*)(l+g_1st_t_int_dn), 1, field_time_slice_cont, g_nb_t_dn, 81, g_cart_grid, &requests[ireq]);
    MPI_Irecv((void*)(l+g_1st_t_ext_up), 1, field_time_slice_cont, g_nb_t_up, 81, g_cart_grid, &requests[ireq+1]);
    ireq=ireq+4;
#    endif

#    if (defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT || defined PARALLELX || defined PARALLELXY || defined PARALLELXYZ )
    /* send the data to the neighbour on the left in x direction */
    /* recieve the data from the neighbour on the right in x direction */
    MPI_Isend((void*)(l+g_1st_x_int_dn), 1, field_x_slice_gath, g_nb_x_dn, 91, g_cart_grid, &requests[ireq]);
    MPI_Irecv((void*)(l+g_1st_x_ext_up), 1, field_x_slice_cont, g_nb_x_up, 91, g_cart_grid, &requests[ireq+1]);
    ireq=ireq+4;
#    endif
    
#    if (defined PARALLELXYT || defined PARALLELXYZT || defined PARALLELXY || defined PARALLELXYZ )
    /* send the data to the neighbour on the left in y direction */
    /* recieve the data from the neighbour on the right in y direction */
    MPI_Isend((void*)(l+g_1st_y_int_dn), 1, field_y_slice_gath, g_nb_y_dn, 101, g_cart_grid, &requests[ireq]);
    MPI_Irecv((void*)(l+g_1st_y_ext_up), 1, field_y_slice_cont, g_nb_y_up, 101, g_cart_grid, &requests[ireq+1]);
    ireq=ireq+4;
#    endif
    
#    if ( defined PARALLELXYZT  || defined PARALLELXYZ )
    /* send the data to the neighbour on the right in z direction */
    /* recieve the data from the neighbour on the left in z direction */  
    if(ieo == 1) { 
      MPI_Isend((void*)(l+g_1st_z_int_up),1,field_z_slice_even_up,g_nb_z_up,504,g_cart_grid,&requests[ireq]);
      MPI_Irecv((void*)(l+g_1st_z_ext_dn),1,field_z_slice_cont,g_nb_z_dn,504,g_cart_grid,&requests[ireq+1]);
    } else { 
      MPI_Isend((void*)(l+g_1st_z_int_up),1,field_z_slice_odd_up,g_nb_z_up,504,g_cart_grid,&requests[ireq]);
      MPI_Irecv((void*)(l+g_1st_z_ext_dn),1,field_z_slice_cont,g_nb_z_dn,504,g_cart_grid,&requests[ireq+1]);
    } 
#    endif

  }
  MPI_Waitall(reqcount, requests, status);


#  endif /* MPI */
  return;
#ifdef _KOJAK_INST
#pragma pomp inst end(xchangefield)
#endif
}


# else /* _INDEX_INDEP_GEOM */

void xchange_field(spinor * const l, const int ieo) {

#ifdef MPI
  MPI_Request requests[16];
  MPI_Status status[16];
#endif
#  ifdef PARALLELT
  int reqcount = 4;
#  elif defined PARALLELXT
  int reqcount = 8;
#  elif defined PARALLELXYT
  int reqcount = 12;
#  elif defined PARALLELXYZT
  int ix=0;
  int reqcount = 16;
#  endif

#ifdef _KOJAK_INST
#pragma pomp inst begin(xchangefield)
#endif
#  if (defined BGL && defined XLC)
#    ifdef PARALLELXYZT
  __alignx(16, field_buffer_z);
  __alignx(16, field_buffer_z2);
#    endif
  __alignx(16, l);
#  endif

#  ifdef MPI


  /* In 4 dimensions there are two processors sharing the   */
  /* communication bandwidth. So the first should start     */
  /* in forward direction, the second in backward direction */
  /* This might only work if the third direction is         */
  /* parallelised only on the node                          */
  if(g_proc_coords[3]%2 == 0) {
    /* send the data to the neighbour on the left */
    /* recieve the data from the neighbour on the right */
    MPI_Isend((void*)l, 1, field_time_slice_cont, g_nb_t_dn, 81, g_cart_grid, &requests[0]);
    MPI_Irecv((void*)(l+T*LX*LY*LZ/2), 1, field_time_slice_cont, g_nb_t_up, 81, g_cart_grid, &requests[1]);
#    if (defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT)
    /* send the data to the neighbour on the left in x direction */
    /* recieve the data from the neighbour on the right in x direction */
    MPI_Isend((void*)l, 1, field_x_slice_gath, g_nb_x_dn, 91, g_cart_grid,  &requests[4]);
    MPI_Irecv((void*)(l+(T+2)*LX*LY*LZ/2), 1, field_x_slice_cont, g_nb_x_up, 91, g_cart_grid, &requests[5]);
#    endif
    
#    if (defined PARALLELXYT || defined PARALLELXYZT)
    /* send the data to the neighbour on the left in y direction */
    /* recieve the data from the neighbour on the right in y direction */
    MPI_Isend((void*)l, 1, field_y_slice_gath, g_nb_y_dn, 101, g_cart_grid, &requests[8]);
    MPI_Irecv((void*)(l+((T+2)*LX*LY*LZ + 2*T*LY*LZ)/2), 1, field_y_slice_cont, g_nb_y_up, 101, g_cart_grid, &requests[9]);
#    endif
    
#    if (defined PARALLELXYZT)
    /* fill buffer ! */
    /* This is now depending on whether the field is */
    /* even or odd */
    if(ieo == 1) { 
      for(ix = 0; ix < T*LX*LY/2; ix++) { 
 	field_buffer_z[ix] = l[ g_field_z_ipt_even[ix] ];  
      } 
    } 
    else { 
      for(ix = 0; ix < T*LX*LY/2; ix++) { 
 	field_buffer_z[ix] = l[ g_field_z_ipt_odd[ix] ];  
      } 
    } 

    /* send the data to the neighbour on the left in z direction */
    /* recieve the data from the neighbour on the right in z direction */
    MPI_Isend((void*)field_buffer_z, 12*T*LX*LY, MPI_DOUBLE, g_nb_z_dn, 503, g_cart_grid, &requests[12]);
    MPI_Irecv((void*)(l+(VOLUME/2 + LX*LY*LZ + T*LY*LZ +T*LX*LZ)), 12*T*LX*LY, MPI_DOUBLE, g_nb_z_up, 503, g_cart_grid, &requests[13]); 

#    endif
    /* send the data to the neighbour on the right */
    /* recieve the data from the neighbour on the left */
    MPI_Isend((void*)(l+(T-1)*LX*LY*LZ/2), 1, field_time_slice_cont, g_nb_t_up, 82, g_cart_grid, &requests[2]);
    MPI_Irecv((void*)(l+(T+1)*LX*LY*LZ/2), 1, field_time_slice_cont, g_nb_t_dn, 82, g_cart_grid, &requests[3]);
    
#    if (defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT)
    /* send the data to the neighbour on the right in x direction */
    /* recieve the data from the neighbour on the left in x direction */  
    MPI_Isend((void*)(l+(LX-1)*LY*LZ/2), 1, field_x_slice_gath, g_nb_x_up, 92, g_cart_grid, &requests[6]);
    MPI_Irecv((void*)(l+((T+2)*LX*LY*LZ + T*LY*LZ)/2), 1, field_x_slice_cont, g_nb_x_dn, 92, g_cart_grid, &requests[7]);
#    endif
    
#    if (defined PARALLELXYT || defined PARALLELXYZT)
    /* send the data to the neighbour on the right in y direction */
    /* recieve the data from the neighbour on the left in y direction */  
    MPI_Isend((void*)(l+(LY-1)*LZ/2), 1, field_y_slice_gath, g_nb_y_up, 102, g_cart_grid, &requests[10]);
    MPI_Irecv((void*)(l+((T+2)*LX*LY*LZ + 2*T*LY*LZ + T*LX*LZ)/2), 1, field_y_slice_cont, g_nb_y_dn, 102, g_cart_grid, &requests[11]);
#    endif
    
#    if defined PARALLELXYZT
    if(ieo == 1) { 
      for(ix = T*LX*LY/2; ix < T*LX*LY; ix++) { 
	field_buffer_z2[ix-T*LX*LY/2] = l[ g_field_z_ipt_even[ix] ];  
      } 
    } 
    else { 
      for(ix = T*LX*LY/2; ix < T*LX*LY; ix++) { 
 	field_buffer_z2[ix-T*LX*LY/2] = l[ g_field_z_ipt_odd[ix] ];  
      } 
    } 
  /* send the data to the neighbour on the right in y direction */
  /* recieve the data from the neighbour on the left in y direction */  
    MPI_Isend((void*)field_buffer_z2, 12*T*LX*LY, MPI_DOUBLE, g_nb_z_up, 504, g_cart_grid, &requests[14]);
    MPI_Irecv((void*)(l+(VOLUME + 2*LX*LY*LZ + 2*T*LY*LZ + 2*T*LX*LZ + T*LX*LY)/2), 12*T*LX*LY, MPI_DOUBLE, g_nb_z_dn, 504, g_cart_grid, &requests[15]); 
#    endif

  }
  else {
    /* send the data to the neighbour on the right */
    /* recieve the data from the neighbour on the left */
    MPI_Isend((void*)(l+(T-1)*LX*LY*LZ/2), 1, field_time_slice_cont, g_nb_t_up, 82, g_cart_grid, &requests[0]);
    MPI_Irecv((void*)(l+(T+1)*LX*LY*LZ/2), 1, field_time_slice_cont, g_nb_t_dn, 82, g_cart_grid, &requests[1]);
#    if (defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT)
    /* send the data to the neighbour on the right in x direction */
    /* recieve the data from the neighbour on the left in x direction */  
    MPI_Isend((void*)(l+(LX-1)*LY*LZ/2), 1, field_x_slice_gath, g_nb_x_up, 92, g_cart_grid, &requests[4]);
    MPI_Irecv((void*)(l+((T+2)*LX*LY*LZ + T*LY*LZ)/2), 1, field_x_slice_cont, g_nb_x_dn, 92, g_cart_grid, &requests[5]);
#    endif
    
#    if (defined PARALLELXYT || defined PARALLELXYZT)
    /* send the data to the neighbour on the right in y direction */
    /* recieve the data from the neighbour on the left in y direction */  
    MPI_Isend((void*)(l+(LY-1)*LZ/2), 1, field_y_slice_gath, g_nb_y_up, 102, g_cart_grid, &requests[8]);
    MPI_Irecv((void*)(l+((T+2)*LX*LY*LZ + 2*T*LY*LZ + T*LX*LZ)/2), 1, field_y_slice_cont, g_nb_y_dn, 102, g_cart_grid, &requests[9]);
#    endif

#    if (defined PARALLELXYZT)
    /* fill buffer ! */
    /* This is now depending on whether the field is */
    /* even or odd */
    if(ieo == 1) { 
      for(ix = 0; ix < T*LX*LY/2; ix++) { 
 	field_buffer_z[ix] = l[ g_field_z_ipt_even[ix] ];  
      } 
    } 
    else { 
      for(ix = 0; ix < T*LX*LY/2; ix++) { 
 	field_buffer_z[ix] = l[ g_field_z_ipt_odd[ix] ];  
      } 
    } 
    /* send the data to the neighbour on the left in z direction */
    /* recieve the data from the neighbour on the right in z direction */
    MPI_Isend((void*)field_buffer_z, 12*T*LX*LY, MPI_DOUBLE, g_nb_z_dn, 503, g_cart_grid, &requests[12]);
    MPI_Irecv((void*)(l+(VOLUME/2 + LX*LY*LZ + T*LY*LZ +T*LX*LZ)), 12*T*LX*LY, MPI_DOUBLE, g_nb_z_up, 503, g_cart_grid, &requests[13]); 
#    endif    

    /* send the data to the neighbour on the left */
    /* recieve the data from the neighbour on the right */
    MPI_Isend((void*)l, 1, field_time_slice_cont, g_nb_t_dn, 81, g_cart_grid, &requests[2]);
    MPI_Irecv((void*)(l+T*LX*LY*LZ/2), 1, field_time_slice_cont, g_nb_t_up, 81, g_cart_grid, &requests[3]);
#    if (defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT)
    /* send the data to the neighbour on the left in x direction */
    /* recieve the data from the neighbour on the right in x direction */
    MPI_Isend((void*)l, 1, field_x_slice_gath, g_nb_x_dn, 91, g_cart_grid,  &requests[6]);
    MPI_Irecv((void*)(l+(T+2)*LX*LY*LZ/2), 1, field_x_slice_cont, g_nb_x_up, 91, g_cart_grid, &requests[7]);
#    endif
    
#    if (defined PARALLELXYT || defined PARALLELXYZT)
    /* send the data to the neighbour on the left in y direction */
    /* recieve the data from the neighbour on the right in y direction */
    MPI_Isend((void*)l, 1, field_y_slice_gath, g_nb_y_dn, 101, g_cart_grid, &requests[10]);
    MPI_Irecv((void*)(l+((T+2)*LX*LY*LZ + 2*T*LY*LZ)/2), 1, field_y_slice_cont, g_nb_y_up, 101, g_cart_grid, &requests[11]);
#    endif
    
#    if defined PARALLELXYZT
    if(ieo == 1) { 
      for(ix = T*LX*LY/2; ix < T*LX*LY; ix++) { 
 	field_buffer_z2[ix-T*LX*LY/2] = l[ g_field_z_ipt_even[ix] ];  
      } 
    } 
    else { 
      for(ix = T*LX*LY/2; ix < T*LX*LY; ix++) { 
 	field_buffer_z2[ix-T*LX*LY/2] = l[ g_field_z_ipt_odd[ix] ];  
      } 
    } 
    /* send the data to the neighbour on the right in y direction */
    /* recieve the data from the neighbour on the left in y direction */  
    MPI_Isend((void*)field_buffer_z2, 12*T*LX*LY, MPI_DOUBLE, g_nb_z_up, 504, g_cart_grid, &requests[14]);
    MPI_Irecv((void*)(l+(VOLUME + 2*LX*LY*LZ + 2*T*LY*LZ + 2*T*LX*LZ + T*LX*LY)/2), 12*T*LX*LY, MPI_DOUBLE, g_nb_z_dn, 504, g_cart_grid, &requests[15]); 
#    endif

  }
  MPI_Waitall(reqcount, requests, status);

#  endif


  return;
#ifdef _KOJAK_INST
#pragma pomp inst end(xchangefield)
#endif
}

# endif /* _INDEX_INDEP_GEOM */

#elif (defined _USE_SHMEM) /* _NON_BLOCKING */

/* Here comes the version with shared memory */
/* exchanges the field  l */
void xchange_field(spinor * const l, const int ieo) {

#  ifdef MPI
  int i,ix, mu, x0, x1, x2, x3, k;

#ifdef _KOJAK_INST
#pragma pomp inst begin(xchangefield)
#endif

  shmem_barrier_all();

  shmem_double_put((double*)(l+T*LX*LY*LZ/2), (double*)l,
                   (LX*LY*LZ*12), g_nb_t_dn);
  shmem_double_put((double*)(l+(T+1)*LX*LY*LZ/2), (double*)(l+(T-1)*LX*LY*LZ/2),
                   (LX*LY*LZ*12), g_nb_t_up);

#    if (defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT)
  k = (T+2)*LX*LY*LZ/2;
  for(x0 = 0; x0 < T; x0++) {
    shmem_double_put((double*)(l + k),
                     (double*)(l + g_lexic2eo[g_ipt[x0][0][0][0]]),
                     12*LZ*LY, g_nb_x_dn);
    k+=LZ*LY;
  }
  k = ((T+2)*LX*LY*LZ + T*LY*LZ)/2;
  for(x0 = 0; x0 < T; x0++) {
    shmem_double_put((double*)(l + k),
                     (double*)(l + g_lexic2eo[g_ipt[x0][LX-1][0][0]]),
                     12*LZ*LY, g_nb_x_up);
    k+=LZ*LY;
  }
#    endif

#    if (defined PARALLELXYT || defined PARALLELXYZT)
  k = ((T+2)*LX*LY*LZ + 2*T*LY*LZ)/2;
  for(x0 = 0; x0 < T; x0++) {
    for(x1 = 0; x1 < LX; x1++) {
      shmem_double_put((double*)(l + k),
                       (double*)(l + g_lexic2eo[g_ipt[x0][x1][0][0]]),
                       12*LZ, g_nb_y_dn);
      k+=LZ;
    }
  }
  k = ((T+2)*LX*LY*LZ + 2*T*LY*LZ + T*LX*LZ)/2;
  for(x0 = 0; x0 < T; x0++) {
    for(x1 = 0; x1 < LX; x1++) {
      shmem_double_put((double*)(l + k),
                       (double*)(l + g_lexic2eo[g_ipt[x0][x1][LY-1][0]]),
                       12*LZ, g_nb_y_up);
      k+=LZ;
    }
  }
#    endif

#    if (defined PARALLELXYZT)
  x0 = (VOLUME/2 + LX*LY*LZ + T*LY*LZ +T*LX*LZ);
  if(ieo == 1) {
    for(k = 0; k < T*LX*LY/2; k++) {
      shmem_double_put((double*)(l + x0),
                       (double*)(l + g_field_z_ipt_even[k]),
                       24, g_nb_z_dn);
      x0++;
    }
  }
  else {
    for(k = 0; k < T*LX*LY/2; k++) {
      shmem_double_put((double*)(l + x0),
                       (double*)(l + g_field_z_ipt_odd[k]),
                       24, g_nb_z_dn);
      x0++;
    }
  }
  x0 = (VOLUME/2 + LX*LY*LZ + T*LY*LZ + T*LX*LZ + T*LX*LY/2);
  if(ieo == 1) {
    for(k = T*LX*LY/2; k < T*LX*LY; k++) {
      shmem_double_put((double*)(l + x0),
                       (double*)(l + g_field_z_ipt_even[k]),
                       24, g_nb_z_up);
      x0++;
    }
  }
  else {
    for(k = T*LX*LY/2; k < T*LX*LY; k++) {
      shmem_double_put((double*)(l + x0),
                       (double*)(l + g_field_z_ipt_even[k]),
                       24, g_nb_z_up);
      x0++;
    }
  }
#    endif

  shmem_barrier_all();
#  endif // MPI
  return;
#ifdef _KOJAK_INST
#pragma pomp inst end(xchangefield)
#endif
}


/* Here comes the naive version */  
/* Using MPI_Sendrecv */
#else  /* _NON_BLOCKING _USE_SHMEM */


/* this is the version independent of the content of the function Index  */
# ifdef  _INDEX_INDEP_GEOM

/* exchanges the field  l */
void xchange_field(spinor * const l, const int ieo) {
  
#  ifdef PARALLELXYZT
  int x0=0, x1=0, x2=0, ix=0;
#  endif
#ifdef _KOJAK_INST
#pragma pomp inst begin(xchangefield)
#endif

#  ifdef MPI
    
  MPI_Status status;

#    if (defined PARALLELT || defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT )
  /* send the data to the neighbour on the left */
  /* recieve the data from the neighbour on the right */
  MPI_Sendrecv((void*)(l+g_1st_t_int_dn),  1, field_time_slice_cont, g_nb_t_dn, 81,
	       (void*)(l+g_1st_t_ext_up), 1, field_time_slice_cont, g_nb_t_up, 81,
	       g_cart_grid, &status);
    
  /* send the data to the neighbour on the right */
  /* recieve the data from the neighbour on the left */
  MPI_Sendrecv((void*)(l+g_1st_t_int_up), 1, field_time_slice_cont, g_nb_t_up, 82,
	       (void*)(l+g_1st_t_ext_dn), 1, field_time_slice_cont, g_nb_t_dn, 82,
	       g_cart_grid, &status);
#    endif    

#    if (defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT || defined PARALLELX || defined PARALLELXY || defined PARALLELXYZ )
  /* send the data to the neighbour on the left in x direction */
  /* recieve the data from the neighbour on the right in x direction */
  MPI_Sendrecv((void*)(l+g_1st_x_int_dn), 1, field_x_slice_gath, g_nb_x_dn, 91, 
	       (void*)(l+g_1st_x_ext_up), 1, field_x_slice_cont, g_nb_x_up, 91,
	       g_cart_grid, &status);
    
  /* send the data to the neighbour on the right in x direction */
  /* recieve the data from the neighbour on the left in x direction */  
  MPI_Sendrecv((void*)(l+g_1st_x_int_up), 1, field_x_slice_gath, g_nb_x_up, 92, 
	       (void*)(l+g_1st_x_ext_dn), 1, field_x_slice_cont, g_nb_x_dn, 92,
	       g_cart_grid, &status);
    
#    endif
    
#    if (defined PARALLELXYT || defined PARALLELXYZT || defined PARALLELXY || defined PARALLELXYZ )
  /* send the data to the neighbour on the left in y direction */
  /* recieve the data from the neighbour on the right in y direction */
  MPI_Sendrecv((void*)(l+g_1st_y_int_dn), 1, field_y_slice_gath, g_nb_y_dn, 101, 
	       (void*)(l+g_1st_y_ext_up), 1, field_y_slice_cont, g_nb_y_up, 101,
	       g_cart_grid, &status);
    
  /* send the data to the neighbour on the right in y direction */
  /* recieve the data from the neighbour on the left in y direction */  
  MPI_Sendrecv((void*)(l+g_1st_y_int_up), 1, field_y_slice_gath, g_nb_y_up, 102, 
	       (void*)(l+g_1st_y_ext_dn), 1, field_y_slice_cont, g_nb_y_dn, 102,
	       g_cart_grid, &status);
    
#    endif
    
#    if (defined PARALLELXYZT || defined PARALLELXYZ )
  /* send the data to the neighbour on the left in z direction */
  /* recieve the data from the neighbour on the right in z direction */
  if(ieo == 1){
    MPI_Sendrecv((void*)(l+g_1st_z_int_dn),1,field_z_slice_even_dn, g_nb_z_dn, 503,  
		 (void*)(l+g_1st_z_ext_up),1,field_z_slice_cont, g_nb_z_up, 503, 
		 g_cart_grid, &status); 
  } else {
    MPI_Sendrecv((void*)(l+g_1st_z_int_dn),1,field_z_slice_odd_dn, g_nb_z_dn, 503,  
		 (void*)(l+g_1st_z_ext_up),1,field_z_slice_cont, g_nb_z_up, 503, 
		 g_cart_grid, &status); 
  }
    
  /* send the data to the neighbour on the right in z direction */
  /* recieve the data from the neighbour on the left in z direction */  
  if(ieo == 1){
    MPI_Sendrecv((void*)(l+g_1st_z_int_up),1,field_z_slice_even_up, g_nb_z_up, 504, 
		 (void*)(l+g_1st_z_ext_dn),1,field_z_slice_cont, g_nb_z_dn, 504, 
		 g_cart_grid, &status); 
  } else {
    MPI_Sendrecv((void*)(l+g_1st_z_int_up),1,field_z_slice_odd_up, g_nb_z_up, 504, 
		 (void*)(l+g_1st_z_ext_dn),1,field_z_slice_cont, g_nb_z_dn, 504, 
		 g_cart_grid, &status); 
  }

#    endif
#  endif // MPI
  return;
#ifdef _KOJAK_INST
#pragma pomp inst end(xchangefield)
#endif
}


# else  /* _INDEX_INDEP_GEOM */

/* exchanges the field  l */
void xchange_field(spinor * const l, const int ieo) {
  
#  ifdef PARALLELXYZT
  int x0=0, x1=0, x2=0, ix=0;
#  endif
#ifdef _KOJAK_INST
#pragma pomp inst begin(xchangefield)
#endif

#  ifdef MPI
    
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
    
#    if (defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT)
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
    
#    endif
    
#    if (defined PARALLELXYT || defined PARALLELXYZT)
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
    
#    endif
    
#    if (defined PARALLELXYZT)
  /* fill buffer ! */
  /* This is now depending on whether the field is */
  /* even or odd */
  if(ieo == 1) {
    for(ix = 0; ix < T*LX*LY/2; ix++) {
      field_buffer_z[ix] = l[ g_field_z_ipt_even[ix] ]; 
    }
  }
  else {
    for(ix = 0; ix < T*LX*LY/2; ix++) {
      field_buffer_z[ix] = l[ g_field_z_ipt_odd[ix] ]; 
    }
  }
  /* send the data to the neighbour on the left in z direction */
  /* recieve the data from the neighbour on the right in z direction */
  MPI_Sendrecv((void*)field_buffer_z, 
	       12*T*LX*LY, MPI_DOUBLE, g_nb_z_dn, 503,  
	       (void*)(l+(VOLUME/2 + LX*LY*LZ + T*LY*LZ +T*LX*LZ)),  
	       12*T*LX*LY, MPI_DOUBLE, g_nb_z_up, 503, 
	       g_cart_grid, &status); 
    
  if(ieo == 1) {
    for(ix = T*LX*LY/2; ix < T*LX*LY; ix++) {
      field_buffer_z[ix-T*LX*LY/2] = l[ g_field_z_ipt_even[ix] ]; 
    }
  }
  else {
    for(ix = T*LX*LY/2; ix < T*LX*LY; ix++) {
      field_buffer_z[ix-T*LX*LY/2] = l[ g_field_z_ipt_odd[ix] ]; 
    }
  }
  /* send the data to the neighbour on the right in y direction */
  /* recieve the data from the neighbour on the left in y direction */  
  MPI_Sendrecv((void*)field_buffer_z,  
	       12*T*LX*LY, MPI_DOUBLE, g_nb_z_up, 504, 
	       (void*)(l+(VOLUME + 2*LX*LY*LZ + 2*T*LY*LZ + 2*T*LX*LZ + T*LX*LY)/2),  
	       12*T*LX*LY, MPI_DOUBLE, g_nb_z_dn, 504, 
	       g_cart_grid, &status); 
    
#    endif
#  endif // MPI
  return;
#ifdef _KOJAK_INST
#pragma pomp inst end(xchangefield)
#endif
}


# endif /* _INDEX_INDEP_GEOM */

#endif /* _NON_BLOCKING */














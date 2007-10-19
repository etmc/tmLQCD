/* $Id$ */

/**********************************************************
 * 
 * exchange routines for 2 spinor fields at once
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
#include "xchange_2fields.h"

#if (defined _NON_BLOCKING)

#if ((defined XLC) && (defined PARALLELXYZT))
#pragma disjoint(*field_buffer_z2, *field_buffer_z, *field_buffer_z3, *field_buffer_z4)
#endif

/* this version uses non-blocking MPI calls */

void xchange_2fields(spinor * const l, spinor * const k, const int ieo) {

  MPI_Request requests[32];
  MPI_Status status[32];
  int reqcount = 0;
#if defined PARALLELXYZT
  int ix=0;
#endif

#ifdef _KOJAK_INST
#pragma pomp inst begin(xchange2fields)
#endif

#  ifdef MPI

#  if (defined BGL && defined XLC)
#    ifdef PARALLELXYZT
  __alignx(16, field_buffer_z);
  __alignx(16, field_buffer_z2);
  __alignx(16, field_buffer_z3);
  __alignx(16, field_buffer_z4);
#    endif
  __alignx(16, l);
#  endif

  /* send the data to the neighbour on the left */
  /* recieve the data from the neighbour on the right */
  MPI_Isend((void*)l, 1, field_time_slice_cont, g_nb_t_dn, 81, g_cart_grid, &requests[reqcount]);
  reqcount++;
  MPI_Irecv((void*)(l+T*LX*LY*LZ/2), 1, field_time_slice_cont, g_nb_t_up, 81, g_cart_grid, &requests[reqcount]);
  reqcount++;
  
  /* send the data to the neighbour on the right */
  /* recieve the data from the neighbour on the left */
  MPI_Isend((void*)(l+(T-1)*LX*LY*LZ/2), 1, field_time_slice_cont, g_nb_t_up, 82, g_cart_grid, &requests[reqcount]);
  reqcount++;
  MPI_Irecv((void*)(l+(T+1)*LX*LY*LZ/2), 1, field_time_slice_cont, g_nb_t_dn, 82, g_cart_grid, &requests[reqcount]);
  reqcount++;
  
  /* send the data to the neighbour on the left */
  /* recieve the data from the neighbour on the right */
  MPI_Isend((void*)k, 1, field_time_slice_cont, g_nb_t_dn, 83, g_cart_grid, &requests[reqcount]);
  reqcount++;
  MPI_Irecv((void*)(k+T*LX*LY*LZ/2), 1, field_time_slice_cont, g_nb_t_up, 83, g_cart_grid, &requests[reqcount]);
  reqcount++;
  
  /* send the data to the neighbour on the right */
  /* recieve the data from the neighbour on the left */
  MPI_Isend((void*)(k+(T-1)*LX*LY*LZ/2), 1, field_time_slice_cont, g_nb_t_up, 84, g_cart_grid, &requests[reqcount]);
  reqcount++;
  MPI_Irecv((void*)(k+(T+1)*LX*LY*LZ/2), 1, field_time_slice_cont, g_nb_t_dn, 84, g_cart_grid, &requests[reqcount]);
  reqcount++;

  
#    if (defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT)
  /* send the data to the neighbour on the left in x direction */
  /* recieve the data from the neighbour on the right in x direction */
  MPI_Isend((void*)l, 1, field_x_slice_gath, g_nb_x_dn, 91, g_cart_grid,  &requests[reqcount]);
  reqcount++;
  MPI_Irecv((void*)(l+(T+2)*LX*LY*LZ/2), 1, field_x_slice_cont, g_nb_x_up, 91, g_cart_grid, &requests[reqcount]);
  reqcount++;
  
    /* send the data to the neighbour on the right in x direction */
    /* recieve the data from the neighbour on the left in x direction */  
  MPI_Isend((void*)(l+(LX-1)*LY*LZ/2), 1, field_x_slice_gath, g_nb_x_up, 92, g_cart_grid, &requests[reqcount]);
  reqcount++;
  MPI_Irecv((void*)(l+((T+2)*LX*LY*LZ + T*LY*LZ)/2), 1, field_x_slice_cont, g_nb_x_dn, 92, g_cart_grid, &requests[reqcount]);
  reqcount++;
  
  /* send the data to the neighbour on the left in x direction */
  /* recieve the data from the neighbour on the right in x direction */
  MPI_Isend((void*)k, 1, field_x_slice_gath, g_nb_x_dn, 93, g_cart_grid,  &requests[reqcount]);
  reqcount++;
  MPI_Irecv((void*)(k+(T+2)*LX*LY*LZ/2), 1, field_x_slice_cont, g_nb_x_up, 93, g_cart_grid, &requests[reqcount]);
  reqcount++;
  
  /* send the data to the neighbour on the right in x direction */
  /* recieve the data from the neighbour on the left in x direction */  
  MPI_Isend((void*)(k+(LX-1)*LY*LZ/2), 1, field_x_slice_gath, g_nb_x_up, 94, g_cart_grid, &requests[reqcount]);
  reqcount++;
  MPI_Irecv((void*)(k+((T+2)*LX*LY*LZ + T*LY*LZ)/2), 1, field_x_slice_cont, g_nb_x_dn, 94, g_cart_grid, &requests[reqcount]);
  reqcount++;
#    endif
  
#    if (defined PARALLELXYT || defined PARALLELXYZT)
  /* send the data to the neighbour on the left in y direction */
  /* recieve the data from the neighbour on the right in y direction */
  MPI_Isend((void*)l, 1, field_y_slice_gath, g_nb_y_dn, 101, g_cart_grid, &requests[reqcount]);
  reqcount++;
  MPI_Irecv((void*)(l+((T+2)*LX*LY*LZ + 2*T*LY*LZ)/2), 1, field_y_slice_cont, g_nb_y_up, 101, g_cart_grid, &requests[reqcount]);
  reqcount++;
  
  /* send the data to the neighbour on the right in y direction */
  /* recieve the data from the neighbour on the left in y direction */  
  MPI_Isend((void*)(l+(LY-1)*LZ/2), 1, field_y_slice_gath, g_nb_y_up, 102, g_cart_grid, &requests[reqcount]);
  reqcount++;
  MPI_Irecv((void*)(l+((T+2)*LX*LY*LZ + 2*T*LY*LZ + T*LX*LZ)/2), 1, field_y_slice_cont, g_nb_y_dn, 102, g_cart_grid, &requests[reqcount]);
  reqcount++;
  
  /* send the data to the neighbour on the left in y direction */
  /* recieve the data from the neighbour on the right in y direction */
  MPI_Isend((void*)k, 1, field_y_slice_gath, g_nb_y_dn, 103, g_cart_grid, &requests[reqcount]);
  reqcount++;
  MPI_Irecv((void*)(k+((T+2)*LX*LY*LZ + 2*T*LY*LZ)/2), 1, field_y_slice_cont, g_nb_y_up, 103, g_cart_grid, &requests[reqcount]);
  reqcount++;
  
  /* send the data to the neighbour on the right in y direction */
  /* recieve the data from the neighbour on the left in y direction */  
  MPI_Isend((void*)(k+(LY-1)*LZ/2), 1, field_y_slice_gath, g_nb_y_up, 104, g_cart_grid, &requests[reqcount]);
  reqcount++;
  MPI_Irecv((void*)(k+((T+2)*LX*LY*LZ + 2*T*LY*LZ + T*LX*LZ)/2), 1, field_y_slice_cont, g_nb_y_dn, 104, g_cart_grid, &requests[reqcount]);
  reqcount++;
  
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
  /* send the data to the neighbour on the left in z direction */
  /* recieve the data from the neighbour on the right in z direction */
  MPI_Isend((void*)field_buffer_z, 12*T*LX*LY, MPI_DOUBLE, g_nb_z_dn, 503, g_cart_grid, &requests[reqcount]); 
  reqcount++;
  MPI_Irecv((void*)(l+(VOLUME/2 + LX*LY*LZ + T*LY*LZ +T*LX*LZ)), 12*T*LX*LY, MPI_DOUBLE, g_nb_z_up, 503, g_cart_grid, &requests[reqcount]); 
  reqcount++;
  
  /* send the data to the neighbour on the right in y direction */
  /* recieve the data from the neighbour on the left in y direction */  
  MPI_Isend((void*)field_buffer_z2, 12*T*LX*LY, MPI_DOUBLE, g_nb_z_up, 504, g_cart_grid, &requests[reqcount]); 
  reqcount++;
  MPI_Irecv((void*)(l+(VOLUME + 2*LX*LY*LZ + 2*T*LY*LZ + 2*T*LX*LZ + T*LX*LY)/2), 12*T*LX*LY, MPI_DOUBLE, g_nb_z_dn, 504, g_cart_grid, &requests[reqcount]); 
  reqcount++;

  /* fill buffer ! */
  /* This is now depending on whether the field is */
  /* even or odd */
  if(ieo == 0) { 
    for(ix = 0; ix < T*LX*LY/2; ix++) { 
      field_buffer_z3[ix] = k[ g_field_z_ipt_even[ix] ];  
    } 
  } 
  else { 
    for(ix = 0; ix < T*LX*LY/2; ix++) { 
      field_buffer_z3[ix] = k[ g_field_z_ipt_odd[ix] ];  
    } 
  } 
  if(ieo == 0) { 
    for(ix = T*LX*LY/2; ix < T*LX*LY; ix++) { 
      field_buffer_z4[ix-T*LX*LY/2] = k[ g_field_z_ipt_even[ix] ];  
    } 
  } 
  else { 
    for(ix = T*LX*LY/2; ix < T*LX*LY; ix++) { 
      field_buffer_z4[ix-T*LX*LY/2] = k[ g_field_z_ipt_odd[ix] ];  
    } 
  } 
  /* send the data to the neighbour on the left in z direction */
  /* recieve the data from the neighbour on the right in z direction */
  MPI_Isend((void*)field_buffer_z3, 12*T*LX*LY, MPI_DOUBLE, g_nb_z_dn, 505, g_cart_grid, &requests[reqcount]); 
  reqcount++;
  MPI_Irecv((void*)(k+(VOLUME/2 + LX*LY*LZ + T*LY*LZ +T*LX*LZ)), 12*T*LX*LY, MPI_DOUBLE, g_nb_z_up, 505, g_cart_grid, &requests[reqcount]); 
  reqcount++;
  
  /* send the data to the neighbour on the right in y direction */
  /* recieve the data from the neighbour on the left in y direction */  
  MPI_Isend((void*)field_buffer_z4, 12*T*LX*LY, MPI_DOUBLE, g_nb_z_up, 506, g_cart_grid, &requests[reqcount]); 
  reqcount++;
  MPI_Irecv((void*)(k+(VOLUME + 2*LX*LY*LZ + 2*T*LY*LZ + 2*T*LX*LZ + T*LX*LY)/2), 12*T*LX*LY, MPI_DOUBLE, g_nb_z_dn, 506, g_cart_grid, &requests[reqcount]); 
  reqcount++;

  
#    endif


  MPI_Waitall(reqcount, requests, status);
#  endif
  return;
#ifdef _KOJAK_INST
#pragma pomp inst end(xchange2fields)
#endif
}

#  endif

static char const rcsid[] = "$Id$";


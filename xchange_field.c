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
#ifdef BGL
#  include "bgl.h"
#endif
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
/*   double _Complex reg00, reg01, reg02, reg03, reg04, reg05; */
/*   spinor *s ALIGN; */
/*   spinor *sp ALIGN; */
/*   spinor *r ALIGN; */
  __alignx(16, l);
  __alignx(16, field_buffer_z);
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

#  endif

#endif
  return;
}

static char const rcsid[] = "$Id$";

/* Sending only half the amount?? */
#if defined BLA
  if(ieo == 1) {
    sp = &l[ g_field_z_ipt_even[0] ];
    for(ix = 0; ix < T*LX*LY/2-1; ix++) {
/*       field_buffer_z[ix] = l[ g_field_z_ipt_even[ix] ]; */
      s = sp;
      sp = &l[ g_field_z_ipt_even[ix+1] ];
      _prefetch_spinor(sp); 
      _bgl_load((*s).s0);
      _bgl_load_up((*s).s2);
      _bgl_vector_i_mul_add();
      _bgl_store(field_buffer_z[ix].s0);
      
      _bgl_load((*sp).s1);
      _bgl_load_up((*sp).s3);
      _bgl_vector_i_mul_sub();
      _bgl_store(field_buffer_z[ix].s1);
    }
    s = sp;
    _bgl_load((*s).s0);
    _bgl_load_up((*s).s2);
    _bgl_vector_i_mul_add();
    _bgl_store(field_buffer_z[ix].s0);
    
    _bgl_load((*sp).s1);
    _bgl_load_up((*sp).s3);
    _bgl_vector_i_mul_sub();
    _bgl_store(field_buffer_z[ix].s1);
  }
  else {
    sp = &l[ g_field_z_ipt_odd[0] ];
    for(ix = 0; ix < T*LX*LY/2-1; ix++) {
      /*       field_buffer_z[ix] = l[ g_field_z_ipt_odd[ix] ]; */
      s = sp;
      sp = &l[ g_field_z_ipt_odd[ix+1] ];
      _prefetch_spinor(sp); 
      _bgl_load((*s).s0);
      _bgl_load_up((*s).s2);
      _bgl_vector_i_mul_add();
      _bgl_store(field_buffer_z[ix].s0);

      _bgl_load((*sp).s1);
      _bgl_load_up((*sp).s3);
      _bgl_vector_i_mul_sub();
      _bgl_store(field_buffer_z[ix].s1);
    }
    s = sp;
    _bgl_load((*s).s0);
    _bgl_load_up((*s).s2);
    _bgl_vector_i_mul_add();
    _bgl_store(field_buffer_z[ix].s0);
    
    _bgl_load((*sp).s1);
    _bgl_load_up((*sp).s3);
    _bgl_vector_i_mul_sub();
    _bgl_store(field_buffer_z[ix].s1);
  }
  /* send the data to the neighbour on the left in z direction */
  /* recieve the data from the neighbour on the right in z direction */
/*   MPI_Sendrecv((void*)field_buffer_z, */
/* 	       12*T*LX*LY, MPI_DOUBLE, g_nb_z_dn, 503,  */
/*  	       (void*)(l+(VOLUME/2 + LX*LY*LZ + T*LY*LZ +T*LX*LZ)),  */
/* 	       12*T*LX*LY, MPI_DOUBLE, g_nb_z_up, 503, */
/* 	       g_cart_grid, &status); */

  MPI_Sendrecv((void*)field_buffer_z,
 	       1, field_z_slice_half, g_nb_z_dn, 503,  
 	       (void*)(l+(VOLUME/2 + LX*LY*LZ + T*LY*LZ +T*LX*LZ)), 
 	       1, field_z_slice_half, g_nb_z_up, 503, 
	       g_cart_grid, &status);

  if(ieo == 1) {
    sp = &l[ g_field_z_ipt_even[T*LX*LY] ];
    for(ix = T*LX*LY/2; ix < T*LX*LY-1; ix++) {
      /*       field_buffer_z[ix-T*LX*LY/2] = l[ g_field_z_ipt_even[ix] ]; */
      s = sp;
      sp =  &l[ g_field_z_ipt_even[ix+1] ];
      _prefetch_spinor(sp);
      _bgl_load((*s).s0);
      _bgl_load_up((*s).s2);
      _bgl_vector_i_mul_sub();
      _bgl_store(field_buffer_z[ix-T*LX*LY/2].s0);
      
      _bgl_load((*s).s1);
      _bgl_load_up((*s).s3);
      _bgl_vector_i_mul_add();
      _bgl_store(field_buffer_z[ix-T*LX*LY/2].s1);
    }
    s = sp;
    _bgl_load((*s).s0);
    _bgl_load_up((*s).s2);
    _bgl_vector_i_mul_sub();
    _bgl_store(field_buffer_z[ix-T*LX*LY/2].s0);
    
    _bgl_load((*s).s1);
    _bgl_load_up((*s).s3);
    _bgl_vector_i_mul_add();
    _bgl_store(field_buffer_z[ix-T*LX*LY/2].s1);
  }
  else {
    sp = &l[ g_field_z_ipt_odd[T*LX*LY] ];
    for(ix = T*LX*LY/2; ix < T*LX*LY-1; ix++) {
/*       field_buffer_z[ix-T*LX*LY/2] = l[ g_field_z_ipt_odd[ix] ]; */
      s = sp;
      sp =  &l[ g_field_z_ipt_odd[ix+1] ];
      _prefetch_spinor(sp);
      _bgl_load((*s).s0);
      _bgl_load_up((*s).s2);
      _bgl_vector_i_mul_sub();
      _bgl_store(field_buffer_z[ix-T*LX*LY/2].s0);
      
      _bgl_load((*s).s1);
      _bgl_load_up((*s).s3);
      _bgl_vector_i_mul_add();
      _bgl_store(field_buffer_z[ix-T*LX*LY/2].s1);
    }
    s = sp;
    _bgl_load((*s).s0);
    _bgl_load_up((*s).s2);
    _bgl_vector_i_mul_sub();
    _bgl_store(field_buffer_z[ix-T*LX*LY/2].s0);
    
    _bgl_load((*s).s1);
    _bgl_load_up((*s).s3);
    _bgl_vector_i_mul_add();
    _bgl_store(field_buffer_z[ix-T*LX*LY/2].s1);
  }
  /* send the data to the neighbour on the right in y direction */
  /* recieve the data from the neighbour on the left in y direction */  
/*   MPI_Sendrecv((void*)field_buffer_z,  */
/* 	       12*T*LX*LY, MPI_DOUBLE, g_nb_z_up, 504, */
/* 	       (void*)(l+(VOLUME + 2*LX*LY*LZ + 2*T*LY*LZ + 2*T*LX*LZ + T*LX*LY)/2),  */
/* 	       12*T*LX*LY, MPI_DOUBLE, g_nb_z_dn, 504, */
/* 	       g_cart_grid, &status); */

  MPI_Sendrecv((void*)field_buffer_z, 
 	       1, field_z_slice_half, g_nb_z_up, 504,  
	       (void*)(l+(VOLUME + 2*LX*LY*LZ + 2*T*LY*LZ + 2*T*LX*LZ + T*LX*LY)/2), 
 	       1, field_z_slice_half, g_nb_z_dn, 504, 
	       g_cart_grid, &status);
#endif

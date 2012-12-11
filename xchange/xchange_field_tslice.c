/**********************************************************
 * 
 * exchange routines for the borders of a timeslice of spinor fields
 *
 * Author: Luigi Scorzato
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
#include "xchange_field_tslice.h"

#ifdef MPI
# ifdef _USE_TSPLITPAR
void xchange_field_open(spinor * const l, const int ieo, const int x0, MPI_Request * requests, 
			MPI_Status * status) {
  
#ifdef _KOJAK_INST
#pragma pomp inst begin(xchangetslicefield)
#endif
#  if (defined BGL && defined XLC)
  __alignx(16, l); /* ?!? */
#  endif

#  ifdef MPI

#    if (defined PARALLELX || defined PARALLELXY || defined PARALLELXYZ )
  /* send the data to the neighbour on the left in x direction */
  /* recieve the data from the neighbour on the right in x direction */
  MPI_Isend((void*)(l+g_1st_xt_int_dn[x0]), 1, field_xt_slice_int, g_nb_x_dn, 91, g_cart_grid,  &requests[0]);
  MPI_Irecv((void*)(l+g_1st_xt_ext_up[x0]), 1, field_xt_slice_ext, g_nb_x_up, 91, g_cart_grid, &requests[1]);
#    endif
    
#    if (defined PARALLELXY || defined PARALLELXYZ )
  /* send the data to the neighbour on the left in y direction */
  /* recieve the data from the neighbour on the right in y direction */
  MPI_Isend((void*)(l+g_1st_yt_int_dn[x0]), 1, field_yt_slice_int, g_nb_y_dn, 101, g_cart_grid, &requests[4]);
  MPI_Irecv((void*)(l+g_1st_yt_ext_up[x0]), 1, field_yt_slice_ext, g_nb_y_up, 101, g_cart_grid, &requests[5]);
#    endif

#    if (defined PARALLELXYZ)
  /* send the data to the neighbour on the left in z direction */
  /* recieve the data from the neighbour on the right in z direction */
  if(ieo == 1){
   if(x0 % 2 == 0) {
    MPI_Isend((void*)(l+g_1st_zt_int_dn[x0]),1,field_zt_slice_even_dn_et,g_nb_z_dn,111,g_cart_grid,&requests[8]);
    MPI_Irecv((void*)(l+g_1st_zt_ext_up[x0]),1 , field_zt_slice_ext_L, g_nb_z_up, 111, g_cart_grid, &requests[9]);
   } else {
    MPI_Isend((void*)(l+g_1st_zt_int_dn[x0]),1,field_zt_slice_even_dn_ot,g_nb_z_dn,111,g_cart_grid,&requests[8]);
    MPI_Irecv((void*)(l+g_1st_zt_ext_up[x0]),1 , field_zt_slice_ext_S, g_nb_z_up, 111, g_cart_grid, &requests[9]);
   }
  } else {
   if(x0 % 2 == 0) {
    MPI_Isend((void*)(l+g_1st_zt_int_dn[x0]),1,field_zt_slice_odd_dn_et,g_nb_z_dn,111,g_cart_grid,&requests[8]);
    MPI_Irecv((void*)(l+g_1st_zt_ext_up[x0]),1 , field_zt_slice_ext_S, g_nb_z_up, 111, g_cart_grid, &requests[9]);
   } else {
    MPI_Isend((void*)(l+g_1st_zt_int_dn[x0]),1,field_zt_slice_odd_dn_ot,g_nb_z_dn,111,g_cart_grid,&requests[8]);
    MPI_Irecv((void*)(l+g_1st_zt_ext_up[x0]),1 , field_zt_slice_ext_L, g_nb_z_up, 111, g_cart_grid, &requests[9]);
   }
  }
#    endif
    
#    if (defined PARALLELX || defined PARALLELXY || defined PARALLELXYZ )
  /* send the data to the neighbour on the right in x direction */
  /* recieve the data from the neighbour on the left in x direction */  
  MPI_Isend((void*)(l+g_1st_xt_int_up[x0]), 1, field_xt_slice_int, g_nb_x_up, 92, g_cart_grid, &requests[2]);
  MPI_Irecv((void*)(l+g_1st_xt_ext_dn[x0]), 1, field_xt_slice_ext, g_nb_x_dn, 92, g_cart_grid, &requests[3]);
#    endif
    
#    if (defined PARALLELXY || defined PARALLELXYZ )
  /* send the data to the neighbour on the right in y direction */
  /* recieve the data from the neighbour on the left in y direction */  
  MPI_Isend((void*)(l+g_1st_yt_int_up[x0]), 1, field_yt_slice_int, g_nb_y_up, 102, g_cart_grid, &requests[6]);
  MPI_Irecv((void*)(l+g_1st_yt_ext_dn[x0]), 1, field_yt_slice_ext, g_nb_y_dn, 102, g_cart_grid, &requests[7]);
#    endif
    
#    if (defined PARALLELXYZ)
  /* send the data to the neighbour on the right in z direction */
  /* recieve the data from the neighbour on the left in z direction */
  if(ieo == 1){
    if(x0 % 2 == 0) {
      MPI_Isend((void*)(l+g_1st_zt_int_up[x0]),1,field_zt_slice_even_up_et,g_nb_z_up,112,g_cart_grid,&requests[10]);
      MPI_Irecv((void*)(l+g_1st_zt_ext_dn[x0]), 1, field_zt_slice_ext_S, g_nb_z_dn, 112, g_cart_grid, &requests[11]);
    } else {
      MPI_Isend((void*)(l+g_1st_zt_int_up[x0]),1,field_zt_slice_even_up_ot,g_nb_z_up,112,g_cart_grid,&requests[10]);
      MPI_Irecv((void*)(l+g_1st_zt_ext_dn[x0]), 1, field_zt_slice_ext_L, g_nb_z_dn, 112, g_cart_grid, &requests[11]);
    }
  } else {
    if(x0 % 2 == 0) {
      MPI_Isend((void*)(l+g_1st_zt_int_up[x0]),1,field_zt_slice_odd_up_et,g_nb_z_up,112,g_cart_grid,&requests[10]);
      MPI_Irecv((void*)(l+g_1st_zt_ext_dn[x0]), 1, field_zt_slice_ext_L, g_nb_z_dn, 112, g_cart_grid, &requests[11]);
    } else {
      MPI_Isend((void*)(l+g_1st_zt_int_up[x0]),1,field_zt_slice_odd_up_ot,g_nb_z_up,112,g_cart_grid,&requests[10]);
      MPI_Irecv((void*)(l+g_1st_zt_ext_dn[x0]), 1, field_zt_slice_ext_S, g_nb_z_dn, 112, g_cart_grid, &requests[11]);
    }
  }
#    endif

#  endif /* MPI */
  return;
#ifdef _KOJAK_INST
#pragma pomp inst end(xchangetslicefield)
#endif
}


void xchange_field_close(MPI_Request * requests, MPI_Status * status, int reqcount) {

#ifdef _KOJAK_INST
#pragma pomp inst begin(xchangetslicefieldclose)
#endif

  MPI_Waitall(reqcount, requests, status);

#ifdef _KOJAK_INST
#pragma pomp inst end(xchangetslicefieldclose)
#endif

}

void xchange_field_slice(spinor * const l, const int ieo, const int x0) {

#ifdef _KOJAK_INST
#pragma pomp inst begin(xchangetslicefield)
#endif
#  if (defined BGL && defined XLC)
  __alignx(16, l); /* ?!? */
#  endif

#  ifdef MPI

    MPI_Status status;

#    if (defined PARALLELX || defined PARALLELXY || defined PARALLELXYZ )
  /* send the data to the neighbour on the left in x direction */
  /* recieve the data from the neighbour on the right in x direction */
  MPI_Sendrecv((void*)(l+g_1st_xt_int_dn[x0]), 1, field_xt_slice_int, g_nb_x_dn, 91,
	       (void*)(l+g_1st_xt_ext_up[x0]), 1, field_xt_slice_ext, g_nb_x_up, 91,  g_cart_grid, &status);
#    endif
    
#    if (defined PARALLELXY || defined PARALLELXYZ )
  /* send the data to the neighbour on the left in y direction */
  /* recieve the data from the neighbour on the right in y direction */
  MPI_Sendrecv((void*)(l+g_1st_yt_int_dn[x0]), 1, field_yt_slice_int, g_nb_y_dn, 101,
	       (void*)(l+g_1st_yt_ext_up[x0]), 1, field_yt_slice_ext, g_nb_y_up, 101, g_cart_grid, &status);
#    endif

#    if (defined PARALLELXYZ)
  /* send the data to the neighbour on the left in z direction */
  /* recieve the data from the neighbour on the right in z direction */
  if(ieo == 1){
   if(x0 % 2 == 0) {
    MPI_Sendrecv((void*)(l+g_1st_zt_int_dn[x0]),1,field_zt_slice_even_dn_et,g_nb_z_dn,111,
		 (void*)(l+g_1st_zt_ext_up[x0]),1 , field_zt_slice_ext_L, g_nb_z_up, 111, g_cart_grid, &status);
   } else {
    MPI_Sendrecv((void*)(l+g_1st_zt_int_dn[x0]),1,field_zt_slice_even_dn_ot,g_nb_z_dn,111,
		 (void*)(l+g_1st_zt_ext_up[x0]),1 , field_zt_slice_ext_S, g_nb_z_up, 111, g_cart_grid, &status);
   }
  } else {
   if(x0 % 2 == 0) {
    MPI_Sendrecv((void*)(l+g_1st_zt_int_dn[x0]),1,field_zt_slice_odd_dn_et,g_nb_z_dn,111,
		 (void*)(l+g_1st_zt_ext_up[x0]),1 , field_zt_slice_ext_S, g_nb_z_up, 111, g_cart_grid, &status);
   } else {
    MPI_Sendrecv((void*)(l+g_1st_zt_int_dn[x0]),1,field_zt_slice_odd_dn_ot,g_nb_z_dn,111,
		 (void*)(l+g_1st_zt_ext_up[x0]),1 , field_zt_slice_ext_L, g_nb_z_up, 111, g_cart_grid, &status);
   }
  }
#    endif
    
#    if (defined PARALLELX || defined PARALLELXY || defined PARALLELXYZ )
  /* send the data to the neighbour on the right in x direction */
  /* recieve the data from the neighbour on the left in x direction */  
  MPI_Sendrecv((void*)(l+g_1st_xt_int_up[x0]), 1, field_xt_slice_int, g_nb_x_up, 92,
	       (void*)(l+g_1st_xt_ext_dn[x0]), 1, field_xt_slice_ext, g_nb_x_dn, 92, g_cart_grid, &status);
#    endif
    
#    if (defined PARALLELXY || defined PARALLELXYZ )
  /* send the data to the neighbour on the right in y direction */
  /* recieve the data from the neighbour on the left in y direction */  
  MPI_Sendrecv((void*)(l+g_1st_yt_int_up[x0]), 1, field_yt_slice_int, g_nb_y_up, 102,
	       (void*)(l+g_1st_yt_ext_dn[x0]), 1, field_yt_slice_ext, g_nb_y_dn, 102, g_cart_grid, &status);
#    endif
    
#    if (defined PARALLELXYZ)
  /* send the data to the neighbour on the right in z direction */
  /* recieve the data from the neighbour on the left in z direction */
  if(ieo == 1){
    if(x0 % 2 == 0) {
      MPI_Sendrecv((void*)(l+g_1st_zt_int_up[x0]),1,field_zt_slice_even_up_et,g_nb_z_up,112,
		   (void*)(l+g_1st_zt_ext_dn[x0]), 1, field_zt_slice_ext_S, g_nb_z_dn, 112, g_cart_grid, &status);
    } else {
      MPI_Sendrecv((void*)(l+g_1st_zt_int_up[x0]),1,field_zt_slice_even_up_ot,g_nb_z_up,112,
		   (void*)(l+g_1st_zt_ext_dn[x0]), 1, field_zt_slice_ext_L, g_nb_z_dn, 112, g_cart_grid, &status);
    }
  } else {
    if(x0 % 2 == 0) {
      MPI_Sendrecv((void*)(l+g_1st_zt_int_up[x0]),1,field_zt_slice_odd_up_et,g_nb_z_up,112,
		   (void*)(l+g_1st_zt_ext_dn[x0]), 1, field_zt_slice_ext_L, g_nb_z_dn, 112, g_cart_grid, &status);
    } else {
      MPI_Sendrecv((void*)(l+g_1st_zt_int_up[x0]),1,field_zt_slice_odd_up_ot,g_nb_z_up,112,
		   (void*)(l+g_1st_zt_ext_dn[x0]), 1, field_zt_slice_ext_S, g_nb_z_dn, 112, g_cart_grid, &status);
    }
  }
#    endif

#  endif /* MPI */
  return;
#ifdef _KOJAK_INST
#pragma pomp inst end(xchangetslicefield)
#endif
}

# endif // _USE_TSPLITPAR
#endif // MPI









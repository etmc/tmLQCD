/* $Id$ */

/**********************************************************
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
#include "xchange_gauge.h"

#ifdef BGL

void xchange_gauge() {
  MPI_Request request[105];
  MPI_Status status[105];
  int cntr=0;
#  ifdef MPI

  /* send the data to the neighbour on the left */
  /* recieve the data from the neighbour on the right */
  
  MPI_Isend(g_gauge_field[0],          1, gauge_time_slice_cont, g_nb_t_dn, 83,
	    g_cart_grid, &request[cntr]);
  cntr++;
  MPI_Irecv(g_gauge_field[VOLUME], 1, gauge_time_slice_cont, g_nb_t_up, 83, 
	    g_cart_grid, &request[cntr]);
  cntr++;

  /* send the data to the neighbour on the right */
  /* recieve the data from the neighbour on the left */
  MPI_Isend(g_gauge_field[(T-1)*LX*LY*LZ], 1, gauge_time_slice_cont, g_nb_t_up, 84,
	    g_cart_grid, &request[cntr]);
  cntr++;
  MPI_Irecv(g_gauge_field[(T+1)*LX*LY*LZ], 1, gauge_time_slice_cont, g_nb_t_dn, 84, 
	    g_cart_grid, &request[cntr]);
  cntr++;

  if(g_dbw2rand > 0) {
    /* send the data to the neighbour on the left */
    /* recieve the data from the neighbour on the right */
    /* t2-Rand */
    MPI_Isend(g_gauge_field[1*LX*LY*LZ],     1, gauge_time_slice_cont, g_nb_t_dn, 85, 
	      g_cart_grid, &request[cntr]);
    cntr++;
    MPI_Irecv(g_gauge_field[VOLUMEPLUSRAND], 1, gauge_time_slice_cont, g_nb_t_up, 85, 
	      g_cart_grid, &request[cntr]);
    cntr++;

    /* send the data to the neighbour on the right */
    /* recieve the data from the neighbour on the left */
    /* t2-Rand */
    MPI_Isend(g_gauge_field[(T-2)*LX*LY*LZ],          1, gauge_time_slice_cont, g_nb_t_up, 86, 
	      g_cart_grid, &request[cntr]);
    cntr++;
    MPI_Irecv(g_gauge_field[VOLUMEPLUSRAND+LX*LY*LZ], 1, gauge_time_slice_cont, g_nb_t_dn, 86, 
	      g_cart_grid, &request[cntr]);
    cntr++;
  }
  
#    if (defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT)
  /* send the data to the neighbour on the left in x direction */
  /* recieve the data from the neighbour on the right in x direction */
  MPI_Isend(g_gauge_field[0],              1, gauge_x_slice_gath, g_nb_x_dn, 87,
	    g_cart_grid, &request[cntr]);
  cntr++;
  MPI_Irecv(g_gauge_field[(T+2)*LX*LY*LZ], 1, gauge_x_slice_cont, g_nb_x_up, 87, 
	    g_cart_grid, &request[cntr]);
  cntr++;

  /* send the data to the neighbour on the right in x direction */
  /* recieve the data from the neighbour on the left in x direction */
  /* x2-Rand */
  MPI_Isend(g_gauge_field[(LX-1)*LY*LZ],             1, gauge_x_slice_gath, g_nb_x_up, 88,
	    g_cart_grid, &request[cntr]);
  cntr++;
  MPI_Irecv(g_gauge_field[(T+2)*LX*LY*LZ + T*LY*LZ], 1, gauge_x_slice_cont, g_nb_x_dn, 88,
	    g_cart_grid, &request[cntr]);
  cntr++;

  if(g_dbw2rand > 0) {
    /* send the data to the neighbour on the left in x direction */
    /* recieve the data from the neighbour on the right in x direction */
    /* x2-Rand */
    MPI_Isend(g_gauge_field[LY*LZ],                     1, gauge_x_slice_gath, g_nb_x_dn, 89,
	     g_cart_grid, &request[cntr]);
    cntr++;
    MPI_Irecv(g_gauge_field[VOLUMEPLUSRAND+2*LX*LY*LZ], 1, gauge_x_slice_cont, g_nb_x_up, 89, 
	      g_cart_grid, &request[cntr]);
    cntr++;
    /* send the data to the neighbour on the right in x direction */
    /* recieve the data from the neighbour on the left in x direction */
    /* x2-Rand */
    MPI_Isend(g_gauge_field[(LX-2)*LY*LZ],                        1, gauge_x_slice_gath, g_nb_x_up, 90,
	      g_cart_grid, &request[cntr]);
    cntr++;
    MPI_Irecv(g_gauge_field[VOLUMEPLUSRAND+2*LX*LY*LZ + T*LY*LZ], 1, gauge_x_slice_cont, g_nb_x_dn, 90,
	      g_cart_grid, &request[cntr]);
    cntr++;
  }
  MPI_Waitall(cntr, request, status);
  cntr=0;
  /* The edges */

  /* send the data to the neighbour on the left in t direction */
  /* recieve the data from the neighbour on the right in t direction */
  /* is on the x-Rand: xt-edge */
  MPI_Isend(g_gauge_field[(T+2)*LX*LY*LZ], 1, gauge_xt_edge_gath, g_nb_t_dn, 100,
	    g_cart_grid, &request[cntr]);
  cntr++;
  MPI_Irecv(g_gauge_field[VOLUME + RAND],  1, gauge_xt_edge_cont, g_nb_t_up, 100, 
	    g_cart_grid, &request[cntr]);
  cntr++;

  /* send the data to the neighbour on the right in t direction */
  /* recieve the data from the neighbour on the left in t direction */
  /* xt-edge */
  MPI_Isend(g_gauge_field[(T+2)*LX*LY*LZ + (T-1)*LY*LZ], 1, gauge_xt_edge_gath, g_nb_t_up, 101,
	    g_cart_grid, &request[cntr]);
  cntr++;
  MPI_Irecv(g_gauge_field[VOLUME + RAND + 2*LY*LZ],      1, gauge_xt_edge_cont, g_nb_t_dn, 101,
	    g_cart_grid, &request[cntr]);
  cntr++;

  if(g_dbw2rand > 0) {
    /* send the data to the neighbour on the left in t direction */
    /* recieve the data from the neighbour on the right in t direction */
    /* t2x-edge */
    MPI_Isend(g_gauge_field[(T+2)*LX*LY*LZ + LY*LZ],
	      1, gauge_xt_edge_gath, g_nb_t_dn, 102,
	      g_cart_grid, &request[cntr]);
    cntr++;
    MPI_Irecv(g_gauge_field[VOLUMEPLUSRAND + RAND],
	      1, gauge_xt_edge_cont, g_nb_t_up, 102, 
	      g_cart_grid, &request[cntr]);
    cntr++;
    
    /* send the data to the neighbour on the right in t direction */
    /* recieve the data from the neighbour on the left in t direction */
    /* t2x-edge */
    MPI_Isend(g_gauge_field[(T+2)*LX*LY*LZ + (T-2)*LY*LZ],
	      1, gauge_xt_edge_gath, g_nb_t_up, 103,
	      g_cart_grid, &request[cntr]);
    cntr++;
    MPI_Irecv(g_gauge_field[VOLUMEPLUSRAND + RAND + 2*LY*LZ],
	      1, gauge_xt_edge_cont, g_nb_t_dn, 103,
	      g_cart_grid, &request[cntr]);
    cntr++;

    /* send the data to the neighbour on the left in t direction */
    /* recieve the data from the neighbour on the right in t direction */
    /* x2t-edge */
    MPI_Isend(g_gauge_field[VOLUMEPLUSRAND + 2*LX*LY*LZ],
	      1, gauge_xt_edge_gath, g_nb_t_dn, 104,
	      g_cart_grid, &request[cntr]);
    cntr++;
    MPI_Irecv(g_gauge_field[VOLUMEPLUSRAND + RAND + 4*LY*LZ],
	      1, gauge_xt_edge_cont, g_nb_t_up, 104, 
	      g_cart_grid, &request[cntr]);
    cntr++;
    
    /* send the data to the neighbour on the right in t direction */
    /* recieve the data from the neighbour on the left in t direction */
    /* x2t-edge */
    MPI_Isend(g_gauge_field[VOLUMEPLUSRAND + 2*LX*LY*LZ + (T-1)*LY*LZ],
	      1, gauge_xt_edge_gath, g_nb_t_up, 105,
	      g_cart_grid, &request[cntr]);
    cntr++;
    MPI_Irecv(g_gauge_field[VOLUMEPLUSRAND + RAND + 6*LY*LZ],
	      1, gauge_xt_edge_cont, g_nb_t_dn, 105,
	      g_cart_grid, &request[cntr]);
    cntr++;
  }
  /* end of if defined PARALLELXT || PARALLELXYT || PARALLELXYZT*/
#    endif

#    if (defined PARALLELXYT || defined PARALLELXYZT)
  /* send the data to the neighbour on the left in y direction */
  /* recieve the data from the neighbour on the right in y direction */
  MPI_Isend(g_gauge_field[0],                            1, gauge_y_slice_gath, g_nb_y_dn, 106,
	    g_cart_grid, &request[cntr]);
  cntr++;
  MPI_Irecv(g_gauge_field[VOLUME + 2*LZ*(LX*LY + T*LY)], 1, gauge_y_slice_cont, g_nb_y_up, 106, 
	    g_cart_grid, &request[cntr]);
  cntr++;
  /* send the data to the neighbour on the right in y direction */
  /* recieve the data from the neighbour on the left in y direction */
  MPI_Isend(g_gauge_field[(LY-1)*LZ],                              1, gauge_y_slice_gath, g_nb_y_up, 107,
	    g_cart_grid, &request[cntr]);
  cntr++;
  MPI_Irecv(g_gauge_field[VOLUME + 2*LZ*(LX*LY + T*LY) + T*LX*LZ], 1, gauge_y_slice_cont, g_nb_y_dn, 107,
	    g_cart_grid, &request[cntr]);
  cntr++;

  if(g_dbw2rand > 0) {
    /* send the data to the neighbour on the left in x direction */
    /* recieve the data from the neighbour on the right in x direction */
    /* y2-Rand */
    MPI_Isend(g_gauge_field[LZ],                              1, gauge_y_slice_gath, g_nb_y_dn, 108,
	      g_cart_grid, &request[cntr]);
    cntr++;
    MPI_Irecv(g_gauge_field[VOLUMEPLUSRAND+(2*LX+2*T)*LY*LZ], 1, gauge_y_slice_cont, g_nb_y_up, 108, 
	      g_cart_grid, &request[cntr]);
    cntr++;
    /* send the data to the neighbour on the right in x direction */
    /* recieve the data from the neighbour on the left in x direction */
    /* y2-Rand */
    MPI_Isend(g_gauge_field[(LY-2)*LZ],                                 1, gauge_y_slice_gath, g_nb_y_up, 109,
	      g_cart_grid, &request[cntr]);
    cntr++;
    MPI_Irecv(g_gauge_field[VOLUMEPLUSRAND+(2*LX+2*T)*LY*LZ + T*LX*LZ], 1, gauge_y_slice_cont, g_nb_y_dn, 109,
	      g_cart_grid, &request[cntr]);
    cntr++;
  }

  MPI_Waitall(cntr, request, status);
  cntr=0;
  /* jetzt wirds richtig eklig ... */

  /* edges */

  /* send the data to the neighbour on the left in x direction */
  /* recieve the data from the neighbour on the right in x direction */
  /* is on the y-Rand -> yx-edge*/
  MPI_Isend(g_gauge_field[VOLUME + 2*LX*LY*LZ + 2*T*LY*LZ], 1, gauge_yx_edge_gath, g_nb_x_dn, 110,
	    g_cart_grid, &request[cntr]);
  cntr++;
  MPI_Irecv(g_gauge_field[VOLUME + RAND + 4*LY*LZ],      1, gauge_yx_edge_cont, g_nb_x_up, 110, 
	    g_cart_grid, &request[cntr]);
  cntr++;

  /* send the data to the neighbour on the right in x direction */
  /* recieve the data from the neighbour on the left in x direction */
  /* yx-edge */
  MPI_Isend(g_gauge_field[VOLUME + 2*LZ*(LX*LY + T*LY) + (LX-1)*LZ], 1, gauge_yx_edge_gath, g_nb_x_up, 111,
	    g_cart_grid, &request[cntr]);
  cntr++;
  MPI_Irecv(g_gauge_field[VOLUME + RAND + 4*LY*LZ + 2*T*LZ],         1, gauge_yx_edge_cont, g_nb_x_dn, 111,
	    g_cart_grid, &request[cntr]);
  cntr++;

  /* send the data to the neighbour on the left in y direction */
  /* recieve the data from the neighbour on the right in y direction */
  /* is on the t-Rand -> ty-edge*/
  MPI_Isend(g_gauge_field[VOLUME],                           1, gauge_ty_edge_gath, g_nb_y_dn, 112,
	    g_cart_grid, &request[cntr]);
  cntr++;
  MPI_Irecv(g_gauge_field[VOLUME + RAND + 4*LY*LZ + 4*T*LZ], 1, gauge_ty_edge_cont, g_nb_y_up, 112, 
	    g_cart_grid, &request[cntr]);
  cntr++;

  /* send the data to the neighbour on the right in y direction */
  /* recieve the data from the neighbour on the left in y direction */
  /* ty-edge */
  MPI_Isend(g_gauge_field[VOLUME + (LY-1)*LZ],                         1, gauge_ty_edge_gath, g_nb_y_up, 113,
	    g_cart_grid, &request[cntr]);
  cntr++;
  MPI_Irecv(g_gauge_field[VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 2*LX*LZ], 1, gauge_ty_edge_cont, g_nb_y_dn, 113,
	    g_cart_grid, &request[cntr]);
  cntr++;



  if(g_dbw2rand > 0) {
    /* send the data to the neighbour on the left in y direction */
    /* recieve the data from the neighbour on the right in y direction */
    /* x2y edge */
    MPI_Isend(g_gauge_field[VOLUME + 2*LX*LY*LZ + 2*T*LY*LZ + LZ],
	      1, gauge_yx_edge_gath, g_nb_x_dn, 114,
	      g_cart_grid, &request[cntr]);
    cntr++;
    MPI_Irecv(g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ],
	      1, gauge_yx_edge_cont, g_nb_x_up, 114, 
	      g_cart_grid, &request[cntr]);
    cntr++;
    
    /* send the data to the neighbour on the right in y direction */
    /* recieve the data from the neighbour on the left in y direction */
    /* x2y-edge */
    MPI_Isend(g_gauge_field[VOLUME + 2*LX*LY*LZ + 2*T*LY*LZ + (LX-2)*LZ],
	      1, gauge_yx_edge_gath, g_nb_x_up, 115,
	      g_cart_grid, &request[cntr]);
    cntr++;
    MPI_Irecv(g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 2*T*LZ],
	      1, gauge_yx_edge_cont, g_nb_x_dn, 115,
	      g_cart_grid, &request[cntr]);
    cntr++;


    /* send the data to the neighbour on the left in y direction */
    /* recieve the data from the neighbour on the right in y direction */
    /* y2x -edge */
    MPI_Isend(g_gauge_field[VOLUMEPLUSRAND + 2*LX*LY*LZ + 2*T*LY*LZ],
	      1, gauge_yx_edge_gath, g_nb_x_dn, 116,
	      g_cart_grid, &request[cntr]);
    cntr++;
    MPI_Irecv(g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 4*T*LZ],
	      1, gauge_yx_edge_cont, g_nb_x_up, 116, 
	      g_cart_grid, &request[cntr]);
    cntr++;
    
    /* send the data to the neighbour on the right in y direction */
    /* recieve the data from the neighbour on the left in y direction */
    /* y2x edge */
    MPI_Isend(g_gauge_field[VOLUMEPLUSRAND + 2*LX*LY*LZ + 2*T*LY*LZ + (LX-1)*LZ],
	      1, gauge_yx_edge_gath, g_nb_x_up, 117,
	      g_cart_grid, &request[cntr]);
    cntr++;
    MPI_Irecv(g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 6*T*LZ],
	      1, gauge_yx_edge_cont, g_nb_x_dn, 117,
	      g_cart_grid, &request[cntr]);
    cntr++;

    /* send the data to the neighbour on the left in y direction */
    /* recieve the data from the neighbour on the right in y direction */
    /* t2y-edge */
    MPI_Isend(g_gauge_field[VOLUMEPLUSRAND],
	      1, gauge_ty_edge_gath, g_nb_y_dn, 118,
	      g_cart_grid, &request[cntr]);
    cntr++;
    MPI_Irecv(g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ],
	      1, gauge_ty_edge_cont, g_nb_y_up, 118, 
	      g_cart_grid, &request[cntr]);
    cntr++;
    
    /* send the data to the neighbour on the right in y direction */
    /* recieve the data from the neighbour on the left in y direction */
    /* t2y edge */
    MPI_Isend(g_gauge_field[VOLUMEPLUSRAND + (LY-1)*LZ], 
	      1, gauge_ty_edge_gath, g_nb_y_up, 119,
	      g_cart_grid, &request[cntr]);
    cntr++;
    MPI_Irecv(g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 2*LX*LZ],
	      1, gauge_ty_edge_cont, g_nb_y_dn, 119,
	      g_cart_grid, &request[cntr]);
    cntr++;

    /* send the data to the neighbour on the left in y direction */
    /* recieve the data from the neighbour on the right in y direction */
    /* y2t edge */
    MPI_Isend(g_gauge_field[VOLUME + LZ],
	      1, gauge_ty_edge_gath, g_nb_y_dn, 120,
	      g_cart_grid, &request[cntr]);
    cntr++;
    MPI_Irecv(g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 4*LX*LZ],
	      1, gauge_ty_edge_cont, g_nb_y_up, 120, 
	      g_cart_grid, &request[cntr]);
    cntr++;
    
    /* send the data to the neighbour on the right in y direction */
    /* recieve the data from the neighbour on the left in y direction */
    /* y2t-edge */
    MPI_Isend(g_gauge_field[VOLUME + (LY-2)*LZ],
	      1, gauge_ty_edge_gath, g_nb_y_up, 121,
	      g_cart_grid, &request[cntr]);
    cntr++;
    MPI_Irecv(g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 6*LX*LZ],
	      1, gauge_ty_edge_cont, g_nb_y_dn, 121,
	      g_cart_grid, &request[cntr]);
    cntr++;
  }

  /* end of if defined PARALLELXYT || PARALLELXYZT */
#    endif
#    if defined PARALLELXYZT
  /* z-Rand */
  /* send the data to the neighbour on the left in z direction */
  /* recieve the data from the neighbour on the right in z direction */
  MPI_Isend(g_gauge_field[0],
	    1, gauge_z_slice_gath, g_nb_z_dn, 122,
	    g_cart_grid, &request[cntr]);
  cntr++;
  MPI_Irecv(g_gauge_field[VOLUME + 2*LZ*(LX*LY + T*LY) + 2*LZ*T*LX], 
	    1, gauge_z_slice_cont, g_nb_z_up, 122, 
	    g_cart_grid, &request[cntr]);
  cntr++;
  /* send the data to the neighbour on the right in z direction */
  /* recieve the data from the neighbour on the left in z direction */
  MPI_Isend(g_gauge_field[LZ-1],
	    1, gauge_z_slice_gath, g_nb_z_up, 123,
	    g_cart_grid, &request[cntr]);
  cntr++;
  MPI_Irecv(g_gauge_field[VOLUME + 2*LZ*(LX*LY + T*LY) + 2*T*LX*LZ + T*LX*LY],
	    1, gauge_z_slice_cont, g_nb_z_dn, 123,
	    g_cart_grid, &request[cntr]);
  cntr++;

  if(g_dbw2rand > 0) {
    /* send the data to the neighbour on the left in z direction */
    /* recieve the data from the neighbour on the right in z direction */
    /* z2-Rand */
    MPI_Isend(g_gauge_field[1],
	      1, gauge_z_slice_gath, g_nb_z_dn, 124,
	      g_cart_grid, &request[cntr]);
    cntr++;
    MPI_Irecv(g_gauge_field[VOLUMEPLUSRAND+(2*LX+2*T)*LY*LZ + 2*T*LX*LZ],
	      1, gauge_z_slice_cont, g_nb_z_up, 124, 
	      g_cart_grid, &request[cntr]);
    cntr++;
    /* send the data to the neighbour on the right in z direction */
    /* recieve the data from the neighbour on the left in z direction */
    /* z2-Rand */
    MPI_Isend(g_gauge_field[LZ-2],
	      1, gauge_z_slice_gath, g_nb_z_up, 125,
	      g_cart_grid, &request[cntr]);
    cntr++;
    MPI_Irecv(g_gauge_field[VOLUMEPLUSRAND+(2*LX+2*T)*LY*LZ + 2*T*LX*LZ + T*LX*LY],
	      1, gauge_z_slice_cont, g_nb_z_dn, 125,
	      g_cart_grid, &request[cntr]);
    cntr++;
  }
  MPI_Waitall(cntr, request, status);
  cntr=0;
  /* edges */

  /* send the data to the neighbour on the left in x direction */
  /* recieve the data from the neighbour on the right in x direction */
  /* is on the z-Rand -> zx-edge*/
  MPI_Isend(g_gauge_field[VOLUME + 2*LX*LY*LZ + 2*T*LY*LZ + 2*T*LX*LZ], 
	    1, gauge_zx_edge_gath, g_nb_x_dn, 126,
	    g_cart_grid, &request[cntr]);
  cntr++;
  MPI_Irecv(g_gauge_field[VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 4*LX*LZ],
	    1, gauge_zx_edge_cont, g_nb_x_up, 126, 
	    g_cart_grid, &request[cntr]);
  cntr++;

  /* send the data to the neighbour on the right in x direction */
  /* recieve the data from the neighbour on the left in x direction */
  /* zx-edge */
  MPI_Isend(g_gauge_field[VOLUME + 2*LZ*(LX*LY + T*LY) + 2*T*LX*LZ + (LX-1)*LY], 
	    1, gauge_zx_edge_gath, g_nb_x_up, 127,
	    g_cart_grid, &request[cntr]);
  cntr++;
  MPI_Irecv(g_gauge_field[VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 4*LX*LZ + 2*T*LY],
	    1, gauge_zx_edge_cont, g_nb_x_dn, 127,
	    g_cart_grid, &request[cntr]);
  cntr++;

  /* send the data to the neighbour on the left in z direction */
  /* recieve the data from the neighbour on the right in z direction */
  /* is on the t-Rand -> tz-edge*/
  MPI_Isend(g_gauge_field[VOLUME],
	    1, gauge_tz_edge_gath, g_nb_z_dn, 128,
	    g_cart_grid, &request[cntr]);
  cntr++;
  MPI_Irecv(g_gauge_field[VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 4*LX*LZ + 4*T*LY], 
	    1, gauge_tz_edge_cont, g_nb_z_up, 128, 
	    g_cart_grid, &request[cntr]);
  cntr++;
  
  /* send the data to the neighbour on the right in z direction */
  /* recieve the data from the neighbour on the left in z direction */
  /* tz-edge */
  MPI_Isend(g_gauge_field[VOLUME + (LZ-1)],
	    1, gauge_tz_edge_gath, g_nb_z_up, 129,
	    g_cart_grid, &request[cntr]);
  cntr++;
  MPI_Irecv(g_gauge_field[VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 4*LX*LZ + 4*T*LY + 2*LX*LY], 
	    1, gauge_tz_edge_cont, g_nb_z_dn, 129,
	    g_cart_grid, &request[cntr]);
  cntr++;

  /* send the data to the neighbour on the left in y direction */
  /* recieve the data from the neighbour on the right in y direction */
  /* is on the z-Rand -> zy-edge*/
  MPI_Isend(g_gauge_field[VOLUME + 2*LZ*(LX*LY + T*LY) + 2*T*LX*LZ],
	    1, gauge_zy_edge_gath, g_nb_y_dn, 130,
	    g_cart_grid, &request[cntr]);
  cntr++;
  MPI_Irecv(g_gauge_field[VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 4*LX*LZ + 4*T*LY + 4*LX*LY], 
	    1, gauge_zy_edge_cont, g_nb_y_up, 130, 
	    g_cart_grid, &request[cntr]);
  cntr++;

  /* send the data to the neighbour on the right in y direction */
  /* recieve the data from the neighbour on the left in y direction */
  /* zy-edge */
  MPI_Isend(g_gauge_field[VOLUME + 2*LZ*(LX*LY + T*LY) + 2*T*LX*LZ + (LY-1)],
	    1, gauge_zy_edge_gath, g_nb_y_up, 131,
	    g_cart_grid, &request[cntr]);
  cntr++;
  MPI_Irecv(g_gauge_field[VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 4*LX*LZ + 4*T*LY + 4*LX*LY + 2*T*LX], 
	    1, gauge_zy_edge_cont, g_nb_y_dn, 131,
	    g_cart_grid, &request[cntr]);
  cntr++;

  /* rectangular gauge action Stuff! */
  if(g_dbw2rand > 0) {
    /* send the data to the neighbour on the left in z direction */
    /* recieve the data from the neighbour on the right in z direction */
    /* t2z edge */
    MPI_Isend(g_gauge_field[VOLUMEPLUSRAND],
	      1, gauge_tz_edge_gath, g_nb_z_dn, 132,
	      g_cart_grid, &request[cntr]);
    cntr++;
    MPI_Irecv(g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ],
	      1, gauge_tz_edge_cont, g_nb_z_up, 132, 
	      g_cart_grid, &request[cntr]);
    cntr++;
    
    /* send the data to the neighbour on the right in z direction */
    /* recieve the data from the neighbour on the left in z direction */
    /* t2z-edge */
    MPI_Isend(g_gauge_field[VOLUMEPLUSRAND + (LZ-1)],
	      1, gauge_tz_edge_gath, g_nb_z_up, 133,
	      g_cart_grid, &request[cntr]);
    cntr++;
    MPI_Irecv(g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 2*LX*LY],
	      1, gauge_tz_edge_cont, g_nb_z_dn, 133,
	      g_cart_grid, &request[cntr]);
    cntr++;


    /* send the data to the neighbour on the left in z direction */
    /* recieve the data from the neighbour on the right in z direction */
    /* z2t -edge */
    MPI_Isend(g_gauge_field[VOLUME + 1],
	      1, gauge_tz_edge_gath, g_nb_z_dn, 134,
	      g_cart_grid, &request[cntr]);
    cntr++;
    MPI_Irecv(g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 4*LX*LY],
	      1, gauge_tz_edge_cont, g_nb_z_up, 134, 
	      g_cart_grid, &request[cntr]);
    cntr++;
    
    /* send the data to the neighbour on the right in z direction */
    /* recieve the data from the neighbour on the left in z direction */
    /* z2t edge */
    MPI_Isend(g_gauge_field[VOLUME + (LZ-2)],
	      1, gauge_tz_edge_gath, g_nb_z_up, 135,
	      g_cart_grid, &request[cntr]);
    cntr++;
    MPI_Irecv(g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 6*LX*LY],
	      1, gauge_tz_edge_cont, g_nb_z_dn, 135,
	      g_cart_grid, &request[cntr]);
    cntr++;

    /* send the data to the neighbour on the left in x direction */
    /* recieve the data from the neighbour on the right in x direction */
    /* z2x-edge */
    MPI_Isend(g_gauge_field[VOLUMEPLUSRAND + 2*LX*LY*LZ + 2*T*LY*LZ + 2*T*LX*LZ],
	      1, gauge_zx_edge_gath, g_nb_x_dn, 136,
	      g_cart_grid, &request[cntr]);
    cntr++;
    MPI_Irecv(g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY],
	      1, gauge_zx_edge_cont, g_nb_x_up, 136, 
	      g_cart_grid, &request[cntr]);
    cntr++;
    
    /* send the data to the neighbour on the right in x direction */
    /* recieve the data from the neighbour on the left in x direction */
    /* z2x edge */
    MPI_Isend(g_gauge_field[VOLUMEPLUSRAND + 2*LX*LY*LZ + 2*T*LY*LZ + 2*T*LX*LZ + (LX-1)*LY], 
	      1, gauge_zx_edge_gath, g_nb_x_up, 137,
	      g_cart_grid, &request[cntr]);
    cntr++;
    MPI_Irecv(g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 2*T*LY],
	      1, gauge_zx_edge_cont, g_nb_x_dn, 137,
	      g_cart_grid, &request[cntr]);
    cntr++;

    /* send the data to the neighbour on the left in x direction */
    /* recieve the data from the neighbour on the right in x direction */
    /* x2z edge */
    MPI_Isend(g_gauge_field[VOLUME + 2*LX*LY*LZ + 2*T*LY*LZ + 2*T*LX*LZ + LY],
	      1, gauge_zx_edge_gath, g_nb_x_dn, 138,
	      g_cart_grid, &request[cntr]);
    cntr++;
    MPI_Irecv(g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 4*T*LY],
	      1, gauge_zx_edge_cont, g_nb_x_up, 138, 
	      g_cart_grid, &request[cntr]);
    cntr++;
    
    /* send the data to the neighbour on the right in x direction */
    /* recieve the data from the neighbour on the left in x direction */
    /* x2z-edge */
    MPI_Isend(g_gauge_field[VOLUME + 2*LX*LY*LZ + 2*T*LY*LZ + 2*T*LX*LZ + (LX-2)*LY],
	      1, gauge_zx_edge_gath, g_nb_x_up, 139,
	      g_cart_grid, &request[cntr]);
    cntr++;
    MPI_Irecv(g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 6*T*LY],
	      1, gauge_zx_edge_cont, g_nb_x_dn, 139,
	      g_cart_grid, &request[cntr]);
    cntr++;

    /* send the data to the neighbour on the left in y direction */
    /* recieve the data from the neighbour on the right in y direction */
    /* z2y-edge */
    MPI_Isend(g_gauge_field[VOLUMEPLUSRAND + 2*LX*LY*LZ + 2*T*LY*LZ + 2*T*LX*LZ],
	      1, gauge_zy_edge_gath, g_nb_y_dn, 140,
	      g_cart_grid, &request[cntr]);
    cntr++;
    MPI_Irecv(g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 8*T*LY],
	      1, gauge_zy_edge_cont, g_nb_y_up, 140, 
	      g_cart_grid, &request[cntr]);
    cntr++;
    
    /* send the data to the neighbour on the right in y direction */
    /* recieve the data from the neighbour on the left in y direction */
    /* z2y edge */
    MPI_Isend(g_gauge_field[VOLUMEPLUSRAND + 2*LX*LY*LZ + 2*T*LY*LZ + 2*T*LX*LZ + (LY-1)], 
	      1, gauge_zy_edge_gath, g_nb_y_up, 141,
	      g_cart_grid, &request[cntr]);
    cntr++;
    MPI_Irecv(g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 8*T*LY + 2*T*LX],
	      1, gauge_zy_edge_cont, g_nb_y_dn, 141,
	      g_cart_grid, &request[cntr]);
    cntr++;

    /* send the data to the neighbour on the left in y direction */
    /* recieve the data from the neighbour on the right in y direction */
    /* y2z edge */
    MPI_Isend(g_gauge_field[VOLUME + 2*LX*LY*LZ + 2*T*LY*LZ + 2*T*LX*LZ + 1],
	      1, gauge_zy_edge_gath, g_nb_y_dn, 142,
	      g_cart_grid, &request[cntr]);
    cntr++;
    MPI_Irecv(g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 8*T*LY + 4*T*LX],
	      1, gauge_zy_edge_cont, g_nb_y_up, 142, 
	      g_cart_grid, &request[cntr]);
    cntr++;
    
    /* send the data to the neighbour on the right in y direction */
    /* recieve the data from the neighbour on the left in y direction */
    /* y2z-edge */
    MPI_Isend(g_gauge_field[VOLUME + 2*LX*LY*LZ + 2*T*LY*LZ + 2*T*LX*LZ + (LY-2)],
	      1, gauge_zy_edge_gath, g_nb_y_up, 143,
	      g_cart_grid, &request[cntr]);
    cntr++;
    MPI_Irecv(g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 8*T*LY + 6*T*LX],
	      1, gauge_zy_edge_cont, g_nb_y_dn, 143,
	      g_cart_grid, &request[cntr]);
    cntr++;
  }
  MPI_Waitall(cntr, request, status);

  /* end of if defined PARALLELXYZT */
#    endif
#  endif
  return;
}
/* BGL */
#else

void xchange_gauge() {

#ifdef MPI

  MPI_Status status;
  /* send the data to the neighbour on the left */
  /* recieve the data from the neighbour on the right */
  MPI_Sendrecv(g_gauge_field[0],          1, gauge_time_slice_cont, g_nb_t_dn, 83, 
	       g_gauge_field[VOLUME], 1, gauge_time_slice_cont, g_nb_t_up, 83, 
	       g_cart_grid, &status);

  /* send the data to the neighbour on the right */
  /* recieve the data from the neighbour on the left */
  MPI_Sendrecv(g_gauge_field[(T-1)*LX*LY*LZ], 1, gauge_time_slice_cont, g_nb_t_up, 84, 
	       g_gauge_field[(T+1)*LX*LY*LZ], 1, gauge_time_slice_cont, g_nb_t_dn, 84, 
	       g_cart_grid, &status);

  if(g_dbw2rand > 0) {
    /* send the data to the neighbour on the left */
    /* recieve the data from the neighbour on the right */
    /* t2-Rand */
    MPI_Sendrecv(g_gauge_field[1*LX*LY*LZ],     1, gauge_time_slice_cont, g_nb_t_dn, 85, 
		 g_gauge_field[VOLUMEPLUSRAND], 1, gauge_time_slice_cont, g_nb_t_up, 85, 
		 g_cart_grid, &status);

    /* send the data to the neighbour on the right */
    /* recieve the data from the neighbour on the left */
    /* t2-Rand */
    MPI_Sendrecv(g_gauge_field[(T-2)*LX*LY*LZ],          1, gauge_time_slice_cont, g_nb_t_up, 86, 
		 g_gauge_field[VOLUMEPLUSRAND+LX*LY*LZ], 1, gauge_time_slice_cont, g_nb_t_dn, 86, 
		 g_cart_grid, &status);
  }
  
#  if (defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT)
  /* send the data to the neighbour on the left in x direction */
  /* recieve the data from the neighbour on the right in x direction */
  MPI_Sendrecv(g_gauge_field[0],              1, gauge_x_slice_gath, g_nb_x_dn, 93,
	       g_gauge_field[(T+2)*LX*LY*LZ], 1, gauge_x_slice_cont, g_nb_x_up, 93, 
	       g_cart_grid, &status);

  /* send the data to the neighbour on the right in x direction */
  /* recieve the data from the neighbour on the left in x direction */
  /* x2-Rand */
  MPI_Sendrecv(g_gauge_field[(LX-1)*LY*LZ],             1, gauge_x_slice_gath, g_nb_x_up, 94,
	       g_gauge_field[(T+2)*LX*LY*LZ + T*LY*LZ], 1, gauge_x_slice_cont, g_nb_x_dn, 94,
	       g_cart_grid, &status);

  if(g_dbw2rand > 0) {
    /* send the data to the neighbour on the left in x direction */
    /* recieve the data from the neighbour on the right in x direction */
    /* x2-Rand */
    MPI_Sendrecv(g_gauge_field[LY*LZ],                     1, gauge_x_slice_gath, g_nb_x_dn, 95,
		 g_gauge_field[VOLUMEPLUSRAND+2*LX*LY*LZ], 1, gauge_x_slice_cont, g_nb_x_up, 95, 
		 g_cart_grid, &status);
    /* send the data to the neighbour on the right in x direction */
    /* recieve the data from the neighbour on the left in x direction */
    /* x2-Rand */
    MPI_Sendrecv(g_gauge_field[(LX-2)*LY*LZ],                        1, gauge_x_slice_gath, g_nb_x_up, 96,
		 g_gauge_field[VOLUMEPLUSRAND+2*LX*LY*LZ + T*LY*LZ], 1, gauge_x_slice_cont, g_nb_x_dn, 96,
		 g_cart_grid, &status);
  }

  /* The edges */

  /* send the data to the neighbour on the left in t direction */
  /* recieve the data from the neighbour on the right in t direction */
  /* is on the x-Rand: xt-edge */
  MPI_Sendrecv(g_gauge_field[(T+2)*LX*LY*LZ], 1, gauge_xt_edge_gath, g_nb_t_dn, 95,
	       g_gauge_field[VOLUME + RAND],  1, gauge_xt_edge_cont, g_nb_t_up, 95, 
	       g_cart_grid, &status);

  /* send the data to the neighbour on the right in t direction */
  /* recieve the data from the neighbour on the left in t direction */
  /* xt-edge */
  MPI_Sendrecv(g_gauge_field[(T+2)*LX*LY*LZ + (T-1)*LY*LZ], 1, gauge_xt_edge_gath, g_nb_t_up, 96,
	       g_gauge_field[VOLUME + RAND + 2*LY*LZ],      1, gauge_xt_edge_cont, g_nb_t_dn, 96,
	       g_cart_grid, &status);

  if(g_dbw2rand > 0) {
    /* send the data to the neighbour on the left in t direction */
    /* recieve the data from the neighbour on the right in t direction */
    /* t2x-edge */
    MPI_Sendrecv(g_gauge_field[(T+2)*LX*LY*LZ + LY*LZ],
		 1, gauge_xt_edge_gath, g_nb_t_dn, 97,
		 g_gauge_field[VOLUMEPLUSRAND + RAND],
		 1, gauge_xt_edge_cont, g_nb_t_up, 97, 
		 g_cart_grid, &status);
    
    /* send the data to the neighbour on the right in t direction */
    /* recieve the data from the neighbour on the left in t direction */
    /* t2x-edge */
    MPI_Sendrecv(g_gauge_field[(T+2)*LX*LY*LZ + (T-2)*LY*LZ],
		 1, gauge_xt_edge_gath, g_nb_t_up, 98,
		 g_gauge_field[VOLUMEPLUSRAND + RAND + 2*LY*LZ],
		 1, gauge_xt_edge_cont, g_nb_t_dn, 98,
		 g_cart_grid, &status);

    /* send the data to the neighbour on the left in t direction */
    /* recieve the data from the neighbour on the right in t direction */
    /* x2t-edge */
    MPI_Sendrecv(g_gauge_field[VOLUMEPLUSRAND + 2*LX*LY*LZ],
		 1, gauge_xt_edge_gath, g_nb_t_dn, 97,
		 g_gauge_field[VOLUMEPLUSRAND + RAND + 4*LY*LZ],
		 1, gauge_xt_edge_cont, g_nb_t_up, 97, 
		 g_cart_grid, &status);
    
    /* send the data to the neighbour on the right in t direction */
    /* recieve the data from the neighbour on the left in t direction */
    /* x2t-edge */
    MPI_Sendrecv(g_gauge_field[VOLUMEPLUSRAND + 2*LX*LY*LZ + (T-1)*LY*LZ],
		 1, gauge_xt_edge_gath, g_nb_t_up, 98,
		 g_gauge_field[VOLUMEPLUSRAND + RAND + 6*LY*LZ],
		 1, gauge_xt_edge_cont, g_nb_t_dn, 98,
		 g_cart_grid, &status);
  }
  /* end of if defined PARALLELXT || PARALLELXYT || PARALLELXYZT*/
#  endif

#  if (defined PARALLELXYT || defined PARALLELXYZT)
  /* send the data to the neighbour on the left in y direction */
  /* recieve the data from the neighbour on the right in y direction */
  MPI_Sendrecv(g_gauge_field[0],                            1, gauge_y_slice_gath, g_nb_y_dn, 103,
	       g_gauge_field[VOLUME + 2*LZ*(LX*LY + T*LY)], 1, gauge_y_slice_cont, g_nb_y_up, 103, 
	       g_cart_grid, &status);
  /* send the data to the neighbour on the right in y direction */
  /* recieve the data from the neighbour on the left in y direction */
  MPI_Sendrecv(g_gauge_field[(LY-1)*LZ],                              1, gauge_y_slice_gath, g_nb_y_up, 104,
	       g_gauge_field[VOLUME + 2*LZ*(LX*LY + T*LY) + T*LX*LZ], 1, gauge_y_slice_cont, g_nb_y_dn, 104,
	       g_cart_grid, &status);

  if(g_dbw2rand > 0) {
    /* send the data to the neighbour on the left in x direction */
    /* recieve the data from the neighbour on the right in x direction */
    /* y2-Rand */
    MPI_Sendrecv(g_gauge_field[LZ],                              1, gauge_y_slice_gath, g_nb_y_dn, 105,
		 g_gauge_field[VOLUMEPLUSRAND+(2*LX+2*T)*LY*LZ], 1, gauge_y_slice_cont, g_nb_y_up, 105, 
		 g_cart_grid, &status);
    /* send the data to the neighbour on the right in x direction */
    /* recieve the data from the neighbour on the left in x direction */
    /* y2-Rand */
    MPI_Sendrecv(g_gauge_field[(LY-2)*LZ],                                 1, gauge_y_slice_gath, g_nb_y_up, 106,
		 g_gauge_field[VOLUMEPLUSRAND+(2*LX+2*T)*LY*LZ + T*LX*LZ], 1, gauge_y_slice_cont, g_nb_y_dn, 106,
		 g_cart_grid, &status);
  }

  /* jetzt wirds richtig eklig ... */

  /* edges */

  /* send the data to the neighbour on the left in x direction */
  /* recieve the data from the neighbour on the right in x direction */
  /* is on the y-Rand -> yx-edge*/
  MPI_Sendrecv(g_gauge_field[VOLUME + 2*LX*LY*LZ + 2*T*LY*LZ], 1, gauge_yx_edge_gath, g_nb_x_dn, 107,
	       g_gauge_field[VOLUME + RAND + 4*LY*LZ],      1, gauge_yx_edge_cont, g_nb_x_up, 107, 
	       g_cart_grid, &status);

  /* send the data to the neighbour on the right in x direction */
  /* recieve the data from the neighbour on the left in x direction */
  /* yx-edge */
  MPI_Sendrecv(g_gauge_field[VOLUME + 2*LZ*(LX*LY + T*LY) + (LX-1)*LZ], 1, gauge_yx_edge_gath, g_nb_x_up, 108,
	       g_gauge_field[VOLUME + RAND + 4*LY*LZ + 2*T*LZ],         1, gauge_yx_edge_cont, g_nb_x_dn, 108,
	       g_cart_grid, &status);

  /* send the data to the neighbour on the left in y direction */
  /* recieve the data from the neighbour on the right in y direction */
  /* is on the t-Rand -> ty-edge*/
  MPI_Sendrecv(g_gauge_field[VOLUME],                           1, gauge_ty_edge_gath, g_nb_y_dn, 109,
	       g_gauge_field[VOLUME + RAND + 4*LY*LZ + 4*T*LZ], 1, gauge_ty_edge_cont, g_nb_y_up, 109, 
	       g_cart_grid, &status);

  /* send the data to the neighbour on the right in y direction */
  /* recieve the data from the neighbour on the left in y direction */
  /* ty-edge */
  MPI_Sendrecv(g_gauge_field[VOLUME + (LY-1)*LZ],                         1, gauge_ty_edge_gath, g_nb_y_up, 110,
	       g_gauge_field[VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 2*LX*LZ], 1, gauge_ty_edge_cont, g_nb_y_dn, 110,
	       g_cart_grid, &status);



  if(g_dbw2rand > 0) {
    /* send the data to the neighbour on the left in y direction */
    /* recieve the data from the neighbour on the right in y direction */
    /* x2y edge */
    MPI_Sendrecv(g_gauge_field[VOLUME + 2*LX*LY*LZ + 2*T*LY*LZ + LZ],
		 1, gauge_yx_edge_gath, g_nb_x_dn, 97,
		 g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ],
		 1, gauge_yx_edge_cont, g_nb_x_up, 97, 
		 g_cart_grid, &status);
    
    /* send the data to the neighbour on the right in y direction */
    /* recieve the data from the neighbour on the left in y direction */
    /* x2y-edge */
    MPI_Sendrecv(g_gauge_field[VOLUME + 2*LX*LY*LZ + 2*T*LY*LZ + (LX-2)*LZ],
		 1, gauge_yx_edge_gath, g_nb_x_up, 98,
		 g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 2*T*LZ],
		 1, gauge_yx_edge_cont, g_nb_x_dn, 98,
		 g_cart_grid, &status);


    /* send the data to the neighbour on the left in y direction */
    /* recieve the data from the neighbour on the right in y direction */
    /* y2x -edge */
    MPI_Sendrecv(g_gauge_field[VOLUMEPLUSRAND + 2*LX*LY*LZ + 2*T*LY*LZ],
		 1, gauge_yx_edge_gath, g_nb_x_dn, 97,
		 g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 4*T*LZ],
		 1, gauge_yx_edge_cont, g_nb_x_up, 97, 
		 g_cart_grid, &status);
    
    /* send the data to the neighbour on the right in y direction */
    /* recieve the data from the neighbour on the left in y direction */
    /* y2x edge */
    MPI_Sendrecv(g_gauge_field[VOLUMEPLUSRAND + 2*LX*LY*LZ + 2*T*LY*LZ + (LX-1)*LZ],
		 1, gauge_yx_edge_gath, g_nb_x_up, 98,
		 g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 6*T*LZ],
		 1, gauge_yx_edge_cont, g_nb_x_dn, 98,
		 g_cart_grid, &status);

    /* send the data to the neighbour on the left in y direction */
    /* recieve the data from the neighbour on the right in y direction */
    /* t2y-edge */
    MPI_Sendrecv(g_gauge_field[VOLUMEPLUSRAND],
		 1, gauge_ty_edge_gath, g_nb_y_dn, 197,
		 g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ],
		 1, gauge_ty_edge_cont, g_nb_y_up, 197, 
		 g_cart_grid, &status);
    
    /* send the data to the neighbour on the right in y direction */
    /* recieve the data from the neighbour on the left in y direction */
    /* t2y edge */
    MPI_Sendrecv(g_gauge_field[VOLUMEPLUSRAND + (LY-1)*LZ], 
		 1, gauge_ty_edge_gath, g_nb_y_up, 198,
		 g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 2*LX*LZ],
		 1, gauge_ty_edge_cont, g_nb_y_dn, 198,
		 g_cart_grid, &status);

    /* send the data to the neighbour on the left in y direction */
    /* recieve the data from the neighbour on the right in y direction */
    /* y2t edge */
    MPI_Sendrecv(g_gauge_field[VOLUME + LZ],
		 1, gauge_ty_edge_gath, g_nb_y_dn, 297,
		 g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 4*LX*LZ],
		 1, gauge_ty_edge_cont, g_nb_y_up, 297, 
		 g_cart_grid, &status);
    
    /* send the data to the neighbour on the right in y direction */
    /* recieve the data from the neighbour on the left in y direction */
    /* y2t-edge */
    MPI_Sendrecv(g_gauge_field[VOLUME + (LY-2)*LZ],
		 1, gauge_ty_edge_gath, g_nb_y_up, 298,
		 g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 6*LX*LZ],
		 1, gauge_ty_edge_cont, g_nb_y_dn, 298,
		 g_cart_grid, &status);
  }

  /* end of if defined PARALLELXYT || PARALLELXYZT */
#  endif
#  if defined PARALLELXYZT
  /* z-Rand */
  /* send the data to the neighbour on the left in z direction */
  /* recieve the data from the neighbour on the right in z direction */
  MPI_Sendrecv(g_gauge_field[0],
	       1, gauge_z_slice_gath, g_nb_z_dn, 303,
	       g_gauge_field[VOLUME + 2*LZ*(LX*LY + T*LY) + 2*LZ*T*LX], 
	       1, gauge_z_slice_cont, g_nb_z_up, 303, 
	       g_cart_grid, &status);
  /* send the data to the neighbour on the right in z direction */
  /* recieve the data from the neighbour on the left in z direction */
  MPI_Sendrecv(g_gauge_field[LZ-1],
	       1, gauge_z_slice_gath, g_nb_z_up, 304,
	       g_gauge_field[VOLUME + 2*LZ*(LX*LY + T*LY) + 2*T*LX*LZ + T*LX*LY],
	       1, gauge_z_slice_cont, g_nb_z_dn, 304,
	       g_cart_grid, &status);

  if(g_dbw2rand > 0) {
    /* send the data to the neighbour on the left in z direction */
    /* recieve the data from the neighbour on the right in z direction */
    /* z2-Rand */
    MPI_Sendrecv(g_gauge_field[1],
		 1, gauge_z_slice_gath, g_nb_z_dn, 305,
		 g_gauge_field[VOLUMEPLUSRAND+(2*LX+2*T)*LY*LZ + 2*T*LX*LZ],
		 1, gauge_z_slice_cont, g_nb_z_up, 305, 
		 g_cart_grid, &status);
    /* send the data to the neighbour on the right in z direction */
    /* recieve the data from the neighbour on the left in z direction */
    /* z2-Rand */
    MPI_Sendrecv(g_gauge_field[LZ-2],
		 1, gauge_z_slice_gath, g_nb_z_up, 306,
		 g_gauge_field[VOLUMEPLUSRAND+(2*LX+2*T)*LY*LZ + 2*T*LX*LZ + T*LX*LY],
		 1, gauge_z_slice_cont, g_nb_z_dn, 306,
		 g_cart_grid, &status);
  }

  /* edges */

  /* send the data to the neighbour on the left in x direction */
  /* recieve the data from the neighbour on the right in x direction */
  /* is on the z-Rand -> zx-edge*/
  MPI_Sendrecv(g_gauge_field[VOLUME + 2*LX*LY*LZ + 2*T*LY*LZ + 2*T*LX*LZ], 
	       1, gauge_zx_edge_gath, g_nb_x_dn, 307,
	       g_gauge_field[VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 4*LX*LZ],
	       1, gauge_zx_edge_cont, g_nb_x_up, 307, 
	       g_cart_grid, &status);

  /* send the data to the neighbour on the right in x direction */
  /* recieve the data from the neighbour on the left in x direction */
  /* zx-edge */
  MPI_Sendrecv(g_gauge_field[VOLUME + 2*LZ*(LX*LY + T*LY) + 2*T*LX*LZ + (LX-1)*LY], 
	       1, gauge_zx_edge_gath, g_nb_x_up, 308,
	       g_gauge_field[VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 4*LX*LZ + 2*T*LY],
	       1, gauge_zx_edge_cont, g_nb_x_dn, 308,
	       g_cart_grid, &status);

  /* send the data to the neighbour on the left in z direction */
  /* recieve the data from the neighbour on the right in z direction */
  /* is on the t-Rand -> tz-edge*/
  MPI_Sendrecv(g_gauge_field[VOLUME],
	       1, gauge_tz_edge_gath, g_nb_z_dn, 309,
	       g_gauge_field[VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 4*LX*LZ + 4*T*LY], 
	       1, gauge_tz_edge_cont, g_nb_z_up, 309, 
	       g_cart_grid, &status);

  /* send the data to the neighbour on the right in z direction */
  /* recieve the data from the neighbour on the left in z direction */
  /* tz-edge */
  MPI_Sendrecv(g_gauge_field[VOLUME + (LZ-1)],
	       1, gauge_tz_edge_gath, g_nb_z_up, 310,
	       g_gauge_field[VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 4*LX*LZ + 4*T*LY + 2*LX*LY], 
	       1, gauge_tz_edge_cont, g_nb_z_dn, 310,
	       g_cart_grid, &status);

  /* send the data to the neighbour on the left in y direction */
  /* recieve the data from the neighbour on the right in y direction */
  /* is on the z-Rand -> zy-edge*/
  MPI_Sendrecv(g_gauge_field[VOLUME + 2*LZ*(LX*LY + T*LY) + 2*T*LX*LZ],
	       1, gauge_zy_edge_gath, g_nb_y_dn, 310,
	       g_gauge_field[VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 4*LX*LZ + 4*T*LY + 4*LX*LY], 
	       1, gauge_zy_edge_cont, g_nb_y_up, 310, 
	       g_cart_grid, &status);

  /* send the data to the neighbour on the right in y direction */
  /* recieve the data from the neighbour on the left in y direction */
  /* zy-edge */
  MPI_Sendrecv(g_gauge_field[VOLUME + 2*LZ*(LX*LY + T*LY) + 2*T*LX*LZ + (LY-1)],
	       1, gauge_zy_edge_gath, g_nb_y_up, 310,
	       g_gauge_field[VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 4*LX*LZ + 4*T*LY + 4*LX*LY + 2*T*LX], 
	       1, gauge_zy_edge_cont, g_nb_y_dn, 310,
	       g_cart_grid, &status);

  /* rectangular gauge action Stuff! */
  if(g_dbw2rand > 0) {
    /* send the data to the neighbour on the left in z direction */
    /* recieve the data from the neighbour on the right in z direction */
    /* t2z edge */
    MPI_Sendrecv(g_gauge_field[VOLUMEPLUSRAND],
		 1, gauge_tz_edge_gath, g_nb_z_dn, 500,
		 g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ],
		 1, gauge_tz_edge_cont, g_nb_z_up, 500, 
		 g_cart_grid, &status);
    
    /* send the data to the neighbour on the right in z direction */
    /* recieve the data from the neighbour on the left in z direction */
    /* t2z-edge */
    MPI_Sendrecv(g_gauge_field[VOLUMEPLUSRAND + (LZ-1)],
		 1, gauge_tz_edge_gath, g_nb_z_up, 501,
		 g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 2*LX*LY],
		 1, gauge_tz_edge_cont, g_nb_z_dn, 501,
		 g_cart_grid, &status);


    /* send the data to the neighbour on the left in z direction */
    /* recieve the data from the neighbour on the right in z direction */
    /* z2t -edge */
    MPI_Sendrecv(g_gauge_field[VOLUME + 1],
		 1, gauge_tz_edge_gath, g_nb_z_dn, 502,
		 g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 4*LX*LY],
		 1, gauge_tz_edge_cont, g_nb_z_up, 502, 
		 g_cart_grid, &status);
    
    /* send the data to the neighbour on the right in z direction */
    /* recieve the data from the neighbour on the left in z direction */
    /* z2t edge */
    MPI_Sendrecv(g_gauge_field[VOLUME + (LZ-2)],
		 1, gauge_tz_edge_gath, g_nb_z_up, 503,
		 g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 6*LX*LY],
		 1, gauge_tz_edge_cont, g_nb_z_dn, 503,
		 g_cart_grid, &status);

    /* send the data to the neighbour on the left in x direction */
    /* recieve the data from the neighbour on the right in x direction */
    /* z2x-edge */
    MPI_Sendrecv(g_gauge_field[VOLUMEPLUSRAND + 2*LX*LY*LZ + 2*T*LY*LZ + 2*T*LX*LZ],
		 1, gauge_zx_edge_gath, g_nb_x_dn, 504,
		 g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY],
		 1, gauge_zx_edge_cont, g_nb_x_up, 504, 
		 g_cart_grid, &status);
    
    /* send the data to the neighbour on the right in x direction */
    /* recieve the data from the neighbour on the left in x direction */
    /* z2x edge */
    MPI_Sendrecv(g_gauge_field[VOLUMEPLUSRAND + 2*LX*LY*LZ + 2*T*LY*LZ + 2*T*LX*LZ + (LX-1)*LY], 
		 1, gauge_zx_edge_gath, g_nb_x_up, 504,
		 g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 2*T*LY],
		 1, gauge_zx_edge_cont, g_nb_x_dn, 504,
		 g_cart_grid, &status);

    /* send the data to the neighbour on the left in x direction */
    /* recieve the data from the neighbour on the right in x direction */
    /* x2z edge */
    MPI_Sendrecv(g_gauge_field[VOLUME + 2*LX*LY*LZ + 2*T*LY*LZ + 2*T*LX*LZ + LY],
		 1, gauge_zx_edge_gath, g_nb_x_dn, 505,
		 g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 4*T*LY],
		 1, gauge_zx_edge_cont, g_nb_x_up, 505, 
		 g_cart_grid, &status);
    
    /* send the data to the neighbour on the right in x direction */
    /* recieve the data from the neighbour on the left in x direction */
    /* x2z-edge */
    MPI_Sendrecv(g_gauge_field[VOLUME + 2*LX*LY*LZ + 2*T*LY*LZ + 2*T*LX*LZ + (LX-2)*LY],
		 1, gauge_zx_edge_gath, g_nb_x_up, 506,
		 g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 6*T*LY],
		 1, gauge_zx_edge_cont, g_nb_x_dn, 506,
		 g_cart_grid, &status);

    /* send the data to the neighbour on the left in y direction */
    /* recieve the data from the neighbour on the right in y direction */
    /* z2y-edge */
    MPI_Sendrecv(g_gauge_field[VOLUMEPLUSRAND + 2*LX*LY*LZ + 2*T*LY*LZ + 2*T*LX*LZ],
		 1, gauge_zy_edge_gath, g_nb_y_dn, 507,
		 g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 8*T*LY],
		 1, gauge_zy_edge_cont, g_nb_y_up, 507, 
		 g_cart_grid, &status);
    
    /* send the data to the neighbour on the right in y direction */
    /* recieve the data from the neighbour on the left in y direction */
    /* z2y edge */
    MPI_Sendrecv(g_gauge_field[VOLUMEPLUSRAND + 2*LX*LY*LZ + 2*T*LY*LZ + 2*T*LX*LZ + (LY-1)], 
		 1, gauge_zy_edge_gath, g_nb_y_up, 508,
		 g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 8*T*LY + 2*T*LX],
		 1, gauge_zy_edge_cont, g_nb_y_dn, 508,
		 g_cart_grid, &status);

    /* send the data to the neighbour on the left in y direction */
    /* recieve the data from the neighbour on the right in y direction */
    /* y2z edge */
    MPI_Sendrecv(g_gauge_field[VOLUME + 2*LX*LY*LZ + 2*T*LY*LZ + 2*T*LX*LZ + 1],
		 1, gauge_zy_edge_gath, g_nb_y_dn, 509,
		 g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 8*T*LY + 4*T*LX],
		 1, gauge_zy_edge_cont, g_nb_y_up, 509, 
		 g_cart_grid, &status);
    
    /* send the data to the neighbour on the right in y direction */
    /* recieve the data from the neighbour on the left in y direction */
    /* y2z-edge */
    MPI_Sendrecv(g_gauge_field[VOLUME + 2*LX*LY*LZ + 2*T*LY*LZ + 2*T*LX*LZ + (LY-2)],
		 1, gauge_zy_edge_gath, g_nb_y_up, 510,
		 g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 8*T*LY + 6*T*LX],
		 1, gauge_zy_edge_cont, g_nb_y_dn, 510,
		 g_cart_grid, &status);

  }


  /* end of if defined PARALLELXYZT */
#  endif
#endif
  return;
}
/* BGL */
#endif 
static char const rcsid[] = "$Id$";

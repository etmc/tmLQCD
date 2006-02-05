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

void xchange_gauge() {

#ifdef MPI

  MPI_Status status;
  /* send the data to the neighbour on the left */
  /* recieve the data from the neighbour on the right */
#  ifndef _NEW_GEOMETRY
  MPI_Sendrecv(g_gauge_field[0],          1, gauge_time_slice_cont, g_nb_t_dn, 83, 
	       g_gauge_field[VOLUME], 1, gauge_time_slice_cont, g_nb_t_up, 83, 
	       g_cart_grid, &status);
#  else
  MPI_Sendrecv(g_gauge_field[0],          1, gauge_time_slice_split, g_nb_t_dn, 83,
	       g_gauge_field[VOLUME], 1, gauge_time_slice_cont , g_nb_t_up, 83,
	       g_cart_grid, &status);
#  endif
  /* send the data to the neighbour on the right */
  /* recieve the data from the neighbour on the left */
# ifndef _NEW_GEOMETRY
  MPI_Sendrecv(g_gauge_field[(T-1)*LX*LY*LZ], 1, gauge_time_slice_cont, g_nb_t_up, 84, 
	       g_gauge_field[(T+1)*LX*LY*LZ], 1, gauge_time_slice_cont, g_nb_t_dn, 84, 
	       g_cart_grid, &status);
#  else
  MPI_Sendrecv(g_gauge_field[(T-1)*LX*LY*LZ/2], 1, gauge_time_slice_split, g_nb_t_up, 84, 
	       g_gauge_field[(T+1)*LX*LY*LZ],   1, gauge_time_slice_cont , g_nb_t_dn, 84,
	       g_cart_grid, &status);
#  endif

  if(g_dbw2rand > 0) {
    /* send the data to the neighbour on the left */
    /* recieve the data from the neighbour on the right */
    /* t2-Rand */
#  ifndef _NEW_GEOMETRY
    MPI_Sendrecv(g_gauge_field[1*LX*LY*LZ],     1, gauge_time_slice_cont, g_nb_t_dn, 85, 
		 g_gauge_field[VOLUMEPLUSRAND], 1, gauge_time_slice_cont, g_nb_t_up, 85, 
		 g_cart_grid, &status);
#  else
    MPI_Sendrecv(g_gauge_field[LX*LY*LZ/2],     1, gauge_time_slice_split, g_nb_t_dn, 85,
		 g_gauge_field[VOLUMEPLUSRAND], 1, gauge_time_slice_cont , g_nb_t_up, 85,
		 g_cart_grid, &status);
#  endif
    /* send the data to the neighbour on the right */
    /* recieve the data from the neighbour on the left */
    /* t2-Rand */
#  ifndef _NEW_GEOMETRY
    MPI_Sendrecv(g_gauge_field[(T-2)*LX*LY*LZ],          1, gauge_time_slice_cont, g_nb_t_up, 86, 
		 g_gauge_field[VOLUMEPLUSRAND+LX*LY*LZ], 1, gauge_time_slice_cont, g_nb_t_dn, 86, 
		 g_cart_grid, &status);
#  else
    MPI_Sendrecv(g_gauge_field[(T-2)*LX*LY*LZ/2],        1, gauge_time_slice_split, g_nb_t_up, 86, 
		 g_gauge_field[VOLUMEPLUSRAND+LX*LY*LZ], 1, gauge_time_slice_cont , g_nb_t_dn, 86,
		 g_cart_grid, &status);
#  endif
  }
  
#  if (defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT)
  /* send the data to the neighbour on the left in x direction */
  /* recieve the data from the neighbour on the right in x direction */
#    ifndef _NEW_GEOMETRY
  MPI_Sendrecv(g_gauge_field[0],              1, gauge_x_slice_gath, g_nb_x_dn, 93,
	       g_gauge_field[(T+2)*LX*LY*LZ], 1, gauge_x_slice_cont, g_nb_x_up, 93, 
	       g_cart_grid, &status);
#    else
  MPI_Sendrecv(g_gauge_field[0],              1, gauge_x_slice_gath_split, g_nb_x_dn, 93,
	       g_gauge_field[(T+2)*LX*LY*LZ], 1, gauge_x_slice_cont      , g_nb_x_up, 93, 
	       g_cart_grid, &status);
#    endif
  /* send the data to the neighbour on the right in x direction */
  /* recieve the data from the neighbour on the left in x direction */
  /* x2-Rand */
#    ifndef _NEW_GEOMETRY
  MPI_Sendrecv(g_gauge_field[(LX-1)*LY*LZ],             1, gauge_x_slice_gath, g_nb_x_up, 94,
	       g_gauge_field[(T+2)*LX*LY*LZ + T*LY*LZ], 1, gauge_x_slice_cont, g_nb_x_dn, 94,
	       g_cart_grid, &status);
#    else
  MPI_Sendrecv(g_gauge_field[(LX-1)*LY*LZ/2],           1, gauge_x_slice_gath_split, g_nb_x_up, 94,
	       g_gauge_field[(T+2)*LX*LY*LZ + T*LY*LZ], 1, gauge_x_slice_cont      , g_nb_x_dn, 94,
	       g_cart_grid, &status);
#    endif
  if(g_dbw2rand > 0) {
    /* send the data to the neighbour on the left in x direction */
    /* recieve the data from the neighbour on the right in x direction */
    /* x2-Rand */
#    ifndef _NEW_GEOMETRY
    MPI_Sendrecv(g_gauge_field[LY*LZ],                     1, gauge_x_slice_gath, g_nb_x_dn, 95,
		 g_gauge_field[VOLUMEPLUSRAND+2*LX*LY*LZ], 1, gauge_x_slice_cont, g_nb_x_up, 95, 
		 g_cart_grid, &status);
#    else
    MPI_Sendrecv(g_gauge_field[LY*LZ/2],                   1, gauge_x_slice_gath_split, g_nb_x_dn, 95,
		 g_gauge_field[VOLUMEPLUSRAND+2*LX*LY*LZ], 1, gauge_x_slice_cont      , g_nb_x_up, 95, 
		 g_cart_grid, &status);
#    endif
    /* send the data to the neighbour on the right in x direction */
    /* recieve the data from the neighbour on the left in x direction */
    /* x2-Rand */
#    ifndef _NEW_GEOMETRY
    MPI_Sendrecv(g_gauge_field[(LX-2)*LY*LZ],                        1, gauge_x_slice_gath, g_nb_x_up, 96,
		 g_gauge_field[VOLUMEPLUSRAND+2*LX*LY*LZ + T*LY*LZ], 1, gauge_x_slice_cont, g_nb_x_dn, 96,
		 g_cart_grid, &status);
#    else
    MPI_Sendrecv(g_gauge_field[(LX-2)*LY*LZ/2],                      1, gauge_x_slice_gath_split, g_nb_x_up, 96,
		 g_gauge_field[VOLUMEPLUSRAND+2*LX*LY*LZ + T*LY*LZ], 1, gauge_x_slice_cont      , g_nb_x_dn, 96,
		 g_cart_grid, &status);
#    endif
  }

  /* The edges */

  /* send the data to the neighbour on the left in t direction */
  /* recieve the data from the neighbour on the right in t direction */
  /* is on the x-Rand: xt-edge */
#    ifndef _NEW_GEOMETRY
  MPI_Sendrecv(g_gauge_field[(T+2)*LX*LY*LZ], 1, gauge_xt_edge_gath, g_nb_t_dn, 95,
	       g_gauge_field[VOLUME + RAND],  1, gauge_xt_edge_cont, g_nb_t_up, 95, 
	       g_cart_grid, &status);
#    else
  MPI_Sendrecv(g_gauge_field[(T+2)*LX*LY*LZ], 1, gauge_xt_edge_gath_split, g_nb_t_dn, 95,
	       g_gauge_field[VOLUME + RAND],  1, gauge_xt_edge_cont      , g_nb_t_up, 95, 
	       g_cart_grid, &status);
#    endif

  /* send the data to the neighbour on the right in t direction */
  /* recieve the data from the neighbour on the left in t direction */
  /* xt-edge */
#    ifndef _NEW_GEOMETRY
  MPI_Sendrecv(g_gauge_field[(T+2)*LX*LY*LZ + (T-1)*LY*LZ], 1, gauge_xt_edge_gath, g_nb_t_up, 96,
	       g_gauge_field[VOLUME + RAND + 2*LY*LZ],      1, gauge_xt_edge_cont, g_nb_t_dn, 96,
	       g_cart_grid, &status);
#    else
  MPI_Sendrecv(g_gauge_field[(T+2)*LX*LY*LZ + (T-1)*LY*LZ/2], 1, gauge_xt_edge_gath_split, g_nb_t_up, 96,
	       g_gauge_field[VOLUME + RAND + 2*LY*LZ],        1, gauge_xt_edge_cont      , g_nb_t_dn, 96,
	       g_cart_grid, &status);
#    endif

  if(g_dbw2rand > 0) {
    /* send the data to the neighbour on the left in t direction */
    /* recieve the data from the neighbour on the right in t direction */
    /* t2x-edge */
#    ifndef _NEW_GEOMETRY
    MPI_Sendrecv(g_gauge_field[(T+2)*LX*LY*LZ + LY*LZ],
		 1, gauge_xt_edge_gath, g_nb_t_dn, 97,
		 g_gauge_field[VOLUMEPLUSRAND + RAND],
		 1, gauge_xt_edge_cont, g_nb_t_up, 97, 
		 g_cart_grid, &status);
#    else
    MPI_Sendrecv(g_gauge_field[(T+2)*LX*LY*LZ + LY*LZ/2],
		 1, gauge_xt_edge_gath_split, g_nb_t_dn, 97,
		 g_gauge_field[VOLUMEPLUSRAND + RAND],
		 1, gauge_xt_edge_cont      , g_nb_t_up, 97, 
		 g_cart_grid, &status);
#    endif
    
    /* send the data to the neighbour on the right in t direction */
    /* recieve the data from the neighbour on the left in t direction */
    /* t2x-edge */
#    ifndef _NEW_GEOMETRY
    MPI_Sendrecv(g_gauge_field[(T+2)*LX*LY*LZ + (T-2)*LY*LZ],
		 1, gauge_xt_edge_gath, g_nb_t_up, 98,
		 g_gauge_field[VOLUMEPLUSRAND + RAND + 2*LY*LZ],
		 1, gauge_xt_edge_cont, g_nb_t_dn, 98,
		 g_cart_grid, &status);
#    else
    MPI_Sendrecv(g_gauge_field[(T+2)*LX*LY*LZ + (T-2)*LY*LZ/2],
		 1, gauge_xt_edge_gath_split, g_nb_t_up, 98,
		 g_gauge_field[VOLUMEPLUSRAND + RAND + 2*LY*LZ],
		 1, gauge_xt_edge_cont      , g_nb_t_dn, 98,
		 g_cart_grid, &status);
#    endif    

    /* send the data to the neighbour on the left in t direction */
    /* recieve the data from the neighbour on the right in t direction */
    /* x2t-edge */
#    ifndef _NEW_GEOMETRY
    MPI_Sendrecv(g_gauge_field[VOLUMEPLUSRAND + 2*LX*LY*LZ],
		 1, gauge_xt_edge_gath, g_nb_t_dn, 97,
		 g_gauge_field[VOLUMEPLUSRAND + RAND + 4*LY*LZ],
		 1, gauge_xt_edge_cont, g_nb_t_up, 97, 
		 g_cart_grid, &status);
#    else
    MPI_Sendrecv(g_gauge_field[VOLUMEPLUSRAND + 2*LX*LY*LZ],
		 1, gauge_xt_edge_gath_split, g_nb_t_dn, 97,
		 g_gauge_field[VOLUMEPLUSRAND + RAND + 4*LY*LZ],
		 1, gauge_xt_edge_cont      , g_nb_t_up, 97, 
		 g_cart_grid, &status);
#    endif
    
    /* send the data to the neighbour on the right in t direction */
    /* recieve the data from the neighbour on the left in t direction */
    /* x2t-edge */
#    ifndef _NEW_GEOMETRY
    MPI_Sendrecv(g_gauge_field[VOLUMEPLUSRAND + 2*LX*LY*LZ + (T-1)*LY*LZ],
		 1, gauge_xt_edge_gath, g_nb_t_up, 98,
		 g_gauge_field[VOLUMEPLUSRAND + RAND + 6*LY*LZ],
		 1, gauge_xt_edge_cont, g_nb_t_dn, 98,
		 g_cart_grid, &status);
#    else
    MPI_Sendrecv(g_gauge_field[VOLUMEPLUSRAND + 2*LX*LY*LZ + (T-1)*LY*LZ/2],
		 1, gauge_xt_edge_gath_split, g_nb_t_up, 98,
		 g_gauge_field[VOLUMEPLUSRAND + RAND + 6*LY*LZ],
		 1, gauge_xt_edge_cont      , g_nb_t_dn, 98,
		 g_cart_grid, &status);
#    endif    
  }
  /* end of if defined PARALLELXT || PARALLELXYT || PARALLELXYZT*/
#  endif

#  if (defined PARALLELXYT || defined PARALLELXYZT)
  /* send the data to the neighbour on the left in y direction */
  /* recieve the data from the neighbour on the right in y direction */
#    ifndef _NEW_GEOMETRY
  MPI_Sendrecv(g_gauge_field[0],                            1, gauge_y_slice_gath, g_nb_y_dn, 103,
	       g_gauge_field[VOLUME + 2*LZ*(LX*LY + T*LY)], 1, gauge_y_slice_cont, g_nb_y_up, 103, 
	       g_cart_grid, &status);
#    else
  MPI_Sendrecv(g_gauge_field[0],                            1, gauge_y_slice_gath_split, g_nb_y_dn, 103,
	       g_gauge_field[VOLUME + 2*LZ*(LX*LY + T*LY)], 1, gauge_y_slice_cont      , g_nb_y_up, 103, 
	       g_cart_grid, &status);
#    endif
  /* send the data to the neighbour on the right in y direction */
  /* recieve the data from the neighbour on the left in y direction */
#    ifndef _NEW_GEOMETRY
  MPI_Sendrecv(g_gauge_field[(LY-1)*LZ],                              1, gauge_y_slice_gath, g_nb_y_up, 104,
	       g_gauge_field[VOLUME + 2*LZ*(LX*LY + T*LY) + T*LX*LZ], 1, gauge_y_slice_cont, g_nb_y_dn, 104,
	       g_cart_grid, &status);
#    else
  MPI_Sendrecv(g_gauge_field[(LY-1)*LZ/2],                            1, gauge_y_slice_gath_split, g_nb_y_up, 104,
	       g_gauge_field[VOLUME + 2*LZ*(LX*LY + T*LY) + T*LX*LZ], 1, gauge_y_slice_cont      , g_nb_y_dn, 104,
	       g_cart_grid, &status);
#    endif

  if(g_dbw2rand > 0) {
    /* send the data to the neighbour on the left in x direction */
    /* recieve the data from the neighbour on the right in x direction */
    /* y2-Rand */
#    ifndef _NEW_GEOMETRY
    MPI_Sendrecv(g_gauge_field[LZ],                              1, gauge_y_slice_gath, g_nb_y_dn, 105,
		 g_gauge_field[VOLUMEPLUSRAND+(2*LX+2*T)*LY*LZ], 1, gauge_y_slice_cont, g_nb_y_up, 105, 
		 g_cart_grid, &status);
#    else
    MPI_Sendrecv(g_gauge_field[LZ/2],                            1, gauge_y_slice_gath_split, g_nb_y_dn, 105,
		 g_gauge_field[VOLUMEPLUSRAND+(2*LX+2*T)*LY*LZ], 1, gauge_y_slice_cont      , g_nb_y_up, 105, 
		 g_cart_grid, &status);
#    endif
    /* send the data to the neighbour on the right in x direction */
    /* recieve the data from the neighbour on the left in x direction */
    /* y2-Rand */
#    ifndef _NEW_GEOMETRY
    MPI_Sendrecv(g_gauge_field[(LY-2)*LZ],                                 1, gauge_y_slice_gath, g_nb_y_up, 106,
		 g_gauge_field[VOLUMEPLUSRAND+(2*LX+2*T)*LY*LZ + T*LX*LZ], 1, gauge_y_slice_cont, g_nb_y_dn, 106,
		 g_cart_grid, &status);
#    else
    MPI_Sendrecv(g_gauge_field[(LY-2)*LZ/2],                               1, gauge_y_slice_gath_split, g_nb_y_up, 106,
		 g_gauge_field[VOLUMEPLUSRAND+(2*LX+2*T)*LY*LZ + T*LX*LZ], 1, gauge_y_slice_cont      , g_nb_y_dn, 106,
		 g_cart_grid, &status);
#    endif
  }

  /* jetzt wirds richtig eklig ... */

  /* edges */

  /* send the data to the neighbour on the left in x direction */
  /* recieve the data from the neighbour on the right in x direction */
  /* is on the y-Rand -> yx-edge*/
#    ifndef _NEW_GEOMETRY
  MPI_Sendrecv(g_gauge_field[VOLUME + 2*LX*LY*LZ + 2*T*LY*LZ], 1, gauge_yx_edge_gath, g_nb_x_dn, 107,
	       g_gauge_field[VOLUME + RAND + 4*LY*LZ],      1, gauge_yx_edge_cont, g_nb_x_up, 107, 
	       g_cart_grid, &status);
#    else
  MPI_Sendrecv(g_gauge_field[VOLUME + 2*LX*LY*LZ + 2*T*LY*LZ], 1, gauge_yx_edge_gath_split, g_nb_x_dn, 107,
	       g_gauge_field[VOLUME + RAND + 4*LY*LZ],       1, gauge_yx_edge_cont      , g_nb_x_up, 107, 
	       g_cart_grid, &status);
#    endif

  /* send the data to the neighbour on the right in x direction */
  /* recieve the data from the neighbour on the left in x direction */
  /* yx-edge */
#    ifndef _NEW_GEOMETRY
  MPI_Sendrecv(g_gauge_field[VOLUME + 2*LZ*(LX*LY + T*LY) + (LX-1)*LZ], 1, gauge_yx_edge_gath, g_nb_x_up, 108,
	       g_gauge_field[VOLUME + RAND + 4*LY*LZ + 2*T*LZ],         1, gauge_yx_edge_cont, g_nb_x_dn, 108,
	       g_cart_grid, &status);
#    else
  MPI_Sendrecv(g_gauge_field[VOLUME + 2*LZ*(LX*LY + T*LY) + (LX-1)*LZ/2], 1, gauge_yx_edge_gath_split, g_nb_x_up, 108,
	       g_gauge_field[VOLUME + RAND + 4*LY*LZ + 2*T*LZ],           1, gauge_yx_edge_cont      , g_nb_x_dn, 108,
	       g_cart_grid, &status);
#    endif

  /* send the data to the neighbour on the left in y direction */
  /* recieve the data from the neighbour on the right in y direction */
  /* is on the t-Rand -> ty-edge*/
#    ifndef _NEW_GEOMETRY
  MPI_Sendrecv(g_gauge_field[VOLUME],                           1, gauge_ty_edge_gath, g_nb_y_dn, 109,
	       g_gauge_field[VOLUME + RAND + 4*LY*LZ + 4*T*LZ], 1, gauge_ty_edge_cont, g_nb_y_up, 109, 
	       g_cart_grid, &status);
#    else
  MPI_Sendrecv(g_gauge_field[VOLUME],                           1, gauge_ty_edge_gath_split, g_nb_y_dn, 109,
	       g_gauge_field[VOLUME + RAND + 4*LY*LZ + 4*T*LZ], 1, gauge_ty_edge_cont      , g_nb_y_up, 109, 
	       g_cart_grid, &status);
#    endif

  /* send the data to the neighbour on the right in y direction */
  /* recieve the data from the neighbour on the left in y direction */
  /* ty-edge */
#    ifndef _NEW_GEOMETRY
  MPI_Sendrecv(g_gauge_field[VOLUME + (LY-1)*LZ],                         1, gauge_ty_edge_gath, g_nb_y_up, 110,
	       g_gauge_field[VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 2*LX*LZ], 1, gauge_ty_edge_cont, g_nb_y_dn, 110,
	       g_cart_grid, &status);
#    else
  MPI_Sendrecv(g_gauge_field[VOLUME + (LY-1)*LZ/2],                      1, gauge_ty_edge_gath_split, g_nb_y_up, 110,
	       g_gauge_field[VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 2*LX*LZ], 1, gauge_ty_edge_cont      , g_nb_y_dn, 110,
	       g_cart_grid, &status);
#    endif



  if(g_dbw2rand > 0) {
    /* send the data to the neighbour on the left in y direction */
    /* recieve the data from the neighbour on the right in y direction */
    /* x2y edge */
#    ifndef _NEW_GEOMETRY
    MPI_Sendrecv(g_gauge_field[VOLUME + 2*LX*LY*LZ + 2*T*LY*LZ + LZ],
		 1, gauge_yx_edge_gath, g_nb_x_dn, 97,
		 g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ],
		 1, gauge_yx_edge_cont, g_nb_x_up, 97, 
		 g_cart_grid, &status);
#    else
    MPI_Sendrecv(g_gauge_field[VOLUME + 2*LX*LY*LZ + 2*T*LY*LZ + LZ],
		 1, gauge_yx_edge_gath_split, g_nb_x_dn, 97,
		 g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ],
		 1, gauge_yx_edge_cont      , g_nb_x_up, 97, 
		 g_cart_grid, &status);
#    endif
    
    /* send the data to the neighbour on the right in y direction */
    /* recieve the data from the neighbour on the left in y direction */
    /* x2y-edge */
#    ifndef _NEW_GEOMETRY
    MPI_Sendrecv(g_gauge_field[VOLUME + 2*LX*LY*LZ + 2*T*LY*LZ + (LX-2)*LZ],
		 1, gauge_yx_edge_gath, g_nb_x_up, 98,
		 g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 2*T*LZ],
		 1, gauge_yx_edge_cont, g_nb_x_dn, 98,
		 g_cart_grid, &status);
#    else
    MPI_Sendrecv(g_gauge_field[VOLUME + 2*LX*LY*LZ + 2*T*LY*LZ + (LX-2)*LZ/2],
		 1, gauge_yx_edge_gath_split, g_nb_x_up, 98,
		 g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 2*T*LZ],
		 1, gauge_yx_edge_cont      , g_nb_x_dn, 98,
		 g_cart_grid, &status);
#    endif


    /* send the data to the neighbour on the left in y direction */
    /* recieve the data from the neighbour on the right in y direction */
    /* y2x -edge */
#    ifndef _NEW_GEOMETRY
    MPI_Sendrecv(g_gauge_field[VOLUMEPLUSRAND + 2*LX*LY*LZ + 2*T*LY*LZ],
		 1, gauge_yx_edge_gath, g_nb_x_dn, 97,
		 g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 4*T*LZ],
		 1, gauge_yx_edge_cont, g_nb_x_up, 97, 
		 g_cart_grid, &status);
#    else
    MPI_Sendrecv(g_gauge_field[VOLUMEPLUSRAND + 2*LX*LY*LZ + 2*T*LY*LZ],
		 1, gauge_yx_edge_gath_split, g_nb_x_dn, 97,
		 g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 4*T*LZ],
		 1, gauge_yx_edge_cont      , g_nb_x_up, 97, 
		 g_cart_grid, &status);
#    endif
    
    /* send the data to the neighbour on the right in y direction */
    /* recieve the data from the neighbour on the left in y direction */
    /* y2x edge */
#    ifndef _NEW_GEOMETRY
    MPI_Sendrecv(g_gauge_field[VOLUMEPLUSRAND + 2*LX*LY*LZ + 2*T*LY*LZ + (LX-1)*LZ],
		 1, gauge_yx_edge_gath, g_nb_x_up, 98,
		 g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 6*T*LZ],
		 1, gauge_yx_edge_cont, g_nb_x_dn, 98,
		 g_cart_grid, &status);
#    else
    MPI_Sendrecv(g_gauge_field[VOLUMEPLUSRAND + 2*LX*LY*LZ + 2*T*LY*LZ + (LX-1)*LZ/2],
		 1, gauge_yx_edge_gath_split, g_nb_x_up, 98,
		 g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 6*T*LZ],
		 1, gauge_yx_edge_cont      , g_nb_x_dn, 98,
		 g_cart_grid, &status);
#    endif    


    /* send the data to the neighbour on the left in y direction */
    /* recieve the data from the neighbour on the right in y direction */
    /* t2y-edge */
#    ifndef _NEW_GEOMETRY
    MPI_Sendrecv(g_gauge_field[VOLUMEPLUSRAND],
		 1, gauge_ty_edge_gath, g_nb_y_dn, 197,
		 g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ],
		 1, gauge_ty_edge_cont, g_nb_y_up, 197, 
		 g_cart_grid, &status);
#    else
    MPI_Sendrecv(g_gauge_field[VOLUMEPLUSRAND],
		 1, gauge_ty_edge_gath_split, g_nb_y_dn, 197,
		 g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ],
		 1, gauge_ty_edge_cont      , g_nb_y_up, 197, 
		 g_cart_grid, &status);
#    endif
    
    /* send the data to the neighbour on the right in y direction */
    /* recieve the data from the neighbour on the left in y direction */
    /* t2y edge */
#    ifndef _NEW_GEOMETRY
    MPI_Sendrecv(g_gauge_field[VOLUMEPLUSRAND + (LY-1)*LZ], 
		 1, gauge_ty_edge_gath, g_nb_y_up, 198,
		 g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 2*LX*LZ],
		 1, gauge_ty_edge_cont, g_nb_y_dn, 198,
		 g_cart_grid, &status);
#    else
    MPI_Sendrecv(g_gauge_field[VOLUMEPLUSRAND + (LY-1)*LZ/2],
		 1, gauge_ty_edge_gath_split, g_nb_y_up, 198,
		 g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 2*LX*LZ],
		 1, gauge_ty_edge_cont      , g_nb_y_dn, 198,
		 g_cart_grid, &status);
#    endif    

    /* send the data to the neighbour on the left in y direction */
    /* recieve the data from the neighbour on the right in y direction */
    /* y2t edge */
#    ifndef _NEW_GEOMETRY
    MPI_Sendrecv(g_gauge_field[VOLUME + LZ],
		 1, gauge_ty_edge_gath, g_nb_y_dn, 297,
		 g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 4*LX*LZ],
		 1, gauge_ty_edge_cont, g_nb_y_up, 297, 
		 g_cart_grid, &status);
#    else
    MPI_Sendrecv(g_gauge_field[VOLUME + LZ/2],
		 1, gauge_ty_edge_gath_split, g_nb_y_dn, 297,
		 g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 4*LX*LZ],
		 1, gauge_ty_edge_cont      , g_nb_y_up, 297, 
		 g_cart_grid, &status);
#    endif
    
    /* send the data to the neighbour on the right in y direction */
    /* recieve the data from the neighbour on the left in y direction */
    /* y2t-edge */
#    ifndef _NEW_GEOMETRY
    MPI_Sendrecv(g_gauge_field[VOLUME + (LY-2)*LZ],
		 1, gauge_ty_edge_gath, g_nb_y_up, 298,
		 g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 6*LX*LZ],
		 1, gauge_ty_edge_cont, g_nb_y_dn, 298,
		 g_cart_grid, &status);
#    else
    MPI_Sendrecv(g_gauge_field[VOLUME + (LY-2)*LZ/2],
		 1, gauge_ty_edge_gath_split, g_nb_y_up, 298,
		 g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 6*LX*LZ],
		 1, gauge_ty_edge_cont      , g_nb_y_dn, 298,
		 g_cart_grid, &status);
#    endif
  }

  /* end of if defined PARALLELXYT || PARALLELXYZT */
#  endif
#  if defined PARALLELXYZT
  /* z-Rand */
  /* send the data to the neighbour on the left in z direction */
  /* recieve the data from the neighbour on the right in z direction */
#    ifndef _NEW_GEOMETRY
  MPI_Sendrecv(g_gauge_field[0],
	       1, gauge_z_slice_gath, g_nb_z_dn, 303,
	       g_gauge_field[VOLUME + 2*LZ*(LX*LY + T*LY) + 2*LZ*T*LX], 
	       1, gauge_z_slice_cont, g_nb_z_up, 303, 
	       g_cart_grid, &status);
#    else
  MPI_Sendrecv(g_gauge_field[0],
	       1, gauge_z_slice_gath_split, g_nb_z_dn, 303,
	       g_gauge_field[VOLUME + 2*LZ*(LX*LY + T*LY) + 2*LZ*T*LX],
	       1, gauge_z_slice_cont      , g_nb_z_up, 303, 
	       g_cart_grid, &status);
#    endif
  /* send the data to the neighbour on the right in z direction */
  /* recieve the data from the neighbour on the left in z direction */
#    ifndef _NEW_GEOMETRY
  MPI_Sendrecv(g_gauge_field[LZ-1],
	       1, gauge_z_slice_gath, g_nb_z_up, 304,
	       g_gauge_field[VOLUME + 2*LZ*(LX*LY + T*LY) + 2*T*LX*LZ + T*LX*LY],
	       1, gauge_z_slice_cont, g_nb_z_dn, 304,
	       g_cart_grid, &status);
#    else
  MPI_Sendrecv(g_gauge_field[(LZ-1)/2],
	       1, gauge_z_slice_gath_split, g_nb_z_up, 304,
	       g_gauge_field[VOLUME + 2*LZ*(LX*LY + T*LY) + 2*T*LX*LZ + T*LX*LY],
	       1, gauge_z_slice_cont      , g_nb_z_dn, 304,
	       g_cart_grid, &status);
#    endif

  if(g_dbw2rand > 0) {
    /* send the data to the neighbour on the left in z direction */
    /* recieve the data from the neighbour on the right in z direction */
    /* z2-Rand */
#    ifndef _NEW_GEOMETRY
    MPI_Sendrecv(g_gauge_field[1],
		 1, gauge_z_slice_gath, g_nb_z_dn, 305,
		 g_gauge_field[VOLUMEPLUSRAND+(2*LX+2*T)*LY*LZ + 2*T*LX*LZ],
		 1, gauge_z_slice_cont, g_nb_z_up, 305, 
		 g_cart_grid, &status);
#    else
    MPI_Sendrecv(g_gauge_field[0],
		 1, gauge_z_slice_gath_split, g_nb_z_dn, 305,
		 g_gauge_field[VOLUMEPLUSRAND+(2*LX+2*T)*LY*LZ + 2*T*LX*LZ],
		 1, gauge_z_slice_cont      , g_nb_z_up, 305, 
		 g_cart_grid, &status);
#    endif
    /* send the data to the neighbour on the right in z direction */
    /* recieve the data from the neighbour on the left in z direction */
    /* z2-Rand */
#    ifndef _NEW_GEOMETRY
    MPI_Sendrecv(g_gauge_field[LZ-2],
		 1, gauge_z_slice_gath, g_nb_z_up, 306,
		 g_gauge_field[VOLUMEPLUSRAND+(2*LX+2*T)*LY*LZ + 2*T*LX*LZ + T*LX*LY],
		 1, gauge_z_slice_cont, g_nb_z_dn, 306,
		 g_cart_grid, &status);
#    else
    MPI_Sendrecv(g_gauge_field[(LZ-2)/2],
		 1, gauge_z_slice_gath_split, g_nb_z_up, 306,
		 g_gauge_field[VOLUMEPLUSRAND+(2*LX+2*T)*LY*LZ + 2*T*LX*LZ + T*LX*LY],
		 1, gauge_z_slice_cont      , g_nb_z_dn, 306,
		 g_cart_grid, &status);
#    endif
  }

  /* edges */

  /* send the data to the neighbour on the left in x direction */
  /* recieve the data from the neighbour on the right in x direction */
  /* is on the z-Rand -> zx-edge*/
#    ifndef _NEW_GEOMETRY
  MPI_Sendrecv(g_gauge_field[VOLUME + 2*LX*LY*LZ + 2*T*LY*LZ + 2*T*LX*LZ], 
	       1, gauge_zx_edge_gath, g_nb_x_dn, 307,
	       g_gauge_field[VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 4*LX*LZ],
	       1, gauge_zx_edge_cont, g_nb_x_up, 307, 
	       g_cart_grid, &status);
#    else
  MPI_Sendrecv(g_gauge_field[VOLUME + 2*LX*LY*LZ + 2*T*LY*LZ + 2*T*LX*LZ], 
	       1, gauge_zx_edge_gath_split, g_nb_x_dn, 307,
	       g_gauge_field[VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 4*LX*LZ],
	       1, gauge_zx_edge_cont      , g_nb_x_up, 307, 
	       g_cart_grid, &status);
#    endif

  /* send the data to the neighbour on the right in x direction */
  /* recieve the data from the neighbour on the left in x direction */
  /* zx-edge */
#    ifndef _NEW_GEOMETRY
  MPI_Sendrecv(g_gauge_field[VOLUME + 2*LZ*(LX*LY + T*LY) + 2*T*LX*LZ + (LX-1)*LY], 
	       1, gauge_zx_edge_gath, g_nb_x_up, 308,
	       g_gauge_field[VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 4*LX*LZ + 2*T*LY],
	       1, gauge_zx_edge_cont, g_nb_x_dn, 308,
	       g_cart_grid, &status);
#    else
  MPI_Sendrecv(g_gauge_field[VOLUME + 2*LZ*(LX*LY + T*LY) + 2*T*LX*LZ + (LX-1)*LY/2], 
	       1, gauge_zx_edge_gath_split, g_nb_x_up, 308,
	       g_gauge_field[VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 4*LX*LZ + 2*T*LY],
	       1, gauge_zx_edge_cont      , g_nb_x_dn, 308,
	       g_cart_grid, &status);
#    endif

  /* send the data to the neighbour on the left in z direction */
  /* recieve the data from the neighbour on the right in z direction */
  /* is on the t-Rand -> tz-edge*/
#    ifndef _NEW_GEOMETRY
  MPI_Sendrecv(g_gauge_field[VOLUME],
	       1, gauge_tz_edge_gath, g_nb_z_dn, 309,
	       g_gauge_field[VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 4*LX*LZ + 4*T*LY], 
	       1, gauge_tz_edge_cont, g_nb_z_up, 309, 
	       g_cart_grid, &status);
#    else
  MPI_Sendrecv(g_gauge_field[VOLUME],
	       1, gauge_tz_edge_gath_split, g_nb_z_dn, 309,
	       g_gauge_field[VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 4*LX*LZ + 4*T*LY], 
	       1, gauge_tz_edge_cont      , g_nb_z_up, 309, 
	       g_cart_grid, &status);
#    endif

  /* send the data to the neighbour on the right in z direction */
  /* recieve the data from the neighbour on the left in z direction */
  /* tz-edge */
#    ifndef _NEW_GEOMETRY
  MPI_Sendrecv(g_gauge_field[VOLUME + (LZ-1)],
	       1, gauge_tz_edge_gath, g_nb_z_up, 310,
	       g_gauge_field[VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 4*LX*LZ + 4*T*LY + 2*LX*LY], 
	       1, gauge_tz_edge_cont, g_nb_z_dn, 310,
	       g_cart_grid, &status);
#    else
  MPI_Sendrecv(g_gauge_field[VOLUME + (LZ-1)/2],
	       1, gauge_tz_edge_gath_split, g_nb_z_up, 310,
	       g_gauge_field[VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 4*LX*LZ + 4*T*LY + 2*LX*LY], 
	       1, gauge_tz_edge_cont      , g_nb_z_dn, 110,
	       g_cart_grid, &status);
#    endif

  /* send the data to the neighbour on the left in y direction */
  /* recieve the data from the neighbour on the right in y direction */
  /* is on the z-Rand -> zy-edge*/
#    ifndef _NEW_GEOMETRY
  MPI_Sendrecv(g_gauge_field[VOLUME + 2*LZ*(LX*LY + T*LY) + 2*T*LX*LZ],
	       1, gauge_zy_edge_gath, g_nb_y_dn, 310,
	       g_gauge_field[VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 4*LX*LZ + 4*T*LY + 4*LX*LY], 
	       1, gauge_zy_edge_cont, g_nb_y_up, 310, 
	       g_cart_grid, &status);
#    else
  MPI_Sendrecv(g_gauge_field[VOLUME + 2*LZ*(LX*LY + T*LY) + 2*T*LX*LZ],
	       1, gauge_zy_edge_gath_split, g_nb_y_dn, 310,
	       g_gauge_field[VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 4*LX*LZ + 4*T*LY + 4*LX*LY], 
	       1, gauge_zy_edge_cont      , g_nb_y_up, 310, 
	       g_cart_grid, &status);
#    endif

  /* send the data to the neighbour on the right in y direction */
  /* recieve the data from the neighbour on the left in y direction */
  /* zy-edge */
#    ifndef _NEW_GEOMETRY
  MPI_Sendrecv(g_gauge_field[VOLUME + 2*LZ*(LX*LY + T*LY) + 2*T*LX*LZ + (LY-1)],
	       1, gauge_zy_edge_gath, g_nb_y_up, 310,
	       g_gauge_field[VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 4*LX*LZ + 4*T*LY + 4*LX*LY + 2*T*LX], 
	       1, gauge_zy_edge_cont, g_nb_y_dn, 310,
	       g_cart_grid, &status);
#    else
  MPI_Sendrecv(g_gauge_field[VOLUME + 2*LZ*(LX*LY + T*LY) + 2*T*LX*LZ + (LY-1)/2],
	       1, gauge_zy_edge_gath_split, g_nb_y_up, 310,
	       g_gauge_field[VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 4*LX*LZ + 4*T*LY + 2*LX*LY + 2*T*LX], 
	       1, gauge_zy_edge_cont      , g_nb_y_dn, 110,
	       g_cart_grid, &status);
#    endif


  /* end of if defined PARALLELXYZT */
#  endif
#endif
  return;
}

static char const rcsid[] = "$Id$";

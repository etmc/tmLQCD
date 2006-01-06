/* $Id$ */

/**********************************************************
 * 
 * exchange routines for gauge fields, spinor fields
 * and the derivatives
 *
 * Autor of the first version: Martin Hasenbusch
 *
 * extended by Carsten Urbach to 2 dimensions and
 * various versions of the geometry
 *
 **********************************************************/



#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include "global.h"
#include "mpi_init.h"
#include "sse.h"
#include "su3.h"
#include "su3adj.h"
#include "xchange.h"

#ifdef MPI
/* exchanges the field  l */
void xchange_field(spinor * const l) {

  MPI_Status status;
  /* send the data to the neighbour on the left */
  /* recieve the data from the neighbour on the right */
  MPI_Sendrecv((void*)l,              1, field_time_slice_cont, g_nb_t_dn, 81,
	       (void*)(l+T*LX*LY*LZ/2), 1, field_time_slice_cont, g_nb_t_up, 81,
	       g_cart_grid, &status);

  /* send the data to the neighbour on the right */
  /* recieve the data from the neighbour on the left */
  MPI_Sendrecv((void*)(l+(T-1)*LX*LY*LZ/2), 1, field_time_slice_cont, g_nb_t_up, 82,
	       (void*)(l+(T+1)*LX*LY*LZ/2),   1, field_time_slice_cont, g_nb_t_dn, 82,
	       g_cart_grid, &status);

#ifdef PARALLELXT
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


#endif

}

/********************************************************************
 * 
 * This version of the xchange routine seems to 
 * overcome a bug on the iwarp cluster: It seems to be
 * not possible to use field_x_slice_gath 
 * consistently on this pc cluster
 *
 ********************************************************************/

void xchange_field2(spinor * const l) {

#if (defined PARALLELXT) 
  int x0, x2, ix = 0, iz=0;
  spinor * buffer_x;
/*   static spinor buffer_x[T*LY*LZ]; */
#endif  

  MPI_Status status;
  /* send the data to the neighbour on the left */
  /* recieve the data from the neighbour on the right */
  MPI_Sendrecv((void*)l,              1, field_time_slice_cont, g_nb_t_dn, 81,
	       (void*)(l+T*LX*LY*LZ/2), 1, field_time_slice_cont, g_nb_t_up, 81,
	       g_cart_grid, &status);

  /* send the data to the neighbour on the right */
  /* recieve the data from the neighbour on the left */
  MPI_Sendrecv((void*)(l+(T-1)*LX*LY*LZ/2), 1, field_time_slice_cont, g_nb_t_up, 82,
	       (void*)(l+(T+1)*LX*LY*LZ/2), 1, field_time_slice_cont, g_nb_t_dn, 82,
	       g_cart_grid, &status);

#ifdef PARALLELXT
  iz = 0;
  for(x0 = 0; x0 < T; x0++) {
    for(x2 = 0; x2 < LY*LZ/2; x2++) {
      buffer_x[ix] = *(l+iz);
      ix++;
      iz++;
    }
    iz += (LX-1)*LY*LZ/2; 
  }

  iz = (LX-1)*LY*LZ/2;
  for(x0 = 0; x0 < T; x0++) {
    for(x2 = 0; x2 < LY*LZ/2; x2++) {
      buffer_x[ix] = *(l+iz);
      ix++;
      iz++;
    }
    iz += (LX-1)*LZ*LY/2;
  }

  /* send the data to the neighbour on the left in x direction */
  /* recieve the data from the neighbour on the right in x direction */
  MPI_Sendrecv(&buffer_x[0]               , 1, field_x_slice_cont, g_nb_x_dn, 91, 
 	       (void*)(l+(T+2)*LX*LY*LZ/2), 1, field_x_slice_cont, g_nb_x_up, 91,
	       g_cart_grid, &status);

  /* send the data to the neighbour on the right in x direction */
  /* recieve the data from the neighbour on the left in x direction */  
  MPI_Sendrecv(&buffer_x[T*LY*LZ/2]                   , 1, field_x_slice_cont, g_nb_x_up, 92,
	       (void*)(l+((T+2)*LX*LY*LZ + T*LY*LZ)/2), 1, field_x_slice_cont, g_nb_x_dn, 92,
	       g_cart_grid, &status);


#endif

}



void xchange_gauge() {

  MPI_Status status;
  /* send the data to the neighbour on the left */
  /* recieve the data from the neighbour on the right */
#ifndef _NEW_GEOMETRY
  MPI_Sendrecv(g_gauge_field[0],        1, gauge_time_slice_cont, g_nb_t_dn, 83, 
	       g_gauge_field[T*LX*LY*LZ], 1, gauge_time_slice_cont, g_nb_t_up, 83, 
	       g_cart_grid, &status);
#else
  MPI_Sendrecv(g_gauge_field[0],        1, gauge_time_slice_split, g_nb_t_dn, 83,
	       g_gauge_field[T*LX*LY*LZ], 1, gauge_time_slice_cont , g_nb_t_up, 83,
	       g_cart_grid, &status);
#endif
  /* send the data to the neighbour on the right */
  /* recieve the data from the neighbour on the left */
#ifndef _NEW_GEOMETRY
  MPI_Sendrecv(g_gauge_field[(T-1)*LX*LY*LZ], 1, gauge_time_slice_cont, g_nb_t_up, 84, 
	       g_gauge_field[(T+1)*LX*LY*LZ], 1, gauge_time_slice_cont, g_nb_t_dn, 84, 
	       g_cart_grid, &status);
#else
  MPI_Sendrecv(g_gauge_field[(T-1)*LX*LY*LZ/2], 1, gauge_time_slice_split, g_nb_t_up, 84, 
	       g_gauge_field[(T+1)*LX*LY*LZ],   1, gauge_time_slice_cont , g_nb_t_dn, 84,
	       g_cart_grid, &status);
#endif

  if(g_dbw2rand > 0) {
    /* send the data to the neighbour on the left */
    /* recieve the data from the neighbour on the right */
#ifndef _NEW_GEOMETRY
    MPI_Sendrecv(g_gauge_field[1*LX*LY*LZ],     1, gauge_time_slice_cont, g_nb_t_dn, 85, 
		 g_gauge_field[VOLUMEPLUSRAND], 1, gauge_time_slice_cont, g_nb_t_up, 85, 
		 g_cart_grid, &status);
#else
    MPI_Sendrecv(g_gauge_field[LX*LY*LZ/2],     1, gauge_time_slice_split, g_nb_t_dn, 85,
		 g_gauge_field[VOLUMEPLUSRAND], 1, gauge_time_slice_cont , g_nb_t_up, 85,
		 g_cart_grid, &status);
#endif
    /* send the data to the neighbour on the right */
    /* recieve the data from the neighbour on the left */
#ifndef _NEW_GEOMETRY
    MPI_Sendrecv(g_gauge_field[(T-2)*LX*LY*LZ],          1, gauge_time_slice_cont, g_nb_t_up, 86, 
		 g_gauge_field[VOLUMEPLUSRAND+LX*LY*LZ], 1, gauge_time_slice_cont, g_nb_t_dn, 86, 
		 g_cart_grid, &status);
#else
    MPI_Sendrecv(g_gauge_field[(T-2)*LX*LY*LZ/2],            1, gauge_time_slice_split, g_nb_t_up, 86, 
		 g_gauge_field[VOLUMEPLUSRAND+LX*LY*LZ], 1, gauge_time_slice_cont , g_nb_t_dn, 86,
		 g_cart_grid, &status);
#endif
  }
  
#ifdef PARALLELXT
  /* send the data to the neighbour on the left in x direction */
  /* recieve the data from the neighbour on the right in x direction */
#ifndef _NEW_GEOMETRY
  MPI_Sendrecv(g_gauge_field[0],              1, gauge_x_slice_gath, g_nb_x_dn, 93,
	       g_gauge_field[(T+2)*LX*LY*LZ], 1, gauge_x_slice_cont, g_nb_x_up, 93, 
	       g_cart_grid, &status);
#else
  MPI_Sendrecv(g_gauge_field[0],              1, gauge_x_slice_gath_split, g_nb_x_dn, 93,
	       g_gauge_field[(T+2)*LX*LY*LZ], 1, gauge_x_slice_cont      , g_nb_x_up, 93, 
	       g_cart_grid, &status);
#endif
  /* send the data to the neighbour on the right in x direction */
  /* recieve the data from the neighbour on the left in x direction */
#ifndef _NEW_GEOMETRY
  MPI_Sendrecv(g_gauge_field[(LX-1)*LY*LZ],             1, gauge_x_slice_gath, g_nb_x_up, 94,
	       g_gauge_field[(T+2)*LX*LY*LZ + T*LY*LZ], 1, gauge_x_slice_cont, g_nb_x_dn, 94,
	       g_cart_grid, &status);
#else
  MPI_Sendrecv(g_gauge_field[(LX-1)*LY*LZ/2],           1, gauge_x_slice_gath_split, g_nb_x_up, 94,
	       g_gauge_field[(T+2)*LX*LY*LZ + T*LY*LZ], 1, gauge_x_slice_cont      , g_nb_x_dn, 94,
	       g_cart_grid, &status);
#endif
  if(g_dbw2rand > 0) {
    /* send the data to the neighbour on the left in x direction */
    /* recieve the data from the neighbour on the right in x direction */
#ifndef _NEW_GEOMETRY
    MPI_Sendrecv(g_gauge_field[LY*LZ],                     1, gauge_x_slice_gath, g_nb_x_dn, 95,
		 g_gauge_field[VOLUMEPLUSRAND+2*LX*LY*LZ], 1, gauge_x_slice_cont, g_nb_x_up, 95, 
		 g_cart_grid, &status);
#else
    MPI_Sendrecv(g_gauge_field[LY*LZ/2],                   1, gauge_x_slice_gath_split, g_nb_x_dn, 95,
		 g_gauge_field[VOLUMEPLUSRAND+2*LX*LY*LZ], 1, gauge_x_slice_cont      , g_nb_x_up, 95, 
		 g_cart_grid, &status);
#endif
    /* send the data to the neighbour on the right in x direction */
    /* recieve the data from the neighbour on the left in x direction */
#ifndef _NEW_GEOMETRY
    MPI_Sendrecv(g_gauge_field[(LX-2)*LY*LZ],                        1, gauge_x_slice_gath, g_nb_x_up, 96,
		 g_gauge_field[VOLUMEPLUSRAND+2*LX*LY*LZ + T*LY*LZ], 1, gauge_x_slice_cont, g_nb_x_dn, 96,
		 g_cart_grid, &status);
#else
    MPI_Sendrecv(g_gauge_field[(LX-2)*LY*LZ/2],                      1, gauge_x_slice_gath_split, g_nb_x_up, 96,
		 g_gauge_field[VOLUMEPLUSRAND+2*LX*LY*LZ + T*LY*LZ], 1, gauge_x_slice_cont      , g_nb_x_dn, 96,
		 g_cart_grid, &status);
#endif
  }

  /* The edges */

  /* send the data to the neighbour on the left in t direction */
  /* recieve the data from the neighbour on the right in t direction */
#ifndef _NEW_GEOMETRY
  MPI_Sendrecv(&g_gauge_field[(T+2)*LX*LY*LZ][0], 1, gauge_yz_edge_gath, g_nb_t_dn, 95,
	       &g_gauge_field[VOLUME + RAND][0],  1, gauge_yz_edge_cont, g_nb_t_up, 95, 
	       g_cart_grid, &status);
#else
  MPI_Sendrecv(&g_gauge_field[(T+2)*LX*LY*LZ][0], 1, gauge_yz_edge_gath_split, g_nb_t_dn, 95,
	       &g_gauge_field[VOLUME + RAND][0],  1, gauge_yz_edge_cont      , g_nb_t_up, 95, 
	       g_cart_grid, &status);
#endif

  /* send the data to the neighbour on the right in t direction */
  /* recieve the data from the neighbour on the left in t direction */
#ifndef _NEW_GEOMETRY
  MPI_Sendrecv(&g_gauge_field[(T+2)*LX*LY*LZ + (T-1)*LY*LZ][0],  1, gauge_yz_edge_gath, g_nb_t_up, 96,
	       &g_gauge_field[VOLUME + RAND + 2*LY*LZ][0],       1, gauge_yz_edge_cont, g_nb_t_dn, 96,
	       g_cart_grid, &status);
#else
  MPI_Sendrecv(&g_gauge_field[(T+2)*LX*LY*LZ + (T-1)*LY*LZ/2][0],  1, gauge_yz_edge_gath_split, g_nb_t_up, 96,
	       &g_gauge_field[VOLUME + RAND + 2*LY*LZ][0],         1, gauge_yz_edge_cont      , g_nb_t_dn, 96,
	       g_cart_grid, &status);
#endif

  if(g_dbw2rand > 0) {
    /* send the data to the neighbour on the left in t direction */
    /* recieve the data from the neighbour on the right in t direction */
#ifndef _NEW_GEOMETRY
    MPI_Sendrecv(g_gauge_field[(T+2)*LX*LY*LZ + LY*LZ],
		 1, gauge_yz_edge_gath, g_nb_t_dn, 97,
		 g_gauge_field[VOLUMEPLUSRAND + 2*LX*LY*LZ + 2*T*LY*LZ],
		 1, gauge_yz_edge_cont, g_nb_t_up, 97, 
		 g_cart_grid, &status);
#else
    MPI_Sendrecv(g_gauge_field[(T+2)*LX*LY*LZ + LY*LZ/2],
		 1, gauge_yz_edge_gath_split, g_nb_t_dn, 97,
		 g_gauge_field[VOLUMEPLUSRAND + 2*LX*LY*LZ + 2*T*LY*LZ],
		 1, gauge_yz_edge_cont      , g_nb_t_up, 97, 
		 g_cart_grid, &status);
#endif
    
    /* send the data to the neighbour on the right in t direction */
    /* recieve the data from the neighbour on the left in t direction */
#ifndef _NEW_GEOMETRY
    MPI_Sendrecv(g_gauge_field[(T+2)*LX*LY*LZ + (T-2)*LY*LZ],
		 1, gauge_yz_edge_gath, g_nb_t_up, 98,
		 g_gauge_field[VOLUMEPLUSRAND + 2*LX*LY*LZ + 2*T*LY*LZ + 2*LY*LZ],
		 1, gauge_yz_edge_cont, g_nb_t_dn, 98,
		 g_cart_grid, &status);
#else
    MPI_Sendrecv(g_gauge_field[(T+2)*LX*LY*LZ + (T-2)*LY*LZ/2],
		 1, gauge_yz_edge_gath_split, g_nb_t_up, 98,
		 g_gauge_field[VOLUMEPLUSRAND + 2*LX*LY*LZ + 2*T*LY*LZ + 2*LY*LZ],
		 1, gauge_yz_edge_cont      , g_nb_t_dn, 98,
		 g_cart_grid, &status);
#endif    

    /* send the data to the neighbour on the left in t direction */
    /* recieve the data from the neighbour on the right in t direction */
#ifndef _NEW_GEOMETRY
    MPI_Sendrecv(g_gauge_field[VOLUMEPLUSRAND + 2*LX*LY*LZ],
		 1, gauge_yz_edge_gath, g_nb_t_dn, 97,
		 g_gauge_field[VOLUMEPLUSRAND + 2*LX*LY*LZ + 2*T*LY*LZ + 4*LY*LZ],
		 1, gauge_yz_edge_cont, g_nb_t_up, 97, 
		 g_cart_grid, &status);
#else
    MPI_Sendrecv(g_gauge_field[VOLUMEPLUSRAND + 2*LX*LY*LZ],
		 1, gauge_yz_edge_gath_split, g_nb_t_dn, 97,
		 g_gauge_field[VOLUMEPLUSRAND + 2*LX*LY*LZ + 2*T*LY*LZ + 4*LY*LZ],
		 1, gauge_yz_edge_cont      , g_nb_t_up, 97, 
		 g_cart_grid, &status);
#endif
    
    /* send the data to the neighbour on the right in t direction */
    /* recieve the data from the neighbour on the left in t direction */
#ifndef _NEW_GEOMETRY
    MPI_Sendrecv(g_gauge_field[VOLUMEPLUSRAND + 2*LX*LY*LZ + (T-1)*LY*LZ],
		 1, gauge_yz_edge_gath, g_nb_t_up, 98,
		 g_gauge_field[VOLUMEPLUSRAND + 2*LX*LY*LZ + 2*T*LY*LZ + 6*LY*LZ],
		 1, gauge_yz_edge_cont, g_nb_t_dn, 98,
		 g_cart_grid, &status);
#else
    MPI_Sendrecv(g_gauge_field[VOLUMEPLUSRAND + 2*LX*LY*LZ + (T-1)*LY*LZ/2],
		 1, gauge_yz_edge_gath_split, g_nb_t_up, 98,
		 g_gauge_field[VOLUMEPLUSRAND + 2*LX*LY*LZ + 2*T*LY*LZ + 6*LY*LZ],
		 1, gauge_yz_edge_cont      , g_nb_t_dn, 98,
		 g_cart_grid, &status);
#endif    
  }
  /* end of if defined PARALLELXT */
#endif
}

/********************************************************************
 * 
 * This version of the xchange routine seems to 
 * overcome a bug on the iwarp cluster: It seems to be
 * not possible to use field_x_slice_gath 
 * consistently on this pc cluster
 *
 ********************************************************************/

void xchange_gauge2() {
#if (defined PARALLELXT) 
  int x0, x2, ix = 0, iz=0, mu;
  su3 ** buffer_x=NULL, * buffer_x_=NULL; 

#endif
  MPI_Status status;
  /* send the data to the neighbour on the left */
  /* recieve the data from the neighbour on the right */
  MPI_Sendrecv(&g_gauge_field[0][0].c00.re,        1, gauge_time_slice_cont, g_nb_t_dn, 83, 
		    &g_gauge_field[T*LX*LY*LZ][0].c00.re, 1, gauge_time_slice_cont, g_nb_t_up, 83, 
		    g_cart_grid, &status);

  /* send the data to the neighbour on the right */
  /* recieve the data from the neighbour on the left */
  MPI_Sendrecv(&g_gauge_field[(T-1)*LX*LY*LZ][0].c00.re, 1, gauge_time_slice_cont, g_nb_t_up, 84, 
		    &g_gauge_field[(T+1)*LX*LY*LZ][0].c00.re, 1, gauge_time_slice_cont, g_nb_t_dn, 84, 
		    g_cart_grid, &status);

#ifdef PARALLELXT
  iz = 0;
  for(x0 = 0; x0 < T; x0++) {
    for(x2 = 0; x2 < LY*LZ; x2++) {
      for(mu = 0; mu < 4; mu++) {
	buffer_x[ix][mu] = g_gauge_field[ iz ][mu];
      }
      ix++;
      iz++;
    }
    iz += (LX-1)*LY*LZ; 
  }

  iz = (LX-1)*LY*LZ;
  for(x0 = 0; x0 < T; x0++) {
    for(x2 = 0; x2 < LY*LZ; x2++) {
      for(mu = 0; mu < 4; mu++) {
	buffer_x[ix][mu] = g_gauge_field[ iz ][mu];
      }
      ix++;
      iz++;
    }
    iz += (LX-1)*LZ*LY;
  }
  /* send the data to the neighbour on the left in x direction */
  /* recieve the data from the neighbour on the right in x direction */
  MPI_Sendrecv(&buffer_x[0][0],                   1, gauge_x_slice_cont, g_nb_x_dn, 93,
	       &g_gauge_field[(T+2)*LX*LY*LZ][0], 1, gauge_x_slice_cont, g_nb_x_up, 93, 
	       g_cart_grid, &status);

  /* send the data to the neighbour on the right in x direction */
  /* recieve the data from the neighbour on the left in x direction */
  MPI_Sendrecv(&buffer_x[T*LY*LZ][0],                            1, gauge_x_slice_cont, g_nb_x_up, 94,
		    &g_gauge_field[(T+2)*LX*LY*LZ + T*LY*LZ][0], 1, gauge_x_slice_cont, g_nb_x_dn, 94,
		    g_cart_grid, &status);

  /* The edges */

  /* send the data to the neighbour on the left in t direction */
  /* recieve the data from the neighbour on the right in t direction */
  MPI_Sendrecv(&g_gauge_field[(T+2)*LX*LY*LZ][0], 1, gauge_yz_edge_gath, g_nb_t_dn, 93,
		    &g_gauge_field[VOLUME + RAND][0],  1, gauge_yz_edge_cont, g_nb_t_up, 93, 
		    g_cart_grid, &status);

  /* send the data to the neighbour on the right in x direction */
  /* recieve the data from the neighbour on the left in x direction */
  MPI_Sendrecv(&g_gauge_field[(T+2)*LX*LY*LZ + (T-1)*LY*LZ][0],  1, gauge_yz_edge_gath, g_nb_t_up, 94,
		    &g_gauge_field[VOLUME + RAND + 2*LY*LZ][0],       1, gauge_yz_edge_cont, g_nb_t_dn, 94,
		    g_cart_grid, &status);
#endif

}

void xchange_deri()
{
  int ix,mu, t, y, z;
  MPI_Status status;
  /* send the data to the neighbour on the left in time direction */
  /* recieve the data from the neighbour on the right in time direction */
#ifndef _NEW_GEOMETRY
  MPI_Sendrecv(&df0[(T+1)*LX*LY*LZ][0].d1,        1, deri_time_slice_cont, g_nb_t_dn, 43,
	       &ddummy[(T-1)*LX*LY*LZ][0].d1,     1, deri_time_slice_cont, g_nb_t_up, 43,
	       g_cart_grid, &status);

  /* add ddummy to df0 */
  for(ix=(T-1)*LX*LY*LZ;ix < T*LX*LY*LZ; ix++){
    for(mu=0;mu<4;mu++){
      df0[ix][mu].d1 += ddummy[ix][mu].d1;
      df0[ix][mu].d2 += ddummy[ix][mu].d2;
      df0[ix][mu].d3 += ddummy[ix][mu].d3;
      df0[ix][mu].d4 += ddummy[ix][mu].d4;
      df0[ix][mu].d5 += ddummy[ix][mu].d5;
      df0[ix][mu].d6 += ddummy[ix][mu].d6;
      df0[ix][mu].d7 += ddummy[ix][mu].d7;
      df0[ix][mu].d8 += ddummy[ix][mu].d8;
    }
  }
#else
  MPI_Sendrecv(&df0[(T+1)*LX*LY*LZ][0].d1,          1, deri_time_slice_cont , g_nb_t_dn, 43, 
	       &ddummy[(T-1)*LX*LY*LZ/2][0].d1,     1, deri_time_slice_split, g_nb_t_up, 43,
	       g_cart_grid, &status);
  /* add ddummy to df0 */
  for(ix=(T-1)*LX*LY*LZ/2;ix < T*LX*LY*LZ/2; ix++){
    for(mu=0;mu<4;mu++){
      df0[ix][mu].d1 += ddummy[ix][mu].d1;
      df0[ix][mu].d2 += ddummy[ix][mu].d2;
      df0[ix][mu].d3 += ddummy[ix][mu].d3;
      df0[ix][mu].d4 += ddummy[ix][mu].d4;
      df0[ix][mu].d5 += ddummy[ix][mu].d5;
      df0[ix][mu].d6 += ddummy[ix][mu].d6;
      df0[ix][mu].d7 += ddummy[ix][mu].d7;
      df0[ix][mu].d8 += ddummy[ix][mu].d8;
    }
  }
  for(ix=VOLUME-LX*LY*LZ/2;ix < VOLUME; ix++){
    for(mu=0;mu<4;mu++){
      df0[ix][mu].d1 += ddummy[ix][mu].d1;
      df0[ix][mu].d2 += ddummy[ix][mu].d2;
      df0[ix][mu].d3 += ddummy[ix][mu].d3;
      df0[ix][mu].d4 += ddummy[ix][mu].d4;
      df0[ix][mu].d5 += ddummy[ix][mu].d5;
      df0[ix][mu].d6 += ddummy[ix][mu].d6;
      df0[ix][mu].d7 += ddummy[ix][mu].d7;
      df0[ix][mu].d8 += ddummy[ix][mu].d8;
    }
  }
#endif
  /* send the data to the neighbour on the right is not needed*/

#ifdef PARALLELXT

  /* send the data to the neighbour on the left in x direction */
  /* recieve the data from the neighbour on the right in x direction */
#ifndef _NEW_GEOMETRY
  MPI_Sendrecv(&df0[(T+2)*LX*LY*LZ + T*LY*LZ][0],    1, deri_x_slice_cont, g_nb_x_dn, 44,
	       &ddummy[(LX-1)*LY*LZ][0],             1, deri_x_slice_gath, g_nb_x_up, 44,
	       g_cart_grid, &status);
#else
  MPI_Sendrecv(&df0[(T+2)*LX*LY*LZ + T*LY*LZ][0],    1, deri_x_slice_cont      , g_nb_x_dn, 44,
	       &ddummy[(LX-1)*LY*LZ/2][0],           1, deri_x_slice_gath_split, g_nb_x_up, 44,
	       g_cart_grid, &status);
#endif
  /* add ddummy to df0 */
  for(t = 0; t < T; t++) {
    for(y = 0; y < LY; y++) {
      for(z = 0; z < LZ; z++) {
	ix = g_ipt[t][LX-1][y][z];
	for(mu=0;mu<4;mu++){
	  df0[ix][mu].d1 += ddummy[ix][mu].d1;
	  df0[ix][mu].d2 += ddummy[ix][mu].d2;
	  df0[ix][mu].d3 += ddummy[ix][mu].d3;
	  df0[ix][mu].d4 += ddummy[ix][mu].d4;
	  df0[ix][mu].d5 += ddummy[ix][mu].d5;
	  df0[ix][mu].d6 += ddummy[ix][mu].d6;
	  df0[ix][mu].d7 += ddummy[ix][mu].d7;
	  df0[ix][mu].d8 += ddummy[ix][mu].d8;
	}
      }
    }
  }
  /* send the data to the neighbour on the right is not needed*/  

#endif


  /* For dclover */

  if(g_use_clover_flag == 1) {
    
    /* send the data to the neighbour on the right */
    /* recieve the data from the neighbour on the left */
    MPI_Sendrecv(&dclover[T*LX*LY*LZ][0].d1, 1, deri_time_slice_cont, g_nb_t_up, 53,
		 &ddummy[0][0].d1,         1, deri_time_slice_cont, g_nb_t_dn, 53,
		 g_cart_grid, &status);
    /* add ddummy to dclover */
    
    for(ix=0;ix < LX*LY*LZ; ix++){
      for(mu=0;mu<4;mu++){
	dclover[ix][mu].d1 += ddummy[ix][mu].d1;
	dclover[ix][mu].d2 += ddummy[ix][mu].d2;
	dclover[ix][mu].d3 += ddummy[ix][mu].d3;
	dclover[ix][mu].d4 += ddummy[ix][mu].d4;
	dclover[ix][mu].d5 += ddummy[ix][mu].d5;
	dclover[ix][mu].d6 += ddummy[ix][mu].d6;
	dclover[ix][mu].d7 += ddummy[ix][mu].d7;
	dclover[ix][mu].d8 += ddummy[ix][mu].d8;
      }
    }
    
    /* send the data to the neighbour on the left */
    /* recieve the data from the neighbour on the right */
    MPI_Sendrecv(&dclover[(T+1)*LX*LY*LZ][0].d1, 1, deri_time_slice_cont, g_nb_t_dn, 54,
		 &ddummy[(T-1)*LX*LY*LZ][0].d1,  1, deri_time_slice_cont, g_nb_t_up, 54,
		 g_cart_grid, &status);
    /* add ddummy to dclover */
    
    for(ix=(T-1)*LX*LY*LZ; ix < T*LX*LY*LZ; ix++){
      for(mu=0;mu<4;mu++){
	dclover[ix][mu].d1 += ddummy[ix][mu].d1;
	dclover[ix][mu].d2 += ddummy[ix][mu].d2;
	dclover[ix][mu].d3 += ddummy[ix][mu].d3;
	dclover[ix][mu].d4 += ddummy[ix][mu].d4;
	dclover[ix][mu].d5 += ddummy[ix][mu].d5;
	dclover[ix][mu].d6 += ddummy[ix][mu].d6;
	dclover[ix][mu].d7 += ddummy[ix][mu].d7;
	dclover[ix][mu].d8 += ddummy[ix][mu].d8;
      }
    }
  }
}



#ifdef OlD

void xchange_gauge()
{
#ifdef MPI
  MPI_Status status;
  /* send the data to the neighbour on the left */
  /* recieve the data from the neighbour on the right */
  MPI_Sendrecv(&g_gauge_field[0][0].c00.re,               1, gauge_time_slice_split, g_nb_t_dn, 83,
	       &g_gauge_field[T*LX*LY*LZ][0].c00.re,        1, gauge_time_slice_cont , g_nb_t_up, 83,
	       g_cart_grid, &status);
  
  /* send the data to the neighbour on the right */
  /* recieve the data from the neighbour on the left */
  MPI_Sendrecv(&g_gauge_field[(T-1)*LX*LY*LX/2][0].c00.re,      1, gauge_time_slice_split, g_nb_t_up, 84, 
	       &g_gauge_field[(T+1)*LX*LY*LZ][0].c00.re,        1, gauge_time_slice_cont , g_nb_t_dn, 84,
	       g_cart_grid, &status);

#endif
}


void xchange_deri()
{
#ifdef MPI
  int ix,mu;
  MPI_Status status;
  /* send the data to the neighbour on the left */
  /* recieve the data from the neighbour on the right */
  MPI_Sendrecv(&df0[(T+1)*LX*LY*LZ][0].d1,          1, deri_time_slice_cont , g_nb_t_dn, 43,
	       &ddummy[(T-1)*LX*LY*LZ/2][0].d1,     1, deri_time_slice_split, g_nb_t_up, 43,
	       g_cart_grid, &status);
  /* add ddummy to df0 */
  for(ix=(T-1)*LX*LY*LZ/2;ix < T*LX*LY*LZ/2; ix++){
    for(mu=0;mu<4;mu++){
      df0[ix][mu].d1 += ddummy[ix][mu].d1;
      df0[ix][mu].d2 += ddummy[ix][mu].d2;
      df0[ix][mu].d3 += ddummy[ix][mu].d3;
      df0[ix][mu].d4 += ddummy[ix][mu].d4;
      df0[ix][mu].d5 += ddummy[ix][mu].d5;
      df0[ix][mu].d6 += ddummy[ix][mu].d6;
      df0[ix][mu].d7 += ddummy[ix][mu].d7;
      df0[ix][mu].d8 += ddummy[ix][mu].d8;
    }
  }
  for(ix=VOLUME-LX*LY*LZ/2;ix < VOLUME; ix++){
    for(mu=0;mu<4;mu++){
      df0[ix][mu].d1 += ddummy[ix][mu].d1;
      df0[ix][mu].d2 += ddummy[ix][mu].d2;
      df0[ix][mu].d3 += ddummy[ix][mu].d3;
      df0[ix][mu].d4 += ddummy[ix][mu].d4;
      df0[ix][mu].d5 += ddummy[ix][mu].d5;
      df0[ix][mu].d6 += ddummy[ix][mu].d6;
      df0[ix][mu].d7 += ddummy[ix][mu].d7;
      df0[ix][mu].d8 += ddummy[ix][mu].d8;
    }
  }
  /* send the data to the neighbour on the right is not needed*/

#endif
}

#endif


static char const rcsid[] = "$Id$";

#endif

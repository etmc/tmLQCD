/* $Id$ */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#ifdef MPI
#include <mpi.h>
#endif
#include "sse.h"
#include "su3.h"
#include "su3adj.h"
#include "global.h"
#include "xchange.h"

#define OlD 

/* exchanges the field  l */
void xchange_field(int l)
{
#ifdef MPI
  MPI_Status status;
  /* send the data to the neighbour on the left */
  /* recieve the data from the neighbour on the right */
  MPI_Sendrecv(&spinor_field[l][0].c1.c1.re,         1, g_field_time_slice_cont, g_nb_t_dn, 81,
	       &spinor_field[l][T*LX*L*L/2].c1.c1.re, 1, g_field_time_slice_cont, g_nb_t_up, 81,
	       g_cart_grid, &status);

  /* send the data to the neighbour on the right */
  /* recieve the data from the neighbour on the left */
  MPI_Sendrecv(&spinor_field[l][(T-1)*LX*L*L/2].c1.c1.re, 1, g_field_time_slice_cont, g_nb_t_up, 82,
	       &spinor_field[l][(T+1)*LX*L*L/2].c1.c1.re, 1, g_field_time_slice_cont, g_nb_t_dn, 82,
	       g_cart_grid, &status);

#ifdef PARALLELXT
  /* send the data to the neighbour on the left in x direction */
  /* recieve the data from the neighbour on the right in x direction */
  MPI_Sendrecv(&spinor_field[l][0],                1, g_field_x_slice_gath, g_nb_x_dn, 91,
	       &spinor_field[l][(T+2)*LX*LY*LZ/2], 1, g_field_x_slice_cont, g_nb_x_up, 91,
	       g_cart_grid, &status);

  /* send the data to the neighbour on the right in x direction */
  /* recieve the data from the neighbour on the left in x direction */  
  MPI_Sendrecv(&spinor_field[l][LX/2-1],                       1, g_field_x_slice_gath, g_nb_x_up, 92,
	       &spinor_field[l][((T+2)*LX*LY*LZ + T*LY*LZ)/2], 1, g_field_x_slice_cont, g_nb_x_dn, 92,
	       g_cart_grid, &status);


#endif

#endif
}


#ifdef OlD

void xchange_gauge()
{
#ifdef MPI

  MPI_Status status;
  /* send the data to the neighbour on the left */
  /* recieve the data from the neighbour on the right */
  MPI_Sendrecv(&g_gauge_field[0][0].c11.re,        1, g_gauge_time_slice_cont, g_nb_t_dn, 83, 
 	       &g_gauge_field[T*LX*L*L][0].c11.re, 1, g_gauge_time_slice_cont, g_nb_t_up, 83, 
	       g_cart_grid, &status);

  /* send the data to the neighbour on the right */
  /* recieve the data from the neighbour on the left */
  MPI_Sendrecv(&g_gauge_field[(T-1)*LX*L*L][0].c11.re, 1, g_gauge_time_slice_cont, g_nb_t_up, 84, 
 	       &g_gauge_field[(T+1)*LX*L*L][0].c11.re, 1, g_gauge_time_slice_cont, g_nb_t_dn, 84, 
	       g_cart_grid, &status);

#ifdef PARALLELXT
  /* send the data to the neighbour on the left in x direction */
  /* recieve the data from the neighbour on the right in x direction */
  MPI_Sendrecv(&g_gauge_field[0][0],              1, g_gauge_x_slice_gath, g_nb_x_dn, 93,
	       &g_gauge_field[(T+2)*LX*LY*LZ][0], 1, g_gauge_x_slice_cont, g_nb_x_up, 93,
	       g_cart_grid, &status);

  /* send the data to the neighbour on the right in x direction */
  /* recieve the data from the neighbour on the left in x direction */
  MPI_Sendrecv(&g_gauge_field[LX-1][0],                     1, g_gauge_x_slice_gath, g_nb_x_up, 94,
	       &g_gauge_field[(T+2)*LX*LY*LZ + T*LY*LZ][0], 1, g_gauge_x_slice_cont, g_nb_x_dn, 94,
	       g_cart_grid, &status);
#endif

#endif
}

void xchange_deri()
{
#ifdef MPI
  int ix,mu, t, y, z;
  MPI_Status status;
  /* send the data to the neighbour on the left in time direction */
  /* recieve the data from the neighbour on the right in time direction */
  MPI_Sendrecv(&df0[(T+1)*LX*L*L][0].d1,        1, g_deri_time_slice_cont, g_nb_t_dn, 43,
	       &ddummy[(T-1)*LX*L*L][0].d1,     1, g_deri_time_slice_cont, g_nb_t_up, 43,
	       g_cart_grid, &status);
  /* add ddummy to df0 */
  for(ix=(T-1)*LX*L*L;ix < T*LX*L*L; ix++){
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

#ifdef PARALLELXT

  /* send the data to the neighbour on the left in x direction */
  /* recieve the data from the neighbour on the right in x direction */
  MPI_Sendrecv(&df0[(T+2)*LX*LY*LZ + T*LY*LZ][0],    1, g_deri_x_slice_cont, g_nb_x_dn, 44,
	       &ddummy[LX-1][0],                     1, g_deri_x_slice_gath, g_nb_x_up, 44,
	       g_cart_grid, &status);
  /* add ddummy to df0 */
  for(t = 0; t < 0; t++) {
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
    MPI_Sendrecv(&dclover[T*LX*L*L][0].d1, 1, g_deri_time_slice_cont, g_nb_t_up, 53,
		 &ddummy[0][0].d1,         1, g_deri_time_slice_cont, g_nb_t_dn, 53,
		 g_cart_grid, &status);
    /* add ddummy to dclover */
    
    for(ix=0;ix < LX*L*L; ix++){
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
    MPI_Sendrecv(&dclover[(T+1)*LX*L*L][0].d1, 1, g_deri_time_slice_cont, g_nb_t_dn, 54,
		 &ddummy[(T-1)*LX*L*L][0].d1,  1, g_deri_time_slice_cont, g_nb_t_up, 54,
		 g_cart_grid, &status);
    /* add ddummy to dclover */
    
    for(ix=(T-1)*LX*L*L; ix < T*LX*L*L; ix++){
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
#endif
}


#endif


#ifndef OlD

void xchange_gauge()
{
#ifdef MPI
  MPI_Status status;
  /* send the data to the neighbour on the left */
  /* recieve the data from the neighbour on the right */
  MPI_Sendrecv(&g_gauge_field[0][0].c11.re,               1, g_gauge_time_slice_split, g_nb_t_dn, 83,
	       &g_gauge_field[T*LX*L*L][0].c11.re,        1, g_gauge_time_slice_cont , g_nb_t_up, 83,
	       g_cart_grid, &status);
  
  /* send the data to the neighbour on the right */
  /* recieve the data from the neighbour on the left */
  MPI_Sendrecv(&g_gauge_field[(T-1)*LX*L*L/2][0].c11.re,      1, g_gauge_time_slice_split, g_nb_t_up, 84, 
	       &g_gauge_field[(T+1)*LX*L*L][0].c11.re,        1, g_gauge_time_slice_cont , g_nb_t_dn, 84,
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
  MPI_Sendrecv(&df0[(T+1)*LX*L*L][0].d1,          1, g_deri_time_slice_cont , g_nb_t_dn, 43,
	       &ddummy[(T-1)*LX*L*L/2][0].d1,     1, g_deri_time_slice_split, g_nb_t_up, 43,
	       g_cart_grid, &status);
  /* add ddummy to df0 */
  for(ix=(T-1)*LX*L*L/2;ix < T*LX*L*L/2; ix++){
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
  for(ix=VOLUME-L*L*L/2;ix < VOLUME; ix++){
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





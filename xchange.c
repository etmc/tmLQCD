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

/* exchanges the field  l */
void xchange_field(int l)
{
#ifdef MPI
  int i,k;
  MPI_Status status;
  /* send the data to the neighbour on the left */
  i=g_proc_id-1; 
  if(g_proc_id==0){
    i=g_nproc-1;
  }
  /* recieve the data from the neighbour on the right */
  k=g_proc_id+1; 
  if(g_proc_id==g_nproc-1){
    k=0;
  }
  MPI_Sendrecv(&spinor_field[l][0].c1.c1.re,       12*L*L*L, MPI_DOUBLE, i, 81,
	       &spinor_field[l][T*L*L*L/2].c1.c1.re, 12*L*L*L, MPI_DOUBLE, k, 81,
	       MPI_COMM_WORLD, &status);

  /* send the data to the neighbour on the right */
  i=g_proc_id+1; 
  if(g_proc_id==g_nproc-1){
    i=0;
  }
  /* recieve the data from the neighbour on the left */
  k=g_proc_id-1; 
  if(g_proc_id==0){
    k=g_nproc-1;
  }
  MPI_Sendrecv(&spinor_field[l][(T-1)*L*L*L/2].c1.c1.re, 12*L*L*L, MPI_DOUBLE, i, 82,
	       &spinor_field[l][(T+1)*L*L*L/2].c1.c1.re, 12*L*L*L, MPI_DOUBLE, k, 82,
	       MPI_COMM_WORLD, &status);
#endif
}

void xchange_gauge()
{
#ifdef MPI
  int i,k;
  MPI_Status status;
  /* send the data to the neighbour on the left */
  i=g_proc_id-1; 
  if(g_proc_id==0){
    i=g_nproc-1;
  }
  /* recieve the data from the neighbour on the right */
  k=g_proc_id+1; 
  if(g_proc_id==g_nproc-1){
    k=0;
  }
  MPI_Sendrecv(&g_gauge_field[0][0].c11.re,       72*L*L*L, MPI_DOUBLE, i, 83,
	       &g_gauge_field[T*L*L*L][0].c11.re, 72*L*L*L, MPI_DOUBLE, k, 83,
	       MPI_COMM_WORLD, &status);

  /* send the data to the neighbour on the right */
  i=g_proc_id+1; 
  if(g_proc_id==g_nproc-1){
    i=0;
  }
  /* recieve the data from the neighbour on the left */
  k=g_proc_id-1; 
  if(g_proc_id==0){
    k=g_nproc-1;
  }
  MPI_Sendrecv(&g_gauge_field[(T-1)*L*L*L][0].c11.re, 72*L*L*L, MPI_DOUBLE, i, 84,
	       &g_gauge_field[(T+1)*L*L*L][0].c11.re, 72*L*L*L, MPI_DOUBLE, k, 84,
	       MPI_COMM_WORLD, &status);
#endif
}

void xchange_deri()
{
#ifdef MPI
  int i,k;
  int ix,mu;
  MPI_Status status;
  /* send the data to the neighbour on the left */
  i=g_proc_id-1; 
  if(g_proc_id==0){
    i=g_nproc-1;
  }
  /* recieve the data from the neighbour on the right */
  k=g_proc_id+1; 
  if(g_proc_id==g_nproc-1){
    k=0;
  }
  MPI_Sendrecv(&df0[(T+1)*L*L*L][0].d1,        32*L*L*L, MPI_DOUBLE, i, 43,
	       &ddummy[(T-1)*L*L*L][0].d1,     32*L*L*L, MPI_DOUBLE, k, 43,
	       MPI_COMM_WORLD, &status);
  /* add ddummy to df0 */
  for(ix=(T-1)*L*L*L;ix < T*L*L*L; ix++){
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

  /* For dclover */

  /* send the data to the neighbour on the right */
  i=g_proc_id+1; 
  if(g_proc_id==g_nproc-1){
    i=0;
  }
  /* recieve the data from the neighbour on the left */
  k=g_proc_id-1; 
  if(g_proc_id==0){
    k=g_nproc-1;
  }
  MPI_Sendrecv(&dclover[T*L*L*L][0].d1, 32*L*L*L, MPI_DOUBLE, i, 53,
	       &ddummy[0][0].d1,        32*L*L*L, MPI_DOUBLE, k, 53,
	       MPI_COMM_WORLD, &status);
  /* add ddummy to dclover */

  for(ix=0;ix < L*L*L; ix++){
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
  i=g_proc_id-1; 
  if(g_proc_id==0){
    i=g_nproc-1;
  }
  /* recieve the data from the neighbour on the right */
  k=g_proc_id+1; 
  if(g_proc_id==g_nproc-1){
    k=0;
  }
  MPI_Sendrecv(&dclover[(T+1)*L*L*L][0].d1, 32*L*L*L, MPI_DOUBLE, i, 54,
	       &ddummy[(T-1)*L*L*L][0].d1,  32*L*L*L, MPI_DOUBLE, k, 54,
	       MPI_COMM_WORLD, &status);
  /* add ddummy to dclover */

  for(ix=(T-1)*L*L*L; ix < T*L*L*L; ix++){
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
#endif
}

/* $Id$ */
#ifndef _MPI_INIT_H
#define _MPI_INIT_H

#ifdef MPI
#include <mpi.h>


/* Datatypes for the data exchange */
extern MPI_Datatype gauge_time_slice_cont;
extern MPI_Datatype gauge_time_slice_split;
extern MPI_Datatype deri_time_slice_cont;
extern MPI_Datatype deri_time_slice_split;
extern MPI_Datatype field_time_slice_cont;
extern MPI_Datatype gauge_x_slice_cont;
extern MPI_Datatype gauge_x_slice_gath;
extern MPI_Datatype gauge_x_slice_gath_split;
extern MPI_Datatype field_x_slice_cont;
extern MPI_Datatype field_x_slice_gath;
extern MPI_Datatype deri_x_slice_cont;
extern MPI_Datatype deri_x_slice_gath;
extern MPI_Datatype deri_x_slice_gath_split;
extern MPI_Datatype gauge_xy_edge_cont;
extern MPI_Datatype gauge_xy_edge_gath;
extern MPI_Datatype gauge_xy_edge_gath_split;
#endif

void mpi_init(int argc, char *argv[]);

#endif

/* $Id$ */
#ifndef _MPI_INIT_H
#define _MPI_INIT_H

#include <mpi.h>

/* Datatypes for the data exchange */
extern MPI_Datatype gauge_point, gauge_time_slice_cont;
extern MPI_Datatype gauge_time_slice_split;
extern MPI_Datatype deri_point, deri_time_slice_cont;
extern MPI_Datatype deri_time_slice_split;
extern MPI_Datatype field_point, field_time_slice_cont;
extern MPI_Datatype gauge_x_slice_cont;
extern MPI_Datatype gauge_yz_subslice;
extern MPI_Datatype gauge_x_slice_gath;
extern MPI_Datatype field_x_slice_cont;
extern MPI_Datatype field_yz_subslice;
extern MPI_Datatype field_x_slice_gath;
extern MPI_Datatype deri_x_slice_cont;
extern MPI_Datatype deri_yz_subslice;
extern MPI_Datatype deri_x_slice_gath;
extern MPI_Datatype gauge_xy_edge_cont;
extern MPI_Datatype gauge_xy_edge_gath;

void mpi_init(int argc,char *argv[]);

#endif

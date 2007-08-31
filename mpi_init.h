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
extern MPI_Datatype lfield_time_slice_cont;
extern MPI_Datatype gauge_x_slice_cont;
extern MPI_Datatype gauge_x_slice_gath;
extern MPI_Datatype field_x_slice_cont;
extern MPI_Datatype field_x_slice_gath;
extern MPI_Datatype lfield_x_slice_cont;
extern MPI_Datatype lfield_x_slice_gath;
extern MPI_Datatype deri_x_slice_cont;
extern MPI_Datatype deri_x_slice_gath;
extern MPI_Datatype gauge_xt_edge_cont;
extern MPI_Datatype gauge_xt_edge_gath;

extern MPI_Datatype gauge_yx_edge_cont;
extern MPI_Datatype gauge_yx_edge_gath;

extern MPI_Datatype gauge_ty_edge_cont;
extern MPI_Datatype gauge_ty_edge_gath;

extern MPI_Datatype gauge_zx_edge_cont;
extern MPI_Datatype gauge_zx_edge_gath;

extern MPI_Datatype gauge_tz_edge_cont;
extern MPI_Datatype gauge_tz_edge_gath;

extern MPI_Datatype gauge_zy_edge_cont;
extern MPI_Datatype gauge_zy_edge_gath;

extern MPI_Datatype gauge_y_slice_cont;
extern MPI_Datatype gauge_y_slice_gath;
extern MPI_Datatype field_y_slice_cont;
extern MPI_Datatype field_y_slice_gath;
extern MPI_Datatype lfield_y_slice_cont;
extern MPI_Datatype lfield_y_slice_gath;
extern MPI_Datatype deri_y_slice_cont;
extern MPI_Datatype deri_y_slice_gath;

extern MPI_Datatype deri_z_slice_cont;
extern MPI_Datatype deri_z_slice_gath;

extern MPI_Datatype gauge_z_slice_gath;
extern MPI_Datatype gauge_z_slice_cont;

extern MPI_Datatype field_z_slice_cont;
extern MPI_Datatype lfield_z_slice_cont;
extern MPI_Datatype lfield_z_slice_gath;
extern MPI_Datatype field_z_slice_half;

extern MPI_Datatype halffield_point;
extern MPI_Datatype halffield_time_slice_cont;
extern MPI_Datatype halffield_x_slice_cont;
extern MPI_Datatype halffield_x_slice_gath;
extern MPI_Datatype halffield_y_slice_cont;
extern MPI_Datatype halffield_y_slice_gath;
extern MPI_Datatype halffield_z_slice_cont;

#ifdef PARALLELXYZT
extern spinor * field_buffer_z ALIGN;
extern spinor * field_buffer_z2 ALIGN;
extern halfspinor * halffield_buffer_z ALIGN;
extern halfspinor * halffield_buffer_z2 ALIGN;
#endif

extern MPI_Comm mpi_time_slices;
#endif

extern int mpi_time_rank;

void mpi_init(int argc, char *argv[]);

#endif

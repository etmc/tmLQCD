/***********************************************************************
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008 Carsten Urbach
 *
 * This file is part of tmLQCD.
 *
 * tmLQCD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * tmLQCD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with tmLQCD.  If not, see <http://www.gnu.org/licenses/>.
 ***********************************************************************/
#ifndef _MPI_INIT_H
#define _MPI_INIT_H

#ifdef MPI
#include <mpi.h>


/* Datatypes for the data exchange */
extern MPI_Datatype mpi_su3;
extern MPI_Datatype field_point;
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
extern MPI_Datatype deri_xt_edge_cont;
extern MPI_Datatype deri_xt_edge_gath;

extern MPI_Datatype gauge_yx_edge_cont;
extern MPI_Datatype gauge_yx_edge_gath;
extern MPI_Datatype deri_yx_edge_cont;

extern MPI_Datatype gauge_ty_edge_cont;
extern MPI_Datatype gauge_ty_edge_gath;
extern MPI_Datatype deri_ty_edge_cont;

extern MPI_Datatype gauge_zx_edge_cont;
extern MPI_Datatype gauge_zx_edge_gath;
extern MPI_Datatype deri_zx_edge_cont;

extern MPI_Datatype gauge_tz_edge_cont;
extern MPI_Datatype gauge_tz_edge_gath;
extern MPI_Datatype deri_tz_edge_cont;

extern MPI_Datatype gauge_zy_edge_cont;
extern MPI_Datatype gauge_zy_edge_gath;
extern MPI_Datatype deri_zy_edge_cont;

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
extern MPI_Datatype field_z_slice_gath;
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

#ifdef _USE_TSPLITPAR
extern MPI_Datatype field_xt_slice_int;
extern MPI_Datatype field_yt_slice_int;
extern MPI_Datatype field_xt_slice_ext;
extern MPI_Datatype field_yt_slice_ext;
# ifdef PARALLELXYZ
extern MPI_Datatype field_zt_slice_ext_L;
extern MPI_Datatype field_zt_slice_ext_S;
extern MPI_Datatype field_zt_slice_even_dn_et;
extern MPI_Datatype field_zt_slice_even_up_et;
extern MPI_Datatype field_zt_slice_odd_dn_et;
extern MPI_Datatype field_zt_slice_odd_up_et;
extern MPI_Datatype field_zt_slice_even_dn_ot;
extern MPI_Datatype field_zt_slice_even_up_ot;
extern MPI_Datatype field_zt_slice_odd_dn_ot;
extern MPI_Datatype field_zt_slice_odd_up_ot;
# endif
#endif
#ifdef WITHLAPH
extern MPI_Datatype su3vect_point;
extern MPI_Datatype jfield_x_slice_cont;
extern MPI_Datatype jfield_y_slice_cont;
extern MPI_Datatype jfield_z_slice_cont;
extern MPI_Datatype jfield_x_slice_gath;
extern MPI_Datatype jfield_y_slice_gath;
extern MPI_Datatype jfield_z_slice_gath;
extern MPI_Datatype jfield_y_subslice;
#endif

#if ( defined PARALLELXYZT || defined PARALLELXYZ )
extern MPI_Datatype field_z_slice_even_dn;
extern MPI_Datatype field_z_slice_even_up;
extern MPI_Datatype field_z_slice_odd_dn;
extern MPI_Datatype field_z_slice_odd_up;

# if (!defined _INDEX_INDEP_GEOM)
extern spinor * field_buffer_z ALIGN;
extern spinor * field_buffer_z2 ALIGN;
extern spinor * field_buffer_z3 ALIGN;
extern spinor * field_buffer_z4 ALIGN;
extern halfspinor * halffield_buffer_z ALIGN;
extern halfspinor * halffield_buffer_z2 ALIGN;
# endif
#endif

extern MPI_Op mpi_reduce_su3_ray;
void reduce_su3_ray(void *u_i, void *u_io, int *len, MPI_Datatype *dt);

#endif

void tmlqcd_mpi_init(int argc, char *argv[]);

#endif

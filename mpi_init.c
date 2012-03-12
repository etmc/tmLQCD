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

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#ifdef MPI
# include <mpi.h>
#endif
#ifdef _USE_SHMEM
# include <mpp/shmem.h>
#endif
#include "global.h"
#include "read_input.h"
#include "mpi_init.h"

#ifdef MPI
/* Datatypes for the data exchange */
MPI_Datatype mpi_su3;
MPI_Datatype gauge_point;
MPI_Datatype gauge_time_slice_cont;
MPI_Datatype gauge_time_slice_split;
MPI_Datatype deri_point;
MPI_Datatype deri_time_slice_cont;
MPI_Datatype deri_time_slice_split;

MPI_Datatype field_point;
MPI_Datatype field_time_slice_cont;
MPI_Datatype lfield_time_slice_cont;
MPI_Datatype gauge_x_slice_cont;
MPI_Datatype gauge_x_subslice;
MPI_Datatype gauge_x_slice_gath;
MPI_Datatype field_x_slice_cont;
MPI_Datatype field_x_subslice;
MPI_Datatype field_x_slice_gath;
MPI_Datatype lfield_x_slice_cont;
MPI_Datatype lfield_x_subslice;
MPI_Datatype lfield_x_slice_gath;
MPI_Datatype deri_x_slice_cont;
MPI_Datatype deri_x_subslice;
MPI_Datatype deri_x_slice_gath;
MPI_Datatype gauge_xt_edge_cont;
MPI_Datatype gauge_xt_edge_gath;
MPI_Datatype deri_xt_edge_cont;

MPI_Datatype gauge_y_slice_gath;
MPI_Datatype gauge_y_slice_cont;
MPI_Datatype gauge_y_subslice;

MPI_Datatype field_y_slice_gath;
MPI_Datatype field_y_slice_cont;
MPI_Datatype field_y_subslice;
MPI_Datatype lfield_y_slice_gath;
MPI_Datatype lfield_y_slice_cont;
MPI_Datatype lfield_y_subslice;

MPI_Datatype field_z_slice_gath;
MPI_Datatype field_z_subslice;
MPI_Datatype field_z_slice_cont;
MPI_Datatype lfield_z_slice_gath;
MPI_Datatype lfield_z_slice_cont;
MPI_Datatype field_z_slice_half;

MPI_Datatype deri_y_slice_cont;
MPI_Datatype deri_y_subslice;
MPI_Datatype deri_y_slice_gath;

MPI_Datatype gauge_yx_edge_cont;
MPI_Datatype gauge_yx_edge_gath;
MPI_Datatype deri_yx_edge_cont;

MPI_Datatype gauge_ty_edge_cont;
MPI_Datatype gauge_ty_edge_gath;
MPI_Datatype deri_ty_edge_cont;

MPI_Datatype gauge_z_slice_gath;
MPI_Datatype gauge_z_slice_cont;
MPI_Datatype gauge_z_subslice;

MPI_Datatype deri_z_slice_cont;
MPI_Datatype deri_z_subslice;
MPI_Datatype deri_z_slice_gath;

MPI_Datatype gauge_zx_edge_cont;
MPI_Datatype gauge_zx_edge_gath;
MPI_Datatype deri_zx_edge_cont;

MPI_Datatype gauge_tz_edge_cont;
MPI_Datatype gauge_tz_edge_gath;
MPI_Datatype deri_tz_edge_cont;

MPI_Datatype gauge_zy_edge_cont;
MPI_Datatype gauge_zy_edge_gath;
MPI_Datatype deri_zy_edge_cont;

MPI_Datatype halffield_point;
MPI_Datatype halffield_time_slice_cont;

MPI_Datatype halffield_x_slice_cont;
MPI_Datatype halffield_x_subslice;
MPI_Datatype halffield_x_slice_gath;

MPI_Datatype halffield_y_slice_cont;
MPI_Datatype halffield_y_subslice;
MPI_Datatype halffield_y_slice_gath;

MPI_Datatype halffield_z_slice_cont;


#ifdef _USE_TSPLITPAR
MPI_Datatype field_xt_slice_int;
MPI_Datatype field_xt_slice_ext;
MPI_Datatype field_yt_slice_int;
MPI_Datatype field_yt_slice_ext;
# ifdef PARALLELXYZ
MPI_Datatype field_zt_slice_ext_L;
MPI_Datatype field_zt_slice_ext_S;
MPI_Datatype field_zt_slice_even_dn_et;
MPI_Datatype field_zt_slice_even_up_et;
MPI_Datatype field_zt_slice_odd_dn_et;
MPI_Datatype field_zt_slice_odd_up_et;
MPI_Datatype field_zt_slice_even_dn_ot;
MPI_Datatype field_zt_slice_even_up_ot;
MPI_Datatype field_zt_slice_odd_dn_ot;
MPI_Datatype field_zt_slice_odd_up_ot;
# endif
#endif
#ifdef WITHLAPH
MPI_Datatype su3vect_point;
MPI_Datatype jfield_x_slice_cont;
MPI_Datatype jfield_y_slice_cont;
MPI_Datatype jfield_z_slice_cont;
MPI_Datatype jfield_x_slice_gath;
MPI_Datatype jfield_y_slice_gath;
MPI_Datatype jfield_z_slice_gath;
MPI_Datatype jfield_y_subslice;
#endif

#if ( defined PARALLELXYZT || defined PARALLELXYZ )
MPI_Datatype field_z_slice_even_dn;
MPI_Datatype field_z_slice_even_up;
MPI_Datatype field_z_slice_odd_dn;
MPI_Datatype field_z_slice_odd_up;

# if (!defined _INDEX_INDEP_GEOM)
spinor * field_buffer_z ALIGN;
spinor * field_buffer_z2 ALIGN;
spinor * field_buffer_z3 ALIGN;
spinor * field_buffer_z4 ALIGN;
halfspinor * halffield_buffer_z ALIGN;
halfspinor * halffield_buffer_z2 ALIGN;
# endif
#endif

MPI_Op mpi_reduce_su3_ray;

void reduce_su3_ray(
    void *u_i        /* in */,
    void *u_io       /* in/out */,
    int *len         /* in */,
    MPI_Datatype *dt /* in */) {

   int n;
   su3 *u, *v, tmp;
   u = (su3 *)u_i;
   v = (su3 *)u_io;

   if(*dt != mpi_su3) {
      fprintf(stderr, "\nInvalid datatype for reduce_su3_ray(); abort.\n");
      MPI_Abort(MPI_COMM_WORLD, 1);
   }
   for(n=0; n<*len; n++) {
      _su3_times_su3(tmp,*(u+n),*(v+n))
      _su3_assign(*(v+n),tmp)
   }
}

#endif


void tmlqcd_mpi_init(int argc,char *argv[]) {
  int i;
#ifdef MPI
  int periods[] = {1,1,1,1};
  int dims[] = {0,0,0,0};
  int ndims = 0;
  int nalldims = 4;
  int reorder = 1, namelen;
  char processor_name[MPI_MAX_PROCESSOR_NAME];
#endif
  g_proc_coords[0] = 0;
  g_proc_coords[1] = 0;
  g_proc_coords[2] = 0;
  g_proc_coords[3] = 0;
  for(i = 0; i < 8; i++) {
    g_nb_list[i] = 0;
  }


#ifdef MPI
#  ifdef _USE_SHMEM
  /* we need that the PE number in MPI_COMM_WORL  */
  /* exactly correspond to the one in g_cart_grid */
  reorder = 0;
#  endif

#    ifndef FIXEDVOLUME
  N_PROC_T=0; /* the other N_PROC_? are read from input, if not constraint below */
              /* N_PROC_T will be set by MPI_Dims_create, if not constraint below */
#    endif

#  if defined PARALLELT
  ndims = 1;
#    ifndef FIXEDVOLUME
  N_PROC_X = 1;
  N_PROC_Y = 1;
  N_PROC_Z = 1;
#    endif
#  endif
#  if defined PARALLELX
  ndims = 1;
#    ifndef FIXEDVOLUME
  N_PROC_T = 1;
  N_PROC_Y = 1;
  N_PROC_Z = 1;
#    endif
#  endif
#  if defined PARALLELXT
  ndims = 2;
#    ifndef FIXEDVOLUME
  N_PROC_Y = 1;
  N_PROC_Z = 1;
#    endif
#  endif
#  if defined PARALLELXY
  ndims = 2;
#    ifndef FIXEDVOLUME
  N_PROC_T = 1;
  N_PROC_Z = 1;
#    endif
#  endif
#  if defined PARALLELXYT
  ndims = 3;
#    ifndef FIXEDVOLUME
  N_PROC_Z = 1;
#    endif
#  endif
#  if defined PARALLELXYZ
  ndims = 3;
#    ifndef FIXEDVOLUME
  N_PROC_T = 1;
#    endif
#  endif
#  if defined PARALLELXYZT
  ndims = 4;
#  endif
  dims[0] = N_PROC_T;
  dims[1] = N_PROC_X;
  dims[2] = N_PROC_Y;
  dims[3] = N_PROC_Z;


  MPI_Comm_size(MPI_COMM_WORLD, &g_nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &g_proc_id);
  MPI_Get_processor_name(processor_name, &namelen);
  MPI_Dims_create(g_nproc, nalldims, dims);
  if(g_proc_id == 0){
    printf("# Creating the following cartesian grid for a %d dimensional parallelisation:\n# %d x %d x %d x %d\n"
            , ndims, dims[0], dims[1], dims[2], dims[3]);
  }

  g_nproc_t = dims[0];
  g_nproc_x = dims[1];
  g_nproc_y = dims[2];
  g_nproc_z = dims[3];

  if( (g_nproc_t < 1 || g_nproc_x < 1 || g_nproc_y < 1 || g_nproc_z < 1) ||
      (LX%g_nproc_x != 0 || LY%g_nproc_y != 0 || LZ%g_nproc_z != 0 || T_global%g_nproc_t != 0) ) {
    if(g_proc_id == 0) {
      fprintf(stderr, "The lattice cannot be properly mapped on the processor grid\n");
      fprintf(stderr, "Please check your number of processors and the Nr?Procs input variables\n");
      fprintf(stderr, "Aborting...!\n");
    }
    MPI_Abort(MPI_COMM_WORLD, 1);
    MPI_Finalize();
    exit(-1);
  }

#  ifndef FIXEDVOLUME
  N_PROC_T = g_nproc_t;
  N_PROC_X = g_nproc_x;
  N_PROC_Y = g_nproc_y;
  N_PROC_Z = g_nproc_z;
  T = T_global/g_nproc_t;
  LX = LX/g_nproc_x;
  LY = LY/g_nproc_y;
  LZ = LZ/g_nproc_z;
  VOLUME = (T*LX*LY*LZ);
  SPACEVOLUME = VOLUME/T;
#    ifdef _USE_TSPLITPAR
  TEOSLICE = (LX*LY*LZ)/2;
#    endif
#    ifdef PARALLELT  
  RAND = (2*LX*LY*LZ);
  EDGES = 0;
#    elif defined PARALLELX  
  RAND = (2*T*LY*LZ);
  EDGES = 0;
#    elif defined PARALLELXT
  RAND = 2*LZ*(LY*LX + T*LY);
  EDGES = 4*LZ*LY;
#    elif defined PARALLELXY
  RAND = 2*LZ*T*(LX + LY);
  EDGES = 4*LZ*T;
#    elif defined PARALLELXYT
  RAND = 2*LZ*(LY*LX + T*LY + T*LX);
  EDGES = 4*LZ*(LY + T + LX);
#    elif defined PARALLELXYZ
  RAND = 2*T*(LY*LZ + LX*LZ + LX*LY);
  EDGES = 4*T*(LX + LY + LZ);
#    elif defined PARALLELXYZT
  RAND = 2*LZ*LY*LX + 2*LZ*T*LY + 2*LZ*T*LX + 2*T*LX*LY;
  EDGES = 4*LZ*LY + 4*LZ*T + 4*LZ*LX + 4*LY*T + 4*LY*LX + 4*T*LX;
#    else /* ifdef PARALLELT */
  RAND = 0;
  EDGES = 0;
#    endif /* ifdef PARALLELT */
  /* Note that VOLUMEPLUSRAND is not always equal to VOLUME+RAND */
  /* VOLUMEPLUSRAND rather includes the edges */
  VOLUMEPLUSRAND = VOLUME + RAND + EDGES;
  SPACERAND=RAND/T;
#  endif /* ifndef FIXEDVOLUME */
  g_dbw2rand = (RAND + 2*EDGES);

#  if (!defined _INDEX_INDEP_GEOM)
#   if ( defined PARALLELXYZT ||  defined PARALLELXYZ )
  field_buffer_z = (spinor*)malloc(T*LX*LY/2*sizeof(spinor));
  field_buffer_z2 = (spinor*)malloc(T*LX*LY/2*sizeof(spinor));
#    ifdef _NON_BLOCKING
  field_buffer_z3 = (spinor*)malloc(T*LX*LY/2*sizeof(spinor));
  field_buffer_z4 = (spinor*)malloc(T*LX*LY/2*sizeof(spinor));
#    endif
  halffield_buffer_z = (halfspinor*)malloc(T*LX*LY/2*sizeof(halfspinor));
  halffield_buffer_z2 = (halfspinor*)malloc(T*LX*LY/2*sizeof(halfspinor));
#  endif
# endif

  MPI_Cart_create(MPI_COMM_WORLD, nalldims, dims, periods, reorder, &g_cart_grid);
  MPI_Comm_rank(g_cart_grid, &g_cart_id);
  MPI_Cart_coords(g_cart_grid, g_cart_id, nalldims, g_proc_coords);
  if (g_debug_level > 1) {
    fprintf(stdout,"# Process %d of %d on %s: cart_id %d, coordinates (%d %d %d %d)\n",
            g_proc_id, g_nproc, processor_name, g_cart_id, 
            g_proc_coords[0], g_proc_coords[1], g_proc_coords[2], g_proc_coords[3]);
    fflush(stdout);
  }
  if(g_stdio_proc == -1){
    g_stdio_proc = g_proc_id;
  }
  for(i = 0; i < 8; i++) {
    g_nb_list[i] = g_cart_id;
  }
#  if (defined PARALLELT || defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT)
  MPI_Cart_shift(g_cart_grid, 0, 1, &g_nb_t_dn, &g_nb_t_up);
  g_nb_list[0] = g_nb_t_up;  
  g_nb_list[1] = g_nb_t_dn;
#  endif
#  if (defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT || defined PARALLELX || defined PARALLELXY || defined PARALLELXYZ )
  MPI_Cart_shift(g_cart_grid, 1, 1, &g_nb_x_dn, &g_nb_x_up);
  g_nb_list[2] = g_nb_x_up;  
  g_nb_list[3] = g_nb_x_dn;
#  endif
#  if (defined PARALLELXYT  || defined PARALLELXYZT || defined PARALLELXY || defined PARALLELXYZ )
  MPI_Cart_shift(g_cart_grid, 2, 1, &g_nb_y_dn, &g_nb_y_up);
  g_nb_list[4] = g_nb_y_up;  
  g_nb_list[5] = g_nb_y_dn;
#  endif
#  if (defined PARALLELXYZT || defined PARALLELXYZ ) 
  MPI_Cart_shift(g_cart_grid, 3, 1, &g_nb_z_dn, &g_nb_z_up);
  g_nb_list[6] = g_nb_z_up;  
  g_nb_list[7] = g_nb_z_dn;
#  endif


#  if ((defined _INDEX_INDEP_GEOM) && (defined _USE_HALFSPINOR))
#   if (defined PARALLELT || defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT)
  g_HS_shift_t = 0;
  g_HS_shift_x = LX*LY*LZ;
  g_HS_shift_y = LX*LY*LZ + T*LY*LZ;
  g_HS_shift_z = LX*LY*LZ + T*LY*LZ + T*LX*LZ;
#   endif
#   if (defined PARALLELX || defined PARALLELXY || defined PARALLELXYZ )
  g_HS_shift_t = 0;
  g_HS_shift_x = 0;
  g_HS_shift_y = T*LY*LZ;
  g_HS_shift_z = T*LY*LZ + T*LX*LZ;
#   endif
#  endif

  /* With internal boundary we mean the fields that are send */
  /* to another processor. It is located wihtin the local    */
  /* volume, whereas the external boundary is the boundary   */
  /* received from another processor lying on the RAND.      */
  /* In general the external bondaries are continuous in     */
  /* memory, while this is not always true for the internal  */
  /* one. */

  /* first the gauge fields */
  MPI_Type_contiguous(18, MPI_DOUBLE, &mpi_su3);
  MPI_Type_commit(&mpi_su3);
  /* This is a gauge field on one space-time point */
  MPI_Type_contiguous(4, mpi_su3, &gauge_point);
  /* This is a type for one gauge time slice continuous */ 
  MPI_Type_contiguous(LX*LY*LZ, gauge_point, &gauge_time_slice_cont);
  /* This is a type for one gauge time slice dis-continuous -> NEW_GEOMETRY */
  /* This are 2 continuous ensembles of gauge_points of length LX*LY*LZ/2   */
  /* separated in memory by (VOLUME)/2 gauge_points                         */
  MPI_Type_vector(2, LX*LY*LZ/2, (VOLUME)/2, gauge_point, &gauge_time_slice_split);
  /* Commit the new types */
  MPI_Type_commit(&gauge_time_slice_split);
  MPI_Type_commit(&gauge_time_slice_cont);

  /* Continuous x-slice as it is found in the external memory.*/
  MPI_Type_contiguous(T*LY*LZ, gauge_point, &gauge_x_slice_cont);
  /* this is a continuous gauge xt-slice */
  MPI_Type_contiguous(LY*LZ, gauge_point, &gauge_x_subslice);
  /* Put T of the latter together, each of which has length 1 (in units */
  /* of gauge_yy_subslice). They are separated by LX of those.          */
  /* This is as the gauge fields are located in the internal memory     */
  MPI_Type_vector(T, 1, LX, gauge_x_subslice, &gauge_x_slice_gath);
  MPI_Type_commit(&gauge_x_slice_gath);
  MPI_Type_commit(&gauge_x_slice_cont);

  /* Continuous y-slice as it is found in the external memory.*/
  MPI_Type_contiguous(T*LX*LZ, gauge_point, &gauge_y_slice_cont);
  /* this is a continuous gauge xyt-slice */
  MPI_Type_contiguous(LZ, gauge_point, &gauge_y_subslice);
  /* Put T*LX together, separated by LY of those */
  MPI_Type_vector(T*LX, 1, LY, gauge_y_subslice, &gauge_y_slice_gath);
  MPI_Type_commit(&gauge_y_slice_cont);
  MPI_Type_commit(&gauge_y_slice_gath);

  /* Continuous z-slice as it is found in the external memory.*/
  MPI_Type_contiguous(T*LX*LY, gauge_point, &gauge_z_slice_cont);
  /* Put T*LX*LY gauge-points together, separated by LZ of those */
  MPI_Type_vector(T*LX*LY, 1, LZ, gauge_point, &gauge_z_slice_gath);
  MPI_Type_commit(&gauge_z_slice_cont);
  MPI_Type_commit(&gauge_z_slice_gath);


  /* external edges: on x-Rand send in t-direction*/
  MPI_Type_contiguous(2*LY*LZ ,gauge_point, &gauge_xt_edge_cont);
  MPI_Type_commit(&gauge_xt_edge_cont);
  /* internal edges, lying in memory nevertheless in the boundary */
  MPI_Type_vector(2, 1, T, gauge_x_subslice, &gauge_xt_edge_gath);
  MPI_Type_commit(&gauge_xt_edge_gath);

  /* external edges: y-Rand send in x-direction */
  MPI_Type_contiguous(2*T*LZ ,gauge_point, &gauge_yx_edge_cont);
  MPI_Type_commit(&gauge_yx_edge_cont);
  /* internal edges */
  MPI_Type_vector(2*T, LZ, LX*LZ, gauge_point, &gauge_yx_edge_gath);
  MPI_Type_commit(&gauge_yx_edge_gath);

  /* external edges: t-Rand send in y-direction */
  MPI_Type_contiguous(2*LX*LZ ,gauge_point, &gauge_ty_edge_cont);
  MPI_Type_commit(&gauge_ty_edge_cont);
  /* internal edges */
  MPI_Type_vector(2*LX, LZ, LY*LZ, gauge_point, &gauge_ty_edge_gath);
  MPI_Type_commit(&gauge_ty_edge_gath);

  /* external edges: z-Rand send in x-direction */
  /* zx-edge */
  MPI_Type_contiguous(2*T*LY ,gauge_point, &gauge_zx_edge_cont);
  MPI_Type_commit(&gauge_zx_edge_cont);
  /* internal edges */
  MPI_Type_vector(2*T, LY, LY*LX, gauge_point, &gauge_zx_edge_gath);
  MPI_Type_commit(&gauge_zx_edge_gath);

  /* external edges: t-Rand send in z-direction */
  /* tz-edge */
  MPI_Type_contiguous(2*LX*LY ,gauge_point, &gauge_tz_edge_cont);
  MPI_Type_commit(&gauge_tz_edge_cont);
  /* internal edges */
  MPI_Type_vector(2*LX*LY, 1, LZ, gauge_point, &gauge_tz_edge_gath);
  MPI_Type_commit(&gauge_tz_edge_gath);

  /* external edges: z-Rand send in y-direction */
  /* zy-edge */
  MPI_Type_contiguous(2*T*LX ,gauge_point, &gauge_zy_edge_cont);
  MPI_Type_commit(&gauge_zy_edge_cont);
  /* internal edges */
  MPI_Type_vector(2*T*LX, 1, LY, gauge_point, &gauge_zy_edge_gath);
  MPI_Type_commit(&gauge_zy_edge_gath);

  /* The spinor fields */
  /* this is a single spinor field on one space-time point */
  MPI_Type_contiguous(24, MPI_DOUBLE, &field_point);
  /* Tis is an even or odd spinor field time slice, continuous */
/*   MPI_Type_contiguous(LX*LY*LZ/2, field_point, &field_time_slice_cont);  */
  MPI_Type_contiguous(LX*LY*LZ*12, MPI_DOUBLE, &field_time_slice_cont); 
  /* Commit the new types */
  MPI_Type_commit(&field_time_slice_cont);

  /* this is the not even/odd field */
  MPI_Type_contiguous(LX*LY*LZ, field_point, &lfield_time_slice_cont);
  MPI_Type_commit(&lfield_time_slice_cont);


  /* This is an even or odd continuous spinor field x-slice */
  MPI_Type_contiguous(T*LY*LZ/2, field_point, &field_x_slice_cont); 
/*   MPI_Type_contiguous(12*T*LY*LZ, MPI_DOUBLE, &field_x_slice_cont); */
  /* this is an even or odd continuous spinor field xt-slice */
  MPI_Type_contiguous(LY*LZ/2, field_point, &field_x_subslice);
  /* this type puts T xt-slices together being the internal x-boundary in */
  /* even/odd ordered spinor fields */
  MPI_Type_vector(T, 1, LX, field_x_subslice, &field_x_slice_gath); 
/*   MPI_Type_vector(T, 12*LY*LZ, 12*LX*LY*LZ, MPI_DOUBLE, &field_x_slice_gath); */
  MPI_Type_commit(&field_x_slice_gath);
  MPI_Type_commit(&field_x_slice_cont);

  /* this is the not even/odd field */
  MPI_Type_contiguous(T*LY*LZ, field_point, &lfield_x_slice_cont);
  MPI_Type_contiguous(LY*LZ, field_point, &lfield_x_subslice);
  MPI_Type_vector(T, 1, LX, lfield_x_subslice, &lfield_x_slice_gath);
  MPI_Type_commit(&lfield_x_slice_gath);
  MPI_Type_commit(&lfield_x_slice_cont);

  /* This is an even or odd continuous spinor field y-slice */
  MPI_Type_contiguous(T*LX*LZ/2, field_point, &field_y_slice_cont); 
/*   MPI_Type_contiguous(12*T*LX*LZ, MPI_DOUBLE, &field_y_slice_cont); */
  /* this is an even or odd continuous spinor field txy-slice */
  MPI_Type_contiguous(LZ/2, field_point, &field_y_subslice);
  /* this type puts T*LX xt-slices together being the internal y-boundary in */
  /* even/odd ordered spinor fields */
  MPI_Type_vector(T*LX, 1, LY, field_y_subslice, &field_y_slice_gath); 
/*   MPI_Type_vector(T*LX, 12*LZ, 12*LY*LZ, MPI_DOUBLE, &field_y_slice_gath); */
  MPI_Type_commit(&field_y_slice_gath);
  MPI_Type_commit(&field_y_slice_cont);

  /* this is the not even/odd field */
  MPI_Type_contiguous(T*LX*LZ, field_point, &lfield_y_slice_cont);
  MPI_Type_contiguous(LZ, field_point, &lfield_y_subslice);
  MPI_Type_vector(T*LX, 1, LY, lfield_y_subslice, &lfield_y_slice_gath);
  MPI_Type_commit(&lfield_y_slice_cont);
  MPI_Type_commit(&lfield_y_slice_gath);

  /* If z-dir is parallelized, I have assumed that both LZ and T*LX*LY are even */
  /* This is an even or odd continuous spinor field z-slice */
  MPI_Type_contiguous(T*LX*LY/2, field_point, &field_z_slice_cont);

  /* this type puts T*LX*LY field_point together being the internal z-boundary in */
  /* even/odd ordered spinor fields */
  MPI_Type_vector(T*LX*LY/2, 12, 24, MPI_DOUBLE, &field_z_slice_half); /* this is ?!? (Not used) */
  MPI_Type_commit(&field_z_slice_half);
  MPI_Type_commit(&field_z_slice_cont);

  /* this is the not even/odd field */
  MPI_Type_contiguous(T*LX*LY, field_point, &lfield_z_slice_cont);
  MPI_Type_vector(T*LX*LY, 1, LZ, field_point, &lfield_z_slice_gath);
  MPI_Type_commit(&lfield_z_slice_cont);
  MPI_Type_commit(&lfield_z_slice_gath);

#ifdef _USE_TSPLITPAR
  /* here I construct the xt yt zt edges for use in _USE_TSPLITPAR  */
  MPI_Type_contiguous(LY*LZ/2, field_point, &field_xt_slice_int); /* OK */
  MPI_Type_vector(LX, LZ/2, LY*LZ/2, field_point, &field_yt_slice_int);  /* OK */
  MPI_Type_contiguous(LY*LZ/2, field_point, &field_xt_slice_ext);  /* OK */
  MPI_Type_contiguous(LX*LZ/2, field_point, &field_yt_slice_ext);  /* OK */
  MPI_Type_commit(&field_xt_slice_int);
  MPI_Type_commit(&field_xt_slice_ext);
  MPI_Type_commit(&field_yt_slice_int);
  MPI_Type_commit(&field_yt_slice_ext);
# ifdef PARALLELXYZ
  MPI_Type_contiguous((LX*LY+1)/2, field_point, &field_zt_slice_ext_L);  /* OK */
  MPI_Type_contiguous(LX*LY/2, field_point, &field_zt_slice_ext_S);  /* OK */
  MPI_Type_commit(&field_zt_slice_ext_L);
  MPI_Type_commit(&field_zt_slice_ext_S);
# endif
#endif

#ifdef WITHLAPH
  MPI_Type_contiguous(6, MPI_DOUBLE, &su3vect_point); 

  MPI_Type_contiguous(LY*LZ, su3vect_point, &jfield_x_slice_cont);
  MPI_Type_contiguous(LX*LZ, su3vect_point, &jfield_y_slice_cont);
  MPI_Type_contiguous(LX*LY, su3vect_point, &jfield_z_slice_cont);
  MPI_Type_contiguous(LY*LZ, su3vect_point, &jfield_x_slice_gath);
  MPI_Type_contiguous(LZ, su3vect_point, &jfield_y_subslice);
  MPI_Type_vector(LX, 1, LY, jfield_y_subslice, &jfield_y_slice_gath);
  MPI_Type_vector(LX*LY, 1, LZ, su3vect_point, &jfield_z_slice_gath);
  MPI_Type_commit(&jfield_x_slice_gath);
  MPI_Type_commit(&jfield_x_slice_cont);
  MPI_Type_commit(&jfield_y_slice_cont);
  MPI_Type_commit(&jfield_y_slice_gath);
  MPI_Type_commit(&jfield_z_slice_cont);
  MPI_Type_commit(&jfield_z_slice_gath);
#endif

  /* The internal z_ and zt_ slices are constructed in geometry() with MPI_Type_indexed() */

  /* Now the derivative fields */
  /* this is a derivative field on one space-time point */
  MPI_Type_contiguous(32, MPI_DOUBLE, &deri_point);
  /* This is a type for one derivative time slice continuous */
  MPI_Type_contiguous(LX*LY*LZ, deri_point, &deri_time_slice_cont);
  /* This is a type for one derivative time slice dis-continuous -> NEW_GEOMETRY */
  MPI_Type_vector(2, LX*LY*LZ/2, VOLUME/2, deri_point, &deri_time_slice_split);
  /* Commit the new types */
  MPI_Type_commit(&deri_time_slice_split);
  MPI_Type_commit(&deri_time_slice_cont);

  MPI_Type_contiguous(T*LY*LZ, deri_point, &deri_x_slice_cont);
  MPI_Type_contiguous(LY*LZ, deri_point, &deri_x_subslice);
  MPI_Type_vector(T, 1, LX, deri_x_subslice, &deri_x_slice_gath);
  MPI_Type_commit(&deri_x_slice_gath);
  MPI_Type_commit(&deri_x_slice_cont);

  /* external edges: on x-boundary send in t-direction first */
  MPI_Type_contiguous(2*LY*LZ ,deri_point, &deri_xt_edge_cont);
  MPI_Type_commit(&deri_xt_edge_cont);
  /* external edges: y-boundary send in x-direction */
  MPI_Type_contiguous(2*T*LZ ,deri_point, &deri_yx_edge_cont);
  MPI_Type_commit(&deri_yx_edge_cont);
  /* external edges: t-boundary send in y-direction */
  MPI_Type_contiguous(2*LX*LZ ,deri_point, &deri_ty_edge_cont);
  MPI_Type_commit(&deri_ty_edge_cont);
  /* external edges: z-boundary send in x-direction */
  MPI_Type_contiguous(2*T*LY ,deri_point, &deri_zx_edge_cont);
  MPI_Type_commit(&deri_zx_edge_cont);
  /* external edges: t-boundary send in z-direction */
  MPI_Type_contiguous(2*LX*LY ,deri_point, &deri_tz_edge_cont);
  MPI_Type_commit(&deri_tz_edge_cont);
  /* external edges: z-boundary send in y-direction */
  MPI_Type_contiguous(2*T*LX ,deri_point, &deri_zy_edge_cont);
  MPI_Type_commit(&deri_zy_edge_cont);

  MPI_Type_contiguous(T*LX*LZ, deri_point, &deri_y_slice_cont);
  MPI_Type_contiguous(LZ, deri_point, &deri_y_subslice);
  MPI_Type_vector(T*LX, 1, LY, deri_y_subslice, &deri_y_slice_gath);
  MPI_Type_commit(&deri_y_slice_gath);
  MPI_Type_commit(&deri_y_slice_cont);

  MPI_Type_contiguous(T*LX*LY, deri_point, &deri_z_slice_cont);
  MPI_Type_vector(T*LX*LY, 1, LZ, deri_point, &deri_z_slice_gath);
  MPI_Type_commit(&deri_z_slice_gath);
  MPI_Type_commit(&deri_z_slice_cont);

  /* this is a single halfspinor field on one space-time point */
  MPI_Type_contiguous(12, MPI_DOUBLE, &halffield_point);
  MPI_Type_vector(LX*LY*LZ/2, 1, 8, halffield_point, &halffield_time_slice_cont); 

  /* Commit the new types */
  MPI_Type_commit(&halffield_time_slice_cont);

  MPI_Type_vector(LY*LZ/2, 1, 8, halffield_point, &halffield_x_subslice);
  MPI_Type_vector(T, 1, LX, halffield_x_subslice, &halffield_x_slice_gath); 
  MPI_Type_commit(&halffield_x_slice_gath);

  MPI_Type_vector(LZ/2, 1, 8, halffield_point, &halffield_y_subslice);
  MPI_Type_vector(T*LX, 1, LY, halffield_y_subslice, &halffield_y_slice_gath); 
  MPI_Type_commit(&halffield_y_slice_gath);

  /* For observables we need communicators for Cartesian time slices */
  MPI_Comm_split(g_cart_grid, g_proc_coords[0], g_cart_id, &g_mpi_time_slices);
  MPI_Comm_rank(g_mpi_time_slices, &g_mpi_time_rank);
  if(g_debug_level > 4) {
    fprintf(stdout, "# My mpi_time_rank = %d, g_proc_coords = (%d,%d,%d,%d), g_cart_id = %d\n", 
	    g_mpi_time_rank, g_proc_coords[0], g_proc_coords[1], g_proc_coords[2], g_proc_coords[3],
	    g_cart_id);
  }

  /* and communicators for Cartesian z-slices */
  MPI_Comm_split(g_cart_grid, g_proc_coords[3], g_cart_id, &g_mpi_z_slices);
  MPI_Comm_rank(g_mpi_z_slices, &g_mpi_z_rank);
  if(g_debug_level > 4) {
    fprintf(stdout, "# My mpi_z_rank = %d, g_proc_coords = (%d,%d,%d,%d), g_cart_id = %d\n", 
	    g_mpi_z_rank, g_proc_coords[0], g_proc_coords[1], g_proc_coords[2], g_proc_coords[3],
	    g_cart_id);
  }

  /* and spatial volume slices */
  MPI_Comm_split(g_cart_grid, g_mpi_time_rank, g_proc_coords[0], &g_mpi_SV_slices);
  MPI_Comm_rank(g_mpi_SV_slices, &g_mpi_SV_rank);
  if(g_debug_level > 4) {
    fprintf(stdout, "# My mpi_SV_rank = %d, g_proc_coords = (%d,%d,%d,%d), g_cart_id = %d\n", 
	    g_mpi_SV_rank, g_proc_coords[0], g_proc_coords[1], g_proc_coords[2], g_proc_coords[3],
	    g_cart_id);
  }

  /* and tim-volume slices orthogonal to the z-direction */
  MPI_Comm_split(g_cart_grid, g_mpi_z_rank, g_proc_coords[3], &g_mpi_ST_slices);
  MPI_Comm_rank(g_mpi_ST_slices, &g_mpi_ST_rank);
  if(g_debug_level > 4) {
    fprintf(stdout, "# My mpi_ST_rank = %d, g_proc_coords = (%d,%d,%d,%d), g_cart_id = %d\n", 
	    g_mpi_ST_rank, g_proc_coords[0], g_proc_coords[1], g_proc_coords[2], g_proc_coords[3],
	    g_cart_id);
  }

  MPI_Op_create(reduce_su3_ray, 0, &mpi_reduce_su3_ray);

#else /*ifdef MPI */
  g_nproc = 1;
  g_proc_id = 0;
  g_nproc_x = 1;
  g_nproc_y = 1;
  g_nproc_z = 1;
  g_nproc_t = 1;
  g_cart_id = 0;
  g_mpi_time_rank = 0;
  g_mpi_z_rank = 0;
  g_mpi_SV_rank = 0;
  g_mpi_ST_rank = 0;
  g_stdio_proc = 0;

#  ifndef FIXEDVOLUME
  T = T_global;
  VOLUME = (T*LX*LY*LZ);
  SPACEVOLUME = VOLUME/T;
#    ifdef _USE_TSPLITPAR
  TEOSLICE = (LX*LY*LZ)/2;
#    endif
  RAND = 0;
  EDGES = 0;
  VOLUMEPLUSRAND = VOLUME;
  SPACERAND=0;
  N_PROC_T = 1;
  N_PROC_X = 1;
  N_PROC_Y = 1;
  N_PROC_Z = 1;
#  endif
  g_dbw2rand = 0;
#endif   /*ifdef MPI */

  /* Here we perform some checks in order not to */
  /* run into trouble later                      */
#if (defined PARALLELXYZT || defined PARALLELXYZ )
  if((T*LX*LY)%2 != 0 && even_odd_flag == 1) {
    fprintf(stderr, "T*LX*LY must be even!\nAborting prgram...\n");
#  ifdef MPI 
    MPI_Finalize();
#  endif
    exit(-1);
  }
#endif

  if(LZ%2 != 0 && even_odd_flag == 1) {
    fprintf(stderr, "LZ must be even!\nAborting prgram...\n");
#ifdef MPI
    MPI_Finalize();
#endif
    exit(-1);
  }
}


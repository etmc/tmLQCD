/* $Id$ */
#include <stdlib.h>
#include <stdio.h>
#ifdef MPI
#include <mpi.h>
#endif
#include "global.h"
#include "mpi_init.h"


void mpi_init(int argc,char *argv[]) {
#ifdef MPI
  int periods[] = {1,1};
  int dims[] = {0,0};
  int ndims = 0;
  int reorder = 1;
  int  namelen;
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  
  g_proc_coords[0] = 0;
  g_proc_coords[1] = 0;
  dims[1] = N_PROC_X;

#ifdef PARALLELT
  ndims = 1;
#endif
#if defined PARALLELXT
  ndims = 2;
#endif

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &g_nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &g_proc_id);
  MPI_Get_processor_name(processor_name, &namelen);
  MPI_Dims_create(g_nproc, ndims, dims);
  if(g_proc_id == 0){
    printf("Creating the following cartesian grid for a %d dimensional parallelisation:\n%d x %d \n"
	   , ndims, dims[0], dims[1]);
  }
  fprintf(stdout,"Process %d of %d on %s\n",
	  g_proc_id, g_nproc, processor_name);
  fflush(stdout);
  g_nproc_t = dims[0];
  MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorder, &g_cart_grid);
  MPI_Comm_rank(g_cart_grid, &g_proc_id);
  MPI_Cart_coords(g_cart_grid, g_proc_id, ndims, g_proc_coords);

  MPI_Comm_rank(MPI_COMM_WORLD, &g_proc_id);
  if(g_stdio_proc == -1){
    g_stdio_proc = g_proc_id;
  }

  MPI_Cart_shift(g_cart_grid, 0, 1, &g_nb_t_dn, &g_nb_t_up);
#ifdef PARALLELXT
  MPI_Cart_shift(g_cart_grid, 1, 1, &g_nb_x_dn, &g_nb_x_up);
#endif

  MPI_Type_contiguous(72, MPI_DOUBLE, &g_gauge_point);
  MPI_Type_commit(&g_gauge_point);
  MPI_Type_contiguous(LX*L*L, g_gauge_point, &g_gauge_time_slice_cont);
  MPI_Type_vector(2, LX*L*L/2, (VOLUME)/2, g_gauge_point, &g_gauge_time_slice_split);
  MPI_Type_commit(&g_gauge_time_slice_split);
  MPI_Type_commit(&g_gauge_time_slice_cont);

  MPI_Type_contiguous(32, MPI_DOUBLE, &g_deri_point);
  MPI_Type_commit(&g_deri_point);
  MPI_Type_vector(2, LX*L*L/2, VOLUME/2, g_deri_point, &g_deri_time_slice_split);
  MPI_Type_contiguous(LX*L*L, g_deri_point, &g_deri_time_slice_cont);
  MPI_Type_commit(&g_deri_time_slice_split);
  MPI_Type_commit(&g_deri_time_slice_cont);

  MPI_Type_contiguous(24, MPI_DOUBLE, &g_field_point);
  MPI_Type_commit(&g_field_point);
  MPI_Type_contiguous(LX*L*L/2, g_field_point, &g_field_time_slice_cont);
  MPI_Type_commit(&g_field_time_slice_cont);

  MPI_Type_contiguous(T*LY*LZ, g_gauge_point, &g_gauge_x_slice_cont);
  MPI_Type_contiguous(LY*LZ, g_gauge_point, &g_gauge_yz_subslice);
  MPI_Type_commit(&g_gauge_yz_subslice);
  MPI_Type_vector(T, 1, LX, g_gauge_yz_subslice, &g_gauge_x_slice_gath);
  MPI_Type_commit(&g_gauge_x_slice_gath);
  MPI_Type_commit(&g_gauge_x_slice_cont);

  MPI_Type_contiguous(T*LY*LZ/2, g_field_point, &g_field_x_slice_cont);
  MPI_Type_commit(&g_field_x_slice_cont);
  MPI_Type_contiguous(LY*LZ/2, g_field_point, &g_field_yz_subslice);
  MPI_Type_commit(&g_field_yz_subslice); 
  MPI_Type_vector(T, 1, LX, g_field_yz_subslice, &g_field_x_slice_gath);
  MPI_Type_commit(&g_field_x_slice_gath);

  MPI_Type_contiguous(T*LY*LZ, g_deri_point, &g_deri_x_slice_cont);
  MPI_Type_contiguous(LY*LZ, g_deri_point, &g_deri_yz_subslice);
  MPI_Type_commit(&g_deri_yz_subslice);
  MPI_Type_vector(T, 1, LX, g_deri_yz_subslice, &g_deri_x_slice_gath);
  MPI_Type_commit(&g_deri_x_slice_gath);
  MPI_Type_commit(&g_deri_x_slice_cont);

#else
  g_nproc = 1;
  g_proc_id = 0;

  g_proc_coords[0] = 0;
  g_proc_coords[1] = 0;
#endif
}

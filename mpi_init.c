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
  int  namelen;
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &g_nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &g_proc_id);
  MPI_Get_processor_name(processor_name, &namelen);
  g_cart_grid = MPI_COMM_WORLD;

  fprintf(stdout,"Process %d of %d on %s\n",
	  g_proc_id, g_nproc, processor_name);
  fflush(stdout);
#else
  g_nproc = 1;
  g_proc_id = 0;
#endif
}

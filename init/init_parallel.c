/***********************************************************************
 *
 * Copyright (C) 2017  Bartosz Kostrzewa, Carsten Urbach
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
 *
 *******************************************************************************/

#ifdef HAVE_CONFIG_H
#include "tmlqcd_config.h"
#endif
#ifdef TM_USE_MPI
#include <mpi.h>
#endif
#ifdef TM_USE_OMP
#include <omp.h>
#include "init/init_openmp.h"
#endif
#ifdef TM_USE_QPHIX
#include "qphix/qphix_config.h"
#endif
#ifdef QPHIX_QMP_COMMS
#include <qmp.h>
#endif

#include "init_parallel.h"
#include "global.h"
#include "read_input.h"

void init_parallel_and_read_input(int argc, char *argv[], char input_filename[]) {
#ifdef QPHIX_QMP_COMMS
  // Initialize QMP
  QMP_thread_level_t prv;
  if (QMP_init_msg_passing(&argc, &argv, QMP_THREAD_SINGLE, &prv) != QMP_SUCCESS) {
    QMP_error("Failed to initialize QMP\n");
    abort();
  }
  if (QMP_is_primary_node()) {
    printf("QMP IS INITIALIZED\n");
  }
#elif defined(TM_USE_MPI) && !defined(QPHIX_QMP_COMMS)
#ifdef TM_USE_OMP
  int mpi_thread_provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_SERIALIZED, &mpi_thread_provided);
#else
  MPI_Init(&argc, &argv);
#endif
#endif  // QPHIX_QMP_COMMS

#if defined(TM_USE_MPI) || defined(QPHIX_QMP_COMMS)
  MPI_Comm_rank(MPI_COMM_WORLD, &g_proc_id);
#else
  g_proc_id = 0;
#endif

// Read the input file
int status = read_input(input_filename);
if (status != 0) {
  fprintf(stderr, "Could not find input file: %s\nAborting...\n", input_filename);
  exit(-1);
}

#ifdef TM_USE_OMP
  init_openmp();
#endif
}

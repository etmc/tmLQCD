/***********************************************************************
 *
 * Copyright (C) 2026 JingJing Li
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
 *******************************************************************************/

#include <stdbool.h>
#include <sys/types.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include "ptbc.h"
#include "fatal_error.h"


/**
 * @brief Use at the end of a ptbc swap to change working directory to the new subdirectory
 *
 */
void ptbc_chdir_instance(void) {
  char subdir[1024];
  char logfile[1024];
  int local_rank;

  fflush(stdout);  /* flush into the OLD instance's log before moving away */

  /* release the old log; nothing is printed between here and the reopen below */
  if (freopen("/dev/null", "w", stdout) == NULL)
    fatal_error("could not detach stdout from the old instance logfile", "ptbc_chdir_instance");

  /* every rank has now closed its previous log, so no file has a writer left */
  MPI_Barrier(app()->mpi.world_comm);

  if (chdir("..") != 0)
    fatal_error("could not chdir back to the run directory", "ptbc_chdir_instance");

  snprintf(subdir, sizeof(subdir), "instance_%.2d", app()->ptbc.instance_id);
  if (chdir(subdir) != 0)
    fatal_error("could not chdir into the instance directory", "ptbc_chdir_instance");

  /* append: the log of this parameter set carries on from the group that was here before */
  MPI_Comm_rank(app()->mpi.comm, &local_rank);
  snprintf(logfile, sizeof(logfile), "hmc_rank%.2d.log", local_rank);
  if (freopen(logfile, "a", stdout) == NULL)
    fatal_error("could not reopen stdout on the instance logfile", "ptbc_chdir_instance");
  setvbuf(stdout, NULL, _IOLBF, 0);
}

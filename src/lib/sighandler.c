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

/************************************************************
 *
 * Routines to handle system signals
 *
 * void catch_ill_inst(int s)
 *
 * catches illegal instructions signal
 * and writes an error indication to
 * stdout.
 *
 * input:
 *  int s: signal number (not needed)
 *
 ************************************************************/

#ifdef HAVE_CONFIG_H
#include <tmlqcd_config.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#ifdef TM_USE_MPI
#include <mpi.h>
#endif

/* Catch an illegal instruction in order */
/* to give the user a hint what was wrong */
void catch_ill_inst(int s) {
  fprintf(stderr, "An illegal instruction occured!\n");
  fprintf(stderr, "Probably this caused the exception.\n");
  fprintf(stderr, "Aborting...\n");
  fflush(stdout);
#ifdef TM_USE_MPI
  MPI_Abort(MPI_COMM_WORLD, 1);
  MPI_Finalize();
#endif
  exit(0);
}

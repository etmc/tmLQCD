/***********************************************************************
 *  
 * Copyright (C) 2012 Carsten Urbach
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
#include <config.h>
#endif

#include <stdio.h>
#include <global.h>

#ifdef MPI
#include <mpi.h>
#endif

#include "fatal_error.h"

void fatal_error(char const *error, char const *function)
{
  if (error != NULL)
  {
    fprintf(stderr, "FATAL ERROR\n");
    if (function != NULL)
    {
#ifdef MPI
      fprintf(stderr, "  Within %s (reported by node %d):\n", function, g_proc_id);
#else
      fprintf(stderr, "  Within %s:\n", function);
#endif
    }
    fprintf(stderr, "    %s\n", error);
    fflush(stderr);
  }
  
#ifdef MPI
  MPI_Abort(MPI_COMM_WORLD, 1);
  MPI_Finalize();
#endif

  exit(500);
}

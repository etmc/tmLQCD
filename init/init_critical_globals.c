/***********************************************************************
 *
 * Copyright (C) 2018  Bartosz Kostrzewa
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

#include "misc_types.h"
#include "global.h"
#include "read_input.h"
#include "init/init_global_states.h"

void init_critical_globals(const tm_ProgramId_t program_id)
{
  verbose = 1;
  g_mpi_thread_level = TM_MPI_THREAD_SINGLE;
  g_timers.lvl = -1;
  
  init_global_states();

  /* further, program-specific initialisations might go here */
}


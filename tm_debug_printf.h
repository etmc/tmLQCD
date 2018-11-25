/***********************************************************************
 *  
 * Copyright (C) 2018 Bartosz Kostrzewa
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

#ifndef TM_DEBUG_PRINTF_H
#define TM_DEBUG_PRINTF_H

/* Function along the lines of printf which produces output on a single
 * or all MPI tasks (unordered) when g_debug_level is at or
 * above the provided threshold 
 * to have output by all MPI tasks, simply pass g_proc_id for proc_id */

void tm_debug_printf(const int proc_id,
                     const int dbg_level_threshold,
                     const char * format,
                     ...);

#endif

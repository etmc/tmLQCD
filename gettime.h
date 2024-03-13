/***********************************************************************                                                             
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008 Carsten Urbach
 * Copyright (C) 2012 Bartosz Kostrzewa (gettime.[c,h])
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

#ifndef _GETTIME_H
#define _GETTIME_H

#define TM_TIMING_MAX_LEVELS 15
#define TM_TIMING_STACK_PATH_LENGTH 500
#define TM_TIMING_NAME_LENGTH 50

/* gettime provides a time measurement with the BGL real time ticker,
   MPI_Wtime, clock_gettime and clock in decreasing order of preference
   depending on availability. Except for clock(), all these measurements
   are good representations of walltime */

#ifdef __cplusplus
extern "C" {
#endif

double gettime(void);

typedef struct tm_timers_s {
  int lvl;
  double t[TM_TIMING_MAX_LEVELS];
  char callstack[TM_TIMING_MAX_LEVELS][TM_TIMING_STACK_PATH_LENGTH];
  char name[TM_TIMING_MAX_LEVELS][TM_TIMING_NAME_LENGTH];
} tm_timers_t;

// tm_stopwatch_push will increase the timing level and perform a gettime()
// measurement
// the 'name' and context arguments are used to build a hierarchical tree
// where 'name' must always be specified and 'group' can be used
// for disambiguiation
// internally, if not at level 0, the context of the level below
// is prepended:
//
//   callstack[lvl-1]/[group:]name
//
void tm_stopwatch_push(tm_timers_t * const timers, const char * const name,
    const char * const group);

// tm_stopwatch_pop will output and decrease the timing level 
//
// # %s: Time for %s : %e s level: %d g_proc_id: %d %s \n
//
// with the prefix, the name and gettime()-startime inserted if the g_proc_id
// matches proc_id and the g_debug_level is equal or higher than
// dbg_level_threshold
// The call stack is given by the last field.
void tm_stopwatch_pop(tm_timers_t * const timers,
    const int proc_id, const int dbg_level_threshold,
    const char * const prefix);

#ifdef __cplusplus
}
#endif

#endif /* _GETTIME_H */


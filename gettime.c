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


#ifdef HAVE_CONFIG_H
# include<tmlqcd_config.h>
#endif
#ifdef HAVE_CLOCK_GETTIME
#  ifndef _POSIX_C_SOURCE
#    define _POSIX_C_SOURCE 199309L
#  endif
#  include <sys/time.h>
//#  include <bits/time.h>
#endif
#include <time.h>
#if (defined BGL && !defined BGP)
#  include <rts.h>
#endif
#ifdef TM_USE_MPI
#include <mpi.h>
#endif

#include "gettime.h"
#include "fatal_error.h"
#include "global.h"

#include <string.h>

double gettime(void) {
  double t;
#if (defined BGL && !defined BGP)

  const double clockspeed=1.0e-6/700.0;
  t = rts_get_timebase() * clockspeed;

#elif defined TM_USE_MPI

  t = MPI_Wtime();

  /* clock_gettime is detected on BGL/BGP but it is an unsupported system call so we can't use it! */
#elif (defined HAVE_CLOCK_GETTIME && !defined BGL)

  struct timespec ts;

  /* on the BGQ the monotonic clock is directly connected to the hardware counters
  and reports process CPU time, that is not a good measurement for threaded applications */
#  ifdef BGQ
  clock_gettime(CLOCK_REALTIME,&ts);
#  else
  clock_gettime(CLOCK_MONOTONIC,&ts);
#  endif
  t = ts.tv_sec + 1.0e-9*ts.tv_nsec;

#else
  /* This number is completely unreliable because the operating system and other processes
     make the clock tick too. This is especially true with multiple threads where the number
     of clock ticks will be multiplied by roughly the number of threads, but not quite, making
     the measurement useless!  */

  t = (double)clock()/(CLOCKS_PER_SEC);

#endif

  return t;
}

void tm_stopwatch_push(tm_timers_t * const timers, const char * const name, const char * group){
  if( timers->lvl+1 >= TM_TIMING_MAX_LEVELS ){
    fatal_error("attempted push would let timers->lvl exceed TM_TIMING_MAX_LEVELS", "tm_stopwatch_push");
  }
  timers->t[++(timers->lvl)] = gettime();
  if( strcmp(name, "") == 0 ){
    char msg[200];
    snprintf(msg, 200, "name may not be empty! (group: %s)", group);
    fatal_error(msg, __func__);
  } else {
    const int empty_group = strcmp(group, "") == 0;
    snprintf(timers->callstack[timers->lvl], TM_TIMING_STACK_PATH_LENGTH,
             "%s/%s%s%s",
             timers->lvl-1 >= 0 ? timers->callstack[timers->lvl-1] : "",
             empty_group ? "" : group,
             empty_group ? "" : ":",
             name);
    snprintf(timers->name[timers->lvl], TM_TIMING_NAME_LENGTH, "%s", name);
  }
}

void tm_stopwatch_pop(tm_timers_t * const timers,
    const int proc_id, const int dbg_level_threshold,
    const char * const prefix){
  if( (timers->lvl-1) < -1 ){
    fatal_error("attempted pop would let timers->lvl drop below the lowest possible value", "tm_stopwatch_pop");
  }
  if( g_proc_id == proc_id && g_debug_level >= dbg_level_threshold ){
    printf("# %s: Time for %s %e s level: %d proc_id: %d %s\n", 
           prefix,
           timers->name[timers->lvl], 
           gettime()-timers->t[timers->lvl], 
           timers->lvl, g_proc_id,
           timers->callstack[timers->lvl]);
  }
  (timers->lvl)--;
}


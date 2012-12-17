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
# include<config.h>
#endif
#include <time.h>
#ifdef HAVE_CLOCK_GETTIME
#include <sys/time.h>
#include <bits/time.h>
#endif
#if (defined BGL && !defined BGP)
#  include <rts.h>
#endif
#ifdef MPI
# include <mpi.h>
#endif

#include "gettime.h"

double gettime(void) {
  double t;
#if (defined BGL && !defined BGP)

  const double clockspeed=1.0e-6/700.0;
  t = rts_get_timebase() * clockspeed;

#elif defined MPI

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


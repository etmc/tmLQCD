/***********************************************************************
 * Copyright (C) 2013 Bartosz Kostrzewa
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
#ifdef OMP
#include <omp.h>
#include "init_omp_accumulators.h"
#endif
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include "global.h"

void init_openmp(void) {
#ifdef OMP  
  if(omp_num_threads > 0) 
  {
     omp_set_num_threads(omp_num_threads);
     if( g_debug_level > 0 && g_proc_id == 0 ) {
       printf("# Instructing OpenMP to use %d threads.\n",omp_num_threads);
     }
  }
  else {
    if( g_proc_id == 0 )
      printf("# No value provided for OmpNumThreads, running in single-threaded mode!\n");

    omp_num_threads = 1;
    omp_set_num_threads(omp_num_threads);
  }

  init_omp_accumulators(omp_num_threads);
#endif
  return;
}


/***********************************************************************
 * Copyright (C) 2012 Bartosz Kostrzewa
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
#endif
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include "global.h"
#include "init_omp_kahan_arrays.h"

int init_omp_kahan_arrays(const int num) {
  g_omp_ks=NULL;
  g_omp_kc=NULL;

  if((void*)(g_omp_ks = (double*)calloc(num, sizeof(double))) == NULL) {
    printf ("init_omp_kahan_arrays malloc errno : %d\n",errno); 
    errno = 0;
    return(1);
  }
  if((void*)(g_omp_kc = (double*)calloc(num, sizeof(double))) == NULL) {
    printf ("init_omp_kahan_arrays malloc errno : %d\n",errno); 
    errno = 0;
    return(2);
  }

  return(0);
}

void free_omp_kahan_arrays() {
  free(g_omp_kc);
  free(g_omp_ks);
}

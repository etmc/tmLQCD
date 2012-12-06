/***********************************************************************
 *
 * Copyright (C) 1995 Ulli Wolff, Stefan Sint
 *               2001,2005 Martin Hasenbusch
 *               2011,2012 Carsten Urbach
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
#ifdef SSE
# undef SSE
#endif
#ifdef SSE2
# undef SSE2
#endif
#ifdef SSE3
# undef SSE3
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <time.h>
#ifdef MPI
# include <mpi.h>
#endif
#ifdef OMP
# include <omp.h>
#endif
#include "global.h"
#include "su3.h"
#include "sse.h"
#include "su3adj.h"
#include "operator/clovertm_operators.h"
#include "operator/clover_leaf.h"

const double tiny_t = 1.0e-20;

su3 ** swm, ** swp;

void mult_6x6(_Complex double a[6][6], _Complex double b[6][6], _Complex double d[6][6]) {

  for(int i = 0; i < 6; i++) {
    for(int j = 0; j < 6; j++) {
      a[i][j] = 0.;
      for(int k = 0; k < 6; k++) {
	a[i][j] += b[i][k] * d[k][j];
      }
    }
  }
  return;
}

void add_6x6(_Complex double a[6][6], _Complex double b[6][6], _Complex double d[6][6]) {

  for(int i = 0; i < 6; i++) {
    for(int j = 0; j < 6; j++) {
      a[i][j] = b[i][j] + d[i][j];
    }
  }
  return;
}

void sub_6x6(_Complex double a[6][6], _Complex double b[6][6], _Complex double d[6][6]) {

  for(int i = 0; i < 6; i++) {
    for(int j = 0; j < 6; j++) {
      a[i][j] = b[i][j] - d[i][j];
    }
  }
  return;
}

void copy_6x6(_Complex double a[6][6], const _Complex double b[6][6]) {
  for(int i = 0; i < 6; i++) {
    for(int j = 0; j < 6; j++) {
      a[i][j] = b[i][j];
    }
  }
  return;
}








su3 * _swp;

int init_swpm(const int V) {
  int i=0;
  static int swpm_init=0;

  if(!swpm_init) {
    if((void*)(swp = (su3**)calloc(V, sizeof(su3*))) == NULL) {
      printf ("malloc errno : %d\n",errno); 
      errno = 0;
      return(1);
    }
    if((void*)(swm = (su3**)calloc(V, sizeof(su3*))) == NULL) {
      printf ("malloc errno : %d\n",errno); 
      errno = 0;
      return(1);
    }
    if((void*)(_swp = (su3*)calloc(2*4*V+1, sizeof(su3))) == NULL) {
      printf ("malloc errno : %d\n",errno); 
      errno = 0;
      return(2);
    }
#if (defined SSE || defined SSE2 || defined SSE3)
    swp[0] = (su3*)(((unsigned long int)(_swp)+ALIGN_BASE)&~ALIGN_BASE);
#else
    swp[0] = _swp;
#endif
    swm[0] = swp[0] + 4*V;
    for(i = 1; i < V; i++){
      swp[i] = swp[i-1]+4;
      swm[i] = swm[i-1]+4;
    }
    swpm_init = 1;
  }
  return(0);
}

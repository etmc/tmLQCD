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

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include "global.h"
#include "su3.h"
#include "sse.h"
#include "init_gauge_tmp.h"

su3 * gauge_tmp_ = NULL;
su3 ** gauge_tmp = NULL;

int init_gauge_tmp(const int V) {
  int i=0;

  if((void*)(gauge_tmp = (su3**)calloc(V, sizeof(su3*))) == NULL) {
    fprintf(stderr, "malloc errno : %d\n", errno);
    errno = 0;
    return(1);
  }
  if((void*)(gauge_tmp_ = (su3*)calloc(4*V+1, sizeof(su3))) == NULL) {
    fprintf(stderr, "malloc errno : %d\n", errno);
    errno = 0;
    return(1);
  }
#if (defined SSE || defined SSE2 || defined SSE3)
  gauge_tmp[0] = (su3*)(((unsigned long int)(gauge_tmp_)+ALIGN_BASE)&~ALIGN_BASE);
#else
  gauge_tmp[0] = gauge_tmp_;
#endif
  for(i = 1; i < V; i++){
    gauge_tmp[i] = gauge_tmp[i-1]+4;
  }
  return(0);
}

void free_gauge_tmp() {
  free(gauge_tmp_);
  free(gauge_tmp);
}

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
#include "init_gauge_fg.h"

su3 * gauge_fg_ = NULL;
su3 ** gauge_fg = NULL;

int init_gauge_fg(const int V) {
  int i=0;

  if((void*)(gauge_fg = (su3**)calloc(V, sizeof(su3*))) == NULL) {
    fprintf(stderr, "malloc errno : %d\n", errno);
    errno = 0;
    return(1);
  }
  if((void*)(gauge_fg_ = (su3*)calloc(4*V+1, sizeof(su3))) == NULL) {
    fprintf(stderr, "malloc errno : %d\n", errno);
    errno = 0;
    return(1);
  }
#if (defined SSE || defined SSE2 || defined SSE3)
  gauge_fg[0] = (su3*)(((unsigned long int)(gauge_fg_)+ALIGN_BASE)&~ALIGN_BASE);
#else
  gauge_fg[0] = gauge_fg_;
#endif
  for(i = 1; i < V; i++){
    gauge_fg[i] = gauge_fg[i-1]+4;
  }
  return(0);
}

void free_gauge_fg() {
  free(gauge_fg_);
  free(gauge_fg);
}

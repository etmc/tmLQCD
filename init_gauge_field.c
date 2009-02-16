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
/* $Id$ */

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include "global.h"
#include "su3.h"
#include "sse.h"
#include "init_gauge_field.h"

su3 * gauge_field = NULL;
su3 * gauge_field_copy = NULL;

int init_gauge_field(const int V, const int back) {
  int i=0;

  g_gauge_field_copy = NULL;
  if((void*)(g_gauge_field = (su3**)calloc(V, sizeof(su3*))) == NULL) {
    printf ("malloc errno : %d\n",errno); 
    errno = 0;
    return(1);
  }
  if((void*)(gauge_field = (su3*)calloc(4*V+1, sizeof(su3))) == NULL) {
    printf ("malloc errno : %d\n",errno); 
    errno = 0;
    return(2);
  }
#if (defined SSE || defined SSE2 || defined SSE3)
  g_gauge_field[0] = (su3*)(((unsigned long int)(gauge_field)+ALIGN_BASE)&~ALIGN_BASE);
#else
  g_gauge_field[0] = gauge_field;
#endif
  for(i = 1; i < V; i++){
    g_gauge_field[i] = g_gauge_field[i-1]+4;
  }

#  if defined _USE_HALFSPINOR
  if(back == 1) {
    /*
      g_gauge_field_copy[ieo][PM][sites/2][mu]
    */
    if((void*)(g_gauge_field_copy = (su3***)calloc(2, sizeof(su3**))) == NULL) {
      printf ("malloc errno : %d\n",errno); 
      errno = 0;
      return(3);
    }
    if((void*)(g_gauge_field_copy[0] = (su3**)calloc(VOLUME, sizeof(su3*))) == NULL) {
      printf ("malloc errno : %d\n",errno); 
      errno = 0;
      return(3);
    }
    g_gauge_field_copy[1] = g_gauge_field_copy[0] + (VOLUME)/2;
    if((void*)(gauge_field_copy = (su3*)calloc(4*(VOLUME)+1, sizeof(su3))) == NULL) {
      printf ("malloc errno : %d\n",errno); 
      errno = 0;
      return(4);
    }
#    if (defined SSE || defined SSE2 || defined SSE3)
    g_gauge_field_copy[0][0] = (su3*)(((unsigned long int)(gauge_field_copy)+ALIGN_BASE)&~ALIGN_BASE);
#    else
    g_gauge_field_copy[0][0] = gauge_field_copy;
#    endif
    for(i = 1; i < (VOLUME)/2; i++) {
      g_gauge_field_copy[0][i] = g_gauge_field_copy[0][i-1]+4;
    }
    g_gauge_field_copy[1][0] = g_gauge_field_copy[0][0] + 2*VOLUME; 
    for(i = 1; i < (VOLUME)/2; i++) {
      g_gauge_field_copy[1][i] = g_gauge_field_copy[1][i-1]+4;
    }
  }
#  else
  if(back == 1) {
    if((void*)(g_gauge_field_copy = (su3**)calloc((VOLUME+RAND), sizeof(su3*))) == NULL) {
      printf ("malloc errno : %d\n",errno); 
      errno = 0;
      return(3);
    }
    if((void*)(gauge_field_copy = (su3*)calloc(8*(VOLUME+RAND)+1, sizeof(su3))) == NULL) {
      printf ("malloc errno : %d\n",errno); 
      errno = 0;
      return(4);
    }
#  if (defined SSE || defined SSE2 || defined SSE3)
    g_gauge_field_copy[0] = (su3*)(((unsigned long int)(gauge_field_copy)+ALIGN_BASE)&~ALIGN_BASE);
#  else
    g_gauge_field_copy[0] = gauge_field_copy;
#  endif
    for(i = 1; i < (VOLUME+RAND); i++) {
      g_gauge_field_copy[i] = g_gauge_field_copy[i-1]+8;
    }
  }
#  endif
  g_update_gauge_copy = 1;
  g_update_gauge_energy = 1;
  g_update_rectangle_energy = 1;
  return(0);
}

void free_gauge_field() {
  free(gauge_field);
  free(g_gauge_field);
  free(gauge_field_copy);
}

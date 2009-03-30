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

spinor * sp_cup = NULL;
spinor * sp_cdn = NULL;

int init_chi_up_copy(const int V) {

  if((void*)(sp_cup = (spinor*)calloc(V+1, sizeof(spinor))) == NULL) {
    printf ("malloc errno : %d\n",errno); 
    errno = 0;
    return(1);
  }
  if((void*)(g_chi_up_copy = malloc(sizeof(spinor*))) == NULL) {
    printf ("malloc errno : %d\n",errno); 
    errno = 0;
    return(2);
  }
#if ( defined SSE || defined SSE2 || defined SSE3)
  g_chi_up_copy = (spinor*)(((unsigned long int)(sp_cup)+ALIGN_BASE)&~ALIGN_BASE);
#else
  g_chi_up_copy = sp_cup;
#endif
  
  return(0);
}

void free_chi_up_copy() {

  free(sp_cup);
}


int init_chi_dn_copy(const int V) {

  if((void*)(sp_cdn = (spinor*)calloc(V+1, sizeof(spinor))) == NULL) {
    printf ("malloc errno : %d\n",errno); 
    errno = 0;
    return(1);
  }
  if((void*)(g_chi_dn_copy = malloc(sizeof(spinor*))) == NULL) {
    printf ("malloc errno : %d\n",errno); 
    errno = 0;
    return(2);
  }
#if ( defined SSE || defined SSE2 || defined SSE3)
  g_chi_dn_copy = (spinor*)(((unsigned long int)(sp_cdn)+ALIGN_BASE)&~ALIGN_BASE);
#else
  g_chi_dn_copy = sp_cdn;
#endif

  return(0);
}

void free_chi_dn_copy() {

  free(sp_cdn);
}




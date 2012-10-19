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
#ifdef _USE_SHMEM
# include <mpp/shmem.h>
#endif
#include "global.h"
#include "su3.h"
#include "sse.h"
#include "monomial/monomial.h"

spinor * sp = NULL;
spinor * sp_csg = NULL;
spinor * sp_tbuff = NULL;

int init_spinor_field(const int V, const int nr) {
  int i = 0;

#if (defined _USE_SHMEM && !(defined _USE_HALFSPINOR))
  if((void*)(sp = (spinor*)shmalloc((nr*V+1)*sizeof(spinor))) == NULL) {
    printf ("malloc errno : %d\n",errno); 
    errno = 0;
    return(1);
  }
#else
  if((void*)(sp = (spinor*)calloc(nr*V+1, sizeof(spinor))) == NULL) {
    printf ("malloc errno : %d\n",errno); 
    errno = 0;
    return(1);
  }
#endif
  if((void*)(g_spinor_field = (spinor**)malloc(nr*sizeof(spinor*))) == NULL) {
    printf ("malloc errno : %d\n",errno); 
    errno = 0;
    return(2);
  }
#if ( defined SSE || defined SSE2 || defined SSE3)
  g_spinor_field[0] = (spinor*)(((unsigned long int)(sp)+ALIGN_BASE)&~ALIGN_BASE);
#else
  g_spinor_field[0] = sp;
#endif
  
  for(i = 1; i < nr; i++){
    g_spinor_field[i] = g_spinor_field[i-1]+V;
  }

  return(0);
}

void free_spinor_field() {
#if (defined _USE_SHMEM && !(defined _USE_HALFSPINOR))
  shfree(sp);
  shfree(sp_csg);
#else
  free(sp);
  free(sp_csg);
#endif
}


/** 
 * costumized spinor allocation routines
 */
int allocate_spinor_field_array(spinor ***spinors,spinor **sp,const int V, const int nr) {
  int i = 0;

#if (defined _USE_SHMEM && !(defined _USE_HALFSPINOR))
  if((void*)((*sp) = (spinor*)shmalloc((nr*V+1)*sizeof(spinor))) == NULL) {
    printf ("malloc errno : %d\n",errno); 
    errno = 0;
    return(1);
  }
#else
  if((void*)((*sp) = (spinor*)calloc(nr*V+1, sizeof(spinor))) == NULL) {
    printf ("malloc errno : %d\n",errno); 
    errno = 0;
    return(1);
  }
#endif
  if((void*)((*spinors) = (spinor**)malloc(nr*sizeof(spinor*))) == NULL) {
    printf ("malloc errno : %d\n",errno); 
    errno = 0;
    return(2);
  }
#if ( defined SSE || defined SSE2 || defined SSE3)
  (*spinors)[0] = (spinor*)(((unsigned long int)(*sp)+ALIGN_BASE)&~ALIGN_BASE);
#else
  (*spinors)[0] = *sp;
#endif
  
  for(i = 1; i < nr; i++){
    (*spinors)[i] = (*spinors)[i-1]+V;
  }

  return(0);
}

void free_spinor_field_array(spinor** sp) {
#if (defined _USE_SHMEM && !(defined _USE_HALFSPINOR))
  shfree(*sp);
#else
  free(*sp);
#endif
}



#ifndef _BENCH_ONLY
int init_csg_field(const int V) {
  int i = 0, j = 0, sum = 0;
  spinor * s;
  for(i = 0; i < no_monomials; i++) {
    sum += monomial_list[i].csg_N;
    sum += monomial_list[i].csg_N2;
  }

  /* if all histories are zero, we do not need initialisation */
  if(sum != 0) {
#if (defined _USE_SHMEM && !(defined _USE_HALFSPINOR))
    sp_csg = (spinor*)shmalloc((sum*V+1)*sizeof(spinor));
#else
    sp_csg = (spinor*)calloc(sum*V+1, sizeof(spinor));
#endif
    if(errno == ENOMEM) {
      return(1);
    }
    for(i = 0; i < no_monomials; i++) {
      monomial_list[i].csg_field = malloc((monomial_list[i].csg_N+1)*sizeof(spinor*));
      if(errno == ENOMEM) {
	return(2);
      }
      monomial_list[i].csg_field2 = malloc(monomial_list[i].csg_N2*sizeof(spinor*));
      if(errno == ENOMEM) {
	return(2);
      }
    }
#if ( defined SSE || defined SSE2 || defined SSE3)
    s = (spinor*)(((unsigned long int)(sp_csg)+ALIGN_BASE)&~ALIGN_BASE);
#else
    s = sp_csg;
#endif
    for(j = 0; j < no_monomials; j++) {
      if(monomial_list[j].csg_N != 0) {
	for(i = 0; i < monomial_list[j].csg_N; i++) {
	  monomial_list[j].csg_field[i] = s;
	  s = s + V;
	}
      }
    }
    for(j = 0; j < no_monomials; j++) {
      if(monomial_list[j].csg_N2 != 0) {
	for(i = 0; i < monomial_list[j].csg_N2; i++) {
	  monomial_list[j].csg_field2[i] = s;
	  s = s + V;
	}
      }
    }
    
    monomial_list[0].csg_index_array = (int*) malloc(sum*sizeof(int));
    for(i = 1; i < no_monomials; i++) {
      monomial_list[i].csg_index_array = monomial_list[i-1].csg_index_array + monomial_list[i-1].csg_N;
    }
    monomial_list[0].csg_index_array2 = monomial_list[no_monomials-1].csg_index_array 
      + monomial_list[no_monomials-1].csg_N;
    for(i = 1; i < no_monomials; i++) {
      monomial_list[i].csg_index_array2 = monomial_list[i-1].csg_index_array2 + monomial_list[i-1].csg_N2;
    }
  }
  return(0);
}

#endif

int init_timslice_buffer_field(const int t_slice) {
  
  if((void*)(sp_tbuff = (spinor*)calloc(t_slice+1, sizeof(spinor))) == NULL) {
    printf ("malloc errno : %d\n",errno);
    errno = 0;
    return(3);
  }

#if (( defined SSE || defined SSE2 || defined SSE3) && defined _USE_TSPLITPAR )
  g_tbuff = (spinor*)(((unsigned long int)(sp_tbuff)+ALIGN_BASE)&~ALIGN_BASE);
#else 
  g_tbuff = sp_tbuff;
#endif
  
  return(0);
}

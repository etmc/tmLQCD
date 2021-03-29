/***********************************************************************                                                             
 * Copyright (C) 2015 Bartosz Kostrzewa
 *               2016 Carsten Urbach
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

#if HAVE_CONFIG_H
#include "tmlqcd_config.h"
#endif

#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include "aligned_malloc.h"
#include "su3.h"
#include "su3adj.h"

#include "fatal_error.h"
 
void *aligned_malloc(size_t const size) {
  void *mem = malloc(size+ALIGN_BASE+sizeof(void*));
  void ** ptr;

  if(mem == NULL) {
    return(mem);
  }

  ptr = (void**)(((uintptr_t)mem+(uintptr_t)ALIGN_BASE+sizeof(void*)) & ~ (uintptr_t)(ALIGN_BASE));
  ptr[-1] = mem;
  
  return ptr;
}

void *aligned_malloc_zero(size_t const size) {
  void *mem = malloc(size+ALIGN_BASE+sizeof(void*));
  void ** ptr;

  if(mem == NULL) {
    return(mem);
  }

  ptr = (void**)(((uintptr_t)mem+(uintptr_t)ALIGN_BASE+sizeof(void*)) & ~ (uintptr_t)(ALIGN_BASE));
  ptr[-1] = mem;
  memset(ptr, 0, size);

  return ptr;
}


void aligned_free(void *ptr) {
  free(((void**)ptr)[-1]);
}

aligned_su3_field_t aligned_su3_field_alloc(const unsigned int V) {
  aligned_su3_field_t f_struct;

  su3** field = (su3**) aligned_malloc(V*sizeof(su3*));
  su3* mem = (su3*)aligned_malloc((4*V+1)*sizeof(su3));

  if( (void*)field == (void*)NULL || (void*)mem == (void*)NULL ) {
    fatal_error("Memory allocation error!","aligned_su3_field_alloc");
  }

  field[0] = mem;
  for(int i = 1; i < V; ++i) {
    field[i] = field[i-1]+4;
  }

  f_struct.field = field;
  f_struct.mem = mem;

  return(f_struct);
}

aligned_su3adj_field_t aligned_su3adj_field_alloc(const unsigned int V) {
  aligned_su3adj_field_t f_struct;
  su3adj** field = (su3adj**) aligned_malloc(V*sizeof(su3adj*));
  su3adj* mem = (su3adj*)aligned_malloc((4*V+1)*sizeof(su3adj));

  if( (void*)field == (void*)NULL || (void*)mem == (void*)NULL ) {
    fatal_error("Memory allocation error!","aligned_su3_field_alloc");
  }

  field[0] = mem;
  for(int i = 1; i < V; ++i) {
    field[i] = field[i-1]+4;
  }

  f_struct.field = field;
  f_struct.mem = mem;

  return(f_struct);
}

void aligned_su3_field_free(const aligned_su3_field_t* f_struct) {
  aligned_free((void*)f_struct->field);
  aligned_free((void*)f_struct->mem);
}

void aligned_su3adj_field_free(const aligned_su3adj_field_t* f_struct) {
  aligned_free((void*)f_struct->field);
  aligned_free((void*)f_struct->mem);
}


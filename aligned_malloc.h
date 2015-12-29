/***********************************************************************                                                                                     
 * Copyright (C) 2015 Bartosz Kostrzewa
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

#ifndef _ALIGNED_MALLOC_H
#define _ALIGNED_MALLOC_H

#include "su3.h"
#include "su3adj.h"

typedef struct {
  su3** field;
  su3* mem;
} aligned_su3_field_t;

typedef struct {
  su3adj** field;
  su3adj* mem;
} aligned_su3adj_field_t;

aligned_su3_field_t aligned_su3_field_alloc(const unsigned int V);
aligned_su3adj_field_t aligned_su3adj_field_alloc(const unsigned int V); 

void aligned_su3_field_free(const aligned_su3_field_t* f_struct);
void aligned_su3adj_field_free(const aligned_su3adj_field_t* f_stuct);

#endif

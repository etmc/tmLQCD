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

#ifndef _ASSIGN_H
#define _ASSIGN_H

#include "su3.h"

/* Assign (*R) = (*S) */
void assign(spinor * const R, spinor * const S, const int N);
void assign_32(spinor32 * const R, spinor32 * const S, const int N);
void assign_su3vect(su3_vector * const R, su3_vector * const S, const int N);

void assign_complex_to_spinor(spinor * const R, _Complex double * const S, const int N); /* N is the size of S */
void assign_spinor_to_complex(_Complex double * const R, spinor * const S, const int N); /* N is the size of S */

void assign_complex_to_spinor_32(spinor32 * const R, _Complex float* const S, const int N);
void assign_spinor_to_complex_32(_Complex float* const R, spinor32 * const S, const int N); /* N is the size of S */

#endif

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

/************************************************
 *
 * Typedefinition of the pointer to the function
 * which contains the matrix multiplication.
 *
 ************************************************/

#ifndef _MATRIX_MULT_TYPEDEF_H
#define _MATRIX_MULT_TYPEDEF_H

typedef void (*matrix_mult)(spinor * const, spinor * const);
typedef void (*matrix_mult_blk)(spinor * const, spinor * const, const int);
typedef void (*matrix_mult_clover)(spinor * const, spinor * const, const double);
typedef void (*c_matrix_mult)(_Complex double * const, _Complex double * const);
typedef void (*matrix_mult_su3vect)(su3_vector * const, su3_vector * const, const int);

#endif

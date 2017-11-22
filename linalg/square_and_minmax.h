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

#ifndef _SQUARE_AND_MAX_H
#define _SQUARE_AND_MAX_H

#include "su3.h"

/* double square_and_minmax(spinor * const P )
 *     Returns the square norm of *P and the local minimal/maximal norm */

/* double square_and_minmax(spinor * const P, spinor * const Q )
 *     Returns the square norm of *P/*Q (locally) and the local minimal/maximal norm */

void square_and_minmax(double * const sum, double * const min, double * const max, const spinor * const P, const int N);
void square_and_minmax_rel(double * const sum, double * const min, double * const max, const spinor * const P,  const spinor * const Q, const int N);
void square_and_minmax_abs(double * const sum, double * const min, double * const max, double * const min_abs, double * const max_abs, const spinor * const P, const int N);
void square_and_minmax_rel_abs(double * const sum, double * const min, double * const max, double * const min_abs, double * const max_abs, const spinor * const P,  const spinor * const Q, const int N);



#endif




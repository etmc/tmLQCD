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
#ifndef _CHEBYSHEV_POLYNOMIAL_ND_H
#define _CHEBYSHEV_POLYNOMIAL_ND_H

#include "solver/matrix_mult_typedef_nd.h"

double func(double u, double exponent);

void chebyshev_coefs(double a, double b, double c[], int n, double exponent);

double cheb_eval(int M, double *c, double s);

void degree_of_polynomial_nd(int * _degree_of_p, double ** coefs,
			     const double EVMin, const double EVMax,
			     matrix_mult_nd Qsq, const int repro);

#endif

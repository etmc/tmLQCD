/***********************************************************************
 *
 * Copyright (C) 2006,2007,2008 Karl Jansen, Carsten Urbach 
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

#ifndef _CHEBYSHEV_POLYNOMIAL_H
#define _CHEBYSHEV_POLYNOMIAL_H

extern double cheb_evmin, cheb_evmax;
extern int dop_n_cheby;
extern double * dop_cheby_coef;


double func(double u, double exponent);
void chebyshev_polynomial(double a, double b, double c[], int n, double exponent);

void QdaggerQ_power(spinor *R_s, spinor *R_c, double *c, int n, spinor *S_s, spinor *S_c);

void degree_of_polynomial(const int repro);

#endif

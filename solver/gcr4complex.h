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

#ifndef _GCR4COMPLEX_H
#define _GCR4COMPLEX_H

#include"solver/matrix_mult_typedef.h"
#include"su3.h"

void ldiff(complex * Q, complex * const R, complex * const S, const int N);
void ladd(complex * Q, complex * const R, complex * const S, const int N);
double lsquare_norm(complex * const Q, const int N, const int parallel);
complex lscalar_prod(complex * const R, complex * const S, const int N, const int parallel);
void lmul_r(complex * const R, const double c, complex * const S, const int N);
void lmul(complex * const R, const complex c, complex * const S, const int N);
void lassign_diff_mul(complex * const R, complex * const S, const complex c, const int N);
void lassign_add_mul(complex * const R, complex * const S, const complex c, const int N);
void ldiff_assign(complex * const Q, complex * const S, 
		  const int N);
void ladd_assign(complex * const Q, complex * const S, 
		  const int N);


int gcr4complex(complex * const P, complex * const Q, 
		const int m, const int max_restarts,
		const double eps_sq, const int rel_prec,
		const int N, const int parallel,
		const int lda, c_matrix_mult f);


#endif

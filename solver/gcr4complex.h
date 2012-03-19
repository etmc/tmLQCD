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

void ldiff(_Complex double * Q, _Complex double * const R, _Complex double * const S, const int N);
void ladd(_Complex double * Q, _Complex double * const R, _Complex double * const S, const int N);
double lsquare_norm(_Complex double * const Q, const int N, const int parallel);
_Complex double lscalar_prod(_Complex double * const R, _Complex double * const S, const int N, const int parallel);
void lmul_r(_Complex double * const R, const double c, _Complex double * const S, const int N);
void lmul(_Complex double * const R, const _Complex double c, _Complex double * const S, const int N);
void lassign_diff_mul(_Complex double * const R, _Complex double * const S, const _Complex double c, const int N);
void lassign_add_mul(_Complex double * const R, _Complex double * const S, const _Complex double c, const int N);
void ldiff_assign(_Complex double * const Q, _Complex double * const S, 
		  const int N);
void ladd_assign(_Complex double * const Q, _Complex double * const S, 
		  const int N);


int gcr4complex(_Complex double * const P, _Complex double * const Q, 
		const int m, const int max_restarts,
		const double eps_sq, const int rel_prec,
		const int N, const int parallel,
		const int lda, c_matrix_mult f);


#endif

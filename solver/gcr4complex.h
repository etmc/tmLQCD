/* $Id$ */

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

/* $Id$ */
#ifndef _LITTLE_D_H
#define _LITTLE_D_H

#include "complex.h"

typedef void (*c_matrix_mult) (complex * const, complex * const);
void little_D(complex * v, complex *w);
void free_lgcr();
int lgcr(complex * const P, complex * const Q, 
	 const int m, const int max_restarts,
	 const double eps_sq, const int rel_prec,
	 const int N, const int lda, c_matrix_mult f);

#endif

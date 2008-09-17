/* $Id$ */

#ifndef _GCR4COMPLEX_H
#define _GCR4COMPLEX_H

#include"solver/matrix_mult_typedef.h"
#include"su3.h"

int gcr4complex(complex * const P, complex * const Q, 
		const int m, const int max_restarts,
		const double eps_sq, const int rel_prec,
		const int N, const int parallel,
		const int lda, c_matrix_mult f);


#endif

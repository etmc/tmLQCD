#ifndef _MR4COMPLEX_H
#define _MR4COMPLEX_H

#include "solver/matrix_mult_typedef.h"
#include "solver/gcr4complex.h"

#include <complex.h>

int mr4complex(_Complex double * const P, _Complex double * const Q,
	       const int max_iter, const double eps_sq,
	       const int rel_prec, const int N, 
	       const int parallel, const int lda,
	       c_matrix_mult f);

#endif

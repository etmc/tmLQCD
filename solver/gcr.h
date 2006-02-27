/* $Id$ */

#ifndef _GCR_H
#define _GCR_H

#include"solver/matrix_mult_typedef.h"
#include"su3.h"

int gcr(spinor * const P, spinor * const Q, 
	const int m, const int max_restarts,
	const double eps_sq, const int rel_prec,
	const int N, matrix_mult f);

#endif

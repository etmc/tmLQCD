#ifndef _CG_HER_ND_H
#define _CG_HER_ND_H

#include"solver/matrix_mult_typedef_nd.h"
#include"su3.h"

int cg_her_nd(spinor * const, spinor * const,spinor * const, spinor * const, const int max_iter, double eps_sq, const int rel_prec,
	   const int N, matrix_mult_nd f, const int, const int);

#endif

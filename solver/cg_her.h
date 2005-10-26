#ifndef _CG_HER_H
#define _CG_HER_H

#include"solver/matrix_mult_typedef.h"
#include"su3.h"

int cg_her(spinor * const, spinor * const, const int max_iter, double eps_sq, const int rel_prec,
	   const int N, matrix_mult f, const int, const int);

#endif

#ifndef _PCG_HER_H
#define _PCG_HER_H

#include"solver/matrix_mult_typedef.h"
#include"su3.h"

int pcg_her(spinor * const, spinor * const, const int max_iter, double eps_sq, const int rel_prec,
	   const int N, matrix_mult f);

#endif

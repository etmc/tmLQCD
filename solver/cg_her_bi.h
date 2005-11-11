#ifndef _CG_HER_BI_H
#define _CG_HER_BI_H

#include"solver/matrix_mult_typedef.h"
#include"su3.h"

#include"solver/matrix_mult_typedef_bi.h"

int cg_her_bi(bispinor * const, bispinor * const, const int max_iter, double eps_sq, const int rel_prec, const int N, matrix_mult_bi f, const int, const int);

#endif

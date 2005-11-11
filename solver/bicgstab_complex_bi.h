/* $Id$ */

#ifndef _BICGSTAB_COMPLEX_BI_H
#define _BICGSTAB_COMPLEX_BI_H

#include"solver/matrix_mult_typedef.h"
#include"su3.h"

#include"solver/matrix_mult_typedef_bi.h"

int bicgstab_complex_bi(bispinor * const, bispinor * const, const int max_iter, double eps_sq, const int rel_prec, const int N, matrix_mult_bi f);

#endif

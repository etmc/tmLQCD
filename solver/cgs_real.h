

#ifndef _CGS_REAL_H
#define _CGS_REAL_H

#include"solver/matrix_mult_typedef.h"
#include"su3.h"

int cgs_real(spinor * const, spinor * const, const int max_iter, double eps_sq, matrix_mult f);

#endif

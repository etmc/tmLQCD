/* $Id$ */

#ifndef _BICGSTAB_COMPLEX_H
#define _BICGSTAB_COMPLEX_H

#include"solver/matrix_mult_typedef.h"
#include"su3.h"

int bicgstab_complex(spinor * const, spinor * const, const int max_iter, double eps_sq, matrix_mult f);

#endif

/* $Id$ */
#ifndef _LITTLE_D_H
#define _LITTLE_D_H

#include "complex.h"

typedef void (*c_matrix_mult) (complex * const, complex * const);
void little_D(complex * v, complex *w);

#endif

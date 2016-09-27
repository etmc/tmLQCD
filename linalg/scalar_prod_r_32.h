#ifndef _SCALAR_PROD_R_32_H
#define _SCALAR_PROD_R_32_H

#include "su3.h"

/* Returns the real part of the scalar product (*R,*S) */
float scalar_prod_r_32(const spinor32 * const S, const spinor32 * const R, const int N, const int parallel);

#endif
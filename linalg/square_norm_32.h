#ifndef _SQUARE_NORM_32_H
#define _SQUARE_NORM_32_H

#include "su3.h"

/* double square_norm(spinor * const P )
 *     Returns the square norm of *P */

float square_norm_32(const spinor32 * const P, const int N, const int parallel);
float square_norm_ts_32(const spinor32 * const P, const int N, const int parallel);
#endif

/* $Id$  */

#ifndef _SCALAR_PROD_I_H
#define _SCALAR_PROD_I_H

#include "su3.h"

/* Returns the imaginary part of the scalar product (*R,*S) */
double scalar_prod_i(spinor * const S,spinor * const R, const int N, const int parallel);

#endif

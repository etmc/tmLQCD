/* $Id$  */

#ifndef _SCALAR_PROD_R_H
#define _SCALAR_PROD_R_H

#include "su3.h"

/* Returns the real part of the scalar product (*R,*S) */
double scalar_prod_r(spinor * const S,spinor * const R, const int N, const int parallel);

#endif

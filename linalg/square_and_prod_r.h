/* $Id$  */

#ifndef _SQUARE_AND_PROD_R_H
#define _SQUARE_AND_PROD_R_H

#include "su3.h"

/*  Returns the real part of (*R,*S) and the square norm of *S 
 *  It's faster than using "scalar_prod_r" and "square_norm"  */
void square_and_prod_r(double * const x1, double * const x2, spinor * const S, spinor * const R, const int N);

#endif

/* $Id$  */

#ifndef _DIFF_H
#define _DIFF_H

#include "su3.h"

/* Makes the difference (*Q) = (*R) - (*S) */
void diff(spinor * const Q, spinor * const R, spinor * const S, const int N);


#endif

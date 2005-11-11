/* $Id$  */

#ifndef _DIFF_BI_H
#define _DIFF_BI_H

#include "su3.h"

/* Makes the difference (*Q) = (*R) - (*S) */
void diff_bi(bispinor * const Q, bispinor * const R, bispinor * const S, const int N);


#endif

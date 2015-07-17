#ifndef _ADDTO_32_H
#define _ADDTO_32_H

#include "su3.h"

/* Makes the sum (*Q) = (*Q) + (*S) */
void addto_32(spinor * const Q, const spinor32 * const R, const int N);


#endif
/* $Id$*/

#ifndef _ASSIGN_ADD_MUL_ADD_MUL_BI_H
#define _ASSIGN_ADD_MUL_ADD_MUL_BI_H

#include "su3.h"

/* (*R) = (*R) + c1*(*S) + c2*(*U) */
void assign_add_mul_add_mul_bi(bispinor * const R, bispinor * const S, bispinor * const U, const complex c1, const complex c2, const int N);


#endif

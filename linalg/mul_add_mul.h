/* $Id$*/

#ifndef _MUL_ADD_MUL_H
#define _MUL_ADD_MUL_H

#include "su3.h"

/* Makes (*R)=c1*(*S)+c2*(*U) , c1 and c2 are complex constants */
void mul_add_mul(spinor * const R,spinor * const S,spinor * const U,const complex c1,const complex c2, const int N);


#endif

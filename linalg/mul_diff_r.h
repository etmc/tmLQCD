/* $Id$*/

#ifndef _MUL_DIFF_R_H
#define _MUL_DIFF_R_H

#include "su3.h"

/* Makes (*R)=c1*(*S) - (*U) , c1 is a real constant */
void mul_diff_r(spinor * const R,spinor * const S,spinor * const U,const double c1, const int N);


#endif

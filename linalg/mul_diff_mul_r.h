/* $Id$*/

#ifndef _MUL_DIFF_MUL_R_H
#define _MUL_DIFF_MUL_R_H

#include "su3.h"

/* Makes (*R)=c1*(*S)-c2*(*U) , c1 and c2 are complex constants */
void mul_diff_mul_r(spinor * const R, spinor * const S, spinor * const U,
		    const double c1, const double c2, const int N);


#endif

/* $Id$*/

#ifndef _MUL_ADD_MUL_R_H
#define _MUL_ADD_MUL_R_H

#include "su3.h"

/* Makes (*R)=c1*(*S)+c2*(*U) , c1 and c2 are real constants */
void mul_add_mul_r(spinor * const R, spinor * const S, spinor * const U,
		   const double c1,const double c2, const int N);


#endif

/* $Id$ */

#ifndef _ASSIGN_ADD_MUL_R_ADD_MUL_H
#define _ASSIGN_ADD_MUL_R_ADD_MUL_H

#include "su3.h"

/* (*R) = (*R) + c1*(*S) + c2*(*U) */
void assign_add_mul_r_add_mul(spinor * const R, spinor * const S, spinor * const U,
			    const double c1,const double c2, const int N);


#endif

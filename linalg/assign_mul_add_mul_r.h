/* $Id$*/

#ifndef _ASSIGN_MUL_ADD_MUL_R_H
#define _ASSIGN_MUL_ADD_MUL_R_H

#include "su3.h"

/* Makes (*R)=c1*(*R)+c2*(*S) , c1 and c2 are real constants */
void assign_mul_add_mul_r(spinor * const R,spinor * const S, 
			  const double c1, const double c2,
			  const int N);

#endif

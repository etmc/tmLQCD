/* $Id$  */

#ifndef _ASSIGN_MUL_BRA_ADD_MUL_KET_ADD_R_H
#define _ASSIGN_MUL_BRA_ADD_MUL_KET_ADD_R_H

#include "su3.h"

/* (*R) =  c2*(*R + c1*(*S)) + (*U) */
void assign_mul_bra_add_mul_ket_add_r(spinor * const R,spinor * const S,spinor * const U,const double c1,const double c2, const int N);

#endif

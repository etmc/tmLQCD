/* $Id$  */

#ifndef _ASSIGN_MUL_BRA_ADD_MUL_KET_ADD_H
#define _ASSIGN_MUL_BRA_ADD_MUL_KET_ADD_H

#include "su3.h"

/* (*R) =  c2*(*R + c1*(*S)) + (*U) */
void assign_mul_bra_add_mul_ket_add(spinor * const R,spinor * const S,spinor * const U,const complex c1,const complex c2);

#endif

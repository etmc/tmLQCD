/* $Id$  */

#ifndef _ASSIGN_MUL_BRA_ADD_MUL_R_H
#define _ASSIGN_MUL_BRA_ADD_MUL_R_H

#include "su3.h"

/*  (*R) = c0*(*R + c*(*S))  */
void assign_mul_bra_add_mul_r(spinor * const R,const double c0, const double c,spinor * const S, const int N);

#endif

/* $Id$  */

#ifndef _MUL_R_H
#define _MUL_R_H

#include "su3.h"

/*   Makes (*R) = c*(*S)   c is a real constant*/
void mul_r(spinor * const R, const double c, spinor * const S, const int N);

#endif

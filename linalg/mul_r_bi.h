/* $Id$  */

#ifndef _MUL_R_BI_H
#define _MUL_R_BI_H

#include "su3.h"

/*   Makes (*R) = c*(*S)   c is a real constant*/
void mul_r_bi(bispinor * const R, const double cup, const double cdn, bispinor * const S, const int N);

#endif

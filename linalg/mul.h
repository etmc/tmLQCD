/* $Id$  */

#ifndef _MUL_H
#define _MUL_H

#include "su3.h"

/*   Makes (*R) = c*(*S) */
void mul(spinor * const R, const complex c, spinor * const S, const int N);

#endif

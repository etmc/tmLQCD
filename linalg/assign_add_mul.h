/* $Id$  */

#ifndef _ASSIGN_ADD_MUL_H
#define _ASSIGN_ADD_MUL_H

#include "su3.h"

/*   (*P) = (*P) + c(*Q)        c is a complex constant   */
void assign_add_mul(spinor * const P, spinor * const Q, const complex c, const int N);

#endif

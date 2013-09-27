#ifndef _ASSIGN_ADD_MUL_32_H
#define _ASSIGN_ADD_MUL_32_H

#include "su3.h"

/*   (*P) = (*P) + c(*Q)        c is a complex constant   */
void assign_add_mul_r_32(spinor32 * const R, spinor32 * const S, const float c, const int N);

#endif

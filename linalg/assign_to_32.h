#ifndef _ASSIGN_TO_32_H
#define _ASSIGN_TO_32_H

#include "su3.h"
void assign_to_32(spinor32 * const R, spinor * const S, const int N);
void assign_to_64(spinor * const R, spinor32 * const S, const int N);
#endif


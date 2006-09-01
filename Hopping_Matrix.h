/* $Id$ */

#ifndef _HOPPING_MATRIX_H
#define _HOPPING_MATRIX_H

#define EO 0
#define OE 1
#define OO 1
#define EE 0

#include "su3.h"

void Hopping_Matrix(const int ieo, spinor * const l, spinor * const k);

#endif

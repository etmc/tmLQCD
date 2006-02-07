/* $Id$  */

#ifndef _COMP_DECOMP_H
#define _COMP_DECOMP_H

#include "su3.h"

/* Build bispinor out of spinors :  (*R) = ((*S), (*T)) */
void compact(bispinor * const R, spinor * const S, spinor * const P);

/* Splits bispinor into spinors :  (*S) = top (*R) ; (*T) = bottom (*R) */
void decompact(spinor * const S, spinor * const P, bispinor * const R);

#endif

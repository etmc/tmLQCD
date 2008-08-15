/* $Id$ */
#ifndef _D_PSI_H
#define _D_PSI_H

#include "deflation/deflation_block.h"

void D_psi(spinor * const P, spinor * const Q);
void Block_D_psi(deflation_block * blk, spinor * const rr, spinor * const s);

#endif

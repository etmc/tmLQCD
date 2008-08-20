/* $Id$ */
#ifndef _D_PSI_H
#define _D_PSI_H

#include "block.h"

void D_psi(spinor * const P, spinor * const Q);
void Block_D_psi(block * blk, spinor * const rr, spinor * const s);

void boundary_D_0(spinor * const r, spinor * const s, su3 *u);
void boundary_D_1(spinor * const r, spinor * const s, su3 *u);
void boundary_D_2(spinor * const r, spinor * const s, su3 *u);
void boundary_D_3(spinor * const r, spinor * const s, su3 *u);
void boundary_D_4(spinor * const r, spinor * const s, su3 *u);
void boundary_D_5(spinor * const r, spinor * const s, su3 *u);
void boundary_D_6(spinor * const r, spinor * const s, su3 *u);
void boundary_D_7(spinor * const r, spinor * const s, su3 *u);


extern void (*boundary_D[])
  (spinor * const r, spinor * const s, su3 * restrict u);

#endif

/* $Id$ */

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "global.h"
#include "su3.h"
#include "start.h"
#include "linalg_eo.h"
#include "tm_operators.h"
#include "boundary.h"
#include "solver.h"
#include "block.h"
#include "D_psi.h"


void dummy_D0(spinor * const P, spinor * const Q) {
  Block_D_psi(&block_list[0], P, Q);
  return;
}

void dummy_D1(spinor * const P, spinor * const Q) {
  Block_D_psi(&block_list[1], P, Q);
  return;
}

void Msap(spinor * const P, spinor * const Q, const int Ncy) {
  int blk, ncy, eo, vol, eolist[2];
  spinor * r, * a, * b;

  r = g_spinor_field[DUM_SOLVER+5];
  a = g_spinor_field[DUM_SOLVER+6];
  b = &g_spinor_field[DUM_SOLVER+6][block_list[0].volume + block_list[0].spinpad];

  if(block_list[0].evenodd == 0) {
    eolist[0] = 0;
    eolist[1] = 1;
  }
  else {
    eolist[0] = 1;
    eolist[1] = 0;
  }
  for(ncy = 0; ncy < Ncy; ncy++) {
    /* even sides first */
    for(eo = 0; eo < 2; eo++) {
      /* compute the global residue */
      D_psi(r, P);
      diff(r, Q, r, VOLUME);
      /* choose the even (odd) block */
      blk = eolist[eo];
      vol = block_list[blk].volume;
      /* get part of r corresponding to block blk into b */
      copy_global_to_upperlower(b, r, blk);
      /* then invert on b block local*/
      if(eolist[eo] == 0) {
	mr(a, b, 100, 1.e-31, 1, vol, 0, &dummy_D0);
      }
      else {
	mr(a, b, 100, 1.e-31, 1, vol, 0, &dummy_D1);
      }
      /* scale a back to full spinor P */
      copy_upperlower_to_global(P, a, blk);
    }
  }
  return;
}

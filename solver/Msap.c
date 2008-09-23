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


void Msap(spinor * const P, spinor * const Q, const int Ncy) {
  int i, j, blk, ncy, eo, vol;
  spinor * r, * tmp, * a, * b;

  r = g_spinor_field[DUM_SOLVER+5];
  for(ncy = 0; ncy < Ncy; ncy++) {
    /* even sides first */
    for(eo = 0; eo < 2; eo++) {
      D_psi(r, P);
      diff(r, Q, r, VOLUME);
      for(blk = 0; blk < 2; blk++) {
	if(block_list[blk].evenodd == eo) {
	  vol = block_list[blk].volume;
	  /* get part of r corresponding to block blk into b */
	  /* to be done */
	  /* then invert on b block local*/
	  mr(a, b, 100, 1.e-31, 1, vol, 0, &Block_D_psi);
	  /* scale a back to full spinor */
	  /* to be done */
	}
      }
    }
  }
  return;
}

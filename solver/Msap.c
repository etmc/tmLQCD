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
#include "gmres.h"
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
  double nrm;

  r = g_spinor_field[DUM_SOLVER+5];
  a = g_spinor_field[DUM_SOLVER+6];
/*   b = g_spinor_field[DUM_SOLVER+7]; */
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
      /* compute the global residue        */
      /* this can be done more efficiently */
      /* here only a naive implementation  */
      D_psi(r, P);
      diff(r, Q, r, VOLUME);
      nrm = square_norm(r, VOLUME, 1);
      if(g_proc_id == 0 && eo == 0) {
	printf("Msap: %d %1.3e\n", ncy, nrm);
      }
      /* choose the even (odd) block */
      blk = eolist[eo];
      vol = block_list[blk].volume;
      /* get part of r corresponding to block blk into b */
      copy_global_to_block(b, r, blk);
      /* then invert on b block local                    */
      /* maybe use a polynomial instead of gmres ?       */
      /* mr does not work, why?                          */
      if(eolist[eo] == 0) {
     	gmres(a, b, 4, 1, 1.e-31, 1, vol, 0, &dummy_D0); 
	/* 	mr(a, b, 4, 1.e-31, 1, vol, 0, &dummy_D0); */
      }
      else {
     	gmres(a, b, 4, 1, 1.e-31, 1, vol, 0, &dummy_D1);
	/* 	mr(a, b, 4, 1.e-31, 1, vol, 0, &dummy_D1); */
      }
      /* add a up to full spinor P */
      add_block_to_global(P, a, blk);
    }
  }
  return;
}

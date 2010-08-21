/***********************************************************************
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008 Carsten Urbach
 *
 * This file is part of tmLQCD.
 *
 * tmLQCD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * tmLQCD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with tmLQCD.  If not, see <http://www.gnu.org/licenses/>.
 ***********************************************************************/
/* $Id$ */

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
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

void dummy_Di(spinor * const P, spinor * const Q, const int i) {
  Block_D_psi(&block_list[i], P, Q);
  return;
}

void dummy_D0(spinor * const P, spinor * const Q) {
  Block_D_psi(&block_list[0], P, Q);
  return;
}

void dummy_D1(spinor * const P, spinor * const Q) {
  Block_D_psi(&block_list[1], P, Q);
  return;
}

void Msap(spinor * const P, spinor * const Q, const int Ncy) {
  int blk, ncy = 0, eo, vol;
  spinor * r, * a, * b;
  double nrm;

  r = g_spinor_field[DUM_SOLVER+5];
  a = g_spinor_field[DUM_SOLVER+6];
  b = g_spinor_field[DUM_SOLVER+7];

  for(ncy = 0; ncy < Ncy; ncy++) {
    /* even sides first */
    /*   for(eo = 0; eo < 2; eo++) { */
    /* compute the global residue        */
    /* this can be done more efficiently */
    /* here only a naive implementation  */
    for(eo = 0; eo < 2; eo++) {
      D_psi(r, P);
      diff(r, Q, r, VOLUME);
      nrm = square_norm(r, VOLUME, 1);
      if(g_proc_id == 0 && g_debug_level > 1) {
	printf("Msap: %d %1.3e\n", ncy, nrm);
      }
      /* choose the even (odd) block */
      
      /*blk = eolist[eo];*/
      
      for (blk = 0; blk < nb_blocks; blk++) {
      	if(block_list[blk].evenodd == eo) {
	  vol = block_list[blk].volume;
	  
	  /* get part of r corresponding to block blk into b */
	  copy_global_to_block(b, r, blk);
	  
	  mrblk(a, b, 16, 1.e-31, 1, vol, &dummy_Di, blk);
	  
	  /* add a up to full spinor P */
	  add_block_to_global(P, a, blk);
	}
      }
    }
  }
  return;
}

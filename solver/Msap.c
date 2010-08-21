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
#include "Hopping_Matrix.h"
#include "D_psi.h"

void dummy_Di(spinor * const P, spinor * const Q, const int i) {
  Block_D_psi(&block_list[i], P, Q);
  return;
}


void Mtm_plus_block_psi(spinor * const l, spinor * const k, const int i) {
  block * blk = &block_list[i];
  int vol = (*blk).volume/2;
  Block_H_psi(blk, g_spinor_field[DUM_MATRIX+1], k, EO);
  mul_one_pm_imu_inv(g_spinor_field[DUM_MATRIX+1], +1., vol);
  Block_H_psi(blk, g_spinor_field[DUM_MATRIX+1], k, OE);
  mul_one_pm_imu_sub_mul(l, k, g_spinor_field[DUM_MATRIX], +1., vol);
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

  if(blk_gauge_eo) {
    init_blocks_gaugefield();
  }

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


void Msap_eo(spinor * const P, spinor * const Q, const int Ncy) {
  int blk, ncy = 0, eo, vol;
  spinor * r, * a, * b;
  double nrm;
  spinor * b_even, * b_odd, * a_even, * a_odd;
  r = g_spinor_field[DUM_SOLVER+5];
  a = g_spinor_field[DUM_SOLVER+6];
  b = g_spinor_field[DUM_SOLVER+7];

  vol = block_list[0].volume/2;
  b_even = b;
  b_odd = b + vol + 1;
  a_even = a;
  a_odd = a + vol + 1;
  if(!blk_gauge_eo) {
    init_blocks_eo_gaugefield();
  }

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
      
      for (blk = 0; blk < nb_blocks; blk++) {
      	if(block_list[blk].evenodd == eo) {
	  /* get part of r corresponding to block blk into b_even and b_odd */
	  copy_global_to_block_eo(b_even, b_odd, r, blk);
	  
	  assign_mul_one_pm_imu_inv(a_even, b_even, +1., vol);
	  Block_H_psi(&block_list[blk], a_odd, a_even, OE);
	  /* a_odd = a_odd + b_odd */
	  assign_mul_add_r(a_odd, +1., b_odd, vol);
	  
	  mrblk(b_odd, a_odd, 3, 1.e-31, 1, vol, &Mtm_plus_block_psi, blk);

	  Block_H_psi(&block_list[blk], b_even, b_odd, EO);
	  mul_one_pm_imu_inv(b_even, +1., vol);
	  /* The sign is plus, since in Hopping_Matrix */
	  /* the minus is missing                      */
	  /* a_even = a_even + b_even */
	  assign_add_mul_r(a_even, b_even, +1., vol);

	  /* add even and odd part up to full spinor P */
	  add_eo_block_to_global(P, a_even, b_odd, blk);
	}
      }
    }
  }
  return;
}

/***********************************************************************
 * $Id$
 * 
 * Copyright (C) 2008 Albert Deuzeman, Siebren Reker, Carsten Urbach
 *               2010 Claude Tadonki
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

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "global.h"
#include "D_psi.h"
#include "linalg_eo.h"
#include "start.h"
#include "xchange_lexicfield.h"
#include "block.h"
#include "su3.h"

#define CALLOC_ERROR_CRASH {printf ("calloc errno : %d\n", errno); errno = 0; return 1;}


int init_blocks_gaugefield();
int init_blocks_geometry();

int **** block_ipt;
int *** bipt__;
int ** bipt_;
int * bipt;
complex * little_A = NULL;

enum{
  NONE = 0,
    T_UP = 1,
    T_DN = 2,
    X_UP = 3,
    X_DN = 4,
    Y_UP = 5,
    Y_DN = 6,
    Z_UP = 7,
    Z_DN = 8
    } Direction;

static void (*boundary_D[8])(spinor * const r, spinor * const s, su3 *u) =
{boundary_D_0, boundary_D_1, boundary_D_2, boundary_D_3, boundary_D_4, boundary_D_5, boundary_D_6, boundary_D_7};



block * block_list = NULL;
static spinor * basis = NULL;
static su3 * u = NULL;
const int spinpad = 1;
static int block_init = 0;

int index_b(int t, int x, int y, int z){
  /* Provides the lexicographic index of (t, x, y, z)
     Claude Tadonki
  */
  return ((t*LX + x)*LY + y)*(LZ) + z;

}

int init_blocks() {
  /*CT: Initialization of global variables for blocks */
  nb_blocks = 1; 
  nblks_t = dfl_nblock_t;
  nblks_x = dfl_nblock_x;
  nblks_y = dfl_nblock_y;
  nblks_z = dfl_nblock_z;
  nblks_dir[0] = nblks_t;
  nblks_dir[1] = nblks_x;
  nblks_dir[2] = nblks_y;
  nblks_dir[3] = nblks_z;
  nb_blocks = nblks_t*nblks_x*nblks_y*nblks_z;
  if(g_debug_level > 0 && g_proc_id == 0) {
    printf("xxxxxxxxxxxxxxxxxxxxxxxxxxx\n");
    printf("Number of deflation blocks = %d\n  n_block_t = %d\n  n_block_x = %d\n  n_block_y = %d\n  n_block_z = %d\n",nb_blocks,nblks_t,nblks_x,nblks_y,nblks_z);
    printf("Number of iteration with the polynomial preconditioner = %d \n",dfl_field_iter);
    printf("Number of iteration in the polynomial preconditioner   = %d \n",dfl_poly_iter);
    printf("xxxxxxxxxxxxxxxxxxxxxxxxxxx\n");
  }

  int i,j;
  free_blocks();
  block_init = 1;
  block_list = calloc(nb_blocks, sizeof(block));
  if((void*)(basis = (spinor*)calloc((nb_blocks + 1) * g_N_s * (VOLUME/nb_blocks + spinpad) + 1, sizeof(spinor))) == NULL) {
    CALLOC_ERROR_CRASH;
  }
  if((void*)(u = (su3*)calloc(1+8*VOLUME, sizeof(su3))) == NULL) {
    CALLOC_ERROR_CRASH;
  }
  for(i = 0; i < nb_blocks; i++) {
    block_list[i].basis = (spinor**)calloc(g_N_s, sizeof(spinor*));
  }

#if ( defined SSE || defined SSE2 || defined SSE3)
  block_list[0].basis[0] = (spinor*)(((unsigned long int)(basis)+ALIGN_BASE)&~ALIGN_BASE);
  block_list[0].u = (su3*)(((unsigned long int)(u)+ALIGN_BASE)&~ALIGN_BASE);
#else
  block_list[0].basis[0] = basis;
  block_list[0].u = u;
#endif
  for(j = 1; j < nb_blocks; j++) {
    block_list[j].basis[0] = block_list[j-1].basis[0] + g_N_s*((VOLUME/nb_blocks) + spinpad) ;
    block_list[j].u = block_list[j-1].u + 8*(VOLUME/nb_blocks);
  }
  for(j = 0; j < nb_blocks; j++) {
    for(i = 1 ; i < g_N_s ; i ++ ) { 
      block_list[j].basis[i] = block_list[j].basis[i-1] + (VOLUME/nb_blocks + spinpad);
    }
  }

  if((void*)(block_ipt = (int****)calloc(T/nblks_t+2,sizeof(int*))) == NULL) return(5);
  if((void*)(bipt__ = (int***)calloc ((T/nblks_t+2)*(LX/nblks_x+2), sizeof(int*))) == NULL) return(4);
  if((void*)(bipt_ = (int**)calloc((T/nblks_t+2)*(LX/nblks_x+2)*(LY/nblks_y+2), sizeof(int*))) == NULL) return(3);
  if((void*)(bipt = (int*)calloc((T/nblks_t+2)*(LX/nblks_x+2)*(LY/nblks_y+2)*(LZ/nblks_z+2), sizeof(int))) == NULL) return(8);
  bipt_[0] = bipt;
  bipt__[0] = bipt_;
  block_ipt[0] = bipt__;
  for(i = 1; i < (T/nblks_t+2)*(LX/nblks_x+2)*(LY/nblks_y+2); i++){
    bipt_[i] = bipt_[i-1]+(LZ/nblks_z+2);
  }
  for(i = 1; i < (T/nblks_t+2)*(LX/nblks_x+2); i++){
    bipt__[i] = bipt__[i-1]+(LY/nblks_y+2);
  }
  for(i = 1; i < (T/nblks_t+2); i++){
    block_ipt[i] = block_ipt[i-1]+(LX/nblks_x+2);
  }

  for (i = 0; i < nb_blocks; ++i) {
    block_list[i].id = i;
    block_list[i].volume = VOLUME/nb_blocks;
    block_list[i].LX = LX/nblks_x;
    block_list[i].LY = LY/nblks_y;
    block_list[i].LZ = LZ/nblks_z;
    block_list[i].T  = T/nblks_t;
    block_list[i].ns = g_N_s;
    block_list[i].spinpad = spinpad;

    for (j = 0 ; j < 6; ++j) {
#ifdef MPI
      block_list[i].mpilocal_neighbour[j] = (g_nb_list[j] == g_cart_id) ? i : -1;
#else
      block_list[i].mpilocal_neighbour[j] = i;
#endif
    }
#ifdef MPI
    block_list[i].mpilocal_neighbour[6] = (i == 0 ? 1 : (g_nb_list[j] == g_cart_id) ? 0 : -1);
    block_list[i].mpilocal_neighbour[7] = (i == 1 ? 0 : (g_nb_list[j] == g_cart_id) ? 1 : -1);
#else
    block_list[i].mpilocal_neighbour[6] = (i == 0 ? 1 : 0);
    block_list[i].mpilocal_neighbour[7] = (i == 0 ? 1 : 0);
#endif
    if(g_debug_level > 4 && g_proc_id == 0) {
      for(j = 0; j < 8; j++) {
	printf("block %d mpilocal_neighbour[%d] = %d\n", i, j, block_list[i].mpilocal_neighbour[j]);
      }
    }

    memcpy(block_list[i].mpilocal_coordinate, g_proc_coords, 4*sizeof(int));
    memcpy(block_list[i].coordinate, g_proc_coords, 3*sizeof(int));
    block_list[i].coordinate[3] = nb_blocks * g_proc_coords[3] + i;
    block_list[i].evenodd = (block_list[i].coordinate[0] + block_list[i].coordinate[1] + 
			     block_list[i].coordinate[2] + block_list[i].coordinate[3]) % 2;

    if ((void*)(block_list[i].idx = calloc(8 * (VOLUME/nb_blocks), sizeof(int))) == NULL)
      CALLOC_ERROR_CRASH;

    for (j = 0; j < g_N_s; j++){ /* write a zero element at the end of every spinor */
      _spinor_null(block_list[i].basis[j][VOLUME/nb_blocks]);
    }

    if ((void*)(block_list[i].little_dirac_operator = calloc(9 * g_N_s * g_N_s, sizeof(complex))) == NULL)
      CALLOC_ERROR_CRASH;
    if ((void*)(block_list[i].little_dirac_operator32 = calloc(9 * g_N_s * g_N_s, sizeof(complex32))) == NULL)
      CALLOC_ERROR_CRASH;

    for (j = 0; j < 9 * g_N_s * g_N_s; ++j){
      _complex_zero(block_list[i].little_dirac_operator[j]);
    }
  }
  
  init_blocks_geometry();
  init_blocks_gaugefield();
  return 0;
}

int free_blocks() {
  int i;
  if(block_init == 1) {
    for(i = 0; i < nb_blocks; ++i) {
      free(block_list[i].basis);
      free(block_list[i].little_dirac_operator);
      free(block_list[i].little_dirac_operator32);
    }
    free(block_ipt);
    free(bipt__);
    free(bipt_);
    free(bipt);
    free(u);
    free(basis);
    free(block_list);
    block_init = 0;
  }
  return 0;
}

int init_blocks_gaugefield() {
  /* Copies the existing gauge field on the processor into the two separate blocks in a form
     that is readable by the block Dirac operator. Specifically, in consecutive memory
     now +t,-t,+x,-x,+y,-y,+z,-z gauge links are stored. This requires double the storage in
     memory. */

  int i, x, y, z, t, ix, ix_new = 0;
  int bx, by, bz, bt;

  for (t = 0; t < T/nblks_t; t++) {
    for (x = 0; x < LX/nblks_x; x++) {
      for (y = 0; y < LY/nblks_y; y++) {
	for (z = 0; z < LZ/nblks_z; z++) {
	  i = 0;
	  for(bt = 0; bt < nblks_t; bt ++) {
	    for(bx = 0; bx < nblks_x; bx ++) {
	      for(by = 0; by < nblks_y; by ++) {
		for(bz = 0; bz < nblks_z; bz ++) {
		  ix = g_ipt[t + bt*(T/nblks_t)][x + bx*(LX/nblks_x)][y + by*(LY/nblks_y)][z + bz*(LZ/nblks_z)];
		  memcpy(block_list[i].u + ix_new,     &g_gauge_field[ ix           ][0], sizeof(su3));
		  memcpy(block_list[i].u + ix_new + 1, &g_gauge_field[ g_idn[ix][0] ][0], sizeof(su3));
		  memcpy(block_list[i].u + ix_new + 2, &g_gauge_field[ ix           ][1], sizeof(su3));
		  memcpy(block_list[i].u + ix_new + 3, &g_gauge_field[ g_idn[ix][1] ][1], sizeof(su3));
		  memcpy(block_list[i].u + ix_new + 4, &g_gauge_field[ ix           ][2], sizeof(su3));
		  memcpy(block_list[i].u + ix_new + 5, &g_gauge_field[ g_idn[ix][2] ][2], sizeof(su3));
		  memcpy(block_list[i].u + ix_new + 6, &g_gauge_field[ ix           ][3], sizeof(su3));
		  memcpy(block_list[i].u + ix_new + 7, &g_gauge_field[ g_idn[ix][3] ][3], sizeof(su3));
		  i++;
		}
	      }
	    }
	  }
	  ix_new += 8;
	}
      }
    }
  }
  return(0);
}

int check_blocks_geometry(block * blk) {
  int i, k=0, x, y, z, t;
  int * itest;
  int * ipt;
  ipt = blk->idx;
  itest = (int*)calloc(blk->volume + blk->spinpad, sizeof(int));
  for(i = 0; i < 8*blk->volume; i++) {
    if(*ipt > blk->volume + blk->spinpad-1 || *ipt < 0) {
      if(g_proc_id == 0) {
        printf("error in block geometry! ipt = %d dir = %d i = %d of %d\n",
               (*ipt), i%8, i/8, blk->volume + blk->spinpad);
      }
    }
    
    itest[*(ipt++)]++;
  }
  
  for(i = 0; i < blk->volume; i++) {
    k += itest[i];
    if(itest[i] < 1 || itest[i] > 8) {
      if(g_proc_id == 0) {
        printf("error in block geometry, itest[%d] = %d\n", i, itest[i]);
      }
    }
  }
  
  if(itest[blk->volume + blk->spinpad-1] != 2*(blk->LX*blk->LY*blk->LZ+blk->T*blk->LX*blk->LY+blk->T*blk->LY*blk->LZ+blk->T*blk->LX*blk->LZ)) {
    if(g_proc_id == 0){
      printf("error in block geometry, boundary points wrong %d != %d\n",
             itest[blk->volume + blk->spinpad-1], 2*(blk->LX*blk->LY*blk->LZ+blk->T*blk->LX*blk->LY+blk->T*blk->LY*blk->LZ+blk->T*blk->LX*blk->LZ));
    }
  }
  k+= itest[blk->volume + blk->spinpad-1];
  if(k != 8*blk->volume) {
    if(g_proc_id == 0){
      printf("error in block geometry, total number of points wrong %d != %d\n",
             k, 8*blk->volume);
    }
  }

  ipt = blk->idx;
  for(t = 0; t < T/nblks_t; t++) {
    for(x = 0; x < LX/nblks_x; x++) {
      for(y = 0; y < LY/nblks_y; y++) {
        for(z = 0; z < LZ/nblks_z; z++) {
          i = block_ipt[t][x][y][z];
          if(t != T/nblks_t-1) {
            if(*ipt != block_ipt[t+1][x][y][z] && g_proc_id == 0)
              printf("Shit +t! %d %d %d %d %d != %d at %d\n",
                     t, x, y, z, *ipt, block_ipt[t+1][x][y][z], i);
          }
          else if(*ipt != VOLUME/nb_blocks)
            printf("Shit +t! %d %d %d %d %d != %d at %d\n",
                   t, x, y, z, *ipt, VOLUME/nb_blocks, i);
          ipt++;
          if(t != 0) {
            if(*ipt != block_ipt[t-1][x][y][z] && g_proc_id == 0)
              printf("Shit -t! %d %d %d %d %d != %d at %d\n",
                     t, x, y, z, *ipt, block_ipt[t+1][x][y][z], i);
          }
          else if(*ipt != VOLUME/nb_blocks)
            printf("Shit -t! %d %d %d %d %d != %d at %d\n",
                   t, x, y, z, *ipt, VOLUME/nb_blocks, i);
          ipt++;
          if(x != LX/nblks_x-1) {
            if(*ipt != block_ipt[t][x+1][y][z] && g_proc_id == 0)
              printf("Shit +x! %d %d %d %d %d != %d at %d\n",
                     t, x, y, z, *ipt, block_ipt[t][x+1][y][z], i);
          }
          else if(*ipt != VOLUME/nb_blocks)
            printf("Shit +x! %d %d %d %d %d != %d at %d\n",
                   t, x, y, z, *ipt, VOLUME/nb_blocks, i);
          ipt++;
          if(x != 0) {
            if(*ipt != block_ipt[t][x-1][y][z] && g_proc_id == 0)
              printf("Shit -x! %d %d %d %d %d != %d at %d\n",
                     t, x, y, z, *ipt, block_ipt[t][x-1][y][z], i);
          }
          else if(*ipt != VOLUME/nb_blocks)
            printf("Shit -x! %d %d %d %d %d != %d at %d\n",
                   t, x, y, z, *ipt, VOLUME/nb_blocks, i);
          ipt++;
          if(y != LY/nblks_y-1) {
            if(*ipt != block_ipt[t][x][y+1][z] && g_proc_id == 0)
              printf("Shit +y! %d %d %d %d %d != %d at %d\n",
                     t, x, y, z, *ipt, block_ipt[t][x][y+1][z], i);
          }
          else if(*ipt != VOLUME/nb_blocks)
            printf("Shit +y! %d %d %d %d %d != %d at %d\n",
                   t, x, y, z, *ipt, VOLUME/nb_blocks, i);
          ipt++;
          if(y != 0) {
            if(*ipt != block_ipt[t][x][y-1][z] && g_proc_id == 0)
              printf("Shit -y! %d %d %d %d %d != %d at %d\n",
                     t, x, y, z, *ipt, block_ipt[t][x][y-1][z], i);
          }
          else if(*ipt != VOLUME/nb_blocks)
            printf("Shit -y! %d %d %d %d %d != %d at %d\n",
                   t, x, y, z, *ipt, VOLUME/nb_blocks, i);
          ipt++;
          if(z != LZ/nblks_z-1) {
            if(*ipt != block_ipt[t][x][y][z+1] && g_proc_id == 0)
              printf("Shit +z! %d %d %d %d %d != %d at %d\n",
                     t, x, y, z, *ipt, block_ipt[t][x][y][z+1], i);
          }
          else if(*ipt != VOLUME/nb_blocks)
            printf("Shit +z! %d %d %d %d %d != %d at %d\n",
                   t, x, y, z, *ipt, VOLUME/nb_blocks, i);
          ipt++;
          if(z != 0) {
            if(*ipt != block_ipt[t][x][y][z-1] && g_proc_id == 0)
              printf("Shit -z! %d %d %d %d %d != %d at %d\n",
                     t, x, y, z, *ipt, block_ipt[t][x][y][z-1], i);
          }
          else if(*ipt != VOLUME/nb_blocks)
            printf("Shit -z! %d %d %d %d %d != %d at %d\n",
                   t, x, y, z, *ipt, VOLUME/nb_blocks, i);
          ipt++;
        }
      }
    }
  }

  free(itest);
  if(g_proc_id == 0) {
    if(g_debug_level > 1) printf("# block geometry checked successfully for block %d !\n", blk->id);
  }
  return(0);
}

int init_blocks_geometry() {
  int i, ix, x, y, z, t;
  int zstride = 1;
  int ystride = LZ/nblks_z;
  int xstride = (LY/nblks_y) * (LZ/nblks_z);
  int tstride = (LX/nblks_x) * (LY/nblks_y) * (LZ/nblks_z);
  int boundidx = VOLUME/nb_blocks;
  for (ix = 0; ix < VOLUME/nb_blocks; ++ix) {
    block_list[0].idx[8 * ix + 0] = (ix           >= VOLUME/nb_blocks - tstride ? boundidx : ix + tstride);/* +t */
    block_list[0].idx[8 * ix + 1] = (ix           <  tstride          ? boundidx : ix - tstride);/* -t */
    block_list[0].idx[8 * ix + 2] = (ix % tstride >= (LZ/nblks_z)*(LY/nblks_y)*((LX/nblks_x)-1)   ? boundidx : ix + xstride);/* +x */
    block_list[0].idx[8 * ix + 3] = (ix % tstride <  (LZ/nblks_z)*(LY/nblks_y)          ? boundidx : ix - xstride);/* -x */
    block_list[0].idx[8 * ix + 4] = (ix % xstride >= (LZ/nblks_z)*((LY/nblks_y) - 1)  ? boundidx : ix + ystride);/* +y */
    block_list[0].idx[8 * ix + 5] = (ix % xstride <  (LZ/nblks_z)             ? boundidx : ix - ystride);/* -y */
    block_list[0].idx[8 * ix + 6] = (ix % ystride == (LZ/nblks_z) - 1         ? boundidx : ix + zstride);/* +z */
    block_list[0].idx[8 * ix + 7] = (ix % ystride == 0                ? boundidx : ix - zstride);/* -z */
  }
  for(i=1;i<nb_blocks;i++) memcpy(block_list[i].idx, block_list[0].idx, 8 * (VOLUME/nb_blocks) * sizeof(int));
  ix = 0;
  for(t = 0; t < T/nblks_t; t++) {
    for(x = 0; x < LX/nblks_x; x++) {
      for(y = 0; y < LY/nblks_y; y++) {
        for(z = 0; z < LZ/nblks_z; z++) {
          block_ipt[t][x][y][z] = ix;
          ix++;
        }
      }
    }
  }
  for(ix = 0; ix < nb_blocks; ix++) {
    zstride = check_blocks_geometry(&block_list[ix]);
  }

  return 0;
}

/* Uses a Modified Gram-Schmidt algorithm to orthonormalize little basis vectors */
void block_orthonormalize(block *parent) {
  int i, j;
  complex coeff;
  double scale;

  for(i = 0; i < g_N_s; ++i) {
    /* rescale the current vector */
    /*     scale = 1. / sqrt(block_two_norm(parent->basis[i], parent->volume)); */
    scale = 1. / sqrt(square_norm(parent->basis[i], parent->volume, 0));
    mul_r(parent->basis[i], scale, parent->basis[i], parent->volume);

    /* rescaling done, now subtract this direction from all vectors that follow */
    for(j = i + 1; j < g_N_s; ++j) {
      /*       coeff = block_scalar_prod(parent->basis[j], parent->basis[i], parent->volume); */
      coeff = scalar_prod(parent->basis[i], parent->basis[j], parent->volume, 0);
      assign_diff_mul(parent->basis[j], parent->basis[i], coeff, parent->volume);
    }
  }

  if(g_debug_level > 4) {
    for(i = 0; i < g_N_s; i++) {
      for(j = 0; j < g_N_s; j++) {
	/*         coeff = block_scalar_prod(parent->basis[j], parent->basis[i], parent->volume); */
        coeff = scalar_prod(parent->basis[i], parent->basis[j], parent->volume, 0);
        if(g_proc_id == 0) printf("basis id = %d <%d, %d> = %1.3e +i %1.3e\n", parent->id, j, i, coeff.re, coeff.im);
      }
    }
  }
  return;
}

void block_orthonormalize_free(block *parent) {
  int i, j;
  complex coeff;
  double scale;

  for(i = 0; i < 12; i++){  /* CHECK THIS !!!!!! 12 */
    /* rescale the current vector */
    constant_spinor_field(parent->basis[i], i, parent->volume);
    scale = 1. / sqrt(square_norm(parent->basis[i], parent->volume, 0));
    mul_r(parent->basis[i], scale, parent->basis[i], parent->volume);
  }

  if(g_debug_level > 4 && g_proc_id == 0) {
    for(i = 0; i < g_N_s; i++) {
      for(j = 0; j < g_N_s; j++) {
	/*         coeff = block_scalar_prod(parent->basis[j], parent->basis[i], parent->volume); */
        coeff = scalar_prod(parent->basis[i], parent->basis[j], parent->volume, 0);
        if(g_proc_id == 0) printf("basis id = %d <%d, %d> = %1.3e +i %1.3e\n", parent->id, j, i, coeff.re, coeff.im);
      }
    }
  }
  return;
}



/* the following 2 functions are reference functions for computing little_d */
/* but much slower than block_compute_little_D_diagonal and                 */
/* block_compute_little_D_offdiagonal                                       */
void block_contract_basis(int const idx, int const vecnum, int const dir, spinor * const psi){
  int l;
  for(l = 0; l < g_N_s; ++l){
    block_list[idx].little_dirac_operator[dir * g_N_s * g_N_s + vecnum * g_N_s + l] =
      /*       block_scalar_prod(psi + idx * VOLUME/2, block_list[idx].basis[l], VOLUME/2); */
      scalar_prod(block_list[idx].basis[l], psi + idx * (VOLUME/nb_blocks+1), VOLUME/nb_blocks, 0);
  }
}

void alt_block_compute_little_D() {
  int i, j, k, l;
  spinor *_rec, *rec, *_app, *app, *zero;
  spinor *psi, **psi_blocks;


  _rec = calloc(VOLUMEPLUSRAND+1, sizeof(spinor));
#if ( defined SSE || defined SSE2 || defined SSE3)
  rec = (spinor*)(((unsigned long int)(_rec)+ALIGN_BASE)&~ALIGN_BASE);
#else
  rec = _rec;
#endif  
  _app = calloc(VOLUMEPLUSRAND+1, sizeof(spinor));
#if ( defined SSE || defined SSE2 || defined SSE3)
  app = (spinor*)(((unsigned long int)(_app)+ALIGN_BASE)&~ALIGN_BASE);
#else
  app = _app;
#endif  
  zero = calloc(VOLUMEPLUSRAND, sizeof(spinor));
  psi = calloc(VOLUME+nb_blocks, sizeof(spinor));
  psi_blocks = (spinor**)calloc(nb_blocks, sizeof(spinor*));
  for(i=0;i<nb_blocks;i++) psi_blocks[i] = psi + i*(VOLUME/nb_blocks + 1);

  for (j = 0; j < VOLUMEPLUSRAND; ++j){
    _spinor_null(zero[j]);
  }

  for (k = 0; k < g_nproc; ++k) {
    for (i = 0; i < g_N_s; ++i) {
      for(l = 0; l < nb_blocks; l++) {
	/* Lower Z block */
	for (j = 0; j < VOLUMEPLUSRAND; ++j){
	  _spinor_null(rec[j]);
	}

	if (g_cart_id == k){
	  reconstruct_global_field_GEN_ID(rec, block_list, i, nb_blocks);
	}

	D_psi(app, rec);

	split_global_field_GEN(psi_blocks, app, nb_blocks);

	if (g_cart_id == k){
	  block_contract_basis(0, i, NONE, psi);
	  block_contract_basis(1, i, Z_DN, psi);
	}
#ifdef MPI
	else if (k == g_nb_t_up){
	  block_contract_basis(0, i, T_UP, psi);
	}
	else if (k == g_nb_t_dn){
	  block_contract_basis(0, i, T_DN, psi);
	}
	else if (k == g_nb_x_up){
	  block_contract_basis(0, i, X_UP, psi);
	}
	else if (k == g_nb_x_dn){
	  block_contract_basis(0, i, X_DN, psi);
	}
	else if (k == g_nb_y_up){
	  block_contract_basis(0, i, Y_UP, psi);
	}
	else if (k == g_nb_y_dn){
	  block_contract_basis(0, i, Y_DN, psi);
	}
	else if (k == g_nb_z_up){
	  block_contract_basis(1, i, Z_UP, psi);
	}
#endif
      }
      /* Upper Z block */
      /*      for (j = 0; j < VOLUMEPLUSRAND; ++j){
	      _spinor_null(rec[j]);
	      }

	      if (g_cart_id == k){
	      reconstruct_global_field(rec, zero, block_list[nb_blocks-1].basis[i]);
	      }

	      D_psi(app, rec);

	      split_global_field(psi_blocks, app);
	      if (g_cart_id == k){
	      block_contract_basis(0, i, Z_UP, psi);
	      block_contract_basis(1, i, NONE, psi);
	      }
	      #ifdef MPI
	      else if (k == g_nb_t_up){
	      block_contract_basis(1, i, T_UP, psi);
	      }
	      else if (k == g_nb_t_dn){
	      block_contract_basis(1, i, T_DN, psi);
	      }
	      else if (k == g_nb_x_up){
	      block_contract_basis(1, i, X_UP, psi);
	      }
	      else if (k == g_nb_x_dn){
	      block_contract_basis(1, i, X_DN, psi);
	      }
	      else if (k == g_nb_y_up){
	      block_contract_basis(1, i, Y_UP, psi);
	      }
	      else if (k == g_nb_y_dn){
	      block_contract_basis(1, i, Y_DN, psi);
	      }
	      else if (k == g_nb_z_dn){
	      block_contract_basis(0, i, Z_DN, psi);
	      }

	      MPI_Barrier(MPI_COMM_WORLD);
	      #endif */
    }
  }

  if(g_debug_level > -1) {
    if (g_N_s <= 5 && g_cart_id == 0){
      printf("\n\n  *** CHECKING LITTLE D ***\n");
      printf("\n  ** node 0, lower block **\n");
      for (i = 0*g_N_s; i < 9 * g_N_s; ++i){
        printf(" [ ");
        for (j = 0; j < g_N_s; ++j){
          printf("%s%1.3e %s %1.3e i", block_list[0].little_dirac_operator[i * g_N_s + j].re >= 0 ? "  " : "- ", block_list[0].little_dirac_operator[i * g_N_s + j].re >= 0 ? block_list[0].little_dirac_operator[i * g_N_s + j].re : -block_list[0].little_dirac_operator[i * g_N_s + j].re, block_list[0].little_dirac_operator[i * g_N_s + j].im >= 0 ? "+" : "-", block_list[0].little_dirac_operator[i * g_N_s + j].im >= 0 ? block_list[0].little_dirac_operator[i * g_N_s + j].im : -block_list[0].little_dirac_operator[i * g_N_s + j].im);
          if (j != g_N_s - 1){
            printf(",\t");
          }
        }
        printf(" ]\n");
        if ((i % g_N_s) == (g_N_s - 1))
          printf("\n");
      }

      printf("\n\n  *** CHECKING LITTLE D ***\n");
      printf("\n  ** node 0, upper block **\n");
      for (i = 0*g_N_s; i < 9 * g_N_s; ++i){
        printf(" [ ");
        for (j = 0; j < g_N_s; ++j){
          printf("%s%1.3e %s %1.3e i", block_list[1].little_dirac_operator[i * g_N_s + j].re >= 0 ? "  " : "- ", block_list[1].little_dirac_operator[i * g_N_s + j].re >= 0 ? block_list[1].little_dirac_operator[i * g_N_s + j].re : -block_list[1].little_dirac_operator[i * g_N_s + j].re, block_list[1].little_dirac_operator[i * g_N_s + j].im >= 0 ? "+" : "-", block_list[1].little_dirac_operator[i * g_N_s + j].im >= 0 ? block_list[1].little_dirac_operator[i * g_N_s + j].im : -block_list[1].little_dirac_operator[i * g_N_s + j].im);
          if (j != g_N_s - 1){
            printf(",\t");
          }
        }
        printf(" ]\n");
        if ((i % g_N_s) == (g_N_s - 1))
          printf("\n");
      }
    }
  }

  free(_rec);
  free(_app);
  free(zero);
  free(psi);
}


/* checked CU */
void compute_little_D_diagonal() {
  int i,j, blk;
  spinor * tmp, * _tmp;
  complex * M;
  _tmp = calloc( block_list[0].volume + block_list[0].spinpad + 1, sizeof(spinor));
#if ( defined SSE || defined SSE2 || defined SSE3)
  tmp = (spinor*)(((unsigned long int)(_tmp)+ALIGN_BASE)&~ALIGN_BASE);
#else
  tmp = _tmp;
#endif  

  for(blk = 0; blk < nb_blocks; blk++) {
    M = block_list[blk].little_dirac_operator;
    for(i = 0; i < g_N_s; i++) {
      Block_D_psi(&block_list[blk], tmp, block_list[blk].basis[i]);
      for(j = 0; j < g_N_s; j++) {
	M[i * g_N_s + j]  = scalar_prod(block_list[blk].basis[j], tmp, block_list[blk].volume, 0);
      }
    }
  }
  free(_tmp);
  return;
}



/* what happens if this routine is called in a one dimensional parallelisation? */
/* or even serially ?                                                           */
/* checked CU */
void compute_little_D(){
  /* this still needs adaption to the new block structure!! */

  spinor *scratch, * temp, *_scratch;
  spinor *r, *s;
  su3 * u;
  int x, y, z, t, ix, iy, i, j, k, pm, mu, blk;
  complex c, *M;

  int bx, by, bz, bt, block_id;
  int dT, dX, dY, dZ;
  dT = T/nblks_t; dX = LX/nblks_x; dY = LY/nblks_y; dZ = LZ/nblks_z;

  /* for a full spinor field we need VOLUMEPLUSRAND                 */
  /* because we use the same geometry as for the                    */
  /* gauge field                                                    */
  /* It is VOLUME + 2*LZ*(LY*LX + T*LY + T*LX) + 4*LZ*(LY + T + LX) */
  _scratch = calloc(2*VOLUMEPLUSRAND+1, sizeof(spinor));
#if ( defined SSE || defined SSE2 || defined SSE3)
  scratch = (spinor*)(((unsigned long int)(_scratch)+ALIGN_BASE)&~ALIGN_BASE);
#else
  scratch = _scratch;
#endif
  temp = scratch + VOLUMEPLUSRAND;

  for(blk = 0; blk < nb_blocks; blk++) {
    M = block_list[blk].little_dirac_operator;
    for(i = 0; i < g_N_s; i++) {
      Block_D_psi(&block_list[blk], scratch, block_list[blk].basis[i]);
      for(j = 0; j < g_N_s; j++) {
	M[i * g_N_s + j]  = scalar_prod(block_list[blk].basis[j], scratch, block_list[blk].volume, 0);
      }
    }
  }


  for (i = 0; i < g_N_s; i++){
    reconstruct_global_field_GEN_ID(scratch, block_list,i, nb_blocks);
#ifdef MPI
    xchange_lexicfield(scratch);
#endif
    /* this sets all elements of scratch to zero, but the boundaries */
    /* and this is what we need here!!                               */
    zero_spinor_field(scratch, VOLUME);
    /* +- t */
    mu = 0;
    for(pm = 0; pm < 2; pm++) {
      if(pm == 0) t = T/nblks_t-1;
      else t = 0;

      r = temp;
      for(x = 0; x < LX; x++) {
	for(y = 0; y < LY; y++) {
	  for(z = 0; z < LZ; z++) {
	    ix = g_ipt[t][x][y][z];
	    if(pm == 0) {
	      s = &scratch[ g_iup[ ix ][mu] ];
	      u = &g_gauge_field[ ix ][mu];
	    }
	    else {
	      s = &scratch[ g_idn[ ix ][mu] ];
	      u = &g_gauge_field[ g_idn[ix][mu] ][mu];
	    }
	    boundary_D[pm](r, s, u);
	    r++;
	  }
	}
      }

      /* now all the scalar products */
      for(j = 0; j < g_N_s; j++) {
	iy = i * g_N_s + j  + (pm + 1) * g_N_s * g_N_s;
	block_id = 0;
	for(bt = 0; bt < nblks_t; bt++) {
	  for(bx = 0; bx < nblks_x; bx++) {
            for(by = 0; by < nblks_y; by++) {
	      for(bz = 0; bz < nblks_z; bz++) {
		_complex_zero(block_list[block_id].little_dirac_operator[ iy ]);
		/* We need to contract g_N_s times with the same set of fields, right? */
		for(x = 0; x < dX; x++) {
		  for(y = 0; y < dY; y++) {
		    for(z = 0; z < dZ; z++) {
		      r = temp + index_b(dT*bt + t, dX*bx + x, dY*by + y, dZ*bz + z);/* CT: TO BE INLINED */
		      ix = block_ipt[t][x][y][z];
		      s = &block_list[block_id].basis[j][ ix ];
		      c = scalar_prod(s, r, 1, 0);/* CT: TO BE INLINED */
		      block_list[block_id].little_dirac_operator[ iy ].re += c.re;
		      block_list[block_id].little_dirac_operator[ iy ].im += c.im;
		    }
		  }
		}
                block_id++;
	      }
	    }
	  }
	}
      }
    }

    /* +- x */
    mu = 1;
    for(pm = 2; pm < 4; pm++) {
      if(pm == 2) x = LX/nblks_x-1;
      else x = 0;

      r = temp;
      for(t = 0; t < T; t++) {
        for(y = 0; y < LY; y++) {
          for(z = 0; z < LZ; z++) {
            ix = g_ipt[t][x][y][z];
            if(pm == 2) {
              s = &scratch[ g_iup[ ix ][mu] ];
              u = &g_gauge_field[ ix ][mu];
            }
            else {
              s = &scratch[ g_idn[ ix ][mu] ];
              u = &g_gauge_field[ g_idn[ix][mu] ][mu];
            }
            boundary_D[pm](r, s, u);
            r++;
          }
        }
      }

      /* now all the scalar products */
      for(j = 0; j < g_N_s; j++) {
	iy = i * g_N_s + j  + (pm + 1) * g_N_s * g_N_s;
	block_id = 0;
	for(bt = 0; bt < nblks_t; bt++) {
	  for(bx = 0; bx < nblks_x; bx++) {
            for(by = 0; by < nblks_y; by++) {
	      for(bz = 0; bz < nblks_z; bz++) {
		_complex_zero(block_list[block_id].little_dirac_operator[ iy ]);
		/* We need to contract g_N_s times with the same set of fields, right? */
		for(t = 0; t < dT; t++) {
		  for(y = 0; y < dY; y++) {
		    for(z = 0; z < dZ; z++) {
		      r = temp + index_b(dT*bt + t, dX*bx + x, dY*by + y, dZ*bz + z);/* CT: TO BE INLINED */
		      ix = block_ipt[t][x][y][z];
		      s = &block_list[block_id].basis[j][ ix ];
		      c = scalar_prod(s, r, 1, 0);/* CT: TO BE INLINED */
		      block_list[block_id].little_dirac_operator[ iy ].re += c.re;
		      block_list[block_id].little_dirac_operator[ iy ].im += c.im;
		    }
		  }
		}
                block_id++;
	      }
	    }
	  }
	}
      }
    }

    /* +- y */
    mu = 2;
    for(pm = 4; pm < 6; pm++) {
      if(pm == 4) y = LY/nblks_y-1;
      else y = 0;

      r = temp;
      for(t = 0; t < T; t++) {
        for(x = 0; x < LX; x++) {
          for(z = 0; z < LZ; z++) {
            ix = g_ipt[t][x][y][z];
            if(pm == 4) {
              s = &scratch[ g_iup[ ix ][mu] ];
              u = &g_gauge_field[ ix ][mu];
            }
            else {
              s = &scratch[ g_idn[ ix ][mu] ];
              u = &g_gauge_field[ g_idn[ix][mu] ][mu];
            }
            /* boundary_D[pm](r, s, u); */
            r++;
          }
        }
      }

      /* now all the scalar products */
      for(j = 0; j < g_N_s; j++) {
	iy = i * g_N_s + j  + (pm + 1) * g_N_s * g_N_s;
	block_id = 0;
	for(bt = 0; bt < nblks_t; bt++) {
	  for(bx = 0; bx < nblks_x; bx++) {
            for(by = 0; by < nblks_y; by++) {
	      for(bz = 0; bz < nblks_z; bz++) {
		_complex_zero(block_list[block_id].little_dirac_operator[ iy ]);
		/* We need to contract g_N_s times with the same set of fields, right? */
		for(t = 0; t < dT; t++) {
		  for(x = 0; x < dX; x++) {
		    for(z = 0; z < dZ; z++) {
		      r = temp + index_b(dT*bt + t, dX*bx + x, dY*by + y, dZ*bz + z);/*  CT: TO BE INLINED */
		      ix = block_ipt[t][x][y][z];
		      s = &block_list[block_id].basis[j][ ix ];
		      c = scalar_prod(s, r, 1, 0);/* CT: TO BE INLINED */
		      block_list[block_id].little_dirac_operator[ iy ].re += c.re;
		      block_list[block_id].little_dirac_operator[ iy ].im += c.im;
		    }
		  }
		}
                block_id++;
	      }
	    }
	  }
	}
      }
    }

    /* z is different */
    /* +-z */
    mu = 3;
    for(pm = 6; pm < 8; pm++) {
      if(pm == 6) z = LZ/nblks_z - 1;
      else z = 0;

      r = temp;
      for(t = 0; t < T; ++t) {
        for(x = 0; x < LX; ++x) {
          for(y = 0; y < LY; ++y) {
            ix = g_ipt[t][x][y][z];
            if(pm == 6) {
              s = &scratch[ g_iup[ ix ][mu] ];
              u = &g_gauge_field[ ix ][mu];
            }
            else {
              s = &scratch[ g_idn[ ix ][mu] ];
              u = &g_gauge_field[ g_idn[ix][mu] ][mu];
            }
            boundary_D[pm](r, s, u);
            r++;
          }
        }
      }
      /*if(pm == 6) z = LZ/nb_blocks-1; */
      /* now for all the MPI scalar products (outer edges) */
      /* this is block 0 -z and block 1 +z CT: should be generalize to blocks (...,0,...) and (...,nb_blocks-1,...)*/
      for(j = 0; j < g_N_s; j++) {
	iy = i * g_N_s + j + (pm + 1) * g_N_s * g_N_s;

	_complex_zero(block_list[(pm+1) % nb_blocks].little_dirac_operator[ iy ]);
	r = temp;
	for(t = 0; t < T; t++) {
	  for(x = 0; x < LX; x++) {
	    for(y = 0; y < LY; y++){
	      ix = block_ipt[t][x][y][z];
	      s = &block_list[(pm+1) % nb_blocks].basis[j][ ix ];
	      c = scalar_prod(s, r, 1, 0);
	      block_list[(pm+1) % nb_blocks].little_dirac_operator[ iy ].re += c.re;
	      block_list[(pm+1) % nb_blocks].little_dirac_operator[ iy ].im += c.im;
	      r++;
	    }
	  }
	}
      }
    }

    pm = 6;
    /* we are on block 0 now and compute direction +z */
    r = temp;
    z = 0;
    for(t = 0; t < T; t++) {
      for(x = 0; x < LX; x++) {
	for(y = 0; y < LY; y++) {
	  ix = g_ipt[t][x][y][LZ/nblks_z-1];
	  iy = block_ipt[t][x][y][z]; /* lowest edge of upper block needed */
	  s = &block_list[ (pm + 1) % nb_blocks ].basis[ i ][ iy ];
	  u = &g_gauge_field[ ix ][mu];
	  boundary_D[pm](r, s, u);
	  r++;
	}
      }
    }
    
    for(j = 0; j < g_N_s; ++j) {
      iy = i * g_N_s + j + (pm + 1) * g_N_s * g_N_s;
      _complex_zero(block_list[(pm) % nb_blocks].little_dirac_operator[ iy ]);
      r = temp;
      for(t = 0; t < T; ++t) {
	for(x = 0; x < LX; ++x) {
	  for(y = 0; y < LY; ++y){
	    ix = block_ipt[t][x][y][LZ/nblks_z-1];
	    s = &block_list[(pm)%nb_blocks].basis[j][ ix ];
	    _add_complex(block_list[(pm) % nb_blocks].little_dirac_operator[ iy ], scalar_prod(s, r, 1, 0));
	    r++;
	  }
	}
      }
    }
    
    /* now we are on block 1 and compute direction -z */
    pm = 7;
    z =  LZ/nblks_z -1;
    r = temp;
    for(t = 0; t < T; ++t) {
      for(x = 0; x < LX; ++x) {
	for(y = 0; y < LY; ++y) {
	  ix = g_ipt[t][x][y][LZ/nblks_z];
	  iy = block_ipt[t][x][y][z];  /* highest edge of lower block needed */
	  s = &block_list[(pm+1) % nb_blocks].basis[ i ][ iy ];
	  u = &g_gauge_field[ g_idn[ ix ][ mu ] ][mu];
	  boundary_D[pm](r, s, u);
	  r++;
	}
      }
    }
      
    /* Now contract with the highest edge of the low block and store */
    for(j = 0; j < g_N_s; ++j) {
      iy = i * g_N_s + j + (pm + 1) * g_N_s * g_N_s;
      _complex_zero(block_list[(pm) % nb_blocks].little_dirac_operator[ iy ]);
      r = temp;
      for(t = 0; t < T; ++t) {
	for(x = 0; x < LX; ++x) {
	  for(y = 0; y < LY; ++y) {
	    ix = block_ipt[t][x][y][0]; 
	    s = &block_list[pm % nb_blocks].basis[j][ ix ];
	    _add_complex(block_list[(pm) % nb_blocks].little_dirac_operator[ iy ], scalar_prod(s, r, 1, 0));
	    r++;
	  }
	}
      }
    }
  }

  if(g_debug_level > 3) {
    if (g_N_s <= 5 && !g_cart_id){
      printf("\n\n  *** CHECKING LITTLE D ***\n");
      printf("\n  ** node 0, lower block **\n");
      for (i = 0*g_N_s; i < 9 * g_N_s; ++i){
        printf(" [ ");
        for (j = 0; j < g_N_s; ++j) {
          printf("%s%1.3e %s %1.3e i", block_list[0].little_dirac_operator[i * g_N_s + j].re >= 0 ? "  " : "- ", block_list[0].little_dirac_operator[i * g_N_s + j].re >= 0 ? block_list[0].little_dirac_operator[i * g_N_s + j].re : -block_list[0].little_dirac_operator[i * g_N_s + j].re, block_list[0].little_dirac_operator[i * g_N_s + j].im >= 0 ? "+" : "-", block_list[0].little_dirac_operator[i * g_N_s + j].im >= 0 ? block_list[0].little_dirac_operator[i * g_N_s + j].im : -block_list[0].little_dirac_operator[i * g_N_s + j].im);
          if (j != g_N_s - 1){
            printf(",\t");
          }
        }
        printf(" ]\n");
        if ((i % g_N_s) == (g_N_s - 1)) {
          printf("\n");
	}
      }
      
      printf("\n\n  *** CHECKING LITTLE D ***\n");
      printf("\n  ** node 0, upper block **\n");
      for (i = 0*g_N_s; i < 9 * g_N_s; ++i) {
        printf(" [ ");
        for (j = 0; j < g_N_s; ++j) {
          printf("%s%1.3e %s %1.3e i", block_list[1].little_dirac_operator[i * g_N_s + j].re >= 0 ? "  " : "- ", block_list[1].little_dirac_operator[i * g_N_s + j].re >= 0 ? block_list[1].little_dirac_operator[i * g_N_s + j].re : -block_list[1].little_dirac_operator[i * g_N_s + j].re, block_list[1].little_dirac_operator[i * g_N_s + j].im >= 0 ? "+" : "-", block_list[1].little_dirac_operator[i * g_N_s + j].im >= 0 ? block_list[1].little_dirac_operator[i * g_N_s + j].im : -block_list[1].little_dirac_operator[i * g_N_s + j].im);
          if (j != g_N_s - 1) {
            printf(",\t");
          }
        }
        printf(" ]\n");
        if ((i % g_N_s) == (g_N_s - 1)) {
          printf("\n");
	}
      }
    }
  }
  
  free(_scratch);
  return;
}


int split_global_field(spinor * const block_low, spinor * const block_high, spinor * const field) {
  int ctr_t;

  for (ctr_t = 0; ctr_t < (VOLUME / LZ); ctr_t++) {
    memcpy(block_low + ctr_t * LZ / 2, field + (2 * ctr_t) * LZ / 2, LZ / 2 * sizeof(spinor));
    memcpy(block_high + ctr_t * LZ / 2, field + (2 * ctr_t + 1) * LZ / 2, LZ / 2 * sizeof(spinor));
  }
  if(g_proc_id == 0 && g_debug_level > 8) {
    printf("lower basis norm = %1.3e\n", square_norm(block_low,  VOLUME / 2, 0));
    printf("upper basis norm = %1.3e\n", square_norm(block_high, VOLUME / 2, 0));
  }

  return 0;
}

int split_global_field_GEN(spinor ** const psi, spinor * const field, int nb_blocks) {
  int j,ctr_t=0;
  int i, x, y, z, t, ix, ix_new = 0;
  int bx, by, bz, bt, block_id;
  int dT, dX, dY, dZ;
  dT = T/nblks_t; dX = LX/nblks_x; dY = LY/nblks_y; dZ = LZ/nblks_z;
  for (t = 0; t < dT; t++) {
    for (x = 0; x < dX; x++) {
      for (y = 0; y < dY; y++) {
	for (z = 0; z < dZ; z++) {
	  block_id = 0;
	  for(bt = 0; bt < nblks_t; bt++) {
	    for(bx = 0; bx < nblks_x; bx++) {
	      for(by = 0; by < nblks_y; by++) {
		for(bz = 0; bz < nblks_z; bz++) {
		  _spinor_assign(*(psi[block_id] + ctr_t), *(field + index_b(dT*bt + t, dX*bx + x, dY*by + y, dZ*bz + z)));
		  block_id++;
		}
	      }
	    }
	  }
	  ctr_t++;
	}
      }
    }
  }

  if(g_proc_id == 0 && g_debug_level > 8) {
    for(j = 0;j < nb_blocks;j++)
      printf("Basis norm %2d = %1.3e\n", j, square_norm(psi[j],  VOLUME / nb_blocks, 0));
  }
  return 0;
}

int split_global_field_GEN_ID(block * const block_list, int id, spinor * const field, int nb_blocks){
  int j,ctr_t=0;
  int i, x, y, z, t, ix, ix_new = 0;
  int bx, by, bz, bt, block_id;
  int dT, dX, dY, dZ;
  dT = T/nblks_t; dX = LX/nblks_x; dY = LY/nblks_y; dZ = LZ/nblks_z;
  for (t = 0; t < dT; t++) {
    for (x = 0; x < dX; x++) {
      for (y = 0; y < dY; y++) {
	for (z = 0; z < dZ; z++) {
	  block_id = 0;
	  for(bt = 0; bt < nblks_t; bt++) {
	    for(bx = 0; bx < nblks_x; bx++) {
	      for(by = 0; by < nblks_y; by++) {
		for(bz = 0; bz < nblks_z; bz++) {
		  _spinor_assign(*(block_list[block_id].basis[id] + ctr_t), *(field + index_b(dT*bt + t, dX*bx + x, dY*by + y, dZ*bz + z)));
		  block_id++;
		}
	      }
	    }
	  }
	  ctr_t++;
	}
      }
    }
  }

  if(g_proc_id == 0 && g_debug_level > 8) {
    for(j = 0; j < nb_blocks; j++) {
      printf("Basis norm %2d = %1.3e\n", j, square_norm(block_list[j].basis[id],  VOLUME / nb_blocks, 0));
    }
  }
  return 0;
}

/* copies the part of globalfields corresponding to block blk */
/* to the block field blockfield                              */
void copy_global_to_block(spinor * const blockfield, spinor * const globalfield, const int blk) {
  int i;
  /*CT: This procedure should changed if multi-dimensional split is considered */
  for (i = 0; i < VOLUME/LZ; i++) {
    memcpy(blockfield + i*(LZ/nb_blocks), globalfield + (nb_blocks * i + blk)*(LZ/nb_blocks), (LZ/nb_blocks)*sizeof(spinor));
  }
  return;
}


/* reconstructs the parts of globalfield corresponding to block blk */
/* from block field blockfield                                      */
void copy_block_to_global(spinor * const globalfield, spinor * const blockfield, const int blk) {
  int i, vol = block_list[blk].volume;
  /*CT: This procedure should changed if multi-dimensional split is considered */
  for (i = 0; i < (vol / (LZ / nb_blocks)); ++i) {
    memcpy(globalfield + (nb_blocks*i + blk)*(LZ/nb_blocks), blockfield + i*(LZ/nb_blocks), (LZ/nb_blocks)*sizeof(spinor));
  }
  return;
}


/* Reconstructs a global field from the little basis of two blocks */
void reconstruct_global_field(spinor * const rec_field, spinor * const block_low, spinor * const block_high) {
  int ctr_t;
  for (ctr_t = 0; ctr_t < ((VOLUME/2) / (LZ / 2)); ++ctr_t) {
    memcpy(rec_field + (2 * ctr_t) * LZ / 2, block_low + ctr_t * LZ / 2, LZ / 2 * sizeof(spinor));
    memcpy(rec_field + (2 * ctr_t + 1) * LZ / 2, block_high + ctr_t * LZ / 2, LZ / 2 * sizeof(spinor));
  }
  return;
}

/* Reconstructs a global field from the little basis of nb_blocks blocks */
void reconstruct_global_field_GEN(spinor * const rec_field, spinor ** const psi, int nb_blocks) {
  int j,ctr_t=0;
  int i, x, y, z, t, ix, ix_new = 0;
  int bx, by, bz, bt, block_id;
  int dT, dX, dY, dZ;
  dT = T/nblks_t; dX = LX/nblks_x; dY = LY/nblks_y; dZ = LZ/nblks_z;
  for (t = 0; t < dT; t++) {
    for (x = 0; x < dX; x++) {
      for (y = 0; y < dY; y++) {
	for (z = 0; z < dZ; z++) {
	  block_id = 0;
	  for(bt = 0; bt < nblks_t; bt++) {
	    for(bx = 0; bx < nblks_x; bx++) {
	      for(by = 0; by < nblks_y; by++) {
		for(bz = 0; bz < nblks_z; bz++) {
		  _spinor_assign(*(rec_field + index_b(dT*bt + t, dX*bx + x, dY*by + y, dZ*bz + z)), *(psi[block_id] + ctr_t));
		  block_id++;
		}
	      }
	    }
	  }
	  ctr_t++;
	}
      }
    }
  }
  return;
}

/* Reconstructs a global field from the little basis of nb_blocks blocks taken from block_list[*].basis[id] */
void reconstruct_global_field_GEN_ID(spinor * const rec_field, block * const block_list, int id, int nb_blocks) {
  int j,ctr_t=0;
  int i, x, y, z, t, ix, ix_new = 0;
  int bx, by, bz, bt, block_id;
  int dT, dX, dY, dZ;
  dT = T/nblks_t; dX = LX/nblks_x; dY = LY/nblks_y; dZ = LZ/nblks_z;
  for (t = 0; t < dT; t++) {
    for (x = 0; x < dX; x++) {
      for (y = 0; y < dY; y++) {
	for (z = 0; z < dZ; z++) {
	  block_id = 0;
	  for(bt = 0; bt < nblks_t; bt++) {
	    for(bx = 0; bx < nblks_x; bx++) {
	      for(by = 0; by < nblks_y; by++) {
		for(bz = 0; bz < nblks_z; bz++) {
		  _spinor_assign(*(rec_field + index_b(dT*bt + t, dX*bx + x, dY*by + y, dZ*bz + z)), *(block_list[block_id].basis[id] + ctr_t));
		  block_id++;
		}
	      }
	    }
	  }
	  ctr_t++;
	}
      }
    }
  }
  return;
}

void add_block_to_global(spinor * const globalfield, spinor * const blockfield, const int blk) {
  int i, vol = block_list[blk].volume;
  spinor * r, * s;
  /*CT: This procedure should changed if multi-dimensional split is considered */
  for (i = 0; i < (vol / (LZ / nb_blocks)); i++) {
    r = globalfield + (nb_blocks*i + blk)*(LZ/nb_blocks);
    s = blockfield + i*(LZ/nb_blocks);
    add(r, r, s, LZ/nb_blocks);
  }  
  return;
}

/***********************************************************************
 * 
 * Copyright (C) 2008 Albert Deuzeman, Siebren Reker, Carsten Urbach
 *               2010 Claude Tadonki, Carsten Urbach
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
#include <complex.h>
#include "global.h"
#include "aligned_malloc.h"
#include "boundary.h"
#include "gettime.h"
#include "read_input.h"
#include "operator/D_psi.h"
#include "linalg_eo.h"
#include "start.h"
#include "gamma.h"
#include "xchange/xchange.h"
#include "block.h"
#include "solver/lu_solve.h"
#include "su3.h"

#define CALLOC_ERROR_CRASH {printf ("calloc errno : %d\n", errno); errno = 0; return 1;}


int init_blocks_geometry();

int **** block_ipt;
int *** bipt__;
int ** bipt_;
int * bipt;
_Complex double * little_A = NULL;
_Complex float * little_A32 = NULL;
_Complex double * little_A_eo = NULL;
_Complex float * little_A_eo_32 = NULL;
int * block_idx;
int * block_evenidx;
int * block_oddidx;
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
static su3_32 * u_32 = NULL;
const int spinpad = 1;
static int block_init = 0;

int dT, dX, dY, dZ; /* Block dimension */
  

inline int index_a(int t, int x, int y, int z){
  /* Provides the absolute lexicographic index of (t, x, y, z)
     Useful to walk over the blocks, maybe could be just g_ipt[t][x][y][z]
     Claude Tadonki (claude.tadonki@u-psud.fr)
  */
  return ((t*LX + x)*LY + y)*(LZ) + z;
}
inline int index_b(int t, int x, int y, int z){
  /* Provides the block lexicographic index of (t, x, y, z)
     Useful to walk inside a block
     Claude Tadonki (claude.tadonki@u-psud.fr)
  */
  return ((t*dX + x)*dY + y)*(dZ) + z;
}
inline int block_index(int t, int x, int y, int z){
  /* Provides the lexicographic index of the block (t, x, y, z)
     Useful to walk over the blocks
     Claude Tadonki (claude.tadonki@u-psud.fr)
  */
  return ((t*nblks_x + x)*nblks_y + y)*(nblks_z) + z;
}

int init_blocks(const int nt, const int nx, const int ny, const int nz) {
  int i,j;
  double atime, etime;
  /* Initialization of block-global variables for blocks */
  atime = gettime();
  nb_blocks = 1; 
  nblks_t = nt;
  nblks_x = nx;
  nblks_y = ny;
  nblks_z = nz;
  blk_gauge_eo = -1;
  nblks_dir[0] = nblks_t;
  nblks_dir[1] = nblks_x;
  nblks_dir[2] = nblks_y;
  nblks_dir[3] = nblks_z;
  nb_blocks = nblks_t*nblks_x*nblks_y*nblks_z;
  //if(nblks_t%2 == 1 || nblks_x%2 == 1 || nblks_y%2 == 1 || nblks_z%2 == 1 ) {
  //    fprintf(stderr, "no of blocks in all directions must be even! Aborting...\n");
  if(nblks_z%2 == 1) {
    fprintf(stderr, "no of MPI local blocks in z direction must be even! Aborting...\n");
    exit(0);
  }
  dT = T/nblks_t; 
  dX = LX/nblks_x; 
  dY = LY/nblks_y; 
  dZ = LZ/nblks_z;
  if(g_proc_id == 0 && g_debug_level > 0) {
    printf("# Number of deflation blocks per MPI process = %d\n  n_block_t = %d\n  n_block_x = %d\n  n_block_y = %d\n  n_block_z = %d\n",
           nb_blocks, nblks_t, nblks_x, nblks_y, nblks_z);
    printf("# Block size: %d x %d x %d x %d\n", dT, dX, dY, dZ);
  }
  
  free_blocks();
  block_init = 1;
  block_list = calloc(nb_blocks, sizeof(block));
  if((void*)(basis = (spinor*)aligned_malloc_zero((nb_blocks + 1) * g_N_s * (VOLUME/nb_blocks + spinpad)*sizeof(spinor))) == NULL) {
    CALLOC_ERROR_CRASH;
  }
  if((void*)(u = (su3*)aligned_malloc_zero(8 * VOLUME * sizeof(su3))) == 0) {
    CALLOC_ERROR_CRASH;
  }
  if((void*)(u_32 = (su3_32*)aligned_malloc_zero(8 * VOLUME * sizeof(su3_32))) == 0) {
    CALLOC_ERROR_CRASH;
  }
  for(i = 0; i < nb_blocks; i++) {
    block_list[i].basis = (spinor**)calloc(g_N_s, sizeof(spinor*));
  }
  
  block_list[0].basis[0] = basis;
  block_list[0].u = u;
  block_list[0].u_32 = u_32;

  for(j = 1; j < nb_blocks; j++) { 
    block_list[j].basis[0] = block_list[j-1].basis[0] + g_N_s*((VOLUME/nb_blocks) + spinpad) ;
    block_list[j].u = block_list[j-1].u + 8*(VOLUME/nb_blocks);
    block_list[j].u_32 = block_list[j-1].u_32 + 8*(VOLUME/nb_blocks);
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
  if((void*)(index_block_eo = (int*)calloc(nblks_t*nblks_x*nblks_y*nblks_z, sizeof(int))) == NULL) return(8);
  bipt_[0] = bipt;
  bipt__[0] = bipt_;
  block_ipt[0] = bipt__;
  for(i = 1; i < (T/nblks_t+2)*(LX/nblks_x+2)*(LY/nblks_y+2); i++) {
    bipt_[i] = bipt_[i-1]+(LZ/nblks_z+2);
  }
  for(i = 1; i < (T/nblks_t+2)*(LX/nblks_x+2); i++) {
    bipt__[i] = bipt__[i-1]+(LY/nblks_y+2);
  }
  for(i = 1; i < (T/nblks_t+2); i++) {
    block_ipt[i] = block_ipt[i-1]+(LX/nblks_x+2);
  }
  
  for (i = 0; i < nb_blocks; ++i) {
    block_list[i].id = i;
    block_list[i].volume = VOLUME/nb_blocks;
    block_list[i].BLX = LX/nblks_x;
    block_list[i].BLY = LY/nblks_y;
    block_list[i].BLZ = LZ/nblks_z;
    block_list[i].BT  = T/nblks_t;
    block_list[i].ns = g_N_s;
    block_list[i].spinpad = spinpad;

    // block coordinate on the mpilocal processor
    block_list[i].mpilocal_coordinate[0] = (i / (nblks_x * nblks_y * nblks_z));
    block_list[i].mpilocal_coordinate[1] = (i / (nblks_y * nblks_z)) % nblks_x;
    block_list[i].mpilocal_coordinate[2] = (i / (nblks_z)) % nblks_y;
    block_list[i].mpilocal_coordinate[3] = i % nblks_z;

    // global block coordinate
    for(j = 0; j < 4; j++) {
      block_list[i].coordinate[j] = nblks_dir[j] * g_proc_coords[j] + block_list[i].mpilocal_coordinate[j];
    }
    // even/odd id of block coordinate
    // using the global coordinates here should allow for
    // odd (MPI) local no of blocks 
    // -> of course the global number of blocks in each direction must be even.
    // and at least the local no of blocks in z-direction must be even.
    block_list[i].evenodd = (block_list[i].coordinate[0] + block_list[i].coordinate[1] + 
                             block_list[i].coordinate[2] + block_list[i].coordinate[3]) % 2;

    block_list[i].evenodd_id = block_list[i].id / 2;
    if(block_list[i].evenodd) block_list[i].evenodd_id += nb_blocks/2;

    /* block_list[i].evenodd = i % 2; */
    if(g_proc_id == 0 && g_debug_level > 4) {
      printf("# Block id %d even odd id %d coordinate (%d %d %d %d)\n", 
             i, block_list[i].evenodd, block_list[i].coordinate[0], block_list[i].coordinate[1], 
             block_list[i].coordinate[2], block_list[i].coordinate[3]);
    }
    for (j = 0; j < g_N_s; j++) { 
      // write a zero element at the end of every spinor
      // this we need for boundary points, which we treat like this
      _spinor_null(block_list[i].basis[j][VOLUME/nb_blocks]);
    }

    if ((void*)(block_list[i].little_dirac_operator =       aligned_malloc_zero(9 * g_N_s * g_N_s * sizeof(_Complex double))) == NULL)
      CALLOC_ERROR_CRASH;
    if ((void*)(block_list[i].little_dirac_operator_32 =    aligned_malloc_zero(9 * g_N_s * g_N_s * sizeof(_Complex float))) == NULL)
      CALLOC_ERROR_CRASH;
    if ((void*)(block_list[i].little_dirac_operator_eo =    aligned_malloc_zero(9 * g_N_s * g_N_s * sizeof(_Complex double))) == NULL)
      CALLOC_ERROR_CRASH;
    if ((void*)(block_list[i].little_dirac_operator_eo_32 = aligned_malloc_zero(9 * g_N_s * g_N_s * sizeof(_Complex float))) == NULL)
      CALLOC_ERROR_CRASH;
  }
  if ((void*)(block_idx = calloc(8 * (VOLUME/nb_blocks), sizeof(int))) == NULL)
    CALLOC_ERROR_CRASH;
  
  if ((void*)(block_evenidx = calloc(8 * (VOLUME/nb_blocks/2), sizeof(int))) == NULL)
    CALLOC_ERROR_CRASH;
  
  if ((void*)(block_oddidx = calloc(8 * (VOLUME/nb_blocks/2), sizeof(int))) == NULL)
    CALLOC_ERROR_CRASH;
  
  init_blocks_geometry();
  init_blocks_gaugefield();
  etime = gettime();
  if(g_proc_id == 0 && g_debug_level > 0) {
    printf("# time for block initialisation %e s\n", etime - atime);
    fflush(stdout);
  }
  return 0;
}

int free_blocks() {
  int i;
  if(block_init == 1) {
    for(i = 0; i < nb_blocks; ++i) {
      free(block_list[i].basis);
      aligned_free(block_list[i].little_dirac_operator);
      aligned_free(block_list[i].little_dirac_operator_32);
      aligned_free(block_list[i].little_dirac_operator_eo);
      aligned_free(block_list[i].little_dirac_operator_eo_32);
    }
    free(block_ipt);
    free(bipt__);
    free(bipt_);
    free(bipt);
    free(index_block_eo);
    free(block_idx);
    free(block_evenidx);
    free(block_oddidx);
    aligned_free(u);
    aligned_free(u_32);
    aligned_free(basis);
    free(block_list);
    block_init = 0;
  }
  return 0;
}

int init_blocks_gaugefield() {
  /* 
     Copies the existing gauge field on the processor into the separate blocks in a form
     that is readable by the block Dirac operator. Specifically, in consecutive memory
     now +t,-t,+x,-x,+y,-y,+z,-z gauge links are stored. This requires double the storage in
     memory. 
  */

  int i, ix, ix_new = 0;

  for(int t = 0; t < dT;  t++) {
    for(int x = 0; x < dX; x++) {
      for(int y = 0; y < dY; y++) {
        for(int z = 0; z < dZ; z++) {
          i = 0;
          for(int bt = 0; bt < nblks_t; bt ++) {
            for(int bx = 0; bx < nblks_x; bx ++) {
              for(int by = 0; by < nblks_y; by ++) {
                for(int bz = 0; bz < nblks_z; bz ++) {
                  ix = g_ipt[t + bt*dT][x + bx*dX][y + by*dY][z + bz*dZ];
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
  blk_gauge_eo = 0;
  return(0);
}

int init_blocks_gaugefield_32() {
  /* 
     Copies the existing gauge field on the processor into the separate blocks in a form
     that is readable by the block Dirac operator. Specifically, in consecutive memory
     now +t,-t,+x,-x,+y,-y,+z,-z gauge links are stored. This requires double the storage in
     memory. 
  */

  int i, ix, ix_new = 0;

  for(int t = 0; t < dT;  t++) {
    for(int x = 0; x < dX; x++) {
      for(int y = 0; y < dY; y++) {
        for(int z = 0; z < dZ; z++) {
          i = 0;
          for(int bt = 0; bt < nblks_t; bt ++) {
            for(int bx = 0; bx < nblks_x; bx ++) {
              for(int by = 0; by < nblks_y; by ++) {
                for(int bz = 0; bz < nblks_z; bz ++) {
                  ix = g_ipt[t + bt*dT][x + bx*dX][y + by*dY][z + bz*dZ];
                  memcpy(block_list[i].u_32 + ix_new,     &g_gauge_field_32[ ix           ][0], sizeof(su3_32));
                  memcpy(block_list[i].u_32 + ix_new + 1, &g_gauge_field_32[ g_idn[ix][0] ][0], sizeof(su3_32));
                  memcpy(block_list[i].u_32 + ix_new + 2, &g_gauge_field_32[ ix           ][1], sizeof(su3_32));
                  memcpy(block_list[i].u_32 + ix_new + 3, &g_gauge_field_32[ g_idn[ix][1] ][1], sizeof(su3_32));
                  memcpy(block_list[i].u_32 + ix_new + 4, &g_gauge_field_32[ ix           ][2], sizeof(su3_32));
                  memcpy(block_list[i].u_32 + ix_new + 5, &g_gauge_field_32[ g_idn[ix][2] ][2], sizeof(su3_32));
                  memcpy(block_list[i].u_32 + ix_new + 6, &g_gauge_field_32[ ix           ][3], sizeof(su3_32));
                  memcpy(block_list[i].u_32 + ix_new + 7, &g_gauge_field_32[ g_idn[ix][3] ][3], sizeof(su3_32));
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
  blk_gauge_eo = 0;
  return(0);
}

int init_blocks_eo_gaugefield() {
  /* 
     Copies the existing gauge field on the processor into the separate blocks in a form
     that is readable by the block Hopping matrix. Specifically, in consecutive memory
     now +t,-t,+x,-x,+y,-y,+z,-z gauge links are stored. This requires double the storage in
     memory. 
  */

  int i, ix, ix_even = 0, ix_odd = (dT*dX*dY*dZ*8)/2, ixeo;

  for (int t = 0; t < dT;  t++) {
    for (int x = 0; x < dX; x++) {
      for (int y = 0; y < dY; y++) {
        for (int z = 0; z < dZ; z++) {
          if((t+x+y+z)%2 == 0) {
            ixeo = ix_even;
            ix_even += 8;
          }
          else {
            ixeo = ix_odd;
            ix_odd += 8;
          }
          i = 0;
          for(int bt = 0; bt < nblks_t; bt ++) {
            for(int bx = 0; bx < nblks_x; bx ++) {
              for(int by = 0; by < nblks_y; by ++) {
                for(int bz = 0; bz < nblks_z; bz ++) {
                  ix = g_ipt[t + bt*dT][x + bx*dX][y + by*dY][z + bz*dZ];
                  memcpy(block_list[i].u + ixeo,     &g_gauge_field[ ix           ][0], sizeof(su3));
                  memcpy(block_list[i].u + ixeo + 1, &g_gauge_field[ g_idn[ix][0] ][0], sizeof(su3));
                  memcpy(block_list[i].u + ixeo + 2, &g_gauge_field[ ix           ][1], sizeof(su3));
                  memcpy(block_list[i].u + ixeo + 3, &g_gauge_field[ g_idn[ix][1] ][1], sizeof(su3));
                  memcpy(block_list[i].u + ixeo + 4, &g_gauge_field[ ix           ][2], sizeof(su3));
                  memcpy(block_list[i].u + ixeo + 5, &g_gauge_field[ g_idn[ix][2] ][2], sizeof(su3));
                  memcpy(block_list[i].u + ixeo + 6, &g_gauge_field[ ix           ][3], sizeof(su3));
                  memcpy(block_list[i].u + ixeo + 7, &g_gauge_field[ g_idn[ix][3] ][3], sizeof(su3));
                  i++;
                }
              }
            }
          }
        }
      }
    }
  }
  blk_gauge_eo = 1;
  return(0);
}

int init_blocks_eo_gaugefield_32() {
  /* 
     Copies the existing gauge field on the processor into the separate blocks in a form
     that is readable by the block Hopping matrix. Specifically, in consecutive memory
     now +t,-t,+x,-x,+y,-y,+z,-z gauge links are stored. This requires double the storage in
     memory. 
  */

  int i, ix, ix_even = 0, ix_odd = (dT*dX*dY*dZ*8)/2, ixeo;

  for(int t = 0; t < dT;  t++) {
    for(int x = 0; x < dX; x++) {
      for(int y = 0; y < dY; y++) {
        for(int z = 0; z < dZ; z++) {
          if((t+x+y+z)%2 == 0) {
            ixeo = ix_even;
            ix_even += 8;
          }
          else {
            ixeo = ix_odd;
            ix_odd += 8;
          }
          i = 0;
          for(int bt = 0; bt < nblks_t; bt ++) {
            for(int bx = 0; bx < nblks_x; bx ++) {
              for(int by = 0; by < nblks_y; by ++) {
                for(int bz = 0; bz < nblks_z; bz ++) {
                  ix = g_ipt[t + bt*dT][x + bx*dX][y + by*dY][z + bz*dZ];
                  memcpy(block_list[i].u_32 + ixeo,     &g_gauge_field_32[ ix           ][0], sizeof(su3_32));
                  memcpy(block_list[i].u_32 + ixeo + 1, &g_gauge_field_32[ g_idn[ix][0] ][0], sizeof(su3_32));
                  memcpy(block_list[i].u_32 + ixeo + 2, &g_gauge_field_32[ ix           ][1], sizeof(su3_32));
                  memcpy(block_list[i].u_32 + ixeo + 3, &g_gauge_field_32[ g_idn[ix][1] ][1], sizeof(su3_32));
                  memcpy(block_list[i].u_32 + ixeo + 4, &g_gauge_field_32[ ix           ][2], sizeof(su3_32));
                  memcpy(block_list[i].u_32 + ixeo + 5, &g_gauge_field_32[ g_idn[ix][2] ][2], sizeof(su3_32));
                  memcpy(block_list[i].u_32 + ixeo + 6, &g_gauge_field_32[ ix           ][3], sizeof(su3_32));
                  memcpy(block_list[i].u_32 + ixeo + 7, &g_gauge_field_32[ g_idn[ix][3] ][3], sizeof(su3_32));
                  i++;
                }
              }
            }
          }
        }
      }
    }
  }
  blk_gauge_eo = 1;
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
  
  if(itest[blk->volume + blk->spinpad-1] != 2*(blk->BLX*blk->BLY*blk->BLZ+blk->BT*blk->BLX*blk->BLY+blk->BT*blk->BLY*blk->BLZ+blk->BT*blk->BLX*blk->BLZ)) {
    if(g_proc_id == 0){
      printf("error in block geometry, boundary points wrong %d != %d\n",
             itest[blk->volume + blk->spinpad-1], 2*(blk->BLX*blk->BLY*blk->BLZ+blk->BT*blk->BLX*blk->BLY+blk->BT*blk->BLY*blk->BLZ+blk->BT*blk->BLX*blk->BLZ));
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

  if(g_proc_id == 0 && g_debug_level > 4) {
    printf("# block geometry checked successfully for block %d !\n", blk->id);
  }
  for(i = 0; i < blk->volume; i++) {
    itest[i] = 0;
  }
  ipt = blk->evenidx;
  for(i = 0; i < 8*blk->volume/2; i++) {
    if(*ipt > (blk->volume/2 + blk->spinpad)-1 || *ipt < 0) {
      if(g_proc_id == 0) {
        printf("error in block eo geometry! ipt = %d dir = %d i = %d of %d\n",
               (*ipt), i%8, i/8, (blk->volume/2 + blk->spinpad));
      }
    }
    
    itest[*(ipt++)]++;
  }

  k = 0;
  for(i = 0; i < blk->volume/2; i++) {
    k += itest[i];
    if(itest[i] < 1 || itest[i] > 8) {
      if(g_proc_id == 0) {
        printf("error in block eo geometry, itest[%d] = %d\n", i, itest[i]);
      }
    }
  }
  k += itest[blk->volume/2 + blk->spinpad-1];
  if(k != 8*blk->volume/2) {
    if(g_proc_id == 0) {
      printf("error in block eo geometry, total number of points wrong %d != %d\n",
             k, 8*blk->volume/2);
    }
  }

  ipt = blk->evenidx;
  for(t = 0; t < T/nblks_t; t++) {
    for(x = 0; x < LX/nblks_x; x++) {
      for(y = 0; y < LY/nblks_y; y++) {
        for(z = 0; z < LZ/nblks_z; z++) {
          if((x + y + z + t)%2 == 0) {
            i = block_ipt[t][x][y][z]/2;
            if(t != T/nblks_t-1) {
              if(*ipt != block_ipt[t+1][x][y][z]/2 && g_proc_id == 0)
                printf("Shit +t! (%d %d %d %d): %d != %d at %d\n",
                       t, x, y, z, *ipt, block_ipt[t+1][x][y][z]/2, i);
            }
            else if(*ipt != VOLUME/nb_blocks/2)
              printf("Shit +t! (%d %d %d %d): %d != %d at %d\n",
                   t, x, y, z, *ipt, VOLUME/nb_blocks/2, i);
            ipt++;
            if(t != 0) {
              if(*ipt != block_ipt[t-1][x][y][z]/2 && g_proc_id == 0)
                printf("Shit -t! (%d %d %d %d): %d != %d at %d\n",
                       t, x, y, z, *ipt, block_ipt[t+1][x][y][z]/2, i);
            }
            else if(*ipt != VOLUME/nb_blocks/2)
              printf("Shit -t! (%d %d %d %d): %d != %d at %d\n",
                     t, x, y, z, *ipt, VOLUME/nb_blocks/2, i);
            ipt++;
            if(x != LX/nblks_x-1) {
              if(*ipt != block_ipt[t][x+1][y][z]/2 && g_proc_id == 0)
                printf("Shit +x! (%d %d %d %d): %d != %d at %d\n",
                       t, x, y, z, *ipt, block_ipt[t][x+1][y][z]/2, i);
            }
            else if(*ipt != VOLUME/nb_blocks/2)
              printf("Shit +x! (%d %d %d %d): %d != %d at %d\n",
                     t, x, y, z, *ipt, VOLUME/nb_blocks/2, i);
            ipt++;
            if(x != 0) {
              if(*ipt != block_ipt[t][x-1][y][z]/2 && g_proc_id == 0)
                printf("Shit -x! (%d %d %d %d): %d != %d at %d\n",
                       t, x, y, z, *ipt, block_ipt[t][x-1][y][z]/2, i);
            }
            else if(*ipt != VOLUME/nb_blocks/2)
              printf("Shit -x! (%d %d %d %d): %d != %d at %d\n",
                     t, x, y, z, *ipt, VOLUME/nb_blocks, i);
            ipt++;
            if(y != LY/nblks_y-1) {
              if(*ipt != block_ipt[t][x][y+1][z]/2 && g_proc_id == 0)
                printf("Shit +y! (%d %d %d %d): %d != %d at %d\n",
                       t, x, y, z, *ipt, block_ipt[t][x][y+1][z]/2, i);
            }
            else if(*ipt != VOLUME/nb_blocks/2)
              printf("Shit +y! (%d %d %d %d): %d != %d at %d\n",
                     t, x, y, z, *ipt, VOLUME/nb_blocks/2, i);
            ipt++;
            if(y != 0) {
              if(*ipt != block_ipt[t][x][y-1][z]/2 && g_proc_id == 0)
                printf("Shit -y! (%d %d %d %d): %d != %d at %d\n",
                       t, x, y, z, *ipt, block_ipt[t][x][y-1][z]/2, i);
            }
            else if(*ipt != VOLUME/nb_blocks/2)
              printf("Shit -y! (%d %d %d %d): %d != %d at %d\n",
                     t, x, y, z, *ipt, VOLUME/nb_blocks/2, i);
            ipt++;
            if(z != LZ/nblks_z-1) {
              if(*ipt != block_ipt[t][x][y][z+1]/2 && g_proc_id == 0)
                printf("Shit +z! (%d %d %d %d): %d != %d at %d\n",
                       t, x, y, z, *ipt, block_ipt[t][x][y][z+1]/2, i);
            }
            else if(*ipt != VOLUME/nb_blocks/2)
              printf("Shit +z! (%d %d %d %d): %d != %d at %d\n",
                     t, x, y, z, *ipt, VOLUME/nb_blocks/2, i);
            ipt++;
            if(z != 0) {
              if(*ipt != block_ipt[t][x][y][z-1]/2 && g_proc_id == 0)
                printf("Shit -z! (%d %d %d %d): %d != %d at %d\n",
                       t, x, y, z, *ipt, block_ipt[t][x][y][z-1]/2, i);
            }
            else if(*ipt != VOLUME/nb_blocks/2)
              printf("Shit -z! (%d %d %d %d): %d != %d at %d\n",
                     t, x, y, z, *ipt, VOLUME/nb_blocks/2, i);
            ipt++;
          }
        }
      }
    }
  }

  if(g_proc_id == 0 && g_debug_level > 4) {
    printf("# block eo geometry checked successfully for block %d !\n", blk->id);
  }

  free(itest);
  return(0);
}

int init_blocks_geometry() {
  int i_even, i_odd;
  int zstride = 1;
  int ystride = dZ;
  int xstride = dY * dZ;
  int tstride = dX * dY * dZ;
  int boundidx = VOLUME/nb_blocks;
  for (int ix = 0; ix < VOLUME/nb_blocks; ++ix) {
    block_idx[8 * ix + 0] = ix           >= VOLUME/nb_blocks - tstride  ? boundidx : ix + tstride;/* +t */
    block_idx[8 * ix + 1] = ix           <  tstride                     ? boundidx : ix - tstride;/* -t */
    block_idx[8 * ix + 2] = ix % tstride >= dZ * dY * (dX - 1)          ? boundidx : ix + xstride;/* +x */
    block_idx[8 * ix + 3] = ix % tstride <  dZ * dY                     ? boundidx : ix - xstride;/* -x */
    block_idx[8 * ix + 4] = ix % xstride >= dZ * (dY - 1)               ? boundidx : ix + ystride;/* +y */
    block_idx[8 * ix + 5] = ix % xstride <  dZ                          ? boundidx : ix - ystride;/* -y */
    block_idx[8 * ix + 6] = ix % ystride == dZ - 1                      ? boundidx : ix + zstride;/* +z */
    block_idx[8 * ix + 7] = ix % ystride == 0                           ? boundidx : ix - zstride;/* -z */
    /* Assume that all directions have even extension */
    /* even and odd versions should be equal          */
    int eo = ((ix%dZ)+(ix/ystride)%dY+(ix/(xstride))%dX
          +ix/(tstride))%2;
    if(eo == 0) {
      block_evenidx[8*(ix/2) + 0] = block_idx[8 * ix + 0] / 2;
      block_evenidx[8*(ix/2) + 1] = block_idx[8 * ix + 1] / 2;
      block_evenidx[8*(ix/2) + 2] = block_idx[8 * ix + 2] / 2;
      block_evenidx[8*(ix/2) + 3] = block_idx[8 * ix + 3] / 2;
      block_evenidx[8*(ix/2) + 4] = block_idx[8 * ix + 4] / 2;
      block_evenidx[8*(ix/2) + 5] = block_idx[8 * ix + 5] / 2;
      block_evenidx[8*(ix/2) + 6] = block_idx[8 * ix + 6] / 2;
      block_evenidx[8*(ix/2) + 7] = block_idx[8 * ix + 7] / 2;
    }
    else {
      block_oddidx[8*(ix/2) + 0] = block_idx[8 * ix + 0] / 2;
      block_oddidx[8*(ix/2) + 1] = block_idx[8 * ix + 1] / 2;
      block_oddidx[8*(ix/2) + 2] = block_idx[8 * ix + 2] / 2;
      block_oddidx[8*(ix/2) + 3] = block_idx[8 * ix + 3] / 2;
      block_oddidx[8*(ix/2) + 4] = block_idx[8 * ix + 4] / 2;
      block_oddidx[8*(ix/2) + 5] = block_idx[8 * ix + 5] / 2;
      block_oddidx[8*(ix/2) + 6] = block_idx[8 * ix + 6] / 2;
      block_oddidx[8*(ix/2) + 7] = block_idx[8 * ix + 7] / 2;
    }
  }
  for(int i = 0; i < nb_blocks; i++) {
    block_list[i].idx = block_idx;
    block_list[i].evenidx = block_evenidx;
    block_list[i].oddidx = block_oddidx;
  }
  int ix = 0;
  for(int t = 0; t < dT; t++) {
    for(int x = 0; x < dX; x++) {
      for(int y = 0; y < dY; y++) {
        for(int z = 0; z < dZ; z++) {
          block_ipt[t][x][y][z] = ix;
          ix++;
        }
      }
    }
  }

  i_even = 0;
  i_odd = 0;
  for (int t = 0; t < nblks_t; t++) {
    for (int x = 0; x < nblks_x; x++) {
      for (int y = 0; y < nblks_y; y++) {
        for (int z = 0; z < nblks_z; z++) {
          ix = block_index(t,x,y,z);
          if (block_list[ix].evenodd == 0) {
            index_block_eo[ix] = i_even;
            i_even++;
          }
          else {
            index_block_eo[ix] = i_odd;
            i_odd++;
          }
        }
      }
    }
  }

  for(int ix = 0; ix < nb_blocks; ix++) {
    zstride = check_blocks_geometry(&block_list[ix]);
  }

  return 0;
}

/* Uses a Modified Gram-Schmidt algorithm to orthonormalize little basis vectors */
void block_orthonormalize(block *parent) {
  int i, j;
  _Complex double coeff;
  double scale;

  for(i = 0; i < g_N_s; ++i){
    /* rescale the current vector */
    scale = 1. / sqrt(square_norm(parent->basis[i], parent->volume, 0));
    mul_r(parent->basis[i], scale, parent->basis[i], parent->volume);

    /* rescaling done, now subtract this direction from all vectors that follow */
    for(j = i + 1; j < g_N_s; ++j){
      coeff = scalar_prod(parent->basis[i], parent->basis[j], parent->volume, 0);
      assign_diff_mul(parent->basis[j], parent->basis[i], coeff, parent->volume);
    }
  }

  if(g_debug_level > 4) {
    for(i = 0; i < g_N_s; i++) {
      for(j = 0; j < g_N_s; j++) {
        coeff = scalar_prod(parent->basis[i], parent->basis[j], parent->volume, 0);
        if(g_proc_id == 0) printf("basis id = %d <%d, %d> = %1.3e +i %1.3e\n", parent->id, j, i, creal(coeff), cimag(coeff));
      }
    }
  }
  return;
}

void block_orthonormalize_free(block *parent) {
  int i, j;
  _Complex double coeff;
  double scale;

  for(i = 0; i < 12; i++){  // CHECK THIS !!!!!! 12
    /* rescale the current vector */
    constant_spinor_field(parent->basis[i], i, parent->volume);
    scale = 1. / sqrt(square_norm(parent->basis[i], parent->volume, 0));
    mul_r(parent->basis[i], scale, parent->basis[i], parent->volume);
  }

  if(g_debug_level > 4 && g_proc_id == 0) {
    for(i = 0; i < g_N_s; i++) {
      for(j = 0; j < g_N_s; j++) {
        coeff = scalar_prod(parent->basis[i], parent->basis[j], parent->volume, 0);
        if(g_proc_id == 0) printf("basis id = %d <%d, %d> = %1.3e +i %1.3e\n", parent->id, j, i, creal(coeff), cimag(coeff));
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
  for(l = 0; l < g_N_s; ++l) {
    block_list[idx].little_dirac_operator[dir * g_N_s * g_N_s + vecnum * g_N_s + l] =
      scalar_prod(block_list[idx].basis[l], psi + idx * (VOLUME/nb_blocks+1), VOLUME/nb_blocks, 0);
  }
}



/* checked CU */
void compute_little_D_diagonal(const int mul_g5) {
  int i,j, blk;
  spinor * tmp, * _tmp;
  _Complex double * M;
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
      if(mul_g5) gamma5(tmp, tmp, block_list[blk].volume);
      for(j = 0; j < g_N_s; j++) {
        M[i * g_N_s + j]  = scalar_prod(block_list[blk].basis[j], tmp, block_list[blk].volume, 0);
        block_list[blk].little_dirac_operator_32[i*g_N_s + j] = (_Complex float)M[i * g_N_s + j];
      }
    }
  }

  if(g_debug_level > 2) {
    if (g_N_s <= 5 && !g_cart_id){
      printf("\n\n  *** CHECKING LITTLE D ***\n");
      printf("\n  ** node 0, lower block **\n");
      for (i = 0*g_N_s; i < 9 * g_N_s; ++i){
        printf(" [ ");
        for (j = 0; j < g_N_s; ++j){
          printf("%s%1.3e %s %1.3e i", creal(block_list[0].little_dirac_operator[i * g_N_s + j]) >= 0 ? "  " : "- ", 
                 creal(block_list[0].little_dirac_operator[i * g_N_s + j]) >= 0 ? creal(block_list[0].little_dirac_operator[i * g_N_s + j]) : -creal(block_list[0].little_dirac_operator[i * g_N_s + j]), 
                 cimag(block_list[0].little_dirac_operator[i * g_N_s + j]) >= 0 ? "+" : "-", 
                 cimag(block_list[0].little_dirac_operator[i * g_N_s + j]) >= 0 ? cimag(block_list[0].little_dirac_operator[i * g_N_s + j]) : -cimag(block_list[0].little_dirac_operator[i * g_N_s + j]));
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
          printf("%s%1.3e %s %1.3e i", creal(block_list[1].little_dirac_operator[i * g_N_s + j]) >= 0 ? "  " : "- ", 
                 creal(block_list[1].little_dirac_operator[i * g_N_s + j]) >= 0 ? creal(block_list[1].little_dirac_operator[i * g_N_s + j]) : -creal(block_list[1].little_dirac_operator[i * g_N_s + j]), 
                 cimag(block_list[1].little_dirac_operator[i * g_N_s + j]) >= 0 ? "+" : "-", 
                 cimag(block_list[1].little_dirac_operator[i * g_N_s + j]) >= 0 ? cimag(block_list[1].little_dirac_operator[i * g_N_s + j]) : -cimag(block_list[1].little_dirac_operator[i * g_N_s + j]));
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


  free(_tmp);
  return;
}


/* checked CU in 4d parallel case */
void compute_little_D(const int mul_g5) {
  /* 
     This is the little dirac routine rewritten according to multidimensional blocking
     Adaptation by Claude Tadonki (claude.tadonki@u-psud.fr)
     Date: May 2010
  */
  spinor *scratch, *temp, *_scratch;
  int mu=0;
  double atime, etime;
  // the block volume
  int bvol = block_list[1].volume;
  int t_start, t_end, x_start, x_end, y_start, y_end, z_start, z_end;
  //_Complex double c;
  //int count=0;
  int is_up = 0;
  int dT, dX, dY, dZ;
  double musave = g_mu;
  double kappasave = g_kappa;
  if(kappa_dfl > 0) {
    g_kappa = kappa_dflgen;
  }
  if(mu_dfl > -10) {
    g_mu = mu_dfl;
    if(g_mu*musave < 0) g_mu *= -1.;
  }
  boundary(g_kappa);

  dT = T/nblks_t; 
  dX = LX/nblks_x; 
  dY = LY/nblks_y; 
  dZ = LZ/nblks_z;

  if(g_proc_id == 0 && g_debug_level > 1) {
    printf("# compute_little_D called with mul_g5 = %d\n", mul_g5);
    printf("# compute_little_D parameters mu= %e, kappa= %e\n", g_mu/2./g_kappa, g_kappa);
  }

  /* for a full spinor field we need VOLUMEPLUSRAND                 */
  /* because we use the same geometry as for the                    */
  /* gauge field                                                    */
  /* It is VOLUME + 2*LZ*(LY*LX + T*LY + T*LX) + 4*LZ*(LY + T + LX) */
  if( (_scratch = calloc(2*VOLUMEPLUSRAND+1, sizeof(spinor))) == NULL) {
    fprintf(stderr, "not enough memory for scratch in compute_little_D! Aborting...\n");
    exit(-1);
  }
  scratch = (spinor*)(((unsigned long int)(_scratch)+ALIGN_BASE)&~ALIGN_BASE);
  // not needed?
  temp = scratch + VOLUMEPLUSRAND;
  // NEEDs TO BE REWRITTEN
  atime = gettime();
#ifdef OMP
#pragma omp parallel
  {
#endif
    spinor * bscratch;
    _Complex double * M;
#ifdef OMP
#pragma omp for
#endif
  for(int blk = 0; blk < nb_blocks; blk++) {
    bscratch = scratch + blk*bvol;
    M = block_list[blk].little_dirac_operator;
    for(int i = 0; i < g_N_s; i++) {
      Block_D_psi(&block_list[blk], bscratch, block_list[blk].basis[i]);
      if(mul_g5) gamma5(bscratch, bscratch, bvol);
      for(int j = 0; j < g_N_s; j++) {
        M[i * g_N_s + j]  = scalar_prod_ts(block_list[blk].basis[j], bscratch, bvol, 0);
        block_list[blk].little_dirac_operator_32[i * g_N_s + j] = (_Complex float)M[i * g_N_s + j];
        
        block_list[block_list[blk].evenodd_id].little_dirac_operator_eo[i * g_N_s + j] = M[i * g_N_s + j];
      }
    }
  }
  /* computation of little_Dhat^{-1}_ee */
#ifdef OMP
#pragma omp for
#endif
  for(int blk = 0; blk < nb_blocks/2; blk++) {
    LUInvert(g_N_s, block_list[blk].little_dirac_operator_eo, g_N_s);
  }
#ifdef OMP
  } /* OpenMP closing brace */
#endif
  for (int i = 0; i < g_N_s; i++) {
    reconstruct_global_field_GEN_ID(scratch, block_list, i , nb_blocks);
    
#ifdef MPI
    xchange_lexicfield(scratch);
#endif
    
    /* the initialisation causes troubles on a single processor */
    if(g_nproc == -1) zero_spinor_field(scratch, VOLUME);
    /* +-t +-x +-y +-z */

    for(int pm = 0; pm < 8; pm++) {
      /* We set up the generic bounds */
      t_start = 0; t_end = dT;
      x_start = 0; x_end = dX;
      y_start = 0; y_end = dY;
      z_start = 0; z_end = dZ;
      switch(pm) { 
      case 0: t_start = dT - 1; t_end = t_start + 1; mu = 0; is_up = 1; break; // Boundary in dir +t
      case 1: t_start = 0;      t_end = t_start + 1; mu = 0; is_up = 0; break; // Boundary in dir -t
      case 2: x_start = dX - 1; x_end = x_start + 1; mu = 1; is_up = 1; break; // Boundary in dir +x
      case 3: x_start = 0;      x_end = x_start + 1; mu = 1; is_up = 0; break; // Boundary in dir -x
      case 4: y_start = dY - 1; y_end = y_start + 1; mu = 2; is_up = 1; break; // Boundary in dir +y
      case 5: y_start = 0;      y_end = y_start + 1; mu = 2; is_up = 0; break; // Boundary in dir -y
      case 6: z_start = dZ - 1; z_end = z_start + 1; mu = 3; is_up = 1; break; // Boundary in dir +z
      case 7: z_start = 0;      z_end = z_start + 1; mu = 3; is_up = 0; break; // Boundary in dir -z
      default: ;
      }
      /* Dirac operator on the boundaries */
#ifdef OMP
#pragma omp parallel 
      {
#endif
        spinor r;
        spinor * s = NULL;
        su3 * u;
        int ib, iy, ix, iz;
        int bx, by, bz, bt;
	int block_id_eo;
        _Complex double c;
#ifdef OMP
#pragma omp for
#endif
        for(int block_id = 0; block_id < nb_blocks; block_id++) {
          block_id_eo =  block_list[block_id].evenodd_id;

          for(int j = 0; j < g_N_s; j++) {
            iz = i * g_N_s + j  + (pm + 1) * g_N_s * g_N_s;
            block_list[block_id].little_dirac_operator[ iz ] = 0.0;
            block_list[block_id_eo].little_dirac_operator_eo[ iz ] = 0.0;
          }

          s = NULL;
          iy = 0; ix = 0;
          bz = block_id % nblks_z;
          by = ((block_id - bz) % (nblks_y*nblks_z)) / nblks_z;
          bx = ((block_id - bz - by*nblks_z) % (nblks_x*nblks_y*nblks_z)) / ( nblks_z*nblks_y);
          bt = (block_id - bz - by*nblks_z - bx*nblks_z*nblks_y) / ( nblks_z*nblks_y*nblks_x);
          for(int t = t_start; t < t_end; t++) {
            for(int x = x_start; x < x_end; x++) {
              for(int y = y_start; y < y_end; y++) {
                for(int z = z_start; z < z_end; z++) {
                  // We treat the case when we need to cross between blocks
                  // We are in block (bt, bx, by, bz) and compute direction pm
                  // We check inner block statement by ( b_ > 0 )&&( b_ < nblks_ - 1 )
                  // Other cases are treated in a standard way using the boundary of the scracth array 
                  ib = -1; // ib is the index of the selected block if any
                  if((pm==0)&&(bt<nblks_t-1)&&(t==t_end-1)){ //direction +t
                    iy = index_b(0, x, y, z); /* lowest edge of upper block needed */
                    ib = block_index(bt+1, bx, by, bz);
                  }
                  else if((pm==1)&&(bt>0)&&(t==0)){ //direction -t
                    iy = index_b(dT - 1, x, y, z); /* highest edge of lower block needed */
                    ib = block_index(bt-1, bx, by, bz);
                  }
                  else if((pm==2)&&(bx<nblks_x-1)&&(x==x_end-1)){ //direction +x
                    iy = index_b(t, 0, y, z); /* lowest edge of upper block needed */
                    ib = block_index(bt, bx+1, by, bz);
                  }
                  else if((pm==3)&&(bx>0)&&(x==0)){ //direction -x
                    iy = index_b(t, dX - 1, y, z); /* highest edge of lower block needed */
                    ib = block_index(bt, bx-1, by, bz);
                  }
                  else if((pm==4)&&(by<nblks_y-1)&&(y==y_end-1)){ //direction +y
                    iy = index_b(t, x, 0, z); /* lowest edge of upper block needed */
                    ib = block_index(bt, bx, by+1, bz);
                  }
                  else if((pm==5)&&(by>0)&&(y==0)){ //direction -y
                    iy = index_b(t, x, dY - 1, z); /* highest edge of lower block needed */
                    ib = block_index(bt, bx, by-1, bz);
                  }
                  else if((pm==6)&&(bz<nblks_z-1)&&(z==z_end-1)){ //direction +z
                    iy = index_b(t, x, y, 0); /* lowest edge of upper block needed */
                    ib = block_index(bt, bx, by, bz+1);
                  }
                  else if((pm==7)&&(bz>0)&&(z==0)){ //direction -z
                    iy = index_b(t, x, y, dZ - 1); /* highest edge of lower block needed */
                    ib = block_index(bt, bx, by, bz-1);
                  }
                  ix = index_a(dT*bt + t, dX*bx + x, dY*by + y, dZ*bz + z);// GAFFE ICI
                  if(is_up == 1) {
                    s = &scratch[ g_iup[ ix ][mu] ]; 
                    u = &g_gauge_field[ ix ][mu];
                  }
                  else {
                    s = &scratch[ g_idn[ ix ][mu] ];
                    u = &g_gauge_field[ g_idn[ix][mu] ][mu];
                  }
                  if(ib >= 0) s = &block_list[ib].basis[ i ][ iy ] ; 
                  boundary_D[pm](&r, s, u);
                  if(mul_g5) gamma5(&r, &r, 1);
                  // now all the scalar products
                  for(int j = 0; j < g_N_s; j++) {
                    iz = i * g_N_s + j  + (pm + 1) * g_N_s * g_N_s;
                    ix = index_b(t, x, y, z); 
                    s = &block_list[block_id].basis[j][ ix ];
                    c = scalar_prod_ts(s, &r, 1, 0);
                    block_list[block_id].little_dirac_operator[ iz ] += c;
                    block_list[block_id_eo].little_dirac_operator_eo[ iz ] += c;
                  }
                }
              }
            }
          }
        }
#ifdef OMP
      } // OMP closing brace
#endif
    }
  }
#ifdef OMP
#pragma omp parallel for
#endif
  for(int ij = 0; ij < nb_blocks*9*g_N_s*g_N_s; ij++) {
    int i = ij / (9*g_N_s*g_N_s);
    int j = ij % (9*g_N_s*g_N_s);
    block_list[i].little_dirac_operator_32[ j ] 
      = (_Complex float)block_list[i].little_dirac_operator[ j ];
    block_list[i].little_dirac_operator_eo_32[ j ] 
      = (_Complex float)block_list[i].little_dirac_operator_eo[ j ];
  }
  etime = gettime();
  if(g_debug_level > 2 && g_proc_id == 0) {
    printf("# time for compute_little_D: %e\n", etime-atime);
  }

  if(g_debug_level > 2) {
    if (g_N_s <= 5 && !g_cart_id){
      printf("\n\n  *** CHECKING LITTLE D ***\n");
      printf("\n  ** node 0, lower block **\n");
      for (int i = 0*g_N_s; i < 9 * g_N_s; ++i){
        printf(" [ ");
        for (int j = 0; j < g_N_s; ++j){
          printf("%s%1.3e %s %1.3e i", 
                 creal(block_list[0].little_dirac_operator[i * g_N_s + j]) >= 0 ? "  " : "- ", 
                 creal(block_list[0].little_dirac_operator[i * g_N_s + j]) >= 0 ? creal(block_list[0].little_dirac_operator[i * g_N_s + j]) : -creal(block_list[0].little_dirac_operator[i * g_N_s + j]), 
                 cimag(block_list[0].little_dirac_operator[i * g_N_s + j]) >= 0 ? "+" : "-", 
                 cimag(block_list[0].little_dirac_operator[i * g_N_s + j]) >= 0 ? cimag(block_list[0].little_dirac_operator[i * g_N_s + j]) : -cimag(block_list[0].little_dirac_operator[i * g_N_s + j]));
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
      for (int i = 0*g_N_s; i < 9 * g_N_s; ++i){
        printf(" [ ");
        for (int j = 0; j < g_N_s; ++j){
          printf("%s%1.3e %s %1.3e i", 
                 creal(block_list[1].little_dirac_operator[i * g_N_s + j]) >= 0 ? "  " : "- ", 
                 creal(block_list[1].little_dirac_operator[i * g_N_s + j]) >= 0 ? creal(block_list[1].little_dirac_operator[i * g_N_s + j]) : -creal(block_list[1].little_dirac_operator[i * g_N_s + j]), 
                 cimag(block_list[1].little_dirac_operator[i * g_N_s + j]) >= 0 ? "+" : "-", 
                 cimag(block_list[1].little_dirac_operator[i * g_N_s + j]) >= 0 ? cimag(block_list[1].little_dirac_operator[i * g_N_s + j]) : -cimag(block_list[1].little_dirac_operator[i * g_N_s + j]));
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
  g_mu = musave;
  g_kappa = kappasave;
  boundary(g_kappa);
  
  free(_scratch);
  return;
}


int split_global_field_GEN(spinor ** const psi, spinor * const field, const int nb_blocks) {
  int j,ctr_t=0;
  int x, y, z, t;
  int bx, by, bz, bt, block_id;
  for (t = 0; t < dT; t++) {
    for (x = 0; x < dX; x++) {
      for (y = 0; y < dY; y++) {
        for (z = 0; z < dZ; z++) {
          block_id = 0;
          for(bt = 0; bt < nblks_t; bt++) {
            for(bx = 0; bx < nblks_x; bx++) {
              for(by = 0; by < nblks_y; by++) {
                for(bz = 0; bz < nblks_z; bz++) {
                  _spinor_assign(*(psi[block_id] + ctr_t), 
                                 *(field + index_a(dT*bt + t, dX*bx + x, dY*by + y, dZ*bz + z)));
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
    for(j = 0; j < nb_blocks; j++)
      printf("Basis norm %2d = %1.3e\n", j, square_norm(psi[j],  VOLUME / nb_blocks, 0));
  }
  return 0;
}

int split_global_field_GEN_ID(block * const block_list, const int id, spinor * const field, const int nb_blocks){
  int j,ctr_t=0;
  int x, y, z, t;
  int bx, by, bz, bt, block_id;
  for (t = 0; t < dT; t++) {
    for (x = 0; x < dX; x++) {
      for (y = 0; y < dY; y++) {
        for (z = 0; z < dZ; z++) {
          block_id = 0;
          for(bt = 0; bt < nblks_t; bt++) {
            for(bx = 0; bx < nblks_x; bx++) {
              for(by = 0; by < nblks_y; by++) {
                for(bz = 0; bz < nblks_z; bz++) {
                  _spinor_assign(*(block_list[block_id].basis[id] + ctr_t), 
                                 *(field + index_a(dT*bt + t, dX*bx + x, dY*by + y, dZ*bz + z)));
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
    for(j = 0; j < nb_blocks; j++)
      printf("Basis norm %2d = %1.3e\n", j, square_norm(block_list[j].basis[id],  VOLUME / nb_blocks, 0));
  }
  return 0;
}

/* copies the part of globalfields corresponding to block blk */
/* to the block field blockfield                              */
void copy_global_to_block(spinor * const blockfield, spinor * const globalfield, const int blk) {
  int it,ix,iy,iz;
  int ibt,ibx,iby,ibz;
  int itb,ixb,iyb,izb;
  int ixcurrent;
  
  ibz = blk%nblks_z;
  iby = (blk / nblks_z)%nblks_y;
  ibx = (blk / (nblks_y * nblks_z))%nblks_x;
  ibt = blk / (nblks_x * nblks_y*nblks_z);

  ixcurrent=0;
  // FIXME: here we should better run only through the block volume!
  for (int i = 0; i < VOLUME; i++) {

    /* global coordinates */
    iz = i%LZ;
    iy = (i / LZ)%LY;
    ix = (i / (LY * LZ))%LX;
    it = i / (LX * LY * LZ);

    /* block coordinates */
    izb = iz / block_list[blk].BLZ;
    iyb = iy / block_list[blk].BLY;
    ixb = ix / block_list[blk].BLX;
    itb = it / block_list[blk].BT;
    
    if ((ibz == izb) && (iby == iyb) && (ibx == ixb) && (ibt==itb)) {
      memcpy(blockfield+ixcurrent, globalfield+i, sizeof(spinor));
      ixcurrent++;
    }
  }
  return;
}

/* copies the part of globalfields corresponding to block blk */
/* to the even and odd block fields                           */
void copy_global_to_block_eo(spinor * const beven, spinor * const bodd, spinor * const globalfield, const int blk) {
  int t, x, y, z;
  int i,it,ix,iy,iz;
  int even = 0, odd = 0;
  
  for(t = 0; t < block_list[blk].BT; t++) {
    it = t + block_list[blk].mpilocal_coordinate[0]*block_list[blk].BT;
    for(x = 0; x < block_list[blk].BLX; x++) {
      ix = x +  block_list[blk].mpilocal_coordinate[1]*block_list[blk].BLX;
      for(y = 0; y < block_list[blk].BLY; y++) {
        iy = y +  block_list[blk].mpilocal_coordinate[2]*block_list[blk].BLY;
        for(z = 0; z < block_list[blk].BLZ; z++) {
          iz = z +  block_list[blk].mpilocal_coordinate[3]*block_list[blk].BLZ;
          i = g_ipt[it][ix][iy][iz];
          if((t+x+y+z)%2 == 0) {
            memcpy(beven + even, globalfield + i, sizeof(spinor));
            even++;
          }
          else {
            memcpy(bodd + odd, globalfield + i, sizeof(spinor));
            odd++;
          }
        }
      }
    }
  }
  return;
}

void copy_global_to_block_eo_32(spinor32 * const beven, spinor32 * const bodd, 
                                spinor * const globalfield, const int blk) {
  int i,it,ix,iy,iz;
  int even = 0, odd = 0;
  _Complex float * to = NULL;
  _Complex double * from = NULL;
  spinor32 * tmp = NULL;
  
  for(int t = 0; t < block_list[blk].BT; t++) {
    it = t + block_list[blk].mpilocal_coordinate[0]*block_list[blk].BT;
    for(int x = 0; x < block_list[blk].BLX; x++) {
      ix = x +  block_list[blk].mpilocal_coordinate[1]*block_list[blk].BLX;
      for(int y = 0; y < block_list[blk].BLY; y++) {
        iy = y +  block_list[blk].mpilocal_coordinate[2]*block_list[blk].BLY;
        for(int z = 0; z < block_list[blk].BLZ; z++) {
          iz = z +  block_list[blk].mpilocal_coordinate[3]*block_list[blk].BLZ;
          i = g_ipt[it][ix][iy][iz];

          if((t+x+y+z)%2 == 0) {
            tmp = beven + even;
            even++;
          }
          else {
            tmp = bodd + odd;
            odd++;
          }

          to = (_Complex float*) tmp;
          from = (_Complex double*) (globalfield + i);
          for(int k = 0; k < 12; k++) {
            (*to) = (_Complex float) (*from);
            to++;
            from++;
          }
        }
      }
    }
  }
  return;
}

/* reverts copy_global_to_block_eo */
void copy_block_eo_to_global(spinor * const globalfield, spinor * const beven, spinor * const bodd, const int blk) {
  int t, x, y, z;
  int i,it,ix,iy,iz;
  int even = 0, odd = 0;
  
  for(t = 0; t < block_list[blk].BT; t++) {
    it = t + block_list[blk].mpilocal_coordinate[0]*block_list[blk].BT;
    for(x = 0; x < block_list[blk].BLX; x++) {
      ix = x +  block_list[blk].mpilocal_coordinate[1]*block_list[blk].BLX;
      for(y = 0; y < block_list[blk].BLY; y++) {
        iy = y +  block_list[blk].mpilocal_coordinate[2]*block_list[blk].BLY;
        for(z = 0; z < block_list[blk].BLZ; z++) {
          iz = z +  block_list[blk].mpilocal_coordinate[3]*block_list[blk].BLZ;
          i = g_ipt[it][ix][iy][iz];
          if((t+x+y+z)%2 == 0) {
            memcpy(globalfield + i, beven + even, sizeof(spinor));
            even++;
          }
          else {
            memcpy(globalfield + i, bodd + odd, sizeof(spinor));
            odd++;
          }
        }
      }
    }
  }
  return;
}


/* reconstructs the parts of globalfield corresponding to block blk */
/* from block field blockfield                                      */
void copy_block_to_global(spinor * const globalfield, spinor * const blockfield, const int blk) {
  int it,ix,iy,iz;
  int ibt,ibx,iby,ibz;
  int itb,ixb,iyb,izb;
  int ixcurrent;
  
  ibz = blk%nblks_z;
  iby = (blk / nblks_z)%nblks_y;
  ibx = (blk / (nblks_y * nblks_z))%nblks_x;
  ibt = blk / (nblks_x * nblks_y*nblks_z);

  ixcurrent=0;
  // FIXME: here we should better run only through the block volume!
  for (int i = 0; i < VOLUME; i++) {

    /* global coordinates */
    iz = i%LZ;
    iy = (i / LZ)%LY;
    ix = (i / (LY * LZ))%LX;
    it = i / (LX * LY * LZ);

    /* block coordinates */
    izb = iz / block_list[blk].BLZ;
    iyb = iy / block_list[blk].BLY;
    ixb = ix / block_list[blk].BLX;
    itb = it / block_list[blk].BT;
    
    if ((ibz == izb) && (iby == iyb) && (ibx == ixb) && (ibt==itb)) {
      memcpy(globalfield+i, blockfield+ixcurrent, sizeof(spinor));
      ixcurrent++;
    }
  }

  return;
}




/* Reconstructs a global field from the little basis of nb_blocks blocks */
void reconstruct_global_field_GEN(spinor * const rec_field, spinor ** const psi, const int nb_blocks) {

#ifdef OMP
#pragma omp parallel for
#endif
  for(int ix = 0; ix < dT*dX*dY*dZ; ix++) {
    int z = ix % dZ;
    int y = ((ix - z) % (dY*dZ)) / dZ;
    int x = (ix - z - y*dZ) % (dX*dY*dZ) / (dZ*dY);
    int t = (ix - z - y*dZ - x*dZ*dY) / (dZ*dY*dX);
    for(int block_id = 0; block_id < nb_blocks; block_id++) {
      int bz = block_id % nblks_z;
      int by = ((block_id - bz) % (nblks_y*nblks_z)) / nblks_z;
      int bx = ((block_id - bz - by*nblks_z) % (nblks_x*nblks_y*nblks_z)) / ( nblks_z*nblks_y);
      int bt = (block_id - bz - by*nblks_z - bx*nblks_z*nblks_y) / ( nblks_z*nblks_y*nblks_x);
      _spinor_assign(*(rec_field + index_a(dT*bt + t, dX*bx + x, dY*by + y, dZ*bz + z)), 
                     *(psi[block_id] + ix));
    }
  }
  return;
}

/* Reconstructs a global field from the little basis of nb_blocks blocks taken from block_list[*].basis[id] */
void reconstruct_global_field_GEN_ID(spinor * const rec_field, block * const block_list, const int id, const int nb_blocks) {
#ifdef OMP
#pragma omp parallel for
#endif
  for(int ix = 0; ix < dT*dX*dY*dZ; ix++) {
    int z = ix % dZ;
    int y = ((ix - z) % (dY*dZ)) / dZ;
    int x = (ix - z - y*dZ) % (dX*dY*dZ) / (dZ*dY);
    int t = (ix - z - y*dZ - x*dZ*dY) / (dZ*dY*dX);
    for(int block_id = 0; block_id < nb_blocks; block_id++) {
      int bz = block_id % nblks_z;
      int by = ((block_id - bz) % (nblks_y*nblks_z)) / nblks_z;
      int bx = ((block_id - bz - by*nblks_z) % (nblks_x*nblks_y*nblks_z)) / ( nblks_z*nblks_y);
      int bt = (block_id - bz - by*nblks_z - bx*nblks_z*nblks_y) / ( nblks_z*nblks_y*nblks_x);
      _spinor_assign(*(rec_field + index_a(dT*bt + t, dX*bx + x, dY*by + y, dZ*bz + z)), 
                     *(block_list[block_id].basis[id] + ix));
    }
  }
  return;
}

void add_eo_block_to_global(spinor * const globalfield, spinor * const beven, spinor * const bodd, const int blk) {
  int t, x, y, z;
  int i,it,ix,iy,iz;
  int even = 0, odd = 0;

  for(t = 0; t < block_list[blk].BT; t++) {
    it = t + block_list[blk].mpilocal_coordinate[0]*block_list[blk].BT;
    for(x = 0; x < block_list[blk].BLX; x++) {
      ix = x +  block_list[blk].mpilocal_coordinate[1]*block_list[blk].BLX;
      for(y = 0; y < block_list[blk].BLY; y++) {
        iy = y +  block_list[blk].mpilocal_coordinate[2]*block_list[blk].BLY;
        for(z = 0; z < block_list[blk].BLZ; z++) {
          iz = z +  block_list[blk].mpilocal_coordinate[3]*block_list[blk].BLZ;
          i = g_ipt[it][ix][iy][iz];
          if((t+x+y+z)%2 == 0) {
            add(globalfield + i, globalfield + i, beven + even, 1);
            even++;
          }
          else {
            add(globalfield + i, globalfield + i, bodd + odd, 1);
            odd++;
          }
        }
      }
    }
  }
  return;
}

void add_eo_block_32_to_global(spinor * const globalfield, 
                               spinor32 * const beven, spinor32 * const bodd, const int blk) {
  int i,it,ix,iy,iz;
  int even = 0, odd = 0;
  spinor32 * tmp = NULL;
  _Complex double * to = NULL;
  _Complex float * from = NULL;

  for(int t = 0; t < block_list[blk].BT; t++) {
    it = t + block_list[blk].mpilocal_coordinate[0]*block_list[blk].BT;
    for(int x = 0; x < block_list[blk].BLX; x++) {
      ix = x +  block_list[blk].mpilocal_coordinate[1]*block_list[blk].BLX;
      for(int y = 0; y < block_list[blk].BLY; y++) {
        iy = y +  block_list[blk].mpilocal_coordinate[2]*block_list[blk].BLY;
        for(int z = 0; z < block_list[blk].BLZ; z++) {
          iz = z +  block_list[blk].mpilocal_coordinate[3]*block_list[blk].BLZ;
          i = g_ipt[it][ix][iy][iz];
          if((t+x+y+z)%2 == 0) {
            tmp = beven + even;
            even++;
          }
          else {
            tmp = bodd + odd;
            odd++;
          }
          to = (_Complex double*) (globalfield + i);
          from = (_Complex float*) tmp;
          for(int k = 0; k < 12; k++) {
            (*to) += (_Complex double) (*from);
            to++;
            from++;
          }

        }
      }
    }
  }
  return;
}


void add_block_to_global(spinor * const globalfield, spinor * const blockfield, const int blk) {
  int i;
  spinor * r, * s;
  int it,ix,iy,iz;
  int ibt,ibx,iby,ibz;
  int itb,ixb,iyb,izb;
  int ixcurrent;
  
  ibz = blk%nblks_z;
  iby = (blk / nblks_z)%nblks_y;
  ibx = (blk / (nblks_y * nblks_z))%nblks_x;
  ibt = blk / (nblks_x * nblks_y * nblks_z);
  
  ixcurrent = 0;
  for (i = 0; i < VOLUME; i++) {
   
    iz = i%LZ;
    iy = (i / LZ)%LY;
    ix = (i / (LY * LZ))%LX;
    it = i / (LX * LY * LZ);
    

    izb = iz / block_list[blk].BLZ;
    iyb = iy / block_list[blk].BLY;
    ixb = ix / block_list[blk].BLX;
    itb = it / block_list[blk].BT;

    if ((ibz == izb) && (iby == iyb) && (ibx == ixb) && (ibt == itb)) {
      r = globalfield + i;
      s = blockfield + ixcurrent;
      add(r, r, s, 1);
      ixcurrent++;
    }
  }
  return;
}

/*      eo -> lexic
 *      P: new spinor with full volume 
 *      s: source spinor even 
 *      r: source spinor odd 
 */
void block_convert_eo_to_lexic(spinor * const P, spinor * const s, spinor * const r) {
  int x, y, z, t, i, ix;
  spinor * p = NULL;

  for(x = 0; x < dX; x++) {
    for(y = 0; y < dY; y++) {
      for(z = 0; z < dZ; z++) {
        for(t = 0; t < dT; t++) {
          ix = block_ipt[t][x][y][z];
          i = ix / 2;
          if((x + y + z + t)%2 == 0) {
            p = s;
          }
          else {
            p = r;
          }
          memcpy((P+ix), (p+i), sizeof(spinor));
        }
      }
    }
  }
  return;
}

/*      lexic -> eo
 *      P: source spinor with full volume 
 *      s: new spinor even 
 *      r: new spinor odd 
 */
void block_convert_lexic_to_eo(spinor * const s, spinor * const r, spinor * const P) {
  int x, y, z, t, i, ix;
  spinor * p = NULL;

  for(x = 0; x < dX; x++) {
    for(y = 0; y < dY; y++) {
      for(z = 0; z < dZ; z++) {
        for(t = 0; t < dT; t++) {
          ix = block_ipt[t][x][y][z];
          i = ix / 2;
          if((x + y + z + t)%2 == 0) {
            p = s;
          }
          else {
            p = r;
          }
          memcpy((p+i), (P+ix), sizeof(spinor));
        }
      }
    }
  }
  return;
}

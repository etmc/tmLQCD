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


#ifndef _BLOCK_H
#define _BLOCK_H

#include "su3.h"
#include "su3spinor.h"

_Complex double * little_A;
_Complex float * little_A32;
_Complex double * little_A_eo;
_Complex float * little_A32_eo;


typedef struct {
  /**** Data members ****/
  int volume;                   /* the block local 4 volume */
  int id;                       /* mpilocal block id */
  int BLX, BLY, BLZ, BT;            /* block local sizes */
  int ns;                       /* the number of basis fields, which is needed almost everywhere */
  int coordinate[4];            /* global block coordinate */
  int mpilocal_coordinate[4];   /* mpi process local coordinate */
  int mpilocal_neighbour[8];    /* contains the block id of mpilocal neighbours, or -1 if non-mpilocal */
  int *idx;                     /* provides the next neighbours for spinors on the block */
  int *evenidx;                 /* provides the next neighbours for spinors on the block even/odd case */
  int *oddidx;                 /* provides the next neighbours for spinors on the block even/odd case */
  spinor **basis;               /* generated orthonormal basis for little D [Ns x local_volume] */
  su3 * u;                      /* block local gauge field, for use in D */
  int spinpad;                  /* number of elements needed to store the boundaries of the spinor */
  int evenodd;                  /* block even or odd (0 or 1) */

  /* storage will be g_Ns x (9 * g_Ns)                 */
  /* build_little_diraclocal g_Ns x g_Ns block first (the diagonal part) */
  /* then +t, -t, +x, -x, +y, -y, +z, -z               */
  _Complex double    *little_dirac_operator;  /* full dense representation of the little D */
  _Complex float  *little_dirac_operator32;
  _Complex double    *little_dirac_operator_eo;  /* full dense representation of the little D in e/o order */
} block;

int init_blocks(const int nt, const int nx, const int ny, const int nz);
int free_blocks();

int init_blocks_gaugefield();
int init_blocks_eo_gaugefield();

void copy_global_to_block(spinor * const blockfield, spinor * const globalfield, const int blk);
void copy_block_to_global(spinor * const globalfield, spinor * const blockfield, const int blk);
void copy_global_to_block_eo(spinor * const beven, spinor * const bodd, spinor * const globalfield, const int blk);
void copy_block_eo_to_global(spinor * const globalfield, spinor * const beven, spinor * const bodd, const int blk);
void add_block_to_global(spinor * const globalfield, spinor * const blockfield, const int blk);
void add_eo_block_to_global(spinor * const globalfield, spinor * const beven, spinor * const bodd, const int blk);

void block_convert_lexic_to_eo(spinor * const s, spinor * const r, spinor * const P);
void block_convert_eo_to_lexic(spinor * const P, spinor * const s, spinor * const r);

void block_orthonormalize(block *parent);
void block_orthonormalize_free(block *parent);

void compute_little_D();
void compute_little_D_diagonal();
void alt_block_compute_little_D();

extern int dfl_field_iter;
extern int dfl_poly_iter;

int nb_blocks;
int nblks_t;
int nblks_x;
int nblks_y;
int nblks_z;
int nblks_dir[4];
int blk_gauge_eo;
void reconstruct_global_field_GEN(spinor * const rec_field, spinor ** const psi, int nb_blocks);
void reconstruct_global_field_GEN_ID(spinor * const rec_field, block * const block_list, const int id, const int nb_blocks);
int split_global_field_GEN(spinor ** const psi, spinor * const field, int nb_blocks);
int split_global_field_GEN_ID(block * const block_list, const int id, spinor * const field, const int nb_blocks);

/* Functions for index manipulation related to blocks, C. Tadonki */
int index_a(int t, int x, int y, int z);
int index_b(int t, int x, int y, int z);
int block_index(int t, int x, int y, int z);

extern block * block_list;

#endif

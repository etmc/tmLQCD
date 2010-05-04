/***********************************************************************
 * $Id$ 
 *
 * Copyright (C) 2008 Albert Deuzeman, Siebren Reker, Carsten Urbach
 *               2010 Claude Tadonki
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

complex * little_A;

typedef struct {
  /**** Data members ****/
  int volume;                   /* the block local 4 volume */
  int id;                       /* mpilocal block id */
  int LX, LY, LZ, T;            /* block local sizes */
  int ns;                       /* the number of basis fields, which is needed almost everywhere */
  int coordinate[4];            /* global block coordinate */
  int mpilocal_coordinate[4];   /* mpi process local coordinate */
  int mpilocal_neighbour[8];    /* contains the block id of mpilocal neighbours, or -1 if non-mpilocal */
  int *idx;                     /* provides the next neighbours for spinors on the block */
  spinor **basis;               /* generated orthonormal basis for little D [Ns x local_volume] */
  su3 * u;                      /* block local gauge field, for use in D */
  int spinpad;                  /* number of elements needed to store the boundaries of the spinor */
  int evenodd;                  /* block even or odd (0 or 1) */

  /* storage will be g_Ns x (9 * g_Ns)                 */
  /* build_little_diraclocal g_Ns x g_Ns block first (the diagonal part) */
  /* then +t, -t, +x, -x, +y, -y, +z, -z               */
  complex    *little_dirac_operator;  /* full dense representation of the little D */

} block;

int init_blocks();
int free_blocks();

int split_global_field(spinor * const block_low, spinor * const block_high, spinor * const field);
void reconstruct_global_field(spinor * const rec_field, spinor * const block_low, spinor * const block_high);
void copy_global_to_block(spinor * const blockfield, spinor * const globalfield, const int blk);
void copy_block_to_global(spinor * const globalfield, spinor * const blockfield, const int blk);
void add_block_to_global(spinor * const globalfield, spinor * const blockfield, const int blk);

void block_orthonormalize(block *parent);
void block_orthonormalize_free(block *parent);

void compute_little_D();
void compute_little_D_diagonal();
void alt_block_compute_little_D();

/* CT:
The parameter "DflFieldIter" (default 80) refers to what is called
GSL in Gilbert reports, and "DflPolyIter" (default 20) is the
number of iterations in the polynomial preconditioner. 
*/
extern int dfl_field_iter;
extern int dfl_poly_iter;

extern int dfl_nblock_t;
extern int dfl_nblock_x;
extern int dfl_nblock_y;
extern int dfl_nblock_z;
int nb_blocks;
int nblks_t;
int nblks_x;
int nblks_y;
int nblks_z;
int nblks_dir[4];
void reconstruct_global_field_GEN(spinor * const rec_field, spinor ** const psi, int nb_blocks);
void reconstruct_global_field_GEN_ID(spinor * const rec_field, block * const block_list, int id, int nb_blocks);
int split_global_field_GEN(spinor ** const psi, spinor * const field, int nb_blocks);
int split_global_field_GEN_ID(block * const block_list, int id, spinor * const field, int nb_blocks);

extern block * block_list;

#endif

/* $Id$ */

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

  /* storage will be g_Ns x (9 * g_Ns)                 */
  /* build_little_diraclocal g_Ns x g_Ns block first (the diagonal part) */
  /* then +t, -t, +x, -x, +y, -y, +z, -z               */
  complex    *little_dirac_operator;  /* full dense representation of the little D */

} block;

int init_blocks();
int free_blocks();

int split_global_field(spinor * const block_low, spinor * const block_high, spinor * const field);
void reconstruct_global_field(spinor * const rec_field, spinor * const block_low, spinor * const block_high);

void block_orthonormalize(block *parent);
void block_orthonormalize_free(block *parent);

void block_compute_little_D_offdiagonal();
void block_compute_little_D_diagonal();
void alt_block_compute_little_D();

extern block * block_list;

#endif

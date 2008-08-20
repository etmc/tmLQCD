#ifndef _BLOCK_H
#define _BLOCK_H

#include "su3.h"
#include "su3spinor.h"

typedef struct {
  /**** Data members ****/
  int volume;                   /* the block local 4 volume */
  int LX, LY, LZ, T;            /* block local sizes */
  int ns;                       /* the number of basis fields, which is needed almost everywhere */
  int coordinate[4];            /* global block coordinate */
  int mpilocal_coordinate[4];   /* mpi process local coordinate */
  int mpilocal_neighbour[8];    /* contains the block id of mpilocal neighbours, or -1 if non-mpilocal */
  int *idx;                     /* provides the next neighbours for spinors on the block */
  spinor *basis;                /* generated orthonormal basis for little D [Ns x local_volume] */
  spinor **neighbour_edges;     /* boundary terms of the basis of the neighbours [8 x Ns x surface_term] */
  su3 * u;                      /* block local gauge field, for use in D */
  int spinpad;                  /* number of elements needed to store the boundaries of the spinor */

  /* storage will be g_Ns x (9 * g_Ns)                 */
  /* build_little_diraclocal g_Ns x g_Ns block first (the diagonal part) */
  /* then +t, -t, +x, -x, +y, -y, +z, -z               */
  complex    *little_dirac_operator;  /* full dense representation of the little D */

} block;

int init_blocks();
int free_blocks();

int add_basis_field(int const index, spinor const *field);
int init_gauge_blocks(su3 ** const field);
int init_geom_blocks();

void block_orthonormalize(block *parent);
void block_reconstruct_global_field(const int index, spinor * const reconstructed_field);

void block_compute_little_D_offdiagonal(block *parent);
void block_compute_little_D_diagonal(block *parent);

extern block * block_list;

#endif

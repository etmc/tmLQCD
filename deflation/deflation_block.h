#ifndef _DEFLATION_BLOCK_H
#define _DEFLATION_BLOCK_H


#include "../su3.h"
#include "../su3spinor.h"

extern int LITTLE_BASIS_SIZE;

typedef struct {
  /**** Data members ****/
  int local_volume;                   /* the block local 4 volume */
  int little_basis_size;              /* the value of Ns, which is needed almost everywhere */
  int coordinate[4];                  /* global block coordinate */
  int mpilocal_coordinate[4];         /* mpi process local coordinate */
  spinor *little_basis;               /* generated orthonormal basis for little D [Ns x local_volume] */
  spinor **little_neighbour_edges;    /* boundary terms of the basis of the neighbours [8 x Ns x surface_term] */
  su3 * u;                            /* block local gauge field, for use in D */

  /* storage will be g_Ns x (9 * g_Ns)                 */
  /* build_little_diraclocal g_Ns x g_Ns block first (the diagonal part) */
  /* then +t, -t, +x, -x, +y, -y, +z, -z               */
  complex    *little_dirac_operator;  /* full dense representation of the litle D */
  complex    *local_little_field;     /* memory reserved for a spinor in little space [Ns] */

  int * iup, * idn;                   /* provides the next neigbours on the block */

  /**** 'Member' functions ****/
  void (*orthonormalize)(void *parent);
  spinor *(*reconstruct_global_field)(void *parent, int const index, spinor *reconstructed_field);
  void (*)(void *parent);
} deflation_block;

int init_deflation_blocks();
int free_deflation_blocks();

int add_basis_field(int const index, spinor const *field);
int copy_block_gauge(su3 const *field);

void block_orthonormalize(void *parent);
spinor *block_reconstruct_global_field(void *parent, int const index, spinor *reconstructed_field);
void block_build_little_dirac(void *parent);

#endif

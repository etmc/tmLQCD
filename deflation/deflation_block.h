#ifndef _DEFLATION_BLOCK_H
#define _DEFLATION_BLOCK_H


#include "su3.h"
#include "../su3spinor.h"

extern int LITTLE_BASIS_SIZE;

typedef struct
{
  /* Data members */
  /* Include, possibly, E/O boolean, storage for 'big' fields */
  /* the block local 4 volume */
  int volume;
  /* global block coordinate */
  int coordinate[4];
  /* mpi process local coordintae */
  int mpilocal_coordinate[4];
  spinor **little_basis;
  spinor **little_neighbour_edges;
  /* block local gauge field */
  su3 * u;

  /* storage will be g_Ns x 9*g_Ns                     */
  /* local g_Ns x g_Ns block first (the diagonal part) */
  /* then +t, -t, +x, -x, +y, -y, +z, -z               */
  complex    *little_dirac_operator;
  complex    *local_little_field;
  
  /* define the next neigbours on the block */
  int * iup, * idn;

  /* 'Member' functions */
  void (*orthonormalize)(void);
  void (*build_little_dirac)(void);
  
} deflation_block;

int init_deflation_blocks();
int free_deflation_blocks();

int add_basis_field(int const index, spinor const *field);

#endif

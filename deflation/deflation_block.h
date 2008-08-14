#ifndef _DEFLATION_BLOCK_H
#define _DEFLATION_BLOCK_H

#include "../su3spinor.h"

extern int LITTLE_BASIS_SIZE;

typedef struct
{
  /* Data members */
  /* Include, possibly, E/O boolean, storage for 'big' fields */
  su3spinor **little_basis;
  su3spinor **little_neighbour_edges;
  complex    *little_dirac_operator;

  complex    *local_little_field;
  
  /* 'Member' functions */
  void (*orthonormalize)(void);
  void (*build_little_dirac)(void);
  
} deflation_block;

int init_deflation_blocks(int const nblocks);
int free_deflation_blocks();

void add_basis_field(int const index, su3spinor const *field);

#endif

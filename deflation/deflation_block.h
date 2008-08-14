#ifndef _DEFLATION_BLOCK_H
#define _DEFLATION_BLOCK_H

#include "../su3spinor.h"

extern int LITTLE_BASIS_SIZE;

typedef struct
{
  /* Data members */
  /* Include, possibly, E/O boolean, storage for 'big' fields */
  spinor **little_basis;
  spinor **little_neighbour_edges;
  complex    *little_dirac_operator;

  complex    *local_little_field;
  
  /* 'Member' functions */
  void (*orthonormalize)(void);
  void (*build_little_dirac)(void);
  
} deflation_block;

int init_deflation_blocks();
int free_deflation_blocks();

int add_basis_field(int const index, spinor const *field);

#endif

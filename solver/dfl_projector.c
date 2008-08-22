#include "dfl_projector.h"
#include "../blocks.h"

/* Break up full volume spinor to blocks
 * loop over block.basis
 * compute inner product and store as complex vector
 * compute A^-1 * complex vector
 * loop over block.basis
 * compute sum of basis vectors times complex element
 * create global vector */

void project(spinor * const out, spinor * const in) {
  int i,j,ctr_t;
  spinor **psi;
  int contig_block = LZ / 2;
  complex *inprod;
  complex *invvec;

  psi = calloc(2, sizeof(spinor*)); /*block local version of global spinor */
  inprod = calloc(2 * 9 * g_N_s, sizeof(complex)); /*inner product of spinors with bases */
  invvec = calloc(2 * 9 * g_N_s, sizeof(complex)); /*inner product of spinors with bases */
  
  /* no loop below because further down we also don't take this cleanly into account */
  psi[0] = calloc(VOLUME / 2, sizeof(spinor));
  psi[1] = calloc(VOLUME / 2, sizeof(spinor));

  /*initialize the local (block) parts of the spinor*/
  for (ctr_t = 0; ctr_t < (VOLUME / LZ); ++ctr_t)
  {
    memcpy(psi[0], in + (2 * ctr_t) * contig_block, contig_block * sizeof(spinor));
    memcpy(psi[1], in + (2 * ctr_t + 1) * contig_block, contig_block * sizeof(spinor));
  }
  for (i = 0; i < 2; ++i) {/* loop over blocks */
    /* compute inner product */
    for (j = 0; j < g_N_s; ++j) {/*loop over block.basis */
      inprod[j + i * 2 * 9 * g_N_s] = block_scalar_prod(block_list[i].basis[j], psi[i], block_list[i].volume);
    }
  }
  lgcr(invvec, inprod, 10, 1000, 1.e-15, 1, 2 * g_N_s, 2 * 9 * g_N_s, &little_D);

}
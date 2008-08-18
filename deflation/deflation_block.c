#include "../global.h"
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "deflation_block.h"
#include "../linalg/diff.h"

#define CALLOC_ERROR_CRASH {printf ("calloc errno : %d\n", errno); errno = 0; return 1;}

extern deflation_block *g_deflation_blocks;
extern int LX,LY,LZ,T,VOLUME;
extern int LITTLE_BASIS_SIZE;

int init_deflation_blocks()
{
  int i;

  g_deflation_blocks = calloc(2, sizeof(deflation_block));
  for (i = 0; i < 2; ++i) {
    if ((void*)(g_deflation_blocks[i].little_basis = calloc(LITTLE_BASIS_SIZE, sizeof(spinor *))) == NULL)
      CALLOC_ERROR_CRASH;

    if ((void*)(g_deflation_blocks[i].little_basis = calloc(LITTLE_BASIS_SIZE * VOLUME, sizeof(spinor))) == NULL)
      CALLOC_ERROR_CRASH;

    if ((void*)(g_deflation_blocks[i].little_neighbour_edges = calloc(8, sizeof(spinor *))) == NULL)
      CALLOC_ERROR_CRASH;

    if ((void*)(g_deflation_blocks[i].little_neighbour_edges[0] = calloc(VOLUME / (2 * LX), sizeof(spinor))) == NULL)
      CALLOC_ERROR_CRASH;
    if ((void*)(g_deflation_blocks[i].little_neighbour_edges[1] = calloc(VOLUME / (2 * LX), sizeof(spinor))) == NULL)
      CALLOC_ERROR_CRASH;
    if ((void*)(g_deflation_blocks[i].little_neighbour_edges[2] = calloc(VOLUME / (2 * LY), sizeof(spinor))) == NULL)
      CALLOC_ERROR_CRASH;
    if ((void*)(g_deflation_blocks[i].little_neighbour_edges[3] = calloc(VOLUME / (2 * LY), sizeof(spinor))) == NULL)
      CALLOC_ERROR_CRASH;
    if ((void*)(g_deflation_blocks[i].little_neighbour_edges[4] = calloc(VOLUME / LZ, sizeof(spinor))) == NULL)
      CALLOC_ERROR_CRASH;
    if ((void*)(g_deflation_blocks[i].little_neighbour_edges[5] = calloc(VOLUME / LZ, sizeof(spinor))) == NULL)
      CALLOC_ERROR_CRASH;
    if ((void*)(g_deflation_blocks[i].little_neighbour_edges[6] = calloc(VOLUME / (2 * T), sizeof(spinor))) == NULL)
      CALLOC_ERROR_CRASH;
    if ((void*)(g_deflation_blocks[i].little_neighbour_edges[7] = calloc(VOLUME / (2 * T), sizeof(spinor))) == NULL)
      CALLOC_ERROR_CRASH;

    if ((void*)(g_deflation_blocks[i].little_dirac_operator = calloc(9 * LITTLE_BASIS_SIZE * LITTLE_BASIS_SIZE, sizeof(complex))) == NULL)
      CALLOC_ERROR_CRASH;

    if ((void*)(g_deflation_blocks[i].local_little_field = calloc(LITTLE_BASIS_SIZE, sizeof(complex))) == NULL)
      CALLOC_ERROR_CRASH;

    g_deflation_blocks[i].orthonormalize = &block_orthonormalize;
    g_deflation_blocks[i].reconstruct_global_field = &block_reconstruct_global_field;

  }
  return 0;
}

int free_deflation_blocks()
{
  int i, j;

  for(i = 0; i < 2; ++i){
    free(g_deflation_blocks[i].little_basis);

    for (j = 0; j < 8; ++j)
      free(g_deflation_blocks[i].little_neighbour_edges[j]);
    free(g_deflation_blocks[i].little_neighbour_edges);

    free(g_deflation_blocks[i].little_dirac_operator);
    free(g_deflation_blocks[i].local_little_field);
  }

  free(g_deflation_blocks);
  return 0;
}

int add_basis_field(int const index, spinor const *field)
{
  int ctr_t;
  int contig_block = LZ / 2;
  for (ctr_t = 0; ctr_t < (2 * VOLUME / LZ); ++ctr_t)
  {
    memcpy(g_deflation_blocks[0].little_basis + index * VOLUME, field + (2 * ctr_t) * contig_block, contig_block * sizeof(spinor));
    memcpy(g_deflation_blocks[1].little_basis + index * VOLUME, field + (2 * ctr_t + 1) * contig_block, contig_block * sizeof(spinor));
  }
  return 0;
}

int copy_block_gauge(su3 const *field)
{
  /* Pseudocode:
  The purpose of this function is to store locally all forward and backward gauge links in a block.
  All positive directions should be available on this processor already, so they can be copied in immediately.
  A subset of the backward links can be done directly, so do them next.
  Then in every direction a set of links needs to be copied from a neighbouring processor and plugged in in the proper location.
  Hermitian conjugation for the backwards part is handled in the Dirac operator explicitly through the _su3_inverse_multiply macro.
  */
}

/* the following should be somewhere else ... */

complex block_scalar_prod_Ugamma(spinor * const r, spinor * const s, 
				const int mu, const int N) {
  complex c;

  return(c);
}

complex block_scalar_prod(spinor * const R, spinor * const S, const int N) {
  int ix;
  static double ks,kc,ds,tr,ts,tt;
  complex c;
  spinor * r;
  spinor * s;
  /* Real Part */

  ks=0.0;
  kc=0.0;
#if (defined BGL && defined XLC)
  __alignx(16, S);
  __alignx(16, R);
#endif  
  for (ix = 0; ix < N; ix++){
    s=(spinor *) S + ix;
    r=(spinor *) R + ix;
    
    ds=(*r).s0.c0.re*(*s).s0.c0.re+(*r).s0.c0.im*(*s).s0.c0.im+
       (*r).s0.c1.re*(*s).s0.c1.re+(*r).s0.c1.im*(*s).s0.c1.im+
       (*r).s0.c2.re*(*s).s0.c2.re+(*r).s0.c2.im*(*s).s0.c2.im+
       (*r).s1.c0.re*(*s).s1.c0.re+(*r).s1.c0.im*(*s).s1.c0.im+
       (*r).s1.c1.re*(*s).s1.c1.re+(*r).s1.c1.im*(*s).s1.c1.im+
       (*r).s1.c2.re*(*s).s1.c2.re+(*r).s1.c2.im*(*s).s1.c2.im+
       (*r).s2.c0.re*(*s).s2.c0.re+(*r).s2.c0.im*(*s).s2.c0.im+
       (*r).s2.c1.re*(*s).s2.c1.re+(*r).s2.c1.im*(*s).s2.c1.im+
       (*r).s2.c2.re*(*s).s2.c2.re+(*r).s2.c2.im*(*s).s2.c2.im+
       (*r).s3.c0.re*(*s).s3.c0.re+(*r).s3.c0.im*(*s).s3.c0.im+
       (*r).s3.c1.re*(*s).s3.c1.re+(*r).s3.c1.im*(*s).s3.c1.im+
       (*r).s3.c2.re*(*s).s3.c2.re+(*r).s3.c2.im*(*s).s3.c2.im;

    /* Kahan Summation */    
    tr=ds+kc;
    ts=tr+ks;
    tt=ts-ks;
    ks=ts;
    kc=tr-tt;
  }
  c.re = ks+kc;

  /* Imaginary Part */

  ks=0.0;
  kc=0.0;
  
  for (ix=0;ix<N;ix++){
    s=(spinor *) S + ix;
    r=(spinor *) R + ix;
    
    ds=-(*r).s0.c0.re*(*s).s0.c0.im+(*r).s0.c0.im*(*s).s0.c0.re-
      (*r).s0.c1.re*(*s).s0.c1.im+(*r).s0.c1.im*(*s).s0.c1.re-
      (*r).s0.c2.re*(*s).s0.c2.im+(*r).s0.c2.im*(*s).s0.c2.re-
      (*r).s1.c0.re*(*s).s1.c0.im+(*r).s1.c0.im*(*s).s1.c0.re-
      (*r).s1.c1.re*(*s).s1.c1.im+(*r).s1.c1.im*(*s).s1.c1.re-
      (*r).s1.c2.re*(*s).s1.c2.im+(*r).s1.c2.im*(*s).s1.c2.re-
      (*r).s2.c0.re*(*s).s2.c0.im+(*r).s2.c0.im*(*s).s2.c0.re-
      (*r).s2.c1.re*(*s).s2.c1.im+(*r).s2.c1.im*(*s).s2.c1.re-
      (*r).s2.c2.re*(*s).s2.c2.im+(*r).s2.c2.im*(*s).s2.c2.re-
      (*r).s3.c0.re*(*s).s3.c0.im+(*r).s3.c0.im*(*s).s3.c0.re-
      (*r).s3.c1.re*(*s).s3.c1.im+(*r).s3.c1.im*(*s).s3.c1.re-
      (*r).s3.c2.re*(*s).s3.c2.im+(*r).s3.c2.im*(*s).s3.c2.re;
    
    tr=ds+kc;
    ts=tr+ks;
    tt=ts-ks;
    ks=ts;
    kc=tr-tt;
  }
  c.im = ks+kc;
  return(c);
}

double block_two_norm(spinor * const R, const int N) {
  int ix;
  static double ks,kc,ds,tr,ts,tt;
  double norm;
  spinor * r;
  /* Real Part */

  ks=0.0;
  kc=0.0;
#if (defined BGL && defined XLC)
  __alignx(16, S);
  __alignx(16, R);
#endif  
  for (ix = 0; ix < N; ix++){
    r=(spinor *) R + ix;
    
    ds=(*r).s0.c0.re*(*r).s0.c0.re+(*r).s0.c0.im*(*r).s0.c0.im+
       (*r).s0.c1.re*(*r).s0.c1.re+(*r).s0.c1.im*(*r).s0.c1.im+
       (*r).s0.c2.re*(*r).s0.c2.re+(*r).s0.c2.im*(*r).s0.c2.im+
       (*r).s1.c0.re*(*r).s1.c0.re+(*r).s1.c0.im*(*r).s1.c0.im+
       (*r).s1.c1.re*(*r).s1.c1.re+(*r).s1.c1.im*(*r).s1.c1.im+
       (*r).s1.c2.re*(*r).s1.c2.re+(*r).s1.c2.im*(*r).s1.c2.im+
       (*r).s2.c0.re*(*r).s2.c0.re+(*r).s2.c0.im*(*r).s2.c0.im+
       (*r).s2.c1.re*(*r).s2.c1.re+(*r).s2.c1.im*(*r).s2.c1.im+
       (*r).s2.c2.re*(*r).s2.c2.re+(*r).s2.c2.im*(*r).s2.c2.im+
       (*r).s3.c0.re*(*r).s3.c0.re+(*r).s3.c0.im*(*r).s3.c0.im+
       (*r).s3.c1.re*(*r).s3.c1.re+(*r).s3.c1.im*(*r).s3.c1.im+
       (*r).s3.c2.re*(*r).s3.c2.re+(*r).s3.c2.im*(*r).s3.c2.im;

    /* Kahan Summation */    
    tr=ds+kc;
    ts=tr+ks;
    tt=ts-ks;
    ks=ts;
    kc=tr-tt;
  }
  norm = ks+kc;
  return(norm);
}

void compute_little_D_diagonal(void *parent) {
  int i,j;
  /* we need working space, where do we best allocate it? */
  spinor * tmp;
  deflation_block *this = (deflation_block*)parent;
  complex * M = this->little_dirac_operator;

  
  for(i = 0; i < this->little_basis_size; i++) {
    Block_D_psi(tmp, this->little_basis[i]);
    for(j = 0; j < this->little_basis_size; j++) {
      /* order correct ? */
      M[i*this->little_basis_size + j] = block_scalar_prod(this->little_basis + j * this->local_volume, tmp, this->local_volume);
    }
  }
  return;
}

void block_compute_little_D_offdiagonal(void *parent) {
/*   Here we need to multiply the boundary with the corresponding  */
/*   U and gamma_i and take the scalar product then */
/*   NOTE I will assume the boundaries are available, for the time being */
  deflation_block *this = (deflation_block*)parent;
  int i,j;
  spinor * tmp; 
  /* Start by going for the first block (t-up) */
  complex * M = this->little_dirac_operator + this->little_basis_size * this->little_basis_size;
  spinor *neighbour = this->little_neighbour_edges[0];
  for(i = 0; i < this->little_basis_size; i++) {
    boundary_D(tmp, this->little_basis[i]); /* NOTE Syntax! */
    for(j = 0; j < this->little_basis_size; j++) {
      /* order correct ? */
      M[i*this->little_basis_size + j] = block_scalar_prod(this->little_basis + j * this->local_volume, tmp, this->local_volume / T);
    }
  }
  /* Start by going for the first block (t-down) */
  M += this->little_basis_size * this->little_basis_size;
  neighbour = this->little_neighbour_edges[1];
  for(i = 0; i < this->little_basis_size; i++) {
    boundary_D(tmp, this->little_basis[i]); /* NOTE Syntax! */
    for(j = 0; j < this->little_basis_size; j++) {
      /* order correct ? */
      M[i*this->little_basis_size + j] = block_scalar_prod(this->little_basis + j * this->local_volume, tmp, this->local_volume / T);
    }
  }
  /* Start by going for the first block (x-up) */
  M += this->little_basis_size * this->little_basis_size;
  neighbour = this->little_neighbour_edges[2];
  for(i = 0; i < this->little_basis_size; i++) {
    boundary_D(tmp, this->little_basis[i]); /* NOTE Syntax! */
    for(j = 0; j < this->little_basis_size; j++) {
      /* order correct ? */
      M[i*this->little_basis_size + j] = block_scalar_prod(this->little_basis + j * this->local_volume, tmp, this->local_volume / LX);
    }
  }
  /* Start by going for the first block (x-down) */
  M += this->little_basis_size * this->little_basis_size;
  neighbour = this->little_neighbour_edges[3];
  for(i = 0; i < this->little_basis_size; i++) {
    boundary_D(tmp, this->little_basis[i]); /* NOTE Syntax! */
    for(j = 0; j < this->little_basis_size; j++) {
      /* order correct ? */
      M[i*this->little_basis_size + j] = block_scalar_prod(this->little_basis + j * this->local_volume, tmp, this->local_volume / LX);
    }
  }
  /* Start by going for the first block (y-up) */
  M += this->little_basis_size * this->little_basis_size;
  neighbour = this->little_neighbour_edges[4];
  for(i = 0; i < this->little_basis_size; i++) {
    boundary_D(tmp, this->little_basis[i]); /* NOTE Syntax! */
    for(j = 0; j < this->little_basis_size; j++) {
      /* order correct ? */
      M[i*this->little_basis_size + j] = block_scalar_prod(this->little_basis + j * this->local_volume, tmp, this->local_volume / LY);
    }
  }
  /* Start by going for the first block (y-down) */
  M += this->little_basis_size * this->little_basis_size;
  neighbour = this->little_neighbour_edges[5];
  for(i = 0; i < this->little_basis_size; i++) {
    boundary_D(tmp, this->little_basis[i]); /* NOTE Syntax! */
    for(j = 0; j < this->little_basis_size; j++) {
      /* order correct ? */
      M[i*this->little_basis_size + j] = block_scalar_prod(this->little_basis + j * this->local_volume, tmp, this->local_volume / LY);
    }
  }
  /* Start by going for the first block (z-up) */
  M += this->little_basis_size * this->little_basis_size;
  neighbour = this->little_neighbour_edges[6];
  for(i = 0; i < this->little_basis_size; i++) {
    boundary_D(tmp, this->little_basis[i]); /* NOTE Syntax! */
    for(j = 0; j < this->little_basis_size; j++) {
      /* order correct ? */
      M[i*this->little_basis_size + j] = block_scalar_prod(this->little_basis + j * this->local_volume, tmp, this->local_volume / LZ);
    }
  }
  /* Start by going for the first block (z-down) */
  M += this->little_basis_size * this->little_basis_size;
  neighbour = this->little_neighbour_edges[7];
  for(i = 0; i < this->little_basis_size; i++) {
    boundary_D(tmp, this->little_basis[i]); /* NOTE Syntax! */
    for(j = 0; j < this->little_basis_size; j++) {
      /* order correct ? */
      M[i*this->little_basis_size + j] = block_scalar_prod(this->little_basis + j * this->local_volume, tmp, this->local_volume / LZ);
    }
  }
}

/* Uses a Modified Gram-Schmidt algorithm to orthonormalize little basis vectors */
void block_orthonormalize(void *parent){
  int i, j, k;
  spinor *current, *next, *iter;
  spinor orig, resbasis;
  complex coeff;
  complex scale;
  deflation_block *this = (deflation_block*)parent;

  scale.im = 0;
  for(i = 0; i < this->little_basis_size; ++i){
    /* rescale the current vector */
    current = this->little_basis + i * this->local_volume;
    scale.re = 1 / sqrt(block_two_norm(current, this->local_volume));
    for(iter = current; iter < current + this->local_volume; ++iter){
      orig = *iter; /* we can't alias */
      _spinor_mul_complex(*iter, scale, orig);
    }
    /* rescaling done, now subtract this direction from all vectors that follow */
    for(j = i + 1; j < this->little_basis_size; ++j){
      next = this->little_basis + j * this->local_volume;
      coeff = block_scalar_prod(current, next, this->local_volume);
      for(k = 0; k < this->local_volume; ++k){
        _spinor_mul_complex(resbasis, coeff, *(current + k));
        orig = *(next + k);
        diff(next + k, &orig, &resbasis, 1); /* NOTE We could also subtract the complete spinors from eachother */
      }
    }
  }
}

/* Reconstructs a global field from the little basis of two blocks */
spinor *block_reconstruct_global_basis(void *parent, int const index, spinor *reconstructed_field)
{
  int ctr_t;
  int contig_block = LZ / 2;
  for (ctr_t = 0; ctr_t < (2 * VOLUME / LZ); ++ctr_t)
  {
    memcpy(reconstructed_field + (2 * ctr_t) * contig_block, g_deflation_blocks[0].little_basis + index * VOLUME,contig_block * sizeof(spinor));
    memcpy(reconstructed_field + (2 * ctr_t + 1) * contig_block, g_deflation_blocks[1].little_basis + index * VOLUME, contig_block * sizeof(spinor));
  }
  return reconstructed_field;
}

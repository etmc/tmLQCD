#include "../global.h"
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "deflation_block.h"
#include "../linalg/diff.h"

#define CALLOC_ERROR_CRASH {printf ("calloc errno : %d\n", errno); errno = 0; return 1;}

extern block *g_blocks;
extern int LX,LY,LZ,T,VOLUME;
extern int g_proc_coords[4]; /* NOTE may need some ifdef guard for MPI availability, though I hope our deflation code will never be made serially capable */
extern int ****g_ipt;        /* Used in the block gauge field initialization */
extern int **g_iup, **g_idn; /* Used in the block gauge field initialization */

int init_blocks()
{
  int i,j;
  int g_N_s = 20;  /* NOTE hardcoded by hand here until we come up with an input way of defining it */
  g_blocks = calloc(2, sizeof(block));
  for (i = 0; i < 2; ++i) {
    g_blocks[i].volume = VOLUME/2;
    g_blocks[i].LX = LX;
    g_blocks[i].LY = LY;
    g_blocks[i].LZ = LZ/2;
    g_blocks[i].T = T;
    g_blocks[i].ns = g_N_s;
    g_blocks[i].spinpad = 1;
    memcpy(g_blocks[i].mpilocal_coordinate, g_proc_coords, 4);
    memcpy(g_blocks[i].coordinate, g_proc_coords, 3);
    g_blocks[i].coordinate[3] = 2 * g_proc_coords[3] + i;

    if ((void*)(g_blocks[i].idx = calloc(8 * VOLUME/2, sizeof(int))) == NULL)
      CALLOC_ERROR_CRASH;

    if ((void*)(g_blocks[i].basis = calloc(g_N_s * (VOLUME/2 + g_blocks[i].spinpad), sizeof(spinor))) == NULL)
      CALLOC_ERROR_CRASH; /* block volume is half of processor volume, add one spinor for the zero element */
    for (j = 0; j < g_N_s; ++j) /* write a zero element at the end of every spinor */
      _spinor_null(g_blocks[i].basis[j * (VOLUME/2 + g_blocks[i].spinpad)]);

    if ((void*)(g_blocks[i].neighbour_edges = calloc(8, sizeof(spinor *))) == NULL)
      CALLOC_ERROR_CRASH;

    if ((void*)(g_blocks[i].neighbour_edges[0] = calloc(VOLUME / (2 * LX), sizeof(spinor))) == NULL)
      CALLOC_ERROR_CRASH;
    if ((void*)(g_blocks[i].neighbour_edges[1] = calloc(VOLUME / (2 * LX), sizeof(spinor))) == NULL)
      CALLOC_ERROR_CRASH;
    if ((void*)(g_blocks[i].neighbour_edges[2] = calloc(VOLUME / (2 * LY), sizeof(spinor))) == NULL)
      CALLOC_ERROR_CRASH;
    if ((void*)(g_blocks[i].neighbour_edges[3] = calloc(VOLUME / (2 * LY), sizeof(spinor))) == NULL)
      CALLOC_ERROR_CRASH;
    if ((void*)(g_blocks[i].neighbour_edges[4] = calloc(VOLUME / LZ, sizeof(spinor))) == NULL)
      CALLOC_ERROR_CRASH;
    if ((void*)(g_blocks[i].neighbour_edges[5] = calloc(VOLUME / LZ, sizeof(spinor))) == NULL)
      CALLOC_ERROR_CRASH;
    if ((void*)(g_blocks[i].neighbour_edges[6] = calloc(VOLUME / (2 * T), sizeof(spinor))) == NULL)
      CALLOC_ERROR_CRASH;
    if ((void*)(g_blocks[i].neighbour_edges[7] = calloc(VOLUME / (2 * T), sizeof(spinor))) == NULL)
      CALLOC_ERROR_CRASH;

    if ((void*)(g_blocks[i].little_dirac_operator = calloc(9 * g_N_s * g_N_s, sizeof(complex))) == NULL)
      CALLOC_ERROR_CRASH;

    if ((void*)(g_blocks[i].u = calloc(8 * VOLUME, sizeof(su3))) == NULL)
      CALLOC_ERROR_CRASH;
  }
  return 0;
}

int free_blocks()
{
  int i, j;

  for(i = 0; i < 2; ++i) {
    free(g_blocks[i].basis);
    free(g_blocks[i].mpilocal_coordinate);
    free(g_blocks[i].coordinate);

    for (j = 0; j < 8; ++j)
      free(g_blocks[i].neighbour_edges[j]);
    free(g_blocks[i].neighbour_edges);

    free(g_blocks[i].little_dirac_operator);
  }

  free(g_blocks);
  return 0;
}

int add_basis_field(int const index, spinor const *field)
{
  int ctr_t;
  int contig_block = LZ / 2;
  for (ctr_t = 0; ctr_t < (2 * VOLUME / LZ); ++ctr_t)
  {
    memcpy(g_blocks[0].basis + index * (VOLUME + g_blocks[0].spinpad), field + (2 * ctr_t) * contig_block, contig_block * sizeof(spinor));
    memcpy(g_blocks[1].basis + index * (VOLUME + g_blocks[1].spinpad), field + (2 * ctr_t + 1) * contig_block, contig_block * sizeof(spinor));
  }
  return 0;
}

int init_gauge_blocks(su3 const *field)
{
  /* Copies the existing gauge field on the processor into the two separate blocks in a form
  that is readable by the by the block Dirac operator. Specifically, in consecutive memory
  now +t,-t,+x,-x,+y,-y,+z,-z gauge links are stored. This requires double the storage in
  memory. */

  int x, y, z, t, ix, ix_new;
  for (t = 0; t < T; ++t) {
    for (x = 0; x < LX; ++x) {
      for (y = 0; y < LY; ++y) {
        for (z = 0; z < LZ/2; ++z) {
          ix = g_ipt[t][x][y][z];
          ix_new = 8 * (z + y * LZ/2 + x * LY * LZ/2 + t * LX * LY * LZ/2); //su3 index on this block
          memcpy(g_blocks[0].u + ix_new, field + g_iup[ix][0], sizeof(su3));
          memcpy(g_blocks[0].u + ix_new + 1, field + g_idn[ix][0], sizeof(su3));
          memcpy(g_blocks[0].u + ix_new + 2, field + g_iup[ix][1], sizeof(su3));
          memcpy(g_blocks[0].u + ix_new + 3, field + g_idn[ix][1], sizeof(su3));
          memcpy(g_blocks[0].u + ix_new + 4, field + g_iup[ix][2], sizeof(su3));
          memcpy(g_blocks[0].u + ix_new + 5, field + g_idn[ix][2], sizeof(su3));
          memcpy(g_blocks[0].u + ix_new + 6, field + g_iup[ix][3], sizeof(su3));
          memcpy(g_blocks[0].u + ix_new + 7, field + g_idn[ix][3], sizeof(su3));
        }
        for (z = LZ/2; z < LZ; ++z) {
          ix = g_ipt[t][x][y][z];
          ix_new = 8 * (z - LZ/2 + y * LZ/2 + x * LY * LZ/2 + t * LX * LY * LZ/2); //su3 index on this block, count anew in the z direction
          memcpy(g_blocks[1].u + ix_new, field + g_iup[ix][0], sizeof(su3));
          memcpy(g_blocks[1].u + ix_new + 1, field + g_idn[ix][0], sizeof(su3));
          memcpy(g_blocks[1].u + ix_new + 2, field + g_iup[ix][1], sizeof(su3));
          memcpy(g_blocks[1].u + ix_new + 3, field + g_idn[ix][1], sizeof(su3));
          memcpy(g_blocks[1].u + ix_new + 4, field + g_iup[ix][2], sizeof(su3));
          memcpy(g_blocks[1].u + ix_new + 5, field + g_idn[ix][2], sizeof(su3));
          memcpy(g_blocks[1].u + ix_new + 6, field + g_iup[ix][3], sizeof(su3));
          memcpy(g_blocks[1].u + ix_new + 7, field + g_idn[ix][3], sizeof(su3));
        }
      }
    }
  }
  return 0;
}

int init_geom_blocks()
{
  int ix;
  int zstride = 1;
  int ystride = LZ/2;
  int xstride = LY * LZ/2;
  int tstride = LX * LY * LZ/2;
  int boundidx = 4 * VOLUME;
  for (ix = 0; ix < VOLUME/2; ++ix) {
    g_blocks[0].idx[8 * ix + 0] = (ix           >= VOLUME - tstride ? boundidx : ix + tstride);//+t
    g_blocks[0].idx[8 * ix + 1] = (ix           <  tstride          ? boundidx : ix - tstride);//-t
    g_blocks[0].idx[8 * ix + 2] = (ix % tstride >= LZ/2*LY*(LX-1)   ? boundidx : ix + xstride);//+x
    g_blocks[0].idx[8 * ix + 3] = (ix % tstride <  LZ/2*LY          ? boundidx : ix - xstride);//-x
    g_blocks[0].idx[8 * ix + 4] = (ix % xstride >= LZ/2 * (LY - 1)  ? boundidx : ix + ystride);//+y
    g_blocks[0].idx[8 * ix + 5] = (ix % xstride <  LZ/2             ? boundidx : ix - ystride);//-y
    g_blocks[0].idx[8 * ix + 6] = (ix % ystride == LZ/2 - 1         ? boundidx : ix + zstride);//+z
    g_blocks[0].idx[8 * ix + 7] = (ix % ystride == 0                ? boundidx : ix - zstride);//-z
  }
  memcpy(g_blocks[1].idx, g_blocks[0].idx, 8 * VOLUME/2);
  return 0;
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
  spinor *r, *s;
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
  spinor *r;
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
  block *this = (block*)parent;
  complex * M = this->little_dirac_operator;

  
  for(i = 0; i < this->ns; i++){
    Block_D_psi(tmp, this->basis + i * (this->volume + this->spinpad));
    for(j = 0; j < this->ns; j++){
      /* order correct ? */
      M[i*this->ns + j] = block_scalar_prod(this->basis + j * (this->volume + this->spinpad), tmp, this->volume);
    }
  }
  return;
}

void block_compute_little_D_offdiagonal(void *parent) {
/*   Here we need to multiply the boundary with the corresponding  */
/*   U and gamma_i and take the scalar product then */
/*   NOTE I will assume the boundaries are available, for the time being */
  block *this = (block*)parent;
  int i, j, vec_ctr, dir, surface;
  int start, stride;
  spinor tmp;
  spinor **neighbours = this->neighbour_edges;
  spinor *basis;
  complex *M = this->little_dirac_operator + this->ns * this->ns;
  complex result, aggregate[this->ns];
  su3 *gauge;

  /* +T direction (+0) */
  for(vec_ctr = 0; vec_ctr < this->ns; ++vec_ctr){
    
    stride = 1;
    surface = this->volume / T;
    start = this->basis + vec_ctr * this->volume;
    gauge_start = u;
    
    aggregate.re = 0;
    aggregate.im = 0;
    
    for(i = 0; i < surface; ++i){
      boundary_D_0(tmp, start + i * stride, gauge_start + 8 * i * stride);
      for(j = 0; j < this->ns; ++j){
        result = block_scalar_prod(neighbours[0][i + j * surface], tmp, 1);
        aggregate[j].re += result.re;
        aggregate[j].im += result.im;
      }
    }
    memcpy(M, aggregate, this->ns); /* NOTE I'm maybe mucking up colums/rows here */
  }

  /* -T direction (-0) */
  M += this->ns * this->ns;
  for(vec_ctr = 0; vec_ctr < this->ns; ++vec_ctr){   
    
    start = this->basis + volume - T + vec_ctr * this->volume;
    gauge_start = u + 8 * (volume - T) + 1;
    
    aggregate.re = 0;
    aggregate.im = 0;
    
    for(i = 0; i < surface; ++i){
      boundary_D_1(tmp, start + i * stride, gauge_start + 8 * i * stride);
      for(j = 0; j < this->ns; ++j){
        result = block_scalar_prod(neighbours[1][i + j * surface], tmp, 1);
        aggregate[j].re += result.re;
        aggregate[j].im += result.im;
      }
    }
    memcpy(M, aggregate, this->ns); /* NOTE I'm maybe mucking up colums/rows here */
  }

  /* +X direction (+1) */
  /* NOTE As above, but make sure to include also readsize to be able to fully define the scan */
}

/* Uses a Modified Gram-Schmidt algorithm to orthonormalize little basis vectors */
void block_orthonormalize(void *parent){
  int i, j, k;
  spinor *current, *next, *iter;
  spinor orig, resbasis;
  complex coeff;
  complex scale;
  block *this = (block*)parent;

  scale.im = 0;
  for(i = 0; i < this->ns; ++i){
    /* rescale the current vector */
    current = this->basis + i * (this->volume + this->spinpad);
    scale.re = 1 / sqrt(block_two_norm(current, this->volume));
    for(iter = current; iter < current + this->volume; ++iter){
      orig = *iter; /* we can't alias */
      _spinor_mul_complex(*iter, scale, orig);
    }
    /* rescaling done, now subtract this direction from all vectors that follow */
    for(j = i + 1; j < this->ns; ++j){
      next = this->basis + j * (this->volume + this->spinpad);
      coeff = block_scalar_prod(current, next, this->volume);
      for(k = 0; k < this->volume; ++k){
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
    memcpy(reconstructed_field + (2 * ctr_t) * contig_block, g_blocks[0].basis + index * (VOLUME + g_blocks[0].spinpad),contig_block * sizeof(spinor));
    memcpy(reconstructed_field + (2 * ctr_t + 1) * contig_block, g_blocks[1].basis + index * (VOLUME + g_blocks[1].spinpad), contig_block * sizeof(spinor));
  }
  return reconstructed_field;
}

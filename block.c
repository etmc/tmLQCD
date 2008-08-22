/* $Id$ */
#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "global.h"
#include "D_psi.h"
#include "linalg/diff.h"
#include "xchange_lexicfield.h"
#include "block.h"

#define CALLOC_ERROR_CRASH {printf ("calloc errno : %d\n", errno); errno = 0; return 1;}

block * block_list = NULL;
static spinor * basis = NULL;
spinor * _edges = NULL;
static spinor * edges = NULL;
static su3 * u = NULL;
const int spinpad = 1;

int init_blocks() {
  int i,j;
  block_list = calloc(2, sizeof(block));
  if((void*)(basis = (spinor*)calloc(2 * g_N_s * (VOLUME / 2 + spinpad) + 2, sizeof(spinor))) == NULL) {
    CALLOC_ERROR_CRASH;
  }
  if((void*)(_edges = (spinor*)calloc(g_N_s * (1+2*VOLUME/T + 2*VOLUME/LX + 2*VOLUME/LY + 4*VOLUME/LZ), sizeof(spinor))) == NULL) {
    CALLOC_ERROR_CRASH;
  }
  if((void*)(u = (su3*)calloc(1+8*VOLUME, sizeof(su3))) == NULL) {
    CALLOC_ERROR_CRASH;
  }

#if ( defined SSE || defined SSE2 || defined SSE3)
  block_list[0].basis = (spinor*)(((unsigned long int)(basis)+ALIGN_BASE)&~ALIGN_BASE);
  edges = (spinor*)(((unsigned long int)(_edges)+ALIGN_BASE)&~ALIGN_BASE);
  block_list[0].u = (su3*)(((unsigned long int)(u)+ALIGN_BASE)&~ALIGN_BASE);
#else
  block_list[0].basis = basis;
  edges = _edges;
  block_list[0].u = u;
#endif
  block_list[1].basis = block_list[0].basis + g_N_s * (VOLUME / 2 + spinpad);
  block_list[1].u = block_list[0].u + 4*VOLUME;

  for (i = 0; i < 2; ++i) {
    block_list[i].volume = VOLUME/2;
    block_list[i].LX = LX;
    block_list[i].LY = LY;
    block_list[i].LZ = LZ/2;
    block_list[i].T = T;
    block_list[i].spinpad = spinpad;
    memcpy(block_list[i].mpilocal_coordinate, g_proc_coords, 4*sizeof(int));
    memcpy(block_list[i].coordinate, g_proc_coords, 3*sizeof(int));
    block_list[i].coordinate[3] = 2 * g_proc_coords[3] + i;

    if ((void*)(block_list[i].idx = calloc(8 * VOLUME/2, sizeof(int))) == NULL)
      CALLOC_ERROR_CRASH;

    for (j = 1; j < g_N_s+1; j++){ /* write a zero element at the end of every spinor */
      _spinor_null(block_list[i].basis[j * (VOLUME/2 + block_list[i].spinpad)]);
    }

    if ((void*)(block_list[i].neighbour_edges = calloc(8, sizeof(spinor *))) == NULL)
      CALLOC_ERROR_CRASH;

    block_list[i].neighbour_edges[0] = edges + i*(VOLUME/T + VOLUME/LX + VOLUME/LY + 2*VOLUME/LZ);
    block_list[i].neighbour_edges[1] = block_list[i].neighbour_edges[0] + VOLUME / (2 * T);
    block_list[i].neighbour_edges[2] = block_list[i].neighbour_edges[1] + VOLUME / (2 * T);
    block_list[i].neighbour_edges[3] = block_list[i].neighbour_edges[2] + VOLUME / (2 * LX);
    block_list[i].neighbour_edges[4] = block_list[i].neighbour_edges[3] + VOLUME / (2 * LX);
    block_list[i].neighbour_edges[5] = block_list[i].neighbour_edges[4] + VOLUME / (2 * LY);
    block_list[i].neighbour_edges[6] = block_list[i].neighbour_edges[5] + VOLUME / (2 * LY);
    block_list[i].neighbour_edges[7] = block_list[i].neighbour_edges[6] + VOLUME / (LZ);

    if ((void*)(block_list[i].little_dirac_operator = calloc(9 * g_N_s * g_N_s, sizeof(complex))) == NULL)
      CALLOC_ERROR_CRASH;
  }

  return 0;
}

int free_blocks() {
  int i;
  for(i = 0; i < 2; ++i) {
    free(block_list[i].neighbour_edges);
    free(block_list[i].little_dirac_operator);
  }
  free(u);
  free(_edges);
  free(basis);
  free(block_list);
  return 0;
}

int add_basis_field(int const index, spinor const *field) {
  int ctr_t;
  int contig_block = LZ / 2;
  for (ctr_t = 0; ctr_t < (VOLUME / LZ); ++ctr_t)
  {
    memcpy(block_list[0].basis + index * (VOLUME/2 + block_list[0].spinpad), field + (2 * ctr_t) * contig_block, contig_block * sizeof(spinor));
    memcpy(block_list[1].basis + index * (VOLUME/2 + block_list[1].spinpad), field + (2 * ctr_t + 1) * contig_block, contig_block * sizeof(spinor));
  }
  return 0;
}

int init_gauge_blocks() {
  /* Copies the existing gauge field on the processor into the two separate blocks in a form
  that is readable by the block Dirac operator. Specifically, in consecutive memory
  now +t,-t,+x,-x,+y,-y,+z,-z gauge links are stored. This requires double the storage in
  memory. */

  int x, y, z, t, ix, ix_new;
  for (t = 0; t < T; ++t) {
    for (x = 0; x < LX; ++x) {
      for (y = 0; y < LY; ++y) {
        for (z = 0; z < LZ/2; ++z) {
          ix = g_ipt[t][x][y][z];
          ix_new = 8 * (z + y * LZ/2 + x * LY * LZ/2 + t * LX * LY * LZ/2); /*su3 index on this block*/
          memcpy(block_list[0].u + ix_new, *g_gauge_field + g_iup[ix][0], sizeof(su3));
          memcpy(block_list[0].u + ix_new + 1, *g_gauge_field + g_idn[ix][0], sizeof(su3));
          memcpy(block_list[0].u + ix_new + 2, *g_gauge_field + g_iup[ix][1], sizeof(su3));
          memcpy(block_list[0].u + ix_new + 3, *g_gauge_field + g_idn[ix][1], sizeof(su3));
          memcpy(block_list[0].u + ix_new + 4, *g_gauge_field + g_iup[ix][2], sizeof(su3));
          memcpy(block_list[0].u + ix_new + 5, *g_gauge_field + g_idn[ix][2], sizeof(su3));
          memcpy(block_list[0].u + ix_new + 6, *g_gauge_field + g_iup[ix][3], sizeof(su3));
          memcpy(block_list[0].u + ix_new + 7, *g_gauge_field + g_idn[ix][3], sizeof(su3));
        }
        for (z = LZ/2; z < LZ; ++z) {
          ix = g_ipt[t][x][y][z];
          ix_new = 8 * (z - LZ/2 + y * LZ/2 + x * LY * LZ/2 + t * LX * LY * LZ/2); /*su3 index on this block, count anew in the z direction*/
          memcpy(block_list[1].u + ix_new, *g_gauge_field + g_iup[ix][0], sizeof(su3));
          memcpy(block_list[1].u + ix_new + 1, *g_gauge_field + g_idn[ix][0], sizeof(su3));
          memcpy(block_list[1].u + ix_new + 2, *g_gauge_field + g_iup[ix][1], sizeof(su3));
          memcpy(block_list[1].u + ix_new + 3, *g_gauge_field + g_idn[ix][1], sizeof(su3));
          memcpy(block_list[1].u + ix_new + 4, *g_gauge_field + g_iup[ix][2], sizeof(su3));
          memcpy(block_list[1].u + ix_new + 5, *g_gauge_field + g_idn[ix][2], sizeof(su3));
          memcpy(block_list[1].u + ix_new + 6, *g_gauge_field + g_iup[ix][3], sizeof(su3));
          memcpy(block_list[1].u + ix_new + 7, *g_gauge_field + g_idn[ix][3], sizeof(su3));
        }
      }
    }
  }
  return(0);
}

int init_geom_blocks() {
  int ix;
  int zstride = 1;
  int ystride = LZ/2;
  int xstride = LY * LZ/2;
  int tstride = LX * LY * LZ/2;
  int boundidx = 4 * VOLUME;
  for (ix = 0; ix < VOLUME/2; ++ix) {
    block_list[0].idx[8 * ix + 0] = (ix           >= VOLUME - tstride ? boundidx : ix + tstride);/* +t */
    block_list[0].idx[8 * ix + 1] = (ix           <  tstride          ? boundidx : ix - tstride);/* -t */
    block_list[0].idx[8 * ix + 2] = (ix % tstride >= LZ/2*LY*(LX-1)   ? boundidx : ix + xstride);/* +x */
    block_list[0].idx[8 * ix + 3] = (ix % tstride <  LZ/2*LY          ? boundidx : ix - xstride);/* -x */
    block_list[0].idx[8 * ix + 4] = (ix % xstride >= LZ/2 * (LY - 1)  ? boundidx : ix + ystride);/* +y */
    block_list[0].idx[8 * ix + 5] = (ix % xstride <  LZ/2             ? boundidx : ix - ystride);/* -y */
    block_list[0].idx[8 * ix + 6] = (ix % ystride == LZ/2 - 1         ? boundidx : ix + zstride);/* +z */
    block_list[0].idx[8 * ix + 7] = (ix % ystride == 0                ? boundidx : ix - zstride);/* -z */
  }
  memcpy(block_list[1].idx, block_list[0].idx, 8 * VOLUME/2 * sizeof(int));
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

void block_compute_little_D_diagonal(block *parent) {
  int i,j;
  /* we need working space, where do we best allocate it? */
  spinor * tmp;
  complex * M;
/*   tmp = (spinor *)calloc((VOLUME/2 + parent->spinpad), sizeof(spinor)); */
  tmp = g_spinor_field[DUM_SOLVER];
  M = parent->little_dirac_operator;

  for(i = 0; i < g_N_s; i++){
    Block_D_psi(parent, tmp, parent->basis + i * (parent->volume + parent->spinpad));
    for(j = 0; j < g_N_s; j++){
      /* order correct ? */
      M[i * g_N_s + j] = block_scalar_prod(parent->basis + j * (parent->volume + parent->spinpad), tmp, parent->volume);
    }
  }
/*   free(tmp); */
}

void surface_D_apply_contract(block *parent, int surface, spinor* offset, int stride, 
                              int step_size, su3 *gauge_start, int direction){
  int i, j, k;
  complex *iter;
  complex * aggregate;
  complex result;
  static spinor tmp;
  void (*boundary_D[8])(spinor * const r, spinor * const s, su3 *u) =
    {boundary_D_0, boundary_D_1, boundary_D_2, boundary_D_3, boundary_D_4, boundary_D_5, boundary_D_6, boundary_D_7};

  aggregate = (complex*)malloc(g_N_s*sizeof(complex));
  for (iter = aggregate; iter != aggregate + g_N_s; ++iter){
    iter->re = 0;
    iter->im = 0;
  }

  for(i = 0; i < (surface / step_size); ++i){
    for(j = 0; j < step_size; ++j){
      boundary_D[direction](&tmp, offset + i * stride + j, gauge_start + 8 * (i * stride + j) + direction);
      for(k = 0; k < g_N_s; ++k){
        result = block_scalar_prod(&(parent->neighbour_edges[direction][i * step_size + j + k * surface]), &tmp, 1);
        aggregate[k].re += result.re;
        aggregate[k].im += result.im;
      }
    }
  }
  memcpy(parent->little_dirac_operator + (direction + 1) * (g_N_s * g_N_s) , aggregate, g_N_s * sizeof(complex));
  /* NOTE Is this the correct location within M eventually? */
  free(aggregate);
}

/* what happens if this routine is called in a one dimensional parallelisation? */
/* or even serially ?                                                           */
void block_compute_little_D_offdiagonal(block *parent) {
/*   Here we need to multiply the boundary with the corresponding  */
/*   U and gamma_i and take the scalar product then */
/*   NOTE I assume the boundaries are available in parent->neighbour_edges*/
  int vec_ctr, surface;
  spinor *offset;
  int stride, step_size;
  su3 *gauge_start;

  /* +T direction (+0) */
  for(vec_ctr = 0; vec_ctr < g_N_s; ++vec_ctr){
    surface   = LX * LY * LZ / 2; /* Corrected for snipped Z direction */
    stride    = surface;
    step_size = surface;

    offset = parent->basis + (T - 1) * stride + vec_ctr * (parent->volume + parent->spinpad);
    gauge_start = parent->u + 8 * (parent->volume - surface);

    surface_D_apply_contract(parent, surface, offset, stride, step_size, gauge_start, 0);
  }

  /* -T direction (-0) */
  for(vec_ctr = 0; vec_ctr < g_N_s; ++vec_ctr){
    offset = parent->basis + vec_ctr * (parent->volume + parent->spinpad);
    gauge_start = parent->u;

    surface_D_apply_contract(parent, surface, offset, stride, step_size, gauge_start, 1);
  }

  /* +X direction (+1) */
  for(vec_ctr = 0; vec_ctr < g_N_s; ++vec_ctr){
    surface =   LZ * LY * T  / 2;
    stride =    LZ * LY * LX / 2;
    step_size = LZ * LY      / 2;

    offset = parent->basis + (LX - 1) * stride + vec_ctr * (parent->volume + parent->spinpad);
    gauge_start = parent->u + (LX - 1) * stride;  /* Pick up -X gauge throughout */

    surface_D_apply_contract(parent, surface, offset, stride, step_size, gauge_start, 2);
  }

  /* -X direction (-1) */
  for(vec_ctr = 0; vec_ctr < g_N_s; ++vec_ctr){
    offset = parent->basis + vec_ctr * (parent->volume + parent->spinpad);
    gauge_start = parent->u;

    surface_D_apply_contract(parent, surface, offset, stride, step_size, gauge_start, 3);
  }

  /* +Y direction (+2) */
  for(vec_ctr = 0; vec_ctr < g_N_s; ++vec_ctr){
    surface =   LZ * LX * T / 2;
    stride =    LZ * LY     / 2;
    step_size = LZ          / 2;

    offset = parent->basis + (LY - 1) * stride + vec_ctr * (parent->volume + parent->spinpad);
    gauge_start = parent->u + (LY - 1) * stride;

    surface_D_apply_contract(parent, surface, offset, stride, step_size, gauge_start, 4);
  }

  /* -Y direction (-2) */
  for(vec_ctr = 0; vec_ctr < g_N_s; ++vec_ctr){
    offset = parent->basis + vec_ctr * (parent->volume + parent->spinpad);
    gauge_start = parent->u;

    surface_D_apply_contract(parent, surface, offset, stride, step_size, gauge_start, 5);
  }

  /* +Z direction (+3) */
  for(vec_ctr = 0; vec_ctr < g_N_s; ++vec_ctr){
    surface =   LY * LX * T   ;
    stride =    LZ         / 2;
    step_size = 1             ;

    offset = parent->basis + (LZ - 1) * stride + vec_ctr * (parent->volume + parent->spinpad);
    gauge_start = parent->u + (LZ - 1) * stride;

    surface_D_apply_contract(parent, surface, offset, stride, step_size, gauge_start, 6);
  }

  /* -Z direction (-3) */
  for(vec_ctr = 0; vec_ctr < g_N_s; ++vec_ctr){
    offset = parent->basis + vec_ctr * (parent->volume + parent->spinpad);
    gauge_start = parent->u;

    surface_D_apply_contract(parent, surface, offset, stride, step_size, gauge_start, 7);
  }
}

/* Uses a Modified Gram-Schmidt algorithm to orthonormalize little basis vectors */
void block_orthonormalize(block *parent) {
  int i, j, k;
  spinor *current, *next, *iter;
  spinor orig, resbasis;
  complex coeff;
  complex scale;

  scale.im = 0;
  for(i = 0; i < parent->ns; ++i){
    /* rescale the current vector */
    current = parent->basis + i * (parent->volume + parent->spinpad);
    scale.re = 1 / sqrt(block_two_norm(current, parent->volume));
    for(iter = current; iter < current + parent->volume; ++iter){
      orig = *iter; /* we can't alias */
      _spinor_mul_complex(*iter, scale, orig);
    }
    /* rescaling done, now subtract this direction from all vectors that follow */
    for(j = i + 1; j < parent->ns; ++j){
      next = parent->basis + j * (parent->volume + parent->spinpad);
      coeff = block_scalar_prod(current, next, parent->volume);
      for(k = 0; k < parent->volume; ++k){
        _spinor_mul_complex(resbasis, coeff, *(current + k));
        orig = *(next + k);
        diff(next + k, &orig, &resbasis, 1); /* NOTE We could also subtract the complete spinors from each other */
      }
    }
  }
  return;
}

/* what happens if this routine is called in a one dimensional parallelisation? */
/* or even serially ?                                                           */
/* it shouldn't be */
void blocks_exchange_edges() {
  int bl_cnt, vec_cnt, div_cnt;
  int div_size = LZ / 2;
  spinor *scratch, *offset_up, *offset_dn;

  int offsets[8];  int surfaces[4];
  offsets[0] = VOLUME;
  offsets[1] = (T + 2) * LX * LY * LZ;
  offsets[2] = VOLUME + 2 * LZ * (LX * LY + T * LY);
  offsets[3] = VOLUME + 2 * LZ * (LX * LY + T * LY + T * LX);
  offsets[4] = (T + 1) * LX * LY * LZ;
  offsets[5] = (T + 2) * LX * LY * LZ + T * LY * LZ;
  offsets[6] = VOLUME + 2 * LZ * (LX * LY + T * LY) + T * LX * LZ;
  offsets[7] = VOLUME + 2 * LZ * (LX * LY + T * LY + T * LX) + T * LX * LY;

  surfaces[0] = VOLUME / T;
  surfaces[1] = VOLUME / LX;
  surfaces[2] = VOLUME / LY;
  surfaces[3] = VOLUME / LZ;

  /* for a full spinor field we need VOLUMEPLUSRAND                 */
  /* because we use the same geometry as for the                    */
  /* gauge field                                                    */
  /* It is VOLUME + 2*LZ*(LY*LX + T*LY + T*LX) + 4*LZ*(LY + T + LX) */
  scratch = calloc(VOLUMEPLUSRAND, sizeof(spinor));
  for (vec_cnt = 0; vec_cnt < g_N_s; ++vec_cnt){
    block_reconstruct_global_field(vec_cnt, scratch);
    xchange_lexicfield(scratch);

    for (bl_cnt = 0; bl_cnt < 6; ++bl_cnt){
      for (div_cnt = 0; div_cnt < (surfaces[bl_cnt / 2] / (2 * div_size)); ++div_cnt) {
        memcpy(block_list[0].neighbour_edges[bl_cnt] + vec_cnt * surfaces[bl_cnt / 2] + 2 * div_cnt * div_size,
               scratch + offsets[bl_cnt] + 2 * div_cnt * div_size, div_size * sizeof(spinor));
        memcpy(block_list[1].neighbour_edges[bl_cnt] + vec_cnt * surfaces[bl_cnt / 2] + (2 * div_cnt + 1) * div_size,
               scratch + offsets[bl_cnt]+ div_cnt * (2 * div_cnt + 1) * div_size, div_size * sizeof(spinor));
      }
    }

   memcpy(block_list[1].neighbour_edges[6] + vec_cnt * surfaces[3], scratch + offsets[6], surfaces[3] * sizeof(spinor));
   if (!g_proc_id){printf("[DEBUG] vec_cnt %u, %u\n", vec_cnt, vec_cnt * surfaces[3]);}
    memcpy(block_list[0].neighbour_edges[7] + vec_cnt * surfaces[3], scratch + offsets[7], surfaces[3] * sizeof(spinor));

    offset_up = block_list[1].basis + (LZ - 1) * LZ + vec_cnt * (block_list[1].volume + block_list[1].spinpad);
    offset_dn = block_list[0].basis + vec_cnt * (block_list[0].volume + block_list[0].spinpad);
    for(div_cnt = 0; div_cnt < surfaces[3]; ++div_cnt){
      memcpy(block_list[0].neighbour_edges[6] + div_cnt, offset_up + div_cnt * LZ, sizeof(spinor));
      memcpy(block_list[1].neighbour_edges[7] + div_cnt, offset_dn + div_cnt * LZ, sizeof(spinor));
    }
  }

  free(scratch);
}

/* Reconstructs a global field from the little basis of two blocks */
void block_reconstruct_global_field(const int index, spinor * const reconstructed_field) {
  int ctr_t;
  int contig_block = LZ / 2;
  for (ctr_t = 0; ctr_t < (VOLUME / contig_block); ++ctr_t) {
    memcpy(reconstructed_field + (2 * ctr_t) * contig_block, 
           block_list[0].basis + index * (VOLUME + block_list[0].spinpad) + ctr_t * contig_block,
	   contig_block * sizeof(spinor));
    memcpy(reconstructed_field + (2 * ctr_t + 1) * contig_block, 
           block_list[1].basis + index * (VOLUME + block_list[1].spinpad) + ctr_t * contig_block, 
	   contig_block * sizeof(spinor));
  }
  return;
}

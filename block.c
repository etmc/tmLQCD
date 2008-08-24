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
#include "linalg_eo.h"
#include "xchange_lexicfield.h"
#include "block.h"

#define CALLOC_ERROR_CRASH {printf ("calloc errno : %d\n", errno); errno = 0; return 1;}

int init_blocks_gaugefield();
int init_blocks_geometry();
double block_two_norm(spinor * const R, const int N);

int **** block_ipt;
int *** bipt__;
int ** bipt_;
int * bipt;

static void (*boundary_D[8])(spinor * const r, spinor * const s, su3 *u) =
{boundary_D_0, boundary_D_1, boundary_D_2, boundary_D_3, boundary_D_4, boundary_D_5, boundary_D_6, boundary_D_7};


block * block_list = NULL;
static spinor * basis = NULL;
static spinor * _edges = NULL;
static spinor * edges = NULL;
static su3 * u = NULL;
const int spinpad = 1;
static int block_init = 0;

int init_blocks() {
  int i,j;
  free_blocks();
  block_init = 1;
  block_list = calloc(2, sizeof(block));
  if((void*)(basis = (spinor*)calloc(2 * g_N_s * (VOLUME / 2 + spinpad) + 1, sizeof(spinor))) == NULL) {
    CALLOC_ERROR_CRASH;
  }
  if((void*)(_edges = (spinor*)calloc(g_N_s * (1+2*VOLUME/T + 2*VOLUME/LX + 2*VOLUME/LY + 4*VOLUME/LZ), sizeof(spinor))) == NULL) {
    CALLOC_ERROR_CRASH;
  }
  if((void*)(u = (su3*)calloc(1+8*VOLUME, sizeof(su3))) == NULL) {
    CALLOC_ERROR_CRASH;
  }
  for(i = 0; i < 2; i++) {
    block_list[i].basis = (spinor**)calloc(g_N_s, sizeof(spinor*));
  }

#if ( defined SSE || defined SSE2 || defined SSE3)
  block_list[0].basis[0] = (spinor*)(((unsigned long int)(basis)+ALIGN_BASE)&~ALIGN_BASE);
  edges = (spinor*)(((unsigned long int)(_edges)+ALIGN_BASE)&~ALIGN_BASE);
  block_list[0].u = (su3*)(((unsigned long int)(u)+ALIGN_BASE)&~ALIGN_BASE);
#else
  block_list[0].basis[0] = basis;
  edges = _edges;
  block_list[0].u = u;
#endif
  block_list[1].basis[0] = block_list[0].basis[0] + g_N_s * (VOLUME / 2 + spinpad);
  for(i = 1; i < g_N_s; i++) {
    block_list[0].basis[i] = block_list[0].basis[i-1] + VOLUME / 2 + spinpad;
    block_list[1].basis[i] = block_list[1].basis[i-1] + VOLUME / 2 + spinpad;
  }
  block_list[1].u = block_list[0].u + 4*VOLUME;

  if((void*)(block_ipt = (int****)calloc(T+4,sizeof(int*))) == NULL) return(5);
  if((void*)(bipt__ = (int***)calloc ((T+4)*(LX+4), sizeof(int*))) == NULL) return(4);
  if((void*)(bipt_ = (int**)calloc((T+4)*(LX+4)*(LY+4), sizeof(int*))) == NULL) return(3);
  if((void*)(bipt = (int*)calloc((T+4)*(LX+4)*(LY+4)*(LZ/2+4), sizeof(int))) == NULL) return(8);
  bipt_[0] = bipt;
  bipt__[0] = bipt_;
  block_ipt[0] = bipt__;
  for(i = 1; i < (T+4)*(LX+4)*(LY+4); i++){
    bipt_[i] = bipt_[i-1]+(LZ/2+4);
  }
  for(i = 1; i < (T+4)*(LX+4); i++){
    bipt__[i] = bipt__[i-1]+(LY+4);
  }
  for(i = 1; i < (T+4); i++){
    block_ipt[i] = block_ipt[i-1]+(LX+4);
  }

  for (i = 0; i < 2; ++i) {
    block_list[i].id = i;
    block_list[i].volume = VOLUME/2;
    block_list[i].LX = LX;
    block_list[i].LY = LY;
    block_list[i].LZ = LZ/2;
    block_list[i].T = T;
    block_list[i].ns = g_N_s;
    block_list[i].spinpad = spinpad;
    for (j = 0 ; j < 6; ++j) {
      #ifdef MPI
        block_list[i].mpilocal_neighbour[j] = (g_nb_list[j] == g_cart_id) ? i : -1;
      #else
        block_list[i].mpilocal_neighbour[j] = i;
      #endif
    }
#ifdef MPI
    block_list[i].mpilocal_neighbour[6] = (i == 0 ? 1 : (g_nb_list[j] == g_cart_id) ? 0 : -1);
    block_list[i].mpilocal_neighbour[7] = (i == 1 ? 0 : (g_nb_list[j] == g_cart_id) ? 1 : -1);
#else
    block_list[i].mpilocal_neighbour[6] = (i == 0 ? 1 : 0);
    block_list[i].mpilocal_neighbour[7] = (i == 0 ? 1 : 0);
#endif
    memcpy(block_list[i].mpilocal_coordinate, g_proc_coords, 4*sizeof(int));
    memcpy(block_list[i].coordinate, g_proc_coords, 3*sizeof(int));
    block_list[i].coordinate[3] = 2 * g_proc_coords[3] + i;

    if ((void*)(block_list[i].idx = calloc(8 * VOLUME/2, sizeof(int))) == NULL)
      CALLOC_ERROR_CRASH;

    for (j = 0; j < g_N_s; j++){ /* write a zero element at the end of every spinor */
      _spinor_null(block_list[i].basis[j][VOLUME/2]);
    }

    if ((void*)(block_list[i].little_dirac_operator = calloc(9 * g_N_s * g_N_s, sizeof(complex))) == NULL)
      CALLOC_ERROR_CRASH;
  }
  init_blocks_geometry();
  init_blocks_gaugefield();

  return 0;
}

int free_blocks() {
  int i;
  if(block_init == 1) {
    for(i = 0; i < 2; ++i) {
      free(block_list[i].basis);
      free(block_list[i].little_dirac_operator);
    }
    free(block_ipt);
    free(bipt__);
    free(bipt_);
    free(bipt);
    free(u);
    free(_edges);
    free(basis);
    free(block_list);
    block_init = 0;
  }
  return 0;
}

int add_basis_field(int const index, spinor const *field) {
  int ctr_t;
  int contig_block = LZ / 2;
  for (ctr_t = 0; ctr_t < (VOLUME / LZ); ctr_t++) {
    memcpy(block_list[0].basis[index] + ctr_t*contig_block, 
	   field + (2 * ctr_t) * contig_block, contig_block * sizeof(spinor));
    memcpy(block_list[1].basis[index] + ctr_t*contig_block, 
	   field + (2 * ctr_t + 1) * contig_block, contig_block * sizeof(spinor));
  }
  if(g_proc_id == 0 && g_debug_level > 4) {
    printf("basis norm = %1.3e\n", block_two_norm(block_list[1].basis[index], block_list[0].volume));
  }
      
  return 0;
}

int init_blocks_gaugefield() {
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

int check_blocks_geometry(block * blk) {
  int i, k=0;
  int * itest;
  int * ipt;
  
  ipt = blk->idx;
  itest = (int*)calloc(blk->volume + blk->spinpad, sizeof(int));
  for(i = 0; i < 8*blk->volume; i++) {
    if(*ipt > blk->volume + blk->spinpad-1 || *ipt < 0) {
      if(g_proc_id == 0) {
	printf("error in block geometry! ipt = %d dir = %d i = %d of %d\n", 
	       (*ipt), i%8, i/8, blk->volume + blk->spinpad);
      }
    }
    
    itest[*(ipt++)]++;
  }
  
  for(i = 0; i < blk->volume; i++) {
    k += itest[i];
    if(itest[i] < 1 || itest[i] > 8) {
      if(g_proc_id == 0) {
	printf("error in block geometry, itest[%d] = %d\n", i, itest[i]);      
      }
    }
  }
  
  if(itest[blk->volume + blk->spinpad-1] != 2*(blk->LX*blk->LY*blk->LZ+blk->T*blk->LX*blk->LY+blk->T*blk->LY*blk->LZ+blk->T*blk->LX*blk->LZ)) {
    if(g_proc_id == 0){
      printf("error in block geometry, boundary points wrong %d != %d\n", 
	    itest[blk->volume + blk->spinpad-1], 2*(blk->LX*blk->LY*blk->LZ+blk->T*blk->LX*blk->LY+blk->T*blk->LY*blk->LZ+blk->T*blk->LX*blk->LZ));
    }
  }
  k+= itest[blk->volume + blk->spinpad-1];
  if(k != 8*blk->volume) {
    if(g_proc_id == 0){
      printf("error in block geometry, total number of points wrong %d != %d\n", 
	    k, 8*blk->volume);
    }
  }
  free(itest);
  if(g_proc_id == 0) {
    printf("# block geometry checked successfully for block %d !\n", blk->id);
  }
  return(0);
}

int init_blocks_geometry() {
  int ix, x, y, z, t;
  int zstride = 1;
  int ystride = LZ/2;
  int xstride = LY * LZ/2;
  int tstride = LX * LY * LZ/2;
  int boundidx = VOLUME/2;
  for (ix = 0; ix < VOLUME/2; ++ix) {
    block_list[0].idx[8 * ix + 0] = (ix           >= VOLUME/2 - tstride ? boundidx : ix + tstride);/* +t */
    block_list[0].idx[8 * ix + 1] = (ix           <  tstride          ? boundidx : ix - tstride);/* -t */
    block_list[0].idx[8 * ix + 2] = (ix % tstride >= LZ/2*LY*(LX-1)   ? boundidx : ix + xstride);/* +x */
    block_list[0].idx[8 * ix + 3] = (ix % tstride <  LZ/2*LY          ? boundidx : ix - xstride);/* -x */
    block_list[0].idx[8 * ix + 4] = (ix % xstride >= LZ/2 * (LY - 1)  ? boundidx : ix + ystride);/* +y */
    block_list[0].idx[8 * ix + 5] = (ix % xstride <  LZ/2             ? boundidx : ix - ystride);/* -y */
    block_list[0].idx[8 * ix + 6] = (ix % ystride == LZ/2 - 1         ? boundidx : ix + zstride);/* +z */
    block_list[0].idx[8 * ix + 7] = (ix % ystride == 0                ? boundidx : ix - zstride);/* -z */
  }
  memcpy(block_list[1].idx, block_list[0].idx, 8 * VOLUME/2 * sizeof(int));
  for(ix = 0; ix < 2; ix++) {
    zstride = check_blocks_geometry(&block_list[ix]);
  }
  ix = 0;
  for(t = 0; t < T; t++) {
    for(x = 0; x < LX; x++) {
      for(y = 0; y < LY; y++) {
	for(z = 0; z < LZ/2; z++) {
	  block_ipt[t][x][y][z] = ix;
	  ix++;
	}
      }
    }
  }

  return 0;
}

/* the following should be somewhere else ... */


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
    Block_D_psi(parent, tmp, parent->basis[i]);
    for(j = 0; j < g_N_s; j++){
      /* order correct ? */
/*       M[i * g_N_s + j] = block_scalar_prod(parent->basis[j], parent->basis[i], parent->volume); */
      M[i * g_N_s + j] = block_scalar_prod(parent->basis[j], tmp, parent->volume);
    }
  }
/*   free(tmp); */
}

/* Uses a Modified Gram-Schmidt algorithm to orthonormalize little basis vectors */
void block_orthonormalize(block *parent) {
  int i, j;
  complex coeff;
  double scale;

  for(i = 0; i < g_N_s; ++i){
    /* rescale the current vector */
    scale = 1. / sqrt(block_two_norm(parent->basis[i], parent->volume));
    mul_r(parent->basis[i], scale, parent->basis[i], parent->volume);

    /* rescaling done, now subtract this direction from all vectors that follow */
    for(j = i + 1; j < g_N_s; ++j){
      coeff = block_scalar_prod(parent->basis[j], parent->basis[i], parent->volume);
      assign_diff_mul(parent->basis[j], parent->basis[i], coeff, parent->volume);
    }
  }

  if(g_debug_level > 4) {
    for(i = 0; i < g_N_s; i++) {
      for(j = 0; j < g_N_s; j++) {
        coeff = block_scalar_prod(parent->basis[j], parent->basis[i], parent->volume);
        if(g_proc_id == 0) printf("basis id = %d <%d, %d> = %1.3e +i %1.3e\n", 
				  parent->id, j, i, coeff.re, coeff.im);
      }
    }
  }
  return;
}

/* what happens if this routine is called in a one dimensional parallelisation? */
/* or even serially ?                                                           */
void block_compute_little_D_offdiagonal(){
  spinor *scratch, * temp;
  spinor *r, *s;
  su3 * u;
  int x, y, z, t, ix, iy, i, j, k, pm, mu;

  /* for a full spinor field we need VOLUMEPLUSRAND                 */
  /* because we use the same geometry as for the                    */
  /* gauge field                                                    */
  /* It is VOLUME + 2*LZ*(LY*LX + T*LY + T*LX) + 4*LZ*(LY + T + LX) */
  scratch = calloc(VOLUMEPLUSRAND, sizeof(spinor));
  temp =  calloc(VOLUMEPLUSRAND - VOLUME, sizeof(spinor));
  for (i = 0; i < g_N_s; i++){
    block_reconstruct_global_field(i, scratch);

#ifdef MPI
    xchange_lexicfield(scratch);
#endif

    /* +- t */
    for(pm = 0; pm < 2; pm++) {
      if(pm == 0) t = T-1;
      else t = 0;
      mu = 0;

      r = temp;
      for(x = 0; x < LX; x++) {
	for(y = 0; y < LY; y++) {
	  for(z = 0; z < LZ; z++) {
	    ix = g_ipt[t][x][y][z];
	    if(pm == 0) {
	      s = &scratch[ g_iup[ ix ][mu] ];
	      u = &g_gauge_field[ ix ][mu];
	    }
	    else {
	      s = &scratch[ g_idn[ ix ][mu] ];
	      u = &g_gauge_field[ g_idn[ix][mu] ][mu];
	    }
	    boundary_D[pm](r, s, u);
	    r++;
	  }
	}
      }
      r = temp;
      /* now all the scalar products */
      for(j = 0; j < g_N_s; j++) {
	for(x = 0; x < LX; x++) {
	  for(y = 0; y < LY; y++) {
	    for(k = 0; k < 2; k++) {
	      iy = j+i*g_N_s + (pm+1)*g_N_s*g_N_s;
	      _complex_zero(block_list[k].little_dirac_operator[ iy ]);
	      for(z = 0; z < LZ/2; z++) {
		ix = block_ipt[t][x][y][z];
		s = &block_list[k].basis[j][ ix ];
		_add_complex(block_list[k].little_dirac_operator[ iy ], block_scalar_prod(s, r, 1));
		r++;
	      }
	    }
	  }
	}
      }
    }

    /* +- x */
    for(pm = 2; pm < 4; pm++) {
      if(pm == 2) x = LX-1;
      else x = 0;
      mu = 1;

      r = temp;
      for(t = 0; t < t; t++) {
	for(y = 0; y < LY; y++) {
	  for(z = 0; z < LZ; z++) {
	    ix = g_ipt[t][x][y][z];
	    if(pm == 2) {
	      s = &scratch[ g_iup[ ix ][mu] ];
	      u = &g_gauge_field[ ix ][mu];
	    }
	    else {
	      s = &scratch[ g_idn[ ix ][mu] ];
	      u = &g_gauge_field[ g_idn[ix][mu] ][mu];
	    }
	    boundary_D[pm](r, s, u);
	    r++;
	  }
	}
      }
      r = temp;
      /* now all the scalar products */
      for(j = 0; j < g_N_s; j++) {
	for(t = 0; t < T; t++) {
	  for(y = 0; y < LY; y++) {
	    for(k = 0; k < 2; k++) {
	      iy = j+i*g_N_s + (pm+1)*g_N_s*g_N_s;
	      _complex_zero(block_list[k].little_dirac_operator[ iy ]);
	      for(z = 0; z < LZ/2; z++) {
		ix = block_ipt[t][x][y][z];
		s = &block_list[k].basis[j][ ix ];
		_add_complex(block_list[k].little_dirac_operator[ iy ], block_scalar_prod(s, r, 1));
		r++;
	      }
	    }
	  }
	}
      }
    }

    /* +- y */
    for(pm = 4; pm < 6; pm++) {
      if(pm == 4) y = LY-1;
      else y = 0;
      mu = 2;

      r = temp;
      for(t = 0; t < t; t++) {
	for(x = 0; x < LX; x++) {
	  for(z = 0; z < LZ; z++) {
	    ix = g_ipt[t][x][y][z];
	    if(pm == 4) {
	      s = &scratch[ g_iup[ ix ][mu] ];
	      u = &g_gauge_field[ ix ][mu];
	    }
	    else {
	      s = &scratch[ g_idn[ ix ][mu] ];
	      u = &g_gauge_field[ g_idn[ix][mu] ][mu];
	    }
	    boundary_D[pm](r, s, u);
	    r++;
	  }
	}
      }
      r = temp;
      /* now all the scalar products */
      for(j = 0; j < g_N_s; j++) {
	for(t = 0; t < T; t++) {
	  for(x = 0; x < LX; x++) {
	    for(k = 0; k < 2; k++) {
	      iy = j+i*g_N_s + (pm+1)*g_N_s*g_N_s;
	      _complex_zero(block_list[k].little_dirac_operator[ iy ]);
	      for(z = 0; z < LZ/2; z++) {
		ix = block_ipt[t][x][y][z];
		s = &block_list[k].basis[j][ ix ];
		_add_complex(block_list[k].little_dirac_operator[ iy ], block_scalar_prod(s, r, 1));
		r++;
	      }
	    }
	  }
	}
      }
    }

    /* z is different */
    /* +-z */
    for(pm = 6; pm < 8; pm++) {
      z = (pm == 6) ? LZ - 1 : 0;
      mu = 3;

      r = temp;
      for(t = 0; t < t; ++t) {
        for(x = 0; x < LX; ++x) {
          for(y = 0; y < LY; ++y) {
            ix = g_ipt[t][x][y][z];
            if(pm == 6) {
              s = &scratch[ g_iup[ ix ][mu] ];
              u = &g_gauge_field[ ix ][mu];
            }
            else {
              s = &scratch[ g_idn[ ix ][mu] ];
              u = &g_gauge_field[ g_idn[ix][mu] ][mu];
            }
            boundary_D[pm](r, s, u);
            r++;
          }
        }
      }
      r = temp;
      /* now for all the MPI scalar products (outer edges) */
      for(j = 0; j < g_N_s; ++j) {
        for(t = 0; t < T; ++t) {
          for(x = 0; x < LX; ++x) {
            for(y = 0; y < LY; ++y){
              iy = j + i * g_N_s + (pm + 1) * g_N_s * g_N_s;
              _complex_zero(block_list[pm % 2].little_dirac_operator[ iy ]);
              ix = block_ipt[t][x][y][z];
              s = &block_list[pm % 2].basis[j][ ix ];
              _add_complex(block_list[pm % 2].little_dirac_operator[ iy ], block_scalar_prod(s, r, 1));
              r++;
            }
          }
        }
      }
    }

    /* and finally the residual inner edges - new D calculations needed here */
    /* this is confusing enough as is, so I unrolled the pm loop for this part */
    z =  LZ / 2 - 1; /* pm = 6 */
    r = temp;
    for(t = 0; t < t; ++t) {
      for(x = 0; x < LX; ++x) {
        for(y = 0; y < LY; ++y) {
          ix = g_ipt[t][x][y][z];
          iy = block_ipt[t][x][y][z]; /* highest edge of low block needed */
          s = &block_list[ 0 ].basis[ i ][ iy ];
          u = &g_gauge_field[ ix ][mu];
          boundary_D[6](r, s, u);
          r++;
        }
      }
    }
    r = temp;

    /* Now contract with the lowest edge of the high block and store */
    for(j = 0; j < g_N_s; ++j) {
      for(t = 0; t < T; ++t) {
        for(x = 0; x < LX; ++x) {
          for(y = 0; y < LY; ++y){
            iy = j + i * g_N_s + (6 + 1) * g_N_s * g_N_s;
            _complex_zero(block_list[1].little_dirac_operator[ iy ]);
            ix = block_ipt[t][x][y][0];
            s = &block_list[1].basis[j][ ix ];
            _add_complex(block_list[1].little_dirac_operator[ iy ], block_scalar_prod(s, r, 1));
            r++;
          }
        }
      }
    }

    z =  LZ / 2; /* pm = 7 */
    r = temp;
    for(t = 0; t < t; ++t) {
      for(x = 0; x < LX; ++x) {
        for(y = 0; y < LY; ++y) {
          ix = g_ipt[t][x][y][z];
          iy = block_ipt[t][x][y][0];  /* lowest edge of high block needed */
          s = &block_list[1].basis[ i ][ iy ];
          u = &g_gauge_field[ g_idn[ ix ][ mu ] ][mu];
          boundary_D[7](r, s, u);
          r++;
        }
      }
    }
    r = temp;

    /* Now contract with the highest edge of the low block and store */
    for(j = 0; j < g_N_s; ++j) {
      for(t = 0; t < T; ++t) {
        for(x = 0; x < LX; ++x) {
          for(y = 0; y < LY; ++y){
            iy = j + i * g_N_s + (7 + 1) * g_N_s * g_N_s;
            _complex_zero(block_list[0].little_dirac_operator[ iy ]);
            ix = block_ipt[t][x][y][z - 1]; /* z - 1 = LZ / 2 - 1 */
            s = &block_list[0].basis[j][ ix ];
            _add_complex(block_list[0].little_dirac_operator[ iy ], block_scalar_prod(s, r, 1));
            r++;
          }
        }
      }
    }
  }

  free(scratch);
  free(temp);

  if(g_debug_level > 4) {
    if (g_N_s <= 5 && !g_cart_id){
      printf("\n  *** CHECKING LITTLE D ***\n");
      printf("\n  ** node 0, lower block **\n");
      for (i = 0; i < 9 * g_N_s; ++i){
        printf(" [ ");
        for (j = 0; j < g_N_s; ++j){
          printf(" %2.2E + %2.2f i", block_list->little_dirac_operator[i * g_N_s + j].re,  block_list->little_dirac_operator[i * g_N_s + j].im);
          if (j != g_N_s){
            printf(", ");
          }
        }
        printf(" ]\n\n");
      }
    }
  }

  return;
}

/* Reconstructs a global field from the little basis of two blocks */
void block_reconstruct_global_field(const int index, spinor * const reconstructed_field) {
  int ctr_t;
  int contig_block = LZ / 2;
  for (ctr_t = 0; ctr_t < (block_list[0].volume / contig_block); ++ctr_t) {
    memcpy(reconstructed_field + (2 * ctr_t) * contig_block, 
           block_list[0].basis[index] + ctr_t * contig_block,
	   contig_block * sizeof(spinor));
    memcpy(reconstructed_field + (2 * ctr_t + 1) * contig_block, 
           block_list[1].basis[index] + ctr_t * contig_block, 
	   contig_block * sizeof(spinor));
  }
  return;
}

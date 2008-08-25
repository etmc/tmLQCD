/* $Id$ */
#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#ifdef MPI
# include <mpi.h>
#endif
#include "global.h"
#include "start.h"
#include "complex.h"
#include "block.h"
#include "linalg/blas.h"
#include "D_psi.h"
#include "little_D.h"
#include "block.h"
#include "linalg_eo.h"
#include "dfl_projector.h"


/* Break up full volume spinor to blocks
 * loop over block.basis
 * compute inner product and store as complex vector
 * compute A^-1 * complex vector
 * loop over block.basis
 * compute sum of basis vectors times complex element
 * create global vector */

/* this is phi_k A^{-1}_{kl} (phi_k, in) */
void project(spinor * const out, spinor * const in) {
  int i,j,ctr_t, iter;
  spinor **psi;
  int contig_block = LZ / 2;
  int vol = block_list[0].volume;
  complex *inprod;
  complex *invvec;

  psi = calloc(4, sizeof(spinor*)); /*block local version of global spinor */
  inprod = calloc(2 * 9 * g_N_s, sizeof(complex)); /*inner product of spinors with bases */
  invvec = calloc(2 * 9 * g_N_s, sizeof(complex)); /*inner product of spinors with bases */

  /* no loop below because further down we also don't take this cleanly into account */
  psi[0] = calloc(VOLUME, sizeof(spinor));
  psi[1] = psi[0] + VOLUME/2;

  /*initialize the local (block) parts of the spinor*/
  for (ctr_t = 0; ctr_t < (VOLUME / LZ); ctr_t++) {
    memcpy(psi[0] + ctr_t * contig_block, in + (2 * ctr_t) * contig_block, contig_block * sizeof(spinor));
    memcpy(psi[1] + ctr_t * contig_block, in + (2 * ctr_t + 1) * contig_block, contig_block * sizeof(spinor));
  }
  for (i = 0; i < 2; i++) {/* loop over blocks */
    /* compute inner product */
    for (j = 0; j < g_N_s; j++) {/*loop over block.basis */
      inprod[j + i*g_N_s] = block_scalar_prod(block_list[i].basis[j], psi[i], vol);
    }
  }

  iter = lgcr(invvec, inprod, 10, 100, 1.e-10, 0, 2 * g_N_s, 2 * 9 * g_N_s, &little_D);

  /* sum up */
  mul(psi[0], invvec[0], block_list[0].basis[0], vol);
  mul(psi[1], invvec[g_N_s], block_list[1].basis[0], vol);
  for(ctr_t = 1; ctr_t < g_N_s; ctr_t++) {
    assign_add_mul(psi[0], block_list[0].basis[ctr_t], invvec[ctr_t], vol);
    assign_add_mul(psi[1], block_list[1].basis[ctr_t], invvec[g_N_s+ctr_t], vol);
  }

  reconstruct_global_field(out, psi[0], psi[1]);

  free(*psi);
  free(psi);
  free(invvec);
  free(inprod);
  return;
}


void project2(spinor * const out, spinor * const in) {
  int i,j,ctr_t;
  spinor **psi;
  int contig_block = LZ / 2;
  int vol = block_list[0].volume;
  complex *inprod;

  psi = calloc(4, sizeof(spinor*)); /*block local version of global spinor */
  inprod = calloc(2 * 9 * g_N_s, sizeof(complex)); /*inner product of spinors with bases */

  /* no loop below because further down we also don't take this cleanly into account */
  psi[0] = calloc(VOLUME, sizeof(spinor));
  psi[1] = psi[0] + VOLUME/2;

  /*initialize the local (block) parts of the spinor*/
  for (ctr_t = 0; ctr_t < (VOLUME / LZ); ctr_t++) {
    memcpy(psi[0] + ctr_t*contig_block, in + (2 * ctr_t) * contig_block, contig_block * sizeof(spinor));
    memcpy(psi[1] + ctr_t*contig_block, in + (2 * ctr_t + 1) * contig_block, contig_block * sizeof(spinor));
  }
  for (i = 0; i < 2; i++) {/* loop over blocks */
    /* compute inner product */
    for (j = 0; j < g_N_s; j++) {/*loop over block.basis */
      inprod[j + i*g_N_s] = block_scalar_prod(block_list[i].basis[j], psi[i], vol);
    }
  }

  /* sum up */
  mul(psi[0], inprod[0], block_list[0].basis[0], vol);
  mul(psi[1], inprod[g_N_s], block_list[1].basis[0], vol);
  for(ctr_t = 1; ctr_t < g_N_s; ctr_t++) {
    assign_add_mul(psi[0], block_list[0].basis[ctr_t], inprod[ctr_t], vol);
    assign_add_mul(psi[1], block_list[1].basis[ctr_t], inprod[g_N_s+ctr_t], vol);
  }

  /* reconstruct global field */
  for (ctr_t = 0; ctr_t < (VOLUME / LZ); ctr_t++) {
    memcpy(out + (2 * ctr_t) * contig_block, 
           psi[0] + ctr_t * contig_block,
	   contig_block * sizeof(spinor));
    memcpy(out + (2 * ctr_t + 1) * contig_block, 
           psi[1] + ctr_t * contig_block, 
	   contig_block * sizeof(spinor));
  }

  free(psi[0]);
  free(psi);
  free(inprod);
  return;
}

void project_left(spinor * const out, spinor * const in) {
  /* out = P_L in = in - D proj in */ 

  project(out, in);
  D_psi(g_spinor_field[DUM_MATRIX], out);
  diff(out, in, g_spinor_field[DUM_MATRIX], VOLUME);
  return;
}

void project_right(spinor * const out, spinor * const in) {
  /* out = P_R in = in - proj D in */

  D_psi(out, in);
  project(g_spinor_field[DUM_MATRIX], out);
  diff(out, in, g_spinor_field[DUM_MATRIX], VOLUME);
  return;
}

void project_left_D(spinor * const out, spinor * const in) {
  /* out = P_L D in  = in - D proj D in*/

  D_psi(g_spinor_field[DUM_SOLVER+3], in);
  project_left(out, g_spinor_field[DUM_SOLVER+3]);
  return;
}

void D_project_right(spinor * const out, spinor * const in) {
  project_right(g_spinor_field[DUM_SOLVER+3], in);
  D_psi(out, g_spinor_field[DUM_SOLVER+3]);
  return;
}

int check_projectors() {
  double nrm = 0.;

  random_spinor_field(g_spinor_field[DUM_SOLVER], VOLUME, 1);

  project_left_D(g_spinor_field[DUM_SOLVER+1], g_spinor_field[DUM_SOLVER]);
  D_project_right(g_spinor_field[DUM_SOLVER+2], g_spinor_field[DUM_SOLVER]);
  diff(g_spinor_field[DUM_SOLVER+3], g_spinor_field[DUM_SOLVER+2], g_spinor_field[DUM_SOLVER+1], VOLUME);
  nrm = square_norm(g_spinor_field[DUM_SOLVER+3], VOLUME);
  if(g_proc_id == 0) {
    printf("||P_L D psi - D P_R psi|| = %1.5e\n", sqrt(nrm));
  }


  project_left(g_spinor_field[DUM_SOLVER+1], g_spinor_field[DUM_SOLVER]);
  project_left(g_spinor_field[DUM_SOLVER+2], g_spinor_field[DUM_SOLVER+1]);
  diff(g_spinor_field[DUM_SOLVER+3], g_spinor_field[DUM_SOLVER+2], g_spinor_field[DUM_SOLVER+1], VOLUME);
  nrm = square_norm(g_spinor_field[DUM_SOLVER+3], VOLUME);
  if(g_proc_id == 0) {
    printf("||P_L^2 psi - P_L psi|| = %1.5e\n", sqrt(nrm));
  }

  project_right(g_spinor_field[DUM_SOLVER+1], g_spinor_field[DUM_SOLVER]);
  project_right(g_spinor_field[DUM_SOLVER+2], g_spinor_field[DUM_SOLVER+1]);
  diff(g_spinor_field[DUM_SOLVER+3], g_spinor_field[DUM_SOLVER+2], g_spinor_field[DUM_SOLVER+1], VOLUME);
  nrm = square_norm(g_spinor_field[DUM_SOLVER+3], VOLUME);
  if(g_proc_id == 0) {
    printf("||P_R^2 psi - P_R psi|| = %1.5e\n", sqrt(nrm));
  }

  project_left(g_spinor_field[DUM_SOLVER+1], g_spinor_field[DUM_SOLVER]);
  project2(g_spinor_field[DUM_SOLVER+2], g_spinor_field[DUM_SOLVER+1]);
  nrm = square_norm(g_spinor_field[DUM_SOLVER+2], VOLUME);
  if(g_proc_id == 0) {
    printf("||P P_L psi|| = %1.5e\n", sqrt(nrm));
  }

  project2(g_spinor_field[DUM_SOLVER+1], g_spinor_field[DUM_SOLVER]);
  project_right(g_spinor_field[DUM_SOLVER+2], g_spinor_field[DUM_SOLVER+1]);
  nrm = square_norm(g_spinor_field[DUM_SOLVER+2], VOLUME);
  if(g_proc_id == 0) {
    printf("||P_R P psi|| = %1.5e\n", sqrt(nrm));
  }

  project2(g_spinor_field[DUM_SOLVER+1], g_spinor_field[DUM_SOLVER]);
  project(g_spinor_field[DUM_SOLVER+2], g_spinor_field[DUM_SOLVER+1]);
  D_psi(g_spinor_field[DUM_SOLVER+3], g_spinor_field[DUM_SOLVER+2]);
  project2(g_spinor_field[DUM_SOLVER+2], g_spinor_field[DUM_SOLVER+3]);
  diff(g_spinor_field[DUM_SOLVER+3], g_spinor_field[DUM_SOLVER+2], g_spinor_field[DUM_SOLVER+1], VOLUME);
  nrm = square_norm(g_spinor_field[DUM_SOLVER+3], VOLUME);
  if(g_proc_id == 0) {
    printf("||P D A^-1 psi - P psi|| = %1.5e\n", sqrt(nrm));
  }

  return(0);
}

void check_little_D_inversion() {
  int i,j,ctr_t;
  spinor **psi;
  int contig_block = LZ / 2;
  int vol = block_list[0].volume;
  complex *inprod;
  complex *result;
  complex *invvec;
  double dif;

  random_spinor_field(g_spinor_field[DUM_SOLVER], VOLUME, 1);
  psi = calloc(4, sizeof(spinor*)); /*block local version of global spinor */
  result = calloc(2 * 9 * g_N_s, sizeof(complex)); /*inner product of spinors with bases */
  inprod = calloc(2 * 9 * g_N_s, sizeof(complex)); /*inner product of spinors with bases */
  invvec = calloc(2 * 9 * g_N_s, sizeof(complex)); /*inner product of spinors with bases */

  /* no loop below because further down we also don't take this cleanly into account */
  psi[0] = calloc(VOLUME, sizeof(spinor));
  psi[1] = psi[0] + VOLUME/2;

  /*initialize the local (block) parts of the spinor*/
  for (ctr_t = 0; ctr_t < (VOLUME / LZ); ++ctr_t) {
    memcpy(psi[0] + ctr_t * contig_block, g_spinor_field[DUM_SOLVER] + (2 * ctr_t) * contig_block, contig_block * sizeof(spinor));
    memcpy(psi[1] + ctr_t * contig_block, g_spinor_field[DUM_SOLVER] + (2 * ctr_t + 1) * contig_block, contig_block * sizeof(spinor));
  }
  for (i = 0; i < 2; ++i) {/* loop over blocks */
    /* compute inner product */
    for (j = 0; j < g_N_s; ++j) {/*loop over block.basis */
      inprod[j + i*g_N_s] = block_scalar_prod(block_list[i].basis[j], psi[i], vol);
    }
  }

  lgcr(invvec, inprod, 10, 100, 1.e-12, 0, 2 * g_N_s, 2 * 9 * g_N_s, &little_D);
  little_D(result, invvec); /* This should be a proper inverse now */

  dif = 0.0;
  for(ctr_t = 0; ctr_t < 2 * g_N_s; ++ctr_t){
    dif += (inprod[ctr_t].re - result[ctr_t].re) * (inprod[ctr_t].re - result[ctr_t].re);
    dif += (inprod[ctr_t].im - result[ctr_t].im) * (inprod[ctr_t].im - result[ctr_t].im);
  }
  dif = sqrt(dif);

  if (dif > 1e-8 * VOLUME){
    printf("[WARNING] check_little_D_inversion: deviation found of size %1.5e!\n", dif);
  }
  MPI_Barrier(MPI_COMM_WORLD);

  if ((g_debug_level > -1) && !g_proc_id){
    printf("Inversion check on little_D\nStart:\n");
    for(ctr_t = 0; ctr_t < 2 * g_N_s; ++ctr_t){
      printf("%1.5e + %1.5e I   ", inprod[ctr_t].re, inprod[ctr_t].im);
      if (ctr_t == g_N_s - 1)
        printf("\n");
    }
    printf("\nInverted:\n");
    for(ctr_t = 0; ctr_t < 2 * g_N_s; ++ctr_t){
      printf("%1.5e + %1.5e I   ", invvec[ctr_t].re, invvec[ctr_t].im);
      if (ctr_t == g_N_s - 1 )
        printf("\n");
    }
    printf("\nResult:\n");
    for(ctr_t = 0; ctr_t < 2 * g_N_s; ++ctr_t){
      printf("%1.5e + %1.5e I   ", result[ctr_t].re, result[ctr_t].im);
      if (ctr_t == g_N_s - 1)
        printf("\n");
    }
  }


  free(psi[0]);
  free(psi);
  free(result);
  free(invvec);
  free(inprod);
  return;
}

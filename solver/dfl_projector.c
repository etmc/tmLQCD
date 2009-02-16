/***********************************************************************
 * $Id$ 
 *
 * Copyright (C) 2008 Alber Deuzeman, Siebren Recker, Carsten Urbach
 *
 * This file is part of tmLQCD.
 *
 * tmLQCD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * tmLQCD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with tmLQCD.  If not, see <http://www.gnu.org/licenses/>.
 ***********************************************************************/

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
#include "gcr4complex.h"
#include "generate_dfl_subspace.h"
#include "dfl_projector.h"

double dfl_little_D_prec = 1.e-24;
int dfl_sloppy_prec = 0;
int init_dfl_projector = 0;
spinor **psi;
complex *inprod;
complex *invvec;
complex *work_block;
const int dfl_work_size = 13;
complex *work[13];

static void alloc_dfl_projector();

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
  int vol = block_list[0].volume;
  complex * v, * w;
  double prec;

  if(init_dfl_projector == 0) {
    alloc_dfl_projector();
  }
  v = work[0];
  w = work[1];
  /*initialize the local (block) parts of the spinor*/
  split_global_field(psi[0],psi[1], in);

  for (j = 0; j < g_N_s; j++) {/*loop over block.basis */
    inprod[j]         = scalar_prod(block_list[0].basis[j], psi[0], vol, 0);
    inprod[j + g_N_s] = scalar_prod(block_list[1].basis[j], psi[1], vol, 0);
  }

  if(dfl_sloppy_prec) prec = dfl_little_D_prec;
  else prec = 1.e-24;

  if(1) {
    iter = gcr4complex(invvec, inprod, 10, 1000, prec, 1, 2 * g_N_s, 1, 2 * 9 * g_N_s, &little_D);
    if(g_proc_id == 0 && g_debug_level > -1) {
      printf("lgcr number of iterations %d (no P_L)\n", iter);
    }
  }
  else {
    little_P_L(v, inprod);
    iter = gcr4complex(w, v, 10, 1000, prec, 1, 2 * g_N_s, 1, 2 * 9 * g_N_s, &little_P_L_D);
    little_P_R(v, w);
    little_project(w, inprod, g_N_s);
    for(i = 0; i < 2*g_N_s; i++) {
      invvec[i].re = w[i].re + v[i].re;
      invvec[i].im = w[i].im + v[i].im;
    }
    if(g_proc_id == 0 && g_debug_level > -1) {
      printf("lgcr number of iterations %d (using P_L)\n", iter);
    }
  }

  /* sum up */
  mul(psi[0], invvec[0], block_list[0].basis[0], vol);
  mul(psi[1], invvec[g_N_s], block_list[1].basis[0], vol);
  for(ctr_t = 1; ctr_t < g_N_s; ctr_t++) {
    assign_add_mul(psi[0], block_list[0].basis[ctr_t], invvec[ctr_t], vol);
    assign_add_mul(psi[1], block_list[1].basis[ctr_t], invvec[g_N_s + ctr_t], vol);
  }

  reconstruct_global_field(out, psi[0], psi[1]);
  return;
}

static void alloc_dfl_projector() {
  int i;
  
  psi = calloc(4, sizeof(spinor*)); /*block local version of global spinor */
  inprod = calloc(2 * 9 * g_N_s, sizeof(complex)); /*inner product of spinors with bases */
  invvec = calloc(2 * 9 * g_N_s, sizeof(complex)); /*inner product of spinors with bases */
  work_block = calloc(dfl_work_size * 2 * 9 * g_N_s, sizeof(complex));
  for(i = 0; i < dfl_work_size; ++i){
    work[i] = work_block + i * 2 * 9 * g_N_s;
  }
  
  /* no loop below because further down we also don't take this cleanly into account */
  psi[0] = calloc(2*(block_list[0].volume + block_list[0].spinpad), sizeof(spinor));
  psi[1] = psi[0] + (block_list[0].volume + block_list[0].spinpad);
  init_dfl_projector = 1;
  return;
}


void free_dfl_projector() {
  free(*psi);
  free(psi);
  free(invvec);
  free(inprod);
  free(work_block);
  init_dfl_projector = 0;
  return;
}

/* this is phi_k (phi_k, in) */
void project2(spinor * const out, spinor * const in) {
  int j;
  int vol = block_list[0].volume;

  if(init_dfl_projector == 0) {
    alloc_dfl_projector();
  }
  /*initialize the local (block) parts of the spinor*/
  split_global_field(psi[0],psi[1], in);

  /* compute inner product */
  for (j = 0; j < g_N_s; j++) {/*loop over block.basis */
    inprod[j]         = scalar_prod(block_list[0].basis[j], psi[0], vol, 0);
    inprod[j + g_N_s] = scalar_prod(block_list[1].basis[j], psi[1], vol, 0);
  }

  /* sum up */
  mul(psi[0], inprod[0], block_list[0].basis[0], vol);
  mul(psi[1], inprod[g_N_s], block_list[1].basis[0], vol);
  for(j = 1; j < g_N_s; j++) {
    assign_add_mul(psi[0], block_list[0].basis[j], inprod[j], vol);
    assign_add_mul(psi[1], block_list[1].basis[j], inprod[g_N_s + j], vol);
  }

  /* reconstruct global field */
  reconstruct_global_field(out, psi[0], psi[1]);
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
  /* out = P_L D in  = D in - D proj D in*/

  D_psi(g_spinor_field[DUM_MATRIX+1], in);
  project_left(out, g_spinor_field[DUM_MATRIX+1]);
  return;
}

void D_project_right(spinor * const out, spinor * const in) {
  project_right(g_spinor_field[DUM_MATRIX+1], in);
  D_psi(out, g_spinor_field[DUM_MATRIX+1]);
  return;
}


void little_project(complex * const out, complex * const in, const int  N) {
  int i, j;
  static complex *phi;
  static complex *psi;

  if(init_dfl_projector == 0) {
    alloc_dfl_projector();
  }

  phi = work[2];
  psi = work[3];
  
  /* NOTE IS THIS REALLY NECESSARY/CORRECT? */
  for(i = 0; i < N; i++) {
    phi[i] = lscalar_prod(little_dfl_fields[i], in, 2*g_N_s, 0);
  }

#ifdef MPI
  MPI_Allreduce(phi, psi, g_N_s, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
#else
  memcpy(psi, phi, g_N_s*sizeof(complex));
#endif
  
  /* apply inverse of little_A */
  for(i = 0; i < g_N_s; i++) {
    _complex_zero(phi[i]);
    for(j = 0; j < g_N_s; j++) {
      _add_assign_complex(phi[i], little_A[j*N + i], psi[j]);
    }
  }

  lmul(out, phi[0], little_dfl_fields[0], 2*g_N_s);
  for(i = 1; i < N; i++) {
    lassign_add_mul(out, little_dfl_fields[i], phi[i], 2*g_N_s);
  }
  return;
}

void little_project2(complex * const out, complex * const in, const int  N) {
  int i;
  static complex *phi;
  static complex *psi;
  
  if(init_dfl_projector == 0) {alloc_dfl_projector();}
  phi = work[4];
  psi = work[5];

  for(i = 0; i < N; i++) {
    phi[i] = lscalar_prod(little_dfl_fields[i], in, 2*g_N_s, 0);
  }
#ifdef MPI
  MPI_Allreduce(phi, psi, g_N_s, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
#else
  memcpy(psi, phi, g_N_s*sizeof(complex));
#endif
  
  lmul(out, psi[0], little_dfl_fields[0], 2*g_N_s);
  for(i = 1; i < N; i++) {
    lassign_add_mul(out, little_dfl_fields[i], psi[i], 2*g_N_s);
  }

  return;
}


void little_P_L(complex * const out, complex * const in) {
  if(init_dfl_projector == 0) {alloc_dfl_projector();}
  little_project(out, in, g_N_s);
  little_D(work[6], out);
  ldiff(out, in, work[6], 2*g_N_s);
  return;
}

void little_P_R(complex * const out, complex * const in) {
  if(init_dfl_projector == 0) {alloc_dfl_projector();}
  little_D(out, in);
  little_project(work[7], out, g_N_s);
  ldiff(out, in, work[7], 2*g_N_s);
  return;
}

void little_P_L_D(complex * const out, complex * const in) {
  if(init_dfl_projector == 0) {alloc_dfl_projector();}
  little_D(work[8], in);
  little_P_L(out, work[8]);
  return;
}

void little_D_P_R(complex * const out, complex * const in) {
  if(init_dfl_projector == 0) {alloc_dfl_projector();}
  little_P_R(work[9], in);
  little_D(out, work[9]);
  return;
}


int check_projectors() {
  double nrm = 0.;
  int j;
  spinor *phi[2];
  complex *v;

  random_spinor_field(g_spinor_field[DUM_SOLVER], VOLUME, 1);
  nrm = square_norm(g_spinor_field[DUM_SOLVER], VOLUME, 1);
  if(g_cart_id == 0) {
    printf("\nNow we check the DFL projection routines!\n\n");
    printf("||psi|| = %1.5e\n", sqrt(nrm));
  }

  split_global_field(g_spinor_field[DUM_SOLVER+1], g_spinor_field[DUM_SOLVER+2], g_spinor_field[DUM_SOLVER]);
  reconstruct_global_field(g_spinor_field[DUM_SOLVER+3], g_spinor_field[DUM_SOLVER+1], g_spinor_field[DUM_SOLVER+2]);
  diff(g_spinor_field[DUM_SOLVER+1], g_spinor_field[DUM_SOLVER+3], g_spinor_field[DUM_SOLVER], VOLUME);
  nrm = square_norm(g_spinor_field[DUM_SOLVER+1], VOLUME, 1);
  if(g_cart_id == 0) {
    printf("||psi_orig - psi_recon|| = %1.5e\n", sqrt(nrm));
    fflush(stdout);
  }

  project2(g_spinor_field[DUM_SOLVER+1], g_spinor_field[DUM_SOLVER]);
  project2(g_spinor_field[DUM_SOLVER+2], g_spinor_field[DUM_SOLVER+1]);
  diff(g_spinor_field[DUM_SOLVER+3], g_spinor_field[DUM_SOLVER+1], g_spinor_field[DUM_SOLVER+2], VOLUME);
  nrm = square_norm(g_spinor_field[DUM_SOLVER+3], VOLUME, 1);
  if(g_cart_id == 0) {
    printf("||P psi - P P psi|| = %1.5e\n", sqrt(nrm));
    fflush(stdout);
  }

  project_left_D(g_spinor_field[DUM_SOLVER+1], g_spinor_field[DUM_SOLVER]);
  D_project_right(g_spinor_field[DUM_SOLVER+2], g_spinor_field[DUM_SOLVER]);
  diff(g_spinor_field[DUM_SOLVER+3], g_spinor_field[DUM_SOLVER+2], g_spinor_field[DUM_SOLVER+1], VOLUME);
  nrm = square_norm(g_spinor_field[DUM_SOLVER+3], VOLUME, 1);
  if(g_cart_id == 0) {
    printf("||P_L D psi - D P_R psi|| = %1.5e\n", sqrt(nrm));
    fflush(stdout);
  }

  project_left(g_spinor_field[DUM_SOLVER+1], g_spinor_field[DUM_SOLVER]);
  project_left(g_spinor_field[DUM_SOLVER+2], g_spinor_field[DUM_SOLVER+1]);
  diff(g_spinor_field[DUM_SOLVER+3], g_spinor_field[DUM_SOLVER+2], g_spinor_field[DUM_SOLVER+1], VOLUME);
  nrm = square_norm(g_spinor_field[DUM_SOLVER+3], VOLUME, 1);
  if(g_cart_id == 0) {
    printf("||P_L^2 psi - P_L psi|| = %1.5e\n", sqrt(nrm));
    fflush(stdout);
  }

  project_right(g_spinor_field[DUM_SOLVER+1], g_spinor_field[DUM_SOLVER]);
  project_right(g_spinor_field[DUM_SOLVER+2], g_spinor_field[DUM_SOLVER+1]);
  diff(g_spinor_field[DUM_SOLVER+3], g_spinor_field[DUM_SOLVER+2], g_spinor_field[DUM_SOLVER+1], VOLUME);
  nrm = square_norm(g_spinor_field[DUM_SOLVER+3], VOLUME, 1);
  if(g_cart_id == 0) {
    printf("||P_R^2 psi - P_R psi|| = %1.5e\n", sqrt(nrm));
    fflush(stdout);
  }

  project_left(g_spinor_field[DUM_SOLVER+1], g_spinor_field[DUM_SOLVER]);
  project2(g_spinor_field[DUM_SOLVER+2], g_spinor_field[DUM_SOLVER+1]);
  nrm = square_norm(g_spinor_field[DUM_SOLVER+2], VOLUME, 1);
  if(g_cart_id == 0) {
    printf("||P P_L psi|| = %1.5e\n", sqrt(nrm));
    fflush(stdout);
  }

  project2(g_spinor_field[DUM_SOLVER+1], g_spinor_field[DUM_SOLVER]);
  project_right(g_spinor_field[DUM_SOLVER+2], g_spinor_field[DUM_SOLVER+1]);
  nrm = square_norm(g_spinor_field[DUM_SOLVER+2], VOLUME, 1);
  if(g_cart_id == 0) {
    printf("||P_R P psi|| = %1.5e\n", sqrt(nrm));
    fflush(stdout);
  }

  project2(g_spinor_field[DUM_SOLVER+1], g_spinor_field[DUM_SOLVER]);
  project(g_spinor_field[DUM_SOLVER+2], g_spinor_field[DUM_SOLVER+1]);
  D_psi(g_spinor_field[DUM_SOLVER+3], g_spinor_field[DUM_SOLVER+2]);
  project2(g_spinor_field[DUM_SOLVER+2], g_spinor_field[DUM_SOLVER+3]);
  diff(g_spinor_field[DUM_SOLVER+3], g_spinor_field[DUM_SOLVER+2], g_spinor_field[DUM_SOLVER+1], VOLUME);
  nrm = square_norm(g_spinor_field[DUM_SOLVER+3], VOLUME, 1);
  if(g_cart_id == 0) {
    printf("||P D A^-1 P psi - P psi|| = %1.5e\n", sqrt(nrm));
    fflush(stdout);
  }

  project2(g_spinor_field[DUM_SOLVER+1], g_spinor_field[DUM_SOLVER]);
  D_psi(g_spinor_field[DUM_SOLVER+2], g_spinor_field[DUM_SOLVER+1]);
  project(g_spinor_field[DUM_SOLVER+3], g_spinor_field[DUM_SOLVER+2]);
  project2(g_spinor_field[DUM_SOLVER+2], g_spinor_field[DUM_SOLVER+3]);
  diff(g_spinor_field[DUM_SOLVER+3], g_spinor_field[DUM_SOLVER+2], g_spinor_field[DUM_SOLVER+1], VOLUME);
  nrm = square_norm(g_spinor_field[DUM_SOLVER+3], VOLUME, 1);
  if(g_cart_id == 0) {
    printf("||P A^-1 D P psi - P psi|| = %1.5e\n", sqrt(nrm));
    fflush(stdout);
  }

  invert_little_D_spinor(g_spinor_field[DUM_SOLVER+1], g_spinor_field[DUM_SOLVER]);
  project2(g_spinor_field[DUM_SOLVER+2], g_spinor_field[DUM_SOLVER+1]);
  D_psi(g_spinor_field[DUM_SOLVER+3], g_spinor_field[DUM_SOLVER+2]);
  project2(g_spinor_field[DUM_SOLVER+2], g_spinor_field[DUM_SOLVER+3]);
  project2(g_spinor_field[DUM_SOLVER+1], g_spinor_field[DUM_SOLVER]);
  diff(g_spinor_field[DUM_SOLVER+3], g_spinor_field[DUM_SOLVER+2], g_spinor_field[DUM_SOLVER+1], VOLUME);
  nrm = square_norm(g_spinor_field[DUM_SOLVER+3], VOLUME, 1);
  if(g_cart_id == 0) {
    printf("||P D P (P D P)^-1 psi - P psi|| = %1.5e\n", sqrt(nrm));
    fflush(stdout);
  }

  invert_little_D_spinor(g_spinor_field[DUM_SOLVER+1], g_spinor_field[DUM_SOLVER]);
  apply_little_D_spinor(g_spinor_field[DUM_SOLVER+2], g_spinor_field[DUM_SOLVER+1]);
  project2(g_spinor_field[DUM_SOLVER+3], g_spinor_field[DUM_SOLVER]);
  diff(g_spinor_field[DUM_SOLVER+1], g_spinor_field[DUM_SOLVER+3], g_spinor_field[DUM_SOLVER+2], VOLUME);
  nrm = square_norm(g_spinor_field[DUM_SOLVER+1], VOLUME, 1);
  if(g_cart_id == 0) {
    printf("||A A^-1 psi - P psi|| = %1.5e\n", sqrt(nrm));
    fflush(stdout);
  }

  invert_little_D_spinor(g_spinor_field[DUM_SOLVER+1], g_spinor_field[DUM_SOLVER]);
  apply_little_D_spinor(g_spinor_field[DUM_SOLVER+2], g_spinor_field[DUM_SOLVER+1]);
  project2(g_spinor_field[DUM_SOLVER+3], g_spinor_field[DUM_SOLVER]);
  project2(g_spinor_field[DUM_SOLVER+1], g_spinor_field[DUM_SOLVER+2]);
  diff(g_spinor_field[DUM_SOLVER+2], g_spinor_field[DUM_SOLVER+3], g_spinor_field[DUM_SOLVER+1], VOLUME);
  nrm = square_norm(g_spinor_field[DUM_SOLVER+2], VOLUME, 1);
  if(g_cart_id == 0) {
    printf("||P A A^-1 psi - P psi|| = %1.5e\n", sqrt(nrm));
    fflush(stdout);
  }

  /* Different flavours for kappa != 0. First project to only a single block */
  for (j = 0; j < (VOLUME * sizeof(spinor) / sizeof(complex)); ++j){
    _complex_zero(((complex*)g_spinor_field[DUM_SOLVER+1])[j]);
    _complex_zero(((complex*)g_spinor_field[DUM_SOLVER+2])[j]);
  }

  if (!g_cart_id){
    reconstruct_global_field(g_spinor_field[DUM_SOLVER+1], block_list[0].basis[0], g_spinor_field[DUM_SOLVER+2]);
  }
  apply_little_D_spinor(g_spinor_field[DUM_SOLVER+3], g_spinor_field[DUM_SOLVER+1]);
  D_psi(g_spinor_field[DUM_SOLVER+2], g_spinor_field[DUM_SOLVER+1]);
  v = calloc(2 * 9 * g_N_s, sizeof(complex));
  phi[0] = calloc(VOLUME + 2, sizeof(spinor));
  phi[1] = phi[0] + VOLUME / 2 + 1;
  split_global_field(phi[0], phi[1], g_spinor_field[DUM_SOLVER+2]);
  if (g_cart_id == 0 && g_debug_level > 4){
    for (j = 0; j < g_N_s; ++j) {
/*       v[j]         = block_scalar_prod(phi[0], block_list[0].basis[j], VOLUME/2); */
/*       v[j + g_N_s] = block_scalar_prod(phi[1], block_list[1].basis[j], VOLUME/2); */
      v[j]         = scalar_prod(block_list[0].basis[j], phi[0], VOLUME/2, 0);
      v[j + g_N_s] = scalar_prod(block_list[1].basis[j], phi[1], VOLUME/2, 0);
    }
    for (j = 0; j < 2* g_N_s; ++j) {
      printf("AFTER D: w[%u] = %1.5e + %1.5e i\n", j, v[j].re, v[j].im);
    }
  }

  project2(g_spinor_field[DUM_SOLVER+1], g_spinor_field[DUM_SOLVER+2]);
  free(v);
  free(phi[0]);

  diff(g_spinor_field[DUM_SOLVER+2], g_spinor_field[DUM_SOLVER+3], g_spinor_field[DUM_SOLVER+1], VOLUME);
  nrm = square_norm(g_spinor_field[DUM_SOLVER+2], VOLUME, 1);
  if(g_proc_id == 0) {
    printf("||(P D - A) phi_i || = %1.5e\n", sqrt(nrm));
    fflush(stdout);
  }

  reconstruct_global_field(g_spinor_field[DUM_SOLVER+1], block_list[0].basis[0], block_list[1].basis[0]);
  apply_little_D_spinor(g_spinor_field[DUM_SOLVER+3], g_spinor_field[DUM_SOLVER+1]);
  D_psi(g_spinor_field[DUM_SOLVER+2], g_spinor_field[DUM_SOLVER+1]);
  v = calloc(2 * 9 * g_N_s, sizeof(complex));
  phi[0] = calloc(VOLUME + 2, sizeof(spinor));
  phi[1] = phi[0] + VOLUME / 2 + 1;
  split_global_field(phi[0], phi[1], g_spinor_field[DUM_SOLVER+2]);
  if (!g_proc_id && g_debug_level > 4){
    for (j = 0; j < g_N_s; ++j) {
/*       v[j]         = block_scalar_prod(phi[0], block_list[0].basis[j], VOLUME/2); */
/*       v[j + g_N_s] = block_scalar_prod(phi[1], block_list[1].basis[j], VOLUME/2); */
      v[j]         = scalar_prod(block_list[0].basis[j], phi[0], VOLUME/2, 0);
      v[j + g_N_s] = scalar_prod(block_list[1].basis[j], phi[1], VOLUME/2, 0);
    }
    for (j = 0; j < 2* g_N_s; ++j) {
      printf("AFTER D: w[%u] = %1.5e + %1.5e i\n", j, v[j].re, v[j].im);
    }
  }
  project2(g_spinor_field[DUM_SOLVER+1], g_spinor_field[DUM_SOLVER+2]);
  split_global_field(phi[0], phi[1], g_spinor_field[DUM_SOLVER+1]);
  free(v);
  free(phi[0]);


  diff(g_spinor_field[DUM_SOLVER+2], g_spinor_field[DUM_SOLVER+3], g_spinor_field[DUM_SOLVER+1], VOLUME);
  nrm = square_norm(g_spinor_field[DUM_SOLVER+2], VOLUME, 1);
  if(g_proc_id == 0) {
    printf("||(P D - A) phi || = %1.5e\n", sqrt(nrm));
    fflush(stdout);
  }

  apply_little_D_spinor(g_spinor_field[DUM_SOLVER+3], g_spinor_field[DUM_SOLVER]);
  project2(g_spinor_field[DUM_SOLVER+1], g_spinor_field[DUM_SOLVER]);
  D_psi(g_spinor_field[DUM_SOLVER+2], g_spinor_field[DUM_SOLVER+1]);
  project2(g_spinor_field[DUM_SOLVER+1], g_spinor_field[DUM_SOLVER+2]);
  diff(g_spinor_field[DUM_SOLVER+2], g_spinor_field[DUM_SOLVER+3], g_spinor_field[DUM_SOLVER+1], VOLUME);
  nrm = square_norm(g_spinor_field[DUM_SOLVER+2], VOLUME, 1);
  if(g_proc_id == 0 && g_debug_level > 4) {
    printf("||P D P psi - A psi|| = %1.5e\n", sqrt(nrm));
    printf("\n*** Comparison of the leading spinor components ***\n");
    printf("%1.5e\t%1.5e\n", g_spinor_field[DUM_SOLVER+1]->s0.c0.re, g_spinor_field[DUM_SOLVER+3]->s0.c0.re);
    printf("%1.5e\t%1.5e\n", g_spinor_field[DUM_SOLVER+1]->s0.c0.im, g_spinor_field[DUM_SOLVER+3]->s0.c0.im);
    printf("%1.5e\t%1.5e\n", g_spinor_field[DUM_SOLVER+1]->s0.c1.re, g_spinor_field[DUM_SOLVER+3]->s0.c1.re);
    printf("%1.5e\t%1.5e\n", g_spinor_field[DUM_SOLVER+1]->s0.c1.im, g_spinor_field[DUM_SOLVER+3]->s0.c1.im);
    printf("%1.5e\t%1.5e\n", g_spinor_field[DUM_SOLVER+1]->s0.c2.re, g_spinor_field[DUM_SOLVER+3]->s0.c2.re);
    printf("%1.5e\t%1.5e\n", g_spinor_field[DUM_SOLVER+1]->s0.c2.im, g_spinor_field[DUM_SOLVER+3]->s0.c2.im);
    printf("%1.5e\t%1.5e\n", g_spinor_field[DUM_SOLVER+1]->s1.c0.re, g_spinor_field[DUM_SOLVER+3]->s1.c0.re);
    printf("%1.5e\t%1.5e\n", g_spinor_field[DUM_SOLVER+1]->s1.c0.im, g_spinor_field[DUM_SOLVER+3]->s1.c0.im);
    printf("%1.5e\t%1.5e\n", g_spinor_field[DUM_SOLVER+1]->s1.c1.re, g_spinor_field[DUM_SOLVER+3]->s1.c1.re);
    printf("%1.5e\t%1.5e\n", g_spinor_field[DUM_SOLVER+1]->s1.c1.im, g_spinor_field[DUM_SOLVER+3]->s1.c1.im);
    printf("%1.5e\t%1.5e\n", g_spinor_field[DUM_SOLVER+1]->s1.c2.re, g_spinor_field[DUM_SOLVER+3]->s1.c2.re);
    printf("%1.5e\t%1.5e\n", g_spinor_field[DUM_SOLVER+1]->s1.c2.im, g_spinor_field[DUM_SOLVER+3]->s1.c2.im);
    printf("%1.5e\t%1.5e\n", g_spinor_field[DUM_SOLVER+1]->s2.c0.re, g_spinor_field[DUM_SOLVER+3]->s2.c0.re);
    printf("%1.5e\t%1.5e\n", g_spinor_field[DUM_SOLVER+1]->s2.c0.im, g_spinor_field[DUM_SOLVER+3]->s2.c0.im);
    printf("%1.5e\t%1.5e\n", g_spinor_field[DUM_SOLVER+1]->s2.c1.re, g_spinor_field[DUM_SOLVER+3]->s2.c1.re);
    printf("%1.5e\t%1.5e\n", g_spinor_field[DUM_SOLVER+1]->s2.c1.im, g_spinor_field[DUM_SOLVER+3]->s2.c1.im);
    printf("%1.5e\t%1.5e\n", g_spinor_field[DUM_SOLVER+1]->s2.c2.re, g_spinor_field[DUM_SOLVER+3]->s2.c2.re);
    printf("%1.5e\t%1.5e\n", g_spinor_field[DUM_SOLVER+1]->s2.c2.im, g_spinor_field[DUM_SOLVER+3]->s2.c2.im);
    printf("*** End of dump ***\n\n");
    fflush(stdout);
  }

  /* check little projectors now */
  if(g_cart_id == 0) {
    printf("\nNow the little little projection routines\n\n");
  }
  if(init_dfl_projector == 0) {alloc_dfl_projector();}
  
  memcpy(work[10], g_spinor_field[DUM_SOLVER], 2*g_N_s*sizeof(complex));
  little_project2(work[11], work[10], g_N_s);
  little_project2(work[12], work[11], g_N_s);
  ldiff(work[12], work[12], work[11], 2*g_N_s);
  nrm = lsquare_norm(work[12], 2*g_N_s, 1);
  if(g_cart_id == 0) {
    printf("||lP v - lP lP v|| = %1.5e\n", sqrt(nrm));
    fflush(stdout);
  }

  little_P_L_D(work[11], work[10]);
  little_P_L_D(work[12], work[10]);
  ldiff(work[12], work[12], work[11], 2*g_N_s);
  nrm = lsquare_norm(work[12], 2*g_N_s, 1);
  if(g_cart_id == 0) {
    printf("||lP_L lD v - lP_L lD v|| = %1.5e\n", sqrt(nrm));
    fflush(stdout);
  }  
  
  little_P_L_D(work[11], work[10]);
  little_D_P_R(work[12], work[10]);
  ldiff(work[12], work[12], work[11], 2*g_N_s);
  nrm = lsquare_norm(work[12], 2*g_N_s, 1);
  if(g_cart_id == 0) {
    printf("||lP_L lD v - lD lP_R v|| = %1.5e\n", sqrt(nrm));
    fflush(stdout);
  }

  little_P_R(work[11], work[10]);
  little_P_R(work[12], work[11]);
  ldiff(work[12], work[12], work[11], 2*g_N_s);
  nrm = lsquare_norm(work[12], 2*g_N_s, 1);
  if(g_cart_id == 0) {
    printf("||lP_R^2 v - lP_R v|| = %1.5e\n", sqrt(nrm));
    fflush(stdout);
  }

  little_P_L(work[11], work[10]);
  little_P_L(work[12], work[11]);
  ldiff(work[12], work[12], work[11], 2*g_N_s);
  nrm = lsquare_norm(work[12], 2*g_N_s, 1);
  if(g_cart_id == 0) {
    printf("||lP_L^2 v - lP_L v|| = %1.5e\n", sqrt(nrm));
    fflush(stdout);
  }
  return(0);
}

void check_little_D_inversion() {
  int i,j,ctr_t;
  int contig_block = LZ / 2;
  int vol = block_list[0].volume;
  complex *result, *v, *w;
  double dif;

  random_spinor_field(g_spinor_field[DUM_SOLVER], VOLUME, 1);
  if(init_dfl_projector == 0) {
    alloc_dfl_projector();
  }
  v = work[11];
  w = work[12];

  result = calloc(2 * 9 * g_N_s, sizeof(complex)); /*inner product of spinors with bases */

  /* no loop below because further down we also don't take this cleanly into account */

  /*initialize the local (block) parts of the spinor*/
  for (ctr_t = 0; ctr_t < (VOLUME / LZ); ++ctr_t) {
    memcpy(psi[0] + ctr_t * contig_block, g_spinor_field[DUM_SOLVER] + (2 * ctr_t) * contig_block, contig_block * sizeof(spinor));
    memcpy(psi[1] + ctr_t * contig_block, g_spinor_field[DUM_SOLVER] + (2 * ctr_t + 1) * contig_block, contig_block * sizeof(spinor));
  }
  for (i = 0; i < 2; ++i) {/* loop over blocks */
    /* compute inner product */
    for (j = 0; j < g_N_s; ++j) {/*loop over block.basis */
/*       inprod[j + i*g_N_s] = block_scalar_prod(block_list[i].basis[j], psi[i], vol); */
      inprod[j + i*g_N_s] = scalar_prod(psi[i], block_list[i].basis[j], vol, 0);
    }
  }

  if(1) {
    gcr4complex(invvec, inprod, 10, 1000, 1.e-31, 0, 2 * g_N_s, 1, 2 * 9 * g_N_s, &little_D);
  }
  else {
    little_P_L(v, inprod);
    gcr4complex(w, v, 10, 1000, 1.e-31, 1, 2 * g_N_s, 1, 2 * 9 * g_N_s, &little_P_L_D);
    little_P_R(v, w);
    little_project(w, inprod, g_N_s);
    for(i = 0; i < 2*g_N_s; i++) {
      invvec[i].re = w[i].re + v[i].re;
      invvec[i].im = w[i].im + v[i].im;
    }
  }
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
#ifdef MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  if ((g_debug_level > 2) && !g_proc_id){
    printf("Inversion check on little_D\nStart:\n");
    for(ctr_t = 0; ctr_t < 2 * g_N_s; ++ctr_t){
      printf("%1.9e + %1.9e I   ", inprod[ctr_t].re, inprod[ctr_t].im);
      if (ctr_t == g_N_s - 1)
        printf("\n");
    }
    printf("\nInverted:\n");
    for(ctr_t = 0; ctr_t < 2 * g_N_s; ++ctr_t){
      printf("%1.9e + %19e I   ", invvec[ctr_t].re, invvec[ctr_t].im);
      if (ctr_t == g_N_s - 1 )
        printf("\n");
    }
    printf("\nResult:\n");
    for(ctr_t = 0; ctr_t < 2 * g_N_s; ++ctr_t){
      printf("%1.9e + %1.9e I   ", result[ctr_t].re, result[ctr_t].im);
      if (ctr_t == g_N_s - 1)
        printf("\n");
    }
    printf("\n");
  }


  free(result);
  return;
}

void check_local_D() /* Should work for kappa = 0 */
{
  int j;
  double nrm;

  for (j = 0; j < (VOLUME * sizeof(spinor) / sizeof(complex)); ++j){
    _complex_zero(((complex*)g_spinor_field[DUM_SOLVER])[j]);
  }

  if (!g_proc_id){
    reconstruct_global_field(g_spinor_field[DUM_SOLVER], block_list[0].basis[0], block_list[1].basis[0]);
  }

  Block_D_psi(block_list, g_spinor_field[DUM_SOLVER + 1], block_list[0].basis[0]);
  Block_D_psi(block_list + 1, g_spinor_field[DUM_SOLVER + 2], block_list[1].basis[0]);
  D_psi(g_spinor_field[DUM_SOLVER + 3], g_spinor_field[DUM_SOLVER]);
  reconstruct_global_field(g_spinor_field[DUM_SOLVER], g_spinor_field[DUM_SOLVER + 1], g_spinor_field[DUM_SOLVER + 2]);
  diff(g_spinor_field[DUM_SOLVER + 1], g_spinor_field[DUM_SOLVER + 3], g_spinor_field[DUM_SOLVER], VOLUME);
  nrm = square_norm(g_spinor_field[DUM_SOLVER + 1], VOLUME, 1);
  if(g_proc_id == 0) {
    printf("\n\nCheck local D (trust for kappa=0 only): %1.5e\n", sqrt(nrm));
  }
}



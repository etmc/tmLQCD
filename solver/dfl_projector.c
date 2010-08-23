/***********************************************************************
 * $Id$ 
 *
 * Copyright (C) 2008 Alber Deuzeman, Siebren Reker, Carsten Urbach
 *               2010 Claude Tadonki
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
#include "Hopping_Matrix.h"
#include "little_D.h"
#include "block.h"
#include "linalg_eo.h"
#include "gcr4complex.h"
#include "generate_dfl_subspace.h"
#include "tm_operators.h"
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
  int i,j, iter;
  int vol = block_list[0].volume;
  complex * v, * w;
  double prec;

  if(init_dfl_projector == 0) {
    alloc_dfl_projector();
  }
  v = work[0];
  w = work[1]; 
  /*initialize the local (block) parts of the spinor*/
  split_global_field_GEN(psi, in, nb_blocks);
  for (j = 0; j < g_N_s; j++) {/*loop over block.basis */
    for(i = 0; i < nb_blocks; i++){
      inprod[j + i*g_N_s]  = scalar_prod(block_list[i].basis[j], psi[i], vol, 0);
    }
  }
  /* if(dfl_sloppy_prec) prec = dfl_little_D_prec; */
  if(dfl_sloppy_prec) prec = 1.e-12;
  else prec = 1.e-24;
  if(0) {
    iter = gcr4complex(invvec, inprod, 10, 50000, prec, 1, nb_blocks * g_N_s, 1, nb_blocks * 9 * g_N_s, &little_D);
    if(g_proc_id == 0 && g_debug_level > 0) {/*CT: was "g_debug_level > -1" */
      printf("lgcr number of iterations %d (no P_L)\n", iter);
    }
  }
  else {
    little_P_L(v, inprod);
    iter = gcr4complex(w, v, 20, 1000, prec, 1, nb_blocks * g_N_s, 1, nb_blocks * 9 * g_N_s, &little_P_L_D);
    little_P_R(v, w);
    little_project(w, inprod, g_N_s);
    for(i = 0; i < nb_blocks*g_N_s; i++) {
      invvec[i].re = w[i].re + v[i].re;
      invvec[i].im = w[i].im + v[i].im;
    }
    if(g_proc_id == 0 && g_debug_level > 0) {/*CT: was "g_debug_level > -1" */
      printf("lgcr number of iterations %d (using P_L)\n", iter);
    }
  }
  /* sum up */
  for(i = 0 ; i < nb_blocks ; i++) mul(psi[i], invvec[i*g_N_s], block_list[i].basis[0], vol);
  for(j = 1; j < g_N_s; j++) {
    for(i = 0 ; i < nb_blocks ; i++) assign_add_mul(psi[i], block_list[i].basis[j], invvec[i*g_N_s + j], vol);
  }

  /* reconstruct global field */
  reconstruct_global_field_GEN(out, psi, nb_blocks);
  return;
}

static void alloc_dfl_projector() {
  int i;
  
  psi = calloc(2*nb_blocks, sizeof(spinor*)); /*block local version of global spinor */
  inprod = calloc(nb_blocks * 9 * g_N_s, sizeof(complex)); /*inner product of spinors with bases */
  invvec = calloc(nb_blocks * 9 * g_N_s, sizeof(complex)); /*inner product of spinors with bases */
  work_block = calloc(dfl_work_size * nb_blocks * 9 * g_N_s, sizeof(complex));
  for(i = 0; i < dfl_work_size; ++i){
    work[i] = work_block + i * nb_blocks * 9 * g_N_s;
  }
  
  /* no loop below because further down we also don't take this cleanly into account */
  psi[0] = calloc(nb_blocks*(block_list[0].volume + block_list[0].spinpad), sizeof(spinor));
  for(i = 1 ; i < nb_blocks ; i++) psi[i] = psi[i-1] + (block_list[0].volume + block_list[0].spinpad);
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
  int i, j; 
  int vol = block_list[0].volume;

  if(init_dfl_projector == 0) {
    alloc_dfl_projector();
  }
  /*initialize the local (block) parts of the spinor*/
  split_global_field_GEN(psi, in, nb_blocks);

  /* compute inner product */
  for (j = 0; j < g_N_s; j++) {
    /*loop over block.basis */
    for(i = 0 ; i < nb_blocks ; i++)  inprod[j + i*g_N_s]  = scalar_prod(block_list[i].basis[j], psi[i], vol, 0);
  }
  
  /* sum up */
  for(i = 0 ; i < nb_blocks ; i++) mul(psi[i], inprod[i*g_N_s], block_list[i].basis[0], vol);
  for(j = 1; j < g_N_s; j++) {
    for(i = 0 ; i < nb_blocks ; i++) assign_add_mul(psi[i], block_list[i].basis[j], inprod[i*g_N_s + j], vol);
  }

  /* reconstruct global field */
  reconstruct_global_field_GEN(out, psi, nb_blocks);
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


/* out = |phi_k> A^{-1}_kl <phi_l|in> */
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
    phi[i] = lscalar_prod(little_dfl_fields[i], in, nb_blocks*N, 0);
  }

#ifdef MPI
  MPI_Allreduce(phi, psi, N, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
#else
  memcpy(psi, phi, N*sizeof(complex));
#endif
  
  /* apply inverse of little_A */
  for(i = 0; i < N; i++) {
    _mult_assign_complex(phi[i], little_A[i], psi[0]);
  }
  for(j = 1; j < N; j++) {
    for(i = 0; i < N; i++) {
      _add_assign_complex(phi[i], little_A[j*N + i], psi[j]);
    }
  }
  /* for(i = 0; i < N; i++) { */
  /*   _complex_zero(phi[i]); */
  /*   for(j = 0; j < N; j++) { */
  /*     _add_assign_complex(phi[i], little_A[j*N + i], psi[j]); */
  /*   } */
  /* } */

  lmul(out, phi[0], little_dfl_fields[0], nb_blocks*N);
  for(i = 1; i < N; i++) {
    lassign_add_mul(out, little_dfl_fields[i], phi[i], nb_blocks*N);
  }
  return;
}

/* out = |phi_k> delta_kl <phi_l|in> */
void little_project2(complex * const out, complex * const in, const int  N) {
  int i;
  static complex *phi;
  static complex *psi;
  
  if(init_dfl_projector == 0) {
    alloc_dfl_projector();
  }
  phi = work[4];
  psi = work[5];
  /* |phi> = <little_fields|in> */ 
  for(i = 0; i < N; i++) {
    phi[i] = lscalar_prod(little_dfl_fields[i], in, nb_blocks*g_N_s, 0);
  }
#ifdef MPI
  MPI_Allreduce(phi, psi, g_N_s, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
#else
  memcpy(psi, phi, g_N_s*sizeof(complex));
#endif
  
  lmul(out, psi[0], little_dfl_fields[0], nb_blocks*g_N_s);
  for(i = 1; i < N; i++) {
    lassign_add_mul(out, little_dfl_fields[i], psi[i], nb_blocks*g_N_s);
  }

  return;
}


void little_P_L(complex * const out, complex * const in) {
  if(init_dfl_projector == 0) {alloc_dfl_projector();}
  little_project(out, in, g_N_s);
  little_D(work[6], out);
  ldiff(out, in, work[6], nb_blocks*g_N_s);
  return;
}

void little_P_R(complex * const out, complex * const in) {
  if(init_dfl_projector == 0) {alloc_dfl_projector();}
  little_D(out, in);
  little_project(work[7], out, g_N_s);
  ldiff(out, in, work[7], nb_blocks*g_N_s);
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
  int i,j;
  spinor **phi;
  spinor **wphi;
  complex *v;
  phi = malloc(nb_blocks*sizeof(spinor *));
  wphi = malloc(nb_blocks*sizeof(spinor *));

  random_spinor_field(g_spinor_field[DUM_SOLVER], VOLUME, 1);
  nrm = square_norm(g_spinor_field[DUM_SOLVER], VOLUME, 1);
  if(g_cart_id == 0) {
    printf("\nNow we check the DFL projection routines!\n\n");
    printf("||psi|| = %1.5e\n", sqrt(nrm));
  }



  /* Check generalized split/reconstruct */
  phi[0] = calloc(VOLUME + nb_blocks, sizeof(spinor));
  for(j = 1; j < nb_blocks; j++) {
    phi[j] = phi[j-1] + (VOLUME/nb_blocks + 1);
  }
  split_global_field_GEN(phi, g_spinor_field[DUM_SOLVER],nb_blocks);
  reconstruct_global_field_GEN(g_spinor_field[DUM_SOLVER+1],phi,nb_blocks);
  diff(g_spinor_field[DUM_SOLVER+2], g_spinor_field[DUM_SOLVER], g_spinor_field[DUM_SOLVER+1], VOLUME);
  nrm = square_norm(g_spinor_field[DUM_SOLVER+2], VOLUME, 1);
  if(g_cart_id == 0) {
    printf("||psi_orig - psi_recon|| = %1.5e\n", sqrt(nrm));
    fflush(stdout);
  }
  /* Check even/odd split reconstruct   */
  assign(g_spinor_field[DUM_SOLVER+3], g_spinor_field[DUM_SOLVER], VOLUME);
  copy_global_to_block_eo(g_spinor_field[DUM_SOLVER+1], g_spinor_field[DUM_SOLVER+2], g_spinor_field[DUM_SOLVER], 0);
  copy_block_eo_to_global(g_spinor_field[DUM_SOLVER+3], g_spinor_field[DUM_SOLVER+1], g_spinor_field[DUM_SOLVER+2], 0);
  diff(g_spinor_field[DUM_SOLVER+2], g_spinor_field[DUM_SOLVER], g_spinor_field[DUM_SOLVER+3], VOLUME);
  nrm = square_norm(g_spinor_field[DUM_SOLVER+2], VOLUME, 1);
  if(g_cart_id == 0) {
    printf("even/odd split: ||psi_orig - psi_recon|| = %1.5e\n", sqrt(nrm));
    fflush(stdout);
  }

  project2(g_spinor_field[DUM_SOLVER+1], g_spinor_field[DUM_SOLVER]);
  project2(g_spinor_field[DUM_SOLVER+2], g_spinor_field[DUM_SOLVER+1]);
  diff(g_spinor_field[DUM_SOLVER+3], g_spinor_field[DUM_SOLVER+1], g_spinor_field[DUM_SOLVER+2], VOLUME);
  nrm = square_norm(g_spinor_field[DUM_SOLVER+3], VOLUME, 1);
  if(g_cart_id == 0) {
    printf("||P2 psi - P2 P2 psi|| = %1.5e\n", sqrt(nrm));
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
    wphi[0] = block_list[0].basis[0];
    for(i = 1; i< nb_blocks; i++) wphi[i] = g_spinor_field[DUM_SOLVER+2];
    reconstruct_global_field_GEN(g_spinor_field[DUM_SOLVER+1], wphi, nb_blocks);
  }
  apply_little_D_spinor(g_spinor_field[DUM_SOLVER+3], g_spinor_field[DUM_SOLVER+1]);
  D_psi(g_spinor_field[DUM_SOLVER+2], g_spinor_field[DUM_SOLVER+1]);
  
  if (g_cart_id == 0 && g_debug_level > 4){
    v = calloc(nb_blocks * 9 * g_N_s, sizeof(complex));
    split_global_field_GEN(phi, g_spinor_field[DUM_SOLVER+2], nb_blocks);
    
    for (j = 0; j < g_N_s; ++j) 
      for(i = 0; i < nb_blocks; i++)
	v[j + i*g_N_s] = scalar_prod(block_list[i].basis[j], phi[i], VOLUME/nb_blocks, 0);
    
    for (j = 0; j < nb_blocks* g_N_s; ++j) {
      printf("AFTER D: w[%u] = %1.5e + %1.5e i\n", j, v[j].re, v[j].im);
    }
    free(v);
  }
  
  project2(g_spinor_field[DUM_SOLVER+1], g_spinor_field[DUM_SOLVER+2]);
  
  
  diff(g_spinor_field[DUM_SOLVER+2], g_spinor_field[DUM_SOLVER+3], g_spinor_field[DUM_SOLVER+1], VOLUME);
  nrm = square_norm(g_spinor_field[DUM_SOLVER+2], VOLUME, 1);
  if(g_proc_id == 0) {
    printf("||(P D - A) phi_i || = %1.5e\n", sqrt(nrm));
    fflush(stdout);
  }
  
  reconstruct_global_field_GEN_ID(g_spinor_field[DUM_SOLVER+1], block_list, 0, nb_blocks);
  apply_little_D_spinor(g_spinor_field[DUM_SOLVER+3], g_spinor_field[DUM_SOLVER+1]);
  D_psi(g_spinor_field[DUM_SOLVER+2], g_spinor_field[DUM_SOLVER+1]);
  if (!g_proc_id && g_debug_level > 4){
    v = calloc(nb_blocks * 9 * g_N_s, sizeof(complex));
    split_global_field_GEN(phi, g_spinor_field[DUM_SOLVER+2],nb_blocks);
    for (j = 0; j < g_N_s; ++j) 
      for(i = 0; i < nb_blocks; i++)
	v[j + i*g_N_s] = scalar_prod(block_list[i].basis[j], phi[i], VOLUME/nb_blocks, 0);
    for (j = 0; j < nb_blocks* g_N_s; ++j) {
      printf("AFTER D: w[%u] = %1.5e + %1.5e i\n", j, v[j].re, v[j].im);
    }
    free(v);
  }
  project2(g_spinor_field[DUM_SOLVER+1], g_spinor_field[DUM_SOLVER+2]);
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
  if(init_dfl_projector == 0) {
    alloc_dfl_projector();
  }
  
  memcpy(work[10], g_spinor_field[DUM_SOLVER], nb_blocks*g_N_s*sizeof(complex));
  little_project2(work[11], work[10], g_N_s);
  little_project2(work[12], work[11], g_N_s);
  ldiff(work[12], work[12], work[11], nb_blocks*g_N_s);
  nrm = lsquare_norm(work[12], nb_blocks*g_N_s, 1);
  if(g_cart_id == 0) {
    printf("||lP2 v - lP2 lP2 v|| = %1.5e\n", sqrt(nrm));
    fflush(stdout);
  }

  little_P_L_D(work[11], work[10]);
  little_P_L_D(work[12], work[10]);
  ldiff(work[12], work[12], work[11], nb_blocks*g_N_s);
  nrm = lsquare_norm(work[12], nb_blocks*g_N_s, 1);
  if(g_cart_id == 0) {
    printf("||lP_L lD v - lP_L lD v|| = %1.5e\n", sqrt(nrm));
    fflush(stdout);
  }  
  
  little_P_L_D(work[11], work[10]);
  little_D_P_R(work[12], work[10]);
  ldiff(work[12], work[12], work[11], nb_blocks*g_N_s);
  nrm = lsquare_norm(work[12], nb_blocks*g_N_s, 1);
  if(g_cart_id == 0) {
    printf("||lP_L lD v - lD lP_R v|| = %1.5e\n", sqrt(nrm));
    fflush(stdout);
  }

  little_P_R(work[11], work[10]);
  little_P_R(work[12], work[11]);
  ldiff(work[12], work[12], work[11], nb_blocks*g_N_s);
  nrm = lsquare_norm(work[12], nb_blocks*g_N_s, 1);
  if(g_cart_id == 0) {
    printf("||lP_R^2 v - lP_R v|| = %1.5e\n", sqrt(nrm));
    fflush(stdout);
  }

  little_P_L(work[11], work[10]);
  little_P_L(work[12], work[11]);
  ldiff(work[12], work[12], work[11], nb_blocks*g_N_s);
  nrm = lsquare_norm(work[12], nb_blocks*g_N_s, 1);
  if(g_cart_id == 0) {
    printf("||lP_L^2 v - lP_L v|| = %1.5e\n", sqrt(nrm));
    fflush(stdout);
  }

  free(phi[0]);
  free(phi);
  free(wphi);

  return(0);
}

void check_little_D_inversion() {
  int i,j,ctr_t;
  int contig_block = LZ / nb_blocks;
  int vol = block_list[0].volume;
  complex *result, *v, *w;
  double dif;

  random_spinor_field(g_spinor_field[DUM_SOLVER], VOLUME, 1);
  if(init_dfl_projector == 0) {
    alloc_dfl_projector();
  }
  v = work[11];
  w = work[12];

  result = calloc(nb_blocks * 9 * g_N_s, sizeof(complex)); /*inner product of spinors with bases */

  /* no loop below because further down we also don't take this cleanly into account */

  /*initialize the local (block) parts of the spinor*/
  for (ctr_t = 0; ctr_t < (VOLUME / LZ); ++ctr_t) {
    for(i=0; i< nb_blocks; i++)
      memcpy(psi[i] + ctr_t * contig_block, g_spinor_field[DUM_SOLVER] + (nb_blocks * ctr_t + i) * contig_block, contig_block * sizeof(spinor));
  }
  for (i = 0; i < nb_blocks; ++i) {/* loop over blocks */
    /* compute inner product */
    for (j = 0; j < g_N_s; ++j) {/*loop over block.basis */
      /*       inprod[j + i*g_N_s] = block_scalar_prod(block_list[i].basis[j], psi[i], vol); */
      inprod[j + i*g_N_s] = scalar_prod(psi[i], block_list[i].basis[j], vol, 0);
    }
  }

  if(1) {
    gcr4complex(invvec, inprod, 10, 1000, 1.e-31, 0, nb_blocks * g_N_s, 1, nb_blocks * 9 * g_N_s, &little_D);
  }
  else {
    little_P_L(v, inprod);
    gcr4complex(w, v, 10, 1000, 1.e-31, 1, nb_blocks * g_N_s, 1, nb_blocks * 9 * g_N_s, &little_P_L_D);
    little_P_R(v, w);
    little_project(w, inprod, g_N_s);
    for(i = 0; i < nb_blocks*g_N_s; i++) {
      invvec[i].re = w[i].re + v[i].re;
      invvec[i].im = w[i].im + v[i].im;
    }
  }
  little_D(result, invvec); /* This should be a proper inverse now */

  dif = 0.0;
  for(ctr_t = 0; ctr_t < nb_blocks * g_N_s; ++ctr_t){
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
    for(ctr_t = 0; ctr_t < nb_blocks * g_N_s; ++ctr_t){
      printf("%1.9e + %1.9e I   ", inprod[ctr_t].re, inprod[ctr_t].im);
      if (ctr_t == g_N_s - 1)
        printf("\n");
    }
    printf("\nInverted:\n");
    for(ctr_t = 0; ctr_t < nb_blocks * g_N_s; ++ctr_t){
      printf("%1.9e + %19e I   ", invvec[ctr_t].re, invvec[ctr_t].im);
      if (ctr_t == g_N_s - 1 )
        printf("\n");
    }
    printf("\nResult:\n");
    for(ctr_t = 0; ctr_t < nb_blocks * g_N_s; ++ctr_t){
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
  int j, vol = block_list[0].volume/2, i;
  double nrm;
  block_convert_lexic_to_eo(g_spinor_field[DUM_SOLVER], g_spinor_field[DUM_SOLVER+1], block_list[0].basis[0]);
  block_convert_eo_to_lexic(g_spinor_field[DUM_SOLVER+2], g_spinor_field[DUM_SOLVER], g_spinor_field[DUM_SOLVER+1]);
  diff(g_spinor_field[DUM_SOLVER], g_spinor_field[DUM_SOLVER+2], block_list[0].basis[0], block_list[0].volume);
  nrm = square_norm(g_spinor_field[DUM_SOLVER], block_list[0].volume, 0);
  if(g_proc_id == 0) {
    printf("\nblock even/odd: ||psi - psi_recon|| = %1.5e\n", sqrt(nrm));
  }

  for(j = 0; j < nb_blocks; j++) {
    zero_spinor_field(g_spinor_field[DUM_SOLVER], VOLUME);
    Block_D_psi(&block_list[j], g_spinor_field[DUM_SOLVER+6], block_list[j].basis[0]);

    /* Now test the block hopping matrix */
    /* split into even/odd sites         */
    block_convert_lexic_to_eo(g_spinor_field[DUM_SOLVER], g_spinor_field[DUM_SOLVER+1], block_list[j].basis[0]);
  
    /* Even sites */
    Block_H_psi(&block_list[j], g_spinor_field[DUM_DERI], g_spinor_field[DUM_SOLVER+1], EO);
    assign_mul_one_pm_imu(g_spinor_field[DUM_SOLVER+2], g_spinor_field[DUM_SOLVER], 1., vol); 
    assign_add_mul_r(g_spinor_field[DUM_SOLVER+2], g_spinor_field[DUM_DERI], 1., vol);

    /* Odd sites */
    Block_H_psi(&block_list[j], g_spinor_field[DUM_DERI], g_spinor_field[DUM_SOLVER], OE);
    assign_mul_one_pm_imu(g_spinor_field[DUM_SOLVER+3], g_spinor_field[DUM_SOLVER+1], 1., vol); 
    assign_add_mul_r(g_spinor_field[DUM_SOLVER+3], g_spinor_field[DUM_DERI], 1., vol);

    /* convert back to block spinor */
    block_convert_eo_to_lexic(g_spinor_field[DUM_SOLVER+5], g_spinor_field[DUM_SOLVER+2], g_spinor_field[DUM_SOLVER+3]);

    if(g_proc_id == 0 && g_debug_level > 5) {
      for(i = 0; i < block_list[0].volume; i++) {
	if(fabs(g_spinor_field[DUM_SOLVER+6][i].s0.c0.re) > 1.e-15 || fabs(g_spinor_field[DUM_SOLVER+5][i].s0.c0.re) > 1.e-15) {
	  printf("%d %e %d\n", i, g_spinor_field[DUM_SOLVER+6][i].s0.c0.re, block_list[0].volume);
	  printf("%d %e\n", i, g_spinor_field[DUM_SOLVER+5][i].s0.c0.re);
	}
      }
    }

    diff(g_spinor_field[DUM_SOLVER + 4], g_spinor_field[DUM_SOLVER + 5], g_spinor_field[DUM_SOLVER+6], block_list[0].volume);
    nrm = square_norm(g_spinor_field[DUM_SOLVER + 4], block_list[0].volume, 0);
    if(g_proc_id == 0) {
      printf("Check local D against Hopping Matrix: %1.5e block %d\n", sqrt(nrm), j);
    }
  }
  return;
}



/***********************************************************************
 *
 * Copyright (C) 2008 Alber Deuzeman, Siebren Reker, Carsten Urbach
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
#include <complex.h>
#include "block.h"
#include "linalg/blas.h"
#include "operator/D_psi.h"
#include "operator/Hopping_Matrix.h"
#include "little_D.h"
#include "block.h"
#include "linalg_eo.h"
#include "gcr4complex.h"
#include "generate_dfl_subspace.h"
#include "operator/tm_operators.h"
#include "boundary.h"
#include "Msap.h"
#include "mr.h"
#include "solver_field.h"
#include "dfl_projector.h"

double dfl_little_D_prec = 1.e-24;
int dfl_sloppy_prec = 0;
int init_dfl_projector = 0;
spinor **psi;
_Complex double *inprod;
_Complex double *inprod_eo;
_Complex double *inprod_o;
_Complex double *inprod_e;
_Complex double *invvec;
_Complex double *invvec_eo;
_Complex double *ctmp;
_Complex double *work_block;
const int dfl_work_size = 16;
_Complex double *work[16];

static void alloc_dfl_projector();

/* Break up full volume spinor to blocks
 * loop over block.basis
 * compute inner product and store as _Complex double vector
 * compute A^-1 * _Complex double vector
 * loop over block.basis
 * compute sum of basis vectors times _Complex double element
 * create global vector */

/* this is phi_k A^{-1}_{kl} (phi_k, in) */
void project(spinor * const out, spinor * const in) {
  int i,j, i_e, i_o, iter;
  int evenodd = 0;
  int usePL = 0;
  int vol = block_list[0].volume;
  _Complex double * v, * w;
  double prec;
  
  if(init_dfl_projector == 0) {
    alloc_dfl_projector();
  }
  v = work[0];
  w = work[1]; 
  /*initialize the local (block) parts of the spinor*/
  split_global_field_GEN(psi, in, nb_blocks);
  
  for (j = 0; j < g_N_s*nb_blocks*9; j++) {
    (inprod[j]) = 0.0;
    (inprod_o[j]) = 0.0;
    (inprod_eo[j]) = 0.0;
    (inprod_e[j]) = 0.0;
    (invvec[j]) = 0.0;
    (invvec_eo[j]) = 0.0;
    (ctmp[j]) = 0.0;
    (w[j]) = 0.0;
    (v[j]) = 0.0;
  }
  
  for (j = 0; j < g_N_s; j++) {/*loop over block.basis */
    i_o=0;
    i_e=0;
    for(i = 0; i < nb_blocks; i++) {
      inprod[j + i*g_N_s]  = scalar_prod(block_list[i].basis[j], psi[i], vol, 0);
      if(evenodd) {
	if (block_list[i].evenodd==0) {
	  inprod_eo[j + i_e*g_N_s] = inprod[j + i*g_N_s];
	  i_e++;
	}
	if (block_list[i].evenodd==1) {
	  inprod_eo[j + nb_blocks*g_N_s/2+i_o*g_N_s] = inprod[j + i*g_N_s];
	  i_o++;
	}
      }
    }
  }
  
  if(evenodd) {
    little_D_ee_inv(inprod_e,inprod_eo);
    little_D_hop(1,inprod_o, inprod_e);
    little_Dhat_rhs(1,inprod_o,-1,inprod_eo);
  }
  
  
  /* if(dfl_sloppy_prec) prec = dfl_little_D_prec;*/
  if(dfl_sloppy_prec) prec = 1.e-12;
  else prec = 1.e-24;
  
  
  
  if(!usePL) {
    if(evenodd) {
      iter = gcr4complex(invvec_eo,inprod_o,10,1000,prec,1,nb_blocks*g_N_s,1,nb_blocks*9*g_N_s,&little_D_sym);
      
      little_D_hop(0,ctmp, invvec_eo);
      little_D_ee_inv(invvec_eo,ctmp);
      little_Dhat_rhs(0,invvec_eo, -1., inprod_e);
    
      for (j = 0; j < g_N_s; j++) {
	i_o=0;
	i_e=0;
	for(i = 0; i < nb_blocks; i++) {
	  if (block_list[i].evenodd==0) {
	    invvec[j + i*g_N_s] = invvec_eo[j + i_e*g_N_s];
	    i_e++;
	  }
	  if (block_list[i].evenodd==1) {
	    invvec[j + i*g_N_s] = invvec_eo[j + nb_blocks*g_N_s/2+i_o*g_N_s];
	    i_o++;
	  }
	}
      }
      if(g_proc_id == 0 && g_debug_level > 0) {/*CT: was "g_debug_level > -1" */
	printf("lgcr evenodd number of iterations %d (no P_L)\n", iter);
      }
    }
    else {
      iter = gcr4complex(invvec, inprod, 10, 1000, prec, 1, nb_blocks * g_N_s, 1, nb_blocks * 9 * g_N_s, &little_D);
      if(g_proc_id == 0 && g_debug_level > 0) {/*CT: was "g_debug_level > -1" */
	printf("lgcr number of iterations %d (no P_L)\n", iter);
      }
    }
  }
  else {
    if(evenodd) {
      little_P_L_sym(v, inprod_o);
      iter = gcr4complex(w, v, 10, 1000, prec, 1, nb_blocks * g_N_s, 1, nb_blocks * 9 * g_N_s, &little_P_L_D_sym);
      little_P_R_sym(v, w);
/*      little_project(w, inprod_o, g_N_s);*/
      little_project_eo(w,inprod_o,g_N_s);
      for(i = 0; i < nb_blocks*g_N_s; ++i)
	invvec_eo[i] = w[i] + v[i];
      little_D_hop(0,ctmp, invvec_eo);
      little_D_ee_inv(invvec_eo,ctmp);
      little_Dhat_rhs(0,invvec_eo, -1., inprod_e);
      for (j = 0; j < g_N_s; j++) {
	i_o=0;
	i_e=0;
	for(i = 0; i < nb_blocks; i++){
	  if (block_list[i].evenodd==0) {
	    invvec[j + i*g_N_s] = invvec_eo[j + i_e*g_N_s];
	    i_e++;
	  }
	  if (block_list[i].evenodd==1) {
	    invvec[j + i*g_N_s] = invvec_eo[j + nb_blocks*g_N_s/2+i_o*g_N_s];
	    i_o++;
	  }
	}
      } 
      if(g_proc_id == 0 && g_debug_level > 0) {/*CT: was "g_debug_level > -1" */
	printf("lgcr even/odd number of iterations %d (using P_L)\n", iter);
      }
    }
    else {
      little_P_L(v, inprod);
      iter = gcr4complex(w, v, 10, 1000, prec, 1, nb_blocks * g_N_s, 1, nb_blocks * 9 * g_N_s, &little_P_L_D);
      little_P_R(v, w);
      little_project(w, inprod, g_N_s);
      for(i = 0; i < nb_blocks*g_N_s; ++i)
	invvec[i] = w[i] + v[i];
      if(g_proc_id == 0 && g_debug_level > 0) {/*CT: was "g_debug_level > -1" */
	printf("lgcr number of iterations %d (using P_L)\n", iter);
      }
    }    
  }
  /* sum up */
  for(i = 0 ; i < nb_blocks ; i++) {
    mul(psi[i], invvec[i*g_N_s], block_list[i].basis[0], vol);
  }
  for(j = 1; j < g_N_s; j++) {
    for(i = 0 ; i < nb_blocks ; i++) {
      assign_add_mul(psi[i], block_list[i].basis[j], invvec[i*g_N_s + j], vol);
    }
  }
  
  /* reconstruct global field */
  reconstruct_global_field_GEN(out, psi, nb_blocks);
  free_dfl_projector();
  return;
}

static void alloc_dfl_projector() {
  int i;
  
  psi = calloc(2*nb_blocks, sizeof(spinor*)); /*block local version of global spinor */
  inprod = calloc(nb_blocks * 9 * g_N_s, sizeof(_Complex double)); /*inner product of spinors with bases */
  inprod_eo = calloc(nb_blocks * 9 * g_N_s, sizeof(_Complex double)); /*inner product of spinors with bases */
  inprod_o = calloc(nb_blocks * 9 * g_N_s, sizeof(_Complex double)); /*inner product of spinors with bases */
  inprod_e = calloc(nb_blocks * 9 * g_N_s, sizeof(_Complex double)); /*inner product of spinors with bases */
  ctmp = calloc(nb_blocks * 9 * g_N_s, sizeof(_Complex double)); /*inner product of spinors with bases */
  invvec = calloc(nb_blocks * 9 * g_N_s, sizeof(_Complex double)); /*inner product of spinors with bases */
  invvec_eo = calloc(nb_blocks * 9 * g_N_s, sizeof(_Complex double)); /*inner product of spinors with bases */
  work_block = calloc(dfl_work_size * nb_blocks * 9 * g_N_s, sizeof(_Complex double));
  for(i = 0; i < dfl_work_size; ++i) {
    work[i] = work_block + i * nb_blocks * 9 * g_N_s;
  }
  
  /* no loop below because further down we also don't take this cleanly into account */
  psi[0] = calloc(nb_blocks*(block_list[0].volume + block_list[0].spinpad), sizeof(spinor));
  for(i = 1 ;i < nb_blocks ;i++) {
    psi[i] = psi[i-1] + (block_list[0].volume + block_list[0].spinpad);
  }
  init_dfl_projector = 1;
  return;
}


void free_dfl_projector() {
  free(*psi);
  free(psi);
  free(invvec);
  free(invvec_eo);
  free(inprod);
  free(inprod_eo);
  free(inprod_e);
  free(inprod_o);
  free(ctmp);
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
void little_project(_Complex double * const out, _Complex double * const in, const int  N) {
  int i, j;
  static _Complex double *phi;
  static _Complex double *psi;

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
  memcpy(psi, phi, N*sizeof(_Complex double));
#endif
  
  /* apply inverse of little_A */
  for(i = 0; i < N; i++) {
    (phi[i]) = 0.0;
    for(j = 0; j < N; j++) {
      (phi[i]) += (little_A[j*N + i]) * (psi[j]);
    }
  }

  lmul(out, phi[0], little_dfl_fields[0], nb_blocks*N);
  for(i = 1; i < N; i++) {
    lassign_add_mul(out, little_dfl_fields[i], phi[i], nb_blocks*N);
  }
  return;
}

void little_project_eo(_Complex double * const out, _Complex double * const in, const int  N) {
  int i, j;
  static _Complex double *phi;
  static _Complex double *psi;
  
  if(init_dfl_projector == 0) {
    alloc_dfl_projector();
  }
      
  phi = work[2];
  psi = work[3];

  /* NOTE IS THIS REALLY NECESSARY/CORRECT? */
  for(i = 0; i < N; i++) {
    phi[i] = lscalar_prod(little_dfl_fields_eo[i], in, nb_blocks*N, 0);
  }
  
#ifdef MPI
  MPI_Allreduce(phi, psi, N, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
#else
  memcpy(psi, phi, N*sizeof(_Complex double));
#endif

  /* apply inverse of little_A_eo */
  for(i = 0; i < N; i++) {
    (phi[i]) = 0.0;
    for(j = 0; j < N; j++) {
      (phi[i]) += (little_A_eo[j*N + i]) * (psi[j]);
    }
  }
  
  lmul(out, phi[0], little_dfl_fields_eo[0], nb_blocks*N);
  for(i = 1; i < N; i++) {
    lassign_add_mul(out, little_dfl_fields_eo[i], phi[i], nb_blocks*N);
  }
  return;
}


void little_project2(_Complex double * const out, _Complex double * const in, const int  N) {
  int i;
  static _Complex double *phi;
  static _Complex double *psi;
  
  if(init_dfl_projector == 0) {alloc_dfl_projector();}
  phi = work[4];
  psi = work[5];

  for(i = 0; i < N; i++) {
    phi[i] = lscalar_prod(little_dfl_fields[i], in, nb_blocks*N, 0);
  }
#ifdef MPI
  MPI_Allreduce(phi, psi, g_N_s, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
#else
  memcpy(psi, phi, g_N_s*sizeof(_Complex double));
#endif
  
  lmul(out, psi[0], little_dfl_fields[0], nb_blocks*g_N_s);
  for(i = 1; i < N; i++) {
    lassign_add_mul(out, little_dfl_fields[i], psi[i], nb_blocks*g_N_s);
  }

  return;
}


void little_P_L(_Complex double * const out, _Complex double * const in) {
  if(init_dfl_projector == 0) {alloc_dfl_projector();}
  little_project(out, in, g_N_s);
  little_D(work[6], out);
  ldiff(out, in, work[6], nb_blocks*g_N_s);
  return;
}

void little_P_R(_Complex double * const out, _Complex double * const in) {
  if(init_dfl_projector == 0) {alloc_dfl_projector();}
  little_D(out, in);
  little_project(work[7], out, g_N_s);
  ldiff(out, in, work[7], nb_blocks*g_N_s);
  return;
}

void little_P_L_sym(_Complex double * const out, _Complex double * const in) {
  if(init_dfl_projector == 0) {alloc_dfl_projector();}
/*  little_project(out, in, g_N_s);*/
  little_project_eo(out,in,g_N_s);
  little_D_sym(work[13], out);
  ldiff(out, in, work[13], nb_blocks*g_N_s);
  return;
}

void little_P_R_sym(_Complex double * const out, _Complex double * const in) {
  if(init_dfl_projector == 0) {alloc_dfl_projector();}
  little_D_sym(out, in);
/*  little_project(work[14], out, g_N_s);*/
  little_project_eo(work[14],out,g_N_s);
  ldiff(out, in, work[14], nb_blocks*g_N_s);
  return;
}

void little_P_L_D(_Complex double * const out, _Complex double * const in) {
  if(init_dfl_projector == 0) {alloc_dfl_projector();}
  little_D(work[8], in);
  little_P_L(out, work[8]);
  return;
}

void little_P_L_D_sym(_Complex double * const out, _Complex double * const in) {
  if(init_dfl_projector == 0) {alloc_dfl_projector();}
  little_D_sym(work[15], in);
  little_P_L_sym(out, work[15]);
  return;
}

void little_D_P_R(_Complex double * const out, _Complex double * const in) {
  if(init_dfl_projector == 0) {alloc_dfl_projector();}
  little_P_R(work[9], in);
  little_D(out, work[9]);
  return;
}


int check_projectors(const int repro) {
  double nrm = 0.;
  int i,j;
  spinor **phi;
  spinor **wphi;
  _Complex double *v;
  spinor ** work_fields = NULL;
  const int nr_wf = 4;

  init_solver_field(&work_fields, VOLUMEPLUSRAND, nr_wf);
  phi = malloc(nb_blocks*sizeof(spinor *));
  wphi = malloc(nb_blocks*sizeof(spinor *));

  random_spinor_field_lexic(work_fields[0], repro, RN_GAUSS);
  nrm = square_norm(work_fields[0], VOLUME, 1);
  if(g_cart_id == 0) {
    printf("\nNow we check the DFL projection routines!\n\n");
    printf("||psi|| = %1.5e\n", sqrt(nrm));
  }



  /* Check generalized split/reconstruct */
  phi[0] = calloc(VOLUME + nb_blocks, sizeof(spinor));
  for(j = 1; j < nb_blocks; j++) {
    phi[j] = phi[j-1] + (VOLUME/nb_blocks + 1);
  }
  split_global_field_GEN(phi, work_fields[0],nb_blocks);
  reconstruct_global_field_GEN(work_fields[1],phi,nb_blocks);
  diff(work_fields[2], work_fields[0], work_fields[1], VOLUME);
  nrm = square_norm(work_fields[2], VOLUME, 1);
  if(g_cart_id == 0) {
    printf("||psi_orig - psi_recon|| = %1.5e\n", sqrt(nrm));
    fflush(stdout);
  }
  /* Check even/odd split reconstruct   */
  assign(work_fields[3], work_fields[0], VOLUME);
  copy_global_to_block_eo(work_fields[1], work_fields[2], work_fields[0], 0);
  copy_block_eo_to_global(work_fields[3], work_fields[1], work_fields[2], 0);
  diff(work_fields[2], work_fields[0], work_fields[3], VOLUME);
  nrm = square_norm(work_fields[2], VOLUME, 1);
  if(g_cart_id == 0) {
    printf("even/odd split: ||psi_orig - psi_recon|| = %1.5e\n", sqrt(nrm));
    fflush(stdout);
  }

  project2(work_fields[1], work_fields[0]);
  project2(work_fields[2], work_fields[1]);
  diff(work_fields[3], work_fields[1], work_fields[2], VOLUME);
  nrm = square_norm(work_fields[3], VOLUME, 1);
  if(g_cart_id == 0) {
    printf("||P2 psi - P2 P2 psi|| = %1.5e\n", sqrt(nrm));
    fflush(stdout);
  }

  project_left_D(work_fields[1], work_fields[0]);
  D_project_right(work_fields[2], work_fields[0]);
  diff(work_fields[3], work_fields[2], work_fields[1], VOLUME);
  nrm = square_norm(work_fields[3], VOLUME, 1);
  if(g_cart_id == 0) {
    printf("||P_L D psi - D P_R psi|| = %1.5e\n", sqrt(nrm));
    fflush(stdout);
  }

  project_left(work_fields[1], work_fields[0]);
  project_left(work_fields[2], work_fields[1]);
  diff(work_fields[3], work_fields[2], work_fields[1], VOLUME);
  nrm = square_norm(work_fields[3], VOLUME, 1);
  if(g_cart_id == 0) {
    printf("||P_L^2 psi - P_L psi|| = %1.5e\n", sqrt(nrm));
    fflush(stdout);
  }

  project_right(work_fields[1], work_fields[0]);
  project_right(work_fields[2], work_fields[1]);
  diff(work_fields[3], work_fields[2], work_fields[1], VOLUME);
  nrm = square_norm(work_fields[3], VOLUME, 1);
  if(g_cart_id == 0) {
    printf("||P_R^2 psi - P_R psi|| = %1.5e\n", sqrt(nrm));
    fflush(stdout);
  }

  project_left(work_fields[1], work_fields[0]);
  project2(work_fields[2], work_fields[1]);
  nrm = square_norm(work_fields[2], VOLUME, 1);
  if(g_cart_id == 0) {
    printf("||P P_L psi|| = %1.5e\n", sqrt(nrm));
    fflush(stdout);
  }

  project2(work_fields[1], work_fields[0]);
  project_right(work_fields[2], work_fields[1]);
  nrm = square_norm(work_fields[2], VOLUME, 1);
  if(g_cart_id == 0) {
    printf("||P_R P psi|| = %1.5e\n", sqrt(nrm));
    fflush(stdout);
  }

  project2(work_fields[1], work_fields[0]);
  project(work_fields[2], work_fields[1]);
  D_psi(work_fields[3], work_fields[2]);
  project2(work_fields[2], work_fields[3]);
  diff(work_fields[3], work_fields[2], work_fields[1], VOLUME);
  nrm = square_norm(work_fields[3], VOLUME, 1);
  if(g_cart_id == 0) {
    printf("||P D A^-1 P psi - P psi|| = %1.5e\n", sqrt(nrm));
    fflush(stdout);
  }

  project2(work_fields[1], work_fields[0]);
  D_psi(work_fields[2], work_fields[1]);
  project(work_fields[3], work_fields[2]);
  project2(work_fields[2], work_fields[3]);
  diff(work_fields[3], work_fields[2], work_fields[1], VOLUME);
  nrm = square_norm(work_fields[3], VOLUME, 1);
  if(g_cart_id == 0) {
    printf("||P A^-1 D P psi - P psi|| = %1.5e\n", sqrt(nrm));
    fflush(stdout);
  }

  invert_little_D_spinor(work_fields[1], work_fields[0]);
  project2(work_fields[2], work_fields[1]);
  D_psi(work_fields[3], work_fields[2]);
  project2(work_fields[2], work_fields[3]);
  project2(work_fields[1], work_fields[0]);
  diff(work_fields[3], work_fields[2], work_fields[1], VOLUME);
  nrm = square_norm(work_fields[3], VOLUME, 1);
  if(g_cart_id == 0) {
    printf("||P D P (P D P)^-1 psi - P psi|| = %1.5e\n", sqrt(nrm));
    fflush(stdout);
  }

  
  invert_little_D_spinor(work_fields[1], work_fields[0]);
  invert_little_D_eo_spinor(work_fields[2], work_fields[0]);
  diff(work_fields[3], work_fields[1], work_fields[2], VOLUME);
  nrm = square_norm(work_fields[3], VOLUME, 1);
  if(g_cart_id == 0) {
    printf("||A^-1 psi - A^-1_eo psi|| = %1.5e\n", sqrt(nrm));
    fflush(stdout);
  }


  invert_little_D_spinor(work_fields[1], work_fields[0]);
  apply_little_D_spinor(work_fields[2], work_fields[1]);
  project2(work_fields[3], work_fields[0]);
  diff(work_fields[1], work_fields[3], work_fields[2], VOLUME);
  nrm = square_norm(work_fields[1], VOLUME, 1);
  if(g_cart_id == 0) {
    printf("||A A^-1 psi - P psi|| = %1.5e\n", sqrt(nrm));
    fflush(stdout);
  }

  invert_little_D_spinor(work_fields[1], work_fields[0]);
  apply_little_D_spinor(work_fields[2], work_fields[1]);
  project2(work_fields[3], work_fields[0]);
  project2(work_fields[1], work_fields[2]);
  diff(work_fields[2], work_fields[3], work_fields[1], VOLUME);
  nrm = square_norm(work_fields[2], VOLUME, 1);
  if(g_cart_id == 0) {
    printf("||P A A^-1 psi - P psi|| = %1.5e\n", sqrt(nrm));
    fflush(stdout);
  }


  /* Different flavours for kappa != 0. First project to only a single block */
  for (j = 0; j < (VOLUME * sizeof(spinor) / sizeof(_Complex double)); ++j){
    ((_Complex double*)work_fields[1])[j] = 0.;
    ((_Complex double*)work_fields[2])[j] = 0.;
  }

  if (!g_cart_id){
    wphi[0] = block_list[0].basis[0];
    for(i = 1; i< nb_blocks; i++) {
      wphi[i] = work_fields[2];
    }
    reconstruct_global_field_GEN(work_fields[1], wphi, nb_blocks);
  }
  apply_little_D_spinor(work_fields[3], work_fields[1]);
  D_psi(work_fields[2], work_fields[1]);
  
  if (g_cart_id == 0 && g_debug_level > 4){
    v = calloc(nb_blocks * 9 * g_N_s, sizeof(_Complex double));
    split_global_field_GEN(phi, work_fields[2], nb_blocks);

    for (j = 0; j < g_N_s; ++j) {
      for(i = 0; i < nb_blocks; i++) {
	v[j + i*g_N_s] = scalar_prod(block_list[i].basis[j], phi[i], VOLUME/nb_blocks, 0);
      }
    }

    for (j = 0; j < nb_blocks* g_N_s; ++j) {
      printf("AFTER D: w[%u] = %1.5e + %1.5e i\n", j, creal(v[j]), cimag(v[j]));
    }
    free(v);
  }

  project2(work_fields[1], work_fields[2]);
  

  diff(work_fields[2], work_fields[3], work_fields[1], VOLUME);
  nrm = square_norm(work_fields[2], VOLUME, 1);
  if(g_proc_id == 0) {
    printf("||(P D - A) phi_i || = %1.5e\n", sqrt(nrm));
    fflush(stdout);
  }

  reconstruct_global_field_GEN_ID(work_fields[1], block_list, 0, nb_blocks);
  apply_little_D_spinor(work_fields[3], work_fields[1]);
  D_psi(work_fields[2], work_fields[1]);
  if (!g_proc_id && g_debug_level > 4){
    v = calloc(nb_blocks * 9 * g_N_s, sizeof(_Complex double));
    split_global_field_GEN(phi, work_fields[2],nb_blocks);
    for (j = 0; j < g_N_s; ++j) 
      for(i = 0; i < nb_blocks; i++)
	v[j + i*g_N_s] = scalar_prod(block_list[i].basis[j], phi[i], VOLUME/nb_blocks, 0);
    for (j = 0; j < nb_blocks* g_N_s; ++j) {
      printf("AFTER D: w[%u] = %1.5e + %1.5e i\n", j, creal(v[j]), cimag(v[j]));
    }
    free(v);
  }
  project2(work_fields[1], work_fields[2]);
  diff(work_fields[2], work_fields[3], work_fields[1], VOLUME);
  nrm = square_norm(work_fields[2], VOLUME, 1);
  if(g_proc_id == 0) {
    printf("||(P D - A) phi || = %1.5e\n", sqrt(nrm));
    fflush(stdout);
  }

  apply_little_D_spinor(work_fields[3], work_fields[0]);
  project2(work_fields[1], work_fields[0]);
  D_psi(work_fields[2], work_fields[1]);
  project2(work_fields[1], work_fields[2]);
  diff(work_fields[2], work_fields[3], work_fields[1], VOLUME);
  nrm = square_norm(work_fields[2], VOLUME, 1);
  if(g_proc_id == 0 && g_debug_level > 4) {
    printf("||P D P psi - A psi|| = %1.5e\n", sqrt(nrm));
    printf("\n*** Comparison of the leading spinor components ***\n");
    printf("%1.5e\t%1.5e\n", creal(work_fields[1]->s0.c0), creal(work_fields[3]->s0.c0));
    printf("%1.5e\t%1.5e\n", cimag(work_fields[1]->s0.c0), cimag(work_fields[3]->s0.c0));
    printf("%1.5e\t%1.5e\n", creal(work_fields[1]->s0.c1), creal(work_fields[3]->s0.c1));
    printf("%1.5e\t%1.5e\n", cimag(work_fields[1]->s0.c1), cimag(work_fields[3]->s0.c1));
    printf("%1.5e\t%1.5e\n", creal(work_fields[1]->s0.c2), creal(work_fields[3]->s0.c2));
    printf("%1.5e\t%1.5e\n", cimag(work_fields[1]->s0.c2), cimag(work_fields[3]->s0.c2));
    printf("%1.5e\t%1.5e\n", creal(work_fields[1]->s1.c0), creal(work_fields[3]->s1.c0));
    printf("%1.5e\t%1.5e\n", cimag(work_fields[1]->s1.c0), cimag(work_fields[3]->s1.c0));
    printf("%1.5e\t%1.5e\n", creal(work_fields[1]->s1.c1), creal(work_fields[3]->s1.c1));
    printf("%1.5e\t%1.5e\n", cimag(work_fields[1]->s1.c1), cimag(work_fields[3]->s1.c1));
    printf("%1.5e\t%1.5e\n", creal(work_fields[1]->s1.c2), creal(work_fields[3]->s1.c2));
    printf("%1.5e\t%1.5e\n", cimag(work_fields[1]->s1.c2), cimag(work_fields[3]->s1.c2));
    printf("%1.5e\t%1.5e\n", creal(work_fields[1]->s2.c0), creal(work_fields[3]->s2.c0));
    printf("%1.5e\t%1.5e\n", cimag(work_fields[1]->s2.c0), cimag(work_fields[3]->s2.c0));
    printf("%1.5e\t%1.5e\n", creal(work_fields[1]->s2.c1), creal(work_fields[3]->s2.c1));
    printf("%1.5e\t%1.5e\n", cimag(work_fields[1]->s2.c1), cimag(work_fields[3]->s2.c1));
    printf("%1.5e\t%1.5e\n", creal(work_fields[1]->s2.c2), creal(work_fields[3]->s2.c2));
    printf("%1.5e\t%1.5e\n", cimag(work_fields[1]->s2.c2), cimag(work_fields[3]->s2.c2));
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
  
  memcpy(work[10], work_fields[0], nb_blocks*g_N_s*sizeof(_Complex double));
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
  finalize_solver(work_fields, nr_wf);
  return(0);
}

void check_little_D_inversion(const int repro) {
  int i,j,ctr_t;
  int contig_block = LZ / nb_blocks;
  int vol = block_list[0].volume;
  _Complex double *result, *v, *w;
  double dif;
  spinor ** work_fields = NULL;
  const int nr_wf = 1;

  init_solver_field(&work_fields, VOLUMEPLUSRAND, nr_wf);
  random_spinor_field_lexic(work_fields[0], repro, RN_GAUSS);
  if(init_dfl_projector == 0) {
    alloc_dfl_projector();
  }
  v = work[11];
  w = work[12];

  result = calloc(nb_blocks * 9 * g_N_s, sizeof(_Complex double)); /*inner product of spinors with bases */

  /* no loop below because further down we also don't take this cleanly into account */

  /*initialize the local (block) parts of the spinor*/
  for (ctr_t = 0; ctr_t < (VOLUME / LZ); ++ctr_t) {
    for(i=0; i< nb_blocks; i++) {
      memcpy(psi[i] + ctr_t * contig_block, work_fields[0] + (nb_blocks * ctr_t + i) * contig_block, contig_block * sizeof(spinor));
    }
  }
  for (i = 0; i < nb_blocks; ++i) {/* loop over blocks */
    /* compute inner product */
    for (j = 0; j < g_N_s; ++j) {/*loop over block.basis */
      /*       inprod[j + i*g_N_s] = block_scalar_prod(block_list[i].basis[j], psi[i], vol); */
      inprod[j + i*g_N_s] = scalar_prod(psi[i], block_list[i].basis[j], vol, 0);
    }
  }

  if(1) {
    gcr4complex(invvec, inprod, 10, 1000, 1.e-24, 0, nb_blocks * g_N_s, 1, nb_blocks * 9 * g_N_s, &little_D);
  }
  else {
    little_P_L(v, inprod);
    gcr4complex(w, v, 10, 1000, 1.e-24, 1, nb_blocks * g_N_s, 1, nb_blocks * 9 * g_N_s, &little_P_L_D);
    little_P_R(v, w);
    little_project(w, inprod, g_N_s);
    for(i = 0; i < nb_blocks*g_N_s; ++i)
      invvec[i] = w[i] + v[i];
  }
  little_D(result, invvec); /* This should be a proper inverse now */

  dif = 0.0;
  for(ctr_t = 0; ctr_t < nb_blocks * g_N_s; ++ctr_t){
    dif += (creal(inprod[ctr_t]) - creal(result[ctr_t])) * (creal(inprod[ctr_t]) - creal(result[ctr_t]));
    dif += (cimag(inprod[ctr_t]) - cimag(result[ctr_t])) * (cimag(inprod[ctr_t]) - cimag(result[ctr_t]));
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
      printf("%1.9e + %1.9e I   ", creal(inprod[ctr_t]), cimag(inprod[ctr_t]));
      if (ctr_t == g_N_s - 1)
        printf("\n");
    }
    printf("\nInverted:\n");
    for(ctr_t = 0; ctr_t < nb_blocks * g_N_s; ++ctr_t){
      printf("%1.9e + %19e I   ", creal(invvec[ctr_t]), cimag(invvec[ctr_t]));
      if (ctr_t == g_N_s - 1 )
        printf("\n");
    }
    printf("\nResult:\n");
    for(ctr_t = 0; ctr_t < nb_blocks * g_N_s; ++ctr_t){
      printf("%1.9e + %1.9e I   ", creal(result[ctr_t]), cimag(result[ctr_t]));
      if (ctr_t == g_N_s - 1)
        printf("\n");
    }
    printf("\n");
  }

  finalize_solver(work_fields, nr_wf);
  free(result);
  return;
}

void check_local_D(const int repro)
{
  spinor * r[8];
  int j, vol = block_list[0].volume/2, i;
  double nrm;
  spinor ** work_fields = NULL;
  const int nr_wf = 7;

  init_solver_field(&work_fields, VOLUMEPLUSRAND, nr_wf);
  block_convert_lexic_to_eo(work_fields[0], work_fields[1], block_list[0].basis[0]);
  block_convert_eo_to_lexic(work_fields[2], work_fields[0], work_fields[1]);
  diff(work_fields[0], work_fields[2], block_list[0].basis[0], block_list[0].volume);
  nrm = square_norm(work_fields[0], block_list[0].volume, 0);
  if(g_proc_id == 0) {
    printf("\nblock even/odd: ||psi - psi_recon|| = %1.5e\n", sqrt(nrm));
    fflush(stdout);
  }

  for(j = 0; j < nb_blocks; j++) {
    zero_spinor_field(work_fields[0], VOLUME);
    Block_D_psi(&block_list[j], work_fields[6], block_list[j].basis[0]);

    /* Now test the block hopping matrix */
    /* split into even/odd sites         */
    block_convert_lexic_to_eo(work_fields[0], work_fields[1], block_list[j].basis[0]);
  
    /* Even sites */
    Block_H_psi(&block_list[j], g_spinor_field[DUM_DERI], work_fields[1], EO);
    assign_mul_one_pm_imu(work_fields[2], work_fields[0], 1., vol); 
    assign_add_mul_r(work_fields[2], g_spinor_field[DUM_DERI], 1., vol);

    /* Odd sites */
    Block_H_psi(&block_list[j], g_spinor_field[DUM_DERI], work_fields[0], OE);
    assign_mul_one_pm_imu(work_fields[3], work_fields[1], 1., vol); 
    assign_add_mul_r(work_fields[3], g_spinor_field[DUM_DERI], 1., vol);

    /* convert back to block spinor */
    block_convert_eo_to_lexic(work_fields[5], work_fields[2], work_fields[3]);

    if(g_proc_id == 0 && g_debug_level > 5) {
      for(i = 0; i < block_list[0].volume; i++) {
	if(fabs(creal(work_fields[6][i].s0.c0)) > 1.e-15 || fabs(creal(work_fields[5][i].s0.c0)) > 1.e-15) {
	  printf("%d %e %d\n", i, creal(work_fields[6][i].s0.c0), block_list[0].volume);
	  printf("%d %e\n", i, creal(work_fields[5][i].s0.c0));
	}
      }
    }

    diff(work_fields[4], work_fields[5], work_fields[6], block_list[0].volume);
    nrm = square_norm(work_fields[4], block_list[0].volume, 0);
    if(sqrt(nrm) > 1.e-12) {
      printf("Check failed for local D against Hopping Matrix: ||delta|| = %1.5e block %d process %d\n", sqrt(nrm), j, g_proc_id);
    }
  }
  /* check Msap and Msap_eo on a radom vector */
  random_spinor_field_lexic(work_fields[0], repro, RN_GAUSS);
  zero_spinor_field(work_fields[1], VOLUME);
  Msap(work_fields[1], work_fields[0], 2);
  D_psi(work_fields[2], work_fields[1]);
  diff(work_fields[3], work_fields[2], work_fields[0], VOLUME);
  nrm = square_norm(work_fields[3], VOLUME, 1);
  if(g_proc_id == 0) {
    printf("Msap relaxed the residue to ||r||^2 = %1.5e\n", nrm);
  }

  zero_spinor_field(work_fields[1], VOLUME);
  Msap_eo(work_fields[1], work_fields[0], 2);
  D_psi(work_fields[2], work_fields[1]);
  diff(work_fields[3], work_fields[2], work_fields[0], VOLUME);
  nrm = square_norm(work_fields[3], VOLUME, 1);
  if(g_proc_id == 0) {
    printf("Msap_eo relaxed the residue to ||r||^2 = %1.5e\n", nrm);
  }

  for(j = 0; j < 6; j++) {
    r[j] = work_fields[j];
  }
  for(j = 0; j < nb_blocks; j++) {
    
    block_convert_lexic_to_eo(r[0], r[1], block_list[j].basis[0]);
    /* check even/odd inversion for Block_D_psi*/
    /* varphi_e in r[2] */
    assign_mul_one_pm_imu_inv(r[2], r[0], +1., vol);
    Block_H_psi(&block_list[j], r[3], r[2], OE);
    /* a_odd = a_odd + b_odd */
    /* varphi_o in r[3] */
    assign_mul_add_r(r[3], -1., r[1], vol);
    /* psi_o in r[1] */
    mrblk(r[1], r[3], 3, 1.e-31, 1, vol, &Mtm_plus_block_psi, j);
    
    Block_H_psi(&block_list[j], r[0], r[1], EO);
    mul_one_pm_imu_inv(r[0], +1., vol);
    /* a_even = a_even + b_even */
    /* check this sign +1 seems to be right in Msap_eo */
    assign_add_mul_r(r[2], r[0], -1., vol);
    
    block_convert_eo_to_lexic(r[4], r[2], r[1]);
    
    Block_D_psi(&block_list[j], r[5], r[4]);
    diff(r[0], block_list[j].basis[0], r[5], block_list[j].volume);
    nrm = square_norm(r[0], block_list[j].volume, 0);
    if(g_proc_id == 0) {
      printf("mr_eo, block=%d: ||r||^2 = %1.5e\n", j, nrm);
    }
  }
  for(j = 0; j < nb_blocks; j++) {
    block_convert_lexic_to_eo(r[0], r[1], block_list[j].basis[0]);
    /* check even/odd inversion for Block_D_psi*/
    /* varphi_e in r[2] */
    assign_mul_one_pm_imu_inv(r[2], r[0], +1., vol);
    Block_H_psi(&block_list[j], r[3], r[2], OE);
    /* a_odd = a_odd + b_odd */
    /* varphi_o in r[3] */
    assign_mul_add_r(r[3], -1., r[1], vol);
    /* psi_o in r[1] */
    mul_one_pm_imu_inv(r[3], +1., vol); 
    mrblk(r[1], r[3], 3, 1.e-31, 1, vol, &Mtm_plus_sym_block_psi, j);
    
    Block_H_psi(&block_list[j], r[0], r[1], EO);
    mul_one_pm_imu_inv(r[0], +1., vol);
    /* a_even = a_even + b_even */
    /* check this sign +1 seems to be right in Msap_eo */
    assign_add_mul_r(r[2], r[0], -1., vol);
    
    block_convert_eo_to_lexic(r[4], r[2], r[1]);
    
    Block_D_psi(&block_list[j], r[5], r[4]);
    diff(r[0], block_list[j].basis[0], r[5], block_list[j].volume);
    nrm = square_norm(r[0], block_list[j].volume, 0);
    if(g_proc_id == 0) {
      printf("mr_eo (symmetric eo), block=%d: ||r||^2 = %1.5e\n", j, nrm);
    }
  }
  finalize_solver(work_fields, nr_wf);
  return;
}




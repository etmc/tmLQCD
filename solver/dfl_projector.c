/***********************************************************************
 *
 * Copyright (C) 2008 Albert Deuzeman, Siebren Reker, Carsten Urbach
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
#include "mcr4complex.h"
#include "mr4complex.h"
#include "cgne4complex.h"
#include "generate_dfl_subspace.h"
#include "operator/tm_operators.h"
#include "operator/clovertm_operators.h"
#include "boundary.h"
#include "Msap.h"
#include "mr.h"
#include "solver_field.h"
#include "solver.h"
#include "dfl_projector.h"

int dfl_sloppy_prec = 1;
int init_dfl_projector = 0;
spinor **psi;
_Complex double *inprod;
_Complex float  *inprod32;
_Complex double *inprod_eo;
_Complex double *inprod_o;
_Complex float *inprod_o32;
_Complex double *inprod_e;
_Complex double *invvec;
_Complex float  *invvec32;
_Complex double *invvec_eo;
_Complex float *invvec_eo32;
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
  int i_e, i_o, iter;
  int evenodd = 1;
  int gcr32 = 1;
  int little_m = little_gmres_m_parameter;
  int little_max_iter = little_solver_max_iter;
  int vol = block_list[0].volume;
  _Complex double * v, * w;
  double prec;
  evenodd = little_evenodd;
  if(init_dfl_projector == 0) {
    alloc_dfl_projector();
  }
  v = work[0];
  w = work[1]; 
  /*initialize the local (block) parts of the spinor*/
  split_global_field_GEN(psi, in, nb_blocks);

  for (int j = 0; j < g_N_s*nb_blocks*9; j++) {
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

  for (int j = 0; j < g_N_s; j++) {/*loop over block.basis */
    i_o = 0;
    i_e = 0;
    for(int i = 0; i < nb_blocks; i++) {
      inprod[j + i*g_N_s]  = scalar_prod(block_list[i].basis[j], psi[i], vol, 0);
      inprod32[j + i*g_N_s]  = (_Complex float)inprod[j + i*g_N_s];
      if(evenodd) {
        if (block_list[i].evenodd == 0) {
          inprod_eo[j + i_e*g_N_s] = inprod[j + i*g_N_s];
          i_e++;
        }
        if (block_list[i].evenodd == 1) {
          inprod_eo[j + nb_blocks*g_N_s/2+i_o*g_N_s] = inprod[j + i*g_N_s];
          i_o++;
        }
      }
    }
  }

  if(evenodd) {
    little_D_ee_inv(inprod_e, inprod_eo);
    little_D_hop(1, inprod_o, inprod_e);
    little_Dhat_rhs(1, inprod_o, -1, inprod_eo);
  }


  if(!dfl_sloppy_prec) prec = little_solver_high_prec;
  else prec = little_solver_low_prec;

  if(!usePL) {
    if(evenodd) {
      if(gcr32) {
        for (int j = 0; j < g_N_s*nb_blocks*9; j++) {
          inprod_o32[j] = (_Complex float) inprod_o[j];
        }
	iter = gcr4complex32(invvec_eo32, inprod_o32, little_m, little_max_iter, prec, 1, 
			     nb_blocks*g_N_s, 1, nb_blocks*9*g_N_s, 0, &little_D_sym32);
	// we could do more in 32bit precision!?
        for (int j = 0; j < g_N_s*nb_blocks*9; j++) {
          invvec_eo[j] = (_Complex double) invvec_eo32[j];
        }
      }
      else {
	iter = gcr4complex(invvec_eo, inprod_o, little_m, little_max_iter, prec, 1, 
			   nb_blocks*g_N_s, 1, nb_blocks*9*g_N_s, 0, &little_D_sym);
      }

      little_D_hop(0, ctmp, invvec_eo);
      little_D_ee_inv(invvec_eo, ctmp);
      little_Dhat_rhs(0,invvec_eo, -1., inprod_e);

      for (int j = 0; j < g_N_s; j++) {
        i_o=0;
        i_e=0;
        for(int i = 0; i < nb_blocks; i++) {
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
      if(g_proc_id == 0 && g_debug_level > 2) {
        printf("lgcr (even/odd) number of iterations %d (no LittleLittleD)\n", iter);
      }
    }
    else {
      if(gcr32) {
        iter = gcr4complex32(invvec32, inprod32, little_m, little_max_iter, prec, 1, 
                             nb_blocks * g_N_s, 1, nb_blocks * 9 * g_N_s, 0, &little_D32);
        
        for (int j = 0; j < g_N_s*nb_blocks*9; j++) {
          invvec[j] = (_Complex double) invvec32[j];
        }
      }
      else {
        iter = gcr4complex(invvec, inprod, little_m, little_max_iter, prec, 1, 
                           nb_blocks * g_N_s, 1, nb_blocks * 9 * g_N_s, 0, &little_D);
      }
      if(g_proc_id == 0 && g_debug_level > 2) {
        printf("lgcr number of iterations %d (no LittleLittleD)\n", iter);
      }       
    }
  }
  else { // usePL = true
    if(evenodd) {
      // this is in adaptive MG style
      if(gcr32) {
        for (int j = 0; j < g_N_s*nb_blocks*9; j++) {
          inprod_o32[j] = (_Complex float) inprod_o[j];
        }
	iter = gcr4complex32(invvec_eo32, inprod_o32, little_m, little_max_iter, prec, 1, 
			     nb_blocks*g_N_s, 1, nb_blocks*9*g_N_s, 0, &little_D_sym32);
	// we could do more in 32bit precision!?
        for (int j = 0; j < g_N_s*nb_blocks*9; j++) {
          invvec_eo[j] = (_Complex double) invvec_eo32[j];
        }
      }
      else {
	iter = gcr4complex(invvec_eo, inprod_o, little_m, little_max_iter, prec, 1, 
			   nb_blocks * g_N_s, 1, nb_blocks * 9 * g_N_s, 1, &little_D_sym);
      }
      little_D_hop(0,ctmp, invvec_eo);
      little_D_ee_inv(invvec_eo,ctmp);
      little_Dhat_rhs(0,invvec_eo, -1., inprod_e);
      for (int j = 0; j < g_N_s; j++) {
        i_o=0;
        i_e=0;
        for(int i = 0; i < nb_blocks; i++){
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
      if(g_proc_id == 0 && g_debug_level > 0) {
        printf("lgcr (even/odd) number of iterations %d (using LittleLittleD)\n", iter);
      }
    }
    else {
      little_P_L(v, inprod);
      iter = gcr4complex(w, v, little_m, little_max_iter, prec, 1, 
                         nb_blocks * g_N_s, 1, nb_blocks * 9 * g_N_s, 0, &little_P_L_D);
      little_P_R(v, w);
      little_project(w, inprod, g_N_s);
      for(int i = 0; i < nb_blocks*g_N_s; ++i)
        invvec[i] = w[i] + v[i];
      if(g_proc_id == 0 && g_debug_level > 0) {
        printf("lgcr number of iterations %d (using LittleLittleD)\n", iter);
      }
    }    
  }
  /* sum up */
  for(int i = 0 ; i < nb_blocks ; i++) {
    mul(psi[i], invvec[i*g_N_s], block_list[i].basis[0], vol);
  }
  for(int j = 1; j < g_N_s; j++) {
    for(int i = 0 ; i < nb_blocks ; i++) {
      assign_add_mul(psi[i], block_list[i].basis[j], invvec[i*g_N_s + j], vol);
    }
  }

  /* reconstruct global field */
  reconstruct_global_field_GEN(out, psi, nb_blocks);
  free_dfl_projector();
  return;
}

static void alloc_dfl_projector() {
  if(init_dfl_projector == 0) {
    
    psi = calloc(2*nb_blocks, sizeof(spinor*)); /*block local version of global spinor */
    inprod = calloc(nb_blocks * 9 * g_N_s, sizeof(_Complex double)); /*inner product of spinors with bases */
    inprod32 = calloc(nb_blocks * 9 * g_N_s, sizeof(_Complex float)); /*inner product of spinors with bases */
    inprod_eo = calloc(nb_blocks * 9 * g_N_s, sizeof(_Complex double)); /*inner product of spinors with bases */
    inprod_o = calloc(nb_blocks * 9 * g_N_s, sizeof(_Complex double)); /*inner product of spinors with bases */
    inprod_o32 = calloc(nb_blocks * 9 * g_N_s, sizeof(_Complex float)); /*inner product of spinors with bases */
    inprod_e = calloc(nb_blocks * 9 * g_N_s, sizeof(_Complex double)); /*inner product of spinors with bases */
    ctmp = calloc(nb_blocks * 9 * g_N_s, sizeof(_Complex double)); /*inner product of spinors with bases */
    invvec = calloc(nb_blocks * 9 * g_N_s, sizeof(_Complex double)); /*inner product of spinors with bases */
    invvec32 = calloc(nb_blocks * 9 * g_N_s, sizeof(_Complex float)); /*inner product of spinors with bases */
    invvec_eo = calloc(nb_blocks * 9 * g_N_s, sizeof(_Complex double)); /*inner product of spinors with bases */
    invvec_eo32 = calloc(nb_blocks * 9 * g_N_s, sizeof(_Complex float)); /*inner product of spinors with bases */
    work_block = calloc(dfl_work_size * nb_blocks * 9 * g_N_s, sizeof(_Complex double));
    for(int i = 0; i < dfl_work_size; ++i) {
      work[i] = work_block + i * nb_blocks * 9 * g_N_s;
    }
    
    /* no loop below because further down we also don't take this cleanly into account */
    psi[0] = calloc(nb_blocks*(block_list[0].volume + block_list[0].spinpad), sizeof(spinor));
    for(int i = 1 ;i < 2*nb_blocks ;i++) {
      psi[i] = psi[i-1] + (block_list[0].volume + block_list[0].spinpad);
    }
    init_dfl_projector = 1;
  }
  return;
}


void free_dfl_projector() {
  if(init_dfl_projector) {
    free(*psi);
    free(psi);
    free(invvec);
    free(invvec32);
    free(invvec_eo);
    free(invvec_eo32);
    free(inprod);
    free(inprod32);
    free(inprod_eo);
    free(inprod_e);
    free(inprod_o);
    free(inprod_o32);
    free(ctmp);
    free(work_block);
    init_dfl_projector = 0;
  }
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

// This is a preconditioner for D in Multi-Grid spirit
// following equation (4.2) in arXiv:1303.1377
// C^(nu) psi = M_sap[(psi - D phi) + phi]
// with approximately P A P phi = psi   (A = little D)
// and nu the M_sap cycles here called Ncy

void mg_precon(spinor * const out, spinor * const in) {
  // phi = PD_c^{-1} P^dagger in
  project(out, in);
  // in - D*phi 
  // need to DUM_MATRIX+2,3 because in Msap_eo DUM_MATRIX+0,1 is used
  D_psi(g_spinor_field[DUM_MATRIX+2], out);
  diff(g_spinor_field[DUM_MATRIX+2], in, g_spinor_field[DUM_MATRIX+2], VOLUME);
  // apply M_SAP
  zero_spinor_field(g_spinor_field[DUM_MATRIX+3], VOLUME);
  Msap_eo(g_spinor_field[DUM_MATRIX+3], g_spinor_field[DUM_MATRIX+2], NcycleMsap, NiterMsap);
  // sum with phi
  add(out, g_spinor_field[DUM_MATRIX+3], out, VOLUME);
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
  little_project_eo(out,in,g_N_s);
  little_D_sym(work[13], out);
  ldiff(out, in, work[13], nb_blocks*g_N_s);
  return;
}

void little_P_R_sym(_Complex double * const out, _Complex double * const in) {
  if(init_dfl_projector == 0) {alloc_dfl_projector();}
  little_D_sym(out, in);
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

void little_mg_precon(_Complex double * const out, _Complex double * const in) {
  // phi = PD_c^{-1} P^dagger in
  little_project_eo(out, in, g_N_s);
  // in - D*phi
  little_D_sym(work[2], out);
  ldiff(work[3], in, work[2], nb_blocks*g_N_s);
  // sum with phi
  ladd(out, work[3], out, nb_blocks*g_N_s);
  return;
}

// little_P_L_D_sym * psi = (1 - PA^-1P little_D_sym) * little_D_sym * psi
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
  spinor ** work_fields = NULL;
  const int nr_wf = 5;
  const double eps = 1.e-8;
  double savelittle_solver_high_prec = little_solver_high_prec;
  little_solver_high_prec = eps*eps/10.;

  init_solver_field(&work_fields, VOLUMEPLUSRAND, nr_wf);
  phi = malloc(nb_blocks * sizeof(spinor *));
  wphi = malloc(nb_blocks * sizeof(spinor *));

  random_spinor_field_lexic(work_fields[0], repro, RN_GAUSS);
  nrm = square_norm(work_fields[0], VOLUME, 1);
  if(g_cart_id == 0) {
    printf("\n######################\n");
    printf("# Now we check the DFL projection routines!\n\n");
    printf("# ||psi|| = %1.5e\n", sqrt(nrm));
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
    printf("# ||psi_orig - psi_recon|| = %1.5e ", sqrt(nrm));
    if(sqrt(nrm) < eps) printf("#  -> passed\n\n");
    else printf("#  -> FAILED!\n\n");
    fflush(stdout);
  }
  /* Check even/odd split reconstruct   */
  assign(work_fields[3], work_fields[0], VOLUME);
  copy_global_to_block_eo(work_fields[1], work_fields[2], work_fields[0], 0);
  copy_block_eo_to_global(work_fields[3], work_fields[1], work_fields[2], 0);
  diff(work_fields[2], work_fields[0], work_fields[3], VOLUME);
  nrm = square_norm(work_fields[2], VOLUME, 1);
  if(g_cart_id == 0) {
    printf("# even/odd split: ||psi_orig - psi_recon|| = %1.5e ", sqrt(nrm));
    if(sqrt(nrm) < eps) printf("#  -> passed\n\n");
    else printf("#  -> FAILED!\n\n");
    fflush(stdout);
  }

  // check assign_mul_one_sw_pm_imu_inv_block for clover case
  if(g_c_sw > 0) {
    for (int blk = 0; blk < nb_blocks; blk++) {
      copy_global_to_block_eo(phi[0], phi[1], work_fields[0], blk);
      assign_mul_one_sw_pm_imu_inv_block(EE, phi[2], phi[0], g_mu, &block_list[blk]);
      
      copy_block_eo_to_global(work_fields[1], phi[2], phi[1], blk);      
    }
    convert_lexic_to_eo(work_fields[2], work_fields[3], work_fields[0]);
    assign_mul_one_sw_pm_imu_inv(EE, work_fields[5], work_fields[2], g_mu);
    convert_eo_to_lexic(work_fields[2], work_fields[5], work_fields[3]);
    diff(work_fields[0], work_fields[1], work_fields[2], VOLUME);
    nrm = square_norm(work_fields[0], VOLUME, 1);
    if(g_cart_id == 0) {
      printf("# assign_mul_one_sw_pm_imu_inv: ||psi_orig - psi_block|| = %1.5e ", sqrt(nrm));
      if(sqrt(nrm) < eps) printf("#  -> passed\n\n");
      else printf("#  -> FAILED!\n\n");
      fflush(stdout);
    }
  }
  
  project2(work_fields[1], work_fields[0]);
  project2(work_fields[2], work_fields[1]);
  diff(work_fields[3], work_fields[1], work_fields[2], VOLUME);
  nrm = square_norm(work_fields[3], VOLUME, 1);
  if(g_cart_id == 0) {
    printf("# ||P2 psi - P2 P2 psi|| = %1.5e ", sqrt(nrm));
    if(sqrt(nrm) < eps) printf("#  -> passed\n\n");
    else printf("#  -> FAILED!\n\n");
    fflush(stdout);
  }

  project_left_D(work_fields[1], work_fields[0]);
  D_project_right(work_fields[2], work_fields[0]);
  diff(work_fields[3], work_fields[2], work_fields[1], VOLUME);
  nrm = square_norm(work_fields[3], VOLUME, 1);
  if(g_cart_id == 0) {
    printf("# ||P_L D psi - D P_R psi|| = %1.5e ", sqrt(nrm));
    if(sqrt(nrm) < eps) printf("#  -> passed\n\n");
    else printf("#  -> FAILED!\n\n");
    fflush(stdout);
    printf("\n######################\n");
    printf("# The following tests are only meaningful up to the precision little_D can be inverted for.\n");
    printf("# They might, therefore, be only useful in a small volume and/or a small condition number of little_D\n");
    printf("# The inversion precision (squared) is set to %e\n", little_solver_high_prec);
    printf("# So don't expect a precision much better than %1.2e in the following tests\n\n", eps*10);
  }

  dfl_sloppy_prec = 0;
  project_left(work_fields[1], work_fields[0]);
  project_left(work_fields[2], work_fields[1]);
  diff(work_fields[3], work_fields[2], work_fields[1], VOLUME);
  nrm = square_norm(work_fields[3], VOLUME, 1);
  if(g_cart_id == 0) {
    printf("# ||P_L^2 psi - P_L psi|| = %1.5e ", sqrt(nrm));
    if(sqrt(nrm) < eps*10) printf("#  -> passed\n\n");
    else printf("#  -> FAILED!\n\n");
    fflush(stdout);
  }

  project_right(work_fields[1], work_fields[0]);
  project_right(work_fields[2], work_fields[1]);
  diff(work_fields[3], work_fields[2], work_fields[1], VOLUME);
  nrm = square_norm(work_fields[3], VOLUME, 1);
  if(g_cart_id == 0) {
    printf("# ||P_R^2 psi - P_R psi|| = %1.5e ", sqrt(nrm));
    if(sqrt(nrm) < eps*10) printf("#  -> passed\n\n");
    else printf("#  -> FAILED!\n\n");
    fflush(stdout);
  }

  project_left(work_fields[1], work_fields[0]);
  project2(work_fields[2], work_fields[1]);
  nrm = square_norm(work_fields[2], VOLUME, 1);
  if(g_cart_id == 0) {
    printf("# ||P P_L psi|| = %1.5e ", sqrt(nrm));
    if(sqrt(nrm) < eps*10) printf("#  -> passed\n\n");
    else printf("#  -> FAILED!\n\n");
    fflush(stdout);
  }

  project2(work_fields[1], work_fields[0]);
  project_right(work_fields[2], work_fields[1]);
  nrm = square_norm(work_fields[2], VOLUME, 1);
  if(g_cart_id == 0) {
    printf("# ||P_R P psi|| = %1.5e ", sqrt(nrm));
    if(sqrt(nrm) < eps*10) printf("#  -> passed\n\n");
    else printf("#  -> FAILED!\n\n");
    fflush(stdout);
  }

  project2(work_fields[1], work_fields[0]);
  project(work_fields[2], work_fields[1]);
  D_psi(work_fields[3], work_fields[2]);
  project2(work_fields[2], work_fields[3]);
  diff(work_fields[3], work_fields[2], work_fields[1], VOLUME);
  nrm = square_norm(work_fields[3], VOLUME, 1);
  if(g_cart_id == 0) {
    printf("# ||P D A^-1 P psi - P psi|| = %1.5e ", sqrt(nrm));
    if(sqrt(nrm) < eps*10) printf("#  -> passed\n\n");
    else printf("#  -> FAILED!\n\n");
    fflush(stdout);
  }

  project2(work_fields[1], work_fields[0]);
  D_psi(work_fields[2], work_fields[1]);
  project(work_fields[3], work_fields[2]);
  project2(work_fields[2], work_fields[3]);
  diff(work_fields[3], work_fields[2], work_fields[1], VOLUME);
  nrm = square_norm(work_fields[3], VOLUME, 1);
  if(g_cart_id == 0) {
    printf("# ||P A^-1 D P psi - P psi|| = %1.5e ", sqrt(nrm));
    if(sqrt(nrm) < eps*10) printf("#  -> passed\n\n");
    else printf("#  -> FAILED!\n\n");
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
    printf("# ||P D P (P D P)^-1 psi - P psi|| = %1.5e ", sqrt(nrm));
    if(sqrt(nrm) < eps*10) printf("#  -> passed\n\n");
    else printf("#  -> FAILED!\n\n");
    fflush(stdout);
  }


  invert_little_D_spinor(work_fields[1], work_fields[0]);
  invert_little_D_eo_spinor(work_fields[2], work_fields[0]);
  diff(work_fields[3], work_fields[1], work_fields[2], VOLUME);
  nrm = square_norm(work_fields[3], VOLUME, 1);
  if(g_cart_id == 0) {
    printf("# ||A^-1 psi - A^-1_eo psi|| = %1.5e ", sqrt(nrm));
    if(sqrt(nrm) < eps*100) printf("#  -> passed\n\n");
    else printf("#  -> FAILED!\n\n");
    fflush(stdout);
  }


  invert_little_D_spinor(work_fields[1], work_fields[0]);
  apply_little_D_spinor(work_fields[2], work_fields[1]);
  project2(work_fields[3], work_fields[0]);
  diff(work_fields[1], work_fields[3], work_fields[2], VOLUME);
  nrm = square_norm(work_fields[1], VOLUME, 1);
  if(g_cart_id == 0) {
    printf("# ||A A^-1 psi - P psi|| = %1.5e ", sqrt(nrm));
    if(sqrt(nrm) < eps*10) printf("#  -> passed\n\n");
    else printf("#  -> FAILED!\n\n");
    fflush(stdout);
  }

  invert_little_D_spinor(work_fields[1], work_fields[0]);
  apply_little_D_spinor(work_fields[2], work_fields[1]);
  project2(work_fields[3], work_fields[0]);
  project2(work_fields[1], work_fields[2]);
  diff(work_fields[2], work_fields[3], work_fields[1], VOLUME);
  nrm = square_norm(work_fields[2], VOLUME, 1);
  if(g_cart_id == 0) {
    printf("# ||P A A^-1 psi - P psi|| = %1.5e", sqrt(nrm));
    if(sqrt(nrm) < eps*10) printf("#  -> passed\n\n");
    else printf("#  -> FAILED!\n\n");
    printf("\n######################\n");
    printf("# The following tests should be again fulfilled up to machine precision\n\n");
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

  project2(work_fields[1], work_fields[2]);
  diff(work_fields[2], work_fields[3], work_fields[1], VOLUME);
  nrm = square_norm(work_fields[2], VOLUME, 1);
  if(g_proc_id == 0) {
    printf("# ||(P D - A) phi_i || = %1.5e ", sqrt(nrm));
    if(sqrt(nrm) < eps*10) printf("#  -> passed\n\n");
    else printf("#  -> FAILED!\n\n");
    fflush(stdout);
  }

  reconstruct_global_field_GEN_ID(work_fields[1], block_list, 0, nb_blocks);
  apply_little_D_spinor(work_fields[3], work_fields[1]);
  D_psi(work_fields[2], work_fields[1]);

  project2(work_fields[1], work_fields[2]);
  diff(work_fields[2], work_fields[3], work_fields[1], VOLUME);
  nrm = square_norm(work_fields[2], VOLUME, 1);
  if(g_proc_id == 0) {
    printf("# ||(P D - A) phi || = %1.5e ", sqrt(nrm));
    if(sqrt(nrm) < eps*10) printf("#  -> passed\n\n");
    else printf("#  -> FAILED!\n\n");
    fflush(stdout);
  }

  apply_little_D_spinor(work_fields[3], work_fields[0]);
  project2(work_fields[1], work_fields[0]);
  D_psi(work_fields[2], work_fields[1]);

  project2(work_fields[1], work_fields[2]);
  diff(work_fields[2], work_fields[3], work_fields[1], VOLUME);
  nrm = square_norm(work_fields[2], VOLUME, 1);
  if(g_proc_id == 0) {
    printf("# ||P D P psi - A psi|| = %1.5e ", sqrt(nrm));
    if(sqrt(nrm) < eps*10) printf("#  -> passed\n\n");
    else printf("#  -> FAILED!\n\n");
    fflush(stdout);
  }

  /* check little projectors now */
  if(g_cart_id == 0) {
    printf("\n######################\n");
    printf("# Now we check the little little projection routines\n\n");
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
    printf("# ||lP2 v - lP2 lP2 v|| = %1.5e ", sqrt(nrm));
    if(sqrt(nrm) < eps*10) printf("#  -> passed\n\n");
    else printf("#  -> FAILED!\n\n");
    fflush(stdout);
  }

  little_P_L_D(work[11], work[10]);
  little_P_L_D(work[12], work[10]);
  ldiff(work[12], work[12], work[11], nb_blocks*g_N_s);
  nrm = lsquare_norm(work[12], nb_blocks*g_N_s, 1);
  if(g_cart_id == 0) {
    printf("# ||lP_L lD v - lP_L lD v|| = %1.5e ", sqrt(nrm));
    if(sqrt(nrm) < eps*10) printf("#  -> passed\n\n");
    else printf("#  -> FAILED!\n\n");
    fflush(stdout);
  }  

  little_P_L_D(work[11], work[10]);
  little_D_P_R(work[12], work[10]);
  ldiff(work[12], work[12], work[11], nb_blocks*g_N_s);
  nrm = lsquare_norm(work[12], nb_blocks*g_N_s, 1);
  if(g_cart_id == 0) {
    printf("# ||lP_L lD v - lD lP_R v|| = %1.5e ", sqrt(nrm));
    if(sqrt(nrm) < eps*10) printf("#  -> passed\n\n");
    else printf("#  -> FAILED!\n\n");
    fflush(stdout);
  }

  little_P_R(work[11], work[10]);
  little_P_R(work[12], work[11]);
  ldiff(work[12], work[12], work[11], nb_blocks*g_N_s);
  nrm = lsquare_norm(work[12], nb_blocks*g_N_s, 1);
  if(g_cart_id == 0) {
    printf("# ||lP_R^2 v - lP_R v|| = %1.5e ", sqrt(nrm));
    if(sqrt(nrm) < eps*10) printf("#  -> passed\n\n");
    else printf("#  -> FAILED!\n\n");
    fflush(stdout);
  }

  little_P_L(work[11], work[10]);
  little_P_L(work[12], work[11]);
  ldiff(work[12], work[12], work[11], nb_blocks*g_N_s);
  nrm = lsquare_norm(work[12], nb_blocks*g_N_s, 1);
  if(g_cart_id == 0) {
    printf("# ||lP_L^2 v - lP_L v|| = %1.5e ", sqrt(nrm));
    if(sqrt(nrm) < eps*10) printf("#  -> passed\n\n");
    else printf("#  -> FAILED!\n\n");
    fflush(stdout);
  }

  little_solver_high_prec = savelittle_solver_high_prec;
  free(phi[0]);
  free(phi);
  free(wphi);
  finalize_solver(work_fields, nr_wf);
  return(0);
}

void check_little_D_inversion(const int repro) {
  int i, j;
  int vol = block_list[0].volume;
  _Complex double *result;
  double dif;
  spinor ** work_fields = NULL;
  const int nr_wf = 1;

  init_solver_field(&work_fields, VOLUMEPLUSRAND, nr_wf);
  random_spinor_field_lexic(work_fields[0], repro, RN_GAUSS);
  if(init_dfl_projector == 0) {
    alloc_dfl_projector();
  }

  if(g_proc_id == 0) {
    printf("# Perform a test inversion of little_D using lGCR\n");
    printf("# This test might be only meaningful for not too large condition number of little_D\n\n");
  }

  result = calloc(nb_blocks * 9 * g_N_s, sizeof(_Complex double)); /*inner product of spinors with bases */

  /*initialize the local (block) parts of the spinor*/
  split_global_field_GEN(psi, work_fields[0], nb_blocks);

  for (i = 0; i < nb_blocks; ++i) {/* loop over blocks */
    /* compute inner product */
    for (j = 0; j < g_N_s; ++j) {/*loop over block.basis */
      /*       inprod[j + i*g_N_s] = block_scalar_prod(block_list[i].basis[j], psi[i], vol); */
      inprod[j + i*g_N_s] = scalar_prod(psi[i], block_list[i].basis[j], vol, 0);
      invvec[j + i*g_N_s] = 0.;
    }
  }


  gcr4complex(invvec, inprod, little_gmres_m_parameter, 1000, 1.e-24, 0, nb_blocks * g_N_s, 1, nb_blocks * 9 * g_N_s, 0, &little_D);

  little_D(result, invvec); /* This should be a proper inverse now */

#ifdef MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  ldiff(invvec, result, inprod, nb_blocks*g_N_s);
  dif = lsquare_norm(invvec, nb_blocks*g_N_s, 1);
  for (i = 0; i < nb_blocks; ++i) {/* loop over blocks */
    for (j = 0; j < 9*g_N_s; ++j) {/*loop over block.basis */
      invvec[j + i*g_N_s] = 0.;
      inprod[j + i*g_N_s] = 0.;
    }
  }

  if(g_proc_id == g_stdio_proc) {
    printf("# # check_little_D_inversion: squared residue found of size %1.5e! ", dif);
    if(dif < 1.e-24) printf("#  -> passed\n\n");
    else printf("#  -> FAILED\n\n");
  }

  finalize_solver(work_fields, nr_wf);
  free(result);
  return;
}

void check_local_D(const int repro)
{
  spinor * r[8];
  int j, vol = block_list[0].volume/2;
  double nrm;
  spinor ** work_fields = NULL;
  const int nr_wf = 7;

  init_solver_field(&work_fields, VOLUMEPLUSRAND, nr_wf);
  block_convert_lexic_to_eo(work_fields[0], work_fields[1], block_list[0].basis[0]);
  block_convert_eo_to_lexic(work_fields[2], work_fields[0], work_fields[1]);
  diff(work_fields[0], work_fields[2], block_list[0].basis[0], block_list[0].volume);
  nrm = square_norm(work_fields[0], block_list[0].volume, 0);
  if(g_proc_id == 0) {
    printf("# \nblock even/odd: ||psi - psi_recon|| = %1.5e\n", sqrt(nrm));
    printf("# next we compare local D against the Hopping matrix\n");
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
    if(g_c_sw > 0)
      assign_mul_one_sw_pm_imu_block(EE, work_fields[2], work_fields[0], g_mu, &block_list[j]);
    else 
      assign_mul_one_pm_imu(work_fields[2], work_fields[0], 1., vol);
    
    assign_add_mul_r(work_fields[2], g_spinor_field[DUM_DERI], 1., vol);

    /* Odd sites */
    Block_H_psi(&block_list[j], g_spinor_field[DUM_DERI], work_fields[0], OE);
    if(g_c_sw > 0)
      assign_mul_one_sw_pm_imu_block(OO,work_fields[3], work_fields[1], g_mu, &block_list[j]);
    else 
      assign_mul_one_pm_imu(work_fields[3], work_fields[1], 1., vol);
 
    assign_add_mul_r(work_fields[3], g_spinor_field[DUM_DERI], 1., vol);

    /* convert back to block spinor */
    block_convert_eo_to_lexic(work_fields[5], work_fields[2], work_fields[3]);

    diff(work_fields[4], work_fields[5], work_fields[6], block_list[0].volume);
    nrm = square_norm(work_fields[4], block_list[0].volume, 0);
    if(sqrt(nrm) > 1.e-12) {
      printf("# Check failed for local D against Hopping Matrix: ||delta|| = %1.5e block %d process %d\n", sqrt(nrm), j, g_proc_id);
    }
  }
  
  if(g_proc_id == 0) {
    printf("# ...done\n");
    printf("# Test Msap and Msap_eo to reduce residue\n");
    printf("# Expect something around 5.e-2 for the relative reduction with Ncycle=4, Niter=4\n\n");
  }
  /* check Msap and Msap_eo on a radom vector */
  random_spinor_field_lexic(work_fields[0], repro, RN_GAUSS);
  double nrm2 =  square_norm(work_fields[0], VOLUME, 1);
  if(g_proc_id == 0) {
    printf("# Initial residue ||r||^2 = %1.5e\n", nrm2);
  }
  zero_spinor_field(work_fields[1], VOLUME);
  Msap(work_fields[1], work_fields[0], 5,3);
  D_psi(work_fields[2], work_fields[1]);
  diff(work_fields[3], work_fields[2], work_fields[0], VOLUME);
  nrm = square_norm(work_fields[3], VOLUME, 1);
  if(g_proc_id == 0) {
    printf("# Msap relaxed the residue to ||r||^2 = %1.5e relative reduction %1.5e\n\n", nrm, nrm/nrm2);
  }

  zero_spinor_field(work_fields[1], VOLUME);
  Msap_eo(work_fields[1], work_fields[0], 5,3);
  D_psi(work_fields[2], work_fields[1]);
  diff(work_fields[3], work_fields[2], work_fields[0], VOLUME);
  nrm = square_norm(work_fields[3], VOLUME, 1);
  if(g_proc_id == 0) {
    printf("# Msap_eo relaxed the residue to ||r||^2 = %1.5e relative reduction %1.5e\n\n", nrm, nrm/nrm2);
    printf("# Now we test the block MR with even/odd\n");
    printf("# Expect 1.e-3 or so on each block\n\n");
  }

  for(j = 0; j < 6; j++) {
    r[j] = work_fields[j];
  }
  for(j = 0; j < nb_blocks; j++) {

    block_convert_lexic_to_eo(r[0], r[1], block_list[j].basis[0]);
    /* check even/odd inversion for Block_D_psi*/
    /* varphi_e in r[2] */
    if(g_c_sw > 0)
      assign_mul_one_sw_pm_imu_inv_block(EE,r[2], r[0], g_mu, &block_list[j]);
    else
      assign_mul_one_pm_imu_inv(r[2], r[0], +1., vol);

    Block_H_psi(&block_list[j], r[3], r[2], OE);
    /* a_odd = a_odd + b_odd */
    /* varphi_o in r[3] */
    assign_mul_add_r(r[3], -1., r[1], vol);
    /* psi_o in r[1] */
    if(g_c_sw > 0) {
      mrblk(r[1], r[3], r[4], 3, 1.e-31, 1, vol, &Msw_plus_block_psi, j);
    }
    else {
      mrblk(r[1], r[3], r[4], 3, 1.e-31, 1, vol, &Mtm_plus_block_psi, j);
    }

    Block_H_psi(&block_list[j], r[0], r[1], EO);
    assign(r[5],r[0],VOLUMEPLUSRAND);
    if(g_c_sw > 0)
      assign_mul_one_sw_pm_imu_inv_block(EE, r[0], r[5], g_mu, &block_list[j]);
    else
      mul_one_pm_imu_inv(r[0], +1., vol);
    /* a_even = a_even + b_even */
    /* check this sign +1 seems to be right in Msap_eo */
    assign_add_mul_r(r[2], r[0], -1., vol);

    block_convert_eo_to_lexic(r[4], r[2], r[1]);

    Block_D_psi(&block_list[j], r[5], r[4]);
    diff(r[0], block_list[j].basis[0], r[5], block_list[j].volume);
    nrm = square_norm(r[0], block_list[j].volume, 0);
    if(g_proc_id == 0) {
      printf("# mr_eo, block=%d: ||r||^2 = %1.5e\n", j, nrm);
    }
  }
  if( g_c_sw <= 0 ) {
    for(j = 0; j < nb_blocks; j++) {
      block_convert_lexic_to_eo(r[0], r[1], block_list[j].basis[0]);
      /* check even/odd inversion for Block_D_psi*/
      /* varphi_e in r[2] */
      if(g_c_sw > 0)
        assign_mul_one_sw_pm_imu_inv_block(EE,r[2],r[0], g_mu, &block_list[j]);
      else
        assign_mul_one_pm_imu_inv(r[2], r[0], +1., vol);
      
      Block_H_psi(&block_list[j], r[3], r[2], OE);
      /* a_odd = a_odd + b_odd */
      /* varphi_o in r[3] */
      assign_mul_add_r(r[3], -1., r[1], vol);
      /* psi_o in r[1] */
      if(g_c_sw > 0) {
        // FIXME: this cannot be correct!
        assign(r[5],r[3],VOLUMEPLUSRAND);
        assign_mul_one_sw_pm_imu_inv_block(OO, r[3], r[5], g_mu, &block_list[j]);
      }
      else {
        mul_one_pm_imu_inv(r[3], +1., vol);
      }
      
      if(g_c_sw > 0)
        mrblk(r[1], r[3], r[4], 3, 1.e-31, 1, vol, &Msw_plus_sym_block_psi, j);
      else
        mrblk(r[1], r[3], r[4], 3, 1.e-31, 1, vol, &Mtm_plus_sym_block_psi, j);
      
      Block_H_psi(&block_list[j], r[0], r[1], EO);
      
      if(g_c_sw > 0){
        assign(r[5],r[0],VOLUMEPLUSRAND);
        assign_mul_one_sw_pm_imu_inv_block(EE, r[0], r[5], g_mu, &block_list[j]);
      }
      else{
        mul_one_pm_imu_inv(r[0], +1., vol);}
      
      /* a_even = a_even + b_even */
      /* check this sign +1 seems to be right in Msap_eo */
      assign_add_mul_r(r[2], r[0], -1., vol);
      
      block_convert_eo_to_lexic(r[4], r[2], r[1]);
      
      Block_D_psi(&block_list[j], r[5], r[4]);
      diff(r[0], block_list[j].basis[0], r[5], block_list[j].volume);
      nrm = square_norm(r[0], block_list[j].volume, 0);
      if(g_proc_id == 0) {
        printf("# mr_eo (symmetric eo), block=%d: ||r||^2 = %1.5e\n", j, nrm);
      }
    }
  }
  finalize_solver(work_fields, nr_wf);
  return;
}




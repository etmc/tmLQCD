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
#include <complex.h>
#include "block.h"
#include "linalg/blas.h"
#include "solver/gcr4complex.h"
#include "solver/generate_dfl_subspace.h"
#include "xchange/little_field_gather.h"
#include "gamma.h"
#include "linalg_eo.h"
#include "little_D.h"


/* assume we have a little field w                       */
/* which has length 9*nb_blocks*N_s                      */
/* with usual order in space                             */
/*                                                       */
/* block[0], block[1], block[0], block[1], block[0]  ... */
/* local             , +t                , -t        ... */
/*                                                       */
/* block[0], block[1], block[0], block[1]                */
/* +z                , -z                                */
/* wasting some memory here...                           */

int dfl_subspace_updated = 1;

/* some lapack related stuff */
static int ONE = 1;
static _Complex double CONE, CZERO, CMONE;
static _Complex float CONE32, CZERO32, CMONE32;

enum{
  NONE = 0,
  T_UP = 1,
  T_DN = 2,
  X_UP = 3,
  X_DN = 4,
  Y_UP = 5,
  Y_DN = 6,
  Z_UP = 7,
  Z_DN = 8
} Direction;

void init_little_field_exchange(_Complex double * w);
void wait_little_field_exchange(const int mu);

void unit_little_D(_Complex double *v, _Complex double *w) {
  memcpy(v, w, nb_blocks*g_N_s*sizeof(_Complex double));

  return;
}

/** ANOTHER TESTING FUNCTION */
void invert_little_D_spinor(spinor *r, spinor *s){
  int i, j;
  spinor **psi;
  _Complex double *v, *w;
  psi = calloc(nb_blocks, sizeof(spinor));
  v = calloc(nb_blocks * 9 * g_N_s, sizeof(_Complex double));
  w = calloc(nb_blocks * 9 * g_N_s, sizeof(_Complex double));
  psi[0] = calloc(VOLUME+nb_blocks, sizeof(spinor));
  for(i = 1; i < nb_blocks; i++) {
    psi[i] = psi[i-1] + (VOLUME / nb_blocks) +1;
  }
  split_global_field_GEN(psi, s, nb_blocks); // ADAPT THIS
     
  for (j = 0; j < g_N_s; ++j) {/*loop over block.basis */
    for(i=0;i<nb_blocks;i++) {
      v[j + i*g_N_s] = scalar_prod(block_list[i].basis[j], psi[i], VOLUME/nb_blocks, 0);
    }
  }

  i = gcr4complex(w, v, 10, 100, little_solver_high_prec, 1, nb_blocks * g_N_s, 1, nb_blocks * 9 * g_N_s, 0, &little_D);
  if(g_proc_id == 0 && g_debug_level > 0) {
    printf("lgcr: %d iterations in invert_little_D_spinor\n", i);
  }
  
  for(i = 0; i < nb_blocks; i++) {
    mul(psi[i], w[i*g_N_s], block_list[i].basis[0], VOLUME/nb_blocks);
  }
  for(j = 1; j < g_N_s; ++j) {
    for(i = 0; i < nb_blocks; i++) {
      assign_add_mul(psi[i], block_list[i].basis[j], w[j+i*g_N_s], VOLUME/nb_blocks);
    }
  }
  reconstruct_global_field_GEN(r, psi, nb_blocks); // ADAPT THIS
      
  free(v);
  free(w);
  free(psi[0]);
  free(psi);
}


/** ANOTHER TESTING FUNCTION */
void invert_little_D_eo_spinor(spinor *r, spinor *s){
  int i, j, iter,i_o, i_e;
  spinor **psi;
  _Complex double *v, *w, *v_o, *v_e, * v_eo, * w_eo, * ctmp2;
  psi = calloc(nb_blocks, sizeof(spinor));
  v = calloc(nb_blocks * 9 * g_N_s, sizeof(_Complex double));
  w = calloc(nb_blocks * 9 * g_N_s, sizeof(_Complex double));
  v_e = calloc(nb_blocks * 9 * g_N_s, sizeof(_Complex double));
  v_o = calloc(nb_blocks * 9 * g_N_s, sizeof(_Complex double));
  v_eo = calloc(nb_blocks * 9 * g_N_s, sizeof(_Complex double));
  w_eo = calloc(nb_blocks * 9 * g_N_s, sizeof(_Complex double));
  ctmp2 = calloc(nb_blocks * 9 * g_N_s, sizeof(_Complex double));
  psi[0] = calloc(VOLUME+nb_blocks, sizeof(spinor));
  for(i = 1; i < nb_blocks; i++) {
    psi[i] = psi[i-1] + (VOLUME / nb_blocks) +1;
  }
  split_global_field_GEN(psi, s, nb_blocks); // ADAPT THIS 
  
  for (j = 0; j < g_N_s; ++j) {/*loop over block.basis */
    i_e=0;
    i_o=0;
    for(i=0;i<nb_blocks;i++) {
      v[j + i*g_N_s] = scalar_prod(block_list[i].basis[j], psi[i], VOLUME/nb_blocks, 0);
      if (block_list[i].evenodd==0) {
	v_eo[j+i_e*g_N_s] = v[j+i*g_N_s];
	i_e++;
      }
      if (block_list[i].evenodd==1) {
	v_eo[j+nb_blocks*g_N_s/2+i_o*g_N_s] = v[j+i*g_N_s];
	i_o++; 
      }
    }
  }
  
  little_D_ee_inv(v_e,v_eo);
  little_D_hop(1,v_o,v_e);
  little_Dhat_rhs(1,v_o,-1,v_eo);
  
  iter = gcr4complex(w_eo, v_o, 10, 100, 1e-31, 1, nb_blocks * g_N_s, 1, nb_blocks * 9 * g_N_s, 0, &little_D_sym);

  little_D_hop(0,ctmp2, w_eo);   
  little_D_ee_inv(w_eo,ctmp2);
  little_Dhat_rhs(0,w_eo, -1., v_e);
            
  for (j = 0; j < g_N_s; j++) {
    i_o=0;
    i_e=0;
    for(i = 0; i < nb_blocks; i++) {
      if (block_list[i].evenodd==0) {
	w[j + i*g_N_s] = w_eo[j + i_e*g_N_s];
	i_e++;
      }
      if (block_list[i].evenodd==1) {
	w[j + i*g_N_s] = w_eo[j + nb_blocks*g_N_s/2+i_o*g_N_s];
	i_o++;
      }
    }
  }

  if(g_proc_id == 0 && g_debug_level > 0) {
    printf("lgcr: %d iterations in invert_little_D_eo_spinor\n", iter);
  }
  
  for(i = 0; i < nb_blocks; i++) {
    mul(psi[i], w[i*g_N_s], block_list[i].basis[0], VOLUME/nb_blocks);
  }
  for(j = 1; j < g_N_s; ++j) {
    for(i = 0; i < nb_blocks; i++) {
      assign_add_mul(psi[i], block_list[i].basis[j], w[j+i*g_N_s], VOLUME/nb_blocks);
    }
  }
  reconstruct_global_field_GEN(r, psi, nb_blocks); // ADAPT THIS

  free(v);
  free(w);
  free(w_eo);
  free(v_eo);
  free(v_o);
  free(v_e);
  free(ctmp2);
  free(psi[0]);
  free(psi);
}


void project2(spinor * const out, spinor * const in);

/** ANOTHER TESTING FUNCTION */
void apply_little_D_spinor(spinor *r, spinor *s){
  int i,j, k;
  spinor **psi;
  _Complex double *v, *w;

  psi = (spinor **)calloc(nb_blocks, sizeof(spinor *));
  v = calloc(nb_blocks * 9 * g_N_s, sizeof(_Complex double));
  w = calloc(nb_blocks * 9 * g_N_s, sizeof(_Complex double));
  psi[0] = calloc(VOLUME + nb_blocks, sizeof(spinor));
  for(i = 1; i < nb_blocks; i++) {
    psi[i] = psi[i-1] + (VOLUME / nb_blocks) + 1;
  }
  split_global_field_GEN(psi, s, nb_blocks);  

  for (j = 0; j < g_N_s; ++j) {
    for(i = 0; i < nb_blocks; i++) v[j + i*g_N_s] = scalar_prod(block_list[i].basis[j], psi[i], VOLUME/nb_blocks, 0);
  }

  little_D(w, v);

  for(i = 0; i < nb_blocks; i++) {
    mul(psi[i], w[i*g_N_s], block_list[i].basis[0], VOLUME/nb_blocks);
  }
  for(j = 1; j < g_N_s; ++j) {
    for(i = 0; i < nb_blocks; i++){
      assign_add_mul(psi[i], block_list[i].basis[j], w[i*g_N_s + j], VOLUME/nb_blocks);
    }
  }
  reconstruct_global_field_GEN(r, psi, nb_blocks);
  
  free(v);
  free(w);
  free(psi[0]);
  free(psi);
}


#define _PSWITCH(s) s
#if (defined NOF77UNDERSCORE || defined NOF77_)
#define _MV(x) zgemv
#else
#define _MV(x) zgemv_
#endif
#define _C_TYPE _Complex double

#include"little_D_body.c"

#undef _C_TYPE
#undef _PSWITCH
#undef _MV

#define _PSWITCH(s) s ## 32
#if (defined NOF77UNDERSCORE || defined NOF77_)
#define _MV(x) cgemv
#else
#define _MV(x) cgemv_
#endif
#define _C_TYPE _Complex float

#include"little_D_body.c"

#undef _C_TYPE
#undef _PSWITCH
#undef _MV



void init_little_field_exchange(_Complex double * w) {
#ifdef MPI
  int i = 0;
#  if (defined PARALLELT || defined PARALLELX)
  int no_dirs = 2;
#  elif (defined PARALLELXT || defined PARALLELXY || defined PARALLELXYZ)
  int no_dirs = 4;
#  elif (defined PARALLELXYT || defined PARALLELXYZT)
  int no_dirs = 6;
#  endif
  if(waitcount != 0) {
    if(g_proc_id == 0) {
      fprintf(stderr, "last little_field_exchange not finished! Aborting...\n");
    }
    exit(-1);
  }
  /* z-dir requires special treatment! */
  for(i = 0; i < no_dirs; i+=nb_blocks) {
    /* send to the right, receive from the left */
    MPI_Isend((void*)w, nb_blocks*g_N_s, MPI_DOUBLE_COMPLEX, g_nb_list[i], 
              i, g_cart_grid, &lrequests[2*i]);
    MPI_Irecv((void*)(w + nb_blocks*(i+2)*g_N_s), nb_blocks*g_N_s, MPI_DOUBLE_COMPLEX, g_nb_list[i+1], 
              i, g_cart_grid, &lrequests[2*i+1]);
    
    /* send to the left, receive from the right */
    MPI_Isend((void*)w, nb_blocks*g_N_s, MPI_DOUBLE_COMPLEX, g_nb_list[i+1], 
              i+1, g_cart_grid, &lrequests[2*i+2]);
    MPI_Irecv((void*)(w + nb_blocks*(i+1)*g_N_s), nb_blocks*g_N_s, MPI_DOUBLE_COMPLEX, g_nb_list[i], 
              i+1, g_cart_grid, &lrequests[2*i+3]);
    waitcount += 4;
  }
#  ifdef PARALLELXYZT
  /* send to the right, receive from the left */
  i = 6;
  MPI_Isend((void*)(w + g_N_s), g_N_s, MPI_DOUBLE_COMPLEX, g_nb_list[i], 
            i, g_cart_grid, &lrequests[2*i]);
  MPI_Irecv((void*)(w + (nb_blocks*(i+1)+1)*g_N_s), g_N_s, MPI_DOUBLE_COMPLEX, g_nb_list[i+1], 
            i, g_cart_grid, &lrequests[2*i+1]);
  
  /* send to the left, receive from the right */
  MPI_Isend((void*)w, g_N_s, MPI_DOUBLE_COMPLEX, g_nb_list[i+1], 
            i+1, g_cart_grid, &lrequests[2*i+2]);
  MPI_Irecv((void*)(w + nb_blocks*(i+1)*g_N_s), g_N_s, MPI_DOUBLE_COMPLEX, g_nb_list[i], 
            i+1, g_cart_grid, &lrequests[2*i+3]);
  waitcount += 4;
#  endif
#endif
  return;
}

void wait_little_field_exchange(const int mu) {
#ifdef MPI
  int err;
  err = MPI_Waitall(2, &lrequests[2*mu], &lstatus[2*mu]);
  waitcount -= 2;
#endif
  return;
}



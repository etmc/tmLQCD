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
#include "block.h"
#include "linalg_eo.h"
#include "little_D.h"


/* assume we have a little field w                       */
/* which has length 9*nb_blocks*N_s                      */
/* with usual order in space                             */
/* nb_blocks = 2 currently fixed                         */
/* and blocks devide z-direction by 2                    */
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

  i = gcr4complex(w, v, 10, 100, 1e-31, 1, nb_blocks * g_N_s, 1, nb_blocks * 9 * g_N_s, &little_D);
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
  
  iter = gcr4complex(w_eo, v_o, 10, 100, 1e-31, 1, nb_blocks * g_N_s, 1, nb_blocks * 9 * g_N_s, &little_D_sym);

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

  if (g_debug_level > 2){
    if (!g_cart_id) {
      for (j = 0; j < nb_blocks* g_N_s; ++j) {
        printf("LITTLE_D for 0: v[%u] = %1.5e + %1.5e i\n", j, creal(v[j]), cimag(v[j]));
      }
    }
#ifdef MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
  }

  if (g_debug_level > 4) {
    for (k = 1; k < 16; ++k) {
      if (g_cart_id == k) {
        for (j = 0; j < nb_blocks* g_N_s; ++j) {
          printf("LITTLE_D for %u: v[%u] = %1.5e + %1.5e i\n", k, j, creal(v[j]), cimag(v[j]));
        }
      }
#ifdef MPI
      MPI_Barrier(MPI_COMM_WORLD);
#endif
    }
  }

  little_D(w, v);

  if (g_debug_level > 2){
    if (!g_cart_id){
      for (j = 0; j < nb_blocks * g_N_s; ++j) {
        printf("LITTLE_D for 0: w[%u] = %1.5e + %1.5e i\n", j, creal(w[j]), cimag(w[j]));
      }
    }
#ifdef MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
  }

  if (g_debug_level > 4) {
    for (k = 1; k < 16; ++k) {
      if (g_cart_id == k) {
        for (j = 0; j < nb_blocks* g_N_s; ++j) {
          printf("LITTLE_D for %u: w[%u] = %1.5e + %1.5e i\n", k, j, creal(w[j]), cimag(w[j]));
        }
      }
#ifdef MPI
      MPI_Barrier(MPI_COMM_WORLD);
#endif
    }
  }
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


void alt_little_field_gather(_Complex double * w) {
#ifdef MPI
  MPI_Status status;
  int size = 25 * g_N_s * sizeof(_Complex double);
  _Complex double *buf = malloc(size);
  MPI_Buffer_attach((void*)buf, size);

  /* LOWER BLOCK */

  /* Send t up */
  MPI_Bsend(w, g_N_s, MPI_DOUBLE_COMPLEX, g_nb_t_up, T_UP, g_cart_grid);
  MPI_Recv(w + 4 * g_N_s, g_N_s, MPI_DOUBLE_COMPLEX, g_nb_t_dn, T_UP, g_cart_grid, &status);

  /* Send t down */
  MPI_Bsend(w, g_N_s, MPI_DOUBLE_COMPLEX, g_nb_t_dn, T_DN, g_cart_grid);
  MPI_Recv(w + 2 * g_N_s, g_N_s, MPI_DOUBLE_COMPLEX, g_nb_t_up, T_DN, g_cart_grid, &status);

  /* Send x up */
  MPI_Bsend(w, g_N_s, MPI_DOUBLE_COMPLEX, g_nb_x_up, X_UP, g_cart_grid);
  MPI_Recv(w + 8 * g_N_s, g_N_s, MPI_DOUBLE_COMPLEX, g_nb_x_dn, X_UP, g_cart_grid, &status);

  /* Send x down */
  MPI_Bsend(w, g_N_s, MPI_DOUBLE_COMPLEX, g_nb_x_dn, X_DN, g_cart_grid);
  MPI_Recv(w + 6 * g_N_s, g_N_s, MPI_DOUBLE_COMPLEX, g_nb_x_up, X_DN, g_cart_grid, &status);

  /* Send y up */
  MPI_Bsend(w, g_N_s, MPI_DOUBLE_COMPLEX, g_nb_y_up, Y_UP, g_cart_grid);
  MPI_Recv(w + 12 * g_N_s, g_N_s, MPI_DOUBLE_COMPLEX, g_nb_y_dn, Y_UP, g_cart_grid, &status);

  /* Send y down */
  MPI_Bsend(w, g_N_s, MPI_DOUBLE_COMPLEX, g_nb_y_dn, Y_DN, g_cart_grid);
  MPI_Recv(w + 10 * g_N_s, g_N_s, MPI_DOUBLE_COMPLEX, g_nb_y_up, Y_DN, g_cart_grid, &status);

  /* Send z up */
  memcpy(w + 17 * g_N_s, w, g_N_s * sizeof(_Complex double));

  /* Send z down */
  MPI_Bsend(w, g_N_s, MPI_DOUBLE_COMPLEX, g_nb_z_dn, Z_DN, g_cart_grid);
  MPI_Recv(w + 15 * g_N_s, g_N_s, MPI_DOUBLE_COMPLEX, g_nb_z_up, Z_DN, g_cart_grid, &status);

  /* END LOWER BLOCK */

  MPI_Barrier(MPI_COMM_WORLD);

  /* UPPER BLOCK */

  /* Send t up */
  MPI_Bsend(w + g_N_s, g_N_s, MPI_DOUBLE_COMPLEX, g_nb_t_up, T_UP, g_cart_grid);
  MPI_Recv(w + 5 * g_N_s, g_N_s, MPI_DOUBLE_COMPLEX, g_nb_t_dn, T_UP, g_cart_grid, &status);

  /* Send t down */
  MPI_Bsend(w + g_N_s, g_N_s, MPI_DOUBLE_COMPLEX, g_nb_t_dn, T_DN, g_cart_grid);
  MPI_Recv(w + 3 * g_N_s, g_N_s, MPI_DOUBLE_COMPLEX, g_nb_t_up, T_DN, g_cart_grid, &status);

  /* Send x up */
  MPI_Bsend(w + g_N_s, g_N_s, MPI_DOUBLE_COMPLEX, g_nb_x_up, X_UP, g_cart_grid);
  MPI_Recv(w + 9 * g_N_s, g_N_s, MPI_DOUBLE_COMPLEX, g_nb_x_dn, X_UP, g_cart_grid, &status);

  /* Send x down */
  MPI_Bsend(w + g_N_s, g_N_s, MPI_DOUBLE_COMPLEX, g_nb_x_dn, X_DN, g_cart_grid);
  MPI_Recv(w + 7 * g_N_s, g_N_s, MPI_DOUBLE_COMPLEX, g_nb_x_up, X_DN, g_cart_grid, &status);

  /* Send y up */
  MPI_Bsend(w + g_N_s, g_N_s, MPI_DOUBLE_COMPLEX, g_nb_y_up, Y_UP, g_cart_grid);
  MPI_Recv(w + 13 * g_N_s, g_N_s, MPI_DOUBLE_COMPLEX, g_nb_y_dn, Y_UP, g_cart_grid, &status);

  /* Send y down */
  MPI_Bsend(w + g_N_s, g_N_s, MPI_DOUBLE_COMPLEX, g_nb_y_dn, Y_DN, g_cart_grid);
  MPI_Recv(w + 11 * g_N_s, g_N_s, MPI_DOUBLE_COMPLEX, g_nb_y_up, Y_DN, g_cart_grid, &status);

  /* Send z up */
  MPI_Bsend(w + g_N_s, g_N_s, MPI_DOUBLE_COMPLEX, g_nb_z_up, Z_UP, g_cart_grid);
  MPI_Recv(w + 16 * g_N_s, g_N_s, MPI_DOUBLE_COMPLEX, g_nb_z_dn, Z_UP, g_cart_grid, &status);

  /* Send z down */
  memcpy(w + 14 * g_N_s, w + g_N_s, g_N_s * sizeof(_Complex double));

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Buffer_detach((void*)buf, &size);

  free(buf);
#endif
  return;
}

#ifdef MPI
MPI_Request lrequests[16];
MPI_Status lstatus[16];
int waitcount = 0;
#endif


void little_field_gather(_Complex double * w) {
#ifdef MPI
  int err, bt, bx, by, bz, pm, ib;
  _Complex double *wt, *wx, *wy, *wz;
  _Complex double *wt_buf, *wx_buf, *wy_buf, *wz_buf, *w_buf, *w_source, *w_dest;
  /************************************************************************/
  /* This routine has been extended for multi_dimensional blocking        */
  /* by Claude Tadonki (claude.tadonki@u-psud.fr) from PetaQCD project    */
  /* June 2010                                                            */
  /************************************************************************/

  w_buf = calloc(8 * nb_blocks * g_N_s, sizeof(_Complex double)); // +-t +-x +-y +-z

  wt = w + ( 0*(2*nb_blocks) + nb_blocks ) * g_N_s; // Were data in the direction t starts
  wx = w + ( 1*(2*nb_blocks) + nb_blocks ) * g_N_s; // Were data in the direction x starts
  wy = w + ( 2*(2*nb_blocks) + nb_blocks ) * g_N_s; // Were data in the direction y starts
  wz = w + ( 3*(2*nb_blocks) + nb_blocks ) * g_N_s; // Were data in the direction z starts

  wt_buf = w_buf + ( 0*(2*nb_blocks)) * g_N_s; // Were data in the direction t starts
  wx_buf = w_buf + ( 1*(2*nb_blocks)) * g_N_s; // Were data in the direction x starts
  wy_buf = w_buf + ( 2*(2*nb_blocks)) * g_N_s; // Were data in the direction y starts
  wz_buf = w_buf + ( 3*(2*nb_blocks)) * g_N_s; // Were data in the direction z starts

  /* We first exchange the fields regardless of block considerations                   */
  /* The data need to be received in an intermediate buffer because of later shuffling */

  /* Send t up */
  MPI_Isend(w, nb_blocks * g_N_s, MPI_DOUBLE_COMPLEX, g_nb_t_up, T_UP, g_cart_grid, &lrequests[0]);
  MPI_Irecv(wt_buf + nb_blocks * g_N_s, nb_blocks * g_N_s, MPI_DOUBLE_COMPLEX, g_nb_t_dn, T_UP, g_cart_grid, &lrequests[1]);

  /* Send t down */
  MPI_Isend(w, nb_blocks * g_N_s, MPI_DOUBLE_COMPLEX, g_nb_t_dn, T_DN, g_cart_grid, &lrequests[2]);
  MPI_Irecv(wt_buf, nb_blocks * g_N_s, MPI_DOUBLE_COMPLEX, g_nb_t_up, T_DN, g_cart_grid, &lrequests[3]);

  /* Send x up */
  MPI_Isend(w, nb_blocks * g_N_s, MPI_DOUBLE_COMPLEX, g_nb_x_up, X_UP, g_cart_grid, &lrequests[4]);
  MPI_Irecv(wx_buf + nb_blocks * g_N_s, nb_blocks * g_N_s, MPI_DOUBLE_COMPLEX, g_nb_x_dn, X_UP, g_cart_grid, &lrequests[5]);

  /* Send x down */
  MPI_Isend(w, nb_blocks * g_N_s, MPI_DOUBLE_COMPLEX, g_nb_x_dn, X_DN, g_cart_grid, &lrequests[6]);
  MPI_Irecv(wx_buf, nb_blocks * g_N_s, MPI_DOUBLE_COMPLEX, g_nb_x_up, X_DN, g_cart_grid, &lrequests[7]);

  /* Send y up */
  MPI_Isend(w, nb_blocks * g_N_s, MPI_DOUBLE_COMPLEX, g_nb_y_up, Y_UP, g_cart_grid, &lrequests[8]);
  MPI_Irecv(wy_buf + nb_blocks * g_N_s, nb_blocks * g_N_s, MPI_DOUBLE_COMPLEX, g_nb_y_dn, Y_UP, g_cart_grid, &lrequests[9]);

  /* Send y down */
  MPI_Isend(w, nb_blocks * g_N_s, MPI_DOUBLE_COMPLEX, g_nb_y_dn, Y_DN, g_cart_grid, &lrequests[10]);
  MPI_Irecv(wy_buf, nb_blocks * g_N_s, MPI_DOUBLE_COMPLEX, g_nb_y_up, Y_DN, g_cart_grid, &lrequests[11]);

  /* Send z up */
  MPI_Isend(w, nb_blocks * g_N_s, MPI_DOUBLE_COMPLEX, g_nb_z_up, Z_UP, g_cart_grid, &lrequests[12]);
  MPI_Irecv(wz_buf + nb_blocks * g_N_s, nb_blocks * g_N_s, MPI_DOUBLE_COMPLEX, g_nb_z_dn, Z_UP, g_cart_grid, &lrequests[13]);

  /* Send z down */
  MPI_Isend(w, nb_blocks * g_N_s, MPI_DOUBLE_COMPLEX, g_nb_z_dn, Z_DN, g_cart_grid, &lrequests[14]);
  MPI_Irecv(wz_buf, nb_blocks * g_N_s, MPI_DOUBLE_COMPLEX, g_nb_z_up, Z_DN, g_cart_grid, &lrequests[15]);
  
  err = MPI_Waitall(16, lrequests, lstatus);
  
  /* We now correct the field according to block partitionning               */
  /* We could have avoid the previous corresponding MPI communication        */
  /* We proceed like this for code simplicity, maybe will be optimized later */
  
  for(pm = 0; pm < 8; pm++) {
    for(bt = 0; bt < nblks_t; bt++) {
      for(bx = 0; bx < nblks_x; bx++) {
	for(by = 0; by < nblks_y; by++) {
	  for(bz = 0; bz < nblks_z; bz++) {
	    ib = block_index(bt, bx, by, bz) * g_N_s;
	    switch(pm){ 
	    case 0: /* Direction +t */
	      w_dest = wt + ib;
	      if( bt == nblks_t - 1 ) {ib = block_index(0, bx, by, bz) * g_N_s; w_source = wt_buf + ib;}					 // got it from the MPI exchange
	      else  {ib = block_index(bt + 1, bx, by, bz) * g_N_s; w_source = w + ib;}										 // got it from the diagonal block
	      break; 
	    case 1: /* Direction -t */
	      w_dest = wt + ib + nb_blocks * g_N_s;
	      if( bt == 0 ) {ib = block_index(nblks_t - 1, bx, by, bz) * g_N_s; w_source = wt_buf + ib + nb_blocks * g_N_s;} // got it from the MPI exchange
	      else  {ib = block_index(bt - 1, bx, by, bz) * g_N_s;w_source = w + ib;}										 // got it from the diagonal block
	      break; 
	    case 2: /* Direction +x */
	      w_dest = wx + ib;
	      if( bx == nblks_x - 1 ) {ib = block_index(bt, 0, by, bz) * g_N_s; w_source = wx_buf + ib;}					 // got it from the MPI exchange
	      else  {ib = block_index(bt, bx + 1, by, bz) * g_N_s; w_source = w + ib;}									     // got it from the diagonal block
	      break; 
	    case 3: /* Direction -x */
	      w_dest = wx + ib + nb_blocks * g_N_s;
	      if( bx == 0 ) {ib = block_index(bt, nblks_x - 1, by, bz) * g_N_s; w_source = wx_buf + ib + nb_blocks * g_N_s;} // got it from the MPI exchange
	      else  {ib = block_index(bt, bx - 1, by, bz) * g_N_s;w_source = w + ib;}									     // got it from the diagonal block
	      break; 
	    case 4: /* Direction +y */
	      w_dest = wy + ib;
	      if( by == nblks_y - 1 ) {ib = block_index(bt, bx, 0, bz) * g_N_s; w_source = wy_buf + ib;}			         // got it from the MPI exchange
	      else  {ib = block_index(bt, bx, by + 1, bz) * g_N_s; w_source = w + ib;}									     // got it from the diagonal block
	      break; 
	    case 5: /* Direction -y */
	      w_dest = wy + ib + nb_blocks * g_N_s;
	      if( by == 0 ) {ib = block_index(bt, bx, nblks_y - 1, bz) * g_N_s; w_source = wy_buf + ib + nb_blocks * g_N_s;} // got it from the MPI exchange
	      else  {ib = block_index(bt, bx, by - 1, bz) * g_N_s;w_source = w + ib;}									     // got it from the diagonal block
	      break; 
	    case 6: /* Direction +z */
	      w_dest = wz + ib;
	      if( bz == nblks_z - 1 ) {ib = block_index(bt, bx, by, 0) * g_N_s; w_source = wz_buf + ib;	}		             // got it from the MPI exchange
	      else  {ib = block_index(bt, bx, by, bz + 1) * g_N_s; w_source = w + ib;	}						             // got it from the diagonal block
	      break; 
	    case 7: /* Direction -z */
	      w_dest = wz + ib + nb_blocks * g_N_s;
	      if( bz == 0 ) {ib = block_index(bt, bx, by, nblks_z - 1) * g_N_s; w_source = wz_buf + ib + nb_blocks * g_N_s;} // got it from the MPI exchange
	      else  {ib = block_index(bt, bx, by, bz - 1) * g_N_s; w_source = w + ib; }                                      // got it from the diagonal block
	      break; 
	      
	    default: 
	      w_dest = NULL;
	      w_source = NULL;
	    }
	    memcpy(w_dest, w_source, g_N_s * sizeof(_Complex double));
	  }
	}
      }
    }
  }
  
  free(w_buf);
  
#endif
  return;
}

void little_field_gather_eo(int eo, _Complex double * w) {
#ifdef MPI
  int err, bt, bx, by, bz, pm, ib,ib2;
  _Complex double *wt, *wx, *wy, *wz;
  _Complex double *wt_buf, *wx_buf, *wy_buf, *wz_buf, *w_buf, *w_source, *w_dest;

  w_buf = calloc(8 * nb_blocks * g_N_s, sizeof(_Complex double)); // +-t +-x +-y +-z
  
  wt = w + ( 0*(2*nb_blocks) + nb_blocks ) * g_N_s; // Were data in the direction t starts
  wx = w + ( 1*(2*nb_blocks) + nb_blocks ) * g_N_s; // Were data in the direction x starts
  wy = w + ( 2*(2*nb_blocks) + nb_blocks ) * g_N_s; // Were data in the direction y starts
  wz = w + ( 3*(2*nb_blocks) + nb_blocks ) * g_N_s; // Were data in the direction z starts

  wt_buf = w_buf + ( 0*(2*nb_blocks)) * g_N_s; // Were data in the direction t starts
  wx_buf = w_buf + ( 1*(2*nb_blocks)) * g_N_s; // Were data in the direction x starts
  wy_buf = w_buf + ( 2*(2*nb_blocks)) * g_N_s; // Were data in the direction y starts
  wz_buf = w_buf + ( 3*(2*nb_blocks)) * g_N_s; // Were data in the direction z starts

  /* We first exchange the fields regardless of block considerations                   */
  /* The data need to be received in an intermediate buffer because of later shuffling */

  /* Send t up */
  MPI_Isend(w, nb_blocks * g_N_s, MPI_DOUBLE_COMPLEX, g_nb_t_up, T_UP, g_cart_grid, &lrequests[0]);
  MPI_Irecv(wt_buf + nb_blocks * g_N_s, nb_blocks * g_N_s, MPI_DOUBLE_COMPLEX, g_nb_t_dn, T_UP, g_cart_grid, &lrequests[1]);

  /* Send t down */
  MPI_Isend(w, nb_blocks * g_N_s, MPI_DOUBLE_COMPLEX, g_nb_t_dn, T_DN, g_cart_grid, &lrequests[2]);
  MPI_Irecv(wt_buf, nb_blocks * g_N_s, MPI_DOUBLE_COMPLEX, g_nb_t_up, T_DN, g_cart_grid, &lrequests[3]);

  /* Send x up */
  MPI_Isend(w, nb_blocks * g_N_s, MPI_DOUBLE_COMPLEX, g_nb_x_up, X_UP, g_cart_grid, &lrequests[4]);
  MPI_Irecv(wx_buf + nb_blocks * g_N_s, nb_blocks * g_N_s, MPI_DOUBLE_COMPLEX, g_nb_x_dn, X_UP, g_cart_grid, &lrequests[5]);

  /* Send x down */
  MPI_Isend(w, nb_blocks * g_N_s, MPI_DOUBLE_COMPLEX, g_nb_x_dn, X_DN, g_cart_grid, &lrequests[6]);
  MPI_Irecv(wx_buf, nb_blocks * g_N_s, MPI_DOUBLE_COMPLEX, g_nb_x_up, X_DN, g_cart_grid, &lrequests[7]);

  /* Send y up */
  MPI_Isend(w, nb_blocks * g_N_s, MPI_DOUBLE_COMPLEX, g_nb_y_up, Y_UP, g_cart_grid, &lrequests[8]);
  MPI_Irecv(wy_buf + nb_blocks * g_N_s, nb_blocks * g_N_s, MPI_DOUBLE_COMPLEX, g_nb_y_dn, Y_UP, g_cart_grid, &lrequests[9]);

  /* Send y down */
  MPI_Isend(w, nb_blocks * g_N_s, MPI_DOUBLE_COMPLEX, g_nb_y_dn, Y_DN, g_cart_grid, &lrequests[10]);
  MPI_Irecv(wy_buf, nb_blocks * g_N_s, MPI_DOUBLE_COMPLEX, g_nb_y_up, Y_DN, g_cart_grid, &lrequests[11]);

  /* Send z up */
  MPI_Isend(w, nb_blocks * g_N_s, MPI_DOUBLE_COMPLEX, g_nb_z_up, Z_UP, g_cart_grid, &lrequests[12]);
  MPI_Irecv(wz_buf + nb_blocks * g_N_s, nb_blocks * g_N_s, MPI_DOUBLE_COMPLEX, g_nb_z_dn, Z_UP, g_cart_grid, &lrequests[13]);

  /* Send z down */
  MPI_Isend(w, nb_blocks * g_N_s, MPI_DOUBLE_COMPLEX, g_nb_z_dn, Z_DN, g_cart_grid, &lrequests[14]);
  MPI_Irecv(wz_buf, nb_blocks * g_N_s, MPI_DOUBLE_COMPLEX, g_nb_z_up, Z_DN, g_cart_grid, &lrequests[15]);
  
  err = MPI_Waitall(16, lrequests, lstatus);
  
  /* We now correct the field according to block partitionning               */
  /* We could have avoid the previous corresponding MPI communication        */
  /* We proceed like this for code simplicity, maybe will be optimized later */

  for(pm = 0; pm < 8; pm++) {
    ib2=0;
    for(bt = 0; bt < nblks_t; bt++) {
      for(bx = 0; bx < nblks_x; bx++) {
	for(by = 0; by < nblks_y; by++) {
	  for(bz = 0; bz < nblks_z; bz++) {
	    if ((bt+bx+by+bz)%2==eo) {
	      ib2 = index_block_eo[block_index(bt, bx, by, bz)] * g_N_s;
	      
	      switch(pm){ 
	      case 0: /* Direction +t */
		w_dest = wt + ib2;
		if( bt == nblks_t - 1 ) {ib = index_block_eo[block_index(0,bx, by,bz)] * g_N_s; w_source = wt_buf + ib;	}		             // got it from the MPI exchange
		else  {ib = index_block_eo[block_index(bt+1, bx, by, bz)] * g_N_s; w_source = w + ib;	}						             // got it from the diagonal block
		break; 
	      case 1: /* Direction -t */
		w_dest = wt + ib2 + nb_blocks * g_N_s;
		if( bt == 0) {ib = index_block_eo[block_index(nblks_t-1, bx,by,bz)] * g_N_s; w_source = wt_buf + ib + nb_blocks * g_N_s;} // got it from the MPI exchange
		else  {ib = index_block_eo[block_index(bt-1,bx, by, bz)] * g_N_s; w_source = w + ib; }                                      // got it from the diagonal block
		break; 
	      case 2: /* Direction +x */
		w_dest = wx + ib2;
		if( bx == nblks_x - 1 ) {ib = index_block_eo[block_index(bt, 0, by,bz)] * g_N_s; w_source = wx_buf + ib;	}		             // got it from the MPI exchange
		else  {ib = index_block_eo[block_index(bt, bx+1, by, bz)] * g_N_s; w_source = w + ib;	}						             // got it from the diagonal block
		break; 
	      case 3: /* Direction -x */
		w_dest = wx + ib2 + nb_blocks * g_N_s;
		if( bx == 0) {ib = index_block_eo[block_index(bt, nblks_x-1, by,bz)] * g_N_s; w_source = wx_buf + ib + nb_blocks * g_N_s;} // got it from the MPI exchange
		else  {ib = index_block_eo[block_index(bt, bx-1, by, bz)] * g_N_s; w_source = w + ib; }                                      // got it from the diagonal block
		break; 
	      case 4: /* Direction +y */
		w_dest = wy + ib2;
		if( by == nblks_y - 1 ) {ib = index_block_eo[block_index(bt, bx, 0,bz)] * g_N_s; w_source = wy_buf + ib;	}		             // got it from the MPI exchange
		else  {ib = index_block_eo[block_index(bt, bx, by+1, bz)] * g_N_s; w_source = w + ib;	}						             // got it from the diagonal block
		break; 
	      case 5: /* Direction -y */
		w_dest = wy + ib2 + nb_blocks * g_N_s;
		if( by == 0) {ib = index_block_eo[block_index(bt, bx, nblks_y-1, bz)] * g_N_s; w_source = wy_buf + ib + nb_blocks * g_N_s;} // got it from the MPI exchange
		else  {ib = index_block_eo[block_index(bt, bx, by-1, bz)] * g_N_s; w_source = w + ib; }                                      // got it from the diagonal block
		break; 
	      case 6: /* Direction +z */
		w_dest = wz + ib2;
		if( bz == nblks_z - 1 ) {ib = index_block_eo[block_index(bt, bx, by, 0)] * g_N_s; w_source = wz_buf + ib;	}		             // got it from the MPI exchange
		else  {ib = index_block_eo[block_index(bt, bx, by, bz + 1)] * g_N_s; w_source = w + ib;	}						             // got it from the diagonal block
		break; 
	      case 7: /* Direction -z */
		w_dest = wz + ib2 + nb_blocks * g_N_s;
		if( bz == 0) {ib = index_block_eo[block_index(bt, bx, by, nblks_z - 1)] * g_N_s; w_source = wz_buf + ib + nb_blocks * g_N_s;} // got it from the MPI exchange
		else  {ib = index_block_eo[block_index(bt, bx, by, bz - 1)] * g_N_s; w_source = w + ib; }                                      // got it from the diagonal block
		break; 
	      default:
		w_dest = NULL;
		w_source = NULL;
	      }
	      memcpy(w_dest, w_source, g_N_s * sizeof(_Complex double));
	    }
	  }
	}
      }
    }
  }
  free(w_buf);
#endif
  return;
}



void little_D(_Complex double * v, _Complex double *w) {
  int i, j, sq = g_N_s*g_N_s;
  CONE = 1.0;
  CMONE = -1.0;
  CZERO = 0.0;

  if(dfl_subspace_updated) {
    compute_little_D();
    dfl_subspace_updated = 0;
  }
  
#ifdef MPI
  /*init_little_field_exchange(w);*/
  little_field_gather(w);
#endif
  
  /* all the mpilocal stuff first */
  for(i = 0; i < nb_blocks; i++) {
    /* diagonal term */
    _FT(zgemv)("N", &g_N_s, &g_N_s, &CONE, block_list[i].little_dirac_operator,
               &g_N_s, w + i * g_N_s, &ONE, &CZERO, v + i * g_N_s, &ONE, 1);
    
    /* offdiagonal terms */
    for(j = 1; j < 9; j++) {
      _FT(zgemv)("N", &g_N_s, &g_N_s, &CONE, block_list[i].little_dirac_operator + j * sq,
		 &g_N_s, w + (nb_blocks * j + i) * g_N_s, &ONE, &CONE, v + i * g_N_s, &ONE, 1);
    }
  }
  return;
}


void little_D_sym(_Complex double * v, _Complex double *w) {
  
  _Complex double* tmpc1, * tmpc2, * tmpc3;
  tmpc1 = calloc(nb_blocks * 9 * g_N_s, sizeof(_Complex double));
  tmpc2 = calloc(nb_blocks * 9 * g_N_s, sizeof(_Complex double));
  tmpc3 = calloc(nb_blocks * 9 * g_N_s, sizeof(_Complex double));
  
  if(dfl_subspace_updated) {
    compute_little_D();
    dfl_subspace_updated = 0;
  }
  
  little_D_hop(0,tmpc1, w);
  little_D_ee_inv(tmpc2,tmpc1);
  little_D_hop(1,tmpc3, tmpc2);
  little_Dhat_lhs(v, w,tmpc3);
  
  free(tmpc1);
  free(tmpc2);
  free(tmpc3);
  return;
}


void little_D_ee_inv(_Complex double * v, _Complex double *w) {
  int i;
  CONE = 1.0;
  CMONE = -1.0;
  CZERO = 0.0;
  
  for(i = 0; i < nb_blocks/2; i++) {
    _FT(zgemv)("N", &g_N_s, &g_N_s, &CONE, block_list[i].little_dirac_operator_eo,
               &g_N_s, w + i * g_N_s, &ONE, &CZERO, v + i * g_N_s, &ONE, 1);
  }
  return;
}


void little_D_hop(int eo,_Complex double * v, _Complex double *w) {
  int i, j, i_eo,sq = g_N_s*g_N_s;
  CONE = 1.0;
  CMONE = -1.0;
  CZERO = 0.0;

  i_eo=(eo+1)%2;
  
#ifdef MPI
  /*init_little_field_exchange(w);*/
  little_field_gather_eo(eo,w+i_eo*nb_blocks*g_N_s/2);
#endif
  
  for(i = 0; i < nb_blocks/2; i++) {
    for(j = 1; j < 9; j++) {
      _FT(zgemv)("N", &g_N_s, &g_N_s, &CONE, block_list[eo*(nb_blocks/2)+i].little_dirac_operator_eo + j * sq,
		 &g_N_s, w + (nb_blocks * j + (nb_blocks/2)*i_eo+i) * g_N_s, &ONE, &CONE, v + (eo*nb_blocks/2+i) * g_N_s, &ONE, 1);
    } 
  }
  return;
}

void little_Dhat_lhs(_Complex double * v, _Complex double *w, _Complex double *u) {
  int i,j;
  CONE = 1.0;
  CMONE = -1.0;
  CZERO = 0.0;


  for(i = nb_blocks/2; i < nb_blocks; i++) {
    _FT(zgemv)("N", &g_N_s, &g_N_s, &CONE, block_list[i].little_dirac_operator_eo,
               &g_N_s, w + i * g_N_s, &ONE, &CZERO, v + i * g_N_s, &ONE, 1);
  }
  
  for (i=nb_blocks/2; i < nb_blocks; i++) {
    for (j=0;j<g_N_s;j++) {
      *(v+ i * g_N_s+ j) = *(v+ i * g_N_s+ j) - *(u+ i * g_N_s+ j);
    }
  }
  return;
}



void little_Dhat_rhs(int eo, _Complex double * v, double r, _Complex double *w) {
  int i, j;
  
  for(i = 0; i < nb_blocks/2; i++) {
    for (j=0;j<g_N_s;j++) {
      *(v+eo*nb_blocks*g_N_s/2+i*g_N_s+j) = *(w+eo*nb_blocks*g_N_s/2+i*g_N_s+j) + r * *(v+eo*nb_blocks*g_N_s/2+i*g_N_s+j);
    }
  }
  return;
}


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



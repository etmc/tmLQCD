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
#include "complex.h"
#include "block.h"
#include "linalg/blas.h"
#include "solver/gcr4complex.h"
#include "solver/generate_dfl_subspace.h"
#include "block.h"
#include "linalg_eo.h"
#include "little_D.h"


/* assume we have a little field w                       */
/* which has length 9*no_blocks*N_s                      */
/* with usual order in space                             */
/* no_blocks = 2 currently fixed                         */
/* and blocks devide z-direction by 2                    */
/*                                                       */
/* block[0], block[1], block[0], block[1], block[0]  ... */
/* local             , +t                , -t        ... */
/*                                                       */
/* block[0], block[1], block[0], block[1]                */
/* +z                , -z                                */
/* wasting some memory here...                           */

const int no_blocks = 2;
const int nblks_t = 1;
const int nblks_x = 1;
const int nblks_y = 1;
const int nblks_z = 2;
int nblks_dir[4] = {1,1,1,2};
int dfl_subspace_updated = 1;

/* some lapack related stuff */
static int ONE = 1;
static complex CONE, CZERO, CMONE;

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

void init_little_field_exchange(complex * w);
void wait_little_field_exchange(const int mu);

void unit_little_D(complex *v, complex *w) {
  memcpy(v, w, 2*g_N_s*sizeof(complex));

  return;
}

/** ANOTHER TESTING FUNCTION */
void invert_little_D_spinor(spinor *r, spinor *s){
  int j;
  spinor *psi[2];
  complex *v, *w;

  v = calloc(2 * 9 * g_N_s, sizeof(complex));
  w = calloc(2 * 9 * g_N_s, sizeof(complex));
  psi[0] = calloc(VOLUME, sizeof(spinor));
  psi[1] = psi[0] + VOLUME / 2;
  split_global_field(psi[0], psi[1], s);

  for (j = 0; j < g_N_s; ++j) {/*loop over block.basis */
/*     v[j]         = block_scalar_prod(psi[0], block_list[0].basis[j], VOLUME/2); */
/*     v[j + g_N_s] = block_scalar_prod(psi[1], block_list[1].basis[j], VOLUME/2); */
    v[j]         = scalar_prod(block_list[0].basis[j], psi[0], VOLUME/2, 0);
    v[j + g_N_s] = scalar_prod(block_list[1].basis[j], psi[1], VOLUME/2, 0);
  }

  gcr4complex(w, v, 10, 100, 1e-31, 1, 2 * g_N_s, 1, 2 * 9 * g_N_s, &little_D);

  mul(psi[0], w[0], block_list[0].basis[0], VOLUME/2);
  mul(psi[1], w[g_N_s], block_list[1].basis[0], VOLUME/2);
  for(j = 1; j < g_N_s; ++j) {
    assign_add_mul(psi[0], block_list[0].basis[j], w[j], VOLUME/2);
    assign_add_mul(psi[1], block_list[1].basis[j], w[g_N_s+j], VOLUME/2);
  }
  reconstruct_global_field(r, psi[0], psi[1]);

  free(v);
  free(w);
  free(psi[0]);
}


void project2(spinor * const out, spinor * const in);

/** ANOTHER TESTING FUNCTION */
void apply_little_D_spinor(spinor *r, spinor *s){
  int j, k;
  spinor *psi[2];
  complex *v, *w;

  v = calloc(2 * 9 * g_N_s, sizeof(complex));
  w = calloc(2 * 9 * g_N_s, sizeof(complex));
  psi[0] = calloc(VOLUME, sizeof(spinor));
  psi[1] = psi[0] + VOLUME / 2;
  split_global_field(psi[0], psi[1], s);

  for (j = 0; j < g_N_s; ++j) {
/*     v[j]         = block_scalar_prod(psi[0], block_list[0].basis[j], VOLUME/2); */
    v[j]         = scalar_prod(block_list[0].basis[j], psi[0], VOLUME/2, 0);
/*     v[j + g_N_s] = block_scalar_prod(psi[1], block_list[1].basis[j], VOLUME/2); */
    v[j + g_N_s] = scalar_prod(block_list[1].basis[j], psi[1], VOLUME/2, 0);
  }

  if (g_debug_level > 2){
    if (!g_cart_id){
      for (j = 0; j < 2* g_N_s; ++j) {
        printf("LITTLE_D for 0: v[%u] = %1.5e + %1.5e i\n", j, v[j].re, v[j].im);
      }
    }
#ifdef MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
  }

  if (g_debug_level > 4){
    for (k = 1; k < 16; ++k){
      if (g_cart_id == k){
        for (j = 0; j < 2* g_N_s; ++j) {
          printf("LITTLE_D for %u: v[%u] = %1.5e + %1.5e i\n", k, j, v[j].re, v[j].im);
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
      for (j = 0; j < 2 * g_N_s; ++j) {
        printf("LITTLE_D for 0: w[%u] = %1.5e + %1.5e i\n", j, w[j].re, w[j].im);
      }
    }
#ifdef MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
  }

  if (g_debug_level > 4)
  {
    for (k = 1; k < 16; ++k){
      if (g_cart_id == k){
        for (j = 0; j < 2* g_N_s; ++j) {
          printf("LITTLE_D for %u: w[%u] = %1.5e + %1.5e i\n", k, j, w[j].re, w[j].im);
        }
      }
#ifdef MPI
      MPI_Barrier(MPI_COMM_WORLD);
#endif
    }
  }

  mul(psi[0], w[0], block_list[0].basis[0], VOLUME/2);
  mul(psi[1], w[g_N_s], block_list[1].basis[0], VOLUME/2);
  for(j = 1; j < g_N_s; ++j) {
    assign_add_mul(psi[0], block_list[0].basis[j], w[j], VOLUME/2);
    assign_add_mul(psi[1], block_list[1].basis[j], w[g_N_s+j], VOLUME/2);
  }
  reconstruct_global_field(r, psi[0], psi[1]);

  free(v);
  free(w);
  free(psi[0]);
}

void alt_little_field_gather(complex * w) {
#ifdef MPI
  MPI_Status status;
  int size = 25 * g_N_s * sizeof(complex);
  complex *buf = malloc(size);
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
  memcpy(w + 17 * g_N_s, w, g_N_s * sizeof(complex));

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
  memcpy(w + 14 * g_N_s, w + g_N_s, g_N_s * sizeof(complex));

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


void little_field_gather(complex * w) {
#ifdef MPI
  int err;
  
  /* LOWER BLOCK */

  /* Send t up */
  MPI_Isend(w, 2*g_N_s, MPI_DOUBLE_COMPLEX, g_nb_t_up, T_UP, g_cart_grid, &lrequests[0]);
  MPI_Irecv(w + 4 * g_N_s, 2*g_N_s, MPI_DOUBLE_COMPLEX, g_nb_t_dn, T_UP, g_cart_grid, &lrequests[1]);

  /* Send t down */
  MPI_Isend(w, 2*g_N_s, MPI_DOUBLE_COMPLEX, g_nb_t_dn, T_DN, g_cart_grid, &lrequests[2]);
  MPI_Irecv(w + 2 * g_N_s, 2*g_N_s, MPI_DOUBLE_COMPLEX, g_nb_t_up, T_DN, g_cart_grid, &lrequests[3]);

  /* Send x up */
  MPI_Isend(w, 2*g_N_s, MPI_DOUBLE_COMPLEX, g_nb_x_up, X_UP, g_cart_grid, &lrequests[4]);
  MPI_Irecv(w + 8 * g_N_s, 2*g_N_s, MPI_DOUBLE_COMPLEX, g_nb_x_dn, X_UP, g_cart_grid, &lrequests[5]);

  /* Send x down */
  MPI_Isend(w, 2*g_N_s, MPI_DOUBLE_COMPLEX, g_nb_x_dn, X_DN, g_cart_grid, &lrequests[6]);
  MPI_Irecv(w + 6 * g_N_s, 2*g_N_s, MPI_DOUBLE_COMPLEX, g_nb_x_up, X_DN, g_cart_grid, &lrequests[7]);

  /* Send y up */
  MPI_Isend(w, 2*g_N_s, MPI_DOUBLE_COMPLEX, g_nb_y_up, Y_UP, g_cart_grid, &lrequests[8]);
  MPI_Irecv(w + 12 * g_N_s, 2*g_N_s, MPI_DOUBLE_COMPLEX, g_nb_y_dn, Y_UP, g_cart_grid, &lrequests[9]);

  /* Send y down */
  MPI_Isend(w, 2*g_N_s, MPI_DOUBLE_COMPLEX, g_nb_y_dn, Y_DN, g_cart_grid, &lrequests[10]);
  MPI_Irecv(w + 10 * g_N_s, 2*g_N_s, MPI_DOUBLE_COMPLEX, g_nb_y_up, Y_DN, g_cart_grid, &lrequests[11]);


  /* Send z down */
  MPI_Isend(w, g_N_s, MPI_DOUBLE_COMPLEX, g_nb_z_dn, Z_DN, g_cart_grid, &lrequests[12]);
  MPI_Irecv(w + 15 * g_N_s, g_N_s, MPI_DOUBLE_COMPLEX, g_nb_z_up, Z_DN, g_cart_grid, &lrequests[13]);


  /* Send z up */
  MPI_Isend(w + g_N_s, g_N_s, MPI_DOUBLE_COMPLEX, g_nb_z_up, Z_UP, g_cart_grid, &lrequests[14]);
  MPI_Irecv(w + 16 * g_N_s, g_N_s, MPI_DOUBLE_COMPLEX, g_nb_z_dn, Z_UP, g_cart_grid, &lrequests[15]);

  /* Send z up */
  memcpy(w + 17 * g_N_s, w, g_N_s * sizeof(complex));

  /* Send z down */
  memcpy(w + 14 * g_N_s, w + g_N_s, g_N_s * sizeof(complex));

  err = MPI_Waitall(16, lrequests, lstatus);
#endif
  return;
}

void little_D(complex * v, complex *w) {
  int i, j, sq = g_N_s*g_N_s;
  CONE.re = 1.;
  CONE.im = 0.;
  CMONE.re = -1.;
  CMONE.im = 0.;
  CZERO.re = 0.;
  CZERO.im = 0.;

  if(dfl_subspace_updated) {
    compute_little_D_diagonal();
    compute_little_D_offdiagonal();
    dfl_subspace_updated = 0;
  }

#ifdef MPI
  /*init_little_field_exchange(w);*/
  little_field_gather(w);
#endif

  /* all the mpilocal stuff first */
  for(i = 0; i < no_blocks; i++) {
    /* diagonal term */
    _FT(zgemv)("N", &g_N_s, &g_N_s, &CONE, block_list[i].little_dirac_operator,
               &g_N_s, w + i * g_N_s, &ONE, &CZERO, v + i * g_N_s, &ONE, 1);

    /* offdiagonal terms */
    for(j = 1; j < 9; j++) {
      _FT(zgemv)("N", &g_N_s, &g_N_s, &CONE, block_list[i].little_dirac_operator + j * sq,
		 &g_N_s, w + (2 * j + i) * g_N_s, &ONE, &CONE, v + i * g_N_s, &ONE, 1);
    }
  }
  return;
}


void init_little_field_exchange(complex * w) {
#ifdef MPI
  int i = 0;
#  ifdef PARALLELT
  int no_dirs = 2;
#  elif defined PARALLELXT
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
  for(i = 0; i < no_dirs; i+=2) {
    /* send to the right, receive from the left */
    MPI_Isend((void*)w, no_blocks*g_N_s, MPI_DOUBLE_COMPLEX, g_nb_list[i], 
              i, g_cart_grid, &lrequests[2*i]);
    MPI_Irecv((void*)(w + no_blocks*(i+2)*g_N_s), no_blocks*g_N_s, MPI_DOUBLE_COMPLEX, g_nb_list[i+1], 
              i, g_cart_grid, &lrequests[2*i+1]);
    
    /* send to the left, receive from the right */
    MPI_Isend((void*)w, no_blocks*g_N_s, MPI_DOUBLE_COMPLEX, g_nb_list[i+1], 
              i+1, g_cart_grid, &lrequests[2*i+2]);
    MPI_Irecv((void*)(w + no_blocks*(i+1)*g_N_s), no_blocks*g_N_s, MPI_DOUBLE_COMPLEX, g_nb_list[i], 
              i+1, g_cart_grid, &lrequests[2*i+3]);
    waitcount += 4;
  }
#  ifdef PARALLELXYZT
  /* send to the right, receive from the left */
  i = 6;
  MPI_Isend((void*)(w + g_N_s), g_N_s, MPI_DOUBLE_COMPLEX, g_nb_list[i], 
            i, g_cart_grid, &lrequests[2*i]);
  MPI_Irecv((void*)(w + (no_blocks*(i+1)+1)*g_N_s), g_N_s, MPI_DOUBLE_COMPLEX, g_nb_list[i+1], 
            i, g_cart_grid, &lrequests[2*i+1]);
  
  /* send to the left, receive from the right */
  MPI_Isend((void*)w, g_N_s, MPI_DOUBLE_COMPLEX, g_nb_list[i+1], 
            i+1, g_cart_grid, &lrequests[2*i+2]);
  MPI_Irecv((void*)(w + no_blocks*(i+1)*g_N_s), g_N_s, MPI_DOUBLE_COMPLEX, g_nb_list[i], 
            i+1, g_cart_grid, &lrequests[2*i+3]);
  waitcount += 4;
#  endif
#endif
  return;
}

void wait_little_field_exchange(const int mu) {
  int err;
#ifdef MPI
  err = MPI_Waitall(2, &lrequests[2*mu], &lstatus[2*mu]);
  waitcount -= 2;
#endif
  return;
}



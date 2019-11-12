/***********************************************************************
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008 Carsten Urbach
 *               2010 claude Tadonki
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
# include<tmlqcd_config.h>
#endif
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<string.h>
#include"global.h"
#include"su3.h"
#include"linalg_eo.h"
#include"mcr4complex.h"
#include "time.h"

static void init_lmcr(const int _M, const int _V);
static void free_lmcr();
static _Complex double * chi = NULL;
static _Complex double * xi = NULL;
static _Complex double * Achi = NULL;
static _Complex double * Axi = NULL;
static _Complex double * tmp = NULL;
static int lmcr_init = 0;

int mcr4complex(_Complex double * const P, _Complex double * const Q, 
                const int m, const int max_restarts,
                const double eps_sq, const int rel_prec,
                const int N, const int parallel, 
                const int lda, c_matrix_mult f) {

  int k, l, restart, i, p=0;
  double norm_sq, norm,err;
  _Complex double ctmp;
  _Complex double alpha,beta;
  _Complex double one = 1.0;


  double atime, etime;
  init_lmcr(m, lda);

  norm_sq = lsquare_norm(Q, N, parallel);
  if(norm_sq < 1.e-20) {
    norm_sq = 1.;
  }

#ifdef TM_USE_MPI
  atime = MPI_Wtime();
#else
  atime = ((double)clock())/((double)(CLOCKS_PER_SEC));
#endif


  for(restart = 0; restart < max_restarts; restart++) {

    f(tmp, P);
    ldiff(chi, Q, tmp, N);
    memcpy(xi, chi, N*sizeof(_Complex double));

    f(Axi,xi);

    err = lsquare_norm(chi, N, parallel);
    if(g_proc_id == g_stdio_proc && g_debug_level > 1){
      printf("lPHCR: %d\t%g true residue\n", p, err); 
      fflush(stdout);
    }

    if(((err <= eps_sq) && (rel_prec == 0)) || ((err <= eps_sq * norm_sq) && (rel_prec == 1))) {
      return (p);
    }

    for(k = 0;k<m ; k++) {
      alpha = lscalar_prod(Axi, chi, N, parallel);
      norm = lsquare_norm(Axi, N, parallel);

      alpha /= norm;

      lassign_add_mul(P,xi, alpha, N);
      lassign_diff_mul(chi, Axi, alpha, N);

      err = lsquare_norm(chi, N, parallel);
      p++;

      if(g_proc_id == g_stdio_proc && g_debug_level > 1){
        printf("mCR: %d\t%g iterated residue\n", p, err); 
        fflush(stdout);
      }
      /* Precision reached? */
      if((k == m-1) || ((err <= eps_sq) && (rel_prec == 0)) || ((err <= eps_sq*norm_sq) && (rel_prec == 1))) {
        break;
      }

      f(Achi, chi);

      beta = lscalar_prod(Axi, Achi, N, 1);
      beta /= -norm;

      lmul_add_mul(tmp,chi,xi,one,beta,N);
      memcpy(xi, tmp, N*sizeof(_Complex double));

      lmul_add_mul(tmp,Achi,Axi,one,beta,N);
      memcpy(Axi, tmp, N*sizeof(_Complex double));

    }

  }

  /* check if it converges in the last restart cycle */
  if (restart == max_restarts) {
    f(tmp, P);
    ldiff(chi, Q, tmp, N);
    memcpy(xi, chi, N*sizeof(_Complex double));

    f(Axi,xi);

    err = lsquare_norm(chi, N, parallel);

#ifdef TM_USE_MPI
    etime = MPI_Wtime();
#else
    etime = ((double)clock())/((double)(CLOCKS_PER_SEC));
#endif
    if(g_proc_id == g_stdio_proc && g_debug_level > 1){
      printf("lPHCR: %d\t%g true residue, time spent %f s\n", p, err, (etime - atime)); 
      fflush(stdout);
    }
    if(((err <= eps_sq) && (rel_prec == 0)) || ((err <= eps_sq * norm_sq) && (rel_prec == 1))) {
      return (p);
    }
  }


  if(g_proc_id == 0 && g_debug_level > 1) printf("lPHcr: for -1 %d %e\n", p, err);
  return(-1);
}

static void init_lmcr(const int _M, const int _V){
  static int Vo = -1;
  static int M = -1;

  int i;
  if((M != _M)||(lmcr_init == 0)||(Vo != _V)){
    if(lmcr_init == 1) free_lmcr();
    Vo = _V;
    M = _M;
    chi = calloc(Vo, sizeof(_Complex double));
    xi = calloc(Vo, sizeof(_Complex double));
    Achi = calloc(Vo, sizeof(_Complex double));
    Axi = calloc(Vo, sizeof(_Complex double));
    tmp = calloc(Vo, sizeof(_Complex double));

    lmcr_init = 1;
  }
}

static void free_lmcr() 
{
  lmcr_init = 0;
  free(chi);
  free(xi);
  free(Achi);
  free(Axi);
  free(tmp);
  return;
}



void lmul_add_mul(_Complex double * const R, _Complex double * const S, _Complex double * const T,const _Complex double c, const _Complex double d, const int N) 
{
  int i;
  for(i = 0; i < N; i++) {
    R[i] = c * S[i] + d * T[i];
  }
  return;
}


void lmul_diff_mul(_Complex double * const R, _Complex double * const S, _Complex double * const T,const _Complex double c, const _Complex double d, const int N) 
{
  int i;
  for(i = 0; i < N; i++) {
    R[i] = c * S[i] - d * T[i];
  }
  return;
}


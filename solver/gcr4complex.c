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
# include<config.h>
#endif
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<string.h>
#include"global.h"
#include"su3.h"
#include"linalg_eo.h"
#include"gcr4complex.h"

static void init_lgcr(const int _M, const int _V);
static void free_lgcr();
static _Complex double ** a = NULL; 
static _Complex double * _a = NULL;
static double * b = NULL;
static _Complex double * c = NULL;
static _Complex double ** chi = NULL;
static _Complex double * _chi = NULL;
static _Complex double ** xi = NULL;
static _Complex double * _xi = NULL;
static _Complex double * alpha = NULL;
static _Complex double * tmp = NULL;
static _Complex double * rho = NULL;
static int lgcr_init = 0;

int gcr4complex(_Complex double * const P, _Complex double * const Q, 
		const int m, const int max_restarts,
		const double eps_sq, const int rel_prec,
		const int N, const int parallel, 
		const int lda, c_matrix_mult f) {
  
  int k, l, restart, i, p=0;
  double norm_sq, err;
  _Complex double ctmp;

  init_lgcr(m, lda);

  norm_sq = lsquare_norm(Q, N, parallel);
  if(norm_sq < 1.e-20) {
    norm_sq = 1.;
  }
  for(restart = 0; restart < max_restarts; restart++) {
    f(tmp, P);
    ldiff(rho, Q, tmp, N);
    err = lsquare_norm(rho, N, parallel);
    if(g_proc_id == g_stdio_proc && g_debug_level > 1){/*CT: was "g_debug_level > 0" */
      printf("lGCR: %d\t%g true residue %1.3e\n", restart * m, err, norm_sq); 
      fflush(stdout);
    }
    if(((err <= eps_sq) && (rel_prec == 0)) || ((err <= eps_sq * norm_sq) && (rel_prec == 1))) {
      if(g_proc_id == 0 && g_debug_level > 1) printf("lgcr: %d %e %e %e %e\n", p, err, norm_sq, err/norm_sq, eps_sq);
      return (p);
    }
    for(k = 0; ; k++) {
      memcpy(xi[k], rho, N*sizeof(_Complex double));
      /* here we could put in a preconditioner */
      f(tmp, xi[k]); 
      /* tmp will become chi[k] */
      for(l = 0; l < k; l++) {
        a[l][k] = lscalar_prod(chi[l], tmp, N, parallel);
        lassign_diff_mul(tmp, chi[l], a[l][k], N);
      }
      b[k] = sqrt(lsquare_norm(tmp, N, parallel));
      lmul_r(chi[k], 1./b[k], tmp, N);
      c[k] = lscalar_prod(chi[k], rho, N, parallel);
      lassign_diff_mul(rho, chi[k], c[k], N);
      err = lsquare_norm(rho, N, parallel);
      if(g_proc_id == g_stdio_proc && g_debug_level > 1){
        printf("lGCR: %d\t%g iterated residue\n", restart*m+k, err); 
        fflush(stdout);
      }
      p++;
      /* Precision reached? */
      if((k == m-1) || ((err <= eps_sq) && (rel_prec == 0)) || ((err <= eps_sq*norm_sq) && (rel_prec == 1))) {
	break;
      }
    }
    /* prepare for restart */
    c[k] /= b[k];
    lassign_add_mul(P, xi[k], c[k], N);
    for(l = k-1; l >= 0; --l)
    {
      for(i = l+1; i <= k; ++i)
      {
        ctmp  = a[l][i] * c[i];
        c[l] -= ctmp;
      }
      c[l] /= b[l];
      lassign_add_mul(P, xi[l], c[l], N);
    }
  }
  if(g_proc_id == 0 && g_debug_level > 1) printf("lgcr: for -1 %d %e %e %e %e\n", p, err, norm_sq, err/norm_sq, eps_sq);
  return(-1);
}

static void init_lgcr(const int _M, const int _V){
  static int Vo = -1;
  static int M = -1;

  int i;
  if((M != _M)||(lgcr_init == 0)||(Vo != _V)){
    if(lgcr_init == 1) free_lgcr();
    Vo = _V;
    M = _M;
    a = calloc(M+1, sizeof(_Complex double *));
    chi = calloc(M, sizeof(_Complex double *));
    xi = calloc(M, sizeof(_Complex double *));
    tmp = calloc(Vo, sizeof(_Complex double));
    rho = calloc(Vo, sizeof(_Complex double));
    _a = calloc((M+1)*M, sizeof(_Complex double));
    a[0] = _a;
    _chi = calloc(M*Vo, sizeof(_Complex double));
    chi[0] = _chi;
    _xi = calloc(M*Vo, sizeof(_Complex double));
    xi[0] = _xi;

    b = calloc(M, sizeof(double));
    c = calloc(M, sizeof(_Complex double));
    alpha = calloc(M+1, sizeof(_Complex double));
    for(i = 1; i < M; i++) { 
      chi[i] = chi[i-1] + Vo;
      xi[i] = xi[i-1] + Vo;
      a[i] = a[i-1] + M;
    }
    a[M] = a[M-1] + M;
    lgcr_init = 1;
  }
}

static void free_lgcr() 
{
  lgcr_init = 0;
  free(a);
  free(chi);
  free(_a);
  free(_chi);
  free(alpha);
  free(c);
  free(_xi);
  free(xi);
  free(rho);
  free(tmp);
  return;
}


void ldiff(_Complex double * const Q, _Complex double * const R, _Complex double * const S, const int N) 
{
  for(int i = 0; i < N; ++i)
    Q[i] = R[i] - S[i];
  return;
}

void ldiff_assign(_Complex double * const Q, _Complex double * const S, const int N) 
{
  for(int i = 0; i < N; ++i)
    Q[i] -= S[i];
  return;
}

void ladd(_Complex double * const Q, _Complex double * const R, _Complex double * const S, const int N) 
{
  for(int i = 0; i < N; ++i)
    Q[i] = R[i] + S[i];
  return;
}

void ladd_assign(_Complex double * const Q, _Complex double * const S, const int N) 
{
  for(int i = 0; i < N; ++i)
    Q[i] += S[i];
  return;
}

double lsquare_norm(_Complex double * const Q, const int N, const int parallel) 
{
  double nrm = 0.0;

  for(int i = 0; i < N; ++i)
    
    nrm += conj(Q[i]) * Q[i];
#ifdef MPI
  if(parallel)
  {
    double nrm2 = nrm;
    MPI_Allreduce(&nrm2, &nrm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  }
#endif

  return(nrm);
}

_Complex double lscalar_prod(_Complex double * const R, _Complex double * const S, const int N, const int parallel) 
{
  _Complex double res = 0.0;

  for(int i = 0; i < N; ++i)
    res += conj(R[i]) * S[i];
  
#ifdef MPI
  if(parallel)
  {
    _Complex double res2 = res;
    MPI_Allreduce(&res2, &res, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
  }
#endif

  return(res);
}

void lmul_r(_Complex double * const R, const double c, _Complex double * const S, const int N) 
{
  for(int i = 0; i < N; ++i)
    R[i] = c * S[i];
}

void lmul(_Complex double * const R, const _Complex double c, _Complex double * const S, const int N) 
{
  for(int i = 0; i < N; ++i)
    R[i] = c * S[i];
}

void lassign_add_mul(_Complex double * const R, _Complex double * const S, const _Complex double c, const int N)
{
  for(int i = 0; i < N; ++i)
    R[i] += c * S[i];
}

void lassign_diff_mul(_Complex double * const R, _Complex double * const S, const _Complex double c, const int N) 
{
  for(int i = 0; i < N; i++)
    R[i] -= c * S[i];
}

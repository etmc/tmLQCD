/***********************************************************************
 *
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008 Carsten Urbach
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
 *
 * Generalized minimal residual (GMRES) with deflated restarting (Morgan)
 *
 * This requires LAPACK to run...
 *
 * Inout:                                                                      
 *  spinor * P       : guess for the solving spinor                                             
 * Input:                                                                      
 *  spinor * Q       : source spinor
 *  int m            : Maximal dimension of Krylov subspace                                     
 *  int nr_ev        : number of eigenvectors to be deflated
 *  int max_restarts : maximal number of restarts                                   
 *  double eps       : stopping criterium                                                     
 *  matrix_mult f    : pointer to a function containing the matrix mult
 *                     for type matrix_mult see matrix_mult_typedef.h
 *
 * Autor: Carsten Urbach <urbach@ifh.de>
 ********************************************************************************/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include"global.h"
#include <complex.h>
#include"su3.h"
#include"linalg_eo.h"
#include"diagonalise_general_matrix.h"
#include"quicksort.h"
#include"linalg/lapack.h"
#include"linalg/blas.h"
#include"solver/gram-schmidt.h"
#include"solver/gmres.h"
#include "solver/solver_field.h"
#include"gmres_dr.h"

#ifndef HAVE_LAPACK
/* In case there is no lapack use normal gmres */
int gmres_dr(spinor * const P,spinor * const Q, 
	  const int m, const int nr_ev, const int max_restarts,
	  const double eps_sq, const int rel_prec,
	  const int N, matrix_mult f){
  return(gmres(P, Q, m, max_restarts, eps_sq, rel_prec, N, 1, f));
}

#else

static void init_gmres_dr(const int _M, const int _V);
_Complex double short_scalar_prod(_Complex double * const x, _Complex double * const y, const int N);
void short_ModifiedGS(_Complex double v[], int n, int m, _Complex double A[], int lda);

static _Complex double ** work;
static _Complex double * _work;
static _Complex double ** work2;
static _Complex double * _work2;
static _Complex double ** H;
static _Complex double ** G;
static _Complex double * alpha;
static _Complex double * c;
static double * s;
static spinor ** V;
static spinor * _v;
static spinor ** Z;
static spinor * _z;
static _Complex double * _h;
static _Complex double * _g;
static _Complex double * alpha;
static _Complex double * c;
static double * s;
static _Complex double * evalues;
static double * sortarray;
static int * idx;
static int one = 1;
static _Complex double cmone;
static _Complex double cpone;
static _Complex double czero;

int gmres_dr(spinor * const P,spinor * const Q, 
	  const int m, const int nr_ev, const int max_restarts,
	  const double eps_sq, const int rel_prec,
	  const int N, matrix_mult f){

  int restart=0, i, j, k, l;
  double beta, eps, norm, beta2=0.;
  _Complex double *lswork = NULL;
  int lwork;
  _Complex double tmp1, tmp2;
  int info=0;
  int _m = m, mp1 = m+1, np1 = nr_ev+1, ne = nr_ev, V2 = 12*(VOLUMEPLUSRAND)/2, _N = 12*N;
  spinor ** solver_field = NULL;
  const int nr_sf = 3;

  if(N == VOLUME) {
    init_solver_field(&solver_field, VOLUMEPLUSRAND, nr_sf);
  }
  else {
    init_solver_field(&solver_field, VOLUMEPLUSRAND/2, nr_sf);
  }
  double err=0.;
  spinor * r0, * x0;

  cmone = -1.;
  cpone = 1.;
  czero = 0.;
  
  r0 = solver_field[0];
  x0 = solver_field[2];
  eps=sqrt(eps_sq);  
  init_gmres_dr(m, (VOLUMEPLUSRAND));
  norm = sqrt(square_norm(Q, N, 1));

  assign(x0, P, N);

  /* first normal GMRES cycle */
  /* r_0=Q-AP  (b=Q, x+0=P) */
  f(r0, x0);
  diff(r0, Q, r0, N);
  
  /* v_0=r_0/||r_0|| */
  alpha[0] = sqrt(square_norm(r0, N, 1));
  err = alpha[0];
  
  if(g_proc_id == g_stdio_proc && g_debug_level > 0){
    printf("%d\t%e true residue\n", restart * m, creal(alpha[0])*creal(alpha[0])); 
    fflush(stdout);
  }
  
  if(creal(alpha[0])==0.){
    assign(P, x0, N);
    finalize_solver(solver_field, nr_sf);
    return(restart*m);
  }
  
  mul_r(V[0], 1./creal(alpha[0]), r0, N);
  
  for(j = 0; j < m; j++){
    /* solver_field[0]=A*v_j */

    /* Set h_ij and omega_j */
    /* solver_field[1] <- omega_j */    
    f(solver_field[1], V[j]);
/*     assign(solver_field[1], solver_field[0], N); */
    for(i = 0; i <= j; i++){
      H[i][j] = scalar_prod(V[i], solver_field[1], N, 1);
      /* G, work and work2 are in Fortran storage: columns first */
      G[j][i] = H[i][j];
      work2[j][i] = H[i][j];
      work[i][j] = conj(H[i][j]);
      assign_diff_mul(solver_field[1], V[i], H[i][j], N);
    }
    
    H[j+1][j] = sqrt(square_norm(solver_field[1], N, 1));
    G[j][j+1] = H[j+1][j];
    work2[j][j+1] = H[j+1][j];
    work[j+1][j] = conj(H[j+1][j]);
    beta2 = creal(H[j+1][j])*creal(H[j+1][j]); 
    for(i = 0; i < j; i++){
      tmp1 = H[i][j];
      tmp2 = H[i+1][j];
      H[i][j] = (tmp2) * (s[i]);
      H[i][j] += conj(c[i]) * (tmp1);
      H[i+1][j] = (tmp1) * (s[i]);
      H[i+1][j] -= (c[i]) * (tmp2);
    }
    
    /* Set beta, s, c, alpha[j],[j+1] */
    beta = sqrt(creal(H[j][j] * conj(H[j][j])) + creal(H[j+1][j] * conj(H[j+1][j])));
    s[j] = creal(H[j+1][j]) / beta;
    c[j] = (H[j][j]) / beta;
    H[j][j] = beta;
    alpha[j+1] = (alpha[j]) * (s[j]);
    tmp1 = alpha[j];
    alpha[j] = conj(c[j]) * (tmp1);
    
    /* precision reached? */
    if(g_proc_id == g_stdio_proc && g_debug_level > 0){
      printf("%d\t%e residue\n", restart*m+j, creal(alpha[j+1])*creal(alpha[j+1])); 
      fflush(stdout);
    }
    if(((creal(alpha[j+1]) <= eps) && (rel_prec == 0)) || ((creal(alpha[j+1]) <= eps*norm) && (rel_prec == 1))){
      alpha[j] = (alpha[j]) * (1./creal(H[j][j]));
      assign_add_mul(x0, V[j], alpha[j], N);
      for(i = j-1; i >= 0; i--){
	for(k = i+1; k <= j; k++){
	  tmp1 = (H[i][k]) * (alpha[k]); 
	  /* alpha[i] -= tmp1 */
	  alpha[i] -= tmp1;
	}
	alpha[i] = (alpha[i]) * (1./creal(H[i][i]));
	assign_add_mul(x0, V[i], alpha[i], N);
      }
      for(i = 0; i < m; i++){
	alpha[i] = creal(alpha[i]);
      }
      assign(P, x0, N);
      finalize_solver(solver_field, nr_sf);
      return(restart*m+j);
    }
    /* if not */
    else {
      mul_r(V[(j+1)], 1./creal(H[j+1][j]), solver_field[1], N); 
    }
    
  }
  j=m-1;
  /* prepare for restart */
  alpha[j] = (alpha[j]) * (1./creal(H[j][j]));
  assign_add_mul(x0, V[j], alpha[j], N);
  if(g_proc_id == 0 && g_debug_level > 3) {
    printf("alpha: %e %e\n", creal(alpha[j]), cimag(alpha[j]));
  }
  for(i = j-1; i >= 0; i--){
    for(k = i+1; k <= j; k++){
      tmp1 = (H[i][k]) * (alpha[k]);
      alpha[i] -= tmp1;
    }
    alpha[i] = (alpha[i]) * (1./creal(H[i][i]));
    if(g_proc_id == 0 && g_debug_level > 3) {
      printf("alpha: %e %e\n", creal(alpha[i]), cimag(alpha[i]));
    }
    assign_add_mul(x0, V[i], alpha[i], N);
  }

  /* This produces c=V_m+1*r0 */
  for(i = 0; i < mp1; i++) {
    c[i] = scalar_prod(V[i], r0, N, 1); 
    if(g_proc_id == 0 && g_debug_level > 3) {
      printf("c: %e %e err = %e\n", creal(c[i]), cimag(c[i]), err);
    }
  }

  for(restart = 1; restart < max_restarts; restart++) {  

    /* compute c-\bar H \alpha */
    _FT(zgemv)("N", &mp1, &_m, &cmone, G[0], &mp1, alpha, &one, &cpone, c, &one, 1);
    err = creal(sqrt(short_scalar_prod(c, c, mp1)));
    if(g_proc_id == 0 && g_debug_level > 0) {
      printf("%d\t %e short residue\n", m*restart, err*err);
    } 
    
    /* Compute new residual r0 */
    /* r_0=Q-AP  (b=Q, x+0=P) */
    if(g_debug_level > 0) {
      f(r0, x0);
      diff(r0, Q, r0, N);
      tmp1 = sqrt(square_norm(r0, N, 1)) * I;
      if(g_proc_id == g_stdio_proc)
      {
	printf("%d\t%e true residue\n", m*restart, cimag(tmp1)*cimag(tmp1)); 
	fflush(stdout);
      }
    }
    mul(r0, c[0], V[0], N);
    for(i = 1; i < mp1; i++) {
      assign_add_mul(r0, V[i], c[i], N);
    } 
    if(g_debug_level > 3) {
      tmp1 = sqrt(square_norm(r0, N, 1)) * I;
      if(g_proc_id == g_stdio_proc){
	printf("%d\t%e residue\n", m*restart, cimag(tmp1)*cimag(tmp1)); 
	fflush(stdout);
      }
    }
    /* Stop if satisfied */
    if(err < eps){
      assign(P, x0, N);
      finalize_solver(solver_field, nr_sf);
      return(restart*m);
    }

    /* Prepare to compute harmonic Ritz pairs */
    for(i = 0; i < m-1; i++){
      alpha[i] = 0.;
    }
    alpha[m-1] = 1.;
    _FT(zgesv)(&_m, &one, work[0], &mp1, idx, alpha, &_m, &info); 
    for(i = 0; i < m; i++) {
      G[m-1][i] += beta2*alpha[idx[i]-1];
    }
    if(g_proc_id == 0 && g_debug_level > 3){
      printf("zgesv returned info = %d, c[m-1]= %e, %e , idx[m-1]=%d\n", 
	     info, creal(alpha[idx[m-1]-1]), cimag(alpha[idx[m-1]-1]), idx[m-1]);
    }
    /* c - \bar H * d -> c */
    /* G contains H + \beta^2 H^-He_n e_n^H */

    /* Compute harmonic Ritz pairs */
    diagonalise_general_matrix(m, G[0], mp1, alpha, evalues);
    for(i = 0; i < m; i++) {
      sortarray[i] = creal(evalues[i] * conj(evalues[i]));
      idx[i] = i;
    }
    quicksort(m, sortarray, idx);
    if(g_proc_id == g_stdio_proc && g_debug_level > 1) {
      for(i = 0; i < m; i++)
	printf("# Evalues %d %e  %e \n", i, creal(evalues[idx[i]]), cimag(evalues[idx[i]]));
      fflush(stdout);
    }
    
    /* Copy the first nr_ev eigenvectors to work */
    for(i = 0; i < ne; i++) {
      for(l = 0; l < m; l++) {
	work[i][l] = G[idx[i]][l];
      }
    }
    /* Orthonormalize them */
    for(i = 0; i < ne; i++) {
      work[i][m] = 0.;
      short_ModifiedGS(work[i], m, i, work[0], mp1); 
    }
    /* Orthonormalize c - \bar H d to work */
    short_ModifiedGS(c, m+1, ne, work[0], mp1);
    for(i = 0; i < mp1; i++) {
      work[nr_ev][i] = c[i];
    }
    /* Now compute \bar H = P^T_k+1 \bar H_m P_k */
    for(i = 0; i < mp1; i++) {
      for(l = 0; l < mp1; l++) {
	H[i][l] = 0.;
      }
    }    

    _FT(zgemm)("N", "N", &mp1, &ne, &_m, &cpone, work2[0], &mp1, work[0], &mp1, &czero, G[0], &mp1, 1, 1); 
    _FT(zgemm)("C", "N", &np1, &ne , &mp1, &cpone, work[0], &mp1, G[0], &mp1, &czero, H[0], &mp1, 1, 1);

    if(g_debug_level > 3) {
      for(i = 0; i < ne+1; i++) {
	for(l = 0; l < ne+1; l++) {
	  if(g_proc_id == 0) {
	    printf("(g[%d], g[%d]) = %e, %e\n", i, l, creal(short_scalar_prod(work[i], work[l], m+1)), 
		   creal(short_scalar_prod(work[i], work[l], m+1)));
	    printf("(g[%d], g[%d]) = %e, %e\n", l, i, creal(short_scalar_prod(work[l], work[i], m+1)), 
		   creal(short_scalar_prod(work[l], work[i], m+1)));
	  }
	}
      }
    }
    /* V_k+1 = V_m+1 P_k+1 */
/*     _FT(zgemm)("N", "N", &_N, &np1, &mp1, &cpone, (_Complex double*)V[0], &V2, work[0], &mp1, &czero, (_Complex double*)Z[0], &V2, 1, 1);  */
    for(l = 0; l < np1; l++) {
      mul(Z[l], work[l][0], V[0], N);
      for(i = 1; i < mp1; i++) {
	assign_add_mul(Z[l], V[i], work[l][i], N);
      }
    }
    /* copy back to V */
    for(i = 0; i < np1; i++) {
      assign(V[i], Z[i], N); 
    }
    /* Reorthogonalise v_nr_ev */
    ModifiedGS((_Complex double*)V[nr_ev], _N, nr_ev, (_Complex double*)V[0], V2);  
    if(g_debug_level > 3) {
      for(i = 0; i < np1; i++) {
	for(l = 0; l < np1; l++) {
	  tmp1 = scalar_prod(V[l], V[i], N, 1);
	  if(g_proc_id == 0) {
	    printf("(V[%d], V[%d]) = %e %e %d %d %d %d %d %d %e %e\n", l, i, creal(tmp1), cimag(tmp1), np1, mp1, ne, _m, _N, V2, creal(H[l][i]), cimag(H[l][i]));
	  }
	}
      }
    }
    /* Copy the content of H to work, work2 and G */
    for(i=0; i < mp1; i++) { 
      for(l = 0; l < mp1; l++) { 
 	G[i][l] = H[i][l];
	work2[i][l] = H[i][l];
	work[l][i] = conj(H[i][l]);
      }
    }

    for(j = ne; j < m; j++) {
      /* solver_field[0]=A*v_j */
      f(solver_field[1], V[j]);
      
      /* Set h_ij and omega_j */
      /* solver_field[1] <- omega_j */
/*       assign(solver_field[1], solver_field[0], N); */
      for(i = 0; i <= j; i++){
	H[j][i] = scalar_prod(V[i], solver_field[1], N, 1);  
	/* H, G, work and work2 are now all in Fortran storage: columns first */
	G[j][i] = H[j][i];
	work2[j][i] = H[j][i];
	work[i][j] = conj(H[j][i]);
	assign_diff_mul(solver_field[1], V[i], H[j][i], N);
      }
      beta2 = square_norm(solver_field[1], N, 1);
      H[j][j+1] = sqrt(beta2);
      G[j][j+1] = H[j][j+1];
      work2[j][j+1] = H[j][j+1];
      work[j+1][j] = conj(H[j][j+1]);
      mul_r(V[(j+1)], 1./creal(H[j][j+1]), solver_field[1], N);
    }

    /* Solve the least square problem for alpha*/
    /* This produces c=V_m+1*r0 */
    for(i = 0; i < mp1; i++) {      
      c[i] = scalar_prod(V[i], r0, N, 1);  
      alpha[i] = c[i];
      if(g_proc_id == 0 && g_debug_level > 3) {
	printf("c: %e %e err = %e\n", creal(c[i]), cimag(c[i]), err);
      }
    }
    if(lswork == NULL) {
      lwork = -1;
      _FT(zgels)("N", &mp1, &_m, &one, H[0], &mp1, alpha, &mp1, &tmp1, &lwork, &info, 1);
      lwork = (int)creal(tmp1);
      lswork = (_Complex double*)malloc(lwork*sizeof(_Complex double));
    }
    _FT(zgels)("N", &mp1, &_m, &one, H[0], &mp1, alpha, &mp1, lswork, &lwork, &info, 1);
    if(g_proc_id == 0 && g_debug_level > 3) {
      printf("zgels returned info = %d\n", info);
      fflush(stdout);
    }
    /* Compute the new solution vector */
    for(i = 0; i < m; i++){
      if(g_proc_id == 0 && g_debug_level > 3) {
	printf("alpha: %e %e\n", creal(alpha[i]), cimag(alpha[i]));
      }
      assign_add_mul(x0, V[i], alpha[i], N);
    }
  }


  /* If maximal number of restart is reached */
  assign(P, x0, N);
  finalize_solver(solver_field, nr_sf);
  return(-1);
}

_Complex double short_scalar_prod(_Complex double * const y, _Complex double * const x, const int N)
{
  _Complex double res = 0.0;

  for (int ix = 0; ix < N; ++ix)
    res += conj(y[ix]) * x[ix];
  return(res);

}

void short_ModifiedGS(_Complex double v[], int n, int m, _Complex double A[], int lda)
{
  double r;
  for (int i = 0; i < m; ++i)
  {
     _Complex double s = -short_scalar_prod(A + i * lda, v, n); 
    _FT(zaxpy)(&n, &s, A+i*lda, &one, v, &one); 
  }
  
  r = creal(sqrt(short_scalar_prod(v, v, n)));
  for(int i = 0; i < n; ++i)
    v[i] /= r;
}

static void init_gmres_dr(const int _M, const int _V){
  static int Vo = -1;
  static int M = -1;
  static int init = 0;
  int i;

  if((M != _M)||(init == 0)||(Vo != _V)){
    if(init == 1){
      free(Z);
      free(_z);
      free(H);
      free(G);
      free(V);
      free(_h);
      free(_g);
      free(_v);
      free(alpha);
      free(c);
      free(s);
      free(evalues);
      free(work);
      free(_work);
      free(work2);
      free(_work2);
    }
    Vo = _V;
    M = _M;
    H = calloc(M+1, sizeof(_Complex double *));
    Z = calloc(M+1, sizeof(spinor *));
    G = calloc(M+1, sizeof(_Complex double *));
    V = calloc(M+1, sizeof(spinor *));
    work = calloc(M+1, sizeof(_Complex double *));
    work2 = calloc(M+1, sizeof(_Complex double *));
#if (defined SSE || defined SSE2)
    _h = calloc((M+2)*(M+1), sizeof(_Complex double));
    H[0] = (_Complex double *)(((unsigned long int)(_h)+ALIGN_BASE)&~ALIGN_BASE); 
    _work = calloc((M+2)*(M+1), sizeof(_Complex double));
    work[0] = (_Complex double *)(((unsigned long int)(_work)+ALIGN_BASE)&~ALIGN_BASE); 
    _work2 = calloc((M+2)*(M+1), sizeof(_Complex double));
    work2[0] = (_Complex double *)(((unsigned long int)(_work2)+ALIGN_BASE)&~ALIGN_BASE); 
    _g = calloc((M+2)*(M+1), sizeof(_Complex double));
    G[0] = (_Complex double *)(((unsigned long int)(_g)+ALIGN_BASE)&~ALIGN_BASE); 
    _v = calloc((M+1)*Vo+1, sizeof(spinor));
    V[0] = (spinor *)(((unsigned long int)(_v)+ALIGN_BASE)&~ALIGN_BASE);
    _z = calloc((M+1)*Vo+1, sizeof(spinor));
    Z[0] = (spinor *)(((unsigned long int)(_z)+ALIGN_BASE)&~ALIGN_BASE);
#else
    _h = calloc((M+1)*(M+1), sizeof(_Complex double));
    H[0] = _h;
    _work = calloc((M+1)*(M+1), sizeof(_Complex double));
    work[0] = _work;
    _work2 = calloc((M+1)*(M+1), sizeof(_Complex double));
    work2[0] = _work2;
    _g = calloc((M+1)*(M+1), sizeof(_Complex double));
    G[0] = _g;
    _v = calloc((M+1)*Vo, sizeof(spinor));
    V[0] = _v;
    _z = calloc((M+1)*Vo, sizeof(spinor));
    Z[0] = _z;
#endif
    s = calloc(M, sizeof(double));
    c = calloc(M+1, sizeof(_Complex double));
    alpha = calloc(M+1, sizeof(_Complex double));
    evalues = calloc(M+1, sizeof(_Complex double));
    sortarray = calloc(M+1, sizeof(double));
    idx = calloc(M+1, sizeof(int));
    for(i = 1; i < M; i++){
      V[i] = V[i-1] + Vo;
      H[i] = H[i-1] + M+1;
      Z[i] = Z[i-1] + Vo;
      G[i] = G[i-1] + M+1;
      work[i] = work[i-1] + M+1;
      work2[i] = work2[i-1] + M+1;
    }
    work[M] = work[M-1] + M+1;
    work2[M] = work2[M-1] + M+1;
    H[M] = H[M-1] + M+1;
    G[M] = G[M-1] + M+1;
    V[M] = V[M-1] + Vo;
    init = 1;
  }
}
#endif

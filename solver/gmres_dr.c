/* $Id$ */

/*******************************************************************************
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
#include"complex.h"
#include"su3.h"
#include"linalg_eo.h"
#include"diagonalise_general_matrix.h"
#include"quicksort.h"
#include"linalg/lapack.h"
#include"linalg/blas.h"
#include"solver/gram-schmidt.h"
#include"solver/gmres.h"
#include"gmres_dr.h"

#ifndef HAVE_LAPACK
/* In case there is no lapack use normal gmres */
int gmres_dr(spinor * const P,spinor * const Q, 
	  const int m, const int nr_ev, const int max_restarts,
	  const double eps_sq, const int rel_prec,
	  const int N, matrix_mult f){
  return(gmres(P, Q, m, max_restarts, eps_sq, rel_prec, N, f));
}

#else

static void init_gmres_dr(const int _M, const int _V);
complex short_scalar_prod(complex * const x, complex * const y, const int N);
void short_ModifiedGS(complex v[], int n, int m, complex A[], int lda);

static complex ** work;
static complex * _work;
static complex ** work2;
static complex * _work2;
static complex ** H;
static complex ** G;
static complex * alpha;
static complex * c;
static double * s;
static spinor ** V;
static spinor * _v;
static spinor ** Z;
static spinor * _z;
static complex * _h;
static complex * _g;
static complex * alpha;
static complex * c;
static double * s;
static complex * evalues;
static double * sortarray;
static int * idx;
static int one = 1;
static double mone = -1.;
static double pone = 1.;
static complex cmone;
static complex cpone;
static complex czero;

int gmres_dr(spinor * const P,spinor * const Q, 
	  const int m, const int nr_ev, const int max_restarts,
	  const double eps_sq, const int rel_prec,
	  const int N, matrix_mult f){

  int restart=0, i, j, k, l;
  double beta, eps, norm, beta2;
  complex *lswork = NULL;
  int lwork;
  complex tmp1, tmp2;
  int info=0;
  int _m = m, mp1 = m+1, np1 = nr_ev+1, ne = nr_ev, V2 = 12*(VOLUMEPLUSRAND)/2, _N = 12*N;
/*   init_solver_field(3); */
  double err=0.;
  spinor * r0, * x0;

  cmone.re = -1.; cmone.im=0.;
  cpone.re = 1.; cpone.im=0.;
  czero.re = 0.; czero.im = 0.;
  
  r0 = g_spinor_field[DUM_SOLVER];
  x0 = g_spinor_field[DUM_SOLVER+2];
  eps=sqrt(eps_sq);  
  init_gmres_dr(m, (VOLUMEPLUSRAND));
  norm = sqrt(square_norm(Q, N));

  assign(x0, P, N);

  /* first normal GMRES cycle */
  /* r_0=Q-AP  (b=Q, x+0=P) */
  f(r0, x0);
  diff(r0, Q, r0, N);
  
  /* v_0=r_0/||r_0|| */
  alpha[0].re=sqrt(square_norm(r0, N));
  err = alpha[0].re;
  
  if(g_proc_id == g_stdio_proc && g_debug_level > 0){
    printf("%d\t%e true residue\n", restart*m, alpha[0].re*alpha[0].re); 
    fflush(stdout);
  }
  
  if(alpha[0].re==0.){
    assign(P, x0, N);
    return(restart*m);
  }
  
  mul_r(V[0], 1./alpha[0].re, r0, N);
  
  for(j = 0; j < m; j++){
    /* g_spinor_field[DUM_SOLVER]=A*v_j */

    /* Set h_ij and omega_j */
    /* g_spinor_field[DUM_SOLVER+1] <- omega_j */    
    f(g_spinor_field[DUM_SOLVER+1], V[j]);
/*     assign(g_spinor_field[DUM_SOLVER+1], g_spinor_field[DUM_SOLVER], N); */
    for(i = 0; i <= j; i++){
      H[i][j] = scalar_prod(V[i], g_spinor_field[DUM_SOLVER+1], N);
      /* G, work and work2 are in Fortran storage: columns first */
      G[j][i] = H[i][j];
      work2[j][i] = H[i][j];
      work[i][j].re = H[i][j].re;
      work[i][j].im = -H[i][j].im;
      assign_diff_mul(g_spinor_field[DUM_SOLVER+1], V[i], H[i][j], N);
    }
    
    _complex_set(H[j+1][j], sqrt(square_norm(g_spinor_field[DUM_SOLVER+1], N)), 0.);
    G[j][j+1] = H[j+1][j];
    work2[j][j+1] = H[j+1][j];
    work[j+1][j].re =  H[j+1][j].re;
    work[j+1][j].im =  -H[j+1][j].im;
    beta2 = H[j+1][j].re*H[j+1][j].re; 
    for(i = 0; i < j; i++){
      tmp1 = H[i][j];
      tmp2 = H[i+1][j];
      _mult_real(H[i][j], tmp2, s[i]);
      _add_assign_complex_conj(H[i][j], c[i], tmp1);
      _mult_real(H[i+1][j], tmp1, s[i]);
      _diff_assign_complex(H[i+1][j], c[i], tmp2);
    }
    
    /* Set beta, s, c, alpha[j],[j+1] */
    beta = sqrt(_complex_square_norm(H[j][j]) + _complex_square_norm(H[j+1][j]));
    s[j] = H[j+1][j].re / beta;
    _mult_real(c[j], H[j][j], 1./beta);
    _complex_set(H[j][j], beta, 0.);
    _mult_real(alpha[j+1], alpha[j], s[j]);
    tmp1 = alpha[j];
    _mult_assign_complex_conj(alpha[j], c[j], tmp1);
    
    /* precision reached? */
    if(g_proc_id == g_stdio_proc && g_debug_level > 0){
      printf("%d\t%e residue\n", restart*m+j, alpha[j+1].re*alpha[j+1].re); 
      fflush(stdout);
    }
    if(((alpha[j+1].re <= eps) && (rel_prec == 0)) || ((alpha[j+1].re <= eps*norm) && (rel_prec == 1))){
      _mult_real(alpha[j], alpha[j], 1./H[j][j].re);
      assign_add_mul(x0, V[j], alpha[j], N);
      for(i = j-1; i >= 0; i--){
	for(k = i+1; k <= j; k++){
	  _mult_assign_complex(tmp1, H[i][k], alpha[k]); 
	  /* alpha[i] -= tmp1 */
	  _diff_complex(alpha[i], tmp1);
	}
	_mult_real(alpha[i], alpha[i], 1./H[i][i].re);
	assign_add_mul(x0, V[i], alpha[i], N);
      }
      for(i = 0; i < m; i++){
	alpha[i].im = 0.;
      }
      assign(P, x0, N);
      return(restart*m+j);
    }
    /* if not */
    else {
      mul_r(V[(j+1)], 1./H[j+1][j].re, g_spinor_field[DUM_SOLVER+1], N); 
    }
    
  }
  j=m-1;
  /* prepare for restart */
  _mult_real(alpha[j], alpha[j], 1./H[j][j].re);
  assign_add_mul(x0, V[j], alpha[j], N);
  if(g_proc_id == 0 && g_debug_level > 3) {
    printf("alpha: %e %e\n", alpha[j].re, alpha[j].im);
  }
  for(i = j-1; i >= 0; i--){
    for(k = i+1; k <= j; k++){
      _mult_assign_complex(tmp1, H[i][k], alpha[k]);
      _diff_complex(alpha[i], tmp1);
    }
    _mult_real(alpha[i], alpha[i], 1./H[i][i].re);
    if(g_proc_id == 0 && g_debug_level > 3) {
      printf("alpha: %e %e\n", alpha[i].re, alpha[i].im);
    }
    assign_add_mul(x0, V[i], alpha[i], N);
  }

  /* This produces c=V_m+1*r0 */
  for(i = 0; i < mp1; i++) {
    c[i] = scalar_prod(V[i], r0, N); 
    if(g_proc_id == 0 && g_debug_level > 3) {
      printf("c: %e %e err = %e\n", c[i].re, c[i].im, err);
    }
  }

  for(restart = 1; restart < max_restarts; restart++) {  

    /* compute c-\bar H \alpha */
    _FT(zgemv) ("N", &mp1, &_m, &cmone, G[0], &mp1, alpha, &one, &cpone, c, &one, 1);
    err = sqrt(short_scalar_prod(c, c, mp1).re);
    if(g_proc_id == 0 && g_debug_level > 0) {
      printf("%d\t %e short residue\n", m*restart, err*err);
    } 
    
    /* Compute new residual r0 */
    /* r_0=Q-AP  (b=Q, x+0=P) */
    if(g_debug_level > 0) {
      f(r0, x0);
      diff(r0, Q, r0, N);
      tmp1.im=sqrt(square_norm(r0, N));
      if(g_proc_id == g_stdio_proc){
	printf("%d\t%e true residue\n", m*restart, tmp1.im*tmp1.im); 
	fflush(stdout);
      }
    }
    mul(r0, c[0], V[0], N);
    for(i = 1; i < mp1; i++) {
      assign_add_mul(r0, V[i], c[i], N);
    } 
    if(g_debug_level > 3) {
      tmp1.im=sqrt(square_norm(r0, N));
      if(g_proc_id == g_stdio_proc){
	printf("%d\t%e residue\n", m*restart, tmp1.im*tmp1.im); 
	fflush(stdout);
      }
    }
    /* Stop if satisfied */
    if(err < eps){
      assign(P, x0, N);
      return(restart*m);
    }

    /* Prepare to compute harmonic Ritz pairs */
    for(i = 0; i < m-1; i++){
      alpha[i].re = 0.;
      alpha[i].im = 0.;
    }
    alpha[m-1].re = 1.;
    alpha[m-1].im = 0.;
    _FT(zgesv) (&_m, &one, work[0], &mp1, idx, alpha, &_m, &info); 
    for(i = 0; i < m; i++) {
      G[m-1][i].re += (beta2*alpha[idx[i]-1].re);
      G[m-1][i].im += (beta2*alpha[idx[i]-1].im);
    }
    if(g_proc_id == 0 && g_debug_level > 3){
      printf("zgesv returned info = %d, c[m-1]= %e, %e , idx[m-1]=%d\n", 
	     info, alpha[idx[m-1]-1].re, alpha[idx[m-1]-1].im, idx[m-1]);
    }
    /* c - \bar H * d -> c */
    /* G contains H + \beta^2 H^-He_n e_n^H */

    /* Compute harmonic Ritz pairs */
    diagonalise_general_matrix(m, G[0], mp1, alpha, evalues);
    for(i = 0; i < m; i++) {
      sortarray[i] = _complex_square_norm(evalues[i]);
      idx[i] = i;
    }
    quicksort(m, sortarray, idx);
    if(g_proc_id == g_stdio_proc && g_debug_level > 1) {
      for(i = 0; i < m; i++) {
	printf("# Evalues %d %e  %e \n", i, evalues[idx[i]].re, evalues[idx[i]].im);
      }
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
      work[i][m].re = 0.;
      work[i][m].im = 0.;
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
	H[i][l].re = 0.;
	H[i][l].im = 0.;
      }
    }    

    _FT(zgemm) ("N", "N", &mp1, &ne, &_m, &cpone, work2[0], &mp1, work[0], &mp1, &czero, G[0], &mp1, 1, 1); 
    _FT(zgemm) ("C", "N", &np1, &ne , &mp1, &cpone, work[0], &mp1, G[0], &mp1, &czero, H[0], &mp1, 1, 1);

    if(g_debug_level > 3) {
      for(i = 0; i < ne+1; i++) {
	for(l = 0; l < ne+1; l++) {
	  if(g_proc_id == 0) {
	    printf("(g[%d], g[%d]) = %e, %e\n", i, l, short_scalar_prod(work[i], work[l], m+1).re, 
		   short_scalar_prod(work[i], work[l], m+1).im);
	    printf("(g[%d], g[%d]) = %e, %e\n", l, i, short_scalar_prod(work[l], work[i], m+1).re, 
		   short_scalar_prod(work[l], work[i], m+1).im);
	  }
	}
      }
    }
    /* V_k+1 = V_m+1 P_k+1 */
/*     _FT(zgemm) ("N", "N", &_N, &np1, &mp1, &cpone, (complex*)V[0], &V2, work[0], &mp1, &czero, (complex*)Z[0], &V2, 1, 1);  */
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
    pModifiedGS((complex*)V[nr_ev], _N, nr_ev, (complex*)V[0], V2);  
    if(g_debug_level > 3) {
      for(i = 0; i < np1; i++) {
	for(l = 0; l < np1; l++) {
	  tmp1 = scalar_prod(V[l], V[i], N);
	  if(g_proc_id == 0) {
	    printf("(V[%d], V[%d]) = %e %e %d %d %d %d %d %d %e %e\n", l, i, tmp1.re, tmp1.im, np1, mp1, ne, _m, _N, V2, H[l][i].re, H[l][i].im);
	  }
	}
      }
    }
    /* Copy the content of H to work, work2 and G */
    for(i=0; i < mp1; i++) { 
      for(l = 0; l < mp1; l++) { 
 	G[i][l] = H[i][l];
	work2[i][l] = H[i][l];
	work[l][i].re = H[i][l].re;
	work[l][i].im = -H[i][l].im;
      }
    }

    for(j = ne; j < m; j++) {
      /* g_spinor_field[DUM_SOLVER]=A*v_j */
      f(g_spinor_field[DUM_SOLVER+1], V[j]);
      
      /* Set h_ij and omega_j */
      /* g_spinor_field[DUM_SOLVER+1] <- omega_j */
/*       assign(g_spinor_field[DUM_SOLVER+1], g_spinor_field[DUM_SOLVER], N); */
      for(i = 0; i <= j; i++){
	H[j][i] = scalar_prod(V[i], g_spinor_field[DUM_SOLVER+1], N);  
	/* H, G, work and work2 are now all in Fortran storage: columns first */
	G[j][i] = H[j][i];
	work2[j][i] = H[j][i];
	work[i][j].re = H[j][i].re;
	work[i][j].im = -H[j][i].im;
	assign_diff_mul(g_spinor_field[DUM_SOLVER+1], V[i], H[j][i], N);
      }
      beta2 = square_norm(g_spinor_field[DUM_SOLVER+1], N);
      _complex_set(H[j][j+1], sqrt(beta2), 0.);
      G[j][j+1] = H[j][j+1];
      work2[j][j+1] = H[j][j+1];
      work[j+1][j].re =  H[j][j+1].re;
      work[j+1][j].im =  -H[j][j+1].im;
      mul_r(V[(j+1)], 1./H[j][j+1].re, g_spinor_field[DUM_SOLVER+1], N);
    }

    /* Solve the least square problem for alpha*/
    /* This produces c=V_m+1*r0 */
    for(i = 0; i < mp1; i++) {      
      c[i] = scalar_prod(V[i], r0, N);  
      alpha[i] = c[i];
      if(g_proc_id == 0 && g_debug_level > 3) {
	printf("c: %e %e err = %e\n", c[i].re, c[i].im, err);
      }
    }
    if(lswork == NULL) {
      lwork = -1;
      _FT(zgels) ("N", &mp1, &_m, &one, H[0], &mp1, alpha, &mp1, &tmp1, &lwork, &info, 1);
      lwork = (int)tmp1.re;
      lswork = (complex*)malloc(lwork*sizeof(complex));
    }
    _FT(zgels) ("N", &mp1, &_m, &one, H[0], &mp1, alpha, &mp1, lswork, &lwork, &info, 1);
    if(g_proc_id == 0 && g_debug_level > 3) {
      printf("zgels returned info = %d\n", info);
      fflush(stdout);
    }
    /* Compute the new solution vector */
    for(i = 0; i < m; i++){
      if(g_proc_id == 0 && g_debug_level > 3) {
	printf("alpha: %e %e\n", alpha[i].re, alpha[i].im);
      }
      assign_add_mul(x0, V[i], alpha[i], N);
    }
  }


  /* If maximal number of restart is reached */
  assign(P, x0, N);

  return(-1);
}

complex short_scalar_prod(complex * const y, complex * const x, const int N) {
  complex res;
  int ix;

  res.re = 0.;
  res.im = 0.;
  for (ix = 0; ix < N; ix++){
    res.re += +x[ix].re*y[ix].re + x[ix].im*y[ix].im;
    res.im += -x[ix].re*y[ix].im + x[ix].im*y[ix].re;
  }
  return(res);

}

void short_ModifiedGS(complex v[], int n, int m, complex A[], int lda) {

  int i;
  complex s;

  for (i = 0; i < m; i++) {
    s = short_scalar_prod(A+i*lda, v, n); 
    s.re = -s.re; s.im = -s.im;
    _FT(zaxpy)(&n, &s, A+i*lda, &one, v, &one); 
  }
  s.re = sqrt(short_scalar_prod(v, v, n).re);
  for(i = 0; i < n; i++) {
    v[i].re /= s.re;
    v[i].im /= s.re;
  }
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
    H = calloc(M+1, sizeof(complex *));
    Z = calloc(M+1, sizeof(spinor *));
    G = calloc(M+1, sizeof(complex *));
    V = calloc(M+1, sizeof(spinor *));
    work = calloc(M+1, sizeof(complex *));
    work2 = calloc(M+1, sizeof(complex *));
#if (defined SSE || defined SSE2)
    _h = calloc((M+2)*(M+1), sizeof(complex));
    H[0] = (complex *)(((unsigned int)(_h)+ALIGN_BASE)&~ALIGN_BASE); 
    _work = calloc((M+2)*(M+1), sizeof(complex));
    work[0] = (complex *)(((unsigned int)(_work)+ALIGN_BASE)&~ALIGN_BASE); 
    _work2 = calloc((M+2)*(M+1), sizeof(complex));
    work2[0] = (complex *)(((unsigned int)(_work2)+ALIGN_BASE)&~ALIGN_BASE); 
    _g = calloc((M+2)*(M+1), sizeof(complex));
    G[0] = (complex *)(((unsigned int)(_g)+ALIGN_BASE)&~ALIGN_BASE); 
    _v = calloc((M+1)*Vo+1, sizeof(spinor));
    V[0] = (spinor *)(((unsigned int)(_v)+ALIGN_BASE)&~ALIGN_BASE);
    _z = calloc((M+1)*Vo+1, sizeof(spinor));
    Z[0] = (spinor *)(((unsigned int)(_z)+ALIGN_BASE)&~ALIGN_BASE);
#else
    _h = calloc((M+1)*(M+1), sizeof(complex));
    H[0] = _h;
    _work = calloc((M+1)*(M+1), sizeof(complex));
    work[0] = _work;
    _work2 = calloc((M+1)*(M+1), sizeof(complex));
    work2[0] = _work2;
    _g = calloc((M+1)*(M+1), sizeof(complex));
    G[0] = _g;
    _v = calloc((M+1)*Vo, sizeof(spinor));
    V[0] = _v;
    _z = calloc((M+1)*Vo, sizeof(spinor));
    Z[0] = _z;
#endif
    s = calloc(M, sizeof(double));
    c = calloc(M+1, sizeof(complex));
    alpha = calloc(M+1, sizeof(complex));
    evalues = calloc(M+1, sizeof(complex));
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

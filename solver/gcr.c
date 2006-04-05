/* $Id$ */

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include"global.h"
#include"su3.h"
#include"linalg_eo.h"
#include"solver/gmres_precon.h"
/* #include"solver/mr_precon.h" */
#include"tm_operators.h"
#include"gcr.h"

static void init_gcr(const int _M, const int _V);

static complex ** a; 
static complex * _a;
static double * b;
static complex * c;
static spinor ** chi;
static spinor * _chi;
static spinor ** xi;
static spinor * _xi;
static complex * alpha;


int gcr(spinor * const P, spinor * const Q, 
	const int m, const int max_restarts,
	const double eps_sq, const int rel_prec,
	const int N, matrix_mult f) {

  int k, l, restart, i;
  double norm_sq, err;
  spinor * rho, * tmp;
  complex ctmp;

  rho = g_spinor_field[DUM_SOLVER];
  tmp = g_spinor_field[DUM_SOLVER+1];

  init_gcr(m, (VOLUME + RAND)/2);

  norm_sq = square_norm(Q, N);
  
  for(restart = 0; restart < max_restarts; restart++) {
    f(tmp, P);
    diff(rho, Q, tmp, N);
    err = square_norm(rho, N);
    if(g_proc_id == g_stdio_proc && g_debug_level > 0){
      printf("%d\t%g true GCR residue\n", restart*m, err); 
      fflush(stdout);
    }
    if(((err <= eps_sq) && (rel_prec == 0)) || ((err <= eps_sq*norm_sq) && (rel_prec == 1))) {
      return(restart*m);
    }
    for(k = 0; k < m; k++) {

      assign(xi[k], rho, N);  
/*       gmres_precon(xi[k], rho, 5, 2,   */
/* 		1.e-13, rel_prec, N, &Mtm_plus_psi_nocom);  */
      f(tmp, xi[k]); 
      /* tmp will become chi[k] */
      for(l = 0; l < k; l++) {
	a[l][k] = scalar_prod(chi[l], tmp, N);
	assign_diff_mul(tmp, chi[l], a[l][k], N);
      }
      b[k] = sqrt(square_norm(tmp, N));
      mul_r(chi[k], 1./b[k], tmp, N);
      c[k] = scalar_prod(chi[k], rho, N);
      assign_diff_mul(rho, chi[k], c[k], N);
      err = square_norm(rho, N);
      if(g_proc_id == g_stdio_proc && g_debug_level > 0){
	printf("%d\t%g GCR residue\n", restart*m+k, err); 
	fflush(stdout);
      }
      /* Precision reached? */
      if(((err <= eps_sq) && (rel_prec == 0)) || ((err <= eps_sq*norm_sq) && (rel_prec == 1))) {
	_mult_real(c[k], c[k], 1./b[k]);
	assign_add_mul(P, xi[k], c[k], N);
	for(l = k-1; l >= 0; l--) {
	  for(i = l+1; i <= k; i++) {
	    _mult_assign_complex(ctmp, a[l][i], c[i]);
	    /* c[l] -= ctmp */
	    _diff_complex(c[l], ctmp);
	  }
	  _mult_real(c[l], c[l], 1./b[l]);
	  assign_add_mul(P, xi[l], c[l], N);
	}
	return(restart*m+k);
      }
    }
    /* prepare for restart */
    k--;
    _mult_real(c[k], c[k], 1./b[k]);
    assign_add_mul(P, xi[k], c[k], N);
    for(l = k-1; l >= 0; l--) {
      for(i = l+1; i <= k; i++) {
	_mult_assign_complex(ctmp, a[l][i], c[i]);
	/* c[l] -= ctmp */
	_diff_complex(c[l], ctmp);
      }
      _mult_real(c[l], c[l], 1./b[l]);
      assign_add_mul(P, xi[l], c[l], N);
    }
  }
  return(-1);
}

static void init_gcr(const int _M, const int _V){
  static int Vo = -1;
  static int M = -1;
  static int init = 0;
  int i;
  if((M != _M)||(init == 0)||(Vo != _V)){
    if(init == 1){
      free(a);
      free(chi);
      free(_a);
      free(_chi);
      free(alpha);
      free(c);
      free(_xi);
      free(xi);
    }
    Vo = _V;
    M = _M;
    a = calloc(M+1, sizeof(complex *));
    chi = calloc(M, sizeof(spinor *));
    xi = calloc(M, sizeof(spinor *));
#if (defined SSE || defined SSE2)
    _a = calloc((M+2)*M, sizeof(complex));
    a[0] = (complex *)(((unsigned int)(_a)+ALIGN_BASE)&~ALIGN_BASE); 
    _chi = calloc(M*Vo+1, sizeof(spinor));
    chi[0] = (spinor *)(((unsigned int)(_chi)+ALIGN_BASE)&~ALIGN_BASE);
    _xi = calloc(M*Vo+1, sizeof(spinor));
    xi[0] = (spinor *)(((unsigned int)(_xi)+ALIGN_BASE)&~ALIGN_BASE);
#else
    _a = calloc((M+1)*M, sizeof(complex));
    a[0] = _a;
    _chi = calloc(M*Vo, sizeof(spinor));
    chi[0] = _chi;
    _xi = calloc(M*Vo, sizeof(spinor));
    xi[0] = _xi;
#endif
    b = calloc(M, sizeof(double));
    c = calloc(M, sizeof(complex));
    alpha = calloc(M+1, sizeof(complex));
    for(i = 1; i < M; i++){
      chi[i] = chi[i-1] + Vo;
      xi[i] = xi[i-1] + Vo;
      a[i] = a[i-1] + M;
    }
    a[M] = a[M-1] + M;
    init = 1;
  }
}

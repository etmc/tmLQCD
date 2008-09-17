/* $Id$ */

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

void ldiff(complex * Q, complex * const R, complex * const S, const int N);
double lsquare_norm(complex * const Q, const int N, const int parallel);
complex lscalar_prod(complex * const R, complex * const S, const int N, const int parallel);
void lmul_r(complex * const R, const double c, complex * const S, const int N);
void lassign_diff_mul(complex * const R, complex * const S, const complex c, const int N);
void lassign_add_mul(complex * const R, complex * const S, const complex c, const int N);

static void init_lgcr(const int _M, const int _V);
static void free_lgcr();
static complex ** a = NULL; 
static complex * _a = NULL;
static double * b = NULL;
static complex * c = NULL;
static complex ** chi = NULL;
static complex * _chi = NULL;
static complex ** xi = NULL;
static complex * _xi = NULL;
static complex * alpha = NULL;
static complex * tmp = NULL;
static complex * rho = NULL;
static int lgcr_init = 0;

int gcr4complex(complex * const P, complex * const Q, 
		const int m, const int max_restarts,
		const double eps_sq, const int rel_prec,
		const int N, const int parallel, 
		const int lda, c_matrix_mult f) {
  
  int k, l, restart, i;
  double norm_sq, err;
  complex ctmp;

  init_lgcr(m, lda);

  norm_sq = lsquare_norm(Q, N, parallel);

  for(restart = 0; restart < max_restarts; restart++) {
    f(tmp, P);
    ldiff(rho, Q, tmp, N);
    err = lsquare_norm(rho, N, parallel);
    if(g_proc_id == g_stdio_proc && g_debug_level > 1){
      printf("lGCR: %d\t%g true residue\n", restart * m, err); 
      fflush(stdout);
    }
    if(((err <= eps_sq) && (rel_prec == 0)) || ((err <= eps_sq * norm_sq) && (rel_prec == 1))) {
      return (restart * m);
    }
    for(k = 0; ; k++) {
      memcpy(xi[k], rho, N*sizeof(complex));
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
      /* Precision reached? */
      if((k == m-1) || ((err <= eps_sq) && (rel_prec == 0)) || ((err <= eps_sq*norm_sq) && (rel_prec == 1))) {
	break;
      }
    }
    /* prepare for restart */
    _mult_real(c[k], c[k], 1./b[k]);
    lassign_add_mul(P, xi[k], c[k], N);
    for(l = k-1; l >= 0; l--) {
      for(i = l+1; i <= k; i++) {
        _mult_assign_complex(ctmp, a[l][i], c[i]);
        /* c[l] -= ctmp */
        _diff_complex(c[l], ctmp);
      }
      _mult_real(c[l], c[l], 1./b[l]);
      lassign_add_mul(P, xi[l], c[l], N);
    }
  }
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
    a = calloc(M+1, sizeof(complex *));
    chi = calloc(M, sizeof(complex *));
    xi = calloc(M, sizeof(complex *));
    tmp = calloc(Vo, sizeof(complex));
    rho = calloc(Vo, sizeof(complex));
    _a = calloc((M+1)*M, sizeof(complex));
    a[0] = _a;
    _chi = calloc(M*Vo, sizeof(complex));
    chi[0] = _chi;
    _xi = calloc(M*Vo, sizeof(complex));
    xi[0] = _xi;

    b = calloc(M, sizeof(double));
    c = calloc(M, sizeof(complex));
    alpha = calloc(M+1, sizeof(complex));
    for(i = 1; i < M; i++){
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


void ldiff(complex * const Q, complex * const R, complex * const S, 
	   const int N) 
{
  int i;
  for(i = 0; i < N; i++) {
    Q[i].re = R[i].re - S[i].re;
    Q[i].im = R[i].im - S[i].im;
  }
  return;
}

double lsquare_norm(complex * const Q, const int N,
		    const int parallel) 
{
  int i;
  double nrm=0., nrm2=0.;
  
  for(i = 0; i < N; i++) {
    nrm += _complex_square_norm(Q[i]);
  }
#ifdef MPI
  if(parallel) {
    nrm2 = nrm;
    MPI_Allreduce(&nrm2, &nrm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  }
#endif
  return(nrm);
}

complex lscalar_prod(complex * const R, complex * const S, 
		     const int N, const int parallel) 
{
  complex res, res2;
  int i;
  res.re = 0.;
  res.im = 0.;
  res2 = res;
  for(i = 0; i < N; i++) {
    _add_assign_complex_conj(res, R[i], S[i]);
  }
#ifdef MPI
  if(parallel) {
    res2 = res;
    MPI_Allreduce(&res2, &res, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
  }
#endif
  return(res);
}

void lmul_r(complex * const R, const double c, complex * const S, const int N) 
{
  int i;
  for(i = 0; i < N; i++) {
    _mult_real(R[i], S[i], c);
  }
  return;
}


void lassign_add_mul(complex * const R, complex * const S, const complex c, const int N)
{
  int i;
  for(i = 0; i < N; i++) {
    _add_assign_complex(R[i], c, S[i]);
  }
  return;
}

void lassign_diff_mul(complex * const R, complex * const S, const complex c, const int N) 
{
  int i;
  complex d;
  d.re = -c.re;
  d.im = -c.im;
  for(i = 0; i < N; i++) {
    _add_assign_complex(R[i], d, S[i]);
  }
  return;
}


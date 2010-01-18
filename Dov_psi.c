/***********************************************************************
 * $Id$
 *
 * Copyright (C) 2003 Ines Wetzorke
 *               2006 Urs Wenger
 *               2004, 2009 Carsten Urbach
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
 * Action of the overlap Dirac operator D on a given spinor field
 *
 * This software is distributed under the terms of the GNU General Public
 * License (GPL)
 *
 * The externally accessible function is
 *
 *   void Dov_psi(spinor * const P, spinor * const S)
 *     Action of the overlap operator Dov on a given spinor field
 *     Dov = (1+s-m0/2){1+gamma5 Q/sqrt(Q^2)} + m0
 *     with Q = gamma5*(-(1+s)+D_W)
 *
 *   void Qov_psi(spinor * const P, spinor * const S) 
 *     Action of the hermitian overlap operator Dov on a given spinor field
 *     i.e. Qov = gamma_5 * Dov
 * 
 *************************************************************************/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "global.h"
#include "su3.h"
#include "linalg_eo.h"
#include "start.h"
#include "D_psi.h"
#include "gamma.h"
#include "chebyshev_polynomial_nd.h"
#include "solver/eigenvalues.h"
#include "solver/sub_low_ev.h"
#include "Dov_psi.h"

void addproj_q_invsqrt(spinor * const Q, spinor * const P, const int n, const int N);
/* |R>=rnorm^2 Q^2 |S> */
void norm_Q_sqr_psi(spinor * const R, spinor * const S, 
		    const double rnorm);
/* void norm_Q_n_psi(spinor *R, spinor *S, double m, int n, double rnorm)  */
/*  norm_Q_n_psi makes multiplication of any power of                      */
/*  Q== gamma5*D_W, initial vector S, final R, finally the                 */
/* vector R is multiplied by a factor rnorm^n                              */
/* |R>=rnorm^n Q^n |S>  where m is a mass                                  */
void norm_Q_n_psi(spinor * const R, spinor * const S, 
		  const int n, const double rnorm);
/* this is Q/sqrt(Q^2) */
void Q_over_sqrt_Q_sqr(spinor * const R, double * const c, 
		       const int n, spinor * const S,
		       const double rnorm, const double minev);

double ov_s = 0.6;
double m_ov = 0.;
extern const int ov_n_cheby;
extern double * ov_cheby_coef;

void Dov_psi(spinor * const P, spinor * const S) {

  int i;
  double c0,c1;
  spinor **s, *s_;
  static int n_cheby = 0;
  static int rec_coefs = 1;

  ov_s = 0.5*(1./g_kappa - 8.) - 1.;

  if(n_cheby != ov_n_cheby || rec_coefs) {
    if(ov_cheby_coef != NULL) free(ov_cheby_coef);
    ov_cheby_coef = (double*)malloc(ov_n_cheby*sizeof(double));
    chebyshev_coefs(ev_minev, 1., ov_cheby_coef, ov_n_cheby, -0.5);
    rec_coefs = 0;
    n_cheby = ov_n_cheby;
  }

  s_ = calloc(2*VOLUMEPLUSRAND+1, sizeof(spinor));
  s  = calloc(2, sizeof(spinor*));

  for(i = 0; i < 2; i++) {
#if (defined SSE3 || defined SSE2 || defined SSE)
    s[i] = (spinor*)(((unsigned long int)(s_)+ALIGN_BASE)&~ALIGN_BASE)+i*VOLUMEPLUSRAND;
#else
    s[i] = s_+i*VOLUMEPLUSRAND;
#endif
  }

  /* here we do with M = 1 + s */
  /* M + m_ov/2 + (M - m_ov/2) \gamma_5 sign(Q(-M)) */
  c0 = 1.0 + ov_s - 0.5*m_ov;
  c1 = 1.0 + ov_s + 0.5*m_ov;

  Q_over_sqrt_Q_sqr(s[0], ov_cheby_coef, ov_n_cheby, S, ev_qnorm, ev_minev);
  gamma5(s[1], s[0], VOLUME);
  assign_mul_add_mul_r(s[1], S, c0, c1, VOLUME);
  assign(P, s[1], VOLUME);

  free(s);
  free(s_);
  return;
}

void Qov_psi(spinor * const P, spinor * const S) {
  Dov_psi(P, S);
  gamma5(P, P, VOLUME);
  return;
}


void addproj_q_invsqrt(spinor * const Q, spinor * const P, const int n, const int N) {
  
  int j;
  spinor *aux_ = NULL, *aux;
  complex cnorm, lambda;
  static double save_ev[2]={-1.,-1.};
  static int * ev_sign = NULL;
  
  if(eigenvls[0] != save_ev[0] && eigenvls[1] != save_ev[1] ) {
    if(g_proc_id == 0 && g_debug_level > 1) {
      printf("# Recomputing eigenvalue signs!\n");
      fflush(stdout);
    }
    for(j = 0; j < 2; j++) {
      save_ev[j] = eigenvls[j];
    }
    free(ev_sign);
    ev_sign = (int*) malloc(n * sizeof(int));
#if ( defined SSE || defined SSE2 || defined SSE3)
    aux_=calloc(VOLUMEPLUSRAND+1, sizeof(spinor));
    aux = (spinor *)(((unsigned long int)(aux_)+ALIGN_BASE)&~ALIGN_BASE);
#else
    aux_=calloc(VOLUMEPLUSRAND, sizeof(spinor));
    aux = aux_;
#endif
    for(j=0; j < n; j++) {
      D_psi(aux, &(eigenvectors[j*evlength]));
      gamma5(aux, aux, N);
      
      lambda = scalar_prod(&(eigenvectors[j*evlength]), aux, N, 1);
      if (lambda.re < 0) {
	ev_sign[j] = -1;
      }
      else {
	ev_sign[j] = 1;
      }
    }
    free(aux_);
  }

  for(j = 0; j < n-1; j++) {
    cnorm = scalar_prod(&(eigenvectors[j*evlength]), P, N, 1);
    
    cnorm.re *= (double)ev_sign[j];
    cnorm.im *= (double)ev_sign[j];
    
    assign_add_mul(Q, &(eigenvectors[j*evlength]), cnorm, N);
  }
  return;
}


/* |R>=rnorm^2 Q^2 |S> */
void norm_Q_sqr_psi(spinor * const R, spinor * const S, 
		    const double rnorm) { 

  spinor *aux_,*aux;
#if ( defined SSE || defined SSE2 || defined SSE3 )
  aux_=calloc(VOLUMEPLUSRAND+1, sizeof(spinor));
  aux = (spinor *)(((unsigned long int)(aux_)+ALIGN_BASE)&~ALIGN_BASE);
#else
  aux_=calloc(VOLUMEPLUSRAND, sizeof(spinor));
  aux = aux_;
#endif

  /* Term -1-s is done in D_psi! does this comment make sense for HMC? */
  /* no, it doesn't, we do have to work on this */
  /* here we need to set kappa = 1./(2 (-1-s) + 8) */
  D_psi(R, S);
  gamma5(aux, R, VOLUME);
  D_psi(R, aux);
  gamma5(R, R, VOLUME);
  mul_r(R, rnorm*rnorm, R, VOLUME);

  free(aux_);
  return;
}

/* void norm_Q_n_psi(spinor *R, spinor *S, double m, int n, double rnorm)  */
/*  norm_Q_n_psi makes multiplication of any power of                      */
/*  Q== gamma5*D_W, initial vector S, final R, finally the                 */
/* vector R is multiplied by a factor rnorm^n                              */
/* |R>=rnorm^n Q^n |S>                                                     */
void norm_Q_n_psi(spinor * const R, spinor * const S, 
		  const int n, const double rnorm) { 

  int i;
  double s_par, npar = 1.;
  spinor *aux_,*aux;
#if (defined SSE || defined SSE2 || defined SSE3)
  aux_=calloc(VOLUMEPLUSRAND+1, sizeof(spinor));
  aux = (spinor *)(((unsigned long int)(aux_)+ALIGN_BASE)&~ALIGN_BASE);
#else
  aux_=calloc(VOLUMEPLUSRAND, sizeof(spinor));
  aux = aux_;
#endif
  assign(aux, S, VOLUME);
  
  s_par=-(ov_s+1);
  
  for(i=0; i < n; i++){
    D_psi(R, aux);
    /* Term -1-s is done in D_psi! does this comment make sense for HMC? */
    gamma5(aux, R, VOLUME);
    npar *= rnorm;
  }
  mul_r(R, npar, aux, VOLUME);
  free(aux_);
  return;
}

void Q_over_sqrt_Q_sqr(spinor * const R, double * const c, 
		       const int n, spinor * const S,
		       const double rnorm, const double minev) {
  
  int j;
  double fact1, fact2, temp1, temp2, temp3, temp4, maxev, tnorm;
  spinor *sv_, *sv, *d_, *d, *dd_, *dd, *aux_, *aux, *aux3_, *aux3;
  double ap_eps_sq = 0.;
  
#if ( defined SSE || defined SSE2 || defined SSE3)
  sv_  = calloc(VOLUMEPLUSRAND+1, sizeof(spinor));
  sv   = (spinor *)(((unsigned long int)(sv_)+ALIGN_BASE)&~ALIGN_BASE);
  d_   = calloc(VOLUMEPLUSRAND+1, sizeof(spinor));
  d    = (spinor *)(((unsigned long int)(d_)+ALIGN_BASE)&~ALIGN_BASE);
  dd_  = calloc(VOLUMEPLUSRAND+1, sizeof(spinor));
  dd   = (spinor *)(((unsigned long int)(dd_)+ALIGN_BASE)&~ALIGN_BASE);
  aux_ = calloc(VOLUMEPLUSRAND+1, sizeof(spinor));
  aux  = (spinor *)(((unsigned long int)(aux_)+ALIGN_BASE)&~ALIGN_BASE);
  aux3_= calloc(VOLUMEPLUSRAND+1, sizeof(spinor));
  aux3 = (spinor *)(((unsigned long int)(aux3_)+ALIGN_BASE)&~ALIGN_BASE);
#else
  sv_=calloc(VOLUMEPLUSRAND, sizeof(spinor));
  sv = sv_;
  d_=calloc(VOLUMEPLUSRAND, sizeof(spinor));
  d = d_;
  dd_=calloc(VOLUMEPLUSRAND, sizeof(spinor));
  dd = dd_;
  aux_=calloc(VOLUMEPLUSRAND, sizeof(spinor));
  aux = aux_;
  aux3_=calloc(VOLUMEPLUSRAND, sizeof(spinor));
  aux3 = aux3_;
#endif

  eigenvalues_for_cg_computed = no_eigenvalues - 1;
  if(eigenvalues_for_cg_computed < 0) eigenvalues_for_cg_computed = 0;
  maxev=1.0;
  
  fact1=4/(maxev-minev);
  fact2=-2*(maxev+minev)/(maxev-minev);
  
  zero_spinor_field(d, VOLUME);
  zero_spinor_field(dd, VOLUME); 
  
  if(1) assign_sub_lowest_eigenvalues(aux3, S, no_eigenvalues-1, VOLUME);
  else assign(aux3, S, VOLUME);
  
  /* Check whether switch for adaptive precision is on */
  /* this might be implemented again in the future */
  /* Use the 'old' version using Clenshaw's recursion for the 
     Chebysheff polynomial 
  */
  if(1) {
    for (j = n-1; j >= 1; j--) {
      assign(sv, d, VOLUME); 
      
      if ( (j%10) == 0 ) {
	assign_sub_lowest_eigenvalues(aux, d, no_eigenvalues-1, VOLUME);
      }
      else {
	assign(aux, d, VOLUME);
      }
      
      norm_Q_sqr_psi(R, aux, rnorm);
      printf("%e %e\n", R[0].s0.c0.re, R[0].s0.c0.im);
      printf("%e %e\n", R[0].s1.c0.re, R[0].s1.c0.im);
      temp1=-1.0;
      temp2=c[j];
      assign_mul_add_mul_add_mul_add_mul_r(d, R, dd, aux3, fact2, fact1, temp1, temp2, VOLUME);
      assign(dd, sv, VOLUME);
    } 
    
    if(1) assign_sub_lowest_eigenvalues(R, d, no_eigenvalues-1, VOLUME);
    else assign(R, d, VOLUME);
    
    norm_Q_sqr_psi(aux, R, rnorm);
    temp1=-1.0;
    temp2=c[0]/2.;
    temp3=fact1/2.;
    temp4=fact2/2.;
    assign_mul_add_mul_add_mul_add_mul_r(aux, d, dd, aux3, temp3, temp4, temp1, temp2, VOLUME);
    norm_Q_n_psi(R, aux, 1, rnorm);
  }
  else {
    /* Use the adaptive precision version using the forward recursion 
       for the Chebysheff polynomial 
    */
    
    /* d = T_0(Q^2) */
    assign(d, aux3, VOLUME);
    /* dd = T_1(Q^2) */
    norm_Q_sqr_psi(dd, d, rnorm);
    temp3 = fact1/2.;
    temp4 = fact2/2.;  
    assign_mul_add_mul_r(dd, d, temp3, temp4, VOLUME);
    /* r = c_1 T_1(Q^2) + 1./2 c_0 */
    temp1 = c[1];
    temp2 = c[0]/2.;
    mul_add_mul_r(R, dd, d, temp1, temp2, VOLUME);
    
    temp1=-1.0;
    for (j = 2; j <= n-1; j++) {
      /* aux = T_j(Q^2) = 2 Q^2 T_{j-1}(Q^2) - T_{j-2}(Q^2) */
      norm_Q_sqr_psi(aux, dd, rnorm);
      assign_mul_add_mul_add_mul_r(aux, dd, d, fact1, fact2, temp1, VOLUME);
      /* r = r + c_j T_j(Q^2) */
      temp2 = c[j];
      assign_add_mul_r(R, aux, temp2, VOLUME);
      /* The stoppping criterio tnorm = |T_j(Q^2)| */
      tnorm=square_norm(aux, VOLUME, 1);
      tnorm*=(temp2*temp2);
      
      /*
	auxnorm=square_norm(R);
	if(g_proc_id == g_stdio_proc){printf("j= %d\t|c T|^2= %g\t c_j= %g\t|r|^2= %g\n",j,tnorm,temp2,auxnorm); fflush( stdout);};
      */
      
      if(tnorm < ap_eps_sq) break; 
       /* d = T_{j-1}(Q^2) */
      assign(d, dd, VOLUME);
      /* dd = T_{j}(Q^2) */
      assign(dd, aux, VOLUME);
    }
    if(g_proc_id == g_stdio_proc && g_debug_level > 0) {
      printf("Order of Chebysheff approximation = %d\n",j); 
      fflush( stdout);
    }
     
    /* r = Q r */
    assign(aux, R, VOLUME); 
    norm_Q_n_psi(R, aux, 1, rnorm);

  }
  /* add in piece from projected subspace */
  addproj_q_invsqrt(R, S, no_eigenvalues-1, VOLUME);
  
  free(sv_);
  free(d_);
  free(dd_);
  free(aux_);
  free(aux3_);
  return;
}

void CheckApproximation(spinor * const P, spinor * const S) {

  int i;
  double c0,c1;
  spinor **s, *s_;
  static int n_cheby = 0;
  static int rec_coefs = 1;

  ov_s = 0.5*(1./g_kappa - 8.) - 1.;

  if(n_cheby != ov_n_cheby || rec_coefs) {
    if(ov_cheby_coef != NULL) free(ov_cheby_coef);
    ov_cheby_coef = (double*)malloc(ov_n_cheby*sizeof(double));
    chebyshev_coefs(ev_minev, 1., ov_cheby_coef, ov_n_cheby, -0.5);
    rec_coefs = 0;
    n_cheby = ov_n_cheby;
  }

  s_ = calloc(2*VOLUMEPLUSRAND+1, sizeof(spinor));
  s  = calloc(2, sizeof(spinor*));

  for(i = 0; i < 2; i++) {
#if (defined SSE3 || defined SSE2 || defined SSE)
    s[i] = (spinor*)(((unsigned long int)(s_)+ALIGN_BASE)&~ALIGN_BASE)+i*VOLUMEPLUSRAND;
#else
    s[i] = s_+i*VOLUMEPLUSRAND;
#endif
  }

  /* here we do with M = 1 + s */
  /* M + m_ov/2 + (M - m_ov/2) \gamma_5 sign(Q(-M)) */
//  c0 = 1.0 + ov_s - 0.5*m_ov;
//  c1 = 1.0 + ov_s + 0.5*m_ov;

  Q_over_sqrt_Q_sqr(s[0], ov_cheby_coef, ov_n_cheby, S, ev_qnorm, ev_minev);
  Q_over_sqrt_Q_sqr(s[1], ov_cheby_coef, ov_n_cheby, s[0], ev_qnorm, ev_minev);

//  gamma5(s[1], s[0], VOLUME);
  //assign_mul_add_mul_r(s[1], S, 1.0, -1., VOLUME);
  assign(P, s[1], VOLUME);

  free(s);
  free(s_);
  return;
}


static char const rcsid[] = "$Id$";

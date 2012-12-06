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
 * Generalized minimal residual (GMRES) with a maximal number of restarts.    
 * Solves Q=AP for complex regular matrices A. This is a special version for 
 * precondition.
 *           
 * int gmres(spinor * const P,spinor * const Q, 
 *	   const int m, const int max_restarts,
 *	   const double eps_sq, matrix_mult f)
 *                                                                 
 * Returns the number of iterations needed or -1 if maximal number of restarts  
 * has been reached.                                                           
 *
 * Inout:                                                                      
 *  spinor * P       : guess for the solving spinor                                             
 * Input:                                                                      
 *  spinor * Q       : source spinor
 *  int m            : Maximal dimension of Krylov subspace                                     
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
#include"su3.h"
#include"linalg_eo.h"
#include"start.h"
#include "solver_field.h"
#include"gmres_precon.h"

#ifdef _SOLVER_OUTPUT
#define _SO(x) x
#else
#define _SO(x)
#endif

static void init_pgmres(const int _M, const int _V);


static complex ** H;
static complex * alpha;
static complex * c;
static double * s;
static spinor ** V;
static spinor * _v;
static complex * _h;
static complex * alpha;
static complex * c;
static double * s;

int gmres_precon(spinor * const P,spinor * const Q, 
	  const int m, const int max_restarts,
	  const double eps_sq, const int rel_prec,
	  const int N, matrix_mult f){

  int restart, i, j, k;
  double beta, eps, norm;
  complex tmp1, tmp2;
  spinor ** solver_field = NULL;
  const int nr_sf = 2;
  
  if(N == VOLUME) {
    init_solver_field(&solver_field, VOLUMEPLUSRAND, nr_sf);
  }
  else {
    init_solver_field(&solver_field, VOLUMEPLUSRAND/2, nr_sf);
  }
  eps=sqrt(eps_sq);
  init_pgmres(m, VOLUMEPLUSRAND);

  norm = sqrt(square_norm(Q, N, 1));

/*   assign(solver_field[1], P, N); */
  zero_spinor_field(solver_field[1], N);
  for(restart = 0; restart < max_restarts; restart++){
    /* r_0=Q-AP  (b=Q, x+0=P) */
    f(solver_field[1], solver_field[1]);
    diff(solver_field[1], Q, solver_field[1], N);
    /* v_0=r_0/||r_0|| */
    alpha[0].re=sqrt(square_norm(solver_field[1], N, 1));

/*     if(alpha[0].re == 0.){  */
/*        assign(P, solver_field[1], N);  */
/*        return(restart*m);  */
/*     }  */

    if(alpha[0].re != 0.) {
      mul_r(V[0], 1./alpha[0].re, solver_field[1], N);

      for(j = 0; j < m; j++){
	/* solver_field[1]=A*v_j */
	
	f(solver_field[1], V[j]);
	
	/* Set h_ij and omega_j */
	/* solver_field[0] <- omega_j */
	assign(solver_field[0], solver_field[1], N);
	for(i = 0; i <= j; i++){
	  H[i][j] = scalar_prod(V[i], solver_field[0], N, 1);
	  assign_diff_mul(solver_field[0], V[i], H[i][j], N);
	}
	
	_complex_set(H[j+1][j], sqrt(square_norm(solver_field[0], N, 1)), 0.);
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
	  printf("gmres precon\t%d\t%g residue\n", restart*m+j, alpha[j+1].re*alpha[j+1].re); 
	  fflush(stdout);
	}
	if(((alpha[j+1].re <= eps) && (rel_prec == 0)) || ((alpha[j+1].re <= eps*norm) && (rel_prec == 1))){
	  _mult_real(alpha[j], alpha[j], 1./H[j][j].re);
	  assign_add_mul(solver_field[1], V[j], alpha[j], N);
	  for(i = j-1; i >= 0; i--){
	    for(k = i+1; k <= j; k++){
	      _mult_assign_complex(tmp1, H[i][k], alpha[k]); 
	      _diff_complex(alpha[i], tmp1);
	    }
	    _mult_real(alpha[i], alpha[i], 1./H[i][i].re);
	    assign_add_mul(solver_field[1], V[i], alpha[i], N);
	  }
	  for(i = 0; i < m; i++){
	    alpha[i].im = 0.;
	  }
	  assign(P, solver_field[1], N);
	  finalize_solver(solver_field, nr_sf);
	  return(restart*m+j);
	}
	/* if not */
	else{
	  if(j != m-1){
	    mul_r(V[(j+1)], 1./H[j+1][j].re, solver_field[0], N);
	  }
	}
      }
    
      j=m-1;
      /* prepare for restart */
      _mult_real(alpha[j], alpha[j], 1./H[j][j].re);
      assign_add_mul(solver_field[1], V[j], alpha[j], N);
      for(i = j-1; i >= 0; i--){
	for(k = i+1; k <= j; k++){
	  _mult_assign_complex(tmp1, H[i][k], alpha[k]);
	  _diff_complex(alpha[i], tmp1);
	}
	_mult_real(alpha[i], alpha[i], 1./H[i][i].re);
	assign_add_mul(solver_field[1], V[i], alpha[i], N);
      }
      for(i = 0; i < m; i++){
	alpha[i].im = 0.;
      }
    }
  }

  /* If maximal number of restarts is reached */
  assign(P, solver_field[1], N);
  
  finalize_solver(solver_field, nr_sf);
  return(-1);
}

static void init_pgmres(const int _M, const int _V){
  static int Vo = -1;
  static int M = -1;
  static int pinit = 0;
  int i;
  if((M != _M)||(pinit == 0)||(Vo != _V)){
    if(pinit == 1){
      free(H);
      free(V);
      free(_h);
      free(_v);
      free(alpha);
      free(c);
      free(s);
    }
    Vo = _V;
    M = _M;
    H = calloc(M+1, sizeof(complex *));
    V = calloc(M, sizeof(spinor *));
#if (defined SSE || defined SSE2)
    _h = calloc((M+2)*M, sizeof(complex));
    H[0] = (complex *)(((unsigned long int)(_h)+ALIGN_BASE)&~ALIGN_BASE); 
    _v = calloc(M*Vo+1, sizeof(spinor));
    V[0] = (spinor *)(((unsigned long int)(_v)+ALIGN_BASE)&~ALIGN_BASE);
#else
    _h = calloc((M+1)*M, sizeof(complex));
    H[0] = _h;
    _v = calloc(M*Vo, sizeof(spinor));
    V[0] = _v;
#endif
    s = calloc(M, sizeof(double));
    c = calloc(M, sizeof(complex));
    alpha = calloc(M+1, sizeof(complex));
    for(i = 1; i < M; i++){
      V[i] = V[i-1] + Vo;
      H[i] = H[i-1] + M;
    }
    H[M] = H[M-1] + M;
    pinit = 1;
  }
}



complex scalar_prod_nocom(spinor * const S,spinor * const R, const int N){
  int ix;
  static double ks,kc,ds,tr,ts,tt;
  spinor *s,*r;
  complex c;
  
  /* Real Part */

  ks=0.0;
  kc=0.0;
  
  for (ix = 0; ix < N; ix++){
    s=(spinor *) S + ix;
    r=(spinor *) R + ix;
    
    ds=(*r).s0.c0.re*(*s).s0.c0.re+(*r).s0.c0.im*(*s).s0.c0.im+
       (*r).s0.c1.re*(*s).s0.c1.re+(*r).s0.c1.im*(*s).s0.c1.im+
       (*r).s0.c2.re*(*s).s0.c2.re+(*r).s0.c2.im*(*s).s0.c2.im+
       (*r).s1.c0.re*(*s).s1.c0.re+(*r).s1.c0.im*(*s).s1.c0.im+
       (*r).s1.c1.re*(*s).s1.c1.re+(*r).s1.c1.im*(*s).s1.c1.im+
       (*r).s1.c2.re*(*s).s1.c2.re+(*r).s1.c2.im*(*s).s1.c2.im+
       (*r).s2.c0.re*(*s).s2.c0.re+(*r).s2.c0.im*(*s).s2.c0.im+
       (*r).s2.c1.re*(*s).s2.c1.re+(*r).s2.c1.im*(*s).s2.c1.im+
       (*r).s2.c2.re*(*s).s2.c2.re+(*r).s2.c2.im*(*s).s2.c2.im+
       (*r).s3.c0.re*(*s).s3.c0.re+(*r).s3.c0.im*(*s).s3.c0.im+
       (*r).s3.c1.re*(*s).s3.c1.re+(*r).s3.c1.im*(*s).s3.c1.im+
       (*r).s3.c2.re*(*s).s3.c2.re+(*r).s3.c2.im*(*s).s3.c2.im;

    /* Kahan Summation */    
    tr=ds+kc;
    ts=tr+ks;
    tt=ts-ks;
    ks=ts;
    kc=tr-tt;
  }
  kc=ks+kc;

  c.re = kc;

  /* Imaginary Part */

  ks=0.0;
  kc=0.0;
  
  for (ix=0;ix<N;ix++){
    s=(spinor *) S + ix;
    r=(spinor *) R + ix;
    
    ds=-(*r).s0.c0.re*(*s).s0.c0.im+(*r).s0.c0.im*(*s).s0.c0.re-
       (*r).s0.c1.re*(*s).s0.c1.im+(*r).s0.c1.im*(*s).s0.c1.re-
       (*r).s0.c2.re*(*s).s0.c2.im+(*r).s0.c2.im*(*s).s0.c2.re-
       (*r).s1.c0.re*(*s).s1.c0.im+(*r).s1.c0.im*(*s).s1.c0.re-
       (*r).s1.c1.re*(*s).s1.c1.im+(*r).s1.c1.im*(*s).s1.c1.re-
       (*r).s1.c2.re*(*s).s1.c2.im+(*r).s1.c2.im*(*s).s1.c2.re-
       (*r).s2.c0.re*(*s).s2.c0.im+(*r).s2.c0.im*(*s).s2.c0.re-
       (*r).s2.c1.re*(*s).s2.c1.im+(*r).s2.c1.im*(*s).s2.c1.re-
       (*r).s2.c2.re*(*s).s2.c2.im+(*r).s2.c2.im*(*s).s2.c2.re-
       (*r).s3.c0.re*(*s).s3.c0.im+(*r).s3.c0.im*(*s).s3.c0.re-
       (*r).s3.c1.re*(*s).s3.c1.im+(*r).s3.c1.im*(*s).s3.c1.re-
       (*r).s3.c2.re*(*s).s3.c2.im+(*r).s3.c2.im*(*s).s3.c2.re;
    
    tr=ds+kc;
    ts=tr+ks;
    tt=ts-ks;
    ks=ts;
    kc=tr-tt;
  }
  kc=ks+kc;

  c.im = kc;
  return(c);
}

double square_norm_nocom(spinor * const P, const int N) {
  int ix;
  static double ks,kc,ds,tr,ts,tt;
  spinor *s;
  
  ks = 0.0;
  kc = 0.0;
  
  /* Change due to even-odd preconditioning : VOLUME   to VOLUME/2 */   
  for (ix  =  0; ix < N; ix++) {
    s = P + ix;
    
    ds = (*s).s0.c0.re*(*s).s0.c0.re + (*s).s0.c0.im*(*s).s0.c0.im + 
         (*s).s0.c1.re*(*s).s0.c1.re + (*s).s0.c1.im*(*s).s0.c1.im + 
         (*s).s0.c2.re*(*s).s0.c2.re + (*s).s0.c2.im*(*s).s0.c2.im + 
         (*s).s1.c0.re*(*s).s1.c0.re + (*s).s1.c0.im*(*s).s1.c0.im + 
         (*s).s1.c1.re*(*s).s1.c1.re + (*s).s1.c1.im*(*s).s1.c1.im + 
         (*s).s1.c2.re*(*s).s1.c2.re + (*s).s1.c2.im*(*s).s1.c2.im + 
         (*s).s2.c0.re*(*s).s2.c0.re + (*s).s2.c0.im*(*s).s2.c0.im + 
         (*s).s2.c1.re*(*s).s2.c1.re + (*s).s2.c1.im*(*s).s2.c1.im + 
         (*s).s2.c2.re*(*s).s2.c2.re + (*s).s2.c2.im*(*s).s2.c2.im + 
         (*s).s3.c0.re*(*s).s3.c0.re + (*s).s3.c0.im*(*s).s3.c0.im + 
         (*s).s3.c1.re*(*s).s3.c1.re + (*s).s3.c1.im*(*s).s3.c1.im + 
         (*s).s3.c2.re*(*s).s3.c2.re + (*s).s3.c2.im*(*s).s3.c2.im;
    
    tr = ds + kc;
    ts = tr + ks;
    tt = ts-ks;
    ks = ts;
    kc = tr-tt;
  }
  kc = ks + kc;
  return kc;
}

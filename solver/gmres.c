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
 * Solves Q=AP for _Complex double regular matrices A.
 * For details see: Andreas Meister, Numerik linearer Gleichungssysteme        
 *   or the original citation:                                                 
 * Y. Saad, M.H.Schultz in GMRES: A generalized minimal residual algorithm    
 *                         for solving nonsymmetric linear systems.            
 * 			SIAM J. Sci. Stat. Comput., 7: 856-869, 1986           
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
#include "solver_field.h"
#include"gmres.h"


static void init_gmres(const int _M, const int _V);

static _Complex double ** H;
static _Complex double * alpha;
static _Complex double * c;
static double * s;
static spinor ** V;
static spinor * _v;
static _Complex double * _h;
static _Complex double * alpha;
static _Complex double * c;
static double * s;

int gmres(spinor * const P,spinor * const Q, 
	  const int m, const int max_restarts,
	  const double eps_sq, const int rel_prec,
	  const int N, const int parallel, matrix_mult f){

  int restart, i, j, k;
  double beta, eps, norm;
  _Complex double tmp1, tmp2;
  spinor ** solver_field = NULL;
  const int nr_sf = 3;

  if(N == VOLUME) {
    init_solver_field(&solver_field, VOLUMEPLUSRAND, nr_sf);
  }
  else {
    init_solver_field(&solver_field, VOLUMEPLUSRAND/2, nr_sf);
  }

  eps=sqrt(eps_sq);
  init_gmres(m, VOLUMEPLUSRAND);

  norm = sqrt(square_norm(Q, N, parallel));

  assign(solver_field[2], P, N);
  for(restart = 0; restart < max_restarts; restart++){
    /* r_0=Q-AP  (b=Q, x+0=P) */
    f(solver_field[0], solver_field[2]);
    diff(solver_field[0], Q, solver_field[0], N);

    /* v_0=r_0/||r_0|| */
    alpha[0] = sqrt(square_norm(solver_field[0], N, parallel));

    if(g_proc_id == g_stdio_proc && g_debug_level > 1){
      printf("%d\t%g true residue\n", restart*m, creal(alpha[0])*creal(alpha[0])); 
      fflush(stdout);
    }

    if(creal(alpha[0])==0.){
      assign(P, solver_field[2], N);
      finalize_solver(solver_field, nr_sf);
      return(restart*m);
    }

    mul_r(V[0], 1./creal(alpha[0]), solver_field[0], N);

    for(j = 0; j < m; j++){
      /* solver_field[0]=A*v_j */

      f(solver_field[0], V[j]);

      /* Set h_ij and omega_j */
      /* solver_field[1] <- omega_j */
      assign(solver_field[1], solver_field[0], N);
      for(i = 0; i <= j; i++){
	H[i][j] = scalar_prod(V[i], solver_field[1], N, parallel);
	assign_diff_mul(solver_field[1], V[i], H[i][j], N);
      }

      H[j+1][j] = sqrt(square_norm(solver_field[1], N, parallel));
      for(i = 0; i < j; i++){
	tmp1 = H[i][j];
	tmp2 = H[i+1][j];
	(H[i][j]) = (tmp2) * (s[i]);
	(H[i][j]) += conj(c[i]) * (tmp1);
	(H[i+1][j]) = (tmp1) * (s[i]);
	(H[i+1][j]) -= (c[i]) * (tmp2);
      }

      /* Set beta, s, c, alpha[j],[j+1] */
      beta = sqrt(creal(H[j][j] * conj(H[j][j])) + creal(H[j+1][j] * conj(H[j+1][j])));
      s[j] = creal(H[j+1][j]) / beta;
      (c[j]) = (H[j][j]) / beta;
      (H[j][j]) = beta;
      (alpha[j+1]) = (alpha[j]) * (s[j]);
      tmp1 = alpha[j];
      (alpha[j]) = conj(c[j]) * (tmp1);

      /* precision reached? */
      if(g_proc_id == g_stdio_proc && g_debug_level > 1){
	printf("%d\t%g residue\n", restart*m+j, creal(alpha[j+1])*creal(alpha[j+1])); 
	fflush(stdout);
      }
      if(((creal(alpha[j+1]) <= eps) && (rel_prec == 0)) || ((creal(alpha[j+1]) <= eps*norm) && (rel_prec == 1))){
	(alpha[j]) = (alpha[j]) * (1./creal(H[j][j]));
	assign_add_mul(solver_field[2], V[j], alpha[j], N);
	for(i = j-1; i >= 0; i--){
	  for(k = i+1; k <= j; k++){
 	    (tmp1) = (H[i][k]) * (alpha[k]); 
	    (alpha[i]) -= tmp1;
	  }
	  (alpha[i]) = (alpha[i]) * (1./creal(H[i][i]));
	  assign_add_mul(solver_field[2], V[i], alpha[i], N);
	}
	for(i = 0; i < m; i++){
	  alpha[i] = creal(alpha[i]);
	}
	assign(P, solver_field[2], N);
	finalize_solver(solver_field, nr_sf);
	return(restart*m+j);
      }
      /* if not */
      else{
	if(j != m-1){
	  mul_r(V[(j+1)], 1./creal(H[j+1][j]), solver_field[1], N);
	}
      }

    }
    j=m-1;
    /* prepare for restart */
    (alpha[j]) = (alpha[j]) * (1./creal(H[j][j]));
    assign_add_mul(solver_field[2], V[j], alpha[j], N);
    for(i = j-1; i >= 0; i--){
      for(k = i+1; k <= j; k++){
	(tmp1) = (H[i][k]) * (alpha[k]);
	(alpha[i]) -= tmp1;
      }
      (alpha[i]) = (alpha[i]) * (1./creal(H[i][i]));
      assign_add_mul(solver_field[2], V[i], alpha[i], N);
    }
    for(i = 0; i < m; i++){
      alpha[i] = creal(alpha[i]);
    }
  }

  /* If maximal number of restarts is reached */
  assign(P, solver_field[2], N);
  finalize_solver(solver_field, nr_sf);
  return(-1);
}

static void init_gmres(const int _M, const int _V){
  static int Vo = -1;
  static int M = -1;
  static int init = 0;
  int i;
  if((M != _M)||(init == 0)||(Vo != _V)){
    if(init == 1){
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
    H = calloc(M+1, sizeof(_Complex double *));
    V = calloc(M, sizeof(spinor *));
#if (defined SSE || defined SSE2)
    _h = calloc((M+2)*M+8, sizeof(_Complex double));
    H[0] = (_Complex double *)(((unsigned long int)(_h)+ALIGN_BASE)&~ALIGN_BASE); 
    _v = calloc(M*Vo+1, sizeof(spinor));
    V[0] = (spinor *)(((unsigned long int)(_v)+ALIGN_BASE)&~ALIGN_BASE);
#else
    _h = calloc((M+1)*M, sizeof(_Complex double));
    H[0] = _h;
    _v = calloc(M*Vo, sizeof(spinor));
    V[0] = _v;
#endif
    s = calloc(M, sizeof(double));
    c = calloc(M, sizeof(_Complex double));
    alpha = calloc(M+1, sizeof(_Complex double));
    for(i = 1; i < M; i++){
      V[i] = V[i-1] + Vo;
      H[i] = H[i-1] + M;
    }
    H[M] = H[M-1] + M;
    init = 1;
  }
  return;
}

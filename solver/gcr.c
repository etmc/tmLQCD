/***********************************************************************
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
 ***********************************************************************/

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
#include"start.h"
#include"operator/tm_operators.h"
#include"solver/poly_precon.h"
#include"solver/cg_her.h"
#include"operator/D_psi.h"
#include"Msap.h"
#include"dfl_projector.h"
#include "solver_field.h"
#include"gcr.h"

static void init_gcr(const int _M, const int _V);

static _Complex double ** a; 
static _Complex double * _a;
static double * b;
static _Complex double * c;
static spinor ** chi;
static spinor * _chi;
static spinor ** xi;
static spinor * _xi;
static _Complex double * alpha;
extern int dfl_poly_iter;

int gcr(spinor * const P, spinor * const Q, 
	const int m, const int max_restarts,
	const double eps_sq, const int rel_prec,
	const int N, const int precon, matrix_mult f) {

  int k, l, restart, i, iter = 0;
  double norm_sq, err;
  spinor * rho, * tmp;
  _Complex double ctmp;
  spinor ** solver_field = NULL;
  const int nr_sf = 2;

  if(N == VOLUME) {
    init_solver_field(&solver_field, VOLUMEPLUSRAND, nr_sf);
  }
  else {
    init_solver_field(&solver_field, VOLUMEPLUSRAND/2, nr_sf);
  }

  rho = solver_field[0];
  tmp = solver_field[1];

  init_gcr(m, N+RAND);

  norm_sq = square_norm(Q, N, 1);
  if(norm_sq < 1.e-32) {
    norm_sq = 1.;
  }
  
  for(restart = 0; restart < max_restarts; restart++) {
    dfl_sloppy_prec = 0;
    f(tmp, P);
    diff(rho, Q, tmp, N);
    err = square_norm(rho, N, 1);
    if(g_proc_id == g_stdio_proc && g_debug_level > 1){
      printf("GCR: iteration number: %d, true residue: %g\n", iter, err); 
      fflush(stdout);
    }
    if(((err <= eps_sq) && (rel_prec == 0)) || ((err <= eps_sq*norm_sq) && (rel_prec == 1))) {
      finalize_solver(solver_field, nr_sf);
      return(iter);
    }
    for(k = 0; k < m; k++) {
      
      if(precon == 0) {
	assign(xi[k], rho, N);
      }
      else {
        zero_spinor_field(xi[k], N);  
        Msap_eo(xi[k], rho, 6);   
 	/* Msap(xi[k], rho, 8); */
      }
	  
      dfl_sloppy_prec = 1;
      dfl_little_D_prec = 1.e-12;
      f(tmp, xi[k]); 
	  
      /* tmp will become chi[k] */
      for(l = 0; l < k; l++) {
        a[l][k] = scalar_prod(chi[l], tmp, N, 1);
        assign_diff_mul(tmp, chi[l], a[l][k], N);
      }
      b[k] = sqrt(square_norm(tmp, N, 1));
      mul_r(chi[k], 1./b[k], tmp, N);
      c[k] = scalar_prod(chi[k], rho, N, 1);
      assign_diff_mul(rho, chi[k], c[k], N);
      err = square_norm(rho, N, 1);
      iter ++;
      if(g_proc_id == g_stdio_proc && g_debug_level > 2){
        if(rel_prec == 1) printf("# GCR: %d\t%g >= %g iterated residue\n", iter, err, eps_sq*norm_sq); 
        else printf("# GCR: %d\t%g >= %giterated residue\n", iter, err, eps_sq);
        fflush(stdout);
      }
      /* Precision reached? */
      if((k == m-1) || ((err <= eps_sq) && (rel_prec == 0)) || ((err <= eps_sq*norm_sq) && (rel_prec == 1))) {
        break;
      }
    }

    /* prepare for restart */
    c[k] /= b[k];
    assign_add_mul(P, xi[k], c[k], N);
    for(l = k - 1; l >= 0; --l)
    {
      for(i = l+1; i <= k; ++i)
      {
        ctmp = a[l][i] * c[i];
        c[l] -= ctmp;
      }
      c[l] /= b[l];
      assign_add_mul(P, xi[l], c[l], N);
    }
  }
  finalize_solver(solver_field, nr_sf);
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
    a = calloc(M+1, sizeof(_Complex double *));
    chi = calloc(M, sizeof(spinor *));
    xi = calloc(M, sizeof(spinor *));
#if (defined SSE || defined SSE2 || defined SSE3)
    _a = calloc((M+2)*M, sizeof(_Complex double));
    a[0] = (_Complex double *)(((unsigned long int)(_a)+ALIGN_BASE)&~ALIGN_BASE); 
    _chi = calloc(M*Vo+1, sizeof(spinor));
    chi[0] = (spinor *)(((unsigned long int)(_chi)+ALIGN_BASE)&~ALIGN_BASE);
    _xi = calloc(M*Vo+1, sizeof(spinor));
    xi[0] = (spinor *)(((unsigned long int)(_xi)+ALIGN_BASE)&~ALIGN_BASE);
#else
    _a = calloc((M+1)*M, sizeof(_Complex double));
    a[0] = _a;
    _chi = calloc(M*Vo, sizeof(spinor));
    chi[0] = _chi;
    _xi = calloc(M*Vo, sizeof(spinor));
    xi[0] = _xi;
#endif
    if(_xi == NULL) {printf("Unable to allocated space for GCR iterations\n");exit(0);  }
    b = calloc(M, sizeof(double));
    c = calloc(M, sizeof(_Complex double));
    alpha = calloc(M+1, sizeof(_Complex double));
    for(i = 1; i < M; i++){
      chi[i] = chi[i-1] + Vo;
      xi[i] = xi[i-1] + Vo;
      a[i] = a[i-1] + M;
    }
    a[M] = a[M-1] + M;
    init = 1;
  }
}

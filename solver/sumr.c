/***********************************************************************
 * $Id$
 *
 * Copyright (C) 2005 Luigi Scorzato
 *               2009 Carsten Urbach
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
 * File: sumr.c
 *
 *
 * The externally accessible functions are
 *
 *
 *   int sumr(spinor * const P, spinor * const Q, int max_iter, double eps_sq)
 *     Inverter for shifted unitary matrices
 *     [C.F.Jagels L.Reichel, Num. Lin. Alg. with Appl. Vol1(6),555-570 (1994)]
 *     [first applied to the Overlap in hep-lat/0311025]
 *
 * input:
 *   Q: source
 * inout:
 *   P: initial guess and result
 *
 * Author: Luigi.Scorzato@physik.hu-berlin.de
 *
 *******************************************************************************/

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
#include "solver/matrix_mult_typedef.h"
#include "complex.h"
#include "gamma.h"
#include "solver/eigenvalues.h"
#include "solver/sub_low_ev.h"
#include "Dov_psi.h"
#include "solver_field.h"
#include "sumr.h"

#define DEBUG_SUMR 0

/* to be fixed somewhere else */


/* P output = solution , Q input = source */
int sumr(spinor * const P, spinor * const Q, const int max_iter,
	 double eps_sq){
  double sigma, delta, s, rho, zeta, z_r, tmpr, normsp, err, tau_hat;
  double ov_s, m_ov, ap_eps_sq, shift;
  complex phi, phi_hat, tau, lambda, c, gamm, alpha, eta, kappa;
  complex r_hat, r_off, r_off_old, r_diag, r_diag_old, tmpc, tmpc1, tmpc2;
  int iteration;
  /* to be fixed somewhere else */
  int switch_on_adaptive_precision=0;
  const int N=VOLUME;
  spinor *x, *r, *p, *v, *v_til, *w, *u, *b, *tmp, *tmp2;
  spinor ** solver_field = NULL;
  printf("Starting SUMR!\n");
  /*        !!!!   INITIALIZATION    !!!! */
  init_solver_field(solver_field, VOLUMEPLUSRAND, 10);
  x = solver_field[0];
  r = solver_field[1];
  p = solver_field[2];
  v = solver_field[3];
  v_til = solver_field[4];
  w = solver_field[5];
  u = solver_field[6];
  b = solver_field[7];
  tmp = solver_field[8];
  tmp2 = solver_field[9];

  assign(b, Q, N);
  assign(x, P, N);
  normsp = square_norm(P, N, 1);

  ov_s = 0.5 * (1. / g_kappa - 8.) - 1.;
  rho  = ov_s+1 - m_ov / 2.;
  zeta = ov_s+1 + m_ov / 2.;
  z_r  = zeta / rho;

  if(normsp == 0) {
    /* if a starting solution vector equal to zero is chosen */
    delta = sqrt(square_norm(b, N, 1));
    assign(r, b, N);
  }
  else {
    /* if a starting solution vector different from zero is chosen */
    Dov_psi(tmp, x);
    diff(r, b, tmp, N);
    delta = sqrt(square_norm(r, N, 1));
  }

  _complex_set(phi_hat, 1 / delta,0);
  tau_hat = delta / rho;
  zero_spinor_field(p, N);
  _complex_zero(phi);
  s = 0.;
  _complex_zero(lambda);
  _complex_zero(r_off_old);
  _complex_one(r_diag_old);
  _complex_one(gamm);
  sigma = 1.;
  _complex_one(c);
  mul(v_til, phi_hat, r, N);
  assign(v, v_til, N);

#if DEBUG_SUMR ==1
  printf("delta=%g;\t phihat=%g;\t tauhat=%g;\t w=%g;\t p=%g;\t phi=%g;\t s=%g;\t lambda=%g;\t r_off=%g;\t r_off_old=%g;\t r_diag=%g;\t r_diag_old=%g;\t gamm=%g;\t sigma=%g;\t c=%g;\t v=%g;\t v_til=%g;\t ",
	 delta,_complex_norm(phi_hat),tau_hat,square_norm(w),square_norm(p),
	 _complex_norm(phi),s,_complex_norm(lambda),_complex_norm(r_off),_complex_norm(r_off_old),_complex_norm(r_diag),_complex_norm(r_diag_old),_complex_norm(gamm),sigma,_complex_norm(c),
	 square_norm(v),square_norm(v_til));
#endif

  if(switch_on_adaptive_precision == 1) {
    ap_eps_sq = 1.0e-2 * eps_sq;
  }

  if(ov_cheby_coef==NULL) calculateOverlapPolynomial();

  /* main loop */
  for(iteration = 0; iteration < max_iter; iteration++) {
    Q_over_sqrt_Q_sqr(tmp, ov_cheby_coef, ov_n_cheby, v, ev_qnorm, ev_minev);
    gamma5(u, tmp, N);
#if DEBUG_SUMR ==1
    printf("u=%g;\t\n", square_norm(u));
#endif
    gamm = scalar_prod(v_til, u, N, 1);
    _complex_chgsig(gamm, gamm);
#if DEBUG_SUMR ==1
    printf("gamm=%g,%g;\t\n",gamm.re,gamm.im);
#endif
    sigma= sqrt((1 - _complex_norm(gamm))*(1 + _complex_norm(gamm)));
#if DEBUG_SUMR ==1
    printf("sigma=%g;\t\n", sigma);
#endif
    _mult_real(alpha,gamm,-delta);
#if DEBUG_SUMR ==1
    printf("alpha=%g,%g;\t\n",alpha.re,alpha.im);
#endif
    _complex_set(r_off,s*z_r,0);
    _add_assign_complex(r_off,alpha,phi);
#if DEBUG_SUMR ==1
    printf("r_off=%g,%g;\t\n",r_off.re,r_off.im);
#endif
    _complex_conj(tmpc, c);
    _mult_real(r_hat, tmpc, z_r);
    _add_assign_complex(r_hat, alpha, phi_hat);
#if DEBUG_SUMR ==1
    printf("r_hat=%g,%g;\t\n",r_hat.re, r_hat.im);
#endif
    tmpr = 1/(sqrt(_complex_square_norm(r_hat) + (sigma*sigma) ));
    _mult_real(tmpc, r_hat, tmpr);
    _complex_conj(c, tmpc);
#if DEBUG_SUMR ==1
    printf("c=%g,%g;\t\n",c.re,c.im);
#endif
    s=-sigma * tmpr;
#if DEBUG_SUMR ==1
    printf("s=%g;\t\n", s);
#endif
    _complex_set(r_diag, s*sigma, 0.);
    _diff_assign_complex(r_diag, c, r_hat);
#if DEBUG_SUMR ==1
    printf("r_diag=%g,%g;\t\n",r_diag.re,r_diag.im);
#endif
    _mult_real(tau, c, -tau_hat);
#if DEBUG_SUMR ==1
    printf("tau=%g,%g;\t\n",tau.re,tau.im);
#endif
    tau_hat *= s;
#if DEBUG_SUMR ==1
    printf("tau_hat=%g;\t\n", tau_hat);
#endif
    _div_complex(eta, tau, r_diag);
#if DEBUG_SUMR ==1
    printf("eta=%g,%g;\t\n",eta.re,eta.im);
#endif
    _div_complex(kappa, r_off, r_diag_old);
#if DEBUG_SUMR ==1
    printf("kappa=%g,%g;\t\n",kappa.re,kappa.im);
#endif
    zero_spinor_field(w, N);
    assign_add_mul_add_mul(w, p, tmp2, alpha, kappa, N);
#if DEBUG_SUMR ==1
    printf("w=%g;\t\n", square_norm(w));
#endif
    assign_add_mul(p, tmp2, lambda, N);
#if DEBUG_SUMR ==1
    printf("p=%g;\t\n", square_norm(p, N, 1));
#endif
    diff(tmp2, v, w, N);
#if DEBUG_SUMR ==1
    printf("w-v=%g;\t\n", square_norm(tmp2, N, 1));
#endif
    assign_add_mul(x, tmp2, eta, N);
#if DEBUG_SUMR ==1
    printf("x=%g;\t\n", square_norm(x, N, 1));
#endif
    
    if(sigma==0) {
      printf("Exit because Sigma = %g\n",sigma);
      finalize_solver(solver_field, 10);
      return(iteration);
    }
    /* Check whether the precision is reached ... */
    err = tau_hat * tau_hat;

    /* relax ap_eps_sq for adaptive precision (abuse of index_shift as a  "prudence factor") */ 
    if(switch_on_adaptive_precision == 1)    ap_eps_sq = (shift * eps_sq) / err;

#if DEBUG_SUMR ==1 
    tmpr = square_norm(x, N, 1);
    if(g_proc_id == g_stdio_proc) {
      printf("it, tau,sigma, ||x||^2: %d\t%g\t%g\t%g\n",iteration,err,sigma,tmpr); 
      fflush( stdout);
    }
#endif
    if ( (iteration%10) == 0 && g_proc_id == g_stdio_proc ) {
      printf("SUMR iteration= %d\t|res|^2= %g\n",iteration,err); 
      /*      fflush( stdout); */ 
    }

    if (err <= eps_sq) {
      assign(P, x, N);
      finalize_solver(solver_field, 10);
      return(iteration);
    }

    delta = delta * sigma;
#if DEBUG_SUMR ==1
    printf("delta=%g;\t\n", delta);
#endif
    _complex_conj(tmpc, gamm);
    _complex_conj(tmpc1, c);
    _mult_real(phi,tmpc, s / delta);
    _diff_assign_complex(phi, c, phi_hat);
#if DEBUG_SUMR ==1
    printf("phi=%g;\t\n", phi);
#endif

    _div_complex(lambda, phi, r_diag);
#if DEBUG_SUMR ==1
    printf("lambda=%g;\t\n", lambda);
#endif

    _div_real(tmpc1, tmpc1, delta);
    _mult_assign_complex(tmpc2, tmpc1, tmpc);
    _mult_real(phi_hat, phi_hat, s);
    _add_complex(phi_hat, tmpc2);
#if DEBUG_SUMR ==1
    printf("phi_hat=%g;\t\n", phi_hat);
#endif

    assign(tmp, u, N);
    assign_add_mul(tmp, v_til, gamm, N);
    mul_r(v, 1 / sigma,tmp, N);
#if DEBUG_SUMR ==1
    printf("v=%g;\t\n", square_norm(v, N, 1));
#endif

    mul(tmp, tmpc, v, N);
    assign_mul_add_r(v_til, sigma, tmp, N);
#if DEBUG_SUMR ==1
    printf("v_til=%g;\t\n", square_norm(v_til, N, 1));
    printf("############################\n");
#endif

    r_diag_old = r_diag;
  }
  finalize_solver(solver_field, 10);
  return(-1);
}



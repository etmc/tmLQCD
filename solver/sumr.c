/***********************************************************************
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
#include <complex.h>
#include "gamma.h"
#include "solver/eigenvalues.h"
#include "solver/sub_low_ev.h"
#include "operator/Dov_psi.h"
#include "solver_field.h"
#include "sumr.h"

#define DEBUG_SUMR 0

/* to be fixed somewhere else */


/* P output = solution , Q input = source */
int sumr(spinor * const P, spinor * const Q, const int max_iter,
	 double eps_sq){
  double sigma, delta, s, rho, zeta, z_r, tmpr, normsp, err, tau_hat;
  double ov_s, m_ov=0.;
  _Complex double phi, phi_hat, tau, lambda, c, gamm, alpha, eta, kappa;
  _Complex double r_hat, r_off, r_diag, r_diag_old, tmpc, tmpc1, tmpc2;
  int iteration;
  /* to be fixed somewhere else */
  const int N=VOLUME;
  spinor *x, *r, *p, *v, *v_til, *w, *u, *b, *tmp, *tmp2;
  spinor ** solver_field = NULL;
  printf("Starting SUMR!\n");
  /*        !!!!   INITIALIZATION    !!!! */
  init_solver_field(&solver_field, VOLUMEPLUSRAND, 10);
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

  phi_hat = 1 / delta;
  tau_hat = delta / rho;
  zero_spinor_field(p, N);
  phi = 0.0;
  s = 0.;
  lambda = 0.0;
  r_diag_old = 1.0;
  gamm = 1.0;
  sigma = 1.;
  c = 1.0;
  mul(v_til, phi_hat, r, N);
  assign(v, v_til, N);

#if DEBUG_SUMR ==1
  printf("delta=%g;\t phihat=%g;\t tauhat=%g;\t w=%g;\t p=%g;\t phi=%g;\t s=%g;\t lambda=%g;\t r_off=%g;\t r_off_old=%g;\t r_diag=%g;\t r_diag_old=%g;\t gamm=%g;\t sigma=%g;\t c=%g;\t v=%g;\t v_til=%g;\t ",
	 delta,cabs(phi_hat),tau_hat,square_norm(w),square_norm(p),
	 cabs(phi),s,cabs(lambda),cabs(r_off),cabs(r_off_old),cabs(r_diag),cabs(r_diag_old),cabs(gamm),sigma,cabs(c),
	 square_norm(v),square_norm(v_til));
#endif

  if(ov_cheby_coef==NULL) calculateOverlapPolynomial();

  /* main loop */
  for(iteration = 0; iteration < max_iter; iteration++) {
    Q_over_sqrt_Q_sqr(tmp, ov_cheby_coef, ov_n_cheby, v, ev_qnorm, ev_minev);
    gamma5(u, tmp, N);
#if DEBUG_SUMR ==1
    printf("u=%g;\t\n", square_norm(u));
#endif
    gamm = scalar_prod(v_til, u, N, 1);
    gamm = -(gamm);
#if DEBUG_SUMR ==1
    printf("gamm=%g,%g;\t\n",creal(gamm),cimag(gamm));
#endif
    sigma= sqrt((1 - cabs(gamm))*(1 + cabs(gamm)));
#if DEBUG_SUMR ==1
    printf("sigma=%g;\t\n", sigma);
#endif
    alpha = -gamm * delta;
#if DEBUG_SUMR ==1
    printf("alpha=%g,%g;\t\n",creal(alpha),cimag(alpha));
#endif
    r_off = s*z_r;
    r_off += (alpha) * (phi);
#if DEBUG_SUMR ==1
    printf("r_off=%g,%g;\t\n",creal(r_off),cimag(r_off));
#endif
    tmpc = conj(c);
    r_hat = (tmpc) * (z_r);
    r_hat += (alpha) * (phi_hat);
#if DEBUG_SUMR ==1
    printf("r_hat=%g,%g;\t\n",creal(r_hat), cimag(r_hat));
#endif
    tmpr = 1/(sqrt(creal(r_hat * conj(r_hat)) + (sigma*sigma)));
    tmpc = (r_hat) * (tmpr);
    c = conj(tmpc);
#if DEBUG_SUMR ==1
    printf("c=%g,%g;\t\n",creal(c),cimag(c));
#endif
    s=-sigma * tmpr;
#if DEBUG_SUMR ==1
    printf("s=%g;\t\n", s);
#endif
    r_diag = s*sigma;
    r_diag -= c * r_hat;
#if DEBUG_SUMR ==1
    printf("r_diag=%g,%g;\t\n",creal(r_diag),cimag(r_diag));
#endif
    tau = -c * tau_hat;
#if DEBUG_SUMR ==1
    printf("tau=%g,%g;\t\n",creal(tau),cimag(tau));
#endif
    tau_hat *= s;
#if DEBUG_SUMR ==1
    printf("tau_hat=%g;\t\n", tau_hat);
#endif
    eta = tau / r_diag;
#if DEBUG_SUMR ==1
    printf("eta=%g,%g;\t\n",creal(eta),cimag(eta));
#endif
    kappa = r_off / r_diag_old;
#if DEBUG_SUMR ==1
    printf("kappa=%g,%g;\t\n",creal(kappa),cimag(kappa));
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

#if DEBUG_SUMR ==1 
    tmpr = square_norm(x, N, 1);
    if(g_proc_id == g_stdio_proc) {
      printf("it, tau,sigma, ||x||^2: %d\t%g\t%g\t%g\n",iteration,err,sigma,tmpr); 
      fflush( stdout);
    }
#endif
    if ((iteration%10) == 0 && g_proc_id == g_stdio_proc ) {
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
    tmpc = conj(gamm);
    tmpc1 = conj(c);
    phi = (tmpc) * (s / delta);
    phi -= (c) * (phi_hat);
#if DEBUG_SUMR ==1
    printf("phi=%g;\t\n", phi);
#endif

    lambda = (phi) / (r_diag);
#if DEBUG_SUMR ==1
    printf("lambda=%g;\t\n", lambda);
#endif

    tmpc1 = (tmpc1) / (delta);
    tmpc2 = (tmpc1) * (tmpc);
    phi_hat = (phi_hat) * (s);
    phi_hat += tmpc2;
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



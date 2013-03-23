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
 * This is an implementation of bicgstab(l)
 * corresponding to the paper of G. L.G. Sleijpen and
 * D.R. Fokkema
 * Transactions on Numerical Analysis
 * Volume1, pp. 11-32, 1993
 *
 * Author: Carsten Urbach
 *         urbach@physik.fu-berlin.de
 *
 *************************************************************/

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
#include "solver_field.h"
#include "bicgstabell.h"

int bicgstabell(spinor * const x0, spinor * const b, const int max_iter, 
		double eps_sq, const int rel_prec, const int _l, const int N, matrix_mult f) {

  double err;
  int i, j, k, l;
  double rho0, rho1, beta, alpha, omega, gamma0 = 0., squarenorm;
  spinor * r[5], * u[5], * r0_tilde, * x;
  double tau[5][5], gamma[25], gammap[25], gammapp[25], sigma[25];
  spinor ** solver_field = NULL;
  const int nr_sf = 2*(_l+1)+2;

  l = _l;
  k = -l;

  if(N == VOLUME) {
    init_solver_field(&solver_field, VOLUMEPLUSRAND, nr_sf);
  }
  else {
    init_solver_field(&solver_field, VOLUMEPLUSRAND/2, nr_sf);
  }
  r0_tilde = solver_field[0];
  for(i = 0; i <= l; i++){
    r[i] = solver_field[2+2*i];
    u[i] = solver_field[3+2*i];
  }

  x = x0; 
  assign(u[0], b, N);
  f(r0_tilde, x);
  diff(r[0], u[0], r0_tilde, N);
  zero_spinor_field(solver_field[1], N);
  assign(r0_tilde, r[0], N);
  squarenorm = square_norm(b, N, 1);

  rho0 = 1.;
  alpha = 0.;
  omega = 1.;
  err = square_norm(r0_tilde, N, 1);
  while( k < max_iter && (((err > eps_sq) && (rel_prec == 0)) 
			  || ((err > eps_sq*squarenorm) && (rel_prec == 1)) 
			  )) {
    k+=l;

    /* The BiCG part */

    rho0 *= -omega;
    for(j = 0; j < l; j++) {
      rho1 = scalar_prod_r(r[j], r0_tilde, N, 1);
      beta = (rho1/rho0);
      beta *= alpha; 
      rho0 = rho1;
      for(i = 0; i <= j; i++) {
	/* u_i = r_i - \beta u_i */
	assign_mul_add_r(u[i], -beta, r[i], N);
      }
      f(u[j+1], u[j]);
      gamma0 = scalar_prod_r(u[j+1], r0_tilde, N, 1);
      alpha = rho0/gamma0;
      /* r_i = r_i - \alpha u_{i+1} */
      for(i = 0; i <= j; i++) {
	assign_add_mul_r(r[i], u[i+1], -alpha, N);
      }
      f(r[j+1], r[j]);
      /* x = x + \alpha u_0 */
      assign_add_mul_r(x, u[0], alpha, N);
      err = square_norm(r[j+1], N, 1);
      if(g_proc_id == 0 && g_debug_level > 2) {printf("%d %d err = %e\n", k, j, err);fflush(stdout);}
    }

    /* The MR part */

    for(j = 1; j <= l; j++){
      for(i = 1; i < j; i++){
	tau[i][j] = scalar_prod_r(r[j], r[i], N, 1)/sigma[i];
	assign_add_mul_r(r[j], r[i], -tau[i][j], N);
      }
      sigma[j] = scalar_prod_r(r[j], r[j], N, 1);
      gammap[j] = scalar_prod_r(r[0], r[j], N, 1)/sigma[j];
    }
    gamma[l] = gammap[l];
    omega = gamma[l];
    for(j = l-1; j > 0; j--) {
      gamma[j] = gammap[j];
      for(i = j+1; i <= l; i++) {
	gamma[j] -= (tau[j][i]*gamma[i]);
      }
    }
    for(j = 1; j < l; j++) {
      gammapp[j] = gamma[j+1];
      for(i = j+1; i < l; i++){
	gammapp[j] += (tau[j][i]*gamma[i+1]);
      }
    }
    assign_add_mul_r(x, r[0], gamma[1], N);
    assign_add_mul_r(r[0], r[l], -gammap[l], N);
    for(j = 1; j < l; j++){
      assign_add_mul_r(x, r[j], gammapp[j], N);
      assign_add_mul_r(r[0], r[j], -gammap[j], N);
    }
    assign_add_mul_r(u[0], u[l], -gamma[l], N);
    for(j = 1; j < l; j++){
      assign_add_mul_r(u[0], u[j], -gamma[j], N);
    }
    err = square_norm(r[0], N, 1);
    if(g_proc_id == 0 && g_debug_level > 0){
      printf(" BiCGstabell iterated %d %d, %e rho0 = %e, alpha = %e, gamma0= %e\n", l, k, err, rho0, alpha, gamma0);
      fflush( stdout );
    }
  }
  finalize_solver(solver_field, nr_sf);
  if(k == max_iter) return(-1);
  return(k);
}

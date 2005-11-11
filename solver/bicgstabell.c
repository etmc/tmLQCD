/* $Id$ */
/*************************************************************
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


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "global.h"
#include "su3.h"
#include "linalg_eo.h"
#include "start.h"
#include "solver/matrix_mult_typedef.h"
#include "bicgstabell.h"

#ifdef _SOLVER_OUTPUT
#define _SO(x) x
#else
#define _SO(x)
#endif

int bicgstabell(spinor * const x0, spinor * const b, const int max_iter, 
		double eps_sq, const int _l, matrix_mult f) {

  double err;
  int i, j, k, l, N=VOLUME/2;
  double rho0, rho1, beta, alpha, omega, gamma0;
  spinor * r[5], * u[5], * r0_tilde, * u0, * x;
  double tau[5][5], gamma[25], gammap[25], gammapp[25], sigma[25];


  l = _l;
  k = -l;

/*   init_solver_field(2*(l+1)+2); */
  r0_tilde = g_spinor_field[DUM_SOLVER];
  u0 = g_spinor_field[DUM_SOLVER+1];
  for(i = 0; i <= l; i++){
    r[i] = g_spinor_field[DUM_SOLVER+2+2*i];
    u[i] = g_spinor_field[DUM_SOLVER+3+2*i];
  }

  x = x0; 
  assign(u[0], b, N);
  f(r0_tilde, x);
  diff(r[0], u[0], r0_tilde, N);
  zero_spinor_field(g_spinor_field[DUM_SOLVER+1]);
  assign(r0_tilde, r[0], N);

  rho0 = 1.;
  alpha = 0.;
  omega = 1.;
  err = square_norm(r0_tilde, N);

  while( k < max_iter && err > eps_sq) {
    k+=l;

    /* The BiCG part */

    rho0 *= -omega;
    for(j = 0; j < l; j++) {
      rho1 = scalar_prod_r(r[j], r0_tilde, N);
      beta = (rho1/rho0);
      beta *= alpha; 
      rho0 = rho1;
      for(i = 0; i <= j; i++) {
	/* u_i = r_i - \beta u_i */
	assign_mul_add_r(u[i], -beta, r[i], N);
      }
      f(u[j+1], u[j]);
      gamma0 = scalar_prod_r(u[j+1], r0_tilde, N);
      alpha = rho0/gamma0;
      /* r_i = r_i - \alpha u_{i+1} */
      for(i = 0; i <= j; i++) {
	assign_add_mul_r(r[i], u[i+1], -alpha, N);
      }
      f(r[j+1], r[j]);
      /* x = x + \alpha u_0 */
      assign_add_mul_r(x, u[0], alpha, N);
    }

    /* The MR part */

    for(j = 1; j <= l; j++){
      for(i = 1; i < j; i++){
	tau[i][j] = scalar_prod_r(r[j], r[i], N)/sigma[i];
	assign_add_mul_r(r[j], r[i], -tau[i][j], N);
      }
      sigma[j] = scalar_prod_r(r[j], r[j], N);
      gammap[j] = scalar_prod_r(r[0], r[j], N)/sigma[j];
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
    err = square_norm(r[0], N);
    _SO(if(g_proc_id == 0){printf(" Iterated %d %d, %e\n", l, k, err);fflush( stdout );});
  }
  if(k == max_iter) return(-1);
  return(k);
}

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
#include "tm_operators.h"
#include "clover_eo.h"
#include "bicgstabell.h"


#ifdef _SOLVER_OUTPUT
#define _SO(x) x
#else
#define _SO(x)
#endif

int bicgstabell(const int x0, const int b, const int max_iter, 
		double eps_sq, const int _l, const double q_off) {

  double err;
  int i, j, k, l;
  double rho0, rho1, beta, alpha, omega, gamma0;
  int r[5], u[5], r0_tilde, u0, x;
  double tau[5][5], gamma[25], gammap[25], gammapp[25], sigma[25];


  l = _l;
  k = -l;
  if(_l > 4){
    l = 4;
    k = -4;
  }

  r0_tilde = DUM_SOLVER;
  u0 = DUM_SOLVER+1;
  x = DUM_SOLVER+2;
  for(i = 0; i <= l; i++){
    r[i] = DUM_SOLVER+3+2*i;
    u[i] = DUM_SOLVER+4+2*i;
  }


  x = x0; 
  gamma5(u[0], b); 
  if(g_use_clover_flag == 1){
    M_psi(r0_tilde, x, q_off); 
  }
  else {
    Mtm_plus_psi(r0_tilde, x);
  }
  diff(r[0], u[0], r0_tilde, VOLUME/2);
  zero_spinor_field(u[0]);
  assign(r0_tilde, r[0], VOLUME/2);

  rho0 = 1.;
  alpha = 0.;
  omega = 1.;
  err = square_norm(r0_tilde, VOLUME/2);

  while( k < max_iter && err > eps_sq) {
    k+=l;

    /* The BiCG part */

    rho0 *= -omega;
    for(j = 0; j < l; j++) {
      rho1 = scalar_prod_r(r[j], r0_tilde, VOLUME/2);
      beta = (rho1/rho0);
      beta *= alpha; 
      rho0 = rho1;
      for(i = 0; i <= j; i++) {
	/* u_i = r_i - \beta u_i */
	assign_mul_add_r(u[i], -beta, r[i], VOLUME/2);
      }
      if(g_use_clover_flag == 1){
	M_psi(u[j+1], u[j], q_off); 
      }
      else {
	Mtm_plus_psi(u[j+1], u[j]);
      }
      gamma0 = scalar_prod_r(u[j+1], r0_tilde, VOLUME/2);
      alpha = rho0/gamma0;
      /* r_i = r_i - \alpha u_{i+1} */
      for(i = 0; i <= j; i++) {
	assign_add_mul(r[i], -alpha, u[i+1], VOLUME/2);
      }
      if(g_use_clover_flag == 1){
	M_psi(r[j+1], r[j], q_off); 
      }
      else {
	Mtm_plus_psi(r[j+1], r[j]);
      }
      /* x = x + \alpha u_0 */
      assign_add_mul(x, alpha, u[0], VOLUME/2);
    }

    /* The MR part */

    for(j = 1; j <= l; j++){
      for(i = 1; i < j; i++){
	tau[i][j] = scalar_prod_r(r[j], r[i], VOLUME/2)/sigma[i];
	assign_add_mul(r[j], -tau[i][j], r[i], VOLUME/2);
      }
      sigma[j] = scalar_prod_r(r[j], r[j], VOLUME/2);
      gammap[j] = scalar_prod_r(r[0], r[j], VOLUME/2)/sigma[j];
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
    assign_add_mul(x, gamma[1], r[0], VOLUME/2);
    assign_add_mul(r[0], -gammap[l], r[l], VOLUME/2);
    for(j = 1; j < l; j++){
      assign_add_mul(x, gammapp[j], r[j], VOLUME/2);
      assign_add_mul(r[0], -gammap[j], r[j], VOLUME/2);
    }
    assign_add_mul(u[0], -gamma[l], u[l], VOLUME/2);
    for(j = 1; j < l; j++){
      assign_add_mul(u[0], -gamma[j], u[j], VOLUME/2);
    }
    err = square_norm(r[0], VOLUME/2);
    _SO(if(g_proc_id == 0){printf(" Iterated %d %d, %e\n", l, k, err);fflush( stdout );});
  }
  if(k == max_iter) return(-1);
  return(k);
}

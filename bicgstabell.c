/* $Id$ */

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

int bicgstabell(const int P, const int Q, const int max_iter, 
		double eps_sq, const int _l, const double q_off) {

  double err;
  int i, j, k, l;
  double rho0, rho1, beta, alpha, omega, gamma0;
  int r[5], u[5], r0, u0, x;
  double tau[5][5], gamma[25], gammap[25], gammapp[25], sigma[25];


  l = _l;
  k = -l;
  if(_l > 4){
    l = 4;
    k = -4;
  }

  r0 = DUM_SOLVER;
  u0 = DUM_SOLVER+1;
  x = DUM_SOLVER+2;
  for(i = 0; i <= l; i++){
    r[i] = DUM_SOLVER+3+2*i;
    u[i] = DUM_SOLVER+4+2*i;
  }


  assign(x, P, VOLUME/2);
  if(g_use_clover_flag == 1){
    M_psi(r0, P, q_off); 
  }
  else {
    Mtm_plus_psi(r0, P);
  }
  diff(r[0], Q, r0, VOLUME/2);
  zero_spinor_field(u[0]);
  assign(r0, r[0], VOLUME/2);

  rho0 = 1.;
  alpha = 0.;
  omega = 1.;
  err = square_norm(r0, VOLUME/2);

  while( k < max_iter && err > eps_sq) {
    k+=l;
    if(g_proc_id == 0){
      printf(" %d, %e\n", k, err);
    }

    /* The BiCG part */

    rho0 *= -omega;
    for(j = 0; j < l; j++) {
      rho1 = scalar_prod_r(r[j], r0, VOLUME/2);
      beta = alpha*(rho1/rho0);
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
      gamma0 = scalar_prod_r(u[j+1], r0, VOLUME/2);
      alpha = rho0/gamma0;
      for(i = 0; i <= j; i++) {
	assign_add_mul(r[i], -alpha, u[i+1], VOLUME/2);
      }
      if(g_use_clover_flag == 1){
	M_psi(r[j+1], r[j], q_off); 
      }
      else {
	Mtm_plus_psi(r[j+1], r[j]);
      }
      assign_add_mul(x, alpha, u[0], VOLUME/2);
    }

    /* The MR part */

    for(j = 1; j <= l; j++){
      for(i = 1; i < j; i++){
	tau[i][j] = 1./sigma[i]*scalar_prod_r(r[j], r[i], VOLUME/2);
	assign_add_mul(r[j], -tau[i][j], r[i], VOLUME/2);
      }
      sigma[j] = scalar_prod_r(r[j], r[j], VOLUME/2);
      gammap[j] = 1./sigma[j]*scalar_prod_r(r[0], r[j], VOLUME/2);
    }
    gamma[l] = gammap[l];
    omega = gamma[l];
    for(j = l-1; j > 0; j--) {
      gamma[j] = gammap[j];
      for(i = j+1; i < l; i++) {
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
  }
  assign(P, x, VOLUME/2);
  if(k == max_iter) return(-1);
  return(k);
}

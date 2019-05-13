/***********************************************************************
 * Copyright (C) 2016 Bartosz Kostrzewa, Florian Burger
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
 **********************
 * rg_mixed_cg_her_nd *
 **********************
 *
 * Mixed precision solver which uses true reliable updates and has a double
 * precision fail-safe mechanism. The Polak-Ribiere computation of beta is
 * implemented but currently not used because the extra scalar product is
 * more expensive than the gain from the self-stabilisation as far as has 
 * been tested.
 *
 * in:
 *   Q: source
 * input:
 *   P: result
 *
 * POSSIBLE IMPROVEMENTS
 * There are still quite a few things that can be tried to make it better,
 * the most significant of which would be to guide the search direction
 * using the previous one upon restart. However, it seems that for the number
 * non-zero entries in the Dirac operator and usual lattice sizes, the
 * requisite projection 
 *
 *   p' = r - <r,Ap>/<p,Ap> p
 *
 * cannot be computed with sufficient precision in 64 bit arithmetic. It should
 * be noted that for L < 24 in general, this does work and produces
 * a mixed solver which converges at the same rate as a double solver, but it's
 * not generally useable... For point sources, it also works for larger lattice 
 * volumes. Might be introduced as an optional mode in the future with some
 * fail-safe mechanism which detects if the search direction begins to diverge.
 *
 **************************************************************************/

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
#include "operator/tm_operators_32.h"
#include "operator/clovertm_operators_32.h"
#include "solver/matrix_mult_typedef_nd.h"
#include "solver/solver_params.h"
#include "read_input.h"

#include "solver_field.h"
#include "solver/rg_mixed_cg_her.h"
#include "gettime.h"

static void output_flops(const double seconds, const unsigned int N, const unsigned int iter_out, 
                  const unsigned int iter_in_sp, const unsigned int iter_in_dp, const double eps_sq);

static inline unsigned int inner_loop_high(spinor * const x_up, spinor * const x_dn, 
                                           spinor * const p_up, spinor * const p_dn,
                                           spinor * const q_up, spinor * const q_dn,
                                           spinor * const r_up, spinor * const r_dn, 
                                           double * const rho1, const double delta,
                                           matrix_mult_nd f, const double eps_sq, const unsigned int N, const unsigned int iter, const unsigned int max_iter ){

  static double alpha, beta, rho, rhomax;
  unsigned int j = 0;

  rho = *rho1;
  rhomax = *rho1;

  /* break out of inner loop if iterated residual goes below some fraction of the maximum observed
  * iterated residual since the last update or if the target precision has been reached 
  * enforce convergence more strictly by a factor of 1.3 to avoid unnecessary restarts 
  * if the real residual is still a bit too large */
  while( rho > delta*rhomax && j+iter <= max_iter ){
    ++j;
    f(q_up,q_dn,p_up,p_dn);
    alpha = rho/( scalar_prod_r(p_up,q_up,N,1) + scalar_prod_r(p_dn,q_dn,N,1) );
    assign_add_mul_r(x_up, p_up, alpha, N); assign_add_mul_r(x_dn, p_dn, alpha, N);
    assign_add_mul_r(r_up, q_up, -alpha, N); assign_add_mul_r(r_dn, q_dn, -alpha, N);
    rho = ( square_norm(r_up,N,1) + square_norm(r_dn,N,1) );
    beta = rho / *rho1;
    *rho1 = rho;
    assign_mul_add_r(p_up, beta, r_up, N); assign_mul_add_r(p_dn, beta, r_dn, N);
    
    if( 1.3*rho < eps_sq ) break;
    if( rho > rhomax ) rhomax = rho;
    
    if(g_debug_level > 2 && g_proc_id == 0) {
      printf("DP_inner CG: %d res^2 %g\t\n", j+iter, rho);
    }
  }

  return j;
}

static inline unsigned int inner_loop(spinor32 * const x_up, spinor32 * const x_dn, 
                                      spinor32 * const p_up, spinor32 * const p_dn, 
                                      spinor32 * const q_up, spinor32 * const q_dn,
                                      spinor32 * const r_up, spinor32 * const r_dn,
                                      float * const rho1, const float delta,
                                      matrix_mult_nd32 f32, const float eps_sq, const unsigned int N, const unsigned int iter, const unsigned int max_iter,
                                      float alpha, float beta, MCG_PIPELINED_TYPE pipelined, MCG_PR_TYPE pr ){

  static float rho, rhomax, pro;
  unsigned int j = 0;

  rho = *rho1;
  rhomax = *rho1;

  if(pipelined==MCG_NO_PIPELINED){
    /* break out of inner loop if iterated residual goes below some fraction of the maximum observed
    * iterated residual since the last update */ 
    while( rho > delta*rhomax && j+iter <= max_iter ){
      ++j;
      f32(q_up,q_dn,p_up,p_dn);
      pro = ( scalar_prod_r_32(p_up,q_up,N,1) + scalar_prod_r_32(p_dn,q_dn,N,1) );
      alpha = rho/pro;
      assign_add_mul_r_32(x_up, p_up, alpha, N); assign_add_mul_r_32(x_dn, p_dn, alpha, N);
      assign_add_mul_r_32(r_up, q_up, -alpha, N); assign_add_mul_r_32(r_dn, q_dn, -alpha, N);
      rho = ( square_norm_32(r_up,N,1) + square_norm_32(r_dn,N,1) );
      // Polak-Ribiere computation of beta, claimed to be self-stabilising, positive effect so far not observed or required
      if(pr==MCG_PR){
        beta = alpha*(alpha*(square_norm_32(q_up,N,1)+square_norm_32(q_dn,N,1)) - pro) / *rho1;
      }else{
        beta = rho / *rho1;
      }
      *rho1 = rho;
      assign_mul_add_r_32(p_up, beta, r_up, N); assign_mul_add_r_32(p_dn, beta, r_dn, N);
      if(g_debug_level > 2 && g_proc_id == 0) {
        printf("SP_inner CG_ND: %d res^2 %g\t\n", j+iter, rho);
      }
       /* enforce convergence more strictly by a factor of 1.3 to avoid unnecessary restarts 
       * if the real residual is still a bit too large */
      if( 1.3*rho < eps_sq ) break;
      if( rho > rhomax ) rhomax = rho;
    }
  }else{
    // pipelined cg requires one more scalar product but may allow optimisations to be made
    // e.g.: one could do the collective communication for sqrnrm(r) while other stuff is being computed
    // it is also self-initialising (alpha=0, beta=0 will work)
    while( rho > delta*rhomax && j+iter <= max_iter ){
      ++j;
      assign_add_mul_r_32(x_up, p_up, alpha, N); assign_add_mul_r_32(x_dn, p_dn, alpha, N);
      assign_add_mul_r_32(r_up, q_up, -alpha, N); assign_add_mul_r_32(r_dn, q_dn, -alpha, N);
      assign_mul_add_r_32(p_up, beta, r_up, N); assign_mul_add_r_32(p_dn, beta, r_dn, N);
      f32(q_up,q_dn,p_up,p_dn);
  
      rho = ( square_norm_32(r_up,N,1) + square_norm_32(r_dn,N,1) );
      pro = ( scalar_prod_r_32(p_up,q_up,N,1) + scalar_prod_r_32(p_dn,q_dn,N,1) );
      alpha = rho/pro;
      if(pr==MCG_PR){
        beta = alpha*(alpha*(square_norm_32(q_up,N,1)+square_norm_32(q_dn,N,1))-pro)/rho;
      }else{
        beta = rho/ *rho1;
      }
      *rho1=rho;

      if(g_debug_level > 2 && g_proc_id == 0) {
        printf("SP_inner CG_ND: %d res^2 %g\t\n", j+iter, rho);
      }
      if( 1.3*rho < eps_sq ) break;
      if( rho > rhomax ) rhomax = rho;
    }
  }

  return j;
}


/* P output = solution , Q input = source */
int rg_mixed_cg_her_nd(spinor * const P_up, spinor * const P_dn, spinor * const Q_up, spinor * const Q_dn, 
                       solver_params_t solver_params, const int max_iter, const double eps_sq, const int rel_prec,
                       const int N, matrix_mult_nd f, matrix_mult_nd32 f32) {

  int iter_in_sp = 0, iter_in_dp = 0, iter_out = 0;
  float rho_sp, delta = solver_params.mcg_delta;
  double beta_dp, rho_dp;
  double sourcesquarenorm, guesssquarenorm, target_eps_sq;

  spinor *xhigh_up, *xhigh_dn, *rhigh_up, *rhigh_dn, *qhigh_up, *qhigh_dn, *phigh_up, *phigh_dn;
  spinor32 *x_up, *x_dn, *p_up, *p_dn, *q_up, *q_dn, *r_up, *r_dn;

  spinor ** solver_field = NULL;
  spinor32 ** solver_field32 = NULL;  
  const int nr_sf = 8;
  const int nr_sf32 = 8;
  
  int high_control = 0;

  double atime, etime, flops;
  
  if(N == VOLUME) {
    init_solver_field(&solver_field, VOLUMEPLUSRAND, nr_sf);    
    init_solver_field_32(&solver_field32, VOLUMEPLUSRAND, nr_sf32);
  }
  else {
    init_solver_field(&solver_field, VOLUMEPLUSRAND/2, nr_sf);
    init_solver_field_32(&solver_field32, VOLUMEPLUSRAND/2, nr_sf32);    
  }

  atime = gettime();

  // we could get away with using fewer fields, of course
  phigh_up = solver_field[7]; phigh_dn = solver_field[6];
  xhigh_up = solver_field[5]; xhigh_dn = solver_field[4];
  rhigh_up = solver_field[3]; rhigh_dn = solver_field[2];
  qhigh_up = solver_field[1]; qhigh_dn = solver_field[0];

  x_up = solver_field32[7]; x_dn = solver_field32[6];
  r_up = solver_field32[5]; r_dn = solver_field32[4];
  p_up = solver_field32[3]; p_dn = solver_field32[2];
  q_up = solver_field32[1]; q_dn = solver_field32[0];

  // we always want to apply the full precision operator in double
  int save_sloppy = g_sloppy_precision_flag;
  g_sloppy_precision_flag = 0;

  sourcesquarenorm = ( square_norm(Q_up,N,1) + square_norm(Q_dn,N,1) );
  if( rel_prec == 1 ) {
    target_eps_sq = eps_sq*sourcesquarenorm;
    if(g_debug_level > 0 && g_proc_id==0) 
      printf("#RG_Mixed CG_ND: Using relative precision! eps_sq: %.6g target_eps_sq: %.6g \n",eps_sq,target_eps_sq);
  }else{
    target_eps_sq = eps_sq;
  }
  
  // compute the maximum number of outer iterations based on the expected reduction
  // of the residual at each run of the inner solver
  int N_outer = (int)ceil(log10( sourcesquarenorm*delta/target_eps_sq ));
  if(g_debug_level > 0 && g_proc_id==0) 
    printf("#RG_Mixed CG_ND: N_outer: %d \n", N_outer);
  
  zero_spinor_field_32(x_up,N); zero_spinor_field_32(x_dn,N);

  if(solver_params.use_initial_guess == 0) {
    assign(phigh_up,Q_up,N); assign(phigh_dn,Q_dn,N);
    assign(rhigh_up,Q_up,N); assign(rhigh_dn,Q_dn,N);
    rho_dp = sourcesquarenorm;
  } else {
    // computing initial guess
    f(rhigh_up,rhigh_dn,P_up,P_dn);
    diff(rhigh_up,Q_up,rhigh_up,N); diff(rhigh_dn,Q_dn,rhigh_dn,N);
    assign(phigh_up,rhigh_up,N); assign(phigh_dn,rhigh_dn,N);
    rho_dp = ( square_norm(rhigh_up,N,1) + square_norm(rhigh_dn,N,1) );
  }

  assign_to_32(r_up,rhigh_up,N); assign_to_32(r_dn,rhigh_dn,N);
  rho_sp = rho_dp;
  assign_32(p_up,r_up,N); assign_32(p_dn,r_dn,N);

  iter_in_sp += inner_loop(x_up, x_dn, p_up, p_dn, q_up, q_dn, r_up, r_dn, &rho_sp, delta, 
                           f32, (float)target_eps_sq, 
                           N, iter_out+iter_in_sp+iter_in_dp, max_iter, 0.0, 0.0, MCG_NO_PIPELINED, MCG_NO_PR);

  for(iter_out = 1; iter_out < N_outer; ++iter_out ) {

    // prepare for defect correction
    // update high precision solution 
    if(high_control==0) {
      // accumulate solution (sp -> dp) 
      addto_32(P_up,x_up,N); addto_32(P_dn,x_dn,N);
      // compute real residual
      f(qhigh_up,qhigh_dn,P_up,P_dn);
      diff(rhigh_up,Q_up,qhigh_up,N); diff(rhigh_dn,Q_dn,qhigh_dn,N);
      beta_dp = 1/rho_dp;
      rho_dp = ( square_norm(rhigh_up,N,1) + square_norm(rhigh_dn,N,1) );
      beta_dp *= rho_dp;
    }
   
    // the iteration limit was reached in the previous iteration, let's try to save the day using double precision
    if( high_control==1 ) {
      assign(phigh_up,rhigh_up,N); assign(phigh_dn,rhigh_dn,N);
      zero_spinor_field(xhigh_up,N); zero_spinor_field(xhigh_dn,N);
      beta_dp = 1/rho_dp;
      iter_in_dp += inner_loop_high(xhigh_up, xhigh_dn, phigh_up, phigh_dn,
                                    qhigh_up, qhigh_dn, rhigh_up, rhigh_dn, &rho_dp, delta, f, 
                                    target_eps_sq, N, iter_out+iter_in_sp+iter_in_dp, max_iter);
      rho_sp = rho_dp;
      // accumulate solution
      add(P_up,P_up,xhigh_up,N); add(P_dn, P_dn, xhigh_dn, N);
      // compute real residual
      f(qhigh_up, qhigh_dn, P_up, P_dn);
      diff(rhigh_up,Q_up,qhigh_up,N); diff(rhigh_dn,Q_dn,qhigh_dn,N);
      rho_dp = ( square_norm(rhigh_up,N,1) + square_norm(rhigh_dn,N,1) );
      beta_dp *= rho_dp;
    }

    if(g_debug_level > 2 && g_proc_id == 0) {
      printf("RG_mixed CG_ND last inner residue:       %17g\n", rho_sp);
      printf("RG_mixed CG_ND true residue:             %6d %10g\n", iter_in_sp+iter_in_dp+iter_out, rho_dp);
      printf("RG_mixed CG_ND residue reduction factor: %6d %10g\n", iter_in_sp+iter_in_dp+iter_out, beta_dp); fflush(stdout);
    }

    if( rho_dp <= target_eps_sq || (iter_in_sp+iter_in_dp+iter_out) >= max_iter ) {
      etime = gettime();
      output_flops(etime-atime, N, iter_out, iter_in_sp, iter_in_dp, eps_sq);
      
      g_sloppy_precision_flag = save_sloppy;
      finalize_solver(solver_field, nr_sf);
      finalize_solver_32(solver_field32, nr_sf32); 
      if( (iter_in_sp+iter_in_dp+iter_out) >= max_iter ){
        return(-1);
      } else {
        return(iter_in_sp+iter_in_dp+iter_out);
      }
    }

    // if it seems like we're stuck and reaching the iteration limit, we skip this correction and proceed in full precision above
    if( iter_out >= (N_outer-2) ){
      if(g_proc_id==0) printf("RG_mixed CG_ND: Reaching iteration limit, switching to DP!\n");
      high_control = 1;
      continue;
    }else{
      // correct defect
      assign_to_32(r_up,rhigh_up,N); assign_to_32(r_dn,rhigh_dn,N);
      rho_sp = rho_dp; // not sure if it's fine to truncate this or whether one should calculate it in SP directly, it seems to work fine though
      assign_32(p_up,r_up,N); assign_32(p_dn,r_dn,N);
    }

    zero_spinor_field_32(x_up,N); zero_spinor_field_32(x_dn,N);
    iter_in_sp += inner_loop(x_up, x_dn, p_up, p_dn, q_up, q_dn, r_up, r_dn, &rho_sp, delta, f32, (float)target_eps_sq, 
                             N, iter_out+iter_in_sp+iter_in_dp, max_iter, 0.0, 0.0, MCG_NO_PIPELINED, MCG_NO_PR);
  }
  
  // convergence failure...
  g_sloppy_precision_flag = save_sloppy;
  finalize_solver(solver_field, nr_sf);
  finalize_solver_32(solver_field32, nr_sf32);
  return -1; 
}

void output_flops(const double seconds, const unsigned int N, const unsigned int iter_out, const unsigned int iter_in_sp, const unsigned int iter_in_dp, const double eps_sq){
  double flops;
  // TODO: compute real number of flops...
  int total_it = iter_in_sp+iter_in_dp+iter_out;
  if(g_debug_level > 0 && g_proc_id == 0) {
    printf("# RG_mixed CG_ND: iter_out: %d iter_in_sp: %d iter_in_dp: %d\n",iter_out,iter_in_sp,iter_in_dp);
  	if(N != VOLUME){
  	  /* 2 A + 2 Nc Ns + N_Count ( 2 A + 10 Nc Ns ) */
  	  /* 2*1608.0 because the linalg is over VOLUME/2 */
  	  flops = 2*(2*(2*1608.0+2*3*4) + 2*3*4 + total_it*(2.*(2*1608.0+2*3*4) + 10*3*4))*N/1.0e6f;
  	}
  	else{
  	  /* 2 A + 2 Nc Ns + N_Count ( 2 A + 10 Nc Ns ) */
  	  flops = 2*(2*(1608.0+2*3*4) + 2*3*4 + total_it*(2.*(1608.0+2*3*4) + 10*3*4))*N/1.0e6f;      
  	}
  	printf("#RG_mixed CG_ND: iter: %d eps_sq: %1.4e t/s: %1.4e\n", total_it, eps_sq, seconds); 
    printf("# FIXME: note the following flop counts are wrong! Consider only the time to solution!\n");
  	printf("#RG_mixed CG_ND: flopcount (for e/o tmWilson only): t/s: %1.4e mflops_local: %.1f mflops: %.1f\n", 
  	       seconds, flops/(seconds), g_nproc*flops/(seconds));
  }      
}

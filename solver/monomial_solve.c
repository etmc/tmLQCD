/***********************************************************************
 *
 * Copyright (C) 2014 Florian Burger
 *               2017 Bartosz Kostrzewa
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
 *  
 * File: monomial_solve.c
 *
 * solver wrapper for monomials
 *
 * The externally accessible functions are
 *
 *
 *   int solve_degenerate(spinor * const P, spinor * const Q, const int max_iter, 
 *                       double eps_sq, const int rel_prec, const int N, matrix_mult f)
 *
 *   int solve_mms_tm(spinor ** const P, spinor * const Q,
 *                    solver_params_t * solver_params)  
 *
 *   int solve_mms_nd(spinor ** const Pup, spinor ** const Pdn, 
 *                    spinor * const Qup, spinor * const Qdn, 
 *                    solver_params_t * solver_params)  
 *
 **************************************************************************/


#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include "global.h"
#include "start.h"
#include "read_input.h"
#include "default_input_values.h"
#include "linalg/mul_gamma5.h"
#include "linalg/diff.h"
#include "linalg/square_norm.h"
#include "linalg/mul_r_gamma5.h"
#include "gamma.h"
// for the non-degenerate operator normalisation
#include "phmc.h"
#include "solver/solver.h"
#include "solver/solver_field.h"
#include "solver/init_guess.h"
#include "solver/matrix_mult_typedef.h"
#include "solver/solver_types.h"
#include "solver/solver_params.h"
#include "operator/tm_operators.h"
#include "operator/tm_operators_32.h"
#include "operator/tm_operators_nd.h"
#include "operator/tm_operators_nd_32.h"
#include "operator/clovertm_operators.h"
#include "operator/clovertm_operators_32.h"
#include "misc_types.h"
#include "monomial_solve.h"
#include "linalg_eo.h"
#ifdef DDalphaAMG
#include "DDalphaAMG_interface.h"
#endif
#ifdef TM_USE_QPHIX
#include "qphix_interface.h"
#endif

#include <io/params.h>
#include <io/spinor.h>

#ifdef HAVE_GPU
#include"../GPU/cudadefs.h"
extern  int linsolve_eo_gpu (spinor * const P, spinor * const Q, const int max_iter, 
                            double eps, const int rel_prec, const int N, matrix_mult f);
extern int dev_cg_mms_tm_nd(spinor ** const Pup, spinor ** const Pdn, 
		 spinor * const Qup, spinor * const Qdn, 
		 solver_params_t * solver_params);
   #ifdef TEMPORALGAUGE
     #include "../temporalgauge.h" 
   #endif
#include "read_input.h" 
#endif

int solve_degenerate(spinor * const P, spinor * const Q, solver_params_t solver_params,
                     const int max_iter, double eps_sq, const int rel_prec, 
                     const int N, matrix_mult f, int solver_type){
  int iteration_count = 0;

#ifdef TM_USE_QPHIX
  if(solver_params.external_inverter == QPHIX_INVERTER){
    spinor** temp;
    // temporary field required by the QPhiX solve
    init_solver_field(&temp, VOLUMEPLUSRAND/2, 1);
    // using CG for the HMC, we always want to have the solution of (Q Q^dagger) x = b, which is equivalent to
    // gamma_5 (M M^dagger)^{-1} gamma_5 b
    // FIXME: this needs to be adjusted to also support BICGSTAB
    gamma5(temp[0], Q, VOLUME/2);
    iteration_count = invert_eo_qphix_oneflavour(P, temp[0], max_iter, eps_sq, solver_type, 
                                                 rel_prec, solver_params, solver_params.sloppy_precision, solver_params.compression_type);
    mul_gamma5(P, VOLUME/2);
    if(g_debug_level > 0){
      f(temp[0], P);
      diff(temp[0], temp[0], Q, VOLUME/2);
      double diffnorm = square_norm(temp[0], VOLUME/2, 1); 
      if( g_proc_id == 0 ){
        printf("# solve_degenerate residual check: %e\n", diffnorm);
      }
    }
    finalize_solver(temp, 1);
  } else
#endif
  if(solver_type == MIXEDCG || solver_type == RGMIXEDCG){
    // the default mixed solver is rg_mixed_cg_her
    int (*msolver_fp)(spinor * const, spinor * const, solver_params_t, 
                      const int, double, const int, const int, matrix_mult, matrix_mult32) = rg_mixed_cg_her;

    // but it might be necessary at some point to use the old version
    if(solver_type == MIXEDCG){
      msolver_fp = mixed_cg_her;
    }

    // FIXME: this GPU stuff needs to go...
    if(usegpu_flag){   
      #ifdef HAVE_GPU     
        #ifdef TEMPORALGAUGE
          to_temporalgauge(g_gauge_field, Q , P);
        #endif          
        iteration_count = linsolve_eo_gpu(P, Q, max_iter, eps_sq, rel_prec, N, f);
        #ifdef TEMPORALGAUGE
          from_temporalgauge(Q, P);
        #endif
      #endif
      return(iteration_count);
    }
    else{
      if(f==Qtm_pm_psi){   
        iteration_count = msolver_fp(P, Q, solver_params, max_iter, eps_sq, rel_prec, N, f, &Qtm_pm_psi_32);
      } else if(f==Q_pm_psi){     
	iteration_count = msolver_fp(P, Q, solver_params, max_iter, eps_sq, rel_prec, N, f, &Q_pm_psi_32);
      } else if(f==Qsw_pm_psi){
        copy_32_sw_fields();
        iteration_count = msolver_fp(P, Q, solver_params, max_iter, eps_sq, rel_prec, N, f, &Qsw_pm_psi_32);
      } else {
        if(g_proc_id==0) printf("Warning: 32 bit matrix not available. Falling back to CG in 64 bit\n"); 
        solver_type = CG;
      }
    }
  } 
  else if(solver_type == CG){
    iteration_count =  cg_her(P, Q, max_iter, eps_sq, rel_prec, N, f);
  }
  else if(solver_type == BICGSTAB){
     iteration_count =  bicgstab_complex(P, Q, max_iter, eps_sq, rel_prec, N, f);     
  }
#ifdef DDalphaAMG 
  else if (solver_type == MG)
    iteration_count =  MG_solver(P, Q, eps_sq, max_iter,rel_prec, N , g_gauge_field, f);
#endif     
  else{
    if(g_proc_id==0) printf("Error: solver not allowed for degenerate solve. Aborting...\n");
    exit(2);
  }

  return(iteration_count);
}

int solve_mms_tm(spinor ** const P, spinor * const Q,
                 solver_params_t * solver_params){ 
  int iteration_count = 0; 

#ifdef TM_USE_QPHIX
  if( solver_params->external_inverter == QPHIX_INVERTER ){
    spinor ** temp;
    init_solver_field(&temp, VOLUMEPLUSRAND/2, 1);
    gamma5(temp[0], Q, VOLUME/2);
    iteration_count = invert_eo_qphix_oneflavour_mshift(P, temp[0],
                                                        solver_params->max_iter, solver_params->squared_solver_prec,
                                                        solver_params->type, solver_params->rel_prec,
                                                        *solver_params,
                                                        solver_params->sloppy_precision,
                                                        solver_params->compression_type);
    for( int shift = 0; shift < solver_params->no_shifts; shift++){
      mul_gamma5(P[shift], VOLUME/2);
    }
    if(g_debug_level > 0){
      // FIXME: in the shift-by-shift branch, the shifted operator exists explicitly and could be used to 
      // truly check the residual here
      solver_params->M_psi(temp[0], P[0]);
      diff(temp[0], temp[0], Q, VOLUME/2);
      double diffnorm = square_norm(temp[0], VOLUME/2, 1); 
      if( g_proc_id == 0 ){
        printf("# solve_mshift_oneflavour residual check: %e\n", diffnorm);
        printf("# NOTE that this currently repors the residual for the *unishfted* operator!\n");
      }
    }
    finalize_solver(temp, 1);
  } else
#endif // TM_USE_QPHIX
  if (solver_params->type == CGMMS){
    iteration_count = cg_mms_tm(P, Q, solver_params);
  }
#ifdef DDalphaAMG
  else if (solver_params->type == MG) {
    // if the mg_mms_mass is larger than the smallest shift we use MG
    if (mg_no_shifts > 0 || mg_mms_mass >= solver_params->shifts[0]) { 
      // if the mg_mms_mass is smaller than the larger shifts, we use CGMMS for those
      // in case mg_no_shifts is used, then mg_mms_mass = 0
      if(mg_mms_mass >= solver_params->shifts[0]) {
        mg_no_shifts = solver_params->no_shifts;
        while (mg_mms_mass < solver_params->shifts[mg_no_shifts-1]) { mg_no_shifts--; }
      }
      int no_shifts = solver_params->no_shifts;
      if (mg_no_shifts < no_shifts) {
        solver_params->no_shifts = no_shifts - mg_no_shifts;
        solver_params->shifts += mg_no_shifts;
        iteration_count = cg_mms_tm( P+mg_no_shifts, Q, solver_params );
        // Restoring solver_params
        solver_params->no_shifts = no_shifts;
        solver_params->shifts -= mg_no_shifts;
      }

      for(int i = mg_no_shifts-1; i>=0; i--){
        // preparing initial guess                                                                                                                                                                       
        init_guess_mms(P, Q, i, solver_params);
        g_mu3 = solver_params->shifts[i]; 
        iteration_count += MG_solver( P[i], Q, solver_params->squared_solver_prec, solver_params->max_iter,
                                         solver_params->rel_prec, solver_params->sdim, g_gauge_field, solver_params->M_psi );
        g_mu3 = _default_g_mu3;
      }
    } else {
      iteration_count = cg_mms_tm( P, Q, solver_params );
    }
  }
#endif
  else if (solver_params->type == RGMIXEDCG){
    matrix_mult32 f32  = Qtm_pm_psi_32;
    if( solver_params->M_psi == Qsw_pm_psi ){ 
      f32  = Qsw_pm_psi_32;
    }
    iteration_count = 0;
    // solver_params_t struct needs to be passed to all solvers except for cgmms, so we need to construct it here
    // and set the one relevant parameter
    solver_params_t temp_params;
    temp_params.mcg_delta = _default_mixcg_innereps;
    double iter_local = 0;
    for(int i = solver_params->no_shifts-1; i>=0; i--){
      // preparing initial guess
      init_guess_mms(P, Q, i, solver_params);
      
      // inverting
      g_mu3 = solver_params->shifts[i]; 
      iter_local = rg_mixed_cg_her( P[i], Q, temp_params, solver_params->max_iter,
                                    solver_params->squared_solver_prec, solver_params->rel_prec, solver_params->sdim,
                                    solver_params->M_psi, f32);
      g_mu3 = _default_g_mu3;
      if(iter_local == -1){
        return(-1);
      } else {
        iteration_count += iter_local;
      }
    }
  } else {
    if(g_proc_id==0) printf("Error: solver not allowed for TM mms solve. Aborting...\n");
    exit(2);      
  }

  return(iteration_count);
}

int solve_mms_nd(spinor ** const Pup, spinor ** const Pdn, 
                 spinor * const Qup, spinor * const Qdn, 
                 solver_params_t * solver_params){ 
  int iteration_count = 0; 

#ifdef TM_USE_QPHIX
  if(solver_params->external_inverter == QPHIX_INVERTER){
    spinor** temp;
    init_solver_field(&temp, VOLUMEPLUSRAND/2, 2);
    //  gamma5 (M.M^dagger)^{-1} gamma5 = [ Q(+mu,eps) Q(-mu,eps) ]^{-1}
    gamma5(temp[0], Qup, VOLUME/2);
    gamma5(temp[1], Qdn, VOLUME/2);
    iteration_count = invert_eo_qphix_twoflavour_mshift(Pup, Pdn, temp[0], temp[1],
                                                        solver_params->max_iter, solver_params->squared_solver_prec,
                                                        solver_params->type, solver_params->rel_prec,
                                                        *solver_params,
                                                        solver_params->sloppy_precision,
                                                        solver_params->compression_type);
    
    // the tmLQCD ND operator used for HMC is normalised by the inverse of the maximum eigenvalue
    // so the inverse of Q^2 is normalised by the square of the maximum eigenvalue
    // or, equivalently, the square of the inverse of the inverse
    // note that in the QPhiX interface, we also correctly normalise the shifts
    const double maxev_sq = (1.0/phmc_invmaxev)*(1.0/phmc_invmaxev);
    for( int shift = 0; shift < solver_params->no_shifts; shift++){
      mul_r_gamma5(Pup[shift], maxev_sq, VOLUME/2);
      mul_r_gamma5(Pdn[shift], maxev_sq, VOLUME/2);
    }
    if( g_debug_level > 0 ){
      // FIXME: in the shift-by-shift branch, the shifted operator exists explicitly and could be used to 
      // truly check the residual here
      solver_params->M_ndpsi(temp[0], temp[1], Pup[0], Pdn[0]);
      diff(temp[0], temp[0], Qup, VOLUME/2);
      diff(temp[1], temp[1], Qdn, VOLUME/2);
      double diffnorm = square_norm(temp[0], VOLUME/2, 1) + square_norm(temp[1], VOLUME/2, 1); 
      if( g_proc_id == 0 ){
        printf("# solve_mms_nd residual check: %e\n", diffnorm);
        printf("# NOTE that this currently repors the residual for the *unishfted* operator!\n");
      }
    }
    finalize_solver(temp, 2);
  } else
#endif //TM_USE_QPHIX
  if(solver_params->type==MIXEDCGMMSND){
    if(usegpu_flag){
    #ifdef HAVE_GPU      
      #ifdef TEMPORALGAUGE
      to_temporalgauge_mms(g_gauge_field , Qup, Qdn, Pup, Pdn, solver_params->no_shifts);
      #endif        
      iteration_count = dev_cg_mms_tm_nd(Pup, Pdn, Qup, Qdn, solver_params);  
      #ifdef TEMPORALGAUGE
      from_temporalgauge_mms(Qup, Qdn, Pup, Pdn, solver_params->no_shifts);
      #endif 
    #endif
    } else {
      iteration_count = mixed_cg_mms_tm_nd(Pup, Pdn, Qup, Qdn, solver_params);
    }
  } else if (solver_params->type == CGMMSND){
    iteration_count = cg_mms_tm_nd(Pup, Pdn, Qup, Qdn, solver_params);
    //Testing the initial guess
    if( g_debug_level > 3 ){
      for(int i = solver_params->no_shifts-1; i>=0; i--){
        // preparing initial guess
        init_guess_mms_nd(Pup, Pdn, Qup, Qdn, i, solver_params);
      }
    }
  }
#ifdef DDalphaAMG
  else if (solver_params->type == MG) {
    // if the mg_mms_mass is larger than the smallest shift we use MG
    if (mg_no_shifts > 0 || mg_mms_mass >= solver_params->shifts[0]) { 

      // if the mg_mms_mass is smaller than the larger shifts, we use CGMMS for those
      // in case mg_no_shifts is used, then mg_mms_mass = 0
      if(mg_mms_mass >= solver_params->shifts[0]) {
        mg_no_shifts = solver_params->no_shifts;
        while (mg_mms_mass < solver_params->shifts[mg_no_shifts-1]) { mg_no_shifts--; }
      }
      int no_shifts = solver_params->no_shifts;
      if (mg_no_shifts < no_shifts) {
        solver_params->no_shifts = no_shifts - mg_no_shifts;
        solver_params->shifts += mg_no_shifts;
        iteration_count = cg_mms_tm_nd( Pup+mg_no_shifts, Pdn+mg_no_shifts, Qup, Qdn, solver_params );
        // Restoring solver_params
        solver_params->no_shifts = no_shifts;
        solver_params->shifts -= mg_no_shifts;
      }

      matrix_mult_nd f = Qtm_pm_ndpsi_shift;
      if( solver_params->M_ndpsi == Qsw_pm_ndpsi ) 
        f = Qsw_pm_ndpsi_shift;

      for(int i = mg_no_shifts-1; i>=0; i--){
        // preparing initial guess
        init_guess_mms_nd(Pup, Pdn, Qup, Qdn, i, solver_params);
        
        g_shift = solver_params->shifts[i]*solver_params->shifts[i]; 
        iteration_count += MG_solver_nd( Pup[i], Pdn[i], Qup, Qdn, solver_params->squared_solver_prec, solver_params->max_iter,
                                         solver_params->rel_prec, solver_params->sdim, g_gauge_field, f );
        g_shift = _default_g_shift;
      }
    } else {
      iteration_count = cg_mms_tm_nd( Pup, Pdn, Qup, Qdn, solver_params );
    }
  }
#endif
  else if (solver_params->type == RGMIXEDCG){
    matrix_mult_nd   f    = Qtm_pm_ndpsi_shift;
    matrix_mult_nd32 f32  = Qtm_pm_ndpsi_shift_32;
    if( solver_params->M_ndpsi == Qsw_pm_ndpsi ){ 
      f    = Qsw_pm_ndpsi_shift;
      f32  = Qsw_pm_ndpsi_shift_32;
    }
    iteration_count = 0;
    // solver_params_t struct needs to be passed to all solvers except for cgmms, so we need to construct it here
    // and set the one relevant parameter
    solver_params_t temp_params;
    temp_params.mcg_delta = _default_mixcg_innereps;
    double iter_local = 0;
    for(int i = solver_params->no_shifts-1; i>=0; i--){
      // preparing initial guess
      init_guess_mms_nd(Pup, Pdn, Qup, Qdn, i, solver_params);
      
      // inverting
      g_shift = solver_params->shifts[i]*solver_params->shifts[i]; 
      iter_local = rg_mixed_cg_her_nd( Pup[i], Pdn[i], Qup, Qdn, temp_params, solver_params->max_iter,
                                       solver_params->squared_solver_prec, solver_params->rel_prec, solver_params->sdim, f, f32);
      g_shift = _default_g_shift;
      if(iter_local == -1){
        return(-1);
      } else {
        iteration_count += iter_local;
      }
    }
  } else {
    if(g_proc_id==0) printf("Error: solver not allowed for ND mms solve. Aborting...\n");
    exit(2);      
  }

  return(iteration_count);
}

int solve_mms_nd_plus(spinor ** const Pup, spinor ** const Pdn, 
                      spinor * const Qup, spinor * const Qdn, 
                      solver_params_t * solver_params){ 

  int iteration_count = 0; 

#ifdef DDalphaAMG
  // With MG we can solve directly the unsquared operator
  if( solver_params->type == MG ){
    matrix_mult_nd f = Qtm_tau1_ndpsi_add_Ishift;
    if( solver_params->M_ndpsi == Qsw_pm_ndpsi )
      f = Qsw_tau1_ndpsi_add_Ishift;
    for(int i = solver_params->no_shifts-1; i>=0; i--){
      // preparing initial guess
      init_guess_mms_nd_plus(Pup, Pdn, Qup, Qdn, i, solver_params);
  
      // g_shift = shift^2 and then in Qsw_tau1_ndpsi_add_Ishift the square root is taken
      g_shift = solver_params->shifts[i]*solver_params->shifts[i]; 
      iteration_count += MG_solver_nd( Pup[i], Pdn[i], Qup, Qdn, solver_params->squared_solver_prec,
                                       solver_params->max_iter, solver_params->rel_prec, solver_params->sdim,
                                       g_gauge_field, f );
      g_shift = _default_g_shift;
    }
  } else 
#endif
  {
    iteration_count = solve_mms_nd(Pup, Pdn, Qup, Qdn, solver_params);
    
    // apply operator for retrieving unsquared solution
    matrix_mult_nd f = Qtm_tau1_ndpsi_sub_Ishift;
    if( solver_params->M_ndpsi == Qsw_pm_ndpsi )
      f = Qsw_tau1_ndpsi_sub_Ishift;
    spinor** temp;
    init_solver_field(&temp, VOLUMEPLUSRAND/2, 1);
    for(int i = solver_params->no_shifts-1; i>=0; i--){
      g_shift = solver_params->shifts[i]*solver_params->shifts[i]; 
      f(temp[0],temp[1],Pup[i],Pdn[i]);
      assign(Pup[i], temp[0], VOLUME/2);
      assign(Pdn[i], temp[1], VOLUME/2);
      g_shift = _default_g_shift;
    }
    finalize_solver(temp, 2);
  }
  return iteration_count;
}

/***********************************************************************
 *
 * Copyright (C) 2014 Florian Burger
 *               2017,2018,2019,2020,2021 Bartosz Kostrzewa
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
# include<tmlqcd_config.h>
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
#include "gettime.h"
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
#include "fatal_error.h"

#include <io/params.h>
#include <io/spinor.h>

#ifdef TM_USE_QUDA
#  include "quda_interface.h"
#endif

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
  tm_stopwatch_push(&g_timers);
  int iteration_count = 0;

  // temporary field required by the QPhiX solve or by residual check
  int nr_sf = 1;
  if(g_debug_level > 3) nr_sf = 2;
  spinor** temp;
  if(g_debug_level > 2 || g_strict_residual_check || 
     solver_params.external_inverter == QPHIX_INVERTER || solver_params.external_inverter == QUDA_INVERTER ){
    init_solver_field(&temp, VOLUMEPLUSRAND/2, nr_sf);
  }
  
  solver_params.use_initial_guess = 0;


#ifdef TM_USE_QUDA
  if ( solver_params.external_inverter == QUDA_INVERTER){

    int QmQp = (f == Qsw_pm_psi || f == Qtm_pm_psi);

    tm_stopwatch_push(&g_timers);
    gamma5(temp[0], Q, VOLUME/2);
    tm_stopwatch_pop(&g_timers, 0, 1, current_mnl, "solve_degenerate:gamma5");
    iteration_count = invert_eo_MMd_quda(P,   //spinor * const Odd_new,
                                         temp[0],
                                         eps_sq, // Marco: check this:   const double precision, 
                                         max_iter,
                                         solver_type,  rel_prec,
                                         1, // Marco: 0 or 1 ? int even_odd_flag,
                                         solver_params,
                                         solver_params.sloppy_precision,
                                         solver_params.compression_type,
                                         QmQp);
    
    if( !(solver_type == MG || solver_type == BICGSTAB) ){
      tm_stopwatch_push(&g_timers);
      mul_gamma5(P, VOLUME/2);
      tm_stopwatch_pop(&g_timers, 0, 1, current_mnl, "solve_degenerate:gamma5");
    }

    //// //////////////////////////////////////////////////////////////// test to be removed
    //// // try matrix application directly
    //// spinor** tempE;
    //// init_solver_field(&tempE, VOLUMEPLUSRAND/2, 2);
    //// //point like source only if mpi=1
    //// for(int x =0; x < (VOLUMEPLUSRAND/2);x++){
    ////   // tempE[0][x].s0.c0=0.0; 
    ////   // tempE[0][x].s0.c1=0.0; 
    ////   // tempE[0][x].s0.c2=0.0; 
    ////   // tempE[0][x].s1.c0=0.0; 
    ////   // tempE[0][x].s1.c1=0.0; 
    ////   // tempE[0][x].s1.c2=0.0; 
    ////   // tempE[0][x].s2.c0=0.0; 
    ////   // tempE[0][x].s2.c1=0.0; 
    ////   // tempE[0][x].s2.c2=0.0; 
    ////   // tempE[0][x].s3.c0=0.0; 
    ////   // tempE[0][x].s3.c1=0.0; 
    ////   // tempE[0][x].s3.c2=0.0; 
    ////   
    ////   // random
    ////   tempE[0][x].s0.c0=((double)rand())/RAND_MAX;
    ////   tempE[0][x].s0.c1=((double)rand())/RAND_MAX;
    ////   tempE[0][x].s0.c2=((double)rand())/RAND_MAX;
    ////   tempE[0][x].s1.c0=((double)rand())/RAND_MAX;
    ////   tempE[0][x].s1.c1=((double)rand())/RAND_MAX;
    ////   tempE[0][x].s1.c2=((double)rand())/RAND_MAX;
    ////   tempE[0][x].s2.c0=((double)rand())/RAND_MAX;
    ////   tempE[0][x].s2.c1=((double)rand())/RAND_MAX;
    ////   tempE[0][x].s2.c2=((double)rand())/RAND_MAX;
    ////   tempE[0][x].s3.c0=((double)rand())/RAND_MAX;
    ////   tempE[0][x].s3.c1=((double)rand())/RAND_MAX;
    ////   tempE[0][x].s3.c2=((double)rand())/RAND_MAX;
    //// }
    //// // set something other than component (0,0) to 1.0
    //// // tempE[0][0].s0.c1=1.0;
    //// // tempE[0][0].s0.c1=1.0;
    //// // tempE[0][0].s0.c2=1.0;
    //// // tempE[0][0].s1.c1=1.0;
    //// // tempE[0][0].s2.c1=1.0;
    //// // tempE[0][0].s3.c1=1.0;

    //// // just in case: copy the source
    //// assign(tempE[1], tempE[0], VOLUMEPLUSRAND/2);

    //// M_quda(P, tempE[0]); // quda changes the source
    //// if( f == Qtm_pm_psi ){
    ////   Mtm_plus_psi(temp[0], tempE[1]);
    //// } else if ( f == Qsw_pm_psi ){
    ////   Msw_plus_psi(temp[0], tempE[1]);
    //// }

    //// // almost certainly we need to account for the gamma basis
    //// for (int ix=0; ix < (VOLUME/2); ix++){
    ////   spinor *hp=((spinor*)temp[0]) + ix;
    ////   spinor *dp=((spinor*)P) + ix;
    ////   double r=creal((hp)->s0.c0)-creal((dp)->s0.c0);
    ////   printf("ix=%d, r=%.3e\n"
    ////          "re tmLQCD=(%.3e,%.3e,%.3e), (%.3e,%.3e,%.3e), (%.3e,%.3e,%.3e), (%.3e,%.3e,%.3e)\n"
    ////          "re  quda=(%.3e,%.3e,%.3e), (%.3e,%.3e,%.3e), (%.3e,%.3e,%.3e), (%.3e,%.3e,%.3e)\n"
    ////          "im tmLQCD=(%.3e,%.3e,%.3e), (%.3e,%.3e,%.3e), (%.3e,%.3e,%.3e), (%.3e,%.3e,%.3e)\n"
    ////          "im  quda=(%.3e,%.3e,%.3e), (%.3e,%.3e,%.3e), (%.3e,%.3e,%.3e), (%.3e,%.3e,%.3e)\n",
    ////          ix,r,
    ////          creal((hp)->s0.c0), creal((hp)->s0.c1), creal((hp)->s0.c2),
    ////          creal((hp)->s1.c0), creal((hp)->s1.c1), creal((hp)->s1.c2),
    ////          creal((hp)->s2.c0), creal((hp)->s2.c1), creal((hp)->s2.c2),
    ////          creal((hp)->s3.c0), creal((hp)->s3.c1), creal((hp)->s3.c2),

    ////          creal((dp)->s0.c0), creal((dp)->s0.c1), creal((dp)->s0.c2),
    ////          creal((dp)->s1.c0), creal((dp)->s1.c1), creal((dp)->s1.c2),
    ////          creal((dp)->s2.c0), creal((dp)->s2.c1), creal((dp)->s2.c2),
    ////          creal((dp)->s3.c0), creal((dp)->s3.c1), creal((dp)->s3.c2),
    ////          
    ////          cimag((hp)->s0.c0), cimag((hp)->s0.c1), cimag((hp)->s0.c2),
    ////          cimag((hp)->s1.c0), cimag((hp)->s1.c1), cimag((hp)->s1.c2),
    ////          cimag((hp)->s2.c0), cimag((hp)->s2.c1), cimag((hp)->s2.c2),
    ////          cimag((hp)->s3.c0), cimag((hp)->s3.c1), cimag((hp)->s3.c2),

    ////          cimag((dp)->s0.c0), cimag((dp)->s0.c1), cimag((dp)->s0.c2),
    ////          cimag((dp)->s1.c0), cimag((dp)->s1.c1), cimag((dp)->s1.c2),
    ////          cimag((dp)->s2.c0), cimag((dp)->s2.c1), cimag((dp)->s2.c2),
    ////          cimag((dp)->s3.c0), cimag((dp)->s3.c1), cimag((dp)->s3.c2)
    ////         );
    //// }
    //// printf("\n\n\n");
    //// print_spinor_similar_components(temp[0], P, VOLUME/2, 1e-4);
    //// finalize_solver(tempE,2);
    //// exit(1);
    //// //////////////////////////////////////////////////////////// end of the test to be removed

  } else
#endif
#ifdef TM_USE_QPHIX
  if(solver_params.external_inverter == QPHIX_INVERTER){
    // using CG for the HMC, we always want to have the solution of (Q Q^dagger) x = b, which is equivalent to
    // gamma_5 (M M^dagger)^{-1} gamma_5 b
    // FIXME: this needs to be adjusted to also support BICGSTAB
    gamma5(temp[0], Q, VOLUME/2);
    iteration_count = invert_eo_qphix_oneflavour(P, temp[0], max_iter, eps_sq, solver_type, 
                                                 rel_prec, solver_params, solver_params.sloppy_precision, solver_params.compression_type);
    mul_gamma5(P, VOLUME/2);
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
    fatal_error("Error: solver not allowed for degenerate solve. Aborting...\n", "solve_degenerate");
  }

  if(g_debug_level > 2 || g_strict_residual_check){
    if(g_proc_id == 0) printf("# solve_degenerate applying operator\n");
    f(temp[0], P);
    if(g_debug_level > 5){
      ratio(temp[1], temp[0], Q, VOLUME/2);
      if(g_proc_id == 0){
        //print_spinor_similar_components(temp[0], Q, VOLUME/2, 1.0e-8);
        print_spinor(temp[1], VOLUME/2);
      }
    }
    diff(temp[0], temp[0], Q, VOLUME/2);
    double diffnorm = square_norm(temp[0], VOLUME/2, 1); 
    double resnorm  = square_norm(P, VOLUME/2, 1); 
    double rel_nrm = rel_prec ? square_norm(Q, VOLUME/2, 1) : 1.0;
    if( g_proc_id == 0 ){
      // checking the norm of the result to make sure it's not zero
      printf("# solve_degenerate result norm: %e\n", resnorm);
      printf("# solve_degenerate target residual: %e\n", eps_sq * rel_nrm);
      printf("# solve_degenerate residual check norm: %e\n", diffnorm);
      fflush(stdout);
    }
    if( g_strict_residual_check ){
      if( diffnorm > 1.1*(eps_sq * rel_nrm) ){
        fatal_error("Residual norm exceeds target by more than a factor of 1.1!", "solve_degenerate");
      }
    }
  }
  if(g_debug_level > 2 || g_strict_residual_check ||
      solver_params.external_inverter == QPHIX_INVERTER  || solver_params.external_inverter == QUDA_INVERTER){
    finalize_solver(temp, nr_sf);
  }
  
  tm_stopwatch_pop(&g_timers, 0, 1, current_mnl, __func__);
  return(iteration_count);
}

int solve_mms_tm(spinor ** const P, spinor * const Q,
                 solver_params_t * solver_params){
  tm_stopwatch_push(&g_timers); 
  int iteration_count = 0; 

  solver_params->use_initial_guess = 0;

  // temporary field required by the QPhiX solve or by residual check
  spinor ** temp;
  if(g_debug_level > 2 || g_strict_residual_check ||
      (solver_params->external_inverter == QUDA_INVERTER ) ||
      (solver_params->external_inverter == QPHIX_INVERTER  && solver_params->type != MG)){
    init_solver_field(&temp, VOLUMEPLUSRAND/2, 1);
  }

#ifdef TM_USE_QPHIX
  if( solver_params->external_inverter == QPHIX_INVERTER && solver_params->type != MG ){
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
  } else
#endif

#ifdef TM_USE_QUDA
  if ( solver_params->external_inverter == QUDA_INVERTER && solver_params->type != MG){
    tm_stopwatch_push(&g_timers);
    gamma5(temp[0], Q, VOLUME/2);
    tm_stopwatch_pop(&g_timers, 0, 1, current_mnl, "solve_mms_tm:gamma5");
    iteration_count = invert_eo_quda_oneflavour_mshift(P, temp[0],
                                                       solver_params->squared_solver_prec,
                                                       solver_params->max_iter,
                                                       solver_params->type,
                                                       solver_params->rel_prec,
                                                       1,
                                                       *solver_params,
                                                       solver_params->sloppy_precision,
                                                       solver_params->compression_type);
    tm_stopwatch_push(&g_timers);
    for( int shift = 0; shift < solver_params->no_shifts; shift++){
      mul_gamma5(P[shift], VOLUME/2);
    }
    tm_stopwatch_pop(&g_timers, 0, 1, current_mnl, "solve_mms_tm:gamma5");
  } else
#endif
    if (solver_params->type == CGMMS){
      iteration_count = cg_mms_tm(P, Q, solver_params);
    }
#ifdef DDalphaAMG
  else if (solver_params->type == MG) {
    // if the mg_mms_mass is larger than the smallest shift we use MG
    if (mg_no_shifts > 0 || mg_mms_mass >= solver_params->shifts[0]) { 
      int nshifts = solver_params->no_shifts;
      int mg_nshifts = mg_no_shifts > nshifts ? nshifts:mg_no_shifts;
      // if the mg_mms_mass is smaller than the larger shifts, we use CGMMS for those
      // in case mg_no_shifts is used, then mg_mms_mass = 0
      if(mg_mms_mass >= solver_params->shifts[0]) {
        mg_nshifts = solver_params->no_shifts;
        while (mg_mms_mass < solver_params->shifts[mg_nshifts-1]) { mg_nshifts--; }
      }
      // Number of initial guesses provided by gcmms
      // README: tunable value. 1 it's fine for now.
      int  no_cgmms_init_guess = 1;
      if(no_cgmms_init_guess > mg_nshifts) {
        no_cgmms_init_guess = mg_nshifts;
      }
#ifdef TM_USE_QPHIX
      if( solver_params->external_inverter == QPHIX_INVERTER && mg_nshifts < nshifts ) {
        // TODO: no initial guess option with QphiX
        no_cgmms_init_guess = 0;
        spinor ** P_cg = P+(mg_nshifts - no_cgmms_init_guess);
        double * shifts_start = solver_params->shifts;
        solver_params->no_shifts = nshifts - (mg_nshifts - no_cgmms_init_guess);
        solver_params->shifts += (mg_nshifts - no_cgmms_init_guess);
        solver_params->type = CGMMS;
        gamma5(temp[0], Q, VOLUME/2);
        iteration_count = invert_eo_qphix_oneflavour_mshift(P, temp[0],
                                                            solver_params->max_iter, solver_params->squared_solver_prec,
                                                            solver_params->type, solver_params->rel_prec,
                                                            *solver_params,
                                                            solver_params->sloppy_precision,
                                                            solver_params->compression_type);
        for( int shift = 0; shift < solver_params->no_shifts; shift++) {
          mul_gamma5(P[shift], VOLUME/2);
        }
        // Restoring solver_params
        solver_params->no_shifts = nshifts;
        solver_params->shifts = shifts_start;
        solver_params->type = MG;
        } else
#endif // TM_USE_QPHIX  
      if (mg_nshifts < nshifts) {
        spinor ** P_cg = P+(mg_nshifts - no_cgmms_init_guess);
        double * shifts_start = solver_params->shifts;
        solver_params->no_shifts = nshifts - (mg_nshifts - no_cgmms_init_guess);
        solver_params->shifts += (mg_nshifts - no_cgmms_init_guess);
        solver_params->type = CGMMS;
        // switching last shift. We run CGMMS for the shift we want to solve.
        if (no_cgmms_init_guess > 0) {
          SWAP(solver_params->shifts[0], solver_params->shifts[no_cgmms_init_guess]);
          SWAP(P_cg[0], P_cg[no_cgmms_init_guess]);
        }
        iteration_count = solve_mms_tm( P_cg, Q, solver_params );
        // Switching back last shift
        if (no_cgmms_init_guess > 0) {
          SWAP(solver_params->shifts[0], solver_params->shifts[no_cgmms_init_guess]);
          SWAP(P_cg[0], P_cg[no_cgmms_init_guess]);
        }
        // Restoring solver_params
        solver_params->no_shifts = nshifts;
        solver_params->shifts = shifts_start;
        solver_params->type = MG;
      } else {
        no_cgmms_init_guess = 0;
      }

      for(int i = mg_nshifts-1; i>=0; i--){
        // preparing initial guess
        if(i<mg_nshifts-no_cgmms_init_guess)
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
      solver_params->use_initial_guess = 1;
     
      // inverting
      g_mu3 = solver_params->shifts[i]; 
      iter_local = rg_mixed_cg_her( P[i], Q, temp_params, solver_params->max_iter,
                                    solver_params->squared_solver_prec, solver_params->rel_prec, solver_params->sdim,
                                    solver_params->M_psi, f32);
      g_mu3 = _default_g_mu3;
      solver_params->use_initial_guess = 0;
      if(iter_local == -1){
        return(-1);
      } else {
        iteration_count += iter_local;
      }
    }
  } else {
    fatal_error("Error: solver not allowed for TM mms solve. Aborting...\n", "solve_mms_tm");
  }

  if(g_debug_level > 2 || g_strict_residual_check){
    matrix_mult f = Qtm_pm_psi_shift;
    if( solver_params->M_psi == Qsw_pm_psi ) 
      f = Qsw_pm_psi_shift;
    
    int rel_nrm = solver_params->rel_prec && g_strict_residual_check ? square_norm(Q, VOLUME/2, 1) : 1.0;
    int check_fail = 0;
    for( int shift = 0; shift < solver_params->no_shifts; shift++){
      g_shift = solver_params->shifts[shift]*solver_params->shifts[shift];
      f(temp[0], P[shift]);
      diff(temp[0], temp[0], Q, VOLUME/2);
      double diffnorm = square_norm(temp[0], VOLUME/2, 1); 
      if( g_proc_id == 0 ){
        printf("# solve_mms_tm residual check: shift %d (%.6e), res. %e\n", shift, g_shift, diffnorm);
        fflush(stdout);
      }
      g_shift = _default_g_shift;

      if( g_strict_residual_check ){
        check_fail += diffnorm > 1.1*( solver_params->squared_solver_prec * rel_nrm );
      }
    }
    if( g_strict_residual_check && check_fail > 0 ){
      fatal_error("Residual norm for at least one shift exceeds target by more than a factor of 1.1!", "solve_mms_tm");
    }
  }
  if(g_debug_level > 2 || g_strict_residual_check ||
      (solver_params->external_inverter == QUDA_INVERTER) ||
      (solver_params->external_inverter == QPHIX_INVERTER && solver_params->type != MG)){
    finalize_solver(temp, 1);
  }
  tm_stopwatch_pop(&g_timers, 0, 1, current_mnl, __func__); 
  return(iteration_count);
}

int solve_mms_nd(spinor ** const Pup, spinor ** const Pdn, 
                 spinor * const Qup, spinor * const Qdn, 
                 solver_params_t * solver_params){
  tm_stopwatch_push(&g_timers); 
  int iteration_count = 0; 
  solver_params->use_initial_guess = 0;

  // temporary field required by the QPhiX solve or by residual check
  spinor ** temp;
  if(g_debug_level > 2 || g_strict_residual_check ||
     (solver_params->external_inverter == QUDA_INVERTER) ||
     (solver_params->external_inverter == QPHIX_INVERTER && solver_params->type != MG) ){
    init_solver_field(&temp, VOLUMEPLUSRAND/2, 2);
  }

#ifdef TM_USE_QUDA
  if(solver_params->external_inverter == QUDA_INVERTER ){
    //  gamma5 (M.M^dagger)^{-1} gamma5 = [ Q(+mu,eps) Q(-mu,eps) ]^{-1}
    tm_stopwatch_push(&g_timers);
    gamma5(temp[0], Qup, VOLUME/2);
    gamma5(temp[1], Qdn, VOLUME/2);
    tm_stopwatch_pop(&g_timers, 0, 1, current_mnl, "solve_mms_nd:gamma5");
    iteration_count = invert_eo_quda_twoflavour_mshift(Pup, Pdn, temp[0], temp[1],
                                                       solver_params->squared_solver_prec, solver_params->max_iter,
                                                       solver_params->type, solver_params->rel_prec,
                                                       1 /*even-odd-flag*/, *solver_params,
                                                       solver_params->sloppy_precision,
                                                       solver_params->compression_type);
    
    // the tmLQCD ND operator used for HMC is normalised by the inverse of the maximum eigenvalue
    // so the inverse of QpQm is normalised by the square of the maximum eigenvalue
    // or, equivalently, the square of the inverse of the inverse
    // note that in the QUDA interface, we also correctly normalise the shifts
    const double maxev_sq = (1.0/phmc_invmaxev)*(1.0/phmc_invmaxev);
    tm_stopwatch_push(&g_timers);
    for( int shift = 0; shift < solver_params->no_shifts; shift++){
      mul_r_gamma5(Pup[shift], maxev_sq, VOLUME/2);
      mul_r_gamma5(Pdn[shift], maxev_sq, VOLUME/2);
    }
    tm_stopwatch_pop(&g_timers, 0, 1, current_mnl, "solve_mms_nd:mul_r_gamm5");
  } else
#endif //TM_USE_QPHIX
#ifdef TM_USE_QPHIX
  if(solver_params->external_inverter == QPHIX_INVERTER && solver_params->type != MG){
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
  }
#ifdef DDalphaAMG
  else if (solver_params->type == MG) {
    // if the mg_mms_mass is larger than the smallest shift we use MG
    if (mg_no_shifts > 0 || mg_mms_mass >= solver_params->shifts[0]) { 

      int nshifts = solver_params->no_shifts;
      int mg_nshifts = mg_no_shifts > nshifts ? nshifts:mg_no_shifts;
      // if the mg_mms_mass is smaller than the larger shifts, we use CGMMS for those
      // in case mg_no_shifts is used, then mg_mms_mass = 0
      if(mg_mms_mass >= solver_params->shifts[0]) {
        mg_nshifts = nshifts;
        while (mg_mms_mass < solver_params->shifts[mg_nshifts-1]) { mg_nshifts--; }
      }
      // Number of initial guesses provided by gcmms
      // README: tunable value. 2 it's fine for now.
      int no_cgmms_init_guess = 2;
      if(no_cgmms_init_guess > mg_nshifts) {
        no_cgmms_init_guess = mg_nshifts;
      }
#ifdef TM_USE_QPHIX
      if(solver_params->external_inverter == QPHIX_INVERTER && mg_nshifts < nshifts){
        // TODO: no initial guess option with QphiX
        no_cgmms_init_guess = 0;
        spinor ** Pup_cg = Pup+(mg_nshifts - no_cgmms_init_guess);
        spinor ** Pdn_cg = Pdn+(mg_nshifts - no_cgmms_init_guess);
        double * shifts_start = solver_params->shifts;
        solver_params->no_shifts = nshifts - (mg_nshifts - no_cgmms_init_guess);
        solver_params->shifts += (mg_nshifts - no_cgmms_init_guess);
        solver_params-> type = CGMMSND;
        //  gamma5 (M.M^dagger)^{-1} gamma5 = [ Q(+mu,eps) Q(-mu,eps) ]^{-1}
        gamma5(temp[0], Qup, VOLUME/2);
        gamma5(temp[1], Qdn, VOLUME/2);
        iteration_count = invert_eo_qphix_twoflavour_mshift(Pup_cg, Pdn_cg, temp[0], temp[1],
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
        // Restoring solver_params
        solver_params->no_shifts = nshifts;
        solver_params->shifts = shifts_start;
        solver_params-> type = MG;
      } else
#endif //TM_USE_QPHIX
      if (mg_nshifts < nshifts) {
        spinor ** Pup_cg = Pup+(mg_nshifts - no_cgmms_init_guess);
        spinor ** Pdn_cg = Pdn+(mg_nshifts - no_cgmms_init_guess);
        double * shifts_start = solver_params->shifts;
        solver_params->no_shifts = nshifts - (mg_nshifts - no_cgmms_init_guess);
        solver_params->shifts += (mg_nshifts - no_cgmms_init_guess);
        solver_params-> type = CGMMSND;
        if (no_cgmms_init_guess > 0) {
          SWAP(solver_params->shifts[0], solver_params->shifts[no_cgmms_init_guess]);
          SWAP(Pup_cg[0], Pup_cg[no_cgmms_init_guess]);
          SWAP(Pdn_cg[0], Pdn_cg[no_cgmms_init_guess]);
        }
        iteration_count = solve_mms_nd( Pup_cg, Pdn_cg, Qup, Qdn, solver_params );
        // Switching back last shift
        if (no_cgmms_init_guess > 0) {
          SWAP(solver_params->shifts[0], solver_params->shifts[no_cgmms_init_guess]);
          SWAP(Pup_cg[0], Pup_cg[no_cgmms_init_guess]);
          SWAP(Pdn_cg[0], Pdn_cg[no_cgmms_init_guess]);
        }
        // Restoring solver_params
        solver_params->no_shifts = nshifts;
        solver_params->shifts = shifts_start;
        solver_params-> type = MG;
      } else {
        no_cgmms_init_guess = 0;
      }

      matrix_mult_nd f = Qtm_pm_ndpsi_shift;
      if( solver_params->M_ndpsi == Qsw_pm_ndpsi ) 
        f = Qsw_pm_ndpsi_shift;

      for(int i = mg_nshifts-1; i>=0; i--){
        // preparing initial guess
        if(i<mg_nshifts-no_cgmms_init_guess)
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
      solver_params->use_initial_guess = 1;
      
      // inverting
      g_shift = solver_params->shifts[i]*solver_params->shifts[i]; 
      iter_local = rg_mixed_cg_her_nd( Pup[i], Pdn[i], Qup, Qdn, temp_params, solver_params->max_iter,
                                       solver_params->squared_solver_prec, solver_params->rel_prec, solver_params->sdim, f, f32);
      g_shift = _default_g_shift;
      solver_params->use_initial_guess = 0;
      if(iter_local == -1){
        return(-1);
      } else {
        iteration_count += iter_local;
      }
    }
  } else {
    fatal_error("Error: solver not allowed for ND mms solve. Aborting...\n", "solve_mss_nd");
  }
  
  if(g_debug_level > 2 || g_strict_residual_check){
    matrix_mult_nd f = Qtm_pm_ndpsi_shift;
    if( solver_params->M_ndpsi == Qsw_pm_ndpsi ) 
      f = Qsw_pm_ndpsi_shift;
    
    int rel_nrm = solver_params->rel_prec && g_strict_residual_check ? 
                    square_norm(Qup, VOLUME/2, 1) + square_norm(Qdn, VOLUME/2, 1) : 1.0;
    int check_fail = 0;
    for( int shift = 0; shift < solver_params->no_shifts; shift++){
      g_shift = solver_params->shifts[shift]*solver_params->shifts[shift];
      f(temp[0], temp[1], Pup[shift], Pdn[shift]);
      diff(temp[0], temp[0], Qup, VOLUME/2);
      diff(temp[1], temp[1], Qdn, VOLUME/2);
      double diffnorm = square_norm(temp[0], VOLUME/2, 1) + square_norm(temp[1], VOLUME/2, 1);
      if( g_proc_id == 0 ){
        printf("# solve_mms_nd residual check: shift %d (%.6e), res. %e\n", shift, g_shift, diffnorm);
        fflush(stdout);
      }
      g_shift = _default_g_shift;

      if( g_strict_residual_check ){
        check_fail += diffnorm > 1.1*( solver_params->squared_solver_prec * rel_nrm );
      }
    }
    if( g_strict_residual_check && check_fail > 0 ){
      fatal_error("Residual norm for at least one shift exceeds target by more than a factor of 1.1!", "solve_mms_nd");
    }
  }
  if(g_debug_level > 2 || g_strict_residual_check ||
      (solver_params->external_inverter == QUDA_INVERTER) ||
      (solver_params->external_inverter == QPHIX_INVERTER && solver_params->type != MG)){
    finalize_solver(temp, 2);
  }
  tm_stopwatch_pop(&g_timers, 0, 1, current_mnl, __func__);
  return(iteration_count);
}

int solve_mms_nd_plus(spinor ** const Pup, spinor ** const Pdn, 
                      spinor * const Qup, spinor * const Qdn, 
                      solver_params_t * solver_params){
  tm_stopwatch_push(&g_timers); 

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
    init_solver_field(&temp, VOLUMEPLUSRAND/2, 2);
    tm_stopwatch_push(&g_timers);
    for(int i = solver_params->no_shifts-1; i>=0; i--){
      g_shift = solver_params->shifts[i]*solver_params->shifts[i]; 
      f(temp[0],temp[1],Pup[i],Pdn[i]);
      assign(Pup[i], temp[0], VOLUME/2);
      assign(Pdn[i], temp[1], VOLUME/2);
      g_shift = _default_g_shift;
    }
    tm_stopwatch_pop(&g_timers, 0, 1, current_mnl, "solve_mms_nd_plus:f_assign");
    finalize_solver(temp, 2);
  }
  tm_stopwatch_pop(&g_timers, 0, 1, current_mnl, __func__);
  return iteration_count;
}

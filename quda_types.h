/***********************************************************************
 * Copyright (C) 2017 Bartosz Kostrzewa
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
 *******************************************************************************/

#ifndef TM_QUDA_TYPES_H
#define TM_QUDA_TYPES_H

#ifdef HAVE_CONFIG_H
#  include<tmlqcd_config.h>
#endif

#ifdef TM_USE_QUDA
#include <quda.h>
#else
#include "quda_dummy_types.h"
#endif

#include "global.h"
#include <float.h>
#include <math.h>

typedef enum tm_quda_ferm_bc_t {
  TM_QUDA_THETABC = 0,
  TM_QUDA_APBC,
  TM_QUDA_PBC
} tm_quda_ferm_bc_t;

/* tm_QudaParams_t provides an interface between the tmLQCD input file and the
 * available QUDA parameters. At the moment, only the fermionic bounday conditions
 * and the MG parameters are exposed like this, but a further refactoring might
 * turn this into a complete representation of the possible input parameters */
typedef struct tm_QudaParams_t {
  int                  enable_device_memory_pool;
  int                  enable_pinned_memory_pool;

  tm_quda_ferm_bc_t fermionbc;

  int                  pipeline;
  double               reliable_delta;
  int                  gcrNkrylov;

  int                  mg_n_level;
  QudaVerbosity        mg_verbosity[QUDA_MAX_MG_LEVEL];                
  int                  mg_n_vec[QUDA_MAX_MG_LEVEL];
  int                  mg_blocksize[QUDA_MAX_MG_LEVEL][4];
  double               mg_mu_factor[QUDA_MAX_MG_LEVEL];
  QudaInverterType     mg_setup_inv_type;
  double               mg_setup_2kappamu;
  double               mg_setup_tol[QUDA_MAX_MG_LEVEL];
  int                  mg_setup_maxiter[QUDA_MAX_MG_LEVEL];
  QudaInverterType     mg_coarse_solver_type[QUDA_MAX_MG_LEVEL];
  int                  mg_coarse_solver_maxiter[QUDA_MAX_MG_LEVEL];
  double               mg_coarse_solver_tol[QUDA_MAX_MG_LEVEL];
  int                  mg_nu_pre[QUDA_MAX_MG_LEVEL];
  int                  mg_nu_post[QUDA_MAX_MG_LEVEL];
  QudaInverterType     mg_smoother_type[QUDA_MAX_MG_LEVEL];
  double               mg_smoother_tol[QUDA_MAX_MG_LEVEL];
  double               mg_omega[QUDA_MAX_MG_LEVEL];
  int                  mg_run_verify;
  int                  mg_enable_size_three_blocks;
  double               mg_reuse_setup_mu_threshold;
  double               mg_reset_setup_mdu_threshold;
  double               mg_refresh_setup_mdu_threshold;

  int                  mg_setup_maxiter_refresh[QUDA_MAX_MG_LEVEL];
  
  // parameters related to communication-avoiding
  // solvers  
  QudaCABasis          mg_setup_ca_basis[QUDA_MAX_MG_LEVEL];
  int                  mg_setup_ca_basis_size[QUDA_MAX_MG_LEVEL];
  double               mg_setup_ca_lambda_min[QUDA_MAX_MG_LEVEL];
  double               mg_setup_ca_lambda_max[QUDA_MAX_MG_LEVEL];

  QudaCABasis          mg_coarse_solver_ca_basis[QUDA_MAX_MG_LEVEL];
  int                  mg_coarse_solver_ca_basis_size[QUDA_MAX_MG_LEVEL];
  double               mg_coarse_solver_ca_lambda_min[QUDA_MAX_MG_LEVEL];
  double               mg_coarse_solver_ca_lambda_max[QUDA_MAX_MG_LEVEL];

  // parameters related to coarse grid deflation in the MG
  int                  mg_use_eig_solver[QUDA_MAX_MG_LEVEL];
  int                  mg_eig_preserve_deflation;
  int                  mg_eig_nEv[QUDA_MAX_MG_LEVEL];
  int                  mg_eig_nKr[QUDA_MAX_MG_LEVEL];
  int                  mg_eig_require_convergence[QUDA_MAX_MG_LEVEL];
  int                  mg_eig_check_interval[QUDA_MAX_MG_LEVEL];
  int                  mg_eig_max_restarts[QUDA_MAX_MG_LEVEL];
  double               mg_eig_tol[QUDA_MAX_MG_LEVEL];
  int                  mg_eig_use_poly_acc[QUDA_MAX_MG_LEVEL];
  int                  mg_eig_poly_deg[QUDA_MAX_MG_LEVEL];
  double               mg_eig_amin[QUDA_MAX_MG_LEVEL];
  double               mg_eig_amax[QUDA_MAX_MG_LEVEL];
  int                  mg_eig_use_normop[QUDA_MAX_MG_LEVEL];
  int                  mg_eig_use_dagger[QUDA_MAX_MG_LEVEL];
  QudaEigSpectrumType  mg_eig_spectrum[QUDA_MAX_MG_LEVEL];
  QudaEigType          mg_eig_type[QUDA_MAX_MG_LEVEL];
  int                  mg_coarse_guess;

  int                  mg_run_low_mode_check;
  int                  mg_run_oblique_proj_check;

} tm_QudaParams_t;

typedef struct tm_QudaMGSetupState_t {
  double gauge_id;
  double c_sw;
  double kappa;
  double mu;
  int initialised;
  double theta_x;
  double theta_y;
  double theta_z;
  double theta_t;
  int force_refresh; // set to 1 when the MG doesn't converge
} tm_QudaMGSetupState_t;

typedef struct tm_QudaCloverState_t {
  double gauge_id;
  double c_sw;
  double kappa;
  double mu;
  int loaded;
  QudaPrecision prec;
  QudaPrecision prec_sloppy;
  QudaPrecision prec_refinement_sloppy;
  QudaPrecision prec_precondition;
  QudaPrecision prec_eigensolver;
  int mg_needs_update;
} tm_QudaCloverState_t;

typedef struct tm_QudaGaugeState_t {
  double gauge_id;
  int loaded;
  double theta_x;
  double theta_y;
  double theta_z;
  double theta_t;
  QudaPrecision prec;
  QudaPrecision prec_sloppy;
  QudaPrecision prec_refinement_sloppy;
  QudaPrecision prec_precondition;
  QudaPrecision prec_eigensolver;
  int mg_needs_update;
} tm_QudaGaugeState_t;

typedef enum tm_QudaMGSetupState_enum_t {
  TM_QUDA_MG_SETUP_RESET = -1,
  TM_QUDA_MG_SETUP_REFRESH,
  TM_QUDA_MG_SETUP_UPDATE,
  TM_QUDA_MG_SETUP_REUSE
} tm_QudaMGSetupState_enum_t; 


static inline int check_quda_clover_state(const tm_QudaCloverState_t * const quda_clover_state,
                                          const tm_QudaGaugeState_t * const quda_gauge_state,
                                          const QudaInvertParam * const inv_param) {
  return( quda_clover_state->loaded &&
          (quda_clover_state->gauge_id == quda_gauge_state->gauge_id) &&
          (fabs(quda_clover_state->c_sw - g_c_sw) < 2*DBL_EPSILON) &&
          (fabs(quda_clover_state->kappa - g_kappa) < 2*DBL_EPSILON) &&
          (fabs(quda_clover_state->mu - g_mu) < 2*DBL_EPSILON) && 
          (quda_clover_state->prec == inv_param->clover_cuda_prec) &&
          (quda_clover_state->prec_sloppy == inv_param->clover_cuda_prec_sloppy) &&
          (quda_clover_state->prec_refinement_sloppy == inv_param->clover_cuda_prec_refinement_sloppy) &&
          (quda_clover_state->prec_precondition == inv_param->clover_cuda_prec_precondition) &&
          (quda_clover_state->prec_eigensolver == inv_param->clover_cuda_prec_eigensolver) );
}

static inline void set_quda_clover_state(tm_QudaCloverState_t * const quda_clover_state,
                                         const tm_QudaGaugeState_t * const quda_gauge_state,
                                         const QudaInvertParam * const inv_param){
  quda_clover_state->gauge_id = quda_gauge_state->gauge_id;
  quda_clover_state->c_sw = g_c_sw;
  quda_clover_state->kappa = g_kappa;
  quda_clover_state->mu = g_mu;
  quda_clover_state->loaded = 1;
  quda_clover_state->prec = inv_param->clover_cuda_prec;
  quda_clover_state->prec_sloppy = inv_param->clover_cuda_prec_sloppy;
  quda_clover_state->prec_refinement_sloppy = inv_param->clover_cuda_prec_refinement_sloppy;
  quda_clover_state->prec_precondition = inv_param->clover_cuda_prec_precondition;
  quda_clover_state->prec_eigensolver = inv_param->clover_cuda_prec_eigensolver;
  quda_clover_state->mg_needs_update = 1;
}

static inline void reset_quda_clover_state(tm_QudaCloverState_t * const quda_clover_state){
  quda_clover_state->gauge_id = -1;
  quda_clover_state->loaded = 0;
  quda_clover_state->mu = -1.0;
  quda_clover_state->c_sw = -1.0;
  quda_clover_state->mu = -1.0;
  quda_clover_state->prec = QUDA_INVALID_PRECISION;
  quda_clover_state->prec_sloppy = QUDA_INVALID_PRECISION;
  quda_clover_state->prec_refinement_sloppy = QUDA_INVALID_PRECISION;
  quda_clover_state->prec_precondition = QUDA_INVALID_PRECISION;
  quda_clover_state->prec_eigensolver = QUDA_INVALID_PRECISION;
  quda_clover_state->mg_needs_update = 1;
}

static inline int check_quda_gauge_state(const tm_QudaGaugeState_t * const quda_gauge_state,
                                         const double gauge_id,
                                         const double theta_x,
                                         const double theta_y,
                                         const double theta_z,
                                         const double theta_t,
                                         const QudaGaugeParam * const gauge_param){
  return( quda_gauge_state->loaded &&
          (fabs(quda_gauge_state->theta_x - theta_x) < 2*DBL_EPSILON) &&
          (fabs(quda_gauge_state->theta_y - theta_y) < 2*DBL_EPSILON) &&
          (fabs(quda_gauge_state->theta_z - theta_z) < 2*DBL_EPSILON) &&
          (fabs(quda_gauge_state->theta_t - theta_t) < 2*DBL_EPSILON) &&
          (fabs(quda_gauge_state->gauge_id - gauge_id) < 2*DBL_EPSILON) &&
          (quda_gauge_state->prec == gauge_param->cuda_prec) &&
          (quda_gauge_state->prec_sloppy == gauge_param->cuda_prec_sloppy) &&
          (quda_gauge_state->prec_refinement_sloppy == gauge_param->cuda_prec_refinement_sloppy) &&
          (quda_gauge_state->prec_precondition == gauge_param->cuda_prec_precondition) &&
          (quda_gauge_state->prec_eigensolver == gauge_param->cuda_prec_eigensolver) );
}

static inline void set_quda_gauge_state(tm_QudaGaugeState_t * const quda_gauge_state,
                                        const double gauge_id,
                                        const double theta_x,
                                        const double theta_y,
                                        const double theta_z,
                                        const double theta_t,
                                        const QudaGaugeParam * const gauge_param){
  quda_gauge_state->gauge_id = gauge_id;
  quda_gauge_state->loaded = 1;
  quda_gauge_state->theta_x = theta_x;
  quda_gauge_state->theta_y = theta_y;
  quda_gauge_state->theta_z = theta_z;
  quda_gauge_state->theta_t = theta_t;
  quda_gauge_state->prec = gauge_param->cuda_prec;
  quda_gauge_state->prec_sloppy = gauge_param->cuda_prec_sloppy;
  quda_gauge_state->prec_refinement_sloppy = gauge_param->cuda_prec_refinement_sloppy;
  quda_gauge_state->prec_precondition = gauge_param->cuda_prec_precondition;
  quda_gauge_state->prec_eigensolver = gauge_param->cuda_prec_eigensolver;
  quda_gauge_state->mg_needs_update = 1;
}

static inline void reset_quda_gauge_state(tm_QudaGaugeState_t * const quda_gauge_state){
  quda_gauge_state->gauge_id = -1;
  quda_gauge_state->loaded = 0;
  quda_gauge_state->prec = QUDA_INVALID_PRECISION;
  quda_gauge_state->prec_sloppy = QUDA_INVALID_PRECISION;
  quda_gauge_state->prec_refinement_sloppy = QUDA_INVALID_PRECISION;
  quda_gauge_state->prec_precondition = QUDA_INVALID_PRECISION;
  quda_gauge_state->prec_eigensolver = QUDA_INVALID_PRECISION;
  quda_gauge_state->mg_needs_update = 1;
}

static inline int check_quda_mg_setup_state(const tm_QudaMGSetupState_t * const quda_mg_setup_state,
                                            const tm_QudaGaugeState_t * const quda_gauge_state,
                                            const tm_QudaCloverState_t * const quda_clover_state,
                                            const tm_QudaParams_t * const quda_params){
  tm_debug_printf(0, 3, "%s mu: %f, g_mu: %f, deltat: %f, reset: %f, refresh: %f\n",
                        __func__,
                        quda_mg_setup_state->mu,
                        g_mu,
                        fabs(quda_mg_setup_state->gauge_id - quda_gauge_state->gauge_id),
                        quda_params->mg_reset_setup_mdu_threshold,
                        quda_params->mg_refresh_setup_mdu_threshold);

  // when the MG setup has not been initialised or when the "gauge_id" has changed by more
  // than the mg_reset_setup_threhold, we need to (re-)do the setup completely
  // similarly, if the boundary conditions for the gauge field change, we need
  // to redo the setup
  if( (quda_mg_setup_state->initialised != 1) ||
      ( fabs(quda_mg_setup_state->theta_x - quda_gauge_state->theta_x) > 2*DBL_EPSILON ) || 
      ( fabs(quda_mg_setup_state->theta_y - quda_gauge_state->theta_y) > 2*DBL_EPSILON ) || 
      ( fabs(quda_mg_setup_state->theta_z - quda_gauge_state->theta_z) > 2*DBL_EPSILON ) || 
      ( fabs(quda_mg_setup_state->theta_t - quda_gauge_state->theta_t) > 2*DBL_EPSILON ) || 
      ( fabs(quda_mg_setup_state->gauge_id - quda_gauge_state->gauge_id) >= quda_params->mg_reset_setup_mdu_threshold )
    ){
    return TM_QUDA_MG_SETUP_RESET;
  // when in the HMC, we have to refresh the setup at regular intervals specified
  // by mg_refresh_setup_mdu_threshold, which triggers a few setup iterations to be
  // run with the existing null vectors as initial guesses, thus refreshing
  // the MG setup for the evolved gauge
  // we also forcibly refresh the setup when the corresponding flag is set in the
  // quda_mg_setup_state 
  } else if ( ( ( fabs(quda_mg_setup_state->gauge_id - quda_gauge_state->gauge_id) < 
                       quda_params->mg_reset_setup_mdu_threshold ) &&
              ( fabs(quda_mg_setup_state->gauge_id - quda_gauge_state->gauge_id) >= 
                       quda_params->mg_refresh_setup_mdu_threshold ) ) ||
              ( quda_mg_setup_state->force_refresh == 1) 
    ){
    return TM_QUDA_MG_SETUP_REFRESH;
  // in other cases, e.g., when the operator parameters change or if the gauge_id has "moved" only a little,
  // we don't need to redo the setup, we can simply rebuild the coarse operators with the
  // new parameters (within reason).
  // Note that we use 2*DBL_EPSILON to have a little bit more wiggle room in case of badly
  // implemented floating point or something like that...
  // TODO: perhaps introduce thresholds also for c_sw, kappa, which might need some
  // more sophisticated logic tree...
  } else if ( ( fabs(quda_mg_setup_state->gauge_id - quda_gauge_state->gauge_id) < 2*DBL_EPSILON ) &&
              ( fabs(quda_mg_setup_state->c_sw - g_c_sw) < 2*DBL_EPSILON) &&
              ( fabs(quda_mg_setup_state->kappa - g_kappa) < 2*DBL_EPSILON) &&
              ( fabs(quda_mg_setup_state->mu - g_mu) < quda_params->mg_reuse_setup_mu_threshold) &&
              ( !quda_gauge_state->mg_needs_update ) &&
              ( quda_clover_state->loaded ? !(quda_clover_state->mg_needs_update) : 1 ) ) {
    return TM_QUDA_MG_SETUP_REUSE;
  } else {
    return TM_QUDA_MG_SETUP_UPDATE;
  }
}

// when we update only the parameters, we are not allowed to touch
// quda_mg_setup_state->gauge_id
static inline void quda_mg_setup_state_update(tm_QudaMGSetupState_t * const quda_mg_setup_state,
                                              tm_QudaGaugeState_t * const quda_gauge_state,
                                              tm_QudaCloverState_t * const quda_clover_state,
                                              const double mu,
                                              const double kappa,
                                              const double c_sw){
  quda_mg_setup_state->mu = mu;
  quda_mg_setup_state->c_sw = c_sw;
  quda_mg_setup_state->kappa = kappa;
  quda_mg_setup_state->force_refresh = 0;
  quda_gauge_state->mg_needs_update = 0;
  quda_clover_state->mg_needs_update = 0;
}

static inline void set_quda_mg_setup_state(tm_QudaMGSetupState_t * const quda_mg_setup_state,
                                           tm_QudaGaugeState_t * const quda_gauge_state,
                                           tm_QudaCloverState_t * const quda_clover_state) {
  quda_mg_setup_state->gauge_id = quda_gauge_state->gauge_id;
  quda_mg_setup_state->theta_x = quda_gauge_state->theta_x;
  quda_mg_setup_state->theta_y = quda_gauge_state->theta_y;
  quda_mg_setup_state->theta_z = quda_gauge_state->theta_z;
  quda_mg_setup_state->theta_t = quda_gauge_state->theta_t;
  quda_mg_setup_state->c_sw = g_c_sw;
  quda_mg_setup_state->kappa = g_kappa;
  quda_mg_setup_state->mu = g_mu;
  quda_mg_setup_state->initialised = 1;
  quda_mg_setup_state->force_refresh = 0;
  quda_gauge_state->mg_needs_update = 0;
  quda_clover_state->mg_needs_update = 0;
}

static inline void reset_quda_mg_setup_state(tm_QudaMGSetupState_t * const quda_mg_setup_state){
  quda_mg_setup_state->gauge_id = -1;
  quda_mg_setup_state->initialised = 0;
  quda_mg_setup_state->force_refresh = 0;
  quda_mg_setup_state->mu = -1.0;
  quda_mg_setup_state->c_sw = -1.0;
  quda_mg_setup_state->mu = -1.0;
}

#endif // TM_QUDA_TYPES_H

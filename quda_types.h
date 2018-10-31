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

#ifdef TM_USE_QUDA
#include <quda.h>
#else
// some definitions which are found in quda.h must be reproduced in case
// we are compiling without it, such that tm_QudaParams_t can be
// defined properly anyway
#define QUDA_MAX_MG_LEVEL 4
#define QUDA_BOOLEAN_YES 1
#define QUDA_BOOLEAN_NO 0
#include "quda_solver_translate.h"
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
  tm_quda_ferm_bc_t fermionbc;

  int               mg_n_level;
  int               mg_n_vec[QUDA_MAX_MG_LEVEL];
  int               mg_blocksize[QUDA_MAX_MG_LEVEL][4];
  double            mg_mu_factor[QUDA_MAX_MG_LEVEL];
  QudaInverterType  mg_setup_inv_type;
  double            mg_setup_tol;
  int               mg_setup_maxiter;
  int 		    mg_coarse_solver_type[QUDA_MAX_MG_LEVEL];
  int               mg_coarse_solver_maxiter;
  double            mg_coarse_solver_tol;
  int               mg_nu_pre;
  int               mg_nu_post;
  int		    mg_smoother_type[QUDA_MAX_MG_LEVEL];
  double            mg_smoother_tol;
  double            mg_omega;
  int               mg_run_verify;
  int               mg_enable_size_three_blocks;
  double            mg_reset_setup_threshold;
  // TODO: extend by mg_smoother_type[QUDA_MAX_MG_LEVEL] and mg_solver[QUDA_MAX_MG_LEVEL]
  // with correct datatype for the purpose, QudaInverterType
} tm_QudaParams_t;

typedef struct tm_QudaMGSetupState_t {
  int gauge_id;
  double c_sw;
  double kappa;
  double mu;
  int initialised;
  double theta_x;
  double theta_y;
  double theta_z;
  double theta_t;
} tm_QudaMGSetupState_t;

typedef struct tm_QudaCloverState_t {
  int gauge_id;
  double c_sw;
  double kappa;
  double mu;
  int loaded;
} tm_QudaCloverState_t;

typedef struct tm_QudaGaugeState_t {
  int gauge_id;
  int loaded;
  double theta_x;
  double theta_y;
  double theta_z;
  double theta_t;
} tm_QudaGaugeState_t;

typedef enum tm_QudaMGSetupState_enum_t {
  TM_QUDA_MG_SETUP_RESET = -1,
  TM_QUDA_MG_SETUP_UPDATE,
  TM_QUDA_MG_SETUP_REUSE
} tm_QudaMGSetupState_enum_t; 


static inline int check_quda_clover_state(const tm_QudaCloverState_t * const quda_clover_state,
                                          const tm_QudaGaugeState_t * const quda_gauge_state){
  return( quda_clover_state->loaded &&
          (quda_clover_state->gauge_id == quda_gauge_state->gauge_id) &&
          (fabs(quda_clover_state->c_sw - g_c_sw) < 2*DBL_EPSILON) &&
          (fabs(quda_clover_state->kappa - g_kappa) < 2*DBL_EPSILON) &&
          (fabs(quda_clover_state->mu - g_mu) < 2*DBL_EPSILON) );
}

static inline void set_quda_clover_state(tm_QudaCloverState_t * const quda_clover_state,
                                         const tm_QudaGaugeState_t * const quda_gauge_state){
  quda_clover_state->gauge_id = quda_gauge_state->gauge_id;
  quda_clover_state->c_sw = g_c_sw;
  quda_clover_state->kappa = g_kappa;
  quda_clover_state->mu = g_mu;
  quda_clover_state->loaded = 1;
}

static inline void reset_quda_clover_state(tm_QudaCloverState_t * const quda_clover_state){
  quda_clover_state->gauge_id = -1;
  quda_clover_state->loaded = 0;
  quda_clover_state->mu = -1.0;
  quda_clover_state->c_sw = -1.0;
  quda_clover_state->mu = -1.0;
}

static inline int check_quda_gauge_state(const tm_QudaGaugeState_t * const quda_gauge_state,
                                         const int gauge_id,
                                         const double theta_x,
                                         const double theta_y,
                                         const double theta_z,
                                         const double theta_t){
  return( quda_gauge_state->loaded &&
          (fabs(quda_gauge_state->theta_x - theta_x) < 2*DBL_EPSILON) &&
          (fabs(quda_gauge_state->theta_y - theta_y) < 2*DBL_EPSILON) &&
          (fabs(quda_gauge_state->theta_z - theta_z) < 2*DBL_EPSILON) &&
          (fabs(quda_gauge_state->theta_t - theta_t) < 2*DBL_EPSILON) &&
          (quda_gauge_state->gauge_id == gauge_id) );
}

static inline void set_quda_gauge_state(tm_QudaGaugeState_t * const quda_gauge_state,
                                        const int gauge_id,
                                        const double theta_x,
                                        const double theta_y,
                                        const double theta_z,
                                        const double theta_t){
  quda_gauge_state->gauge_id = gauge_id;
  quda_gauge_state->loaded = 1;
  quda_gauge_state->theta_x = theta_x;
  quda_gauge_state->theta_y = theta_y;
  quda_gauge_state->theta_z = theta_z;
  quda_gauge_state->theta_t = theta_t;
}

static inline void reset_quda_gauge_state(tm_QudaGaugeState_t * const quda_gauge_state){
  quda_gauge_state->gauge_id = -1;
  quda_gauge_state->loaded = 0;
}

static inline int check_quda_mg_setup_state(const tm_QudaMGSetupState_t * const quda_mg_setup_state,
                                            const tm_QudaGaugeState_t * const quda_gauge_state,
                                            const tm_QudaParams_t * const quda_params){
  // when the MG setup has not been initialised or when the "gauge_id" has changed by more
  // than the mg_redo_setup_threhold, we need to (re-)do the setup completely
  // similarly, if the boundary conditions for the gauge field change, we need
  // to redo the setup
  if( (quda_mg_setup_state->initialised != 1) ||
      ( fabs(quda_mg_setup_state->theta_x - quda_gauge_state->theta_x) > 2*DBL_EPSILON ) || 
      ( fabs(quda_mg_setup_state->theta_y - quda_gauge_state->theta_y) > 2*DBL_EPSILON ) || 
      ( fabs(quda_mg_setup_state->theta_z - quda_gauge_state->theta_z) > 2*DBL_EPSILON ) || 
      ( fabs(quda_mg_setup_state->theta_t - quda_gauge_state->theta_t) > 2*DBL_EPSILON ) || 
      ( fabs(quda_mg_setup_state->gauge_id - quda_gauge_state->gauge_id) > quda_params->mg_reset_setup_threshold ) ){
    return TM_QUDA_MG_SETUP_RESET;
  // in other cases, e.g., when the operator parameters change or if the gauge_id has "moved" only a little,
  // we don't need to redo the setup, we can simply rebuild the coarse operators with the
  // new parameters (within reason).
  // Note that we use 2*DBL_EPSILON to have a little bit more wiggle room in case of badly
  // implemented floating point or something like that...
  // TODO: perhaps introduce thresholds also for c_sw, kappa and mu, which might need some
  // more sophisticated logic tree...
  } else if( ( fabs(quda_mg_setup_state->gauge_id - quda_gauge_state->gauge_id) < 2*DBL_EPSILON ) &&
             ( fabs(quda_mg_setup_state->c_sw - g_c_sw) < 2*DBL_EPSILON) &&
             ( fabs(quda_mg_setup_state->kappa - g_kappa) < 2*DBL_EPSILON) &&
             ( fabs(quda_mg_setup_state->mu - g_mu) < 2*DBL_EPSILON) ){
    return TM_QUDA_MG_SETUP_REUSE;
  } else {
    return TM_QUDA_MG_SETUP_UPDATE;
  }
}

static inline void set_quda_mg_setup_state(tm_QudaMGSetupState_t * const quda_mg_setup_state,
                                           const tm_QudaGaugeState_t * const quda_gauge_state){
  quda_mg_setup_state->gauge_id = quda_gauge_state->gauge_id;
  quda_mg_setup_state->theta_x = quda_gauge_state->theta_x;
  quda_mg_setup_state->theta_y = quda_gauge_state->theta_y;
  quda_mg_setup_state->theta_z = quda_gauge_state->theta_z;
  quda_mg_setup_state->theta_t = quda_gauge_state->theta_t;
  quda_mg_setup_state->c_sw = g_c_sw;
  quda_mg_setup_state->kappa = g_kappa;
  quda_mg_setup_state->mu = g_mu;
  quda_mg_setup_state->initialised = 1;
}

static inline void reset_quda_mg_setup_state(tm_QudaMGSetupState_t * const quda_mg_setup_state){
  quda_mg_setup_state->gauge_id = -1;
  quda_mg_setup_state->initialised = 0;
  quda_mg_setup_state->mu = -1.0;
  quda_mg_setup_state->c_sw = -1.0;
  quda_mg_setup_state->mu = -1.0;
}

#endif // TM_QUDA_TYPES_H

/***********************************************************************
 *
 * Copyright (C) 2015       Mario Schroeck
 *               2016, 2017 Bartosz Kostrzewa
 *               2018       Bartosz Kostrzewa, Ferenc Pittler
 *               2019, 2020 Bartosz Kostrzewa
 *               2021       Bartosz Kostrzewa, Marco Garofalo, Ferenc Pittler, Simone Bacchio
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
 ***********************************************************************/
/***********************************************************************
*
* File quda_interface.h
*
* Authors: Mario Schroeck <mario.schroeck@roma3.infn.it>
*          Bartosz Kostrzewa <bartosz_kostrzewa@fastmail.com>
* 
* Interface to QUDA for multi-GPU inverters
*
* The externally accessible functions are
*
*   void _initQuda()
*     Initializes the QUDA library. Carries over the lattice size and the
*     MPI process grid and thus must be called after initializing MPI.
*     Currently it is called in init_operators() if optr->use_qudainverter
*     flag is set.
*     Memory for the QUDA gaugefield on the host is allocated but not filled
*     yet (the latter is done in _loadGaugeQuda(), see below).
*     Performance critical settings are done here and can be changed.
*
*   void _endQuda()
*     Finalizes the QUDA library. Call before MPI_Finalize().
*
*   void _loadGaugeQuda()
*     Copies and reorders the gaugefield on the host and copies it to the GPU.
*     Must be called between last changes on the gaugefield (smearing etc.)
*     and first call of the inverter. In particular, 'boundary(const double kappa)'
*     must be called before if nontrivial boundary conditions are to be used since
*     those will be applied directly to the gaugefield. Currently it is called just
*     before the inversion is done (might result in wasted loads...).
*     It checks whether the curently loaded gauge field corresponds to the gauge field
*     about to be loaded and returns with a no-op if they agree.
*
*   void _loadCloverQuda()
*     Wrapper for loadCloverQuda() which checks that the currently loaded gauge field
*     and the clover field about to be constructed agree. If they do, the currently
*     loaded clover field is reused.
*
*   void _setQudaMultigridParam()
*     borrowed from QUDA multigrid_invert_test, sets up the input parameters
*     for running the QUDA-MG implementation
*
*   The functions
*
*     int invert_eo_quda(...);
*     int invert_doublet_eo_quda(...);
*     void M_full_quda(...);
*     void D_psi_quda(...);
*
*   mimic their tmLQCD counterparts in functionality as well as input and
*   output parameters. The invert functions will check the parameters
*   g_mu, g_c_sw do decide which QUDA operator to create.
*
*   To activate those, set "UseQudaInverter = yes" in the operator
*   declaration of the input file. For details see the documentation.
*
*   The function
*
*     int invert_quda_direct(...);
*
*   provides a direct interface to the QUDA solver and is not accessible through
*   the input file.
*
* Notes:
*
* Minimum QUDA version is 0.7.0 (see https://github.com/lattice/quda/issues/151 
* and https://github.com/lattice/quda/issues/157).
*
*
**************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "quda_interface.h"
#include "quda_types.h"
#include "boundary.h"
#include "linalg/convert_eo_to_lexic.h"
#include "linalg/mul_r.h"
#include "solver/solver.h"
#include "solver/solver_field.h"
#include "gettime.h"
#include "boundary.h"
#include "quda.h"
#include "global.h"
#include "operator.h"
#include "tm_debug_printf.h"
#include "phmc.h"
#include "quda_gauge_paths.inc"

// nstore is generally like a gauge id, for measurements it identifies the gauge field
// uniquely 
extern int nstore;

// define order of the spatial indices
// default is LX-LY-LZ-T, see below def. of local lattice size, this is related to
// the gamma basis transformation from tmLQCD -> UKQCD
// for details see https://github.com/lattice/quda/issues/157
#define USE_LZ_LY_LX_T 0

#define MAX(a,b) ((a)>(b)?(a):(b))

tm_QudaMGSetupState_t quda_mg_setup_state;
tm_QudaGaugeState_t quda_gauge_state;
tm_QudaCloverState_t quda_clover_state;

// gauge and invert paramameter structs; init. in _initQuda()
QudaGaugeParam  gauge_param;
QudaInvertParam inv_param;
// params to pass to MG
QudaMultigridParam quda_mg_param;
QudaInvertParam mg_inv_param;
void* quda_mg_preconditioner;
// MGEigSolver params for all levels, these need to be around as they
// will be populated in _setQudaMultigridParam and then assigned
// on a per-level basis to QudaInvertParam members (if the EigSolver is enabled
// on the given level)
QudaEigParam mg_eig_param[QUDA_MAX_MG_LEVEL];

// input params specific to tmLQCD QUDA interface
tm_QudaParams_t quda_input;

// pointer to the QUDA gaugefield
double *gauge_quda[4];

QudaGaugeParam  f_gauge_param;
// pointer to the QUDA momentum field
double *mom_quda[4];
double *mom_quda_reordered[4];

// pointer to a temp. spinor, used for reordering etc.
double *tempSpinor;
  
// function that maps coordinates in the communication grid to MPI ranks
int commsMap(const int *coords, void *fdata) {
#if USE_LZ_LY_LX_T
  int n[4] = {coords[3], coords[2], coords[1], coords[0]};
#else
  int n[4] = {coords[3], coords[0], coords[1], coords[2]};
#endif

  int rank = 0;
#ifdef TM_USE_MPI
  MPI_Cart_rank( g_cart_grid, n, &rank );
#endif

  return rank;
}

// variable to check if quda has been initialized
static int quda_initialized = 0;

void _setQudaMultigridParam(QudaMultigridParam* mg_param);
void _setMGInvertParam(QudaInvertParam * mg_inv_param, const QudaInvertParam * const inv_param);
void _updateQudaMultigridPreconditioner(void);
void _setOneFlavourSolverParam(const double kappa, const double c_sw, const double mu, 
                               const int solver_type, const int even_odd,
                               const double eps_sq, const int maxiter,
                               const int single_parity_solve, const int QpQm);

void set_default_gauge_param(QudaGaugeParam * gauge_param){
  // local lattice size
#if USE_LZ_LY_LX_T
  gauge_param->X[0] = LZ;
  gauge_param->X[1] = LY;
  gauge_param->X[2] = LX;
  gauge_param->X[3] = T;
#else
  gauge_param->X[0] = LX;
  gauge_param->X[1] = LY;
  gauge_param->X[2] = LZ;
  gauge_param->X[3] = T;
#endif

  gauge_param->anisotropy = 1.0;
  gauge_param->type = QUDA_WILSON_LINKS;
  gauge_param->gauge_order = QUDA_QDP_GAUGE_ORDER;

  gauge_param->cpu_prec = QUDA_DOUBLE_PRECISION;
  gauge_param->cuda_prec = QUDA_DOUBLE_PRECISION;
  
  gauge_param->reconstruct = NO_COMPRESSION;
  gauge_param->reconstruct_sloppy = NO_COMPRESSION;
  gauge_param->reconstruct_precondition = NO_COMPRESSION;
  gauge_param->reconstruct_refinement_sloppy = NO_COMPRESSION;
  gauge_param->reconstruct_eigensolver = NO_COMPRESSION;
  
  gauge_param->gauge_fix = QUDA_GAUGE_FIXED_NO;

  // For multi-GPU, ga_pad must be large enough to store a time-slice
  int x_face_size = gauge_param->X[1]*gauge_param->X[2]*gauge_param->X[3]/2;
  int y_face_size = gauge_param->X[0]*gauge_param->X[2]*gauge_param->X[3]/2;
  int z_face_size = gauge_param->X[0]*gauge_param->X[1]*gauge_param->X[3]/2;
  int t_face_size = gauge_param->X[0]*gauge_param->X[1]*gauge_param->X[2]/2;
  int pad_size =MAX(x_face_size, y_face_size);
  pad_size = MAX(pad_size, z_face_size);
  pad_size = MAX(pad_size, t_face_size);
  gauge_param->ga_pad = pad_size;
  
  gauge_param->make_resident_gauge = QUDA_BOOLEAN_NO;
  gauge_param->t_boundary = QUDA_PERIODIC_T;
}

void _setDefaultQudaParam(void){
  reset_quda_gauge_state(&quda_gauge_state);
  reset_quda_clover_state(&quda_clover_state);
  reset_quda_mg_setup_state(&quda_mg_setup_state);

  quda_mg_preconditioner = NULL;

  // *** QUDA parameters begin here (sloppy prec. will be adjusted in invert)
  QudaPrecision cpu_prec  = QUDA_DOUBLE_PRECISION;
  QudaPrecision cuda_prec = QUDA_DOUBLE_PRECISION;
  QudaPrecision cuda_prec_sloppy = QUDA_SINGLE_PRECISION;
  QudaPrecision cuda_prec_precondition = QUDA_HALF_PRECISION;

  QudaTune tune = QUDA_TUNE_YES;

  // *** the remainder should not be changed for this application

  inv_param.Ls = 1;

  set_default_gauge_param(&gauge_param);
  
  gauge_param.cuda_prec_sloppy = cuda_prec_sloppy;
  gauge_param.cuda_prec_refinement_sloppy = cuda_prec_sloppy;
  gauge_param.cuda_prec_precondition = cuda_prec_precondition;
  gauge_param.cuda_prec_eigensolver = cuda_prec_precondition;
  
  inv_param.dagger = QUDA_DAG_NO;
  inv_param.mass_normalization = QUDA_KAPPA_NORMALIZATION;
  inv_param.solver_normalization = QUDA_DEFAULT_NORMALIZATION;

  inv_param.pipeline = quda_input.pipeline;
  inv_param.gcrNkrylov = quda_input.gcrNkrylov;

  inv_param.residual_type = (QudaResidualType)(QUDA_L2_RELATIVE_RESIDUAL);
  inv_param.tol_hq = 0.1;
  inv_param.use_alternative_reliable = 0;
  // anything smaller than this and we break down in double-half in most solves
  inv_param.reliable_delta = 1e-3; // ignored by multi-shift solver
  inv_param.use_sloppy_partial_accumulator = 0;

  // domain decomposition preconditioner parameters
  inv_param.inv_type_precondition = QUDA_CG_INVERTER;
  inv_param.schwarz_type = QUDA_ADDITIVE_SCHWARZ;
  inv_param.precondition_cycle = 1;
  inv_param.tol_precondition = 1e-1;
  inv_param.maxiter_precondition = 10;
  inv_param.verbosity_precondition = QUDA_SILENT;
  if( g_debug_level > 3 )
    inv_param.verbosity_precondition = QUDA_SUMMARIZE;
  if( g_debug_level > 5 )
    inv_param.verbosity_precondition = QUDA_VERBOSE;

  inv_param.omega = 1.0;

  inv_param.cpu_prec = cpu_prec;

  inv_param.cuda_prec = cuda_prec;
  inv_param.cuda_prec_sloppy = cuda_prec_sloppy;
  inv_param.cuda_prec_refinement_sloppy = cuda_prec_sloppy;
  inv_param.cuda_prec_precondition = cuda_prec_precondition;
  inv_param.cuda_prec_eigensolver = cuda_prec_precondition;
  
  inv_param.chrono_precision = cuda_prec_sloppy;

  inv_param.clover_rho = 0.0;
  
  inv_param.clover_cpu_prec = cpu_prec;
  inv_param.clover_cuda_prec = cuda_prec;
  inv_param.clover_cuda_prec_sloppy = cuda_prec_sloppy;
  inv_param.clover_cuda_prec_refinement_sloppy = cuda_prec_sloppy;
  inv_param.clover_cuda_prec_precondition = cuda_prec_precondition;
  inv_param.clover_cuda_prec_eigensolver = cuda_prec_precondition;

  inv_param.return_clover = QUDA_BOOLEAN_FALSE;
  inv_param.return_clover_inverse = QUDA_BOOLEAN_FALSE;

  inv_param.preserve_source = QUDA_PRESERVE_SOURCE_YES;
  inv_param.gamma_basis = QUDA_CHIRAL_GAMMA_BASIS; // CHIRAL -> UKQCD does not seem to be supported right now...
  inv_param.dirac_order = QUDA_DIRAC_ORDER;

  inv_param.clover_location = QUDA_CUDA_FIELD_LOCATION;
  
  inv_param.input_location = QUDA_CPU_FIELD_LOCATION;
  inv_param.output_location = QUDA_CPU_FIELD_LOCATION;

  inv_param.tune = tune ? QUDA_TUNE_YES : QUDA_TUNE_NO;

  inv_param.sp_pad = 0; // 24*24*24/2;
  inv_param.cl_pad = 0; // 24*24*24/2;

  // solver verbosity and general verbosity
  QudaVerbosity gen_verb = QUDA_SUMMARIZE;
  if( g_debug_level == 0 ) {
    inv_param.verbosity = QUDA_SILENT;
    gen_verb = QUDA_SUMMARIZE;
  }
  else if( g_debug_level >= 1 && g_debug_level < 3 ) {
    inv_param.verbosity = QUDA_SUMMARIZE;
  }
  else if( g_debug_level >= 3 && g_debug_level < 5 ) {
    inv_param.verbosity = QUDA_VERBOSE;
    gen_verb = QUDA_VERBOSE;
  }
  else if( g_debug_level >= 5 ) {
    inv_param.verbosity = QUDA_DEBUG_VERBOSE;
    gen_verb = QUDA_DEBUG_VERBOSE;
  }

  // general verbosity
  setVerbosityQuda(gen_verb, "# QUDA: ", stdout);
}

void set_force_gauge_param( QudaGaugeParam * f_gauge_param){
  set_default_gauge_param(f_gauge_param);

  f_gauge_param->t_boundary = QUDA_PERIODIC_T;
  f_gauge_param->ga_pad = 0;

  f_gauge_param->use_resident_gauge = QUDA_BOOLEAN_NO;
  f_gauge_param->make_resident_gauge = QUDA_BOOLEAN_NO;
  f_gauge_param->use_resident_mom = QUDA_BOOLEAN_NO;
  f_gauge_param->make_resident_mom = QUDA_BOOLEAN_NO;
  f_gauge_param->return_result_mom = QUDA_BOOLEAN_YES;
  f_gauge_param->overwrite_mom = QUDA_BOOLEAN_YES;
}

void _initQuda() {
  if( quda_initialized )
    return;

  if( g_debug_level > 0 )
    if(g_proc_id == 0)
      printf("\n# TM_QUDA: Detected QUDA version %d.%d.%d\n\n", QUDA_VERSION_MAJOR, QUDA_VERSION_MINOR, QUDA_VERSION_SUBMINOR);
  if( QUDA_VERSION_MAJOR == 0 && QUDA_VERSION_MINOR < 7) {
    fprintf(stderr, "Error: minimum QUDA version required is 0.7.0 (for support of chiral basis and removal of bug in mass normalization with preconditioning).\n");
    exit(-2);
  }

  gauge_param = newQudaGaugeParam();
  f_gauge_param = newQudaGaugeParam();
  inv_param = newQudaInvertParam();
  mg_inv_param = newQudaInvertParam();
  quda_mg_param = newQudaMultigridParam();
  for( int level = 0; level < QUDA_MAX_MG_LEVEL; ++level ){
    mg_eig_param[level] = newQudaEigParam();
  }

  _setDefaultQudaParam();

  // declare the grid mapping used for communications in a multi-GPU grid
#if USE_LZ_LY_LX_T
  int grid[4] = {g_nproc_z, g_nproc_y, g_nproc_x, g_nproc_t};
#else
  int grid[4] = {g_nproc_x, g_nproc_y, g_nproc_z, g_nproc_t};
#endif

  initCommsGridQuda(4, grid, commsMap, NULL);

  // alloc gauge_quda
  size_t gSize = (gauge_param.cpu_prec == QUDA_DOUBLE_PRECISION) ? sizeof(double) : sizeof(float);

  for (int dir = 0; dir < 4; dir++) {
    gauge_quda[dir] = (double*) malloc(VOLUME*18*gSize);
    if(gauge_quda[dir] == NULL) {
      fprintf(stderr, "_initQuda: malloc for gauge_quda[dir] failed");
      exit(-2);
    }
  }

  // alloc space for a temp. spinor, used throughout this module
  tempSpinor  = (double*)malloc( 2*VOLUME*24*sizeof(double) ); /* factor 2 for doublet */
  if(tempSpinor == NULL) {
    fprintf(stderr, "_initQuda: malloc for tempSpinor failed");
    exit(-2);
  }

  // initialize the QUDA library
#ifdef TM_USE_MPI
  initQuda(-1); //sets device numbers automatically
#else
  // when running in 'subprocess' mode, the external program should have provided us with a unique
  // id in the range 0 to (N-1), where N is the number of NVIDIA devices available (see wrapper/lib_wrapper.c)
  if(subprocess_flag){
    initQuda(g_external_id);
  }else{
    initQuda(0);  //scalar build without subprocess: use device 0
  }
#endif
  quda_initialized = 1;
}

// finalize the QUDA library
void _endQuda() {
  if( quda_initialized ) {
    if( quda_mg_preconditioner != NULL ){
      destroyMultigridQuda(quda_mg_preconditioner);
      quda_mg_preconditioner = NULL;
    }
    for(int dir = 0; dir < 4; dir++){
      if( (void*)gauge_quda[dir] != NULL ) free((void*)gauge_quda[dir]);
      if( (void*)mom_quda[dir] != NULL ) free((void*)mom_quda[dir]);
      if( (void*)mom_quda_reordered[dir] != NULL ) free((void*)mom_quda_reordered[dir]);
    }
    freeGaugeQuda();
    freeCloverQuda(); // this is safe even if there is no Clover field loaded, at least it was in QUDA v0.7.2
    free((void*)tempSpinor);
    endQuda();
  }
}

void _loadCloverQuda(QudaInvertParam* inv_param){
  static int first_call = 1;
  // check if loaded clover and gauge fields agree
  if( check_quda_clover_state(&quda_clover_state, &quda_gauge_state) ){
      tm_debug_printf(0, 0, "# TM_QUDA: Clover field and inverse already loaded for gauge_id: %f\n", quda_gauge_state.gauge_id);
  } else {
    tm_stopwatch_push(&g_timers, "loadCloverQuda", "");
    if(first_call){
      first_call = 1;
    } else {
      freeCloverQuda();
    }
    reset_quda_clover_state(&quda_clover_state);
    loadCloverQuda(NULL, NULL, inv_param);
    set_quda_clover_state(&quda_clover_state, &quda_gauge_state);
    tm_stopwatch_pop(&g_timers, 0, 0, "TM_QUDA");
  }
}

void reorder_gauge_toQuda( const su3 ** const gaugefield, const CompressionType compression ) {
  tm_stopwatch_push(&g_timers, __func__, "");

#ifdef TM_USE_OMP
#pragma omp parallel
  {
#endif
  _Complex double tmpcplx;

  size_t gSize = (gauge_param.cpu_prec == QUDA_DOUBLE_PRECISION) ? sizeof(double) : sizeof(float);
  
  // now copy and reorder
#ifdef TM_USE_OMP
  #pragma omp for collapse(4)
#endif
  for( int x0=0; x0<T; x0++ )
    for( int x1=0; x1<LX; x1++ )
      for( int x2=0; x2<LY; x2++ )
        for( int x3=0; x3<LZ; x3++ ) {
#if USE_LZ_LY_LX_T
          int j = x3 + LZ*x2 + LY*LZ*x1 + LX*LY*LZ*x0;
          int tm_idx = x1 + LX*x2 + LY*LX*x3 + LZ*LY*LX*x0;
#else
          int j = x1 + LX*x2 + LY*LX*x3 + LZ*LY*LX*x0;
          int tm_idx = x3 + LZ*x2 + LY*LZ*x1 + LX*LY*LZ*x0;
#endif
          int oddBit = (x0+x1+x2+x3) & 1;
          int quda_idx = 18*(oddBit*VOLUME/2+j/2);

#if USE_LZ_LY_LX_T
          memcpy( &(gauge_quda[0][quda_idx]), &(gaugefield[tm_idx][3]), 18*gSize);
          memcpy( &(gauge_quda[1][quda_idx]), &(gaugefield[tm_idx][2]), 18*gSize);
          memcpy( &(gauge_quda[2][quda_idx]), &(gaugefield[tm_idx][1]), 18*gSize);
          memcpy( &(gauge_quda[3][quda_idx]), &(gaugefield[tm_idx][0]), 18*gSize);
#else
          memcpy( &(gauge_quda[0][quda_idx]), &(gaugefield[tm_idx][1]), 18*gSize);
          memcpy( &(gauge_quda[1][quda_idx]), &(gaugefield[tm_idx][2]), 18*gSize);
          memcpy( &(gauge_quda[2][quda_idx]), &(gaugefield[tm_idx][3]), 18*gSize);
          memcpy( &(gauge_quda[3][quda_idx]), &(gaugefield[tm_idx][0]), 18*gSize);
#endif
        if( compression == NO_COMPRESSION && quda_input.fermionbc == TM_QUDA_THETABC ) {
          // apply theta boundary conditions if compression is not used
          for( int i=0; i<9; i++ ) {
            tmpcplx = gauge_quda[0][quda_idx+2*i] + I*gauge_quda[0][quda_idx+2*i+1];
            tmpcplx *= -phase_1/g_kappa;
            gauge_quda[0][quda_idx+2*i]   = creal(tmpcplx);
            gauge_quda[0][quda_idx+2*i+1] = cimag(tmpcplx);

            tmpcplx = gauge_quda[1][quda_idx+2*i] + I*gauge_quda[1][quda_idx+2*i+1];
            tmpcplx *= -phase_2/g_kappa;
            gauge_quda[1][quda_idx+2*i]   = creal(tmpcplx);
            gauge_quda[1][quda_idx+2*i+1] = cimag(tmpcplx);

            tmpcplx = gauge_quda[2][quda_idx+2*i] + I*gauge_quda[2][quda_idx+2*i+1];
            tmpcplx *= -phase_3/g_kappa;
            gauge_quda[2][quda_idx+2*i]   = creal(tmpcplx);
            gauge_quda[2][quda_idx+2*i+1] = cimag(tmpcplx);

            tmpcplx = gauge_quda[3][quda_idx+2*i] + I*gauge_quda[3][quda_idx+2*i+1];
            tmpcplx *= -phase_0/g_kappa;
            gauge_quda[3][quda_idx+2*i]   = creal(tmpcplx);
            gauge_quda[3][quda_idx+2*i+1] = cimag(tmpcplx);
          }
          // when compression is not used, we can still force naive anti-periodic boundary conditions
        } else {
          if ( quda_input.fermionbc == TM_QUDA_APBC && x0+g_proc_coords[0]*T == g_nproc_t*T-1 ) {
            for( int i=0; i<18; i++ ) {
              gauge_quda[3][quda_idx+i]   = -gauge_quda[3][quda_idx+i];
            }
          } // quda_input.fermionbc
        } // if(compression & boundary conditions)
      } // volume loop
#ifdef TM_USE_OMP
  } // OpenMP parallel closing brace 
#endif

  tm_stopwatch_pop(&g_timers, 0, 0, "TM_QUDA");
}

void _loadGaugeQuda( const CompressionType compression ) {
  static int first_call = 1;
  // check if the currently loaded gauge field is also the current gauge field
  // and if so, return immediately
  tm_debug_printf(0, 1, "# TM_QUDA: Called _loadGaugeQuda for gauge_id: %f\n", g_gauge_state.gauge_id);
  
  if( inv_param.verbosity > QUDA_SILENT ){
    if(g_proc_id == 0) {
      if( compression == NO_COMPRESSION ){
        if( quda_input.fermionbc == TM_QUDA_THETABC ){
          printf("# TM_QUDA: Theta boundary conditions will be applied to gauge field\n");
        }
      } else {
        if( quda_input.fermionbc == TM_QUDA_APBC ){
          printf("# TM_QUDA: Temporal ABPC will be applied to gauge field\n");
        }
      }
    }
  }

  if( check_quda_gauge_state(&quda_gauge_state, g_gauge_state.gauge_id, X1, X2, X3, X0) ){
    return;
  } else {
    if( first_call ){
      first_call = 0;
    } else {
      freeGaugeQuda();
    }
    reset_quda_gauge_state(&quda_gauge_state);
  }

  reorder_gauge_toQuda(g_gauge_field, compression);

  tm_stopwatch_push(&g_timers, "loadGaugeQuda", "");
  loadGaugeQuda((void*)gauge_quda, &gauge_param);
  tm_stopwatch_pop(&g_timers, 0, 0, "TM_QUDA");

  set_quda_gauge_state(&quda_gauge_state, g_gauge_state.gauge_id, X1, X2, X3, X0);
}

// reorder spinor to QUDA format
void reorder_spinor_toQuda( double* sp, QudaPrecision precision, int doublet ) {
  tm_stopwatch_push(&g_timers, __func__, "");

  memcpy( tempSpinor, sp, (1+doublet)*VOLUME*24*sizeof(double) );

  // now copy and reorder from tempSpinor to spinor
#ifdef TM_USE_OMP
  #pragma omp parallel for collapse(4)
#endif
  for( int x0=0; x0<T; x0++ )
    for( int x1=0; x1<LX; x1++ )
      for( int x2=0; x2<LY; x2++ )
        for( int x3=0; x3<LZ; x3++ ) {
#if USE_LZ_LY_LX_T
          int j = x3 + LZ*x2 + LY*LZ*x1 + LX*LY*LZ*x0;
          int tm_idx = x1 + LX*x2 + LY*LX*x3 + LZ*LY*LX*x0;
#else
          int j = x1 + LX*x2 + LY*LX*x3 + LZ*LY*LX*x0;
          int tm_idx   = x3 + LZ*x2 + LY*LZ*x1 + LX*LY*LZ*x0;
#endif
          int oddBit = (x0+x1+x2+x3) & 1;

          if( doublet ) {
            memcpy( &(sp[24*(oddBit*VOLUME+j/2)]),          &(tempSpinor[24*tm_idx         ]), 24*sizeof(double));
            memcpy( &(sp[24*(oddBit*VOLUME+j/2+VOLUME/2)]), &(tempSpinor[24*(tm_idx+VOLUME)]), 24*sizeof(double));
          }
          else {
            memcpy( &(sp[24*(oddBit*VOLUME/2+j/2)]), &(tempSpinor[24*tm_idx]), 24*sizeof(double));
          }

        }

  tm_stopwatch_pop(&g_timers, 0, 0, "TM_QUDA");
}

void _initMomQuda(void) {
  static int first_call = 1;
  if( first_call ){
    first_call = 0;
    set_force_gauge_param(&f_gauge_param);
    for(int i = 0; i < 4; i++){
      mom_quda[i] = (double*)malloc(VOLUME*10*sizeof(double));
      mom_quda_reordered[i] = (double*)malloc(VOLUME*10*sizeof(double));
      if( (void*)mom_quda[i] == NULL || (void*)mom_quda_reordered[i] == NULL ){
        fatal_error("Memory allocation for host momentum field failed!", __func__);
      }
    }
  }
}

void reorder_mom_fromQuda() {
  // mom_quda -> mom_quda_reordered
  tm_stopwatch_push(&g_timers, __func__, "");

#ifdef TM_USE_OMP
#pragma omp parallel
  {
#endif

#ifdef TM_USE_OMP
  #pragma omp parallel for collapse(4)
#endif
  for( int x0=0; x0<T; x0++ )
    for( int x1=0; x1<LX; x1++ )
      for( int x2=0; x2<LY; x2++ )
        for( int x3=0; x3<LZ; x3++ ) {
#if USE_LZ_LY_LX_T
          int j = x3 + LZ*x2 + LY*LZ*x1 + LX*LY*LZ*x0;
          int tm_idx = x1 + LX*x2 + LY*LX*x3 + LZ*LY*LX*x0;
#else
          int j = x1 + LX*x2 + LY*LX*x3 + LZ*LY*LX*x0;
          int tm_idx   = x3 + LZ*x2 + LY*LZ*x1 + LX*LY*LZ*x0;
#endif
          int oddBit = (x0+x1+x2+x3) & 1;
          int quda_idx = 10*(oddBit*VOLUME/2+j/2);
          tm_idx *= 10;

#if USE_LZ_LY_LX_T
          memcpy( &(mom_quda_reordered[3][tm_idx]), &(mom_quda[0][quda_idx]), 10*sizeof(double));
          memcpy( &(mom_quda_reordered[2][tm_idx]), &(mom_quda[1][quda_idx]), 10*sizeof(double));
          memcpy( &(mom_quda_reordered[1][tm_idx]), &(mom_quda[2][quda_idx]), 10*sizeof(double));
          memcpy( &(mom_quda_reordered[0][tm_idx]), &(mom_quda[3][quda_idx]), 10*sizeof(double));
#else
          memcpy( &(mom_quda_reordered[1][tm_idx]), &(mom_quda[0][quda_idx]), 10*sizeof(double));
          memcpy( &(mom_quda_reordered[2][tm_idx]), &(mom_quda[1][quda_idx]), 10*sizeof(double));
          memcpy( &(mom_quda_reordered[3][tm_idx]), &(mom_quda[2][quda_idx]), 10*sizeof(double));
          memcpy( &(mom_quda_reordered[0][tm_idx]), &(mom_quda[3][quda_idx]), 10*sizeof(double));
#endif
        }
  
#ifdef TM_USE_OMP
  }
#endif
  tm_stopwatch_pop(&g_timers, 0, 0, "TM_QUDA");
}


// reorder spinor from QUDA format
void reorder_spinor_fromQuda( double* sp, QudaPrecision precision, int doublet) {
  tm_stopwatch_push(&g_timers, __func__, "");

  memcpy( tempSpinor, sp, (1+doublet)*VOLUME*24*sizeof(double) );

  // now copy and reorder from tempSpinor to spinor
#ifdef TM_USE_OMP
  #pragma omp parallel for collapse(4)
#endif
  for( int x0=0; x0<T; x0++ )
    for( int x1=0; x1<LX; x1++ )
      for( int x2=0; x2<LY; x2++ )
        for( int x3=0; x3<LZ; x3++ ) {
#if USE_LZ_LY_LX_T
          int j = x3 + LZ*x2 + LY*LZ*x1 + LX*LY*LZ*x0;
          int tm_idx = x1 + LX*x2 + LY*LX*x3 + LZ*LY*LX*x0;
#else
          int j = x1 + LX*x2 + LY*LX*x3 + LZ*LY*LX*x0;
          int tm_idx   = x3 + LZ*x2 + LY*LZ*x1 + LX*LY*LZ*x0;
#endif
          int oddBit = (x0+x1+x2+x3) & 1;

          if( doublet ) {
            memcpy( &(sp[24*tm_idx]),          &(tempSpinor[24*(oddBit*VOLUME+j/2)         ]), 24*sizeof(double));
            memcpy( &(sp[24*(tm_idx+VOLUME)]), &(tempSpinor[24*(oddBit*VOLUME+j/2+VOLUME/2)]), 24*sizeof(double));
          }
          else {
            memcpy( &(sp[24*tm_idx]), &(tempSpinor[24*(oddBit*VOLUME/2+j/2)]), 24*sizeof(double));
          }
        }
  tm_stopwatch_pop(&g_timers, 0, 0, "TM_QUDA");
}


// reorder spinor to QUDA format
void reorder_spinor_eo_toQuda(double* sp, QudaPrecision precision, int doublet, int odd) {
  tm_stopwatch_push(&g_timers, __func__, "");

  static const int change_sign[4] = {-1, 1, 1, -1};
  static const int change_spin[4] = {3, 2, 1, 0};
  const int Vh = VOLUME/2;

  memcpy(tempSpinor, sp, (1+doublet)*Vh*24*sizeof(double) );

  // now copy and reorder from tempSpinor to spinor
#ifdef TM_USE_OMP
  #pragma omp parallel for collapse(4)
#endif
  for( int x0=0; x0<T; x0++ )
    for( int x1=0; x1<LX; x1++ )
      for( int x2=0; x2<LY; x2++ )
        for( int x3=0; x3<LZ; x3++ ) {
#if USE_LZ_LY_LX_T
          const int q_eo_idx = (x3 + LZ*x2 + LY*LZ*x1 + LX*LY*LZ*x0)/2;
          const int tm_eo_idx = (x1 + LX*x2 + LY*LX*x3 + LZ*LY*LX*x0)/2;
#else
          const int q_eo_idx = (x1 + LX*x2 + LY*LX*x3 + LZ*LY*LX*x0)/2;
          const int tm_eo_idx  = (x3 + LZ*x2 + LY*LZ*x1 + LX*LY*LZ*x0)/2;
#endif
          const int oddBit = (x0+x1+x2+x3) & 1;
          if( oddBit == odd ){
            for(int q_spin = 0; q_spin < 4; q_spin++){
              const int tm_spin = change_spin[q_spin];
              for(int col = 0; col < 3; col++){
                for(int reim = 0; reim < 2; reim++){
                  sp[24*q_eo_idx + 6*q_spin + 2*col + reim] = 
                    change_sign[q_spin] * tempSpinor[24*tm_eo_idx + 6*tm_spin + 2*col + reim];
                  if(doublet){
                    sp[24*(q_eo_idx+Vh) + 6*q_spin + 2*col + reim] = 
                      change_sign[q_spin] * tempSpinor[24*(tm_eo_idx+Vh) + 6*tm_spin + 2*col + reim];
                  }
                }
              }
            }
          }
        }
  tm_stopwatch_pop(&g_timers, 0, 0, "TM_QUDA");
}

// reorder spinor from QUDA format
void reorder_spinor_eo_fromQuda( double* sp, QudaPrecision precision, int doublet, int odd) {
  tm_stopwatch_push(&g_timers, __func__, "");

  const int change_sign[4] = {-1, 1, 1, -1};
  const int change_spin[4] = {3, 2, 1, 0};
  const int Vh = VOLUME/2;


  memcpy( tempSpinor, sp, (1+doublet)*(VOLUME/2)*24*sizeof(double) );

  // now copy and reorder from tempSpinor to spinor
#ifdef TM_USE_OMP
  #pragma omp parallel for collapse(4)
#endif
  for( int x0=0; x0<T; x0++ )
    for( int x1=0; x1<LX; x1++ )
      for( int x2=0; x2<LY; x2++ )
        for( int x3=0; x3<LZ; x3++ ) {
#if USE_LZ_LY_LX_T
          const int q_eo_idx = (x3 + LZ*x2 + LY*LZ*x1 + LX*LY*LZ*x0)/2;
          const int tm_eo_idx = (x1 + LX*x2 + LY*LX*x3 + LZ*LY*LX*x0)/2;
#else
          const int q_eo_idx = (x1 + LX*x2 + LY*LX*x3 + LZ*LY*LX*x0)/2;
          const int tm_eo_idx = (x3 + LZ*x2 + LY*LZ*x1 + LX*LY*LZ*x0)/2;
#endif
          const int oddBit = (x0+x1+x2+x3) & 1;
          if( oddBit == odd ){
            for(int q_spin = 0; q_spin < 4; q_spin++){
              const int tm_spin = change_spin[q_spin];
              for(int col = 0; col < 3; col++){
                for(int reim = 0; reim < 2; reim++){
                  sp[24*tm_eo_idx + 6*tm_spin + 2*col + reim] = 
                    change_sign[q_spin] * tempSpinor[24*q_eo_idx + 6*q_spin + 2*col + reim];
                  if(doublet){
                    sp[24*(tm_eo_idx+Vh) + 6*tm_spin + 2*col + reim] = 
                      change_sign[q_spin] * tempSpinor[24*(q_eo_idx+Vh) + 6*q_spin + 2*col + reim];
                  }
                }
              }
            }
          }
        }

  tm_stopwatch_pop(&g_timers, 0, 0, "TM_QUDA");
}


void set_boundary_conditions(CompressionType* compression, QudaGaugeParam * gauge_param) {
  // we can't have compression and theta-BC, but we will support compression
  // for theta_0 = 0.0 or theta_0 = 1.0 (using naive periodic or anti-periodic boundary conditions
  // warning the user that the residual check will fail)
  if( fabs(X1)>0.0 || fabs(X2)>0.0 || fabs(X3)>0.0 || (fabs(X0) > 2*DBL_EPSILON && fabs(fabs(X0)-1.0) > 2*DBL_EPSILON  ) ) {
    if( *compression!=NO_COMPRESSION ) {
      if(g_proc_id == 0) {
        printf("\n# TM_QUDA: WARNING you can't use compression %d with boundary conditions for fermion fields (t,x,y,z)*pi: (%f,%f,%f,%f) \n", *compression,X0,X1,X2,X3);
        printf("# TM_QUDA: disabling compression.\n\n");
      }
      *compression=NO_COMPRESSION;
    }
  }

  QudaReconstructType link_recon;
  QudaReconstructType link_recon_sloppy;

  if( *compression==NO_COMPRESSION ) {
    // without compression, any kind of boundary conditions are supported
    // and will be applied to the gauge field as required
    if( quda_input.fermionbc != TM_QUDA_APBC ){
      gauge_param->t_boundary = QUDA_PERIODIC_T;
    } else {
      gauge_param->t_boundary = QUDA_ANTI_PERIODIC_T;
    }
    link_recon = NO_COMPRESSION;
    link_recon_sloppy = NO_COMPRESSION;
  } else {
    // if we reach this point with compression (see logic above), theta_0 is either 0.0 or 1.0
    // if it is 1.0, we explicitly enabled TM_QUDA_APBC to force simple anti-periodic boundary
    // conditions
    if( fabs(X0)>0.0 ){
      quda_input.fermionbc = TM_QUDA_APBC;
      tm_debug_printf(0, 0,
          "# TM_QUDA: WARNING You have set temporal theta-BC but gauge compression is enabled. "
          "This will be overriden to use naive APBC instead. This works fine, but the residual "
          "check on the host (CPU) will fail.\n");
    }

    if( quda_input.fermionbc == TM_QUDA_APBC ) {
      gauge_param->t_boundary = QUDA_ANTI_PERIODIC_T;
    } else {
      gauge_param->t_boundary = QUDA_PERIODIC_T;
    }

    link_recon = COMPRESSION_12;
    link_recon_sloppy = *compression;

    tm_debug_printf(0, 0, 
        "\n# TM_QUDA: WARNING using %d compression with trivial (A)PBC instead "
        "of theta-BC ((t,x,y,z)*pi: (%f,%f,%f,%f))! This works fine but the residual "
        "check on the host (CPU) will fail.\n",
        *compression,X0,X1,X2,X3);
  }

  gauge_param->reconstruct = link_recon;
  gauge_param->reconstruct_sloppy = link_recon_sloppy;
  gauge_param->reconstruct_precondition = link_recon_sloppy;
}

void set_sloppy_prec(const SloppyPrecision sloppy_precision, QudaGaugeParam * gauge_param, QudaInvertParam * inv_param) {
  // choose sloppy prec.
  QudaPrecision cuda_prec_sloppy;
  if( sloppy_precision==SLOPPY_DOUBLE ) {
    inv_param->reliable_delta = 1e-8;
    cuda_prec_sloppy = QUDA_DOUBLE_PRECISION;
    if(g_proc_id == 0) printf("# TM_QUDA: Using double prec. as sloppy!\n");
  }
  else if( sloppy_precision==SLOPPY_HALF ) {
    // in double-half, we perform many reliable updates
    inv_param->reliable_delta = 1e-2;
    cuda_prec_sloppy = QUDA_HALF_PRECISION;
    if(g_proc_id == 0) printf("# TM_QUDA: Using half prec. as sloppy!\n");
  }
  else {
    inv_param->reliable_delta = 1e-4;
    cuda_prec_sloppy = QUDA_SINGLE_PRECISION;
    if(g_proc_id == 0) printf("# TM_QUDA: Using single prec. as sloppy!\n");
  }
  
  gauge_param->cuda_prec_sloppy = cuda_prec_sloppy;
  gauge_param->cuda_prec_refinement_sloppy = cuda_prec_sloppy;
  
  inv_param->cuda_prec_sloppy = cuda_prec_sloppy;
  inv_param->cuda_prec_refinement_sloppy = cuda_prec_sloppy;
  
  inv_param->clover_cuda_prec_sloppy = cuda_prec_sloppy;
  inv_param->clover_cuda_prec_refinement_sloppy = cuda_prec_sloppy;
}



int invert_quda_direct(double * const propagator, double const * const source,
                       const int op_id) {
  tm_stopwatch_push(&g_timers, __func__, ""); 
  spinor ** solver_field = NULL;
  init_solver_field(&solver_field, VOLUME, 1);

  memcpy((void*)(solver_field[0]), (void*)(source), VOLUME*sizeof(spinor));

  double atime, atotaltime = gettime();
  void *spinorIn  = (void*)solver_field[0]; // source
  void *spinorOut = (void*)propagator; // solution

  operator * optr = &operator_list[op_id];
  // g_kappa is necessary for the gauge field to be correctly translated from tmLQCD to QUDA
  g_kappa = optr->kappa;
  g_c_sw = optr->c_sw;
  g_mu = optr->mu;

  boundary(optr->kappa);
  
  if ( g_relative_precision_flag )
    inv_param.residual_type = QUDA_L2_RELATIVE_RESIDUAL;
  else
    inv_param.residual_type = QUDA_L2_ABSOLUTE_RESIDUAL;
  
  inv_param.kappa = optr->kappa;

  // figure out which BC to use (theta, trivial...)
  set_boundary_conditions(&optr->compression_type, &gauge_param);

  // set the sloppy precision of the mixed prec solver
  set_sloppy_prec(optr->sloppy_precision, &gauge_param, &inv_param);
 
  // load gauge after setting precision, this is a no-op if the current gauge field
  // is already loaded and the boundary conditions have not changed
  atime = gettime();
  _loadGaugeQuda(optr->compression_type);
  if(g_proc_id==0 && g_debug_level > 0 ) printf("# TM_QUDA: Time for loadGaugeQuda: %.4e\n",gettime()-atime);

  // this will also construct the clover field and its inverse, if required
  // it will also run the MG setup
  _setOneFlavourSolverParam(optr->kappa, 
                            optr->c_sw, 
                            optr->mu, 
                            optr->solver,
                            optr->even_odd_flag,
                            optr->eps_sq,
                            optr->maxiter,
                            0, 0);
  
  // reorder spinor
  reorder_spinor_toQuda( (double*)spinorIn, inv_param.cpu_prec, 0 );

  // perform the inversion
  invertQuda(spinorOut, spinorIn, &inv_param);

  if( inv_param.verbosity > QUDA_SILENT )
    if(g_proc_id == 0)
      printf("# TM_QUDA: Done: %i iter / %g secs = %g Gflops\n",
             inv_param.iter, inv_param.secs, inv_param.gflops/inv_param.secs);

  optr->iterations = inv_param.iter;

  // reorder spinor
  reorder_spinor_fromQuda( (double*)spinorOut, inv_param.cpu_prec, 0 );
  // propagator in usual normalisation, this is only necessary in invert_quda_direct
  // since the rescaling is otherwise done in the operator inversion driver
  mul_r((spinor*)spinorOut, (2*optr->kappa), (spinor*)spinorOut, VOLUME );

  finalize_solver(solver_field, 1);

  tm_stopwatch_pop(&g_timers, 0, 0, "TM_QUDA");

  if(optr->iterations >= optr->maxiter)
    return(-1);

  return(optr->iterations);
}

int invert_eo_quda(spinor * const Even_new, spinor * const Odd_new,
                   spinor * const Even, spinor * const Odd,
                   const double precision, const int max_iter,
                   const int solver_flag, const int rel_prec,
                   const int even_odd_flag, solver_params_t solver_params,
                   SloppyPrecision sloppy_precision,
                   CompressionType compression) {

  tm_stopwatch_push(&g_timers, __func__, "");

  spinor ** solver_field = NULL;
  const int nr_sf = 2;
  init_solver_field(&solver_field, VOLUME, nr_sf);

  convert_eo_to_lexic(solver_field[0],  Even, Odd);

// this is basically not necessary, but if we want to use an a nitial guess, it will be
//  convert_eo_to_lexic(solver_field[1], Even_new, Odd_new);

  void *spinorIn  = (void*)solver_field[0]; // source
  void *spinorOut = (void*)solver_field[1]; // solution

  if ( rel_prec )
    inv_param.residual_type = QUDA_L2_RELATIVE_RESIDUAL;
  else
    inv_param.residual_type = QUDA_L2_ABSOLUTE_RESIDUAL;

  inv_param.kappa = g_kappa;

  // figure out which BC to use (theta, trivial...)
  set_boundary_conditions(&compression, &gauge_param);
  // set the sloppy precision of the mixed prec solver
  set_sloppy_prec(sloppy_precision, &gauge_param, &inv_param);
  
  // load gauge after setting precision
  _loadGaugeQuda(compression);

  // this will also construct the clover field and its inverse, if required
  // it will also run the MG setup
  _setOneFlavourSolverParam(g_kappa,
                            g_c_sw,
                            g_mu,
                            solver_flag,
                            even_odd_flag,
                            precision,
                            max_iter,
                            0, 0);

  // reorder spinor
  reorder_spinor_toQuda( (double*)spinorIn, inv_param.cpu_prec, 0 );

  // perform the inversion
  tm_stopwatch_push(&g_timers, "invertQuda", "");
  invertQuda(spinorOut, spinorIn, &inv_param);
  tm_stopwatch_pop(&g_timers, 0, 1, "TM_QUDA");


  if( inv_param.verbosity > QUDA_SILENT )
    if(g_proc_id == 0)
      printf("# TM_QUDA: Done: %i iter / %g secs = %g Gflops\n",
             inv_param.iter, inv_param.secs, inv_param.gflops/inv_param.secs);

  // number of CG iterations
  int iteration = inv_param.iter;

  reorder_spinor_fromQuda( (double*)spinorOut, inv_param.cpu_prec, 0 );
  convert_lexic_to_eo(Even_new, Odd_new, solver_field[1]);

  finalize_solver(solver_field, nr_sf);

  tm_stopwatch_pop(&g_timers, 0, 0, "TM_QUDA");

  if(iteration >= max_iter)
    return(-1);

  return(iteration);
}

int invert_doublet_eo_quda(spinor * const Even_new_s, spinor * const Odd_new_s,
                           spinor * const Even_new_c, spinor * const Odd_new_c,
                           spinor * const Even_s, spinor * const Odd_s,
                           spinor * const Even_c, spinor * const Odd_c,
                           const double precision, const int max_iter,
                           const int solver_flag, const int rel_prec, const int even_odd_flag,
                           const SloppyPrecision sloppy_precision,
                           CompressionType compression) {
  tm_stopwatch_push(&g_timers, __func__, "");

  spinor ** solver_field = NULL;
  const int nr_sf = 2;
  init_solver_field(&solver_field, 2*VOLUME, nr_sf);

  convert_eo_to_lexic(solver_field[0],          Even_s,  Odd_s);
  convert_eo_to_lexic(solver_field[0]+VOLUME,   Even_c,  Odd_c);

  // if we were to use an initial guess, we would need to also prepare the
  // solution spinor here

  void *spinorIn    = (void*)solver_field[0];
  void *spinorOut   = (void*)solver_field[1];

  if ( rel_prec )
    inv_param.residual_type = QUDA_L2_RELATIVE_RESIDUAL;
  else
    inv_param.residual_type = QUDA_L2_ABSOLUTE_RESIDUAL;

  inv_param.kappa = g_kappa;

  // IMPORTANT: use opposite TM mu-flavor since gamma5 -> -gamma5
  inv_param.mu           = -g_mubar /2./g_kappa;
  inv_param.epsilon      =  g_epsbar/2./g_kappa;
  // FIXME: in principle, there is also QUDA_TWIST_DEG_DOUBLET
  inv_param.twist_flavor =  QUDA_TWIST_NONDEG_DOUBLET; 
  inv_param.Ls = 2;

  // figure out which BC to use (theta, trivial...)
  set_boundary_conditions(&compression, &gauge_param);

  // set the sloppy precision of the mixed prec solver
  set_sloppy_prec(sloppy_precision, &gauge_param, &inv_param);

  // load gauge after setting precision
   _loadGaugeQuda(compression);

  // choose dslash type
  if( g_c_sw > 0.0 ) {
    inv_param.dslash_type = QUDA_TWISTED_CLOVER_DSLASH;
    inv_param.matpc_type = QUDA_MATPC_EVEN_EVEN; // FIXME: note sure if this is the correct PC type
    inv_param.solution_type = QUDA_MAT_SOLUTION;
    inv_param.clover_order = QUDA_PACKED_CLOVER_ORDER;
    inv_param.clover_coeff = g_c_sw*g_kappa;
    inv_param.compute_clover = 1;
    inv_param.compute_clover_inverse = 1;
  }
  else {
    inv_param.dslash_type = QUDA_TWISTED_MASS_DSLASH;
    inv_param.matpc_type = QUDA_MATPC_EVEN_EVEN_ASYMMETRIC;
    inv_param.solution_type = QUDA_MAT_SOLUTION;
  }

  // choose solver
  if(solver_flag == BICGSTAB) {
    if(g_proc_id == 0) {printf("# TM_QUDA: Using BiCGstab!\n"); fflush(stdout);}
    inv_param.inv_type = QUDA_BICGSTAB_INVERTER;
  }
  else {
    /* Here we invert the hermitean operator squared */
    inv_param.inv_type = QUDA_CG_INVERTER;
    if(g_proc_id == 0) {
      printf("# TM_QUDA: Using mixed precision CG!\n");
      printf("# TM_QUDA: mu = %.12f, kappa = %.12f\n", g_mu/2./g_kappa, g_kappa);
      fflush(stdout);
    }
  }

  if( even_odd_flag ) {
    inv_param.solve_type = QUDA_NORMERR_PC_SOLVE;
    if(g_proc_id == 0) printf("# TM_QUDA: Using EO preconditioning!\n");
  }
  else {
    inv_param.solve_type = QUDA_NORMERR_SOLVE;
    if(g_proc_id == 0) printf("# TM_QUDA: Not using EO preconditioning!\n");
  }

  inv_param.tol = sqrt(precision);
  inv_param.maxiter = max_iter;

  if( g_c_sw > 0.0 ){
    _loadCloverQuda(&inv_param);
  }

  // reorder spinor
  reorder_spinor_toQuda( (double*)spinorIn,   inv_param.cpu_prec, 1 );

  // perform the inversion
  tm_stopwatch_push(&g_timers, "invertQuda", "");
  invertQuda(spinorOut, spinorIn, &inv_param);
  tm_stopwatch_pop(&g_timers, 0, 0, "TM_QUDA");

  if( inv_param.verbosity > QUDA_SILENT )
    if(g_proc_id == 0)
      printf("# TM_QUDA: Done: %i iter / %g secs = %g Gflops\n",
             inv_param.iter, inv_param.secs, inv_param.gflops/inv_param.secs);

  // number of CG iterations
  int iteration = inv_param.iter;

  reorder_spinor_fromQuda( (double*)spinorOut,   inv_param.cpu_prec, 1 );
  convert_lexic_to_eo(Even_new_s, Odd_new_s, solver_field[1]);
  convert_lexic_to_eo(Even_new_c, Odd_new_c, solver_field[1]+VOLUME);

  finalize_solver(solver_field, nr_sf);

  tm_stopwatch_pop(&g_timers, 0, 0, "TM_QUDA");

  if(iteration >= max_iter)
    return(-1);

  return(iteration);
}

// if even_odd_flag set
void M_full_quda(spinor * const Even_new, spinor * const Odd_new,  spinor * const Even, spinor * const Odd) {
  inv_param.kappa = g_kappa;
  // IMPORTANT: use opposite TM flavor since gamma5 -> -gamma5 (until LXLYLZT prob. resolved)
  inv_param.mu = -g_mu;
  inv_param.epsilon = 0.0;

  inv_param.twist_flavor = QUDA_TWIST_SINGLET;
  inv_param.Ls = (inv_param.twist_flavor == QUDA_TWIST_NONDEG_DOUBLET ) ? 2 : 1;

  void *spinorIn  = (void*)g_spinor_field[DUM_DERI];   // source
  void *spinorOut = (void*)g_spinor_field[DUM_DERI+1]; // solution

  // reorder spinor
  convert_eo_to_lexic( spinorIn, Even, Odd );
  reorder_spinor_toQuda( (double*)spinorIn, inv_param.cpu_prec, 0 );

  // multiply
  inv_param.solution_type = QUDA_MAT_SOLUTION;
  MatQuda( spinorOut, spinorIn, &inv_param);

  // reorder spinor
  reorder_spinor_fromQuda( (double*)spinorOut, inv_param.cpu_prec, 0 );
  convert_lexic_to_eo( Even_new, Odd_new, spinorOut );
}

// no even-odd
void D_psi_quda(spinor * const P, spinor * const Q) {
  inv_param.kappa = g_kappa;
  // IMPORTANT: use opposite TM flavor since gamma5 -> -gamma5 (until LXLYLZT prob. resolved)
  inv_param.mu = -g_mu;
  inv_param.epsilon = 0.0;

  inv_param.twist_flavor = QUDA_TWIST_SINGLET;
  inv_param.Ls = (inv_param.twist_flavor == QUDA_TWIST_NONDEG_DOUBLET ) ? 2 : 1;

  void *spinorIn  = (void*)Q;
  void *spinorOut = (void*)P;

  // reorder spinor
  reorder_spinor_toQuda( (double*)spinorIn, inv_param.cpu_prec, 0 );

  // multiply
  inv_param.solution_type = QUDA_MAT_SOLUTION;
  MatQuda( spinorOut, spinorIn, &inv_param);

  // reorder spinor
  reorder_spinor_fromQuda( (double*)spinorIn,  inv_param.cpu_prec, 0 );
  reorder_spinor_fromQuda( (double*)spinorOut, inv_param.cpu_prec, 0 );
}

// even-odd
void M_quda(spinor * const P, spinor * const Q) {

  _initQuda();

  inv_param.kappa = g_kappa;

  inv_param.twist_flavor = QUDA_TWIST_SINGLET;
  inv_param.preserve_source = QUDA_PRESERVE_SOURCE_YES;
  inv_param.Ls = (inv_param.twist_flavor == QUDA_TWIST_NONDEG_DOUBLET ) ? 2 : 1;
  //inv_param.Ls=1;
  //custom
  _setOneFlavourSolverParam(g_kappa,
                            g_c_sw,
                            g_mu,
                            1,//solver flag
                            1,//even_odd
                            1e-12,
                            1000,
                            1, 0);
  
  inv_param.solution_type = QUDA_MATPC_SOLUTION; 
  inv_param.solve_type = QUDA_NORMOP_PC_SOLVE;
  inv_param.matpc_type = QUDA_MATPC_ODD_ODD_ASYMMETRIC;
  inv_param.dagger = QUDA_DAG_NO;
  inv_param.gamma_basis = QUDA_DEGRAND_ROSSI_GAMMA_BASIS;

  void *spinorIn  = (void*)Q;
  void *spinorOut = (void*)P;

  // reorder spinor
  reorder_spinor_eo_toQuda( (double*)spinorIn, inv_param.cpu_prec, 0 ,1);

  // multiply
  MatQuda(spinorOut, spinorIn, &inv_param);

  // reorder spinor
  reorder_spinor_eo_fromQuda( (double*)spinorIn,  inv_param.cpu_prec, 0, 1);
  reorder_spinor_eo_fromQuda( (double*)spinorOut, inv_param.cpu_prec, 0, 1);
}


void _setOneFlavourSolverParam(const double kappa, const double c_sw, const double mu, 
                               const int solver_type, const int even_odd,
                               const double eps_sq, const int maxiter,
                               const int single_parity_solve,
                               const int QpQm) {

  inv_param.tol = sqrt(eps_sq);
  inv_param.maxiter = maxiter;
  inv_param.Ls = 1;
#ifdef TM_QUDA_EXPERIMENTAL
  inv_param.tm_rho = 0.0;
#endif
  inv_param.clover_rho = 0.0;
  
  // chiral by default, for single-parity, will switch to DEGRAND_ROSSI
  inv_param.gamma_basis = QUDA_CHIRAL_GAMMA_BASIS;

  // choose dslash type
  if( fabs(mu) > 0.0 && c_sw > 0.0 ) {
    inv_param.twist_flavor = QUDA_TWIST_SINGLET;
    inv_param.dslash_type = QUDA_TWISTED_CLOVER_DSLASH;
    inv_param.matpc_type = QUDA_MATPC_EVEN_EVEN;
    inv_param.solution_type = QUDA_MAT_SOLUTION;
    inv_param.clover_order = QUDA_PACKED_CLOVER_ORDER;
    // IMPORTANT: use opposite TM flavor since gamma5 -> -gamma5 (until LXLYLZT prob. resolved)
    inv_param.mu = -mu/2./kappa;
    inv_param.clover_coeff = c_sw*kappa;
    if( fabs(g_mu3) > 2*DBL_EPSILON ){
#ifdef TM_QUDA_EXPERIMENTAL
      inv_param.tm_rho = -g_mu3/2./kappa;
#else
      fatal_error("Attempt to set inv_param.tm_rho but --enable-quda_experimental was not set!", __func__);
#endif
    }

    inv_param.compute_clover_inverse = 1;
    inv_param.compute_clover = 1;
  }
  else if( fabs(mu) > 0.0 ) {
    inv_param.clover_coeff = 0.0;

    inv_param.twist_flavor = QUDA_TWIST_SINGLET;
    inv_param.dslash_type = QUDA_TWISTED_MASS_DSLASH;
    inv_param.matpc_type = QUDA_MATPC_EVEN_EVEN_ASYMMETRIC;
    inv_param.solution_type = QUDA_MAT_SOLUTION;
    // IMPORTANT: use opposite TM flavor since gamma5 -> -gamma5 (until LXLYLZT prob. resolved)
    inv_param.mu = -mu/2./kappa;
  }
  else if( c_sw > 0.0 ) {
    // when we are dealing with the 'rho' mass preconditioning parameter
    // the way this is implemented in QUDA is not consistent with the
    // way it is implemented in tmLQCD
    // To get agreement, we use the twisted clover operator also
    // for Wilson clover fermions in this case.
    if( fabs(g_mu3) > 2*DBL_EPSILON ){
      inv_param.twist_flavor = QUDA_TWIST_SINGLET;
      inv_param.dslash_type = QUDA_TWISTED_CLOVER_DSLASH;
      inv_param.mu = 0.0;
#ifdef TM_QUDA_EXPERIMENTAL
      inv_param.tm_rho = -g_mu3/2./kappa;
#else
    fatal_error("Attempt to set inv_param.tm_rho but --enable-quda_experimental was not set!", __func__);
#endif
    } else {
      inv_param.twist_flavor = QUDA_TWIST_NO;
      inv_param.dslash_type = QUDA_CLOVER_WILSON_DSLASH;
    }
    inv_param.matpc_type = QUDA_MATPC_EVEN_EVEN;
    inv_param.solution_type = QUDA_MAT_SOLUTION;
    inv_param.clover_order = QUDA_PACKED_CLOVER_ORDER;
    inv_param.clover_coeff = c_sw*kappa;
    inv_param.compute_clover_inverse = 1;
    inv_param.compute_clover = 1;
  }
  else {
    inv_param.mu = 0.0;
    inv_param.clover_coeff = 0.0;
    inv_param.twist_flavor = QUDA_TWIST_NO;
    inv_param.dslash_type = QUDA_WILSON_DSLASH;
    if( single_parity_solve ){
      // for single parity solves, we employ QUDA_MATPC_ODD_ODD_ASYMMETRIC below
      // which is not supported for the plain Wilson operator
      // so we work around this by using the twisted mass dslash with zero
      // twisted mass
      inv_param.dslash_type = QUDA_TWISTED_MASS_DSLASH;
      inv_param.twist_flavor = QUDA_TWIST_SINGLET;
    }
    inv_param.matpc_type = QUDA_MATPC_EVEN_EVEN;
    inv_param.solution_type = QUDA_MAT_SOLUTION;
  }
  
  // set up for single parity solves (such as those in the HMC)
  if( single_parity_solve ){
    // when doing single parity, we change to the DeGrand-Rossi basis
    inv_param.gamma_basis = QUDA_DEGRAND_ROSSI_GAMMA_BASIS;
    // and we always want the solution to the asymetrically preconditioned 
    // problem
    inv_param.matpc_type = QUDA_MATPC_ODD_ODD_ASYMMETRIC;

    if( QpQm ){
      if( solver_type == MG || solver_type == BICGSTAB ){
        // when we use MG or BiCGstab to solve the QpQm^{-1} problem, we will
        // perform two solves
        inv_param.solution_type = QUDA_MATPC_SOLUTION;
      } else {
        // otherwise we can get the solution in a single solve
        inv_param.solution_type = QUDA_MATPCDAG_MATPC_SOLUTION;
      }
    } else {
      // when not wanting the inverse of the squred operator explicitly, 
      // we always just need the MATPC solution
      inv_param.solution_type = QUDA_MATPC_SOLUTION;
    }
  }

  // choose solver
  if( solver_type == BICGSTAB ) {
    if(g_proc_id == 0) {printf("# TM_QUDA: Using BiCGstab!\n"); fflush(stdout);}
    inv_param.inv_type = QUDA_BICGSTAB_INVERTER;
  } else if ( solver_type == MG ) {
    if(g_proc_id == 0) {printf("# TM_QUDA: Using MG!\n"); fflush(stdout);}
    inv_param.inv_type = QUDA_GCR_INVERTER;
    inv_param.gcrNkrylov = quda_input.gcrNkrylov;
    inv_param.inv_type_precondition = QUDA_MG_INVERTER;
    inv_param.schwarz_type = QUDA_ADDITIVE_SCHWARZ;
    inv_param.precondition_cycle = 1;
    inv_param.tol_precondition = 1e-1;
    inv_param.maxiter_precondition = 1;
    // this under/overrelaxation parameter is not related to the ones
    // used in the MG 
    inv_param.omega = 1.0;
  } else {
    /* Here we invert the hermitean operator squared */
    inv_param.inv_type = QUDA_CG_INVERTER;
    if(g_proc_id == 0) {
      printf("# TM_QUDA: Using mixed precision CG!\n");
      fflush(stdout);
    }
  }

  // make sure to reset the preconditioner if we've switched from MG to another
  // solver (for example in the HMC or when doing light and heavy inversions)
  if( solver_type != MG ){
    inv_param.inv_type_precondition = QUDA_INVALID_INVERTER;
    inv_param.preconditioner = NULL;
  }

  // direct or norm-op. solve
  if( inv_param.inv_type == QUDA_CG_INVERTER ) {
    if( even_odd || single_parity_solve ) {
      inv_param.solve_type = QUDA_NORMOP_PC_SOLVE;
      if(g_proc_id == 0) printf("# TM_QUDA: Using EO preconditioning!\n");
    }
    else {
      inv_param.solve_type = QUDA_NORMOP_SOLVE;
      if(g_proc_id == 0) printf("# TM_QUDA: Not using EO preconditioning!\n");
    }
  }
  else {
    if( even_odd || single_parity_solve ) {
      inv_param.solve_type = QUDA_DIRECT_PC_SOLVE;
      if(g_proc_id == 0) printf("# TM_QUDA: Using EO preconditioning!\n");
    }
    else {
      inv_param.solve_type = QUDA_DIRECT_SOLVE;
      if(g_proc_id == 0) printf("# TM_QUDA: Not using EO preconditioning!\n");
    }
  }

  // load clover field if required, doing so in this odd place because we need
  // basic stuff to be set in inv_param
  if( c_sw > 0.0 ) {
    _loadCloverQuda(&inv_param);
  }

  if( g_proc_id == 0){
    printf("# TM_QUDA: mu = %.12f, kappa = %.12f, csw = %.12f\n", mu/2./kappa, kappa, c_sw);
  }
  if(g_proc_id == 0 && g_debug_level > 3){
    printf("------------- OUTER SOLVER InvertParam --------------\n");
    printQudaInvertParam(&inv_param);
    printf("----------------------------------------\n");
  }

    // run / refresh / update the MG setup if required
  if( inv_param.inv_type_precondition == QUDA_MG_INVERTER ){
    tm_debug_printf(0,0,"# TM_QUDA: using MG solver to invert operator with 2kappamu = %.12f\n",
                        -g_kappa*2.0*inv_param.mu);
  
    quda_mg_param.invert_param = &mg_inv_param;
    _setQudaMultigridParam(&quda_mg_param);
    _updateQudaMultigridPreconditioner();
  }
}

void _updateQudaMultigridPreconditioner(){

  if( check_quda_mg_setup_state(&quda_mg_setup_state, &quda_gauge_state, &quda_input) == TM_QUDA_MG_SETUP_RESET ){

    tm_stopwatch_push(&g_timers, "MG_Preconditioner_Setup", "");

    if( quda_mg_preconditioner != NULL ){
      tm_debug_printf(0,0,"# TM_QUDA: Destroying MG Preconditioner Setup\n");
      destroyMultigridQuda(quda_mg_preconditioner);
      reset_quda_mg_setup_state(&quda_mg_setup_state);
      quda_mg_preconditioner = NULL;
    }
    tm_debug_printf(0,0,"# TM_QUDA: Performing MG Preconditioner Setup for gauge_id: %f\n", quda_gauge_state.gauge_id);

    int update_required_after_setup = 0;

    // if this call originates from a monomial with tm_rho set, we put this to zero
    // for generting the setup
#ifdef TM_QUDA_EXPERIMENTAL
    double save_tm_rho = 0.0;
    if( fabs(quda_mg_param.invert_param->tm_rho) > 2*DBL_EPSILON ){
      save_tm_rho = quda_mg_param.invert_param->tm_rho;
      quda_mg_param.invert_param->tm_rho = 0.0;
      update_required_after_setup = 1;
    }
#endif

    // if we have set an explicit mu value for the generation of our MG setup,
    // we would like to use it here (but only if it really differs from the
    // target mu)
    double save_mu = 0.0;
    if( fabs(quda_input.mg_setup_2kappamu) > 2*DBL_EPSILON && 
        ( fabs(quda_input.mg_setup_2kappamu) - 2*g_kappa*fabs(quda_mg_param.invert_param->mu) ) > 2*DBL_EPSILON ){

      save_mu = quda_mg_param.invert_param->mu;
      // note the minus sign
      quda_mg_param.invert_param->mu = -quda_input.mg_setup_2kappamu/2.0/g_kappa;
      tm_debug_printf(0,0,"# TM_QUDA: Generating MG Setup with mu = %.12f instead of %.12f\n", 
                          -quda_mg_param.invert_param->mu,
                          -save_mu);

      update_required_after_setup = 1;
    }
   
    // perform the actual setup 
    quda_mg_preconditioner = newMultigridQuda(&quda_mg_param);

#ifdef TM_QUDA_EXPERIMENTAL
    if( fabs(save_tm_rho) > 2*DBL_EPSILON ){
      quda_mg_param.invert_param->tm_rho = save_tm_rho;
    }
#endif
    if( fabs(save_mu) > 2*DBL_EPSILON ){
      quda_mg_param.invert_param->mu = save_mu;
    }
   
    // if either mu or tm_rho were changed during setup, we run an update now 
    if( update_required_after_setup ){
      updateMultigridQuda(quda_mg_preconditioner, &quda_mg_param);
    }

    inv_param.preconditioner = quda_mg_preconditioner;
    set_quda_mg_setup_state(&quda_mg_setup_state, &quda_gauge_state);

    tm_stopwatch_pop(&g_timers, 0, 1, "TM_QUDA");

  } else if ( check_quda_mg_setup_state(&quda_mg_setup_state, &quda_gauge_state, &quda_input) == TM_QUDA_MG_SETUP_REFRESH ) {

    tm_stopwatch_push(&g_timers, "MG_Preconditioner_Setup_Refresh", "");
    tm_debug_printf(0,0,"# TM_QUDA: Refreshing MG Preconditioner Setup for gauge_id: %f\n", quda_gauge_state.gauge_id);
  
    int update_required_after_refresh = 0;

    // the same spiel as we had for mu and tm_rho in the setup phase above
#ifdef TM_QUDA_EXPERIMENTAL
    double save_tm_rho = 0.0;
    if( fabs(quda_mg_param.invert_param->tm_rho) > 2*DBL_EPSILON ){
      save_tm_rho = quda_mg_param.invert_param->tm_rho;
      quda_mg_param.invert_param->tm_rho = 0.0;
      update_required_after_refresh = 1;
    }
#endif // TM_QUDA_EXPERIMENTAL
    
    double save_mu = 0.0;
    if( fabs(quda_input.mg_setup_2kappamu) > 2*DBL_EPSILON && 
        ( fabs(quda_input.mg_setup_2kappamu) - 2*g_kappa*fabs(quda_mg_param.invert_param->mu) ) > 2*DBL_EPSILON ){

      save_mu = quda_mg_param.invert_param->mu;
      // note the minus sign
      quda_mg_param.invert_param->mu = -quda_input.mg_setup_2kappamu/2.0/g_kappa;
      tm_debug_printf(0,0,"# TM_QUDA: Refreshing MG Setup with mu = %.12f instead of %.12f\n", 
                          -quda_mg_param.invert_param->mu,
                          -save_mu);
      update_required_after_refresh = 1;
    }

    for(int level = 0; level < (quda_input.mg_n_level-1); level++){
      quda_mg_param.setup_maxiter_refresh[level] = quda_input.mg_setup_maxiter_refresh[level];
    }

    // update the parameters AND refresh the setup
    updateMultigridQuda(quda_mg_preconditioner, &quda_mg_param);
    // reset refresh iterations to zero such that the next call
    // to updateMultigridQuda only updates parameters and coarse
    // operator(s) (unless another refresh is due)
    for(int level = 0; level < (quda_input.mg_n_level-1); level++){
      quda_mg_param.setup_maxiter_refresh[level] = 0;
    }

#ifdef TM_QUDA_EXPERIMENTAL
    if( fabs(save_tm_rho) > 2*DBL_EPSILON ){
      quda_mg_param.invert_param->tm_rho = save_tm_rho;
    }
#endif // TM_QUDA_EXPERIMENTAL
    if( fabs(save_mu) > 2*DBL_EPSILON ){
      quda_mg_param.invert_param->mu = save_mu;
    }

    // if we used parameters different from the target parameters for the refresh.
    // then we need to run the update again (without the refresh part)
    if( update_required_after_refresh ){
      updateMultigridQuda(quda_mg_preconditioner, &quda_mg_param);
    }

    inv_param.preconditioner = quda_mg_preconditioner;
    set_quda_mg_setup_state(&quda_mg_setup_state, &quda_gauge_state);

    tm_stopwatch_pop(&g_timers, 0, 1, "TM_QUDA");

  } else if ( check_quda_mg_setup_state(&quda_mg_setup_state, &quda_gauge_state, &quda_input) == TM_QUDA_MG_SETUP_UPDATE )  {

    tm_stopwatch_push(&g_timers, "MG_Preconditioner_Setup_Update", "");

    tm_debug_printf(0,0,"# TM_QUDA: Updating MG Preconditioner Setup for gauge_id: %f\n", quda_gauge_state.gauge_id);
#ifdef TM_QUDA_EXPERIMENTAL
    if( quda_input.mg_eig_preserve_deflation == QUDA_BOOLEAN_YES ){
      tm_debug_printf(0,0,"# TM_QUDA: Deflation subspace for gauge_id: %f will be re-used!\n", quda_gauge_state.gauge_id);
    }
#endif

    updateMultigridQuda(quda_mg_preconditioner, &quda_mg_param);

    // if the precondioner was disabled because we switched solvers from MG to some other
    // solver, re-enable it here
    inv_param.preconditioner = quda_mg_preconditioner;

    tm_stopwatch_pop(&g_timers, 0, 1, "TM_QUDA");

  } else {
    // if the precondioner was disabled because we switched solvers from MG to some other
    // solver, re-enable it here
    inv_param.preconditioner = quda_mg_preconditioner;
    tm_debug_printf(0,0,"# TM_QUDA: Reusing MG Preconditioner Setup for gauge_id: %f\n", quda_gauge_state.gauge_id);
  }

  if(g_proc_id == 0 && g_debug_level > 3 && inv_param.inv_type_precondition == QUDA_MG_INVERTER){
    printf("--------------- MG InvertParam ------------------\n");
    printQudaInvertParam(quda_mg_param.invert_param);
    printf("---------------- MG MultigridParam ------------------------\n");
    printQudaMultigridParam(&quda_mg_param);
    printf("----------------------------------------\n");
  }
}

void _setTwoFlavourSolverParam(const double kappa, const double c_sw, const double mu,
                               const double epsilon, const int solver_type, const int even_odd,
                               const double eps_sq, const int maxiter,
                               const int single_parity_solve,
                               const int QpQm) {

  inv_param.tol = sqrt(eps_sq);
  inv_param.maxiter = maxiter;
  inv_param.Ls = 2;
  
  inv_param.twist_flavor = QUDA_TWIST_NONDEG_DOUBLET;

#ifdef TM_QUDA_EXPERIMENTAL
  inv_param.tm_rho = 0.0;
#endif
  inv_param.clover_rho = 0.0;

  // choose dslash type
  if( c_sw > DBL_EPSILON ) {
    inv_param.dslash_type = QUDA_TWISTED_CLOVER_DSLASH;
    inv_param.matpc_type = QUDA_MATPC_EVEN_EVEN;
    inv_param.solution_type = QUDA_MAT_SOLUTION;
    inv_param.clover_order = QUDA_PACKED_CLOVER_ORDER;
    // IMPORTANT: use opposite TM flavor since gamma5 -> -gamma5 (until LXLYLZT prob. resolved)
    inv_param.mu = -mu/2./kappa;
    inv_param.epsilon = epsilon/2./kappa;
    inv_param.clover_coeff = c_sw*kappa;
    inv_param.compute_clover_inverse = 1;
    inv_param.compute_clover = 1;
  } else {
    inv_param.dslash_type = QUDA_TWISTED_MASS_DSLASH;
    inv_param.matpc_type = QUDA_MATPC_EVEN_EVEN_ASYMMETRIC;
    inv_param.solution_type = QUDA_MAT_SOLUTION;
    
    inv_param.clover_coeff = 0.0;

    // IMPORTANT: use opposite TM flavor since gamma5 -> -gamma5 (until LXLYLZT prob. resolved)
    inv_param.mu = -mu/2./kappa;
    inv_param.epsilon = epsilon/2./kappa; 
    inv_param.compute_clover_inverse = 0;
    inv_param.compute_clover = 0;
  }
  
  // set up for single parity solves (such as those in the HMC)
  if( single_parity_solve ){
    // when doing single parity, we change to the DeGrand-Rossi basis
    inv_param.gamma_basis = QUDA_DEGRAND_ROSSI_GAMMA_BASIS;
    // and we always want the solution to the asymetrically preconditioned 
    // problem
    inv_param.matpc_type = QUDA_MATPC_ODD_ODD_ASYMMETRIC;

    if( QpQm ){
      if( solver_type == MG || solver_type == BICGSTAB ){
        // when we use MG or BiCGstab to solve the QpQm^{-1} problem, we will
        // perform two solves
        inv_param.solution_type = QUDA_MATPC_SOLUTION;
      } else {
        // otherwise we can get the solution in a single solve
        inv_param.solution_type = QUDA_MATPCDAG_MATPC_SOLUTION;
      }
    } else {
      // when not wanting the inverse of the squred operator explicitly, 
      // we always just need the MATPC solution
      inv_param.solution_type = QUDA_MATPC_SOLUTION;
    }
  }

  // choose solver
  if( solver_type == BICGSTAB ) {
    if(g_proc_id == 0) {printf("# TM_QUDA: Using BiCGstab!\n"); fflush(stdout);}
    inv_param.inv_type = QUDA_BICGSTAB_INVERTER;
  } else if ( solver_type == MG ) {
    fatal_error("MG unsupported for non-degenerate inversions", "setTwoFlavourSolverParam");
    if(g_proc_id == 0) {printf("# TM_QUDA: Using MG!\n"); fflush(stdout);}
    inv_param.inv_type = QUDA_GCR_INVERTER;
    inv_param.gcrNkrylov = quda_input.gcrNkrylov;
    inv_param.inv_type_precondition = QUDA_MG_INVERTER;
    inv_param.schwarz_type = QUDA_ADDITIVE_SCHWARZ;
    inv_param.precondition_cycle = 1;
    inv_param.tol_precondition = 1e-1;
    inv_param.maxiter_precondition = 1;
    // this under/overrelaxation parameter is not related to the ones
    // used in the MG 
    inv_param.omega = 1.0;
  } else {
    /* Here we invert the hermitean operator squared */
    inv_param.inv_type = QUDA_CG_INVERTER;
    if(g_proc_id == 0) {
      printf("# TM_QUDA: Using mixed precision CG!\n");
      fflush(stdout);
    }
  }

  // make sure to reset the preconditioner if we've switched from MG to another
  // solver (for example in the HMC or when doing light and heavy inversions)
  if( solver_type != MG ){
    inv_param.inv_type_precondition = QUDA_INVALID_INVERTER;
  }

  // direct or norm-op. solve
  if( inv_param.inv_type == QUDA_CG_INVERTER ) {
    if( even_odd || single_parity_solve ) {
      inv_param.solve_type = QUDA_NORMOP_PC_SOLVE;
      if(g_proc_id == 0) printf("# TM_QUDA: Using EO preconditioning!\n");
    }
    else {
      inv_param.solve_type = QUDA_NORMOP_SOLVE;
      if(g_proc_id == 0) printf("# TM_QUDA: Not using EO preconditioning!\n");
    }
  }
  else {
    if( even_odd || single_parity_solve ) {
      inv_param.solve_type = QUDA_DIRECT_PC_SOLVE;
      if(g_proc_id == 0) printf("# TM_QUDA: Using EO preconditioning!\n");
    }
    else {
      inv_param.solve_type = QUDA_DIRECT_SOLVE;
      if(g_proc_id == 0) printf("# TM_QUDA: Not using EO preconditioning!\n");
    }
  }

  // load clover field if required, doing so in this odd place because we need
  // basic stuff to be set in inv_param
  if( c_sw > 0.0 ) {
    _loadCloverQuda(&inv_param);
  }

  if( g_proc_id == 0){
    printf("# TM_QUDA: mu = %.12f, epsilon = %.12f kappa = %.12f, csw = %.12f\n", 
           mu/2./kappa, epsilon/2./kappa, kappa, c_sw);
  }
  if(g_proc_id == 0 && g_debug_level > 3){
    printf("------------- OUTER SOLVER InvertParam --------------\n");
    printQudaInvertParam(&inv_param);
    printf("----------------------------------------\n");
  }

  // run the MG setup if required
  // FIXME: even if this is supported at some point, it will require some work to make
  // sure that we can keep to MG preconditioners around (or perhaps the same one can
  // be used for both operators?)
  if( inv_param.inv_type_precondition == QUDA_MG_INVERTER ){
    tm_debug_printf(0,0,"# TM_QUDA: using MG solver to invert operator with 2kappamu = %.12f 2kappaeps = %.12f\n",
                        mu, epsilon);
    quda_mg_param.invert_param = &mg_inv_param;
    _setQudaMultigridParam(&quda_mg_param);
    _updateQudaMultigridPreconditioner();
  }
}

_setMGInvertParam(QudaInvertParam * mg_inv_param, const QudaInvertParam * const inv_param){
  // reset the mg_inv_param to start from a clean slate
  (*mg_inv_param) = newQudaInvertParam();

  mg_inv_param->Ls = 1;
  mg_inv_param->sp_pad = 0;
  mg_inv_param->cl_pad = 0;

  mg_inv_param->residual_type = QUDA_L2_RELATIVE_RESIDUAL;

  mg_inv_param->preserve_source = QUDA_PRESERVE_SOURCE_YES;
  // the MG internal Gamma basis is always DEGRAND_ROSSI
  mg_inv_param->gamma_basis = QUDA_DEGRAND_ROSSI_GAMMA_BASIS;
  mg_inv_param->dirac_order = QUDA_DIRAC_ORDER;
  
  mg_inv_param->input_location = QUDA_CPU_FIELD_LOCATION;
  mg_inv_param->output_location = QUDA_CPU_FIELD_LOCATION;
  
  // currently, only QUDA_DIRECT_SOLVE is supported for this, thus also QUDA_MAT_SOLUTION 
  mg_inv_param->solve_type = QUDA_DIRECT_SOLVE;
  mg_inv_param->solution_type = QUDA_MAT_SOLUTION;

  mg_inv_param->dagger = QUDA_DAG_NO;

  mg_inv_param->verbosity = QUDA_SUMMARIZE;
  if( g_debug_level >= 3 ){
    mg_inv_param->verbosity = QUDA_VERBOSE;
  }
  mg_inv_param->verbosity_precondition = QUDA_SUMMARIZE;

  // now copy over relevant stuff from the outer solver
  
  // only the symmetric operator may be coarsened, so this
  // is what we use in the MG internal InvertParam
  //
  // when the outer solve is not matpc, we still do matpc internally
  if( inv_param->matpc_type == QUDA_MATPC_ODD_ODD_ASYMMETRIC || inv_param->matpc_type == QUDA_MATPC_ODD_ODD ){
    mg_inv_param->matpc_type = QUDA_MATPC_ODD_ODD;
  } else if( inv_param->matpc_type == QUDA_MATPC_EVEN_EVEN_ASYMMETRIC || inv_param->matpc_type == QUDA_MATPC_EVEN_EVEN ||
             inv_param->matpc_type == QUDA_MATPC_INVALID ) {
    mg_inv_param->matpc_type = QUDA_MATPC_EVEN_EVEN;
  }

  // now set the relevant MG-invernal inverter parameters
  mg_inv_param->inv_type = inv_param->inv_type;
  mg_inv_param->tol = inv_param->tol;
  mg_inv_param->maxiter = inv_param->maxiter;
  mg_inv_param->reliable_delta = inv_param->reliable_delta;
  mg_inv_param->mass_normalization = inv_param->mass_normalization;

  mg_inv_param->cpu_prec = inv_param->cpu_prec;
  mg_inv_param->cuda_prec = inv_param->cuda_prec;
  mg_inv_param->cuda_prec_sloppy = inv_param->cuda_prec_sloppy;
  mg_inv_param->cuda_prec_refinement_sloppy = inv_param->cuda_prec_refinement_sloppy;
  mg_inv_param->cuda_prec_precondition = inv_param->cuda_prec_precondition;
  mg_inv_param->cuda_prec_eigensolver = inv_param->cuda_prec_eigensolver;
  
  mg_inv_param->clover_cpu_prec = inv_param->clover_cpu_prec;
  mg_inv_param->clover_cuda_prec = inv_param->clover_cuda_prec;
  mg_inv_param->clover_cuda_prec_sloppy = inv_param->clover_cuda_prec_sloppy;
  mg_inv_param->clover_cuda_prec_refinement_sloppy = inv_param->clover_cuda_prec_refinement_sloppy;
  mg_inv_param->clover_cuda_prec_precondition = inv_param->clover_cuda_prec_precondition;
  mg_inv_param->clover_cuda_prec_eigensolver = inv_param->clover_cuda_prec_eigensolver;
  
  mg_inv_param->clover_order = inv_param->clover_order;
  mg_inv_param->gcrNkrylov = inv_param->gcrNkrylov;

  mg_inv_param->dslash_type = inv_param->dslash_type;
  mg_inv_param->twist_flavor = inv_param->twist_flavor;
  mg_inv_param->mu = inv_param->mu;
  mg_inv_param->kappa = inv_param->kappa;
  mg_inv_param->clover_coeff = inv_param->clover_coeff;
#ifdef TM_QUDA_EXPERIMENTAL
  mg_inv_param->tm_rho = inv_param->tm_rho;
#endif
}

void _setQudaMultigridParam(QudaMultigridParam* mg_param) {
  QudaInvertParam * mg_inv_param = mg_param->invert_param;
  _setMGInvertParam(mg_inv_param, &inv_param);

#ifdef TM_QUDA_EXPERIMENTAL
  mg_param->preserve_deflation = quda_input.mg_eig_preserve_deflation;
#endif

  mg_param->n_level = quda_input.mg_n_level;
  for (int level=0; level < mg_param->n_level; level++) {
    for (int dim=0; dim<4; dim++) {
      int extent;
      switch(dim){
        case 0:
          extent = LX;
          break;
        case 1:
          extent = LY;
          break;
        case 2:
          extent = LZ;
          break;
        case 3:
        default:
          extent = T;
          break;
      }
      // determine how many lattice sites remain at the current level
      for(int k = level; k > 0; k--) {
        extent = extent/mg_param->geo_block_size[k-1][dim];
      }

      if( level == (quda_input.mg_n_level-1) ){
        // for the coarsest level, the block size is always set to 1
        mg_param->geo_block_size[level][dim] = 1;
      } else if( quda_input.mg_blocksize[level][dim] != 0 ){
        // the block size for this level and dimension has been set non-zero in the input file
        // we respect this no matter what
        mg_param->geo_block_size[level][dim] = quda_input.mg_blocksize[level][dim];
        // otherwise we employ our blocking algorithm
      } else {
        // on all levels, we try to use a block size of 4^4 and compute the
        // number of fine or aggregate lattice sites on a given level,
        // resulting in block sizes of:
        // - 4 if the extent is larger or equal to 16 and
        // - 2 otherwise
        // When an extent is divisible by three, smaller or equal to 24 and when we're
        // not on the finest grid [and the user has explicitly enabled support 
        // for these block lengths  (and therefore also adjusted QUDA to instantiate them)],
        // we use a block length of 3.
        // If aggregation using an even number of lattice points (if size 3 is disabled)
        // is not possible or if the extent is 1 or some divisible only by some prime number
        // other than 3 or 2, we use a block size of 1
        int even_block_size = 4;
        if( extent < 16 ) even_block_size = 2;
     
        // special treatment of size 24 lattice extents on the fine grid
        if ( extent <= 24 && extent % 3 == 0 && quda_input.mg_enable_size_three_blocks ) {
          mg_param->geo_block_size[level][dim] = 3;
        } else if ( extent % even_block_size == 0 ) { 
          mg_param->geo_block_size[level][dim] = even_block_size;
        } else {
          mg_param->geo_block_size[level][dim] = 1;
        }
      }
      
      // this output is only relevant on levels 0, 1, ..., n-2
      if( level < (mg_param->n_level-1) && g_proc_id == 0 && g_debug_level >= 2 ) {
        printf("# TM_QUDA: MG level %d, extent of (xyzt) dim %d: %d\n", level, dim, extent);
        printf("# TM_QUDA: MG aggregation size set to: %d\n", mg_param->geo_block_size[level][dim]);
        fflush(stdout);
      }

      // all lattice extents must be even after blocking on all levels
      if( level < mg_param->n_level-1 && 
          (extent / mg_param->geo_block_size[level][dim]) % 2 != 0 ){
        tm_debug_printf(0, 0,
                        "MG level %d, dim %d (xyzt) has extent %d. Block size of %d would result "
                        "in odd extent on level %d, aborting!\n"
                        "Adjust your block sizes or parallelisation, all local lattice extents on all levels must be even!\n",
                        level, dim, extent, mg_param->geo_block_size[level][dim], level+1);
        fflush(stdout);
        fatal_error("Blocking error.\n", "_setQudaMultigridParam");
      }

    } // for( dim=0 to dim=3 ) (space-time dimensions)
    
    mg_param->verbosity[level] = quda_input.mg_verbosity[level];
    mg_param->precision_null[level] = QUDA_HALF_PRECISION;
    mg_param->setup_inv_type[level] = quda_input.mg_setup_inv_type;
    // Kate says: experimental, leave at 1 (will be used for bootstrap-style setup later)
    mg_param->num_setup_iter[level] = 1;
    mg_param->setup_tol[level] = quda_input.mg_setup_tol[level];
    mg_param->setup_maxiter[level] = quda_input.mg_setup_maxiter[level];
    // If doing twisted mass, we can scale the twisted mass on the coarser grids
    // which significantly increases speed of convergence as a result of making
    // the coarsest grid solve a lot better conditioned.
    // Dean Howarth has some RG arguments on why the coarse mass parameter should be
    // rescaled for the coarse operator to be optimal.
    if( fabs(mg_inv_param->mu) > 2*DBL_EPSILON ) {
      mg_param->mu_factor[level] = quda_input.mg_mu_factor[level];
      if( g_proc_id == 0 && g_debug_level >= 2 ){
        printf("# TM_QUDA: MG setting coarse mu scaling factor on level %d to %lf\n", level, mg_param->mu_factor[level]);
      }
    }
    
    mg_param->coarse_solver[level] = quda_input.mg_coarse_solver_type[level];
    mg_param->coarse_solver_tol[level] = quda_input.mg_coarse_solver_tol[level];
    mg_param->coarse_solver_maxiter[level] = quda_input.mg_coarse_solver_maxiter[level];
    // spin block size on level zero will be reset to 2 below
    mg_param->spin_block_size[level] = 1;
    mg_param->n_vec[level] = quda_input.mg_n_vec[level];
    mg_param->nu_pre[level] = quda_input.mg_nu_pre[level];
    mg_param->nu_post[level] = quda_input.mg_nu_post[level];

    mg_param->cycle_type[level] = QUDA_MG_CYCLE_RECURSIVE;
    mg_param->location[level] = QUDA_CUDA_FIELD_LOCATION;
    mg_param->setup_location[level] = QUDA_CUDA_FIELD_LOCATION;
    
    mg_param->smoother[level] = quda_input.mg_smoother_type[level];
    mg_param->smoother_tol[level] = quda_input.mg_smoother_tol[level];
    // unless the Schwarz-alternating smoother is used, this should be 1
    mg_param->smoother_schwarz_cycle[level] = 1;
    // Kate says this should be EO always for performance
    mg_param->smoother_solve_type[level] = QUDA_DIRECT_PC_SOLVE;
    mg_param->smoother_schwarz_type[level] = QUDA_INVALID_SCHWARZ;
    mg_param->smoother_halo_precision[level] = QUDA_HALF_PRECISION;
   
    // when the Schwarz-alternating smoother is used, this can be set to NO, otherwise it must be YES 
    mg_param->global_reduction[level] = QUDA_BOOLEAN_YES;

    // set to QUDA_MAT_SOLUTION to inject a full field into coarse grid
    // set to QUDA_MATPC_SOLUTION to inject single parity field into coarse grid
    // if we are using an outer even-odd preconditioned solve, then we
    // use single parity injection into the coarse grid
    mg_param->coarse_grid_solution_type[level] = (level == 0 && inv_param.solve_type == QUDA_DIRECT_PC_SOLVE) ? 
                                                    QUDA_MATPC_SOLUTION : 
                                                    QUDA_MAT_SOLUTION;

    mg_param->omega[level] = quda_input.mg_omega[level]; // over/under relaxation factor

    mg_param->location[level] = QUDA_CUDA_FIELD_LOCATION;

    mg_param->setup_ca_basis[level]      = quda_input.mg_setup_ca_basis[level];
    mg_param->setup_ca_basis_size[level] = quda_input.mg_setup_ca_basis_size[level];
    mg_param->setup_ca_lambda_min[level] = quda_input.mg_setup_ca_lambda_min[level];
    mg_param->setup_ca_lambda_max[level] = quda_input.mg_setup_ca_lambda_max[level];

    mg_param->coarse_solver_ca_basis[level]      = quda_input.mg_coarse_solver_ca_basis[level];
    mg_param->coarse_solver_ca_basis_size[level] = quda_input.mg_coarse_solver_ca_basis_size[level];
    mg_param->coarse_solver_ca_lambda_min[level] = quda_input.mg_coarse_solver_ca_lambda_min[level];
    mg_param->coarse_solver_ca_lambda_max[level] = quda_input.mg_coarse_solver_ca_lambda_max[level];
    

    // set the MG EigSolver parameters, almost equivalent to
    // setEigParam from QUDA's multigrid_invert_test, except
    // for cuda_prec_ritz (on 20190822)
    if( quda_input.mg_use_eig_solver[level] == QUDA_BOOLEAN_YES ){
      mg_param->use_eig_solver[level] = QUDA_BOOLEAN_YES;
      mg_eig_param[level].eig_type = quda_input.mg_eig_type[level];
      mg_eig_param[level].spectrum = quda_input.mg_eig_spectrum[level];
      if ((quda_input.mg_eig_type[level] == QUDA_EIG_TR_LANCZOS || 
           quda_input.mg_eig_type[level] == QUDA_EIG_IR_ARNOLDI)
          && !(quda_input.mg_eig_spectrum[level] == QUDA_SPECTRUM_LR_EIG || 
               quda_input.mg_eig_spectrum[level] == QUDA_SPECTRUM_SR_EIG)) {
        tm_debug_printf(0, 0,
                        "ERROR: MG level %d: Only real spectrum type (LR or SR)"
                          "can be passed to the a Lanczos type solver!\n",
                        level);
        fflush(stdout);
        fatal_error("Eigensolver parameter error.\n", "_setQudaMultigridParam");
      }

      mg_eig_param[level].n_ev = quda_input.mg_eig_nEv[level];
      mg_eig_param[level].n_kr = quda_input.mg_eig_nKr[level];
      mg_eig_param[level].n_conv = quda_input.mg_n_vec[level];
      mg_eig_param[level].require_convergence = quda_input.mg_eig_require_convergence[level];

      mg_eig_param[level].tol = quda_input.mg_eig_tol[level];
      mg_eig_param[level].check_interval = quda_input.mg_eig_check_interval[level];
      mg_eig_param[level].max_restarts = quda_input.mg_eig_max_restarts[level];
			// in principle this can be set to a different precision, but we always
      // use double precision in the outer solver
      mg_eig_param[level].cuda_prec_ritz = QUDA_DOUBLE_PRECISION;

      // this seems to be set to NO in multigrid_invert_test
      mg_eig_param[level].compute_svd = QUDA_BOOLEAN_NO;
      mg_eig_param[level].use_norm_op = quda_input.mg_eig_use_normop[level]; 
      mg_eig_param[level].use_dagger = quda_input.mg_eig_use_dagger[level];
      mg_eig_param[level].use_poly_acc = quda_input.mg_eig_use_poly_acc[level]; 
      mg_eig_param[level].poly_deg = quda_input.mg_eig_poly_deg[level];
      mg_eig_param[level].a_min = quda_input.mg_eig_amin[level];
      mg_eig_param[level].a_max = quda_input.mg_eig_amax[level];

      // set file i/o parameters
      // Give empty strings, Multigrid will handle IO.
      strcpy(mg_eig_param[level].vec_infile, "");
      strcpy(mg_eig_param[level].vec_outfile, "");
      strncpy(mg_eig_param[level].QUDA_logfile, "quda_eig.log", 512);

      mg_param->eig_param[level] = &(mg_eig_param[level]);
    } else {
      mg_param->eig_param[level] = NULL;
      mg_param->use_eig_solver[level] = QUDA_BOOLEAN_NO;
    } // end of else branch of if(quda_input.mg_use_eig_solver[level] == QUDA_BOOLEAN_YES) above
  } // for(i=0 to n_level-1)

  // only coarsen the spin on the first restriction
  mg_param->spin_block_size[0] = 2;

  mg_param->compute_null_vector = QUDA_COMPUTE_NULL_VECTOR_YES;
  mg_param->generate_all_levels = QUDA_BOOLEAN_YES;

  mg_param->run_low_mode_check = quda_input.mg_run_low_mode_check;
  mg_param->run_oblique_proj_check = quda_input.mg_run_oblique_proj_check;
  mg_param->run_verify = quda_input.mg_run_verify;

  // set file i/o parameters
  strcpy(mg_param->vec_infile, "");
  strcpy(mg_param->vec_outfile, "");
}

int invert_eo_degenerate_quda(spinor * const out,
                              spinor * const in,
                              const double precision, const int max_iter,
                              const int solver_flag, const int rel_prec,
                              const int even_odd_flag, solver_params_t solver_params,
                              SloppyPrecision sloppy_precision,
                              CompressionType compression,
                              const int QpQm) {
  tm_stopwatch_push(&g_timers, __func__, "");

  int iterations = 0;

  // it returns if quda is already init
  _initQuda();

  void *spinorIn  = (void*)in; // source
  void *spinorOut = (void*)out; // solution

  if ( rel_prec )
    inv_param.residual_type = QUDA_L2_RELATIVE_RESIDUAL;
  else
    inv_param.residual_type = QUDA_L2_ABSOLUTE_RESIDUAL;

  inv_param.kappa = g_kappa;

  // figure out which BC to use (theta, trivial...)
  set_boundary_conditions(&compression, &gauge_param);
  // set the sloppy precision of the mixed prec solver
  set_sloppy_prec(sloppy_precision, &gauge_param, &inv_param);
  
  // load gauge after setting precision
  _loadGaugeQuda(compression);

  // this will also construct the clover field and its inverse, if required
  // it will also run the MG setup
  // internally, it decides about the preconditioning type
  // depending on whether we do a single parity solve
  // and decides further things depending on whether we're inverting
  // ( \hat{Q}^{+} \hat{Q}^{-} ) or just one of \hat{Q}^{+/-} 
  _setOneFlavourSolverParam(g_kappa,
                            g_c_sw,
                            g_mu,
                            solver_flag,
                            even_odd_flag,
                            // when doing a two-step solve, we are more stringent
                            // FIXME: surely this needs to be fine-tuned
                            (solver_flag == MG || solver_flag == BICGSTAB) && QpQm ? precision/10.0 : precision,
                            max_iter,
                            1, QpQm);

  if( solver_flag == MG || solver_flag == BICGSTAB ){
    // for MG and BiCGstab, we solve QpQm in two steps
    // we start with [ \hat{M}^{+} ]^{-1}
    // also in the direct solve of just \hat{Q}^{+/-}, we don't
    // want to invert the daggered operator
    inv_param.dagger = QUDA_DAG_NO; 
    if(solver_flag == MG){
      quda_mg_param.invert_param->gamma_basis = QUDA_DEGRAND_ROSSI_GAMMA_BASIS;
    }
  } else {
    // QUDA applies the MMdag operator, we need QpQm^{-1) in the end
    // so we want QUDA to use the MdagM operator
    inv_param.dagger = QUDA_DAG_YES; 
  }

  reorder_spinor_eo_toQuda( (double*)spinorIn, inv_param.cpu_prec, 0, 1);

  tm_stopwatch_push(&g_timers, "invertQuda", "");
  invertQuda(spinorOut, spinorIn, &inv_param);
  tm_stopwatch_pop(&g_timers, 0, 0, "TM_QUDA");
  
  // the second solve is only necessary in the derivative where we want the inverse of
  // \hat{Q}^{+} \hat{Q}^{-}
  // but we're using solvers that don't operate on the normal system
  if( (solver_flag == MG || solver_flag == BICGSTAB) && QpQm ){
    if(g_proc_id == 0)
      printf("# TM_QUDA: Qp solve done: %i iter / %g secs = %g Gflops\n",
             inv_param.iter, inv_param.secs, inv_param.gflops/inv_param.secs);
    
    iterations += inv_param.iter;

    inv_param.preserve_source = QUDA_PRESERVE_SOURCE_YES;

    reorder_spinor_eo_fromQuda( (double*)spinorOut, inv_param.cpu_prec, 0, 1);
    // we multiply by gamma5 here to obtain the inverse of \hat{Q}^{-}
    // in the next solve
    tm_stopwatch_push(&g_timers, "mul_gamma5", "");
    mul_gamma5((spinor*)spinorOut, VOLUME/2);
    tm_stopwatch_pop(&g_timers, 0, 0, "TM_QUDA");
    reorder_spinor_eo_toQuda( (double*)spinorOut, inv_param.cpu_prec, 0, 1);
   
    // now we invert \hat{M}^{-} to get the inverse of \hat{Q}^{-} in the end
    inv_param.mu = -inv_param.mu;
#ifdef TM_QUDA_EXPERIMENTAL
    inv_param.tm_rho = -inv_param.tm_rho;
#endif
    if(solver_flag == MG){
      // flip the sign of the coarse operator and update the setup
      quda_mg_param.invert_param->mu = -quda_mg_param.invert_param->mu;
      tm_stopwatch_push(&g_timers, "updateMultigridQuda_sign_flip", "");
      updateMultigridQuda(quda_mg_preconditioner, &quda_mg_param);
      // we need to do this to make sure that the MG setup is updated at the next
      // mu flip
      set_quda_mg_setup_mu(&quda_mg_setup_state, -g_mu);
      tm_stopwatch_pop(&g_timers, 0, 1, "TM_QUDA");
    }
    tm_stopwatch_push(&g_timers, "invertQuda", "");
    invertQuda(spinorOut, spinorOut, &inv_param);
    tm_stopwatch_pop(&g_timers, 0, 0, "TM_QUDA");
    reorder_spinor_eo_fromQuda( (double*)spinorOut, inv_param.cpu_prec, 0, 1);
  } else {
    reorder_spinor_eo_fromQuda( (double*)spinorOut, inv_param.cpu_prec, 0, 1);
  }

  if( inv_param.verbosity > QUDA_SILENT )
    if(g_proc_id == 0)
      printf("# TM_QUDA: QpQm solve done: %i iter / %g secs = %g Gflops\n",
             inv_param.iter, inv_param.secs, inv_param.gflops/inv_param.secs);

  iterations += inv_param.iter;

  tm_stopwatch_pop(&g_timers, 0, 0, "TM_QUDA");

  if(iterations >= max_iter)
    return(-1);

  return(iterations);
}

int invert_eo_quda_oneflavour_mshift(spinor ** const out,
                                     spinor * const in,
                                     const double precision, const int max_iter,
                                     const int solver_flag, const int rel_prec,
                                     const int even_odd_flag, solver_params_t solver_params,
                                     SloppyPrecision sloppy_precision,
                                     CompressionType compression){
  tm_stopwatch_push(&g_timers, __func__, "");

  int iterations = 0;

  // it returns if quda is already init
  _initQuda();

  void *spinorIn  = (void*)in; // source
  void **spinorOut = (void**)out; // solution

  if ( rel_prec )
    inv_param.residual_type = QUDA_L2_RELATIVE_RESIDUAL;
  else
    inv_param.residual_type = QUDA_L2_ABSOLUTE_RESIDUAL;

  inv_param.kappa = g_kappa;

  // figure out which BC to use (theta, trivial...)
  set_boundary_conditions(&compression, &gauge_param);
  // set the sloppy precision of the mixed prec solver
  set_sloppy_prec(sloppy_precision, &gauge_param, &inv_param);
  
  // load gauge after setting precision
  _loadGaugeQuda(compression);

  _setOneFlavourSolverParam(g_kappa,
                            g_c_sw,
                            g_mu,
                            solver_flag,
                            even_odd_flag,
                            precision,
                            max_iter,
                            1 /*single_parity_solve */,
                            1 /*always QpQm*/);

  // QUDA applies the MMdag operator, we need QpQm^{-1) in the end
  // so we want QUDA to use the MdagM operator
  inv_param.dagger = QUDA_DAG_YES;

  // just to avoid any issues due to no_shifts being set to zero by accident
  const int num_shifts = solver_params.no_shifts == 0 ? 1 : solver_params.no_shifts;
  inv_param.num_offset = num_shifts;
  for(int shift = 0; shift < num_shifts; shift++){
    inv_param.offset[shift] = solver_params.shifts[shift]*solver_params.shifts[shift];
    inv_param.tol_offset[shift] = sqrt(precision);
  }

  reorder_spinor_eo_toQuda( (double*)spinorIn, inv_param.cpu_prec, 0, 1);

  tm_stopwatch_push(&g_timers, "invertMultiShiftQuda", "");
  invertMultiShiftQuda(spinorOut, spinorIn, &inv_param);
  tm_stopwatch_pop(&g_timers, 0, 0, "TM_QUDA");
  
  for(int shift = 0; shift < num_shifts; shift++){
    reorder_spinor_eo_fromQuda( (double*)spinorOut[shift], inv_param.cpu_prec, 0, 1);
  }

  if( inv_param.verbosity > QUDA_SILENT )
    if(g_proc_id == 0)
      printf("# TM_QUDA: QpQm solve done: %i iter / %g secs = %g Gflops\n",
             inv_param.iter, inv_param.secs, inv_param.gflops/inv_param.secs);

  iterations += inv_param.iter;

  tm_stopwatch_pop(&g_timers, 0, 0, "TM_QUDA");

  if(iterations >= max_iter)
    return(-1);

  return(iterations);
}

int invert_eo_quda_twoflavour_mshift(spinor ** const out_up, spinor ** const out_dn,
                                     spinor * const in_up, spinor * const in_dn,
                                     const double precision, const int max_iter,
                                     const int solver_flag, const int rel_prec,
                                     const int even_odd_flag, solver_params_t solver_params,
                                     SloppyPrecision sloppy_precision,
                                     CompressionType compression){
  tm_stopwatch_push(&g_timers, __func__, "");

  const int Vh = VOLUME/2;

  int iterations = 0;

  const int nr_sf_in = 1;
  spinor ** in;
  spinor ** out;
  init_solver_field(&in, VOLUME, nr_sf_in);

  tm_stopwatch_push(&g_timers, "twoflavour_input_overhead", "");
  memcpy(in[0],      in_up, Vh*24*sizeof(double));
  memcpy(in[0] + Vh, in_dn, Vh*24*sizeof(double));
  tm_stopwatch_pop(&g_timers, 0, 1, "TM_QUDA");

  const int nr_sf_out = solver_params.no_shifts;
  init_solver_field(&out, VOLUME, nr_sf_out);

  // it returns if quda is already init
  _initQuda();

  void *spinorIn  = (void*)in[0]; // source
  void **spinorOut = (void**)out; // solution

  if ( rel_prec )
    inv_param.residual_type = QUDA_L2_RELATIVE_RESIDUAL;
  else
    inv_param.residual_type = QUDA_L2_ABSOLUTE_RESIDUAL;

  inv_param.kappa = g_kappa;

  // figure out which BC to use (theta, trivial...)
  set_boundary_conditions(&compression, &gauge_param);
  // set the sloppy precision of the mixed prec solver
  set_sloppy_prec(sloppy_precision, &gauge_param, &inv_param);
  
  // load gauge after setting precision
  _loadGaugeQuda(compression);

  _setTwoFlavourSolverParam(g_kappa,
                            g_c_sw,
                            g_mubar,
                            g_epsbar,
                            solver_flag,
                            even_odd_flag,
                            precision,
                            max_iter,
                            1 /*single_parity_solve */,
                            1 /*always QpQm*/);

  // QUDA applies the MMdag operator, we need QpQm^{-1) in the end
  // so we want QUDA to use the MdagM operator
  inv_param.dagger = QUDA_DAG_YES;

  // just to avoid any issues due to no_shifts being set to zero by accident
  const int num_shifts = solver_params.no_shifts == 0 ? 1 : solver_params.no_shifts;

  // tmLQCD weights the operator with 1/maxev in the RHMC relative to the shifts 
  // we will do this externally on the shifted inverses (in monomial_solve) and thus need to weight
  // the shifts by maxev^2 to compensate for this
  const double maxev_sq = (1.0/phmc_invmaxev)*(1.0/phmc_invmaxev); 
  inv_param.num_offset = num_shifts;
  for(int shift = 0; shift < num_shifts; shift++){
    inv_param.offset[shift] = maxev_sq*solver_params.shifts[shift]*solver_params.shifts[shift];
    inv_param.tol_offset[shift] = sqrt(precision);
  }

  reorder_spinor_eo_toQuda((double*)spinorIn, inv_param.cpu_prec, 1, 1);

  tm_stopwatch_push(&g_timers, "invertMultiShiftQuda", "");
  invertMultiShiftQuda(spinorOut, spinorIn, &inv_param);
  tm_stopwatch_pop(&g_timers, 0, 0, "TM_QUDA");
 
  for(int shift = 0; shift < num_shifts; shift++){
    reorder_spinor_eo_fromQuda((double*)spinorOut[shift], inv_param.cpu_prec, 1, 1);
    tm_stopwatch_push(&g_timers, "twoflavour_output_overhead", ""); 
    memcpy(out_up[shift], spinorOut[shift],                    24*Vh*sizeof(double));
    memcpy(out_dn[shift], ((double*)spinorOut[shift]) + 24*Vh, 24*Vh*sizeof(double));
    tm_stopwatch_pop(&g_timers, 0, 0, "TM_QUDA");
  }

  finalize_solver(in, nr_sf_in);
  finalize_solver(out, nr_sf_out);

  if( inv_param.verbosity > QUDA_SILENT )
    if(g_proc_id == 0)
      printf("# TM_QUDA: QpQm solve done: %i iter / %g secs = %g Gflops\n",
             inv_param.iter, inv_param.secs, inv_param.gflops/inv_param.secs);

  iterations += inv_param.iter;

  tm_stopwatch_pop(&g_timers, 0, 0, "TM_QUDA");

  if(iterations >= max_iter)
    return(-1);

  return(iterations);
}

void add_mom_to_derivative(su3adj** der) {
  tm_stopwatch_push(&g_timers, __func__, "");

#ifdef TM_USE_OMP
#pragma omp parallel for collapse(2)
#endif
  for( size_t v=0; v<VOLUME; v++ )
    for( int mu=0; mu<4; mu++ ){
      der[v][mu].d1 += mom_quda_reordered[mu][10*v+1]; // imag 01 
      der[v][mu].d2 += mom_quda_reordered[mu][10*v+0]; // real 01
      der[v][mu].d4 += mom_quda_reordered[mu][10*v+3]; // imag 02
      der[v][mu].d5 += mom_quda_reordered[mu][10*v+2]; // real 02
      der[v][mu].d6 += mom_quda_reordered[mu][10*v+5]; // imag 12
      der[v][mu].d7 += mom_quda_reordered[mu][10*v+4]; // real 12 

      double c00 = mom_quda_reordered[mu][10*v+6];
      double c11 = mom_quda_reordered[mu][10*v+7];
      double c22 = mom_quda_reordered[mu][10*v+8];
      der[v][mu].d3 += -(c11-c00)*0.5; // imag 11 - 00
      der[v][mu].d8 += -(2*c22 - c00 - c11) * 0.288675134594813; // imag (2*22 - 00 - 11)/sqrt(3)/2 
                                                                // factor -1/2 required to match tmLQCD results
    }
  
  tm_stopwatch_pop(&g_timers, 0, 0, "TM_QUDA");
}


void compute_gauge_derivative_quda(monomial * const mnl, hamiltonian_field_t * const hf) {
  tm_stopwatch_push(&g_timers, __func__, "");
  
  _initQuda();
  _initMomQuda();

  const int rect = mnl->use_rectangles;

  int * path_length = rect ? plaq_rect_length : plaq_length;

  const int num_paths = rect ? 24 : 6;
  const int max_length = rect ? 5 : 3;

  // can't use stack memory for this so we need to copy to a heap buffer
  int *** path_buf = malloc(4*(num_paths+1)*sizeof(void*));
  for(int i=0; i<4; i++) {
    path_buf[i] = (int**) path_buf+4+i*num_paths;
    for(int j=0; j<num_paths; j++) {
      if(rect) 
        path_buf[i][j] = plaq_rect_path[i][j];
      else
        path_buf[i][j] = plaq_path[i][j];
    }
  }

  double * loop_coeff = malloc(num_paths*sizeof(double));
  for(int i=0; i<num_paths; i++) {
    // Plaq coeffs
    if(i<6)
      loop_coeff[i] = -0.66666666666 * g_beta * mnl->c0; //2/3 conversion factor to match tmLQCD results
    // Rect coeffs
    else
      loop_coeff[i] = -0.66666666666 * g_beta * mnl->c1;
  }
  
  #pragma omp parallel for
  for(int i = 0; i < 4; i++){
    memset(mom_quda[i], 0, VOLUME*10*sizeof(double));
  }

  reorder_gauge_toQuda(hf->gaugefield, NO_COMPRESSION);
  // the reordering above overwrites gauge_quda
  // to make sure that things become synchronised again at the
  // next _loadGaugeQuda, we reset the QUDA gauge state here
  reset_quda_gauge_state(&quda_gauge_state);
  
  tm_stopwatch_push(&g_timers, "computeGaugeForceQuda", ""); 
  computeGaugeForceQuda((void*)mom_quda, 
                        (void*)gauge_quda, 
                        path_buf, path_length, loop_coeff, num_paths, max_length, 1.0,
                        &f_gauge_param);
  tm_stopwatch_pop(&g_timers, 0, 1, "TM_QUDA");

  free(path_buf);
  free(loop_coeff);
  
  reorder_mom_fromQuda(mom_quda);
  add_mom_to_derivative(hf->derivative);

  tm_stopwatch_pop(&g_timers, 0, 1, "TM_QUDA");
}


/***********************************************************************
 *
 * Copyright (C) 2015       Mario Schroeck
 *               2016, 2017 Bartosz Kostrzewa
 *               2018       Bartosz Kostrzewa, Ferenc Pittler
 *               2019       Bartosz Kostrzewa
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
QudaInvertParam inv_mg_param;
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
void _setOneFlavourSolverParam(const double kappa, const double c_sw, const double mu, 
                               const int solver_type, const int even_odd,
                               const double eps_sq, const int maxiter);

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
  // local lattice size
#if USE_LZ_LY_LX_T
  gauge_param.X[0] = LZ;
  gauge_param.X[1] = LY;
  gauge_param.X[2] = LX;
  gauge_param.X[3] = T;
#else
  gauge_param.X[0] = LX;
  gauge_param.X[1] = LY;
  gauge_param.X[2] = LZ;
  gauge_param.X[3] = T;
#endif

  inv_param.Ls = 1;

  gauge_param.anisotropy = 1.0;
  gauge_param.type = QUDA_WILSON_LINKS;
  gauge_param.gauge_order = QUDA_QDP_GAUGE_ORDER;

  gauge_param.cpu_prec = cpu_prec;
  
  gauge_param.cuda_prec = cuda_prec;
  gauge_param.cuda_prec_sloppy = cuda_prec_sloppy;
  gauge_param.cuda_prec_precondition = cuda_prec_precondition;
  gauge_param.cuda_prec_refinement_sloppy = cuda_prec_sloppy;
  gauge_param.cuda_prec_eigensolver = cuda_prec_precondition;
  
  gauge_param.reconstruct = 18;
  gauge_param.reconstruct_sloppy = 18;
  gauge_param.reconstruct_precondition = 18;
  gauge_param.reconstruct_refinement_sloppy = 18;
  gauge_param.reconstruct_eigensolver = 18;
  
  gauge_param.gauge_fix = QUDA_GAUGE_FIXED_NO;

  inv_param.dagger = QUDA_DAG_NO;
  inv_param.mass_normalization = QUDA_KAPPA_NORMALIZATION;
  inv_param.solver_normalization = QUDA_DEFAULT_NORMALIZATION;

  inv_param.pipeline = quda_input.pipeline;
  inv_param.gcrNkrylov = quda_input.gcrNkrylov;

  inv_param.residual_type = (QudaResidualType)(QUDA_L2_RELATIVE_RESIDUAL);
  inv_param.tol_hq = 0.1;
  inv_param.use_alternative_reliable = 0;
  inv_param.reliable_delta = 1e-3; // ignored by multi-shift solver
  inv_param.use_sloppy_partial_accumulator = 0;

  // domain decomposition preconditioner parameters
  inv_param.inv_type_precondition = QUDA_CG_INVERTER;
  inv_param.schwarz_type = QUDA_ADDITIVE_SCHWARZ;
  inv_param.precondition_cycle = 1;
  inv_param.tol_precondition = 1e-1;
  inv_param.maxiter_precondition = 10;
  inv_param.verbosity_precondition = QUDA_SILENT;
  if( g_debug_level >= 5 )
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
  inv_param.clover_cuda_prec_precondition = cuda_prec_precondition;
  inv_param.clover_cuda_prec_refinement_sloppy = cuda_prec_sloppy;
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

  gauge_param.ga_pad = 0; // 24*24*24/2;
  inv_param.sp_pad = 0; // 24*24*24/2;
  inv_param.cl_pad = 0; // 24*24*24/2;

  // For multi-GPU, ga_pad must be large enough to store a time-slice
  int x_face_size = gauge_param.X[1]*gauge_param.X[2]*gauge_param.X[3]/2;
  int y_face_size = gauge_param.X[0]*gauge_param.X[2]*gauge_param.X[3]/2;
  int z_face_size = gauge_param.X[0]*gauge_param.X[1]*gauge_param.X[3]/2;
  int t_face_size = gauge_param.X[0]*gauge_param.X[1]*gauge_param.X[2]/2;
  int pad_size =MAX(x_face_size, y_face_size);
  pad_size = MAX(pad_size, z_face_size);
  pad_size = MAX(pad_size, t_face_size);
  gauge_param.ga_pad = pad_size;

  // solver verbosity
  if( g_debug_level == 0 )
    inv_param.verbosity = QUDA_SILENT;
  else if( g_debug_level >= 1 && g_debug_level < 3 )
    inv_param.verbosity = QUDA_SUMMARIZE;
  else if( g_debug_level >= 3 && g_debug_level < 5 )
    inv_param.verbosity = QUDA_VERBOSE;
  else if( g_debug_level >= 5 )
    inv_param.verbosity = QUDA_DEBUG_VERBOSE;

  // general verbosity
  setVerbosityQuda( QUDA_SUMMARIZE, "# QUDA: ", stdout);
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
  inv_param = newQudaInvertParam();
  inv_mg_param = newQudaInvertParam();
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
    freeGaugeQuda();
    freeCloverQuda(); // this is safe even if there is no Clover field loaded, at least it was in QUDA v0.7.2
    free((void*)tempSpinor);
    endQuda();
  }
}

void _loadCloverQuda(QudaInvertParam* inv_param){
  // check if loaded clover and gauge fields agree
  if( check_quda_clover_state(&quda_clover_state, &quda_gauge_state) ){
    if(g_proc_id==0 && g_debug_level > 0 ) {
      printf("# TM_QUDA: Clover field and inverse already loaded for gauge %d\n", quda_gauge_state.gauge_id);
    }
  } else {
    double atime = gettime();
    freeCloverQuda();
    reset_quda_clover_state(&quda_clover_state);
    loadCloverQuda(NULL, NULL, inv_param);
    set_quda_clover_state(&quda_clover_state, &quda_gauge_state);
    if(g_proc_id==0 && g_debug_level > 0 ) {
      printf("# TM_QUDA: Time for loadCloverQuda: %.4e\n",gettime()-atime);
    }
  }
}

void _loadGaugeQuda( const int compression ) {
  // check if the currently loaded gauge field is also the current gauge field
  // and if so, return immediately
  if( check_quda_gauge_state(&quda_gauge_state, nstore, X1, X2, X3, X0) ){
    return;
  } else {
    freeGaugeQuda();
    reset_quda_gauge_state(&quda_gauge_state);
  }

  if( inv_param.verbosity > QUDA_SILENT ){
    if(g_proc_id == 0) {
      printf("# TM_QUDA: Called _loadGaugeQuda\n");
      if( compression == 18 ){
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

#ifdef TM_USE_OMP
#pragma omp parallel
  {
#endif
  _Complex double tmpcplx;

  size_t gSize = (gauge_param.cpu_prec == QUDA_DOUBLE_PRECISION) ? sizeof(double) : sizeof(float);
  
  // now copy and reorder
#ifdef TM_USE_OMP
  #pragma omp for
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
          memcpy( &(gauge_quda[0][quda_idx]), &(g_gauge_field[tm_idx][3]), 18*gSize);
          memcpy( &(gauge_quda[1][quda_idx]), &(g_gauge_field[tm_idx][2]), 18*gSize);
          memcpy( &(gauge_quda[2][quda_idx]), &(g_gauge_field[tm_idx][1]), 18*gSize);
          memcpy( &(gauge_quda[3][quda_idx]), &(g_gauge_field[tm_idx][0]), 18*gSize);
#else
          memcpy( &(gauge_quda[0][quda_idx]), &(g_gauge_field[tm_idx][1]), 18*gSize);
          memcpy( &(gauge_quda[1][quda_idx]), &(g_gauge_field[tm_idx][2]), 18*gSize);
          memcpy( &(gauge_quda[2][quda_idx]), &(g_gauge_field[tm_idx][3]), 18*gSize);
          memcpy( &(gauge_quda[3][quda_idx]), &(g_gauge_field[tm_idx][0]), 18*gSize);
#endif
        if( compression == 18 && quda_input.fermionbc == TM_QUDA_THETABC ) {
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

  loadGaugeQuda((void*)gauge_quda, &gauge_param);
  set_quda_gauge_state(&quda_gauge_state, nstore, X1, X2, X3, X0);
}

// reorder spinor to QUDA format
void reorder_spinor_toQuda( double* sp, QudaPrecision precision, int doublet ) {
  double startTime = gettime();

  memcpy( tempSpinor, sp, (1+doublet)*VOLUME*24*sizeof(double) );

  // now copy and reorder from tempSpinor to spinor
#ifdef TM_USE_OMP
  #pragma omp parallel for
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

  double endTime = gettime();
  double diffTime = endTime - startTime;
  if(g_proc_id == 0)
    printf("# TM_QUDA: time spent in reorder_spinor_toQuda: %f secs\n", diffTime);
}

// reorder spinor from QUDA format
void reorder_spinor_fromQuda( double* sp, QudaPrecision precision, int doublet) {
  double startTime = gettime();

  memcpy( tempSpinor, sp, (1+doublet)*VOLUME*24*sizeof(double) );

  // now copy and reorder from tempSpinor to spinor
#ifdef TM_USE_OMP
  #pragma omp parallel for
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

  double endTime = gettime();
  double diffTime = endTime - startTime;
  tm_debug_printf(0,0,"# TM_QUDA: time spent in reorder_spinor_fromQuda: %f secs\n", diffTime);
}


// reorder spinor to QUDA format
void reorder_spinor_eo_toQuda(double* sp, QudaPrecision precision, int doublet, int odd) {
  double startTime = gettime();

  memcpy( tempSpinor, sp, (1+doublet)*(VOLUME/2)*24*sizeof(double) );

  // now copy and reorder from tempSpinor to spinor
#ifdef TM_USE_OMP
  #pragma omp parallel for
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
          if( oddBit == odd ){
            if( doublet ) {
              memcpy(&(sp[24*(j/2)]),        &(tempSpinor[24*tm_idx/2         ]), 24*sizeof(double));
              memcpy(&(sp[24*(j+VOLUME)/2]), &(tempSpinor[24*(tm_idx+VOLUME)/2]), 24*sizeof(double));
            }
            else {
              memcpy( &(sp[24*(j/2)]), &(tempSpinor[24*tm_idx/2]), 24*sizeof(double));
            }
          }
        }

  double endTime = gettime();
  double diffTime = endTime - startTime;
  if(g_proc_id == 0)
    printf("# TM_QUDA: time spent in reorder_spinor_toQuda: %f secs\n", diffTime);
}

// reorder spinor from QUDA format
void reorder_spinor_eo_fromQuda( double* sp, QudaPrecision precision, int doublet, int odd) {
  double startTime = gettime();

  memcpy( tempSpinor, sp, (1+doublet)*(VOLUME/2)*24*sizeof(double) );

  // now copy and reorder from tempSpinor to spinor
#ifdef TM_USE_OMP
  #pragma omp parallel for
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
          if( oddBit == odd ){
            if( doublet ) {
              memcpy(&(sp[24*tm_idx/2]),          &(tempSpinor[24*j/2         ]), 24*sizeof(double));
              memcpy(&(sp[24*(tm_idx+VOLUME)/2]), &(tempSpinor[24*(j+VOLUME)/2]), 24*sizeof(double));
            } else {
              memcpy( &(sp[24*tm_idx/2]), &(tempSpinor[24*(j/2)]), 24*sizeof(double));
            }
          }
        }

  double endTime = gettime();
  double diffTime = endTime - startTime;
  tm_debug_printf(0,0,"# TM_QUDA: time spent in reorder_spinor_fromQuda: %f secs\n", diffTime);
}


void set_boundary_conditions( CompressionType* compression ) {
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
      gauge_param.t_boundary = QUDA_PERIODIC_T;
    } else {
      gauge_param.t_boundary = QUDA_ANTI_PERIODIC_T;
    }
    link_recon = 18;
    link_recon_sloppy = 18;
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
      gauge_param.t_boundary = QUDA_ANTI_PERIODIC_T;
    } else {
      gauge_param.t_boundary = QUDA_PERIODIC_T;
    }

    link_recon = 12;
    link_recon_sloppy = *compression;

    tm_debug_printf(0, 0, 
        "\n# TM_QUDA: WARNING using %d compression with trivial (A)PBC instead "
        "of theta-BC ((t,x,y,z)*pi: (%f,%f,%f,%f))! This works fine but the residual "
        "check on the host (CPU) will fail.\n",
        *compression,X0,X1,X2,X3);
  }

  gauge_param.reconstruct = link_recon;
  gauge_param.reconstruct_sloppy = link_recon_sloppy;
  gauge_param.reconstruct_precondition = link_recon_sloppy;
}

void set_sloppy_prec( const SloppyPrecision sloppy_precision ) {
  // choose sloppy prec.
  QudaPrecision cuda_prec_sloppy;
  if( sloppy_precision==SLOPPY_DOUBLE ) {
    cuda_prec_sloppy = QUDA_DOUBLE_PRECISION;
    if(g_proc_id == 0) printf("# TM_QUDA: Using double prec. as sloppy!\n");
  }
  else if( sloppy_precision==SLOPPY_HALF ) {
    cuda_prec_sloppy = QUDA_HALF_PRECISION;
    if(g_proc_id == 0) printf("# TM_QUDA: Using half prec. as sloppy!\n");
  }
  else {
    cuda_prec_sloppy = QUDA_SINGLE_PRECISION;
    if(g_proc_id == 0) printf("# TM_QUDA: Using single prec. as sloppy!\n");
  }
  
  gauge_param.cuda_prec_sloppy = cuda_prec_sloppy;
  gauge_param.cuda_prec_refinement_sloppy = cuda_prec_sloppy;
  
  inv_param.cuda_prec_sloppy = cuda_prec_sloppy;
  inv_param.clover_cuda_prec_sloppy = cuda_prec_sloppy;
  inv_param.clover_cuda_prec_refinement_sloppy = cuda_prec_sloppy;
}



int invert_quda_direct(double * const propagator, double const * const source,
                       const int op_id) {
  
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
  set_boundary_conditions(&optr->compression_type);

  // set the sloppy precision of the mixed prec solver
  set_sloppy_prec(optr->sloppy_precision);
 
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
                            optr->maxiter);
  
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
  
  if( g_proc_id==0 && g_debug_level > 0 )
    printf("# TM_QUDA: Total time for invert_quda_direct: %.4e\n",gettime()-atotaltime); 

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
  set_boundary_conditions(&compression);
  // set the sloppy precision of the mixed prec solver
  set_sloppy_prec(sloppy_precision);
  
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
                            max_iter);

  // reorder spinor
  reorder_spinor_toQuda( (double*)spinorIn, inv_param.cpu_prec, 0 );

  // perform the inversion
  invertQuda(spinorOut, spinorIn, &inv_param);


  if( inv_param.verbosity > QUDA_SILENT )
    if(g_proc_id == 0)
      printf("# TM_QUDA: Done: %i iter / %g secs = %g Gflops\n",
             inv_param.iter, inv_param.secs, inv_param.gflops/inv_param.secs);

  // number of CG iterations
  int iteration = inv_param.iter;

  // reorder spinor
  // BaKo 20170901: not sure why the source was also re-ordered after inversion
  // we leave that commented out for now
  //reorder_spinor_fromQuda( (double*)spinorIn,  inv_param.cpu_prec, 0, NULL );
  //convert_lexic_to_eo(Even,     Odd,     solver_field[0]);
  
  reorder_spinor_fromQuda( (double*)spinorOut, inv_param.cpu_prec, 0 );
  convert_lexic_to_eo(Even_new, Odd_new, solver_field[1]);

  finalize_solver(solver_field, nr_sf);

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
  set_boundary_conditions(&compression);

  // set the sloppy precision of the mixed prec solver
  set_sloppy_prec(sloppy_precision);

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
  invertQuda(spinorOut, spinorIn, &inv_param);

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
  // IMPORTANT: use opposite TM flavor since gamma5 -> -gamma5 (until LXLYLZT prob. resolved)
  inv_param.mu = -g_mu;
  inv_param.epsilon = 0.0;

  inv_param.twist_flavor = QUDA_TWIST_SINGLET;
  inv_param.Ls = (inv_param.twist_flavor == QUDA_TWIST_NONDEG_DOUBLET ) ? 2 : 1;
  inv_param.Ls=1;
  //custom
  _setOneFlavourSolverParam(g_kappa,
                            g_c_sw,
                            g_mu,
                            0,//solver flag
                            1,//even_odd
                            1e-12,
                            1000);

  inv_param.solution_type = QUDA_MATPCDAG_MATPC_SOLUTION; 
  inv_param.solve_type = QUDA_NORMOP_PC_SOLVE;
  inv_param.matpc_type = QUDA_MATPC_ODD_ODD_ASYMMETRIC;

  void *spinorIn  = (void*)Q;
  void *spinorOut = (void*)P;

  // reorder spinor
  reorder_spinor_eo_toQuda( (double*)spinorIn, inv_param.cpu_prec, 0 ,1);

  // multiply

  
  MatQuda( spinorOut, spinorIn, &inv_param);

  // reorder spinor
  reorder_spinor_eo_fromQuda( (double*)spinorIn,  inv_param.cpu_prec, 0, 1);
  reorder_spinor_eo_fromQuda( (double*)spinorOut, inv_param.cpu_prec, 0, 1);
}


void _setOneFlavourSolverParam(const double kappa, const double c_sw, const double mu, 
                               const int solver_type, const int even_odd,
                               const double eps_sq, const int maxiter) {

  inv_param.tol = sqrt(eps_sq);
  inv_param.maxiter = maxiter;
  inv_param.Ls = 1;

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
    inv_param.compute_clover_inverse = 1;
    inv_param.compute_clover = 1;
  }
  else if( fabs(mu) > 0.0 ) {
    inv_param.twist_flavor = QUDA_TWIST_SINGLET;
    inv_param.dslash_type = QUDA_TWISTED_MASS_DSLASH;
    inv_param.matpc_type = QUDA_MATPC_EVEN_EVEN_ASYMMETRIC;
    inv_param.solution_type = QUDA_MAT_SOLUTION;
    // IMPORTANT: use opposite TM flavor since gamma5 -> -gamma5 (until LXLYLZT prob. resolved)
    inv_param.mu = -mu/2./kappa;
  }
  else if( c_sw > 0.0 ) {
    inv_param.twist_flavor = QUDA_TWIST_NO;
    inv_param.dslash_type = QUDA_CLOVER_WILSON_DSLASH;
    inv_param.matpc_type = QUDA_MATPC_EVEN_EVEN;
    inv_param.solution_type = QUDA_MAT_SOLUTION;
    inv_param.clover_order = QUDA_PACKED_CLOVER_ORDER;
    inv_param.clover_coeff = c_sw*kappa;
    inv_param.compute_clover_inverse = 1;
    inv_param.compute_clover = 1;
  }
  else {
    inv_param.twist_flavor = QUDA_TWIST_NO;
    inv_param.dslash_type = QUDA_WILSON_DSLASH;
    inv_param.matpc_type = QUDA_MATPC_EVEN_EVEN;
    inv_param.solution_type = QUDA_MAT_SOLUTION;
  }
  
  // choose solver
  if( solver_type == BICGSTAB ) {
    if(g_proc_id == 0) {printf("# TM_QUDA: Using BiCGstab!\n"); fflush(stdout);}
    inv_param.inv_type = QUDA_BICGSTAB_INVERTER;
  }
  else if ( solver_type == MG ) {
    if(g_proc_id == 0) {printf("# TM_QUDA: Using MG!\n"); fflush(stdout);}
    // coarsening does not support QUDA_MATPC_EVEN_EVEN_ASYMMETRIC
    if( inv_param.matpc_type == QUDA_MATPC_EVEN_EVEN_ASYMMETRIC ) inv_param.matpc_type = QUDA_MATPC_EVEN_EVEN;
    inv_param.inv_type = QUDA_GCR_INVERTER;
    inv_param.gcrNkrylov = quda_input.gcrNkrylov;
    inv_param.inv_type_precondition = QUDA_MG_INVERTER;
    inv_param.schwarz_type = QUDA_ADDITIVE_SCHWARZ;
    inv_param.reliable_delta = quda_input.reliable_delta;
    inv_param.precondition_cycle = 1;
    inv_param.tol_precondition = 1e-1;
    inv_param.maxiter_precondition = 1;
    // this under/overrelaxation parameter is not related to the ones
    // used in the MG 
    inv_param.omega = 1.0;
  }
  else {
    /* Here we invert the hermitean operator squared */
    inv_param.inv_type = QUDA_CG_INVERTER;
    if(g_proc_id == 0) {
      printf("# TM_QUDA: Using mixed precision CG!\n");
      fflush(stdout);
    }
  }

  // direct or norm-op. solve
  if( inv_param.inv_type == QUDA_CG_INVERTER ) {
    if( even_odd ) {
      inv_param.solve_type = QUDA_NORMERR_PC_SOLVE;
      if(g_proc_id == 0) printf("# TM_QUDA: Using EO preconditioning!\n");
    }
    else {
      inv_param.solve_type = QUDA_NORMERR_SOLVE;
      if(g_proc_id == 0) printf("# TM_QUDA: Not using EO preconditioning!\n");
    }
  }
  else {
    if( even_odd ) {
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

    // run the MG setup if required
  if( inv_param.inv_type_precondition == QUDA_MG_INVERTER ){
    // we begin by setting the inverter params for the quda_mg_param struct equal to the outer inv_param
    inv_mg_param = inv_param;
    // when the preconditioner for the outer solver has already been set below, the line just
    // above would set a preconditioner for the MG smoothers, which is not allowed
    // so we set this to NULL explicitly
    inv_mg_param.preconditioner = NULL;
    quda_mg_param.invert_param = &inv_mg_param;
    _setQudaMultigridParam(&quda_mg_param);

    if( check_quda_mg_setup_state(&quda_mg_setup_state, &quda_gauge_state, &quda_input) == TM_QUDA_MG_SETUP_RESET ){
      double atime = gettime();
      if( quda_mg_preconditioner != NULL ){
        tm_debug_printf(0,0,"# TM_QUDA: Destroying MG Preconditioner Setup\n");
        destroyMultigridQuda(quda_mg_preconditioner);
        reset_quda_mg_setup_state(&quda_mg_setup_state);
        quda_mg_preconditioner = NULL;
      }
      tm_debug_printf(0,0,"# TM_QUDA: Performing MG Preconditioner Setup\n");
      quda_mg_preconditioner = newMultigridQuda(&quda_mg_param);
      inv_param.preconditioner = quda_mg_preconditioner;
      set_quda_mg_setup_state(&quda_mg_setup_state, &quda_gauge_state);
      tm_debug_printf(0,1,"# TM_QUDA: MG Preconditioner Setup took %.3f seconds\n", gettime()-atime);
    } else if ( check_quda_mg_setup_state(&quda_mg_setup_state, &quda_gauge_state, &quda_input) == TM_QUDA_MG_SETUP_UPDATE )  {
      tm_debug_printf(0,0,"# TM_QUDA: Updating MG Preconditioner Setup for gauge %d\n", quda_gauge_state.gauge_id);
#ifdef TM_QUDA_EXPERIMENTAL
      if( quda_input.mg_eig_preserve_deflation == QUDA_BOOLEAN_YES ){
        tm_debug_printf(0,0,"# TM_QUDA: Deflation subspace for gauge %d will be re-used!\n", quda_gauge_state.gauge_id);
      }
#endif
      double atime = gettime();
      updateMultigridQuda(quda_mg_preconditioner, &quda_mg_param);
      set_quda_mg_setup_state(&quda_mg_setup_state, &quda_gauge_state);
      tm_debug_printf(0,1,"# TM_QUDA: MG Preconditioner Setup Update took %.3f seconds\n", gettime()-atime);
     } else {
      tm_debug_printf(0,0,"# TM_QUDA: Reusing MG Preconditioner Setup for gauge %d\n", quda_gauge_state.gauge_id);
    }
  }
  
  if(g_proc_id == 0 && g_debug_level > 3 && inv_param.inv_type_precondition == QUDA_MG_INVERTER){
    printf("--------------- MG InvertParam ------------------\n");
    printQudaInvertParam(quda_mg_param.invert_param);
    printf("---------------- MG MultigridParam ------------------------\n");
    printQudaMultigridParam(&quda_mg_param);
    printf("----------------------------------------\n");
  }
}



void _setQudaMultigridParam(QudaMultigridParam* mg_param) {
  QudaInvertParam *mg_inv_param = mg_param->invert_param;

  // FIXME: we also want to do MG for the ND operator, perhaps
  mg_inv_param->Ls = 1;
  mg_inv_param->sp_pad = 0;
  mg_inv_param->cl_pad = 0;

  // in the MG, the residual type should always be relative,
  // otherwisethe solver fails to converge
  // in the outer solver, we are still free to choose
  // absolute or relative
  mg_inv_param->residual_type = QUDA_L2_RELATIVE_RESIDUAL;

  mg_inv_param->preserve_source = QUDA_PRESERVE_SOURCE_NO;
  // the MG internal Gamma basis is always DEGRAND_ROSSI
  mg_inv_param->gamma_basis = QUDA_DEGRAND_ROSSI_GAMMA_BASIS;
  mg_inv_param->dirac_order = QUDA_DIRAC_ORDER;
  // just to be safe, we also set the input and output gamma basis again
  inv_param.gamma_basis = QUDA_CHIRAL_GAMMA_BASIS; // CHIRAL -> UKQCD does not seem to be supported right now...

  mg_inv_param->input_location = QUDA_CPU_FIELD_LOCATION;
  mg_inv_param->output_location = QUDA_CPU_FIELD_LOCATION;
  
  // currently, only QUDA_DIRECT_SOLVE is supported for this, thus also QUDA_MAT_SOLUTION 
  mg_inv_param->solve_type = QUDA_DIRECT_SOLVE;
  mg_inv_param->solution_type = QUDA_MAT_SOLUTION;

  mg_inv_param->dagger = QUDA_DAG_NO;

  mg_param->setup_type = QUDA_NULL_VECTOR_SETUP;
  mg_param->pre_orthonormalize = QUDA_BOOLEAN_NO;
  mg_param->post_orthonormalize = QUDA_BOOLEAN_YES;

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
    mg_param->coarse_grid_solution_type[level] = inv_param.solve_type == QUDA_DIRECT_PC_SOLVE ? QUDA_MATPC_SOLUTION : QUDA_MAT_SOLUTION;

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
    } // if(quda_input.mg_use_eig_solver[level] == QUDA_BOOLEAN_YES)
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

  mg_inv_param->verbosity = QUDA_SUMMARIZE;
  if( g_debug_level >= 3 ){
    mg_inv_param->verbosity = QUDA_VERBOSE;
  }
  mg_inv_param->verbosity_precondition = QUDA_SUMMARIZE;;
}

int invert_eo_MMd_quda(spinor * const out,
                       spinor * const in,
                       const double precision, const int max_iter,
                       const int solver_flag, const int rel_prec,
                       const int even_odd_flag, solver_params_t solver_params,
                       SloppyPrecision sloppy_precision,
                       CompressionType compression) {

  // it returns if quda is already init
  _initQuda();


  spinor ** solver_field = NULL;
  const int nr_sf = 2;
  init_solver_field(&solver_field, VOLUME/2, nr_sf);

  // no need anymore everithing is done in reorder_spinor_eo_toQuda
  //convert_eo_to_lexic(solver_field[0],  Even, Odd);
  // this is basically not necessary, but if we want to use an a nitial guess, it will be
  //  convert_eo_to_lexic(solver_field[1], Even_new, Odd_new);

  void *spinorIn  = (void*)in; // source
  void *spinorOut = (void*)out; // solution

  if ( rel_prec )
    inv_param.residual_type = QUDA_L2_RELATIVE_RESIDUAL;
  else
    inv_param.residual_type = QUDA_L2_ABSOLUTE_RESIDUAL;

  inv_param.kappa = g_kappa;

  // figure out which BC to use (theta, trivial...)
  set_boundary_conditions(&compression);
  // set the sloppy precision of the mixed prec solver
  set_sloppy_prec(sloppy_precision);
  
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
                            max_iter);
  // overwriting  inv_param set by _setOneFlavourSolverParam

  //solution_type    solve_type    Effect
  //  -------------    ----------    ------
  //  MAT              DIRECT        Solve Ax=b
  //  MATDAG_MAT       DIRECT        Solve A^dag y = b, followed by Ax=y
  //  MAT              NORMOP        Solve (A^dag A) x = (A^dag b)
  //  MATDAG_MAT       NORMOP        Solve (A^dag A) x = b
  
  //inv_param.solution_type = QUDA_MAT_SOLUTION;   
  //inv_param.solution_type =QUDA_MATDAG_MAT_SOLUTION  ;  // no
  //inv_param.solution_type = QUDA_MATPC_SOLUTION ; // no
  //inv_param.solution_type = QUDA_MATPC_DAG_SOLUTION ;//no
  inv_param.solution_type = QUDA_MATPCDAG_MATPC_SOLUTION; // # QUDA: ERROR: Source has zero norm
  //inv_param.solution_type = QUDA_MATPCDAG_MATPC_SHIFT_SOLUTION; //# QUDA: ERROR: Solution type 5 not supported 
  //inv_param.solution_type = QUDA_MATDAG_MAT_SOLUTION ;
  //inv_param.solve_type = QUDA_NORMOP_SOLVE ;
  inv_param.solve_type = QUDA_NORMOP_PC_SOLVE;
  //inv_param.solve_type = QUDA_NORMERR_PC_SOLVE ;  // QUDA: ERROR: Unpreconditioned MATDAG_MAT solution_type requires an unpreconditioned solve_type
  inv_param.matpc_type = QUDA_MATPC_ODD_ODD_ASYMMETRIC;
  //inv_param.matpc_type = QUDA_MATPC_ODD_ODD;
  //inv_param.matpc_type = QUDA_MATPC_EVEN_EVEN_ASYMMETRIC;
  //inv_param.matpc_type = QUDA_MATPC_EVEN_EVEN;

  // reorder spinor
  reorder_spinor_eo_toQuda( (double*)spinorIn, inv_param.cpu_prec, 0, 1);

  // perform the inversion
  invertQuda(spinorOut, spinorIn, &inv_param);


  if( inv_param.verbosity > QUDA_SILENT )
    if(g_proc_id == 0)
      printf("# TM_QUDA: Done: %i iter / %g secs = %g Gflops\n",
             inv_param.iter, inv_param.secs, inv_param.gflops/inv_param.secs);

  // number of CG iterations
  int iteration = inv_param.iter;

  // reorder spinor
  // BaKo 20170901: not sure why the source was also re-ordered after inversion
  // we leave that commented out for now
  //reorder_spinor_fromQuda( (double*)spinorIn,  inv_param.cpu_prec, 0, NULL );
  //convert_lexic_to_eo(Even,     Odd,     solver_field[0]);
  
  reorder_spinor_eo_fromQuda( (double*)spinorOut, inv_param.cpu_prec, 0, 1);
  //convert_lexic_to_eo(Even_new, Odd_new, solver_field[1]);

  finalize_solver(solver_field, nr_sf);

  if(iteration >= max_iter)
    return(-1);

  return(iteration);
}


int invert_eo_MMd_quda_ref(spinor * const Even_new, spinor * const Odd_new,
                   spinor * const Even, spinor * const Odd,
                   const double precision, const int max_iter,
                   const int solver_flag, const int rel_prec,
                   const int even_odd_flag, solver_params_t solver_params,
                   SloppyPrecision sloppy_precision,
                   CompressionType compression) {
// it returns if quda is already init
  _initQuda();

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
  set_boundary_conditions(&compression);
  // set the sloppy precision of the mixed prec solver
  set_sloppy_prec(sloppy_precision);
  
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
                            max_iter);
  // overriting parameters set in _setOneFlavourSolverParam
  inv_param.solution_type = QUDA_MATPCDAG_MATPC_SOLUTION; 
  inv_param.solve_type = QUDA_NORMOP_PC_SOLVE;
  inv_param.matpc_type = QUDA_MATPC_ODD_ODD_ASYMMETRIC;
  

  // reorder spinor
  reorder_spinor_toQuda( (double*)spinorIn, inv_param.cpu_prec, 0 );

  // perform the inversion
  invertQuda(spinorOut, spinorIn, &inv_param);


  if( inv_param.verbosity > QUDA_SILENT )
    if(g_proc_id == 0)
      printf("# TM_QUDA: Done: %i iter / %g secs = %g Gflops\n",
             inv_param.iter, inv_param.secs, inv_param.gflops/inv_param.secs);

  // number of CG iterations
  int iteration = inv_param.iter;

  // reorder spinor
  // BaKo 20170901: not sure why the source was also re-ordered after inversion
  // we leave that commented out for now
  //reorder_spinor_fromQuda( (double*)spinorIn,  inv_param.cpu_prec, 0, NULL );
  //convert_lexic_to_eo(Even,     Odd,     solver_field[0]);
  
  reorder_spinor_fromQuda( (double*)spinorOut, inv_param.cpu_prec, 0 );
  convert_lexic_to_eo(Even_new, Odd_new, solver_field[1]);

  finalize_solver(solver_field, nr_sf);

  if(iteration >= max_iter)
    return(-1);

  return(iteration);
}

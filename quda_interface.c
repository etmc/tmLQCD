/***********************************************************************
 *
 * Copyright (C) 2015 Mario Schroeck
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
* Author: Mario Schroeck <mario.schroeck@roma3.infn.it>
* 
* Last changes: 06/2015
*
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
*
* Notes:
*
* Minimum QUDA version is 0.7.0 (see https://github.com/lattice/quda/issues/151 
* and https://github.com/lattice/quda/issues/157).
*
*
**************************************************************************/

#include "quda_interface.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "global.h"
#include "boundary.h"
#include "linalg/convert_eo_to_lexic.h"
#include "solver/solver.h"
#include "solver/solver_field.h"
#include "gettime.h"
#include "quda.h"


// define order of the spatial indices
// default is LX-LY-LZ-T, see below def. of local lattice size, this is related to
// the gamma basis transformation from tmLQCD -> UKQCD
// for details see https://github.com/lattice/quda/issues/157
#define USE_LZ_LY_LX_T 0

// TRIVIAL_BC are trivial (anti-)periodic boundary conditions,
// i.e. 1 or -1 on last timeslice
// tmLQCD uses twisted BC, i.e. phases on all timeslices.
// if using TRIVIAL_BC: can't compare inversion result to tmLQCD
// if not using TRIVIAL_BC: BC will be applied to gauge field,
// can't use 12 parameter reconstruction
#define TRIVIAL_BC 0

#define MAX(a,b) ((a)>(b)?(a):(b))

// gauge and invert paramameter structs; init. in _initQuda()
QudaGaugeParam  gauge_param;
QudaInvertParam inv_param;

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
#ifdef MPI
  MPI_Cart_rank( g_cart_grid, n, &rank );
#endif

  return rank;
}

// variable to check if quda has been initialized
static int quda_initialized = 0;

void _initQuda() {
  if( quda_initialized )
    return;

  if( g_debug_level > 0 )
    if(g_proc_id == 0)
      printf("\n# QUDA: Detected QUDA version %d.%d.%d\n\n", QUDA_VERSION_MAJOR, QUDA_VERSION_MINOR, QUDA_VERSION_SUBMINOR);
  if( QUDA_VERSION_MAJOR == 0 && QUDA_VERSION_MINOR < 7) {
    fprintf(stderr, "Error: minimum QUDA version required is 0.7.0 (for support of chiral basis and removal of bug in mass normalization with preconditioning).\n");
    exit(-2);
  }

  gauge_param = newQudaGaugeParam();
  inv_param = newQudaInvertParam();

  // *** QUDA parameters begin here (may be modified)
  QudaPrecision cpu_prec  = QUDA_DOUBLE_PRECISION;
  QudaPrecision cuda_prec = QUDA_DOUBLE_PRECISION;
  QudaPrecision cuda_prec_sloppy = QUDA_SINGLE_PRECISION;
  QudaPrecision cuda_prec_precondition = QUDA_HALF_PRECISION;
#if TRIVIAL_BC
  gauge_param.t_boundary = ( abs(phase_0)>0.0 ? QUDA_ANTI_PERIODIC_T : QUDA_PERIODIC_T );
  QudaReconstructType link_recon = 12;
  QudaReconstructType link_recon_sloppy = 12;
  if( g_debug_level > 0 )
    if(g_proc_id == 0)
      printf("\n# QUDA: WARNING using TRIVIAL_BC! This works fine but the residual check on the host (CPU) will fail.\n");
#else
  gauge_param.t_boundary = QUDA_PERIODIC_T; // BC will be applied to gaugefield
  QudaReconstructType link_recon = 18;
  QudaReconstructType link_recon_sloppy = 18;
#endif
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
  gauge_param.reconstruct = link_recon;
  gauge_param.cuda_prec_sloppy = cuda_prec_sloppy;
  gauge_param.reconstruct_sloppy = link_recon_sloppy;
  gauge_param.cuda_prec_precondition = cuda_prec_precondition;
  gauge_param.reconstruct_precondition = link_recon_sloppy;
  gauge_param.gauge_fix = QUDA_GAUGE_FIXED_NO;

  inv_param.dagger = QUDA_DAG_NO;
  inv_param.mass_normalization = QUDA_KAPPA_NORMALIZATION;
  inv_param.solver_normalization = QUDA_DEFAULT_NORMALIZATION;

  inv_param.pipeline = 0;
  inv_param.gcrNkrylov = 10;

  // require both L2 relative and heavy quark residual to determine convergence
//  inv_param.residual_type = (QudaResidualType)(QUDA_L2_RELATIVE_RESIDUAL | QUDA_HEAVY_QUARK_RESIDUAL);
  inv_param.tol_hq = 1.0;//1e-3; // specify a tolerance for the residual for heavy quark residual
  inv_param.reliable_delta = 1e-2; // ignored by multi-shift solver

  // domain decomposition preconditioner parameters
  inv_param.inv_type_precondition = QUDA_CG_INVERTER;
  inv_param.schwarz_type = QUDA_ADDITIVE_SCHWARZ;
  inv_param.precondition_cycle = 1;
  inv_param.tol_precondition = 1e-1;
  inv_param.maxiter_precondition = 10;
  inv_param.verbosity_precondition = QUDA_SILENT;
  inv_param.cuda_prec_precondition = cuda_prec_precondition;
  inv_param.omega = 1.0;

  inv_param.cpu_prec = cpu_prec;
  inv_param.cuda_prec = cuda_prec;
  inv_param.cuda_prec_sloppy = cuda_prec_sloppy;

  inv_param.clover_cpu_prec = cpu_prec;
  inv_param.clover_cuda_prec = cuda_prec;
  inv_param.clover_cuda_prec_sloppy = cuda_prec_sloppy;
  inv_param.clover_cuda_prec_precondition = cuda_prec_precondition;

  inv_param.preserve_source = QUDA_PRESERVE_SOURCE_YES;
  inv_param.gamma_basis = QUDA_CHIRAL_GAMMA_BASIS;
  inv_param.dirac_order = QUDA_DIRAC_ORDER;

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
  else if( g_debug_level == 1 )
    inv_param.verbosity = QUDA_SUMMARIZE;
  else
    inv_param.verbosity = QUDA_VERBOSE;

  // general verbosity
  setVerbosityQuda( QUDA_SUMMARIZE, "# QUDA: ", stdout);

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
#ifdef MPI
  initQuda(-1); //sets device numbers automatically
#else
  initQuda(0);  //scalar build: use device 0
#endif
  quda_initialized = 1;
}

// finalize the QUDA library
void _endQuda() {
  if( quda_initialized ) {
    freeGaugeQuda();
    free((void*)tempSpinor);
    endQuda();
  }
}


void _loadGaugeQuda() {
  if( inv_param.verbosity > QUDA_SILENT )
    if(g_proc_id == 0)
      printf("\n# QUDA: Called _loadGaugeQuda\n");

  _Complex double tmpcplx;

  size_t gSize = (gauge_param.cpu_prec == QUDA_DOUBLE_PRECISION) ? sizeof(double) : sizeof(float);

  // now copy and reorder
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

#if !(TRIVIAL_BC)
          // apply boundary conditions
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
#endif

#endif
        }

  loadGaugeQuda((void*)gauge_quda, &gauge_param);
}


// reorder spinor to QUDA format
void reorder_spinor_toQuda( double* spinor, QudaPrecision precision, int doublet, double* spinor2 ) {
  double startTime = gettime();

  if( doublet ) {
    memcpy( tempSpinor,           spinor,  VOLUME*24*sizeof(double) );
    memcpy( tempSpinor+VOLUME*24, spinor2, VOLUME*24*sizeof(double) );
  }
  else {
    memcpy( tempSpinor, spinor, VOLUME*24*sizeof(double) );
  }

  // now copy and reorder from tempSpinor to spinor
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
            memcpy( &(spinor[24*(oddBit*VOLUME+j/2)]),          &(tempSpinor[24* tm_idx        ]), 24*sizeof(double));
            memcpy( &(spinor[24*(oddBit*VOLUME+j/2+VOLUME/2)]), &(tempSpinor[24*(tm_idx+VOLUME)]), 24*sizeof(double));
          }
          else {
            memcpy( &(spinor[24*(oddBit*VOLUME/2+j/2)]), &(tempSpinor[24*tm_idx]), 24*sizeof(double));
          }

        }

  double endTime = gettime();
  double diffTime = endTime - startTime;
  if(g_proc_id == 0)
    printf("# QUDA: time spent in reorder_spinor_toQuda: %f secs\n", diffTime);
}

// reorder spinor from QUDA format
void reorder_spinor_fromQuda( double* spinor, QudaPrecision precision, int doublet, double* spinor2 ) {
  double startTime = gettime();

  if( doublet ) {
    memcpy( tempSpinor, spinor, 2*VOLUME*24*sizeof(double) );
  }
  else {
    memcpy( tempSpinor, spinor, VOLUME*24*sizeof(double) );
  }

  // now copy and reorder from tempSpinor to spinor
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
            memcpy( &(spinor [24* tm_idx]),  &(tempSpinor[24*(oddBit*VOLUME+j/2)         ]), 24*sizeof(double));
            memcpy( &(spinor2[24*(tm_idx)]), &(tempSpinor[24*(oddBit*VOLUME+j/2+VOLUME/2)]), 24*sizeof(double));
          }
          else {
            memcpy( &(spinor[24*tm_idx]), &(tempSpinor[24*(oddBit*VOLUME/2+j/2)]), 24*sizeof(double));
          }
        }

  double endTime = gettime();
  double diffTime = endTime - startTime;
  if(g_proc_id == 0)
    printf("# QUDA: time spent in reorder_spinor_fromQuda: %f secs\n", diffTime);
}

int invert_eo_quda(spinor * const Even_new, spinor * const Odd_new,
                   spinor * const Even, spinor * const Odd,
                   const double precision, const int max_iter,
                   const int solver_flag, const int rel_prec,
                   const int even_odd_flag, solver_params_t solver_params) {
  _loadGaugeQuda();

  spinor ** solver_field = NULL;
  const int nr_sf = 2;
  init_solver_field(&solver_field, VOLUMEPLUSRAND, nr_sf);

  convert_eo_to_lexic(solver_field[0],  Even, Odd);
//  convert_eo_to_lexic(solver_field[1], Even_new, Odd_new);

  void *spinorIn  = (void*)solver_field[0]; // source
  void *spinorOut = (void*)solver_field[1]; // solution

  if ( rel_prec )
    inv_param.residual_type = QUDA_L2_RELATIVE_RESIDUAL;
  else
    inv_param.residual_type = QUDA_L2_ABSOLUTE_RESIDUAL;

  inv_param.kappa = g_kappa;

  // choose dslash type
  if( g_mu != 0.0 && g_c_sw > 0.0 ) {
    inv_param.dslash_type = QUDA_TWISTED_CLOVER_DSLASH;
    inv_param.matpc_type = QUDA_MATPC_EVEN_EVEN;
    inv_param.solution_type = QUDA_MAT_SOLUTION;
    inv_param.clover_order = QUDA_PACKED_CLOVER_ORDER;
    inv_param.mu = fabs(g_mu/2./g_kappa);
    inv_param.clover_coeff = g_c_sw*g_kappa;

  }
  else if( g_mu != 0.0 ) {
    inv_param.dslash_type = QUDA_TWISTED_MASS_DSLASH;
    inv_param.matpc_type = QUDA_MATPC_EVEN_EVEN_ASYMMETRIC;
    inv_param.solution_type = QUDA_MAT_SOLUTION;
    inv_param.mu = fabs(g_mu/2./g_kappa);
  }
  else if( g_c_sw > 0.0 ) {
    inv_param.dslash_type = QUDA_CLOVER_WILSON_DSLASH;
    inv_param.matpc_type = QUDA_MATPC_EVEN_EVEN;
    inv_param.solution_type = QUDA_MAT_SOLUTION;
    inv_param.clover_order = QUDA_PACKED_CLOVER_ORDER;
    inv_param.clover_coeff = g_c_sw*g_kappa;
  }
  else {
    inv_param.dslash_type = QUDA_WILSON_DSLASH;
    inv_param.matpc_type = QUDA_MATPC_EVEN_EVEN;
    inv_param.solution_type = QUDA_MAT_SOLUTION;
  }

  // choose solver
  if(solver_flag == BICGSTAB) {
    if(g_proc_id == 0) {printf("# QUDA: Using BiCGstab!\n"); fflush(stdout);}
    inv_param.inv_type = QUDA_BICGSTAB_INVERTER;
  }
  else {
    /* Here we invert the hermitean operator squared */
    inv_param.inv_type = QUDA_CG_INVERTER;
    if(g_proc_id == 0) {
      printf("# QUDA:  Using mixed precision CG!\n");
      printf("# mu = %f, kappa = %f\n", g_mu/2./g_kappa, g_kappa);
      fflush(stdout);
    }
  }

  // direct or norm-op. solve
  if( inv_param.inv_type == QUDA_CG_INVERTER ) {
    if( even_odd_flag ) {
      inv_param.solve_type = QUDA_NORMOP_PC_SOLVE;
      if(g_proc_id == 0) printf("# QUDA: Using preconditioning!\n");
    }
    else {
      inv_param.solve_type = QUDA_NORMOP_SOLVE;
      if(g_proc_id == 0) printf("# QUDA: Not using preconditioning!\n");
    }
  }
  else {
    if( even_odd_flag ) {
      inv_param.solve_type = QUDA_DIRECT_PC_SOLVE;
      if(g_proc_id == 0) printf("# QUDA: Using preconditioning!\n");
    }
    else {
      inv_param.solve_type = QUDA_DIRECT_SOLVE;
      if(g_proc_id == 0) printf("# QUDA: Not using preconditioning!\n");
    }
  }


  inv_param.tol = sqrt(precision)*0.25;
  inv_param.maxiter = max_iter;

  // IMPORTANT: use opposite TM flavor since gamma5 -> -gamma5 (until LXLYLZT prob. resolved)
  inv_param.twist_flavor = (g_mu < 0.0 ? QUDA_TWIST_PLUS : QUDA_TWIST_MINUS);
  inv_param.Ls = 1;

  // NULL pointers to host fields to force
  // construction instead of download of the clover field:
  if( g_c_sw > 0.0 )
    loadCloverQuda(NULL, NULL, &inv_param);

  // reorder spinor
  reorder_spinor_toQuda( (double*)spinorIn, inv_param.cpu_prec, 0, NULL );

  // perform the inversion
  invertQuda(spinorOut, spinorIn, &inv_param);

  if( inv_param.verbosity == QUDA_VERBOSE )
    if(g_proc_id == 0)
      printf("# QUDA: Device memory used:  Spinor: %f GiB,  Gauge: %f GiB, Clover: %f GiB\n",
   inv_param.spinorGiB, gauge_param.gaugeGiB, inv_param.cloverGiB);
  if( inv_param.verbosity > QUDA_SILENT )
    if(g_proc_id == 0)
      printf("# QUDA: Done: %i iter / %g secs = %g Gflops\n",
   inv_param.iter, inv_param.secs, inv_param.gflops/inv_param.secs);

  // number of CG iterations
  int iteration = inv_param.iter;

  // reorder spinor
  reorder_spinor_fromQuda( (double*)spinorIn,  inv_param.cpu_prec, 0, NULL );
  reorder_spinor_fromQuda( (double*)spinorOut, inv_param.cpu_prec, 0, NULL );
  convert_lexic_to_eo(Even,     Odd,     solver_field[0]);
  convert_lexic_to_eo(Even_new, Odd_new, solver_field[1]);

  finalize_solver(solver_field, nr_sf);
  freeGaugeQuda();

  if(iteration >= max_iter)
    return(-1);

  return(iteration);
}

int invert_doublet_eo_quda(spinor * const Even_new_s, spinor * const Odd_new_s,
                           spinor * const Even_new_c, spinor * const Odd_new_c,
                           spinor * const Even_s, spinor * const Odd_s,
                           spinor * const Even_c, spinor * const Odd_c,
                           const double precision, const int max_iter,
                           const int solver_flag, const int rel_prec, const int even_odd_flag) {
  _loadGaugeQuda();

  spinor ** solver_field = NULL;
  const int nr_sf = 4;
  init_solver_field(&solver_field, VOLUMEPLUSRAND, nr_sf);

  convert_eo_to_lexic(solver_field[0],   Even_s,  Odd_s);
  convert_eo_to_lexic(solver_field[1],   Even_c,  Odd_c);
//  convert_eo_to_lexic(g_spinor_field[DUM_DERI+1], Even_new, Odd_new);

  void *spinorIn    = (void*)solver_field[0]; // source
  void *spinorIn_c  = (void*)solver_field[1]; // charme source
  void *spinorOut   = (void*)solver_field[2]; // solution
  void *spinorOut_c = (void*)solver_field[3]; // charme solution

  if ( rel_prec )
    inv_param.residual_type = QUDA_L2_RELATIVE_RESIDUAL;
  else
    inv_param.residual_type = QUDA_L2_ABSOLUTE_RESIDUAL;

  inv_param.kappa = g_kappa;

  // IMPORTANT: use opposite TM mu-flavor since gamma5 -> -gamma5
  inv_param.mu      = -g_mubar /2./g_kappa;
  inv_param.epsilon =  g_epsbar/2./g_kappa;

  // choose dslash type
  if( g_c_sw > 0.0 ) {
    inv_param.dslash_type = QUDA_TWISTED_CLOVER_DSLASH;
    inv_param.matpc_type = QUDA_MATPC_EVEN_EVEN;
    inv_param.solution_type = QUDA_MAT_SOLUTION;
    inv_param.clover_order = QUDA_PACKED_CLOVER_ORDER;
    inv_param.clover_coeff = g_c_sw*g_kappa;
  }
  else {
    inv_param.dslash_type = QUDA_TWISTED_MASS_DSLASH;
    inv_param.matpc_type = QUDA_MATPC_EVEN_EVEN_ASYMMETRIC;
    inv_param.solution_type = QUDA_MAT_SOLUTION;
  }

  // choose solver
  if(solver_flag == BICGSTAB) {
    if(g_proc_id == 0) {printf("# QUDA: Using BiCGstab!\n"); fflush(stdout);}
    inv_param.inv_type = QUDA_BICGSTAB_INVERTER;
  }
  else {
    /* Here we invert the hermitean operator squared */
    inv_param.inv_type = QUDA_CG_INVERTER;
    if(g_proc_id == 0) {
      printf("# QUDA:  Using mixed precision CG!\n");
      printf("# mu = %f, kappa = %f\n", g_mu/2./g_kappa, g_kappa);
      fflush(stdout);
    }
  }

  if( even_odd_flag ) {
    inv_param.solve_type = QUDA_NORMOP_PC_SOLVE;
    if(g_proc_id == 0) printf("# QUDA: Using preconditioning!\n");
  }
  else {
    inv_param.solve_type = QUDA_NORMOP_SOLVE;
    if(g_proc_id == 0) printf("# QUDA: Not using preconditioning!\n");
  }

  inv_param.tol = sqrt(precision)*0.25;
  inv_param.maxiter = max_iter;

  inv_param.twist_flavor = QUDA_TWIST_NONDEG_DOUBLET;
  inv_param.Ls = 2;

  // NULL pointers to host fields to force
  // construction instead of download of the clover field:
  if( g_c_sw > 0.0 )
    loadCloverQuda(NULL, NULL, &inv_param);

  // reorder spinor
  reorder_spinor_toQuda( (double*)spinorIn,   inv_param.cpu_prec, 1, (double*)spinorIn_c );

  // perform the inversion
  invertQuda(spinorOut, spinorIn, &inv_param);

  if( inv_param.verbosity == QUDA_VERBOSE )
    if(g_proc_id == 0)
      printf("# QUDA: Device memory used:  Spinor: %f GiB,  Gauge: %f GiB, Clover: %f GiB\n",
   inv_param.spinorGiB, gauge_param.gaugeGiB, inv_param.cloverGiB);
  if( inv_param.verbosity > QUDA_SILENT )
    if(g_proc_id == 0)
      printf("# QUDA: Done: %i iter / %g secs = %g Gflops\n",
   inv_param.iter, inv_param.secs, inv_param.gflops/inv_param.secs);

  // number of CG iterations
  int iteration = inv_param.iter;

  // reorder spinor
  reorder_spinor_fromQuda( (double*)spinorIn,    inv_param.cpu_prec, 1, (double*)spinorIn_c );
  reorder_spinor_fromQuda( (double*)spinorOut,   inv_param.cpu_prec, 1, (double*)spinorOut_c );
  convert_lexic_to_eo(Even_s,     Odd_s,     solver_field[0]);
  convert_lexic_to_eo(Even_c,     Odd_c,     solver_field[1]);
  convert_lexic_to_eo(Even_new_s, Odd_new_s, solver_field[2]);
  convert_lexic_to_eo(Even_new_c, Odd_new_c, solver_field[3]);

  finalize_solver(solver_field, nr_sf);
  freeGaugeQuda();

  if(iteration >= max_iter)
    return(-1);

  return(iteration);
}

// if even_odd_flag set
void M_full_quda(spinor * const Even_new, spinor * const Odd_new,  spinor * const Even, spinor * const Odd) {
  inv_param.kappa = g_kappa;
  inv_param.mu = fabs(g_mu);
  inv_param.epsilon = 0.0;

  // IMPORTANT: use opposite TM flavor since gamma5 -> -gamma5 (until LXLYLZT prob. resolved)
  inv_param.twist_flavor = (g_mu < 0.0 ? QUDA_TWIST_PLUS : QUDA_TWIST_MINUS);
  inv_param.Ls = (inv_param.twist_flavor == QUDA_TWIST_NONDEG_DOUBLET ||
       inv_param.twist_flavor == QUDA_TWIST_DEG_DOUBLET ) ? 2 : 1;

  void *spinorIn  = (void*)g_spinor_field[DUM_DERI];   // source
  void *spinorOut = (void*)g_spinor_field[DUM_DERI+1]; // solution

  // reorder spinor
  convert_eo_to_lexic( spinorIn, Even, Odd );
  reorder_spinor_toQuda( (double*)spinorIn, inv_param.cpu_prec, 0, NULL );

  // multiply
  inv_param.solution_type = QUDA_MAT_SOLUTION;
  MatQuda( spinorOut, spinorIn, &inv_param);

  // reorder spinor
  reorder_spinor_fromQuda( (double*)spinorOut, inv_param.cpu_prec, 0, NULL );
  convert_lexic_to_eo( Even_new, Odd_new, spinorOut );
}

// no even-odd
void D_psi_quda(spinor * const P, spinor * const Q) {
  inv_param.kappa = g_kappa;
  inv_param.mu = fabs(g_mu);
  inv_param.epsilon = 0.0;

  // IMPORTANT: use opposite TM flavor since gamma5 -> -gamma5 (until LXLYLZT prob. resolved)
  inv_param.twist_flavor = (g_mu < 0.0 ? QUDA_TWIST_PLUS : QUDA_TWIST_MINUS);
  inv_param.Ls = (inv_param.twist_flavor == QUDA_TWIST_NONDEG_DOUBLET ||
       inv_param.twist_flavor == QUDA_TWIST_DEG_DOUBLET ) ? 2 : 1;

  void *spinorIn  = (void*)Q;
  void *spinorOut = (void*)P;

  // reorder spinor
  reorder_spinor_toQuda( (double*)spinorIn, inv_param.cpu_prec, 0, NULL );

  // multiply
  inv_param.solution_type = QUDA_MAT_SOLUTION;
  MatQuda( spinorOut, spinorIn, &inv_param);

  // reorder spinor
  reorder_spinor_fromQuda( (double*)spinorIn,  inv_param.cpu_prec, 0, NULL );
  reorder_spinor_fromQuda( (double*)spinorOut, inv_param.cpu_prec, 0, NULL );
}


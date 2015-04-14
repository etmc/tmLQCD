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
 ***********************************************************************/
/***********************************************************************
*
* File quda_interface.c
*
* Author: Mario Schroeck <mario.schroeck@roma3.infn.it>
* 
* Last changes: 03/2015
*
*
* Integration of the QUDA inverter for multi-GPU usage
*
* The externally accessible functions are
*
*   void _initQuda()
*     Initializes the QUDA library. Carries over the lattice size and the 
*     MPI process grid and thus must be called after initializing MPI (and 
*     after 'read_infile(argc,argv)').
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
*     those will be applied directly to the gaugefield.
*
*   double tmcgne_quda(int nmx,double res,int k,int l,int *status,int *ifail)
*     The same functionality as 'tmcgne' (see tmcg.c) but inversion is performed on 
*     the GPU using QUDA. Final residuum check is performed on the host (CPU)
*     with the function 'void tmQnohat_dble(int k,int l)' (see tmdirac.c).
*
*   void tmQnohat_quda(int k, int l)
*     The implementation of the QUDA equivalent of 'tmQnohat_dble'. 
*
*
* Notes:
*
* Minimum QUDA version is 0.7.0 (see https://github.com/lattice/quda/issues/151 
* and https://github.com/lattice/quda/issues/157).
*
* To enable compilation of the same code for QUDA usage and standard non-QUDA usage, 
* all calls of these functions should be wrapped in precompiler switches of the form
*
*   #ifdef QUDA
*     ...
*   #endif  
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

// final check of residual with DD functions on the CPU
#define FINAL_RESIDUAL_CHECK_CPU_DD 1


#define MAX(a,b) ((a)>(b)?(a):(b))

// gauge and invert paramameter structs; init. in _initQuda()
QudaGaugeParam  gauge_param;
QudaInvertParam inv_param;

// pointer to a temp. spinor, used for reordering etc.
double *tempSpinor;
// needed if even_odd_flag set
double *fullSpinor1;
double *fullSpinor2;

// function that maps coordinates in the communication grid to MPI ranks
int commsMap(const int *coords, void *fdata)
{
#if USE_LZ_LY_LX_T
  int n[4] = {coords[3], coords[2], coords[1], coords[0]};
#else
  int n[4] = {coords[3], coords[0], coords[1], coords[2]};
#endif

  int rank;

  MPI_Cart_rank( g_cart_grid, n, &rank );

  return rank;
}

void _initQuda()
{
	static initialized = 0;
	if( initialized )
		return;

  if( g_debug_level > 0 )
    printf("\nDetected QUDA version %d.%d.%d\n\n", QUDA_VERSION_MAJOR, QUDA_VERSION_MINOR, QUDA_VERSION_SUBMINOR);
//TODO
//  error_root((QUDA_VERSION_MAJOR == 0 && QUDA_VERSION_MINOR < 7),1,"_initQuda [quda_interface.c]","minimum QUDA version required is 0.7.0 (for support of chiral basis and removal of bug in mass normalization with preconditioning).");
  
    
  
  gauge_param = newQudaGaugeParam();
  inv_param = newQudaInvertParam();
  
  // *** QUDA parameters begin here (may be modified)
  QudaDslashType dslash_type = QUDA_TWISTED_MASS_DSLASH;
  QudaPrecision cpu_prec  = QUDA_DOUBLE_PRECISION;
  QudaPrecision cuda_prec = QUDA_DOUBLE_PRECISION;
  QudaPrecision cuda_prec_sloppy = QUDA_SINGLE_PRECISION;
  QudaPrecision cuda_prec_precondition = QUDA_HALF_PRECISION;
#if TRIVIAL_BC
  gauge_param.t_boundary = ( abs(phase_0)>0.0 ? QUDA_ANTI_PERIODIC_T : QUDA_PERIODIC_T );
  QudaReconstructType link_recon = 12;
  QudaReconstructType link_recon_sloppy = 12;
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

  inv_param.dslash_type = dslash_type;
  

  // offsets used only by multi-shift solver
  inv_param.num_offset = 4;
  double offset[4] = {0.01, 0.02, 0.03, 0.04};
  for (int i=0; i<inv_param.num_offset; i++) inv_param.offset[i] = offset[i];

  inv_param.matpc_type = QUDA_MATPC_EVEN_EVEN_ASYMMETRIC;
  inv_param.solution_type = QUDA_MAT_SOLUTION;

  inv_param.dagger = QUDA_DAG_NO;
  inv_param.mass_normalization = QUDA_KAPPA_NORMALIZATION;
  inv_param.solver_normalization = QUDA_DEFAULT_NORMALIZATION;
  inv_param.solve_type = QUDA_NORMOP_PC_SOLVE;
  inv_param.inv_type = QUDA_CG_INVERTER;

  inv_param.pipeline = 0;
  inv_param.gcrNkrylov = 10;

  // require both L2 relative and heavy quark residual to determine convergence
  inv_param.residual_type = (QudaResidualType)(QUDA_L2_RELATIVE_RESIDUAL | QUDA_HEAVY_QUARK_RESIDUAL);
  inv_param.tol_hq = 1.0;//1e-3; // specify a tolerance for the residual for heavy quark residual
  inv_param.reliable_delta = 1e-2; // ignored by multi-shift solver

  // domain decomposition preconditioner parameters
  inv_param.inv_type_precondition = QUDA_INVALID_INVERTER;
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

  if( g_debug_level == 0 )
    inv_param.verbosity = QUDA_SILENT;
  else if( g_debug_level == 1 )
    inv_param.verbosity = QUDA_SUMMARIZE;
  else
    inv_param.verbosity = QUDA_VERBOSE;

  
  // declare the grid mapping used for communications in a multi-GPU grid
#if USE_LZ_LY_LX_T
  int grid[4] = {NPROC3, NPROC2, NPROC1, NPROC0};
#else
  int grid[4] = {g_nproc_x, g_nproc_y, g_nproc_z, g_nproc_t};
#endif

  initCommsGridQuda(4, grid, commsMap, NULL);

  // alloc gauge_quda
  size_t gSize = (gauge_param.cpu_prec == QUDA_DOUBLE_PRECISION) ? sizeof(double) : sizeof(float);

  for (int dir = 0; dir < 4; dir++) {
    gauge_quda[dir] = (double*) malloc(VOLUME*18*gSize);
//    error_root((gauge_quda[dir] == NULL),1,"_initQuda [quda_interface.c]","malloc for gauge_quda[dir] failed");
  } 
  
  // alloc space for a temp. spinor, used throughout this module
  tempSpinor  = (double*)malloc( VOLUME*24*sizeof(double) );

  // only needed if even_odd_flag set
  fullSpinor1 = (double*)malloc( VOLUME*24*sizeof(double) );
  fullSpinor2 = (double*)malloc( VOLUME*24*sizeof(double) );

//  error_root((tempSpinor == NULL),1,"reorder_spinor_toQuda [quda_interface.c]","malloc for tempSpinor failed");

  // initialize the QUDA library
  initQuda(-1);
  initialized = 1;
}

// finalize the QUDA library
void _endQuda()
{
  freeGaugeQuda();
  free((void*)tempSpinor);
  endQuda();
}


void _loadGaugeQuda()
{
  if( inv_param.verbosity > QUDA_SILENT )
    printf("\nCalled _loadGaugeQuda\n\n");

  // update boundary if necessary  
//   if( query_flags(UDBUF_UP2DATE)!=1 )
//    copy_bnd_ud();

  _Complex double tmpcplx;

  size_t gSize = (gauge_param.cpu_prec == QUDA_DOUBLE_PRECISION) ? sizeof(double) : sizeof(float);

  // now copy and reorder 
  for( int x0=0; x0<T; x0++ )
    for( int x1=0; x1<LX; x1++ )
      for( int x2=0; x2<LY; x2++ )
        for( int x3=0; x3<LZ; x3++ )
        {
          /* ipt[x3+LZ*x2+LY*LZ*x1+LX*LY*LZ*x0] is the index of the
             point on the local lattice with cartesian coordinates
             (x0,x1,x2,x3) */
             
#if USE_LZ_LY_LX_T
          int j = x3 + LZ*x2 + LY*LZ*x1 + LX*LY*LZ*x0;
          int tm_idx   = g_ipt[x0][x1][x2][x3];
#else
          int j = x1 + LX*x2 + LY*LX*x3 + LZ*LY*LX*x0;
          int tm_idx   = x3 + LZ*x2 + LY*LZ*x1 + LX*LY*LZ*x0;//g_ipt[x0][x3][x2][x1];
#endif
          
          int oddBit = (x0+x1+x2+x3) & 1;
          int quda_idx = 18*(oddBit*VOLUME/2+j/2);
          

#if USE_LZ_LY_LX_T
            memcpy( &(gauge_quda[0][quda_idx]), pud[tm_idx][3], 18*gSize);
            memcpy( &(gauge_quda[1][quda_idx]), pud[tm_idx][2], 18*gSize);
            memcpy( &(gauge_quda[2][quda_idx]), pud[tm_idx][1], 18*gSize);
            memcpy( &(gauge_quda[3][quda_idx]), pud[tm_idx][0], 18*gSize);
#else
            memcpy( &(gauge_quda[0][quda_idx]), &(g_gauge_field[tm_idx][1]), 18*gSize);
            memcpy( &(gauge_quda[1][quda_idx]), &(g_gauge_field[tm_idx][2]), 18*gSize);
            memcpy( &(gauge_quda[2][quda_idx]), &(g_gauge_field[tm_idx][3]), 18*gSize);
            memcpy( &(gauge_quda[3][quda_idx]), &(g_gauge_field[tm_idx][0]), 18*gSize);

#if !(TRIVIAL_BC)
            // apply boundary conditions
            for( int i=0; i<9; i++ )
            {
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
void reorder_spinor_toQuda( double* spinor, QudaPrecision precision )
{
  double startTime = MPI_Wtime();
  memcpy( tempSpinor, spinor, VOLUME*24*sizeof(double) );

  // now copy and reorder from tempSpinor to spinor 
  for( int x0=0; x0<T; x0++ )
    for( int x1=0; x1<LX; x1++ )
      for( int x2=0; x2<LY; x2++ )
        for( int x3=0; x3<LZ; x3++ )
        {
#if USE_LZ_LY_LX_T
          int j = x3 + LZ*x2 + LY*LZ*x1 + LX*LY*LZ*x0;
          int tm_idx   = g_ipt[x0][x1][x2][x3];
#else
          int j = x1 + LX*x2 + LY*LX*x3 + LZ*LY*LX*x0;
          int tm_idx   = x3 + LZ*x2 + LY*LZ*x1 + LX*LY*LZ*x0;//g_ipt[x0][x3][x2][x1];
#endif
          
          int oddBit = (x0+x1+x2+x3) & 1;
//          int quda_idx = 18*(oddBit*VOLUME/2+j/2);

          memcpy( &(spinor[24*(oddBit*VOLUME/2+j/2)]), &(tempSpinor[24*tm_idx]), 24*sizeof(double));
        } 
        
  double endTime = MPI_Wtime();
  double diffTime = endTime - startTime;
  printf("time spent in reorder_spinor_toQuda: %f secs\n", diffTime);
}

// reorder spinor from QUDA format
void reorder_spinor_fromQuda( double* spinor, QudaPrecision precision )
{
  double startTime = MPI_Wtime();
  memcpy( tempSpinor, spinor, VOLUME*24*sizeof(double) );

  // now copy and reorder from tempSpinor to spinor 
  for( int x0=0; x0<T; x0++ )
    for( int x1=0; x1<LX; x1++ )
      for( int x2=0; x2<LY; x2++ )
        for( int x3=0; x3<LZ; x3++ )
        {
#if USE_LZ_LY_LX_T
          int j = x3 + LZ*x2 + LY*LZ*x1 + LX*LY*LZ*x0;
          int tm_idx   = g_ipt[x0][x1][x2][x3];
#else
          int j = x1 + LX*x2 + LY*LX*x3 + LZ*LY*LX*x0;
          int tm_idx   = x3 + LZ*x2 + LY*LZ*x1 + LX*LY*LZ*x0;//g_ipt[x0][x3][x2][x1];
#endif
          
          int oddBit = (x0+x1+x2+x3) & 1;
//          int quda_idx = 18*(oddBit*VOLUME/2+j/2);

          memcpy( &(spinor[24*tm_idx]), &(tempSpinor[24*(oddBit*VOLUME/2+j/2)]), 24*sizeof(double));
        }
        
  double endTime = MPI_Wtime();
  double diffTime = endTime - startTime;
  printf("time spent in reorder_spinor_fromQuda: %f secs\n", diffTime);
}

// if even_odd_flag set
void M_full_quda(spinor * const Even_new, spinor * const Odd_new,  spinor * const Even, spinor * const Odd)
{
  inv_param.kappa = g_kappa;
  inv_param.mu = fabs(g_mu);
  inv_param.epsilon = 0.0;

  // IMPORTANT: use opposite TM flavor since gamma5 -> -gamma5 (until LXLYLZT prob. resolved)
  inv_param.twist_flavor = (g_mu < 0.0 ? QUDA_TWIST_PLUS : QUDA_TWIST_MINUS);
  inv_param.Ls = (inv_param.twist_flavor == QUDA_TWIST_NONDEG_DOUBLET ||
       inv_param.twist_flavor == QUDA_TWIST_DEG_DOUBLET ) ? 2 : 1;

  void *spinorIn  = (void*)fullSpinor1;
  void *spinorOut = (void*)fullSpinor2;

  // reorder spinor
  convert_eo_to_lexic( spinorIn, Even, Odd );
  reorder_spinor_toQuda( (double*)spinorIn, inv_param.cpu_prec );

  // multiply
//   inv_param.solution_type = QUDA_MAT_SOLUTION;
  MatQuda( spinorOut, spinorIn, &inv_param);

  // reorder spinor
//  reorder_spinor_fromQuda( (double*)spinorIn,  inv_param.cpu_prec );
//  convert_lexic_to_eo( Even, Odd, spinorIn );

  reorder_spinor_fromQuda( (double*)spinorOut, inv_param.cpu_prec );
  convert_lexic_to_eo( Even_new, Odd_new, spinorOut );
}


// no even-odd
void D_psi_quda(spinor * const P, spinor * const Q)
{
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
  reorder_spinor_toQuda( (double*)spinorIn, inv_param.cpu_prec );

  // multiply
//   inv_param.solution_type = QUDA_MAT_SOLUTION;
  MatQuda( spinorOut, spinorIn, &inv_param);

  // reorder spinor
  reorder_spinor_fromQuda( (double*)spinorIn,  inv_param.cpu_prec );
  reorder_spinor_fromQuda( (double*)spinorOut, inv_param.cpu_prec );
}


// the following function calculates the relative residual with DD functions, no QUDA here.
double getResidualDD( int k, int l, double *nrm2 )
{
//  double startTime = MPI_Wtime();
//  double rn, locnrm2, locsrc2, globnrm2, globsrc2;
//  // copy source to temp. spinor (to recover source later)
//  memcpy( tempSpinor, psd[k][0], VOLUME*24*sizeof(double));
//
//  tmQnohat_dble(l,k); // apply op. on solution(=l) and store in source(=k): x -> A*x
//  complex_dble z;
//  z.re=-1.0;
//	z.im=0.0;
// 	mulg5_dble(VOLUME,psd[k][0]);  // A*x -> gamma_5 * A*x
//	mulc_spinor_add_dble(VOLUME,psd[k][0],(spinor_dble*)tempSpinor,z); // A*x -> A*x-b
//  locnrm2 = norm_square_dble(VOLUME,1,psd[k][0]);      // ||A*x-b||^2
//  locsrc2 = norm_square_dble(VOLUME,1,(spinor_dble*)tempSpinor);// ||b||^2
//
//  memcpy( psd[k][0], tempSpinor, VOLUME*24*sizeof(double));
//
//  MPI_Bcast(&rn,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
//  MPI_Allreduce(&locnrm2, &globnrm2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//  MPI_Allreduce(&locsrc2, &globsrc2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//
//  (*nrm2) = globnrm2;
//  rn = sqrt(globnrm2/globsrc2);
//
//  double endTime = MPI_Wtime();
//  double diffTime = endTime - startTime;
//  message("time spent in getResidualDD: %f secs\n", diffTime);
  
  return 0.0;//rn;
}

int invert_eo_quda(spinor * const Even_new, spinor * const Odd_new,
                   spinor * const Even, spinor * const Odd,
                   const double precision, const int max_iter,
                   const int solver_flag, const int rel_prec,
                   const int sub_evs_flag, const int even_odd_flag,
                   const int no_extra_masses, double * const extra_masses, solver_params_t solver_params,
                   const int id )
{
	convert_eo_to_lexic(g_spinor_field[DUM_DERI],   Even,     Odd);
//	convert_eo_to_lexic(g_spinor_field[DUM_DERI+1], Even_new, Odd_new);

    void *spinorIn  = (void*)g_spinor_field[DUM_DERI];   // source
    void *spinorOut = (void*)g_spinor_field[DUM_DERI+1]; // solution

    inv_param.maxiter = max_iter;

  //ranlxd(spinorIn,VOLUME*24); // a random source for debugging

  // get initial rel. residual (rn) and residual^2 (nrm2) from DD
  int iteration;
//  rn = getResidualDD(k,l,&nrm2);
//  double tol = res*rn;
  inv_param.tol = precision;

  // these can be set individually
  for (int i=0; i<inv_param.num_offset; i++) {
    inv_param.tol_offset[i] = inv_param.tol;
    inv_param.tol_hq_offset[i] = inv_param.tol_hq;
  }
  inv_param.kappa = g_kappa;
  inv_param.mu = fabs(g_mu/2./g_kappa);
  inv_param.epsilon = 0.0;
  // IMPORTANT: use opposite TM flavor since gamma5 -> -gamma5 (until LXLYLZT prob. resolved)
  inv_param.twist_flavor = (g_mu < 0.0 ? QUDA_TWIST_PLUS : QUDA_TWIST_MINUS);
  inv_param.Ls = (inv_param.twist_flavor == QUDA_TWIST_NONDEG_DOUBLET ||
       inv_param.twist_flavor == QUDA_TWIST_DEG_DOUBLET ) ? 2 : 1;

  // reorder spinor
  reorder_spinor_toQuda( (double*)spinorIn, inv_param.cpu_prec );

  // perform the inversion
  invertQuda(spinorOut, spinorIn, &inv_param);

  if( inv_param.verbosity == QUDA_VERBOSE )
    printf("Device memory used:\n   Spinor: %f GiB\n    Gauge: %f GiB\n",
	 inv_param.spinorGiB, gauge_param.gaugeGiB);
	if( inv_param.verbosity > QUDA_SILENT )
    printf("Done: %i iter / %g secs = %g Gflops\n",
	 inv_param.iter, inv_param.secs, inv_param.gflops/inv_param.secs);

  // number of CG iterations
  iteration = inv_param.iter;

  // reorder spinor
  reorder_spinor_fromQuda( (double*)spinorIn,  inv_param.cpu_prec );
  reorder_spinor_fromQuda( (double*)spinorOut, inv_param.cpu_prec );
  convert_lexic_to_eo(Even,     Odd,     g_spinor_field[DUM_DERI]);
  convert_lexic_to_eo(Even_new, Odd_new, g_spinor_field[DUM_DERI+1]);

//#if FINAL_RESIDUAL_CHECK_CPU_DD
//  rn = getResidualDD(k,l,&nrm2);
//  if( inv_param.verbosity > QUDA_SILENT )
//    message("CPU (DD) rel. residual is %e (tol. is %e)\n", rn, tol);
//#else
//    rn = inv_param.true_res;
//    if( inv_param.verbosity > QUDA_SILENT )
//      message("true_res is %e (tol. is %e)\n", rn, tol);
//#endif

  // negative value if the program failed
//  if (((*status)>=nmx)||(rn>tol)) //NOTE: here was AND, should be OR
//    (*status)=-1;
//  else if ((100.0*DBL_EPSILON*sqrt(nrm2))>tol)
//    (*status)=-2;

  if(iteration > max_iter) return(-1);
  	  return(iteration);
}

int invert_quda(spinor * const P, spinor * const Q, const int max_iter, double eps_sq, const int rel_prec )
{
  double startTime = MPI_Wtime();

  if( inv_param.verbosity > QUDA_SILENT )
    printf("\nCalled invert_quda\n\n");

  inv_param.maxiter = max_iter;

  void *spinorIn  = (void*)Q; // source
  void *spinorOut = (void*)P; // solution

  //ranlxd(spinorIn,VOLUME*24); // a random source for debugging

  // get initial rel. residual (rn) and residual^2 (nrm2) from DD
  int iteration;
//  rn = getResidualDD(k,l,&nrm2);
//  double tol = res*rn;
  inv_param.tol = rel_prec;

  // these can be set individually
  for (int i=0; i<inv_param.num_offset; i++) {
    inv_param.tol_offset[i] = inv_param.tol;
    inv_param.tol_hq_offset[i] = inv_param.tol_hq;
  }
  inv_param.kappa = g_kappa;
  inv_param.mu = fabs(g_mu);
  inv_param.epsilon = 0.0;
  // IMPORTANT: use opposite TM flavor since gamma5 -> -gamma5 (until LXLYLZT prob. resolved)
  inv_param.twist_flavor = (g_mu < 0.0 ? QUDA_TWIST_PLUS : QUDA_TWIST_MINUS);
  inv_param.Ls = (inv_param.twist_flavor == QUDA_TWIST_NONDEG_DOUBLET ||
       inv_param.twist_flavor == QUDA_TWIST_DEG_DOUBLET ) ? 2 : 1;

  // reorder spinor
  reorder_spinor_toQuda( (double*)spinorIn, inv_param.cpu_prec );

  // perform the inversion
  invertQuda(spinorOut, spinorIn, &inv_param);

  if( inv_param.verbosity == QUDA_VERBOSE )
    printf("Device memory used:\n   Spinor: %f GiB\n    Gauge: %f GiB\n",
	 inv_param.spinorGiB, gauge_param.gaugeGiB);
	if( inv_param.verbosity > QUDA_SILENT )
    printf("Done: %i iter / %g secs = %g Gflops\n",
	 inv_param.iter, inv_param.secs, inv_param.gflops/inv_param.secs);

  // number of CG iterations
  iteration = inv_param.iter;

  // reorder spinor
  reorder_spinor_fromQuda( (double*)spinorIn,  inv_param.cpu_prec );
  reorder_spinor_fromQuda( (double*)spinorOut, inv_param.cpu_prec );

//#if FINAL_RESIDUAL_CHECK_CPU_DD
//  rn = getResidualDD(k,l,&nrm2);
//  if( inv_param.verbosity > QUDA_SILENT )
//    message("CPU (DD) rel. residual is %e (tol. is %e)\n", rn, tol);
//#else
//    rn = inv_param.true_res;
//    if( inv_param.verbosity > QUDA_SILENT )
//      message("true_res is %e (tol. is %e)\n", rn, tol);
//#endif

  // negative value if the program failed
//  if (((*status)>=nmx)||(rn>tol)) //NOTE: here was AND, should be OR
//    (*status)=-1;
//  else if ((100.0*DBL_EPSILON*sqrt(nrm2))>tol)
//    (*status)=-2;

  double endTime = MPI_Wtime();
  double diffTime = endTime - startTime;
  printf("time spent in tmcgne_quda: %f secs\n", diffTime);
    
  if(iteration > max_iter) return(-1);
  	  return(iteration);
}


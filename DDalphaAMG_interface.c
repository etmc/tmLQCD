/***********************************************************************
 *
 * Copyright (C) 2016 Simone Bacchio, Jacob Finkenrath
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
 * Interface for DDalphaAMG
 *
 *******************************************************************************/

#include "DDalphaAMG_interface.h"

#ifndef DDalphaAMG

int mg_setup_iter;
int mg_coarse_setup_iter;
int mg_update_setup_iter;
int mg_update_gauge;
int mg_omp_num_threads;
int mg_Nvec;
int mg_lvl;
int mg_blk[4];
int mg_mixed_prec;
int mg_setup_mu_set;
int mg_no_shifts = 0;
double mg_mms_mass = 0;
double mg_setup_mu;
double mg_cmu_factor;
double mg_dtau_update;
double mg_rho_update;

void MG_init(void) {
    printf("ERROR: MG_init called but DDalphaAMG library not included.\n");
    exit(1);
}

void MG_update_gauge(double step) {
    printf("ERROR: MG_update_gauge called but DDalphaAMG library not included.\n");
    exit(1);
}

void MG_update_mu(double mu_tmLQCD, double odd_tmLQCD) {
    printf("ERROR: MG_update_mu called but DDalphaAMG library not included.\n");
    exit(1);
}

void MG_reset(void) {
    printf("ERROR: MG_reset called but DDalphaAMG library not included.\n");
    exit(1);
}

void MG_finalize(void) {
    printf("ERROR: MG_finalize called but DDalphaAMG library not included.\n");
    exit(1);
}

int MG_solver(spinor * const phi_new, spinor * const phi_old,
              const double precision, const int max_iter,const int rel_prec,
              const int N, su3 **gf, matrix_mult f) {
    printf("ERROR: MG_solver called but DDalphaAMG library not included.\n");
    exit(1);
}

int MG_solver_eo(spinor * const Even_new, spinor * const Odd_new,
                 spinor * const Even, spinor * const Odd,
                 const double precision, const int max_iter, const int rel_prec,
                 const int N, su3 **gf, matrix_mult_full f_full) {
    printf("ERROR: MG_solver_eo called but DDalphaAMG library not included.\n");
    exit(1);
}

int MG_solver_nd(spinor * const up_new, spinor * const dn_new,
		 spinor * const up_old, spinor * const dn_old,
		 const double precision, const int max_iter, const int rel_prec,
		 const int N, su3 **gf, matrix_mult_nd f) {
    printf("ERROR: MG_solver_nd called but DDalphaAMG library not included.\n");
    exit(1);
}

#else
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "boundary.h"
#include "linalg/convert_eo_to_lexic.h"
#include "solver/solver.h"
#include "solver/solver_field.h"
#include "gettime.h"
#include "read_input.h"
#include "DDalphaAMG.h"
#include "linalg_eo.h"
#include "phmc.h"
#include "operator/D_psi.h"
#include "operator/tm_operators.h"
#include "operator/tm_operators_nd.h"
#include "operator/clovertm_operators.h"
#include "operator/Hopping_Matrix.h"

//Enable variant for shifted operator in the ND sector.
//The variant is used in case of initial guess for the squared operator.
//It is faster and tests prove also to be safe (read Appendix A of arxiv:1801.##### by S.Bacchio et al.)
#define VARIANT_FOR_SHIFTED

DDalphaAMG_init mg_init;
DDalphaAMG_parameters mg_params;
DDalphaAMG_status mg_status;
int mg_do_setup=1; //if one do or redo the setup
int mg_update_gauge=1; //set to zero if gaugefield is up to date, set to one if it has to be updated
int mg_update_setup=0; //Number of additional setup iteration 
int mg_initialized=0;
int mg_setup_iter=5;
int mg_coarse_setup_iter=3;
int mg_update_setup_iter=1;
int mg_omp_num_threads=0;
int mg_Nvec=24;
int mg_lvl=3;
int mg_blk[4] = {0, 0, 0, 0};
int mg_mixed_prec=0;
int mg_setup_mu_set = 0; //flag that enable the use of mg_setup_mu in the setup phase
int mg_no_shifts = 0; // number of shifts to invert with MG. MMS-CG is used for the others at larger mass.
double mg_mms_mass = 0.1; // mass shift value for switching from MMS-CG to MG. MMS-CG is used for larger masses than the value.
double mg_setup_mu = 0.; 
double mg_cmu_factor = 1.0;
double mg_dtau_update = 0.0;
double mg_rho_update = -1.0;
double mg_tau = 0.0;
double gauge_tau = 0.0;

static int Cart_rank(MPI_Comm comm, const int coords[], int *rank) 
{
  int coords_l[4];
  
  coords_l[0]=coords[0];
  coords_l[1]=coords[3];
  coords_l[2]=coords[2];
  coords_l[3]=coords[1];
  
  return MPI_Cart_rank(comm, coords_l, rank);
}

static int Cart_coords(MPI_Comm comm, int rank, int maxdims, int coords[]) 
{
  int stat;
  
  stat=MPI_Cart_coords(comm, rank, maxdims, coords);
  int tmp=coords[1];
  coords[1]=coords[3];
  coords[3]=tmp;

  return stat;
}

static int conf_index_fct(int t, int z, int y, int x, int mu) 
{
  int id;
  
  id=(g_ipt[t][x][y][z])*72; //9*2*4
  id+=18*((mu%2==0)?mu:((mu==1)?3:1));//9*2
  
  return id;
}

static int vector_index_fct(int t, int z, int y, int x )
{
   int id;
   
   id=24*(g_ipt[t][x][y][z]);
   
   return id;
}

static inline int MG_check(spinor * const phi_new, spinor * const phi_old, const int N, const double precision, matrix_mult f) 
{
  double differ[2], residual;
  spinor ** check_vect = NULL;
  double acc_factor = 2;
  
  init_solver_field(&check_vect, VOLUMEPLUSRAND,1);
  f( check_vect[0], phi_new);
  diff( check_vect[0], check_vect[0], phi_old, N);
  differ[0] = sqrt(square_norm(check_vect[0], N, 1));
  differ[1] = sqrt(square_norm(phi_old, N, 1));
  finalize_solver(check_vect, 1);
  
  residual = differ[0]/differ[1];
  
  if( residual > precision && residual < acc_factor*precision ) {
    if(g_proc_id == 0)
      printf("WARNING: solution accepted even if the residual wasn't complitely acceptable (%e > %e) \n", residual, precision);
  } else if( residual > acc_factor*precision ) {
    if(g_proc_id == 0) {
      printf("ERROR: something bad happened... MG converged giving the wrong solution!! Trying to restart... \n");
      printf("ERROR contd: || s - f_{tmLQC} * f_{DDalphaAMG}^{-1} * s || / ||s|| = %e / %e = %e > %e \n", differ[0],differ[1],differ[0]/differ[1],precision);
    }
    return 0;
  } 

  if (g_debug_level > 0 && g_proc_id == 0)
    printf("MGTEST:  || s - f_{tmLQC} * f_{DDalphaAMG}^{-1} * s || / ||s|| = %e / %e = %e \n", differ[0],differ[1],differ[0]/differ[1]);
  
  return 1;
  
}

static inline int MG_check_nd( spinor * const up_new, spinor * const dn_new, spinor * const up_old, spinor * const dn_old,
                               const int N, const double precision, matrix_mult_nd f) 
{
  double differ[2], residual;
  spinor ** check_vect = NULL;
  double acc_factor = 4;
#ifdef VARIANT_FOR_SHIFTED
  if((  f == Qtm_pm_ndpsi_shift ||  // (Gamma5 Dh tau1)^2 - Schur complement squared with csw = 0 and shift
        f == Qsw_pm_ndpsi_shift )   // (Gamma5 Dh tau1)^2 - Schur complement squared with shift
     && g_shift != 0 )
    acc_factor = 1/sqrt(phmc_cheb_evmin/phmc_cheb_evmax + g_shift);
#endif

  init_solver_field(&check_vect, VOLUMEPLUSRAND,2);
  f( check_vect[0], check_vect[1], up_new, dn_new);
  diff( check_vect[0], check_vect[0], up_old, N);
  diff( check_vect[1], check_vect[1], dn_old, N);
  differ[0] = sqrt(square_norm(check_vect[0], N, 1)+square_norm(check_vect[1], N, 1));
  differ[1] = sqrt(square_norm(up_old, N, 1)+square_norm(dn_old, N, 1));
  finalize_solver(check_vect, 2);
  
  residual = differ[0]/differ[1];
  
  if( residual > precision && residual < acc_factor*precision ) {
    if(g_proc_id == 0)
      printf("WARNING: solution accepted even if the residual wasn't complitely acceptable (%e > %e). Max acc. factor %f.\n", residual, precision, acc_factor);
  } else if( residual > acc_factor*precision ) {
    if(g_proc_id == 0) {
      printf("ERROR: something bad happened... MG converged giving the wrong solution!! Trying to restart... \n");
      printf("ERROR contd: || s - f_{tmLQC} * f_{DDalphaAMG}^{-1} * s || / ||s|| = %e / %e = %e > %e \n", differ[0],differ[1],differ[0]/differ[1],precision);
    }
    return 0;
  } 

  if (g_debug_level > 0 && g_proc_id == 0)
    printf("MGTEST:  || s - f_{tmLQC} * f_{DDalphaAMG}^{-1} * s || / ||s|| = %e / %e = %e \n", differ[0],differ[1],differ[0]/differ[1]);
  
  return 1;
  
}

static inline int MG_mms_check_nd( spinor **const up_new, spinor **const dn_new, 
                                   spinor * const up_old, spinor * const dn_old,
                                   const double * shifts, const int no_shifts, 
                                   const int N, double * precision, matrix_mult_nd f) 
{
  double differ[2], residual;
  spinor ** check_vect = NULL;
  double acc_factor = 2;
  
  init_solver_field(&check_vect, VOLUMEPLUSRAND,2);

  for( int i = 0; i < no_shifts; i++ ) {

    g_shift = shifts[i]*shifts[i]; 

    f( check_vect[0], check_vect[1], up_new[i], dn_new[i]);
    diff( check_vect[0], check_vect[0], up_old, N);
    diff( check_vect[1], check_vect[1], dn_old, N);
    differ[0] = sqrt(square_norm(check_vect[0], N, 1)+square_norm(check_vect[1], N, 1));
    differ[1] = sqrt(square_norm(up_old, N, 1)+square_norm(dn_old, N, 1));
  
    residual = differ[0]/differ[1];
    
    if( residual > precision[i] && residual < acc_factor*precision[i] ) {
      if(g_proc_id == 0)
        printf("WARNING: solution accepted even if the residual wasn't complitely acceptable (%e > %e) \n", residual, precision[i]);
    } else if( residual > acc_factor*precision[i] ) {
      if(g_proc_id == 0) {
        printf("ERROR: something bad happened... MG converged giving the wrong solution!! Trying to restart... \n");
        printf("ERROR contd: || s - f_{tmLQC} * f_{DDalphaAMG}^{-1} * s || / ||s|| = %e / %e = %e > %e \n", differ[0],differ[1],differ[0]/differ[1],precision[i]);
      }
      finalize_solver(check_vect, 2);
      return 0;
    } 
    
    if (g_debug_level > 0 && g_proc_id == 0)
      printf("MGTEST:  || s - f_{tmLQC} * f_{DDalphaAMG}^{-1} * s || / ||s|| = %e / %e = %e \n", differ[0],differ[1],differ[0]/differ[1]);
    
  }

  finalize_solver(check_vect, 2);

  return 1;
  
}


static int MG_pre_solve( su3 **gf )
{
  
  double dtau = fabs(mg_tau-gauge_tau);
  // Checking if:
  //  mg_update_setup < mg_update_setup_iter : maybe you want to do more iteration at this run
  //  mg_dtau_update < dtau  : regular condition for update of setup
  //  mg_dtau_update < -dtau : during reversability check dtau is negative!
  //  mg_dtau_update == 0.0  : updating at every change of configuration -> valid as well if configuration changed outside the HMC
  //  mg_rho_update < 0.0    : parameter ignore
  //  mg_rho_update == rho   : updating only if this condition and the others are satisfied
  if ( mg_do_setup == 0 && mg_update_setup < mg_update_setup_iter && 
       ( (mg_dtau_update > 0.0 && dtau > 0.0 && mg_dtau_update < dtau+1e-6) ||
         (mg_dtau_update == 0.0 && mg_update_gauge == 1)         ||
         (mg_rho_update >= 0.0 && mg_rho_update == g_mu3) )) 
    mg_update_setup = mg_update_setup_iter;
  
  if(g_debug_level > 0 && g_proc_id == 0)
    printf("Tau has been increased since last MG setup update of %e\n", dtau);
  
  if (mg_initialized==0) {
    MG_init();
    mg_initialized = 1;
    if (g_proc_id == 0)
      printf("DDalphaAMG initialized\n");
    MPI_Barrier(MPI_COMM_WORLD);
  }
  
  if (mg_update_gauge==1) {
    DDalphaAMG_set_configuration( (double*) &(gf[0][0]), &mg_status );
    mg_update_gauge = 0;
    if (mg_status.success && g_proc_id == 0) 
      printf("DDalphaAMG cnfg set, plaquette %e\n", mg_status.info);
    else if ( g_proc_id == 0)
      printf("ERROR: configuration updating did not run correctly");
  }
  
  if (mg_do_setup==1) {
    if( mg_setup_mu_set ) {
      if (g_proc_id == 0)
        printf("DDalphaAMG using mu=%f during setup\n", mg_setup_mu);
      MG_update_mu(mg_setup_mu, 0); 
    } else
      MG_update_mu(g_mu, 0);
    if (g_proc_id == 0)
      printf("DDalphaAMG running setup\n");
    DDalphaAMG_setup(&mg_status);
    mg_do_setup = 0;
    mg_tau = gauge_tau;
    if (mg_status.success && g_proc_id == 0)
      printf("DDalphaAMG setup ran, time %.2f sec (%.2f %% on coarse grid)\n",
             mg_status.time, 100.*(mg_status.coarse_time/mg_status.time));
    else if ( g_proc_id == 0)
      printf("ERROR: setup procedure did not run correctly");
  }
  
  if (mg_update_setup>0) {
    if( mg_setup_mu_set ) {
      if (g_proc_id == 0)
        printf("DDalphaAMG using mu=%f during setup\n", mg_setup_mu);
      MG_update_mu(mg_setup_mu, 0); 
    } else
      MG_update_mu(g_mu, 0);
    if (g_proc_id == 0)
      printf("DDalphaAMG updating setup\n");
    DDalphaAMG_update_setup(mg_update_setup, &mg_status);
    mg_update_setup = 0;
    mg_tau = gauge_tau;
    if (mg_status.success && g_proc_id == 0)
      printf("DDalphaAMG setup ran, time %.2f sec (%.2f %% on coarse grid)\n",
             mg_status.time, 100.*(mg_status.coarse_time/mg_status.time));
    else if ( g_proc_id == 0)
      printf("ERROR: setup updating did not run correctly");
  }
  
  return mg_status.success;
}

static int MG_solve(spinor * const phi_new, spinor * const phi_old, const double precision,
                    const int N, matrix_mult f)
{
  
  // for rescaling  convention in DDalphaAMG: (4+m)*\delta_{x,y} in tmLQCD: 1*\delta_{x,y} -> rescale by 1/4+m
  double mg_scale=0.5/g_kappa;
  double *old = (double*) phi_old; 
  double *new = (double*) phi_new;
  spinor ** solver_field = NULL;
  
  if( N != VOLUME && N != VOLUME/2 ) {
    if( g_proc_id == 0 )
      printf("ERROR: N = %d in MG_solve. Expettected N == VOLUME (%d) or VOLUME/2 (%d)\n", N, VOLUME, VOLUME/2);
    return 0;
  }

  if (N==VOLUME/2) {
    init_solver_field(&solver_field, VOLUMEPLUSRAND,2);
    old = (double*) solver_field[0];
    new = (double*) solver_field[1];
    convert_odd_to_lexic( (spinor*) old, phi_old);
  }
  
  // Checking if the operator is in the list and compatible with N
  if (      f == Msw_psi ||       //          Schur complement with mu=0 on odd sites
            f == Qsw_psi ||       // Gamma5 - Schur complement with mu=0 on odd sites
            f == Mtm_plus_psi ||  //          Schur complement with plus mu 
            f == Msw_plus_psi ||  //          Schur complement with plus mu
            f == Qtm_plus_psi ||  // Gamma5 - Schur complement with plus mu 
            f == Qsw_plus_psi ||  // Gamma5 - Schur complement with plus mu
            f == Mtm_minus_psi || //          Schur complement with minus mu 
            f == Msw_minus_psi || //          Schur complement with minus mu
            f == Qtm_minus_psi || // Gamma5 - Schur complement with minus mu 
            f == Qsw_minus_psi || // Gamma5 - Schur complement with minus mu
            f == Qtm_pm_psi ||    //          Schur complement squared
            f == Qsw_pm_psi ) {   //          Schur complement squared
    if( N != VOLUME/2 && g_proc_id == 0 )
      printf("WARNING: expected N == VOLUME/2 for the required operator in MG_solve. Continuing with N == VOLUME\n");
  }
  else if ( f == D_psi ||         //          Full operator    with plus mu
            f == Q_plus_psi ||    // Gamma5 - Full operator    with plus mu 
            f == Q_minus_psi ||   // Gamma5 - Full operator    with minus mu
            f == Q_pm_psi ||      //          Full operator    squared
            f == Qsw_full_plus_psi || // Gamma5 - Full operator    with plus mu
            f == Qsw_full_minus_psi|| //Gamma5 - Full operator    with plus mu
            f == Qsw_full_pm_psi   || //          Full operator    squared
            f == Msw_full_minus_psi) {//         Full operator    with minus mu
    if( N != VOLUME && g_proc_id == 0 )
      printf("WARNING: expected N == VOLUME for the required operator in MG_solve. Continuing with N == VOLUME/2\n");
  }
  else if( g_proc_id == 0 )
    printf("WARNING: required operator unknown for MG_solve. Using standard operator: %s.\n",
           N==VOLUME?"D_psi":"Msw_plus_psi");

  // Setting mu
  if (      f == Msw_psi ||       //          Schur complement with mu=0 on odd sites
            f == Qsw_psi )        // Gamma5 - Schur complement with mu=0 on odd sites
    MG_update_mu(g_mu, -g_mu);
  else if ( f == Mtm_minus_psi || //          Schur complement with minus mu 
            f == Msw_minus_psi || //          Schur complement with minus mu
            f == Qtm_minus_psi || // Gamma5 - Schur complement with minus mu 
            f == Qsw_minus_psi || // Gamma5 - Schur complement with minus mu
            f == Qsw_full_minus_psi|| //Gamma5 - Full operator    with plus mu
            f == Msw_full_minus_psi|| //         Full operator    with minus mu
            f == Q_minus_psi )    // Gamma5 - Full operator    with minus mu
    MG_update_mu(-g_mu, -g_mu3);
  else if ( f == Mtm_plus_psi ||  //          Schur complement with plus mu 
            f == Msw_plus_psi ||  //          Schur complement with plus mu
            f == Qtm_plus_psi ||  // Gamma5 - Schur complement with plus mu 
            f == Qsw_plus_psi ||  // Gamma5 - Schur complement with plus mu
            f == D_psi ||         //          Full operator    with plus mu
            f == Q_plus_psi ||    // Gamma5 - Full operator    with plus mu 
            f == Qtm_pm_psi ||    //          Schur complement squared
            f == Qsw_pm_psi ||    //          Schur complement squared
            f == Qsw_full_plus_psi || // Gamma5 - Full operator    with plus mu
            f == Qsw_full_pm_psi   || //          Full operator    squared
            f == Q_pm_psi )       //          Full operator    squared
    MG_update_mu(g_mu, g_mu3); 
  else
    MG_update_mu(g_mu, g_mu3); 

  //Solving
  if (      f == Qtm_plus_psi ||  // Gamma5 - Schur complement with plus mu 
            f == Qsw_plus_psi ||  // Gamma5 - Schur complement with plus mu
            f == Qtm_minus_psi || // Gamma5 - Schur complement with minus mu 
            f == Qsw_minus_psi || // Gamma5 - Schur complement with minus mu 
            f == Qsw_psi ||       // Gamma5 - Schur complement with mu=0 on odd sites
            f == Q_plus_psi ||    // Gamma5 - Full operator    with plus mu 
            f == Q_minus_psi ||   // Gamma5 - Full operator    with minus mu
            f == Qsw_full_plus_psi || // Gamma5 - Full operator    with plus mu
            f == Qsw_full_minus_psi|| //Gamma5 - Full operator    with plus mu
            f == Qsw_full_pm_psi ) {  //          Full operator    squared
    mul_gamma5((spinor *const) old, VOLUME);
    DDalphaAMG_solve( new, old, precision, &mg_status );
    if( N == VOLUME ) // in case of VOLUME/2 old is a just local vector
      mul_gamma5((spinor *const) old, VOLUME);
  }
  else if ( f == Qtm_pm_psi ||    //          Schur complement squared
            f == Qsw_pm_psi ) {   //          Schur complement squared
    mg_scale *= mg_scale;
    DDalphaAMG_solve_squared_odd( new, old, precision, &mg_status );
  }
  else if ( f == Q_pm_psi ) {     //          Full operator    squared
    mg_scale *= mg_scale;
    DDalphaAMG_solve_squared( new, old, precision, &mg_status );
  }
  else if ( f == Mtm_plus_psi ||  //          Schur complement with plus mu 
            f == Msw_plus_psi ||  //          Schur complement with plus mu
            f == Mtm_minus_psi || //          Schur complement with minus mu 
            f == Msw_minus_psi || //          Schur complement with minus mu
            f == Msw_psi ||       //          Schur complement with mu=0 on odd sites
            f == D_psi ||         //          Full operator    with plus mu
            f == Msw_full_minus_psi) {//         Full operator    with minus mu
    DDalphaAMG_solve( new, old, precision, &mg_status );
  }
  else
    DDalphaAMG_solve( new, old, precision, &mg_status );
  
  if (N==VOLUME/2) {
    convert_lexic_to_odd(phi_new, (spinor*) new);
    finalize_solver(solver_field, 2);
  }
  
  mul_r(phi_new ,mg_scale, phi_new, N);

  if (g_proc_id == 0) {
    printf("Solving time %.2f sec (%.1f %% on coarse grid)\n", mg_status.time,
           100.*(mg_status.coarse_time/mg_status.time));
    printf("Total iterations on fine grid %d\n", mg_status.iter_count);
    printf("Total iterations on coarse grids %d\n", mg_status.coarse_iter_count);
    if (!mg_status.success) 
      printf("ERROR: the solver did not converge!\n");
  }
  
  return mg_status.success;
}

static int MG_solve_nd( spinor * up_new, spinor * dn_new, spinor * const up_old, spinor * const dn_old,
                        const double precision, const int N, matrix_mult_nd f)
{
  
  // for rescaling  convention in DDalphaAMG: (4+m)*\delta_{x,y} in tmLQCD: 1*\delta_{x,y} -> rescale by 1/4+m
  // moreover in the nd case, the tmLQCD is multiplied by phmc_invmaxev
  double mg_scale=0.5/g_kappa/phmc_invmaxev;
  double sqnorm;
  int init_guess = 0;
  spinor *old1 = up_old; 
  spinor *old2 = dn_old; 
  spinor *new1 = up_new, *new1tmp;
  spinor *new2 = dn_new, *new2tmp;
  spinor ** solver_field = NULL, ** oe_solver_field = NULL;
  int no_solver_field = 0;

  if( N != VOLUME && N != VOLUME/2 ) {
    if( g_proc_id == 0 )
      printf("ERROR: N = %d in MG_solve. Expettected N == VOLUME (%d) or VOLUME/2 (%d)\n", N, VOLUME, VOLUME/2);
    return 0;
  }

  if (N==VOLUME/2) no_solver_field += 4;

  // Checking if initial guess is given
  sqnorm = square_norm(up_new, N, 1);
  sqnorm += square_norm(dn_new, N, 1);
  if ( sqnorm > 1e-14 ) init_guess = 1;

  // In case of initial guess and squared operator, we do the inversion in two step and we need two more vectors
  if ( init_guess && (
            f == Qtm_pm_ndpsi ||        // (Gamma5 Dh tau1)^2 - Schur complement squared with csw = 0
            f == Qtm_pm_ndpsi_shift ||  // (Gamma5 Dh tau1)^2 - Schur complement squared with csw = 0 and shift
            f == Qsw_pm_ndpsi ||        // (Gamma5 Dh tau1)^2 - Schur complement squared
            f == Qsw_pm_ndpsi_shift ))  // (Gamma5 Dh tau1)^2 - Schur complement squared with shift
    no_solver_field += 2;

  // Allocating and assigning fields
  if(no_solver_field>0)
    init_solver_field(&solver_field, VOLUMEPLUSRAND,no_solver_field);

  int assign_solver_field = 0;
  if (N==VOLUME/2) {
    old1 = solver_field[assign_solver_field++];
    old2 = solver_field[assign_solver_field++];
    new1 = solver_field[assign_solver_field++];
    new2 = solver_field[assign_solver_field++];
    convert_odd_to_lexic(old1, up_old);
    convert_odd_to_lexic(old2, dn_old);
    set_even_to_zero(old1);
    set_even_to_zero(old2);
  }

  if ( init_guess && (
            f == Qtm_pm_ndpsi ||        // (Gamma5 Dh tau1)^2 - Schur complement squared with csw = 0
            f == Qtm_pm_ndpsi_shift ||  // (Gamma5 Dh tau1)^2 - Schur complement squared with csw = 0 and shift
            f == Qsw_pm_ndpsi ||        // (Gamma5 Dh tau1)^2 - Schur complement squared
            f == Qsw_pm_ndpsi_shift )) {// (Gamma5 Dh tau1)^2 - Schur complement squared with shift
    new1tmp = solver_field[assign_solver_field++];
    new2tmp = solver_field[assign_solver_field++];
  }

  // Reconstructing initial guess in case of oe
  if ( init_guess && N==VOLUME/2 ) {
    init_solver_field(&oe_solver_field, VOLUMEPLUSRAND, 4);
    spinor* tmp11 = oe_solver_field[0];
    spinor* tmp21 = oe_solver_field[1];
    spinor* tmp12 = oe_solver_field[2];
    spinor* tmp22 = oe_solver_field[3];

    if (g_debug_level > 2) {
      double differ[2];
      f( tmp11, tmp12, up_new, dn_new);
      diff( tmp11, tmp11, up_old, N);
      diff( tmp12, tmp12, dn_old, N);
      differ[0] = sqrt(square_norm(tmp11, N, 1)+square_norm(tmp12, N, 1));
      differ[1] = sqrt(square_norm(up_old, N, 1)+square_norm(dn_old, N, 1));
  
      if(g_proc_id == 0)
        printf("MG TEST: using initial guess. Relative residual = %e  \n", differ[0]/differ[1]);
    }

    /* Reconstruct the even sites                */
    if (    f == Qtm_pm_ndpsi       ||  // (Gamma5 Dh tau1)^2 - Schur complement squared with csw = 0
            f == Qsw_pm_ndpsi       ||  // (Gamma5 Dh tau1)^2 - Schur complement squared
            f == Qtm_pm_ndpsi_shift ||  // (Gamma5 Dh tau1)^2 - Schur complement squared with csw = 0 and shift
            f == Qsw_pm_ndpsi_shift ||  // (Gamma5 Dh tau1)^2 - Schur complement squared with shift
            f == Qtm_tau1_ndpsi_add_Ishift || // Gamma5 Dh tau1 - Schur complement with csw = 0 and plus shift
            f == Qtm_tau1_ndpsi_sub_Ishift || // Gamma5 Dh tau1 - Schur complement with csw = 0 and minus shift
            f == Qsw_tau1_ndpsi_add_Ishift || // Gamma5 Dh tau1 - Schur complement with plus shift
            f == Qsw_tau1_ndpsi_sub_Ishift ) {// Gamma5 Dh tau1 - Schur complement with minus shift
#ifdef VARIANT_FOR_SHIFTED
      if((  f == Qtm_pm_ndpsi_shift ||  // (Gamma5 Dh tau1)^2 - Schur complement squared with csw = 0 and shift
            f == Qsw_pm_ndpsi_shift )   // (Gamma5 Dh tau1)^2 - Schur complement squared with shift
         && g_shift != 0 ) {
        if( f == Qtm_pm_ndpsi_shift ) { // (Gamma5 Dh tau1)^2 - Schur complement squared with csw = 0 and shift
          Qtm_tau1_ndpsi_add_Ishift(tmp12, tmp11, up_new, dn_new); // tau1 exchange tmp11 <-> tmp12  
        } else {                        // (Gamma5 Dh tau1)^2 - Schur complement squared with shift
          Qsw_tau1_ndpsi_add_Ishift(tmp12, tmp11, up_new, dn_new); // tau1 exchange tmp11 <-> tmp12
        }
        // tau1 exchange new1tmp <-> new2tmp
        convert_odd_to_lexic( new2, tmp11);
        convert_odd_to_lexic( new1, tmp12);
        Hopping_Matrix(EO, tmp21, tmp11);
        Hopping_Matrix(EO, tmp22, tmp12);
        Msw_ee_inv_ndpsi(tmp11, tmp12, tmp21, tmp22);
        convert_even_to_lexic(new2, tmp11);
        convert_even_to_lexic(new1, tmp12);
      } else
#endif
      {
        // tau1 exchange tmp11 <-> tmp12
        Hopping_Matrix(EO, tmp12, up_new);
        Hopping_Matrix(EO, tmp11, dn_new);

        Msw_ee_inv_ndpsi(tmp21, tmp22, tmp11, tmp12);

        /* Assigning with plus sign for the even
         * since in Hopping_Matrix the minus is missing
         */
        // tau1 exchange tmp22 <-> tmp21
        convert_eo_to_lexic(new1, tmp22, up_new);
        convert_eo_to_lexic(new2, tmp21, dn_new);
      }
    } else {
      Hopping_Matrix(EO, tmp11, up_new);
      Hopping_Matrix(EO, tmp12, dn_new);

      Msw_ee_inv_ndpsi(tmp21, tmp22, tmp11, tmp12);

      /* Assigning with plus sign for the even
       * since in Hopping_Matrix the minus is missing
       */
      convert_eo_to_lexic(new1, tmp21, up_new);
      convert_eo_to_lexic(new2, tmp22, dn_new);
    }
  
    // if squared obtaining initial guess for Gamma5 Dh
    if (    f == Qtm_pm_ndpsi       ) { // (Gamma5 Dh tau1)^2 - Schur complement squared with csw = 0
      Qtm_dagger_ndpsi(tmp11, tmp12, up_new, dn_new); // tau1 Gamma5 Dh tau1
    }
    else if(f == Qsw_pm_ndpsi       ) { // (Gamma5 Dh tau1)^2 - Schur complement squared
      Qsw_dagger_ndpsi(tmp11, tmp12, up_new, dn_new); // tau1 Gamma5 Dh tau1
    }
    else if(f == Qtm_pm_ndpsi_shift ) { // (Gamma5 Dh tau1)^2 - Schur complement squared with csw = 0 and shift
      Qtm_tau1_ndpsi_sub_Ishift(tmp12, tmp11, up_new, dn_new); // tau1 exchange tmp11 <-> tmp12  
    }
    else if(f == Qsw_pm_ndpsi_shift ) { // (Gamma5 Dh tau1)^2 - Schur complement squared with shift
      Qsw_tau1_ndpsi_sub_Ishift(tmp12, tmp11, up_new, dn_new); // tau1 exchange tmp11 <-> tmp12
    }

    if (    f == Qtm_pm_ndpsi ||        // (Gamma5 Dh tau1)^2 - Schur complement squared with csw = 0
            f == Qtm_pm_ndpsi_shift ||  // (Gamma5 Dh tau1)^2 - Schur complement squared with csw = 0 and shift
            f == Qsw_pm_ndpsi ||        // (Gamma5 Dh tau1)^2 - Schur complement squared
            f == Qsw_pm_ndpsi_shift ){  // (Gamma5 Dh tau1)^2 - Schur complement squared with shift

      // tau1 exchange new1tmp <-> new2tmp
      convert_odd_to_lexic( new2tmp, tmp11);
      convert_odd_to_lexic( new1tmp, tmp12);
      Hopping_Matrix(EO, tmp21, tmp11);
      Hopping_Matrix(EO, tmp22, tmp12);
      Msw_ee_inv_ndpsi(tmp11, tmp12, tmp21, tmp22);
      convert_even_to_lexic(new2tmp, tmp11);
      convert_even_to_lexic(new1tmp, tmp12);
    } 
    finalize_solver(oe_solver_field, 4);
  } 

  // Checking if the operator is in the list and compatible with N
  if (      f == Qtm_ndpsi ||           //  Gamma5 Dh    - Schur complement with csw = 0
            f == Qsw_ndpsi ||           //  Gamma5 Dh    - Schur complement
            f == Qtm_dagger_ndpsi ||    //  Gamma5 Dh    - Schur complement with mu = -mubar and csw = 0
            f == Qsw_dagger_ndpsi ||    //  Gamma5 Dh    - Schur complement with mu = -mubar
            f == Qtm_tau1_ndpsi_add_Ishift || // Gamma5 Dh tau1 - Schur complement with csw = 0 and plus shift
            f == Qtm_tau1_ndpsi_sub_Ishift || // Gamma5 Dh tau1 - Schur complement with csw = 0 and minus shift
            f == Qsw_tau1_ndpsi_add_Ishift || // Gamma5 Dh tau1 - Schur complement with plus shift
            f == Qsw_tau1_ndpsi_sub_Ishift || // Gamma5 Dh tau1 - Schur complement with minus shift
            f == Qtm_pm_ndpsi ||        // (Gamma5 Dh tau1)^2 - Schur complement squared with csw = 0
            f == Qtm_pm_ndpsi_shift ||  // (Gamma5 Dh tau1)^2 - Schur complement squared with csw = 0 and shift
            f == Qsw_pm_ndpsi ||        // (Gamma5 Dh tau1)^2 - Schur complement squared
            f == Qsw_pm_ndpsi_shift ) { // (Gamma5 Dh tau1)^2 - Schur complement squared with shift
    if( N != VOLUME/2 && g_proc_id == 0 )
      printf("WARNING: expected N == VOLUME/2 for the required operator in MG_solve. Continuing with N == VOLUME\n");
  }
  else if ( f == D_ndpsi ) {            //  Dh
    if( N != VOLUME && g_proc_id == 0 )
      printf("WARNING: expected N == VOLUME for the required operator in MG_solve. Continuing with N == VOLUME/2\n");
  }
  else if( g_proc_id == 0 )
    printf("WARNING: required operator unknown for MG_solve. Using standard operator: %s.\n",
           N==VOLUME?"":"Qsw_ndpsi");

  // Setting mu and eps
  if (      f == Qtm_pm_ndpsi_shift ||  // (Gamma5 Dh tau1)^2 - Schur complement squared with csw = 0 and shift
            f == Qsw_pm_ndpsi_shift )   // (Gamma5 Dh tau1)^2 - Schur complement squared with shift
    MG_update_mubar_epsbar( g_mubar, g_epsbar, sqrt(g_shift) );
  else if ( f == Qtm_tau1_ndpsi_add_Ishift || // Gamma5 Dh tau1 - Schur complement with csw = 0 and plus shift
            f == Qsw_tau1_ndpsi_add_Ishift )  // Gamma5 Dh tau1 - Schur complement with plus shift
    MG_update_mubar_epsbar( g_mubar, g_epsbar, sqrt(g_shift) );
  else if ( f == Qtm_tau1_ndpsi_sub_Ishift || // Gamma5 Dh tau1 - Schur complement with csw = 0 and minus shift
            f == Qsw_tau1_ndpsi_sub_Ishift )  // Gamma5 Dh tau1 - Schur complement with minus shift
    MG_update_mubar_epsbar( g_mubar, g_epsbar, -sqrt(g_shift) );
  else if ( f == Qtm_dagger_ndpsi ||    //  Gamma5 Dh    - Schur complement with mu = -mubar csw = 0
            f == Qsw_dagger_ndpsi )     //  Gamma5 Dh    - Schur complement with mu = -mubar
    MG_update_mubar_epsbar( -g_mubar, g_epsbar, 0 );
  else if ( f == Qtm_pm_ndpsi ||        // (Gamma5 Dh tau1)^2 - Schur complement squared with csw = 0
            f == Qsw_pm_ndpsi ||        // (Gamma5 Dh tau1)^2 - Schur complement squared
            f == D_ndpsi )              //  Dh
    MG_update_mubar_epsbar( g_mubar, g_epsbar, 0 );
  else
    MG_update_mubar_epsbar( g_mubar, g_epsbar, 0 );
  
  //Solving
  if (      f == Qtm_ndpsi ||           //  Gamma5 Dh    - Schur complement with csw = 0
            f == Qsw_ndpsi ||           //  Gamma5 Dh    - Schur complement
            f == Qtm_dagger_ndpsi ||    //  Gamma5 Dh    - Schur complement with mu = -mubar csw = 0
            f == Qsw_dagger_ndpsi ) {   //  Gamma5 Dh    - Schur complement with mu = -mubar
    mul_gamma5(old1, VOLUME);
    mul_gamma5(old2, VOLUME);
    if (init_guess) {
      // Removing normalization from initial guess
      mul_r(new1, 1/mg_scale, new1, VOLUME);
      mul_r(new2, 1/mg_scale, new2, VOLUME);
      DDalphaAMG_solve_doublet_with_guess( (double*) new1, (double*) old1, (double*) new2, (double*) old2,
                                           precision, &mg_status );
    } else {
      DDalphaAMG_solve_doublet( (double*) new1, (double*) old1, (double*) new2, (double*) old2, 
                                precision, &mg_status );
    }
    if( N == VOLUME ) { // in case of VOLUME/2 old is a just local vector
      mul_gamma5(old1, VOLUME);
      mul_gamma5(old2, VOLUME);
    }
  }
  else if ( f == Qtm_tau1_ndpsi_add_Ishift || // Gamma5 Dh tau1 - Schur complement with csw = 0 and plus shift
            f == Qtm_tau1_ndpsi_sub_Ishift || // Gamma5 Dh tau1 - Schur complement with csw = 0 and minus shift
            f == Qsw_tau1_ndpsi_add_Ishift || // Gamma5 Dh tau1 - Schur complement with plus shift
            f == Qsw_tau1_ndpsi_sub_Ishift ) {// Gamma5 Dh tau1 - Schur complement with minus shift
    mul_gamma5(old1, VOLUME);
    mul_gamma5(old2, VOLUME);
    // tau1 exchange new1 <-> new2
    if (init_guess) {
      // Removing normalization from initial guess
      mul_r(new1, 1/mg_scale, new1, VOLUME);
      mul_r(new2, 1/mg_scale, new2, VOLUME);
      DDalphaAMG_solve_doublet_with_guess( (double*) new2, (double*) old1, (double*) new1, (double*) old2, 
                                           precision, &mg_status );
    } else {
      DDalphaAMG_solve_doublet( (double*) new2, (double*) old1, (double*) new1, (double*) old2, 
                                precision, &mg_status );
    }
    if( N == VOLUME ) { // in case of VOLUME/2 old is a just local vector
      mul_gamma5(old1, VOLUME);
      mul_gamma5(old2, VOLUME);
    }
  }            
  else if ( f == Qtm_pm_ndpsi ||        // (Gamma5 Dh tau1)^2 - Schur complement squared with csw = 0
            f == Qtm_pm_ndpsi_shift ||  // (Gamma5 Dh tau1)^2 - Schur complement squared with csw = 0 and shift
            f == Qsw_pm_ndpsi ||        // (Gamma5 Dh tau1)^2 - Schur complement squared
            f == Qsw_pm_ndpsi_shift ) { // (Gamma5 Dh tau1)^2 - Schur complement squared with shift
    // DDalphaAMG: tau1 gamma5 Dh tau1 gamma5 Dh
    // tmLQCD:          gamma5 Dh tau1 gamma5 Dh tau1
    if (init_guess) {
      mul_gamma5(old1, VOLUME);
      mul_gamma5(old2, VOLUME);
      // Removing normalization from initial guess
      mul_r(new1tmp, 1/mg_scale, new1tmp, VOLUME);
      mul_r(new2tmp, 1/mg_scale, new2tmp, VOLUME);
      DDalphaAMG_solve_doublet_with_guess( (double*) new2tmp, (double*) old1, (double*) new1tmp, (double*) old2,
                                           precision/2, &mg_status );
#ifdef VARIANT_FOR_SHIFTED
      if((  f == Qtm_pm_ndpsi_shift ||  // (Gamma5 Dh tau1)^2 - Schur complement squared with csw = 0 and shift
            f == Qsw_pm_ndpsi_shift )   // (Gamma5 Dh tau1)^2 - Schur complement squared with shift
         && g_shift != 0 ) {
        // Removing normalization from initial guess
        mul_r(new1, 1/mg_scale, new1, VOLUME);
        mul_r(new2, 1/mg_scale, new2, VOLUME);
        MG_update_mubar_epsbar( g_mubar, g_epsbar, -sqrt(g_shift) );
        DDalphaAMG_solve_doublet_with_guess( (double*) new2, (double*) old1, (double*) new1, (double*) old2,
                                             precision/2, &mg_status );
        assign_mul_add_mul(new1, -_Complex_I/2./sqrt(g_shift), new1tmp, _Complex_I/2./sqrt(g_shift), VOLUME);
        assign_mul_add_mul(new2, -_Complex_I/2./sqrt(g_shift), new2tmp, _Complex_I/2./sqrt(g_shift), VOLUME);
      } else 
#endif
      {
        mul_gamma5(new1tmp, VOLUME);
        mul_gamma5(new2tmp, VOLUME);
        set_even_to_zero(new1tmp);
        set_even_to_zero(new2tmp);
        // Removing normalization from initial guess
        mg_scale *= mg_scale;
        mul_r(new1, 1/mg_scale, new1, VOLUME);
        mul_r(new2, 1/mg_scale, new2, VOLUME);
        if (      f == Qtm_pm_ndpsi_shift ||  // (Gamma5 Dh tau1)^2 - Schur complement squared with csw = 0 and shift
                  f == Qsw_pm_ndpsi_shift )   // (Gamma5 Dh tau1)^2 - Schur complement squared with shift
          MG_update_mubar_epsbar( g_mubar, g_epsbar, -sqrt(g_shift) );
        DDalphaAMG_solve_doublet_with_guess( (double*) new2, (double*) new1tmp, (double*) new1, (double*) new2tmp,
                                             precision/2, &mg_status );      
      }
      if( N == VOLUME ) { // in case of VOLUME/2 old is a just local vector
        mul_gamma5(old1, VOLUME);
        mul_gamma5(old2, VOLUME);
      }
    } else {
      mg_scale *= mg_scale;
      DDalphaAMG_solve_doublet_squared_odd( (double*) new2, (double*) old2, (double*) new1, (double*) old1,
                                            precision, &mg_status );
    }
  }
  else if ( f == D_ndpsi ) {            //  Dh
    if (init_guess) {
      // Removing normalization from initial guess
      mul_r(new1, 1/mg_scale, new1, VOLUME);
      mul_r(new2, 1/mg_scale, new2, VOLUME);
      DDalphaAMG_solve_doublet_with_guess( (double*) new1, (double*) old1, (double*) new2, (double*) old2,
                                           precision, &mg_status );
    } else {
      DDalphaAMG_solve_doublet( (double*) new1, (double*) old1, (double*) new2, (double*) old2,
                                precision, &mg_status );
    }
  } else {
    if (init_guess) {
      // Removing normalization from initial guess
      mul_r(new1, 1/mg_scale, new1, VOLUME);
      mul_r(new2, 1/mg_scale, new2, VOLUME);
      DDalphaAMG_solve_doublet_with_guess( (double*) new1, (double*) old1, (double*) new2, (double*) old2,
                                           precision, &mg_status );
    } else {
      DDalphaAMG_solve_doublet( (double*) new1, (double*) old1, (double*) new2, (double*) old2,
                                precision, &mg_status );
    }
  }
  if (N==VOLUME/2) {
    convert_lexic_to_odd(up_new, new1);
    convert_lexic_to_odd(dn_new, new2);
  }
  if (no_solver_field>0)
    finalize_solver(solver_field, no_solver_field);
  mul_r(up_new ,mg_scale, up_new, N);
  mul_r(dn_new ,mg_scale, dn_new, N);
  
  if (g_proc_id == 0) {
    printf("Solving time %.2f sec (%.1f %% on coarse grid)\n", mg_status.time,
           100.*(mg_status.coarse_time/mg_status.time));
    printf("Total iterations on fine grid %d\n", mg_status.iter_count);
    printf("Total iterations on coarse grids %d\n", mg_status.coarse_iter_count);
    if (!mg_status.success) 
      printf("ERROR: the solver did not converge!\n");
  }
  
  return mg_status.success;
}

static int MG_mms_solve_nd( spinor **const up_new, spinor **const dn_new, 
                            spinor * const up_old, spinor * const dn_old,
                            const double * shifts, const int no_shifts,
                            double * precision, const int N, matrix_mult_nd f)
{
  
  // for rescaling  convention in DDalphaAMG: (4+m)*\delta_{x,y} in tmLQCD: 1*\delta_{x,y} -> rescale by 1/4+m
  // moreover in the nd case, the tmLQCD is multiplied by phmc_invmaxev
  double mg_scale=0.5/g_kappa/phmc_invmaxev;
  double *old1 = (double*) up_old; 
  double *old2 = (double*) dn_old; 
  double **new1, **new2, *mg_odd_shifts, *mg_even_shifts;
  spinor ** solver_field = NULL;

  //  if( N != VOLUME && N != VOLUME/2 ) {
  if( N != VOLUME/2 ) { // no full VOLUME functions implemented at the moment 
    if( g_proc_id == 0 )
      printf("ERROR: N = %d in MG_solve. Expettected N == VOLUME (%d) or VOLUME/2 (%d)\n", N, VOLUME, VOLUME/2);
    return 0;
  }

  new1 = (double**) malloc(no_shifts*sizeof(double*));
  new2 = (double**) malloc(no_shifts*sizeof(double*));
  mg_odd_shifts  = (double*) malloc(no_shifts*sizeof(double));
  mg_even_shifts = (double*) malloc(no_shifts*sizeof(double));

  if( N==VOLUME/2 ) {
    init_solver_field(&solver_field, VOLUMEPLUSRAND,2+2*no_shifts);
    old1 = (double*) solver_field[0];
    old2 = (double*) solver_field[1];
    convert_odd_to_lexic( (spinor*) old1, up_old);
    convert_odd_to_lexic( (spinor*) old2, dn_old);

    for( int i = 0; i < no_shifts; i++ ) {
      new1[i] = (double*) solver_field[2+2*i];
      new2[i] = (double*) solver_field[3+2*i];
    }
  } else {
    for( int i = 0; i < no_shifts; i++ ) {
      new1[i] = (double*) up_new[i];
      new2[i] = (double*) dn_new[i];
    }
  }

  // Checking if the operator is in the list and compatible with N
  if (      f == Qtm_tau1_ndpsi_add_Ishift || // Gamma5 Dh tau1 - Schur complement with csw = 0 and plus shift
            f == Qtm_tau1_ndpsi_sub_Ishift || // Gamma5 Dh tau1 - Schur complement with csw = 0 and minus shift
            f == Qsw_tau1_ndpsi_add_Ishift || // Gamma5 Dh tau1 - Schur complement with plus shift
            f == Qsw_tau1_ndpsi_sub_Ishift || // Gamma5 Dh tau1 - Schur complement with minus shift
            f == Qtm_pm_ndpsi_shift ||  // (Gamma5 Dh tau1)^2 - Schur complement squared with csw = 0 and shift
            f == Qsw_pm_ndpsi_shift ) { // (Gamma5 Dh tau1)^2 - Schur complement squared with shift
    if( N != VOLUME/2 ) {
      if( g_proc_id == 0 )
        printf("ERROR: expected N == VOLUME/2 for the required operator in MG_mms_solve_nd.\n");
      return 0;
    }
  }  else if( g_proc_id == 0 )
    printf("WARNING: required operator unknown for MG_solve. Using standard operator: %s.\n",
           N==VOLUME?"":"Qsw_pm_ndpsi_shift");

  // Setting mubar, epsbar and shifts
  if (      f == Qtm_tau1_ndpsi_add_Ishift || // Gamma5 Dh tau1 - Schur complement with csw = 0 and plus shift
            f == Qsw_tau1_ndpsi_add_Ishift || // Gamma5 Dh tau1 - Schur complement with plus shift
            f == Qtm_pm_ndpsi_shift ||  // (Gamma5 Dh tau1)^2 - Schur complement squared with csw = 0 and shift
            f == Qsw_pm_ndpsi_shift ) { // (Gamma5 Dh tau1)^2 - Schur complement squared with shift
    MG_update_mubar_epsbar( g_mubar, g_epsbar, shifts[0] );
    for( int i = 0; i < no_shifts; i++ ) {
      mg_odd_shifts[i]  = shifts[i]*mg_scale;
      mg_even_shifts[i] = 0;
    }
  }
  else if ( f == Qtm_tau1_ndpsi_sub_Ishift || // Gamma5 Dh tau1 - Schur complement with csw = 0 and minus shift
            f == Qsw_tau1_ndpsi_sub_Ishift ) {// Gamma5 Dh tau1 - Schur complement with minus shift
    MG_update_mubar_epsbar( g_mubar, g_epsbar, -shifts[0] );
    for( int i = 0; i < no_shifts; i++ ) {
      mg_odd_shifts[i]  = -shifts[i]*mg_scale;
      mg_even_shifts[i] = 0;
    }
  }

  //Solving
  if (      f == Qtm_tau1_ndpsi_add_Ishift || // Gamma5 Dh tau1 - Schur complement with csw = 0 and plus shift
            f == Qtm_tau1_ndpsi_sub_Ishift || // Gamma5 Dh tau1 - Schur complement with csw = 0 and minus shift
            f == Qsw_tau1_ndpsi_add_Ishift || // Gamma5 Dh tau1 - Schur complement with plus shift
            f == Qsw_tau1_ndpsi_sub_Ishift ) {// Gamma5 Dh tau1 - Schur complement with minus shift
    mul_gamma5((spinor *const) old1, VOLUME);
    mul_gamma5((spinor *const) old2, VOLUME);
    // tau1 exchange new1 <-> new2
    DDalphaAMG_solve_ms_doublet( new2, old1, new1, old2, mg_even_shifts, mg_odd_shifts, no_shifts, 
                                 precision, &mg_status );
    if( N == VOLUME ) { // in case of VOLUME/2 old is a just local vector
      mul_gamma5((spinor *const) old1, VOLUME);
      mul_gamma5((spinor *const) old2, VOLUME);
    }
  }            
  else if ( f == Qtm_pm_ndpsi_shift ||  // (Gamma5 Dh tau1)^2 - Schur complement squared with csw = 0 and shift
            f == Qsw_pm_ndpsi_shift ) { // (Gamma5 Dh tau1)^2 - Schur complement squared with shift
    mg_scale *= mg_scale;
    // DDalphaAMG: tau1 gamma5 Dh tau1 gamma5 Dh
    // tmLQCD:          gamma5 Dh tau1 gamma5 Dh tau1
    DDalphaAMG_solve_ms_doublet_squared_odd( new2, old2, new1, old1, mg_even_shifts, mg_odd_shifts, no_shifts,
                                             precision, &mg_status );
  }
  else
    DDalphaAMG_solve_ms_doublet( new1, old1, new2, old2, mg_even_shifts, mg_odd_shifts, no_shifts, 
                                 precision, &mg_status );

  if (N==VOLUME/2) {
    for( int i = 0; i < no_shifts; i++ ) {
      convert_lexic_to_odd(up_new[i], (spinor*) new1[i]);
      convert_lexic_to_odd(dn_new[i], (spinor*) new2[i]);
    }
    finalize_solver(solver_field, 2+2*no_shifts);
  }

  for( int i = 0; i < no_shifts; i++ ) {
    mul_r(up_new[i], mg_scale, up_new[i], N);
    mul_r(dn_new[i], mg_scale, dn_new[i], N);
  }

  if (g_proc_id == 0) {
    printf("Solving time %.2f sec (%.1f %% on coarse grid)\n", mg_status.time,
           100.*(mg_status.coarse_time/mg_status.time));
    printf("Total iterations on fine grid %d\n", mg_status.iter_count);
    printf("Total iterations on coarse grids %d\n", mg_status.coarse_iter_count);
    if (!mg_status.success) 
      printf("ERROR: the solver did not converge!\n");
  }

  free(new1);
  free(new2);
  free(mg_odd_shifts);
  free(mg_even_shifts);
  
  return mg_status.success;
}

void MG_init()
{
  mg_init.comm_cart=g_cart_grid;
  mg_init.Cart_rank=Cart_rank;
  mg_init.Cart_coords=Cart_coords;
  
  mg_init.global_lattice[0]=T*g_nproc_t;
  mg_init.global_lattice[1]=LZ*g_nproc_z;
  mg_init.global_lattice[2]=LY*g_nproc_y;
  mg_init.global_lattice[3]=LX*g_nproc_x;
  
  mg_init.procs[0]=g_nproc_t;
  mg_init.procs[1]=g_nproc_z;
  mg_init.procs[2]=g_nproc_y;
  mg_init.procs[3]=g_nproc_x;
  
  for(int i = 0; i<4; i++)
    if(mg_blk[i]==0)
      mg_blk[i]=(((L/g_nproc_x)%2==0)?(((L/g_nproc_x)%4==0)?4:2):
                 (((L/g_nproc_x)%3==0)?3:1));
  
  mg_init.block_lattice[0]=mg_blk[0];
  mg_init.block_lattice[1]=mg_blk[1];
  mg_init.block_lattice[2]=mg_blk[2];
  mg_init.block_lattice[3]=mg_blk[3];
  
  if (X0==0 && X1==0 && X2==0 && X3==0)
    mg_init.bc=0;
  else
    mg_init.bc=2;
  mg_init.theta[0] = X0;
  mg_init.theta[1] = X3;
  mg_init.theta[2] = X2;
  mg_init.theta[3] = X1;

  mg_init.number_of_levels=mg_lvl;
#ifdef TM_USE_OMP
  if(mg_omp_num_threads<=0)
      mg_init.number_openmp_threads=omp_num_threads;
  else
      mg_init.number_openmp_threads=mg_omp_num_threads;
#else
  mg_init.number_openmp_threads=1;
#endif   
  mg_init.kappa=g_kappa;
  mg_init.mu=0.5*g_mu/g_kappa;
  
  if (g_c_sw<0.00)
    mg_init.csw=0.0;
  else
    mg_init.csw=g_c_sw;

  if (reproduce_randomnumber_flag) {
    mg_init.rnd_seeds = (unsigned int *) malloc(g_nproc*sizeof(unsigned int));
    for (int i=0; i<g_nproc; i++)
      mg_init.rnd_seeds[i] = random_seed + i*1000;   
  }
  else
    mg_init.rnd_seeds = NULL;
  
  DDalphaAMG_initialize( &mg_init, &mg_params, &mg_status);

  if (reproduce_randomnumber_flag)
    free(mg_init.rnd_seeds);
  
  if (mg_status.success!=mg_lvl) {
      if (g_proc_id == 0) {
          printf("MG WARNING: %d level initialized instead of %d\n",mg_status.success,mg_lvl);
          printf("MG WARNING: parameter: mg_lvl is changed to %d\n\n",mg_status.success);
      }
      mg_lvl=mg_status.success;
  }
  
  mg_params.conf_index_fct=conf_index_fct;
  mg_params.vector_index_fct=vector_index_fct;
  
  /* in DDalphaAMG
   * Printing level:
   *  -1: silent (errors or warnings)
   *   0: minimal //default
   *   1: g_debug_level > 0
   */
  if(g_debug_level > 0) {
    mg_params.print=1;
  }
  else
    mg_params.print=0;
  
  mg_params.mu_factor[mg_lvl-1]=mg_cmu_factor; // input param mg_cmu_factor
  mg_params.mg_basis_vectors[0]=mg_Nvec;
  for (int j=1;j < mg_lvl-1; j++)
    mg_params.mg_basis_vectors[j]=fmax(28,mg_params.mg_basis_vectors[j-1]);
  
  mg_params.setup_iterations[0]=mg_setup_iter;
  mg_params.setup_iterations[1]=mg_coarse_setup_iter;
 
  // with mixed_precision = 2 the library adapt the solving precision according to the vector components
  if(mg_mixed_prec)
    mg_params.mixed_precision = 2;
  else
    mg_params.mixed_precision = 1;

  DDalphaAMG_update_parameters(&mg_params, &mg_status);
  
}

void MG_update_gauge(double step)
{
  gauge_tau += step;
  mg_update_gauge = 1;
}

void MG_update_mu(double mu_tmLQCD, double shift_tmLQCD)
{
  double mu, shift;
  mu    = 0.5 * mu_tmLQCD   /g_kappa;
  shift = 0.5 * shift_tmLQCD/g_kappa;
  
  DDalphaAMG_get_parameters(&mg_params);
  
  if (mu != mg_params.mu || shift != mg_params.mu_odd_shift || mg_params.mu_even_shift != 0.0 || mg_params.smoother_iterations != 4 ) {
    //Taking advantage of this function for updating printing in HMC
    if(g_debug_level > 0) 
      mg_params.print=1;
    else
      mg_params.print=0;

    mg_params.mu = mu;
    mg_params.mu_even_shift = 0.0;
    mg_params.mu_odd_shift = shift;
    mg_params.mu_factor[mg_lvl-1] = mg_cmu_factor;
    mg_params.epsbar = 0.0;
    mg_params.epsbar_ig5_even_shift = 0.0;
    mg_params.epsbar_ig5_odd_shift = 0.0;
    mg_params.smoother_iterations = 4;
    DDalphaAMG_update_parameters(&mg_params, &mg_status);
  }         
}

void MG_update_mubar_epsbar(double mubar_tmLQCD, double epsbar_tmLQCD, double shift_tmLQCD)
{
  double mubar, epsbar, shift;
  mubar  = 0.5 * mubar_tmLQCD /g_kappa;
  epsbar = 0.5 * epsbar_tmLQCD/g_kappa;
  shift  = 0.5 * shift_tmLQCD/g_kappa/phmc_invmaxev;
  
  DDalphaAMG_get_parameters(&mg_params);
  
  if ( mubar != mg_params.mu || mg_params.mu_odd_shift != 0.0 || mg_params.mu_even_shift != 0.0 ||
       epsbar != mg_params.epsbar || shift != mg_params.epsbar_ig5_odd_shift || mg_params.epsbar_ig5_even_shift != 0.0 || mg_params.smoother_iterations != 2 ) {
    //Taking advantage of this function for updating printing in HMC
    if(g_debug_level > 0) 
      mg_params.print=1;
    else
      mg_params.print=0;

    mg_params.mu = mubar;
    mg_params.mu_even_shift = 0.0;
    mg_params.mu_odd_shift = 0.0;
    mg_params.mu_factor[mg_lvl-1] = 1.0;
    mg_params.epsbar = epsbar;
    mg_params.epsbar_ig5_even_shift = 0.0;
    mg_params.epsbar_ig5_odd_shift = shift;
    mg_params.smoother_iterations = 2;
    DDalphaAMG_update_parameters(&mg_params, &mg_status);
  }         
}

void MG_reset() {

  if(mg_do_setup == 0)
    DDalphaAMG_free();
  
  mg_update_gauge = 1;
  mg_do_setup = 1;
  mg_update_setup = 0;
  mg_tau = 0.0;
  gauge_tau = 0.0;
}

void MG_finalize()
{
  DDalphaAMG_finalize();
}


int MG_solver(spinor * const phi_new, spinor * const phi_old,
              const double precision, const int max_iter,const int rel_prec,
              const int N, su3 **gf, matrix_mult f)
{
  
  int success=0;
  double mg_prec = rel_prec?sqrt(precision):sqrt(precision/square_norm(phi_old, N, 1));
  
  MG_pre_solve(gf);

  success = MG_solve( phi_new, phi_old, mg_prec, N, f );

  if(success && g_debug_level > 2) 
    success = MG_check( phi_new, phi_old, N, mg_prec, f );
  
  if(!success) {
    MG_reset();
    MG_pre_solve(gf);
    success = MG_solve( phi_new, phi_old, mg_prec, N, f);
    
    if(success && g_debug_level > 2) 
      success = MG_check( phi_new, phi_old, N, mg_prec, f );
  }
  
  if(!success) {
    if( g_proc_id == 0 )
      printf("ERROR: solver didn't converge after two trials!! Aborting... \n");
    //TODO: handle abort
    DDalphaAMG_finalize();
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    exit(1);
  } 
  // mg_status should have been used last time for the inversion.
  return mg_status.iter_count;
}

int MG_solver_eo(spinor * const Even_new, spinor * const Odd_new,
                 spinor * const Even, spinor * const Odd,
                 const double precision, const int max_iter, const int rel_prec,
                 const int N, su3 **gf, matrix_mult_full f_full)
{
  
  int iter_count;
  spinor ** solver_field = NULL;
  matrix_mult f;
  
  init_solver_field(&solver_field, VOLUMEPLUSRAND, 2);
  convert_eo_to_lexic(solver_field[0], Even, Odd);
  
  if (f_full == M_full)
    f=&D_psi;
  else if (f_full == Q_full)
    f=&Q_plus_psi;
  else if (f_full == Msw_full)
    f=&D_psi;
  else if (f_full == Qsw_full)
    f=&Qsw_full_plus_psi;
  else {
    f=&D_psi;
    if( g_proc_id == 0 )
      printf("WARNING: required operator unknown for MG_solver_eo. Using standard operator.\n");
  }

  iter_count = MG_solver( solver_field[1], solver_field[0], precision, max_iter, rel_prec, VOLUME, gf, f );
  
  convert_lexic_to_eo(Even_new, Odd_new, solver_field[1]);
  finalize_solver(solver_field, 2);
  
  return iter_count;
}

int MG_solver_nd(spinor * const up_new, spinor * const dn_new,
                 spinor * const up_old, spinor * const dn_old,
                 const double precision, const int max_iter, const int rel_prec,
                 const int N, su3 **gf, matrix_mult_nd f)
{
  
  int success=0;
  double mg_prec = rel_prec?sqrt(precision):sqrt(precision/(square_norm(up_old, N, 1)+square_norm(dn_old, N, 1)));
  
  MG_pre_solve(gf);

  success = MG_solve_nd( up_new, dn_new, up_old, dn_old, mg_prec, N, f );
  
  if(success && g_debug_level > 2) {
    success = MG_check_nd( up_new, dn_new, up_old, dn_old, N, mg_prec, f );

    if(!success) {
      success = MG_solve_nd( up_new, dn_new, up_old, dn_old, mg_prec, N, f);
    
      if(success) 
        success = MG_check_nd( up_new, dn_new, up_old, dn_old, N, mg_prec, f );
    }
  }
  
  if(!success) {
    MG_reset();
    MG_pre_solve(gf);
    success = MG_solve_nd( up_new, dn_new, up_old, dn_old, mg_prec, N, f);
    
    if(success && g_debug_level > 2) 
      success = MG_check_nd( up_new, dn_new, up_old, dn_old, N, mg_prec, f );
  }
  
  if(!success) {
    if( g_proc_id == 0 )
      printf("ERROR: solver didn't converge after two trials!! Aborting... \n");
    //TODO: handle abort
    DDalphaAMG_finalize();
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    exit(1);
  } 
  // mg_status should have been used last time for the inversion.
  return mg_status.iter_count;
}

int MG_solver_nd_eo(spinor * const Even_new_up, spinor * const Odd_new_up, 
                    spinor * const Even_new_dn, spinor * const Odd_new_dn,
                    spinor * const Even_up, spinor * const Odd_up,
                    spinor * const Even_dn, spinor * const Odd_dn,
                    const double precision, const int max_iter, const int rel_prec,
                    const int N, su3 **gf, matrix_mult_full_nd f_full)
{
  
  int iter_count;
  spinor ** solver_field = NULL;
  matrix_mult_nd f;
  
  init_solver_field(&solver_field, VOLUMEPLUSRAND, 4);
  convert_eo_to_lexic(solver_field[0], Even_up, Odd_up);
  convert_eo_to_lexic(solver_field[1], Even_dn, Odd_dn);
  
  if (f_full == M_full_ndpsi)
    f=&D_ndpsi;
  else if (f_full == Msw_full_ndpsi)
    f=&D_ndpsi;
  else {
    f=&D_ndpsi;
    if( g_proc_id == 0 )
      printf("WARNING: required operator unknown for MG_solver_eo. Using standard operator.\n");
  }

  iter_count = MG_solver_nd( solver_field[2], solver_field[3], solver_field[0], solver_field[1], precision, max_iter,
                             rel_prec, VOLUME, gf, f );
  
  convert_lexic_to_eo(Even_new_up, Odd_new_up, solver_field[2]);
  convert_lexic_to_eo(Even_new_dn, Odd_new_dn, solver_field[3]);
  finalize_solver(solver_field, 4);
  
  return iter_count;
}

int MG_mms_solver_nd(spinor **const up_new, spinor **const dn_new,
                     spinor * const up_old, spinor * const dn_old,
                     const double * shifts, const int no_shifts,
                     const double * precision, const int max_iter, const int rel_prec,
                     const int N, su3 **gf, matrix_mult_nd f)
{
  
  int success=0;
  double mg_prec[no_shifts];
  if(rel_prec) {
    for(int i=0; i<no_shifts; i++)
      mg_prec[i] = sqrt(precision[i]);
  } else {
    double nrhs = square_norm(up_old, N, 1)+square_norm(dn_old, N, 1);
    for(int i=0; i<no_shifts; i++)
      mg_prec[i] = sqrt(precision[i]/nrhs);
  }  

  MG_pre_solve(gf);

  success = MG_mms_solve_nd( up_new, dn_new, up_old, dn_old, shifts, no_shifts, mg_prec, N, f );
  
  if(success && g_debug_level > 2) 
    success = MG_mms_check_nd( up_new, dn_new, up_old, dn_old, shifts, no_shifts, N, mg_prec, f );
  
  if(!success) {
    MG_reset();
    MG_pre_solve(gf);
    success = MG_mms_solve_nd( up_new, dn_new, up_old, dn_old, shifts, no_shifts, mg_prec, N, f);
    
    if(success && g_debug_level > 2) 
      success = MG_mms_check_nd( up_new, dn_new, up_old, dn_old, shifts, no_shifts, N, mg_prec, f );
  }
  
  if(!success) {
    if( g_proc_id == 0 )
      printf("ERROR: solver didn't converge after two trials!! Aborting... \n");
    //TODO: handle abort
    DDalphaAMG_finalize();
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    exit(1);
  } 
  // mg_status should have been used last time for the inversion.
  return mg_status.iter_count;
}

#endif

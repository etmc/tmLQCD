/***********************************************************************
 * Copyright (C) 2012 Carsten Urbach
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
 * invert_clover_eo makes an inversion with EO preconditioned
 * clover tm Operator
 *
 * Even and Odd are the numbers of spinor_field that contain
 * the even and the odd sites of the source. The result is stored
 * int Even_new and Odd_new.
 *
 * invert_clover_eo returns the number of iterations needed or -1 if the 
 * solver did not converge.
 *
 * Author: Carsten Urbach
 *         urbach@physik.fu-berlin.de
 *
 ***********************************************************************/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include<stdlib.h>
#include"global.h"
#include"su3.h"
#include"linalg_eo.h"
#include"operator/tm_operators.h"
#include"operator/Hopping_Matrix.h"
#include"operator/clovertm_operators.h"
#include"operator/clovertm_operators_32.h"
#include"operator/D_psi.h"
#include"gamma.h"
#include"read_input.h"
#include"solver/solver.h"
#include"solver/solver_params.h"
#include"invert_clover_eo.h"
#include "solver/dirac_operator_eigenvectors.h"
#include "solver/dfl_projector.h"
#ifdef TM_USE_QUDA
#  include "quda_interface.h"
#endif
#ifdef DDalphaAMG
#  include "DDalphaAMG_interface.h"
#endif
#ifdef TM_USE_QPHIX
#  include "qphix_interface.h"
#endif

int invert_clover_eo(spinor * const Even_new, spinor * const Odd_new, 
                     spinor * const Even, spinor * const Odd,
                     const double precision, const int max_iter,
                     const int solver_flag, const int rel_prec, const int even_odd_flag,
		     solver_params_t solver_params,
                     su3 *** gf, matrix_mult Qsq, matrix_mult Qm,
                     const ExternalInverter external_inverter, const SloppyPrecision sloppy, const CompressionType compression) {
  int iter;

  if(even_odd_flag) {  
    if(g_proc_id == 0 && g_debug_level > 0) {
      printf("# Using even/odd preconditioning!\n"); fflush(stdout);
    }
    
#ifdef TM_USE_QUDA
    if( external_inverter==QUDA_INVERTER ) {
      return invert_eo_quda(Even_new, Odd_new, Even, Odd,
                            precision, max_iter,
                            solver_flag, rel_prec,
                            even_odd_flag, solver_params,
                            sloppy, compression);
    }
#endif
    
#ifdef DDalphaAMG
     if ( solver_flag == MG )
    {
      return MG_solver_eo(Even_new, Odd_new, Even, Odd, precision, max_iter,
                          rel_prec, VOLUME/2, gf[0], &Msw_full);
    }
#endif

    if(g_proc_id == 0) {
      printf("# mu = %.12f, kappa = %.12f, csw = %.12f\n", 
             g_mu/2./g_kappa, g_kappa, g_c_sw);
      fflush(stdout);
    }

    assign_mul_one_sw_pm_imu_inv(EE, Even_new, Even, +g_mu);
    
    Hopping_Matrix(OE, g_spinor_field[DUM_DERI], Even_new); 
    /* The sign is plus, since in Hopping_Matrix */
    /* the minus is missing                      */
    assign_mul_add_r(g_spinor_field[DUM_DERI], +1., Odd, VOLUME/2);
    /* Do the inversion with the preconditioned  */
    /* matrix to get the odd sites               */
    
    /* Here we invert the hermitean operator squared */
#ifdef TM_USE_QPHIX
    if( external_inverter==QPHIX_INVERTER ) {
      // QPhiX inverts M(mu)M(mu)^dag or M(mu), no gamma_5 multiplication required
      iter = invert_eo_qphix_oneflavour(Odd_new, g_spinor_field[DUM_DERI],
                                        max_iter, precision,
                                        solver_flag, rel_prec,
                                        solver_params,
                                        sloppy,
                                        compression);
      // for solver_params.solution_type == TM_SOLUTION_M (the default)
      // QPhiX applies M(mu)^dag internally for normal equation solves, no call to tmLQCD operaor required
    } else
#endif    
    if(solver_flag == CG) {
      if(g_proc_id == 0) {printf("# Using CG!\n"); fflush(stdout);}
      gamma5(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI], VOLUME/2);
      iter = cg_her(Odd_new, g_spinor_field[DUM_DERI], max_iter, 
                    precision, rel_prec, 
                    VOLUME/2, Qsq);
      Qm(Odd_new, Odd_new);
    }
    else if(solver_flag == INCREIGCG){
      gamma5(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI], VOLUME/2);
      if(g_proc_id == 0) {printf("# Using Incremental Eig-CG!\n"); fflush(stdout);}
      iter = incr_eigcg(VOLUME/2,solver_params.eigcg_nrhs, solver_params.eigcg_nrhs1, Odd_new, g_spinor_field[DUM_DERI], solver_params.eigcg_ldh, Qsq,
                        solver_params.eigcg_tolsq1, solver_params.eigcg_tolsq, solver_params.eigcg_restolsq , solver_params.eigcg_rand_guess_opt, 
                        rel_prec, max_iter, solver_params.eigcg_nev, solver_params.eigcg_vmax);
      Qm(Odd_new, Odd_new);
    }
    else if(solver_flag == MIXEDCG){
      gamma5(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI], VOLUME/2);
      iter = mixed_cg_her(Odd_new, g_spinor_field[DUM_DERI], solver_params, 
			  max_iter, precision, rel_prec, 
                          VOLUME/2, &Qsw_pm_psi, &Qsw_pm_psi_32);
      Qm(Odd_new, Odd_new);
    }
    else if(solver_flag == RGMIXEDCG){
      gamma5(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI], VOLUME/2);
      iter = rg_mixed_cg_her(Odd_new, g_spinor_field[DUM_DERI], solver_params, max_iter, precision, rel_prec,
			                     VOLUME/2, &Qsw_pm_psi, &Qsw_pm_psi_32);
      Qm(Odd_new, Odd_new);
    }
    else{
      if(g_proc_id == 0) {printf("# This solver is not available for this operator. Exiting!\n"); fflush(stdout);}
      return 0;
    }
    
    /* Reconstruct the even sites                */
    Hopping_Matrix(EO, g_spinor_field[DUM_DERI], Odd_new);
    clover_inv(g_spinor_field[DUM_DERI], +1, g_mu);
    /* The sign is plus, since in Hopping_Matrix */
    /* the minus is missing                      */
    assign_add_mul_r(Even_new, g_spinor_field[DUM_DERI], +1., VOLUME/2);
  }

  else {
    if(g_proc_id == 0) {
      printf("# Not using even/odd preconditioning!\n"); fflush(stdout);
    }
#ifdef TM_USE_QUDA
    if( external_inverter==QUDA_INVERTER ) {
      return invert_eo_quda(Even_new, Odd_new, Even, Odd,
                            precision, max_iter,
                            solver_flag, rel_prec,
                            even_odd_flag, solver_params,
                            sloppy, compression);
    }
#endif
    convert_eo_to_lexic(g_spinor_field[DUM_DERI], Even, Odd);

    if(solver_flag == DFLGCR) {
      if(g_proc_id == 0) {printf("# Using deflated FGMRES solver! m = %d\n", gmres_m_parameter); fflush(stdout);}
      iter = gcr(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI], gmres_m_parameter, 
                 max_iter/gmres_m_parameter, precision, rel_prec, VOLUME, 2, &D_psi);
    }
    else if (solver_flag == DFLFGMRES) {
      if(g_proc_id == 0) {printf("# Using deflated FGMRES solver! m = %d\n", gmres_m_parameter); fflush(stdout);}
      iter = fgmres(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI], gmres_m_parameter, 
                    max_iter/gmres_m_parameter, precision, rel_prec, VOLUME, 2, &D_psi);
    }
    else if(solver_flag == CG){
      if(g_proc_id == 0) {
	printf("# Using CG!\n"); fflush(stdout);
      }
      gamma5(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI], VOLUME);
      iter = cg_her(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+1], max_iter, precision, 
		    rel_prec, VOLUME, Qsq);
      Qm(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI]);
    }
#ifdef DDalphaAMG
    else if ( solver_flag == MG )
    {
      return MG_solver_eo(Even_new, Odd_new, Even, Odd, precision, max_iter,
                          rel_prec, VOLUME/2, gf[0], &Msw_full);
    }
#endif
    convert_lexic_to_eo(Even_new, Odd_new, g_spinor_field[DUM_DERI+1]);
  }
  return(iter);
}


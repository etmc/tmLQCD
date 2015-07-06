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
#include"operator/D_psi.h"
#include"gamma.h"
#include"solver/solver.h"
#include"invert_clover_eo.h"
#include "solver/dirac_operator_eigenvectors.h"
#ifdef QUDA
#  include "quda_interface.h"
#endif


int invert_clover_eo(spinor * const Even_new, spinor * const Odd_new, 
                     spinor * const Even, spinor * const Odd,
                     const double precision, const int max_iter,
                     const int solver_flag, const int rel_prec,solver_params_t solver_params,
                     su3 *** gf, matrix_mult Qsq, matrix_mult Qm,
                     const ExternalInverter inverter, const SloppyPrecision sloppy, const CompressionType compression) {
  int iter;

#ifdef QUDA
  if( inverter==QUDA_INVERTER ) {
    return invert_eo_quda(Even_new, Odd_new, Even, Odd,
                                  precision, max_iter,
                                  solver_flag, rel_prec,
                                  1, solver_params,
                                  sloppy, compression);
  }
#endif

  if(g_proc_id == 0 && g_debug_level > 0) {
    printf("# Using even/odd preconditioning!\n"); fflush(stdout);
  }

  assign_mul_one_sw_pm_imu_inv(EE, Even_new, Even, +g_mu);
    
  Hopping_Matrix(OE, g_spinor_field[DUM_DERI], Even_new); 
  /* The sign is plus, since in Hopping_Matrix */
  /* the minus is missing                      */
  assign_mul_add_r(g_spinor_field[DUM_DERI], +1., Odd, VOLUME/2);
  /* Do the inversion with the preconditioned  */
  /* matrix to get the odd sites               */

  /* Here we invert the hermitean operator squared */
  gamma5(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI], VOLUME/2);
  if(g_proc_id == 0) {
    //printf("# Using CG!\n"); 
    printf("# mu = %f, kappa = %f, csw = %f\n", 
	   g_mu/2./g_kappa, g_kappa, g_c_sw);
    fflush(stdout);
  }
  
  if(solver_flag == CG){
    if(g_proc_id == 0) {printf("# Using CG!\n"); fflush(stdout);}
    iter = cg_her(Odd_new, g_spinor_field[DUM_DERI], max_iter, 
		 precision, rel_prec, 
		 VOLUME/2, Qsq);
    Qm(Odd_new, Odd_new);
    }else if(solver_flag == INCREIGCG){

       if(g_proc_id == 0) {printf("# Using Incremental Eig-CG!\n"); fflush(stdout);}
       iter = incr_eigcg(VOLUME/2,solver_params.eigcg_nrhs, solver_params.eigcg_nrhs1, Odd_new, g_spinor_field[DUM_DERI], solver_params.eigcg_ldh, Qsq,
 		    	            solver_params.eigcg_tolsq1, solver_params.eigcg_tolsq, solver_params.eigcg_restolsq , solver_params.eigcg_rand_guess_opt, 
                                    rel_prec, max_iter, solver_params.eigcg_nev, solver_params.eigcg_vmax);
       Qm(Odd_new, Odd_new);

   }else{
    if(g_proc_id == 0) {printf("# This solver is not available for this operator. Exisiting!\n"); fflush(stdout);}
    return 0;
  }


  /* Reconstruct the even sites                */
  Hopping_Matrix(EO, g_spinor_field[DUM_DERI], Odd_new);
  clover_inv(g_spinor_field[DUM_DERI], +1, g_mu);
  /* The sign is plus, since in Hopping_Matrix */
  /* the minus is missing                      */
  assign_add_mul_r(Even_new, g_spinor_field[DUM_DERI], +1., VOLUME/2);

  return(iter);
}

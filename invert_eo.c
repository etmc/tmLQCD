/***********************************************************************
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008 Carsten Urbach
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
 * invert_eo makes an inversion with EO preconditioned
 * tm Operator
 *
 * Even and Odd are the numbers of spinor_field that contain
 * the even and the odd sites of the source. The result is stored
 * int Even_new and Odd_new.
 *
 * invert_eo returns the number of iterations needed or -1 if the 
 * solver did not converge.
 *
 * Author: Carsten Urbach
 *         urbach@physik.fu-berlin.de
 *
 ****************************************************************/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include<stdlib.h>
#include"global.h"
#include"linalg_eo.h"
#include"operator/tm_operators.h"
#include"operator/Hopping_Matrix.h"
#include"operator/D_psi.h"
#include"operator/tm_operators_32.h"
#include"gamma.h"
#include"solver/solver.h"
#include"read_input.h"
#include"xchange/xchange.h"
#include"solver/poly_precon.h"
#include"solver/dfl_projector.h"
#include"invert_eo.h"
#include "solver/dirac_operator_eigenvectors.h"
/* FIXME temporary includes and declarations until IO and interface for invert and CGMMS are generelized */
#include "init/init_spinor_field.h"
#include <io/params.h>
#include <io/spinor.h>
#ifdef TM_USE_QUDA
#  include "quda_interface.h"
#endif
#ifdef TM_USE_QPHIX
#  include "qphix_interface.h"
#endif
#ifdef DDalphaAMG
#  include "DDalphaAMG_interface.h"
#endif

static double cgmms_reached_prec = 0.0; 
static void cgmms_write_props(spinor ** const P, double const * const extra_masses, const int no_extra_masses, const int id, const int iteration);

#ifdef HAVE_GPU
#include"GPU/cudadefs.h"
#include"temporalgauge.h"
#include"measure_gauge_action.h"

extern int mixed_solve (spinor * const P, spinor * const Q, const int max_iter, 
                        double eps, const int rel_prec,const int N);
extern  int mixed_solve_eo (spinor * const P, spinor * const Q, const int max_iter, 
                            double eps, const int rel_prec, const int N);
#ifdef TEMPORALGAUGE
extern su3* g_trafo;
#endif
#endif

int invert_eo(spinor * const Even_new, spinor * const Odd_new, 
              spinor * const Even, spinor * const Odd,
              const double precision, const int max_iter,
              const int solver_flag, const int rel_prec,
              const int sub_evs_flag, const int even_odd_flag,
              const int no_extra_masses, double * const extra_masses, solver_params_t solver_params, const int id,
              const ExternalInverter external_inverter, const SloppyPrecision sloppy, const CompressionType compression )  {

  int iter = 0;

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
  if ( solver_flag==MG )
    return MG_solver_eo(Even_new, Odd_new, Even, Odd, precision, max_iter,
			rel_prec, VOLUME/2, g_gauge_field, &M_full);
#endif

  /* here comes the inversion using even/odd preconditioning */
  if(even_odd_flag) {
    if(g_proc_id == 0) {printf("# Using even/odd preconditioning!\n"); fflush(stdout);}
    
#ifdef HAVE_GPU
#ifdef TEMPORALGAUGE
    /* initialize temporal gauge here */
    int retval;
    double dret;
    double plaquette = 0.0;

    if(usegpu_flag){
    
      /* need VOLUME here (not N=VOLUME/2)*/
      if((retval=init_temporalgauge_trafo(VOLUME, g_gauge_field)) !=0){
        if(g_proc_id == 0) printf("Error while gauge fixing to temporal gauge. Aborting...\n");   
        exit(200);
      }
      plaquette = measure_plaquette(g_gauge_field);
      if(g_proc_id == 0) printf("Plaquette before gauge fixing: %.16e\n", plaquette/6./VOLUME);
      /* do trafo */
      apply_gtrafo(g_gauge_field, g_trafo);
      plaquette = measure_plaquette(g_gauge_field);
      if(g_proc_id == 0) printf("Plaquette after gauge fixing: %.16e\n", plaquette/6./VOLUME);
    
      /* do trafo to odd part of source */
      dret = square_norm(Odd, VOLUME/2 , 1);
      if(g_proc_id == 0) printf("square norm before gauge fixing: %.16e\n", dret); 
      apply_gtrafo_spinor_odd(Odd, g_trafo);
      dret = square_norm(Odd, VOLUME/2, 1);
      if(g_proc_id == 0) printf("square norm after gauge fixing: %.16e\n", dret);       
    
      /* do trafo to even part of source */
      dret = square_norm(Even, VOLUME/2 , 1);
      if(g_proc_id == 0) printf("square norm before gauge fixing: %.16e\n", dret); 
      apply_gtrafo_spinor_even(Even, g_trafo);
      dret = square_norm(Even, VOLUME/2, 1);
      if(g_proc_id == 0) printf("square norm after gauge fixing: %.16e\n", dret);      
    } 
#endif  
#endif /* HAVE_GPU*/    
    
 
    assign_mul_one_pm_imu_inv(Even_new, Even, +1., VOLUME/2);
    
    Hopping_Matrix(OE, g_spinor_field[DUM_DERI], Even_new);
    /* The sign is plus, since in Hopping_Matrix */
    /* the minus is missing                      */
    assign_mul_add_r(g_spinor_field[DUM_DERI], +1., Odd, VOLUME/2);
    /* Do the inversion with the preconditioned  */
    /* matrix to get the odd sites               */
    
#ifdef TM_USE_QPHIX
    if( external_inverter==QPHIX_INVERTER ) {
      // QPhiX inverts M(mu)M(mu)^dag or M(mu), no gamma_5 source multiplication required
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
    if(solver_flag == BICGSTAB) {
      if(g_proc_id == 0) {printf("# Using BiCGstab!\n"); fflush(stdout);}
      mul_one_pm_imu_inv(g_spinor_field[DUM_DERI], +1., VOLUME/2); 
      iter = bicgstab_complex(Odd_new, g_spinor_field[DUM_DERI], max_iter, precision, rel_prec, VOLUME/2, &Mtm_plus_sym_psi);
    }
    else if(solver_flag == BICG) {
      if(g_proc_id == 0) {printf("# Using BiCG!\n"); fflush(stdout);}
      mul_one_pm_imu_inv(g_spinor_field[DUM_DERI], +1., VOLUME/2);
      iter = bicg_complex(Odd_new, g_spinor_field[DUM_DERI], max_iter, precision, rel_prec, VOLUME/2, &Mtm_plus_sym_psi, &Mtm_plus_sym_dagg_psi);
    }

    else if(solver_flag == GMRES) {
      if(g_proc_id == 0) {printf("# Using GMRES! m = %d\n", gmres_m_parameter); fflush(stdout);}
      mul_one_pm_imu_inv(g_spinor_field[DUM_DERI], +1., VOLUME/2);
      iter = gmres(Odd_new, g_spinor_field[DUM_DERI], gmres_m_parameter, max_iter/gmres_m_parameter, precision, rel_prec, VOLUME/2, 1, &Mtm_plus_sym_psi);
    }
    else if(solver_flag == GCR) {
      if(g_proc_id == 0) {printf("# Using GCR! m = %d\n", gmres_m_parameter); fflush(stdout);}
      mul_one_pm_imu_inv(g_spinor_field[DUM_DERI], +1., VOLUME/2);
      iter = gcr(Odd_new, g_spinor_field[DUM_DERI], gmres_m_parameter, max_iter/gmres_m_parameter, precision, rel_prec, VOLUME/2, 0, &Mtm_plus_sym_psi);
    }
    else if (solver_flag == MCR) {
      if(g_proc_id == 0) {printf("# Using MCR! m = %d\n", gmres_m_parameter); fflush(stdout);}
      gamma5(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI], VOLUME/2);
      iter = mcr(Odd_new, g_spinor_field[DUM_DERI], gmres_m_parameter, max_iter/gmres_m_parameter, precision, rel_prec, VOLUME/2, 0, &Qtm_pm_psi);
      Qtm_minus_psi(Odd_new, Odd_new);
    }
    else if(solver_flag == GMRESDR) {
      if(g_proc_id == 0) {printf("# Using GMRES-DR! m = %d, NrEv = %d\n", 
                                 gmres_m_parameter, gmresdr_nr_ev); fflush(stdout);}
      mul_one_pm_imu_inv(g_spinor_field[DUM_DERI], +1., VOLUME/2);
      iter = gmres_dr(Odd_new, g_spinor_field[DUM_DERI], gmres_m_parameter, gmresdr_nr_ev, max_iter/gmres_m_parameter, precision, rel_prec, VOLUME/2, &Mtm_plus_sym_psi);
    }
    else if(solver_flag == FGMRES) {
      if(g_proc_id == 0) {printf("# Using FGMRES!\n"); fflush(stdout);}
      iter = fgmres(Odd_new, g_spinor_field[DUM_DERI], gmres_m_parameter, max_iter/gmres_m_parameter, precision, rel_prec, VOLUME/2, 0, &Qtm_pm_psi);
      gamma5(Odd_new, Odd_new, VOLUME/2);
      Qtm_minus_psi(Odd_new, Odd_new);
    }
    else if(solver_flag == BICGSTABELL) {
      if(g_proc_id == 0) {printf("# Using BiCGstab2!\n"); fflush(stdout);}
      mul_one_pm_imu_inv(g_spinor_field[DUM_DERI], +1., VOLUME/2); 
      iter = bicgstabell(Odd_new, g_spinor_field[DUM_DERI], max_iter, precision, rel_prec, 3, VOLUME/2, &Mtm_plus_sym_psi);
    }
    else if(solver_flag == PCG) {
      /* Here we invert the hermitean operator squared */
      gamma5(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI], VOLUME/2);  
      if(g_proc_id == 0) {printf("# Using PCG!\n"); fflush(stdout);}
      iter = pcg_her(Odd_new, g_spinor_field[DUM_DERI], max_iter, precision, rel_prec, VOLUME/2, &Qtm_pm_psi);
      Qtm_minus_psi(Odd_new, Odd_new);
    }
    else if(solver_flag == INCREIGCG) {
      /* Here we invert the hermitean operator squared */
      gamma5(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI], VOLUME/2);  
      if(g_proc_id == 0) {printf("# Using Incremental Eig-CG!\n"); fflush(stdout);}
      iter = incr_eigcg(VOLUME/2,solver_params.eigcg_nrhs,solver_params.eigcg_nrhs1, Odd_new, g_spinor_field[DUM_DERI], solver_params.eigcg_ldh, &Qtm_pm_psi,
                        solver_params.eigcg_tolsq1, solver_params.eigcg_tolsq, solver_params.eigcg_restolsq , solver_params.eigcg_rand_guess_opt,
                        rel_prec, max_iter, solver_params.eigcg_nev, solver_params.eigcg_vmax);
      Qtm_minus_psi(Odd_new, Odd_new);
    }
    else if(solver_flag == MIXEDCG) {
      /* Here we invert the hermitean operator squared */
      gamma5(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI], VOLUME/2);
      if(g_proc_id == 0) {printf("# Using Mixed Precision CG!\n"); fflush(stdout);}
      iter = mixed_cg_her(Odd_new, g_spinor_field[DUM_DERI], solver_params, max_iter, precision, rel_prec, 
                          VOLUME/2, &Qtm_pm_psi, &Qtm_pm_psi_32);
      Qtm_minus_psi(Odd_new, Odd_new);
    }
    else if(solver_flag == RGMIXEDCG) {
      /* Here we invert the hermitean operator squared */
      gamma5(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI], VOLUME/2);
      if(g_proc_id == 0) {printf("# Using Mixed Precision CG!\n"); fflush(stdout);}
      iter = rg_mixed_cg_her(Odd_new, g_spinor_field[DUM_DERI], solver_params, max_iter, precision, rel_prec,
                             VOLUME/2, &Qtm_pm_psi, &Qtm_pm_psi_32);
      Qtm_minus_psi(Odd_new, Odd_new);
    }
    else if(solver_flag == CG) {
      /* Here we invert the hermitean operator squared */
      gamma5(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI], VOLUME/2);
      if(g_proc_id == 0) {
        printf("# Using CG!\n"); 
        printf("# mu = %.12f, kappa = %.12f\n", g_mu/2./g_kappa, g_kappa);
        fflush(stdout);
      }
#ifdef HAVE_GPU
      if(usegpu_flag){
        if(g_proc_id == 0) printf("Using GPU for inversion\n");
        iter = mixed_solve_eo(Odd_new, g_spinor_field[DUM_DERI], max_iter,   precision, rel_prec, VOLUME/2);
      }
      else {
        iter = cg_her(Odd_new, g_spinor_field[DUM_DERI], max_iter, precision, rel_prec, VOLUME/2, &Qtm_pm_psi);
        Qtm_minus_psi(Odd_new, Odd_new);
      }
#else        
      iter = cg_her(Odd_new, g_spinor_field[DUM_DERI], max_iter, precision, rel_prec, 
                    VOLUME/2, &Qtm_pm_psi);
      Qtm_minus_psi(Odd_new, Odd_new);
#endif /*HAVE_GPU*/
    }
    else if(solver_flag == MR) {
      if(g_proc_id == 0) {printf("# Using MR!\n"); fflush(stdout);}
      iter = mr(Odd_new, g_spinor_field[DUM_DERI], max_iter, precision, rel_prec, VOLUME/2, 1, &Mtm_plus_psi);
    }
    else if(solver_flag == CGS) {
      if(g_proc_id == 0) {printf("# Using CGS!\n"); fflush(stdout);}
      mul_one_pm_imu_inv(g_spinor_field[DUM_DERI], +1., VOLUME/2); 
      iter = cgs_real(Odd_new, g_spinor_field[DUM_DERI], max_iter, precision, rel_prec, VOLUME/2, &Mtm_plus_sym_psi);
    }
    else {
      if(g_proc_id == 0) {printf("# Using CG as default solver!\n"); fflush(stdout);}
      gamma5(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI], VOLUME/2);
#ifdef HAVE_GPU
      if(g_proc_id == 0) {printf("Using GPU for inversion\n");
        fflush(stdout);}
      iter = mixed_solve_eo(Odd_new, g_spinor_field[DUM_DERI], max_iter,   precision, rel_prec, VOLUME/2);
#else
      iter = cg_her(Odd_new, g_spinor_field[DUM_DERI], max_iter, precision, rel_prec, VOLUME/2, &Qtm_pm_psi);
      Qtm_minus_psi(Odd_new, Odd_new);
#endif
    }
    
    /* In case of failure, redo with CG */
    if(iter == -1 && solver_flag !=CG) {
      /* Here we invert the hermitean operator squared */
      mul_one_pm_imu(g_spinor_field[DUM_DERI], +1.); 
      gamma5(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI], VOLUME/2);  
      if(g_proc_id == 0) {printf("# Redoing it with CG!\n"); fflush(stdout);}
      iter = cg_her(Odd_new, g_spinor_field[DUM_DERI], max_iter, precision, rel_prec, VOLUME/2, &Qtm_pm_psi);
      Qtm_minus_psi(Odd_new, Odd_new);
    }
    
    /* Reconstruct the even sites                */
    Hopping_Matrix(EO, g_spinor_field[DUM_DERI], Odd_new);
    mul_one_pm_imu_inv(g_spinor_field[DUM_DERI], +1., VOLUME/2);
    /* The sign is plus, since in Hopping_Matrix */
    /* the minus is missing                      */
    assign_add_mul_r(Even_new, g_spinor_field[DUM_DERI], +1., VOLUME/2);
 
#ifdef HAVE_GPU  
    /* return from temporal gauge again */
#ifdef TEMPORALGAUGE
    if(usegpu_flag){ 
      plaquette = measure_plaquette(g_gauge_field);
      if(g_proc_id == 0) printf("Plaquette before inverse gauge fixing: %.16e\n", plaquette/6./VOLUME);
    
      /* undo trafo */
    
      /*apply_inv_gtrafo(g_gauge_field, g_trafo);*/
      /* copy back the saved original field located in g_tempgauge_field -> update necessary*/
      copy_gauge_field(g_gauge_field, g_tempgauge_field);
      g_update_gauge_copy = 1;
    
    
      plaquette = measure_plaquette(g_gauge_field);
      if(g_proc_id == 0) printf("Plaquette after inverse gauge fixing: %.16e\n", plaquette/6./VOLUME);
   
      /* undo trafo to source (Even, Odd) */
      dret = square_norm(Even, VOLUME/2 , 1);
      if(g_proc_id == 0) printf("square norm before gauge fixing: %.16e\n", dret); 
      apply_inv_gtrafo_spinor_even(Even, g_trafo);
      dret = square_norm(Even, VOLUME/2, 1);
      if(g_proc_id == 0) printf("square norm after gauge fixing: %.16e\n", dret);  
      dret = square_norm(Odd, VOLUME/2 , 1);
      if(g_proc_id == 0) printf("square norm before gauge fixing: %.16e\n", dret); 
      apply_inv_gtrafo_spinor_odd(Odd, g_trafo);
      dret = square_norm(Odd, VOLUME/2, 1);
      if(g_proc_id == 0) printf("square norm after gauge fixing: %.16e\n", dret); 
    
    
      dret = square_norm(Even_new, VOLUME/2 , 1);
      if(g_proc_id == 0) printf("square norm before gauge fixing: %.16e\n", dret); 
      apply_inv_gtrafo_spinor_even(Even_new, g_trafo);
      dret = square_norm(Even_new, VOLUME/2, 1);
      if(g_proc_id == 0) printf("square norm after gauge fixing: %.16e\n", dret);  

      dret = square_norm(Odd_new, VOLUME/2 , 1);
      if(g_proc_id == 0) printf("square norm before gauge fixing: %.16e\n", dret); 
      apply_inv_gtrafo_spinor_odd(Odd_new, g_trafo);
      dret = square_norm(Odd_new, VOLUME/2, 1);
      if(g_proc_id == 0) printf("square norm after gauge fixing: %.16e\n", dret); 
  
    
      finalize_temporalgauge();
    }
#endif
#endif     
  
  
  }

  else {
    /* here comes the inversion not using even/odd preconditioning */
    if(g_proc_id == 0) {printf("# Not using even/odd preconditioning!\n"); fflush(stdout);}
    convert_eo_to_lexic(g_spinor_field[DUM_DERI], Even, Odd);
    convert_eo_to_lexic(g_spinor_field[DUM_DERI+1], Even_new, Odd_new);
    
    if(solver_flag == BICGSTAB) {
      if(g_proc_id == 0) {printf("# Using BiCGstab!\n"); fflush(stdout);}
      if(use_preconditioning==1 && g_precWS!=NULL){
        if(g_proc_id == 0) {printf("# Using preconditioning (which one?)!\n");}
        iter = bicgstab_complex(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI], max_iter, precision, rel_prec, VOLUME, &D_psi_prec);
      } else {
        if(g_proc_id == 0) {printf("# Not using preconditioning (which one?)!\n");}
        iter = bicgstab_complex(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI], max_iter, precision, rel_prec, VOLUME, &D_psi);
      }
    }
    else if(solver_flag == BICG) {
      if(g_proc_id == 0) {printf("# Using BiCG!\n"); fflush(stdout);}
      if(use_preconditioning==1 && g_precWS!=NULL){
        //if(g_proc_id == 0) {printf("# Using preconditioning (which one?)!\n");}
        //iter = bicg_complex(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI], max_iter, precision, rel_prec, VOLUME, &D_psi_prec, &D_dagg_psi_prec);
        
        if(g_proc_id == 0) {printf("# Not using preconditioning (which one?)!\n");}
        iter = bicg_complex(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI], max_iter, precision, rel_prec, VOLUME, &D_psi, &D_dagg_psi);
      } else {
        if(g_proc_id == 0) {printf("# Not using preconditioning (which one?)!\n");}
        iter = bicg_complex(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI], max_iter, precision, rel_prec, VOLUME, &D_psi, &D_dagg_psi);
      }
    }
    else if(solver_flag == CGS) {
      if(g_proc_id == 0) {printf("# Using CGS!\n"); fflush(stdout);}

      if(use_preconditioning==1 && g_precWS!=NULL){
        if(g_proc_id == 0) {printf("# Using preconditioning (which one?)!\n");}
        iter = cgs_real(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI], max_iter, precision, rel_prec, VOLUME, &D_psi_prec);
      } else {
        if(g_proc_id == 0) {printf("# Not using preconditioning (which one?)!\n");}
        iter = cgs_real(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI], max_iter, precision, rel_prec, VOLUME, &D_psi);
      }


    }    
    else if(solver_flag == GMRES) {
      if(g_proc_id == 0) {printf("# Using GMRES! m = %d\n", gmres_m_parameter); fflush(stdout);}

      if(use_preconditioning==1 && g_precWS!=NULL){
        if(g_proc_id == 0) {printf("# Using preconditioning (which one?)!\n");}
        iter = gmres(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI], gmres_m_parameter, max_iter/gmres_m_parameter, precision, rel_prec, VOLUME, 1, &D_psi_prec);
      } else {
        if(g_proc_id == 0) {printf("# not using preconditioning (which one?)!\n");}
        iter = gmres(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI], gmres_m_parameter, max_iter/gmres_m_parameter, precision, rel_prec, VOLUME, 1, &D_psi);
      }
    }
    else if(solver_flag == MIXEDCG) {
      if(g_proc_id == 0) {printf("# Using MIXEDCG!\n"); fflush(stdout);}
      gamma5(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI], VOLUME);
      iter = mixed_cg_her(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+1], solver_params, max_iter, 
                          precision, rel_prec, VOLUME, &Q_pm_psi, &Q_pm_psi_32);
      Q_minus_psi(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI]);
    } else if(solver_flag == RGMIXEDCG) {
      if(g_proc_id == 0) {printf("# Using MIXEDCG!\n"); fflush(stdout);}
      gamma5(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI], VOLUME);
      iter = rg_mixed_cg_her(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+1], solver_params, max_iter, 
                             precision, rel_prec, VOLUME, &Q_pm_psi, &Q_pm_psi_32);
      Q_minus_psi(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI]);
    } else if(solver_flag == FGMRES) {
      if(g_proc_id == 0) {printf("# Using FGMRES! m = %d\n", gmres_m_parameter); fflush(stdout);}
      iter = fgmres(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI], gmres_m_parameter, max_iter/gmres_m_parameter, precision, rel_prec, VOLUME, Msap_precon, &D_psi); 
    }
    else if(solver_flag == GCR) {
      if(g_proc_id == 0) {printf("# Using GCR! m = %d\n", gmres_m_parameter); fflush(stdout);}
      iter = gcr(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI], gmres_m_parameter, max_iter/gmres_m_parameter, precision, rel_prec, VOLUME, Msap_precon, &D_psi); 
    }
    else if(solver_flag == MCR) {
      if(g_proc_id == 0) {printf("# Using mCR! m = %d\n", gmres_m_parameter); fflush(stdout);}
      /* iter = mcr(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI], gmres_m_parameter, max_iter/gmres_m_parameter, precision, rel_prec, VOLUME, 0, &D_psi); */
      gamma5(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI], VOLUME);
      iter = mcr(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI], gmres_m_parameter, max_iter/gmres_m_parameter, precision, rel_prec, VOLUME, 0, &Q_plus_psi);
    }
    else if(solver_flag == CR) {
      if(g_proc_id == 0) {printf("# Using CR and iQ!\n"); fflush(stdout);}
      gamma5(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI], VOLUME);
      iter = cr(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI], gmres_m_parameter, max_iter/gmres_m_parameter, precision, rel_prec, VOLUME, 0, &Q_plus_psi);
      /* Solve DdaggD */
      /* gamma5(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI], VOLUME);
         iter = cr(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+1], gmres_m_parameter, max_iter/gmres_m_parameter, precision, rel_prec, VOLUME, 0, &Q_pm_psi);
         Q_minus_psi(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI]);
      */  
    }
    else if (solver_flag == DFLGCR) {
      if(g_proc_id == 0) {printf("# Using deflated GCR solver! m = %d\n", gmres_m_parameter); fflush(stdout);}
      iter = gcr(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI], gmres_m_parameter, 
                 max_iter/gmres_m_parameter, precision, rel_prec, VOLUME, 2, &D_psi);
    }
    else if (solver_flag == DFLFGMRES) {
      if(g_proc_id == 0) {printf("# Using deflated FGMRES solver! m = %d\n", gmres_m_parameter); fflush(stdout);}
      iter = fgmres(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI], gmres_m_parameter, 
                    max_iter/gmres_m_parameter, precision, rel_prec, VOLUME, 2, &D_psi);
    }
    else if (solver_flag == CGMMS) {
      /* FIXME temporary workaround for the multiple masses interface */
      double * shifts = (double*)calloc(no_extra_masses+1,sizeof(double));
      shifts[0]=g_mu;
      for(int i = 0; i < no_extra_masses; ++i)
        shifts[i+1] = extra_masses[i];
      g_mu = 0;
      solver_params_t solver_params;
      solver_params.shifts = shifts;
      solver_params.no_shifts = no_extra_masses+1;
      solver_params.rel_prec = rel_prec;
      solver_params.max_iter = max_iter;
      solver_params.squared_solver_prec = precision;
      solver_params.sdim = VOLUME;
      solver_params.M_psi = &Q_pm_psi;
      solver_params.type = solver_flag; 

      /* FIXME temporary workaround for the multiple shift solver interface and integration of IO */
      spinor * P_memory;
      spinor ** P;
      allocate_spinor_field_array(&P,&P_memory,VOLUME,no_extra_masses+1);

      if(g_proc_id == 0) {printf("# Using multi mass CG!\n"); fflush(stdout);}
      
      gamma5(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI], VOLUME);
      iter = cg_mms_tm(P, g_spinor_field[DUM_DERI+1],&solver_params,&cgmms_reached_prec);
      g_mu = shifts[0];
      Q_minus_psi(g_spinor_field[DUM_DERI+1], P[0]);
      
      cgmms_write_props(P,shifts,no_extra_masses+1,id,iter);
    
      free_spinor_field_array(&P_memory);
      free(P);
      free(shifts);
    }
    else if(solver_flag == PCG) {
      if(g_proc_id == 0) {printf("# Using PCG!\n"); fflush(stdout);}
      gamma5(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI], VOLUME);
      iter = pcg_her(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+1], max_iter, precision, 
                     rel_prec, VOLUME, &Q_pm_psi);
      Q_minus_psi(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI]);
    }
    else {
      if(g_proc_id == 0) {printf("# Using CG!\n"); fflush(stdout);}
#ifdef HAVE_GPU 
      if(usegpu_flag){
        if(g_proc_id == 0) printf("# Using GPU for inversion\n");
        iter = mixed_solve(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI], max_iter, precision, rel_prec, VOLUME);
      }
      else{
        gamma5(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI], VOLUME);
        iter = cg_her(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+1], max_iter, precision, 
                      rel_prec, VOLUME, &Q_pm_psi);
        Q_minus_psi(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI]);
      }
#else
      gamma5(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI], VOLUME);

      if(use_preconditioning==1 && g_precWS!=NULL){
        spinorPrecWS *ws=(spinorPrecWS*)g_precWS;
        static _Complex double alpha = 0.0;
        if(g_proc_id==0) {printf("# Using preconditioning (which one?)!\n");}

        if(g_prec_sequence_d_dagger_d[2] != 0.0){
          alpha = g_prec_sequence_d_dagger_d[2];
          spinorPrecondition(g_spinor_field[DUM_DERI+1],g_spinor_field[DUM_DERI+1],ws,T,L,alpha,0,1);
        }

        iter = cg_her(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+1], max_iter, precision, 
                      rel_prec, VOLUME, &Q_pm_psi_prec);

        if(g_prec_sequence_d_dagger_d[0] != 0.0){
          alpha = g_prec_sequence_d_dagger_d[0];
          spinorPrecondition(g_spinor_field[DUM_DERI],g_spinor_field[DUM_DERI],ws,T,L,alpha,0,1);
        }

      } else {
        if(g_proc_id==0) {printf("# Not using preconditioning!\n");}
        iter = cg_her(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+1], max_iter, precision, 
                      rel_prec, VOLUME, &Q_pm_psi);
      }


      Q_minus_psi(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI]);

      if(use_preconditioning==1 && g_precWS!=NULL){
        spinorPrecWS *ws=(spinorPrecWS*)g_precWS;
        static _Complex double alpha = 0.0;
        if(g_prec_sequence_d_dagger_d[1] != 0.0){
          alpha = g_prec_sequence_d_dagger_d[1];
          spinorPrecondition(g_spinor_field[DUM_DERI+1],g_spinor_field[DUM_DERI+1],ws,T,L,alpha,0,1);
        }
      }
#endif
    }
    convert_lexic_to_eo(Even_new, Odd_new, g_spinor_field[DUM_DERI+1]);
  }
  return(iter);
}

/* FIXME temporary solution for the writing of CGMMS propagators until the input/output interface for
   invert_eo has been generalized
   NOTE that no_shifts = no_extra_masses+1 */
static void cgmms_write_props(spinor ** const P, double const * const shifts, const int no_shifts, const int id, const int iteration) {
  int append = 0;
  char filename[300];
  WRITER * writer = NULL;
  paramsInverterInfo *inverterInfo = NULL;
  paramsPropagatorFormat *propagatorFormat = NULL;
  
  spinor * temp_eo_spinors_memory;
  spinor ** temp_eo_spinors;
  
  allocate_spinor_field_array(&temp_eo_spinors, &temp_eo_spinors_memory, VOLUME/2, 2);

  /* save all the results of (Q^dagger Q)^(-1) \gamma_5 \phi */
  for(int im = 0; im < no_shifts; im++) {
    if(SourceInfo.type != SRC_TYPE_VOL) {
      if (PropInfo.splitted) {
        if(T_global > 99) sprintf(filename, "%s.%.2d.%.4d.%.3d.%.2d.cgmms.%.2d.inverted", PropInfo.basename, id, SourceInfo.nstore, SourceInfo.t, SourceInfo.ix, im);
        else sprintf(filename, "%s.%.2d.%.4d.%.2d.%.2d.cgmms.%.2d.inverted", PropInfo.basename, id, SourceInfo.nstore, SourceInfo.t, SourceInfo.ix, im);
      } else {
        sprintf(filename, "%s.%.2d.%.4d.%.2d.cgmms.%.2d.inverted", PropInfo.basename, id, SourceInfo.nstore, SourceInfo.t, im);
      }
    } else {
      sprintf(filename, "%s.%.2d.%.4d.%.5d.cgmms.%.2d.0", PropInfo.basename, id, SourceInfo.nstore, SourceInfo.sample, im);
    }
    
    if(g_kappa != 0) {
      mul_r(P[im], (2*g_kappa)*(2*g_kappa), P[im], VOLUME);
    }

    append = !PropInfo.splitted;

    construct_writer(&writer, filename, append);

    if (PropInfo.splitted || SourceInfo.ix == index_start) {
      //Create the inverter info NOTE: always set to TWILSON=12 and 1 flavour (to be adjusted)
      inverterInfo = construct_paramsInverterInfo(cgmms_reached_prec, iteration+1, 12, 1);
      inverterInfo->cgmms_mass = shifts[im]/(2 * inverterInfo->kappa);
      write_spinor_info(writer, PropInfo.format, inverterInfo, append);
      //Create the propagatorFormat NOTE: always set to 1 flavour (to be adjusted)
      propagatorFormat = construct_paramsPropagatorFormat(PropInfo.precision, 1);
      write_propagator_format(writer, propagatorFormat);
      free(inverterInfo);
      free(propagatorFormat);
    }
    convert_lexic_to_eo(temp_eo_spinors[1], temp_eo_spinors[0], P[im]);
    write_spinor(writer, &temp_eo_spinors[1], &temp_eo_spinors[0], 1, PropInfo.precision);
    destruct_writer(writer);
  }
  free_spinor_field_array(&temp_eo_spinors_memory);
  free(temp_eo_spinors);
}


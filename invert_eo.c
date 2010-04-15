/***********************************************************************
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008 Carsten Urbach
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
 ***********************************************************************/
/* $Id$ */

/****************************************************************
 *
 * invert_eo makes an inversion with EO precoditioned
 * tm Operator
 *
 * Even and Odd are the numbers of spinor_field that contain
 * the even and the odd sites of the source. The result is stored
 * int Even_new and Odd_new.
 *
 * invert_eo returns the number of iterations neede or -1 if the 
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
#include"tm_operators.h"
#include"Hopping_Matrix.h"
#include"D_psi.h"
#include"linsolve.h"
#include"gamma.h"
#include"solver/solver.h"
#include"read_input.h"
#include"xchange.h"
#include"solver/poly_precon.h"
#include"solver/dfl_projector.h"
#include"invert_eo.h"


#ifdef HAVE_GPU
extern int mixed_solve (spinor * const P, spinor * const Q, const int max_iter, 
           double eps, const int rel_prec,const int N);
extern  int mixed_solve_eo (spinor * const P, spinor * const Q, const int max_iter, 
           double eps, const int rel_prec, const int N);
#endif



int invert_eo(spinor * const Even_new, spinor * const Odd_new, 
	      spinor * const Even, spinor * const Odd,
	      const double precision, const int max_iter,
	      const int solver_flag, const int rel_prec,
	      const int sub_evs_flag, const int even_odd_flag) {

  int iter = 0;
  /* here comes the inversion using even/odd preconditioning */
  if(even_odd_flag) {
    if(g_proc_id == 0) {printf("# Using Even/Odd preconditioning!\n"); fflush(stdout);}
    assign_mul_one_pm_imu_inv(Even_new, Even, +1.);
    
    Hopping_Matrix(OE, g_spinor_field[DUM_DERI], Even_new); 
    /* The sign is plus, since in Hopping_Matrix */
    /* the minus is missing                      */
    assign_mul_add_r(g_spinor_field[DUM_DERI], +1., Odd, VOLUME/2);
    /* Do the inversion with the preconditioned  */
    /* matrix to get the odd sites               */
    
    if(solver_flag == BICGSTAB) {
      if(g_proc_id == 0) {printf("# Using BiCGstab!\n"); fflush(stdout);}
      mul_one_pm_imu_inv(g_spinor_field[DUM_DERI], +1.); 
      iter = bicgstab_complex(Odd_new, g_spinor_field[DUM_DERI], max_iter, precision, rel_prec, VOLUME/2, &Mtm_plus_sym_psi);
    }
    else if(solver_flag == GMRES) {
      if(g_proc_id == 0) {printf("# Using GMRES! m = %d\n", gmres_m_parameter); fflush(stdout);}
      mul_one_pm_imu_inv(g_spinor_field[DUM_DERI], +1.);
      iter = gmres(Odd_new, g_spinor_field[DUM_DERI], gmres_m_parameter, max_iter/gmres_m_parameter, precision, rel_prec, VOLUME/2, 1, &Mtm_plus_sym_psi);
    }
    else if(solver_flag == GCR) {
      if(g_proc_id == 0) {printf("# Using GCR! m = %d\n", gmres_m_parameter); fflush(stdout);}
      mul_one_pm_imu_inv(g_spinor_field[DUM_DERI], +1.);
      iter = gcr(Odd_new, g_spinor_field[DUM_DERI], gmres_m_parameter, max_iter/gmres_m_parameter, precision, rel_prec, VOLUME/2, 0, &Mtm_plus_sym_psi);
    }
    else if(solver_flag == GMRESDR) {
      if(g_proc_id == 0) {printf("# Using GMRES-DR! m = %d, NrEv = %d\n", 
				 gmres_m_parameter, gmresdr_nr_ev); fflush(stdout);}
      mul_one_pm_imu_inv(g_spinor_field[DUM_DERI], +1.);
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
      mul_one_pm_imu_inv(g_spinor_field[DUM_DERI], +1.); 
      iter = bicgstabell(Odd_new, g_spinor_field[DUM_DERI], max_iter, precision, rel_prec, 3, VOLUME/2, &Mtm_plus_sym_psi);
    }
    else if(solver_flag == PCG) {
      /* Here we invert the hermitean operator squared */
      gamma5(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI], VOLUME/2);  
      if(g_proc_id == 0) {printf("# Using PCG!\n"); fflush(stdout);}
      iter = pcg_her(Odd_new, g_spinor_field[DUM_DERI], max_iter, precision, rel_prec, VOLUME/2, &Qtm_pm_psi);
      Qtm_minus_psi(Odd_new, Odd_new);
    }
    else if(solver_flag == MIXEDCG) {
      /* Here we invert the hermitean operator squared */
      gamma5(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI], VOLUME/2);
      if(g_proc_id == 0) {printf("# Using Mixed Precision CG!\n"); fflush(stdout);}
      iter = mixed_cg_her(Odd_new, g_spinor_field[DUM_DERI], max_iter, precision, rel_prec, 
			  VOLUME/2, &Qtm_pm_psi);
      Qtm_minus_psi(Odd_new, Odd_new);
    }
    else if(solver_flag == CG) {
      /* Here we invert the hermitean operator squared */
      gamma5(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI], VOLUME/2);
      if(g_proc_id == 0) {printf("# Using CG!\n"); fflush(stdout);}
#ifdef HAVE_GPU  
      if(g_proc_id == 0) {printf("Using GPU for inversion\n");
      fflush(stdout);}
      iter = mixed_solve_eo(Odd_new, g_spinor_field[DUM_DERI], max_iter,   precision, rel_prec, VOLUME/2); 
#else        
      iter = cg_her(Odd_new, g_spinor_field[DUM_DERI], max_iter, precision, rel_prec, 
		    VOLUME/2, &Qtm_pm_psi, sub_evs_flag, 1000);
      Qtm_minus_psi(Odd_new, Odd_new);
#endif
    }
    else if(solver_flag == MR) {
      if(g_proc_id == 0) {printf("# Using MR!\n"); fflush(stdout);}
      iter = mr(Odd_new, g_spinor_field[DUM_DERI], max_iter, precision, rel_prec, VOLUME/2, 1, &Mtm_plus_psi);
    }
    else if(solver_flag == CGS) {
      if(g_proc_id == 0) {printf("# Using CGS!\n"); fflush(stdout);}
      mul_one_pm_imu_inv(g_spinor_field[DUM_DERI], +1.); 
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
      iter = cg_her(Odd_new, g_spinor_field[DUM_DERI], max_iter, precision, rel_prec, VOLUME/2, &Qtm_pm_psi, 0, 0);
      Qtm_minus_psi(Odd_new, Odd_new);
#endif
    }
    
    /* In case of failure, redo with CG */
    if(iter == -1 && solver_flag !=CG) {
      /* Here we invert the hermitean operator squared */
      mul_one_pm_imu(g_spinor_field[DUM_DERI], +1.); 
      gamma5(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI], VOLUME/2);  
      if(g_proc_id == 0) {printf("# Redoing it with CG!\n"); fflush(stdout);}
      iter = cg_her(Odd_new, g_spinor_field[DUM_DERI], max_iter, precision, rel_prec, VOLUME/2, &Qtm_pm_psi, 0, 0.);
      Qtm_minus_psi(Odd_new, Odd_new);
    }
    
    /* Reconstruct the even sites                */
    Hopping_Matrix(EO, g_spinor_field[DUM_DERI], Odd_new);
    mul_one_pm_imu_inv(g_spinor_field[DUM_DERI], +1.);
    /* The sign is plus, since in Hopping_Matrix */
    /* the minus is missing                      */
    assign_add_mul_r(Even_new, g_spinor_field[DUM_DERI], +1., VOLUME/2);
  }

  else {
    /* here comes the inversion not using even/odd preconditioning */
    if(g_proc_id == 0) {printf("# Not using Even/Odd preconditioning!\n"); fflush(stdout);}
    convert_eo_to_lexic(g_spinor_field[DUM_DERI], Even, Odd);
    convert_eo_to_lexic(g_spinor_field[DUM_DERI+1], Even_new, Odd_new);
    
    if(solver_flag == BICGSTAB) {
      if(g_proc_id == 0) {printf("# Using BiCGstab!\n"); fflush(stdout);}
      iter = bicgstab_complex(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI], max_iter, precision, rel_prec, VOLUME, &D_psi);
    }
    else if(solver_flag == CGS) {
      if(g_proc_id == 0) {printf("# Using CGS!\n"); fflush(stdout);}
      iter = cgs_real(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI], max_iter, precision, rel_prec, VOLUME, &D_psi);
    }    
    else if(solver_flag == GMRES) {
      if(g_proc_id == 0) {printf("# Using GMRES! m = %d\n", gmres_m_parameter); fflush(stdout);}
      iter = gmres(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI], gmres_m_parameter, max_iter/gmres_m_parameter, precision, rel_prec, VOLUME, 1, &D_psi);
    }
    else if(solver_flag == FGMRES) {
      if(g_proc_id == 0) {printf("# Using FGMRES! m = %d\n", gmres_m_parameter); fflush(stdout);}
      iter = fgmres(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI], gmres_m_parameter, max_iter/gmres_m_parameter, precision, rel_prec, VOLUME, 1, &D_psi); 
/*       gamma5(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI], VOLUME); */
/*       iter = fgmres(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+1], gmres_m_parameter, max_iter/gmres_m_parameter, precision, rel_prec, VOLUME, &Q_pm_psi);  */
/*       Q_minus_psi(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI]); */
    }
    else if(solver_flag == GCR) {
      if(g_proc_id == 0) {printf("# Using GCR! m = %d\n", gmres_m_parameter); fflush(stdout);}
      iter = gcr(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI], gmres_m_parameter, max_iter/gmres_m_parameter, precision, rel_prec, VOLUME, 1, &D_psi); 
/*       gamma5(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI], VOLUME); */
/*       iter = gcr(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+1], gmres_m_parameter, max_iter/gmres_m_parameter, precision, rel_prec, VOLUME, &Q_pm_psi); */
/*       Q_minus_psi(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI]); */
    }
    else if(solver_flag == DFLGCR || solver_flag == DFLFGMRES) {
      if(g_proc_id == 0) {printf("# Using deflated solver! m = %d\n", gmres_m_parameter); fflush(stdout);}
      /* apply P_L to source           */
      project_left(g_spinor_field[DUM_DERI+2], g_spinor_field[DUM_DERI]);
      if(g_proc_id == 0) printf("Applied P_L to source\n");
      /* invert P_L D on source -> chi */
      if(solver_flag == DFLGCR) {
	iter = gcr(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI+2], gmres_m_parameter, 
		   max_iter/gmres_m_parameter, precision, rel_prec, VOLUME, 1, &project_left_D);
      }
      else {
	iter = fgmres(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI+2], gmres_m_parameter, 
		      max_iter/gmres_m_parameter, precision, rel_prec, VOLUME, 1, &project_left_D);
      }
      /* apply P_R to chi              */
      project_right(g_spinor_field[DUM_DERI+2], g_spinor_field[DUM_DERI+1]);
      if(g_proc_id == 0) printf("Applied P_R to solution\n");
      /* reconstruct solution          */
      project(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI]);
      add(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI+2], VOLUME);
    }
    else if (solver_flag == CGMMS) {
      if(g_proc_id == 0) {printf("# Using multi mass CG!\n"); fflush(stdout);}
      gamma5(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI], VOLUME);
      iter = cg_mms_tm(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+1], 
		       max_iter, precision, rel_prec, VOLUME, &Q_pm_psi);
      Q_minus_psi(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI]);
    }
    else {
      if(g_proc_id == 0) {printf("# Using CG!\n"); fflush(stdout);}
#ifdef HAVE_GPU
      if(g_proc_id == 0) {printf("Using GPU for inversion\n");
      fflush(stdout);}      
      iter = mixed_solve(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI], max_iter, precision, rel_prec, VOLUME);      
#else
      gamma5(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI], VOLUME);
      iter = cg_her(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+1], max_iter, precision, rel_prec, VOLUME, &Q_pm_psi, 0, 0);
      Q_minus_psi(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI]);
#endif      
    }
    convert_lexic_to_eo(Even_new, Odd_new, g_spinor_field[DUM_DERI+1]);
  }
  return(iter);
}

void M_full(spinor * const Even_new, spinor * const Odd_new, 
	    spinor * const Even, spinor * const Odd) {
  /* Even sites */
  Hopping_Matrix(EO, g_spinor_field[DUM_DERI], Odd);
  assign_mul_one_pm_imu(Even_new, Even, 1.); 
  assign_add_mul_r(Even_new, g_spinor_field[DUM_DERI], -1., VOLUME/2);

  /* Odd sites */
  Hopping_Matrix(OE, g_spinor_field[DUM_DERI], Even);
  assign_mul_one_pm_imu(Odd_new, Odd, 1.); 
  assign_add_mul_r(Odd_new, g_spinor_field[DUM_DERI], -1., VOLUME/2);
}

void Q_full(spinor * const Even_new, spinor * const Odd_new, 
	    spinor * const Even, spinor * const Odd) {
  /* Even sites */
  Hopping_Matrix(EO, g_spinor_field[DUM_DERI], Odd);
  assign_mul_one_pm_imu(Even_new, Even, 1.); 
  assign_add_mul_r(Even_new, g_spinor_field[DUM_DERI], -1., VOLUME/2);
  gamma5(Even_new, Even_new, VOLUME/2);

  /* Odd sites */
  Hopping_Matrix(OE, g_spinor_field[DUM_DERI], Even);
  assign_mul_one_pm_imu(Odd_new, Odd, 1.); 
  assign_add_mul_r(Odd_new, g_spinor_field[DUM_DERI], -1., VOLUME/2);
  gamma5(Odd_new, Odd_new, VOLUME/2);
}

void M_minus_1_timesC(spinor * const Even_new, spinor * const Odd_new, 
		      spinor * const Even, spinor * const Odd) {
  /* Even sites */
  Hopping_Matrix(EO, Even_new, Odd);
  mul_one_pm_imu_inv(Even_new, 1.); 

  /* Odd sites */
  Hopping_Matrix(OE, Odd_new, Even);
  mul_one_pm_imu_inv(Odd_new, 1.); 
}

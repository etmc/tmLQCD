/***********************************************************************
 *
 * Copyright (C) 2014 Florian Burger
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
 * File: monomial_solve.c
 *
 * solver wrapper for monomials
 *
 * The externally accessible functions are
 *
 *
 *   int solve_degenerate(spinor * const P, spinor * const Q, const int max_iter, 
           double eps_sq, const int rel_prec, const int N, matrix_mult f)
 *   int solve_mms_nd(spinor ** const Pup, spinor ** const Pdn, 
 *                    spinor * const Qup, spinor * const Qdn, 
 *                    solver_pm_t * solver_pm)  
 *
 **************************************************************************/


#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include "global.h"
#include "read_input.h"
#include "solver/solver.h"
#include "solver/matrix_mult_typedef.h"
#include "solver/solver_types.h"
#include "solver/solver_params.h"
#include "operator/tm_operators.h"
#include "operator/tm_operators_32.h"
#include "operator/tm_operators_nd.h"
#include "operator/tm_operators_nd_32.h"
#include "operator/clovertm_operators.h"
#include "operator/clovertm_operators_32.h"
#include "monomial_solve.h"
#ifdef DDalphaAMG
#include "DDalphaAMG_interface.h"
#endif


#ifdef HAVE_GPU
#include"../GPU/cudadefs.h"
extern  int linsolve_eo_gpu (spinor * const P, spinor * const Q, const int max_iter, 
                            double eps, const int rel_prec, const int N, matrix_mult f);
extern int dev_cg_mms_tm_nd(spinor ** const Pup, spinor ** const Pdn, 
		 spinor * const Qup, spinor * const Qdn, 
		 solver_pm_t * solver_pm);
   #ifdef TEMPORALGAUGE
     #include "../temporalgauge.h" 
   #endif
#include "read_input.h" 
#endif

int solve_degenerate(spinor * const P, spinor * const Q, solver_params_t solver_params,
                     const int max_iter, double eps_sq, const int rel_prec, 
                     const int N, matrix_mult f, int solver_type){
  int iteration_count = 0;
  int use_solver = solver_type;
  
  if(use_solver == MIXEDCG || use_solver == RGMIXEDCG){
    // the default mixed solver is rg_mixed_cg_her
    int (*msolver_fp)(spinor * const, spinor * const, solver_params_t, 
                      const int, double, const int, const int, matrix_mult, matrix_mult32) = rg_mixed_cg_her;

    // but it might be necessary at some point to use the old version
    if(use_solver == MIXEDCG){
      msolver_fp = mixed_cg_her;
    }

    if(usegpu_flag){   
      #ifdef HAVE_GPU     
	      #ifdef TEMPORALGAUGE
          to_temporalgauge(g_gauge_field, Q , P);
        #endif          
        iteration_count = linsolve_eo_gpu(P, Q, max_iter, eps_sq, rel_prec, N, f);			     
        #ifdef TEMPORALGAUGE
          from_temporalgauge(Q, P);
        #endif
      #endif
      return(iteration_count);
    }
    else{
      if(f==Qtm_pm_psi){   
        iteration_count =  msolver_fp(P, Q, solver_params, max_iter, eps_sq, rel_prec, N, f, &Qtm_pm_psi_32);
        return(iteration_count);
      }
      else if(f==Q_pm_psi){     
	iteration_count =  msolver_fp(P, Q, solver_params, max_iter, eps_sq, rel_prec, N, f, &Q_pm_psi_32);
	return(iteration_count);      
      } else if(f==Qsw_pm_psi){
        copy_32_sw_fields();
        iteration_count = msolver_fp(P, Q, solver_params, max_iter, eps_sq, rel_prec, N, f, &Qsw_pm_psi_32);
        return(iteration_count);
      } else {
        if(g_proc_id==0) printf("Warning: 32 bit matrix not available. Falling back to CG in 64 bit\n"); 
        use_solver = CG;
      }
    }
  } 
  if(use_solver == CG){
     iteration_count =  cg_her(P, Q, max_iter, eps_sq, rel_prec, N, f);   
  }
  else if(use_solver == BICGSTAB){
     iteration_count =  bicgstab_complex(P, Q, max_iter, eps_sq, rel_prec, N, f);     
  }
#ifdef DDalphaAMG 
  else if (use_solver == MG)
    iteration_count =  MG_solver(P, Q, eps_sq, max_iter,rel_prec, N , g_gauge_field, f);
#endif     
  else{
    if(g_proc_id==0) printf("Error: solver not allowed for degenerate solve. Aborting...\n");
    exit(2);
  }
  return(iteration_count);
}


int solve_mms_nd(spinor ** const Pup, spinor ** const Pdn, 
                 spinor * const Qup, spinor * const Qdn, 
                 solver_pm_t * solver_pm){ 
  int iteration_count = 0; 
    if(solver_pm->type==MIXEDCGMMSND){
      if(usegpu_flag){
	#ifdef HAVE_GPU      
	  #ifdef TEMPORALGAUGE
	    to_temporalgauge_mms(g_gauge_field , Qup, Qdn, Pup, Pdn, solver_pm->no_shifts);
	  #endif        
	  iteration_count = dev_cg_mms_tm_nd(Pup, Pdn, Qup, Qdn, solver_pm);  
	  #ifdef TEMPORALGAUGE
	    from_temporalgauge_mms(Qup, Qdn, Pup, Pdn, solver_pm->no_shifts);
	  #endif 
	#endif
      }
      else{
	iteration_count = mixed_cg_mms_tm_nd(Pup, Pdn, Qup, Qdn, solver_pm);
      }
    }
    else if (solver_pm->type==CGMMSND){
      iteration_count = cg_mms_tm_nd(Pup, Pdn, Qup, Qdn, solver_pm);
    }
    else{
      if(g_proc_id==0) printf("Error: solver not allowed for ND mms solve. Aborting...\n");
      exit(2);      
    }
  return(iteration_count);
}

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
#include "solver/solver.h"
#include "solver/matrix_mult_typedef.h"
#include "monomial_solve.h"

#ifdef HAVE_GPU
#include"../GPU/cudadefs.h"
extern int mixed_solve (spinor * const P, spinor * const Q, const int max_iter, 
                        double eps, const int rel_prec,const int N);
extern  int linsolve_eo_gpu (spinor * const P, spinor * const Q, const int max_iter, 
                            double eps, const int rel_prec, const int N);
extern int dev_cg_mms_tm_nd(spinor ** const Pup, spinor ** const Pdn, 
		 spinor * const Qup, spinor * const Qdn, 
		 solver_pm_t * solver_pm);
   #ifdef TEMPORALGAUGE
     #include "../temporalgauge.h" 
   #endif
#include "read_input.h" 
#endif




int solve_degenerate(spinor * const P, spinor * const Q, const int max_iter, 
           double eps_sq, const int rel_prec, const int N, matrix_mult f){
  int iteration_count;
  if(usegpu_flag){     
    #ifdef TEMPORALGAUGE
      to_temporalgauge(g_gauge_field, Q , P);
    #endif          
    iteration_count = linsolve_eo_gpu(P, Q, max_iter, eps_sq, rel_prec, N);			     
    #ifdef TEMPORALGAUGE
      from_temporalgauge(Q, P);
    #endif     
  }
  else{
    iteration_count =  cg_her(P, Q, max_iter, eps_sq, rel_prec, N, f);
  }
  return(iteration_count);
}


int solve_mms_nd(spinor ** const Pup, spinor ** const Pdn, 
                 spinor * const Qup, spinor * const Qdn, 
                 solver_pm_t * solver_pm){ 
  int iteration_count;  
  if(usegpu_flag){
    #ifdef TEMPORALGAUGE
      to_temporalgauge_mms(g_gauge_field , Qup, Qdn, Pup, Pdn, solver_pm->no_shifts);
    #endif        
    iteration_count = dev_cg_mms_tm_nd(Pup, Pdn, Qup, Qdn, solver_pm);  
    #ifdef TEMPORALGAUGE
      from_temporalgauge_mms(Qup, Qdn, Pup, Pdn, solver_pm->no_shifts);
    #endif 			     
  }
  else{
    iteration_count = cg_mms_tm_nd(Pup, Pdn, Qup, Qdn, solver_pm);
  }
  return(iteration_count);
}

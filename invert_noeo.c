/***********************************************************************
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
 * invert_noeo makes an inversion without EO preconditioned
 * tm Operator
 *
 *
 * invert_noeo returns the number of iterations needed or -1 if the 
 * solver did not converge.
 *
 * Author: Florian Burger
 *         burger@physik.hu-berlin.de
 *
 ****************************************************************/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include<stdlib.h>
#include"global.h"
#include"linalg_eo.h"
#include"operator/tm_operators.h"
#include"operator/D_psi.h"
#include"gamma.h"
#include"solver/solver.h"
#include"read_input.h"
#include"xchange/xchange.h"
#include"invert_noeo.h"


#ifdef HAVE_GPU
#include"GPU/cudadefs.h"
#include"temporalgauge.h"


extern int mixed_solve (spinor * const P, spinor * const Q, const int max_iter, 
			double eps, const int rel_prec,const int N);
#endif

int invert_noeo(spinor * const Spin_new, 
	      spinor * const Spin, 
	      const double precision, const int max_iter,
	      const int solver_flag, const int rel_prec,
              const int id )  {

  int iter = 0;

    if(g_proc_id == 0) {printf("#NOT using even/odd preconditioning!\n"); fflush(stdout);}
    
#ifdef HAVE_GPU
#ifdef TEMPORALGAUGE
    /* initialize temporal gauge here */
    if(usegpu_flag){
      //to_temporalgauge_invert_noeo(g_gauge_field, Spin);   
    } 
#endif  
#endif /* HAVE_GPU*/   
    
    if(solver_flag != CG){
      if(g_proc_id == 0) {printf("Error: Only CG ist implemented as inverter for genuinely non-EO operator! Stopping here.\n"); fflush(stdout);}
      return(-1);
    }
    
    if(g_proc_id == 0) {printf("# Using CG!\n"); fflush(stdout);}
#ifdef HAVE_GPU 
      if(usegpu_flag){
	if(g_proc_id == 0) printf("# Using GPU for inversion\n");
	iter = mixed_solve(Spin_new, Spin, max_iter, precision, rel_prec, VOLUME);
      }
      else{
	gamma5(g_spinor_field[DUM_DERI+1], Spin, VOLUME);
	iter = cg_her(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+1], max_iter, precision, 
		      rel_prec, VOLUME, &Q_pm_psi);
	Q_minus_psi(Spin_new, g_spinor_field[DUM_DERI]);
      }
#else
      gamma5(g_spinor_field[DUM_DERI+1], Spin, VOLUME);
      iter = cg_her(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+1], max_iter, precision, 
		      rel_prec, VOLUME, &Q_pm_psi);
      Q_minus_psi(Spin_new, g_spinor_field[DUM_DERI]);

#endif
    
    
#ifdef HAVE_GPU  
    /* return from temporal gauge again */
#ifdef TEMPORALGAUGE
    if(usegpu_flag){ 
      //from_temporalgauge_invert_noeo(Spin, Spin_new);
    }
#endif
#endif       
  
  
  return(iter);
}


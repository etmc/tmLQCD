/***********************************************************************
 *
 * Copyright (C) 2009 Carsten Urbach
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

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include "global.h"
#include "solver/sumr.h"
#include "solver/cgs_real.h"
#include "operator.h"
#include "invert_overlap.h"
#include "operator/Dov_psi.h"
#include "linalg_eo.h"
#include "read_input.h"
#include "operator/tm_operators.h"
#include "gamma.h"
#include "solver/cg_her.h"


void invert_overlap(const int op_id, const int index_start) {
  operator * optr;
  void (*op)(spinor*,spinor*);
  static _Complex double alpha = 0.;
  spinorPrecWS *ws;
  optr = &operator_list[op_id];
  op=&Dov_psi;

  /* here we need to (re)compute the kernel eigenvectors */
  /* for new gauge fields                                */

  if(g_proc_id == 0) {printf("# Not using even/odd preconditioning!\n"); fflush(stdout);}
  convert_eo_to_lexic(g_spinor_field[DUM_DERI], optr->sr0, optr->sr1);
  convert_eo_to_lexic(g_spinor_field[DUM_DERI+1], optr->prop0, optr->prop1);

  if(optr->solver == 13 ){
    optr->iterations = sumr(g_spinor_field[DUM_DERI+1],g_spinor_field[DUM_DERI] , optr->maxiter, optr->eps_sq);
  } 
  else if(optr->solver == 1 /* CG */) {

    gamma5(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI], VOLUME);
  
    if(use_preconditioning==1 && g_precWS!=NULL){
      ws=(spinorPrecWS*)g_precWS;
      printf("# Using preconditioning (which one?)!\n");
    
      alpha = ws->precExpo[2];
      spinorPrecondition(g_spinor_field[DUM_DERI+1],g_spinor_field[DUM_DERI+1],ws,T,L,alpha,0,1);

      /* 	iter = cg_her(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+1], max_iter, precision,  */
      /* 		    rel_prec, VOLUME, &Q_pm_psi_prec); */
      optr->iterations = cg_her(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+1], optr->maxiter, optr->eps_sq,
				optr->rel_prec, VOLUME, &Qov_sq_psi_prec);
    
      alpha = ws->precExpo[0];
      spinorPrecondition(g_spinor_field[DUM_DERI],g_spinor_field[DUM_DERI],ws,T,L,alpha,0,1);
    
    } 
    else {
      printf("# Not using preconditioning (which one?)!\n");
      /* 	iter = cg_her(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+1], max_iter, precision,  */
      /* 		      rel_prec, VOLUME, &Q_pm_psi); */
      optr->iterations = cg_her(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+1], optr->maxiter, optr->eps_sq,
				optr->rel_prec, VOLUME, &Qov_sq_psi);
    }
  
  
    Qov_psi(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI]);
  
    if(use_preconditioning == 1 && g_precWS!=NULL){
      ws=(spinorPrecWS*)g_precWS;
      alpha = ws->precExpo[1];
      spinorPrecondition(g_spinor_field[DUM_DERI+1],g_spinor_field[DUM_DERI+1],ws,T,L,alpha,0,1);
    }
  
  }
  
  op(g_spinor_field[4],g_spinor_field[DUM_DERI+1]);

  convert_eo_to_lexic(g_spinor_field[DUM_DERI], optr->sr0, optr->sr1);

  optr->reached_prec=diff_and_square_norm(g_spinor_field[4],g_spinor_field[DUM_DERI],VOLUME);
  
  convert_lexic_to_eo(optr->prop0, optr->prop1 , g_spinor_field[DUM_DERI+1]);

  return;
}

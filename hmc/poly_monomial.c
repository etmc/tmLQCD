/***********************************************************************
 * 
 * Copyright (C) 2010 Andreas Nube
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
#include <stdio.h>
#include <math.h>
/*
.
.
.
*/

#include "global.h"
#include "start.h"
#include "read_input.h"
#include "monomial.h"
#include "poly_monomial.h"
#include "boundary.h"
#include "linalg/square_norm.h"
#include "linalg/assign.h"
#include "linalg/mul_r.h"
#include "linalg/diff.h"
#include "linalg_eo.h"
#include "deriv_Sb.h"
#include "linsolve.h"
#include "tm_operators.h"
#include "solver/solver.h"
#include "solver/chrono_guess.h"
#include "solver/eigenvalues.h"
#include "Nondegenerate_Matrix.h"
#include "Hopping_Matrix.h"
#include "phmc.h"



inline void setPhmcVars(monomial *mnl){
  phmc_invmaxev=1.0/mnl->MDPolyLmax;
  phmc_dop_n_cheby=(mnl->MDPolyDegree/2)+1;
  phmc_Cpol=mnl->MDPolyLocNormConst;
  phmc_root=mnl->MDPolyRoots;
}

void poly_derivative(const int id){

  monomial * mnl = &monomial_list[id];
  int k,j;
  int degreehalf=mnl->MDPolyDegree/2;

  spinor** chi_spinor_field=mnl->MDPoly_chi_spinor_fields;


  (*mnl).forcefactor = -2.*mnl->MDPolyLocNormConst/mnl->MDPolyLmax;


  /* push and set phmc vars */
  pushPhmcVars();
  setPhmcVars(mnl);


  if(mnl->even_odd_flag){



    if(mnl->MDPolyDetRatio==1){
      g_mu=mnl->mu2;
      boundary(mnl->kappa2);
      Qtm_plus_psi(chi_spinor_field[0],mnl->pf);
    } else {
      assign(chi_spinor_field[0],mnl->pf,VOLUME/2);
    }


    g_mu=mnl->mu;
    boundary(mnl->kappa);

    
    /* Here comes the definitions for the chi_j fields */
    /* from  j=0  (chi_0 = phi)  .....  to j = n-1 */
    for(k = 0; k < degreehalf-1 ; k++) {
      Qtm_pm_min_cconst_nrm(chi_spinor_field[k+1],
				 chi_spinor_field[k], 
				 mnl->MDPolyRoots[k]);
    }



    assign(chi_spinor_field[degreehalf+1],
	   chi_spinor_field[degreehalf-1], VOLUME/2);
  
    /* loop over monoms */
    for(j=degreehalf; j>=1; j--) {

      assign(chi_spinor_field[degreehalf],
	     chi_spinor_field[degreehalf+1], VOLUME/2);
      
      Qtm_pm_min_cconst_nrm(chi_spinor_field[degreehalf+1], 
			    chi_spinor_field[degreehalf],
			    mnl->MDPolyRoots[mnl->MDPolyDegree-(j+1)]);
      

      Qtm_minus_psi(g_spinor_field[DUM_DERI+3],chi_spinor_field[j-1]); 
      
      H_eo_tm_inv_psi(g_spinor_field[DUM_DERI+2], chi_spinor_field[degreehalf+1], EO, -1.);
      deriv_Sb(OE, g_spinor_field[DUM_DERI+3], g_spinor_field[DUM_DERI+2]); 
      
      H_eo_tm_inv_psi(g_spinor_field[DUM_DERI+2], g_spinor_field[DUM_DERI+3], EO, 1.); 
      deriv_Sb(EO, g_spinor_field[DUM_DERI+2], chi_spinor_field[degreehalf+1]);
      
    
      Qtm_minus_psi(g_spinor_field[DUM_DERI+3],chi_spinor_field[degreehalf+1]); 

      H_eo_tm_inv_psi(g_spinor_field[DUM_DERI+2],g_spinor_field[DUM_DERI+3], EO, +1.);
      deriv_Sb(OE, chi_spinor_field[j-1] , g_spinor_field[DUM_DERI+2]); 
      
      H_eo_tm_inv_psi(g_spinor_field[DUM_DERI+2], chi_spinor_field[j-1], EO, -1.); 
      deriv_Sb(EO, g_spinor_field[DUM_DERI+2], g_spinor_field[DUM_DERI+3]);
      
    }


    if(mnl->MDPolyDetRatio==1){
      /******************************************
      * multiply with the last missing monomial *
      * such that we get an evaluation of P     *
      ******************************************/
      Qtm_pm_min_cconst_nrm(chi_spinor_field[degreehalf], 
			    chi_spinor_field[degreehalf+1],
			    mnl->MDPolyRoots[mnl->MDPolyDegree-1]);
      
    /* devide by this factor cause its multiplied again in update_fermion_momenta see comment below */
    mul_r(chi_spinor_field[degreehalf],
	  1./mnl->MDPolyLocNormConst*mnl->MDPolyLmax,
	  chi_spinor_field[degreehalf],
	  VOLUME/2);


      g_mu=mnl->mu2;
      boundary(mnl->kappa2);

      H_eo_tm_inv_psi(g_spinor_field[DUM_DERI+2],chi_spinor_field[degreehalf], EO, -1.);
      deriv_Sb(OE, mnl->pf , g_spinor_field[DUM_DERI+2]);
      
      H_eo_tm_inv_psi(g_spinor_field[DUM_DERI+2], mnl->pf, EO, +1.);
      deriv_Sb(EO, g_spinor_field[DUM_DERI+2], chi_spinor_field[degreehalf]);



    }




  } else {
    if(g_proc_id == 0) {
      fprintf(stderr,"Error: PHMC for light quarks not implementeted for non even/odd preconditioning\n");
    }

    g_mu = g_mu1;
    boundary(g_kappa);
    popPhmcVars();

    return;
  }

  /* restore all changed global vars */
  g_mu = g_mu1;
  boundary(g_kappa);
  popPhmcVars();


}

double poly_acc(const int id){

  monomial * mnl = &monomial_list[id];
  int j;
  spinor* spinor1=g_spinor_field[2];
  spinor* spinor2=g_spinor_field[3];
  double diff;
  int no_eigenvalues=-1;




  if(mnl->even_odd_flag){

    if(mnl->MDPolyDetRatio==1){
      g_mu = mnl->mu2;
      boundary(mnl->kappa2);

      Qtm_plus_psi(spinor2,mnl->pf);

    } else {
      assign(spinor2,mnl->pf,VOLUME/2);
    }

    g_mu = mnl->mu;
    boundary(mnl->kappa);

    /* push and set phmc vars */
    pushPhmcVars();
    setPhmcVars(mnl);

    /* apply B */
    for(j = 0; j < mnl->MDPolyDegree/2; j++){
      assign(spinor1, spinor2, VOLUME/2);
      Qtm_pm_min_cconst_nrm(spinor2,
			    spinor1,
			    mnl->MDPolyRoots[j]);
    }

    mnl->energy1 =  square_norm(spinor2, VOLUME/2,1);

    /* calculate evs */
    if (compute_evs != 0) {
      no_eigenvalues=10;
      eigenvalues(&no_eigenvalues, mnl->maxiter, eigenvalue_precision,
                  0/* compute minimal evs*/, 0/*dont write evecs*/, nstore, mnl->even_odd_flag);

      no_eigenvalues=1;
      eigenvalues(&no_eigenvalues, mnl->maxiter, eigenvalue_precision,
                  1/* compute maximal evs*/, 0/*dont write evecs*/, nstore, mnl->even_odd_flag);
    }


    /* restore global phmc vars */
    popPhmcVars();

    
    /* return the energy differnce */
    g_mu = g_mu1;
    boundary(g_kappa);


  
    if(g_proc_id == 0 && g_debug_level > 3) {
      fprintf(stderr," Poly energy1     = %e \n" , mnl->energy1);
      fprintf(stderr," Poly energy0     = %e \n" , mnl->energy0);
      diff = mnl->energy1 - mnl->energy0;
      fprintf(stderr," Poly energy diff = %e \n" , diff);
    }
    
    return (mnl->energy1 - mnl->energy0);
  } else {
    if(g_proc_id == 0) {
      fprintf(stderr,"Error: PHMC for light quarks not implementeted for non even/odd preconditioning\n");
    }

    g_mu = g_mu1;
    boundary(g_kappa);

    return NAN;
  }


}

void poly_heatbath(const int id){
  monomial * mnl = &monomial_list[id];
  int j;
  spinor* spinor1=g_spinor_field[2];
  spinor* spinor2=g_spinor_field[3];

  mnl->csg_n = 0;
  mnl->csg_n2 = 0;
  mnl->iter0 = 0;
  mnl->iter1 = 0;

  g_mu = mnl->mu;
  boundary(mnl->kappa);
  
  /* push and set phmc vars */
  pushPhmcVars();
  setPhmcVars(mnl);

  if(mnl->even_odd_flag) {


    random_spinor_field(spinor1, VOLUME/2, mnl->rngrepro);
    mnl->energy0 = square_norm(spinor1, VOLUME/2, 1);

    if(g_proc_id == 0 && g_debug_level > 3) {
      fprintf(stderr," Poly energy0     = %e \n" , mnl->energy0);
    }

    /* calculate the phmc hamiltonian */
    Qtm_pm_psi(spinor2, spinor1);

    /* solve (Q+)*(Q-)*P((Q+)*(Q-)) *x=y */
    cg_her(spinor1, spinor2,
	   1000,mnl->accprec,g_relative_precision_flag,VOLUME/2, Qtm_pm_Ptm_pm_psi);
    
    /*  phi= Bdagger phi  */
    for(j = 0; j < (mnl->MDPolyDegree/2); j++){
      assign(spinor2, spinor1, VOLUME/2);
      Qtm_pm_min_cconst_nrm(spinor1,
				 spinor2,
				 mnl->MDPolyRoots[mnl->MDPolyDegree/2+j]);
    }


    if(mnl->MDPolyDetRatio==1){
      g_mu = mnl->mu2;
      boundary(mnl->kappa2);
      zero_spinor_field(mnl->pf,VOLUME/2);
      if(mnl->solver == CG) ITER_MAX_BCG = 0;
      ITER_MAX_CG = mnl->maxiter;
      mnl->iter0 += bicg(mnl->pf, spinor1, mnl->accprec, g_relative_precision_flag);
      
      chrono_add_solution(mnl->pf, mnl->csg_field, mnl->csg_index_array,
			  mnl->csg_N, &mnl->csg_n, VOLUME/2);
      
      if(mnl->solver != CG) {
	chrono_add_solution(mnl->pf, mnl->csg_field2, mnl->csg_index_array2,
			    mnl->csg_N2, &mnl->csg_n2, VOLUME/2);
      }
    } else {
      /* store constructed phi field */
      assign(mnl->pf, spinor1, VOLUME/2);
    }
    
  }
  else {
    /* not implemented */
    fprintf(stderr,"Error: non even odd preconditioned \"light\" phmc not implemented \n");
  }

  g_mu = g_mu1;
  boundary(g_kappa);
  popPhmcVars();

  if(g_proc_id == 0 && g_debug_level > 3) {
    printf("called poly_heatbath for id %d %d\n", id, mnl->even_odd_flag);
  }


  return;
}

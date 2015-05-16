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
 *
 * invert_doublet_eo makes an inversion with EO precoditioned
 * tm Operator with a nondegenerate doublet
 *
 * Even and Odd are the numbers of spinor_field that contain
 * the even and the odd sites of the source. The result is stored
 * int Even_new and Odd_new.
 *
 * invert_doublet_eo returns the number of iterations neede or -1 if the 
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
#include"gamma.h"
#include"solver/solver.h"
#include"read_input.h"
#include"xchange/xchange.h"
#include"operator/tm_operators_nd.h"
#include"invert_doublet_eo.h"


#ifdef HAVE_GPU
#  include"GPU/cudadefs.h"
#  include"temporalgauge.h"
#  include"measure_gauge_action.h"
int mixedsolve_eo_nd (spinor *, spinor *, spinor *, spinor *, double, int, double, int, int, matrix_mult_nd);
#  ifdef TEMPORALGAUGE
extern su3* g_trafo;
#  endif
#endif


int invert_doublet_eo(spinor * const Even_new_s, spinor * const Odd_new_s, 
		      spinor * const Even_new_c, spinor * const Odd_new_c, 
		      spinor * const Even_s, spinor * const Odd_s,
		      spinor * const Even_c, spinor * const Odd_c,
		      const double precision, const int max_iter,
		      const int solver_flag, const int rel_prec, const int even_odd_flag) {

  int iter = 0;

  if(even_odd_flag){

#ifdef HAVE_GPU
  /* initialize temporal gauge here */
  int retval;
  double dret1, dret2;
  double plaquette1 = 0.0;
  double plaquette2 = 0.0;
  
  if (usegpu_flag) {
    #ifdef TEMPORALGAUGE 
      to_temporalgauge_invert_doublet_eo(g_gauge_field, Even_s, Odd_s, Even_c, Odd_c);
    #endif
  } 
#endif  



  /* here comes the inversion using even/odd preconditioning */
  if(g_proc_id == 0) {printf("# Using even/odd preconditioning!\n"); fflush(stdout);}
  M_ee_inv_ndpsi(Even_new_s, Even_new_c, 
		 Even_s, Even_c,
		 g_mubar, g_epsbar);
  Hopping_Matrix(OE, g_spinor_field[DUM_DERI], Even_new_s);
  Hopping_Matrix(OE, g_spinor_field[DUM_DERI+1], Even_new_c);
  
  /* The sign is plus, since in Hopping_Matrix */
  /* the minus is missing                      */
  assign_mul_add_r(g_spinor_field[DUM_DERI], +1., Odd_s, VOLUME/2);
  assign_mul_add_r(g_spinor_field[DUM_DERI+1], +1., Odd_c, VOLUME/2);
  
  /* Do the inversion with the preconditioned  */
  /* matrix to get the odd sites               */
  
  /* Here we invert the hermitean operator squared */
  
  if(g_proc_id == 0) {
    printf("# Using CG for TMWILSON flavour doublet!\n"); 
    fflush(stdout);
  }
  gamma5(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI], VOLUME/2);
  gamma5(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI+1], VOLUME/2);
  
  
  if (usegpu_flag) {	// GPU, mixed precision solver, shift==0
   #ifdef HAVE_GPU
    iter = mixedsolve_eo_nd(Odd_new_s, Odd_new_c, g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+1],
			    0.0, max_iter, precision, rel_prec, even_odd_flag, &Qtm_pm_ndpsi);
   #endif
  }
  else {
    iter = cg_her_nd(Odd_new_s, Odd_new_c, g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+1],
		    max_iter, precision, rel_prec, 
		    VOLUME/2, &Qtm_pm_ndpsi);
  }
  
  
  Qtm_dagger_ndpsi(Odd_new_s, Odd_new_c,
		   Odd_new_s, Odd_new_c);

  /* Reconstruct the even sites                */
  Hopping_Matrix(EO, g_spinor_field[DUM_DERI], Odd_new_s);
  Hopping_Matrix(EO, g_spinor_field[DUM_DERI+1], Odd_new_c);
  M_ee_inv_ndpsi(g_spinor_field[DUM_DERI+2], g_spinor_field[DUM_DERI+3],
		 g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+1],
		 g_mubar, g_epsbar);
  
  /* The sign is plus, since in Hopping_Matrix */
  /* the minus is missing                      */
  assign_add_mul_r(Even_new_s, g_spinor_field[DUM_DERI+2], +1., VOLUME/2);
  assign_add_mul_r(Even_new_c, g_spinor_field[DUM_DERI+3], +1., VOLUME/2);
  
  
  }
  else{
    
    /* here comes the inversion not using even/odd preconditioning */
    if(g_proc_id == 0) {printf("# Not using even/odd preconditioning!\n"); fflush(stdout);}
    
    //if odd set to NULL (e.g. if global T is odd) then do not 
    //use EO management of fields but use the Even and Even_new pointers
    //exclusively, supposed to hold the complete field
    if((Odd_s == NULL) || (Odd_c == NULL)){
      assign(g_spinor_field[DUM_DERI],   Even_s, VOLUME);
      assign(g_spinor_field[DUM_DERI+1], Even_c, VOLUME);
      assign(g_spinor_field[DUM_DERI+2], Even_new_s, VOLUME);
      assign(g_spinor_field[DUM_DERI+3], Even_new_c, VOLUME);
    }
    else{
      convert_eo_to_lexic(g_spinor_field[DUM_DERI], Even_s, Odd_s);
      convert_eo_to_lexic(g_spinor_field[DUM_DERI+1], Even_c, Odd_c);
      convert_eo_to_lexic(g_spinor_field[DUM_DERI+2], Even_new_s, Odd_new_s);
      convert_eo_to_lexic(g_spinor_field[DUM_DERI+3], Even_new_c, Odd_new_c); 
    }
   
    
    gamma5(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI], VOLUME);
    gamma5(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI+1], VOLUME);   
    
    if (usegpu_flag) {	// GPU, mixed precision solver, shift==0
      #ifdef HAVE_GPU
	iter = mixedsolve_eo_nd(g_spinor_field[DUM_DERI+2], g_spinor_field[DUM_DERI+3], g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+1],
				0.0, max_iter, precision, rel_prec, even_odd_flag, &Q_pm_ndpsi);
      #endif
    }
    else{
      iter = cg_her_nd(g_spinor_field[DUM_DERI+2], g_spinor_field[DUM_DERI+3], g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+1],
		     max_iter, precision, rel_prec, 
		     VOLUME, &Q_pm_ndpsi);
    }
    Q_minus_ndpsi(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI+2], g_spinor_field[DUM_DERI+3]);
    
    if((Odd_s == NULL) || (Odd_c == NULL)){
      assign(Even_new_s, g_spinor_field[DUM_DERI], VOLUME);
      assign(Even_new_c, g_spinor_field[DUM_DERI+1], VOLUME);
    }
    else{
      convert_lexic_to_eo(Even_new_s, Odd_new_s, g_spinor_field[DUM_DERI]);
      convert_lexic_to_eo(Even_new_c, Odd_new_c, g_spinor_field[DUM_DERI+1]);
    }
  }
  
  if (usegpu_flag) {
    #ifdef HAVE_GPU
      /* return from temporal gauge again */
      #ifdef TEMPORALGAUGE
	from_temporalgauge_invert_doublet_eo(Even_s, Odd_s, Even_new_s, Odd_new_s,
					  Even_c, Odd_c, Even_new_c, Odd_new_c);
      #endif  
    #endif  
  }  

  return(iter);
}


int invert_cloverdoublet_eo(spinor * const Even_new_s, spinor * const Odd_new_s, 
			    spinor * const Even_new_c, spinor * const Odd_new_c, 
			    spinor * const Even_s, spinor * const Odd_s,
			    spinor * const Even_c, spinor * const Odd_c,
			    const double precision, const int max_iter,
			    const int solver_flag, const int rel_prec) {
  
  int iter = 0;
  
  
  /* here comes the inversion using even/odd preconditioning */
  if(g_proc_id == 0) {printf("# Using even/odd preconditioning!\n"); fflush(stdout);}
  Msw_ee_inv_ndpsi(Even_new_s, Even_new_c, 
		   Even_s, Even_c);
  Hopping_Matrix(OE, g_spinor_field[DUM_DERI], Even_new_s);
  Hopping_Matrix(OE, g_spinor_field[DUM_DERI+1], Even_new_c);
  
  /* The sign is plus, since in Hopping_Matrix */
  /* the minus is missing                      */
  assign_mul_add_r(g_spinor_field[DUM_DERI], +1., Odd_s, VOLUME/2);
  assign_mul_add_r(g_spinor_field[DUM_DERI+1], +1., Odd_c, VOLUME/2);
  
  /* Do the inversion with the preconditioned  */
  /* matrix to get the odd sites               */
  
  /* Here we invert the hermitean operator squared */
  
  if(g_proc_id == 0) {
    printf("# Using CG for TMWILSON flavour doublet!\n"); 
    fflush(stdout);
  }
  gamma5(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI], VOLUME/2);
  gamma5(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI+1], VOLUME/2);
  
  iter = cg_her_nd(Odd_new_s, Odd_new_c, g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+1],
		   max_iter, precision, rel_prec, 
		   VOLUME/2, &Qsw_pm_ndpsi);
  
  
  Qsw_dagger_ndpsi(Odd_new_s, Odd_new_c,
		   Odd_new_s, Odd_new_c);
  
  /* Reconstruct the even sites                */
  Hopping_Matrix(EO, g_spinor_field[DUM_DERI], Odd_new_s);
  Hopping_Matrix(EO, g_spinor_field[DUM_DERI+1], Odd_new_c);
  Msw_ee_inv_ndpsi(g_spinor_field[DUM_DERI+2], g_spinor_field[DUM_DERI+3],
		   g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+1]);
  
  /* The sign is plus, since in Hopping_Matrix */
  /* the minus is missing                      */
  assign_add_mul_r(Even_new_s, g_spinor_field[DUM_DERI+2], +1., VOLUME/2);
  assign_add_mul_r(Even_new_c, g_spinor_field[DUM_DERI+3], +1., VOLUME/2);
  
  return(iter);
}


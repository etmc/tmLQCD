/***********************************************************************
 *
 * Copyright (C) 2006,2007,2008 Karl Jansen, Thomas Chiarappa, 
 *                              Carsten Urbach
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
 * This file contains operators for twisted mass Wilson QCD   
 * to construct a multiplication with a non-degenerate        
 * flavour matrix                                             
 *                                                            
 *                                                            
 ***********************************************************************/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "global.h"
#include "su3.h"
#include "operator/Hopping_Matrix.h"
#include "phmc.h"
#include "gamma.h"
#include "linalg_eo.h"
#include "operator/D_psi.h"
#include "operator/tm_operators.h"
#include "operator/clovertm_operators.h"
#include "operator/tm_operators_nd.h"


void mul_one_pm_iconst(spinor * const l, spinor * const k, 
		       const double mu_, const int sign_);

void M_oo_sub_g5_ndpsi(spinor * const l_s, spinor * const l_c, 
		       spinor * const k_s, spinor * const k_c,
		       spinor * const j_s, spinor * const j_c,
		       const double mu, const double eps);

/* external functions */


/******************************************
 *
 * This is the implementation of
 *
 *  M_full_ndpsi = D_w I_f + i gamma5 mubar tau3 - epsbar tau1
 *  the full operator done for testing purpose
 ******************************************/
void M_full_ndpsi(spinor * const Even_new_s, spinor * const Odd_new_s, 
                  spinor * const Even_new_c, spinor * const Odd_new_c, 
                  spinor * const Even_s, spinor * const Odd_s,
                  spinor * const Even_c, spinor * const Odd_c) {
  
  double mu = g_mu;
  g_mu = g_mubar;
  M_full(Even_new_s, Odd_new_s, Even_s, Odd_s);

  assign_add_mul_r(Even_new_s, Even_c, -g_epsbar, VOLUME/2);
  assign_add_mul_r(Odd_new_s, Odd_c, -g_epsbar, VOLUME/2);
  
  g_mu = -g_mu;
  M_full(Even_new_c, Odd_new_c, Even_c, Odd_c);
  
  assign_add_mul_r(Even_new_c, Even_s, -g_epsbar, VOLUME/2);
  assign_add_mul_r(Odd_new_c, Odd_s, -g_epsbar, VOLUME/2);

  g_mu = mu;
}

void Msw_full_ndpsi(spinor * const Even_new_s, spinor * const Odd_new_s, 
                    spinor * const Even_new_c, spinor * const Odd_new_c, 
                    spinor * const Even_s, spinor * const Odd_s,
                    spinor * const Even_c, spinor * const Odd_c) {

  double mu = g_mu;
  g_mu = g_mubar;
  Msw_full(Even_new_s, Odd_new_s, Even_s, Odd_s);

  assign_add_mul_r(Even_new_s, Even_c, -g_epsbar, VOLUME/2);
  assign_add_mul_r(Odd_new_s, Odd_c, -g_epsbar, VOLUME/2);
  
  g_mu = -g_mu;
  Msw_full(Even_new_c, Odd_new_c, Even_c, Odd_c);
  
  assign_add_mul_r(Even_new_c, Even_s, -g_epsbar, VOLUME/2);
  assign_add_mul_r(Odd_new_c, Odd_s, -g_epsbar, VOLUME/2);

  g_mu = mu;
}

// full VOLUME operator; it used D_psi which works with tm and tm+clover
void D_ndpsi(spinor * const l_strange, spinor * const l_charm,
             spinor * const k_strange, spinor * const k_charm) {

  double mu = g_mu;
  g_mu = g_mubar;
  D_psi(l_strange,k_strange);

  assign_add_mul_r(l_strange, k_charm, -g_epsbar, VOLUME);
  
  g_mu = -g_mu;
  D_psi(l_charm,k_charm);
  
  assign_add_mul_r(l_charm, k_strange, -g_epsbar, VOLUME);

  g_mu = mu;
}

/******************************************
 *
 * This is the implementation of
 *
 * Qhat(2x2) = gamma_5 * [ M_oo - M_oe M_ee^-1 M_eo ]
 *
 * see documentation for details
 * k_charm and k_strange are the input fields
 * l_* the output fields
 *
 * it acts only on the odd part or only 
 * on a half spinor
 ******************************************/
void Qtm_ndpsi(spinor * const l_strange, spinor * const l_charm,
	       spinor * const k_strange, spinor * const k_charm){

  /* Here the  M_oe Mee^-1 M_eo  implementation  */
  Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX], k_strange);
  Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX+1], k_charm);

  M_ee_inv_ndpsi(g_spinor_field[DUM_MATRIX+3], g_spinor_field[DUM_MATRIX+2],
		 g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+1],
		 g_mubar, g_epsbar);
  
  Hopping_Matrix(OE, l_strange, g_spinor_field[DUM_MATRIX+3]);
  Hopping_Matrix(OE, l_charm, g_spinor_field[DUM_MATRIX+2]);

  /* Here the M_oo  implementation  */
  M_oo_sub_g5_ndpsi(g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+1], k_strange, k_charm,
  		    l_strange, l_charm,
  		    -g_mubar, -g_epsbar);
  /* At the end, the normalisation by the max. eigenvalue  */
  mul_r(l_strange, phmc_invmaxev, g_spinor_field[DUM_MATRIX], VOLUME/2);
  mul_r(l_charm, phmc_invmaxev, g_spinor_field[DUM_MATRIX+1], VOLUME/2);
}

void Qsw_ndpsi(spinor * const l_strange, spinor * const l_charm,
		spinor * const k_strange, spinor * const k_charm) {

  Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX], k_charm);
  Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX+1], k_strange);

  assign_mul_one_sw_pm_imu_eps(EE, g_spinor_field[DUM_MATRIX+2], g_spinor_field[DUM_MATRIX+3], 
			       g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+1], g_mubar, g_epsbar);
  clover_inv_nd(EE, g_spinor_field[DUM_MATRIX+2], g_spinor_field[DUM_MATRIX+3]);

  Hopping_Matrix(OE, g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+2]);
  Hopping_Matrix(OE, g_spinor_field[DUM_MATRIX+1], g_spinor_field[DUM_MATRIX+3]);

  clover_gamma5_nd(OO, g_spinor_field[DUM_MATRIX+2], g_spinor_field[DUM_MATRIX+3], 
  		   k_charm, k_strange,
  		   g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+1],
  		   g_mubar, -g_epsbar);
  mul_r(l_charm, phmc_invmaxev, g_spinor_field[DUM_MATRIX+2], VOLUME/2);
  mul_r(l_strange, phmc_invmaxev, g_spinor_field[DUM_MATRIX+3], VOLUME/2);
  return;
}

/******************************************
 *
 * This is the implementation of 
 *
 *  Q_tau1_ndpsi_add/sub_Ishift =  ( M +/- I z_k )
 *
 *  with M = Qhat(2x2) tau_1   and z_k from sqrt(g_shift) 
 *
 *
 *  needed in the evaluation of the heatbath when 
 *  the Rational approximation is used
 *
 *
 * For details, see documentation and comments of the
 * above mentioned routines
 *
 * k_charm and k_strange are the input fields
 * l_* the output fields
 *
 * it acts only on the odd part or only
 * on a half spinor
 ******************************************/


void Qtm_tau1_ndpsi_add_Ishift(spinor * const l_strange, spinor * const l_charm,
                               spinor * const k_strange, spinor * const k_charm) {

  Q_tau1_sub_const_ndpsi(l_strange,l_charm,k_strange,k_charm,-I*sqrt(g_shift),1.,phmc_invmaxev);

  return;
}

void Qtm_tau1_ndpsi_sub_Ishift(spinor * const l_strange, spinor * const l_charm,
                               spinor * const k_strange, spinor * const k_charm) {

  Q_tau1_sub_const_ndpsi(l_strange,l_charm,k_strange,k_charm, I*sqrt(g_shift),1.,phmc_invmaxev);

  return;
}

void Qsw_tau1_ndpsi_add_Ishift(spinor * const l_strange, spinor * const l_charm,
                               spinor * const k_strange, spinor * const k_charm) {

  Qsw_tau1_sub_const_ndpsi(l_strange,l_charm,k_strange,k_charm,-I*sqrt(g_shift),1.,phmc_invmaxev);

  return;
}

void Qsw_tau1_ndpsi_sub_Ishift(spinor * const l_strange, spinor * const l_charm,
                               spinor * const k_strange, spinor * const k_charm) {

  Qsw_tau1_sub_const_ndpsi(l_strange,l_charm,k_strange,k_charm, I*sqrt(g_shift),1.,phmc_invmaxev);

  return;
}


/******************************************
 *
 * This is the implementation of
 * 
 * Qhat(2x2)^dagger = tau_1  Qhat(2x2) tau_1 =
 *
 *  = Qhat(2x2)  with   g_mubar  ->  - g_mubar
 *
 * With respect to Qtm_ndpsi the role of charme and strange fields
 * are interchenged, since Qdagger=tau_1 Q tau_1
 * see documentation for details
 * k_charm and k_strange are the input fields
 * l_* the output fields
 *
 * it acts only on the odd part or only
 * on a half spinor
 ******************************************/
void Qtm_dagger_ndpsi(spinor * const l_strange, spinor * const l_charm,
		      spinor * const k_strange, spinor * const k_charm) {

  /* Here the  M_oe Mee^-1 M_eo  implementation  */
  Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX], k_charm);
  Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX+1], k_strange);

  M_ee_inv_ndpsi(g_spinor_field[DUM_MATRIX+2], g_spinor_field[DUM_MATRIX+3],
		 g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+1],
		 g_mubar, g_epsbar);
  
  Hopping_Matrix(OE, g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+2]);
  Hopping_Matrix(OE, g_spinor_field[DUM_MATRIX+1], g_spinor_field[DUM_MATRIX+3]);

  /* Here the M_oo  implementation  */
  M_oo_sub_g5_ndpsi(l_strange, l_charm, k_strange, k_charm,
  		    g_spinor_field[DUM_MATRIX+1], g_spinor_field[DUM_MATRIX],
  		    g_mubar, -g_epsbar);
  /* At the end, the normalisation by the max. eigenvalue  */
  mul_r(l_charm, phmc_invmaxev, l_charm, VOLUME/2);
  mul_r(l_strange, phmc_invmaxev, l_strange, VOLUME/2);

}

void Qsw_dagger_ndpsi(spinor * const l_strange, spinor * const l_charm,
		      spinor * const k_strange, spinor * const k_charm) {

  Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX], k_charm);
  Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX+1], k_strange);

  assign_mul_one_sw_pm_imu_eps(EE, g_spinor_field[DUM_MATRIX+2], g_spinor_field[DUM_MATRIX+3], 
			       g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+1], -g_mubar, g_epsbar);
  clover_inv_nd(EE, g_spinor_field[DUM_MATRIX+2], g_spinor_field[DUM_MATRIX+3]);

  Hopping_Matrix(OE, g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+2]);
  Hopping_Matrix(OE, g_spinor_field[DUM_MATRIX+1], g_spinor_field[DUM_MATRIX+3]);

  clover_gamma5_nd(OO, g_spinor_field[DUM_MATRIX+2], g_spinor_field[DUM_MATRIX+3], 
  		   k_charm, k_strange,
  		   g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+1],
  		   -g_mubar, -g_epsbar);
  mul_r(l_charm, phmc_invmaxev, g_spinor_field[DUM_MATRIX+2], VOLUME/2);
  mul_r(l_strange, phmc_invmaxev, g_spinor_field[DUM_MATRIX+3], VOLUME/2);
  return;
}


/******************************************
 *
 * This is the implementation of
 * 
 * Qhat(2x2) Qhat(2x2)^dagger 
 *
 *
 * For details, see documentation and comments of the 
 * above mentioned routines 
 *
 * k_charm and k_strange are the input fields
 * l_* the output fields
 *
 * l_ and k_ can be identical
 *
 * it acts only on the odd part or only
 * on a half spinor
 ******************************************/
void Qtm_pm_ndpsi(spinor * const l_strange, spinor * const l_charm,
		  spinor * const k_strange, spinor * const k_charm){

  /* first the  Qhat(2x2)^dagger  PART*/
  /* Here the  M_oe Mee^-1 M_eo  implementation  */
  Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX], k_charm);
  Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX+1], k_strange);

  M_ee_inv_ndpsi(g_spinor_field[DUM_MATRIX+2], g_spinor_field[DUM_MATRIX+3],
		 g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+1],
		 g_mubar, g_epsbar);

  Hopping_Matrix(OE, g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+2]);
  Hopping_Matrix(OE, g_spinor_field[DUM_MATRIX+1], g_spinor_field[DUM_MATRIX+3]);

  /* Here the M_oo  implementation  */
  M_oo_sub_g5_ndpsi(g_spinor_field[DUM_MATRIX+2], g_spinor_field[DUM_MATRIX+3], k_charm, k_strange,
  		    g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+1],
  		    -g_mubar, -g_epsbar);
  /* We have to reassigin as follows to avoid overwriting */
  /* Recall in fact that   Q^hat = tau_1 Q tau_1  , hence  */
  /* and then the  Qhat(2x2)  PART */

  /* Here the  M_oe Mee^-1 M_eo  implementation  */
  Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+3]);
  Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX+1], g_spinor_field[DUM_MATRIX+2]);

  M_ee_inv_ndpsi(g_spinor_field[DUM_MATRIX+5], g_spinor_field[DUM_MATRIX+4],
		 g_spinor_field[DUM_MATRIX+1], g_spinor_field[DUM_MATRIX],
		 -g_mubar, g_epsbar);

  Hopping_Matrix(OE, l_strange, g_spinor_field[DUM_MATRIX+4]);
  Hopping_Matrix(OE, l_charm, g_spinor_field[DUM_MATRIX+5]);

  /* Here the M_oo  implementation  */
  M_oo_sub_g5_ndpsi(l_strange, l_charm, g_spinor_field[DUM_MATRIX+3], g_spinor_field[DUM_MATRIX+2],
		    l_strange, l_charm,
  		    -g_mubar, -g_epsbar);
  /* At the end, the normalisation by the max. eigenvalue  */ 
  /* Twice  phmc_invmaxev  since we consider here  D Ddag  !!! */
  mul_r(l_charm, phmc_invmaxev*phmc_invmaxev, l_charm, VOLUME/2);
  mul_r(l_strange, phmc_invmaxev*phmc_invmaxev, l_strange, VOLUME/2);
  return;
}

void Qtm_pm_ndpsi_shift(spinor * const l_strange, spinor * const l_charm,
                       spinor * const k_strange, spinor * const k_charm) {
  Qtm_pm_ndpsi(l_strange,l_charm,k_strange,k_charm);  
  assign_add_mul_r( l_strange, k_strange, g_shift, VOLUME/2 );
  assign_add_mul_r( l_charm, k_charm, g_shift, VOLUME/2 );
  return;
}

void Qsw_pm_ndpsi(spinor * const l_strange, spinor * const l_charm,
		  spinor * const k_strange, spinor * const k_charm) {

  /* FIRST THE  Qhat(2x2)^dagger  PART*/
  /* Here the  M_oe Mee^-1 M_eo  implementation  */
  Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX], k_charm);
  Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX+1], k_strange);

  assign_mul_one_sw_pm_imu_eps(EE, g_spinor_field[DUM_MATRIX+2], g_spinor_field[DUM_MATRIX+3], 
			       g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+1], -g_mubar, g_epsbar);
  clover_inv_nd(EE, g_spinor_field[DUM_MATRIX+2], g_spinor_field[DUM_MATRIX+3]);

  Hopping_Matrix(OE, g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+2]);
  Hopping_Matrix(OE, g_spinor_field[DUM_MATRIX+1], g_spinor_field[DUM_MATRIX+3]);

  // Here the M_oo  implementation  
  clover_gamma5_nd(OO, g_spinor_field[DUM_MATRIX+2], g_spinor_field[DUM_MATRIX+3], 
  		   k_charm, k_strange,
  		   g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+1],
  		   -g_mubar, -g_epsbar);

  // and then the  Qhat(2x2)  PART 
  // Recall in fact that   Q^hat = tau_1 Q tau_1  
  // Here the  M_oe Mee^-1 M_eo  implementation  
  // the re-ordering in s and c components is due to tau_1
  Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+3]);
  Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX+1], g_spinor_field[DUM_MATRIX+2]);

  assign_mul_one_sw_pm_imu_eps(EE, g_spinor_field[DUM_MATRIX+7], g_spinor_field[DUM_MATRIX+6], 
			       g_spinor_field[DUM_MATRIX+1], g_spinor_field[DUM_MATRIX], g_mubar, g_epsbar);
  clover_inv_nd(EE, g_spinor_field[DUM_MATRIX+6], g_spinor_field[DUM_MATRIX+7]);

  Hopping_Matrix(OE, g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+6]);
  Hopping_Matrix(OE, g_spinor_field[DUM_MATRIX+1], g_spinor_field[DUM_MATRIX+7]);

  clover_gamma5_nd(OO, l_charm, l_strange,
  		   g_spinor_field[DUM_MATRIX+2], g_spinor_field[DUM_MATRIX+3],
  		   g_spinor_field[DUM_MATRIX+1], g_spinor_field[DUM_MATRIX],
  		   g_mubar, -g_epsbar);

  /* At the end, the normalisation by the max. eigenvalue  */ 
  /* Twice  phmc_invmaxev  since we consider here  D Ddag  !!! */
  mul_r(l_charm, phmc_invmaxev*phmc_invmaxev, l_charm, VOLUME/2);
  mul_r(l_strange, phmc_invmaxev*phmc_invmaxev, l_strange, VOLUME/2);
  return;
}

void Qsw_pm_ndpsi_shift(spinor * const l_strange, spinor * const l_charm,
                       spinor * const k_strange, spinor * const k_charm) {
  Qsw_pm_ndpsi(l_strange,l_charm,k_strange,k_charm);
  
  assign_add_mul_r( l_strange, k_strange, g_shift, VOLUME/2 );
  assign_add_mul_r( l_charm, k_charm, g_shift, VOLUME/2 );

  return;
}


/******************************************
 *
 * This is the implementation of 
 *
 *  Q_tau1_sub_const_ndpsi =  Cpol*( M - z_k )
 *
 *  with M = Qhat(2x2) tau_1   and z_k \in Complex
 *
 *
 *  needed in the evaluation of the forces when 
 *  the Polynomial approximation is used
 *
 *
 * For details, see documentation and comments of the
 * above mentioned routines
 *
 * k_charm and k_strange are the input fields
 * l_* the output fields
 *
 * it acts only on the odd part or only
 * on a half spinor
 ******************************************/
void Q_tau1_sub_const_ndpsi(spinor * const l_strange, spinor * const l_charm,
			    spinor * const k_strange, spinor * const k_charm, 
			    const _Complex double z, const double Cpol, const double invev) {

  spinor *r, *s;
  su3_vector ALIGN phi1;

  /*   tau_1   inverts the   k_charm  <->  k_strange   spinors */
  /*  Apply first  Qhat(2x2)  and finally substract the constant  */

  /* Here the  M_oe Mee^-1 M_eo  implementation  */

  Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX], k_charm);
  Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX+1], k_strange);

  M_ee_inv_ndpsi(g_spinor_field[DUM_MATRIX+3], g_spinor_field[DUM_MATRIX+2],
				g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+1],
				g_mubar, g_epsbar);

  Hopping_Matrix(OE, l_strange, g_spinor_field[DUM_MATRIX+3]);
  Hopping_Matrix(OE, l_charm, g_spinor_field[DUM_MATRIX+2]);

  /* Here the M_oo  implementation  */
  M_oo_sub_g5_ndpsi(g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+1], k_charm, k_strange,
  		    l_strange, l_charm,
  		    -g_mubar, -g_epsbar);

  /* At the end, the normalisation by the max. eigenvalue  */
  mul_r(l_strange, Cpol*invev, g_spinor_field[DUM_MATRIX], VOLUME/2);
  mul_r(l_charm, Cpol*invev, g_spinor_field[DUM_MATRIX+1], VOLUME/2);

  /* Finally, we add k to l and multiply all */
  /* by the constant  phmc_Cpol  */
  /* which renders the polynomial in monomials  */
  /* identical to the polynomial a la clenshaw */;
#ifdef TM_USE_OMP
#pragma omp parallel for private(r) private(s) private(phi1)
#endif
  for(int ix = 0; ix < (VOLUME/2); ix++){

    r=l_strange + ix;
    s=k_strange + ix;
    
    _complex_times_vector(phi1, Cpol*z, s->s0);
    _vector_sub_assign(r->s0, phi1);
    _complex_times_vector(phi1, Cpol*z, s->s1);
    _vector_sub_assign(r->s1, phi1);
    _complex_times_vector(phi1, Cpol*z, s->s2);
    _vector_sub_assign(r->s2, phi1);
    _complex_times_vector(phi1, Cpol*z, s->s3);
    _vector_sub_assign(r->s3, phi1);

    r=l_charm + ix;
    s=k_charm + ix;
    
    _complex_times_vector(phi1, Cpol*z, s->s0);
    _vector_sub_assign(r->s0, phi1);
    _complex_times_vector(phi1, Cpol*z, s->s1);
    _vector_sub_assign(r->s1, phi1);
    _complex_times_vector(phi1, Cpol*z, s->s2);
    _vector_sub_assign(r->s2, phi1);
    _complex_times_vector(phi1, Cpol*z, s->s3);
    _vector_sub_assign(r->s3, phi1);    
  }
  return;
}

void Qsw_tau1_sub_const_ndpsi(spinor * const l_strange, spinor * const l_charm,
			      spinor * const k_strange, spinor * const k_charm, 
			      const _Complex double z, const double Cpol, const double invev) {

  spinor *r, *s;
  su3_vector ALIGN phi1;

  /*   tau_1   inverts the   k_charm  <->  k_strange   spinors */
  /*  Apply first  Qhat(2x2)  and finally substract the constant  */

  /* Here the  M_oe Mee^-1 M_eo  implementation  */

  Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX], k_charm);
  Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX+1], k_strange);

  assign_mul_one_sw_pm_imu_eps(EE, g_spinor_field[DUM_MATRIX+3], g_spinor_field[DUM_MATRIX+2], 
			       g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+1], -g_mubar, g_epsbar);
  clover_inv_nd(EE, g_spinor_field[DUM_MATRIX+2], g_spinor_field[DUM_MATRIX+3]);

  Hopping_Matrix(OE, l_strange, g_spinor_field[DUM_MATRIX+3]);
  Hopping_Matrix(OE, l_charm, g_spinor_field[DUM_MATRIX+2]);

  /* Here the M_oo  implementation  */
  clover_gamma5_nd(OO, g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+1], 
  		   k_charm, k_strange,
  		   l_strange, l_charm,
  		   -g_mubar, -g_epsbar);

  /* At the end, the normalisation by the max. eigenvalue  */
  mul_r(l_strange, Cpol*invev, g_spinor_field[DUM_MATRIX], VOLUME/2);
  mul_r(l_charm, Cpol*invev, g_spinor_field[DUM_MATRIX+1], VOLUME/2);

  /* Finally, we add k to l and multiply all */
  /* by the constant  phmc_Cpol  */
  /* which renders the polynomial in monomials  */
  /* identical to the polynomial a la clenshaw */;
#ifdef TM_USE_OMP
#pragma omp parallel for private(r) private(s) private(phi1)
#endif
  for(int ix = 0; ix < (VOLUME/2); ix++){

    r=l_strange + ix;
    s=k_strange + ix;
    
    _complex_times_vector(phi1, Cpol*z, s->s0);
    _vector_sub_assign(r->s0, phi1);
    _complex_times_vector(phi1, Cpol*z, s->s1);
    _vector_sub_assign(r->s1, phi1);
    _complex_times_vector(phi1, Cpol*z, s->s2);
    _vector_sub_assign(r->s2, phi1);
    _complex_times_vector(phi1, Cpol*z, s->s3);
    _vector_sub_assign(r->s3, phi1);

    r=l_charm + ix;
    s=k_charm + ix;
    
    _complex_times_vector(phi1, Cpol*z, s->s0);
    _vector_sub_assign(r->s0, phi1);
    _complex_times_vector(phi1, Cpol*z, s->s1);
    _vector_sub_assign(r->s1, phi1);
    _complex_times_vector(phi1, Cpol*z, s->s2);
    _vector_sub_assign(r->s2, phi1);
    _complex_times_vector(phi1, Cpol*z, s->s3);
    _vector_sub_assign(r->s3, phi1);    
  }
  return;
}




/******************************************
 *
 * This is the same implementation as above of
 * 
 * Qhat(2x2) Qhat(2x2)^dagger 
 *
 *
 *  but now input and output are bispinors !!!!
 *
 * For details, see documentation and comments of the 
 * above mentioned routines 
 *
 * k_charm and k_strange are the input fields
 * l_* the output fields
 *
 * it acts only on the odd part or only
 * on a half spinor
 ******************************************/
void Qtm_pm_ndbipsi(bispinor * const bisp_l, bispinor * const bisp_k) {

  /*  create 2 spinors out of 1 (input) bispinor  */
  decompact(g_spinor_field[DUM_MATRIX+6], g_spinor_field[DUM_MATRIX+7], bisp_k);

  Qtm_pm_ndpsi(g_spinor_field[DUM_MATRIX+6], g_spinor_field[DUM_MATRIX+7],
	       g_spinor_field[DUM_MATRIX+6], g_spinor_field[DUM_MATRIX+7]);

  /*  create 1 (output) bispinor out of 2 spinors  */
  compact(bisp_l, g_spinor_field[DUM_MATRIX+6], g_spinor_field[DUM_MATRIX+7]);
  return;
}

void Qsw_pm_ndbipsi(bispinor * const bisp_l, bispinor * const bisp_k) {

  /*  create 2 spinors out of 1 (input) bispinor  */
  decompact(g_spinor_field[DUM_MATRIX+6], g_spinor_field[DUM_MATRIX+7], bisp_k);

  Qsw_pm_ndpsi(g_spinor_field[DUM_MATRIX+6], g_spinor_field[DUM_MATRIX+7],
	       g_spinor_field[DUM_MATRIX+6], g_spinor_field[DUM_MATRIX+7]);

  /*  create 1 (output) bispinor out of 2 spinors  */
  compact(bisp_l, g_spinor_field[DUM_MATRIX+6], g_spinor_field[DUM_MATRIX+7]);
  return;
}


/******************************************
 *
 * This is the implementation of
 *
 * (M_{ee}^\pm)^{-1}M_{eo} tau^1
 *
 * see documentation for details
 * k is the number of the input field
 * l is the number of the output field
 *
 * it acts only on the odd part or only 
 * on a half spinor
 ******************************************/

void H_eo_tm_ndpsi(spinor * const l_strange, spinor * const l_charm, 
             spinor * const k_strange, spinor * const k_charm, 
	     const int ieo) {
  /* recall:   strange <-> up    while    charm <-> dn   */
  Hopping_Matrix(ieo, g_spinor_field[DUM_MATRIX], k_strange);
  Hopping_Matrix(ieo, g_spinor_field[DUM_MATRIX+1], k_charm);

  M_ee_inv_ndpsi(l_charm, l_strange,
		 g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+1],
		 -g_mubar, g_epsbar);
  return;
}

void H_eo_sw_ndpsi(spinor * const l_strange, spinor * const l_charm, 
		   spinor * const k_strange, spinor * const k_charm) {

  /* recall:   strange <-> up    while    charm <-> dn   */
  Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX], k_strange);
  Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX+1], k_charm);
  
  assign_mul_one_sw_pm_imu_eps(EE, l_charm, l_strange,
			       g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+1], 
			       g_mubar, g_epsbar);
  // here the order doesn't matter
  clover_inv_nd(EE, l_strange, l_charm);

  return;
}

// for this routine we need to have sw_invert_nd and sw_term called before hand
// and the clover term must be initialised
void Msw_ee_inv_ndpsi(spinor * const l_strange, spinor * const l_charm, 
		      spinor * const k_strange, spinor * const k_charm) {
  

  /* recall:   strange <-> up    while    charm <-> dn   */

  assign_mul_one_sw_pm_imu_eps(EE, l_strange, l_charm, k_strange, k_charm, -g_mubar, g_epsbar);

  clover_inv_nd(EE, l_strange, l_charm);
  return;
}



void Q_test_epsilon(spinor * const l_strange, spinor * const l_charm,
                           spinor * const k_strange, spinor * const k_charm){

  double nrm = 1./(1.+g_mubar*g_mubar-g_epsbar*g_epsbar);

  /* Here the  M_oe Mee^-1 M_eo  implementation  */
  Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX], k_strange);
  Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX+1], k_charm);

  Hopping_Matrix(OE, g_spinor_field[DUM_MATRIX+2], g_spinor_field[DUM_MATRIX]);
  Hopping_Matrix(OE, g_spinor_field[DUM_MATRIX+3], g_spinor_field[DUM_MATRIX+1]);

  assign_add_mul_r(k_strange, g_spinor_field[DUM_MATRIX+2], nrm, VOLUME/2);
  assign_add_mul_r(k_charm, g_spinor_field[DUM_MATRIX+3], nrm, VOLUME/2);

  mul_r(l_strange, -2, k_strange, VOLUME/2);
  mul_r(l_charm, -2, k_charm, VOLUME/2);

  /* and finally the  gamma_5  multiplication  */
  gamma5(l_strange, l_strange, VOLUME/2);
  gamma5(l_charm, l_charm, VOLUME/2);

  /* At the end, the normalisation by the max. eigenvalue  */
  mul_r(l_charm, phmc_invmaxev, l_charm, VOLUME/2);
  mul_r(l_strange, phmc_invmaxev, l_strange, VOLUME/2);
  return;
}


void mul_one_pm_itau2(spinor * const p, spinor * const q,
		      spinor * const r, spinor * const s,
		      const double sign, const int N) {
  double fac = 1./sqrt(2.);

  if(sign > 0) {
    add(p, r, s, N);
    diff(q, s, r, N);
  }
  else {
    diff(p, r, s, N);
    add(q, r, s, N);
  }
  mul_r(p, fac, p, N);
  mul_r(q, fac, q, N);
}

void mul_one_pm_iconst(spinor * const l, spinor * const k, 
		       const double mu_, const int sign_) {
#ifdef TM_USE_OMP
#pragma omp parallel
  {
#endif

  spinor *r, *s;
  su3_vector ALIGN phi1;
  double mu = mu_;
  if(sign_ < 0) {
    mu = -mu_;
  }

  /************ loop over all lattice sites ************/
#ifdef TM_USE_OMP
#pragma omp for
#endif
  for(unsigned int ix = 0; ix < (VOLUME/2); ++ix){
    r=l + ix;
    s=k + ix;
    /* Multiply the spinorfield with 1+imu\gamma_5 */
    _complex_times_vector(phi1, (1. + mu * I), s->s0);
    _vector_assign(r->s0, phi1);
    _complex_times_vector(phi1, (1. + mu * I), s->s1);
    _vector_assign(r->s1, phi1);
    _complex_times_vector(phi1, (1. - mu * I), s->s2);
    _vector_assign(r->s2, phi1);
    _complex_times_vector(phi1, (1. - mu * I), s->s3);
    _vector_assign(r->s3, phi1);
  }

#ifdef TM_USE_OMP
  } /* OpenMP closing brace */
#endif

  return;
}

// l_ and k_ are allowed to be the same spinors
void M_ee_inv_ndpsi(spinor * const l_s, spinor * const l_c, 
		    spinor * const k_s, spinor * const k_c,
		    const double mu, const double eps) {
#ifdef TM_USE_OMP
#pragma omp parallel
  {
#endif
  double nrm = 1./(1.+ mu*mu - eps*eps);
  spinor *r_s, *r_c, *s_s, *s_c;
  su3_vector ALIGN phi1, phi2;

#ifdef TM_USE_OMP
#pragma omp for
#endif
  for(unsigned int ix = 0; ix < (VOLUME/2); ++ix){
    r_s = l_s + ix;
    r_c = l_c + ix;
    s_s = k_s + ix;
    s_c = k_c + ix;

    _complex_times_vector(phi1, (1. - mu * I), s_s->s0);
    _vector_add_mul(phi1, eps, s_c->s0);
    _complex_times_vector(phi2, (1. + mu * I), s_c->s0);
    _vector_add_mul(phi2, eps, s_s->s0);
    _vector_mul(r_s->s0, nrm, phi1);
    _vector_mul(r_c->s0, nrm, phi2);

    _complex_times_vector(phi1, (1. - mu * I), s_s->s1);
    _vector_add_mul(phi1, eps, s_c->s1);
    _complex_times_vector(phi2, (1. + mu * I), s_c->s1);
    _vector_add_mul(phi2, eps, s_s->s1);
    _vector_mul(r_s->s1, nrm, phi1);
    _vector_mul(r_c->s1, nrm, phi2);

    _complex_times_vector(phi1, (1. + mu * I), s_s->s2);
    _vector_add_mul(phi1, eps, s_c->s2);
    _complex_times_vector(phi2, (1. - mu * I), s_c->s2);
    _vector_add_mul(phi2, eps, s_s->s2);
    _vector_mul(r_s->s2, nrm, phi1);
    _vector_mul(r_c->s2, nrm, phi2);

    _complex_times_vector(phi1, (1. + mu * I), s_s->s3);
    _vector_add_mul(phi1, eps, s_c->s3);
    _complex_times_vector(phi2, (1. - mu * I), s_c->s3);
    _vector_add_mul(phi2, eps, s_s->s3);
    _vector_mul(r_s->s3, nrm, phi1);
    _vector_mul(r_c->s3, nrm, phi2);

  }

#ifdef TM_USE_OMP
  } /* OpenMP closing brace */
#endif

  return;
}


// l_ and k_ are allowed to be the same spinors
void M_oo_sub_g5_ndpsi(spinor * const l_s, spinor * const l_c, 
		       spinor * const k_s, spinor * const k_c,
		       spinor * const j_s, spinor * const j_c,
		       const double mu, const double eps) {
#ifdef TM_USE_OMP
#pragma omp parallel
  {
#endif
  spinor *r_s, *r_c, *s_s, *s_c, *t_s, *t_c;
  su3_vector ALIGN phi1, phi2;

#ifdef TM_USE_OMP
#pragma omp for
#endif
  for(unsigned int ix = 0; ix < (VOLUME/2); ++ix){
    r_s = l_s + ix;
    r_c = l_c + ix;
    s_s = k_s + ix;
    s_c = k_c + ix;
    t_s = j_s + ix;
    t_c = j_c + ix;

    _complex_times_vector(phi1, (1. - mu * I), s_s->s0);
    _vector_add_mul(phi1, eps, s_c->s0);
    _complex_times_vector(phi2, (1. + mu * I), s_c->s0);
    _vector_add_mul(phi2, eps, s_s->s0);
    _vector_sub(r_s->s0, phi1, t_s->s0);
    _vector_sub(r_c->s0, phi2, t_c->s0);

    _complex_times_vector(phi1, (1. - mu * I), s_s->s1);
    _vector_add_mul(phi1, eps, s_c->s1);
    _complex_times_vector(phi2, (1. + mu * I), s_c->s1);
    _vector_add_mul(phi2, eps, s_s->s1);
    _vector_sub(r_s->s1, phi1, t_s->s1);
    _vector_sub(r_c->s1, phi2, t_c->s1);

    _complex_times_vector(phi1, (1. + mu * I), s_s->s2);
    _vector_add_mul(phi1, eps, s_c->s2);
    _complex_times_vector(phi2, (1. - mu * I), s_c->s2);
    _vector_add_mul(phi2, eps, s_s->s2);
    _vector_sub(r_s->s2, t_s->s2, phi1);
    _vector_sub(r_c->s2, t_c->s2, phi2);

    _complex_times_vector(phi1, (1. + mu * I), s_s->s3);
    _vector_add_mul(phi1, eps, s_c->s3);
    _complex_times_vector(phi2, (1. - mu * I), s_c->s3);
    _vector_add_mul(phi2, eps, s_s->s3);
    _vector_sub(r_s->s3, t_s->s3, phi1);
    _vector_sub(r_c->s3, t_c->s3, phi2);
  }

#ifdef TM_USE_OMP
  } /* OpenMP closing brace */
#endif

  return;
}


/*  calculates P(Q Q^dagger) for the nondegenerate case */

void P_ndpsi(spinor * const l_strange, spinor * const l_charm,
	     spinor * const k_strange, spinor * const k_charm){
  
  
  
  int j;
  spinor *dum_up,*dum_dn;
  dum_up=g_chi_up_spinor_field[DUM_MATRIX];
  dum_dn=g_chi_dn_spinor_field[DUM_MATRIX];
  
  assign(dum_up,k_strange,VOLUME/2);
  assign(dum_dn,k_charm,VOLUME/2);
    
  for(j = 0; j < (2*phmc_dop_n_cheby -2); j++) {
    if(j>0) {
      assign(dum_up,l_strange,VOLUME/2);
      assign(dum_dn,l_charm,VOLUME/2);
    }
    
    Q_tau1_sub_const_ndpsi(l_strange, l_charm,
			   dum_up, dum_dn,
			   phmc_root[j], phmc_Cpol, phmc_invmaxev);
  }
  return;
}


/* calculates  Q * \tau^1  for the nondegenerate case */
void Qtau1_P_ndpsi(spinor * const l_strange, spinor * const l_charm,
		   spinor * const k_strange, spinor * const k_charm){
  
  
  spinor * dum_up,* dum_dn;
  dum_up = g_chi_up_spinor_field[DUM_MATRIX+1];
  dum_dn = g_chi_dn_spinor_field[DUM_MATRIX+1];
  
  P_ndpsi(l_strange, l_charm, k_strange, k_charm);
  
  assign(dum_up, l_strange, VOLUME/2);
  assign(dum_dn, l_charm, VOLUME/2);
  
  Qtm_ndpsi(l_strange, l_charm, dum_dn, dum_up);
  return;
}



/* this is neccessary for the calculation of the polynomial */

void Qtm_pm_sub_const_nrm_psi(spinor * const l, spinor * const k,
			      const _Complex double z){
  su3_vector ALIGN phi1;
  spinor *r,*s;
  int ix;

  Qtm_pm_psi(l, k);
  mul_r(l, phmc_invmaxev, l, VOLUME/2);

  /*  AND FINALLY WE SUBSTRACT THE C-CONSTANT  */


  /************ loop over all lattice sites ************/
#ifdef TM_USE_OMP
#pragma omp parallel for private(ix) private(r) private(s) private(phi1)
#endif
  for(ix = 0; ix < (VOLUME/2); ix++){

    r=l + ix;
    s=k + ix;

    _complex_times_vector(phi1, z, s->s0);
    _vector_sub_assign(r->s0, phi1);
    _complex_times_vector(phi1, z, s->s1);
    _vector_sub_assign(r->s1, phi1);
    _complex_times_vector(phi1, z, s->s2);
    _vector_sub_assign(r->s2, phi1);
    _complex_times_vector(phi1, z, s->s3);
    _vector_sub_assign(r->s3, phi1);
  }

  mul_r(l, phmc_Cpol, l, VOLUME/2);
  return;
}

/* calculate a polynomial in (Q+)*(Q-) */


void Ptm_pm_psi(spinor * const l, spinor * const k){
  
  int j;
  spinor *spinDum;
  spinDum=g_spinor_field[DUM_MATRIX+2];
  
  assign(spinDum,k,VOLUME/2);
  
  
  for(j=0; j<(2*phmc_dop_n_cheby -2); j++){
    if(j>0) {
      assign(spinDum,l,VOLUME/2);
    }
    
    Qtm_pm_sub_const_nrm_psi(l,spinDum,phmc_root[j]);
  }
  return;
}

/* **********************************************
 * Qpm * P(Qpm)
 * this operator is neccessary for the inverter 
 ************************************************/

void Qtm_pm_Ptm_pm_psi(spinor * const l, spinor * const k){
  spinor * spinDum;
  
  spinDum=g_spinor_field[DUM_MATRIX+3]; 
  Ptm_pm_psi(l,k);
  assign(spinDum,l,VOLUME/2);
  Qtm_pm_psi(l,spinDum);
  return;
}


/* ************************************************
 * for noise reduction 
 * this implements
 * a = B^dagger H b
 * 
 * with Hopping matrix H and
 *
 * B = (1-i\g5\tau^1\musigma-\tau^3\mudelta)/c
 * where
 * c = 1+\musigma^2-\mudelta^2
 *
 * so it is in the convention of hep-lat/0606011
 * not in the internal one, see documentation
 * 
 **************************************************/

void red_noise_nd(spinor * const lse, spinor * const lso, 
		  spinor * const lce, spinor * const lco) 
{

  double nrm0 = (1.-g_epsbar)/(1+g_mubar*g_mubar-g_epsbar*g_epsbar);
  double nrm1 = (1.+g_epsbar)/(1+g_mubar*g_mubar-g_epsbar*g_epsbar);
  _Complex double z;
  int ix, i;
  su3_vector ALIGN phi;
  spinor * r, * s;
  
  /* need B^\dagger, so change sign of g_mubar */
  z = (g_mubar / (1 + g_mubar * g_mubar - g_epsbar * g_epsbar)) * I;

  /* first multiply with Hopping matrix */
  Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX], lso);
  Hopping_Matrix(OE, g_spinor_field[DUM_MATRIX+1], lse);

  Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX+2], lco);
  Hopping_Matrix(OE, g_spinor_field[DUM_MATRIX+3], lce);
  
  /* now with A^{-1}*/
  mul_r(lse, nrm0, g_spinor_field[DUM_MATRIX], VOLUME/2);
  mul_r(lso, nrm0, g_spinor_field[DUM_MATRIX+1], VOLUME/2);

  mul_r(lce, nrm1, g_spinor_field[DUM_MATRIX+2], VOLUME/2);
  mul_r(lco, nrm1, g_spinor_field[DUM_MATRIX+3], VOLUME/2);

  /************ loop over all lattice sites ************/
  for(i = 0; i < 4; i++) {
    if(i == 0) {
      r = lse, s = g_spinor_field[DUM_MATRIX];
    }
    else if(i == 1) {
      r = lso, s = g_spinor_field[DUM_MATRIX+1];
    }
    else if(i == 2) {
      r = lce, s = g_spinor_field[DUM_MATRIX+2];
    }
    else {
      r = lco, s = g_spinor_field[DUM_MATRIX+3];
    }
    for(ix = 0; ix < (VOLUME/2); ix++){
      /* Multiply the spinorfield with (i epsbar \gamma_5)/c */
      /* and add it to */
      _complex_times_vector(phi, z, s->s0);
      _vector_add_assign(r->s0, phi);
      _complex_times_vector(phi, z, s->s1);
      _vector_add_assign(r->s1, phi);
      _complex_times_vector(phi, -z, s->s2);
      _vector_add_assign(r->s2, phi);
      _complex_times_vector(phi, -z, s->s3);
      _vector_add_assign(r->s3, phi);
      r++; s++;
    }  
  }
  return;
}







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
#include "Hopping_Matrix.h"
#include "phmc.h"
#include "gamma.h"
#include "linsolve.h"
#include "linalg_eo.h"
#include "Nondegenerate_Matrix.h"


void mul_one_minus_imubar(spinor * const l, spinor * const k);
/******************************************
 * mul_one_plus_imubar_inv computes
 * l = [(1-i\mubar\gamma_5) * l
 *
*/

void mul_one_plus_imubar(spinor * const l, spinor * const k);
/******************************************
 * mul_one_plus_imubar_inv computes
 * l = [(1+i\mubar\gamma_5) * l
 *
*/

void Qtm_pm_psi(spinor *l,spinor *k);

/* external functions */

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
void QNon_degenerate(spinor * const l_strange, spinor * const l_charm,
                     spinor * const k_strange, spinor * const k_charm){

  double nrm = 1./(1.+g_mubar*g_mubar-g_epsbar*g_epsbar);

  /* Here the  M_oe Mee^-1 M_eo  implementation  */
  Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX], k_strange);
  Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX+1], k_charm);

  mul_one_minus_imubar(g_spinor_field[DUM_MATRIX+4], g_spinor_field[DUM_MATRIX]);
  mul_one_plus_imubar(g_spinor_field[DUM_MATRIX+3], g_spinor_field[DUM_MATRIX+1]);

  assign_add_mul_r(g_spinor_field[DUM_MATRIX+4], g_spinor_field[DUM_MATRIX+1], g_epsbar, VOLUME/2);
  assign_add_mul_r(g_spinor_field[DUM_MATRIX+3], g_spinor_field[DUM_MATRIX], g_epsbar, VOLUME/2);

  mul_r(g_spinor_field[DUM_MATRIX+4], nrm, g_spinor_field[DUM_MATRIX+4], VOLUME/2);
  mul_r(g_spinor_field[DUM_MATRIX+3], nrm, g_spinor_field[DUM_MATRIX+3], VOLUME/2);
  /* where nrm (= 1/(1+mu^2 -eps^2)) has been defined at the beginning of 
     the subroutine */
  
  Hopping_Matrix(OE, l_strange, g_spinor_field[DUM_MATRIX+4]);
  Hopping_Matrix(OE, l_charm, g_spinor_field[DUM_MATRIX+3]);

  /* Here the M_oo  implementation  */
  mul_one_plus_imubar(g_spinor_field[DUM_MATRIX], k_strange);
  mul_one_minus_imubar(g_spinor_field[DUM_MATRIX+1], k_charm);

  assign_add_mul_r(g_spinor_field[DUM_MATRIX], k_charm, -g_epsbar, VOLUME/2);
  assign_add_mul_r(g_spinor_field[DUM_MATRIX+1], k_strange, -g_epsbar, VOLUME/2);
   
  diff(l_strange, g_spinor_field[DUM_MATRIX], l_strange, VOLUME/2);
  diff(l_charm, g_spinor_field[DUM_MATRIX+1], l_charm, VOLUME/2);

  /* and finally the  gamma_5  multiplication  */
  gamma5(l_strange, l_strange, VOLUME/2);
  gamma5(l_charm, l_charm, VOLUME/2);

  /* At the end, the normalisation by the max. eigenvalue  */
  mul_r(l_strange, phmc_invmaxev, l_strange, VOLUME/2);
  mul_r(l_charm, phmc_invmaxev, l_charm, VOLUME/2);
}

/******************************************
 *
 * This is the implementation of
 * 
 * Qhat(2x2)^dagger = tau_1  Qhat(2x2) tau_1 =
 *
 *  = Qhat(2x2)  with   g_mubar  ->  - g_mubar
 *
 * With respect to QNon_degenerate the role of charme and strange fields
 * are interchenged, since Qdagger=tau_1 Q tau_1
 * see documentation for details
 * k_charm and k_strange are the input fields
 * l_* the output fields
 *
 * it acts only on the odd part or only
 * on a half spinor
 ******************************************/
void QdaggerNon_degenerate(spinor * const l_strange, spinor * const l_charm,
                           spinor * const k_strange, spinor * const k_charm){

  double nrm = 1./(1.+g_mubar*g_mubar-g_epsbar*g_epsbar);

  /* Here the  M_oe Mee^-1 M_eo  implementation  */
  Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX], k_charm);
  Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX+1], k_strange);

  mul_one_minus_imubar(g_spinor_field[DUM_MATRIX+2], g_spinor_field[DUM_MATRIX]);
  mul_one_plus_imubar(g_spinor_field[DUM_MATRIX+3], g_spinor_field[DUM_MATRIX+1]);

  assign_add_mul_r(g_spinor_field[DUM_MATRIX+2], g_spinor_field[DUM_MATRIX+1], g_epsbar, VOLUME/2);
  assign_add_mul_r(g_spinor_field[DUM_MATRIX+3], g_spinor_field[DUM_MATRIX], g_epsbar, VOLUME/2);

  mul_r(g_spinor_field[DUM_MATRIX+2], nrm, g_spinor_field[DUM_MATRIX+2], VOLUME/2);
  mul_r(g_spinor_field[DUM_MATRIX+3], nrm, g_spinor_field[DUM_MATRIX+3], VOLUME/2);
  /* where nrm (= 1/(1+mu^2 -eps^2)) has been defined at the beginning of 
     the subroutine */
  
  Hopping_Matrix(OE, g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+2]);
  Hopping_Matrix(OE, g_spinor_field[DUM_MATRIX+1], g_spinor_field[DUM_MATRIX+3]);

  /* Here the M_oo  implementation  */
  mul_one_plus_imubar(g_spinor_field[DUM_MATRIX+2], k_charm);
  mul_one_minus_imubar(g_spinor_field[DUM_MATRIX+3], k_strange);

  assign_add_mul_r(g_spinor_field[DUM_MATRIX+2], k_strange, -g_epsbar, VOLUME/2);
  assign_add_mul_r(g_spinor_field[DUM_MATRIX+3], k_charm, -g_epsbar, VOLUME/2);
   
  diff(l_charm, g_spinor_field[DUM_MATRIX+2], g_spinor_field[DUM_MATRIX], VOLUME/2);
  diff(l_strange, g_spinor_field[DUM_MATRIX+3], g_spinor_field[DUM_MATRIX+1], VOLUME/2);

  /* and finally the  gamma_5  multiplication  */
  gamma5(l_charm, l_charm, VOLUME/2);
  gamma5(l_strange, l_strange, VOLUME/2);

  /* At the end, the normalisation by the max. eigenvalue  */
  mul_r(l_charm, phmc_invmaxev, l_charm, VOLUME/2);
  mul_r(l_strange, phmc_invmaxev, l_strange, VOLUME/2);

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
 * it acts only on the odd part or only
 * on a half spinor
 ******************************************/
void Q_Qdagger_ND(spinor * const l_strange, spinor * const l_charm,
                           spinor * const k_strange, spinor * const k_charm){

  double nrm = 1./(1.+g_mubar*g_mubar-g_epsbar*g_epsbar);

  /* FIRST THE  Qhat(2x2)^dagger  PART*/
  /* Here the  M_oe Mee^-1 M_eo  implementation  */
  Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX], k_charm);
  Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX+1], k_strange);

  mul_one_minus_imubar(g_spinor_field[DUM_MATRIX+2], g_spinor_field[DUM_MATRIX]);
  mul_one_plus_imubar(g_spinor_field[DUM_MATRIX+3], g_spinor_field[DUM_MATRIX+1]);

  assign_add_mul_r(g_spinor_field[DUM_MATRIX+2], g_spinor_field[DUM_MATRIX+1], g_epsbar, VOLUME/2);
  assign_add_mul_r(g_spinor_field[DUM_MATRIX+3], g_spinor_field[DUM_MATRIX], g_epsbar, VOLUME/2);

  mul_r(g_spinor_field[DUM_MATRIX+2], nrm, g_spinor_field[DUM_MATRIX+2], VOLUME/2);
  mul_r(g_spinor_field[DUM_MATRIX+3], nrm, g_spinor_field[DUM_MATRIX+3], VOLUME/2);
  
  Hopping_Matrix(OE, g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+2]);
  Hopping_Matrix(OE, g_spinor_field[DUM_MATRIX+1], g_spinor_field[DUM_MATRIX+3]);

  /* Here the M_oo  implementation  */
  mul_one_plus_imubar(g_spinor_field[DUM_MATRIX+2], k_charm);
  mul_one_minus_imubar(g_spinor_field[DUM_MATRIX+3], k_strange);

  assign_add_mul_r(g_spinor_field[DUM_MATRIX+2], k_strange, -g_epsbar, VOLUME/2);
  assign_add_mul_r(g_spinor_field[DUM_MATRIX+3], k_charm, -g_epsbar, VOLUME/2);
   
  diff(g_spinor_field[DUM_MATRIX+4], g_spinor_field[DUM_MATRIX+2], g_spinor_field[DUM_MATRIX], VOLUME/2);
  diff(g_spinor_field[DUM_MATRIX+5], g_spinor_field[DUM_MATRIX+3], g_spinor_field[DUM_MATRIX+1], VOLUME/2);

  /* and finally the  gamma_5  multiplication  */
  gamma5(g_spinor_field[DUM_MATRIX+2], g_spinor_field[DUM_MATRIX+4], VOLUME/2);
  gamma5(g_spinor_field[DUM_MATRIX+3], g_spinor_field[DUM_MATRIX+5], VOLUME/2);

  /* The normalisation by the max. eigenvalue  is done twice at the end */


  /* We have to reassigin as follows to avoid overwriting */
  /* Recall in fact that   Q^hat = tau_1 Q tau_1  , hence  */

  /*  ABOVE: dum_matrix+2  is  l_charm   goes to  dum_matrix+6 :BELOW */
  /*  ABOVE: dum_matrix+3  is  l_strange   goes to  dum_matrix+7 :BELOW */
  assign(g_spinor_field[DUM_MATRIX+6], g_spinor_field[DUM_MATRIX+2], VOLUME/2);
  assign(g_spinor_field[DUM_MATRIX+7], g_spinor_field[DUM_MATRIX+3], VOLUME/2);


  /* AND THEN THE  Qhat(2x2)  PART */

  /* Here the  M_oe Mee^-1 M_eo  implementation  */
  Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+7]);
  Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX+1], g_spinor_field[DUM_MATRIX+6]);

  mul_one_minus_imubar(g_spinor_field[DUM_MATRIX+2], g_spinor_field[DUM_MATRIX]);
  mul_one_plus_imubar(g_spinor_field[DUM_MATRIX+3], g_spinor_field[DUM_MATRIX+1]);

  assign_add_mul_r(g_spinor_field[DUM_MATRIX+2], g_spinor_field[DUM_MATRIX+1], g_epsbar, VOLUME/2);
  assign_add_mul_r(g_spinor_field[DUM_MATRIX+3], g_spinor_field[DUM_MATRIX], g_epsbar, VOLUME/2);

  mul_r(g_spinor_field[DUM_MATRIX+2], nrm, g_spinor_field[DUM_MATRIX+2], VOLUME/2);
  mul_r(g_spinor_field[DUM_MATRIX+3], nrm, g_spinor_field[DUM_MATRIX+3], VOLUME/2);
 
  Hopping_Matrix(OE, l_strange, g_spinor_field[DUM_MATRIX+2]);
  Hopping_Matrix(OE, l_charm, g_spinor_field[DUM_MATRIX+3]);

  /* Here the M_oo  implementation  */
  mul_one_plus_imubar(g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+7]);
  mul_one_minus_imubar(g_spinor_field[DUM_MATRIX+1], g_spinor_field[DUM_MATRIX+6]);

  assign_add_mul_r(g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+6], -g_epsbar, VOLUME/2);
  assign_add_mul_r(g_spinor_field[DUM_MATRIX+1], g_spinor_field[DUM_MATRIX+7], -g_epsbar, VOLUME/2);
   
  diff(l_strange, g_spinor_field[DUM_MATRIX], l_strange, VOLUME/2);
  diff(l_charm, g_spinor_field[DUM_MATRIX+1], l_charm, VOLUME/2);

  /* and finally the  gamma_5  multiplication  */
  gamma5(l_strange, l_strange, VOLUME/2);
  gamma5(l_charm, l_charm, VOLUME/2);


  /* At the end, the normalisation by the max. eigenvalue  */ 
  /* Twice  phmc_invmaxev  since we consider here  D Ddag  !!! */
  mul_r(l_charm, phmc_invmaxev*phmc_invmaxev, l_charm, VOLUME/2);
  mul_r(l_strange, phmc_invmaxev*phmc_invmaxev, l_strange, VOLUME/2);
  return;
}


/******************************************
 *
 * This is the implementation of 
 *
 *  Q_tau1_min_cconst_ND =  M - z_k 
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
void Q_tau1_min_cconst_ND(spinor * const l_strange, spinor * const l_charm,
                     spinor * const k_strange, spinor * const k_charm, const _Complex double z){


  int ix;
  spinor *r, *s;
  static su3_vector phi1;

  double nrm = 1./(1.+g_mubar*g_mubar-g_epsbar*g_epsbar);


  /*   tau_1   inverts the   k_charm  <->  k_strange   spinors */
  /*  Apply first  Qhat(2x2)  and finally substract the constant  */

  /* Here the  M_oe Mee^-1 M_eo  implementation  */

  Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX], k_charm);
  Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX+1], k_strange);


  mul_one_minus_imubar(g_spinor_field[DUM_MATRIX+4], g_spinor_field[DUM_MATRIX]);
  mul_one_plus_imubar(g_spinor_field[DUM_MATRIX+3], g_spinor_field[DUM_MATRIX+1]);


  assign_add_mul_r(g_spinor_field[DUM_MATRIX+4], g_spinor_field[DUM_MATRIX+1], g_epsbar, VOLUME/2);
  assign_add_mul_r(g_spinor_field[DUM_MATRIX+3], g_spinor_field[DUM_MATRIX], g_epsbar, VOLUME/2);

  mul_r(g_spinor_field[DUM_MATRIX+4], nrm, g_spinor_field[DUM_MATRIX+4], VOLUME/2);
  mul_r(g_spinor_field[DUM_MATRIX+3], nrm, g_spinor_field[DUM_MATRIX+3], VOLUME/2);
  /* where nrm (= 1/(1+mu^2 -eps^2)) has been defined at the beginning of
     the subroutine */


  Hopping_Matrix(OE, l_strange, g_spinor_field[DUM_MATRIX+4]);
  Hopping_Matrix(OE, l_charm, g_spinor_field[DUM_MATRIX+3]);


  /* Here the M_oo  implementation  */
  mul_one_plus_imubar(g_spinor_field[DUM_MATRIX], k_charm);
  mul_one_minus_imubar(g_spinor_field[DUM_MATRIX+1], k_strange);


  assign_add_mul_r(g_spinor_field[DUM_MATRIX], k_strange, -g_epsbar, VOLUME/2);
  assign_add_mul_r(g_spinor_field[DUM_MATRIX+1], k_charm, -g_epsbar, VOLUME/2);


  diff(l_strange, g_spinor_field[DUM_MATRIX], l_strange, VOLUME/2);
  diff(l_charm, g_spinor_field[DUM_MATRIX+1], l_charm, VOLUME/2);


  /* and finally the  gamma_5  multiplication  */
  gamma5(l_strange, l_strange, VOLUME/2);
  gamma5(l_charm, l_charm, VOLUME/2);

  /* At the end, the normalisation by the max. eigenvalue  */
  mul_r(l_strange, phmc_invmaxev, l_strange, VOLUME/2);
  mul_r(l_charm, phmc_invmaxev, l_charm, VOLUME/2);

  /*     
  printf(" IN UP: %f %f \n", l_strange[0].creal(s2.c1), l_strange[0].cimag(s2.c1));
  printf(" IN DN: %f %f \n", l_charm[0].creal(s2.c1), l_charm[0].cimag(s2.c1));
  */

  /*  AND FINALLY WE SUBSTRACT THE C-CONSTANT  */


  /************ loop over all lattice sites ************/
  for(ix = 0; ix < (VOLUME/2); ix++){

    r=l_strange + ix;
    s=k_strange + ix;

    _complex_times_vector(phi1, z, s->s0);
    _vector_sub_assign(r->s0, phi1);
    _complex_times_vector(phi1, z, s->s1);
    _vector_sub_assign(r->s1, phi1);
    _complex_times_vector(phi1, z, s->s2);
    _vector_sub_assign(r->s2, phi1);
    _complex_times_vector(phi1, z, s->s3);
    _vector_sub_assign(r->s3, phi1);
    
    r=l_charm + ix;
    s=k_charm + ix;

    _complex_times_vector(phi1, z, s->s0);
    _vector_sub_assign(r->s0, phi1);
    _complex_times_vector(phi1, z, s->s1);
    _vector_sub_assign(r->s1, phi1);
    _complex_times_vector(phi1, z, s->s2);
    _vector_sub_assign(r->s2, phi1);
    _complex_times_vector(phi1, z, s->s3);
    _vector_sub_assign(r->s3, phi1);    
  }

  /* Finally, we multiply by the constant  phmc_Cpol  */
  /* which renders the polynomial in monomials  */
  /* identical to the polynomial a la clenshaw */;
  mul_r(l_strange, phmc_Cpol, l_strange, VOLUME/2);
  mul_r(l_charm, phmc_Cpol, l_charm, VOLUME/2);

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
void Q_Qdagger_ND_BI(bispinor * const bisp_l, bispinor * const bisp_k){

  double nrm = 1./(1.+g_mubar*g_mubar-g_epsbar*g_epsbar);
  static int memalloc = 0;

  static spinor *k_strange, *k_charm;
  static spinor *l_strange, *l_charm;

#if ( defined SSE || defined SSE2 || defined SSE3)
  static spinor *k_strange_, *k_charm_;
  static spinor *l_strange_, *l_charm_;
  if(memalloc == 0) {
    memalloc = 1;
    k_strange_ = (spinor*)calloc(VOLUMEPLUSRAND/2+1, sizeof(spinor));
    k_strange  = (spinor *)(((unsigned long int)(k_strange_)+ALIGN_BASE)&~ALIGN_BASE);
    k_charm_   = (spinor*)calloc(VOLUMEPLUSRAND/2+1, sizeof(spinor));
    k_charm    = (spinor *)(((unsigned long int)(k_charm_)+ALIGN_BASE)&~ALIGN_BASE);

    l_strange_ = (spinor*)calloc(VOLUMEPLUSRAND/2+1, sizeof(spinor));
    l_strange  = (spinor *)(((unsigned long int)(l_strange_)+ALIGN_BASE)&~ALIGN_BASE);
    l_charm_   = (spinor*)calloc(VOLUMEPLUSRAND/2+1, sizeof(spinor));
    l_charm    = (spinor *)(((unsigned long int)(l_charm_)+ALIGN_BASE)&~ALIGN_BASE);
  }
#else
  if(memalloc == 0) {
    memalloc = 1;
    k_strange  = (spinor*)calloc(VOLUMEPLUSRAND/2, sizeof(spinor));
    k_charm    = (spinor*)calloc(VOLUMEPLUSRAND/2, sizeof(spinor));
    
    l_strange  = (spinor*)calloc(VOLUMEPLUSRAND/2, sizeof(spinor));
    l_charm    = (spinor*)calloc(VOLUMEPLUSRAND/2, sizeof(spinor));
  }
#endif

  /*  CREATE 2 SPINORS OUT OF 1 (INPUT) BISPINOR  */
  decompact(k_strange, k_charm, bisp_k);

  /* FIRST THE  Qhat(2x2)^dagger  PART*/

  /* Here the  M_oe Mee^-1 M_eo  implementation  */
  Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX], k_charm);
  Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX+1], k_strange);

  mul_one_minus_imubar(g_spinor_field[DUM_MATRIX+2], g_spinor_field[DUM_MATRIX]);
  mul_one_plus_imubar(g_spinor_field[DUM_MATRIX+3], g_spinor_field[DUM_MATRIX+1]);

  assign_add_mul_r(g_spinor_field[DUM_MATRIX+2], g_spinor_field[DUM_MATRIX+1], g_epsbar, VOLUME/2);
  assign_add_mul_r(g_spinor_field[DUM_MATRIX+3], g_spinor_field[DUM_MATRIX], g_epsbar, VOLUME/2);

  mul_r(g_spinor_field[DUM_MATRIX+2], nrm, g_spinor_field[DUM_MATRIX+2], VOLUME/2);
  mul_r(g_spinor_field[DUM_MATRIX+3], nrm, g_spinor_field[DUM_MATRIX+3], VOLUME/2);
  /* where nrm (= 1/(1+mu^2 -eps^2)) has been defined at the beginning of 
     the subroutine */

  Hopping_Matrix(OE, g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+2]);
  Hopping_Matrix(OE, g_spinor_field[DUM_MATRIX+1], g_spinor_field[DUM_MATRIX+3]);

  /* Here the M_oo  implementation  */
  mul_one_plus_imubar(g_spinor_field[DUM_MATRIX+2], k_charm);
  mul_one_minus_imubar(g_spinor_field[DUM_MATRIX+3], k_strange);

  assign_add_mul_r(g_spinor_field[DUM_MATRIX+2], k_strange, -g_epsbar, VOLUME/2);
  assign_add_mul_r(g_spinor_field[DUM_MATRIX+3], k_charm, -g_epsbar, VOLUME/2);
   
  diff(g_spinor_field[DUM_MATRIX+4], g_spinor_field[DUM_MATRIX+2], g_spinor_field[DUM_MATRIX], VOLUME/2);
  diff(g_spinor_field[DUM_MATRIX+5], g_spinor_field[DUM_MATRIX+3], g_spinor_field[DUM_MATRIX+1], VOLUME/2);

  /* and finally the  gamma_5  multiplication  */
  gamma5(g_spinor_field[DUM_MATRIX+2], g_spinor_field[DUM_MATRIX+4], VOLUME/2);
  gamma5(g_spinor_field[DUM_MATRIX+3], g_spinor_field[DUM_MATRIX+5], VOLUME/2);

  /* The normalisation by the max. eigenvalue  is done twice at the end */


  /* We have to reassigin as follows to avoid overwriting */
  /* Recall in fact that   Q^hat = tau_1 Q tau_1  , hence  */

  /*  ABOVE: dum_matrix+2  is  l_charm   goes to  dum_matrix+6 :BELOW */
  /*  ABOVE: dum_matrix+3  is  l_strange   goes to  dum_matrix+7 :BELOW */
  assign(g_spinor_field[DUM_MATRIX+6], g_spinor_field[DUM_MATRIX+2], VOLUME/2);
  assign(g_spinor_field[DUM_MATRIX+7], g_spinor_field[DUM_MATRIX+3], VOLUME/2);


  /* AND THEN THE  Qhat(2x2)  PART */

  /* Here the  M_oe Mee^-1 M_eo  implementation  */
  Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+7]);
  Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX+1], g_spinor_field[DUM_MATRIX+6]);

  mul_one_minus_imubar(g_spinor_field[DUM_MATRIX+2], g_spinor_field[DUM_MATRIX]);
  mul_one_plus_imubar(g_spinor_field[DUM_MATRIX+3], g_spinor_field[DUM_MATRIX+1]);

  assign_add_mul_r(g_spinor_field[DUM_MATRIX+2], g_spinor_field[DUM_MATRIX+1], g_epsbar, VOLUME/2);
  assign_add_mul_r(g_spinor_field[DUM_MATRIX+3], g_spinor_field[DUM_MATRIX], g_epsbar, VOLUME/2);

  mul_r(g_spinor_field[DUM_MATRIX+2], nrm, g_spinor_field[DUM_MATRIX+2], VOLUME/2);
  mul_r(g_spinor_field[DUM_MATRIX+3], nrm, g_spinor_field[DUM_MATRIX+3], VOLUME/2);
  /* where nrm (= 1/(1+mu^2 -eps^2)) has been defined at the beginning of 
     the subroutine */

  Hopping_Matrix(OE, l_strange, g_spinor_field[DUM_MATRIX+2]);
  Hopping_Matrix(OE, l_charm, g_spinor_field[DUM_MATRIX+3]);

  /* Here the M_oo  implementation  */
  mul_one_plus_imubar(g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+7]);
  mul_one_minus_imubar(g_spinor_field[DUM_MATRIX+1], g_spinor_field[DUM_MATRIX+6]);

  assign_add_mul_r(g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+6], -g_epsbar, VOLUME/2);
  assign_add_mul_r(g_spinor_field[DUM_MATRIX+1], g_spinor_field[DUM_MATRIX+7], -g_epsbar, VOLUME/2);
   
  diff(l_strange, g_spinor_field[DUM_MATRIX], l_strange, VOLUME/2);
  diff(l_charm, g_spinor_field[DUM_MATRIX+1], l_charm, VOLUME/2);

  /* and finally the  gamma_5  multiplication  */
  gamma5(l_strange, l_strange, VOLUME/2);
  gamma5(l_charm, l_charm, VOLUME/2);

  /* At the end, the normalisation by the max. eigenvalue  */
  /* Twice  phmc_invmaxev  since we consider here  D Ddag  !!! */
  mul_r(l_charm, phmc_invmaxev*phmc_invmaxev, l_charm, VOLUME/2);
  mul_r(l_strange, phmc_invmaxev*phmc_invmaxev, l_strange, VOLUME/2);


  /*  CREATE 1 (OUTPUT) BISPINOR OUT OF 2 SPINORS  */
  compact(bisp_l, l_strange, l_charm);
}


/******************************************
 *
 * This is the implementation of
 *
 * (M_{ee}^\pm)^{-1}M_{eo}
 *
 * see documentation for details
 * k is the number of the input field
 * l is the number of the output field
 *
 * it acts only on the odd part or only 
 * on a half spinor
 ******************************************/
void H_eo_ND(spinor * const l_strange, spinor * const l_charm, 
             spinor * const k_strange, spinor * const k_charm, 
	     const int ieo) {

  double nrm = 1./(1.+g_mubar*g_mubar-g_epsbar*g_epsbar);


  /* recall:   strange <-> up    while    charm <-> dn   */
  Hopping_Matrix(ieo, g_spinor_field[DUM_MATRIX], k_strange);
  Hopping_Matrix(ieo, g_spinor_field[DUM_MATRIX+1], k_charm);

  mul_one_minus_imubar(l_strange, g_spinor_field[DUM_MATRIX+1]);
  mul_one_plus_imubar(l_charm, g_spinor_field[DUM_MATRIX]);

  assign_add_mul_r(l_strange, g_spinor_field[DUM_MATRIX], g_epsbar, VOLUME/2);
  assign_add_mul_r(l_charm, g_spinor_field[DUM_MATRIX+1], g_epsbar, VOLUME/2);

  mul_r(l_strange, nrm, l_strange, VOLUME/2);
  mul_r(l_charm, nrm, l_charm, VOLUME/2);

}

void M_ee_inv_ND(spinor * const l_strange, spinor * const l_charm, 
		 spinor * const k_strange, spinor * const k_charm) {
  
  double nrm = 1./(1.+g_mubar*g_mubar-g_epsbar*g_epsbar);


  /* recall:   strange <-> up    while    charm <-> dn   */

  mul_one_minus_imubar(l_strange, k_strange);
  mul_one_plus_imubar(l_charm, k_charm);

  assign_add_mul_r(l_strange, k_charm, g_epsbar, VOLUME/2);
  assign_add_mul_r(l_charm, k_strange, g_epsbar, VOLUME/2);

  mul_r(l_strange, nrm, l_strange, VOLUME/2);
  mul_r(l_charm, nrm, l_charm, VOLUME/2);

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


void mul_one_minus_imubar(spinor * const l, spinor * const k)
{
  spinor *r, *s;
  static su3_vector phi1;

  /************ loop over all lattice sites ************/
  for(int ix = 0; ix < (VOLUME/2); ++ix){
    r=l + ix;
    s=k + ix;
    /* Multiply the spinorfield with the inverse of 1+imu\gamma_5 */
    _complex_times_vector(phi1, (1. - g_mubar * I), s->s0);
    _vector_assign(r->s0, phi1);
    _complex_times_vector(phi1, (1. - g_mubar * I), s->s1);
    _vector_assign(r->s1, phi1);
    _complex_times_vector(phi1, (1. + g_mubar * I), s->s2);
    _vector_assign(r->s2, phi1);
    _complex_times_vector(phi1, (1. + g_mubar * I), s->s3);
    _vector_assign(r->s3, phi1);
  }
}


void mul_one_plus_imubar(spinor * const l, spinor * const k){
  spinor *r, *s;
  static su3_vector phi1;

  /************ loop over all lattice sites ************/
  for(int ix = 0; ix < (VOLUME/2); ++ix){
    r=l + ix;
    s=k + ix;
    /* Multiply the spinorfield with the inverse of 1+imu\gamma_5 */
    _complex_times_vector(phi1, (1. + g_mubar * I), s->s0);
    _vector_assign(r->s0, phi1);
    _complex_times_vector(phi1, (1. + g_mubar * I), s->s1);
    _vector_assign(r->s1, phi1);
    _complex_times_vector(phi1, (1. - g_mubar * I), s->s2);
    _vector_assign(r->s2, phi1);
    _complex_times_vector(phi1, (1. - g_mubar * I), s->s3);
    _vector_assign(r->s3, phi1);
  }
  return;
}


/*  calculates P(Q Q^dagger) for the nondegenerate case */

void P_ND(spinor * const l_strange, spinor * const l_charm,
	  spinor * const k_strange, spinor * const k_charm){
  
  
  
  int j;
  spinor *dum_up,*dum_dn;
  dum_up=g_chi_up_spinor_field[DUM_MATRIX];
  dum_dn=g_chi_dn_spinor_field[DUM_MATRIX];
  
  assign(dum_up,k_strange,VOLUME/2);
  assign(dum_dn,k_charm,VOLUME/2);
  
  
  for(j=0; j<(2*phmc_dop_n_cheby -2); j++){
    if(j>0) {
      assign(dum_up,l_strange,VOLUME/2);
      assign(dum_dn,l_charm,VOLUME/2);
    }
    
    Q_tau1_min_cconst_ND(l_strange, l_charm,
			 dum_up, dum_dn,
			 phmc_root[j]);
  }
  return;
}


/* calculates  Q * \tau^1  for the nondegenerate case */
void Qtau1_P_ND(spinor * const l_strange, spinor * const l_charm,
		spinor * const k_strange, spinor * const k_charm){
  
  
  spinor * dum_up,* dum_dn;
  dum_up=g_chi_up_spinor_field[DUM_MATRIX+1];
  dum_dn=g_chi_dn_spinor_field[DUM_MATRIX+1];
  
  P_ND(l_strange, l_charm,k_strange,k_charm);
  
  assign(dum_up,l_strange,VOLUME/2);
  assign(dum_dn,l_charm,VOLUME/2);
  
  QNon_degenerate(l_strange,l_charm,dum_dn,dum_up);
  return;
}



/* this is neccessary for the calculation of the polynomial */

void Qtm_pm_min_cconst_nrm(spinor * const l, spinor * const k,
			   const _Complex double z){
  static su3_vector phi1;
  spinor *r,*s;
  int ix;
  Qtm_pm_psi(l,k);
  mul_r(l, phmc_invmaxev, l, VOLUME/2);

  /*  AND FINALLY WE SUBSTRACT THE C-CONSTANT  */


  /************ loop over all lattice sites ************/
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
    
    Qtm_pm_min_cconst_nrm(l,spinDum,phmc_root[j]);
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
  static su3_vector phi;
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







/**************************************************************
 * $Id$ *
 *                                                            *
 * This file contains operators for twisted mass Wilson QCD   *
 * to construct a multiplication with a non-degenerate        *
 * flavour matrix                                             *
 *                                                            *
 * see notes of R. Frezzoti and T. Chiarappa for more details *
 * Author: Karl Jansen                                        *
 *         Karl.Jansen@desy.de                                *
 **************************************************************/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include "global.h"
#include "su3.h"
#include "Hopping_Matrix.h"
/* #include "sse.h" */
#include "gamma.h"
/* in piu` */
#include "linsolve.h"
/* fine in piu` */
#include "linalg_eo.h"
#include "Nondegenerate_Matrix.h"

/* internal */


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

  /*  double invmaxev=1./20.;*/
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
  /*
  mul_r(l_strange, invmaxev, l_strange, VOLUME/2);
  mul_r(l_charm, invmaxev, l_charm, VOLUME/2);
  */
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

  /*  double invmaxev=1./20.;*/
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
  /*
  mul_r(l_charm, invmaxev, l_charm, VOLUME/2);
  mul_r(l_strange, invmaxev, l_strange, VOLUME/2);
  */

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

  /*  double invmaxev=1./20.;*/
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

  /* At the end, the normalisation by the max. eigenvalue  */
  /*
  mul_r(k_charm, invmaxev, g_spinor_field[DUM_MATRIX+2], VOLUME/2);
  mul_r(k_strange, invmaxev, g_spinor_field[DUM_MATRIX+3], VOLUME/2);
  */


  /* AND THEN THE  Qhat(2x2)  PART */

  /* Here the  M_oe Mee^-1 M_eo  implementation  */
  Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX], k_strange);
  Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX+1], k_charm);

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
  /*
  mul_r(l_charm, invmaxev, l_charm, VOLUME/2);
  mul_r(l_strange, invmaxev, l_strange, VOLUME/2);
  */

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

  /*  double invmaxev=1./20.;*/
  double nrm = 1./(1.+g_mubar*g_mubar-g_epsbar*g_epsbar);


  spinor *k_strange, *k_charm, *k_strange_, *k_charm_;
  spinor *l_strange, *l_charm, *l_strange_, *l_charm_;

#if ( defined SSE || defined SSE2 || defined SSE3)
  k_strange_ = calloc(VOLUMEPLUSRAND/2+1, sizeof(spinor));
  k_strange  = (spinor *)(((unsigned int)(k_strange_)+ALIGN_BASE)&~ALIGN_BASE);
  k_charm_   = calloc(VOLUMEPLUSRAND/2+1, sizeof(spinor));
  k_charm    = (spinor *)(((unsigned int)(k_charm_)+ALIGN_BASE)&~ALIGN_BASE);

  l_strange_ = calloc(VOLUMEPLUSRAND/2+1, sizeof(spinor));
  l_strange  = (spinor *)(((unsigned int)(l_strange_)+ALIGN_BASE)&~ALIGN_BASE);
  l_charm_   = calloc(VOLUMEPLUSRAND/2+1, sizeof(spinor));
  l_charm    = (spinor *)(((unsigned int)(l_charm_)+ALIGN_BASE)&~ALIGN_BASE);
#else
  k_strange  =calloc(VOLUMEPLUSRAND/2, sizeof(spinor));
  k_charm    =calloc(VOLUMEPLUSRAND/2, sizeof(spinor));

  l_strange  =calloc(VOLUMEPLUSRAND/2, sizeof(spinor));
  l_charm    =calloc(VOLUMEPLUSRAND/2, sizeof(spinor));
#endif

  /*  CREATE 2 SPINORS OUT OF 1 (INPUT) BISPINOR  */
  decompact(&k_strange[0], &k_charm[0], &bisp_k[0]);


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

  /* At the end, the normalisation by the max. eigenvalue  */
  /*  
  mul_r(k_charm, invmaxev, g_spinor_field[DUM_MATRIX+2], VOLUME/2);
  mul_r(k_strange, invmaxev, g_spinor_field[DUM_MATRIX+3], VOLUME/2);
  */
  

  /*    !!!  I HAVE REPLACED THE ABOVE LINES BY ASSIGNMENTS  !!!     */
  assign(k_charm, g_spinor_field[DUM_MATRIX+2], VOLUME/2);
  assign(k_strange, g_spinor_field[DUM_MATRIX+3], VOLUME/2);


  /* AND THEN THE  Qhat(2x2)  PART */

  /* Here the  M_oe Mee^-1 M_eo  implementation  */
  Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX], k_strange);
  Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX+1], k_charm);

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
  /*
  mul_r(l_charm, invmaxev, l_charm, VOLUME/2);
  mul_r(l_strange, invmaxev, l_strange, VOLUME/2);
  */  

  /*    !!!  I HAVE REPLACED THE ABOVE LINES BY ASSIGNMENTS  !!!     */


  /*  CREATE 1 (OUTPUT) BISPINOR OUT OF 2 SPINORS  */
  compact(&bisp_l[0], &l_strange[0], &l_charm[0]);


#if ( defined SSE || defined SSE2 || defined SSE3)
  free(k_strange_);
  free(k_charm_);
  free(l_strange_);
  free(l_charm_);
#else
  free(k_strange);
  free(k_charm);
  free(l_strange);
  free(l_charm);
#endif  


}


void Q_test_epsilon(spinor * const l_strange, spinor * const l_charm,
                           spinor * const k_strange, spinor * const k_charm){

  /*  double invmaxev=1./20.; */
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
  /*
  mul_r(l_charm, invmaxev, l_charm, VOLUME/2);
  mul_r(l_strange, invmaxev, l_strange, VOLUME/2);
  */

}


void mul_one_minus_imubar(spinor * const l, spinor * const k){
  complex z,w;
  int ix;
  spinor *r, *s;
  static su3_vector phi1;

  z.re = 1.;
  z.im = - g_mubar;
  w.re = 1.;
  w.im = - z.im; 

  /************ loop over all lattice sites ************/
  for(ix = 0; ix < (VOLUME/2); ix++){
    r=l + ix;
    s=k + ix;
    /* Multiply the spinorfield with the inverse of 1+imu\gamma_5 */
    _complex_times_vector(phi1, z, (*s).s0);
    _vector_assign((*r).s0, phi1);
    _complex_times_vector(phi1, z, (*s).s1);
    _vector_assign((*r).s1, phi1);
    _complex_times_vector(phi1, w, (*s).s2);
    _vector_assign((*r).s2, phi1);
    _complex_times_vector(phi1, w, (*s).s3);
    _vector_assign((*r).s3, phi1);
  }
}


void mul_one_plus_imubar(spinor * const l, spinor * const k){
  complex z,w;
  int ix;
  spinor *r, *s;
  static su3_vector phi1;

  z.re = 1.;
  z.im = g_mubar;
  w.re = 1.;
  w.im =  -z.im; 

  /************ loop over all lattice sites ************/
  for(ix = 0; ix < (VOLUME/2); ix++){
    r=l + ix;
    s=k + ix;
    /* Multiply the spinorfield with the inverse of 1+imu\gamma_5 */
    _complex_times_vector(phi1, z, (*s).s0);
    _vector_assign((*r).s0, phi1);
    _complex_times_vector(phi1, z, (*s).s1);
    _vector_assign((*r).s1, phi1);
    _complex_times_vector(phi1, w, (*s).s2);
    _vector_assign((*r).s2, phi1);
    _complex_times_vector(phi1, w, (*s).s3);
    _vector_assign((*r).s3, phi1);
  }
}



static char const rcsid[] = "$Id$";

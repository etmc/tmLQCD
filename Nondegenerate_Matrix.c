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

#include <stdlib.h>
#include <stdio.h>
#include "global.h"
#include "su3.h"
#include "Hopping_Matrix.h"
/* #include "sse.h" */
#include "gamma.h"
#include "linalg_eo.h"
#include "Nondegenerate_Matrix.h"

/* internal */

/******************************************
 * mul_one_minus_imubar_inv computes
 * l = [(1-i\mubar\gamma_5)/(1+mubar^2-epsbar^2)] * l
 *
 ******************************************/
void mul_one_minus_imubar_inv(spinor * const l, spinor * const k);
/******************************************
 * mul_one_plus_imubar_inv computes
 * l = [(1+i\mubar\gamma_5)/(1+mubar^2-epsbar^2)] * l
 *
 ******************************************/
void mul_one_plus_imubar_inv(spinor * const l, spinor * const k);
/******************************************
 * mul_one_plus_imubar_inv computes
 * l = [(1+i\mubar\gamma_5)/(1+mubar^2-epsbar^2)] * l
 *
*/


/* "mul_one_minus_imubar_inv"  and  "mul_one_minus_imubar_inv" 
       NO LONGER USED 

      REPLACED BY    "mul_one_minus_imubar"  and  "mul_one_minus_imubar"  
*/


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

  double temp;
  double invmaxev=1./20.; 
  double nrm = 1./(1.+g_mubar*g_mubar-g_epsbar*g_epsbar);

  Hopping_Matrix(EO, spinor_field[DUM_MATRIX], k_strange);
  Hopping_Matrix(EO, spinor_field[DUM_MATRIX+1], k_charm);


  /* Even if in double precision has almost no meaning, according to 
     Roberto`s suggestion (as far as he said about its and Karl`s 
     experience on PHMC), in order to avoid possible roundoff it should be 
     better to postpone the division by  (1+mu^2 -eps^2)  of M_ee^1 later*/

  /* Therefore, I commented out the following 2 lines (your own)*/

  /*  mul_one_minus_imubar_inv(spinor_field[DUM_MATRIX+2], spinor_field[DUM_MATRIX]); */
  /*  mul_one_plus_imubar_inv(spinor_field[DUM_MATRIX+3], spinor_field[DUM_MATRIX+1]);*/

  /* and include the following 4 (2 here and 2 at line 102) new lines */

    /* DUM+2 Substitued by DUM+4 */
  mul_one_minus_imubar(spinor_field[DUM_MATRIX+4], spinor_field[DUM_MATRIX]);
  mul_one_plus_imubar(spinor_field[DUM_MATRIX+3], spinor_field[DUM_MATRIX+1]);

    /* DUM+2 Substitued by DUM+4 */
  assign_add_mul_r(spinor_field[DUM_MATRIX+4], spinor_field[DUM_MATRIX+1], g_epsbar, VOLUME/2);
  assign_add_mul_r(spinor_field[DUM_MATRIX+3], spinor_field[DUM_MATRIX], g_epsbar, VOLUME/2);

    /* DUM+2 Substitued by DUM+4 */
  mul_r(spinor_field[DUM_MATRIX+4], nrm, spinor_field[DUM_MATRIX+4], VOLUME/2);
  mul_r(spinor_field[DUM_MATRIX+3], nrm, spinor_field[DUM_MATRIX+3], VOLUME/2);
  /* where nrm (= 1/(1+mu^2 -eps^2)) has been defined at the beginning of 
     the subroutine */
  
    /* DUM+2 Substitued by DUM+4 */
  Hopping_Matrix(OE, l_strange, spinor_field[DUM_MATRIX+4]);
  Hopping_Matrix(OE, l_charm, spinor_field[DUM_MATRIX+3]);


  /* Here in my opinion  M_oo  is not implemented correctly.
     I propose the following  */

  mul_one_plus_imubar(spinor_field[DUM_MATRIX], k_strange);
  mul_one_minus_imubar(spinor_field[DUM_MATRIX+1], k_charm);

  assign_add_mul_r(spinor_field[DUM_MATRIX], k_charm, -g_epsbar, VOLUME/2);
  assign_add_mul_r(spinor_field[DUM_MATRIX+1], k_strange, -g_epsbar, VOLUME/2);
   
 
  /* Therefore I have also slightly modified the  diff  command*/
  diff(l_strange, spinor_field[DUM_MATRIX], l_strange, VOLUME/2);
  diff(l_charm, spinor_field[DUM_MATRIX+1], l_charm, VOLUME/2);

  gamma5(l_strange, l_strange, VOLUME/2);
  gamma5(l_charm, l_charm, VOLUME/2);

  mul_r(l_strange, invmaxev, l_strange, VOLUME/2);
  mul_r(l_charm, invmaxev, l_charm, VOLUME/2);


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

  double invmaxev=1./20.; 

  double nrm = 1./(1.+g_mubar*g_mubar-g_epsbar*g_epsbar);

  Hopping_Matrix(EO, spinor_field[DUM_MATRIX], k_charm);
  Hopping_Matrix(EO, spinor_field[DUM_MATRIX+1], k_strange);


    /* DUM+2 Substitued by DUM+4 */
  mul_one_minus_imubar(spinor_field[DUM_MATRIX+4], spinor_field[DUM_MATRIX]);
  mul_one_plus_imubar(spinor_field[DUM_MATRIX+3], spinor_field[DUM_MATRIX+1]);

    /* DUM+2 Substitued by DUM+4 */
  assign_add_mul_r(spinor_field[DUM_MATRIX+4], spinor_field[DUM_MATRIX+1], g_epsbar, VOLUME/2);
  assign_add_mul_r(spinor_field[DUM_MATRIX+3], spinor_field[DUM_MATRIX], g_epsbar, VOLUME/2);

    /* DUM+2 Substitued by DUM+4 */
  mul_r(spinor_field[DUM_MATRIX+4], nrm, spinor_field[DUM_MATRIX+4], VOLUME/2);
  mul_r(spinor_field[DUM_MATRIX+3], nrm, spinor_field[DUM_MATRIX+3], VOLUME/2);
  /* where nrm (= 1/(1+mu^2 -eps^2)) has been defined at the beginning of 
     the subroutine */
  

    /* DUM+2 Substitued by DUM+4 */
  Hopping_Matrix(OE, l_strange, spinor_field[DUM_MATRIX+4]);
  Hopping_Matrix(OE, l_charm, spinor_field[DUM_MATRIX+3]);

  /* Here in my opinion  M_oo  is not implemented correctly.
     I propose the following  */
  mul_one_plus_imubar(spinor_field[DUM_MATRIX], k_charm);
  mul_one_minus_imubar(spinor_field[DUM_MATRIX+1], k_strange);

  assign_add_mul_r(spinor_field[DUM_MATRIX], k_strange, -g_epsbar, VOLUME/2);
  assign_add_mul_r(spinor_field[DUM_MATRIX+1], k_charm, -g_epsbar, VOLUME/2);
   
 
  /* Therefore I have also slightly modified the  diff  command*/
  diff(spinor_field[DUM_MATRIX+3], spinor_field[DUM_MATRIX], l_strange, VOLUME/2);
  diff(spinor_field[DUM_MATRIX+4], spinor_field[DUM_MATRIX+1], l_charm, VOLUME/2);


    /* DUM+2 Substitued by DUM+4 */
  gamma5(l_charm, spinor_field[DUM_MATRIX+3], VOLUME/2);
  gamma5(l_strange, spinor_field[DUM_MATRIX+4], VOLUME/2);

    /* DUM+2 Substitued by DUM+4 */
  mul_r(l_strange, invmaxev, l_strange, VOLUME/2);
  mul_r(l_charm, invmaxev, l_charm, VOLUME/2);
  
}



/*      subroutines: 
   "mul_one_minus_imubar_inv"  and  "mul_one_minus_imubar_inv" 
       NO LONGER USED  ....   REPLACED BY    
   "mul_one_minus_imubar"  and  "mul_one_minus_imubar"  
*/



void mul_one_minus_imubar(spinor * const l, spinor * const k){
  complex z,w;
  int ix;
  spinor *r, *s;
  static su3_vector phi1;

  z.re = 1.;
  z.im = - g_mubar;
  w.re = 1.;
  w.im = z.im; 

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



void mul_one_minus_imubar_inv(spinor * const l, spinor * const k){
  complex z,w;
  int ix;
  spinor *r, *s;
  static su3_vector phi1;
  double nrm = 1./(1.+g_mubar*g_mubar-g_epsbar*g_epsbar);

  z.re = nrm;
  z.im =  -nrm * g_mubar;
  w.re = nrm;
  w.im =  z.im; 


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

void mul_one_plus_imubar_inv(spinor * const l, spinor * const k){
  complex z,w;
  int ix;
  spinor *r, *s;
  static su3_vector phi1;
  double nrm = 1./(1.+g_mubar*g_mubar-g_epsbar*g_epsbar);

  z.re = nrm;
  z.im =   nrm * g_mubar;
  w.re = nrm;
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
/******************************************
 *
 * This is the implementation of
 *
 *  [ M_ee^{-1} M_eo ]
 *
 * to reconstruct the even sites for the force
 ******************************************/
void QNon_degenerate_eo(spinor * const l_strange, spinor * const l_charm,
                        spinor * const k_strange, spinor * const k_charm){

  double nrm = 1./(1.+g_mubar*g_mubar-g_epsbar*g_epsbar);


  Hopping_Matrix(EO, spinor_field[DUM_MATRIX], k_strange);
  Hopping_Matrix(EO, spinor_field[DUM_MATRIX+1], k_charm);

  /* inverse? Carsten thinks _inv */
  mul_one_minus_imubar_inv(spinor_field[DUM_MATRIX+4], spinor_field[DUM_MATRIX]);
  mul_one_plus_imubar_inv(spinor_field[DUM_MATRIX+3], spinor_field[DUM_MATRIX+1]);

  assign_add_mul_r(spinor_field[DUM_MATRIX+4], spinor_field[DUM_MATRIX+1], g_epsbar, VOLUME/2);
  assign_add_mul_r(spinor_field[DUM_MATRIX+3], spinor_field[DUM_MATRIX], g_epsbar, VOLUME/2);

  mul_r(spinor_field[DUM_MATRIX+4], nrm, spinor_field[DUM_MATRIX+4], VOLUME/2);
  mul_r(spinor_field[DUM_MATRIX+3], nrm, spinor_field[DUM_MATRIX+3], VOLUME/2);

  assign(l_strange, spinor_field[DUM_MATRIX+4], VOLUME/2);
  assign(l_charm, spinor_field[DUM_MATRIX+3], VOLUME/2);
}


static char const rcsid[] = "$Id$";

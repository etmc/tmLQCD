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
void mul_one_minus_imubar_inv(spinor * const l);
/******************************************
 * mul_one_plus_imubar_inv computes
 * l = [(1+i\mubar\gamma_5)/(1+mubar^2-epsbar^2)] * l
 *
 ******************************************/
void mul_one_minus_imubar_inv(spinor * const l);
/******************************************
 * mul_one_plus_imubar_inv computes
 * l = [(1+i\mubar\gamma_5)/(1+mubar^2-epsbar^2)] * l
 *


/* external functions */

/******************************************
 *
 * This is the implementation of
 *
 * MISSING
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

  Hopping_Matrix(EO, spinor_field[DUM_MATRIX], k_strange);
  Hopping_Matrix(EO, spinor_field[DUM_MATRIX+1], k_charm);

  mul_one_minus_imubar_inv(spinor_field[DUM_MATRIX+2], spinor_field[DUM_MATRIX]);
  mul_one_plus_imubar_inv(spinor_field[DUM_MATRIX+3], spinor_field[DUM_MATRIX+1]);

  assign_add_mul_r(spinor_field[DUM_MATRIX+2], spinor_field[DUM_MATRIX+1], g_epsbar, VOLUME/2);
  assign_add_mul_r(spinor_field[DUM_MATRIX+3], spinor_field[DUM_MATRIX], g_epsbar, VOLUME/2);

  Hopping_Matrix(OE, l_strange, spinor_field[DUM_MATRIX+2]);
  Hopping_Matrix(OE, l_charm, spinor_field[DUM_MATRIX+3]);

  gamma5(l_strange, l_strange, VOLUME/2);
  gamma5(l_charm, l_charm, VOLUME/2);

}

/******************************************
 *
 * This is the implementation of
 *
 * MISSING
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

  Hopping_Matrix(EO, spinor_field[DUM_MATRIX], k_charm);
  Hopping_Matrix(EO, spinor_field[DUM_MATRIX+1], k_strange);

  mul_one_minus_imubar_inv(spinor_field[DUM_MATRIX+2], spinor_field[DUM_MATRIX]);
  mul_one_plus_imubar_inv(spinor_field[DUM_MATRIX+3], spinor_field[DUM_MATRIX+1]);

  assign_add_mul_r(spinor_field[DUM_MATRIX+2], spinor_field[DUM_MATRIX+1], g_epsbar, VOLUME/2);
  assign_add_mul_r(spinor_field[DUM_MATRIX+3], spinor_field[DUM_MATRIX], g_epsbar, VOLUME/2);

  Hopping_Matrix(OE, spinor_field[DUM_MATRIX], spinor_field[DUM_MATRIX+2]);
  Hopping_Matrix(OE, spinor_field[DUM_MATRIX+1], spinor_field[DUM_MATRIX+3]);

  gamma5(l_charm, spinor_field[DUM_MATRIX], VOLUME/2);
  gamma5(l_stange, spinor_field[DUM_MATRIX+1], VOLUME/2);

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

static char const rcsid[] = "$Id$";

/**************************************************************
 * $Id$ *
 *                                                            *
 * This file contains operators for twisted mass Wilson QCD   *
 * prepared for even odd preconditioning                      *
 *                                                            *
 * see documentation for details                              *
 * Author: Carsten Urbach                                     *
 *         urbach@physik.fu-berlin.de                         *
 **************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include "global.h"
#include "su3.h"
#include "Hopping_Matrix.h"
#include "sse.h"
#include "tm_operators.h"

/* internal */

/******************************************
 * mul_one_pm_imu_inv computes
 * l = (1\pm i\mu\gamma_5)^{-1} * l
 *
 * sign is the sign used in 
 *      1\pm i\mu\gamma_5
 * l is number of input and output field
 *
 ******************************************/
void mul_one_pm_imu_inv(const int l, const double _sign);
void mul_one_pm_imu(const int l, const double _sign);
/******************************************
 * mul_one_pm_imu_sub_mul_gamma5 computes
 * l = gamma_5*((1\pm i\mu\gamma_5)*k - j)
 *
 * l is the number of the output field
 * k and j the numbers of the input fields
 *
 * sign indicates which sign should be used
 * in 1\pm i\mu\gamma_5
 ******************************************/
void mul_one_pm_imu_sub_mul_gamma5(const int l, const int k, 
				   const int j, const double _sign);

/******************************************
 * mul_one_pm_imu_sub_mul computes
 * l = ((1\pm i\mu\gamma_5)*k - j)
 *
 * l is the number of the output field
 * k and j the numbers of the input fields
 *
 * sign indicates which sign should be used
 * in 1\pm i\mu\gamma_5
 ******************************************/
void mul_one_pm_imu_sub_mul(const int l, const int k,
			    const int j, const double _sign);

/* external functions */

/******************************************
 *
 * This is the implementation of
 *
 * \hat Q_{+} =
 * \gamma_5(M_{oo}^+ - M_{oe}(M_{ee}^+ )^{-1}M_{eo})
 *
 * see documentation for details
 * k is the number of the input field
 * l is the number of the output field
 *
 * it acts only on the odd part or only 
 * on a half spinor
 ******************************************/
void Qtm_plus_psi(const int l, const int k){
  Hopping_Matrix(EO, DUM_MATRIX+1, k);
  mul_one_pm_imu_inv(DUM_MATRIX+1, +1.);
  Hopping_Matrix(OE, DUM_MATRIX, DUM_MATRIX+1);
  mul_one_pm_imu_sub_mul_gamma5(l, k, DUM_MATRIX, +1.);
}

/******************************************
 *
 * This is the implementation of
 *
 * \hat Q_{-} =
 * \gamma_5(M_{oo}^- - M_{oe}(M_{ee}^- )^{-1}M_{eo})
 *
 * see documentation for details
 * k is the number of the input field
 * l is the number of the output field
 *
 * it acts only on the odd part or only 
 * on a half spinor
 ******************************************/
void Qtm_minus_psi(const int l, const int k){
  Hopping_Matrix(EO, DUM_MATRIX+1, k);
  mul_one_pm_imu_inv(DUM_MATRIX+1, -1.);
  Hopping_Matrix(OE, DUM_MATRIX, DUM_MATRIX+1);
  mul_one_pm_imu_sub_mul_gamma5(l, k, DUM_MATRIX, -1.);
}

/******************************************
 *
 * This is the implementation of
 *
 * \gamma_5 \hat Q_{+} =
 * (M_{oo}^+ - M_{oe}(M_{ee}^+ )^{-1}M_{eo})
 *
 * see documentation for details
 * k is the number of the input field
 * l is the number of the output field
 *
 * it acts only on the odd part or only 
 * on a half spinor
 ******************************************/
void Mtm_plus_psi(const int l, const int k){
  Hopping_Matrix(EO, DUM_MATRIX+1, k);
  mul_one_pm_imu_inv(DUM_MATRIX+1, +1.);
  Hopping_Matrix(OE, DUM_MATRIX, DUM_MATRIX+1);
  mul_one_pm_imu_sub_mul(l, k, DUM_MATRIX, +1.);
}

/******************************************
 *
 * This is the implementation of
 *
 * \gamma_5 \hat Q_{-} =
 * (M_{oo}^- - M_{oe}(M_{ee}^- )^{-1}M_{eo})
 *
 * see documentation for details
 * k is the number of the input field
 * l is the number of the output field
 *
 * it acts only on the odd part or only 
 * on a half spinor
 ******************************************/
void Mtm_minus_psi(const int l, const int k){
  Hopping_Matrix(EO, DUM_MATRIX+1, k);
  mul_one_pm_imu_inv(DUM_MATRIX+1, -1.);
  Hopping_Matrix(OE, DUM_MATRIX, DUM_MATRIX+1);
  mul_one_pm_imu_sub_mul(l, k, DUM_MATRIX, -1.);
}

/******************************************
 *
 * This is the implementation of
 *
 * \hat Q_{+} \hat Q_{-} 
 *
 * see documentation for details
 * k is the number of the input field
 * l is the number of the output field
 *
 * it acts only on the odd part or only 
 * on a half spinor
 ******************************************/
void Qtm_pm_psi(const int l, const int k){
  /* Q_{-} */
  Hopping_Matrix(EO, DUM_MATRIX+1, k);
  mul_one_pm_imu_inv(DUM_MATRIX+1, -1.);
  Hopping_Matrix(OE, DUM_MATRIX, DUM_MATRIX+1);
  mul_one_pm_imu_sub_mul_gamma5(DUM_MATRIX, k, DUM_MATRIX, -1.);
  /* Q_{+} */
  Hopping_Matrix(EO, l, DUM_MATRIX);
  mul_one_pm_imu_inv(l, +1.);
  Hopping_Matrix(OE, DUM_MATRIX+1, l);
  mul_one_pm_imu_sub_mul_gamma5(l, DUM_MATRIX, DUM_MATRIX+1, +1.);
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
void H_eo_tm_inv_psi(const int l, const int k, 
		     const int ieo, const double sign){
  Hopping_Matrix(ieo, l, k);
  mul_one_pm_imu_inv(l, sign);
}

/**********************************************
 *
 * All the results are only stored in the first
 * half of the spinor fields, they have anly
 * lenght VOLUME/2
 * 
 * That's why mul_... do not need a iput
 * parameter ieo.
 *
 * the next functions are internal and you 
 * can find comments above at the declaration 
 *
 **********************************************/

void mul_one_pm_imu_inv(const int l, const double _sign){
  complex z,w;
  int ix;
  double sign=-1.; 
  spinor *r;
  static su3_vector phi1;
  double nrm = 1./(1.+g_mu*g_mu);

  if(_sign < 0.){
    sign = 1.; 
  }

  z.re = nrm;
  z.im =  sign * nrm * g_mu;
  w.re = nrm;
  w.im = -z.im; /*-sign * nrm * g_mu;*/

  /************ loop over all lattice sites ************/
  for(ix = 0; ix < (VOLUME/2); ix++){
    r=&spinor_field[l][ix];
    /* Multiply the spinorfield with the inverse of 1+imu\gamma_5 */
    _complex_times_vector(phi1, z, (*r).c1);
    _vector_assign((*r).c1, phi1);
    _complex_times_vector(phi1, z, (*r).c2);
    _vector_assign((*r).c2, phi1);
    _complex_times_vector(phi1, w, (*r).c3);
    _vector_assign((*r).c3, phi1);
    _complex_times_vector(phi1, w, (*r).c4);
    _vector_assign((*r).c4, phi1);
  }
}

void assign_mul_one_pm_imu_inv(const int l, const int k, const double _sign){
  complex z,w;
  int ix;
  double sign=-1.; 
  spinor *r, *s;
  double nrm = 1./(1.+g_mu*g_mu);

  if(_sign < 0.){
    sign = 1.; 
  }

  z.re = nrm;
  z.im =  sign * nrm * g_mu;
  w.re = nrm;
  w.im = -z.im; /*-sign * nrm * g_mu;*/

  /************ loop over all lattice sites ************/
  for(ix = 0; ix < (VOLUME/2); ix++){
    r=&spinor_field[k][ix];
    s=&spinor_field[l][ix];
    /* Multiply the spinorfield with the inverse of 1+imu\gamma_5 */
    _complex_times_vector((*s).c1, z, (*r).c1);
    _complex_times_vector((*s).c2, z, (*r).c2);
    _complex_times_vector((*s).c3, w, (*r).c3);
    _complex_times_vector((*s).c4, w, (*r).c4);
  }
}

void mul_one_pm_imu(const int l, const double _sign){
  complex z,w;
  int ix;
  double sign = 1.; 
  spinor *r;
  static su3_vector phi1;

  if(_sign < 0.){
    sign = -1.; 
  }

  z.re = 1.;
  z.im =  sign * g_mu;
  w.re = 1.;
  w.im = -z.im; /*-sign * nrm * g_mu;*/

  /************ loop over all lattice sites ************/
  for(ix = 0; ix < (VOLUME/2); ix++){
    r=&spinor_field[l][ix];
    /* Multiply the spinorfield with 1+imu\gamma_5 */
    _complex_times_vector(phi1, z, (*r).c1);
    _vector_assign((*r).c1, phi1);
    _complex_times_vector(phi1, z, (*r).c2);
    _vector_assign((*r).c2, phi1);
    _complex_times_vector(phi1, w, (*r).c3);
    _vector_assign((*r).c3, phi1);
    _complex_times_vector(phi1, w, (*r).c4);
    _vector_assign((*r).c4, phi1);
  }
}

void assign_mul_one_pm_imu(const int l, const int k, const double _sign){
  complex z,w;
  int ix;
  double sign = 1.; 
  spinor *r, *s;

  if(_sign < 0.){
    sign = -1.; 
  }

  z.re = 1.;
  z.im =  sign * g_mu;
  w.re = 1.;
  w.im = -z.im; /*-sign * nrm * g_mu;*/

  /************ loop over all lattice sites ************/
  for(ix = 0; ix < (VOLUME/2); ix++){
    s=&spinor_field[l][ix];
    r=&spinor_field[k][ix];

    /* Multiply the spinorfield with of 1+imu\gamma_5 */
    _complex_times_vector((*s).c1, z, (*r).c1);
    _complex_times_vector((*s).c2, z, (*r).c2);
    _complex_times_vector((*s).c3, w, (*r).c3);
    _complex_times_vector((*s).c4, w, (*r).c4);
  }
}

void mul_one_pm_imu_sub_mul_gamma5(const int l, const int k, 
				   const int j, const double _sign){
  complex z,w;
  int ix;
  double sign=1.;
  spinor *r, *s, *t;
  static su3_vector phi1, phi2, phi3, phi4;

  if(_sign < 0.){
    sign = -1.;
  }

  z.re = 1.;
  z.im =  sign * g_mu;
  w.re = 1.;
  w.im = -sign * g_mu;

  /************ loop over all lattice sites ************/
  for(ix = 0; ix < (VOLUME/2); ix++){
    r = &spinor_field[k][ix];
    s = &spinor_field[j][ix];
    t = &spinor_field[l][ix];
    /* Multiply the spinorfield with 1+imu\gamma_5 */
    _complex_times_vector(phi1, z, (*r).c1);
    _complex_times_vector(phi2, z, (*r).c2);
    _complex_times_vector(phi3, w, (*r).c3);
    _complex_times_vector(phi4, w, (*r).c4);
    /* Subtract s and store the result in t */
    /* multiply with  gamma5 included by    */
    /* reversed order of s and phi3|4       */
    _vector_sub((*t).c1, phi1, (*s).c1);
    _vector_sub((*t).c2, phi2, (*s).c2);
    _vector_sub((*t).c3, (*s).c3, phi3);
    _vector_sub((*t).c4, (*s).c4, phi4);
  }
}

void mul_one_pm_imu_sub_mul(const int l, const int k, 
			    const int j, const double _sign){
  complex z,w;
  int ix;
  double sign=1.;
  spinor *r, *s, *t;
  static su3_vector phi1, phi2, phi3, phi4;

  if(_sign < 0.){
    sign = -1.;
  }

  z.re = 1.;
  z.im =  sign * g_mu;
  w.re = 1.;
  w.im = -sign * g_mu;

  /************ loop over all lattice sites ************/
  for(ix = 0; ix < (VOLUME/2); ix++){
    r = &spinor_field[k][ix];
    s = &spinor_field[j][ix];
    t = &spinor_field[l][ix];
    /* Multiply the spinorfield with 1+imu\gamma_5 */
    _complex_times_vector(phi1, z, (*r).c1);
    _complex_times_vector(phi2, z, (*r).c2);
    _complex_times_vector(phi3, w, (*r).c3);
    _complex_times_vector(phi4, w, (*r).c4);
    /* Subtract s and store the result in t */
    _vector_sub((*t).c1, phi1, (*s).c1);
    _vector_sub((*t).c2, phi2, (*s).c2);
    _vector_sub((*t).c3, phi3, (*s).c3);
    _vector_sub((*t).c4, phi4, (*s).c4);
  }
}

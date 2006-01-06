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

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include "global.h"
#include "su3.h"
#include "Hopping_Matrix.h"
#include "Hopping_Matrix_nocom.h"
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
void mul_one_pm_imu_inv(spinor * const l, const double _sign);
void mul_one_pm_imu(spinor * const l, const double _sign);
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
void mul_one_pm_imu_sub_mul_gamma5(spinor * const l, spinor * const k, 
				   spinor * const j, const double _sign);

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
void mul_one_pm_imu_sub_mul(spinor * const l, spinor * const k,
			    spinor * const j, const double _sign);

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
void Qtm_plus_psi(spinor * const l, spinor * const k){
  Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX+1], k);
  mul_one_pm_imu_inv(g_spinor_field[DUM_MATRIX+1], +1.);
  Hopping_Matrix(OE, g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+1]);
  mul_one_pm_imu_sub_mul_gamma5(l, k, g_spinor_field[DUM_MATRIX], +1.);
}

void Qtm_plus_psi_nocom(spinor * const l, spinor * const k){
  Hopping_Matrix_nocom(EO, g_spinor_field[DUM_MATRIX+1], k);
  mul_one_pm_imu_inv(g_spinor_field[DUM_MATRIX+1], +1.);
  Hopping_Matrix_nocom(OE, g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+1]);
  mul_one_pm_imu_sub_mul_gamma5(l, k, g_spinor_field[DUM_MATRIX], +1.);
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
void Qtm_minus_psi(spinor * const l, spinor * const k){
  Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX+1], k);
  mul_one_pm_imu_inv(g_spinor_field[DUM_MATRIX+1], -1.);
  Hopping_Matrix(OE, g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+1]);
  mul_one_pm_imu_sub_mul_gamma5(l, k, g_spinor_field[DUM_MATRIX], -1.);
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
void Mtm_plus_psi(spinor * const l, spinor * const k){
  Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX+1], k);
  mul_one_pm_imu_inv(g_spinor_field[DUM_MATRIX+1], +1.);
  Hopping_Matrix(OE, g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+1]);
  mul_one_pm_imu_sub_mul(l, k, g_spinor_field[DUM_MATRIX], +1.);
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
void Mtm_minus_psi(spinor * const l, spinor * const k){
  Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX+1], k);
  mul_one_pm_imu_inv(g_spinor_field[DUM_MATRIX+1], -1.);
  Hopping_Matrix(OE, g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+1]);
  mul_one_pm_imu_sub_mul(l, k, g_spinor_field[DUM_MATRIX], -1.);
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
void Qtm_pm_psi(spinor * const l, spinor * const k){
  /* Q_{-} */
  Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX+1], k);
  mul_one_pm_imu_inv(g_spinor_field[DUM_MATRIX+1], -1.);
  Hopping_Matrix(OE, g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+1]);
  mul_one_pm_imu_sub_mul_gamma5(g_spinor_field[DUM_MATRIX], k, g_spinor_field[DUM_MATRIX], -1.);
  /* Q_{+} */
  Hopping_Matrix(EO, l, g_spinor_field[DUM_MATRIX]);
  mul_one_pm_imu_inv(l, +1.);
  Hopping_Matrix(OE, g_spinor_field[DUM_MATRIX+1], l);
  mul_one_pm_imu_sub_mul_gamma5(l, g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+1], +1.);
}

void Qtm_pm_psi_nocom(spinor * const l, spinor * const k){
  /* Q_{-} */
  Hopping_Matrix_nocom(EO, g_spinor_field[DUM_MATRIX+1], k);
  mul_one_pm_imu_inv(g_spinor_field[DUM_MATRIX+1], -1.);
  Hopping_Matrix_nocom(OE, g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+1]);
  mul_one_pm_imu_sub_mul_gamma5(g_spinor_field[DUM_MATRIX], k, g_spinor_field[DUM_MATRIX], -1.);
  /* Q_{+} */
  Hopping_Matrix_nocom(EO, l, g_spinor_field[DUM_MATRIX]);
  mul_one_pm_imu_inv(l, +1.);
  Hopping_Matrix_nocom(OE, g_spinor_field[DUM_MATRIX+1], l);
  mul_one_pm_imu_sub_mul_gamma5(l, g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+1], +1.);
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
void H_eo_tm_inv_psi(spinor * const l, spinor * const k, 
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

void mul_one_pm_imu_inv(spinor * const l, const double _sign){
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
    r=l + ix;
    /* Multiply the spinorfield with the inverse of 1+imu\gamma_5 */
    _complex_times_vector(phi1, z, (*r).s0);
    _vector_assign((*r).s0, phi1);
    _complex_times_vector(phi1, z, (*r).s1);
    _vector_assign((*r).s1, phi1);
    _complex_times_vector(phi1, w, (*r).s2);
    _vector_assign((*r).s2, phi1);
    _complex_times_vector(phi1, w, (*r).s3);
    _vector_assign((*r).s3, phi1);
  }
}

void assign_mul_one_pm_imu_inv(spinor * const l, spinor * const k, const double _sign){
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
    r=k+ix;
    s=l+ix;
    /* Multiply the spinorfield with the inverse of 1+imu\gamma_5 */
    _complex_times_vector((*s).s0, z, (*r).s0);
    _complex_times_vector((*s).s1, z, (*r).s1);
    _complex_times_vector((*s).s2, w, (*r).s2);
    _complex_times_vector((*s).s3, w, (*r).s3);
  }
}

void mul_one_pm_imu(spinor * const l, const double _sign){
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
    r=l+ix;
    /* Multiply the spinorfield with 1+imu\gamma_5 */
    _complex_times_vector(phi1, z, (*r).s0);
    _vector_assign((*r).s0, phi1);
    _complex_times_vector(phi1, z, (*r).s1);
    _vector_assign((*r).s1, phi1);
    _complex_times_vector(phi1, w, (*r).s2);
    _vector_assign((*r).s2, phi1);
    _complex_times_vector(phi1, w, (*r).s3);
    _vector_assign((*r).s3, phi1);
  }
}

void assign_mul_one_pm_imu(spinor * const l, spinor * const k, const double _sign){
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
    s=l+ix;
    r=k+ix;

    /* Multiply the spinorfield with of 1+imu\gamma_5 */
    _complex_times_vector((*s).s0, z, (*r).s0);
    _complex_times_vector((*s).s1, z, (*r).s1);
    _complex_times_vector((*s).s2, w, (*r).s2);
    _complex_times_vector((*s).s3, w, (*r).s3);
  }
}

void mul_one_pm_imu_sub_mul_gamma5(spinor * const l, spinor * const k, 
				   spinor * const j, const double _sign){
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
    r = k+ix;
    s = j+ix;
    t = l+ix;
    /* Multiply the spinorfield with 1+imu\gamma_5 */
    _complex_times_vector(phi1, z, (*r).s0);
    _complex_times_vector(phi2, z, (*r).s1);
    _complex_times_vector(phi3, w, (*r).s2);
    _complex_times_vector(phi4, w, (*r).s3);
    /* Subtract s and store the result in t */
    /* multiply with  gamma5 included by    */
    /* reversed order of s and phi3|4       */
    _vector_sub((*t).s0, phi1, (*s).s0);
    _vector_sub((*t).s1, phi2, (*s).s1);
    _vector_sub((*t).s2, (*s).s2, phi3);
    _vector_sub((*t).s3, (*s).s3, phi4);
  }
}

void mul_one_pm_imu_sub_mul(spinor * const l, spinor * const k, 
			    spinor * const j, const double _sign){
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
    r = k+ix;
    s = j+ix;
    t = l+ix;
    /* Multiply the spinorfield with 1+imu\gamma_5 */
    _complex_times_vector(phi1, z, (*r).s0);
    _complex_times_vector(phi2, z, (*r).s1);
    _complex_times_vector(phi3, w, (*r).s2);
    _complex_times_vector(phi4, w, (*r).s3);
    /* Subtract s and store the result in t */
    _vector_sub((*t).s0, phi1, (*s).s0);
    _vector_sub((*t).s1, phi2, (*s).s1);
    _vector_sub((*t).s2, phi3, (*s).s2);
    _vector_sub((*t).s3, phi4, (*s).s3);
  }
}

static char const rcsid[] = "$Id$";

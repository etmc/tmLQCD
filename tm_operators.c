/* $Id$ */
#include <stdlib.h>
#include <stdio.h>
#include "global.h"
#include "su3.h"
#include "H_eo.h"
#include "tm_operators.h"

/* internal */
void mul_one_pm_imu_inv(const int l, const double _sign);
void mul_one_pm_imu_sub_mul_gamma5(int l, int k, int j, const double _sign);

/* external functions */


void Qtm_plus_psi(const int l, const int k){
  H_eo(1, DUM_MATRIX+1, k);
  mul_one_pm_imu_inv(DUM_MATRIX+1, +1.);
  H_eo(0, DUM_MATRIX, DUM_MATRIX+1);
  mul_one_pm_imu_sub_mul_gamma5(l, k, DUM_MATRIX, +1.);
}

void Qtm_minus_psi(const int l, const int k){
  H_eo(1, DUM_MATRIX+1, k);
  mul_one_pm_imu_inv(DUM_MATRIX+1, -1.);
  H_eo(0, DUM_MATRIX, DUM_MATRIX+1);
  mul_one_pm_imu_sub_mul_gamma5(l, k, DUM_MATRIX, -1.);
}

void Qtm_pm_psi(const int l, const int k){
  /* Q_{-} */
  H_eo(1, DUM_MATRIX+1, k);
  mul_one_pm_imu_inv(DUM_MATRIX+1, -1.);
  H_eo(0, DUM_MATRIX, DUM_MATRIX+1);
  mul_one_pm_imu_sub_mul_gamma5(DUM_MATRIX, k, DUM_MATRIX, -1.);
  /* Q_{+} */
  H_eo(1, l, DUM_MATRIX);
  mul_one_pm_imu_inv(l, +1.);
  H_eo(0, DUM_MATRIX+1, l);
  mul_one_pm_imu_sub_mul_gamma5(l, DUM_MATRIX, DUM_MATRIX+1, +1.);
}

void H_eo_tm_inv_psi(const int l, const int k, const int ieo, const double sign){
  H_eo(ieo, l, k);
  mul_one_pm_imu_inv(l, sign);
}

/* !!
 * All the results are only stored in the first
 * half of the spinor fields, they have anly
 * lenght VOLUME/2
 * 
 * That's why mul_... do not need a iput
 * parameter ieo.
 *
 * !!
 */

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
  w.im = -sign * nrm * g_mu;

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

void mul_one_pm_imu_sub_mul_gamma5(int l, int k, int j, const double _sign){
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

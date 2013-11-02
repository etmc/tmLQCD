#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include "global.h"
#include "su3.h"
#include "operator/Hopping_Matrix.h"
#include "linalg_eo.h"
#include "gamma.h"

#include "tm_operators_32.h"



void mul_one_pm_imu_inv_32(spinor32 * const l, const float _sign, const int N){
  _Complex float ALIGN z,w;
  int ix;
  float sign=-1.; 
  spinor32 *r;

  su3_vector32 ALIGN phi1;

  double ALIGN nrm = 1./(1.+g_mu*g_mu);

  if(_sign < 0.){
    sign = 1.; 
  }

  z = nrm + (sign * nrm * g_mu) * I;
  w = conj(z);
  /************ loop over all lattice sites ************/
#ifdef OMP
#pragma omp for
#endif
  for(ix = 0; ix < N; ix++){
    r=l + ix;
    /* Multiply the spinorfield with the inverse of 1+imu\gamma_5 */
    _complex_times_vector(phi1, z, r->s0);
    _vector_assign(r->s0, phi1);
    _complex_times_vector(phi1, z, r->s1);
    _vector_assign(r->s1, phi1);
    _complex_times_vector(phi1, w, r->s2);
    _vector_assign(r->s2, phi1);
    _complex_times_vector(phi1, w, r->s3);
    _vector_assign(r->s3, phi1);
  }


}




void mul_one_pm_imu_sub_mul_gamma5_32(spinor32 * const l, spinor32 * const k, 
				   spinor32 * const j, const float _sign){

  _Complex float z,w;
  int ix;
  float sign=1.;
  spinor32 *r, *s, *t;

  su3_vector32 ALIGN phi1, phi2, phi3, phi4;

  if(_sign < 0.){
    sign = -1.;
  }

  z = 1. + (sign * g_mu) * I;
  w = conj(z);
  
  /************ loop over all lattice sites ************/
#ifdef OMP
#pragma omp for
#endif
  for(ix = 0; ix < (VOLUME/2); ix++){
    r = k+ix;
    s = j+ix;
    t = l+ix;
    /* Multiply the spinorfield with 1+imu\gamma_5 */
    _complex_times_vector(phi1, z, r->s0);
    _complex_times_vector(phi2, z, r->s1);
    _complex_times_vector(phi3, w, r->s2);
    _complex_times_vector(phi4, w, r->s3);
    /* Subtract s and store the result in t */
    /* multiply with  gamma5 included by    */
    /* reversed order of s and phi3|4       */
    _vector_sub(t->s0, phi1, s->s0);
    _vector_sub(t->s1, phi2, s->s1);
    _vector_sub(t->s2, s->s2, phi3);
    _vector_sub(t->s3, s->s3, phi4);
  }


}




void Qtm_pm_psi_32(spinor32 * const l, spinor32 * const k){
  /* Q_{-} */
#ifdef OMP
#pragma omp parallel
  {
#endif  
  Hopping_Matrix_32(EO, g_spinor_field32[1], k);
  mul_one_pm_imu_inv_32(g_spinor_field32[1], -1., VOLUME/2);
  Hopping_Matrix_32(OE, g_spinor_field32[0], g_spinor_field32[1]);
  mul_one_pm_imu_sub_mul_gamma5_32(g_spinor_field32[0], k, g_spinor_field32[0], -1.);
  /* Q_{+} */
  Hopping_Matrix_32(EO, l, g_spinor_field32[0]);
  mul_one_pm_imu_inv_32(l, +1., VOLUME/2);
  Hopping_Matrix_32(OE, g_spinor_field32[1], l);
  mul_one_pm_imu_sub_mul_gamma5_32(l, g_spinor_field32[0], g_spinor_field32[1], +1.);
#ifdef OMP
  } /* OpenMP closing brace */
#endif  
}








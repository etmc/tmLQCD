void _PSWITCH(mul_one_pm_imu_inv)(_PTSWITCH(spinor) * const l, const double _sign, const int N){
#ifdef TM_USE_OMP
#pragma omp parallel
  {
#endif
  _C_TYPE ALIGN z,w;
  int ix;
  double sign=-1.; 
  _PTSWITCH(spinor) *r;

  _PTSWITCH(su3_vector) ALIGN phi1;

  _F_TYPE ALIGN nrm = 1./(1.+g_mu*g_mu);

  if(_sign < 0.){
    sign = 1.; 
  }

  z = nrm + (sign * nrm * g_mu) * I;
  w = conj(z);
  /************ loop over all lattice sites ************/
#ifdef TM_USE_OMP
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

#ifdef TM_USE_OMP
  } /* OpenMP closing brace */
#endif

}

void _PSWITCH(assign_mul_one_pm_imu_inv)(_PTSWITCH(spinor) * const l, _PTSWITCH(spinor) * const k, 
					 const double _sign, const int N){
#ifdef TM_USE_OMP
#pragma omp parallel
  {
#endif
  _C_TYPE z,w;
  int ix;
  double sign=-1.; 
  _PTSWITCH(spinor) *r, *s;
  _F_TYPE nrm = (_F_TYPE) 1./(1.+g_mu*g_mu);

  if(_sign < 0.){
    sign = 1.; 
  }

  z = nrm + (_F_TYPE) (sign * nrm * g_mu) * I;
  w = conj(z);

  /************ loop over all lattice sites ************/
#ifdef TM_USE_OMP
#pragma omp for
#endif
  for(ix = 0; ix < N; ix++){
    r=k+ix;
    s=l+ix;
    /* Multiply the spinorfield with the inverse of 1+imu\gamma_5 */
    _complex_times_vector(s->s0, z, r->s0);
    _complex_times_vector(s->s1, z, r->s1);
    _complex_times_vector(s->s2, w, r->s2);
    _complex_times_vector(s->s3, w, r->s3);
  }

#ifdef TM_USE_OMP
  } /* OpenMP closing brace */
#endif
}

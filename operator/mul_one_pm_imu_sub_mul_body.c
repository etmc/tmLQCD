void _PSWITCH(mul_one_pm_imu_sub_mul)(_PTSWITCH(spinor) * const l, _PTSWITCH(spinor) * const k, 
				      _PTSWITCH(spinor) * const j, const double _sign, const int N){
#ifdef TM_USE_OMP
#pragma omp parallel
  {
#endif
  _C_TYPE z,w;
  int ix;
  double sign=1.;
  _PTSWITCH(spinor) *r, *s, *t;

#if (!defined SSE2 && !defined SSE3)

  _PTSWITCH(su3_vector) ALIGN phi1, phi2, phi3, phi4;
  
#endif

  if(_sign < 0.){
    sign = -1.;
  }

  z = 1. + (sign * g_mu) * I;
  w = conj(z);
  /************ loop over all lattice sites ************/
#ifdef TM_USE_OMP
#pragma omp for
#endif
  for(ix = 0; ix < N; ix++){
    r = k+ix;
    s = j+ix;
    t = l+ix;
    /* Multiply the spinorfield with 1+imu\gamma_5 */
    _complex_times_vector(phi1, z, r->s0);
    _complex_times_vector(phi2, z, r->s1);
    _complex_times_vector(phi3, w, r->s2);
    _complex_times_vector(phi4, w, r->s3);
    /* Subtract s and store the result in t */
    _vector_sub(t->s0, phi1, s->s0);
    _vector_sub(t->s1, phi2, s->s1);
    _vector_sub(t->s2, phi3, s->s2);
    _vector_sub(t->s3, phi4, s->s3);
  }

#ifdef TM_USE_OMP
  } /* OpenMP closing brace */
#endif
}


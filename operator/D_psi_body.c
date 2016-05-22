static inline void _PTSWITCH(p0add)(_PTSWITCH(spinor) * restrict const tmpr , _PTSWITCH(spinor) const * restrict const s, 
				    _PSWITCH(su3) const * restrict const u, const _C_TYPE phase) {
#ifdef TM_USE_OMP
#define static
#endif
  static _PTSWITCH(su3_vector) chi, psi;
#ifdef TM_USE_OMP
#undef static
#endif

  _vector_add(psi,s->s0, s->s2);
  _su3_multiply(chi, (*u), psi);

  _complex_times_vector(psi, phase, chi);
  _vector_add_assign(tmpr->s0, psi);
  _vector_add_assign(tmpr->s2, psi);

  _vector_add(psi, s->s1, s->s3);
  _su3_multiply(chi, (*u), psi);

  _complex_times_vector(psi, phase, chi);
  _vector_add_assign(tmpr->s1, psi);
  _vector_add_assign(tmpr->s3, psi);

  return;
}


static inline void _PTSWITCH(m0add)(_PTSWITCH(spinor) * restrict const tmpr, _PTSWITCH(spinor) const * restrict const s, 
				    _PSWITCH(su3) const * restrict const u, const _C_TYPE phase) {
#ifdef TM_USE_OMP
#define static
#endif
  static _PTSWITCH(su3_vector) chi, psi;
#ifdef TM_USE_OMP
#undef static
#endif

  _vector_sub(psi, s->s0, s->s2);
  _su3_inverse_multiply(chi, (*u), psi);

  _complexcjg_times_vector(psi, phase, chi);
  _vector_add_assign(tmpr->s0, psi);
  _vector_sub_assign(tmpr->s2, psi);

  _vector_sub(psi, s->s1, s->s3);
  _su3_inverse_multiply(chi, (*u), psi);

  _complexcjg_times_vector(psi, phase, chi);
  _vector_add_assign(tmpr->s1, psi);
  _vector_sub_assign(tmpr->s3, psi);

  return;
}

static inline void _PTSWITCH(p1add)(_PTSWITCH(spinor) * restrict const tmpr, _PTSWITCH(spinor) const * restrict const s, 
				    _PSWITCH(su3) const * restrict const u, const _C_TYPE phase) {
#ifdef TM_USE_OMP
#define static
#endif
  static _PTSWITCH(su3_vector) chi, psi;
#ifdef TM_USE_OMP
#undef static
#endif

  _vector_i_add(psi,s->s0,s->s3);
  _su3_multiply(chi,(*u),psi);

  _complex_times_vector(psi, phase, chi);
  _vector_add_assign(tmpr->s0, psi);
  _vector_i_sub_assign(tmpr->s3, psi);
 
  _vector_i_add(psi, s->s1, s->s2);
  _su3_multiply(chi, (*u), psi);

  _complex_times_vector(psi, phase, chi);
  _vector_add_assign(tmpr->s1, psi);
  _vector_i_sub_assign(tmpr->s2, psi);

  return;
}

static inline void _PTSWITCH(m1add)(_PTSWITCH(spinor) * restrict const tmpr, _PTSWITCH(spinor) const * restrict const s, 
				    _PSWITCH(su3) const * restrict const u, const _C_TYPE phase) {
#ifdef TM_USE_OMP
#define static
#endif
  static _PTSWITCH(su3_vector) chi, psi;
#ifdef TM_USE_OMP
#undef static
#endif

  _vector_i_sub(psi,s->s0, s->s3);
  _su3_inverse_multiply(chi,(*u), psi);

  _complexcjg_times_vector(psi, phase, chi);
  _vector_add_assign(tmpr->s0, psi);
  _vector_i_add_assign(tmpr->s3, psi);

  _vector_i_sub(psi, s->s1, s->s2);
  _su3_inverse_multiply(chi, (*u), psi);

  _complexcjg_times_vector(psi, phase, chi);
  _vector_add_assign(tmpr->s1, psi);
  _vector_i_add_assign(tmpr->s2, psi);

  return;
}

static inline void _PTSWITCH(p2add)(_PTSWITCH(spinor) * restrict const tmpr, _PTSWITCH(spinor) const * restrict const s, 
				    _PSWITCH(su3) const * restrict const u, const _C_TYPE phase) {
#ifdef TM_USE_OMP
#define static
#endif
  static _PTSWITCH(su3_vector) chi, psi;
#ifdef TM_USE_OMP
#undef static
#endif

  _vector_add(psi,s->s0,s->s3);
  _su3_multiply(chi, (*u), psi);

  _complex_times_vector(psi, phase, chi);
  _vector_add_assign(tmpr->s0, psi);
  _vector_add_assign(tmpr->s3, psi);

  _vector_sub(psi,s->s1,s->s2);
  _su3_multiply(chi, (*u), psi);

  _complex_times_vector(psi, phase, chi);
  _vector_add_assign(tmpr->s1, psi);
  _vector_sub_assign(tmpr->s2, psi);


  return;
}

static inline void _PTSWITCH(m2add)(_PTSWITCH(spinor) * restrict const tmpr, _PTSWITCH(spinor) const * restrict const s, 
				    _PSWITCH(su3) const * restrict const u, const _C_TYPE phase) {
#ifdef TM_USE_OMP
#define static
#endif
  static _PTSWITCH(su3_vector) chi, psi;
#ifdef TM_USE_OMP
#undef static
#endif

  _vector_sub(psi, s->s0, s->s3);
  _su3_inverse_multiply(chi, (*u), psi);

  _complexcjg_times_vector(psi, phase, chi);
  _vector_add_assign(tmpr->s0, psi);
  _vector_sub_assign(tmpr->s3, psi);

  _vector_add(psi, s->s1, s->s2);
  _su3_inverse_multiply(chi, (*u),psi);

  _complexcjg_times_vector(psi, phase, chi);
  _vector_add_assign(tmpr->s1, psi);
  _vector_add_assign(tmpr->s2, psi);

  return;
}

static inline void _PTSWITCH(p3add)(_PTSWITCH(spinor) * restrict const tmpr, _PTSWITCH(spinor) const * restrict const s, 
				    _PSWITCH(su3) const * restrict const u, const _C_TYPE phase) {
#ifdef TM_USE_OMP
#define static
#endif
  static _PTSWITCH(su3_vector) chi, psi;
#ifdef TM_USE_OMP
#undef static
#endif

  _vector_i_add(psi, s->s0, s->s2);
  _su3_multiply(chi, (*u), psi);

  _complex_times_vector(psi, phase, chi);
  _vector_add_assign(tmpr->s0, psi);
  _vector_i_sub_assign(tmpr->s2, psi);

  _vector_i_sub(psi,s->s1, s->s3);
  _su3_multiply(chi, (*u), psi);

  _complex_times_vector(psi, phase, chi);
  _vector_add_assign(tmpr->s1, psi);
  _vector_i_add_assign(tmpr->s3, psi);

  return;
}

static inline void _PTSWITCH(m3addandstore)(_PTSWITCH(spinor) * restrict const r, _PTSWITCH(spinor) const * restrict const s, 
					    _PSWITCH(su3) const * restrict const u, const _C_TYPE phase,
					    _PTSWITCH(spinor) const * restrict const tmpr) {
#ifdef TM_USE_OMP
#define static
#endif
  static _PTSWITCH(su3_vector) chi, psi;
#ifdef TM_USE_OMP
#undef static
#endif

  _vector_i_sub(psi,s->s0, s->s2);
  _su3_inverse_multiply(chi, (*u), psi);

  _complexcjg_times_vector(psi, phase, chi);
  _vector_add(r->s0, tmpr->s0, psi);
  _vector_i_add(r->s2, tmpr->s2, psi);

  _vector_i_add(psi, s->s1, s->s3);
  _su3_inverse_multiply(chi, (*u), psi);

  _complexcjg_times_vector(psi, phase, chi);
  _vector_add(r->s1, tmpr->s1, psi);
  _vector_i_sub(r->s3, tmpr->s3, psi);

  return;
}

/* this is the hopping part only */
static inline void _PSWITCH(local_H)(_PTSWITCH(spinor) * const rr, _PTSWITCH(spinor) const * const s, 
				      _PSWITCH(su3) const * restrict u, 
				      int * _idx, _PTSWITCH(spinor) * const restrict tmpr) {
  // convert phases to _C_TYPE locally
  _C_TYPE ALIGN32 phase_0l = (_C_TYPE) phase_0;
  _C_TYPE ALIGN32 phase_1l = (_C_TYPE) phase_1;
  _C_TYPE ALIGN32 phase_2l = (_C_TYPE) phase_2;
  _C_TYPE ALIGN32 phase_3l = (_C_TYPE) phase_3;  

  int * idx = _idx;

  /****** direction +0 ******/
  _PTSWITCH(p0add)(tmpr, s + (*idx), u, phase_0l);
  u++;
  idx++;
  /****** direction -0 ******/
  _PTSWITCH(m0add)(tmpr, s + (*idx), u, phase_0l);
  u++;
  idx++;
  /****** direction +1 ******/
  _PTSWITCH(p1add)(tmpr, s + (*idx), u, phase_1l);
  u++;
  idx++;
  /****** direction -1 ******/
  _PTSWITCH(m1add)(tmpr, s + (*idx), u, phase_1l);
  u++;
  idx++;
  /****** direction +2 ******/
  _PTSWITCH(p2add)(tmpr, s + (*idx), u, phase_2l);
  u++;
  idx++;
  /****** direction -2 ******/
  _PTSWITCH(m2add)(tmpr, s + (*idx), u, phase_2l);
  u++;
  idx++;
  /****** direction +3 ******/
  _PTSWITCH(p3add)(tmpr, s + (*idx), u, phase_3l);
  u++;
  idx++;
  /****** direction -3 ******/
  _PTSWITCH(m3addandstore)(rr, s + (*idx), u, phase_3l, tmpr);

  return;
}

void _PSWITCH(D_psi)(_PTSWITCH(spinor) * const P, _PTSWITCH(spinor) * const Q){
  if(P==Q){
    printf("Error in D_psi (operator.c):\n");
    printf("Arguments must be different spinor fields\n");
    printf("Program aborted\n");
    exit(1);
  }
  //convert phases to float locally
  _C_TYPE ALIGN32 phase_0l = (_C_TYPE) phase_0;
  _C_TYPE ALIGN32 phase_1l = (_C_TYPE) phase_1;
  _C_TYPE ALIGN32 phase_2l = (_C_TYPE) phase_2;
  _C_TYPE ALIGN32 phase_3l = (_C_TYPE) phase_3;  

#ifdef _GAUGE_COPY
  if(_PSWITCH(g_update_gauge_copy)) {
    _PSWITCH(update_backward_gauge)(_PSWITCH(g_gauge_field));
  }
#endif
# if defined MPI
  _PTSWITCH(xchange_lexicfield)(Q);
# endif

#ifdef TM_USE_OMP
#pragma omp parallel
  {
#endif

  int ix,iy;
  _PSWITCH(su3) * restrict up,* restrict um;
  _PTSWITCH(spinor) * restrict rr; 
  _PTSWITCH(spinor) const * restrict s;
  _PTSWITCH(spinor) const * restrict sp;
  _PTSWITCH(spinor) const * restrict sm;
  _C_TYPE rho1, rho2;
  _PTSWITCH(spinor) tmpr;

  rho1 = (_F_TYPE)1. + (_F_TYPE) g_mu * I;
  rho2 = conj(rho1);

  /************************ loop over all lattice sites *************************/

#ifdef TM_USE_OMP
#pragma omp for
#endif
  for (ix = 0; ix < VOLUME; ix++) {
    rr  = (_PTSWITCH(spinor) *) P +ix;
    s  = (_PTSWITCH(spinor) *) Q +ix;

    if(g_c_sw > 0) {
      _PSWITCH(assign_mul_one_sw_pm_imu_site_lexic)(ix, &tmpr, s, (_F_TYPE) g_mu);
    }
    else {
      _complex_times_vector(tmpr.s0, rho1, s->s0);
      _complex_times_vector(tmpr.s1, rho1, s->s1);
      _complex_times_vector(tmpr.s2, rho2, s->s2);
      _complex_times_vector(tmpr.s3, rho2, s->s3);
    }

    /******************************* direction +0 *********************************/
    iy=g_iup[ix][0];
    sp = (_PTSWITCH(spinor) *) Q +iy;
    up=&_PSWITCH(g_gauge_field)[ix][0];
    _PTSWITCH(p0add)(&tmpr, sp, up, phase_0l);

    /******************************* direction -0 *********************************/
    iy=g_idn[ix][0];
    sm  = (_PTSWITCH(spinor) *) Q +iy;
    um=&_PSWITCH(g_gauge_field)[iy][0];
    _PTSWITCH(m0add)(&tmpr, sm, um, phase_0l);

    /******************************* direction +1 *********************************/
    iy=g_iup[ix][1];
    sp = (_PTSWITCH(spinor) *) Q +iy;
    up=&_PSWITCH(g_gauge_field)[ix][1];
    _PTSWITCH(p1add)(&tmpr, sp, up, phase_1l);

    /******************************* direction -1 *********************************/
    iy=g_idn[ix][1];
    sm = (_PTSWITCH(spinor) *) Q +iy;
    um=&_PSWITCH(g_gauge_field)[iy][1];
    _PTSWITCH(m1add)(&tmpr, sm, um, phase_1l);

    /******************************* direction +2 *********************************/
    iy=g_iup[ix][2];
    sp = (_PTSWITCH(spinor) *) Q +iy;
    up=&_PSWITCH(g_gauge_field)[ix][2];
    _PTSWITCH(p2add)(&tmpr, sp, up, phase_2l);

    /******************************* direction -2 *********************************/
    iy=g_idn[ix][2];
    sm = (_PTSWITCH(spinor) *) Q +iy;
    um=&_PSWITCH(g_gauge_field)[iy][2];
    _PTSWITCH(m2add)(&tmpr, sm, um, phase_2l);

    /******************************* direction +3 *********************************/
    iy=g_iup[ix][3];
    sp = (_PTSWITCH(spinor) *) Q +iy;
    up=&_PSWITCH(g_gauge_field)[ix][3];
    _PTSWITCH(p3add)(&tmpr, sp, up, phase_3l);

    /******************************* direction -3 *********************************/
    iy=g_idn[ix][3];
    sm = (_PTSWITCH(spinor) *) Q +iy;
    um=&_PSWITCH(g_gauge_field)[iy][3];
    _PTSWITCH(m3addandstore)(rr, sm, um, phase_3l, &tmpr);
  }
#ifdef TM_USE_OMP
  } /* OpenMP closing brace */
#endif
}

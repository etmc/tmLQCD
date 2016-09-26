void _PSWITCH(assign_mul_one_sw_pm_imu_site_lexic)(const int ix,
						    _PTSWITCH(spinor) * const k, const _PTSWITCH(spinor) * const l, 
						    const _F_TYPE mu) {

  _PTSWITCH(su3_vector) ALIGN chi, psi1, psi2;
  const _PSWITCH(su3) *w1, *w2, *w3;

  // upper two spin components first
  w1 = &_PSWITCH(sw)[ix][0][0];
  w2 = w1 + 2; /*&sw[ix][1][0];*/
  w3 = w1 + 4; /*&sw[ix][2][0];*/
  _su3_multiply(psi1, *w1, (*l).s0);
  _su3_multiply(chi, *w2, (*l).s1);
  _vector_add_assign(psi1, chi);
  _su3_inverse_multiply(psi2, *w2, (*l).s0);
  _su3_multiply(chi, *w3, (*l).s1);
  _vector_add_assign(psi2, chi);

  // add in the twisted mass term (plus in the upper components)
  _vector_add_i_mul(psi1, mu, (*l).s0);
  _vector_add_i_mul(psi2, mu, (*l).s1);

  _vector_assign((*k).s0, psi1);
  _vector_assign((*k).s1, psi2);

  // now lower to spin components
  w1++; /*=&sw[ix][0][1];*/
  w2++; /*=&sw[ix][1][1];*/
  w3++; /*=&sw[ix][2][1];*/
  _su3_multiply(psi1, *w1, (*l).s2);
  _su3_multiply(chi, *w2, (*l).s3);
  _vector_add_assign(psi1, chi);
  _su3_inverse_multiply(psi2, *w2, (*l).s2);
  _su3_multiply(chi, *w3, (*l).s3);
  _vector_add_assign(psi2, chi);

  // add in the twisted mass term (minus from g5 in the lower components)
  _vector_add_i_mul(psi1, -mu, (*l).s2);
  _vector_add_i_mul(psi2, -mu, (*l).s3);

  _vector_assign((*k).s2, psi1);
  _vector_assign((*k).s3, psi2);
  return;
}

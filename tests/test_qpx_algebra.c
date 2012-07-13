#include <stdio.h>
#include <config.h>
#include <complex.h>
#include <cu/cu.h>
#if (defined SSE || defined SSE2 || defined SSE3)
# include "sse.h"
#endif
#include "../su3.h"
#include "../su3adj.h"
#include "../expo.h"
#if (!defined BGQ && defined XLC) 
# include "../bgq.h"
#endif

#define EPS 5e-16

TEST(qpx_algebra) {
  int test = 0;
#if (!defined BGQ && defined XLC)
  complex double ca, cb, cc, cd, k ALIGN;
  vector4double a, b, c, d;
  su3 u ALIGN;
  su3_vector phi0, phi1, phi2, phi3, phi4, phi5 ALIGN;
  spinor s, temp, temp2 ALIGN;
  spinor * sp;
  vector4double rs[12];
  vector4double r[12];
  vector4double U[9];

  k = 0.3-0.7*I;
  phi0.c0 = 0.9+2*I;
  phi0.c1 = 0.6+1.3*I;
  phi0.c2 = 0.6-.2*I;

  phi1.c0 = -0.1+1*I;
  phi1.c1 = 0.1+0.3*I;
  phi1.c2 = -0.1-.12*I;

  u.c00 = +0.3391 -0.1635*I;
  u.c01 = -0.2357 +0.5203*I;
  u.c02 = +0.5609 +0.4663*I;
  u.c10 = -0.0740 -0.4204*I;
  u.c11 = -0.7706 -0.1863*I;
  u.c12 = +0.1191 -0.4185*I;
  u.c20 = +0.5351 -0.6243*I;
  u.c21 = +0.1825 +0.1089*I;
  u.c22 = -0.5279 -0.0022*I;

  _su3_multiply(phi2, u, phi0);
  _su3_multiply(phi3, u, phi1);

  vec_load2(r, &phi0);
  vec_load2(&r[3], &phi1);
  vec_su3_multiply_double2(&u, U, r);
  vec_store2(&phi4, &r[6]);
  vec_store2(&phi5, &r[9]);
  if( cabs(phi4.c0 - phi2.c0) > EPS || 
      cabs(phi4.c1 - phi2.c1) > EPS || 
      cabs(phi4.c2 - phi2.c2) > EPS || 
      cabs(phi5.c0 - phi3.c0) > EPS || 
      cabs(phi5.c1 - phi3.c1) > EPS || 
      cabs(phi5.c2 - phi3.c2) > EPS )
    test = 1;

  assertFalseM(test, "vec_su3_multiply_double2 failed\n");
  test = 0;

  _su3_inverse_multiply(phi2, u, phi0);
  _su3_inverse_multiply(phi3, u, phi1);

  vec_load2(r, &phi0);
  vec_load2(&r[3], &phi1);
  vec_su3_inverse_multiply_double2(&u, U, r);
  vec_store2(&phi4, &r[6]);
  vec_store2(&phi5, &r[9]);
  if( cabs(phi4.c0 - phi2.c0) > EPS || 
      cabs(phi4.c1 - phi2.c1) > EPS || 
      cabs(phi4.c2 - phi2.c2) > EPS || 
      cabs(phi5.c0 - phi3.c0) > EPS || 
      cabs(phi5.c1 - phi3.c1) > EPS || 
      cabs(phi5.c2 - phi3.c2) > EPS )
    test = 1;

  assertFalseM(test, "vec_su3_inverse_multiply_double2 failed\n");
  test = 0;

  vec_load2(r, &phi0);
  vec_load2(r+3, &phi1);
  vec_add2(r, r+3);
  vec_store2(&phi4, r);

  if( cabs(phi4.c0 - phi0.c0 - phi1.c0) > EPS || 
      cabs(phi4.c1 - phi0.c1 - phi1.c1) > EPS || 
      cabs(phi4.c2 - phi0.c2 - phi1.c2) > EPS )
    test = 1;

  assertFalseM(test, "vec_add2 failed\n");
  test = 0;

  vec_load2(r, &phi0);
  vec_load2(r+3, &phi1);
  vec_sub2(r, r+3);
  vec_store2(&phi4, r);

  if( cabs(phi4.c0 - phi0.c0 + phi1.c0) > EPS || 
      cabs(phi4.c1 - phi0.c1 + phi1.c1) > EPS || 
      cabs(phi4.c2 - phi0.c2 + phi1.c2) > EPS )
    test = 1;

  assertFalseM(test, "vec_sub2 failed\n");
  test = 0;

  vec_load2(r, &phi0);
  vec_load2(r+3, &phi1);
  vec_cmplx_mul_double2(&r[6], r, U, &k);
  vec_store2(&phi4, r+6);
  vec_store2(&phi5, r+9);

  if( cabs(phi4.c0 - k*phi0.c0) > EPS || 
      cabs(phi4.c1 - k*phi0.c1) > EPS || 
      cabs(phi4.c2 - k*phi0.c2) > EPS ||
      cabs(phi5.c0 - k*phi1.c0) > EPS || 
      cabs(phi5.c1 - k*phi1.c1) > EPS || 
      cabs(phi5.c2 - k*phi1.c2) > EPS )
    test = 1;

  assertFalseM(test, "vec_cmplx_mul_double2 failed\n");
  test = 0;

  vec_load2(r, &phi0);
  vec_load2(r+3, &phi1);
  vec_cmplxcg_mul_double2(r+6, r, U, &k);
  vec_store2(&phi4, r+6);
  vec_store2(&phi5, r+9);

  if( cabs(phi4.c0 - conj(k)*phi0.c0) > EPS || 
      cabs(phi4.c1 - conj(k)*phi0.c1) > EPS || 
      cabs(phi4.c2 - conj(k)*phi0.c2) > EPS ||
      cabs(phi5.c0 - conj(k)*phi1.c0) > EPS || 
      cabs(phi5.c1 - conj(k)*phi1.c1) > EPS || 
      cabs(phi5.c2 - conj(k)*phi1.c2) > EPS )
    test = 1;

  assertFalseM(test, "vec_cmplxcg_mul_double2 failed\n");
  test = 0;

  vec_load2(r, &phi0);
  vec_load2(r+3, &phi1);
  vec_i_mul_add2(r, r+3, U);
  vec_store2(&phi4, r);

  if( cabs(phi4.c0 - phi0.c0 - I*phi1.c0) > EPS || 
      cabs(phi4.c1 - phi0.c1 - I*phi1.c1) > EPS || 
      cabs(phi4.c2 - phi0.c2 - I*phi1.c2) > EPS )
    test = 1;

  assertFalseM(test, "vec_i_mul_add2 failed\n");
  test = 0;

  vec_load2(r, &phi0);
  vec_load2(r+3, &phi1);
  vec_i_mul_sub2(r, r+3, U);
  vec_store2(&phi4, r);

  if( cabs(phi4.c0 - phi0.c0 + I*phi1.c0) > EPS || 
      cabs(phi4.c1 - phi0.c1 + I*phi1.c1) > EPS || 
      cabs(phi4.c2 - phi0.c2 + I*phi1.c2) > EPS )
    test = 1;

  assertFalseM(test, "vec_i_mul_sub2 failed\n");
  test = 0;

  phi2.c0 = 0.9+2*I;
  phi2.c1 = 0.6+1.3*I;
  phi2.c2 = 0.6-.2*I;

  phi3.c0 = -0.4+0.3*I;
  phi3.c1 = 0.5-1.3*I;
  phi3.c2 = -0.1-3.12*I;
  
  s.s0 = phi0;
  s.s1 = phi1;
  s.s2 = phi2;
  s.s3 = phi3;
  sp = &s;
  _vector_add(psi,sp->s0,sp->s2);
  _su3_multiply(chi,u,psi);
  _complex_times_vector(psi,k,chi);
  _vector_assign(temp.s0,psi);
  _vector_assign(temp.s2,psi);
  _vector_add(psi, sp->s1, sp->s3);
  _su3_multiply(chi,u,psi);
  _complex_times_vector(psi,k,chi);
  _vector_assign(temp.s1,psi);
  _vector_assign(temp.s3,psi);

  vec_load2(r, &sp->s0);
  vec_load2(r+3, &sp->s1);
  vec_load2(r+6, &sp->s2);
  vec_load2(r+9, &sp->s3);
  // s0 + s2 and s1 + s3
  vec_add_double2(r, &r[6]);
  // result is now in r[0-5] 
  vec_su3_multiply_double2(&u, U, r);
  // result is now in r[6-11]
  // mult with ka0 and store in rs
  vec_cmplx_mul_double2(rs, &r[6], U, &k);
  rs[6] = rs[0]; rs[7] = rs[1]; rs[8] = rs[2];
  rs[9] = rs[3]; rs[10]= rs[4]; rs[11]= rs[5];
  vec_store2(&temp2.s0, rs);
  vec_store2(&temp2.s1, rs+3);
  vec_store2(&temp2.s2, rs+6);
  vec_store2(&temp2.s3, rs+9);

  if( cabs(temp.s0.c0 - temp2.s0.c0) > EPS || 
      cabs(temp.s0.c1 - temp2.s0.c1) > EPS || 
      cabs(temp.s0.c2 - temp2.s0.c2) > EPS ||
      cabs(temp.s1.c0 - temp2.s1.c0) > EPS || 
      cabs(temp.s1.c1 - temp2.s1.c1) > EPS || 
      cabs(temp.s1.c2 - temp2.s1.c2) > EPS ||
      cabs(temp.s2.c0 - temp2.s2.c0) > EPS || 
      cabs(temp.s2.c1 - temp2.s2.c1) > EPS || 
      cabs(temp.s2.c2 - temp2.s2.c2) > EPS ||
      cabs(temp.s3.c0 - temp2.s3.c0) > EPS || 
      cabs(temp.s3.c1 - temp2.s3.c1) > EPS || 
      cabs(temp.s3.c2 - temp2.s3.c2) > EPS )
    test = 1;

  assertFalseM(test, "D_t+ failed\n");
  test = 0;

  sm = sp;
  _vector_sub(psi,sm->s0,sm->s2);
  _su3_inverse_multiply(chi,u,psi);
  _complexcjg_times_vector(psi,k,chi);
  _vector_add_assign(temp.s0,psi);
  _vector_sub_assign(temp.s2,psi);
  _vector_sub(psi,sm->s1,sm->s3);
  _su3_inverse_multiply(chi,u,psi);
  _complexcjg_times_vector(psi,k,chi);
  _vector_add_assign(temp.s1,psi);
  _vector_sub_assign(temp.s3,psi);

  vec_load2(r, &sm->s0);
  vec_load2(r+3, &sm->s1);
  vec_load2(r+6, &sm->s2);
  vec_load2(r+9, &sm->s3);
  // s0 - s2 and s1 - s3
  vec_sub_double2(r, &r[6]);
  // result is now in r[0-5]
  vec_su3_inverse_multiply_double2(&u, U, r);
  // result is now in r[6-11]
  // mult with k0
  vec_cmplxcg_mul_double2(r, &r[6], U, &k);
  // result in r[0-5] now
  vec_add_double2(rs, r);
  vec_sub_double2(&rs[6], r);

  vec_store2(&temp2.s0, rs);
  vec_store2(&temp2.s1, rs+3);
  vec_store2(&temp2.s2, rs+6);
  vec_store2(&temp2.s3, rs+9);

  if( cabs(temp.s0.c0 - temp2.s0.c0) > EPS || 
      cabs(temp.s0.c1 - temp2.s0.c1) > EPS || 
      cabs(temp.s0.c2 - temp2.s0.c2) > EPS ||
      cabs(temp.s1.c0 - temp2.s1.c0) > EPS || 
      cabs(temp.s1.c1 - temp2.s1.c1) > EPS || 
      cabs(temp.s1.c2 - temp2.s1.c2) > EPS ||
      cabs(temp.s2.c0 - temp2.s2.c0) > EPS || 
      cabs(temp.s2.c1 - temp2.s2.c1) > EPS || 
      cabs(temp.s2.c2 - temp2.s2.c2) > EPS ||
      cabs(temp.s3.c0 - temp2.s3.c0) > EPS || 
      cabs(temp.s3.c1 - temp2.s3.c1) > EPS || 
      cabs(temp.s3.c2 - temp2.s3.c2) > EPS )
    test = 1;

  assertFalseM(test, "D_t- failed\n");
  test = 0;

  _vector_i_add(psi,sp->s0,sp->s3);
  _su3_multiply(chi,u,psi);
  _complex_times_vector(psi,k,chi);
  _vector_add_assign(temp.s0,psi);
  _vector_i_sub_assign(temp.s3,psi);
  _vector_i_add(psi,sp->s1,sp->s2);
  _su3_multiply(chi,u,psi);
  _complex_times_vector(psi,k,chi);
  _vector_add_assign(temp.s1,psi);
  _vector_i_sub_assign(temp.s2,psi);

  vec_load2(r, &sp->s0);
  vec_load2(r+3, &sp->s1);
  vec_load2(r+9, &sp->s2);
  vec_load2(r+6, &sp->s3);
  vec_i_mul_add_double2(r, &r[6], U);
  vec_su3_multiply_double2(&u, U, r);
  vec_cmplx_mul_double2(r, &r[6], U, &k);
  vec_add_double2(rs, r);
  vec_i_mul_sub2(&rs[6], &r[3], U);
  vec_i_mul_sub2(&rs[9], &r[0], U);

  vec_store2(&temp2.s0, rs);
  vec_store2(&temp2.s1, rs+3);
  vec_store2(&temp2.s2, rs+6);
  vec_store2(&temp2.s3, rs+9);

  if( cabs(temp.s0.c0 - temp2.s0.c0) > EPS || 
      cabs(temp.s0.c1 - temp2.s0.c1) > EPS || 
      cabs(temp.s0.c2 - temp2.s0.c2) > EPS ||
      cabs(temp.s1.c0 - temp2.s1.c0) > EPS || 
      cabs(temp.s1.c1 - temp2.s1.c1) > EPS || 
      cabs(temp.s1.c2 - temp2.s1.c2) > EPS ||
      cabs(temp.s2.c0 - temp2.s2.c0) > EPS || 
      cabs(temp.s2.c1 - temp2.s2.c1) > EPS || 
      cabs(temp.s2.c2 - temp2.s2.c2) > EPS ||
      cabs(temp.s3.c0 - temp2.s3.c0) > EPS || 
      cabs(temp.s3.c1 - temp2.s3.c1) > EPS || 
      cabs(temp.s3.c2 - temp2.s3.c2) > EPS )
    test = 1;

  assertFalseM(test, "D_x+ failed\n");
  test = 0;

  _vector_i_sub(psi,sm->s0,sm->s3);
  _su3_inverse_multiply(chi,u,psi);
  _complexcjg_times_vector(psi,k,chi);
  _vector_add_assign(temp.s0,psi);
  _vector_i_add_assign(temp.s3,psi);
  _vector_i_sub(psi,sm->s1,sm->s2);
  _su3_inverse_multiply(chi,u,psi);
  _complexcjg_times_vector(psi,k,chi);
  _vector_add_assign(temp.s1,psi);
  _vector_i_add_assign(temp.s2,psi);

  vec_load2(r, &sm->s0);
  vec_load2(r+3, &sm->s1);
  vec_load2(r+9, &sm->s2);
  vec_load2(r+6, &sm->s3);
  vec_i_mul_sub_double2(r, &r[6], U);
  vec_su3_inverse_multiply_double2(&u, U, r);
  vec_cmplxcg_mul_double2(r, &r[6], U, &k);
  vec_add_double2(rs, r);
  vec_i_mul_add2(&rs[6], &r[3], U);
  vec_i_mul_add2(&rs[9], &r[0], U);

  vec_store2(&temp2.s0, rs);
  vec_store2(&temp2.s1, rs+3);
  vec_store2(&temp2.s2, rs+6);
  vec_store2(&temp2.s3, rs+9);

  if( cabs(temp.s0.c0 - temp2.s0.c0) > EPS || 
      cabs(temp.s0.c1 - temp2.s0.c1) > EPS || 
      cabs(temp.s0.c2 - temp2.s0.c2) > EPS ||
      cabs(temp.s1.c0 - temp2.s1.c0) > EPS || 
      cabs(temp.s1.c1 - temp2.s1.c1) > EPS || 
      cabs(temp.s1.c2 - temp2.s1.c2) > EPS ||
      cabs(temp.s2.c0 - temp2.s2.c0) > EPS || 
      cabs(temp.s2.c1 - temp2.s2.c1) > EPS || 
      cabs(temp.s2.c2 - temp2.s2.c2) > EPS ||
      cabs(temp.s3.c0 - temp2.s3.c0) > EPS || 
      cabs(temp.s3.c1 - temp2.s3.c1) > EPS || 
      cabs(temp.s3.c2 - temp2.s3.c2) > EPS )
    test = 1;

  assertFalseM(test, "D_x- failed\n");
  test = 0;

  _vector_add(psi,sp->s0,sp->s3);
  _su3_multiply(chi,u,psi);
  _complex_times_vector(psi,k,chi);
  _vector_add_assign(temp.s0,psi);
  _vector_add_assign(temp.s3,psi);
  _vector_sub(psi,sp->s1,sp->s2);
  _su3_multiply(chi,u,psi);
  _complex_times_vector(psi,k,chi);
  _vector_add_assign(temp.s1,psi);
  _vector_sub_assign(temp.s2,psi);

  vec_load2(r, &sp->s0);
  vec_load2(r+3, &sp->s1);
  vec_load2(r+9, &sp->s2);
  vec_load2(r+6, &sp->s3);
  vec_add2(r, &r[6]);
  vec_sub2(r+3, &r[9]);
  vec_su3_multiply_double2(&u, U, r);
  vec_cmplx_mul_double2(r, &r[6], U, &k);
  vec_add_double2(rs, r);
  vec_sub2(&rs[6], &r[3]);
  vec_add2(&rs[9], r);

  vec_store2(&temp2.s0, rs);
  vec_store2(&temp2.s1, rs+3);
  vec_store2(&temp2.s2, rs+6);
  vec_store2(&temp2.s3, rs+9);

  if( cabs(temp.s0.c0 - temp2.s0.c0) > EPS || 
      cabs(temp.s0.c1 - temp2.s0.c1) > EPS || 
      cabs(temp.s0.c2 - temp2.s0.c2) > EPS ||
      cabs(temp.s1.c0 - temp2.s1.c0) > EPS || 
      cabs(temp.s1.c1 - temp2.s1.c1) > EPS || 
      cabs(temp.s1.c2 - temp2.s1.c2) > EPS ||
      cabs(temp.s2.c0 - temp2.s2.c0) > EPS || 
      cabs(temp.s2.c1 - temp2.s2.c1) > EPS || 
      cabs(temp.s2.c2 - temp2.s2.c2) > EPS ||
      cabs(temp.s3.c0 - temp2.s3.c0) > EPS || 
      cabs(temp.s3.c1 - temp2.s3.c1) > EPS || 
      cabs(temp.s3.c2 - temp2.s3.c2) > EPS )
    test = 1;

  assertFalseM(test, "D_y+ failed\n");
  test = 0;

  _vector_sub(psi,sm->s0,sm->s3);
  _su3_inverse_multiply(chi,u,psi);
  _complexcjg_times_vector(psi,k,chi);
  _vector_add_assign(temp.s0,psi);
  _vector_sub_assign(temp.s3,psi);
  _vector_add(psi,sm->s1,sm->s2);
  _su3_inverse_multiply(chi,u,psi);
  _complexcjg_times_vector(psi,k,chi);
  _vector_add_assign(temp.s1,psi);
  _vector_add_assign(temp.s2,psi);

  vec_load2(r, &sm->s0);
  vec_load2(r+3, &sm->s1);
  vec_load2(r+9, &sm->s2);
  vec_load2(r+6, &sm->s3);
  vec_sub2(r, r+6);
  vec_add2(r+3, r+9);
  vec_su3_inverse_multiply_double2(&u, U, r);
  vec_cmplxcg_mul_double2(r, &r[6], U, &k);
  vec_add_double2(rs, r);
  vec_add2(rs+6, r+3);
  vec_sub2(rs+9, r);

  vec_store2(&temp2.s0, rs);
  vec_store2(&temp2.s1, rs+3);
  vec_store2(&temp2.s2, rs+6);
  vec_store2(&temp2.s3, rs+9);

  if( cabs(temp.s0.c0 - temp2.s0.c0) > EPS || 
      cabs(temp.s0.c1 - temp2.s0.c1) > EPS || 
      cabs(temp.s0.c2 - temp2.s0.c2) > EPS ||
      cabs(temp.s1.c0 - temp2.s1.c0) > EPS || 
      cabs(temp.s1.c1 - temp2.s1.c1) > EPS || 
      cabs(temp.s1.c2 - temp2.s1.c2) > EPS ||
      cabs(temp.s2.c0 - temp2.s2.c0) > EPS || 
      cabs(temp.s2.c1 - temp2.s2.c1) > EPS || 
      cabs(temp.s2.c2 - temp2.s2.c2) > EPS ||
      cabs(temp.s3.c0 - temp2.s3.c0) > EPS || 
      cabs(temp.s3.c1 - temp2.s3.c1) > EPS || 
      cabs(temp.s3.c2 - temp2.s3.c2) > EPS )
    test = 1;

  assertFalseM(test, "D_y- failed\n");
  test = 0;

  _vector_i_add(psi,sp->s0,sp->s2);
  _su3_multiply(chi,u,psi);
  _complex_times_vector(psi,k,chi);
  _vector_add_assign(temp.s0,psi);
  _vector_i_sub_assign(temp.s2,psi);
  _vector_i_sub(psi,sp->s1,sp->s3);
  _su3_multiply(chi,u,psi);
  _complex_times_vector(psi,k,chi);
  _vector_add_assign(temp.s1,psi);
  _vector_i_add_assign(temp.s3,psi);

  vec_load2(r, &sp->s0);
  vec_load2(r+3, &sp->s1);
  vec_load2(r+6, &sp->s2);
  vec_load2(r+9, &sp->s3);
  vec_i_mul_add2(r, r+6, U);
  vec_i_mul_sub2(r+3, r+9, U);
  vec_su3_multiply_double2(&u, U, r);
  vec_cmplx_mul_double2(r, &r[6], U, &k);
  vec_add_double2(rs, r);
  vec_i_mul_sub2(rs+6, r, U);
  vec_i_mul_add2(rs+9, r+3, U);

  vec_store2(&temp2.s0, rs);
  vec_store2(&temp2.s1, rs+3);
  vec_store2(&temp2.s2, rs+6);
  vec_store2(&temp2.s3, rs+9);

  if( cabs(temp.s0.c0 - temp2.s0.c0) > EPS || 
      cabs(temp.s0.c1 - temp2.s0.c1) > EPS || 
      cabs(temp.s0.c2 - temp2.s0.c2) > EPS ||
      cabs(temp.s1.c0 - temp2.s1.c0) > EPS || 
      cabs(temp.s1.c1 - temp2.s1.c1) > EPS || 
      cabs(temp.s1.c2 - temp2.s1.c2) > EPS ||
      cabs(temp.s2.c0 - temp2.s2.c0) > EPS || 
      cabs(temp.s2.c1 - temp2.s2.c1) > EPS || 
      cabs(temp.s2.c2 - temp2.s2.c2) > EPS ||
      cabs(temp.s3.c0 - temp2.s3.c0) > EPS || 
      cabs(temp.s3.c1 - temp2.s3.c1) > EPS || 
      cabs(temp.s3.c2 - temp2.s3.c2) > EPS )
    test = 1;

  assertFalseM(test, "D_z+ failed\n");
  test = 0;

  _vector_i_sub(psi,sm->s0,sm->s2);
  _su3_inverse_multiply(chi,u,psi);
  _complexcjg_times_vector(psi,k,chi);
  _vector_add_assign(temp.s0, psi);
  _vector_i_add_assign(temp.s2, psi);
  _vector_i_add(psi,sm->s1,sm->s3);
  _su3_inverse_multiply(chi,u,psi);
  _complexcjg_times_vector(psi,k,chi);
  _vector_add_assign(temp.s1, psi);
  _vector_i_sub_assign(temp.s3, psi);

  vec_load2(r, &sm->s0);
  vec_load2(r+3, &sm->s1);
  vec_load2(r+6, &sm->s2);
  vec_load2(r+9, &sm->s3);
  vec_i_mul_sub2(r, r+6, U);
  vec_i_mul_add2(r+3, r+9, U);
  vec_su3_inverse_multiply_double2(&u, U, r);
  vec_cmplxcg_mul_double2(r, &r[6], U, &k);
  vec_add_double2(rs, r);
  vec_store2(&temp2.s0, rs);
  vec_store2(&temp2.s1, rs+3);
  vec_i_mul_add2(rs+6, r, U);
  vec_store2(&temp2.s2, rs+6);
  vec_i_mul_sub2(rs+9, r+3, U);
  vec_store2(&temp2.s3, rs+9);

  vec_store2(&temp2.s0, rs);
  vec_store2(&temp2.s1, rs+3);
  vec_store2(&temp2.s2, rs+6);
  vec_store2(&temp2.s3, rs+9);

  if( cabs(temp.s0.c0 - temp2.s0.c0) > EPS || 
      cabs(temp.s0.c1 - temp2.s0.c1) > EPS || 
      cabs(temp.s0.c2 - temp2.s0.c2) > EPS ||
      cabs(temp.s1.c0 - temp2.s1.c0) > EPS || 
      cabs(temp.s1.c1 - temp2.s1.c1) > EPS || 
      cabs(temp.s1.c2 - temp2.s1.c2) > EPS ||
      cabs(temp.s2.c0 - temp2.s2.c0) > EPS || 
      cabs(temp.s2.c1 - temp2.s2.c1) > EPS || 
      cabs(temp.s2.c2 - temp2.s2.c2) > EPS ||
      cabs(temp.s3.c0 - temp2.s3.c0) > EPS || 
      cabs(temp.s3.c1 - temp2.s3.c1) > EPS || 
      cabs(temp.s3.c2 - temp2.s3.c2) > EPS )
    test = 1;

  assertFalseM(test, "D_z- failed\n");
  test = 0;

#else
  test = 1;
  assertFalseM(test, "not on a BG/Q or not compiling with XLC\n");
#endif
}

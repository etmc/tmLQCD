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
  complex double ca, cb, cc, cd, k;
  static vector4double a, b, c, d;
  static su3 u;
  static su3_vector phi0, phi1, phi2, phi3, phi4, phi5;
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

#else
  test = 1;
  assertFalseM(test, "not on a BG/Q or not compiling with XLC\n");
#endif
}

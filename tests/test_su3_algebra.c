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

#define EPS 5e-16

TEST(su3_multiply) {
  su3 u;
  su3_vector phi0, phi1, phi2, phi3;
  int test = 0;
  
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

  if( (creal(phi2.c0) - 2.441800e-01) > EPS || (cimag(phi2.c0) - 7.044200e-01) > EPS ||
      (creal(phi2.c1) - 5.417900e-01) > EPS || (cimag(phi2.c1) + 1.914840e+00) > EPS ||
      (creal(phi2.c2) - 1.380940e+00) > EPS || (cimag(phi2.c2) - 9.151800e-01) > EPS ||
      (creal(phi3.c0) + 5.020400e-02) > EPS || (cimag(phi3.c0) - 2.228320e-01) > EPS ||
      (creal(phi3.c1) - 3.445000e-01) > EPS || (cimag(phi3.c1) + 2.542120e-01) > EPS ||
      (creal(phi3.c2) - 6.088960e-01) > EPS || (cimag(phi3.c2) - 7.267380e-01) > EPS)
    test = 1;
  assertFalseM(test, "_su3_multiply failed\n.");
}

TEST(su3_inverse_multiply) {
  su3 u;
  su3_vector phi0, phi1, phi2, phi3;
  int test = 0;
  
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

  _su3_inverse_multiply(phi2, u, phi0);
  _su3_inverse_multiply(phi3, u, phi1);

  if( (creal(phi2.c0) + 1.668100e-01) > EPS || (cimag(phi2.c0) - 1.248950e+00) > EPS ||
      (creal(phi2.c1) - 2.116400e-01) > EPS || (cimag(phi2.c1) + 1.931510e+00) > EPS ||
      (creal(phi2.c2) - 6.485200e-01) > EPS || (cimag(phi2.c2) - 1.214960e+00) > EPS ||
      (creal(phi3.c0) + 3.095240e-01) > EPS || (cimag(phi3.c0) - 2.159480e-01) > EPS ||
      (creal(phi3.c1) - 3.796020e-01) > EPS || (cimag(phi3.c1) + 4.072300e-01) > EPS ||
      (creal(phi3.c2) - 3.496240e-01) > EPS || (cimag(phi3.c2) - 7.482380e-01) > EPS)
    test = 1;
  assertFalseM(test, "_su3_multiply failed\n.");
}

TEST(vector_add) {
  su3_vector phi0, phi1, phi2;
  int test = 0;
  
  phi0.c0 = 0.9+2*I;
  phi0.c1 = 0.6+1.3*I;
  phi0.c2 = 0.6-.2*I;

  phi1.c0 = -0.1+1*I;
  phi1.c1 = 0.1+0.3*I;
  phi1.c2 = -0.1-.12*I;

  _vector_add(phi2, phi0, phi1);

  if( cabs(phi2.c0 - phi0.c0 - phi1.c0) > EPS ||
      cabs(phi2.c1 - phi0.c1 - phi1.c1) > EPS ||
      cabs(phi2.c2 - phi0.c2 - phi1.c2) > EPS )
    test = 1;
  assertFalseM(test, "_vector_add failed\n.");
}

TEST(vector_sub) {
  su3_vector phi0, phi1, phi2;
  int test = 0;
  
  phi0.c0 = 0.9+2*I;
  phi0.c1 = 0.6+1.3*I;
  phi0.c2 = 0.6-.2*I;

  phi1.c0 = -0.1+1*I;
  phi1.c1 = 0.1+0.3*I;
  phi1.c2 = -0.1-.12*I;

  _vector_sub(phi2, phi0, phi1);

  if( cabs(phi2.c0 - phi0.c0 + phi1.c0) > EPS ||
      cabs(phi2.c1 - phi0.c1 + phi1.c1) > EPS ||
      cabs(phi2.c2 - phi0.c2 + phi1.c2) > EPS )
    test = 1;
  assertFalseM(test, "_vector_sub failed\n.");
}

TEST(vector_i_add) {
  su3_vector phi0, phi1, phi2;
  int test = 0;
  
  phi0.c0 = 0.9+2*I;
  phi0.c1 = 0.6+1.3*I;
  phi0.c2 = 0.6-.2*I;

  phi1.c0 = -0.1+1*I;
  phi1.c1 = 0.1+0.3*I;
  phi1.c2 = -0.1-.12*I;

  _vector_i_add(phi2, phi0, phi1);

  if( cabs(phi2.c0 - phi0.c0 - I*phi1.c0) > EPS ||
      cabs(phi2.c1 - phi0.c1 - I*phi1.c1) > EPS ||
      cabs(phi2.c2 - phi0.c2 - I*phi1.c2) > EPS)
    test = 1;
  assertFalseM(test, "_vector_i_add failed\n.");
}

TEST(vector_i_sub) {
  su3_vector phi0, phi1, phi2;
  int test = 0;
  
  phi0.c0 = 0.9+2*I;
  phi0.c1 = 0.6+1.3*I;
  phi0.c2 = 0.6-.2*I;

  phi1.c0 = -0.1+1*I;
  phi1.c1 = 0.1+0.3*I;
  phi1.c2 = -0.1-.12*I;

  _vector_i_sub(phi2, phi0, phi1);

  if( cabs(phi2.c0 - phi0.c0 + I*phi1.c0) > EPS ||
      cabs(phi2.c1 - phi0.c1 + I*phi1.c1) > EPS ||
      cabs(phi2.c2 - phi0.c2 + I*phi1.c2) > EPS)
    test = 1;
  assertFalseM(test, "_vector_i_sub failed\n.");
}

TEST(cmplx_times_vector) {
  su3_vector phi0, phi1;
  complex double c = 3.7655 - 0.3*I;
  int test = 0;
    
  phi0.c0 = 0.9+2*I;
  phi0.c1 = 0.6+1.3*I;
  phi0.c2 = 0.6-.2*I;

  _complex_times_vector(phi1, c, phi0);

  if( cabs(phi1.c0 - c*phi0.c0) > EPS ||
      cabs(phi1.c1 - c*phi0.c1) > EPS ||
      cabs(phi1.c2 - c*phi0.c2) > EPS)
    test = 1;
  assertFalseM(test, "_complex_times_vector failed\n.");
}

TEST(cmplxcjg_times_vector) {
  su3_vector phi0, phi1;
  complex double c = 3.7655 - 0.3*I;
  int test = 0;
    
  phi0.c0 = 0.9+2*I;
  phi0.c1 = 0.6+1.3*I;
  phi0.c2 = 0.6-.2*I;

  _complexcjg_times_vector(phi1, c, phi0);

  if( cabs(phi1.c0 - conj(c)*phi0.c0) > EPS ||
      cabs(phi1.c1 - conj(c)*phi0.c1) > EPS ||
      cabs(phi1.c2 - conj(c)*phi0.c2) > EPS)
    test = 1;
  assertFalseM(test, "_complexcjg_times_vector failed\n.");
}


TEST(su3_assign) {
  su3 m1,m2;
  
  int test = 0;
  m1.c00 = 1 + 1.*I; m1.c01 = 0.; m1.c02 = 0.;
  m1.c10 = 0.; m1.c11 = 1 + 1.*I; m1.c12 = 0.;
  m1.c20 = 0.; m1.c21 = 0.; m1.c22 = 1 + 1.*I;

  _su3_assign(m2,m1);

  if( creal(m2.c00) == 1 && cimag(m2.c00) == 1 && creal(m2.c01) == 0 && cimag(m2.c01) == 0 && creal(m2.c02) == 0 && cimag(m2.c02) == 0 &&
      creal(m2.c10) == 0 && cimag(m2.c10) == 0 && creal(m2.c11) == 1 && cimag(m2.c11) == 1 && creal(m2.c12) == 0 && cimag(m2.c12) == 0 &&
      creal(m2.c20) == 0 && cimag(m2.c20) == 0 && creal(m2.c21) == 0 && cimag(m2.c21) == 0 && creal(m2.c22) == 1 && cimag(m2.c22) == 1 )
    test = 1;
  
  assertTrueM(test,"The SU3 assignment operator does not work correctly!\n");
}

TEST(su3_expo_positivedet) {
  su3 Q, U;
  su3adj T;

  int test = 0;

  /* Positive determinant */
  Q.c00 = -0.2994;
  Q.c01 =  0.5952 + 1.3123*I;
  Q.c02 = -0.7943 + 0.0913*I;
  Q.c11 = -1.1430;
  Q.c12 = -2.0025 + 0.2978*I;
  Q.c22 = +1.4424;
  Q.c10 = conj(Q.c01);
  Q.c20 = conj(Q.c02);
  Q.c21 = conj(Q.c12);
  
  /* Matlab's solution for U = exp(i * Q) */
  U.c00 = +0.3391 -0.1635*I;
  U.c01 = -0.2357 +0.5203*I;
  U.c02 = +0.5609 +0.4663*I;
  U.c10 = -0.0740 -0.4204*I;
  U.c11 = -0.7706 -0.1863*I;
  U.c12 = +0.1191 -0.4185*I;
  U.c20 = +0.5351 -0.6243*I;
  U.c21 = +0.1825 +0.1089*I;
  U.c22 = -0.5279 -0.0022*I;

  _trace_lambda(T,Q);
  exposu3(&Q,&T);

  if( creal(Q.c00 - U.c00) > EPS && creal(Q.c01 - U.c01) > EPS && creal(Q.c02 - U.c02) > EPS &&
      creal(Q.c10 - U.c10) > EPS && creal(Q.c11 - U.c11) > EPS && creal(Q.c12 - U.c12) > EPS &&   
      creal(Q.c20 - U.c20) > EPS && creal(Q.c21 - U.c21) > EPS && creal(Q.c22 - U.c22) > EPS &&
      cimag(Q.c00 - U.c00) > EPS && cimag(Q.c01 - U.c01) > EPS && cimag(Q.c02 - U.c02) > EPS &&
      cimag(Q.c10 - U.c10) > EPS && cimag(Q.c11 - U.c11) > EPS && cimag(Q.c12 - U.c12) > EPS &&   
      cimag(Q.c20 - U.c20) > EPS && cimag(Q.c21 - U.c21) > EPS && cimag(Q.c22 - U.c22) > EPS )
    test = 1;

  assertFalseM(test,"The exponentation of Q with a positive determinant failed.\n");
}

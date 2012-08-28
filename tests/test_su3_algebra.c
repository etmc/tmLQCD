#include <config.h>
#include <complex.h>
#include <cu/cu.h>

#include "../su3.h"
#include "../su3adj.h"
#include "../expo.h"

#define EPS 5e-16

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
  Q = exposu3(T);

  if( creal(Q.c00 - U.c00) > EPS && creal(Q.c01 - U.c01) > EPS && creal(Q.c02 - U.c02) > EPS &&
      creal(Q.c10 - U.c10) > EPS && creal(Q.c11 - U.c11) > EPS && creal(Q.c12 - U.c12) > EPS &&   
      creal(Q.c20 - U.c20) > EPS && creal(Q.c21 - U.c21) > EPS && creal(Q.c22 - U.c22) > EPS &&
      cimag(Q.c00 - U.c00) > EPS && cimag(Q.c01 - U.c01) > EPS && cimag(Q.c02 - U.c02) > EPS &&
      cimag(Q.c10 - U.c10) > EPS && cimag(Q.c11 - U.c11) > EPS && cimag(Q.c12 - U.c12) > EPS &&   
      cimag(Q.c20 - U.c20) > EPS && cimag(Q.c21 - U.c21) > EPS && cimag(Q.c22 - U.c22) > EPS )
    test = 1;

  assertFalseM(test,"The exponentation of Q with a positive determinant failed.\n");
}

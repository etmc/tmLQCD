#include <config.h>

#include <cu/cu.h>

#include "../su3.h"
#include "../su3adj.h"
#include "../expo.h"

#define EPS 5e-16

TEST(su3_assign) {
  su3 m1,m2;
  
  int test = 0;

  m1.c00.re = 1; m1.c00.im = 1;  m1.c01.re = 0; m1.c01.im = 0;  m1.c02.re=0; m1.c02.im=0;
  m1.c10.re = 0; m1.c10.im = 0;  m1.c11.re = 1; m1.c11.im = 1;  m1.c12.re=0; m1.c12.im=0;
  m1.c20.re = 0; m1.c20.im = 0;  m1.c21.re = 0; m1.c21.im = 0;  m1.c22.re=1; m1.c22.im=1;

  _su3_assign(m2,m1);

  if( m2.c00.re == 1 && m2.c00.im == 1 && m2.c01.re == 0 && m2.c01.im == 0 && m2.c02.re == 0 && m2.c02.im == 0 &&
 m2.c10.re == 0 && m2.c10.im == 0 && m2.c11.re == 1 && m2.c11.im == 1 && m2.c12.re == 0 && m2.c12.im == 0 &&
 m2.c20.re == 0 && m2.c20.im == 0 && m2.c21.re == 0 && m2.c21.im == 0 && m2.c22.re == 1 && m2.c22.im == 1 )
    test = 1;

  assertTrueM(test,"The SU3 assignment operator does not work correctly!\n");
}

TEST(su3_expo_positivedet) {
  su3 Q, U;
  su3adj T;

  int test = 0;

  /* Positive determinant */
  Q.c00.re = -0.2994; Q.c00.im = 0.0;
  Q.c01.re =  0.5952; Q.c01.im = 1.3123;
  Q.c02.re = -0.7943; Q.c02.im = 0.0913;
  Q.c11.re = -1.1430; Q.c11.im = 0.0;
  Q.c12.re = -2.0025; Q.c12.im = 0.2978;
  Q.c22.re = +1.4424; Q.c22.im = 0.0;
  Q.c10.re = Q.c01.re; Q.c10.im = -Q.c01.im;
  Q.c20.re = Q.c02.re; Q.c20.im = -Q.c02.im;
  Q.c21.re = Q.c12.re; Q.c21.im = -Q.c12.im;
  
  /* Matlab's solution for U = exp(i * Q) */
  U.c00.re = +0.3391; U.c00.im = -0.1635;
  U.c01.re = -0.2357; U.c01.im = +0.5203;
  U.c02.re = +0.5609; U.c02.im = +0.4663;
  U.c10.re = -0.0740; U.c10.im = -0.4204;
  U.c11.re = -0.7706; U.c11.im = -0.1863;
  U.c12.re = +0.1191; U.c12.im = -0.4185;
  U.c20.re = +0.5351; U.c20.im = -0.6243;
  U.c21.re = +0.1825; U.c21.im = +0.1089;
  U.c22.re = -0.5279; U.c22.im = -0.0022;

  _trace_lambda(T,Q);
  Q = exposu3(T);

  if( Q.c00.re - U.c00.re > EPS &&  Q.c01.re - U.c01.re > EPS && Q.c02.re - U.c02.re > EPS &&
  Q.c10.re - U.c10.re > EPS && Q.c11.re - U.c11.re > EPS && Q.c12.re - U.c12.re > EPS &&   
  Q.c20.re - U.c20.re > EPS && Q.c21.re - U.c21.re > EPS && Q.c22.re - U.c22.re > EPS &&
  Q.c00.im - U.c00.im > EPS && Q.c01.im - U.c01.im > EPS && Q.c02.im - U.c02.im > EPS &&
  Q.c10.im - U.c10.im > EPS && Q.c11.im - U.c11.im > EPS && Q.c12.im - U.c12.im > EPS &&   
  Q.c20.im - U.c20.im > EPS && Q.c21.im - U.c21.im > EPS && Q.c22.im - U.c22.im > EPS )
    test = 1;

  assertFalseM(test,"The exponentation of Q with a positive determinant failed.\n");
}

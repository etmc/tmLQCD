#include <config.h>

#include "../cu/cu.h"

#include "../su3.h"


TEST(su3_assign) {
  su3 m1,m2;
  
  int test = 0;

  m1.c00.re = 1; m1.c00.im = 1;  m1.c01.re = 0; m1.c01.im = 0;  m1.c02.re=0; m1.c02.im=0;
  m1.c10.re = 0; m1.c10.im = 1;  m1.c11.re = 1; m1.c11.im = 1;  m1.c12.re=0; m1.c12.im=0;
  m1.c20.re = 0; m1.c20.im = 0;  m1.c21.re = 0; m1.c21.im = 0;  m1.c22.re=1; m1.c22.im=1;

  _su3_assign(m2,m1);

  if( m2.c00.re == 1 && m2.c00.im == 1 && m2.c01.re == 0 && m2.c01.im == 0 && m2.c02.re == 0 && m2.c02.im == 0 &&
 m2.c10.re == 1 && m2.c10.im == 1 && m2.c11.re == 0 && m2.c11.im == 0 && m2.c12.re == 0 && m2.c12.im == 0 &&
 m2.c20.re == 1 && m2.c20.im == 1 && m2.c21.re == 0 && m2.c21.im == 0 && m2.c22.re == 0 && m2.c22.im == 0 )
    test = 1;

  assertFalseM(test,"The SU3 assignment operator does not work correctly!\n");
}


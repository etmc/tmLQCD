#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <config.h>
#include <complex.h>
#include <time.h>
#include <cu/cu.h>
#include "../su3.h"
#include "../linalg_eo.h"

#define EPS 1e-15

TEST(scalar_prod_real) {
  const int N = 1000;
  int test = 0;
  double snrm = 0., atime=0., etime=0.;
  double *s, *r;
  spinor R[N] ALIGN;
  spinor S[N] ALIGN;

  R[0].s0.c0 =-1.1+0.7*I;
  R[0].s0.c1 =-0.0-0.7*I;
  R[0].s0.c2 =-1.3+1.9*I;
  R[0].s1.c0 = 4.4-4.0*I;
  R[0].s1.c1 =-1.5-3.1*I;
  R[0].s1.c2 = 6.6+2.3*I;
  R[0].s2.c0 =-5.4-0.7*I;
  R[0].s2.c1 =-7.8-6.3*I;
  R[0].s2.c2 = 1.3+3.7*I;
  R[0].s3.c0 =-8.3-4.6*I;
  R[0].s3.c1 =-1.3+4.5*I;
  R[0].s3.c2 = 9.3-2.3*I;

  S[0].s0.c0 =+1.1-0.7*I;
  S[0].s0.c1 =-0.0-1.7*I;
  S[0].s0.c2 =-1.3+1.9*I;
  S[0].s1.c0 =-4.4-4.0*I;
  S[0].s1.c1 =-1.5-3.1*I;
  S[0].s1.c2 = 3.4+0.3*I;
  S[0].s2.c0 =-5.4-0.7*I;
  S[0].s2.c1 =+7.4-6.3*I;
  S[0].s2.c2 = 1.6+3.7*I;
  S[0].s3.c0 =+8.3-4.6*I;
  S[0].s3.c1 =-1.8-0.5*I;
  S[0].s3.c2 =-9.0+0.3*I;

  R[1].s0.c0 =-1.1+0.7*I;
  R[1].s0.c1 = 0.0-0.7*I;
  R[1].s0.c2 =-1.3+1.9*I;
  R[1].s1.c0 = 1.4-4.0*I;
  R[1].s1.c1 =-1.0-3.1*I;
  R[1].s1.c2 = 6.6+3.3*I;
  R[1].s2.c0 =-5.4+0.7*I;
  R[1].s2.c1 =+7.8-4.3*I;
  R[1].s2.c2 =-1.3+3.7*I;
  R[1].s3.c0 =-8.3-4.6*I;
  R[1].s3.c1 =-1.3+4.5*I;
  R[1].s3.c2 = 9.3-2.3*I;

  S[1].s0.c0 =+1.1-0.7*I;
  S[1].s0.c1 =+0.0-1.7*I;
  S[1].s0.c2 = 1.2+4.1*I;
  S[1].s1.c0 =-4.4-4.7*I;
  S[1].s1.c1 =-1.5-2.1*I;
  S[1].s1.c2 = 3.4-0.3*I;
  S[1].s2.c0 =-5.4-0.7*I;
  S[1].s2.c1 = 7.4+6.3*I;
  S[1].s2.c2 =-1.6+3.7*I;
  S[1].s3.c0 =+1.3-4.6*I;
  S[1].s3.c1 =-1.8-0.5*I;
  S[1].s3.c2 =-2.0+0.3*I;

  snrm = scalar_prod_r(R, S, 2, 0);
  if( (snrm - 4.584000e+01) > EPS) test = 1;
  assertFalseM(test, "scalar_prod_r failed\n.");
  for(int i = 0; i < N; i++) {
    s = (double*)(S+i);
    r = (double*)(R+i);
    for(int j = 0; j < 24; j++) {
      s[j] = (double)random()/(double)RAND_MAX;
      r[j] = (double)random()/(double)RAND_MAX;
    }
  }

  atime = (double)clock()/(double)(CLOCKS_PER_SEC);
  for(int i = 0; i < 10000; i++) {
    snrm += scalar_prod_r(R, S, N, 0);
  }
  etime = (double)clock()/(double)(CLOCKS_PER_SEC);
  printf("res %e\n\n", snrm);
  printf("time = %e\n", etime-atime);
}


TEST(snorm) {
  const int N = 1000;
  int test = 0;
  double snrm = 0., atime=0., etime=0.;
  double *s, *r;
  spinor R[N] ALIGN;
  spinor S[N] ALIGN;

  R[0].s0.c0 =-1.1+0.7*I;
  R[0].s0.c1 =-0.0-0.7*I;
  R[0].s0.c2 =-1.3+1.9*I;
  R[0].s1.c0 = 4.4-4.0*I;
  R[0].s1.c1 =-1.5-3.1*I;
  R[0].s1.c2 = 6.6+2.3*I;
  R[0].s2.c0 =-5.4-0.7*I;
  R[0].s2.c1 =-7.8-6.3*I;
  R[0].s2.c2 = 1.3+3.7*I;
  R[0].s3.c0 =-8.3-4.6*I;
  R[0].s3.c1 =-1.3+4.5*I;
  R[0].s3.c2 = 9.3-2.3*I;

  S[0].s0.c0 =+1.1-0.7*I;
  S[0].s0.c1 =-0.0-1.7*I;
  S[0].s0.c2 =-1.3+1.9*I;
  S[0].s1.c0 =-4.4-4.0*I;
  S[0].s1.c1 =-1.5-3.1*I;
  S[0].s1.c2 = 3.4+0.3*I;
  S[0].s2.c0 =-5.4-0.7*I;
  S[0].s2.c1 =+7.4-6.3*I;
  S[0].s2.c2 = 1.6+3.7*I;
  S[0].s3.c0 =+8.3-4.6*I;
  S[0].s3.c1 =-1.8-0.5*I;
  S[0].s3.c2 =-9.0+0.3*I;

  R[1].s0.c0 =-1.1+0.7*I;
  R[1].s0.c1 = 0.0-0.7*I;
  R[1].s0.c2 =-1.3+1.9*I;
  R[1].s1.c0 = 1.4-4.0*I;
  R[1].s1.c1 =-1.0-3.1*I;
  R[1].s1.c2 = 6.6+3.3*I;
  R[1].s2.c0 =-5.4+0.7*I;
  R[1].s2.c1 =+7.8-4.3*I;
  R[1].s2.c2 =-1.3+3.7*I;
  R[1].s3.c0 =-8.3-4.6*I;
  R[1].s3.c1 =-1.3+4.5*I;
  R[1].s3.c2 = 9.3-2.3*I;

  S[1].s0.c0 =+1.1-0.7*I;
  S[1].s0.c1 =+0.0-1.7*I;
  S[1].s0.c2 = 1.2+4.1*I;
  S[1].s1.c0 =-4.4-4.7*I;
  S[1].s1.c1 =-1.5-2.1*I;
  S[1].s1.c2 = 3.4-0.3*I;
  S[1].s2.c0 =-5.4-0.7*I;
  S[1].s2.c1 = 7.4+6.3*I;
  S[1].s2.c2 =-1.6+3.7*I;
  S[1].s3.c0 =+1.3-4.6*I;
  S[1].s3.c1 =-1.8-0.5*I;
  S[1].s3.c2 =-2.0+0.3*I;

  snrm = square_norm(R, 2, 0);
  printf("square norm = %.16e\n", snrm);
  if( (snrm - 8.7152999999999997e+02 ) > EPS) test = 1;
  assertFalseM(test, "square_norm failed\n.");
  for(int i = 0; i < N; i++) {
    s = (double*)(S+i);
    r = (double*)(R+i);
    for(int j = 0; j < 24; j++) {
      s[j] = (double)random()/(double)RAND_MAX;
      r[j] = (double)random()/(double)RAND_MAX;
    }
  }

  atime = (double)clock()/(double)(CLOCKS_PER_SEC);
  for(int i = 0; i < 10000; i++) {
    snrm += square_norm(R, N, 0);
  }
  etime = (double)clock()/(double)(CLOCKS_PER_SEC);
  printf("res %e\n\n", snrm);
  printf("time = %e\n", etime-atime);

}

TEST(sdiff) {
  const int N = 1000;
  int test = 0;
  double snrm = 0., atime=0., etime=0.;
  double *s, *r;
  spinor R[N] ALIGN;
  spinor S[N] ALIGN;
  spinor Q[N] ALIGN;

  R[0].s0.c0 =-1.1+0.7*I;
  R[0].s0.c1 =-0.0-0.7*I;
  R[0].s0.c2 =-1.3+1.9*I;
  R[0].s1.c0 = 4.4-4.0*I;
  R[0].s1.c1 =-1.5-3.1*I;
  R[0].s1.c2 = 6.6+2.3*I;
  R[0].s2.c0 =-5.4-0.7*I;
  R[0].s2.c1 =-7.8-6.3*I;
  R[0].s2.c2 = 1.3+3.7*I;
  R[0].s3.c0 =-8.3-4.6*I;
  R[0].s3.c1 =-1.3+4.5*I;
  R[0].s3.c2 = 9.3-2.3*I;

  S[0].s0.c0 =+1.1-0.7*I;
  S[0].s0.c1 =-0.0-1.7*I;
  S[0].s0.c2 =-1.3+1.9*I;
  S[0].s1.c0 =-4.4-4.0*I;
  S[0].s1.c1 =-1.5-3.1*I;
  S[0].s1.c2 = 3.4+0.3*I;
  S[0].s2.c0 =-5.4-0.7*I;
  S[0].s2.c1 =+7.4-6.3*I;
  S[0].s2.c2 = 1.6+3.7*I;
  S[0].s3.c0 =+8.3-4.6*I;
  S[0].s3.c1 =-1.8-0.5*I;
  S[0].s3.c2 =-9.0+0.3*I;

  R[1].s0.c0 =-1.1+0.7*I;
  R[1].s0.c1 = 0.0-0.7*I;
  R[1].s0.c2 =-1.3+1.9*I;
  R[1].s1.c0 = 1.4-4.0*I;
  R[1].s1.c1 =-1.0-3.1*I;
  R[1].s1.c2 = 6.6+3.3*I;
  R[1].s2.c0 =-5.4+0.7*I;
  R[1].s2.c1 =+7.8-4.3*I;
  R[1].s2.c2 =-1.3+3.7*I;
  R[1].s3.c0 =-8.3-4.6*I;
  R[1].s3.c1 =-1.3+4.5*I;
  R[1].s3.c2 = 9.3-2.3*I;

  S[1].s0.c0 =+1.1-0.7*I;
  S[1].s0.c1 =+0.0-1.7*I;
  S[1].s0.c2 = 1.2+4.1*I;
  S[1].s1.c0 =-4.4-4.7*I;
  S[1].s1.c1 =-1.5-2.1*I;
  S[1].s1.c2 = 3.4-0.3*I;
  S[1].s2.c0 =-5.4-0.7*I;
  S[1].s2.c1 = 7.4+6.3*I;
  S[1].s2.c2 =-1.6+3.7*I;
  S[1].s3.c0 =+1.3-4.6*I;
  S[1].s3.c1 =-1.8-0.5*I;
  S[1].s3.c2 =-2.0+0.3*I;

  diff(Q, R, S, 2);
  snrm = square_norm(Q, 2, 0);
  printf("diff %.16e %.16e\n", creal(Q[0].s0.c0), cimag(Q[1].s2.c1));
  printf("square_norm Q=R-S  = %.16e\n\n", snrm);
  if( (snrm - 1.4169700000000000e+03 ) > EPS || 
      (creal(Q[0].s0.c0) +2.2 )>EPS || (cimag(Q[1].s2.c1 + 10.6) > EPS)) test = 1;
  assertFalseM(test, "diff failed\n.");
  for(int i = 0; i < N; i++) {
    s = (double*)(S+i);
    r = (double*)(R+i);
    for(int j = 0; j < 24; j++) {
      s[j] = (double)random()/(double)RAND_MAX;
      r[j] = (double)random()/(double)RAND_MAX;
    }
  }

  atime = (double)clock()/(double)(CLOCKS_PER_SEC);
  for(int i = 0; i < 10000; i++) {
    diff(Q, R, S, N);
    diff(R, S, Q, N);
  }
  etime = (double)clock()/(double)(CLOCKS_PER_SEC);
  printf("time = %e\n", etime-atime);

}

TEST(aaddm_r) {
  const int N = 1000;
  int test = 0;
  double c = 0.756;
  double snrm = 0., atime=0., etime=0.;
  double *s, *r;
  spinor R[N] ALIGN;
  spinor S[N] ALIGN;

  R[0].s0.c0 =-1.1+0.7*I;
  R[0].s0.c1 =-0.0-0.7*I;
  R[0].s0.c2 =-1.3+1.9*I;
  R[0].s1.c0 = 4.4-4.0*I;
  R[0].s1.c1 =-1.5-3.1*I;
  R[0].s1.c2 = 6.6+2.3*I;
  R[0].s2.c0 =-5.4-0.7*I;
  R[0].s2.c1 =-7.8-6.3*I;
  R[0].s2.c2 = 1.3+3.7*I;
  R[0].s3.c0 =-8.3-4.6*I;
  R[0].s3.c1 =-1.3+4.5*I;
  R[0].s3.c2 = 9.3-2.3*I;

  S[0].s0.c0 =+1.1-0.7*I;
  S[0].s0.c1 =-0.0-1.7*I;
  S[0].s0.c2 =-1.3+1.9*I;
  S[0].s1.c0 =-4.4-4.0*I;
  S[0].s1.c1 =-1.5-3.1*I;
  S[0].s1.c2 = 3.4+0.3*I;
  S[0].s2.c0 =-5.4-0.7*I;
  S[0].s2.c1 =+7.4-6.3*I;
  S[0].s2.c2 = 1.6+3.7*I;
  S[0].s3.c0 =+8.3-4.6*I;
  S[0].s3.c1 =-1.8-0.5*I;
  S[0].s3.c2 =-9.0+0.3*I;

  R[1].s0.c0 =-1.1+0.7*I;
  R[1].s0.c1 = 0.0-0.7*I;
  R[1].s0.c2 =-1.3+1.9*I;
  R[1].s1.c0 = 1.4-4.0*I;
  R[1].s1.c1 =-1.0-3.1*I;
  R[1].s1.c2 = 6.6+3.3*I;
  R[1].s2.c0 =-5.4+0.7*I;
  R[1].s2.c1 =+7.8-4.3*I;
  R[1].s2.c2 =-1.3+3.7*I;
  R[1].s3.c0 =-8.3-4.6*I;
  R[1].s3.c1 =-1.3+4.5*I;
  R[1].s3.c2 = 9.3-2.3*I;

  S[1].s0.c0 =+1.1-0.7*I;
  S[1].s0.c1 =+0.0-1.7*I;
  S[1].s0.c2 = 1.2+4.1*I;
  S[1].s1.c0 =-4.4-4.7*I;
  S[1].s1.c1 =-1.5-2.1*I;
  S[1].s1.c2 = 3.4-0.3*I;
  S[1].s2.c0 =-5.4-0.7*I;
  S[1].s2.c1 = 7.4+6.3*I;
  S[1].s2.c2 =-1.6+3.7*I;
  S[1].s3.c0 =+1.3-4.6*I;
  S[1].s3.c1 =-1.8-0.5*I;
  S[1].s3.c2 =-2.0+0.3*I;

  assign_add_mul_r(R, S, c, 2);
  snrm = square_norm(R, 2, 0);
  printf("single parts %.16e %.16e\n", creal(R[0].s0.c0), cimag(R[1].s2.c1));
  printf("assign_add_mul_r  = %.16e\n\n", snrm);
  if( (snrm/1000. - 1.3049770963199999 ) > EPS || 
      (creal(R[0].s0.c0)*10 + 2.6840000000000003 ) > EPS || (cimag(R[1].s2.c1)*10 - 4.6280000000000010) > EPS) test = 1;
  assertFalseM(test, "assign_add_mul_r failed\n.");

  for(int i = 0; i < N; i++) {
    s = (double*)(S+i);
    r = (double*)(R+i);
    for(int j = 0; j < 24; j++) {
      s[j] = (double)random()/(double)RAND_MAX;
      r[j] = (double)random()/(double)RAND_MAX;
    }
  }

  atime = (double)clock()/(double)(CLOCKS_PER_SEC);
  for(int i = 0; i < 10000; i++) {
    assign_add_mul_r(R, S, c, N);
  }
  etime = (double)clock()/(double)(CLOCKS_PER_SEC);
  printf("time assign_add_mul_r = %e\n", etime-atime);
}

TEST(amadd_r) {
  const int N = 1000;
  int test = 0;
  double c = 0.756;
  double snrm = 0., atime=0., etime=0.;
  double *s, *r;
  spinor R[N] ALIGN;
  spinor S[N] ALIGN;

  R[0].s0.c0 =-1.1+0.7*I;
  R[0].s0.c1 =-0.0-0.7*I;
  R[0].s0.c2 =-1.3+1.9*I;
  R[0].s1.c0 = 4.4-4.0*I;
  R[0].s1.c1 =-1.5-3.1*I;
  R[0].s1.c2 = 6.6+2.3*I;
  R[0].s2.c0 =-5.4-0.7*I;
  R[0].s2.c1 =-7.8-6.3*I;
  R[0].s2.c2 = 1.3+3.7*I;
  R[0].s3.c0 =-8.3-4.6*I;
  R[0].s3.c1 =-1.3+4.5*I;
  R[0].s3.c2 = 9.3-2.3*I;

  S[0].s0.c0 =+1.1-0.7*I;
  S[0].s0.c1 =-0.0-1.7*I;
  S[0].s0.c2 =-1.3+1.9*I;
  S[0].s1.c0 =-4.4-4.0*I;
  S[0].s1.c1 =-1.5-3.1*I;
  S[0].s1.c2 = 3.4+0.3*I;
  S[0].s2.c0 =-5.4-0.7*I;
  S[0].s2.c1 =+7.4-6.3*I;
  S[0].s2.c2 = 1.6+3.7*I;
  S[0].s3.c0 =+8.3-4.6*I;
  S[0].s3.c1 =-1.8-0.5*I;
  S[0].s3.c2 =-9.0+0.3*I;

  R[1].s0.c0 =-1.1+0.7*I;
  R[1].s0.c1 = 0.0-0.7*I;
  R[1].s0.c2 =-1.3+1.9*I;
  R[1].s1.c0 = 1.4-4.0*I;
  R[1].s1.c1 =-1.0-3.1*I;
  R[1].s1.c2 = 6.6+3.3*I;
  R[1].s2.c0 =-5.4+0.7*I;
  R[1].s2.c1 =+7.8-4.3*I;
  R[1].s2.c2 =-1.3+3.7*I;
  R[1].s3.c0 =-8.3-4.6*I;
  R[1].s3.c1 =-1.3+4.5*I;
  R[1].s3.c2 = 9.3-2.3*I;

  S[1].s0.c0 =+1.1-0.7*I;
  S[1].s0.c1 =+0.0-1.7*I;
  S[1].s0.c2 = 1.2+4.1*I;
  S[1].s1.c0 =-4.4-4.7*I;
  S[1].s1.c1 =-1.5-2.1*I;
  S[1].s1.c2 = 3.4-0.3*I;
  S[1].s2.c0 =-5.4-0.7*I;
  S[1].s2.c1 = 7.4+6.3*I;
  S[1].s2.c2 =-1.6+3.7*I;
  S[1].s3.c0 =+1.3-4.6*I;
  S[1].s3.c1 =-1.8-0.5*I;
  S[1].s3.c2 =-2.0+0.3*I;

  assign_mul_add_r(R, c, S, 2);
  snrm = square_norm(R, 2, 0);
  printf("single parts %.16e %.16e\n", creal(R[0].s0.c0), cimag(R[1].s2.c1));
  printf("assign_mul_add_r  = %.16e\n\n", snrm);
  if( (snrm/1000. - 1.2045408500799999 ) > EPS || 
      (creal(R[0].s0.c0)*10 - 2.6840000000000003 ) > EPS || (cimag(R[1].s2.c1) - 3.0491999999999999) > EPS) test = 1;
  assertFalseM(test, "assign_mul_add_r failed\n.");

  for(int i = 0; i < N; i++) {
    s = (double*)(S+i);
    r = (double*)(R+i);
    for(int j = 0; j < 24; j++) {
      s[j] = (double)random()/(double)RAND_MAX;
      r[j] = (double)random()/(double)RAND_MAX;
    }
  }

  atime = (double)clock()/(double)(CLOCKS_PER_SEC);
  for(int i = 0; i < 10000; i++) {
    assign_mul_add_r(R, c, S, N);
  }
  etime = (double)clock()/(double)(CLOCKS_PER_SEC);
  printf("time assign_mul_add_r = %e\n", etime-atime);
}

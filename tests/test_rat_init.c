#include <stdlib.h>
#include <stdio.h>
#include <config.h>
#include <complex.h>
#include <cu/cu.h>
#include <complex.h>
#include <math.h>
#include "../rational/rational.h"

#define EPS 1e-8

// apply rational approximation as partial fraction
double apply_R(const int order, const double y, const double A, double * rs, double * as) {
  double x = 1.;
  
  for(int i = 0; i < order; i++) {
    x += rs[i]/(y+as[i]);
  }
  return(A*x);
}

_Complex double apply_C(const int order, const double y, const double A, double * rs, double * as) {
  _Complex double x = 1.;
  for(int i = 0; i < order; i++) {
    x += I*rs[i]/(sqrt(y)+I*as[i]);
  }
  return(x/sqrt(A));
}

_Complex double apply_Cdag(const int order, const double y, const double A, double * rs, double * as) {
  _Complex double x = 1.;
  for(int i = 0; i < order; i++) {
    x -= I*rs[i]/(sqrt(y)-I*as[i]);
  }
  return(x/sqrt(A));
}

TEST(rat_init) {
  int t = 0, ret=0;
  int order = 10;
  double eps = 1.e-4;
  double ra = eps, rb = 1.;
  rational_t rat;
  double * ar = malloc(order*sizeof(double));
  rat.order = order;
  rat.range[0] = ra;
  rat.range[1] = rb;
  rat.crange[0] = 0;
  rat.crange[1] = order-1;
 
  ret = init_rational(&rat);
  assertFalseM(ret, "rat_init failed\n");

  for(int i = 0; i < order; i++) {
    ar[i] = (rat.mu[i])*(rat.mu[i]);
  }
  
  for(double y = eps; y < 1.; y += eps) {
    double x = apply_R(rat.order, y, rat.A, rat.rmu, ar);
    if(fabs(1 - x*sqrt(y)) > rat.delta + EPS) {
      t++;
      printf("%e %e %e %e %e\n", y, x, 1./sqrt(y), fabs(1 - x*sqrt(y)), rat.delta);
    }
  }
  assertFalseM(t, "rational approximation not as accurate as expected\n.");
  t = 0;

  for(double y = eps; y < 1.; y += eps) {
    _Complex double c0 = apply_C(rat.order, y, rat.A, rat.rnu, rat.nu);
    _Complex double c1 = apply_Cdag(rat.order, y, rat.A, rat.rnu, rat.nu);
    double x = apply_R(rat.order, y, rat.A, rat.rmu, ar);
    if((fabs(1-creal(c0*x*c1)) > rat.delta + EPS) || cimag(c0*c1) > EPS ) {
      t++;
      printf("res: %e %e %e %e (%e, %e) (%e, %e)\n", y, sqrt(y), x, creal(c0*c1), creal(c0), cimag(c0), creal(c1), cimag(c1));
    }
  }
  assertFalseM(t, "C^dagger C approximation not as accurate as expected\n.");
  t = 0;

}

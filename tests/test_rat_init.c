#include <stdlib.h>
#include <stdio.h>
#include <config.h>
#include <complex.h>
#include <cu/cu.h>
#include <math.h>
#include "../rational/rational.h"

#define EPS 1e-7

// apply rational approximation as partial fraction
double apply_R(const int order, const double y, const double A, double * rs, double * as) {
  double x = 1.;
  
  for(int i = 0; i < order; i++) {
    x += rs[i]/(y+as[i]);
  }
  return(A*x);
}


TEST(rat_init) {
  int t = 0, ret=0;
  int order = 10;
  double eps = 1.e-4;
  double ra = eps, rb = 1.;
  rational_t rat;
  double * ar = malloc(order*sizeof(double));
  
  ret = init_rational(&rat, order, ra, rb, 0, order-1);
  assertFalseM(ret, "rat_init failed\n");
  printf("%d\n", ret);

  for(int i = 0; i < order; i++) {
    ar[i] = (rat.mu[i])*(rat.mu[i])/rb;
  }
  
  for(double y = eps; y < 1.; y += eps) {
    double x = apply_R(rat.order, y, rat.A, rat.rmu, ar);
    if(fabs(fabs(x - 1./sqrt(y)) - rat.delta) > EPS) {
      t++;
      printf("%e %e %e %e %e\n", y, x, 1./sqrt(y), fabs(fabs(x - 1./sqrt(y)) - rat.delta), rat.delta);
    }
  }
  assertFalseM(t, "rational approximation not as accurate as expected\n.");
}

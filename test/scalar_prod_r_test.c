/* #ifndef apenext */
/* #include <complex.h> */
/* #endif */

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#ifndef _STD_C99_COMPLEX
#include "complex.h"
#endif

#include "qdran64.h"

#include "su3.h"
#include "linalg/scalar_prod_r.c"


#define N 16

int main(void) {
  double s;
  int i;
  spinor a[N], b[N];

  qdran64_init(42,13);

  for(i=0; i<N; i++) {
#ifndef _STD_C99_COMPLEX
    qdran64_2d(&(a[i].s0.c0.re),&(a[i].s0.c0.im));
    qdran64_2d(&(a[i].s0.c1.re),&(a[i].s0.c1.im));
    qdran64_2d(&(a[i].s0.c2.re),&(a[i].s0.c2.im));

    qdran64_2d(&(a[i].s1.c0.re),&(a[i].s1.c0.im));
    qdran64_2d(&(a[i].s1.c1.re),&(a[i].s1.c1.im));
    qdran64_2d(&(a[i].s1.c2.re),&(a[i].s1.c2.im));

    qdran64_2d(&(a[i].s2.c0.re),&(a[i].s2.c0.im));
    qdran64_2d(&(a[i].s2.c1.re),&(a[i].s2.c1.im));
    qdran64_2d(&(a[i].s2.c2.re),&(a[i].s2.c2.im));

    qdran64_2d(&(a[i].s3.c0.re),&(a[i].s3.c0.im));
    qdran64_2d(&(a[i].s3.c1.re),&(a[i].s3.c1.im));
    qdran64_2d(&(a[i].s3.c2.re),&(a[i].s3.c2.im));

    qdran64_2d(&(b[i].s0.c0.re),&(b[i].s0.c0.im));
    qdran64_2d(&(b[i].s0.c1.re),&(b[i].s0.c1.im));
    qdran64_2d(&(b[i].s0.c2.re),&(b[i].s0.c2.im));

    qdran64_2d(&(b[i].s1.c0.re),&(b[i].s1.c0.im));
    qdran64_2d(&(b[i].s1.c1.re),&(b[i].s1.c1.im));
    qdran64_2d(&(b[i].s1.c2.re),&(b[i].s1.c2.im));

    qdran64_2d(&(b[i].s2.c0.re),&(b[i].s2.c0.im));
    qdran64_2d(&(b[i].s2.c1.re),&(b[i].s2.c1.im));
    qdran64_2d(&(b[i].s2.c2.re),&(b[i].s2.c2.im));

    qdran64_2d(&(b[i].s3.c0.re),&(b[i].s3.c0.im));
    qdran64_2d(&(b[i].s3.c1.re),&(b[i].s3.c1.im));
    qdran64_2d(&(b[i].s3.c2.re),&(b[i].s3.c2.im));
#else
    qdran64z(&a[i].s0.c0);
    qdran64z(&a[i].s0.c1);
    qdran64z(&a[i].s0.c2);

    qdran64z(&a[i].s1.c0);
    qdran64z(&a[i].s1.c1);
    qdran64z(&a[i].s1.c2);

    qdran64z(&a[i].s2.c0);
    qdran64z(&a[i].s2.c1);
    qdran64z(&a[i].s2.c2);

    qdran64z(&a[i].s3.c0);
    qdran64z(&a[i].s3.c1);
    qdran64z(&a[i].s3.c2);

    qdran64z(&b[i].s0.c0);
    qdran64z(&b[i].s0.c1);
    qdran64z(&b[i].s0.c2);

    qdran64z(&b[i].s1.c0);
    qdran64z(&b[i].s1.c1);
    qdran64z(&b[i].s1.c2);

    qdran64z(&b[i].s2.c0);
    qdran64z(&b[i].s2.c1);
    qdran64z(&b[i].s2.c2);

    qdran64z(&b[i].s3.c0);
    qdran64z(&b[i].s3.c1);
    qdran64z(&b[i].s3.c2);
#endif
 }

#ifndef _STD_C99_COMPLEX
  printf("%e %e\n",a[0].s0.c0.re,a[0].s0.c0.im);
  printf("%e %e\n",b[0].s0.c0.re,b[0].s0.c0.im);
  printf("%e %e\n",a[N-1].s3.c2.re,a[N-1].s3.c2.im);
  printf("%e %e\n",b[N-1].s3.c2.re,b[N-1].s3.c2.im);
#else
  printf("%e %e\n",creal(a[0].s0.c0),cimag(a[0].s0.c0));
  printf("%e %e\n",creal(b[0].s0.c0),cimag(b[0].s0.c0));
  printf("%e %e\n",creal(a[N-1].s3.c2),cimag(a[N-1].s3.c2));
  printf("%e %e\n",creal(b[N-1].s3.c2),cimag(b[N-1].s3.c2));
#endif
  s=scalar_prod_r(a,b,N, 1);
  printf("scalar_prod_r(a,b,%d)=%1.16e\n",N,s);

  return 0;
}

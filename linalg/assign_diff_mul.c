/* $Id$ */

#include <stdlib.h>
#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include "su3.h"
#include "assign_diff_mul.h"

/* S=S-c*Q */
void assign_diff_mul(spinor * const S,spinor * const R, const complex c, const int N){
  int ix;
  double a,b;
  spinor *r, *s;

  a = - c.re;
  b = - c.im;

  for (ix=0;ix<N;ix++){
    s = (spinor *) S + ix;
    r = (spinor *) R + ix;

    (*s).s0.c0.re += a*(*r).s0.c0.re - b*(*r).s0.c0.im;
    (*s).s0.c0.im += a*(*r).s0.c0.im + b*(*r).s0.c0.re;
    (*s).s0.c1.re += a*(*r).s0.c1.re - b*(*r).s0.c1.im;
    (*s).s0.c1.im += a*(*r).s0.c1.im + b*(*r).s0.c1.re;
    (*s).s0.c2.re += a*(*r).s0.c2.re - b*(*r).s0.c2.im;
    (*s).s0.c2.im += a*(*r).s0.c2.im + b*(*r).s0.c2.re;

    (*s).s1.c0.re += a*(*r).s1.c0.re - b*(*r).s1.c0.im;
    (*s).s1.c0.im += a*(*r).s1.c0.im + b*(*r).s1.c0.re;
    (*s).s1.c1.re += a*(*r).s1.c1.re - b*(*r).s1.c1.im;
    (*s).s1.c1.im += a*(*r).s1.c1.im + b*(*r).s1.c1.re;
    (*s).s1.c2.re += a*(*r).s1.c2.re - b*(*r).s1.c2.im;
    (*s).s1.c2.im += a*(*r).s1.c2.im + b*(*r).s1.c2.re;

    (*s).s2.c0.re += a*(*r).s2.c0.re - b*(*r).s2.c0.im;
    (*s).s2.c0.im += a*(*r).s2.c0.im + b*(*r).s2.c0.re;
    (*s).s2.c1.re += a*(*r).s2.c1.re - b*(*r).s2.c1.im;
    (*s).s2.c1.im += a*(*r).s2.c1.im + b*(*r).s2.c1.re;
    (*s).s2.c2.re += a*(*r).s2.c2.re - b*(*r).s2.c2.im;
    (*s).s2.c2.im += a*(*r).s2.c2.im + b*(*r).s2.c2.re;
     
    (*s).s3.c0.re += a*(*r).s3.c0.re - b*(*r).s3.c0.im;
    (*s).s3.c0.im += a*(*r).s3.c0.im + b*(*r).s3.c0.re;
    (*s).s3.c1.re += a*(*r).s3.c1.re - b*(*r).s3.c1.im;
    (*s).s3.c1.im += a*(*r).s3.c1.im + b*(*r).s3.c1.re;
    (*s).s3.c2.re += a*(*r).s3.c2.re - b*(*r).s3.c2.im;
    (*s).s3.c2.im += a*(*r).s3.c2.im + b*(*r).s3.c2.re;
  }
}

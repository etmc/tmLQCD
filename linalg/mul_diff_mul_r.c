/* $Id$*/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "su3.h"
#include "mul_diff_mul_r.h"


/* Makes (*R)=c1*(*S)-c2*(*U) , c1 and c2 are complex constants */
void mul_diff_mul_r(spinor * const R, spinor * const S,spinor * const U,
		    const double c1, const double c2, const int N){
  int ix;
  spinor *r,*s,*u;
  
  for (ix=0;ix<N;ix++){
    r=(spinor *) R + ix;
    s=(spinor *) S + ix;
    u=(spinor *) U + ix;
    
    (*r).s0.c0.re = c1*(*s).s0.c0.re - c2*(*u).s0.c0.re;
    (*r).s0.c0.im = c1*(*s).s0.c0.im - c2*(*u).s0.c0.im;
    (*r).s0.c1.re = c1*(*s).s0.c1.re - c2*(*u).s0.c1.re;
    (*r).s0.c1.im = c1*(*s).s0.c1.im - c2*(*u).s0.c1.im;
    (*r).s0.c2.re = c1*(*s).s0.c2.re - c2*(*u).s0.c2.re;
    (*r).s0.c2.im = c1*(*s).s0.c2.im - c2*(*u).s0.c2.im;
    (*r).s1.c0.re = c1*(*s).s1.c0.re - c2*(*u).s1.c0.re;
    (*r).s1.c0.im = c1*(*s).s1.c0.im - c2*(*u).s1.c0.im;
    (*r).s1.c1.re = c1*(*s).s1.c1.re - c2*(*u).s1.c1.re;
    (*r).s1.c1.im = c1*(*s).s1.c1.im - c2*(*u).s1.c1.im;
    (*r).s1.c2.re = c1*(*s).s1.c2.re - c2*(*u).s1.c2.re;
    (*r).s1.c2.im = c1*(*s).s1.c2.im - c2*(*u).s1.c2.im;
    (*r).s2.c0.re = c1*(*s).s2.c0.re - c2*(*u).s2.c0.re;
    (*r).s2.c0.im = c1*(*s).s2.c0.im - c2*(*u).s2.c0.im;
    (*r).s2.c1.re = c1*(*s).s2.c1.re - c2*(*u).s2.c1.re;
    (*r).s2.c1.im = c1*(*s).s2.c1.im - c2*(*u).s2.c1.im;
    (*r).s2.c2.re = c1*(*s).s2.c2.re - c2*(*u).s2.c2.re;
    (*r).s2.c2.im = c1*(*s).s2.c2.im - c2*(*u).s2.c2.im;
    (*r).s3.c0.re = c1*(*s).s3.c0.re - c2*(*u).s3.c0.re;
    (*r).s3.c0.im = c1*(*s).s3.c0.im - c2*(*u).s3.c0.im;
    (*r).s3.c1.re = c1*(*s).s3.c1.re - c2*(*u).s3.c1.re;
    (*r).s3.c1.im = c1*(*s).s3.c1.im - c2*(*u).s3.c1.im;
    (*r).s3.c2.re = c1*(*s).s3.c2.re - c2*(*u).s3.c2.re;
    (*r).s3.c2.im = c1*(*s).s3.c2.im - c2*(*u).s3.c2.im;
  }
}




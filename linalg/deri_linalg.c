/* $Id$ */

#include <stdlib.h>
#include "su3.h"
#include "deri_linalg.h"

void deri_linalg(spinor * const R, const double c1, spinor * const S, const double c2, spinor * const U, const int N) {
  int ix;
  spinor *r,*s,*t;
  
  /* Change due to even-odd preconditioning : VOLUME   to VOLUME/2 */   
  for (ix = 0; ix < N; ix++) {
    r=R+ix;
    s=S+ix;
    t=U+ix;
    
    (*r).s0.c0.re = c1*(*s).s0.c0.re + c2*(*t).s0.c0.re;
    (*r).s0.c0.im = c1*(*s).s0.c0.im + c2*(*t).s0.c0.im;
    (*r).s0.c1.re = c1*(*s).s0.c1.re + c2*(*t).s0.c1.re;
    (*r).s0.c1.im = c1*(*s).s0.c1.im + c2*(*t).s0.c1.im;
    (*r).s0.c2.re = c1*(*s).s0.c2.re + c2*(*t).s0.c2.re;
    (*r).s0.c2.im = c1*(*s).s0.c2.im + c2*(*t).s0.c2.im;
    
    (*r).s1.c0.re = c1*(*s).s1.c0.re + c2*(*t).s1.c0.re;
    (*r).s1.c0.im = c1*(*s).s1.c0.im + c2*(*t).s1.c0.im;
    (*r).s1.c1.re = c1*(*s).s1.c1.re + c2*(*t).s1.c1.re;
    (*r).s1.c1.im = c1*(*s).s1.c1.im + c2*(*t).s1.c1.im;
    (*r).s1.c2.re = c1*(*s).s1.c2.re + c2*(*t).s1.c2.re;
    (*r).s1.c2.im = c1*(*s).s1.c2.im + c2*(*t).s1.c2.im;
    
    (*r).s2.c0.re = c1*(*s).s2.c0.re - c2*(*t).s2.c0.re;
    (*r).s2.c0.im = c1*(*s).s2.c0.im - c2*(*t).s2.c0.im;
    (*r).s2.c1.re = c1*(*s).s2.c1.re - c2*(*t).s2.c1.re;
    (*r).s2.c1.im = c1*(*s).s2.c1.im - c2*(*t).s2.c1.im;
    (*r).s2.c2.re = c1*(*s).s2.c2.re - c2*(*t).s2.c2.re;
    (*r).s2.c2.im = c1*(*s).s2.c2.im - c2*(*t).s2.c2.im;
    
    (*r).s3.c0.re = c1*(*s).s3.c0.re - c2*(*t).s3.c0.re;
    (*r).s3.c0.im = c1*(*s).s3.c0.im - c2*(*t).s3.c0.im;
    (*r).s3.c1.re = c1*(*s).s3.c1.re - c2*(*t).s3.c1.re;
    (*r).s3.c1.im = c1*(*s).s3.c1.im - c2*(*t).s3.c1.im;
    (*r).s3.c2.re = c1*(*s).s3.c2.re - c2*(*t).s3.c2.re;
    (*r).s3.c2.im = c1*(*s).s3.c2.im - c2*(*t).s3.c2.im;
  }
}

/* $Id$ */

#include <stdlib.h>
#include "su3.h"
#include "sse.h"
#include "assign_mul_add_r.h"

#if defined SSE2
/* k input , l output*/
void assign_mul_add_r(spinor * const S, const double c, spinor * const R, const int N) {

  int ix;
  su3_vector *s,*r;
  __asm__ __volatile__ ("movsd %0, %%xmm7 \n\t"
			"unpcklpd %%xmm7, %%xmm7"
			:
			:
			"m" (c));
  s=&S[0].s0;
  r=&R[0].s0;
  for (ix=0;ix<4*N;ix++) {
    _sse_load(*s);
    __asm__ __volatile__ ("mulpd %%xmm7, %%xmm0 \n\t"
			  "mulpd %%xmm7, %%xmm1 \n\t"
			  "mulpd %%xmm7, %%xmm2"
			  :
			  :);
    _sse_load_up(*r);
    __asm__ __volatile__ ("addpd %%xmm3, %%xmm0 \n\t"
			  "addpd %%xmm4, %%xmm1 \n\t"
			  "addpd %%xmm5, %%xmm2"
			  :
			  :);
    _sse_store(*s);
    s++; r++;
  }
}
#else
/* k input , l output*/
void assign_mul_add_r(spinor * const R, const double c, spinor * const S, const int N) {

  int ix;
  static double fact;
  spinor *r,*s;
  
  fact=c;
  
  /* Change due to even-odd preconditioning : VOLUME   to VOLUME/2 */   
  for (ix = 0; ix < N; ix++) {
    r = R + ix;
    s = S + ix;
    
    (*r).s0.c0.re = fact*(*r).s0.c0.re + (*s).s0.c0.re;
    (*r).s0.c0.im = fact*(*r).s0.c0.im + (*s).s0.c0.im;
    (*r).s0.c1.re = fact*(*r).s0.c1.re + (*s).s0.c1.re;
    (*r).s0.c1.im = fact*(*r).s0.c1.im + (*s).s0.c1.im;
    (*r).s0.c2.re = fact*(*r).s0.c2.re + (*s).s0.c2.re;
    (*r).s0.c2.im = fact*(*r).s0.c2.im + (*s).s0.c2.im;
    
    (*r).s1.c0.re = fact*(*r).s1.c0.re + (*s).s1.c0.re;
    (*r).s1.c0.im = fact*(*r).s1.c0.im + (*s).s1.c0.im;
    (*r).s1.c1.re = fact*(*r).s1.c1.re + (*s).s1.c1.re;
    (*r).s1.c1.im = fact*(*r).s1.c1.im + (*s).s1.c1.im;
    (*r).s1.c2.re = fact*(*r).s1.c2.re + (*s).s1.c2.re;
    (*r).s1.c2.im = fact*(*r).s1.c2.im + (*s).s1.c2.im;         
    
    (*r).s2.c0.re = fact*(*r).s2.c0.re + (*s).s2.c0.re;
    (*r).s2.c0.im = fact*(*r).s2.c0.im + (*s).s2.c0.im;
    (*r).s2.c1.re = fact*(*r).s2.c1.re + (*s).s2.c1.re;
    (*r).s2.c1.im = fact*(*r).s2.c1.im + (*s).s2.c1.im;
    (*r).s2.c2.re = fact*(*r).s2.c2.re + (*s).s2.c2.re;
    (*r).s2.c2.im = fact*(*r).s2.c2.im + (*s).s2.c2.im;
    
    (*r).s3.c0.re = fact*(*r).s3.c0.re + (*s).s3.c0.re;
    (*r).s3.c0.im = fact*(*r).s3.c0.im + (*s).s3.c0.im;
    (*r).s3.c1.re = fact*(*r).s3.c1.re + (*s).s3.c1.re;
    (*r).s3.c1.im = fact*(*r).s3.c1.im + (*s).s3.c1.im;
    (*r).s3.c2.re = fact*(*r).s3.c2.re + (*s).s3.c2.re;
    (*r).s3.c2.im = fact*(*r).s3.c2.im + (*s).s3.c2.im;
  }
}
#endif

/* $Id$ */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#ifdef MPI
#include <mpi.h>
#endif
#include "sse.h"
#include "su3.h"
#include "assign_add_mul_r.h"


#if defined SSE2
/*  k input, l output */
void assign_add_mul_r(spinor * const P, spinor * const Q, const double c, const int N) {
  int ix;
  su3_vector *s,*r;
  __asm__ __volatile__ ("movsd %0, %%xmm7 \n\t"
			"unpcklpd %%xmm7, %%xmm7"
			:
			:
			"m" (c));
  s=(su3_vector*)P;
  r=(su3_vector*)Q;
  for (ix = 0;ix < 4*N; ix++) {
    _sse_load_up(*r);
    __asm__ __volatile__ ("mulpd %%xmm7, %%xmm3 \n\t"
			  "mulpd %%xmm7, %%xmm4 \n\t"
			  "mulpd %%xmm7, %%xmm5"
			  :
			  :);
    _sse_load(*s);
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
/*  k input, l output */
void assign_add_mul_r(spinor * const P, spinor * const Q, const double c, const int N) {
  int ix;
  static double fact;
  spinor *r,*s;

  fact=c;

  /* Change due to even-odd preconditioning : VOLUME   to VOLUME/2 */   
  for (ix = 0; ix < N; ix++) {
    r=P+ix;      
    s=Q+ix;
    
    (*r).c1.c1.re+=fact*(*s).s0.c1.re;
    (*r).c1.c1.im+=fact*(*s).s0.c1.im;
    (*r).c1.c2.re+=fact*(*s).s0.c2.re;
    (*r).c1.c2.im+=fact*(*s).s0.c2.im;
    (*r).c1.c3.re+=fact*(*s).s0.c3.re;
    (*r).c1.c3.im+=fact*(*s).s0.c3.im;
    
    (*r).s1.c1.re+=fact*(*s).s1.c1.re;
    (*r).s1.c1.im+=fact*(*s).s1.c1.im;
    (*r).s1.c2.re+=fact*(*s).s1.c2.re;
    (*r).s1.c2.im+=fact*(*s).s1.c2.im;
    (*r).s1.c3.re+=fact*(*s).s1.c3.re;
    (*r).s1.c3.im+=fact*(*s).s1.c3.im;         
    
    (*r).s2.c1.re+=fact*(*s).s2.c1.re;
    (*r).s2.c1.im+=fact*(*s).s2.c1.im;
    (*r).s2.c2.re+=fact*(*s).s2.c2.re;
    (*r).s2.c2.im+=fact*(*s).s2.c2.im;
    (*r).s2.c3.re+=fact*(*s).s2.c3.re;
    (*r).s2.c3.im+=fact*(*s).s2.c3.im;         
    
    (*r).s3.c1.re+=fact*(*s).s3.c1.re;
    (*r).s3.c1.im+=fact*(*s).s3.c1.im;
    (*r).s3.c2.re+=fact*(*s).s3.c2.re;
    (*r).s3.c2.im+=fact*(*s).s3.c2.im;
    (*r).s3.c3.re+=fact*(*s).s3.c3.re;
    (*r).s3.c3.im+=fact*(*s).s3.c3.im;
  }
}
#endif

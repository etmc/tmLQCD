/* $Id$ */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#ifdef MPI
#include <mpi.h>
#endif
#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include "sse.h"
#include "su3.h"
#include "su3adj.h"
#include "assign_add_mul_r_add_mul.h"

#if defined SSE2
void assign_add_mul_r_add_mul(spinor * const R, spinor * const S, spinor * const U,
			      const double c1,const double c2, const int N) {

  int ix;
  su3_vector *s,*r,*t;
  r=&R[0].s0;
  s=&S[0].s0;
  t=&U[0].s0;
  __asm__ __volatile__ ("movsd %0, %%xmm6 \n\t"
			"unpcklpd %%xmm6, %%xmm6"
			:
			:
			"m" (c1));
  __asm__ __volatile__ ("movsd %0, %%xmm7 \n\t"
			"unpcklpd %%xmm7, %%xmm7"
			:
			:
			"m" (c2));
  
  for (ix = 0; ix < 4*N; ix++) {
    _sse_load_up(*s);
    __asm__ __volatile__ ("mulpd %%xmm6, %%xmm3 \n\t"
			  "mulpd %%xmm6, %%xmm4 \n\t"
			  "mulpd %%xmm6, %%xmm5"
			  :
			  :);
    _sse_load(*r);
    __asm__ __volatile__ ("addpd %%xmm3, %%xmm0 \n\t"
			  "addpd %%xmm4, %%xmm1 \n\t"
			  "addpd %%xmm5, %%xmm2"
			  :
			  :);
    _sse_load_up(*t);
    __asm__ __volatile__ ("mulpd %%xmm7, %%xmm3 \n\t"
			  "mulpd %%xmm7, %%xmm4 \n\t"
			  "mulpd %%xmm7, %%xmm5"
			  :
			  :);
    __asm__ __volatile__ ("addpd %%xmm3, %%xmm0 \n\t"
			  "addpd %%xmm4, %%xmm1 \n\t"
			  "addpd %%xmm5, %%xmm2"
			  :
			  :);
    _sse_store(*r);
    r++; s++; t++;
  }
}
#else
/* j, k input, l output */
void assign_add_mul_r_add_mul(spinor * const R, spinor * const S, spinor * const U,
			      const double c1,const double c2, const int N) {

  int ix;
   spinor *r,*s,*t;

   for (ix = 0; ix < N; ix++) {
     r=R+ix;      
     s=S+ix;
     t=U+ix;
     
     (*r).s0.c0.re+=c1*(*s).s0.c0.re+c2*(*t).s0.c0.re;
     (*r).s0.c0.im+=c1*(*s).s0.c0.im+c2*(*t).s0.c0.im;
     (*r).s0.c1.re+=c1*(*s).s0.c1.re+c2*(*t).s0.c1.re;
     (*r).s0.c1.im+=c1*(*s).s0.c1.im+c2*(*t).s0.c1.im;
     (*r).s0.c2.re+=c1*(*s).s0.c2.re+c2*(*t).s0.c2.re;
     (*r).s0.c2.im+=c1*(*s).s0.c2.im+c2*(*t).s0.c2.im;
     
     (*r).s1.c0.re+=c1*(*s).s1.c0.re+c2*(*t).s1.c0.re;
     (*r).s1.c0.im+=c1*(*s).s1.c0.im+c2*(*t).s1.c0.im;
     (*r).s1.c1.re+=c1*(*s).s1.c1.re+c2*(*t).s1.c1.re;
     (*r).s1.c1.im+=c1*(*s).s1.c1.im+c2*(*t).s1.c1.im;
     (*r).s1.c2.re+=c1*(*s).s1.c2.re+c2*(*t).s1.c2.re;
     (*r).s1.c2.im+=c1*(*s).s1.c2.im+c2*(*t).s1.c2.im;         
     
     (*r).s2.c0.re+=c1*(*s).s2.c0.re+c2*(*t).s2.c0.re;
     (*r).s2.c0.im+=c1*(*s).s2.c0.im+c2*(*t).s2.c0.im;
     (*r).s2.c1.re+=c1*(*s).s2.c1.re+c2*(*t).s2.c1.re;
     (*r).s2.c1.im+=c1*(*s).s2.c1.im+c2*(*t).s2.c1.im;
     (*r).s2.c2.re+=c1*(*s).s2.c2.re+c2*(*t).s2.c2.re;
     (*r).s2.c2.im+=c1*(*s).s2.c2.im+c2*(*t).s2.c2.im;         
     
     (*r).s3.c0.re+=c1*(*s).s3.c0.re+c2*(*t).s3.c0.re;
     (*r).s3.c0.im+=c1*(*s).s3.c0.im+c2*(*t).s3.c0.im;
     (*r).s3.c1.re+=c1*(*s).s3.c1.re+c2*(*t).s3.c1.re;
     (*r).s3.c1.im+=c1*(*s).s3.c1.im+c2*(*t).s3.c1.im;
     (*r).s3.c2.re+=c1*(*s).s3.c2.re+c2*(*t).s3.c2.re;
     (*r).s3.c2.im+=c1*(*s).s3.c2.im+c2*(*t).s3.c2.im;
   }
}
#endif

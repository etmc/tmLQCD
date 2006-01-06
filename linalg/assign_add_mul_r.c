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
#include "global.h"

#ifdef _STD_C99_COMPLEX_CHECKED
#include <complex.h>
#endif

#ifdef apenext
#include <topology.h>
#include <queue.h>
#endif

#if defined SSE2

/*   (*P) = (*P) + c(*Q)        c is a complex constant   */

void assign_add_mul_r(spinor * const P, spinor * const Q, const double c, const int N) {
  int ix;
  su3_vector *s,*r;
  __asm__ __volatile__ ("movsd %0, %%xmm7 \n\t"
			"unpcklpd %%xmm7, %%xmm7"
			:
			:
			"m" (c));
  s=&P[0].s0;
  r=&Q[0].s0;
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

#endif

#if ((!defined _STD_C99_COMPLEX_CHECKED) && (!defined apenext))

/*   (*P) = (*P) + c(*Q)        c is a complex constant   */

void assign_add_mul_r(spinor * const P, spinor * const Q, const double c, const int N) {
  int ix;
  static double fact;
  spinor *r,*s;

  fact=c;

  /* Change due to even-odd preconditioning : VOLUME   to VOLUME/2 */   
  for (ix = 0; ix < N; ix++) {
    r=P+ix;      
    s=Q+ix;
    
    (*r).s0.c0.re+=fact*(*s).s0.c0.re;
    (*r).s0.c0.im+=fact*(*s).s0.c0.im;
    (*r).s0.c1.re+=fact*(*s).s0.c1.re;
    (*r).s0.c1.im+=fact*(*s).s0.c1.im;
    (*r).s0.c2.re+=fact*(*s).s0.c2.re;
    (*r).s0.c2.im+=fact*(*s).s0.c2.im;
    
    (*r).s1.c0.re+=fact*(*s).s1.c0.re;
    (*r).s1.c0.im+=fact*(*s).s1.c0.im;
    (*r).s1.c1.re+=fact*(*s).s1.c1.re;
    (*r).s1.c1.im+=fact*(*s).s1.c1.im;
    (*r).s1.c2.re+=fact*(*s).s1.c2.re;
    (*r).s1.c2.im+=fact*(*s).s1.c2.im;         
    
    (*r).s2.c0.re+=fact*(*s).s2.c0.re;
    (*r).s2.c0.im+=fact*(*s).s2.c0.im;
    (*r).s2.c1.re+=fact*(*s).s2.c1.re;
    (*r).s2.c1.im+=fact*(*s).s2.c1.im;
    (*r).s2.c2.re+=fact*(*s).s2.c2.re;
    (*r).s2.c2.im+=fact*(*s).s2.c2.im;         
    
    (*r).s3.c0.re+=fact*(*s).s3.c0.re;
    (*r).s3.c0.im+=fact*(*s).s3.c0.im;
    (*r).s3.c1.re+=fact*(*s).s3.c1.re;
    (*r).s3.c1.im+=fact*(*s).s3.c1.im;
    (*r).s3.c2.re+=fact*(*s).s3.c2.re;
    (*r).s3.c2.im+=fact*(*s).s3.c2.im;
  }
}
#endif

#if ((defined _STD_C99_COMPLEX_CHECKED) && (!defined apenext))

/*   (*P) = (*P) + c(*Q)        c is a complex constant   */

void assign_add_mul_r(spinor * const P, spinor * const Q, const double c, const int N) {
  register int ix=0;
  register double fact;
  register spinor *rPointer,*sPointer;

  fact=c;

  rPointer = P;
  sPointer = Q;

  /* Change due to even-odd preconditioning : VOLUME   to VOLUME/2 */   
  do {
    register spinor s, r;
    ix+=1;
    
    s = *(sPointer);
    r = *(rPointer);

    r.s0.c0+=fact*s.s0.c0;
    r.s0.c1+=fact*s.s0.c1;
    r.s0.c2+=fact*s.s0.c2;

    r.s1.c0+=fact*s.s1.c0;
    r.s1.c1+=fact*s.s1.c1;
    r.s1.c2+=fact*s.s1.c2;
    
    r.s2.c0+=fact*s.s2.c0;
    r.s2.c1+=fact*s.s2.c1;
    r.s2.c2+=fact*s.s2.c2;

    r.s3.c0+=fact*s.s3.c0;
    r.s3.c1+=fact*s.s3.c1;
    r.s3.c2+=fact*s.s3.c2;

    sPointer+=1;
    rPointer+=1;

  } while (ix<N);
}
#endif

#ifdef apenext

#define NOWHERE_COND(condition) ((condition) ? 0x0 : NOWHERE )

/*   (*P) = (*P) + c(*Q)        c is a complex constant   */

void assign_add_mul_r(spinor * const P, spinor * const Q, const double c, const int N) {
  register int ix=N;
  register double fact;
  register spinor *rPointer,*sPointer;

  fact=c;

  rPointer = P;
  sPointer = Q;

  prefetch(*(rPointer));
  prefetch(*(sPointer));
  {
#pragma cache
  /* Change due to even-odd preconditioning : VOLUME   to VOLUME/2 */   
    do {
      register spinor s, r;
      ix--;
    
      rPointer++;
      sPointer++;

      fetch(r);
      fetch(s);      

      prefetch (*(rPointer+NOWHERE_COND(ix)));
      prefetch (*(sPointer+NOWHERE_COND(ix)));

      r.s0.c0+=fact*s.s0.c0;
      r.s0.c1+=fact*s.s0.c1;
      r.s0.c2+=fact*s.s0.c2;

      r.s1.c0+=fact*s.s1.c0;
      r.s1.c1+=fact*s.s1.c1;
      r.s1.c2+=fact*s.s1.c2;
    
      r.s2.c0+=fact*s.s2.c0;
      r.s2.c1+=fact*s.s2.c1;
      r.s2.c2+=fact*s.s2.c2;

      r.s3.c0+=fact*s.s3.c0;
      r.s3.c1+=fact*s.s3.c1;
      r.s3.c2+=fact*s.s3.c2;

      *(rPointer-1) = r;

    } while (ix>0);
  }
}
#endif

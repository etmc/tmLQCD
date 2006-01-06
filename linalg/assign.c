/*******************************************************************************
 * $Id$
 *
 * File assign.c
 *
 *   void assign(spinor * const R, spinor * const S)
 *     Assign (*R) = (*S)
 *
 *******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "su3.h"
#include "assign.h"
#include "global.h"

#ifdef _STD_C99_COMPLEX_CHECKED
#include <complex.h>
#endif

#ifdef apenext
#include <topology.h>
#include <queue.h>
#endif

/* S input, R output */

#if ((!defined _STD_C99_COMPLEX_CHECKED) && (!defined apenext))

void assign(spinor * const R, spinor * const S, const int N){
  int ix;
  spinor *r,*s;
  
  for (ix = 0; ix < N; ix++){
    r=(spinor *) R + ix;
    s=(spinor *) S + ix;
    
    (*r).s0.c0.re = (*s).s0.c0.re;
    (*r).s0.c0.im = (*s).s0.c0.im;
    (*r).s0.c1.re = (*s).s0.c1.re;
    (*r).s0.c1.im = (*s).s0.c1.im;
    (*r).s0.c2.re = (*s).s0.c2.re;
    (*r).s0.c2.im = (*s).s0.c2.im;
    
    (*r).s1.c0.re = (*s).s1.c0.re;
    (*r).s1.c0.im = (*s).s1.c0.im;
    (*r).s1.c1.re = (*s).s1.c1.re;
    (*r).s1.c1.im = (*s).s1.c1.im;
    (*r).s1.c2.re = (*s).s1.c2.re;
    (*r).s1.c2.im = (*s).s1.c2.im;         
    
    (*r).s2.c0.re = (*s).s2.c0.re;
    (*r).s2.c0.im = (*s).s2.c0.im;
    (*r).s2.c1.re = (*s).s2.c1.re;
    (*r).s2.c1.im = (*s).s2.c1.im;
    (*r).s2.c2.re = (*s).s2.c2.re;
    (*r).s2.c2.im = (*s).s2.c2.im;         
    
    (*r).s3.c0.re = (*s).s3.c0.re;
    (*r).s3.c0.im = (*s).s3.c0.im;
    (*r).s3.c1.re = (*s).s3.c1.re;
    (*r).s3.c1.im = (*s).s3.c1.im;
    (*r).s3.c2.re = (*s).s3.c2.re;
    (*r).s3.c2.im = (*s).s3.c2.im;
  }
}
#endif

#if ((defined _STD_C99_COMPLEX_CHECKED) && (!defined apenext))

void assign(spinor * const R, spinor * const S, const int N){
  register int ix=0;
  register spinor *rPointer,*sPointer;

  rPointer = R;
  sPointer = S;

  do {
    register spinor s, r;
    ix+=1;

    s = *(sPointer);

    r.s0.c0 = s.s0.c0;
    r.s0.c1 = s.s0.c1;
    r.s0.c2 = s.s0.c2;

    r.s1.c0 = s.s1.c0;
    r.s1.c1 = s.s1.c1;
    r.s1.c2 = s.s1.c2;

    r.s2.c0 = s.s2.c0;
    r.s2.c1 = s.s2.c1;
    r.s2.c2 = s.s2.c2;
    
    r.s3.c0 = s.s3.c0;
    r.s3.c1 = s.s3.c1;
    r.s3.c2 = s.s3.c2;

    *(rPointer) = r;

    rPointer+=1;
    sPointer+=1;

 } while (ix<N);
}
#endif


#ifdef apenext

#define NOWHERE_COND(condition) ((condition) ? 0x0 : NOWHERE )

void assign(spinor * const R, spinor * const S, const int N){
  register int ix=N;
  register spinor *rPointer,*sPointer;

  rPointer = R;
  sPointer = S;

  prefetch (*(sPointer));
  {
#pragma cache

    do {
      register spinor s, r;
      ix-=1;

      sPointer+=1;

      fetch(s);

      prefetch (*(sPointer+NOWHERE_COND(ix)));
      
      r.s0.c0 = s.s0.c0;
      r.s0.c1 = s.s0.c1;
      r.s0.c2 = s.s0.c2;

      r.s1.c0 = s.s1.c0;
      r.s1.c1 = s.s1.c1;
      r.s1.c2 = s.s1.c2;

      r.s2.c0 = s.s2.c0;
      r.s2.c1 = s.s2.c1;
      r.s2.c2 = s.s2.c2;
    
      r.s3.c0 = s.s3.c0;
      r.s3.c1 = s.s3.c1;
      r.s3.c2 = s.s3.c2;

      *(rPointer) = r;

      rPointer+=1;

    } while (ix>0);
  }
}
#endif

/*******************************************************************************
 * $Id$
 *
 *   void diff(spinor * const Q,spinor * const R,spinor * const S)
 *     Makes the difference (*Q) = (*R) - (*S)
 *******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "su3.h"
#include "global.h"
#include "diff.h"

#ifdef _STD_C99_COMPLEX_CHECKED
#include <complex.h>
#endif

#ifdef apenext
#include <topology.h>
#include <queue.h>
#endif

#if ((!defined _STD_C99_COMPLEX_CHECKED) && (!defined apenext))

void diff(spinor * const Q,spinor * const R,spinor * const S, const int N)
{
   int ix;
   spinor *q,*r,*s;


/* Change due to even-odd preconditioning : VOLUME   to VOLUME/2 */   
   for (ix = 0; ix < N; ix++) {
     q=(spinor *) Q + ix;
     r=(spinor *) R + ix;
     s=(spinor *) S + ix;
     
     (*q).s0.c0.re=(*r).s0.c0.re-(*s).s0.c0.re;
     (*q).s0.c0.im=(*r).s0.c0.im-(*s).s0.c0.im;
     (*q).s0.c1.re=(*r).s0.c1.re-(*s).s0.c1.re;
     (*q).s0.c1.im=(*r).s0.c1.im-(*s).s0.c1.im;
     (*q).s0.c2.re=(*r).s0.c2.re-(*s).s0.c2.re;
     (*q).s0.c2.im=(*r).s0.c2.im-(*s).s0.c2.im;
     
     (*q).s1.c0.re=(*r).s1.c0.re-(*s).s1.c0.re;
     (*q).s1.c0.im=(*r).s1.c0.im-(*s).s1.c0.im;
     (*q).s1.c1.re=(*r).s1.c1.re-(*s).s1.c1.re;
     (*q).s1.c1.im=(*r).s1.c1.im-(*s).s1.c1.im;
     (*q).s1.c2.re=(*r).s1.c2.re-(*s).s1.c2.re;
     (*q).s1.c2.im=(*r).s1.c2.im-(*s).s1.c2.im;         
     
     (*q).s2.c0.re=(*r).s2.c0.re-(*s).s2.c0.re;
     (*q).s2.c0.im=(*r).s2.c0.im-(*s).s2.c0.im;
     (*q).s2.c1.re=(*r).s2.c1.re-(*s).s2.c1.re;
     (*q).s2.c1.im=(*r).s2.c1.im-(*s).s2.c1.im;
     (*q).s2.c2.re=(*r).s2.c2.re-(*s).s2.c2.re;
     (*q).s2.c2.im=(*r).s2.c2.im-(*s).s2.c2.im;         
     
     (*q).s3.c0.re=(*r).s3.c0.re-(*s).s3.c0.re;
     (*q).s3.c0.im=(*r).s3.c0.im-(*s).s3.c0.im;
     (*q).s3.c1.re=(*r).s3.c1.re-(*s).s3.c1.re;
     (*q).s3.c1.im=(*r).s3.c1.im-(*s).s3.c1.im;
     (*q).s3.c2.re=(*r).s3.c2.re-(*s).s3.c2.re;
     (*q).s3.c2.im=(*r).s3.c2.im-(*s).s3.c2.im;
   }
}

#endif

#if ((defined _STD_C99_COMPLEX_CHECKED) && (!defined apenext))

void diff(spinor * const Q,spinor * const R,spinor * const S, const int N){
  register int ix=0;
  register spinor *qPointer,*rPointer,*sPointer;

  qPointer = Q;
  rPointer = R;
  sPointer = S;

  do {
    register spinor q,r,s;
    ix+=1;
/* Change due to even-odd preconditioning : VOLUME   to VOLUME/2 */   

    q = *(qPointer);
    r = *(rPointer);
    s = *(sPointer);

    q.s0.c0 = r.s0.c0 - s.s0.c0;
    q.s0.c1 = r.s0.c1 - s.s0.c1;
    q.s0.c2 = r.s0.c2 - s.s0.c2;

    q.s1.c0 = r.s1.c0 - s.s1.c0;
    q.s1.c1 = r.s1.c1 - s.s1.c1;
    q.s1.c2 = r.s1.c2 - s.s1.c2;

    q.s2.c0 = r.s2.c0 - s.s2.c0;
    q.s2.c1 = r.s2.c1 - s.s2.c1;
    q.s2.c2 = r.s2.c2 - s.s2.c2;

    q.s3.c0 = r.s3.c0 - s.s3.c0;
    q.s3.c1 = r.s3.c1 - s.s3.c1;
    q.s3.c2 = r.s3.c2 - s.s3.c2;

    qPointer+=1;
    rPointer+=1;
    sPointer+=1;
    
  } while (ix<N);
}

#endif

#ifdef apenext

#define NOWHERE_COND(condition) ((condition) ? 0x0 : NOWHERE )

void diff(spinor * const Q,spinor * const R,spinor * const S, const int N){
  register int ix=N;
  register spinor *qPointer,*rPointer,*sPointer;

  qPointer = Q;
  rPointer = R;
  sPointer = S;

  prefetch (*(rPointer));
  prefetch (*(sPointer));
  
  {
#pragma cache

    do {
      register spinor q,r,s;
      ix-=1;

      rPointer+=1;
      sPointer+=1;

      fetch(r);
      fetch(s);

      prefetch (*(rPointer+NOWHERE_COND(ix)));
      prefetch (*(sPointer+NOWHERE_COND(ix)));

      q.s0.c0 = r.s0.c0 - s.s0.c0;
      q.s0.c1 = r.s0.c1 - s.s0.c1;
      q.s0.c2 = r.s0.c2 - s.s0.c2;

      q.s1.c0 = r.s1.c0 - s.s1.c0;
      q.s1.c1 = r.s1.c1 - s.s1.c1;
      q.s1.c2 = r.s1.c2 - s.s1.c2;

      q.s2.c0 = r.s2.c0 - s.s2.c0;
      q.s2.c1 = r.s2.c1 - s.s2.c1;
      q.s2.c2 = r.s2.c2 - s.s2.c2;

      q.s3.c0 = r.s3.c0 - s.s3.c0;
      q.s3.c1 = r.s3.c1 - s.s3.c1;
      q.s3.c2 = r.s3.c2 - s.s3.c2;

      *(qPointer) = q;

      qPointer+=1;
    
    } while (ix>0);
  }
}

#endif

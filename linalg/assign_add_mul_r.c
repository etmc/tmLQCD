/* $Id$ */

#ifdef _STD_C99_COMPLEX_CHECKED
# include <complex.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#ifdef MPI
# include <mpi.h>
#endif
#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#ifdef apenext
#include <topology.h>
#include <queue.h>
#endif
#include "global.h"
#include "sse.h"
#include "su3.h"
#include "assign_add_mul_r.h"


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

#elif ((defined BGL) && (defined XLC))

#  include"bgl.h"

void assign_add_mul_r(spinor * const R, spinor * const S, const double c, const int N) {
  int ix = 1;
  double *s ALIGN;
  double *sp ALIGN;
  double *r ALIGN;
  double *rp ALIGN;
  double _Complex x00, x01, x02, x03, x04, x05, x06, x07, 
    x08, x09, x10, x11;
  double _Complex y00, y01, y02, y03, y04, y05, y06, y07, 
    y08, y09, y10, y11;
  double _Complex a;

#pragma disjoint(*S, *R)
  a = __cmplx(c, c);
  __alignx(16, S);
  __alignx(16, R);
  s = (double*) S;
  r = (double*) R;
  rp = r + 24;
  sp = s + 24;
  _prefetch_spinor(rp);
  _prefetch_spinor(sp);
  x00 = __lfpd(s);    
  x01 = __lfpd(s+2);  
  x02 = __lfpd(s+4);  
  x03 = __lfpd(s+6);  
  x04 = __lfpd(s+8);  
  x05 = __lfpd(s+10); 
  x06 = __lfpd(s+12); 
  x07 = __lfpd(s+14); 
  x08 = __lfpd(s+16); 
  x09 = __lfpd(s+18); 
  x10 = __lfpd(s+20); 
  x11 = __lfpd(s+22); 
  y00 = __lfpd(r);   
  y01 = __lfpd(r+2); 
  y02 = __lfpd(r+4); 
  y03 = __lfpd(r+6); 
  y04 = __lfpd(r+8); 
  y05 = __lfpd(r+10);
  y06 = __lfpd(r+12);
  y07 = __lfpd(r+14);
  y08 = __lfpd(r+16);
  y09 = __lfpd(r+18);
  y10 = __lfpd(r+20);
  y11 = __lfpd(r+22);

  y00 = __fpmadd(y00, x00, a);
  y01 = __fpmadd(y01, x01, a);
  y02 = __fpmadd(y02, x02, a);
  y03 = __fpmadd(y03, x03, a);
  y04 = __fpmadd(y04, x04, a);
  y05 = __fpmadd(y05, x05, a);
  y06 = __fpmadd(y06, x06, a);
  y07 = __fpmadd(y07, x07, a);
  y08 = __fpmadd(y08, x08, a);
  y09 = __fpmadd(y09, x09, a);
  y10 = __fpmadd(y10, x10, a);
  y11 = __fpmadd(y11, x11, a);
  __stfpd(r, y00);
  __stfpd(r+2, y01);
  __stfpd(r+4, y02);
  __stfpd(r+6, y03);
  __stfpd(r+8, y04);
  __stfpd(r+10, y05);
  __stfpd(r+12, y06);
  __stfpd(r+14, y07);
  __stfpd(r+16, y08);
  __stfpd(r+18, y09);
  __stfpd(r+20, y10);
  __stfpd(r+22, y11);
  s = sp;
  r = rp;

#pragma unroll(12)
  for(ix = 1; ix < N-1; ix++) {
    rp += 24;
    sp += 24;
    _prefetch_spinor(rp);
    _prefetch_spinor(sp);
    x00 = __lfpd(s);    
    x01 = __lfpd(s+2);  
    x02 = __lfpd(s+4);  
    x03 = __lfpd(s+6);  
    x04 = __lfpd(s+8);  
    x05 = __lfpd(s+10); 
    x06 = __lfpd(s+12); 
    x07 = __lfpd(s+14); 
    x08 = __lfpd(s+16); 
    x09 = __lfpd(s+18); 
    x10 = __lfpd(s+20); 
    x11 = __lfpd(s+22); 
    y00 = __lfpd(r);   
    y01 = __lfpd(r+2); 
    y02 = __lfpd(r+4); 
    y03 = __lfpd(r+6); 
    y04 = __lfpd(r+8); 
    y05 = __lfpd(r+10);
    y06 = __lfpd(r+12);
    y07 = __lfpd(r+14);
    y08 = __lfpd(r+16);
    y09 = __lfpd(r+18);
    y10 = __lfpd(r+20);
    y11 = __lfpd(r+22);

    y00 = __fpmadd(y00, x00, a);
    y01 = __fpmadd(y01, x01, a);
    y02 = __fpmadd(y02, x02, a);
    y03 = __fpmadd(y03, x03, a);
    y04 = __fpmadd(y04, x04, a);
    y05 = __fpmadd(y05, x05, a);
    y06 = __fpmadd(y06, x06, a);
    y07 = __fpmadd(y07, x07, a);
    y08 = __fpmadd(y08, x08, a);
    y09 = __fpmadd(y09, x09, a);
    y10 = __fpmadd(y10, x10, a);
    y11 = __fpmadd(y11, x11, a);
    __stfpd(r, y00);
    __stfpd(r+2, y01);
    __stfpd(r+4, y02);
    __stfpd(r+6, y03);
    __stfpd(r+8, y04);
    __stfpd(r+10, y05);
    __stfpd(r+12, y06);
    __stfpd(r+14, y07);
    __stfpd(r+16, y08);
    __stfpd(r+18, y09);
    __stfpd(r+20, y10);
    __stfpd(r+22, y11);
    s = sp;
    r = rp;

  }
  x00 = __lfpd(s);    
  x01 = __lfpd(s+2);  
  x02 = __lfpd(s+4);  
  x03 = __lfpd(s+6);  
  x04 = __lfpd(s+8);  
  x05 = __lfpd(s+10); 
  x06 = __lfpd(s+12); 
  x07 = __lfpd(s+14); 
  x08 = __lfpd(s+16); 
  x09 = __lfpd(s+18); 
  x10 = __lfpd(s+20); 
  x11 = __lfpd(s+22); 
  y00 = __lfpd(r);   
  y01 = __lfpd(r+2); 
  y02 = __lfpd(r+4); 
  y03 = __lfpd(r+6); 
  y04 = __lfpd(r+8); 
  y05 = __lfpd(r+10);
  y06 = __lfpd(r+12);
  y07 = __lfpd(r+14);
  y08 = __lfpd(r+16);
  y09 = __lfpd(r+18);
  y10 = __lfpd(r+20);
  y11 = __lfpd(r+22);

  y00 = __fpmadd(y00, x00, a);
  y01 = __fpmadd(y01, x01, a);
  y02 = __fpmadd(y02, x02, a);
  y03 = __fpmadd(y03, x03, a);
  y04 = __fpmadd(y04, x04, a);
  y05 = __fpmadd(y05, x05, a);
  y06 = __fpmadd(y06, x06, a);
  y07 = __fpmadd(y07, x07, a);
  y08 = __fpmadd(y08, x08, a);
  y09 = __fpmadd(y09, x09, a);
  y10 = __fpmadd(y10, x10, a);
  y11 = __fpmadd(y11, x11, a);
  __stfpd(r, y00);
  __stfpd(r+2, y01);
  __stfpd(r+4, y02);
  __stfpd(r+6, y03);
  __stfpd(r+8, y04);
  __stfpd(r+10, y05);
  __stfpd(r+12, y06);
  __stfpd(r+14, y07);
  __stfpd(r+16, y08);
  __stfpd(r+18, y09);
  __stfpd(r+20, y10);
  __stfpd(r+22, y11);

  return;
}

#elif ((!defined _STD_C99_COMPLEX_CHECKED) && (!defined apenext))

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
#elif ((defined _STD_C99_COMPLEX_CHECKED) && (!defined apenext))

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
#elif defined apenext

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

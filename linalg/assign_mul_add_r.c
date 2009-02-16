/***********************************************************************
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008 Carsten Urbach
 *
 * This file is part of tmLQCD.
 *
 * tmLQCD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * tmLQCD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with tmLQCD.  If not, see <http://www.gnu.org/licenses/>.
 ***********************************************************************/
/* $Id$ */

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
# ifdef _STD_C99_COMPLEX_CHECKED
#include <complex.h>
#endif
#ifdef apenext
# include <topology.h>
# include <queue.h>
#endif
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

#elif ((defined BGL) && (defined XLC))

#  include"bgl.h"

void assign_mul_add_r(spinor * const R, const double c, spinor * const S, const int N) {
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
  x00 = __lfpd(r);    
  x01 = __lfpd(r+2);  
  x02 = __lfpd(r+4);  
  x03 = __lfpd(r+6);  
  x04 = __lfpd(r+8);  
  x05 = __lfpd(r+10); 
  x06 = __lfpd(r+12); 
  x07 = __lfpd(r+14); 
  x08 = __lfpd(r+16); 
  x09 = __lfpd(r+18); 
  x10 = __lfpd(r+20); 
  x11 = __lfpd(r+22); 
  y00 = __lfpd(s);   
  y01 = __lfpd(s+2); 
  y02 = __lfpd(s+4); 
  y03 = __lfpd(s+6); 
  y04 = __lfpd(s+8); 
  y05 = __lfpd(s+10);
  y06 = __lfpd(s+12);
  y07 = __lfpd(s+14);
  y08 = __lfpd(s+16);
  y09 = __lfpd(s+18);
  y10 = __lfpd(s+20);
  y11 = __lfpd(s+22);

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
    x00 = __lfpd(r);    
    x01 = __lfpd(r+2);  
    x02 = __lfpd(r+4);  
    x03 = __lfpd(r+6);  
    x04 = __lfpd(r+8);  
    x05 = __lfpd(r+10); 
    x06 = __lfpd(r+12); 
    x07 = __lfpd(r+14); 
    x08 = __lfpd(r+16); 
    x09 = __lfpd(r+18); 
    x10 = __lfpd(r+20); 
    x11 = __lfpd(r+22); 
    y00 = __lfpd(s);   
    y01 = __lfpd(s+2); 
    y02 = __lfpd(s+4); 
    y03 = __lfpd(s+6); 
    y04 = __lfpd(s+8); 
    y05 = __lfpd(s+10);
    y06 = __lfpd(s+12);
    y07 = __lfpd(s+14);
    y08 = __lfpd(s+16);
    y09 = __lfpd(s+18);
    y10 = __lfpd(s+20);
    y11 = __lfpd(s+22);

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
  x00 = __lfpd(r);    
  x01 = __lfpd(r+2);  
  x02 = __lfpd(r+4);  
  x03 = __lfpd(r+6);  
  x04 = __lfpd(r+8);  
  x05 = __lfpd(r+10); 
  x06 = __lfpd(r+12); 
  x07 = __lfpd(r+14); 
  x08 = __lfpd(r+16); 
  x09 = __lfpd(r+18); 
  x10 = __lfpd(r+20); 
  x11 = __lfpd(r+22); 
  y00 = __lfpd(s);   
  y01 = __lfpd(s+2); 
  y02 = __lfpd(s+4); 
  y03 = __lfpd(s+6); 
  y04 = __lfpd(s+8); 
  y05 = __lfpd(s+10);
  y06 = __lfpd(s+12);
  y07 = __lfpd(s+14);
  y08 = __lfpd(s+16);
  y09 = __lfpd(s+18);
  y10 = __lfpd(s+20);
  y11 = __lfpd(s+22);

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

/* R inoutput , c,S input*/
/*   (*R) = c*(*R) + (*S)        c is a real constant   */

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

/*       printf("%1.16e %1.16e\n",(*r).s0.c0.re,(*r).s0.c0.im); */
/*       printf("%1.16e %1.16e\n",(*r).s0.c1.re,(*r).s0.c1.im); */
/*       printf("%1.16e %1.16e\n",(*r).s0.c2.re,(*r).s0.c2.im); */

/*       printf("%1.16e %1.16e\n",(*r).s1.c0.re,(*r).s1.c0.im); */
/*       printf("%1.16e %1.16e\n",(*r).s1.c1.re,(*r).s1.c1.im); */
/*       printf("%1.16e %1.16e\n",(*r).s1.c2.re,(*r).s1.c2.im); */

/*       printf("%1.16e %1.16e\n",(*r).s2.c0.re,(*r).s2.c0.im); */
/*       printf("%1.16e %1.16e\n",(*r).s2.c1.re,(*r).s2.c1.im); */
/*       printf("%1.16e %1.16e\n",(*r).s2.c2.re,(*r).s2.c2.im); */

/*       printf("%1.16e %1.16e\n",(*r).s3.c0.re,(*r).s3.c0.im); */
/*       printf("%1.16e %1.16e\n",(*r).s3.c1.re,(*r).s3.c1.im); */
/*       printf("%1.16e %1.16e\n",(*r).s3.c2.re,(*r).s3.c2.im); */
  }
}

#elif ((defined _STD_C99_COMPLEX_CHECKED) && (!defined apenext))

/*   (*R) = c*(*R) + (*S)        c is a real constant   */

void assign_mul_add_r(spinor * const R, const double c, spinor * const S, const int N) {

  register int ix=0;
  register double fact;
  register spinor *rPointer,*sPointer;

  fact=c;

  rPointer = R;
  sPointer = S;

  /* Change due to even-odd preconditioning : VOLUME   to VOLUME/2 */   
  do {
    register spinor s, r;
    ix+=1;
    
    s = *(sPointer);
    r = *(rPointer);

    r.s0.c0 =fact*r.s0.c0+s.s0.c0;
    r.s0.c1 =fact*r.s0.c1+s.s0.c1;
    r.s0.c2 =fact*r.s0.c2+s.s0.c2;
   	   	 	   
    r.s1.c0 =fact*r.s1.c0+s.s1.c0;
    r.s1.c1 =fact*r.s1.c1+s.s1.c1;
    r.s1.c2 =fact*r.s1.c2+s.s1.c2;
   	    	 	   
    r.s2.c0 =fact*r.s2.c0+s.s2.c0;
    r.s2.c1 =fact*r.s2.c1+s.s2.c1;
    r.s2.c2 =fact*r.s2.c2+s.s2.c2;
   	   	 	   
    r.s3.c0 =fact*r.s3.c0+s.s3.c0;
    r.s3.c1 =fact*r.s3.c1+s.s3.c1;
    r.s3.c2 =fact*r.s3.c2+s.s3.c2;

    *(rPointer) = r;

    sPointer+=1;
    rPointer+=1;

  } while (ix<N);
}

#elif defined apenext

#define NOWHERE_COND(condition) ((condition) ? 0x0 : NOWHERE ) 

/*   (*R) = c*(*R) + (*S)        c is a real constant   */

void assign_mul_add_r(spinor * const R, const double c, spinor * const S, const int N) {

  register int ix=N;
  register double fact;
  register spinor *rPointer,*sPointer;

  fact=c;

  rPointer = R;
  sPointer = S;

  prefetch(*(rPointer));
  prefetch(*(sPointer));

  {
#pragma cache

    do {
      register spinor s, r;
      register spinor *aux;
      ix--;
    
      rPointer++;
      sPointer++;

      fetch(r);
      fetch(s);


      prefetch (*(rPointer+NOWHERE_COND(ix)));
      prefetch (*(sPointer+NOWHERE_COND(ix)));


      r.s0.c0 = fact*r.s0.c0+s.s0.c0;
      r.s0.c1 = fact*r.s0.c1+s.s0.c1;
      r.s0.c2 = fact*r.s0.c2+s.s0.c2;
   	        	 	   
      r.s1.c0 = fact*r.s1.c0+s.s1.c0;
      r.s1.c1 = fact*r.s1.c1+s.s1.c1;
      r.s1.c2 = fact*r.s1.c2+s.s1.c2;
   	        	 	   
      r.s2.c0 = fact*r.s2.c0+s.s2.c0;
      r.s2.c1 = fact*r.s2.c1+s.s2.c1;
      r.s2.c2 = fact*r.s2.c2+s.s2.c2;
   	        	 	   
      r.s3.c0 = fact*r.s3.c0+s.s3.c0;
      r.s3.c1 = fact*r.s3.c1+s.s3.c1;
      r.s3.c2 = fact*r.s3.c2+s.s3.c2;

      *(rPointer-1) = r;

    } while (ix>0);
  }
}
#endif



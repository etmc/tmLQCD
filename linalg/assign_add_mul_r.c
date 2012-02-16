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

#include <complex.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#ifdef MPI
# include <mpi.h>
#endif
#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include "global.h"
#include "sse.h"
#include "su3.h"
#include "assign_add_mul_r.h"


#if defined SSE2

/*   (*P) = (*P) + c(*Q)        c is a complex constant   */

void assign_add_mul_r(spinor * const P, spinor * const Q, const double c, const int N)
{
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
    __stfpd(r+10, y05);r
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

#else

/*   (*P) = (*P) + c(*Q)        c is a complex constant   */

void assign_add_mul_r(spinor * const P, spinor * const Q, const double c, const int N)
{
  register spinor *p;
  register spinor *q;

  /* Change due to even-odd preconditioning : VOLUME   to VOLUME/2 */   
  for (int ix = 0; ix < N; ++ix)
  {
    p = P + ix;
    q = Q + ix;
    p->s0.c0 += c * q->s0.c0;
    p->s0.c1 += c * q->s0.c1;
    p->s0.c2 += c * q->s0.c2;

    p->s1.c0 += c * q->s1.c0;
    p->s1.c1 += c * q->s1.c1;
    p->s1.c2 += c * q->s1.c2;
    
    p->s2.c0 += c * q->s2.c0;
    p->s2.c1 += c * q->s2.c1;
    p->s2.c2 += c * q->s2.c2;

    p->s3.c0 += c * q->s3.c0;
    p->s3.c1 += c * q->s3.c1;
    p->s3.c2 += c * q->s3.c2;
  }
}
#endif

#ifdef WITHLAPH
void assign_add_mul_r_su3vect(su3_vector * const P, su3_vector * const Q, const double c, const int N)
{
  su3_vector *p,*q;

  for (int ix = 0; ix < N; ++ix) 
  {
    p = P + ix;      
    q = Q + ix;
    
    p->c0 += c * q->c0;
    p->c1 += c * q->c1;
    p->c2 += c * q->c2;
  }
}
#endif

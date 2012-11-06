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
/*******************************************************************************
 *
 *   void diff(spinor * const Q,spinor * const R,spinor * const S)
 *     Makes the difference (*Q) = (*R) - (*S)
 *******************************************************************************/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#ifdef OMP
# include <omp.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include "su3.h"
#include "diff.h"

#if ((defined BGL) && (defined XLC))

/***************************************
 *
 * diff with intrinsics
 *
 * Carsten.Urbach@liverpool.ac.uk
 *
 ***************************************/

#  include"bgl.h"

void diff(spinor * const Q,spinor * const R,spinor * const S, const int N)
{
  int ix = 1;
  double *s ALIGN;
  double *sp ALIGN;
  double *r ALIGN;
  double *rp ALIGN;
  double *q ALIGN;
  double _Complex x00, x01, x02, x03, x04, x05, x06, x07, 
    x08, x09, x10, x11;
  double _Complex y00, y01, y02, y03, y04, y05, y06, y07, 
    y08, y09, y10, y11;
#pragma disjoint(*R, *S)

  __alignx(16, Q);
  __alignx(16, R);
  __alignx(16, S);
  r = (double*) R;
  s = (double*) S;
  q = (double*) Q;
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

  __stfpd(q, __fpsub(x00, y00));
  __stfpd(q+2, __fpsub(x01, y01));
  __stfpd(q+4, __fpsub(x02, y02));
  __stfpd(q+6, __fpsub(x03, y03));
  __stfpd(q+8, __fpsub(x04, y04));
  __stfpd(q+10, __fpsub(x05, y05));
  __stfpd(q+12, __fpsub(x06, y06));
  __stfpd(q+14, __fpsub(x07, y07));
  __stfpd(q+16, __fpsub(x08, y08));
  __stfpd(q+18, __fpsub(x09, y09));
  __stfpd(q+20, __fpsub(x10, y10));
  __stfpd(q+22, __fpsub(x11, y11));
  s = sp;
  r = rp;
  q+=24;
#pragma unroll(12)
  for(ix = 1; ix < N-1; ix++) {
    rp+=24;
    sp+=24;
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
    
    __stfpd(q, __fpsub(x00, y00));
    __stfpd(q+2, __fpsub(x01, y01));
    __stfpd(q+4, __fpsub(x02, y02));
    __stfpd(q+6, __fpsub(x03, y03));
    __stfpd(q+8, __fpsub(x04, y04));
    __stfpd(q+10, __fpsub(x05, y05));
    __stfpd(q+12, __fpsub(x06, y06));
    __stfpd(q+14, __fpsub(x07, y07));
    __stfpd(q+16, __fpsub(x08, y08));
    __stfpd(q+18, __fpsub(x09, y09));
    __stfpd(q+20, __fpsub(x10, y10));
    __stfpd(q+22, __fpsub(x11, y11));
    s = sp;
    r = rp;
    q+=24;
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

  __stfpd(q, __fpsub(x00, y00));
  __stfpd(q+2, __fpsub(x01, y01));
  __stfpd(q+4, __fpsub(x02, y02));
  __stfpd(q+6, __fpsub(x03, y03));
  __stfpd(q+8, __fpsub(x04, y04));
  __stfpd(q+10, __fpsub(x05, y05));
  __stfpd(q+12, __fpsub(x06, y06));
  __stfpd(q+14, __fpsub(x07, y07));
  __stfpd(q+16, __fpsub(x08, y08));
  __stfpd(q+18, __fpsub(x09, y09));
  __stfpd(q+20, __fpsub(x10, y10));
  __stfpd(q+22, __fpsub(x11, y11));

  return;
}

#elif (defined BGQ && defined XLC)

void diff(spinor * const Q,const spinor * const R,const spinor * const S, const int N) {
#ifdef OMP
#pragma omp parallel
  {
#endif
  vector4double x0, x1, x2, x3, x4, x5, y0, y1, y2, y3, y4, y5;
  vector4double z0, z1, z2, z3, z4, z5;
  double *s, *r, *q;

  __alignx(32, s);
  __alignx(32, r);
  __alignx(32, q);
  __alignx(32, S);
  __alignx(32, R);

  __prefetch_by_load(S);
  __prefetch_by_load(R);
  __prefetch_by_load(Q);

#ifndef OMP
#pragma unroll(2)
#else
#pragma omp for
#endif
  for (int ix = 0; ix < N; ++ix) {
    s=(double*)((spinor *) S + ix);
    r=(double*)((spinor *) R + ix);
    q=(double*)((spinor *) Q + ix);
    __prefetch_by_load(S + ix + 1);
    __prefetch_by_load(R + ix + 1);
    __prefetch_by_stream(1, Q + ix + 1);
    x0 = vec_ld(0, r);
    x1 = vec_ld(0, r+4);
    x2 = vec_ld(0, r+8);
    x3 = vec_ld(0, r+12);
    x4 = vec_ld(0, r+16);
    x5 = vec_ld(0, r+20);
    y0 = vec_ld(0, s);
    y1 = vec_ld(0, s+4);
    y2 = vec_ld(0, s+8);
    y3 = vec_ld(0, s+12);
    y4 = vec_ld(0, s+16);
    y5 = vec_ld(0, s+20);
    z0 = vec_sub(x0, y0);
    z1 = vec_sub(x1, y1);
    z2 = vec_sub(x2, y2);
    z3 = vec_sub(x3, y3);
    z4 = vec_sub(x4, y4);
    z5 = vec_sub(x5, y5);
    vec_st(z0, 0, q);
    vec_st(z1, 0, q+4);
    vec_st(z2, 0, q+8);
    vec_st(z3, 0, q+12);
    vec_st(z4, 0, q+16);
    vec_st(z5, 0, q+20);
  }

#ifdef OMP
  } /* OpenMP parallel closing brace */
#endif

  return;
}

#else

void diff(spinor * const Q, const spinor * const R, const spinor * const S, const int N)
{
#ifdef OMP
#pragma omp parallel
  {
#endif

   spinor *q;
   const spinor *r,*s;

/* Change due to even-odd preconditioning : VOLUME   to VOLUME/2 */   
#ifdef OMP
#pragma omp for
#endif
   for (int ix = 0; ix < N; ix++)
   {
     q=(spinor *) Q + ix;
     r=(spinor *) R + ix;
     s=(spinor *) S + ix;
     
     q->s0.c0 = r->s0.c0 - s->s0.c0;
     q->s0.c1 = r->s0.c1 - s->s0.c1;
     q->s0.c2 = r->s0.c2 - s->s0.c2;

     q->s1.c0 = r->s1.c0 - s->s1.c0;
     q->s1.c1 = r->s1.c1 - s->s1.c1;
     q->s1.c2 = r->s1.c2 - s->s1.c2;

     q->s2.c0 = r->s2.c0 - s->s2.c0;
     q->s2.c1 = r->s2.c1 - s->s2.c1;
     q->s2.c2 = r->s2.c2 - s->s2.c2;

     q->s3.c0 = r->s3.c0 - s->s3.c0;
     q->s3.c1 = r->s3.c1 - s->s3.c1;
     q->s3.c2 = r->s3.c2 - s->s3.c2;
   }
#ifdef OMP
  } /* OpenMP closing brace */
#endif
}

#endif

#ifdef WITHLAPH
void diff_su3vect(su3_vector * const Q,su3_vector * const R,su3_vector * const S, const int N)
{
#ifdef OMP
#pragma omp parallel
  {
#endif

  su3_vector *q,*r,*s;

#ifdef OMP
#pragma omp for
#endif
  for (int ix = 0; ix < N; ++ix) 
  {
    q=(su3_vector *) Q + ix;
    r=(su3_vector *) R + ix;
    s=(su3_vector *) S + ix;
     
    q->c0 = r->c0 - s->c0;
    q->c1 = r->c1 - s->c1;
    q->c2 = r->c2 - s->c2;
  } 
#ifdef OMP
  } /* OpenMP closing brace */
#endif
}
#endif

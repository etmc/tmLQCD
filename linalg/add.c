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
 *   void add(spinor * const Q,spinor * const R,spinor * const S)
 *     Makes the sum (*Q) = (*R) + (*S)
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
#include "su3.h"
#include "add.h"

#if (defined BGQ && defined XLC)

void add(spinor * const Q,const spinor * const R,const spinor * const S, const int N) {
#ifdef OMP
#pragma omp parallel
  {
#endif
  vector4double x0, x1, x2, x3, x4, x5, y0, y1, y2, y3, y4, y5;
  vector4double z0, z1, z2, z3, z4, z5;
  double *q;
  double *r,*s;

  __alignx(32, s);
  __alignx(32, r);
  __alignx(32, q);
  __alignx(32, S);
  __alignx(32, R);

  __prefetch_by_load(S);
  __prefetch_by_load(R);
  __prefetch_by_stream(1, Q);

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
    z0 = vec_add(x0, y0);
    z1 = vec_add(x1, y1);
    z2 = vec_add(x2, y2);
    z3 = vec_add(x3, y3);
    z4 = vec_add(x4, y4);
    z5 = vec_add(x5, y5);
    vec_st(z0, 0, q);
    vec_st(z1, 0, q+4);
    vec_st(z2, 0, q+8);
    vec_st(z3, 0, q+12);
    vec_st(z4, 0, q+16);
    vec_st(z5, 0, q+20);
  }

#ifdef OMP
  } /*OpenMP parallel closing brace */
#endif
  return;
}

#else

/* Q output, R input, S input */
void add(spinor * const Q,const spinor * const R,const spinor * const S, const int N){
#ifdef OMP
#pragma omp parallel
  {
#endif

  int ix;
  spinor *q,*r,*s;
  
#ifdef OMP
#pragma omp for
#endif
  for (ix = 0; ix < N; ix++){
    q=(spinor *) Q + ix;
    r=(spinor *) R + ix;
    s=(spinor *) S + ix;
    
    q->s0.c0 = r->s0.c0 + s->s0.c0;
    q->s0.c1 = r->s0.c1 + s->s0.c1;
    q->s0.c2 = r->s0.c2 + s->s0.c2;
    
    q->s1.c0 = r->s1.c0 + s->s1.c0;
    q->s1.c1 = r->s1.c1 + s->s1.c1;
    q->s1.c2 = r->s1.c2 + s->s1.c2;
    
    q->s2.c0 = r->s2.c0 + s->s2.c0;
    q->s2.c1 = r->s2.c1 + s->s2.c1;
    q->s2.c2 = r->s2.c2 + s->s2.c2;
    
    q->s3.c0 = r->s3.c0 + s->s3.c0;
    q->s3.c1 = r->s3.c1 + s->s3.c1;
    q->s3.c2 = r->s3.c2 + s->s3.c2;
  }

#ifdef OMP
  } /* OpenMP closing brace */
#endif

}

#endif

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
 * File assign_add_mul.c 
 *
 *   void assign_add_mul(spinor * const P, spinor * const Q, const complex c)
 *     (*P) = (*P) + c(*Q)        c is a complex constant
 *
 *******************************************************************************/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#ifdef TM_USE_OMP
# include <omp.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "su3.h"
#include "assign_add_mul_r_32.h"

#if (defined BGQ && defined XLC)
void assign_add_mul_r_32_orphaned(spinor32 * const R, spinor32 * const S, const float c, const int N) {
  vector4double x0, x1, x2, x3, x4, x5, y0, y1, y2, y3, y4, y5;
  vector4double z0, z1, z2, z3, z4, z5, k;
  float *s, *r;
  float ALIGN32 _c;
  _c = c;
  __prefetch_by_load(S);
  __prefetch_by_load(R);

  k = vec_splats((double)_c);
  __alignx(16, s);
  __alignx(16, r);
  __alignx(16, S);
  __alignx(16, R);

#ifdef TM_USE_OMP
#pragma omp for
#else
#pragma unroll(2)
#endif
  for(int i = 0; i < N; i++) {
    s=(float*)((spinor32 *) S + i);
    r=(float*)((spinor32 *) R + i);
    __prefetch_by_load(S + i + 1);
    __prefetch_by_stream(1, R + i + 1);
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
    z0 = vec_madd(k, y0, x0);
    z1 = vec_madd(k, y1, x1);
    z2 = vec_madd(k, y2, x2);
    z3 = vec_madd(k, y3, x3);
    z4 = vec_madd(k, y4, x4);
    z5 = vec_madd(k, y5, x5);
    vec_st(z0, 0, r);
    vec_st(z1, 0, r+4);
    vec_st(z2, 0, r+8);
    vec_st(z3, 0, r+12);
    vec_st(z4, 0, r+16);
    vec_st(z5, 0, r+20);
  }
  return;
}

#else

void assign_add_mul_r_32_orphaned(spinor32 * const R, spinor32 * const S, const float c, const int N)
{
  spinor32 *r,*s;

#ifdef TM_USE_OMP
#pragma omp for
#endif
  for (int ix=0; ix<N; ix++)
  {
    r=(spinor32 *) R + ix;
    s=(spinor32 *) S + ix;

    r->s0.c0 += c * s->s0.c0;
    r->s0.c1 += c * s->s0.c1;
    r->s0.c2 += c * s->s0.c2;

    r->s1.c0 += c * s->s1.c0;
    r->s1.c1 += c * s->s1.c1;
    r->s1.c2 += c * s->s1.c2;

    r->s2.c0 += c * s->s2.c0;
    r->s2.c1 += c * s->s2.c1;
    r->s2.c2 += c * s->s2.c2;

    r->s3.c0 += c * s->s3.c0;
    r->s3.c1 += c * s->s3.c1;
    r->s3.c2 += c * s->s3.c2;
  }

}

#endif

void assign_add_mul_r_32(spinor32 * const R, spinor32 * const S, const float c, const int N)
{
#ifdef TM_USE_OMP
#pragma omp parallel
  {
#endif
  assign_add_mul_r_32_orphaned(R,S,c,N);
#ifdef TM_USE_OMP
  } /* OpenMP closing brace */
#endif
return;
}


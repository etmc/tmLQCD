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
#include <tmlqcd_config.h>
#endif
#ifdef TM_USE_OMP
#include <omp.h>
#endif
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "assign_add_mul_r_32.h"
#include "su3.h"

inline void assign_add_mul_r_32_orphaned(spinor32 *const R, spinor32 *const S, const float c,
                                         const int N) {
#ifdef TM_USE_OMP
#pragma omp parallel for
#endif
  for (int ix = 0; ix < N; ix++) {
    spinor32 *const r = (spinor32 *)R + ix;
    spinor32 *const s = (spinor32 *)S + ix;

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

void assign_add_mul_r_32(spinor32 *const R, spinor32 *const S, const float c, const int N) {
  assign_add_mul_r_32_orphaned(R, S, c, N);
  return;
}

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
#include <tmlqcd_config.h>
#endif
#ifdef TM_USE_OMP
#include <omp.h>
#endif
#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "diff.h"
#include "su3.h"

void diff(spinor *const Q, const spinor *const R, const spinor *const S, const int N) {
  /* Change due to even-odd preconditioning : VOLUME   to VOLUME/2 */

#ifdef TM_USE_OMP
#pragma omp parallel for
#endif
  for (int ix = 0; ix < N; ix++) {
    spinor *q = (spinor *)Q + ix;
    const spinor *r = (spinor *)R + ix;
    const spinor *s = (spinor *)S + ix;

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
}

void diff_ts(spinor *const Q, const spinor *const R, const spinor *const S, const int N) {
  for (int ix = 0; ix < N; ix++) {
    spinor *q = (spinor *)Q + ix;
    const spinor *r = (spinor *)R + ix;
    const spinor *s = (spinor *)S + ix;

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
}

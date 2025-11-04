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

#ifdef HAVE_CONFIG_H
#include <tmlqcd_config.h>
#endif
#include <complex.h>
#include <stdlib.h>
#ifdef TM_USE_OMP
#include <omp.h>
#endif
#include "assign_mul_add_r.h"
#include "su3.h"

/* R inoutput , c,S input*/
/*   (*R) = c*(*R) + (*S)        c is a real constant   */

void assign_mul_add_r(spinor *const R, const double c, const spinor *const S, const int N) {
  /* Change due to even-odd preconditioning : VOLUME   to VOLUME/2 */

#ifdef TM_USE_OMP
#pragma omp parallel for
#endif
  for (int ix = 0; ix < N; ++ix) {
    spinor *r = R + ix;
    const spinor *s = S + ix;

    r->s0.c0 = c * r->s0.c0 + s->s0.c0;
    r->s0.c1 = c * r->s0.c1 + s->s0.c1;
    r->s0.c2 = c * r->s0.c2 + s->s0.c2;

    r->s1.c0 = c * r->s1.c0 + s->s1.c0;
    r->s1.c1 = c * r->s1.c1 + s->s1.c1;
    r->s1.c2 = c * r->s1.c2 + s->s1.c2;

    r->s2.c0 = c * r->s2.c0 + s->s2.c0;
    r->s2.c1 = c * r->s2.c1 + s->s2.c1;
    r->s2.c2 = c * r->s2.c2 + s->s2.c2;

    r->s3.c0 = c * r->s3.c0 + s->s3.c0;
    r->s3.c1 = c * r->s3.c1 + s->s3.c1;
    r->s3.c2 = c * r->s3.c2 + s->s3.c2;
  }
}

#ifdef WITHLAPH
void assign_mul_add_r_su3vect(su3_vector *const R, const double c, su3_vector *const S,
                              const int N) {
#ifdef TM_USE_OMP
#pragma omp parallel for
#endif
  for (int ix = 0; ix < N; ++ix) {
    su3_vector *const r = R + ix;
    su3_vector *const s = S + ix;
    r->c0 = c * r->c0 + s->c0;
    r->c1 = c * r->c1 + s->c1;
    r->c2 = c * r->c2 + s->c2;
  }
}
#endif

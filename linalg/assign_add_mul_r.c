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
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#ifdef TM_USE_MPI
#include <mpi.h>
#endif
#ifdef TM_USE_OMP
#include <omp.h>
#endif
#ifdef HAVE_CONFIG_H
#include <tmlqcd_config.h>
#endif
#include "assign_add_mul_r.h"
#include "su3.h"

/*   (*P) = (*P) + c(*Q)        c is a real constant   */

void assign_add_mul_r(spinor *const P, spinor *const Q, const double c, const int N) {
#ifdef TM_USE_OMP
#pragma omp parallel for
#endif
  for (int ix = 0; ix < N; ++ix) {
    register spinor *p = P + ix;
    register spinor *q = Q + ix;
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

#ifdef WITHLAPH
void assign_add_mul_r_su3vect(su3_vector *const P, su3_vector *const Q, const double c,
                              const int N) {
#ifdef TM_USE_OMP
#pragma omp parallel for
#endif
  for (int ix = 0; ix < N; ++ix) {
    su3_vector *p = P + ix;
    su3_vector *q = Q + ix;

    p->c0 += c * q->c0;
    p->c1 += c * q->c1;
    p->c2 += c * q->c2;
  }
}
#endif

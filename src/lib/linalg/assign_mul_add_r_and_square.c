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
#ifdef TM_USE_MPI
#include <mpi.h>
#endif
#include <complex.h>
#include <stdlib.h>
#ifdef TM_USE_OMP
#include <global.h>
#include <omp.h>
#endif
#include "assign_mul_add_r_and_square.h"
#include "su3.h"

/* R inoutput , c,S input*/
/*   (*R) = c*(*R) + (*S)        c is a real constant   */

double assign_mul_add_r_and_square(spinor *const R, const double c, const spinor *const S,
                                   const int N, const int parallel) {
  double ALIGN res = 0.0;
#ifdef TM_USE_MPI
  double ALIGN mres;
#endif

#ifdef TM_USE_OMP
#pragma omp parallel
  {
    int thread_num = omp_get_thread_num();
#endif
    double ALIGN ds = 0.0;

    /* Change due to even-odd preconditioning : VOLUME   to VOLUME/2 */
#ifdef TM_USE_OMP
#pragma omp for
#endif
    for (int ix = 0; ix < N; ++ix) {
      spinor *r = R + ix;
      const spinor *s = S + ix;

      r->s0.c0 = c * r->s0.c0 + s->s0.c0;
      ds += creal(r->s0.c0) * creal(r->s0.c0) + cimag(r->s0.c0) * cimag(r->s0.c0);
      r->s0.c1 = c * r->s0.c1 + s->s0.c1;
      ds += creal(r->s0.c1) * creal(r->s0.c1) + cimag(r->s0.c1) * cimag(r->s0.c1);
      r->s0.c2 = c * r->s0.c2 + s->s0.c2;
      ds += creal(r->s0.c2) * creal(r->s0.c2) + cimag(r->s0.c2) * cimag(r->s0.c2);

      r->s1.c0 = c * r->s1.c0 + s->s1.c0;
      ds += creal(r->s1.c0) * creal(r->s1.c0) + cimag(r->s1.c0) * cimag(r->s1.c0);
      r->s1.c1 = c * r->s1.c1 + s->s1.c1;
      ds += creal(r->s1.c1) * creal(r->s1.c1) + cimag(r->s1.c1) * cimag(r->s1.c1);
      r->s1.c2 = c * r->s1.c2 + s->s1.c2;
      ds += creal(r->s1.c2) * creal(r->s1.c2) + cimag(r->s1.c2) * cimag(r->s1.c2);

      r->s2.c0 = c * r->s2.c0 + s->s2.c0;
      ds += creal(r->s2.c0) * creal(r->s2.c0) + cimag(r->s2.c0) * cimag(r->s2.c0);
      r->s2.c1 = c * r->s2.c1 + s->s2.c1;
      ds += creal(r->s2.c1) * creal(r->s2.c1) + cimag(r->s2.c1) * cimag(r->s2.c1);
      r->s2.c2 = c * r->s2.c2 + s->s2.c2;
      ds += creal(r->s2.c2) * creal(r->s2.c2) + cimag(r->s2.c2) * cimag(r->s2.c2);

      r->s3.c0 = c * r->s3.c0 + s->s3.c0;
      ds += creal(r->s3.c0) * creal(r->s3.c0) + cimag(r->s3.c0) * cimag(r->s3.c0);
      r->s3.c1 = c * r->s3.c1 + s->s3.c1;
      ds += creal(r->s3.c1) * creal(r->s3.c1) + cimag(r->s3.c1) * cimag(r->s3.c1);
      r->s3.c2 = c * r->s3.c2 + s->s3.c2;
      ds += creal(r->s3.c2) * creal(r->s3.c2) + cimag(r->s3.c2) * cimag(r->s3.c2);
    }

#ifdef TM_USE_OMP
    g_omp_acc_re[thread_num] = ds;
  } /* OpenMP closing brace */

  for (int i = 0; i < omp_num_threads; ++i) {
    res += g_omp_acc_re[i];
  }
#else
  res = ds;
#endif

#ifdef TM_USE_MPI
  if (parallel) {
    MPI_Allreduce(&res, &mres, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return (mres);
  }
#endif
  return (res);
}

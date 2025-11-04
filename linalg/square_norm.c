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
 *
 * File square_norm.c
 *
 *   double square_norm(spinor * const P )
 *     Returns the square norm of *P
 *
 *******************************************************************************/

#ifdef HAVE_CONFIG_H
#include <tmlqcd_config.h>
#endif
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#ifdef TM_USE_MPI
#include <mpi.h>
#endif
#ifdef TM_USE_OMP
#include <omp.h>
#include "global.h"
#endif
#include <complex.h>
#include "square_norm.h"
#include "su3.h"

double square_norm(const spinor *const P, const int N, const int parallel) {
  double ALIGN res = 0.0;
#ifdef TM_USE_MPI
  double ALIGN mres;
#endif

#ifdef TM_USE_OMP
#pragma omp parallel
  {
    int thread_num = omp_get_thread_num();
    g_omp_acc_re[thread_num] = 0.0;
#endif
    double ALIGN ks, kc, ds, tr, ts, tt;
    const spinor *s;

    ks = 0.0;
    kc = 0.0;

#ifdef TM_USE_OMP
#pragma omp for
#endif
    for (int ix = 0; ix < N; ix++) {
      s = P + ix;

      ds = conj(s->s0.c0) * s->s0.c0 + conj(s->s0.c1) * s->s0.c1 + conj(s->s0.c2) * s->s0.c2 +
           conj(s->s1.c0) * s->s1.c0 + conj(s->s1.c1) * s->s1.c1 + conj(s->s1.c2) * s->s1.c2 +
           conj(s->s2.c0) * s->s2.c0 + conj(s->s2.c1) * s->s2.c1 + conj(s->s2.c2) * s->s2.c2 +
           conj(s->s3.c0) * s->s3.c0 + conj(s->s3.c1) * s->s3.c1 + conj(s->s3.c2) * s->s3.c2;

      tr = ds + kc;
      ts = tr + ks;
      tt = ts - ks;
      ks = ts;
      kc = tr - tt;
    }
    kc = ks + kc;

#ifdef TM_USE_OMP
    g_omp_acc_re[thread_num] = kc;

  } /* OpenMP closing brace */

  /* having left the parallel section, we can now sum up the Kahan
     corrected sums from each thread into kc */
  for (int i = 0; i < omp_num_threads; ++i) res += g_omp_acc_re[i];
#else
  res = kc;
#endif

#ifdef TM_USE_MPI
  if (parallel) {
    MPI_Allreduce(&res, &mres, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return mres;
  }
#endif

  return res;
}

// threadsafe version

double square_norm_ts(const spinor *const P, const int N, const int parallel) {
  double ALIGN res = 0.0;
#ifdef TM_USE_MPI
  double ALIGN mres;
#endif

#ifdef TM_USE_OMP2
#pragma omp parallel reduction(+ : res)
  {
#endif
    double ALIGN ks, kc, ds, tr, ts, tt;
    const spinor *s;

    ks = 0.0;
    kc = 0.0;

#ifdef TM_USE_OMP
#pragma omp for
#endif
    for (int ix = 0; ix < N; ix++) {
      s = P + ix;

      ds = conj(s->s0.c0) * s->s0.c0 + conj(s->s0.c1) * s->s0.c1 + conj(s->s0.c2) * s->s0.c2 +
           conj(s->s1.c0) * s->s1.c0 + conj(s->s1.c1) * s->s1.c1 + conj(s->s1.c2) * s->s1.c2 +
           conj(s->s2.c0) * s->s2.c0 + conj(s->s2.c1) * s->s2.c1 + conj(s->s2.c2) * s->s2.c2 +
           conj(s->s3.c0) * s->s3.c0 + conj(s->s3.c1) * s->s3.c1 + conj(s->s3.c2) * s->s3.c2;

      tr = ds + kc;
      ts = tr + ks;
      tt = ts - ks;
      ks = ts;
      kc = tr - tt;
    }
    res = ks + kc;

#ifdef TM_USE_OMP2
  } /* OpenMP closing brace */
#endif

#ifdef TM_USE_MPI
  if (parallel) {
    MPI_Allreduce(&res, &mres, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return mres;
  }
#endif

  return res;
}

#ifdef WITHLAPH
double square_norm_su3vect(su3_vector *const P, const int N, const int parallel) {
  int ix;
  double ALIGN ks, kc, ds, tr, ts, tt;
  su3_vector *s;

  ks = 0.0;
  kc = 0.0;

  for (ix = 0; ix < N; ix++) {
    s = P + ix;

    ds = creal(s->c0) * creal(s->c0) + cimag(s->c0) * cimag(s->c0) + creal(s->c1) * creal(s->c1) +
         cimag(s->c1) * cimag(s->c1) + creal(s->c2) * creal(s->c2) + cimag(s->c2) * cimag(s->c2);

    tr = ds + kc;
    ts = tr + ks;
    tt = ts - ks;
    ks = ts;
    kc = tr - tt;
  }
  kc = ks + kc;
#ifdef TM_USE_MPI
  if (parallel) {
    MPI_Allreduce(&kc, &ks, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return ks;
  }
#endif
  return kc;
}
#endif

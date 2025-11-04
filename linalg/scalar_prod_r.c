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
 * File scalar_prod_r.c
 *
 *   double scalar_prod_r(spinor * const S,spinor * const R, const int N)
 *     Returns the real part of the scalar product (*R,*S)
 *
 *******************************************************************************/

#include <complex.h>

#ifdef HAVE_CONFIG_H
#include "tmlqcd_config.h"
#endif
#ifdef TM_USE_MPI
#include <mpi.h>
#endif

#ifdef TM_USE_OMP
#include <global.h>
#include <omp.h>
#endif
#include "scalar_prod_r.h"
#include "su3.h"

/*  R input, S input */
double scalar_prod_r(const spinor *const S, const spinor *const R, const int N,
                     const int parallel) {
  double ALIGN res = 0.0;
#ifdef TM_USE_MPI
  double ALIGN mres;
#endif

#ifdef TM_USE_OMP
#pragma omp parallel
  {
    int thread_num = omp_get_thread_num();
#endif
    double ALIGN kc, ks, ds, tr, ts, tt;

    ks = 0.0;
    kc = 0.0;

#ifdef TM_USE_OMP
#pragma omp for
#endif
    for (int ix = 0; ix < N; ++ix) {
      const spinor *const s = S + ix;
      const spinor *const r = R + ix;

      ds = creal(r->s0.c0 * conj(s->s0.c0)) + creal(r->s0.c1 * conj(s->s0.c1)) +
           creal(r->s0.c2 * conj(s->s0.c2)) + creal(r->s1.c0 * conj(s->s1.c0)) +
           creal(r->s1.c1 * conj(s->s1.c1)) + creal(r->s1.c2 * conj(s->s1.c2)) +
           creal(r->s2.c0 * conj(s->s2.c0)) + creal(r->s2.c1 * conj(s->s2.c1)) +
           creal(r->s2.c2 * conj(s->s2.c2)) + creal(r->s3.c0 * conj(s->s3.c0)) +
           creal(r->s3.c1 * conj(s->s3.c1)) + creal(r->s3.c2 * conj(s->s3.c2));

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

  for (int i = 0; i < omp_num_threads; ++i) res += g_omp_acc_re[i];
#else
  res = kc;
#endif

#if defined TM_USE_MPI
  if (parallel) {
    MPI_Allreduce(&res, &mres, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return mres;
  }
#endif
  return res;
}

#ifdef WITHLAPH
double scalar_prod_r_su3vect(su3_vector *const S, su3_vector *const R, const int N,
                             const int parallel) {
  double ALIGN ks, kc, ds, tr, ts, tt;
  ks = 0.0;
  kc = 0.0;
  for (int ix = 0; ix < N; ++ix) {
    su3_vector *const s = (su3_vector *)S + ix;
    su3_vector *const r = (su3_vector *)R + ix;

    ds = creal(r->c0) * creal(s->c0) + cimag(r->c0) * cimag(s->c0) + creal(r->c1) * creal(s->c1) +
         cimag(r->c1) * cimag(s->c1) + creal(r->c2) * creal(s->c2) + cimag(r->c2) * cimag(s->c2);

    tr = ds + kc;
    ts = tr + ks;
    tt = ts - ks;
    ks = ts;
    kc = tr - tt;
  }
  kc = ks + kc;
#if defined TM_USE_MPI
  if (parallel) {
    MPI_Allreduce(&kc, &ks, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return ks;
  }
#endif
  return kc;
}

#endif  // WITHLAPH

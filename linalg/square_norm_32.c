#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#ifdef MPI
# include <mpi.h>
#endif
#ifdef OMP
# include <omp.h>
# include "global.h"
#endif
#include <complex.h>
#include "su3.h"
#include "square_norm_32.h"


float square_norm_32(const spinor32 * const P, const int N, const int parallel)
{
  float ALIGN res = 0.0;
#ifdef MPI
  float ALIGN mres;
#endif

#ifdef OMP
#pragma omp parallel
  {
    int thread_num = omp_get_thread_num();
    g_omp_acc_re[thread_num] = 0.0;
#endif
  float ALIGN ks,kc,ds,tr,ts,tt;
  const spinor32 *s;
  
  ks = 0.0;
  kc = 0.0;
  
#ifdef OMP
#pragma omp for
#endif    
  for (int ix  =  0; ix < N; ix++) {
    s = P + ix;
    
    ds = conj(s->s0.c0) * s->s0.c0 +
         conj(s->s0.c1) * s->s0.c1 +
         conj(s->s0.c2) * s->s0.c2 +
         conj(s->s1.c0) * s->s1.c0 +
         conj(s->s1.c1) * s->s1.c1 +
         conj(s->s1.c2) * s->s1.c2 +
         conj(s->s2.c0) * s->s2.c0 +
         conj(s->s2.c1) * s->s2.c1 +
         conj(s->s2.c2) * s->s2.c2 +
         conj(s->s3.c0) * s->s3.c0 +
         conj(s->s3.c1) * s->s3.c1 +
         conj(s->s3.c2) * s->s3.c2;

    tr = ds + kc;
    ts = tr + ks;
    tt = ts-ks;
    ks = ts;
    kc = tr-tt;
  }
  kc=ks+kc;

#ifdef OMP
  g_omp_acc_re[thread_num] = kc;

  } /* OpenMP closing brace */

  /* having left the parallel section, we can now sum up the Kahan
     corrected sums from each thread into kc */
  for(int i = 0; i < omp_num_threads; ++i)
    res += g_omp_acc_re[i];
#else
  res = kc;
#endif

#  ifdef MPI
  if(parallel) {
    MPI_Allreduce(&res, &mres, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    return mres;
  }
#endif

  return res;
}


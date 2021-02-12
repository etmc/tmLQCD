#ifdef HAVE_CONFIG_H
#include "tmlqcd_config.h"
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#ifdef TM_USE_MPI
#include <mpi.h>
#endif
#ifdef TM_USE_OMP
#include <omp.h>
#include "global.h"
#endif
#include <complex.h>
#include "su3.h"
#include "square_norm_32.h"

#if (defined BGQ && defined XLC)

float square_norm_32(spinor32 * const P, const int N, const int parallel) {
  float ALIGN32 res = 0.0;
#ifdef TM_USE_MPI
  float ALIGN32 mres;
#endif

#ifdef TM_USE_OMP
#pragma omp parallel
  {
    int thread_num = omp_get_thread_num();
#endif
  vector4double x0, x1, x2, x3, x4, x5, y0, y1, y2, y3, y4, y5;
  vector4double ds,tt,tr,ts,kc,ks,buffer;
  float *s ALIGN32;

  ks = vec_splats(0.);
  kc = vec_splats(0.);

#ifndef TM_USE_OMP
#pragma unroll(4)
#else
#pragma omp for
#endif
  for(int i = 0; i < N; i++) {
    s = (float*)((spinor32*) P+i);
    __prefetch_by_load(P+i+1);
    x0 = vec_ld(0, s);
    x1 = vec_ld(0, s+4);
    x2 = vec_ld(0, s+8);
    x3 = vec_ld(0, s+12);
    x4 = vec_ld(0, s+16);
    x5 = vec_ld(0, s+20);
    y0 = vec_mul(x0, x0);
    y1 = vec_mul(x1, x1);
    y2 = vec_mul(x2, x2);
    y3 = vec_mul(x3, x3);
    y4 = vec_mul(x4, x4);
    y5 = vec_mul(x5, x5);

    x0 = vec_add(y0, y1);
    x1 = vec_add(y2, y3);
    x2 = vec_add(y4, y5);
    x3 = vec_add(x0, x1);
    ds = vec_add(x2, x3);

    tr = vec_add(ds, kc);
    ts = vec_add(tr, ks);
    tt = vec_sub(ts, ks);
    ks = ts;
    kc = vec_sub(tr, tt);
  }
  buffer = vec_add(kc,ks);

#ifdef TM_USE_OMP
  g_omp_acc_re[thread_num] = buffer[0] + buffer[1] + buffer[2] + buffer[3];
  } /* OpenMP closing brace */

  for(int i = 0; i < omp_num_threads; ++i)
    res += g_omp_acc_re[i];
#else
  res = buffer[0] + buffer[1] + buffer[2] + buffer[3];
#endif

#  ifdef TM_USE_MPI
  if(parallel) {
    MPI_Allreduce(&res, &mres, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    return mres;
  }
#  endif

  return res;
}


#else 
float square_norm_32(const spinor32 * const P, const int N, const int parallel)
{
  float ALIGN32 res = 0.0;
#ifdef TM_USE_MPI
  float ALIGN32 mres;
#endif

#ifdef TM_USE_OMP
#pragma omp parallel
  {
    int thread_num = omp_get_thread_num();
    g_omp_acc_re[thread_num] = 0.0;
#endif
  float ALIGN32 ks,kc,ds,tr,ts,tt;
  const spinor32 *s;
  
  ks = 0.0;
  kc = 0.0;
  
#ifdef TM_USE_OMP
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

#ifdef TM_USE_OMP
  g_omp_acc_re[thread_num] = kc;

  } /* OpenMP closing brace */

  /* having left the parallel section, we can now sum up the Kahan
     corrected sums from each thread into kc */
  for(int i = 0; i < omp_num_threads; ++i)
    res += g_omp_acc_re[i];
#else
  res = kc;
#endif

#  ifdef TM_USE_MPI
  if(parallel) {
    MPI_Allreduce(&res, &mres, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    return mres;
  }
#endif

  return res;
}

#endif

// threadsafe version

float square_norm_ts_32(const spinor32 * const P, const int N, const int parallel)
{
  float ALIGN32 res = 0.0;
#ifdef TM_USE_MPI
  float ALIGN32 mres;
#endif

#ifdef TM_USE_OMP2
#pragma omp parallel reduction(+:res)
  {
#endif
  float ALIGN32 ks,kc,ds,tr,ts,tt;
  const spinor32 *s;
  
  ks = 0.0;
  kc = 0.0;
  
#ifdef TM_USE_OMP2
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
  res=ks+kc;
#ifdef TM_USE_OMP2
  } /* OpenMP closing brace */
#endif

#  ifdef TM_USE_MPI
  if(parallel) {
    MPI_Allreduce(&res, &mres, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    return mres;
  }
#endif

  return res;
}

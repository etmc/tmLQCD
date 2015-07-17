#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#ifdef _USE_MPI
# include <mpi.h>
#endif
#ifdef OMP
# include <omp.h>
# include <global.h>
#endif
#include "su3.h"
#include "scalar_prod_r_32.h"

/*  R input, S input */

#include <complex.h>
#if (defined BGQ && defined XLC)

float scalar_prod_r_32(const spinor32 * const S, const spinor32 * const R, const int N, const int parallel) {
  float ALIGN32 res = 0.0;
#ifdef _USE_MPI
  float ALIGN32 mres;
#endif

#ifdef OMP
#pragma omp parallel
  {
  int thread_num = omp_get_thread_num();
#endif
  vector4double ks, kc, ds, tr, ts, tt;
  vector4double x0, x1, x2, x3, x4, x5, y0, y1, y2, y3, y4, y5;
  vector4double z0, z1, z2, z3, z4, z5;
  float *s, *r;
  vector4double buffer;
  __alignx(16, s);
  __alignx(16, r);
  __alignx(16, S);
  __alignx(16, R);

  __prefetch_by_load(S);
  __prefetch_by_load(R);

  ks = vec_splats(0.0);
  kc = vec_splats(0.0);

#ifndef OMP
#pragma unroll(2)
#else
#pragma omp for
#endif
  for (int ix = 0; ix < N; ++ix) {
    s=(float*)((spinor32 *) S + ix);
    r=(float*)((spinor32 *) R + ix);
    __prefetch_by_load(S + ix + 1);
    __prefetch_by_load(R + ix + 1);
    x0 = vec_ld(0, s);
    x1 = vec_ld(0, s+4);
    x2 = vec_ld(0, s+8);
    x3 = vec_ld(0, s+12);
    x4 = vec_ld(0, s+16);
    x5 = vec_ld(0, s+20);
    y0 = vec_ld(0, r);
    y1 = vec_ld(0, r+4);
    y2 = vec_ld(0, r+8);
    y3 = vec_ld(0, r+12);
    y4 = vec_ld(0, r+16);
    y5 = vec_ld(0, r+20);
    z0 = vec_mul(x0, y0);
    z1 = vec_mul(x1, y1);
    z2 = vec_mul(x2, y2);
    z3 = vec_mul(x3, y3);
    z4 = vec_mul(x4, y4);
    z5 = vec_mul(x5, y5);
    x0 = vec_add(z0, z1);
    x1 = vec_add(z2, z3);
    x2 = vec_add(z4, z5);
    x3 = vec_add(x0, x1);
    ds = vec_add(x2, x3);

    tr = vec_add(ds, kc);
    ts = vec_add(tr, ks);
    tt = vec_sub(ts, ks);
    ks = ts;
    kc = vec_sub(tr, tt);
  }
  buffer = vec_add(kc, ks);

#ifdef OMP
  g_omp_acc_re[thread_num] = buffer[0] + buffer[1] + buffer[2] + buffer[3];
  } /* OpenMP parallel closing brace */
  for( int i = 0; i < omp_num_threads; ++i)
    res += g_omp_acc_re[i];
#else
  res = buffer[0] + buffer[1] + buffer[2] + buffer[3]; 
#endif

#if defined _USE_MPI
  if(parallel) {
    MPI_Allreduce(&res, &mres, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    return(mres);
  }
#endif

  return (res);
}

#else

float scalar_prod_r_32(const spinor32 * const S, const spinor32 * const R, const int N, const int parallel)
{
  float ALIGN32 res = 0.0;
#ifdef _USE_MPI
  float ALIGN32 mres;
#endif

#ifdef OMP
#pragma omp parallel
  {
  int thread_num = omp_get_thread_num();
#endif
  float ALIGN32 kc,ks,ds,tr,ts,tt;
  const spinor32 *s,*r;

  ks = 0.0;
  kc = 0.0;

#ifdef OMP
#pragma omp for
#endif
  for (int ix = 0; ix < N; ++ix) {
    s = S + ix;
    r = R + ix;
    
    ds = creal(r->s0.c0 * conj(s->s0.c0)) + creal(r->s0.c1 * conj(s->s0.c1)) + creal(r->s0.c2 * conj(s->s0.c2)) +
      creal(r->s1.c0 * conj(s->s1.c0)) + creal(r->s1.c1 * conj(s->s1.c1)) + creal(r->s1.c2 * conj(s->s1.c2)) +
      creal(r->s2.c0 * conj(s->s2.c0)) + creal(r->s2.c1 * conj(s->s2.c1)) + creal(r->s2.c2 * conj(s->s2.c2)) +
      creal(r->s3.c0 * conj(s->s3.c0)) + creal(r->s3.c1 * conj(s->s3.c1)) + creal(r->s3.c2 * conj(s->s3.c2));    
    
    tr=ds+kc;
    ts=tr+ks;
    tt=ts-ks;
    ks=ts;
    kc=tr-tt;
  }
  kc=ks+kc;

#ifdef OMP
  g_omp_acc_re[thread_num] = kc;

  } /* OpenMP closing brace */

  for(int i = 0; i < omp_num_threads; ++i)
    res += g_omp_acc_re[i];
#else
  res = kc;
#endif

#if defined _USE_MPI
  if(parallel)
  {
    MPI_Allreduce(&res, &mres, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    return mres;
  }
#endif
  return res;
}

#endif

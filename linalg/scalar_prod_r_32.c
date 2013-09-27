#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#ifdef MPI
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


float scalar_prod_r_32(const spinor32 * const S, const spinor32 * const R, const int N, const int parallel)
{
  float ALIGN res = 0.0;
#ifdef MPI
  float ALIGN mres;
#endif

#ifdef OMP
#pragma omp parallel
  {
  int thread_num = omp_get_thread_num();
#endif
  float ALIGN kc,ks,ds,tr,ts,tt;
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

#if defined MPI
  if(parallel)
  {
    MPI_Allreduce(&res, &mres, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    return mres;
  }
#endif
  return res;
}



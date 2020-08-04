#ifdef HAVE_CONFIG_H
#include "tmlqcd_config.h"
#endif
#include <stdlib.h>
#include <complex.h>
#ifdef TM_USE_OMP
#include <omp.h>
#endif
#include "su3.h"
#include "assign_mul_add_r_32.h"


/* R inoutput , c,S input*/
/*   (*R) = c*(*R) + (*S)        c is a real constant   */

#if (defined BGQ && defined XLC)

void assign_mul_add_r_32(spinor32 * const R, const float c, const spinor32 * const S, const int N) {
#ifdef TM_USE_OMP
#pragma omp parallel
  {
#endif

  vector4double x0, x1, x2, x3, x4, x5, y0, y1, y2, y3, y4, y5;
  vector4double z0, z1, z2, z3, z4, z5, k;
  float *s, *r;
  float ALIGN32 _c;
  _c = c;
  __prefetch_by_load(S);
  __prefetch_by_load(R);

  k = vec_splats((double)_c);
  __alignx(16, s);
  __alignx(16, r);
  __alignx(16, S);
  __alignx(16, R);

#ifdef TM_USE_OMP
#pragma omp for
#else
#pragma unroll(4)
#endif
  for(int i = 0; i < N; i++) {
    s=(float*)((spinor32 *) S + i);
    r=(float*)((spinor32 *) R + i);
    __prefetch_by_load(S + i + 1);
    __prefetch_by_stream(1, R + i + 1);
    x0 = vec_ld(0, r);
    x1 = vec_ld(0, r+4);
    x2 = vec_ld(0, r+8);
    x3 = vec_ld(0, r+12);
    x4 = vec_ld(0, r+16);
    x5 = vec_ld(0, r+20);
    y0 = vec_ld(0, s);
    y1 = vec_ld(0, s+4);
    y2 = vec_ld(0, s+8);
    y3 = vec_ld(0, s+12);
    y4 = vec_ld(0, s+16);
    y5 = vec_ld(0, s+20);
    z0 = vec_madd(k, x0, y0);
    z1 = vec_madd(k, x1, y1);
    z2 = vec_madd(k, x2, y2);
    z3 = vec_madd(k, x3, y3);
    z4 = vec_madd(k, x4, y4);
    z5 = vec_madd(k, x5, y5);
    vec_st(z0, 0, r);
    vec_st(z1, 0, r+4);
    vec_st(z2, 0, r+8);
    vec_st(z3, 0, r+12);
    vec_st(z4, 0, r+16);
    vec_st(z5, 0, r+20);
  }
#ifdef TM_USE_OMP
  } /* OpenMP closing brace */
#endif  
  return;
}

#else

void assign_mul_add_r_32(spinor32 * const R, const float c, const spinor32 * const S, const int N)
{
#ifdef TM_USE_OMP
#pragma omp parallel
  {
#endif
  spinor32 *r;
  const spinor32 *s;
  
  /* Change due to even-odd preconditioning : VOLUME   to VOLUME/2 */   
#ifdef TM_USE_OMP
#pragma omp for
#endif
  for (int ix = 0; ix < N; ++ix)
  {
    r = R + ix;
    s = S + ix;
    
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
#ifdef TM_USE_OMP
  } /* OpenMP closing brace */
#endif
}


#endif

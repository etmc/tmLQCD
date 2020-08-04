#ifdef HAVE_CONFIG_H
#include "tmlqcd_config.h"
#endif
#ifdef TM_USE_OMP
#include <omp.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "su3.h"
#include "addto_32.h"



/* Q output, R input, S input */
void addto_32(spinor * const Q, const spinor32 * const R, const int N)
{
#ifdef TM_USE_OMP
#pragma omp parallel
  {
#endif

  int ix;
  spinor *q;
  spinor32 * r;
#ifdef TM_USE_OMP
#pragma omp for
#endif
  for (ix = 0; ix < N; ix++){
    q=(spinor *) Q + ix;
    r=(spinor32 *) R + ix;

    
    q->s0.c0 += r->s0.c0;
    q->s0.c1 += r->s0.c1;
    q->s0.c2 += r->s0.c2;
    
    q->s1.c0 += r->s1.c0;
    q->s1.c1 += r->s1.c1;
    q->s1.c2 += r->s1.c2;
    
    q->s2.c0 += r->s2.c0;
    q->s2.c1 += r->s2.c1;
    q->s2.c2 += r->s2.c2;
    
    q->s3.c0 += r->s3.c0;
    q->s3.c1 += r->s3.c1;
    q->s3.c2 += r->s3.c2;
  }

#ifdef TM_USE_OMP
  } /* OpenMP closing brace */
#endif

}

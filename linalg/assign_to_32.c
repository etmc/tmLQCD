#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#ifdef OMP
# include <omp.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <math.h>
#include <string.h>
#include "su3.h"
#include "assign_to_32.h"


/* S input, R output        */
/* S and R must not overlap */
void assign_to_32(spinor32 * const R, spinor * const S, const int N)
{
#ifdef OMP
#pragma omp parallel
  {
#endif
  spinor32 *r;
  spinor *s;

#ifdef OMP
#pragma omp for
#endif
  for (int ix=0; ix<N; ix++)
  {
    r=(spinor32 *) R + ix;
    s=(spinor*) S + ix;

    r->s0.c0 = s->s0.c0;
    r->s0.c1 = s->s0.c1;
    r->s0.c2 = s->s0.c2;

    r->s1.c0 = s->s1.c0;
    r->s1.c1 = s->s1.c1;
    r->s1.c2 = s->s1.c2;

    r->s2.c0 = s->s2.c0;
    r->s2.c1 = s->s2.c1;
    r->s2.c2 = s->s2.c2;

    r->s3.c0 = s->s3.c0;
    r->s3.c1 = s->s3.c1;
    r->s3.c2 = s->s3.c2;
  }

#ifdef OMP
  } /* OpenMP closing brace */
#endif
 
  return;
}





/* S input, R output        */
/* S and R must not overlap */
void assign_to_64(spinor * const R, spinor32 * const S, const int N)
{
#ifdef OMP
#pragma omp parallel
  {
#endif
  spinor *r;
  spinor32 *s;

#ifdef OMP
#pragma omp for
#endif
  for (int ix=0; ix<N; ix++)
  {
    r=(spinor *) R + ix;
    s=(spinor32 *) S + ix;

    r->s0.c0 = s->s0.c0;
    r->s0.c1 = s->s0.c1;
    r->s0.c2 = s->s0.c2;

    r->s1.c0 = s->s1.c0;
    r->s1.c1 = s->s1.c1;
    r->s1.c2 = s->s1.c2;

    r->s2.c0 = s->s2.c0;
    r->s2.c1 = s->s2.c1;
    r->s2.c2 = s->s2.c2;

    r->s3.c0 = s->s3.c0;
    r->s3.c1 = s->s3.c1;
    r->s3.c2 = s->s3.c2;
  }

#ifdef OMP
  } /* OpenMP closing brace */
#endif
 
  return;
}










#include "su3.h"

/* Q output, R input, S input */
void ratio(spinor * const Q,const spinor * const R,const spinor * const S, const int N){
#ifdef TM_USE_OMP
#pragma omp parallel
  {
#endif

  int ix;
  spinor *q,*r,*s;
  
#ifdef TM_USE_OMP
#pragma omp for
#endif
  for (ix = 0; ix < N; ix++){
    q=(spinor *) Q + ix;
    r=(spinor *) R + ix;
    s=(spinor *) S + ix;
    
    q->s0.c0 = creal(r->s0.c0) / creal(s->s0.c0) + I*cimag(r->s0.c0) / cimag(s->s0.c0);
    q->s0.c1 = creal(r->s0.c1) / creal(s->s0.c1) + I*cimag(r->s0.c1) / cimag(s->s0.c1);
    q->s0.c2 = creal(r->s0.c2) / creal(s->s0.c2) + I*cimag(r->s0.c2) / cimag(s->s0.c2);
    
    q->s1.c0 = creal(r->s1.c0) / creal(s->s1.c0) + I*cimag(r->s1.c0) / cimag(s->s1.c0);
    q->s1.c1 = creal(r->s1.c1) / creal(s->s1.c1) + I*cimag(r->s1.c1) / cimag(s->s1.c1);
    q->s1.c2 = creal(r->s1.c2) / creal(s->s1.c2) + I*cimag(r->s1.c2) / cimag(s->s1.c2);
    
    q->s2.c0 = creal(r->s2.c0) / creal(s->s2.c0) + I*cimag(r->s2.c0) / cimag(s->s2.c0);
    q->s2.c1 = creal(r->s2.c1) / creal(s->s2.c1) + I*cimag(r->s2.c1) / cimag(s->s2.c1);
    q->s2.c2 = creal(r->s2.c2) / creal(s->s2.c2) + I*cimag(r->s2.c2) / cimag(s->s2.c2);
    
    q->s3.c0 = creal(r->s3.c0) / creal(s->s3.c0) + I*cimag(r->s3.c0) / cimag(s->s3.c0);
    q->s3.c1 = creal(r->s3.c1) / creal(s->s3.c1) + I*cimag(r->s3.c1) / cimag(s->s3.c1);
    q->s3.c2 = creal(r->s3.c2) / creal(s->s3.c2) + I*cimag(r->s3.c2) / cimag(s->s3.c2);
  }

#ifdef TM_USE_OMP
  } /* OpenMP closing brace */
#endif

}


#include "su3.h"
#include <stdio.h>

/* Q input */
void print_spinor(spinor const * const Q, const int N){
  int ix;
  spinor *q;
  for (ix = 0; ix < N; ix++){
    q=(spinor *) Q + ix;
    
    printf("ix: %d s0: (%f,%f)  (%f,%f) (%f,%f)\n",   ix, creal(q->s0.c0), cimag(q->s0.c0), 
        creal(q->s0.c1), cimag(q->s0.c1), creal(q->s0.c2), cimag(q->s0.c2));
    printf("ix: %d s1: (%f,%f)  (%f,%f) (%f,%f)\n",   ix, creal(q->s1.c0), cimag(q->s1.c0), 
        creal(q->s1.c1), cimag(q->s1.c1), creal(q->s1.c2), cimag(q->s1.c2));
    printf("ix: %d s2: (%f,%f)  (%f,%f) (%f,%f)\n",   ix, creal(q->s2.c0), cimag(q->s2.c0), 
        creal(q->s2.c1), cimag(q->s2.c1), creal(q->s2.c2), cimag(q->s2.c2));
    printf("ix: %d s3: (%f,%f)  (%f,%f) (%f,%f)\n\n", ix, creal(q->s3.c0), cimag(q->s3.c0), 
        creal(q->s3.c1), cimag(q->s3.c1), creal(q->s3.c2), cimag(q->s3.c2));
  }
}


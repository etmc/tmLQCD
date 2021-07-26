#include "su3.h"
#include <stdio.h>

/* Q input */
void print_spinor_similar_components(spinor const * const Q, spinor const * const P, const int N, const double thresh){
  spinor *q, *p;
  for (int ix = 0; ix < N; ix++){
    for (int iy = 0; iy < N; iy++){
      q=(spinor *) Q + ix;
      p=(spinor *) P + iy;
      
      if( fabs( (creal(q->s0.c0) - creal(p->s0.c0)) / creal(q->s0.c0) )  < thresh  &&
          fabs( (cimag(q->s0.c0) - cimag(p->s0.c0)) / cimag(q->s0.c0) )  < thresh  &&
          fabs( (creal(q->s0.c1) - creal(p->s0.c1)) / creal(q->s0.c1) )  < thresh  &&
          fabs( (cimag(q->s0.c1) - cimag(p->s0.c1)) / cimag(q->s0.c1) )  < thresh  &&
          fabs( (creal(q->s0.c2) - creal(p->s0.c2)) / creal(q->s0.c2) )  < thresh  &&
          fabs( (cimag(q->s0.c2) - cimag(p->s0.c2)) / cimag(q->s0.c2) )  < thresh  &&
          fabs( (creal(q->s1.c0) - creal(p->s1.c0)) / creal(q->s1.c0) )  < thresh  &&
          fabs( (cimag(q->s1.c0) - cimag(p->s1.c0)) / cimag(q->s1.c0) )  < thresh  &&
          fabs( (creal(q->s1.c1) - creal(p->s1.c1)) / creal(q->s1.c1) )  < thresh  &&
          fabs( (cimag(q->s1.c1) - cimag(p->s1.c1)) / cimag(q->s1.c1) )  < thresh  &&
          fabs( (creal(q->s1.c2) - creal(p->s1.c2)) / creal(q->s1.c2) )  < thresh  && 
          fabs( (cimag(q->s1.c2) - cimag(p->s1.c2)) / cimag(q->s1.c2) )  < thresh  &&
          fabs( (creal(q->s2.c0) - creal(p->s2.c0)) / creal(q->s2.c0) )  < thresh  &&
          fabs( (cimag(q->s2.c0) - cimag(p->s2.c0)) / cimag(q->s2.c0) )  < thresh  &&
          fabs( (creal(q->s2.c1) - creal(p->s2.c1)) / creal(q->s2.c1) )  < thresh  &&
          fabs( (cimag(q->s2.c1) - cimag(p->s2.c1)) / cimag(q->s2.c1) )  < thresh  &&
          fabs( (creal(q->s2.c2) - creal(p->s2.c2)) / creal(q->s2.c2) )  < thresh  &&
          fabs( (cimag(q->s2.c2) - cimag(p->s2.c2)) / cimag(q->s2.c2) )  < thresh  &&
          fabs( (creal(q->s3.c0) - creal(p->s3.c0)) / creal(q->s3.c0) )  < thresh  &&
          fabs( (cimag(q->s3.c0) - cimag(p->s3.c0)) / cimag(q->s3.c0) )  < thresh  &&
          fabs( (creal(q->s3.c1) - creal(p->s3.c1)) / creal(q->s3.c1) )  < thresh  &&
          fabs( (cimag(q->s3.c1) - cimag(p->s3.c1)) / cimag(q->s3.c1) )  < thresh  &&
          fabs( (creal(q->s3.c2) - creal(p->s3.c2)) / creal(q->s3.c2) )  < thresh  &&
          fabs( (cimag(q->s3.c2) - cimag(p->s3.c2)) / cimag(q->s3.c2) )  < thresh  ){

        printf("ix: %d iy: %d s0.c0: [Q,P]:(%f,%f)(%f,%f)\n",ix,iy,creal(q->s0.c0),cimag(q->s0.c0),creal(p->s0.c0),cimag(p->s0.c0)); 
        printf("ix: %d iy: %d s0.c1: [Q,P]:(%f,%f)(%f,%f)\n",ix,iy,creal(q->s0.c0),cimag(q->s0.c0),creal(p->s0.c0),cimag(p->s0.c0)); 
        printf("ix: %d iy: %d s0.c2: [Q,P]:(%f,%f)(%f,%f)\n",ix,iy,creal(q->s0.c0),cimag(q->s0.c0),creal(p->s0.c0),cimag(p->s0.c0)); 

        printf("ix: %d iy: %d s1.c0: [Q,P]:(%f,%f)(%f,%f)\n",ix,iy,creal(q->s1.c0),cimag(q->s1.c0),creal(p->s1.c0),cimag(p->s1.c0)); 
        printf("ix: %d iy: %d s1.c1: [Q,P]:(%f,%f)(%f,%f)\n",ix,iy,creal(q->s1.c1),cimag(q->s1.c1),creal(p->s1.c1),cimag(p->s1.c1)); 
        printf("ix: %d iy: %d s1.c2: [Q,P]:(%f,%f)(%f,%f)\n",ix,iy,creal(q->s1.c2),cimag(q->s1.c2),creal(p->s1.c2),cimag(p->s1.c2)); 

        printf("ix: %d iy: %d s2.c0: [Q,P]:(%f,%f)(%f,%f)\n",ix,iy,creal(q->s2.c0),cimag(q->s2.c0),creal(p->s2.c0),cimag(p->s2.c0)); 
        printf("ix: %d iy: %d s2.c1: [Q,P]:(%f,%f)(%f,%f)\n",ix,iy,creal(q->s2.c1),cimag(q->s2.c1),creal(p->s2.c1),cimag(p->s2.c1)); 
        printf("ix: %d iy: %d s2.c2: [Q,P]:(%f,%f)(%f,%f)\n",ix,iy,creal(q->s2.c2),cimag(q->s2.c2),creal(p->s2.c2),cimag(p->s2.c2)); 

        printf("ix: %d iy: %d s3.c0: [Q,P]:(%f,%f)(%f,%f)\n",ix,iy,creal(q->s3.c0),cimag(q->s3.c0),creal(p->s3.c0),cimag(p->s3.c0)); 
        printf("ix: %d iy: %d s3.c1: [Q,P]:(%f,%f)(%f,%f)\n",ix,iy,creal(q->s3.c1),cimag(q->s3.c1),creal(p->s3.c1),cimag(p->s3.c1)); 
        printf("ix: %d iy: %d s3.c2: [Q,P]:(%f,%f)(%f,%f)\n",ix,iy,creal(q->s3.c2),cimag(q->s3.c2),creal(p->s3.c2),cimag(p->s3.c2)); 
      }
    }
  }
}


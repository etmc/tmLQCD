#include "su3.h"
#include <stdio.h>
#include <float.h>

int static inline match(const complex double a, const complex double b, const double thresh){
  // either a and b are zero or very close to zero or they satisfy our threshold
  // note that in the second expression we may have divide by zero
  // so we evaluate it in this order
  return ( ( fabs(creal(a) - creal(b)) < 4*DBL_EPSILON &&
             fabs(cimag(a) - cimag(b)) < 4*DBL_EPSILON ) ||
           ( fabs( (creal(a) - creal(b)) / creal(b) ) < thresh &&
             fabs( (cimag(a) - cimag(b)) / cimag(b) ) < thresh ) );
}

            

/* Q input */
void print_spinor_similar_components(spinor const * const Q, spinor const * const P, const int N, const double thresh){
  spinor *q, *p;
  int n_matches = 0;
  for (int ix = 0; ix < N; ix++){
    for (int iy = 0; iy < N; iy++){
      q=(spinor *) Q + ix;
      p=(spinor *) P + iy;
      
      if( match(q->s0.c0, p->s0.c0, thresh) &&
          match(q->s0.c1, p->s0.c1, thresh) &&
          match(q->s0.c2, p->s0.c2, thresh) &&
          match(q->s1.c0, p->s1.c0, thresh) &&
          match(q->s1.c1, p->s1.c1, thresh) &&
          match(q->s1.c2, p->s1.c2, thresh) && 
          match(q->s2.c0, p->s2.c0, thresh) &&
          match(q->s2.c1, p->s2.c1, thresh) &&
          match(q->s2.c2, p->s2.c2, thresh) &&
          match(q->s3.c0, p->s3.c0, thresh) &&
          match(q->s3.c1, p->s3.c1, thresh) &&
          match(q->s3.c2, p->s3.c2, thresh) ){

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
        
        n_matches++;
        break; 
      }
    }
  }
  printf("%d out of %d match\n\n", n_matches, N);
}


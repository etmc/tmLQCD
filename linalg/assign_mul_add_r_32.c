#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <complex.h>
#ifdef OMP
# include <omp.h>
#endif
#include "su3.h"
#include "assign_mul_add_r_32.h"


/* R inoutput , c,S input*/
/*   (*R) = c*(*R) + (*S)        c is a real constant   */

void assign_mul_add_r_32(spinor32 * const R, const float c, const spinor32 * const S, const int N)
{
#ifdef OMP
#pragma omp parallel
  {
#endif
  spinor32 *r;
  const spinor32 *s;
  
  /* Change due to even-odd preconditioning : VOLUME   to VOLUME/2 */   
#ifdef OMP
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
#ifdef OMP
  } /* OpenMP closing brace */
#endif
}

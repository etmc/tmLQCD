/*******************************************************************************
 * $Id$
 *
 * File assign_mul_bra_add_mul_ket_add.c
 *
 *   void assign_mul_bra_add_mul_ket_add
 *   (spinor * const R,spinor * const S,spinor * const U,const double c1,const double c2)
 *     (*R) = c2*(*R + c1*(*S)) + (*U)  with c1 and c2 complex variables
 *
 *******************************************************************************/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "su3.h" 
#include "assign_mul_bra_add_mul_ket_add_r.h"

void assign_mul_bra_add_mul_ket_add_r(spinor * const R, spinor * const S, spinor * const U, 
				      const double c1, const double c2, const int N) {
   int ix;
   spinor *r,*s,*t;

   for (ix = 0; ix < N; ix++) {
     r=R+ix;
     s=S+ix;
     t=U+ix;
     (*r).s0.c0.re=c2*( (*r).s0.c0.re+c1*(*s).s0.c0.re )+(*t).s0.c0.re;
     (*r).s0.c0.im=c2*( (*r).s0.c0.im+c1*(*s).s0.c0.im )+(*t).s0.c0.im;
     (*r).s0.c1.re=c2*( (*r).s0.c1.re+c1*(*s).s0.c1.re )+(*t).s0.c1.re;
     (*r).s0.c1.im=c2*( (*r).s0.c1.im+c1*(*s).s0.c1.im )+(*t).s0.c1.im;
     (*r).s0.c2.re=c2*( (*r).s0.c2.re+c1*(*s).s0.c2.re )+(*t).s0.c2.re;
     (*r).s0.c2.im=c2*( (*r).s0.c2.im+c1*(*s).s0.c2.im )+(*t).s0.c2.im;
     
     (*r).s1.c0.re=c2*( (*r).s1.c0.re+c1*(*s).s1.c0.re )+(*t).s1.c0.re;
     (*r).s1.c0.im=c2*( (*r).s1.c0.im+c1*(*s).s1.c0.im )+(*t).s1.c0.im;
     (*r).s1.c1.re=c2*( (*r).s1.c1.re+c1*(*s).s1.c1.re )+(*t).s1.c1.re;
     (*r).s1.c1.im=c2*( (*r).s1.c1.im+c1*(*s).s1.c1.im )+(*t).s1.c1.im;
     (*r).s1.c2.re=c2*( (*r).s1.c2.re+c1*(*s).s1.c2.re )+(*t).s1.c2.re;
     (*r).s1.c2.im=c2*( (*r).s1.c2.im+c1*(*s).s1.c2.im )+(*t).s1.c2.im;
     
     (*r).s2.c0.re=c2*( (*r).s2.c0.re+c1*(*s).s2.c0.re )+(*t).s2.c0.re;
     (*r).s2.c0.im=c2*( (*r).s2.c0.im+c1*(*s).s2.c0.im )+(*t).s2.c0.im;
     (*r).s2.c1.re=c2*( (*r).s2.c1.re+c1*(*s).s2.c1.re )+(*t).s2.c1.re;
     (*r).s2.c1.im=c2*( (*r).s2.c1.im+c1*(*s).s2.c1.im )+(*t).s2.c1.im;
     (*r).s2.c2.re=c2*( (*r).s2.c2.re+c1*(*s).s2.c2.re )+(*t).s2.c2.re;
     (*r).s2.c2.im=c2*( (*r).s2.c2.im+c1*(*s).s2.c2.im )+(*t).s2.c2.im;
     
     (*r).s3.c0.re=c2*( (*r).s3.c0.re+c1*(*s).s3.c0.re )+(*t).s3.c0.re;
     (*r).s3.c0.im=c2*( (*r).s3.c0.im+c1*(*s).s3.c0.im )+(*t).s3.c0.im;
     (*r).s3.c1.re=c2*( (*r).s3.c1.re+c1*(*s).s3.c1.re )+(*t).s3.c1.re;
     (*r).s3.c1.im=c2*( (*r).s3.c1.im+c1*(*s).s3.c1.im )+(*t).s3.c1.im;
     (*r).s3.c2.re=c2*( (*r).s3.c2.re+c1*(*s).s3.c2.re )+(*t).s3.c2.re;
     (*r).s3.c2.im=c2*( (*r).s3.c2.im+c1*(*s).s3.c2.im )+(*t).s3.c2.im;
   }
}

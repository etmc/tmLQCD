/*******************************************************************************
 * $Id$
 *
 * File assign_mul_bra_add_mul_r.c
 *
 *   void assign_mul_bra_add_mul_r(spinor * const R,const double c0, const double c,spinor * const S)
 *     (*R) = c0*(*R + c*(*S))
 *
 *******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "su3.h"
#include "assign_mul_bra_add_mul_r.h"

/*  R output, S input, c0 input, c input */
void assign_mul_bra_add_mul_r(spinor * const R,const double c0, const double c,spinor * const S, const int N){
  int ix;
  static double fact0,fact;
  spinor *r,*s;
  
  fact=c;
  fact0=c0;

  for (ix = 0;ix < N; ix++){
    r=(spinor *) R + ix;
    s=(spinor *) S + ix;
    
    (*r).s0.c0.re=fact0*((*r).s0.c0.re+fact*(*s).s0.c0.re);
    (*r).s0.c0.im=fact0*((*r).s0.c0.im+fact*(*s).s0.c0.im);
    (*r).s0.c1.re=fact0*((*r).s0.c1.re+fact*(*s).s0.c1.re);
    (*r).s0.c1.im=fact0*((*r).s0.c1.im+fact*(*s).s0.c1.im);
    (*r).s0.c2.re=fact0*((*r).s0.c2.re+fact*(*s).s0.c2.re);
    (*r).s0.c2.im=fact0*((*r).s0.c2.im+fact*(*s).s0.c2.im);
    
    (*r).s1.c0.re=fact0*((*r).s1.c0.re+fact*(*s).s1.c0.re);
    (*r).s1.c0.im=fact0*((*r).s1.c0.im+fact*(*s).s1.c0.im);
    (*r).s1.c1.re=fact0*((*r).s1.c1.re+fact*(*s).s1.c1.re);
    (*r).s1.c1.im=fact0*((*r).s1.c1.im+fact*(*s).s1.c1.im);
    (*r).s1.c2.re=fact0*((*r).s1.c2.re+fact*(*s).s1.c2.re);
    (*r).s1.c2.im=fact0*((*r).s1.c2.im+fact*(*s).s1.c2.im);
    
    (*r).s2.c0.re=fact0*((*r).s2.c0.re+fact*(*s).s2.c0.re);
    (*r).s2.c0.im=fact0*((*r).s2.c0.im+fact*(*s).s2.c0.im);
    (*r).s2.c1.re=fact0*((*r).s2.c1.re+fact*(*s).s2.c1.re);
    (*r).s2.c1.im=fact0*((*r).s2.c1.im+fact*(*s).s2.c1.im);
    (*r).s2.c2.re=fact0*((*r).s2.c2.re+fact*(*s).s2.c2.re);
    (*r).s2.c2.im=fact0*((*r).s2.c2.im+fact*(*s).s2.c2.im);
    
    (*r).s3.c0.re=fact0*((*r).s3.c0.re+fact*(*s).s3.c0.re);
    (*r).s3.c0.im=fact0*((*r).s3.c0.im+fact*(*s).s3.c0.im);
    (*r).s3.c1.re=fact0*((*r).s3.c1.re+fact*(*s).s3.c1.re);
    (*r).s3.c1.im=fact0*((*r).s3.c1.im+fact*(*s).s3.c1.im);
    (*r).s3.c2.re=fact0*((*r).s3.c2.re+fact*(*s).s3.c2.re);
    (*r).s3.c2.im=fact0*((*r).s3.c2.im+fact*(*s).s3.c2.im);
  }
}

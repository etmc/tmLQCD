/*******************************************************************************
 * $Id$ 
 * Makes (*R)=c1*(*R)+c2*(*S) , c1 and c2 are real constants 
 *******************************************************************************/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "su3.h"
#include "assign_mul_add_mul_r.h"


/* S input, R inoutput, c1,c2 input */
void assign_mul_add_mul_r(spinor * const R,spinor * const S, 
			  const double c1, const double c2,
			  const int N) {
  int ix;
  spinor *r,*s;
  
  for (ix = 0; ix < N; ix++){
    r=(spinor *) R + ix;
    s=(spinor *) S + ix;

    (*r).s0.c0.re=c1*(*r).s0.c0.re+c2*(*s).s0.c0.re;
    (*r).s0.c0.im=c1*(*r).s0.c0.im+c2*(*s).s0.c0.im;
    (*r).s0.c1.re=c1*(*r).s0.c1.re+c2*(*s).s0.c1.re;
    (*r).s0.c1.im=c1*(*r).s0.c1.im+c2*(*s).s0.c1.im;
    (*r).s0.c2.re=c1*(*r).s0.c2.re+c2*(*s).s0.c2.re;
    (*r).s0.c2.im=c1*(*r).s0.c2.im+c2*(*s).s0.c2.im;
    
    (*r).s1.c0.re=c1*(*r).s1.c0.re+c2*(*s).s1.c0.re;
    (*r).s1.c0.im=c1*(*r).s1.c0.im+c2*(*s).s1.c0.im;
    (*r).s1.c1.re=c1*(*r).s1.c1.re+c2*(*s).s1.c1.re;
    (*r).s1.c1.im=c1*(*r).s1.c1.im+c2*(*s).s1.c1.im;
    (*r).s1.c2.re=c1*(*r).s1.c2.re+c2*(*s).s1.c2.re;
    (*r).s1.c2.im=c1*(*r).s1.c2.im+c2*(*s).s1.c2.im;         
    
    (*r).s2.c0.re=c1*(*r).s2.c0.re+c2*(*s).s2.c0.re;
    (*r).s2.c0.im=c1*(*r).s2.c0.im+c2*(*s).s2.c0.im;
    (*r).s2.c1.re=c1*(*r).s2.c1.re+c2*(*s).s2.c1.re;
    (*r).s2.c1.im=c1*(*r).s2.c1.im+c2*(*s).s2.c1.im;
    (*r).s2.c2.re=c1*(*r).s2.c2.re+c2*(*s).s2.c2.re;
    (*r).s2.c2.im=c1*(*r).s2.c2.im+c2*(*s).s2.c2.im;         
    
    (*r).s3.c0.re=c1*(*r).s3.c0.re+c2*(*s).s3.c0.re;
    (*r).s3.c0.im=c1*(*r).s3.c0.im+c2*(*s).s3.c0.im;
    (*r).s3.c1.re=c1*(*r).s3.c1.re+c2*(*s).s3.c1.re;
    (*r).s3.c1.im=c1*(*r).s3.c1.im+c2*(*s).s3.c1.im;
    (*r).s3.c2.re=c1*(*r).s3.c2.re+c2*(*s).s3.c2.re;
    (*r).s3.c2.im=c1*(*r).s3.c2.im+c2*(*s).s3.c2.im;
  }
}






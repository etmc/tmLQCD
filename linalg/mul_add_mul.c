/* $Id$*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "su3.h"
#include "mul_add_mul.h"


/* Makes (*R)=c1*(*S)+c2*(*U) , c1 and c2 are complex constants */
void mul_add_mul(spinor * const R,spinor * const S,spinor * const U,const complex c1,const complex c2, const int N){
  int ix;
  spinor *r,*s,*u;
  
  for (ix=0;ix<N;ix++){
    r=(spinor *) R + ix;
    s=(spinor *) S + ix;
    u=(spinor *) U + ix;
    
    (*r).s0.c0.re=c1.re*(*s).s0.c0.re-c1.im*(*s).s0.c0.im+c2.re*(*u).s0.c0.re-c1.im*(*u).s0.c0.im;
    (*r).s0.c0.im=c1.re*(*s).s0.c0.im+c1.im*(*s).s0.c0.re+c2.re*(*u).s0.c0.im+c2.im*(*u).s0.c0.re;
    (*r).s0.c1.re=c1.re*(*s).s0.c1.re-c1.im*(*s).s0.c1.im+c2.re*(*u).s0.c1.re-c1.im*(*u).s0.c1.im;
    (*r).s0.c1.im=c1.re*(*s).s0.c1.im+c1.im*(*s).s0.c1.re+c2.re*(*u).s0.c1.im+c2.im*(*u).s0.c1.re;
    (*r).s0.c2.re=c1.re*(*s).s0.c2.re-c1.im*(*s).s0.c2.im+c2.re*(*u).s0.c2.re-c1.im*(*u).s0.c2.im;
    (*r).s0.c2.im=c1.re*(*s).s0.c2.im+c1.im*(*s).s0.c2.re+c2.re*(*u).s0.c2.im+c2.im*(*u).s0.c2.re;

    (*r).s1.c0.re=c1.re*(*s).s1.c0.re-c1.im*(*s).s1.c0.im+c2.re*(*u).s1.c0.re-c1.im*(*u).s1.c0.im;
    (*r).s1.c0.im=c1.re*(*s).s1.c0.im+c1.im*(*s).s1.c0.re+c2.re*(*u).s1.c0.im+c2.im*(*u).s1.c0.re;
    (*r).s1.c1.re=c1.re*(*s).s1.c1.re-c1.im*(*s).s1.c1.im+c2.re*(*u).s1.c1.re-c1.im*(*u).s1.c1.im;
    (*r).s1.c1.im=c1.re*(*s).s1.c1.im+c1.im*(*s).s1.c1.re+c2.re*(*u).s1.c1.im+c2.im*(*u).s1.c1.re;
    (*r).s1.c2.re=c1.re*(*s).s1.c2.re-c1.im*(*s).s1.c2.im+c2.re*(*u).s1.c2.re-c1.im*(*u).s1.c2.im;
    (*r).s1.c2.im=c1.re*(*s).s1.c2.im+c1.im*(*s).s1.c2.re+c2.re*(*u).s1.c2.im+c2.im*(*u).s1.c2.re;

    (*r).s2.c0.re=c1.re*(*s).s2.c0.re-c1.im*(*s).s2.c0.im+c2.re*(*u).s2.c0.re-c1.im*(*u).s2.c0.im;
    (*r).s2.c0.im=c1.re*(*s).s2.c0.im+c1.im*(*s).s2.c0.re+c2.re*(*u).s2.c0.im+c2.im*(*u).s2.c0.re;
    (*r).s2.c1.re=c1.re*(*s).s2.c1.re-c1.im*(*s).s2.c1.im+c2.re*(*u).s2.c1.re-c1.im*(*u).s2.c1.im;
    (*r).s2.c1.im=c1.re*(*s).s2.c1.im+c1.im*(*s).s2.c1.re+c2.re*(*u).s2.c1.im+c2.im*(*u).s2.c1.re;
    (*r).s2.c2.re=c1.re*(*s).s2.c2.re-c1.im*(*s).s2.c2.im+c2.re*(*u).s2.c2.re-c1.im*(*u).s2.c2.im;
    (*r).s2.c2.im=c1.re*(*s).s2.c2.im+c1.im*(*s).s2.c2.re+c2.re*(*u).s2.c2.im+c2.im*(*u).s2.c2.re;

    (*r).s3.c0.re=c1.re*(*s).s3.c0.re-c1.im*(*s).s3.c0.im+c2.re*(*u).s3.c0.re-c1.im*(*u).s3.c0.im;
    (*r).s3.c0.im=c1.re*(*s).s3.c0.im+c1.im*(*s).s3.c0.re+c2.re*(*u).s3.c0.im+c2.im*(*u).s3.c0.re;
    (*r).s3.c1.re=c1.re*(*s).s3.c1.re-c1.im*(*s).s3.c1.im+c2.re*(*u).s3.c1.re-c1.im*(*u).s3.c1.im;
    (*r).s3.c1.im=c1.re*(*s).s3.c1.im+c1.im*(*s).s3.c1.re+c2.re*(*u).s3.c1.im+c2.im*(*u).s3.c1.re;
    (*r).s3.c2.re=c1.re*(*s).s3.c2.re-c1.im*(*s).s3.c2.im+c2.re*(*u).s3.c2.re-c1.im*(*u).s3.c2.im;
    (*r).s3.c2.im=c1.re*(*s).s3.c2.im+c1.im*(*s).s3.c2.re+c2.re*(*u).s3.c2.im+c2.im*(*u).s3.c2.re;
  }
}




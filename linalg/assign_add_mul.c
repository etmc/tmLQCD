/*******************************************************************************
 * $Id$
 *
 * File assign_add_mul.c 
 *
 *   void assign_add_mul(spinor * const P, spinor * const Q, const complex c)
 *     (*P) = (*P) + c(*Q)        c is a complex constant
 *
 *******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "su3.h"
#include "assign_add_mul.h"


void assign_add_mul(spinor * const P, spinor * const Q, const complex c, const int N){
  int ix;
  static double a,b;
  spinor *r,*s;

  a=c.re;
  b=c.im;
   
  for (ix=0;ix<N;ix++){
    r=(spinor *) P + ix;
    s=(spinor *) Q + ix;

    (*r).s0.c0.re+=a*(*s).s0.c0.re - b*(*s).s0.c0.im;
    (*r).s0.c0.im+=a*(*s).s0.c0.im + b*(*s).s0.c0.re;
    (*r).s0.c1.re+=a*(*s).s0.c1.re - b*(*s).s0.c1.im;
    (*r).s0.c1.im+=a*(*s).s0.c1.im + b*(*s).s0.c1.re;
    (*r).s0.c2.re+=a*(*s).s0.c2.re - b*(*s).s0.c2.im;
    (*r).s0.c2.im+=a*(*s).s0.c2.im + b*(*s).s0.c2.re;

    (*r).s1.c0.re+=a*(*s).s1.c0.re - b*(*s).s1.c0.im;
    (*r).s1.c0.im+=a*(*s).s1.c0.im + b*(*s).s1.c0.re;
    (*r).s1.c1.re+=a*(*s).s1.c1.re - b*(*s).s1.c1.im;
    (*r).s1.c1.im+=a*(*s).s1.c1.im + b*(*s).s1.c1.re;
    (*r).s1.c2.re+=a*(*s).s1.c2.re - b*(*s).s1.c2.im;
    (*r).s1.c2.im+=a*(*s).s1.c2.im + b*(*s).s1.c2.re;

    (*r).s2.c0.re+=a*(*s).s2.c0.re - b*(*s).s2.c0.im;
    (*r).s2.c0.im+=a*(*s).s2.c0.im + b*(*s).s2.c0.re;
    (*r).s2.c1.re+=a*(*s).s2.c1.re - b*(*s).s2.c1.im;
    (*r).s2.c1.im+=a*(*s).s2.c1.im + b*(*s).s2.c1.re;
    (*r).s2.c2.re+=a*(*s).s2.c2.re - b*(*s).s2.c2.im;
    (*r).s2.c2.im+=a*(*s).s2.c2.im + b*(*s).s2.c2.re;

    (*r).s3.c0.re+=a*(*s).s3.c0.re - b*(*s).s3.c0.im;
    (*r).s3.c0.im+=a*(*s).s3.c0.im + b*(*s).s3.c0.re;
    (*r).s3.c1.re+=a*(*s).s3.c1.re - b*(*s).s3.c1.im;
    (*r).s3.c1.im+=a*(*s).s3.c1.im + b*(*s).s3.c1.re;
    (*r).s3.c2.re+=a*(*s).s3.c2.re - b*(*s).s3.c2.im;
    (*r).s3.c2.im+=a*(*s).s3.c2.im + b*(*s).s3.c2.re;

  }
}



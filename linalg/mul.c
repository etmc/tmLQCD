/*******************************************************************************
 * $Id$
 *
 * File mul.c
 *
 *   void mul(spinor * const R, const complex c, spinor * const S){
 *     Makes (*R) = c*(*S)
 *       
 *******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include "su3.h"
#include "mul.h"

void mul(spinor * const R, const complex c, spinor * const S, const int N){
  int ix;
  spinor *r,*s;
  
  for (ix=0;ix<N;ix++){
    r=(spinor *) R + ix;
    s=(spinor *) S + ix;
    
    (*r).s0.c0.re=c.re*(*s).s0.c0.re-c.im*(*s).s0.c0.im;
    (*r).s0.c0.im=c.re*(*s).s0.c0.im+c.im*(*s).s0.c0.re;
    (*r).s0.c1.re=c.re*(*s).s0.c1.re-c.im*(*s).s0.c1.im;
    (*r).s0.c1.im=c.re*(*s).s0.c1.im+c.im*(*s).s0.c1.re;
    (*r).s0.c2.re=c.re*(*s).s0.c2.re-c.im*(*s).s0.c2.im;
    (*r).s0.c2.im=c.re*(*s).s0.c2.im+c.im*(*s).s0.c2.re;
    
    (*r).s1.c0.re=c.re*(*s).s1.c0.re-c.im*(*s).s1.c0.im;
    (*r).s1.c0.im=c.re*(*s).s1.c0.im+c.im*(*s).s1.c0.re;
    (*r).s1.c1.re=c.re*(*s).s1.c1.re-c.im*(*s).s1.c1.im;
    (*r).s1.c1.im=c.re*(*s).s1.c1.im+c.im*(*s).s1.c1.re;
    (*r).s1.c2.re=c.re*(*s).s1.c2.re-c.im*(*s).s1.c2.im;
    (*r).s1.c2.im=c.re*(*s).s1.c2.im+c.im*(*s).s1.c2.re;

    (*r).s2.c0.re=c.re*(*s).s2.c0.re-c.im*(*s).s2.c0.im;
    (*r).s2.c0.im=c.re*(*s).s2.c0.im+c.im*(*s).s2.c0.re;
    (*r).s2.c1.re=c.re*(*s).s2.c1.re-c.im*(*s).s2.c1.im;
    (*r).s2.c1.im=c.re*(*s).s2.c1.im+c.im*(*s).s2.c1.re;
    (*r).s2.c2.re=c.re*(*s).s2.c2.re-c.im*(*s).s2.c2.im;
    (*r).s2.c2.im=c.re*(*s).s2.c2.im+c.im*(*s).s2.c2.re;

    (*r).s3.c0.re=c.re*(*s).s3.c0.re-c.im*(*s).s3.c0.im;
    (*r).s3.c0.im=c.re*(*s).s3.c0.im+c.im*(*s).s3.c0.re;
    (*r).s3.c1.re=c.re*(*s).s3.c1.re-c.im*(*s).s3.c1.im;
    (*r).s3.c1.im=c.re*(*s).s3.c1.im+c.im*(*s).s3.c1.re;
    (*r).s3.c2.re=c.re*(*s).s3.c2.re-c.im*(*s).s3.c2.im;
    (*r).s3.c2.im=c.re*(*s).s3.c2.im+c.im*(*s).s3.c2.re;
  }
}

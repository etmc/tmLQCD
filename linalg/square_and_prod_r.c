/*******************************************************************************
 * $Id$
 *
 * File square_and_prod_r.c
 *
 *   void square_and_prod_r(double * const x1, double * const x2, spinor * const S, spinor * const R)
 *     Returns the real part of (*R,*S) and the square norm of *S
 *     It's faster than using "scalar_prod_r" and "square_norm"
 *       
 *******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#ifdef MPI
# include <mpi.h>
#endif
#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include "su3.h"
#include "square_and_prod_r.h"

void square_and_prod_r(double * const x1, double * const x2, spinor * const S, spinor * const R, const int N){
  int ix;
  static double ks,kc,ds,tr,ts,tt;
  static double xks,xkc,xds,xtr,xts,xtt;
  spinor *s,*r;
  
  ks=0.0;
  kc=0.0;
  
  xks=0.0;
  xkc=0.0;
  
  for (ix = 0; ix < N; ix++){
    s=(spinor *) S + ix;
    r=(spinor *) R + ix;

    ds=(*r).s0.c0.re*(*s).s0.c0.re + (*r).s0.c0.im*(*s).s0.c0.im + 
       (*r).s0.c1.re*(*s).s0.c1.re + (*r).s0.c1.im*(*s).s0.c1.im + 
       (*r).s0.c2.re*(*s).s0.c2.re + (*r).s0.c2.im*(*s).s0.c2.im + 
       (*r).s1.c0.re*(*s).s1.c0.re + (*r).s1.c0.im*(*s).s1.c0.im + 
       (*r).s1.c1.re*(*s).s1.c1.re + (*r).s1.c1.im*(*s).s1.c1.im + 
       (*r).s1.c2.re*(*s).s1.c2.re + (*r).s1.c2.im*(*s).s1.c2.im + 
       (*r).s2.c0.re*(*s).s2.c0.re + (*r).s2.c0.im*(*s).s2.c0.im + 
       (*r).s2.c1.re*(*s).s2.c1.re + (*r).s2.c1.im*(*s).s2.c1.im + 
       (*r).s2.c2.re*(*s).s2.c2.re + (*r).s2.c2.im*(*s).s2.c2.im + 
       (*r).s3.c0.re*(*s).s3.c0.re + (*r).s3.c0.im*(*s).s3.c0.im + 
       (*r).s3.c1.re*(*s).s3.c1.re + (*r).s3.c1.im*(*s).s3.c1.im + 
       (*r).s3.c2.re*(*s).s3.c2.re + (*r).s3.c2.im*(*s).s3.c2.im;
     
    xds=(*s).s0.c0.re*(*s).s0.c0.re + (*s).s0.c0.im*(*s).s0.c0.im + 
        (*s).s0.c1.re*(*s).s0.c1.re + (*s).s0.c1.im*(*s).s0.c1.im + 
        (*s).s0.c2.re*(*s).s0.c2.re + (*s).s0.c2.im*(*s).s0.c2.im + 
        (*s).s1.c0.re*(*s).s1.c0.re + (*s).s1.c0.im*(*s).s1.c0.im + 
        (*s).s1.c1.re*(*s).s1.c1.re + (*s).s1.c1.im*(*s).s1.c1.im + 
        (*s).s1.c2.re*(*s).s1.c2.re + (*s).s1.c2.im*(*s).s1.c2.im + 
        (*s).s2.c0.re*(*s).s2.c0.re + (*s).s2.c0.im*(*s).s2.c0.im + 
        (*s).s2.c1.re*(*s).s2.c1.re + (*s).s2.c1.im*(*s).s2.c1.im + 
        (*s).s2.c2.re*(*s).s2.c2.re + (*s).s2.c2.im*(*s).s2.c2.im + 
        (*s).s3.c0.re*(*s).s3.c0.re + (*s).s3.c0.im*(*s).s3.c0.im + 
        (*s).s3.c1.re*(*s).s3.c1.re + (*s).s3.c1.im*(*s).s3.c1.im + 
        (*s).s3.c2.re*(*s).s3.c2.re + (*s).s3.c2.im*(*s).s3.c2.im;
    
    tr=ds + kc;
    ts=tr + ks;
    tt=ts-ks;
    ks=ts;
    kc=tr-tt;
    
    xtr=xds + xkc;
    xts=xtr + xks;
    xtt=xts-xks;
    xks=xts;
    xkc=xtr-xtt;
  }
  xkc=xks + xkc;
  *x1=xkc;

#if defined MPI

  MPI_Allreduce(&xkc, x1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

#endif
  kc=ks + kc;
  *x2=kc;

#if defined MPI

    MPI_Allreduce(&kc, x2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

#endif
}


/*******************************************************************************
 * $Id$
 *
 * File scalar_prod_r.c
 *
 *   double scalar_prod_r(spinor * const S,spinor * const R, const int N)
 *     Returns the real part of the scalar product (*R,*S)
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
#ifdef _STD_C99_COMPLEX_CHECKED
# include <complex.h>
#endif
#ifdef apenext
# include <topology.h>
# include <queue.h>
#endif
#include "su3.h"
#include "scalar_prod_r.h"

/*  R input, S input */

#if ((!defined _STD_C99_COMPLEX_CHECKED) && (!defined apenext))

double scalar_prod_r(spinor * const S,spinor * const R, const int N){
  int ix;
  static double ks,kc,ds,tr,ts,tt;
  spinor *s,*r;
  ks=0.0;
  kc=0.0;
  
  for (ix=0;ix<N;ix++){
    s = (spinor *) S + ix;
    r = (spinor *) R + ix;
    
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
    
    tr = ds + kc;
    ts = tr + ks;
    tt = ts-ks;
    ks = ts;
    kc = tr-tt;
  }
  kc = ks + kc;

#if defined MPI
  
  MPI_Allreduce(&kc, &ks, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return ks;

#endif

  return kc;

}
#endif

#if ((defined _STD_C99_COMPLEX_CHECKED) && (!defined apenext))

double scalar_prod_r(spinor * const S,spinor * const R, const int N){
  register int ix=0;
  register complex ds00,ds01,ds02,ds10,ds11,ds12,ds20,ds21,ds22,ds30,ds31,ds32;
  register double kc;
  register spinor *sPointer, *rPointer;

  sPointer=S;
  rPointer=R;
  ds00=0.; ds01=0.; ds02=0.;
  ds10=0.; ds11=0.; ds12=0.;
  ds20=0.; ds21=0.; ds22=0.;
  ds30=0.; ds31=0.; ds32=0.;

  do {
    register spinor s, r;
    ix+=1;
    
    s = *(sPointer);
    r = *(rPointer);
    ds00+=conj(r.s0.c0)*s.s0.c0;
    ds01+=conj(r.s0.c1)*s.s0.c1;
    ds02+=conj(r.s0.c2)*s.s0.c2;
    ds10+=conj(r.s1.c0)*s.s1.c0;
    ds11+=conj(r.s1.c1)*s.s1.c1;
    ds12+=conj(r.s1.c2)*s.s1.c2;
    ds20+=conj(r.s2.c0)*s.s2.c0;
    ds21+=conj(r.s2.c1)*s.s2.c1;
    ds22+=conj(r.s2.c2)*s.s2.c2;
    ds30+=conj(r.s3.c0)*s.s3.c0;
    ds31+=conj(r.s3.c1)*s.s3.c1;
    ds32+=conj(r.s3.c2)*s.s3.c2;

    sPointer+=1;
    rPointer+=1;
      
  } while (ix<N);

 ds00=ds00+ds01;
 ds02=ds02+ds10;
 ds11=ds11+ds12;
 ds20=ds20+ds21;
 ds22=ds22+ds30;
 ds31=ds31+ds32;
 
 ds00=ds00+ds02;
 ds01=ds11+ds20;
 ds02=ds22+ds31;
 
 ds10=ds00+ds01;
 kc=creal(ds10+ds02);

#if defined MPI
  
  MPI_Allreduce(&kc, &ks, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return ks;

#endif

  return kc;

}
#endif


#ifdef apenext

#define NOWHERE_COND(condition) ((condition) ? 0x0 : NOWHERE )

double scalar_prod_r(spinor * const S,spinor * const R, const int N){
  register int ix=N;
  register complex ds00,ds01,ds02,ds10,ds11,ds12,ds20,ds21,ds22,ds30,ds31,ds32;
  register double kc;
  register spinor *sPointer, *rPointer;

  sPointer=S;
  rPointer=R;
  ds00=0.; ds01=0.; ds02=0.;
  ds10=0.; ds11=0.; ds12=0.;
  ds20=0.; ds21=0.; ds22=0.;
  ds30=0.; ds31=0.; ds32=0.;

  prefetch (*(sPointer));
  prefetch (*(rPointer));
  {
#pragma cache

    do {
      register spinor s, r;
      ix-=1;
      
      sPointer+=1;
      rPointer+=1;
      fetch(s);
      fetch(r);      
      prefetch (*(sPointer+NOWHERE_COND(ix)));
      prefetch (*(rPointer+NOWHERE_COND(ix)));

      ds00+=conj(r.s0.c0)*s.s0.c0;
      ds01+=conj(r.s0.c1)*s.s0.c1;
      ds02+=conj(r.s0.c2)*s.s0.c2;
      ds10+=conj(r.s1.c0)*s.s1.c0;
      ds11+=conj(r.s1.c1)*s.s1.c1;
      ds12+=conj(r.s1.c2)*s.s1.c2;
      ds20+=conj(r.s2.c0)*s.s2.c0;
      ds21+=conj(r.s2.c1)*s.s2.c1;
      ds22+=conj(r.s2.c2)*s.s2.c2;
      ds30+=conj(r.s3.c0)*s.s3.c0;
      ds31+=conj(r.s3.c1)*s.s3.c1;
      ds32+=conj(r.s3.c2)*s.s3.c2;

    } while (ix>0);
  }
  ds00=ds00+ds01;
  ds02=ds02+ds10;
  ds11=ds11+ds12;
  ds20=ds20+ds21;
  ds22=ds22+ds30;
  ds31=ds31+ds32;
 
  ds00=ds00+ds02;
  ds01=ds11+ds20;
  ds02=ds22+ds31;
  
  ds10=ds00+ds01;
  kc=creal(ds10+ds02);

#if defined MPI
  
  MPI_Allreduce(&kc, &ks, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return ks;

#endif

  return kc;

}
#endif

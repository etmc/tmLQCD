/*******************************************************************************
 * $Id$
 *
 * File square_norm.c
 *
 *   double square_norm(spinor * const P )
 *     Returns the square norm of *P
 *
 *******************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#ifdef MPI
#include <mpi.h>
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
#include "sse.h"
#include "square_norm.h"

#if ((!defined _STD_C99_COMPLEX_CHECKED) && (!defined apenext))

double square_norm(spinor * const P, const int N) {
  int ix;
  static double ks,kc,ds,tr,ts,tt;
  spinor *s;
  
  ks = 0.0;
  kc = 0.0;
  
  /* Change due to even-odd preconditioning : VOLUME   to VOLUME/2 */   
  for (ix  =  0; ix < N; ix++) {
    s = P + ix;
    
    ds = (*s).s0.c0.re*(*s).s0.c0.re + (*s).s0.c0.im*(*s).s0.c0.im + 
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
    
    tr = ds + kc;
    ts = tr + ks;
    tt = ts-ks;
    ks = ts;
    kc = tr-tt;
  }
  kc = ks + kc;
#ifdef MPI
  MPI_Allreduce(&kc, &ks, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return ks;
#else
  return kc;
#endif
}
#endif

#if ((defined _STD_C99_COMPLEX_CHECKED) && (!defined apenext))

double square_norm(spinor * const P, const int N) {

  register int ix=0;
  register complex ds00,ds01,ds02,ds10,ds11,ds12,ds20,ds21,ds22,ds30,ds31,ds32;
  register double kc;
  register spinor *sPointer;

  sPointer=P;

  ds00=0.; ds01=0.; ds02=0.;
  ds10=0.; ds11=0.; ds12=0.;
  ds20=0.; ds21=0.; ds22=0.;
  ds30=0.; ds31=0.; ds32=0.;

  do {
    register spinor s;
    ix+=1;
    
    s = *(sPointer);

    ds00+=conj(s.s0.c0)*s.s0.c0;
    ds01+=conj(s.s0.c1)*s.s0.c1;
    ds02+=conj(s.s0.c2)*s.s0.c2;
    ds10+=conj(s.s1.c0)*s.s1.c0;
    ds11+=conj(s.s1.c1)*s.s1.c1;
    ds12+=conj(s.s1.c2)*s.s1.c2;
    ds20+=conj(s.s2.c0)*s.s2.c0;
    ds21+=conj(s.s2.c1)*s.s2.c1;
    ds22+=conj(s.s2.c2)*s.s2.c2;
    ds30+=conj(s.s3.c0)*s.s3.c0;
    ds31+=conj(s.s3.c1)*s.s3.c1;
    ds32+=conj(s.s3.c2)*s.s3.c2;

    sPointer+=1;
      
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

double square_norm(spinor * const P, const int N) {
  register int ix=N;
  register complex ds00,ds01,ds02,ds10,ds11,ds12,ds20,ds21,ds22,ds30,ds31,ds32;
  register double kc;
  register spinor *sPointer;

  sPointer=P;

  ds00=0.; ds01=0.; ds02=0.;
  ds10=0.; ds11=0.; ds12=0.;
  ds20=0.; ds21=0.; ds22=0.;
  ds30=0.; ds31=0.; ds32=0.;

  prefetch (*(sPointer));

  {
#pragma cache

    do {
      register spinor s;
      ix-=1;
      
      sPointer+=1;

      fetch(s);

      prefetch (*(sPointer+NOWHERE_COND(ix)));


      ds00+=conj(s.s0.c0)*s.s0.c0;
      ds01+=conj(s.s0.c1)*s.s0.c1;
      ds02+=conj(s.s0.c2)*s.s0.c2;
      ds10+=conj(s.s1.c0)*s.s1.c0;
      ds11+=conj(s.s1.c1)*s.s1.c1;
      ds12+=conj(s.s1.c2)*s.s1.c2;
      ds20+=conj(s.s2.c0)*s.s2.c0;
      ds21+=conj(s.s2.c1)*s.s2.c1;
      ds22+=conj(s.s2.c2)*s.s2.c2;
      ds30+=conj(s.s3.c0)*s.s3.c0;
      ds31+=conj(s.s3.c1)*s.s3.c1;
      ds32+=conj(s.s3.c2)*s.s3.c2;

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

/*******************************************************************************
 * $Id$
 *
 * File square_norm.c
 *
 *   double square_norm(spinor * const P )
 *     Returns the square norm of *P
 *
 *******************************************************************************/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#ifdef MPI
# include <mpi.h>
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

#if ((!defined _STD_C99_COMPLEX_CHECKED) && (!defined apenext) && (!defined BGL))

double square_norm(spinor * const P, const int N, const int parallel) {
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
#  ifdef MPI
  if(parallel) {
    MPI_Allreduce(&kc, &ks, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return ks;
  }
#endif
  return kc;
}

#elif ((defined BGL) && (defined XLC))

/***************************************
 *
 * square norm with intrinsics
 *
 * Carsten.Urbach@liverpool.ac.uk
 *
 ***************************************/

#  include"bgl.h"
double square_norm(spinor * const P, const int N, const int parallel) {
  int ix=0;
  double res, res2;
  double *s ALIGN;
  double *sp ALIGN;
  double _Complex x00, x01, x02, x03, x04, x05, x06, x07, 
    x08, x09, x10, x11;
  double _Complex y00, y01, y02, y03, y04, y05, y06, y07, 
    y08, y09, y10, y11;

  __alignx(16, P);
  s = (double*)P;
  sp = s+24;
  _prefetch_spinor(sp);
  x00 = __lfpd(s);
  x01 = __lfpd(s+2);
  x02 = __lfpd(s+4);
  x03 = __lfpd(s+6);
  x04 = __lfpd(s+8);
  x05 = __lfpd(s+10);
  x06 = __lfpd(s+12);
  x07 = __lfpd(s+14);
  x08 = __lfpd(s+16);
  x09 = __lfpd(s+18);
  x10 = __lfpd(s+20);
  x11 = __lfpd(s+22);
  
  y00 = __fpmul(x00, x00);
  y01 = __fpmul(x01, x01);
  y02 = __fpmul(x02, x02);
  y03 = __fpmul(x03, x03);
  y04 = __fpmul(x04, x04);
  y05 = __fpmul(x05, x05);
  y06 = __fpmul(x06, x06);
  y07 = __fpmul(x07, x07);
  y08 = __fpmul(x08, x08);
  y09 = __fpmul(x09, x09);
  y10 = __fpmul(x10, x10);
  y11 = __fpmul(x11, x11);
  s = sp;


#pragma unroll(12)
  for(ix = 1; ix < N-1; ix++) {
    sp+=24;;
    _prefetch_spinor(sp);
    x00 = __lfpd(s);   
    x01 = __lfpd(s+2); 
    x02 = __lfpd(s+4); 
    x03 = __lfpd(s+6); 
    x04 = __lfpd(s+8); 
    x05 = __lfpd(s+10);
    x06 = __lfpd(s+12);
    x07 = __lfpd(s+14);
    x08 = __lfpd(s+16);
    x09 = __lfpd(s+18);
    x10 = __lfpd(s+20);
    x11 = __lfpd(s+22);
    y00 = __fpmadd(y00, x00, x00); 
    y01 = __fpmadd(y01, x01, x01); 
    y02 = __fpmadd(y02, x02, x02); 
    y03 = __fpmadd(y03, x03, x03); 
    y04 = __fpmadd(y04, x04, x04); 
    y05 = __fpmadd(y05, x05, x05); 
    y06 = __fpmadd(y06, x06, x06); 
    y07 = __fpmadd(y07, x07, x07); 
    y08 = __fpmadd(y08, x08, x08); 
    y09 = __fpmadd(y09, x09, x09); 
    y10 = __fpmadd(y10, x10, x10); 
    y11 = __fpmadd(y11, x11, x11); 
    s=sp;
  }
  x00 = __lfpd(s);   
  x01 = __lfpd(s+2); 
  x02 = __lfpd(s+4); 
  x03 = __lfpd(s+6); 
  x04 = __lfpd(s+8); 
  x05 = __lfpd(s+10);
  x06 = __lfpd(s+12);
  x07 = __lfpd(s+14);
  x08 = __lfpd(s+16);
  x09 = __lfpd(s+18);
  x10 = __lfpd(s+20);
  x11 = __lfpd(s+22);
  y00 = __fpmadd(y00, x00, x00); 
  y01 = __fpmadd(y01, x01, x01); 
  y02 = __fpmadd(y02, x02, x02); 
  y03 = __fpmadd(y03, x03, x03); 
  y04 = __fpmadd(y04, x04, x04); 
  y05 = __fpmadd(y05, x05, x05); 
  y06 = __fpmadd(y06, x06, x06); 
  y07 = __fpmadd(y07, x07, x07); 
  y08 = __fpmadd(y08, x08, x08); 
  y09 = __fpmadd(y09, x09, x09); 
  y10 = __fpmadd(y10, x10, x10); 
  y11 = __fpmadd(y11, x11, x11); 
  
  y00 = __fpadd(y00, y01);
  y02 = __fpadd(y02, y03);
  y04 = __fpadd(y04, y05);
  y06 = __fpadd(y06, y07);
  y08 = __fpadd(y08, y09);
  y10 = __fpadd(y10, y11);
  y00 = __fpadd(y00, y02);
  y04 = __fpadd(y04, y06);
  y08 = __fpadd(y08, y10);
  y00 = __fpadd(y00, y04);
  y00 = __fpadd(y00, y08);
  res = __creal(y00)+__cimag(y00);
#  ifdef MPI
  if(parallel) {
    MPI_Allreduce(&res, &res2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return res2;
  }
#  endif
  return res;
}

#elif ((defined _STD_C99_COMPLEX_CHECKED) && (!defined apenext))

double square_norm(spinor * const P, const int N, const int parallel) {

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
 if(parallel) {
   MPI_Allreduce(&kc, &ks, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   return ks;
 }
#endif
  return kc;
}
#elif defined apenext

#define NOWHERE_COND(condition) ((condition) ? 0x0 : NOWHERE )

double square_norm(spinor * const P, const int N, const int parallel) {
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
  if(parallel) {
    MPI_Allreduce(&kc, &ks, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return ks;
  }
#endif
  return kc;
}
#endif

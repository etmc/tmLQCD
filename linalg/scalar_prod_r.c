/***********************************************************************
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008 Carsten Urbach
 *
 * This file is part of tmLQCD.
 *
 * tmLQCD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * tmLQCD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with tmLQCD.  If not, see <http://www.gnu.org/licenses/>.
 ***********************************************************************/
/*******************************************************************************
 *
 * File scalar_prod_r.c
 *
 *   double scalar_prod_r(spinor * const S,spinor * const R, const int N)
 *     Returns the real part of the scalar product (*R,*S)
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
#include "scalar_prod_r.h"

/*  R input, S input */

#if ((!defined _STD_C99_COMPLEX_CHECKED) && (!defined apenext))

double scalar_prod_r(spinor * const S,spinor * const R, const int N, const int parallel){
  int ix;
  static double ks,kc,ds,tr,ts,tt;
  spinor *s,*r;
  ks=0.0;
  kc=0.0;

#if (defined BGL && defined XLC)
  __alignx(16, S);
  __alignx(16, R);
#endif
  
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
  if(parallel) {
    MPI_Allreduce(&kc, &ks, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return ks;
  }
#endif

  return kc;

}
#endif

#if ((defined _STD_C99_COMPLEX_CHECKED) && (!defined apenext))

double scalar_prod_r(spinor * const S,spinor * const R, const int N, const int parallel){
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
 if(parallel) {
   MPI_Allreduce(&kc, &ks, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   return ks;
 }

#endif

  return kc;

}
#endif


#ifdef apenext

#define NOWHERE_COND(condition) ((condition) ? 0x0 : NOWHERE )

double scalar_prod_r(spinor * const S,spinor * const R, const int N, const int parallel){
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
  if(parallel) {
    MPI_Allreduce(&kc, &ks, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return ks;
  }
#endif

  return kc;

}
#endif // apenext

#ifdef WITHLAPH
double scalar_prod_r_su3vect(su3_vector * const S,su3_vector * const R, const int N, const int parallel)
{
  int ix;
  static double ks,kc,ds,tr,ts,tt;
  su3_vector *s,*r;

  ks=0.0;
  kc=0.0;
  for (ix=0;ix<N;ix++)
    {
      s = (su3_vector *) S + ix;
      r = (su3_vector *) R + ix;
    
      ds=(*r).c0.re*(*s).c0.re + (*r).c0.im*(*s).c0.im + 
	(*r).c1.re*(*s).c1.re + (*r).c1.im*(*s).c1.im + 
	(*r).c2.re*(*s).c2.re + (*r).c2.im*(*s).c2.im;
    
      tr = ds + kc;
      ts = tr + ks;
      tt = ts-ks;
      ks = ts;
      kc = tr-tt;
    }
  kc = ks + kc;
#if defined MPI
  if(parallel) {
    MPI_Allreduce(&kc, &ks, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return ks;
  }
#endif

  return kc;
}

#endif // WITHLAPH


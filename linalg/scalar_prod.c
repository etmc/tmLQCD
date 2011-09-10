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
/* $Id$ */

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#ifdef MPI
#include <mpi.h>
#endif
#include "su3.h"
#include "scalar_prod.h"

/*  <S,R>=S^* times R */
complex scalar_prod(spinor * const S, spinor * const R, const int N, const int parallel){
  int ix;
  static double ks,kc,ds,tr,ts,tt;
  spinor *s,*r;
  complex c;
#ifdef MPI
  complex d;
#endif
  
  /* Real Part */

  ks=0.0;
  kc=0.0;
#if (defined BGL && defined XLC)
  __alignx(16, S);
  __alignx(16, R);
#endif  
  for (ix = 0; ix < N; ix++){
    s=(spinor *) S + ix;
    r=(spinor *) R + ix;
    
    ds=(*r).s0.c0.re*(*s).s0.c0.re+(*r).s0.c0.im*(*s).s0.c0.im+
      (*r).s0.c1.re*(*s).s0.c1.re+(*r).s0.c1.im*(*s).s0.c1.im+
      (*r).s0.c2.re*(*s).s0.c2.re+(*r).s0.c2.im*(*s).s0.c2.im+
      (*r).s1.c0.re*(*s).s1.c0.re+(*r).s1.c0.im*(*s).s1.c0.im+
      (*r).s1.c1.re*(*s).s1.c1.re+(*r).s1.c1.im*(*s).s1.c1.im+
      (*r).s1.c2.re*(*s).s1.c2.re+(*r).s1.c2.im*(*s).s1.c2.im+
      (*r).s2.c0.re*(*s).s2.c0.re+(*r).s2.c0.im*(*s).s2.c0.im+
      (*r).s2.c1.re*(*s).s2.c1.re+(*r).s2.c1.im*(*s).s2.c1.im+
      (*r).s2.c2.re*(*s).s2.c2.re+(*r).s2.c2.im*(*s).s2.c2.im+
      (*r).s3.c0.re*(*s).s3.c0.re+(*r).s3.c0.im*(*s).s3.c0.im+
      (*r).s3.c1.re*(*s).s3.c1.re+(*r).s3.c1.im*(*s).s3.c1.im+
      (*r).s3.c2.re*(*s).s3.c2.re+(*r).s3.c2.im*(*s).s3.c2.im;

    /* Kahan Summation */    
    tr=ds+kc;
    ts=tr+ks;
    tt=ts-ks;
    ks=ts;
    kc=tr-tt;
  }
  kc=ks+kc;

#if defined MPI0
  if(parallel == 1) {
    MPI_Allreduce(&kc, &ks, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    kc = ks;
  }
#endif

  c.re = kc;

  /* Imaginary Part */

  ks=0.0;
  kc=0.0;
  
  for (ix=0;ix<N;ix++){
    s=(spinor *) S + ix;
    r=(spinor *) R + ix;
    
    ds=-(*r).s0.c0.re*(*s).s0.c0.im+(*r).s0.c0.im*(*s).s0.c0.re-
      (*r).s0.c1.re*(*s).s0.c1.im+(*r).s0.c1.im*(*s).s0.c1.re-
      (*r).s0.c2.re*(*s).s0.c2.im+(*r).s0.c2.im*(*s).s0.c2.re-
      (*r).s1.c0.re*(*s).s1.c0.im+(*r).s1.c0.im*(*s).s1.c0.re-
      (*r).s1.c1.re*(*s).s1.c1.im+(*r).s1.c1.im*(*s).s1.c1.re-
      (*r).s1.c2.re*(*s).s1.c2.im+(*r).s1.c2.im*(*s).s1.c2.re-
      (*r).s2.c0.re*(*s).s2.c0.im+(*r).s2.c0.im*(*s).s2.c0.re-
      (*r).s2.c1.re*(*s).s2.c1.im+(*r).s2.c1.im*(*s).s2.c1.re-
      (*r).s2.c2.re*(*s).s2.c2.im+(*r).s2.c2.im*(*s).s2.c2.re-
      (*r).s3.c0.re*(*s).s3.c0.im+(*r).s3.c0.im*(*s).s3.c0.re-
      (*r).s3.c1.re*(*s).s3.c1.im+(*r).s3.c1.im*(*s).s3.c1.re-
      (*r).s3.c2.re*(*s).s3.c2.im+(*r).s3.c2.im*(*s).s3.c2.re;
    
    tr=ds+kc;
    ts=tr+ks;
    tt=ts-ks;
    ks=ts;
    kc=tr-tt;
  }
  kc=ks+kc;

#if defined MPI0
  if(parallel == 1) {
    MPI_Allreduce(&kc, &ks, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    kc = ks;
  }
#endif

  c.im = kc;
#ifdef MPI
  if(parallel == 1) {
    d = c;
    MPI_Allreduce(&d, &c, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
  }
#endif
  return(c);
}

#ifdef WITHLAPH
complex scalar_prod_su3vect(su3_vector * const S, su3_vector * const R, const int N, const int parallel){
  int ix;
  static double ks,kc,ds,tr,ts,tt;
  su3_vector *s,*r;
  complex c;
#ifdef MPI
  complex d;
#endif

  /* Real Part */

  ks=0.0;
  kc=0.0;
  for (ix = 0; ix < N; ix++)
    {
      s=(su3_vector *) S + ix;
      r=(su3_vector *) R + ix;
    
      ds=(*r).c0.re*(*s).c0.re+(*r).c0.im*(*s).c0.im+
	(*r).c1.re*(*s).c1.re+(*r).c1.im*(*s).c1.im+
	(*r).c2.re*(*s).c2.re+(*r).c2.im*(*s).c2.im;

      /* Kahan Summation */    
      tr=ds+kc;
      ts=tr+ks;
      tt=ts-ks;
      ks=ts;
      kc=tr-tt;
    }
  kc=ks+kc;

#if defined MPI0
  if(parallel == 1) {
    MPI_Allreduce(&kc, &ks, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    kc = ks;
  }
#endif

  c.re = kc;

  /* Imaginary Part */

  ks=0.0;
  kc=0.0;
  
  for (ix=0;ix<N;ix++)
    {
      s=(su3_vector *) S + ix;
      r=(su3_vector *) R + ix;
    
      ds=-(*r).c0.re*(*s).c0.im+(*r).c0.im*(*s).c0.re-
	(*r).c1.re*(*s).c1.im+(*r).c1.im*(*s).c1.re-
	(*r).c2.re*(*s).c2.im+(*r).c2.im*(*s).c2.re;
    
      tr=ds+kc;
      ts=tr+ks;
      tt=ts-ks;
      ks=ts;
      kc=tr-tt;
    }
  kc=ks+kc;

#if defined MPI0
  if(parallel == 1) {
    MPI_Allreduce(&kc, &ks, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    kc = ks;
  }
#endif

  c.im = kc;
#ifdef MPI
  if(parallel == 1) {
    d = c;
    MPI_Allreduce(&d, &c, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
  }
#endif
  return(c);
}
#endif

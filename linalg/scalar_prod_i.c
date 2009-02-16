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
 * $Id$
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
#include "scalar_prod_i.h"

/*  R input, S input */


double scalar_prod_i(spinor * const S,spinor * const R, const int N, const int parallel){
  int ix;
  static double ks,kc,ds,tr,ts,tt;
  spinor *s,*r;
  ks=0.0;
  kc=0.0;

#if (defined BGL && defined XLC)
  __alignx(16, S);
  __alignx(16, R);
#endif
  
  for (ix = 0; ix < N; ix++){
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
  
#if defined MPI
  if(parallel == 1) {
    MPI_Allreduce(&kc, &ks, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    kc = ks;
  }
#endif

  return kc;

}



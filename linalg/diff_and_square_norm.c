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
# include <mpi.h>
#endif
#include "su3.h"
#include "diff_and_square_norm.h"

double diff_and_square_norm(spinor * const Q, spinor * const R, const int N) {
  int ix;
  static double ks,kc,ds,tr,ts,tt;
  spinor *q,*r;
  
  ks=0.0;
  kc=0.0;
  
  /* Change due to even-odd preconditioning : VOLUME   to VOLUME/2 */   
  for (ix = 0; ix < N; ix++) {
    q=Q+ix;
    r=R+ix;
    
    (*q).s0.c0.re = (*r).s0.c0.re-(*q).s0.c0.re;
    (*q).s0.c0.im = (*r).s0.c0.im-(*q).s0.c0.im;
    (*q).s0.c1.re = (*r).s0.c1.re-(*q).s0.c1.re;
    (*q).s0.c1.im = (*r).s0.c1.im-(*q).s0.c1.im;
    (*q).s0.c2.re = (*r).s0.c2.re-(*q).s0.c2.re;
    (*q).s0.c2.im = (*r).s0.c2.im-(*q).s0.c2.im;
    
    ds = 
      (*q).s0.c0.re*(*q).s0.c0.re+(*q).s0.c0.im*(*q).s0.c0.im+
      (*q).s0.c1.re*(*q).s0.c1.re+(*q).s0.c1.im*(*q).s0.c1.im+
      (*q).s0.c2.re*(*q).s0.c2.re+(*q).s0.c2.im*(*q).s0.c2.im;
    
    (*q).s1.c0.re = (*r).s1.c0.re-(*q).s1.c0.re;
    (*q).s1.c0.im = (*r).s1.c0.im-(*q).s1.c0.im;
    (*q).s1.c1.re = (*r).s1.c1.re-(*q).s1.c1.re;
    (*q).s1.c1.im = (*r).s1.c1.im-(*q).s1.c1.im;
    (*q).s1.c2.re = (*r).s1.c2.re-(*q).s1.c2.re;
    (*q).s1.c2.im = (*r).s1.c2.im-(*q).s1.c2.im;         
    
    ds+= 
      (*q).s1.c0.re*(*q).s1.c0.re+(*q).s1.c0.im*(*q).s1.c0.im+
      (*q).s1.c1.re*(*q).s1.c1.re+(*q).s1.c1.im*(*q).s1.c1.im+
      (*q).s1.c2.re*(*q).s1.c2.re+(*q).s1.c2.im*(*q).s1.c2.im;
    
    (*q).s2.c0.re = (*r).s2.c0.re-(*q).s2.c0.re;
    (*q).s2.c0.im = (*r).s2.c0.im-(*q).s2.c0.im;
    (*q).s2.c1.re = (*r).s2.c1.re-(*q).s2.c1.re;
    (*q).s2.c1.im = (*r).s2.c1.im-(*q).s2.c1.im;
    (*q).s2.c2.re = (*r).s2.c2.re-(*q).s2.c2.re;
    (*q).s2.c2.im = (*r).s2.c2.im-(*q).s2.c2.im;         
    
    ds+= 
      (*q).s2.c0.re*(*q).s2.c0.re+(*q).s2.c0.im*(*q).s2.c0.im+
      (*q).s2.c1.re*(*q).s2.c1.re+(*q).s2.c1.im*(*q).s2.c1.im+
      (*q).s2.c2.re*(*q).s2.c2.re+(*q).s2.c2.im*(*q).s2.c2.im;
    
    (*q).s3.c0.re = (*r).s3.c0.re-(*q).s3.c0.re;
    (*q).s3.c0.im = (*r).s3.c0.im-(*q).s3.c0.im;
    (*q).s3.c1.re = (*r).s3.c1.re-(*q).s3.c1.re;
    (*q).s3.c1.im = (*r).s3.c1.im-(*q).s3.c1.im;
    (*q).s3.c2.re = (*r).s3.c2.re-(*q).s3.c2.re;
    (*q).s3.c2.im = (*r).s3.c2.im-(*q).s3.c2.im;
    
    ds+= 
      (*q).s3.c0.re*(*q).s3.c0.re+(*q).s3.c0.im*(*q).s3.c0.im+
      (*q).s3.c1.re*(*q).s3.c1.re+(*q).s3.c1.im*(*q).s3.c1.im+
      (*q).s3.c2.re*(*q).s3.c2.re+(*q).s3.c2.im*(*q).s3.c2.im;
    
    tr = ds+kc;
    ts = tr+ks;
    tt = ts-ks;
    ks = ts;
    kc = tr-tt;
  }
  kc = ks+kc;
#ifdef MPI
  MPI_Allreduce(&kc, &ks, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return ks;
#else
  return kc;
#endif

}

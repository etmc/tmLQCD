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
 * File square_and_prod_r.c
 *
 *   void square_and_prod_r(double * const x1, double * const x2, spinor * const S, spinor * const R)
 *     Returns the real part of (*R,*S) and the square norm of *S
 *     It's faster than using "scalar_prod_r" and "square_norm"
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
#include "su3.h"
#include "square_and_prod_r.h"

void square_and_prod_r(double * const x1, double * const x2, spinor * const S, spinor * const R, const int N)
{
  int ix;
  double ALIGN ks,kc,ds,tr,ts,tt;
  double ALIGN xks,xkc,xds,xtr,xts,xtt;
  spinor *s,*r;
  
  ks=0.0;
  kc=0.0;
  
  xks=0.0;
  xkc=0.0;

#if (defined BGL && defined XLC)
  __alignx(16, S);
  __alignx(16, R);
#endif
  
  for (ix = 0; ix < N; ix++)
  {
    s=(spinor *) S + ix;
    r=(spinor *) R + ix;

    ds= r->s0.c0 * conj(s->s0.c0) + r->s0.c1 * conj(s->s0.c1) + r->s0.c2 * conj(s->s0.c2) +  
        r->s1.c0 * conj(s->s1.c0) + r->s1.c1 * conj(s->s1.c1) + r->s1.c2 * conj(s->s1.c2) +  
        r->s2.c0 * conj(s->s2.c0) + r->s2.c1 * conj(s->s2.c1) + r->s2.c2 * conj(s->s2.c2) + 
        r->s3.c0 * conj(s->s3.c0) + r->s3.c1 * conj(s->s3.c1) + r->s3.c2 * conj(s->s3.c2);

    xds=s->s0.c0 * conj(s->s0.c0) + s->s0.c1 * conj(s->s0.c1) + s->s0.c2 * conj(s->s0.c2) +  
        s->s1.c0 * conj(s->s1.c0) + s->s1.c1 * conj(s->s1.c1) + s->s1.c2 * conj(s->s1.c2) +  
        s->s2.c0 * conj(s->s2.c0) + s->s2.c1 * conj(s->s2.c1) + s->s2.c2 * conj(s->s2.c2) +
        s->s3.c0 * conj(s->s3.c0) + s->s3.c1 * conj(s->s3.c1) + s->s3.c2 * conj(s->s3.c2);
    
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


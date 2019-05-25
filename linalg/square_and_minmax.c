/***********************************************************************
 * copyright (C) 2002,2003,2004,2005,2006,2007,2008 Carsten Urbach
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
 *
 * File square_and_max.c
 *
 *   void square_and_max(spinor * const P )
 *     Returns the square norm and max local deviation of *P
 *
 *******************************************************************************/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#ifdef TM_USE_MPI
# include <mpi.h>
#endif
#ifdef TM_USE_OMP
# include <omp.h>
# include "global.h"
#endif
#include <complex.h>
#include "su3.h"
#include "su3adj.h"
#include "su3spinor.h"
#include "square_and_minmax.h"

void square_and_minmax(double * const sum, double * const min, double * const max, const spinor * const P, const int N)
{
  int ix;
  double ALIGN ks,kc,ds,tr,ts,tt;
  spinor *s;
  
  ks=0.0;
  kc=0.0;
  *max = 0.0;
  *min = -1;

#if (defined BGL && defined XLC)
  __alignx(16, S);
  __alignx(16, R);
#endif
  
  for (ix = 0; ix < N; ix++)
  {
    s=(spinor *) P + ix;

    ds=s->s0.c0 * conj(s->s0.c0) + s->s0.c1 * conj(s->s0.c1) + s->s0.c2 * conj(s->s0.c2) +  
      s->s1.c0 * conj(s->s1.c0) + s->s1.c1 * conj(s->s1.c1) + s->s1.c2 * conj(s->s1.c2) +  
      s->s2.c0 * conj(s->s2.c0) + s->s2.c1 * conj(s->s2.c1) + s->s2.c2 * conj(s->s2.c2) +
      s->s3.c0 * conj(s->s3.c0) + s->s3.c1 * conj(s->s3.c1) + s->s3.c2 * conj(s->s3.c2);
    
    tr=ds + kc;
    ts=tr + ks;
    tt=ts-ks;
    ks=ts;
    kc=tr-tt;

    if(ds > *max) *max = ds;
    if(ds < *min || *min < 0) *min = ds;
  }
  kc=ks + kc;
  *sum=kc;

#if defined TM_USE_MPI

  MPI_Allreduce(&kc, sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  MPI_Allreduce(min, &kc, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  *min = kc;

  MPI_Allreduce(max, &kc, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  *max = kc;

#endif

  return;
}

void square_and_minmax_rel(double * const sum, double * const min, double * const max, const spinor * const P, const spinor * const Q, const int N)
{
  int ix;
  double ALIGN ks,kc,ds,dr,tr,ts,tt;
  spinor *s, *r;
  
  ks=0.0;
  kc=0.0;
  *max = 0.0;
  *min = -1;

#if (defined BGL && defined XLC)
  __alignx(16, S);
  __alignx(16, R);
#endif
  
  for (ix = 0; ix < N; ix++)
  {
    s=(spinor *) P + ix;
    r=(spinor *) Q + ix;

    ds=s->s0.c0 * conj(s->s0.c0) + s->s0.c1 * conj(s->s0.c1) + s->s0.c2 * conj(s->s0.c2) +  
      s->s1.c0 * conj(s->s1.c0) + s->s1.c1 * conj(s->s1.c1) + s->s1.c2 * conj(s->s1.c2) +  
      s->s2.c0 * conj(s->s2.c0) + s->s2.c1 * conj(s->s2.c1) + s->s2.c2 * conj(s->s2.c2) +
      s->s3.c0 * conj(s->s3.c0) + s->s3.c1 * conj(s->s3.c1) + s->s3.c2 * conj(s->s3.c2);

    dr=r->s0.c0 * conj(r->s0.c0) + r->s0.c1 * conj(r->s0.c1) + r->s0.c2 * conj(r->s0.c2) +  
      r->s1.c0 * conj(r->s1.c0) + r->s1.c1 * conj(r->s1.c1) + r->s1.c2 * conj(r->s1.c2) +  
      r->s2.c0 * conj(r->s2.c0) + r->s2.c1 * conj(r->s2.c1) + r->s2.c2 * conj(r->s2.c2) +
      r->s3.c0 * conj(r->s3.c0) + r->s3.c1 * conj(r->s3.c1) + r->s3.c2 * conj(r->s3.c2);
    
    ds = ds/dr;

    tr=ds + kc;
    ts=tr + ks;
    tt=ts-ks;
    ks=ts;
    kc=tr-tt;

    if(ds > *max) *max = ds;
    if(ds < *min || *min < 0) *min = ds;
  }
  kc=ks + kc;
  *sum=kc;

#if defined TM_USE_MPI

  MPI_Allreduce(&kc, sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  MPI_Allreduce(min, &kc, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  *min = kc;

  MPI_Allreduce(max, &kc, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  *max = kc;

#endif

  return;
}

void square_and_minmax_abs(double * const sum, double * const min, double * const max,  double * const min_abs, double * const max_abs, const spinor * const P, const int N)
{
  int ix;
  double ALIGN ks,kc,ds,dds,tr,ts,tt;
  spinor *s;
  
  ks=0.0;
  kc=0.0;
  *max = 0.0;
  *min = -1;
  *max_abs = 0.0;
  *min_abs = -1;

#if (defined BGL && defined XLC)
  __alignx(16, S);
  __alignx(16, R);
#endif
  
  for (ix = 0; ix < N; ix++)
  {
    s=(spinor *) P + ix;

    dds=s->s0.c0 * conj(s->s0.c0);
    ds=dds;
    if(dds > *max_abs) *max_abs = dds;
    else if(dds < *min_abs || *min_abs < 0) *min_abs = dds;

    dds=s->s0.c1 * conj(s->s0.c1);
    ds+=dds;
    if(dds > *max_abs) *max_abs = dds;
    else if(dds < *min_abs || *min_abs < 0) *min_abs = dds;

    dds=s->s0.c2 * conj(s->s0.c2);
    ds+=dds;
    if(dds > *max_abs) *max_abs = dds;
    else if(dds < *min_abs || *min_abs < 0) *min_abs = dds;

    dds=s->s1.c0 * conj(s->s1.c0);
    ds+=dds;
    if(dds > *max_abs) *max_abs = dds;
    else if(dds < *min_abs || *min_abs < 0) *min_abs = dds;

    dds=s->s1.c1 * conj(s->s1.c1);
    ds+=dds;
    if(dds > *max_abs) *max_abs = dds;
    else if(dds < *min_abs || *min_abs < 0) *min_abs = dds;

    dds=s->s1.c2 * conj(s->s1.c2);
    ds+=dds;
    if(dds > *max_abs) *max_abs = dds;
    else if(dds < *min_abs || *min_abs < 0) *min_abs = dds;

    dds=s->s2.c0 * conj(s->s2.c0);
    ds+=dds;
    if(dds > *max_abs) *max_abs = dds;
    else if(dds < *min_abs || *min_abs < 0) *min_abs = dds;

    dds=s->s2.c1 * conj(s->s2.c1);
    ds+=dds;
    if(dds > *max_abs) *max_abs = dds;
    else if(dds < *min_abs || *min_abs < 0) *min_abs = dds;

    dds=s->s2.c2 * conj(s->s2.c2);
    ds+=dds;
    if(dds > *max_abs) *max_abs = dds;
    else if(dds < *min_abs || *min_abs < 0) *min_abs = dds;

    dds=s->s3.c0 * conj(s->s3.c0);
    ds+=dds;
    if(dds > *max_abs) *max_abs = dds;
    else if(dds < *min_abs || *min_abs < 0) *min_abs = dds;

    dds=s->s3.c1 * conj(s->s3.c1);
    ds+=dds;
    if(dds > *max_abs) *max_abs = dds;
    else if(dds < *min_abs || *min_abs < 0) *min_abs = dds;

    dds=s->s3.c2 * conj(s->s3.c2);
    ds+=dds;
    if(dds > *max_abs) *max_abs = dds;
    else if(dds < *min_abs || *min_abs < 0) *min_abs = dds;

    tr=ds + kc;
    ts=tr + ks;
    tt=ts-ks;
    ks=ts;
    kc=tr-tt;

    if(ds > *max) *max = ds;
    if(ds < *min || *min < 0) *min = ds;
  }
  kc=ks + kc;
  *sum=kc;

#if defined TM_USE_MPI

  MPI_Allreduce(&kc, sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  MPI_Allreduce(min, &kc, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  *min = kc;

  MPI_Allreduce(max, &kc, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  *max = kc;

  MPI_Allreduce(min_abs, &kc, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  *min_abs = kc;

  MPI_Allreduce(max_abs, &kc, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  *max_abs = kc;

#endif

  return;
}

void square_and_minmax_rel_abs(double * const sum, double * const min, double * const max, double * const min_abs, double * const max_abs, const spinor * const P, const spinor * const Q, const int N)
{
  int ix;
  double ALIGN ks,kc,ds,dds,dr,ddr,tr,ts,tt;
  spinor *s, *r;
  
  ks=0.0;
  kc=0.0;
  *max = 0.0;
  *min = -1;
  *max_abs = 0.0;
  *min_abs = -1;

#if (defined BGL && defined XLC)
  __alignx(16, S);
  __alignx(16, R);
#endif
  
  for (ix = 0; ix < N; ix++)
  {
    s=(spinor *) P + ix;
    r=(spinor *) Q + ix;

    dds=s->s0.c0 * conj(s->s0.c0);
    ddr=r->s0.c0 * conj(r->s0.c0);
    ds=dds;
    dr=ddr;
    dds/=ddr;
    if(dds > *max_abs) *max_abs = dds;
    else if(dds < *min_abs || *min_abs < 0) *min_abs = dds;

    dds=s->s0.c1 * conj(s->s0.c1);
    ddr=r->s0.c1 * conj(r->s0.c1);
    ds+=dds;
    dr+=ddr;
    dds/=ddr;
    if(dds > *max_abs) *max_abs = dds;
    else if(dds < *min_abs || *min_abs < 0) *min_abs = dds;

    dds=s->s0.c2 * conj(s->s0.c2);
    ddr=r->s0.c2 * conj(r->s0.c2);
    ds+=dds;
    dr+=ddr;
    dds/=ddr;
    if(dds > *max_abs) *max_abs = dds;
    else if(dds < *min_abs || *min_abs < 0) *min_abs = dds;

    dds=s->s1.c0 * conj(s->s1.c0);
    ddr=r->s1.c0 * conj(r->s1.c0);
    ds+=dds;
    dr+=ddr;
    dds/=ddr;
    if(dds > *max_abs) *max_abs = dds;
    else if(dds < *min_abs || *min_abs < 0) *min_abs = dds;

    dds=s->s1.c1 * conj(s->s1.c1);
    ddr=r->s1.c1 * conj(r->s1.c1);
    ds+=dds;
    dr+=ddr;
    dds/=ddr;
    if(dds > *max_abs) *max_abs = dds;
    else if(dds < *min_abs || *min_abs < 0) *min_abs = dds;

    dds=s->s1.c2 * conj(s->s1.c2);
    ddr=r->s1.c2 * conj(r->s1.c2);
    ds+=dds;
    dr+=ddr;
    dds/=ddr;
    if(dds > *max_abs) *max_abs = dds;
    else if(dds < *min_abs || *min_abs < 0) *min_abs = dds;

    dds=s->s2.c0 * conj(s->s2.c0);
    ddr=r->s2.c0 * conj(r->s2.c0);
    ds+=dds;
    dr+=ddr;
    dds/=ddr;
    if(dds > *max_abs) *max_abs = dds;
    else if(dds < *min_abs || *min_abs < 0) *min_abs = dds;

    dds=s->s2.c1 * conj(s->s2.c1);
    ddr=r->s2.c1 * conj(r->s2.c1);
    ds+=dds;
    dr+=ddr;
    dds/=ddr;
    if(dds > *max_abs) *max_abs = dds;
    else if(dds < *min_abs || *min_abs < 0) *min_abs = dds;

    dds=s->s2.c2 * conj(s->s2.c2);
    ddr=r->s2.c2 * conj(r->s2.c2);
    ds+=dds;
    dr+=ddr;
    dds/=ddr;
    if(dds > *max_abs) *max_abs = dds;
    else if(dds < *min_abs || *min_abs < 0) *min_abs = dds;

    dds=s->s3.c0 * conj(s->s3.c0);
    ddr=r->s3.c0 * conj(r->s3.c0);
    ds+=dds;
    dr+=ddr;
    dds/=ddr;
    if(dds > *max_abs) *max_abs = dds;
    else if(dds < *min_abs || *min_abs < 0) *min_abs = dds;

    dds=s->s3.c1 * conj(s->s3.c1);
    ddr=r->s3.c1 * conj(r->s3.c1);
    ds+=dds;
    dr+=ddr;
    dds/=ddr;
    if(dds > *max_abs) *max_abs = dds;
    else if(dds < *min_abs || *min_abs < 0) *min_abs = dds;

    dds=s->s3.c2 * conj(s->s3.c2);
    ddr=r->s3.c2 * conj(r->s3.c2);
    ds+=dds;
    dr+=ddr;
    dds/=ddr;
    if(dds > *max_abs) *max_abs = dds;
    else if(dds < *min_abs || *min_abs < 0) *min_abs = dds;
    
    ds = ds/dr;

    tr=ds + kc;
    ts=tr + ks;
    tt=ts-ks;
    ks=ts;
    kc=tr-tt;

    if(ds > *max) *max = ds;
    if(ds < *min || *min < 0) *min = ds;
  }
  kc=ks + kc;
  *sum=kc;

#if defined TM_USE_MPI

  MPI_Allreduce(&kc, sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  MPI_Allreduce(min, &kc, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  *min = kc;

  MPI_Allreduce(max, &kc, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  *max = kc;

  MPI_Allreduce(min_abs, &kc, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  *min_abs = kc;

  MPI_Allreduce(max_abs, &kc, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  *max_abs = kc;

#endif

  return;
}

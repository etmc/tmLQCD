/***********************************************************************
 *
 * Copyright (C) 1995 Ulli Wolff, Stefan Sint
 *               2001,2005 Martin Hasenbusch
 *               2011,2012 Carsten Urbach
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

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#ifdef SSE
# undef SSE
#endif
#ifdef SSE2
# undef SSE2
#endif
#ifdef SSE3
# undef SSE3
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <time.h>
#ifdef MPI
# include <mpi.h>
#endif
#ifdef OMP
# include <omp.h>
#endif
#include "global.h"
#include "su3.h"
#include "sse.h"
#include "su3adj.h"
#include "operator/clovertm_operators.h"
#include "operator/clover_leaf.h"
#include "operator/clover_inline.h"

#define nm1 5
void six_det(_Complex double* const rval, _Complex double a[6][6])
{
  /* required for thread safety */
  _Complex double ALIGN sigma,z;
  _Complex double ALIGN det;
  double ALIGN p[nm1+1];
  double ALIGN s,q;
  int i,j,k;
  int ifail;
  ifail=0;
  /* compute the determinant:*/
  det = 1.0;
  
  for(k = 0; k < nm1; k++) {
    s=0.0;
    for(j = k+1; j <= nm1; ++j) {
      s += conj(a[j][k]) * a[j][k];
    }
    s = sqrt(1. + s / (conj(a[k][k]) * a[k][k]));
    sigma = s * a[k][k];
    
    /* determinant */
    det *= sigma;
    q   = sigma * conj(sigma);
    if (q < tiny_t)
      ifail++;
    
    a[k][k] += sigma;
    p[k]     = sigma * conj(a[k][k]);
    
    /* reflect all columns to the right */
    for(j = k+1; j <= nm1; j++) {
      z = 0.;
      for(i = k; i <= nm1; i++) {
	z += conj(a[i][k]) * a[i][j];
      }
      z /= p[k];
      for(i = k; i <= nm1; i++) {
	a[i][j] -= z * a[i][k];
      }
    }
  }
  sigma = a[nm1][nm1];
  
  /* determinant */
  det *= sigma;
  q = conj(sigma) * sigma;
  
  if(q < tiny_t) {
    ifail++;
  }
  if(g_proc_id == 0 && ifail > 0) {
    fprintf(stderr, "Warning: ifail = %d > 0 in six_det\n", ifail);
  }
  *rval = det;
}


double sw_trace(const int ieo, const double mu) {
  double ALIGN res = 0.0;
#ifdef MPI
  double ALIGN mres;
#endif

#ifdef OMP
#pragma omp parallel
  {
  int thread_num = omp_get_thread_num();
#endif

  int i,x,ioff;
  _Complex double ALIGN a[6][6];
  double ALIGN tra;
  double ALIGN ks,kc,tr,ts,tt;
  _Complex double ALIGN det;

  ks = 0.0;
  kc = 0.0;

  if(ieo==0) {
    ioff=0;
  } 
  else {
    ioff=(VOLUME+RAND)/2;
  }
  
#ifdef OMP
#pragma omp for
#endif
  for(int icx = ioff; icx < (VOLUME/2+ioff); icx++) {
    x = g_eo2lexic[icx];
    for(i=0;i<2;i++) {
      memcpy(a[0], sw[x][i], 36*sizeof(_Complex double));
      // we add the twisted mass term
      if(i == 0) add_tm(a, mu);
      else add_tm(a, -mu);
      // and compute the tr log (or log det)
      six_det(&det,a);
      tra = log(conj(det)*det);
      // we need to compute only the one with +mu
      // the one with -mu must be the complex conjugate!
      
      tr=tra+kc;
      ts=tr+ks;
      tt=ts-ks;
      ks=ts;
      kc=tr-tt;
    }
  }
  kc=ks+kc;

#ifdef OMP
  g_omp_acc_re[thread_num] = kc;
  } /* OpenMP parallel closing brace */

  for(int i = 0; i < omp_num_threads; ++i) {
    res += g_omp_acc_re[i];
  }
#else
  res=kc;
#endif

#ifdef MPI
  MPI_Allreduce(&res, &mres, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return(mres);
#else
  return(res);
#endif

}


// This function computes the trace-log part of the clover term
// in case of even/odd preconditioning in the nd case
//
// it is expected that sw_term is called beforehand such that
// the array sw is populated properly
//
// it is tested to deliver bit-identical results to sw_trace
// if eps is set to zero

double sw_trace_nd(const int ieo, const double mu, const double eps) {
  double ALIGN res = 0.0;
#ifdef MPI
  double ALIGN mres;
#endif

#ifdef OMP
#pragma omp parallel
  {
  int thread_num = omp_get_thread_num();
#endif

  int x,ioff;
  _Complex double ALIGN a[6][6];
  double ALIGN tra;
  double ALIGN ks,kc,tr,ts,tt;
  _Complex double ALIGN det[2];
  double se = (eps*eps)*(eps*eps)*(eps*eps);
  ks=0.0;
  kc=0.0;

  if(ieo==0) {
    ioff=0;
  } 
  else {
    ioff=(VOLUME+RAND)/2;
  }

#ifdef OMP
#pragma omp for
#endif
  for(unsigned int icx = ioff; icx < (VOLUME/2+ioff); icx++) {
    x = g_eo2lexic[icx];
    for(unsigned int i = 0; i < 2; i++) {
      memcpy(a[0], sw[x][i], 36*sizeof(_Complex double));
      // we add the twisted mass term prop to tau^3
      if(i == 0) add_tm(a, mu);
      else add_tm(a, -mu);
      six_det(&det[i], a);
    }
    // and compute the tr log (or log det)
    // for the 2x2 matrix in flavour space
    // with eps*tau^1 in the off diagonal
    tra = log(conj(det[0])*det[0]*conj(det[1])*det[1] - se*se);

    tr=tra+kc;
    ts=tr+ks;
    tt=ts-ks;
    ks=ts;
    kc=tr-tt;
  }
  kc=ks+kc;
  
#ifdef OMP
  g_omp_acc_re[thread_num] = kc;
  } /* OpenMP parallel closing brace */

  for(int i = 0; i < omp_num_threads; ++i) {
    res += g_omp_acc_re[i];
  }
#else
  res=kc;
#endif

#ifdef MPI
  MPI_Allreduce(&res, &mres, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return(mres);
#else
  return(res);
#endif
}

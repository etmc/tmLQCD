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

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#ifdef MPI
# include<mpi.h>
#endif
#include <stdlib.h>
#include <complex.h>
#ifdef OMP
# include <omp.h>
# include <global.h>
#endif
#include "su3.h"
#include "assign_mul_add_r_and_square.h"


#if (defined BGQ && defined XLC)

double assign_mul_add_r_and_square(spinor * const R, const double c, spinor * const S, 
				   const int N, const int parallel) {
  double ALIGN res = 0.0;
#ifdef MPI
  double ALIGN mres;
#endif

#ifdef OMP
#pragma omp parallel
  {
  int thread_num = omp_get_thread_num();
#endif
  vector4double x0, x1, x2, x3, x4, x5, y0, y1, y2, y3, y4, y5;
  vector4double z0, z1, z2, z3, z4, z5, k;
  vector4double r0, r1, r2, r3, r4, r5;
  double *s, *r;
  double ALIGN _c = c;
  double ALIGN ds = 0.0;
#ifndef OMP
  __prefetch_by_load(S);
  __prefetch_by_load(R);
#endif

  k = vec_splats(_c);
  __alignx(32, s);
  __alignx(32, r);
  __alignx(32, S);
  __alignx(32, R);
  r0 = vec_splats(0.);
  r1 = vec_splats(0.);
  r2 = vec_splats(0.);
  r3 = vec_splats(0.);
  r4 = vec_splats(0.);
  r5 = vec_splats(0.);


#ifdef OMP
#pragma omp for 
#endif
  for(int i = 0; i < N; i++) {
    s=(double*)((spinor *) S + i);
    r=(double*)((spinor *) R + i);
    __prefetch_by_load(S + i + 1);
    __prefetch_by_stream(1, R + i + 1);
    x0 = vec_ld(0, r);
    x1 = vec_ld(0, r+4);
    x2 = vec_ld(0, r+8);
    x3 = vec_ld(0, r+12);
    x4 = vec_ld(0, r+16);
    x5 = vec_ld(0, r+20);
    y0 = vec_ld(0, s);
    y1 = vec_ld(0, s+4);
    y2 = vec_ld(0, s+8);
    y3 = vec_ld(0, s+12);
    y4 = vec_ld(0, s+16);
    y5 = vec_ld(0, s+20);
    z0 = vec_madd(k, x0, y0);
    z1 = vec_madd(k, x1, y1);
    z2 = vec_madd(k, x2, y2);
    z3 = vec_madd(k, x3, y3);
    z4 = vec_madd(k, x4, y4);
    z5 = vec_madd(k, x5, y5);
    vec_st(z0, 0, r);
    vec_st(z1, 0, r+4);
    vec_st(z2, 0, r+8);
    vec_st(z3, 0, r+12);
    vec_st(z4, 0, r+16);
    vec_st(z5, 0, r+20);
    r0 = vec_madd(z0, z0, r0);
    r1 = vec_madd(z1, z1, r1);
    r2 = vec_madd(z2, z2, r2);
    r3 = vec_madd(z3, z3, r3);
    r4 = vec_madd(z4, z4, r4);
    r5 = vec_madd(z5, z5, r5);
  }
  x0 = vec_add(r0, r1);
  x1 = vec_add(r2, r3);
  x2 = vec_add(r4, r5);
  y0 = vec_add(x0, x1);
  y1 = vec_add(x2, y0);
  ds = y1[0] + y1[1] + y1[2] + y1[3];

#ifdef OMP
  g_omp_acc_re[thread_num] = ds;
  } /* OpenMP closing brace */

  for(int i = 0; i < omp_num_threads; ++i) {
    res += g_omp_acc_re[i];
  }
#else
  res = ds;
#endif

#  ifdef MPI
  if(parallel) {
    MPI_Allreduce(&res, &mres, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return(mres);
  }
#endif
  return(res);
}

#else

/* R inoutput , c,S input*/
/*   (*R) = c*(*R) + (*S)        c is a real constant   */

double assign_mul_add_r_and_square(spinor * const R, const double c, const spinor * const S, 
				   const int N, const int parallel) {
  double ALIGN res = 0.0;
#ifdef MPI
  double ALIGN mres;
#endif

#ifdef OMP
#pragma omp parallel
  {
  int thread_num = omp_get_thread_num();
#endif
  spinor *r;
  const spinor *s;
  double ALIGN ds = 0.0;

  /* Change due to even-odd preconditioning : VOLUME   to VOLUME/2 */   
#ifdef OMP
#pragma omp for 
#endif
  for (int ix = 0; ix < N; ++ix) {
    r = R + ix;
    s = S + ix;
    
    r->s0.c0 = c * r->s0.c0 + s->s0.c0;
    ds += creal(r->s0.c0)*creal(r->s0.c0) + cimag(r->s0.c0)*cimag(r->s0.c0);
    r->s0.c1 = c * r->s0.c1 + s->s0.c1;
    ds += creal(r->s0.c1)*creal(r->s0.c1) + cimag(r->s0.c1)*cimag(r->s0.c1);
    r->s0.c2 = c * r->s0.c2 + s->s0.c2;    
    ds += creal(r->s0.c2)*creal(r->s0.c2) + cimag(r->s0.c2)*cimag(r->s0.c2);

    r->s1.c0 = c * r->s1.c0 + s->s1.c0;
    ds += creal(r->s1.c0)*creal(r->s1.c0) + cimag(r->s1.c0)*cimag(r->s1.c0);
    r->s1.c1 = c * r->s1.c1 + s->s1.c1;
    ds += creal(r->s1.c1)*creal(r->s1.c1) + cimag(r->s1.c1)*cimag(r->s1.c1);
    r->s1.c2 = c * r->s1.c2 + s->s1.c2;    
    ds += creal(r->s1.c2)*creal(r->s1.c2) + cimag(r->s1.c2)*cimag(r->s1.c2);

    r->s2.c0 = c * r->s2.c0 + s->s2.c0;
    ds += creal(r->s2.c0)*creal(r->s2.c0) + cimag(r->s2.c0)*cimag(r->s2.c0);
    r->s2.c1 = c * r->s2.c1 + s->s2.c1;
    ds += creal(r->s2.c1)*creal(r->s2.c1) + cimag(r->s2.c1)*cimag(r->s2.c1);
    r->s2.c2 = c * r->s2.c2 + s->s2.c2;    
    ds += creal(r->s2.c2)*creal(r->s2.c2) + cimag(r->s2.c2)*cimag(r->s2.c2);

    r->s3.c0 = c * r->s3.c0 + s->s3.c0;
    ds += creal(r->s3.c0)*creal(r->s3.c0) + cimag(r->s3.c0)*cimag(r->s3.c0);
    r->s3.c1 = c * r->s3.c1 + s->s3.c1;
    ds += creal(r->s3.c1)*creal(r->s3.c1) + cimag(r->s3.c1)*cimag(r->s3.c1);
    r->s3.c2 = c * r->s3.c2 + s->s3.c2;   
    ds += creal(r->s3.c2)*creal(r->s3.c2) + cimag(r->s3.c2)*cimag(r->s3.c2);
  }

#ifdef OMP
  g_omp_acc_re[thread_num] = ds;
  } /* OpenMP closing brace */

  for(int i = 0; i < omp_num_threads; ++i) {
    res += g_omp_acc_re[i];
  }
#else
  res = ds;
#endif

#  ifdef MPI
  if(parallel) {
    MPI_Allreduce(&res, &mres, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return(mres);
  }
#endif
  return(res);
}

#endif


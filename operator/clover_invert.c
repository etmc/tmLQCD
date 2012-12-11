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

/*
  !--------------------------------------------------------------!
  !  The subroutine sw_invert is needed for the                  !
  !  even_odd preconditioned Dirac operator with SW improvement. !
  !  Details can be found in  the notes sw.ps on tsun.desy.de    !
  !  by P. Weisz and U. Wolff.                                   !
  !--------------------------------------------------------------!
  !  inversion in place of complex matrix a without pivoting     !
  !  triangularization by householder reflections                 !
  !  inversion of triangular matrix                              !
  !  inverse reflections                                          !
  !--------------------------------------------------------------!
  !  a square matrix, dimensioned 0:n-1                          !
  !  itrouble is counted up, when a dangerously small diagonal   !
  !  element is encountered in the tringular matrix              !
  !  has to be initialized outside                               !
  !                                                              !
  !  Author: U. Wolff, adapted to fortran90 by S. Sint, 29/10/95 !
  !--------------------------------------------------------------!
  !  ported to C by M.Hasenbusch Wed Oct 24 15:46:46 MEST 2001   !
  !______________________________________________________________!
*/


/* six_invert and six_det are called from multiple threads, they are thus
 * made thread-safe by removing the static keywords but they are NOT
 * parallelised for OpenMP */

#define nm1 5
void six_invert(int* ifail ,_Complex double a[6][6])
{
  /* required for thread safety */
  _Complex double ALIGN d[nm1+1],u[nm1+1];
  _Complex double ALIGN sigma,z;
  double ALIGN p[nm1+1];
  double ALIGN s,q;
  int i,j,k;
  *ifail=0;
  for(k = 0; k < nm1; ++k)
  {
    s=0.0;
    for(j = k+1; j <= nm1; ++j)
      s += conj(a[j][k]) * a[j][k];
    s = sqrt(1. + s / (conj(a[k][k]) * a[k][k]));
    sigma = s * a[k][k];

    a[k][k] += sigma;
    p[k] = conj(sigma) * a[k][k];
    q = conj(sigma) * sigma;
    if (q < tiny_t)
      (*ifail)++;
    d[k] = -conj(sigma) / q;

    /* reflect all columns to the right */
    for(j = k+1; j <= nm1; ++j)
    {
      z = 0.0;
      for(i = k; i <= nm1; ++i)
	z += conj(a[i][k]) * a[i][j];
      z /= p[k];
      for(i = k; i <= nm1; ++i)
	a[i][j] -= z * a[i][k];
    }
  }
  sigma = a[nm1][nm1];
  q = conj(sigma) * sigma;
  if (q < tiny_t)
    (*ifail)++;
  d[nm1] = conj(sigma) / q;

  /*  inversion of upper triangular matrix in place
      (diagonal elements done already): */

  for(k = nm1; k >= 0; k--) {
    for(i = k-1; i >= 0;i--) {
      z = 0.0;
      for(j = i+1; j < k; j++)
	z += a[i][j] * a[j][k];
      z += a[i][k] * d[k];
      a[i][k] = -z * d[i];
    }
  }     
  /* execute reflections in reverse order from the right: */
  
  a[nm1][nm1] = d[nm1];
  for(k = nm1-1; k >= 0; k--)
  {
    for(j=k;j<=nm1;j++)
      u[j] = a[j][k];
    a[k][k] = d[k];
    for(j = k+1; j <= nm1; j++)
      a[j][k] = 0.0;
    for(i = 0; i <= nm1; i++)
    {
      z = 0.0;
      for(j = k; j <= nm1; j++)
        z += a[i][j] * u[j];
      z /= p[k];         /* normalization */
      
      for(j = k; j <= nm1; j++)
        a[i][j] -= conj(u[j]) * z; /* reflection */
    }
  }
}

// This function computes the inverse of
// (1 + T_ee \pm I\mu\gamma_5)
//
// + is stored in sw_inv[0-(VOLUME/2-1)] 
// - is stored in sw_inv[VOLUME/2-(VOLUME-1)]

void sw_invert(const int ieo, const double mu) {
#ifdef OMP
#pragma omp parallel
  {
#endif
  int icy;
  int ioff, err=0;
  int i, x;
  su3 ALIGN v;
  _Complex double ALIGN a[6][6];

  if(ieo==0) {
    ioff=0;
  } 
  else {
    ioff=(VOLUME+RAND)/2;
  }

#ifndef OMP
  icy=0;
#endif

#ifdef OMP
#pragma omp for
#endif
  for(int icx = ioff; icx < (VOLUME/2+ioff); icx++) {
#ifdef OMP
    icy = icx - ioff;
#endif
    x = g_eo2lexic[icx];

    for(i = 0; i < 2; i++) {
      populate_6x6_matrix(a, &sw[x][0][i], 0, 0);
      populate_6x6_matrix(a, &sw[x][1][i], 0, 3);
      _su3_dagger(v, sw[x][1][i]); 
      populate_6x6_matrix(a, &v, 3, 0);
      populate_6x6_matrix(a, &sw[x][2][i], 3, 3);
      // we add the twisted mass term
      if(i == 0) add_tm(a, +mu);
      else add_tm(a, -mu);
      // and invert the resulting matrix

      six_invert(&err,a); 
      // here we need to catch the error! 
      if(err > 0 && g_proc_id == 0) {
	printf("# inversion failed in six_invert code %d\n", err);
	err = 0;
      }

      /*  copy "a" back to sw_inv */
      get_3x3_block_matrix(&sw_inv[icy][0][i], a, 0, 0);
      get_3x3_block_matrix(&sw_inv[icy][1][i], a, 0, 3);
      get_3x3_block_matrix(&sw_inv[icy][2][i], a, 3, 3);
      get_3x3_block_matrix(&sw_inv[icy][3][i], a, 3, 0);
    }

    if(fabs(mu) > 0.) {
      for(i = 0; i < 2; i++) {
	populate_6x6_matrix(a, &sw[x][0][i], 0, 0);
	populate_6x6_matrix(a, &sw[x][1][i], 0, 3);
	_su3_dagger(v, sw[x][1][i]); 
	populate_6x6_matrix(a, &v, 3, 0);
	populate_6x6_matrix(a, &sw[x][2][i], 3, 3);

	// we add the twisted mass term
	if(i == 0) add_tm(a, -mu);
	else add_tm(a, +mu);
	// and invert the resulting matrix
	six_invert(&err,a); 
	// here we need to catch the error! 
	if(err > 0 && g_proc_id == 0) {
	  printf("# %d\n", err);
	  err = 0;
	}

	/*  copy "a" back to sw_inv */
	get_3x3_block_matrix(&sw_inv[icy+VOLUME/2][0][i], a, 0, 0);
	get_3x3_block_matrix(&sw_inv[icy+VOLUME/2][1][i], a, 0, 3);
	get_3x3_block_matrix(&sw_inv[icy+VOLUME/2][2][i], a, 3, 3);
	get_3x3_block_matrix(&sw_inv[icy+VOLUME/2][3][i], a, 3, 0);
      }
    }
#ifndef OMP
    ++icy;
#endif
  }
#ifdef OMP
  } /* OpenMP closing brace */
#endif
  return;
}

// This function computes
//
// 1/((1+T)^2 + barmu^2 - bareps^2)^{-1}
//
// for all even x,
// which is stored in sw_inv[0-(VOLUME/2-1)]
//
// it is the complement of sw_invert for the
// non-degenerate case
// multiplication with
// (1+T - i\bar\mu\gamma_5\tau^3 + \bar\epsion\tau^1)
// must be done elsewhere because of flavour structure

void sw_invert_nd(const double mshift) {
#ifdef OMP
#pragma omp parallel
  {
#endif
  int err=0;
  int i, x;
  su3 ALIGN v;
  _Complex double ALIGN a[6][6], b[6][6];

#ifdef OMP
#pragma omp for
#endif
  for(int icx = 0; icx < (VOLUME/2); icx++) {
    x = g_eo2lexic[icx];

    for(i = 0; i < 2; i++) {
      populate_6x6_matrix(a, &sw[x][0][i], 0, 0);
      populate_6x6_matrix(a, &sw[x][1][i], 0, 3);
      _su3_dagger(v, sw[x][1][i]); 
      populate_6x6_matrix(a, &v, 3, 0);
      populate_6x6_matrix(a, &sw[x][2][i], 3, 3);

      // compute (1+T)^2 and store in b
      mult_6x6(b, a, a);
      // we add the mass shift term, which is a real number
      add_shift_6x6(b, mshift);
      // so b = (1+T)^2 + shift
      // now invert this matrix
      six_invert(&err, b); 
      // here we need to catch the error! 
      if(err > 0 && g_proc_id == 0) {
	printf("# inversion failed in six_invert_nd code %d\n", err);
	err = 0;
      }

      /*  copy "a" back to sw_inv */
      get_3x3_block_matrix(&sw_inv[icx][0][i], b, 0, 0);
      get_3x3_block_matrix(&sw_inv[icx][1][i], b, 0, 3);
      get_3x3_block_matrix(&sw_inv[icx][2][i], b, 3, 3);
      get_3x3_block_matrix(&sw_inv[icx][3][i], b, 3, 0);
    }
  }
#ifdef OMP
  } /* OpenMP closing brace */
#endif
  return;
}

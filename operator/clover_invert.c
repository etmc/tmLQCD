/***********************************************************************
 *
 * Copyright (C) 1995 Ulli Wolff, Stefan Sint
 *               2001,2005 Martin Hasenbusch
 *               2011,2012 Carsten Urbach
 *               2017      Bartosz Kostrzewa
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
# include<tmlqcd_config.h>
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
#ifdef TM_USE_MPI
# include <mpi.h>
#endif
#ifdef TM_USE_OMP
# include <omp.h>
#endif
#include "global.h"
#include "su3.h"
#include "sse.h"
#include "su3adj.h"
#include "operator/clovertm_operators.h"
#include "operator/clover_leaf.h"
#include "operator/clover_inline.h"
#include "gettime.h"

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

// for debugging purposes, the print statements can be enabled
// #define CLOVER_INVERT_DEBUG

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
  double swtime = gettime();
#ifdef TM_USE_OMP
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

#ifndef TM_USE_OMP
  icy=0;
#endif

#ifdef TM_USE_OMP
#pragma omp for
#endif
  for(int icx = ioff; icx < (VOLUME/2+ioff); icx++) {
#ifdef TM_USE_OMP
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
      // here we need to catch the error! probably this check should be
      // performed by all processes?! 
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
#ifndef TM_USE_OMP
    ++icy;
#endif
  }
#ifdef TM_USE_OMP
  } /* OpenMP closing brace */
#endif
  tm_stopwatch(0, 2, "", __func__, swtime);
  return;
}

/* This function computes 
 * 
 *   \bar\epsilon / ( (1 + Tee)^2 + \bar\mu^2 - \bar\epsilon^2 )
 * 
 * for use in the QPhiX packing routine for the non-degenerate
 * clover doublet
 *
 * sw_inv should contain 
 *   1 / ( (1 + Tee)^2 + \bar\mu^2 - \bar\epsilon^2 ) 
 * the last VOLUME/2 elements (which should not be relevant at this stage) 
 * of sw_inv will be overwritten
 */ 

void sw_invert_epsbar(const double epsbar) {
  double swtime = gettime();
#ifdef TM_USE_OMP
#pragma omp parallel
  {
#endif
  int icy, i;
  _Complex double ALIGN a[6][6];

#ifndef TM_USE_OMP
  icy=VOLUME/2;
#endif

#ifdef TM_USE_OMP
#pragma omp for
#endif
  for(int icx = 0; icx < VOLUME/2; icx++) {
#ifdef TM_USE_OMP
    icy = icx + VOLUME/2;
#endif
    for(i = 0; i < 2; i++) {
      // extract 1/((1+Tee)^2 + \bar\mu^2 - \bar\epsilon^2)
      populate_6x6_matrix(a, &sw_inv[icx][0][i], 0, 0);
      populate_6x6_matrix(a, &sw_inv[icx][1][i], 0, 3);
      populate_6x6_matrix(a, &sw_inv[icx][2][i], 3, 3);
      populate_6x6_matrix(a, &sw_inv[icx][3][i], 3, 0);

      // scale by epsbar
      scale_real_6x6(a, epsbar);

#ifdef CLOVER_INVERT_DEBUG
      if(icx==0) print_6x6(a, "sw_invert_epsbar epsilon*sw_inv");
#endif

      /*  and write the result into the last VOLUME/2 elements of sw_inv */
      get_3x3_block_matrix(&sw_inv[icy][0][i], a, 0, 0);
      get_3x3_block_matrix(&sw_inv[icy][1][i], a, 0, 3);
      get_3x3_block_matrix(&sw_inv[icy][2][i], a, 3, 3);
      get_3x3_block_matrix(&sw_inv[icy][3][i], a, 3, 0);
    }
#ifndef TM_USE_OMP
    icy++;
#endif
  }
#ifdef TM_USE_OMP
  } /* OpenMP closing brace */
#endif
  tm_stopwatch(0, 2, "", __func__, swtime);
  return;
}

/* This function computes 
 * 
 *   (1 + Tee - I*\bar\mu\gamma_5\tau3 ) / ( 1 + Tee + \bar\mu^2 - \bar\eps^2 )
 * 
 * for use in the QPhiX packing routine for the non-degenerate
 * clover doublet
 *
 * sw should be populated with (1+Tee) and 
 * sw_inv should contain 1 / ( (1 + Tee)^2 + \bar\mu^2 - \bar\eps^2 )
 *
 * !! all elements of sw_inv will be overwritten !!
 */ 

void sw_invert_mubar(const double mubar) {
  double swtime = gettime();
#ifdef TM_USE_OMP
#pragma omp parallel
  {
#endif
  int err=0;
  int icy, i, x;
  su3 ALIGN v;
  _Complex double ALIGN a[6][6];
  _Complex double ALIGN b[6][6];
  _Complex double ALIGN c[6][6];
  _Complex double ALIGN d[6][6];

#ifndef TM_USE_OMP
  icy=VOLUME/2;
#endif

#ifdef TM_USE_OMP
#pragma omp for
#endif
  for(int icx = 0; icx < VOLUME/2; icx++) {
#ifdef TM_USE_OMP
    icy = icx + VOLUME/2;
#endif
    x = g_eo2lexic[icx];
    
    for(i = 0; i < 2; i++) {
      // extract (1+Tee)
      populate_6x6_matrix(a, &sw[x][0][i], 0, 0);
      populate_6x6_matrix(a, &sw[x][1][i], 0, 3);
      _su3_dagger(v, sw[x][1][i]); 
      populate_6x6_matrix(a, &v, 3, 0);
      populate_6x6_matrix(a, &sw[x][2][i], 3, 3);

      copy_6x6(b,a);

      // we add the twisted quark masses for both the 'up' and the 'down' flavour
      // (note that this is the inverse, so -mu is associated with 'up')
      // the i index denotes the halfspinor block and thus implements gamma5
      add_tm(a, -(i==0?1.0:-1.0)*mubar);
      add_tm(b, +(i==0?1.0:-1.0)*mubar);

#ifdef CLOVER_INVERT_DEBUG
      if(icx==0) {
        print_6x6(a,"sw_invert_mubar sw_up");
        print_6x6(b,"sw_invert_mubar sw_dn");
      }
#endif
  
      // extract 1/((1+Tee)^2 + \bar\mu^2 - \bar\eps^2)
      populate_6x6_matrix(c, &sw_inv[icx][0][i], 0, 0);
      populate_6x6_matrix(c, &sw_inv[icx][1][i], 0, 3);
      populate_6x6_matrix(c, &sw_inv[icx][2][i], 3, 3);
      populate_6x6_matrix(c, &sw_inv[icx][3][i], 3, 0);
  
      // multiply the two together and store in d
      mult_6x6(d, a, c);
  
      /*  and write the result into sw_inv */
      get_3x3_block_matrix(&sw_inv[icx][0][i], d, 0, 0);
      get_3x3_block_matrix(&sw_inv[icx][1][i], d, 0, 3);
      get_3x3_block_matrix(&sw_inv[icx][2][i], d, 3, 3);
      get_3x3_block_matrix(&sw_inv[icx][3][i], d, 3, 0);

#ifdef CLOVER_INVERT_DEBUG
      if(icx==0) print_6x6(d,"sw_invert_mubar sw_inv_up");
#endif

      // and the same for the 'down'
      mult_6x6(d, b, c);
  
      get_3x3_block_matrix(&sw_inv[icy][0][i], d, 0, 0);
      get_3x3_block_matrix(&sw_inv[icy][1][i], d, 0, 3);
      get_3x3_block_matrix(&sw_inv[icy][2][i], d, 3, 3);
      get_3x3_block_matrix(&sw_inv[icy][3][i], d, 3, 0);

#ifdef CLOVER_INVERT_DEBUG
      if(icx==0) print_6x6(d,"sw_invert_mubar sw_inv_dn");
#endif
    }
#ifndef TM_USE_OMP
    icy++;
#endif
  }
#ifdef TM_USE_OMP
  } /* OpenMP closing brace */
#endif
  tm_stopwatch(0, 2, "", __func__, swtime);
  return;
}

// This function computes
//
// 1/((1+T)^2 + \bar\mu^2 - \bar\epsilon^2)
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
  double swtime = gettime();
#ifdef TM_USE_OMP
#pragma omp parallel
  {
#endif
  int err=0;
  int i, x;
  su3 ALIGN v;
  _Complex double ALIGN a[6][6], b[6][6];

#ifdef TM_USE_OMP
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

#ifdef CLOVER_INVERT_DEBUG
      if(icx==0) print_6x6(a, "sw_invert_nd sw");
#endif

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

#ifdef CLOVER_INVERT_DEBUG
      if(icx==0) print_6x6(b, "sw_invert_nd sw_inv");
#endif

      /*  copy "a" back to sw_inv */
      get_3x3_block_matrix(&sw_inv[icx][0][i], b, 0, 0);
      get_3x3_block_matrix(&sw_inv[icx][1][i], b, 0, 3);
      get_3x3_block_matrix(&sw_inv[icx][2][i], b, 3, 3);
      get_3x3_block_matrix(&sw_inv[icx][3][i], b, 3, 0);
    }
  }
#ifdef TM_USE_OMP
  } /* OpenMP closing brace */
#endif
  tm_stopwatch(0, 2, "", __func__, swtime);
  return;
}

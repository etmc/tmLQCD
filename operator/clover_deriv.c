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

// this is (-tr(1+T_ee(+mu)) -tr(1+T_ee(-mu)))      
// (or T_oo of course)
// 
// see equation (24) of hep-lat/9603008             
//
// or in more detail the insertion matrix at even sites
// is computed
// and stored in swm and swp, which are 4 su3 matrices 
// each per site
// refereing to upwards or downwards winding paths  
//
// swm and swp are representing 6x6 complex matrices
// (colour matrices)
//
// this function depends on mu

void sw_deriv(const int ieo, const double mu) {
#ifdef OMP
#pragma omp parallel
  {
#endif
  int icy;
  int ioff;
  int x;
  double fac = 1.0000;
  su3 ALIGN lswp[4], lswm[4];

  /* convention: Tr clover-leaf times insertion */
  if(ieo == 0) {
    ioff=0;
  } 
  else {
    ioff = (VOLUME+RAND)/2;
  }
  if(fabs(mu) > 0.) fac = 0.5;

#ifndef OMP
  icy = 0;
#endif

#ifdef OMP
#pragma omp for
#endif
  for(int icx = ioff; icx < (VOLUME/2+ioff); icx++) {
#ifdef OMP
    icy = icx - ioff;
#endif
    x = g_eo2lexic[icx];
    /* compute the insertion matrix */
    _su3_plus_su3(lswp[0], sw_inv[icy][0][1], sw_inv[icy][0][0]);
    _su3_plus_su3(lswp[1], sw_inv[icy][1][1], sw_inv[icy][1][0]);
    _su3_plus_su3(lswp[2], sw_inv[icy][2][1], sw_inv[icy][2][0]);
    _su3_plus_su3(lswp[3], sw_inv[icy][3][1], sw_inv[icy][3][0]);

    _su3_minus_su3(lswm[0], sw_inv[icy][0][1], sw_inv[icy][0][0]);
    _su3_minus_su3(lswm[1], sw_inv[icy][1][1], sw_inv[icy][1][0]);
    _su3_minus_su3(lswm[2], sw_inv[icy][2][1], sw_inv[icy][2][0]);
    _su3_minus_su3(lswm[3], sw_inv[icy][3][1], sw_inv[icy][3][0]);
    
    /* add up to swm[] and swp[] */
    _su3_refac_acc(swm[x][0], fac, lswm[0]);
    _su3_refac_acc(swm[x][1], fac, lswm[1]);
    _su3_refac_acc(swm[x][2], fac, lswm[2]);
    _su3_refac_acc(swm[x][3], fac, lswm[3]);
    _su3_refac_acc(swp[x][0], fac, lswp[0]);
    _su3_refac_acc(swp[x][1], fac, lswp[1]);
    _su3_refac_acc(swp[x][2], fac, lswp[2]);
    _su3_refac_acc(swp[x][3], fac, lswp[3]);
    if(fabs(mu) > 0.) {
      /* compute the insertion matrix */
      _su3_plus_su3(lswp[0], sw_inv[icy+VOLUME/2][0][1], sw_inv[icy+VOLUME/2][0][0]);
      _su3_plus_su3(lswp[1], sw_inv[icy+VOLUME/2][1][1], sw_inv[icy+VOLUME/2][1][0]);
      _su3_plus_su3(lswp[2], sw_inv[icy+VOLUME/2][2][1], sw_inv[icy+VOLUME/2][2][0]);
      _su3_plus_su3(lswp[3], sw_inv[icy+VOLUME/2][3][1], sw_inv[icy+VOLUME/2][3][0]); 

      _su3_minus_su3(lswm[0], sw_inv[icy+VOLUME/2][0][1], sw_inv[icy+VOLUME/2][0][0]);
      _su3_minus_su3(lswm[1], sw_inv[icy+VOLUME/2][1][1], sw_inv[icy+VOLUME/2][1][0]);
      _su3_minus_su3(lswm[2], sw_inv[icy+VOLUME/2][2][1], sw_inv[icy+VOLUME/2][2][0]);
      _su3_minus_su3(lswm[3], sw_inv[icy+VOLUME/2][3][1], sw_inv[icy+VOLUME/2][3][0]);
      
      /* add up to swm[] and swp[] */
      _su3_refac_acc(swm[x][0], fac, lswm[0]);
      _su3_refac_acc(swm[x][1], fac, lswm[1]);
      _su3_refac_acc(swm[x][2], fac, lswm[2]);
      _su3_refac_acc(swm[x][3], fac, lswm[3]);
      _su3_refac_acc(swp[x][0], fac, lswp[0]);
      _su3_refac_acc(swp[x][1], fac, lswp[1]);
      _su3_refac_acc(swp[x][2], fac, lswp[2]);
      _su3_refac_acc(swp[x][3], fac, lswp[3]);
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

void sw_deriv_nd(const int ieo) {
#ifdef OMP
#pragma omp parallel
  {
#endif
  int icy;
  int ioff;
  int x;
  double fac = 1.0000;
  su3 ALIGN lswp[4], lswm[4], v;
  _Complex double ALIGN a0[6][6], a1[6][6], b[6][6], c[6][6];

  /* convention: Tr clover-leaf times insertion */
  if(ieo == 0) {
    ioff=0;
  } 
  else {
    ioff = (VOLUME+RAND)/2;
  }

#ifndef OMP
  icy = 0;
#endif

#ifdef OMP
#pragma omp for
#endif
  for(int icx = ioff; icx < (VOLUME/2+ioff); icx++) {
#ifdef OMP
    icy = icx - ioff;
#endif
    x = g_eo2lexic[icx];
    /* compute the insertion matrix */
    populate_6x6_matrix(b, &sw[x][0][0], 0, 0);
    populate_6x6_matrix(b, &sw[x][1][0], 0, 3);
    _su3_dagger(v, sw[x][1][0]); 
    populate_6x6_matrix(b, &v, 3, 0);
    populate_6x6_matrix(b, &sw[x][2][0], 3, 3);

    populate_6x6_matrix(c, &sw_inv[icy][0][0], 0, 0);
    populate_6x6_matrix(c, &sw_inv[icy][1][0], 0, 3);
    populate_6x6_matrix(c, &sw_inv[icy][2][0], 3, 3);
    populate_6x6_matrix(c, &sw_inv[icy][3][0], 3, 0);

    mult_6x6(a0, b, c);

    populate_6x6_matrix(b, &sw[x][0][1], 0, 0);
    populate_6x6_matrix(b, &sw[x][1][1], 0, 3);
    _su3_dagger(v, sw[x][1][1]); 
    populate_6x6_matrix(b, &v, 3, 0);
    populate_6x6_matrix(b, &sw[x][2][1], 3, 3);

    populate_6x6_matrix(c, &sw_inv[icy][0][1], 0, 0);
    populate_6x6_matrix(c, &sw_inv[icy][1][1], 0, 3);
    populate_6x6_matrix(c, &sw_inv[icy][2][1], 3, 3);
    populate_6x6_matrix(c, &sw_inv[icy][3][1], 3, 0);

    mult_6x6(a1, b, c);
    add_6x6(b, a1, a0);
    get_3x3_block_matrix(&lswp[0], b, 0, 0);
    get_3x3_block_matrix(&lswp[1], b, 0, 3);
    get_3x3_block_matrix(&lswp[2], b, 3, 3);
    get_3x3_block_matrix(&lswp[3], b, 3, 0);

    sub_6x6(b, a1, a0);
    get_3x3_block_matrix(&lswm[0], b, 0, 0);
    get_3x3_block_matrix(&lswm[1], b, 0, 3);
    get_3x3_block_matrix(&lswm[2], b, 3, 3);
    get_3x3_block_matrix(&lswm[3], b, 3, 0);
    
    /* add up to swm[] and swp[] */
    _su3_refac_acc(swm[x][0], fac, lswm[0]);
    _su3_refac_acc(swm[x][1], fac, lswm[1]);
    _su3_refac_acc(swm[x][2], fac, lswm[2]);
    _su3_refac_acc(swm[x][3], fac, lswm[3]);
    _su3_refac_acc(swp[x][0], fac, lswp[0]);
    _su3_refac_acc(swp[x][1], fac, lswp[1]);
    _su3_refac_acc(swp[x][2], fac, lswp[2]);
    _su3_refac_acc(swp[x][3], fac, lswp[3]);
#ifndef OMP
    ++icy;
#endif
  }
#ifdef OMP
  } /* OpenMP closing brace */
#endif
  return;
}


// direct product of Y_e(o) and X_e(o) in colour space   
// with insertion matrix at site x
// see equation (22) of hep-lat/9603008                  
// result is again stored in swm and swp                 
// includes a gamma5 multiplication for kk

void sw_spinor(const int ieo, const spinor * const kk, const spinor * const ll, 
	       const double fac) {
#ifdef OMP
#pragma omp parallel
  {
#endif

  int ioff;
  int icx;
  int x;
  const spinor *r,*s;
  su3 ALIGN v0,v1,v2,v3;
  su3 ALIGN u0,u1,u2,u3;
  su3 ALIGN lswp[4],lswm[4];

  if(ieo == 0) {
    ioff=0;
  } 
  else {
    ioff=(VOLUME+RAND)/2;
  }
  /************************ loop over half of the lattice sites ***********/

#ifdef OMP
#pragma omp for
#endif  
  for(icx = ioff; icx < (VOLUME/2+ioff); icx++) {
    x = g_eo2lexic[icx];
    r = kk + icx - ioff;
    s = ll + icx - ioff;
    
    _vector_tensor_vector(v0,(*r).s0,(*s).s0);
    _vector_tensor_vector(v1,(*r).s0,(*s).s1);
    _vector_tensor_vector(v2,(*r).s1,(*s).s1);
    _vector_tensor_vector(v3,(*r).s1,(*s).s0);
    // mvector takes g5 into account
    _mvector_tensor_vector(u0,(*r).s2,(*s).s2);
    _mvector_tensor_vector(u1,(*r).s2,(*s).s3);
    _mvector_tensor_vector(u2,(*r).s3,(*s).s3);
    _mvector_tensor_vector(u3,(*r).s3,(*s).s2);
    
    /* compute the insertion matrix */
    _su3_plus_su3(lswp[0],u0,v0);
    _su3_plus_su3(lswp[1],u1,v1);
    _su3_plus_su3(lswp[2],u2,v2);
    _su3_plus_su3(lswp[3],u3,v3);

    _su3_minus_su3(lswm[0],u0,v0);
    _su3_minus_su3(lswm[1],u1,v1);
    _su3_minus_su3(lswm[2],u2,v2);
    _su3_minus_su3(lswm[3],u3,v3);
    
    /* add up to swm[0] and swp[0] */
    _su3_refac_acc(swm[x][0], fac, lswm[0]);
    _su3_refac_acc(swm[x][1], fac, lswm[1]);
    _su3_refac_acc(swm[x][2], fac, lswm[2]);
    _su3_refac_acc(swm[x][3], fac, lswm[3]);
    _su3_refac_acc(swp[x][0], fac, lswp[0]);
    _su3_refac_acc(swp[x][1], fac, lswp[1]);
    _su3_refac_acc(swp[x][2], fac, lswp[2]);
    _su3_refac_acc(swp[x][3], fac, lswp[3]);
  }
#ifdef OMP
  } /* OpenMP closing brace */
#endif
  return;
}



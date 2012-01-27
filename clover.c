/***********************************************************************
 * $Id: operator.c 1763 2011-04-21 11:51:47Z reker $ 
 *
 * Copyright (C) 2005 Martin Hasenbusch
 *               2011 Carsten Urbach
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
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <time.h>
#ifdef MPI
# include <mpi.h>
#endif
#include "global.h"
#include "su3.h"
#include "sse.h"
#include "Hopping_Matrix.h"
#include "clover.h"


su3 *** sw;
su3 *** sw_inv;

void clover_inv(const int ieo, spinor * const l);
void clover_gamma5(const int ieo, 
		   spinor * const l, spinor * const k, spinor * const j,
		   const double q_off);
void clover(const int ieo, 
	    spinor * const l, spinor * const k, spinor * const j,
	    const double q_off);

/*******************************************************************
 *
 *
 * \hat Q_{+} =
 * \gamma_5(M_{oo}^+ - M_{oe}(M_{ee}^+ )^{-1}M_{eo})
 *
 * with clover term!
 * see documentation for details
 * k is the input field
 * l is the output field
 *
 * it acts only on the odd part or only 
 * on a half spinor
 *******************************************************************/


void Qsw_psi(spinor * const l, spinor * const k) {
  Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX+1], k);
  clover_inv(EO, g_spinor_field[DUM_MATRIX+1]);
  Hopping_Matrix(OE, g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+1]);
  /* temporarily set to g_mu -> needs to be fixed */
  clover_gamma5(OE, l, k, g_spinor_field[DUM_MATRIX], g_mu);
}

void Msw_psi(spinor * const l, spinor * const k) {
  Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX+1], k);
  clover_inv(EO, g_spinor_field[DUM_MATRIX+1]);
  Hopping_Matrix(OE, g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+1]);
  /* temporarily set to g_mu -> needs to be fixed */
  clover(OE, l, k, g_spinor_field[DUM_MATRIX], g_mu);
}

/**********************************************************
 *
 * clover_inv applies the inverse of the clover term
 * to spinor field l
 * it is assumed that the corresponding inverted matrices
 * are stored in sw_inv
 *
 * this is needed for even/odd preconditioning
 *
 **********************************************************/

void clover_inv(const int ieo, spinor * const l) {
  int ix;
  int ioff,icx;
  su3 *w1,*w2,*w3;
  spinor *rn;
  static su3_vector psi, chi, phi1, phi3;

  if(ieo == 0) {
    ioff = 0;
  } 
  else {
    ioff = (VOLUME+RAND)/2;
  }
  /************************ loop over all lattice sites *************************/

  for(icx = ioff; icx < (VOLUME/2+ioff); icx++) {
    ix = g_eo2lexic[icx];

    rn = l + icx - ioff;
    /* noch schlimmer murks fuer den clover-term*/
    _vector_assign(phi1,(*rn).s0);
    _vector_assign(phi3,(*rn).s2);
    
    w1=&sw_inv[ix][0][0];
    w2=w1+2;  /* &sw_inv[ix][1][0]; */
    w3=w1+4;  /* &sw_inv[ix][2][0]; */
    _su3_multiply(psi,*w1,phi1); 
    _su3_multiply(chi,*w2,(*rn).s1);
    _vector_add((*rn).s0,psi,chi);
    _su3_inverse_multiply(psi,*w2,phi1); 
    _su3_multiply(chi,*w3,(*rn).s1);
    _vector_add((*rn).s1,psi,chi);
    
    w1++; /* &sw_inv[ix][0][1]; */
    w2++; /* &sw_inv[ix][1][1]; */
    w3++; /* &sw_inv[ix][2][1]; */
    _su3_multiply(psi,*w1,phi3); 
    _su3_multiply(chi,*w2,(*rn).s3);
    _vector_add((*rn).s2,psi,chi);
    _su3_inverse_multiply(psi,*w2,phi3); 
    _su3_multiply(chi,*w3,(*rn).s3);
    _vector_add((*rn).s3,psi,chi);
    /******************************** end of loop *********************************/
  }
  return;
}

/**************************************************************
 *
 * clover_gamma5 applies the clover term to spinor k, adds k 
 * to j then and stores it in l multiplied by gamma_5
 *
 * it is assumed that the clover leaf is computed and stored
 * in sw
 * the corresponding routine can be found in clover_leaf.c
 *
 **************************************************************/

void clover_gamma5(const int ieo, 
		   spinor * const l, spinor * const k, spinor * const j,
		   const double q_off) {
  int ix;
  int ioff,icx;
  su3 *w1,*w2,*w3;
  spinor *r,*s,*t;
  static su3_vector chi, psi1, psi2;

  if(ieo == 0) {
    ioff = 0;
  } 
  else {
    ioff = (VOLUME+RAND)/2;
  }
/************************ loop over all lattice sites *************************/
  for(icx = ioff; icx < (VOLUME/2+ioff); icx++) {
    ix = g_eo2lexic[icx];
    
    r = l + icx-ioff;
    s = k + icx-ioff;
    t = j + icx-ioff;
    
    w1=&sw[ix][0][0];
    w2=w1+2; /*&sw[ix][1][0];*/
    w3=w1+4; /*&sw[ix][2][0];*/
    _su3_multiply(psi1,*w1,(*s).s0); 
    _su3_multiply(chi,*w2,(*s).s1);
    _vector_add_assign(psi1,chi);
    _su3_inverse_multiply(psi2,*w2,(*s).s0); 
    _su3_multiply(chi,*w3,(*s).s1);
    _vector_add_assign(psi2,chi); 
    /*********************** contribution from q_off ***************************/
    _vector_add_mul(psi1,q_off,(*s).s0);
    _vector_add_mul(psi2,q_off,(*s).s1);
    /**************************************************************/
    _vector_sub((*r).s0,psi1,(*t).s0);
    _vector_sub((*r).s1,psi2,(*t).s1);
    
    w1++; /*=&sw[ix][0][1];*/
    w2++; /*=&sw[ix][1][1];*/
    w3++; /*=&sw[ix][2][1];*/
    _su3_multiply(psi1,*w1,(*s).s2); _su3_multiply(chi,*w2,(*s).s3);
    _vector_add_assign(psi1,chi); 
    _su3_inverse_multiply(psi2,*w2,(*s).s2); _su3_multiply(chi,*w3,(*s).s3);
    _vector_add_assign(psi2,chi); 
    /*********************** contribution from q_off ***************************/
    _vector_add_mul(psi1,q_off,(*s).s2);
    _vector_add_mul(psi2,q_off,(*s).s3);
    /**************** multiply with  gamma5 included ******************************/
    _vector_sub((*r).s2,(*t).s2,psi1);
    _vector_sub((*r).s3,(*t).s3,psi2);
    /******************************** end of loop *********************************/
  }
  return;
}

/**************************************************************
 *
 * clover applies the clover term to spinor k, adds k 
 * to j then and stores it in l
 *
 * it is assumed that the clover leaf is computed and stored
 * in sw
 * the corresponding routine can be found in clover_leaf.c
 *
 **************************************************************/


void clover(const int ieo, 
	    spinor * const l, spinor * const k, spinor * const j,
	    const double q_off) {
  int ix;
  int ioff,icx;
  su3 *w1,*w2,*w3;
  spinor *r,*s,*t;
  static su3_vector chi, psi1, psi2;
  
  if(ieo == 0) {
    ioff = 0;
  } 
  else {
    ioff = (VOLUME+RAND)/2;
  }
  /************************ loop over all lattice sites *************************/
  for(icx = ioff; icx < (VOLUME/2+ioff); icx++) {
    ix = g_eo2lexic[icx];
    
    r = l + icx-ioff;
    s = k + icx-ioff;
    t = j + icx-ioff;

    w1=&sw[ix][0][0];
    w2=w1+2; /*&sw[ix][1][0];*/
    w3=w1+4; /*&sw[ix][2][0];*/
    _su3_multiply(psi1,*w1,(*s).s0); 
    _su3_multiply(chi,*w2,(*s).s1);
    _vector_add_assign(psi1,chi);
    _su3_inverse_multiply(psi2,*w2,(*s).s0); 
    _su3_multiply(chi,*w3,(*s).s1);
    _vector_add_assign(psi2,chi); 
    /*********************** contribution from q_off ***************************/
    _vector_add_mul(psi1,q_off,(*s).s0);
    _vector_add_mul(psi2,q_off,(*s).s1);
    /*************************************************/
    _vector_sub((*r).s0,psi1,(*t).s0);
    _vector_sub((*r).s1,psi2,(*t).s1);

    w1++; /*=&sw[ix][0][1];*/
    w2++; /*=&sw[ix][1][1];*/
    w3++; /*=&sw[ix][2][1];*/
    _su3_multiply(psi1,*w1,(*s).s2); 
    _su3_multiply(chi,*w2,(*s).s3);
    _vector_add_assign(psi1,chi); 
    _su3_inverse_multiply(psi2,*w2,(*s).s2); 
    _su3_multiply(chi,*w3,(*s).s3);
    _vector_add_assign(psi2,chi); 
    /*********************** contribution from q_off ***************************/
    _vector_add_mul(psi1,q_off,(*s).s2);
    _vector_add_mul(psi2,q_off,(*s).s3);
    /*************************************************/
    _vector_sub((*r).s2,psi1,(*t).s2);
    _vector_sub((*r).s3,psi2,(*t).s3);
    /******************************** end of loop *********************************/
  }
  return;
}

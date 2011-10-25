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
#include "clover.h"

/****
 *
 *
 *
 ****/

void clover_inv(const int ieo, spionor * const l) {
  int ix;
  int ioff,icx;
  su3 *w1,*w2,*w3;
  spinor *rn;

  if(ieo==0) {
    ioff=0;
  } 
  else {
    ioff=(VOLUME+RAND)/2;
  }
  /************************ loop over all lattice sites *************************/

  for(icx = ioff; icx < (VOLUME/2+ioff); icx++) {
    ix = g_eo2lexic[icx];

    rn = l + icx - ioff;
    /* noch schlimmer murks fuer den clover-term*/
    _vector_assign(phi1,(*rn).c0);
    _vector_assign(phi3,(*rn).c2);
    
    w1=&sw_inv[ix][0][0];
    w2=w1+2;  /* &sw_inv[ix][1][0]; */
    w3=w1+4;  /* &sw_inv[ix][2][0]; */
    _su3_multiply(psi,*w1,phi1); 
    _su3_multiply(chi,*w2,(*rn).c1);
    _vector_add((*rn).c0,psi,chi);
    _su3_inverse_multiply(psi,*w2,phi1); 
    _su3_multiply(chi,*w3,(*rn).c1);
    _vector_add((*rn).c1,psi,chi);
    
    w1++; /* &sw_inv[ix][0][1]; */
    w2++; /* &sw_inv[ix][1][1]; */
    w3++; /* &sw_inv[ix][2][1]; */
    _su3_multiply(psi,*w1,phi3); 
    _su3_multiply(chi,*w2,(*rn).c3);
    _vector_add((*rn).c2,psi,chi);
    _su3_inverse_multiply(psi,*w2,phi3); 
    _su3_multiply(chi,*w3,(*rn).c3);
    _vector_add((*rn).c3,psi,chi);
    /******************************** end of loop *********************************/
  }
  return;
}

void clover_gamma5(const int ieo, 
		   spinor * const l, spinor * const k, spinor * const j,
		   const double q_off) {
  int ix;
  int ioff,icx;
  su3 *w1,*w2,*w3;
  spinor *r,*s,*t;

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
    _su3_multiply(psi1,*w1,(*s).c0); 
    _su3_multiply(chi,*w2,(*s).c1);
    _vector_add_assign(psi1,chi);
    _su3_inverse_multiply(psi2,*w2,(*s).c0); 
    _su3_multiply(chi,*w3,(*s).c1);
    _vector_add_assign(psi2,chi); 
    /*********************** contribution from q_off ***************************/
    _vector_add_mul(psi1,q_off,(*s).c0);
    _vector_add_mul(psi2,q_off,(*s).c1);
    /**************************************************************/
    _vector_sub((*r).c0,psi1,(*t).c0);
    _vector_sub((*r).c1,psi2,(*t).c1);
    
    w1++; /*=&sw[ix][0][1];*/
    w2++; /*=&sw[ix][1][1];*/
    w3++; /*=&sw[ix][2][1];*/
    _su3_multiply(psi1,*w1,(*s).c2); _su3_multiply(chi,*w2,(*s).c3);
    _vector_add_assign(psi1,chi); 
    _su3_inverse_multiply(psi2,*w2,(*s).c2); _su3_multiply(chi,*w3,(*s).c3);
    _vector_add_assign(psi2,chi); 
    /*********************** contribution from q_off ***************************/
    _vector_add_mul(psi1,q_off,(*s).c2);
    _vector_add_mul(psi2,q_off,(*s).c3);
    /**************** multiply with  gamma5 included ******************************/
    _vector_sub((*r).c2,(*t).c2,psi1);
    _vector_sub((*r).c3,(*t).c3,psi2);
    /******************************** end of loop *********************************/
  }
  return;
}

void clover(const int ieo, 
	    spinor * const l, spinor * const k, spinor * const j,
	    const double q_off) {
  int ix;
  int ioff,icx;
  su3 *w1,*w2,*w3;
  spinor *r,*s,*t;
  
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
    _su3_multiply(psi1,*w1,(*s).c0); 
    _su3_multiply(chi,*w2,(*s).c1);
    _vector_add_assign(psi1,chi);
    _su3_inverse_multiply(psi2,*w2,(*s).c0); 
    _su3_multiply(chi,*w3,(*s).c1);
    _vector_add_assign(psi2,chi); 
    /*********************** contribution from q_off ***************************/
    _vector_add_mul(psi1,q_off,(*s).c0);
    _vector_add_mul(psi2,q_off,(*s).c1);
    /*************************************************/
    _vector_sub((*r).c0,psi1,(*t).c0);
    _vector_sub((*r).c1,psi2,(*t).c1);

    w1++; /*=&sw[ix][0][1];*/
    w2++; /*=&sw[ix][1][1];*/
    w3++; /*=&sw[ix][2][1];*/
    _su3_multiply(psi1,*w1,(*s).c2); 
    _su3_multiply(chi,*w2,(*s).c3);
    _vector_add_assign(psi1,chi); 
    _su3_inverse_multiply(psi2,*w2,(*s).c2); 
    _su3_multiply(chi,*w3,(*s).c3);
    _vector_add_assign(psi2,chi); 
    /*********************** contribution from q_off ***************************/
    _vector_add_mul(psi1,q_off,(*s).c2);
    _vector_add_mul(psi2,q_off,(*s).c3);
    /*************************************************/
    _vector_sub((*r).c2,psi1,(*t).c2);
    _vector_sub((*r).c3,psi2,(*t).c3);
    /******************************** end of loop *********************************/
  }
  return;
}

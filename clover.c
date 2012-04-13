/***********************************************************************
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
#include "linalg_eo.h"
#include "Hopping_Matrix.h"
#include "tm_operators.h"
#include "clover.h"


su3 *** sw;
su3 *** sw_inv;

void clover_gamma5(const int ieo, 
		   spinor * const l, spinor * const k, spinor * const j,
		   const double mu);
void clover(const int ieo, 
	    spinor * const l, spinor * const k, spinor * const j,
	    const double mu);

void Msw_full(spinor * const Even_new, spinor * const Odd_new, 
	      spinor * const Even, spinor * const Odd) {
  /* Even sites */
  Hopping_Matrix(EO, g_spinor_field[DUM_DERI], Odd);
  assign_mul_one_sw_pm_imu(EE, Even_new, Even, +g_mu);
  assign_add_mul_r(Even_new, g_spinor_field[DUM_DERI], -1., VOLUME/2);
  
  /* Odd sites */
  Hopping_Matrix(OE, g_spinor_field[DUM_DERI], Even);
  assign_mul_one_sw_pm_imu(OO, Odd_new, Odd, +g_mu);
  assign_add_mul_r(Odd_new, g_spinor_field[DUM_DERI], -1., VOLUME/2);
}


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


// this is the clover Qhat with mu = 0
void Qsw_psi(spinor * const l, spinor * const k) {
  Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX+1], k);
  clover_inv(EE, g_spinor_field[DUM_MATRIX+1], 0.);
  Hopping_Matrix(OE, g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+1]);
  clover_gamma5(OO, l, k, g_spinor_field[DUM_MATRIX], 0.);
}

// this is the twisted clover Qhat with -mu
void Qsw_minus_psi(spinor * const l, spinor * const k) {
  Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX+1], k);
  clover_inv(EE, g_spinor_field[DUM_MATRIX+1], -g_mu);
  Hopping_Matrix(OE, g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+1]);
  clover_gamma5(OO, l, k, g_spinor_field[DUM_MATRIX], -(g_mu + g_mu3));
}

// this is the twisted clover Qhat with +mu
void Qsw_plus_psi(spinor * const l, spinor * const k) {
  Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX+1], k);
  clover_inv(EE, g_spinor_field[DUM_MATRIX+1], +g_mu);
  Hopping_Matrix(OE, g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+1]);
  clover_gamma5(OO, l, k, g_spinor_field[DUM_MATRIX], +(g_mu + g_mu3));
}


void Qsw_sq_psi(spinor * const l, spinor * const k) {
  /* \hat Q_{-} */
  Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX+1], k);
  clover_inv(EE, g_spinor_field[DUM_MATRIX+1], 0.);
  Hopping_Matrix(OE, g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+1]);
  clover_gamma5(OO, g_spinor_field[DUM_MATRIX], k, g_spinor_field[DUM_MATRIX], 0.);
  /* \hat Q_{+} */
  Hopping_Matrix(EO, l, g_spinor_field[DUM_MATRIX]);
  clover_inv(EE, l, 0.); 
  Hopping_Matrix(OE, g_spinor_field[DUM_MATRIX+1], l);
  clover_gamma5(OO, l, g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+1], 0.);
}

void Qsw_pm_psi(spinor * const l, spinor * const k) {
  /* \hat Q_{-} */
  Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX+1], k);
  clover_inv(EE, g_spinor_field[DUM_MATRIX+1], -g_mu);
  Hopping_Matrix(OE, g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+1]);
  clover_gamma5(OO, g_spinor_field[DUM_MATRIX], k, g_spinor_field[DUM_MATRIX], -(g_mu + g_mu3));
  /* \hat Q_{+} */
  Hopping_Matrix(EO, l, g_spinor_field[DUM_MATRIX]);
  clover_inv(EE, l, +g_mu); 
  Hopping_Matrix(OE, g_spinor_field[DUM_MATRIX+1], l);
  clover_gamma5(OO, l, g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+1], +(g_mu + g_mu3));
}

// this is the clover Mhat with mu = 0
void Msw_psi(spinor * const l, spinor * const k) {
  Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX+1], k);
  clover_inv(EE, g_spinor_field[DUM_MATRIX+1], 0.);
  Hopping_Matrix(OE, g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+1]);
  clover(OO, l, k, g_spinor_field[DUM_MATRIX], 0.);
}

void Msw_plus_psi(spinor * const l, spinor * const k) {
  Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX+1], k);
  clover_inv(EE, g_spinor_field[DUM_MATRIX+1], +g_mu);
  Hopping_Matrix(OE, g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+1]);
  clover(OO, l, k, g_spinor_field[DUM_MATRIX], +(g_mu + g_mu3));
}

void Msw_minus_psi(spinor * const l, spinor * const k) {
  Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX+1], k);
  clover_inv(EE, g_spinor_field[DUM_MATRIX+1], -g_mu);
  Hopping_Matrix(OE, g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+1]);
  clover(OO, l, k, g_spinor_field[DUM_MATRIX], -(g_mu + g_mu3));
}


void H_eo_sw_inv_psi(spinor * const l, spinor * const k, const int ieo, const double mu) {
  Hopping_Matrix(ieo, l, k);
  clover_inv(ieo, l, mu);
  return;
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

void clover_inv(const int ieo, spinor * const l, const double mu) {
#ifdef OMP
#pragma omp parallel
  {
  su3_vector psi, chi, phi1, phi3;
  int icy;
#else
  static su3_vector psi, chi, phi1, phi3;
#endif
  int ioff = 0;
  su3 *w1, *w2, *w3, *w4;
  spinor *rn;

  if(mu < 0) ioff = VOLUME/2;
  /************************ loop over all lattice sites *************************/
#ifdef OMP
#pragma omp for
  for(int icx = 0; icx < (VOLUME/2); icx++) {
    icy = ioff + icx;
#else
  for(int icx = 0, icy = ioff; icx < (VOLUME/2); icx++, icy++) {
#endif

    rn = l + icx;
    _vector_assign(phi1,(*rn).s0);
    _vector_assign(phi3,(*rn).s2);

    w1=&sw_inv[icy][0][0];
    w2=w1+2;  /* &sw_inv[icy][1][0]; */
    w3=w1+4;  /* &sw_inv[icy][2][0]; */
    w4=w1+6;  /* &sw_inv[icy][3][0]; */
    _su3_multiply(psi,*w1,phi1); 
    _su3_multiply(chi,*w2,(*rn).s1);
    _vector_add((*rn).s0,psi,chi);
    _su3_multiply(psi,*w4,phi1); 
    _su3_multiply(chi,*w3,(*rn).s1);
    _vector_add((*rn).s1,psi,chi);

    w1++; /* &sw_inv[icy][0][1]; */
    w2++; /* &sw_inv[icy][1][1]; */
    w3++; /* &sw_inv[icy][2][1]; */
    w4++; /* &sw_inv[icy][3][1]; */
    _su3_multiply(psi,*w1,phi3); 
    _su3_multiply(chi,*w2,(*rn).s3);
    _vector_add((*rn).s2,psi,chi);
    _su3_multiply(psi,*w4,phi3); 
    _su3_multiply(chi,*w3,(*rn).s3);
    _vector_add((*rn).s3,psi,chi);

    /******************************** end of loop *********************************/
  }
#ifdef OMP
  } /* OpenMP closing brace */
#endif
  return;
}

/**************************************************************
 *
 * clover_gamma5 applies the clover term to spinor k, adds k 
 * to j then and stores it in l multiplied by gamma_5
 *
 * it is assumed that the clover leaf is computed and stored
 * in sw[VOLUME][3][2]
 * the corresponding routine can be found in clover_leaf.c
 *
 **************************************************************/

void clover_gamma5(const int ieo, 
		   spinor * const l, spinor * const k, spinor * const j,
		   const double mu) {
#ifdef OMP
#pragma omp parallel
  {
  su3_vector chi, psi1, psi2;
#else
  static su3_vector chi, psi1, psi2;
#endif
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
#ifdef OMP
#pragma omp for
#endif
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
    // add in the twisted mass term (plus in the upper components)
    _vector_add_i_mul(psi1, mu, (*s).s0);
    _vector_add_i_mul(psi2, mu, (*s).s1);

    _vector_sub((*r).s0,psi1,(*t).s0);
    _vector_sub((*r).s1,psi2,(*t).s1);
    
    w1++; /*=&sw[ix][0][1];*/
    w2++; /*=&sw[ix][1][1];*/
    w3++; /*=&sw[ix][2][1];*/
    _su3_multiply(psi1,*w1,(*s).s2); _su3_multiply(chi,*w2,(*s).s3);
    _vector_add_assign(psi1,chi); 
    _su3_inverse_multiply(psi2,*w2,(*s).s2); _su3_multiply(chi,*w3,(*s).s3);
    _vector_add_assign(psi2,chi); 
    // add in the twisted mass term (minus from g5 in the lower components)
    _vector_add_i_mul(psi1, -mu, (*s).s2);
    _vector_add_i_mul(psi2, -mu, (*s).s3);

    /**************** multiply with  gamma5 included ******************************/
    _vector_sub((*r).s2,(*t).s2,psi1);
    _vector_sub((*r).s3,(*t).s3,psi2);
    /******************************** end of loop *********************************/
  }
#ifdef OMP
  } /* OMP closing brace */
#endif
  return;
}

/**************************************************************
 *
 * clover applies the clover term to spinor k, adds k 
 * to j then and stores it in l
 *
 * it is assumed that the clover leaf is computed and stored
 * in sw[VOLUME][3][2]
 * the corresponding routine can be found in clover_leaf.c
 *
 **************************************************************/


void clover(const int ieo, 
	    spinor * const l, spinor * const k, spinor * const j,
	    const double mu) {
#ifdef OMP
#pragma omp parallel
  {
  su3_vector chi, psi1, psi2;
#else
  static su3_vector chi, psi1, psi2;
#endif
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
#ifdef OMP
#pragma omp for
#endif
  for(icx = ioff; icx < (VOLUME/2+ioff); icx++) {
    ix = g_eo2lexic[icx];
    
    r = l + icx-ioff;
    s = k + icx-ioff;
    t = j + icx-ioff;

    // upper two spin components first
    w1=&sw[ix][0][0];
    w2=w1+2; /*&sw[ix][1][0];*/
    w3=w1+4; /*&sw[ix][2][0];*/
    _su3_multiply(psi1,*w1,(*s).s0); 
    _su3_multiply(chi,*w2,(*s).s1);
    _vector_add_assign(psi1,chi);
    _su3_inverse_multiply(psi2,*w2,(*s).s0); 
    _su3_multiply(chi,*w3,(*s).s1);
    _vector_add_assign(psi2,chi); 

    // add in the twisted mass term (plus in the upper components)
    _vector_add_i_mul(psi1, mu, (*s).s0);
    _vector_add_i_mul(psi2, mu, (*s).s1);

    _vector_sub((*r).s0,psi1,(*t).s0);
    _vector_sub((*r).s1,psi2,(*t).s1);

    // now lower to spin components
    w1++; /*=&sw[ix][0][1];*/
    w2++; /*=&sw[ix][1][1];*/
    w3++; /*=&sw[ix][2][1];*/
    _su3_multiply(psi1,*w1,(*s).s2); 
    _su3_multiply(chi,*w2,(*s).s3);
    _vector_add_assign(psi1,chi); 
    _su3_inverse_multiply(psi2,*w2,(*s).s2); 
    _su3_multiply(chi,*w3,(*s).s3);
    _vector_add_assign(psi2,chi); 

    // add in the twisted mass term (minus from g5 in the lower components)
    _vector_add_i_mul(psi1, -mu, (*s).s2);
    _vector_add_i_mul(psi2, -mu, (*s).s3);

    _vector_sub((*r).s2,psi1,(*t).s2);
    _vector_sub((*r).s3,psi2,(*t).s3);
  }
#ifdef OMP
  } /* OpenMP closing brace */
#endif
  return;
}

void assign_mul_one_sw_pm_imu(const int ieo, 
			      spinor * const k, spinor * const l,
			      const double mu) {
  int ix;
  int ioff, icx;
  su3 *w1, *w2, *w3;
  spinor *r, *s;
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
    
    r = k + icx-ioff;
    s = l + icx-ioff;

    // upper two spin components first
    w1=&sw[ix][0][0];
    w2=w1+2; /*&sw[ix][1][0];*/
    w3=w1+4; /*&sw[ix][2][0];*/
    _su3_multiply(psi1,*w1,(*s).s0); 
    _su3_multiply(chi,*w2,(*s).s1);
    _vector_add_assign(psi1,chi);
    _su3_inverse_multiply(psi2,*w2,(*s).s0); 
    _su3_multiply(chi,*w3,(*s).s1);
    _vector_add_assign(psi2,chi); 

    // add in the twisted mass term (plus in the upper components)
    _vector_add_i_mul(psi1, mu, (*s).s0);
    _vector_add_i_mul(psi2, mu, (*s).s1);

    _vector_assign((*r).s0, psi1);
    _vector_assign((*r).s1, psi2);

    // now lower to spin components
    w1++; /*=&sw[ix][0][1];*/
    w2++; /*=&sw[ix][1][1];*/
    w3++; /*=&sw[ix][2][1];*/
    _su3_multiply(psi1,*w1,(*s).s2); 
    _su3_multiply(chi,*w2,(*s).s3);
    _vector_add_assign(psi1,chi); 
    _su3_inverse_multiply(psi2,*w2,(*s).s2); 
    _su3_multiply(chi,*w3,(*s).s3);
    _vector_add_assign(psi2,chi); 

    // add in the twisted mass term (minus from g5 in the lower components)
    _vector_add_i_mul(psi1, -mu, (*s).s2);
    _vector_add_i_mul(psi2, -mu, (*s).s3);

    _vector_assign((*r).s2, psi1);
    _vector_assign((*r).s3, psi2);
  }
  return;
}


void assign_mul_one_sw_pm_imu_inv(const int ieo, 
				  spinor * const k, spinor * const l,
				  const double mu) {
  su3 *w1, *w2, *w3, *w4;
  spinor *rn, *s;
  static su3_vector psi, chi, phi1, phi3;

  /************************ loop over all lattice sites *************************/

  for(int icx = 0; icx < (VOLUME/2); icx++) {

    rn = l + icx;
    s = k + icx;
    _vector_assign(phi1,(*rn).s0);
    _vector_assign(phi3,(*rn).s2);

    w1=&sw_inv[icx][0][0];
    w2=w1+2;  /* &sw_inv[icx][1][0]; */
    w3=w1+4;  /* &sw_inv[icx][2][0]; */
    w4=w1+6;  /* &sw_inv[icx][3][0]; */
    _su3_multiply(psi,*w1,phi1); 
    _su3_multiply(chi,*w2,(*rn).s1);
    _vector_add((*s).s0,psi,chi);
    _su3_multiply(psi,*w4,phi1); 
    _su3_multiply(chi,*w3,(*rn).s1);
    _vector_add((*s).s1,psi,chi);

    w1++; /* &sw_inv[icx][0][1]; */
    w2++; /* &sw_inv[icx][1][1]; */
    w3++; /* &sw_inv[icx][2][1]; */
    w4++; /* &sw_inv[icx][3][1]; */
    _su3_multiply(psi,*w1,phi3); 
    _su3_multiply(chi,*w2,(*rn).s3);
    _vector_add((*s).s2,psi,chi);
    _su3_multiply(psi,*w4,phi3); 
    _su3_multiply(chi,*w3,(*rn).s3);
    _vector_add((*s).s3,psi,chi);

    /******************************** end of loop *********************************/
  }
  return;
}


/********
 *
 * temporary initialisation function
 *
 ********/

su3 ** sw1, ** sw_inv1;
su3 * _sw, *_sw_inv;

void init_sw_fields() {
  int V = VOLUME;
  su3 * tmp;
  static int sw_init = 0;

  if(!sw_init) {
    if((void*)(sw = (su3***)calloc(V, sizeof(su3**))) == NULL) {
      fprintf (stderr, "sw malloc err\n"); 
    }
    if((void*)(sw_inv = (su3***)calloc(V, sizeof(su3**))) == NULL) {
      fprintf (stderr, "sw_inv malloc err\n"); 
    }
    if((void*)(sw1 = (su3**)calloc(3*V, sizeof(su3*))) == NULL) {
      fprintf (stderr, "sw1 malloc err\n"); 
    }
    if((void*)(sw_inv1 = (su3**)calloc(4*V, sizeof(su3*))) == NULL) {
      fprintf (stderr, "sw_inv1 malloc err\n"); 
    }
    if((void*)(_sw = (su3*)calloc(3*2*V+1, sizeof(su3))) == NULL) {
      fprintf (stderr, "_sw malloc err\n"); 
    }
    if((void*)(_sw_inv = (su3*)calloc(4*2*V+1, sizeof(su3))) == NULL) {
      fprintf (stderr, "_sw_inv malloc err\n"); 
    }
    sw[0] = sw1;
    sw_inv[0] = sw_inv1;
    for(int i = 1; i < V; i++) {
      sw[i] = sw[i-1]+3;
      sw_inv[i] = sw_inv[i-1]+4;
    }
#    if (defined SSE || defined SSE2 || defined SSE3)
    sw[0][0] = (su3*)(((unsigned long int)(_sw)+ALIGN_BASE)&~ALIGN_BASE);
    sw_inv[0][0] = (su3*)(((unsigned long int)(_sw_inv)+ALIGN_BASE)&~ALIGN_BASE);
#    else
    sw[0][0] = _sw;
    sw_inv[0][0] = _sw_inv;
#    endif
    tmp = sw[0][0];
    for(int i = 0; i < V; i++) {
      for(int j = 0; j < 3; j++) {
	sw[i][j] = tmp;
	tmp = tmp+2;
      }
    }
    
    tmp = sw_inv[0][0];
    for(int i = 0; i < V; i++) {
      for(int j = 0; j < 4; j++) {
	sw_inv[i][j] = tmp;
	tmp = tmp+2;
      }
    }
    sw_init = 1;
  }
  return;
}

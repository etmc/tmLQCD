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
#include "operator/Hopping_Matrix.h"
#include "tm_operators.h"
#include "operator/clovertm_operators.h"


su3 *** sw;
su3 *** sw_inv;

void clover_gamma5(const int ieo, 
		   spinor * const l, const spinor * const k, const spinor * const j,
		   const double mu);
void clover(const int ieo, 
	    spinor * const l, const spinor * const k, const spinor * const j,
	    const double mu);

void Msw_full(spinor * const Even_new, spinor * const Odd_new, 
	      spinor * const Even, spinor * const Odd) {
  /* Even sites */
  Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX], Odd);
  assign_mul_one_sw_pm_imu(EE, Even_new, Even, +g_mu);
  assign_add_mul_r(Even_new, g_spinor_field[DUM_MATRIX], -1., VOLUME/2);
  
  /* Odd sites */
  Hopping_Matrix(OE, g_spinor_field[DUM_MATRIX], Even);
  assign_mul_one_sw_pm_imu(OO, Odd_new, Odd, +g_mu);
  assign_add_mul_r(Odd_new, g_spinor_field[DUM_MATRIX], -1., VOLUME/2);
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
  clover_inv(g_spinor_field[DUM_MATRIX+1], +1, g_mu);
  Hopping_Matrix(OE, g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+1]);
  clover_gamma5(OO, l, k, g_spinor_field[DUM_MATRIX], 0.);
}

// this is the twisted clover Qhat with -mu
void Qsw_minus_psi(spinor * const l, spinor * const k) {
  Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX+1], k);
  clover_inv(g_spinor_field[DUM_MATRIX+1], -1, g_mu);
  Hopping_Matrix(OE, g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+1]);
  clover_gamma5(OO, l, k, g_spinor_field[DUM_MATRIX], -(g_mu + g_mu3));
}

// this is the twisted clover Qhat with +mu
void Qsw_plus_psi(spinor * const l, spinor * const k) {
  Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX+1], k);
  clover_inv(g_spinor_field[DUM_MATRIX+1], +1, g_mu);
  Hopping_Matrix(OE, g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+1]);
  clover_gamma5(OO, l, k, g_spinor_field[DUM_MATRIX], +(g_mu + g_mu3));
}


void Qsw_sq_psi(spinor * const l, spinor * const k) {
  /* \hat Q_{-} */
  Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX+1], k);
  clover_inv(g_spinor_field[DUM_MATRIX+1], +1, g_mu);
  Hopping_Matrix(OE, g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+1]);
  clover_gamma5(OO, g_spinor_field[DUM_MATRIX], k, g_spinor_field[DUM_MATRIX], 0.);
  /* \hat Q_{+} */
  Hopping_Matrix(EO, l, g_spinor_field[DUM_MATRIX]);
  clover_inv(l, +1, g_mu); 
  Hopping_Matrix(OE, g_spinor_field[DUM_MATRIX+1], l);
  clover_gamma5(OO, l, g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+1], 0.);
}

void Qsw_pm_psi(spinor * const l, spinor * const k) {
  /* \hat Q_{-} */
  Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX+1], k);
  clover_inv(g_spinor_field[DUM_MATRIX+1], -1, g_mu);
  Hopping_Matrix(OE, g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+1]);
  clover_gamma5(OO, g_spinor_field[DUM_MATRIX], k, g_spinor_field[DUM_MATRIX], -(g_mu + g_mu3));
  /* \hat Q_{+} */
  Hopping_Matrix(EO, l, g_spinor_field[DUM_MATRIX]);
  clover_inv(l, +1, g_mu); 
  Hopping_Matrix(OE, g_spinor_field[DUM_MATRIX+1], l);
  clover_gamma5(OO, l, g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+1], +(g_mu + g_mu3));
}

// this is the clover Mhat with mu = 0
void Msw_psi(spinor * const l, spinor * const k) {
  Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX+1], k);
  clover_inv(g_spinor_field[DUM_MATRIX+1], +1, g_mu);
  Hopping_Matrix(OE, g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+1]);
  clover(OO, l, k, g_spinor_field[DUM_MATRIX], 0.);
}

void Msw_plus_psi(spinor * const l, spinor * const k) {
  Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX+1], k);
  clover_inv(g_spinor_field[DUM_MATRIX+1], +1, g_mu);
  Hopping_Matrix(OE, g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+1]);
  clover(OO, l, k, g_spinor_field[DUM_MATRIX], +(g_mu + g_mu3));
}

void Msw_minus_psi(spinor * const l, spinor * const k) {
  Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX+1], k);
  clover_inv(g_spinor_field[DUM_MATRIX+1], -1, g_mu);
  Hopping_Matrix(OE, g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+1]);
  clover(OO, l, k, g_spinor_field[DUM_MATRIX], -(g_mu + g_mu3));
}


void H_eo_sw_inv_psi(spinor * const l, spinor * const k, const int ieo, const int tau3sign, const double mu) {
  Hopping_Matrix(ieo, l, k);
  clover_inv(l, tau3sign, mu);
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

void clover_inv(spinor * const l, const int tau3sign, const double mu) {
#ifdef OMP
#pragma omp parallel
  {
#endif
  int icy;
  su3_vector ALIGN psi, chi, phi1, phi3;
  int ioff = 0;
  const su3 *w1, *w2, *w3, *w4;
  spinor *rn;

  if(tau3sign < 0 && fabs(mu) > 0) {
    ioff = VOLUME/2;
  }

#ifndef OMP
  icy = ioff;
#endif
  /************************ loop over all lattice sites *************************/
#ifdef OMP
#pragma omp for
#endif
  for(int icx = 0; icx < (VOLUME/2); icx++) {
#ifdef OMP
    icy = ioff + icx;
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

#ifndef OMP
    ++icy;
#endif

    /******************************** end of loop *********************************/
  }
#ifdef OMP
  } /* OpenMP closing brace */
#endif
  return;
}

void clover_inv_nd(const int ieo, spinor * const l_c, spinor * const l_s) {
#ifdef OMP
#pragma omp parallel
  {
#endif
  int icy;
  su3_vector ALIGN psi, chi, phi1, phi3;
  int ioff = 0;
  const su3 *w1, *w2, *w3, *w4;
  spinor *rn_s, *rn_c;


  if(ieo == 1) ioff = VOLUME/2;

#ifndef OMP
  icy = ioff;
#endif

#ifdef OMP
#pragma omp for
#endif
  for(unsigned int icx = 0; icx < (VOLUME/2); icx++) {
#ifdef OMP
    icy = ioff + icx;
#endif

    rn_s = l_s + icx;
    rn_c = l_c + icx;
    _vector_assign(phi1,(*rn_s).s0);

    w1=&sw_inv[icy][0][0];
    w2=w1+2;  /* &sw_inv[icy][1][0]; */
    w3=w1+4;  /* &sw_inv[icy][2][0]; */
    w4=w1+6;  /* &sw_inv[icy][3][0]; */
    _su3_multiply(psi, *w1, phi1); 
    _su3_multiply(chi, *w2, (*rn_s).s1);
    _vector_add((*rn_s).s0, psi,chi);
    _su3_multiply(psi, *w4, phi1); 
    _su3_multiply(chi, *w3, (*rn_s).s1);
    _vector_add((*rn_s).s1, psi, chi);

    _vector_assign(phi1,(*rn_c).s0);

    _su3_multiply(psi, *w1, phi1); 
    _su3_multiply(chi, *w2, (*rn_c).s1);
    _vector_add((*rn_c).s0, psi,chi);
    _su3_multiply(psi, *w4, phi1); 
    _su3_multiply(chi, *w3, (*rn_c).s1);
    _vector_add((*rn_c).s1, psi, chi);

    _vector_assign(phi3,(*rn_s).s2);

    w1++; /* &sw_inv[icy][0][1]; */
    w2++; /* &sw_inv[icy][1][1]; */
    w3++; /* &sw_inv[icy][2][1]; */
    w4++; /* &sw_inv[icy][3][1]; */
    _su3_multiply(psi, *w1, phi3); 
    _su3_multiply(chi, *w2, (*rn_s).s3);
    _vector_add((*rn_s).s2, psi, chi);
    _su3_multiply(psi, *w4, phi3); 
    _su3_multiply(chi, *w3, (*rn_s).s3);
    _vector_add((*rn_s).s3, psi, chi);

    _vector_assign(phi3,(*rn_c).s2);

    _su3_multiply(psi, *w1, phi3); 
    _su3_multiply(chi, *w2, (*rn_c).s3);
    _vector_add((*rn_c).s2, psi, chi);
    _su3_multiply(psi, *w4, phi3); 
    _su3_multiply(chi, *w3, (*rn_c).s3);
    _vector_add((*rn_c).s3, psi, chi);

#ifndef OMP
    ++icy;
#endif

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
		   spinor * const l, const spinor * const k, const spinor * const j,
		   const double mu) {
#ifdef OMP
#pragma omp parallel
  {
#endif
  su3_vector ALIGN chi, psi1, psi2;
  int ix;
  int ioff,icx;
  const su3 *w1,*w2,*w3;
  spinor *r;
  const spinor *s,*t;

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
 * clover applies (1 + T + imug5) to spinor k, 
 * subtracts j from k and stores in l
 *
 * it is assumed that the clover leaf is computed and stored
 * in sw[VOLUME][3][2]
 * the corresponding routine can be found in clover_leaf.c
 *
 **************************************************************/


void clover(const int ieo, 
	    spinor * const l, const spinor * const k, const spinor * const j,
	    const double mu) {
#ifdef OMP
#pragma omp parallel
  {
#endif
  su3_vector ALIGN chi, psi1, psi2;
  int ix;
  int ioff;
  const su3 *w1,*w2,*w3;
  spinor *r;
  const spinor *s,*t;
  
  if(ieo == 0) {
    ioff = 0;
  } 
  else {
    ioff = (VOLUME+RAND)/2;
  }
#ifdef OMP
#pragma omp for
#endif
  for(unsigned int icx = ioff; icx < (VOLUME/2+ioff); icx++) {
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

/**************************************************************
 *
 * clover_nd applies the clover (1 + T + imug5tau3 + epstau1) 
 * term to spinor k, subtracts j from k and stores in l
 *
 * it is assumed that the clover leaf is computed and stored
 * in sw[VOLUME][3][2]
 * the corresponding routine can be found in clover_leaf.c
 *
 **************************************************************/

void clover_nd(const int ieo, 
	       spinor * const l_c, spinor * const l_s, 
	       const spinor * const k_c, const spinor * const k_s, 
	       const spinor * const j_c, const spinor * const j_s,
	       const double mubar, const double epsbar) {
#ifdef OMP
#pragma omp parallel
  {
#endif
  su3_vector ALIGN chi, psi1, psi2;
  int ix;
  int ioff;
  const su3 *w1,*w2,*w3;
  spinor *r_s, *r_c;
  const spinor *s_s, *s_c, *t_s, *t_c;
  
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
  for(unsigned int icx = ioff; icx < (VOLUME/2+ioff); icx++) {
    ix = g_eo2lexic[icx];
    
    r_s = l_s + icx-ioff;
    r_c = l_c + icx-ioff;
    s_s = k_s + icx-ioff;
    s_c = k_c + icx-ioff;
    t_s = j_s + icx-ioff;
    t_c = j_c + icx-ioff;

    // upper two spin components first
    w1=&sw[ix][0][0];
    w2=w1+2; /*&sw[ix][1][0];*/
    w3=w1+4; /*&sw[ix][2][0];*/
    _su3_multiply(psi1, *w1, (*s_s).s0); 
    _su3_multiply(chi, *w2, (*s_s).s1);
    _vector_add_assign(psi1, chi);
    _su3_inverse_multiply(psi2, *w2, (*s_s).s0); 
    _su3_multiply(chi, *w3, (*s_s).s1);
    _vector_add_assign(psi2, chi); 

    // add in the twisted mass term (plus in the upper components)
    _vector_add_i_mul(psi1, mubar, (*s_s).s0);
    _vector_add_i_mul(psi2, mubar, (*s_s).s1);

    _vector_add_mul(psi1, epsbar, (*s_c).s0);
    _vector_add_mul(psi2, epsbar, (*s_c).s1);

    _vector_sub((*r_s).s0, psi1, (*t_s).s0);
    _vector_sub((*r_s).s1, psi2, (*t_s).s1);

    _su3_multiply(psi1, *w1, (*s_c).s0); 
    _su3_multiply(chi, *w2, (*s_c).s1);
    _vector_add_assign(psi1, chi);
    _su3_inverse_multiply(psi2, *w2, (*s_c).s0); 
    _su3_multiply(chi, *w3, (*s_c).s1);
    _vector_add_assign(psi2, chi); 

    // add in the twisted mass term (plus in the upper components)
    _vector_add_i_mul(psi1, -mubar, (*s_c).s0);
    _vector_add_i_mul(psi2, -mubar, (*s_c).s1);

    _vector_add_mul(psi1, epsbar, (*s_s).s0);
    _vector_add_mul(psi2, epsbar, (*s_s).s1);

    _vector_sub((*r_c).s0, psi1, (*t_c).s0);
    _vector_sub((*r_c).s1, psi2, (*t_c).s1);


    // now lower to spin components
    w1++; /*=&sw[ix][0][1];*/
    w2++; /*=&sw[ix][1][1];*/
    w3++; /*=&sw[ix][2][1];*/
    _su3_multiply(psi1, *w1, (*s_s).s2); 
    _su3_multiply(chi, *w2, (*s_s).s3);
    _vector_add_assign(psi1, chi); 
    _su3_inverse_multiply(psi2, *w2, (*s_s).s2); 
    _su3_multiply(chi, *w3, (*s_s).s3);
    _vector_add_assign(psi2, chi); 

    // add in the twisted mass term (minus from g5 in the lower components)
    _vector_add_i_mul(psi1, -mubar, (*s_s).s2);
    _vector_add_i_mul(psi2, -mubar, (*s_s).s3);

    _vector_add_mul(psi1, epsbar, (*s_c).s2);
    _vector_add_mul(psi2, epsbar, (*s_c).s3);

    _vector_sub((*r_s).s2,psi1,(*t_s).s2);
    _vector_sub((*r_s).s3,psi2,(*t_s).s3);

    _su3_multiply(psi1, *w1, (*s_c).s2); 
    _su3_multiply(chi, *w2, (*s_c).s3);
    _vector_add_assign(psi1, chi); 
    _su3_inverse_multiply(psi2, *w2, (*s_c).s2); 
    _su3_multiply(chi, *w3, (*s_c).s3);
    _vector_add_assign(psi2, chi); 

    // add in the twisted mass term (minus from g5 in the lower components)
    _vector_add_i_mul(psi1, mubar, (*s_c).s2);
    _vector_add_i_mul(psi2, mubar, (*s_c).s3);

    _vector_add_mul(psi1, epsbar, (*s_s).s2);
    _vector_add_mul(psi2, epsbar, (*s_s).s3);

    _vector_sub((*r_c).s2, psi1, (*t_c).s2);
    _vector_sub((*r_c).s3, psi2, (*t_c).s3);
  }
#ifdef OMP
  } /* OpenMP closing brace */
#endif
  return;
}

void clover_gamma5_nd(const int ieo, 
		      spinor * const l_c, spinor * const l_s, 
		      const spinor * const k_c, const spinor * const k_s, 
		      const spinor * const j_c, const spinor * const j_s,
		      const double mubar, const double epsbar) {
#ifdef OMP
#pragma omp parallel
  {
#endif
  su3_vector ALIGN chi, psi1, psi2;
  int ix;
  int ioff;
  const su3 *w1,*w2,*w3;
  spinor *r_s, *r_c;
  const spinor *s_s, *s_c, *t_s, *t_c;
  
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
  for(unsigned int icx = ioff; icx < (VOLUME/2+ioff); icx++) {
    ix = g_eo2lexic[icx];
    
    r_s = l_s + icx-ioff;
    r_c = l_c + icx-ioff;
    s_s = k_s + icx-ioff;
    s_c = k_c + icx-ioff;
    t_s = j_s + icx-ioff;
    t_c = j_c + icx-ioff;

    // upper two spin components first
    w1=&sw[ix][0][0];
    w2=w1+2; /*&sw[ix][1][0];*/
    w3=w1+4; /*&sw[ix][2][0];*/
    _su3_multiply(psi1, *w1, (*s_s).s0); 
    _su3_multiply(chi, *w2, (*s_s).s1);
    _vector_add_assign(psi1, chi);
    _su3_inverse_multiply(psi2, *w2, (*s_s).s0); 
    _su3_multiply(chi, *w3, (*s_s).s1);
    _vector_add_assign(psi2, chi); 

    // add in the twisted mass term (plus in the upper components)
    _vector_add_i_mul(psi1, mubar, (*s_s).s0);
    _vector_add_i_mul(psi2, mubar, (*s_s).s1);

    _vector_add_mul(psi1, epsbar, (*s_c).s0);
    _vector_add_mul(psi2, epsbar, (*s_c).s1);

    _vector_sub((*r_s).s0, psi1, (*t_s).s0);
    _vector_sub((*r_s).s1, psi2, (*t_s).s1);

    _su3_multiply(psi1, *w1, (*s_c).s0); 
    _su3_multiply(chi, *w2, (*s_c).s1);
    _vector_add_assign(psi1, chi);
    _su3_inverse_multiply(psi2, *w2, (*s_c).s0); 
    _su3_multiply(chi, *w3, (*s_c).s1);
    _vector_add_assign(psi2, chi); 

    // add in the twisted mass term (plus in the upper components)
    _vector_add_i_mul(psi1, -mubar, (*s_c).s0);
    _vector_add_i_mul(psi2, -mubar, (*s_c).s1);

    _vector_add_mul(psi1, epsbar, (*s_s).s0);
    _vector_add_mul(psi2, epsbar, (*s_s).s1);

    _vector_sub((*r_c).s0, psi1, (*t_c).s0);
    _vector_sub((*r_c).s1, psi2, (*t_c).s1);


    // now lower to spin components
    w1++; /*=&sw[ix][0][1];*/
    w2++; /*=&sw[ix][1][1];*/
    w3++; /*=&sw[ix][2][1];*/
    _su3_multiply(psi1, *w1, (*s_s).s2); 
    _su3_multiply(chi, *w2, (*s_s).s3);
    _vector_add_assign(psi1, chi); 
    _su3_inverse_multiply(psi2, *w2, (*s_s).s2); 
    _su3_multiply(chi, *w3, (*s_s).s3);
    _vector_add_assign(psi2, chi); 

    // add in the twisted mass term (minus from g5 in the lower components)
    _vector_add_i_mul(psi1, -mubar, (*s_s).s2);
    _vector_add_i_mul(psi2, -mubar, (*s_s).s3);

    _vector_add_mul(psi1, epsbar, (*s_c).s2);
    _vector_add_mul(psi2, epsbar, (*s_c).s3);

    _vector_sub((*r_s).s2, (*t_s).s2, psi1);
    _vector_sub((*r_s).s3, (*t_s).s3, psi2);

    _su3_multiply(psi1, *w1, (*s_c).s2); 
    _su3_multiply(chi, *w2, (*s_c).s3);
    _vector_add_assign(psi1, chi); 
    _su3_inverse_multiply(psi2, *w2, (*s_c).s2); 
    _su3_multiply(chi, *w3, (*s_c).s3);
    _vector_add_assign(psi2, chi); 

    // add in the twisted mass term (minus from g5 in the lower components)
    _vector_add_i_mul(psi1, mubar, (*s_c).s2);
    _vector_add_i_mul(psi2, mubar, (*s_c).s3);

    _vector_add_mul(psi1, epsbar, (*s_s).s2);
    _vector_add_mul(psi2, epsbar, (*s_s).s3);

    _vector_sub((*r_c).s2, (*t_c).s2, psi1);
    _vector_sub((*r_c).s3, (*t_c).s3, psi2);
  }
#ifdef OMP
  } /* OpenMP closing brace */
#endif
  return;
}


/**************************************************************
 *
 * assign_mul_one_sw_pm_imu applies (1 + T + imug5) to spinor l
 * and stores it in k
 *
 * it is assumed that the clover leaf is computed and stored
 * in sw[VOLUME][3][2]
 * the corresponding routine can be found in clover_leaf.c
 *
 **************************************************************/


void assign_mul_one_sw_pm_imu(const int ieo, 
			      spinor * const k, const spinor * const l,
			      const double mu) {
#ifdef OMP
#pragma omp parallel
  {
#endif
  su3_vector ALIGN chi, psi1, psi2;
  int ix;
  int ioff;
  const su3 *w1, *w2, *w3;
  spinor *r;
  const spinor *s;
  
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
  for(unsigned icx = ioff; icx < (VOLUME/2+ioff); icx++) {
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
#ifdef OMP
  } /* OpenMP closing brace */
#endif
  return;
}

/**************************************************************
 *
 * assign_mul_one_sw_pm_imu_eps applies 
 * (1 + T + imug5tau3 + epstau1) to spinor l
 * and stores it in k
 *
 * it is assumed that the clover leaf is computed and stored
 * in sw[VOLUME][3][2]
 * the corresponding routine can be found in clover_leaf.c
 *
 **************************************************************/


void assign_mul_one_sw_pm_imu_eps(const int ieo, 
				  spinor * const k_s, spinor * const k_c, 
				  const spinor * const l_s, const spinor * const l_c,
				  const double mu, const double eps) {
#ifdef OMP
#pragma omp parallel
  {
#endif
  su3_vector ALIGN chi, psi1, psi2;
  int ix;
  int ioff;
  const su3 *w1, *w2, *w3;
  spinor *r_s, *r_c;
  const spinor *s_s, *s_c;
  
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
  for(unsigned int icx = ioff; icx < (VOLUME/2+ioff); icx++) {
    ix = g_eo2lexic[icx];
    
    r_s = k_s + icx-ioff;
    r_c = k_c + icx-ioff;
    s_s = l_s + icx-ioff;
    s_c = l_c + icx-ioff;

    // upper two spin components first
    w1=&sw[ix][0][0];
    w2=w1+2; /*&sw[ix][1][0];*/
    w3=w1+4; /*&sw[ix][2][0];*/
    _su3_multiply(psi1, *w1, (*s_s).s0); 
    _su3_multiply(chi, *w2, (*s_s).s1);
    _vector_add_assign(psi1, chi);
    _su3_inverse_multiply(psi2, *w2, (*s_s).s0); 
    _su3_multiply(chi, *w3, (*s_s).s1);
    _vector_add_assign(psi2, chi); 

    // add in the twisted mass term (plus in the upper components)
    _vector_add_i_mul(psi1, mu, (*s_s).s0);
    _vector_add_i_mul(psi2, mu, (*s_s).s1);

    _vector_add_mul(psi1, eps, (*s_c).s0);
    _vector_add_mul(psi2, eps, (*s_c).s1);

    _vector_assign((*r_s).s0, psi1);
    _vector_assign((*r_s).s1, psi2);

    _su3_multiply(psi1, *w1, (*s_c).s0); 
    _su3_multiply(chi, *w2, (*s_c).s1);
    _vector_add_assign(psi1, chi);
    _su3_inverse_multiply(psi2, *w2, (*s_c).s0); 
    _su3_multiply(chi, *w3, (*s_c).s1);
    _vector_add_assign(psi2, chi); 

    // add in the twisted mass term (plus in the upper components)
    _vector_add_i_mul(psi1, -mu, (*s_c).s0);
    _vector_add_i_mul(psi2, -mu, (*s_c).s1);

    _vector_add_mul(psi1, eps, (*s_s).s0);
    _vector_add_mul(psi2, eps, (*s_s).s1);

    _vector_assign((*r_c).s0, psi1);
    _vector_assign((*r_c).s1, psi2);

    // now lower two spin components
    w1++; /*=&sw[ix][0][1];*/
    w2++; /*=&sw[ix][1][1];*/
    w3++; /*=&sw[ix][2][1];*/
    _su3_multiply(psi1, *w1, (*s_s).s2); 
    _su3_multiply(chi, *w2, (*s_s).s3);
    _vector_add_assign(psi1, chi); 
    _su3_inverse_multiply(psi2, *w2, (*s_s).s2); 
    _su3_multiply(chi, *w3, (*s_s).s3);
    _vector_add_assign(psi2, chi); 

    // add in the twisted mass term (minus from g5 in the lower components)
    _vector_add_i_mul(psi1, -mu, (*s_s).s2);
    _vector_add_i_mul(psi2, -mu, (*s_s).s3);

    _vector_add_mul(psi1, eps, (*s_c).s2);
    _vector_add_mul(psi2, eps, (*s_c).s3);

    _vector_assign((*r_s).s2, psi1);
    _vector_assign((*r_s).s3, psi2);

    _su3_multiply(psi1, *w1, (*s_c).s2); 
    _su3_multiply(chi, *w2, (*s_c).s3);
    _vector_add_assign(psi1, chi); 
    _su3_inverse_multiply(psi2, *w2, (*s_c).s2); 
    _su3_multiply(chi, *w3, (*s_c).s3);
    _vector_add_assign(psi2, chi); 

    // add in the twisted mass term (minus from g5 in the lower components)
    _vector_add_i_mul(psi1, mu, (*s_c).s2);
    _vector_add_i_mul(psi2, mu, (*s_c).s3);

    _vector_add_mul(psi1, eps, (*s_s).s2);
    _vector_add_mul(psi2, eps, (*s_s).s3);

    _vector_assign((*r_c).s2, psi1);
    _vector_assign((*r_c).s3, psi2);

  }
#ifdef OMP
  } /* OpenMP closing brace */
#endif
  return;
}



void assign_mul_one_sw_pm_imu_inv(const int ieo, 
				  spinor * const k, const spinor * const l,
				  const double mu) {
#ifdef OMP
#pragma omp parallel
  {
#endif
  su3_vector ALIGN psi, chi, phi1, phi3;
  const su3 *w1, *w2, *w3, *w4;
  const spinor *rn;
  spinor *s;

  /************************ loop over all lattice sites *************************/
#ifdef OMP
#pragma omp for
#endif
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
#ifdef OMP
  } /* OpenMP closing brace */
#endif
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
    sw[0][0] = (su3*)(((unsigned long int)(_sw)+ALIGN_BASE)&~ALIGN_BASE);
    sw_inv[0][0] = (su3*)(((unsigned long int)(_sw_inv)+ALIGN_BASE)&~ALIGN_BASE);
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

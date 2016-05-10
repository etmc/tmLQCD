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

// work-around for missing single precision implementation of inline SSE
#ifdef SSE
#define REDEFSSE
#undef SSE
#endif

#ifdef SSE2
#define REDEFSSE2
#undef SSE2
#endif

#ifdef SSE3
#define REDEFSSE3
#undef SSE3
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
#include "operator/Hopping_Matrix_32.h"

#include "tm_operators.h"
#include "tm_operators_32.h"

#include "operator/clovertm_operators.h"
#include "operator/clovertm_operators_32.h"


void Qsw_pm_psi_32(spinor32 * const l, spinor32 * const k) {
#ifdef OMP
#pragma omp parallel
  {
#endif
  /* \hat Q_{-} */
  Hopping_Matrix_32_orphaned(EO, g_spinor_field32[1], k);
  clover_inv_32_orphaned(g_spinor_field32[1], -1, g_mu);
  Hopping_Matrix_32_orphaned(OE, g_spinor_field32[0], g_spinor_field32[1]);
  clover_gamma5_32_orphaned(OO, g_spinor_field32[0], k, g_spinor_field32[0], -(g_mu + g_mu3));
  /* \hat Q_{+} */
  Hopping_Matrix_32_orphaned(EO, l, g_spinor_field32[0]);
  clover_inv_32_orphaned(l, +1, g_mu); 
  Hopping_Matrix_32_orphaned(OE, g_spinor_field32[1], l);
  clover_gamma5_32_orphaned(OO, l, g_spinor_field32[0], g_spinor_field32[1], +(g_mu + g_mu3));
#ifdef OMP
  } /* OpenMP parallel closing brace */
#endif
}

void clover_inv_32_orphaned(spinor32 * const l, const int tau3sign, const double mu) {
  int icy;
  su3_vector32 ALIGN32 psi, chi, phi1, phi3;
  int ioff = 0;
  const su3_32 *w1, *w2, *w3, *w4;
  spinor32 *rn;

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

    w1=&sw_inv_32[icy][0][0];
    w2=w1+2;  /* &sw_inv_32[icy][1][0]; */
    w3=w1+4;  /* &sw_inv_32[icy][2][0]; */
    w4=w1+6;  /* &sw_inv_32[icy][3][0]; */
    _su3_multiply(psi,*w1,phi1); 
    _su3_multiply(chi,*w2,(*rn).s1);
    _vector_add((*rn).s0,psi,chi);
    _su3_multiply(psi,*w4,phi1); 
    _su3_multiply(chi,*w3,(*rn).s1);
    _vector_add((*rn).s1,psi,chi);

    w1++; /* &sw_inv_32[icy][0][1]; */
    w2++; /* &sw_inv_32[icy][1][1]; */
    w3++; /* &sw_inv_32[icy][2][1]; */
    w4++; /* &sw_inv_32[icy][3][1]; */
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
}

void clover_inv_32(spinor32 * const l, const int tau3sign, const double mu) {
#ifdef OMP
#pragma omp parallel
  {
#endif
  clover_inv_32_orphaned(l,tau3sign,mu);
#ifdef OMP
  } /* OpenMP closing brace */
#endif
  return;
}

void clover_inv_nd_32_orphaned(const int ieo, spinor32 * const l_c, spinor32 * const l_s) {
  int icy;
  su3_vector32 ALIGN psi, chi, phi1, phi3;
  int ioff = 0;
  const su3_32 *w1, *w2, *w3, *w4;
  spinor32 *rn_s, *rn_c;


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

    w1=&sw_inv_32[icy][0][0];
    w2=w1+2;  /* &sw_inv_32[icy][1][0]; */
    w3=w1+4;  /* &sw_inv_32[icy][2][0]; */
    w4=w1+6;  /* &sw_inv_32[icy][3][0]; */
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

    w1++; /* &sw_inv_32[icy][0][1]; */
    w2++; /* &sw_inv_32[icy][1][1]; */
    w3++; /* &sw_inv_32[icy][2][1]; */
    w4++; /* &sw_inv_32[icy][3][1]; */
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
  return;
}

void clover_inv_nd_32(const int ieo, spinor32 * const l_c, spinor32 * const l_s) {
#ifdef OMP
#pragma omp parallel
  {
#endif
  clover_inv_nd_32_orphaned(ieo,l_c,l_s);
#ifdef OMP
  } /* OpenMP closing brace */
#endif
  return;
}

void clover_gamma5_32_orphaned(const int ieo, 
		   spinor32 * const l, const spinor32 * const k, const spinor32 * const j,
		   const double mu) {

  su3_vector32 ALIGN32 chi, psi1, psi2;
  int ix;
  int ioff,icx;
  const su3_32 *w1,*w2,*w3;
  spinor32 *r;
  const spinor32 *s,*t;

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
    
    w1=&sw_32[ix][0][0];
    w2=w1+2; /*&sw[ix][1][0];*/
    w3=w1+4; /*&sw[ix][2][0];*/
    _su3_multiply(psi1,*w1,(*s).s0); 
    _su3_multiply(chi,*w2,(*s).s1);
    _vector_add_assign(psi1,chi);
    _su3_inverse_multiply(psi2,*w2,(*s).s0); 
    _su3_multiply(chi,*w3,(*s).s1);
    _vector_add_assign(psi2,chi); 
    // add in the twisted mass term (plus in the upper components)
    _vector_add_i_mul(psi1, (float)mu, (*s).s0);
    _vector_add_i_mul(psi2, (float)mu, (*s).s1);

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
}

void clover_gamma5_32(const int ieo, 
		   spinor32 * const l, const spinor32 * const k, const spinor32 * const j,
		   const double mu) {
#ifdef OMP
#pragma omp parallel
  {
#endif
  clover_gamma5_32_orphaned(ieo,l,k,j,mu);
#ifdef OMP
  } /* OMP closing brace */
#endif
  return;
}

void clover_gamma5_nd_32_orphaned(const int ieo, 
          spinor32 * const l_c, spinor32 * const l_s, 
          const spinor32 * const k_c, const spinor32 * const k_s, 
          const spinor32 * const j_c, const spinor32 * const j_s,
          const float mubar, const float epsbar) {
  su3_vector32 ALIGN chi, psi1, psi2;
  int ix;
  int ioff;
  const su3_32 *w1,*w2,*w3;
  spinor32 *r_s, *r_c;
  const spinor32 *s_s, *s_c, *t_s, *t_c;
  
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
    w1=&sw_32[ix][0][0];
    w2=w1+2; /*&sw_32[ix][1][0];*/
    w3=w1+4; /*&sw_32[ix][2][0];*/
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
    w1++; /*=&sw_32[ix][0][1];*/
    w2++; /*=&sw_32[ix][1][1];*/
    w3++; /*=&sw_32[ix][2][1];*/
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
}

void clover_gamma5_nd_32(const int ieo, 
          spinor32 * const l_c, spinor32 * const l_s, 
          const spinor32 * const k_c, const spinor32 * const k_s, 
          const spinor32 * const j_c, const spinor32 * const j_s,
          const float mubar, const float epsbar) {
#ifdef OMP
#pragma omp parallel
  {
#endif
  clover_gamma5_nd_32_orphaned(ieo,l_c,l_s,k_c,k_s,j_c,j_s,mubar,epsbar);
#ifdef OMP
  } /* OpenMP parallel closing brace */
#endif
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

void assign_mul_one_sw_pm_imu_eps_32_orphaned(const int ieo, 
          spinor32 * const k_s, spinor32 * const k_c, 
          const spinor32 * const l_s, const spinor32 * const l_c,
          const float mu, const float eps) {
  su3_vector32 ALIGN chi, psi1, psi2;
  int ix;
  int ioff;
  const su3_32 *w1, *w2, *w3;
  spinor32 *r_s, *r_c;
  const spinor32 *s_s, *s_c;
  
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
    w1=&sw_32[ix][0][0];
    w2=w1+2; /*&sw_32[ix][1][0];*/
    w3=w1+4; /*&sw_32[ix][2][0];*/
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
    w1++; /*=&sw_32[ix][0][1];*/
    w2++; /*=&sw_32[ix][1][1];*/
    w3++; /*=&sw_32[ix][2][1];*/
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
}

void assign_mul_one_sw_pm_imu_eps_32(const int ieo, 
          spinor32 * const k_s, spinor32 * const k_c, 
          const spinor32 * const l_s, const spinor32 * const l_c,
          const float mu, const float eps) {
  #ifdef OMP
  #pragma omp parallel
  {
  #endif
  assign_mul_one_sw_pm_imu_eps_32_orphaned(ieo,k_s,k_c,l_s,l_c,mu,eps);
  #ifdef OMP
  } /* OpenMP parallel closing brace */
  #endif
}

#ifdef REDEFSSE
#undef REDEFSSE
#define SSE
#endif

#ifdef REDEFSSE2
#undef REDEFSSE2
#define SSE2
#endif

#ifdef REDEFSSE3
#undef REDEFSSE3
#define SSE3
#endif

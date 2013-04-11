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
#include "operator/swinv_times_Hopping_Matrix.h"
#include "operator/tm_operators.h"
#include "operator/clover_inline.h"
#include "operator/clovertm_operators.h"


_Complex double *** sw;
_Complex double *** sw_inv;

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
  //Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX+1], k);
  //clover_inv(EE, g_spinor_field[DUM_MATRIX+1], -g_mu);
  swinv_times_Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX+1], k, -g_mu);
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
#endif
  int icy;
  int ioff = 0;
  _Complex double * restrict w, * restrict r;

#if (defined BGQ && defined XLC)
  vector4double r0, r1, r2, r3, r4, r5, r6;
  vector4double w0, w1, w2;
#else 
  _Complex double s[6];
#endif

  if(mu < 0) ioff = VOLUME/2;

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

    r = (_Complex double *) (l + icx);
    w = sw_inv[icy][0];
    _mul_colourmatrix(w, r);

    r += 6;
    w = sw_inv[icy][1];
    _mul_colourmatrix(w, r);

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
  int ioff = 0;
#if (defined BGQ && defined XLC)
  vector4double r0, r1, r2, r3, r4, r5;
  vector4double r6, r7, r8, r9, r10, r11;
  vector4double w0, w1, w2;
#else 
  _Complex double ALIGN s[12];
#endif
  _Complex double * restrict w, * restrict rn_s, * restrict rn_c;


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

    rn_s = (_Complex double *) (l_s + icx);
    rn_c = (_Complex double *) (l_c + icx);
    w = sw_inv[icy][0];
    _mul_colourmatrix_nd(w, rn_s, rn_c);

    rn_s += 6;
    rn_c += 6;
    w = sw_inv[icy][1];
    _mul_colourmatrix_nd(w, rn_s, rn_c);

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
 * in sw[VOLUME][2][36]
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
  int ix;
  int ioff;
#if (defined BGQ && defined XLC)
  vector4double r0, r1, r2, r3, r4, r5;
  vector4double r6, r7, r8, tmp;
  vector4double w0, w1, w2;
#endif
  _Complex double * restrict w;
  _Complex double * restrict r, * restrict s, * restrict t;

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
  for(int icx = ioff; icx < (VOLUME/2+ioff); icx++) {
    ix = g_eo2lexic[icx];
    
    r = (_Complex double*)(l + icx-ioff);
    s = (_Complex double*)(k + icx-ioff);
    t = (_Complex double*)(j + icx-ioff);

    w = sw[ix][0];
    _mul_colourmatrix_add_i_mul_sub_assign(r, w, mu, s, t);

    r += 6;
    s += 6;
    t += 6;

    // includes g5 multiplication
    w = sw[ix][1];
    _mul_colourmatrix_add_i_mul_sub_nassign(r, w, -mu, s, t);

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
  int ix;
  int ioff;
#if (defined BGQ && defined XLC)
  vector4double r0, r1, r2, r3, r4, r5;
  vector4double r6, r7, r8, tmp;
  vector4double w0, w1, w2;
#endif
  _Complex double * restrict w;
  _Complex double * restrict r, * restrict s, * restrict t;

  
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
    
    
    r = (_Complex double*)(l + icx-ioff);
    s = (_Complex double*)(k + icx-ioff);
    t = (_Complex double*)(j + icx-ioff);

    w = sw[ix][0];
    _mul_colourmatrix_add_i_mul_sub_assign(r, w, mu, s, t);

    r += 6;
    s += 6;
    t += 6;

    w = sw[ix][1];
    _mul_colourmatrix_add_i_mul_sub_assign(r, w, -mu, s, t);

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
  int ix;
  int ioff;
#if (defined BGQ && defined XLC)
  vector4double r0, r1, r2, r3, r4, r5;
  vector4double r6, r7, r8, r9, r10, r11;
  vector4double w0, w1, w2;
#endif
  _Complex double * restrict w;
  _Complex double  *r_s, *r_c, *s_s, *s_c, *t_s, *t_c;
  
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
    
    r_s = (_Complex double*)(l_s + icx-ioff);
    r_c = (_Complex double*)(l_c + icx-ioff);
    s_s = (_Complex double*)(k_s + icx-ioff);
    s_c = (_Complex double*)(k_c + icx-ioff);
    t_s = (_Complex double*)(j_s + icx-ioff);
    t_c = (_Complex double*)(j_c + icx-ioff);

    w = sw[ix][0];
    _mul_colourmatrix_assign_nd(r_s, r_c, w, s_s, s_c);
    _add_nd_massterm_sub(r_s, r_c, mubar, epsbar, s_s, s_c, t_s, t_c);

    r_s += 6; r_c += 6;
    s_s += 6; s_c += 6;
    t_s += 6; t_c += 6;

    w = sw[ix][1];
    _mul_colourmatrix_assign_nd(r_s, r_c, w, s_s, s_c);
    _add_nd_massterm_sub(r_s, r_c, -mubar, epsbar, s_s, s_c, t_s, t_c);

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
  int ix;
  int ioff;
#if (defined BGQ && defined XLC)
  vector4double r0, r1, r2, r3, r4, r5;
  vector4double r6, r7, r8, r9, r10, r11;
  vector4double w0, w1, w2;
#endif
  _Complex double * restrict w;
  _Complex double  *r_s, *r_c, *s_s, *s_c, *t_s, *t_c;
  
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
    
    r_s = (_Complex double*)(l_s + icx-ioff);
    r_c = (_Complex double*)(l_c + icx-ioff);
    s_s = (_Complex double*)(k_s + icx-ioff);
    s_c = (_Complex double*)(k_c + icx-ioff);
    t_s = (_Complex double*)(j_s + icx-ioff);
    t_c = (_Complex double*)(j_c + icx-ioff);

    w = sw[ix][0];
    _mul_colourmatrix_assign_nd(r_s, r_c, w, s_s, s_c);
    _add_nd_massterm_sub(r_s, r_c, mubar, epsbar, s_s, s_c, t_s, t_c);

    r_s += 6; r_c += 6;
    s_s += 6; s_c += 6;
    t_s += 6; t_c += 6;

    w = sw[ix][1];
    _mul_colourmatrix_assign_nd(r_s, r_c, w, s_s, s_c);
    _add_nd_massterm_nsub(r_s, r_c, -mubar, epsbar, s_s, s_c, t_s, t_c);

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
  int ix;
  int ioff;
  _Complex double * restrict r, * restrict s, * restrict w; 
  
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
    
    r = (_Complex double*)(k + icx-ioff);
    s = (_Complex double*)(l + icx-ioff);

    w = sw[ix][0];
    _mul_colourmatrix_add_i_mul_assign(r, w, mu, s);

    r += 6;
    s += 6;
    w = sw[ix][1];
    _mul_colourmatrix_add_i_mul_assign(r, w, -mu, s);

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
  int ix;
  int ioff;
#if (defined BGQ && defined XLC)
  vector4double r0, r1, r2, r3, r4, r5;
  vector4double r6, r7, r8, r9, r10, r11;
  vector4double w0, w1, w2;
#endif
  _Complex double * restrict r_s, * restrict r_c, * restrict s_s, * restrict s_c, * restrict w;
  
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
    
    r_s = (_Complex double*)(k_s + icx-ioff);
    r_c = (_Complex double*)(k_c + icx-ioff);
    s_s = (_Complex double*)(l_s + icx-ioff);
    s_c = (_Complex double*)(l_c + icx-ioff);

    w = sw[ix][0];
    _mul_colourmatrix_assign_nd(r_s, r_c, w, s_s, s_c);
    _add_nd_massterm(r_s, r_c, mu, eps, s_s, s_c);

    r_s += 6; r_c += 6;
    s_s += 6; s_c += 6;
    w = sw[ix][1];

    _mul_colourmatrix_assign_nd(r_s, r_c, w, s_s, s_c);
    _add_nd_massterm(r_s, r_c, -mu, eps, s_s, s_c);

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
  _Complex double * restrict w, * restrict r, * restrict s;
#if (defined BGQ && defined XLC)
  vector4double r0, r1, r2, r3, r4, r5, r6;
  vector4double w0, w1, w2;
#endif

#ifdef OMP
#pragma omp for
#endif
  for(int icx = 0; icx < (VOLUME/2); icx++) {

    r = (_Complex double *) (l + icx);
    s = (_Complex double *) (k + icx);
    w = sw_inv[icx][0];
    _mul_colourmatrix_assign(s, w, r);

    r += 6;
    s += 6;
    w = sw_inv[icx][1];
    _mul_colourmatrix_assign(s, w, r);
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

_Complex double ** sw1, * _sw;
_Complex double **  sw_inv1, *_sw_inv;

void init_sw_fields() {
  int V = VOLUME;
  _Complex double * ctmp;
  static int sw_init = 0;

  if(!sw_init) {
    if((void*)(sw = (_Complex double***)calloc(V, sizeof(_Complex double**))) == NULL) {
      fprintf (stderr, "sw malloc err\n"); 
    }
    if((void*)(sw_inv = (_Complex double***)calloc(V, sizeof(_Complex double**))) == NULL) {
      fprintf (stderr, "sw_inv malloc err\n"); 
    }
    if((void*)(sw1 = (_Complex double**)calloc(2*V, sizeof(_Complex double*))) == NULL) {
      fprintf (stderr, "sw1 malloc err\n"); 
    }
    if((void*)(sw_inv1 = (_Complex double**)calloc(2*V, sizeof(_Complex double*))) == NULL) {
      fprintf (stderr, "sw_inv1 malloc err\n"); 
    }
    if((void*)(_sw = (_Complex double*)calloc(36*2*V+1, sizeof(_Complex double))) == NULL) {
      fprintf (stderr, "_sw malloc err\n"); 
    }
    if((void*)(_sw_inv = (_Complex double*)calloc(36*2*V+1, sizeof(_Complex double))) == NULL) {
      fprintf (stderr, "_sw_inv malloc err\n"); 
    }
    sw[0] = sw1;
    sw_inv[0] = sw_inv1;
    for(int i = 1; i < V; i++) {
      sw[i] = sw[i-1]+2;
      sw_inv[i] = sw_inv[i-1]+2;
    }
    sw[0][0] = (_Complex double*)(((unsigned long int)(_sw)+ALIGN_BASE)&~ALIGN_BASE);
    sw_inv[0][0] = (_Complex double*)(((unsigned long int)(_sw_inv)+ALIGN_BASE)&~ALIGN_BASE);
    ctmp = sw[0][0];
    for(int i = 0; i < V; i++) {
      for(int j = 0; j < 2; j++) {
	sw[i][j] = ctmp;
	ctmp = ctmp+36;
      }
    }
    
    ctmp = sw_inv[0][0];
    for(int i = 0; i < V; i++) {
      for(int j = 0; j < 2; j++) {
	sw_inv[i][j] = ctmp;
	ctmp = ctmp + 36;
      }
    }
    sw_init = 1;
  }
  return;
}

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
 *                                                            *
 * This file contains operators for twisted mass Wilson QCD   *
 * prepared for even odd preconditioning                      *
 *                                                            *
 * see documentation for details                              *
 * Author: Carsten Urbach                                     *
 *         urbach@physik.fu-berlin.de                         *
 **************************************************************/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include "global.h"
#include "su3.h"
#include "operator/Hopping_Matrix.h"
#include "operator/Hopping_Matrix_nocom.h"
#include "operator/tm_times_Hopping_Matrix.h"
#include "operator/tm_sub_Hopping_Matrix.h"
#include "sse.h"
#include "linalg_eo.h"
#include "gamma.h"
#include "operator/D_psi.h"
#ifdef BGL
#  include "bgl.h"
#endif
#ifdef BGQ
#  include "bgq.h"
#endif

#include "solver/dirac_operator_eigenvectors.h"

#include "tm_operators.h"

#if (defined SSE2 || defined SSE3 || defined BGL)
const int predist=2;
#endif
/* internal */

/******************************************
 * mul_one_pm_imu_inv computes
 * l = (1\pm i\mu\gamma_5)^{-1} * l
 *
 * sign is the sign used in 
 *      1\pm i\mu\gamma_5
 * l is number of input and output field
 *
 ******************************************/
void mul_one_pm_imu_inv(spinor * const l, const double _sign, const int N);
void mul_one_pm_imu(spinor * const l, const double _sign);
/******************************************
 * mul_one_pm_imu_sub_mul_gamma5 computes
 * l = gamma_5*((1\pm i\mu\gamma_5)*k - j)
 *
 * l is the number of the output field
 * k and j the numbers of the input fields
 *
 * sign indicates which sign should be used
 * in 1\pm i\mu\gamma_5
 ******************************************/
void mul_one_pm_imu_sub_mul_gamma5(spinor * const l, spinor * const k, 
				   spinor * const j, const double _sign);
void mul_one_sub_mul_gamma5(spinor * const l, spinor * const k, 
			    spinor * const j);

/******************************************
 * mul_one_pm_imu_sub_mul computes
 * l = ((1\pm i\mu\gamma_5)*k - j)
 *
 * l is the number of the output field
 * k and j the numbers of the input fields
 *
 * sign indicates which sign should be used
 * in 1\pm i\mu\gamma_5
 ******************************************/
void mul_one_pm_imu_sub_mul(spinor * const l, spinor * const k,
			    spinor * const j, const double _sign, const int N);
void tm_sub_H_eo_gamma5(spinor* const l, spinor * const p, spinor * const k,
			const int ieo, const double _sign);

/* external functions */

void M_full(spinor * const Even_new, spinor * const Odd_new, 
	    spinor * const Even, spinor * const Odd) {
  /* Even sites */
  Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX], Odd);
  assign_mul_one_pm_imu(Even_new, Even, 1., VOLUME/2); 
  assign_add_mul_r(Even_new, g_spinor_field[DUM_MATRIX], -1., VOLUME/2);

  /* Odd sites */
  Hopping_Matrix(OE, g_spinor_field[DUM_MATRIX], Even);
  assign_mul_one_pm_imu(Odd_new, Odd, 1., VOLUME/2); 
  assign_add_mul_r(Odd_new, g_spinor_field[DUM_MATRIX], -1., VOLUME/2);
}

void Q_full(spinor * const Even_new, spinor * const Odd_new, 
	    spinor * const Even, spinor * const Odd) {
  /* Even sites */
  Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX], Odd);
  assign_mul_one_pm_imu(Even_new, Even, 1., VOLUME/2); 
  assign_add_mul_r(Even_new, g_spinor_field[DUM_MATRIX], -1., VOLUME/2);
  gamma5(Even_new, Even_new, VOLUME/2);

  /* Odd sites */
  Hopping_Matrix(OE, g_spinor_field[DUM_MATRIX], Even);
  assign_mul_one_pm_imu(Odd_new, Odd, 1., VOLUME/2); 
  assign_add_mul_r(Odd_new, g_spinor_field[DUM_MATRIX], -1., VOLUME/2);
  gamma5(Odd_new, Odd_new, VOLUME/2);
}

void M_minus_1_timesC(spinor * const Even_new, spinor * const Odd_new, 
		      spinor * const Even, spinor * const Odd) {
  /* Even sites */
  Hopping_Matrix(EO, Even_new, Odd);
  mul_one_pm_imu_inv(Even_new, 1., VOLUME/2); 

  /* Odd sites */
  Hopping_Matrix(OE, Odd_new, Even);
  mul_one_pm_imu_inv(Odd_new, 1., VOLUME/2); 
}



/******************************************
 *
 * This is the implementation of
 *
 * \hat Q_{+} =
 * \gamma_5(M_{oo}^+ - M_{oe}(M_{ee}^+ )^{-1}M_{eo})
 *
 * see documentation for details
 * k is the number of the input field
 * l is the number of the output field
 *
 * it acts only on the odd part or only 
 * on a half spinor
 ******************************************/
void Qtm_plus_psi(spinor * const l, spinor * const k){
  Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX+1], k);
  mul_one_pm_imu_inv(g_spinor_field[DUM_MATRIX+1], +1., VOLUME/2);
  Hopping_Matrix(OE, g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+1]);
  mul_one_pm_imu_sub_mul_gamma5(l, k, g_spinor_field[DUM_MATRIX], +1.);
}

void Qtm_plus_psi_nocom(spinor * const l, spinor * const k){
  Hopping_Matrix_nocom(EO, g_spinor_field[DUM_MATRIX+1], k);
  mul_one_pm_imu_inv(g_spinor_field[DUM_MATRIX+1], +1., VOLUME/2);
  Hopping_Matrix_nocom(OE, g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+1]);
  mul_one_pm_imu_sub_mul_gamma5(l, k, g_spinor_field[DUM_MATRIX], +1.);
}

void Qtm_plus_sym_psi(spinor * const l, spinor * const k){
  Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX+1], k);
  mul_one_pm_imu_inv(g_spinor_field[DUM_MATRIX+1], +1., VOLUME/2);
  Hopping_Matrix(OE, g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+1]);
  mul_one_pm_imu_inv(g_spinor_field[DUM_MATRIX], +1., VOLUME/2);
  mul_one_sub_mul_gamma5(l, k, g_spinor_field[DUM_MATRIX]);
}

void Qtm_plus_sym_psi_nocom(spinor * const l, spinor * const k){
  Hopping_Matrix_nocom(EO, g_spinor_field[DUM_MATRIX+1], k);
  mul_one_pm_imu_inv(g_spinor_field[DUM_MATRIX+1], +1., VOLUME/2);
  Hopping_Matrix_nocom(OE, g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+1]);
  mul_one_pm_imu_inv(g_spinor_field[DUM_MATRIX], +1., VOLUME/2);
  mul_one_sub_mul_gamma5(l, k, g_spinor_field[DUM_MATRIX]);
}

/******************************************
 *
 * This is the implementation of
 *
 * \hat Q_{-} =
 * \gamma_5(M_{oo}^- - M_{oe}(M_{ee}^- )^{-1}M_{eo})
 *
 * see documentation for details
 * k is the number of the input field
 * l is the number of the output field
 *
 * it acts only on the odd part or only 
 * on a half spinor
 ******************************************/
void Qtm_minus_psi(spinor * const l, spinor * const k) {
  H_eo_tm_inv_psi(g_spinor_field[DUM_MATRIX+1], k, EO, -1);
  Hopping_Matrix(OE, g_spinor_field[DUM_MATRIX+2], g_spinor_field[DUM_MATRIX+1]);
  mul_one_pm_imu_sub_mul_gamma5(l, k, g_spinor_field[DUM_MATRIX+2], -1);
  //tm_sub_H_eo_gamma5(l, k, g_spinor_field[DUM_MATRIX+1], OE, -1.);
}

void Qtm_minus_sym_psi(spinor * const l, spinor * const k){
  Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX+1], k);
  mul_one_pm_imu_inv(g_spinor_field[DUM_MATRIX+1], -1., VOLUME/2);
  Hopping_Matrix(OE, g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+1]);
  mul_one_pm_imu_inv(g_spinor_field[DUM_MATRIX], -1., VOLUME/2);
  mul_one_sub_mul_gamma5(l, k, g_spinor_field[DUM_MATRIX]);
}

/******************************************
 *
 * This is the implementation of
 *
 * \gamma_5 \hat Q_{+} =
 * (M_{oo}^+ - M_{oe}(M_{ee}^+ )^{-1}M_{eo})
 *
 * see documentation for details
 * k is the number of the input field
 * l is the number of the output field
 *
 * it acts only on the odd part or only 
 * on a half spinor
 ******************************************/
void Mtm_plus_psi(spinor * const l, spinor * const k){
  Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX+1], k);
  mul_one_pm_imu_inv(g_spinor_field[DUM_MATRIX+1], +1., VOLUME/2);
  Hopping_Matrix(OE, g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+1]);
  mul_one_pm_imu_sub_mul(l, k, g_spinor_field[DUM_MATRIX], +1., VOLUME/2);
}

void Mtm_plus_psi_nocom(spinor * const l, spinor * const k){
  Hopping_Matrix_nocom(EO, g_spinor_field[DUM_MATRIX+1], k);
  mul_one_pm_imu_inv(g_spinor_field[DUM_MATRIX+1], +1., VOLUME/2);
  Hopping_Matrix_nocom(OE, g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+1]);
  mul_one_pm_imu_sub_mul(l, k, g_spinor_field[DUM_MATRIX], +1., VOLUME/2);
}

void Mtm_plus_sym_psi(spinor * const l, spinor * const k){
  Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX+1], k);
  mul_one_pm_imu_inv(g_spinor_field[DUM_MATRIX+1], +1., VOLUME/2);
  Hopping_Matrix(OE, g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+1]);
  mul_one_pm_imu_inv(g_spinor_field[DUM_MATRIX], +1., VOLUME/2);
  diff(l, k, g_spinor_field[DUM_MATRIX], VOLUME/2);
}

void Mtm_plus_sym_psi_nocom(spinor * const l, spinor * const k){
  Hopping_Matrix_nocom(EO, g_spinor_field[DUM_MATRIX+1], k);
  mul_one_pm_imu_inv(g_spinor_field[DUM_MATRIX+1], +1., VOLUME/2);
  Hopping_Matrix_nocom(OE, g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+1]);
  mul_one_pm_imu_inv(g_spinor_field[DUM_MATRIX], +1., VOLUME/2);
  diff(l, k, g_spinor_field[DUM_MATRIX], VOLUME/2);
}

/******************************************
 *
 * This is the implementation of
 *
 * \gamma_5 \hat Q_{-} =
 * (M_{oo}^- - M_{oe}(M_{ee}^- )^{-1}M_{eo})
 *
 * see documentation for details
 * k is the number of the input field
 * l is the number of the output field
 *
 * it acts only on the odd part or only 
 * on a half spinor
 ******************************************/
void Mtm_minus_psi(spinor * const l, spinor * const k) {
  Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX+1], k);
  mul_one_pm_imu_inv(g_spinor_field[DUM_MATRIX+1], -1., VOLUME/2);
  Hopping_Matrix(OE, g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+1]);
  mul_one_pm_imu_sub_mul(l, k, g_spinor_field[DUM_MATRIX], -1., VOLUME/2);
}

void Mtm_minus_sym_psi(spinor * const l, spinor * const k) {
  Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX+1], k);
  mul_one_pm_imu_inv(g_spinor_field[DUM_MATRIX+1], -1., VOLUME/2);
  Hopping_Matrix(OE, g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+1]);
  mul_one_pm_imu_inv(g_spinor_field[DUM_MATRIX], -1., VOLUME/2);
  diff(l, k, g_spinor_field[DUM_MATRIX], VOLUME/2);
}

void Mtm_minus_sym_psi_nocom(spinor * const l, spinor * const k) {
  Hopping_Matrix_nocom(EO, g_spinor_field[DUM_MATRIX+1], k);
  mul_one_pm_imu_inv(g_spinor_field[DUM_MATRIX+1], -1., VOLUME/2);
  Hopping_Matrix_nocom(OE, g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+1]);
  mul_one_pm_imu_inv(g_spinor_field[DUM_MATRIX], -1., VOLUME/2);
  diff(l, k, g_spinor_field[DUM_MATRIX], VOLUME/2);
}

/******************************************
 *
 * This is the implementation of
 *
 * \hat Q_{+} \hat Q_{-} 
 *
 * see documentation for details
 * k is the number of the input field
 * l is the number of the output field
 *
 * it acts only on the odd part or only 
 * on a half spinor
 ******************************************/
void Qtm_pm_psi(spinor * const l, spinor * const k){
  /* Q_{-} */
  H_eo_tm_inv_psi(g_spinor_field[DUM_MATRIX+1], k, EO, -1);
  tm_sub_H_eo_gamma5(g_spinor_field[DUM_MATRIX], k, g_spinor_field[DUM_MATRIX+1], OE, -1);
  /* Q_{+} */
  H_eo_tm_inv_psi(g_spinor_field[DUM_MATRIX+1], g_spinor_field[DUM_MATRIX], EO, +1);
  tm_sub_H_eo_gamma5(l, g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+1], OE, +1);
}

void Qtm_pm_sym_psi(spinor * const l, spinor * const k){
  /* Q_{-} */
  Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX+1], k);
  mul_one_pm_imu_inv(g_spinor_field[DUM_MATRIX+1], -1., VOLUME/2);
  Hopping_Matrix(OE, g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+1]);
  mul_one_pm_imu_inv(g_spinor_field[DUM_MATRIX], -1., VOLUME/2);
  diff(l, k, g_spinor_field[DUM_MATRIX], VOLUME/2);
  gamma5(l, l, VOLUME/2);

  /* Q_{+} */
  Hopping_Matrix(EO, l, g_spinor_field[DUM_MATRIX]);
  mul_one_pm_imu_inv(l, +1., VOLUME/2);
  Hopping_Matrix(OE, g_spinor_field[DUM_MATRIX+1], l);
  mul_one_pm_imu_inv(g_spinor_field[DUM_MATRIX], +1., VOLUME/2);
  diff(l, k, g_spinor_field[DUM_MATRIX], VOLUME/2);
  gamma5(l, l, VOLUME/2);

}

void Qtm_pm_psi_nocom(spinor * const l, spinor * const k){
  /* Q_{-} */
  Hopping_Matrix_nocom(EO, g_spinor_field[DUM_MATRIX+1], k);
  mul_one_pm_imu_inv(g_spinor_field[DUM_MATRIX+1], -1., VOLUME/2);
  Hopping_Matrix_nocom(OE, g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+1]);
  mul_one_pm_imu_sub_mul_gamma5(g_spinor_field[DUM_MATRIX], k, g_spinor_field[DUM_MATRIX], -1.);
  /* Q_{+} */
  Hopping_Matrix_nocom(EO, l, g_spinor_field[DUM_MATRIX]);
  mul_one_pm_imu_inv(l, +1., VOLUME/2);
  Hopping_Matrix_nocom(OE, g_spinor_field[DUM_MATRIX+1], l);
  mul_one_pm_imu_sub_mul_gamma5(l, g_spinor_field[DUM_MATRIX], g_spinor_field[DUM_MATRIX+1], +1.);
}

/* the "full" operators */
void Q_pm_psi(spinor * const l, spinor * const k)
{
  g_mu = -g_mu;
  D_psi(l, k);
  gamma5(g_spinor_field[DUM_MATRIX], l, VOLUME);
  g_mu = -g_mu;
  D_psi(l, g_spinor_field[DUM_MATRIX]);
  gamma5(l, l, VOLUME);
}


/* the "full" operators */
void Q_pm_psi_prec(spinor * const l, spinor * const k)
{
  spinorPrecWS *ws=(spinorPrecWS*)g_precWS;

  _Complex double ALIGN alpha= -1.0;

  if(g_prec_sequence_d_dagger_d[0]!=0.0)
  {
    alpha = g_prec_sequence_d_dagger_d[0];
    spinorPrecondition(l,k,ws,T,L,alpha,0,1);
  } 
  else
    assign(l,k,VOLUME);

  g_mu = -g_mu;
  D_psi(g_spinor_field[DUM_MATRIX], l);
  gamma5(l, g_spinor_field[DUM_MATRIX], VOLUME);
  g_mu = -g_mu;

  if(g_prec_sequence_d_dagger_d[1]!=0.0)
  {
    alpha = g_prec_sequence_d_dagger_d[1];
    spinorPrecondition(l,l,ws,T,L,alpha,0,1);
  }

  D_psi(g_spinor_field[DUM_MATRIX], l);
  gamma5(l, g_spinor_field[DUM_MATRIX], VOLUME);

  if(g_prec_sequence_d_dagger_d[2]!=0.0)
  {
    alpha = g_prec_sequence_d_dagger_d[2]; 
    spinorPrecondition(l,l,ws,T,L,alpha,0,1);
  }

}



/* This is the version for the gpu with interchanged order of gamma5 and D_psi (Florian Burger)*/
void Q_pm_psi_gpu(spinor * const l, spinor * const k)
{
  gamma5(k, k, VOLUME);
  g_mu = -g_mu;
  D_psi(l, k);
  gamma5(g_spinor_field[DUM_MATRIX], l, VOLUME);
  g_mu = -g_mu;
  D_psi(l, g_spinor_field[DUM_MATRIX]);
  
}

/* the "full" operators */
void Q_pm_psi2(spinor * const l, spinor * const k)
{
  g_mu = -10.*g_mu;
  D_psi(l, k);
  gamma5(g_spinor_field[DUM_MATRIX], l, VOLUME);
  g_mu = -g_mu/10.;
  D_psi(l, g_spinor_field[DUM_MATRIX]);
  gamma5(l, l, VOLUME);
}

void Q_minus_psi(spinor * const l, spinor * const k)
{
  g_mu = -g_mu;
  D_psi(l, k);
  g_mu = -g_mu;
  gamma5(l, l, VOLUME);
}

/* This is the version for the gpu (Florian Burger)*/
void Q_minus_psi_gpu(spinor * const l, spinor * const k)
{
  gamma5(k,k,VOLUME);
  g_mu = -g_mu;
  D_psi(l, k);
  g_mu = -g_mu;
  gamma5(l, l, VOLUME);
}

void Q_plus_psi(spinor * const l, spinor * const k)
{
  D_psi(l, k);
  gamma5(l, l, VOLUME);
}

/******************************************
 *
 * This is the implementation of
 *
 * (M_{ee}^\pm)^{-1}M_{eo}
 *
 * see documentation for details
 * k is the number of the input field
 * l is the number of the output field
 *
 * it acts only on the odd part or only 
 * on a half spinor
 ******************************************/

void H_eo_tm_inv_psi(spinor * const l, spinor * const k, 
		     const int ieo, const double _sign) {
#if ((defined BGL && defined XLC) || defined _USE_TSPLITPAR)
  Hopping_Matrix(ieo, l, k);
  mul_one_pm_imu_inv(l, _sign, VOLUME/2);
#else
  double ALIGN nrm = 1./(1.+g_mu*g_mu);
  double sign=-1.; 
  complex double ALIGN z;
  if(_sign < 0.){
    sign = 1.; 
  }

  z = nrm + (sign * nrm * g_mu) * I;
  tm_times_Hopping_Matrix(ieo, l, k, z);
  return;
#endif

}

void tm_sub_H_eo_gamma5(spinor* const l, spinor * const p, spinor * const k,
                       const int ieo, const double _sign) {
#if ((defined BGL && defined XLC) || defined _USE_TSPLITPAR)
  Hopping_Matrix(ieo, g_spinor_field[DUM_MATRIX+2], k);
  mul_one_pm_imu_sub_mul_gamma5(l, p, g_spinor_field[DUM_MATRIX+2], _sign);
#else
  _Complex double ALIGN z;
  double sign=1.;

  if(_sign < 0.){
    sign = -1.;
  }

  z = 1. + (sign * g_mu) * I;
  tm_sub_Hopping_Matrix(ieo, l, p, k, z);
#endif

  return;
}


/**********************************************
 *
 * All the results are only stored in the first
 * half of the spinor fields, they have only
 * length VOLUME/2
 * 
 * That's why mul_... do not need a iput
 * parameter ieo.
 *
 * the next functions are internal and you 
 * can find comments above at the declaration 
 *
 **********************************************/

void mul_one_pm_imu_inv(spinor * const l, const double _sign, const int N){
#ifdef OMP
#pragma omp parallel
  {
#endif
  _Complex double ALIGN z,w;
  int ix;
  double sign=-1.; 
  spinor *r;

  su3_vector ALIGN phi1;

  double ALIGN nrm = 1./(1.+g_mu*g_mu);

  if(_sign < 0.){
    sign = 1.; 
  }

  z = nrm + (sign * nrm * g_mu) * I;
  w = conj(z);
  /************ loop over all lattice sites ************/
#ifdef OMP
#pragma omp for
#endif
  for(ix = 0; ix < N; ix++){
    r=l + ix;
    /* Multiply the spinorfield with the inverse of 1+imu\gamma_5 */
#if ( defined SSE2 || defined SSE3 )
    _prefetch_spinor((r+predist)); 
    _sse_load_up(r->s0);
    _sse_vector_cmplx_mul(z);
    _sse_store_nt_up(r->s0);
    _sse_load_up(r->s1);
    _sse_vector_cmplx_mul_two();
    _sse_store_nt_up(r->s1);
    _sse_load_up(r->s2);
    _sse_vector_cmplx_mul(w);
    _sse_store_nt_up(r->s2);
    _sse_load_up(r->s3);
    _sse_vector_cmplx_mul_two();
    _sse_store_nt_up(r->s3);
#else
    _complex_times_vector(phi1, z, r->s0);
    _vector_assign(r->s0, phi1);
    _complex_times_vector(phi1, z, r->s1);
    _vector_assign(r->s1, phi1);
    _complex_times_vector(phi1, w, r->s2);
    _vector_assign(r->s2, phi1);
    _complex_times_vector(phi1, w, r->s3);
    _vector_assign(r->s3, phi1);
#endif
  }

#ifdef OMP
  } /* OpenMP closing brace */
#endif

}

void assign_mul_one_pm_imu_inv(spinor * const l, spinor * const k, const double _sign, const int N){
#ifdef OMP
#pragma omp parallel
  {
#endif
  _Complex double z,w;
  int ix;
  double sign=-1.; 
  spinor *r, *s;
  double nrm = 1./(1.+g_mu*g_mu);

  if(_sign < 0.){
    sign = 1.; 
  }

  z = nrm + (sign * nrm * g_mu) * I;
  w = conj(z);

  /************ loop over all lattice sites ************/
#ifdef OMP
#pragma omp for
#endif
  for(ix = 0; ix < N; ix++){
    r=k+ix;
    s=l+ix;
    /* Multiply the spinorfield with the inverse of 1+imu\gamma_5 */
    _complex_times_vector(s->s0, z, r->s0);
    _complex_times_vector(s->s1, z, r->s1);
    _complex_times_vector(s->s2, w, r->s2);
    _complex_times_vector(s->s3, w, r->s3);
  }

#ifdef OMP
  } /* OpenMP closing brace */
#endif
}

void mul_one_pm_imu(spinor * const l, const double _sign){
#ifdef OMP
#pragma omp parallel
  {
#endif
  _Complex double z,w;
  int ix;
  double sign = 1.; 
  spinor *r;

  su3_vector ALIGN phi1;

  if(_sign < 0.){
    sign = -1.; 
  }

  z = 1. + (sign * g_mu) * I;
  w = conj(z);

  /************ loop over all lattice sites ************/
#ifdef OMP
#pragma omp for
#endif
  for(ix = 0; ix < (VOLUME/2); ix++){
    r=l+ix;
    /* Multiply the spinorfield with 1+imu\gamma_5 */
    _complex_times_vector(phi1, z, r->s0);
    _vector_assign(r->s0, phi1);
    _complex_times_vector(phi1, z, r->s1);
    _vector_assign(r->s1, phi1);
    _complex_times_vector(phi1, w, r->s2);
    _vector_assign(r->s2, phi1);
    _complex_times_vector(phi1, w, r->s3);
    _vector_assign(r->s3, phi1);
  }

#ifdef OMP
  } /* OpenMP closing brace */
#endif

}

void assign_mul_one_pm_imu(spinor * const l, spinor * const k, const double _sign, const int N){
#ifdef OMP
#pragma omp parallel
  {
#endif
  _Complex double z,w;
  int ix;
  double sign = 1.; 
  spinor *r, *s;

  if(_sign < 0.){
    sign = -1.; 
  }

  z = 1. + (sign * g_mu) * I;
  w = conj(z);

  /************ loop over all lattice sites ************/
#ifdef OMP
#pragma omp for
#endif
  for(ix = 0; ix < N; ix++){
    s=l+ix;
    r=k+ix;

    /* Multiply the spinorfield with of 1+imu\gamma_5 */
#if ( defined SSE2 || defined SSE3 )
    _prefetch_spinor((r+predist));
    _prefetch_spinor((s+predist));
    _sse_load_up(r->s0);
    _sse_vector_cmplx_mul(z);
    _sse_store_nt_up(s->s0);
    _sse_load_up(r->s1);
    _sse_vector_cmplx_mul_two();
    _sse_store_nt_up(s->s1);
    _sse_load_up(r->s2);
    _sse_vector_cmplx_mul(w);
    _sse_store_nt_up(s->s2);
    _sse_load_up(r->s3);
    _sse_vector_cmplx_mul_two();
    _sse_store_nt_up(s->s3);
#else
    _complex_times_vector(s->s0, z, r->s0);
    _complex_times_vector(s->s1, z, r->s1);
    _complex_times_vector(s->s2, w, r->s2);
    _complex_times_vector(s->s3, w, r->s3);
#endif
  }
#ifdef OMP
  } /* OpenMP closing brace */
#endif
}

void mul_one_sub_mul_gamma5(spinor * const l, spinor * const k, 
				   spinor * const j){
#ifdef OMP
#pragma omp parallel
  {
#endif
  spinor *r, *s, *t;

  /************ loop over all lattice sites ************/
#ifdef OMP
#pragma omp for
#endif
  for(int ix = 0; ix < (VOLUME/2); ++ix)
  {
    r = k+ix;
    s = j+ix;
    t = l+ix;
    /* Subtract s and store the result in t */
    /* multiply with  gamma5 included by    */
    /* reversed order of s and r (2&3)       */
    _vector_sub(t->s0, r->s0, s->s0);  
    _vector_sub(t->s1, r->s1, s->s1);  
    _vector_sub(t->s2, s->s2, r->s2);  
    _vector_sub(t->s3, s->s3, r->s3);  
  }

#ifdef OMP
  } /* OpenMP closing brace */
#endif
}


void mul_one_pm_imu_sub_mul_gamma5(spinor * const l, spinor * const k, 
				   spinor * const j, const double _sign){
#ifdef OMP
#pragma omp parallel
  {
#endif
  _Complex double z,w;
  int ix;
  double sign=1.;
  spinor *r, *s, *t;

  su3_vector ALIGN phi1, phi2, phi3, phi4;

  if(_sign < 0.){
    sign = -1.;
  }

  z = 1. + (sign * g_mu) * I;
  w = conj(z);
  
  /************ loop over all lattice sites ************/
#ifdef OMP
#pragma omp for
#endif
  for(ix = 0; ix < (VOLUME/2); ix++){
    r = k+ix;
    s = j+ix;
    t = l+ix;
    /* Multiply the spinorfield with 1+imu\gamma_5 */
    _complex_times_vector(phi1, z, r->s0);
    _complex_times_vector(phi2, z, r->s1);
    _complex_times_vector(phi3, w, r->s2);
    _complex_times_vector(phi4, w, r->s3);
    /* Subtract s and store the result in t */
    /* multiply with  gamma5 included by    */
    /* reversed order of s and phi3|4       */
    _vector_sub(t->s0, phi1, s->s0);
    _vector_sub(t->s1, phi2, s->s1);
    _vector_sub(t->s2, s->s2, phi3);
    _vector_sub(t->s3, s->s3, phi4);
  }

#ifdef OMP
  } /* OpenMP closing brace */
#endif
}

void mul_one_pm_imu_sub_mul(spinor * const l, spinor * const k, 
			    spinor * const j, const double _sign, const int N){
#ifdef OMP
#pragma omp parallel
  {
#endif
  _Complex double z,w;
  int ix;
  double sign=1.;
  spinor *r, *s, *t;

#if (!defined SSE2 && !defined SSE3)

  su3_vector ALIGN phi1, phi2, phi3, phi4;
  
#endif

  if(_sign < 0.){
    sign = -1.;
  }

  z = 1. + (sign * g_mu) * I;
  w = conj(z);
  /************ loop over all lattice sites ************/
#ifdef OMP
#pragma omp for
#endif
  for(ix = 0; ix < N; ix++){
    r = k+ix;
    s = j+ix;
    t = l+ix;
    /* Multiply the spinorfield with 1+imu\gamma_5 */
#if (defined SSE2 || defined SSE3)
    _prefetch_spinor((r+predist));
    _prefetch_spinor((s+predist));
    _sse_load_up(r->s0);
    _sse_vector_cmplx_mul(z);
    _sse_load(s->s0);
    _sse_vector_sub_up();
    _sse_store_nt_up(t->s0);
    _sse_load_up(r->s1);
    _sse_vector_cmplx_mul_two();
    _sse_load(s->s1);
    _sse_vector_sub_up();
    _sse_store_nt_up(t->s1);
    _sse_load_up(r->s2);
    _sse_vector_cmplx_mul(w);
    _sse_load(s->s2);
    _sse_vector_sub_up();
    _sse_store_nt_up(t->s2);
    _sse_load_up(r->s3);
    _sse_vector_cmplx_mul_two();
    _sse_load(s->s3);
    _sse_vector_sub_up();
    _sse_store_nt_up(t->s3);
#else
    _complex_times_vector(phi1, z, r->s0);
    _complex_times_vector(phi2, z, r->s1);
    _complex_times_vector(phi3, w, r->s2);
    _complex_times_vector(phi4, w, r->s3);
    /* Subtract s and store the result in t */
    _vector_sub(t->s0, phi1, s->s0);
    _vector_sub(t->s1, phi2, s->s1);
    _vector_sub(t->s2, phi3, s->s2);
    _vector_sub(t->s3, phi4, s->s3);
#endif
  }

#ifdef OMP
  } /* OpenMP closing brace */
#endif
}


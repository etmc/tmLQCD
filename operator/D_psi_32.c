/***********************************************************************
 *
 * Copyright (C) 2001 Martin Luescher
 * original code 
 * changed and extended for twisted mass 2002 Andrea Shindler
 *               2007,2008 Carsten Urbach
 *
 *  32 bit version 2015 Florian Burger 
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
 *
 * Action of a Dirac operator D (Wilson or twisted) on a given spinor field
 *
 *
 *******************************************************************************/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "global.h"
#include "su3.h"
#include "boundary.h"
#ifdef MPI
# include "xchange/xchange.h"
#endif
#include "update_backward_gauge.h"
#include "operator/D_psi_32.h"



static inline void p0add32(spinor32 * restrict const tmpr , spinor32 const * restrict const s, 
			 su3_32 const * restrict const u, const _Complex float phase) {

#ifdef OMP
#define static
#endif
  static su3_vector32 chi, psi;
#ifdef OMP
#undef static
#endif

  _vector_add(psi,s->s0, s->s2);
  _su3_multiply(chi, (*u), psi);

  _complex_times_vector(psi, phase, chi);
  _vector_add_assign(tmpr->s0, psi);
  _vector_add_assign(tmpr->s2, psi);

  _vector_add(psi, s->s1, s->s3);
  _su3_multiply(chi, (*u), psi);

  _complex_times_vector(psi, phase, chi);
  _vector_add_assign(tmpr->s1, psi);
  _vector_add_assign(tmpr->s3, psi);

  return;
}


static inline void m0add32(spinor32 * restrict const tmpr, spinor32 const * restrict const s, 
			 su3_32 const * restrict const u, const _Complex float phase) {
#ifdef OMP
#define static
#endif
  static su3_vector32 chi, psi;
#ifdef OMP
#undef static
#endif

  _vector_sub(psi, s->s0, s->s2);
  _su3_inverse_multiply(chi, (*u), psi);

  _complexcjg_times_vector(psi, phase, chi);
  _vector_add_assign(tmpr->s0, psi);
  _vector_sub_assign(tmpr->s2, psi);

  _vector_sub(psi, s->s1, s->s3);
  _su3_inverse_multiply(chi, (*u), psi);

  _complexcjg_times_vector(psi, phase, chi);
  _vector_add_assign(tmpr->s1, psi);
  _vector_sub_assign(tmpr->s3, psi);

  return;
}

static inline void p1add32(spinor32 * restrict const tmpr, spinor32 const * restrict const s, 
			 su3_32 const * restrict const u, const _Complex float phase) {
#ifdef OMP
#define static
#endif
  static su3_vector32 chi, psi;
#ifdef OMP
#undef static
#endif

  _vector_i_add(psi,s->s0,s->s3);
  _su3_multiply(chi,(*u),psi);

  _complex_times_vector(psi, phase, chi);
  _vector_add_assign(tmpr->s0, psi);
  _vector_i_sub_assign(tmpr->s3, psi);
 
  _vector_i_add(psi, s->s1, s->s2);
  _su3_multiply(chi, (*u), psi);

  _complex_times_vector(psi, phase, chi);
  _vector_add_assign(tmpr->s1, psi);
  _vector_i_sub_assign(tmpr->s2, psi);

  return;
}

static inline void m1add32(spinor32 * restrict const tmpr, spinor32 const * restrict const s, 
			 su3_32 const * restrict const u, const _Complex float phase) {
#ifdef OMP
#define static
#endif
  static su3_vector32 chi, psi;
#ifdef OMP
#undef static
#endif

  _vector_i_sub(psi,s->s0, s->s3);
  _su3_inverse_multiply(chi,(*u), psi);

  _complexcjg_times_vector(psi, phase, chi);
  _vector_add_assign(tmpr->s0, psi);
  _vector_i_add_assign(tmpr->s3, psi);

  _vector_i_sub(psi, s->s1, s->s2);
  _su3_inverse_multiply(chi, (*u), psi);

  _complexcjg_times_vector(psi, phase, chi);
  _vector_add_assign(tmpr->s1, psi);
  _vector_i_add_assign(tmpr->s2, psi);

  return;
}

static inline void p2add32(spinor32 * restrict const tmpr, spinor32 const * restrict const s, 
			 su3_32 const * restrict const u, const _Complex float phase) {
#ifdef OMP
#define static
#endif
  static su3_vector32 chi, psi;
#ifdef OMP
#undef static
#endif

  _vector_add(psi,s->s0,s->s3);
  _su3_multiply(chi, (*u), psi);

  _complex_times_vector(psi, phase, chi);
  _vector_add_assign(tmpr->s0, psi);
  _vector_add_assign(tmpr->s3, psi);

  _vector_sub(psi,s->s1,s->s2);
  _su3_multiply(chi, (*u), psi);

  _complex_times_vector(psi, phase, chi);
  _vector_add_assign(tmpr->s1, psi);
  _vector_sub_assign(tmpr->s2, psi);


  return;
}

static inline void m2add32(spinor32 * restrict const tmpr, spinor32 const * restrict const s, 
			 su3_32 const * restrict const u, const _Complex float phase) {
#ifdef OMP
#define static
#endif
  static su3_vector32 chi, psi;
#ifdef OMP
#undef static
#endif

  _vector_sub(psi, s->s0, s->s3);
  _su3_inverse_multiply(chi, (*u), psi);

  _complexcjg_times_vector(psi, phase, chi);
  _vector_add_assign(tmpr->s0, psi);
  _vector_sub_assign(tmpr->s3, psi);

  _vector_add(psi, s->s1, s->s2);
  _su3_inverse_multiply(chi, (*u),psi);

  _complexcjg_times_vector(psi, phase, chi);
  _vector_add_assign(tmpr->s1, psi);
  _vector_add_assign(tmpr->s2, psi);

  return;
}

static inline void p3add32(spinor32 * restrict const tmpr, spinor32 const * restrict const s, 
			 su3_32 const * restrict const u, const _Complex float phase) {
#ifdef OMP
#define static
#endif
  static su3_vector32 chi, psi;
#ifdef OMP
#undef static
#endif

  _vector_i_add(psi, s->s0, s->s2);
  _su3_multiply(chi, (*u), psi);

  _complex_times_vector(psi, phase, chi);
  _vector_add_assign(tmpr->s0, psi);
  _vector_i_sub_assign(tmpr->s2, psi);

  _vector_i_sub(psi,s->s1, s->s3);
  _su3_multiply(chi, (*u), psi);

  _complex_times_vector(psi, phase, chi);
  _vector_add_assign(tmpr->s1, psi);
  _vector_i_add_assign(tmpr->s3, psi);

  return;
}

static inline void m3addandstore32(spinor32 * restrict const r, spinor32 const * restrict const s, 
				 su3_32 const * restrict const u, const _Complex float phase,
         spinor32 const * restrict const tmpr) {
#ifdef OMP
#define static
#endif
  static su3_vector32 chi, psi;
#ifdef OMP
#undef static
#endif

  _vector_i_sub(psi,s->s0, s->s2);
  _su3_inverse_multiply(chi, (*u), psi);

  _complexcjg_times_vector(psi, phase, chi);
  _vector_add(r->s0, tmpr->s0, psi);
  _vector_i_add(r->s2, tmpr->s2, psi);

  _vector_i_add(psi, s->s1, s->s3);
  _su3_inverse_multiply(chi, (*u), psi);

  _complexcjg_times_vector(psi, phase, chi);
  _vector_add(r->s1, tmpr->s1, psi);
  _vector_i_sub(r->s3, tmpr->s3, psi);

  return;
}




void D_psi_32(spinor32 * const P, spinor32 * const Q){
  if(P==Q){
    printf("Error in D_psi (operator.c):\n");
    printf("Arguments must be different spinor fields\n");
    printf("Program aborted\n");
    exit(1);
  }
//convert phases to float locally
_Complex float ALIGN32 phase_0_32 = (_Complex float) phase_0;
_Complex float ALIGN32 phase_1_32 = (_Complex float) phase_1;
_Complex float ALIGN32 phase_2_32 = (_Complex float) phase_2;
_Complex float ALIGN32 phase_3_32 = (_Complex float) phase_3;  

#ifdef _GAUGE_COPY
  if(g_update_gauge_copy_32) {
      update_backward_gauge_32(g_gauge_field_32);
  }
#endif
# if defined MPI
  xchange_lexicfield32(Q);
# endif

#ifdef OMP
#pragma omp parallel
  {
#endif

  int ix,iy;
  su3_32 * restrict up,* restrict um;
  spinor32 * restrict rr; 
  spinor32 const * restrict s;
  spinor32 const * restrict sp;
  spinor32 const * restrict sm;
  _Complex float rho1, rho2;
  spinor32 tmpr;

  rho1 = 1.f + (float) g_mu * I;
  rho2 = conj(rho1);

  /************************ loop over all lattice sites *************************/

#ifdef OMP
#pragma omp for
#endif
  for (ix=0;ix<VOLUME;ix++)
  {
    rr  = (spinor32 *) P +ix;
    s  = (spinor32 *) Q +ix;

    _complex_times_vector(tmpr.s0, rho1, s->s0);
    _complex_times_vector(tmpr.s1, rho1, s->s1);
    _complex_times_vector(tmpr.s2, rho2, s->s2);
    _complex_times_vector(tmpr.s3, rho2, s->s3);

    /******************************* direction +0 *********************************/
    iy=g_iup[ix][0];
    sp = (spinor32 *) Q +iy;
    up=&g_gauge_field_32[ix][0];
    p0add32(&tmpr, sp, up, phase_0_32);

    /******************************* direction -0 *********************************/
    iy=g_idn[ix][0];
    sm  = (spinor32 *) Q +iy;
    um=&g_gauge_field_32[iy][0];
    m0add32(&tmpr, sm, um, phase_0_32);

    /******************************* direction +1 *********************************/
    iy=g_iup[ix][1];
    sp = (spinor32 *) Q +iy;
    up=&g_gauge_field_32[ix][1];
    p1add32(&tmpr, sp, up, phase_1_32);

    /******************************* direction -1 *********************************/
    iy=g_idn[ix][1];
    sm = (spinor32 *) Q +iy;
    um=&g_gauge_field_32[iy][1];
    m1add32(&tmpr, sm, um, phase_1_32);

    /******************************* direction +2 *********************************/
    iy=g_iup[ix][2];
    sp = (spinor32 *) Q +iy;
    up=&g_gauge_field_32[ix][2];
    p2add32(&tmpr, sp, up, phase_2_32);

    /******************************* direction -2 *********************************/
    iy=g_idn[ix][2];
    sm = (spinor32 *) Q +iy;
    um=&g_gauge_field_32[iy][2];
    m2add32(&tmpr, sm, um, phase_2_32);

    /******************************* direction +3 *********************************/
    iy=g_iup[ix][3];
    sp = (spinor32 *) Q +iy;
    up=&g_gauge_field_32[ix][3];
    p3add32(&tmpr, sp, up, phase_3_32);

    /******************************* direction -3 *********************************/
    iy=g_idn[ix][3];
    sm = (spinor32 *) Q +iy;
    um=&g_gauge_field_32[iy][3];
    m3addandstore32(rr, sm, um, phase_3_32, &tmpr);
  }
#ifdef OMP
  } /* OpenMP closing brace */
#endif
}


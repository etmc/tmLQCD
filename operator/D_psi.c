/***********************************************************************
 *
 * Copyright (C) 2001 Martin Luescher
 * original code 
 * changed and extended for twisted mass 2002 Andrea Shindler
 *               2007,2008 Carsten Urbach
 *
 * Blue Gene version Copyright (C) 2007 Carsten Urbach 
 * Block Dirac operator Copyright (C) 2008 Carsten Urbach
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
 * various versions including a block version.
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
#include "sse.h"
#include "boundary.h"
#include "gamma.h"
#include "linalg_eo.h"
#ifdef MPI
# include "xchange/xchange.h"
#endif
#include "update_backward_gauge.h"
#include "block.h"
#include "operator/D_psi.h"
#include "operator/clovertm_operators.h"
#include "solver/dirac_operator_eigenvectors.h"

#if (defined SSE23 || defined SSE33)

#else


static inline void p0add(spinor * restrict const tmpr , spinor const * restrict const s, 
                         su3 const * restrict const u, const _Complex double phase) {
#ifdef OMP
#define static
#endif
  static su3_vector chi, psi;
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


static inline void m0add(spinor * restrict const tmpr, spinor const * restrict const s, 
                         su3 const * restrict const u, const _Complex double phase) {
#ifdef OMP
#define static
#endif
  static su3_vector chi, psi;
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

static inline void p1add(spinor * restrict const tmpr, spinor const * restrict const s, 
                         su3 const * restrict const u, const _Complex double phase) {
#ifdef OMP
#define static
#endif
  static su3_vector chi, psi;
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

static inline void m1add(spinor * restrict const tmpr, spinor const * restrict const s, 
                         su3 const * restrict const u, const _Complex double phase) {
#ifdef OMP
#define static
#endif
  static su3_vector chi, psi;
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

static inline void p2add(spinor * restrict const tmpr, spinor const * restrict const s, 
                         su3 const * restrict const u, const _Complex double phase) {
#ifdef OMP
#define static
#endif
  static su3_vector chi, psi;
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

static inline void m2add(spinor * restrict const tmpr, spinor const * restrict const s, 
                         su3 const * restrict const u, const _Complex double phase) {
#ifdef OMP
#define static
#endif
  static su3_vector chi, psi;
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

static inline void p3add(spinor * restrict const tmpr, spinor const * restrict const s, 
                         su3 const * restrict const u, const _Complex double phase) {
#ifdef OMP
#define static
#endif
  static su3_vector chi, psi;
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

static inline void m3addandstore(spinor * restrict const r, spinor const * restrict const s, 
                                 su3 const * restrict const u, const _Complex double phase,
                                 spinor const * restrict const tmpr) {
#ifdef OMP
#define static
#endif
  static su3_vector chi, psi;
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

/* this is the hopping part only */
static inline void local_H(spinor * const rr, spinor const * const s, su3 const * restrict u, 
                           int * _idx, spinor * const restrict tmpr) {

  int * idx = _idx;

  /****** direction +0 ******/
  p0add(tmpr, s + (*idx), u, phase_0);
  u++;
  idx++;
  /****** direction -0 ******/
  m0add(tmpr, s + (*idx), u, phase_0);
  u++;
  idx++;
  /****** direction +1 ******/
  p1add(tmpr, s + (*idx), u, phase_1);
  u++;
  idx++;
  /****** direction -1 ******/
  m1add(tmpr, s + (*idx), u, phase_1);
  u++;
  idx++;
  /****** direction +2 ******/
  p2add(tmpr, s + (*idx), u, phase_2);
  u++;
  idx++;
  /****** direction -2 ******/
  m2add(tmpr, s + (*idx), u, phase_2);
  u++;
  idx++;
  /****** direction +3 ******/
  p3add(tmpr, s + (*idx), u, phase_3);
  u++;
  idx++;
  /****** direction -3 ******/
  m3addandstore(rr, s + (*idx), u, phase_3, tmpr);

  return;
}


#endif

#if (defined SSE2 || defined SSE3)

/* Serially Checked ! */
void Dtm_psi(spinor * const P, spinor * const Q){

  if(P==Q){
    printf("Error in Dtm_psi (D_psi.c):\n");
    printf("Arguments must be differen spinor fields\n");
    printf("Program aborted\n");
    exit(1);
  }

#ifdef _GAUGE_COPY2
  if(g_update_gauge_copy) {
    update_backward_gauge(g_gauge_field);
  }
#endif

# if defined MPI
  xchange_lexicfield(Q);
# endif

#ifdef OMP
#pragma omp parallel
  {
#endif
    int ix,iy,iz;
    su3 *up,*um;
    spinor *s,*sp,*sm,*rn;
    _Complex double fact1, fact2;
    spinor rs __attribute__ ((aligned (16)));

    fact1 = 1. + g_mu * I;
    fact2 = conj(fact1);

#ifndef OMP
    iy=g_iup[0][0];
    sp=(spinor *) Q + iy;
    up=&g_gauge_field[0][0];
#endif

    /************************ loop over all lattice sites *************************/
#ifdef OMP
#pragma omp for
#endif
    for (ix=0;ix<VOLUME;ix++){
#ifdef OMP
      iy=g_iup[ix][0];
      up=&g_gauge_field[ix][0];
      sp=(spinor *) Q + iy;
#endif
      s=(spinor *) Q + ix;
      _prefetch_spinor(s);

      /******************************* direction +0 *********************************/

      iy=g_idn[ix][0];
      
      sm = (spinor *) Q + iy;
      _prefetch_spinor(sm);       

      _sse_load(sp->s0);
      _sse_load_up(sp->s2);
      _sse_vector_add();

      _sse_su3_multiply((*up));
      _sse_vector_cmplx_mul(phase_0);
      _sse_store_up(rs.s2);

      // the diagonal bit
      _sse_load_up(s->s0);
      _sse_vector_cmplx_mul(fact1);
      _sse_load(rs.s2);
      _sse_vector_add();
      _sse_store(rs.s0);

      // g5 in the twisted term
      _sse_load_up(s->s2);
      _sse_vector_cmplx_mul(fact2);
      _sse_load(rs.s2);
      _sse_vector_add();
      _sse_store(rs.s2);      
      
      um=&g_gauge_field[iy][0];
      _prefetch_su3(um);
      
      _sse_load(sp->s1);
      _sse_load_up(sp->s3);
      _sse_vector_add();
      
      _sse_su3_multiply((*up));
      _sse_vector_cmplx_mul(phase_0);
      _sse_store_up(rs.s3);
    
      // the diagonal bit
      _sse_load_up(s->s1);
      _sse_vector_cmplx_mul(fact1);
      _sse_load(rs.s3);
      _sse_vector_add();
      _sse_store(rs.s1);

      // g5 in the twisted term
      _sse_load_up(s->s3);
      _sse_vector_cmplx_mul(fact2);
      _sse_load(rs.s3);
      _sse_vector_add();
      _sse_store(rs.s3); 

      /******************************* direction -0 *********************************/

      iy=g_iup[ix][1];

      sp = (spinor *) Q + iy;
      _prefetch_spinor(sp);

      _sse_load(sm->s0);
      _sse_load_up(sm->s2);
      _sse_vector_sub();
      
      _sse_su3_inverse_multiply((*um));
      _sse_vector_cmplxcg_mul(phase_0);
      _sse_load(rs.s0);
      _sse_vector_add();
      _sse_store(rs.s0);

      _sse_load(rs.s2);
      _sse_vector_sub();
      _sse_store(rs.s2);
      
      up+=1;
      _prefetch_su3(up);
      
      _sse_load(sm->s1);
      _sse_load_up(sm->s3);
      _sse_vector_sub();
      
      _sse_su3_inverse_multiply((*um));
      _sse_vector_cmplxcg_mul(phase_0);
      _sse_load(rs.s1);
      _sse_vector_add();
      _sse_store(rs.s1);

      _sse_load(rs.s3);
      _sse_vector_sub();
      _sse_store(rs.s3);
      
      /******************************* direction +1 *********************************/

      iy=g_idn[ix][1];
      
      sm = (spinor *) Q + iy;
      _prefetch_spinor(sm);

      _sse_load(sp->s0);
      _sse_load_up(sp->s3);
      _sse_vector_i_mul();
      _sse_vector_add();

      _sse_su3_multiply((*up));
      _sse_vector_cmplx_mul(phase_1);
      _sse_load(rs.s0);
      _sse_vector_add();
      _sse_store(rs.s0);

      _sse_load(rs.s3);
      _sse_vector_i_mul();      
      _sse_vector_sub();
      _sse_store(rs.s3); 
      
      um=&g_gauge_field[iy][1];
      _prefetch_su3(um);

      _sse_load(sp->s1);
      _sse_load_up(sp->s2);
      _sse_vector_i_mul();
      _sse_vector_add();

      _sse_su3_multiply((*up));
      _sse_vector_cmplx_mul(phase_1);
      _sse_load(rs.s1);
      _sse_vector_add();
      _sse_store(rs.s1);

      _sse_load(rs.s2);
      _sse_vector_i_mul();      
      _sse_vector_sub();
      _sse_store(rs.s2);       

      /******************************* direction -1 *********************************/

      iy=g_iup[ix][2];

      sp = (spinor *) Q + iy;
      _prefetch_spinor(sp);

      _sse_load(sm->s0);
      _sse_load_up(sm->s3);
      _sse_vector_i_mul();
      _sse_vector_sub();
      
      _sse_su3_inverse_multiply((*um));
      _sse_vector_cmplxcg_mul(phase_1);
      _sse_load(rs.s0);
      _sse_vector_add();
      _sse_store(rs.s0);

      _sse_load(rs.s3);
      _sse_vector_i_mul();      
      _sse_vector_add();
      _sse_store(rs.s3);

      up+=1;
      _prefetch_su3(up);

      _sse_load(sm->s1);
      _sse_load_up(sm->s2);
      _sse_vector_i_mul();
      _sse_vector_sub();
      
      _sse_su3_inverse_multiply((*um));
      _sse_vector_cmplxcg_mul(phase_1);
      _sse_load(rs.s1);
      _sse_vector_add();
      _sse_store(rs.s1);

      _sse_load(rs.s2);
      _sse_vector_i_mul();      
      _sse_vector_add();
      _sse_store(rs.s2);

      /******************************* direction +2 *********************************/

      iy=g_idn[ix][2];

      sm = (spinor *) Q + iy;
      _prefetch_spinor(sm);

      _sse_load(sp->s0);
      _sse_load_up(sp->s3);
      _sse_vector_add();

      _sse_su3_multiply((*up));
      _sse_vector_cmplx_mul(phase_2);
      _sse_load(rs.s0);
      _sse_vector_add();
      _sse_store(rs.s0);

      _sse_load(rs.s3);
      _sse_vector_add();
      _sse_store(rs.s3);
      
      um=&g_gauge_field[iy][2];
      _prefetch_su3(um);

      _sse_load(sp->s1);
      _sse_load_up(sp->s2);
      _sse_vector_sub();

      _sse_su3_multiply((*up));
      _sse_vector_cmplx_mul(phase_2);
      _sse_load(rs.s1);
      _sse_vector_add();
      _sse_store(rs.s1);

      _sse_load(rs.s2);
      _sse_vector_sub();
      _sse_store(rs.s2);      

      /******************************* direction -2 *********************************/

      iy=g_iup[ix][3];

      sp = (spinor *) Q + iy;
      _prefetch_spinor(sp);

      _sse_load(sm->s0);
      _sse_load_up(sm->s3);
      _sse_vector_sub();
      
      _sse_su3_inverse_multiply((*um));
      _sse_vector_cmplxcg_mul(phase_2);
      _sse_load(rs.s0);
      _sse_vector_add();
      _sse_store(rs.s0);

      _sse_load(rs.s3);
      _sse_vector_sub();
      _sse_store(rs.s3);
      
      up+=1;
      _prefetch_su3(up);

      _sse_load(sm->s1);
      _sse_load_up(sm->s2);
      _sse_vector_add();
      
      _sse_su3_inverse_multiply((*um));
      _sse_vector_cmplxcg_mul(phase_2);
      _sse_load(rs.s1);
      _sse_vector_add();
      _sse_store(rs.s1);

      _sse_load(rs.s2);
      _sse_vector_add();
      _sse_store(rs.s2);      
      
      /******************************* direction +3 *********************************/

      iy=g_idn[ix][3];

      sm = (spinor *) Q + iy;
      _prefetch_spinor(sm);

      _sse_load(sp->s0);
      _sse_load_up(sp->s2);
      _sse_vector_i_mul();
      _sse_vector_add();

      _sse_su3_multiply((*up));
      _sse_vector_cmplx_mul(phase_3);
      _sse_load(rs.s0);
      _sse_vector_add();
      _sse_store(rs.s0);

      _sse_load(rs.s2);
      _sse_vector_i_mul();      
      _sse_vector_sub();
      _sse_store(rs.s2);
      
      um=&g_gauge_field[iy][3];
      _prefetch_su3(um);

      _sse_load(sp->s1);
      _sse_load_up(sp->s3);
      _sse_vector_i_mul();
      _sse_vector_sub();

      _sse_su3_multiply((*up));
      _sse_vector_cmplx_mul(phase_3);
      _sse_load(rs.s1);
      _sse_vector_add();
      _sse_store(rs.s1);

      _sse_load(rs.s3);
      _sse_vector_i_mul();      
      _sse_vector_add();
      _sse_store(rs.s3);
      
      /******************************* direction -3 *********************************/

      iz=(ix+1+VOLUME)%VOLUME;

      iy=g_iup[iz][0];
      
      sp = (spinor *) Q + iy;
      _prefetch_spinor(sp);

      _sse_load(sm->s0);
      _sse_load_up(sm->s2);
      _sse_vector_i_mul();
      _sse_vector_sub();
      
      _sse_su3_inverse_multiply((*um));
      _sse_vector_cmplxcg_mul(phase_3);
      rn = (spinor *) P + ix;
      
      _sse_load(rs.s0);
      _sse_vector_add();
      _sse_store_nt(rn->s0);

      _sse_load(rs.s2);
      _sse_vector_i_mul();      
      _sse_vector_add();
      _sse_store_nt(rn->s2);

      up=&g_gauge_field[iz][0];
      _prefetch_su3(up);

      _sse_load(sm->s1);
      _sse_load_up(sm->s3);
      _sse_vector_i_mul();
      _sse_vector_add();
      
      _sse_su3_inverse_multiply((*um));
      _sse_vector_cmplxcg_mul(phase_3);
      _sse_load(rs.s1);
      _sse_vector_add();
      _sse_store_nt(rn->s1);

      _sse_load(rs.s3);
      _sse_vector_i_mul();      
      _sse_vector_sub();
      _sse_store_nt(rn->s3);
      
      /******************************** end of loop *********************************/

    }
#ifdef OMP
  } /* OpenMP closing brace */
#endif
}

void Dsw_psi(spinor * const P, spinor * const Q){

  if(P==Q){
    printf("Error in Dsw_psi (D_psi.c):\n");
    printf("Arguments must be differen spinor fields\n");
    printf("Program aborted\n");
    exit(1);
  }

#ifdef _GAUGE_COPY2
  if(g_update_gauge_copy) {
    update_backward_gauge(g_gauge_field);
  }
#endif

# if defined MPI
  xchange_lexicfield(Q);
# endif

#ifdef OMP
#pragma omp parallel
  {
#endif
    int ix,iy,iz;
    su3 *up,*um;
    spinor *s,*sp,*sm,*rn;
    spinor ALIGN stmp,rs;

#ifndef OMP
    iy=g_iup[0][0];
    sp=(spinor *) Q + iy;
    up=&g_gauge_field[0][0];
#endif

    /************************ loop over all lattice sites *************************/
#ifdef OMP
#pragma omp for
#endif
    for (ix=0;ix<VOLUME;ix++){
#ifdef OMP
      iy=g_iup[ix][0];
      up=&g_gauge_field[ix][0];
      sp=(spinor *) Q + iy;
#endif
      s=(spinor *) Q + ix;
      _prefetch_spinor(s);

      /******************************* direction +0 *********************************/

      iy=g_idn[ix][0];
      
      sm = (spinor *) Q + iy;
      _prefetch_spinor(sm);       

      _sse_load(sp->s0);
      _sse_load_up(sp->s2);
      _sse_vector_add();

      _sse_su3_multiply((*up));
      _sse_vector_cmplx_mul(phase_0);
      _sse_store_up(rs.s2);

      // apply the clover plus twisted term to diagonal bit
      assign_mul_one_sw_pm_imu_site_lexic(ix, &stmp, s, g_mu);

      _sse_load_up(stmp.s0);
      _sse_load(rs.s2);
      _sse_vector_add();
      _sse_store(rs.s0);

      _sse_load_up(stmp.s2);
      _sse_load(rs.s2);
      _sse_vector_add();
      _sse_store(rs.s2);      
      
      um=&g_gauge_field[iy][0];
      _prefetch_su3(um);
      
      _sse_load(sp->s1);
      _sse_load_up(sp->s3);
      _sse_vector_add();
      
      _sse_su3_multiply((*up));
      _sse_vector_cmplx_mul(phase_0);
      _sse_store_up(rs.s3);
    
      _sse_load_up(stmp.s1);
      _sse_load(rs.s3);
      _sse_vector_add();
      _sse_store(rs.s1);

      _sse_load_up(stmp.s3);
      _sse_load(rs.s3);
      _sse_vector_add();
      _sse_store(rs.s3); 

      /******************************* direction -0 *********************************/

      iy=g_iup[ix][1];

      sp = (spinor *) Q + iy;
      _prefetch_spinor(sp);

      _sse_load(sm->s0);
      _sse_load_up(sm->s2);
      _sse_vector_sub();
      
      _sse_su3_inverse_multiply((*um));
      _sse_vector_cmplxcg_mul(phase_0);
      _sse_load(rs.s0);
      _sse_vector_add();
      _sse_store(rs.s0);

      _sse_load(rs.s2);
      _sse_vector_sub();
      _sse_store(rs.s2);
      
      up+=1;
      _prefetch_su3(up);
      
      _sse_load(sm->s1);
      _sse_load_up(sm->s3);
      _sse_vector_sub();
      
      _sse_su3_inverse_multiply((*um));
      _sse_vector_cmplxcg_mul(phase_0);
      _sse_load(rs.s1);
      _sse_vector_add();
      _sse_store(rs.s1);

      _sse_load(rs.s3);
      _sse_vector_sub();
      _sse_store(rs.s3);
      
      /******************************* direction +1 *********************************/

      iy=g_idn[ix][1];
      
      sm = (spinor *) Q + iy;
      _prefetch_spinor(sm);

      _sse_load(sp->s0);
      _sse_load_up(sp->s3);
      _sse_vector_i_mul();
      _sse_vector_add();

      _sse_su3_multiply((*up));
      _sse_vector_cmplx_mul(phase_1);
      _sse_load(rs.s0);
      _sse_vector_add();
      _sse_store(rs.s0);

      _sse_load(rs.s3);
      _sse_vector_i_mul();      
      _sse_vector_sub();
      _sse_store(rs.s3); 
      
      um=&g_gauge_field[iy][1];
      _prefetch_su3(um);

      _sse_load(sp->s1);
      _sse_load_up(sp->s2);
      _sse_vector_i_mul();
      _sse_vector_add();

      _sse_su3_multiply((*up));
      _sse_vector_cmplx_mul(phase_1);
      _sse_load(rs.s1);
      _sse_vector_add();
      _sse_store(rs.s1);

      _sse_load(rs.s2);
      _sse_vector_i_mul();      
      _sse_vector_sub();
      _sse_store(rs.s2);       

      /******************************* direction -1 *********************************/

      iy=g_iup[ix][2];

      sp = (spinor *) Q + iy;
      _prefetch_spinor(sp);

      _sse_load(sm->s0);
      _sse_load_up(sm->s3);
      _sse_vector_i_mul();
      _sse_vector_sub();
      
      _sse_su3_inverse_multiply((*um));
      _sse_vector_cmplxcg_mul(phase_1);
      _sse_load(rs.s0);
      _sse_vector_add();
      _sse_store(rs.s0);

      _sse_load(rs.s3);
      _sse_vector_i_mul();      
      _sse_vector_add();
      _sse_store(rs.s3);

      up+=1;
      _prefetch_su3(up);

      _sse_load(sm->s1);
      _sse_load_up(sm->s2);
      _sse_vector_i_mul();
      _sse_vector_sub();
      
      _sse_su3_inverse_multiply((*um));
      _sse_vector_cmplxcg_mul(phase_1);
      _sse_load(rs.s1);
      _sse_vector_add();
      _sse_store(rs.s1);

      _sse_load(rs.s2);
      _sse_vector_i_mul();      
      _sse_vector_add();
      _sse_store(rs.s2);

      /******************************* direction +2 *********************************/

      iy=g_idn[ix][2];

      sm = (spinor *) Q + iy;
      _prefetch_spinor(sm);

      _sse_load(sp->s0);
      _sse_load_up(sp->s3);
      _sse_vector_add();

      _sse_su3_multiply((*up));
      _sse_vector_cmplx_mul(phase_2);
      _sse_load(rs.s0);
      _sse_vector_add();
      _sse_store(rs.s0);

      _sse_load(rs.s3);
      _sse_vector_add();
      _sse_store(rs.s3);
      
      um=&g_gauge_field[iy][2];
      _prefetch_su3(um);

      _sse_load(sp->s1);
      _sse_load_up(sp->s2);
      _sse_vector_sub();

      _sse_su3_multiply((*up));
      _sse_vector_cmplx_mul(phase_2);
      _sse_load(rs.s1);
      _sse_vector_add();
      _sse_store(rs.s1);

      _sse_load(rs.s2);
      _sse_vector_sub();
      _sse_store(rs.s2);      

      /******************************* direction -2 *********************************/

      iy=g_iup[ix][3];

      sp = (spinor *) Q + iy;
      _prefetch_spinor(sp);

      _sse_load(sm->s0);
      _sse_load_up(sm->s3);
      _sse_vector_sub();
      
      _sse_su3_inverse_multiply((*um));
      _sse_vector_cmplxcg_mul(phase_2);
      _sse_load(rs.s0);
      _sse_vector_add();
      _sse_store(rs.s0);

      _sse_load(rs.s3);
      _sse_vector_sub();
      _sse_store(rs.s3);
      
      up+=1;
      _prefetch_su3(up);

      _sse_load(sm->s1);
      _sse_load_up(sm->s2);
      _sse_vector_add();
      
      _sse_su3_inverse_multiply((*um));
      _sse_vector_cmplxcg_mul(phase_2);
      _sse_load(rs.s1);
      _sse_vector_add();
      _sse_store(rs.s1);

      _sse_load(rs.s2);
      _sse_vector_add();
      _sse_store(rs.s2);      
      
      /******************************* direction +3 *********************************/

      iy=g_idn[ix][3];

      sm = (spinor *) Q + iy;
      _prefetch_spinor(sm);

      _sse_load(sp->s0);
      _sse_load_up(sp->s2);
      _sse_vector_i_mul();
      _sse_vector_add();

      _sse_su3_multiply((*up));
      _sse_vector_cmplx_mul(phase_3);
      _sse_load(rs.s0);
      _sse_vector_add();
      _sse_store(rs.s0);

      _sse_load(rs.s2);
      _sse_vector_i_mul();      
      _sse_vector_sub();
      _sse_store(rs.s2);
      
      um=&g_gauge_field[iy][3];
      _prefetch_su3(um);

      _sse_load(sp->s1);
      _sse_load_up(sp->s3);
      _sse_vector_i_mul();
      _sse_vector_sub();

      _sse_su3_multiply((*up));
      _sse_vector_cmplx_mul(phase_3);
      _sse_load(rs.s1);
      _sse_vector_add();
      _sse_store(rs.s1);

      _sse_load(rs.s3);
      _sse_vector_i_mul();      
      _sse_vector_add();
      _sse_store(rs.s3);
      
      /******************************* direction -3 *********************************/

      iz=(ix+1+VOLUME)%VOLUME;

      iy=g_iup[iz][0];
      
      sp = (spinor *) Q + iy;
      _prefetch_spinor(sp);

      _sse_load(sm->s0);
      _sse_load_up(sm->s2);
      _sse_vector_i_mul();
      _sse_vector_sub();
      
      _sse_su3_inverse_multiply((*um));
      _sse_vector_cmplxcg_mul(phase_3);
      rn = (spinor *) P + ix;
      
      _sse_load(rs.s0);
      _sse_vector_add();
      _sse_store_nt(rn->s0);

      _sse_load(rs.s2);
      _sse_vector_i_mul();      
      _sse_vector_add();
      _sse_store_nt(rn->s2);

      up=&g_gauge_field[iz][0];
      _prefetch_su3(up);

      _sse_load(sm->s1);
      _sse_load_up(sm->s3);
      _sse_vector_i_mul();
      _sse_vector_add();
      
      _sse_su3_inverse_multiply((*um));
      _sse_vector_cmplxcg_mul(phase_3);
      _sse_load(rs.s1);
      _sse_vector_add();
      _sse_store_nt(rn->s1);

      _sse_load(rs.s3);
      _sse_vector_i_mul();      
      _sse_vector_sub();
      _sse_store_nt(rn->s3);
      
      /******************************** end of loop *********************************/

    }
#ifdef OMP
  } /* OpenMP closing brace */
#endif
}


#else

/* Serially Checked ! */

void Dtm_psi(spinor * const P, spinor * const Q){
  if(P==Q){
    printf("Error in Dtm_psi (operator.c):\n");
    printf("Arguments must be different spinor fields\n");
    printf("Program aborted\n");
    exit(1);
  }
#ifdef _GAUGE_COPY
  if(g_update_gauge_copy) {
    update_backward_gauge(g_gauge_field);
  }
#endif
# if defined MPI
  xchange_lexicfield(Q);
# endif

#ifdef OMP
#pragma omp parallel
  {
#endif

    int ix,iy;
    su3 * restrict up,* restrict um;
    spinor * restrict rr; 
    spinor const * restrict s;
    spinor const * restrict sp;
    spinor const * restrict sm;
    _Complex double rho1, rho2;
    spinor tmpr;

    rho1 = 1. + g_mu * I;
    rho2 = conj(rho1);

    /************************ loop over all lattice sites *************************/

#ifdef OMP
#pragma omp for
#endif
    for (ix=0;ix<VOLUME;ix++) {
      rr  = (spinor *) P +ix;
      s  = (spinor *) Q +ix;
      
      _complex_times_vector(tmpr.s0, rho1, s->s0);
      _complex_times_vector(tmpr.s1, rho1, s->s1);
      _complex_times_vector(tmpr.s2, rho2, s->s2);
      _complex_times_vector(tmpr.s3, rho2, s->s3);
      
      /******************************* direction +0 *********************************/
      iy=g_iup[ix][0];
      sp = (spinor *) Q +iy;
      up=&g_gauge_field[ix][0];
      p0add(&tmpr, sp, up, phase_0);
      
      /******************************* direction -0 *********************************/
      iy=g_idn[ix][0];
      sm  = (spinor *) Q +iy;
      um=&g_gauge_field[iy][0];
      m0add(&tmpr, sm, um, phase_0);
      
      /******************************* direction +1 *********************************/
      iy=g_iup[ix][1];
      sp = (spinor *) Q +iy;
      up=&g_gauge_field[ix][1];
      p1add(&tmpr, sp, up, phase_1);
      
      /******************************* direction -1 *********************************/
      iy=g_idn[ix][1];
      sm = (spinor *) Q +iy;
      um=&g_gauge_field[iy][1];
      m1add(&tmpr, sm, um, phase_1);
      
      /******************************* direction +2 *********************************/
      iy=g_iup[ix][2];
      sp = (spinor *) Q +iy;
      up=&g_gauge_field[ix][2];
      p2add(&tmpr, sp, up, phase_2);
      
      /******************************* direction -2 *********************************/
      iy=g_idn[ix][2];
      sm = (spinor *) Q +iy;
      um=&g_gauge_field[iy][2];
      m2add(&tmpr, sm, um, phase_2);
      
      /******************************* direction +3 *********************************/
      iy=g_iup[ix][3];
      sp = (spinor *) Q +iy;
      up=&g_gauge_field[ix][3];
      p3add(&tmpr, sp, up, phase_3);
      
      /******************************* direction -3 *********************************/
      iy=g_idn[ix][3];
      sm = (spinor *) Q +iy;
      um=&g_gauge_field[iy][3];
      m3addandstore(rr, sm, um, phase_3, &tmpr);
      if(ix < 0) {
        ix = 0;
      }
    }
#ifdef OMP
  } /* OpenMP closing brace */
#endif
}

void Dsw_psi(spinor * const P, spinor * const Q){
  if(P==Q) {
    printf("Error in Dsw_psi (D_psi.c):\n");
    printf("Arguments must be different spinor fields\n");
    printf("Program aborted\n");
    exit(1);
  }

#ifdef _GAUGE_COPY
  if(g_update_gauge_copy) {
    update_backward_gauge(g_gauge_field);
  }
#endif
# if defined MPI
  xchange_lexicfield(Q);
# endif

#ifdef OMP
#pragma omp parallel
  {
#endif
    
    int ix,iy;
    su3 * restrict up,* restrict um;
    spinor * restrict rr; 
    spinor const * restrict s;
    spinor const * restrict sp;
    spinor const * restrict sm;
    spinor tmpr;
    
    /************************ loop over all lattice sites *************************/
    
#ifdef OMP
#pragma omp for
#endif
    for (ix=0;ix<VOLUME;ix++)
      {
        rr  = (spinor *) P +ix;
        s  = (spinor *) Q +ix;
        assign_mul_one_sw_pm_imu_site_lexic(ix,&tmpr,s,g_mu);
        
        /******************************* direction +0 *********************************/
        iy=g_iup[ix][0];
        sp = (spinor *) Q +iy;
        up=&g_gauge_field[ix][0];
        p0add(&tmpr, sp, up, phase_0);
        
        /******************************* direction -0 *********************************/
        iy=g_idn[ix][0];
        sm  = (spinor *) Q +iy;
        um=&g_gauge_field[iy][0];
        m0add(&tmpr, sm, um, phase_0);
        
        /******************************* direction +1 *********************************/
        iy=g_iup[ix][1];
        sp = (spinor *) Q +iy;
        up=&g_gauge_field[ix][1];
        p1add(&tmpr, sp, up, phase_1);
        
        /******************************* direction -1 *********************************/
        iy=g_idn[ix][1];
        sm = (spinor *) Q +iy;
        um=&g_gauge_field[iy][1];
        m1add(&tmpr, sm, um, phase_1);
        
        /******************************* direction +2 *********************************/
        iy=g_iup[ix][2];
        sp = (spinor *) Q +iy;
        up=&g_gauge_field[ix][2];
        p2add(&tmpr, sp, up, phase_2);
        
        /******************************* direction -2 *********************************/
        iy=g_idn[ix][2];
        sm = (spinor *) Q +iy;
        um=&g_gauge_field[iy][2];
        m2add(&tmpr, sm, um, phase_2);
        
        /******************************* direction +3 *********************************/
        iy=g_iup[ix][3];
        sp = (spinor *) Q +iy;
        up=&g_gauge_field[ix][3];
        p3add(&tmpr, sp, up, phase_3);
        
        /******************************* direction -3 *********************************/
        iy=g_idn[ix][3];
        sm = (spinor *) Q +iy;
        um=&g_gauge_field[iy][3];
        m3addandstore(rr, sm, um, phase_3, &tmpr);
      }
#ifdef OMP
  } /* OpenMP closing brace */
#endif
}

#endif


void D_psi(spinor * const P, spinor * const Q){
   if(g_c_sw > 0)
     Dsw_psi(P,Q);
   else
     Dtm_psi(P,Q);

   return ;
}



void D_psi_prec(spinor * const P, spinor * const Q){

  /* todo: do preconditioning */
  spinorPrecWS *ws=(spinorPrecWS*)g_precWS;
  static _Complex double alpha = -1.0;

  alpha = -0.5;
  spinorPrecondition(P,Q,ws,T,L,alpha,0,1);
  D_psi(g_spinor_field[DUM_MATRIX],P);
  alpha = -0.5;
  spinorPrecondition(P,g_spinor_field[DUM_MATRIX],ws,T,L,alpha,0,1);
}

/* apply the Dirac operator to the block local spinor field s */
/* and store the result in block local spinor field rr        */
/* for block blk                                              */
/* the block local gauge field is assumed to be in the order  */
/* that is needed int local_D, which means also that it is a  */
/* double copy                                                */
// CU: has problems with SSE2,3
void Block_D_psi(block * blk, spinor * const rr, spinor * const s) {

  if(g_c_sw > 0)
     Block_Dsw_psi(blk,rr,s);
  else 
     Block_Dtm_psi(blk,rr,s);

  return;
}
 
//this version is for c_sw=0
void Block_Dtm_psi(block * blk, spinor * const rr, spinor * const s) {
  int i;
  spinor *r = rr;
  spinor *t = s;
  su3 * u = blk->u;
  int * idx = blk->idx;
  static _Complex double rhoa, rhob;
  spinor ALIGN tmpr;
  if(blk_gauge_eo) {
    init_blocks_gaugefield();
  }
  rhoa = 1.0 + g_mu * I;
  rhob = conj(rhoa);

  /* set the boundary term to zero */
  _spinor_null(rr[blk->volume]);
  _spinor_null(s[blk->volume]);

  for(i = 0; i < blk->volume; i++) {
    _complex_times_vector(tmpr.s0, rhoa, t->s0);
    _complex_times_vector(tmpr.s1, rhoa, t->s1);
    _complex_times_vector(tmpr.s2, rhob, t->s2);
    _complex_times_vector(tmpr.s3, rhob, t->s3);

    local_H(r, s, u, idx, &tmpr);

    r++;
    t++;
    idx += 8;
    u += 8;
  }

  return;
}


/* apply the Dirac operator to the block local spinor field s */
/* and store the result in block local spinor field rr        */
/* for block blk                                              */
/* the block local gauge field is assumed to be in the order  */
/* that is needed int local_D, which means also that it is a  */
/* double copy                                                */
// CU: has problems with SSE2,3
void Block_Dsw_psi(block * blk, spinor * const rr, spinor * const s) {
  int i;
  spinor *r = rr;
  spinor *t = s;
  su3 * u = blk->u;
  int * idx = blk->idx;
  //static _Complex double rhoa, rhob;
  spinor ALIGN tmpr;

  int it,ix,iy,iz; //lexiographic index of the site w.r.t the block
  int bt,bx,by,bz; //block coordinate on the local mpi process
  int dT,dX,dY,dZ; //block size
  int sT,sX,sY,sZ; //constant shifts
  int lx; //lexiographic index of the block site w.r.t the local mpi process

  dT = blk->BT;
  dX = blk->BLX;
  dY = blk->BLY;
  dZ = blk->BLZ;

  bt = blk->mpilocal_coordinate[0];
  bx = blk->mpilocal_coordinate[1];
  by = blk->mpilocal_coordinate[2];
  bz = blk->mpilocal_coordinate[3];

  sT = bt*dT;
  sX = bx*dX;
  sY = by*dY;
  sZ = bz*dZ;
  
  if(blk_gauge_eo) {
    init_blocks_gaugefield();
  }
  //rhoa = 1.0 + g_mu * I;
  //rhob = conj(rhoa);

  /* set the boundary term to zero */
  _spinor_null(rr[blk->volume]);
  _spinor_null(s[blk->volume]);

  for(i = 0; i < blk->volume; i++) {

    iz = i%dZ;
    iy = (i/dZ)%dY;
    ix = (i/(dZ*dY))%dX;
    it = i/(dZ*dY*dX);

    lx = g_ipt[it+sT][ix+sX][iy+sY][iz+sZ];

    assign_mul_one_sw_pm_imu_site_lexic(lx,&tmpr,t,g_mu);

    local_H(r, s, u, idx, &tmpr);

    r++;
    t++;
    idx += 8;
    u += 8;
  }

  return;
}


/* Apply Hopping Matrix to a even(odd) spinor */
void Block_H_psi(block * blk, spinor * const rr, spinor * const s, const int eo) {
  int i;
  spinor *r = rr;
  su3 * u = blk->u;
  int * eoidx = blk->evenidx;
  spinor ALIGN tmpr;

  if(!blk_gauge_eo) {
    init_blocks_eo_gaugefield();
  }

  /* for OE */
  if(eo == 1) {
    u = blk->u + blk->volume*8/2;
    eoidx = blk->oddidx;
  }

  /* set the boundary term to zero */
  _spinor_null(rr[blk->volume/2]);
  _spinor_null(s[blk->volume/2]);
  
  for(i = 0; i < blk->volume/2; i++) {
    _spinor_null(tmpr);

    local_H(r, s, u, eoidx, &tmpr);

    r++;
    eoidx += 8;
    u += 8;
  }

  return;
}

/* direction +t */
void boundary_D_0(spinor * const r, spinor * const s, su3 * const u) {

  static su3_vector chi, psi;

  _vector_add(psi,s->s0,s->s2);

  _su3_multiply(chi,(*u),psi);

  _complex_times_vector(r->s0, phase_0, chi);
  _vector_assign(r->s2,r->s0);

  _vector_add(psi,s->s1,s->s3);

  _su3_multiply(chi,(*u),psi);

  _complex_times_vector(r->s1, phase_0, chi);
  _vector_assign(r->s3, r->s1);

  return;
}

/* direction -t */
void boundary_D_1(spinor * const r, spinor * const s, su3 * restrict u) {

  static su3_vector chi, psi;

  _vector_sub(psi, s->s0, s->s2);

  _su3_inverse_multiply(chi, (*u), psi);

  _complexcjg_times_vector(r->s0, phase_0, chi);
  _vector_minus_assign(r->s2, r->s0);

  _vector_sub(psi,s->s1,s->s3);

  _su3_inverse_multiply(chi,(*u),psi);

  _complexcjg_times_vector(r->s1,phase_0,chi);
  _vector_minus_assign(r->s3, r->s1);

  return;
}

/* direction +x */
void boundary_D_2(spinor * const r, spinor * const s, su3 * restrict u) {

  static su3_vector chi, psi;

  _vector_i_add(psi,s->s0,s->s3);

  _su3_multiply(chi,(*u),psi);

  _complex_times_vector(r->s0, phase_1, chi);
  _vector_null(r->s3);
  _vector_i_sub_assign(r->s3, r->s0);

  _vector_i_add(psi,s->s1,s->s2);

  _su3_multiply(chi,(*u),psi);

  _complex_times_vector(r->s1, phase_1, chi);
  _vector_null(r->s2);
  _vector_i_sub_assign(r->s2, r->s1);

  return;
}

/* direction -x */
void boundary_D_3(spinor * const r, spinor * const s, su3 * restrict u) {

  static su3_vector chi, psi;

  _vector_i_sub(psi,s->s0,s->s3);

  _su3_inverse_multiply(chi,(*u),psi);

  _complexcjg_times_vector(r->s0, phase_1, chi);
  _vector_null(r->s3);
  _vector_i_add_assign(r->s3, r->s0);

  _vector_i_sub(psi,s->s1,s->s2);

  _su3_inverse_multiply(chi,(*u),psi);

  _complexcjg_times_vector(r->s1, phase_1, chi);
  _vector_null(r->s2);
  _vector_i_add_assign(r->s2, r->s1);

  return;
}

/* direction +y */
void boundary_D_4(spinor * const r, spinor * const s, su3 * restrict u) {

  static su3_vector chi, psi;

  _vector_add(psi,s->s0,s->s3);

  _su3_multiply(chi,(*u),psi);

  _complex_times_vector(r->s0, phase_2, chi);
  _vector_assign(r->s3, r->s0);

  _vector_sub(psi,s->s1,s->s2);

  _su3_multiply(chi,(*u),psi);

  _complex_times_vector(r->s1, phase_2, chi);
  _vector_minus_assign(r->s2, r->s1);

  return;
}

/* direction -y */
void boundary_D_5(spinor * const r, spinor * const s, su3 * restrict u) {

  static su3_vector chi, psi;

  _vector_sub(psi,s->s0,s->s3);

  _su3_inverse_multiply(chi,(*u),psi);

  _complexcjg_times_vector(r->s0, phase_2, chi);
  _vector_minus_assign(r->s3, r->s0);

  _vector_add(psi,s->s1,s->s2);

  _su3_inverse_multiply(chi,(*u),psi);

  _complexcjg_times_vector(r->s1, phase_2, chi);
  _vector_assign(r->s2, r->s1);


  return;
}

/* direction +z */
void boundary_D_6(spinor * const r, spinor * const s, su3 * restrict u) {

  static su3_vector chi, psi;

  _vector_i_add(psi,s->s0,s->s2);

  _su3_multiply(chi,(*u),psi);

  _complex_times_vector(r->s0, phase_3, chi);
  _vector_null(r->s2);
  _vector_i_sub_assign(r->s2, r->s0);

  _vector_i_sub(psi,s->s1,s->s3);

  _su3_multiply(chi,(*u),psi);

  _complex_times_vector(r->s1, phase_3, chi);
  _vector_null(r->s3);
  _vector_i_add_assign(r->s3, r->s1);

  return;
}

/* direction -z */
void boundary_D_7(spinor * const r, spinor * const s, su3 * restrict u) {

  static su3_vector chi, psi;

  _vector_i_sub(psi,s->s0,s->s2);

  _su3_inverse_multiply(chi,(*u),psi);

  _complexcjg_times_vector(r->s0, phase_3, chi);
  _vector_null(r->s2);
  _vector_i_add_assign(r->s2, r->s0);

  _vector_i_add(psi,s->s1,s->s3);

  _su3_inverse_multiply(chi,(*u),psi);

  _complexcjg_times_vector(r->s1, phase_3, chi);
  _vector_null(r->s3);
  _vector_i_sub_assign(r->s3, r->s1);

  return;
}


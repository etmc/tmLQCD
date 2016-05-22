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


#define _C_TYPE _Complex float
#define _F_TYPE float
#define _PSWITCH(s) s ## _32
#define _PTSWITCH(s) s ## 32

#include"D_psi_body.c"

#undef _C_TYPE
#undef _F_TYPE
#undef _PSWITCH
#undef _PTSWITCH


#if (!defined SSE && !defined SSE2 && !defined SSE3)

#define _C_TYPE _Complex double
#define _F_TYPE double
#define _PSWITCH(s) s
#define _PTSWITCH(s) s

#include"D_psi_body.c"

#undef _C_TYPE
#undef _F_TYPE
#undef _PSWITCH
#undef _PTSWITCH

#else

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

#ifdef TM_USE_OMP
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

#ifndef TM_USE_OMP
    iy=g_iup[0][0];
    sp=(spinor *) Q + iy;
    up=&g_gauge_field[0][0];
#endif

    /************************ loop over all lattice sites *************************/
#ifdef TM_USE_OMP
#pragma omp for
#endif
    for (ix=0;ix<VOLUME;ix++){
#ifdef TM_USE_OMP
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
#ifdef TM_USE_OMP
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

#ifdef TM_USE_OMP
#pragma omp parallel
  {
#endif
    int ix,iy,iz;
    su3 *up,*um;
    spinor *s,*sp,*sm,*rn;
    spinor ALIGN stmp,rs;

#ifndef TM_USE_OMP
    iy=g_iup[0][0];
    sp=(spinor *) Q + iy;
    up=&g_gauge_field[0][0];
#endif

    /************************ loop over all lattice sites *************************/
#ifdef TM_USE_OMP
#pragma omp for
#endif
    for (ix=0;ix<VOLUME;ix++){
#ifdef TM_USE_OMP
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
#ifdef TM_USE_OMP
  } /* OpenMP closing brace */
#endif
}

void D_psi(spinor * const P, spinor * const Q){
   if(g_c_sw > 0)
     Dsw_psi(P,Q);
   else
     Dtm_psi(P,Q);

   return ;
}

#endif // SSE

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

#define _F_TYPE float
#define _C_TYPE _Complex float
#define _PSWITCH(s) s ## _32
#define _PTSWITCH(s) s ## 32

#include "Block_D_psi_body.c"

#undef _F_TYPE
#undef _C_TYPE
#undef _PSWITCH
#undef _PTSWITCH

#define _F_TYPE double
#define _C_TYPE _Complex double
#define _PSWITCH(s) s
#define _PTSWITCH(s) s

#include "Block_D_psi_body.c"

#undef _F_TYPE
#undef _C_TYPE
#undef _PSWITCH
#undef _PTSWITCH

#ifdef TM_USE_OMP
#define static
#endif

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

#ifdef TM_USE_OMP
#undef static
#endif

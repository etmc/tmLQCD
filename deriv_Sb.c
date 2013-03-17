/***********************************************************************
 *
 * Copyright (C) 2001 Martin Hasenbusch
 *
 * some changes to initial version by Carsten Urbach
 *
 * BG version Copyright (C) 2006, 2007 Carsten Urbach
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
 * deriv_Sb: function to compute the derivative 
 * of the phi^{\dag} Q psi with respect
 * to the generators of the gauge group.
 * without the clover part.
 *
 * Author: Martin Hasenbusch <Martin.Hasenbusch@desy.de>
 * Date: Fri Oct 26 15:06:27 MEST 2001
 *
 *  both l and k are input
 *  for ieo = 0 
 *  l resides on even lattice points and k on odd lattice points
 *  for ieo = 1 
 *  l resides on odd lattice points and k on even lattice points
 *  the output is a su3adj field that is written to df0[][]
 *
 ************************************************************************/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "global.h"
#include "su3.h"
#include "boundary.h"
#include "xchange/xchange.h"
#include "sse.h"
#include "update_backward_gauge.h"
#include "hamiltonian_field.h"
#include "deriv_Sb.h"


#if (defined BGL && defined XLC)

void deriv_Sb(const int ieo, spinor * const l, spinor * const k, 
	      hamiltonian_field_t * const hf, const double factor) {

  int ix,iy, iz;
  int ioff, icx, icy, icz;
  su3 * restrict up ALIGN;
  su3 * restrict um ALIGN;
  su3adj * restrict ddd;
  static su3adj der;
  static su3 v1,v2;
  static su3_vector psia,psib,phia,phib;
  static spinor rr;
  spinor * restrict r ALIGN;
  spinor * restrict sp ALIGN;
  spinor * restrict sm ALIGN;

  /* We have 32 registers available */
  double _Complex reg00, reg01, reg02, reg03, reg04, reg05;
  double _Complex reg10, reg11, reg12, reg13, reg14, reg15;
  /* For su3 matrix, use reg00 for missing register */
  double _Complex v00, v01, v02, v10, v11, v12, v20, v21;
  /* The following contains the left spinor (12 regs) and the final  */
  /* su3 matrix to trace over */
  double _Complex r00, r01, r02, r10, r11, r12, r20, r21, r22, 
    r30, r31, r32;

#ifdef _KOJAK_INST
# pragma pomp inst begin(derivSb)
#endif

#pragma disjoint(*r, *sp, *sm, *up, *um, *ddd)
  __alignx(16, l);
  __alignx(16, k);

  if(ieo==0) {
    ioff=0;
  }
  else {
    ioff=(VOLUME+RAND)/2;
  } 

  /* for parallelization */
#ifdef MPI
  xchange_2fields(k, l, ieo);
#endif
  /************** loop over all lattice sites ****************/

  ix=g_eo2lexic[ioff];
  iy=g_iup[ix][0]; icy=g_lexic2eosub[iy];
  sp = k + icy;
  _prefetch_spinor(sp);
  up=&hf->gaugefield[ix][0];
  _prefetch_su3(up);

  for(icx = ioff; icx < (VOLUME/2+ioff); icx++){
/*     rr = (*(l + (icx-ioff))); */
    /*     rr=g_spinor_field[l][icx-ioff]; */
/*     r=&rr; */

    /* load left vector r and */
    /* multiply with gamma5   */
    r = l + (icx-ioff);
    ix=g_eo2lexic[icx];

    /*********************** direction +0 ********************/

    ddd = &hf->derivative[ix][0];
    _bgl_load_r0((*r).s0);
    _bgl_load_r1((*r).s1);
    _bgl_load_minus_r2((*r).s2);
    _bgl_load_minus_r3((*r).s3);

    _bgl_load_reg0((*sp).s0);
    _bgl_load_reg0_up((*sp).s1);
    _bgl_load_reg1((*sp).s2);
    _bgl_load_reg1_up((*sp).s3);

    _bgl_add_to_reg0_reg1();
    _bgl_add_to_reg0_up_reg1_up();

    _bgl_add_r0_to_r2_reg1();
    _bgl_add_r1_to_r3_reg1_up();

    iy=g_idn[ix][0]; icy=g_lexic2eosub[iy];
    sm = k + icy;
    _prefetch_spinor(sm);
    um=&hf->gaugefield[iy][0];
    _prefetch_su3(um);

    _bgl_tensor_product_and_add();
    /* result in v now */
    /* v is daggered as compared to non-bgl version */
    _bgl_su3_times_v(*up);
    /* result in r now */
    _bgl_complex_times_r(ka0);
    _bgl_trace_lambda_mul_add_assign((*ddd), 2.*factor);


    /************** direction -0 ****************************/

    ddd = &hf->derivative[iy][0];
    _bgl_load_r0(r->s0);
    _bgl_load_r1(r->s1);
    _bgl_load_minus_r2(r->s2);
    _bgl_load_minus_r3(r->s3);

    _bgl_load_reg0(sm->s0);
    _bgl_load_reg0_up(sm->s1);
    _bgl_load_reg1(sm->s2);
    _bgl_load_reg1_up(sm->s3);

    _bgl_sub_from_reg0_reg1();
    _bgl_sub_from_reg0_up_reg1_up();

    _bgl_sub_from_r0_r2_reg1();
    _bgl_sub_from_r1_r3_reg1_up();

    iy=g_iup[ix][1]; icy=g_lexic2eosub[iy];

    sp = k + icy;
    _prefetch_spinor(sp);
    up=&hf->gaugefield[ix][1];      
    _prefetch_su3(up);

    _bgl_tensor_product_and_add_d();
    /* result in v now */
    _bgl_su3_times_v(*um);

    /* result in r now */
    _bgl_complex_times_r(ka0);
    _bgl_trace_lambda_mul_add_assign((*ddd), 2.*factor);

    /*************** direction +1 **************************/

    ddd = &hf->derivative[ix][1];
    _bgl_load_r0(r->s0);
    _bgl_load_r1(r->s1);
    _bgl_load_minus_r2(r->s2);
    _bgl_load_minus_r3(r->s3);

    _bgl_load_reg0(sp->s0);
    _bgl_load_reg0_up(sp->s1);
    _bgl_load_reg1(sp->s2);
    _bgl_load_reg1_up(sp->s3);
    
    _bgl_i_mul_add_to_reg0_reg1_up();
    _bgl_i_mul_add_to_reg0_up_reg1();
      
    _bgl_i_mul_add_r0_to_r3_reg1();
    _bgl_i_mul_add_r1_to_r2_reg1_up();

    iy=g_idn[ix][1]; icy=g_lexic2eosub[iy];

    sm = k + icy;
    _prefetch_spinor(sm);
    um=&hf->gaugefield[iy][1];
    _prefetch_su3(um);

    _bgl_tensor_product_and_add();
    /* result in v now */
    _bgl_su3_times_v(*up);
    /* result in r now */
    _bgl_complex_times_r(ka1);
    _bgl_trace_lambda_mul_add_assign((*ddd), 2.*factor);

    /**************** direction -1 *************************/

    ddd = &hf->derivative[iy][1];
    _bgl_load_r0(r->s0);
    _bgl_load_r1(r->s1);
    _bgl_load_minus_r2(r->s2);
    _bgl_load_minus_r3(r->s3);

    _bgl_load_reg0(sm->s0);
    _bgl_load_reg0_up(sm->s1);
    _bgl_load_reg1(sm->s2);
    _bgl_load_reg1_up(sm->s3);
    
    _bgl_i_mul_sub_from_reg0_reg1_up();
    _bgl_i_mul_sub_from_reg0_up_reg1();
      
    _bgl_i_mul_sub_from_r0_r3_reg1();
    _bgl_i_mul_sub_from_r1_r2_reg1_up();

    iy=g_iup[ix][2]; icy=g_lexic2eosub[iy];

    sp = k + icy;
    _prefetch_spinor(sp);
    up=&hf->gaugefield[ix][2];
    _prefetch_su3(up);

    _bgl_tensor_product_and_add_d();
    /* result in v now */
    _bgl_su3_times_v(*um);
    /* result in r now */
    _bgl_complex_times_r(ka1);
    _bgl_trace_lambda_mul_add_assign((*ddd), 2.*factor);

    /*************** direction +2 **************************/

    ddd = &hf->derivative[ix][2];
    _bgl_load_r0(r->s0);
    _bgl_load_r1(r->s1);
    _bgl_load_minus_r2(r->s2);
    _bgl_load_minus_r3(r->s3);

    _bgl_load_reg0(sp->s0);
    _bgl_load_reg0_up(sp->s1);
    _bgl_load_reg1(sp->s2);
    _bgl_load_reg1_up(sp->s3);

    _bgl_add_to_reg0_reg1_up();
    _bgl_sub_from_reg0_up_reg1();

    _bgl_add_r0_to_r3_reg1();
    _bgl_sub_from_r1_r2_reg1_up();

    iy=g_idn[ix][2]; icy=g_lexic2eosub[iy];

    sm = k + icy;
    _prefetch_spinor(sm);
    um=&hf->gaugefield[iy][2];
    _prefetch_su3(um);

    _bgl_tensor_product_and_add();
    /* result in v now */
    _bgl_su3_times_v(*up);
    /* result in r now */
    _bgl_complex_times_r(ka2);
    _bgl_trace_lambda_mul_add_assign((*ddd), 2.*factor);
      
    /***************** direction -2 ************************/

    ddd = &hf->derivative[iy][2];
    _bgl_load_r0(r->s0);
    _bgl_load_r1(r->s1);
    _bgl_load_minus_r2(r->s2);
    _bgl_load_minus_r3(r->s3);

    _bgl_load_reg0(sm->s0);
    _bgl_load_reg0_up(sm->s1);
    _bgl_load_reg1(sm->s2);
    _bgl_load_reg1_up(sm->s3);

    _bgl_sub_from_reg0_reg1_up();
    _bgl_add_to_reg0_up_reg1();

    _bgl_sub_from_r0_r3_reg1();
    _bgl_add_r1_to_r2_reg1_up();

    iy=g_iup[ix][3]; icy=g_lexic2eosub[iy];

    sp = k + icy;
    _prefetch_spinor(sp);
    up=&hf->gaugefield[ix][3];
    _prefetch_su3(up);

    _bgl_tensor_product_and_add_d();
    /* result in v now */
    _bgl_su3_times_v(*um);
    /* result in r now */
    _bgl_complex_times_r(ka1);
    _bgl_trace_lambda_mul_add_assign(*ddd, 2.*factor);

    /****************** direction +3 ***********************/

    ddd = &hf->derivative[ix][3];
    _bgl_load_r0(r->s0);
    _bgl_load_r1(r->s1);
    _bgl_load_minus_r2(r->s2);
    _bgl_load_minus_r3(r->s3);

    _bgl_load_reg0(sp->s0);
    _bgl_load_reg0_up(sp->s1);
    _bgl_load_reg1(sp->s2);
    _bgl_load_reg1_up(sp->s3);

    _bgl_i_mul_add_to_reg0_reg1();
    _bgl_i_mul_sub_from_reg0_up_reg1_up();

    _bgl_i_mul_add_r0_to_r2_reg1();
    _bgl_i_mul_sub_from_r1_r3_reg1_up();

    iy=g_idn[ix][3]; icy=g_lexic2eosub[iy];

    sm = k + icy;
    _prefetch_spinor(sm);
    um=&hf->gaugefield[iy][3];
    _prefetch_su3(um);

    _bgl_tensor_product_and_add();
    /* result in v now */
    _bgl_su3_times_v(*up);
    /* result in r now */
    _bgl_complex_times_r(ka3);
    _bgl_trace_lambda_mul_add_assign((*ddd), 2.*factor);

    /***************** direction -3 ************************/

    ddd = &hf->derivative[iy][3];
    _bgl_load_r0(r->s0);
    _bgl_load_r1(r->s1);
    _bgl_load_minus_r2(r->s2);
    _bgl_load_minus_r3(r->s3);

    _bgl_load_reg0(sm->s0);
    _bgl_load_reg0_up(sm->s1);
    _bgl_load_reg1(sm->s2);
    _bgl_load_reg1_up(sm->s3);

    _bgl_i_mul_sub_from_reg0_reg1();
    _bgl_i_mul_add_to_reg0_up_reg1_up();

    _bgl_i_mul_sub_from_r0_r2_reg1();
    _bgl_i_mul_add_r1_to_r3_reg1_up();

    /* something wrong here...*/
    icz=icx+1;
    if(icz==((VOLUME+RAND)/2+ioff)) icz=ioff;
    iz=g_eo2lexic[icz];
    iy=g_iup[iz][0]; icy=g_lexic2eosub[iy];

    sp = k + icy;
    _prefetch_spinor(sp);
    up=&hf->gaugefield[iz][0];
    _prefetch_su3(up);

    _bgl_tensor_product_and_add_d();
    /* result in v now */
    _bgl_su3_times_v(*um);
    /* result in r now */
    _bgl_complex_times_r(ka3);
    _bgl_trace_lambda_mul_add_assign((*ddd), 2.*factor);

    /****************** end of loop ************************/
  }
#ifdef _KOJAK_INST
#pragma pomp inst end(derivSb)
#endif
}

#else

void deriv_Sb(const int ieo, spinor * const l, spinor * const k, 
	      hamiltonian_field_t * const hf, const double factor) {

#ifdef _GAUGE_COPY
  if(g_update_gauge_copy) {
    update_backward_gauge(hf->gaugefield);
  }
#endif
  /* for parallelization */
#ifdef MPI
  xchange_2fields(k, l, ieo);
#endif

#ifdef OMP
#define static
#pragma omp parallel
  {
#endif
  int ix,iy;
  int ioff, icx, icy;
  su3 * restrict up ALIGN;
  su3 * restrict um ALIGN;
  static su3 v1,v2;
  static su3_vector psia,psib,phia,phib;
  static spinor rr;
  spinor * restrict sp ALIGN;
  spinor * restrict sm ALIGN;

#ifdef OMP
#undef static
#endif

#ifdef _KOJAK_INST
#pragma pomp inst begin(derivSb)
#endif
#ifdef XLC
#pragma disjoint(*sp, *sm, *up, *um)
#endif

#ifdef BGL
  __alignx(16, l);
  __alignx(16, k);
#endif

  if(ieo==0) {
    ioff=0;
  }
  else {
    ioff=(VOLUME+RAND)/2;
  } 

  /************** loop over all lattice sites ****************/
#ifdef OMP
#pragma omp for
#endif
  for(icx = ioff; icx < (VOLUME/2+ioff); icx++){
    ix=g_eo2lexic[icx];
    rr = (*(l + (icx-ioff)));
    /*     rr=g_spinor_field[l][icx-ioff]; */

    /*multiply the left vector with gamma5*/
    _vector_minus_assign(rr.s2, rr.s2);
    _vector_minus_assign(rr.s3, rr.s3);

    /*********************** direction +0 ********************/

    iy=g_iup[ix][0]; icy=g_lexic2eosub[iy];

    sp = k + icy;
#if (defined _GAUGE_COPY && !defined _USE_HALFSPINOR && !defined  _USE_TSPLITPAR)
    up=&g_gauge_field_copy[icx][0];
#else
    up=&hf->gaugefield[ix][0];
#endif      
    _vector_add(psia,sp->s0,sp->s2);
    _vector_add(psib,sp->s1,sp->s3);
      
    _vector_add(phia,rr.s0,rr.s2);
    _vector_add(phib,rr.s1,rr.s3);

    _vector_tensor_vector_add(v1, phia, psia, phib, psib);
    _su3_times_su3d(v2,*up,v1);
    _complex_times_su3(v1, ka0, v2);
    _trace_lambda_mul_add_assign_nonlocal(hf->derivative[ix][0], 2.*factor, v1);

    /************** direction -0 ****************************/

    iy=g_idn[ix][0]; icy=g_lexic2eosub[iy];

    sm = k + icy;
#if (defined _GAUGE_COPY && !defined _USE_HALFSPINOR && !defined  _USE_TSPLITPAR)
    um = up+1;
#else
    um=&hf->gaugefield[iy][0];
#endif
      
    _vector_sub(psia,sm->s0,sm->s2);
    _vector_sub(psib,sm->s1,sm->s3);

    _vector_sub(phia,rr.s0,rr.s2);
    _vector_sub(phib,rr.s1,rr.s3);


    _vector_tensor_vector_add(v1, psia, phia, psib, phib);
    _su3_times_su3d(v2,*um,v1);
    _complex_times_su3(v1,ka0,v2);
    _trace_lambda_mul_add_assign_nonlocal(hf->derivative[iy][0], 2.*factor, v1);

    /*************** direction +1 **************************/

    iy=g_iup[ix][1]; icy=g_lexic2eosub[iy];

    sp = k + icy;
#if (defined _GAUGE_COPY && !defined _USE_HALFSPINOR && !defined  _USE_TSPLITPAR)
    up=um+1;
#else
    up=&hf->gaugefield[ix][1];      
#endif
    _vector_i_add(psia,sp->s0,sp->s3);
    _vector_i_add(psib,sp->s1,sp->s2);

    _vector_i_add(phia,rr.s0,rr.s3);
    _vector_i_add(phib,rr.s1,rr.s2);

    _vector_tensor_vector_add(v1, phia, psia, phib, psib);
    _su3_times_su3d(v2,*up,v1);
    _complex_times_su3(v1,ka1,v2);
    _trace_lambda_mul_add_assign_nonlocal(hf->derivative[ix][1], 2.*factor, v1);

    /**************** direction -1 *************************/

    iy=g_idn[ix][1]; icy=g_lexic2eosub[iy];

    sm = k + icy;
#if (defined _GAUGE_COPY && !defined _USE_HALFSPINOR && !defined  _USE_TSPLITPAR)
      um=up+1;
#else
    um=&hf->gaugefield[iy][1];
#endif
    _vector_i_sub(psia,sm->s0,sm->s3);
    _vector_i_sub(psib,sm->s1,sm->s2);

    _vector_i_sub(phia,rr.s0,rr.s3);
    _vector_i_sub(phib,rr.s1,rr.s2);

    _vector_tensor_vector_add(v1, psia, phia, psib, phib);
    _su3_times_su3d(v2,*um,v1);
    _complex_times_su3(v1,ka1,v2);
    _trace_lambda_mul_add_assign_nonlocal(hf->derivative[iy][1], 2.*factor, v1);

    /*************** direction +2 **************************/

    iy=g_iup[ix][2]; icy=g_lexic2eosub[iy];

    sp = k + icy;
#if (defined _GAUGE_COPY && !defined _USE_HALFSPINOR && !defined  _USE_TSPLITPAR)
    up=um+1;
#else
    up=&hf->gaugefield[ix][2];
#endif      
    _vector_add(psia,sp->s0,sp->s3);
    _vector_sub(psib,sp->s1,sp->s2);
      
    _vector_add(phia,rr.s0,rr.s3);
    _vector_sub(phib,rr.s1,rr.s2);

    _vector_tensor_vector_add(v1, phia, psia, phib, psib);
    _su3_times_su3d(v2,*up,v1);
    _complex_times_su3(v1,ka2,v2);
    _trace_lambda_mul_add_assign_nonlocal(hf->derivative[ix][2], 2.*factor, v1);

    /***************** direction -2 ************************/

    iy=g_idn[ix][2]; icy=g_lexic2eosub[iy];

    sm = k + icy;
#if (defined _GAUGE_COPY && !defined _USE_HALFSPINOR && !defined  _USE_TSPLITPAR)
      um = up +1;
#else
    um=&hf->gaugefield[iy][2];
#endif
    _vector_sub(psia,sm->s0,sm->s3);
    _vector_add(psib,sm->s1,sm->s2);

    _vector_sub(phia,rr.s0,rr.s3);
    _vector_add(phib,rr.s1,rr.s2);

    _vector_tensor_vector_add(v1, psia, phia, psib, phib);
    _su3_times_su3d(v2,*um,v1);
    _complex_times_su3(v1,ka2,v2);
    _trace_lambda_mul_add_assign_nonlocal(hf->derivative[iy][2], 2.*factor, v1);

    /****************** direction +3 ***********************/

    iy=g_iup[ix][3]; icy=g_lexic2eosub[iy];

    sp = k + icy;
#if (defined _GAUGE_COPY && !defined _USE_HALFSPINOR && !defined  _USE_TSPLITPAR)
    up=um+1;
#else
    up=&hf->gaugefield[ix][3];
#endif      
    _vector_i_add(psia,sp->s0,sp->s2);
    _vector_i_sub(psib,sp->s1,sp->s3);

    _vector_i_add(phia,rr.s0,rr.s2);
    _vector_i_sub(phib,rr.s1,rr.s3);

    _vector_tensor_vector_add(v1, phia, psia, phib, psib);
    _su3_times_su3d(v2,*up,v1);
    _complex_times_su3(v1, ka3, v2);
    _trace_lambda_mul_add_assign_nonlocal(hf->derivative[ix][3], 2.*factor, v1);

    /***************** direction -3 ************************/

    iy=g_idn[ix][3]; icy=g_lexic2eosub[iy];

    sm = k + icy;
#if (defined _GAUGE_COPY && !defined _USE_HALFSPINOR && !defined  _USE_TSPLITPAR)
    um = up+1;
#else
    um=&hf->gaugefield[iy][3];
#endif
    _vector_i_sub(psia,sm->s0,sm->s2);
    _vector_i_add(psib,sm->s1,sm->s3);

    _vector_i_sub(phia,rr.s0,rr.s2);
    _vector_i_add(phib,rr.s1,rr.s3);

    _vector_tensor_vector_add(v1, psia, phia, psib, phib);
    _su3_times_su3d(v2,*um,v1);
    _complex_times_su3(v1,ka3,v2);
    _trace_lambda_mul_add_assign_nonlocal(hf->derivative[iy][3], 2.*factor, v1);
     
    /****************** end of loop ************************/
  }

#ifdef OMP
  } /* OpenMP closing brace */
#endif

#ifdef _KOJAK_INST
#pragma pomp inst end(derivSb)
#endif
}

#endif


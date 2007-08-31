/* $Id$ */

/**************************************************************
 * deriv_Sb: function to compute the derivative 
 * of the phi^{\dag} Q psi with respect
 * to the generators of the gauge group.
 * without the clover part.
 *
 * Version: 0.0
 * Author: Martin Hasenbusch <Martin.Hasenbusch@desy.de>
 * Date: Fri Oct 26 15:06:27 MEST 2001
 *
 **************************************************************/
/*
  both l and k are input
  for ieo = 0 
  l resides on even lattice points and k on odd lattice points
  for ieo = 1 
  l resides on odd lattice points and k on even lattice points
  the output is a su3adj field that is written to df0[][]
*/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "global.h"
#include "su3.h"
#include "boundary.h"
#include "xchange_field.h"
#include "sse.h"
#include "deriv_Sb.h"


#if (defined BGLnotchecked && defined XLC)

void deriv_Sb(const int ieo, spinor * const l, spinor * const k) {

  int ix,iy, iz;
  int ioff,ioff2,icx,icy, icz;
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
  ioff2=(VOLUME+RAND)/2-ioff;

  /* for parallelization */
#ifdef MPI
  xchange_field(k, ieo);
  xchange_field(l, (ieo+1)%2);
#endif
  /************** loop over all lattice sites ****************/

  ix=g_eo2lexic[ioff];
  iy=g_iup[ix][0]; icy=g_lexic2eosub[iy];
  sp = k + icy;
  _prefetch_spinor(sp);
  up=&g_gauge_field[ix][0];
  _prefetch_su3(up);

  for(icx = ioff; icx < (VOLUME/2+ioff); icx++){
/*     rr = (*(l + (icx-ioff))); */
    /*     rr=g_spinor_field[l][icx-ioff]; */
/*     r=&rr; */

    /* load left vector r and */
    /* multiply with gamma5   */
    r = l + (icx-ioff);


    /*********************** direction +0 ********************/

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
    um=&g_gauge_field[iy][0];
    _prefetch_su3(um);

    _bgl_tensor_product_and_add();
    /* result in v now */
    _bgl_su3_times_v_dagger(*up);
    /* result in r now */
    _bgl_complex_times_r(ka0);
    _bgl_trace_lambda_add_assign(df0[ix][0]);

    /************** direction -0 ****************************/

    _bgl_load_r0((*r).s0);
    _bgl_load_r1((*r).s1);
    _bgl_load_minus_r2((*r).s2);
    _bgl_load_minus_r3((*r).s3);

    _bgl_load_reg0((*sm).s0);
    _bgl_load_reg0_up((*sm).s1);
    _bgl_load_reg1((*sm).s2);
    _bgl_load_reg1_up((*sm).s3);

    _bgl_sub_from_reg0_reg1();
    _bgl_sub_from_reg0_up_reg1_up();

    _bgl_sub_from_r0_r2_reg1();
    _bgl_sub_from_r1_r3_reg1_up();

    iy=g_iup[ix][1]; icy=g_lexic2eosub[iy];

    sp = k + icy;
    _prefetch_spinor(sp);
    up=&g_gauge_field[ix][1];      
    _prefetch_su3(up);

    _bgl_tensor_product_and_add();
    /* result in v now */
    _bgl_su3_times_v_dagger(*um);
    /* result in r now */
    _bgl_complex_times_r(ka0);
    _bgl_trace_lambda_add_assign(df0[iy][0]);

    /*************** direction +1 **************************/

    _bgl_load_r0((*r).s0);
    _bgl_load_r1((*r).s1);
    _bgl_load_minus_r2((*r).s2);
    _bgl_load_minus_r3((*r).s3);

    _bgl_load_reg0((*sp).s0);
    _bgl_load_reg0_up((*sp).s1);
    _bgl_load_reg1((*sp).s2);
    _bgl_load_reg1_up((*sp).s3);
    
    _bgl_i_mul_add_to_reg0_reg1_up();
    _bgl_i_mul_add_to_reg0_up_reg1();
      
    _bgl_i_mul_add_r0_to_r3_reg1();
    _bgl_i_mul_add_r1_to_r2_reg1_up();

    iy=g_idn[ix][1]; icy=g_lexic2eosub[iy];

    sm = k + icy;
    _prefetch_spinor(sm);
    um=&g_gauge_field[iy][1];
    _prefetch_su3(um);

    _bgl_tensor_product_and_add();
    /* result in v now */
    _bgl_su3_times_v_dagger(*up);
    /* result in r now */
    _bgl_complex_times_r(ka1);
    _bgl_trace_lambda_add_assign(df0[ix][1]);

    /**************** direction -1 *************************/

    _bgl_load_r0((*r).s0);
    _bgl_load_r1((*r).s1);
    _bgl_load_minus_r2((*r).s2);
    _bgl_load_minus_r3((*r).s3);

    _bgl_load_reg0((*sp).s0);
    _bgl_load_reg0_up((*sp).s1);
    _bgl_load_reg1((*sp).s2);
    _bgl_load_reg1_up((*sp).s3);
    
    _bgl_i_mul_sub_from_reg0_reg1_up();
    _bgl_i_mul_sub_from_reg0_up_reg1();
      
    _bgl_i_mul_sub_from_r0_r3_reg1();
    _bgl_i_mul_sub_from_r1_r2_reg1_up();

    iy=g_iup[ix][2]; icy=g_lexic2eosub[iy];

    sp = k + icy;
    _prefetch_spinor(sp);
    up=&g_gauge_field[ix][2];
    _prefetch_su3(up);

    _bgl_tensor_product_and_add();
    /* result in v now */
    _bgl_su3_times_v_dagger(*um);
    /* result in r now */
    _bgl_complex_times_r(ka1);
    _bgl_trace_lambda_add_assign(df0[iy][1]);

    /*************** direction +2 **************************/

    _bgl_load_r0((*r).s0);
    _bgl_load_r1((*r).s1);
    _bgl_load_minus_r2((*r).s2);
    _bgl_load_minus_r3((*r).s3);

    _bgl_load_reg0((*sp).s0);
    _bgl_load_reg0_up((*sp).s1);
    _bgl_load_reg1((*sp).s2);
    _bgl_load_reg1_up((*sp).s3);
    
    _bgl_add_to_reg0_reg1_up();
    _bgl_sub_from_reg0_up_reg1();
      
    _bgl_add_r0_to_r3_reg1();
    _bgl_sub_from_r1_r2_reg1_up();

    iy=g_idn[ix][2]; icy=g_lexic2eosub[iy];

    sm = k + icy;
    _prefetch_spinor(sm);
    um=&g_gauge_field[iy][2];
    _prefetch_su3(um);

    _bgl_tensor_product_and_add();
    /* result in v now */
    _bgl_su3_times_v_dagger(*up);
    /* result in r now */
    _bgl_complex_times_r(ka2);
    _bgl_trace_lambda_add_assign(df0[ix][2]);
      
    /***************** direction -2 ************************/

    _bgl_load_r0((*r).s0);
    _bgl_load_r1((*r).s1);
    _bgl_load_minus_r2((*r).s2);
    _bgl_load_minus_r3((*r).s3);

    _bgl_load_reg0((*sp).s0);
    _bgl_load_reg0_up((*sp).s1);
    _bgl_load_reg1((*sp).s2);
    _bgl_load_reg1_up((*sp).s3);
    
    _bgl_sub_from_reg0_reg1_up();
    _bgl_add_to_reg0_up_reg1();
      
    _bgl_sub_from_r0_r3_reg1();
    _bgl_add_r1_to_r2_reg1_up();

    iy=g_iup[ix][3]; icy=g_lexic2eosub[iy];

    sp = k + icy;
    _prefetch_spinor(sp);
    up=&g_gauge_field[ix][3];
    _prefetch_su3(up);

    _bgl_tensor_product_and_add();
    /* result in v now */
    _bgl_su3_times_v_dagger(*um);
    /* result in r now */
    _bgl_complex_times_r(ka1);
    _bgl_trace_lambda_add_assign(df0[iy][2]);

    /****************** direction +3 ***********************/

    _bgl_load_r0((*r).s0);
    _bgl_load_r1((*r).s1);
    _bgl_load_minus_r2((*r).s2);
    _bgl_load_minus_r3((*r).s3);

    _bgl_load_reg0((*sp).s0);
    _bgl_load_reg0_up((*sp).s1);
    _bgl_load_reg1((*sp).s2);
    _bgl_load_reg1_up((*sp).s3);
    
    _bgl_i_mul_add_to_reg0_reg1();
    _bgl_i_mul_sub_from_reg0_up_reg1_up();
      
    _bgl_i_mul_add_r0_to_r2_reg1();
    _bgl_i_mul_sub_from_r1_r3_reg1_up();

    iy=g_idn[ix][3]; icy=g_lexic2eosub[iy];

    sm = k + icy;
    _prefetch_spinor(sm);
    um=&g_gauge_field[iy][3];
    _prefetch_su3(um);

    _bgl_tensor_product_and_add();
    /* result in v now */
    _bgl_su3_times_v_dagger(*up);
    /* result in r now */
    _bgl_complex_times_r(ka3);
    _bgl_trace_lambda_add_assign(df0[ix][3]);

    /***************** direction -3 ************************/

    _bgl_load_r0((*r).s0);
    _bgl_load_r1((*r).s1);
    _bgl_load_minus_r2((*r).s2);
    _bgl_load_minus_r3((*r).s3);

    _bgl_load_reg0((*sp).s0);
    _bgl_load_reg0_up((*sp).s1);
    _bgl_load_reg1((*sp).s2);
    _bgl_load_reg1_up((*sp).s3);
    
    _bgl_i_mul_sub_from_reg0_reg1();
    _bgl_i_mul_add_to_reg0_up_reg1_up();
      
    _bgl_i_mul_sub_from_r0_r2_reg1();
    _bgl_i_mul_add_r1_to_r3_reg1_up();

    icz=icx+1;
    if(icz==((VOLUME+RAND)/2+ioff)) icz=ioff;
    iz=g_eo2lexic[icz];
    iy=g_iup[iz][0]; icy=g_lexic2eosub[iy];

    sp = k + icy;
    _prefetch_spinor(sp);
    up=&g_gauge_field[iz][0];
    _prefetch_su3(up);

    _bgl_tensor_product_and_add();
    /* result in v now */
    _bgl_su3_times_v_dagger(*um);
    /* result in r now */
    _bgl_complex_times_r(ka3);
    _bgl_trace_lambda_add_assign(df0[iy][3]);

    /****************** end of loop ************************/
  }
#ifdef _KOJAK_INST
#pragma pomp inst end(derivSb)
#endif
}

#else



void deriv_Sb(const int ieo, spinor * const l, spinor * const k) {
/* const int l, const int k){ */
  int ix,iy;
  int ioff,ioff2,icx,icy;
  su3 * restrict up ALIGN;
  su3 * restrict um ALIGN;
/*   su3adj * restrict ddd; */
/*   static su3adj der; */
  static su3 v1,v2;
  static su3_vector psia,psib,phia,phib;
  static spinor rr;
/*   spinor * restrict r ALIGN; */
  spinor * restrict sp ALIGN;
  spinor * restrict sm ALIGN;

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
  ioff2=(VOLUME+RAND)/2-ioff;

  /* for parallelization */
#ifdef MPI
  xchange_field(k, ieo);
  xchange_field(l, (ieo+1)%2);
#endif
  /************** loop over all lattice sites ****************/

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
/*     sp=&g_spinor_field[k][icy]; */
    up=&g_gauge_field[ix][0];
      
    _vector_add(psia,(*sp).s0,(*sp).s2);
    _vector_add(psib,(*sp).s1,(*sp).s3);
      
    _vector_add(phia,rr.s0,rr.s2);
    _vector_add(phib,rr.s1,rr.s3);

    _vector_tensor_vector_add(v1, phia, psia, phib, psib);
/*     _vector_tensor_vector(v1,phia,psia); */
/*     _vector_tensor_vector(v2,phib,psib); */
/*     _su3_plus_su3(v1,v1,v2); */
    _su3_times_su3d(v2,*up,v1);
    _complex_times_su3(v1,ka0,v2);
    _trace_lambda_add_assign(df0[ix][0], v1);
/*     _trace_lambda(der,v1); */
/*     ddd=&df0[ix][0]; */
/*     _add_su3adj(*ddd,der); */

    /************** direction -0 ****************************/

    iy=g_idn[ix][0]; icy=g_lexic2eosub[iy];

    sm = k + icy;
/*     sm=&g_spinor_field[k][icy]; */
    um=&g_gauge_field[iy][0];
      
    _vector_sub(psia,(*sm).s0,(*sm).s2);
    _vector_sub(psib,(*sm).s1,(*sm).s3);

    _vector_sub(phia,rr.s0,rr.s2);
    _vector_sub(phib,rr.s1,rr.s3);

    _vector_tensor_vector_add(v1, psia, phia, psib, phib);
/*     _vector_tensor_vector(v1,psia,phia); */
/*     _vector_tensor_vector(v2,psib,phib); */
/*     _su3_plus_su3(v1,v1,v2); */
    _su3_times_su3d(v2,*um,v1);
    _complex_times_su3(v1,ka0,v2);
    _trace_lambda_add_assign(df0[iy][0], v1);
/*     _trace_lambda(der,v1); */
/*     ddd=&df0[iy][0]; */
/*     _add_su3adj(*ddd,der); */

    /*************** direction +1 **************************/

    iy=g_iup[ix][1]; icy=g_lexic2eosub[iy];

    sp = k + icy;
    /*     sp=&g_spinor_field[k][icy]; */
    up=&g_gauge_field[ix][1];      

    _vector_i_add(psia,(*sp).s0,(*sp).s3);
    _vector_i_add(psib,(*sp).s1,(*sp).s2);

    _vector_i_add(phia,rr.s0,rr.s3);
    _vector_i_add(phib,rr.s1,rr.s2);

    _vector_tensor_vector_add(v1, phia, psia, phib, psib);
/*     _vector_tensor_vector(v1,phia,psia); */
/*     _vector_tensor_vector(v2,phib,psib); */
/*     _su3_plus_su3(v1,v1,v2); */
    _su3_times_su3d(v2,*up,v1);
    _complex_times_su3(v1,ka1,v2);
    _trace_lambda_add_assign(df0[ix][1], v1);
/*     _trace_lambda(der,v1); */
/*     ddd=&df0[ix][1]; */
/*     _add_su3adj(*ddd,der); */

    /**************** direction -1 *************************/

    iy=g_idn[ix][1]; icy=g_lexic2eosub[iy];

    sm = k + icy;
    /*     sm=&g_spinor_field[k][icy]; */
    um=&g_gauge_field[iy][1];
      
    _vector_i_sub(psia,(*sm).s0,(*sm).s3);
    _vector_i_sub(psib,(*sm).s1,(*sm).s2);

    _vector_i_sub(phia,rr.s0,rr.s3);
    _vector_i_sub(phib,rr.s1,rr.s2);

    _vector_tensor_vector_add(v1, psia, phia, psib, phib);
/*     _vector_tensor_vector(v1,psia,phia); */
/*     _vector_tensor_vector(v2,psib,phib); */
/*     _su3_plus_su3(v1,v1,v2); */
    _su3_times_su3d(v2,*um,v1);
    _complex_times_su3(v1,ka1,v2);
    _trace_lambda_add_assign(df0[iy][1], v1);
/*     _trace_lambda(der,v1); */
/*     ddd=&df0[iy][1]; */
/*     _add_su3adj(*ddd,der); */

    /*************** direction +2 **************************/

    iy=g_iup[ix][2]; icy=g_lexic2eosub[iy];

    sp = k + icy;
    /*     sp=&g_spinor_field[k][icy]; */
    up=&g_gauge_field[ix][2];
      
    _vector_add(psia,(*sp).s0,(*sp).s3);
    _vector_sub(psib,(*sp).s1,(*sp).s2);
      
    _vector_add(phia,rr.s0,rr.s3);
    _vector_sub(phib,rr.s1,rr.s2);

    _vector_tensor_vector_add(v1, phia, psia, phib, psib);
/*     _vector_tensor_vector(v1,phia,psia); */
/*     _vector_tensor_vector(v2,phib,psib); */
/*     _su3_plus_su3(v1,v1,v2); */
    _su3_times_su3d(v2,*up,v1);
    _complex_times_su3(v1,ka2,v2);
    _trace_lambda_add_assign(df0[ix][2], v1);
/*     _trace_lambda(der,v1); */
/*     ddd=&df0[ix][2]; */
/*     _add_su3adj(*ddd,der); */

    /***************** direction -2 ************************/

    iy=g_idn[ix][2]; icy=g_lexic2eosub[iy];

    sm = k + icy;
    /*     sm=&g_spinor_field[k][icy]; */
    um=&g_gauge_field[iy][2];
      
    _vector_sub(psia,(*sm).s0,(*sm).s3);
    _vector_add(psib,(*sm).s1,(*sm).s2);

    _vector_sub(phia,rr.s0,rr.s3);
    _vector_add(phib,rr.s1,rr.s2);

    _vector_tensor_vector_add(v1, psia, phia, psib, phib);
/*     _vector_tensor_vector(v1,psia,phia); */
/*     _vector_tensor_vector(v2,psib,phib); */
/*     _su3_plus_su3(v1,v1,v2); */
    _su3_times_su3d(v2,*um,v1);
    _complex_times_su3(v1,ka2,v2);
    _trace_lambda_add_assign(df0[iy][2], v1);
/*     _trace_lambda(der,v1); */
/*     ddd=&df0[iy][2]; */
/*     _add_su3adj(*ddd,der); */

    /****************** direction +3 ***********************/

    iy=g_iup[ix][3]; icy=g_lexic2eosub[iy];

    sp = k + icy;
    /*     sp=&g_spinor_field[k][icy]; */
    up=&g_gauge_field[ix][3];
      
    _vector_i_add(psia,(*sp).s0,(*sp).s2);
    _vector_i_sub(psib,(*sp).s1,(*sp).s3);

    _vector_i_add(phia,rr.s0,rr.s2);
    _vector_i_sub(phib,rr.s1,rr.s3);

    _vector_tensor_vector_add(v1, phia, psia, phib, psib);
/*     _vector_tensor_vector(v1,phia,psia); */
/*     _vector_tensor_vector(v2,phib,psib); */
/*     _su3_plus_su3(v1,v1,v2); */
    _su3_times_su3d(v2,*up,v1);
    _complex_times_su3(v1,ka3,v2);
    _trace_lambda_add_assign(df0[ix][3], v1);
/*     _trace_lambda(der,v1); */
/*     ddd=&df0[ix][3]; */
/*     _add_su3adj(*ddd,der); */

    /***************** direction -3 ************************/

    iy=g_idn[ix][3]; icy=g_lexic2eosub[iy];

    sm = k + icy;
    /*     sm=&g_spinor_field[k][icy]; */
    um=&g_gauge_field[iy][3];
      
    _vector_i_sub(psia,(*sm).s0,(*sm).s2);
    _vector_i_add(psib,(*sm).s1,(*sm).s3);

    _vector_i_sub(phia,rr.s0,rr.s2);
    _vector_i_add(phib,rr.s1,rr.s3);

    _vector_tensor_vector_add(v1, psia, phia, psib, phib);
/*     _vector_tensor_vector(v1,psia,phia); */
/*     _vector_tensor_vector(v2,psib,phib); */
/*     _su3_plus_su3(v1,v1,v2); */
    _su3_times_su3d(v2,*um,v1);
    _complex_times_su3(v1,ka3,v2);
    _trace_lambda_add_assign(df0[iy][3], v1);
/*     _trace_lambda(der,v1); */
/*     ddd=&df0[iy][3]; */
/*     _add_su3adj(*ddd,der); */
     
    /****************** end of loop ************************/
  }
#ifdef _KOJAK_INST
#pragma pomp inst end(derivSb)
#endif
}

#endif

static char const rcsid[] = "$Id$";

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

#include <stdlib.h>
#include <stdio.h>
#include "global.h"
#include "su3.h"
#include "boundary.h"
#include "xchange.h"
#include "deriv_Sb.h"

static su3_vector psi1, psi2, psi, chi, phi1, phi3;

void deriv_Sb(const int ieo, const int l, const int k){
  int ix,iy;
  int ioff,ioff2,icx,icy;
  su3 *up,*um;
  su3adj *ddd;
  static su3adj der;
  static su3 v1,v2;
  static su3_vector psia,psib,phia,phib;
  static spinor rr;
  spinor *r,*sp,*sm;
  if(ieo==0){
    ioff=0;
  } 
  else{
    ioff=(VOLUME+RAND)/2;
  } 
  ioff2=(VOLUME+RAND)/2-ioff;

  /* for parallelization */
#ifdef MPI
  xchange_field(k);
  xchange_field(l);
#endif
  /************** loop over all lattice sites ****************/

  for(icx = ioff; icx < (VOLUME/2+ioff); icx++){
    ix=trans2[icx];
    rr=spinor_field[l][icx-ioff];
    r=&rr;

    /*multiply the left vector with gamma5*/
    _vector_minus_assign((*r).c3, (*r).c3);
    _vector_minus_assign((*r).c4, (*r).c4);

    /*********************** direction +0 ********************/

    iy=g_iup[ix][0]; icy=trans1[iy]-ioff2;


    sp=&spinor_field[k][icy];
    up=&g_gauge_field[ix][0];
      
    _vector_add(psia,(*sp).c1,(*sp).c3);
    _vector_add(psib,(*sp).c2,(*sp).c4);
      
    _vector_add(phia,(*r).c1,(*r).c3);
    _vector_add(phib,(*r).c2,(*r).c4);

    _vector_tensor_vector(v1,phia,psia);
    _vector_tensor_vector(v2,phib,psib);
    _su3_plus_su3(v1,v1,v2);
    _su3_times_su3d(v2,*up,v1);
    _complex_times_su3(v1,ka0,v2);
    _trace_lambda(der,v1);
    ddd=&df0[ix][0];
    _add_su3adj(*ddd,der);
    /************** direction -0 ****************************/

    iy=g_idn[ix][0]; icy=trans1[iy]-ioff2;

    sm=&spinor_field[k][icy];
    um=&g_gauge_field[iy][0];
      
    _vector_sub(psia,(*sm).c1,(*sm).c3);
    _vector_sub(psib,(*sm).c2,(*sm).c4);

    _vector_sub(phia,(*r).c1,(*r).c3);
    _vector_sub(phib,(*r).c2,(*r).c4);

    _vector_tensor_vector(v1,psia,phia);
    _vector_tensor_vector(v2,psib,phib);
    _su3_plus_su3(v1,v1,v2);
    _su3_times_su3d(v2,*um,v1);
    _complex_times_su3(v1,ka0,v2);
    _trace_lambda(der,v1);
    ddd=&df0[iy][0];
    _add_su3adj(*ddd,der);
    /*************** direction +1 **************************/

    iy=g_iup[ix][1]; icy=trans1[iy]-ioff2;

    sp=&spinor_field[k][icy];
    up=&g_gauge_field[ix][1];      

    _vector_i_add(psia,(*sp).c1,(*sp).c4);
    _vector_i_add(psib,(*sp).c2,(*sp).c3);

    _vector_i_add(phia,(*r).c1,(*r).c4);
    _vector_i_add(phib,(*r).c2,(*r).c3);

    _vector_tensor_vector(v1,phia,psia);
    _vector_tensor_vector(v2,phib,psib);
    _su3_plus_su3(v1,v1,v2);
    _su3_times_su3d(v2,*up,v1);
    _complex_times_su3(v1,ka1,v2);
    _trace_lambda(der,v1);
    ddd=&df0[ix][1];
    _add_su3adj(*ddd,der);
    /**************** direction -1 *************************/

    iy=g_idn[ix][1]; icy=trans1[iy]-ioff2;

    sm=&spinor_field[k][icy];
    um=&g_gauge_field[iy][1];
      
    _vector_i_sub(psia,(*sm).c1,(*sm).c4);
    _vector_i_sub(psib,(*sm).c2,(*sm).c3);

    _vector_i_sub(phia,(*r).c1,(*r).c4);
    _vector_i_sub(phib,(*r).c2,(*r).c3);

    _vector_tensor_vector(v1,psia,phia);
    _vector_tensor_vector(v2,psib,phib);
    _su3_plus_su3(v1,v1,v2);
    _su3_times_su3d(v2,*um,v1);
    _complex_times_su3(v1,ka1,v2);
    _trace_lambda(der,v1);
    ddd=&df0[iy][1];
    _add_su3adj(*ddd,der);

    /*************** direction +2 **************************/

    iy=g_iup[ix][2]; icy=trans1[iy]-ioff2;

    sp=&spinor_field[k][icy];
    up=&g_gauge_field[ix][2];
      
    _vector_add(psia,(*sp).c1,(*sp).c4);
    _vector_sub(psib,(*sp).c2,(*sp).c3);
      
    _vector_add(phia,(*r).c1,(*r).c4);
    _vector_sub(phib,(*r).c2,(*r).c3);

    _vector_tensor_vector(v1,phia,psia);
    _vector_tensor_vector(v2,phib,psib);
    _su3_plus_su3(v1,v1,v2);
    _su3_times_su3d(v2,*up,v1);
    _complex_times_su3(v1,ka2,v2);
    _trace_lambda(der,v1);
    ddd=&df0[ix][2];
    _add_su3adj(*ddd,der);

    /***************** direction -2 ************************/

    iy=g_idn[ix][2]; icy=trans1[iy]-ioff2;

    sm=&spinor_field[k][icy];
    um=&g_gauge_field[iy][2];
      
    _vector_sub(psia,(*sm).c1,(*sm).c4);
    _vector_add(psib,(*sm).c2,(*sm).c3);

    _vector_sub(phia,(*r).c1,(*r).c4);
    _vector_add(phib,(*r).c2,(*r).c3);

    _vector_tensor_vector(v1,psia,phia);
    _vector_tensor_vector(v2,psib,phib);

    _su3_plus_su3(v1,v1,v2);
    _su3_times_su3d(v2,*um,v1);
    _complex_times_su3(v1,ka2,v2);
    _trace_lambda(der,v1);
    ddd=&df0[iy][2];
    _add_su3adj(*ddd,der);

    /****************** direction +3 ***********************/

    iy=g_iup[ix][3]; icy=trans1[iy]-ioff2;

    sp=&spinor_field[k][icy];
    up=&g_gauge_field[ix][3];
      
    _vector_i_add(psia,(*sp).c1,(*sp).c3);
    _vector_i_sub(psib,(*sp).c2,(*sp).c4);

    _vector_i_add(phia,(*r).c1,(*r).c3);
    _vector_i_sub(phib,(*r).c2,(*r).c4);

    _vector_tensor_vector(v1,phia,psia);
    _vector_tensor_vector(v2,phib,psib);

    _su3_plus_su3(v1,v1,v2);
    _su3_times_su3d(v2,*up,v1);
    _complex_times_su3(v1,ka3,v2);
    _trace_lambda(der,v1);
    ddd=&df0[ix][3];
    _add_su3adj(*ddd,der);

    /***************** direction -3 ************************/

    iy=g_idn[ix][3]; icy=trans1[iy]-ioff2;

    sm=&spinor_field[k][icy];
    um=&g_gauge_field[iy][3];
      
    _vector_i_sub(psia,(*sm).c1,(*sm).c3);
    _vector_i_add(psib,(*sm).c2,(*sm).c4);

    _vector_i_sub(phia,(*r).c1,(*r).c3);
    _vector_i_add(phib,(*r).c2,(*r).c4);

    _vector_tensor_vector(v1,psia,phia);
    _vector_tensor_vector(v2,psib,phib);
    _su3_plus_su3(v1,v1,v2);
    _su3_times_su3d(v2,*um,v1);
    _complex_times_su3(v1,ka3,v2);
    _trace_lambda(der,v1);
    ddd=&df0[iy][3];
    _add_su3adj(*ddd,der);
     
    /****************** end of loop ************************/
  }
}

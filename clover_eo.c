#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "sse.h"
#include "su3.h"
#include "su3adj.h"
#include "global.h"
#include "xchange.h"
#include "clover_eo.h"

static complex ka0, ka1, ka2, ka3;
static su3_vector psi1, psi2, psi, chi, phi1, phi3;
const double PI_ = 3.14159265358979;

/* declaration of internal functions */

/******************************************
 * H_eo is the conventional Wilson 
 * hopping matrix
 *
 * for ieo == 1 H_eo connects even with odd
 * sides, and for ieo == 0 H_eo connects
 * odd with even sides.
 *
 * l is the number of the output field
 * k is the number of the input field
 *
 ******************************************/

void H_eo(int ieo, int l, int k);

/******************************************
 * clover_inv is (1+T_{ee|oo})^{-1}
 *
 * for ieo == 1 clover_inv is for even even
 * for ieo == 0 it is for odd odd
 *
 * l is number of input and output field
 *
 ******************************************/

void clover_inv(int ieo, int l);

/******************************************
 * clover computes
 * l = ((q_off*1+T_{ee|oo})*k-j)
 *
 * for ieo == 0 clover is for T_{oo}
 *     but uses only even sides of l
 * for iwo == 0 it is for T_{ee}
 *     and uses only odd sides of l
 *
 * l is number of the output field
 * k and j the numbers of the input fields
 *
 ******************************************/

void clover(int ieo, int l, int k, int j, double q_off);

/******************************************
 * clover_gamma5 computes
 * l = gamma_5*((q_off*1+T_{ee|oo})*k-j)
 *
 * l is the number of the output field
 * k and j the numbers of the input fields
 *
 * for ieo == 0 it is T_{oo} but the 
 *     but the it works only on even sides
 *     of the output field l
 * for ieo == 1 it is T_{ee}
 *     and it works on odd sides of l.
 *
 ******************************************/
void clover_gamma5(int ieo, int l, int k, int j, double q_off);

/******************************************
 * gamma5 multiplies k with \gamma_5
 * and stores the result in l
 * 
 * l is the number of the output field
 * k is the number of the input field 
 *
 ******************************************/

void gamma5(int l,int k);

/* definition of external accessable functions */

/***********************************************
 * This function defines the boundary cond.
 * 
 * Here is not really the place for this
 * function, but ka* are needed in this
 * file
 *
 ***********************************************/

void boundary(){
  double x0,x1,x2,x3;
  /* anti-periodic in time */
  x0 = X0 * PI_/((T)*g_nproc);
  x1 = X1 * PI_/(L);
  x2 = X2 * PI_/(L);
  x3 = X3 * PI_/(L);
  ka0.re = g_kappa * cos(x0); 
  ka0.im = g_kappa * sin(x0);
  ka1.re = g_kappa * cos(x1); 
  ka1.im = g_kappa * sin(x1);
  ka2.re = g_kappa * cos(x2); 
  ka2.im = g_kappa * sin(x2);
  ka3.re = g_kappa * cos(x3); 
  ka3.im = g_kappa * sin(x3);
}

/******************************************
 * Q_psi is the same as \hat Q
 * applied to a field \psi (see hep-lat/9603008)
 *
 * \hat Q = \gamma_5(q_off*1+T_{ee}
 *          - M_{oe}(1+T_{ee})^(-1)M_{eo})
 *
 * l is the number of the output field 
 * k is the number of the input field
 * q_off
 *
 ******************************************/

void Q_psi(int l, int k, double q_off){
  H_eo(1, DUM_MATRIX+1, k);
  clover_inv(1, DUM_MATRIX+1);
  H_eo(0, DUM_MATRIX, DUM_MATRIX+1);
  clover_gamma5(0, l, k, DUM_MATRIX, q_off);
}

/******************************************
 * This is the fermion matrix:
 * \gamma_5 \hat Q or simply
 * \hat Q without \gamma_5
 *
 * for input output see Q_psi
 *
 ******************************************/

void M_psi(int l, int k, double q_off){
  H_eo(1, DUM_MATRIX+1, k);
  clover_inv(1, DUM_MATRIX+1);
  H_eo(0, DUM_MATRIX, DUM_MATRIX+1);
  clover(0, l, k, DUM_MATRIX, q_off);
}

void H_eo_psi(int ieo, int l, int k){
  H_eo(ieo, l, k);
  clover_inv(ieo, l);
}

/**************************************************************
 * derivf: function to compute the derivative 
 * of the phi^{\dag} Q psi with respect
 * to the generators of the gauge group.
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

void derivf(int ieo, int l,int k){
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

/* definitions of internal functions */

#if defined SSE2

/* input on k; output on l */
void H_eo(int ieo, int l, int k){
  int icx,icy,icz,ioff,ioff2;
  int ix,iy,iz;
  su3 *up,*um;
  spinor *sp,*sm,*rn;
  static spinor rs;

  /* for parallelization */
#ifdef MPI
  xchange_field(k);
#endif

  if(k == l){
    printf("Error in subroutine D_psi: improper arguments\n");
    printf("Program aborted\n");
    exit(1);
  }
  if(ieo == 0){
    ioff = 0;
  } 
  else{
    ioff = (VOLUME+RAND)/2;
  }
  ioff2 = (VOLUME+RAND)/2-ioff;

  ix=trans2[ioff];
  iy=g_iup[ix][0]; 
  icy=trans1[iy]-ioff2;

  sp=&spinor_field[k][icy];
  up=&g_gauge_field[ix][0];
   
  /**************** loop over all lattice sites ******************/
  for(icx = ioff; icx < (VOLUME/2+ioff); icx++){
    ix=trans2[icx];
    /*********************** direction +0 ************************/

    iy=g_idn[ix][0]; icy=trans1[iy]-ioff2;

    um=&g_gauge_field[iy][0];
    _prefetch_su3(um);
      
    sm=&spinor_field[k][icy];
    _prefetch_spinor(sm);

    _sse_load((*sp).c1);
    _sse_load_up((*sp).c3);
    _sse_vector_add();

    _sse_su3_multiply((*up));
    _sse_vector_cmplx_mul(ka0);
    _sse_store_up(rs.c1);
    _sse_store_up(rs.c3);      
      
    _sse_load((*sp).c2);
    _sse_load_up((*sp).c4);
    _sse_vector_add();
      
    _sse_su3_multiply((*up));
    _sse_vector_cmplx_mul(ka0);
    _sse_store_up(rs.c2);
    _sse_store_up(rs.c4); 

    /*********************** direction -0 ************************/

    iy=g_iup[ix][1]; icy=trans1[iy]-ioff2;

    up+=1;
    _prefetch_su3(up);
      
    sp=&spinor_field[k][icy];
    _prefetch_spinor(sp);

    _sse_load((*sm).c1);
    _sse_load_up((*sm).c3);
    _sse_vector_sub();
      
    _sse_su3_inverse_multiply((*um));
    _sse_vector_cmplxcg_mul(ka0);
      
    _sse_load(rs.c1);
    _sse_vector_add();
    _sse_store(rs.c1);

    _sse_load(rs.c3);
    _sse_vector_sub();
    _sse_store(rs.c3);
      
    _sse_load((*sm).c2);
    _sse_load_up((*sm).c4);
    _sse_vector_sub();
      
    _sse_su3_inverse_multiply((*um));
    _sse_vector_cmplxcg_mul(ka0);
      
    _sse_load(rs.c2);
    _sse_vector_add();
    _sse_store(rs.c2);

    _sse_load(rs.c4);
    _sse_vector_sub();
    _sse_store(rs.c4);
      
    /*********************** direction +1 ************************/

    iy=g_idn[ix][1]; icy=trans1[iy]-ioff2;

    um=&g_gauge_field[iy][1];
    _prefetch_su3(um);

    sm=&spinor_field[k][icy];
    _prefetch_spinor(sm);

    _sse_load((*sp).c1);
    _sse_load_up((*sp).c4);
    _sse_vector_i_mul();
    _sse_vector_add();

    _sse_su3_multiply((*up));
    _sse_vector_cmplx_mul(ka1);

    _sse_load(rs.c1);
    _sse_vector_add();
    _sse_store(rs.c1);

    _sse_load(rs.c4);
    _sse_vector_i_mul();      
    _sse_vector_sub();
    _sse_store(rs.c4); 
      
    _sse_load((*sp).c2);
    _sse_load_up((*sp).c3);
    _sse_vector_i_mul();
    _sse_vector_add();

    _sse_su3_multiply((*up));
    _sse_vector_cmplx_mul(ka1);

    _sse_load(rs.c2);
    _sse_vector_add();
    _sse_store(rs.c2);

    _sse_load(rs.c3);
    _sse_vector_i_mul();      
    _sse_vector_sub();
    _sse_store(rs.c3);       

    /*********************** direction -1 ************************/

    iy=g_iup[ix][2]; icy=trans1[iy]-ioff2;

    up+=1;
    _prefetch_su3(up);

    sp=&spinor_field[k][icy];
    _prefetch_spinor(sp);

    _sse_load((*sm).c1);
    _sse_load_up((*sm).c4);
    _sse_vector_i_mul();
    _sse_vector_sub();
      
    _sse_su3_inverse_multiply((*um));
    _sse_vector_cmplxcg_mul(ka1);
      
    _sse_load(rs.c1);
    _sse_vector_add();
    _sse_store(rs.c1);

    _sse_load(rs.c4);
    _sse_vector_i_mul();      
    _sse_vector_add();
    _sse_store(rs.c4);

    _sse_load((*sm).c2);
    _sse_load_up((*sm).c3);
    _sse_vector_i_mul();
    _sse_vector_sub();
      
    _sse_su3_inverse_multiply((*um));
    _sse_vector_cmplxcg_mul(ka1);
      
    _sse_load(rs.c2);
    _sse_vector_add();
    _sse_store(rs.c2);

    _sse_load(rs.c3);
    _sse_vector_i_mul();      
    _sse_vector_add();
    _sse_store(rs.c3);

    /*********************** direction +2 ************************/

    iy=g_idn[ix][2]; icy=trans1[iy]-ioff2;

    um=&g_gauge_field[iy][2];
    _prefetch_su3(um);

    sm=&spinor_field[k][icy];
    _prefetch_spinor(sm);

    _sse_load((*sp).c1);
    _sse_load_up((*sp).c4);
    _sse_vector_add();

    _sse_su3_multiply((*up));
    _sse_vector_cmplx_mul(ka2);

    _sse_load(rs.c1);
    _sse_vector_add();
    _sse_store(rs.c1);

    _sse_load(rs.c4);
    _sse_vector_add();
    _sse_store(rs.c4);
      
    _sse_load((*sp).c2);
    _sse_load_up((*sp).c3);
    _sse_vector_sub();

    _sse_su3_multiply((*up));
    _sse_vector_cmplx_mul(ka2);

    _sse_load(rs.c2);
    _sse_vector_add();
    _sse_store(rs.c2);

    _sse_load(rs.c3);
    _sse_vector_sub();
    _sse_store(rs.c3);      

    /*********************** direction -2 ************************/

    iy=g_iup[ix][3]; icy=trans1[iy]-ioff2;

    up+=1;
    _prefetch_su3(up);

    sp=&spinor_field[k][icy];
    _prefetch_spinor(sp);

    _sse_load((*sm).c1);
    _sse_load_up((*sm).c4);
    _sse_vector_sub();
      
    _sse_su3_inverse_multiply((*um));
    _sse_vector_cmplxcg_mul(ka2);
      
    _sse_load(rs.c1);
    _sse_vector_add();
    _sse_store(rs.c1);

    _sse_load(rs.c4);
    _sse_vector_sub();
    _sse_store(rs.c4);
      
    _sse_load((*sm).c2);
    _sse_load_up((*sm).c3);
    _sse_vector_add();
      
    _sse_su3_inverse_multiply((*um));
    _sse_vector_cmplxcg_mul(ka2);
      
    _sse_load(rs.c2);
    _sse_vector_add();
    _sse_store(rs.c2);

    _sse_load(rs.c3);
    _sse_vector_add();
    _sse_store(rs.c3);      
      
    /*********************** direction +3 ************************/

    iy=g_idn[ix][3]; icy=trans1[iy]-ioff2;

    um=&g_gauge_field[iy][3];
    _prefetch_su3(um);

    sm=&spinor_field[k][icy];
    _prefetch_spinor(sm);

    _sse_load((*sp).c1);
    _sse_load_up((*sp).c3);
    _sse_vector_i_mul();
    _sse_vector_add();

    _sse_su3_multiply((*up));
    _sse_vector_cmplx_mul(ka3);

    _sse_load(rs.c1);
    _sse_vector_add();
    _sse_store(rs.c1);

    _sse_load(rs.c3);
    _sse_vector_i_mul();      
    _sse_vector_sub();
    _sse_store(rs.c3);
      
    _sse_load((*sp).c2);
    _sse_load_up((*sp).c4);
    _sse_vector_i_mul();
    _sse_vector_sub();

    _sse_su3_multiply((*up));
    _sse_vector_cmplx_mul(ka3);

    _sse_load(rs.c2);
    _sse_vector_add();
    _sse_store(rs.c2);

    _sse_load(rs.c4);
    _sse_vector_i_mul();      
    _sse_vector_add();
    _sse_store(rs.c4);
      
    /*********************** direction -3 ************************/

    icz=icx+1;
    if(icz==((VOLUME+RAND)/2+ioff)) icz=ioff;
    iz=trans2[icz];
    iy=g_iup[iz][0]; icy=trans1[iy]-ioff2;

    up=&g_gauge_field[iz][0];
    _prefetch_su3(up);
      
    sp=&spinor_field[k][icy];
    _prefetch_spinor(sp);

    _sse_load((*sm).c1);
    _sse_load_up((*sm).c3);
    _sse_vector_i_mul();
    _sse_vector_sub();
      
    _sse_su3_inverse_multiply((*um));
    _sse_vector_cmplxcg_mul(ka3);

    rn=&spinor_field[l][icx-ioff];
      
    _sse_load(rs.c1);
    _sse_vector_add();
    _sse_store((*rn).c1);

    _sse_load(rs.c3);
    _sse_vector_i_mul();      
    _sse_vector_add();
    _sse_store((*rn).c3);

    _sse_load((*sm).c2);
    _sse_load_up((*sm).c4);
    _sse_vector_i_mul();
    _sse_vector_add();
      
    _sse_su3_inverse_multiply((*um));
    _sse_vector_cmplxcg_mul(ka3);

    _sse_load(rs.c2);
    _sse_vector_add();
    _sse_store((*rn).c2);

    _sse_load(rs.c4);
    _sse_vector_i_mul();      
    _sse_vector_sub();
    _sse_store((*rn).c4);
  }
}

/* If defined SSE2 */
#else

/* l output , k input*/
/* for ieo=0, k resides on  odd sites and l on even sites */
void H_eo(int ieo, int l, int k){
  int ix,iy;
  int ioff,ioff2,icx,icy;
  su3 *up,*um;
  spinor *r,*sp,*sm;

  /* for parallelization */
#ifdef MPI
  xchange_field(k);
#endif

  if(k == l){
    printf("Error in H_psi (simple.c):\n");
    printf("Arguments k and l must be different\n");
    printf("Program aborted\n");
    exit(1);
  }
  if(ieo == 0){
    ioff = 0;
  } 
  else{
    ioff = (VOLUME+RAND)/2;
  } 
  ioff2 = (VOLUME+RAND)/2-ioff;
  /**************** loop over all lattice sites ****************/

  for (icx = ioff; icx < (VOLUME/2 + ioff); icx++){
    ix=trans2[icx];

    r=&spinor_field[l][icx-ioff];

    /*********************** direction +0 ************************/

    iy=g_iup[ix][0]; icy=trans1[iy]-ioff2;


    sp=&spinor_field[k][icy];
    up=&g_gauge_field[ix][0];
      
    _vector_add(psi,(*sp).c1,(*sp).c3);

    _su3_multiply(chi,(*up),psi);
    _complex_times_vector(psi,ka0,chi);
      
    _vector_assign((*r).c1,psi);
    _vector_assign((*r).c3,psi);

    _vector_add(psi,(*sp).c2,(*sp).c4);

    _su3_multiply(chi,(*up),psi);
    _complex_times_vector(psi,ka0,chi);
            
    _vector_assign((*r).c2,psi);
    _vector_assign((*r).c4,psi);

    /*********************** direction -0 ************************/

    iy=g_idn[ix][0]; icy=trans1[iy]-ioff2;

    sm=&spinor_field[k][icy];
    um=&g_gauge_field[iy][0];
      
    _vector_sub(psi,(*sm).c1,(*sm).c3);

    _su3_inverse_multiply(chi,(*um),psi);
    _complexcjg_times_vector(psi,ka0,chi);

    _vector_add_assign((*r).c1,psi);
    _vector_sub_assign((*r).c3,psi);

    _vector_sub(psi,(*sm).c2,(*sm).c4);

    _su3_inverse_multiply(chi,(*um),psi);
    _complexcjg_times_vector(psi,ka0,chi);
      
    _vector_add_assign((*r).c2,psi);
    _vector_sub_assign((*r).c4,psi);

    /*********************** direction +1 ************************/

    iy=g_iup[ix][1]; icy=trans1[iy]-ioff2;

    sp=&spinor_field[k][icy];
    up=&g_gauge_field[ix][1];      
      
    _vector_i_add(psi,(*sp).c1,(*sp).c4);

    _su3_multiply(chi,(*up),psi);
    _complex_times_vector(psi,ka1,chi);

    _vector_add_assign((*r).c1,psi);
    _vector_i_sub_assign((*r).c4,psi);

    _vector_i_add(psi,(*sp).c2,(*sp).c3);

    _su3_multiply(chi,(*up),psi);
    _complex_times_vector(psi,ka1,chi);

    _vector_add_assign((*r).c2,psi);
    _vector_i_sub_assign((*r).c3,psi);

    /*********************** direction -1 ************************/

    iy=g_idn[ix][1]; icy=trans1[iy]-ioff2;

    sm=&spinor_field[k][icy];
    um=&g_gauge_field[iy][1];
      
    _vector_i_sub(psi,(*sm).c1,(*sm).c4);

    _su3_inverse_multiply(chi,(*um),psi);
    _complexcjg_times_vector(psi,ka1,chi);

    _vector_add_assign((*r).c1,psi);
    _vector_i_add_assign((*r).c4,psi);

    _vector_i_sub(psi,(*sm).c2,(*sm).c3);

    _su3_inverse_multiply(chi,(*um),psi);
    _complexcjg_times_vector(psi,ka1,chi);

    _vector_add_assign((*r).c2,psi);
    _vector_i_add_assign((*r).c3,psi);

    /*********************** direction +2 ************************/

    iy=g_iup[ix][2]; icy=trans1[iy]-ioff2;

    sp=&spinor_field[k][icy];
    up=&g_gauge_field[ix][2];
      
    _vector_add(psi,(*sp).c1,(*sp).c4);

    _su3_multiply(chi,(*up),psi);
    _complex_times_vector(psi,ka2,chi);

    _vector_add_assign((*r).c1,psi);
    _vector_add_assign((*r).c4,psi);

    _vector_sub(psi,(*sp).c2,(*sp).c3);

    _su3_multiply(chi,(*up),psi);
    _complex_times_vector(psi,ka2,chi);
      
    _vector_add_assign((*r).c2,psi);
    _vector_sub_assign((*r).c3,psi);

    /*********************** direction -2 ************************/

    iy=g_idn[ix][2]; icy=trans1[iy]-ioff2;

    sm=&spinor_field[k][icy];
    um=&g_gauge_field[iy][2];
      
    _vector_sub(psi,(*sm).c1,(*sm).c4);

    _su3_inverse_multiply(chi,(*um),psi);
    _complexcjg_times_vector(psi,ka2,chi);

    _vector_add_assign((*r).c1,psi);
    _vector_sub_assign((*r).c4,psi);

    _vector_add(psi,(*sm).c2,(*sm).c3);

    _su3_inverse_multiply(chi,(*um),psi);
    _complexcjg_times_vector(psi,ka2,chi);
      
    _vector_add_assign((*r).c2,psi);
    _vector_add_assign((*r).c3,psi);

    /*********************** direction +3 ************************/

    iy=g_iup[ix][3]; icy=trans1[iy]-ioff2;

    sp=&spinor_field[k][icy];
    up=&g_gauge_field[ix][3];
      
    _vector_i_add(psi,(*sp).c1,(*sp).c3);
      
    _su3_multiply(chi,(*up),psi);
    _complex_times_vector(psi,ka3,chi);

    _vector_add_assign((*r).c1,psi);
    _vector_i_sub_assign((*r).c3,psi);

    _vector_i_sub(psi,(*sp).c2,(*sp).c4);

    _su3_multiply(chi,(*up),psi);
    _complex_times_vector(psi,ka3,chi);

    _vector_add_assign((*r).c2,psi);
    _vector_i_add_assign((*r).c4,psi);

    /*********************** direction -3 ************************/

    iy=g_idn[ix][3]; icy=trans1[iy]-ioff2;

    sm=&spinor_field[k][icy];
    um=&g_gauge_field[iy][3];
      
    _vector_i_sub(psi,(*sm).c1,(*sm).c3);

    _su3_inverse_multiply(chi,(*um),psi);
    _complexcjg_times_vector(psi,ka3,chi);
      
    _vector_add_assign((*r).c1,psi);
    _vector_i_add_assign((*r).c3,psi);

    _vector_i_add(psi,(*sm).c2,(*sm).c4);

    _su3_inverse_multiply(chi,(*um),psi);
    _complexcjg_times_vector(psi,ka3,chi);

    _vector_add_assign((*r).c2,psi);
    _vector_i_sub_assign((*r).c4,psi);
    /************************ end of loop ************************/
  }
}
/* If defined SSE2 */
#endif

void clover_inv(int ieo, int l){
  int ix;
  int ioff,icx;
  su3 *w1,*w2,*w3;
  spinor *rn;

  /* Set the offset to work only */
  /* on even or odd sides */
  if(ieo == 0){
    ioff = 0;
  } 
  else{
    ioff = (VOLUME+RAND)/2;
  }
  /************ loop over all lattice sites ************/
  for(icx = ioff; icx < (VOLUME/2+ioff); icx++){
    ix=trans2[icx];
    rn=&spinor_field[l][icx-ioff];
    /* noch schlimmer murks fuer den clover-term */
    _vector_assign(phi1, (*rn).c1);
    _vector_assign(phi3, (*rn).c3);

    w1=&sw_inv[ix][0][0];
    w2=w1+2;  /* &sw_inv[ix][1][0]; */
    w3=w1+4;  /* &sw_inv[ix][2][0]; */
    _su3_multiply(psi, *w1, phi1); 
    _su3_multiply(chi, *w2, (*rn).c2);
    _vector_add((*rn).c1, psi,chi);
    _su3_inverse_multiply(psi, *w2,phi1); 
    _su3_multiply(chi, *w3, (*rn).c2);
    _vector_add((*rn).c2, psi, chi);

    w1++; /* &sw_inv[ix][0][1]; */
    w2++; /* &sw_inv[ix][1][1]; */
    w3++; /* &sw_inv[ix][2][1]; */
    _su3_multiply(psi, *w1, phi3); 
    _su3_multiply(chi, *w2, (*rn).c4);
    _vector_add((*rn).c3, psi, chi);
    _su3_inverse_multiply(psi, *w2, phi3); 
    _su3_multiply(chi,*w3,(*rn).c4);
    _vector_add((*rn).c4, psi, chi);
    /****************** end of loop *******************/
  }
}
/*
 * l is the number of the output field
 * k,j are the numbers of input fields
 *
 * ieo == 0 indicates to use T_{oo}
 *          but to work only on even sides of
 *          l
 * ieo == 1 the opposite
 * 
 *
 * computes l = gamma_5*((q_off*1+T_{ee|oo})*k-j)
 */

void clover_gamma5(int ieo, int l, int k, int j, double q_off){
  int ix;
  int ioff,icx;
  su3 *w1,*w2,*w3;
  spinor *r,*s,*t;

  /* Set the offset to work only */
  /* on even or odd sides */
  if(ieo==0){
    ioff = 0;
  } 
  else{
    ioff = (VOLUME+RAND)/2;
  }
  /******* loop over all lattice sites ********/
  for(icx = ioff; icx < (VOLUME/2+ioff); icx++){
    ix=trans2[icx];
    r=&spinor_field[l][icx-ioff];
    s=&spinor_field[k][icx-ioff];
    t=&spinor_field[j][icx-ioff];

    /* upper left part */
    w1=&sw[ix][0][0];
    w2=w1+2; /*&sw[ix][1][0];*/
    w3=w1+4; /*&sw[ix][2][0];*/
    _su3_multiply(psi1, *w1, (*s).c1); 
    _su3_multiply(chi, *w2, (*s).c2);
    _vector_add_assign(psi1, chi);
    _su3_inverse_multiply(psi2, *w2, (*s).c1); 
    _su3_multiply(chi, *w3, (*s).c2);
    _vector_add_assign(psi2, chi); 

    /*************** contribution from q_off ******************/
    _vector_add_mul(psi1, q_off, (*s).c1);
    _vector_add_mul(psi2, q_off,(*s).c2);

    /**********************************************************/
    /* Subtract t and store the result in r */
    _vector_sub((*r).c1, psi1, (*t).c1);
    _vector_sub((*r).c2, psi2, (*t).c2);

    /* lower right part */
    w1++; /*=&sw[ix][0][1];*/
    w2++; /*=&sw[ix][1][1];*/
    w3++; /*=&sw[ix][2][1];*/
    _su3_multiply(psi1, *w1, (*s).c3);
    _su3_multiply(chi, *w2, (*s).c4);
    _vector_add_assign(psi1, chi); 
    _su3_inverse_multiply(psi2, *w2, (*s).c3); 
    _su3_multiply(chi, *w3, (*s).c4);
    _vector_add_assign(psi2, chi); 

    /*************** contribution from q_off ******************/
    _vector_add_mul(psi1, q_off, (*s).c3);
    _vector_add_mul(psi2, q_off, (*s).c4);

    /* Subtract t and store the result in r */
    /* multiply with  gamma5 included by    */
    /* reversed order of t and psi1|2       */
    _vector_sub((*r).c3, (*t).c3, psi1);
    _vector_sub((*r).c4, (*t).c4, psi2);

    /************************ end of loop *********************/
  }
}

void clover(int ieo, int l, int k, int j, double q_off){
  int ix;
  int ioff,icx;
  su3 *w1,*w2,*w3;
  spinor *r,*s,*t;

  if(ieo==0){
    ioff = 0;
  } 
  else{
    ioff = (VOLUME+RAND)/2;
  }
  /**************** loop over all lattice sites ************/
  for(icx = ioff; icx < (VOLUME/2+ioff); icx++){
    ix=trans2[icx];
    r=&spinor_field[l][icx-ioff];
    s=&spinor_field[k][icx-ioff];
    t=&spinor_field[j][icx-ioff];

    w1=&sw[ix][0][0];
    w2=w1+2; /*&sw[ix][1][0];*/
    w3=w1+4; /*&sw[ix][2][0];*/
    _su3_multiply(psi1, *w1, (*s).c1); 
    _su3_multiply(chi, *w2, (*s).c2);
    _vector_add_assign(psi1, chi);
    _su3_inverse_multiply(psi2, *w2, (*s).c1); 
    _su3_multiply(chi, *w3, (*s).c2);
    _vector_add_assign(psi2, chi); 

    /*************** contribution from q_off *********/
    _vector_add_mul(psi1, q_off, (*s).c1);
    _vector_add_mul(psi2, q_off, (*s).c2);

    /*************************************************/
    /* Subtract t and store the result in r */
    _vector_sub((*r).c1, psi1, (*t).c1);
    _vector_sub((*r).c2, psi2, (*t).c2);

    w1++; /*=&sw[ix][0][1];*/
    w2++; /*=&sw[ix][1][1];*/
    w3++; /*=&sw[ix][2][1];*/
    _su3_multiply(psi1, *w1, (*s).c3); 
    _su3_multiply(chi, *w2, (*s).c4);
    _vector_add_assign(psi1, chi); 
    _su3_inverse_multiply(psi2, *w2, (*s).c3); 
    _su3_multiply(chi, *w3, (*s).c4);
    _vector_add_assign(psi2, chi); 

    /************** contribution from q_off **********/
    _vector_add_mul(psi1, q_off, (*s).c3);
    _vector_add_mul(psi2, q_off, (*s).c4);

    /*************************************************/
    /* Subtract t and store the result in r */
    _vector_sub((*r).c3, psi1, (*t).c3);
    _vector_sub((*r).c4, psi2, (*t).c4);

    /******************* end of loop *****************/
  }
}

/* output l , input k; k and l identical is permitted*/






void gamma5(int l,int k){
  int ix;
  spinor *r,*s;
  for (ix=0;ix<VOLUME/2;ix++){
    r=&spinor_field[l][ix];
    s=&spinor_field[k][ix];
    _vector_assign((*r).c1,(*s).c1);
    _vector_assign((*r).c2,(*s).c2);
    _vector_minus_assign((*r).c3,(*s).c3);
    _vector_minus_assign((*r).c4,(*s).c4);
  }
}

/* $Id$ */

/******************************************
 * Hopping_Matrix is the conventional Wilson 
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

#include <stdlib.h>
#include <stdio.h>
#include "global.h"
#include "su3.h"
#include "sse.h"
#ifdef MPI
#include "xchange.h"
#endif
#include "boundary.h"
#include "Hopping_Matrix.h"

static su3_vector psi1, psi2, psi, chi, phi1, phi3;

#if defined SSE2

/* input on k; output on l */
void Hopping_Matrix(const int ieo, const int l, const int k){
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

/* else of If defined SSE2 */
#else

/* l output , k input*/
/* for ieo=0, k resides on  odd sites and l on even sites */
void Hopping_Matrix(int ieo, int l, int k){
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
/* end of If defined SSE2 */
#endif

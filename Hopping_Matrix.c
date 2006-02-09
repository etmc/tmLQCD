/* $Id$ */

/******************************************
 * Hopping_Matrix is the conventional Wilson 
 * hopping matrix
 *
 * \kappa\sum_{\pm\mu}(r+\gamma_\mu)U_{x,\mu}
 *
 * for ieo = 0 this is M_{eo}, for ieo = 1
 * it is M_{oe}
 *
 * l is the number of the output field
 * k is the number of the input field
 *
 ******************************************/

/* #define BGL */

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include "global.h"
#include "su3.h"
#include "sse.h"
#ifdef MPI
#  include "xchange_field.h"
#endif
#include "boundary.h"
#include "Hopping_Matrix.h"

void xchange_field(spinor * const l, const int even);

#if ((defined SSE2)||(defined SSE3))

/* input on k; output on l */
void Hopping_Matrix(const int ieo, spinor * const l, spinor * const k){
  int icx,icy,icz,ioff,ioff2;
  int ix,iy,iz;
  su3 *up,*um;
  spinor *sp,*sm,*rn;
  static spinor rs;
  
  /* for parallelization */
#ifdef MPI
  xchange_field(k, ieo);
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

  ix=g_eo2lexic[ioff];
  iy=g_iup[ix][0]; 
  icy=g_lexic2eosub[iy];

  sp=k+icy;
  up=&g_gauge_field[ix][0];
   
  /**************** loop over all lattice sites ******************/
  for(icx = ioff; icx < (VOLUME/2+ioff); icx++){
    ix=g_eo2lexic[icx];
    /*********************** direction +0 ************************/

    iy=g_idn[ix][0]; icy=g_lexic2eosub[iy];

#ifndef _GAUGE_COPY
    um=&g_gauge_field[iy][0]; 
#else
    um=&g_gauge_field_back[ix][0]; 
#endif
    _prefetch_su3(um);
      
    sm=k+icy;
    _prefetch_spinor(sm);

    _sse_load((*sp).s0);
    _sse_load_up((*sp).s2);
    _sse_vector_add();

    _sse_su3_multiply((*up));
    _sse_vector_cmplx_mul(ka0);
    _sse_store_up(rs.s0);
    _sse_store_up(rs.s2);      
      
    _sse_load((*sp).s1);
    _sse_load_up((*sp).s3);
    _sse_vector_add();
      
    _sse_su3_multiply((*up));
    _sse_vector_cmplx_mul(ka0);
    _sse_store_up(rs.s1);
    _sse_store_up(rs.s3); 

    /*********************** direction -0 ************************/

    iy=g_iup[ix][1]; icy=g_lexic2eosub[iy];

    up+=1;
    _prefetch_su3(up);
      
    sp=k+icy;
    _prefetch_spinor(sp);

    _sse_load((*sm).s0);
    _sse_load_up((*sm).s2);
    _sse_vector_sub();
      
    _sse_su3_inverse_multiply((*um));
    _sse_vector_cmplxcg_mul(ka0);
      
    _sse_load(rs.s0);
    _sse_vector_add();
    _sse_store(rs.s0);

    _sse_load(rs.s2);
    _sse_vector_sub();
    _sse_store(rs.s2);
      
    _sse_load((*sm).s1);
    _sse_load_up((*sm).s3);
    _sse_vector_sub();
      
    _sse_su3_inverse_multiply((*um));
    _sse_vector_cmplxcg_mul(ka0);
      
    _sse_load(rs.s1);
    _sse_vector_add();
    _sse_store(rs.s1);

    _sse_load(rs.s3);
    _sse_vector_sub();
    _sse_store(rs.s3);
      
    /*********************** direction +1 ************************/

    iy=g_idn[ix][1]; icy=g_lexic2eosub[iy];

#ifndef _GAUGE_COPY
    um=&g_gauge_field[iy][1]; 
#else
    um+=1;
#endif
    _prefetch_su3(um);

    sm=k+icy;
    _prefetch_spinor(sm);

    _sse_load((*sp).s0);
    _sse_load_up((*sp).s3);
    _sse_vector_i_mul();
    _sse_vector_add();

    _sse_su3_multiply((*up));
    _sse_vector_cmplx_mul(ka1);

    _sse_load(rs.s0);
    _sse_vector_add();
    _sse_store(rs.s0);

    _sse_load(rs.s3);
    _sse_vector_i_mul();      
    _sse_vector_sub();
    _sse_store(rs.s3); 
      
    _sse_load((*sp).s1);
    _sse_load_up((*sp).s2);
    _sse_vector_i_mul();
    _sse_vector_add();

    _sse_su3_multiply((*up));
    _sse_vector_cmplx_mul(ka1);

    _sse_load(rs.s1);
    _sse_vector_add();
    _sse_store(rs.s1);

    _sse_load(rs.s2);
    _sse_vector_i_mul();      
    _sse_vector_sub();
    _sse_store(rs.s2);       

    /*********************** direction -1 ************************/

    iy=g_iup[ix][2]; icy=g_lexic2eosub[iy];

    up+=1;
    _prefetch_su3(up);

    sp=k+icy;
    _prefetch_spinor(sp);

    _sse_load((*sm).s0);
    _sse_load_up((*sm).s3);
    _sse_vector_i_mul();
    _sse_vector_sub();
      
    _sse_su3_inverse_multiply((*um));
    _sse_vector_cmplxcg_mul(ka1);
      
    _sse_load(rs.s0);
    _sse_vector_add();
    _sse_store(rs.s0);

    _sse_load(rs.s3);
    _sse_vector_i_mul();      
    _sse_vector_add();
    _sse_store(rs.s3);

    _sse_load((*sm).s1);
    _sse_load_up((*sm).s2);
    _sse_vector_i_mul();
    _sse_vector_sub();
      
    _sse_su3_inverse_multiply((*um));
    _sse_vector_cmplxcg_mul(ka1);
      
    _sse_load(rs.s1);
    _sse_vector_add();
    _sse_store(rs.s1);

    _sse_load(rs.s2);
    _sse_vector_i_mul();      
    _sse_vector_add();
    _sse_store(rs.s2);

    /*********************** direction +2 ************************/

    iy=g_idn[ix][2]; icy=g_lexic2eosub[iy];

#ifndef _GAUGE_COPY
    um=&g_gauge_field[iy][2]; 
#else
    um+=1;
#endif
    _prefetch_su3(um);

    sm=k+icy;
    _prefetch_spinor(sm);

    _sse_load((*sp).s0);
    _sse_load_up((*sp).s3);
    _sse_vector_add();

    _sse_su3_multiply((*up));
    _sse_vector_cmplx_mul(ka2);

    _sse_load(rs.s0);
    _sse_vector_add();
    _sse_store(rs.s0);

    _sse_load(rs.s3);
    _sse_vector_add();
    _sse_store(rs.s3);
      
    _sse_load((*sp).s1);
    _sse_load_up((*sp).s2);
    _sse_vector_sub();

    _sse_su3_multiply((*up));
    _sse_vector_cmplx_mul(ka2);

    _sse_load(rs.s1);
    _sse_vector_add();
    _sse_store(rs.s1);

    _sse_load(rs.s2);
    _sse_vector_sub();
    _sse_store(rs.s2);      

    /*********************** direction -2 ************************/

    iy=g_iup[ix][3]; icy=g_lexic2eosub[iy];

    up+=1;
    _prefetch_su3(up);

    sp=k+icy;
    _prefetch_spinor(sp);

    _sse_load((*sm).s0);
    _sse_load_up((*sm).s3);
    _sse_vector_sub();
      
    _sse_su3_inverse_multiply((*um));
    _sse_vector_cmplxcg_mul(ka2);
      
    _sse_load(rs.s0);
    _sse_vector_add();
    _sse_store(rs.s0);

    _sse_load(rs.s3);
    _sse_vector_sub();
    _sse_store(rs.s3);
      
    _sse_load((*sm).s1);
    _sse_load_up((*sm).s2);
    _sse_vector_add();
      
    _sse_su3_inverse_multiply((*um));
    _sse_vector_cmplxcg_mul(ka2);
      
    _sse_load(rs.s1);
    _sse_vector_add();
    _sse_store(rs.s1);

    _sse_load(rs.s2);
    _sse_vector_add();
    _sse_store(rs.s2);      
      
    /*********************** direction +3 ************************/

    iy=g_idn[ix][3]; icy=g_lexic2eosub[iy];

#ifndef _GAUGE_COPY
    um=&g_gauge_field[iy][3]; 
#else
    um+=1;
#endif
    _prefetch_su3(um);

    sm=k+icy;
    _prefetch_spinor(sm);

    _sse_load((*sp).s0);
    _sse_load_up((*sp).s2);
    _sse_vector_i_mul();
    _sse_vector_add();

    _sse_su3_multiply((*up));
    _sse_vector_cmplx_mul(ka3);

    _sse_load(rs.s0);
    _sse_vector_add();
    _sse_store(rs.s0);

    _sse_load(rs.s2);
    _sse_vector_i_mul();      
    _sse_vector_sub();
    _sse_store(rs.s2);
      
    _sse_load((*sp).s1);
    _sse_load_up((*sp).s3);
    _sse_vector_i_mul();
    _sse_vector_sub();

    _sse_su3_multiply((*up));
    _sse_vector_cmplx_mul(ka3);

    _sse_load(rs.s1);
    _sse_vector_add();
    _sse_store(rs.s1);

    _sse_load(rs.s3);
    _sse_vector_i_mul();      
    _sse_vector_add();
    _sse_store(rs.s3);
      
    /*********************** direction -3 ************************/

    icz=icx+1;
    if(icz==((VOLUME+RAND)/2+ioff)) icz=ioff;
    iz=g_eo2lexic[icz];
    iy=g_iup[iz][0]; icy=g_lexic2eosub[iy];

    up=&g_gauge_field[iz][0];
    _prefetch_su3(up);
      
    sp=k+icy;
    _prefetch_spinor(sp);

    _sse_load((*sm).s0);
    _sse_load_up((*sm).s2);
    _sse_vector_i_mul();
    _sse_vector_sub();
      
    _sse_su3_inverse_multiply((*um));
    _sse_vector_cmplxcg_mul(ka3);

    rn=l+(icx-ioff);
      
    _sse_load(rs.s0);
    _sse_vector_add();
    _sse_store((*rn).s0);

    _sse_load(rs.s2);
    _sse_vector_i_mul();      
    _sse_vector_add();
    _sse_store((*rn).s2);

    _sse_load((*sm).s1);
    _sse_load_up((*sm).s3);
    _sse_vector_i_mul();
    _sse_vector_add();
      
    _sse_su3_inverse_multiply((*um));
    _sse_vector_cmplxcg_mul(ka3);

    _sse_load(rs.s1);
    _sse_vector_add();
    _sse_store((*rn).s1);

    _sse_load(rs.s3);
    _sse_vector_i_mul();      
    _sse_vector_sub();
    _sse_store((*rn).s3);
  }
}

#elif (defined BGL && defined XLC)

#include "bgl.h"

#undef ALIGN
#define ALIGN __attribute__ ((aligned (16)))

void Hopping_Matrix(const int ieo, spinor * const l, spinor * const k){
  int icx,icy,icz,ioff,ioff2;
  int ix,iy,iz;
  su3 *up ALIGN;
  su3 *um ALIGN;
  spinor *sp ALIGN;
  spinor *sm ALIGN;
  spinor *rn ALIGN;
  double _Complex reg00, reg01, reg02, reg03, reg04, reg05, reg06, reg07;
  double _Complex reg10, reg20;
  double _Complex rs00, rs01, rs02, rs10, rs11, rs12, rs20, rs21, rs22, 
    rs30, rs31, rs32;

#pragma disjoint(*sp, *sm, *rn, *up, *um, *l, *k)

  __alignx(16,l);
  __alignx(16,k);

#ifdef MPI
  xchange_field(k, ieo);
#endif

  if(ieo == 0){
    ioff = 0;
  } 
  else{
    ioff = (VOLUME+RAND)/2;
  }
  ioff2 = (VOLUME+RAND)/2-ioff;

  ix=g_eo2lexic[ioff];
  iy=g_iup[ix][0]; 
  icy=g_lexic2eosub[iy];

  sp=k+icy;
  up=&g_gauge_field[ix][0];
   
  /**************** loop over all lattice sites ******************/
  for(icx = ioff; icx < (VOLUME/2+ioff); icx++){
    ix=g_eo2lexic[icx];
    /*********************** direction +0 ************************/
    iy=g_idn[ix][0]; 
    icy=g_lexic2eosub[iy];
#ifndef _GAUGE_COPY
    um=&g_gauge_field[iy][0]; 
#else
    um=&g_gauge_field_back[ix][0]; 
#endif
    _prefetch_su3(um); 
    sm=k+icy;
    _prefetch_spinor(sm); 

    _bgl_load((*sp).s0);
    _bgl_load_up((*sp).s2);
    _bgl_vector_add();

    _bgl_su3_multiply((*up));
    _bgl_vector_cmplx_mul(ka0);
    _bgl_store_up_rs0();
    _bgl_store_up_rs1();

    _bgl_load((*sp).s1);
    _bgl_load_up((*sp).s3);
    _bgl_vector_add();
      
    _bgl_su3_multiply((*up));
    _bgl_vector_cmplx_mul(ka0);
    _bgl_store_up_rs1();
    _bgl_store_up_rs3();


    /*********************** direction -0 ************************/

    iy=g_iup[ix][1]; 
    icy=g_lexic2eosub[iy];

    up+=1;
    _prefetch_su3(up); 
    sp=k+icy;
    _prefetch_spinor(sp); 

    _bgl_load((*sm).s0);
    _bgl_load_up((*sm).s2);
    _bgl_vector_sub();
      
    _bgl_su3_inverse_multiply((*um));
    _bgl_vector_cmplxcg_mul(ka0);

    _bgl_add_to_rs0();
    _bgl_sub_from_rs2();
      
    _bgl_load((*sm).s1);
    _bgl_load_up((*sm).s3);
    _bgl_vector_sub();
      
    _bgl_su3_inverse_multiply((*um));
    _bgl_vector_cmplxcg_mul(ka0);

    _bgl_add_to_rs1();
    _bgl_sub_from_rs3();

    /*********************** direction +1 ************************/

    iy=g_idn[ix][1]; 
    icy=g_lexic2eosub[iy];

#ifndef _GAUGE_COPY
    um=&g_gauge_field[iy][1]; 
#else
    um+=1;
#endif
    _prefetch_su3(um); 
    sm=k+icy;
    _prefetch_spinor(sm); 

    _bgl_load((*sp).s0);
    _bgl_load_up((*sp).s3);
/*     _bgl_vector_i_mul(); */
    _bgl_vector_i_mul_add();

    _bgl_su3_multiply((*up));
    _bgl_vector_cmplx_mul(ka1);

    _bgl_add_to_rs0();
/*     _bgl_vector_i_mul(); */
    _bgl_i_mul_sub_from_rs3();

    _bgl_load((*sp).s1);
    _bgl_load_up((*sp).s2);
/*     _bgl_vector_i_mul(); */
    _bgl_vector_i_mul_add();

    _bgl_su3_multiply((*up));
    _bgl_vector_cmplx_mul(ka1);

    _bgl_add_to_rs1();
    
/*     _bgl_vector_i_mul(); */
    _bgl_i_mul_sub_from_rs2()

    /*********************** direction -1 ************************/

    iy=g_iup[ix][2]; 
    icy=g_lexic2eosub[iy];

    up+=1;
    _prefetch_su3(up); 
    sp=k+icy;
    _prefetch_spinor(sp); 

    _bgl_load((*sm).s0);
    _bgl_load_up((*sm).s3);
/*     _bgl_vector_i_mul(); */
    _bgl_vector_i_mul_sub();
      
    _bgl_su3_inverse_multiply((*um));
    _bgl_vector_cmplxcg_mul(ka1);
      
    _bgl_add_to_rs0();
/*     _bgl_vector_i_mul();       */
    _bgl_i_mul_add_to_rs3();

    _bgl_load((*sm).s1);
    _bgl_load_up((*sm).s2);
/*     _bgl_vector_i_mul(); */
    _bgl_vector_i_mul_sub();
      
    _bgl_su3_inverse_multiply((*um));
    _bgl_vector_cmplxcg_mul(ka1);
      
    _bgl_add_to_rs1();

/*     _bgl_vector_i_mul(); */
    _bgl_i_mul_add_to_rs2();      

    /*********************** direction +2 ************************/

    iy=g_idn[ix][2]; 
    icy=g_lexic2eosub[iy];

#ifndef _GAUGE_COPY
    um=&g_gauge_field[iy][2]; 
#else
    um+=1;
#endif
    _prefetch_su3(um); 
    sm=k+icy;
    _prefetch_spinor(sm); 

    _bgl_load((*sp).s0);
    _bgl_load_up((*sp).s3);
    _bgl_vector_add();

    _bgl_su3_multiply((*up));
    _bgl_vector_cmplx_mul(ka2);

    _bgl_add_to_rs0();
    _bgl_add_to_rs3();

    _bgl_load((*sp).s1);
    _bgl_load_up((*sp).s2);
    _bgl_vector_sub();

    _bgl_su3_multiply((*up));
    _bgl_vector_cmplx_mul(ka2);

    _bgl_add_to_rs1();
    _bgl_sub_from_rs2();


    /*********************** direction -2 ************************/

    iy=g_iup[ix][3]; 
    icy=g_lexic2eosub[iy];

    up+=1;
    _prefetch_su3(up); 
    sp=k+icy;
    _prefetch_spinor(sp); 

    _bgl_load((*sm).s0);
    _bgl_load_up((*sm).s3);
    _bgl_vector_sub();
      
    _bgl_su3_inverse_multiply((*um));
    _bgl_vector_cmplxcg_mul(ka2);
      
    _bgl_add_to_rs0();
    _bgl_sub_from_rs3();

    _bgl_load((*sm).s1);
    _bgl_load_up((*sm).s2);
    _bgl_vector_add();
      
    _bgl_su3_inverse_multiply((*um));
    _bgl_vector_cmplxcg_mul(ka2);
      
    _bgl_add_to_rs1();

    _bgl_add_to_rs2();
      
    /*********************** direction +3 ************************/

    iy=g_idn[ix][3]; 
    icy=g_lexic2eosub[iy];

#ifndef _GAUGE_COPY
    um=&g_gauge_field[iy][3]; 
#else
    um+=1;
#endif
    _prefetch_su3(um); 
    sm=k+icy;
    _prefetch_spinor(sm); 

    _bgl_load((*sp).s0);
    _bgl_load_up((*sp).s2);
/*     _bgl_vector_i_mul(); */
    _bgl_vector_i_mul_add();

    _bgl_su3_multiply((*up));
    _bgl_vector_cmplx_mul(ka3);

    _bgl_add_to_rs0();

/*     _bgl_vector_i_mul();       */
    _bgl_i_mul_sub_from_rs3();

    _bgl_load((*sp).s1);
    _bgl_load_up((*sp).s3);
/*     _bgl_vector_i_mul(); */
    _bgl_vector_i_mul_sub();

    _bgl_su3_multiply((*up));
    _bgl_vector_cmplx_mul(ka3);

    _bgl_add_to_rs1();

/*     _bgl_vector_i_mul();       */
    _bgl_i_mul_add_to_rs3();
      
    /*********************** direction -3 ************************/

    icz=icx+1;
    if(icz==((VOLUME+RAND)/2+ioff)) icz=ioff;
    iz=g_eo2lexic[icz];
    iy=g_iup[iz][0]; icy=g_lexic2eosub[iy];

    up=&g_gauge_field[iz][0];
    _prefetch_su3(up); 
    sp=k+icy;
    _prefetch_spinor(sp); 

    _bgl_load((*sm).s0);
    _bgl_load_up((*sm).s2);
/*     _bgl_vector_i_mul(); */
    _bgl_vector_i_mul_sub();
      
    _bgl_su3_inverse_multiply((*um));
    _bgl_vector_cmplxcg_mul(ka3);

    rn=l+(icx-ioff);
      
    _bgl_add_to_rs0();
    _bgl_store_rs0((*rn).s0);

/*     _bgl_vector_i_mul();      */
    _bgl_i_mul_add_to_rs2();
    _bgl_store_rs2((*rn).s2);

    _bgl_load((*sm).s1);
    _bgl_load_up((*sm).s3);
/*     _bgl_vector_i_mul(); */
    _bgl_vector_i_mul_add();
      
    _bgl_su3_inverse_multiply((*um));
    _bgl_vector_cmplxcg_mul(ka3);

    _bgl_add_to_rs1();
    _bgl_store_rs1((*rn).s1);

/*     _bgl_vector_i_mul();       */
    _bgl_i_mul_sub_from_rs3();
    _bgl_store_rs3((*rn).s3);

  }
}
#elif defined XLC

static su3_vector psi1, psi2, psi, chi, phi1, phi3;

#include"xlc_prefetch.h"

/* l output , k input*/
/* for ieo=0, k resides on  odd sites and l on even sites */
void Hopping_Matrix(int ieo, spinor * const l, spinor * const k){
  int ix, iy, iz;
  int ioff, ioff2, icx, icy, icz;
  su3 *up,*um;
  spinor *r,*sp,*sm;
  void *dd;
  spinor temp;
#pragma disjoint(temp, *sp, *sm, *r, *up, *um)

  double _Complex reg0, reg1, reg2, reg3, reg4, reg5, reg6, reg7;

  /* for parallelization */
#ifdef MPI
  xchange_field(k, ieo);
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

  ix = g_eo2lexic[ioff];

  iy=g_iup[ix][0]; icy=g_lexic2eosub[iy];

  sp=k+icy;
  up=&g_gauge_field[ix][0];

  for (icx = ioff; icx < (VOLUME/2 + ioff); icx++) {
    ix = g_eo2lexic[icx];
    /*********************** direction +0 ************************/
    
    iy=g_idn[ix][0]; icy=g_lexic2eosub[iy];
    
    sm=k+icy;
    dd=&(*sm).s2.c2.re;
    _prefetch_spinor_dcbt((void*)sm, dd);

#ifndef _GAUGE_COPY
    um=&g_gauge_field[iy][0];
    dd=&g_gauge_field[iy][0].c22.re;
#else
    um = &g_gauge_field_back[ix][0];
    dd = &g_gauge_field_back[ix][0].c22.re;
#endif
    _prefetch_su3_dcbt((void*)um, dd);

    _vector_add(psi,(*sp).s0,(*sp).s2);

    _su3_multiply(chi,(*up),psi);
    _complex_times_vector(psi,ka0,chi);

    _vector_assign(temp.s0,psi);
    _vector_assign(temp.s2,psi);

    _vector_add(psi,(*sp).s1,(*sp).s3);

    _su3_multiply(chi,(*up),psi);
    _complex_times_vector(psi,ka0,chi);

    _vector_assign(temp.s1,psi);
    _vector_assign(temp.s3,psi);

    /*********************** direction -0 ************************/

    iy=g_iup[ix][1]; icy=g_lexic2eosub[iy];

    sp=k+icy;
    dd=&(*sp).s2.c2.re;
    _prefetch_spinor_dcbt((void*)sp, dd);

    up=&g_gauge_field[ix][1];     
    dd=&g_gauge_field[iy][1].c22.re;
    _prefetch_su3_dcbt((void*)up, dd); 
      
    _vector_sub(psi,(*sm).s0,(*sm).s2);

    _su3_inverse_multiply(chi,(*um),psi);
    _complexcjg_times_vector(psi,ka0,chi);

    _vector_add_assign(temp.s0,psi);
    _vector_sub_assign(temp.s2,psi);

    _vector_sub(psi,(*sm).s1,(*sm).s3);

    _su3_inverse_multiply(chi,(*um),psi);
    _complexcjg_times_vector(psi,ka0,chi);

    _vector_add_assign(temp.s1,psi);
    _vector_sub_assign(temp.s3,psi);

    /*********************** direction +1 ************************/

    iy=g_idn[ix][1]; icy=g_lexic2eosub[iy];

    sm=k+icy;
    dd=&(*sm).s2.c2.re;
    _prefetch_spinor_dcbt((void*)sm, dd);

#ifndef _GAUGE_COPY
    um=&g_gauge_field[iy][1];
    dd=&g_gauge_field[iy][1].c22.re;
#else
    um=&g_gauge_field_back[ix][1];
    dd=&g_gauge_field_back[ix][1].c22.re;
#endif
    _prefetch_su3_dcbt((void*)um, dd);
      
    _vector_i_add(psi,(*sp).s0,(*sp).s3);

    _su3_multiply(chi,(*up),psi);
    _complex_times_vector(psi,ka1,chi);

    _vector_add_assign(temp.s0,psi);
    _vector_i_sub_assign(temp.s3,psi);

    _vector_i_add(psi,(*sp).s1,(*sp).s2);

    _su3_multiply(chi,(*up),psi);
    _complex_times_vector(psi,ka1,chi);

    _vector_add_assign(temp.s1,psi);
    _vector_i_sub_assign(temp.s2,psi);

    /*********************** direction -1 ************************/

    iy=g_iup[ix][2]; icy=g_lexic2eosub[iy];

    sp=k+icy;
    dd=&(*sp).s2.c2.re;
    _prefetch_spinor_dcbt((void*)sp, dd);

    up=&g_gauge_field[ix][2];
    dd=&g_gauge_field[iy][2].c22.re;
    _prefetch_su3_dcbt((void*)up, dd);

    _vector_i_sub(psi,(*sm).s0,(*sm).s3);

    _su3_inverse_multiply(chi,(*um),psi);
    _complexcjg_times_vector(psi,ka1,chi);

    _vector_add_assign(temp.s0,psi);
    _vector_i_add_assign(temp.s3,psi);

    _vector_i_sub(psi,(*sm).s1,(*sm).s2);

    _su3_inverse_multiply(chi,(*um),psi);
    _complexcjg_times_vector(psi,ka1,chi);

    _vector_add_assign(temp.s1,psi);
    _vector_i_add_assign(temp.s2,psi);

    /*********************** direction +2 ************************/

    iy=g_idn[ix][2]; icy=g_lexic2eosub[iy];

    sm=k+icy;
    dd=&(*sm).s2.c2.re;
    _prefetch_spinor_dcbt((void*)sm, dd);

#ifndef _GAUGE_COPY
    um=&g_gauge_field[iy][2];
    dd=&g_gauge_field[iy][2].c22.re;
#else
    um=&g_gauge_field_back[ix][2];
    dd=&g_gauge_field_back[ix][2].c22.re;
#endif
    _prefetch_su3_dcbt((void*)um, dd);

    _vector_add(psi,(*sp).s0,(*sp).s3);

    _su3_multiply(chi,(*up),psi);
    _complex_times_vector(psi,ka2,chi);

    _vector_add_assign(temp.s0,psi);
    _vector_add_assign(temp.s3,psi);

    _vector_sub(psi,(*sp).s1,(*sp).s2);

    _su3_multiply(chi,(*up),psi);
    _complex_times_vector(psi,ka2,chi);

    _vector_add_assign(temp.s1,psi);
    _vector_sub_assign(temp.s2,psi);

    /*********************** direction -2 ************************/

    iy=g_iup[ix][3]; icy=g_lexic2eosub[iy];

    sp=k+icy;
    dd=&(*sp).s2.c2.re;
    _prefetch_spinor_dcbt((void*)sp, dd);

    up=&g_gauge_field[ix][3];
    dd=&g_gauge_field[iy][3].c22.re;
    _prefetch_su3_dcbt((void*)up, dd);

    _vector_sub(psi,(*sm).s0,(*sm).s3);

    _su3_inverse_multiply(chi,(*um),psi);
    _complexcjg_times_vector(psi,ka2,chi);

    _vector_add_assign(temp.s0,psi);
    _vector_sub_assign(temp.s3,psi);

    _vector_add(psi,(*sm).s1,(*sm).s2);

    _su3_inverse_multiply(chi,(*um),psi);
    _complexcjg_times_vector(psi,ka2,chi);

    _vector_add_assign(temp.s1,psi);
    _vector_add_assign(temp.s2,psi);

    /*********************** direction +3 ************************/

    iy=g_idn[ix][3]; icy=g_lexic2eosub[iy];

    sm=k+icy;
    dd=&(*sm).s2.c2.re;
    _prefetch_spinor_dcbt((void*)sm, dd);

#ifndef _GAUGE_COPY
    um=&g_gauge_field[iy][3];
    dd=&g_gauge_field[iy][3].c22.re;
#else
    um=&g_gauge_field_back[ix][3];
    dd=&g_gauge_field_back[ix][3].c22.re;
#endif
    _prefetch_su3_dcbt((void*)um, dd);

    _vector_i_add(psi,(*sp).s0,(*sp).s2);
      
    _su3_multiply(chi,(*up),psi);
    _complex_times_vector(psi,ka3,chi);

    _vector_add_assign(temp.s0,psi);
    _vector_i_sub_assign(temp.s2,psi);

    _vector_i_sub(psi,(*sp).s1,(*sp).s3);

    _su3_multiply(chi,(*up),psi);
    _complex_times_vector(psi,ka3,chi);

    _vector_add_assign(temp.s1,psi); 
    _vector_i_add_assign(temp.s3,psi); 

    /*********************** direction -3 ************************/

    r=l+(icx - ioff);
    dd=&(*r).s2.c2.re;
    _prefetch_spinor_dcbt((void*)r, dd);

    icz = icx + 1;
    if(icz==((VOLUME+RAND)/2+ioff)) icz=ioff;
    iz=g_eo2lexic[icz];
    iy=g_iup[iz][0]; icy=g_lexic2eosub[iy];

    up=&g_gauge_field[iz][0];
    dd=&g_gauge_field[iz][0].c22.re;
    _prefetch_su3_dcbt((void*)up, dd);

    sp=k+icy;
    dd=&(*sp).s2.c2.re;
    _prefetch_spinor_dcbt((void*)sp, dd);

    _vector_i_sub(psi,(*sm).s0,(*sm).s2);

    _su3_inverse_multiply(chi,(*um),psi);
    _complexcjg_times_vector(psi,ka3,chi);
      
    _vector_add((*r).s0, temp.s0, psi);
    _vector_i_add((*r).s2, temp.s2, psi);

    _vector_i_add(psi,(*sm).s1,(*sm).s3);

    _su3_inverse_multiply(chi,(*um),psi);
    _complexcjg_times_vector(psi,ka3,chi);

    _vector_add((*r).s1, temp.s1, psi);
    _vector_i_sub((*r).s3, temp.s3, psi);

    /************************ end of loop ************************/
  }
}


/* else of If defined SSE2  and if defined XLC */
#else

static su3_vector psi1, psi2, psi, chi, phi1, phi3;

/* l output , k input*/
/* for ieo=0, k resides on  odd sites and l on even sites */
void Hopping_Matrix(int ieo, spinor * const l, spinor * const k){
  int ix,iy;
  int ioff,ioff2,icx,icy;
  su3 *up,*um;
  spinor *r,*sp,*sm;
  spinor temp;

  /* for parallelization */
#ifdef MPI
  xchange_field(k, ieo);
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
    ix=g_eo2lexic[icx];

    r=l+(icx-ioff);

    /*********************** direction +0 ************************/

    iy=g_iup[ix][0]; icy=g_lexic2eosub[iy];


    sp=k+icy;
    up=&g_gauge_field[ix][0];
      
    _vector_add(psi,(*sp).s0,(*sp).s2);

    _su3_multiply(chi,(*up),psi);
    _complex_times_vector(psi,ka0,chi);
      
    _vector_assign(temp.s0,psi);
    _vector_assign(temp.s2,psi);

    _vector_add(psi,(*sp).s1,(*sp).s3);

    _su3_multiply(chi,(*up),psi);
    _complex_times_vector(psi,ka0,chi);
            
    _vector_assign(temp.s1,psi);
    _vector_assign(temp.s3,psi);

    /*********************** direction -0 ************************/

    iy=g_idn[ix][0]; icy=g_lexic2eosub[iy];

    sm=k+icy;
#ifndef _GAUGE_COPY
    um=&g_gauge_field[iy][0];
#else
    um = &g_gauge_field_back[ix][0];
#endif

    _vector_sub(psi,(*sm).s0,(*sm).s2);

    _su3_inverse_multiply(chi,(*um),psi);
    _complexcjg_times_vector(psi,ka0,chi);

    _vector_add_assign(temp.s0,psi);
    _vector_sub_assign(temp.s2,psi);

    _vector_sub(psi,(*sm).s1,(*sm).s3);

    _su3_inverse_multiply(chi,(*um),psi);
    _complexcjg_times_vector(psi,ka0,chi);
      
    _vector_add_assign(temp.s1,psi);
    _vector_sub_assign(temp.s3,psi);

    /*********************** direction +1 ************************/

    iy=g_iup[ix][1]; icy=g_lexic2eosub[iy];

    sp=k+icy;
    up=&g_gauge_field[ix][1];      
      
    _vector_i_add(psi,(*sp).s0,(*sp).s3);

    _su3_multiply(chi,(*up),psi);
    _complex_times_vector(psi,ka1,chi);

    _vector_add_assign(temp.s0,psi);
    _vector_i_sub_assign(temp.s3,psi);

    _vector_i_add(psi,(*sp).s1,(*sp).s2);

    _su3_multiply(chi,(*up),psi);
    _complex_times_vector(psi,ka1,chi);

    _vector_add_assign(temp.s1,psi);
    _vector_i_sub_assign(temp.s2,psi);

    /*********************** direction -1 ************************/

    iy=g_idn[ix][1]; icy=g_lexic2eosub[iy];

    sm=k+icy;
#ifndef _GAUGE_COPY
    um=&g_gauge_field[iy][1];
#else
    um = &g_gauge_field_back[ix][1];
#endif

    _vector_i_sub(psi,(*sm).s0,(*sm).s3);

    _su3_inverse_multiply(chi,(*um),psi);
    _complexcjg_times_vector(psi,ka1,chi);

    _vector_add_assign(temp.s0,psi);
    _vector_i_add_assign(temp.s3,psi);

    _vector_i_sub(psi,(*sm).s1,(*sm).s2);

    _su3_inverse_multiply(chi,(*um),psi);
    _complexcjg_times_vector(psi,ka1,chi);

    _vector_add_assign(temp.s1,psi);
    _vector_i_add_assign(temp.s2,psi);

    /*********************** direction +2 ************************/

    iy=g_iup[ix][2]; icy=g_lexic2eosub[iy];

    sp=k+icy;
    up=&g_gauge_field[ix][2];
      
    _vector_add(psi,(*sp).s0,(*sp).s3);

    _su3_multiply(chi,(*up),psi);
    _complex_times_vector(psi,ka2,chi);

    _vector_add_assign(temp.s0,psi);
    _vector_add_assign(temp.s3,psi);

    _vector_sub(psi,(*sp).s1,(*sp).s2);

    _su3_multiply(chi,(*up),psi);
    _complex_times_vector(psi,ka2,chi);
      
    _vector_add_assign(temp.s1,psi);
    _vector_sub_assign(temp.s2,psi);


    /*********************** direction -2 ************************/

    iy=g_idn[ix][2]; icy=g_lexic2eosub[iy];

    sm=k+icy;
#ifndef _GAUGE_COPY
    um = &g_gauge_field[iy][2];
#else
    um = &g_gauge_field_back[ix][2];
#endif

    _vector_sub(psi,(*sm).s0,(*sm).s3);

    _su3_inverse_multiply(chi,(*um),psi);
    _complexcjg_times_vector(psi,ka2,chi);

    _vector_add_assign(temp.s0,psi);
    _vector_sub_assign(temp.s3,psi);

    _vector_add(psi,(*sm).s1,(*sm).s2);

    _su3_inverse_multiply(chi,(*um),psi);
    _complexcjg_times_vector(psi,ka2,chi);
      
    _vector_add_assign(temp.s1,psi);
    _vector_add_assign(temp.s2,psi);

    /*********************** direction +3 ************************/

    iy=g_iup[ix][3]; icy=g_lexic2eosub[iy];

    sp=k+icy;
    up=&g_gauge_field[ix][3];
      
    _vector_i_add(psi,(*sp).s0,(*sp).s2);
      
    _su3_multiply(chi,(*up),psi);
    _complex_times_vector(psi,ka3,chi);

    _vector_add_assign(temp.s0,psi);
    _vector_i_sub_assign(temp.s2,psi);

    _vector_i_sub(psi,(*sp).s1,(*sp).s3);

    _su3_multiply(chi,(*up),psi);
    _complex_times_vector(psi,ka3,chi);

    _vector_add_assign(temp.s1,psi);
    _vector_i_add_assign(temp.s3,psi);

    /*********************** direction -3 ************************/

    iy=g_idn[ix][3]; icy=g_lexic2eosub[iy];

    sm=k+icy;
#ifndef _GAUGE_COPY
    um = &g_gauge_field[iy][3];
#else
    um = &g_gauge_field_back[ix][3];
#endif

    _vector_i_sub(psi,(*sm).s0,(*sm).s2);

    _su3_inverse_multiply(chi,(*um),psi);
    _complexcjg_times_vector(psi,ka3,chi);
      
#ifdef OlD
    _vector_add_assign((*r).s0,psi);
    _vector_i_add_assign((*r).s2,psi);
#else
    _vector_add((*r).s0, temp.s0, psi);
    _vector_i_add((*r).s2, temp.s2, psi);
#endif

    _vector_i_add(psi,(*sm).s1,(*sm).s3);

    _su3_inverse_multiply(chi,(*um),psi);
    _complexcjg_times_vector(psi,ka3,chi);

#ifdef OlD
    _vector_add_assign((*r).s1,psi);
    _vector_i_sub_assign((*r).s3,psi);
#else
    _vector_add((*r).s1, temp.s1, psi);
    _vector_i_sub((*r).s3, temp.s3, psi);
#endif
    /************************ end of loop ************************/
  }
}
/* end of If defined SSE2 */
#endif

static char const rcsid[] = "$Id$";

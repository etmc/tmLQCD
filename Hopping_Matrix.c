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
#  include "xchange_halffield.h"
#endif
#include "boundary.h"
#include "init_dirac_halfspinor.h"
#include "Hopping_Matrix.h"

void xchange_field(spinor * const l, const int even);


#if ((defined SSE2)||(defined SSE3))

/* input on k; output on l */
void Hopping_Matrix(const int ieo, spinor * const l, spinor * const k){
  int icx,icy,icz,ioff;
  int ix,iy,iz;
  su3 *up,*um;
  static spinor rs;
  spinor *s, *sp, *rn;
  halfspinor *r0, *r1, *r2, *r3, *rp, *rm;
  
  /* for parallelization */

  if(k == l){
    printf("Error in subroutine D_psi: improper arguments\n");
    printf("Program aborted\n");
    exit(1);
  }
  if(ieo == 1){
    ioff = 0;
  } 
  else{
    ioff = (VOLUME+RAND)/2;
  }

  ix=g_eo2lexic[ioff];
  iy=g_idn[ix][0]; 
  icy=g_lexic2eosub[iy];
  /* s contains the source vector */
  /* r0,r1,r2,r3 contain the intermediate half spinor */
  r0 = &halfs[0][ g_halfpt[icy][0] ];

#if ((defined _GAUGE_COPY))
  up=&g_gauge_field_copy[ioff][0];
#else
  up=&g_gauge_field[iy][0];
#endif

  /**************** loop over all lattice sites ******************/
  for(icx = ioff; icx < (VOLUME/2+ioff); icx++){
    s = k+(icx-ioff);
    ix=g_eo2lexic[icx];
    /*********************** direction +0 ************************/
    iy=g_idn[ix][1]; 
    icy=g_lexic2eosub[iy];

#if ((defined _GAUGE_COPY))
    um=up+1;
#else
    um=&g_gauge_field[iy][1]; 
#endif
    _prefetch_su3(um);

    _sse_load((*s).s0);
    _sse_load_up((*s).s2);
    _sse_vector_add();

    _sse_su3_multiply((*up));
    _sse_vector_cmplx_mul(ka0);
    _sse_store_nt_up((*r0).s0);

    r1 =&halfs[0][ g_halfpt[icy][1] ];      
/*     _prefetch_halfspinor(r1);  */

    _sse_load((*s).s1);
    _sse_load_up((*s).s3);
    _sse_vector_add();
      
    _sse_su3_multiply((*up));
    _sse_vector_cmplx_mul(ka0);
    _sse_store_nt_up((*r0).s1);

    /*********************** direction +1 ************************/

    iy=g_idn[ix][2]; 
    icy=g_lexic2eosub[iy];

#ifndef _GAUGE_COPY
    up=&g_gauge_field[iy][2]; 
#else
    up=um+1;
#endif
    _prefetch_su3(up);

    _sse_load((*s).s0);
    _sse_load_up((*s).s3);
    _sse_vector_i_mul();
    _sse_vector_add();

    _sse_su3_multiply((*um));
    _sse_vector_cmplx_mul(ka1);
    _sse_store_nt_up((*r1).s0);

    r2=&halfs[0][ g_halfpt[icy][2] ];
/*     _prefetch_halfspinor(r2); */

    _sse_load((*s).s1);
    _sse_load_up((*s).s2);
    _sse_vector_i_mul();
    _sse_vector_add();

    _sse_su3_multiply((*um));
    _sse_vector_cmplx_mul(ka1);
    _sse_store_nt_up((*r1).s1);

    /*********************** direction +2 ************************/

    iy=g_idn[ix][3]; 
    icy=g_lexic2eosub[iy];

#ifndef _GAUGE_COPY
    um=&g_gauge_field[iy][3];
#else
    um=up+1;
#endif
    _prefetch_su3(um);

    _sse_load((*s).s0);
    _sse_load_up((*s).s3);
    _sse_vector_add();

    _sse_su3_multiply((*up));
    _sse_vector_cmplx_mul(ka2);
    _sse_store_nt_up((*r2).s0);

    r3=&halfs[0][ g_halfpt[icy][3] ];
/*     _prefetch_halfspinor(r3);  */

    _sse_load((*s).s1);
    _sse_load_up((*s).s2);
    _sse_vector_sub();

    _sse_su3_multiply((*up));
    _sse_vector_cmplx_mul(ka2);
    _sse_store_nt_up((*r2).s1);

    /*********************** direction +3 ************************/
    icz=icx+1;
    if(icz==((VOLUME+RAND)/2+ioff)) icz=ioff;
    iz=g_eo2lexic[icz];
    iy=g_idn[iz][0]; icy=g_lexic2eosub[iy];

#ifndef _GAUGE_COPY
    up=&g_gauge_field[iy][0];
#else
    up=um+1;
#endif
    _prefetch_su3(up);
    _prefetch_spinor(s+1);

    _sse_load((*s).s0);
    _sse_load_up((*s).s2);
    _sse_vector_i_mul();
    _sse_vector_add();

    _sse_su3_multiply((*um));
    _sse_vector_cmplx_mul(ka3);
    _sse_store_nt_up((*r3).s0);

    r0=&halfs[0][ g_halfpt[icy][0] ];
/*     _prefetch_halfspinor(r0);  */

    _sse_load((*s).s1);
    _sse_load_up((*s).s3);
    _sse_vector_i_mul();
    _sse_vector_sub();

    _sse_su3_multiply((*um));
    _sse_vector_cmplx_mul(ka3);
    _sse_store_nt_up((*r3).s1);

  }

  ix=g_eo2lexic[ioff];
  iy=g_iup[ix][0]; 
  icy=g_lexic2eosub[iy];

  /* s contains the source vector */
  s = k;
  /* rp and rm contain the intermediate half spinor */
  r0 = &halfs[1][ g_halfpt[icy][0] ];

#if ((defined _GAUGE_COPY))
/*   up=g_gauge_field_copy[ioff]; */
#else
  up=g_gauge_field[ix];
#endif

  /* Here we need to do some exchange for plus direction */
#  ifdef MPI

#endif

  /* Now all the minus directions */
  for(icx = ioff; icx < (VOLUME/2+ioff); icx++){
    s = k+(icx-ioff);
    ix=g_eo2lexic[icx];
    /*********************** direction -0 ************************/

    iy=g_iup[ix][1]; 
    icy=g_lexic2eosub[iy];

    um = up + 1;
    _prefetch_su3(um);
      
    _sse_load((*s).s0);
    _sse_load_up((*s).s2);
    _sse_vector_sub();
      
    _sse_su3_inverse_multiply((*up));
    _sse_vector_cmplxcg_mul(ka0);
    _sse_store_nt_up((*r0).s0);
      
    r1=&halfs[1][ g_halfpt[icy][1] ];
/*     _prefetch_halfspinor(r1);  */

    _sse_load((*s).s1);
    _sse_load_up((*s).s3);
    _sse_vector_sub();
      
    _sse_su3_inverse_multiply((*up));
    _sse_vector_cmplxcg_mul(ka0);
    _sse_store_nt_up((*r0).s1);

    /*********************** direction -1 ************************/

    iy=g_iup[ix][2]; 
    icy=g_lexic2eosub[iy];

    up = um + 1;
    _prefetch_su3(up);

    _sse_load((*s).s0);
    _sse_load_up((*s).s3);
    _sse_vector_i_mul();
    _sse_vector_sub();
      
    _sse_su3_inverse_multiply((*um));
    _sse_vector_cmplxcg_mul(ka1);
    _sse_store_nt_up((*r1).s0);
      
    r2=&halfs[1][ g_halfpt[icy][2] ];
/*     _prefetch_halfspinor(r2);  */

    _sse_load((*s).s1);
    _sse_load_up((*s).s2);
    _sse_vector_i_mul();
    _sse_vector_sub();
      
    _sse_su3_inverse_multiply((*um));
    _sse_vector_cmplxcg_mul(ka1);
    _sse_store_nt_up((*r1).s1);

    /*********************** direction -2 ************************/

    iy=g_iup[ix][3]; 
    icy=g_lexic2eosub[iy];

    um = up + 1;
    _prefetch_su3(um);

    _sse_load((*s).s0);
    _sse_load_up((*s).s3);
    _sse_vector_sub();
      
    _sse_su3_inverse_multiply((*up));
    _sse_vector_cmplxcg_mul(ka2);
    _sse_store_nt_up((*r2).s0);
      
    r3=&halfs[1][ g_halfpt[icy][3] ];
/*     _prefetch_halfspinor(r3);  */

    _sse_load((*s).s1);
    _sse_load_up((*s).s2);
    _sse_vector_add();
      
    _sse_su3_inverse_multiply((*up));
    _sse_vector_cmplxcg_mul(ka2);
    _sse_store_nt_up((*r2).s1);
      
    /*********************** direction -3 ************************/
    icz=icx+1;
    if(icz==((VOLUME+RAND)/2+ioff)) icz=ioff;
    iz=g_eo2lexic[icz];
    iy=g_iup[iz][0]; 
    icy=g_lexic2eosub[iy];

#  if (defined _GAUGE_COPY)
    up=um+1;
#else
    up=&g_gauge_field[iz][0];
#endif

    _prefetch_su3(up); 
    _prefetch_spinor(s+1);

    _sse_load((*s).s0);
    _sse_load_up((*s).s2);
    _sse_vector_i_mul();
    _sse_vector_sub();
      
    _sse_su3_inverse_multiply((*um));
    _sse_vector_cmplxcg_mul(ka3);
    _sse_store_nt_up((*r3).s0);

    r0=&halfs[1][ g_halfpt[icy][0] ];
/*     _prefetch_halfspinor(r0);  */

    _sse_load((*s).s1);
    _sse_load_up((*s).s3);
    _sse_vector_i_mul();
    _sse_vector_add();
      
    _sse_su3_inverse_multiply((*um));
    _sse_vector_cmplxcg_mul(ka3);
    _sse_store_nt_up((*r3).s1);

  }
  rp = halfs[0];
  _prefetch_nta_halfspinor(rp);
  rm = halfs[1]; 
  _prefetch_nta_halfspinor(rm); 
#  ifdef MPI

#endif
  /* Now we sum up and expand to full spinor */
  for(icx = 0; icx < (VOLUME)/2; icx++) {
    /*********************** direction +0 and -0 *****************/
    _sse_load((*rp).s0);
    _sse_load_up((*rm).s0);
    _sse_vector_add();
    _sse_store(rs.s0);

    _sse_load((*rp).s0);
    _sse_vector_sub();
    _sse_store(rs.s2);

    _sse_load_up((*rm).s1);
    _sse_load((*rp).s1);
    _sse_vector_add();
    _sse_store(rs.s1);

    _sse_load((*rp).s1);
    _sse_vector_sub();
    _sse_store(rs.s3);
    _prefetch_nta_halfspinor(rp+4);
    _prefetch_nta_halfspinor(rm+4);
    rp++;
    rm++;
    /*********************** direction +1 ************************/
/*     _prefetch_nta_halfspinor(rp+1); */
    _sse_load_up((*rp).s0);
    _sse_load(rs.s0);
    _sse_vector_add();
    _sse_store(rs.s0);

    _sse_load(rs.s3);
    _sse_vector_i_mul();      
    _sse_vector_sub();
    _sse_store(rs.s3); 

    _sse_load_up((*rp).s1);
    _sse_load(rs.s1);
    _sse_vector_add();
    _sse_store(rs.s1);

    _sse_load(rs.s2);
    _sse_vector_i_mul();      
    _sse_vector_sub();
    _sse_store(rs.s2);
    _prefetch_nta_halfspinor(rp+4);
    rp++;
    /*********************** direction -1 ************************/
/*     _prefetch_nta_halfspinor(rm+1); */
    _sse_load_up((*rm).s0);
    _sse_load(rs.s0);
    _sse_vector_add();
    _sse_store(rs.s0);

    _sse_load(rs.s3);
    _sse_vector_i_mul();      
    _sse_vector_add();
    _sse_store(rs.s3);
/*     _prefetch_spinor(rp+1); */
    _sse_load_up((*rm).s1);
    _sse_load(rs.s1);
    _sse_vector_add();
    _sse_store(rs.s1);

    _sse_load(rs.s2);
    _sse_vector_i_mul();      
    _sse_vector_add();
    _sse_store(rs.s2);
    _prefetch_nta_halfspinor(rm+4);
    rm++;
    /*********************** direction +4 ************************/
/*     _prefetch_nta_halfspinor(rp+1); */

    _sse_load_up((*rp).s0);
    _sse_load(rs.s0);
    _sse_vector_add();
    _sse_store(rs.s0);

    _sse_load(rs.s3);
    _sse_vector_add();
    _sse_store(rs.s3);

    _sse_load_up((*rp).s1);
    _sse_load(rs.s1);
    _sse_vector_add();
    _sse_store(rs.s1);

    _sse_load(rs.s2);
    _sse_vector_sub();
    _sse_store(rs.s2);
    _prefetch_nta_halfspinor(rp+4);
    rp++; 
    /*********************** direction -2 ************************/
/*     _prefetch_nta_halfspinor(rm+1); */
    _sse_load_up((*rm).s0);
    _sse_load(rs.s0);
    _sse_vector_add();
    _sse_store(rs.s0);

    _sse_load(rs.s3);
    _sse_vector_sub();
    _sse_store(rs.s3);

    _sse_load_up((*rm).s1);
    _sse_load(rs.s1);
    _sse_vector_add();
    _sse_store(rs.s1);

    _sse_load(rs.s2);
    _sse_vector_add();
    _sse_store(rs.s2);
    _prefetch_nta_halfspinor(rm+4);
    rm++;
    /*********************** direction +4 ************************/

    _sse_load_up((*rp).s0);
    _sse_load(rs.s0);
    _sse_vector_add();
    _sse_store(rs.s0);

    _sse_load(rs.s2);
    _sse_vector_i_mul();      
    _sse_vector_sub();
    _sse_store(rs.s2);

    rn = l + icx;
    _prefetch_spinor(rn);  

/*     _prefetch_nta_halfspinor(rp+1);  */

    _sse_load_up((*rp).s1);
    _sse_load(rs.s1);
    _sse_vector_add();
    _sse_store(rs.s1);

    _sse_load(rs.s3);
    _sse_vector_i_mul();      
    _sse_vector_add();
    _sse_store(rs.s3);
    _prefetch_nta_halfspinor(rp+4);
    rp++;
    /*********************** direction -3 ************************/
/*     _prefetch_nta_halfspinor(rm+1); */
    _sse_load_up((*rm).s0);
    _sse_load(rs.s0);
    _sse_vector_add();
    _sse_store((*rn).s0);

    _sse_load(rs.s2);
    _sse_vector_i_mul();      
    _sse_vector_add();
    _sse_store((*rn).s2);

    _sse_load_up((*rm).s1);
    _sse_load(rs.s1);
    _sse_vector_add();
    _sse_store((*rn).s1);

    _sse_load(rs.s3);
    _sse_vector_i_mul();      
    _sse_vector_sub();
    _sse_store((*rn).s3);
    _prefetch_nta_halfspinor(rm+4);
    rm++;
  }
}

#elif (defined BGL && defined XLC)

/**********************************
 *
 * Blue Gene/L Version
 *
 * Author: Carsten Urbach
 *
 **********************************/


void Hopping_Matrix(const int ieo_, spinor * const l, spinor * const k){
  int icx,icy,icz,ieo;
  int i;
  int ix,iy,iz;
  su3 *U ALIGN;
  spinor *s ALIGN;
  halfspinor ** phi ALIGN;
  halfspinor *r ALIGN;
  /* We have 32 registers available */
  double _Complex reg00, reg01, reg02, reg03, reg04, reg05;
  double _Complex reg10, reg11, reg12, reg13, reg14, reg15;
  /* For the gauge field, reuse the first three!*/
  double _Complex u00, u01, u02, u10, u11, u12;
  double _Complex reg20, reg21;
  /* The following contains the result spinor (12 regs) */
  double _Complex rs00, rs01, rs02, rs10, rs11, rs12, rs20, rs21, rs22, 
    rs30, rs31, rs32;

#pragma disjoint(*s, *U, *l, *k, *r, *s)

  __alignx(16, l);
  __alignx(16, k);
  __alignx(16, HalfSpinor);
  if(ieo_ == 0) ieo = 0;
  else ieo = 1;
  /* We will run through the source vector now */
  /* instead of the solution vector            */
  s = k;
  _prefetch_spinor(s); 

  /* s contains the source vector */
  /* r0,r1,r2,r3 contain the intermediate half spinor */

  if(ieo == 0) {
    U = g_gauge_field_copy[0][0][0];
  }
  else {
    U = g_gauge_field_copy[0][1][0];
  }
  phi = NBPointer[ieo];

  _prefetch_su3(U);
  /**************** loop over all lattice sites ******************/
  ix=0;
  for(i = 0; i < (VOLUME)/2; i++){

    _bgl_load_rs0((*s).s0);
    _bgl_load_rs1((*s).s1);
    _bgl_load_rs2((*s).s2);
    _bgl_load_rs3((*s).s3);
    s++; 
    _prefetch_spinor(s); 
    /*********************** direction +0 ************************/
    _prefetch_su3(U+1);

    _bgl_vector_add_rs2_to_rs0_reg0();
    _bgl_vector_add_rs3_to_rs1_reg1();

    _bgl_su3_multiply_double((*U));
    _bgl_vector_cmplx_mul_double(ka0);
    /* result is now in regx0, regx1, regx2 , x=0,1 */

    _bgl_store_reg0_up((*phi[ix]).s0);
    _bgl_store_reg1_up((*phi[ix]).s1);
    U++;
    ix++;

    /*********************** direction -0 ************************/
    _bgl_vector_sub_rs2_from_rs0_reg0();
    _bgl_vector_sub_rs3_from_rs1_reg1();

    _bgl_store_reg0((*phi[ix]).s0);
    _bgl_store_reg1((*phi[ix]).s1);
    ix++;

    /*********************** direction +1 ************************/
    _prefetch_su3(U+1);
    _bgl_vector_i_mul_add_rs3_to_rs0_reg0();
    _bgl_vector_i_mul_add_rs2_to_rs1_reg1();

    _bgl_su3_multiply_double((*U));
    _bgl_vector_cmplx_mul_double(ka1);

    _bgl_store_reg0_up((*phi[ix]).s0);
    _bgl_store_reg1_up((*phi[ix]).s1);
    ix++;
    U++;

    /*********************** direction -1 ************************/
    _bgl_vector_i_mul_sub_rs3_from_rs0_reg0();
    _bgl_vector_i_mul_sub_rs2_from_rs1_reg1();

    _bgl_store_reg0((*phi[ix]).s0);
    _bgl_store_reg1((*phi[ix]).s1);
    ix++;


    /*********************** direction +2 ************************/
    _prefetch_su3(U+1);

    _bgl_vector_add_rs3_to_rs0_reg0();
    _bgl_vector_sub_rs2_from_rs1_reg1();

    _bgl_su3_multiply_double((*U));
    _bgl_vector_cmplx_mul_double(ka2);

    _bgl_store_reg0_up((*phi[ix]).s0);
    _bgl_store_reg1_up((*phi[ix]).s1);
    ix++;
    U++;

    /*********************** direction -2 ************************/
    _bgl_vector_sub_rs3_from_rs0_reg0();
    _bgl_vector_add_rs2_to_rs1_reg1();

    _bgl_store_reg0((*phi[ix]).s0);
    _bgl_store_reg1((*phi[ix]).s1);
    ix++;

    /*********************** direction +3 ************************/
    _prefetch_su3(U+1); 

    _bgl_vector_i_mul_add_rs2_to_rs0_reg0();
    _bgl_vector_i_mul_sub_rs3_from_rs1_reg1();

    _bgl_su3_multiply_double((*U));
    _bgl_vector_cmplx_mul_double(ka3);

    _bgl_store_reg0_up((*phi[ix]).s0);
    _bgl_store_reg1_up((*phi[ix]).s1);
    ix++;
    U++;

    /*********************** direction -3 ************************/
    _bgl_vector_i_mul_sub_rs2_from_rs0_reg0();
    _bgl_vector_i_mul_add_rs3_to_rs1_reg1();

    _bgl_store_reg0((*phi[ix]).s0);
    _bgl_store_reg1((*phi[ix]).s1);
    ix++;

    /************************ end of loop ************************/
  }

#  if (defined MPI && !defined _NO_COMM)
  xchange_halffield(ieo); 
#  endif
  s = l;
  phi = NBPointer[2 + ieo];
  if(ieo == 0) {
    U = g_gauge_field_copy[0][1][0];
  }
  else {
    U = g_gauge_field_copy[0][0][0];
  }
  _prefetch_halfspinor(phi[0]);
  _prefetch_su3(U);
  
  /* Now we sum up and expand to a full spinor */
  ix = 0;
  _prefetch_spinor_for_store(s); 
  for(i = 0; i < (VOLUME)/2; i++){
    _prefetch_spinor_for_store(s);
    _prefetch_halfspinor(phi[ix+1]);
    /*********************** direction +0 ************************/
    _bgl_load_rs0((*phi[ix]).s0);
    rs20 = rs00;
    rs21 = rs01;
    rs22 = rs02;
    _bgl_load_rs1((*phi[ix]).s1);
    rs30 = rs10;
    rs31 = rs11;
    rs32 = rs12;
    if(g_proc_id == -1) {
      printf("%d dir +0 %e %e\n", i, __creal(rs00), __cimag(rs00));
      printf("%d dir +0 %e %e\n", i, __creal(rs10), __cimag(rs10));
      printf("%d dir +0 %e %e\n", i, __creal(rs20), __cimag(rs20));
      printf("%d dir +0 %e %e\n", i, __creal(rs30), __cimag(rs30));
    }
    ix++;
    /*********************** direction -0 ************************/
    _prefetch_su3(U+1);

    _bgl_load_reg0((*phi[ix]).s0);
    _bgl_load_reg1((*phi[ix]).s1);

    _bgl_su3_inverse_multiply_double((*U));
    _bgl_vector_cmplxcg_mul_double(ka0);
    /* result is in the upper parts of reg0? and reg1? */

    _bgl_add_to_rs0_reg0();
    _bgl_sub_from_rs2_reg0();
    _bgl_add_to_rs1_reg1();
    _bgl_sub_from_rs3_reg1();
    if(g_proc_id == -1) {
      printf("%d dir -0 %e %e\n", i, __creal(reg03), __cimag(reg03));
      printf("%d dir -0 %e %e\n", i, __creal(rs00), __cimag(rs00));
      printf("%d dir +0 %e %e\n", i, __creal(rs10), __cimag(rs10));
      printf("%d dir +0 %e %e\n", i, __creal(rs20), __cimag(rs20));
      printf("%d dir +0 %e %e\n", i, __creal(rs30), __cimag(rs30));
    }
    U++;
    ix++;
    /*********************** direction +1 ************************/
    _bgl_load_reg0_up((*phi[ix]).s0);
    _bgl_load_reg1_up((*phi[ix]).s1);

    _bgl_add_to_rs0_reg0();
    _bgl_i_mul_sub_from_rs3_reg0();
    _bgl_add_to_rs1_reg1();
    _bgl_i_mul_sub_from_rs2_reg1();
    if(g_proc_id == -1) {
      printf("%d dir +1 %e %e\n", i, __creal(reg03), __cimag(reg03));
      printf("%d dir +1 %e %e\n", i, __creal(rs00), __cimag(rs00));
      printf("%d dir +0 %e %e\n", i, __creal(rs10), __cimag(rs10));
      printf("%d dir +0 %e %e\n", i, __creal(rs20), __cimag(rs20));
      printf("%d dir +0 %e %e\n", i, __creal(rs30), __cimag(rs30));
    }
    ix++;
    /*********************** direction -1 ************************/
    _prefetch_su3(U+1);

    _bgl_load_reg0((*phi[ix]).s0);
    _bgl_load_reg1((*phi[ix]).s1);

    _bgl_su3_inverse_multiply_double((*U));
    _bgl_vector_cmplxcg_mul_double(ka1);

    _bgl_add_to_rs0_reg0();
    _bgl_add_to_rs1_reg1();
    _bgl_i_mul_add_to_rs3_reg0();
    _bgl_i_mul_add_to_rs2_reg1();      
    if(g_proc_id == -1) {
      printf("%d dir -1 %e %e\n", i, __creal(reg03), __cimag(reg03));
      printf("%d dir -1 %e %e\n", i, __creal(rs00), __cimag(rs00));
      printf("%d dir +0 %e %e\n", i, __creal(rs10), __cimag(rs10));
      printf("%d dir +0 %e %e\n", i, __creal(rs20), __cimag(rs20));
      printf("%d dir +0 %e %e\n", i, __creal(rs30), __cimag(rs30));
    }
    U++;
    ix++;
    /*********************** direction +2 ************************/
    _bgl_load_reg0_up((*phi[ix]).s0);
    _bgl_load_reg1_up((*phi[ix]).s1);

    _bgl_add_to_rs0_reg0();
    _bgl_add_to_rs1_reg1();
    _bgl_sub_from_rs2_reg1();
    _bgl_add_to_rs3_reg0();
    if(g_proc_id == -1) {
      printf("%d dir +2 %e %e\n", i, __creal(reg03), __cimag(reg03));
      printf("%d dir +2 %e %e\n", i, __creal(rs00), __cimag(rs00));
      printf("%d dir +0 %e %e\n", i, __creal(rs10), __cimag(rs10));
      printf("%d dir +0 %e %e\n", i, __creal(rs20), __cimag(rs20));
      printf("%d dir +0 %e %e\n", i, __creal(rs30), __cimag(rs30));
    }
    ix++;
    /*********************** direction -2 ************************/
    _prefetch_su3(U+1);

    _bgl_load_reg0((*phi[ix]).s0);
    _bgl_load_reg1((*phi[ix]).s1);

    _bgl_su3_inverse_multiply_double((*U));
    _bgl_vector_cmplxcg_mul_double(ka2);

    _bgl_add_to_rs0_reg0();
    _bgl_add_to_rs1_reg1();
    _bgl_add_to_rs2_reg1();
    _bgl_sub_from_rs3_reg0();
    if(g_proc_id == -1) {
      printf("%d dir -2 %e %e\n", i, __creal(reg03), __cimag(reg03));
      printf("%d dir -2 %e %e\n", i, __creal(rs00), __cimag(rs00));
      printf("%d dir +0 %e %e\n", i, __creal(rs10), __cimag(rs10));
      printf("%d dir +0 %e %e\n", i, __creal(rs20), __cimag(rs20));
      printf("%d dir +0 %e %e\n", i, __creal(rs30), __cimag(rs30));
    }
    U++;
    ix++;
    /*********************** direction +3 ************************/
    _bgl_load_reg0_up((*phi[ix]).s0);
    _bgl_load_reg1_up((*phi[ix]).s1);

    _bgl_add_to_rs0_reg0();
    _bgl_add_to_rs1_reg1();
    _bgl_i_mul_sub_from_rs2_reg0();
    _bgl_i_mul_add_to_rs3_reg1();
    if(g_proc_id == -1) {
      printf("%d dir +3 %e %e\n", i, __creal(reg03), __cimag(reg03));
      printf("%d dir +3 %e %e\n", i, __creal(rs00), __cimag(rs00));
      printf("%d dir +0 %e %e\n", i, __creal(rs10), __cimag(rs10));
      printf("%d dir +0 %e %e\n", i, __creal(rs20), __cimag(rs20));
      printf("%d dir +0 %e %e\n", i, __creal(rs30), __cimag(rs30));
    }
    ix++;
    /*********************** direction -3 ************************/
    _prefetch_su3(U+1);

    _bgl_load_reg0((*phi[ix]).s0);
    _bgl_load_reg1((*phi[ix]).s1);

    _bgl_su3_inverse_multiply_double((*U));
    _bgl_vector_cmplxcg_mul_double(ka3);

    _bgl_add_to_rs0_reg0();
    _bgl_store_rs0((*s).s0);
    _bgl_i_mul_add_to_rs2_reg0();
    _bgl_store_rs2((*s).s2);

    _bgl_add_to_rs1_reg1();
    _bgl_store_rs1((*s).s1);
    _bgl_i_mul_sub_from_rs3_reg1();
    _bgl_store_rs3((*s).s3);
    if(g_proc_id == -1) {
      printf("%d dir -3 %e %e\n", i, __creal(reg03), __cimag(reg03));
      printf("%d dir -3 %e %e\n", i, __creal(rs00), __cimag(rs00));
      printf("%d dir +0 %e %e\n", i, __creal(rs10), __cimag(rs10));
      printf("%d dir +0 %e %e\n", i, __creal(rs20), __cimag(rs20));
      printf("%d dir +0 %e %e\n", i, __creal(rs30), __cimag(rs30));
      printf("%e %e %e %e \n", (*s).s0.c0.re, (*s).s0.c0.im, (*s).s0.c1.re, (*s).s0.c1.im);
      printf("%e %e %e %e \n", (*s).s0.c2.re, (*s).s0.c2.im, (*s).s1.c0.re, (*s).s1.c0.im);
      printf("%e %e %e %e \n", (*s).s1.c1.re, (*s).s1.c1.im, (*s).s1.c2.re, (*s).s1.c2.im);
      printf("%e %e %e %e \n", (*s).s2.c0.re, (*s).s2.c0.im, (*s).s2.c1.re, (*s).s2.c1.im);
      printf("%e %e %e %e \n", (*s).s2.c2.re, (*s).s2.c2.im, (*s).s3.c0.re, (*s).s3.c0.im);
      printf("%e %e %e %e \n", (*s).s3.c1.re, (*s).s3.c1.im, (*s).s3.c2.re, (*s).s3.c2.im);
    }
    U++;
    ix++;
    s++;
  }
/*   if(ieo_ == 1) { */
/*     exit(0); */
/*   } */
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
  spinor temp;
#pragma disjoint(temp, *sp, *sm, *r, *up, *um)

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
#if ((defined _GAUGE_COPY))
  up=&g_gauge_field_copy[ioff][0];
#else
  up=&g_gauge_field[ix][0];
#endif
  for (icx = ioff; icx < (VOLUME/2 + ioff); icx++) {
    ix = g_eo2lexic[icx];
    /*********************** direction +0 ************************/
    
    iy=g_idn[ix][0]; icy=g_lexic2eosub[iy];
    
    sm=k+icy;
    _prefetch_spinor((void*)sm);

#ifndef _GAUGE_COPY
    um=&g_gauge_field[iy][0];
#else
    um=up+1;
#endif
    _prefetch_su3((void*)um);

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
    _prefetch_spinor((void*)sp);
#if ((defined _GAUGE_COPY))
    up=um+1;
#else
    up+=1;
#endif
    _prefetch_su3((void*)up); 
      
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
    _prefetch_spinor((void*)sm);

#ifndef _GAUGE_COPY
    um=&g_gauge_field[iy][1];
#else
    um=up+1;
#endif
    _prefetch_su3((void*)um);
      
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
    _prefetch_spinor((void*)sp);
#if ((defined _GAUGE_COPY))
    up=um+1;
#else
    up+=1;
#endif
    _prefetch_su3((void*)up);

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
    _prefetch_spinor((void*)sm);

#ifndef _GAUGE_COPY
    um=&g_gauge_field[iy][2];
#else
    um=up+1;
#endif
    _prefetch_su3((void*)um);

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
    _prefetch_spinor((void*)sp);
#if ((defined _GAUGE_COPY))
    up=um+1;
#else
    up+=1;
#endif
    _prefetch_su3((void*)up);

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
    _prefetch_spinor((void*)sm);

#ifndef _GAUGE_COPY
    um=&g_gauge_field[iy][3];
#else
    um=up+1;
#endif
    _prefetch_su3((void*)um);

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
    _prefetch_spinor((void*)r);

    icz = icx + 1;
    if(icz==((VOLUME+RAND)/2+ioff)) icz=ioff;
    iz=g_eo2lexic[icz];
    iy=g_iup[iz][0]; icy=g_lexic2eosub[iy];
#if ((defined _GAUGE_COPY))
    up=&g_gauge_field_copy[icz][0];
#else
    up=&g_gauge_field[iz][0];
#endif
    _prefetch_su3((void*)up);

    sp=k+icy;
    _prefetch_spinor((void*)sp);

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

/* l output , k input*/
/* for ieo=0, k resides on  odd sites and l on even sites */
void Hopping_Matrix(const int ieo, spinor * const l, spinor * const k){
  int i,ix;
  su3 * U ALIGN;
  spinor *s ALIGN;
  spinor rs;
  su3_vector psi, chi;
  halfspinor ** phi ALIGN;


  if(k == l){
    printf("Error in H_psi (simple.c):\n");
    printf("Arguments k and l must be different\n");
    printf("Program aborted\n");
    exit(1);
  }
  s = k;

  if(ieo == 0) {
    U = g_gauge_field_copy[0][0][0];
  }
  else {
    U = g_gauge_field_copy[0][1][0];
  }
  phi = NBPointer[ieo];

  /**************** loop over all lattice sites ****************/
  ix=0;
  for(i = 0; i < (VOLUME)/2; i++){
    _vector_assign(rs.s0, (*s).s0);
    _vector_assign(rs.s1, (*s).s1);
    _vector_assign(rs.s2, (*s).s2);
    _vector_assign(rs.s3, (*s).s3);
    s++;
    /*********************** direction +0 ************************/
      
    _vector_add(psi, rs.s0, rs.s2);

    _su3_multiply(chi,(*U),psi);
    _complex_times_vector((*phi[ix]).s0, ka0, chi);
      
    _vector_add(psi, rs.s1, rs.s3);

    _su3_multiply(chi,(*U),psi);
    _complex_times_vector((*phi[ix]).s1, ka0, chi);
            
    U++;
    ix++;
    
    /*********************** direction -0 ************************/

    _vector_sub((*phi[ix]).s0, rs.s0, rs.s2);
    _vector_sub((*phi[ix]).s1, rs.s1, rs.s3);

    ix++;

    /*********************** direction +1 ************************/

    _vector_i_add(psi, rs.s0, rs.s3);

    _su3_multiply(chi, (*U), psi);
    _complex_times_vector((*phi[ix]).s0, ka1, chi);

    _vector_i_add(psi, rs.s1, rs.s2);

    _su3_multiply(chi, (*U), psi);
    _complex_times_vector((*phi[ix]).s1, ka1, chi);

    U++;
    ix++;

    /*********************** direction -1 ************************/

    _vector_i_sub((*phi[ix]).s0, rs.s0, rs.s3);
    _vector_i_sub((*phi[ix]).s1, rs.s1, rs.s2);

    ix++;
    /*********************** direction +2 ************************/

    _vector_add(psi, rs.s0, rs.s3);

    _su3_multiply(chi,(*U),psi);
    _complex_times_vector((*phi[ix]).s0, ka2, chi);

    _vector_sub(psi, rs.s1, rs.s2);

    _su3_multiply(chi,(*U),psi);
    _complex_times_vector((*phi[ix]).s1, ka2, chi);
      
    U++;
    ix++;

    /*********************** direction -2 ************************/

    _vector_sub((*phi[ix]).s0, rs.s0, rs.s3);
    _vector_add((*phi[ix]).s1, rs.s1, rs.s2);
    ix++;

    /*********************** direction +3 ************************/

    _vector_i_add(psi, rs.s0, rs.s2);
      
    _su3_multiply(chi, (*U), psi);
    _complex_times_vector((*phi[ix]).s0, ka3, chi);


    _vector_i_sub(psi, rs.s1, rs.s3);

    _su3_multiply(chi,(*U),psi);
    _complex_times_vector((*phi[ix]).s1, ka3, chi);

    U++;
    ix++;
    /*********************** direction -3 ************************/

    _vector_i_sub((*phi[ix]).s0, rs.s0, rs.s2);
    _vector_i_add((*phi[ix]).s1, rs.s1, rs.s3);

    ix++;
    /************************ end of loop ************************/
  }
#  if (defined MPI && !defined _NO_COMM)
  xchange_halffield(ieo); 
#  endif
  s = l;
  phi = NBPointer[2 + ieo];
  if(ieo == 0) {
    U = g_gauge_field_copy[0][1][0];
  }
  else {
    U = g_gauge_field_copy[0][0][0];
  }

  ix = 0;
  for(i = 0; i < (VOLUME)/2; i++){
    /*********************** direction +0 ************************/
    _vector_assign(rs.s0, (*phi[ix]).s0);
    _vector_assign(rs.s2, (*phi[ix]).s0);
    _vector_assign(rs.s1, (*phi[ix]).s1);
    _vector_assign(rs.s3, (*phi[ix]).s1);
    ix++;
    /*********************** direction -0 ************************/
    _su3_inverse_multiply(chi,(*U),(*phi[ix]).s0);
    _complexcjg_times_vector(psi,ka0,chi);

    _vector_add_assign(rs.s0, psi);
    _vector_sub_assign(rs.s2, psi);

    _su3_inverse_multiply(chi,(*U),(*phi[ix]).s1);
    _complexcjg_times_vector(psi,ka0,chi);
      
    _vector_add_assign(rs.s1, psi);
    _vector_sub_assign(rs.s3, psi);
    ix++;
    U++;
    /*********************** direction +1 ************************/

    _vector_add_assign(rs.s0, (*phi[ix]).s0);
    _vector_i_sub_assign(rs.s3, (*phi[ix]).s0);

    _vector_add_assign(rs.s1, (*phi[ix]).s1);
    _vector_i_sub_assign(rs.s2, (*phi[ix]).s1);
    
    ix++;
    /*********************** direction -1 ************************/

    _su3_inverse_multiply(chi,(*U), (*phi[ix]).s0);
    _complexcjg_times_vector(psi,ka1,chi);

    _vector_add_assign(rs.s0, psi);
    _vector_i_add_assign(rs.s3, psi);

    _su3_inverse_multiply(chi,(*U), (*phi[ix]).s1);
    _complexcjg_times_vector(psi,ka1,chi);

    _vector_add_assign(rs.s1, psi);
    _vector_i_add_assign(rs.s2, psi);

    U++;
    ix++;

    /*********************** direction +2 ************************/

    _vector_add_assign(rs.s0, (*phi[ix]).s0);
    _vector_add_assign(rs.s3, (*phi[ix]).s0);

    _vector_add_assign(rs.s1, (*phi[ix]).s1);
    _vector_sub_assign(rs.s2, (*phi[ix]).s1);
    
    ix++;
    /*********************** direction -2 ************************/

    _su3_inverse_multiply(chi,(*U), (*phi[ix]).s0);
    _complexcjg_times_vector(psi,ka2,chi);

    _vector_add_assign(rs.s0, psi);
    _vector_sub_assign(rs.s3, psi);

    _su3_inverse_multiply(chi, (*U), (*phi[ix]).s1);
    _complexcjg_times_vector(psi,ka2,chi);
      
    _vector_add_assign(rs.s1, psi);
    _vector_add_assign(rs.s2, psi);

    U++;
    ix++;
    /*********************** direction +3 ************************/

    _vector_add_assign(rs.s0, (*phi[ix]).s0);
    _vector_i_sub_assign(rs.s2, (*phi[ix]).s0);

    _vector_add_assign(rs.s1, (*phi[ix]).s1);
    _vector_i_add_assign(rs.s3, (*phi[ix]).s1);

    ix++;

    /*********************** direction -3 ************************/

    _su3_inverse_multiply(chi,(*U), (*phi[ix]).s0);
    _complexcjg_times_vector(psi,ka3,chi);
      
    _vector_add((*s).s0, rs.s0, psi);
    _vector_i_add((*s).s2, rs.s2, psi);

    _su3_inverse_multiply(chi,(*U), (*phi[ix]).s1);
    _complexcjg_times_vector(psi,ka3,chi);

    _vector_add((*s).s1, rs.s1, psi);
    _vector_i_sub((*s).s3, rs.s3, psi);

    U++;
    ix++;
    s++;
  }
}
/* end of If defined SSE2 */
#endif

static char const rcsid[] = "$Id$";

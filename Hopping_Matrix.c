/**********************************************************************
 *
 *
 * Copyright (C) 2001 Martin Luescher
 *               2002 Martin Hasenbusch
 *               2003, 2004, 2005, 2006, 2007, 2008 Carsten Urbach
 *
 * BG and halfspinor versions (C) 2007, 2008 Carsten Urbach
 *
 * This file is based on an implementation of the Dirac operator 
 * written by Martin Luescher, modified by Martin Hasenbusch in 2002 
 * and modified and extended by Carsten Urbach from 2003-2008
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
 *
 * Hopping_Matrix is the conventional Wilson 
 * hopping matrix
 *
 * \kappa\sum_{\pm\mu}(r+\gamma_\mu)U_{x,\mu}
 *
 * for ieo = 0 this is M_{eo}, for ieo = 1
 * it is M_{oe}
 *
 * l is the output, k the input field
 *
 *  Structure of top level precompiler directives 
 *
 * - defining _USE_HALFSPINOR implies that we also use
 *   a "gauge copy"
 *
 * - such that we are checking for the _USE_GAUGECOPY feature seperatly in the 
 *   ELSE branch of the "if defined _USE_HALFSPINOR" statement
 *
 * - the numbers represent the certain implementation of Hopping_Matrix
 *
 * #if defined _USE_HALFSPINOR
 *   #if ((defined SSE2)||(defined SSE3))
 *     1. 
 *    #elif (defined BGL && defined XLC)
 *     2.
 *   #else
 *     3.
 *   #endif
 * #else * thats _USE_HALFSPINOR *
 *   #if ((defined SSE2)||(defined SSE3))
 *     #if (defined _USE_TSPLITPAR)
 *     4.
 *     #else
 *     5.
 *     #endif
 *   #elif (defined BGL && defined XLC)
 *     6.
 *   #elif defined XLC
 *     7.
 *   * else of If defined SSE2  and if defined XLC *
 *   #else
 *     8.
 *   #endif
 * #endif * thats _USE_HALFSPINOR *
 *
 ****************************************************************/

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
#  ifdef _USE_TSPLITPAR
#    include "xchange_field_tslice.h"
#  endif
#  if defined _USE_HALFSPINOR
#    include "xchange_halffield.h"
#  endif
#endif
#include "boundary.h"
#include "init_dirac_halfspinor.h"
#include "update_backward_gauge.h"
#include "Hopping_Matrix.h"

#if defined _USE_HALFSPINOR
#  if ((defined SSE2)||(defined SSE3))

/* 1. */
/* input on k; output on l */
void Hopping_Matrix(const int ieo, spinor * const l, spinor * const k){
  int ix, i;
  su3 * restrict U ALIGN;
  static spinor rs;
  spinor * restrict s ALIGN;
  halfspinor ** phi ALIGN;
#if defined OPTERON
  const int predist=2;
#else
  const int predist=1;
#endif
#ifdef _KOJAK_INST
#pragma pomp inst begin(hoppingmatrix)
#endif

#ifdef _GAUGE_COPY
  if(g_update_gauge_copy) {
    update_backward_gauge();
  }
#endif
  /* We will run through the source vector now */
  /* instead of the solution vector            */
  s = k;
  _prefetch_spinor(s);

  if(ieo == 0) {
    U = g_gauge_field_copy[0][0];
  }
  else {
    U = g_gauge_field_copy[1][0];
  }
  phi = NBPointer[ieo];

  _prefetch_su3(U);
  /**************** loop over all lattice sites ******************/
  ix=0;
  for(i = 0; i < (VOLUME)/2; i++){

    /*********************** direction +0 ************************/
    _prefetch_su3(U+predist);

    _sse_load((*s).s0);
    _sse_load_up((*s).s2);
    _sse_vector_add();

    _sse_su3_multiply((*U));
    _sse_vector_cmplx_mul(ka0);
    _sse_store_nt_up((*phi[ix]).s0);

    _sse_load((*s).s1);
    _sse_load_up((*s).s3);
    _sse_vector_add();
      
    _sse_su3_multiply((*U));
    _sse_vector_cmplx_mul(ka0);
    _sse_store_nt_up((*phi[ix]).s1);
    U++;
    ix++;
    /*********************** direction -0 ************************/
    _sse_load((*s).s0);
    _sse_load_up((*s).s2);
    _sse_vector_sub();
    _sse_store_nt((*phi[ix]).s0);

    _sse_load((*s).s1);
    _sse_load_up((*s).s3);
    _sse_vector_sub();
    _sse_store_nt((*phi[ix]).s1);
    ix++;

    /*********************** direction +1 ************************/
    _prefetch_su3(U+predist);

    _sse_load((*s).s0);
    /*next not needed?*/
    _sse_load_up((*s).s3);
    _sse_vector_i_mul();
    _sse_vector_add();

    _sse_su3_multiply((*U));
    _sse_vector_cmplx_mul(ka1);
    _sse_store_nt_up((*phi[ix]).s0);

    _sse_load((*s).s1);
    _sse_load_up((*s).s2);
    _sse_vector_i_mul();
    _sse_vector_add();

    _sse_su3_multiply((*U));
    _sse_vector_cmplx_mul(ka1);
    _sse_store_nt_up((*phi[ix]).s1);
    ix++;
    U++;

    /*********************** direction -1 ************************/
    _sse_load((*s).s0);
    _sse_load_up((*s).s3);
    _sse_vector_i_mul();
    _sse_vector_sub();
    _sse_store_nt((*phi[ix]).s0);

    _sse_load((*s).s1);
    _sse_load_up((*s).s2);
    _sse_vector_i_mul();
    _sse_vector_sub();
    _sse_store_nt((*phi[ix]).s1);
    ix++;

    /*********************** direction +2 ************************/
    _prefetch_su3(U+predist);

    _sse_load((*s).s0);
    _sse_load_up((*s).s3);
    _sse_vector_add();

    _sse_su3_multiply((*U));
    _sse_vector_cmplx_mul(ka2);
    _sse_store_nt_up((*phi[ix]).s0);

    _sse_load((*s).s1);
    _sse_load_up((*s).s2);
    _sse_vector_sub();

    _sse_su3_multiply((*U));
    _sse_vector_cmplx_mul(ka2);
    _sse_store_nt_up((*phi[ix]).s1);
    ix++;
    U++;
    /*********************** direction -2 ************************/
    _sse_load((*s).s0);
    _sse_load_up((*s).s3);
    _sse_vector_sub();
    _sse_store_nt((*phi[ix]).s0);

    _sse_load((*s).s1);
    _sse_load_up((*s).s2);
    _sse_vector_add();
    _sse_store_nt((*phi[ix]).s1);
    ix++;

    /*********************** direction +3 ************************/
    _prefetch_su3(U+predist);
    _prefetch_spinor(s+1);

    _sse_load((*s).s0);
    _sse_load_up((*s).s2);
    _sse_vector_i_mul();
    _sse_vector_add();

    _sse_su3_multiply((*U));
    _sse_vector_cmplx_mul(ka3);
    _sse_store_nt_up((*phi[ix]).s0);

    _sse_load((*s).s1);
    _sse_load_up((*s).s3);
    _sse_vector_i_mul();
    _sse_vector_sub();

    _sse_su3_multiply((*U));
    _sse_vector_cmplx_mul(ka3);
    _sse_store_nt_up((*phi[ix]).s1);
    ix++;
    U++;

    /*********************** direction -3 ************************/
    _sse_load((*s).s0);
    _sse_load_up((*s).s2);
    _sse_vector_i_mul();
    _sse_vector_sub();
    _sse_store_nt((*phi[ix]).s0);

    _sse_load((*s).s1);
    _sse_load_up((*s).s3);
    _sse_vector_i_mul();
    _sse_vector_add();
    _sse_store_nt((*phi[ix]).s1);
    ix++;
    s++;
  }

#    if (defined MPI && !defined _NO_COMM)
  xchange_halffield(); 
#    endif
  s = l;
  phi = NBPointer[2 + ieo];
  if(ieo == 0) {
    U = g_gauge_field_copy[1][0];
  }
  else {
    U = g_gauge_field_copy[0][0];
  }
  _prefetch_su3(U);
  
  /* Now we sum up and expand to a full spinor */
  ix = 0;
  for(i = 0; i < (VOLUME)/2; i++){
    /*********************** direction +0 ************************/
    _vector_assign(rs.s0, (*phi[ix]).s0);
    _vector_assign(rs.s2, (*phi[ix]).s0);
    _vector_assign(rs.s1, (*phi[ix]).s1);
    _vector_assign(rs.s3, (*phi[ix]).s1);
    ix++;

    /*********************** direction -0 ************************/
    _prefetch_su3(U+predist);
      
    _sse_load((*phi[ix]).s0);
    _sse_su3_inverse_multiply((*U));
    _sse_vector_cmplxcg_mul(ka0);
    _sse_load(rs.s0);
    _sse_vector_add();
    _sse_store(rs.s0);

    _sse_load(rs.s2);
    _sse_vector_sub();
    _sse_store(rs.s2);

    _sse_load((*phi[ix]).s1);
    _sse_su3_inverse_multiply((*U));
    _sse_vector_cmplxcg_mul(ka0);

    _sse_load(rs.s1);
    _sse_vector_add();
    _sse_store(rs.s1);

    _sse_load(rs.s3);
    _sse_vector_sub();
    _sse_store(rs.s3);

    ix++;
    U++;
    /*********************** direction +1 ************************/
    _sse_load_up((*phi[ix]).s0);
    _sse_load(rs.s0);
    _sse_vector_add();
    _sse_store(rs.s0);

    _sse_load(rs.s3);
    _sse_vector_i_mul();      
    _sse_vector_sub();
    _sse_store(rs.s3); 

    _sse_load_up((*phi[ix]).s1);
    _sse_load(rs.s1);
    _sse_vector_add();
    _sse_store(rs.s1);

    _sse_load(rs.s2);
    _sse_vector_i_mul();      
    _sse_vector_sub();
    _sse_store(rs.s2);       
    ix++;

    /*********************** direction -1 ************************/

    _prefetch_su3(U+predist);

    _sse_load((*phi[ix]).s0);
    _sse_su3_inverse_multiply((*U));
    _sse_vector_cmplxcg_mul(ka1);

    _sse_load(rs.s0);
    _sse_vector_add();
    _sse_store(rs.s0);

    _sse_load(rs.s3);
    _sse_vector_i_mul();      
    _sse_vector_add();
    _sse_store(rs.s3);

    _sse_load((*phi[ix]).s1);
      
    _sse_su3_inverse_multiply((*U));
    _sse_vector_cmplxcg_mul(ka1);

    _sse_load(rs.s1);
    _sse_vector_add();
    _sse_store(rs.s1);

    _sse_load(rs.s2);
    _sse_vector_i_mul();      
    _sse_vector_add();
    _sse_store(rs.s2);
    ix++;
    U++;

    /*********************** direction +2 ************************/
    _sse_load_up((*phi[ix]).s0);
    _sse_load(rs.s0);
    _sse_vector_add();
    _sse_store(rs.s0);

    _sse_load(rs.s3);
    _sse_vector_add();
    _sse_store(rs.s3);

    _sse_load_up((*phi[ix]).s1);
    _sse_load(rs.s1);
    _sse_vector_add();
    _sse_store(rs.s1);

    _sse_load(rs.s2);
    _sse_vector_sub();
    _sse_store(rs.s2);      
    ix++;

    /*********************** direction -2 ************************/

    _prefetch_su3(U+predist);

    _sse_load((*phi[ix]).s0);
    _sse_su3_inverse_multiply((*U));
    _sse_vector_cmplxcg_mul(ka2);

    _sse_load(rs.s0);
    _sse_vector_add();
    _sse_store(rs.s0);

    _sse_load(rs.s3);
    _sse_vector_sub();
    _sse_store(rs.s3);

    _sse_load((*phi[ix]).s1);
      
    _sse_su3_inverse_multiply((*U));
    _sse_vector_cmplxcg_mul(ka2);

    _sse_load(rs.s1);
    _sse_vector_add();
    _sse_store(rs.s1);

    _sse_load(rs.s2);
    _sse_vector_add();
    _sse_store(rs.s2);      
    ix++;
    U++;
    /*********************** direction +3 ************************/
    _sse_load_up((*phi[ix]).s0);

    _sse_load(rs.s0);
    _sse_vector_add();
    _sse_store(rs.s0);

    _sse_load(rs.s2);
    _sse_vector_i_mul();      
    _sse_vector_sub();
    _sse_store(rs.s2);

    _sse_load_up((*phi[ix]).s1);

    _sse_load(rs.s1);
    _sse_vector_add();
    _sse_store(rs.s1);

    _sse_load(rs.s3);
    _sse_vector_i_mul();      
    _sse_vector_add();
    _sse_store(rs.s3);

    ix++;
    /*********************** direction -3 ************************/

    _prefetch_su3(U+predist); 
    _prefetch_spinor(s+1);

    _sse_load((*phi[ix]).s0);
      
    _sse_su3_inverse_multiply((*U));
    _sse_vector_cmplxcg_mul(ka3);

    _sse_load(rs.s0);
    _sse_vector_add();
    _sse_store_nt((*s).s0);

    _sse_load(rs.s2);
    _sse_vector_i_mul();      
    _sse_vector_add();
    _sse_store_nt((*s).s2);

    _sse_load((*phi[ix]).s1);
      
    _sse_su3_inverse_multiply((*U));
    _sse_vector_cmplxcg_mul(ka3);

    _sse_load(rs.s1);
    _sse_vector_add();
    _sse_store_nt((*s).s1);

    _sse_load(rs.s3);
    _sse_vector_i_mul();      
    _sse_vector_sub();
    _sse_store_nt((*s).s3);
    ix++;
    U++;
    s++;
  }
#ifdef _KOJAK_INST
#pragma pomp inst end(hoppingmatrix)
#endif
}

#  elif (defined BGL && defined XLC)

/**********************************
 *
 * Blue Gene/L Version
 *
 * Author: Carsten Urbach
 *
 **********************************/


/* 2. */
void Hopping_Matrix(const int ieo, spinor * const l, spinor * const k){
  int i, ix;
  su3 * restrict U ALIGN;
  spinor * restrict s ALIGN;
  halfspinor * restrict * phi ALIGN;
  halfspinor32 * restrict * phi32 ALIGN;
  /* We have 32 registers available */
  double _Complex reg00, reg01, reg02, reg03, reg04, reg05;
  double _Complex reg10, reg11, reg12, reg13, reg14, reg15;
  /* For the gauge field, reuse the first three!*/
  double _Complex u00, u01, u02, u10, u11, u12;
  double _Complex reg20, reg21;
  /* The following contains the result spinor (12 regs) */
  double _Complex rs00, rs01, rs02, rs10, rs11, rs12, rs20, rs21, rs22, 
    rs30, rs31, rs32;
#ifdef _KOJAK_INST
#pragma pomp inst begin(hoppingmatrix)
#endif
#pragma disjoint(*s, *U)

#ifdef _GAUGE_COPY
  if(g_update_gauge_copy) {
    update_backward_gauge();
  }
#endif

  __alignx(16, l);
  __alignx(16, k);
  if(g_sloppy_precision == 1 && g_sloppy_precision_flag == 1) {
    /* Take the 64 Bit precision part and replace */
    /* _bgl_load_reg0|1 with _bgl_load_reg0|1_32 */
    /* _bgl_load_rs0|1|2|3 with _bgl_load_rs0|1|2|3_32*/
    /* phi with phi32*/
    /* _bgl_store_reg0|1 with _bgl_store_reg0|1_32 */
    /* _bgl_store_reg0|1_up with _bgl_store_reg0|1_up_32 */
    /* HalfSpinor with Halfspinor32 */
    /* _bgl_load_rs0|1 with _bgl_load_rs0|1_32*/
    /* xchange_halffield with xchange_halffield_32 */
    __alignx(16, HalfSpinor32);
    /* We will run through the source vector now */
    /* instead of the solution vector            */
    s = k;
    _prefetch_spinor(s); 

    /* s contains the source vector */

    if(ieo == 0) {
      U = g_gauge_field_copy[0][0];
    }
    else {
      U = g_gauge_field_copy[1][0];
    }
    phi32 = NBPointer32[ieo];

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

      _bgl_store_reg0_up_32((*phi32[ix]).s0);
      _bgl_store_reg1_up_32((*phi32[ix]).s1);
      U++;
      ix++;

      /*********************** direction -0 ************************/
      _bgl_vector_sub_rs2_from_rs0_reg0();
      _bgl_vector_sub_rs3_from_rs1_reg1();

      _bgl_store_reg0_32((*phi32[ix]).s0);
      _bgl_store_reg1_32((*phi32[ix]).s1);
      ix++;

      /*********************** direction +1 ************************/
      _prefetch_su3(U+1);
      _bgl_vector_i_mul_add_rs3_to_rs0_reg0();
      _bgl_vector_i_mul_add_rs2_to_rs1_reg1();

      _bgl_su3_multiply_double((*U));
      _bgl_vector_cmplx_mul_double(ka1);

      _bgl_store_reg0_up_32((*phi32[ix]).s0);
      _bgl_store_reg1_up_32((*phi32[ix]).s1);
      ix++;
      U++;

      /*********************** direction -1 ************************/
      _bgl_vector_i_mul_sub_rs3_from_rs0_reg0();
      _bgl_vector_i_mul_sub_rs2_from_rs1_reg1();

      _bgl_store_reg0_32((*phi32[ix]).s0);
      _bgl_store_reg1_32((*phi32[ix]).s1);
      ix++;


      /*********************** direction +2 ************************/
      _prefetch_su3(U+1);

      _bgl_vector_add_rs3_to_rs0_reg0();
      _bgl_vector_sub_rs2_from_rs1_reg1();

      _bgl_su3_multiply_double((*U));
      _bgl_vector_cmplx_mul_double(ka2);

      _bgl_store_reg0_up_32((*phi32[ix]).s0);
      _bgl_store_reg1_up_32((*phi32[ix]).s1);
      ix++;
      U++;

      /*********************** direction -2 ************************/
      _bgl_vector_sub_rs3_from_rs0_reg0();
      _bgl_vector_add_rs2_to_rs1_reg1();

      _bgl_store_reg0_32((*phi32[ix]).s0);
      _bgl_store_reg1_32((*phi32[ix]).s1);
      ix++;

      /*********************** direction +3 ************************/
      _prefetch_su3(U+1); 

      _bgl_vector_i_mul_add_rs2_to_rs0_reg0();
      _bgl_vector_i_mul_sub_rs3_from_rs1_reg1();

      _bgl_su3_multiply_double((*U));
      _bgl_vector_cmplx_mul_double(ka3);

      _bgl_store_reg0_up_32((*phi32[ix]).s0);
      _bgl_store_reg1_up_32((*phi32[ix]).s1);
      ix++;
      U++;

      /*********************** direction -3 ************************/
      _bgl_vector_i_mul_sub_rs2_from_rs0_reg0();
      _bgl_vector_i_mul_add_rs3_to_rs1_reg1();

      _bgl_store_reg0_32((*phi32[ix]).s0);
      _bgl_store_reg1_32((*phi32[ix]).s1);
      ix++;

      /************************ end of loop ************************/
    }

#    if (defined MPI && !defined _NO_COMM)
    xchange_halffield32(); 
#    endif
    s = l;
    phi32 = NBPointer32[2 + ieo];
    if(ieo == 0) {
      U = g_gauge_field_copy[1][0];
    }
    else {
      U = g_gauge_field_copy[0][0];
    }
    _prefetch_halfspinor(phi32[0]);
    _prefetch_su3(U);
  
    /* Now we sum up and expand to a full spinor */
    ix = 0;
    /*   _prefetch_spinor_for_store(s); */
    for(i = 0; i < (VOLUME)/2; i++){
      /* This causes a lot of trouble, do we understand this? */
      /*     _prefetch_spinor_for_store(s); */
      _prefetch_halfspinor(phi32[ix+1]);
      /*********************** direction +0 ************************/
      _bgl_load_rs0_32((*phi32[ix]).s0);
      rs20 = rs00;
      rs21 = rs01;
      rs22 = rs02;
      _bgl_load_rs1_32((*phi32[ix]).s1);
      rs30 = rs10;
      rs31 = rs11;
      rs32 = rs12;
      ix++;
      /*********************** direction -0 ************************/
      _prefetch_su3(U+1);

      _bgl_load_reg0_32((*phi32[ix]).s0);
      _bgl_load_reg1_32((*phi32[ix]).s1);

      _bgl_su3_inverse_multiply_double((*U));
      _bgl_vector_cmplxcg_mul_double(ka0);
      /* result is in the upper parts of reg0? and reg1? */

      _bgl_add_to_rs0_reg0();
      _bgl_sub_from_rs2_reg0();
      _bgl_add_to_rs1_reg1();
      _bgl_sub_from_rs3_reg1();
      U++;
      ix++;
      /*********************** direction +1 ************************/
      _bgl_load_reg0_up_32((*phi32[ix]).s0);
      _bgl_load_reg1_up_32((*phi32[ix]).s1);

      _bgl_add_to_rs0_reg0();
      _bgl_i_mul_sub_from_rs3_reg0();
      _bgl_add_to_rs1_reg1();
      _bgl_i_mul_sub_from_rs2_reg1();
      ix++;
      /*********************** direction -1 ************************/
      _prefetch_su3(U+1);

      _bgl_load_reg0_32((*phi32[ix]).s0);
      _bgl_load_reg1_32((*phi32[ix]).s1);

      _bgl_su3_inverse_multiply_double((*U));
      _bgl_vector_cmplxcg_mul_double(ka1);

      _bgl_add_to_rs0_reg0();
      _bgl_add_to_rs1_reg1();
      _bgl_i_mul_add_to_rs3_reg0();
      _bgl_i_mul_add_to_rs2_reg1();      
      U++;
      ix++;
      /*********************** direction +2 ************************/
      _bgl_load_reg0_up_32((*phi32[ix]).s0);
      _bgl_load_reg1_up_32((*phi32[ix]).s1);

      _bgl_add_to_rs0_reg0();
      _bgl_add_to_rs1_reg1();
      _bgl_sub_from_rs2_reg1();
      _bgl_add_to_rs3_reg0();
      ix++;
      /*********************** direction -2 ************************/
      _prefetch_su3(U+1);

      _bgl_load_reg0_32((*phi32[ix]).s0);
      _bgl_load_reg1_32((*phi32[ix]).s1);

      _bgl_su3_inverse_multiply_double((*U));
      _bgl_vector_cmplxcg_mul_double(ka2);

      _bgl_add_to_rs0_reg0();
      _bgl_add_to_rs1_reg1();
      _bgl_add_to_rs2_reg1();
      _bgl_sub_from_rs3_reg0();
      U++;
      ix++;
      /*********************** direction +3 ************************/
      _bgl_load_reg0_up_32((*phi32[ix]).s0);
      _bgl_load_reg1_up_32((*phi32[ix]).s1);

      _bgl_add_to_rs0_reg0();
      _bgl_add_to_rs1_reg1();
      _bgl_i_mul_sub_from_rs2_reg0();
      _bgl_i_mul_add_to_rs3_reg1();
      ix++;
      /*********************** direction -3 ************************/
      _prefetch_su3(U+1);

      _bgl_load_reg0_32((*phi32[ix]).s0);
      _bgl_load_reg1_32((*phi32[ix]).s1);

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

      U++;
      ix++;
      s++;
    }
  }
  else {
    __alignx(16, HalfSpinor);
    /* We will run through the source vector now */
    /* instead of the solution vector            */
    s = k;
    _prefetch_spinor(s); 

    /* s contains the source vector */

    if(ieo == 0) {
      U = g_gauge_field_copy[0][0];
    }
    else {
      U = g_gauge_field_copy[1][0];
    }
    phi = NBPointer[ieo];

    _prefetch_su3(U);
    /**************** loop over all lattice sites ******************/
    ix=0;
    for(i = 0; i < (VOLUME)/2; i++){
      _prefetch_halfspinor(phi[ix+4]);
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
      _prefetch_halfspinor(phi[ix+4]);
      _bgl_vector_sub_rs2_from_rs0_reg0();
      _bgl_vector_sub_rs3_from_rs1_reg1();

      _bgl_store_reg0((*phi[ix]).s0);
      _bgl_store_reg1((*phi[ix]).s1);
      ix++;

      /*********************** direction +1 ************************/
      _prefetch_halfspinor(phi[ix+4]);
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
      _prefetch_halfspinor(phi[ix+4]);
      _bgl_vector_i_mul_sub_rs3_from_rs0_reg0();
      _bgl_vector_i_mul_sub_rs2_from_rs1_reg1();

      _bgl_store_reg0((*phi[ix]).s0);
      _bgl_store_reg1((*phi[ix]).s1);
      ix++;


      /*********************** direction +2 ************************/
      _prefetch_halfspinor(phi[ix+4]);
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
      _prefetch_halfspinor(phi[ix+4]);
      _bgl_vector_sub_rs3_from_rs0_reg0();
      _bgl_vector_add_rs2_to_rs1_reg1();

      _bgl_store_reg0((*phi[ix]).s0);
      _bgl_store_reg1((*phi[ix]).s1);
      ix++;

      /*********************** direction +3 ************************/
      _prefetch_halfspinor(phi[ix+4]);
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
      _prefetch_halfspinor(phi[ix+4]);
      _bgl_vector_i_mul_sub_rs2_from_rs0_reg0();
      _bgl_vector_i_mul_add_rs3_to_rs1_reg1();

      _bgl_store_reg0((*phi[ix]).s0);
      _bgl_store_reg1((*phi[ix]).s1);
      ix++;

      /************************ end of loop ************************/
    }

#    if (defined MPI && !defined _NO_COMM)
    xchange_halffield(); 
#    endif
    s = l;
    phi = NBPointer[2 + ieo];
    _prefetch_halfspinor(phi[0]);
    if(ieo == 0) {
      U = g_gauge_field_copy[1][0];
    }
    else {
      U = g_gauge_field_copy[0][0];
    }
    _prefetch_su3(U);
  
    /* Now we sum up and expand to a full spinor */
    ix = 0;
    /*   _prefetch_spinor_for_store(s); */
    for(i = 0; i < (VOLUME)/2; i++){
      /* This causes a lot of trouble, do we understand this? */
      _prefetch_halfspinor(phi[ix+3]);
      /*********************** direction +0 ************************/
      _bgl_load_rs0((*phi[ix]).s0);
      rs20 = rs00;
      rs21 = rs01;
      rs22 = rs02;
      _bgl_load_rs1((*phi[ix]).s1);
      rs30 = rs10;
      rs31 = rs11;
      rs32 = rs12;
      ix++;
      /*********************** direction -0 ************************/
      _prefetch_halfspinor(phi[ix+3]);
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
      U++;
      ix++;
      /*********************** direction +1 ************************/
      _prefetch_halfspinor(phi[ix+3]);
      _bgl_load_reg0_up((*phi[ix]).s0);
      _bgl_load_reg1_up((*phi[ix]).s1);

      _bgl_add_to_rs0_reg0();
      _bgl_i_mul_sub_from_rs3_reg0();
      _bgl_add_to_rs1_reg1();
      _bgl_i_mul_sub_from_rs2_reg1();
      ix++;
      /*********************** direction -1 ************************/
      _prefetch_halfspinor(phi[ix+3]);
      _prefetch_su3(U+1);

      _bgl_load_reg0((*phi[ix]).s0);
      _bgl_load_reg1((*phi[ix]).s1);

      _bgl_su3_inverse_multiply_double((*U));
      _bgl_vector_cmplxcg_mul_double(ka1);

      _bgl_add_to_rs0_reg0();
      _bgl_add_to_rs1_reg1();
      _bgl_i_mul_add_to_rs3_reg0();
      _bgl_i_mul_add_to_rs2_reg1();      
      U++;
      ix++;
      /*********************** direction +2 ************************/
      _prefetch_halfspinor(phi[ix+3]);
      _bgl_load_reg0_up((*phi[ix]).s0);
      _bgl_load_reg1_up((*phi[ix]).s1);

      _bgl_add_to_rs0_reg0();
      _bgl_add_to_rs1_reg1();
      _bgl_sub_from_rs2_reg1();
      _bgl_add_to_rs3_reg0();
      ix++;
      /*********************** direction -2 ************************/

      _prefetch_halfspinor(phi[ix+3]);
      _prefetch_su3(U+1);

      _bgl_load_reg0((*phi[ix]).s0);
      _bgl_load_reg1((*phi[ix]).s1);

      _bgl_su3_inverse_multiply_double((*U));
      _bgl_vector_cmplxcg_mul_double(ka2);

      _bgl_add_to_rs0_reg0();
      _bgl_add_to_rs1_reg1();
      _bgl_add_to_rs2_reg1();
      _bgl_sub_from_rs3_reg0();
      U++;
      ix++;
      /*********************** direction +3 ************************/

      _prefetch_halfspinor(phi[ix+3]);
      _bgl_load_reg0_up((*phi[ix]).s0);
      _bgl_load_reg1_up((*phi[ix]).s1);

      _bgl_add_to_rs0_reg0();
      _bgl_add_to_rs1_reg1();
      _bgl_i_mul_sub_from_rs2_reg0();
      _bgl_i_mul_add_to_rs3_reg1();
      ix++;
      /*********************** direction -3 ************************/
      _prefetch_spinor(s);
      _prefetch_halfspinor(phi[ix+3]);
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

      U++;
      ix++;
      s++;
    }
  }
#ifdef _KOJAK_INST
#pragma pomp inst end(hoppingmatrix)
#endif
}

#  else

/* 3. */
/* l output , k input*/
/* for ieo=0, k resides on  odd sites and l on even sites */
void Hopping_Matrix(const int ieo, spinor * const l, spinor * const k){
  int i,ix;
  su3 * restrict U ALIGN;
  spinor * restrict s ALIGN;
  spinor rs;
  static su3_vector psi, chi, psi2, chi2;
  halfspinor * restrict * phi ALIGN;
  halfspinor32 * restrict * phi32 ALIGN;
#ifdef _KOJAK_INST
#pragma pomp inst begin(hoppingmatrix)
#endif
#ifdef XLC
#pragma disjoint(*l, *k, *U, *s)
#endif

#ifdef _GAUGE_COPY
  if(g_update_gauge_copy) {
    update_backward_gauge();
  }
#endif

  if(k == l){
    printf("Error in H_psi (simple.c):\n");
    printf("Arguments k and l must be different\n");
    printf("Program aborted\n");
    exit(1);
  }
  s = k;

  if(ieo == 0) {
    U = g_gauge_field_copy[0][0];
  }
  else {
    U = g_gauge_field_copy[1][0];
  }
  if(g_sloppy_precision == 1 && g_sloppy_precision_flag == 1) {
    phi32 = NBPointer32[ieo];
      
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
      _complex_times_vector32((*phi32[ix]).s0, ka0, chi);
      
      _vector_add(psi, rs.s1, rs.s3);

      _su3_multiply(chi,(*U),psi);
      _complex_times_vector32((*phi32[ix]).s1, ka0, chi);
            
      U++;
      ix++;
    
      /*********************** direction -0 ************************/

      _vector_sub32((*phi32[ix]).s0, rs.s0, rs.s2);
      _vector_sub32((*phi32[ix]).s1, rs.s1, rs.s3);

      ix++;

      /*********************** direction +1 ************************/

      _vector_i_add(psi, rs.s0, rs.s3);

      _su3_multiply(chi, (*U), psi);
      _complex_times_vector32((*phi32[ix]).s0, ka1, chi);

      _vector_i_add(psi, rs.s1, rs.s2);

      _su3_multiply(chi, (*U), psi);
      _complex_times_vector32((*phi32[ix]).s1, ka1, chi);

      U++;
      ix++;

      /*********************** direction -1 ************************/

      _vector_i_sub32((*phi32[ix]).s0, rs.s0, rs.s3);
      _vector_i_sub32((*phi32[ix]).s1, rs.s1, rs.s2);

      ix++;
      /*********************** direction +2 ************************/

      _vector_add(psi, rs.s0, rs.s3);

      _su3_multiply(chi,(*U),psi);
      _complex_times_vector32((*phi32[ix]).s0, ka2, chi);

      _vector_sub(psi, rs.s1, rs.s2);

      _su3_multiply(chi,(*U),psi);
      _complex_times_vector32((*phi32[ix]).s1, ka2, chi);
      
      U++;
      ix++;

      /*********************** direction -2 ************************/

      _vector_sub32((*phi32[ix]).s0, rs.s0, rs.s3);
      _vector_add32((*phi32[ix]).s1, rs.s1, rs.s2);
      ix++;

      /*********************** direction +3 ************************/

      _vector_i_add(psi, rs.s0, rs.s2);
      
      _su3_multiply(chi, (*U), psi);
      _complex_times_vector32((*phi32[ix]).s0, ka3, chi);


      _vector_i_sub(psi, rs.s1, rs.s3);

      _su3_multiply(chi,(*U),psi);
      _complex_times_vector32((*phi32[ix]).s1, ka3, chi);

      U++;
      ix++;
      /*********************** direction -3 ************************/

      _vector_i_sub32((*phi32[ix]).s0, rs.s0, rs.s2);
      _vector_i_add32((*phi32[ix]).s1, rs.s1, rs.s3);

      ix++;
      /************************ end of loop ************************/
    }
#    if (defined MPI && !defined _NO_COMM)
    xchange_halffield32(); 
#    endif
    s = l;
    phi32 = NBPointer32[2 + ieo];
    if(ieo == 0) {
      U = g_gauge_field_copy[1][0];
    }
    else {
      U = g_gauge_field_copy[0][0];
    }

    ix = 0;
    for(i = 0; i < (VOLUME)/2; i++){
      /*********************** direction +0 ************************/
      _vector_assign32(rs.s0, (*phi32[ix]).s0);
      _vector_assign32(rs.s2, (*phi32[ix]).s0);
      _vector_assign32(rs.s1, (*phi32[ix]).s1);
      _vector_assign32(rs.s3, (*phi32[ix]).s1);
      ix++;
      /*********************** direction -0 ************************/
      _vector_assign32(psi, (*phi32[ix]).s0);
      _su3_inverse_multiply(chi,(*U), psi);
      _complexcjg_times_vector(psi,ka0,chi);

      _vector_add_assign(rs.s0, psi);
      _vector_sub_assign(rs.s2, psi);

      _vector_assign32(psi, (*phi32[ix]).s1);
      _su3_inverse_multiply(chi,(*U), psi);
      _complexcjg_times_vector(psi,ka0,chi);
      
      _vector_add_assign(rs.s1, psi);
      _vector_sub_assign(rs.s3, psi);
      ix++;
      U++;
      /*********************** direction +1 ************************/

      _vector_add_assign32(rs.s0, (*phi32[ix]).s0);
      _vector_i_sub_assign32(rs.s3, (*phi32[ix]).s0);

      _vector_add_assign32(rs.s1, (*phi32[ix]).s1);
      _vector_i_sub_assign32(rs.s2, (*phi32[ix]).s1);
    
      ix++;
      /*********************** direction -1 ************************/
      _vector_assign32(psi, (*phi32[ix]).s0);
      _su3_inverse_multiply(chi,(*U), psi);
      _complexcjg_times_vector(psi,ka1,chi);

      _vector_add_assign(rs.s0, psi);
      _vector_i_add_assign(rs.s3, psi);

      _vector_assign32(psi, (*phi32[ix]).s1);
      _su3_inverse_multiply(chi,(*U), psi);
      _complexcjg_times_vector(psi,ka1,chi);

      _vector_add_assign(rs.s1, psi);
      _vector_i_add_assign(rs.s2, psi);

      U++;
      ix++;

      /*********************** direction +2 ************************/

      _vector_add_assign32(rs.s0, (*phi32[ix]).s0);
      _vector_add_assign32(rs.s3, (*phi32[ix]).s0);

      _vector_add_assign32(rs.s1, (*phi32[ix]).s1);
      _vector_sub_assign32(rs.s2, (*phi32[ix]).s1);
    
      ix++;
      /*********************** direction -2 ************************/

      _vector_assign32(psi, (*phi32[ix]).s0);
      _su3_inverse_multiply(chi,(*U), psi);
      _complexcjg_times_vector(psi,ka2,chi);

      _vector_add_assign(rs.s0, psi);
      _vector_sub_assign(rs.s3, psi);

      _vector_assign32(psi, (*phi32[ix]).s1);
      _su3_inverse_multiply(chi, (*U), psi);
      _complexcjg_times_vector(psi,ka2,chi);
      
      _vector_add_assign(rs.s1, psi);
      _vector_add_assign(rs.s2, psi);

      U++;
      ix++;
      /*********************** direction +3 ************************/

      _vector_add_assign32(rs.s0, (*phi32[ix]).s0);
      _vector_i_sub_assign32(rs.s2, (*phi32[ix]).s0);

      _vector_add_assign32(rs.s1, (*phi32[ix]).s1);
      _vector_i_add_assign32(rs.s3, (*phi32[ix]).s1);

      ix++;

      /*********************** direction -3 ************************/

      _vector_assign32(psi, (*phi32[ix]).s0);
      _su3_inverse_multiply(chi,(*U), psi);
      _complexcjg_times_vector(psi,ka3,chi);
      
      _vector_add((*s).s0, rs.s0, psi);
      _vector_i_add((*s).s2, rs.s2, psi);

      _vector_assign32(psi, (*phi32[ix]).s1);
      _su3_inverse_multiply(chi,(*U), psi);
      _complexcjg_times_vector(psi,ka3,chi);

      _vector_add((*s).s1, rs.s1, psi);
      _vector_i_sub((*s).s3, rs.s3, psi);

      U++;
      ix++;
      s++;
    }
  }
  else {
    phi = NBPointer[ieo];
      
    /**************** loop over all lattice sites ****************/
    ix=0;
    /* #pragma ivdep*/
    for(i = 0; i < (VOLUME)/2; i++){
      _vector_assign(rs.s0, (*s).s0);
      _vector_assign(rs.s1, (*s).s1);
      _vector_assign(rs.s2, (*s).s2);
      _vector_assign(rs.s3, (*s).s3);
      s++;
      /*********************** direction +0 ************************/
      
      _vector_add(psi, rs.s0, rs.s2);
      _vector_add(psi2, rs.s1, rs.s3);
      _su3_multiply(chi,(*U),psi);
      _su3_multiply(chi2,(*U),psi2);
      _complex_times_vector((*phi[ix]).s0, ka0, chi);
      _complex_times_vector((*phi[ix]).s1, ka0, chi2);
            
      U++;
      ix++;
    
      /*********************** direction -0 ************************/

      _vector_sub((*phi[ix]).s0, rs.s0, rs.s2);
      _vector_sub((*phi[ix]).s1, rs.s1, rs.s3);

      ix++;

      /*********************** direction +1 ************************/

      _vector_i_add(psi, rs.s0, rs.s3);
      _vector_i_add(psi2, rs.s1, rs.s2);
      _su3_multiply(chi, (*U), psi);
      _su3_multiply(chi2, (*U), psi2);
      _complex_times_vector((*phi[ix]).s0, ka1, chi);
      _complex_times_vector((*phi[ix]).s1, ka1, chi2);

      U++;
      ix++;

      /*********************** direction -1 ************************/

      _vector_i_sub((*phi[ix]).s0, rs.s0, rs.s3);
      _vector_i_sub((*phi[ix]).s1, rs.s1, rs.s2);

      ix++;
      /*********************** direction +2 ************************/

      _vector_add(psi, rs.s0, rs.s3);
      _vector_sub(psi2, rs.s1, rs.s2);
      _su3_multiply(chi,(*U),psi);
      _su3_multiply(chi2,(*U),psi2);
      _complex_times_vector((*phi[ix]).s0, ka2, chi);
      _complex_times_vector((*phi[ix]).s1, ka2, chi2);
      
      U++;
      ix++;

      /*********************** direction -2 ************************/

      _vector_sub((*phi[ix]).s0, rs.s0, rs.s3);
      _vector_add((*phi[ix]).s1, rs.s1, rs.s2);
      ix++;

      /*********************** direction +3 ************************/

      _vector_i_add(psi, rs.s0, rs.s2);
      _vector_i_sub(psi2, rs.s1, rs.s3);      
      _su3_multiply(chi, (*U), psi);
      _su3_multiply(chi2,(*U),psi2);
      _complex_times_vector((*phi[ix]).s0, ka3, chi);
      _complex_times_vector((*phi[ix]).s1, ka3, chi2);

      U++;
      ix++;
      /*********************** direction -3 ************************/

      _vector_i_sub((*phi[ix]).s0, rs.s0, rs.s2);
      _vector_i_add((*phi[ix]).s1, rs.s1, rs.s3);

      ix++;
      /************************ end of loop ************************/
    }
#    if (defined MPI && !defined _NO_COMM)
    xchange_halffield(); 
#    endif
    s = l;
    phi = NBPointer[2 + ieo];
    if(ieo == 0) {
      U = g_gauge_field_copy[1][0];
    }
    else {
      U = g_gauge_field_copy[0][0];
    }

    ix = 0;
    /* #pragma ivdep */
    for(i = 0; i < (VOLUME)/2; i++){
      /*********************** direction +0 ************************/
      _vector_assign(rs.s0, (*phi[ix]).s0);
      _vector_assign(rs.s2, (*phi[ix]).s0);
      _vector_assign(rs.s1, (*phi[ix]).s1);
      _vector_assign(rs.s3, (*phi[ix]).s1);
      ix++;
      /*********************** direction -0 ************************/
      _su3_inverse_multiply(chi,(*U),(*phi[ix]).s0);
      _su3_inverse_multiply(chi2,(*U),(*phi[ix]).s1);
      _complexcjg_times_vector(psi,ka0,chi);
      _complexcjg_times_vector(psi2,ka0,chi2);
      _vector_add_assign(rs.s0, psi);
      _vector_sub_assign(rs.s2, psi);
      _vector_add_assign(rs.s1, psi2);
      _vector_sub_assign(rs.s3, psi2);
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
      _su3_inverse_multiply(chi2, (*U), (*phi[ix]).s1);
      _complexcjg_times_vector(psi,ka1,chi);
      _complexcjg_times_vector(psi2,ka1,chi2);
      _vector_add_assign(rs.s0, psi);
      _vector_i_add_assign(rs.s3, psi);
      _vector_add_assign(rs.s1, psi2);
      _vector_i_add_assign(rs.s2, psi2);

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
      _su3_inverse_multiply(chi2, (*U), (*phi[ix]).s1);
      _complexcjg_times_vector(psi,ka2,chi);
      _complexcjg_times_vector(psi2,ka2,chi2);
      _vector_add_assign(rs.s0, psi);
      _vector_sub_assign(rs.s3, psi);
      _vector_add_assign(rs.s1, psi2);
      _vector_add_assign(rs.s2, psi2);

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
      _su3_inverse_multiply(chi2, (*U), (*phi[ix]).s1);
      _complexcjg_times_vector(psi,ka3,chi);
      _complexcjg_times_vector(psi2,ka3,chi2);      
      _vector_add((*s).s0, rs.s0, psi);
      _vector_i_add((*s).s2, rs.s2, psi);
      _vector_add((*s).s1, rs.s1, psi2);
      _vector_i_sub((*s).s3, rs.s3, psi2);

      U++;
      ix++;
      s++;
    }
  }
#ifdef _KOJAK_INST
#pragma pomp inst end(hoppingmatrix)
#endif
}
#  endif

#else /* thats _USE_HALFSPINOR */

#  if ((defined SSE2)||(defined SSE3))

#    if (defined _USE_TSPLITPAR) /* needs also SSE */

/***********************************
 * 
 * Aurora version
 * Author: Luigi Scorzato (scorzato@ect.it)
 * (last modified 20.4.2009)
 * The strategy of the code is explained in the file Strategy.txt
 *
 ************************************/

/* 4. */
/* input on k; output on l */
void Hopping_Matrix(const int ieo, spinor * const l, spinor * const k){
  int icx,icz,ioff;
  int ix,iz;
  int x0,icx0,jj;
  su3 *restrict up;
  su3 * restrict um;
  spinor * restrict sp;
  spinor * restrict sm;
  spinor * restrict rn;

# if (defined MPI)
#  ifdef PARALLELX
#   define  REQC 4
#  elif defined PARALLELXY
#   define  REQC 8
#  elif defined PARALLELXYZ
#   define  REQC 12
#  endif
  MPI_Request requests[REQC];
  MPI_Status status[REQC];
# endif

#ifdef _GAUGE_COPY
  if(g_update_gauge_copy) {
    update_backward_gauge();
  }
#endif

  if(ieo == 0){ /* even out - odd in */
    ioff = 0;
  } 
  else{ /* odd out - even in */
    ioff = (VOLUME+RAND)/2;
  }

  /* Loop over time direction. This is the outmost loop */
  for(x0=0;x0<T;x0++){

    /* start the communication of the timslice borders (non-blocking send and receive)*/
#    if (defined MPI && !defined _NO_COMM)
   xchange_field_open(k, ieo, x0, requests, status);
#    endif
    

  /* loop over timeslice. At: contribution of timelike links  */
   icx0=g_1st_eot[x0][ieo];
   jj =0;
   um=&g_gauge_field_copyt[icx0][0]-1; /* allowed? */
   for(icx = icx0; icx < icx0+TEOSLICE; icx++){
     rn=l+(icx-ioff);
    /*********************** direction +0 ************************/

    sp=k+g_iup_eo[icx][0]; /* all sp,sm,up,um could be moved up */
    up=um+1;

    _sse_load((*sp).s0);
    _sse_load_up((*sp).s2);
    _sse_vector_add();

    _sse_su3_multiply((*up));
    _sse_vector_cmplx_mul(ka0);
    _sse_store_up((*rn).s0);
    _sse_store_up((*rn).s2);      
      
    _sse_load((*sp).s1);
    _sse_load_up((*sp).s3);
    _sse_vector_add();
      
    _sse_su3_multiply((*up));
    _sse_vector_cmplx_mul(ka0);
    _sse_store_up((*rn).s1);
    _sse_store_up((*rn).s3); 

    /*********************** direction -0 ************************/

    sm=k+g_idn_eo[icx][0];
    um=up+1;

    _sse_load((*sm).s0);
    _sse_load_up((*sm).s2);
    _sse_vector_sub();
      
    _sse_su3_inverse_multiply((*um));
    _sse_vector_cmplxcg_mul(ka0);
      
    _sse_load((*rn).s0);
    _sse_vector_add();
    _sse_store((*rn).s0);

    _sse_load((*rn).s2);
    _sse_vector_sub();
    _sse_store((*rn).s2);
      
    _sse_load((*sm).s1);
    _sse_load_up((*sm).s3);
    _sse_vector_sub();
      
    _sse_su3_inverse_multiply((*um));
    _sse_vector_cmplxcg_mul(ka0);
      
    _sse_load((*rn).s1);
    _sse_vector_add();
    _sse_store((*rn).s1);

    _sse_load((*rn).s3);
    _sse_vector_sub();
    _sse_store((*rn).s3);
    jj++;
   } /* end of loop over timeslice (At)*/

       
   /* complete the communication of the timslice borders (and wait) */
#if (defined MPI && !defined _NO_COMM)
   xchange_field_close(requests, status, REQC); /*    MPI_Waitall */
#endif

   /* loop over timeslice. Bt: contribution of spacelike links  */
   um=&g_gauge_field_copys[icx0][0]-1;
   for(icx = icx0; icx < icx0+TEOSLICE; icx++){
    ix=g_eo2lexic[icx];
    rn=l+(icx-ioff);
    /*********************** direction +1 ************************/

    sp=k+g_iup_eo[icx][1];
    up=um+1;

    _sse_load((*sp).s0);
    _sse_load_up((*sp).s3);
    _sse_vector_i_mul();
    _sse_vector_add();

    _sse_su3_multiply((*up));
    _sse_vector_cmplx_mul(ka1);

    _sse_load((*rn).s0);
    _sse_vector_add();
    _sse_store((*rn).s0);

    _sse_load((*rn).s3);
    _sse_vector_i_mul();      
    _sse_vector_sub();
    _sse_store((*rn).s3); 
      
    _sse_load((*sp).s1);
    _sse_load_up((*sp).s2);
    _sse_vector_i_mul();
    _sse_vector_add();

    _sse_su3_multiply((*up));
    _sse_vector_cmplx_mul(ka1);

    _sse_load((*rn).s1);
    _sse_vector_add();
    _sse_store((*rn).s1);

    _sse_load((*rn).s2);
    _sse_vector_i_mul();      
    _sse_vector_sub();
    _sse_store((*rn).s2);       

    /*********************** direction -1 ************************/

    sm=k+g_idn_eo[icx][1];
    um=up+1;

    _sse_load((*sm).s0);
    _sse_load_up((*sm).s3);
    _sse_vector_i_mul();
    _sse_vector_sub();
      
    _sse_su3_inverse_multiply((*um));
    _sse_vector_cmplxcg_mul(ka1);
      
    _sse_load((*rn).s0);
    _sse_vector_add();
    _sse_store((*rn).s0);

    _sse_load((*rn).s3);
    _sse_vector_i_mul();      
    _sse_vector_add();
    _sse_store((*rn).s3);

    _sse_load((*sm).s1);
    _sse_load_up((*sm).s2);
    _sse_vector_i_mul();
    _sse_vector_sub();
      
    _sse_su3_inverse_multiply((*um));
    _sse_vector_cmplxcg_mul(ka1);
      
    _sse_load((*rn).s1);
    _sse_vector_add();
    _sse_store((*rn).s1);

    _sse_load((*rn).s2);
    _sse_vector_i_mul();      
    _sse_vector_add();
    _sse_store((*rn).s2);

    /*********************** direction +2 ************************/

    sp=k+g_iup_eo[icx][2];
    up=um+1;

    _sse_load((*sp).s0);
    _sse_load_up((*sp).s3);
    _sse_vector_add();

    _sse_su3_multiply((*up));
    _sse_vector_cmplx_mul(ka2);

    _sse_load((*rn).s0);
    _sse_vector_add();
    _sse_store((*rn).s0);

    _sse_load((*rn).s3);
    _sse_vector_add();
    _sse_store((*rn).s3);
      
    _sse_load((*sp).s1);
    _sse_load_up((*sp).s2);
    _sse_vector_sub();

    _sse_su3_multiply((*up));
    _sse_vector_cmplx_mul(ka2);

    _sse_load((*rn).s1);
    _sse_vector_add();
    _sse_store((*rn).s1);

    _sse_load((*rn).s2);
    _sse_vector_sub();
    _sse_store((*rn).s2);      

    /*********************** direction -2 ************************/

    sm=k+g_idn_eo[icx][2];
    um=up+1;

    _sse_load((*sm).s0);
    _sse_load_up((*sm).s3);
    _sse_vector_sub();
      
    _sse_su3_inverse_multiply((*um));
    _sse_vector_cmplxcg_mul(ka2);
      
    _sse_load((*rn).s0);
    _sse_vector_add();
    _sse_store((*rn).s0);

    _sse_load((*rn).s3);
    _sse_vector_sub();
    _sse_store((*rn).s3);
      
    _sse_load((*sm).s1);
    _sse_load_up((*sm).s2);
    _sse_vector_add();
      
    _sse_su3_inverse_multiply((*um));
    _sse_vector_cmplxcg_mul(ka2);
      
    _sse_load((*rn).s1);
    _sse_vector_add();
    _sse_store((*rn).s1);

    _sse_load((*rn).s2);
    _sse_vector_add();
    _sse_store((*rn).s2);      
      
    /*********************** direction +3 ************************/

    sp=k+g_iup_eo[icx][3];
    up=um+1;

    _sse_load((*sp).s0);
    _sse_load_up((*sp).s2);
    _sse_vector_i_mul();
    _sse_vector_add();

    _sse_su3_multiply((*up));
    _sse_vector_cmplx_mul(ka3);

    _sse_load((*rn).s0);
    _sse_vector_add();
    _sse_store((*rn).s0);

    _sse_load((*rn).s2);
    _sse_vector_i_mul();      
    _sse_vector_sub();
    _sse_store((*rn).s2);
      
    _sse_load((*sp).s1);
    _sse_load_up((*sp).s3);
    _sse_vector_i_mul();
    _sse_vector_sub();

    _sse_su3_multiply((*up));
    _sse_vector_cmplx_mul(ka3);

    _sse_load((*rn).s1);
    _sse_vector_add();
    _sse_store((*rn).s1);

    _sse_load((*rn).s3);
    _sse_vector_i_mul();      
    _sse_vector_add();
    _sse_store((*rn).s3);
      
    /*********************** direction -3 ************************/

    sm=k+g_idn_eo[icx][3];
    um=up+1;

    _sse_load((*sm).s0);
    _sse_load_up((*sm).s2);
    _sse_vector_i_mul();
    _sse_vector_sub();
      
    _sse_su3_inverse_multiply((*um));
    _sse_vector_cmplxcg_mul(ka3);

    _sse_load((*rn).s0);
    _sse_vector_add();
    _sse_store((*rn).s0);

    _sse_load((*rn).s2);
    _sse_vector_i_mul();
    _sse_vector_add();
    _sse_store((*rn).s2);

    _sse_load((*sm).s1);
    _sse_load_up((*sm).s3);
    _sse_vector_i_mul();
    _sse_vector_add();
      
    _sse_su3_inverse_multiply((*um));
    _sse_vector_cmplxcg_mul(ka3);

    _sse_load((*rn).s1);
    _sse_vector_add();
    _sse_store((*rn).s1);

    _sse_load((*rn).s3);
    _sse_vector_i_mul();      
    _sse_vector_sub();
    _sse_store((*rn).s3);
   }  /* end of loop over timeslice (Bt)*/
  } /* x0=0; x0<T */
}

#    else /* _USE_TSPLITPAR */

/* 5. */
/* input on k; output on l */
void Hopping_Matrix(const int ieo, spinor * const l, spinor * const k){
  int icx,icy,icz,ioff,ioff2;
  int ix,iy,iz;
  su3 *restrict up;
  su3 * restrict um;
  spinor * restrict sp;
  spinor * restrict sm;
  spinor * restrict rn;
  static spinor rs;
  
  /* for parallelization */
#ifdef _GAUGE_COPY
  if(g_update_gauge_copy) {
    update_backward_gauge();
  }
#endif

#    if (defined MPI && !defined _NO_COMM)
  xchange_field(k, ieo);
#    endif

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
#    if ((defined _GAUGE_COPY))
  up=&g_gauge_field_copy[ioff][0];
#    else
  up=&g_gauge_field[ix][0];
#    endif
   
  /**************** loop over all lattice sites ******************/
  for(icx = ioff; icx < (VOLUME/2+ioff); icx++){
    ix=g_eo2lexic[icx];
    /*********************** direction +0 ************************/

    iy=g_idn[ix][0]; icy=g_lexic2eosub[iy];

    sm=k+icy;
    _prefetch_spinor(sm);

#    if ((defined _GAUGE_COPY))
    um=up+1;
#    else
    um=&g_gauge_field[iy][0]; 
#    endif
    _prefetch_su3(um);
      
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
      
    sp=k+icy;
    _prefetch_spinor(sp);

#    if ((defined _GAUGE_COPY))
    up = um + 1;
#    else
    up+=1;
#    endif
    _prefetch_su3(up);

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

    sm=k+icy;
    _prefetch_spinor(sm);

#    ifndef _GAUGE_COPY
    um=&g_gauge_field[iy][1]; 
#    else
    um=up+1;
#    endif
    _prefetch_su3(um);

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

    sp=k+icy;
    _prefetch_spinor(sp);

#    if ((defined _GAUGE_COPY))
    up = um + 1;
#    else
    up+=1;
#    endif
    _prefetch_su3(up);

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

    sm=k+icy;
    _prefetch_spinor(sm);

#    ifndef _GAUGE_COPY
    um=&g_gauge_field[iy][2]; 
#    else
    um=up+1;
#    endif
    _prefetch_su3(um);

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

    sp=k+icy;
    _prefetch_spinor(sp);

#    if ((defined _GAUGE_COPY))
    up = um + 1;
#    else
    up+=1;
#    endif
    _prefetch_su3(up);

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

    sm=k+icy;
    _prefetch_spinor(sm);

#    ifndef _GAUGE_COPY
    um=&g_gauge_field[iy][3]; 
#    else
    um=up+1;
#    endif
    _prefetch_su3(um);

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
      
    sp=k+icy;
    _prefetch_spinor(sp);

#    if ((defined _GAUGE_COPY))
    up=&g_gauge_field_copy[icz][0];
#    else
    up=&g_gauge_field[iz][0];
#    endif
    _prefetch_su3(up);

    _sse_load((*sm).s0);
    _sse_load_up((*sm).s2);
    _sse_vector_i_mul();
    _sse_vector_sub();
      
    _sse_su3_inverse_multiply((*um));
    _sse_vector_cmplxcg_mul(ka3);

    rn=l+(icx-ioff);
      
    _sse_load(rs.s0);
    _sse_vector_add();
    _sse_store_nt((*rn).s0);

    _sse_load(rs.s2);
    _sse_vector_i_mul();      
    _sse_vector_add();
    _sse_store_nt((*rn).s2);

    _sse_load((*sm).s1);
    _sse_load_up((*sm).s3);
    _sse_vector_i_mul();
    _sse_vector_add();
      
    _sse_su3_inverse_multiply((*um));
    _sse_vector_cmplxcg_mul(ka3);

    _sse_load(rs.s1);
    _sse_vector_add();
    _sse_store_nt((*rn).s1);

    _sse_load(rs.s3);
    _sse_vector_i_mul();      
    _sse_vector_sub();
    _sse_store_nt((*rn).s3);
  }
}
#    endif /* _USE_TSPLITPAR */

#  elif (defined BGL && defined XLC)

/**********************************
 *
 * Blue Gene/L Version
 *
 * Author: Carsten Urbach
 *
 **********************************/

/* 6. */
void Hopping_Matrix(const int ieo, spinor * const l, spinor * const k){
  int icx,icy,icz,ioff,ioff2;
  int ix,iy,iz;
  su3 * restrict up ALIGN;
  su3 * restrict um ALIGN;
  spinor * restrict sp ALIGN;
  spinor * restrict sm ALIGN;
  spinor * restrict rn ALIGN;
  /* We have 32 registers available */
  double _Complex reg00, reg01, reg02, reg03, reg04, reg05;
  double _Complex reg10, reg11, reg12, reg13, reg14, reg15;
  /* For the gauge field, reuse the first three!*/
  double _Complex u00, u01, u02, u10, u11, u12;
  double _Complex reg20, reg21;
  /* The following contains the result spinor (12 regs) */
  double _Complex rs00, rs01, rs02, rs10, rs11, rs12, rs20, rs21, rs22, 
    rs30, rs31, rs32;

#pragma disjoint(*sp, *sm, *rn, *up, *um, *l, *k)

  __alignx(16,l);
  __alignx(16,k);

#ifdef _GAUGE_COPY
  if(g_update_gauge_copy) {
    update_backward_gauge();
  }
#endif

#    if (defined MPI && !(defined _NO_COMM))
  xchange_field(k, ieo);
#    endif

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

#    if ((defined _GAUGE_COPY))
  up=&g_gauge_field_copy[ioff][0];
#    else
  up=&g_gauge_field[ix][0];
#    endif
  /**************** loop over all lattice sites ******************/
  for(icx = ioff; icx < (VOLUME/2+ioff); icx++){
    rn=l+(icx-ioff);
    ix=g_eo2lexic[icx];
    /*********************** direction +0 ************************/
    iy=g_idn[ix][0]; 
    icy=g_lexic2eosub[iy];
#    if (!defined _GAUGE_COPY)
    um=&g_gauge_field[iy][0]; 
#    else
    um=up+1;
#    endif
    _prefetch_su3(um); 
    sm=k+icy;
    _prefetch_spinor(sm); 

    _bgl_load_reg0((*sp).s0);
    _bgl_load_reg1((*sp).s1);
    _bgl_load_reg0_up((*sp).s2);
    _bgl_load_reg1_up((*sp).s3);
    _bgl_vector_add_reg0();
    _bgl_vector_add_reg1();
    /* result is now in regx0, regx1, regx2 x = 0,1 */

    _bgl_su3_multiply_double((*up));
    _bgl_vector_cmplx_mul_double(ka0);
    _bgl_store_reg0_up_rs0();
    _bgl_store_reg0_up_rs2();
    _bgl_store_reg1_up_rs1();
    _bgl_store_reg1_up_rs3();

    /*********************** direction -0 ************************/

    iy=g_iup[ix][1]; 
    icy=g_lexic2eosub[iy];

#    if ((defined _GAUGE_COPY))
    up=um+1;
#    else
    up+=1;
#    endif
    _prefetch_su3(up); 
    sp=k+icy;
    _prefetch_spinor(sp); 

    _bgl_load_reg0((*sm).s0);
    _bgl_load_reg1((*sm).s1);
    _bgl_load_reg0_up((*sm).s2);
    _bgl_load_reg1_up((*sm).s3);
    _bgl_vector_sub_reg0();
    _bgl_vector_sub_reg1();

    _bgl_su3_inverse_multiply_double((*um));
    _bgl_vector_cmplxcg_mul_double(ka0);

    _bgl_add_to_rs0_reg0();
    _bgl_sub_from_rs2_reg0();
    _bgl_add_to_rs1_reg1();
    _bgl_sub_from_rs3_reg1();

    /*********************** direction +1 ************************/

    iy=g_idn[ix][1]; 
    icy=g_lexic2eosub[iy];

#    ifndef _GAUGE_COPY
    um=&g_gauge_field[iy][1]; 
#    else
    um = up+1;
#    endif
    _prefetch_su3(um); 
    sm=k+icy;
    _prefetch_spinor(sm); 

    _bgl_load_reg0((*sp).s0);
    _bgl_load_reg1((*sp).s1);
    _bgl_load_reg0_up((*sp).s3);
    _bgl_load_reg1_up((*sp).s2);
    _bgl_vector_i_mul_add_reg0();
    _bgl_vector_i_mul_add_reg1();

    _bgl_su3_multiply_double((*up));
    _bgl_vector_cmplx_mul_double(ka1);

    _bgl_add_to_rs0_reg0();
    _bgl_i_mul_sub_from_rs3_reg0();
    _bgl_add_to_rs1_reg1();
    _bgl_i_mul_sub_from_rs2_reg1();

    /*********************** direction -1 ************************/

    iy=g_iup[ix][2]; 
    icy=g_lexic2eosub[iy];

#    if ((defined _GAUGE_COPY))
    up=um+1;
#    else
    up+=1;
#    endif
    _prefetch_su3(up); 
    sp=k+icy;
    _prefetch_spinor(sp); 

    _bgl_load_reg0((*sm).s0);
    _bgl_load_reg1((*sm).s1);
    _bgl_load_reg0_up((*sm).s3);
    _bgl_load_reg1_up((*sm).s2);
    _bgl_vector_i_mul_sub_reg0();
    _bgl_vector_i_mul_sub_reg1();
      
    _bgl_su3_inverse_multiply_double((*um));
    _bgl_vector_cmplxcg_mul_double(ka1);
      
    _bgl_add_to_rs0_reg0();
    _bgl_add_to_rs1_reg1();
    _bgl_i_mul_add_to_rs3_reg0();
    _bgl_i_mul_add_to_rs2_reg1();      

    /*********************** direction +2 ************************/

    iy=g_idn[ix][2]; 
    icy=g_lexic2eosub[iy];

#    ifndef _GAUGE_COPY
    um=&g_gauge_field[iy][2]; 
#    else
    um= up+1;
#    endif
    _prefetch_su3(um); 
    sm=k+icy;
    _prefetch_spinor(sm); 

    _bgl_load_reg0((*sp).s0);
    _bgl_load_reg1((*sp).s1);
    _bgl_load_reg1_up((*sp).s2);
    _bgl_load_reg0_up((*sp).s3);
    _bgl_vector_add_reg0();
    _bgl_vector_sub_reg1();

    _bgl_su3_multiply_double((*up));
    _bgl_vector_cmplx_mul_double(ka2);

    _bgl_add_to_rs0_reg0();
    _bgl_add_to_rs1_reg1();
    _bgl_sub_from_rs2_reg1();
    _bgl_add_to_rs3_reg0();


    /*********************** direction -2 ************************/

    iy=g_iup[ix][3]; 
    icy=g_lexic2eosub[iy];

#    if ((defined _GAUGE_COPY))
    up=um+1;
#    else
    up+=1;
#    endif
    _prefetch_su3(up); 
    sp=k+icy;
    _prefetch_spinor(sp); 

    _bgl_load_reg0((*sm).s0);
    _bgl_load_reg1((*sm).s1);
    _bgl_load_reg1_up((*sm).s2);
    _bgl_load_reg0_up((*sm).s3);
    _bgl_vector_sub_reg0();
    _bgl_vector_add_reg1();
      
    _bgl_su3_inverse_multiply_double((*um));
    _bgl_vector_cmplxcg_mul_double(ka2);
      
    _bgl_add_to_rs0_reg0();
    _bgl_add_to_rs1_reg1();
    _bgl_add_to_rs2_reg1();
    _bgl_sub_from_rs3_reg0();

    /*********************** direction +3 ************************/

    iy=g_idn[ix][3]; 
    icy=g_lexic2eosub[iy];

#    ifndef _GAUGE_COPY
    um=&g_gauge_field[iy][3]; 
#    else
    um=up+1;
#    endif
    _prefetch_su3(um); 
    sm=k+icy;
    _prefetch_spinor(sm); 

    _bgl_load_reg0((*sp).s0);
    _bgl_load_reg1((*sp).s1);
    _bgl_load_reg0_up((*sp).s2);
    _bgl_load_reg1_up((*sp).s3);
    _bgl_vector_i_mul_add_reg0();
    _bgl_vector_i_mul_sub_reg1();

    _bgl_su3_multiply_double((*up));
    _bgl_vector_cmplx_mul_double(ka3);

    _bgl_add_to_rs0_reg0();
    _bgl_add_to_rs1_reg1();
    _bgl_i_mul_sub_from_rs2_reg0();
    _bgl_i_mul_add_to_rs3_reg1();

    /*********************** direction -3 ************************/

    icz=icx+1;
    if(icz==((VOLUME+RAND)/2+ioff)) icz=ioff;
    iz=g_eo2lexic[icz];
    iy=g_iup[iz][0]; icy=g_lexic2eosub[iy];



#    if ((defined _GAUGE_COPY))
    up=um+1;
#    else
    up=&g_gauge_field[iz][0];
#    endif
    _prefetch_su3(up); 
    sp=k+icy;
    _prefetch_spinor(sp); 

    _bgl_load_reg0((*sm).s0);
    _bgl_load_reg1((*sm).s1);
    _bgl_load_reg0_up((*sm).s2);
    _bgl_load_reg1_up((*sm).s3);
    _bgl_vector_i_mul_sub_reg0();
    _bgl_vector_i_mul_add_reg1();
      
    _bgl_su3_inverse_multiply_double((*um));
    _bgl_vector_cmplxcg_mul_double(ka3);


      
    _bgl_add_to_rs0_reg0();
    _bgl_store_rs0((*rn).s0);
    _bgl_i_mul_add_to_rs2_reg0();
    _bgl_store_rs2((*rn).s2);

    _bgl_add_to_rs1_reg1();
    _bgl_store_rs1((*rn).s1);
    _bgl_i_mul_sub_from_rs3_reg1();
    _bgl_store_rs3((*rn).s3);

    /************************ end of loop ************************/
  }
}
#  elif defined XLC

static su3_vector psi1, psi2, psi, chi, phi1, phi3;

#include"xlc_prefetch.h"

/* 7. */
/* l output , k input*/
/* for ieo=0, k resides on  odd sites and l on even sites */
void Hopping_Matrix(int ieo, spinor * const l, spinor * const k){
  int ix, iy, iz;
  int ioff, ioff2, icx, icy, icz;
  su3 * restrict up, * restrict um;
  spinor * restrict r,* restrict sp,* restrict sm;
  spinor temp;
#pragma disjoint(temp, *sp, *sm, *r, *up, *um)

#ifdef _GAUGE_COPY
  if(g_update_gauge_copy) {
    update_backward_gauge();
  }
#endif

  /* for parallelization */
#    if (defined MPI && !(defined _NO_COMM))
  xchange_field(k, ieo);
#    endif

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
#    if ((defined _GAUGE_COPY))
  up=&g_gauge_field_copy[ioff][0];
#    else
  up=&g_gauge_field[ix][0];
#    endif
  for (icx = ioff; icx < (VOLUME/2 + ioff); icx++) {
    ix = g_eo2lexic[icx];
    /*********************** direction +0 ************************/
    
    iy=g_idn[ix][0]; icy=g_lexic2eosub[iy];
    
    sm=k+icy;
    _prefetch_spinor((void*)sm);

#    ifndef _GAUGE_COPY
    um=&g_gauge_field[iy][0];
#    else
    um=up+1;
#    endif
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
#    if ((defined _GAUGE_COPY))
    up=um+1;
#    else
    up+=1;
#    endif
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

#    ifndef _GAUGE_COPY
    um=&g_gauge_field[iy][1];
#    else
    um=up+1;
#    endif
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
#    if ((defined _GAUGE_COPY))
    up=um+1;
#    else
    up+=1;
#    endif
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

#    ifndef _GAUGE_COPY
    um=&g_gauge_field[iy][2];
#    else
    um=up+1;
#    endif
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
#    if ((defined _GAUGE_COPY))
    up=um+1;
#    else
    up+=1;
#    endif
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

#    ifndef _GAUGE_COPY
    um=&g_gauge_field[iy][3];
#    else
    um=up+1;
#    endif
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
#    if ((defined _GAUGE_COPY))
    up=&g_gauge_field_copy[icz][0];
#    else
    up=&g_gauge_field[iz][0];
#    endif
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
#  else

static su3_vector psi1, psi2, psi, chi, phi1, phi3;

/* 8. */
/* l output , k input*/
/* for ieo=0, k resides on  odd sites and l on even sites */
void Hopping_Matrix(int ieo, spinor * const l, spinor * const k){
  int ix,iy;
  int ioff,ioff2,icx,icy;
  su3 * restrict up, * restrict um;
  spinor * restrict r, * restrict sp, * restrict sm;
  spinor temp;

#ifdef _GAUGE_COPY
  if(g_update_gauge_copy) {
    update_backward_gauge();
  }
#endif

  /* for parallelization */
#    if (defined MPI && !(defined _NO_COMM))
  xchange_field(k, ieo);
#    endif

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
#    if ((defined _GAUGE_COPY))
    up=&g_gauge_field_copy[icx][0];
#    else
    up=&g_gauge_field[ix][0];
#    endif
      
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
#    if ((defined _GAUGE_COPY))
    um = up+1;
#    else
    um=&g_gauge_field[iy][0];
#    endif

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

#    if ((defined _GAUGE_COPY))
    up=um+1;
#    else
    up+=1;
#    endif
      
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
#    ifndef _GAUGE_COPY
    um=&g_gauge_field[iy][1];
#    else
    um=up+1;
#    endif

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
#    if ((defined _GAUGE_COPY))
    up=um+1;
#    else
    up+=1;
#    endif 
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
#    ifndef _GAUGE_COPY
    um = &g_gauge_field[iy][2];
#    else
    um = up +1;
#    endif

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
#    if ((defined _GAUGE_COPY))
    up=um+1;
#    else
    up+=1;
#    endif 
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
#    ifndef _GAUGE_COPY
    um = &g_gauge_field[iy][3];
#    else
    um = up+1;
#    endif

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
/* end of If defined SSE2 */
#  endif


#endif /* thats _USE_HALFSPINOR */


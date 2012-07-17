/**********************************************************************
 *
 *
 * Copyright (C) 2012 Carsten Urbach
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
 **********************************************************************/


void Hopping_Matrix(const int ieo, spinor * const l, spinor * const k){
  int ix;
  su3 * restrict ALIGN U;
  spinor * restrict ALIGN s;
  halfspinor * restrict * phi ALIGN;
  halfspinor32 * restrict * phi32 ALIGN;
  /* We have 32 registers available */
  vector4double ALIGN r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11;
  vector4double ALIGN rs0, rs1, rs2, rs3, rs4, rs5, rs6, rs7, rs8, rs9, rs10, rs11;
  vector4double ALIGN U0, U1, U2, U3, U4, U6, U7;
  vector4double ALIGN rtmp;

#ifdef _KOJAK_INST
#pragma pomp inst begin(hoppingmatrix)
#endif
#pragma disjoint(*s, *U)

#ifdef _GAUGE_COPY
  if(g_update_gauge_copy) {
    update_backward_gauge(g_gauge_field);
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
    for(int i = 0; i < (VOLUME)/2; i++){

      _vec_load2(rs0, rs1, rs2, s->s0);
      _vec_load2(rs3, rs4, rs5, s->s1);
      _vec_load2(rs6, rs7, rs8, s->s2);
      _vec_load2(rs9, rs10, rs11, s->s3);
      s++; 
      _prefetch_spinor(s); 
      /*********************** direction +0 ************************/
      _prefetch_su3(U+1);

      _vec_add_to2(r0, r1, r2, rs0, rs1, rs2, rs6, rs7, rs8);
      _vec_add_to2(r3, r4, r5, rs3, rs4, rs5, rs9, rs10, rs11);
      _vec_su3_multiply_double2((*U));
      _vec_cmplx_mul_double2(r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, U0,ka0);
      
      _vec_store2_32(phi32[ix]->s0, r6, r7, r8);
      _vec_store2_32(phi32[ix]->s1, r9, r10, r11);
      U++;
      ix++;

      /*********************** direction -0 ************************/
      _vec_sub_to2(r0, r1, r2, rs0, rs1, rs2, rs6, rs7, rs8);
      _vec_sub_to2(r3, r4, r5, rs3, rs4, rs5, rs9, rs10, rs11);

      _vec_store2_32(phi32[ix]->s0, r0, r1, r2);
      _vec_store2_32(phi32[ix]->s1, r3, r4, r5);
      ix++;

      /*********************** direction +1 ************************/
      _prefetch_su3(U+1);
      _vec_i_mul_add_to2(r0, r1, r2, rs0, rs1, rs2, rs9, rs10, rs11);
      _vec_i_mul_add_to2(r3, r4, r5, rs3, rs4, rs5, rs6, rs7, rs8);

      _vec_su3_multiply_double2((*U));
      _vec_cmplx_mul_double2(r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, U0,ka1);
      
      _vec_store2_32(phi32[ix]->s0, r6, r7, r8);
      _vec_store2_32(phi32[ix]->s1, r9, r10, r11);
      ix++;
      U++;

      /*********************** direction -1 ************************/
      _vec_i_mul_sub_to2(r0, r1, r2, rs0, rs1, rs2, rs9, rs10, rs11);
      _vec_i_mul_sub_to2(r3, r4, r5, rs3, rs4, rs5, rs6, rs7, rs8);

      _vec_store2_32(phi32[ix]->s0, r0, r1, r2);
      _vec_store2_32(phi32[ix]->s1, r3, r4, r5);
      ix++;


      /*********************** direction +2 ************************/
      _prefetch_su3(U+1);

      _vec_add_to2(r0, r1, r2, rs0, rs1, rs2, rs9, rs10, rs11);
      _vec_sub_to2(r3, r4, r5, rs3, rs4, rs5, rs6, rs7, rs8);
      _vec_su3_multiply_double2((*U));
      _vec_cmplx_mul_double2(r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, U0,ka2);

      _vec_store2_32(phi32[ix]->s0, r6, r7, r8);
      _vec_store2_32(phi32[ix]->s1, r9, r10, r11);
      ix++;
      U++;

      /*********************** direction -2 ************************/
      _vec_sub_to2(r0, r1, r2, rs0, rs1, rs2, rs9, rs10, rs11);
      _vec_sub_to2(r3, r4, r5, rs3, rs4, rs5, rs6, rs7, rs8);

      _vec_store2_32(phi32[ix]->s0, r0, r1, r2);
      _vec_store2_32(phi32[ix]->s1, r3, r4, r5);
      ix++;

      /*********************** direction +3 ************************/
      _prefetch_su3(U+1); 
      _vec_i_mul_add_to2(r0, r1, r2, rs0, rs1, rs2, rs6, rs7, rs8);
      _vec_i_mul_sub_to2(r3, r4, r5, rs3, rs4, rs5, rs9, rs10, rs11);
      _vec_su3_multiply_double2((*U));
      _vec_cmplx_mul_double2(r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, U0,ka3);

      _vec_store2_32(phi32[ix]->s0, r6, r7, r8);
      _vec_store2_32(phi32[ix]->s1, r9, r10, r11);
      ix++;
      U++;

      /*********************** direction -3 ************************/
      _vec_i_mul_sub_to2(r0, r1, r2, rs0, rs1, rs2, rs6, rs7, rs8);
      _vec_i_mul_add_to2(r3, r4, r5, rs3, rs4, rs5, rs9, rs10, rs11);

      _vec_store2_32(phi32[ix]->s0, r0, r1, r2);
      _vec_store2_32(phi32[ix]->s1, r3, r4, r5);
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
    for(int i = 0; i < (VOLUME)/2; i++){
      /* This causes a lot of trouble, do we understand this? */
      /*     _prefetch_spinor_for_store(s); */
      _prefetch_halfspinor(phi32[ix+1]);
      /*********************** direction +0 ************************/
      _vec_load2_32(rs0, rs1, rs2, phi32[ix]->s0);
      rs6 = rs0;
      rs7 = rs1;
      rs8 = rs2;
      _vec_load2_32(rs3, rs4, rs5, phi32[ix]->s1);
      rs9 = rs3;
      rs10= rs4;
      rs11= rs5;
      ix++;
      /*********************** direction -0 ************************/
      _prefetch_su3(U+1);

      _vec_load2_32(r0, r1, r2, phi32[ix]->s0);
      _vec_load2_32(r3, r4, r5, phi32[ix]->s1);
      _vec_su3_inverse_multiply_double2((*U));
      _vec_cmplxcg_mul_double2(r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, U0,ka0);

      _vec_add2(rs0, rs1, rs2, r0, r1, r2);
      _vec_sub2(rs6, rs7, rs8, r0, r1, r2);
      _vec_add2(rs3, rs4, rs5, r3, r4, r5);
      _vec_sub2(rs9, rs10, rs11, r3, r4, r5);
      U++;
      ix++;
      /*********************** direction +1 ************************/
      _vec_load2_32(r0, r1, r2, phi32[ix]->s0);
      _vec_load2_32(r3, r4, r5, phi32[ix]->s1);

      _vec_add2(rs0, rs1, rs2, r0, r1, r2);
      _vec_i_mul_sub2(rs9, rs10, rs11, r0, r1, r2);
      _vec_add2(rs3, rs4, rs5, r3, r4, r5);
      _vec_i_mul_sub2(rs6, rs7, rs8, r3, r4, r5);
      ix++;
      /*********************** direction -1 ************************/
      _prefetch_su3(U+1);
      _vec_load2_32(r0, r1, r2, phi32[ix]->s0);
      _vec_load2_32(r3, r4, r5, phi32[ix]->s1);
      _vec_su3_inverse_multiply_double2((*U));
      _vec_cmplxcg_mul_double2(r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, U0,ka1);

      _vec_add_double2(rs0, rs1, rs2, rs3, rs4, rs5, r0, r1, r2, r3, r4, r5);
      _vec_i_mul_add_double2(rs6, rs7, rs8, rs9, rs10, rs11, r0, r1, r2, r3, r4, r5);
      U++;
      ix++;
      /*********************** direction +2 ************************/
      _vec_load2_32(r0, r1, r2, phi32[ix]->s0);
      _vec_load2_32(r3, r4, r5, phi32[ix]->s1);

      _vec_add_double2(rs0, rs1, rs2, rs3, rs4, rs5, r0, r1, r2, r3, r4, r5);
      _vec_sub2(rs6, rs7, rs8, r3, r4, r5);
      _vec_add2(rs9, rs10, rs11, r0, r1, r2);
      ix++;
      /*********************** direction -2 ************************/
      _prefetch_su3(U+1);

      _vec_load2_32(r0, r1, r2, phi32[ix]->s0);
      _vec_load2_32(r3, r4, r5, phi32[ix]->s1);
      _vec_su3_inverse_multiply_double2((*U));
      _vec_cmplxcg_mul_double2(r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, U0,ka2);

      _vec_add_double2(rs0, rs1, rs2, rs3, rs4, rs5, r0, r1, r2, r3, r4, r5);
      _vec_add2(rs6, rs7, rs8, r3, r4, r5);
      _vec_sub2(rs9, rs10, rs11, r0, r1, r2);
      U++;
      ix++;
      /*********************** direction +3 ************************/
      _vec_load2_32(r0, r1, r2, phi32[ix]->s0);
      _vec_load2_32(r3, r4, r5, phi32[ix]->s1);

      _vec_add_double2(rs0, rs1, rs2, rs3, rs4, rs5, r0, r1, r2, r3, r4, r5);
      _vec_i_mul_sub2(rs6, rs7, rs8, r0, r1, r2);
      _vec_i_mul_add2(rs9, rs10, rs11, r3, r4, r5);
      ix++;
      /*********************** direction -3 ************************/
      _prefetch_su3(U+1);

      _vec_load2_32(r0, r1, r2, phi32[ix]->s0);
      _vec_load2_32(r3, r4, r5, phi32[ix]->s1);
      _vec_su3_inverse_multiply_double2((*U));
      _vec_cmplxcg_mul_double2(r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, U0,ka2);

      _vec_add2(rs0, rs1, rs2, r0, r1, r2);
      _vec_store2(s->s0, rs0, rs1, rs2);
      _vec_i_mul_add2(rs6, rs7, rs8, r0, r1, r2);
      _vec_store2(s->s2, rs6, rs7, rs8);

      _vec_add2(rs3, rs4, rs5, r3, r4, r5);
      _vec_store2(s->s1, rs3, rs4, rs5);
      _vec_i_mul_sub2(rs9, rs10, rs11, r3, r4, r5);
      _vec_store2(s->s3, rs9, rs10, rs11);
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
    for(int i = 0; i < (VOLUME)/2; i++){
      _vec_load2(rs0, rs1, rs2, s->s0);
      _vec_load2(rs3, rs4, rs5, s->s1);
      _vec_load2(rs6, rs7, rs8, s->s2);
      _vec_load2(rs9, rs10, rs11, s->s3);
      s++; 
      _prefetch_spinor(s); 
      /*********************** direction +0 ************************/
      _prefetch_su3(U+1);

      _vec_add_to2(r0, r1, r2, rs0, rs1, rs2, rs6, rs7, rs8);
      _vec_add_to2(r3, r4, r5, rs3, rs4, rs5, rs9, rs10, rs11);
      _vec_su3_multiply_double2((*U));
      _vec_cmplx_mul_double2(r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, U0,ka0);
      
      _vec_store2(phi[ix]->s0, r6, r7, r8);
      _vec_store2(phi[ix]->s1, r9, r10, r11);
      U++;
      ix++;

      /*********************** direction -0 ************************/
      _vec_sub_to2(r0, r1, r2, rs0, rs1, rs2, rs6, rs7, rs8);
      _vec_sub_to2(r3, r4, r5, rs3, rs4, rs5, rs9, rs10, rs11);

      _vec_store2(phi[ix]->s0, r0, r1, r2);
      _vec_store2(phi[ix]->s1, r3, r4, r5);
      ix++;

      /*********************** direction +1 ************************/
      _prefetch_su3(U+1);
      _vec_i_mul_add_to2(r0, r1, r2, rs0, rs1, rs2, rs9, rs10, rs11);
      _vec_i_mul_add_to2(r3, r4, r5, rs3, rs4, rs5, rs6, rs7, rs8);

      _vec_su3_multiply_double2((*U));
      _vec_cmplx_mul_double2(r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, U0,ka1);
      
      _vec_store2(phi[ix]->s0, r6, r7, r8);
      _vec_store2(phi[ix]->s1, r9, r10, r11);
      ix++;
      U++;

      /*********************** direction -1 ************************/
      _vec_i_mul_sub_to2(r0, r1, r2, rs0, rs1, rs2, rs9, rs10, rs11);
      _vec_i_mul_sub_to2(r3, r4, r5, rs3, rs4, rs5, rs6, rs7, rs8);

      _vec_store2(phi[ix]->s0, r0, r1, r2);
      _vec_store2(phi[ix]->s1, r3, r4, r5);
      ix++;


      /*********************** direction +2 ************************/
      _prefetch_su3(U+1);

      _vec_add_to2(r0, r1, r2, rs0, rs1, rs2, rs9, rs10, rs11);
      _vec_sub_to2(r3, r4, r5, rs3, rs4, rs5, rs6, rs7, rs8);
      _vec_su3_multiply_double2((*U));
      _vec_cmplx_mul_double2(r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, U0,ka2);

      _vec_store2(phi[ix]->s0, r6, r7, r8);
      _vec_store2(phi[ix]->s1, r9, r10, r11);
      ix++;
      U++;

      /*********************** direction -2 ************************/
      _vec_sub_to2(r0, r1, r2, rs0, rs1, rs2, rs9, rs10, rs11);
      _vec_sub_to2(r3, r4, r5, rs3, rs4, rs5, rs6, rs7, rs8);

      _vec_store2(phi[ix]->s0, r0, r1, r2);
      _vec_store2(phi[ix]->s1, r3, r4, r5);
      ix++;

      /*********************** direction +3 ************************/
      _prefetch_su3(U+1); 
      _vec_i_mul_add_to2(r0, r1, r2, rs0, rs1, rs2, rs6, rs7, rs8);
      _vec_i_mul_sub_to2(r3, r4, r5, rs3, rs4, rs5, rs9, rs10, rs11);
      _vec_su3_multiply_double2((*U));
      _vec_cmplx_mul_double2(r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, U0,ka3);

      _vec_store2(phi[ix]->s0, r6, r7, r8);
      _vec_store2(phi[ix]->s1, r9, r10, r11);
      ix++;
      U++;

      /*********************** direction -3 ************************/
      _vec_i_mul_sub_to2(r0, r1, r2, rs0, rs1, rs2, rs6, rs7, rs8);
      _vec_i_mul_add_to2(r3, r4, r5, rs3, rs4, rs5, rs9, rs10, rs11);

      _vec_store2(phi[ix]->s0, r0, r1, r2);
      _vec_store2(phi[ix]->s1, r3, r4, r5);
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
    for(int i = 0; i < (VOLUME)/2; i++){
      /* This causes a lot of trouble, do we understand this? */
      /*     _prefetch_spinor_for_store(s); */
      _prefetch_halfspinor(phi[ix+1]);
      /*********************** direction +0 ************************/
      _vec_load2(rs0, rs1, rs2, phi[ix]->s0);
      rs6 = rs0;
      rs7 = rs1;
      rs8 = rs2;
      _vec_load2(rs3, rs4, rs5, phi[ix]->s1);
      rs9 = rs3;
      rs10= rs4;
      rs11= rs5;
      ix++;
      /*********************** direction -0 ************************/
      _prefetch_su3(U+1);

      _vec_load2(r0, r1, r2, phi[ix]->s0);
      _vec_load2(r3, r4, r5, phi[ix]->s1);
      _vec_su3_inverse_multiply_double2((*U));
      _vec_cmplxcg_mul_double2(r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, U0,ka0);

      _vec_add2(rs0, rs1, rs2, r0, r1, r2);
      _vec_sub2(rs6, rs7, rs8, r0, r1, r2);
      _vec_add2(rs3, rs4, rs5, r3, r4, r5);
      _vec_sub2(rs9, rs10, rs11, r3, r4, r5);
      U++;
      ix++;
      /*********************** direction +1 ************************/
      _vec_load2(r0, r1, r2, phi[ix]->s0);
      _vec_load2(r3, r4, r5, phi[ix]->s1);

      _vec_add2(rs0, rs1, rs2, r0, r1, r2);
      _vec_i_mul_sub2(rs9, rs10, rs11, r0, r1, r2);
      _vec_add2(rs3, rs4, rs5, r3, r4, r5);
      _vec_i_mul_sub2(rs6, rs7, rs8, r3, r4, r5);
      ix++;
      /*********************** direction -1 ************************/
      _prefetch_su3(U+1);
      _vec_load2(r0, r1, r2, phi[ix]->s0);
      _vec_load2(r3, r4, r5, phi[ix]->s1);
      _vec_su3_inverse_multiply_double2((*U));
      _vec_cmplxcg_mul_double2(r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, U0,ka1);

      _vec_add_double2(rs0, rs1, rs2, rs3, rs4, rs5, r0, r1, r2, r3, r4, r5);
      _vec_i_mul_add_double2(rs6, rs7, rs8, rs9, rs10, rs11, r0, r1, r2, r3, r4, r5);
      U++;
      ix++;
      /*********************** direction +2 ************************/
      _vec_load2(r0, r1, r2, phi[ix]->s0);
      _vec_load2(r3, r4, r5, phi[ix]->s1);

      _vec_add_double2(rs0, rs1, rs2, rs3, rs4, rs5, r0, r1, r2, r3, r4, r5);
      _vec_sub2(rs6, rs7, rs8, r3, r4, r5);
      _vec_add2(rs9, rs10, rs11, r0, r1, r2);
      ix++;
      /*********************** direction -2 ************************/
      _prefetch_su3(U+1);

      _vec_load2(r0, r1, r2, phi[ix]->s0);
      _vec_load2(r3, r4, r5, phi[ix]->s1);
      _vec_su3_inverse_multiply_double2((*U));
      _vec_cmplxcg_mul_double2(r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, U0,ka2);

      _vec_add_double2(rs0, rs1, rs2, rs3, rs4, rs5, r0, r1, r2, r3, r4, r5);
      _vec_add2(rs6, rs7, rs8, r3, r4, r5);
      _vec_sub2(rs9, rs10, rs11, r0, r1, r2);
      U++;
      ix++;
      /*********************** direction +3 ************************/
      _vec_load2(r0, r1, r2, phi[ix]->s0);
      _vec_load2(r3, r4, r5, phi[ix]->s1);

      _vec_add_double2(rs0, rs1, rs2, rs3, rs4, rs5, r0, r1, r2, r3, r4, r5);
      _vec_i_mul_sub2(rs6, rs7, rs8, r0, r1, r2);
      _vec_i_mul_add2(rs9, rs10, rs11, r3, r4, r5);
      ix++;
      /*********************** direction -3 ************************/
      _prefetch_su3(U+1);

      _vec_load2(r0, r1, r2, phi[ix]->s0);
      _vec_load2(r3, r4, r5, phi[ix]->s1);
      _vec_su3_inverse_multiply_double2((*U));
      _vec_cmplxcg_mul_double2(r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, U0,ka2);

      _vec_add2(rs0, rs1, rs2, r0, r1, r2);
      _vec_store2(s->s0, rs0, rs1, rs2);
      _vec_i_mul_add2(rs6, rs7, rs8, r0, r1, r2);
      _vec_store2(s->s2, rs6, rs7, rs8);

      _vec_add2(rs3, rs4, rs5, r3, r4, r5);
      _vec_store2(s->s1, rs3, rs4, rs5);
      _vec_i_mul_sub2(rs9, rs10, rs11, r3, r4, r5);
      _vec_store2(s->s3, rs9, rs10, rs11);
      U++;
      ix++;
      s++;
    }
  }
#ifdef _KOJAK_INST
#pragma pomp inst end(hoppingmatrix)
#endif
}

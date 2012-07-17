/**********************************************************************
 *
 *
 * Copyright (C) 2003, 2004, 2005, 2006, 2007, 2008 Carsten Urbach
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
    for(i = 0; i < (VOLUME)/2; i++){

      _bgl_load_rs0(s->s0);
      _bgl_load_rs1(s->s1);
      _bgl_load_rs2(s->s2);
      _bgl_load_rs3(s->s3);
      s++; 
      _prefetch_spinor(s); 
      /*********************** direction +0 ************************/
      _prefetch_su3(U+1);

      _bgl_vector_add_rs2_to_rs0_reg0();
      _bgl_vector_add_rs3_to_rs1_reg1();

      _bgl_su3_multiply_double((*U));
      _bgl_vector_cmplx_mul_double(ka0);
      /* result is now in regx0, regx1, regx2 , x=0,1 */

      _bgl_store_reg0_up_32(phi32[ix]->s0);
      _bgl_store_reg1_up_32(phi32[ix]->s1);
      U++;
      ix++;

      /*********************** direction -0 ************************/
      _bgl_vector_sub_rs2_from_rs0_reg0();
      _bgl_vector_sub_rs3_from_rs1_reg1();

      _bgl_store_reg0_32(phi32[ix]->s0);
      _bgl_store_reg1_32(phi32[ix]->s1);
      ix++;

      /*********************** direction +1 ************************/
      _prefetch_su3(U+1);
      _bgl_vector_i_mul_add_rs3_to_rs0_reg0();
      _bgl_vector_i_mul_add_rs2_to_rs1_reg1();

      _bgl_su3_multiply_double((*U));
      _bgl_vector_cmplx_mul_double(ka1);

      _bgl_store_reg0_up_32(phi32[ix]->s0);
      _bgl_store_reg1_up_32(phi32[ix]->s1);
      ix++;
      U++;

      /*********************** direction -1 ************************/
      _bgl_vector_i_mul_sub_rs3_from_rs0_reg0();
      _bgl_vector_i_mul_sub_rs2_from_rs1_reg1();

      _bgl_store_reg0_32(phi32[ix]->s0);
      _bgl_store_reg1_32(phi32[ix]->s1);
      ix++;


      /*********************** direction +2 ************************/
      _prefetch_su3(U+1);

      _bgl_vector_add_rs3_to_rs0_reg0();
      _bgl_vector_sub_rs2_from_rs1_reg1();

      _bgl_su3_multiply_double((*U));
      _bgl_vector_cmplx_mul_double(ka2);

      _bgl_store_reg0_up_32(phi32[ix]->s0);
      _bgl_store_reg1_up_32(phi32[ix]->s1);
      ix++;
      U++;

      /*********************** direction -2 ************************/
      _bgl_vector_sub_rs3_from_rs0_reg0();
      _bgl_vector_add_rs2_to_rs1_reg1();

      _bgl_store_reg0_32(phi32[ix]->s0);
      _bgl_store_reg1_32(phi32[ix]->s1);
      ix++;

      /*********************** direction +3 ************************/
      _prefetch_su3(U+1); 

      _bgl_vector_i_mul_add_rs2_to_rs0_reg0();
      _bgl_vector_i_mul_sub_rs3_from_rs1_reg1();

      _bgl_su3_multiply_double((*U));
      _bgl_vector_cmplx_mul_double(ka3);

      _bgl_store_reg0_up_32(phi32[ix]->s0);
      _bgl_store_reg1_up_32(phi32[ix]->s1);
      ix++;
      U++;

      /*********************** direction -3 ************************/
      _bgl_vector_i_mul_sub_rs2_from_rs0_reg0();
      _bgl_vector_i_mul_add_rs3_to_rs1_reg1();

      _bgl_store_reg0_32(phi32[ix]->s0);
      _bgl_store_reg1_32(phi32[ix]->s1);
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
      _bgl_load_rs0_32(phi32[ix]->s0);
      rs20 = rs00;
      rs21 = rs01;
      rs22 = rs02;
      _bgl_load_rs1_32(phi32[ix]->s1);
      rs30 = rs10;
      rs31 = rs11;
      rs32 = rs12;
      ix++;
      /*********************** direction -0 ************************/
      _prefetch_su3(U+1);

      _bgl_load_reg0_32(phi32[ix]->s0);
      _bgl_load_reg1_32(phi32[ix]->s1);

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
      _bgl_load_reg0_up_32(phi32[ix]->s0);
      _bgl_load_reg1_up_32(phi32[ix]->s1);

      _bgl_add_to_rs0_reg0();
      _bgl_i_mul_sub_from_rs3_reg0();
      _bgl_add_to_rs1_reg1();
      _bgl_i_mul_sub_from_rs2_reg1();
      ix++;
      /*********************** direction -1 ************************/
      _prefetch_su3(U+1);

      _bgl_load_reg0_32(phi32[ix]->s0);
      _bgl_load_reg1_32(phi32[ix]->s1);

      _bgl_su3_inverse_multiply_double((*U));
      _bgl_vector_cmplxcg_mul_double(ka1);

      _bgl_add_to_rs0_reg0();
      _bgl_add_to_rs1_reg1();
      _bgl_i_mul_add_to_rs3_reg0();
      _bgl_i_mul_add_to_rs2_reg1();      
      U++;
      ix++;
      /*********************** direction +2 ************************/
      _bgl_load_reg0_up_32(phi32[ix]->s0);
      _bgl_load_reg1_up_32(phi32[ix]->s1);

      _bgl_add_to_rs0_reg0();
      _bgl_add_to_rs1_reg1();
      _bgl_sub_from_rs2_reg1();
      _bgl_add_to_rs3_reg0();
      ix++;
      /*********************** direction -2 ************************/
      _prefetch_su3(U+1);

      _bgl_load_reg0_32(phi32[ix]->s0);
      _bgl_load_reg1_32(phi32[ix]->s1);

      _bgl_su3_inverse_multiply_double((*U));
      _bgl_vector_cmplxcg_mul_double(ka2);

      _bgl_add_to_rs0_reg0();
      _bgl_add_to_rs1_reg1();
      _bgl_add_to_rs2_reg1();
      _bgl_sub_from_rs3_reg0();
      U++;
      ix++;
      /*********************** direction +3 ************************/
      _bgl_load_reg0_up_32(phi32[ix]->s0);
      _bgl_load_reg1_up_32(phi32[ix]->s1);

      _bgl_add_to_rs0_reg0();
      _bgl_add_to_rs1_reg1();
      _bgl_i_mul_sub_from_rs2_reg0();
      _bgl_i_mul_add_to_rs3_reg1();
      ix++;
      /*********************** direction -3 ************************/
      _prefetch_su3(U+1);

      _bgl_load_reg0_32(phi32[ix]->s0);
      _bgl_load_reg1_32(phi32[ix]->s1);

      _bgl_su3_inverse_multiply_double((*U));
      _bgl_vector_cmplxcg_mul_double(ka3);

      _bgl_add_to_rs0_reg0();
      _bgl_store_rs0(s->s0);
      _bgl_i_mul_add_to_rs2_reg0();
      _bgl_store_rs2(s->s2);

      _bgl_add_to_rs1_reg1();
      _bgl_store_rs1(s->s1);
      _bgl_i_mul_sub_from_rs3_reg1();
      _bgl_store_rs3(s->s3);

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
      _bgl_load_rs0(s->s0);
      _bgl_load_rs1(s->s1);
      _bgl_load_rs2(s->s2);
      _bgl_load_rs3(s->s3);
      s++; 
      _prefetch_spinor(s); 
      /*********************** direction +0 ************************/
      _prefetch_su3(U+1);

      _bgl_vector_add_rs2_to_rs0_reg0();
      _bgl_vector_add_rs3_to_rs1_reg1();

      _bgl_su3_multiply_double((*U));
      _bgl_vector_cmplx_mul_double(ka0);
      /* result is now in regx0, regx1, regx2 , x=0,1 */

      _bgl_store_reg0_up(phi[ix]->s0);
      _bgl_store_reg1_up(phi[ix]->s1);
      U++;
      ix++;

      /*********************** direction -0 ************************/
      _prefetch_halfspinor(phi[ix+4]);
      _bgl_vector_sub_rs2_from_rs0_reg0();
      _bgl_vector_sub_rs3_from_rs1_reg1();

      _bgl_store_reg0(phi[ix]->s0);
      _bgl_store_reg1(phi[ix]->s1);
      ix++;

      /*********************** direction +1 ************************/
      _prefetch_halfspinor(phi[ix+4]);
      _prefetch_su3(U+1);
      _bgl_vector_i_mul_add_rs3_to_rs0_reg0();
      _bgl_vector_i_mul_add_rs2_to_rs1_reg1();

      _bgl_su3_multiply_double((*U));
      _bgl_vector_cmplx_mul_double(ka1);

      _bgl_store_reg0_up(phi[ix]->s0);
      _bgl_store_reg1_up(phi[ix]->s1);
      ix++;
      U++;

      /*********************** direction -1 ************************/
      _prefetch_halfspinor(phi[ix+4]);
      _bgl_vector_i_mul_sub_rs3_from_rs0_reg0();
      _bgl_vector_i_mul_sub_rs2_from_rs1_reg1();

      _bgl_store_reg0(phi[ix]->s0);
      _bgl_store_reg1(phi[ix]->s1);
      ix++;


      /*********************** direction +2 ************************/
      _prefetch_halfspinor(phi[ix+4]);
      _prefetch_su3(U+1);

      _bgl_vector_add_rs3_to_rs0_reg0();
      _bgl_vector_sub_rs2_from_rs1_reg1();

      _bgl_su3_multiply_double((*U));
      _bgl_vector_cmplx_mul_double(ka2);

      _bgl_store_reg0_up(phi[ix]->s0);
      _bgl_store_reg1_up(phi[ix]->s1);
      ix++;
      U++;

      /*********************** direction -2 ************************/
      _prefetch_halfspinor(phi[ix+4]);
      _bgl_vector_sub_rs3_from_rs0_reg0();
      _bgl_vector_add_rs2_to_rs1_reg1();

      _bgl_store_reg0(phi[ix]->s0);
      _bgl_store_reg1(phi[ix]->s1);
      ix++;

      /*********************** direction +3 ************************/
      _prefetch_halfspinor(phi[ix+4]);
      _prefetch_su3(U+1); 

      _bgl_vector_i_mul_add_rs2_to_rs0_reg0();
      _bgl_vector_i_mul_sub_rs3_from_rs1_reg1();

      _bgl_su3_multiply_double((*U));
      _bgl_vector_cmplx_mul_double(ka3);

      _bgl_store_reg0_up(phi[ix]->s0);
      _bgl_store_reg1_up(phi[ix]->s1);
      ix++;
      U++;

      /*********************** direction -3 ************************/
      _prefetch_halfspinor(phi[ix+4]);
      _bgl_vector_i_mul_sub_rs2_from_rs0_reg0();
      _bgl_vector_i_mul_add_rs3_to_rs1_reg1();

      _bgl_store_reg0(phi[ix]->s0);
      _bgl_store_reg1(phi[ix]->s1);
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
      _bgl_load_rs0(phi[ix]->s0);
      rs20 = rs00;
      rs21 = rs01;
      rs22 = rs02;
      _bgl_load_rs1(phi[ix]->s1);
      rs30 = rs10;
      rs31 = rs11;
      rs32 = rs12;
      ix++;
      /*********************** direction -0 ************************/
      _prefetch_halfspinor(phi[ix+3]);
      _prefetch_su3(U+1);

      _bgl_load_reg0(phi[ix]->s0);
      _bgl_load_reg1(phi[ix]->s1);

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
      _bgl_load_reg0_up(phi[ix]->s0);
      _bgl_load_reg1_up(phi[ix]->s1);

      _bgl_add_to_rs0_reg0();
      _bgl_i_mul_sub_from_rs3_reg0();
      _bgl_add_to_rs1_reg1();
      _bgl_i_mul_sub_from_rs2_reg1();
      ix++;
      /*********************** direction -1 ************************/
      _prefetch_halfspinor(phi[ix+3]);
      _prefetch_su3(U+1);

      _bgl_load_reg0(phi[ix]->s0);
      _bgl_load_reg1(phi[ix]->s1);

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
      _bgl_load_reg0_up(phi[ix]->s0);
      _bgl_load_reg1_up(phi[ix]->s1);

      _bgl_add_to_rs0_reg0();
      _bgl_add_to_rs1_reg1();
      _bgl_sub_from_rs2_reg1();
      _bgl_add_to_rs3_reg0();
      ix++;
      /*********************** direction -2 ************************/

      _prefetch_halfspinor(phi[ix+3]);
      _prefetch_su3(U+1);

      _bgl_load_reg0(phi[ix]->s0);
      _bgl_load_reg1(phi[ix]->s1);

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
      _bgl_load_reg0_up(phi[ix]->s0);
      _bgl_load_reg1_up(phi[ix]->s1);

      _bgl_add_to_rs0_reg0();
      _bgl_add_to_rs1_reg1();
      _bgl_i_mul_sub_from_rs2_reg0();
      _bgl_i_mul_add_to_rs3_reg1();
      ix++;
      /*********************** direction -3 ************************/
      _prefetch_spinor(s);
      _prefetch_halfspinor(phi[ix+3]);
      _prefetch_su3(U+1);

      _bgl_load_reg0(phi[ix]->s0);
      _bgl_load_reg1(phi[ix]->s1);

      _bgl_su3_inverse_multiply_double((*U));
      _bgl_vector_cmplxcg_mul_double(ka3);

      _bgl_add_to_rs0_reg0();
      _bgl_store_rs0(s->s0);
      _bgl_i_mul_add_to_rs2_reg0();
      _bgl_store_rs2(s->s2);

      _bgl_add_to_rs1_reg1();
      _bgl_store_rs1(s->s1);
      _bgl_i_mul_sub_from_rs3_reg1();
      _bgl_store_rs3(s->s3);

      U++;
      ix++;
      s++;
    }
  }
#ifdef _KOJAK_INST
#pragma pomp inst end(hoppingmatrix)
#endif
}


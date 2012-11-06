/**********************************************************************
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

#ifndef _HOPPING_H
#define _HOPPING_H

#  if (defined BGQ && defined XLC)

/* We have 32 registers available */
#define _declare_regs()							\
  vector4double ALIGN r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11; \
  vector4double ALIGN rs0, rs1, rs2, rs3, rs4, rs5, rs6, rs7, rs8, rs9, rs10, rs11; \
  vector4double ALIGN U0, U1, U2, U3, U4, U6, U7;			\
  vector4double ALIGN rtmp;							\
  __alignx(16,l);							\
  __alignx(16,k);

#define _hop_t_p()							\
  _vec_load_spinor(r4, r5, r6, r7, r8, r9, sp->s0);			\
  _vec_add_ul_spinor(r0, r1, r2, r4, r5, r6, r7, r8, r9);		\
  _vec_su3_multiply_double2ct(up);					\
  rtmp = vec_ld2(0, (double*) &ka0);					\
  _vec_cmplx_mul_double2c(rs0, rs1, rs2, r4, r5, r6, rtmp);		\
  _vec_unfuse(rs0, rs1, rs2, rs3, rs4, rs5);				\
  rs6 = rs0; rs7 = rs1; rs8 = rs2;					\
  rs9 = rs3; rs10= rs4; rs11= rs5;

#define _hop_t_m()							\
  _vec_load_spinor(r4, r5, r6, r7, r8, r9, sm->s0);			\
  _vec_sub_ul_spinor(r0, r1, r2, r4, r5, r6, r7, r8, r9);		\
  _vec_su3_inverse_multiply_double2ct(um);				\
  _vec_cmplxcg_mul_double2c(r0, r1, r2, r4, r5, r6, rtmp);		\
  _vec_unfuse(r0, r1, r2, r3, r4, r5);					\
  _vec_add_double2(rs0, rs1, rs2, rs3, rs4, rs5, r0, r1, r2, r3, r4, r5); \
  _vec_sub_double2(rs6, rs7, rs8, rs9, rs10, rs11, r0, r1, r2, r3, r4, r5);

#define _hop_x_p()							\
  _vec_load(r4, r5, sp->s0);						\
  _vec_load16(r6, r7, sp->s1, U0);					\
  _vec_load(r10, r11, sp->s2);						\
  _vec_load16(r8, r9, sp->s3, U0);					\
  _vec_i_mul_add(r0, r1, r4, r5, r8, r9, U0);				\
  _vec_i_mul_add(r2, r3, r6, r7, r10, r11, U0);				\
  _vec_su3_multiply_double2c(up);					\
  rtmp = vec_ld2(0, (double*) &ka1);					\
  _vec_cmplx_mul_double2c(r0, r1, r2, r4, r5, r6, rtmp);		\
  _vec_unfuse(r0, r1, r2, r3, r4, r5);					\
  _vec_add_double2(rs0, rs1, rs2, rs3, rs4, rs5, r0, r1, r2, r3, r4, r5); \
  _vec_i_mul_sub2(rs6, rs7, rs8, r3, r4, r5, U0);			\
  _vec_i_mul_sub2(rs9, rs10, rs11, r0, r1, r2, U1);

#define _hop_x_m()							\
  _vec_load(r4, r5, sm->s0);						\
  _vec_load16(r6, r7, sm->s1, U0);					\
  _vec_load(r10, r11, sm->s2);						\
  _vec_load16(r8, r9, sm->s3, U0);					\
  _vec_i_mul_sub(r0, r1, r4, r5, r8, r9, U0);				\
  _vec_i_mul_sub(r2, r3, r6, r7, r10, r11, U0);				\
  _vec_su3_inverse_multiply_double2c(um);				\
  _vec_cmplxcg_mul_double2c(r0, r1, r2, r4, r5, r6, rtmp);		\
  _vec_unfuse(r0, r1, r2, r3, r4, r5);					\
  _vec_add_double2(rs0, rs1, rs2, rs3, rs4, rs5, r0, r1, r2, r3, r4, r5); \
  _vec_i_mul_add2(rs6, rs7, rs8, r3, r4, r5, U0);			\
  _vec_i_mul_add2(rs9, rs10, rs11, r0, r1, r2, U1);

#define _hop_y_p()							\
  _vec_load(r4, r5, sp->s0);						\
  _vec_load16(r6, r7, sp->s1, U0);					\
  _vec_load(r10, r11, sp->s2);						\
  _vec_load16(r8, r9, sp->s3, U0);					\
  _vec_add(r0, r1, r4, r5, r8, r9);					\
  _vec_sub(r2, r3, r6, r7, r10, r11);					\
  _vec_su3_multiply_double2c(up);					\
  rtmp = vec_ld2(0, (double*) &ka2);					\
  _vec_cmplx_mul_double2c(r0, r1, r2, r4, r5, r6, rtmp);		\
  _vec_unfuse(r0, r1, r2, r3, r4, r5);					\
  _vec_add_double2(rs0, rs1, rs2, rs3, rs4, rs5, r0, r1, r2, r3, r4, r5);	\
  _vec_sub2(rs6, rs7, rs8, r3, r4, r5);					\
  _vec_add2(rs9, rs10, rs11, r0, r1, r2);

#define _hop_y_m()							\
  _vec_load(r4, r5, sm->s0);						\
  _vec_load16(r6, r7, sm->s1, U0);					\
  _vec_load(r10, r11, sm->s2);						\
  _vec_load16(r8, r9, sm->s3, U0);					\
  _vec_sub(r0, r1, r4, r5, r8, r9);					\
  _vec_add(r2, r3, r6, r7, r10, r11);					\
  _vec_su3_inverse_multiply_double2c(um);				\
  _vec_cmplxcg_mul_double2c(r0, r1, r2, r4, r5, r6, rtmp);		\
  _vec_unfuse(r0, r1, r2, r3, r4, r5);					\
  _vec_add_double2(rs0, rs1, rs2, rs3, rs4, rs5, r0, r1, r2, r3, r4, r5); \
  _vec_add2(rs6, rs7, rs8, r3, r4, r5);					\
  _vec_sub2(rs9, rs10, rs11, r0, r1, r2);

#define _hop_z_p()							\
  _vec_load(r4, r5, sp->s0);						\
  _vec_load16(r6, r7, sp->s1, U0);					\
  _vec_load(r8, r9, sp->s2);						\
  _vec_load16(r10, r11, sp->s3, U0);					\
  _vec_i_mul_add(r0, r1, r4, r5, r8, r9, U0);				\
  _vec_i_mul_sub(r2, r3, r6, r7, r10, r11, U1);				\
  _vec_su3_multiply_double2c(up);					\
  rtmp = vec_ld2(0, (double*) &ka3);					\
  _vec_cmplx_mul_double2c(r0, r1, r2, r4, r5, r6, rtmp);		\
  _vec_unfuse(r0, r1, r2, r3, r4, r5);					\
  _vec_add_double2(rs0, rs1, rs2, rs3, rs4, rs5, r0, r1, r2, r3, r4, r5); \
  _vec_i_mul_sub2(rs6, rs7, rs8, r0, r1, r2, U0);			\
  _vec_i_mul_add2(rs9, rs10, rs11, r3, r4, r5, U1);

#define _hop_z_m()							\
  _vec_load(r4, r5, sm->s0);						\
  _vec_load16(r6, r7, sm->s1, U0);					\
  _vec_load(r8, r9, sm->s2);						\
  _vec_load16(r10, r11, sm->s3, U0);					\
  _vec_i_mul_sub(r0, r1, r4, r5, r8, r9, U0);				\
  _vec_i_mul_add(r2, r3, r6, r7, r10, r11, U1);				\
  _vec_su3_inverse_multiply_double2c(um);				\
  _vec_cmplxcg_mul_double2c(r0, r1, r2, r4, r5, r6, rtmp);		\
  _vec_unfuse(r0, r1, r2, r3, r4, r5);					\
  _vec_add_double2(rs0, rs1, rs2, rs3, rs4, rs5, r0, r1, r2, r3, r4, r5); \
  _vec_i_mul_add2(rs6, rs7, rs8, r0, r1, r2, U0);			\
  _vec_i_mul_sub2(rs9, rs10, rs11, r3, r4, r5, U1);

#define _hop_mul_g5_cmplx_and_store()					\
  _vec_cmplx_mul_double2(r0, r1, r2, r3, r4, r5, rs0, rs1, rs2, rs3, rs4, rs5, cf); \
  _vec_cmplxcg_mul_double2(r6, r7, r8, r9, r10, r11, rs6, rs7, rs8, rs9, rs10, rs11, cf); \
  _vec_store2(rn->s0, r0, r1, r2);					\
  _vec_store2(rn->s1, r3, r4, r5);					\
  _vec_store2(rn->s2, r6, r7, r8);					\
  _vec_store2(rn->s3, r9, r10, r11);

#define _g5_cmplx_sub_hop_and_g5store()					\
  _vec_load_halfspinor(r3, r4, r5, pn->s0);				\
  _vec_cmplx_mul_double2c(r0, r1, r2, r3, r4, r5, cf);			\
  _vec_unfuse(r0, r1, r2, r3, r4, r5);					\
  _vec_sub_double2(r0, r3, r1, r4, r2, r5, rs0, rs1, rs2, rs3, rs4, rs5); \
  _vec_store2(rn->s0, r0, r3, r1);					\
  _vec_store2(rn->s1, r4, r2, r5);					\
  _vec_load_halfspinor(r3, r4, r5, pn->s2);				\
  _vec_cmplxcg_mul_double2c(r0, r1, r2, r3, r4, r5, cf);		\
  _vec_unfuse(r0, r1, r2, r3, r4, r5);					\
  _vec_sub_double2(rs6, rs7, rs8, rs9, rs10, rs11, r0, r3, r1, r4, r2, r5); \
  _vec_store2(rn->s2, rs6, rs7, rs8);					\
  _vec_store2(rn->s3, rs9, rs10, rs11);
  

#define _store_res()				\
  _vec_store2(rn->s0, rs0, rs1, rs2);		\
  _vec_store2(rn->s1, rs3, rs4, rs5);		\
  _vec_store2(rn->s2, rs6, rs7, rs8);		\
  _vec_store2(rn->s3, rs9, rs10, rs11);

#  elif (defined BGL && defined XLC)

#define _declare_regs()						\
  double _Complex reg00, reg01, reg02, reg03, reg04, reg05;	\
  double _Complex reg10, reg11, reg12, reg13, reg14, reg15;	\
  double _Complex u00, u01, u02, u10, u11, u12;			\
  double _Complex reg20, reg21;						\
  double _Complex rs00, rs01, rs02, rs10, rs11, rs12, rs20, rs21, rs22, \
    rs30, rs31, rs32;

#define _hop_t_p()				\
  _prefetch_su3(um);				\
  _prefetch_spinor(sm);				\
  _bgl_load_reg0(sp->s0);			\
  _bgl_load_reg1(sp->s1);			\
  _bgl_load_reg0_up(sp->s2);			\
  _bgl_load_reg1_up(sp->s3);			\
  _bgl_vector_add_reg0();			\
  _bgl_vector_add_reg1();			\
  _bgl_su3_multiply_double((*up));		\
  _bgl_vector_cmplx_mul_double(ka0);		\
  _bgl_store_reg0_up_rs0();			\
  _bgl_store_reg0_up_rs2();			\
  _bgl_store_reg1_up_rs1();			\
  _bgl_store_reg1_up_rs3();

#define _hop_t_m()				\
  _prefetch_su3(up);				\
  _prefetch_spinor(sp);				\
  _bgl_load_reg0(sm->s0);			\
  _bgl_load_reg1(sm->s1);			\
  _bgl_load_reg0_up(sm->s2);			\
  _bgl_load_reg1_up(sm->s3);			\
  _bgl_vector_sub_reg0();			\
  _bgl_vector_sub_reg1();			\
  _bgl_su3_inverse_multiply_double((*um));	\
  _bgl_vector_cmplxcg_mul_double(ka0);		\
  _bgl_add_to_rs0_reg0();			\
  _bgl_sub_from_rs2_reg0();			\
  _bgl_add_to_rs1_reg1();			\
  _bgl_sub_from_rs3_reg1();

#define _hop_x_p()				\
  _prefetch_su3(um);				\
  _prefetch_spinor(sm);				\
  _bgl_load_reg0(sp->s0);			\
  _bgl_load_reg1(sp->s1);			\
  _bgl_load_reg0_up(sp->s3);			\
  _bgl_load_reg1_up(sp->s2);			\
  _bgl_vector_i_mul_add_reg0();			\
  _bgl_vector_i_mul_add_reg1();			\
  _bgl_su3_multiply_double((*up));		\
  _bgl_vector_cmplx_mul_double(ka1);		\
  _bgl_add_to_rs0_reg0();			\
  _bgl_i_mul_sub_from_rs3_reg0();		\
  _bgl_add_to_rs1_reg1();			\
  _bgl_i_mul_sub_from_rs2_reg1();

#define _hop_x_m()				\
  _prefetch_su3(up);				\
  _prefetch_spinor(sp);				\
  _bgl_load_reg0(sm->s0);			\
  _bgl_load_reg1(sm->s1);			\
  _bgl_load_reg0_up(sm->s3);			\
  _bgl_load_reg1_up(sm->s2);			\
  _bgl_vector_i_mul_sub_reg0();			\
  _bgl_vector_i_mul_sub_reg1();			\
  _bgl_su3_inverse_multiply_double((*um));	\
  _bgl_vector_cmplxcg_mul_double(ka1);		\
  _bgl_add_to_rs0_reg0();			\
  _bgl_add_to_rs1_reg1();			\
  _bgl_i_mul_add_to_rs3_reg0();			\
  _bgl_i_mul_add_to_rs2_reg1();      

#define _hop_y_p()				\
  _prefetch_su3(um);				\
  _prefetch_spinor(sm);				\
  _bgl_load_reg0(sp->s0);			\
  _bgl_load_reg1(sp->s1);			\
  _bgl_load_reg1_up(sp->s2);			\
  _bgl_load_reg0_up(sp->s3);			\
  _bgl_vector_add_reg0();			\
  _bgl_vector_sub_reg1();			\
  _bgl_su3_multiply_double((*up));		\
  _bgl_vector_cmplx_mul_double(ka2);		\
  _bgl_add_to_rs0_reg0();			\
  _bgl_add_to_rs1_reg1();			\
  _bgl_sub_from_rs2_reg1();			\
  _bgl_add_to_rs3_reg0();
  
#define _hop_y_m()				\
  _prefetch_su3(up);				\
  _prefetch_spinor(sp);				\
  _bgl_load_reg0(sm->s0);			\
  _bgl_load_reg1(sm->s1);			\
  _bgl_load_reg1_up(sm->s2);			\
  _bgl_load_reg0_up(sm->s3);			\
  _bgl_vector_sub_reg0();			\
  _bgl_vector_add_reg1();			\
  _bgl_su3_inverse_multiply_double((*um));	\
  _bgl_vector_cmplxcg_mul_double(ka2);		\
  _bgl_add_to_rs0_reg0();			\
  _bgl_add_to_rs1_reg1();			\
  _bgl_add_to_rs2_reg1();			\
  _bgl_sub_from_rs3_reg0();

#define _hop_z_p()				\
  _prefetch_su3(um);				\
  _prefetch_spinor(sm);				\
  _bgl_load_reg0(sp->s0);			\
  _bgl_load_reg1(sp->s1);			\
  _bgl_load_reg0_up(sp->s2);			\
  _bgl_load_reg1_up(sp->s3);			\
  _bgl_vector_i_mul_add_reg0();			\
  _bgl_vector_i_mul_sub_reg1();			\
  _bgl_su3_multiply_double((*up));		\
  _bgl_vector_cmplx_mul_double(ka3);		\
  _bgl_add_to_rs0_reg0();			\
  _bgl_add_to_rs1_reg1();			\
  _bgl_i_mul_sub_from_rs2_reg0();		\
  _bgl_i_mul_add_to_rs3_reg1();

#define _hop_z_m()				\
  _prefetch_su3(up);				\
  _prefetch_spinor(sp);				\
  _bgl_load_reg0(sm->s0);			\
  _bgl_load_reg1(sm->s1);			\
  _bgl_load_reg0_up(sm->s2);			\
  _bgl_load_reg1_up(sm->s3);			\
  _bgl_vector_i_mul_sub_reg0();			\
  _bgl_vector_i_mul_add_reg1();			\
  _bgl_su3_inverse_multiply_double((*um));	\
  _bgl_vector_cmplxcg_mul_double(ka3);		\
  _bgl_add_to_rs0_reg0();			\
  _bgl_i_mul_add_to_rs2_reg0();			\
  _bgl_add_to_rs1_reg1();			\
  _bgl_i_mul_sub_from_rs3_reg1();

#define _store_res()				\
  _bgl_store_rs0(rn->s0);			\
  _bgl_store_rs1(rn->s1);			\
  _bgl_store_rs2(rn->s2);			\
  _bgl_store_rs3(rn->s3);

#  elif (defined SSE2 || defined SSE3)

#define _declare_regs()				\
  spinor ALIGN rs;

#define _hop_t_p()				\
  _prefetch_su3(um);				\
  _sse_load(sp->s0);				\
  _sse_load_up(sp->s2);				\
  _sse_vector_add();				\
  _sse_su3_multiply((*up));			\
  _sse_vector_cmplx_mul(ka0);			\
  _sse_store_up(rs.s0);				\
  _sse_store_up(rs.s2);				\
  _sse_load(sp->s1);				\
  _sse_load_up(sp->s3);				\
  _sse_vector_add();				\
  _sse_su3_multiply((*up));			\
  _sse_vector_cmplx_mul(ka0);			\
  _sse_store_up(rs.s1);				\
  _sse_store_up(rs.s3);

#define _hop_t_m()				\
  _prefetch_su3(up);				\
  _sse_load(sm->s0);				\
  _sse_load_up(sm->s2);				\
  _sse_vector_sub();				\
  _sse_su3_inverse_multiply((*um));		\
  _sse_vector_cmplxcg_mul(ka0);			\
  _sse_load(rs.s0);				\
  _sse_vector_add();				\
  _sse_store(rs.s0);				\
  _sse_load(rs.s2);				\
  _sse_vector_sub();				\
  _sse_store(rs.s2);				\
  _sse_load(sm->s1);				\
  _sse_load_up(sm->s3);				\
  _sse_vector_sub();				\
  _sse_su3_inverse_multiply((*um));		\
  _sse_vector_cmplxcg_mul(ka0);			\
  _sse_load(rs.s1);				\
  _sse_vector_add();				\
  _sse_store(rs.s1);				\
  _sse_load(rs.s3);				\
  _sse_vector_sub();				\
  _sse_store(rs.s3);

#define _hop_x_p()				\
  _prefetch_su3(um);				\
  _sse_load(sp->s0);				\
  _sse_load_up(sp->s3);				\
  _sse_vector_i_mul();				\
  _sse_vector_add();				\
  _sse_su3_multiply((*up));			\
  _sse_vector_cmplx_mul(ka1);			\
  _sse_load(rs.s0);				\
  _sse_vector_add();				\
  _sse_store(rs.s0);				\
  _sse_load(rs.s3);				\
  _sse_vector_i_mul();				\
  _sse_vector_sub();				\
  _sse_store(rs.s3);				\
  _sse_load(sp->s1);				\
  _sse_load_up(sp->s2);				\
  _sse_vector_i_mul();				\
  _sse_vector_add();				\
  _sse_su3_multiply((*up));			\
  _sse_vector_cmplx_mul(ka1);			\
  _sse_load(rs.s1);				\
  _sse_vector_add();				\
  _sse_store(rs.s1);				\
  _sse_load(rs.s2);				\
  _sse_vector_i_mul();				\
  _sse_vector_sub();				\
  _sse_store(rs.s2);

#define _hop_x_m()				\
  _prefetch_su3(up);				\
  _sse_load(sm->s0);				\
  _sse_load_up(sm->s3);				\
  _sse_vector_i_mul();				\
  _sse_vector_sub();				\
  _sse_su3_inverse_multiply((*um));		\
  _sse_vector_cmplxcg_mul(ka1);			\
  _sse_load(rs.s0);				\
  _sse_vector_add();				\
  _sse_store(rs.s0);				\
  _sse_load(rs.s3);				\
  _sse_vector_i_mul();				\
  _sse_vector_add();				\
  _sse_store(rs.s3);				\
  _sse_load(sm->s1);				\
  _sse_load_up(sm->s2);				\
  _sse_vector_i_mul();				\
  _sse_vector_sub();				\
  _sse_su3_inverse_multiply((*um));		\
  _sse_vector_cmplxcg_mul(ka1);			\
  _sse_load(rs.s1);				\
  _sse_vector_add();				\
  _sse_store(rs.s1);				\
  _sse_load(rs.s2);				\
  _sse_vector_i_mul();				\
  _sse_vector_add();				\
  _sse_store(rs.s2);

#define _hop_y_p()				\
  _prefetch_su3(um);				\
  _sse_load(sp->s0);				\
  _sse_load_up(sp->s3);				\
  _sse_vector_add();				\
  _sse_su3_multiply((*up));			\
  _sse_vector_cmplx_mul(ka2);			\
  _sse_load(rs.s0);				\
  _sse_vector_add();				\
  _sse_store(rs.s0);				\
  _sse_load(rs.s3);				\
  _sse_vector_add();				\
  _sse_store(rs.s3);				\
  _sse_load(sp->s1);				\
  _sse_load_up(sp->s2);				\
  _sse_vector_sub();				\
  _sse_su3_multiply((*up));			\
  _sse_vector_cmplx_mul(ka2);			\
  _sse_load(rs.s1);				\
  _sse_vector_add();				\
  _sse_store(rs.s1);				\
  _sse_load(rs.s2);				\
  _sse_vector_sub();				\
  _sse_store(rs.s2);      

#define _hop_y_m()				\
  _prefetch_su3(up);				\
  _sse_load(sm->s0);				\
  _sse_load_up(sm->s3);				\
  _sse_vector_sub();				\
  _sse_su3_inverse_multiply((*um));		\
  _sse_vector_cmplxcg_mul(ka2);			\
  _sse_load(rs.s0);				\
  _sse_vector_add();				\
  _sse_store(rs.s0);				\
  _sse_load(rs.s3);				\
  _sse_vector_sub();				\
  _sse_store(rs.s3);				\
  _sse_load(sm->s1);				\
  _sse_load_up(sm->s2);				\
  _sse_vector_add();				\
  _sse_su3_inverse_multiply((*um));		\
  _sse_vector_cmplxcg_mul(ka2);			\
  _sse_load(rs.s1);				\
  _sse_vector_add();				\
  _sse_store(rs.s1);				\
  _sse_load(rs.s2);				\
  _sse_vector_add();				\
  _sse_store(rs.s2);

#define _hop_z_p()				\
  _prefetch_su3(um);				\
  _sse_load(sp->s0);				\
  _sse_load_up(sp->s2);				\
  _sse_vector_i_mul();				\
  _sse_vector_add();				\
  _sse_su3_multiply((*up));			\
  _sse_vector_cmplx_mul(ka3);			\
  _sse_load(rs.s0);				\
  _sse_vector_add();				\
  _sse_store(rs.s0);				\
  _sse_load(rs.s2);				\
  _sse_vector_i_mul();				\
  _sse_vector_sub();				\
  _sse_store(rs.s2);				\
  _sse_load(sp->s1);				\
  _sse_load_up(sp->s3);				\
  _sse_vector_i_mul();				\
  _sse_vector_sub();				\
  _sse_su3_multiply((*up));			\
  _sse_vector_cmplx_mul(ka3);			\
  _sse_load(rs.s1);				\
  _sse_vector_add();				\
  _sse_store(rs.s1);				\
  _sse_load(rs.s3);				\
  _sse_vector_i_mul();				\
  _sse_vector_add();				\
  _sse_store(rs.s3);

#define _hop_z_m()				\
  _prefetch_su3(up);				\
  _sse_load(sm->s0);				\
  _sse_load_up(sm->s2);				\
  _sse_vector_i_mul();				\
  _sse_vector_sub();				\
  _sse_su3_inverse_multiply((*um));		\
  _sse_vector_cmplxcg_mul(ka3);			\
  _sse_load(rs.s0);				\
  _sse_vector_add();				\
  _sse_store_nt(rn->s0);			\
  _sse_load(rs.s2);				\
  _sse_vector_i_mul();				\
  _sse_vector_add();				\
  _sse_store_nt(rn->s2);			\
  _sse_load(sm->s1);				\
  _sse_load_up(sm->s3);				\
  _sse_vector_i_mul();				\
  _sse_vector_add();				\
  _sse_su3_inverse_multiply((*um));		\
  _sse_vector_cmplxcg_mul(ka3);			\
  _sse_load(rs.s1);				\
  _sse_vector_add();				\
  _sse_store_nt(rn->s1);			\
  _sse_load(rs.s3);				\
  _sse_vector_i_mul();				\
  _sse_vector_sub();				\
  _sse_store_nt(rn->s3);

#define _hop_mul_g5_cmplx_and_store()			\
  _sse_load_up(rn->s0);					\
  _sse_vector_cmplx_mul(cf);				\
  _sse_store_nt_up(rn->s0);				\
  _sse_load_up(rn->s1);					\
  _sse_vector_cmplx_mul(cf);				\
  _sse_store_nt_up(rn->s1);				\
  _sse_load_up(rn->s2);					\
  _sse_vector_cmplxcg_mul(cf);				\
  _sse_store_nt_up(rn->s2);				\
  _sse_load_up(rn->s3);					\
  _sse_vector_cmplxcg_mul(cf);				\
  _sse_store_nt_up(rn->s3);

#define _g5_cmplx_sub_hop_and_g5store()			\
  _sse_load_up(pn->s0);					\
  _sse_vector_cmplx_mul(cf);				\
  _sse_load(rn->s0);					\
  _sse_vector_sub_up();					\
  _sse_store_nt_up(rn->s0);				\
  _sse_load_up(pn->s1);					\
  _sse_vector_cmplx_mul(cf);				\
  _sse_load(rn->s1);					\
  _sse_vector_sub_up();					\
  _sse_store_nt_up(rn->s1);				\
  _sse_load_up(pn->s2);					\
  _sse_vector_cmplxcg_mul(cf);				\
  _sse_load(rn->s2);					\
  _sse_vector_sub();					\
  _sse_store_nt(rn->s2);				\
  _sse_load_up(pn->s3);					\
  _sse_vector_cmplxcg_mul(cf);				\
  _sse_load(rn->s3);					\
  _sse_vector_sub();					\
  _sse_store_nt(rn->s3);
  
#define _store_res()

#  else

#define _declare_regs()	     \
  su3_vector ALIGN psi, chi; \
  spinor ALIGN temp;

#define _hop_t_p()				\
  _vector_add(psi,sp->s0,sp->s2);		\
  _su3_multiply(chi,(*up),psi);			\
  _complex_times_vector(psi,ka0,chi);		\
  _vector_assign(temp.s0,psi);			\
  _vector_assign(temp.s2,psi);			\
  _vector_add(psi,sp->s1,sp->s3);		\
  _su3_multiply(chi,(*up),psi);			\
  _complex_times_vector(psi,ka0,chi);		\
  _vector_assign(temp.s1,psi);			\
  _vector_assign(temp.s3,psi);
  
#define _hop_t_m()				\
  _vector_sub(psi,sm->s0,sm->s2);		\
  _su3_inverse_multiply(chi,(*um),psi);		\
  _complexcjg_times_vector(psi,ka0,chi);	\
  _vector_add_assign(temp.s0,psi);		\
  _vector_sub_assign(temp.s2,psi);		\
  _vector_sub(psi,sm->s1,sm->s3);		\
  _su3_inverse_multiply(chi,(*um),psi);		\
  _complexcjg_times_vector(psi,ka0,chi);	\
  _vector_add_assign(temp.s1,psi);		\
  _vector_sub_assign(temp.s3,psi);

#define _hop_x_p()				\
  _vector_i_add(psi,sp->s0,sp->s3);		\
  _su3_multiply(chi,(*up),psi);			\
  _complex_times_vector(psi,ka1,chi);		\
  _vector_add_assign(temp.s0,psi);		\
  _vector_i_sub_assign(temp.s3,psi);		\
  _vector_i_add(psi,sp->s1,sp->s2);		\
  _su3_multiply(chi,(*up),psi);			\
  _complex_times_vector(psi,ka1,chi);		\
  _vector_add_assign(temp.s1,psi);		\
  _vector_i_sub_assign(temp.s2,psi);

#define _hop_x_m()				\
  _vector_i_sub(psi,sm->s0,sm->s3);		\
  _su3_inverse_multiply(chi,(*um),psi);		\
  _complexcjg_times_vector(psi,ka1,chi);	\
  _vector_add_assign(temp.s0,psi);		\
  _vector_i_add_assign(temp.s3,psi);		\
  _vector_i_sub(psi,sm->s1,sm->s2);		\
  _su3_inverse_multiply(chi,(*um),psi);		\
  _complexcjg_times_vector(psi,ka1,chi);	\
  _vector_add_assign(temp.s1,psi);		\
  _vector_i_add_assign(temp.s2,psi);

#define _hop_y_p()				\
  _vector_add(psi,sp->s0,sp->s3);		\
  _su3_multiply(chi,(*up),psi);			\
  _complex_times_vector(psi,ka2,chi);		\
  _vector_add_assign(temp.s0,psi);		\
  _vector_add_assign(temp.s3,psi);		\
  _vector_sub(psi,sp->s1,sp->s2);		\
  _su3_multiply(chi,(*up),psi);			\
  _complex_times_vector(psi,ka2,chi);		\
  _vector_add_assign(temp.s1,psi);		\
  _vector_sub_assign(temp.s2,psi);

#define _hop_y_m()				\
  _vector_sub(psi,sm->s0,sm->s3);		\
  _su3_inverse_multiply(chi,(*um),psi);		\
  _complexcjg_times_vector(psi,ka2,chi);	\
  _vector_add_assign(temp.s0,psi);		\
  _vector_sub_assign(temp.s3,psi);		\
  _vector_add(psi,sm->s1,sm->s2);		\
  _su3_inverse_multiply(chi,(*um),psi);		\
  _complexcjg_times_vector(psi,ka2,chi);	\
  _vector_add_assign(temp.s1,psi);		\
  _vector_add_assign(temp.s2,psi);

#define _hop_z_p()				\
  _vector_i_add(psi,sp->s0,sp->s2);		\
  _su3_multiply(chi,(*up),psi);			\
  _complex_times_vector(psi,ka3,chi);		\
  _vector_add_assign(temp.s0,psi);		\
  _vector_i_sub_assign(temp.s2,psi);		\
  _vector_i_sub(psi,sp->s1,sp->s3);		\
  _su3_multiply(chi,(*up),psi);			\
  _complex_times_vector(psi,ka3,chi);		\
  _vector_add_assign(temp.s1,psi);		\
  _vector_i_add_assign(temp.s3,psi);

#define _hop_z_m()				\
  _vector_i_sub(psi,sm->s0,sm->s2);		\
  _su3_inverse_multiply(chi,(*um),psi);		\
  _complexcjg_times_vector(psi,ka3,chi);	\
  _vector_add_assign(temp.s0, psi);		\
  _vector_i_add_assign(temp.s2, psi);		\
  _vector_i_add(psi,sm->s1,sm->s3);		\
  _su3_inverse_multiply(chi,(*um),psi);		\
  _complexcjg_times_vector(psi,ka3,chi);	\
  _vector_add_assign(temp.s1, psi);		\
  _vector_i_sub_assign(temp.s3, psi);

#define _hop_mul_g5_cmplx_and_store()			\
  _complex_times_vector(rn->s0, cfactor, temp.s0);	\
  _complex_times_vector(rn->s1, cfactor, temp.s1);	\
  _complexcjg_times_vector(rn->s2, cfactor, temp.s2);	\
  _complexcjg_times_vector(rn->s3, cfactor, temp.s3);

#define _g5_cmplx_sub_hop_and_g5store()			\
  _complex_times_vector(psi, cfactor, pn->s0);		\
  _vector_sub(rn->s0, psi, temp.s0);			\
  _complex_times_vector(chi, cfactor, pn->s1);		\
  _vector_sub(rn->s1, chi, temp.s1);			\
  _complexcjg_times_vector(psi, cfactor, pn->s2);	\
  _vector_sub(rn->s2, temp.s2, psi);			\
  _complexcjg_times_vector(chi, cfactor, pn->s3);	\
  _vector_sub(rn->s3, temp.s3, chi);

#define _store_res()				\
  _vector_assign(rn->s0, temp.s0);		\
  _vector_assign(rn->s1, temp.s1);		\
  _vector_assign(rn->s2, temp.s2);		\
  _vector_assign(rn->s3, temp.s3);

#  endif

#endif

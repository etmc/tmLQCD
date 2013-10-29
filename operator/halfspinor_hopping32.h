/**********************************************************************
 *
 * Copyright (C) 2013  Florian Burger
 *
 * A 32-bit version of the Half-spinor implementation by Carsten Urbach
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

#ifndef _HALFSPINOR_HOPPING32_H
#define _HALFSPINOR_HOPPING32_H

#if (defined BGQ && defined XLC)

#define o_hop_t_p_pre32()					\
  _vec_load2_32(rs0, rs1, rs2, s->s0);				\
  _vec_load2_32(rs3, rs4, rs5, s->s1);				\
  _vec_load2_32(rs6, rs7, rs8, s->s2);				\
  _vec_load2_32(rs9, rs10, rs11, s->s3);				\
  _prefetch_spinor_32(s+1);					\
  _prefetch_su3_32(U+1);						\
  _vec_add_to2(r0, r1, r2, rs0, rs1, rs2, rs6, rs7, rs8);	\
  _vec_add_to2(r3, r4, r5, rs3, rs4, rs5, rs9, rs10, rs11);	\
  rtmp = vec_ld2a(0L, (float*) &ka0_32);					\
  _vec_su3_multiply_double2_32(U);						\
  _vec_cmplx_mul_double2(r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, rtmp); \
  _vec_store2_32(phi2[ix]->s0, r0, r1, r2);				\
  _vec_store2_32(phi2[ix]->s1, r3, r4, r5);

#define o_hop_t_m_pre32()					\
  _vec_sub_to2(r0, r1, r2, rs0, rs1, rs2, rs6, rs7, rs8);	\
  _vec_sub_to2(r3, r4, r5, rs3, rs4, rs5, rs9, rs10, rs11);	\
  _vec_store2_32(phi2[ix]->s0, r0, r1, r2);			\
  _vec_store2_32(phi2[ix]->s1, r3, r4, r5);

#define o_hop_x_p_pre32()						\
  _prefetch_su3_32(U+1);							\
  _vec_i_mul_add_to2(r0, r1, r2, rs0, rs1, rs2, rs9, rs10, rs11, U0);	\
  _vec_i_mul_add_to2(r3, r4, r5, rs3, rs4, rs5, rs6, rs7, rs8, U0);	\
  rtmp = vec_ld2a(0L, (float*) &ka1_32);					\
  _vec_su3_multiply_double2_32(U);						\
  _vec_cmplx_mul_double2(r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, rtmp); \
  _vec_store2_32(phi2[ix]->s0, r0, r1, r2);				\
  _vec_store2_32(phi2[ix]->s1, r3, r4, r5);

#define o_hop_x_m_pre32()						\
  _vec_i_mul_sub_to2(r0, r1, r2, rs0, rs1, rs2, rs9, rs10, rs11, U0);	\
  _vec_i_mul_sub_to2(r3, r4, r5, rs3, rs4, rs5, rs6, rs7, rs8, U0);	\
  _vec_store2_32(phi2[ix]->s0, r0, r1, r2);				\
  _vec_store2_32(phi2[ix]->s1, r3, r4, r5);

#define o_hop_y_p_pre32()					\
  _prefetch_su3_32(U+1);						\
  _vec_add_to2(r0, r1, r2, rs0, rs1, rs2, rs9, rs10, rs11);	\
  _vec_sub_to2(r3, r4, r5, rs3, rs4, rs5, rs6, rs7, rs8);	\
  rtmp = vec_ld2a(0L, (float*) &ka2_32);					\
  _vec_su3_multiply_double2_32(U);						\
  _vec_cmplx_mul_double2(r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, rtmp); \
  _vec_store2_32(phi2[ix]->s0, r0, r1, r2);				\
  _vec_store2_32(phi2[ix]->s1, r3, r4, r5);

#define o_hop_y_m_pre32()					\
  _vec_sub_to2(r0, r1, r2, rs0, rs1, rs2, rs9, rs10, rs11);	\
  _vec_add_to2(r3, r4, r5, rs3, rs4, rs5, rs6, rs7, rs8);	\
  _vec_store2_32(phi2[ix]->s0, r0, r1, r2);			\
  _vec_store2_32(phi2[ix]->s1, r3, r4, r5);

#define o_hop_z_p_pre32()						\
  _prefetch_su3_32(U+1);							\
  _vec_i_mul_add_to2(r0, r1, r2, rs0, rs1, rs2, rs6, rs7, rs8, U0);	\
  _vec_i_mul_sub_to2(r3, r4, r5, rs3, rs4, rs5, rs9, rs10, rs11, U0);	\
  rtmp = vec_ld2a(0L, (float*) &ka3_32);					\
  _vec_su3_multiply_double2_32(U);						\
  _vec_cmplx_mul_double2(r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, rtmp); \
  _vec_store2_32(phi2[ix]->s0, r0, r1, r2);				\
  _vec_store2_32(phi2[ix]->s1, r3, r4, r5);

#define o_hop_z_m_pre32()						\
  _vec_i_mul_sub_to2(r0, r1, r2, rs0, rs1, rs2, rs6, rs7, rs8, U0);	\
  _vec_i_mul_add_to2(r3, r4, r5, rs3, rs4, rs5, rs9, rs10, rs11, U0);	\
  _vec_store2_32(phi2[ix]->s0, r0, r1, r2);				\
  _vec_store2_32(phi2[ix]->s1, r3, r4, r5);

#define o_hop_t_p_post32()				\
  _vec_load2_32(rs0, rs1, rs2, phi2[ix]->s0);		\
  rs6 = rs0;						\
  rs7 = rs1;						\
  rs8 = rs2;						\
  _vec_load2_32(rs3, rs4, rs5, phi2[ix]->s1);		\
  rs9 = rs3;						\
  rs10= rs4;						\
  rs11= rs5;

#define o_hop_t_m_post32()						\
  _prefetch_su3_32(U+1);							\
  _vec_load2_32(r0, r1, r2, phi2[ix]->s0);				\
  _vec_load2_32(r3, r4, r5, phi2[ix]->s1);				\
  rtmp = vec_ld2a(0L, (float*) &ka0_32);					\
  _vec_su3_inverse_multiply_double2_32(U);					\
  _vec_cmplxcg_mul_double2(r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, rtmp); \
  _vec_add2(rs0, rs1, rs2, r0, r1, r2);					\
  _vec_sub2(rs6, rs7, rs8, r0, r1, r2);					\
  _vec_add2(rs3, rs4, rs5, r3, r4, r5);					\
  _vec_sub2(rs9, rs10, rs11, r3, r4, r5);

#define o_hop_x_p_post32()						\
  _vec_load2_32(r0, r1, r2, phi2[ix]->s0);				\
  _vec_load2_32(r3, r4, r5, phi2[ix]->s1);				\
  _vec_add2(rs0, rs1, rs2, r0, r1, r2);					\
  _vec_i_mul_sub2(rs9, rs10, rs11, r0, r1, r2, U0);			\
  _vec_add2(rs3, rs4, rs5, r3, r4, r5);					\
  _vec_i_mul_sub2(rs6, rs7, rs8, r3, r4, r5, U0);

#define o_hop_x_m_post32()						\
  _prefetch_su3_32(U+1);							\
  _vec_load2_32(r0, r1, r2, phi2[ix]->s0);				\
  _vec_load2_32(r3, r4, r5, phi2[ix]->s1);				\
  rtmp = vec_ld2a(0L, (float*) &ka1_32);					\
  _vec_su3_inverse_multiply_double2_32(U);					\
  _vec_cmplxcg_mul_double2(r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, rtmp); \
  _vec_add_double2(rs0, rs1, rs2, rs3, rs4, rs5, r0, r1, r2, r3, r4, r5); \
  _vec_i_mul_add_double2(rs9, rs10, rs11, rs6, rs7, rs8, r0, r1, r2, r3, r4, r5, U0);

#define o_hop_y_p_post32()						\
  _vec_load2_32(r0, r1, r2, phi2[ix]->s0);				\
  _vec_load2_32(r3, r4, r5, phi2[ix]->s1);				\
  _vec_add_double2(rs0, rs1, rs2, rs3, rs4, rs5, r0, r1, r2, r3, r4, r5); \
  _vec_sub2(rs6, rs7, rs8, r3, r4, r5);					\
  _vec_add2(rs9, rs10, rs11, r0, r1, r2);

#define o_hop_y_m_post32()						\
  _prefetch_su3_32(U+1);							\
  _vec_load2_32(r0, r1, r2, phi2[ix]->s0);				\
  _vec_load2_32(r3, r4, r5, phi2[ix]->s1);				\
  rtmp = vec_ld2a(0L, (float*) &ka2_32);					\
  _vec_su3_inverse_multiply_double2_32(U);					\
  _vec_cmplxcg_mul_double2(r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, rtmp); \
  _vec_add_double2(rs0, rs1, rs2, rs3, rs4, rs5, r0, r1, r2, r3, r4, r5); \
  _vec_add2(rs6, rs7, rs8, r3, r4, r5);					\
  _vec_sub2(rs9, rs10, rs11, r0, r1, r2);

#define o_hop_z_p_post32()						\
  _vec_load2_32(r0, r1, r2, phi2[ix]->s0);				\
  _vec_load2_32(r3, r4, r5, phi2[ix]->s1);				\
  _vec_add_double2(rs0, rs1, rs2, rs3, rs4, rs5, r0, r1, r2, r3, r4, r5); \
  _vec_i_mul_sub2(rs6, rs7, rs8, r0, r1, r2, U0);			\
  _vec_i_mul_add2(rs9, rs10, rs11, r3, r4, r5, U0);

#define o_hop_z_m_post32()						\
  _prefetch_su3_32(U+1);							\
  _vec_load2_32(r0, r1, r2, phi2[ix]->s0);				\
  _vec_load2_32(r3, r4, r5, phi2[ix]->s1);				\
  rtmp = vec_ld2a(0L, (float*) &ka3_32);					\
  _vec_su3_inverse_multiply_double2_32(U);					\
  _vec_cmplxcg_mul_double2(r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, rtmp); \
  _vec_add2(rs0, rs1, rs2, r0, r1, r2);					\
  _vec_i_mul_add2(rs6, rs7, rs8, r0, r1, r2, U0);			\
  _vec_add2(rs3, rs4, rs5, r3, r4, r5);					\
  _vec_i_mul_sub2(rs9, rs10, rs11, r3, r4, r5, U0);

  
  
  
  
  
  
  
  
  
  
//new versions





#define _hop_t_p_pre32()							\
  _vec_load_32(rs0, rs1, s->s0);						\
  _vec_load16_32(rs2, rs3, s->s1, rtmp);					\
  _vec_load_32(rs4, rs5, s->s2);						\
  _vec_load16_32(rs6, rs7, s->s3, rtmp);					\
  _prefetch_spinor_32(s+1);						\
  _prefetch_su3_32(U+1);							\
  _vec_add(r0, r1, rs0, rs1, rs4, rs5);					\
  _vec_add(r2, r3, rs2, rs3, rs6, rs7);					\
  _vec_su3_multiply_double2c_32(U);					\
  rtmp = vec_ld2(0, (float*) &ka0_32);					\
  _vec_cmplx_mul_double2c(r0, r1, r2, r4, r5, r6, rtmp);		\
  _vec_store_halfspinor_32(phi2[ix]->s0, r0, r1, r2);



#define _hop_t_m_pre32()						\
  _vec_sub(r0, r1, rs0, rs1, rs4, rs5);				\
  _vec_sub(r2, r3, rs2, rs3, rs6, rs7);				\
  _vec_store_32(phi2[ix]->s0, r0, r1);				\
  _vec_store16_32(phi2[ix]->s1, r2, r3, U0);


#define _hop_x_p_pre32()						\
  _prefetch_su3_32(U+1);						\
  _vec_i_mul_add(r0, r1, rs0, rs1, rs6, rs7, U0);		\
  _vec_i_mul_add(r2, r3, rs2, rs3, rs4, rs5, U0);		\
  rtmp = vec_ld2(0, (float*) &ka1_32);				\
  _vec_su3_multiply_double2c_32(U);				\
  _vec_cmplx_mul_double2c(r0, r1, r2, r4, r5, r6, rtmp);	\
  _vec_store_halfspinor_32(phi2[ix]->s0, r0, r1, r2);


#define _hop_x_m_pre32()					\
  _vec_i_mul_sub(r0, r1, rs0, rs1, rs6, rs7, U0);	\
  _vec_i_mul_sub(r2, r3, rs2, rs3, rs4, rs5, U0);	\
  _vec_store_32(phi2[ix]->s0, r0, r1);			\
  _vec_store16_32(phi2[ix]->s1, r2, r3, U0);


#define _hop_y_p_pre32()						\
  _prefetch_su3_32(U+1);						\
  _vec_add(r0, r1, rs0, rs1, rs6, rs7);				\
  _vec_sub(r2, r3, rs2, rs3, rs4, rs5);				\
  rtmp = vec_ld2(0, (float*) &ka2_32);				\
  _vec_su3_multiply_double2c_32(U);				\
  _vec_cmplx_mul_double2c(r0, r1, r2, r4, r5, r6, rtmp);	\
  _vec_store_halfspinor_32(phi2[ix]->s0, r0, r1, r2);



#define _hop_y_m_pre32()				\
  _vec_sub(r0, r1, rs0, rs1, rs6, rs7);		\
  _vec_add(r2, r3, rs2, rs3, rs4, rs5);		\
  _vec_store_32(phi2[ix]->s0, r0, r1);		\
  _vec_store16_32(phi2[ix]->s1, r2, r3, U0);

  
#define _hop_z_p_pre32()						\
  _prefetch_su3_32(U+1);						\
  _vec_i_mul_add(r0, r1, rs0, rs1, rs4, rs5, U0);		\
  _vec_i_mul_sub(r2, r3, rs2, rs3, rs6, rs7, U0);		\
  rtmp = vec_ld2(0, (float*) &ka3_32);				\
  _vec_su3_multiply_double2c_32(U);				\
  _vec_cmplx_mul_double2c(r0, r1, r2, r4, r5, r6, rtmp);	\
  _vec_store_halfspinor_32(phi2[ix]->s0, r0, r1, r2);


#define _hop_z_m_pre32()					\
  _vec_i_mul_sub(r0, r1, rs0, rs1, rs4, rs5, U0);	\
  _vec_i_mul_add(r2, r3, rs2, rs3, rs6, rs7, U0);	\
  _vec_store_32(phi2[ix]->s0, r0, r1);			\
  _vec_store16_32(phi2[ix]->s1, r2, r3, U0);

  
#define _hop_t_p_post32()				\
  _vec_load_halfspinor_32(rs0, rs1, rs2, phi2[ix]->s0);	\
  _vec_unfuse(rs0, rs1, rs2, rs3, rs4, rs5);		\
  rs6 = rs0; rs7 = rs1; rs8 = rs2;			\
  rs9 = rs3; rs10= rs4; rs11= rs5;


#define _hop_t_m_post32()						\
  _prefetch_su3_32(U+1);							\
  _vec_load_32(r0, r1, phi2[ix]->s0);					\
  _vec_load16_32(r2, r3, phi2[ix]->s1, rtmp);				\
  rtmp = vec_ld2(0, (float*) &ka0_32);					\
  _vec_su3_inverse_multiply_double2c_32(U);				\
  _vec_cmplxcg_mul_double2c(r0, r1, r2, r4, r5, r6, rtmp);		\
  _vec_unfuse(r0, r1, r2, r3, r4, r5);					\
  _vec_add_double2(rs0, rs1, rs2, rs3, rs4, rs5, r0, r1, r2, r3, r4, r5); \
  _vec_sub_double2(rs6, rs7, rs8, rs9, rs10, rs11, r0, r1, r2, r3, r4, r5);


#define _hop_x_p_post32()						\
  _vec_load_halfspinor_32(r0, r1, r2, phi2[ix]->s0);			\
  _vec_unfuse(r0, r1, r2, r3, r4, r5);				\
  _vec_add_double2(rs0, rs1, rs2, rs3, rs4, rs5, r0, r1, r2, r3, r4, r5); \
  _vec_i_mul_sub2(rs6, rs7, rs8, r3, r4, r5, U0);			\
  _vec_i_mul_sub2(rs9, rs10, rs11, r0, r1, r2, U1);
  

#define _hop_x_m_post32()						\
  _prefetch_su3_32(U+1);							\
  _vec_load_32(r0, r1, phi2[ix]->s0);					\
  _vec_load16_32(r2, r3, phi2[ix]->s1, rtmp);				\
  rtmp = vec_ld2(0, (float*) &ka1_32);					\
  _vec_su3_inverse_multiply_double2c_32(U);				\
  _vec_cmplxcg_mul_double2c(r0, r1, r2, r4, r5, r6, rtmp);		\
  _vec_unfuse(r0, r1, r2, r3, r4, r5);					\
  _vec_add_double2(rs0, rs1, rs2, rs3, rs4, rs5, r0, r1, r2, r3, r4, r5); \
  _vec_i_mul_add_double2(rs9, rs10, rs11, rs6, rs7, rs8, r0, r1, r2, r3, r4, r5, U0);



#define _hop_y_p_post32()						\
  _vec_load_halfspinor_32(r0, r1, r2, phi2[ix]->s0);			\
  _vec_unfuse(r0, r1, r2, r3, r4, r5);				\
  _vec_add_double2(rs0, rs1, rs2, rs3, rs4, rs5, r0, r1, r2, r3, r4, r5);	\
  _vec_sub2(rs6, rs7, rs8, r3, r4, r5);					\
  _vec_add2(rs9, rs10, rs11, r0, r1, r2);


#define _hop_y_m_post32()						\
  _prefetch_su3_32(U+1);							\
  _vec_load_32(r0, r1, phi2[ix]->s0);					\
  _vec_load16_32(r2, r3, phi2[ix]->s1, rtmp);				\
  rtmp = vec_ld2(0, (float*) &ka2_32);					\
  _vec_su3_inverse_multiply_double2c_32(U);				\
  _vec_cmplxcg_mul_double2c(r0, r1, r2, r4, r5, r6, rtmp);		\
  _vec_unfuse(r0, r1, r2, r3, r4, r5);					\
  _vec_add_double2(rs0, rs1, rs2, rs3, rs4, rs5, r0, r1, r2, r3, r4, r5); \
  _vec_add2(rs6, rs7, rs8, r3, r4, r5);					\
  _vec_sub2(rs9, rs10, rs11, r0, r1, r2);


#define _hop_z_p_post32()						\
  _vec_load_halfspinor_32(r0, r1, r2, phi2[ix]->s0);			\
  _vec_unfuse(r0, r1, r2, r3, r4, r5);				\
  _vec_add_double2(rs0, rs1, rs2, rs3, rs4, rs5, r0, r1, r2, r3, r4, r5); \
  _vec_i_mul_sub2(rs6, rs7, rs8, r0, r1, r2, U0);			\
  _vec_i_mul_add2(rs9, rs10, rs11, r3, r4, r5, U1);


#define _hop_z_m_post32()						\
  _prefetch_su3_32(U+1);							\
  _vec_load_32(r0, r1, phi2[ix]->s0);					\
  _vec_load16_32(r2, r3, phi2[ix]->s1, rtmp);				\
  rtmp = vec_ld2(0, (float*) &ka3_32);					\
  _vec_su3_inverse_multiply_double2c_32(U);				\
  _vec_cmplxcg_mul_double2c(r0, r1, r2, r4, r5, r6, rtmp);		\
  _vec_unfuse(r0, r1, r2, r3, r4, r5);					\
  _vec_add_double2(rs0, rs1, rs2, rs3, rs4, rs5, r0, r1, r2, r3, r4, r5); \
  _vec_i_mul_add2(rs6, rs7, rs8, r0, r1, r2, U0);			\
  _vec_i_mul_sub2(rs9, rs10, rs11, r3, r4, r5, U1);

  
  
  
//end new versions
  
  
  
  
  
#define _hop_mul_g5_cmplx_and_store32(res)					\
  _vec_cmplx_mul_double2(r0, r1, r2, r3, r4, r5, rs0, rs1, rs2, rs3, rs4, rs5, cf); \
  _vec_cmplxcg_mul_double2(r6, r7, r8, r9, r10, r11, rs6, rs7, rs8, rs9, rs10, rs11, cf); \
  _vec_store2_32((res)->s0, r0, r1, r2);					\
  _vec_store2_32((res)->s1, r3, r4, r5);					\
  _vec_store2_32((res)->s2, r6, r7, r8);					\
  _vec_store2_32((res)->s3, r9, r10, r11);

#define _g5_cmplx_sub_hop_and_g5store32(res)					\
  _vec_load_halfspinor_32(r3, r4, r5, pn->s0);				\
  _vec_cmplx_mul_double2c_32(r0, r1, r2, r3, r4, r5, cf);			\
  _vec_unfuse(r0, r1, r2, r3, r4, r5);					\
  _vec_sub_double2(r0, r3, r1, r4, r2, r5, rs0, rs1, rs2, rs3, rs4, rs5); \
  _vec_store2_32((res)->s0, r0, r3, r1);					\
  _vec_store2_32((res)->s1, r4, r2, r5);					\
  _vec_load_halfspinor_32(r3, r4, r5, pn->s2);				\
  _vec_cmplxcg_mul_double2c(r0, r1, r2, r3, r4, r5, cf);		\
  _vec_unfuse(r0, r1, r2, r3, r4, r5);					\
  _vec_sub_double2(rs6, rs7, rs8, rs9, rs10, rs11, r0, r3, r1, r4, r2, r5); \
  _vec_store2_32((res)->s2, rs6, rs7, rs8);					\
  _vec_store2_32((res)->s3, rs9, rs10, rs11);

#define _hop_store_post32(res)		\
  _vec_store2_32((res)->s0, rs0, rs1, rs2);	\
  _vec_store2_32((res)->s1, rs3, rs4, rs5);	\
  _vec_store2_32((res)->s2, rs6, rs7, rs8);	\
  _vec_store2_32((res)->s3, rs9, rs10, rs11);


#define _declare_hregs()						\
  vector4double ALIGN r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11;	\
  vector4double ALIGN rs0, rs1, rs2, rs3, rs4, rs5, rs6, rs7, rs8, rs9, rs10, rs11; \
  vector4double ALIGN U0, U1, U2, U3, U4, U6, U7;			\
  vector4double ALIGN rtmp;

#else

#define _prefetch_spinor(s)
#define _prefetch_halfspinor(hs)
#define _prefetch_spinor_32(s)
#define _prefetch_su3_32(U)


#define _hop_t_p_pre32()				\
  _vector_assign(rs.s0, s->s0);				\
  _vector_assign(rs.s1, s->s1);				\
  _vector_assign(rs.s2, s->s2);				\
  _vector_assign(rs.s3, s->s3);				\
  _vector_add(psi, rs.s0, rs.s2);			\
  _su3_multiply(chi,(*U),psi);				\
  _complex_times_vector(phi2[ix]->s0, ka0_32, chi);	\
  _vector_add(psi, rs.s1, rs.s3);			\
  _su3_multiply(chi,(*U),psi);				\
  _complex_times_vector(phi2[ix]->s1, ka0_32, chi);

#define _hop_t_m_pre32()				\
  _vector_sub(phi2[ix]->s0, rs.s0, rs.s2);		\
  _vector_sub(phi2[ix]->s1, rs.s1, rs.s3);

#define _hop_x_p_pre32()				\
  _vector_i_add(psi, rs.s0, rs.s3);			\
  _su3_multiply(chi, (*U), psi);			\
  _complex_times_vector(phi2[ix]->s0, ka1_32, chi);	\
  _vector_i_add(psi, rs.s1, rs.s2);			\
  _su3_multiply(chi, (*U), psi);			\
  _complex_times_vector(phi2[ix]->s1, ka1_32, chi);

#define _hop_x_m_pre32()				\
  _vector_i_sub(phi2[ix]->s0, rs.s0, rs.s3);		\
  _vector_i_sub(phi2[ix]->s1, rs.s1, rs.s2);

#define _hop_y_p_pre32()				\
  _vector_add(psi, rs.s0, rs.s3);			\
  _su3_multiply(chi,(*U),psi);				\
  _complex_times_vector(phi2[ix]->s0, ka2_32, chi);	\
  _vector_sub(psi, rs.s1, rs.s2);			\
  _su3_multiply(chi,(*U),psi);				\
  _complex_times_vector(phi2[ix]->s1, ka2_32, chi);

#define _hop_y_m_pre32()			\
  _vector_sub(phi2[ix]->s0, rs.s0, rs.s3);	\
  _vector_add(phi2[ix]->s1, rs.s1, rs.s2);

#define _hop_z_p_pre32()				\
  _vector_i_add(psi, rs.s0, rs.s2);			\
  _su3_multiply(chi, (*U), psi);			\
  _complex_times_vector(phi2[ix]->s0, ka3_32, chi);	\
  _vector_i_sub(psi, rs.s1, rs.s3);			\
  _su3_multiply(chi,(*U),psi);				\
  _complex_times_vector(phi2[ix]->s1, ka3_32, chi);

#define _hop_z_m_pre32()			\
  _vector_i_sub(phi2[ix]->s0, rs.s0, rs.s2);	\
  _vector_i_add(phi2[ix]->s1, rs.s1, rs.s3);

#define _hop_t_p_post32();			\
  _vector_assign(rs.s0, phi2[ix]->s0);		\
  _vector_assign(rs.s2, phi2[ix]->s0);		\
  _vector_assign(rs.s1, phi2[ix]->s1);		\
  _vector_assign(rs.s3, phi2[ix]->s1);		\

#define _hop_t_m_post32();			\
  _vector_assign(psi, phi2[ix]->s0);		\
  _su3_inverse_multiply(chi,(*U), psi);		\
  _complexcjg_times_vector(psi,ka0_32,chi);	\
  _vector_add_assign(rs.s0, psi);		\
  _vector_sub_assign(rs.s2, psi);		\
  _vector_assign(psi, phi2[ix]->s1);		\
  _su3_inverse_multiply(chi,(*U), psi);		\
  _complexcjg_times_vector(psi,ka0_32,chi);	\
  _vector_add_assign(rs.s1, psi);		\
  _vector_sub_assign(rs.s3, psi);

#define _hop_x_p_post32();				\
  _vector_add_assign(rs.s0, phi2[ix]->s0);		\
  _vector_i_sub_assign(rs.s3, phi2[ix]->s0);		\
  _vector_add_assign(rs.s1, phi2[ix]->s1);		\
  _vector_i_sub_assign(rs.s2, phi2[ix]->s1);

#define _hop_x_m_post32();			\
  _vector_assign(psi, phi2[ix]->s0);		\
  _su3_inverse_multiply(chi,(*U), psi);		\
  _complexcjg_times_vector(psi,ka1_32,chi);	\
  _vector_add_assign(rs.s0, psi);		\
  _vector_i_add_assign(rs.s3, psi);		\
  _vector_assign(psi, phi2[ix]->s1);		\
  _su3_inverse_multiply(chi,(*U), psi);		\
  _complexcjg_times_vector(psi,ka1_32,chi);	\
  _vector_add_assign(rs.s1, psi);		\
  _vector_i_add_assign(rs.s2, psi);

#define _hop_y_p_post32();			\
  _vector_add_assign(rs.s0, phi2[ix]->s0);	\
  _vector_add_assign(rs.s3, phi2[ix]->s0);	\
  _vector_add_assign(rs.s1, phi2[ix]->s1);	\
  _vector_sub_assign(rs.s2, phi2[ix]->s1);

#define _hop_y_m_post32();			\
  _vector_assign(psi, phi2[ix]->s0);		\
  _su3_inverse_multiply(chi,(*U), psi);		\
  _complexcjg_times_vector(psi,ka2_32,chi);	\
  _vector_add_assign(rs.s0, psi);		\
  _vector_sub_assign(rs.s3, psi);		\
  _vector_assign(psi, phi2[ix]->s1);		\
  _su3_inverse_multiply(chi, (*U), psi);	\
  _complexcjg_times_vector(psi,ka2_32,chi);	\
  _vector_add_assign(rs.s1, psi);		\
  _vector_add_assign(rs.s2, psi);

#define _hop_z_p_post32();			\
  _vector_add_assign(rs.s0, phi2[ix]->s0);	\
  _vector_i_sub_assign(rs.s2, phi2[ix]->s0);	\
  _vector_add_assign(rs.s1, phi2[ix]->s1);	\
  _vector_i_add_assign(rs.s3, phi2[ix]->s1);

#define _hop_z_m_post32();			\
  _vector_assign(psi, phi2[ix]->s0);		\
  _su3_inverse_multiply(chi,(*U), psi);		\
  _complexcjg_times_vector(psi,ka3_32,chi);	\
  _vector_add_assign(rs.s0, psi);		\
  _vector_i_add_assign(rs.s2, psi);		\
  _vector_assign(psi, phi2[ix]->s1);		\
  _su3_inverse_multiply(chi,(*U), psi);		\
  _complexcjg_times_vector(psi,ka3_32,chi);	\
  _vector_add_assign(rs.s1, psi);		\
  _vector_i_sub_assign(rs.s3, psi);

#define _hop_mul_g5_cmplx_and_store32(res)			\
  _complex_times_vector((res)->s0, cfactor, rs.s0);		\
  _complex_times_vector((res)->s1, cfactor, rs.s1);		\
  _complexcjg_times_vector((res)->s2, cfactor, rs.s2);	\
  _complexcjg_times_vector((res)->s3, cfactor, rs.s3);

#define _g5_cmplx_sub_hop_and_g5store32(res)		\
  _complex_times_vector(psi, cfactor, pn->s0);		\
  _vector_sub((res)->s0, psi, rs.s0);			\
  _complex_times_vector(psi2, cfactor, pn->s1);		\
  _vector_sub((res)->s1, psi2, rs.s1);			\
  _complexcjg_times_vector(psi, cfactor, pn->s2);	\
  _vector_sub((res)->s2, rs.s2, psi);			\
  _complexcjg_times_vector(psi2, cfactor, pn->s3);	\
  _vector_sub((res)->s3, rs.s3, psi2);


#define _hop_store_post32(res)		\
  _vector_assign(res->s0, rs.s0);	\
  _vector_assign(res->s1, rs.s1);	\
  _vector_assign(res->s2, rs.s2);	\
  _vector_assign(res->s3, rs.s3);


#define _declare_hregs()				\
  spinor32 ALIGN32 rs;					\
  su3_vector32 ALIGN32 psi, chi, psi2, chi2;

#endif

#endif


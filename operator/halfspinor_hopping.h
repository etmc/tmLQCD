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

#ifndef _HALFSPINOR_HOPPING_H
#define _HALFSPINOR_HOPPING_H

#if (defined SSE2 || defined SSE3)

#define _hop_t_p_pre32()
#define _hop_t_m_pre32()
#define _hop_x_p_pre32()
#define _hop_x_m_pre32()
#define _hop_y_p_pre32()
#define _hop_y_m_pre32()
#define _hop_z_p_pre32()
#define _hop_z_m_pre32()
#define _hop_t_p_post32()
#define _hop_t_m_post32()
#define _hop_x_p_post32()
#define _hop_x_m_post32()
#define _hop_y_p_post32()
#define _hop_y_m_post32()
#define _hop_z_p_post32()
#define _hop_z_m_post32()

#define _hop_t_p_pre()					\
  _prefetch_su3(U+predist);				\
  _sse_load(s->s0);					\
  _sse_load_up(s->s2);					\
  _sse_vector_add();					\
  _sse_su3_multiply((*U));				\
  _sse_vector_cmplx_mul(ka0);				\
  _sse_store_nt_up(phi[ix]->s0);			\
  _sse_load(s->s1);					\
  _sse_load_up(s->s3);					\
  _sse_vector_add();					\
  _sse_su3_multiply((*U));				\
  _sse_vector_cmplx_mul(ka0);				\
  _sse_store_nt_up(phi[ix]->s1);

#define _hop_t_m_pre()				\
  _sse_load(s->s0);				\
  _sse_load_up(s->s2);				\
  _sse_vector_sub();				\
  _sse_store_nt(phi[ix]->s0);			\
  _sse_load(s->s1);				\
  _sse_load_up(s->s3);				\
  _sse_vector_sub();				\
  _sse_store_nt(phi[ix]->s1);

#define _hop_x_p_pre()					\
  _prefetch_su3(U+predist);				\
  _sse_load(s->s0);					\
  _sse_load_up(s->s3);					\
  _sse_vector_i_mul();					\
  _sse_vector_add();					\
  _sse_su3_multiply((*U));				\
  _sse_vector_cmplx_mul(ka1);				\
  _sse_store_nt_up(phi[ix]->s0);			\
  _sse_load(s->s1);					\
  _sse_load_up(s->s2);					\
  _sse_vector_i_mul();					\
  _sse_vector_add();					\
  _sse_su3_multiply((*U));				\
  _sse_vector_cmplx_mul(ka1);				\
  _sse_store_nt_up(phi[ix]->s1);

#define _hop_x_m_pre()				\
  _sse_load(s->s0);				\
  _sse_load_up(s->s3);				\
  _sse_vector_i_mul();				\
  _sse_vector_sub();				\
  _sse_store_nt(phi[ix]->s0);			\
  _sse_load(s->s1);				\
  _sse_load_up(s->s2);				\
  _sse_vector_i_mul();				\
  _sse_vector_sub();				\
  _sse_store_nt(phi[ix]->s1);

#define _hop_y_p_pre()					\
  _prefetch_su3(U+predist);				\
  _sse_load(s->s0);					\
  _sse_load_up(s->s3);					\
  _sse_vector_add();					\
  _sse_su3_multiply((*U));				\
  _sse_vector_cmplx_mul(ka2);				\
  _sse_store_nt_up(phi[ix]->s0);			\
  _sse_load(s->s1);					\
  _sse_load_up(s->s2);					\
  _sse_vector_sub();					\
  _sse_su3_multiply((*U));				\
  _sse_vector_cmplx_mul(ka2);				\
  _sse_store_nt_up(phi[ix]->s1);

#define _hop_y_m_pre()				\
  _sse_load(s->s0);				\
  _sse_load_up(s->s3);				\
  _sse_vector_sub();				\
  _sse_store_nt(phi[ix]->s0);			\
  _sse_load(s->s1);				\
  _sse_load_up(s->s2);				\
  _sse_vector_add();				\
  _sse_store_nt(phi[ix]->s1);

#define _hop_z_p_pre()					\
  _prefetch_su3(U+predist);				\
  _prefetch_spinor(s+1);				\
  _sse_load(s->s0);					\
  _sse_load_up(s->s2);					\
  _sse_vector_i_mul();					\
  _sse_vector_add();					\
  _sse_su3_multiply((*U));				\
  _sse_vector_cmplx_mul(ka3);				\
  _sse_store_nt_up(phi[ix]->s0);			\
  _sse_load(s->s1);					\
  _sse_load_up(s->s3);					\
  _sse_vector_i_mul();					\
  _sse_vector_sub();					\
  _sse_su3_multiply((*U));				\
  _sse_vector_cmplx_mul(ka3);				\
  _sse_store_nt_up(phi[ix]->s1);			\

#define _hop_z_m_pre()				\
  _sse_load(s->s0);				\
  _sse_load_up(s->s2);				\
  _sse_vector_i_mul();				\
  _sse_vector_sub();				\
  _sse_store_nt(phi[ix]->s0);			\
  _sse_load(s->s1);				\
  _sse_load_up(s->s3);				\
  _sse_vector_i_mul();				\
  _sse_vector_add();				\
  _sse_store_nt(phi[ix]->s1);

#define _hop_t_p_post()					\
  _vector_assign(rs.s0, phi[ix]->s0);			\
  _vector_assign(rs.s2, phi[ix]->s0);			\
  _vector_assign(rs.s1, phi[ix]->s1);			\
  _vector_assign(rs.s3, phi[ix]->s1);

#define _hop_t_m_post()				\
  _prefetch_su3(U+predist);			\
  _sse_load(phi[ix]->s0);			\
  _sse_su3_inverse_multiply((*U));		\
  _sse_vector_cmplxcg_mul(ka0);			\
  _sse_load(rs.s0);				\
  _sse_vector_add();				\
  _sse_store(rs.s0);				\
  _sse_load(rs.s2);				\
  _sse_vector_sub();				\
  _sse_store(rs.s2);				\
  _sse_load(phi[ix]->s1);			\
  _sse_su3_inverse_multiply((*U));		\
  _sse_vector_cmplxcg_mul(ka0);			\
  _sse_load(rs.s1);				\
  _sse_vector_add();				\
  _sse_store(rs.s1);				\
  _sse_load(rs.s3);				\
  _sse_vector_sub();				\
  _sse_store(rs.s3);

#define _hop_x_p_post()					\
  _sse_load_up(phi[ix]->s0);				\
  _sse_load(rs.s0);					\
  _sse_vector_add();					\
  _sse_store(rs.s0);					\
  _sse_load(rs.s3);					\
  _sse_vector_i_mul();					\
  _sse_vector_sub();					\
  _sse_store(rs.s3);					\
  _sse_load_up(phi[ix]->s1);				\
  _sse_load(rs.s1);					\
  _sse_vector_add();					\
  _sse_store(rs.s1);					\
  _sse_load(rs.s2);					\
  _sse_vector_i_mul();					\
  _sse_vector_sub();					\
  _sse_store(rs.s2);       

#define _hop_x_m_post()				\
  _prefetch_su3(U+predist);			\
  _sse_load(phi[ix]->s0);			\
  _sse_su3_inverse_multiply((*U));		\
  _sse_vector_cmplxcg_mul(ka1);			\
  _sse_load(rs.s0);				\
  _sse_vector_add();				\
  _sse_store(rs.s0);				\
  _sse_load(rs.s3);				\
  _sse_vector_i_mul();				\
  _sse_vector_add();				\
  _sse_store(rs.s3);				\
  _sse_load(phi[ix]->s1);			\
  _sse_su3_inverse_multiply((*U));		\
  _sse_vector_cmplxcg_mul(ka1);			\
  _sse_load(rs.s1);				\
  _sse_vector_add();				\
  _sse_store(rs.s1);				\
  _sse_load(rs.s2);				\
  _sse_vector_i_mul();				\
  _sse_vector_add();				\
  _sse_store(rs.s2);

#define _hop_y_p_post()					\
  _sse_load_up(phi[ix]->s0);				\
  _sse_load(rs.s0);					\
  _sse_vector_add();					\
  _sse_store(rs.s0);					\
  _sse_load(rs.s3);					\
  _sse_vector_add();					\
  _sse_store(rs.s3);					\
  _sse_load_up(phi[ix]->s1);				\
  _sse_load(rs.s1);					\
  _sse_vector_add();					\
  _sse_store(rs.s1);					\
  _sse_load(rs.s2);					\
  _sse_vector_sub();					\
  _sse_store(rs.s2);      

#define _hop_y_m_post()				\
  _prefetch_su3(U+predist);			\
  _sse_load(phi[ix]->s0);			\
  _sse_su3_inverse_multiply((*U));		\
  _sse_vector_cmplxcg_mul(ka2);			\
  _sse_load(rs.s0);				\
  _sse_vector_add();				\
  _sse_store(rs.s0);				\
  _sse_load(rs.s3);				\
  _sse_vector_sub();				\
  _sse_store(rs.s3);				\
  _sse_load(phi[ix]->s1);			\
  _sse_su3_inverse_multiply((*U));		\
  _sse_vector_cmplxcg_mul(ka2);			\
  _sse_load(rs.s1);				\
  _sse_vector_add();				\
  _sse_store(rs.s1);				\
  _sse_load(rs.s2);				\
  _sse_vector_add();				\
  _sse_store(rs.s2);

#define _hop_z_p_post()					\
  _sse_load_up(phi[ix]->s0);				\
  _sse_load(rs.s0);					\
  _sse_vector_add();					\
  _sse_store(rs.s0);					\
  _sse_load(rs.s2);					\
  _sse_vector_i_mul();					\
  _sse_vector_sub();					\
  _sse_store(rs.s2);					\
  _sse_load_up(phi[ix]->s1);				\
  _sse_load(rs.s1);					\
  _sse_vector_add();					\
  _sse_store(rs.s1);					\
  _sse_load(rs.s3);					\
  _sse_vector_i_mul();					\
  _sse_vector_add();					\
  _sse_store(rs.s3);

#define _hop_z_m_post()				\
  _prefetch_su3(U+predist);			\
  _prefetch_spinor(s+1);			\
  _sse_load(phi[ix]->s0);			\
  _sse_su3_inverse_multiply((*U));		\
  _sse_vector_cmplxcg_mul(ka3);			\
  _sse_load(rs.s0);				\
  _sse_vector_add();				\
  _sse_store_nt(s->s0);				\
  _sse_load(rs.s2);				\
  _sse_vector_i_mul();				\
  _sse_vector_add();				\
  _sse_store_nt(s->s2);				\
  _sse_load(phi[ix]->s1);			\
  _sse_su3_inverse_multiply((*U));		\
  _sse_vector_cmplxcg_mul(ka3);			\
  _sse_load(rs.s1);				\
  _sse_vector_add();				\
  _sse_store_nt(s->s1);				\
  _sse_load(rs.s3);				\
  _sse_vector_i_mul();				\
  _sse_vector_sub();				\
  _sse_store_nt(s->s3);

#define _hop_mul_g5_cmplx_and_store(res)	\
  _sse_load_up((res)->s0);			\
  _sse_vector_cmplx_mul(cf);			\
  _sse_store_nt_up((res)->s0);			\
  _sse_load_up((res)->s1);			\
  _sse_vector_cmplx_mul(cf);			\
  _sse_store_nt_up((res)->s1);			\
  _sse_load_up((res)->s2);			\
  _sse_vector_cmplxcg_mul(cf);			\
  _sse_store_nt_up((res)->s2);			\
  _sse_load_up((res)->s3);			\
  _sse_vector_cmplxcg_mul(cf);			\
  _sse_store_nt_up((res)->s3);

#define _g5_cmplx_sub_hop_and_g5store(res)		\
  _sse_load_up(pn->s0);					\
  _sse_vector_cmplx_mul(cf);				\
  _sse_load((res)->s0);					\
  _sse_vector_sub_up();					\
  _sse_store_nt_up((res)->s0);				\
  _sse_load_up(pn->s1);					\
  _sse_vector_cmplx_mul(cf);				\
  _sse_load((res)->s1);					\
  _sse_vector_sub_up();					\
  _sse_store_nt_up((res)->s1);				\
  _sse_load_up(pn->s2);					\
  _sse_vector_cmplxcg_mul(cf);				\
  _sse_load((res)->s2);					\
  _sse_vector_sub();					\
  _sse_store_nt((res)->s2);				\
  _sse_load_up(pn->s3);					\
  _sse_vector_cmplxcg_mul(cf);				\
  _sse_load((res)->s3);					\
  _sse_vector_sub();					\
  _sse_store_nt((res)->s3);

#define _hop_store_post(res)

#if defined OPTERON
#  define _declare_hregs()			\
  spinor rs ALIGN;				\
  const int predist=2;
#else
#  define _declare_hregs()			\
  spinor rs ALIGN;				\
  const int predist=1;
#endif

#elif (defined BGL && defined XLC)

#define _declare_hregs()					\
  double _Complex reg00, reg01, reg02, reg03, reg04, reg05;	\
  double _Complex reg10, reg11, reg12, reg13, reg14, reg15;	\
  double _Complex u00, u01, u02, u10, u11, u12;			\
  double _Complex reg20, reg21;					\
  double _Complex rs00, rs01, rs02, rs10, rs11, rs12, rs20, rs21, rs22, \
    rs30, rs31, rs32;

#define _hop_t_p_pre32()					\
  _bgl_load_rs0(s->s0);						\
  _bgl_load_rs1(s->s1);						\
  _bgl_load_rs2(s->s2);						\
  _bgl_load_rs3(s->s3);						\
  _prefetch_spinor(s+1);					\
  _prefetch_su3(U+1);						\
  _bgl_vector_add_rs2_to_rs0_reg0();				\
  _bgl_vector_add_rs3_to_rs1_reg1();				\
  _bgl_su3_multiply_double((*U));				\
  _bgl_vector_cmplx_mul_double(ka0);				\
  _bgl_store_reg0_up_32(phi32[ix]->s0);				\
  _bgl_store_reg1_up_32(phi32[ix]->s1);

#define _hop_t_m_pre32()					\
  _bgl_vector_sub_rs2_from_rs0_reg0();				\
  _bgl_vector_sub_rs3_from_rs1_reg1();				\
  _bgl_store_reg0_32(phi32[ix]->s0);				\
  _bgl_store_reg1_32(phi32[ix]->s1);

#define _hop_x_p_pre32()					\
  _prefetch_su3(U+1);						\
  _bgl_vector_i_mul_add_rs3_to_rs0_reg0();			\
  _bgl_vector_i_mul_add_rs2_to_rs1_reg1();			\
  _bgl_su3_multiply_double((*U));				\
  _bgl_vector_cmplx_mul_double(ka1);				\
  _bgl_store_reg0_up_32(phi32[ix]->s0);				\
  _bgl_store_reg1_up_32(phi32[ix]->s1);

#define _hop_x_m_pre32()					\
  _bgl_vector_i_mul_sub_rs3_from_rs0_reg0();			\
  _bgl_vector_i_mul_sub_rs2_from_rs1_reg1();			\
  _bgl_store_reg0_32(phi32[ix]->s0);				\
  _bgl_store_reg1_32(phi32[ix]->s1);

#define _hop_y_p_pre32()					\
  _prefetch_su3(U+1);						\
  _bgl_vector_add_rs3_to_rs0_reg0();				\
  _bgl_vector_sub_rs2_from_rs1_reg1();				\
  _bgl_su3_multiply_double((*U));				\
  _bgl_vector_cmplx_mul_double(ka2);				\
  _bgl_store_reg0_up_32(phi32[ix]->s0);				\
  _bgl_store_reg1_up_32(phi32[ix]->s1);

#define _hop_y_m_pre32()					\
  _bgl_vector_sub_rs3_from_rs0_reg0();				\
  _bgl_vector_add_rs2_to_rs1_reg1();				\
  _bgl_store_reg0_32(phi32[ix]->s0);				\
  _bgl_store_reg1_32(phi32[ix]->s1);

#define _hop_z_p_pre32()					\
  _bgl_vector_i_mul_add_rs2_to_rs0_reg0();			\
  _bgl_vector_i_mul_sub_rs3_from_rs1_reg1();			\
  _bgl_su3_multiply_double((*U));				\
  _bgl_vector_cmplx_mul_double(ka3);				\
  _bgl_store_reg0_up_32(phi32[ix]->s0);				\
  _bgl_store_reg1_up_32(phi32[ix]->s1);

#define _hop_z_m_pre32()					\
  _bgl_vector_i_mul_sub_rs2_from_rs0_reg0();			\
  _bgl_vector_i_mul_add_rs3_to_rs1_reg1();			\
  _bgl_store_reg0_32(phi32[ix]->s0);				\
  _bgl_store_reg1_32(phi32[ix]->s1);

#define _hop_t_p_post32()					\
  _bgl_load_rs0_32(phi32[ix]->s0);				\
  rs20 = rs00;							\
  rs21 = rs01;							\
  rs22 = rs02;							\
  _bgl_load_rs1_32(phi32[ix]->s1);				\
  rs30 = rs10;							\
  rs31 = rs11;							\
  rs32 = rs12;

#define _hop_t_m_post32()					\
  _prefetch_su3(U+1);						\
  _bgl_load_reg0_32(phi32[ix]->s0);				\
  _bgl_load_reg1_32(phi32[ix]->s1);				\
  _bgl_su3_inverse_multiply_double((*U));			\
  _bgl_vector_cmplxcg_mul_double(ka0);				\
  _bgl_add_to_rs0_reg0();					\
  _bgl_sub_from_rs2_reg0();					\
  _bgl_add_to_rs1_reg1();					\
  _bgl_sub_from_rs3_reg1();

#define _hop_x_p_post32()					\
  _bgl_load_reg0_up_32(phi32[ix]->s0);				\
  _bgl_load_reg1_up_32(phi32[ix]->s1);				\
  _bgl_add_to_rs0_reg0();					\
  _bgl_i_mul_sub_from_rs3_reg0();				\
  _bgl_add_to_rs1_reg1();					\
  _bgl_i_mul_sub_from_rs2_reg1();

#define _hop_x_m_post32()					\
  _prefetch_su3(U+1);						\
  _bgl_load_reg0_32(phi32[ix]->s0);				\
  _bgl_load_reg1_32(phi32[ix]->s1);				\
  _bgl_su3_inverse_multiply_double((*U));			\
  _bgl_vector_cmplxcg_mul_double(ka1);				\
  _bgl_add_to_rs0_reg0();					\
  _bgl_add_to_rs1_reg1();					\
  _bgl_i_mul_add_to_rs3_reg0();					\
  _bgl_i_mul_add_to_rs2_reg1();

#define _hop_y_p_post32()					\
  _bgl_load_reg0_up_32(phi32[ix]->s0);				\
  _bgl_load_reg1_up_32(phi32[ix]->s1);				\
  _bgl_add_to_rs0_reg0();					\
  _bgl_add_to_rs1_reg1();					\
  _bgl_sub_from_rs2_reg1();					\
  _bgl_add_to_rs3_reg0();

#define _hop_y_m_post32()					\
  _prefetch_su3(U+1);						\
  _bgl_load_reg0_32(phi32[ix]->s0);				\
  _bgl_load_reg1_32(phi32[ix]->s1);				\
  _bgl_su3_inverse_multiply_double((*U));			\
  _bgl_vector_cmplxcg_mul_double(ka2);				\
  _bgl_add_to_rs0_reg0();					\
  _bgl_add_to_rs1_reg1();					\
  _bgl_add_to_rs2_reg1();					\
  _bgl_sub_from_rs3_reg0();

#define _hop_z_p_post32()					\
  _bgl_load_reg0_up_32(phi32[ix]->s0);				\
  _bgl_load_reg1_up_32(phi32[ix]->s1);				\
  _bgl_add_to_rs0_reg0();					\
  _bgl_add_to_rs1_reg1();					\
  _bgl_i_mul_sub_from_rs2_reg0();				\
  _bgl_i_mul_add_to_rs3_reg1();

#define _hop_z_m_post32()					\
  _prefetch_su3(U+1);						\
  _bgl_load_reg0_32(phi32[ix]->s0);				\
  _bgl_load_reg1_32(phi32[ix]->s1);				\
  _bgl_su3_inverse_multiply_double((*U));			\
  _bgl_vector_cmplxcg_mul_double(ka3);				\
  _bgl_add_to_rs0_reg0();					\
  _bgl_i_mul_add_to_rs2_reg0();					\
  _bgl_add_to_rs1_reg1();					\
  _bgl_i_mul_sub_from_rs3_reg1();

#define _hop_t_p_pre()				\
  _prefetch_halfspinor(phi[ix+4]);		\
  _bgl_load_rs0(s->s0);				\
  _bgl_load_rs1(s->s1);				\
  _bgl_load_rs2(s->s2);				\
  _bgl_load_rs3(s->s3);				\
  _prefetch_spinor(s+1);			\
  _prefetch_su3(U+1);				\
  _bgl_vector_add_rs2_to_rs0_reg0();		\
  _bgl_vector_add_rs3_to_rs1_reg1();		\
  _bgl_su3_multiply_double((*U));		\
  _bgl_vector_cmplx_mul_double(ka0);		\
  _bgl_store_reg0_up(phi[ix]->s0);		\
  _bgl_store_reg1_up(phi[ix]->s1);

#define _hop_t_m_pre()				\
  _prefetch_halfspinor(phi[ix+4]);		\
  _bgl_vector_sub_rs2_from_rs0_reg0();		\
  _bgl_vector_sub_rs3_from_rs1_reg1();		\
  _bgl_store_reg0(phi[ix]->s0);			\
  _bgl_store_reg1(phi[ix]->s1);

#define _hop_x_p_pre()				\
  _prefetch_halfspinor(phi[ix+4]);		\
  _prefetch_su3(U+1);				\
  _bgl_vector_i_mul_add_rs3_to_rs0_reg0();	\
  _bgl_vector_i_mul_add_rs2_to_rs1_reg1();	\
  _bgl_su3_multiply_double((*U));		\
  _bgl_vector_cmplx_mul_double(ka1);		\
  _bgl_store_reg0_up(phi[ix]->s0);		\
  _bgl_store_reg1_up(phi[ix]->s1);

#define _hop_x_m_pre()				\
  _prefetch_halfspinor(phi[ix+4]);			\
  _bgl_vector_i_mul_sub_rs3_from_rs0_reg0();		\
  _bgl_vector_i_mul_sub_rs2_from_rs1_reg1();		\
  _bgl_store_reg0(phi[ix]->s0);				\
  _bgl_store_reg1(phi[ix]->s1);

#define _hop_y_p_pre()				\
  _prefetch_halfspinor(phi[ix+4]);		\
  _prefetch_su3(U+1);				\
  _bgl_vector_add_rs3_to_rs0_reg0();		\
  _bgl_vector_sub_rs2_from_rs1_reg1();		\
  _bgl_su3_multiply_double((*U));		\
  _bgl_vector_cmplx_mul_double(ka2);		\
  _bgl_store_reg0_up(phi[ix]->s0);		\
  _bgl_store_reg1_up(phi[ix]->s1);

#define _hop_y_m_pre()				\
  _prefetch_halfspinor(phi[ix+4]);		\
  _bgl_vector_sub_rs3_from_rs0_reg0();		\
  _bgl_vector_add_rs2_to_rs1_reg1();		\
  _bgl_store_reg0(phi[ix]->s0);			\
  _bgl_store_reg1(phi[ix]->s1);

#define _hop_z_p_pre()				\
  _prefetch_halfspinor(phi[ix+4]);		\
  _prefetch_su3(U+1);				\
  _bgl_vector_i_mul_add_rs2_to_rs0_reg0();		\
  _bgl_vector_i_mul_sub_rs3_from_rs1_reg1();		\
  _bgl_su3_multiply_double((*U));			\
  _bgl_vector_cmplx_mul_double(ka3);			\
  _bgl_store_reg0_up(phi[ix]->s0);			\
  _bgl_store_reg1_up(phi[ix]->s1);

#define _hop_z_m_pre()				\
  _prefetch_halfspinor(phi[ix+4]);			\
  _bgl_vector_i_mul_sub_rs2_from_rs0_reg0();		\
  _bgl_vector_i_mul_add_rs3_to_rs1_reg1();		\
  _bgl_store_reg0(phi[ix]->s0);				\
  _bgl_store_reg1(phi[ix]->s1);

#define _hop_t_p_post();			\
  _prefetch_halfspinor(phi[ix+3]);		\
  _bgl_load_rs0(phi[ix]->s0);			\
  rs20 = rs00;					\
  rs21 = rs01;					\
  rs22 = rs02;					\
  _bgl_load_rs1(phi[ix]->s1);			\
  rs30 = rs10;					\
  rs31 = rs11;					\
  rs32 = rs12;

#define _hop_t_m_post();			\
  _prefetch_halfspinor(phi[ix+3]);		\
  _prefetch_su3(U+1);				\
  _bgl_load_reg0(phi[ix]->s0);			\
  _bgl_load_reg1(phi[ix]->s1);			\
  _bgl_su3_inverse_multiply_double((*U));	\
  _bgl_vector_cmplxcg_mul_double(ka0);		\
  _bgl_add_to_rs0_reg0();			\
  _bgl_sub_from_rs2_reg0();			\
  _bgl_add_to_rs1_reg1();			\
  _bgl_sub_from_rs3_reg1();

#define _hop_x_p_post();			\
  _prefetch_halfspinor(phi[ix+3]);		\
  _bgl_load_reg0_up(phi[ix]->s0);		\
  _bgl_load_reg1_up(phi[ix]->s1);		\
  _bgl_add_to_rs0_reg0();			\
  _bgl_i_mul_sub_from_rs3_reg0();		\
  _bgl_add_to_rs1_reg1();			\
  _bgl_i_mul_sub_from_rs2_reg1();

#define _hop_x_m_post();			\
  _prefetch_halfspinor(phi[ix+3]);		\
  _prefetch_su3(U+1);				\
  _bgl_load_reg0(phi[ix]->s0);			\
  _bgl_load_reg1(phi[ix]->s1);			\
  _bgl_su3_inverse_multiply_double((*U));	\
  _bgl_vector_cmplxcg_mul_double(ka1);		\
  _bgl_add_to_rs0_reg0();			\
  _bgl_add_to_rs1_reg1();			\
  _bgl_i_mul_add_to_rs3_reg0();			\
  _bgl_i_mul_add_to_rs2_reg1();      

#define _hop_y_p_post();			\
  _prefetch_halfspinor(phi[ix+3]);		\
  _bgl_load_reg0_up(phi[ix]->s0);		\
  _bgl_load_reg1_up(phi[ix]->s1);		\
  _bgl_add_to_rs0_reg0();			\
  _bgl_add_to_rs1_reg1();			\
  _bgl_sub_from_rs2_reg1();			\
  _bgl_add_to_rs3_reg0();

#define _hop_y_m_post();			\
  _prefetch_halfspinor(phi[ix+3]);		\
  _prefetch_su3(U+1);				\
  _bgl_load_reg0(phi[ix]->s0);			\
  _bgl_load_reg1(phi[ix]->s1);			\
  _bgl_su3_inverse_multiply_double((*U));	\
  _bgl_vector_cmplxcg_mul_double(ka2);		\
  _bgl_add_to_rs0_reg0();			\
  _bgl_add_to_rs1_reg1();			\
  _bgl_add_to_rs2_reg1();			\
  _bgl_sub_from_rs3_reg0();

#define _hop_z_p_post();			\
  _prefetch_halfspinor(phi[ix+3]);		\
  _bgl_load_reg0_up(phi[ix]->s0);		\
  _bgl_load_reg1_up(phi[ix]->s1);		\
  _bgl_add_to_rs0_reg0();			\
  _bgl_add_to_rs1_reg1();			\
  _bgl_i_mul_sub_from_rs2_reg0();		\
  _bgl_i_mul_add_to_rs3_reg1();

#define _hop_z_m_post();			\
  _prefetch_spinor(s);				\
  _prefetch_halfspinor(phi[ix+3]);		\
  _prefetch_su3(U+1);				\
  _bgl_load_reg0(phi[ix]->s0);			\
  _bgl_load_reg1(phi[ix]->s1);			\
  _bgl_su3_inverse_multiply_double((*U));	\
  _bgl_vector_cmplxcg_mul_double(ka3);		\
  _bgl_add_to_rs0_reg0();			\
  _bgl_i_mul_add_to_rs2_reg0();			\
  _bgl_add_to_rs1_reg1();			\
  _bgl_i_mul_sub_from_rs3_reg1();



#define _hop_store_post(res)			\
  _bgl_store_rs0((res)->s0);			\
  _bgl_store_rs1((res)->s1);			\
  _bgl_store_rs2((res)->s2);			\
  _bgl_store_rs3((res)->s3);


#elif (defined BGQ && defined XLC)

#define _hop_t_p_pre32()					\
  _vec_load2(rs0, rs1, rs2, s->s0);				\
  _vec_load2(rs3, rs4, rs5, s->s1);				\
  _vec_load2(rs6, rs7, rs8, s->s2);				\
  _vec_load2(rs9, rs10, rs11, s->s3);				\
  _prefetch_spinor(s+1);					\
  _prefetch_su3(U+1);						\
  _vec_add_to2(r0, r1, r2, rs0, rs1, rs2, rs6, rs7, rs8);	\
  _vec_add_to2(r3, r4, r5, rs3, rs4, rs5, rs9, rs10, rs11);	\
  rtmp = vec_ld2(0, (double*) &ka0);					\
  _vec_su3_multiply_double2(U);						\
  _vec_cmplx_mul_double2(r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, rtmp); \
  _vec_store2_32(phi32[ix]->s0, r0, r1, r2);				\
  _vec_store2_32(phi32[ix]->s1, r3, r4, r5);

#define _hop_t_m_pre32()					\
  _vec_sub_to2(r0, r1, r2, rs0, rs1, rs2, rs6, rs7, rs8);	\
  _vec_sub_to2(r3, r4, r5, rs3, rs4, rs5, rs9, rs10, rs11);	\
  _vec_store2_32(phi32[ix]->s0, r0, r1, r2);			\
  _vec_store2_32(phi32[ix]->s1, r3, r4, r5);

#define _hop_x_p_pre32()						\
  _prefetch_su3(U+1);							\
  _vec_i_mul_add_to2(r0, r1, r2, rs0, rs1, rs2, rs9, rs10, rs11, U0);	\
  _vec_i_mul_add_to2(r3, r4, r5, rs3, rs4, rs5, rs6, rs7, rs8, U0);	\
  rtmp = vec_ld2(0, (double*) &ka1);					\
  _vec_su3_multiply_double2(U);						\
  _vec_cmplx_mul_double2(r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, rtmp); \
  _vec_store2_32(phi32[ix]->s0, r0, r1, r2);				\
  _vec_store2_32(phi32[ix]->s1, r3, r4, r5);

#define _hop_x_m_pre32()						\
  _vec_i_mul_sub_to2(r0, r1, r2, rs0, rs1, rs2, rs9, rs10, rs11, U0);	\
  _vec_i_mul_sub_to2(r3, r4, r5, rs3, rs4, rs5, rs6, rs7, rs8, U0);	\
  _vec_store2_32(phi32[ix]->s0, r0, r1, r2);				\
  _vec_store2_32(phi32[ix]->s1, r3, r4, r5);

#define _hop_y_p_pre32()					\
  _prefetch_su3(U+1);						\
  _vec_add_to2(r0, r1, r2, rs0, rs1, rs2, rs9, rs10, rs11);	\
  _vec_sub_to2(r3, r4, r5, rs3, rs4, rs5, rs6, rs7, rs8);	\
  rtmp = vec_ld2(0, (double*) &ka2);					\
  _vec_su3_multiply_double2(U);						\
  _vec_cmplx_mul_double2(r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, rtmp); \
  _vec_store2_32(phi32[ix]->s0, r0, r1, r2);				\
  _vec_store2_32(phi32[ix]->s1, r3, r4, r5);

#define _hop_y_m_pre32()					\
  _vec_sub_to2(r0, r1, r2, rs0, rs1, rs2, rs9, rs10, rs11);	\
  _vec_add_to2(r3, r4, r5, rs3, rs4, rs5, rs6, rs7, rs8);	\
  _vec_store2_32(phi32[ix]->s0, r0, r1, r2);			\
  _vec_store2_32(phi32[ix]->s1, r3, r4, r5);

#define _hop_z_p_pre32()						\
  _prefetch_su3(U+1);							\
  _vec_i_mul_add_to2(r0, r1, r2, rs0, rs1, rs2, rs6, rs7, rs8, U0);	\
  _vec_i_mul_sub_to2(r3, r4, r5, rs3, rs4, rs5, rs9, rs10, rs11, U0);	\
  rtmp = vec_ld2(0, (double*) &ka3);					\
  _vec_su3_multiply_double2(U);						\
  _vec_cmplx_mul_double2(r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, rtmp); \
  _vec_store2_32(phi32[ix]->s0, r0, r1, r2);				\
  _vec_store2_32(phi32[ix]->s1, r3, r4, r5);

#define _hop_z_m_pre32()						\
  _vec_i_mul_sub_to2(r0, r1, r2, rs0, rs1, rs2, rs6, rs7, rs8, U0);	\
  _vec_i_mul_add_to2(r3, r4, r5, rs3, rs4, rs5, rs9, rs10, rs11, U0);	\
  _vec_store2_32(phi32[ix]->s0, r0, r1, r2);				\
  _vec_store2_32(phi32[ix]->s1, r3, r4, r5);

#define _hop_t_p_post32()				\
  _vec_load2_32(rs0, rs1, rs2, phi32[ix]->s0);		\
  rs6 = rs0;						\
  rs7 = rs1;						\
  rs8 = rs2;						\
  _vec_load2_32(rs3, rs4, rs5, phi32[ix]->s1);		\
  rs9 = rs3;						\
  rs10= rs4;						\
  rs11= rs5;

#define _hop_t_m_post32()						\
  _prefetch_su3(U+1);							\
  _vec_load2_32(r0, r1, r2, phi32[ix]->s0);				\
  _vec_load2_32(r3, r4, r5, phi32[ix]->s1);				\
  rtmp = vec_ld2(0, (double*) &ka0);					\
  _vec_su3_inverse_multiply_double2(U);					\
  _vec_cmplxcg_mul_double2(r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, rtmp); \
  _vec_add2(rs0, rs1, rs2, r0, r1, r2);					\
  _vec_sub2(rs6, rs7, rs8, r0, r1, r2);					\
  _vec_add2(rs3, rs4, rs5, r3, r4, r5);					\
  _vec_sub2(rs9, rs10, rs11, r3, r4, r5);

#define _hop_x_p_post32()						\
  _vec_load2_32(r0, r1, r2, phi32[ix]->s0);				\
  _vec_load2_32(r3, r4, r5, phi32[ix]->s1);				\
  _vec_add2(rs0, rs1, rs2, r0, r1, r2);					\
  _vec_i_mul_sub2(rs9, rs10, rs11, r0, r1, r2, U0);			\
  _vec_add2(rs3, rs4, rs5, r3, r4, r5);					\
  _vec_i_mul_sub2(rs6, rs7, rs8, r3, r4, r5, U0);

#define _hop_x_m_post32()						\
  _prefetch_su3(U+1);							\
  _vec_load2_32(r0, r1, r2, phi32[ix]->s0);				\
  _vec_load2_32(r3, r4, r5, phi32[ix]->s1);				\
  rtmp = vec_ld2(0, (double*) &ka1);					\
  _vec_su3_inverse_multiply_double2(U);					\
  _vec_cmplxcg_mul_double2(r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, rtmp); \
  _vec_add_double2(rs0, rs1, rs2, rs3, rs4, rs5, r0, r1, r2, r3, r4, r5); \
  _vec_i_mul_add_double2(rs9, rs10, rs11, rs6, rs7, rs8, r0, r1, r2, r3, r4, r5, U0);

#define _hop_y_p_post32()						\
  _vec_load2_32(r0, r1, r2, phi32[ix]->s0);				\
  _vec_load2_32(r3, r4, r5, phi32[ix]->s1);				\
  _vec_add_double2(rs0, rs1, rs2, rs3, rs4, rs5, r0, r1, r2, r3, r4, r5); \
  _vec_sub2(rs6, rs7, rs8, r3, r4, r5);					\
  _vec_add2(rs9, rs10, rs11, r0, r1, r2);

#define _hop_y_m_post32()						\
  _prefetch_su3(U+1);							\
  _vec_load2_32(r0, r1, r2, phi32[ix]->s0);				\
  _vec_load2_32(r3, r4, r5, phi32[ix]->s1);				\
  rtmp = vec_ld2(0, (double*) &ka2);					\
  _vec_su3_inverse_multiply_double2(U);					\
  _vec_cmplxcg_mul_double2(r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, rtmp); \
  _vec_add_double2(rs0, rs1, rs2, rs3, rs4, rs5, r0, r1, r2, r3, r4, r5); \
  _vec_add2(rs6, rs7, rs8, r3, r4, r5);					\
  _vec_sub2(rs9, rs10, rs11, r0, r1, r2);

#define _hop_z_p_post32()						\
  _vec_load2_32(r0, r1, r2, phi32[ix]->s0);				\
  _vec_load2_32(r3, r4, r5, phi32[ix]->s1);				\
  _vec_add_double2(rs0, rs1, rs2, rs3, rs4, rs5, r0, r1, r2, r3, r4, r5); \
  _vec_i_mul_sub2(rs6, rs7, rs8, r0, r1, r2, U0);			\
  _vec_i_mul_add2(rs9, rs10, rs11, r3, r4, r5, U0);

#define _hop_z_m_post32()						\
  _prefetch_su3(U+1);							\
  _vec_load2_32(r0, r1, r2, phi32[ix]->s0);				\
  _vec_load2_32(r3, r4, r5, phi32[ix]->s1);				\
  rtmp = vec_ld2(0, (double*) &ka3);					\
  _vec_su3_inverse_multiply_double2(U);					\
  _vec_cmplxcg_mul_double2(r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, rtmp); \
  _vec_add2(rs0, rs1, rs2, r0, r1, r2);					\
  _vec_i_mul_add2(rs6, rs7, rs8, r0, r1, r2, U0);			\
  _vec_add2(rs3, rs4, rs5, r3, r4, r5);					\
  _vec_i_mul_sub2(rs9, rs10, rs11, r3, r4, r5, U0);

#define _hop_t_p_pre2()						\
  _vec_load2(rs0, rs1, rs2, s->s0);				\
  _vec_load2(rs3, rs4, rs5, s->s1);				\
  _vec_load2(rs6, rs7, rs8, s->s2);				\
  _vec_load2(rs9, rs10, rs11, s->s3);				\
  _prefetch_spinor(s+1);					\
  _prefetch_su3(U+1);						\
  _vec_add_to2(r0, r1, r2, rs0, rs1, rs2, rs6, rs7, rs8);	\
  _vec_add_to2(r3, r4, r5, rs3, rs4, rs5, rs9, rs10, rs11);	\
  rtmp = vec_ld2(0, (double*) &ka0);					\
  _vec_su3_multiply_double2(U);						\
  _vec_cmplx_mul_double2(r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, rtmp); \
  _vec_store2(phi[ix]->s0, r0, r1, r2);					\
  _vec_store2(phi[ix]->s1, r3, r4, r5);

// note that this _hop_t_p_pre stores the quadwords in phi[ix] 
// in a different order than expected, but this is taken care of in 
// the correspondin _hop_t_p_post version
//
// it might be good to check whether unfusing is better done here
// instead of in the corresponding post version!?
#define _hop_t_p_pre()							\
  _vec_load(rs0, rs1, s->s0);						\
  _vec_load16(rs2, rs3, s->s1, rtmp);					\
  _vec_load(rs4, rs5, s->s2);						\
  _vec_load16(rs6, rs7, s->s3, rtmp);					\
  _prefetch_spinor(s+1);						\
  _prefetch_su3(U+1);							\
  _vec_add(r0, r1, rs0, rs1, rs4, rs5);					\
  _vec_add(r2, r3, rs2, rs3, rs6, rs7);					\
  _vec_su3_multiply_double2c(U);					\
  rtmp = vec_ld2(0, (double*) &ka0);					\
  _vec_cmplx_mul_double2c(r0, r1, r2, r4, r5, r6, rtmp);		\
  _vec_store_halfspinor(phi[ix]->s0, r0, r1, r2);

#define _hop_t_m_pre2()						\
  _vec_sub_to2(r0, r1, r2, rs0, rs1, rs2, rs6, rs7, rs8);	\
  _vec_sub_to2(r3, r4, r5, rs3, rs4, rs5, rs9, rs10, rs11);	\
  _vec_store2(phi[ix]->s0, r0, r1, r2);				\
  _vec_store2(phi[ix]->s1, r3, r4, r5);

#define _hop_t_m_pre()						\
  _vec_sub(r0, r1, rs0, rs1, rs4, rs5);				\
  _vec_sub(r2, r3, rs2, rs3, rs6, rs7);				\
  _vec_store(phi[ix]->s0, r0, r1);				\
  _vec_store16(phi[ix]->s1, r2, r3, U0);

#define _hop_x_p_pre2()							\
  _prefetch_su3(U+1);							\
  _vec_i_mul_add_to2(r0, r1, r2, rs0, rs1, rs2, rs9, rs10, rs11, U0);	\
  _vec_i_mul_add_to2(r3, r4, r5, rs3, rs4, rs5, rs6, rs7, rs8, U0);	\
  rtmp = vec_ld2(0, (double*) &ka1);					\
  _vec_su3_multiply_double2(U);						\
  _vec_cmplx_mul_double2(r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, rtmp); \
  _vec_store2(phi[ix]->s0, r0, r1, r2);					\
  _vec_store2(phi[ix]->s1, r3, r4, r5);

#define _hop_x_p_pre()						\
  _prefetch_su3(U+1);						\
  _vec_i_mul_add(r0, r1, rs0, rs1, rs6, rs7, U0);		\
  _vec_i_mul_add(r2, r3, rs2, rs3, rs4, rs5, U0);		\
  rtmp = vec_ld2(0, (double*) &ka1);				\
  _vec_su3_multiply_double2c(U);				\
  _vec_cmplx_mul_double2c(r0, r1, r2, r4, r5, r6, rtmp);	\
  _vec_store_halfspinor(phi[ix]->s0, r0, r1, r2);

#define _hop_x_m_pre2()							\
  _vec_i_mul_sub_to2(r0, r1, r2, rs0, rs1, rs2, rs9, rs10, rs11, U0);	\
  _vec_i_mul_sub_to2(r3, r4, r5, rs3, rs4, rs5, rs6, rs7, rs8, U0);	\
  _vec_store2(phi[ix]->s0, r0, r1, r2);					\
  _vec_store2(phi[ix]->s1, r3, r4, r5);

#define _hop_x_m_pre()					\
  _vec_i_mul_sub(r0, r1, rs0, rs1, rs6, rs7, U0);	\
  _vec_i_mul_sub(r2, r3, rs2, rs3, rs4, rs5, U0);	\
  _vec_store(phi[ix]->s0, r0, r1);			\
  _vec_store16(phi[ix]->s1, r2, r3, U0);

#define _hop_y_p_pre2()						\
  _prefetch_su3(U+1);						\
  _vec_add_to2(r0, r1, r2, rs0, rs1, rs2, rs9, rs10, rs11);	\
  _vec_sub_to2(r3, r4, r5, rs3, rs4, rs5, rs6, rs7, rs8);	\
  rtmp = vec_ld2(0, (double*) &ka2);					\
  _vec_su3_multiply_double2(U);						\
  _vec_cmplx_mul_double2(r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, rtmp); \
  _vec_store2(phi[ix]->s0, r0, r1, r2);					\
  _vec_store2(phi[ix]->s1, r3, r4, r5);

#define _hop_y_p_pre()						\
  _prefetch_su3(U+1);						\
  _vec_add(r0, r1, rs0, rs1, rs6, rs7);				\
  _vec_sub(r2, r3, rs2, rs3, rs4, rs5);				\
  rtmp = vec_ld2(0, (double*) &ka2);				\
  _vec_su3_multiply_double2c(U);				\
  _vec_cmplx_mul_double2c(r0, r1, r2, r4, r5, r6, rtmp);	\
  _vec_store_halfspinor(phi[ix]->s0, r0, r1, r2);

#define _hop_y_m_pre2()						\
  _vec_sub_to2(r0, r1, r2, rs0, rs1, rs2, rs9, rs10, rs11);	\
  _vec_add_to2(r3, r4, r5, rs3, rs4, rs5, rs6, rs7, rs8);	\
  _vec_store2(phi[ix]->s0, r0, r1, r2);				\
  _vec_store2(phi[ix]->s1, r3, r4, r5);

#define _hop_y_m_pre()				\
  _vec_sub(r0, r1, rs0, rs1, rs6, rs7);		\
  _vec_add(r2, r3, rs2, rs3, rs4, rs5);		\
  _vec_store(phi[ix]->s0, r0, r1);		\
  _vec_store16(phi[ix]->s1, r2, r3, U0);

#define _hop_z_p_pre2()							\
  _prefetch_su3(U+1);							\
  _vec_i_mul_add_to2(r0, r1, r2, rs0, rs1, rs2, rs6, rs7, rs8, U0);	\
  _vec_i_mul_sub_to2(r3, r4, r5, rs3, rs4, rs5, rs9, rs10, rs11, U0);	\
  rtmp = vec_ld2(0, (double*) &ka3);					\
  _vec_su3_multiply_double2(U);						\
  _vec_cmplx_mul_double2(r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, rtmp); \
  _vec_store2(phi[ix]->s0, r0, r1, r2);					\
  _vec_store2(phi[ix]->s1, r3, r4, r5);

#define _hop_z_p_pre()						\
  _prefetch_su3(U+1);						\
  _vec_i_mul_add(r0, r1, rs0, rs1, rs4, rs5, U0);		\
  _vec_i_mul_sub(r2, r3, rs2, rs3, rs6, rs7, U0);		\
  rtmp = vec_ld2(0, (double*) &ka3);				\
  _vec_su3_multiply_double2c(U);				\
  _vec_cmplx_mul_double2c(r0, r1, r2, r4, r5, r6, rtmp);	\
  _vec_store_halfspinor(phi[ix]->s0, r0, r1, r2);

#define _hop_z_m_pre2()							\
  _vec_i_mul_sub_to2(r0, r1, r2, rs0, rs1, rs2, rs6, rs7, rs8, U0);	\
  _vec_i_mul_add_to2(r3, r4, r5, rs3, rs4, rs5, rs9, rs10, rs11, U0);	\
  _vec_store2(phi[ix]->s0, r0, r1, r2);					\
  _vec_store2(phi[ix]->s1, r3, r4, r5);

#define _hop_z_m_pre()					\
  _vec_i_mul_sub(r0, r1, rs0, rs1, rs4, rs5, U0);	\
  _vec_i_mul_add(r2, r3, rs2, rs3, rs6, rs7, U0);	\
  _vec_store(phi[ix]->s0, r0, r1);			\
  _vec_store16(phi[ix]->s1, r2, r3, U0);

#define _hop_t_p_post2()				\
  _vec_load2(rs0, rs1, rs2, phi[ix]->s0);	\
  rs6 = rs0;					\
  rs7 = rs1;					\
  rs8 = rs2;					\
  _vec_load2(rs3, rs4, rs5, phi[ix]->s1);	\
  rs9 = rs3;					\
  rs10= rs4;					\
  rs11= rs5;

#define _hop_t_p_post()				\
  _vec_load_halfspinor(rs0, rs1, rs2, phi[ix]->s0);	\
  _vec_unfuse(rs0, rs1, rs2, rs3, rs4, rs5);		\
  rs6 = rs0; rs7 = rs1; rs8 = rs2;			\
  rs9 = rs3; rs10= rs4; rs11= rs5;

#define _hop_t_m_post2()							\
  _prefetch_su3(U+1);							\
  _vec_load2(r0, r1, r2, phi[ix]->s0);					\
  _vec_load2(r3, r4, r5, phi[ix]->s1);					\
  rtmp = vec_ld2(0, (double*) &ka0);					\
  _vec_su3_inverse_multiply_double2(U);					\
  _vec_cmplxcg_mul_double2(r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, rtmp); \
  _vec_add2(rs0, rs1, rs2, r0, r1, r2);					\
  _vec_sub2(rs6, rs7, rs8, r0, r1, r2);					\
  _vec_add2(rs3, rs4, rs5, r3, r4, r5);					\
  _vec_sub2(rs9, rs10, rs11, r3, r4, r5);


#define _hop_t_m_post()						\
  _prefetch_su3(U+1);							\
  _vec_load(r0, r1, phi[ix]->s0);					\
  _vec_load16(r2, r3, phi[ix]->s1, rtmp);				\
  rtmp = vec_ld2(0, (double*) &ka0);					\
  _vec_su3_inverse_multiply_double2c(U);				\
  _vec_cmplxcg_mul_double2c(r0, r1, r2, r4, r5, r6, rtmp);		\
  _vec_unfuse(r0, r1, r2, r3, r4, r5);					\
  _vec_add_double2(rs0, rs1, rs2, rs3, rs4, rs5, r0, r1, r2, r3, r4, r5); \
  _vec_sub_double2(rs6, rs7, rs8, rs9, rs10, rs11, r0, r1, r2, r3, r4, r5);

#define _hop_x_p_post2()					\
  _vec_load2(r0, r1, r2, phi[ix]->s0);			\
  _vec_load2(r3, r4, r5, phi[ix]->s1);			\
  _vec_add2(rs0, rs1, rs2, r0, r1, r2);			\
  _vec_i_mul_sub2(rs9, rs10, rs11, r0, r1, r2, U0);	\
  _vec_add2(rs3, rs4, rs5, r3, r4, r5);			\
  _vec_i_mul_sub2(rs6, rs7, rs8, r3, r4, r5, U0);

#define _hop_x_p_post()						\
  _vec_load_halfspinor(r0, r1, r2, phi[ix]->s0);			\
  _vec_unfuse(r0, r1, r2, r3, r4, r5);				\
  _vec_add_double2(rs0, rs1, rs2, rs3, rs4, rs5, r0, r1, r2, r3, r4, r5); \
  _vec_i_mul_sub2(rs6, rs7, rs8, r3, r4, r5, U0);			\
  _vec_i_mul_sub2(rs9, rs10, rs11, r0, r1, r2, U1);
  
#define _hop_x_m_post2()							\
  _prefetch_su3(U+1);							\
  _vec_load2(r0, r1, r2, phi[ix]->s0);					\
  _vec_load2(r3, r4, r5, phi[ix]->s1);					\
  rtmp = vec_ld2(0, (double*) &ka1);					\
  _vec_su3_inverse_multiply_double2(U);					\
  _vec_cmplxcg_mul_double2(r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, rtmp); \
  _vec_add_double2(rs0, rs1, rs2, rs3, rs4, rs5, r0, r1, r2, r3, r4, r5); \
  _vec_i_mul_add_double2(rs9, rs10, rs11, rs6, rs7, rs8, r0, r1, r2, r3, r4, r5, U0);

#define _hop_x_m_post()						\
  _prefetch_su3(U+1);							\
  _vec_load(r0, r1, phi[ix]->s0);					\
  _vec_load16(r2, r3, phi[ix]->s1, rtmp);				\
  rtmp = vec_ld2(0, (double*) &ka1);					\
  _vec_su3_inverse_multiply_double2c(U);				\
  _vec_cmplxcg_mul_double2c(r0, r1, r2, r4, r5, r6, rtmp);		\
  _vec_unfuse(r0, r1, r2, r3, r4, r5);					\
  _vec_add_double2(rs0, rs1, rs2, rs3, rs4, rs5, r0, r1, r2, r3, r4, r5); \
  _vec_i_mul_add_double2(rs9, rs10, rs11, rs6, rs7, rs8, r0, r1, r2, r3, r4, r5, U0);

#define _hop_y_p_post2()							\
  _vec_load2(r0, r1, r2, phi[ix]->s0);					\
  _vec_load2(r3, r4, r5, phi[ix]->s1);					\
  _vec_add_double2(rs0, rs1, rs2, rs3, rs4, rs5, r0, r1, r2, r3, r4, r5); \
  _vec_sub2(rs6, rs7, rs8, r3, r4, r5);					\
  _vec_add2(rs9, rs10, rs11, r0, r1, r2);

#define _hop_y_p_post()						\
  _vec_load_halfspinor(r0, r1, r2, phi[ix]->s0);			\
  _vec_unfuse(r0, r1, r2, r3, r4, r5);				\
  _vec_add_double2(rs0, rs1, rs2, rs3, rs4, rs5, r0, r1, r2, r3, r4, r5);	\
  _vec_sub2(rs6, rs7, rs8, r3, r4, r5);					\
  _vec_add2(rs9, rs10, rs11, r0, r1, r2);

#define _hop_y_m_post2()							\
  _prefetch_su3(U+1);							\
  _vec_load2(r0, r1, r2, phi[ix]->s0);					\
  _vec_load2(r3, r4, r5, phi[ix]->s1);					\
  rtmp = vec_ld2(0, (double*) &ka2);					\
  _vec_su3_inverse_multiply_double2(U);					\
  _vec_cmplxcg_mul_double2(r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, rtmp); \
  _vec_add_double2(rs0, rs1, rs2, rs3, rs4, rs5, r0, r1, r2, r3, r4, r5); \
  _vec_add2(rs6, rs7, rs8, r3, r4, r5);					\
  _vec_sub2(rs9, rs10, rs11, r0, r1, r2);

#define _hop_y_m_post()						\
  _prefetch_su3(U+1);							\
  _vec_load(r0, r1, phi[ix]->s0);					\
  _vec_load16(r2, r3, phi[ix]->s1, rtmp);				\
  rtmp = vec_ld2(0, (double*) &ka2);					\
  _vec_su3_inverse_multiply_double2c(U);				\
  _vec_cmplxcg_mul_double2c(r0, r1, r2, r4, r5, r6, rtmp);		\
  _vec_unfuse(r0, r1, r2, r3, r4, r5);					\
  _vec_add_double2(rs0, rs1, rs2, rs3, rs4, rs5, r0, r1, r2, r3, r4, r5); \
  _vec_add2(rs6, rs7, rs8, r3, r4, r5);					\
  _vec_sub2(rs9, rs10, rs11, r0, r1, r2);

#define _hop_z_p_post2()				\
  _vec_load2(r0, r1, r2, phi[ix]->s0);		\
  _vec_load2(r3, r4, r5, phi[ix]->s1);		\
  _vec_add_double2(rs0, rs1, rs2, rs3, rs4, rs5, r0, r1, r2, r3, r4, r5); \
  _vec_i_mul_sub2(rs6, rs7, rs8, r0, r1, r2, U0);			\
  _vec_i_mul_add2(rs9, rs10, rs11, r3, r4, r5, U0);

#define _hop_z_p_post()						\
  _vec_load_halfspinor(r0, r1, r2, phi[ix]->s0);			\
  _vec_unfuse(r0, r1, r2, r3, r4, r5);				\
  _vec_add_double2(rs0, rs1, rs2, rs3, rs4, rs5, r0, r1, r2, r3, r4, r5); \
  _vec_i_mul_sub2(rs6, rs7, rs8, r0, r1, r2, U0);			\
  _vec_i_mul_add2(rs9, rs10, rs11, r3, r4, r5, U1);

#define _hop_z_m_post2()							\
  _prefetch_su3(U+1);							\
  _vec_load2(r0, r1, r2, phi[ix]->s0);					\
  _vec_load2(r3, r4, r5, phi[ix]->s1);					\
  rtmp = vec_ld2(0, (double*) &ka3);					\
  _vec_su3_inverse_multiply_double2(U);					\
  _vec_cmplxcg_mul_double2(r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, rtmp); \
  _vec_add2(rs0, rs1, rs2, r0, r1, r2);					\
  _vec_add2(rs3, rs4, rs5, r3, r4, r5);					\
  _vec_i_mul_add2(rs6, rs7, rs8, r0, r1, r2, U0);			\
  _vec_i_mul_sub2(rs9, rs10, rs11, r3, r4, r5, U0);

#define _hop_z_m_post()						\
  _prefetch_su3(U+1);							\
  _vec_load(r0, r1, phi[ix]->s0);					\
  _vec_load16(r2, r3, phi[ix]->s1, rtmp);				\
  rtmp = vec_ld2(0, (double*) &ka3);					\
  _vec_su3_inverse_multiply_double2c(U);				\
  _vec_cmplxcg_mul_double2c(r0, r1, r2, r4, r5, r6, rtmp);		\
  _vec_unfuse(r0, r1, r2, r3, r4, r5);					\
  _vec_add_double2(rs0, rs1, rs2, rs3, rs4, rs5, r0, r1, r2, r3, r4, r5); \
  _vec_i_mul_add2(rs6, rs7, rs8, r0, r1, r2, U0);			\
  _vec_i_mul_sub2(rs9, rs10, rs11, r3, r4, r5, U1);

#define _hop_mul_g5_cmplx_and_store(res)					\
  _vec_cmplx_mul_double2(r0, r1, r2, r3, r4, r5, rs0, rs1, rs2, rs3, rs4, rs5, cf); \
  _vec_cmplxcg_mul_double2(r6, r7, r8, r9, r10, r11, rs6, rs7, rs8, rs9, rs10, rs11, cf); \
  _vec_store2((res)->s0, r0, r1, r2);					\
  _vec_store2((res)->s1, r3, r4, r5);					\
  _vec_store2((res)->s2, r6, r7, r8);					\
  _vec_store2((res)->s3, r9, r10, r11);

#define _g5_cmplx_sub_hop_and_g5store(res)					\
  _vec_load_halfspinor(r3, r4, r5, pn->s0);				\
  _vec_cmplx_mul_double2c(r0, r1, r2, r3, r4, r5, cf);			\
  _vec_unfuse(r0, r1, r2, r3, r4, r5);					\
  _vec_sub_double2(r0, r3, r1, r4, r2, r5, rs0, rs1, rs2, rs3, rs4, rs5); \
  _vec_store2((res)->s0, r0, r3, r1);					\
  _vec_store2((res)->s1, r4, r2, r5);					\
  _vec_load_halfspinor(r3, r4, r5, pn->s2);				\
  _vec_cmplxcg_mul_double2c(r0, r1, r2, r3, r4, r5, cf);		\
  _vec_unfuse(r0, r1, r2, r3, r4, r5);					\
  _vec_sub_double2(rs6, rs7, rs8, rs9, rs10, rs11, r0, r3, r1, r4, r2, r5); \
  _vec_store2((res)->s2, rs6, rs7, rs8);					\
  _vec_store2((res)->s3, rs9, rs10, rs11);

#define _hop_store_post(res)		\
  _vec_store2((res)->s0, rs0, rs1, rs2);	\
  _vec_store2((res)->s1, rs3, rs4, rs5);	\
  _vec_store2((res)->s2, rs6, rs7, rs8);	\
  _vec_store2((res)->s3, rs9, rs10, rs11);


#define _declare_hregs()						\
  vector4double ALIGN r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11;	\
  vector4double ALIGN rs0, rs1, rs2, rs3, rs4, rs5, rs6, rs7, rs8, rs9, rs10, rs11; \
  vector4double ALIGN U0, U1, U2, U3, U4, U6, U7;			\
  vector4double ALIGN rtmp;

#else

#define _prefetch_spinor(s)
#define _prefetch_halfspinor(hs)
#define _prefetch_su3(U)

#define _hop_t_p_pre32()				\
  _vector_assign(rs.s0, s->s0);				\
  _vector_assign(rs.s1, s->s1);				\
  _vector_assign(rs.s2, s->s2);				\
  _vector_assign(rs.s3, s->s3);				\
  _vector_add(psi, rs.s0, rs.s2);			\
  _su3_multiply(chi,(*U),psi);				\
  _complex_times_vector(phi32[ix]->s0, ka0, chi);	\
  _vector_add(psi, rs.s1, rs.s3);			\
  _su3_multiply(chi,(*U),psi);				\
  _complex_times_vector(phi32[ix]->s1, ka0, chi);

#define _hop_t_m_pre32()				\
  _vector_sub(phi32[ix]->s0, rs.s0, rs.s2);		\
  _vector_sub(phi32[ix]->s1, rs.s1, rs.s3);

#define _hop_x_p_pre32()				\
  _vector_i_add(psi, rs.s0, rs.s3);			\
  _su3_multiply(chi, (*U), psi);			\
  _complex_times_vector(phi32[ix]->s0, ka1, chi);	\
  _vector_i_add(psi, rs.s1, rs.s2);			\
  _su3_multiply(chi, (*U), psi);			\
  _complex_times_vector(phi32[ix]->s1, ka1, chi);

#define _hop_x_m_pre32()				\
  _vector_i_sub(phi32[ix]->s0, rs.s0, rs.s3);		\
  _vector_i_sub(phi32[ix]->s1, rs.s1, rs.s2);

#define _hop_y_p_pre32()				\
  _vector_add(psi, rs.s0, rs.s3);			\
  _su3_multiply(chi,(*U),psi);				\
  _complex_times_vector(phi32[ix]->s0, ka2, chi);	\
  _vector_sub(psi, rs.s1, rs.s2);			\
  _su3_multiply(chi,(*U),psi);				\
  _complex_times_vector(phi32[ix]->s1, ka2, chi);

#define _hop_y_m_pre32()			\
  _vector_sub(phi32[ix]->s0, rs.s0, rs.s3);	\
  _vector_add(phi32[ix]->s1, rs.s1, rs.s2);

#define _hop_z_p_pre32()				\
  _vector_i_add(psi, rs.s0, rs.s2);			\
  _su3_multiply(chi, (*U), psi);			\
  _complex_times_vector(phi32[ix]->s0, ka3, chi);	\
  _vector_i_sub(psi, rs.s1, rs.s3);			\
  _su3_multiply(chi,(*U),psi);				\
  _complex_times_vector(phi32[ix]->s1, ka3, chi);

#define _hop_z_m_pre32()			\
  _vector_i_sub(phi32[ix]->s0, rs.s0, rs.s2);	\
  _vector_i_add(phi32[ix]->s1, rs.s1, rs.s3);

#define _hop_t_p_post32();			\
  _vector_assign(rs.s0, phi32[ix]->s0);		\
  _vector_assign(rs.s2, phi32[ix]->s0);		\
  _vector_assign(rs.s1, phi32[ix]->s1);		\
  _vector_assign(rs.s3, phi32[ix]->s1);		\

#define _hop_t_m_post32();			\
  _vector_assign(psi, phi32[ix]->s0);		\
  _su3_inverse_multiply(chi,(*U), psi);		\
  _complexcjg_times_vector(psi,ka0,chi);	\
  _vector_add_assign(rs.s0, psi);		\
  _vector_sub_assign(rs.s2, psi);		\
  _vector_assign(psi, phi32[ix]->s1);		\
  _su3_inverse_multiply(chi,(*U), psi);		\
  _complexcjg_times_vector(psi,ka0,chi);	\
  _vector_add_assign(rs.s1, psi);		\
  _vector_sub_assign(rs.s3, psi);

#define _hop_x_p_post32();				\
  _vector_add_assign(rs.s0, phi32[ix]->s0);		\
  _vector_i_sub_assign(rs.s3, phi32[ix]->s0);		\
  _vector_add_assign(rs.s1, phi32[ix]->s1);		\
  _vector_i_sub_assign(rs.s2, phi32[ix]->s1);

#define _hop_x_m_post32();			\
  _vector_assign(psi, phi32[ix]->s0);		\
  _su3_inverse_multiply(chi,(*U), psi);		\
  _complexcjg_times_vector(psi,ka1,chi);	\
  _vector_add_assign(rs.s0, psi);		\
  _vector_i_add_assign(rs.s3, psi);		\
  _vector_assign(psi, phi32[ix]->s1);		\
  _su3_inverse_multiply(chi,(*U), psi);		\
  _complexcjg_times_vector(psi,ka1,chi);	\
  _vector_add_assign(rs.s1, psi);		\
  _vector_i_add_assign(rs.s2, psi);

#define _hop_y_p_post32();			\
  _vector_add_assign(rs.s0, phi32[ix]->s0);	\
  _vector_add_assign(rs.s3, phi32[ix]->s0);	\
  _vector_add_assign(rs.s1, phi32[ix]->s1);	\
  _vector_sub_assign(rs.s2, phi32[ix]->s1);

#define _hop_y_m_post32();			\
  _vector_assign(psi, phi32[ix]->s0);		\
  _su3_inverse_multiply(chi,(*U), psi);		\
  _complexcjg_times_vector(psi,ka2,chi);	\
  _vector_add_assign(rs.s0, psi);		\
  _vector_sub_assign(rs.s3, psi);		\
  _vector_assign(psi, phi32[ix]->s1);		\
  _su3_inverse_multiply(chi, (*U), psi);	\
  _complexcjg_times_vector(psi,ka2,chi);	\
  _vector_add_assign(rs.s1, psi);		\
  _vector_add_assign(rs.s2, psi);

#define _hop_z_p_post32();			\
  _vector_add_assign(rs.s0, phi32[ix]->s0);	\
  _vector_i_sub_assign(rs.s2, phi32[ix]->s0);	\
  _vector_add_assign(rs.s1, phi32[ix]->s1);	\
  _vector_i_add_assign(rs.s3, phi32[ix]->s1);

#define _hop_z_m_post32();			\
  _vector_assign(psi, phi32[ix]->s0);		\
  _su3_inverse_multiply(chi,(*U), psi);		\
  _complexcjg_times_vector(psi,ka3,chi);	\
  _vector_add_assign(rs.s0, psi);		\
  _vector_i_add_assign(rs.s2, psi);		\
  _vector_assign(psi, phi32[ix]->s1);		\
  _su3_inverse_multiply(chi,(*U), psi);		\
  _complexcjg_times_vector(psi,ka3,chi);	\
  _vector_add_assign(rs.s1, psi);		\
  _vector_i_sub_assign(rs.s3, psi);

#define _hop_t_p_pre()					\
  _vector_assign(rs.s0, s->s0);				\
  _vector_assign(rs.s1, s->s1);				\
  _vector_assign(rs.s2, s->s2);				\
  _vector_assign(rs.s3, s->s3);				\
  _vector_add(psi, rs.s0, rs.s2);			\
  _vector_add(psi2, rs.s1, rs.s3);			\
  _su3_multiply(chi,(*U),psi);				\
  _su3_multiply(chi2,(*U),psi2);			\
  _complex_times_vector(phi[ix]->s0, ka0, chi);		\
  _complex_times_vector(phi[ix]->s1, ka0, chi2);

#define _hop_t_m_pre()				\
  _vector_sub(phi[ix]->s0, rs.s0, rs.s2);	\
  _vector_sub(phi[ix]->s1, rs.s1, rs.s3);

#define _hop_x_p_pre()					\
  _vector_i_add(psi, rs.s0, rs.s3);			\
  _vector_i_add(psi2, rs.s1, rs.s2);			\
  _su3_multiply(chi, (*U), psi);			\
  _su3_multiply(chi2, (*U), psi2);			\
  _complex_times_vector(phi[ix]->s0, ka1, chi);		\
  _complex_times_vector(phi[ix]->s1, ka1, chi2);

#define _hop_x_m_pre()					\
  _vector_i_sub(phi[ix]->s0, rs.s0, rs.s3);		\
  _vector_i_sub(phi[ix]->s1, rs.s1, rs.s2);

#define _hop_y_p_pre()					\
  _vector_add(psi, rs.s0, rs.s3);			\
  _vector_sub(psi2, rs.s1, rs.s2);			\
  _su3_multiply(chi,(*U),psi);				\
  _su3_multiply(chi2,(*U),psi2);			\
  _complex_times_vector(phi[ix]->s0, ka2, chi);		\
  _complex_times_vector(phi[ix]->s1, ka2, chi2);

#define _hop_y_m_pre()					\
  _vector_sub(phi[ix]->s0, rs.s0, rs.s3);		\
  _vector_add(phi[ix]->s1, rs.s1, rs.s2);

#define _hop_z_p_pre()					\
  _vector_i_add(psi, rs.s0, rs.s2);			\
  _vector_i_sub(psi2, rs.s1, rs.s3);			\
  _su3_multiply(chi, (*U), psi);			\
  _su3_multiply(chi2,(*U),psi2);			\
  _complex_times_vector(phi[ix]->s0, ka3, chi);		\
  _complex_times_vector(phi[ix]->s1, ka3, chi2);

#define _hop_z_m_pre()					\
  _vector_i_sub(phi[ix]->s0, rs.s0, rs.s2);		\
  _vector_i_add(phi[ix]->s1, rs.s1, rs.s3);

#define _hop_t_p_post()					\
  _vector_assign(rs.s0, phi[ix]->s0);			\
  _vector_assign(rs.s2, phi[ix]->s0);			\
  _vector_assign(rs.s1, phi[ix]->s1);			\
  _vector_assign(rs.s3, phi[ix]->s1);

#define _hop_t_m_post()					\
  _su3_inverse_multiply(chi,(*U),phi[ix]->s0);		\
  _su3_inverse_multiply(chi2,(*U),phi[ix]->s1);		\
  _complexcjg_times_vector(psi,ka0,chi);		\
  _complexcjg_times_vector(psi2,ka0,chi2);		\
  _vector_add_assign(rs.s0, psi);			\
  _vector_sub_assign(rs.s2, psi);			\
  _vector_add_assign(rs.s1, psi2);			\
  _vector_sub_assign(rs.s3, psi2);

#define _hop_x_p_post()					\
  _vector_add_assign(rs.s0, phi[ix]->s0);		\
  _vector_i_sub_assign(rs.s3, phi[ix]->s0);		\
  _vector_add_assign(rs.s1, phi[ix]->s1);		\
  _vector_i_sub_assign(rs.s2, phi[ix]->s1);

#define _hop_x_m_post()					\
  _su3_inverse_multiply(chi,(*U), phi[ix]->s0);		\
  _su3_inverse_multiply(chi2, (*U), phi[ix]->s1);	\
  _complexcjg_times_vector(psi,ka1,chi);		\
  _complexcjg_times_vector(psi2,ka1,chi2);		\
  _vector_add_assign(rs.s0, psi);			\
  _vector_i_add_assign(rs.s3, psi);			\
  _vector_add_assign(rs.s1, psi2);			\
  _vector_i_add_assign(rs.s2, psi2);

#define _hop_y_p_post()					\
  _vector_add_assign(rs.s0, phi[ix]->s0);		\
  _vector_add_assign(rs.s3, phi[ix]->s0);		\
  _vector_add_assign(rs.s1, phi[ix]->s1);		\
  _vector_sub_assign(rs.s2, phi[ix]->s1);

#define _hop_y_m_post()					\
  _su3_inverse_multiply(chi,(*U), phi[ix]->s0);		\
  _su3_inverse_multiply(chi2, (*U), phi[ix]->s1);	\
  _complexcjg_times_vector(psi,ka2,chi);		\
  _complexcjg_times_vector(psi2,ka2,chi2);		\
  _vector_add_assign(rs.s0, psi);			\
  _vector_sub_assign(rs.s3, psi);			\
  _vector_add_assign(rs.s1, psi2);			\
  _vector_add_assign(rs.s2, psi2);

#define _hop_z_p_post()					\
  _vector_add_assign(rs.s0, phi[ix]->s0);		\
  _vector_i_sub_assign(rs.s2, phi[ix]->s0);		\
  _vector_add_assign(rs.s1, phi[ix]->s1);		\
  _vector_i_add_assign(rs.s3, phi[ix]->s1);

#define _hop_z_m_post()					\
  _su3_inverse_multiply(chi,(*U), phi[ix]->s0);		\
  _su3_inverse_multiply(chi2, (*U), phi[ix]->s1);	\
  _complexcjg_times_vector(psi,ka3,chi);		\
  _complexcjg_times_vector(psi2,ka3,chi2);		\
  _vector_add_assign(rs.s0, psi);			\
  _vector_add_assign(rs.s1, psi2);			\
  _vector_i_add_assign(rs.s2, psi);			\
  _vector_i_sub_assign(rs.s3, psi2);

#define _hop_mul_g5_cmplx_and_store(res)			\
  _complex_times_vector((res)->s0, cfactor, rs.s0);		\
  _complex_times_vector((res)->s1, cfactor, rs.s1);		\
  _complexcjg_times_vector((res)->s2, cfactor, rs.s2);	\
  _complexcjg_times_vector((res)->s3, cfactor, rs.s3);

#define _g5_cmplx_sub_hop_and_g5store(res)		\
  _complex_times_vector(psi, cfactor, pn->s0);		\
  _vector_sub((res)->s0, psi, rs.s0);			\
  _complex_times_vector(psi2, cfactor, pn->s1);		\
  _vector_sub((res)->s1, psi2, rs.s1);			\
  _complexcjg_times_vector(psi, cfactor, pn->s2);	\
  _vector_sub((res)->s2, rs.s2, psi);			\
  _complexcjg_times_vector(psi2, cfactor, pn->s3);	\
  _vector_sub((res)->s3, rs.s3, psi2);


#define _hop_store_post(res)		\
  _vector_assign(res->s0, rs.s0);	\
  _vector_assign(res->s1, rs.s1);	\
  _vector_assign(res->s2, rs.s2);	\
  _vector_assign(res->s3, rs.s3);


#define _declare_hregs()				\
  spinor ALIGN rs;					\
  su3_vector ALIGN psi, chi, psi2, chi2;

#endif

#endif

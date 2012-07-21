/**********************************************************************
 *
 * Copyright (C) 2012 Carsten Urbach
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
 * These are routines and macros for using half of the available
 * four floating point units of the BG/Q processor
 *
 **********************************************************************/

#ifndef _BGQ2_H
#define _BGQ2_H

//#define regtype vector4double

#define _vec_load2(r0, r1, r2, phi)	  \
  (r0) = vec_ld2(0L, (double*) &(phi).c0); \
  (r1) = vec_ld2(0L, (double*) &(phi).c1); \
  (r2) = vec_ld2(0L, (double*) &(phi).c2);

#define _vec_load2_32(r0, r1, r2, phi)	 \
  (r0) = vec_ld2(0L, (float*) &(phi).c0); \
  (r1) = vec_ld2(0L, (float*) &(phi).c1); \
  (r2) = vec_ld2(0L, (float*) &(phi).c2);

#define _vec_load2c(r0, r1, phi)	 \
  r0 = vec_ld(0L, (double*) &(phi).c0);	 \
  r1 = vec_ld2(0L, (double*) &(phi).c2); 


inline void vec_load2(vector4double * r, su3_vector * const phi) {
#pragma disjoint(*r,*phi)
  r[0] = vec_ld2(0, (double*) &phi->c0); 
  r[1] = vec_ld2(0, (double*) &phi->c1);
  r[2] = vec_ld2(0, (double*) &phi->c2);
  return;
}

#define _vec_store2(phi, r0, r1, r2)		\
  vec_st2((r0), 0, (double*) &phi.c0);		\
  vec_st2((r1), 0, (double*) &phi.c1);		\
  vec_st2((r2), 0, (double*) &phi.c2);

#define _vec_store2_32(phi, r0, r1, r2)		\
  vec_st2((r0), 0, (float*) &phi.c0);		\
  vec_st2((r1), 0, (float*) &phi.c1);		\
  vec_st2((r2), 0, (float*) &phi.c2);

inline void vec_store2(su3_vector * const phi, vector4double * r) {
#pragma disjoint(*r,*phi)
  vec_st2(r[0], 0, (double*) &phi->c0);
  vec_st2(r[1], 0, (double*) &phi->c1);
  vec_st2(r[2], 0, (double*) &phi->c2);
  return;
}

// r = r + s
#define _vec_add2(r0, r1, r2, s0, s1, s2) \
  (r0) = vec_add((r0), (s0));		  \
  (r1) = vec_add((r1), (s1));		  \
  (r2) = vec_add((r2), (s2));

#define _vec_add_to2(rs0, rs1, rs2, r0, r1, r2, s0, s1, s2)	\
  (rs0) = vec_add((r0), (s0));		  \
  (rs1) = vec_add((r1), (s1));		  \
  (rs2) = vec_add((r2), (s2));

inline void vec_add2(vector4double * restrict r, vector4double * restrict s) {
#pragma disjoint(*s, *r)
  r[0] = vec_add(r[0], s[0]);
  r[1] = vec_add(r[1], s[1]);
  r[2] = vec_add(r[2], s[2]);
  return;
}

// r = r + s
#define _vec_add_double2(r0, r1, r2, r3, r4, r5, s0, s1, s2, s3, s4, s5) \
  (r0) = vec_add((r0), (s0));						\
  (r1) = vec_add((r1), (s1));						\
  (r2) = vec_add((r2), (s2));						\
  (r3) = vec_add((r3), (s3));						\
  (r4) = vec_add((r4), (s4));						\
  (r5) = vec_add((r5), (s5));						

#define _vec_add_double_to2(rs0, rs1, rs2, rs3, rs4, rs5, r0, r1, r2, r3, r4, r5, s0, s1, s2, s3, s4, s5) \
  (rs0) = vec_add((r0), (s0));						\
  (rs1) = vec_add((r1), (s1));						\
  (rs2) = vec_add((r2), (s2));						\
  (rs3) = vec_add((r3), (s3));						\
  (rs4) = vec_add((r4), (s4));						\
  (rs5) = vec_add((r5), (s5));						

inline void vec_add_double2(vector4double * restrict r, vector4double * restrict s) {
#pragma disjoint(*s, *r)
  r[0] = vec_add(r[0], s[0]);
  r[1] = vec_add(r[1], s[1]);
  r[2] = vec_add(r[2], s[2]);
  r[3] = vec_add(r[3], s[3]);
  r[4] = vec_add(r[4], s[4]);
  r[5] = vec_add(r[5], s[5]);

  return;
}

// r = r - s
#define _vec_sub2(r0, r1, r2, s0, s1, s2) \
  (r0) = vec_sub((r0), (s0));		  \
  (r1) = vec_sub((r1), (s1));		  \
  (r2) = vec_sub((r2), (s2));

#define _vec_sub_to2(rs0, rs1, rs2, r0, r1, r2, s0, s1, s2)	\
  (rs0) = vec_sub((r0), (s0));		  \
  (rs1) = vec_sub((r1), (s1));		  \
  (rs2) = vec_sub((r2), (s2));


inline void vec_sub2(vector4double * restrict r, vector4double * restrict s) {
#pragma disjoint(*s, *r)
  r[0] = vec_sub(r[0], s[0]);
  r[1] = vec_sub(r[1], s[1]);
  r[2] = vec_sub(r[2], s[2]);
  return;
}

// r = r - s
#define _vec_sub_double2(r0, r1, r2, r3, r4, r5, s0, s1, s2, s3, s4, s5) \
  (r0) = vec_sub((r0), (s0));						\
  (r1) = vec_sub((r1), (s1));						\
  (r2) = vec_sub((r2), (s2));						\
  (r3) = vec_sub((r3), (s3));						\
  (r4) = vec_sub((r4), (s4));						\
  (r5) = vec_sub((r5), (s5));						

#define _vec_sub_to_double2(rs0, rs1, rs2, rs3, rs4, rs5, r0, r1, r2, r3, r4, r5, s0, s1, s2, s3, s4, s5) \
  (rs0) = vec_sub((r0), (s0));						\
  (rs1) = vec_sub((r1), (s1));						\
  (rs2) = vec_sub((r2), (s2));						\
  (rs3) = vec_sub((r3), (s3));						\
  (rs4) = vec_sub((r4), (s4));						\
  (rs5) = vec_sub((r5), (s5));						


inline void vec_sub_double2(vector4double * restrict r, vector4double * restrict s) {
#pragma disjoint(*s, *r)
  r[0] = vec_sub(r[0], s[0]);
  r[1] = vec_sub(r[1], s[1]);
  r[2] = vec_sub(r[2], s[2]);
  r[3] = vec_sub(r[3], s[3]);
  r[4] = vec_sub(r[4], s[4]);
  r[5] = vec_sub(r[5], s[5]);

  return;
}

// r = r + i*s
#define _vec_i_mul_add2(r0, r1, r2, s0, s1, s2, tmp) \
  tmp = vec_splats(1.);				     \
  r0 = vec_xxnpmadd(s0, tmp, r0);		     \
  r1 = vec_xxnpmadd(s1, tmp, r1);		     \
  r2 = vec_xxnpmadd(s2, tmp, r2);

#define _vec_i_mul_add_to2(rs0, rs1, rs2, r0, r1, r2, s0, s1, s2, tmp)	\
  tmp = vec_splats(1.);				     \
  rs0 = vec_xxnpmadd(s0, tmp, r0);		     \
  rs1 = vec_xxnpmadd(s1, tmp, r1);		     \
  rs2 = vec_xxnpmadd(s2, tmp, r2);

// tmp, temporary memory
inline void vec_i_mul_add2(vector4double * restrict r, vector4double * restrict s, 
			   vector4double * restrict tmp) {
#pragma disjoint(*s, *r)
#pragma disjoint(*s, *tmp)
#pragma disjoint(*r, *tmp)
  tmp[0] = vec_splats(1.);
  r[0] = vec_xxnpmadd(s[0], tmp[0], r[0]);
  r[1] = vec_xxnpmadd(s[1], tmp[0], r[1]);
  r[2] = vec_xxnpmadd(s[2], tmp[0], r[2]);
  return;
}

// r = r + i*s
#define _vec_i_mul_add_double2(r0, r1, r2, r3, r4, r5, s0, s1, s2, s3, s4, s5, tmp) \
  tmp = vec_splats(1.);				     \
  r0 = vec_xxnpmadd(s0, tmp, r0);		     \
  r1 = vec_xxnpmadd(s1, tmp, r1);		     \
  r2 = vec_xxnpmadd(s2, tmp, r2);		     \
  r3 = vec_xxnpmadd(s3, tmp, r3);		     \
  r4 = vec_xxnpmadd(s4, tmp, r4);		     \
  r5 = vec_xxnpmadd(s5, tmp, r5);

#define _vec_i_mul_add_double_to2(rs0, rs1, rs2, rs3, rs4, rs5, r0, r1, r2, r3, r4, r5, s0, s1, s2, s3, s4, s5, tmp) \
  tmp = vec_splats(1.);				     \
  rs0 = vec_xxnpmadd(s0, tmp, r0);		     \
  rs1 = vec_xxnpmadd(s1, tmp, r1);		     \
  rs2 = vec_xxnpmadd(s2, tmp, r2);		     \
  rs3 = vec_xxnpmadd(s3, tmp, r3);		     \
  rs4 = vec_xxnpmadd(s4, tmp, r4);		     \
  rs5 = vec_xxnpmadd(s5, tmp, r5);


// tmp, temporary memory
inline void vec_i_mul_add_double2(vector4double * restrict r, vector4double * restrict s, 
				  vector4double * restrict tmp) {
#pragma disjoint(*s, *r)
#pragma disjoint(*s, *tmp)
#pragma disjoint(*r, *tmp)
  tmp[0] = vec_splats(1.);
  r[0] = vec_xxnpmadd(s[0], tmp[0], r[0]);
  r[1] = vec_xxnpmadd(s[1], tmp[0], r[1]);
  r[2] = vec_xxnpmadd(s[2], tmp[0], r[2]);
  r[3] = vec_xxnpmadd(s[3], tmp[0], r[3]);
  r[4] = vec_xxnpmadd(s[4], tmp[0], r[4]);
  r[5] = vec_xxnpmadd(s[5], tmp[0], r[5]);
  return;
}

// r = r - i*s
#define _vec_i_mul_sub2(r0, r1, r2, s0, s1, s2, tmp) \
  tmp = vec_splats(-1.);				     \
  r0 = vec_xxnpmadd(s0, tmp, r0);		     \
  r1 = vec_xxnpmadd(s1, tmp, r1);		     \
  r2 = vec_xxnpmadd(s2, tmp, r2);

#define _vec_i_mul_sub_to2(rs0, rs1, rs2, r0, r1, r2, s0, s1, s2, tmp)	\
  tmp = vec_splats(-1.);				     \
  rs0 = vec_xxnpmadd(s0, tmp, r0);		     \
  rs1 = vec_xxnpmadd(s1, tmp, r1);		     \
  rs2 = vec_xxnpmadd(s2, tmp, r2);


// tmp, temporary memory
inline void vec_i_mul_sub2(vector4double * restrict r, vector4double * restrict s, 
			   vector4double * restrict tmp) {
#pragma disjoint(*s, *r)
#pragma disjoint(*s, *tmp)
#pragma disjoint(*r, *tmp)
  tmp[0] = vec_splats(-1.);
  r[0] = vec_xxnpmadd(s[0], tmp[0], r[0]);
  r[1] = vec_xxnpmadd(s[1], tmp[0], r[1]);
  r[2] = vec_xxnpmadd(s[2], tmp[0], r[2]);
  return;
}

// r = r - i*s
#define _vec_i_mul_sub_double2(r0, r1, r2, r3, r4, r5, s0, s1, s2, s3, s4, s5, tmp) \
  tmp = vec_splats(-1.);				     \
  r0 = vec_xxnpmadd(s0, tmp, r0);		     \
  r1 = vec_xxnpmadd(s1, tmp, r1);		     \
  r2 = vec_xxnpmadd(s2, tmp, r2);		     \
  r3 = vec_xxnpmadd(s3, tmp, r3);		     \
  r4 = vec_xxnpmadd(s4, tmp, r4);		     \
  r5 = vec_xxnpmadd(s5, tmp, r5);

#define _vec_i_mul_sub_double_to2(rs0, rs1, rs2, rs3, rs4, rs5, r0, r1, r2, r3, r4, r5, s0, s1, s2, s3, s4, s5, tmp) \
  tmp = vec_splats(-1.);				     \
  rs0 = vec_xxnpmadd(s0, tmp, r0);		     \
  rs1 = vec_xxnpmadd(s1, tmp, r1);		     \
  rs2 = vec_xxnpmadd(s2, tmp, r2);		     \
  rs3 = vec_xxnpmadd(s3, tmp, r3);		     \
  rs4 = vec_xxnpmadd(s4, tmp, r4);		     \
  rs5 = vec_xxnpmadd(s5, tmp, r5);


// tmp, temporary memory
inline void vec_i_mul_sub_double2(vector4double * restrict r, vector4double * restrict s, 
				  vector4double * restrict tmp) {
#pragma disjoint(*s, *r)
#pragma disjoint(*s, *tmp)
#pragma disjoint(*r, *tmp)
  tmp[0] = vec_splats(-1.);
  r[0] = vec_xxnpmadd(s[0], tmp[0], r[0]);
  r[1] = vec_xxnpmadd(s[1], tmp[0], r[1]);
  r[2] = vec_xxnpmadd(s[2], tmp[0], r[2]);
  r[3] = vec_xxnpmadd(s[3], tmp[0], r[3]);
  r[4] = vec_xxnpmadd(s[4], tmp[0], r[4]);
  r[5] = vec_xxnpmadd(s[5], tmp[0], r[5]);
  return;
}

#define _vec_cmplx_mul_double2(rs0, rs1, rs2, rs3, rs4, rs5, r0, r1, r2, r3, r4, r5, rtmp) \
  rs0 = vec_xmul(r0, rtmp);						\
  rs1 = vec_xmul(r1, rtmp);						\
  rs2 = vec_xmul(r2, rtmp);						\
  rs3 = vec_xmul(r3, rtmp);						\
  rs4 = vec_xmul(r4, rtmp);						\
  rs5 = vec_xmul(r5, rtmp);						\
  rs0 = vec_xxnpmadd(rtmp, r0, rs0);					\
  rs1 = vec_xxnpmadd(rtmp, r1, rs1);					\
  rs2 = vec_xxnpmadd(rtmp, r2, rs2);					\
  rs3 = vec_xxnpmadd(rtmp, r3, rs3);					\
  rs4 = vec_xxnpmadd(rtmp, r4, rs4);					\
  rs5 = vec_xxnpmadd(rtmp, r5, rs5);

#define _vec_cmplx_mul_double2c(rs0, rs1, rs2, r0, r1, r2, rtmp) \
  rs0 = vec_xmul(r0, rtmp);						\
  rs1 = vec_xmul(r1, rtmp);						\
  rs2 = vec_xmul(r2, rtmp);						\
  rs0 = vec_xxnpmadd(rtmp, r0, rs0);					\
  rs1 = vec_xxnpmadd(rtmp, r1, rs1);					\
  rs2 = vec_xxnpmadd(rtmp, r2, rs2);


inline void vec_cmplx_mul_double2(vector4double * restrict rs, vector4double * restrict r,
				  vector4double * tmp, complex double * c) {
#pragma disjoint(*rs, *r)
#pragma disjoint(*r, *tmp)
#pragma disjoint(*tmp, rs)
  __alignx(32, c);
  tmp[0] = vec_ld2(0, (double*) c);
  rs[0] = vec_xmul(r[0], tmp[0]);
  rs[1] = vec_xmul(r[1], tmp[0]);
  rs[2] = vec_xmul(r[2], tmp[0]);
  rs[3] = vec_xmul(r[3], tmp[0]);
  rs[4] = vec_xmul(r[4], tmp[0]);
  rs[5] = vec_xmul(r[5], tmp[0]);
  rs[0] = vec_xxnpmadd(tmp[0], r[0], rs[0]);
  rs[1] = vec_xxnpmadd(tmp[0], r[1], rs[1]);
  rs[2] = vec_xxnpmadd(tmp[0], r[2], rs[2]);
  rs[3] = vec_xxnpmadd(tmp[0], r[3], rs[3]);
  rs[4] = vec_xxnpmadd(tmp[0], r[4], rs[4]);
  rs[5] = vec_xxnpmadd(tmp[0], r[5], rs[5]);
  return;
}

#define _vec_cmplxcg_mul_double2(rs0, rs1, rs2, rs3, rs4, rs5, r0, r1, r2, r3, r4, r5, tmp) \
  rs0 = vec_xmul(rtmp, r0);						\
  rs1 = vec_xmul(rtmp, r1);						\
  rs2 = vec_xmul(rtmp, r2);						\
  rs3 = vec_xmul(rtmp, r3);						\
  rs4 = vec_xmul(rtmp, r4);						\
  rs5 = vec_xmul(rtmp, r5);						\
  rs0 = vec_xxcpnmadd(r0, rtmp, rs0);					\
  rs1 = vec_xxcpnmadd(r1, rtmp, rs1);					\
  rs2 = vec_xxcpnmadd(r2, rtmp, rs2);					\
  rs3 = vec_xxcpnmadd(r3, rtmp, rs3);					\
  rs4 = vec_xxcpnmadd(r4, rtmp, rs4);					\
  rs5 = vec_xxcpnmadd(r5, rtmp, rs5);					\

#define _vec_cmplxcg_mul_double2c(rs0, rs1, rs2, r0, r1, r2, tmp) \
  rs0 = vec_xmul(rtmp, r0);						\
  rs1 = vec_xmul(rtmp, r1);						\
  rs2 = vec_xmul(rtmp, r2);						\
  rs0 = vec_xxcpnmadd(r0, rtmp, rs0);					\
  rs1 = vec_xxcpnmadd(r1, rtmp, rs1);					\
  rs2 = vec_xxcpnmadd(r2, rtmp, rs2);

#define _vec_unfuse(r0, r1, r2, r3, r4, r5)	\
  r3 = vec_sldw(r0, r0, 2);			\
  r4 = vec_sldw(r1, r1, 2);			\
  r5 = vec_sldw(r2, r2, 2);

  
  
inline void vec_cmplxcg_mul_double2(vector4double * restrict rs, vector4double * restrict r,
				    vector4double * tmp, complex double * c) {
#pragma disjoint(*rs, *r)
#pragma disjoint(*r, *tmp)
#pragma disjoint(*tmp, rs)
  __alignx(32, c);
  tmp[0] = vec_ld2(0, (double*) c);
  rs[0] = vec_xmul(tmp[0], r[0]);
  rs[1] = vec_xmul(tmp[0], r[1]);
  rs[2] = vec_xmul(tmp[0], r[2]);
  rs[3] = vec_xmul(tmp[0], r[3]);
  rs[4] = vec_xmul(tmp[0], r[4]);
  rs[5] = vec_xmul(tmp[0], r[5]);
  rs[0] = vec_xxcpnmadd(r[0], tmp[0], rs[0]);
  rs[1] = vec_xxcpnmadd(r[1], tmp[0], rs[1]);
  rs[2] = vec_xxcpnmadd(r[2], tmp[0], rs[2]);
  rs[3] = vec_xxcpnmadd(r[3], tmp[0], rs[3]);
  rs[4] = vec_xxcpnmadd(r[4], tmp[0], rs[4]);
  rs[5] = vec_xxcpnmadd(r[5], tmp[0], rs[5]);
  return;
}


// multiplies one su3 matrix with two su3_vectors
// the first of which stored in r[0-2]
// and the second one in r[3-5]
//
// the resulting two vectors are stored in
// r[6-11]
//
// this routine uses only half of the 4 doubles in vector4double
#define _vec_su3_multiply_double2b(u)		\
  U[0] = vec_ld2(0, (double*) &(u)->c00);	\
  U[3] = vec_ld2(0, (double*) &(u)->c01);	\
  U[6] = vec_ld2(0, (double*) &(u)->c02);	\
  U[1] = vec_ld2(0, (double*) &(u)->c10);	\
  U[4] = vec_ld2(0, (double*) &(u)->c11);	\
  U[7] = vec_ld2(0, (double*) &(u)->c12);	\
  U[2] = vec_ld2(0, (double*) &(u)->c20);	\
  r[6] = vec_xmul(r[0], U[0]);			\
  r[7] = vec_xmul(r[0], U[1]);			\
  r[8] = vec_xmul(r[0], U[2]);			\
  r[9] = vec_xmul(r[3], U[0]);			\
  r[10] = vec_xmul(r[3], U[1]);			\
  r[11] = vec_xmul(r[3], U[2]);			\
						\
  r[6] = vec_xxnpmadd(U[0], r[0], r[6]);	\
  r[7] = vec_xxnpmadd(U[1], r[0], r[7]);	\
  r[8] = vec_xxnpmadd(U[2], r[0], r[8]);	\
  r[9] = vec_xxnpmadd(U[0], r[3], r[9]);	\
  r[10] = vec_xxnpmadd(U[1], r[3], r[10]);	\
  r[11] = vec_xxnpmadd(U[2], r[3], r[11]);	\
  U[5] = vec_ld2(0, (double*) &(u)->c21);	\
  						\
  r[6] = vec_xmadd(r[1], U[3], r[6]);		\
  r[7] = vec_xmadd(r[1], U[4], r[7]);		\
  r[8] = vec_xmadd(r[1], U[5], r[8]);		\
  r[9] = vec_xmadd(r[4], U[3], r[9]);		\
  r[10] = vec_xmadd(r[4], U[4], r[10]);		\
  r[11] = vec_xmadd(r[4], U[5], r[11]);		\
  						\
  r[6] = vec_xxnpmadd(U[3], r[1], r[6]);	\
  r[7] = vec_xxnpmadd(U[4], r[1], r[7]);	\
  r[8] = vec_xxnpmadd(U[5], r[1], r[8]);	\
  r[9] = vec_xxnpmadd(U[3], r[4], r[9]);	\
  r[10] = vec_xxnpmadd(U[4], r[4], r[10]);	\
  r[11] = vec_xxnpmadd(U[5], r[4], r[11]);	\
  U[8] = vec_ld2(0, (double*) &(u)->c22);	\
						\
  r[6] = vec_xmadd(r[2], U[6], r[6]);		\
  r[7] = vec_xmadd(r[2], U[7], r[7]);		\
  r[8] = vec_xmadd(r[2], U[8], r[8]);		\
  r[9] = vec_xmadd(r[5], U[6], r[9]);		\
  r[10] = vec_xmadd(r[5], U[7], r[10]);		\
  r[11] = vec_xmadd(r[5], U[8], r[11]);		\
						\
  r[6] = vec_xxnpmadd(U[6], r[2], r[6]);	\
  r[7] = vec_xxnpmadd(U[7], r[2], r[7]);	\
  r[8] = vec_xxnpmadd(U[8], r[2], r[8]);	\
  r[9] = vec_xxnpmadd(U[6], r[5], r[9]);	\
  r[10] = vec_xxnpmadd(U[7], r[5], r[10]);	\
  r[11] = vec_xxnpmadd(U[8], r[5], r[11]); 

#define _vec_su3_multiply_double2(u)		\
  U0 = vec_ld2(0, (double*) &(u)->c00);	\
  U3 = vec_ld2(0, (double*) &(u)->c01);	\
  U6 = vec_ld2(0, (double*) &(u)->c02);	\
  U1 = vec_ld2(0, (double*) &(u)->c10);	\
  U4 = vec_ld2(0, (double*) &(u)->c11);	\
  U7 = vec_ld2(0, (double*) &(u)->c12);	\
  U2 = vec_ld2(0, (double*) &(u)->c20);	\
  r6 = vec_xmul(r0, U0);			\
  r7 = vec_xmul(r0, U1);			\
  r8 = vec_xmul(r0, U2);			\
  r9 = vec_xmul(r3, U0);			\
  r10= vec_xmul(r3, U1);			\
  r11= vec_xmul(r3, U2);			\
						\
  r6 = vec_xxnpmadd(U0, r0, r6);	\
  r7 = vec_xxnpmadd(U1, r0, r7);	\
  r8 = vec_xxnpmadd(U2, r0, r8);	\
  r9 = vec_xxnpmadd(U0, r3, r9);	\
  r10= vec_xxnpmadd(U1, r3, r10);	\
  r11= vec_xxnpmadd(U2, r3, r11);	\
  U0 = vec_ld2(0, (double*) &(u)->c21);	\
  						\
  r6 = vec_xmadd(r1, U3, r6);		\
  r7 = vec_xmadd(r1, U4, r7);		\
  r8 = vec_xmadd(r1, U0, r8);		\
  r9 = vec_xmadd(r4, U3, r9);		\
  r10= vec_xmadd(r4, U4, r10);		\
  r11= vec_xmadd(r4, U0, r11);		\
       					\
  r6 = vec_xxnpmadd(U3, r1, r6);	\
  r7 = vec_xxnpmadd(U4, r1, r7);	\
  r8 = vec_xxnpmadd(U0, r1, r8);	\
  r9 = vec_xxnpmadd(U3, r4, r9);	\
  r10= vec_xxnpmadd(U4, r4, r10);	\
  r11= vec_xxnpmadd(U0, r4, r11);	\
  U1 = vec_ld2(0, (double*) &(u)->c22);	\
       					\
  r6 = vec_xmadd(r2, U6, r6);		\
  r7 = vec_xmadd(r2, U7, r7);		\
  r8 = vec_xmadd(r2, U1, r8);		\
  r9 = vec_xmadd(r5, U6, r9);		\
  r10= vec_xmadd(r5, U7, r10);		\
  r11= vec_xmadd(r5, U1, r11);		\
       					\
  r6 = vec_xxnpmadd(U6, r2, r6);	\
  r7 = vec_xxnpmadd(U7, r2, r7);	\
  r8 = vec_xxnpmadd(U1, r2, r8);	\
  r9 = vec_xxnpmadd(U6, r5, r9);	\
  r10= vec_xxnpmadd(U7, r5, r10);	\
  r11= vec_xxnpmadd(U1, r5, r11); 

#define _vec_su3_multiply_double2c(u)		\
  r8 = vec_gpci(00145);				\
  r9 = vec_gpci(02367);				\
  U0 = vec_ld2(0, (double*) &(u)->c00);		\
  U3 = vec_ld2(0, (double*) &(u)->c01);		\
  U6 = vec_ld2(0, (double*) &(u)->c02);		\
  U1 = vec_ld2(0, (double*) &(u)->c10);		\
  r7 = vec_perm(r0, r2, r8);			\
  U4 = vec_ld2(0, (double*) &(u)->c11);		\
  U7 = vec_ld2(0, (double*) &(u)->c12);		\
  U2 = vec_ld2(0, (double*) &(u)->c20);		\
  r4 = vec_xmul(r7, U0);			\
  r5 = vec_xmul(r7, U1);			\
  r6 = vec_xmul(r7, U2);			\
						\
  r4 = vec_xxnpmadd(U0, r7, r4);		\
  r5 = vec_xxnpmadd(U1, r7, r5);		\
  r6 = vec_xxnpmadd(U2, r7, r6);		\
  r7 = vec_perm(r0, r2, r9);			\
  U0 = vec_ld2(0, (double*) &(u)->c21);		\
						\
  r4 = vec_xmadd(r7, U3, r4);			\
  r5 = vec_xmadd(r7, U4, r5);			\
  r6 = vec_xmadd(r7, U0, r6);			\
  						\
  r4 = vec_xxnpmadd(U3, r7, r4);		\
  r5 = vec_xxnpmadd(U4, r7, r5);		\
  r6 = vec_xxnpmadd(U0, r7, r6);		\
  r7 = vec_perm(r1, r3, r8);			\
  U1 = vec_ld2(0, (double*) &(u)->c22);		\
						\
  r4 = vec_xmadd(r7, U6, r4);			\
  r5 = vec_xmadd(r7, U7, r5);			\
  r6 = vec_xmadd(r7, U1, r6);			\
  						\
  r4 = vec_xxnpmadd(U6, r7, r4);		\
  r5 = vec_xxnpmadd(U7, r7, r5);		\
  r6 = vec_xxnpmadd(U1, r7, r6);

#define _vec_su3_multiply_double2ct(u)		\
  r8 = vec_gpci(00167);				\
  U0 = vec_ld2(0, (double*) &(u)->c00);		\
  U3 = vec_ld2(0, (double*) &(u)->c01);		\
  U6 = vec_ld2(0, (double*) &(u)->c02);		\
  U1 = vec_ld2(0, (double*) &(u)->c10);		\
  r7 = vec_perm(r0, r1, r8);			\
  U4 = vec_ld2(0, (double*) &(u)->c11);		\
  U7 = vec_ld2(0, (double*) &(u)->c12);		\
  U2 = vec_ld2(0, (double*) &(u)->c20);		\
  r4 = vec_xmul(r7, U0);			\
  r5 = vec_xmul(r7, U1);			\
  r6 = vec_xmul(r7, U2);			\
						\
  r4 = vec_xxnpmadd(U0, r7, r4);		\
  r5 = vec_xxnpmadd(U1, r7, r5);		\
  r6 = vec_xxnpmadd(U2, r7, r6);		\
  r7 = vec_sldw(r0, r2, 2);			\
  U0 = vec_ld2(0, (double*) &(u)->c21);		\
						\
  r4 = vec_xmadd(r7, U3, r4);			\
  r5 = vec_xmadd(r7, U4, r5);			\
  r6 = vec_xmadd(r7, U0, r6);			\
  						\
  r4 = vec_xxnpmadd(U3, r7, r4);		\
  r5 = vec_xxnpmadd(U4, r7, r5);		\
  r6 = vec_xxnpmadd(U0, r7, r6);		\
  r7 = vec_perm(r1, r2, r8);			\
  U1 = vec_ld2(0, (double*) &(u)->c22);		\
						\
  r4 = vec_xmadd(r7, U6, r4);			\
  r5 = vec_xmadd(r7, U7, r5);			\
  r6 = vec_xmadd(r7, U1, r6);			\
  						\
  r4 = vec_xxnpmadd(U6, r7, r4);		\
  r5 = vec_xxnpmadd(U7, r7, r5);		\
  r6 = vec_xxnpmadd(U1, r7, r6);


inline void vec_su3_multiply_double2(su3 * const restrict u, vector4double * restrict U, 
				     vector4double * restrict r) {
#pragma disjoint(*U, *r)
  __alignx(32, u);
  __alignx(32, U);
  __alignx(32, r);

  U[0] = vec_ld2(0, (double*) &u->c00);
  U[3] = vec_ld2(0, (double*) &u->c01);
  U[6] = vec_ld2(0, (double*) &u->c02);
  U[1] = vec_ld2(0, (double*) &u->c10);
  U[4] = vec_ld2(0, (double*) &u->c11);
  U[7] = vec_ld2(0, (double*) &u->c12);
  U[2] = vec_ld2(0, (double*) &u->c20);
  r[6] = vec_xmul(r[0], U[0]);
  r[7] = vec_xmul(r[0], U[1]);
  r[8] = vec_xmul(r[0], U[2]);
  r[9] = vec_xmul(r[3], U[0]);
  r[10] = vec_xmul(r[3], U[1]);
  r[11] = vec_xmul(r[3], U[2]);

  r[6] = vec_xxnpmadd(U[0], r[0], r[6]);
  r[7] = vec_xxnpmadd(U[1], r[0], r[7]);
  r[8] = vec_xxnpmadd(U[2], r[0], r[8]);
  r[9] = vec_xxnpmadd(U[0], r[3], r[9]);
  r[10] = vec_xxnpmadd(U[1], r[3], r[10]);
  r[11] = vec_xxnpmadd(U[2], r[3], r[11]);
  U[5] = vec_ld2(0, (double*) &u->c21);

  r[6] = vec_xmadd(r[1], U[3], r[6]);
  r[7] = vec_xmadd(r[1], U[4], r[7]);
  r[8] = vec_xmadd(r[1], U[5], r[8]);
  r[9] = vec_xmadd(r[4], U[3], r[9]);
  r[10] = vec_xmadd(r[4], U[4], r[10]);
  r[11] = vec_xmadd(r[4], U[5], r[11]);

  r[6] = vec_xxnpmadd(U[3], r[1], r[6]);
  r[7] = vec_xxnpmadd(U[4], r[1], r[7]);
  r[8] = vec_xxnpmadd(U[5], r[1], r[8]);
  r[9] = vec_xxnpmadd(U[3], r[4], r[9]);
  r[10] = vec_xxnpmadd(U[4], r[4], r[10]);
  r[11] = vec_xxnpmadd(U[5], r[4], r[11]);
  U[8] = vec_ld2(0, (double*) &u->c22);

  r[6] = vec_xmadd(r[2], U[6], r[6]);
  r[7] = vec_xmadd(r[2], U[7], r[7]);
  r[8] = vec_xmadd(r[2], U[8], r[8]);
  r[9] = vec_xmadd(r[5], U[6], r[9]);
  r[10] = vec_xmadd(r[5], U[7], r[10]);
  r[11] = vec_xmadd(r[5], U[8], r[11]);

  r[6] = vec_xxnpmadd(U[6], r[2], r[6]);
  r[7] = vec_xxnpmadd(U[7], r[2], r[7]);
  r[8] = vec_xxnpmadd(U[8], r[2], r[8]);
  r[9] = vec_xxnpmadd(U[6], r[5], r[9]);
  r[10] = vec_xxnpmadd(U[7], r[5], r[10]);
  r[11] = vec_xxnpmadd(U[8], r[5], r[11]);
  return;
}

// multiplies the inverse of one su3 matrix with two su3_vectors
// the first of which stored in r[0-2]
// and the second one in r[3-5]
//
// the resulting two vectors are stored in
// r[6-11]
//
// this routine uses only half of the 4 doubles in vector4double
#define _vec_su3_inverse_multiply_double2(u)    \
  U0 = vec_ld2(0, (double*) &(u)->c00);		\
  U1 = vec_ld2(0, (double*) &(u)->c01);		\
  U2 = vec_ld2(0, (double*) &(u)->c02);		\
  						\
  r6 = vec_xmul(U0, r0);                        \
  r7 = vec_xmul(U1, r0);                        \
  r8 = vec_xmul(U2, r0);                        \
  r9 = vec_xmul(U0, r3);                        \
  r10= vec_xmul(U1, r3);                        \
  r11= vec_xmul(U2, r3);                        \
  						\
  r6 = vec_xxcpnmadd(r0, U0, r6);		\
  r7 = vec_xxcpnmadd(r0, U1, r7);		\
  r8 = vec_xxcpnmadd(r0, U2, r8);		\
  r9 = vec_xxcpnmadd(r3, U0, r9);		\
  r10= vec_xxcpnmadd(r3, U1, r10);		\
  r11= vec_xxcpnmadd(r3, U2, r11);		\
  						\
  U3 = vec_ld2(0, (double*) &(u)->c10);		\
  U4 = vec_ld2(0, (double*) &(u)->c11);		\
  U6 = vec_ld2(0, (double*) &(u)->c12);		\
  						\
  r6 = vec_xmadd(U3, r1, r6);			\
  r7 = vec_xmadd(U4, r1, r7);			\
  r8 = vec_xmadd(U6, r1, r8);			\
  r9 = vec_xmadd(U3, r4, r9);			\
  r10= vec_xmadd(U4, r4, r10);			\
  r11= vec_xmadd(U6, r4, r11);			\
  						\
  r6 = vec_xxcpnmadd(r1, U3, r6);		\
  r7 = vec_xxcpnmadd(r1, U4, r7);		\
  r8 = vec_xxcpnmadd(r1, U6, r8);		\
  r9 = vec_xxcpnmadd(r4, U3, r9);		\
  r10= vec_xxcpnmadd(r4, U4, r10);		\
  r11= vec_xxcpnmadd(r4, U6, r11);		\
  						\
  U0 = vec_ld2(0, (double*) &(u)->c20);		\
  U1 = vec_ld2(0, (double*) &(u)->c21);		\
  U2 = vec_ld2(0, (double*) &(u)->c22);		\
  						\
  r6 = vec_xmadd(U0, r2, r6);			\
  r7 = vec_xmadd(U1, r2, r7);			\
  r8 = vec_xmadd(U2, r2, r8);			\
  r9 = vec_xmadd(U0, r5, r9);			\
  r10= vec_xmadd(U1, r5, r10);			\
  r11= vec_xmadd(U2, r5, r11);			\
  						\
  r6 = vec_xxcpnmadd(r2, U0, r6);		\
  r7 = vec_xxcpnmadd(r2, U1, r7);		\
  r8 = vec_xxcpnmadd(r2, U2, r8);		\
  r9 = vec_xxcpnmadd(r5, U0, r9);		\
  r10= vec_xxcpnmadd(r5, U1, r10);		\
  r11= vec_xxcpnmadd(r5, U2, r11);


#define _vec_su3_inverse_multiply_double2c(u)	\
  U0 = vec_ld2(0, (double*) &(u)->c00);		\
  r8 = vec_gpci(00145);				\
  r9 = vec_gpci(02367);				\
  U1 = vec_ld2(0, (double*) &(u)->c01);		\
  r7 = vec_perm(r0, r2, r8);			\
  U2 = vec_ld2(0, (double*) &(u)->c02);		\
						\
  r4 = vec_xmul(U0, r7);			\
  r5 = vec_xmul(U1, r7);			\
  r6 = vec_xmul(U2, r7);			\
						\
  r4 = vec_xxcpnmadd(r7, U0, r4);		\
  r5 = vec_xxcpnmadd(r7, U1, r5);		\
  r6 = vec_xxcpnmadd(r7, U2, r6);		\
						\
  r7 = vec_perm(r0, r2, r9);			\
  U3 = vec_ld2(0, (double*) &(u)->c10);		\
  U4 = vec_ld2(0, (double*) &(u)->c11);		\
  U6 = vec_ld2(0, (double*) &(u)->c12);		\
  						\
  r4 = vec_xmadd(U3, r7, r4);			\
  r5 = vec_xmadd(U4, r7, r5);			\
  r6 = vec_xmadd(U6, r7, r6);			\
  						\
  r4 = vec_xxcpnmadd(r7, U3, r4);		\
  r5 = vec_xxcpnmadd(r7, U4, r5);		\
  r6 = vec_xxcpnmadd(r7, U6, r6);		\
						\
  r7 = vec_perm(r1, r3, r8);			\
  U0 = vec_ld2(0, (double*) &(u)->c20);		\
  U1 = vec_ld2(0, (double*) &(u)->c21);		\
  U2 = vec_ld2(0, (double*) &(u)->c22);		\
  						\
  r4 = vec_xmadd(U0, r7, r4);			\
  r5 = vec_xmadd(U1, r7, r5);			\
  r6 = vec_xmadd(U2, r7, r6);			\
  						\
  r4 = vec_xxcpnmadd(r7, U0, r4);		\
  r5 = vec_xxcpnmadd(r7, U1, r5);		\
  r6 = vec_xxcpnmadd(r7, U2, r6);

#define _vec_su3_inverse_multiply_double2ct(u)	\
  U0 = vec_ld2(0, (double*) &(u)->c00);		\
  r8 = vec_gpci(00167);				\
  U1 = vec_ld2(0, (double*) &(u)->c01);		\
  r7 = vec_perm(r0, r1, r8);			\
  U2 = vec_ld2(0, (double*) &(u)->c02);		\
						\
  r4 = vec_xmul(U0, r7);			\
  r5 = vec_xmul(U1, r7);			\
  r6 = vec_xmul(U2, r7);			\
						\
  r4 = vec_xxcpnmadd(r7, U0, r4);		\
  r5 = vec_xxcpnmadd(r7, U1, r5);		\
  r6 = vec_xxcpnmadd(r7, U2, r6);		\
						\
  r7 = vec_sldw(r0, r2, 2);			\
  U3 = vec_ld2(0, (double*) &(u)->c10);		\
  U4 = vec_ld2(0, (double*) &(u)->c11);		\
  U6 = vec_ld2(0, (double*) &(u)->c12);		\
  						\
  r4 = vec_xmadd(U3, r7, r4);			\
  r5 = vec_xmadd(U4, r7, r5);			\
  r6 = vec_xmadd(U6, r7, r6);			\
  						\
  r4 = vec_xxcpnmadd(r7, U3, r4);		\
  r5 = vec_xxcpnmadd(r7, U4, r5);		\
  r6 = vec_xxcpnmadd(r7, U6, r6);		\
						\
  r7 = vec_perm(r1, r2, r8);			\
  U0 = vec_ld2(0, (double*) &(u)->c20);		\
  U1 = vec_ld2(0, (double*) &(u)->c21);		\
  U2 = vec_ld2(0, (double*) &(u)->c22);		\
  						\
  r4 = vec_xmadd(U0, r7, r4);			\
  r5 = vec_xmadd(U1, r7, r5);			\
  r6 = vec_xmadd(U2, r7, r6);			\
  						\
  r4 = vec_xxcpnmadd(r7, U0, r4);		\
  r5 = vec_xxcpnmadd(r7, U1, r5);		\
  r6 = vec_xxcpnmadd(r7, U2, r6);


inline void vec_su3_inverse_multiply_double2(su3 * const restrict u, vector4double * restrict U, 
					     vector4double * restrict r) {
#pragma disjoint(*U, *r)
  __alignx(32, u);
  __alignx(32, U);
  __alignx(32, r);

  U[0] = vec_ld2(0, (double*) &u->c00);
  U[1] = vec_ld2(0, (double*) &u->c01);
  U[2] = vec_ld2(0, (double*) &u->c02);

  r[6] = vec_xmul(U[0], r[0]);
  r[7] = vec_xmul(U[1], r[0]);
  r[8] = vec_xmul(U[2], r[0]);
  r[9] = vec_xmul(U[0], r[3]);
  r[10] = vec_xmul(U[1], r[3]);
  r[11] = vec_xmul(U[2], r[3]);

  r[6] = vec_xxcpnmadd(r[0], U[0], r[6]);
  r[7] = vec_xxcpnmadd(r[0], U[1], r[7]);
  r[8] = vec_xxcpnmadd(r[0], U[2], r[8]);
  r[9] = vec_xxcpnmadd(r[3], U[0], r[9]);
  r[10] = vec_xxcpnmadd(r[3], U[1], r[10]);
  r[11] = vec_xxcpnmadd(r[3], U[2], r[11]);

  U[0] = vec_ld2(0, (double*) &u->c10);
  U[1] = vec_ld2(0, (double*) &u->c11);
  U[2] = vec_ld2(0, (double*) &u->c12);

  r[6] = vec_xmadd(U[0], r[1], r[6]);
  r[7] = vec_xmadd(U[1], r[1], r[7]);
  r[8] = vec_xmadd(U[2], r[1], r[8]);
  r[9] = vec_xmadd(U[0], r[4], r[9]);
  r[10] = vec_xmadd(U[1], r[4], r[10]);
  r[11] = vec_xmadd(U[2], r[4], r[11]);
  
  r[6] = vec_xxcpnmadd(r[1], U[0], r[6]);
  r[7] = vec_xxcpnmadd(r[1], U[1], r[7]);
  r[8] = vec_xxcpnmadd(r[1], U[2], r[8]);
  r[9] = vec_xxcpnmadd(r[4], U[0], r[9]);
  r[10] = vec_xxcpnmadd(r[4], U[1], r[10]);
  r[11] = vec_xxcpnmadd(r[4], U[2], r[11]);

  U[0] = vec_ld2(0, (double*) &u->c20);
  U[1] = vec_ld2(0, (double*) &u->c21);
  U[2] = vec_ld2(0, (double*) &u->c22);

  r[6] = vec_xmadd(U[0], r[2], r[6]);
  r[7] = vec_xmadd(U[1], r[2], r[7]);
  r[8] = vec_xmadd(U[2], r[2], r[8]);
  r[9] = vec_xmadd(U[0], r[5], r[9]);
  r[10] = vec_xmadd(U[1], r[5], r[10]);
  r[11] = vec_xmadd(U[2], r[5], r[11]);

  r[6] = vec_xxcpnmadd(r[2], U[0], r[6]);
  r[7] = vec_xxcpnmadd(r[2], U[1], r[7]);
  r[8] = vec_xxcpnmadd(r[2], U[2], r[8]);
  r[9] = vec_xxcpnmadd(r[5], U[0], r[9]);
  r[10] = vec_xxcpnmadd(r[5], U[1], r[10]);
  r[11] = vec_xxcpnmadd(r[5], U[2], r[11]);
  return;
}


// alternative implementation
//
// might not be optimal for pipeline as result is 
// re-used in the next line.
inline void vec_su3_multiply_double2b(su3 * const u, vector4double * U, vector4double * r) {
  __alignx(32, u);
  __alignx(32, U);
  __alignx(32, r);
  U[0] = vec_ld2(0, (double*) &u->c00);
  U[1] = vec_ld2(0, (double*) &u->c01);
  U[2] = vec_ld2(0, (double*) &u->c02);

  r[6] = vec_xmul(r[0], U[0]);
  r[6] = vec_xxnpmadd(U[0], r[0], r[6]);
#pragma unroll(2)
  for(int i = 1; i < 3; i++) {
    r[6] = vec_xmadd(r[i], U[i], r[6]);
    r[6] = vec_xxnpmadd(U[i], r[i], r[6]);
  }
  r[9] = vec_xmul(r[3], U[0]);
  r[9] = vec_xxnpmadd(U[0], r[3], r[9]);
#pragma unroll(2)
  for(int i = 1; i < 3; i++) {
    r[9] = vec_xmadd(r[3+i], U[i], r[9]);
    r[9] = vec_xxnpmadd(U[i], r[3+i], r[9]);
  }

  U[0] = vec_ld2(0, (double*) &u->c10);
  U[1] = vec_ld2(0, (double*) &u->c11);
  U[2] = vec_ld2(0, (double*) &u->c12);

  r[7] = vec_xmul(r[0], U[0]);
  r[7] = vec_xxnpmadd(U[0], r[0], r[7]);
#pragma unroll(2)
  for(int i = 1; i < 3; i++) {
    r[7] = vec_xmadd(r[i], U[i], r[7]);
    r[7] = vec_xxnpmadd(U[i], r[i], r[7]);
  }
  r[10] = vec_xmul(r[3], U[0]);
  r[10] = vec_xxnpmadd(U[0], r[3], r[10]);
#pragma unroll(2)
  for(int i = 1; i < 3; i++) {
    r[10] = vec_xmadd(r[3+i], U[i], r[10]);
    r[10] = vec_xxnpmadd(U[i], r[3+i], r[10]);
  }

  U[0] = vec_ld2(0, (double*) &u->c20);
  U[1] = vec_ld2(0, (double*) &u->c21);
  U[2] = vec_ld2(0, (double*) &u->c22);

  r[8] = vec_xmul(r[0], U[0]);
  r[8] = vec_xxnpmadd(U[0], r[0], r[8]);
#pragma unroll(2)
  for(int i = 1; i < 3; i++) {
    r[8] = vec_xmadd(r[i], U[i], r[8]);
    r[8] = vec_xxnpmadd(U[i], r[i], r[8]);
  }
  r[11] = vec_xmul(r[3], U[0]);
  r[11] = vec_xxnpmadd(U[0], r[3], r[11]);
#pragma unroll(2)
  for(int i = 1; i < 3; i++) {
    r[11] = vec_xmadd(r[3+i], U[i], r[11]);
    r[11] = vec_xxnpmadd(U[i], r[3+i], r[11]);
  }
  return;
}


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
  _vec_load_spinor(r4, r5, r6, r7, r8, r9, sp->s0);			\
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
  _vec_load(r4, r5, sp->s0);						\
  _vec_load16(r6, r7, sp->s1, U0);					\
  _vec_load(r10, r11, sp->s2);						\
  _vec_load16(r8, r9, sp->s3, U0);					\
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
  _vec_load(r4, r5, sp->s0);						\
  _vec_load16(r6, r7, sp->s1, U0);					\
  _vec_load(r10, r11, sp->s2);						\
  _vec_load16(r8, r9, sp->s3, U0);					\
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
  _vec_load(r4, r5, sp->s0);						\
  _vec_load16(r6, r7, sp->s1, U0);					\
  _vec_load(r8, r9, sp->s2);						\
  _vec_load16(r10, r11, sp->s3, U0);					\
  _vec_i_mul_sub(r0, r1, r4, r5, r8, r9, U0);				\
  _vec_i_mul_add(r2, r3, r6, r7, r10, r11, U1);				\
  _vec_su3_inverse_multiply_double2c(um);				\
  _vec_cmplxcg_mul_double2c(r0, r1, r2, r4, r5, r6, rtmp);		\
  _vec_unfuse(r0, r1, r2, r3, r4, r5);					\
  _vec_add_double2(rs0, rs1, rs2, rs3, rs4, rs5, r0, r1, r2, r3, r4, r5); \
  _vec_i_mul_add2(rs6, rs7, rs8, r0, r1, r2, U0);			\
  _vec_i_mul_sub2(rs9, rs10, rs11, r3, r4, r5, U1);


#endif

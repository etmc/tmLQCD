#ifndef _BGQ_H
#define _BGQ_H

#include "bgq2.h"

// requires 32 byte alignment of phi
#define _vec_load(r0, r1, phi)			\
  r0 = vec_ld(0, (double*) &(phi).c0);		\
  r1 = vec_ld2(0, (double*) &(phi).c2); 

// works also with 16 byte alignement of phi
#define _vec_load16(r0, r1, phi, tmp)		\
  r0 = vec_ld2(0, (double*) &(phi).c0);		\
  r1 = vec_ld(0, (double*) &(phi).c1);		\
  tmp = vec_gpci(00145);			\
  r0 = vec_perm(r0, r1, tmp);			\
  tmp = vec_gpci(02301);			\
  r1 = vec_perm(r1, r0, tmp);

// alternative
#define _vec_load16c(r0, r1, phi, tmp)		\
  r0 = vec_ld2(0, (double*) &(phi).c0);		\
  r1 = vec_ld(0, (double*) &(phi).c1);		\
  tmp = vec_gpci(00145);			\
  r0 = vec_perm(r0, r1, tmp);			\
  r1 = vec_ld2(0, (double*) &(phi).c2);

// requires 32 byte alignment of phi
#define _vec_store(phi, r0, r1)				\
  vec_st((r0), 0, (double*) &(phi).c0);			\
  vec_st2((r1), 0, (double*) &(phi).c2);

#define _vec_add(rs0, rs1, r0, r1, s0, s1) \
  rs0 = vec_add(r0, s0);		   \
  rs1 = vec_add(r1, s1);

#define _vec_sub(rs0, rs1, r0, r1, s0, s1) \
  rs0 = vec_sub(r0, s0);		   \
  rs1 = vec_sub(r1, s1);

#define _vec_i_mul_add(rs0, rs1, r0, r1, s0, s1, tmp)	\
  tmp = vec_splats(1.);					\
  rs0 = vec_xxnpmadd(s0, tmp, r0);			\
  rs1 = vec_xxnpmadd(s1, tmp, r1);

#define _vec_i_mul_sub(rs0, rs1, r0, r1, s0, s1, tmp)	\
  tmp = vec_splats(-1.);					\
  rs0 = vec_xxnpmadd(s0, tmp, r0);			\
  rs1 = vec_xxnpmadd(s1, tmp, r1);

#define _vec_cmplx_mul_double(rs0, rs1, rs2, rs3, r0, r1, r2, r3, tmp, c) \
  tmp = vec_ld2(0, (double*) &(c));					\
  rs0 = vec_xmul(r0, tmp);						\
  rs1 = vec_xmul(r1, tmp);						\
  rs2 = vec_xmul(r2, tmp);						\
  rs3 = vec_xmul(r3, tmp);						\
  rs0 = vec_xxnpmadd(tmp, r0, rs0);					\
  rs1 = vec_xxnpmadd(tmp, r1, rs1);					\
  rs2 = vec_xxnpmadd(tmp, r2, rs2);					\
  rs3 = vec_xxnpmadd(tmp, r3, rs3);

#define _vec_cmplxcg_mul_double(rs0, rs1, rs2, rs3, r0, r1, r2, r3, tmp, c) \
  tmp = vec_ld2(0, (double*) &(c));					\
  rs0 = vec_xmul(tmp, r0);						\
  rs1 = vec_xmul(tmp, r1);						\
  rs2 = vec_xmul(tmp, r2);						\
  rs3 = vec_xmul(tmp, r3);						\
  rs0 = vec_xxcpnmadd(r0, tmp, rs0);					\
  rs1 = vec_xxcpnmadd(r1, tmp, rs1);					\
  rs2 = vec_xxcpnmadd(r2, tmp, rs2);					\
  rs3 = vec_xxcpnmadd(r3, tmp, rs3);

#endif

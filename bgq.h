#ifndef _BGQ_H
#define _BGQ_H

#include "bgq2.h"

#define _vec_load_spinor(r0, r1, r2, r3, r4, r5, phi)	\
  r0 = vec_ld(0L, (double*) &(phi).c0);			\
  r1 = vec_ld(32L, (double*) &(phi).c0);		\
  r2 = vec_ld(64L, (double*) &(phi).c0);		\
  r3 = vec_ld(96L, (double*) &(phi).c0);		\
  r4 = vec_ld(128L, (double*) &(phi).c0);		\
  r5 = vec_ld(160L, (double*) &(phi).c0);

#define _vec_load_halfspinor(r0, r1, r2, phi)	\
  r0 = vec_ld(0L, (double*) &(phi).c0);			\
  r1 = vec_ld(32L, (double*) &(phi).c0);		\
  r2 = vec_ld(64L, (double*) &(phi).c0);


#define _vec_store_spinor(phi, r0, r1, r2, r3, r4, r5) \
  vec_st(r0, 0L, (double*) &(phi).c0);		       \
  vec_st(r1, 32L, (double*) &(phi).c0);		       \
  vec_st(r2, 64L, (double*) &(phi).c0);		       \
  vec_st(r3, 96L, (double*) &(phi).c0);		       \
  vec_st(r4, 128L, (double*) &(phi).c0);	       \
  vec_st(r5, 160L, (double*) &(phi).c0);

#define _vec_add_ul_spinor(rs0, rs1, rs2, r0, r1, r2, r3, r4, r5) \
  rs0 = vec_add(r0, r3);					  \
  rs1 = vec_add(r1, r4);					  \
  rs2 = vec_add(r2, r5);

#define _vec_sub_ul_spinor(rs0, rs1, rs2, r0, r1, r2, r3, r4, r5) \
  rs0 = vec_sub(r0, r3);					  \
  rs1 = vec_sub(r1, r4);					  \
  rs2 = vec_sub(r2, r5);

// requires 32 byte alignment of phi
#define _vec_load(r0, r1, phi)			\
  r0 = vec_ld(0L, (double*) &(phi).c0);		\
  r1 = vec_ld2(0L, (double*) &(phi).c2); 

// works also with 16 byte alignement of phi
#define _vec_load16(r0, r1, phi, tmp)		\
  r0 = vec_ld2(0L, (double*) &(phi).c0);	\
  r1 = vec_ld(0L, (double*) &(phi).c1);		\
  tmp = vec_gpci(00145);			\
  r0 = vec_perm(r0, r1, tmp);			\
  tmp = vec_gpci(02301);			\
  r1 = vec_perm(r1, r0, tmp);

// alternative
#define _vec_load16c(r0, r1, phi, tmp)		\
  r0 = vec_ld2(0L, (double*) &(phi).c0);	\
  r1 = vec_ld(0L, (double*) &(phi).c1);		\
  tmp = vec_gpci(00145);			\
  r0 = vec_perm(r0, r1, tmp);			\
  r1 = vec_ld2(0L, (double*) &(phi).c2);

// requires 32 byte alignment of phi
#define _vec_store(phi, r0, r1)			\
  vec_st((r0), 0L, (double*) &(phi).c0);	\
  vec_st2((r1), 0L, (double*) &(phi).c2);

// requires 16 (and must not be 32) byte alignment of phi
#define _vec_store16(phi, r0, r1, tmp)		\
  vec_st2((r0), 0L, (double*) &(phi).c0);	\
  tmp = vec_gpci(02345);			\
  r0 = vec_perm(r0, r1, tmp);			\
  vec_st((r0), 0L, (double *) &(phi).c1);

// requires 32 byte alignment of phi
#define _vec_store_halfspinor(phi, r0, r1, r2)	\
  vec_st((r0), 0L, (double*) &(phi).c0);	\
  vec_st((r1), 32L, (double*) &(phi).c0);	\
  vec_st((r2), 64L, (double*) &(phi).c0);

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

#endif

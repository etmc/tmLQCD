#ifndef _SU3_H
#define _SU3_H

/***********************************************************************
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008 Carsten Urbach
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
 * File su3.h
 *
 * Type definitions and macros for SU(3) matrices and spinors  
 *
 * Version: 1.0
 * Author: Martin Luescher <luscher@mail.desy.de>
 * Date: 24.10.2000
 *
 * Extended by Martin Hasenbusch 2001.  <Martin.Hasenbusch@desy.de>
 * Rewritten for C99 complex by Albert Deuzeman 2012 <deuzeman@itp.unibe.ch>
 *
 *******************************************************************************/

#include <complex.h>
#if (defined XLC && defined BGL)
# include "bgl.h"
#endif

typedef struct 
{
   _Complex double c00, c01, c02, c10, c11, c12, c20, c21, c22;
} su3;

typedef struct
{
   _Complex double c0,c1,c2;
} su3_vector;

typedef struct
{
    _Complex float c0,c1,c2;
} su3_vector32;

typedef struct
{
   su3_vector s0,s1,s2,s3;
} spinor;

typedef struct
{
  su3_vector s0, s1;
} halfspinor;

typedef struct
{
   su3_vector32 s0, s1;
} halfspinor32;

typedef struct
{
   spinor sp_up,sp_dn;
} bispinor;

typedef struct
{
  _Complex double s00,s01,s02,s03,s10,s11,s12,s13,s20,s21,s22,s23,s30,s31,s32,s33;
} spinor_matrix;

typedef struct
{
  _Complex double sc0,sc1,sc2,sc3;
} complex_spinor;


/*******************************************************************************
*
* Macros for SU(3) vectors
*
* Arguments are variables of type su3_vector (or su3 in the case of the
* matrix times vector multiplication macros)
*
*******************************************************************************/
/* M. Hasenbusch Mon Sep 24
* r.c1=0
* r.c2=0
* r.c3=0
*/
#define _vector_null(r) \
  (r).c0 = 0.;		\
  (r).c1 = 0.;		\
  (r).c2 = 0.;

/* M. Hasenbusch Mon Sep 24
* r.c1=s.c1
* r.c2=s.c2
* r.c3=s.c3
*/

/* The following should be taken care of by the compiler, actually... Just redefine as _vector_assign */
#define _vector_assign32(r,s) _vector_assign(r,s)

#define _vector_assign(r,s)			\
  (r).c0 = (s).c0;				\
  (r).c1 = (s).c1;				\
  (r).c2 = (s).c2;

   
#define _spinor_assign(r,s) \
   _vector_assign((r).s0,(s).s0);\
   _vector_assign((r).s1,(s).s1);\
   _vector_assign((r).s2,(s).s2);\
   _vector_assign((r).s3,(s).s3);

#define _vector_norm_square(r) \
   conj((r).c0) * (r).c0 + conj((r).c1) * (r).c1 + conj((r).c2) * (r).c2
   
#define _vector_minus_assign(r,s)		\
  (r).c0 = -(s).c0;				\
  (r).c1 = -(s).c1;				\
  (r).c2 = -(s).c2;

#define _vector_mul(r,c,s)			\
  (r).c0 = (c) * (s).c0;			\
  (r).c1 = (c) * (s).c1;			\
  (r).c2 = (c) * (s).c2;

#define _vector_add_mul(r,c,s)			\
  (r).c0 += (c) * (s).c0;			\
  (r).c1 += (c) * (s).c1;			\
  (r).c2 += (c) * (s).c2;


/* r += I * c * s (c real) */
#define _vector_add_i_mul(r,c,s) \
  (r).c0 += I*(c)*(s).c0; \
  (r).c1 += I*(c)*(s).c1; \
  (r).c2 += I*(c)*(s).c2;

#if ((defined SSE2)||(defined SSE3))

#define _vector_add(r,s1,s2) \
_sse_load(s1); \
_sse_load_up(s2); \
_sse_vector_add(); \
_sse_store(r);

#define _vector_sub(r,s1,s2) \
_sse_load(s1); \
_sse_load_up(s2); \
_sse_vector_sub(); \
_sse_store(r);

#elif (defined XLC && defined BGLNOTCHECKED)

#define _vector_add(r,s1,s2) \
  _bgl_load(s1);	     \
  _bgl_load_up(s2);	     \
  _bgl_vector_add();	     \
  _bgl_store(r);

#define _vector_sub(r,s1,s2) \
  _bgl_load(s1);	     \
  _bgl_load_up(s2);	     \
  _bgl_vector_sub();	     \
  _bgl_store(r);

#else

#define _vector_add(r,s1,s2)		\
  (r).c0 = (s1).c0 + (s2).c0;		\
  (r).c1 = (s1).c1 + (s2).c1;		\
  (r).c2 = (s1).c2 + (s2).c2;

#define _vector_sub(r,s1,s2)		\
  (r).c0 = (s1).c0 - (s2).c0;		\
  (r).c1 = (s1).c1 - (s2).c1;		\
  (r).c2 = (s1).c2 - (s2).c2;
#endif


#define _vector_i_add(r,s1,s2)		\
  (r).c0 = (s1).c0 + I * (s2).c0;	\
  (r).c1 = (s1).c1 + I * (s2).c1;	\
  (r).c2 = (s1).c2 + I * (s2).c2;

#define _vector_i_sub(r,s1,s2)		\
  (r).c0 = (s1).c0 - I * (s2).c0;	\
  (r).c1 = (s1).c1 - I * (s2).c1;	\
  (r).c2 = (s1).c2 - I * (s2).c2;

#define _vector_combined_add_i_add(r1, s1, r2, s2, s)	\
  (r1).c0 = (s1).c0 + (s).c0;				\
  (r2).c0 = (s2).c0 + I * (s).c0;			\
  (r1).c1 = (s1).c1 + (s).c1;				\
  (r2).c1 = (s2).c1 + I * (s).c1;			\
  (r1).c2 = (s1).c2 + (s).c2;				\
  (r2).c2 = (s2).c2 + I * (s).c2;			\

#if ((defined SSE2) || (defined SSE3))

#define _vector_add_assign(r,s) \
_sse_load(r); \
_sse_load_up(s); \
_sse_vector_add(); \
_sse_store(r);

#define _vector_sub_assign(r,s) \
_sse_load(r); \
_sse_load_up(s); \
_sse_vector_sub(); \
_sse_store(r);

#else

#define _vector_add_assign(r,s)			\
  (r).c0 += (s).c0;				\
  (r).c1 += (s).c1;				\
  (r).c2 += (s).c2;

#define _vector_sub_assign(r,s)			\
  (r).c0 -= (s).c0;				\
  (r).c1 -= (s).c1;				\
  (r).c2 -= (s).c2;

#endif 

#define _vector_i_add_assign(r,s)		\
  (r).c0 += I * (s).c0;			\
  (r).c1 += I * (s).c1;			\
  (r).c2 += I * (s).c2;


#define _vector_i_sub_assign(r,s)		\
  (r).c0 -= I * (s).c0;				\
  (r).c1 -= I * (s).c1;				\
  (r).c2 -= I * (s).c2;

#define complex_times_vector(r,c,s)		\
  (r).c0 = (c) * (s).c0;			\
  (r).c1 = (c) * (s).c1;			\
  (r).c2 = (c) * (s).c2;

/* M.Hasenbusch */
#define _complexcjg_times_vector(r,c,s)		\
  (r).c0 = conj(c) * (s).c0;			\
  (r).c1 = conj(c) * (s).c1;			\
  (r).c2 = conj(c) * (s).c2;

#define _vector_project(r,z,s)			\
  (r).c0 -= z * (s).c0;				\
  (r).c1 -= z * (s).c1;		        	\
  (r).c2 -= z * (s).c2;

#if ((defined SSE2) || (defined SSE3))

#define _su3_multiply(r,u,s) \
_sse_load(s); \
_sse_su3_multiply(u); \
_sse_store_up(r);

#define _su3_inverse_multiply(r,u,s) \
_sse_load(s); \
_sse_su3_inverse_multiply(u); \
_sse_store_up(r);

#elif (defined XLC && defined BGLNOTCHECKED)

#define _su3_multiply(r,u,s) \
  _bgl_load(s); \
  _bgl_su3_multiply(u); \
  _bgl_store_up(r);

#define _su3_inverse_multiply(r,u,s)		\
  _bgl_load(s); \
  _bgl_su3_inverse_multiply(u); \
  _bgl_store_up(r);

#else

#define _su3_multiply(r,u,s)					        \
  (r).c0 = (u).c00 * (s).c0 + (u).c01 * (s).c1 + (u).c02 * (s).c2;	\
  (r).c1 = (u).c10 * (s).c0 + (u).c11 * (s).c1 + (u).c12 * (s).c2;	\
  (r).c2 = (u).c20 * (s).c0 + (u).c21 * (s).c1 + (u).c22 * (s).c2;

#define _su3_inverse_multiply(r,u,s)		\
(r).c0 = conj((u).c00) * (s).c0 + conj((u).c10) * (s).c1 + conj((u).c20) * (s).c2;	\
(r).c1 = conj((u).c01) * (s).c0 + conj((u).c11) * (s).c1 + conj((u).c21) * (s).c2;	\
(r).c2 = conj((u).c02) * (s).c0 + conj((u).c12) * (s).c1 + conj((u).c22) * (s).c2;   

#endif

/*******************************************************************************
*
* Macros for SU(3) matrices
*
* Arguments are variables of type su3
*
*******************************************************************************/

#define _su3_norm_sq(x,u)			\
  x = (u).c00 * conj((u).c00) + (u).c01 * conj((u).c01) + (u).c02 * conj((u).c02) \
      (u).c10 * conj((u).c10) + (u).c11 * conj((u).c11) + (u).c12 * conj((u).c12) \
      (u).c20 * conj((u).c20) + (u).c21 * conj((u).c21) + (u).c22 * conj((u).c22); 

#define _su3_one(u)	\
  (u).c00 = 1.;		\
  (u).c01 = 0.;		\
  (u).c02 = 0.;		\
  (u).c10 = 0.;		\
  (u).c11 = 1.;		\
  (u).c12 = 0.;		\
  (u).c20 = 0.;		\
  (u).c21 = 0.;		\
  (u).c22 = 1.;

#define _su3_zero(u)		\
  (u).c00 = 0.;			\
  (u).c01 = 0.;			\
  (u).c02 = 0.;			\
  (u).c10 = 0.;			\
  (u).c11 = 0.;			\
  (u).c12 = 0.;			\
  (u).c20 = 0.;			\
  (u).c21 = 0.;			\
  (u).c22 = 0.;

#define _su3_assign(u,v)			\
  (u).c00 = (v).c00;				\
  (u).c01 = (v).c01;				\
  (u).c02 = (v).c02;				\
  (u).c10 = (v).c10;				\
  (u).c11 = (v).c11;				\
  (u).c12 = (v).c12;				\
  (u).c20 = (v).c20;				\
  (u).c21 = (v).c21;				\
  (u).c22 = (v).c22;

#define _su3_minus_assign(u,v)			\
  (u).c00 = -(v).c00;				\
  (u).c01 = -(v).c01;				\
  (u).c02 = -(v).c02;				\
  (u).c10 = -(v).c10;				\
  (u).c11 = -(v).c11;				\
  (u).c12 = -(v).c12;				\
  (u).c20 = -(v).c20;				\
  (u).c21 = -(v).c21;				\
  (u).c22 = -(v).c22;

#define _su3_dagger(u,v)			\
  (u).c00 = conj((v).c00);			\
  (u).c01 = conj((v).c10);			\
  (u).c02 = conj((v).c20);			\
  (u).c10 = conj((v).c01);			\
  (u).c11 = conj((v).c11);			\
  (u).c12 = conj((v).c21);			\
  (u).c20 = conj((v).c02);			\
  (u).c21 = conj((v).c12);			\
  (u).c22 = conj((v).c22); 

#define _itimes_su3(u,v)			\
  (u).c00 = I * (v).c00;			\
  (u).c01 = I * (v).c01;			\
  (u).c02 = I * (v).c02;			\
  (u).c10 = I * (v).c10;			\
  (u).c11 = I * (v).c11;			\
  (u).c12 = I * (v).c12;			\
  (u).c20 = I * (v).c20;			\
  (u).c21 = I * (v).c21;			\
  (u).c22 = I * (v).c22;

#define _real_times_su3(u,a,v)			\
  (u).c00 = (a) * (v).c00;			\
  (u).c01 = (a) * (v).c01;			\
  (u).c02 = (a) * (v).c02;			\
  (u).c10 = (a) * (v).c10;			\
  (u).c11 = (a) * (v).c11;			\
  (u).c12 = (a) * (v).c12;			\
  (u).c20 = (a) * (v).c20;			\
  (u).c21 = (a) * (v).c21;			\
  (u).c22 = (a) * (v).c22;

#define _real_times_su3_plus_real_times_su3(u, a, v, b, w) \
  (u).c00 = (a)*(v).c00 + (b)*(w).c00; \
  (u).c01 = (a)*(v).c01 + (b)*(w).c01; \
  (u).c02 = (a)*(v).c02 + (b)*(w).c02; \
  (u).c10 = (a)*(v).c10 + (b)*(w).c10; \
  (u).c11 = (a)*(v).c11 + (b)*(w).c11; \
  (u).c12 = (a)*(v).c12 + (b)*(w).c12; \
  (u).c20 = (a)*(v).c20 + (b)*(w).c20; \
  (u).c21 = (a)*(v).c21 + (b)*(w).c21; \
  (u).c22 = (a)*(v).c22 + (b)*(w).c22;

#define _su3_minus_su3(u,v,w) \
   (u).c00 = (v).c00 - (w).c00; \
   (u).c01 = (v).c01 - (w).c01; \
   (u).c02 = (v).c02 - (w).c02; \
   (u).c10 = (v).c10 - (w).c10; \
   (u).c11 = (v).c11 - (w).c11; \
   (u).c12 = (v).c12 - (w).c12; \
   (u).c20 = (v).c20 - (w).c20; \
   (u).c21 = (v).c21 - (w).c21; \
   (u).c22 = (v).c22 - (w).c22; \

#define _itimes_su3_minus_su3(u,v,w)  \
   (u).c00 = I * ((v).c00 - (w).c00); \
   (u).c01 = I * ((v).c01 - (w).c01); \
   (u).c02 = I * ((v).c02 - (w).c02); \
   (u).c10 = I * ((v).c10 - (w).c10); \
   (u).c11 = I * ((v).c11 - (w).c11); \
   (u).c12 = I * ((v).c12 - (w).c12); \
   (u).c20 = I * ((v).c20 - (w).c20); \
   (u).c21 = I * ((v).c21 - (w).c21); \
   (u).c22 = I * ((v).c22 - (w).c22);

#define _su3_plus_su3(u,v,w) \
   (u).c00 = (v).c00 + (w).c00; \
   (u).c01 = (v).c01 + (w).c01; \
   (u).c02 = (v).c02 + (w).c02; \
   (u).c10 = (v).c10 + (w).c10; \
   (u).c11 = (v).c11 + (w).c11; \
   (u).c12 = (v).c12 + (w).c12; \
   (u).c20 = (v).c20 + (w).c20; \
   (u).c21 = (v).c21 + (w).c21; \
   (u).c22 = (v).c22 + (w).c22;

#define _minus_su3_plus_su3(u,v,w) \
   (u).c00 = -((v).c00 + (w).c00); \
   (u).c01 = -((v).c01 + (w).c01); \
   (u).c02 = -((v).c02 + (w).c02); \
   (u).c10 = -((v).c10 + (w).c10); \
   (u).c11 = -((v).c11 + (w).c11); \
   (u).c12 = -((v).c12 + (w).c12); \
   (u).c20 = -((v).c20 + (w).c20); \
   (u).c21 = -((v).c21 + (w).c21); \
   (u).c22 = -((v).c22 + (w).c22);

#define _itimes_su3_plus_su3(u,v,w) \
   (u).c00 = I * ((v).c00 + (w).c00); \
   (u).c01 = I * ((v).c01 + (w).c01); \
   (u).c02 = I * ((v).c02 + (w).c02); \
   (u).c10 = I * ((v).c10 + (w).c10); \
   (u).c11 = I * ((v).c11 + (w).c11); \
   (u).c12 = I * ((v).c12 + (w).c12); \
   (u).c20 = I * ((v).c20 + (w).c20); \
   (u).c21 = I * ((v).c21 + (w).c21); \
   (u).c22 = I * ((v).c22 + (w).c22);

#define _minus_itimes_su3_plus_su3(u,v,w) \
  (u).c00 = -I * ((v).c00 + (w).c00); \
  (u).c01 = -I * ((v).c01 + (w).c01); \
  (u).c02 = -I * ((v).c02 + (w).c02); \
  (u).c10 = -I * ((v).c10 + (w).c10); \
  (u).c11 = -I * ((v).c11 + (w).c11); \
  (u).c12 = -I * ((v).c12 + (w).c12); \
  (u).c20 = -I * ((v).c20 + (w).c20); \
  (u).c21 = -I * ((v).c21 + (w).c21); \
  (u).c22 = -I * ((v).c22 + (w).c22);

#define _complex_times_su3(r,c,s) \
   (r).c00 = (c) * (s).c00; \
   (r).c01 = (c) * (s).c01; \
   (r).c02 = (c) * (s).c02; \
   (r).c10 = (c) * (s).c10; \
   (r).c11 = (c) * (s).c11; \
   (r).c12 = (c) * (s).c12; \
   (r).c20 = (c) * (s).c20; \
   (r).c21 = (c) * (s).c21; \
   (r).c22 = (c) * (s).c22;

#define _complexcjg_times_su3(r,c,s) \
   (r).c00 = conj(c) * (s).c00; \
   (r).c01 = conj(c) * (s).c01; \
   (r).c02 = conj(c) * (s).c02; \
   (r).c10 = conj(c) * (s).c10; \
   (r).c11 = conj(c) * (s).c11; \
   (r).c12 = conj(c) * (s).c12; \
   (r).c20 = conj(c) * (s).c20; \
   (r).c21 = conj(c) * (s).c21; \
   (r).c22 = conj(c) * (s).c22;


/* M. Hasenbusch
* su3_acc
*/
#if ((defined SSE2) || (defined SSE3))
#define _su3_acc(u,v) _sse_su3_acc(u,v) 
#else
#define _su3_acc(u,v) \
   (u).c00 += (v).c00; \
   (u).c01 += (v).c01; \
   (u).c02 += (v).c02; \
   (u).c10 += (v).c10; \
   (u).c11 += (v).c11; \
   (u).c12 += (v).c12; \
   (u).c20 += (v).c20; \
   (u).c21 += (v).c21; \
   (u).c22 += (v).c22;
#endif

/*
* su3_refac_acc
*/
#define _su3_refac_acc(u,a,v) \
   (u).c00 += a * (v).c00; \
   (u).c01 += a * (v).c01; \
   (u).c02 += a * (v).c02; \
   (u).c10 += a * (v).c10; \
   (u).c11 += a * (v).c11; \
   (u).c12 += a * (v).c12; \
   (u).c20 += a * (v).c20; \
   (u).c21 += a * (v).c21; \
   (u).c22 += a * (v).c22;

/*
* su3_imfac_acc
*/
#define _su3_imfac_acc(u,a,v) \
   (u).c00 += I * a * (v).c00; \
   (u).c01 += I * a * (v).c01; \
   (u).c02 += I * a * (v).c02; \
   (u).c10 += I * a * (v).c10; \
   (u).c11 += I * a * (v).c11; \
   (u).c12 += I * a * (v).c12; \
   (u).c20 += I * a * (v).c20; \
   (u).c21 += I * a * (v).c21; \
   (u).c22 += I * a * (v).c22;

#define _su3_square_norm(s, v) \
  s = conj(v.c00) * (v.c00) + conj(v.c01) * (v.c01) + conj(v.c02) * (v.c02) + \
    conj(v.c10) * (v.c10) + conj(v.c11) * (v.c11) + conj(v.c12) * (v.c12) + \
    conj(v.c20) * (v.c20) + conj(v.c21) * (v.c21) + conj(v.c22) * (v.c22);

#if ((defined SSE2) || (defined SSE3))

#define _su3_times_su3(u,v,w) _sse_su3_times_su3(u,v,w)
#define _su3_times_su3_acc(u,v,w) _sse_su3_times_su3_acc(u,v,w)
#define _su3d_times_su3(u,v,w) _sse_su3d_times_su3(u,v,w)
#define _su3d_times_su3_acc(u,v,w) _sse_su3d_times_su3_acc(u,v,w)
#define _su3_times_su3d(u,v,w) _sse_su3_times_su3d(u,v,w)
#define _su3_times_su3d_acc(u,v,w) _sse_su3_times_su3d_acc(u,v,w)  

#else

#define _su3_times_su3(u,v,w)					\
  (u).c00 = (v).c00 * (w).c00 + (v).c01 * (w).c10 + (v).c02*(w).c20;	\
  (u).c01 = (v).c00 * (w).c01 + (v).c01 * (w).c11 + (v).c02*(w).c21;	\
  (u).c02 = (v).c00 * (w).c02 + (v).c01 * (w).c12 + (v).c02*(w).c22;	\
  (u).c10 = (v).c10 * (w).c00 + (v).c11 * (w).c10 + (v).c12*(w).c20;	\
  (u).c11 = (v).c10 * (w).c01 + (v).c11 * (w).c11 + (v).c12*(w).c21;	\
  (u).c12 = (v).c10 * (w).c02 + (v).c11 * (w).c12 + (v).c12*(w).c22;	\
  (u).c20 = (v).c20 * (w).c00 + (v).c21 * (w).c10 + (v).c22*(w).c20;	\
  (u).c21 = (v).c20 * (w).c01 + (v).c21 * (w).c11 + (v).c22*(w).c21;	\
  (u).c22 = (v).c20 * (w).c02 + (v).c21 * (w).c12 + (v).c22*(w).c22;	\

#define _su3_times_su3_acc(u,v,w)					\
  (u).c00 += (v).c00 * (w).c00 + (v).c01*(w).c10 + (v).c02*(w).c20;	\
  (u).c01 += (v).c00 * (w).c01 + (v).c01*(w).c11 + (v).c02*(w).c21;	\
  (u).c02 += (v).c00 * (w).c02 + (v).c01*(w).c12 + (v).c02*(w).c22;	\
  (u).c10 += (v).c10 * (w).c00 + (v).c11*(w).c10 + (v).c12*(w).c20;	\
  (u).c11 += (v).c10 * (w).c01 + (v).c11*(w).c11 + (v).c12*(w).c21;	\
  (u).c12 += (v).c10 * (w).c02 + (v).c11*(w).c12 + (v).c12*(w).c22;	\
  (u).c20 += (v).c20 * (w).c00 + (v).c21*(w).c10 + (v).c22*(w).c20;	\
  (u).c21 += (v).c20 * (w).c01 + (v).c21*(w).c11 + (v).c22*(w).c21;	\
  (u).c22 += (v).c20 * (w).c02 + (v).c21*(w).c12 + (v).c22*(w).c22;

#define _su3_times_su3d(u,v,w)						\
  (u).c00 =  (v).c00 * conj((w).c00) + (v).c01 * conj((w).c01) + (v).c02 * conj((w).c02); \
  (u).c01 =  (v).c00 * conj((w).c10) + (v).c01 * conj((w).c11) + (v).c02 * conj((w).c12); \
  (u).c02 =  (v).c00 * conj((w).c20) + (v).c01 * conj((w).c21) + (v).c02 * conj((w).c22); \
  (u).c10 =  (v).c10 * conj((w).c00) + (v).c11 * conj((w).c01) + (v).c12 * conj((w).c02); \
  (u).c11 =  (v).c10 * conj((w).c10) + (v).c11 * conj((w).c11) + (v).c12 * conj((w).c12); \
  (u).c12 =  (v).c10 * conj((w).c20) + (v).c11 * conj((w).c21) + (v).c12 * conj((w).c22); \
  (u).c20 =  (v).c20 * conj((w).c00) + (v).c21 * conj((w).c01) + (v).c22 * conj((w).c02); \
  (u).c21 =  (v).c20 * conj((w).c10) + (v).c21 * conj((w).c11) + (v).c22 * conj((w).c12); \
  (u).c22 =  (v).c20 * conj((w).c20) + (v).c21 * conj((w).c21) + (v).c22 * conj((w).c22);

#define _su3_times_su3d_acc(u,v,w)						\
  (u).c00 += (v).c00 * conj((w).c00) + (v).c01 * conj((w).c01)  + (v).c02 * conj((w).c02); \
  (u).c01 += (v).c00 * conj((w).c10) + (v).c01 * conj((w).c11)  + (v).c02 * conj((w).c12); \
  (u).c02 += (v).c00 * conj((w).c20) + (v).c01 * conj((w).c21)  + (v).c02 * conj((w).c22); \
  (u).c10 += (v).c10 * conj((w).c00) + (v).c11 * conj((w).c01)  + (v).c12 * conj((w).c02); \
  (u).c11 += (v).c10 * conj((w).c10) + (v).c11 * conj((w).c11)  + (v).c12 * conj((w).c12); \
  (u).c12 += (v).c10 * conj((w).c20) + (v).c11 * conj((w).c21)  + (v).c12 * conj((w).c22); \
  (u).c20 += (v).c20 * conj((w).c00) + (v).c21 * conj((w).c01)  + (v).c22 * conj((w).c02); \
  (u).c21 += (v).c20 * conj((w).c10) + (v).c21 * conj((w).c11)  + (v).c22 * conj((w).c12); \
  (u).c22 += (v).c20 * conj((w).c20) + (v).c21 * conj((w).c21)  + (v).c22 * conj((w).c22);

#define _su3d_times_su3(u,v,w)			\
  (u).c00 = conj((v).c00) * (w).c00 + conj((v).c10) * (w).c10 + conj((v).c20) * (w).c20; \
  (u).c01 = conj((v).c00) * (w).c01 + conj((v).c10) * (w).c11 + conj((v).c20) * (w).c21; \
  (u).c02 = conj((v).c00) * (w).c02 + conj((v).c10) * (w).c12 + conj((v).c20) * (w).c22; \
  (u).c10 = conj((v).c01) * (w).c00 + conj((v).c11) * (w).c10 + conj((v).c21) * (w).c20; \
  (u).c11 = conj((v).c01) * (w).c01 + conj((v).c11) * (w).c11 + conj((v).c21) * (w).c21; \
  (u).c12 = conj((v).c01) * (w).c02 + conj((v).c11) * (w).c12 + conj((v).c21) * (w).c22; \
  (u).c20 = conj((v).c02) * (w).c00 + conj((v).c12) * (w).c10 + conj((v).c22) * (w).c20; \
  (u).c21 = conj((v).c02) * (w).c01 + conj((v).c12) * (w).c11 + conj((v).c22) * (w).c21; \
  (u).c22 = conj((v).c02) * (w).c02 + conj((v).c12) * (w).c12 + conj((v).c22) * (w).c22;

#define _su3d_times_su3_acc(u,v,w)						\
  (u).c00 += conj((v).c00) * (w).c00 + conj((v).c10) * (w).c10 + conj((v).c20) * (w).c20; \
  (u).c01 += conj((v).c00) * (w).c01 + conj((v).c10) * (w).c11 + conj((v).c20) * (w).c21; \
  (u).c02 += conj((v).c00) * (w).c02 + conj((v).c10) * (w).c12 + conj((v).c20) * (w).c22; \
  (u).c10 += conj((v).c01) * (w).c00 + conj((v).c11) * (w).c10 + conj((v).c21) * (w).c20; \
  (u).c11 += conj((v).c01) * (w).c01 + conj((v).c11) * (w).c11 + conj((v).c21) * (w).c21; \
  (u).c12 += conj((v).c01) * (w).c02 + conj((v).c11) * (w).c12 + conj((v).c21) * (w).c22; \
  (u).c20 += conj((v).c02) * (w).c00 + conj((v).c12) * (w).c10 + conj((v).c22) * (w).c20; \
  (u).c21 += conj((v).c02) * (w).c01 + conj((v).c12) * (w).c11 + conj((v).c22) * (w).c21; \
  (u).c22 += conj((v).c02) * (w).c02 + conj((v).c12) * (w).c12 + conj((v).c22) * (w).c22;

#endif

#define _trace_su3_times_su3d(x,v,w)	\
  x =   (v).c00 * conj((w).c00)		\
      + (v).c01 * conj((w).c01)		\
      + (v).c02 * conj((w).c02)		\
      + (v).c10 * conj((w).c10)		\
      + (v).c11 * conj((w).c11)		\
      + (v).c12 * conj((w).c12)		\
      + (v).c20 * conj((w).c20)		\
      + (v).c21 * conj((w).c21)		\
      + (v).c22 * conj((w).c22);

#define _trace_su3_times_su3(x,v,w)	\
  x = (v).c00 * (w).c00			\
    + (v).c01 * (w).c10			\
    + (v).c02 * (w).c20			\
    + (v).c10 * (w).c01			\
    + (v).c11 * (w).c11			\
    + (v).c12 * (w).c21			\
    + (v).c20 * (w).c02			\
    + (v).c21 * (w).c12			\
    + (v).c22 * (w).c22;

#define _complex_times_vector(x, c, y)	\
   x.c0 = (c) * (y).c0;			\
   x.c1 = (c) * (y).c1;			\
   x.c2 = (c) * (y).c2;
    
#define _vector_tensor_vector(t,u,v)	\
  (t).c00 = (u).c0 * conj((v).c0);	\
  (t).c01 = (u).c0 * conj((v).c1);	\
  (t).c02 = (u).c0 * conj((v).c2);	\
  (t).c10 = (u).c1 * conj((v).c0);	\
  (t).c11 = (u).c1 * conj((v).c1);	\
  (t).c12 = (u).c1 * conj((v).c2);	\
  (t).c20 = (u).c2 * conj((v).c0);	\
  (t).c21 = (u).c2 * conj((v).c1);	\
  (t).c22 = (u).c2 * conj((v).c2);

#define _mvector_tensor_vector(t,u,v)	\
  (t).c00 = -(u).c0 * conj((v).c0);	\
  (t).c01 = -(u).c0 * conj((v).c1);	\
  (t).c02 = -(u).c0 * conj((v).c2);	\
  (t).c10 = -(u).c1 * conj((v).c0);	\
  (t).c11 = -(u).c1 * conj((v).c1);	\
  (t).c12 = -(u).c1 * conj((v).c2);	\
  (t).c20 = -(u).c2 * conj((v).c0);	\
  (t).c21 = -(u).c2 * conj((v).c1);	\
  (t).c22 = -(u).c2 * conj((v).c2);


#define _vector_tensor_vector_add(t, u, v, w, z) \
  (t).c00 = (u).c0 * conj((v).c0) + (w).c0 * conj((z).c0) ;	\
  (t).c01 = (u).c0 * conj((v).c1) + (w).c0 * conj((z).c1);	\
  (t).c02 = (u).c0 * conj((v).c2) + (w).c0 * conj((z).c2);	\
  (t).c10 = (u).c1 * conj((v).c0) + (w).c1 * conj((z).c0);	\
  (t).c11 = (u).c1 * conj((v).c1) + (w).c1 * conj((z).c1);	\
  (t).c12 = (u).c1 * conj((v).c2) + (w).c1 * conj((z).c2);	\
  (t).c20 = (u).c2 * conj((v).c0) + (w).c2 * conj((z).c0);	\
  (t).c21 = (u).c2 * conj((v).c1) + (w).c2 * conj((z).c1);	\
  (t).c22 = (u).c2 * conj((v).c2) + (w).c2 * conj((z).c2);



#endif

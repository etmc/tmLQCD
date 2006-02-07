/* $Id$ */
#ifndef _SU3_H
#define _SU3_H

/*******************************************************************************
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
*
*
*******************************************************************************/

#ifdef _STD_C99_COMPLEX_CHECKED
#include <complex.h>
#else
#include "complex.h"
#endif
#if (defined XLC && defined BGL)
# include "bgl.h"
#endif

typedef struct 
{
   complex c00,c01,c02,c10,c11,c12,c20,c21,c22;
} su3;

typedef struct
{
   complex c0,c1,c2;
} su3_vector;

typedef struct
{
   su3_vector s0,s1,s2,s3;
} spinor;

typedef struct
{
   spinor sp_up,sp_dn;
} bispinor;

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
#ifdef _STD_C99_COMPLEX
#define _vector_null(r) \
  (r).c0 = 0.;		\
  (r).c1 = 0.;		\
  (r).c2 = 0.;

#else 
#define _vector_null(r) \
   (r).c0.re=0.0; \
   (r).c0.im=0.0; \
   (r).c1.re=0.0; \
   (r).c1.im=0.0; \
   (r).c2.re=0.0; \
   (r).c2.im=0.0;
#endif

/* M. Hasenbusch Mon Sep 24
* r.c1=s.c1
* r.c2=s.c2
* r.c3=s.c3
*/
#ifdef _STD_C99_COMPLEX_CHECKED
#define _vector_assign(r,s)			\
  (r).c0 = (s).c0;				\
  (r).c1 = (s).c1;				\
  (r).c2 = (s).c2;

#else
#define _vector_assign(r,s) \
   (r).c0.re=(s).c0.re; \
   (r).c0.im=(s).c0.im; \
   (r).c1.re=(s).c1.re; \
   (r).c1.im=(s).c1.im; \
   (r).c2.re=(s).c2.re; \
   (r).c2.im=(s).c2.im;

#endif

/* M. Hasenbusch Mon Sep 24
* r.c1=-s.c1
* r.c2=-s.c2
* r.c3=-s.c3
*/
#ifdef _STD_C99_COMPLEX
#define _vector_minus_assign(r,s)		\
  (r).c0 = -(s).c0;				\
  (r).c1 = -(s).c1;				\
  (r).c2 = -(s).c2;

#else
#define _vector_minus_assign(r,s) \
   (r).c0.re=-(s).c0.re; \
   (r).c0.im=-(s).c0.im; \
   (r).c1.re=-(s).c1.re; \
   (r).c1.im=-(s).c1.im; \
   (r).c2.re=-(s).c2.re; \
   (r).c2.im=-(s).c2.im;

#endif

/*
* r.c1=c*s.c1 (c real)
* r.c2=c*s.c2
* r.c3=c*s.c3
*/
#ifdef _STD_C99_COMPLEX_CHECKED
#define _vector_mul(r,c,s)			\
  (r).c0 = (c)*(s).c0;				\
  (r).c1 = (c)*(s).c1;				\
  (r).c2 = (c)*(s).c2;

#else
#define _vector_mul(r,c,s) \
   (r).c0.re=(c)*(s).c0.re; \
   (r).c0.im=(c)*(s).c0.im; \
   (r).c1.re=(c)*(s).c1.re; \
   (r).c1.im=(c)*(s).c1.im; \
   (r).c2.re=(c)*(s).c2.re; \
   (r).c2.im=(c)*(s).c2.im;

#endif

#ifdef _STD_C99_COMPLEX
#define _vector_add_mul(r,c,s)			\
  (r).c0 += ((c)*(s).c0);			\
  (r).c1 += ((c)*(s).c1);			\
  (r).c2 += ((c)*(s).c2);

#elif defined XLC3
#define _vector_add_mul(r,c,s)			\
  (r).c0 = __fmadd(c, (s).c0, (r).c0);		\
  (r).c1 = __fmadd(c, (s).c1, (r).c1);		\
  (r).c2 = __fmadd(c, (s).c2, (r).c2); 

#else
#define _vector_add_mul(r,c,s) \
   (r).c0.re+=(c)*(s).c0.re; \
   (r).c0.im+=(c)*(s).c0.im; \
   (r).c1.re+=(c)*(s).c1.re; \
   (r).c1.im+=(c)*(s).c1.im; \
   (r).c2.re+=(c)*(s).c2.re; \
   (r).c2.im+=(c)*(s).c2.im;

#endif

/*
* r.c1=s1.c1+s2.c1
* r.c2=s1.c2+s2.c2
* r.c3=s1.c3+s2.c3
*/

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

#elif (defined XLC && defined BGL)

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

#elif defined _STD_C99_COMPLEX_CHECKED

#define _vector_add(r,s1,s2)			\
  (r).c0 = (s1).c0 + (s2).c0;			\
  (r).c1 = (s1).c1 + (s2).c1;			\
  (r).c2 = (s1).c2 + (s2).c2;

#define _vector_sub(r,s1,s2)			\
  (r).c0 = (s1).c0 - (s2).c0;			\
  (r).c1 = (s1).c1 - (s2).c1;			\
  (r).c2 = (s1).c2 - (s2).c2;

#else

#define _vector_add(r,s1,s2) \
   (r).c0.re=(s1).c0.re+(s2).c0.re; \
   (r).c0.im=(s1).c0.im+(s2).c0.im; \
   (r).c1.re=(s1).c1.re+(s2).c1.re; \
   (r).c1.im=(s1).c1.im+(s2).c1.im; \
   (r).c2.re=(s1).c2.re+(s2).c2.re; \
   (r).c2.im=(s1).c2.im+(s2).c2.im;

/*
* r.c1=s1.c1-s2.c1
* r.c2=s1.c2-s2.c2
* r.c3=s1.c3-s2.c3
*/

#define _vector_sub(r,s1,s2) \
   (r).c0.re=(s1).c0.re-(s2).c0.re; \
   (r).c0.im=(s1).c0.im-(s2).c0.im; \
   (r).c1.re=(s1).c1.re-(s2).c1.re; \
   (r).c1.im=(s1).c1.im-(s2).c1.im; \
   (r).c2.re=(s1).c2.re-(s2).c2.re; \
   (r).c2.im=(s1).c2.im-(s2).c2.im;

#endif


/*
* r.c1=s1.c1+i*s2.c1
* r.c2=s1.c2+i*s2.c2
* r.c3=s1.c3+i*s2.c3
*/

#ifdef _STD_C99_COMPLEX_CHECKED
#define _vector_i_add(r,s1,s2)			\
  (r).c0 = (s1).c0 + I*(s2).c0;			\
  (r).c1 = (s1).c1 + I*(s2).c1;			\
  (r).c2 = (s1).c2 + I*(s2).c2;

#else
#define _vector_i_add(r,s1,s2) \
   (r).c0.re=(s1).c0.re-(s2).c0.im; \
   (r).c0.im=(s1).c0.im+(s2).c0.re; \
   (r).c1.re=(s1).c1.re-(s2).c1.im; \
   (r).c1.im=(s1).c1.im+(s2).c1.re; \
   (r).c2.re=(s1).c2.re-(s2).c2.im; \
   (r).c2.im=(s1).c2.im+(s2).c2.re;

#endif

/*
* r.c1=s1.c1-i*s2.c1
* r.c2=s1.c2-i*s2.c2
* r.c3=s1.c3-i*s2.c3
*/

#ifdef _STD_C99_COMPLEX_CHECKED
#define _vector_i_sub(r,s1,s2)			\
  (r).c0 = (s1).c0 - I*(s2).c0;			\
  (r).c1 = (s1).c1 - I*(s2).c1;			\
  (r).c2 = (s1).c2 - I*(s2).c2;

#else
#define _vector_i_sub(r,s1,s2)	    \
   (r).c0.re=(s1).c0.re+(s2).c0.im; \
   (r).c0.im=(s1).c0.im-(s2).c0.re; \
   (r).c1.re=(s1).c1.re+(s2).c1.im; \
   (r).c1.im=(s1).c1.im-(s2).c1.re; \
   (r).c2.re=(s1).c2.re+(s2).c2.im; \
   (r).c2.im=(s1).c2.im-(s2).c2.re;

#endif

#ifdef _STD_C99_COMPLEX
#define _vector_combined_add_i_add(r1, s1, r2, s2, s)	\
  (r1).c0 = (s1).c0 + (s).c0;				\
  (r2).c0 = (s2).c0 + I*(s).c0;				\
  (r1).c1 = (s1).c1 + (s).c1;				\
  (r2).c1 = (s2).c1 + I*(s).c1;				\
  (r1).c2 = (s1).c2 + (s).c2;				\
  (r2).c2 = (s2).c2 + I*(s).c2;				\

#else
#define _vector_combined_add_i_add(r1, s1, r2, s2, s) \
   (r1).c0.re=(s1).c0.re+(s).c0.re; \
   (r1).c0.im=(s1).c0.im+(s).c0.im; \
   (r2).c0.re=(s2).c0.re-(s).c0.im; \
   (r2).c0.im=(s2).c0.im+(s).c0.re; \
   (r1).c1.re=(s1).c1.re+(s).c1.re; \
   (r1).c1.im=(s1).c1.im+(s).c1.im; \
   (r2).c1.re=(s2).c1.re-(s).c1.im; \
   (r2).c1.im=(s2).c1.im+(s).c1.re; \
   (r1).c2.re=(s1).c2.re+(s).c2.re; \
   (r1).c2.im=(s1).c2.im+(s).c2.im; \
   (r2).c2.re=(s2).c2.re-(s).c2.im; \
   (r2).c2.im=(s2).c2.im+(s).c2.re;

#endif

/*
* r.c1+=s.c1
* r.c2+=s.c2
* r.c3+=s.c3
*/

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

#elif defined _STD_C99_COMPLEX_CHECKED

#define _vector_add_assign(r,s)			\
  (r).c0 += (s).c0;				\
  (r).c1 += (s).c1;				\
  (r).c2 += (s).c2;

#define _vector_sub_assign(r,s)			\
  (r).c0 -= (s).c0;				\
  (r).c1 -= (s).c1;				\
  (r).c2 -= (s).c2;

#else

#define _vector_add_assign(r,s) \
   (r).c0.re+=(s).c0.re; \
   (r).c0.im+=(s).c0.im; \
   (r).c1.re+=(s).c1.re; \
   (r).c1.im+=(s).c1.im; \
   (r).c2.re+=(s).c2.re; \
   (r).c2.im+=(s).c2.im;

/*
* r.c1-=s.c1
* r.c2-=s.c2
* r.c3-=s.c3
*/

#define _vector_sub_assign(r,s) \
   (r).c0.re-=(s).c0.re; \
   (r).c0.im-=(s).c0.im; \
   (r).c1.re-=(s).c1.re; \
   (r).c1.im-=(s).c1.im; \
   (r).c2.re-=(s).c2.re; \
   (r).c2.im-=(s).c2.im;

#endif 



/*
* r.c1+=i*s.c1
* r.c2+=i*s.c2
* r.c3+=i*s.c3
*/

#ifdef _STD_C99_COMPLEX_CHECKED
#define _vector_i_add_assign(r,s)		\
  (r).c0 += (I*(s).c0);				\
  (r).c1 += (I*(s).c1);				\
  (r).c2 += (I*(s).c2);

#else
#define _vector_i_add_assign(r,s) \
   (r).c0.re-=(s).c0.im; \
   (r).c0.im+=(s).c0.re; \
   (r).c1.re-=(s).c1.im; \
   (r).c1.im+=(s).c1.re; \
   (r).c2.re-=(s).c2.im; \
   (r).c2.im+=(s).c2.re;

#endif

/*
* r.c1-=i*s.c1
* r.c2-=i*s.c2
* r.c3-=i*s.c3
*/

#ifdef _STD_C99_COMPLEX_CHECKED
#define _vector_i_sub_assign(r,s)		\
  (r).c0 -= (I*(s).c0);				\
  (r).c1 -= (I*(s).c1);				\
  (r).c2 -= (I*(s).c2);

#else
#define _vector_i_sub_assign(r,s) \
   (r).c0.re+=(s).c0.im; \
   (r).c0.im-=(s).c0.re; \
   (r).c1.re+=(s).c1.im; \
   (r).c1.im-=(s).c1.re; \
   (r).c2.re+=(s).c2.im; \
   (r).c2.im-=(s).c2.re;

#endif

/* M.Hasenbusch 
* r.c1=c*s.c1
* r.c2=c*s.c2
* r.c3=c*s.c3
*
* c complex
*/

#ifdef _STD_C99_COMPLEX_CHECKED
#define _complex_times_vector(r,c,s)		\
  (r).c0 = (c)*(s).c0;				\
  (r).c1 = (c)*(s).c1;				\
  (r).c2 = (c)*(s).c2;

#else
#define _complex_times_vector(r,c,s) \
   (r).c0.re=(c).re*(s).c0.re-(c).im*(s).c0.im; \
   (r).c0.im=(c).re*(s).c0.im+(c).im*(s).c0.re; \
   (r).c1.re=(c).re*(s).c1.re-(c).im*(s).c1.im; \
   (r).c1.im=(c).re*(s).c1.im+(c).im*(s).c1.re; \
   (r).c2.re=(c).re*(s).c2.re-(c).im*(s).c2.im; \
   (r).c2.im=(c).re*(s).c2.im+(c).im*(s).c2.re;

#endif

/* M.Hasenbusch */
#ifdef _STD_C99_COMPLEX_CHECKED
#define _complexcjg_times_vector(r,c,s)		\
  (r).c0 = conj(c)*(s).c0;			\
  (r).c1 = conj(c)*(s).c1;			\
  (r).c2 = conj(c)*(s).c2;

#else
#define _complexcjg_times_vector(r,c,s) \
   (r).c0.re=(c).re*(s).c0.re+(c).im*(s).c0.im; \
   (r).c0.im=(c).re*(s).c0.im-(c).im*(s).c0.re; \
   (r).c1.re=(c).re*(s).c1.re+(c).im*(s).c1.im; \
   (r).c1.im=(c).re*(s).c1.im-(c).im*(s).c1.re; \
   (r).c2.re=(c).re*(s).c2.re+(c).im*(s).c2.im; \
   (r).c2.im=(c).re*(s).c2.im-(c).im*(s).c2.re;

#endif

/*
* Real part of the scalar product (r,s)
*/

#ifdef _STD_C99_COMPLEX_CHECKED
#define _vector_prod_re(r,s)						\
  creal(conj((r).c0)*(s).c0 + conj((r).c1)*(s).c1 + conj((r).c2)*(s).c2);

#else
#define _vector_prod_re(r,s)		    \
  (r).c0.re*(s).c0.re+(r).c0.im*(s).c0.im+  \
  (r).c1.re*(s).c1.re+(r).c1.im*(s).c1.im+  \
  (r).c2.re*(s).c2.re+(r).c2.im*(s).c2.im;

#endif
/*
 * Imaginary part of the scalar product (r,s)
 */
#ifdef _STD_C99_COMPLEX_CHECKED
#define _vector_prod_im(r,s)			\
  cimag(conj((r).c0)*(s).c0 + conj((r).c1)*(s).c1 + conj((r).c2)*(s).c2);

#else
#define _vector_prod_im(r,s)		    \
  (r).c0.re*(s).c0.im-(r).c0.im*(s).c0.re+  \
  (r).c1.re*(s).c1.im-(r).c1.im*(s).c1.re+  \
  (r).c2.re*(s).c2.im-(r).c2.im*(s).c2.re; 
#endif

/*
* r.c1-=z*s.c1 (z of type complex)
* r.c2-=z*s.c2
* r.c3-=z*s.c3
*/

#ifdef _STD_C99_COMPLEX_CHECKED
#define _vector_project(r,z,s)			\
  (r).c0 -= (z*(s).c0);				\
  (r).c1 -= (z*(s).c1);				\
  (r).c2 -= (z*(s).c2);
#else
#define _vector_project(r,z,s) \
   (r).c0.re-=((z).re*(s).c0.re-(z).im*(s).c0.im); \
   (r).c0.im-=((z).re*(s).c0.im+(z).im*(s).c0.re); \
   (r).c1.re-=((z).re*(s).c1.re-(z).im*(s).c1.im); \
   (r).c1.im-=((z).re*(s).c1.im+(z).im*(s).c1.re); \
   (r).c2.re-=((z).re*(s).c2.re-(z).im*(s).c2.im); \
   (r).c2.im-=((z).re*(s).c2.im+(z).im*(s).c2.re);
#endif

/*
* SU(3) matrix u times SU(3) vector s
*  
* r.c1=(u*s).c1
* r.c2=(u*s).c2
* r.c3=(u*s).c3
*/

#if ((defined SSE2) || (defined SSE3))

#define _su3_multiply(r,u,s) \
_sse_load(s); \
_sse_su3_multiply(u); \
_sse_store_up(r);

#define _su3_inverse_multiply(r,u,s) \
_sse_load(s); \
_sse_su3_inverse_multiply(u); \
_sse_store_up(r);

#elif (defined XLC && defined BGL)

#define _su3_multiply(r,u,s) \
  _bgl_load(s); \
  _bgl_su3_multiply(u); \
  _bgl_store_up(r);

#define _su3_inverse_multiply(r,u,s)		\
  _bgl_load(s); \
  _bgl_su3_inverse_multiply(u); \
  _bgl_store_up(r);

#elif defined _STD_C99_COMPLEX_CHECKED

#define _su3_multiply(r,u,s)					\
  (r).c0 = (u).c00*(s).c0 + (u).c01*(s).c1 + (u).c02*(s).c2;	\
  (r).c1 = (u).c10*(s).c0 + (u).c11*(s).c1 + (u).c12*(s).c2;	\
  (r).c2 = (u).c20*(s).c0 + (u).c21*(s).c1 + (u).c22*(s).c2;  

#define _su3_inverse_multiply(r,u,s)		\
(r).c0 = conj((u).c00)*(s).c0 + conj((u).c10)*(s).c1 + conj((u).c20)*(s).c2;	\
(r).c1 = conj((u).c01)*(s).c0 + conj((u).c11)*(s).c1 + conj((u).c21)*(s).c2;	\
(r).c2 = conj((u).c02)*(s).c0 + conj((u).c12)*(s).c1 + conj((u).c22)*(s).c2;   

#else

#define _su3_multiply(r,u,s) \
   (r).c0.re= (u).c00.re*(s).c0.re-(u).c00.im*(s).c0.im  \
             +(u).c01.re*(s).c1.re-(u).c01.im*(s).c1.im  \
             +(u).c02.re*(s).c2.re-(u).c02.im*(s).c2.im; \
   (r).c0.im= (u).c00.re*(s).c0.im+(u).c00.im*(s).c0.re  \
             +(u).c01.re*(s).c1.im+(u).c01.im*(s).c1.re  \
             +(u).c02.re*(s).c2.im+(u).c02.im*(s).c2.re; \
   (r).c1.re= (u).c10.re*(s).c0.re-(u).c10.im*(s).c0.im  \
             +(u).c11.re*(s).c1.re-(u).c11.im*(s).c1.im  \
             +(u).c12.re*(s).c2.re-(u).c12.im*(s).c2.im; \
   (r).c1.im= (u).c10.re*(s).c0.im+(u).c10.im*(s).c0.re  \
             +(u).c11.re*(s).c1.im+(u).c11.im*(s).c1.re  \
             +(u).c12.re*(s).c2.im+(u).c12.im*(s).c2.re; \
   (r).c2.re= (u).c20.re*(s).c0.re-(u).c20.im*(s).c0.im  \
             +(u).c21.re*(s).c1.re-(u).c21.im*(s).c1.im  \
             +(u).c22.re*(s).c2.re-(u).c22.im*(s).c2.im; \
   (r).c2.im= (u).c20.re*(s).c0.im+(u).c20.im*(s).c0.re  \
             +(u).c21.re*(s).c1.im+(u).c21.im*(s).c1.re  \
             +(u).c22.re*(s).c2.im+(u).c22.im*(s).c2.re;

/*
* SU(3) matrix u^dagger times SU(3) vector s
*  
* r.c1=(u^dagger*s).c1
* r.c2=(u^dagger*s).c2
* r.c3=(u^dagger*s).c3
*/

#define _su3_inverse_multiply(r,u,s) \
   (r).c0.re= (u).c00.re*(s).c0.re+(u).c00.im*(s).c0.im  \
             +(u).c10.re*(s).c1.re+(u).c10.im*(s).c1.im  \
             +(u).c20.re*(s).c2.re+(u).c20.im*(s).c2.im; \
   (r).c0.im= (u).c00.re*(s).c0.im-(u).c00.im*(s).c0.re  \
             +(u).c10.re*(s).c1.im-(u).c10.im*(s).c1.re  \
             +(u).c20.re*(s).c2.im-(u).c20.im*(s).c2.re; \
   (r).c1.re= (u).c01.re*(s).c0.re+(u).c01.im*(s).c0.im  \
             +(u).c11.re*(s).c1.re+(u).c11.im*(s).c1.im  \
             +(u).c21.re*(s).c2.re+(u).c21.im*(s).c2.im; \
   (r).c1.im= (u).c01.re*(s).c0.im-(u).c01.im*(s).c0.re  \
             +(u).c11.re*(s).c1.im-(u).c11.im*(s).c1.re  \
             +(u).c21.re*(s).c2.im-(u).c21.im*(s).c2.re; \
   (r).c2.re= (u).c02.re*(s).c0.re+(u).c02.im*(s).c0.im  \
             +(u).c12.re*(s).c1.re+(u).c12.im*(s).c1.im  \
             +(u).c22.re*(s).c2.re+(u).c22.im*(s).c2.im; \
   (r).c2.im= (u).c02.re*(s).c0.im-(u).c02.im*(s).c0.re  \
             +(u).c12.re*(s).c1.im-(u).c12.im*(s).c1.re  \
             +(u).c22.re*(s).c2.im-(u).c22.im*(s).c2.re;
#endif



/*******************************************************************************
*
* Macros for SU(3) matrices
*
* Arguments are variables of type su3
*
*******************************************************************************/

/* M. Hasenbusch */
#ifdef _STD_C99_COMPLEX
#define _su3_norm_sq(x,u)			\
  x = (u).c00*conj((u).c00) + (u).c01*conj((u).c01) + (u).c02*conj((u).c02) \
    (u).c10*conj((u).c10) + (u).c11*conj((u).c11) + (u).c12*conj((u).c12) \
    (u).c20*conj((u).c20) + (u).c21*conj((u).c21) + (u).c22*conj((u).c22); 

#else

#define _su3_norm_sq(x,u) \
x = (u).c00.re*(u).c00.re + (u).c00.im*(u).c00.im \
   +(u).c01.re*(u).c01.re + (u).c01.im*(u).c01.im \
   +(u).c02.re*(u).c02.re + (u).c02.im*(u).c02.im \
   +(u).c10.re*(u).c10.re + (u).c10.im*(u).c10.im \
   +(u).c11.re*(u).c11.re + (u).c11.im*(u).c11.im \
   +(u).c12.re*(u).c12.re + (u).c12.im*(u).c12.im \
   +(u).c20.re*(u).c20.re + (u).c20.im*(u).c20.im \
   +(u).c21.re*(u).c21.re + (u).c21.im*(u).c21.im \
   +(u).c22.re*(u).c22.re + (u).c22.im*(u).c22.im; 

#endif

/*
 u=1 
 added by M.Hasenbusch Thu Aug  9 10:27:28 MEST 2001 */

#ifdef _STD_C99_COMPLEX
#define _su2_one(u)				\
  (u).c00 = 1.;					\
  (u).c01 = 0.;					\
  (u).c02 = 0.;					\
  (u).c10 = 0.;					\
  (u).c11 = 1.;					\
  (u).c12 = 0.;					\
  (u).c20 = 0.;					\
  (u).c21 = 0.;					\
  (u).c22 = 1.;

#else
#define _su3_one(u) \
   (u).c00.re=1.0; \
   (u).c00.im=0.0; \
   (u).c01.re=0.0; \
   (u).c01.im=0.0; \
   (u).c02.re=0.0; \
   (u).c02.im=0.0; \
   (u).c10.re=0.0; \
   (u).c10.im=0.0; \
   (u).c11.re=1.0; \
   (u).c11.im=0.0; \
   (u).c12.re=0.0; \
   (u).c12.im=0.0; \
   (u).c20.re=0.0; \
   (u).c20.im=0.0; \
   (u).c21.re=0.0; \
   (u).c21.im=0.0; \
   (u).c22.re=1.0; \
   (u).c22.im=0.0;

#endif
/*
 u=0 
 added by M.Hasenbusch Thu Aug  9 10:27:28 MEST 2001 */
#ifdef _STD_C99_COMPLEX
#define _su2_zero(u)				\
  (u).c00 = 0.;					\
  (u).c01 = 0.;					\
  (u).c02 = 0.;					\
  (u).c10 = 0.;					\
  (u).c11 = 0.;					\
  (u).c12 = 0.;					\
  (u).c20 = 0.;					\
  (u).c21 = 0.;					\
  (u).c22 = 0.;
#else
#define _su3_zero(u) \
   (u).c00.re=0.0; \
   (u).c00.im=0.0; \
   (u).c01.re=0.0; \
   (u).c01.im=0.0; \
   (u).c02.re=0.0; \
   (u).c02.im=0.0; \
   (u).c10.re=0.0; \
   (u).c10.im=0.0; \
   (u).c11.re=0.0; \
   (u).c11.im=0.0; \
   (u).c12.re=0.0; \
   (u).c12.im=0.0; \
   (u).c20.re=0.0; \
   (u).c20.im=0.0; \
   (u).c21.re=0.0; \
   (u).c21.im=0.0; \
   (u).c22.re=0.0; \
   (u).c22.im=0.0;

#endif

/* M. Hasenbusch
* u=v
*/

#ifdef _STD_C99_COMPLEX
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
#else
#define _su3_assign(u,v) \
   (u).c00.re= (v).c00.re; \
   (u).c00.im= (v).c00.im; \
   (u).c01.re= (v).c01.re; \
   (u).c01.im= (v).c01.im; \
   (u).c02.re= (v).c02.re; \
   (u).c02.im= (v).c02.im; \
   (u).c10.re= (v).c10.re; \
   (u).c10.im= (v).c10.im; \
   (u).c11.re= (v).c11.re; \
   (u).c11.im= (v).c11.im; \
   (u).c12.re= (v).c12.re; \
   (u).c12.im= (v).c12.im; \
   (u).c20.re= (v).c20.re; \
   (u).c20.im= (v).c20.im; \
   (u).c21.re= (v).c21.re; \
   (u).c21.im= (v).c21.im; \
   (u).c22.re= (v).c22.re; \
   (u).c22.im= (v).c22.im;
#endif

/* M. Hasenbusch
* u=-v
*/
#ifdef _STD_C99_COMPLEX
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
#else
#define _su3_minus_assign(u,v) \
   (u).c00.re= -(v).c00.re; \
   (u).c00.im= -(v).c00.im; \
   (u).c01.re= -(v).c01.re; \
   (u).c01.im= -(v).c01.im; \
   (u).c02.re= -(v).c02.re; \
   (u).c02.im= -(v).c02.im; \
   (u).c10.re= -(v).c10.re; \
   (u).c10.im= -(v).c10.im; \
   (u).c11.re= -(v).c11.re; \
   (u).c11.im= -(v).c11.im; \
   (u).c12.re= -(v).c12.re; \
   (u).c12.im= -(v).c12.im; \
   (u).c20.re= -(v).c20.re; \
   (u).c20.im= -(v).c20.im; \
   (u).c21.re= -(v).c21.re; \
   (u).c21.im= -(v).c21.im; \
   (u).c22.re= -(v).c22.re; \
   (u).c22.im= -(v).c22.im;
#endif
/*
* u=v^dagger
*/

#ifdef _STD_C99_COMPLEX
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
#else
#define _su3_dagger(u,v) \
   (u).c00.re= (v).c00.re; \
   (u).c00.im=-(v).c00.im; \
   (u).c01.re= (v).c10.re; \
   (u).c01.im=-(v).c10.im; \
   (u).c02.re= (v).c20.re; \
   (u).c02.im=-(v).c20.im; \
   (u).c10.re= (v).c01.re; \
   (u).c10.im=-(v).c01.im; \
   (u).c11.re= (v).c11.re; \
   (u).c11.im=-(v).c11.im; \
   (u).c12.re= (v).c21.re; \
   (u).c12.im=-(v).c21.im; \
   (u).c20.re= (v).c02.re; \
   (u).c20.im=-(v).c02.im; \
   (u).c21.re= (v).c12.re; \
   (u).c21.im=-(v).c12.im; \
   (u).c22.re= (v).c22.re; \
   (u).c22.im=-(v).c22.im;
#endif

#ifdef _STD_C99_COMPLEX
#define _itimes_su3(u,v)			\
  (u).c00 = I*(v).c00;				\
  (u).c01 = I*(v).c01;				\
  (u).c02 = I*(v).c02;				\
  (u).c10 = I*(v).c10;				\
  (u).c11 = I*(v).c11;				\
  (u).c12 = I*(v).c12;				\
  (u).c20 = I*(v).c20;				\
  (u).c21 = I*(v).c21;				\
  (u).c22 = I*(v).c22;

#else
/* M.Hasenbusch */
#define _itimes_su3(u,v) \
   (u).c00.re=-(v).c00.im; \
   (u).c00.im= (v).c00.re; \
   (u).c01.re=-(v).c01.im; \
   (u).c01.im= (v).c01.re; \
   (u).c02.re=-(v).c02.im; \
   (u).c02.im= (v).c02.re; \
   (u).c10.re=-(v).c10.im; \
   (u).c10.im= (v).c10.re; \
   (u).c11.re=-(v).c11.im; \
   (u).c11.im= (v).c11.re; \
   (u).c12.re=-(v).c12.im; \
   (u).c12.im= (v).c12.re; \
   (u).c20.re=-(v).c20.im; \
   (u).c20.im= (v).c20.re; \
   (u).c21.re=-(v).c21.im; \
   (u).c21.im= (v).c21.re; \
   (u).c22.re=-(v).c22.im; \
   (u).c22.im= (v).c22.re;
#endif

/* M. Hasenbusch
* u=c*v 
* c real
*/

#ifdef _STD_C99_COMPLEX
#define _real_times_su3(u,a,v)				\
  (u).c00 = (a)*(v).c00;				\
  (u).c01 = (a)*(v).c01;				\
  (u).c02 = (a)*(v).c02;				\
  (u).c10 = (a)*(v).c10;				\
  (u).c11 = (a)*(v).c11;				\
  (u).c12 = (a)*(v).c12;				\
  (u).c20 = (a)*(v).c20;				\
  (u).c21 = (a)*(v).c21;				\
  (u).c22 = (a)*(v).c22;

#else 
#define _real_times_su3(u,a,v) \
   (u).c00.re= (a)*(v).c00.re; \
   (u).c00.im= (a)*(v).c00.im; \
   (u).c01.re= (a)*(v).c01.re; \
   (u).c01.im= (a)*(v).c01.im; \
   (u).c02.re= (a)*(v).c02.re; \
   (u).c02.im= (a)*(v).c02.im; \
   (u).c10.re= (a)*(v).c10.re; \
   (u).c10.im= (a)*(v).c10.im; \
   (u).c11.re= (a)*(v).c11.re; \
   (u).c11.im= (a)*(v).c11.im; \
   (u).c12.re= (a)*(v).c12.re; \
   (u).c12.im= (a)*(v).c12.im; \
   (u).c20.re= (a)*(v).c20.re; \
   (u).c20.im= (a)*(v).c20.im; \
   (u).c21.re= (a)*(v).c21.re; \
   (u).c21.im= (a)*(v).c21.im; \
   (u).c22.re= (a)*(v).c22.re; \
   (u).c22.im= (a)*(v).c22.im;
#endif

/* M. Hasenbusch
* u=v-w
*/

#define _su3_minus_su3(u,v,w) \
   (u).c00.re= (v).c00.re-(w).c00.re; \
   (u).c00.im= (v).c00.im-(w).c00.im; \
   (u).c01.re= (v).c01.re-(w).c01.re; \
   (u).c01.im= (v).c01.im-(w).c01.im; \
   (u).c02.re= (v).c02.re-(w).c02.re; \
   (u).c02.im= (v).c02.im-(w).c02.im; \
   (u).c10.re= (v).c10.re-(w).c10.re; \
   (u).c10.im= (v).c10.im-(w).c10.im; \
   (u).c11.re= (v).c11.re-(w).c11.re; \
   (u).c11.im= (v).c11.im-(w).c11.im; \
   (u).c12.re= (v).c12.re-(w).c12.re; \
   (u).c12.im= (v).c12.im-(w).c12.im; \
   (u).c20.re= (v).c20.re-(w).c20.re; \
   (u).c20.im= (v).c20.im-(w).c20.im; \
   (u).c21.re= (v).c21.re-(w).c21.re; \
   (u).c21.im= (v).c21.im-(w).c21.im; \
   (u).c22.re= (v).c22.re-(w).c22.re; \
   (u).c22.im= (v).c22.im-(w).c22.im;

/* M. Hasenbusch
* u=i*(v-w)
*/

#define _itimes_su3_minus_su3(u,v,w) \
   (u).c00.im= (v).c00.re-(w).c00.re; \
   (u).c00.re= (w).c00.im-(v).c00.im; \
   (u).c01.im= (v).c01.re-(w).c01.re; \
   (u).c01.re= (w).c01.im-(v).c01.im; \
   (u).c02.im= (v).c02.re-(w).c02.re; \
   (u).c02.re= (w).c02.im-(v).c02.im; \
   (u).c10.im= (v).c10.re-(w).c10.re; \
   (u).c10.re= (w).c10.im-(v).c10.im; \
   (u).c11.im= (v).c11.re-(w).c11.re; \
   (u).c11.re= (w).c11.im-(v).c11.im; \
   (u).c12.im= (v).c12.re-(w).c12.re; \
   (u).c12.re= (w).c12.im-(v).c12.im; \
   (u).c20.im= (v).c20.re-(w).c20.re; \
   (u).c20.re= (w).c20.im-(v).c20.im; \
   (u).c21.im= (v).c21.re-(w).c21.re; \
   (u).c21.re= (w).c21.im-(v).c21.im; \
   (u).c22.im= (v).c22.re-(w).c22.re; \
   (u).c22.re= (w).c22.im-(v).c22.im;

/* M. Hasenbusch
* u=v+w
*/

#define _su3_plus_su3(u,v,w) \
   (u).c00.re= (v).c00.re+(w).c00.re; \
   (u).c00.im= (v).c00.im+(w).c00.im; \
   (u).c01.re= (v).c01.re+(w).c01.re; \
   (u).c01.im= (v).c01.im+(w).c01.im; \
   (u).c02.re= (v).c02.re+(w).c02.re; \
   (u).c02.im= (v).c02.im+(w).c02.im; \
   (u).c10.re= (v).c10.re+(w).c10.re; \
   (u).c10.im= (v).c10.im+(w).c10.im; \
   (u).c11.re= (v).c11.re+(w).c11.re; \
   (u).c11.im= (v).c11.im+(w).c11.im; \
   (u).c12.re= (v).c12.re+(w).c12.re; \
   (u).c12.im= (v).c12.im+(w).c12.im; \
   (u).c20.re= (v).c20.re+(w).c20.re; \
   (u).c20.im= (v).c20.im+(w).c20.im; \
   (u).c21.re= (v).c21.re+(w).c21.re; \
   (u).c21.im= (v).c21.im+(w).c21.im; \
   (u).c22.re= (v).c22.re+(w).c22.re; \
   (u).c22.im= (v).c22.im+(w).c22.im;

/* M. Hasenbusch
* u=-(v+w)
*/

#define _minus_su3_plus_su3(u,v,w) \
   (u).c00.re=-(v).c00.re-(w).c00.re; \
   (u).c00.im=-(v).c00.im-(w).c00.im; \
   (u).c01.re=-(v).c01.re-(w).c01.re; \
   (u).c01.im=-(v).c01.im-(w).c01.im; \
   (u).c02.re=-(v).c02.re-(w).c02.re; \
   (u).c02.im=-(v).c02.im-(w).c02.im; \
   (u).c10.re=-(v).c10.re-(w).c10.re; \
   (u).c10.im=-(v).c10.im-(w).c10.im; \
   (u).c11.re=-(v).c11.re-(w).c11.re; \
   (u).c11.im=-(v).c11.im-(w).c11.im; \
   (u).c12.re=-(v).c12.re-(w).c12.re; \
   (u).c12.im=-(v).c12.im-(w).c12.im; \
   (u).c20.re=-(v).c20.re-(w).c20.re; \
   (u).c20.im=-(v).c20.im-(w).c20.im; \
   (u).c21.re=-(v).c21.re-(w).c21.re; \
   (u).c21.im=-(v).c21.im-(w).c21.im; \
   (u).c22.re=-(v).c22.re-(w).c22.re; \
   (u).c22.im=-(v).c22.im-(w).c22.im;

/* M. Hasenbusch
* u=i*(v+w)
*/

#define _itimes_su3_plus_su3(u,v,w) \
   (u).c00.im= (v).c00.re+(w).c00.re; \
   (u).c00.re=-(v).c00.im-(w).c00.im; \
   (u).c01.im= (v).c01.re+(w).c01.re; \
   (u).c01.re=-(v).c01.im-(w).c01.im; \
   (u).c02.im= (v).c02.re+(w).c02.re; \
   (u).c02.re=-(v).c02.im-(w).c02.im; \
   (u).c10.im= (v).c10.re+(w).c10.re; \
   (u).c10.re=-(v).c10.im-(w).c10.im; \
   (u).c11.im= (v).c11.re+(w).c11.re; \
   (u).c11.re=-(v).c11.im-(w).c11.im; \
   (u).c12.im= (v).c12.re+(w).c12.re; \
   (u).c12.re=-(v).c12.im-(w).c12.im; \
   (u).c20.im= (v).c20.re+(w).c20.re; \
   (u).c20.re=-(v).c20.im-(w).c20.im; \
   (u).c21.im= (v).c21.re+(w).c21.re; \
   (u).c21.re=-(v).c21.im-(w).c21.im; \
   (u).c22.im= (v).c22.re+(w).c22.re; \
   (u).c22.re=-(v).c22.im-(w).c22.im;

/* M. Hasenbusch
* u=-i*(v+w)
*/

#define _minus_itimes_su3_plus_su3(u,v,w) \
   (u).c00.im=-(v).c00.re-(w).c00.re; \
   (u).c00.re= (v).c00.im+(w).c00.im; \
   (u).c01.im=-(v).c01.re-(w).c01.re; \
   (u).c01.re= (v).c01.im+(w).c01.im; \
   (u).c02.im=-(v).c02.re-(w).c02.re; \
   (u).c02.re= (v).c02.im+(w).c02.im; \
   (u).c10.im=-(v).c10.re-(w).c10.re; \
   (u).c10.re= (v).c10.im+(w).c10.im; \
   (u).c11.im=-(v).c11.re-(w).c11.re; \
   (u).c11.re= (v).c11.im+(w).c11.im; \
   (u).c12.im=-(v).c12.re-(w).c12.re; \
   (u).c12.re= (v).c12.im+(w).c12.im; \
   (u).c20.im=-(v).c20.re-(w).c20.re; \
   (u).c20.re= (v).c20.im+(w).c20.im; \
   (u).c21.im=-(v).c21.re-(w).c21.re; \
   (u).c21.re= (v).c21.im+(w).c21.im; \
   (u).c22.im=-(v).c22.re-(w).c22.re; \
   (u).c22.re= (v).c22.im+(w).c22.im;

/* M.Hasenbusch */
#define _complex_times_su3(r,c,s) \
   (r).c00.re=(c).re*(s).c00.re-(c).im*(s).c00.im; \
   (r).c00.im=(c).re*(s).c00.im+(c).im*(s).c00.re; \
   (r).c01.re=(c).re*(s).c01.re-(c).im*(s).c01.im; \
   (r).c01.im=(c).re*(s).c01.im+(c).im*(s).c01.re; \
   (r).c02.re=(c).re*(s).c02.re-(c).im*(s).c02.im; \
   (r).c02.im=(c).re*(s).c02.im+(c).im*(s).c02.re; \
   (r).c10.re=(c).re*(s).c10.re-(c).im*(s).c10.im; \
   (r).c10.im=(c).re*(s).c10.im+(c).im*(s).c10.re; \
   (r).c11.re=(c).re*(s).c11.re-(c).im*(s).c11.im; \
   (r).c11.im=(c).re*(s).c11.im+(c).im*(s).c11.re; \
   (r).c12.re=(c).re*(s).c12.re-(c).im*(s).c12.im; \
   (r).c12.im=(c).re*(s).c12.im+(c).im*(s).c12.re; \
   (r).c20.re=(c).re*(s).c20.re-(c).im*(s).c20.im; \
   (r).c20.im=(c).re*(s).c20.im+(c).im*(s).c20.re; \
   (r).c21.re=(c).re*(s).c21.re-(c).im*(s).c21.im; \
   (r).c21.im=(c).re*(s).c21.im+(c).im*(s).c21.re; \
   (r).c22.re=(c).re*(s).c22.re-(c).im*(s).c22.im; \
   (r).c22.im=(c).re*(s).c22.im+(c).im*(s).c22.re; 

/* M.Hasenbusch */
#define _complexcjg_times_su3(r,c,s) \
   (r).c00.re=(c).re*(s).c00.re+(c).im*(s).c00.im; \
   (r).c00.im=(c).re*(s).c00.im-(c).im*(s).c00.re; \
   (r).c01.re=(c).re*(s).c01.re+(c).im*(s).c01.im; \
   (r).c01.im=(c).re*(s).c01.im-(c).im*(s).c01.re; \
   (r).c02.re=(c).re*(s).c02.re+(c).im*(s).c02.im; \
   (r).c02.im=(c).re*(s).c02.im-(c).im*(s).c02.re; \
   (r).c10.re=(c).re*(s).c10.re+(c).im*(s).c10.im; \
   (r).c10.im=(c).re*(s).c10.im-(c).im*(s).c10.re; \
   (r).c11.re=(c).re*(s).c11.re+(c).im*(s).c11.im; \
   (r).c11.im=(c).re*(s).c11.im-(c).im*(s).c11.re; \
   (r).c12.re=(c).re*(s).c12.re+(c).im*(s).c12.im; \
   (r).c12.im=(c).re*(s).c12.im-(c).im*(s).c12.re; \
   (r).c20.re=(c).re*(s).c20.re+(c).im*(s).c20.im; \
   (r).c20.im=(c).re*(s).c20.im-(c).im*(s).c20.re; \
   (r).c21.re=(c).re*(s).c21.re+(c).im*(s).c21.im; \
   (r).c21.im=(c).re*(s).c21.im-(c).im*(s).c21.re; \
   (r).c22.re=(c).re*(s).c22.re+(c).im*(s).c22.im; \
   (r).c22.im=(c).re*(s).c22.im-(c).im*(s).c22.re;


/* M. Hasenbusch
* su3_acc
*/

#if ((defined SSE2) || (defined SSE3))
#define _su3_acc(u,v) _sse_su3_acc(u,v) 
#else
#define _su3_acc(u,v) \
   (u).c00.re+=(v).c00.re; \
   (u).c00.im+=(v).c00.im; \
   (u).c01.re+=(v).c01.re; \
   (u).c01.im+=(v).c01.im; \
   (u).c02.re+=(v).c02.re; \
   (u).c02.im+=(v).c02.im; \
   (u).c10.re+=(v).c10.re; \
   (u).c10.im+=(v).c10.im; \
   (u).c11.re+=(v).c11.re; \
   (u).c11.im+=(v).c11.im; \
   (u).c12.re+=(v).c12.re; \
   (u).c12.im+=(v).c12.im; \
   (u).c20.re+=(v).c20.re; \
   (u).c20.im+=(v).c20.im; \
   (u).c21.re+=(v).c21.re; \
   (u).c21.im+=(v).c21.im; \
   (u).c22.re+=(v).c22.re; \
   (u).c22.im+=(v).c22.im;
#endif

/*
* su3_refac_acc
*/

#define _su3_refac_acc(u,a,v) \
   (u).c00.re+=a*(v).c00.re; \
   (u).c00.im+=a*(v).c00.im; \
   (u).c01.re+=a*(v).c01.re; \
   (u).c01.im+=a*(v).c01.im; \
   (u).c02.re+=a*(v).c02.re; \
   (u).c02.im+=a*(v).c02.im; \
   (u).c10.re+=a*(v).c10.re; \
   (u).c10.im+=a*(v).c10.im; \
   (u).c11.re+=a*(v).c11.re; \
   (u).c11.im+=a*(v).c11.im; \
   (u).c12.re+=a*(v).c12.re; \
   (u).c12.im+=a*(v).c12.im; \
   (u).c20.re+=a*(v).c20.re; \
   (u).c20.im+=a*(v).c20.im; \
   (u).c21.re+=a*(v).c21.re; \
   (u).c21.im+=a*(v).c21.im; \
   (u).c22.re+=a*(v).c22.re; \
   (u).c22.im+=a*(v).c22.im;
/*
* su3_imfac_acc
*/

#define _su3_imfac_acc(u,a,v) \
   (u).c00.re-=a*(v).c00.im; \
   (u).c00.im+=a*(v).c00.re; \
   (u).c01.re-=a*(v).c01.im; \
   (u).c01.im+=a*(v).c01.re; \
   (u).c02.re-=a*(v).c02.im; \
   (u).c02.im+=a*(v).c02.re; \
   (u).c10.re-=a*(v).c10.im; \
   (u).c10.im+=a*(v).c10.re; \
   (u).c11.re-=a*(v).c11.im; \
   (u).c11.im+=a*(v).c11.re; \
   (u).c12.re-=a*(v).c12.im; \
   (u).c12.im+=a*(v).c12.re; \
   (u).c20.re-=a*(v).c20.im; \
   (u).c20.im+=a*(v).c20.re; \
   (u).c21.re-=a*(v).c21.im; \
   (u).c21.im+=a*(v).c21.re; \
   (u).c22.re-=a*(v).c22.im; \
   (u).c22.im+=a*(v).c22.re;


#if ((defined SSE2) || (defined SSE3))
#define _su3_times_su3(u,v,w) _sse_su3_times_su3(u,v,w)
#define _su3_times_su3_acc(u,v,w) _sse_su3_times_su3_acc(u,v,w)
#define _su3d_times_su3(u,v,w) _sse_su3d_times_su3(u,v,w)
#define _su3d_times_su3_acc(u,v,w) _sse_su3d_times_su3_acc(u,v,w)
#define _su3_times_su3d(u,v,w) _sse_su3_times_su3d(u,v,w)
#define _su3_times_su3d_acc(u,v,w) _sse_su3_times_su3d_acc(u,v,w)  

#elif defined _STD_C99_COMPLEX_CHECKED
#define _su3_times_su3(u,v,w)					\
  (u).c00 = (v).c00*(w).c00 + (v).c01*(w).c10  + (v).c02*(w).c20;	\
  (u).c01 = (v).c00*(w).c01 + (v).c01*(w).c11  + (v).c02*(w).c21;	\
  (u).c02 = (v).c00*(w).c02 + (v).c01*(w).c12  + (v).c02*(w).c22;	\
  (u).c10 = (v).c10*(w).c00 + (v).c11*(w).c10	+ (v).c12*(w).c20;	\
  (u).c11 = (v).c10*(w).c01 + (v).c11*(w).c11	+ (v).c12*(w).c21;	\
  (u).c12 = (v).c10*(w).c02 + (v).c11*(w).c12	+ (v).c12*(w).c22;	\
  (u).c20 = (v).c20*(w).c00 + (v).c21*(w).c10	+ (v).c22*(w).c20;	\
  (u).c21 = (v).c20*(w).c01 + (v).c21*(w).c11	+ (v).c22*(w).c21;	\
  (u).c22 = (v).c20*(w).c02 + (v).c21*(w).c12	+ (v).c22*(w).c22;

/* C.Urbach u=u + v * w */
#define _su3_times_su3_acc(u,v,w)					\
  (u).c00 += (v).c00*(w).c00 + (v).c01*(w).c10  + (v).c02*(w).c20;	\
  (u).c01 += (v).c00*(w).c01 + (v).c01*(w).c11  + (v).c02*(w).c21;	\
  (u).c02 += (v).c00*(w).c02 + (v).c01*(w).c12  + (v).c02*(w).c22;	\
  (u).c10 += (v).c10*(w).c00 + (v).c11*(w).c10	+ (v).c12*(w).c20;	\
  (u).c11 += (v).c10*(w).c01 + (v).c11*(w).c11	+ (v).c12*(w).c21;	\
  (u).c12 += (v).c10*(w).c02 + (v).c11*(w).c12	+ (v).c12*(w).c22;	\
  (u).c20 += (v).c20*(w).c00 + (v).c21*(w).c10	+ (v).c22*(w).c20;	\
  (u).c21 += (v).c20*(w).c01 + (v).c21*(w).c11	+ (v).c22*(w).c21;	\
  (u).c22 += (v).c20*(w).c02 + (v).c21*(w).c12	+ (v).c22*(w).c22;

/* C.Urbach u=v * w^{dag} */
#define _su3_times_su3d(u,v,w)						\
  (u).c00 =  (v).c00*conj((w).c00) + (v).c01*conj((w).c01)  + (v).c02*conj((w).c02); \
  (u).c01 =  (v).c00*conj((w).c10) + (v).c01*conj((w).c11)  + (v).c02*conj((w).c12); \
  (u).c02 =  (v).c00*conj((w).c20) + (v).c01*conj((w).c21)  + (v).c02*conj((w).c22); \
  (u).c10 =  (v).c10*conj((w).c00) + (v).c11*conj((w).c01)  + (v).c12*conj((w).c02); \
  (u).c11 =  (v).c10*conj((w).c10) + (v).c11*conj((w).c11)  + (v).c12*conj((w).c12); \
  (u).c12 =  (v).c10*conj((w).c20) + (v).c11*conj((w).c21)  + (v).c12*conj((w).c22); \
  (u).c20 =  (v).c20*conj((w).c00) + (v).c21*conj((w).c01)  + (v).c22*conj((w).c02); \
  (u).c21 =  (v).c20*conj((w).c10) + (v).c21*conj((w).c11)  + (v).c22*conj((w).c12); \
  (u).c22 =  (v).c20*conj((w).c20) + (v).c21*conj((w).c21)  + (v).c22*conj((w).c22);

/* C.Urbach u=u + v * w^{dag} */
#define _su3_times_su3d_acc(u,v,w)						\
  (u).c00 +=  (v).c00*conj((w).c00) + (v).c01*conj((w).c01)  + (v).c02*conj((w).c02); \
  (u).c01 +=  (v).c00*conj((w).c10) + (v).c01*conj((w).c11)  + (v).c02*conj((w).c12); \
  (u).c02 +=  (v).c00*conj((w).c20) + (v).c01*conj((w).c21)  + (v).c02*conj((w).c22); \
  (u).c10 +=  (v).c10*conj((w).c00) + (v).c11*conj((w).c01)  + (v).c12*conj((w).c02); \
  (u).c11 +=  (v).c10*conj((w).c10) + (v).c11*conj((w).c11)  + (v).c12*conj((w).c12); \
  (u).c12 +=  (v).c10*conj((w).c20) + (v).c11*conj((w).c21)  + (v).c12*conj((w).c22); \
  (u).c20 +=  (v).c20*conj((w).c00) + (v).c21*conj((w).c01)  + (v).c22*conj((w).c02); \
  (u).c21 +=  (v).c20*conj((w).c10) + (v).c21*conj((w).c11)  + (v).c22*conj((w).c12); \
  (u).c22 +=  (v).c20*conj((w).c20) + (v).c21*conj((w).c21)  + (v).c22*conj((w).c22);

/* C.Urbach u = v^dag * w */
#define _su3d_times_su3(u,v,w)			\
  (u).c00 = conj((v).c00)*(w).c00 + conj((v).c10)*(w).c10 + conj((v).c20)*(w).c20; \
  (u).c01 = conj((v).c00)*(w).c01 + conj((v).c10)*(w).c11 + conj((v).c20)*(w).c21; \
  (u).c02 = conj((v).c00)*(w).c02 + conj((v).c10)*(w).c12 + conj((v).c20)*(w).c22; \
  (u).c10 = conj((v).c01)*(w).c00 + conj((v).c11)*(w).c10 + conj((v).c21)*(w).c20; \
  (u).c11 = conj((v).c01)*(w).c01 + conj((v).c11)*(w).c11 + conj((v).c21)*(w).c21; \
  (u).c12 = conj((v).c01)*(w).c02 + conj((v).c11)*(w).c12 + conj((v).c21)*(w).c22; \
  (u).c20 = conj((v).c02)*(w).c00 + conj((v).c12)*(w).c10 + conj((v).c22)*(w).c20; \
  (u).c21 = conj((v).c02)*(w).c01 + conj((v).c12)*(w).c11 + conj((v).c22)*(w).c21; \
  (u).c22 = conj((v).c02)*(w).c02 + conj((v).c12)*(w).c12 + conj((v).c22)*(w).c22;

/* C.Urbach u = u + v^dag * w */
#define _su3d_times_su3_acc(u,v,w)						\
  (u).c00 += conj((v).c00)*(w).c00 + conj((v).c10)*(w).c10 + conj((v).c20)*(w).c20; \
  (u).c01 += conj((v).c00)*(w).c01 + conj((v).c10)*(w).c11 + conj((v).c20)*(w).c21; \
  (u).c02 += conj((v).c00)*(w).c02 + conj((v).c10)*(w).c12 + conj((v).c20)*(w).c22; \
  (u).c10 += conj((v).c01)*(w).c00 + conj((v).c11)*(w).c10 + conj((v).c21)*(w).c20; \
  (u).c11 += conj((v).c01)*(w).c01 + conj((v).c11)*(w).c11 + conj((v).c21)*(w).c21; \
  (u).c12 += conj((v).c01)*(w).c02 + conj((v).c11)*(w).c12 + conj((v).c21)*(w).c22; \
  (u).c20 += conj((v).c02)*(w).c00 + conj((v).c12)*(w).c10 + conj((v).c22)*(w).c20; \
  (u).c21 += conj((v).c02)*(w).c01 + conj((v).c12)*(w).c11 + conj((v).c22)*(w).c21; \
  (u).c22 += conj((v).c02)*(w).c02 + conj((v).c12)*(w).c12 + conj((v).c22)*(w).c22;

#else

/*
* u=v*w
*/

#define _su3_times_su3(u,v,w) \
   (u).c00.re= (v).c00.re*(w).c00.re-(v).c00.im*(w).c00.im  \
              +(v).c01.re*(w).c10.re-(v).c01.im*(w).c10.im  \
              +(v).c02.re*(w).c20.re-(v).c02.im*(w).c20.im; \
   (u).c00.im= (v).c00.re*(w).c00.im+(v).c00.im*(w).c00.re  \
              +(v).c01.re*(w).c10.im+(v).c01.im*(w).c10.re  \
              +(v).c02.re*(w).c20.im+(v).c02.im*(w).c20.re; \
   (u).c01.re= (v).c00.re*(w).c01.re-(v).c00.im*(w).c01.im  \
              +(v).c01.re*(w).c11.re-(v).c01.im*(w).c11.im  \
              +(v).c02.re*(w).c21.re-(v).c02.im*(w).c21.im; \
   (u).c01.im= (v).c00.re*(w).c01.im+(v).c00.im*(w).c01.re  \
              +(v).c01.re*(w).c11.im+(v).c01.im*(w).c11.re  \
              +(v).c02.re*(w).c21.im+(v).c02.im*(w).c21.re; \
   (u).c02.re= (v).c00.re*(w).c02.re-(v).c00.im*(w).c02.im  \
              +(v).c01.re*(w).c12.re-(v).c01.im*(w).c12.im  \
              +(v).c02.re*(w).c22.re-(v).c02.im*(w).c22.im; \
   (u).c02.im= (v).c00.re*(w).c02.im+(v).c00.im*(w).c02.re  \
              +(v).c01.re*(w).c12.im+(v).c01.im*(w).c12.re  \
              +(v).c02.re*(w).c22.im+(v).c02.im*(w).c22.re; \
   (u).c10.re= (v).c10.re*(w).c00.re-(v).c10.im*(w).c00.im  \
              +(v).c11.re*(w).c10.re-(v).c11.im*(w).c10.im  \
              +(v).c12.re*(w).c20.re-(v).c12.im*(w).c20.im; \
   (u).c10.im= (v).c10.re*(w).c00.im+(v).c10.im*(w).c00.re  \
              +(v).c11.re*(w).c10.im+(v).c11.im*(w).c10.re  \
              +(v).c12.re*(w).c20.im+(v).c12.im*(w).c20.re; \
   (u).c11.re= (v).c10.re*(w).c01.re-(v).c10.im*(w).c01.im  \
              +(v).c11.re*(w).c11.re-(v).c11.im*(w).c11.im  \
              +(v).c12.re*(w).c21.re-(v).c12.im*(w).c21.im; \
   (u).c11.im= (v).c10.re*(w).c01.im+(v).c10.im*(w).c01.re  \
              +(v).c11.re*(w).c11.im+(v).c11.im*(w).c11.re  \
              +(v).c12.re*(w).c21.im+(v).c12.im*(w).c21.re; \
   (u).c12.re= (v).c10.re*(w).c02.re-(v).c10.im*(w).c02.im  \
              +(v).c11.re*(w).c12.re-(v).c11.im*(w).c12.im  \
              +(v).c12.re*(w).c22.re-(v).c12.im*(w).c22.im; \
   (u).c12.im= (v).c10.re*(w).c02.im+(v).c10.im*(w).c02.re  \
              +(v).c11.re*(w).c12.im+(v).c11.im*(w).c12.re  \
              +(v).c12.re*(w).c22.im+(v).c12.im*(w).c22.re; \
   (u).c20.re= (v).c20.re*(w).c00.re-(v).c20.im*(w).c00.im  \
              +(v).c21.re*(w).c10.re-(v).c21.im*(w).c10.im  \
              +(v).c22.re*(w).c20.re-(v).c22.im*(w).c20.im; \
   (u).c20.im= (v).c20.re*(w).c00.im+(v).c20.im*(w).c00.re  \
              +(v).c21.re*(w).c10.im+(v).c21.im*(w).c10.re  \
              +(v).c22.re*(w).c20.im+(v).c22.im*(w).c20.re; \
   (u).c21.re= (v).c20.re*(w).c01.re-(v).c20.im*(w).c01.im  \
              +(v).c21.re*(w).c11.re-(v).c21.im*(w).c11.im  \
              +(v).c22.re*(w).c21.re-(v).c22.im*(w).c21.im; \
   (u).c21.im= (v).c20.re*(w).c01.im+(v).c20.im*(w).c01.re  \
              +(v).c21.re*(w).c11.im+(v).c21.im*(w).c11.re  \
              +(v).c22.re*(w).c21.im+(v).c22.im*(w).c21.re; \
   (u).c22.re= (v).c20.re*(w).c02.re-(v).c20.im*(w).c02.im  \
              +(v).c21.re*(w).c12.re-(v).c21.im*(w).c12.im  \
              +(v).c22.re*(w).c22.re-(v).c22.im*(w).c22.im; \
   (u).c22.im= (v).c20.re*(w).c02.im+(v).c20.im*(w).c02.re  \
              +(v).c21.re*(w).c12.im+(v).c21.im*(w).c12.re  \
              +(v).c22.re*(w).c22.im+(v).c22.im*(w).c22.re;

/* M.Hasenbusch u = u + v * w */

#define _su3_times_su3_acc(u,v,w) \
   (u).c00.re+=(v).c00.re*(w).c00.re-(v).c00.im*(w).c00.im  \
              +(v).c01.re*(w).c10.re-(v).c01.im*(w).c10.im  \
              +(v).c02.re*(w).c20.re-(v).c02.im*(w).c20.im; \
   (u).c00.im+=(v).c00.re*(w).c00.im+(v).c00.im*(w).c00.re  \
              +(v).c01.re*(w).c10.im+(v).c01.im*(w).c10.re  \
              +(v).c02.re*(w).c20.im+(v).c02.im*(w).c20.re; \
   (u).c01.re+=(v).c00.re*(w).c01.re-(v).c00.im*(w).c01.im  \
              +(v).c01.re*(w).c11.re-(v).c01.im*(w).c11.im  \
              +(v).c02.re*(w).c21.re-(v).c02.im*(w).c21.im; \
   (u).c01.im+=(v).c00.re*(w).c01.im+(v).c00.im*(w).c01.re  \
              +(v).c01.re*(w).c11.im+(v).c01.im*(w).c11.re  \
              +(v).c02.re*(w).c21.im+(v).c02.im*(w).c21.re; \
   (u).c02.re+=(v).c00.re*(w).c02.re-(v).c00.im*(w).c02.im  \
              +(v).c01.re*(w).c12.re-(v).c01.im*(w).c12.im  \
              +(v).c02.re*(w).c22.re-(v).c02.im*(w).c22.im; \
   (u).c02.im+=(v).c00.re*(w).c02.im+(v).c00.im*(w).c02.re  \
              +(v).c01.re*(w).c12.im+(v).c01.im*(w).c12.re  \
              +(v).c02.re*(w).c22.im+(v).c02.im*(w).c22.re; \
   (u).c10.re+=(v).c10.re*(w).c00.re-(v).c10.im*(w).c00.im  \
              +(v).c11.re*(w).c10.re-(v).c11.im*(w).c10.im  \
              +(v).c12.re*(w).c20.re-(v).c12.im*(w).c20.im; \
   (u).c10.im+=(v).c10.re*(w).c00.im+(v).c10.im*(w).c00.re  \
              +(v).c11.re*(w).c10.im+(v).c11.im*(w).c10.re  \
              +(v).c12.re*(w).c20.im+(v).c12.im*(w).c20.re; \
   (u).c11.re+=(v).c10.re*(w).c01.re-(v).c10.im*(w).c01.im  \
              +(v).c11.re*(w).c11.re-(v).c11.im*(w).c11.im  \
              +(v).c12.re*(w).c21.re-(v).c12.im*(w).c21.im; \
   (u).c11.im+=(v).c10.re*(w).c01.im+(v).c10.im*(w).c01.re  \
              +(v).c11.re*(w).c11.im+(v).c11.im*(w).c11.re  \
              +(v).c12.re*(w).c21.im+(v).c12.im*(w).c21.re; \
   (u).c12.re+=(v).c10.re*(w).c02.re-(v).c10.im*(w).c02.im  \
              +(v).c11.re*(w).c12.re-(v).c11.im*(w).c12.im  \
              +(v).c12.re*(w).c22.re-(v).c12.im*(w).c22.im; \
   (u).c12.im+=(v).c10.re*(w).c02.im+(v).c10.im*(w).c02.re  \
              +(v).c11.re*(w).c12.im+(v).c11.im*(w).c12.re  \
              +(v).c12.re*(w).c22.im+(v).c12.im*(w).c22.re; \
   (u).c20.re+=(v).c20.re*(w).c00.re-(v).c20.im*(w).c00.im  \
              +(v).c21.re*(w).c10.re-(v).c21.im*(w).c10.im  \
              +(v).c22.re*(w).c20.re-(v).c22.im*(w).c20.im; \
   (u).c20.im+=(v).c20.re*(w).c00.im+(v).c20.im*(w).c00.re  \
              +(v).c21.re*(w).c10.im+(v).c21.im*(w).c10.re  \
              +(v).c22.re*(w).c20.im+(v).c22.im*(w).c20.re; \
   (u).c21.re+=(v).c20.re*(w).c01.re-(v).c20.im*(w).c01.im  \
              +(v).c21.re*(w).c11.re-(v).c21.im*(w).c11.im  \
              +(v).c22.re*(w).c21.re-(v).c22.im*(w).c21.im; \
   (u).c21.im+=(v).c20.re*(w).c01.im+(v).c20.im*(w).c01.re  \
              +(v).c21.re*(w).c11.im+(v).c21.im*(w).c11.re  \
              +(v).c22.re*(w).c21.im+(v).c22.im*(w).c21.re; \
   (u).c22.re+=(v).c20.re*(w).c02.re-(v).c20.im*(w).c02.im  \
              +(v).c21.re*(w).c12.re-(v).c21.im*(w).c12.im  \
              +(v).c22.re*(w).c22.re-(v).c22.im*(w).c22.im; \
   (u).c22.im+=(v).c20.re*(w).c02.im+(v).c20.im*(w).c02.re  \
              +(v).c21.re*(w).c12.im+(v).c21.im*(w).c12.re  \
              +(v).c22.re*(w).c22.im+(v).c22.im*(w).c22.re;

/* M.Hasenbusch u=v * w^{dag} */

#define _su3_times_su3d(u,v,w) \
   (u).c00.re= (v).c00.re*(w).c00.re+(v).c00.im*(w).c00.im  \
              +(v).c01.re*(w).c01.re+(v).c01.im*(w).c01.im  \
              +(v).c02.re*(w).c02.re+(v).c02.im*(w).c02.im; \
   (u).c00.im=-(v).c00.re*(w).c00.im+(v).c00.im*(w).c00.re  \
              -(v).c01.re*(w).c01.im+(v).c01.im*(w).c01.re  \
              -(v).c02.re*(w).c02.im+(v).c02.im*(w).c02.re; \
   (u).c01.re= (v).c00.re*(w).c10.re+(v).c00.im*(w).c10.im  \
              +(v).c01.re*(w).c11.re+(v).c01.im*(w).c11.im  \
              +(v).c02.re*(w).c12.re+(v).c02.im*(w).c12.im; \
   (u).c01.im=-(v).c00.re*(w).c10.im+(v).c00.im*(w).c10.re  \
              -(v).c01.re*(w).c11.im+(v).c01.im*(w).c11.re  \
              -(v).c02.re*(w).c12.im+(v).c02.im*(w).c12.re; \
   (u).c02.re= (v).c00.re*(w).c20.re+(v).c00.im*(w).c20.im  \
              +(v).c01.re*(w).c21.re+(v).c01.im*(w).c21.im  \
              +(v).c02.re*(w).c22.re+(v).c02.im*(w).c22.im; \
   (u).c02.im=-(v).c00.re*(w).c20.im+(v).c00.im*(w).c20.re  \
              -(v).c01.re*(w).c21.im+(v).c01.im*(w).c21.re  \
              -(v).c02.re*(w).c22.im+(v).c02.im*(w).c22.re; \
   (u).c10.re= (v).c10.re*(w).c00.re+(v).c10.im*(w).c00.im  \
              +(v).c11.re*(w).c01.re+(v).c11.im*(w).c01.im  \
              +(v).c12.re*(w).c02.re+(v).c12.im*(w).c02.im; \
   (u).c10.im=-(v).c10.re*(w).c00.im+(v).c10.im*(w).c00.re  \
              -(v).c11.re*(w).c01.im+(v).c11.im*(w).c01.re  \
              -(v).c12.re*(w).c02.im+(v).c12.im*(w).c02.re; \
   (u).c11.re= (v).c10.re*(w).c10.re+(v).c10.im*(w).c10.im  \
              +(v).c11.re*(w).c11.re+(v).c11.im*(w).c11.im  \
              +(v).c12.re*(w).c12.re+(v).c12.im*(w).c12.im; \
   (u).c11.im=-(v).c10.re*(w).c10.im+(v).c10.im*(w).c10.re  \
              -(v).c11.re*(w).c11.im+(v).c11.im*(w).c11.re  \
              -(v).c12.re*(w).c12.im+(v).c12.im*(w).c12.re; \
   (u).c12.re= (v).c10.re*(w).c20.re+(v).c10.im*(w).c20.im  \
              +(v).c11.re*(w).c21.re+(v).c11.im*(w).c21.im  \
              +(v).c12.re*(w).c22.re+(v).c12.im*(w).c22.im; \
   (u).c12.im=-(v).c10.re*(w).c20.im+(v).c10.im*(w).c20.re  \
              -(v).c11.re*(w).c21.im+(v).c11.im*(w).c21.re  \
              -(v).c12.re*(w).c22.im+(v).c12.im*(w).c22.re; \
   (u).c20.re= (v).c20.re*(w).c00.re+(v).c20.im*(w).c00.im  \
              +(v).c21.re*(w).c01.re+(v).c21.im*(w).c01.im  \
              +(v).c22.re*(w).c02.re+(v).c22.im*(w).c02.im; \
   (u).c20.im=-(v).c20.re*(w).c00.im+(v).c20.im*(w).c00.re  \
              -(v).c21.re*(w).c01.im+(v).c21.im*(w).c01.re  \
              -(v).c22.re*(w).c02.im+(v).c22.im*(w).c02.re; \
   (u).c21.re= (v).c20.re*(w).c10.re+(v).c20.im*(w).c10.im  \
              +(v).c21.re*(w).c11.re+(v).c21.im*(w).c11.im  \
              +(v).c22.re*(w).c12.re+(v).c22.im*(w).c12.im; \
   (u).c21.im=-(v).c20.re*(w).c10.im+(v).c20.im*(w).c10.re  \
              -(v).c21.re*(w).c11.im+(v).c21.im*(w).c11.re  \
              -(v).c22.re*(w).c12.im+(v).c22.im*(w).c12.re; \
   (u).c22.re= (v).c20.re*(w).c20.re+(v).c20.im*(w).c20.im  \
              +(v).c21.re*(w).c21.re+(v).c21.im*(w).c21.im  \
              +(v).c22.re*(w).c22.re+(v).c22.im*(w).c22.im; \
   (u).c22.im=-(v).c20.re*(w).c20.im+(v).c20.im*(w).c20.re  \
              -(v).c21.re*(w).c21.im+(v).c21.im*(w).c21.re  \
              -(v).c22.re*(w).c22.im+(v).c22.im*(w).c22.re;

/* C.Urbach u=u + v * w^{dag} */

#define _su3_times_su3d_acc(u,v,w) \
   (u).c00.re+=(v).c00.re*(w).c00.re+(v).c00.im*(w).c00.im  \
              +(v).c01.re*(w).c01.re+(v).c01.im*(w).c01.im  \
              +(v).c02.re*(w).c02.re+(v).c02.im*(w).c02.im; \
   (u).c00.im+=-(v).c00.re*(w).c00.im+(v).c00.im*(w).c00.re  \
              -(v).c01.re*(w).c01.im+(v).c01.im*(w).c01.re  \
              -(v).c02.re*(w).c02.im+(v).c02.im*(w).c02.re; \
   (u).c01.re+=(v).c00.re*(w).c10.re+(v).c00.im*(w).c10.im  \
              +(v).c01.re*(w).c11.re+(v).c01.im*(w).c11.im  \
              +(v).c02.re*(w).c12.re+(v).c02.im*(w).c12.im; \
   (u).c01.im+=-(v).c00.re*(w).c10.im+(v).c00.im*(w).c10.re  \
              -(v).c01.re*(w).c11.im+(v).c01.im*(w).c11.re  \
              -(v).c02.re*(w).c12.im+(v).c02.im*(w).c12.re; \
   (u).c02.re+=(v).c00.re*(w).c20.re+(v).c00.im*(w).c20.im  \
              +(v).c01.re*(w).c21.re+(v).c01.im*(w).c21.im  \
              +(v).c02.re*(w).c22.re+(v).c02.im*(w).c22.im; \
   (u).c02.im+=-(v).c00.re*(w).c20.im+(v).c00.im*(w).c20.re  \
              -(v).c01.re*(w).c21.im+(v).c01.im*(w).c21.re  \
              -(v).c02.re*(w).c22.im+(v).c02.im*(w).c22.re; \
   (u).c10.re+=(v).c10.re*(w).c00.re+(v).c10.im*(w).c00.im  \
              +(v).c11.re*(w).c01.re+(v).c11.im*(w).c01.im  \
              +(v).c12.re*(w).c02.re+(v).c12.im*(w).c02.im; \
   (u).c10.im+=-(v).c10.re*(w).c00.im+(v).c10.im*(w).c00.re  \
              -(v).c11.re*(w).c01.im+(v).c11.im*(w).c01.re  \
              -(v).c12.re*(w).c02.im+(v).c12.im*(w).c02.re; \
   (u).c11.re+=(v).c10.re*(w).c10.re+(v).c10.im*(w).c10.im  \
              +(v).c11.re*(w).c11.re+(v).c11.im*(w).c11.im  \
              +(v).c12.re*(w).c12.re+(v).c12.im*(w).c12.im; \
   (u).c11.im+=-(v).c10.re*(w).c10.im+(v).c10.im*(w).c10.re  \
              -(v).c11.re*(w).c11.im+(v).c11.im*(w).c11.re  \
              -(v).c12.re*(w).c12.im+(v).c12.im*(w).c12.re; \
   (u).c12.re+=(v).c10.re*(w).c20.re+(v).c10.im*(w).c20.im  \
              +(v).c11.re*(w).c21.re+(v).c11.im*(w).c21.im  \
              +(v).c12.re*(w).c22.re+(v).c12.im*(w).c22.im; \
   (u).c12.im+=-(v).c10.re*(w).c20.im+(v).c10.im*(w).c20.re  \
              -(v).c11.re*(w).c21.im+(v).c11.im*(w).c21.re  \
              -(v).c12.re*(w).c22.im+(v).c12.im*(w).c22.re; \
   (u).c20.re+=(v).c20.re*(w).c00.re+(v).c20.im*(w).c00.im  \
              +(v).c21.re*(w).c01.re+(v).c21.im*(w).c01.im  \
              +(v).c22.re*(w).c02.re+(v).c22.im*(w).c02.im; \
   (u).c20.im+=-(v).c20.re*(w).c00.im+(v).c20.im*(w).c00.re  \
              -(v).c21.re*(w).c01.im+(v).c21.im*(w).c01.re  \
              -(v).c22.re*(w).c02.im+(v).c22.im*(w).c02.re; \
   (u).c21.re+=(v).c20.re*(w).c10.re+(v).c20.im*(w).c10.im  \
              +(v).c21.re*(w).c11.re+(v).c21.im*(w).c11.im  \
              +(v).c22.re*(w).c12.re+(v).c22.im*(w).c12.im; \
   (u).c21.im+=-(v).c20.re*(w).c10.im+(v).c20.im*(w).c10.re  \
              -(v).c21.re*(w).c11.im+(v).c21.im*(w).c11.re  \
              -(v).c22.re*(w).c12.im+(v).c22.im*(w).c12.re; \
   (u).c22.re+=(v).c20.re*(w).c20.re+(v).c20.im*(w).c20.im  \
              +(v).c21.re*(w).c21.re+(v).c21.im*(w).c21.im  \
              +(v).c22.re*(w).c22.re+(v).c22.im*(w).c22.im; \
   (u).c22.im+=-(v).c20.re*(w).c20.im+(v).c20.im*(w).c20.re  \
              -(v).c21.re*(w).c21.im+(v).c21.im*(w).c21.re  \
              -(v).c22.re*(w).c22.im+(v).c22.im*(w).c22.re;

/* M.Hasenbusch u=v^{dag} w */

#define _su3d_times_su3(u,v,w) \
   (u).c00.re= (v).c00.re*(w).c00.re+(v).c00.im*(w).c00.im  \
              +(v).c10.re*(w).c10.re+(v).c10.im*(w).c10.im  \
              +(v).c20.re*(w).c20.re+(v).c20.im*(w).c20.im; \
   (u).c00.im= (v).c00.re*(w).c00.im-(v).c00.im*(w).c00.re  \
              +(v).c10.re*(w).c10.im-(v).c10.im*(w).c10.re  \
              +(v).c20.re*(w).c20.im-(v).c20.im*(w).c20.re; \
   (u).c01.re= (v).c00.re*(w).c01.re+(v).c00.im*(w).c01.im  \
              +(v).c10.re*(w).c11.re+(v).c10.im*(w).c11.im  \
              +(v).c20.re*(w).c21.re+(v).c20.im*(w).c21.im; \
   (u).c01.im= (v).c00.re*(w).c01.im-(v).c00.im*(w).c01.re  \
              +(v).c10.re*(w).c11.im-(v).c10.im*(w).c11.re  \
              +(v).c20.re*(w).c21.im-(v).c20.im*(w).c21.re; \
   (u).c02.re= (v).c00.re*(w).c02.re+(v).c00.im*(w).c02.im  \
              +(v).c10.re*(w).c12.re+(v).c10.im*(w).c12.im  \
              +(v).c20.re*(w).c22.re+(v).c20.im*(w).c22.im; \
   (u).c02.im= (v).c00.re*(w).c02.im-(v).c00.im*(w).c02.re  \
              +(v).c10.re*(w).c12.im-(v).c10.im*(w).c12.re  \
              +(v).c20.re*(w).c22.im-(v).c20.im*(w).c22.re; \
   (u).c10.re= (v).c01.re*(w).c00.re+(v).c01.im*(w).c00.im  \
              +(v).c11.re*(w).c10.re+(v).c11.im*(w).c10.im  \
              +(v).c21.re*(w).c20.re+(v).c21.im*(w).c20.im; \
   (u).c10.im= (v).c01.re*(w).c00.im-(v).c01.im*(w).c00.re  \
              +(v).c11.re*(w).c10.im-(v).c11.im*(w).c10.re  \
              +(v).c21.re*(w).c20.im-(v).c21.im*(w).c20.re; \
   (u).c11.re= (v).c01.re*(w).c01.re+(v).c01.im*(w).c01.im  \
              +(v).c11.re*(w).c11.re+(v).c11.im*(w).c11.im  \
              +(v).c21.re*(w).c21.re+(v).c21.im*(w).c21.im; \
   (u).c11.im= (v).c01.re*(w).c01.im-(v).c01.im*(w).c01.re  \
              +(v).c11.re*(w).c11.im-(v).c11.im*(w).c11.re  \
              +(v).c21.re*(w).c21.im-(v).c21.im*(w).c21.re; \
   (u).c12.re= (v).c01.re*(w).c02.re+(v).c01.im*(w).c02.im  \
              +(v).c11.re*(w).c12.re+(v).c11.im*(w).c12.im  \
              +(v).c21.re*(w).c22.re+(v).c21.im*(w).c22.im; \
   (u).c12.im= (v).c01.re*(w).c02.im-(v).c01.im*(w).c02.re  \
              +(v).c11.re*(w).c12.im-(v).c11.im*(w).c12.re  \
              +(v).c21.re*(w).c22.im-(v).c21.im*(w).c22.re; \
   (u).c20.re= (v).c02.re*(w).c00.re+(v).c02.im*(w).c00.im  \
              +(v).c12.re*(w).c10.re+(v).c12.im*(w).c10.im  \
              +(v).c22.re*(w).c20.re+(v).c22.im*(w).c20.im; \
   (u).c20.im= (v).c02.re*(w).c00.im-(v).c02.im*(w).c00.re  \
              +(v).c12.re*(w).c10.im-(v).c12.im*(w).c10.re  \
              +(v).c22.re*(w).c20.im-(v).c22.im*(w).c20.re; \
   (u).c21.re= (v).c02.re*(w).c01.re+(v).c02.im*(w).c01.im  \
              +(v).c12.re*(w).c11.re+(v).c12.im*(w).c11.im  \
              +(v).c22.re*(w).c21.re+(v).c22.im*(w).c21.im; \
   (u).c21.im= (v).c02.re*(w).c01.im-(v).c02.im*(w).c01.re  \
              +(v).c12.re*(w).c11.im-(v).c12.im*(w).c11.re  \
              +(v).c22.re*(w).c21.im-(v).c22.im*(w).c21.re; \
   (u).c22.re= (v).c02.re*(w).c02.re+(v).c02.im*(w).c02.im  \
              +(v).c12.re*(w).c12.re+(v).c12.im*(w).c12.im  \
              +(v).c22.re*(w).c22.re+(v).c22.im*(w).c22.im; \
   (u).c22.im= (v).c02.re*(w).c02.im-(v).c02.im*(w).c02.re  \
              +(v).c12.re*(w).c12.im-(v).c12.im*(w).c12.re  \
              +(v).c22.re*(w).c22.im-(v).c22.im*(w).c22.re;

#define _su3d_times_su3_acc(u,v,w) \
   (u).c00.re+=(v).c00.re*(w).c00.re+(v).c00.im*(w).c00.im  \
              +(v).c10.re*(w).c10.re+(v).c10.im*(w).c10.im  \
              +(v).c20.re*(w).c20.re+(v).c20.im*(w).c20.im; \
   (u).c00.im+=(v).c00.re*(w).c00.im-(v).c00.im*(w).c00.re  \
              +(v).c10.re*(w).c10.im-(v).c10.im*(w).c10.re  \
              +(v).c20.re*(w).c20.im-(v).c20.im*(w).c20.re; \
   (u).c01.re+=(v).c00.re*(w).c01.re+(v).c00.im*(w).c01.im  \
              +(v).c10.re*(w).c11.re+(v).c10.im*(w).c11.im  \
              +(v).c20.re*(w).c21.re+(v).c20.im*(w).c21.im; \
   (u).c01.im+=(v).c00.re*(w).c01.im-(v).c00.im*(w).c01.re  \
              +(v).c10.re*(w).c11.im-(v).c10.im*(w).c11.re  \
              +(v).c20.re*(w).c21.im-(v).c20.im*(w).c21.re; \
   (u).c02.re+=(v).c00.re*(w).c02.re+(v).c00.im*(w).c02.im  \
              +(v).c10.re*(w).c12.re+(v).c10.im*(w).c12.im  \
              +(v).c20.re*(w).c22.re+(v).c20.im*(w).c22.im; \
   (u).c02.im+=(v).c00.re*(w).c02.im-(v).c00.im*(w).c02.re  \
              +(v).c10.re*(w).c12.im-(v).c10.im*(w).c12.re  \
              +(v).c20.re*(w).c22.im-(v).c20.im*(w).c22.re; \
   (u).c10.re+=(v).c01.re*(w).c00.re+(v).c01.im*(w).c00.im  \
              +(v).c11.re*(w).c10.re+(v).c11.im*(w).c10.im  \
              +(v).c21.re*(w).c20.re+(v).c21.im*(w).c20.im; \
   (u).c10.im+=(v).c01.re*(w).c00.im-(v).c01.im*(w).c00.re  \
              +(v).c11.re*(w).c10.im-(v).c11.im*(w).c10.re  \
              +(v).c21.re*(w).c20.im-(v).c21.im*(w).c20.re; \
   (u).c11.re+=(v).c01.re*(w).c01.re+(v).c01.im*(w).c01.im  \
              +(v).c11.re*(w).c11.re+(v).c11.im*(w).c11.im  \
              +(v).c21.re*(w).c21.re+(v).c21.im*(w).c21.im; \
   (u).c11.im+=(v).c01.re*(w).c01.im-(v).c01.im*(w).c01.re  \
              +(v).c11.re*(w).c11.im-(v).c11.im*(w).c11.re  \
              +(v).c21.re*(w).c21.im-(v).c21.im*(w).c21.re; \
   (u).c12.re+=(v).c01.re*(w).c02.re+(v).c01.im*(w).c02.im  \
              +(v).c11.re*(w).c12.re+(v).c11.im*(w).c12.im  \
              +(v).c21.re*(w).c22.re+(v).c21.im*(w).c22.im; \
   (u).c12.im+=(v).c01.re*(w).c02.im-(v).c01.im*(w).c02.re  \
              +(v).c11.re*(w).c12.im-(v).c11.im*(w).c12.re  \
              +(v).c21.re*(w).c22.im-(v).c21.im*(w).c22.re; \
   (u).c20.re+=(v).c02.re*(w).c00.re+(v).c02.im*(w).c00.im  \
              +(v).c12.re*(w).c10.re+(v).c12.im*(w).c10.im  \
              +(v).c22.re*(w).c20.re+(v).c22.im*(w).c20.im; \
   (u).c20.im+=(v).c02.re*(w).c00.im-(v).c02.im*(w).c00.re  \
              +(v).c12.re*(w).c10.im-(v).c12.im*(w).c10.re  \
              +(v).c22.re*(w).c20.im-(v).c22.im*(w).c20.re; \
   (u).c21.re+=(v).c02.re*(w).c01.re+(v).c02.im*(w).c01.im  \
              +(v).c12.re*(w).c11.re+(v).c12.im*(w).c11.im  \
              +(v).c22.re*(w).c21.re+(v).c22.im*(w).c21.im; \
   (u).c21.im+=(v).c02.re*(w).c01.im-(v).c02.im*(w).c01.re  \
              +(v).c12.re*(w).c11.im-(v).c12.im*(w).c11.re  \
              +(v).c22.re*(w).c21.im-(v).c22.im*(w).c21.re; \
   (u).c22.re+=(v).c02.re*(w).c02.re+(v).c02.im*(w).c02.im  \
              +(v).c12.re*(w).c12.re+(v).c12.im*(w).c12.im  \
              +(v).c22.re*(w).c22.re+(v).c22.im*(w).c22.im; \
   (u).c22.im+=(v).c02.re*(w).c02.im-(v).c02.im*(w).c02.re  \
              +(v).c12.re*(w).c12.im-(v).c12.im*(w).c12.re  \
              +(v).c22.re*(w).c22.im-(v).c22.im*(w).c22.re;

#endif

/* M. Hasenbusch x=Re Tr (v * w^{\dag}) */
#ifdef _STD_C99_COMPLEX
#define _trace_su3_times_su3d(x,v,w)		    \
  x = (v).c00*conj((w).c00)			    \
    +(v).c01*conj((w).c01)			    \
    +(v).c02*conj((w).c02)			    \
    +(v).c10*conj((w).c10)			    \
    +(v).c11*conj((w).c11)			    \
    +(v).c12*conj((w).c12)			    \
    +(v).c20*conj((w).c20)			    \
    +(v).c21*conj((w).c21)			    \
    +(v).c22*conj((w).c22);

#else

#define _trace_su3_times_su3d(x,v,w) \
   x = (v).c00.re*(w).c00.re+(v).c00.im*(w).c00.im  \
      +(v).c01.re*(w).c01.re+(v).c01.im*(w).c01.im  \
      +(v).c02.re*(w).c02.re+(v).c02.im*(w).c02.im  \
      +(v).c10.re*(w).c10.re+(v).c10.im*(w).c10.im  \
      +(v).c11.re*(w).c11.re+(v).c11.im*(w).c11.im  \
      +(v).c12.re*(w).c12.re+(v).c12.im*(w).c12.im  \
      +(v).c20.re*(w).c20.re+(v).c20.im*(w).c20.im  \
      +(v).c21.re*(w).c21.re+(v).c21.im*(w).c21.im  \
      +(v).c22.re*(w).c22.re+(v).c22.im*(w).c22.im;
#endif

/* U. Wenger x=Tr (v * w) */

#ifdef _STD_C99_COMPLEX
#define _trace_su3_times_su3(x,v,w)		\
  x = (v).c00*(w).c00				\
    + (v).c01*(w).c10				\
    + (v).c02*(w).c20				\
    + (v).c10*(w).c01				\
    + (v).c11*(w).c11				\
    + (v).c12*(w).c21				\
    + (v).c20*(w).c02				\
    + (v).c21*(w).c12				\
    + (v).c22*(w).c22;
#else
#define _trace_su3_times_su3(x,v,w) \
   x.re = (v).c00.re*(w).c00.re-(v).c00.im*(w).c00.im  \
         +(v).c01.re*(w).c10.re-(v).c01.im*(w).c10.im  \
         +(v).c02.re*(w).c20.re-(v).c02.im*(w).c20.im  \
         +(v).c10.re*(w).c01.re-(v).c10.im*(w).c01.im  \
         +(v).c11.re*(w).c11.re-(v).c11.im*(w).c11.im  \
         +(v).c12.re*(w).c21.re-(v).c12.im*(w).c21.im  \
         +(v).c20.re*(w).c02.re-(v).c20.im*(w).c02.im  \
         +(v).c21.re*(w).c12.re-(v).c21.im*(w).c12.im  \
         +(v).c22.re*(w).c22.re-(v).c22.im*(w).c22.im; \
   x.im = (v).c00.re*(w).c00.im+(v).c00.im*(w).c00.re  \
         +(v).c01.re*(w).c10.im+(v).c01.im*(w).c10.re  \
         +(v).c02.re*(w).c20.im+(v).c02.im*(w).c20.re  \
         +(v).c10.re*(w).c01.im+(v).c10.im*(w).c01.re  \
         +(v).c11.re*(w).c11.im+(v).c11.im*(w).c11.re  \
         +(v).c12.re*(w).c21.im+(v).c12.im*(w).c21.re  \
         +(v).c20.re*(w).c02.im+(v).c20.im*(w).c02.re  \
         +(v).c21.re*(w).c12.im+(v).c21.im*(w).c12.re  \
         +(v).c22.re*(w).c22.im+(v).c22.im*(w).c22.re; 
#endif

/* M. Hasenbusch t =u tensor v^dag */
#ifdef _STD_C99_COMPLEX
#define _vector_tensor_vector(t,u,v)	\
  (t).c00 = (u).c0*conj((v).c0);	\
  (t).c01 = (u).c0*conj((v).c1);	\
  (t).c02 = (u).c0*conj((v).c2);	\
  (t).c10 = (u).c1*conj((v).c0);	\
  (t).c11 = (u).c1*conj((v).c1);	\
  (t).c12 = (u).c1*conj((v).c2);	\
  (t).c20 = (u).c2*conj((v).c0);	\
  (t).c21 = (u).c2*conj((v).c1);	\
  (t).c22 = (u).c2*conj((v).c2);
#else
#define _vector_tensor_vector(t,u,v) \
   (t).c00.re=(u).c0.re*(v).c0.re+(u).c0.im*(v).c0.im; \
   (t).c00.im=(u).c0.im*(v).c0.re-(u).c0.re*(v).c0.im; \
   (t).c01.re=(u).c0.re*(v).c1.re+(u).c0.im*(v).c1.im; \
   (t).c01.im=(u).c0.im*(v).c1.re-(u).c0.re*(v).c1.im; \
   (t).c02.re=(u).c0.re*(v).c2.re+(u).c0.im*(v).c2.im; \
   (t).c02.im=(u).c0.im*(v).c2.re-(u).c0.re*(v).c2.im; \
   (t).c10.re=(u).c1.re*(v).c0.re+(u).c1.im*(v).c0.im; \
   (t).c10.im=(u).c1.im*(v).c0.re-(u).c1.re*(v).c0.im; \
   (t).c11.re=(u).c1.re*(v).c1.re+(u).c1.im*(v).c1.im; \
   (t).c11.im=(u).c1.im*(v).c1.re-(u).c1.re*(v).c1.im; \
   (t).c12.re=(u).c1.re*(v).c2.re+(u).c1.im*(v).c2.im; \
   (t).c12.im=(u).c1.im*(v).c2.re-(u).c1.re*(v).c2.im; \
   (t).c20.re=(u).c2.re*(v).c0.re+(u).c2.im*(v).c0.im; \
   (t).c20.im=(u).c2.im*(v).c0.re-(u).c2.re*(v).c0.im; \
   (t).c21.re=(u).c2.re*(v).c1.re+(u).c2.im*(v).c1.im; \
   (t).c21.im=(u).c2.im*(v).c1.re-(u).c2.re*(v).c1.im; \
   (t).c22.re=(u).c2.re*(v).c2.re+(u).c2.im*(v).c2.im; \
   (t).c22.im=(u).c2.im*(v).c2.re-(u).c2.re*(v).c2.im; 
#endif

#ifdef MAIN_PROGRAM
static char const su3rcsid[] = "$Id$";
#endif

#endif

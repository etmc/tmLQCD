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

#include "complex.h"

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
   (r).c0.re=0.0; \
   (r).c0.im=0.0; \
   (r).c1.re=0.0; \
   (r).c1.im=0.0; \
   (r).c2.re=0.0; \
   (r).c2.im=0.0;

/* M. Hasenbusch Mon Sep 24
* r.c1=s.c1
* r.c2=s.c2
* r.c3=s.c3
*/

#define _vector_assign(r,s) \
   (r).c0.re=(s).c0.re; \
   (r).c0.im=(s).c0.im; \
   (r).c1.re=(s).c1.re; \
   (r).c1.im=(s).c1.im; \
   (r).c2.re=(s).c2.re; \
   (r).c2.im=(s).c2.im;

/* M. Hasenbusch Mon Sep 24
* r.c1=-s.c1
* r.c2=-s.c2
* r.c3=-s.c3
*/

#define _vector_minus_assign(r,s) \
   (r).c0.re=-(s).c0.re; \
   (r).c0.im=-(s).c0.im; \
   (r).c1.re=-(s).c1.re; \
   (r).c1.im=-(s).c1.im; \
   (r).c2.re=-(s).c2.re; \
   (r).c2.im=-(s).c2.im;

/*
* r.c1=c*s.c1 (c real)
* r.c2=c*s.c2
* r.c3=c*s.c3
*/

#define _vector_mul(r,c,s) \
   (r).c0.re=(c)*(s).c0.re; \
   (r).c0.im=(c)*(s).c0.im; \
   (r).c1.re=(c)*(s).c1.re; \
   (r).c1.im=(c)*(s).c1.im; \
   (r).c2.re=(c)*(s).c2.re; \
   (r).c2.im=(c)*(s).c2.im;

#define _vector_add_mul(r,c,s) \
   (r).c0.re+=(c)*(s).c0.re; \
   (r).c0.im+=(c)*(s).c0.im; \
   (r).c1.re+=(c)*(s).c1.re; \
   (r).c1.im+=(c)*(s).c1.im; \
   (r).c2.re+=(c)*(s).c2.re; \
   (r).c2.im+=(c)*(s).c2.im;

/*
* r.c1=s1.c1+s2.c1
* r.c2=s1.c2+s2.c2
* r.c3=s1.c3+s2.c3
*/

#if defined SSE2

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

#define _vector_i_add(r,s1,s2) \
   (r).c0.re=(s1).c0.re-(s2).c0.im; \
   (r).c0.im=(s1).c0.im+(s2).c0.re; \
   (r).c1.re=(s1).c1.re-(s2).c1.im; \
   (r).c1.im=(s1).c1.im+(s2).c1.re; \
   (r).c2.re=(s1).c2.re-(s2).c2.im; \
   (r).c2.im=(s1).c2.im+(s2).c2.re;

/*
* r.c1=s1.c1+i*s2.c1
* r.c2=s1.c2+i*s2.c2
* r.c3=s1.c3+i*s2.c3
*/

#define _vector_i_sub(r,s1,s2) \
   (r).c0.re=(s1).c0.re+(s2).c0.im; \
   (r).c0.im=(s1).c0.im-(s2).c0.re; \
   (r).c1.re=(s1).c1.re+(s2).c1.im; \
   (r).c1.im=(s1).c1.im-(s2).c1.re; \
   (r).c2.re=(s1).c2.re+(s2).c2.im; \
   (r).c2.im=(s1).c2.im-(s2).c2.re;

/*
* r.c1+=s.c1
* r.c2+=s.c2
* r.c3+=s.c3
*/

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



#if defined SSE2
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

#define _vector_i_add_assign(r,s) \
   (r).c0.re-=(s).c0.im; \
   (r).c0.im+=(s).c0.re; \
   (r).c1.re-=(s).c1.im; \
   (r).c1.im+=(s).c1.re; \
   (r).c2.re-=(s).c2.im; \
   (r).c2.im+=(s).c2.re;

/*
* r.c1-=i*s.c1
* r.c2-=i*s.c2
* r.c3-=i*s.c3
*/

#define _vector_i_sub_assign(r,s) \
   (r).c0.re+=(s).c0.im; \
   (r).c0.im-=(s).c0.re; \
   (r).c1.re+=(s).c1.im; \
   (r).c1.im-=(s).c1.re; \
   (r).c2.re+=(s).c2.im; \
   (r).c2.im-=(s).c2.re;

/* M.Hasenbusch 
* r.c1=c*s.c1
* r.c2=c*s.c2
* r.c3=c*s.c3
*
* c complex
*/

#define _complex_times_vector(r,c,s) \
   (r).c0.re=(c).re*(s).c0.re-(c).im*(s).c0.im; \
   (r).c0.im=(c).re*(s).c0.im+(c).im*(s).c0.re; \
   (r).c1.re=(c).re*(s).c1.re-(c).im*(s).c1.im; \
   (r).c1.im=(c).re*(s).c1.im+(c).im*(s).c1.re; \
   (r).c2.re=(c).re*(s).c2.re-(c).im*(s).c2.im; \
   (r).c2.im=(c).re*(s).c2.im+(c).im*(s).c2.re;

/* M.Hasenbusch */
#define _complexcjg_times_vector(r,c,s) \
   (r).c0.re=(c).re*(s).c0.re+(c).im*(s).c0.im; \
   (r).c0.im=(c).re*(s).c0.im-(c).im*(s).c0.re; \
   (r).c1.re=(c).re*(s).c1.re+(c).im*(s).c1.im; \
   (r).c1.im=(c).re*(s).c1.im-(c).im*(s).c1.re; \
   (r).c2.re=(c).re*(s).c2.re+(c).im*(s).c2.im; \
   (r).c2.im=(c).re*(s).c2.im-(c).im*(s).c2.re;

/*
* Real part of the scalar product (r,s)
*/

#define _vector_prod_re(r,s) \
   (r).c0.re*(s).c0.re+(r).c0.im*(s).c0.im+ \
   (r).c1.re*(s).c1.re+(r).c1.im*(s).c1.im+ \
   (r).c2.re*(s).c2.re+(r).c2.im*(s).c2.im;

/*
* Imaginary part of the scalar product (r,s)
*/

#define _vector_prod_im(r,s) \
   (r).c0.re*(s).c0.im-(r).c0.im*(s).c0.re+ \
   (r).c1.re*(s).c1.im-(r).c1.im*(s).c1.re+ \
   (r).c2.re*(s).c2.im-(r).c2.im*(s).c2.re; 

/*
* r.c1-=z*s.c1 (z of type complex)
* r.c2-=z*s.c2
* r.c3-=z*s.c3
*/

#define _vector_project(r,z,s) \
   (r).c0.re-=((z).re*(s).c0.re-(z).im*(s).c0.im); \
   (r).c0.im-=((z).re*(s).c0.im+(z).im*(s).c0.re); \
   (r).c1.re-=((z).re*(s).c1.re-(z).im*(s).c1.im); \
   (r).c1.im-=((z).re*(s).c1.im+(z).im*(s).c1.re); \
   (r).c2.re-=((z).re*(s).c2.re-(z).im*(s).c2.im); \
   (r).c2.im-=((z).re*(s).c2.im+(z).im*(s).c2.re);

/*
* SU(3) matrix u times SU(3) vector s
*  
* r.c1=(u*s).c1
* r.c2=(u*s).c2
* r.c3=(u*s).c3
*/

#if defined SSE2

#define _su3_multiply(r,u,s) \
_sse_load(s); \
_sse_su3_multiply(u); \
_sse_store_up(r);

#define _su3_inverse_multiply(r,u,s) \
_sse_load(s); \
_sse_su3_inverse_multiply(u); \
_sse_store_up(r);

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

/*
 u=1 
 added by M.Hasenbusch Thu Aug  9 10:27:28 MEST 2001 */

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

/*
 u=0 
 added by M.Hasenbusch Thu Aug  9 10:27:28 MEST 2001 */
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

/* M. Hasenbusch
* u=v
*/

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

/* M. Hasenbusch
* u=-v
*/

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

/*
* u=v^dagger
*/

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

/* M. Hasenbusch
* u=c*v 
* c real
*/

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

#if defined SSE2
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


#if defined SSE2
#define _su3_times_su3(u,v,w) _sse_su3_times_su3(u,v,w)
#define _su3_times_su3_acc(u,v,w) _sse_su3_times_su3_acc(u,v,w)
#define _su3d_times_su3(u,v,w) _sse_su3d_times_su3(u,v,w)
#define _su3d_times_su3_acc(u,v,w) _sse_su3d_times_su3_acc(u,v,w)
#define _su3_times_su3d(u,v,w) _sse_su3_times_su3d(u,v,w)
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

/* M.Hasenbusch u = u + v * v^{dag} */

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

/* M.Hasenbusch u=v * v^{dag} */

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

/* M. Hasenbusch t =u tensor v^dag */
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

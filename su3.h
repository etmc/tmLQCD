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

typedef struct
{
   double re,im;
} complex;

typedef struct 
{
   complex c11,c12,c13,c21,c22,c23,c31,c32,c33;
} su3;

typedef struct
{
   complex c1,c2,c3;
} su3_vector;

typedef struct
{
   su3_vector c1,c2,c3,c4;
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
   (r).c1.re=0.0; \
   (r).c1.im=0.0; \
   (r).c2.re=0.0; \
   (r).c2.im=0.0; \
   (r).c3.re=0.0; \
   (r).c3.im=0.0;

/* M. Hasenbusch Mon Sep 24
* r.c1=s.c1
* r.c2=s.c2
* r.c3=s.c3
*/

#define _vector_assign(r,s) \
   (r).c1.re=(s).c1.re; \
   (r).c1.im=(s).c1.im; \
   (r).c2.re=(s).c2.re; \
   (r).c2.im=(s).c2.im; \
   (r).c3.re=(s).c3.re; \
   (r).c3.im=(s).c3.im;

/* M. Hasenbusch Mon Sep 24
* r.c1=-s.c1
* r.c2=-s.c2
* r.c3=-s.c3
*/

#define _vector_minus_assign(r,s) \
   (r).c1.re=-(s).c1.re; \
   (r).c1.im=-(s).c1.im; \
   (r).c2.re=-(s).c2.re; \
   (r).c2.im=-(s).c2.im; \
   (r).c3.re=-(s).c3.re; \
   (r).c3.im=-(s).c3.im;

/*
* r.c1=c*s.c1 (c real)
* r.c2=c*s.c2
* r.c3=c*s.c3
*/

#define _vector_mul(r,c,s) \
   (r).c1.re=(c)*(s).c1.re; \
   (r).c1.im=(c)*(s).c1.im; \
   (r).c2.re=(c)*(s).c2.re; \
   (r).c2.im=(c)*(s).c2.im; \
   (r).c3.re=(c)*(s).c3.re; \
   (r).c3.im=(c)*(s).c3.im;

#define _vector_add_mul(r,c,s) \
   (r).c1.re+=(c)*(s).c1.re; \
   (r).c1.im+=(c)*(s).c1.im; \
   (r).c2.re+=(c)*(s).c2.re; \
   (r).c2.im+=(c)*(s).c2.im; \
   (r).c3.re+=(c)*(s).c3.re; \
   (r).c3.im+=(c)*(s).c3.im;

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
   (r).c1.re=(s1).c1.re+(s2).c1.re; \
   (r).c1.im=(s1).c1.im+(s2).c1.im; \
   (r).c2.re=(s1).c2.re+(s2).c2.re; \
   (r).c2.im=(s1).c2.im+(s2).c2.im; \
   (r).c3.re=(s1).c3.re+(s2).c3.re; \
   (r).c3.im=(s1).c3.im+(s2).c3.im;

/*
* r.c1=s1.c1-s2.c1
* r.c2=s1.c2-s2.c2
* r.c3=s1.c3-s2.c3
*/

#define _vector_sub(r,s1,s2) \
   (r).c1.re=(s1).c1.re-(s2).c1.re; \
   (r).c1.im=(s1).c1.im-(s2).c1.im; \
   (r).c2.re=(s1).c2.re-(s2).c2.re; \
   (r).c2.im=(s1).c2.im-(s2).c2.im; \
   (r).c3.re=(s1).c3.re-(s2).c3.re; \
   (r).c3.im=(s1).c3.im-(s2).c3.im;

#endif

/*
* r.c1=s1.c1+i*s2.c1
* r.c2=s1.c2+i*s2.c2
* r.c3=s1.c3+i*s2.c3
*/

#define _vector_i_add(r,s1,s2) \
   (r).c1.re=(s1).c1.re-(s2).c1.im; \
   (r).c1.im=(s1).c1.im+(s2).c1.re; \
   (r).c2.re=(s1).c2.re-(s2).c2.im; \
   (r).c2.im=(s1).c2.im+(s2).c2.re; \
   (r).c3.re=(s1).c3.re-(s2).c3.im; \
   (r).c3.im=(s1).c3.im+(s2).c3.re;

/*
* r.c1=s1.c1+i*s2.c1
* r.c2=s1.c2+i*s2.c2
* r.c3=s1.c3+i*s2.c3
*/

#define _vector_i_sub(r,s1,s2) \
   (r).c1.re=(s1).c1.re+(s2).c1.im; \
   (r).c1.im=(s1).c1.im-(s2).c1.re; \
   (r).c2.re=(s1).c2.re+(s2).c2.im; \
   (r).c2.im=(s1).c2.im-(s2).c2.re; \
   (r).c3.re=(s1).c3.re+(s2).c3.im; \
   (r).c3.im=(s1).c3.im-(s2).c3.re;

/*
* r.c1+=s.c1
* r.c2+=s.c2
* r.c3+=s.c3
*/

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
   (r).c1.re+=(s).c1.re; \
   (r).c1.im+=(s).c1.im; \
   (r).c2.re+=(s).c2.re; \
   (r).c2.im+=(s).c2.im; \
   (r).c3.re+=(s).c3.re; \
   (r).c3.im+=(s).c3.im;

/*
* r.c1-=s.c1
* r.c2-=s.c2
* r.c3-=s.c3
*/

#define _vector_sub_assign(r,s) \
   (r).c1.re-=(s).c1.re; \
   (r).c1.im-=(s).c1.im; \
   (r).c2.re-=(s).c2.re; \
   (r).c2.im-=(s).c2.im; \
   (r).c3.re-=(s).c3.re; \
   (r).c3.im-=(s).c3.im;

#endif 

/*
* r.c1+=i*s.c1
* r.c2+=i*s.c2
* r.c3+=i*s.c3
*/

#define _vector_i_add_assign(r,s) \
   (r).c1.re-=(s).c1.im; \
   (r).c1.im+=(s).c1.re; \
   (r).c2.re-=(s).c2.im; \
   (r).c2.im+=(s).c2.re; \
   (r).c3.re-=(s).c3.im; \
   (r).c3.im+=(s).c3.re;

/*
* r.c1-=i*s.c1
* r.c2-=i*s.c2
* r.c3-=i*s.c3
*/

#define _vector_i_sub_assign(r,s) \
   (r).c1.re+=(s).c1.im; \
   (r).c1.im-=(s).c1.re; \
   (r).c2.re+=(s).c2.im; \
   (r).c2.im-=(s).c2.re; \
   (r).c3.re+=(s).c3.im; \
   (r).c3.im-=(s).c3.re;

/* M.Hasenbusch 
* r.c1=c*s.c1
* r.c2=c*s.c2
* r.c3=c*s.c3
*
* c complex
*/

#define _complex_times_vector(r,c,s) \
   (r).c1.re=(c).re*(s).c1.re-(c).im*(s).c1.im; \
   (r).c1.im=(c).re*(s).c1.im+(c).im*(s).c1.re; \
   (r).c2.re=(c).re*(s).c2.re-(c).im*(s).c2.im; \
   (r).c2.im=(c).re*(s).c2.im+(c).im*(s).c2.re; \
   (r).c3.re=(c).re*(s).c3.re-(c).im*(s).c3.im; \
   (r).c3.im=(c).re*(s).c3.im+(c).im*(s).c3.re;

/* M.Hasenbusch */
#define _complexcjg_times_vector(r,c,s) \
   (r).c1.re=(c).re*(s).c1.re+(c).im*(s).c1.im; \
   (r).c1.im=(c).re*(s).c1.im-(c).im*(s).c1.re; \
   (r).c2.re=(c).re*(s).c2.re+(c).im*(s).c2.im; \
   (r).c2.im=(c).re*(s).c2.im-(c).im*(s).c2.re; \
   (r).c3.re=(c).re*(s).c3.re+(c).im*(s).c3.im; \
   (r).c3.im=(c).re*(s).c3.im-(c).im*(s).c3.re;

/*
* Real part of the scalar product (r,s)
*/

#define _vector_prod_re(r,s) \
   (r).c1.re*(s).c1.re+(r).c1.im*(s).c1.im+ \
   (r).c2.re*(s).c2.re+(r).c2.im*(s).c2.im+ \
   (r).c3.re*(s).c3.re+(r).c3.im*(s).c3.im;

/*
* Imaginary part of the scalar product (r,s)
*/

#define _vector_prod_im(r,s) \
   (r).c1.re*(s).c1.im-(r).c1.im*(s).c1.re+ \
   (r).c2.re*(s).c2.im-(r).c2.im*(s).c2.re+ \
   (r).c3.re*(s).c3.im-(r).c3.im*(s).c3.re; 

/*
* r.c1-=z*s.c1 (z of type complex)
* r.c2-=z*s.c2
* r.c3-=z*s.c3
*/

#define _vector_project(r,z,s) \
   (r).c1.re-=((z).re*(s).c1.re-(z).im*(s).c1.im); \
   (r).c1.im-=((z).re*(s).c1.im+(z).im*(s).c1.re); \
   (r).c2.re-=((z).re*(s).c2.re-(z).im*(s).c2.im); \
   (r).c2.im-=((z).re*(s).c2.im+(z).im*(s).c2.re); \
   (r).c3.re-=((z).re*(s).c3.re-(z).im*(s).c3.im); \
   (r).c3.im-=((z).re*(s).c3.im+(z).im*(s).c3.re);

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
   (r).c1.re= (u).c11.re*(s).c1.re-(u).c11.im*(s).c1.im  \
             +(u).c12.re*(s).c2.re-(u).c12.im*(s).c2.im  \
             +(u).c13.re*(s).c3.re-(u).c13.im*(s).c3.im; \
   (r).c1.im= (u).c11.re*(s).c1.im+(u).c11.im*(s).c1.re  \
             +(u).c12.re*(s).c2.im+(u).c12.im*(s).c2.re  \
             +(u).c13.re*(s).c3.im+(u).c13.im*(s).c3.re; \
   (r).c2.re= (u).c21.re*(s).c1.re-(u).c21.im*(s).c1.im  \
             +(u).c22.re*(s).c2.re-(u).c22.im*(s).c2.im  \
             +(u).c23.re*(s).c3.re-(u).c23.im*(s).c3.im; \
   (r).c2.im= (u).c21.re*(s).c1.im+(u).c21.im*(s).c1.re  \
             +(u).c22.re*(s).c2.im+(u).c22.im*(s).c2.re  \
             +(u).c23.re*(s).c3.im+(u).c23.im*(s).c3.re; \
   (r).c3.re= (u).c31.re*(s).c1.re-(u).c31.im*(s).c1.im  \
             +(u).c32.re*(s).c2.re-(u).c32.im*(s).c2.im  \
             +(u).c33.re*(s).c3.re-(u).c33.im*(s).c3.im; \
   (r).c3.im= (u).c31.re*(s).c1.im+(u).c31.im*(s).c1.re  \
             +(u).c32.re*(s).c2.im+(u).c32.im*(s).c2.re  \
             +(u).c33.re*(s).c3.im+(u).c33.im*(s).c3.re;

/*
* SU(3) matrix u^dagger times SU(3) vector s
*  
* r.c1=(u^dagger*s).c1
* r.c2=(u^dagger*s).c2
* r.c3=(u^dagger*s).c3
*/

#define _su3_inverse_multiply(r,u,s) \
   (r).c1.re= (u).c11.re*(s).c1.re+(u).c11.im*(s).c1.im  \
             +(u).c21.re*(s).c2.re+(u).c21.im*(s).c2.im  \
             +(u).c31.re*(s).c3.re+(u).c31.im*(s).c3.im; \
   (r).c1.im= (u).c11.re*(s).c1.im-(u).c11.im*(s).c1.re  \
             +(u).c21.re*(s).c2.im-(u).c21.im*(s).c2.re  \
             +(u).c31.re*(s).c3.im-(u).c31.im*(s).c3.re; \
   (r).c2.re= (u).c12.re*(s).c1.re+(u).c12.im*(s).c1.im  \
             +(u).c22.re*(s).c2.re+(u).c22.im*(s).c2.im  \
             +(u).c32.re*(s).c3.re+(u).c32.im*(s).c3.im; \
   (r).c2.im= (u).c12.re*(s).c1.im-(u).c12.im*(s).c1.re  \
             +(u).c22.re*(s).c2.im-(u).c22.im*(s).c2.re  \
             +(u).c32.re*(s).c3.im-(u).c32.im*(s).c3.re; \
   (r).c3.re= (u).c13.re*(s).c1.re+(u).c13.im*(s).c1.im  \
             +(u).c23.re*(s).c2.re+(u).c23.im*(s).c2.im  \
             +(u).c33.re*(s).c3.re+(u).c33.im*(s).c3.im; \
   (r).c3.im= (u).c13.re*(s).c1.im-(u).c13.im*(s).c1.re  \
             +(u).c23.re*(s).c2.im-(u).c23.im*(s).c2.re  \
             +(u).c33.re*(s).c3.im-(u).c33.im*(s).c3.re;
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
x = (u).c11.re*(u).c11.re + (u).c11.im*(u).c11.im \
   +(u).c12.re*(u).c12.re + (u).c12.im*(u).c12.im \
   +(u).c13.re*(u).c13.re + (u).c13.im*(u).c13.im \
   +(u).c21.re*(u).c21.re + (u).c21.im*(u).c21.im \
   +(u).c22.re*(u).c22.re + (u).c22.im*(u).c22.im \
   +(u).c23.re*(u).c23.re + (u).c23.im*(u).c23.im \
   +(u).c31.re*(u).c31.re + (u).c31.im*(u).c31.im \
   +(u).c32.re*(u).c32.re + (u).c32.im*(u).c32.im \
   +(u).c33.re*(u).c33.re + (u).c33.im*(u).c33.im; 

/*
 u=1 
 added by M.Hasenbusch Thu Aug  9 10:27:28 MEST 2001 */

#define _su3_one(u) \
   (u).c11.re=1.0; \
   (u).c11.im=0.0; \
   (u).c12.re=0.0; \
   (u).c12.im=0.0; \
   (u).c13.re=0.0; \
   (u).c13.im=0.0; \
   (u).c21.re=0.0; \
   (u).c21.im=0.0; \
   (u).c22.re=1.0; \
   (u).c22.im=0.0; \
   (u).c23.re=0.0; \
   (u).c23.im=0.0; \
   (u).c31.re=0.0; \
   (u).c31.im=0.0; \
   (u).c32.re=0.0; \
   (u).c32.im=0.0; \
   (u).c33.re=1.0; \
   (u).c33.im=0.0;

/*
 u=0 
 added by M.Hasenbusch Thu Aug  9 10:27:28 MEST 2001 */
#define _su3_zero(u) \
   (u).c11.re=0.0; \
   (u).c11.im=0.0; \
   (u).c12.re=0.0; \
   (u).c12.im=0.0; \
   (u).c13.re=0.0; \
   (u).c13.im=0.0; \
   (u).c21.re=0.0; \
   (u).c21.im=0.0; \
   (u).c22.re=0.0; \
   (u).c22.im=0.0; \
   (u).c23.re=0.0; \
   (u).c23.im=0.0; \
   (u).c31.re=0.0; \
   (u).c31.im=0.0; \
   (u).c32.re=0.0; \
   (u).c32.im=0.0; \
   (u).c33.re=0.0; \
   (u).c33.im=0.0;

/* M. Hasenbusch
* u=v
*/

#define _su3_assign(u,v) \
   (u).c11.re= (v).c11.re; \
   (u).c11.im= (v).c11.im; \
   (u).c12.re= (v).c12.re; \
   (u).c12.im= (v).c12.im; \
   (u).c13.re= (v).c13.re; \
   (u).c13.im= (v).c13.im; \
   (u).c21.re= (v).c21.re; \
   (u).c21.im= (v).c21.im; \
   (u).c22.re= (v).c22.re; \
   (u).c22.im= (v).c22.im; \
   (u).c23.re= (v).c23.re; \
   (u).c23.im= (v).c23.im; \
   (u).c31.re= (v).c31.re; \
   (u).c31.im= (v).c31.im; \
   (u).c32.re= (v).c32.re; \
   (u).c32.im= (v).c32.im; \
   (u).c33.re= (v).c33.re; \
   (u).c33.im= (v).c33.im;

/* M. Hasenbusch
* u=-v
*/

#define _su3_minus_assign(u,v) \
   (u).c11.re= -(v).c11.re; \
   (u).c11.im= -(v).c11.im; \
   (u).c12.re= -(v).c12.re; \
   (u).c12.im= -(v).c12.im; \
   (u).c13.re= -(v).c13.re; \
   (u).c13.im= -(v).c13.im; \
   (u).c21.re= -(v).c21.re; \
   (u).c21.im= -(v).c21.im; \
   (u).c22.re= -(v).c22.re; \
   (u).c22.im= -(v).c22.im; \
   (u).c23.re= -(v).c23.re; \
   (u).c23.im= -(v).c23.im; \
   (u).c31.re= -(v).c31.re; \
   (u).c31.im= -(v).c31.im; \
   (u).c32.re= -(v).c32.re; \
   (u).c32.im= -(v).c32.im; \
   (u).c33.re= -(v).c33.re; \
   (u).c33.im= -(v).c33.im;

/*
* u=v^dagger
*/

#define _su3_dagger(u,v) \
   (u).c11.re= (v).c11.re; \
   (u).c11.im=-(v).c11.im; \
   (u).c12.re= (v).c21.re; \
   (u).c12.im=-(v).c21.im; \
   (u).c13.re= (v).c31.re; \
   (u).c13.im=-(v).c31.im; \
   (u).c21.re= (v).c12.re; \
   (u).c21.im=-(v).c12.im; \
   (u).c22.re= (v).c22.re; \
   (u).c22.im=-(v).c22.im; \
   (u).c23.re= (v).c32.re; \
   (u).c23.im=-(v).c32.im; \
   (u).c31.re= (v).c13.re; \
   (u).c31.im=-(v).c13.im; \
   (u).c32.re= (v).c23.re; \
   (u).c32.im=-(v).c23.im; \
   (u).c33.re= (v).c33.re; \
   (u).c33.im=-(v).c33.im;

/* M.Hasenbusch */
#define _itimes_su3(u,v) \
   (u).c11.re=-(v).c11.im; \
   (u).c11.im= (v).c11.re; \
   (u).c12.re=-(v).c12.im; \
   (u).c12.im= (v).c12.re; \
   (u).c13.re=-(v).c13.im; \
   (u).c13.im= (v).c13.re; \
   (u).c21.re=-(v).c21.im; \
   (u).c21.im= (v).c21.re; \
   (u).c22.re=-(v).c22.im; \
   (u).c22.im= (v).c22.re; \
   (u).c23.re=-(v).c23.im; \
   (u).c23.im= (v).c23.re; \
   (u).c31.re=-(v).c31.im; \
   (u).c31.im= (v).c31.re; \
   (u).c32.re=-(v).c32.im; \
   (u).c32.im= (v).c32.re; \
   (u).c33.re=-(v).c33.im; \
   (u).c33.im= (v).c33.re;

/* M. Hasenbusch
* u=c*v 
* c real
*/

#define _real_times_su3(u,a,v) \
   (u).c11.re= (a)*(v).c11.re; \
   (u).c11.im= (a)*(v).c11.im; \
   (u).c12.re= (a)*(v).c12.re; \
   (u).c12.im= (a)*(v).c12.im; \
   (u).c13.re= (a)*(v).c13.re; \
   (u).c13.im= (a)*(v).c13.im; \
   (u).c21.re= (a)*(v).c21.re; \
   (u).c21.im= (a)*(v).c21.im; \
   (u).c22.re= (a)*(v).c22.re; \
   (u).c22.im= (a)*(v).c22.im; \
   (u).c23.re= (a)*(v).c23.re; \
   (u).c23.im= (a)*(v).c23.im; \
   (u).c31.re= (a)*(v).c31.re; \
   (u).c31.im= (a)*(v).c31.im; \
   (u).c32.re= (a)*(v).c32.re; \
   (u).c32.im= (a)*(v).c32.im; \
   (u).c33.re= (a)*(v).c33.re; \
   (u).c33.im= (a)*(v).c33.im;

/* M. Hasenbusch
* u=v-w
*/

#define _su3_minus_su3(u,v,w) \
   (u).c11.re= (v).c11.re-(w).c11.re; \
   (u).c11.im= (v).c11.im-(w).c11.im; \
   (u).c12.re= (v).c12.re-(w).c12.re; \
   (u).c12.im= (v).c12.im-(w).c12.im; \
   (u).c13.re= (v).c13.re-(w).c13.re; \
   (u).c13.im= (v).c13.im-(w).c13.im; \
   (u).c21.re= (v).c21.re-(w).c21.re; \
   (u).c21.im= (v).c21.im-(w).c21.im; \
   (u).c22.re= (v).c22.re-(w).c22.re; \
   (u).c22.im= (v).c22.im-(w).c22.im; \
   (u).c23.re= (v).c23.re-(w).c23.re; \
   (u).c23.im= (v).c23.im-(w).c23.im; \
   (u).c31.re= (v).c31.re-(w).c31.re; \
   (u).c31.im= (v).c31.im-(w).c31.im; \
   (u).c32.re= (v).c32.re-(w).c32.re; \
   (u).c32.im= (v).c32.im-(w).c32.im; \
   (u).c33.re= (v).c33.re-(w).c33.re; \
   (u).c33.im= (v).c33.im-(w).c33.im;

/* M. Hasenbusch
* u=i*(v-w)
*/

#define _itimes_su3_minus_su3(u,v,w) \
   (u).c11.im= (v).c11.re-(w).c11.re; \
   (u).c11.re= (w).c11.im-(v).c11.im; \
   (u).c12.im= (v).c12.re-(w).c12.re; \
   (u).c12.re= (w).c12.im-(v).c12.im; \
   (u).c13.im= (v).c13.re-(w).c13.re; \
   (u).c13.re= (w).c13.im-(v).c13.im; \
   (u).c21.im= (v).c21.re-(w).c21.re; \
   (u).c21.re= (w).c21.im-(v).c21.im; \
   (u).c22.im= (v).c22.re-(w).c22.re; \
   (u).c22.re= (w).c22.im-(v).c22.im; \
   (u).c23.im= (v).c23.re-(w).c23.re; \
   (u).c23.re= (w).c23.im-(v).c23.im; \
   (u).c31.im= (v).c31.re-(w).c31.re; \
   (u).c31.re= (w).c31.im-(v).c31.im; \
   (u).c32.im= (v).c32.re-(w).c32.re; \
   (u).c32.re= (w).c32.im-(v).c32.im; \
   (u).c33.im= (v).c33.re-(w).c33.re; \
   (u).c33.re= (w).c33.im-(v).c33.im;

/* M. Hasenbusch
* u=v+w
*/

#define _su3_plus_su3(u,v,w) \
   (u).c11.re= (v).c11.re+(w).c11.re; \
   (u).c11.im= (v).c11.im+(w).c11.im; \
   (u).c12.re= (v).c12.re+(w).c12.re; \
   (u).c12.im= (v).c12.im+(w).c12.im; \
   (u).c13.re= (v).c13.re+(w).c13.re; \
   (u).c13.im= (v).c13.im+(w).c13.im; \
   (u).c21.re= (v).c21.re+(w).c21.re; \
   (u).c21.im= (v).c21.im+(w).c21.im; \
   (u).c22.re= (v).c22.re+(w).c22.re; \
   (u).c22.im= (v).c22.im+(w).c22.im; \
   (u).c23.re= (v).c23.re+(w).c23.re; \
   (u).c23.im= (v).c23.im+(w).c23.im; \
   (u).c31.re= (v).c31.re+(w).c31.re; \
   (u).c31.im= (v).c31.im+(w).c31.im; \
   (u).c32.re= (v).c32.re+(w).c32.re; \
   (u).c32.im= (v).c32.im+(w).c32.im; \
   (u).c33.re= (v).c33.re+(w).c33.re; \
   (u).c33.im= (v).c33.im+(w).c33.im;

/* M. Hasenbusch
* u=-(v+w)
*/

#define _minus_su3_plus_su3(u,v,w) \
   (u).c11.re=-(v).c11.re-(w).c11.re; \
   (u).c11.im=-(v).c11.im-(w).c11.im; \
   (u).c12.re=-(v).c12.re-(w).c12.re; \
   (u).c12.im=-(v).c12.im-(w).c12.im; \
   (u).c13.re=-(v).c13.re-(w).c13.re; \
   (u).c13.im=-(v).c13.im-(w).c13.im; \
   (u).c21.re=-(v).c21.re-(w).c21.re; \
   (u).c21.im=-(v).c21.im-(w).c21.im; \
   (u).c22.re=-(v).c22.re-(w).c22.re; \
   (u).c22.im=-(v).c22.im-(w).c22.im; \
   (u).c23.re=-(v).c23.re-(w).c23.re; \
   (u).c23.im=-(v).c23.im-(w).c23.im; \
   (u).c31.re=-(v).c31.re-(w).c31.re; \
   (u).c31.im=-(v).c31.im-(w).c31.im; \
   (u).c32.re=-(v).c32.re-(w).c32.re; \
   (u).c32.im=-(v).c32.im-(w).c32.im; \
   (u).c33.re=-(v).c33.re-(w).c33.re; \
   (u).c33.im=-(v).c33.im-(w).c33.im;

/* M. Hasenbusch
* u=i*(v+w)
*/

#define _itimes_su3_plus_su3(u,v,w) \
   (u).c11.im= (v).c11.re+(w).c11.re; \
   (u).c11.re=-(v).c11.im-(w).c11.im; \
   (u).c12.im= (v).c12.re+(w).c12.re; \
   (u).c12.re=-(v).c12.im-(w).c12.im; \
   (u).c13.im= (v).c13.re+(w).c13.re; \
   (u).c13.re=-(v).c13.im-(w).c13.im; \
   (u).c21.im= (v).c21.re+(w).c21.re; \
   (u).c21.re=-(v).c21.im-(w).c21.im; \
   (u).c22.im= (v).c22.re+(w).c22.re; \
   (u).c22.re=-(v).c22.im-(w).c22.im; \
   (u).c23.im= (v).c23.re+(w).c23.re; \
   (u).c23.re=-(v).c23.im-(w).c23.im; \
   (u).c31.im= (v).c31.re+(w).c31.re; \
   (u).c31.re=-(v).c31.im-(w).c31.im; \
   (u).c32.im= (v).c32.re+(w).c32.re; \
   (u).c32.re=-(v).c32.im-(w).c32.im; \
   (u).c33.im= (v).c33.re+(w).c33.re; \
   (u).c33.re=-(v).c33.im-(w).c33.im;

/* M. Hasenbusch
* u=-i*(v+w)
*/

#define _minus_itimes_su3_plus_su3(u,v,w) \
   (u).c11.im=-(v).c11.re-(w).c11.re; \
   (u).c11.re= (v).c11.im+(w).c11.im; \
   (u).c12.im=-(v).c12.re-(w).c12.re; \
   (u).c12.re= (v).c12.im+(w).c12.im; \
   (u).c13.im=-(v).c13.re-(w).c13.re; \
   (u).c13.re= (v).c13.im+(w).c13.im; \
   (u).c21.im=-(v).c21.re-(w).c21.re; \
   (u).c21.re= (v).c21.im+(w).c21.im; \
   (u).c22.im=-(v).c22.re-(w).c22.re; \
   (u).c22.re= (v).c22.im+(w).c22.im; \
   (u).c23.im=-(v).c23.re-(w).c23.re; \
   (u).c23.re= (v).c23.im+(w).c23.im; \
   (u).c31.im=-(v).c31.re-(w).c31.re; \
   (u).c31.re= (v).c31.im+(w).c31.im; \
   (u).c32.im=-(v).c32.re-(w).c32.re; \
   (u).c32.re= (v).c32.im+(w).c32.im; \
   (u).c33.im=-(v).c33.re-(w).c33.re; \
   (u).c33.re= (v).c33.im+(w).c33.im;

/* M.Hasenbusch */
#define _complex_times_su3(r,c,s) \
   (r).c11.re=(c).re*(s).c11.re-(c).im*(s).c11.im; \
   (r).c11.im=(c).re*(s).c11.im+(c).im*(s).c11.re; \
   (r).c12.re=(c).re*(s).c12.re-(c).im*(s).c12.im; \
   (r).c12.im=(c).re*(s).c12.im+(c).im*(s).c12.re; \
   (r).c13.re=(c).re*(s).c13.re-(c).im*(s).c13.im; \
   (r).c13.im=(c).re*(s).c13.im+(c).im*(s).c13.re; \
   (r).c21.re=(c).re*(s).c21.re-(c).im*(s).c21.im; \
   (r).c21.im=(c).re*(s).c21.im+(c).im*(s).c21.re; \
   (r).c22.re=(c).re*(s).c22.re-(c).im*(s).c22.im; \
   (r).c22.im=(c).re*(s).c22.im+(c).im*(s).c22.re; \
   (r).c23.re=(c).re*(s).c23.re-(c).im*(s).c23.im; \
   (r).c23.im=(c).re*(s).c23.im+(c).im*(s).c23.re; \
   (r).c31.re=(c).re*(s).c31.re-(c).im*(s).c31.im; \
   (r).c31.im=(c).re*(s).c31.im+(c).im*(s).c31.re; \
   (r).c32.re=(c).re*(s).c32.re-(c).im*(s).c32.im; \
   (r).c32.im=(c).re*(s).c32.im+(c).im*(s).c32.re; \
   (r).c33.re=(c).re*(s).c33.re-(c).im*(s).c33.im; \
   (r).c33.im=(c).re*(s).c33.im+(c).im*(s).c33.re; 

/* M.Hasenbusch */
#define _complexcjg_times_su3(r,c,s) \
   (r).c11.re=(c).re*(s).c11.re+(c).im*(s).c11.im; \
   (r).c11.im=(c).re*(s).c11.im-(c).im*(s).c11.re; \
   (r).c12.re=(c).re*(s).c12.re+(c).im*(s).c12.im; \
   (r).c12.im=(c).re*(s).c12.im-(c).im*(s).c12.re; \
   (r).c13.re=(c).re*(s).c13.re+(c).im*(s).c13.im; \
   (r).c13.im=(c).re*(s).c13.im-(c).im*(s).c13.re; \
   (r).c21.re=(c).re*(s).c21.re+(c).im*(s).c21.im; \
   (r).c21.im=(c).re*(s).c21.im-(c).im*(s).c21.re; \
   (r).c22.re=(c).re*(s).c22.re+(c).im*(s).c22.im; \
   (r).c22.im=(c).re*(s).c22.im-(c).im*(s).c22.re; \
   (r).c23.re=(c).re*(s).c23.re+(c).im*(s).c23.im; \
   (r).c23.im=(c).re*(s).c23.im-(c).im*(s).c23.re; \
   (r).c31.re=(c).re*(s).c31.re+(c).im*(s).c31.im; \
   (r).c31.im=(c).re*(s).c31.im-(c).im*(s).c31.re; \
   (r).c32.re=(c).re*(s).c32.re+(c).im*(s).c32.im; \
   (r).c32.im=(c).re*(s).c32.im-(c).im*(s).c32.re; \
   (r).c33.re=(c).re*(s).c33.re+(c).im*(s).c33.im; \
   (r).c33.im=(c).re*(s).c33.im-(c).im*(s).c33.re;


/* M. Hasenbusch
* su3_acc
*/

#if defined SSE2
#define _su3_acc(u,v) _sse_su3_acc(u,v) 
#else
#define _su3_acc(u,v) \
   (u).c11.re+=(v).c11.re; \
   (u).c11.im+=(v).c11.im; \
   (u).c12.re+=(v).c12.re; \
   (u).c12.im+=(v).c12.im; \
   (u).c13.re+=(v).c13.re; \
   (u).c13.im+=(v).c13.im; \
   (u).c21.re+=(v).c21.re; \
   (u).c21.im+=(v).c21.im; \
   (u).c22.re+=(v).c22.re; \
   (u).c22.im+=(v).c22.im; \
   (u).c23.re+=(v).c23.re; \
   (u).c23.im+=(v).c23.im; \
   (u).c31.re+=(v).c31.re; \
   (u).c31.im+=(v).c31.im; \
   (u).c32.re+=(v).c32.re; \
   (u).c32.im+=(v).c32.im; \
   (u).c33.re+=(v).c33.re; \
   (u).c33.im+=(v).c33.im;
#endif

/*
* su3_refac_acc
*/

#define _su3_refac_acc(u,a,v) \
   (u).c11.re+=a*(v).c11.re; \
   (u).c11.im+=a*(v).c11.im; \
   (u).c12.re+=a*(v).c12.re; \
   (u).c12.im+=a*(v).c12.im; \
   (u).c13.re+=a*(v).c13.re; \
   (u).c13.im+=a*(v).c13.im; \
   (u).c21.re+=a*(v).c21.re; \
   (u).c21.im+=a*(v).c21.im; \
   (u).c22.re+=a*(v).c22.re; \
   (u).c22.im+=a*(v).c22.im; \
   (u).c23.re+=a*(v).c23.re; \
   (u).c23.im+=a*(v).c23.im; \
   (u).c31.re+=a*(v).c31.re; \
   (u).c31.im+=a*(v).c31.im; \
   (u).c32.re+=a*(v).c32.re; \
   (u).c32.im+=a*(v).c32.im; \
   (u).c33.re+=a*(v).c33.re; \
   (u).c33.im+=a*(v).c33.im;
/*
* su3_imfac_acc
*/

#define _su3_imfac_acc(u,a,v) \
   (u).c11.re-=a*(v).c11.im; \
   (u).c11.im+=a*(v).c11.re; \
   (u).c12.re-=a*(v).c12.im; \
   (u).c12.im+=a*(v).c12.re; \
   (u).c13.re-=a*(v).c13.im; \
   (u).c13.im+=a*(v).c13.re; \
   (u).c21.re-=a*(v).c21.im; \
   (u).c21.im+=a*(v).c21.re; \
   (u).c22.re-=a*(v).c22.im; \
   (u).c22.im+=a*(v).c22.re; \
   (u).c23.re-=a*(v).c23.im; \
   (u).c23.im+=a*(v).c23.re; \
   (u).c31.re-=a*(v).c31.im; \
   (u).c31.im+=a*(v).c31.re; \
   (u).c32.re-=a*(v).c32.im; \
   (u).c32.im+=a*(v).c32.re; \
   (u).c33.re-=a*(v).c33.im; \
   (u).c33.im+=a*(v).c33.re;


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
   (u).c11.re= (v).c11.re*(w).c11.re-(v).c11.im*(w).c11.im  \
              +(v).c12.re*(w).c21.re-(v).c12.im*(w).c21.im  \
              +(v).c13.re*(w).c31.re-(v).c13.im*(w).c31.im; \
   (u).c11.im= (v).c11.re*(w).c11.im+(v).c11.im*(w).c11.re  \
              +(v).c12.re*(w).c21.im+(v).c12.im*(w).c21.re  \
              +(v).c13.re*(w).c31.im+(v).c13.im*(w).c31.re; \
   (u).c12.re= (v).c11.re*(w).c12.re-(v).c11.im*(w).c12.im  \
              +(v).c12.re*(w).c22.re-(v).c12.im*(w).c22.im  \
              +(v).c13.re*(w).c32.re-(v).c13.im*(w).c32.im; \
   (u).c12.im= (v).c11.re*(w).c12.im+(v).c11.im*(w).c12.re  \
              +(v).c12.re*(w).c22.im+(v).c12.im*(w).c22.re  \
              +(v).c13.re*(w).c32.im+(v).c13.im*(w).c32.re; \
   (u).c13.re= (v).c11.re*(w).c13.re-(v).c11.im*(w).c13.im  \
              +(v).c12.re*(w).c23.re-(v).c12.im*(w).c23.im  \
              +(v).c13.re*(w).c33.re-(v).c13.im*(w).c33.im; \
   (u).c13.im= (v).c11.re*(w).c13.im+(v).c11.im*(w).c13.re  \
              +(v).c12.re*(w).c23.im+(v).c12.im*(w).c23.re  \
              +(v).c13.re*(w).c33.im+(v).c13.im*(w).c33.re; \
   (u).c21.re= (v).c21.re*(w).c11.re-(v).c21.im*(w).c11.im  \
              +(v).c22.re*(w).c21.re-(v).c22.im*(w).c21.im  \
              +(v).c23.re*(w).c31.re-(v).c23.im*(w).c31.im; \
   (u).c21.im= (v).c21.re*(w).c11.im+(v).c21.im*(w).c11.re  \
              +(v).c22.re*(w).c21.im+(v).c22.im*(w).c21.re  \
              +(v).c23.re*(w).c31.im+(v).c23.im*(w).c31.re; \
   (u).c22.re= (v).c21.re*(w).c12.re-(v).c21.im*(w).c12.im  \
              +(v).c22.re*(w).c22.re-(v).c22.im*(w).c22.im  \
              +(v).c23.re*(w).c32.re-(v).c23.im*(w).c32.im; \
   (u).c22.im= (v).c21.re*(w).c12.im+(v).c21.im*(w).c12.re  \
              +(v).c22.re*(w).c22.im+(v).c22.im*(w).c22.re  \
              +(v).c23.re*(w).c32.im+(v).c23.im*(w).c32.re; \
   (u).c23.re= (v).c21.re*(w).c13.re-(v).c21.im*(w).c13.im  \
              +(v).c22.re*(w).c23.re-(v).c22.im*(w).c23.im  \
              +(v).c23.re*(w).c33.re-(v).c23.im*(w).c33.im; \
   (u).c23.im= (v).c21.re*(w).c13.im+(v).c21.im*(w).c13.re  \
              +(v).c22.re*(w).c23.im+(v).c22.im*(w).c23.re  \
              +(v).c23.re*(w).c33.im+(v).c23.im*(w).c33.re; \
   (u).c31.re= (v).c31.re*(w).c11.re-(v).c31.im*(w).c11.im  \
              +(v).c32.re*(w).c21.re-(v).c32.im*(w).c21.im  \
              +(v).c33.re*(w).c31.re-(v).c33.im*(w).c31.im; \
   (u).c31.im= (v).c31.re*(w).c11.im+(v).c31.im*(w).c11.re  \
              +(v).c32.re*(w).c21.im+(v).c32.im*(w).c21.re  \
              +(v).c33.re*(w).c31.im+(v).c33.im*(w).c31.re; \
   (u).c32.re= (v).c31.re*(w).c12.re-(v).c31.im*(w).c12.im  \
              +(v).c32.re*(w).c22.re-(v).c32.im*(w).c22.im  \
              +(v).c33.re*(w).c32.re-(v).c33.im*(w).c32.im; \
   (u).c32.im= (v).c31.re*(w).c12.im+(v).c31.im*(w).c12.re  \
              +(v).c32.re*(w).c22.im+(v).c32.im*(w).c22.re  \
              +(v).c33.re*(w).c32.im+(v).c33.im*(w).c32.re; \
   (u).c33.re= (v).c31.re*(w).c13.re-(v).c31.im*(w).c13.im  \
              +(v).c32.re*(w).c23.re-(v).c32.im*(w).c23.im  \
              +(v).c33.re*(w).c33.re-(v).c33.im*(w).c33.im; \
   (u).c33.im= (v).c31.re*(w).c13.im+(v).c31.im*(w).c13.re  \
              +(v).c32.re*(w).c23.im+(v).c32.im*(w).c23.re  \
              +(v).c33.re*(w).c33.im+(v).c33.im*(w).c33.re;

/* M.Hasenbusch u = u + v * v^{dag} */

#define _su3_times_su3_acc(u,v,w) \
   (u).c11.re+=(v).c11.re*(w).c11.re-(v).c11.im*(w).c11.im  \
              +(v).c12.re*(w).c21.re-(v).c12.im*(w).c21.im  \
              +(v).c13.re*(w).c31.re-(v).c13.im*(w).c31.im; \
   (u).c11.im+=(v).c11.re*(w).c11.im+(v).c11.im*(w).c11.re  \
              +(v).c12.re*(w).c21.im+(v).c12.im*(w).c21.re  \
              +(v).c13.re*(w).c31.im+(v).c13.im*(w).c31.re; \
   (u).c12.re+=(v).c11.re*(w).c12.re-(v).c11.im*(w).c12.im  \
              +(v).c12.re*(w).c22.re-(v).c12.im*(w).c22.im  \
              +(v).c13.re*(w).c32.re-(v).c13.im*(w).c32.im; \
   (u).c12.im+=(v).c11.re*(w).c12.im+(v).c11.im*(w).c12.re  \
              +(v).c12.re*(w).c22.im+(v).c12.im*(w).c22.re  \
              +(v).c13.re*(w).c32.im+(v).c13.im*(w).c32.re; \
   (u).c13.re+=(v).c11.re*(w).c13.re-(v).c11.im*(w).c13.im  \
              +(v).c12.re*(w).c23.re-(v).c12.im*(w).c23.im  \
              +(v).c13.re*(w).c33.re-(v).c13.im*(w).c33.im; \
   (u).c13.im+=(v).c11.re*(w).c13.im+(v).c11.im*(w).c13.re  \
              +(v).c12.re*(w).c23.im+(v).c12.im*(w).c23.re  \
              +(v).c13.re*(w).c33.im+(v).c13.im*(w).c33.re; \
   (u).c21.re+=(v).c21.re*(w).c11.re-(v).c21.im*(w).c11.im  \
              +(v).c22.re*(w).c21.re-(v).c22.im*(w).c21.im  \
              +(v).c23.re*(w).c31.re-(v).c23.im*(w).c31.im; \
   (u).c21.im+=(v).c21.re*(w).c11.im+(v).c21.im*(w).c11.re  \
              +(v).c22.re*(w).c21.im+(v).c22.im*(w).c21.re  \
              +(v).c23.re*(w).c31.im+(v).c23.im*(w).c31.re; \
   (u).c22.re+=(v).c21.re*(w).c12.re-(v).c21.im*(w).c12.im  \
              +(v).c22.re*(w).c22.re-(v).c22.im*(w).c22.im  \
              +(v).c23.re*(w).c32.re-(v).c23.im*(w).c32.im; \
   (u).c22.im+=(v).c21.re*(w).c12.im+(v).c21.im*(w).c12.re  \
              +(v).c22.re*(w).c22.im+(v).c22.im*(w).c22.re  \
              +(v).c23.re*(w).c32.im+(v).c23.im*(w).c32.re; \
   (u).c23.re+=(v).c21.re*(w).c13.re-(v).c21.im*(w).c13.im  \
              +(v).c22.re*(w).c23.re-(v).c22.im*(w).c23.im  \
              +(v).c23.re*(w).c33.re-(v).c23.im*(w).c33.im; \
   (u).c23.im+=(v).c21.re*(w).c13.im+(v).c21.im*(w).c13.re  \
              +(v).c22.re*(w).c23.im+(v).c22.im*(w).c23.re  \
              +(v).c23.re*(w).c33.im+(v).c23.im*(w).c33.re; \
   (u).c31.re+=(v).c31.re*(w).c11.re-(v).c31.im*(w).c11.im  \
              +(v).c32.re*(w).c21.re-(v).c32.im*(w).c21.im  \
              +(v).c33.re*(w).c31.re-(v).c33.im*(w).c31.im; \
   (u).c31.im+=(v).c31.re*(w).c11.im+(v).c31.im*(w).c11.re  \
              +(v).c32.re*(w).c21.im+(v).c32.im*(w).c21.re  \
              +(v).c33.re*(w).c31.im+(v).c33.im*(w).c31.re; \
   (u).c32.re+=(v).c31.re*(w).c12.re-(v).c31.im*(w).c12.im  \
              +(v).c32.re*(w).c22.re-(v).c32.im*(w).c22.im  \
              +(v).c33.re*(w).c32.re-(v).c33.im*(w).c32.im; \
   (u).c32.im+=(v).c31.re*(w).c12.im+(v).c31.im*(w).c12.re  \
              +(v).c32.re*(w).c22.im+(v).c32.im*(w).c22.re  \
              +(v).c33.re*(w).c32.im+(v).c33.im*(w).c32.re; \
   (u).c33.re+=(v).c31.re*(w).c13.re-(v).c31.im*(w).c13.im  \
              +(v).c32.re*(w).c23.re-(v).c32.im*(w).c23.im  \
              +(v).c33.re*(w).c33.re-(v).c33.im*(w).c33.im; \
   (u).c33.im+=(v).c31.re*(w).c13.im+(v).c31.im*(w).c13.re  \
              +(v).c32.re*(w).c23.im+(v).c32.im*(w).c23.re  \
              +(v).c33.re*(w).c33.im+(v).c33.im*(w).c33.re;

/* M.Hasenbusch u=v * v^{dag} */

#define _su3_times_su3d(u,v,w) \
   (u).c11.re= (v).c11.re*(w).c11.re+(v).c11.im*(w).c11.im  \
              +(v).c12.re*(w).c12.re+(v).c12.im*(w).c12.im  \
              +(v).c13.re*(w).c13.re+(v).c13.im*(w).c13.im; \
   (u).c11.im=-(v).c11.re*(w).c11.im+(v).c11.im*(w).c11.re  \
              -(v).c12.re*(w).c12.im+(v).c12.im*(w).c12.re  \
              -(v).c13.re*(w).c13.im+(v).c13.im*(w).c13.re; \
   (u).c12.re= (v).c11.re*(w).c21.re+(v).c11.im*(w).c21.im  \
              +(v).c12.re*(w).c22.re+(v).c12.im*(w).c22.im  \
              +(v).c13.re*(w).c23.re+(v).c13.im*(w).c23.im; \
   (u).c12.im=-(v).c11.re*(w).c21.im+(v).c11.im*(w).c21.re  \
              -(v).c12.re*(w).c22.im+(v).c12.im*(w).c22.re  \
              -(v).c13.re*(w).c23.im+(v).c13.im*(w).c23.re; \
   (u).c13.re= (v).c11.re*(w).c31.re+(v).c11.im*(w).c31.im  \
              +(v).c12.re*(w).c32.re+(v).c12.im*(w).c32.im  \
              +(v).c13.re*(w).c33.re+(v).c13.im*(w).c33.im; \
   (u).c13.im=-(v).c11.re*(w).c31.im+(v).c11.im*(w).c31.re  \
              -(v).c12.re*(w).c32.im+(v).c12.im*(w).c32.re  \
              -(v).c13.re*(w).c33.im+(v).c13.im*(w).c33.re; \
   (u).c21.re= (v).c21.re*(w).c11.re+(v).c21.im*(w).c11.im  \
              +(v).c22.re*(w).c12.re+(v).c22.im*(w).c12.im  \
              +(v).c23.re*(w).c13.re+(v).c23.im*(w).c13.im; \
   (u).c21.im=-(v).c21.re*(w).c11.im+(v).c21.im*(w).c11.re  \
              -(v).c22.re*(w).c12.im+(v).c22.im*(w).c12.re  \
              -(v).c23.re*(w).c13.im+(v).c23.im*(w).c13.re; \
   (u).c22.re= (v).c21.re*(w).c21.re+(v).c21.im*(w).c21.im  \
              +(v).c22.re*(w).c22.re+(v).c22.im*(w).c22.im  \
              +(v).c23.re*(w).c23.re+(v).c23.im*(w).c23.im; \
   (u).c22.im=-(v).c21.re*(w).c21.im+(v).c21.im*(w).c21.re  \
              -(v).c22.re*(w).c22.im+(v).c22.im*(w).c22.re  \
              -(v).c23.re*(w).c23.im+(v).c23.im*(w).c23.re; \
   (u).c23.re= (v).c21.re*(w).c31.re+(v).c21.im*(w).c31.im  \
              +(v).c22.re*(w).c32.re+(v).c22.im*(w).c32.im  \
              +(v).c23.re*(w).c33.re+(v).c23.im*(w).c33.im; \
   (u).c23.im=-(v).c21.re*(w).c31.im+(v).c21.im*(w).c31.re  \
              -(v).c22.re*(w).c32.im+(v).c22.im*(w).c32.re  \
              -(v).c23.re*(w).c33.im+(v).c23.im*(w).c33.re; \
   (u).c31.re= (v).c31.re*(w).c11.re+(v).c31.im*(w).c11.im  \
              +(v).c32.re*(w).c12.re+(v).c32.im*(w).c12.im  \
              +(v).c33.re*(w).c13.re+(v).c33.im*(w).c13.im; \
   (u).c31.im=-(v).c31.re*(w).c11.im+(v).c31.im*(w).c11.re  \
              -(v).c32.re*(w).c12.im+(v).c32.im*(w).c12.re  \
              -(v).c33.re*(w).c13.im+(v).c33.im*(w).c13.re; \
   (u).c32.re= (v).c31.re*(w).c21.re+(v).c31.im*(w).c21.im  \
              +(v).c32.re*(w).c22.re+(v).c32.im*(w).c22.im  \
              +(v).c33.re*(w).c23.re+(v).c33.im*(w).c23.im; \
   (u).c32.im=-(v).c31.re*(w).c21.im+(v).c31.im*(w).c21.re  \
              -(v).c32.re*(w).c22.im+(v).c32.im*(w).c22.re  \
              -(v).c33.re*(w).c23.im+(v).c33.im*(w).c23.re; \
   (u).c33.re= (v).c31.re*(w).c31.re+(v).c31.im*(w).c31.im  \
              +(v).c32.re*(w).c32.re+(v).c32.im*(w).c32.im  \
              +(v).c33.re*(w).c33.re+(v).c33.im*(w).c33.im; \
   (u).c33.im=-(v).c31.re*(w).c31.im+(v).c31.im*(w).c31.re  \
              -(v).c32.re*(w).c32.im+(v).c32.im*(w).c32.re  \
              -(v).c33.re*(w).c33.im+(v).c33.im*(w).c33.re;

/* M.Hasenbusch u=v^{dag} w */

#define _su3d_times_su3(u,v,w) \
   (u).c11.re= (v).c11.re*(w).c11.re+(v).c11.im*(w).c11.im  \
              +(v).c21.re*(w).c21.re+(v).c21.im*(w).c21.im  \
              +(v).c31.re*(w).c31.re+(v).c31.im*(w).c31.im; \
   (u).c11.im= (v).c11.re*(w).c11.im-(v).c11.im*(w).c11.re  \
              +(v).c21.re*(w).c21.im-(v).c21.im*(w).c21.re  \
              +(v).c31.re*(w).c31.im-(v).c31.im*(w).c31.re; \
   (u).c12.re= (v).c11.re*(w).c12.re+(v).c11.im*(w).c12.im  \
              +(v).c21.re*(w).c22.re+(v).c21.im*(w).c22.im  \
              +(v).c31.re*(w).c32.re+(v).c31.im*(w).c32.im; \
   (u).c12.im= (v).c11.re*(w).c12.im-(v).c11.im*(w).c12.re  \
              +(v).c21.re*(w).c22.im-(v).c21.im*(w).c22.re  \
              +(v).c31.re*(w).c32.im-(v).c31.im*(w).c32.re; \
   (u).c13.re= (v).c11.re*(w).c13.re+(v).c11.im*(w).c13.im  \
              +(v).c21.re*(w).c23.re+(v).c21.im*(w).c23.im  \
              +(v).c31.re*(w).c33.re+(v).c31.im*(w).c33.im; \
   (u).c13.im= (v).c11.re*(w).c13.im-(v).c11.im*(w).c13.re  \
              +(v).c21.re*(w).c23.im-(v).c21.im*(w).c23.re  \
              +(v).c31.re*(w).c33.im-(v).c31.im*(w).c33.re; \
   (u).c21.re= (v).c12.re*(w).c11.re+(v).c12.im*(w).c11.im  \
              +(v).c22.re*(w).c21.re+(v).c22.im*(w).c21.im  \
              +(v).c32.re*(w).c31.re+(v).c32.im*(w).c31.im; \
   (u).c21.im= (v).c12.re*(w).c11.im-(v).c12.im*(w).c11.re  \
              +(v).c22.re*(w).c21.im-(v).c22.im*(w).c21.re  \
              +(v).c32.re*(w).c31.im-(v).c32.im*(w).c31.re; \
   (u).c22.re= (v).c12.re*(w).c12.re+(v).c12.im*(w).c12.im  \
              +(v).c22.re*(w).c22.re+(v).c22.im*(w).c22.im  \
              +(v).c32.re*(w).c32.re+(v).c32.im*(w).c32.im; \
   (u).c22.im= (v).c12.re*(w).c12.im-(v).c12.im*(w).c12.re  \
              +(v).c22.re*(w).c22.im-(v).c22.im*(w).c22.re  \
              +(v).c32.re*(w).c32.im-(v).c32.im*(w).c32.re; \
   (u).c23.re= (v).c12.re*(w).c13.re+(v).c12.im*(w).c13.im  \
              +(v).c22.re*(w).c23.re+(v).c22.im*(w).c23.im  \
              +(v).c32.re*(w).c33.re+(v).c32.im*(w).c33.im; \
   (u).c23.im= (v).c12.re*(w).c13.im-(v).c12.im*(w).c13.re  \
              +(v).c22.re*(w).c23.im-(v).c22.im*(w).c23.re  \
              +(v).c32.re*(w).c33.im-(v).c32.im*(w).c33.re; \
   (u).c31.re= (v).c13.re*(w).c11.re+(v).c13.im*(w).c11.im  \
              +(v).c23.re*(w).c21.re+(v).c23.im*(w).c21.im  \
              +(v).c33.re*(w).c31.re+(v).c33.im*(w).c31.im; \
   (u).c31.im= (v).c13.re*(w).c11.im-(v).c13.im*(w).c11.re  \
              +(v).c23.re*(w).c21.im-(v).c23.im*(w).c21.re  \
              +(v).c33.re*(w).c31.im-(v).c33.im*(w).c31.re; \
   (u).c32.re= (v).c13.re*(w).c12.re+(v).c13.im*(w).c12.im  \
              +(v).c23.re*(w).c22.re+(v).c23.im*(w).c22.im  \
              +(v).c33.re*(w).c32.re+(v).c33.im*(w).c32.im; \
   (u).c32.im= (v).c13.re*(w).c12.im-(v).c13.im*(w).c12.re  \
              +(v).c23.re*(w).c22.im-(v).c23.im*(w).c22.re  \
              +(v).c33.re*(w).c32.im-(v).c33.im*(w).c32.re; \
   (u).c33.re= (v).c13.re*(w).c13.re+(v).c13.im*(w).c13.im  \
              +(v).c23.re*(w).c23.re+(v).c23.im*(w).c23.im  \
              +(v).c33.re*(w).c33.re+(v).c33.im*(w).c33.im; \
   (u).c33.im= (v).c13.re*(w).c13.im-(v).c13.im*(w).c13.re  \
              +(v).c23.re*(w).c23.im-(v).c23.im*(w).c23.re  \
              +(v).c33.re*(w).c33.im-(v).c33.im*(w).c33.re;

#define _su3d_times_su3_acc(u,v,w) \
   (u).c11.re+=(v).c11.re*(w).c11.re+(v).c11.im*(w).c11.im  \
              +(v).c21.re*(w).c21.re+(v).c21.im*(w).c21.im  \
              +(v).c31.re*(w).c31.re+(v).c31.im*(w).c31.im; \
   (u).c11.im+=(v).c11.re*(w).c11.im-(v).c11.im*(w).c11.re  \
              +(v).c21.re*(w).c21.im-(v).c21.im*(w).c21.re  \
              +(v).c31.re*(w).c31.im-(v).c31.im*(w).c31.re; \
   (u).c12.re+=(v).c11.re*(w).c12.re+(v).c11.im*(w).c12.im  \
              +(v).c21.re*(w).c22.re+(v).c21.im*(w).c22.im  \
              +(v).c31.re*(w).c32.re+(v).c31.im*(w).c32.im; \
   (u).c12.im+=(v).c11.re*(w).c12.im-(v).c11.im*(w).c12.re  \
              +(v).c21.re*(w).c22.im-(v).c21.im*(w).c22.re  \
              +(v).c31.re*(w).c32.im-(v).c31.im*(w).c32.re; \
   (u).c13.re+=(v).c11.re*(w).c13.re+(v).c11.im*(w).c13.im  \
              +(v).c21.re*(w).c23.re+(v).c21.im*(w).c23.im  \
              +(v).c31.re*(w).c33.re+(v).c31.im*(w).c33.im; \
   (u).c13.im+=(v).c11.re*(w).c13.im-(v).c11.im*(w).c13.re  \
              +(v).c21.re*(w).c23.im-(v).c21.im*(w).c23.re  \
              +(v).c31.re*(w).c33.im-(v).c31.im*(w).c33.re; \
   (u).c21.re+=(v).c12.re*(w).c11.re+(v).c12.im*(w).c11.im  \
              +(v).c22.re*(w).c21.re+(v).c22.im*(w).c21.im  \
              +(v).c32.re*(w).c31.re+(v).c32.im*(w).c31.im; \
   (u).c21.im+=(v).c12.re*(w).c11.im-(v).c12.im*(w).c11.re  \
              +(v).c22.re*(w).c21.im-(v).c22.im*(w).c21.re  \
              +(v).c32.re*(w).c31.im-(v).c32.im*(w).c31.re; \
   (u).c22.re+=(v).c12.re*(w).c12.re+(v).c12.im*(w).c12.im  \
              +(v).c22.re*(w).c22.re+(v).c22.im*(w).c22.im  \
              +(v).c32.re*(w).c32.re+(v).c32.im*(w).c32.im; \
   (u).c22.im+=(v).c12.re*(w).c12.im-(v).c12.im*(w).c12.re  \
              +(v).c22.re*(w).c22.im-(v).c22.im*(w).c22.re  \
              +(v).c32.re*(w).c32.im-(v).c32.im*(w).c32.re; \
   (u).c23.re+=(v).c12.re*(w).c13.re+(v).c12.im*(w).c13.im  \
              +(v).c22.re*(w).c23.re+(v).c22.im*(w).c23.im  \
              +(v).c32.re*(w).c33.re+(v).c32.im*(w).c33.im; \
   (u).c23.im+=(v).c12.re*(w).c13.im-(v).c12.im*(w).c13.re  \
              +(v).c22.re*(w).c23.im-(v).c22.im*(w).c23.re  \
              +(v).c32.re*(w).c33.im-(v).c32.im*(w).c33.re; \
   (u).c31.re+=(v).c13.re*(w).c11.re+(v).c13.im*(w).c11.im  \
              +(v).c23.re*(w).c21.re+(v).c23.im*(w).c21.im  \
              +(v).c33.re*(w).c31.re+(v).c33.im*(w).c31.im; \
   (u).c31.im+=(v).c13.re*(w).c11.im-(v).c13.im*(w).c11.re  \
              +(v).c23.re*(w).c21.im-(v).c23.im*(w).c21.re  \
              +(v).c33.re*(w).c31.im-(v).c33.im*(w).c31.re; \
   (u).c32.re+=(v).c13.re*(w).c12.re+(v).c13.im*(w).c12.im  \
              +(v).c23.re*(w).c22.re+(v).c23.im*(w).c22.im  \
              +(v).c33.re*(w).c32.re+(v).c33.im*(w).c32.im; \
   (u).c32.im+=(v).c13.re*(w).c12.im-(v).c13.im*(w).c12.re  \
              +(v).c23.re*(w).c22.im-(v).c23.im*(w).c22.re  \
              +(v).c33.re*(w).c32.im-(v).c33.im*(w).c32.re; \
   (u).c33.re+=(v).c13.re*(w).c13.re+(v).c13.im*(w).c13.im  \
              +(v).c23.re*(w).c23.re+(v).c23.im*(w).c23.im  \
              +(v).c33.re*(w).c33.re+(v).c33.im*(w).c33.im; \
   (u).c33.im+=(v).c13.re*(w).c13.im-(v).c13.im*(w).c13.re  \
              +(v).c23.re*(w).c23.im-(v).c23.im*(w).c23.re  \
              +(v).c33.re*(w).c33.im-(v).c33.im*(w).c33.re;

#endif

/* M. Hasenbusch x=Re Tr (v * w^{\dag}) */

#define _trace_su3_times_su3d(x,v,w) \
   x = (v).c11.re*(w).c11.re+(v).c11.im*(w).c11.im  \
      +(v).c12.re*(w).c12.re+(v).c12.im*(w).c12.im  \
      +(v).c13.re*(w).c13.re+(v).c13.im*(w).c13.im  \
      +(v).c21.re*(w).c21.re+(v).c21.im*(w).c21.im  \
      +(v).c22.re*(w).c22.re+(v).c22.im*(w).c22.im  \
      +(v).c23.re*(w).c23.re+(v).c23.im*(w).c23.im  \
      +(v).c31.re*(w).c31.re+(v).c31.im*(w).c31.im  \
      +(v).c32.re*(w).c32.re+(v).c32.im*(w).c32.im  \
      +(v).c33.re*(w).c33.re+(v).c33.im*(w).c33.im;

/* M. Hasenbusch t =u tensor v^dag */
#define _vector_tensor_vector(t,u,v) \
   (t).c11.re=(u).c1.re*(v).c1.re+(u).c1.im*(v).c1.im; \
   (t).c11.im=(u).c1.im*(v).c1.re-(u).c1.re*(v).c1.im; \
   (t).c12.re=(u).c1.re*(v).c2.re+(u).c1.im*(v).c2.im; \
   (t).c12.im=(u).c1.im*(v).c2.re-(u).c1.re*(v).c2.im; \
   (t).c13.re=(u).c1.re*(v).c3.re+(u).c1.im*(v).c3.im; \
   (t).c13.im=(u).c1.im*(v).c3.re-(u).c1.re*(v).c3.im; \
   (t).c21.re=(u).c2.re*(v).c1.re+(u).c2.im*(v).c1.im; \
   (t).c21.im=(u).c2.im*(v).c1.re-(u).c2.re*(v).c1.im; \
   (t).c22.re=(u).c2.re*(v).c2.re+(u).c2.im*(v).c2.im; \
   (t).c22.im=(u).c2.im*(v).c2.re-(u).c2.re*(v).c2.im; \
   (t).c23.re=(u).c2.re*(v).c3.re+(u).c2.im*(v).c3.im; \
   (t).c23.im=(u).c2.im*(v).c3.re-(u).c2.re*(v).c3.im; \
   (t).c31.re=(u).c3.re*(v).c1.re+(u).c3.im*(v).c1.im; \
   (t).c31.im=(u).c3.im*(v).c1.re-(u).c3.re*(v).c1.im; \
   (t).c32.re=(u).c3.re*(v).c2.re+(u).c3.im*(v).c2.im; \
   (t).c32.im=(u).c3.im*(v).c2.re-(u).c3.re*(v).c2.im; \
   (t).c33.re=(u).c3.re*(v).c3.re+(u).c3.im*(v).c3.im; \
   (t).c33.im=(u).c3.im*(v).c3.re-(u).c3.re*(v).c3.im; 

#endif

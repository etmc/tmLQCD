/* $Id$ */
/*******************************************************************************
*
* File linalg.c
*
* Linear combination and norm of spinor fields
*
* The externally accessible function are
*
*   void assign_add_mul(double c,int k,int l)
*     Adds c*psi[k] to psi[l]
*
*   double square_norm(int k)
*     Returns the square norm of psi[k]
*
* Based on
* Version: 2.0
* Author: Martin Luescher <luscher@mail.desy.de>
* Date: 21.03.2001
* Present version: M.Hasenbusch Thu Oct 25 10:58:45 METDST 2001
* Extended and adapted for even-odd precondtioning
*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#ifdef MPI
#include <mpi.h>
#endif
#include "sse.h"
#include "su3.h"
#include "su3adj.h"
#include "global.h"
#include "linalg_eo.h"

/* j output k input , l input */
void diff(const int j, const int k, const int l, const int N)
{
   int ix;
   spinor *q,*r,*s;


/* Change due to even-odd preconditioning : VOLUME   to VOLUME/2 */   
   for (ix = 0; ix < N; ix++)
   {
      q=&spinor_field[j][ix];      
      r=&spinor_field[k][ix];      
      s=&spinor_field[l][ix];

      (*q).c1.c1.re=(*r).c1.c1.re-(*s).c1.c1.re;
      (*q).c1.c1.im=(*r).c1.c1.im-(*s).c1.c1.im;
      (*q).c1.c2.re=(*r).c1.c2.re-(*s).c1.c2.re;
      (*q).c1.c2.im=(*r).c1.c2.im-(*s).c1.c2.im;
      (*q).c1.c3.re=(*r).c1.c3.re-(*s).c1.c3.re;
      (*q).c1.c3.im=(*r).c1.c3.im-(*s).c1.c3.im;

      (*q).c2.c1.re=(*r).c2.c1.re-(*s).c2.c1.re;
      (*q).c2.c1.im=(*r).c2.c1.im-(*s).c2.c1.im;
      (*q).c2.c2.re=(*r).c2.c2.re-(*s).c2.c2.re;
      (*q).c2.c2.im=(*r).c2.c2.im-(*s).c2.c2.im;
      (*q).c2.c3.re=(*r).c2.c3.re-(*s).c2.c3.re;
      (*q).c2.c3.im=(*r).c2.c3.im-(*s).c2.c3.im;         

      (*q).c3.c1.re=(*r).c3.c1.re-(*s).c3.c1.re;
      (*q).c3.c1.im=(*r).c3.c1.im-(*s).c3.c1.im;
      (*q).c3.c2.re=(*r).c3.c2.re-(*s).c3.c2.re;
      (*q).c3.c2.im=(*r).c3.c2.im-(*s).c3.c2.im;
      (*q).c3.c3.re=(*r).c3.c3.re-(*s).c3.c3.re;
      (*q).c3.c3.im=(*r).c3.c3.im-(*s).c3.c3.im;         
      
      (*q).c4.c1.re=(*r).c4.c1.re-(*s).c4.c1.re;
      (*q).c4.c1.im=(*r).c4.c1.im-(*s).c4.c1.im;
      (*q).c4.c2.re=(*r).c4.c2.re-(*s).c4.c2.re;
      (*q).c4.c2.im=(*r).c4.c2.im-(*s).c4.c2.im;
      (*q).c4.c3.re=(*r).c4.c3.re-(*s).c4.c3.re;
      (*q).c4.c3.im=(*r).c4.c3.im-(*s).c4.c3.im;
   }
}
#if defined SSE2
/*  k input, l output */
void assign_add_mul(const int l, const double c, const int k, const int N) {
  int ix;
  su3_vector *s,*r;
  __asm__ __volatile__ ("movsd %0, %%xmm7 \n\t"
			"unpcklpd %%xmm7, %%xmm7"
			:
			:
			"m" (c));
  s=&spinor_field[l][0].c1;
  r=&spinor_field[k][0].c1;
  for (ix = 0;ix < 4*N; ix++) {
    _sse_load_up(*r);
    __asm__ __volatile__ ("mulpd %%xmm7, %%xmm3 \n\t"
			  "mulpd %%xmm7, %%xmm4 \n\t"
			  "mulpd %%xmm7, %%xmm5"
			  :
			  :);
    _sse_load(*s);
    __asm__ __volatile__ ("addpd %%xmm3, %%xmm0 \n\t"
			  "addpd %%xmm4, %%xmm1 \n\t"
			  "addpd %%xmm5, %%xmm2"
			  :
			  :);
    _sse_store(*s);
    s++; r++;
  }
}

#else
/*  k input, l output */
void assign_add_mul(const int l, const double c, const int k, const int N) {
  int ix;
  static double fact;
  spinor *r,*s;

  fact=c;

  /* Change due to even-odd preconditioning : VOLUME   to VOLUME/2 */   
  for (ix = 0; ix < N; ix++) {
    r=&spinor_field[l][ix];      
    s=&spinor_field[k][ix];
    
    (*r).c1.c1.re+=fact*(*s).c1.c1.re;
    (*r).c1.c1.im+=fact*(*s).c1.c1.im;
    (*r).c1.c2.re+=fact*(*s).c1.c2.re;
    (*r).c1.c2.im+=fact*(*s).c1.c2.im;
    (*r).c1.c3.re+=fact*(*s).c1.c3.re;
    (*r).c1.c3.im+=fact*(*s).c1.c3.im;
    
    (*r).c2.c1.re+=fact*(*s).c2.c1.re;
    (*r).c2.c1.im+=fact*(*s).c2.c1.im;
    (*r).c2.c2.re+=fact*(*s).c2.c2.re;
    (*r).c2.c2.im+=fact*(*s).c2.c2.im;
    (*r).c2.c3.re+=fact*(*s).c2.c3.re;
    (*r).c2.c3.im+=fact*(*s).c2.c3.im;         
    
    (*r).c3.c1.re+=fact*(*s).c3.c1.re;
    (*r).c3.c1.im+=fact*(*s).c3.c1.im;
    (*r).c3.c2.re+=fact*(*s).c3.c2.re;
    (*r).c3.c2.im+=fact*(*s).c3.c2.im;
    (*r).c3.c3.re+=fact*(*s).c3.c3.re;
    (*r).c3.c3.im+=fact*(*s).c3.c3.im;         
    
    (*r).c4.c1.re+=fact*(*s).c4.c1.re;
    (*r).c4.c1.im+=fact*(*s).c4.c1.im;
    (*r).c4.c2.re+=fact*(*s).c4.c2.re;
    (*r).c4.c2.im+=fact*(*s).c4.c2.im;
    (*r).c4.c3.re+=fact*(*s).c4.c3.re;
    (*r).c4.c3.im+=fact*(*s).c4.c3.im;
  }
}
#endif
#if defined SSE2
void assign_add_mul_add_mul(const int l, const double c1, const int k, const double c2, const int j, const int N) {

  int ix;
  su3_vector *s,*r,*t;
  r=&spinor_field[l][0].c1;
  s=&spinor_field[k][0].c1;
  t=&spinor_field[j][0].c1;
  __asm__ __volatile__ ("movsd %0, %%xmm6 \n\t"
			"unpcklpd %%xmm6, %%xmm6"
			:
			:
			"m" (c1));
  __asm__ __volatile__ ("movsd %0, %%xmm7 \n\t"
			"unpcklpd %%xmm7, %%xmm7"
			:
			:
			"m" (c2));
  
  for (ix = 0; ix < 4*N; ix++) {
    _sse_load_up(*s);
    __asm__ __volatile__ ("mulpd %%xmm6, %%xmm3 \n\t"
			  "mulpd %%xmm6, %%xmm4 \n\t"
			  "mulpd %%xmm6, %%xmm5"
			  :
			  :);
    _sse_load(*r);
    __asm__ __volatile__ ("addpd %%xmm3, %%xmm0 \n\t"
			  "addpd %%xmm4, %%xmm1 \n\t"
			  "addpd %%xmm5, %%xmm2"
			  :
			  :);
    _sse_load_up(*t);
    __asm__ __volatile__ ("mulpd %%xmm7, %%xmm3 \n\t"
			  "mulpd %%xmm7, %%xmm4 \n\t"
			  "mulpd %%xmm7, %%xmm5"
			  :
			  :);
    __asm__ __volatile__ ("addpd %%xmm3, %%xmm0 \n\t"
			  "addpd %%xmm4, %%xmm1 \n\t"
			  "addpd %%xmm5, %%xmm2"
			  :
			  :);
    _sse_store(*r);
    r++; s++; t++;
  }
}
#else
/* j, k input, l output */
void assign_add_mul_add_mul(const int l, const double c1, const int k, const double c2, const int j, const int N) {

  int ix;
   spinor *r,*s,*t;

/* Change due to even-odd preconditioning : VOLUME   to VOLUME/2 */   
   for (ix = 0; ix < N; ix++) {
     r=&spinor_field[l][ix];      
     s=&spinor_field[k][ix];
     t=&spinor_field[j][ix];
     
     (*r).c1.c1.re+=c1*(*s).c1.c1.re+c2*(*t).c1.c1.re;
     (*r).c1.c1.im+=c1*(*s).c1.c1.im+c2*(*t).c1.c1.im;
     (*r).c1.c2.re+=c1*(*s).c1.c2.re+c2*(*t).c1.c2.re;
     (*r).c1.c2.im+=c1*(*s).c1.c2.im+c2*(*t).c1.c2.im;
     (*r).c1.c3.re+=c1*(*s).c1.c3.re+c2*(*t).c1.c3.re;
     (*r).c1.c3.im+=c1*(*s).c1.c3.im+c2*(*t).c1.c3.im;
     
     (*r).c2.c1.re+=c1*(*s).c2.c1.re+c2*(*t).c2.c1.re;
     (*r).c2.c1.im+=c1*(*s).c2.c1.im+c2*(*t).c2.c1.im;
     (*r).c2.c2.re+=c1*(*s).c2.c2.re+c2*(*t).c2.c2.re;
     (*r).c2.c2.im+=c1*(*s).c2.c2.im+c2*(*t).c2.c2.im;
     (*r).c2.c3.re+=c1*(*s).c2.c3.re+c2*(*t).c2.c3.re;
     (*r).c2.c3.im+=c1*(*s).c2.c3.im+c2*(*t).c2.c3.im;         
     
     (*r).c3.c1.re+=c1*(*s).c3.c1.re+c2*(*t).c3.c1.re;
     (*r).c3.c1.im+=c1*(*s).c3.c1.im+c2*(*t).c3.c1.im;
     (*r).c3.c2.re+=c1*(*s).c3.c2.re+c2*(*t).c3.c2.re;
     (*r).c3.c2.im+=c1*(*s).c3.c2.im+c2*(*t).c3.c2.im;
     (*r).c3.c3.re+=c1*(*s).c3.c3.re+c2*(*t).c3.c3.re;
     (*r).c3.c3.im+=c1*(*s).c3.c3.im+c2*(*t).c3.c3.im;         
     
     (*r).c4.c1.re+=c1*(*s).c4.c1.re+c2*(*t).c4.c1.re;
     (*r).c4.c1.im+=c1*(*s).c4.c1.im+c2*(*t).c4.c1.im;
     (*r).c4.c2.re+=c1*(*s).c4.c2.re+c2*(*t).c4.c2.re;
     (*r).c4.c2.im+=c1*(*s).c4.c2.im+c2*(*t).c4.c2.im;
     (*r).c4.c3.re+=c1*(*s).c4.c3.re+c2*(*t).c4.c3.re;
     (*r).c4.c3.im+=c1*(*s).c4.c3.im+c2*(*t).c4.c3.im;
   }
}
#endif

void assign_mul_bra_add_mul_ket_add(const int l,const double c1,const int k,
				    const double c2, const int j, const int N) {
   int ix;
   spinor *r,*s,*t;

/* Change due to even-odd preconditioning : VOLUME   to VOLUME/2 */   
   for (ix = 0; ix < N; ix++) {
     r=&spinor_field[l][ix];      
     s=&spinor_field[k][ix];
     t=&spinor_field[j][ix];
     (*r).c1.c1.re=c2*( (*r).c1.c1.re+c1*(*s).c1.c1.re )+(*t).c1.c1.re;
     (*r).c1.c1.im=c2*( (*r).c1.c1.im+c1*(*s).c1.c1.im )+(*t).c1.c1.im;
     (*r).c1.c2.re=c2*( (*r).c1.c2.re+c1*(*s).c1.c2.re )+(*t).c1.c2.re;
     (*r).c1.c2.im=c2*( (*r).c1.c2.im+c1*(*s).c1.c2.im )+(*t).c1.c2.im;
     (*r).c1.c3.re=c2*( (*r).c1.c3.re+c1*(*s).c1.c3.re )+(*t).c1.c3.re;
     (*r).c1.c3.im=c2*( (*r).c1.c3.im+c1*(*s).c1.c3.im )+(*t).c1.c3.im;
     
     (*r).c2.c1.re=c2*( (*r).c2.c1.re+c1*(*s).c2.c1.re )+(*t).c2.c1.re;
     (*r).c2.c1.im=c2*( (*r).c2.c1.im+c1*(*s).c2.c1.im )+(*t).c2.c1.im;
     (*r).c2.c2.re=c2*( (*r).c2.c2.re+c1*(*s).c2.c2.re )+(*t).c2.c2.re;
     (*r).c2.c2.im=c2*( (*r).c2.c2.im+c1*(*s).c2.c2.im )+(*t).c2.c2.im;
     (*r).c2.c3.re=c2*( (*r).c2.c3.re+c1*(*s).c2.c3.re )+(*t).c2.c3.re;
     (*r).c2.c3.im=c2*( (*r).c2.c3.im+c1*(*s).c2.c3.im )+(*t).c2.c3.im;
     
     (*r).c3.c1.re=c2*( (*r).c3.c1.re+c1*(*s).c3.c1.re )+(*t).c3.c1.re;
     (*r).c3.c1.im=c2*( (*r).c3.c1.im+c1*(*s).c3.c1.im )+(*t).c3.c1.im;
     (*r).c3.c2.re=c2*( (*r).c3.c2.re+c1*(*s).c3.c2.re )+(*t).c3.c2.re;
     (*r).c3.c2.im=c2*( (*r).c3.c2.im+c1*(*s).c3.c2.im )+(*t).c3.c2.im;
     (*r).c3.c3.re=c2*( (*r).c3.c3.re+c1*(*s).c3.c3.re )+(*t).c3.c3.re;
     (*r).c3.c3.im=c2*( (*r).c3.c3.im+c1*(*s).c3.c3.im )+(*t).c3.c3.im;
     
     (*r).c4.c1.re=c2*( (*r).c4.c1.re+c1*(*s).c4.c1.re )+(*t).c4.c1.re;
     (*r).c4.c1.im=c2*( (*r).c4.c1.im+c1*(*s).c4.c1.im )+(*t).c4.c1.im;
     (*r).c4.c2.re=c2*( (*r).c4.c2.re+c1*(*s).c4.c2.re )+(*t).c4.c2.re;
     (*r).c4.c2.im=c2*( (*r).c4.c2.im+c1*(*s).c4.c2.im )+(*t).c4.c2.im;
     (*r).c4.c3.re=c2*( (*r).c4.c3.re+c1*(*s).c4.c3.re )+(*t).c4.c3.re;
     (*r).c4.c3.im=c2*( (*r).c4.c3.im+c1*(*s).c4.c3.im )+(*t).c4.c3.im;
   }
}

/* j, k input, l output */
void deri_linalg(const int l, const double c1, const int k, const double c2, const int j, const int N) {
  int ix;
  spinor *r,*s,*t;
  
  /* Change due to even-odd preconditioning : VOLUME   to VOLUME/2 */   
  for (ix = 0; ix < N; ix++) {
    r=&spinor_field[l][ix];
    s=&spinor_field[k][ix];
    t=&spinor_field[j][ix];
    
    (*r).c1.c1.re=c1*(*s).c1.c1.re+c2*(*t).c1.c1.re;
    (*r).c1.c1.im=c1*(*s).c1.c1.im+c2*(*t).c1.c1.im;
    (*r).c1.c2.re=c1*(*s).c1.c2.re+c2*(*t).c1.c2.re;
    (*r).c1.c2.im=c1*(*s).c1.c2.im+c2*(*t).c1.c2.im;
    (*r).c1.c3.re=c1*(*s).c1.c3.re+c2*(*t).c1.c3.re;
    (*r).c1.c3.im=c1*(*s).c1.c3.im+c2*(*t).c1.c3.im;
    
    (*r).c2.c1.re=c1*(*s).c2.c1.re+c2*(*t).c2.c1.re;
    (*r).c2.c1.im=c1*(*s).c2.c1.im+c2*(*t).c2.c1.im;
    (*r).c2.c2.re=c1*(*s).c2.c2.re+c2*(*t).c2.c2.re;
    (*r).c2.c2.im=c1*(*s).c2.c2.im+c2*(*t).c2.c2.im;
    (*r).c2.c3.re=c1*(*s).c2.c3.re+c2*(*t).c2.c3.re;
    (*r).c2.c3.im=c1*(*s).c2.c3.im+c2*(*t).c2.c3.im;
    
    (*r).c3.c1.re=c1*(*s).c3.c1.re-c2*(*t).c3.c1.re;
    (*r).c3.c1.im=c1*(*s).c3.c1.im-c2*(*t).c3.c1.im;
    (*r).c3.c2.re=c1*(*s).c3.c2.re-c2*(*t).c3.c2.re;
    (*r).c3.c2.im=c1*(*s).c3.c2.im-c2*(*t).c3.c2.im;
    (*r).c3.c3.re=c1*(*s).c3.c3.re-c2*(*t).c3.c3.re;
    (*r).c3.c3.im=c1*(*s).c3.c3.im-c2*(*t).c3.c3.im;
    
    (*r).c4.c1.re=c1*(*s).c4.c1.re-c2*(*t).c4.c1.re;
    (*r).c4.c1.im=c1*(*s).c4.c1.im-c2*(*t).c4.c1.im;
    (*r).c4.c2.re=c1*(*s).c4.c2.re-c2*(*t).c4.c2.re;
    (*r).c4.c2.im=c1*(*s).c4.c2.im-c2*(*t).c4.c2.im;
    (*r).c4.c3.re=c1*(*s).c4.c3.re-c2*(*t).c4.c3.re;
    (*r).c4.c3.im=c1*(*s).c4.c3.im-c2*(*t).c4.c3.im;
  }
}

/* k input, l output */
void assign(const int l, const int k, const int N) {

  int ix;
  spinor *r,*s;
  
  /* Change due to even-odd preconditioning : VOLUME   to VOLUME/2 */   
  for (ix = 0; ix < N; ix++) {
    r=&spinor_field[l][ix];      
    s=&spinor_field[k][ix];
    
    (*r).c1.c1.re=(*s).c1.c1.re;
    (*r).c1.c1.im=(*s).c1.c1.im;
    (*r).c1.c2.re=(*s).c1.c2.re;
    (*r).c1.c2.im=(*s).c1.c2.im;
    (*r).c1.c3.re=(*s).c1.c3.re;
    (*r).c1.c3.im=(*s).c1.c3.im;
    
    (*r).c2.c1.re=(*s).c2.c1.re;
    (*r).c2.c1.im=(*s).c2.c1.im;
    (*r).c2.c2.re=(*s).c2.c2.re;
    (*r).c2.c2.im=(*s).c2.c2.im;
    (*r).c2.c3.re=(*s).c2.c3.re;
    (*r).c2.c3.im=(*s).c2.c3.im;         
    
    (*r).c3.c1.re=(*s).c3.c1.re;
    (*r).c3.c1.im=(*s).c3.c1.im;
    (*r).c3.c2.re=(*s).c3.c2.re;
    (*r).c3.c2.im=(*s).c3.c2.im;
    (*r).c3.c3.re=(*s).c3.c3.re;
    (*r).c3.c3.im=(*s).c3.c3.im;         
    
    (*r).c4.c1.re=(*s).c4.c1.re;
    (*r).c4.c1.im=(*s).c4.c1.im;
    (*r).c4.c2.re=(*s).c4.c2.re;
    (*r).c4.c2.im=(*s).c4.c2.im;
    (*r).c4.c3.re=(*s).c4.c3.re;
    (*r).c4.c3.im=(*s).c4.c3.im;
  }
}

#if defined SSE2
/* k input , l output*/
void assign_mul_add_r(const int l, const double c, const int k, const int N) {

  int ix;
  su3_vector *s,*r;
  __asm__ __volatile__ ("movsd %0, %%xmm7 \n\t"
			"unpcklpd %%xmm7, %%xmm7"
			:
			:
			"m" (c));
  s=&spinor_field[l][0].c1;
  r=&spinor_field[k][0].c1;
  for (ix=0;ix<4*N;ix++) {
    _sse_load(*s);
    __asm__ __volatile__ ("mulpd %%xmm7, %%xmm0 \n\t"
			  "mulpd %%xmm7, %%xmm1 \n\t"
			  "mulpd %%xmm7, %%xmm2"
			  :
			  :);
    _sse_load_up(*r);
    __asm__ __volatile__ ("addpd %%xmm3, %%xmm0 \n\t"
			  "addpd %%xmm4, %%xmm1 \n\t"
			  "addpd %%xmm5, %%xmm2"
			  :
			  :);
    _sse_store(*s);
    s++; r++;
  }
}
#else
/* k input , l output*/
void assign_mul_add_r(const int l, const double c, const int k, const int N) {

  int ix;
  static double fact;
  spinor *r,*s;
  
  fact=c;
  
  /* Change due to even-odd preconditioning : VOLUME   to VOLUME/2 */   
  for (ix = 0; ix < N; ix++) {
    r=&spinor_field[l][ix];      
    s=&spinor_field[k][ix];
    
    (*r).c1.c1.re=fact*(*r).c1.c1.re+(*s).c1.c1.re;
    (*r).c1.c1.im=fact*(*r).c1.c1.im+(*s).c1.c1.im;
    (*r).c1.c2.re=fact*(*r).c1.c2.re+(*s).c1.c2.re;
    (*r).c1.c2.im=fact*(*r).c1.c2.im+(*s).c1.c2.im;
    (*r).c1.c3.re=fact*(*r).c1.c3.re+(*s).c1.c3.re;
    (*r).c1.c3.im=fact*(*r).c1.c3.im+(*s).c1.c3.im;
    
    (*r).c2.c1.re=fact*(*r).c2.c1.re+(*s).c2.c1.re;
    (*r).c2.c1.im=fact*(*r).c2.c1.im+(*s).c2.c1.im;
    (*r).c2.c2.re=fact*(*r).c2.c2.re+(*s).c2.c2.re;
    (*r).c2.c2.im=fact*(*r).c2.c2.im+(*s).c2.c2.im;
    (*r).c2.c3.re=fact*(*r).c2.c3.re+(*s).c2.c3.re;
    (*r).c2.c3.im=fact*(*r).c2.c3.im+(*s).c2.c3.im;         
    
    (*r).c3.c1.re=fact*(*r).c3.c1.re+(*s).c3.c1.re;
    (*r).c3.c1.im=fact*(*r).c3.c1.im+(*s).c3.c1.im;
    (*r).c3.c2.re=fact*(*r).c3.c2.re+(*s).c3.c2.re;
    (*r).c3.c2.im=fact*(*r).c3.c2.im+(*s).c3.c2.im;
    (*r).c3.c3.re=fact*(*r).c3.c3.re+(*s).c3.c3.re;
    (*r).c3.c3.im=fact*(*r).c3.c3.im+(*s).c3.c3.im;
    
    (*r).c4.c1.re=fact*(*r).c4.c1.re+(*s).c4.c1.re;
    (*r).c4.c1.im=fact*(*r).c4.c1.im+(*s).c4.c1.im;
    (*r).c4.c2.re=fact*(*r).c4.c2.re+(*s).c4.c2.re;
    (*r).c4.c2.im=fact*(*r).c4.c2.im+(*s).c4.c2.im;
    (*r).c4.c3.re=fact*(*r).c4.c3.re+(*s).c4.c3.re;
    (*r).c4.c3.im=fact*(*r).c4.c3.im+(*s).c4.c3.im;
  }
}
#endif

double square_norm(const int k, const int N) {
  int ix;
  static double ks,kc,ds,tr,ts,tt;
  spinor *s;
  
  ks=0.0;
  kc=0.0;
  
  /* Change due to even-odd preconditioning : VOLUME   to VOLUME/2 */   
  for (ix = 0; ix < N; ix++) {
    s=&spinor_field[k][ix];
    
    ds=(*s).c1.c1.re*(*s).c1.c1.re+(*s).c1.c1.im*(*s).c1.c1.im+
      (*s).c1.c2.re*(*s).c1.c2.re+(*s).c1.c2.im*(*s).c1.c2.im+
      (*s).c1.c3.re*(*s).c1.c3.re+(*s).c1.c3.im*(*s).c1.c3.im+
      (*s).c2.c1.re*(*s).c2.c1.re+(*s).c2.c1.im*(*s).c2.c1.im+
      (*s).c2.c2.re*(*s).c2.c2.re+(*s).c2.c2.im*(*s).c2.c2.im+
      (*s).c2.c3.re*(*s).c2.c3.re+(*s).c2.c3.im*(*s).c2.c3.im+
      (*s).c3.c1.re*(*s).c3.c1.re+(*s).c3.c1.im*(*s).c3.c1.im+
      (*s).c3.c2.re*(*s).c3.c2.re+(*s).c3.c2.im*(*s).c3.c2.im+
      (*s).c3.c3.re*(*s).c3.c3.re+(*s).c3.c3.im*(*s).c3.c3.im+
      (*s).c4.c1.re*(*s).c4.c1.re+(*s).c4.c1.im*(*s).c4.c1.im+
      (*s).c4.c2.re*(*s).c4.c2.re+(*s).c4.c2.im*(*s).c4.c2.im+
      (*s).c4.c3.re*(*s).c4.c3.re+(*s).c4.c3.im*(*s).c4.c3.im;
    
    tr=ds+kc;
    ts=tr+ks;
    tt=ts-ks;
    ks=ts;
    kc=tr-tt;
  }
  kc=ks+kc;
#ifdef MPI
  MPI_Allreduce(&kc, &ks, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return ks;
#else
  return kc;
#endif
}

double scalar_prod_r(const int k, const int l, const int N) {
  int ix;
  static double ks,kc,ds,tr,ts,tt;
  spinor *s, *r;
  
  ks=0.0;
  kc=0.0;
  
  /* Change due to even-odd preconditioning : VOLUME   to VOLUME/2 */   
  for (ix = 0; ix < N; ix++) {
    s=&spinor_field[k][ix];
    r=&spinor_field[l][ix];
    
    ds=(*r).c1.c1.re*(*s).c1.c1.re+(*r).c1.c1.im*(*s).c1.c1.im+
      (*r).c1.c2.re*(*s).c1.c2.re+(*r).c1.c2.im*(*s).c1.c2.im+
      (*r).c1.c3.re*(*s).c1.c3.re+(*r).c1.c3.im*(*s).c1.c3.im+
      (*r).c2.c1.re*(*s).c2.c1.re+(*r).c2.c1.im*(*s).c2.c1.im+
      (*r).c2.c2.re*(*s).c2.c2.re+(*r).c2.c2.im*(*s).c2.c2.im+
      (*r).c2.c3.re*(*s).c2.c3.re+(*r).c2.c3.im*(*s).c2.c3.im+
      (*r).c3.c1.re*(*s).c3.c1.re+(*r).c3.c1.im*(*s).c3.c1.im+
      (*r).c3.c2.re*(*s).c3.c2.re+(*r).c3.c2.im*(*s).c3.c2.im+
      (*r).c3.c3.re*(*s).c3.c3.re+(*r).c3.c3.im*(*s).c3.c3.im+
      (*r).c4.c1.re*(*s).c4.c1.re+(*r).c4.c1.im*(*s).c4.c1.im+
      (*r).c4.c2.re*(*s).c4.c2.re+(*r).c4.c2.im*(*s).c4.c2.im+
      (*r).c4.c3.re*(*s).c4.c3.re+(*r).c4.c3.im*(*s).c4.c3.im;
    
    tr=ds+kc;
    ts=tr+ks;
    tt=ts-ks;
    ks=ts;
    kc=tr-tt;
   }
  kc=ks+kc;
#ifdef MPI
   MPI_Allreduce(&kc, &ks, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   return ks;
#else
   return kc;
#endif
}

void square_and_prod_r(double *x1, double *x2, const int k, const int l, const int N) {
  int ix;
  static double ks,kc,ds,tr,ts,tt;
  static double xks,xkc,xds,xtr,xts,xtt;
  spinor *s, *r;
  
  ks=0.0;
  kc=0.0;
  
  xks=0.0;
  xkc=0.0;
  
  /* Change due to even-odd preconditioning : VOLUME   to VOLUME/2 */   
  for (ix = 0; ix < N; ix++) {
    s=&spinor_field[k][ix];
    r=&spinor_field[l][ix];
    
    ds=(*r).c1.c1.re*(*s).c1.c1.re+(*r).c1.c1.im*(*s).c1.c1.im+
      (*r).c1.c2.re*(*s).c1.c2.re+(*r).c1.c2.im*(*s).c1.c2.im+
      (*r).c1.c3.re*(*s).c1.c3.re+(*r).c1.c3.im*(*s).c1.c3.im+
      (*r).c2.c1.re*(*s).c2.c1.re+(*r).c2.c1.im*(*s).c2.c1.im+
      (*r).c2.c2.re*(*s).c2.c2.re+(*r).c2.c2.im*(*s).c2.c2.im+
      (*r).c2.c3.re*(*s).c2.c3.re+(*r).c2.c3.im*(*s).c2.c3.im+
      (*r).c3.c1.re*(*s).c3.c1.re+(*r).c3.c1.im*(*s).c3.c1.im+
      (*r).c3.c2.re*(*s).c3.c2.re+(*r).c3.c2.im*(*s).c3.c2.im+
      (*r).c3.c3.re*(*s).c3.c3.re+(*r).c3.c3.im*(*s).c3.c3.im+
      (*r).c4.c1.re*(*s).c4.c1.re+(*r).c4.c1.im*(*s).c4.c1.im+
      (*r).c4.c2.re*(*s).c4.c2.re+(*r).c4.c2.im*(*s).c4.c2.im+
      (*r).c4.c3.re*(*s).c4.c3.re+(*r).c4.c3.im*(*s).c4.c3.im;
    
    xds=(*s).c1.c1.re*(*s).c1.c1.re+(*s).c1.c1.im*(*s).c1.c1.im+
      (*s).c1.c2.re*(*s).c1.c2.re+(*s).c1.c2.im*(*s).c1.c2.im+
      (*s).c1.c3.re*(*s).c1.c3.re+(*s).c1.c3.im*(*s).c1.c3.im+
      (*s).c2.c1.re*(*s).c2.c1.re+(*s).c2.c1.im*(*s).c2.c1.im+
      (*s).c2.c2.re*(*s).c2.c2.re+(*s).c2.c2.im*(*s).c2.c2.im+
      (*s).c2.c3.re*(*s).c2.c3.re+(*s).c2.c3.im*(*s).c2.c3.im+
      (*s).c3.c1.re*(*s).c3.c1.re+(*s).c3.c1.im*(*s).c3.c1.im+
      (*s).c3.c2.re*(*s).c3.c2.re+(*s).c3.c2.im*(*s).c3.c2.im+
      (*s).c3.c3.re*(*s).c3.c3.re+(*s).c3.c3.im*(*s).c3.c3.im+
      (*s).c4.c1.re*(*s).c4.c1.re+(*s).c4.c1.im*(*s).c4.c1.im+
      (*s).c4.c2.re*(*s).c4.c2.re+(*s).c4.c2.im*(*s).c4.c2.im+
      (*s).c4.c3.re*(*s).c4.c3.re+(*s).c4.c3.im*(*s).c4.c3.im;
    
    tr=ds+kc;
    ts=tr+ks;
    tt=ts-ks;
    ks=ts;
    kc=tr-tt;
    
    xtr=xds+xkc;
    xts=xtr+xks;
    xtt=xts-xks;
    xks=xts;
    xkc=xtr-xtt;
  }
  xkc=xks+xkc;
#ifdef MPI
  MPI_Allreduce(&xkc, x1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
  (*x1) =  xkc;
#endif
  kc=ks+kc;
#ifdef MPI
  MPI_Allreduce(&kc, x2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
  (*x2) = kc;
#endif
}

/*  k input, l output */
void assign_mul_bra_add_mul_r(const int l, const double c0, const double c, const int k, const int N)
{
  int ix;
  static double fact0,fact;
  spinor *r,*s;

  fact=c;
  fact0=c0;

/* Change due to even-odd preconditioning : VOLUME   to VOLUME/2 */   
  for (ix = 0; ix < N; ix++) {
    r=&spinor_field[l][ix];
    s=&spinor_field[k][ix];
    
    (*r).c1.c1.re=fact0*((*r).c1.c1.re+fact*(*s).c1.c1.re);
    (*r).c1.c1.im=fact0*((*r).c1.c1.im+fact*(*s).c1.c1.im);
    (*r).c1.c2.re=fact0*((*r).c1.c2.re+fact*(*s).c1.c2.re);
    (*r).c1.c2.im=fact0*((*r).c1.c2.im+fact*(*s).c1.c2.im);
    (*r).c1.c3.re=fact0*((*r).c1.c3.re+fact*(*s).c1.c3.re);
    (*r).c1.c3.im=fact0*((*r).c1.c3.im+fact*(*s).c1.c3.im);

    (*r).c2.c1.re=fact0*((*r).c2.c1.re+fact*(*s).c2.c1.re);
    (*r).c2.c1.im=fact0*((*r).c2.c1.im+fact*(*s).c2.c1.im);
    (*r).c2.c2.re=fact0*((*r).c2.c2.re+fact*(*s).c2.c2.re);
    (*r).c2.c2.im=fact0*((*r).c2.c2.im+fact*(*s).c2.c2.im);
    (*r).c2.c3.re=fact0*((*r).c2.c3.re+fact*(*s).c2.c3.re);
    (*r).c2.c3.im=fact0*((*r).c2.c3.im+fact*(*s).c2.c3.im);

    (*r).c3.c1.re=fact0*((*r).c3.c1.re+fact*(*s).c3.c1.re);
    (*r).c3.c1.im=fact0*((*r).c3.c1.im+fact*(*s).c3.c1.im);
    (*r).c3.c2.re=fact0*((*r).c3.c2.re+fact*(*s).c3.c2.re);
    (*r).c3.c2.im=fact0*((*r).c3.c2.im+fact*(*s).c3.c2.im);
    (*r).c3.c3.re=fact0*((*r).c3.c3.re+fact*(*s).c3.c3.re);
    (*r).c3.c3.im=fact0*((*r).c3.c3.im+fact*(*s).c3.c3.im);

    (*r).c4.c1.re=fact0*((*r).c4.c1.re+fact*(*s).c4.c1.re);
    (*r).c4.c1.im=fact0*((*r).c4.c1.im+fact*(*s).c4.c1.im);
    (*r).c4.c2.re=fact0*((*r).c4.c2.re+fact*(*s).c4.c2.re);
    (*r).c4.c2.im=fact0*((*r).c4.c2.im+fact*(*s).c4.c2.im);
    (*r).c4.c3.re=fact0*((*r).c4.c3.re+fact*(*s).c4.c3.re);
    (*r).c4.c3.im=fact0*((*r).c4.c3.im+fact*(*s).c4.c3.im);
  }
}

/* j output k input , l input */
double diff_and_square_norm(const int j, const int k, const int N) {
  int ix;
  static double ks,kc,ds,tr,ts,tt;
  spinor *q,*r;
  
  ks=0.0;
  kc=0.0;
  
  /* Change due to even-odd preconditioning : VOLUME   to VOLUME/2 */   
  for (ix = 0; ix < N; ix++) {
    q=&spinor_field[j][ix];      
    r=&spinor_field[k][ix];      
    
    (*q).c1.c1.re=(*r).c1.c1.re-(*q).c1.c1.re;
    (*q).c1.c1.im=(*r).c1.c1.im-(*q).c1.c1.im;
    (*q).c1.c2.re=(*r).c1.c2.re-(*q).c1.c2.re;
    (*q).c1.c2.im=(*r).c1.c2.im-(*q).c1.c2.im;
    (*q).c1.c3.re=(*r).c1.c3.re-(*q).c1.c3.re;
    (*q).c1.c3.im=(*r).c1.c3.im-(*q).c1.c3.im;
    
    ds=(*q).c1.c1.re*(*q).c1.c1.re+(*q).c1.c1.im*(*q).c1.c1.im+
      (*q).c1.c2.re*(*q).c1.c2.re+(*q).c1.c2.im*(*q).c1.c2.im+
      (*q).c1.c3.re*(*q).c1.c3.re+(*q).c1.c3.im*(*q).c1.c3.im;
    
    (*q).c2.c1.re=(*r).c2.c1.re-(*q).c2.c1.re;
    (*q).c2.c1.im=(*r).c2.c1.im-(*q).c2.c1.im;
    (*q).c2.c2.re=(*r).c2.c2.re-(*q).c2.c2.re;
    (*q).c2.c2.im=(*r).c2.c2.im-(*q).c2.c2.im;
    (*q).c2.c3.re=(*r).c2.c3.re-(*q).c2.c3.re;
    (*q).c2.c3.im=(*r).c2.c3.im-(*q).c2.c3.im;         
    
    ds+=
      (*q).c2.c1.re*(*q).c2.c1.re+(*q).c2.c1.im*(*q).c2.c1.im+
      (*q).c2.c2.re*(*q).c2.c2.re+(*q).c2.c2.im*(*q).c2.c2.im+
      (*q).c2.c3.re*(*q).c2.c3.re+(*q).c2.c3.im*(*q).c2.c3.im;
    
    (*q).c3.c1.re=(*r).c3.c1.re-(*q).c3.c1.re;
    (*q).c3.c1.im=(*r).c3.c1.im-(*q).c3.c1.im;
    (*q).c3.c2.re=(*r).c3.c2.re-(*q).c3.c2.re;
    (*q).c3.c2.im=(*r).c3.c2.im-(*q).c3.c2.im;
    (*q).c3.c3.re=(*r).c3.c3.re-(*q).c3.c3.re;
    (*q).c3.c3.im=(*r).c3.c3.im-(*q).c3.c3.im;         
    
    ds+=
      (*q).c3.c1.re*(*q).c3.c1.re+(*q).c3.c1.im*(*q).c3.c1.im+
      (*q).c3.c2.re*(*q).c3.c2.re+(*q).c3.c2.im*(*q).c3.c2.im+
      (*q).c3.c3.re*(*q).c3.c3.re+(*q).c3.c3.im*(*q).c3.c3.im;
    
    (*q).c4.c1.re=(*r).c4.c1.re-(*q).c4.c1.re;
    (*q).c4.c1.im=(*r).c4.c1.im-(*q).c4.c1.im;
    (*q).c4.c2.re=(*r).c4.c2.re-(*q).c4.c2.re;
    (*q).c4.c2.im=(*r).c4.c2.im-(*q).c4.c2.im;
    (*q).c4.c3.re=(*r).c4.c3.re-(*q).c4.c3.re;
    (*q).c4.c3.im=(*r).c4.c3.im-(*q).c4.c3.im;
    
    ds+=
      (*q).c4.c1.re*(*q).c4.c1.re+(*q).c4.c1.im*(*q).c4.c1.im+
      (*q).c4.c2.re*(*q).c4.c2.re+(*q).c4.c2.im*(*q).c4.c2.im+
      (*q).c4.c3.re*(*q).c4.c3.re+(*q).c4.c3.im*(*q).c4.c3.im;
    
    tr=ds+kc;
    ts=tr+ks;
    tt=ts-ks;
    ks=ts;
    kc=tr-tt;
  }
  kc=ks+kc;
#ifdef MPI
  MPI_Allreduce(&kc, &ks, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return ks;
#else
  return kc;
#endif

}


void mul_r(const int R, const double c, const int S, const int N){
  int ix;
  spinor *r,*s;
  
  for (ix = 0; ix < N; ix++){
    r = &spinor_field[R][ix];
    s = &spinor_field[S][ix];
    
    (*r).c1.c1.re=c*(*s).c1.c1.re;
    (*r).c1.c1.im=c*(*s).c1.c1.im;
    (*r).c1.c2.re=c*(*s).c1.c2.re;
    (*r).c1.c2.im=c*(*s).c1.c2.im;
    (*r).c1.c3.re=c*(*s).c1.c3.re;
    (*r).c1.c3.im=c*(*s).c1.c3.im;
    
    (*r).c2.c1.re=c*(*s).c2.c1.re;
    (*r).c2.c1.im=c*(*s).c2.c1.im;
    (*r).c2.c2.re=c*(*s).c2.c2.re;
    (*r).c2.c2.im=c*(*s).c2.c2.im;
    (*r).c2.c3.re=c*(*s).c2.c3.re;
    (*r).c2.c3.im=c*(*s).c2.c3.im;
    
    (*r).c3.c1.re=c*(*s).c3.c1.re;
    (*r).c3.c1.im=c*(*s).c3.c1.im;
    (*r).c3.c2.re=c*(*s).c3.c2.re;
    (*r).c3.c2.im=c*(*s).c3.c2.im;
    (*r).c3.c3.re=c*(*s).c3.c3.re;
    (*r).c3.c3.im=c*(*s).c3.c3.im;
    
    (*r).c4.c1.re=c*(*s).c4.c1.re;
    (*r).c4.c1.im=c*(*s).c4.c1.im;
    (*r).c4.c2.re=c*(*s).c4.c2.re;
    (*r).c4.c2.im=c*(*s).c4.c2.im;
    (*r).c4.c3.re=c*(*s).c4.c3.re;
    (*r).c4.c3.im=c*(*s).c4.c3.im;
  }
}


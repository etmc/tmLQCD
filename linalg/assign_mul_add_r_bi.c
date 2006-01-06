/* $Id$ */


/************************************************************************
 *
 *      Adpated routine evaluating the P=c*P+Q where P,Q are bispinors
 *
 * Author: Thomas Chiarappa
 *         Thomas.Chiarappa@mib.infn.it
 *
 ************************************************************************/


#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include "su3.h"
#include "sse.h"
#include "assign_mul_add_r_bi.h"


#if defined SSE2
/* k input , l output*/
void assign_mul_add_r_bi(bispinor * const S, const double c, bispinor * const R, const int N) {

  int ix;
  su3_vector *s,*r;
  /*
  su3_vector *t,*u;
  */

  __asm__ __volatile__ ("movsd %0, %%xmm7 \n\t"
			"unpcklpd %%xmm7, %%xmm7"
			:
			:
			"m" (c));
  
  
  s=(su3_vector *) &S[0].sp_up.s0;
  r=(su3_vector *) &R[0].sp_up.s0;

/*  for (ix=0;ix<4*N;ix++) { */
  for (ix=0;ix<2*4*N;ix++) {
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
void assign_mul_add_r_bi(bispinor * const R, const double c, bispinor * const S, const int N) {

  int ix;
  static double fact;
  spinor *r,*s;
  
  fact=c;
  
  /* Change due to even-odd preconditioning : VOLUME   to VOLUME/2 */   
  for (ix = 0; ix < N; ix++) {

    r = (spinor *) &R[ix].sp_up;
    s = (spinor *) &S[ix].sp_up;
    
    (*r).s0.c0.re = fact*(*r).s0.c0.re + (*s).s0.c0.re;
    (*r).s0.c0.im = fact*(*r).s0.c0.im + (*s).s0.c0.im;
    (*r).s0.c1.re = fact*(*r).s0.c1.re + (*s).s0.c1.re;
    (*r).s0.c1.im = fact*(*r).s0.c1.im + (*s).s0.c1.im;
    (*r).s0.c2.re = fact*(*r).s0.c2.re + (*s).s0.c2.re;
    (*r).s0.c2.im = fact*(*r).s0.c2.im + (*s).s0.c2.im;
    
    (*r).s1.c0.re = fact*(*r).s1.c0.re + (*s).s1.c0.re;
    (*r).s1.c0.im = fact*(*r).s1.c0.im + (*s).s1.c0.im;
    (*r).s1.c1.re = fact*(*r).s1.c1.re + (*s).s1.c1.re;
    (*r).s1.c1.im = fact*(*r).s1.c1.im + (*s).s1.c1.im;
    (*r).s1.c2.re = fact*(*r).s1.c2.re + (*s).s1.c2.re;
    (*r).s1.c2.im = fact*(*r).s1.c2.im + (*s).s1.c2.im;         
    
    (*r).s2.c0.re = fact*(*r).s2.c0.re + (*s).s2.c0.re;
    (*r).s2.c0.im = fact*(*r).s2.c0.im + (*s).s2.c0.im;
    (*r).s2.c1.re = fact*(*r).s2.c1.re + (*s).s2.c1.re;
    (*r).s2.c1.im = fact*(*r).s2.c1.im + (*s).s2.c1.im;
    (*r).s2.c2.re = fact*(*r).s2.c2.re + (*s).s2.c2.re;
    (*r).s2.c2.im = fact*(*r).s2.c2.im + (*s).s2.c2.im;
    
    (*r).s3.c0.re = fact*(*r).s3.c0.re + (*s).s3.c0.re;
    (*r).s3.c0.im = fact*(*r).s3.c0.im + (*s).s3.c0.im;
    (*r).s3.c1.re = fact*(*r).s3.c1.re + (*s).s3.c1.re;
    (*r).s3.c1.im = fact*(*r).s3.c1.im + (*s).s3.c1.im;
    (*r).s3.c2.re = fact*(*r).s3.c2.re + (*s).s3.c2.re;
    (*r).s3.c2.im = fact*(*r).s3.c2.im + (*s).s3.c2.im;


    r = (spinor *) &R[ix].sp_dn;
    s = (spinor *) &S[ix].sp_dn;
    
    (*r).s0.c0.re = fact*(*r).s0.c0.re + (*s).s0.c0.re;
    (*r).s0.c0.im = fact*(*r).s0.c0.im + (*s).s0.c0.im;
    (*r).s0.c1.re = fact*(*r).s0.c1.re + (*s).s0.c1.re;
    (*r).s0.c1.im = fact*(*r).s0.c1.im + (*s).s0.c1.im;
    (*r).s0.c2.re = fact*(*r).s0.c2.re + (*s).s0.c2.re;
    (*r).s0.c2.im = fact*(*r).s0.c2.im + (*s).s0.c2.im;
    
    (*r).s1.c0.re = fact*(*r).s1.c0.re + (*s).s1.c0.re;
    (*r).s1.c0.im = fact*(*r).s1.c0.im + (*s).s1.c0.im;
    (*r).s1.c1.re = fact*(*r).s1.c1.re + (*s).s1.c1.re;
    (*r).s1.c1.im = fact*(*r).s1.c1.im + (*s).s1.c1.im;
    (*r).s1.c2.re = fact*(*r).s1.c2.re + (*s).s1.c2.re;
    (*r).s1.c2.im = fact*(*r).s1.c2.im + (*s).s1.c2.im;         
    
    (*r).s2.c0.re = fact*(*r).s2.c0.re + (*s).s2.c0.re;
    (*r).s2.c0.im = fact*(*r).s2.c0.im + (*s).s2.c0.im;
    (*r).s2.c1.re = fact*(*r).s2.c1.re + (*s).s2.c1.re;
    (*r).s2.c1.im = fact*(*r).s2.c1.im + (*s).s2.c1.im;
    (*r).s2.c2.re = fact*(*r).s2.c2.re + (*s).s2.c2.re;
    (*r).s2.c2.im = fact*(*r).s2.c2.im + (*s).s2.c2.im;
    
    (*r).s3.c0.re = fact*(*r).s3.c0.re + (*s).s3.c0.re;
    (*r).s3.c0.im = fact*(*r).s3.c0.im + (*s).s3.c0.im;
    (*r).s3.c1.re = fact*(*r).s3.c1.re + (*s).s3.c1.re;
    (*r).s3.c1.im = fact*(*r).s3.c1.im + (*s).s3.c1.im;
    (*r).s3.c2.re = fact*(*r).s3.c2.re + (*s).s3.c2.re;
    (*r).s3.c2.im = fact*(*r).s3.c2.im + (*s).s3.c2.im;
  }
}

#endif




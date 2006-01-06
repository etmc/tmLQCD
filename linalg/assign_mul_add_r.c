/* $Id$ */

#include <stdlib.h>
#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include "su3.h"
#include "sse.h"
#include "global.h"
#include "assign_mul_add_r.h"

#ifdef _STD_C99_COMPLEX_CHECKED
#include <complex.h>
#endif

#ifdef apenext
#include <topology.h>
#include <queue.h>
#endif


#if defined SSE2
/* k input , l output*/
void assign_mul_add_r(spinor * const S, const double c, spinor * const R, const int N) {

  int ix;
  su3_vector *s,*r;
  __asm__ __volatile__ ("movsd %0, %%xmm7 \n\t"
			"unpcklpd %%xmm7, %%xmm7"
			:
			:
			"m" (c));
  s=&S[0].s0;
  r=&R[0].s0;
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

#elif ((!defined _STD_C99_COMPLEX_CHECKED) && (!defined apenext))

/* R inoutput , c,S input*/
/*   (*R) = c*(*R) + (*S)        c is a real constant   */

void assign_mul_add_r(spinor * const R, const double c, spinor * const S, const int N) {

  int ix;
  static double fact;
  spinor *r,*s;
  
  fact=c;
  
  /* Change due to even-odd preconditioning : VOLUME   to VOLUME/2 */   
  for (ix = 0; ix < N; ix++) {
    r = R + ix;
    s = S + ix;
    
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

/*       printf("%1.16e %1.16e\n",(*r).s0.c0.re,(*r).s0.c0.im); */
/*       printf("%1.16e %1.16e\n",(*r).s0.c1.re,(*r).s0.c1.im); */
/*       printf("%1.16e %1.16e\n",(*r).s0.c2.re,(*r).s0.c2.im); */

/*       printf("%1.16e %1.16e\n",(*r).s1.c0.re,(*r).s1.c0.im); */
/*       printf("%1.16e %1.16e\n",(*r).s1.c1.re,(*r).s1.c1.im); */
/*       printf("%1.16e %1.16e\n",(*r).s1.c2.re,(*r).s1.c2.im); */

/*       printf("%1.16e %1.16e\n",(*r).s2.c0.re,(*r).s2.c0.im); */
/*       printf("%1.16e %1.16e\n",(*r).s2.c1.re,(*r).s2.c1.im); */
/*       printf("%1.16e %1.16e\n",(*r).s2.c2.re,(*r).s2.c2.im); */

/*       printf("%1.16e %1.16e\n",(*r).s3.c0.re,(*r).s3.c0.im); */
/*       printf("%1.16e %1.16e\n",(*r).s3.c1.re,(*r).s3.c1.im); */
/*       printf("%1.16e %1.16e\n",(*r).s3.c2.re,(*r).s3.c2.im); */
  }
}

#elif ((defined _STD_C99_COMPLEX_CHECKED) && (!defined apenext))

/*   (*R) = c*(*R) + (*S)        c is a real constant   */

void assign_mul_add_r(spinor * const R, const double c, spinor * const S, const int N) {

  register int ix=0;
  register double fact;
  register spinor *rPointer,*sPointer;

  fact=c;

  rPointer = R;
  sPointer = S;

  /* Change due to even-odd preconditioning : VOLUME   to VOLUME/2 */   
  do {
    register spinor s, r;
    ix+=1;
    
    s = *(sPointer);
    r = *(rPointer);

    r.s0.c0 =fact*r.s0.c0+s.s0.c0;
    r.s0.c1 =fact*r.s0.c1+s.s0.c1;
    r.s0.c2 =fact*r.s0.c2+s.s0.c2;
   	   	 	   
    r.s1.c0 =fact*r.s1.c0+s.s1.c0;
    r.s1.c1 =fact*r.s1.c1+s.s1.c1;
    r.s1.c2 =fact*r.s1.c2+s.s1.c2;
   	    	 	   
    r.s2.c0 =fact*r.s2.c0+s.s2.c0;
    r.s2.c1 =fact*r.s2.c1+s.s2.c1;
    r.s2.c2 =fact*r.s2.c2+s.s2.c2;
   	   	 	   
    r.s3.c0 =fact*r.s3.c0+s.s3.c0;
    r.s3.c1 =fact*r.s3.c1+s.s3.c1;
    r.s3.c2 =fact*r.s3.c2+s.s3.c2;

    *(rPointer) = r;

    sPointer+=1;
    rPointer+=1;

  } while (ix<N);
}

#elif defined apenext

#define NOWHERE_COND(condition) ((condition) ? 0x0 : NOWHERE ) 

/*   (*R) = c*(*R) + (*S)        c is a real constant   */

void assign_mul_add_r(spinor * const R, const double c, spinor * const S, const int N) {

  register int ix=N;
  register double fact;
  register spinor *rPointer,*sPointer;

  fact=c;

  rPointer = R;
  sPointer = S;

  prefetch(*(rPointer));
  prefetch(*(sPointer));

  {
#pragma cache

    do {
      register spinor s, r;
      register spinor *aux;
      ix--;
    
      rPointer++;
      sPointer++;

      fetch(r);
      fetch(s);


      prefetch (*(rPointer+NOWHERE_COND(ix)));
      prefetch (*(sPointer+NOWHERE_COND(ix)));


      r.s0.c0 = fact*r.s0.c0+s.s0.c0;
      r.s0.c1 = fact*r.s0.c1+s.s0.c1;
      r.s0.c2 = fact*r.s0.c2+s.s0.c2;
   	        	 	   
      r.s1.c0 = fact*r.s1.c0+s.s1.c0;
      r.s1.c1 = fact*r.s1.c1+s.s1.c1;
      r.s1.c2 = fact*r.s1.c2+s.s1.c2;
   	        	 	   
      r.s2.c0 = fact*r.s2.c0+s.s2.c0;
      r.s2.c1 = fact*r.s2.c1+s.s2.c1;
      r.s2.c2 = fact*r.s2.c2+s.s2.c2;
   	        	 	   
      r.s3.c0 = fact*r.s3.c0+s.s3.c0;
      r.s3.c1 = fact*r.s3.c1+s.s3.c1;
      r.s3.c2 = fact*r.s3.c2+s.s3.c2;

      *(rPointer-1) = r;

    } while (ix>0);
  }
}
#endif



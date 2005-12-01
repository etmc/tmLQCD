/*******************************************************************************
 * $Id$
 *
 * File assign_mul_bra_add_mul_ket_add.c
 *
 *   void assign_mul_bra_add_mul_ket_add
 *   (spinor * const R,spinor * const S,spinor * const U,const double c1,const double c2)
 *     (*R) = c2*(*R + c1*(*S)) + (*U)  with c1 and c2 complex variables
 *
 * Author: Thomas Chiarappa
 *         Thomas.Chiarappa@mib.infn.it
 *
 *******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "global.h"
#include "su3.h"

#include "assign_mul_bra_add_mul_ket_add_bi.h"



/* R inoutput, S input, U input, c1 input, c2 input */
void assign_mul_bra_add_mul_ket_add_bi(bispinor * const R, bispinor * const S, bispinor * const U, const complex c1, const complex c2, const int N) { 

  int ix;
  spinor *r,*s,*u,*w,*w_;
  double a1,b1,a2,b2;
#if ( defined SSE || defined SSE2 )
   w_=calloc(N+1, sizeof(spinor));
   w = (spinor *)(((unsigned long int)(w_)+ALIGN_BASE)&~ALIGN_BASE);
#else
   w_=calloc(N, sizeof(spinor));
   w = w_;
#endif

  a1 = c1.re;
  b1 = c1.im;

  a2 = c2.re;
  b2 = c2.im;
  
  
  for (ix = 0;ix < N; ix++){

    r=(spinor *) &R[ix].sp_up;
    s=(spinor *) &S[ix].sp_up;
    u=(spinor *) &U[ix].sp_up;
    /* w=(spinor *) U + ix; */

    /* *W = *R + c1*(*S) */

    (*w).s0.c0.re=(*r).s0.c0.re + a1*(*s).s0.c0.re - b1*(*s).s0.c0.im;
    (*w).s0.c0.im=(*r).s0.c0.im + a1*(*s).s0.c0.im + b1*(*s).s0.c0.re;
    (*w).s0.c1.re=(*r).s0.c1.re + a1*(*s).s0.c1.re - b1*(*s).s0.c1.im;
    (*w).s0.c1.im=(*r).s0.c1.im + a1*(*s).s0.c1.im + b1*(*s).s0.c1.re;
    (*w).s0.c2.re=(*r).s0.c2.re + a1*(*s).s0.c2.re - b1*(*s).s0.c2.im;
    (*w).s0.c2.im=(*r).s0.c2.im + a1*(*s).s0.c2.im + b1*(*s).s0.c2.re;

    (*w).s1.c0.re=(*r).s1.c0.re + a1*(*s).s1.c0.re - b1*(*s).s1.c0.im;
    (*w).s1.c0.im=(*r).s1.c0.im + a1*(*s).s1.c0.im + b1*(*s).s1.c0.re;
    (*w).s1.c1.re=(*r).s1.c1.re + a1*(*s).s1.c1.re - b1*(*s).s1.c1.im;
    (*w).s1.c1.im=(*r).s1.c1.im + a1*(*s).s1.c1.im + b1*(*s).s1.c1.re;
    (*w).s1.c2.re=(*r).s1.c2.re + a1*(*s).s1.c2.re - b1*(*s).s1.c2.im;
    (*w).s1.c2.im=(*r).s1.c2.im + a1*(*s).s1.c2.im + b1*(*s).s1.c2.re;
		 
    (*w).s2.c0.re=(*r).s2.c0.re + a1*(*s).s2.c0.re - b1*(*s).s2.c0.im;
    (*w).s2.c0.im=(*r).s2.c0.im + a1*(*s).s2.c0.im + b1*(*s).s2.c0.re;
    (*w).s2.c1.re=(*r).s2.c1.re + a1*(*s).s2.c1.re - b1*(*s).s2.c1.im;
    (*w).s2.c1.im=(*r).s2.c1.im + a1*(*s).s2.c1.im + b1*(*s).s2.c1.re;
    (*w).s2.c2.re=(*r).s2.c2.re + a1*(*s).s2.c2.re - b1*(*s).s2.c2.im;
    (*w).s2.c2.im=(*r).s2.c2.im + a1*(*s).s2.c2.im + b1*(*s).s2.c2.re;
		 
    (*w).s3.c0.re=(*r).s3.c0.re + a1*(*s).s3.c0.re - b1*(*s).s3.c0.im;
    (*w).s3.c0.im=(*r).s3.c0.im + a1*(*s).s3.c0.im + b1*(*s).s3.c0.re;
    (*w).s3.c1.re=(*r).s3.c1.re + a1*(*s).s3.c1.re - b1*(*s).s3.c1.im;
    (*w).s3.c1.im=(*r).s3.c1.im + a1*(*s).s3.c1.im + b1*(*s).s3.c1.re;
    (*w).s3.c2.re=(*r).s3.c2.re + a1*(*s).s3.c2.re - b1*(*s).s3.c2.im;
    (*w).s3.c2.im=(*r).s3.c2.im + a1*(*s).s3.c2.im + b1*(*s).s3.c2.re;

    /* *R = c2*(*W) + *U  */

    (*r).s0.c0.re=a2*(*w).s0.c0.re - b2*(*w).s0.c0.im + (*u).s0.c0.re;
    (*r).s0.c1.re=a2*(*w).s0.c1.re - b2*(*w).s0.c1.im + (*u).s0.c1.re;
    (*r).s0.c2.re=a2*(*w).s0.c2.re - b2*(*w).s0.c2.im + (*u).s0.c2.re;

    (*r).s1.c0.re=a2*(*w).s1.c0.re - b2*(*w).s1.c0.im + (*u).s1.c0.re;
    (*r).s1.c1.re=a2*(*w).s1.c1.re - b2*(*w).s1.c1.im + (*u).s1.c1.re;
    (*r).s1.c2.re=a2*(*w).s1.c2.re - b2*(*w).s1.c2.im + (*u).s1.c2.re;

    (*r).s2.c0.re=a2*(*w).s2.c0.re - b2*(*w).s2.c0.im + (*u).s2.c0.re;
    (*r).s2.c1.re=a2*(*w).s2.c1.re - b2*(*w).s2.c1.im + (*u).s2.c1.re;
    (*r).s2.c2.re=a2*(*w).s2.c2.re - b2*(*w).s2.c2.im + (*u).s2.c2.re;

    (*r).s3.c0.re=a2*(*w).s3.c0.re - b2*(*w).s3.c0.im + (*u).s3.c0.re;
    (*r).s3.c1.re=a2*(*w).s3.c1.re - b2*(*w).s3.c1.im + (*u).s3.c1.re;
    (*r).s3.c2.re=a2*(*w).s3.c2.re - b2*(*w).s3.c2.im + (*u).s3.c2.re;


    (*r).s0.c0.im=b2*(*w).s0.c0.re + a2*(*w).s0.c0.im + (*u).s0.c0.im;
    (*r).s0.c1.im=b2*(*w).s0.c1.re + a2*(*w).s0.c1.im + (*u).s0.c1.im;
    (*r).s0.c2.im=b2*(*w).s0.c2.re + a2*(*w).s0.c2.im + (*u).s0.c2.im;

    (*r).s1.c0.im=b2*(*w).s1.c0.re + a2*(*w).s1.c0.im + (*u).s1.c0.im;
    (*r).s1.c1.im=b2*(*w).s1.c1.re + a2*(*w).s1.c1.im + (*u).s1.c1.im;
    (*r).s1.c2.im=b2*(*w).s1.c2.re + a2*(*w).s1.c2.im + (*u).s1.c2.im;

    (*r).s2.c0.im=b2*(*w).s2.c0.re + a2*(*w).s2.c0.im + (*u).s2.c0.im;
    (*r).s2.c1.im=b2*(*w).s2.c1.re + a2*(*w).s2.c1.im + (*u).s2.c1.im;
    (*r).s2.c2.im=b2*(*w).s2.c2.re + a2*(*w).s2.c2.im + (*u).s2.c2.im;

    (*r).s3.c0.im=b2*(*w).s3.c0.re + a2*(*w).s3.c0.im + (*u).s3.c0.im;
    (*r).s3.c1.im=b2*(*w).s3.c1.re + a2*(*w).s3.c1.im + (*u).s3.c1.im;
    (*r).s3.c2.im=b2*(*w).s3.c2.re + a2*(*w).s3.c2.im + (*u).s3.c2.im;


    r=(spinor *) &R[ix].sp_dn;
    s=(spinor *) &S[ix].sp_dn;
    u=(spinor *) &U[ix].sp_dn;
    /* w=(spinor *) U + ix; */

    /* *W = *R + c1*(*S) */

    (*w).s0.c0.re=(*r).s0.c0.re + a1*(*s).s0.c0.re - b1*(*s).s0.c0.im;
    (*w).s0.c0.im=(*r).s0.c0.im + a1*(*s).s0.c0.im + b1*(*s).s0.c0.re;
    (*w).s0.c1.re=(*r).s0.c1.re + a1*(*s).s0.c1.re - b1*(*s).s0.c1.im;
    (*w).s0.c1.im=(*r).s0.c1.im + a1*(*s).s0.c1.im + b1*(*s).s0.c1.re;
    (*w).s0.c2.re=(*r).s0.c2.re + a1*(*s).s0.c2.re - b1*(*s).s0.c2.im;
    (*w).s0.c2.im=(*r).s0.c2.im + a1*(*s).s0.c2.im + b1*(*s).s0.c2.re;

    (*w).s1.c0.re=(*r).s1.c0.re + a1*(*s).s1.c0.re - b1*(*s).s1.c0.im;
    (*w).s1.c0.im=(*r).s1.c0.im + a1*(*s).s1.c0.im + b1*(*s).s1.c0.re;
    (*w).s1.c1.re=(*r).s1.c1.re + a1*(*s).s1.c1.re - b1*(*s).s1.c1.im;
    (*w).s1.c1.im=(*r).s1.c1.im + a1*(*s).s1.c1.im + b1*(*s).s1.c1.re;
    (*w).s1.c2.re=(*r).s1.c2.re + a1*(*s).s1.c2.re - b1*(*s).s1.c2.im;
    (*w).s1.c2.im=(*r).s1.c2.im + a1*(*s).s1.c2.im + b1*(*s).s1.c2.re;
		 
    (*w).s2.c0.re=(*r).s2.c0.re + a1*(*s).s2.c0.re - b1*(*s).s2.c0.im;
    (*w).s2.c0.im=(*r).s2.c0.im + a1*(*s).s2.c0.im + b1*(*s).s2.c0.re;
    (*w).s2.c1.re=(*r).s2.c1.re + a1*(*s).s2.c1.re - b1*(*s).s2.c1.im;
    (*w).s2.c1.im=(*r).s2.c1.im + a1*(*s).s2.c1.im + b1*(*s).s2.c1.re;
    (*w).s2.c2.re=(*r).s2.c2.re + a1*(*s).s2.c2.re - b1*(*s).s2.c2.im;
    (*w).s2.c2.im=(*r).s2.c2.im + a1*(*s).s2.c2.im + b1*(*s).s2.c2.re;
		 
    (*w).s3.c0.re=(*r).s3.c0.re + a1*(*s).s3.c0.re - b1*(*s).s3.c0.im;
    (*w).s3.c0.im=(*r).s3.c0.im + a1*(*s).s3.c0.im + b1*(*s).s3.c0.re;
    (*w).s3.c1.re=(*r).s3.c1.re + a1*(*s).s3.c1.re - b1*(*s).s3.c1.im;
    (*w).s3.c1.im=(*r).s3.c1.im + a1*(*s).s3.c1.im + b1*(*s).s3.c1.re;
    (*w).s3.c2.re=(*r).s3.c2.re + a1*(*s).s3.c2.re - b1*(*s).s3.c2.im;
    (*w).s3.c2.im=(*r).s3.c2.im + a1*(*s).s3.c2.im + b1*(*s).s3.c2.re;

    /* *R = c2*(*W) + *U  */

    (*r).s0.c0.re=a2*(*w).s0.c0.re - b2*(*w).s0.c0.im + (*u).s0.c0.re;
    (*r).s0.c1.re=a2*(*w).s0.c1.re - b2*(*w).s0.c1.im + (*u).s0.c1.re;
    (*r).s0.c2.re=a2*(*w).s0.c2.re - b2*(*w).s0.c2.im + (*u).s0.c2.re;

    (*r).s1.c0.re=a2*(*w).s1.c0.re - b2*(*w).s1.c0.im + (*u).s1.c0.re;
    (*r).s1.c1.re=a2*(*w).s1.c1.re - b2*(*w).s1.c1.im + (*u).s1.c1.re;
    (*r).s1.c2.re=a2*(*w).s1.c2.re - b2*(*w).s1.c2.im + (*u).s1.c2.re;

    (*r).s2.c0.re=a2*(*w).s2.c0.re - b2*(*w).s2.c0.im + (*u).s2.c0.re;
    (*r).s2.c1.re=a2*(*w).s2.c1.re - b2*(*w).s2.c1.im + (*u).s2.c1.re;
    (*r).s2.c2.re=a2*(*w).s2.c2.re - b2*(*w).s2.c2.im + (*u).s2.c2.re;

    (*r).s3.c0.re=a2*(*w).s3.c0.re - b2*(*w).s3.c0.im + (*u).s3.c0.re;
    (*r).s3.c1.re=a2*(*w).s3.c1.re - b2*(*w).s3.c1.im + (*u).s3.c1.re;
    (*r).s3.c2.re=a2*(*w).s3.c2.re - b2*(*w).s3.c2.im + (*u).s3.c2.re;


    (*r).s0.c0.im=b2*(*w).s0.c0.re + a2*(*w).s0.c0.im + (*u).s0.c0.im;
    (*r).s0.c1.im=b2*(*w).s0.c1.re + a2*(*w).s0.c1.im + (*u).s0.c1.im;
    (*r).s0.c2.im=b2*(*w).s0.c2.re + a2*(*w).s0.c2.im + (*u).s0.c2.im;

    (*r).s1.c0.im=b2*(*w).s1.c0.re + a2*(*w).s1.c0.im + (*u).s1.c0.im;
    (*r).s1.c1.im=b2*(*w).s1.c1.re + a2*(*w).s1.c1.im + (*u).s1.c1.im;
    (*r).s1.c2.im=b2*(*w).s1.c2.re + a2*(*w).s1.c2.im + (*u).s1.c2.im;

    (*r).s2.c0.im=b2*(*w).s2.c0.re + a2*(*w).s2.c0.im + (*u).s2.c0.im;
    (*r).s2.c1.im=b2*(*w).s2.c1.re + a2*(*w).s2.c1.im + (*u).s2.c1.im;
    (*r).s2.c2.im=b2*(*w).s2.c2.re + a2*(*w).s2.c2.im + (*u).s2.c2.im;

    (*r).s3.c0.im=b2*(*w).s3.c0.re + a2*(*w).s3.c0.im + (*u).s3.c0.im;
    (*r).s3.c1.im=b2*(*w).s3.c1.re + a2*(*w).s3.c1.im + (*u).s3.c1.im;
    (*r).s3.c2.im=b2*(*w).s3.c2.re + a2*(*w).s3.c2.im + (*u).s3.c2.im;


  }
  free(w_);

}

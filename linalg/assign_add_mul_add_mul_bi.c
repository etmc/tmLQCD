/*******************************************************************************
 * $Id$
 *
 * File assign_add_mul_add_mul.c
 *
 *   void assign_add_mul_add_mul(spinor * const R,spinor * const S,spinor * const U,const complex c1,const complex c2)
 *     (*R) = (*R) + c1*(*S) + c2*(*U) with c1 and c2 complex variables
 *
 * Author: Thomas Chiarappa
 *         Thomas.Chiarappa@mib.infn.it
 * 
 *******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "su3.h"
#include "global.h"

#include "assign_add_mul_add_mul_bi.h"


/* S,U input, R inoutput, c1,c2 input */
void assign_add_mul_add_mul_bi(bispinor * const R, bispinor * const S, bispinor * const U, const complex c1, const complex c2, const int N){


  int ix;
  spinor *r,*s,*u;
  double a1,b1,a2,b2;

  a1 = c1.re;
  b1 = c1.im;

  a2 = c2.re;
  b2 = c2.im;

  for (ix=0;ix<N;ix++){

    r=(spinor *) &R[ix].sp_up;
    s=(spinor *) &S[ix].sp_up;
    u=(spinor *) &U[ix].sp_up;
    
    (*r).s0.c0.re+=a1*(*s).s0.c0.re - b1*(*s).s0.c0.im + a2*(*u).s0.c0.re - b2*(*u).s0.c0.im;
    (*r).s0.c1.re+=a1*(*s).s0.c1.re - b1*(*s).s0.c1.im + a2*(*u).s0.c1.re - b2*(*u).s0.c1.im;
    (*r).s0.c2.re+=a1*(*s).s0.c2.re - b1*(*s).s0.c2.im + a2*(*u).s0.c2.re - b2*(*u).s0.c2.im;

    (*r).s1.c0.re+=a1*(*s).s1.c0.re - b1*(*s).s1.c0.im + a2*(*u).s1.c0.re - b2*(*u).s1.c0.im;
    (*r).s1.c1.re+=a1*(*s).s1.c1.re - b1*(*s).s1.c1.im + a2*(*u).s1.c1.re - b2*(*u).s1.c1.im;
    (*r).s1.c2.re+=a1*(*s).s1.c2.re - b1*(*s).s1.c2.im + a2*(*u).s1.c2.re - b2*(*u).s1.c2.im;

    (*r).s2.c0.re+=a1*(*s).s2.c0.re - b1*(*s).s2.c0.im + a2*(*u).s2.c0.re - b2*(*u).s2.c0.im;
    (*r).s2.c1.re+=a1*(*s).s2.c1.re - b1*(*s).s2.c1.im + a2*(*u).s2.c1.re - b2*(*u).s2.c1.im;
    (*r).s2.c2.re+=a1*(*s).s2.c2.re - b1*(*s).s2.c2.im + a2*(*u).s2.c2.re - b2*(*u).s2.c2.im;

    (*r).s3.c0.re+=a1*(*s).s3.c0.re - b1*(*s).s3.c0.im + a2*(*u).s3.c0.re - b2*(*u).s3.c0.im;
    (*r).s3.c1.re+=a1*(*s).s3.c1.re - b1*(*s).s3.c1.im + a2*(*u).s3.c1.re - b2*(*u).s3.c1.im;
    (*r).s3.c2.re+=a1*(*s).s3.c2.re - b1*(*s).s3.c2.im + a2*(*u).s3.c2.re - b2*(*u).s3.c2.im;



    (*r).s0.c0.im+=a1*(*s).s0.c0.im + b1*(*s).s0.c0.re + a2*(*u).s0.c0.im + b2*(*u).s0.c0.re;
    (*r).s0.c1.im+=a1*(*s).s0.c1.im + b1*(*s).s0.c1.re + a2*(*u).s0.c1.im + b2*(*u).s0.c1.re;
    (*r).s0.c2.im+=a1*(*s).s0.c2.im + b1*(*s).s0.c2.re + a2*(*u).s0.c2.im + b2*(*u).s0.c2.re;

    (*r).s1.c0.im+=a1*(*s).s1.c0.im + b1*(*s).s1.c0.re + a2*(*u).s1.c0.im + b2*(*u).s1.c0.re;
    (*r).s1.c1.im+=a1*(*s).s1.c1.im + b1*(*s).s1.c1.re + a2*(*u).s1.c1.im + b2*(*u).s1.c1.re;
    (*r).s1.c2.im+=a1*(*s).s1.c2.im + b1*(*s).s1.c2.re + a2*(*u).s1.c2.im + b2*(*u).s1.c2.re;

    (*r).s2.c0.im+=a1*(*s).s2.c0.im + b1*(*s).s2.c0.re + a2*(*u).s2.c0.im + b2*(*u).s2.c0.re;
    (*r).s2.c1.im+=a1*(*s).s2.c1.im + b1*(*s).s2.c1.re + a2*(*u).s2.c1.im + b2*(*u).s2.c1.re;
    (*r).s2.c2.im+=a1*(*s).s2.c2.im + b1*(*s).s2.c2.re + a2*(*u).s2.c2.im + b2*(*u).s2.c2.re;

    (*r).s3.c0.im+=a1*(*s).s3.c0.im + b1*(*s).s3.c0.re + a2*(*u).s3.c0.im + b2*(*u).s3.c0.re;
    (*r).s3.c1.im+=a1*(*s).s3.c1.im + b1*(*s).s3.c1.re + a2*(*u).s3.c1.im + b2*(*u).s3.c1.re;
    (*r).s3.c2.im+=a1*(*s).s3.c2.im + b1*(*s).s3.c2.re + a2*(*u).s3.c2.im + b2*(*u).s3.c2.re;



    r=(spinor *) &R[ix].sp_dn;
    s=(spinor *) &S[ix].sp_dn;
    u=(spinor *) &U[ix].sp_dn;
    
    (*r).s0.c0.re+=a1*(*s).s0.c0.re - b1*(*s).s0.c0.im + a2*(*u).s0.c0.re - b2*(*u).s0.c0.im;
    (*r).s0.c1.re+=a1*(*s).s0.c1.re - b1*(*s).s0.c1.im + a2*(*u).s0.c1.re - b2*(*u).s0.c1.im;
    (*r).s0.c2.re+=a1*(*s).s0.c2.re - b1*(*s).s0.c2.im + a2*(*u).s0.c2.re - b2*(*u).s0.c2.im;

    (*r).s1.c0.re+=a1*(*s).s1.c0.re - b1*(*s).s1.c0.im + a2*(*u).s1.c0.re - b2*(*u).s1.c0.im;
    (*r).s1.c1.re+=a1*(*s).s1.c1.re - b1*(*s).s1.c1.im + a2*(*u).s1.c1.re - b2*(*u).s1.c1.im;
    (*r).s1.c2.re+=a1*(*s).s1.c2.re - b1*(*s).s1.c2.im + a2*(*u).s1.c2.re - b2*(*u).s1.c2.im;

    (*r).s2.c0.re+=a1*(*s).s2.c0.re - b1*(*s).s2.c0.im + a2*(*u).s2.c0.re - b2*(*u).s2.c0.im;
    (*r).s2.c1.re+=a1*(*s).s2.c1.re - b1*(*s).s2.c1.im + a2*(*u).s2.c1.re - b2*(*u).s2.c1.im;
    (*r).s2.c2.re+=a1*(*s).s2.c2.re - b1*(*s).s2.c2.im + a2*(*u).s2.c2.re - b2*(*u).s2.c2.im;

    (*r).s3.c0.re+=a1*(*s).s3.c0.re - b1*(*s).s3.c0.im + a2*(*u).s3.c0.re - b2*(*u).s3.c0.im;
    (*r).s3.c1.re+=a1*(*s).s3.c1.re - b1*(*s).s3.c1.im + a2*(*u).s3.c1.re - b2*(*u).s3.c1.im;
    (*r).s3.c2.re+=a1*(*s).s3.c2.re - b1*(*s).s3.c2.im + a2*(*u).s3.c2.re - b2*(*u).s3.c2.im;



    (*r).s0.c0.im+=a1*(*s).s0.c0.im + b1*(*s).s0.c0.re + a2*(*u).s0.c0.im + b2*(*u).s0.c0.re;
    (*r).s0.c1.im+=a1*(*s).s0.c1.im + b1*(*s).s0.c1.re + a2*(*u).s0.c1.im + b2*(*u).s0.c1.re;
    (*r).s0.c2.im+=a1*(*s).s0.c2.im + b1*(*s).s0.c2.re + a2*(*u).s0.c2.im + b2*(*u).s0.c2.re;

    (*r).s1.c0.im+=a1*(*s).s1.c0.im + b1*(*s).s1.c0.re + a2*(*u).s1.c0.im + b2*(*u).s1.c0.re;
    (*r).s1.c1.im+=a1*(*s).s1.c1.im + b1*(*s).s1.c1.re + a2*(*u).s1.c1.im + b2*(*u).s1.c1.re;
    (*r).s1.c2.im+=a1*(*s).s1.c2.im + b1*(*s).s1.c2.re + a2*(*u).s1.c2.im + b2*(*u).s1.c2.re;

    (*r).s2.c0.im+=a1*(*s).s2.c0.im + b1*(*s).s2.c0.re + a2*(*u).s2.c0.im + b2*(*u).s2.c0.re;
    (*r).s2.c1.im+=a1*(*s).s2.c1.im + b1*(*s).s2.c1.re + a2*(*u).s2.c1.im + b2*(*u).s2.c1.re;
    (*r).s2.c2.im+=a1*(*s).s2.c2.im + b1*(*s).s2.c2.re + a2*(*u).s2.c2.im + b2*(*u).s2.c2.re;

    (*r).s3.c0.im+=a1*(*s).s3.c0.im + b1*(*s).s3.c0.re + a2*(*u).s3.c0.im + b2*(*u).s3.c0.re;
    (*r).s3.c1.im+=a1*(*s).s3.c1.im + b1*(*s).s3.c1.re + a2*(*u).s3.c1.im + b2*(*u).s3.c1.re;
    (*r).s3.c2.im+=a1*(*s).s3.c2.im + b1*(*s).s3.c2.re + a2*(*u).s3.c2.im + b2*(*u).s3.c2.re;



  }

}

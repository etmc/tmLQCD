/*******************************************************************************
 * $Id$
 *
 * File square_norm.c
 *
 *   double square_norm(spinor * const P )
 *     Returns the square norm of *P
 *
 *******************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "su3.h"
#include "sse.h"
#include "square_norm.h"
#ifdef MPI
#include <mpi.h>
#endif

double square_norm(spinor * const P, const int N) {
  int ix;
  static double ks,kc,ds,tr,ts,tt;
  spinor *s;
  
  ks = 0.0;
  kc = 0.0;
  
  /* Change due to even-odd preconditioning : VOLUME   to VOLUME/2 */   
  for (ix  =  0; ix < N; ix++) {
    s = P + ix;
    
    ds = (*s).s0.c0.re*(*s).s0.c0.re + (*s).s0.c0.im*(*s).s0.c0.im + 
         (*s).s0.c1.re*(*s).s0.c1.re + (*s).s0.c1.im*(*s).s0.c1.im + 
         (*s).s0.c2.re*(*s).s0.c2.re + (*s).s0.c2.im*(*s).s0.c2.im + 
         (*s).s1.c0.re*(*s).s1.c0.re + (*s).s1.c0.im*(*s).s1.c0.im + 
         (*s).s1.c1.re*(*s).s1.c1.re + (*s).s1.c1.im*(*s).s1.c1.im + 
         (*s).s1.c2.re*(*s).s1.c2.re + (*s).s1.c2.im*(*s).s1.c2.im + 
         (*s).s2.c0.re*(*s).s2.c0.re + (*s).s2.c0.im*(*s).s2.c0.im + 
         (*s).s2.c1.re*(*s).s2.c1.re + (*s).s2.c1.im*(*s).s2.c1.im + 
         (*s).s2.c2.re*(*s).s2.c2.re + (*s).s2.c2.im*(*s).s2.c2.im + 
         (*s).s3.c0.re*(*s).s3.c0.re + (*s).s3.c0.im*(*s).s3.c0.im + 
         (*s).s3.c1.re*(*s).s3.c1.re + (*s).s3.c1.im*(*s).s3.c1.im + 
         (*s).s3.c2.re*(*s).s3.c2.re + (*s).s3.c2.im*(*s).s3.c2.im;
    
    tr = ds + kc;
    ts = tr + ks;
    tt = ts-ks;
    ks = ts;
    kc = tr-tt;
  }
  kc = ks + kc;
#ifdef MPI
  MPI_Allreduce(&kc, &ks, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return ks;
#else
  return kc;
#endif
}

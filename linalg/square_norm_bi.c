/*******************************************************************************
 * $Id$
 *
 * File square_norm_bi.c
 *
 *   double square_norm_bi(bispinor * const P )
 *     Returns the square norm of *P
 *
 * Author: Thomas Chiarappa
 *         Thomas.Chiarappa@mib.infn.it
 * 
 *******************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#ifdef MPI
# include <mpi.h>
#endif
#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include "su3.h"
#include "sse.h"
#include "square_norm_bi.h"

double square_norm_bi(bispinor * const P, const int N) {
  int ix;
  static double ks,kc,ds,tr,ts,tt;
  spinor *s, *t;
  
  ks = 0.0;
  kc = 0.0;
  
  /* Change due to even-odd preconditioning : VOLUME   to VOLUME/2 */   
  for (ix  =  0; ix < N; ix++) {
    s = &P[ix].sp_up;
    t = &P[ix].sp_dn;

    
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
         (*s).s3.c2.re*(*s).s3.c2.re + (*s).s3.c2.im*(*s).s3.c2.im +
         (*t).s0.c0.re*(*t).s0.c0.re + (*t).s0.c0.im*(*t).s0.c0.im + 
         (*t).s0.c1.re*(*t).s0.c1.re + (*t).s0.c1.im*(*t).s0.c1.im + 
         (*t).s0.c2.re*(*t).s0.c2.re + (*t).s0.c2.im*(*t).s0.c2.im + 
         (*t).s1.c0.re*(*t).s1.c0.re + (*t).s1.c0.im*(*t).s1.c0.im + 
         (*t).s1.c1.re*(*t).s1.c1.re + (*t).s1.c1.im*(*t).s1.c1.im + 
         (*t).s1.c2.re*(*t).s1.c2.re + (*t).s1.c2.im*(*t).s1.c2.im + 
         (*t).s2.c0.re*(*t).s2.c0.re + (*t).s2.c0.im*(*t).s2.c0.im + 
         (*t).s2.c1.re*(*t).s2.c1.re + (*t).s2.c1.im*(*t).s2.c1.im + 
         (*t).s2.c2.re*(*t).s2.c2.re + (*t).s2.c2.im*(*t).s2.c2.im + 
         (*t).s3.c0.re*(*t).s3.c0.re + (*t).s3.c0.im*(*t).s3.c0.im + 
         (*t).s3.c1.re*(*t).s3.c1.re + (*t).s3.c1.im*(*t).s3.c1.im + 
         (*t).s3.c2.re*(*t).s3.c2.re + (*t).s3.c2.im*(*t).s3.c2.im;
    
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

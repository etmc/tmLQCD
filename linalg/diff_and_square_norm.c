/* $Id$ */

#include <stdlib.h>
#include "su3.h"
#ifdef MPI
#include <mpi.h>
#endif
#include "diff_and_square_norm.h"

double diff_and_square_norm(spinor * const Q, spinor * const R, const int N) {
  int ix;
  static double ks,kc,ds,tr,ts,tt;
  spinor *q,*r;
  
  ks=0.0;
  kc=0.0;
  
  /* Change due to even-odd preconditioning : VOLUME   to VOLUME/2 */   
  for (ix = 0; ix < N; ix++) {
    q=Q+ix;
    r=R+ix;
    
    (*q).s0.c0.re = (*r).s0.c0.re-(*q).s0.c0.re;
    (*q).s0.c0.im = (*r).s0.c0.im-(*q).s0.c0.im;
    (*q).s0.c1.re = (*r).s0.c1.re-(*q).s0.c1.re;
    (*q).s0.c1.im = (*r).s0.c1.im-(*q).s0.c1.im;
    (*q).s0.c2.re = (*r).s0.c2.re-(*q).s0.c2.re;
    (*q).s0.c2.im = (*r).s0.c2.im-(*q).s0.c2.im;
    
    ds = 
      (*q).s0.c0.re*(*q).s0.c0.re+(*q).s0.c0.im*(*q).s0.c0.im+
      (*q).s0.c1.re*(*q).s0.c1.re+(*q).s0.c1.im*(*q).s0.c1.im+
      (*q).s0.c2.re*(*q).s0.c2.re+(*q).s0.c2.im*(*q).s0.c2.im;
    
    (*q).s1.c0.re = (*r).s1.c0.re-(*q).s1.c0.re;
    (*q).s1.c0.im = (*r).s1.c0.im-(*q).s1.c0.im;
    (*q).s1.c1.re = (*r).s1.c1.re-(*q).s1.c1.re;
    (*q).s1.c1.im = (*r).s1.c1.im-(*q).s1.c1.im;
    (*q).s1.c2.re = (*r).s1.c2.re-(*q).s1.c2.re;
    (*q).s1.c2.im = (*r).s1.c2.im-(*q).s1.c2.im;         
    
    ds+= 
      (*q).s1.c0.re*(*q).s1.c0.re+(*q).s1.c0.im*(*q).s1.c0.im+
      (*q).s1.c1.re*(*q).s1.c1.re+(*q).s1.c1.im*(*q).s1.c1.im+
      (*q).s1.c2.re*(*q).s1.c2.re+(*q).s1.c2.im*(*q).s1.c2.im;
    
    (*q).s2.c0.re = (*r).s2.c0.re-(*q).s2.c0.re;
    (*q).s2.c0.im = (*r).s2.c0.im-(*q).s2.c0.im;
    (*q).s2.c1.re = (*r).s2.c1.re-(*q).s2.c1.re;
    (*q).s2.c1.im = (*r).s2.c1.im-(*q).s2.c1.im;
    (*q).s2.c2.re = (*r).s2.c2.re-(*q).s2.c2.re;
    (*q).s2.c2.im = (*r).s2.c2.im-(*q).s2.c2.im;         
    
    ds+ = 
      (*q).s2.c0.re*(*q).s2.c0.re+(*q).s2.c0.im*(*q).s2.c0.im+
      (*q).s2.c1.re*(*q).s2.c1.re+(*q).s2.c1.im*(*q).s2.c1.im+
      (*q).s2.c2.re*(*q).s2.c2.re+(*q).s2.c2.im*(*q).s2.c2.im;
    
    (*q).s3.c0.re = (*r).s3.c0.re-(*q).s3.c0.re;
    (*q).s3.c0.im = (*r).s3.c0.im-(*q).s3.c0.im;
    (*q).s3.c1.re = (*r).s3.c1.re-(*q).s3.c1.re;
    (*q).s3.c1.im = (*r).s3.c1.im-(*q).s3.c1.im;
    (*q).s3.c2.re = (*r).s3.c2.re-(*q).s3.c2.re;
    (*q).s3.c2.im = (*r).s3.c2.im-(*q).s3.c2.im;
    
    ds+ = 
      (*q).s3.c0.re*(*q).s3.c0.re+(*q).s3.c0.im*(*q).s3.c0.im+
      (*q).s3.c1.re*(*q).s3.c1.re+(*q).s3.c1.im*(*q).s3.c1.im+
      (*q).s3.c2.re*(*q).s3.c2.re+(*q).s3.c2.im*(*q).s3.c2.im;
    
    tr = ds+kc;
    ts = tr+ks;
    tt = ts-ks;
    ks = ts;
    kc = tr-tt;
  }
  kc = ks+kc;
#ifdef MPI
  MPI_Allreduce(&kc, &ks, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return ks;
#else
  return kc;
#endif

}

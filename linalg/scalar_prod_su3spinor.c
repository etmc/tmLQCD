/***********************************************************************
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008 Carsten Urbach
 *
 * This file is part of tmLQCD.
 *
 * tmLQCD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * tmLQCD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with tmLQCD.  If not, see <http://www.gnu.org/licenses/>.
 ***********************************************************************/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#ifdef MPI
#include <mpi.h>
#endif
#include "su3.h"
#include "scalar_prod_su3spinor.h"

#ifdef WITHLAPH
complex_spinor scalar_prod_su3spinor(su3_vector * const S, spinor * const R, const int N, const int parallel){
  int ix;
  static double ks,kc,ds,tr,ts,tt;
  su3_vector *s,*r;
  complex_spinor c;
#ifdef MPI
  complex_spinor d;
#endif

  /* sc0.re */

  ks=0.0;
  kc=0.0;
  for (ix = 0; ix < N; ix++)
    {
      s=(su3_vector *) S + ix;
      r=&(R[ix].s0);
    
      ds=(*r).c0.re*(*s).c0.re+(*r).c0.im*(*s).c0.im+
	(*r).c1.re*(*s).c1.re+(*r).c1.im*(*s).c1.im+
	(*r).c2.re*(*s).c2.re+(*r).c2.im*(*s).c2.im;

      /* Kahan Summation */    
      tr=ds+kc;
      ts=tr+ks;
      tt=ts-ks;
      ks=ts;
      kc=tr-tt;
    }
  kc=ks+kc;
  c.sc0.re = kc;

  /* sc0.im */

  ks=0.0;
  kc=0.0;
  for (ix=0;ix<N;ix++)
    {
      s=(su3_vector *) S + ix;
      r=&(R[ix].s0);
    
      ds=-(*r).c0.re*(*s).c0.im+(*r).c0.im*(*s).c0.re-
	(*r).c1.re*(*s).c1.im+(*r).c1.im*(*s).c1.re-
	(*r).c2.re*(*s).c2.im+(*r).c2.im*(*s).c2.re;
    
      tr=ds+kc;
      ts=tr+ks;
      tt=ts-ks;
      ks=ts;
      kc=tr-tt;
    }
  kc=ks+kc;
  c.sc0.im = kc;

  /* sc1.re */

  ks=0.0;
  kc=0.0;
  for (ix = 0; ix < N; ix++)
    {
      s=(su3_vector *) S + ix;
      r=&(R[ix].s0);
    
      ds=(*r).c0.re*(*s).c0.re+(*r).c0.im*(*s).c0.im+
	(*r).c1.re*(*s).c1.re+(*r).c1.im*(*s).c1.im+
	(*r).c2.re*(*s).c2.re+(*r).c2.im*(*s).c2.im;

      /* Kahan Summation */    
      tr=ds+kc;
      ts=tr+ks;
      tt=ts-ks;
      ks=ts;
      kc=tr-tt;
    }
  kc=ks+kc;
  c.sc1.re = kc;

  /* sc1.im */

  ks=0.0;
  kc=0.0;
  for (ix=0;ix<N;ix++)
    {
      s=(su3_vector *) S + ix;
      r=&(R[ix].s1);
    
      ds=-(*r).c0.re*(*s).c0.im+(*r).c0.im*(*s).c0.re-
	(*r).c1.re*(*s).c1.im+(*r).c1.im*(*s).c1.re-
	(*r).c2.re*(*s).c2.im+(*r).c2.im*(*s).c2.re;
    
      tr=ds+kc;
      ts=tr+ks;
      tt=ts-ks;
      ks=ts;
      kc=tr-tt;
    }
  kc=ks+kc;
  c.sc1.im = kc;

  /* sc2.re */

  ks=0.0;
  kc=0.0;
  for (ix = 0; ix < N; ix++)
    {
      s=(su3_vector *) S + ix;
      r=&(R[ix].s2);
    
      ds=(*r).c0.re*(*s).c0.re+(*r).c0.im*(*s).c0.im+
	(*r).c1.re*(*s).c1.re+(*r).c1.im*(*s).c1.im+
	(*r).c2.re*(*s).c2.re+(*r).c2.im*(*s).c2.im;

      /* Kahan Summation */    
      tr=ds+kc;
      ts=tr+ks;
      tt=ts-ks;
      ks=ts;
      kc=tr-tt;
    }
  kc=ks+kc;
  c.sc2.re = kc;

  /* sc2.im */

  ks=0.0;
  kc=0.0;
  for (ix=0;ix<N;ix++)
    {
      s=(su3_vector *) S + ix;
      r=&(R[ix].s2);
    
      ds=-(*r).c0.re*(*s).c0.im+(*r).c0.im*(*s).c0.re-
	(*r).c1.re*(*s).c1.im+(*r).c1.im*(*s).c1.re-
	(*r).c2.re*(*s).c2.im+(*r).c2.im*(*s).c2.re;
    
      tr=ds+kc;
      ts=tr+ks;
      tt=ts-ks;
      ks=ts;
      kc=tr-tt;
    }
  kc=ks+kc;
  c.sc2.im = kc;

  /* sc3.re */

  ks=0.0;
  kc=0.0;
  for (ix = 0; ix < N; ix++)
    {
      s=(su3_vector *) S + ix;
      r=&(R[ix].s3);
    
      ds=(*r).c0.re*(*s).c0.re+(*r).c0.im*(*s).c0.im+
	(*r).c1.re*(*s).c1.re+(*r).c1.im*(*s).c1.im+
	(*r).c2.re*(*s).c2.re+(*r).c2.im*(*s).c2.im;

      /* Kahan Summation */    
      tr=ds+kc;
      ts=tr+ks;
      tt=ts-ks;
      ks=ts;
      kc=tr-tt;
    }
  kc=ks+kc;
  c.sc3.re = kc;

  /* sc3.im */

  ks=0.0;
  kc=0.0;
  for (ix=0;ix<N;ix++)
    {
      s=(su3_vector *) S + ix;
      r=&(R[ix].s3);
    
      ds=-(*r).c0.re*(*s).c0.im+(*r).c0.im*(*s).c0.re-
	(*r).c1.re*(*s).c1.im+(*r).c1.im*(*s).c1.re-
	(*r).c2.re*(*s).c2.im+(*r).c2.im*(*s).c2.re;
    
      tr=ds+kc;
      ts=tr+ks;
      tt=ts-ks;
      ks=ts;
      kc=tr-tt;
    }
  kc=ks+kc;
  c.sc3.im = kc;

#ifdef MPI
  if(parallel == 1) {
    d = c;
    MPI_Allreduce(&d, &c, 4, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD); //???
  }
#endif

  return(c);
}
#endif // WITHLAPH

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
#include "scalar_prod.h"

/*  <S,R>=S^* times R */
_Complex double scalar_prod(spinor * const S, spinor * const R, const int N, const int parallel) {
  int ix;
  _Complex double ALIGN ds,ks,kc,tr,ts,tt;
  spinor *s,*r;
  _Complex double c;
#ifdef MPI
  _Complex double d;
#endif
  
  ks=0.0;
  kc=0.0;
#if (defined BGL && defined XLC)
  __alignx(16, S);
  __alignx(16, R);
#endif  
  for (ix = 0; ix < N; ix++)
  {
    s=(spinor *) S + ix;
    r=(spinor *) R + ix;
    
    ds = r->s0.c0 * conj(s->s0.c0) + r->s0.c1 * conj(s->s0.c1) + r->s0.c2 * conj(s->s0.c2) +
         r->s1.c0 * conj(s->s1.c0) + r->s1.c1 * conj(s->s1.c1) + r->s1.c2 * conj(s->s1.c2) +
	 r->s2.c0 * conj(s->s2.c0) + r->s2.c1 * conj(s->s2.c1) + r->s2.c2 * conj(s->s2.c2) + 
         r->s3.c0 * conj(s->s3.c0) + r->s3.c1 * conj(s->s3.c1) + r->s3.c2 * conj(s->s3.c2);

    /* Kahan Summation */    
    tr=ds+kc;
    ts=tr+ks;
    tt=ts-ks;
    ks=ts;
    kc=tr-tt;
  }
  kc=ks+kc;

  c = kc;

#ifdef MPI
  if(parallel == 1)
  {
    d = c;
    MPI_Allreduce(&d, &c, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
  }
#endif
  return(c);
}

#ifdef WITHLAPH
_Complex double scalar_prod_su3vect(su3_vector * const S, su3_vector * const R, const int N, const int parallel)
{
  double ALIGN ks, ds, tr, ts, tt;
  su3_vector *s, *r;
  _Complex double c;
#ifdef MPI
  _Complex double d;
#endif

  /* Real Part */

  ks = 0.0;
  c  = 0.0;
  for (int ix = 0; ix < N; ++ix)
    {
      s = (su3_vector *) S + ix;
      r = (su3_vector *) R + ix;
    
      ds = r->c0 * conj(s->c0) + r->c1 * conj(s->c1) + r->c2 * conj(s->c2);

      /* Kahan Summation */    
      tr = ds + c;
      ts = tr + ks;
      tt = ts - ks;
      ks = ts;
      c  = tr - tt;
    }
  c = ks + c;

#ifdef MPI
  if(parallel == 1)
  {
    d = c;
    MPI_Allreduce(&d, &c, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
  }
#endif
  return(c);
}
#endif

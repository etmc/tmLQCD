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
  static _Complex double ks, kc, ds, tr, ts, tt;
  su3_vector *s, *r;
  complex_spinor c;
#ifdef MPI
  complex_spinor d;
#endif

  /* sc0 */
  ks = 0.0;
  kc = 0.0;
  for (ix = 0; ix < N; ix++)
  {
    s = (su3_vector *) S + ix;
    r = &(R[ix].s0);
  
    ds = r->c0 * conj(s->c0) + r->c1 * conj(s->c1) + r->c2 * conj(s->c2);

    /* Kahan Summation */    
    tr = ds + kc;
    ts = tr + ks;
    tt = ts - ks;
    ks = ts;
    kc = tr - tt;
  }
  kc = ks + kc;
  c.sc0 = kc;

  /* sc1 */
  ks = 0.0;
  kc = 0.0;
  for (ix = 0; ix < N; ix++)
  {
    s = (su3_vector *) S + ix;
    r = &(R[ix].s1);
  
    ds = r->c0 * conj(s->c0) + r->c1 * conj(s->c1) + r->c2 * conj(s->c2);

    /* Kahan Summation */    
    tr = ds + kc;
    ts = tr + ks;
    tt = ts - ks;
    ks = ts;
    kc = tr - tt;
  }
  kc = ks + kc;
  c.sc1 = kc;

  /* sc2 */
  ks = 0.0;
  kc = 0.0;
  for (ix = 0; ix < N; ix++)
  {
    s = (su3_vector *) S + ix;
    r = &(R[ix].s2);
  
    ds = r->c0 * conj(s->c0) + r->c1 * conj(s->c1) + r->c2 * conj(s->c2);

    /* Kahan Summation */    
    tr = ds + kc;
    ts = tr + ks;
    tt = ts - ks;
    ks = ts;
    kc = tr - tt;
  }
  kc = ks + kc;
  c.sc2 = kc;

  /* sc3 */
  ks = 0.0;
  kc = 0.0;
  for (ix = 0; ix < N; ix++)
  {
    s = (su3_vector *) S + ix;
    r = &(R[ix].s3);
  
    ds = r->c0 * conj(s->c0) + r->c1 * conj(s->c1) + r->c2 * conj(s->c2);

    /* Kahan Summation */    
    tr = ds + kc;
    ts = tr + ks;
    tt = ts - ks;
    ks = ts;
    kc = tr - tt;
  }
  kc = ks + kc;
  c.sc3 = kc;

#ifdef MPI
  if(parallel == 1) {
    d = c;
    MPI_Allreduce(&d, &c, 4, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD); //???
  }
#endif

  return(c);
}
#endif // WITHLAPH

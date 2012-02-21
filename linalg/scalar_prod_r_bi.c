/***********************************************************************
 * Copyright (C) 2006 Thomas Chiarappa
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
/*******************************************************************************
 *
 * File scalar_prod_r_bi.c
 *
 *   double scalar_prod_r_bi(bispinor * const S,bispinor * const R, const int N)
 *     Returns the real part of the scalar product (*R,*S)
 *
 * 
 * Author: Thomas Chiarappa
 *         Thomas.Chiarappa@mib.infn.it
 * 
 *******************************************************************************/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#ifdef MPI
# include <mpi.h>
#endif
#include "su3.h"
#include "scalar_prod_r_bi.h"


/*  R input, S input */
double scalar_prod_r_bi(bispinor * const S,bispinor * const R, const int N){


  int ix;
  static double ks,kc,ds,tr,ts,tt;
  spinor *s,*r,*t,*u;
  
  ks=0.0;
  kc=0.0;
  
  for (ix=0;ix<N;ix++){

    s = (spinor *) &S[ix].sp_up;
    r = (spinor *) &R[ix].sp_up;
    t = (spinor *) &S[ix].sp_dn;
    u = (spinor *) &R[ix].sp_dn;
    
    ds = r->s0.c0 * conj(s->s0.c0) + r->s0.c1 * conj(s->s0.c1) + r->s0.c2 * conj(s->s0.c2) +
    r->s1.c0 * conj(s->s1.c0) + r->s1.c1 * conj(s->s1.c1) + r->s1.c2 * conj(s->s1.c2) +
    r->s2.c0 * conj(s->s2.c0) + r->s2.c1 * conj(s->s2.c1) + r->s2.c2 * conj(s->s2.c2) +
    r->s3.c0 * conj(s->s3.c0) + r->s3.c1 * conj(s->s3.c1) + r->s3.c2 * conj(s->s3.c2) +
    u->t0.c0 * conj(t->t0.c0) + u->t0.c1 * conj(t->t0.c1) + u->t0.c2 * conj(t->t0.c2) +
    u->t1.c0 * conj(t->t1.c0) + u->t1.c1 * conj(t->t1.c1) + u->t1.c2 * conj(t->t1.c2) +
    u->t2.c0 * conj(t->t2.c0) + u->t2.c1 * conj(t->t2.c1) + u->t2.c2 * conj(t->t2.c2) +
    u->t3.c0 * conj(t->t3.c0) + u->t3.c1 * conj(t->t3.c1) + u->t3.c2 * conj(t->t3.c2);
    
    tr = ds + kc;
    ts = tr + ks;
    tt = ts-ks;
    ks = ts;
    kc = tr-tt;
  }
  kc = ks + kc;

#if defined MPI
  
  MPI_Allreduce(&kc, &ks, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return ks;

#endif
  
  return kc;
}
 


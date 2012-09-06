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
/*******************************************************************************
 *
 * File assign_mul_bra_add_mul_r.c
 *
 *   void assign_mul_bra_add_mul_r(spinor * const R,const double c0, const double c,spinor * const S)
 *     (*R) = c0*(*R + c*(*S))
 *
 *******************************************************************************/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#ifdef OMP
# include <omp.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "su3.h"
#include "assign_mul_bra_add_mul_r.h"

/*  R output, S input, c0 input, c input */
void assign_mul_bra_add_mul_r(spinor * const R,const double c0, const double c,spinor * const S, const int N){
#ifdef OMP
#pragma omp parallel
  {
#endif

  int ix;
  double ALIGN fact0,fact;
  spinor *r,*s;
  
  fact=c;
  fact0=c0;


#ifdef OMP
#pragma omp for
#endif
  for (ix = 0;ix < N; ++ix)
  {
    r=(spinor *) R + ix;
    s=(spinor *) S + ix;
    
    r->s0.c0 = fact0 * (r->s0.c0 + fact * s->s0.c0);
    r->s0.c1 = fact0 * (r->s0.c1 + fact * s->s0.c1);
    r->s0.c2 = fact0 * (r->s0.c2 + fact * s->s0.c2);
    
    r->s1.c0 = fact0 * (r->s1.c0 + fact * s->s1.c0);
    r->s1.c1 = fact0 * (r->s1.c1 + fact * s->s1.c1);
    r->s1.c2 = fact0 * (r->s1.c2 + fact * s->s1.c2);

    r->s2.c0 = fact0 * (r->s2.c0 + fact * s->s2.c0);
    r->s2.c1 = fact0 * (r->s2.c1 + fact * s->s2.c1);
    r->s2.c2 = fact0 * (r->s2.c2 + fact * s->s2.c2);

    r->s3.c0 = fact0 * (r->s3.c0 + fact * s->s3.c0);
    r->s3.c1 = fact0 * (r->s3.c1 + fact * s->s3.c1);
    r->s3.c2 = fact0 * (r->s3.c2 + fact * s->s3.c2);
  }
#ifdef OMP
  } /* OpenMP closing brace */
#endif
}

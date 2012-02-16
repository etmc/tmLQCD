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
#include "complex.h"
#include "mattimesvec.h"

/* v = M*w                         */
/* v,w complex vectors of length N */
/* M a NxN complex matrix with     */
/* leading dimension ldM >= N      */
/* we should provide special SSE2  */
/* and BG/P versions               */

void mattimesvec(_Complex double * const v, _Complex double * const M, _Complex double * const w, 
		 const int N, const int ldM) 
{
  for(int i = 0; i < N; ++i)
  {
    v[i] = M[i*ldM] * w[0];
    for(int j = 1; j < N; ++j)
      v[i] += M[i*ldM + j] * w[j];
  }
  return;
}

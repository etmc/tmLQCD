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
#include <stdio.h>
#include <string.h>
#ifdef MPI
# include <mpi.h>
#endif
#ifdef OMP
# include <omp.h>
#endif
#include "global.h"
#include "su3.h"
#include "convert_eo_to_lexic.h"

void convert_eo_to_lexic(spinor * const P, spinor * const s, spinor * const r) {
#ifdef OMP
#pragma omp parallel
  {
#endif

  int x, y, z, t, i, ix;
  spinor * p = NULL;

#ifdef OMP
#pragma omp for
#endif
  for(x = 0; x < LX; x++) {
    for(y = 0; y < LY; y++) {
      for(z = 0; z < LZ; z++) {
	for(t = 0; t < T; t++) {
	  ix = g_ipt[t][x][y][z];
	  i = g_lexic2eosub[ ix ];
	  if((t+x+y+z+g_proc_coords[3]*LZ+g_proc_coords[2]*LY 
	      + g_proc_coords[0]*T+g_proc_coords[1]*LX)%2 == 0) {
	    p = s;
	  }
	  else {
	    p = r;
	  }
	  memcpy((P+ix), (p+i), sizeof(spinor));
	}
      }
    }
  }

#ifdef OMP
  } /*OpenMP closing brace */
#endif

  return;
}

/*
 *      P: spinor with full volume 
 *      s: new spinor even 
 *      r: new spinor odd 
 */
void convert_lexic_to_eo(spinor * const s, spinor * const r, spinor * const P) {
#ifdef OMP
#pragma omp parallel
  {
#endif

  int x, y, z, t, i, ix;
  spinor * p = NULL;

#ifdef OMP
#pragma omp for
#endif
  for(x = 0; x < LX; x++) {
    for(y = 0; y < LY; y++) {
      for(z = 0; z < LZ; z++) {
	for(t = 0; t < T; t++) {
	  ix = g_ipt[t][x][y][z];
	  i = g_lexic2eosub[ ix ];
	  if((t+x+y+z+g_proc_coords[3]*LZ+g_proc_coords[2]*LY 
	      + g_proc_coords[0]*T+g_proc_coords[1]*LX)%2 == 0) {
	    p = s;
	  }
	  else {
	    p = r;
	  }
	  memcpy((p+i), (P+ix), sizeof(spinor));
	}
      }
    }
  }

#ifdef OMP
  } /* OpenMP closing brace */
#endif

  return;
}

/***********************************************************************
* Copyright (C) 2012 Carsten Urbach
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
#include "global.h"
#ifdef MPI
# include <mpi.h>
#endif
#include "su3.h"
#include "io/spinor_write_stdout.h"


void spinor_write_stdout(spinor * const s) {
  int X, Y, Z, t0, id = 0, ix, iy;
  int coords[4];

  for(int t = 0; t < g_nproc_t*T; t++) {
    t0 = t - g_proc_coords[0]*T;
    coords[0] = t / T;
    for(int x = 0; x < g_nproc_x*LX; x++) {
      X = x - g_proc_coords[1]*LX; 
      coords[1] = x / LX;
      for(int y = 0; y < g_nproc_y*LY; y++) {
	Y = y - g_proc_coords[2]*LY;
	coords[2] = y / LY;
	for(int z = 0; z < g_nproc_z*LZ; z++) {
	  Z = z - g_proc_coords[3]*LZ;
	  coords[3] = z / LZ;
#ifdef MPI
	  MPI_Cart_rank(g_cart_grid, coords, &id);
#endif
	  if((t+x+y+z)%2 == 0 && g_cart_id == id) {
	    ix = g_lexic2eosub[ g_ipt[t0][X][Y][Z] ];
	    iy = t*(g_nproc_x*LX*g_nproc_y*LY*g_nproc_z*LZ) +
	      x*(g_nproc_y*LY*g_nproc_z*LZ) +
	      y*(g_nproc_z*LZ) + z;
	    printf(" %d %d %d %d %d, %d %d %d %d: %e %e sp\n",
		   iy, t, x, y, z, t0, X, Y, Z, 
		   creal(s[ix].s0.c0), cimag(s[ix].s0.c0));
	    fflush(stdout);
	  }
#ifdef MPI
	  MPI_Barrier(MPI_COMM_WORLD);
#endif
	}
      }
    }
  }
  return;
}

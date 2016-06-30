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

#include "gauge.ih"
#include <stdio.h>

// not parallelised!

int write_luscher_gauge_binary(const double plaq, char* filename, su3 ** const gf) {
#ifdef TM_USE_MPI
  fprintf(stdout, "Luescher DD-HMC format for gauges not implemented for MPI! Not writing anything!\n");
#else
  FILE * ofs;
  int ix;

  ofs = fopen(filename, "w");

  fwrite(&T, sizeof(int), 1, ofs);
  fwrite(&L, sizeof(int), 1, ofs);
  fwrite(&L, sizeof(int), 1, ofs);
  fwrite(&L, sizeof(int), 1, ofs);
  fwrite(&plaq, sizeof(double), 1, ofs);

  for(int t = 0; t < T; t++) {
    for(int x = 0; x < LX; x++) {
      for(int y = 0; y < LY; y++) {
	for(int z = 0; z < LZ; z++) {
	  ix = g_ipt[t][x][y][z];
	  // if odd
	  if((t + x + y + z)%2 == 1) {
	    for(int mu =0; mu < 4; mu++) {
	      // forward direction
	      fwrite(&gf[ix][mu], sizeof(double), 18, ofs);
	      // backward direction
	      fwrite(&gf[ g_idn[ix][mu] ][mu], sizeof(double), 18, ofs);
	    }
	  }
	}
      }
    }
  }
  fclose(ofs);

#endif
  return(0);
}

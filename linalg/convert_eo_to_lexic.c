/* $Id$ */

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#ifdef MPI
# include <mpi.h>
#endif
#include "global.h"
#include "su3.h"
#include "convert_eo_to_lexic.h"

void convert_eo_to_lexic(spinor * const P, spinor * const s, spinor * const r) {
  int x, y, z, t, i, ix;
  spinor * p = NULL;

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
  return;
}

void convert_lexic_to_eo(spinor * const s, spinor * const r, spinor * const P) {
  int x, y, z, t, i, ix;
  spinor * p = NULL;

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
  return;
}

/* $Id$ */

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "global.h"
#include "start.h"
#include "ranlxd.h"
#include "source_generation.h"

#ifndef M_PI
# define M_PI           3.14159265358979323846
#endif

void source_generation_nucleon(spinor * const P, spinor * const Q, 
		       const int is, const int ic,
		       const int t, const int nt, const int nx, 
		       const int sample, const int nstore) {

  double rnumber, si=0., co=0.;
  int rlxd_state[105];
  int reset = 0, seed, r, tt, lt, xx, lx, yy, ly, zz, lz, ix;
  int coords[4], id=0, i;
  complex * p = NULL;
  const double s0=0.;
  const double c0=1.;
  const double s1=sin(2.*M_PI/3.);
  const double c1=cos(2.*M_PI/3.);
  const double s2=sin(4.*M_PI/3.);
  const double c2=cos(4.*M_PI/3.);

  zero_spinor_field(P,VOLUME/2);
  zero_spinor_field(Q,VOLUME/2);

  /* save the ranlxd_state if neccessary */
  if(ranlxd_init == 1) {
    rlxd_get(rlxd_state);
  }

  /* Compute the seed */
  seed =(int) abs(1+sample + t*10*97 + nstore*100*53);

  rlxd_init(1, seed);

  for(tt = t; tt < T*g_nproc_t; tt+=nt) {
    lt = tt - g_proc_coords[0]*T;
    coords[0] = tt / T;
    ix = 0;
    for(xx = 0; xx < LX*g_nproc_x; xx++) {
      lx = xx - g_proc_coords[1]*LX;
      coords[1] = xx / LX;
      for(yy = 0; yy < LY*g_nproc_y; yy++) {
	ly = yy - g_proc_coords[2]*LY;
	coords[2] = yy / LY;
	for(zz = 0; zz < LZ*g_nproc_z; zz++) {
	  lz = zz - g_proc_coords[3]*LZ;
	  coords[3] = zz / LZ;
	  if(ix%nx == 0) {
#ifdef MPI
	    MPI_Cart_rank(g_cart_grid, coords, &id);
#endif
	    ranlxd(&rnumber, 1);
	    if(g_cart_id  == id) {
	      r = (int)floor(3.*rnumber);
	      if(r == 0) {
		si = s0;
		co = c0;
	      }
	      else if(r == 1) {
		si = s1;
		co = c1;
	      }
	      else {
		si = s2;
		co = c2;
	      }

	      i = g_lexic2eosub[ g_ipt[lt][lx][ly][lz] ];
	      if((lt+lx+ly+lz+g_proc_coords[3]*LZ+g_proc_coords[2]*LY 
		  + g_proc_coords[0]*T+g_proc_coords[1]*LX)%2 == 0) {
		p = (complex*)(P + i);
	      }
	      else {
		p = (complex*)(Q + i);
	      }

	      (*(p+3*is+ic)).re = co;
 	      (*(p+3*is+ic)).im = si;
	    }
	  }
	  ix++;
	}
      }
    }
  }

  /* reset the ranlxd if neccessary */
  if(reset) {
    rlxd_reset(rlxd_state);
  }
  return;
}

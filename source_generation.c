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

void extended_pion_source(spinor * const P, spinor * const Q,
			  spinor * const S, spinor * const R,
			  const int t,
			  const double px, const double py, const double pz) {
  int lt, lx, ly, lz, i, x, y, z, id=0;
  int coords[4];
  spinor * p, * q, r;
  complex efac;

  zero_spinor_field(P,VOLUME/2);
  zero_spinor_field(Q,VOLUME/2);

  lt = t - g_proc_coords[0]*T;
  coords[0] = t / T;
  for(x = 0; x < LX*g_nproc_x; x++) {
    lx = x - g_proc_coords[1]*LX;
    coords[1] = x / LX;
    for(y = 0; y < LY*g_nproc_y; y++) {
      ly = y - g_proc_coords[2]*LY;
      coords[2] = y / LY;
      for(z = 0; z < LZ*g_nproc_z; z++) {
	lz = z - g_proc_coords[3]*LZ;
	coords[3] = z / LZ;
#ifdef MPI
	MPI_Cart_rank(g_cart_grid, coords, &id);
#endif
	if(g_cart_id == id) {
	  efac.re= cos(px*x + py*y + pz*z);
	  efac.im=-sin(px*x + py*y + pz*z);

	  i = g_lexic2eosub[ g_ipt[lt][lx][ly][lz] ];
	  if((lt+lx+ly+lz+g_proc_coords[3]*LZ+g_proc_coords[2]*LY 
	      + g_proc_coords[0]*T+g_proc_coords[1]*LX)%2 == 0) {
	    p = (P + i);
	    q = (R + i);
	  }
	  else {
	    p = (Q + i);
	    q = (S + i);
	  }
	  _gamma5(r, (*q));
	  _spinor_mul_complex((*p),efac,r);
	}
      }
    }
  }
  return;
}

void source_generation_pion_only(spinor * const P, spinor * const Q,
				 const int t,
				 const int sample, const int nstore) {

  int reset = 0, i, x, y, z, is, ic, lt, lx, ly, lz, id=0;
  int coords[4], seed, r;
  double rnumber, si=0., co=0.;
  int rlxd_state[105];
  const double sqr2 = 1./sqrt(2.);
  complex * p = NULL;
  
  zero_spinor_field(P,VOLUME/2);
  zero_spinor_field(Q,VOLUME/2);

  /* save the ranlxd_state if neccessary */
  if(ranlxd_init == 1) {
    rlxd_get(rlxd_state);
    reset = 1;
  }

  /* Compute the seed */
  seed =(int) abs(1 + sample + t*10*97 + nstore*100*53);

  rlxd_init(1, seed);

  lt = t - g_proc_coords[0]*T;
  coords[0] = t / T;
  for(x = 0; x < LX*g_nproc_x; x++) {
    lx = x - g_proc_coords[1]*LX;
    coords[1] = x / LX;
    for(y = 0; y < LY*g_nproc_y; y++) {
      ly = y - g_proc_coords[2]*LY;
      coords[2] = y / LY;
      for(z = 0; z < LZ*g_nproc_z; z++) {
	lz = z - g_proc_coords[3]*LZ;
	coords[3] = z / LZ;
#ifdef MPI
	MPI_Cart_rank(g_cart_grid, coords, &id);
#endif
	for(is = 0; is < 4; is++) {
	  for(ic = 0; ic < 3; ic++) {
	    ranlxd(&rnumber, 1);
	    if(g_cart_id  == id) {
	      r = (int)floor(4.*rnumber);
	      if(r == 0) {
		si = sqr2;
		co = sqr2;
	      }
	      else if(r == 1) {
		si = -sqr2;
		co = sqr2;
	      }
	      else if(r==2) {
		si = sqr2;
		co = -sqr2;
	      }
	      else {
		si = -sqr2;
		co = -sqr2;
	      }
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
      }
    }
  }
	    
  /* reset the ranlxd if neccessary */
  if(reset) {
    rlxd_reset(rlxd_state);
  }
  return;
}

void source_generation_nucleon(spinor * const P, spinor * const Q, 
			       const int is, const int ic,
			       const int t, const int nt, const int nx, 
			       const int sample, const int nstore, 
			       const int meson) {

  double rnumber, si=0., co=0., sqr2;
  int rlxd_state[105];
  int reset = 0, seed, r, tt, lt, xx, lx, yy, ly, zz, lz;
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

  sqr2 = 1./sqrt(2.);
  /* save the ranlxd_state if neccessary */
  if(ranlxd_init == 1) {
    rlxd_get(rlxd_state);
    reset = 1;
  }

  /* Compute the seed */
  seed =(int) abs(1 + sample + t*10*97 + nstore*100*53);

  rlxd_init(1, seed);

  for(tt = t; tt < T*g_nproc_t; tt+=nt) {
    lt = tt - g_proc_coords[0]*T;
    coords[0] = tt / T;
    for(xx = 0; xx < LX*g_nproc_x; xx+=nx) {
      lx = xx - g_proc_coords[1]*LX;
      coords[1] = xx / LX;
      for(yy = 0; yy < LY*g_nproc_y; yy+=nx) {
	ly = yy - g_proc_coords[2]*LY;
	coords[2] = yy / LY;
	for(zz = 0; zz < LZ*g_nproc_z; zz+=nx) {
	  lz = zz - g_proc_coords[3]*LZ;
	  coords[3] = zz / LZ;
#ifdef MPI
	  MPI_Cart_rank(g_cart_grid, coords, &id);
#endif
	  ranlxd(&rnumber, 1);
	  if(g_cart_id  == id) {
	    if(meson) {
	      r = (int)floor(4.*rnumber);
	      if(r == 0) {
		si = sqr2;
		co = sqr2;
	      }
	      else if(r == 1) {
		si = -sqr2;
		co = sqr2;
	      }
	      else if(r==2) {
		si = sqr2;
		co = -sqr2;
	      }
	      else {
		si = -sqr2;
		co = -sqr2;
	      }
	    }
	    else {
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
      }
    }
  }

  /* reset the ranlxd if neccessary */
  if(reset) {
    rlxd_reset(rlxd_state);
  }
  return;
}

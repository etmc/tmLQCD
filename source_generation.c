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
# include<tmlqcd_config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "global.h"
#include "start.h"
#include "ranlxd.h"
#include "su3spinor.h"
#include "source_generation.h"
#include "solver/solver_field.h"
#include "linalg/convert_eo_to_lexic.c"

#ifndef M_PI
# define M_PI           3.14159265358979323846
#endif

/* Generates normal distributed random numbers */
/* using the box-muller method                 */
/* this is even standard normal distributed    */
/* so mean = 0, sd = 1                         */
void rnormal(double * r, const int n) 
{
  double u[2], s, l;
  int i;

  /* basic form, but trig. functions needed */
/*   for(i = 0; i < n; i+=2) { */
/*     ranlxd(u, 2); */
/*     l = sqrt(-2*log(u[0])); */
/*     r[i] = l*cos(2*M_PI*u[1]); */
/*     r[i+1] = l*sin(2*M_PI*u[1]); */
/*     printf("%f\n", r[i]); */
/*     printf("%f\n", r[i+1]); */
/*   } */
/*   return; */
  /* polar form, no trig. functions, but more random numbers */
  /* which one is faster? */
  for(i = 0; i < n; i += 2) {
    ranlxd(u, 2);
    u[0] = 2.*u[0] - 1.;
    u[1] = 2.*u[1] - 1.;
    s = u[0]*u[0]+u[1]*u[1];
    while(s == 0. || s > 1.) {
      ranlxd(u, 2);
      u[0] = 2.*u[0] - 1.;
      u[1] = 2.*u[1] - 1.;
      s = u[0]*u[0]+u[1]*u[1];
    }
    l = sqrt(-2.*log(s)/s);
    r[i] = u[0]*l;
    r[i+1] = u[1]*l;
  }
  return;
}

/* Generates a volume source with gaussian noise */
/* in all real and imaginary elements            */
/*                                               */
/* i.e. xi*.xi = 2                               */
/* is the normalisation                          */
/* this is corrected for in the contraction      */
/* codes                                         */
void gaussian_volume_source(spinor * const P, spinor * const Q,
			    const int sample, const int nstore, const int f) 
{
  int x, y, z, t, i, reset = 0, seed; 
  int rlxd_state[105];
  spinor * p;

  /* save the ranlxd_state if neccessary */
  if(ranlxd_init == 1) {
    rlxd_get(rlxd_state);
    reset = 1;
  }

  /* Compute the seed */
  seed =(int) abs(1 + sample + f*10*97 + nstore*100*53 + g_cart_id*13);

  rlxd_init(2, seed);

  for(t = 0; t < T; t++) {
    for(x = 0; x < LX; x++) {
      for(y =0; y < LY; y++) {
	for(z = 0; z < LZ; z++) {
	  i = g_lexic2eosub[ g_ipt[t][x][y][z] ];
	  if((t+x+y+z+g_proc_coords[3]*LZ+g_proc_coords[2]*LY 
	      + g_proc_coords[0]*T+g_proc_coords[1]*LX)%2 == 0) {
	    p = P + i;
	  }
	  else {
	    p = Q + i;
	  }
	  rnormal((double*)p, 24);
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

void extended_pion_source(spinor * const P, spinor * const Q,
			  spinor * const R, spinor * const S,
			  const int t0, const int ts,
			  const double px, const double py, const double pz) {
  int lt, lx, ly, lz, i, x, y, z, id=0, t;
  int coords[4];
  spinor * p, * q, r;
  _Complex double efac;

  zero_spinor_field(P,VOLUME/2);
  zero_spinor_field(Q,VOLUME/2);
  
  t=(ts + t0)%(g_nproc_t*T);
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
#ifdef TM_USE_MPI
	MPI_Cart_rank(g_cart_grid, coords, &id);
#endif
	if(g_cart_id == id) {
	  efac = cexp(-(px * x + py * y + pz * z) * I);

	  i = g_lexic2eosub[ g_ipt[lt][lx][ly][lz] ];
	  if((lt+lx+ly+lz+g_proc_coords[3]*LZ+g_proc_coords[2]*LY 
	      + g_proc_coords[0]*T+g_proc_coords[1]*LX)%2 == 0) {
	    p = P + i;
	    q = R + i;
	  }
	  else {
	    p = Q + i;
	    q = S + i;
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
				 const int t, const int sample, 
                                 const int nstore, const unsigned int _seed) {

  int reset = 0, i, x, y, z, is, ic, lt, lx, ly, lz, id=0;
  int coords[4], seed, r;
  double rnumber, si=0., co=0.;
  int rlxd_state[105];
  const double sqr2 = 1./sqrt(2.);
  _Complex double * p = NULL;
  
  zero_spinor_field(P,VOLUME/2);
  zero_spinor_field(Q,VOLUME/2);

  /* save the ranlxd_state if neccessary */
  if(ranlxd_init == 1) {
    rlxd_get(rlxd_state);
    reset = 1;
  }

  /* Compute the seed */
  seed =(int) abs(_seed + sample + t*10*97 + nstore*100*53);

  rlxd_init(2, seed);

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
#ifdef TM_USE_MPI
	MPI_Cart_rank(g_cart_grid, coords, &id);
#endif
	for(is = 0; is < 4; is++) {
	  for(ic = 0; ic < 3; ic++) {
	    ranlxd(&rnumber, 1);
	    if(g_cart_id  == id) {
	      r = (int)floor(4.*rnumber);
	      if(r == 0)
	      {
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
	    
	      i = g_lexic2eosub[ g_ipt[lt][lx][ly][lz] ];
	      if((lt+lx+ly+lz+g_proc_coords[3]*LZ+g_proc_coords[2]*LY 
		  + g_proc_coords[0]*T+g_proc_coords[1]*LX)%2 == 0) {
		p = (_Complex double*)(P + i);
	      }
	      else {
		p = (_Complex double*)(Q + i);
	      }
	      
	      (*(p+3*is+ic)) = co + si * I;
	    }
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

/* Florian Burger 4.11.2009 */
void source_generation_pion_zdir(spinor * const P, spinor * const Q,
                                 const int z,
                                 const int sample, const int nstore) {

  int reset = 0, i, x, y, t, is, ic, lt, lx, ly, lz, id=0;
  int coords[4], seed, r;
  double rnumber, si=0., co=0.;
  int rlxd_state[105];
  const double sqr2 = 1./sqrt(2.);
  _Complex double * p = NULL;
  
  zero_spinor_field(P,VOLUME/2);
  zero_spinor_field(Q,VOLUME/2);

  /* save the ranlxd_state if neccessary */
  if(ranlxd_init == 1) {
    rlxd_get(rlxd_state);
    reset = 1;
  }

  /* Compute the seed */
  seed =(int) abs(1 + sample + z*10*97 + nstore*100*53);

  rlxd_init(2, seed);
  lz = z - g_proc_coords[3]*LZ;
  coords[3] = z / LZ;
 for(t = 0; t < T*g_nproc_t; t++) {
   lt = t - g_proc_coords[0]*T;
   coords[0] = t / T;  
   for(x = 0; x < LX*g_nproc_x; x++) {
    lx = x - g_proc_coords[1]*LX;
    coords[1] = x / LX;
    for(y = 0; y < LY*g_nproc_y; y++) {
      ly = y - g_proc_coords[2]*LY;
      coords[2] = y / LY;

#ifdef TM_USE_MPI
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
            
              i = g_lexic2eosub[ g_ipt[lt][lx][ly][lz] ];
              if((lt+lx+ly+lz+g_proc_coords[3]*LZ+g_proc_coords[2]*LY 
                  + g_proc_coords[0]*T+g_proc_coords[1]*LX)%2 == 0) {
                p = (_Complex double*)(P + i);
              }
              else {
                p = (_Complex double*)(Q + i);
              }
              
	      (*(p+3*is+ic)) = co + si * I;
            }
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

/* end Florian Burger 4.11.2009 */





void source_generation_nucleon(spinor * const P, spinor * const Q, 
			       const int is, const int ic,
			       const int t, const int nt, const int nx, 
			       const int sample, const int nstore, 
			       const int meson) {

  double rnumber, si=0., co=0., sqr2;
  int rlxd_state[105];
  int reset = 0, seed, r, tt, lt, xx, lx, yy, ly, zz, lz;
  int coords[4], id=0, i;
  _Complex double * p = NULL;
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

  rlxd_init(2, seed);

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
#ifdef TM_USE_MPI
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
	      p = (_Complex double*)(P + i);
	    }
	    else {
	      p = (_Complex double*)(Q + i);
	    }

	    (*(p+3*is+ic)) = co + si * I;
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

void full_source_spinor_field_point(spinor * const full_spinor,
                                    const int is, const int ic,
                                    const int * const global_txyz_src_pos){
  spinor ** temp_cb_spinors;
  init_solver_field(&temp_cb_spinors, VOLUME/2, 2);
  eo_source_spinor_field_point(temp_cb_spinors[0], temp_cb_spinors[1], is, ic, global_txyz_src_pos);
  convert_eo_to_lexic(full_spinor, temp_cb_spinors[0], temp_cb_spinors[1]);
  finalize_solver(temp_cb_spinors, 2);
}

// create a point source at the lattice-global coordinates specified by the
// four-element array global_txyz_src_pos
void eo_source_spinor_field_point(spinor * const even_cb_spinor, 
                                  spinor * const odd_cb_spinor,
                                  const int is, const int ic,
                                  const int * const global_txyz_src_pos){
  zero_spinor_field(even_cb_spinor, VOLUME/2);
  zero_spinor_field(odd_cb_spinor, VOLUME/2);

  if( (global_txyz_src_pos[0] >= g_nproc_t * T || global_txyz_src_pos[0] < 0) ||
      (global_txyz_src_pos[1] >= g_nproc_x * LX || global_txyz_src_pos[1] < 0) ||
      (global_txyz_src_pos[2] >= g_nproc_y * LY || global_txyz_src_pos[2] < 0) ||
      (global_txyz_src_pos[3] >= g_nproc_z * LZ || global_txyz_src_pos[3] < 0) ){
    char error_message[500];
    snprintf(error_message,
             500,
             "Coordinates t=%d x=%d y=%d z=%d are outside the global lattice extent!\n",
             global_txyz_src_pos[0],
             global_txyz_src_pos[1],
             global_txyz_src_pos[2],
             global_txyz_src_pos[3]);
    fatal_error(error_message, "source_spinor_field_point");
  }

  if( (g_proc_coords[0] == global_txyz_src_pos[0] / T) &&
      (g_proc_coords[1] == global_txyz_src_pos[1] / LX) &&
      (g_proc_coords[2] == global_txyz_src_pos[2] / LY) &&
      (g_proc_coords[3] == global_txyz_src_pos[3] / LZ) ) {
    
    const int idx = g_ipt[ global_txyz_src_pos[0] % T ]
                           [ global_txyz_src_pos[1] % LX ]
                           [ global_txyz_src_pos[2] % LY ]
                           [ global_txyz_src_pos[3] % LZ ];

    const int eo_idx = g_lexic2eosub[ idx ];

    spinor * s = (g_lexic2eo[idx] < VOLUME/2 ? even_cb_spinor : odd_cb_spinor) +
                 eo_idx;

    /* put source to 1.0 */
    if (is==0){
      if      (ic==0) s->s0.c0 = 1.0;
      else if (ic==1) s->s0.c1 = 1.0;
      else if (ic==2) s->s0.c2 = 1.0;
    }
    else if (is==1){
      if      (ic==0) s->s1.c0 = 1.0;
      else if (ic==1) s->s1.c1 = 1.0;
      else if (ic==2) s->s1.c2 = 1.0;
    }
    else if (is==2){
      if      (ic==0) s->s2.c0 = 1.0;
      else if (ic==1) s->s2.c1 = 1.0;
      else if (ic==2) s->s2.c2 = 1.0;
    }
    else if (is==3){
      if      (ic==0) s->s3.c0 = 1.0;
      else if (ic==1) s->s3.c1 = 1.0;
      else if (ic==2) s->s3.c2 = 1.0;
    }
  }
}

void full_source_spinor_field_spin_diluted_oet_ts(spinor * const full_spinor,
                                                  const int src_ts,
                                                  const int src_d,
                                                  const int sample,
                                                  const int nstore,
                                                  const unsigned int oet_seed){
  spinor ** temp_cb_spinors;
  init_solver_field(&temp_cb_spinors, VOLUME/2, 2);
  eo_source_spinor_field_spin_diluted_oet_ts(temp_cb_spinors[0], 
                                             temp_cb_spinors[1],
                                             src_ts,
                                             src_d,
                                             sample,
                                             nstore,
                                             oet_seed);
  convert_eo_to_lexic(full_spinor, temp_cb_spinors[0], temp_cb_spinors[1]);
  finalize_solver(temp_cb_spinors, 2);
}

void eo_source_spinor_field_spin_diluted_oet_ts(spinor * const even_cb_spinor,
                                                spinor * const odd_cb_spinor,
                                                const int src_ts,
                                                const int src_d,
                                                const int sample,
                                                const int nstore,
                                                const unsigned int oet_seed){
  int reset = 0, i, x, y, z, ic, lt, lx, ly, lz, id=0;
  int coords[4], seed, r;
  double rnumber, si=0., co=0.;
  int rlxd_state[105];
  const double sqr2 = 1./sqrt(2.);
  _Complex double * p = NULL;
  
  zero_spinor_field(even_cb_spinor,VOLUME/2);
  zero_spinor_field(odd_cb_spinor,VOLUME/2);

  /* save the ranlxd_state if neccessary */
  if(ranlxd_init == 1) {
    rlxd_get(rlxd_state);
    reset = 1;
  }

  seed =(int) abs(oet_seed + sample + src_ts*10*97 + nstore);

  rlxd_init(2, seed);

  lt = src_ts - g_proc_coords[0]*T;
  coords[0] = src_ts / T;
  for(x = 0; x < LX*g_nproc_x; x++) {
    lx = x - g_proc_coords[1]*LX;
    coords[1] = x / LX;
    for(y = 0; y < LY*g_nproc_y; y++) {
      ly = y - g_proc_coords[2]*LY;
      coords[2] = y / LY;
      for(z = 0; z < LZ*g_nproc_z; z++) {
	      lz = z - g_proc_coords[3]*LZ;
	      coords[3] = z / LZ;
#ifdef TM_USE_MPI
	MPI_Cart_rank(g_cart_grid, coords, &id);
#endif
	    for(ic = 0; ic < 3; ic++) {
	      ranlxd(&rnumber, 1);
	      if(g_cart_id  == id) {
	        r = (int)floor(4.*rnumber);
	        if(r == 0){
	          si = sqr2;
	          co = sqr2;
	        } else if(r == 1) {
	          si = -sqr2;
	          co = sqr2;
	        } else if(r==2) {
	          si = sqr2;
	          co = -sqr2;
	        } else {
	          si = -sqr2;
	          co = -sqr2;
	        }
	        
	        i = g_lexic2eosub[ g_ipt[lt][lx][ly][lz] ];
	        if((lt+lx+ly+lz+g_proc_coords[3]*LZ+g_proc_coords[2]*LY 
	            + g_proc_coords[0]*T+g_proc_coords[1]*LX)%2 == 0) {
	          p = (_Complex double*)(even_cb_spinor + i);
	        } else {
	          p = (_Complex double*)(odd_cb_spinor + i);
	        }  
	          (*(p+3*src_d+ic)) = co + si * I;
	        }
	      } // ic
      } // z
    } // y
  } // x
	    
  /* reset the ranlxd if neccessary */
  if(reset) {
    rlxd_reset(rlxd_state);
  }
  return;

} 

/***********************************************************************
 *
 * Copyright (C) 2001 Martin Hasebusch
 *
 * some changes by C. Urbach 2002-2008
 *
 * Modified by Jenifer Gonzalez Lopez for the Schroedinger Functional
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
#include <math.h>
#include <time.h>
#include "global.h"
#include "su3.h"
#include "su3adj.h"
#include "su3spinor.h"
#include "expo.h"
#include "sse.h"
#include "xchange/xchange.h"
#include "get_rectangle_staples.h"
#include "gamma.h"
#include "get_staples.h"
#include "read_input.h"
#include "smearing/stout.h"

#include "ranlxd.h"
#include "start.h"
#include "phmc.h"
#include "hybrid_update.h"




/*----------------------------------------------------------------------------*/

/*******************************************
 *
 * This computes the contribution to
 * the Hamiltonian coming from the momenta
 *
 *******************************************/
double moment_energy(su3adj ** const momenta) {

  su3adj *xm;
  int i,mu;
  static double tt,tr,ts,kc,ks,sum;
  kc=0.; ks=0.;
  
  for(i=0;i<VOLUME;i++){
    for(mu=0;mu<4;mu++){
      xm=&momenta[i][mu];
      sum=(*xm).d1*(*xm).d1
	+(*xm).d2*(*xm).d2
	+(*xm).d3*(*xm).d3
	+(*xm).d4*(*xm).d4
	+(*xm).d5*(*xm).d5
	+(*xm).d6*(*xm).d6
	+(*xm).d7*(*xm).d7
	+(*xm).d8*(*xm).d8;
      tr=sum+kc;
      ts=tr+ks;
      tt=ts-ks;
      ks=ts;
      kc=tr-tt;
    }
  }
  /* from the loop I got: p^2 */
  /* the contribution to the E is however (p^2)/2: */
  kc=0.5*(ks+kc);
#ifdef MPI
  MPI_Allreduce(&kc, &ks, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return ks;
#else
  return kc;
#endif
}

/*----------------------------------------------------------------------------*/

/**************************************
 *
 * Initialises the momenta
 * with the gaussian distribution
 *
 **************************************/
double init_momenta(const int repro, su3adj ** const momenta) {
  
  su3adj *xm;
  int i, mu, t0, x, y, z, X, Y, Z, t, id = 0;
  int coords[4];
#ifdef MPI
  int k;
  int rlxd_state[105];
#endif
  double ALIGN yy[8];
  double ALIGN tt, tr, ts, kc = 0., ks = 0., sum;
  
  if(repro) {
#ifdef MPI
    if(g_proc_id == 0) {
      rlxd_get(rlxd_state);
    }
    MPI_Bcast(rlxd_state, 105, MPI_INT, 0, MPI_COMM_WORLD);
    rlxd_reset(rlxd_state);
#endif
    for(t0 = 0; t0 < g_nproc_t*T; t0++) {
      t = t0 - T*g_proc_coords[0];
      coords[0] = t0 / T;
      for(x = 0; x < g_nproc_x*LX; x++) {
	X = x - g_proc_coords[1]*LX;
	coords[1] = x / LX;
	for(y = 0; y < g_nproc_y*LY; y++) {
	  Y = y - g_proc_coords[2]*LY;
	  coords[2] = y / LY;
	  for(z = 0; z < g_nproc_z*LZ; z++) {
	    Z = z - g_proc_coords[3]*LZ;
	    coords[3] = z / LZ;
#ifdef MPI
	    MPI_Cart_rank(g_cart_grid, coords, &id);
#endif
	    if(g_cart_id == id) i = g_ipt[t][X][Y][Z];
	    for(mu = 0; mu < 4; mu++) {
	      gauss_vector(yy,8);
	      if(g_cart_id == id) {
		sum = 0.;
		xm = &momenta[i][mu];
		(*xm).d1 = 1.4142135623731*yy[0];
		(*xm).d2 = 1.4142135623731*yy[1];
		sum += (*xm).d1*(*xm).d1+(*xm).d2*(*xm).d2;
		(*xm).d3 = 1.4142135623731*yy[2];
		(*xm).d4 = 1.4142135623731*yy[3];
		sum += (*xm).d3*(*xm).d3+(*xm).d4*(*xm).d4;
		(*xm).d5 = 1.4142135623731*yy[4];
		(*xm).d6 = 1.4142135623731*yy[5];
		sum += (*xm).d5*(*xm).d5+(*xm).d6*(*xm).d6;
		(*xm).d7 = 1.4142135623731*yy[6];
		(*xm).d8 = 1.4142135623731*yy[7];
		sum += (*xm).d7*(*xm).d7+(*xm).d8*(*xm).d8;
		tr = sum+kc;
		ts = tr+ks;
		tt = ts-ks;
		ks = ts;
		kc = tr-tt;
	      }
	    }
	  }
	}
      }
    }
    kc=0.5*(ks+kc);
  }
  else {
    for(i = 0; i < VOLUME; i++) { 
      for(mu = 0; mu < 4; mu++) {
	sum=0.;
	xm=&momenta[i][mu];
	gauss_vector(yy,8);
	(*xm).d1=1.4142135623731*yy[0];
	(*xm).d2=1.4142135623731*yy[1];
	sum+=(*xm).d1*(*xm).d1+(*xm).d2*(*xm).d2;
	(*xm).d3=1.4142135623731*yy[2];
	(*xm).d4=1.4142135623731*yy[3];
	sum+=(*xm).d3*(*xm).d3+(*xm).d4*(*xm).d4;
	(*xm).d5=1.4142135623731*yy[4];
	(*xm).d6=1.4142135623731*yy[5];
	sum+=(*xm).d5*(*xm).d5+(*xm).d6*(*xm).d6;
	(*xm).d7=1.4142135623731*yy[6];
	(*xm).d8=1.4142135623731*yy[7];
	sum+=(*xm).d7*(*xm).d7+(*xm).d8*(*xm).d8;
	tr=sum+kc;
	ts=tr+ks;
	tt=ts-ks;
	ks=ts;
	kc=tr-tt;
      }
    }
    kc=0.5*(ks+kc);
  }
#ifdef MPI
  MPI_Allreduce(&kc, &ks, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return ks;
#endif
  return kc;
}


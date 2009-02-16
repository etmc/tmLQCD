/***********************************************************************
 * $Id$
 *
 * Copyright (C) 2001 Martin Hasebusch
 *               2002,2003,2004,2005,2006,2007,2008 Carsten Urbach
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
#include "xchange.h"
#include "get_rectangle_staples.h"
#include "gamma.h"
#include "get_staples.h"
#include "read_input.h"
#include "stout_smear.h"
#include "stout_smear_force.h"
#include "ranlxd.h"
#include "start.h"
#include "phmc.h"
#include "hybrid_update.h"


/***********************************
 *
 * Computes the new gauge field
 *
 ***********************************/



/*----------------------------------------------------------------------------*/

void update_gauge(double step) {

  int i,mu;
  static su3 v,w;
  su3 *z;
  static su3adj deriv;
  su3adj *xm;
#ifdef _KOJAK_INST
#pragma pomp inst begin(updategauge)
#endif

  for(i = 0; i < VOLUME; i++) { 
    for(mu = 0; mu < 4; mu++){
      xm=&moment[i][mu];
      z=&g_gauge_field[i][mu];
      _assign_const_times_mom(deriv, step, *xm);
      v=restoresu3( exposu3(deriv) );
      _su3_times_su3(w, v, *z);
      _su3_assign(*z, w);
    }
  }

#ifdef MPI
  /* for parallelization */
  xchange_gauge();
#endif
  /*
   * The backward copy of the gauge field
   * is not updated here!
   */
  g_update_gauge_copy = 1;
  g_update_gauge_energy = 1;
  g_update_rectangle_energy = 1;
  return;
#ifdef _KOJAK_INST
#pragma pomp inst end(updategauge)
#endif
}


/*----------------------------------------------------------------------------*/

/*******************************************
 *
 * This computes the contribution to
 * the Hamiltonian coming from the momenta
 *
 *******************************************/
double moment_energy() {

  su3adj *xm;
  int i,mu;
  static double tt,tr,ts,kc,ks,sum;
  kc=0.; ks=0.;
  
  for(i=0;i<VOLUME;i++){
    for(mu=0;mu<4;mu++){
      xm=&moment[i][mu];
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
double ini_momenta(const int repro) {
  
  su3adj *xm;
  int i,mu,k;
  int rlxd_state[105];
  static double y[8];
  static double tt,tr,ts,kc,ks,sum;
  
  if(repro == 1) {
    if(g_proc_id==0){
      kc=0.; 
      ks=0.;
      for(i=0;i<VOLUME;i++){ 
	for(mu=0;mu<4;mu++){
	  sum=0.;
	  xm=&moment[i][mu];
	  gauss_vector(y,8);
	  (*xm).d1=1.4142135623731*y[0];
	  (*xm).d2=1.4142135623731*y[1];
	  sum+=(*xm).d1*(*xm).d1+(*xm).d2*(*xm).d2;
	  (*xm).d3=1.4142135623731*y[2];
	  (*xm).d4=1.4142135623731*y[3];
	  sum+=(*xm).d3*(*xm).d3+(*xm).d4*(*xm).d4;
	  (*xm).d5=1.4142135623731*y[4];
	  (*xm).d6=1.4142135623731*y[5];
	  sum+=(*xm).d5*(*xm).d5+(*xm).d6*(*xm).d6;
	  (*xm).d7=1.4142135623731*y[6];
	  (*xm).d8=1.4142135623731*y[7];
	  sum+=(*xm).d7*(*xm).d7+(*xm).d8*(*xm).d8;
	  tr=sum+kc;
	  ts=tr+ks;
	  tt=ts-ks;
	  ks=ts;
	  kc=tr-tt;
	}
      }
#ifdef MPI
      /* send the state for the random-number generator to 1 */
      rlxd_get(rlxd_state);
      MPI_Send(&rlxd_state[0], 105, MPI_INT, 1, 101, MPI_COMM_WORLD);
#endif
    }
    
#ifdef MPI
    if(g_proc_id != 0){
      MPI_Recv(&rlxd_state[0], 105, MPI_INT, g_proc_id-1, 101, MPI_COMM_WORLD, &status);
      rlxd_reset(rlxd_state);
      kc=0.; ks=0.;
      for(i=0;i<VOLUME;i++){ 
	for(mu=0;mu<4;mu++){
	  sum=0.;
	  xm=&moment[i][mu];
	  gauss_vector(y,8);
	  (*xm).d1=1.4142135623731*y[0];
	  (*xm).d2=1.4142135623731*y[1];
	  sum+=(*xm).d1*(*xm).d1+(*xm).d2*(*xm).d2;
	  (*xm).d3=1.4142135623731*y[2];
	  (*xm).d4=1.4142135623731*y[3];
	  sum+=(*xm).d3*(*xm).d3+(*xm).d4*(*xm).d4;
	  (*xm).d5=1.4142135623731*y[4];
	  (*xm).d6=1.4142135623731*y[5];
	  sum+=(*xm).d5*(*xm).d5+(*xm).d6*(*xm).d6;
	  (*xm).d7=1.4142135623731*y[6];
	  (*xm).d8=1.4142135623731*y[7];
	  sum+=(*xm).d7*(*xm).d7+(*xm).d8*(*xm).d8;
	  tr=sum+kc;
	  ts=tr+ks;
	  tt=ts-ks;
	  ks=ts;
	  kc=tr-tt;
	}
      }
      /* send the state fo the random-number 
	 generator to next processor */
      
      k=g_proc_id+1; 
      if(k==g_nproc){ 
	k=0;
      }
      rlxd_get(rlxd_state);
      MPI_Send(&rlxd_state[0], 105, MPI_INT, k, 101, MPI_COMM_WORLD);
    }
#endif
    kc=0.5*(ks+kc);
    
#ifdef MPI
    if(g_proc_id == 0){
      MPI_Recv(&rlxd_state[0], 105, MPI_INT, g_nproc-1, 101, MPI_COMM_WORLD, &status);
      rlxd_reset(rlxd_state);
    }
#endif
  }
  else {
    kc=0.; 
    ks=0.;
    for(i=0;i<VOLUME;i++){ 
      for(mu=0;mu<4;mu++){
	sum=0.;
	xm=&moment[i][mu];
	gauss_vector(y,8);
	(*xm).d1=1.4142135623731*y[0];
	(*xm).d2=1.4142135623731*y[1];
	sum+=(*xm).d1*(*xm).d1+(*xm).d2*(*xm).d2;
	(*xm).d3=1.4142135623731*y[2];
	(*xm).d4=1.4142135623731*y[3];
	sum+=(*xm).d3*(*xm).d3+(*xm).d4*(*xm).d4;
	(*xm).d5=1.4142135623731*y[4];
	(*xm).d6=1.4142135623731*y[5];
	sum+=(*xm).d5*(*xm).d5+(*xm).d6*(*xm).d6;
	(*xm).d7=1.4142135623731*y[6];
	(*xm).d8=1.4142135623731*y[7];
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

static char const rcsid[] = "$Id$";

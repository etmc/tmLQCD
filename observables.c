/*******************************************************************************
*
* File observables.c
*
*
* The externally accessible functions are
*
*   double measure_gauge_action(void)
*     Returns the value of the action
*
* Author: Martin Hasenbusch <martin.hasenbusch@desy.de>
* Date: Thu, Aug 9, 2001 
*
*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "sse.h"
#include "su3.h"
#include "su3adj.h"
#include "geometry_eo.h"
#include "global.h"
#include "observables.h"

double measure_gauge_action(void)
{
  int ix,ix1,ix2,mu1,mu2;
  static su3 pr1,pr2; 
  su3 *v,*w;
  static double ga,ac,gas; 
  static double ks,kc,tr,ts,tt;
  kc=0.0; ks=0.0;
  for (ix=0;ix<VOLUME;ix++){
    for (mu1=0;mu1<3;mu1++){ 
      ix1=g_iup[ix][mu1];
      for (mu2=mu1+1;mu2<4;mu2++){ 
	ix2=g_iup[ix][mu2];
	v=&g_gauge_field[ix][mu1];
	w=&g_gauge_field[ix1][mu2];
	_su3_times_su3(pr1,*v,*w);
	v=&g_gauge_field[ix][mu2];
	w=&g_gauge_field[ix2][mu1];
	_su3_times_su3(pr2,*v,*w);
	_trace_su3_times_su3d(ac,pr1,pr2);
	tr=ac+kc;
	ts=tr+ks;
	tt=ts-ks;
	ks=ts;
	kc=tr-tt;
      }
    }
  }
  ga=(kc+ks)/3.0;
#ifdef MPI
  MPI_Allreduce(&ga, &gas, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return gas;
#else
  return ga;
#endif
}

/* $Id$ */

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
#include "ranlxd.h"
#include "sse.h"
#include "linalg_eo.h"
#include "start.h"
#include "linsolve.h"
#include "xchange.h"
#include "deriv_Sb.h"
#include "Hopping_Matrix.h"
#include "tm_operators.h"
#include "solver/solver.h"
#include "get_rectangle_staples.h"
#include "derivative_psf.h"
#include "gamma.h"
#include "get_staples.h"
#include "update_backward_gauge.h"
#include "read_input.h"
#include "stout_smear.h"
#include "stout_smear_force.h"
#include "phmc.h"
#include "derivative_nondegenerate.h"
#include "hybrid_update.h"

extern int ITER_MAX_BCG;
extern int ITER_MAX_CG;


/***********************************
 *
 * Computes the new gauge momenta
 *
 ***********************************/

void gauge_momenta(double step) 
{

  int i, j, mu;
  static su3 v, w;
  su3 *z, force;
  static su3adj deriv;
  su3adj *xm;
  static double st, st1;
  double sum=0., sum1=0., max=0., max2=0.;
  double sum2=0.;
  double atime=0., etime=0.;

  extern su3 * stout_force_field;
  extern su3 ** g_stout_force_field;
#ifdef _KOJAK_INST
#pragma pomp inst begin(gaugemomenta)
#endif
/*   int x0, x1, x2, x3; */

  st  = -step * g_rgi_C0 * g_beta/3.0;
  st1 = -step * g_rgi_C1 * g_beta/3.0;

#ifdef MPI
  atime = MPI_Wtime();
#else
  atime = (double)clock()/(double)(CLOCKS_PER_SEC);
#endif
    for(i = 0; i < VOLUME; i++)
    { 
      for(mu=0;mu<4;mu++)
      {
        z=&g_gauge_field[i][mu];
        xm=&moment[i][mu];
        v=get_staples(i,mu, g_gauge_field); 
        _su3_times_su3d(w,*z,v);
        _trace_lambda(deriv,w);
        if(g_debug_level > 0) 
        {
          sum2 = _su3adj_square_norm(deriv);
          sum+= sum2;
          if(sum2 > max) 
            max = sum2;
        }
        _minus_const_times_mom(*xm,st,deriv);

        if(g_rgi_C1 > 0. || g_rgi_C1 < 0.) 
        {
          get_rectangle_staples(&v, i, mu);
          _su3_times_su3d(w, *z, v);
          _trace_lambda(deriv, w);
          if(g_debug_level > 0) 
          {
            sum2 =_su3adj_square_norm(deriv);
            sum1+= sum2;
            if(sum2 > max2) 
              max2 = sum2;
          }
          _minus_const_times_mom(*xm, st1, deriv);
        }
      }
    }
#ifdef MPI
  etime = MPI_Wtime();
#else
  etime = (double)clock()/(double)(CLOCKS_PER_SEC);
#endif

  if(g_debug_level > 0) {
#ifdef MPI
    MPI_Reduce(&sum, &sum2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    sum = sum2;
    MPI_Reduce(&max, &sum2, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    max = sum2;
#endif
    if(g_rgi_C1 > 0. || g_rgi_C1 < 0.) {
#ifdef MPI
      MPI_Reduce(&sum1, &sum2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      sum1 = sum2;
      MPI_Reduce(&max2, &sum2, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
      max2 = sum2;
#endif
      if(g_proc_id == 0) {
	/* Here is a factor 1/4 missing            */
	printf("gaugerec %e max %e factor %e\n", sum1/((double)(VOLUME*g_nproc))/4., max2,
	       g_rgi_C1*g_rgi_C1*g_beta*g_beta/9);
	fflush(stdout);
      }
    }

    if(g_proc_id == 0) {
      /* Here is a factor 1/4 missing            */
      printf("gaugeplaq %e max %e factor %e time/s %f\n", sum/((double)(VOLUME*g_nproc))/4., max, 
	     g_rgi_C0*g_rgi_C0*g_beta*g_beta/9, etime-atime);
      fflush(stdout);
    }
  }
#ifdef _KOJAK_INST
#pragma pomp inst end(gaugemomenta)
#endif
}

/*----------------------------------------------------------------------------*/

void update_fermion_momenta(double step, const int S, 
			    const int do_all) {
  int i,mu;
  double tmp;
  su3adj *xm,*deriv;

  double sum=0., max=0.;
  double sum2=0.;
  static int co = 0;
  co++;
  double atime=0., etime=0.;

  if(!(g_running_phmc && phmc_no_flavours == 2)) {
#ifdef MPI
    atime = MPI_Wtime();
#else
    atime = (double)clock()/(double)(CLOCKS_PER_SEC);
#endif
    if(do_all == 1) {
      /* set deriv to zero here */
      derivative_psf(0, 1);
      for(i = 1; i < g_nr_of_psf; i++) {
        /* and add all the rest */
	derivative_psf(i, 0);
      }
    }
    else {
      derivative_psf(S, 1);
    }
#ifdef MPI
    xchange_deri();
#endif
    for(i = 0; i < VOLUME; i++) {
      for(mu=0;mu<4;mu++){
	xm=&moment[i][mu];
	deriv=&df0[i][mu];
	if(g_debug_level > 0) {
	  sum2 = _su3adj_square_norm(*deriv); 
	  sum+= sum2;
	  if(g_debug_level > 2 && co < 20 && g_proc_id == 0) {
	    printf("histogram%d %e\n",S,sum2);
	    fflush(stdout);
	  }
	  if(sum2 > max) max = sum2;
	}
	/* This 2* is coming from what?             */
	/* From a missing factor 2 in trace_lambda? */
	tmp = 2.*step;
	_minus_const_times_mom(*xm,tmp,*deriv); 
      }
    }
#ifdef MPI
    etime = MPI_Wtime();
#else
    etime = (double)clock()/(double)(CLOCKS_PER_SEC);
#endif
    if(g_debug_level > 0) {
#ifdef MPI
      MPI_Reduce(&sum, &sum2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      sum = sum2;
      MPI_Reduce(&max, &sum2, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
      max = sum2;
#endif
      if(g_proc_id == 0) {
	/* The factor 2 from above is missing here */
	printf("fermionforce%d %e max %e time/s %f\n", S, sum/((double)(VOLUME*g_nproc))/4., max, etime-atime);
	fflush(stdout);
      }
    }
  }

  /* the 1+1 part */
  if(g_running_phmc && (do_all || (phmc_heavy_timescale == S))) {
#ifdef MPI
    atime = MPI_Wtime();
#endif
    sum = 0.;
    max = 0.;
    derivative_nondegenerate(); 
    
#ifdef MPI
    xchange_deri();
#endif
    for(i = 0; i < VOLUME; i++) {
      for(mu=0;mu<4;mu++){
	xm=&moment[i][mu];
	deriv=&df0[i][mu];

	if(g_debug_level > 0) {
	  sum2 = _su3adj_square_norm(*deriv); 
	  sum+= sum2;
	  if(sum2 > max) max = sum2;
	}

	tmp = -2.*step*phmc_Cpol*phmc_invmaxev;
	_minus_const_times_mom(*xm,tmp,*deriv); 
      }
    }
#ifdef MPI
    etime = MPI_Wtime();
#endif
    if(g_debug_level > 0) {
#ifdef MPI
      MPI_Reduce(&sum, &sum2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      sum = sum2;
      MPI_Reduce(&max, &sum2, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
      max = sum2;
#endif
      if(g_proc_id == 0) {
	printf("fermionforceheavydoubletts%d %e max %e time/s %f\n", S, 
	       sum/((double)(VOLUME*g_nproc))/4.*phmc_Cpol*phmc_invmaxev, 
	       max*phmc_Cpol*phmc_invmaxev, etime-atime);
	fflush(stdout);
      }
    }
  }
}

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
   *
   * The backward copy of the gauge field
   * is not updated here!
   *
   */
#ifdef _KOJAK_INST
#pragma pomp inst end(updategauge)
#endif
}

/*----------------------------------------------------------------------------*/

void leap_frog(double step, int m, int nsmall) {
  int i,j;
  double smallstep;

  /* initialize the counter for the inverter */
  count00=0; count01=0; count10=0; count11=0; count20=0; count21=0;
  /* adjust the step-size to standard convention */
  /*   step*=0.7071067811865; */
  smallstep=step/nsmall;

#ifdef _GAUGE_COPY
  update_backward_gauge();
#endif
  update_fermion_momenta(0.5*step, 0, 1);
  gauge_momenta(0.5*smallstep);
  for(i=1;i<m;i++){
    for(j=0;j<nsmall;j++){
      update_gauge(smallstep); 
      gauge_momenta(smallstep);
    }
#ifdef _GAUGE_COPY
    update_backward_gauge();
#endif
    update_fermion_momenta(step, 0, 1);
  }
  for(j=1;j<nsmall;j++){
    update_gauge(smallstep); 
    gauge_momenta(smallstep);
  }
  update_gauge(smallstep); 
  gauge_momenta(0.5*smallstep);
#ifdef _GAUGE_COPY
  update_backward_gauge();
#endif
  update_fermion_momenta(0.5*step, 0, 1);
}

/*----------------------------------------------------------------------------*/

void sexton(double step, int m, int nsmall) {
  int i,j;
  /*   int ev = 10; */
  double smallstep;
  /* initialize the counter for the inverter */
  count00=0; count01=0; count10=0; count11=0; count20=0; count21=0;
  /* adjust the step-size to standard convention */
  /*   step*=0.7071067811865; */
  smallstep=step/nsmall;

#ifdef _GAUGE_COPY
  update_backward_gauge();
#endif
  update_fermion_momenta(step/6., 0, 1);
  gauge_momenta(smallstep/12.);
  for(i=1;i<m;i++){
    for(j=0;j<nsmall;j++){
      update_gauge(smallstep/4.);
      gauge_momenta(smallstep/3.);
      update_gauge(smallstep/4.);
      gauge_momenta(smallstep/6.);
    }

#ifdef _GAUGE_COPY
    update_backward_gauge();
#endif
    update_fermion_momenta(2.*step/3., 0, 1);
    for(j=0;j<nsmall;j++) {
      update_gauge(smallstep/4.);
      gauge_momenta(smallstep/3.);
      update_gauge(smallstep/4.);
      gauge_momenta(smallstep/6.);
    }
#ifdef _GAUGE_COPY
    update_backward_gauge();
#endif
    update_fermion_momenta(step/3., 0, 1);
  }
  for(j=0;j<nsmall;j++){
    update_gauge(smallstep/4.);
    gauge_momenta(smallstep/3.);
    update_gauge(smallstep/4.);
    gauge_momenta(smallstep/6.);
  }
#ifdef _GAUGE_COPY
  update_backward_gauge();
#endif
  update_fermion_momenta(2*step/3., 0, 1);
  for(j=1;j<nsmall;j++){
    update_gauge(smallstep/4.);
    gauge_momenta(smallstep/3.);
    update_gauge(smallstep/4.);
    gauge_momenta(smallstep/6.);
  }
  update_gauge(smallstep/4.);
  gauge_momenta(smallstep/3.);
  update_gauge(smallstep/4.);
  gauge_momenta(smallstep/12.);
#ifdef _GAUGE_COPY
  update_backward_gauge();
#endif
  update_fermion_momenta(step/6., 0, 1);
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
double ini_momenta() {
  
  su3adj *xm;
  int i,mu,k;
  int rlxd_state[105];
  static double y[8];
  static double tt,tr,ts,kc,ks,sum;

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

  MPI_Allreduce(&kc, &ks, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return ks;
#else
  return kc;
#endif
}

static char const rcsid[] = "$Id$";

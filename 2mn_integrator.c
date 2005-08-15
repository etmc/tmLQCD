/* $Id$ */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "su3.h"
#include "su3adj.h"
#include "expo.h"
#include "ranlxd.h"
#include "sse.h"
#include "global.h"
#include "linalg_eo.h"
#include "clover_eo.h"
#include "start.h"
#include "sw.h"
#include "linsolve.h"
#include "xchange.h"
#include "deriv_Sb.h"
#include "Hopping_Matrix.h"
#include "tm_operators.h"
#include "hybrid_update.h"
#include "derivative_psf.h"
#include "2mn_integrator.h"

extern int count00, count01, count10, count11, count20, count21;


/* double dtmp[5]; */

void mn2_integrator(int * const n_int, const double tau, 
		    const int S, const int halfstep, double * const lambda) {
  int i,j;
  double eps, eps0,
    oneminus2lambda = (1.-2.*lambda[S]), oneminus2lambda0;
  
  if(S == g_nr_of_psf) {
    /* initialize the counter for the inverter */
    count00=0; count01=0; count10=0; count11=0; count20=0; count21=0;
/*     for(i=0;i<5;i++) { */
/*       dtmp[i] = 0.; */
/*     } */

    eps = tau/((double)n_int[S]);

#ifdef _GAUGE_COPY
    update_backward_gauge();
#endif
    for(i = 1; i < S; i++) {
      update_fermion_momenta(lambda[S-i+1]*eps, S-i);
/*       dtmp[S-i+1] += lambda[S-i+1]*eps; */
      eps /= ((double)n_int[S-i])*2;
    }
    update_fermion_momenta(lambda[1]*eps, 0);
/*     dtmp[1] += lambda[1]*eps; */
    eps /= ((double)n_int[0])*2;
    gauge_momenta(lambda[0]*eps);
/*     dtmp[0] += lambda[0]*eps; */
    
  }
  
  eps = tau/((double)n_int[S]);
  if(S == 1) {
    eps0 = 0.5*eps/((double)n_int[0]);
    oneminus2lambda0 = (1.-2.*lambda[0]);

    for(i = 1; i < n_int[S]; i++) {
      for(j = 0; j < n_int[0]; j++) {
	update_gauge(0.5*eps0);
	gauge_momenta(oneminus2lambda0*eps0);
/* 	dtmp[0] += oneminus2lambda0*eps0; */
	update_gauge(0.5*eps0);
	gauge_momenta(2*lambda[0]*eps0);
/* 	dtmp[0] += 2*lambda[0]*eps0; */
      }
#ifdef _GAUGE_COPY
      update_backward_gauge();
#endif
      update_fermion_momenta(oneminus2lambda*eps, S-1);
/*       dtmp[S] += oneminus2lambda*eps; */
      for(j = 0; j < n_int[0]; j++) {
	update_gauge(0.5*eps0);
	gauge_momenta(oneminus2lambda0*eps0);
/* 	dtmp[0] += oneminus2lambda0*eps0; */
	update_gauge(0.5*eps0);
	gauge_momenta(2*lambda[0]*eps0);
/* 	dtmp[0] += 2*lambda[0]*eps0; */
      }
#ifdef _GAUGE_COPY
      update_backward_gauge();
#endif
      update_fermion_momenta(2*lambda[S]*eps, S-1);
/*       dtmp[S] += 2*lambda[S]*eps; */
    }
    for(j = 0; j < n_int[0]; j++) {
      update_gauge(0.5*eps0);
      gauge_momenta(oneminus2lambda0*eps0);
/*       dtmp[0] += oneminus2lambda0*eps0; */
      update_gauge(0.5*eps0);
      gauge_momenta(2*lambda[0]*eps0);
/*       dtmp[0] += 2*lambda[0]*eps0; */
    }
#ifdef _GAUGE_COPY
    update_backward_gauge();
#endif
    update_fermion_momenta(oneminus2lambda*eps, S-1);
/*     dtmp[S] += oneminus2lambda*eps; */
    for(j = 1; j < n_int[0]; j++) {
      update_gauge(0.5*eps0);
      gauge_momenta(oneminus2lambda0*eps0);
/*       dtmp[0] += oneminus2lambda0*eps0; */
      update_gauge(0.5*eps0);
      gauge_momenta(2*lambda[0]*eps0);
/*       dtmp[0] += 2*lambda[0]*eps0; */
    }
    update_gauge(0.5*eps0);
    gauge_momenta(oneminus2lambda0*eps0);
/*     dtmp[0] += oneminus2lambda0*eps0; */
    update_gauge(0.5*eps0);
#ifdef _GAUGE_COPY
    update_backward_gauge();
#endif
    if(halfstep != 1) {
      gauge_momenta(2.*lambda[0]*eps0);
/*       dtmp[0] += 2.*lambda[0]*eps0; */
      update_fermion_momenta(2.*lambda[S]*eps, S-1);
/*       dtmp[S] += 2.*lambda[S]*eps; */
    }
  }
  else {
    for(i = 1; i < n_int[S]; i++){
      mn2_integrator(n_int, eps/2., S-1, 0, lambda);
      update_fermion_momenta(oneminus2lambda*eps, S-1);
      mn2_integrator(n_int, eps/2., S-1, 0, lambda);
      update_fermion_momenta(2*lambda[S]*eps, S-1);
    }
    mn2_integrator(n_int, eps/2., S-1, 0, lambda);
    update_fermion_momenta(oneminus2lambda*eps, S-1);
    if(S == g_nr_of_psf) {
      mn2_integrator(n_int, eps/2., S-1, 1, lambda);
    }
    else mn2_integrator(n_int, eps/2., S-1, halfstep, lambda);
    if(halfstep != 1 && S != g_nr_of_psf) {
      update_fermion_momenta(2*lambda[S]*eps, S-1);
    }
  }

  if(S == g_nr_of_psf) {
    for(i = 1; i < S; i++) { 
      update_fermion_momenta(lambda[S-i+1]*eps, S-i); 
      eps /= ((double)n_int[S-i])*2; 
    } 
    update_fermion_momenta(lambda[1]*eps, 0); 
    eps /= ((double)n_int[0])*2; 
    gauge_momenta(lambda[0]*eps); 

/*     for(i=0;i<5;i++) { */
/*       printf("%d %e %e\n", i, dtmp[i], lambda[i]); */
/*     } */

  }
}


void mn2p_integrator(int * const n_int, const double tau, 
		    const int S, double * const lambda) {
  int i,j;
  double eps = tau/((double)n_int[S]);
  double oneminus2lambda = (1.-2.*lambda[S]);
  
  if(S == g_nr_of_psf) {
    /* initialize the counter for the inverter */
    count00=0; count01=0; count10=0; count11=0; count20=0; count21=0;
/*     for(i=0;i<5;i++) { */
/*       dtmp[i] = 0.; */
/*     } */
  }
  
  if(S == 0) {
    for(i = 0; i < n_int[0]; i++) {
      update_gauge(lambda[0]*eps);
      gauge_momenta(0.5*eps);
/*       dtmp[0] += 0.5*eps; */
      update_gauge(oneminus2lambda*eps);
      gauge_momenta(0.5*eps);
/*       dtmp[0] += 0.5*eps; */
      update_gauge(lambda[0]*eps);
    }
#ifdef _GAUGE_COPY
    update_backward_gauge();
#endif
  }
  else {
    for(i = 0; i < n_int[S]; i++) {
      mn2p_integrator(n_int, lambda[S]*eps, S-1, lambda);
      update_fermion_momenta(0.5*eps, S-1); 
/*       dtmp[S] += 0.5*eps; */
      mn2p_integrator(n_int, oneminus2lambda*eps, S-1, lambda);
      update_fermion_momenta(0.5*eps, S-1);
/*       dtmp[S] += 0.5*eps; */
      mn2p_integrator(n_int, lambda[S]*eps, S-1, lambda);
    }
  }

/*   if(S == g_nr_of_psf) { */
/*     for(i=0;i<5;i++) { */
/*       printf("%d %e %e\n", i, dtmp[i], lambda[i]); */
/*     } */
/*   } */
}

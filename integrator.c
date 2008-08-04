/* $Id$ */

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
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
#include "start.h"
#include "linsolve.h"
#include "xchange.h"
#include "deriv_Sb.h"
#include "Hopping_Matrix.h"
#include "tm_operators.h"
#include "hybrid_update.h"
#include "update_backward_gauge.h"
#include "monomial.h"
#include "update_momenta.h"
#include "integrator.h"

extern int count00, count01, count10, count11, count20, count21;

integrator Integrator;

/* double dtmp[5]; */

int init_integrator() {
  int i, ts;
  for(i = 0; i < 10; i++) {
    Integrator.no_mnls_per_ts[i] = 0;
  }

  for(i = 0; i < no_monomials; i++) {
    ts = monomial_list[i].timescale;
    if(ts < Integrator.no_timescales) {
      Integrator.mnls_per_ts[ ts ][ Integrator.no_mnls_per_ts[ts] ] = monomial_list[i].id;
      Integrator.no_mnls_per_ts[ ts ]++;
    }
    else {
      if(g_proc_id == 0) {
	fprintf(stderr, "Warning: monomial %d is not on a valid timescale and discarded\n", i);
      }
    }
  }
  for(i = 0; i < Integrator.no_timescales; i++) {
    if(Integrator.no_mnls_per_ts[ ts ] < 1) {
      fprintf(stderr, "Error, no monomial on timescale %d!\nAborting...\n", i);
      exit(-1);
    }
  }
  return(0);
}

void integrate_2mn(integrator * itgr, const double tau, 
		   const int S, const int halfstep);
void integrate_2mnp(integrator *itgr, const double tau, const int S);
void integrate_leap_frog(integrator * itgr, const double tau, const int S, const int halfstep);

void integrate_md(integrator *itgr, const int forward) {
  int halfstep = 0;
  double tau = itgr->tau;
  if(!forward) tau = -tau;
  if(itgr->no_timescales == 2) halfstep = 1;

  if(itgr->type == 5) {
    integrate_leap_frog(itgr, tau, itgr->no_timescales-1, halfstep);
  }
  else if(itgr->type == 6) {
    integrate_2mn(itgr, tau, itgr->no_timescales-1, halfstep);
  }
  else if(itgr->type == 7) {
    integrate_2mnp(itgr, tau, itgr->no_timescales-1);
  }
  else {
    fprintf(stderr, "Unkown integrator type...!\nAborting...!\n");
    exit(1);
  }
}

void integrate_2mn(integrator * itgr, const double tau, 
		    const int S, const int halfstep) {
  int i,j=0;
  double eps, eps0,
    oneminus2lambda = (1.-2.*itgr->lambda[S]), oneminus2lambda0;

  if(S == itgr->no_timescales-1) {

    eps = tau/((double)itgr->n_int[S]);
#ifdef _GAUGE_COPY
    update_backward_gauge();
#endif
    for(i = S; i > 0; i--) {
      update_momenta(itgr->mnls_per_ts[i], itgr->lambda[i]*eps, itgr->no_mnls_per_ts[i]);
      eps /= ((double)itgr->n_int[i-1])*2;
    }
    update_momenta(itgr->mnls_per_ts[0], itgr->lambda[0]*eps, itgr->no_mnls_per_ts[0]);
  }
  
  eps = tau/((double)itgr->n_int[S]);
  if(S == 1) {
    eps0 = 0.5*eps/((double)itgr->n_int[0]);
    oneminus2lambda0 = (1.-2.*itgr->lambda[0]);

    for(i = 1; i < itgr->n_int[S]; i++) {
      for(j = 0; j < itgr->n_int[0]; j++) {
	update_gauge(0.5*eps0);
	update_momenta(itgr->mnls_per_ts[0], oneminus2lambda0*eps0, itgr->no_mnls_per_ts[0]);
	update_gauge(0.5*eps0);
	update_momenta(itgr->mnls_per_ts[0], 2.*itgr->lambda[0]*eps0, itgr->no_mnls_per_ts[0]);
      }
#ifdef _GAUGE_COPY
      update_backward_gauge();
#endif
      update_momenta(itgr->mnls_per_ts[S], oneminus2lambda*eps, itgr->no_mnls_per_ts[S]);
      for(j = 0; j < itgr->n_int[0]; j++) {
	update_gauge(0.5*eps0);
	update_momenta(itgr->mnls_per_ts[0], oneminus2lambda0*eps0, itgr->no_mnls_per_ts[0]);
	update_gauge(0.5*eps0);
	update_momenta(itgr->mnls_per_ts[0], 2.*itgr->lambda[0]*eps0, itgr->no_mnls_per_ts[0]);
      }
#ifdef _GAUGE_COPY
      update_backward_gauge();
#endif
      update_momenta(itgr->mnls_per_ts[S], 2*itgr->lambda[S]*eps, itgr->no_mnls_per_ts[S]);
    }
    for(j = 0; j < itgr->n_int[0]; j++) {
      update_gauge(0.5*eps0);
      update_momenta(itgr->mnls_per_ts[0], oneminus2lambda0*eps0, itgr->no_mnls_per_ts[0]);
      update_gauge(0.5*eps0);
      update_momenta(itgr->mnls_per_ts[0], 2.*itgr->lambda[0]*eps0, itgr->no_mnls_per_ts[0]);
    }
#ifdef _GAUGE_COPY
    update_backward_gauge();
#endif
    update_momenta(itgr->mnls_per_ts[S], oneminus2lambda*eps, itgr->no_mnls_per_ts[S]);
    for(j = 1; j < itgr->n_int[0]; j++) {
      update_gauge(0.5*eps0);
      update_momenta(itgr->mnls_per_ts[0], oneminus2lambda0*eps0, itgr->no_mnls_per_ts[0]);
      update_gauge(0.5*eps0);
      update_momenta(itgr->mnls_per_ts[0], 2.*itgr->lambda[0]*eps0, itgr->no_mnls_per_ts[0]);
    }
    update_gauge(0.5*eps0);
    update_momenta(itgr->mnls_per_ts[0], oneminus2lambda0*eps0, itgr->no_mnls_per_ts[0]);
    update_gauge(0.5*eps0);
#ifdef _GAUGE_COPY
    update_backward_gauge();
#endif
    if(halfstep != 1) {
      update_momenta(itgr->mnls_per_ts[0], 2*itgr->lambda[0]*eps0, itgr->no_mnls_per_ts[0]);
      update_momenta(itgr->mnls_per_ts[S], 2*itgr->lambda[S]*eps, itgr->no_mnls_per_ts[S]);
    }
  }
  else {
    for(i = 1; i < itgr->n_int[S]; i++){
      integrate_2mn(itgr, eps/2., S-1, 0);
      update_momenta(itgr->mnls_per_ts[S], oneminus2lambda*eps, itgr->no_mnls_per_ts[S]);
      integrate_2mn(itgr, eps/2., S-1, 0);
      update_momenta(itgr->mnls_per_ts[S], 2*itgr->lambda[S]*eps, itgr->no_mnls_per_ts[S]);
    }
    integrate_2mn(itgr, eps/2., S-1, 0);
    update_momenta(itgr->mnls_per_ts[S], oneminus2lambda*eps, itgr->no_mnls_per_ts[S]);
    if(S == itgr->no_timescales-1) {
      integrate_2mn(itgr, eps/2., S-1, 1);
    }
    else integrate_2mn(itgr, eps/2., S-1, halfstep);
    if(halfstep != 1 && S != itgr->no_timescales-1) {
      update_momenta(itgr->mnls_per_ts[S], 2*itgr->lambda[S]*eps, itgr->no_mnls_per_ts[S]);
    }
  }

  if(S == itgr->no_timescales-1) {
    for(i = S; i > 0; i--) {
      update_momenta(itgr->mnls_per_ts[i], itgr->lambda[i]*eps, itgr->no_mnls_per_ts[i]);
      eps /= ((double)itgr->n_int[i-1])*2;
    }
    update_momenta(itgr->mnls_per_ts[0], itgr->lambda[0]*eps, itgr->no_mnls_per_ts[0]);
  }
}

void integrate_2mnp(integrator *itgr, const double tau, const int S) {
  int i;
  double eps = tau/((double)itgr->n_int[S]);
  double oneminus2lambda = (1.-2.*itgr->lambda[S]);
  
  if(S == 0) {
    update_gauge(itgr->lambda[0]*eps);
    for(i = 1; i < itgr->n_int[0]; i++) {
      update_momenta(itgr->mnls_per_ts[0], 0.5*eps, itgr->no_mnls_per_ts[0]);
/*       dtmp[0] += 0.5*eps; */
      update_gauge(oneminus2lambda*eps);
      update_momenta(itgr->mnls_per_ts[0], 0.5*eps, itgr->no_mnls_per_ts[0]);
/*       dtmp[0] += 0.5*eps; */
      update_gauge(2*itgr->lambda[0]*eps);
    }
    update_momenta(itgr->mnls_per_ts[0], 0.5*eps, itgr->no_mnls_per_ts[0]);
/*       dtmp[0] += 0.5*eps; */
    update_gauge(oneminus2lambda*eps);
    update_momenta(itgr->mnls_per_ts[0], 0.5*eps, itgr->no_mnls_per_ts[0]);
/*       dtmp[0] += 0.5*eps; */
    update_gauge(itgr->lambda[0]*eps);
#ifdef _GAUGE_COPY
    update_backward_gauge();
#endif
  }
  else {
    for(i = 0; i < itgr->n_int[S]; i++) {
      integrate_2mnp(itgr, itgr->lambda[S]*eps, S-1);
      update_momenta(itgr->mnls_per_ts[S], 0.5*eps, itgr->no_mnls_per_ts[S]);

      integrate_2mnp(itgr, oneminus2lambda*eps, S-1);
      update_momenta(itgr->mnls_per_ts[S], 0.5*eps, itgr->no_mnls_per_ts[S]);

      integrate_2mnp(itgr, itgr->lambda[S]*eps, S-1);
    }
  }
}


void integrate_leap_frog(integrator * itgr, const double tau, const int S, const int halfstep) {
  int i,j;
  double eps, eps0;

  if(S == itgr->no_timescales-1) {
    /* initialize the counter for the inverter */
    eps = tau/((double)itgr->n_int[S]);

#ifdef _GAUGE_COPY
    update_backward_gauge();
#endif
    for(i = S; i > 0; i--) {
      update_momenta(itgr->mnls_per_ts[S], 0.5*eps, itgr->no_mnls_per_ts[S]);
      eps /= ((double)itgr->n_int[i-1]);
    }
    update_momenta(itgr->mnls_per_ts[0], 0.5*eps, itgr->no_mnls_per_ts[0]);
  }

  eps = tau/((double)itgr->n_int[S]);
  if(S == 1) {
    eps0 = eps/((double)itgr->n_int[0]);
    for(j = 1; j < itgr->n_int[S]; j++) {
      for(i = 0; i < itgr->n_int[0]; i++) {
	update_gauge(eps0); 
	update_momenta(itgr->mnls_per_ts[0], eps0, itgr->no_mnls_per_ts[0]);
      }
#ifdef _GAUGE_COPY
      update_backward_gauge();
#endif
      update_momenta(itgr->mnls_per_ts[S], eps, itgr->no_mnls_per_ts[S]);
    }
    for(i = 1; i < itgr->n_int[0]; i++) {
      update_gauge(eps0); 
      update_momenta(itgr->mnls_per_ts[0], eps0, itgr->no_mnls_per_ts[0]);
    }
    update_gauge(eps0); 
#ifdef _GAUGE_COPY
    update_backward_gauge();
#endif
    if(halfstep != 1) {
      gauge_momenta(eps0);
      update_momenta(itgr->mnls_per_ts[S], eps, itgr->no_mnls_per_ts[S]);
    }
  }
  else {
    for(i = 1; i < itgr->n_int[S]; i++){
      integrate_leap_frog(itgr, eps, S-1, 0);
      update_momenta(itgr->mnls_per_ts[S], eps, itgr->no_mnls_per_ts[S]);
    }
    if(S == itgr->no_timescales-1) {
      integrate_leap_frog(itgr, eps, S-1, 1);
    }
    else integrate_leap_frog(itgr, eps, S-1, halfstep);
    if(halfstep != 1 && S != itgr->no_timescales-1) {
      update_momenta(itgr->mnls_per_ts[S], eps, itgr->no_mnls_per_ts[S]);
    }
  }

  if(S == itgr->no_timescales-1) {
    for(i = S; i > 0; i--) {
      update_momenta(itgr->mnls_per_ts[S], 0.5*eps, itgr->no_mnls_per_ts[S]);
      eps /= ((double)itgr->n_int[i-1]);
    }
    update_momenta(itgr->mnls_per_ts[0], 0.5*eps, itgr->no_mnls_per_ts[0]);
  }
}

/* $Id$ */

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "global.h"
#include "su3.h"
#include "su3adj.h"
#include "ranlxd.h"
#include "sse.h"
#include "start.h"
#include "get_rectangle_staples.h"
#include "gamma.h"
#include "get_staples.h"
#include "read_input.h"
#include "observables.h"
#include "measure_rectangles.h"
#include "monomial.h"
#include "gauge_monomial.h"


void gauge_derivative(const int id) {

  int i, mu;
  static su3 v, w;
  su3 *z;
  su3adj *xm;
  monomial * mnl = &monomial_list[id];

  if(mnl->use_rectangles) {
    mnl->forcefactor = -mnl->c0 * g_beta/3.0/2.;
  }
  else {
    mnl->forcefactor = -1. * g_beta/3.0/2.;
  }

  for(i = 0; i < VOLUME; i++) { 
    for(mu=0;mu<4;mu++) {
      z=&g_gauge_field[i][mu];
      xm=&df0[i][mu];
      v=get_staples(i,mu, g_gauge_field); 
      _su3_times_su3d(w,*z,v);
      _add_trace_lambda((*xm),w);

      if(mnl->use_rectangles) {
	get_rectangle_staples(&v, i, mu);
	_su3_times_su3d(w, *z, v);
 	_mul_add_trace_lambda((*xm), w, mnl->c1/mnl->c0);
      }
    }
  }
  return;
}

void gauge_heatbath(const int id) {
  monomial * mnl = &monomial_list[id];
  if(mnl->use_rectangles) mnl->c0 = 1. - 8.*mnl->c1;

  mnl->energy0 = mnl->c0 * measure_gauge_action();
  if(mnl->use_rectangles) {
    mnl->energy0 += mnl->c1 * measure_rectangles();
  }
  if(g_proc_id == 0 && g_debug_level > 3) {
    printf("called gauge_heatbath for id %d %d\n", id, mnl->even_odd_flag);
  }
}

double gauge_acc(const int id) {
  monomial * mnl = &monomial_list[id];

  mnl->energy1 = mnl->c0 * measure_gauge_action();
  if(mnl->use_rectangles) {
    mnl->energy1 += mnl->c1 * measure_rectangles();
  }
  if(g_proc_id == 0 && g_debug_level > 3) {
    printf("called gauge_acc for id %d %d\n", id, mnl->even_odd_flag);
  }
  return(g_beta*(mnl->energy0 - mnl->energy1));
}

/***********************************************************************
 * $Id: sf_gauge_monomial.c,v 1.6 2009/03/27
 *
 * Jenifer Gonzalez Lopez
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
#include "sf_gauge_monomial.h"
#include "sf_utils.h"

void sf_gauge_derivative(const int id) {

  int i, mu;
  static su3 v, w;
  su3 *z;
  su3adj *xm;
  monomial * mnl = &monomial_list[id];

  printf ("hola");

  if(mnl->use_rectangles) {
    mnl->forcefactor = -mnl->c0 * g_beta/3.0;
  }
  else {
    mnl->forcefactor = -1. * g_beta/3.0;
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

void sf_gauge_heatbath( const int id )
{
  monomial* mnl = &(monomial_list[id]);

  if( mnl->use_rectangles ){ mnl->c0 = 1. - 8.*mnl->c1; }

  mnl->energy0 = g_beta * ( mnl->c0 * measure_gauge_action() );

  if(mnl->use_rectangles) {
    mnl->energy0 += g_beta*(mnl->c1 * measure_rectangles());
  }
  if(g_proc_id == 0 && g_debug_level > 3) {
    printf("called gauge_heatbath for id %d %d\n", id, mnl->even_odd_flag);
  }
}

double sf_gauge_acc( const int id )
{
  monomial* mnl = &(monomial_list[id]);
  double sq_plaq = 0;
  double sq_bulk_plaq = 0;
  double sq_boundary_space_space_plaq = 0;
  double sq_boundary_space_time_plaq = 0;
  double sq_wrapped_plaq = 0;

  double rect_plaq = 0;

  sq_plaq = calc_sq_plaq();
  sq_bulk_plaq = calc_bulk_sq_plaq();
  sq_boundary_space_space_plaq = calc_boundary_space_space_sq_plaq();
  sq_boundary_space_time_plaq = calc_boundary_space_time_sq_plaq();
  sq_wrapped_plaq = calc_wrapped_sq_plaq();

  rect_plaq = calc_rect_plaq();

  #if 1
  {
    fprintf( stderr, "sq_plaq = %e\n", sq_plaq );
    fprintf( stderr, "beta * c0 * sq_plaq = %e\n", g_beta * mnl->c0 * sq_plaq );

    fprintf( stderr, "sq_bulk_plaq = %e\n", sq_bulk_plaq );
    fprintf( stderr, "beta * c0 * sq_bulk_plaq = %e\n", g_beta * mnl->c0 * sq_bulk_plaq );

    fprintf( stderr, "sq_wrapped_plaq = %e\n", sq_wrapped_plaq );
    fprintf( stderr, "beta * c0 * sq_wrapped_plaq = %e\n", g_beta * mnl->c0 * sq_wrapped_plaq );

    fprintf( stderr, "rect_plaq = %e\n", rect_plaq );
    fprintf( stderr, "beta * c1 * rect_plaq = %e\n", g_beta * mnl->c1 * rect_plaq );

    fprintf( stderr, "bulk + bound(ss) + bound(st) + wrapped = %e + %e + %e + %e = %e =?= %e = total\n",
             sq_bulk_plaq, sq_boundary_space_space_plaq, sq_boundary_space_time_plaq, sq_wrapped_plaq,
             sq_bulk_plaq + sq_boundary_space_space_plaq + sq_boundary_space_time_plaq + sq_wrapped_plaq, sq_plaq );

    fprintf( stderr, "my energy =    %e\n", g_beta * ( mnl->c0 * sq_plaq + mnl->c1 * rect_plaq ) );
  }
  #endif

  /*mnl->energy1 = g_beta*( mnl->c0 * measure_gauge_action() );*/

  /* The bulk contribution is the same. */
  mnl->energy1  = g_beta * mnl->c0 * sq_bulk_plaq;

  /* The space-time boundary contribution must be weighted differently. */
  fprintf( stderr, "mnl->ct = %e\n", mnl->ct );
  mnl->energy1 += g_beta * mnl->c0 * mnl->ct * sq_boundary_space_time_plaq;  

  /* The space-space boundary contribution must be weighted differently. */
  fprintf( stderr, "mnl->cs = %e\n", mnl->cs );
  mnl->energy1 += g_beta * mnl->c0 * mnl->cs * sq_boundary_space_space_plaq;  

  /* Include the missing plaquettes if requested. */
  if( g_sf_inc_wrap_sq == 1 ){ mnl->energy1 += g_beta * mnl->c0 * sq_wrapped_plaq; }

  if( mnl->use_rectangles )
  {
    mnl->energy1 += g_beta*( mnl->c1 * measure_rectangles() );
  }
  fprintf( stderr, "mnl->energy1 = %e\n", mnl->energy1 );

  if( ( g_proc_id == 0 ) & ( g_debug_level > 3 ) )
  {
    printf( "called sf_gauge_acc for id %d %d dH = %1.4e\n", 
	    id, mnl->even_odd_flag, mnl->energy0 - mnl->energy1 );
  }

  return ( mnl->energy0 - mnl->energy1 );
}

/***********************************************************************
 *
 * Copyright (C) 2001 Martin Hasebusch
 *
 * some changes by C. Urbach 2002-2008,2012
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
 *
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
#include "hamiltonian_field.h"
#include "update_gauge.h"


/*******************************************************
 *
 * Updates the gauge field corresponding to the momenta
 *
 *******************************************************/


void update_gauge(const double step, hamiltonian_field_t * const hf) {

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
      /* moment[i][mu] = h_{i,mu}^{alpha} */
      xm = &hf->momenta[i][mu];
      z = &hf->gaugefield[i][mu];
      _assign_const_times_mom(deriv, step, *xm);
      v = restoresu3( exposu3(deriv) );
      _su3_times_su3(w, v, *z);
      _su3_assign(*z, w);
    }
  }
  
#ifdef MPI
  /* for parallelization */
  xchange_gauge(hf->gaugefield);
#endif
  /*
   * The backward copy of the gauge field
   * is not updated here!
   */
  hf->update_gauge_copy = 1;
  g_update_gauge_copy = 1;
  hf->update_gauge_energy = 1;
  g_update_gauge_energy = 1;
  hf->update_rectangle_energy = 1;
  g_update_rectangle_energy = 1;
  return;
#ifdef _KOJAK_INST
#pragma pomp inst end(updategauge)
#endif
}

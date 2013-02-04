/***********************************************************************
 *
 * Copyright (C) 2003 Mauro Papinutto
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
 * perform a random gauge transformation
 *
 *
 *******************************************************************************/

#if HAVE_CONFIG_H
#include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include "global.h"
#include "su3.h"
#include "start.h"
#include "rnd_gauge_trafo.h"

void rnd_gauge_trafo(const int repro, su3 ** const gf){
  int ix,iy,mu;
  static su3 u,v,w,x,y;
  su3 * _gauge_trafo = NULL;
  su3 * gauge_trafo = NULL;

  if((_gauge_trafo = calloc(VOLUMEPLUSRAND+1, sizeof(su3))) == NULL) {
    fprintf(stderr, "Could not allocate memory in rnd_gauge_trafo. Exiting!\n");
    exit(0);
  }
  gauge_trafo = (su3*)(((unsigned long int)(gauge_trafo)+ALIGN_BASE)&~ALIGN_BASE);

  random_gauge_field(repro, gauge_trafo);

#ifdef MPI
  xchange_gauge(gauge_trafo);
#endif

  for (ix=0;ix<VOLUME;ix++){

    u=gauge_trafo[ix];

    for (mu=0;mu<4;mu++){
      iy=g_iup[ix][mu];
      w=gauge_trafo[iy];
      _su3_dagger(v,w);
      w=g_gauge_field[ix][mu];

      _su3_times_su3(x,w,v);
      _su3_times_su3(y,u,x);

      gf[ix][mu]=y;
    }
  }
      
  free(_gauge_trafo);
}


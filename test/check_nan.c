/***********************************************************************
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008 Carsten Urbach
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
#include "global.h"
#include "su3adj.h"
#include "check_nan.h"

int check_nan() {
  int i, mu;
  su3 * um;
  
  um = &g_gauge_field[0][0];
  for(i = 0; i < VOLUMEPLUSRAND; i++) {
    for(mu = 0; mu < 4; mu++) {
      if(__isnan(creal(um->c00))|| __isnan(cimag(um->c00)) || __isnan(creal(um->c01)) || __isnan(cimag(um->c01)) ||
	 __isnan(creal(um->c02)) || __isnan(cimag(um->c02)) || __isnan(creal(um->c10)) || __isnan(cimag(um->c10)) ||
	 __isnan(creal(um->c11)) || __isnan(cimag(um->c11)) || __isnan(creal(um->c12)) || __isnan(cimag(um->c12)) ||
	 __isnan(creal(um->c20)) || __isnan(cimag(um->c20)) || __isnan(creal(um->c21)) || __isnan(cimag(um->c21)) ||
	 __isnan(creal(um->c22)) || __isnan(cimag(um->c22))) {
	return(i);
      }
      um++;
    }
  }
  return(-1);
}

int check_greater(const double a) {
  int i, mu;
  su3 * um;
  
  um = &g_gauge_field[0][0];
  for(i = 0; i < VOLUMEPLUSRAND; i++) {
    for(mu = 0; mu < 4; mu++) {
      if((creal(um->c00) > a)|| (cimag(um->c00) > a) || (creal(um->c01) > a) || (cimag(um->c01) > a) ||
	 (creal(um->c02) > a) || (cimag(um->c02) > a) || (creal(um->c10) > a) || (cimag(um->c10) > a) ||
	 (creal(um->c11) > a) || (cimag(um->c11) > a) || (creal(um->c12) > a) || (cimag(um->c12) > a) ||
	 (creal(um->c20) > a) || (cimag(um->c20) > a) || (creal(um->c21) > a) || (cimag(um->c21) > a) ||
	 (creal(um->c22) > a) || (cimag(um->c22) > a)) {
	return(i);
      }
      um++;
    }
  }
  return(-1);
}

int check_nan_gauge(const int i, const int mu) {
  su3 * um;
  
  um = &g_gauge_field[i][mu];
  if(__isnan(creal(um->c00))|| __isnan(cimag(um->c00)) || __isnan(creal(um->c01)) || __isnan(cimag(um->c01)) ||
     __isnan(creal(um->c02)) || __isnan(cimag(um->c02)) || __isnan(creal(um->c10)) || __isnan(cimag(um->c10)) ||
     __isnan(creal(um->c11)) || __isnan(cimag(um->c11)) || __isnan(creal(um->c12)) || __isnan(cimag(um->c12)) ||
     __isnan(creal(um->c20)) || __isnan(cimag(um->c20)) || __isnan(creal(um->c21)) || __isnan(cimag(um->c21)) ||
     __isnan(creal(um->c22)) || __isnan(cimag(um->c22))) {
    return(i);
  }
  return(-1);
}

int check_su3adj(su3adj * s, const double a) {

  if(s->d1 > a || s->d2 > a || s->d3 > a || s->d4 > a || 
     s->d5 > a || s->d6 > a || s->d7 > a || s->d8 > a ) {
    return(1);
  }
  return(0);
}


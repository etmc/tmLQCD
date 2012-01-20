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
      if(__isnan((*um).c00.re)|| __isnan((*um).c00.im) || __isnan((*um).c01.re) || __isnan((*um).c01.im) ||
	 __isnan((*um).c02.re) || __isnan((*um).c02.im) || __isnan((*um).c10.re) || __isnan((*um).c10.im) ||
	 __isnan((*um).c11.re) || __isnan((*um).c11.im) || __isnan((*um).c12.re) || __isnan((*um).c12.im) ||
	 __isnan((*um).c20.re) || __isnan((*um).c20.im) || __isnan((*um).c21.re) || __isnan((*um).c21.im) ||
	 __isnan((*um).c22.re) || __isnan((*um).c22.im)) {
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
      if(((*um).c00.re > a)|| ((*um).c00.im > a) || ((*um).c01.re > a) || ((*um).c01.im > a) ||
	 ((*um).c02.re > a) || ((*um).c02.im > a) || ((*um).c10.re > a) || ((*um).c10.im > a) ||
	 ((*um).c11.re > a) || ((*um).c11.im > a) || ((*um).c12.re > a) || ((*um).c12.im > a) ||
	 ((*um).c20.re > a) || ((*um).c20.im > a) || ((*um).c21.re > a) || ((*um).c21.im > a) ||
	 ((*um).c22.re > a) || ((*um).c22.im > a)) {
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
  if(__isnan((*um).c00.re)|| __isnan((*um).c00.im) || __isnan((*um).c01.re) || __isnan((*um).c01.im) ||
     __isnan((*um).c02.re) || __isnan((*um).c02.im) || __isnan((*um).c10.re) || __isnan((*um).c10.im) ||
     __isnan((*um).c11.re) || __isnan((*um).c11.im) || __isnan((*um).c12.re) || __isnan((*um).c12.im) ||
     __isnan((*um).c20.re) || __isnan((*um).c20.im) || __isnan((*um).c21.re) || __isnan((*um).c21.im) ||
     __isnan((*um).c22.re) || __isnan((*um).c22.im)) {
    return(i);
  }
  return(-1);
}

int check_su3adj(su3adj * s, const double a) {

  if((*s).d1 > a || (*s).d2 > a || (*s).d3 > a || (*s).d4 > a || 
     (*s).d5 > a || (*s).d6 > a || (*s).d7 > a || (*s).d8 > a ) {
    return(1);
  }
  return(0);
}


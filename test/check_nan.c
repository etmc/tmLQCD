/* $Id$ */

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
      if(__isnan((*um).c11.re)|| __isnan((*um).c11.im) || __isnan((*um).c12.re) || __isnan((*um).c12.im) ||
	 __isnan((*um).c13.re) || __isnan((*um).c13.im) || __isnan((*um).c21.re) || __isnan((*um).c21.im) ||
	 __isnan((*um).c22.re) || __isnan((*um).c22.im) || __isnan((*um).c23.re) || __isnan((*um).c23.im) ||
	 __isnan((*um).c31.re) || __isnan((*um).c31.im) || __isnan((*um).c32.re) || __isnan((*um).c32.im) ||
	 __isnan((*um).c33.re) || __isnan((*um).c33.im)) {
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
      if(((*um).c11.re > a)|| ((*um).c11.im > a) || ((*um).c12.re > a) || ((*um).c12.im > a) ||
	 ((*um).c13.re > a) || ((*um).c13.im > a) || ((*um).c21.re > a) || ((*um).c21.im > a) ||
	 ((*um).c22.re > a) || ((*um).c22.im > a) || ((*um).c23.re > a) || ((*um).c23.im > a) ||
	 ((*um).c31.re > a) || ((*um).c31.im > a) || ((*um).c32.re > a) || ((*um).c32.im > a) ||
	 ((*um).c33.re > a) || ((*um).c33.im > a)) {
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
  if(__isnan((*um).c11.re)|| __isnan((*um).c11.im) || __isnan((*um).c12.re) || __isnan((*um).c12.im) ||
     __isnan((*um).c13.re) || __isnan((*um).c13.im) || __isnan((*um).c21.re) || __isnan((*um).c21.im) ||
     __isnan((*um).c22.re) || __isnan((*um).c22.im) || __isnan((*um).c23.re) || __isnan((*um).c23.im) ||
     __isnan((*um).c31.re) || __isnan((*um).c31.im) || __isnan((*um).c32.re) || __isnan((*um).c32.im) ||
     __isnan((*um).c33.re) || __isnan((*um).c33.im)) {
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


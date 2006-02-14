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
#include "start.h"
#include "get_staples.h"

su3 get_staples(int x, int mu) {

  int k,iy;
  static su3 v,st;
  su3 *w1,*w2,*w3;
  
  _su3_zero(v);
  for(k=0;k<4;k++) {
    if(k!=mu){
      w1=&g_gauge_field[x][k];
      w2=&g_gauge_field[g_iup[x][k]][mu];
      w3=&g_gauge_field[g_iup[x][mu]][k];
      /* st = w2 * w3^d */
      _su3_times_su3d(st,*w2,*w3);
      /* v = v + w1 * st */
      _su3_times_su3_acc(v,*w1,st); 

      iy=g_idn[x][k];
      w1=&g_gauge_field[iy][k];
      w2=&g_gauge_field[iy][mu];
      w3=&g_gauge_field[g_iup[iy][mu]][k];
      /* st = w2 * w3 */
      _su3_times_su3(st,*w2,*w3);
      /* v = v + w1^d * st */
      _su3d_times_su3_acc(v,*w1,st);
    }
  }
  return v;
}

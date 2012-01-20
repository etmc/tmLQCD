/***********************************************************************
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008 Carsten Urbach
 *
 * Jenifer Gonzalez Lopez
 * (SF piece of the code)
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
#include <time.h>
#include "global.h"
#include "su3.h"
#include "su3adj.h"
#include "start.h"
#include "sf_get_staples.h"

su3 sf_get_staples(int x, int mu, su3 ** in_gauge_field) {
  
  int k, iy, flag1, flag2;
  static su3 v, st, cst;
  su3 *w1, *w2, *w3;
#ifdef _KOJAK_INST
#pragma pomp inst begin(staples)
#endif
  flag1 = -1;
  flag2 = -1;
  _su3_zero(v);
  for (k = 0; k < 4; k++) {
    if (k != mu) {
      if (g_t[x] > 1 && g_t[x] < (g_Tbsf - 1)) {
        flag1 = 0;
        flag2 = 0;
      } else if (g_t[x] == 0 && mu == 0) {
        flag1 = 1;  
        flag2 = 1;
      } else if (g_t[x] == 1 && mu == 0) {
        flag1 = 0;
        flag2 = 0;
      } else if (g_t[x] == 1 && mu != 0) {
        if (k != 0) {
          flag1 = 0;
          flag2 = 0;
        } else if (k == 0) {
          flag1 = 0;  
          flag2 = 1;
        }
      } else if (g_t[x] == (g_Tbsf - 1) && mu == 0) {
        flag1 = 1;  
        flag2 = 1;
      } else if (g_t[x] == (g_Tbsf - 1) && mu != 0) {
        if (k != 0) {
          flag1 = 0;
          flag2 = 0;
        } else if (k == 0) {
          flag1 = 1;  
          flag2 = 0;
        }
      }
      if (flag1 < 0 || flag2 < 0)
        exit(0);
      w1 = &in_gauge_field[x][k];
      w2 = &in_gauge_field[g_iup[x][k]][mu];
      w3 = &in_gauge_field[g_iup[x][mu]][k]; 
      _su3_times_su3d(st, *w2, *w3);
      _real_times_su3(cst, (flag1 == 0 ? g_rgi_C0 : g_Ct) ,st); /* specific to SF */
      _su3_times_su3_acc(v, *w1, cst); 
      iy = g_idn[x][k];
      w1 = &in_gauge_field[iy][k];
      w2 = &in_gauge_field[iy][mu];
      w3 = &in_gauge_field[g_iup[iy][mu]][k];
      _su3_times_su3(st, *w2, *w3);
      _real_times_su3(cst, (flag2 == 0 ? g_rgi_C0 : g_Ct) ,st); /* specific to SF */
      _su3d_times_su3_acc(v, *w1, cst);
    }
  }

  return v;
#ifdef _KOJAK_INST
#pragma pomp inst end(staples)
#endif 
}



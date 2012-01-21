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
#undef SSE
#undef SSE2
#undef SSE3
#include "global.h"
#include "su3.h"
#include "sf_get_rectangle_staples.h"

/* this function is valid ONLY IF Nt => 6 */
void sf_get_rectangle_staples(su3 * const v, const int x, const int mu) {
  static su3 tmp1, tmp2;
  int y, z, nu;
  su3 * a, * b, * c, * d, * e;
#ifdef _KOJAK_INST
#pragma pomp inst begin(rectstaples)
#endif
#ifdef XLC
#pragma disjoint(*v, tmp1, tmp2, *a, *b, *c, *d, *e)
#endif
  _su3_zero((*v));
  for(nu = 0; nu < 4; nu++) {
    if(mu != nu) {
      if (g_t[x] > 2 && g_t[x] < (g_Tbsf-2)) {
        /* first contr. starting from x 
          * a b c e^+ d^+
          *   c
          *   _
          * b| |e
          * a| |d
          */
        a = &g_gauge_field[x][nu];
        y = g_iup[x][nu];
        b = &g_gauge_field[y][nu];
        _su3_times_su3(tmp1, *a, *b);
        z = g_iup[y][nu];
        c = &g_gauge_field[z][mu];
        _su3_times_su3(tmp2, tmp1, *c);
        
        y = g_iup[x][mu];
        d = &g_gauge_field[y][nu];
        z = g_iup[y][nu];
        e = &g_gauge_field[z][nu];
        _su3_times_su3(tmp1, *d, *e);
        /* tmp1 = c1 * tmp1 */
        _real_times_su3(tmp1,g_rgi_C1,tmp1); /* that is the new thing specific of SF */
        _su3_times_su3d_acc((*v), tmp2, tmp1);
        
        /* 1 contr. starting idn[idn[x][nu]][nu] 
          * e^+ d^+ a b c
          *
          *e| |c
          *d|_|b
          *  a
          */
        y = g_idn[x][nu];
        z = g_idn[y][nu];
        d = &g_gauge_field[z][nu];
        a = &g_gauge_field[z][mu];
        _su3d_times_su3(tmp1, *d, *a);
        e = &g_gauge_field[y][nu];
        _su3d_times_su3(tmp2, *e, tmp1);
        
        y = g_iup[z][mu];
        b = &g_gauge_field[y][nu];
        z = g_iup[y][nu];
        c = &g_gauge_field[z][nu];
        _su3_times_su3(tmp1, *b, *c);
        /* tmp1 = c1 * tmp1 */
        _real_times_su3(tmp1,g_rgi_C1,tmp1); /* that is the new thing specific of SF */
        _su3_times_su3_acc((*v), tmp2, tmp1);
        
        /* second contr. starting from x 
          * a b c e^+ d^+
          *
          *   bc
          *   __
          * a| _|e
          *    d
          */
        a = &g_gauge_field[x][nu];
        y = g_iup[x][nu];
        b = &g_gauge_field[y][mu];
        _su3_times_su3(tmp1, *a, *b);
        z = g_iup[y][mu];
        c = &g_gauge_field[z][mu];
        _su3_times_su3(tmp2, tmp1, *c);
        
        y = g_iup[x][mu];
        d = &g_gauge_field[y][mu];
        z = g_iup[y][mu];
        e = &g_gauge_field[z][nu];
        _su3_times_su3(tmp1, *d, *e);
        /* tmp1 = c1 * tmp1 */
        _real_times_su3(tmp1,g_rgi_C1,tmp1); /* that is the new thing specific of SF */
        _su3_times_su3d_acc((*v), tmp2, tmp1);
        
        /* 1 contr. starting idn[x][nu] 
          * d^+ a b c e^+
          *
          *    e
          *    _
          * d|__|c
          *   ab
          */
        y = g_idn[x][nu];
        d = &g_gauge_field[y][nu];
        a = &g_gauge_field[y][mu];
        _su3d_times_su3(tmp1, *d, *a);
        z = g_iup[y][mu];
        b = &g_gauge_field[z][mu];
        _su3_times_su3(tmp2, tmp1, *b);
        
        y = g_iup[z][mu];
        c = &g_gauge_field[y][nu];
        z = g_iup[x][mu];
        e = &g_gauge_field[z][mu];
        _su3_times_su3d(tmp1, *c, *e);
        /* tmp1 = c1 * tmp1 */
        _real_times_su3(tmp1,g_rgi_C1,tmp1); /* that is the new thing specific of SF */
        _su3_times_su3_acc((*v), tmp2, tmp1);
        
        /* 1 contr. starting idn[idn[x][mu]][nu] 
          *  e^+ d^+ a b c
          *
          *  e
          *  _
          *d|__|c
          *  ab
          */
        y = g_idn[x][mu];
        z = g_idn[y][nu];
        d = &g_gauge_field[z][nu];
        a = &g_gauge_field[z][mu];
        _su3d_times_su3(tmp1, *d, *a);
        e = &g_gauge_field[y][mu];
        _su3d_times_su3(tmp2, *e, tmp1);
        
        y = g_idn[x][nu];
        b = &g_gauge_field[y][mu];
        z = g_iup[y][mu];
        c = &g_gauge_field[z][nu];
        _su3_times_su3(tmp1, *b, *c);
        /* tmp1 = c1 * tmp1 */
        _real_times_su3(tmp1,g_rgi_C1,tmp1); /* that is the new thing specific of SF */
        _su3_times_su3_acc((*v), tmp2, tmp1); 
        
        /* 1 contr. starting idn[x][mu] 
          * d^+ a b c e^+
          *
          *  bc
          *  __
          *a|_ |e
          *  d
          */
        y = g_idn[x][mu];
        d = &g_gauge_field[y][mu];
        z = g_iup[y][nu];
        a = &g_gauge_field[y][nu];
        _su3d_times_su3(tmp1, *d, *a);
        b = &g_gauge_field[z][mu];
        _su3_times_su3(tmp2, tmp1, *b);
        
        y = g_iup[x][mu];
        e = &g_gauge_field[y][nu];
        z = g_iup[x][nu];
        c = &g_gauge_field[z][mu];
        _su3_times_su3d(tmp1, *c, *e);
        /* tmp1 = c1 * tmp1 */
        _real_times_su3(tmp1,g_rgi_C1,tmp1); /* that is the new thing specific of SF */
        _su3_times_su3_acc((*v), tmp2, tmp1);
      }
      else if (g_t[x] == 2) {
        if (mu == 0 || (mu != 0 && nu != 0)){
          
          /* first contr. starting from x 
            * a b c e^+ d^+
            *   c
            *   _
            * b| |e
            * a| |d
            */
          a = &g_gauge_field[x][nu];
          y = g_iup[x][nu];
          b = &g_gauge_field[y][nu];
          _su3_times_su3(tmp1, *a, *b);
          z = g_iup[y][nu];
          c = &g_gauge_field[z][mu];
          _su3_times_su3(tmp2, tmp1, *c);
          
          y = g_iup[x][mu];
          d = &g_gauge_field[y][nu];
          z = g_iup[y][nu];
          e = &g_gauge_field[z][nu];
          _su3_times_su3(tmp1, *d, *e);
          /* tmp1 = c1 * tmp1 */
          _real_times_su3(tmp1,g_rgi_C1,tmp1); /* that is the new thing specific of SF */
          _su3_times_su3d_acc((*v), tmp2, tmp1);
          
          /* 1 contr. starting idn[idn[x][nu]][nu] 
            * e^+ d^+ a b c
            *
            *e| |c
            *d|_|b
            *  a
            */
          y = g_idn[x][nu];
          z = g_idn[y][nu];
          d = &g_gauge_field[z][nu];
          a = &g_gauge_field[z][mu];
          _su3d_times_su3(tmp1, *d, *a);
          e = &g_gauge_field[y][nu];
          _su3d_times_su3(tmp2, *e, tmp1);
          
          y = g_iup[z][mu];
          b = &g_gauge_field[y][nu];
          z = g_iup[y][nu];
          c = &g_gauge_field[z][nu];
          _su3_times_su3(tmp1, *b, *c);
          /* tmp1 = c1 * tmp1 */
          _real_times_su3(tmp1,g_rgi_C1,tmp1); /* that is the new thing specific of SF */
          _su3_times_su3_acc((*v), tmp2, tmp1);
          
      #if 1 /* PROBLEMS HERE for T=4 */
          
          /* second contr. starting from x 
            * a b c e^+ d^+
            *
            *   bc
            *   __
            * a| _|e
            *    d
            */
          a = &g_gauge_field[x][nu];
          y = g_iup[x][nu];
          b = &g_gauge_field[y][mu];
          _su3_times_su3(tmp1, *a, *b);
          z = g_iup[y][mu];
          c = &g_gauge_field[z][mu];
          _su3_times_su3(tmp2, tmp1, *c);
          
          y = g_iup[x][mu];
          d = &g_gauge_field[y][mu];
          z = g_iup[y][mu];
          e = &g_gauge_field[z][nu];
          _su3_times_su3(tmp1, *d, *e);
          /* tmp1 = c1 * tmp1 */
          _real_times_su3(tmp1,g_rgi_C1,tmp1); /* that is the new thing specific of SF */
          _su3_times_su3d_acc((*v), tmp2, tmp1);
      #endif
          
      #if 1 /*PROBLEMS HERE for T=4 */
          
          /* 1 contr. starting idn[x][nu] 
            * d^+ a b c e^+
            *
            *    e
            *    _
            * d|__|c
            *   ab
            */
          y = g_idn[x][nu];
          d = &g_gauge_field[y][nu];
          a = &g_gauge_field[y][mu];
          _su3d_times_su3(tmp1, *d, *a);
          z = g_iup[y][mu];
          b = &g_gauge_field[z][mu];
          _su3_times_su3(tmp2, tmp1, *b);
          
          y = g_iup[z][mu];
          c = &g_gauge_field[y][nu];
          z = g_iup[x][mu];
          e = &g_gauge_field[z][mu];
          _su3_times_su3d(tmp1, *c, *e);
          /* tmp1 = c1 * tmp1 */
          _real_times_su3(tmp1,g_rgi_C1,tmp1); /* that is the new thing specific of SF */
          _su3_times_su3_acc((*v), tmp2, tmp1);
      #endif
          
          
          /* 1 contr. starting idn[idn[x][mu]][nu] 
            *  e^+ d^+ a b c
            *
            *  e
            *  _
            *d|__|c
            *  ab
            */
          y = g_idn[x][mu];
          z = g_idn[y][nu];
          d = &g_gauge_field[z][nu];
          a = &g_gauge_field[z][mu];
          _su3d_times_su3(tmp1, *d, *a);
          e = &g_gauge_field[y][mu];
          _su3d_times_su3(tmp2, *e, tmp1);
          
          y = g_idn[x][nu];
          b = &g_gauge_field[y][mu];
          z = g_iup[y][mu];
          c = &g_gauge_field[z][nu];
          _su3_times_su3(tmp1, *b, *c);
          /* tmp1 = c1 * tmp1 */
          _real_times_su3(tmp1,g_rgi_C1,tmp1); /* that is the new thing specific of SF */
          _su3_times_su3_acc((*v), tmp2, tmp1); 
          
          /* 1 contr. starting idn[x][mu] 
            * d^+ a b c e^+
            *
            *  bc
            *  __
            *a|_ |e
            *  d
            */
          y = g_idn[x][mu];
          d = &g_gauge_field[y][mu];
          z = g_iup[y][nu];
          a = &g_gauge_field[y][nu];
          _su3d_times_su3(tmp1, *d, *a);
          b = &g_gauge_field[z][mu];
          _su3_times_su3(tmp2, tmp1, *b);
          
          y = g_iup[x][mu];
          e = &g_gauge_field[y][nu];
          z = g_iup[x][nu];
          c = &g_gauge_field[z][mu];
          _su3_times_su3d(tmp1, *c, *e);
          /* tmp1 = c1 * tmp1 */
          _real_times_su3(tmp1,g_rgi_C1,tmp1); /* that is the new thing specific of SF */
          _su3_times_su3_acc((*v), tmp2, tmp1);
        }
        else if (mu != 0 && nu == 0){
          
      #if 1 /* PROBLEMS HERE for T=4 */
          
          /* first contr. starting from x 
            * a b c e^+ d^+
            *   c
            *   _
            * b| |e
            * a| |d
            */
          a = &g_gauge_field[x][nu];
          y = g_iup[x][nu];
          b = &g_gauge_field[y][nu];
          _su3_times_su3(tmp1, *a, *b);
          z = g_iup[y][nu];
          c = &g_gauge_field[z][mu];
          _su3_times_su3(tmp2, tmp1, *c);
          
          y = g_iup[x][mu];
          d = &g_gauge_field[y][nu];
          z = g_iup[y][nu];
          e = &g_gauge_field[z][nu];
          _su3_times_su3(tmp1, *d, *e);
          /* tmp1 = c1 * tmp1 */
          _real_times_su3(tmp1,g_rgi_C1,tmp1); /* that is the new thing specific of SF */
          _su3_times_su3d_acc((*v), tmp2, tmp1);
      #endif
          
          /* 1 contr. starting idn[idn[x][nu]][nu] 
            * e^+ d^+ a b c
            *
            *e| |c
            *d|_|b
            *  a
            */
          y = g_idn[x][nu];
          z = g_idn[y][nu];
          d = &g_gauge_field[z][nu];
          a = &g_gauge_field[z][mu];
          _su3d_times_su3(tmp1, *d, *a);
          e = &g_gauge_field[y][nu];
          _su3d_times_su3(tmp2, *e, tmp1);
          
          y = g_iup[z][mu];
          b = &g_gauge_field[y][nu];
          z = g_iup[y][nu];
          c = &g_gauge_field[z][nu];
          _su3_times_su3(tmp1, *b, *c);
          /* tmp1 = c1tts * tmp1 */
          _real_times_su3(tmp1,g_C1tts,tmp1); /* that is the new thing specific of SF */
          _su3_times_su3_acc((*v), tmp2, tmp1);
          
          
          /* second contr. starting from x 
            * a b c e^+ d^+
            *
            *   bc
            *   __
            * a| _|e
            *    d
            */
          a = &g_gauge_field[x][nu];
          y = g_iup[x][nu];
          b = &g_gauge_field[y][mu];
          _su3_times_su3(tmp1, *a, *b);
          z = g_iup[y][mu];
          c = &g_gauge_field[z][mu];
          _su3_times_su3(tmp2, tmp1, *c);
          
          y = g_iup[x][mu];
          d = &g_gauge_field[y][mu];
          z = g_iup[y][mu];
          e = &g_gauge_field[z][nu];
          _su3_times_su3(tmp1, *d, *e);
          /* tmp1 = c1 * tmp1 */
          _real_times_su3(tmp1,g_rgi_C1,tmp1); /* that is the new thing specific of SF */
          _su3_times_su3d_acc((*v), tmp2, tmp1);
          
          /* 1 contr. starting idn[x][nu] 
            * d^+ a b c e^+
            *
            *    e
            *    _
            * d|__|c
            *   ab
            */
          y = g_idn[x][nu];
          d = &g_gauge_field[y][nu];
          a = &g_gauge_field[y][mu];
          _su3d_times_su3(tmp1, *d, *a);
          z = g_iup[y][mu];
          b = &g_gauge_field[z][mu];
          _su3_times_su3(tmp2, tmp1, *b);
          
          y = g_iup[z][mu];
          c = &g_gauge_field[y][nu];
          z = g_iup[x][mu];
          e = &g_gauge_field[z][mu];
          _su3_times_su3d(tmp1, *c, *e);
          /* tmp1 = c1 * tmp1 */
          _real_times_su3(tmp1,g_rgi_C1,tmp1); /* that is the new thing specific of SF */
          _su3_times_su3_acc((*v), tmp2, tmp1);
          
          /* 1 contr. starting idn[idn[x][mu]][nu] 
            *  e^+ d^+ a b c
            *
            *  e
            *  _
            *d|__|c
            *  ab
            */
          y = g_idn[x][mu];
          z = g_idn[y][nu];
          d = &g_gauge_field[z][nu];
          a = &g_gauge_field[z][mu];
          _su3d_times_su3(tmp1, *d, *a);
          e = &g_gauge_field[y][mu];
          _su3d_times_su3(tmp2, *e, tmp1);
          
          y = g_idn[x][nu];
          b = &g_gauge_field[y][mu];
          z = g_iup[y][mu];
          c = &g_gauge_field[z][nu];
          _su3_times_su3(tmp1, *b, *c);
          /* tmp1 = c1 * tmp1 */
          _real_times_su3(tmp1,g_rgi_C1,tmp1); /* that is the new thing specific of SF */
          _su3_times_su3_acc((*v), tmp2, tmp1); 
        
          /* 1 contr. starting idn[x][mu] 
            * d^+ a b c e^+
            *
            *  bc
            *  __
            *a|_ |e
            *  d
            */
          y = g_idn[x][mu];
          d = &g_gauge_field[y][mu];
          z = g_iup[y][nu];
          a = &g_gauge_field[y][nu];
          _su3d_times_su3(tmp1, *d, *a);
          b = &g_gauge_field[z][mu];
          _su3_times_su3(tmp2, tmp1, *b);
          
          y = g_iup[x][mu];
          e = &g_gauge_field[y][nu];
          z = g_iup[x][nu];
          c = &g_gauge_field[z][mu];
          _su3_times_su3d(tmp1, *c, *e);
          /* tmp1 = c1 * tmp1 */
          _real_times_su3(tmp1,g_rgi_C1,tmp1); /* that is the new thing specific of SF */
          _su3_times_su3_acc((*v), tmp2, tmp1);
        }	
      }
      else if (g_t[x] == (g_Tbsf-2)) {
        if (mu == 0){
          /* first contr. starting from x 
            * a b c e^+ d^+
            *   c
            *   _
            * b| |e
            * a| |d
            */
          a = &g_gauge_field[x][nu];
          y = g_iup[x][nu];
          b = &g_gauge_field[y][nu];
          _su3_times_su3(tmp1, *a, *b);
          z = g_iup[y][nu];
          c = &g_gauge_field[z][mu];
          _su3_times_su3(tmp2, tmp1, *c);
          
          y = g_iup[x][mu];
          d = &g_gauge_field[y][nu];
          z = g_iup[y][nu];
          e = &g_gauge_field[z][nu];
          _su3_times_su3(tmp1, *d, *e);
          /* tmp1 = c1 * tmp1 */
          _real_times_su3(tmp1,g_rgi_C1,tmp1); /* that is the new thing specific of SF */
          _su3_times_su3d_acc((*v), tmp2, tmp1);
          
          /* 1 contr. starting idn[idn[x][nu]][nu] 
            * e^+ d^+ a b c
            *
            *e| |c
            *d|_|b
            *  a
            */
          y = g_idn[x][nu];
          z = g_idn[y][nu];
          d = &g_gauge_field[z][nu];
          a = &g_gauge_field[z][mu];
          _su3d_times_su3(tmp1, *d, *a);
          e = &g_gauge_field[y][nu];
          _su3d_times_su3(tmp2, *e, tmp1);
          
          y = g_iup[z][mu];
          b = &g_gauge_field[y][nu];
          z = g_iup[y][nu];
          c = &g_gauge_field[z][nu];
          _su3_times_su3(tmp1, *b, *c);
          /* tmp1 = c1 * tmp1 */
          _real_times_su3(tmp1,g_rgi_C1,tmp1); /* that is the new thing specific of SF */
          _su3_times_su3_acc((*v), tmp2, tmp1);
          
          /* second contr. starting from x 
            * a b c e^+ d^+
            *
            *   bc
            *   __
            * a| _|e
            *    d
            */
          a = &g_gauge_field[x][nu];
          y = g_iup[x][nu];
          b = &g_gauge_field[y][mu];
          _su3_times_su3(tmp1, *a, *b);
          z = g_iup[y][mu];
          c = &g_gauge_field[z][mu];
          _su3_times_su3(tmp2, tmp1, *c);
          
          y = g_iup[x][mu];
          d = &g_gauge_field[y][mu];
          z = g_iup[y][mu];
          e = &g_gauge_field[z][nu];
          _su3_times_su3(tmp1, *d, *e);
          /* tmp1 = c1tts * tmp1 */
          _real_times_su3(tmp1,g_C1tts,tmp1); /* that is the new thing specific of SF */
          _su3_times_su3d_acc((*v), tmp2, tmp1);
          
          /* 1 contr. starting idn[x][nu] 
            * d^+ a b c e^+
            *
            *    e
            *    _
            * d|__|c
            *   ab
            */
          y = g_idn[x][nu];
          d = &g_gauge_field[y][nu];
          a = &g_gauge_field[y][mu];
          _su3d_times_su3(tmp1, *d, *a);
          z = g_iup[y][mu];
          b = &g_gauge_field[z][mu];
          _su3_times_su3(tmp2, tmp1, *b);
          
          y = g_iup[z][mu];
          c = &g_gauge_field[y][nu];
          z = g_iup[x][mu];
          e = &g_gauge_field[z][mu];
          _su3_times_su3d(tmp1, *c, *e);
          /* tmp1 = c1tts * tmp1 */
          _real_times_su3(tmp1,g_C1tts,tmp1); /* that is the new thing specific of SF */
          _su3_times_su3_acc((*v), tmp2, tmp1);
          
          /* 1 contr. starting idn[idn[x][mu]][nu] 
            *  e^+ d^+ a b c
            *
            *  e
            *  _
            *d|__|c
            *  ab
            */
          y = g_idn[x][mu];
          z = g_idn[y][nu];
          d = &g_gauge_field[z][nu];
          a = &g_gauge_field[z][mu];
          _su3d_times_su3(tmp1, *d, *a);
          e = &g_gauge_field[y][mu];
          _su3d_times_su3(tmp2, *e, tmp1);
          
          y = g_idn[x][nu];
          b = &g_gauge_field[y][mu];
          z = g_iup[y][mu];
          c = &g_gauge_field[z][nu];
          _su3_times_su3(tmp1, *b, *c);
          /* tmp1 = c1 * tmp1 */
          _real_times_su3(tmp1,g_rgi_C1,tmp1); /* that is the new thing specific of SF */
          _su3_times_su3_acc((*v), tmp2, tmp1); 
          
          /* 1 contr. starting idn[x][mu] 
            * d^+ a b c e^+
            *
            *  bc
            *  __
            *a|_ |e
            *  d
            */
          y = g_idn[x][mu];
          d = &g_gauge_field[y][mu];
          z = g_iup[y][nu];
          a = &g_gauge_field[y][nu];
          _su3d_times_su3(tmp1, *d, *a);
          b = &g_gauge_field[z][mu];
          _su3_times_su3(tmp2, tmp1, *b);
        
          y = g_iup[x][mu];
          e = &g_gauge_field[y][nu];
          z = g_iup[x][nu];
          c = &g_gauge_field[z][mu];
          _su3_times_su3d(tmp1, *c, *e);
          /* tmp1 = c1 * tmp1 */
          _real_times_su3(tmp1,g_rgi_C1,tmp1); /* that is the new thing specific of SF */
          _su3_times_su3_acc((*v), tmp2, tmp1);
        }
        
        else if (mu != 0 && nu == 0){
          /* first contr. starting from x 
            * a b c e^+ d^+
            *   c
            *   _
            * b| |e
            * a| |d
            */
          a = &g_gauge_field[x][nu];
          y = g_iup[x][nu];
          b = &g_gauge_field[y][nu];
          _su3_times_su3(tmp1, *a, *b);
          z = g_iup[y][nu];
          c = &g_gauge_field[z][mu];
          _su3_times_su3(tmp2, tmp1, *c);
          
          y = g_iup[x][mu];
          d = &g_gauge_field[y][nu];
          z = g_iup[y][nu];
          e = &g_gauge_field[z][nu];
          _su3_times_su3(tmp1, *d, *e);
          /* tmp1 = c1tts * tmp1 */
          _real_times_su3(tmp1,g_C1tts,tmp1); /* that is the new thing specific of SF */
          _su3_times_su3d_acc((*v), tmp2, tmp1);
          
      #if 1 /* PROBLEMS HERE for T=4 */
          
          /* 1 contr. starting idn[idn[x][nu]][nu] 
            * e^+ d^+ a b c
            *
            *e| |c
            *d|_|b
            *  a
            */
          y = g_idn[x][nu];
          z = g_idn[y][nu];
          d = &g_gauge_field[z][nu];
          a = &g_gauge_field[z][mu];
          _su3d_times_su3(tmp1, *d, *a);
          e = &g_gauge_field[y][nu];
          _su3d_times_su3(tmp2, *e, tmp1);
          
          y = g_iup[z][mu];
          b = &g_gauge_field[y][nu];
          z = g_iup[y][nu];
          c = &g_gauge_field[z][nu];
          _su3_times_su3(tmp1, *b, *c);
          /* tmp1 = c1 * tmp1 */
          _real_times_su3(tmp1,g_rgi_C1,tmp1); /* that is the new thing specific of SF */
          _su3_times_su3_acc((*v), tmp2, tmp1);
      #endif
          
          
          /* second contr. starting from x 
            * a b c e^+ d^+
            *
            *   bc
            *   __
            * a| _|e
            *    d
            */
          a = &g_gauge_field[x][nu];
          y = g_iup[x][nu];
          b = &g_gauge_field[y][mu];
          _su3_times_su3(tmp1, *a, *b);
          z = g_iup[y][mu];
          c = &g_gauge_field[z][mu];
          _su3_times_su3(tmp2, tmp1, *c);
          
          y = g_iup[x][mu];
          d = &g_gauge_field[y][mu];
          z = g_iup[y][mu];
          e = &g_gauge_field[z][nu];
          _su3_times_su3(tmp1, *d, *e);
          /* tmp1 = c1 * tmp1 */
          _real_times_su3(tmp1,g_rgi_C1,tmp1); /* that is the new thing specific of SF */
          _su3_times_su3d_acc((*v), tmp2, tmp1);
          
          /* 1 contr. starting idn[x][nu] 
            * d^+ a b c e^+
            *
            *    e
            *    _
            * d|__|c
            *   ab
            */
          y = g_idn[x][nu];
          d = &g_gauge_field[y][nu];
          a = &g_gauge_field[y][mu];
          _su3d_times_su3(tmp1, *d, *a);
          z = g_iup[y][mu];
          b = &g_gauge_field[z][mu];
          _su3_times_su3(tmp2, tmp1, *b);
          
          y = g_iup[z][mu];
          c = &g_gauge_field[y][nu];
          z = g_iup[x][mu];
          e = &g_gauge_field[z][mu];
          _su3_times_su3d(tmp1, *c, *e);
          /* tmp1 = c1 * tmp1 */
          _real_times_su3(tmp1,g_rgi_C1,tmp1); /* that is the new thing specific of SF */
          _su3_times_su3_acc((*v), tmp2, tmp1);
          
          /* 1 contr. starting idn[idn[x][mu]][nu] 
            *  e^+ d^+ a b c
            *
            *  e
            *  _
            *d|__|c
            *  ab
            */
          y = g_idn[x][mu];
          z = g_idn[y][nu];
          d = &g_gauge_field[z][nu];
          a = &g_gauge_field[z][mu];
          _su3d_times_su3(tmp1, *d, *a);
          e = &g_gauge_field[y][mu];
          _su3d_times_su3(tmp2, *e, tmp1);
          
          y = g_idn[x][nu];
          b = &g_gauge_field[y][mu];
          z = g_iup[y][mu];
          c = &g_gauge_field[z][nu];
          _su3_times_su3(tmp1, *b, *c);
          /* tmp1 = c1 * tmp1 */
          _real_times_su3(tmp1,g_rgi_C1,tmp1); /* that is the new thing specific of SF */
          _su3_times_su3_acc((*v), tmp2, tmp1); 
          
          /* 1 contr. starting idn[x][mu] 
            * d^+ a b c e^+
            *
            *  bc
            *  __
            *a|_ |e
            *  d
            */
          y = g_idn[x][mu];
          d = &g_gauge_field[y][mu];
          z = g_iup[y][nu];
          a = &g_gauge_field[y][nu];
          _su3d_times_su3(tmp1, *d, *a);
          b = &g_gauge_field[z][mu];
          _su3_times_su3(tmp2, tmp1, *b);
          
          y = g_iup[x][mu];
          e = &g_gauge_field[y][nu];
          z = g_iup[x][nu];
          c = &g_gauge_field[z][mu];
          _su3_times_su3d(tmp1, *c, *e);
          /* tmp1 = c1 * tmp1 */
          _real_times_su3(tmp1,g_rgi_C1,tmp1); /* that is the new thing specific of SF */
          _su3_times_su3_acc((*v), tmp2, tmp1);
        }
        else if (mu != 0 && nu != 0){
          /* first contr. starting from x 
            * a b c e^+ d^+
            *   c
            *   _
            * b| |e
            * a| |d
        */
          a = &g_gauge_field[x][nu];
          y = g_iup[x][nu];
          b = &g_gauge_field[y][nu];
          _su3_times_su3(tmp1, *a, *b);
          z = g_iup[y][nu];
          c = &g_gauge_field[z][mu];
          _su3_times_su3(tmp2, tmp1, *c);
          
          y = g_iup[x][mu];
          d = &g_gauge_field[y][nu];
          z = g_iup[y][nu];
          e = &g_gauge_field[z][nu];
          _su3_times_su3(tmp1, *d, *e);
          /* tmp1 = c1 * tmp1 */
          _real_times_su3(tmp1,g_rgi_C1,tmp1); /* that is the new thing specific of SF */
          _su3_times_su3d_acc((*v), tmp2, tmp1);
          
          /* 1 contr. starting idn[idn[x][nu]][nu] 
            * e^+ d^+ a b c
            *
            *e| |c
            *d|_|b
            *  a
            */
          y = g_idn[x][nu];
          z = g_idn[y][nu];
          d = &g_gauge_field[z][nu];
          a = &g_gauge_field[z][mu];
          _su3d_times_su3(tmp1, *d, *a);
          e = &g_gauge_field[y][nu];
          _su3d_times_su3(tmp2, *e, tmp1);
          
          y = g_iup[z][mu];
          b = &g_gauge_field[y][nu];
          z = g_iup[y][nu];
          c = &g_gauge_field[z][nu];
          _su3_times_su3(tmp1, *b, *c);
          /* tmp1 = c1 * tmp1 */
          _real_times_su3(tmp1,g_rgi_C1,tmp1); /* that is the new thing specific of SF */
          _su3_times_su3_acc((*v), tmp2, tmp1);
          
          /* second contr. starting from x 
            * a b c e^+ d^+
            *
            *   bc
            *   __
            * a| _|e
            *    d
            */
          a = &g_gauge_field[x][nu];
          y = g_iup[x][nu];
          b = &g_gauge_field[y][mu];
          _su3_times_su3(tmp1, *a, *b);
          z = g_iup[y][mu];
          c = &g_gauge_field[z][mu];
          _su3_times_su3(tmp2, tmp1, *c);
          
          y = g_iup[x][mu];
          d = &g_gauge_field[y][mu];
          z = g_iup[y][mu];
          e = &g_gauge_field[z][nu];
          _su3_times_su3(tmp1, *d, *e);
          /* tmp1 = c1 * tmp1 */
          _real_times_su3(tmp1,g_rgi_C1,tmp1); /* that is the new thing specific of SF */
          _su3_times_su3d_acc((*v), tmp2, tmp1);
          
          /* 1 contr. starting idn[x][nu] 
            * d^+ a b c e^+
            *
            *    e
            *    _
            * d|__|c
            *   ab
            */
          y = g_idn[x][nu];
          d = &g_gauge_field[y][nu];
          a = &g_gauge_field[y][mu];
          _su3d_times_su3(tmp1, *d, *a);
          z = g_iup[y][mu];
          b = &g_gauge_field[z][mu];
          _su3_times_su3(tmp2, tmp1, *b);
          
          y = g_iup[z][mu];
          c = &g_gauge_field[y][nu];
          z = g_iup[x][mu];
          e = &g_gauge_field[z][mu];
          _su3_times_su3d(tmp1, *c, *e);
          /* tmp1 = c1 * tmp1 */
          _real_times_su3(tmp1,g_rgi_C1,tmp1); /* that is the new thing specific of SF */
          _su3_times_su3_acc((*v), tmp2, tmp1);
          
          /* 1 contr. starting idn[idn[x][mu]][nu] 
            *  e^+ d^+ a b c
            *
            *  e
            *  _
            *d|__|c
            *  ab
            */
          y = g_idn[x][mu];
          z = g_idn[y][nu];
          d = &g_gauge_field[z][nu];
          a = &g_gauge_field[z][mu];
          _su3d_times_su3(tmp1, *d, *a);
          e = &g_gauge_field[y][mu];
          _su3d_times_su3(tmp2, *e, tmp1);
          
          y = g_idn[x][nu];
          b = &g_gauge_field[y][mu];
          z = g_iup[y][mu];
          c = &g_gauge_field[z][nu];
          _su3_times_su3(tmp1, *b, *c);
          /* tmp1 = c1 * tmp1 */
          _real_times_su3(tmp1,g_rgi_C1,tmp1); /* that is the new thing specific of SF */
          _su3_times_su3_acc((*v), tmp2, tmp1); 
          
          /* 1 contr. starting idn[x][mu] 
            * d^+ a b c e^+
            *
            *  bc
            *  __
            *a|_ |e
            *  d
            */
          y = g_idn[x][mu];
          d = &g_gauge_field[y][mu];
          z = g_iup[y][nu];
          a = &g_gauge_field[y][nu];
          _su3d_times_su3(tmp1, *d, *a);
          b = &g_gauge_field[z][mu];
          _su3_times_su3(tmp2, tmp1, *b);
          
          y = g_iup[x][mu];
          e = &g_gauge_field[y][nu];
          z = g_iup[x][nu];
          c = &g_gauge_field[z][mu];
          _su3_times_su3d(tmp1, *c, *e);
          /* tmp1 = c1 * tmp1 */
          _real_times_su3(tmp1,g_rgi_C1,tmp1); /* that is the new thing specific of SF */
          _su3_times_su3_acc((*v), tmp2, tmp1);
        }
      }
      else if (g_t[x] == 1) {
        if (mu == 0) {
          /* first contr. starting from x 
            * a b c e^+ d^+
            *   c
            *   _
            * b| |e
            * a| |d
            */
          a = &g_gauge_field[x][nu];
          y = g_iup[x][nu];
          b = &g_gauge_field[y][nu];
          _su3_times_su3(tmp1, *a, *b);
          z = g_iup[y][nu];
          c = &g_gauge_field[z][mu];
          _su3_times_su3(tmp2, tmp1, *c);
          
          y = g_iup[x][mu];
          d = &g_gauge_field[y][nu];
          z = g_iup[y][nu];
          e = &g_gauge_field[z][nu];
          _su3_times_su3(tmp1, *d, *e);
          /* tmp1 = c1 * tmp1 */
          _real_times_su3(tmp1,g_rgi_C1,tmp1); /* that is the new thing specific of SF */
          _su3_times_su3d_acc((*v), tmp2, tmp1);
          
          /* 1 contr. starting idn[idn[x][nu]][nu] 
            * e^+ d^+ a b c
            *
            *e| |c
            *d|_|b
            *  a
            */
          y = g_idn[x][nu];
          z = g_idn[y][nu];
          d = &g_gauge_field[z][nu];
          a = &g_gauge_field[z][mu];
          _su3d_times_su3(tmp1, *d, *a);
          e = &g_gauge_field[y][nu];
          _su3d_times_su3(tmp2, *e, tmp1);
          
          y = g_iup[z][mu];
          b = &g_gauge_field[y][nu];
          z = g_iup[y][nu];
          c = &g_gauge_field[z][nu];
          _su3_times_su3(tmp1, *b, *c);
          /* tmp1 = c1 * tmp1 */
          _real_times_su3(tmp1,g_rgi_C1,tmp1); /* that is the new thing specific of SF */
          _su3_times_su3_acc((*v), tmp2, tmp1);
          
          /* second contr. starting from x 
            * a b c e^+ d^+
            *
            *   bc
            *   __
            * a| _|e
            *    d
            */
          a = &g_gauge_field[x][nu];
          y = g_iup[x][nu];
          b = &g_gauge_field[y][mu];
          _su3_times_su3(tmp1, *a, *b);
          z = g_iup[y][mu];
          c = &g_gauge_field[z][mu];
          _su3_times_su3(tmp2, tmp1, *c);
          
          y = g_iup[x][mu];
          d = &g_gauge_field[y][mu];
          z = g_iup[y][mu];
          e = &g_gauge_field[z][nu];
          _su3_times_su3(tmp1, *d, *e);
          /* tmp1 = c1 * tmp1 */
          _real_times_su3(tmp1,g_rgi_C1,tmp1); /* that is the new thing specific of SF */
          _su3_times_su3d_acc((*v), tmp2, tmp1);
          
          /* 1 contr. starting idn[x][nu] 
            * d^+ a b c e^+
            *
            *    e
            *    _
            * d|__|c
            *   ab
            */
          y = g_idn[x][nu];
          d = &g_gauge_field[y][nu];
          a = &g_gauge_field[y][mu];
          _su3d_times_su3(tmp1, *d, *a);
          z = g_iup[y][mu];
          b = &g_gauge_field[z][mu];
          _su3_times_su3(tmp2, tmp1, *b);
          
          y = g_iup[z][mu];
          c = &g_gauge_field[y][nu];
          z = g_iup[x][mu];
          e = &g_gauge_field[z][mu];
          _su3_times_su3d(tmp1, *c, *e);
          /* tmp1 = c1 * tmp1 */
          _real_times_su3(tmp1,g_rgi_C1,tmp1); /* that is the new thing specific of SF */
          _su3_times_su3_acc((*v), tmp2, tmp1);
          
          /* 1 contr. starting idn[idn[x][mu]][nu] 
            *  e^+ d^+ a b c
            *
            *  e
            *  _
            *d|__|c
            *  ab
            */
          y = g_idn[x][mu];
          z = g_idn[y][nu];
          d = &g_gauge_field[z][nu];
          a = &g_gauge_field[z][mu];
          _su3d_times_su3(tmp1, *d, *a);
          e = &g_gauge_field[y][mu];
          _su3d_times_su3(tmp2, *e, tmp1);
          
          y = g_idn[x][nu];
          b = &g_gauge_field[y][mu];
          z = g_iup[y][mu];
          c = &g_gauge_field[z][nu];
          _su3_times_su3(tmp1, *b, *c);
          /* tmp1 = c1tts * tmp1 */
          _real_times_su3(tmp1,g_C1tts,tmp1); /* that is the new thing specific of SF */
          _su3_times_su3_acc((*v), tmp2, tmp1); 
          
          /* 1 contr. starting idn[x][mu] 
            * d^+ a b c e^+
            *
            *  bc
            *  __
            *a|_ |e
            *  d
            */
          y = g_idn[x][mu];
          d = &g_gauge_field[y][mu];
          z = g_iup[y][nu];
          a = &g_gauge_field[y][nu];
          _su3d_times_su3(tmp1, *d, *a);
          b = &g_gauge_field[z][mu];
          _su3_times_su3(tmp2, tmp1, *b);
          
          y = g_iup[x][mu];
          e = &g_gauge_field[y][nu];
          z = g_iup[x][nu];
          c = &g_gauge_field[z][mu];
          _su3_times_su3d(tmp1, *c, *e);
          /* tmp1 = c1tts * tmp1 */
          _real_times_su3(tmp1,g_C1tts,tmp1); /* that is the new thing specific of SF */
          _su3_times_su3_acc((*v), tmp2, tmp1);
        }
        else if (mu != 0 && nu == 0) {
          
          /* first contr. starting from x 
            * a b c e^+ d^+
            *   c
            *   _
            * b| |e
            * a| |d
            */
          a = &g_gauge_field[x][nu];
          y = g_iup[x][nu];
          b = &g_gauge_field[y][nu];
          _su3_times_su3(tmp1, *a, *b);
          z = g_iup[y][nu];
          c = &g_gauge_field[z][mu];
          _su3_times_su3(tmp2, tmp1, *c);
          
          y = g_iup[x][mu];
          d = &g_gauge_field[y][nu];
          z = g_iup[y][nu];
          e = &g_gauge_field[z][nu];
          _su3_times_su3(tmp1, *d, *e);
          /* tmp1 = c1 * tmp1 */
          _real_times_su3(tmp1,g_rgi_C1,tmp1); /* that is the new thing specific of SF */
          _su3_times_su3d_acc((*v), tmp2, tmp1);
          
          /* second contr. starting from x 
            * a b c e^+ d^+
            *
            *   bc
            *   __
            * a| _|e
            *    d
            */
          a = &g_gauge_field[x][nu];
          y = g_iup[x][nu];
          b = &g_gauge_field[y][mu];
          _su3_times_su3(tmp1, *a, *b);
          z = g_iup[y][mu];
          c = &g_gauge_field[z][mu];
          _su3_times_su3(tmp2, tmp1, *c);
          
          y = g_iup[x][mu];
          d = &g_gauge_field[y][mu];
          z = g_iup[y][mu];
          e = &g_gauge_field[z][nu];
          _su3_times_su3(tmp1, *d, *e);
          /* tmp1 = c1 * tmp1 */
          _real_times_su3(tmp1,g_rgi_C1,tmp1); /* that is the new thing specific of SF */
          _su3_times_su3d_acc((*v), tmp2, tmp1);
          
          /* 1 contr. starting idn[x][nu] 
            * d^+ a b c e^+
            *
            *    e
            *    _
            * d|__|c
            *   ab
            */
          y = g_idn[x][nu];
          d = &g_gauge_field[y][nu];
          a = &g_gauge_field[y][mu];
          _su3d_times_su3(tmp1, *d, *a);
          z = g_iup[y][mu];
          b = &g_gauge_field[z][mu];
          _su3_times_su3(tmp2, tmp1, *b);
          
          y = g_iup[z][mu];
          c = &g_gauge_field[y][nu];
          z = g_iup[x][mu];
          e = &g_gauge_field[z][mu];
          _su3_times_su3d(tmp1, *c, *e);
          /* tmp1 = c1tss * tmp1 */
          _real_times_su3(tmp1,g_C1tss,tmp1); /* that is the new thing specific of SF */
          _su3_times_su3_acc((*v), tmp2, tmp1);
          
          /* 1 contr. starting idn[idn[x][mu]][nu] 
            *  e^+ d^+ a b c
            *
            *  e
            *  _
            *d|__|c
            *  ab
            */
          y = g_idn[x][mu];
          z = g_idn[y][nu];
          d = &g_gauge_field[z][nu];
          a = &g_gauge_field[z][mu];
          _su3d_times_su3(tmp1, *d, *a);
          e = &g_gauge_field[y][mu];
          _su3d_times_su3(tmp2, *e, tmp1);
          
          y = g_idn[x][nu];
          b = &g_gauge_field[y][mu];
          z = g_iup[y][mu];
          c = &g_gauge_field[z][nu];
          _su3_times_su3(tmp1, *b, *c);
          /* tmp1 = c1tss * tmp1 */
          _real_times_su3(tmp1,g_C1tss,tmp1); /* that is the new thing specific of SF */
          _su3_times_su3_acc((*v), tmp2, tmp1); 
          
          /* 1 contr. starting idn[x][mu] 
            * d^+ a b c e^+
            *
            *  bc
            *  __
            *a|_ |e
            *  d
            */
          y = g_idn[x][mu];
          d = &g_gauge_field[y][mu];
          z = g_iup[y][nu];
          a = &g_gauge_field[y][nu];
          _su3d_times_su3(tmp1, *d, *a);
          b = &g_gauge_field[z][mu];
          _su3_times_su3(tmp2, tmp1, *b);
          
          y = g_iup[x][mu];
          e = &g_gauge_field[y][nu];
          z = g_iup[x][nu];
          c = &g_gauge_field[z][mu];
          _su3_times_su3d(tmp1, *c, *e);
          /* tmp1 = c1 * tmp1 */
          _real_times_su3(tmp1,g_rgi_C1,tmp1); /* that is the new thing specific of SF */
          _su3_times_su3_acc((*v), tmp2, tmp1);
        }
        else if (mu != 0 && nu != 0) {
          /* first contr. starting from x 
            * a b c e^+ d^+
            *   c
            *   _
            * b| |e
            * a| |d
            */
          a = &g_gauge_field[x][nu];
          y = g_iup[x][nu];
          b = &g_gauge_field[y][nu];
          _su3_times_su3(tmp1, *a, *b);
          z = g_iup[y][nu];
          c = &g_gauge_field[z][mu];
          _su3_times_su3(tmp2, tmp1, *c);
          
          y = g_iup[x][mu];
          d = &g_gauge_field[y][nu];
          z = g_iup[y][nu];
          e = &g_gauge_field[z][nu];
          _su3_times_su3(tmp1, *d, *e);
          /* tmp1 = c1 * tmp1 */
          _real_times_su3(tmp1,g_rgi_C1,tmp1); /* that is the new thing specific of SF */
          _su3_times_su3d_acc((*v), tmp2, tmp1);
          
          /* 1 contr. starting idn[idn[x][nu]][nu] 
            * e^+ d^+ a b c
            *
            *e| |c
            *d|_|b
            *  a
            */
          y = g_idn[x][nu];
          z = g_idn[y][nu];
          d = &g_gauge_field[z][nu];
          a = &g_gauge_field[z][mu];
          _su3d_times_su3(tmp1, *d, *a);
          e = &g_gauge_field[y][nu];
          _su3d_times_su3(tmp2, *e, tmp1);
          
          y = g_iup[z][mu];
          b = &g_gauge_field[y][nu];
          z = g_iup[y][nu];
          c = &g_gauge_field[z][nu];
          _su3_times_su3(tmp1, *b, *c);
          /* tmp1 = c1 * tmp1 */
          _real_times_su3(tmp1,g_rgi_C1,tmp1); /* that is the new thing specific of SF */
          _su3_times_su3_acc((*v), tmp2, tmp1);
          
          /* second contr. starting from x 
            * a b c e^+ d^+
            *
            *   bc
            *   __
            * a| _|e
            *    d
            */
          a = &g_gauge_field[x][nu];
          y = g_iup[x][nu];
          b = &g_gauge_field[y][mu];
          _su3_times_su3(tmp1, *a, *b);
          z = g_iup[y][mu];
          c = &g_gauge_field[z][mu];
          _su3_times_su3(tmp2, tmp1, *c);
        
          y = g_iup[x][mu];
          d = &g_gauge_field[y][mu];
          z = g_iup[y][mu];
          e = &g_gauge_field[z][nu];
          _su3_times_su3(tmp1, *d, *e);
          /* tmp1 = c1 * tmp1 */
          _real_times_su3(tmp1,g_rgi_C1,tmp1); /* that is the new thing specific of SF */
          _su3_times_su3d_acc((*v), tmp2, tmp1);
          
          /* 1 contr. starting idn[x][nu] 
            * d^+ a b c e^+
            *
            *    e
            *    _
            * d|__|c
            *   ab
            */
          y = g_idn[x][nu];
          d = &g_gauge_field[y][nu];
          a = &g_gauge_field[y][mu];
          _su3d_times_su3(tmp1, *d, *a);
          z = g_iup[y][mu];
          b = &g_gauge_field[z][mu];
          _su3_times_su3(tmp2, tmp1, *b);
          
          y = g_iup[z][mu];
          c = &g_gauge_field[y][nu];
          z = g_iup[x][mu];
          e = &g_gauge_field[z][mu];
          _su3_times_su3d(tmp1, *c, *e);
          /* tmp1 = c1 * tmp1 */
          _real_times_su3(tmp1,g_rgi_C1,tmp1); /* that is the new thing specific of SF */
          _su3_times_su3_acc((*v), tmp2, tmp1);
          
          /* 1 contr. starting idn[idn[x][mu]][nu] 
            *  e^+ d^+ a b c
            *
            *  e
            *  _
            *d|__|c
            *  ab
            */
          y = g_idn[x][mu];
          z = g_idn[y][nu];
          d = &g_gauge_field[z][nu];
          a = &g_gauge_field[z][mu];
          _su3d_times_su3(tmp1, *d, *a);
          e = &g_gauge_field[y][mu];
          _su3d_times_su3(tmp2, *e, tmp1);
          
          y = g_idn[x][nu];
          b = &g_gauge_field[y][mu];
          z = g_iup[y][mu];
          c = &g_gauge_field[z][nu];
          _su3_times_su3(tmp1, *b, *c);
          /* tmp1 = c1 * tmp1 */
          _real_times_su3(tmp1,g_rgi_C1,tmp1); /* that is the new thing specific of SF */
          _su3_times_su3_acc((*v), tmp2, tmp1); 
          
          /* 1 contr. starting idn[x][mu] 
            * d^+ a b c e^+
            *
            *  bc
            *  __
            *a|_ |e
            *  d
            */
          y = g_idn[x][mu];
          d = &g_gauge_field[y][mu];
          z = g_iup[y][nu];
          a = &g_gauge_field[y][nu];
          _su3d_times_su3(tmp1, *d, *a);
          b = &g_gauge_field[z][mu];
          _su3_times_su3(tmp2, tmp1, *b);
          
          y = g_iup[x][mu];
          e = &g_gauge_field[y][nu];
          z = g_iup[x][nu];
          c = &g_gauge_field[z][mu];
          _su3_times_su3d(tmp1, *c, *e);
          /* tmp1 = c1 * tmp1 */
          _real_times_su3(tmp1,g_rgi_C1,tmp1); /* that is the new thing specific of SF */
          _su3_times_su3_acc((*v), tmp2, tmp1);
        }
      }
      else if (g_t[x] == (g_Tbsf-1)) {
        if (mu == 0) {
          /* first contr. starting from x 
            * a b c e^+ d^+
            *   c
            *   _
            * b| |e
            * a| |d
            */
          a = &g_gauge_field[x][nu];
          y = g_iup[x][nu];
          b = &g_gauge_field[y][nu];
          _su3_times_su3(tmp1, *a, *b);
          z = g_iup[y][nu];
          c = &g_gauge_field[z][mu];
          _su3_times_su3(tmp2, tmp1, *c);
          
          y = g_iup[x][mu];
          d = &g_gauge_field[y][nu];
          z = g_iup[y][nu];
          e = &g_gauge_field[z][nu];
          _su3_times_su3(tmp1, *d, *e);
          /* tmp1 = c1tss * tmp1 */
          _real_times_su3(tmp1,g_C1tss,tmp1); /* that is the new thing specific of SF */
          _su3_times_su3d_acc((*v), tmp2, tmp1);
          
          /* 1 contr. starting idn[idn[x][nu]][nu] 
            * e^+ d^+ a b c
            *
            *e| |c
            *d|_|b
            *  a
            */
          y = g_idn[x][nu];
          z = g_idn[y][nu];
          d = &g_gauge_field[z][nu];
          a = &g_gauge_field[z][mu];
          _su3d_times_su3(tmp1, *d, *a);
          e = &g_gauge_field[y][nu];
          _su3d_times_su3(tmp2, *e, tmp1);
          
          y = g_iup[z][mu];
          b = &g_gauge_field[y][nu];
          z = g_iup[y][nu];
          c = &g_gauge_field[z][nu];
          _su3_times_su3(tmp1, *b, *c);
          /* tmp1 = c1tss * tmp1 */
          _real_times_su3(tmp1,g_C1tss,tmp1); /* that is the new thing specific of SF */
          _su3_times_su3_acc((*v), tmp2, tmp1);
          
          /* 1 contr. starting idn[idn[x][mu]][nu] 
            *  e^+ d^+ a b c
            *
            *  e
            *  _
            *d|__|c
            *  ab
            */
          y = g_idn[x][mu];
          z = g_idn[y][nu];
          d = &g_gauge_field[z][nu];
          a = &g_gauge_field[z][mu];
          _su3d_times_su3(tmp1, *d, *a);
          e = &g_gauge_field[y][mu];
          _su3d_times_su3(tmp2, *e, tmp1);
          
          y = g_idn[x][nu];
          b = &g_gauge_field[y][mu];
          z = g_iup[y][mu];
          c = &g_gauge_field[z][nu];
          _su3_times_su3(tmp1, *b, *c);
          /* tmp1 = c1tts * tmp1 */
          _real_times_su3(tmp1,g_C1tts,tmp1); /* that is the new thing specific of SF */
          _su3_times_su3_acc((*v), tmp2, tmp1); 
          
          /* 1 contr. starting idn[x][mu] 
            * d^+ a b c e^+
            *
            *  bc
            *  __
            *a|_ |e
            *  d
            */
          y = g_idn[x][mu];
          d = &g_gauge_field[y][mu];
          z = g_iup[y][nu];
          a = &g_gauge_field[y][nu];
          _su3d_times_su3(tmp1, *d, *a);
          b = &g_gauge_field[z][mu];
          _su3_times_su3(tmp2, tmp1, *b);
          
          y = g_iup[x][mu];
          e = &g_gauge_field[y][nu];
          z = g_iup[x][nu];
          c = &g_gauge_field[z][mu];
          _su3_times_su3d(tmp1, *c, *e);
          /* tmp1 = c1tts * tmp1 */
          _real_times_su3(tmp1,g_C1tts,tmp1); /* that is the new thing specific of SF */
          _su3_times_su3_acc((*v), tmp2, tmp1);
        }
        else if (mu != 0 && nu == 0) {
          
          /* 1 contr. starting idn[idn[x][nu]][nu] 
            * e^+ d^+ a b c
            *
            *e| |c
            *d|_|b
            *  a
            */
          y = g_idn[x][nu];
          z = g_idn[y][nu];
          d = &g_gauge_field[z][nu];
          a = &g_gauge_field[z][mu];
          _su3d_times_su3(tmp1, *d, *a);
          e = &g_gauge_field[y][nu];
          _su3d_times_su3(tmp2, *e, tmp1);
          
          y = g_iup[z][mu];
          b = &g_gauge_field[y][nu];
          z = g_iup[y][nu];
          c = &g_gauge_field[z][nu];
          _su3_times_su3(tmp1, *b, *c);
          /* tmp1 = c1 * tmp1 */
          _real_times_su3(tmp1,g_rgi_C1,tmp1); /* that is the new thing specific of SF */
          _su3_times_su3_acc((*v), tmp2, tmp1);
          
          /* second contr. starting from x 
            * a b c e^+ d^+
            *
            *   bc
            *   __
            * a| _|e
            *    d
            */
          a = &g_gauge_field[x][nu];
          y = g_iup[x][nu];
          b = &g_gauge_field[y][mu];
          _su3_times_su3(tmp1, *a, *b);
          z = g_iup[y][mu];
          c = &g_gauge_field[z][mu];
          _su3_times_su3(tmp2, tmp1, *c);
          
          y = g_iup[x][mu];
          d = &g_gauge_field[y][mu];
          z = g_iup[y][mu];
          e = &g_gauge_field[z][nu];
          _su3_times_su3(tmp1, *d, *e);
          /* tmp1 = c1tss * tmp1 */
          _real_times_su3(tmp1,g_C1tss,tmp1); /* that is the new thing specific of SF */
          _su3_times_su3d_acc((*v), tmp2, tmp1);
          
          /* 1 contr. starting idn[x][nu] 
            * d^+ a b c e^+
            *
            *    e
            *    _
            * d|__|c
            *   ab
            */
          y = g_idn[x][nu];
          d = &g_gauge_field[y][nu];
          a = &g_gauge_field[y][mu];
          _su3d_times_su3(tmp1, *d, *a);
          z = g_iup[y][mu];
          b = &g_gauge_field[z][mu];
          _su3_times_su3(tmp2, tmp1, *b);
          
          y = g_iup[z][mu];
          c = &g_gauge_field[y][nu];
          z = g_iup[x][mu];
          e = &g_gauge_field[z][mu];
          _su3_times_su3d(tmp1, *c, *e);
          /* tmp1 = c1 * tmp1 */
          _real_times_su3(tmp1,g_rgi_C1,tmp1); /* that is the new thing specific of SF */
          _su3_times_su3_acc((*v), tmp2, tmp1);
          
          /* 1 contr. starting idn[idn[x][mu]][nu] 
            *  e^+ d^+ a b c
            *
            *  e
            *  _
            *d|__|c
            *  ab
            */
          y = g_idn[x][mu];
          z = g_idn[y][nu];
          d = &g_gauge_field[z][nu];
          a = &g_gauge_field[z][mu];
          _su3d_times_su3(tmp1, *d, *a);
          e = &g_gauge_field[y][mu];
          _su3d_times_su3(tmp2, *e, tmp1);
          
          y = g_idn[x][nu];
          b = &g_gauge_field[y][mu];
          z = g_iup[y][mu];
          c = &g_gauge_field[z][nu];
          _su3_times_su3(tmp1, *b, *c);
          /* tmp1 = c1 * tmp1 */
          _real_times_su3(tmp1,g_rgi_C1,tmp1); /* that is the new thing specific of SF */
          _su3_times_su3_acc((*v), tmp2, tmp1); 
          
          /* 1 contr. starting idn[x][mu] 
            * d^+ a b c e^+
            *
            *  bc
            *  __
            *a|_ |e
            *  d
            */
          y = g_idn[x][mu];
          d = &g_gauge_field[y][mu];
          z = g_iup[y][nu];
          a = &g_gauge_field[y][nu];
          _su3d_times_su3(tmp1, *d, *a);
          b = &g_gauge_field[z][mu];
          _su3_times_su3(tmp2, tmp1, *b);
          
          y = g_iup[x][mu];
          e = &g_gauge_field[y][nu];
          z = g_iup[x][nu];
          c = &g_gauge_field[z][mu];
          _su3_times_su3d(tmp1, *c, *e);
          /* tmp1 = c1tss * tmp1 */
          _real_times_su3(tmp1,g_C1tss,tmp1); /* that is the new thing specific of SF */
          _su3_times_su3_acc((*v), tmp2, tmp1);
        }
        else if (mu != 0 && nu != 0) {
          /* first contr. starting from x 
            * a b c e^+ d^+
        *   c
        *   _
        * b| |e
        * a| |d
        */
          a = &g_gauge_field[x][nu];
          y = g_iup[x][nu];
          b = &g_gauge_field[y][nu];
          _su3_times_su3(tmp1, *a, *b);
          z = g_iup[y][nu];
          c = &g_gauge_field[z][mu];
          _su3_times_su3(tmp2, tmp1, *c);
          
          y = g_iup[x][mu];
          d = &g_gauge_field[y][nu];
          z = g_iup[y][nu];
          e = &g_gauge_field[z][nu];
          _su3_times_su3(tmp1, *d, *e);
        /* tmp1 = c1 * tmp1 */
          _real_times_su3(tmp1,g_rgi_C1,tmp1); /* that is the new thing specific of SF */
          _su3_times_su3d_acc((*v), tmp2, tmp1);
          
          /* 1 contr. starting idn[idn[x][nu]][nu] 
            * e^+ d^+ a b c
            *
            *e| |c
            *d|_|b
            *  a
            */
          y = g_idn[x][nu];
          z = g_idn[y][nu];
          d = &g_gauge_field[z][nu];
          a = &g_gauge_field[z][mu];
          _su3d_times_su3(tmp1, *d, *a);
          e = &g_gauge_field[y][nu];
          _su3d_times_su3(tmp2, *e, tmp1);
          
          y = g_iup[z][mu];
          b = &g_gauge_field[y][nu];
          z = g_iup[y][nu];
          c = &g_gauge_field[z][nu];
          _su3_times_su3(tmp1, *b, *c);
          /* tmp1 = c1 * tmp1 */
          _real_times_su3(tmp1,g_rgi_C1,tmp1); /* that is the new thing specific of SF */
          _su3_times_su3_acc((*v), tmp2, tmp1);
          
          /* second contr. starting from x 
            * a b c e^+ d^+
            *
            *   bc
            *   __
            * a| _|e
            *    d
            */
          a = &g_gauge_field[x][nu];
          y = g_iup[x][nu];
          b = &g_gauge_field[y][mu];
          _su3_times_su3(tmp1, *a, *b);
          z = g_iup[y][mu];
          c = &g_gauge_field[z][mu];
          _su3_times_su3(tmp2, tmp1, *c);
          
          y = g_iup[x][mu];
          d = &g_gauge_field[y][mu];
          z = g_iup[y][mu];
          e = &g_gauge_field[z][nu];
          _su3_times_su3(tmp1, *d, *e);
          /* tmp1 = c1 * tmp1 */
          _real_times_su3(tmp1,g_rgi_C1,tmp1); /* that is the new thing specific of SF */
          _su3_times_su3d_acc((*v), tmp2, tmp1);
        
          /* 1 contr. starting idn[x][nu] 
            * d^+ a b c e^+
            *
            *    e
            *    _
            * d|__|c
            *   ab
            */
          y = g_idn[x][nu];
          d = &g_gauge_field[y][nu];
          a = &g_gauge_field[y][mu];
          _su3d_times_su3(tmp1, *d, *a);
          z = g_iup[y][mu];
          b = &g_gauge_field[z][mu];
          _su3_times_su3(tmp2, tmp1, *b);
          
          y = g_iup[z][mu];
          c = &g_gauge_field[y][nu];
          z = g_iup[x][mu];
          e = &g_gauge_field[z][mu];
          _su3_times_su3d(tmp1, *c, *e);
          /* tmp1 = c1 * tmp1 */
          _real_times_su3(tmp1,g_rgi_C1,tmp1); /* that is the new thing specific of SF */
          _su3_times_su3_acc((*v), tmp2, tmp1);
          
          /* 1 contr. starting idn[idn[x][mu]][nu] 
            *  e^+ d^+ a b c
            *
            *  e
            *  _
            *d|__|c
            *  ab
            */
          y = g_idn[x][mu];
          z = g_idn[y][nu];
          d = &g_gauge_field[z][nu];
          a = &g_gauge_field[z][mu];
          _su3d_times_su3(tmp1, *d, *a);
          e = &g_gauge_field[y][mu];
          _su3d_times_su3(tmp2, *e, tmp1);
          
          y = g_idn[x][nu];
          b = &g_gauge_field[y][mu];
          z = g_iup[y][mu];
          c = &g_gauge_field[z][nu];
          _su3_times_su3(tmp1, *b, *c);
          /* tmp1 = c1 * tmp1 */
          _real_times_su3(tmp1,g_rgi_C1,tmp1); /* that is the new thing specific of SF */
          _su3_times_su3_acc((*v), tmp2, tmp1); 
          
          /* 1 contr. starting idn[x][mu] 
            * d^+ a b c e^+
            *
            *  bc
            *  __
            *a|_ |e
            *  d
            */
          y = g_idn[x][mu];
          d = &g_gauge_field[y][mu];
          z = g_iup[y][nu];
          a = &g_gauge_field[y][nu];
          _su3d_times_su3(tmp1, *d, *a);
          b = &g_gauge_field[z][mu];
          _su3_times_su3(tmp2, tmp1, *b);
          
          y = g_iup[x][mu];
          e = &g_gauge_field[y][nu];
          z = g_iup[x][nu];
          c = &g_gauge_field[z][mu];
          _su3_times_su3d(tmp1, *c, *e);
          /* tmp1 = c1 * tmp1 */
          _real_times_su3(tmp1,g_rgi_C1,tmp1); /* that is the new thing specific of SF */
          _su3_times_su3_acc((*v), tmp2, tmp1);
        }
      }
      else if (g_t[x] == 0) {
        /* first contr. starting from x 
          * a b c e^+ d^+
          *   c
          *   _
          * b| |e
          * a| |d
          */
        a = &g_gauge_field[x][nu];
        y = g_iup[x][nu];
        b = &g_gauge_field[y][nu];
        _su3_times_su3(tmp1, *a, *b);
        z = g_iup[y][nu];
        c = &g_gauge_field[z][mu];
        _su3_times_su3(tmp2, tmp1, *c);
        
        y = g_iup[x][mu];
        d = &g_gauge_field[y][nu];
        z = g_iup[y][nu];
        e = &g_gauge_field[z][nu];
        _su3_times_su3(tmp1, *d, *e);
        /* tmp1 = c1tss * tmp1 */
        _real_times_su3(tmp1,g_C1tss,tmp1); /* that is the new thing specific of SF */
        _su3_times_su3d_acc((*v), tmp2, tmp1);
        
        /* 1 contr. starting idn[idn[x][nu]][nu] 
          * e^+ d^+ a b c
          *
          *e| |c
          *d|_|b
          *  a
          */
        y = g_idn[x][nu];
        z = g_idn[y][nu];
        d = &g_gauge_field[z][nu];
        a = &g_gauge_field[z][mu];
        _su3d_times_su3(tmp1, *d, *a);
        e = &g_gauge_field[y][nu];
        _su3d_times_su3(tmp2, *e, tmp1);

        y = g_iup[z][mu];
        b = &g_gauge_field[y][nu];
        z = g_iup[y][nu];
        c = &g_gauge_field[z][nu];
        _su3_times_su3(tmp1, *b, *c);
        /* tmp1 = c1tss * tmp1 */
        _real_times_su3(tmp1,g_C1tss,tmp1); /* that is the new thing specific of SF */
        _su3_times_su3_acc((*v), tmp2, tmp1);
        
        /* second contr. starting from x 
          * a b c e^+ d^+
          *
          *   bc
          *   __
          * a| _|e
          *    d
          */
        a = &g_gauge_field[x][nu];
        y = g_iup[x][nu];
        b = &g_gauge_field[y][mu];
        _su3_times_su3(tmp1, *a, *b);
        z = g_iup[y][mu];
        c = &g_gauge_field[z][mu];
        _su3_times_su3(tmp2, tmp1, *c);
        
        y = g_iup[x][mu];
        d = &g_gauge_field[y][mu];
        z = g_iup[y][mu];
        e = &g_gauge_field[z][nu];
        _su3_times_su3(tmp1, *d, *e);
        /* tmp1 = c1tts * tmp1 */
        _real_times_su3(tmp1,g_C1tts,tmp1); /* that is the new thing specific of SF */
        _su3_times_su3d_acc((*v), tmp2, tmp1);
        
        /* 1 contr. starting idn[x][nu] 
          * d^+ a b c e^+
          *
          *    e
          *    _
          * d|__|c
          *   ab
          */
        y = g_idn[x][nu];
        d = &g_gauge_field[y][nu];
        a = &g_gauge_field[y][mu];
        _su3d_times_su3(tmp1, *d, *a);
        z = g_iup[y][mu];
        b = &g_gauge_field[z][mu];
        _su3_times_su3(tmp2, tmp1, *b);
        
        y = g_iup[z][mu];
        c = &g_gauge_field[y][nu];
        z = g_iup[x][mu];
        e = &g_gauge_field[z][mu];
        _su3_times_su3d(tmp1, *c, *e);
        /* tmp1 = c1tts * tmp1 */
        _real_times_su3(tmp1,g_C1tts,tmp1); /* that is the new thing specific of SF */
        _su3_times_su3_acc((*v), tmp2, tmp1);

      }
    }
  }
#ifdef _KOJAK_INST
#pragma pomp inst end(rectstaples)
#endif
}

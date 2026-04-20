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
#include <tmlqcd_config.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include "get_rectangle_staples.h"
#include "global.h"
#include "su3.h"
#include "ptbc.h"

void get_rectangle_staples(su3 *const v, const int x, const int mu) {
  get_rectangle_staples_general(v, x, mu, (const su3 *const *const)g_gauge_field);
}

void get_rectangle_staples_general(su3 *const v, const int x, const int mu,
                                   const su3 *const *const gf) {
  su3 ALIGN tmp1, tmp2;
  const su3 *a, *b, *c, *d, *e;
  _su3_zero((*v));
  for (int nu = 0; nu < 4; nu++) {
    if (mu != nu) {
      int y, z;
      /* first contr. starting from x
       * a b c e^+ d^+
       *   c
       *   _
       * b| |e
       * a| |d
       */
      a = &gf[x][nu];
      y = g_iup[x][nu];
      b = &gf[y][nu];
      _su3_times_su3(tmp1, *a, *b);
      z = g_iup[y][nu];
      c = &gf[z][mu];
      double const ptbc_fac0 = get_ptbc_coeff(x, nu) * get_ptbc_coeff(y, nu) * get_ptbc_coeff(z, mu);
      _su3_times_su3(tmp2, tmp1, *c);

      y = g_iup[x][mu];
      d = &gf[y][nu];
      z = g_iup[y][nu];
      e = &gf[z][nu];
      double const ptbc_fac1 = get_ptbc_coeff(y, nu) * get_ptbc_coeff(z, nu);
      _su3_times_su3(tmp1, *d, *e);
      //_su3_times_su3d_acc((*v), tmp2, tmp1);
      _real_times_su3_times_su3d_acc((*v), tmp2, tmp1, ptbc_fac0 * ptbc_fac1);

      /* 1 contr. starting idn[idn[x][nu]][nu]
       * e^+ d^+ a b c
       *
       *e| |c
       *d|_|b
       *  a
       */
      y = g_idn[x][nu];
      z = g_idn[y][nu];
      d = &gf[z][nu];
      a = &gf[z][mu];
      _su3d_times_su3(tmp1, *d, *a);
      e = &gf[y][nu];
      double const ptbc_fac2 = get_ptbc_coeff(z, nu) * get_ptbc_coeff(z, mu) * get_ptbc_coeff(y, nu);
      _su3d_times_su3(tmp2, *e, tmp1);

      y = g_iup[z][mu];
      b = &gf[y][nu];
      z = g_iup[y][nu];
      c = &gf[z][nu];
      double const ptbc_fac3 = get_ptbc_coeff(y, nu) * get_ptbc_coeff(z, nu);
      _su3_times_su3(tmp1, *b, *c);
      //_su3_times_su3_acc((*v), tmp2, tmp1);
      _real_times_su3_times_su3_acc((*v), tmp2, tmp1, ptbc_fac2 * ptbc_fac3);

      /* second contr. starting from x
       * a b c e^+ d^+
       *
       *   bc
       *   __
       * a| _|e
       *    d
       */
      a = &gf[x][nu];
      y = g_iup[x][nu];
      b = &gf[y][mu];
      _su3_times_su3(tmp1, *a, *b);
      z = g_iup[y][mu];
      c = &gf[z][mu];
      double const ptbc_fac4 = get_ptbc_coeff(x, nu) * get_ptbc_coeff(y, mu) * get_ptbc_coeff(z, mu);
      _su3_times_su3(tmp2, tmp1, *c);

      y = g_iup[x][mu];
      d = &gf[y][mu];
      z = g_iup[y][mu];
      e = &gf[z][nu];
      double const ptbc_fac5 = get_ptbc_coeff(y, mu) * get_ptbc_coeff(z, nu);
      _su3_times_su3(tmp1, *d, *e);
      //_su3_times_su3d_acc((*v), tmp2, tmp1);
      _real_times_su3_times_su3d_acc((*v), tmp2, tmp1, ptbc_fac4 * ptbc_fac5);

      /* 1 contr. starting idn[x][nu]
       * d^+ a b c e^+
       *
       *    e
       *    _
       * d|__|c
       *   ab
       */
      y = g_idn[x][nu];
      d = &gf[y][nu];
      a = &gf[y][mu];
      _su3d_times_su3(tmp1, *d, *a);
      z = g_iup[y][mu];
      b = &gf[z][mu];
      double const ptbc_fac6 = get_ptbc_coeff(y, nu) * get_ptbc_coeff(y, mu) * get_ptbc_coeff(z, mu);
      _su3_times_su3(tmp2, tmp1, *b);

      y = g_iup[z][mu];
      c = &gf[y][nu];
      z = g_iup[x][mu];
      e = &gf[z][mu];
      double const ptbc_fac7 = get_ptbc_coeff(y, nu) * get_ptbc_coeff(z, mu);
      _su3_times_su3d(tmp1, *c, *e);
      //_su3_times_su3_acc((*v), tmp2, tmp1);
      _real_times_su3_times_su3_acc((*v), tmp2, tmp1, ptbc_fac6 * ptbc_fac7);

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
      d = &gf[z][nu];
      a = &gf[z][mu];
      _su3d_times_su3(tmp1, *d, *a);
      e = &gf[y][mu];
      double const ptbc_fac8 = get_ptbc_coeff(z, nu) * get_ptbc_coeff(z, mu) * get_ptbc_coeff(y, mu);
      _su3d_times_su3(tmp2, *e, tmp1);

      y = g_idn[x][nu];
      b = &gf[y][mu];
      z = g_iup[y][mu];
      c = &gf[z][nu];
      double const ptbc_fac9 = get_ptbc_coeff(y, mu) * get_ptbc_coeff(z, nu);
      _su3_times_su3(tmp1, *b, *c);
      //_su3_times_su3_acc((*v), tmp2, tmp1);
      _real_times_su3_times_su3_acc((*v), tmp2, tmp1, ptbc_fac8 * ptbc_fac9);

      /* 1 contr. starting idn[x][mu]
       * d^+ a b c e^+
       *
       *  bc
       *  __
       *a|_ |e
       *  d
       */
      y = g_idn[x][mu];
      d = &gf[y][mu];
      z = g_iup[y][nu];
      a = &gf[y][nu];
      _su3d_times_su3(tmp1, *d, *a);
      b = &gf[z][mu];
      double const ptbc_fac10 = get_ptbc_coeff(y, mu) * get_ptbc_coeff(y, nu) * get_ptbc_coeff(z, mu);
      _su3_times_su3(tmp2, tmp1, *b);

      y = g_iup[x][mu];
      e = &gf[y][nu];
      z = g_iup[x][nu];
      c = &gf[z][mu];
      double const ptbc_fac11 = get_ptbc_coeff(y, nu) * get_ptbc_coeff(z, mu);
      _su3_times_su3d(tmp1, *c, *e);
      //_su3_times_su3_acc((*v), tmp2, tmp1);
      _real_times_su3_times_su3_acc((*v), tmp2, tmp1, ptbc_fac10 * ptbc_fac11);
    }
  }
}

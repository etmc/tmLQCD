/***********************************************************************
 * Copyright (C) 2001 Martin Hasenbusch
 * Copyright (C) 2006,2007,2008 Carsten Urbach
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
 * This function defines the boundary cond.
 * with arbitrary angle in all directions
 *
 ***********************************************************************/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <math.h>
#include "global.h"
#include "su3.h"
#include "boundary.h"
  
complex ka0, ka1, ka2, ka3;
complex phase_0, phase_1, phase_2, phase_3;
const double PI_ = 3.14159265358979;
double X0, X1, X2, X3;

void boundary(const double kappa) {
  double x0,x1,x2,x3;
  x0 = X0 * PI_/((T)*g_nproc_t);
  x1 = X1 * PI_/((LX)*g_nproc_x);
  x2 = X2 * PI_/((LY)*g_nproc_y);
  x3 = X3 * PI_/((LZ)*g_nproc_z);
  ka0.re = kappa * cos(x0); 
  ka0.im = kappa * sin(x0);
  ka1.re = kappa * cos(x1); 
  ka1.im = kappa * sin(x1);
  ka2.re = kappa * cos(x2); 
  ka2.im = kappa * sin(x2);
  ka3.re = kappa * cos(x3); 
  ka3.im = kappa * sin(x3);
  phase_0.re = -ka0.re;
  phase_1.re = -ka1.re;
  phase_2.re = -ka2.re;
  phase_3.re = -ka3.re;
  phase_0.im = -ka0.im;
  phase_1.im = -ka1.im;
  phase_2.im = -ka2.im;
  phase_3.im = -ka3.im;
}

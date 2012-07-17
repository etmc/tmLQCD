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
#include <complex.h>
#include "global.h"
#include "su3.h"
#include "boundary.h"
  
_Complex double ALIGN ka0, ka1, ka2, ka3;
_Complex double ALIGN phase_0, phase_1, phase_2, phase_3;
const double PI_ = 3.14159265358979;
double X0, X1, X2, X3;

void boundary(const double kappa)
{
  double x0,x1,x2,x3;
  x0 = X0 * PI_/((T)*g_nproc_t);
  x1 = X1 * PI_/((LX)*g_nproc_x);
  x2 = X2 * PI_/((LY)*g_nproc_y);
  x3 = X3 * PI_/((LZ)*g_nproc_z);
  ka0 = kappa * cexp(x0 * I);
  ka1 = kappa * cexp(x1 * I);
  ka2 = kappa * cexp(x2 * I);
  ka3 = kappa * cexp(x3 * I);
  phase_0 = -ka0;
  phase_1 = -ka1;
  phase_2 = -ka2;
  phase_3 = -ka3;
}

/***********************************************************************
 *
 * Copyright (C) 2011 Elena Garcia-Ramos
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

#ifndef _P_M_ETA_H
#define _P_M_ETA_H

#include "su3.h"

extern int x_n_cheby;
extern double * x_cheby_coef;

void norm_X_sqr_psi(spinor * const R, spinor * const S, double const mstar);

void norm_X_n_psi(spinor * const R, spinor * const S, const int n, double const mstar);

void X_over_sqrt_X_sqr(spinor * const R, double * const c, const int n, spinor * const S, const double minev, double const mstar);

void h_X_sqr_eta(spinor * const R1,spinor * const R2,spinor * const S, double const mstar);

void h_X_eta(spinor * const R,spinor * const S, double const mstar);

void h_X_4_eta(spinor * const R1, spinor * const R2, spinor * const S, double const mstar);

void Check_Approximation(double const mstar, const int repro);

#endif


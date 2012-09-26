/***********************************************************************
 *
 * Copyright (C) 2005 Martin Hasenbusch
 *               2009 Carsten Urbach
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

#ifndef _CLOVER_H
#define _CLOVER_H

#include "su3.h"

extern su3 *** sw;
extern su3 *** sw_inv;
extern su3 ** swm, ** swp;

void assign_mul_one_sw_pm_imu(const int ieo, spinor * const k, const spinor * const l, const double mu);
void assign_mul_one_sw_pm_imu_inv(const int ieo, spinor * const k, const spinor * const l, const double mu);
void Msw_full(spinor * const Even_new, spinor * const Odd_new, 
	      spinor * const Even, spinor * const Odd);
void clover_inv(const int ieo, spinor * const l, const double mu);
void Qsw_psi(spinor * const l, spinor * const k);
void Qsw_plus_psi(spinor * const l, spinor * const k);
void Qsw_minus_psi(spinor * const l, spinor * const k);
void Qsw_sq_psi(spinor * const l, spinor * const k);
void Qsw_pm_psi(spinor * const l, spinor * const k);
void Msw_psi(spinor * const l, spinor * const k);
void Msw_plus_psi(spinor * const l, spinor * const k);
void Msw_minus_psi(spinor * const l, spinor * const k);
void H_eo_sw_inv_psi(spinor * const l, spinor * const k, const int ieo, const double mu);
void init_sw_fields();

#endif

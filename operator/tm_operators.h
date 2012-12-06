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

#ifndef _TM_OPERATORS_H
#define _TM_OPERATORS_H

#include "su3.h"

/* This is the full matrix multiplication */
void M_full(spinor * const Even_new, spinor * const Odd_new, 
	    spinor * const Even, spinor * const Odd);
void Q_full(spinor * const Even_new, spinor * const Odd_new, 
	    spinor * const Even, spinor * const Odd);
void M_minus_1_timesC(spinor * const Even_new, spinor * const Odd_new, 
		      spinor * const Even, spinor * const Odd);

void Qtm_plus_psi(spinor * const l, spinor * const k);
void Qtm_plus_psi_nocom(spinor * const l, spinor * const k);
void Qtm_minus_psi(spinor * const l, spinor * const k);
void Mtm_plus_psi(spinor * const l, spinor * const k);
void Mtm_plus_psi_nocom(spinor * const l, spinor * const k);
void Mtm_minus_psi(spinor * const l, spinor * const k);
void Qtm_pm_psi(spinor * const l, spinor * const k);
void Qtm_pm_psi_nocom(spinor * const l, spinor * const k);
void H_eo_tm_inv_psi(spinor * const l, spinor * const k, const int ieo, const double sign);
void mul_one_pm_imu_inv(spinor * const l, const double _sign, const int N);
void assign_mul_one_pm_imu_inv(spinor * const l, spinor * const k, const double _sign, const int N);
void assign_mul_one_pm_imu(spinor * const l, spinor * const k, const double _sign, const int N);
void mul_one_pm_imu(spinor * const l, const double _sign);
void mul_one_pm_imu_sub_mul(spinor * const l, spinor * const k,
			    spinor * const j, const double _sign, const int N);

void Qtm_plus_sym_psi(spinor * const l, spinor * const k);
void Qtm_plus_sym_psi_nocom(spinor * const l, spinor * const k);
void Qtm_minus_sym_psi(spinor * const l, spinor * const k);
void Mtm_plus_sym_psi(spinor * const l, spinor * const k);
void Mtm_minus_sym_psi(spinor * const l, spinor * const k);
void Mtm_plus_sym_psi_nocom(spinor * const l, spinor * const k);
void Mtm_minus_sym_psi_nocom(spinor * const l, spinor * const k);
void Qtm_pm_sym_psi(spinor * const l, spinor * const k);

void Q_pm_psi(spinor * const l, spinor * const k);
void Q_pm_psi_prec(spinor * const l, spinor * const k);
void Q_pm_psi_gpu(spinor * const l, spinor * const k);
void Q_pm_psi2(spinor * const l, spinor * const k);
void Q_minus_psi(spinor * const l, spinor * const k);
void Q_minus_psi_gpu(spinor * const l, spinor * const k);
void Q_plus_psi(spinor * const l, spinor * const k);

#endif

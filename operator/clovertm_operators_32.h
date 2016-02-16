/***********************************************************************
 *
 * Copyright (C) 2005 Martin Hasenbusch
 *               2009 Carsten Urbach
 *               2012 Carsten Urbach
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

#ifndef _CLOVERTM_OPERATORS_32_H
#define _CLOVERTM_OPERATORS_32_H

#include "su3.h"

extern su3 *** sw;
extern su3 *** sw_inv;
extern su3_32 *** sw_32;
extern su3_32 *** sw_inv_32;
extern su3 ** swm, ** swp;

void clover_inv_32_orphaned(spinor32 * const l, const int tau3sign, const double mu);
void clover_inv_32(spinor32 * const l, const int tau3sign, const double mu);
void Qsw_pm_psi_32(spinor32 * const l, spinor32 * const k);
void clover_gamma5_32_orphaned(const int ieo, 
		   spinor32 * const l, const spinor32 * const k, const spinor32 * const j,
		   const double mu);
void clover_gamma5_32(const int ieo, 
		   spinor32 * const l, const spinor32 * const k, const spinor32 * const j,
		   const double mu);

void assign_mul_one_sw_pm_imu_eps_32(const int ieo, 
          spinor32 * const k_s, spinor32 * const k_c, 
          const spinor32 * const l_s, const spinor32 * const l_c,
          const float mu, const float eps);
void assign_mul_one_sw_pm_imu_eps_32_orphaned(const int ieo, 
          spinor32 * const k_s, spinor32 * const k_c, 
          const spinor32 * const l_s, const spinor32 * const l_c,
          const float mu, const float eps);

void clover_gamma5_nd_32(const int ieo,
          spinor32 * const l_c, spinor32 * const l_s,
          const spinor32 * const k_c, const spinor32 * const k_s,
          const spinor32 * const j_c, const spinor32 * const j_s,
          const float mubar, const float epsbar);
void clover_gamma5_nd_32_orphaned(const int ieo,
          spinor32 * const l_c, spinor32 * const l_s,
          const spinor32 * const k_c, const spinor32 * const k_s,
          const spinor32 * const j_c, const spinor32 * const j_s,
          const float mubar, const float epsbar);

void clover_inv_nd_32(const int ieo, spinor32 * const l_c, spinor32 * const l_s);
void clover_inv_nd_32_orphaned(const int ieo, spinor32 * const l_c, spinor32 * const l_s);

#endif


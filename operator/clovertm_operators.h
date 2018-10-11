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

#ifndef _CLOVERTM_OPERATORS_H
#define _CLOVERTM_OPERATORS_H

#include "su3.h"
#include "block.h"

extern su3 *** sw;
extern su3 *** sw_inv;
extern su3_32 *** sw_32;
extern su3_32 *** sw_inv_32;
extern su3 ** swm, ** swp;

void assign_mul_one_sw_pm_imu_site_lexic(const int ix, spinor * const k,  const spinor * const l, const double mu);
void assign_mul_one_sw_pm_imu_site_lexic_32(const int ix, spinor32 * const k,  const spinor32 * const l, const float mu);

void Qsw_full(spinor * const Even_new, spinor * const Odd_new,
              spinor * const Even, spinor * const Odd);
void Qsw_full_plus_psi(spinor * const l, spinor * const k);
void Qsw_full_minus_psi(spinor * const l, spinor * const k);
void Qsw_full_pm_psi(spinor * const l, spinor * const k);
void Msw_full_minus_psi(spinor * const l, spinor * const k);

void assign_mul_one_sw_pm_imu(const int ieo, spinor * const k, spinor * const l, const double mu);
void assign_mul_one_sw_pm_imu_32(const int ieo, spinor32 * const k, spinor32 * const l, const float mu);
void assign_mul_one_sw_pm_imu_block(const int ieo, spinor * const k, spinor * const l, const double mu, block *blk);
void assign_mul_one_sw_pm_imu_block_32(const int ieo, spinor32 * const k, spinor32 * const l, const float mu, block *blk);
void assign_mul_one_sw_pm_imu_inv(const int ieo, spinor * const k, spinor * const l, const double mu);
void assign_mul_one_sw_pm_imu_inv_32(const int ieo, spinor32 * const k, spinor32 * const l, const float mu);
void assign_mul_one_sw_pm_imu_inv_block(const int ieo, spinor * const k, spinor * const l, const double mu, block *blk);
void assign_mul_one_sw_pm_imu_inv_block_32(const int ieo, spinor32 * const k, spinor32 * const l, const float mu, block *blk);

void Mee_sw_psi(spinor * const l, spinor * const k, const double mu);
void Mee_sw_inv_psi(spinor * const k, spinor * const l, const double mu);
void Msw_full(spinor * const Even_new, spinor * const Odd_new, 
	      spinor * const Even, spinor * const Odd);

void clover_inv(spinor * const l, const int tau3sign, const double mu);
void Qsw_psi(spinor * const l, spinor * const k);
void Qsw_plus_psi(spinor * const l, spinor * const k);
void Qsw_minus_psi(spinor * const l, spinor * const k);
void Qsw_sq_psi(spinor * const l, spinor * const k);
void Qsw_pm_psi(spinor * const l, spinor * const k);
void Msw_psi(spinor * const l, spinor * const k);
void Msw_plus_psi(spinor * const l, spinor * const k);
void Msw_minus_psi(spinor * const l, spinor * const k);
void H_eo_sw_inv_psi(spinor * const l, spinor * const k, const int ieo, const int tau3sign, const double mu);
void init_sw_fields();
void copy_32_sw_fields();

void clover_nd(const int ieo, 
	       spinor * const l_s, spinor * const l_c, 
	       const spinor * const k_s, const spinor * const k_c, 
	       const spinor * const j_s, const spinor * const j_c,
	       const double mubar, const double epsbar);

void clover_gamma5_nd(const int ieo, 
		      spinor * const l_s, spinor * const l_c, 
		      const spinor * const k_s, const spinor * const k_c, 
		      const spinor * const j_s, const spinor * const j_c,
		      const double mubar, const double epsbar);
void clover_inv_nd(const int ieo, spinor * const l_s, spinor * const l_c);

void assign_mul_one_sw_pm_imu_eps(const int ieo, 
				  spinor * const k_s, spinor * const k_c, 
				  const spinor * const l_s, const spinor * const l_c,
				  const double mu, const double eps);
#endif

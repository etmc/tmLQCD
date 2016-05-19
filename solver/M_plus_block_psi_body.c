/***********************************************************************
 * Copyright (C) 2016 Carsten Urbach
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


void _PSWITCH(Mtm_plus_block_psi)(_PTSWITCH(spinor) * const l, _PTSWITCH(spinor) * const k, const int i) {
  block * blk = &block_list[i];
  int vol = (*blk).volume/2;
  _PSWITCH(Block_H_psi)(blk, &_PTSWITCH(g_spinor_field)[DUM_MATRIX+1][i*vol], k, EO);
  _PSWITCH(mul_one_pm_imu_inv)(&_PTSWITCH(g_spinor_field)[DUM_MATRIX+1][i*vol], +1., vol);
  _PSWITCH(Block_H_psi)(blk, &_PTSWITCH(g_spinor_field)[DUM_MATRIX][i*vol], &_PTSWITCH(g_spinor_field)[DUM_MATRIX+1][i*vol], OE);
  _PSWITCH(mul_one_pm_imu_sub_mul)(l, k, &_PTSWITCH(g_spinor_field)[DUM_MATRIX][i*vol], +1., vol);
  return;
}

void _PSWITCH(Msw_plus_block_psi)(_PTSWITCH(spinor) * l, _PTSWITCH(spinor) *  k, const int i) {
  block * blk = &block_list[i];
  int vol = (*blk).volume/2;
  _PSWITCH(Block_H_psi)(blk, &_PTSWITCH(g_spinor_field)[DUM_MATRIX+1][i*vol], k, EO);
  _PSWITCH(assign_mul_one_sw_pm_imu_inv_block)(EE, &_PTSWITCH(g_spinor_field)[DUM_MATRIX][i*vol],&_PTSWITCH(g_spinor_field)[DUM_MATRIX+1][i*vol], g_mu, blk);
  _PSWITCH(Block_H_psi)(blk, &_PTSWITCH(g_spinor_field)[DUM_MATRIX+1][i*vol], &_PTSWITCH(g_spinor_field)[DUM_MATRIX][i*vol], OE);
  _PSWITCH(assign_mul_one_sw_pm_imu_block)(OO, &_PTSWITCH(g_spinor_field)[DUM_MATRIX][i*vol],k,g_mu,blk);
  _PSWITCH(diff)(l,&_PTSWITCH(g_spinor_field)[DUM_MATRIX][i*vol],&_PTSWITCH(g_spinor_field)[DUM_MATRIX+1][i*vol],vol);
  return;
}

/* $Id$ */

#ifndef _TM_OPERATORS_H
#define _TM_OPERATORS_H

void Qtm_plus_psi(const int l, const int k);
void Qtm_minus_psi(const int l, const int k);
void Mtm_plus_psi(const int l, const int k);
void Mtm_minus_psi(const int l, const int k);
void Qtm_pm_psi(const int l, const int k);
void H_eo_tm_inv_psi(const int l, const int k, const int ieo, const double sign);
void mul_one_pm_imu_inv(const int l, const double _sign);
void assign_mul_one_pm_imu_inv(const int l, const int k, const double _sign);
void assign_mul_one_pm_imu(const int l, const int k, const double _sign);
void mul_one_pm_imu(const int l, const double _sign);

#endif


#ifndef _TM_OPERATORS_32_H
#define _TM_OPERATORS_32_H

void mul_one_pm_imu_inv_32(spinor32 * const l, const float _sign, const int N);
void mul_one_pm_imu_sub_mul_gamma5_32(spinor32 * const l, spinor32 * const k, spinor32 * const j, const float _sign);
void Qtm_pm_psi_32(spinor32 * const l, spinor32 * const k);
void Q_pm_psi_32(spinor32 * const l, spinor32 * const k);

#endif
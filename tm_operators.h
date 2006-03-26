/* $Id$ */

#ifndef _TM_OPERATORS_H
#define _TM_OPERATORS_H

#include "su3.h"

void Qtm_plus_psi(spinor * const l, spinor * const k);
void Qtm_plus_psi_nocom(spinor * const l, spinor * const k);
void Qtm_minus_psi(spinor * const l, spinor * const k);
void Mtm_plus_psi(spinor * const l, spinor * const k);
void Mtm_plus_psi_nocom(spinor * const l, spinor * const k);
void Mtm_minus_psi(spinor * const l, spinor * const k);
void Qtm_pm_psi(spinor * const l, spinor * const k);
void Qtm_pm_psi_nocom(spinor * const l, spinor * const k);
void H_eo_tm_inv_psi(spinor * const l, spinor * const k, const int ieo, const double sign);
void mul_one_pm_imu_inv(spinor * const l, const double _sign);
void assign_mul_one_pm_imu_inv(spinor * const l, spinor * const k, const double _sign);
void assign_mul_one_pm_imu(spinor * const l, spinor * const k, const double _sign);
void mul_one_pm_imu(spinor * const l, const double _sign);

void Qtm_plus_sym_psi(spinor * const l, spinor * const k);
void Qtm_plus_sym_psi_nocom(spinor * const l, spinor * const k);
void Qtm_minus_sym_psi(spinor * const l, spinor * const k);
void Mtm_plus_sym_psi(spinor * const l, spinor * const k);
void Mtm_minus_sym_psi(spinor * const l, spinor * const k);
void Mtm_plus_sym_psi_nocom(spinor * const l, spinor * const k);
void Mtm_minus_sym_psi_nocom(spinor * const l, spinor * const k);
void Qtm_pm_sym_psi(spinor * const l, spinor * const k);


#endif

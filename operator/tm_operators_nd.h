/***********************************************************************
 *
 * Copyright (C) 2006,2007,2008 Karl Jansen, Thomas Chiarappa, 
 *                              Carsten Urbach
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

#ifndef _TM_OPERATTORS_ND_H
#define _TM_OPERATTORS_ND_H

void mul_one_pm_itau2(spinor * const p, spinor * const q,
                      spinor * const r, spinor * const s,
                      const double sign, const int N);

void M_full_ndpsi(spinor * const Even_new_s, spinor * const Odd_new_s, 
                  spinor * const Even_new_c, spinor * const Odd_new_c, 
                  spinor * const Even_s, spinor * const Odd_s,
                  spinor * const Even_c, spinor * const Odd_c);

void Msw_full_ndpsi(spinor * const Even_new_s, spinor * const Odd_new_s, 
                    spinor * const Even_new_c, spinor * const Odd_new_c, 
                    spinor * const Even_s, spinor * const Odd_s,
                    spinor * const Even_c, spinor * const Odd_c);

//This works with tm and tm+clover 
void D_ndpsi(spinor * const l_strange, spinor * const l_charm,
             spinor * const k_strange,  spinor * const k_charm);

void Qtm_ndpsi(spinor * const l_strange, spinor * const l_charm,
               spinor * const k_strange,  spinor * const k_charm);
void Qsw_ndpsi(spinor * const l_strange, spinor * const l_charm,
               spinor * const k_strange, spinor * const k_charm);

void Qtm_tau1_ndpsi_add_Ishift(spinor * const l_strange, spinor * const l_charm,
                               spinor * const k_strange,  spinor * const k_charm);
void Qtm_tau1_ndpsi_sub_Ishift(spinor * const l_strange, spinor * const l_charm,
                               spinor * const k_strange,  spinor * const k_charm);
void Qsw_tau1_ndpsi_add_Ishift(spinor * const l_strange, spinor * const l_charm,
                               spinor * const k_strange,  spinor * const k_charm);
void Qsw_tau1_ndpsi_sub_Ishift(spinor * const l_strange, spinor * const l_charm,
                               spinor * const k_strange,  spinor * const k_charm);


void Qtm_dagger_ndpsi(spinor * const l_strange, spinor * const l_charm,
                      spinor * const k_strange, spinor * const k_charm);
void Qsw_dagger_ndpsi(spinor * const l_strange, spinor * const l_charm,
                      spinor * const k_strange, spinor * const k_charm);

void Qtm_pm_ndpsi(spinor * const l_strange, spinor * const l_charm,
                  spinor * const k_strange, spinor * const k_charm);
void Qtm_pm_ndpsi_shift(spinor * const l_strange, spinor * const l_charm,
                        spinor * const k_strange, spinor * const k_charm);

void Qsw_pm_ndpsi(spinor * const l_strange, spinor * const l_charm,
                  spinor * const k_strange, spinor * const k_charm);
void Qsw_pm_ndpsi_shift(spinor * const l_strange, spinor * const l_charm,
                        spinor * const k_strange, spinor * const k_charm);

void Qtm_pm_ndbipsi(bispinor * const bisp_l, bispinor * const bisp_k);
void Qsw_pm_ndbipsi(bispinor * const bisp_l, bispinor * const bisp_k);

void Q_tau1_sub_const_ndpsi(spinor * const l_strange, spinor * const l_charm,
                            spinor * const k_strange, spinor * const k_charm, 
                            const _Complex double z, const double Cpol, const double invev);
void Qsw_tau1_sub_const_ndpsi(spinor * const l_strange, spinor * const l_charm,
                              spinor * const k_strange, spinor * const k_charm, 
                              const _Complex double z, const double Cpol, const double invev);

void H_eo_tm_ndpsi(spinor * const l_strange, spinor * const l_charm, 
                   spinor * const k_strange, spinor * const k_charm, 
                   const int ieo);
void H_eo_sw_ndpsi(spinor * const l_strange, spinor * const l_charm, 
                   spinor * const k_strange, spinor * const k_charm);


void M_ee_inv_ndpsi(spinor * const l_strange, spinor * const l_charm, 
                    spinor * const k_strange, spinor * const k_charm,
                    const double mu, const double eps);

void Msw_ee_inv_ndpsi(spinor * const l_strange, spinor * const l_charm, 
                      spinor * const k_strange, spinor * const k_charm);

void Q_test_epsilon(spinor * const l_strange, spinor * const l_charm,
                    spinor * const k_strange, spinor * const k_charm);

void Qtau1_P_ndpsi(spinor * const l_strange, spinor * const l_charm,
                   spinor * const k_strange, spinor * const k_charm);

void Qtm_pm_Ptm_pm_psi(spinor * const l, spinor * const k);

void Qtm_pm_sub_const_nrm_psi(spinor * const l, spinor * const k,const _Complex double z);

/* ************************************************
 * for noise reduction 
 * this implements
 * a = B^dagger H b
 * 
 * with Hopping matrix H and
 *
 * B = (1-i\g5\tau^1\musigma-\tau^3\mudelta)/c
 * where
 * c = 1+\musigma^2-\mudelta^2
 *
 **************************************************/

void red_noise_nd(spinor * const lse, spinor * const lso, spinor * const lce, spinor * const lco);


#endif

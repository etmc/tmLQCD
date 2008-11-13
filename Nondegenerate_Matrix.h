/* $Id$ */
#ifndef _NONDEGENRATE_MATRIX_H
#define _NONDEGENRATE_MATRIX_H

void mul_one_minus_imubar(spinor * const l, spinor * const k);
/******************************************
 * mul_one_plus_imubar_inv computes
 * l = [(1-i\mubar\gamma_5) * l
 *
 */

void mul_one_plus_imubar(spinor * const l, spinor * const k);
/******************************************
 * mul_one_plus_imubar_inv computes
 * l = [(1+i\mubar\gamma_5) * l
 *
*/

void mul_one_pm_itau2(spinor * const p, spinor * const q,
		      spinor * const r, spinor * const s,
		      const double sign, const int N);

void QNon_degenerate(spinor * const l_strange, spinor * const l_charm,
                     spinor * const k_strange,  spinor * const k_charm);

void QdaggerNon_degenerate(spinor * const l_strange, spinor * const l_charm,
                           spinor * const k_strange, spinor * const k_charm);

void Q_Qdagger_ND(spinor * const l_strange, spinor * const l_charm,
                  spinor * const k_strange, spinor * const k_charm);

void Q_Qdagger_ND_BI(bispinor * const bisp_l, bispinor * const bisp_k);

void Q_tau1_min_cconst_ND(spinor * const l_strange, spinor * const l_charm,
                       spinor * const k_strange, spinor * const k_charm, 
                       const complex z);

void H_eo_ND(spinor * const l_strange, spinor * const l_charm, 
             spinor * const k_strange, spinor * const k_charm, 
	     const int ieo);

void M_ee_inv_ND(spinor * const l_strange, spinor * const l_charm, 
		 spinor * const k_strange, spinor * const k_charm);

void Q_test_epsilon(spinor * const l_strange, spinor * const l_charm,
                    spinor * const k_strange, spinor * const k_charm);

void Qtau1_P_ND(spinor * const l_strange, spinor * const l_charm,
		spinor * const k_strange, spinor * const k_charm);

void Qtm_pm_Ptm_pm_psi(spinor * const l, spinor * const k);
void Qtm_pm_min_cconst_nrm(spinor * const l, spinor * const k,const complex z);

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

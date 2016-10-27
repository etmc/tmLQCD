/***********************************************************************
 *
 * Copyright (C) 2015 Florian Burger
 * based on the corresponding 64 bit operators in tm_operators_nd.c
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
 *                                                            
 * This file contains operators for twisted mass Wilson QCD   
 * to construct a multiplication with a non-degenerate        
 * flavour matrix                                             
 *                                                            
 *                                                            
 ***********************************************************************/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "global.h"
#include "su3.h"
#include "operator/Hopping_Matrix_32.h"
#include "phmc.h"
#include "gamma.h"
#include "linalg_eo.h"
#include "operator/tm_operators_32.h"
#include "operator/tm_operators_nd.h"
#include "operator/clovertm_operators_32.h"
#include "operator/D_psi.h"
#include "tm_operators_nd_32.h"



void sub_epsbar_tau1_32(spinor32 * const l_strange, spinor32 * const l_charm , spinor32 * const k_strange, spinor32 * const k_charm){
  mul_r_32(g_spinor_field32[2], (float) g_epsbar, k_strange , VOLUME);
  mul_r_32(g_spinor_field32[3], (float) g_epsbar, k_charm, VOLUME);
  diff_32(l_strange, l_strange, g_spinor_field32[3], VOLUME);
  diff_32(l_charm, l_charm, g_spinor_field32[2], VOLUME);  
}


void Q_pm_ndpsi_32(spinor32 * const l_strange, spinor32 * const l_charm, spinor32 * const k_strange, spinor32 * const k_charm)
{

  //D_h^{dagger}
  //tau^1 by s<->c
  
     D_psi_32(l_strange, k_charm);
     g_mu = -g_mu;
     D_psi_32(l_charm, k_strange);     
     g_mu = -g_mu;
     
     sub_epsbar_tau1_32(l_strange, l_charm, k_charm, k_strange);
     
     gamma5_32(g_spinor_field32[0], l_strange, VOLUME);
     gamma5_32(g_spinor_field32[1], l_charm, VOLUME);    
     
    //D_h
    //tau^1 by s<->c     
     D_psi_32(l_strange, g_spinor_field32[1]);
     g_mu = -g_mu;
     D_psi_32(l_charm, g_spinor_field32[0]);         
     g_mu = -g_mu;
     sub_epsbar_tau1_32(l_strange, l_charm, g_spinor_field32[1], g_spinor_field32[0]);
     
     gamma5_32(l_strange, l_strange, VOLUME);      
     gamma5_32(l_charm, l_charm, VOLUME);
     /* At the end, the normalisation by the max. eigenvalue  */ 
     /* Twice  phmc_invmaxev  since we consider here  D Ddag  !!! */
     mul_r_32(l_charm, (float) phmc_invmaxev*phmc_invmaxev, l_charm, VOLUME);
     mul_r_32(l_strange, (float) phmc_invmaxev*phmc_invmaxev, l_strange, VOLUME);     

}

// l_ and k_ are allowed to be the same spinors
void M_ee_inv_ndpsi_32_orphaned(spinor32 * const l_s, spinor32 * const l_c, 
		    spinor32 * const k_s, spinor32 * const k_c,
		    const float mu, const float eps) {
  float nrm = 1./(1.+ mu*mu - eps*eps);
  spinor32 *r_s, *r_c, *s_s, *s_c;
  su3_vector32 ALIGN32 phi1, phi2;

#ifdef TM_USE_OMP
#pragma omp for
#endif
  for(unsigned int ix = 0; ix < (VOLUME/2); ++ix){
    r_s = l_s + ix;
    r_c = l_c + ix;
    s_s = k_s + ix;
    s_c = k_c + ix;

    _complex_times_vector(phi1, (1. - mu * I), s_s->s0);
    _vector_add_mul(phi1, eps, s_c->s0);
    _complex_times_vector(phi2, (1. + mu * I), s_c->s0);
    _vector_add_mul(phi2, eps, s_s->s0);
    _vector_mul(r_s->s0, nrm, phi1);
    _vector_mul(r_c->s0, nrm, phi2);

    _complex_times_vector(phi1, (1. - mu * I), s_s->s1);
    _vector_add_mul(phi1, eps, s_c->s1);
    _complex_times_vector(phi2, (1. + mu * I), s_c->s1);
    _vector_add_mul(phi2, eps, s_s->s1);
    _vector_mul(r_s->s1, nrm, phi1);
    _vector_mul(r_c->s1, nrm, phi2);

    _complex_times_vector(phi1, (1. + mu * I), s_s->s2);
    _vector_add_mul(phi1, eps, s_c->s2);
    _complex_times_vector(phi2, (1. - mu * I), s_c->s2);
    _vector_add_mul(phi2, eps, s_s->s2);
    _vector_mul(r_s->s2, nrm, phi1);
    _vector_mul(r_c->s2, nrm, phi2);

    _complex_times_vector(phi1, (1. + mu * I), s_s->s3);
    _vector_add_mul(phi1, eps, s_c->s3);
    _complex_times_vector(phi2, (1. - mu * I), s_c->s3);
    _vector_add_mul(phi2, eps, s_s->s3);
    _vector_mul(r_s->s3, nrm, phi1);
    _vector_mul(r_c->s3, nrm, phi2);

  }
}

void M_ee_inv_ndpsi_32(spinor32 * const l_s, spinor32 * const l_c, 
		    spinor32 * const k_s, spinor32 * const k_c,
		    const float mu, const float eps) {
#ifdef TM_USE_OMP
#pragma omp parallel
  {
#endif
  M_ee_inv_ndpsi_32_orphaned(l_s, l_c, k_s, k_c, mu, eps);
#ifdef TM_USE_OMP
  } /* OpenMP closing brace */
#endif
  return;
}

// l_ and k_ are allowed to be the same spinors
void M_oo_sub_g5_ndpsi_32_orphaned(spinor32 * const l_s, spinor32 * const l_c, 
		       spinor32 * const k_s, spinor32 * const k_c,
		       spinor32 * const j_s, spinor32 * const j_c,
		       const float mu, const float eps) {
  spinor32 *r_s, *r_c, *s_s, *s_c, *t_s, *t_c;
  su3_vector32 ALIGN32 phi1, phi2;

#ifdef TM_USE_OMP
#pragma omp for
#endif
  for(unsigned int ix = 0; ix < (VOLUME/2); ++ix){
    r_s = l_s + ix;
    r_c = l_c + ix;
    s_s = k_s + ix;
    s_c = k_c + ix;
    t_s = j_s + ix;
    t_c = j_c + ix;

    _complex_times_vector(phi1, (1. - mu * I), s_s->s0);
    _vector_add_mul(phi1, eps, s_c->s0);
    _complex_times_vector(phi2, (1. + mu * I), s_c->s0);
    _vector_add_mul(phi2, eps, s_s->s0);
    _vector_sub(r_s->s0, phi1, t_s->s0);
    _vector_sub(r_c->s0, phi2, t_c->s0);

    _complex_times_vector(phi1, (1. - mu * I), s_s->s1);
    _vector_add_mul(phi1, eps, s_c->s1);
    _complex_times_vector(phi2, (1. + mu * I), s_c->s1);
    _vector_add_mul(phi2, eps, s_s->s1);
    _vector_sub(r_s->s1, phi1, t_s->s1);
    _vector_sub(r_c->s1, phi2, t_c->s1);

    _complex_times_vector(phi1, (1. + mu * I), s_s->s2);
    _vector_add_mul(phi1, eps, s_c->s2);
    _complex_times_vector(phi2, (1. - mu * I), s_c->s2);
    _vector_add_mul(phi2, eps, s_s->s2);
    _vector_sub(r_s->s2, t_s->s2, phi1);
    _vector_sub(r_c->s2, t_c->s2, phi2);

    _complex_times_vector(phi1, (1. + mu * I), s_s->s3);
    _vector_add_mul(phi1, eps, s_c->s3);
    _complex_times_vector(phi2, (1. - mu * I), s_c->s3);
    _vector_add_mul(phi2, eps, s_s->s3);
    _vector_sub(r_s->s3, t_s->s3, phi1);
    _vector_sub(r_c->s3, t_c->s3, phi2);
  }
}

void M_oo_sub_g5_ndpsi_32(spinor32 * const l_s, spinor32 * const l_c, 
		       spinor32 * const k_s, spinor32 * const k_c,
		       spinor32 * const j_s, spinor32 * const j_c,
		       const float mu, const float eps) {
#ifdef TM_USE_OMP
#pragma omp parallel
  {
#endif
  M_oo_sub_g5_ndpsi_32_orphaned(l_s,l_c,k_s,k_c,j_s,j_c,mu,eps);
#ifdef TM_USE_OMP
  } /* OpenMP closing brace */
#endif
  return;
}

void Qtm_pm_ndpsi_32(spinor32 * const l_strange, spinor32 * const l_charm,
		  spinor32 * const k_strange, spinor32 * const k_charm){
#ifdef TM_USE_OMP
#pragma omp parallel
  {
#endif
  /* first the  Qhat(2x2)^dagger  PART*/
  /* Here the  M_oe Mee^-1 M_eo  implementation  */
  Hopping_Matrix_32_orphaned(EO, g_spinor_field32[0], k_charm);
  Hopping_Matrix_32_orphaned(EO, g_spinor_field32[1], k_strange);

  M_ee_inv_ndpsi_32_orphaned(g_spinor_field32[2], g_spinor_field32[3],
		 g_spinor_field32[0], g_spinor_field32[1],
		 (float) g_mubar, (float) g_epsbar);

  Hopping_Matrix_32_orphaned(OE, g_spinor_field32[0], g_spinor_field32[2]);
  Hopping_Matrix_32_orphaned(OE, g_spinor_field32[1], g_spinor_field32[3]);

  /* Here the M_oo  implementation  */
  M_oo_sub_g5_ndpsi_32_orphaned(g_spinor_field32[2], g_spinor_field32[3], k_charm, k_strange,
  		    g_spinor_field32[0], g_spinor_field32[1],
  		    (float)(-g_mubar), (float)(-g_epsbar));
  /* We have to reassigin as follows to avoid overwriting */
  /* Recall in fact that   Q^hat = tau_1 Q tau_1  , hence  */
  /* and then the  Qhat(2x2)  PART */

  /* Here the  M_oe Mee^-1 M_eo  implementation  */
  Hopping_Matrix_32_orphaned(EO, g_spinor_field32[0], g_spinor_field32[3]);
  Hopping_Matrix_32_orphaned(EO, g_spinor_field32[1], g_spinor_field32[2]);

  M_ee_inv_ndpsi_32_orphaned(g_spinor_field32[5], g_spinor_field32[4],
		 g_spinor_field32[1], g_spinor_field32[0],
		 (float)(-g_mubar), (float)g_epsbar);

  Hopping_Matrix_32_orphaned(OE, l_strange, g_spinor_field32[4]);
  Hopping_Matrix_32_orphaned(OE, l_charm, g_spinor_field32[5]);

  /* Here the M_oo  implementation  */
  M_oo_sub_g5_ndpsi_32_orphaned(l_strange, l_charm, g_spinor_field32[3], g_spinor_field32[2],
		    l_strange, l_charm, (float)(-g_mubar), (float)(-g_epsbar));
  /* At the end, the normalisation by the max. eigenvalue  */ 
  /* Twice  phmc_invmaxev  since we consider here  D Ddag  !!! */
  mul_r_32_orphaned(l_charm, (float) phmc_invmaxev*phmc_invmaxev, l_charm, VOLUME/2);
  mul_r_32_orphaned(l_strange, (float) phmc_invmaxev*phmc_invmaxev, l_strange, VOLUME/2);
#ifdef TM_USE_OMP
  } /* OpenMP closing brace */
#endif
  return;
}

void Qtm_pm_ndpsi_shift_32(spinor32 * const l_strange, spinor32 * const l_charm,
          spinor32 * const k_strange, spinor32 * const k_charm){
  Qtm_pm_ndpsi_32(l_strange, l_charm, k_strange, k_charm);
  assign_add_mul_r_32(l_strange, k_strange, (float)g_shift, VOLUME/2 );
  assign_add_mul_r_32(l_charm, k_charm, (float)g_shift, VOLUME/2 );
  return;
}

void Qsw_pm_ndpsi_32(spinor32 * const l_strange, spinor32 * const l_charm,
      spinor32 * const k_strange, spinor32 * const k_charm) {
#ifdef TM_USE_OMP
#pragma omp parallel
  {
#endif
  /* FIRST THE  Qhat(2x2)^dagger  PART*/
  /* Here the  M_oe Mee^-1 M_eo  implementation  */
  Hopping_Matrix_32_orphaned(EO, g_spinor_field32[0], k_charm);
  Hopping_Matrix_32_orphaned(EO, g_spinor_field32[1], k_strange);

  assign_mul_one_sw_pm_imu_eps_32_orphaned(EE, g_spinor_field32[2], g_spinor_field32[3], 
             g_spinor_field32[0], g_spinor_field32[1], -g_mubar, g_epsbar);
  clover_inv_nd_32_orphaned(EE, g_spinor_field32[2], g_spinor_field32[3]);

  Hopping_Matrix_32_orphaned(OE, g_spinor_field32[0], g_spinor_field32[2]);
  Hopping_Matrix_32_orphaned(OE, g_spinor_field32[1], g_spinor_field32[3]);

  // Here the M_oo  implementation  
  clover_gamma5_nd_32_orphaned(OO, g_spinor_field32[2], g_spinor_field32[3], 
         k_charm, k_strange,
         g_spinor_field32[0], g_spinor_field32[1],
         -g_mubar, -g_epsbar);

  // and then the  Qhat(2x2)  PART 
  // Recall in fact that   Q^hat = tau_1 Q tau_1  
  // Here the  M_oe Mee^-1 M_eo  implementation  
  // the re-ordering in s and c components is due to tau_1
  Hopping_Matrix_32_orphaned(EO, g_spinor_field32[0], g_spinor_field32[3]);
  Hopping_Matrix_32_orphaned(EO, g_spinor_field32[1], g_spinor_field32[2]);

  assign_mul_one_sw_pm_imu_eps_32_orphaned(EE, g_spinor_field32[4], g_spinor_field32[5], 
             g_spinor_field32[1], g_spinor_field32[0], g_mubar, g_epsbar);
  clover_inv_nd_32_orphaned(EE, g_spinor_field32[4], g_spinor_field32[5]);

  Hopping_Matrix_32_orphaned(OE, g_spinor_field32[0], g_spinor_field32[5]);
  Hopping_Matrix_32_orphaned(OE, g_spinor_field32[1], g_spinor_field32[4]);

  clover_gamma5_nd_32_orphaned(OO, l_charm, l_strange,
         g_spinor_field32[2], g_spinor_field32[3],
         g_spinor_field32[1], g_spinor_field32[0],
         g_mubar, -g_epsbar);

  /* At the end, the normalisation by the max. eigenvalue  */ 
  /* Twice  phmc_invmaxev  since we consider here  D Ddag  !!! */
  mul_r_32_orphaned(l_charm, phmc_invmaxev*phmc_invmaxev, l_charm, VOLUME/2);
  mul_r_32_orphaned(l_strange, phmc_invmaxev*phmc_invmaxev, l_strange, VOLUME/2);

#ifdef TM_USE_OMP /* OpenMP parallel closing brace */
  }
#endif

  return;
}

void Qsw_pm_ndpsi_shift_32(spinor32* const l_strange, spinor32 * const l_charm,
      spinor32 * const k_strange, spinor32 * const k_charm){
  Qsw_pm_ndpsi_32(l_strange,l_charm,k_strange,k_charm);
  assign_add_mul_r_32(l_strange, k_strange, (float)g_shift, VOLUME/2 );
  assign_add_mul_r_32(l_charm, k_charm, (float)g_shift, VOLUME/2 );
  return;
}



/***********************************************************************
 *
 * Copyright (C) 2016 Simone Bacchio, Jacob Finkenrath
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
 * Interface for DDalphaAMG
 *
 *******************************************************************************/

#ifndef DDalphaAMG_INTERFACE_H_
#define DDalphaAMG_INTERFACE_H_
#include "global.h"
#include "su3.h"
#include"solver/matrix_mult_typedef.h"
#include"solver/matrix_mult_typedef_nd.h"

extern int mg_setup_iter;
extern int mg_coarse_setup_iter;
extern int mg_update_setup_iter;
extern int mg_omp_num_threads;
extern int mg_Nvec;
extern int mg_lvl;
extern int mg_blk[4];
extern int mg_mixed_prec;
extern int mg_setup_mu_set;
extern int mg_no_shifts;
extern double mg_mms_mass;
extern double mg_setup_mu;
extern double mg_cmu_factor;
extern double mg_dtau_update;
extern double mg_rho_update;

void MG_init(void);
void MG_update_gauge(double step);
void MG_update_mu(double mu_tmLQCD, double odd_tmLQCD);
void MG_update_mubar_epsbar(double mubar_tmLQCD, double epsbar_tmLQCD, double shift_tmLQCD);
void MG_reset(void);
void MG_finalize(void);

int MG_solver(spinor * const phi_new, spinor * const phi_old,
	      const double precision, const int max_iter,const int rel_prec,
	      const int N, su3 **gf, matrix_mult f);

int MG_solver_eo(spinor * const Even_new, spinor * const Odd_new,
		 spinor * const Even, spinor * const Odd,
		 const double precision, const int max_iter, const int rel_prec,
		 const int N, su3 **gf, matrix_mult_full f_full);

int MG_solver_nd(spinor * const up_new, spinor * const dn_new,
		 spinor * const up_old, spinor * const dn_old,
		 const double precision, const int max_iter, const int rel_prec,
		 const int N, su3 **gf, matrix_mult_nd f);

int MG_solver_nd_eo(spinor * const Even_new_up, spinor * const Odd_new_up, 
                    spinor * const Even_new_dn, spinor * const Odd_new_dn,
                    spinor * const Even_up, spinor * const Odd_up,
                    spinor * const Even_dn, spinor * const Odd_dn,
                    const double precision, const int max_iter, const int rel_prec,
                    const int N, su3 **gf, matrix_mult_full_nd f_full);

int MG_mms_solver_nd(spinor **const up_new, spinor **const dn_new,
                     spinor * const up_old, spinor * const dn_old,
                     const double * shifts, const int no_shifts,
                     const double * precision, const int max_iter, const int rel_prec,
                     const int N, su3 **gf, matrix_mult_nd f);

#endif /* DDalphaAMG_INTERFACE_H_ */

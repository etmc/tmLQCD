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
 * Interface for MG4QCD
 *
 *******************************************************************************/

#ifndef MG4QCD_INTERFACE_H_
#define MG4QCD_INTERFACE_H_
#include "global.h"
#include "su3.h"
#include"solver/matrix_mult_typedef.h"

extern int mg_setup_iter;
extern int mg_Nvec;
extern int mg_lvl;
extern int mg_blk[4];
extern double mg_cmu_factor;

void MG_init(void);
void MG_update_mu(double mu);
void MG_finalize(void);
int MG_solver(spinor * const Even_new, spinor * const Odd_new,
                   spinor * const Even, spinor * const Odd,
                   const double precision, const int max_iter,
                   const int solver_flag, const int rel_prec,
                   const int even_odd_flag, su3 **gf, matrix_mult f);
#endif /* MG4QCD_INTERFACE_H_ */
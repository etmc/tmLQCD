/***********************************************************************
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008 Carsten Urbach
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

#ifndef _PHMC_H
#define _PHMC_H

#include "solver/matrix_mult_typedef_bi.h"

/* the normalisation constant appearing in the product representation of */
/* the polynomial */
extern double phmc_Cpol;
/* maximal and minimal eigenvalue of the ND operator */
extern double phmc_cheb_evmin, phmc_cheb_evmax;
/* inverse maximal EV, needed for normalisation */
extern double phmc_invmaxev;
/* These are the roots */
extern _Complex double * phmc_root;
/* degree and coefs of P */
extern int phmc_dop_n_cheby;
extern double * phmc_dop_cheby_coef;
/* degree of coefs \tilde P */
extern int phmc_ptilde_n_cheby;
extern double * phmc_ptilde_cheby_coef;
extern int phmc_max_ptilde_degree;

/* structure for holding a set of phmc specific variables*/
typedef struct phmc_vars_ {
  void *previous,*next;
  double invmaxev;
  double Cpol;
  int dop_n_cheby;
  _Complex double *root;
  int stacksize;
} phmc_vars;

/* stack for saving and restoring phmc variables*/
extern phmc_vars *phmc_var_stack;

/* functions for pushing and poping phmc vars */
void pushPhmcVars();
void popPhmcVars();

void phmc_compute_ev(const int trajectory_counter, const int id,
		     matrix_mult_bi Qsq);

#endif

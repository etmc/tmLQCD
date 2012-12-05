/***********************************************************************
 *
 * Copyright (C) 2003 Ines Wetzorke
 *               2009 Carsten Urbach
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
 * Action of the overlap Dirac operator D on a given spinor field
 *
 * This software is distributed under the terms of the GNU General Public
 * License (GPL)
 *
 ************************************************************************/

#ifndef _DOV_PSI_H
#define _DOV_PSI_H

#include "su3.h"


/**
 * this is for bookeeping of auxiliary spinors among the routines
 *                                          index 
 * Dov_psi:           (1 auxiliary spinors)   0
 *
 * addproj_q_invsqrt: (1 auxiliary spinor )   1  > these functions share their contigent
 * norm_Q_sqr_psi   : (1 auxiliary spinor )   1  > as they are not called at the same time
 * norm_Q_n_psi     : (1 auxiliary spinor )   1  > (if this is not the case anymore it will be detected by the code see below)
 *
 * Q_over_sqrt_Q_sqr: (5 auxiliary spinors)  2-6
 *                    --------------------------
 *                     7 auxiliary spinors   0-6
 *
 * for additional safety the Dov_WS struct has a lock_map member tacking track of 
 * locks to specific spinors, if a function requests one spinor before it has been 
 * unlocked by a previously called function lock_Dov_WS_spinor will strike
 */
typedef struct Dov_WS_{
  int n_spinors;
  spinor **dum_spinors;
  spinor *dum_spinors_membuf;
  int *lock_map;
} Dov_WS;

extern double m_ov;
extern int ov_n_cheby;
extern double * ov_cheby_coef;
extern Dov_WS *dov_ws;

void Dov_psi(spinor * const, spinor * const);
void Dov_psi_prec(spinor * const, spinor * const);
void Qov_psi(spinor * const, spinor * const);
void Qov_sq_psi(spinor * const P, spinor * const S);
void Qov_sq_psi_prec(spinor * const P, spinor * const S);

void Q_over_sqrt_Q_sqr(spinor * const R, double * const c, 
		       const int n, spinor * const S,
		       const double rnorm, const double minev);

void calculateOverlapPolynomial();

void free_Dov_WS();
#endif

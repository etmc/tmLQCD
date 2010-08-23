/***********************************************************************
 * $Id$
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

extern double m_ov;
extern int ov_n_cheby;
extern double * ov_cheby_coef;

void Dov_psi(spinor * const, spinor * const);
void Qov_psi(spinor * const, spinor * const);

void Q_over_sqrt_Q_sqr(spinor * const R, double * const c, 
		       const int n, spinor * const S,
		       const double rnorm, const double minev);


#endif

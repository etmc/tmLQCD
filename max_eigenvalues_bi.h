/***********************************************************************
 *
 * Copyright (C) 2006 Thomas Chiarappa
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

#ifndef _MAX_EIGENVALUES_BI_H
#define _MAX_EIGENVALUES_BI_H

extern bispinor * max_evs;
extern double * max_evls;
/*
extern int eigenvalues_for_cg_computed;
*/
double max_eigenvalues_bi(int * nev, const int operator_flag, const int max_iterations, const double prec);

#endif

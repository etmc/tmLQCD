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
#ifndef _EIGENVALUES_H
#define _EIGENVALUES_H

extern spinor * eigenvectors;
extern double * eigenvls;
extern double * inv_eigenvls;
extern int eigenvalues_for_cg_computed;
extern int no_eigenvalues;
extern int evlength;

/* the folowing two are needed for the overlap */
extern double ev_qnorm, ev_minev;

double eigenvalues(int * nr_of_eigenvalues, const int max_iterations, 
		   const double precision, const int maxmin,
		   const int readwrite, const int nstore,
		   const int even_odd_flag);

#endif

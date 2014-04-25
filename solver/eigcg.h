/****************************************************************************
 * Copyright (C) 2008,2009,2010,2011,2012  
 * Andreas Stathopoulos, Kostas Orginos, Abdou M. Abdel-Rehim
 *
 * This program is based on interfacing the eigCG solver to the tmLQCD code.
 * It was written by Abdou M. Abdel-Rehim. The original code was written
 * by Andreas Stathopoulos and Kostas Orginos and integrated in Chroma.
 * In this interface we use functions from tmLQCD.
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
 ******************************************************************************/

#ifndef _EIGCG_H
#define _EIGCG_H

#include "su3.h"
#include "solver/matrix_mult_typedef.h"


void eigcg(int n, int lde, spinor * const x, spinor * const b, double *normb, const double eps_sq, 
           double restart_eps_sq, const int rel_prec, int maxit, int *iter, double *reshist, int *flag,  
           spinor **work, matrix_mult f, int nev, int v_max, spinor *V, int esize, _Complex double *ework);



#endif

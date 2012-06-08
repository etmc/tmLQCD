/***********************************************************************
 *
 *
 * Copyright (C) 2004 Andrea Shindler
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
 ***********************************************************************/

#ifndef _CG_MMS_TM_H
#define _CG_MMS_TM_H

#include"matrix_mult_typedef.h"
#include"su3.h"

int cg_mms_tm(spinor * const P,spinor * const Q, const int max_iter, 
	      double eps_sq, const int rel_prec, const int N, matrix_mult f,
        const int no_extra_masses, double * const extra_masses, 
        const int id );

#endif

/***********************************************************************
 *
 *
 * Copyright (C) 2016 Simone Bacchio
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

#ifndef _INIT_GUESS_H
#define _INIT_GUESS_H

#include"su3.h"
#include"solver.h"

int init_guess_mms(spinor ** const P, spinor * const Q,
                   int shift, solver_params_t * const params);

int init_guess_mms_nd(spinor ** const Pup, spinor ** const Pdn, 
                      spinor * const Qup, spinor * const Qdn, 
                      int shift, solver_params_t * solver_params);

int init_guess_mms_nd_plus(spinor ** const Pup, spinor ** const Pdn, 
                           spinor * const Qup, spinor * const Qdn, 
                           int shift, solver_params_t * solver_params);
#endif

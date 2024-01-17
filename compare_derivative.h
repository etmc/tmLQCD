/***********************************************************************
 *
 * Copyright (C) 2024 Bartosz Kostrzewa
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

#ifndef COMPARE_DERIVATIVE_H
#define COMPARE_DERIVATIVE_H

#include "monomial/monomial.h"
#include "su3adj.h"

void compare_derivative(monomial *mnl, su3adj **ext_lib, su3adj **native, const double threshold, const char * name);

#endif

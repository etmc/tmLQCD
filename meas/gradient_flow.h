/***********************************************************************
 *
 * Copyright (C) 2013 Albert Deuzeman
 *               2015 Bartosz Kostrzewa
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

#ifndef _GRADIENT_FLOW_H
#define _GRADIENT_FLOW_H

#include "su3.h"

void step_gradient_flow(su3 ** vt, su3 ** x1, su3 ** x2, su3 ** z, const unsigned int type, const double eps);
void gradient_flow_measurement(const int traj, const int id, const int ieo);

#endif

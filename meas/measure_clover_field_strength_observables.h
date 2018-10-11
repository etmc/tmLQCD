/***********************************************************************
 *
 * Copyright (C) 2015,2018 Bartosz Kostrzewa
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

#ifndef MEASURE_CLOVER_FIELD_STRENGTH_OBSERVABLES_H
#define MEASURE_CLOVER_FIELD_STRENGTH_OBSERVABLES_H

#include "meas/field_strength_types.h"

void measure_clover_field_strength_observables(const su3 ** const gf, field_strength_obs_t * const ret);

#endif  

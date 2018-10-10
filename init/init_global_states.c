/***********************************************************************
 *
 * Copyright (C) 2018  Bartosz Kostrzewa
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
 *******************************************************************************/

#include "global.h"
#include "misc_types.h"

void init_global_states(void) {
  g_gauge_state = new_tm_GaugeState("g_gauge_field");
  g_gauge_state_32 = new_tm_GaugeState("g_gauge_field_32");
  g_clover_state = new_tm_CloverState();
  g_clover_state_32 = new_tm_CloverState();
  g_clover_inverse_state = new_tm_CloverInverseState();
  g_clover_inverse_state_32 = new_tm_CloverInverseState();
}

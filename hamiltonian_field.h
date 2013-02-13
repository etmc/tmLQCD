/***********************************************************************
 *
 * Copyright (C) 2012 Carsten Urbach
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

#ifndef _HAMILTONIAN_FIELD_H
#define _HAMILTONIAN_FIELD_H

#include <su3.h>
#include <su3adj.h>

typedef struct {
  su3 ** gaugefield;
  su3adj ** momenta;
  su3adj ** derivative;
  int update_gauge_copy;
  int traj_counter;
} hamiltonian_field_t;



#endif

/***********************************************************************
 * Copyright (C) 2010 Florian Burger
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

#ifndef _TEMPORALGAUGE_H
#define _TEMPORALGAUGE_H

int init_temporalgauge_trafo(const int V, su3** gfield);
void apply_gtrafo(su3 ** gfield, su3 * trafofield);
void apply_gtrafo_spinor(spinor * spin, su3 * trafofield);
void apply_inv_gtrafo(su3 ** gfield, su3 * trafofield);
void apply_inv_gtrafo_spinor(spinor * spin, su3 * trafofield);
void finalize_temporalgauge();

void apply_gtrafo_spinor_odd(spinor * spin, su3 * trafofield);
void apply_inv_gtrafo_spinor_odd(spinor * spin, su3 * trafofield);
void apply_gtrafo_spinor_even(spinor * spin, su3 * trafofield);
void apply_inv_gtrafo_spinor_even(spinor * spin, su3 * trafofield);

void copy_gauge_field(su3** to, su3** from);

#endif





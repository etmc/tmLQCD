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
#ifndef _INIT_SPINOR_FIELD_H
#define _INIT_SPINOR_FIELD_H

int init_spinor_field(const int V, const int nr);
int init_spinor_field_32(const int V, const int nr);
int init_csg_field(const int V);

int allocate_spinor_field_array(spinor ***spinors,spinor **sp,const int V, const int nr);
void free_spinor_field_array(spinor** sp);

void free_spinor_field();
void free_spinor_field_32();
int init_timslice_buffer_field(const int t_slice);

#endif

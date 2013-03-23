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
#ifndef _GENERATE_DFL_SUBSPACE
#define _GENERATE_DFL_SUBSPACE

#include "su3.h"
#include <complex.h>

int init_dfl_subspace(const int);
int free_dfl_subspace();
int generate_dfl_subspace(const int Ns, const int N, const int repro);
int generate_dfl_subspace_free(const int Ns, const int N);

extern spinor ** dfl_fields;
extern _Complex double ** little_dfl_fields;
extern _Complex double ** little_dfl_fields_eo;

#endif

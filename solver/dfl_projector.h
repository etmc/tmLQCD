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
#ifndef _DFL_PROJECTOR_H
#define _DFL_PROJECTOR_H

#include "su3spinor.h"

void project(spinor * const out, spinor * const in);
void project_left(spinor * const out, spinor * const in);
void project_right(spinor * const out, spinor * const in);
void project_left_D(spinor * const out, spinor * const in);
void D_project_right(spinor * const out, spinor * const in);
int check_projectors();
void check_little_D_inversion();
void check_local_D();
void free_dfl_projector();

void little_project(complex * const out, complex * const in, const int  N);
void little_P_L_D(complex * const out, complex * const in);
void little_P_L_D_sym(complex * const out, complex * const in);
void little_D_P_R(complex * const out, complex * const in);
void little_P_R(complex * const out, complex * const in);
void little_P_L(complex * const out, complex * const in);
void little_P_R_sym(complex * const out, complex * const in);
void little_P_L_sym(complex * const out, complex * const in);

extern double dfl_little_D_prec;
extern int dfl_sloppy_prec;


#endif

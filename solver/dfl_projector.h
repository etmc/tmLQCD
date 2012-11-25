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
int check_projectors(const int repro);
void check_little_D_inversion(const int repro);
void check_local_D(const int repro);
void free_dfl_projector();

void little_project(_Complex double * const out, _Complex double * const in, const int  N);
void little_project_eo(_Complex double * const out, _Complex double * const in, const int N);
void little_P_L_D(_Complex double * const out, _Complex double * const in);
void little_P_L_D_sym(_Complex double * const out, _Complex double * const in);
void little_D_P_R(_Complex double * const out, _Complex double * const in);
void little_P_R(_Complex double * const out, _Complex double * const in);
void little_P_L(_Complex double * const out, _Complex double * const in);
void little_P_R_sym(_Complex double * const out, _Complex double * const in);
void little_P_L_sym(_Complex double * const out, _Complex double * const in);

extern double dfl_little_D_prec;
extern int dfl_sloppy_prec;


#endif

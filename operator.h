/***********************************************************************
 * $Id$ 
 *
 * Copyright (C) 2009 Carsten Urbach
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

#ifndef _OPERATOR_H
#define _OPERATOR_H

#include "spinor.h"

#define TMWILSON 0
#define OVERLAP 1
#define WILSON 2

#define max_no_operators 1

typedef struct {
  int type;
  int n_cheby;
  int deg_poly;
  int sloppy_precision;
  int even_odd;
  
  double kappa;
  double m;
  double s;

  double * coefs;
  void (*applyD) (spinor * const, spinor * const);
  void (*applyQ) (spinor * const, spinor * const);
  void (*applyQp) (spinor * const, spinor * const);
  void (*applyQm) (spinor * const, spinor * const);
  void (*applyQsq) (spinor * const, spinor * const);
} operator;

extern operator operator_list[max_no_operators];

#endif

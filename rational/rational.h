/***********************************************************************
 *
 * Copyright (C) 2013 Carsten Urbach
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

#ifndef _RATIONAL_H
#define _RATIONAL_H

typedef struct {
  int order, np;
  int crange[2];
  double range[2];
  double eps;
  double A, delta;
  double *mu,*rmu;
  double *nu,*rnu;
} rational_t;

int init_rational(rational_t * rat, const unsigned int scale);
int free_rational(rational_t * rat);

#endif

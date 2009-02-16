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
#ifndef _INTEGRATOR_H
#define _INTEGRATOR_H

#define LEAPFROG 1
#define SEXTON 2
#define EXTLEAPFROG 3
#define EXTSEXTON 4
#define IMPRLEAPFROG 5
#define MN2 6
#define MN2p 7

typedef void (*integratefk)(const double, const int, const int);

typedef struct{
  int type[10];
  int no_timescales;
  int n_int[10];
  double tau;
  double lambda[10];
  int mnls_per_ts[10][10];
  int no_mnls_per_ts[10];
  integratefk integrate[10];
} integrator;

extern integrator Integrator;


int init_integrator();
void integrate_2mn(const double tau, const int S, const int halfstep);
void integrate_2mnp(const double tau, const int S, const int halfstep);
void integrate_leap_frog(const double tau, const int S, const int halfstep);


#endif

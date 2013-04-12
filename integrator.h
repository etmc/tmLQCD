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

#include <su3.h>
#include <su3adj.h>
#include <hamiltonian_field.h>

#define LEAPFROG 1
#define SEXTON 2
#define EXTLEAPFROG 3
#define EXTSEXTON 4
#define IMPRLEAPFROG 5
#define MN2 6
#define MN2p 7
#define OMF4 8

typedef void (*integratefk)(const double, const int, const int);

typedef struct {
  /* gauge, momenta and derivative fields to be used during integration */
  hamiltonian_field_t hf;
  /* list of types of integrators */
  int type[10];
  /* number of timescales */
  int no_timescales;
  /* monitor forces */
  int monitor_forces;
  /* steps per timescale */
  int n_int[10];
  /* trajectory length */
  double tau;
  /* lambda parameter for 2MN integration scheme */
  double lambda[10];
  /* monomials per timescale */
  int mnls_per_ts[10][10];
  /* number of monomials per timescale */
  int no_mnls_per_ts[10];
  /* function pointers to integration scheme functions */
  integratefk integrate[10];
} integrator;

extern integrator Integrator;

/* all following functions are currently defined in integrator.c */
/* function to initialise the integrator, to be called once at the beginning */
int init_integrator();
/* function to set the gauge and momenta fields for the integration */
void integrator_set_fields(hamiltonian_field_t * hf);
/* and unsets again (to NULL pointer ) */
void integrator_unset_fields();

#endif

/***********************************************************************
 *
 * Measurements of the reweighting factors by Georg Bergner 2016
 *
 * Copyright (C) 20016 Georg Bergner
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

#ifndef _REWEIGHTING_MEASUREMENT_H
#define _REWEIGHTING_MEASUREMENT_H
#include <stdlib.h>
#include <stdio.h>

/**
 * This struct collects data for the Chebyshev reweighting factor measurements.
 * These are the coefficients for the log(x) approximation provided from outside.
 */
typedef struct{
	double* el;
	unsigned int s;
} vector_list;

inline void allocate(vector_list* v, const unsigned int length){
  v->s=length;
  v->el = malloc(length * sizeof(double));
}

/**
 * This struct collects data for the Chebyshev reweighting factor measurements.
 * The split_list allows to combine different approximations.
 */
typedef struct {
	unsigned int * ord;
	unsigned int* est;
	unsigned int s;

} split_list;


/**
 * Parameters for the reweighting factor measurements.
 */
typedef struct {
	  int reweighting_operator;
	  int reweighting_number_sources;
	  int use_evenodd;
	  double k2mu0;
	  double kappa0;
	  vector_list kappaarray;
	  double rmu0;
	  double rmu;
	  double minev;
	  double maxev;
	  double testchebconvergence;
	  int interpolationsteps;
	  int estimatorscheb;
	  int cheborder;
	  int use_cheb;
	  int only_cheb;
	  int evest;
	vector_list coeff;
	split_list splitlist;
} reweighting_parameter;

/**
 * Destructor for parameter.
 * @param par
 */
void free_reweighting_parameter(void* par);

/**
 * Constructor for parameter.
 * @param parameter
 */
void initialize_reweighting_parameter(void** parameter);

/**
 * main measurement.
 * @param traj
 * @param t0
 * @param ieo
 */
void reweighting_measurement(const int traj, const int t0, const int ieo);

#endif

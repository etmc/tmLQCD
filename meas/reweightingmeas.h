/***********************************************************************
 *
 * Measurements of the reweighting factors by Georg Bergner 2016
 *
 * Copyright (C) 2008 Carsten Urbach
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
typedef struct{
	  int reweighting_operator;
	  int reweighting_number_sources;
	  int use_evenodd;
	  double k2mu0;
	  double kappa0;
	  int interpolationsteps;
	  double rmu0;
	  double rmu;

} reweighting_parameter;

inline void initialize_reweighting_parameter(void** parameter){
	reweighting_parameter* param;
	if(!(*parameter)){
		(*parameter)=malloc(sizeof(reweighting_parameter));
		param=(reweighting_parameter*)(*parameter);
		param->reweighting_number_sources=0;
		param->reweighting_operator=0;
		param->use_evenodd=0;
		param->k2mu0=0.0;
		param->rmu0=0.0;
		param->rmu=0.0;
		param->kappa0=0.0;
		param->interpolationsteps=1;
	}
}

void reweighting_measurement(const int traj, const int t0, const int ieo);

#endif

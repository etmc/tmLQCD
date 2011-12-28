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

#include "params.ih"
#include "solver/solver.h"

/* This needs fixing */

paramsInverterInfo *construct_paramsInverterInfo(double const epssq, const int iter, 
						 const int solver, const int noflavours) {
  int i;
  struct timeval t1;
  paramsInverterInfo *info = malloc(sizeof(paramsInverterInfo));

  if (info == (paramsInverterInfo*)NULL)
    kill_with_error(NULL, g_cart_id, "Could not allocate paramsInverterInfo.");

  gettimeofday(&t1, NULL);

  info->iter = iter;
  info->epssq = epssq;
  info->noflavours = noflavours;

  info->kappa = g_kappa;
  info->mu = g_mu / 2. / g_kappa;

  strcpy(info->package_version, PACKAGE_VERSION);

  if(noflavours == 2) {
    info->mubar = g_mubar / 2. / g_kappa;
    info->epsbar = g_epsbar / 2. / g_kappa;
  }
  else {
    info->mubar = 0.;
    info->epsbar = 0.;
  }
  strcpy(info->date, ctime(&t1.tv_sec));
  info->mms = 0;
  info->heavy = 0;
  info->cgmms_mass = 0;
  switch (solver) {
  case CG:
    strcpy(info->inverter, "CG");
    break;
  case BICGSTAB:
    strcpy(info->inverter, "BiCGstab");
    break;
  case GMRES:
    strcpy(info->inverter, "GMRES");
    break;
  case CGMMS:
    strcpy(info->inverter, "CGMMS");
    info->mms = 1;
    break;
  case CGS:
    strcpy(info->inverter, "CGS");
    break;
  default:
    strcpy(info->inverter, "other");
    break;
  }
  return(info);
}

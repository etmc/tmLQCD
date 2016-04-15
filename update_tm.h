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
#ifndef _UPDATE_TM_H
#define _UPDATE_TM_H

#ifdef MG4QCD
typedef struct{
   double tauMC;
   int gcopy_up2date;
   int basis_up2date;
   double tau_basis;
}hmc_control_t;

hmc_control_t reset_hmc_control(void);
hmc_control_t update_hmc_control(double dtau);
hmc_control_t get_hmc_control(void);
hmc_control_t set_hmc_control(int gcopy_up2date,int basis_up2date,double tau_basis);
#endif

int update_tm(double *plaquette_energy, double *rectangle_energy, 
	      char * filename, const int return_check, const int acctest, 
	      const int traj_counter);

#endif

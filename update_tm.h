#ifndef _UPDATE_TM_H
#define _UPDATE_TM_H

int update_tm(const int integtyp, double * plaquette_energy, double * rectangle_energy, char * filename, 
	      const double dtau, const int Nsteps, const int nsmall, 
	      const double tau, int * n_int, const int return_check,
	      const double q_off, const double q_off2);

#endif

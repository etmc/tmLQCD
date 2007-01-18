#ifndef _UPDATE_TM_ND_H
#define _UPDATE_TM_ND_H

int update_tm_nd(const int integtyp, double * plaquette_energy, double * rectangle_energy, char * filename, 
	      const double dtau, const int Nsteps, const int nsmall, 
	      const double tau, int * n_int, const int return_check,
	      double * lambda, const int rngrepro);

#endif

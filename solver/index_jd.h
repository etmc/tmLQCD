#ifndef _INDEX_JD_H
#define _INDEX_JD_H

#ifndef RESTART_JACOBI_DAVIDSON
#define RESTART_JACOBI_DAVIDSON 5
#endif
#ifndef RESTART_RITZ_JACOBI
#define RESTART_RITZ_JACOBI 6
#endif

void index_jd(int * nr_of_eigenvalues, 
	      const int max_iterations, const double precision, 
	      char * conf_filename, const int nstore,
	      const int method);

#endif

/*******************************************************************************
 * $Id$
 *
 * Here we compute the nr_of_eigenvalues lowest eigenvalues
 * of (gamma5*D)^2. Therefore we use the arnoldi routines.
 * 
 * The computed eigenvalues are stored in g_eigenvalues
 * and the computed eigenvectors in g_ev
 * 
 * inout:
 *   nr_of_eigenvalues:      input:  Number of eigenvalues to compute
 *                           output: Number of computed eigenvalues
 *
 * Autor: Carsten Urbach <urbach@ifh.de>
 *
 *******************************************************************************/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "global.h"
#include "su3.h"
#include "linalg_eo.h"
#include "start.h"
#include "tm_operators.h"
#include "solver/solver.h"
#include "solver/jdher.h"
#include "eigenvalues.h"


spinor  *eigenvectors = NULL;
double * eigenvls = NULL;
int eigenvalues_for_cg_computed = 0;

double eigenvalues(int * nr_of_eigenvalues, 
		 const int max_iterations, const double precision,const int maxmin) {
#ifdef HAVE_LAPACK
  static spinor * eigenvectors_ = NULL;
  static int allocated = 0;

  /**********************
   * For Jacobi-Davidson 
   **********************/
  int verbosity = g_debug_level, converged = 0, blocksize = 1, blockwise = 0;
  int solver_it_max = 50, j_max, j_min;
  /*int it_max = 10000;*/
  /* complex *eigv_ = NULL, *eigv; */
  double decay_min = 1.7, decay_max = 1.5, prec,
    threshold_min = 1.e-3, threshold_max = 5.e-2,
    startvalue, threshold, decay, returnvalue;
  /* static int v0dim = 0; */
  int v0dim = 0;
  

  /**********************
   * General variables
   **********************/
  int returncode=0;

  if(maxmin == JD_MINIMAL) {
    startvalue = 0.;
    threshold = threshold_min;
    decay = decay_min;
    solver_it_max = 200;
  }
  else {
    startvalue = 50.;
    threshold = threshold_max;
    decay = decay_max;
    solver_it_max = 50;
  }

  if(g_proc_id == g_stdio_proc && g_debug_level > 0) {
    printf("\nNumber of lowest eigenvalues to compute = %d\n\n",(*nr_of_eigenvalues));
    printf("Using Jacobi-Davidson method! \n");
  }
  eigenvalues_for_cg_computed = 1;

  if((*nr_of_eigenvalues) < 8){
    j_max = 15;
    j_min = 8;
  }
  else{
    j_max = 2*(*nr_of_eigenvalues);
    j_min = (*nr_of_eigenvalues);
  }
  if(precision < 1.e-14){
    prec = 1.e-14;
  }
  else{
    prec = precision;
  }

  if(allocated == 0) {
    allocated = 1;
#if (defined SSE || defined SSE2 || defined SSE3)
    eigenvectors_ = calloc((VOLUMEPLUSRAND)/2*(*nr_of_eigenvalues)+1, sizeof(spinor)); 
    eigenvectors = (spinor *)(((unsigned long int)(eigenvectors_)+ALIGN_BASE)&~ALIGN_BASE);
#else
    eigenvectors_= calloc((VOLUMEPLUSRAND)/2*(*nr_of_eigenvalues), sizeof(spinor));
    eigenvectors = eigenvectors_;
#endif
    eigenvls = (double*)malloc((*nr_of_eigenvalues)*sizeof(double));
  }

  /* compute minimal eigenvalues */

  jdher((VOLUME)/2*sizeof(spinor)/sizeof(complex), (VOLUMEPLUSRAND)/2*sizeof(spinor)/sizeof(complex),
	startvalue, prec, 
	(*nr_of_eigenvalues), j_max, j_min, 
	max_iterations, blocksize, blockwise, v0dim, (complex*) eigenvectors,
	BICGSTAB, solver_it_max,
	threshold, decay, verbosity,
	&converged, (complex*) eigenvectors, eigenvls,
	&returncode, maxmin, 1,
	&Qtm_pm_psi);

  (*nr_of_eigenvalues) = converged;
  /* v0dim = converged; */
  returnvalue=eigenvls[0];
#else
  fprintf(stderr, "lapack not available, so JD method not available \n");
#endif
  return(returnvalue);
}

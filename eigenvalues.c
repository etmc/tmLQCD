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
 * input:
 *   crylov_space_dimension: Dimension of crylov space dimension
 *                           to be used in the arnoldi routines
 *   operator_flag:          Choose if we want to use D_Wilson
 *                           or D_Overlap
 *
 * Autor: Carsten Urbach <urbach@ifh.de>
 *
 *******************************************************************************/

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
#ifdef MPI
#include "solver/pjdher.h"
#endif
#include "eigenvalues.h"

/*********************************************************
 *
 * We need here another function Qsqr_psi, representing
 * (gamma5*D)^2, because this is used in the CG solver
 *
 * It is not identical to Q_sqr_psi and not externally
 * accessible.
 *
 *********************************************************/

spinor  *eigenvectors = NULL;
  double * eigenvls = NULL;
int eigenvalues_for_cg_computed = 0;

void eigenvalues(int * nr_of_eigenvalues, const int operator_flag, 
		 const int max_iterations, const double precision) {

  static spinor * eigenvectors_ = NULL;
  static int allocated = 0;
  spinor  *temp_field, *temp_field_ = NULL, *aux_ = NULL;

  /**********************
   * For Jacobi-Davidson 
   **********************/
  int verbosity = 5, converged = 0, blocksize = 1, blockwise = 0;
  int solver_it_max = 50, j_max, j_min;
  /*int it_max = 10000;*/
  complex *eigv_ = NULL, *eigv;
  double decay_min = 1.7, decay_max = 1.5, prec,
    threshold_min = 1.e-3, threshold_max = 5.e-2;
  static int v0dim = 0;

  /**********************
   * General variables
   **********************/
  int returncode=0;

  if(g_proc_id == g_stdio_proc) printf("\nNumber of lowest eigenvalues to compute = %d\n\n",(*nr_of_eigenvalues));
  eigenvalues_for_cg_computed = 1;

  if(g_proc_id == g_stdio_proc) printf("Using Jacobi-Davidson method! \n");
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
  g_mu = 0.00343;
/*   prec = 1.e-10; */
  if(allocated == 0) {
    allocated = 1;
#if (defined SSE || defined SSE2 )
    eigenvectors_ = calloc((VOLUMEPLUSRAND)/2*(*nr_of_eigenvalues)+1, sizeof(spinor)); 
    eigenvectors = (spinor *)(((unsigned long int)(eigenvectors_)+ALIGN_BASE)&~ALIGN_BASE);
    temp_field_ = calloc((VOLUMEPLUSRAND)/2+1, sizeof(spinor));
    temp_field = (spinor *)(((unsigned long int)(temp_field_)+ALIGN_BASE)&~ALIGN_BASE);
#else
    eigenvectors_= calloc((VOLUMEPLUSRAND)/2*(*nr_of_eigenvalues), sizeof(spinor));
    temp_field_ = calloc((VOLUMEPLUSRAND)/2, sizeof(spinor));
    eigenvectors = eigenvectors_;
    temp_field = temp_field_;
#endif
    eigenvls = (double*)malloc((*nr_of_eigenvalues)*sizeof(double));
  }

  /* compute minimal eigenvalues */

#ifdef MPI
  pjdher((VOLUME)/2*sizeof(spinor)/sizeof(complex), (VOLUMEPLUSRAND)/2*sizeof(spinor)/sizeof(complex),
	 0., prec, 
	 (*nr_of_eigenvalues), j_max, j_min, 
	 max_iterations, blocksize, blockwise, v0dim, (complex*) eigenvectors,
	 CG, solver_it_max,
	 threshold_min, decay_min, verbosity,
	 &converged, (complex*) eigenvectors, eigenvls,
	 &returncode, JD_MINIMAL, 1,
	 &Qtm_pm_psi);
#else
  jdher((VOLUME)/2*sizeof(spinor)/sizeof(complex),
	0., prec, 
	(*nr_of_eigenvalues), j_max, j_min, 
	max_iterations, blocksize, blockwise, v0dim, (complex*) eigenvectors,
	BICGSTAB, solver_it_max,
	threshold_min, decay_min, verbosity,
	&converged, (complex*) eigenvectors, eigenvls,
	&returncode, JD_MINIMAL, 1,
	&Qtm_pm_psi);
#endif
  (*nr_of_eigenvalues) = converged;
  v0dim = converged;

  free(eigenvls);
}

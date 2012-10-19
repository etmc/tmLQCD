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
 *
 * Autor: Thomas Chiarappa
 *        Thomas.Chiarappa@mib.infn.it
 *
 *******************************************************************************/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#ifdef MPI
# include <mpi.h>
#endif
#include "global.h"
#include "su3.h"
#include "linalg_eo.h"
#include "start.h"
#include "operator/tm_operators.h"
#include "solver/solver.h"
#include "solver/jdher_bi.h"
#include "solver/matrix_mult_typedef_bi.h"
#include "eigenvalues_bi.h"
#include "operator/tm_operators_nd.h"


double eigenvalues_bi(int * nr_of_eigenvalues,  
		      const int max_iterations, const double precision,
		      const int maxmin, matrix_mult_bi Qsq) {


  static bispinor * eigenvectors_bi_ = NULL;
  static int allocated = 0;
  static bispinor  *eigenvectors_bi = NULL;
  static double * eigenvls_bi = NULL;

  /**********************
   * For Jacobi-Davidson 
   **********************/
  int verbosity = g_debug_level, converged = 0, blocksize = 1, blockwise = 0;
  int solver_it_max = 200, j_max, j_min; 
  double decay_min = 1.7, decay_max = 1.5, prec,
    threshold_min = 1.e-3, threshold_max = 5.e-2, 
    startvalue, threshold, decay, returnvalue;
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

  if(g_proc_id == g_stdio_proc) {
    printf("Number of %s eigenvalues to compute = %d\n",
	   maxmin ? "maximal" : "minimal",(*nr_of_eigenvalues));
    printf("Using Jacobi-Davidson method! \n");
  }

  if((*nr_of_eigenvalues) < 8){
    j_max = 15;
    j_min = 8;
  }
  else{
    j_max = 2*(*nr_of_eigenvalues);
    j_min = *nr_of_eigenvalues;
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
    eigenvectors_bi_ = calloc((VOLUME)/2*(*nr_of_eigenvalues)+1, sizeof(bispinor)); 
    eigenvectors_bi = (bispinor *)(((unsigned long int)(eigenvectors_bi_)+ALIGN_BASE)&~ALIGN_BASE);
#else
    eigenvectors_bi_= calloc((VOLUME)/2*(*nr_of_eigenvalues), sizeof(bispinor));
    eigenvectors_bi = eigenvectors_bi_;
#endif
    eigenvls_bi = (double*)malloc((*nr_of_eigenvalues)*sizeof(double));
  }

  /* compute eigenvalues */

  if((g_proc_id==0) && (g_debug_level > 4)) {
    printf(" Values of   mu = %e     mubar = %e     eps = %e     precision = %e  \n \n", g_mu, g_mubar, g_epsbar, precision);
  }

  /* here n and lda are equal, because Q_Qdagger_ND_BI does an internal */
  /* conversion to non _bi fields which are subject to xchange_fields   */
  /* so _bi fields do not need boundary                                 */
  jdher_bi((VOLUME)/2*sizeof(bispinor)/sizeof(_Complex double), (VOLUME)/2*sizeof(bispinor)/sizeof(_Complex double),
	   startvalue, prec, 
	   (*nr_of_eigenvalues), j_max, j_min, 
	   max_iterations, blocksize, blockwise, v0dim, (_Complex double*) eigenvectors_bi,
	   BICGSTAB, solver_it_max,
	   threshold, decay, verbosity,
	   &converged, (_Complex double*) eigenvectors_bi, eigenvls_bi,
	   &returncode, maxmin, 1,
	   Qsq);
  
  *nr_of_eigenvalues = converged;

  returnvalue = eigenvls_bi[0];
  return(returnvalue);
}

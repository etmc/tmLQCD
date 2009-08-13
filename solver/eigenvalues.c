/***********************************************************************
 * $Id$
 *
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
 *
 * Autor: Carsten Urbach <urbach@ifh.de>
 *
 **************************************************************************/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "global.h"
#include "su3.h"
#include "io.h"
#include "tm_operators.h"
#include "solver/solver.h"
#include "solver/jdher.h"
#include "solver/matrix_mult_typedef.h"
#include "eigenvalues.h"


spinor  *eigenvectors = NULL;
double * eigenvls = NULL;
double max_eigenvalue;
double * inv_eigenvls = NULL;
int eigenvalues_for_cg_computed = 0;
int no_eigenvalues;

/* the folowing two are needed for the overlap */
double ev_minev=-1., ev_qnorm=-1.;

double eigenvalues(int * nr_of_eigenvalues, const int max_iterations, 
		   const double precision, const int maxmin,
		   const int readwrite, const int nstore, 
		   const int even_odd_flag) {
  double returnvalue;

#ifdef HAVE_LAPACK
  static spinor * eigenvectors_ = NULL;
  static int allocated = 0;
  char filename[200];
  FILE * ofs;
#ifdef MPI
  double atime, etime;
#endif

  /**********************
   * For Jacobi-Davidson 
   **********************/
  int verbosity = g_debug_level, converged = 0, blocksize = 1, blockwise = 0;
  int solver_it_max = 50, j_max, j_min;
  /*int it_max = 10000;*/
  /* complex *eigv_ = NULL, *eigv; */
  double decay_min = 1.7, decay_max = 1.5, prec,
    threshold_min = 1.e-3, threshold_max = 5.e-2;

  /* static int v0dim = 0; */
  int v0dim = 0;
  matrix_mult f;
  int N = (VOLUME)/2, N2 = (VOLUMEPLUSRAND)/2;
  spinor * max_eigenvector_ = NULL, * max_eigenvector;

  /**********************
   * General variables
   **********************/
  int returncode=0;
  int returncode2=0;

  no_eigenvalues = *nr_of_eigenvalues;
  if(!even_odd_flag) {
    N = (VOLUME);
    N2 = (VOLUMEPLUSRAND);
    f = &Q_pm_psi;
  }
  else {
    f = &Qtm_pm_psi;
  }

  if(g_proc_id == g_stdio_proc && g_debug_level >0) {
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
    j_min = (*nr_of_eigenvalues);
  }
  if(precision < 1.e-14){
    prec = 1.e-14;
  }
  else{
    prec = precision;
  }
#if (defined SSE || defined SSE2 || defined SSE3)
    max_eigenvector_ = calloc(N2+1, sizeof(spinor));
    max_eigenvector = (spinor *)(((unsigned long int)(max_eigenvector_)+ALIGN_BASE)&~ALIGN_BASE);
#else
    max_eigenvector_= calloc(N2, sizeof(spinor));
    max_eigenvector = max_eigenvectors;
#endif  


  if(allocated == 0) {
    allocated = 1;
#if (defined SSE || defined SSE2 || defined SSE3)
    eigenvectors_ = calloc(N2*(*nr_of_eigenvalues)+1, sizeof(spinor)); 
    eigenvectors = (spinor *)(((unsigned long int)(eigenvectors_)+ALIGN_BASE)&~ALIGN_BASE);
#else
    eigenvectors_= calloc(N2*(*nr_of_eigenvalues), sizeof(spinor));
    eigenvectors = eigenvectors_;
#endif
    eigenvls = (double*)malloc((*nr_of_eigenvalues)*sizeof(double));
    inv_eigenvls = (double*)malloc((*nr_of_eigenvalues)*sizeof(double));
  }

  solver_it_max = 50;
  /* compute the maximal one first */
  jdher(N*sizeof(spinor)/sizeof(complex), N2*sizeof(spinor)/sizeof(complex),
	50., 1.e-12, 
	1, 15, 8, max_iterations, 1, 0, 0, NULL,
	CG, solver_it_max,
	threshold_max, decay_max, verbosity,
	&converged, (complex*) max_eigenvector, (double*) &max_eigenvalue,
	&returncode2, JD_MAXIMAL, 1,
	f);

  if(readwrite && even_odd_flag) {
    for(v0dim = 0; v0dim < (*nr_of_eigenvalues); v0dim++) {
      sprintf(filename, "eigenvector.%s.%.2d.%.4d", maxmin ? "max" : "min", v0dim, nstore);
      if((read_eospinor(&eigenvectors[v0dim*N2], filename)) != 0) {
	break;
      }
    }
  }

  if(readwrite != 2) {
#ifdef MPI
    atime = MPI_Wtime();
#endif
    /* (re-) compute minimal eigenvalues */
    converged = 0;
    solver_it_max = 200;
    jdher(N*sizeof(spinor)/sizeof(complex), N2*sizeof(spinor)/sizeof(complex),
	  0., prec, 
	  (*nr_of_eigenvalues), j_max, j_min, 
	  max_iterations, blocksize, blockwise, v0dim, (complex*) eigenvectors,
	  CG, solver_it_max,
	  threshold_min, decay_min, verbosity,
	  &converged, (complex*) eigenvectors, eigenvls,
	  &returncode, JD_MINIMAL, 1,
	  f);
    
#ifdef MPI
    etime = MPI_Wtime();
    if(g_proc_id == 0) {
      printf("Eigenvalues computed in %e sec. (MPI_Wtime)\n", etime-atime);
    }
#endif
  }
  else {
    sprintf(filename, "eigenvalues.%s.%.4d", maxmin ? "max" : "min", nstore); 
    if((ofs = fopen(filename, "r")) != (FILE*) NULL) {
      for(v0dim = 0; v0dim < (*nr_of_eigenvalues); v0dim++) {
	fscanf(ofs, "%d %lf\n", &v0dim, &eigenvls[v0dim]);
	if(feof(ofs)) break;
	converged = v0dim;
      }
    }
    fclose(ofs);
  }

  (*nr_of_eigenvalues) = converged;
  no_eigenvalues = converged;
  ev_minev = eigenvls[(*nr_of_eigenvalues)-1];
  if(maxmin == JD_MINIMAL) {
    eigenvalues_for_cg_computed = converged;
  }
  if(readwrite == 1 && even_odd_flag) {
    for(v0dim = 0; v0dim < (*nr_of_eigenvalues); v0dim++) {
      sprintf(filename, "eigenvector.%s.%.2d.%.4d", maxmin ? "max" : "min", v0dim, nstore);
      if((write_eospinor(&eigenvectors[v0dim*(VOLUMEPLUSRAND)/2], filename, eigenvls[v0dim], prec, nstore)) != 0) {
	break;
      }
    }
  }
  if(g_proc_id == 0 && readwrite != 2) {
    sprintf(filename, "eigenvalues.%s.%.4d", maxmin ? "max" : "min", nstore); 
    ofs = fopen(filename, "w");
    for(v0dim = 0; v0dim < (*nr_of_eigenvalues); v0dim++) {
      fprintf(ofs, "%d %e\n", v0dim, eigenvls[v0dim]);
    }
    fclose(ofs);
  }
  for(v0dim = 0; v0dim < converged; v0dim++) {
    inv_eigenvls[v0dim] = 1./eigenvls[v0dim];
  }

  ev_qnorm=1.0/(sqrt(max_eigenvalue)+0.1);
  ev_minev*=ev_qnorm*ev_qnorm;

  returnvalue=eigenvls[0];
#else
  fprintf(stderr, "lapack not available, so JD method for EV computation not available \n");
#endif
  free(max_eigenvector_);
  return(returnvalue);
}

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

/* ************************************************************************
 * Main routine for the LapH_ev program: computes eigensystem of the Laplacian operator.
 * Authors: Luigi Scorzato, Marco Cristoforetti
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
#include <io/eospinor.h>
#include <io/params.h>
#include <io/gauge.h>
#include <io/spinor.h>
#include <io/utils.h>
#include "jacobi.h"
#include "solver/solver.h"
#include "solver/jdher_su3vect.h"
#include "solver/matrix_mult_typedef.h"
#include "linalg_eo.h"
#include "eigenvalues_Jacobi.h"
#include "gettime.h"

#ifdef WITHLAPH

su3_vector *eigenvectors_su3v = NULL;
double *eigenvls_su3v = NULL;
double max_eigenvalue_su3v;
double * inv_eigenvls_su3v = NULL;

int eigenvalues_for_cg_computed_su3v = 0;
int evlength_su3v;

double eigenvalues_Jacobi(int * nr_of_eigenvalues, const int max_iterations, 
			  const double precision, const int maxmin,int tslice, 
			  const int nstore) {
  double returnvalue;
  static int allocated = 0;

#ifdef HAVE_LAPACK


  int verbosity = 1, converged = 0, blocksize = 1 , blockwise=0;
  int solver_it_max = 50, j_max, j_min;
  double decay_min = 1.7, decay_max = 1.5, prec, threshold_min = 1.e-3, threshold_max = 5.e-2;
  int v0dim = 0;
  matrix_mult_su3vect f;
  int N=SPACEVOLUME, N2=(SPACEVOLUME + SPACERAND);
  su3_vector * max_eigenvector_ = NULL, *max_eigenvector;
  
  int returncode=0;
  int returncode2=0;
  su3_vector *s;
  double sqnorm;
  
  char filename[200];
  char eigvl_filename[200];
  //  int dims[]={T*g_nproc_t, LX*g_nproc_x, LY*g_nproc_y, LZ*g_nproc_z};
  int dims[]={1, LX*g_nproc_x, LY*g_nproc_y, LZ*g_nproc_z};
  FILE *efp;

#ifdef MPI
  double atime, etime;
  MPI_File fp;
  MPI_Offset siteSize=3*2*sizeof(double);
  LemonRecordHeader *header;
  LemonWriter *writer;
#else
  FILE *fp;
  int siteSize=3*2*sizeof(double);
#endif

  f = &Jacobi;
  evlength_su3v = N2;
  
  if(g_proc_id == g_stdio_proc && g_debug_level >0) 
    {
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
  max_eigenvector_= calloc(N2, sizeof(su3_vector));
  max_eigenvector = max_eigenvector_;
  
  if(allocated == 0) 
    {
      allocated = 1;
      eigenvectors_su3v = calloc(N2*(*nr_of_eigenvalues), sizeof(su3_vector));;
      eigenvls_su3v = (double*)malloc((*nr_of_eigenvalues)*sizeof(double));
      inv_eigenvls_su3v = (double*)malloc((*nr_of_eigenvalues)*sizeof(double));
    }
  
  solver_it_max = 64;
  /* compute the maximal one first */
  /* DEBUG 
  jdher_su3vect(N*sizeof(su3_vector)/sizeof(_Complex double), N2*sizeof(su3_vector)/sizeof(_Complex double),
		50., 1.e-12, 
		1, 15, 8, max_iterations, 1, 0, 0, NULL,
		CG, solver_it_max,
		threshold_max, decay_max, verbosity,
		&converged, (_Complex double*) max_eigenvector, (double*) &max_eigenvalue_su3v,
		&returncode2, JD_MAXIMAL, 1,tslice,f);
  */
  
  atime = gettime();
  
  /* (re-) compute minimal eigenvalues */
  converged = 0;
  solver_it_max = 256;
  
  if(maxmin)
    jdher_su3vect(N*sizeof(su3_vector)/sizeof(_Complex double), N2*sizeof(su3_vector)/sizeof(_Complex double),
		  50., prec, 
		  (*nr_of_eigenvalues), j_max, j_min, 
		  max_iterations, blocksize, blockwise, v0dim, (_Complex double*) eigenvectors_su3v,
		  CG, solver_it_max,
		  threshold_max, decay_max, verbosity,
		  &converged, (_Complex double*) eigenvectors_su3v, eigenvls_su3v,
		  &returncode, JD_MAXIMAL, 1,tslice,
		  f);
  else
    jdher_su3vect(N*sizeof(su3_vector)/sizeof(_Complex double), N2*sizeof(su3_vector)/sizeof(_Complex double),
		  0., prec, 
		  (*nr_of_eigenvalues), j_max, j_min, 
		  max_iterations, blocksize, blockwise, v0dim, (_Complex double*) eigenvectors_su3v,
		  CG, solver_it_max,
		  threshold_min, decay_min, verbosity,
		  &converged, (_Complex double*) eigenvectors_su3v, eigenvls_su3v,
		  &returncode, JD_MINIMAL, 1,tslice,
		  f);
  
  etime = gettime();
  if(g_proc_id == 0) {
    printf("Eigenvalues computed in %e sec. (gettime)\n", etime-atime);
    }

  
  /* Printout eigenvalues.  */
  if(g_proc_id == 0) {
    sprintf(eigvl_filename,"eigenvalues.%.3d.%.4d", tslice, nstore);
    efp=fopen(eigvl_filename,"w");
    for(v0dim = 0; v0dim < (*nr_of_eigenvalues); v0dim++) {
      fprintf(efp,"%e\n",eigenvls_su3v[v0dim]);
    }
    fclose(efp);    
  }

  /* Printout eigenvectors.  */
  for(v0dim = 0; v0dim < (*nr_of_eigenvalues); v0dim++) {
    sprintf(filename, "eigenvector.%.3d.%.3d.%.4d", v0dim, tslice, nstore);
    s=(su3_vector*)&eigenvectors_su3v[v0dim*N2];
#ifdef MPI 
# ifdef HAVE_LIBLEMON
    // SEGNO: dovrebbe stampare 8*2*3*SPACEVOLUME data per file, ma ne stampa 8*2*4n*SPACEVOLUME (n=4-1 per ev 0-3)

    MPI_File_open(g_cart_grid, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fp);
    writer = lemonCreateWriter(&fp, g_cart_grid);
    header = lemonCreateHeader(1 /* MB */, 1 /* ME */, "lattice-su3_vector-data",SPACEVOLUME*3*sizeof(_Complex double));
    lemonWriteRecordHeader(header, writer);
    lemonDestroyHeader(header);
    lemonWriteLatticeParallel(writer, s, siteSize, dims);
    lemonWriterCloseRecord(writer);
    lemonDestroyWriter(writer);
    MPI_File_close(&fp);
# else
  if(g_proc_id == 0) {
    printf("Cannot write eigenvectors: you need LEMON for writing eigenvectors with MPI\n");
    }
# endif
#else
    fp=fopen(filename,"wb");
    fwrite(s,siteSize,SPACEVOLUME,fp);
    fclose(fp);
#endif // MPI
    sqnorm=square_norm_su3vect(s,SPACEVOLUME,1);
    if(g_proc_id == 0) {
      printf("wrote eigenvector | |^2 = %e \n",sqnorm);
    }
  }

  returnvalue=eigenvls_su3v[0];
  free(max_eigenvector_);
#else
  fprintf(stderr, "lapack not available, so JD method for EV computation not available \n");
#endif // LAPACK
  return(returnvalue);
}

#endif // WITHLAPH

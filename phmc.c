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
/* $Id$ */

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "global.h"

#include "read_input.h"
#include "init_bispinor_field.h"
#include "eigenvalues_bi.h"
#include "solver/solver.h"
#include "init_chi_spinor_field.h"
#include "init_chi_copy.h"
#include "chebyshev_polynomial_nd.h"
#include "Ptilde_nd.h"
#include "phmc.h"

double phmc_Cpol;
double phmc_cheb_evmin, phmc_cheb_evmax;
double phmc_invmaxev;
complex * phmc_root;
int phmc_dop_n_cheby;
double * phmc_dop_cheby_coef;
int phmc_ptilde_n_cheby;
double * phmc_ptilde_cheby_coef;

void init_phmc() {
  int max_iter_ev, j, k;
  double temp, temp2;
  FILE *roots;
  char *filename_phmc_root = "Square_root_BR_roots.dat";
  char *filename_phmc_root_oox = "Square_root_BR_roots.dat.oox";
  char title[50];

  FILE *Const;
  char *filename_const = "normierungLocal.dat";
  char *filename_const_oox = "normierungLocal.dat.oox";

  /* START IF PHMC */

  phmc_invmaxev=1.0;

  if(phmc_compute_evs != 0) {
    g_mu = g_mu1;
    max_iter_ev = 1000;
    
    no_eigenvalues = 10;   /* Number of lowest eigenvalues to be computed */
    if(g_epsbar!=0.0)
      phmc_cheb_evmin = eigenvalues_bi(&no_eigenvalues, max_iter_ev, eigenvalue_precision, 0);
    else {
      phmc_cheb_evmin = eigenvalues(&no_eigenvalues, max_iter_ev, eigenvalue_precision, 0, 0, nstore, even_odd_flag);
    }

    no_eigenvalues = 4;   /* Number of highest eigenvalues to be computed */
    if(g_epsbar!=0.0)
      phmc_cheb_evmax = eigenvalues_bi(&no_eigenvalues, max_iter_ev, eigenvalue_precision, 1);
    else
      phmc_cheb_evmax = eigenvalues(&no_eigenvalues, max_iter_ev, eigenvalue_precision, 1, 0, nstore, even_odd_flag);
       
    temp=phmc_cheb_evmin;
    temp2=phmc_cheb_evmax;
    
    if(g_proc_id==0) {
      printf("PHMC: Ev-max = %e \n", phmc_cheb_evmax);
      printf("PHMC: Ev-min = %e \n", phmc_cheb_evmin); 
    }
#ifdef MPI
    MPI_Finalize();
#endif
    exit(0);
  }

  /* This is the epsilon parameter */
  phmc_cheb_evmin = stilde_min/(stilde_max);

  /* In the following there is the  "sqrt"  since the value refers to 
     the hermitian Dirac operator (used in EV-computation), namely 
     S = Q Q^dag         
     When  "S"  is applied, we call  phmc_invmaxev  twice !!! */
  if(g_epsbar!=0.0 || phmc_exact_poly==0) phmc_invmaxev=1./(sqrt(stilde_max));
  else if(g_epsbar==0.0 && phmc_exact_poly==1) phmc_invmaxev=1./stilde_max;
  phmc_cheb_evmax = 1.0;

  /* Here we prepare the less precise polynomial first */
  degree_of_polynomial_nd(degree_of_p);

  if((g_proc_id == 0) && (g_debug_level > 1)) {
    printf("PHMC: interval of approximation [stilde_min, stilde_max] = [%e, %e]\n", stilde_min, stilde_max);
    printf("PHMC: degree for P = %d, epsilont = %e, normalisation = %e", 
	   phmc_dop_n_cheby-1, phmc_cheb_evmin, phmc_invmaxev);
  }

  /* Chi`s-spinors  memory allocation */
  j = init_chi_up_spinor_field(VOLUMEPLUSRAND/2, (phmc_dop_n_cheby+1));
  if ( j!= 0) {
    fprintf(stderr, "Not enough memory for PHMC Chi_up fields! Aborting...\n");
    exit(0);
  }
  j = init_chi_dn_spinor_field(VOLUMEPLUSRAND/2, (phmc_dop_n_cheby+1));
  if ( j!= 0) {
    fprintf(stderr, "Not enough memory for PHMC Chi_dn fields! Aborting...\n");
    exit(0);
  }
  /* End memory allocation */
  /* Here we prepare the precise polynomial */
  degree_of_Ptilde();

  /* THIS IS THE OVERALL CONSTANT */
  /* write phmc_Cpol as the result of the simple-program files (BigC^(1/2))^1/2 
     since  BigC^(1/2)  is the constant appearing in each factor of the 
     multiplication defining the monomial basis representation of the 
     polinomial in s,  while its square phmc_root  (BigC^(1/2))^1/2  is the 
     constant appearing in the multiplication representing the 
     polinomial in  sqrt(s) .
  */
  if(g_epsbar==0.0 && phmc_exact_poly ==1) 
    filename_const=filename_const_oox;
  if((Const=fopen(filename_const,"r")) != (FILE*)NULL) {
    fscanf(Const, " %lf \n", &phmc_Cpol);
    fclose(Const);
  }
  else {
    fprintf(stderr, "File %s is missing! Aborting...\n", filename_const);
#ifdef MPI
    MPI_Finalize();
#endif
    exit(6);
  }
  if(g_epsbar!=0.0 || phmc_exact_poly==0) phmc_Cpol = sqrt(phmc_Cpol);

  phmc_root = calloc((2*phmc_dop_n_cheby-2),sizeof(complex));


  if(g_epsbar==0.0 && phmc_exact_poly == 1) 
    filename_phmc_root=filename_phmc_root_oox;
  if((roots=fopen(filename_phmc_root,"r")) != (FILE*)NULL) {
    fgets(title, 100, roots);
    
    /* Here we read in the 2n roots needed for the polinomial in sqrt(s) */
    for(j=0; j<(2*phmc_dop_n_cheby-2); j++){
      fscanf(roots," %d %lf %lf \n", &k, &phmc_root[j].re, &phmc_root[j].im);
    }
    fclose(roots);
  }
  else {
    fprintf(stderr, "File %s is missing! Aborting ...\n", filename_phmc_root);
#ifdef MPI
    MPI_Finalize();
#endif
    exit(6);
  }
  
  /* END IF PHMC */
  return;
}


void phmc_compute_ev(const int trajectory_counter,
		     const double plaquette_energy) {
  double atime, etime, temp=0., temp2=0.;
  int max_iter_ev, no_eigenvalues;
  char * phmcfilename = "phmc.data";
  FILE * countfile;

#ifdef MPI
  atime = MPI_Wtime();
#else
  atime = (double)clock()/(double)(CLOCKS_PER_SEC);
#endif
  max_iter_ev = 1000;
  g_mu = g_mu1;
  
  no_eigenvalues = 1;

  if(g_epsbar!=0.0)
    temp = eigenvalues_bi(&no_eigenvalues, max_iter_ev, eigenvalue_precision, 0);
  else
    temp = eigenvalues(&no_eigenvalues, max_iter_ev, eigenvalue_precision, 0, 0, nstore, even_odd_flag);
  
  no_eigenvalues = 1;
  if(g_epsbar!=0.0)
    temp2 = eigenvalues_bi(&no_eigenvalues, max_iter_ev, eigenvalue_precision, 1);
  else
    temp2 = eigenvalues(&no_eigenvalues, max_iter_ev, eigenvalue_precision, 1, 0, nstore, even_odd_flag);
  
  if((g_proc_id == 0) && (g_debug_level > 0)) {
    printf("PHMC: lowest eigenvalue end of trajectory %d = %e\n", 
	   trajectory_counter, temp);
    printf("PHMC: maximal eigenvalue end of trajectory %d = %e\n", 
	   trajectory_counter, temp2);
  }
  if(g_proc_id == 0) {
    countfile = fopen(phmcfilename, "a");
    fprintf(countfile, "%d %1.12f %1.5e %1.5e %1.5e %1.5e\n", 
	    trajectory_counter, plaquette_energy/(6.*VOLUME*g_nproc), temp, temp2, stilde_min, stilde_max);
    fclose(countfile);
  }
#ifdef MPI
  etime = MPI_Wtime();
#else
  etime = (double)clock()/(double)(CLOCKS_PER_SEC);
#endif
  if((g_proc_id == 0)) {
    printf("PHMC: time/s for eigenvalue computation %e\n", etime-atime);
  }
}

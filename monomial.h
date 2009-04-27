/***********************************************************************
 * $Id$ 
 *
 * Copyright (C) 2008 Carsten Urbach
 *
 * Modified by Jenifer Gonzalez Lopez 2009/03/27
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
#ifndef _MONOMIAL_H
#define _MONOMIAL_H


#define DET 0
#define DETRATIO 1
#define GAUGE 2
#define POLY 3
#define NDPOLY 4
#define SFGAUGE 5
#define NDDETRATIO 6

#define max_no_monomials 20

typedef struct {
  int type;
  int gtype;
  int initialised;
  int timescale;
  int maxiter;
  int id;
  int even_odd_flag;
  int rngrepro;
  int solver;
  int iter0, iter1, iter2;
  int csg_N, csg_N2;
  int csg_n, csg_n2;
  int use_rectangles;
  int * csg_index_array, *csg_index_array2;
  /* det or detratio related */
  double mu, mu2, kappa, kappa2;
  /* polynomial related, not yet in use */
  double mubar, epsbar, mubar2, epsbar2;
  /* energies at beginning and end of trajectory */
  double energy0; 
  double energy1;
  /* gauge related */
  double c0, c1, beta;
  /* solver related*/
  double epsilon;
  double forceprec;
  double accprec;
  /* force normalisation */
  double forcefactor;
  /* sf */          
  double eta;             
  double ct, cs; 
  double c1ss, c1tss, c1tts;
  /* some book-keeping */
  char name[100];
  /* pseudo fermion field */
  /* second one needed for ND monomials */
  spinor * pf, * pf2;
  /* chronological solver fields */
  spinor ** csg_field;
  spinor ** csg_field2;
  /* functions for the HMC update */
  void (*hbfunction) (const int no);
  double (*accfunction) (const int no);
  void (*derivativefunction) (const int no);
} monomial;

#include "su3.h"
#include "su3spinor.h"
#include "det_monomial.h"
#include "detratio_monomial.h"
#include "ndpoly_monomial.h"
#include "nddetratio_monomial.h"
#include "gauge_monomial.h"
#include "sf_gauge_monomial.h"

/* list of all monomials */
extern monomial monomial_list[max_no_monomials];
/* number of initialised monomials */
extern int no_monomials;
/* number of gauge monomials, currently 0 or 1 */
extern int no_gauge_monomials;
/* number of ndpoly monomials, currently 0 or 1 */
extern int no_ndpoly_monomials;

/* add a new monomial to the list of monomials */
int add_monomial(const int type);
/* initialise all monomials in the list */
int init_monomials(const int V, const int even_odd_flag);
/* free space again */
void free_monomials();

/* some dummy functions */
void dummy_derivative(const int id);
void dummy_heatbath(const int id);
double dummy_acc(const int id);

#endif

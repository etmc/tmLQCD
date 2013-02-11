/***********************************************************************
 *
 * Copyright (C) 2008,2011,2012 Carsten Urbach
 *               2009 Jenifer Gonzalez Lopez
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

#include "su3.h"
#include "su3spinor.h"
#include "hamiltonian_field.h"
#include "rational/rational.h"

#define DET 0
#define DETRATIO 1
#define GAUGE 2
#define POLY 3
#define NDPOLY 4
#define SFGAUGE 5
#define NDDETRATIO 6
#define POLYDETRATIO 7
#define CLOVERTRLOG 8
#define CLOVERDET 9
#define CLOVERDETRATIO 10
#define NDCLOVER 11
#define CLOVERNDTRLOG 12
#define NDRAT 13
#define NDCLOVERRAT 14
#define NDRATCOR 15
#define NDCLOVERRATCOR 16
#define RAT 17
#define RATCOR 18
#define CLOVERRAT 19
#define CLOVERRATCOR 20

#define max_no_monomials 30

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
  /* trlog */
  int trlog;
  int * csg_index_array, *csg_index_array2;
  /* det or detratio related */
  double mu, mu2, kappa, kappa2;
  /* clover coefficient */
  double c_sw, rho, rho2;
  /* polynomial related, not yet in use */
  double mubar, epsbar, mubar2, epsbar2;
  /* energies at beginning and end of trajectory */
  double energy0; 
  double energy1;
  /* gauge related */
  double c0, c1, beta, glambda;
  /* solver related*/
  double epsilon;
  double forceprec;
  double accprec;
  /* force normalisation */
  double forcefactor;
  /* some book-keeping */
  char name[100];
  /* pseudo fermion field */
  /* second one needed for ND monomials */
  spinor * pf, * pf2;
  /* parameters for the POLY Monomial*/
  int rec_ev;
  int MDPolyDegree, MaxPtildeDegree, PtildeDegree;
  double MDPolyLmin, MDPolyLmax;
  char MDPolyRootsFile[256];
  _Complex double *MDPolyRoots;
  spinor **MDPoly_chi_spinor_fields;
  double MDPolyLocNormConst;
  int MDPolyDetRatio;
  int no_wfields;
  double PrecisionPtilde;
  double PrecisionHfinal;
  double StildeMin, StildeMax;
  double EVMin, EVMax, EVMaxInv;
  double * MDPolyCoefs, * PtildeCoefs;
  /* rational approximation */
  rational_t rat;
  /* chronological solver fields */
  spinor ** csg_field;
  spinor ** csg_field2;
  spinor ** w_fields;
  /* functions for the HMC update */
  void (*hbfunction) (const int no, hamiltonian_field_t * const hf);
  double (*accfunction) (const int no, hamiltonian_field_t * const hf);
  void (*derivativefunction) (const int no, hamiltonian_field_t * const hf);
  /* the operator definitions */
  void (*Qsq) (spinor * const, spinor * const);
  void (*Qp) (spinor * const, spinor * const);
  void (*Qm) (spinor * const, spinor * const);
} monomial;

#include "monomial/det_monomial.h"
#include "monomial/detratio_monomial.h"
#include "monomial/poly_monomial.h"
#include "monomial/ndpoly_monomial.h"
#include "monomial/nddetratio_monomial.h"
#include "monomial/gauge_monomial.h"
#include "monomial/sf_gauge_monomial.h"
#include "monomial/clover_trlog_monomial.h"
#include "monomial/clovernd_trlog_monomial.h"
#include "monomial/cloverdet_monomial.h"
#include "monomial/cloverdetratio_monomial.h"
#include "monomial/cloverndpoly_monomial.h"
#include "monomial/ndrat_monomial.h"
#include "monomial/rat_monomial.h"
#include "monomial/ndratcor_monomial.h"
#include "monomial/ratcor_monomial.h"
#include "monomial/moment_energy.h"
#include "monomial/monitor_forces.h"

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

/* initialisation function for a poly monomial */
int init_poly_monomial(const int V,const int id);


/* some dummy functions */
void dummy_derivative(const int id, hamiltonian_field_t * const hf);
void dummy_heatbath(const int id, hamiltonian_field_t * const hf);
double dummy_acc(const int id, hamiltonian_field_t * const hf);

#endif

#ifndef _MONOMIAL_H
#define _MONOMIAL_H


#define DET 0
#define DETRATIO 1
#define GAUGE 2
#define POLY 3
#define NDPOLY 4

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
  double mu;
  double mu2;
  double kappa;
  double kappa2;
  double epsilon;
  double forceprec;
  double accprec;
  double forcefactor;
  double energy0; 
  double energy1;
  double c0, c1;
  double beta;
  char name[100];
  spinor * pf;
  spinor ** csg_field;
  spinor ** csg_field2;
  void (*hbfunction) (const int no);
  double (*accfunction) (const int no);
  void (*derivativefunction) (const int no);
} monomial;

#include "su3.h"
#include "su3spinor.h"
#include "det_monomial.h"
#include "detratio_monomial.h"
#include "ndpoly_monomial.h"
#include "gauge_monomial.h"

extern monomial monomial_list[20];
extern int no_monomials;
extern int no_gauge_monomials;
extern int no_ndpoly_monomials;

int add_monomial(const int type);
int init_monomials(const int V, const int even_odd_flag);
void free_monomials();

void dummy_heatbath(const int id);
double dummy_acc(const int id);

#endif

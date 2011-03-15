
#include "su3.h"

extern int x_n_cheby;
extern double * x_cheby_coef;

void norm_X_sqr_psi(spinor * const R, spinor * const S, double const mstar);

void norm_X_n_psi(spinor * const R, spinor * const S, const int n, double const mstar);

void X_over_sqrt_X_sqr(spinor * const R, double * const c, const int n, spinor * const S, const double minev, double const mstar);

void h_X_sqr_eta(spinor * const R1,spinor * const R2,spinor * const S, double const mstar);

void h_X_eta(spinor * const R,spinor * const S, double const mstar);

void h_X_4_eta(spinor * const R1, spinor * const R2, spinor * const S, double const mstar);

void Check_Approximation(double const mstar);



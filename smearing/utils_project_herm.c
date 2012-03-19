#include "utils.ih"

/* This implements the approach of taking (I/2) * (Omega' - Omega) - (I/6) * Tr(Omega' - Omega) */

void project_herm(su3 *omega)
{
  static const double fac_3 = 1.00 / 3.00;
  double tr_omega = fac_3 * (cimag(omega->c00) + cimag(omega->c11) + cimag(omega->c22));

  omega->c00 = cimag(omega->c00) - tr_omega;
  omega->c11 = cimag(omega->c11) - tr_omega;
  omega->c22 = cimag(omega->c22) - tr_omega;

  omega->c01  = 0.5 * (omega->c10 - conj(omega->c01));
  omega->c10 = conj(omega->c01);

  omega->c02  = 0.5 * (omega->c20 - conj(omega->c02));
  omega->c20 = conj(omega->c02);
  
  omega->c21  = 0.5 * (omega->c12 - conj(omega->c21));
  omega->c12 = conj(omega->c21);
}

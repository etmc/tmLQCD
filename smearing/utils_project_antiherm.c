#include "utils.ih"

void project_antiherm(su3 *omega)
{
  static const double fac_3 = 1.00 / 3.00;
  double tr_omega = creal(-I * fac_3 * (omega->c00 + omega->c11 + omega->c22);

  
  omega->c00 = (cimag(omega->c00) - tr_omega) * I;
  omega->c11 = (cimag(omega->c11) - tr_omega) * I;
  omega->c22 = (cimag(omega->c22) - tr_omega) * I;

  omega->c01 -= conj(omega->c10);
  omega->c01 *= 0.50;
  omega->c10  = -conj(omega->c01);


  omega->c02 -= conj(omega->c20);
  omega->c02 *= 0.50;
  omega->c20  = -conj(omega->c02);

  omega->c12 -= conj(omega->c21);
  omega->c12 *= 0.50;
  omega->c21  = -conj(omega->c12);
}

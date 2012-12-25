#include "utils.ih"

void project_traceless_antiherm(su3 *in)
{
  static const double fac_3 = 1.00 / 3.00;
  double tr_in = fac_3 * (cimag(in->c00) + cimag(in->c11) + cimag(in->c22));
  
  in->c00  = (cimag(in->c00) - tr_in) * I;
  in->c11  = (cimag(in->c11) - tr_in) * I;
  in->c22  = (cimag(in->c22) - tr_in) * I;

  in->c01 -= conj(in->c10);
  in->c01 *= 0.50;
  in->c10  = -conj(in->c01);

  in->c02 -= conj(in->c20);
  in->c02 *= 0.50;
  in->c20  = -conj(in->c02);

  in->c12 -= conj(in->c21);
  in->c12 *= 0.50;
  in->c21  = -conj(in->c12);
}

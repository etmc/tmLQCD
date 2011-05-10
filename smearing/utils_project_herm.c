#include "utils.ih"

/* This implements the approach of taking (I/2) * (Omega' - Omega) - (I/6) * Tr(Omega' - Omega) */

void project_herm(su3 *omega)
{
  static const double fac_3 = 1.00 / 3.00;
  static double tmp;
  double tr_omega = (omega->c00.im + omega->c11.im + omega->c22.im) * fac_3;

  omega->c00.re = omega->c00.im - tr_omega;
  omega->c00.im = 0.00;

  omega->c11.re = omega->c11.im - tr_omega;
  omega->c11.im = 0.00;

  omega->c22.re = omega->c22.im - tr_omega;
  omega->c22.im = 0.00;

  tmp = 0.5 * (omega->c10.im + omega->c01.im);
  omega->c01.im = 0.5 * (omega->c10.re - omega->c01.re);
  omega->c01.re = tmp;
  omega->c10.im = -omega->c01.im;
  omega->c10.re = tmp;

  tmp = 0.5 * (omega->c02.im + omega->c20.im);
  omega->c02.im = 0.5 * (omega->c20.re - omega->c02.re);
  omega->c02.re = tmp;
  omega->c20.im = -omega->c02.im;
  omega->c20.re = tmp;
  
  tmp = 0.5 * (omega->c21.im + omega->c12.im);
  omega->c21.im = 0.5 * (omega->c12.re - omega->c21.re);
  omega->c21.re = tmp;
  omega->c12.im = -omega->c21.im;
  omega->c12.re = tmp;
}

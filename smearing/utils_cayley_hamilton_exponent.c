#include "utils.ih"

void cayley_hamilton_exponent(su3 *U, _Complex double *f0, _Complex double *f1, _Complex double *f2, su3 *Q)
{
  static double const fac_1_3 = 1 / 3.0;
  
  double c0 = Q->c00 * (Q->c11 * Q->c22 - Q->c12 * Q->c21) + 
              Q->c01 * (Q->c12 * Q->c20 - Q->c10 * Q->c22) +
              Q->c02 * (Q->c10 * Q->c21 - Q->c11 * Q->c20)  ;
  double sign_c0 = copysign(1.0, c0);
  c0 = fabs(c0);

  /* We'll need U^2 down the line, so we might as well calculate it now. */
  /* This is also why we need the buffer space -- can't do matrix multiplication in place */
  _su3_times_su3(*U, *Q, *Q);
  double c1 = 0.5 * (U->c00 + U->c11 + U->c22);

  double c0max = 2.0 * pow(fac_1_3 * c1, 1.5);
  double theta_3 = fac_1_3 * acos(c0 / c0max); 

  double u = sqrt(fac_1_3 * c1) * cos(theta_3);
  double w = sqrt(c1) * sin(theta_3);
  
  /* Modification w.r.t. Peardon & Morningstar: w is always positive, so |w| = w */
  double xi0 = (w > 0.05) ? (sin(w) / w) : 1 - 0.16666666666666667 * w * w * (1 - 0.05 * w * w * (1 - 0.023809523809523808 * w * w));
  double divisor = 1.0 / (9.0 * u * u - w * w);
  
  *f0 = divisor * ((u * u - w * w) * cexp(2 * I * u) + cexp(-I * u) * (8 * u * u * cos(w) + 2 * I * u * (3 * u * u + w * w) * xi0));
  *f1 = divisor * (2 * u * cexp(2 * I * u) - cexp(-I * u) * (2 * u * cos(w) - I * (3 * u * u - w * w) * xi0));
  *f2 = divisor * (cexp(2 * I * u) - cexp(-I * u) * (cos(w) + 3 * I * u * xi0));

  if (sign_c0 < 0)
  {
    *f0 =  conj(*f0);
    *f1 = -conj(*f1);
    *f2 =  conj(*f2);
  }

  U->c00 = *f0 + *f1 * Q->c00 + *f2 * U->c00;
  U->c01 =       *f1 * Q->c01 + *f2 * U->c01;
  U->c02 =       *f1 * Q->c02 + *f2 * U->c02;
  U->c10 =       *f1 * Q->c10 + *f2 * U->c10;
  U->c11 = *f0 + *f1 * Q->c11 + *f2 * U->c11;
  U->c12 =       *f1 * Q->c12 + *f2 * U->c12;
  U->c20 =       *f1 * Q->c20 + *f2 * U->c20;
  U->c21 =       *f1 * Q->c21 + *f2 * U->c21;
  U->c22 = *f0 + *f1 * Q->c22 + *f2 * U->c22;
}

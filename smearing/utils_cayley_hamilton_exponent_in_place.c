#include "utils.ih"

void cayley_hamilton_exponent_in_place(su3 *omega, _Complex double *f0, _Complex double *f1, _Complex double *f2)
{
  static double const fac_1_3 = 1 / 3.0;
  
  static su3 tmp;
  
  double c0 = omega->c00 * omega->c11 * omega->c22 - 
              omega->c00 * omega->c12 * omega->c21 + 
              omega->c01 * omega->c12 * omega->c20 - 
              omega->c01 * omega->c10 * omega->c22 +
              omega->c02 * omega->c10 * omega->c21 - 
              omega->c02 * omega->c11 * omega->c20;
  double sign_c0 = copysign(1.0, c0);
  c0 = fabs(c0);

  /* We'll need omega^2 down the line, so we might as well calculate it now. */
  /* This is also why we need the buffer space -- can't do matrix multiplication in place */
  _su3_times_su3(tmp, *omega, *omega);
  double c1 = 0.5 * (tmp.c00 + tmp.c11 + tmp.c22);

  double c0max = 2.0 * pow(fac_1_3 * c1, 1.5);
  double theta_3 = fac_1_3 * acos(c0 / c0max); 

  double u = sqrt(fac_1_3 * c1) * cos(theta_3);
  double w = sqrt(c1) * sin(theta_3);
  
  /* Modification w.r.t. Peardon & Morningstar: w is always positive, so |w| = w */
  double xi0 = (w > 0.05) ? sin(w) / w : 1 - 0.16666666666666667 * w * w * (1 - 0.05 * w * w * (1 - 0.023809523809523808 * w * w));
  double divisor = 1.0 / (9.0 * u * u - w * w);
  
  *f0 = divisor * ((u * u - w * w) * cexp(2 * I * u) + cexp(-I * u) * (8 * u * u * cos(w) + 2 * I * u * (3 * u * u + w * w) * xi0));
  *f1 = divisor * (2 * u * cexp(2 * I * u) + cexp(-I * u) * (2 * u * cos(w) - I * (3 * u * u - w * w) * xi0));
  *f2 = divisor * (cexp(2 * I * u) + cexp(-I * u) * (cos(w) + 3 * I * u * xi0));

  if (sign_c0 < 0)
  {
    *f0 =  conj(*f0);
    *f1 = -conj(*f1);
    *f2 =  conj(*f2);
  }

  omega->c00 = *f0 + *f1 * omega->c00 + *f2 * tmp.c00;
  omega->c01 = *f0 + *f1 * omega->c01 + *f2 * tmp.c01;
  omega->c02 = *f0 + *f1 * omega->c02 + *f2 * tmp.c02;
  omega->c10 = *f0 + *f1 * omega->c10 + *f2 * tmp.c10;
  omega->c11 = *f0 + *f1 * omega->c11 + *f2 * tmp.c11;
  omega->c12 = *f0 + *f1 * omega->c12 + *f2 * tmp.c12;
  omega->c20 = *f0 + *f1 * omega->c20 + *f2 * tmp.c20;
  omega->c21 = *f0 + *f1 * omega->c21 + *f2 * tmp.c21;
  omega->c22 = *f0 + *f1 * omega->c22 + *f2 * tmp.c22;
}

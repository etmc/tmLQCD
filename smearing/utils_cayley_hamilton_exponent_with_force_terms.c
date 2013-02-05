#include "utils.ih"

#include "utils_exponent_from_coefficients.static"

/* This is a convenience function with a fairly ugly interface -- sorry about that. */
void cayley_hamilton_exponent_with_force_terms(su3* expA, su3 *B1, su3 *B2, _Complex double *f1, _Complex double *f2, su3 const *A)
{
  static double const fac_1_3 = 1 / 3.0;
    
  /* The value of f0 is never needed beyond this function, unlike f1 and f2. We make scoped room for all three,
   * in case f1 and f2 are not requested.  NOTE We want to check the performance impact of this -- this function is
   * called on the full volume of a lattice, after all. */
  _Complex double f0[3];
  
  if (!f1)
  {
    f1 = f0 + 1;
    f2 = f0 + 2;
  }
  
  /* NOTE The expressions below are adapted from Peardon-Morningstar. Note that there is a factor -I between A and Q! */
  
  /* c0 = det[A] */
  double c0 = I * (A->c00 * (A->c11 * A->c22 - A->c12 * A->c21) + 
                   A->c01 * (A->c12 * A->c20 - A->c10 * A->c22) +
                   A->c02 * (A->c10 * A->c21 - A->c11 * A->c20)  );
  
  /* c1 = 0.5 * Tr[AA] */
  double c1 = -0.5 * (A->c00 * A->c00 + A->c01 * A->c10 + A->c02 * A->c20 +
                      A->c10 * A->c01 + A->c11 * A->c11 + A->c12 * A->c21 +
                      A->c20 * A->c02 + A->c21 * A->c12 + A->c22 * A->c22  );

  /* There is a special, but common (cold start) case where the given matrix is actually 0!
   * We need to account for it. */
  if (c0 == 0 && c1 == 0) 
  {
    _su3_one(*expA);
    _su3_zero(*B1);
    _su3_zero(*B2); /* FIXME Not quite sure about this one, check the limit for A->0 explicitly. */
    *f1 = I;
    *f2 = -0.5;
    return;
  }

  /* P&M give symmetry relations that can be used when c0 < 0, to avoid the numerically problematic c0 -> -c0_max limit.
     We note the sign here for future reference, then continue with c0 as if it were positive. */
  int c0_negative = (c0 < 0);
  c0 = fabs(c0);

  /* The call to fmin below is needed, because for small deviations alpha from zero -- O(10e-12) -- rounding errors can cause c0 > c0max by epsilon.
     In that case, acos(c0/c0max) will produce NaNs, whereas the mathematically correct conclusion would be that theta is zero to machine precision! 
     Note that this approach will *not* produce identity and zero for all output, but rather the correct answer of order (I + alpha) for exp(iQ). */
  
  double c0max = 2.0 * pow(fac_1_3 * c1, 1.5);
  double theta_3 = fac_1_3 * acos(fmin(c0 / c0max, 1.0));

  double u = sqrt(fac_1_3 * c1) * cos(theta_3);
  double w = sqrt(c1) * sin(theta_3);
  
  /* Calculate and cache some repeating factors. *
   * We can fold in the sign immediately -- c.f. f_j(-c0, c1) = -1^j * conj(f_j(c0, c1)) 
   * This should just amount to potentially adding a minus to all imaginary components and an overall phase for f1. */
  _Complex double ma = cexp(2 * I * u);
  _Complex double mb = cexp(-I * u);
  double cw = cos(w);
  double u2 = u * u;
  double w2 = w * w;
  
  /* Modification w.r.t. Peardon & Morningstar:  w is always positive, so |w| =  w */
  double xi0 = (w > 0.05) ? (sin(w) / w) 
                          : 1 - 0.16666666666666667 *  w2 * (1 - 0.05 *  w2 * (1 - 0.023809523809523808 *  w2));
  double divisor = 1.0 / (9.0 * u2 -  w2);

  *f0 = divisor * (ma * (u * u -  w * w) + mb * (8 * u * u * cw + 2 * I * u * (3 * u * u +  w * w) * xi0));
  *f1 = divisor * (2 * u * ma - mb * (2 * u * cw - I * (3 * u * u -  w * w) * xi0));
  *f2 = divisor * (ma - mb * (cw + 3 * I * u * xi0));

  /* The first point where we use the symmetry relations to calculate the negative c0 possibility */
  if (c0_negative)
  {
    *f0 = conj(*f0);
    *f1 = -conj(*f1);
    *f2 = conj(*f2);
  }
  
  exponent_from_coefficients(expA,  *f0, *f1, *f2, A);
  
  /* To also provide a replacement for plain su3_expo, we add a check on a request for B1 and B2 calculation and quit if they're not needed. */
  if (!B1)
    return;
  
  double xi1 = (w > 0.05) ? (w * cw - sin(w)) / (w * w2) 
                          : -fac_1_3 * (1 - 0.10 * w2 * (1 - 0.03571428571428571 * w2 * (1 - 0.01851851851851852 * w2)));
  
  _Complex double r10 = 2 * (u + I * (u2 - w2)) * ma + 2 * mb * ((4 * u * (2 - I * u) * cw) + I * (9 * u2 + w2 - I * u * (3 * u2 + w2)) * xi0);
  _Complex double r11 = 2 * (1 + 2 * I * u) * ma + mb * (-2 * (1 - I * u) * cw + I * (6 * u + I * (w2 - 3 * u2)) * xi0);
  _Complex double r12 = 2 * I * ma + I * mb * (cw - 3 * (1 - I * u) * xi0);
  _Complex double r20 = - 2 * ma + 2 * I * u * mb * (cw + (1 + 4 * I * u) * xi0 + 3 * u2 * xi1);
  _Complex double r21 = -I * mb * (cw + (1 + 2 * I * u) * xi0 - 3 * u2 * xi1);
  _Complex double r22 = mb * (xi0 - 3 * I * u * xi1);

  divisor = 0.5 * divisor * divisor;
  double mr2 = (3 * u2 - w2);
  double mf  = 2 * (15 * u2 + w2);
  
  _Complex double bn0 = divisor * (2 * u * r10 + mr2 * r20 - mf * *f0);
  _Complex double bn1 = divisor * (2 * u * r11 + mr2 * r21 - mf * *f1);
  _Complex double bn2 = divisor * (2 * u * r12 + mr2 * r22 - mf * *f2);
  
  if (c0_negative)
  {
    bn0 = conj(bn0);
    bn1 = -conj(bn1);
    bn2 = conj(bn2);
  }
  
  exponent_from_coefficients(B1,    bn0, bn1, bn2, A);
  
  bn0 = divisor * (r10 - 3 * u * r20 - 24 * u * *f0);
  bn1 = divisor * (r11 - 3 * u * r21 - 24 * u * *f1);
  bn2 = divisor * (r12 - 3 * u * r22 - 24 * u * *f2);
  
  if (c0_negative)
  {
    bn0 = -conj(bn0);
    bn1 = conj(bn1);
    bn2 = -conj(bn2);
  }
  
  exponent_from_coefficients(B2, bn0, bn1, bn2, A);
}

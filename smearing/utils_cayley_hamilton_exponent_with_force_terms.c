#include "utils.ih"

#include "utils_exponent_from_coefficients.static"

/* This is a convenience function with a fairly ugly interface -- sorry about that. */
void cayley_hamilton_exponent_with_force_terms(su3* expA, su3 *B1, su3 *B2, double *f1, double *f2, su3 const *A)
{
  dump_su3(A, stderr);
  static double const fac_1_3 = 1 / 3.0;
    
  /* The value of f0 is never needed beyond this function, unlike f1 and f2. We make scoped room for all three,
   * in case f1 and f2 are not requested.  NOTE We want to check the performance impact of this -- this function is
   * called on the full volume of a lattice, after all. */
  static double f0[3];
  
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
  double sign_c0 = copysign(1.0, c0);
  
  c0 = fabs(c0);
  
  /* c1 = 0.5 * Tr[AA] */
  double c1 = -0.5 * (A->c00 * A->c00 + A->c01 * A->c10 + A->c02 * A->c20 +
                      A->c10 * A->c01 + A->c11 * A->c11 + A->c12 * A->c21 +
                      A->c20 * A->c02 + A->c21 * A->c12 + A->c22 * A->c22  );

  double c0max = 2.0 * pow(fac_1_3 * c1, 1.5);
  double theta_3 = fac_1_3 * acos(c0 / c0max); 

  double u = sqrt(fac_1_3 * c1) * cos(theta_3);
  double w = sqrt(c1) * sin(theta_3);
  
  /* Modification w.r.t. Peardon & Morningstar:  w is always positive, so |w| =  w */
  double xi0 = (w > 0.05) ? (sin(w) / w) 
                          : 1 - 0.16666666666666667 *  w *  w * (1 - 0.05 *  w *  w * (1 - 0.023809523809523808 *  w * w));
  double xi1 = (w > 0.05) ? (w * cos(w) - sin(w)) / (w * w * w) 
                          : -fac_1_3 * (1 - 0.10 * w * w * (1 - 0.03571428571428571 * w * w * (1 - 0.01851851851851852 * w * w)));
  double divisor = 1.0 / (9.0 * u * u -  w * w);
  
  /* NOTE The following can probably be reorganized for efficiency -- some of the factors keep returning. */
  
  double r10 = 2 * (u + I * (u * u - w * w)) * cexp(2 * I * u) + 2 * cexp(-I * u) * (4 * u * (2 - I * u) * cos(w) + I * (9 * u * u + w * w - I * u * (3 * u * u + w * w)) * xi0);
  double r11 = 2 * (1 + 2 * I * u) * cexp(2 * I * u) + cexp(-I * u) * (-2 * (1 - I * u) * cos(w) + I * (6 * u + I * (w * w - 3 * u * u)) * xi0);
  double r12 = 2 * I * cexp(2 * I * u) + I * cexp(-I * u) * (cos(w) - 3 * (1 - I * u) * xi0);
  double r20 = - 2 * cexp(2 * I * u) + 2 * I * u * cexp(-I * u) * (cos(w) + (1 + 4 * I * u) * xi0 + 3 * u * u * xi1);
  double r21 = -I * cexp(-I * u) * (cos(w) + (1 + 2 * I * u) * xi0 - 3 * u * u * xi1);
  double r22 = cexp(-I * u) * (xi0 - 3 * I * u * xi1);
  
  /* We can fold in the sign immediately -- c.f. f_j(-c0, c1) = -1^j * conj(f_j(c0, c1)) */
  *f0 = divisor * ((u * u -  w * w) * cexp(sign_c0 * 2 * I * u) + cexp(-sign_c0 * I * u) * (8 * u * u * cos(w) + 2 * sign_c0 * I * u * (3 * u * u +  w * w) * xi0));
  *f1 = sign_c0 * divisor * (2 * u * cexp(sign_c0 * 2 * I * u) - cexp(-sign_c0 * I * u) * (2 * u * cos(w) - sign_c0 * I * (3 * u * u -  w * w) * xi0));
  *f2 =  divisor * (cexp(2 * sign_c0 * I * u) - cexp(-sign_c0 * I * u) * (cos(w) + 3 * sign_c0 * I * u * xi0));
  
  exponent_from_coefficients(expA,  *f0, *f1, *f2, A);
  fprintf(stderr, "%8.6f, %8.6f, %8.6f\n",  *f0, *f1, *f2);
  dump_su3(expA, stderr);
  fprintf(stderr, "----------------\n");
  
  /* To also provide a replacement for plain su3_expo, we add a check on a request for B1 and B2 calculation and quit if they're not needed. */
  if (!B1)
    return;
    
  double bdiv = 0.5 * divisor * divisor;
  
  double bn0 = bdiv * (2 * u * r10 + (3 * u * u - w * w) * r20 - 2 * (15 * u * u + w * w) * *f0);
  double bn1 = bdiv * (2 * u * r11 + (3 * u * u - w * w) * r21 - 2 * (15 * u * u + w * w) * *f1);
  double bn2 = bdiv * (2 * u * r12 + (3 * u * u - w * w) * r22 - 2 * (15 * u * u + w * w) * *f2);
  
  exponent_from_coefficients(B1,    bn0, bn1, bn2, A);
  
  bn0 = bdiv * (r10 - 3 * u * r20 - 24 * u * *f0);
  bn1 = bdiv * (r11 - 3 * u * r21 - 24 * u * *f1);
  bn2 = bdiv * (r12 - 3 * u * r22 - 24 * u * *f2);
  
  exponent_from_coefficients(B2,    bn0, bn1, bn2, A);
}

#include "utils.ih"

void cayley_hamilton_exponent_with_force_terms(gauge_field_t U, gauge_field_t Q, real_field_array_t u, complex_field_array_t v,
                                               complex_field_array_t f0, complex_field_array_t f1, complex_field_array_t f2)
{
  static double const fac_1_3 = 1 / 3.0;
  for (unsigned int x = 0; x < VOLUME; ++x)
    for (unsigned int mu = 0; mu < 4; ++mu)
    {
      double c0 = Q.field[x][mu].c00 * (Q.field[x][mu].c11 * Q.field[x][mu].c22 - Q.field[x][mu].c12 * Q.field[x][mu].c21) + 
                  Q.field[x][mu].c01 * (Q.field[x][mu].c12 * Q.field[x][mu].c20 - Q.field[x][mu].c10 * Q.field[x][mu].c22) +
                  Q.field[x][mu].c02 * (Q.field[x][mu].c10 * Q.field[x][mu].c21 - Q.field[x][mu].c11 * Q.field[x][mu].c20)  ;
      double sign_c0 = copysign(1.0, c0);
      c0 = fabs(c0);

      /* We'll need U^2 down the line, so we might as well calculate it now. */
      /* This is also why we need the buffer space -- can't do matrix multiplication in place */
      _su3_times_su3(*U, *Q, *Q);
      double c1 = 0.5 * (U.field[x][mu].c00 + U.field[x][mu].c11 + U.field[x][mu].c22);

      double c0max = 2.0 * pow(fac_1_3 * c1, 1.5);
      double theta_3 = fac_1_3 * acos(c0 / c0max); 

      u.field_array[mu].field[x] = sqrt(fac_1_3 * c1) * cos(theta_3);
      w.field_array[mu].field[x] = sqrt(c1) * sin(theta_3);
      
      /* Modification w.r.t. Peardon & Morningstar:  w is always positive, so |w| =  w */
      double xi0 = (w.field_array[mu].field[x] > 0.05) ? (sin(w.field_array[mu].field[x]) / w.field_array[mu].field[x]) : 1 - 0.16666666666666667 *  w.field_array[mu].field[x] *  w.field_array[mu].field[x] * (1 - 0.05 *  w.field_array[mu].field[x] *  w.field_array[mu].field[x] * (1 - 0.023809523809523808 *  w.field_array[mu].field[x] * w.field_array[mu].field[x]));
      double divisor = 1.0 / (9.0 * u.field_array[mu].field[x] * u.field_array[mu].field[x] -  w.field_array[mu].field[x] * w.field_array[mu].field[x]);

      /* We can fold in the sign immediately -- c.f. f_j(-c0, c1) = -1^j * conj(f_j(c0, c1)) */
      f0.field_array[mu].field[x] =           divisor * ((u.field_array[mu].field[x] * u.field_array[mu].field[x] -  w.field_array[mu].field[x] * w.field_array[mu].field[x]) * cexp(sign_c0 * 2 * I * u.field_array[mu].field[x]) + 
                                                         cexp(-sign_c0 * I * u.field_array[mu].field[x]) * (8 * u.field_array[mu].field[x] * u.field_array[mu].field[x] * cos(w.field_array[mu].field[x]) + 2 * sign_c0 * I * u.field_array[mu].field[x] * (3 * u.field_array[mu].field[x] * u.field_array[mu].field[x] +  w.field_array[mu].field[x] * w.field_array[mu].field[x]) * xi0));
      f1.field_array[mu].field[x] = sign_c0 * divisor * (2 * u.field_array[mu].field[x] * cexp(sign_c0 * 2 * I * u.field_array[mu].field[x]) - 
                                                         cexp(-sign_c0 * I * u.field_array[mu].field[x]) * (2 * u.field_array[mu].field[x] * cos(w.field_array[mu].field[x]) - sign_c0 * I * (3 * u.field_array[mu].field[x] * u.field_array[mu].field[x] -  w.field_array[mu].field[x] * w.field_array[mu].field[x]) * xi0));
      f2.field_array[mu].field[x] =           divisor * (cexp(2 * sign_c0 * I * u.field_array[mu].field[x]) - cexp(-sign_c0 * I * u.field_array[mu].field[x]) * (cos(w.field_array[mu].field[x]) + 3 * sign_c0 * I * u.field_array[mu].field[x] * xi0));

      U.field[x][mu].c00 = f0.field_array[mu].field[x] + f1.field_array[mu].field[x] * Q.field[x][mu].c00 + f2.field_array[mu].field[x] * U.field[x][mu].c00;
      U.field[x][mu].c01 =                               f1.field_array[mu].field[x] * Q.field[x][mu].c01 + f2.field_array[mu].field[x] * U.field[x][mu].c01;
      U.field[x][mu].c02 =                               f1.field_array[mu].field[x] * Q.field[x][mu].c02 + f2.field_array[mu].field[x] * U.field[x][mu].c02;
      U.field[x][mu].c10 =                               f1.field_array[mu].field[x] * Q.field[x][mu].c10 + f2.field_array[mu].field[x] * U.field[x][mu].c10;
      U.field[x][mu].c11 = f0.field_array[mu].field[x] + f1.field_array[mu].field[x] * Q.field[x][mu].c11 + f2.field_array[mu].field[x] * U.field[x][mu].c11;
      U.field[x][mu].c12 =                               f1.field_array[mu].field[x] * Q.field[x][mu].c12 + f2.field_array[mu].field[x] * U.field[x][mu].c12;
      U.field[x][mu].c20 =                               f1.field_array[mu].field[x] * Q.field[x][mu].c20 + f2.field_array[mu].field[x] * U.field[x][mu].c20;
      U.field[x][mu].c21 =                               f1.field_array[mu].field[x] * Q.field[x][mu].c21 + f2.field_array[mu].field[x] * U.field[x][mu].c21;
      U.field[x][mu].c22 = f0.field_array[mu].field[x] + f1.field_array[mu].field[x] * Q.field[x][mu].c22 + f2.field_array[mu].field[x] * U.field[x][mu].c22;
    }
}

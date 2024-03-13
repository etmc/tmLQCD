/***********************************************************************
 *
 * Copyright (C) 2013 Albert Deuzeman 
 *
 * This file is part of tmLQCD.
 *
 * tmLQCD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * tmLQCD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with tmLQCD.  If not, see <http://www.gnu.org/licenses/>.
 *
 ************************************************************************/

#if HAVE_CONFIG_H
#include <tmlqcd_config.h>
#endif
#include <math.h>
#include <complex.h> 

#if (defined SSE || defined SSE2 || defined SSE3)
# include "sse.h"
#endif
#include "su3.h"

#ifndef TM_USE_OMP
static
#endif
void exponent_from_coefficients(su3 *out, _Complex double f0, _Complex double f1, _Complex double f2, su3 const *in)                                  
{
  su3 ALIGN tmp;
  _complex_times_su3(tmp, f2, *in);
  _su3_add_equals_complex_identity(tmp, f1);
  _su3_times_su3(*out, tmp, *in);
  _su3_add_equals_complex_identity(*out, f0);
}

void cayley_hamilton_exponent(su3* expA, su3 const *A)
{
  static double const fac_1_3 = 1 / 3.0;
    
  _Complex double f0,f1,f2;

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
    f1 = I;
    f2 = -0.5;
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

  f0 = divisor * (ma * (u * u -  w * w) + mb * (8 * u * u * cw + 2 * I * u * (3 * u * u +  w * w) * xi0));
  f1 = divisor * (-2 * I * u * ma + mb * (2 * I * u * cw + (3 * u * u -  w * w) * xi0));
  f2 = divisor * (mb * (cw + 3 * I * u * xi0) - ma);

  /* The first point where we use the symmetry relations to calculate the negative c0 possibility */
  if (c0_negative)
  {
    f0 = conj(f0);
    f1 = conj(f1);
    f2 = conj(f2);
  }
  
  exponent_from_coefficients(expA, f0, f1, f2, A);
   
  return;
 }

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


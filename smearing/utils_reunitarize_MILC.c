#include "utils.ih"

/* No reunitarization code seems to be available, so I've adapted (stolen) this routine from the MILC code (who stole it elsewhere, I think ;]) -- AD. */
void reunitarize(su3 *omega)
{
  _Complex double a, bj0, bj1, bj2, t;

  /* first normalize row 0 */
  a = 1.0 / sqrt(conj(omega->c00) * omega->c00 + conj(omega->c01) * omega->c01 +conj(omega->c02) * omega->c02);

  omega->c00 *= a;
  omega->c01 *= a;
  omega->c02 *= a;

  /* now make row 1 orthogonal to row 0 */
  a = conj(omega->c00) * omega->c10 + conj(omega->c01) * omega->c11 + conj(omega->c02) * omega->c12;

  /* row 1 -= a * row 0 */
  omega->c10 -= a * omega->c00;
  omega->c11 -= a * omega->c01;
  omega->c12 -= a * omega->c02;

  /* now normalize row 1 */
  a = 1.0 / sqrt(conj(omega->c10) * omega->c10 + conj(omega->c11) * omega->c11 +conj(omega->c12) * omega->c12);

  omega->c10 *= a;
  omega->c11 *= a;
  omega->c12 *= a;

  /* reconstruct row 2 */
  bj0 = omega->c00;
  bj1 = omega->c01;
  bj2 = omega->c02;

  omega->c20  = bj1 * omega->c12;
  omega->c20 -= bj2 * omega->c11

  omega->c21  = bj2 * omega->c10;
  omega->c21 -= bj0 * omega->c12;

  omega->c22  = bj0 * omega->c11;
  omega->c22 -= bj1r * omega->c10;
}

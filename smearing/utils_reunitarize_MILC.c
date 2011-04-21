#include "utils.ih"

/* No reunitarization code seems to be available, so I've adapted (stolen) this routine from the MILC code (who stole it elsewhere, I think ;]) -- AD. */
void reunitarize(su3 *omega)
{
  register double ar, ai, bj0r, bj0i, bj1r, bj1i, bj2r, bj2i, tr, ti;

  /* first normalize row 0 */
  ar = omega->c00.re * omega->c00.re +    /* sum of squares of row */
	omega->c00.im * omega->c00.im +
	omega->c01.re * omega->c01.re +
	omega->c01.im * omega->c01.im +
	omega->c02.re * omega->c02.re +
	omega->c02.im * omega->c02.im;

  ar = 1.0 / sqrt(ar);	           /* used to normalize row */
  omega->c00.re *= ar;
  omega->c00.im *= ar;
  omega->c01.re *= ar;
  omega->c01.im *= ar;
  omega->c02.re *= ar;
  omega->c02.im *= ar;

  /* now make row 1 orthogonal to row 0 */
  ar = omega->c00.re * omega->c10.re +     /* real part of 0 dot 1 */
	omega->c00.im * omega->c10.im +
	omega->c01.re * omega->c11.re +
	omega->c01.im * omega->c11.im +
	omega->c02.re * omega->c12.re +
	omega->c02.im * omega->c12.im;
  ai = omega->c00.re * omega->c10.im -     /* imag part of 0 dot 1 */
	omega->c00.im * omega->c10.re +
	omega->c01.re * omega->c11.im -
	omega->c01.im * omega->c11.re +
	omega->c02.re * omega->c12.im -
	omega->c02.im * omega->c12.re;

  /* row 1 -= a * row 0 */
  omega->c10.re -= ar * omega->c00.re - ai * omega->c00.im;
  omega->c10.im -= ar * omega->c00.im + ai * omega->c00.re;
  omega->c11.re -= ar * omega->c01.re - ai * omega->c01.im;
  omega->c11.im -= ar * omega->c01.im + ai * omega->c01.re;
  omega->c12.re -= ar * omega->c02.re - ai * omega->c02.im;
  omega->c12.im -= ar * omega->c02.im + ai * omega->c02.re;

  /* now normalize row 1 */
  ar = omega->c10.re * omega->c10.re +    /* sum of squares of row */
	omega->c10.im * omega->c10.im +
	omega->c11.re * omega->c11.re +
	omega->c11.im * omega->c11.im +
	omega->c12.re * omega->c12.re +
	omega->c12.im * omega->c12.im;

  ar = 1.0 / sqrt(ar);	           /* used to normalize row */
  omega->c10.re *= ar;
  omega->c10.im *= ar;
  omega->c11.re *= ar;
  omega->c11.im *= ar;
  omega->c12.re *= ar;
  omega->c12.im *= ar;

  /* reconstruct row 2 */
  bj0r = omega->c00.re;
  bj0i = omega->c00.im;
  bj1r = omega->c01.re;
  bj1i = omega->c01.im;
  bj2r = omega->c02.re;
  bj2i = omega->c02.im;

  ar = omega->c12.re;
  ai = omega->c12.im;
  tr = bj1r * ar - bj1i * ai;
  ti = bj1r * ai + bj1i * ar;
  ar = omega->c11.re;
  ai = omega->c11.im;
  tr = tr - bj2r * ar + bj2i * ai;
  ti = ti - bj2r * ai - bj2i * ar;
  omega->c20.re = tr;
  omega->c20.im = -ti;

  ar = omega->c10.re;
  ai = omega->c10.im;
  tr = bj2r * ar - bj2i * ai;
  ti = bj2r * ai + bj2i * ar;
  ar = omega->c12.re;
  ai = omega->c12.im;
  tr = tr - bj0r * ar + bj0i * ai;
  ti = ti - bj0r * ai - bj0i * ar;
  omega->c21.re = tr;
  omega->c21.im = -ti;

  ar = omega->c11.re;
  ai = omega->c11.im;
  tr = bj0r * ar - bj0i * ai;
  ti = bj0r * ai + bj0i * ar;
  ar = omega->c10.re;
  ai = omega->c10.im;
  tr = tr - bj1r * ar + bj1i * ai;
  ti = ti - bj1r * ai - bj1i * ar;
  omega->c22.re = tr;
  omega->c22.im = -ti;
}

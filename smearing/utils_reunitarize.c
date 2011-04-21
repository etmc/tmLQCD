#include "utils.ih"

/* Method based on Givens' rotations, as used by Urs Wenger */
void reunitarize(su3 *omega)
{
  static su3 w, rot, tmp;
  static double trace_old, trace_new;
  static complex s0, s1;
  static double scale;

  _su3_one(w);
  trace_old = omega->c00.re + omega->c11.re + omega->c22.re;

  for (int iter = 0; iter < 200; ++iter)
  {
    /* Givens' rotation 01 */

      _assign_add_conj(s0, omega->c00, omega->c11);
      _assign_diff_conj(s1, omega->c01, omega->c10);
      scale = 1.0 / sqrt(_complex_square_norm(s0) + _complex_square_norm(s1));
      _mult_assign_real(s0, scale);
      _mult_assign_real(s1, scale);

      /* Projecting */
      _su3_one(rot);
      _assign(rot.c00, s0);
      _complex_conj(rot.c11, s0);
      _assign(rot.c01, s1);
      _complex_conj_chgsig(rot.c10, s1);

      _su3_times_su3(tmp, rot, w);
      _su3_assign(w, tmp);
      _su3_times_su3d(tmp, *omega, rot);
      _su3_assign(*omega, tmp);

    /* Givens' rotation 12 */

      _assign_add_conj(s0, omega->c11, omega->c22);
      _assign_diff_conj(s1, omega->c12, omega->c21);
      scale = 1.0 / sqrt(_complex_square_norm(s0) + _complex_square_norm(s1));
      _mult_assign_real(s0, scale);
      _mult_assign_real(s1, scale);

      /* Projecting */
      _su3_one(rot);
      _assign(rot.c11, s0);
      _complex_conj(rot.c22, s0);
      _assign(rot.c12, s1);
      _complex_conj_chgsig(rot.c21, s1);

      _su3_times_su3(tmp, rot, w);
      _su3_assign(w, tmp);
      _su3_times_su3d(tmp, *omega, rot);
      _su3_assign(*omega, tmp);


    /* Givens' rotation 20 */

      _assign_add_conj(s0, omega->c22, omega->c00);
      _assign_diff_conj(s1, omega->c20, omega->c02);
      scale = 1.0 / sqrt(_complex_square_norm(s0) + _complex_square_norm(s1));
      _mult_assign_real(s0, scale);
      _mult_assign_real(s1, scale);

      /* Projecting */
      _su3_one(rot);
      _assign(rot.c22, s0);
      _complex_conj(rot.c00, s0);
      _assign(rot.c20, s1);
      _complex_conj_chgsig(rot.c02, s1);

      _su3_times_su3(tmp, rot, w);
      _su3_assign(w, tmp);
      _su3_times_su3d(tmp, *omega, rot);
      _su3_assign(*omega, tmp);

    trace_new = omega->c00.re + omega->c11.re + omega->c22.re;

    if (trace_new - trace_old < 1e-15)
      break;
    trace_old = trace_new;
  }
  _su3_assign(*omega, w);
}

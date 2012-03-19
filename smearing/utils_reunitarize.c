#include "utils.ih"

/* Method based on Givens' rotations, as used by Urs Wenger */
void reunitarize(su3 *omega)
{
  static su3 w, rot, tmp;
  static double trace_old, trace_new;
  static _Complex double s0, s1;
  static double scale;

  _su3_one(w);
  trace_old = omega->c00 + omega->c11 + omega->c22;

  for (int iter = 0; iter < 200; ++iter)
  {
    /* Givens' rotation 01 */
      s0 = omega->c00 + conj(omega->c11);
      s1 = omega->c01 - conj(omega->c10);
      scale = 1.0 / sqrt(conj(s0) * s0 + conj(s1) * s1);
      s0 *= scale;
      s1 *= scale;
      
      /* Projecting */
      _su3_one(rot);
      rot.c00 = s0;
      rot.c11 = conj(s0);
      rot.c01 = s1;
      rot.c10 = -conj(s1);

      _su3_times_su3(tmp, rot, w);
      _su3_assign(w, tmp);
      _su3_times_su3d(tmp, *omega, rot);
      _su3_assign(*omega, tmp);

    /* Givens' rotation 12 */
      s0 = omega->c11 + conj(omega->c22);
      s1 = omega->c12 - conj(omega->c21);
      scale = 1.0 / sqrt(conj(s0) * s0 + conj(s1) * s1);
      s0 *= scale;
      s1 *= scale;

      /* Projecting */
      _su3_one(rot);
      rot.c11 = s0;
      rot.c22 = conj(s0);
      rot.c12 = s1;
      rot.c21 = -conj(s1);

      _su3_times_su3(tmp, rot, w);
      _su3_assign(w, tmp);
      _su3_times_su3d(tmp, *omega, rot);
      _su3_assign(*omega, tmp);

    /* Givens' rotation 20 */
      s0 = omega->c22 + conj(omega->c00);
      s1 = omega->c20 - conj(omega->c02);
      scale = 1.0 / sqrt(conj(s0) * s0 + conj(s1) * s1);
      s0 *= scale;
      s1 *= scale;

      /* Projecting */
      _su3_one(rot);
      rot.c22 = s0;
      rot.c00 = conj(s0);
      rot.c20 = s1;
      rot.c02 = -conj(s1);

      _su3_times_su3(tmp, rot, w);
      _su3_assign(w, tmp);
      _su3_times_su3d(tmp, *omega, rot);
      _su3_assign(*omega, tmp);

    trace_new = omega->c00 + omega->c11 + omega->c22;

    if (trace_new - trace_old < 1e-15)
      break;
    trace_old = trace_new;
  }
  _su3_assign(*omega, w);
}

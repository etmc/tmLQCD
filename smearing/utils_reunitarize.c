#include "utils.ih"

/* Method based on Givens' rotations, as used by Urs Wenger */
void reunitarize(su3 *in)
{
  static su3 w, rot, tmp;
  static double trace_old, trace_new;
  static _Complex double s0, s1;
  static double scale;

  _su3_one(w);
  trace_old = in->c00 + in->c11 + in->c22;

  for (int iter = 0; iter < 200; ++iter)
  {
    /* Givens' rotation 01 */
      s0 = in->c00 + conj(in->c11);
      s1 = in->c01 - conj(in->c10);
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
      _su3_times_su3d(tmp, *in, rot);
      _su3_assign(*in, tmp);

    /* Givens' rotation 12 */
      s0 = in->c11 + conj(in->c22);
      s1 = in->c12 - conj(in->c21);
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
      _su3_times_su3d(tmp, *in, rot);
      _su3_assign(*in, tmp);

    /* Givens' rotation 20 */
      s0 = in->c22 + conj(in->c00);
      s1 = in->c20 - conj(in->c02);
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
      _su3_times_su3d(tmp, *in, rot);
      _su3_assign(*in, tmp);

    trace_new = in->c00 + in->c11 + in->c22;

    if (trace_new - trace_old < 1e-15)
      break;
    trace_old = trace_new;
  }
  _su3_assign(*in, w);
}

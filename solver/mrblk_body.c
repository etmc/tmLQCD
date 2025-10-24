/***********************************************************************
 * Copyright (C) 2016 Carsten Urbach
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
 ***********************************************************************/

int _PSWITCH(mrblk)(_PTSWITCH(spinor) * const P, _PTSWITCH(spinor) * const Q,
                    _PTSWITCH(spinor) * const s_, const int max_iter, const double eps_sq,
                    const int rel_prec, const int N, _PTSWITCH(matrix_mult_blk) f, const int blk) {
  int i = 0;
  _F_TYPE norm_r, beta;
  _C_TYPE alpha;
  _PTSWITCH(spinor) * r;
  _PTSWITCH(spinor) * s[3];
  const int parallel = 0;

  s[0] = s_;
  s[1] = s[0] + N + 1;
  s[2] = s[1] + N + 1;

  r = s[0];
  norm_r = _PSWITCH(square_norm_ts)(Q, N, parallel);

  _PSWITCH(zero_spinor_field)(P, N);
  f(s[2], P, blk);
  _PSWITCH(diff_ts)(r, Q, s[2], N);
  norm_r = _PSWITCH(square_norm_ts)(r, N, parallel);
  if (g_proc_id == g_stdio_proc && g_debug_level > 2 && blk == 0) {
    printf("MRblk iteration= %d  |res|^2= %e\n", i, norm_r);
    fflush(stdout);
  }

  while ((norm_r > eps_sq) && (i < max_iter)) {
    i++;
    f(s[1], r, blk);
    alpha = _PSWITCH(scalar_prod_ts)(s[1], r, N, parallel);
    beta = _PSWITCH(square_norm_ts)(s[1], N, parallel);
    alpha /= beta;
    _PSWITCH(assign_add_mul_ts)(P, r, alpha, N);
    if (i % 50 == 0) {
      f(s[2], P, blk);
    } else {
      _PSWITCH(assign_add_mul_ts)(s[2], s[1], alpha, N);
    }
    _PSWITCH(diff_ts)(r, Q, s[2], N);
    norm_r = _PSWITCH(square_norm_ts)(r, N, parallel);
    if (g_proc_id == g_stdio_proc && g_debug_level > 2 && blk == 0) {
      printf("MRblk iteration= %d  |res|^2= %g\n", i, norm_r);
      fflush(stdout);
    }
  }
  if (norm_r > eps_sq) {
    return (-1);
  }
  return (i);
}

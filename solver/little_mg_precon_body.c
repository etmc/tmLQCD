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

void _PSWITCH(little_mg_precon)(_Complex _F_TYPE * const out, _Complex _F_TYPE * const in) {
  // phi = PD_c^{-1} P^dagger in
  _PSWITCH(little_project_eo)(out, in, g_N_s);
  // in - D*phi
  _PSWITCH(little_D_sym)((_Complex _F_TYPE *) work[2], out);
  _PSWITCH(ldiff)((_Complex _F_TYPE *) work[3], in, (_Complex _F_TYPE *) work[2], nb_blocks*g_N_s);
  // sum with phi
  _PSWITCH(ladd)(out, (_Complex _F_TYPE *) work[3], out, nb_blocks*g_N_s);
  return;
}

/***********************************************************************
 *
 * Copyright (C) 2008 Albert Deuzeman, Siebren Reker, Carsten Urbach
 *                    Claude Tadonki
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

void _PSWITCH(little_project_eo)(_Complex _F_TYPE * const out, _Complex _F_TYPE * const in, const int  N) {
  static _Complex _F_TYPE * phi;
  static _Complex _F_TYPE * psi;

  if(init_dfl_projector == 0) {
    alloc_dfl_projector();
  }

  phi = (_Complex _F_TYPE *)work[2];
  psi = (_Complex _F_TYPE *)work[3];

  /* NOTE IS THIS REALLY NECESSARY/CORRECT? */
  for(int i = 0; i < N; i++) {
    phi[i] = _PSWITCH(lscalar_prod)(_PSWITCH(little_dfl_fields_eo)[i], in, nb_blocks*N, 0);
  }

#ifdef MPI
  MPI_Allreduce(phi, psi, N, _MPI_C_TYPE, MPI_SUM, MPI_COMM_WORLD);
#else
  memcpy(psi, phi, N*sizeof(_Complex _F_TYPE));
#endif

  /* apply inverse of little_A_eo */
  for(int i = 0; i < N; i++) {
    (phi[i]) = 0.0;
    for(int j = 0; j < N; j++) {
      (phi[i]) += (_PSWITCH(little_A_eo)[j*N + i]) * (psi[j]);
    }
  }
  _PSWITCH(lmul)(out, phi[0], _PSWITCH(little_dfl_fields_eo)[0], nb_blocks*N);
  for(int i = 1; i < N; i++) {
    _PSWITCH(lassign_add_mul)(out, _PSWITCH(little_dfl_fields_eo)[i], phi[i], nb_blocks*N);
  }

  return;
}

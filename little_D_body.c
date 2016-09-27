/***********************************************************************
 *
 * Copyright (C) 2008 Albert Deuzeman, Siebren Reker, Carsten Urbach
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

void _PSWITCH(little_D)(_C_TYPE * v, _C_TYPE *w) {
  int sq = g_N_s*g_N_s;
  _PSWITCH(CONE) = 1.0;
  _PSWITCH(CMONE) = -1.0;
  _PSWITCH(CZERO) = 0.0;

  if(dfl_subspace_updated) {
    compute_little_D(0);
    dfl_subspace_updated = 0;
  }
  
  _PSWITCH(little_field_gather)(w);
  
#ifdef TM_USE_OMP
#pragma omp parallel for
#endif
  for(int i = 0; i < nb_blocks; i++) {
    for(int j = 0; j < 9; j++) {
      _MV(zgemv)("N", &g_N_s, &g_N_s, &_PSWITCH(CONE), _PSWITCH(block_list[i].little_dirac_operator) + j * sq,
                 &g_N_s, w + (nb_blocks * j + i) * g_N_s, &ONE, &_PSWITCH(CONE), v + i * g_N_s, &ONE, 1);
    }
  }
  return;
}


void _PSWITCH(little_Q_pm)(_C_TYPE * v, _C_TYPE *w) {
  _C_TYPE * tmp = (_C_TYPE*) aligned_malloc_zero(nb_blocks * 9 * g_N_s * sizeof(_C_TYPE));
  double musave= g_mu;
  if(dfl_subspace_updated) {
    g_mu = 0.;
    compute_little_D(1);
    g_mu = musave;
    dfl_subspace_updated = 0;
  }
  _PSWITCH(little_D)(tmp, w);
  _PSWITCH(little_D)(v, tmp);
  aligned_free(tmp);
  _PSWITCH(lassign_add_mul)(v, w, g_mu*g_mu + g_mu2*g_mu2, nb_blocks*g_N_s);
}


void _PSWITCH(little_D_sym)(_C_TYPE * v, _C_TYPE *w) {
  
  _C_TYPE* tmpc1, * tmpc2, * tmpc3;
  tmpc1 = (_C_TYPE*) aligned_malloc_zero(3*nb_blocks * 9 * g_N_s * sizeof(_C_TYPE));
  tmpc2 = tmpc1 + nb_blocks * 9 * g_N_s;
  tmpc3 = tmpc1 + 2*nb_blocks * 9 * g_N_s;
  
  if(dfl_subspace_updated) {
    compute_little_D(0);
    dfl_subspace_updated = 0;
  }
  
  _PSWITCH(little_D_hop)(0, tmpc1, w);
  _PSWITCH(little_D_ee_inv)(tmpc2, tmpc1);
  _PSWITCH(little_D_hop)(1, tmpc3, tmpc2);
  _PSWITCH(little_Dhat_lhs)(v, w, tmpc3);
  
  aligned_free(tmpc1);
  return;
}


void _PSWITCH(little_D_ee_inv)(_C_TYPE * v, _C_TYPE *w) {
  int i;
  _PSWITCH(CONE) = 1.0;
  _PSWITCH(CMONE) = -1.0;
  _PSWITCH(CZERO) = 0.0;
  
  for(i = 0; i < nb_blocks/2; i++) {
    _MV(zgemv)("N", &g_N_s, &g_N_s, &_PSWITCH(CONE), _PSWITCH(block_list[i].little_dirac_operator_eo),
               &g_N_s, w + i * g_N_s, &ONE, &_PSWITCH(CZERO), v + i * g_N_s, &ONE, 1);
  }
  return;
}


void _PSWITCH(little_D_hop)(int eo, _C_TYPE * v, _C_TYPE *w) {
  int i, j, i_eo,sq = g_N_s*g_N_s;
  _PSWITCH(CONE) = 1.0;
  _PSWITCH(CMONE) = -1.0;
  _PSWITCH(CZERO) = 0.0;

  i_eo = (eo+1) % 2;
  
  _PSWITCH(little_field_gather_eo)(eo, w + i_eo*nb_blocks*g_N_s/2);
  
  for(j = 1; j < 9; j++) {
    for(i = 0; i < nb_blocks/2; i++) {
      _MV(zgemv)("N", &g_N_s, &g_N_s, &_PSWITCH(CONE), _PSWITCH(block_list[eo*(nb_blocks/2)+i].little_dirac_operator_eo) + j * sq,
                 &g_N_s, w + (nb_blocks * j + (nb_blocks/2)*i_eo+i) * g_N_s, &ONE, &_PSWITCH(CONE), v + (eo*nb_blocks/2+i) * g_N_s, &ONE, 1);
    } 
  }
  return;
}

void _PSWITCH(little_Dhat_lhs)(_C_TYPE * v, _C_TYPE *w, _C_TYPE *u) {
  int i,j;
  _PSWITCH(CONE) = 1.0;
  _PSWITCH(CMONE) = -1.0;
  _PSWITCH(CZERO) = 0.0;


  for(i = nb_blocks/2; i < nb_blocks; i++) {
    _MV(zgemv)("N", &g_N_s, &g_N_s, &_PSWITCH(CONE), _PSWITCH(block_list[i].little_dirac_operator_eo),
               &g_N_s, w + i * g_N_s, &ONE, &_PSWITCH(CZERO), v + i * g_N_s, &ONE, 1);
  }
  
  for (i=nb_blocks/2; i < nb_blocks; i++) {
    for (j=0;j<g_N_s;j++) {
      *(v+ i * g_N_s+ j) = *(v+ i * g_N_s+ j) - *(u+ i * g_N_s+ j);
    }
  }
  return;
}



void _PSWITCH(little_Dhat_rhs)(int eo, _C_TYPE * v, double r, _C_TYPE *w) {
  int i, j;
  
  for(i = 0; i < nb_blocks/2; i++) {
    for (j=0;j<g_N_s;j++) {
      *(v+eo*nb_blocks*g_N_s/2+i*g_N_s+j) = *(w+eo*nb_blocks*g_N_s/2+i*g_N_s+j) + r * *(v+eo*nb_blocks*g_N_s/2+i*g_N_s+j);
    }
  }
  return;
}

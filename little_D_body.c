
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
  
  /* all the mpilocal stuff first */
#ifdef OMP
#pragma omp parallel for
#endif
  for(int i = 0; i < nb_blocks; i++) {
    /* diagonal term */
    _MV(dummy)("N", &g_N_s, &g_N_s, &_PSWITCH(CONE), _PSWITCH(block_list[i].little_dirac_operator),
               &g_N_s, w + i * g_N_s, &ONE, &_PSWITCH(CZERO), v + i * g_N_s, &ONE, 1);
  }
  /* offdiagonal terms */
#ifdef OMP
#pragma omp parallel for
#endif
  for(int j = 1; j < 9; j++) {
#ifdef OMP
#pragma omp parallel for
#endif
    for(int i = 0; i < nb_blocks; i++) {
      _MV(dummy)("N", &g_N_s, &g_N_s, &_PSWITCH(CONE), _PSWITCH(block_list[i].little_dirac_operator) + j * sq,
                 &g_N_s, w + (nb_blocks * j + i) * g_N_s, &ONE, &_PSWITCH(CONE), v + i * g_N_s, &ONE, 1);
    }
  }
  return;
}


void _PSWITCH(little_Q_pm)(_C_TYPE * v, _C_TYPE *w) {
  _C_TYPE * tmp = calloc(nb_blocks * 9 * g_N_s, sizeof(_C_TYPE));
  double musave= g_mu;
  if(dfl_subspace_updated) {
    g_mu = 0.;
    compute_little_D(1);
    g_mu = musave;
    dfl_subspace_updated = 0;
  }
  _PSWITCH(little_D)(tmp, w);
  _PSWITCH(little_D)(v, tmp);
  free(tmp);
  _PSWITCH(lassign_add_mul)(v, w, g_mu*g_mu + g_mu2*g_mu2, nb_blocks*g_N_s);
  //memcpy(v, w, nb_blocks * g_N_s*sizeof(_C_TYPE));
}


void _PSWITCH(little_D_sym)(_C_TYPE * v, _C_TYPE *w) {
  
  _C_TYPE* tmpc1, * tmpc2, * tmpc3;
  tmpc1 = calloc(3*nb_blocks * 9 * g_N_s, sizeof(_C_TYPE));
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
  
  free(tmpc1);
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


void _PSWITCH(little_D_hop)(int eo,_C_TYPE * v, _C_TYPE *w) {
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

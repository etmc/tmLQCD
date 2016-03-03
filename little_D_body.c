
void _PSWITCH(little_D)(_C_TYPE * v, _C_TYPE *w) {
  int sq = g_N_s*g_N_s;
  _PSWITCH(CONE) = 1.0;
  _PSWITCH(CMONE) = -1.0;
  _PSWITCH(CZERO) = 0.0;

  if(dfl_subspace_updated) {
    compute_little_D(0);
    dfl_subspace_updated = 0;
  }
  
  little_field_gather(w);
  
  /* all the mpilocal stuff first */
  for(int i = 0; i < nb_blocks; i++) {
    /* diagonal term */
    _MV(dummy)("N", &g_N_s, &g_N_s, &CONE, _PSWITCH(block_list[i].little_dirac_operator),
               &g_N_s, w + i * g_N_s, &ONE, &CZERO, v + i * g_N_s, &ONE, 1);
  }
  /* offdiagonal terms */
  for(int j = 1; j < 9; j++) {
    for(int i = 0; i < nb_blocks; i++) {
      _MV(dummy)("N", &g_N_s, &g_N_s, &CONE, _PSWITCH(block_list[i].little_dirac_operator) + j * sq,
                 &g_N_s, w + (nb_blocks * j + i) * g_N_s, &ONE, &CONE, v + i * g_N_s, &ONE, 1);
    }
  }
  return;
}

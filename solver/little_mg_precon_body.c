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

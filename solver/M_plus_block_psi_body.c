void _PSWITCH(Mtm_plus_block_psi)(_PTSWITCH(spinor) * const l, _PTSWITCH(spinor) * const k, const int i) {
  block * blk = &block_list[i];
  int vol = (*blk).volume/2;
  _PSWITCH(Block_H_psi)(blk, &_PTSWITCH(g_spinor_field)[DUM_MATRIX+1][i*vol], k, EO);
  _PSWITCH(mul_one_pm_imu_inv)(&_PTSWITCH(g_spinor_field)[DUM_MATRIX+1][i*vol], +1., vol);
  _PSWITCH(Block_H_psi)(blk, &_PTSWITCH(g_spinor_field)[DUM_MATRIX][i*vol], &_PTSWITCH(g_spinor_field)[DUM_MATRIX+1][i*vol], OE);
  _PSWITCH(mul_one_pm_imu_sub_mul)(l, k, &_PTSWITCH(g_spinor_field)[DUM_MATRIX][i*vol], +1., vol);
  return;
}

void _PSWITCH(Msw_plus_block_psi)(_PTSWITCH(spinor) * l, _PTSWITCH(spinor) *  k, const int i) {
  block * blk = &block_list[i];
  int vol = (*blk).volume/2;
  _PSWITCH(Block_H_psi)(blk, &_PTSWITCH(g_spinor_field)[DUM_MATRIX+1][i*vol], k, EO);
  _PSWITCH(assign_mul_one_sw_pm_imu_inv_block)(EE, &_PTSWITCH(g_spinor_field)[DUM_MATRIX][i*vol],&_PTSWITCH(g_spinor_field)[DUM_MATRIX+1][i*vol], g_mu, blk);
  _PSWITCH(Block_H_psi)(blk, &_PTSWITCH(g_spinor_field)[DUM_MATRIX+1][i*vol], &_PTSWITCH(g_spinor_field)[DUM_MATRIX][i*vol], OE);
  _PSWITCH(assign_mul_one_sw_pm_imu_block)(OO, &_PTSWITCH(g_spinor_field)[DUM_MATRIX][i*vol],k,g_mu,blk);
  _PSWITCH(diff)(l,&_PTSWITCH(g_spinor_field)[DUM_MATRIX][i*vol],&_PTSWITCH(g_spinor_field)[DUM_MATRIX+1][i*vol],vol);
  return;
}

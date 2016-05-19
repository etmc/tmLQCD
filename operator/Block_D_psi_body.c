/* apply the Dirac operator to the block local spinor field s */
/* and store the result in block local spinor field rr        */
/* for block blk                                              */
/* the block local gauge field is assumed to be in the order  */
/* that is needed int local_D, which means also that it is a  */
/* double copy                                                */
// CU: has problems with SSE2,3
void _PSWITCH(Block_D_psi)(block * blk, _PTSWITCH(spinor) * const rr, _PTSWITCH(spinor) * const s) {

  if(g_c_sw > 0)
    _PSWITCH(Block_Dsw_psi)(blk, rr, s);
  else 
    _PSWITCH(Block_Dtm_psi)(blk,rr,s);

  return;
}
 
//this version is for c_sw=0
void _PSWITCH(Block_Dtm_psi)(block * blk, _PTSWITCH(spinor) * const rr, _PTSWITCH(spinor) * const s) {
  int i;
  _PTSWITCH(spinor) *r = rr;
  _PTSWITCH(spinor) *t = s;
  _PSWITCH(su3) * u = blk->_PSWITCH(u);
  int * idx = blk->idx;
  static _C_TYPE rhoa, rhob;
  _PTSWITCH(spinor) ALIGN tmpr;
  if(blk_gauge_eo) {
    _PSWITCH(init_blocks_gaugefield)();
  }
  rhoa = 1.0 + (_F_TYPE)g_mu * I;
  rhob = conj(rhoa);

  /* set the boundary term to zero */
  _spinor_null(rr[blk->volume]);
  _spinor_null(s[blk->volume]);

  for(i = 0; i < blk->volume; i++) {
    _complex_times_vector(tmpr.s0, rhoa, t->s0);
    _complex_times_vector(tmpr.s1, rhoa, t->s1);
    _complex_times_vector(tmpr.s2, rhob, t->s2);
    _complex_times_vector(tmpr.s3, rhob, t->s3);

    _PSWITCH(local_H)(r, s, u, idx, &tmpr);

    r++;
    t++;
    idx += 8;
    u += 8;
  }

  return;
}


/* apply the Dirac operator to the block local spinor field s */
/* and store the result in block local spinor field rr        */
/* for block blk                                              */
/* the block local gauge field is assumed to be in the order  */
/* that is needed int local_D, which means also that it is a  */
/* double copy                                                */
// CU: has problems with SSE2,3
void _PSWITCH(Block_Dsw_psi)(block * blk, _PTSWITCH(spinor) * const rr, _PTSWITCH(spinor) * const s) {
  _PTSWITCH(spinor) *r = rr;
  _PTSWITCH(spinor) *t = s;
  _PSWITCH(su3) * u = blk->_PSWITCH(u);
  int * idx = blk->idx;
  //static _Complex double rhoa, rhob;
  _PTSWITCH(spinor) ALIGN tmpr;

  int it, ix, iy, iz; //lexiographic index of the site w.r.t the block
  int bt, bx, by, bz; //block coordinate on the local mpi process
  int dT, dX, dY, dZ; //block size
  int sT, sX, sY, sZ; //constant shifts
  int lx; //lexiographic index of the block site w.r.t the local mpi process

  dT = blk->BT;
  dX = blk->BLX;
  dY = blk->BLY;
  dZ = blk->BLZ;

  bt = blk->mpilocal_coordinate[0];
  bx = blk->mpilocal_coordinate[1];
  by = blk->mpilocal_coordinate[2];
  bz = blk->mpilocal_coordinate[3];

  sT = bt * dT;
  sX = bx * dX;
  sY = by * dY;
  sZ = bz * dZ;
  
  if(blk_gauge_eo) {
    init_blocks_gaugefield();
  }
  
  /* set the boundary term to zero */
  _spinor_null(rr[blk->volume]);
  _spinor_null(s[blk->volume]);

  for(int i = 0; i < blk->volume; i++) {

    iz = i % dZ;
    iy = (i / dZ) % dY;
    ix = (i / (dZ * dY)) % dX;
    it = i / (dZ * dY * dX);

    lx = g_ipt[it + sT][ix + sX][iy + sY][iz + sZ];

    _PSWITCH(assign_mul_one_sw_pm_imu_site_lexic)(lx, &tmpr, t, (_F_TYPE)g_mu);

    _PSWITCH(local_H)(r, s, u, idx, &tmpr);

    r++;
    t++;
    idx += 8;
    u += 8;
  }

  return;
}


/* Apply Hopping Matrix to a even(odd) spinor */
void _PSWITCH(Block_H_psi)(block * blk, _PTSWITCH(spinor) * const rr, _PTSWITCH(spinor) * const s, 
			   const int eo) {
  int i;
  _PTSWITCH(spinor) *r = rr;
  _PSWITCH(su3) * u = blk->_PSWITCH(u);
  int * eoidx = blk->evenidx;
  _PTSWITCH(spinor) ALIGN tmpr;

  if(!blk_gauge_eo) {
    _PSWITCH(init_blocks_eo_gaugefield)();
  }

  /* for OE */
  if(eo == 1) {
    u = blk->_PSWITCH(u) + blk->volume*8/2;
    eoidx = blk->oddidx;
  }

  /* set the boundary term to zero */
  _spinor_null(rr[blk->volume/2]);
  _spinor_null(s[blk->volume/2]);
  
  for(i = 0; i < blk->volume/2; i++) {
    _spinor_null(tmpr);

    _PSWITCH(local_H)(r, s, u, eoidx, &tmpr);

    r++;
    eoidx += 8;
    u += 8;
  }

  return;
}

void _PSWITCH(assign_mul_one_sw_pm_imu)(const int ieo, 
					_PTSWITCH(spinor) * const k, _PTSWITCH(spinor) * const l,
					const _F_TYPE mu) {
#ifdef OMP
#pragma omp parallel
  {
#endif
    _PTSWITCH(su3_vector) ALIGN chi, psi1, psi2;
    int ix;
    int ioff;
    const _PSWITCH(su3) *w1, *w2, *w3;
    _PTSWITCH(spinor) *r;
    const _PTSWITCH(spinor) *s;
  
    if(ieo == 0) {
      ioff = 0;
    } 
    else {
      ioff = (VOLUME+RAND)/2;
    }
    /************************ loop over all lattice sites *************************/
#ifdef OMP
#pragma omp for
#endif
    for(unsigned icx = ioff; icx < (VOLUME/2+ioff); icx++) {
      ix = g_eo2lexic[icx];
    
      r = k + icx-ioff;
      s = l + icx-ioff;

      // upper two spin components first
      w1=&_PSWITCH(sw)[ix][0][0];
      w2=w1+2; /*&sw[ix][1][0];*/
      w3=w1+4; /*&sw[ix][2][0];*/
      _su3_multiply(psi1,*w1,(*s).s0); 
      _su3_multiply(chi,*w2,(*s).s1);
      _vector_add_assign(psi1,chi);
      _su3_inverse_multiply(psi2,*w2,(*s).s0); 
      _su3_multiply(chi,*w3,(*s).s1);
      _vector_add_assign(psi2,chi); 

      // add in the twisted mass term (plus in the upper components)
      _vector_add_i_mul(psi1, mu, (*s).s0);
      _vector_add_i_mul(psi2, mu, (*s).s1);

      _vector_assign((*r).s0, psi1);
      _vector_assign((*r).s1, psi2);

      // now lower to spin components
      w1++; /*=&sw[ix][0][1];*/
      w2++; /*=&sw[ix][1][1];*/
      w3++; /*=&sw[ix][2][1];*/
      _su3_multiply(psi1,*w1,(*s).s2); 
      _su3_multiply(chi,*w2,(*s).s3);
      _vector_add_assign(psi1,chi); 
      _su3_inverse_multiply(psi2,*w2,(*s).s2); 
      _su3_multiply(chi,*w3,(*s).s3);
      _vector_add_assign(psi2,chi); 

      // add in the twisted mass term (minus from g5 in the lower components)
      _vector_add_i_mul(psi1, -mu, (*s).s2);
      _vector_add_i_mul(psi2, -mu, (*s).s3);

      _vector_assign((*r).s2, psi1);
      _vector_assign((*r).s3, psi2);
    }
#ifdef OMP
  } /* OpenMP closing brace */
#endif
  return;
}

/**************************************************************
 * assign_mul_one_sw_pm_imu applies (1 + T + imug5) to 
 * a block spinor k
 *
 * it is assumed that the clover leaf is computed and stored
 * in sw[VOLUME][3][2]
 * the corresponding routine can be found in clover_leaf.c
 *
 **************************************************************/


void _PSWITCH(assign_mul_one_sw_pm_imu_block)(const int ieo, 
					      _PTSWITCH(spinor) * const k, _PTSWITCH(spinor) * const l,
					      const _F_TYPE mu, block *blk) {
  int jeo;
  _PTSWITCH(spinor) *r,*s;
  _PTSWITCH(spinor)  ALIGN tmpr;

  int t,x,y,z;     //lexiographic index of the site w.r.t the current process
  int bt,bx,by,bz; //block coordinate on the local mpi process
  int dT,dX,dY,dZ; //block size
  int sT,sX,sY,sZ; //constant shifts
  int lx; //lexiographic index of the block site w.r.t the local mpi process

  dT  = blk->BT;
  dX  = blk->BLX;
  dY  = blk->BLY;
  dZ  = blk->BLZ;

  bt = blk->mpilocal_coordinate[0];
  bx = blk->mpilocal_coordinate[1];
  by = blk->mpilocal_coordinate[2];
  bz = blk->mpilocal_coordinate[3];

  sT = bt*dT;
  sX = bx*dX;
  sY = by*dY;
  sZ = bz*dZ;
 

  r = k; 
  s = l;

  for(int it = 0; it < dT; it++)
    for(int ix = 0; ix < dX; ix++)
      for(int iy = 0; iy < dY; iy++)
        for(int iz = 0; iz < dZ; iz++)
          { 
            t = it + sT;
            x = ix + sX;
            y = iy + sY;
            z = iz + sZ;

            lx = g_ipt[t][x][y][z];

            jeo= (t+x+y+z)%2;

            if( ieo == jeo) {
	      _PSWITCH(assign_mul_one_sw_pm_imu_site_lexic)(lx, &tmpr, s, (_F_TYPE) g_mu);
	      _PSWITCH(assign)(r, &tmpr, 1);
	      r++;
	      s++;
	    }
          }

  return;
}



void _PSWITCH(assign_mul_one_sw_pm_imu_inv)(const int ieo, 
					    _PTSWITCH(spinor) * const k, _PTSWITCH(spinor) * const l,
					    const _F_TYPE mu) {
#ifdef OMP
#pragma omp parallel
  {
#endif
    _PTSWITCH(su3_vector) ALIGN psi, chi, phi1, phi3;
    const _PSWITCH(su3) *w1, *w2, *w3, *w4;
    const _PTSWITCH(spinor) *rn;
    _PTSWITCH(spinor) *s;

    /************************ loop over all lattice sites *************************/
#ifdef OMP
#pragma omp for
#endif
    for(int icx = 0; icx < (VOLUME/2); icx++) {

      rn = l + icx;
      s = k + icx;
      _vector_assign(phi1,(*rn).s0);
      _vector_assign(phi3,(*rn).s2);

      w1=&_PSWITCH(sw_inv)[icx][0][0];
      w2=w1+2;  /* &sw_inv[icx][1][0]; */
      w3=w1+4;  /* &sw_inv[icx][2][0]; */
      w4=w1+6;  /* &sw_inv[icx][3][0]; */
      _su3_multiply(psi,*w1,phi1); 
      _su3_multiply(chi,*w2,(*rn).s1);
      _vector_add((*s).s0,psi,chi);
      _su3_multiply(psi,*w4,phi1); 
      _su3_multiply(chi,*w3,(*rn).s1);
      _vector_add((*s).s1,psi,chi);

      w1++; /* &sw_inv[icx][0][1]; */
      w2++; /* &sw_inv[icx][1][1]; */
      w3++; /* &sw_inv[icx][2][1]; */
      w4++; /* &sw_inv[icx][3][1]; */
      _su3_multiply(psi,*w1,phi3); 
      _su3_multiply(chi,*w2,(*rn).s3);
      _vector_add((*s).s2,psi,chi);
      _su3_multiply(psi,*w4,phi3); 
      _su3_multiply(chi,*w3,(*rn).s3);
      _vector_add((*s).s3,psi,chi);

      /******************************** end of loop *********************************/
    }
#ifdef OMP
  } /* OpenMP closing brace */
#endif
  return;
}



void _PSWITCH(assign_mul_one_sw_pm_imu_inv_block)(const int ieo, 
						  _PTSWITCH(spinor) * const k, _PTSWITCH(spinor) * const l, 
						  const _F_TYPE mu, block *blk) {
  int i,it,ix,iy,iz, lxeo;
  _PTSWITCH(su3_vector) ALIGN psi, chi, phi1, phi3;
  const _PSWITCH(su3) *w1, *w2, *w3, *w4;
  const _PTSWITCH(spinor) *rn;
  _PTSWITCH(spinor) *s;

  rn = l; 
  s  = k;

  for(int t = 0; t < blk->BT; t++) {
    it = t + blk->mpilocal_coordinate[0]*blk->BT;
    for(int x = 0; x < blk->BLX; x++) {
      ix = x +  blk->mpilocal_coordinate[1]*blk->BLX;
      for(int y = 0; y < blk->BLY; y++) {
        iy = y +  blk->mpilocal_coordinate[2]*blk->BLY;
        for(int z = 0; z < blk->BLZ; z++) {
          iz = z +  blk->mpilocal_coordinate[3]*blk->BLZ;
          i = g_ipt[it][ix][iy][iz];
          lxeo = g_lexic2eo[i]; 
          if((t+x+y+z)%2 == ieo) {
            _vector_assign(phi1, (*rn).s0);
            _vector_assign(phi3, (*rn).s2);
            
            w1 = &_PSWITCH(sw_inv)[lxeo][0][0];
            w2 = w1 + 2;  /* &sw_inv[lxeo][1][0]; */
            w3 = w1 + 4;  /* &sw_inv[lxeo][2][0]; */
            w4 = w1 + 6;  /* &sw_inv[lxeo][3][0]; */
            
            _su3_multiply(psi, *w1, phi1); 
            _su3_multiply(chi, *w2, (*rn).s1);
            _vector_add((*s).s0, psi, chi);
            _su3_multiply(psi, *w4, phi1); 
            _su3_multiply(chi, *w3, (*rn).s1);
            _vector_add((*s).s1, psi, chi);
            
            w1++; /* &sw_inv[lxeo][0][1]; */
            w2++; /* &sw_inv[lxeo][1][1]; */
            w3++; /* &sw_inv[lxeo][2][1]; */
            w4++; /* &sw_inv[lxeo][3][1]; */
            _su3_multiply(psi, *w1, phi3); 
            _su3_multiply(chi, *w2, (*rn).s3);
            _vector_add((*s).s2, psi, chi);
            _su3_multiply(psi, *w4, phi3); 
            _su3_multiply(chi, *w3, (*rn).s3);
            _vector_add((*s).s3, psi, chi);
            rn++;
            s++;
          }
        }
      }
    }
  }
  return;
}

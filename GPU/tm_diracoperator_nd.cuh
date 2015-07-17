


__global__ void dev_sub_epsbar_tau1(dev_spinor* k_strange, dev_spinor* k_charm, 
                                    dev_spinor* l_strange, dev_spinor* l_charm){
   float4 slocal_strange[6], slocal_charm[6];
   float4 sout_strange[6], sout_charm[6];
   int pos;
   pos= threadIdx.x + blockDim.x*blockIdx.x; 

   if(pos < dev_VOLUME){
     dev_read_spinor(&(slocal_strange[0]), &(k_strange[pos]));
     dev_read_spinor(&(slocal_charm[0]), &(k_charm[pos]));
     dev_read_spinor(&(sout_strange[0]), &(l_strange[pos]));
     dev_read_spinor(&(sout_charm[0]), &(l_charm[pos]));

     sout_strange[0].x -= epsbar*slocal_charm[0].x;
     sout_strange[0].y -= epsbar*slocal_charm[0].y;
     sout_strange[0].z -= epsbar*slocal_charm[0].z;
     sout_strange[0].w -= epsbar*slocal_charm[0].w;

     sout_strange[1].x -= epsbar*slocal_charm[1].x;
     sout_strange[1].y -= epsbar*slocal_charm[1].y;
     sout_strange[1].z -= epsbar*slocal_charm[1].z;
     sout_strange[1].w -= epsbar*slocal_charm[1].w;

     sout_strange[2].x -= epsbar*slocal_charm[2].x;
     sout_strange[2].y -= epsbar*slocal_charm[2].y;
     sout_strange[2].z -= epsbar*slocal_charm[2].z;
     sout_strange[2].w -= epsbar*slocal_charm[2].w;

     sout_strange[3].x -= epsbar*slocal_charm[3].x;
     sout_strange[3].y -= epsbar*slocal_charm[3].y;
     sout_strange[3].z -= epsbar*slocal_charm[3].z;
     sout_strange[3].w -= epsbar*slocal_charm[3].w;

     sout_strange[4].x -= epsbar*slocal_charm[4].x;
     sout_strange[4].y -= epsbar*slocal_charm[4].y;
     sout_strange[4].z -= epsbar*slocal_charm[4].z;
     sout_strange[4].w -= epsbar*slocal_charm[4].w;

     sout_strange[5].x -= epsbar*slocal_charm[5].x;
     sout_strange[5].y -= epsbar*slocal_charm[5].y;
     sout_strange[5].z -= epsbar*slocal_charm[5].z;
     sout_strange[5].w -= epsbar*slocal_charm[5].w;



     sout_charm[0].x -= epsbar*slocal_strange[0].x;
     sout_charm[0].y -= epsbar*slocal_strange[0].y;
     sout_charm[0].z -= epsbar*slocal_strange[0].z;
     sout_charm[0].w -= epsbar*slocal_strange[0].w;

     sout_charm[1].x -= epsbar*slocal_strange[1].x;
     sout_charm[1].y -= epsbar*slocal_strange[1].y;
     sout_charm[1].z -= epsbar*slocal_strange[1].z;
     sout_charm[1].w -= epsbar*slocal_strange[1].w;

     sout_charm[2].x -= epsbar*slocal_strange[2].x;
     sout_charm[2].y -= epsbar*slocal_strange[2].y;
     sout_charm[2].z -= epsbar*slocal_strange[2].z;
     sout_charm[2].w -= epsbar*slocal_strange[2].w;

     sout_charm[3].x -= epsbar*slocal_strange[3].x;
     sout_charm[3].y -= epsbar*slocal_strange[3].y;
     sout_charm[3].z -= epsbar*slocal_strange[3].z;
     sout_charm[3].w -= epsbar*slocal_strange[3].w;

     sout_charm[4].x -= epsbar*slocal_strange[4].x;
     sout_charm[4].y -= epsbar*slocal_strange[4].y;
     sout_charm[4].z -= epsbar*slocal_strange[4].z;
     sout_charm[4].w -= epsbar*slocal_strange[4].w;

     sout_charm[5].x -= epsbar*slocal_strange[5].x;
     sout_charm[5].y -= epsbar*slocal_strange[5].y;
     sout_charm[5].z -= epsbar*slocal_strange[5].z;
     sout_charm[5].w -= epsbar*slocal_strange[5].w;


     dev_write_spinor(&(sout_strange[0]), &(l_strange[pos]));
     dev_write_spinor(&(sout_charm[0]), &(l_charm[pos]));
   }   
}


__global__ void dev_sub_epsbar_tau1_d(dev_spinor_d* k_strange, dev_spinor_d* k_charm, 
                                      dev_spinor_d* l_strange, dev_spinor_d* l_charm){
   double4 slocal_strange[6], slocal_charm[6];
   double4 sout_strange[6], sout_charm[6];
   int pos;
   pos= threadIdx.x + blockDim.x*blockIdx.x; 

   if(pos < dev_VOLUME){
     dev_read_spinor_d(&(slocal_strange[0]), &(k_strange[pos]));
     dev_read_spinor_d(&(slocal_charm[0]), &(k_charm[pos]));
     dev_read_spinor_d(&(sout_strange[0]), &(l_strange[pos]));
     dev_read_spinor_d(&(sout_charm[0]), &(l_charm[pos]));

     sout_strange[0].x -= epsbar_d*slocal_charm[0].x;
     sout_strange[0].y -= epsbar_d*slocal_charm[0].y;
     sout_strange[0].z -= epsbar_d*slocal_charm[0].z;
     sout_strange[0].w -= epsbar_d*slocal_charm[0].w;

     sout_strange[1].x -= epsbar_d*slocal_charm[1].x;
     sout_strange[1].y -= epsbar_d*slocal_charm[1].y;
     sout_strange[1].z -= epsbar_d*slocal_charm[1].z;
     sout_strange[1].w -= epsbar_d*slocal_charm[1].w;

     sout_strange[2].x -= epsbar_d*slocal_charm[2].x;
     sout_strange[2].y -= epsbar_d*slocal_charm[2].y;
     sout_strange[2].z -= epsbar_d*slocal_charm[2].z;
     sout_strange[2].w -= epsbar_d*slocal_charm[2].w;

     sout_strange[3].x -= epsbar_d*slocal_charm[3].x;
     sout_strange[3].y -= epsbar_d*slocal_charm[3].y;
     sout_strange[3].z -= epsbar_d*slocal_charm[3].z;
     sout_strange[3].w -= epsbar_d*slocal_charm[3].w;

     sout_strange[4].x -= epsbar_d*slocal_charm[4].x;
     sout_strange[4].y -= epsbar_d*slocal_charm[4].y;
     sout_strange[4].z -= epsbar_d*slocal_charm[4].z;
     sout_strange[4].w -= epsbar_d*slocal_charm[4].w;

     sout_strange[5].x -= epsbar_d*slocal_charm[5].x;
     sout_strange[5].y -= epsbar_d*slocal_charm[5].y;
     sout_strange[5].z -= epsbar_d*slocal_charm[5].z;
     sout_strange[5].w -= epsbar_d*slocal_charm[5].w;


     sout_charm[0].x -= epsbar_d*slocal_strange[0].x;
     sout_charm[0].y -= epsbar_d*slocal_strange[0].y;
     sout_charm[0].z -= epsbar_d*slocal_strange[0].z;
     sout_charm[0].w -= epsbar_d*slocal_strange[0].w;

     sout_charm[1].x -= epsbar_d*slocal_strange[1].x;
     sout_charm[1].y -= epsbar_d*slocal_strange[1].y;
     sout_charm[1].z -= epsbar_d*slocal_strange[1].z;
     sout_charm[1].w -= epsbar_d*slocal_strange[1].w;

     sout_charm[2].x -= epsbar_d*slocal_strange[2].x;
     sout_charm[2].y -= epsbar_d*slocal_strange[2].y;
     sout_charm[2].z -= epsbar_d*slocal_strange[2].z;
     sout_charm[2].w -= epsbar_d*slocal_strange[2].w;

     sout_charm[3].x -= epsbar_d*slocal_strange[3].x;
     sout_charm[3].y -= epsbar_d*slocal_strange[3].y;
     sout_charm[3].z -= epsbar_d*slocal_strange[3].z;
     sout_charm[3].w -= epsbar_d*slocal_strange[3].w;

     sout_charm[4].x -= epsbar_d*slocal_strange[4].x;
     sout_charm[4].y -= epsbar_d*slocal_strange[4].y;
     sout_charm[4].z -= epsbar_d*slocal_strange[4].z;
     sout_charm[4].w -= epsbar_d*slocal_strange[4].w;

     sout_charm[5].x -= epsbar_d*slocal_strange[5].x;
     sout_charm[5].y -= epsbar_d*slocal_strange[5].y;
     sout_charm[5].z -= epsbar_d*slocal_strange[5].z;
     sout_charm[5].w -= epsbar_d*slocal_strange[5].w;

     dev_write_spinor_d(&(sout_strange[0]), &(l_strange[pos]));
     dev_write_spinor_d(&(sout_charm[0]), &(l_charm[pos]));
   }   
}




void dev_Q_pm_ndpsi(dev_spinor * const l_strange, dev_spinor * const l_charm, 
                    dev_spinor * const k_strange, dev_spinor * const k_charm)
{

  int N_floats = 24*VOLUME;		// #floats

  //D_h^{dagger}
  //tau^1 by s<->c
  
     //D_psi(l_strange, k_charm);
     #ifdef USETEXTURE
       bind_texture_spin(k_charm,1);
     #endif
     dev_tm_dirac_kappa <<<gpu_gd_M, gpu_bd_M >>> (dev_gf, k_charm, l_strange, dev_nn);
     #ifdef USETEXTURE
      unbind_texture_spin(1);
     #endif

     //g_mu = -g_mu;     
     dev_swapmu <<<1,1>>> ();

     //D_psi(l_charm, k_strange);  
     #ifdef USETEXTURE
       bind_texture_spin(k_strange,1);
     #endif
     dev_tm_dirac_kappa <<<gpu_gd_M, gpu_bd_M >>> (dev_gf, k_strange, l_charm, dev_nn);   
     #ifdef USETEXTURE
      unbind_texture_spin(1);
     #endif

     //g_mu = -g_mu;
     dev_swapmu <<<1,1>>> (); 
   
     //sub_epsbar_tau1(l_strange, l_charm, k_charm, k_strange);
     dev_sub_epsbar_tau1<<<gpu_gd_linalg, gpu_bd_linalg >>>(k_charm, k_strange, l_strange, l_charm);
   
     //gamma5(g_spinor_field[DUM_MATRIX], l_strange, VOLUME);
     //gamma5(g_spinor_field[DUM_MATRIX+1], l_charm, VOLUME); 
  #ifdef RELATIVISTIC_BASIS 
     dev_gamma5_rel <<<gpu_gd_blas, gpu_bd_blas >>> (l_strange,dev_spin_eo1_up);
     dev_gamma5_rel <<<gpu_gd_blas, gpu_bd_blas >>> (l_charm,dev_spin_eo1_dn);
  #else
     dev_gamma5 <<<gpu_gd_blas, gpu_bd_blas >>> (l_strange,dev_spin_eo1_up);
     dev_gamma5 <<<gpu_gd_blas, gpu_bd_blas >>> (l_charm,dev_spin_eo1_dn);    
  #endif
   
     
    //D_h
    //tau^1 by s<->c     
     //D_psi(l_strange, g_spinor_field[DUM_MATRIX+1]);
     #ifdef USETEXTURE
       bind_texture_spin(dev_spin_eo1_dn,1);
     #endif
     dev_tm_dirac_kappa <<<gpu_gd_M, gpu_bd_M >>> (dev_gf, dev_spin_eo1_dn, l_strange, dev_nn);   
     #ifdef USETEXTURE
      unbind_texture_spin(1);
     #endif

     //g_mu = -g_mu;
     dev_swapmu <<<1,1>>> (); 

     //D_psi(l_charm, g_spinor_field[DUM_MATRIX]);  
     #ifdef USETEXTURE
       bind_texture_spin(dev_spin_eo1_up,1);
     #endif
     dev_tm_dirac_kappa <<<gpu_gd_M, gpu_bd_M >>> (dev_gf, dev_spin_eo1_up, l_charm, dev_nn);   
     #ifdef USETEXTURE
      unbind_texture_spin(1);
     #endif    

     //g_mu = -g_mu;
     dev_swapmu <<<1,1>>> ();

     //sub_epsbar_tau1(l_strange, l_charm, g_spinor_field[DUM_MATRIX+1], g_spinor_field[DUM_MATRIX]);
     dev_sub_epsbar_tau1<<<gpu_gd_linalg, gpu_bd_linalg >>>(dev_spin_eo1_dn, dev_spin_eo1_up, l_strange, l_charm);
     
     //gamma5(l_strange, l_strange, VOLUME);      
     //gamma5(l_charm, l_charm, VOLUME);
  #ifdef RELATIVISTIC_BASIS 
     dev_gamma5_rel <<<gpu_gd_blas, gpu_bd_blas >>> (l_strange, l_strange);
     dev_gamma5_rel <<<gpu_gd_blas, gpu_bd_blas >>> (l_charm, l_charm);
  #else
     dev_gamma5 <<<gpu_gd_blas, gpu_bd_blas >>> (l_strange, l_strange);
     dev_gamma5 <<<gpu_gd_blas, gpu_bd_blas >>> (l_charm, l_charm);    
  #endif

     /* At the end, the normalisation by the max. eigenvalue  */ 
     /* Twice  phmc_invmaxev  since we consider here  D Ddag  !!! */
     //mul_r(l_charm, phmc_invmaxev*phmc_invmaxev, l_charm, VOLUME);
     //mul_r(l_strange, phmc_invmaxev*phmc_invmaxev, l_strange, VOLUME);     
     cublasSscal (N_floats, phmc_invmaxev*phmc_invmaxev, (float *) l_charm, 1);
     cublasSscal (N_floats, phmc_invmaxev*phmc_invmaxev, (float *) l_strange, 1);
}







void dev_Q_pm_ndpsi_d(dev_spinor_d * const l_strange, dev_spinor_d * const l_charm, 
                                 dev_spinor_d * const k_strange, dev_spinor_d * const k_charm)
{


  //D_h^{dagger}
  //tau^1 by s<->c
  
     //D_psi(l_strange, k_charm);
     dev_tm_dirac_kappa_d <<<gpu_gd_M_d, gpu_bd_M_d >>> (dev_gf_d, k_charm, l_strange, dev_nn);
     //g_mu = -g_mu;     
     dev_swapmu_d <<<1,1>>> ();

     //D_psi(l_charm, k_strange);  
     dev_tm_dirac_kappa_d <<<gpu_gd_M_d, gpu_bd_M_d >>> (dev_gf_d, k_strange, l_charm, dev_nn);   
     //g_mu = -g_mu;
     dev_swapmu_d <<<1,1>>> (); 
   
     //sub_epsbar_tau1(l_strange, l_charm, k_charm, k_strange);
     dev_sub_epsbar_tau1_d<<<gpu_gd_linalg_d, gpu_bd_linalg_d >>>(k_charm, k_strange, l_strange, l_charm);
   
     //gamma5(g_spinor_field[DUM_MATRIX], l_strange, VOLUME);
     //gamma5(g_spinor_field[DUM_MATRIX+1], l_charm, VOLUME); 
  #ifdef RELATIVISTIC_BASIS 
     dev_gamma5_rel_d <<<gpu_gd_blas_d, gpu_bd_blas_d >>> (l_strange,dev_spin_eo1_up_d);
     dev_gamma5_rel_d <<<gpu_gd_blas_d, gpu_bd_blas_d >>> (l_charm,dev_spin_eo1_dn_d);
  #else
     dev_gamma5_d <<<gpu_gd_blas_d, gpu_bd_blas_d >>> (l_strange,dev_spin_eo1_up_d);
     dev_gamma5_d <<<gpu_gd_blas_d, gpu_bd_blas_d >>> (l_charm,dev_spin_eo1_dn_d);    
  #endif
   
     
    //D_h
    //tau^1 by s<->c     
     //D_psi(l_strange, g_spinor_field[DUM_MATRIX+1]);
     dev_tm_dirac_kappa_d <<<gpu_gd_M_d, gpu_bd_M_d >>> (dev_gf_d, dev_spin_eo1_dn_d, l_strange, dev_nn);   

     //g_mu = -g_mu;
     dev_swapmu_d <<<1,1>>> (); 

     //D_psi(l_charm, g_spinor_field[DUM_MATRIX]);  
     dev_tm_dirac_kappa_d <<<gpu_gd_M_d, gpu_bd_M_d >>> (dev_gf_d, dev_spin_eo1_up_d, l_charm, dev_nn);   
       
     //g_mu = -g_mu;
     dev_swapmu_d <<<1,1>>> ();

     //sub_epsbar_tau1(l_strange, l_charm, g_spinor_field[DUM_MATRIX+1], g_spinor_field[DUM_MATRIX]);
     dev_sub_epsbar_tau1_d<<<gpu_gd_linalg_d, gpu_bd_linalg_d >>>(dev_spin_eo1_dn_d, dev_spin_eo1_up_d, l_strange, l_charm);
     
     //gamma5(l_strange, l_strange, VOLUME);      
     //gamma5(l_charm, l_charm, VOLUME);
  #ifdef RELATIVISTIC_BASIS 
     dev_gamma5_rel_d <<<gpu_gd_blas_d, gpu_bd_blas_d >>> (l_strange, l_strange);
     dev_gamma5_rel_d <<<gpu_gd_blas_d, gpu_bd_blas_d >>> (l_charm, l_charm);
  #else
     dev_gamma5_d <<<gpu_gd_blas_d, gpu_bd_blas_d >>> (l_strange, l_strange);
     dev_gamma5_d <<<gpu_gd_blas_d, gpu_bd_blas_d >>> (l_charm, l_charm);    
  #endif

     /* At the end, the normalisation by the max. eigenvalue  */ 
     /* Twice  phmc_invmaxev  since we consider here  D Ddag  !!! */     
     dev_blasscal_d<<<gpu_gd_blas_d, gpu_bd_blas_d>>>(phmc_invmaxev*phmc_invmaxev, l_charm);
     dev_blasscal_d<<<gpu_gd_blas_d, gpu_bd_blas_d>>>(phmc_invmaxev*phmc_invmaxev, l_strange);


}


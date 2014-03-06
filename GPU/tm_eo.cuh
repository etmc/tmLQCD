// computes sout = 1/(1 +- mutilde gamma5) sin = (1 -+ i mutilde gamma5)/(1+mutilde^2) sin
// mutilde = 2 kappa mu
__global__ void dev_mul_one_pm_imu_inv(dev_spinor* sin, dev_spinor* sout, const float sign){
   
   dev_spinor slocal[6];
   //need the inverse sign in the numerator because of inverse
   dev_complex pm_imu = dev_initcomplex(0.0,-1.0*sign*twokappamu);
   
   float one_plus_musquare_inv = 1.0/(1.0 + twokappamu*twokappamu);
   int pos;
   pos= threadIdx.x + blockDim.x*blockIdx.x;  

   if(pos < dev_VOLUME){
     #ifdef RELATIVISTIC_BASIS
       dev_skalarmult_gamma5_globalspinor_rel(&(slocal[0]), pm_imu, &(sin[pos]) );
     #else
       dev_skalarmult_gamma5_globalspinor(&(slocal[0]), pm_imu, &(sin[pos]) );
     #endif
     dev_add_globalspinor_assign(&(slocal[0]), &(sin[pos])); 
     dev_realmult_spinor_assigntoglobal(&(sout[pos]), one_plus_musquare_inv, &(slocal[0]) );
   }
}





// sout = gamma_5*((1\pm i\mutilde \gamma_5)*sin1 - sin2)
__global__ void dev_mul_one_pm_imu_sub_mul_gamma5(dev_spinor* sin1, dev_spinor* sin2, dev_spinor* sout, const float sign){
   dev_spinor slocal[6];
   dev_complex pm_imu = dev_initcomplex(0.0, sign*twokappamu); // i mutilde
   int pos;
   pos= threadIdx.x + blockDim.x*blockIdx.x; 

   if(pos < dev_VOLUME){
     #ifdef RELATIVISTIC_BASIS
       dev_skalarmult_gamma5_globalspinor_rel(&(slocal[0]), pm_imu, &(sin1[pos]) );
     #else
       dev_skalarmult_gamma5_globalspinor(&(slocal[0]),pm_imu,&(sin1[pos]));
     #endif
     dev_add_globalspinor_assign(&(slocal[0]), &(sin1[pos]));
     dev_sub_globalspinor_assign(&(slocal[0]), &(sin2[pos]));
     #ifdef RELATIVISTIC_BASIS
       dev_Gamma5_assigntoglobal_rel(&(sout[pos]), &(slocal[0]));
     #else
       dev_Gamma5_assigntoglobal(&(sout[pos]), &(slocal[0]));
     #endif
   }   
}













// aequivalent to Qtm_pm_psi in tm_operators.c, this is NON-MPI version
extern "C" void dev_Qtm_pm_psi_old(dev_spinor* spinin, dev_spinor* spinout, int gridsize, dim3 blocksize, int gridsize2, int blocksize2){
  //spinin == odd
  //spinout == odd
  
  int VolumeEO = VOLUME/2;
  //Q_{-}
  #ifdef USETEXTURE
    bind_texture_spin(spinin,1);
  #endif
  //bind_texture_nn(dev_nn_eo);
    dev_Hopping_Matrix<<<gridsize, blocksize>>>
             (dev_gf, spinin, dev_spin_eo1, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0, 0, VolumeEO); //dev_spin_eo1 == even -> 0           
  //unbind_texture_nn();           
  #ifdef USETEXTURE
    unbind_texture_spin(1);
  #endif
  dev_mul_one_pm_imu_inv<<<gridsize2, blocksize2>>>(dev_spin_eo1,dev_spin_eo2, -1.);
  

  #ifdef USETEXTURE
    bind_texture_spin(dev_spin_eo2,1);
  #endif
  //bind_texture_nn(dev_nn_oe);
    dev_Hopping_Matrix<<<gridsize, blocksize>>>
            (dev_gf, dev_spin_eo2, dev_spin_eo1, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1, 0, VolumeEO); 
  //unbind_texture_nn();
  #ifdef USETEXTURE
    unbind_texture_spin(1);
  #endif
  dev_mul_one_pm_imu_sub_mul_gamma5<<<gridsize2, blocksize2>>>(spinin, dev_spin_eo1,  dev_spin_eo2, -1.);
  
  
  //Q_{+}

  #ifdef USETEXTURE
    bind_texture_spin(dev_spin_eo2,1);
  #endif
  //bind_texture_nn(dev_nn_eo);
    dev_Hopping_Matrix<<<gridsize, blocksize>>>
          (dev_gf, dev_spin_eo2, dev_spin_eo1, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0, 0, VolumeEO); //dev_spin_eo1 == even -> 0
  //unbind_texture_nn();      
  #ifdef USETEXTURE  
    unbind_texture_spin(1);
  #endif
  dev_mul_one_pm_imu_inv<<<gridsize2, blocksize2>>>(dev_spin_eo1,spinout, +1.);
  

  #ifdef USETEXTURE
    bind_texture_spin(spinout,1);
  #endif
  //bind_texture_nn(dev_nn_oe);
    dev_Hopping_Matrix<<<gridsize, blocksize>>>
             (dev_gf, spinout, dev_spin_eo1, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1, 0, VolumeEO); 
  //unbind_texture_nn();  
  #ifdef USETEXTURE
    unbind_texture_spin(1);
  #endif
  dev_mul_one_pm_imu_sub_mul_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo2, dev_spin_eo1,  spinout , +1.); 
}





// aequivalent to Qtm_pm_psi in tm_operators.c, this is NON-MPI version
// fused hopping and linalg kernels to be more efficient on kepler
extern "C" void dev_Qtm_pm_psi_scalar(dev_spinor* spinin, dev_spinor* spinout, int gridsize, dim3 blocksize, int gridsize2, int blocksize2){
  //spinin == odd
  //spinout == odd
  
  int VolumeEO = VOLUME/2;
  //Q_{-}
  #ifdef USETEXTURE
    bind_texture_spin(spinin,1);
  #endif
    dev_Hopping_Matrix_ext3<<<gridsize, blocksize>>>
             (dev_gf, spinin, spinout, -1.0, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0, 0, VolumeEO); //dev_spin_eo1 == even -> 0           
  //unbind_texture_nn();           
  #ifdef USETEXTURE
    unbind_texture_spin(1);
  #endif

  

  #ifdef USETEXTURE
    bind_texture_spin(spinout,1);
  #endif
    dev_Hopping_Matrix_ext2<<<gridsize, blocksize>>>
            (dev_gf, spinout, dev_spin_eo2, -1.0, spinin, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1, 0, VolumeEO); 
  //unbind_texture_nn();
  #ifdef USETEXTURE
    unbind_texture_spin(1);
  #endif

  
  //Q_{+}

  #ifdef USETEXTURE
    bind_texture_spin(dev_spin_eo2,1);
  #endif
    dev_Hopping_Matrix_ext3<<<gridsize, blocksize>>>
          (dev_gf, dev_spin_eo2, dev_spin_eo1, +1.0, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0, 0, VolumeEO); //dev_spin_eo1 == even -> 0
  //unbind_texture_nn();      
  #ifdef USETEXTURE  
    unbind_texture_spin(1);
  #endif

  

  #ifdef USETEXTURE
    bind_texture_spin(dev_spin_eo1,1);
  #endif
  //bind_texture_nn(dev_nn_oe);
    dev_Hopping_Matrix_ext2<<<gridsize, blocksize>>>
             (dev_gf, dev_spin_eo1, spinout, 1.0, dev_spin_eo2,  dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1, 0, VolumeEO); 
  //unbind_texture_nn();  
  #ifdef USETEXTURE
    unbind_texture_spin(1);
  #endif

}






#ifdef MPI
// aequivalent to Qtm_pm_psi in tm_operators.c
// using HOPPING_ASYNC for mpi
extern "C" void dev_Qtm_pm_psi_mpi(dev_spinor* spinin, dev_spinor* spinout, int gridsize, dim3 blocksize, int gridsize2, int blocksize2){
  //spinin == odd
  //spinout == odd
  
  //Q_{-}


    HOPPING_ASYNC(dev_gf, spinin, dev_spin_eo1, dev_eoidx_even, 
               dev_eoidx_odd, dev_nn_eo, 0,gridsize, blocksize); //dev_spin_eo1 == even -> 0           
          


  dev_mul_one_pm_imu_inv<<<gridsize2, blocksize2>>>(dev_spin_eo1,dev_spin_eo2, -1.);
  



    HOPPING_ASYNC(dev_gf, dev_spin_eo2, dev_spin_eo1, 
          dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1,gridsize, 
          blocksize); 

  dev_mul_one_pm_imu_sub_mul_gamma5<<<gridsize2, blocksize2>>>(spinin, dev_spin_eo1,  dev_spin_eo2, -1.);
  
  
  //Q_{+}

    HOPPING_ASYNC(dev_gf, dev_spin_eo2, dev_spin_eo1, 
         dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0, gridsize, 
         blocksize); //dev_spin_eo1 == even -> 0

  dev_mul_one_pm_imu_inv<<<gridsize2, blocksize2>>>(dev_spin_eo1,spinout, +1.);
  

    HOPPING_ASYNC(dev_gf, spinout, dev_spin_eo1, dev_eoidx_odd, 
           dev_eoidx_even, dev_nn_oe, 1,gridsize, blocksize); 

  dev_mul_one_pm_imu_sub_mul_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo2, dev_spin_eo1,  spinout , +1.); 
}
#endif




extern "C" void dev_Qtm_pm_psi(dev_spinor* spinin, dev_spinor* spinout, int gridsize, dim3 blocksize, 
			       int gridsize2, int blocksize2){

#ifdef MPI
  dev_Qtm_pm_psi_mpi(spinin, spinout, gridsize, blocksize, gridsize2, blocksize2);
#else
  dev_Qtm_pm_psi_scalar(spinin, spinout, gridsize, blocksize, gridsize2, blocksize2);
#endif

}



#ifdef HALF

// computes sout = 1/(1 +- mutilde gamma5) sin = (1 -+ i mutilde gamma5)/(1+mutilde^2) sin
// mutilde = 2 kappa mu
// uses shared local memory for manipulation
__global__ void dev_mul_one_pm_imu_inv_half(dev_spinor_half* sin, float* sin_norm, dev_spinor_half* sout, float* sout_norm, const float sign){
   
   dev_spinor slocal[6];
   dev_spinor s[6];
   float norm;
   
   //need the inverse sign in the numerator because of inverse
   dev_complex pm_imu = dev_initcomplex(0.0,-1.0*sign*twokappamu);
   
   float one_plus_musquare_inv = 1.0/(1.0 + twokappamu*twokappamu);
   int pos;
   pos= threadIdx.x + blockDim.x*blockIdx.x;  
   int ix = threadIdx.x;
   if(pos < dev_VOLUME){
     norm = sin_norm[pos];
     construct_spinor_fromhalf(s, sin, norm, pos);
     #ifdef RELATIVISTIC_BASIS
       dev_skalarmult_gamma5_spinor_rel(&(slocal[0]), pm_imu, &(s[0]) );
     #else
       dev_skalarmult_gamma5_spinor(&(slocal[0]), pm_imu, &(s[0]) );
     #endif
     dev_add_spinor_assign(&(slocal[0]), &(s[0]));
     
     dev_realmult_spinor_assign(&(s[0]), one_plus_musquare_inv, &(slocal[0]) );
     
     dev_write_spinor_half(&(s[0]),&(sout[pos]), &(sout_norm[pos]));
   }
}





// sout = gamma_5*((1\pm i\mutilde \gamma_5)*sin1 - sin2)
// uses shared local memory for manipulation
__global__ void dev_mul_one_pm_imu_sub_mul_gamma5_half(dev_spinor_half* sin1, float* sin1_norm, dev_spinor_half* sin2, float* sin2_norm, dev_spinor_half* sout, float* sout_norm, const float sign){
   dev_spinor slocal[6];
   dev_spinor s1[6];
   dev_spinor s2[6];
   float norm;
   dev_complex pm_imu = dev_initcomplex(0.0, sign*twokappamu); // i mutilde
   int pos;
   pos= threadIdx.x + blockDim.x*blockIdx.x; 
   int ix = threadIdx.x;
   if(pos < dev_VOLUME){
     norm = sin1_norm[pos];
     construct_spinor_fromhalf(s1, sin1,norm, pos);
     norm = sin2_norm[pos];
     construct_spinor_fromhalf(s2, sin2, norm, pos);

     #ifdef RELATIVISTIC_BASIS
       dev_skalarmult_gamma5_spinor_rel(&(slocal[0]),pm_imu,&(s1[0]));
     #else
       dev_skalarmult_gamma5_spinor(&(slocal[0]),pm_imu,&(s1[0]));
     #endif
     
     dev_add_spinor_assign(&(slocal[0]), &(s1[0]));
     dev_sub_spinor_assign(&(slocal[0]), &(s2[0]));
     
     #ifdef RELATIVISTIC_BASIS
       dev_Gamma5_assign_rel(&(s1[0]), &(slocal[0]));
     #else
       dev_Gamma5_assign(&(s1[0]), &(slocal[0]));
     #endif
     
     dev_write_spinor_half(&(s1[0]),&(sout[pos]), &(sout_norm[pos]));
   }   
}





// aequivalent to Qtm_pm_psi in tm_operators.c for half precision
extern "C" void dev_Qtm_pm_psi_half(dev_spinor_half* spinin, float* spinin_norm, dev_spinor_half* spinout, float* spinout_norm, int gridsize, dim3 blocksize, int gridsize2, int blocksize2){
  //spinin == odd
  //spinout == odd
  
  //Q_{-}
  #ifdef USETEXTURE
    bind_halfspinor_texture(spinin, spinin_norm);
  #endif
    dev_Hopping_Matrix_half<<<gridsize, blocksize>>>
             (dev_gf_half, spinin, spinin_norm, dev_spin_eo1_half, dev_spin_eo1_half_norm, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0); //dev_spin_eo1 == even -> 0  
  #ifdef USETEXTURE
    unbind_halfspinor_texture();
  #endif
  dev_mul_one_pm_imu_inv_half<<<gridsize2, blocksize2>>>(dev_spin_eo1_half, dev_spin_eo1_half_norm ,dev_spin_eo2_half, dev_spin_eo2_half_norm, -1.);
  
  #ifdef USETEXTURE
    bind_halfspinor_texture(dev_spin_eo2_half, dev_spin_eo2_half_norm);
  #endif
    dev_Hopping_Matrix_half<<<gridsize, blocksize>>>
            (dev_gf_half, dev_spin_eo2_half, dev_spin_eo2_half_norm, dev_spin_eo1_half, dev_spin_eo1_half_norm, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1); 
  #ifdef USETEXTURE
    unbind_halfspinor_texture();
  #endif
  dev_mul_one_pm_imu_sub_mul_gamma5_half<<<gridsize2, blocksize2>>>(spinin, spinin_norm, dev_spin_eo1_half, dev_spin_eo1_half_norm,  dev_spin_eo2_half, dev_spin_eo2_half_norm, -1.);
  
  //Q_{+}
  #ifdef USETEXTURE
    bind_halfspinor_texture(dev_spin_eo2_half, dev_spin_eo2_half_norm);
  #endif
    dev_Hopping_Matrix_half<<<gridsize, blocksize>>>
          (dev_gf_half, dev_spin_eo2_half, dev_spin_eo2_half_norm, dev_spin_eo1_half, dev_spin_eo1_half_norm, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0); //dev_spin_eo1 == even -> 0    
  #ifdef USETEXTURE  
    unbind_halfspinor_texture();
  #endif
  dev_mul_one_pm_imu_inv_half<<<gridsize2, blocksize2>>>(dev_spin_eo1_half, dev_spin_eo1_half_norm,spinout, spinout_norm, +1.);
  
  #ifdef USETEXTURE
    bind_halfspinor_texture(spinout, spinout_norm);
  #endif
    dev_Hopping_Matrix_half<<<gridsize, blocksize>>>
             (dev_gf_half, spinout, spinout_norm, dev_spin_eo1_half, dev_spin_eo1_half_norm, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1);  
  #ifdef USETEXTURE
    unbind_halfspinor_texture();
  #endif
  dev_mul_one_pm_imu_sub_mul_gamma5_half<<<gridsize2, blocksize2>>>(dev_spin_eo2_half, dev_spin_eo2_half_norm, dev_spin_eo1_half, dev_spin_eo1_half_norm,  spinout, spinout_norm , +1.); 
}


#ifdef MPI

// aequivalent to Qtm_pm_psi in tm_operators.c for half precision
extern "C" void dev_Qtm_pm_psi_half_mpi(dev_spinor_half* spinin, float* spinin_norm, dev_spinor_half* spinout, float* spinout_norm, int gridsize, dim3 blocksize, int gridsize2, int blocksize2){
  //spinin == odd
  //spinout == odd
  
  //Q_{-}
  HOPPING_HALF_ASYNC(dev_gf_half, spinin, spinin_norm, dev_spin_eo1_half, dev_spin_eo1_half_norm, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0,gridsize, blocksize); //dev_spin_eo1 == even -> 0  

  dev_mul_one_pm_imu_inv_half<<<gridsize2, blocksize2>>>(dev_spin_eo1_half, dev_spin_eo1_half_norm ,dev_spin_eo2_half, dev_spin_eo2_half_norm, -1.);
  

    HOPPING_HALF_ASYNC(dev_gf_half, dev_spin_eo2_half, dev_spin_eo2_half_norm, dev_spin_eo1_half, dev_spin_eo1_half_norm, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1,gridsize, blocksize); 

  dev_mul_one_pm_imu_sub_mul_gamma5_half<<<gridsize2, blocksize2>>>(spinin, spinin_norm, dev_spin_eo1_half, dev_spin_eo1_half_norm,  dev_spin_eo2_half, dev_spin_eo2_half_norm, -1.);
  
  //Q_{+}
    HOPPING_HALF_ASYNC (dev_gf_half, dev_spin_eo2_half, dev_spin_eo2_half_norm, dev_spin_eo1_half, dev_spin_eo1_half_norm, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0,gridsize, blocksize); //dev_spin_eo1 == even -> 0    
    
  dev_mul_one_pm_imu_inv_half<<<gridsize2, blocksize2>>>(dev_spin_eo1_half, dev_spin_eo1_half_norm,spinout, spinout_norm, +1.);
  
    HOPPING_HALF_ASYNC (dev_gf_half, spinout, spinout_norm, dev_spin_eo1_half, dev_spin_eo1_half_norm, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1,gridsize, blocksize);  

  dev_mul_one_pm_imu_sub_mul_gamma5_half<<<gridsize2, blocksize2>>>(dev_spin_eo2_half, dev_spin_eo2_half_norm, dev_spin_eo1_half, dev_spin_eo1_half_norm,  spinout, spinout_norm , +1.); 
}
#endif // MPI





/*
extern "C" void dev_Qtm_pm_psi(dev_spinor* spinin, dev_spinor* spinout, int gridsize, int blocksize, int gridsize2, int blocksize2){

  printf("WARNING: dummy function 'dev_Qtm_pm_psi' was called\n");
  
}
*/





#endif //HALF

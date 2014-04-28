
/////////////
// KERNELS //
/////////////

// derived from function  dev_mul_one_pm_imu_inv
//	order of the arguments (spinin, spinout)
// applies (1 +- imubar*gamma5)
// one thread per lattice site
__global__ void dev_mul_one_pm_imubar_gamma5 (dev_spinor * sin,
                                              dev_spinor * sout,
                                              float sign         ) {
   
  dev_spinor slocal[6];									// dev_spinor = float4		// 6*float4 = 24 floats		// auxiliary for each thread
  dev_complex pm_imu = dev_initcomplex(0.0, sign * mubar);				// dev_complex = struct { REAL re; REAL im; }	// pm_imu.re = 0.0																	// pm_imu.im = sign * mubar
  int pos = threadIdx.x + blockDim.x*blockIdx.x;
  
  if (pos < dev_VOLUME) {
    #ifdef RELATIVISTIC_BASIS
      dev_skalarmult_gamma5_globalspinor_rel(&(slocal[0]), pm_imu, &(sin[pos]) );			// slocal  =  pm_imu * (gamma5) * sin
    #else
      dev_skalarmult_gamma5_globalspinor(&(slocal[0]), pm_imu, &(sin[pos]) );			// slocal  =  pm_imu * (gamma5) * sin
    #endif
    dev_add_globalspinor_assign(&(slocal[0]), &(sin[pos]));					// slocal  =  slocal + sin  =  pm_imu * (gamma5) * sin + sin
    dev_realmult_spinor_assigntoglobal(&(sout[pos]), 1.0, &(slocal[0]) );			// sout    =  slocal
  }
}



//s2_up = s2_up -epsbar s3_dn - s1_up
//s2_dn = s2_dn -epsbar s3_up - s1_dn
__global__ void dev_nd_linalg1 (dev_spinor * s1_up, dev_spinor * s1_dn,
				dev_spinor * s3_up, dev_spinor * s3_dn,  
				dev_spinor * s2_up, dev_spinor * s2_dn,	float epsbar			
				) {
  dev_spinor s1[6], s3[6], sout[6];	
  int pos = threadIdx.x + blockDim.x*blockIdx.x;
  
  if (pos < dev_VOLUME) {
    //upper output spinor
    dev_read_spinor(&(sout[0]), &(s2_up[pos]));
    dev_read_spinor(&(s1[0]), &(s1_up[pos]));    
    dev_read_spinor(&(s3[0]), &(s3_up[pos])); 
    
    sout[0].x -= epsbar*s3[0].x;
    sout[0].x -= s1[0].x;
    sout[0].y -= epsbar*s3[0].y;
    sout[0].y -= s1[0].y;    
    sout[0].z -= epsbar*s3[0].z;
    sout[0].z -= s1[0].z;    
    sout[0].w -= epsbar*s3[0].w;
    sout[0].w -= s1[0].w;       
    
    sout[1].x -= epsbar*s3[1].x;
    sout[1].x -= s1[1].x;
    sout[1].y -= epsbar*s3[1].y;
    sout[1].y -= s1[1].y;    
    sout[1].z -= epsbar*s3[1].z;
    sout[1].z -= s1[1].z;    
    sout[1].w -= epsbar*s3[1].w;
    sout[1].w -= s1[1].w;       
    
    sout[2].x -= epsbar*s3[2].x;
    sout[2].x -= s1[2].x;
    sout[2].y -= epsbar*s3[2].y;
    sout[2].y -= s1[2].y;    
    sout[2].z -= epsbar*s3[2].z;
    sout[2].z -= s1[2].z;    
    sout[2].w -= epsbar*s3[2].w;
    sout[2].w -= s1[2].w;       
    
    sout[3].x -= epsbar*s3[3].x;
    sout[3].x -= s1[3].x;
    sout[3].y -= epsbar*s3[3].y;
    sout[3].y -= s1[3].y;    
    sout[3].z -= epsbar*s3[3].z;
    sout[3].z -= s1[3].z;    
    sout[3].w -= epsbar*s3[3].w;
    sout[3].w -= s1[3].w;   
    
    sout[4].x -= epsbar*s3[4].x;
    sout[4].x -= s1[4].x;
    sout[4].y -= epsbar*s3[4].y;
    sout[4].y -= s1[4].y;    
    sout[4].z -= epsbar*s3[4].z;
    sout[4].z -= s1[4].z;    
    sout[4].w -= epsbar*s3[4].w;
    sout[4].w -= s1[4].w;       
    
    sout[5].x -= epsbar*s3[5].x;
    sout[5].x -= s1[5].x;
    sout[5].y -= epsbar*s3[5].y;
    sout[5].y -= s1[5].y;    
    sout[5].z -= epsbar*s3[5].z;
    sout[5].z -= s1[5].z;    
    sout[5].w -= epsbar*s3[5].w;
    sout[5].w -= s1[5].w;       
    
    dev_write_spinor(&(sout[0]),&(s2_up[pos]));    
    
    
    //lower output spinor
    dev_read_spinor(&(sout[0]), &(s2_dn[pos]));
    dev_read_spinor(&(s1[0]), &(s1_dn[pos]));    
    dev_read_spinor(&(s3[0]), &(s3_dn[pos]));    
    

    sout[0].x -= epsbar*s3[0].x;
    sout[0].x -= s1[0].x;
    sout[0].y -= epsbar*s3[0].y;
    sout[0].y -= s1[0].y;    
    sout[0].z -= epsbar*s3[0].z;
    sout[0].z -= s1[0].z;    
    sout[0].w -= epsbar*s3[0].w;
    sout[0].w -= s1[0].w;       
    
    sout[1].x -= epsbar*s3[1].x;
    sout[1].x -= s1[1].x;
    sout[1].y -= epsbar*s3[1].y;
    sout[1].y -= s1[1].y;    
    sout[1].z -= epsbar*s3[1].z;
    sout[1].z -= s1[1].z;    
    sout[1].w -= epsbar*s3[1].w;
    sout[1].w -= s1[1].w;       
    
    sout[2].x -= epsbar*s3[2].x;
    sout[2].x -= s1[2].x;
    sout[2].y -= epsbar*s3[2].y;
    sout[2].y -= s1[2].y;    
    sout[2].z -= epsbar*s3[2].z;
    sout[2].z -= s1[2].z;    
    sout[2].w -= epsbar*s3[2].w;
    sout[2].w -= s1[2].w;       
    
    sout[3].x -= epsbar*s3[3].x;
    sout[3].x -= s1[3].x;
    sout[3].y -= epsbar*s3[3].y;
    sout[3].y -= s1[3].y;    
    sout[3].z -= epsbar*s3[3].z;
    sout[3].z -= s1[3].z;    
    sout[3].w -= epsbar*s3[3].w;
    sout[3].w -= s1[3].w;   
    
    sout[4].x -= epsbar*s3[4].x;
    sout[4].x -= s1[4].x;
    sout[4].y -= epsbar*s3[4].y;
    sout[4].y -= s1[4].y;    
    sout[4].z -= epsbar*s3[4].z;
    sout[4].z -= s1[4].z;    
    sout[4].w -= epsbar*s3[4].w;
    sout[4].w -= s1[4].w;       
    
    sout[5].x -= epsbar*s3[5].x;
    sout[5].x -= s1[5].x;
    sout[5].y -= epsbar*s3[5].y;
    sout[5].y -= s1[5].y;    
    sout[5].z -= epsbar*s3[5].z;
    sout[5].z -= s1[5].z;    
    sout[5].w -= epsbar*s3[5].w;
    sout[5].w -= s1[5].w;         
    
    
    dev_write_spinor(&(sout[0]),&(s2_dn[pos]));       
  }
}








//s2_up = gamma5*(s2_up -epsbar s3_up - s1_up)
//s2_dn = gamma5*(s2_dn -epsbar s3_dn - s1_dn)
__global__ void dev_nd_linalg1_gamma5 (dev_spinor * s1_up, dev_spinor * s1_dn,
				dev_spinor * s3_up, dev_spinor * s3_dn,  
				dev_spinor * s2_up, dev_spinor * s2_dn,	float epsbar			
				) {
  dev_spinor s1[6], s3[6], sout[6];	
  int pos = threadIdx.x + blockDim.x*blockIdx.x;
  
  if (pos < dev_VOLUME) {
    //upper output spinor
    dev_read_spinor(&(sout[0]), &(s2_up[pos]));
    dev_read_spinor(&(s1[0]), &(s1_up[pos]));    
    dev_read_spinor(&(s3[0]), &(s3_up[pos])); 
    
    sout[0].x -= epsbar*s3[0].x;
    sout[0].x -= s1[0].x;
    sout[0].y -= epsbar*s3[0].y;
    sout[0].y -= s1[0].y;    
    sout[0].z -= epsbar*s3[0].z;
    sout[0].z -= s1[0].z;    
    sout[0].w -= epsbar*s3[0].w;
    sout[0].w -= s1[0].w;       
    
    sout[1].x -= epsbar*s3[1].x;
    sout[1].x -= s1[1].x;
    sout[1].y -= epsbar*s3[1].y;
    sout[1].y -= s1[1].y;    
    sout[1].z -= epsbar*s3[1].z;
    sout[1].z -= s1[1].z;    
    sout[1].w -= epsbar*s3[1].w;
    sout[1].w -= s1[1].w;       
    
    sout[2].x -= epsbar*s3[2].x;
    sout[2].x -= s1[2].x;
    sout[2].y -= epsbar*s3[2].y;
    sout[2].y -= s1[2].y;    
    sout[2].z -= epsbar*s3[2].z;
    sout[2].z -= s1[2].z;    
    sout[2].w -= epsbar*s3[2].w;
    sout[2].w -= s1[2].w;       
    
    sout[3].x -= epsbar*s3[3].x;
    sout[3].x -= s1[3].x;
    sout[3].y -= epsbar*s3[3].y;
    sout[3].y -= s1[3].y;    
    sout[3].z -= epsbar*s3[3].z;
    sout[3].z -= s1[3].z;    
    sout[3].w -= epsbar*s3[3].w;
    sout[3].w -= s1[3].w;   
    
    sout[4].x -= epsbar*s3[4].x;
    sout[4].x -= s1[4].x;
    sout[4].y -= epsbar*s3[4].y;
    sout[4].y -= s1[4].y;    
    sout[4].z -= epsbar*s3[4].z;
    sout[4].z -= s1[4].z;    
    sout[4].w -= epsbar*s3[4].w;
    sout[4].w -= s1[4].w;       
    
    sout[5].x -= epsbar*s3[5].x;
    sout[5].x -= s1[5].x;
    sout[5].y -= epsbar*s3[5].y;
    sout[5].y -= s1[5].y;    
    sout[5].z -= epsbar*s3[5].z;
    sout[5].z -= s1[5].z;    
    sout[5].w -= epsbar*s3[5].w;
    sout[5].w -= s1[5].w;       
    
    #ifdef RELATIVISTIC_BASIS
      dev_Gamma5_assign_rel(&(s1[0]), &(sout[0]));
    #else
      dev_Gamma5_assign(&(s1[0]), &(sout[0]));    
    #endif   
    dev_write_spinor(&(s1[0]),&(s2_up[pos]));    
    
    
    //lower output spinor
    dev_read_spinor(&(sout[0]), &(s2_dn[pos]));
    dev_read_spinor(&(s1[0]), &(s1_dn[pos]));    
    dev_read_spinor(&(s3[0]), &(s3_dn[pos]));    
    

     sout[0].x -= epsbar*s3[0].x;
    sout[0].x -= s1[0].x;
    sout[0].y -= epsbar*s3[0].y;
    sout[0].y -= s1[0].y;    
    sout[0].z -= epsbar*s3[0].z;
    sout[0].z -= s1[0].z;    
    sout[0].w -= epsbar*s3[0].w;
    sout[0].w -= s1[0].w;       
    
    sout[1].x -= epsbar*s3[1].x;
    sout[1].x -= s1[1].x;
    sout[1].y -= epsbar*s3[1].y;
    sout[1].y -= s1[1].y;    
    sout[1].z -= epsbar*s3[1].z;
    sout[1].z -= s1[1].z;    
    sout[1].w -= epsbar*s3[1].w;
    sout[1].w -= s1[1].w;       
    
    sout[2].x -= epsbar*s3[2].x;
    sout[2].x -= s1[2].x;
    sout[2].y -= epsbar*s3[2].y;
    sout[2].y -= s1[2].y;    
    sout[2].z -= epsbar*s3[2].z;
    sout[2].z -= s1[2].z;    
    sout[2].w -= epsbar*s3[2].w;
    sout[2].w -= s1[2].w;       
    
    sout[3].x -= epsbar*s3[3].x;
    sout[3].x -= s1[3].x;
    sout[3].y -= epsbar*s3[3].y;
    sout[3].y -= s1[3].y;    
    sout[3].z -= epsbar*s3[3].z;
    sout[3].z -= s1[3].z;    
    sout[3].w -= epsbar*s3[3].w;
    sout[3].w -= s1[3].w;   
    
    sout[4].x -= epsbar*s3[4].x;
    sout[4].x -= s1[4].x;
    sout[4].y -= epsbar*s3[4].y;
    sout[4].y -= s1[4].y;    
    sout[4].z -= epsbar*s3[4].z;
    sout[4].z -= s1[4].z;    
    sout[4].w -= epsbar*s3[4].w;
    sout[4].w -= s1[4].w;       
    
    sout[5].x -= epsbar*s3[5].x;
    sout[5].x -= s1[5].x;
    sout[5].y -= epsbar*s3[5].y;
    sout[5].y -= s1[5].y;    
    sout[5].z -= epsbar*s3[5].z;
    sout[5].z -= s1[5].z;    
    sout[5].w -= epsbar*s3[5].w;
    sout[5].w -= s1[5].w;         
    
    #ifdef RELATIVISTIC_BASIS
      dev_Gamma5_assign_rel(&(s1[0]), &(sout[0]));
    #else
      dev_Gamma5_assign(&(s1[0]), &(sout[0]));    
    #endif    
    dev_write_spinor(&(s1[0]),&(s2_dn[pos]));       
  }
}



//sout_up = gamma5*((1+imubar gamma5) s2_dn - epsbar s2_up - s1_up)
//sout_dn = gamma5*((1+imubar gamma5) s2_up - epsbar s2_dn - s1_dn)
__global__ void dev_nd_linalg3_gamma5 (dev_spinor * s1_up, dev_spinor * s1_dn,
				dev_spinor * s2_up, dev_spinor * s2_dn,  
				dev_spinor * sout_up, dev_spinor * sout_dn, float epsbar, float sign
				) {
  dev_spinor s1[6], s3[6], sout[6], slocal[6];	
  int pos = threadIdx.x + blockDim.x*blockIdx.x;
  dev_complex pm_imu = dev_initcomplex(0.0, sign * mubar);
  
  if (pos < dev_VOLUME) {
    //upper output spinor

    dev_read_spinor(&(s1[0]), &(s1_up[pos]));    

    #ifdef USETEXTURE
      dev_read_spinor_tex_up(&(s3[0]), pos);  
      dev_read_spinor_tex_dn(&(slocal[0]), pos);      
    #else
      dev_read_spinor(&(s3[0]), &(s2_up[pos]));  
      dev_read_spinor(&(slocal[0]), &(s2_dn[pos])); 
    #endif
   
    
    #ifdef RELATIVISTIC_BASIS
      dev_skalarmult_gamma5_spinor_rel(&(sout[0]), pm_imu, &(slocal[0]) );
    #else
      dev_skalarmult_gamma5_spinor(&(sout[0]), pm_imu, &(slocal[0]) );
    #endif
    dev_add_spinor_assign(&(sout[0]), &(slocal[0]));
    
 
    
    sout[0].x -= epsbar*s3[0].x;
    sout[0].x -= s1[0].x;
    sout[0].y -= epsbar*s3[0].y;
    sout[0].y -= s1[0].y;    
    sout[0].z -= epsbar*s3[0].z;
    sout[0].z -= s1[0].z;    
    sout[0].w -= epsbar*s3[0].w;
    sout[0].w -= s1[0].w;       
    
    sout[1].x -= epsbar*s3[1].x;
    sout[1].x -= s1[1].x;
    sout[1].y -= epsbar*s3[1].y;
    sout[1].y -= s1[1].y;    
    sout[1].z -= epsbar*s3[1].z;
    sout[1].z -= s1[1].z;    
    sout[1].w -= epsbar*s3[1].w;
    sout[1].w -= s1[1].w;       
    
    sout[2].x -= epsbar*s3[2].x;
    sout[2].x -= s1[2].x;
    sout[2].y -= epsbar*s3[2].y;
    sout[2].y -= s1[2].y;    
    sout[2].z -= epsbar*s3[2].z;
    sout[2].z -= s1[2].z;    
    sout[2].w -= epsbar*s3[2].w;
    sout[2].w -= s1[2].w;       
    
    sout[3].x -= epsbar*s3[3].x;
    sout[3].x -= s1[3].x;
    sout[3].y -= epsbar*s3[3].y;
    sout[3].y -= s1[3].y;    
    sout[3].z -= epsbar*s3[3].z;
    sout[3].z -= s1[3].z;    
    sout[3].w -= epsbar*s3[3].w;
    sout[3].w -= s1[3].w;   
    
    sout[4].x -= epsbar*s3[4].x;
    sout[4].x -= s1[4].x;
    sout[4].y -= epsbar*s3[4].y;
    sout[4].y -= s1[4].y;    
    sout[4].z -= epsbar*s3[4].z;
    sout[4].z -= s1[4].z;    
    sout[4].w -= epsbar*s3[4].w;
    sout[4].w -= s1[4].w;       
    
    sout[5].x -= epsbar*s3[5].x;
    sout[5].x -= s1[5].x;
    sout[5].y -= epsbar*s3[5].y;
    sout[5].y -= s1[5].y;    
    sout[5].z -= epsbar*s3[5].z;
    sout[5].z -= s1[5].z;    
    sout[5].w -= epsbar*s3[5].w;
    sout[5].w -= s1[5].w;       
    
    #ifdef RELATIVISTIC_BASIS
      dev_Gamma5_assign_rel(&(s1[0]), &(sout[0]));
    #else
      dev_Gamma5_assign(&(s1[0]), &(sout[0]));    
    #endif   
    dev_write_spinor(&(s1[0]),&(sout_up[pos]));    
    
    
    //lower output spinor
    
    pm_imu = dev_initcomplex(0.0, -sign * mubar);

    dev_read_spinor(&(s1[0]), &(s1_dn[pos]));    
    #ifdef USETEXTURE
      dev_read_spinor_tex_dn(&(s3[0]), pos);  
      dev_read_spinor_tex_up(&(slocal[0]), pos);      
    #else
      dev_read_spinor(&(s3[0]), &(s2_dn[pos]));  
      dev_read_spinor(&(slocal[0]), &(s2_up[pos]));    
    #endif
    //restore s2_up/dn from local variables
    //dev_copy_spinor_local(&(save_dn[0]), &(s3[0]));
    //dev_copy_spinor_local(&(save_up[0]), &(slocal[0])); 
    
    #ifdef RELATIVISTIC_BASIS
      dev_skalarmult_gamma5_spinor_rel(&(sout[0]), pm_imu, &(slocal[0]) );
    #else
      dev_skalarmult_gamma5_spinor(&(sout[0]), pm_imu, &(slocal[0]) );
    #endif
    dev_add_spinor_assign(&(sout[0]), &(slocal[0]));
    

    

    sout[0].x -= epsbar*s3[0].x;
    sout[0].x -= s1[0].x;
    sout[0].y -= epsbar*s3[0].y;
    sout[0].y -= s1[0].y;    
    sout[0].z -= epsbar*s3[0].z;
    sout[0].z -= s1[0].z;    
    sout[0].w -= epsbar*s3[0].w;
    sout[0].w -= s1[0].w;       
    
    sout[1].x -= epsbar*s3[1].x;
    sout[1].x -= s1[1].x;
    sout[1].y -= epsbar*s3[1].y;
    sout[1].y -= s1[1].y;    
    sout[1].z -= epsbar*s3[1].z;
    sout[1].z -= s1[1].z;    
    sout[1].w -= epsbar*s3[1].w;
    sout[1].w -= s1[1].w;       
    
    sout[2].x -= epsbar*s3[2].x;
    sout[2].x -= s1[2].x;
    sout[2].y -= epsbar*s3[2].y;
    sout[2].y -= s1[2].y;    
    sout[2].z -= epsbar*s3[2].z;
    sout[2].z -= s1[2].z;    
    sout[2].w -= epsbar*s3[2].w;
    sout[2].w -= s1[2].w;       
    
    sout[3].x -= epsbar*s3[3].x;
    sout[3].x -= s1[3].x;
    sout[3].y -= epsbar*s3[3].y;
    sout[3].y -= s1[3].y;    
    sout[3].z -= epsbar*s3[3].z;
    sout[3].z -= s1[3].z;    
    sout[3].w -= epsbar*s3[3].w;
    sout[3].w -= s1[3].w;   
    
    sout[4].x -= epsbar*s3[4].x;
    sout[4].x -= s1[4].x;
    sout[4].y -= epsbar*s3[4].y;
    sout[4].y -= s1[4].y;    
    sout[4].z -= epsbar*s3[4].z;
    sout[4].z -= s1[4].z;    
    sout[4].w -= epsbar*s3[4].w;
    sout[4].w -= s1[4].w;       
    
    sout[5].x -= epsbar*s3[5].x;
    sout[5].x -= s1[5].x;
    sout[5].y -= epsbar*s3[5].y;
    sout[5].y -= s1[5].y;    
    sout[5].z -= epsbar*s3[5].z;
    sout[5].z -= s1[5].z;    
    sout[5].w -= epsbar*s3[5].w;
    sout[5].w -= s1[5].w;         
    
    #ifdef RELATIVISTIC_BASIS
      dev_Gamma5_assign_rel(&(s1[0]), &(sout[0]));
    #else
      dev_Gamma5_assign(&(s1[0]), &(sout[0]));    
    #endif    
    dev_write_spinor(&(s1[0]),&(sout_dn[pos]));   
    
  }
}






//sout_up = gamma5*((1+imubar gamma5) s2_up - epsbar s2_dn - s1_up)
//sout_dn = gamma5*((1+imubar gamma5) s2_dn - epsbar s2_up - s1_dn)
__global__ void dev_nd_linalg4_gamma5 (dev_spinor * s1_up, dev_spinor * s1_dn,
				dev_spinor * s2_up, dev_spinor * s2_dn,  
				dev_spinor * sout_up, dev_spinor * sout_dn, float epsbar, float sign
				) {
  dev_spinor s1[6], s3[6], sout[6], slocal[6];	
  int pos = threadIdx.x + blockDim.x*blockIdx.x;
  dev_complex pm_imu = dev_initcomplex(0.0, sign * mubar);
  
  if (pos < dev_VOLUME) {
    //upper output spinor

    dev_read_spinor(&(s1[0]), &(s1_up[pos]));    
    #ifdef USETEXTURE
      dev_read_spinor_tex_dn(&(s3[0]), pos);  
      dev_read_spinor_tex_up(&(slocal[0]), pos);      
    #else
      dev_read_spinor(&(slocal[0]), &(s2_up[pos]));    
      dev_read_spinor(&(s3[0]), &(s2_dn[pos]));  
    #endif       
    
    #ifdef RELATIVISTIC_BASIS
      dev_skalarmult_gamma5_spinor_rel(&(sout[0]), pm_imu, &(slocal[0]) );
    #else
      dev_skalarmult_gamma5_spinor(&(sout[0]), pm_imu, &(slocal[0]) );
    #endif
    dev_add_spinor_assign(&(sout[0]), &(slocal[0]));
    
 
    
    sout[0].x -= epsbar*s3[0].x;
    sout[0].x -= s1[0].x;
    sout[0].y -= epsbar*s3[0].y;
    sout[0].y -= s1[0].y;    
    sout[0].z -= epsbar*s3[0].z;
    sout[0].z -= s1[0].z;    
    sout[0].w -= epsbar*s3[0].w;
    sout[0].w -= s1[0].w;       
    
    sout[1].x -= epsbar*s3[1].x;
    sout[1].x -= s1[1].x;
    sout[1].y -= epsbar*s3[1].y;
    sout[1].y -= s1[1].y;    
    sout[1].z -= epsbar*s3[1].z;
    sout[1].z -= s1[1].z;    
    sout[1].w -= epsbar*s3[1].w;
    sout[1].w -= s1[1].w;       
    
    sout[2].x -= epsbar*s3[2].x;
    sout[2].x -= s1[2].x;
    sout[2].y -= epsbar*s3[2].y;
    sout[2].y -= s1[2].y;    
    sout[2].z -= epsbar*s3[2].z;
    sout[2].z -= s1[2].z;    
    sout[2].w -= epsbar*s3[2].w;
    sout[2].w -= s1[2].w;       
    
    sout[3].x -= epsbar*s3[3].x;
    sout[3].x -= s1[3].x;
    sout[3].y -= epsbar*s3[3].y;
    sout[3].y -= s1[3].y;    
    sout[3].z -= epsbar*s3[3].z;
    sout[3].z -= s1[3].z;    
    sout[3].w -= epsbar*s3[3].w;
    sout[3].w -= s1[3].w;   
    
    sout[4].x -= epsbar*s3[4].x;
    sout[4].x -= s1[4].x;
    sout[4].y -= epsbar*s3[4].y;
    sout[4].y -= s1[4].y;    
    sout[4].z -= epsbar*s3[4].z;
    sout[4].z -= s1[4].z;    
    sout[4].w -= epsbar*s3[4].w;
    sout[4].w -= s1[4].w;       
    
    sout[5].x -= epsbar*s3[5].x;
    sout[5].x -= s1[5].x;
    sout[5].y -= epsbar*s3[5].y;
    sout[5].y -= s1[5].y;    
    sout[5].z -= epsbar*s3[5].z;
    sout[5].z -= s1[5].z;    
    sout[5].w -= epsbar*s3[5].w;
    sout[5].w -= s1[5].w;       
    
    #ifdef RELATIVISTIC_BASIS
      dev_Gamma5_assign_rel(&(s1[0]), &(sout[0]));
    #else
      dev_Gamma5_assign(&(s1[0]), &(sout[0]));    
    #endif   
    dev_write_spinor(&(s1[0]),&(sout_up[pos]));    
    
    
    //lower output spinor
    
    pm_imu = dev_initcomplex(0.0, -sign * mubar);

    dev_read_spinor(&(s1[0]), &(s1_dn[pos]));
    #ifdef USETEXTURE
      dev_read_spinor_tex_dn(&(slocal[0]), pos);     
      dev_read_spinor_tex_up(&(s3[0]), pos);  
    #else    
      dev_read_spinor(&(slocal[0]), &(s2_dn[pos]));    
      dev_read_spinor(&(s3[0]), &(s2_up[pos]));  
    #endif  
    
    
    #ifdef RELATIVISTIC_BASIS
      dev_skalarmult_gamma5_spinor_rel(&(sout[0]), pm_imu, &(slocal[0]) );
    #else
      dev_skalarmult_gamma5_spinor(&(sout[0]), pm_imu, &(slocal[0]) );
    #endif
    dev_add_spinor_assign(&(sout[0]), &(slocal[0]));
    

    

    sout[0].x -= epsbar*s3[0].x;
    sout[0].x -= s1[0].x;
    sout[0].y -= epsbar*s3[0].y;
    sout[0].y -= s1[0].y;    
    sout[0].z -= epsbar*s3[0].z;
    sout[0].z -= s1[0].z;    
    sout[0].w -= epsbar*s3[0].w;
    sout[0].w -= s1[0].w;       
    
    sout[1].x -= epsbar*s3[1].x;
    sout[1].x -= s1[1].x;
    sout[1].y -= epsbar*s3[1].y;
    sout[1].y -= s1[1].y;    
    sout[1].z -= epsbar*s3[1].z;
    sout[1].z -= s1[1].z;    
    sout[1].w -= epsbar*s3[1].w;
    sout[1].w -= s1[1].w;       
    
    sout[2].x -= epsbar*s3[2].x;
    sout[2].x -= s1[2].x;
    sout[2].y -= epsbar*s3[2].y;
    sout[2].y -= s1[2].y;    
    sout[2].z -= epsbar*s3[2].z;
    sout[2].z -= s1[2].z;    
    sout[2].w -= epsbar*s3[2].w;
    sout[2].w -= s1[2].w;       
    
    sout[3].x -= epsbar*s3[3].x;
    sout[3].x -= s1[3].x;
    sout[3].y -= epsbar*s3[3].y;
    sout[3].y -= s1[3].y;    
    sout[3].z -= epsbar*s3[3].z;
    sout[3].z -= s1[3].z;    
    sout[3].w -= epsbar*s3[3].w;
    sout[3].w -= s1[3].w;   
    
    sout[4].x -= epsbar*s3[4].x;
    sout[4].x -= s1[4].x;
    sout[4].y -= epsbar*s3[4].y;
    sout[4].y -= s1[4].y;    
    sout[4].z -= epsbar*s3[4].z;
    sout[4].z -= s1[4].z;    
    sout[4].w -= epsbar*s3[4].w;
    sout[4].w -= s1[4].w;       
    
    sout[5].x -= epsbar*s3[5].x;
    sout[5].x -= s1[5].x;
    sout[5].y -= epsbar*s3[5].y;
    sout[5].y -= s1[5].y;    
    sout[5].z -= epsbar*s3[5].z;
    sout[5].z -= s1[5].z;    
    sout[5].w -= epsbar*s3[5].w;
    sout[5].w -= s1[5].w;         
    
    #ifdef RELATIVISTIC_BASIS
      dev_Gamma5_assign_rel(&(s1[0]), &(sout[0]));
    #else
      dev_Gamma5_assign(&(s1[0]), &(sout[0]));    
    #endif    
    dev_write_spinor(&(s1[0]),&(sout_dn[pos]));   
    
  }
}







//s2_up = nrm*(s2_up +epsbar s1_dn)
//s2_dn = nrm*(s2_dn +epsbar s1_up)
// __global__ void dev_nd_linalg2 (dev_spinor * s1_up, dev_spinor * s1_dn, 
// 				dev_spinor * s2_up, dev_spinor * s2_dn,	float epsbar, float nrm			
// 				) {
//   dev_spinor s1[6], sout[6];	
//   int pos = threadIdx.x + blockDim.x*blockIdx.x;
//   
//   if (pos < dev_VOLUME) {
//     //upper output spinor
//     dev_read_spinor(&(sout[0]), &(s2_up[pos]));
//     #ifdef USETEXTURE
//       dev_read_spinor_tex_dn(&(s1[0]), pos);      
//     #else
//       dev_read_spinor(&(s1[0]), &(s1_dn[pos]));    
//     #endif
//     
//     sout[0].x += epsbar*s1[0].x;
//     sout[0].x *= nrm;
//     sout[0].y += epsbar*s1[0].y;
//     sout[0].y *= nrm;       
//     sout[0].z += epsbar*s1[0].z;
//     sout[0].z *= nrm;      
//     sout[0].w += epsbar*s1[0].w;
//     sout[0].w *= nrm;
//        
//     sout[1].x += epsbar*s1[1].x;
//     sout[1].x *= nrm;
//     sout[1].y += epsbar*s1[1].y;
//     sout[1].y *= nrm;       
//     sout[1].z += epsbar*s1[1].z;
//     sout[1].z *= nrm;      
//     sout[1].w += epsbar*s1[1].w;
//     sout[1].w *= nrm;    
//     
//     sout[2].x += epsbar*s1[2].x;
//     sout[2].x *= nrm;
//     sout[2].y += epsbar*s1[2].y;
//     sout[2].y *= nrm;       
//     sout[2].z += epsbar*s1[2].z;
//     sout[2].z *= nrm;      
//     sout[2].w += epsbar*s1[2].w;
//     sout[2].w *= nrm;    
// 
//     sout[3].x += epsbar*s1[3].x;
//     sout[3].x *= nrm;
//     sout[3].y += epsbar*s1[3].y;
//     sout[3].y *= nrm;       
//     sout[3].z += epsbar*s1[3].z;
//     sout[3].z *= nrm;      
//     sout[3].w += epsbar*s1[3].w;
//     sout[3].w *= nrm;    
// 
//     sout[4].x += epsbar*s1[4].x;
//     sout[4].x *= nrm;
//     sout[4].y += epsbar*s1[4].y;
//     sout[4].y *= nrm;       
//     sout[4].z += epsbar*s1[4].z;
//     sout[4].z *= nrm;      
//     sout[4].w += epsbar*s1[4].w;
//     sout[4].w *= nrm;    
//     
//     sout[5].x += epsbar*s1[5].x;
//     sout[5].x *= nrm;
//     sout[5].y += epsbar*s1[5].y;
//     sout[5].y *= nrm;       
//     sout[5].z += epsbar*s1[5].z;
//     sout[5].z *= nrm;      
//     sout[5].w += epsbar*s1[5].w;
//     sout[5].w *= nrm;    
//     
//     dev_write_spinor(&(sout[0]),&(s2_up[pos]));  
//     
//     
//     //upper output spinor
//     dev_read_spinor(&(sout[0]), &(s2_dn[pos]));
//     #ifdef USETEXTURE
//       dev_read_spinor_tex_up(&(s1[0]), pos);      
//     #else
//      dev_read_spinor(&(s1[0]), &(s1_up[pos]));    
//     #endif
//     
//     sout[0].x += epsbar*s1[0].x;
//     sout[0].x *= nrm;
//     sout[0].y += epsbar*s1[0].y;
//     sout[0].y *= nrm;       
//     sout[0].z += epsbar*s1[0].z;
//     sout[0].z *= nrm;      
//     sout[0].w += epsbar*s1[0].w;
//     sout[0].w *= nrm;
//        
//     sout[1].x += epsbar*s1[1].x;
//     sout[1].x *= nrm;
//     sout[1].y += epsbar*s1[1].y;
//     sout[1].y *= nrm;       
//     sout[1].z += epsbar*s1[1].z;
//     sout[1].z *= nrm;      
//     sout[1].w += epsbar*s1[1].w;
//     sout[1].w *= nrm;    
//     
//     sout[2].x += epsbar*s1[2].x;
//     sout[2].x *= nrm;
//     sout[2].y += epsbar*s1[2].y;
//     sout[2].y *= nrm;       
//     sout[2].z += epsbar*s1[2].z;
//     sout[2].z *= nrm;      
//     sout[2].w += epsbar*s1[2].w;
//     sout[2].w *= nrm;    
// 
//     sout[3].x += epsbar*s1[3].x;
//     sout[3].x *= nrm;
//     sout[3].y += epsbar*s1[3].y;
//     sout[3].y *= nrm;       
//     sout[3].z += epsbar*s1[3].z;
//     sout[3].z *= nrm;      
//     sout[3].w += epsbar*s1[3].w;
//     sout[3].w *= nrm;    
// 
//     sout[4].x += epsbar*s1[4].x;
//     sout[4].x *= nrm;
//     sout[4].y += epsbar*s1[4].y;
//     sout[4].y *= nrm;       
//     sout[4].z += epsbar*s1[4].z;
//     sout[4].z *= nrm;      
//     sout[4].w += epsbar*s1[4].w;
//     sout[4].w *= nrm;    
//     
//     sout[5].x += epsbar*s1[5].x;
//     sout[5].x *= nrm;
//     sout[5].y += epsbar*s1[5].y;
//     sout[5].y *= nrm;       
//     sout[5].z += epsbar*s1[5].z;
//     sout[5].z *= nrm;      
//     sout[5].w += epsbar*s1[5].w;
//     sout[5].w *= nrm;    
//     
// /*    s2_dn[pos+0*DEVOFF].x = sout[0].x;
//     s2_dn[pos+0*DEVOFF].y = sout[0].y;
//     s2_dn[pos+0*DEVOFF].z = sout[0].z;
//     s2_dn[pos+0*DEVOFF].w = sout[0].w;   
//     
//     s2_dn[pos+1*DEVOFF].x = sout[1].x;
//     s2_dn[pos+1*DEVOFF].y = sout[1].y;
//     s2_dn[pos+1*DEVOFF].z = sout[1].z;
//     s2_dn[pos+1*DEVOFF].w = sout[1].w;       
//     
//     s2_dn[pos+2*DEVOFF].x = sout[2].x;
//     s2_dn[pos+2*DEVOFF].y = sout[2].y;
//     s2_dn[pos+2*DEVOFF].z = sout[2].z;
//     s2_dn[pos+2*DEVOFF].w = sout[2].w;      
//     
//     s2_dn[pos+3*DEVOFF].x = sout[3].x;
//     s2_dn[pos+3*DEVOFF].y = sout[3].y;
//     s2_dn[pos+3*DEVOFF].z = sout[3].z;
//     s2_dn[pos+3*DEVOFF].w = sout[3].w;      
//     
//     s2_dn[pos+4*DEVOFF].x = sout[4].x;
//     s2_dn[pos+4*DEVOFF].y = sout[4].y;
//     s2_dn[pos+4*DEVOFF].z = sout[4].z;
//     s2_dn[pos+4*DEVOFF].w = sout[4].w;    
//     
//     s2_dn[pos+5*DEVOFF].x = sout[5].x;
//     s2_dn[pos+5*DEVOFF].y = sout[5].y;
//     s2_dn[pos+5*DEVOFF].z = sout[5].z;
//     s2_dn[pos+5*DEVOFF].w = sout[5].w;    */  
//     
//     dev_write_spinor(&(sout[0]),&(s2_dn[pos]));     
//     
//     
//     
//   }
// }
// 


__global__ void dev_nd_linalg2 (dev_spinor * s1_up, dev_spinor * s1_dn, 
				dev_spinor * s2_up, dev_spinor * s2_dn,	float epsbar, float nrm			
				) {
  dev_spinor s1[6], sout[6];	
  int pos = threadIdx.x + blockDim.x*blockIdx.x;
  
  if (pos < dev_VOLUME) {
    //upper output spinor
    dev_read_spinor(&(sout[0]), &(s2_up[pos]));
    #ifdef USETEXTURE
      dev_read_spinor_tex_dn(&(s1[0]), pos);      
    #else
      dev_read_spinor(&(s1[0]), &(s1_dn[pos]));    
    #endif
    
    sout[0].x += epsbar*s1[0].x;
    s2_up[pos+0*DEVOFF].x = sout[0].x*nrm;    
    sout[0].y += epsbar*s1[0].y;
    s2_up[pos+0*DEVOFF].y = sout[0].y*nrm;    
    sout[0].z += epsbar*s1[0].z;
    s2_up[pos+0*DEVOFF].z = sout[0].z*nrm;      
    sout[0].w += epsbar*s1[0].w;
    s2_up[pos+0*DEVOFF].w = sout[0].w*nrm; 
       
    sout[1].x += epsbar*s1[1].x;
    s2_up[pos+1*DEVOFF].x = sout[1].x*nrm; 
    sout[1].y += epsbar*s1[1].y;
    s2_up[pos+1*DEVOFF].y = sout[1].y*nrm;      
    sout[1].z += epsbar*s1[1].z;
    s2_up[pos+1*DEVOFF].z = sout[1].z*nrm;       
    sout[1].w += epsbar*s1[1].w;
    s2_up[pos+1*DEVOFF].w = sout[1].w*nrm;    
    
    sout[2].x += epsbar*s1[2].x;
    s2_up[pos+2*DEVOFF].x = sout[2].x*nrm; 
    sout[2].y += epsbar*s1[2].y;
    s2_up[pos+2*DEVOFF].y = sout[2].y*nrm;  
    sout[2].z += epsbar*s1[2].z;
    s2_up[pos+2*DEVOFF].z = sout[2].z*nrm;    
    sout[2].w += epsbar*s1[2].w;
    s2_up[pos+2*DEVOFF].w = sout[2].w*nrm;   

    sout[3].x += epsbar*s1[3].x;
    s2_up[pos+3*DEVOFF].x = sout[3].x*nrm;
    sout[3].y += epsbar*s1[3].y;
    s2_up[pos+3*DEVOFF].y = sout[3].y*nrm;   
    sout[3].z += epsbar*s1[3].z;
    s2_up[pos+3*DEVOFF].z = sout[3].z*nrm;     
    sout[3].w += epsbar*s1[3].w;
    s2_up[pos+3*DEVOFF].w = sout[3].w*nrm;   

    sout[4].x += epsbar*s1[4].x;
    s2_up[pos+4*DEVOFF].x = sout[4].x*nrm;
    sout[4].y += epsbar*s1[4].y;
    s2_up[pos+4*DEVOFF].y = sout[4].y*nrm;      
    sout[4].z += epsbar*s1[4].z;
    s2_up[pos+4*DEVOFF].z = sout[4].z*nrm;    
    sout[4].w += epsbar*s1[4].w;
    s2_up[pos+4*DEVOFF].w = sout[4].w*nrm;   
    
    sout[5].x += epsbar*s1[5].x;
    s2_up[pos+5*DEVOFF].x = sout[5].x*nrm;
    sout[5].y += epsbar*s1[5].y;
    s2_up[pos+5*DEVOFF].y = sout[5].y*nrm;      
    sout[5].z += epsbar*s1[5].z;
    s2_up[pos+5*DEVOFF].z = sout[5].z*nrm;   
    sout[5].w += epsbar*s1[5].w;
    s2_up[pos+5*DEVOFF].w = sout[5].w*nrm;   
    

    
    
    //upper output spinor
    dev_read_spinor(&(sout[0]), &(s2_dn[pos]));
    #ifdef USETEXTURE
      dev_read_spinor_tex_up(&(s1[0]), pos);      
    #else
     dev_read_spinor(&(s1[0]), &(s1_up[pos]));    
    #endif
    
    sout[0].x += epsbar*s1[0].x;
    s2_dn[pos+0*DEVOFF].x = sout[0].x*nrm;
    sout[0].y += epsbar*s1[0].y;
    s2_dn[pos+0*DEVOFF].y = sout[0].y*nrm;        
    sout[0].z += epsbar*s1[0].z;
    s2_dn[pos+0*DEVOFF].z = sout[0].z*nrm;             
    sout[0].w += epsbar*s1[0].w;
    s2_dn[pos+0*DEVOFF].w = sout[0].w*nrm;  
       
    sout[1].x += epsbar*s1[1].x;
    s2_dn[pos+1*DEVOFF].x = sout[1].x*nrm;    
    sout[1].y += epsbar*s1[1].y;
    s2_dn[pos+1*DEVOFF].y = sout[1].y*nrm; 
    sout[1].z += epsbar*s1[1].z;
    s2_dn[pos+1*DEVOFF].z = sout[1].z*nrm; 
    sout[1].w += epsbar*s1[1].w;
    s2_dn[pos+1*DEVOFF].w = sout[1].w*nrm; 
    
    sout[2].x += epsbar*s1[2].x;
    s2_dn[pos+2*DEVOFF].x = sout[2].x*nrm; 
    sout[2].y += epsbar*s1[2].y;
    s2_dn[pos+2*DEVOFF].y = sout[2].y*nrm;     
    sout[2].z += epsbar*s1[2].z;
    s2_dn[pos+2*DEVOFF].z = sout[2].z*nrm; 
    sout[2].w += epsbar*s1[2].w;
    s2_dn[pos+2*DEVOFF].w = sout[2].w*nrm;  

    sout[3].x += epsbar*s1[3].x;
    s2_dn[pos+3*DEVOFF].x = sout[3].x*nrm; 
    sout[3].y += epsbar*s1[3].y;
    s2_dn[pos+3*DEVOFF].y = sout[3].y*nrm;  
    sout[3].z += epsbar*s1[3].z;
    s2_dn[pos+3*DEVOFF].z = sout[3].z*nrm;      
    sout[3].w += epsbar*s1[3].w;
    s2_dn[pos+3*DEVOFF].w = sout[3].w*nrm;     

    sout[4].x += epsbar*s1[4].x;
    s2_dn[pos+4*DEVOFF].x = sout[4].x*nrm; 
    sout[4].y += epsbar*s1[4].y;
    s2_dn[pos+4*DEVOFF].y = sout[4].y*nrm;    
    sout[4].z += epsbar*s1[4].z;
    s2_dn[pos+4*DEVOFF].z = sout[4].z*nrm;      
    sout[4].w += epsbar*s1[4].w;
    s2_dn[pos+4*DEVOFF].w = sout[4].w*nrm;    
    
    sout[5].x += epsbar*s1[5].x;
    s2_dn[pos+5*DEVOFF].x = sout[5].x*nrm; 
    sout[5].y += epsbar*s1[5].y;
    s2_dn[pos+5*DEVOFF].y = sout[5].y*nrm; 
    sout[5].z += epsbar*s1[5].z;
    s2_dn[pos+5*DEVOFF].z = sout[5].z*nrm;    
    sout[5].w += epsbar*s1[5].w;
    s2_dn[pos+5*DEVOFF].w = sout[5].w*nrm;    

  }
}










///////////////////////////
// MATRIX MULTIPLICATION //
///////////////////////////

// the GPU implementation of  Qtm_pm_ndpsi(...)  from tm_operators_nd.c
// we have the possibility to add a constant shift > 0 
void dev_Qtm_pm_ndpsi_old (dev_spinor * spinout_up, dev_spinor * spinout_dn,
                              dev_spinor * spinin_up , dev_spinor * spinin_dn , 
                              int gridsize1, int blocksize1, int gridsize2, int blocksize2,
                              int gridsize3, int blocksize3, int gridsize4, int blocksize4) {
  
  
  // we will use the auxiliary fields  dev_spin_eo{1,2}_up/dn  for working on and buffering
  // and set  dev_spin_eo2_up/dn  equal  spinout_up/dn
  // spinin_up/dn  have to remain unchanged !!
  // spinout_up/dn  can be freely used
  
  /////////////////////
  // LOCAL VARIABLES //
  /////////////////////
  
  int N_sites  =    VOLUME/2;		// #lattice sites
  
  
  
  
  ////////////////////////////////////
  //   MATCHING with  Qtm_pm_ndpsi  //
  ////////////////////////////////////
  //                                //
  // _strange = _up                 //
  // _charm   = _dn                 //
  //                                //
  // DUM_MATRIX   = dev_spin_eo1_up //
  // DUM_MATRIX+1 = dev_spin_eo1_dn //
  //                                //
  // DUM_MATRIX+2 = dev_spin_eo2_up //
  // DUM_MATRIX+3 = dev_spin_eo2_dn //
  //                                //
  ////////////////////////////////////
  
  
  
  
  ///////////////////////////////////
  // INITIALIZATIONS & ASSIGNMENTS //	// have to use (one) other auxiliary field(s) than the calling function dev_cg_eo_nd
  ///////////////////////////////////
  
  dev_spin_eo2_up = spinout_up;		// need no memory allocated
  dev_spin_eo2_dn = spinout_dn;
  												///////////// THEORY ////////////////////////////////////////////////////////////////
  												//                                                                                 //
  												//  (Q_tilde) = gamma5 * ((M_oo) - (M_oe)(Mee^-1)(M_eo))                           //
  												//  (Q_tilde)(Q_tilde_dagger) * (up,dn) = (Q_tilde) * (b,a)                        //
  ///////////////										//                                                    (a,b) = (Q_tilde) * (dn,up)  //
  // MAIN BODY //										//                                                                                 //
  ///////////////										/////////////////////////////////////////////////////////////////////////////////////
  
  
  double nrm = 1.0 / (1.0 + g_mubar*g_mubar - g_epsbar*g_epsbar);
  
  
  ///////////////////////////////////////							/////////////////////////////////
  //        Q_tilde_dagger(2x2)        //							// (a,b) = (Q_tilde) * (dn,up) //
  ///////////////////////////////////////							/////////////////////////////////
  
  
  ////////////////////////////
  // (M_oe) (Mee^-1) (M_eo) //
  ////////////////////////////
  
  // Hopping:
  #ifdef USETEXTURE
    bind_texture_spin(spinin_dn,1);
  #endif
    dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, spinin_dn, dev_spin_eo1_up, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0,0,N_sites);	// dev_spin_eo1_up = (M_eo) * spinin_dn
  #ifdef USETEXTURE
    unbind_texture_spin(1);
  #endif
  
  #ifdef USETEXTURE
    bind_texture_spin(spinin_up,1);
  #endif
    dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, spinin_up, dev_spin_eo1_dn, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0,0,N_sites);	// dev_spin_eo1_dn = (M_eo) * spinin_up
  #ifdef USETEXTURE
    unbind_texture_spin(1);
  #endif

  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo1_up, dev_spin_eo2_up, -1.0);		// dev_spin_eo2_up  =  (1 - imubar) * dev_spin_eo1_up  =  (1 - imubar)*(M_eo) * spinin_dn
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo1_dn, dev_spin_eo2_dn, +1.0);		// dev_spin_eo2_dn  =  (1 + imubar) * dev_spin_eo1_dn  =  (1 + imubar)*(M_eo) * spinin_up

  dev_nd_linalg2<<<gridsize2, blocksize2>>>(dev_spin_eo1_up, dev_spin_eo1_dn, dev_spin_eo2_up, dev_spin_eo2_dn, g_epsbar, nrm); 


  // Hopping:
  #ifdef USETEXTURE
    bind_texture_spin(dev_spin_eo2_up,1);									// remember: this is ((M_oe)(Mee^-1)(M_eo)) * (spinin_dn, spinin_up):
  #endif
    dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo2_up, dev_spin_eo1_up, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1,0,N_sites);
  #ifdef USETEXTURE
    unbind_texture_spin(1);							// dev_spin_eo1_up  =  (M_oe) * dev_spin_eo2_up  =  (M_oe)*nrm*(1-imubar)*(M_eo) * spinin_dn  +  (M_oe)*nrm*epsbar*(M_eo) * spinin_up
  #endif									
  
  #ifdef USETEXTURE
    bind_texture_spin(dev_spin_eo2_dn,1);
  #endif
    dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo2_dn, dev_spin_eo1_dn, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1,0,N_sites);
  #ifdef USETEXTURE
    unbind_texture_spin(1);							// dev_spin_eo1_dn  =  (M_oe) * dev_spin_eo2_dn  =  (M_oe)*nrm*(1+imubar)*(M_eo) * spinin_up  +  (M_oe)*nrm*epsbar*(M_eo) * spinin_dn
  #endif
  
 
  ////////////
  // (M_oo) //
  ////////////
  
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(spinin_dn, dev_spin_eo2_up, +1.0);			// dev_spin_eo2_up = (1 + imubar) * spinin_dn
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(spinin_up, dev_spin_eo2_dn, -1.0);			// dev_spin_eo2_dn = (1 - imubar) * spinin_up
  
  ///////////////////////////////////////
  // (M_oo)  -  (M_oe) (Mee^-1) (M_eo) //
  ///////////////////////////////////////
   
  dev_nd_linalg1_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo1_up, dev_spin_eo1_dn, spinin_up, spinin_dn, dev_spin_eo2_up, dev_spin_eo2_dn, g_epsbar ); 
    

  ////////////////////
  // (a,b) -> (b,a) //
  ////////////////////
  dev_copy_spinor_field<<<gridsize4, blocksize4>>>(dev_spin_eo2_dn, dev_spin_eo3_up);			// dev_spin_eo3_up = dev_spin_eo2_dn
  dev_copy_spinor_field<<<gridsize4, blocksize4>>>(dev_spin_eo2_up, dev_spin_eo3_dn);			// dev_spin_eo3_dn = dev_spin_eo2_up

  ///////////////////////////////////								///////////////////////
  //        Q_tilde(2x2)           //								// (Q_tilde) * (b,a) //
  ///////////////////////////////////								///////////////////////
  
  ////////////////////////////
  // (M_oe) (Mee^-1) (M_eo) //
  ////////////////////////////
  
  // Hopping:
  #ifdef USETEXTURE
    bind_texture_spin(dev_spin_eo3_up,1);
  #endif
    dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo3_up, dev_spin_eo1_up, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0,0,N_sites);	// dev_spin_eo1_up = (M_eo) * dev_spin_eo3_up
  #ifdef USETEXTURE
    unbind_texture_spin(1);
  #endif
  
  #ifdef USETEXTURE
    bind_texture_spin(dev_spin_eo3_dn,1);
  #endif
    dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo3_dn, dev_spin_eo1_dn, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0,0,N_sites);	// dev_spin_eo1_dn = (M_eo) * dev_spin_eo3_dn
  #ifdef USETEXTURE
    unbind_texture_spin(1);
  #endif
  
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo1_up, dev_spin_eo2_up, -1.0);			// dev_spin_eo2_up  =  (1 - imubar) * dev_spin_eo1_up  =  (1 - imubar)*(M_eo) * dev_spin_eo3_up
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo1_dn, dev_spin_eo2_dn, +1.0);			// dev_spin_eo2_dn  =  (1 + imubar) * dev_spin_eo1_dn  =  (1 + imubar)*(M_eo) * dev_spin_eo3_dn
  dev_nd_linalg2<<<gridsize2, blocksize2>>>(dev_spin_eo1_up, dev_spin_eo1_dn, dev_spin_eo2_up, dev_spin_eo2_dn, g_epsbar, nrm); 
  
  // Hopping:
  #ifdef USETEXTURE
    bind_texture_spin(dev_spin_eo2_up,1);									// remember: this is ((M_oe) (Mee^-1) (M_eo)) * (dev_spin_eo3_up, dev_spin_eo3_dn):
  #endif
    dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo2_up, dev_spin_eo1_up, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1,0,N_sites);
  #ifdef USETEXTURE
    unbind_texture_spin(1);							// dev_spin_eo1_up  =  (M_oe) * dev_spin_eo2_up  =  (M_oe)*nrm*(1-imubar)*(M_eo) * dev_spin_eo3_up  +  (M_oe)*nrm*epsbar*(M_eo) * dev_spin_eo3_dn
  #endif
  
  #ifdef USETEXTURE
    bind_texture_spin(dev_spin_eo2_dn,1);
  #endif
    dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo2_dn, dev_spin_eo1_dn, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1,0,N_sites);
  #ifdef USETEXTURE
    unbind_texture_spin(1);							// dev_spin_eo1_dn  =  (M_oe) * dev_spin_eo2_dn  =  (M_oe)*nrm*(1+imubar)*(M_eo) * dev_spin_eo3_dn  +  (M_oe)*nrm*epsbar*(M_eo) * dev_spin_eo3_up
  #endif
  
  ////////////
  // (M_oo) //
  ////////////
  
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo3_up, dev_spin_eo2_up, +1.0);				// dev_spin_eo2_up = (1 + imubar) * dev_spin_eo3_up
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo3_dn, dev_spin_eo2_dn, -1.0);				// dev_spin_eo2_dn = (1 - imubar) * dev_spin_eo3_dn
  
  //dev_spin_eo3_up <-> dev_spin_eo3_dn here!! 
  dev_nd_linalg1_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo1_up, dev_spin_eo1_dn, dev_spin_eo3_dn, dev_spin_eo3_up, dev_spin_eo2_up, dev_spin_eo2_dn, g_epsbar ); 
  
  
  ////////////
  // output //	output is already done by setting  dev_spin_eo2_up/dn = spinout_up/dn
  ////////////
 
  return;
  
}//dev_Qtm_pm_ndpsi_old()









void dev_Qtm_pm_ndpsi_scalar (dev_spinor * spinout_up, dev_spinor * spinout_dn,
                              dev_spinor * spinin_up , dev_spinor * spinin_dn , 
                              int gridsize1, int blocksize1, int gridsize2, int blocksize2,
                              int gridsize3, int blocksize3, int gridsize4, int blocksize4) {
  
  
  // we will use the auxiliary fields  dev_spin_eo{1,2}_up/dn  for working on and buffering
  // and set  dev_spin_eo2_up/dn  equal  spinout_up/dn
  // spinin_up/dn  have to remain unchanged !!
  // spinout_up/dn  can be freely used
  
  /////////////////////
  // LOCAL VARIABLES //
  /////////////////////
  
  int N_sites  =    VOLUME/2;		// #lattice sites
  int N_floats = 24*VOLUME/2;		// #floats
  
  
  
  
  ////////////////////////////////////
  //   MATCHING with  Qtm_pm_ndpsi  //
  ////////////////////////////////////
  //                                //
  // _strange = _up                 //
  // _charm   = _dn                 //
  //                                //
  // DUM_MATRIX   = dev_spin_eo1_up //
  // DUM_MATRIX+1 = dev_spin_eo1_dn //
  //                                //
  // DUM_MATRIX+2 = dev_spin_eo2_up //
  // DUM_MATRIX+3 = dev_spin_eo2_dn //
  //                                //
  ////////////////////////////////////
  
  
  
  
  ///////////////////////////////////
  // INITIALIZATIONS & ASSIGNMENTS //	// have to use (one) other auxiliary field(s) than the calling function dev_cg_eo_nd
  ///////////////////////////////////
  
  dev_spin_eo2_up = spinout_up;		// need no memory allocated
  dev_spin_eo2_dn = spinout_dn;
  												///////////// THEORY ////////////////////////////////////////////////////////////////
  												//                                                                                 //
  												//  (Q_tilde) = gamma5 * ((M_oo) - (M_oe)(Mee^-1)(M_eo))                           //
  												//  (Q_tilde)(Q_tilde_dagger) * (up,dn) = (Q_tilde) * (b,a)                        //
  ///////////////										//                                                    (a,b) = (Q_tilde) * (dn,up)  //
  // MAIN BODY //										//                                                                                 //
  ///////////////										/////////////////////////////////////////////////////////////////////////////////////
  
  
  double nrm = 1.0 / (1.0 + g_mubar*g_mubar - g_epsbar*g_epsbar);
  
  
  ///////////////////////////////////////							/////////////////////////////////
  //        Q_tilde_dagger(2x2)        //							// (a,b) = (Q_tilde) * (dn,up) //
  ///////////////////////////////////////							/////////////////////////////////
  
  
  ////////////////////////////
  // (M_oe) (Mee^-1) (M_eo) //
  ////////////////////////////
  
  // Hopping:
  #ifdef USETEXTURE
    bind_texture_spin(spinin_dn,1);
  #endif
    dev_Hopping_Matrix_ext<<<gridsize1, blocksize1, 0, stream_nd[0]>>>(dev_gf, spinin_dn, dev_spin_eo1_up, dev_spin_eo2_up, -1.0, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0,0,N_sites);	// dev_spin_eo1_up = (M_eo) * spinin_dn
  #ifdef USETEXTURE
    unbind_texture_spin(1);
  #endif
  
  #ifdef USETEXTURE
    bind_texture_spin(spinin_up,1);
  #endif
    dev_Hopping_Matrix_ext<<<gridsize1, blocksize1, 0, stream_nd[1]>>>(dev_gf, spinin_up, dev_spin_eo1_dn, dev_spin_eo2_dn, 1.0, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0,0,N_sites);	// dev_spin_eo1_dn = (M_eo) * spinin_up
  #ifdef USETEXTURE
    unbind_texture_spin(1);
  #endif

  cudaStreamSynchronize(stream_nd[0]);
  cudaStreamSynchronize(stream_nd[1]);  

  #ifdef USETEXTURE
    bind_texture_spin(dev_spin_eo1_up,1);
    bind_texture_spin_dn(dev_spin_eo1_dn,1);   
  #endif  
  dev_nd_linalg2<<<gridsize2, blocksize2>>>(dev_spin_eo1_up, dev_spin_eo1_dn, dev_spin_eo2_up, dev_spin_eo2_dn, g_epsbar, nrm); 
  #ifdef USETEXTURE
    unbind_texture_spin(1);
    unbind_texture_spin_dn(1);    
  #endif
  

  // Hopping:
  #ifdef USETEXTURE
    bind_texture_spin(dev_spin_eo2_up,1);									// remember: this is ((M_oe)(Mee^-1)(M_eo)) * (spinin_dn, spinin_up):
  #endif
    dev_Hopping_Matrix<<<gridsize1, blocksize1, 0, stream_nd[0]>>>(dev_gf, dev_spin_eo2_up, dev_spin_eo1_up, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1,0,N_sites);
  #ifdef USETEXTURE
    unbind_texture_spin(1);							// dev_spin_eo1_up  =  (M_oe) * dev_spin_eo2_up  =  (M_oe)*nrm*(1-imubar)*(M_eo) * spinin_dn  +  (M_oe)*nrm*epsbar*(M_eo) * spinin_up
  #endif									
  
  #ifdef USETEXTURE
    bind_texture_spin(dev_spin_eo2_dn,1);
  #endif
    dev_Hopping_Matrix<<<gridsize1, blocksize1, 0, stream_nd[1]>>>(dev_gf, dev_spin_eo2_dn, dev_spin_eo1_dn, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1,0,N_sites);
  #ifdef USETEXTURE
    unbind_texture_spin(1);							// dev_spin_eo1_dn  =  (M_oe) * dev_spin_eo2_dn  =  (M_oe)*nrm*(1+imubar)*(M_eo) * spinin_up  +  (M_oe)*nrm*epsbar*(M_eo) * spinin_dn
  #endif
  
 
  ////////////
  // (M_oo) //
  ////////////
  
  ///////////////////////////////////////
  // (M_oo)  -  (M_oe) (Mee^-1) (M_eo) //
  ///////////////////////////////////////  
  
  cudaStreamSynchronize(stream_nd[0]);
  cudaStreamSynchronize(stream_nd[1]);  
  
  #ifdef USETEXTURE
    bind_texture_spin(spinin_up,1);
    bind_texture_spin_dn(spinin_dn,1);   
  #endif     
  dev_nd_linalg3_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo1_up, dev_spin_eo1_dn, spinin_up, spinin_dn, dev_spin_eo2_up, dev_spin_eo2_dn, g_epsbar, +1.0 ); 
  #ifdef USETEXTURE
    unbind_texture_spin(1);
    unbind_texture_spin_dn(1);    
  #endif 
 
  

  ////////////////////
  // (a,b) -> (b,a) //
  ////////////////////
  dev_copy_spinor_field<<<gridsize4, blocksize4, 0, stream_nd[0]>>>(dev_spin_eo2_dn, dev_spin_eo3_up);			// dev_spin_eo3_up = dev_spin_eo2_dn
  dev_copy_spinor_field<<<gridsize4, blocksize4, 0, stream_nd[1]>>>(dev_spin_eo2_up, dev_spin_eo3_dn);			// dev_spin_eo3_dn = dev_spin_eo2_up

  ///////////////////////////////////								///////////////////////
  //        Q_tilde(2x2)           //								// (Q_tilde) * (b,a) //
  ///////////////////////////////////								///////////////////////
  
  ////////////////////////////
  // (M_oe) (Mee^-1) (M_eo) //
  ////////////////////////////
  
  // Hopping:
  #ifdef USETEXTURE
    bind_texture_spin(dev_spin_eo3_up,1);
  #endif
    dev_Hopping_Matrix_ext<<<gridsize1, blocksize1, 0, stream_nd[0]>>>(dev_gf, dev_spin_eo3_up, dev_spin_eo1_up, dev_spin_eo2_up, -1.0, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0,0,N_sites);	// dev_spin_eo1_up = (M_eo) * dev_spin_eo3_up
  #ifdef USETEXTURE
    unbind_texture_spin(1);
  #endif
  
  #ifdef USETEXTURE
    bind_texture_spin(dev_spin_eo3_dn,1);
  #endif
    dev_Hopping_Matrix_ext<<<gridsize1, blocksize1, 0, stream_nd[1]>>>(dev_gf, dev_spin_eo3_dn, dev_spin_eo1_dn, dev_spin_eo2_dn, +1.0, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0,0,N_sites);	// dev_spin_eo1_dn = (M_eo) * dev_spin_eo3_dn
  #ifdef USETEXTURE
    unbind_texture_spin(1);
  #endif
  
  cudaStreamSynchronize(stream_nd[0]);
  cudaStreamSynchronize(stream_nd[1]);  
  
  #ifdef USETEXTURE
    bind_texture_spin(dev_spin_eo1_up,1);
    bind_texture_spin_dn(dev_spin_eo1_dn,1);   
  #endif    
  dev_nd_linalg2<<<gridsize2, blocksize2>>>(dev_spin_eo1_up, dev_spin_eo1_dn, dev_spin_eo2_up, dev_spin_eo2_dn, g_epsbar, nrm); 
  #ifdef USETEXTURE
    unbind_texture_spin(1);
    unbind_texture_spin_dn(1);    
  #endif 
  
  // Hopping:
  #ifdef USETEXTURE
    bind_texture_spin(dev_spin_eo2_up,1);									// remember: this is ((M_oe) (Mee^-1) (M_eo)) * (dev_spin_eo3_up, dev_spin_eo3_dn):
  #endif
    dev_Hopping_Matrix<<<gridsize1, blocksize1, 0, stream_nd[0]>>>(dev_gf, dev_spin_eo2_up, dev_spin_eo1_up, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1,0,N_sites);
  #ifdef USETEXTURE
    unbind_texture_spin(1);							// dev_spin_eo1_up  =  (M_oe) * dev_spin_eo2_up  =  (M_oe)*nrm*(1-imubar)*(M_eo) * dev_spin_eo3_up  +  (M_oe)*nrm*epsbar*(M_eo) * dev_spin_eo3_dn
  #endif
  
  #ifdef USETEXTURE
    bind_texture_spin(dev_spin_eo2_dn,1);
  #endif
    dev_Hopping_Matrix<<<gridsize1, blocksize1, 0, stream_nd[1]>>>(dev_gf, dev_spin_eo2_dn, dev_spin_eo1_dn, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1,0,N_sites);
  #ifdef USETEXTURE
    unbind_texture_spin(1);							// dev_spin_eo1_dn  =  (M_oe) * dev_spin_eo2_dn  =  (M_oe)*nrm*(1+imubar)*(M_eo) * dev_spin_eo3_dn  +  (M_oe)*nrm*epsbar*(M_eo) * dev_spin_eo3_up
  #endif
  
  ////////////
  // (M_oo) //
  ////////////
  cudaStreamSynchronize(stream_nd[0]);
  cudaStreamSynchronize(stream_nd[1]);   

  #ifdef USETEXTURE
    bind_texture_spin(dev_spin_eo3_up,1);
    bind_texture_spin_dn(dev_spin_eo3_dn,1);   
  #endif    
    dev_nd_linalg4_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo1_up, dev_spin_eo1_dn, dev_spin_eo3_up, dev_spin_eo3_dn, dev_spin_eo2_up, dev_spin_eo2_dn, g_epsbar, +1.0 ); 
  #ifdef USETEXTURE
    unbind_texture_spin(1);
    unbind_texture_spin_dn(1);    
  #endif 
  
  cublasSscal (N_floats, phmc_invmaxev*phmc_invmaxev, (float *) dev_spin_eo2_up, 1);
  cublasSscal (N_floats, phmc_invmaxev*phmc_invmaxev, (float *) dev_spin_eo2_dn, 1);
  
  ////////////
  // output //	output is already done by setting  dev_spin_eo2_up/dn = spinout_up/dn
  ////////////
 
  return;
  
}//dev_Qtm_pm_ndpsi()









///////////////////////////
// MATRIX MULTIPLICATION //
///////////////////////////

// the GPU implementation of  Qtm_pm_ndpsi(...)  from tm_operators_nd.c
// we have the possibility to add a constant shift > 0 
void dev_Qtm_pm_ndpsi_updn (dev_spinor * spinout_up, dev_spinor * spinout_dn,
                              dev_spinor * spinin_up , dev_spinor * spinin_dn , 
                              int gridsize1, int blocksize1, int gridsize2, int blocksize2,
                              int gridsize3, int blocksize3, int gridsize4, int blocksize4) {
  
  
  // we will use the auxiliary fields  dev_spin_eo{1,2}_up/dn  for working on and buffering
  // and set  dev_spin_eo2_up/dn  equal  spinout_up/dn
  // spinin_up/dn  have to remain unchanged !!
  // spinout_up/dn  can be freely used
  

  int N_sites  =    VOLUME/2;		// #lattice sites
  int N_floats = 24*VOLUME/2;		// #floats

  dev_spin_eo2_up = spinout_up;		// need no memory allocated
  dev_spin_eo2_dn = spinout_dn;

 
  double nrm = 1.0 / (1.0 + g_mubar*g_mubar - g_epsbar*g_epsbar);
  
  #ifdef USETEXTURE
    bind_texture_spin(spinin_dn,1);
    bind_texture_spin_dn(spinin_up,1);    
  #endif
    dev_Hopping_Matrix_updn<<<gridsize1, blocksize1>>>(dev_gf, spinin_dn, spinin_up, dev_spin_eo1_up, dev_spin_eo1_dn,dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0,0,N_sites);	// dev_spin_eo1_up = (M_eo) * spinin_dn
  #ifdef USETEXTURE
    unbind_texture_spin(1);
    unbind_texture_spin_dn(1);   
  #endif
  
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo1_up, dev_spin_eo2_up, -1.0);		// dev_spin_eo2_up  =  (1 - imubar) * dev_spin_eo1_up  =  (1 - imubar)*(M_eo) * spinin_dn
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo1_dn, dev_spin_eo2_dn, +1.0);		// dev_spin_eo2_dn  =  (1 + imubar) * dev_spin_eo1_dn  =  (1 + imubar)*(M_eo) * spinin_up

  /*
  cublasSaxpy (N_floats, g_epsbar, (float *) dev_spin_eo1_dn, 1, (float *) dev_spin_eo2_up, 1);
  										// dev_spin_eo2_up  =  dev_spin_eo2_up  +  epsbar * dev_spin_eo1_dn  =  (1 - imubar)*(M_eo) * spinin_dn  +  epsbar * (M_eo) * spinin_up
  cublasSaxpy (N_floats, g_epsbar, (float *) dev_spin_eo1_up, 1, (float *) dev_spin_eo2_dn, 1);
  										// dev_spin_eo2_dn  =  dev_spin_eo2_dn  +  epsbar * dev_spin_eo1_up  =  (1 + imubar)*(M_eo) * spinin_up  +  epsbar * (M_eo) * spinin_dn
  

  cublasSscal (N_floats, nrm, (float *) dev_spin_eo2_up, 1);			// dev_spin_eo2_up  =  nrm * dev_spin_eo2_up  =  nrm*(1-imubar)*(M_eo) * spinin_dn  +  nrm*epsbar*(M_eo) * spinin_up
  
  cublasSscal (N_floats, nrm, (float *) dev_spin_eo2_dn, 1);			// dev_spin_eo2_dn  =  nrm * dev_spin_eo2_dn  =  nrm*(1+imubar)*(M_eo) * spinin_up  +  nrm*epsbar*(M_eo) * spinin_dn
  */
  dev_nd_linalg2<<<gridsize2, blocksize2>>>(dev_spin_eo1_up, dev_spin_eo1_dn, dev_spin_eo2_up, dev_spin_eo2_dn, g_epsbar, nrm); 

  #ifdef USETEXTURE
    bind_texture_spin(dev_spin_eo2_up,1);									// remember: this is ((M_oe)(Mee^-1)(M_eo)) * (spinin_dn, spinin_up):
    bind_texture_spin_dn(dev_spin_eo2_dn,1);									// remember: this is ((M_oe)(Mee^-1)(M_eo)) * (spinin_dn, spinin_up):    
  #endif
    dev_Hopping_Matrix_updn<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo2_up, dev_spin_eo2_dn, dev_spin_eo1_up, dev_spin_eo1_dn, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1,0,N_sites);
  #ifdef USETEXTURE
    unbind_texture_spin(1);	
    unbind_texture_spin_dn(1);
  #endif									

  
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(spinin_dn, dev_spin_eo2_up, +1.0);			// dev_spin_eo2_up = (1 + imubar) * spinin_dn
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(spinin_up, dev_spin_eo2_dn, -1.0);			// dev_spin_eo2_dn = (1 - imubar) * spinin_up
  
  /*
  cublasSaxpy (N_floats, -g_epsbar, (float *) spinin_up, 1, (float *) dev_spin_eo2_up, 1);
  												// dev_spin_eo2_up  =  dev_spin_eo2_up - epsbar*spinin_up  =  (1+imubar)*spinin_dn - epsbar*spinin_up
  cublasSaxpy (N_floats, -g_epsbar, (float *) spinin_dn, 1, (float *) dev_spin_eo2_dn, 1);
  												// dev_spin_eo2_dn  =  dev_spin_eo2_dn - epsbar*spinin_dn  =  (1-imubar)*spinin_up - epsbar*spinin_dn
  
  													// this is ((M_oo) - (M_oe)(Mee^-1)(M_eo)) * (spinin_dn, spinin_up):
  cublasSaxpy (N_floats, -1.0, (float *) dev_spin_eo1_up, 1, (float *) dev_spin_eo2_up, 1);
  							// dev_spin_eo2_up  =  dev_spin_eo2_up  -  dev_spin_eo1_up  =                      (1+imubar) * spinin_dn  -                    epsbar * spinin_up
  							//                                                             - (M_oe)*nrm*(1-imubar)*(M_eo) * spinin_dn  -  (M_oe)*nrm*epsbar*(M_eo) * spinin_up
  cublasSaxpy (N_floats, -1.0, (float *) dev_spin_eo1_dn, 1, (float *) dev_spin_eo2_dn, 1);
  							// dev_spin_eo2_dn  =  dev_spin_eo2_dn  -  dev_spin_eo1_dn  =                      (1-imubar) * spinin_up  -                    epsbar * spinin_dn
  							//                                                             - (M_oe)*nrm*(1+imubar)*(M_eo) * spinin_up  -  (M_oe)*nrm*epsbar*(M_eo) * spinin_dn
  */
  
  dev_nd_linalg1<<<gridsize2, blocksize2>>>(dev_spin_eo1_up, dev_spin_eo1_dn, spinin_up, spinin_dn, dev_spin_eo2_up, dev_spin_eo2_dn, g_epsbar ); 
  
  
  #ifdef RELATIVISTIC_BASIS
    dev_gamma5_rel<<<gridsize3, blocksize3>>>(dev_spin_eo2_up, dev_spin_eo2_up);					// dev_spin_eo2_up = gamma5 * dev_spin_eo2_up 
    dev_gamma5_rel<<<gridsize3, blocksize3>>>(dev_spin_eo2_dn, dev_spin_eo2_dn);					// dev_spin_eo2_dn = gamma5 * dev_spin_eo2_dn  
  #else
    dev_gamma5<<<gridsize3, blocksize3>>>(dev_spin_eo2_up, dev_spin_eo2_up);					// dev_spin_eo2_up = gamma5 * dev_spin_eo2_up 
    dev_gamma5<<<gridsize3, blocksize3>>>(dev_spin_eo2_dn, dev_spin_eo2_dn);					// dev_spin_eo2_dn = gamma5 * dev_spin_eo2_dn
  #endif
  
  
  
  ////////////////////
  // (a,b) -> (b,a) //
  ////////////////////
  dev_copy_spinor_field<<<gridsize4, blocksize4>>>(dev_spin_eo2_dn, dev_spin_eo3_up);			// dev_spin_eo3_up = dev_spin_eo2_dn
  dev_copy_spinor_field<<<gridsize4, blocksize4>>>(dev_spin_eo2_up, dev_spin_eo3_dn);			// dev_spin_eo3_dn = dev_spin_eo2_up
  
  
  
  #ifdef USETEXTURE
    bind_texture_spin(dev_spin_eo3_up,1);
    bind_texture_spin_dn(dev_spin_eo3_dn,1);   
  #endif
    dev_Hopping_Matrix_updn<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo3_up, dev_spin_eo3_dn, dev_spin_eo1_up, dev_spin_eo1_dn, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0,0,N_sites);	// dev_spin_eo1_up = (M_eo) * dev_spin_eo3_up
  #ifdef USETEXTURE
    unbind_texture_spin(1);
    unbind_texture_spin_dn(1);    
  #endif
  
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo1_up, dev_spin_eo2_up, -1.0);			// dev_spin_eo2_up  =  (1 - imubar) * dev_spin_eo1_up  =  (1 - imubar)*(M_eo) * dev_spin_eo3_up
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo1_dn, dev_spin_eo2_dn, +1.0);			// dev_spin_eo2_dn  =  (1 + imubar) * dev_spin_eo1_dn  =  (1 + imubar)*(M_eo) * dev_spin_eo3_dn

  /*
  cublasSaxpy (N_floats, g_epsbar, (float *) dev_spin_eo1_dn, 1, (float *) dev_spin_eo2_up, 1);
  										// dev_spin_eo2_up  =  dev_spin_eo2_up  +  epsbar * dev_spin_eo1_dn  =  (1 - imubar)*(M_eo) * dev_spin_eo3_up  +  epsbar * (M_eo) * dev_spin_eo3_dn
  cublasSaxpy (N_floats, g_epsbar, (float *) dev_spin_eo1_up, 1, (float *) dev_spin_eo2_dn, 1);
  										// dev_spin_eo2_dn  =  dev_spin_eo2_dn  +  epsbar * dev_spin_eo1_up  =  (1 + imubar)*(M_eo) * dev_spin_eo3_dn  +  epsbar * (M_eo) * dev_spin_eo3_up
  
  cublasSscal (N_floats, nrm, (float *) dev_spin_eo2_up, 1);			// dev_spin_eo2_up  =  nrm * dev_spin_eo2_up  =  nrm*(1-imubar)*(M_eo) * dev_spin_eo3_up  +  nrm*epsbar*(M_eo) * dev_spin_eo3_dn
  
  cublasSscal (N_floats, nrm, (float *) dev_spin_eo2_dn, 1);			// dev_spin_eo2_dn  =  nrm * dev_spin_eo2_dn  =  nrm*(1+imubar)*(M_eo) * dev_spin_eo3_dn  +  nrm*epsbar*(M_eo) * dev_spin_eo3_up
  */
  
  dev_nd_linalg2<<<gridsize2, blocksize2>>>(dev_spin_eo1_up, dev_spin_eo1_dn, dev_spin_eo2_up, dev_spin_eo2_dn, g_epsbar, nrm); 

  
  #ifdef USETEXTURE
    bind_texture_spin(dev_spin_eo2_up,1);	
    bind_texture_spin_dn(dev_spin_eo2_dn,1);    
  #endif
    dev_Hopping_Matrix_updn<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo2_up, dev_spin_eo2_dn, dev_spin_eo1_up,dev_spin_eo1_dn, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1,0,N_sites);
  #ifdef USETEXTURE
    unbind_texture_spin(1);
    unbind_texture_spin_dn(1);
  #endif
  

  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo3_up, dev_spin_eo2_up, +1.0);				// dev_spin_eo2_up = (1 + imubar) * dev_spin_eo3_up
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo3_dn, dev_spin_eo2_dn, -1.0);				// dev_spin_eo2_dn = (1 - imubar) * dev_spin_eo3_dn
  
  /*											// remember: this is (M_oo) * (dev_spin_eo3_up, dev_spin_eo3_dn):
  cublasSaxpy (N_floats, -g_epsbar, (float *) dev_spin_eo3_dn, 1, (float *) dev_spin_eo2_up, 1);
  												// dev_spin_eo2_up  =  dev_spin_eo2_up - epsbar*dev_spin_eo3_dn  =  (1+imubar)*dev_spin_eo3_up - epsbar*dev_spin_eo3_dn
  cublasSaxpy (N_floats, -g_epsbar, (float *) dev_spin_eo3_up, 1, (float *) dev_spin_eo2_dn, 1);
  												// dev_spin_eo2_dn  =  dev_spin_eo2_dn - epsbar*dev_spin_eo3_up  =  (1-imubar)*dev_spin_eo3_dn - epsbar*dev_spin_eo3_up
  
  

  											// this is ( (M_oo) - (M_oe) (Mee^-1) (M_eo) ) * (dev_spin_eo3_up, dev_spin_eo3_dn)
  cublasSaxpy (N_floats, -1.0, (float *) dev_spin_eo1_up, 1, (float *) dev_spin_eo2_up, 1);
  							// dev_spin_eo2_up  =  dev_spin_eo2_up  -  dev_spin_eo1_up  =                      (1+imubar) * dev_spin_eo3_up  -                    epsbar * dev_spin_eo3_dn
  							//                                                             - (M_oe)*nrm*(1-imubar)*(M_eo) * dev_spin_eo3_up  -  (M_oe)*nrm*epsbar*(M_eo) * dev_spin_eo3_dn
  cublasSaxpy (N_floats, -1.0, (float *) dev_spin_eo1_dn, 1, (float *) dev_spin_eo2_dn, 1);
  							// dev_spin_eo2_dn  =  dev_spin_eo2_dn  -  dev_spin_eo1_dn  =                      (1-imubar) * dev_spin_eo3_dn  -                    epsbar * dev_spin_eo3_up
  							//                                                             - (M_oe)*nrm*(1+imubar)*(M_eo) * dev_spin_eo3_dn  -  (M_oe)*nrm*epsbar*(M_eo) * dev_spin_eo3_up
  */
  dev_nd_linalg1<<<gridsize2, blocksize2>>>(dev_spin_eo1_up, dev_spin_eo1_dn, dev_spin_eo3_dn, dev_spin_eo3_up, dev_spin_eo2_up, dev_spin_eo2_dn, g_epsbar ); 

  
  #ifdef RELATIVISTIC_BASIS
    dev_gamma5_rel<<<gridsize3, blocksize3>>>(dev_spin_eo2_up, dev_spin_eo2_up);					// dev_spin_eo2_up = gamma5 * dev_spin_eo2_up 
    dev_gamma5_rel<<<gridsize3, blocksize3>>>(dev_spin_eo2_dn, dev_spin_eo2_dn);					// dev_spin_eo2_dn = gamma5 * dev_spin_eo2_dn  
  #else
    dev_gamma5<<<gridsize3, blocksize3>>>(dev_spin_eo2_up, dev_spin_eo2_up);					// dev_spin_eo2_up = gamma5 * dev_spin_eo2_up 
    dev_gamma5<<<gridsize3, blocksize3>>>(dev_spin_eo2_dn, dev_spin_eo2_dn);					// dev_spin_eo2_dn = gamma5 * dev_spin_eo2_dn
  #endif
  
  cublasSscal (N_floats, phmc_invmaxev*phmc_invmaxev, (float *) dev_spin_eo2_up, 1);
  cublasSscal (N_floats, phmc_invmaxev*phmc_invmaxev, (float *) dev_spin_eo2_dn, 1);
  
  return;
  
}//dev_Qtm_pm_ndpsi_updn()




#ifdef _USE_MPI

void dev_Qtm_pm_ndpsi_mpi (dev_spinor * spinout_up, dev_spinor * spinout_dn,
                                        dev_spinor * spinin_up , dev_spinor * spinin_dn ,
                                        int gridsize1, int blocksize1, int gridsize2, int blocksize2,
                                        int gridsize3, int blocksize3, int gridsize4, int blocksize4) {
  
  
  // we will use the auxiliary fields  dev_spin_eo{1,2}_up/dn  for working on and buffering
  // and set  dev_spin_eo2_up/dn  equal  spinout_up/dn
  // spinin_up/dn  have to remain unchanged !!
  // spinout_up/dn  can be freely used
  
  
  
  
  /////////////////////
  // LOCAL VARIABLES //
  /////////////////////
  
  int N_floats = 24*VOLUME/2;		// #floats
  
    
  
  ////////////////////////////////////
  //   MATCHING with Q_Qdagger_ND   //
  ////////////////////////////////////
  //                                //
  // _strange = _up                 //
  // _charm   = _dn                 //
  //                                //
  // DUM_MATRIX   = dev_spin_eo1_up //
  // DUM_MATRIX+1 = dev_spin_eo1_dn //
  //                                //
  // DUM_MATRIX+2 = dev_spin_eo2_up //
  // DUM_MATRIX+3 = dev_spin_eo2_dn //
  //                                //
  ////////////////////////////////////
  
  
  
  
  ///////////////////////////////////
  // INITIALIZATIONS & ASSIGNMENTS //	// have to use (one) other auxiliary field(s) than the calling function dev_cg_eo_nd
  ///////////////////////////////////
  
  dev_spin_eo2_up = spinout_up;		// need no memory allocated
  dev_spin_eo2_dn = spinout_dn;
  												///////////// THEORY ////////////////////////////////////////////////////////////////
  												//                                                                                 //
  												//  (Q_tilde) = gamma5 * ((M_oo) - (M_oe)(Mee^-1)(M_eo))                           //
  												//  (Q_tilde)(Q_tilde_dagger) * (up,dn) = (Q_tilde) * (b,a)                        //
  ///////////////										//                                                    (a,b) = (Q_tilde) * (dn,up)  //
  // MAIN BODY //										//                                                                                 //
  ///////////////										/////////////////////////////////////////////////////////////////////////////////////
  
  
  double nrm = 1.0 / (1.0 + g_mubar*g_mubar - g_epsbar*g_epsbar);
  
  //printf("This is matrix_multiplication32_mpi_ASYNC().\n");
  
  
  ///////////////////////////////////////							/////////////////////////////////
  //        Q_tilde_dagger(2x2)        //							// (a,b) = (Q_tilde) * (dn,up) //
  ///////////////////////////////////////							/////////////////////////////////
  
  
  ////////////////////////////
  // (M_oe) (Mee^-1) (M_eo) //
  ////////////////////////////
  
#if ASYNC == 0

  int N_sites  =    VOLUME/2;		// #lattice sites

  ////////////////////////////
  // (M_oe) (Mee^-1) (M_eo) //
  ////////////////////////////
  
  
  		// xchange
  		#ifndef HOPPING_DEBUG
  		  xchange_field_wrapper(spinin_dn, 0);
  		#endif
  
  // hopping matrix
  #ifdef USETEXTURE
    bind_texture_spin(spinin_dn,1);
  #endif
  
    #ifndef HOPPING_DEBUG
      dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, spinin_dn, dev_spin_eo1_up, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0,0,N_sites);	// dev_spin_eo1_up = (M_eo) * spinin_dn
    #else
      Hopping_Matrix_wrapper(0, dev_spin_eo1_up, spinin_dn);
    #endif
  
  #ifdef USETEXTURE
    unbind_texture_spin(1);
  #endif
  
  
  		// xchange
  		#ifndef HOPPING_DEBUG
  		  xchange_field_wrapper(spinin_up, 0);
  		#endif
  
  // hopping matrix
  #ifdef USETEXTURE
    bind_texture_spin(spinin_up,1);
  #endif
  
    #ifndef HOPPING_DEBUG
      dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, spinin_up, dev_spin_eo1_dn, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0,0,N_sites);	// dev_spin_eo1_dn = (M_eo) * spinin_up
    #else
      Hopping_Matrix_wrapper(0, dev_spin_eo1_dn, spinin_up);
    #endif
  
  #ifdef USETEXTURE
    unbind_texture_spin(1);
  #endif
  
  
  // imubar, gamma5
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo1_up, dev_spin_eo2_up, -1.0);		// dev_spin_eo2_up  =  (1 - imubar) * dev_spin_eo1_up  =  (1 - imubar)*(M_eo) * spinin_dn
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo1_dn, dev_spin_eo2_dn, +1.0);		// dev_spin_eo2_dn  =  (1 + imubar) * dev_spin_eo1_dn  =  (1 + imubar)*(M_eo) * spinin_up
  
  
  // linear algebra
  cublasSaxpy_wrapper (N_floats, g_epsbar, (float *) dev_spin_eo1_dn, (float *) dev_spin_eo2_up);
  										// dev_spin_eo2_up  =  dev_spin_eo2_up  +  epsbar * dev_spin_eo1_dn  =  (1 - imubar)*(M_eo) * spinin_dn  +  epsbar * (M_eo) * spinin_up
  cublasSaxpy_wrapper (N_floats, g_epsbar, (float *) dev_spin_eo1_up, (float *) dev_spin_eo2_dn);
  										// dev_spin_eo2_dn  =  dev_spin_eo2_dn  +  epsbar * dev_spin_eo1_up  =  (1 + imubar)*(M_eo) * spinin_up  +  epsbar * (M_eo) * spinin_dn
  
  // linear algebra
  cublasSscal_wrapper (N_floats, nrm, (float *) dev_spin_eo2_up);			// dev_spin_eo2_up  =  nrm * dev_spin_eo2_up  =  nrm*(1-imubar)*(M_eo) * spinin_dn  +  nrm*epsbar*(M_eo) * spinin_up
  
  cublasSscal_wrapper (N_floats, nrm, (float *) dev_spin_eo2_dn);			// dev_spin_eo2_dn  =  nrm * dev_spin_eo2_dn  =  nrm*(1+imubar)*(M_eo) * spinin_up  +  nrm*epsbar*(M_eo) * spinin_dn
  
  
  		// xchange
  		#ifndef HOPPING_DEBUG
  		  xchange_field_wrapper(dev_spin_eo2_up, 1);
  		#endif
  
  // hopping matrix
  #ifdef USETEXTURE
    bind_texture_spin(dev_spin_eo2_up,1);									// remember: this is ((M_oe)(Mee^-1)(M_eo)) * (spinin_dn, spinin_up):
  #endif
  
    #ifndef HOPPING_DEBUG
      dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo2_up, dev_spin_eo1_up, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1,0,N_sites);
    #else
      Hopping_Matrix_wrapper(1, dev_spin_eo1_up, dev_spin_eo2_up);
    #endif
  
  #ifdef USETEXTURE
    unbind_texture_spin(1);							// dev_spin_eo1_up  =  (M_oe) * dev_spin_eo2_up  =  (M_oe)*nrm*(1-imubar)*(M_eo) * spinin_dn  +  (M_oe)*nrm*epsbar*(M_eo) * spinin_up
  #endif									
  
  
  		// xchange
  		#ifndef HOPPING_DEBUG
  		  xchange_field_wrapper(dev_spin_eo2_dn, 1);
  		#endif
  
  // hopping matrix
  #ifdef USETEXTURE
    bind_texture_spin(dev_spin_eo2_dn,1);
  #endif
  
    #ifndef HOPPING_DEBUG
      dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo2_dn, dev_spin_eo1_dn, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1,0,N_sites);
    #else
      Hopping_Matrix_wrapper(1, dev_spin_eo1_dn, dev_spin_eo2_dn);
    #endif
    
  #ifdef USETEXTURE
    unbind_texture_spin(1);							// dev_spin_eo1_dn  =  (M_oe) * dev_spin_eo2_dn  =  (M_oe)*nrm*(1+imubar)*(M_eo) * spinin_up  +  (M_oe)*nrm*epsbar*(M_eo) * spinin_dn
  #endif
  
  
  
  
  ////////////
  // (M_oo) //
  ////////////
  
  
  // imubar, gamma5
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(spinin_dn, dev_spin_eo2_up, +1.0);			// dev_spin_eo2_up = (1 + imubar) * spinin_dn
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(spinin_up, dev_spin_eo2_dn, -1.0);			// dev_spin_eo2_dn = (1 - imubar) * spinin_up
  
  
  // linear algebra													// remember: this is (M_oo) * (spinin_dn, spinin_up):
  cublasSaxpy_wrapper (N_floats, -g_epsbar, (float *) spinin_up, (float *) dev_spin_eo2_up);
  												// dev_spin_eo2_up  =  dev_spin_eo2_up - epsbar*spinin_up  =  (1+imubar)*spinin_dn - epsbar*spinin_up
  cublasSaxpy_wrapper (N_floats, -g_epsbar, (float *) spinin_dn, (float *) dev_spin_eo2_dn);
  												// dev_spin_eo2_dn  =  dev_spin_eo2_dn - epsbar*spinin_dn  =  (1-imubar)*spinin_up - epsbar*spinin_dn
  
  
  
  ///////////////////////////////////////
  // (M_oo)  -  (M_oe) (Mee^-1) (M_eo) //
  ///////////////////////////////////////
  
  // linear algebra													// this is ((M_oo) - (M_oe)(Mee^-1)(M_eo)) * (spinin_dn, spinin_up):
  cublasSaxpy_wrapper(N_floats, -1.0, (float *) dev_spin_eo1_up, (float *) dev_spin_eo2_up);
  							// dev_spin_eo2_up  =  dev_spin_eo2_up  -  dev_spin_eo1_up  =                      (1+imubar) * spinin_dn  -                    epsbar * spinin_up
  							//                                                             - (M_oe)*nrm*(1-imubar)*(M_eo) * spinin_dn  -  (M_oe)*nrm*epsbar*(M_eo) * spinin_up
  cublasSaxpy_wrapper (N_floats, -1.0, (float *) dev_spin_eo1_dn, (float *) dev_spin_eo2_dn);
  							// dev_spin_eo2_dn  =  dev_spin_eo2_dn  -  dev_spin_eo1_dn  =                      (1-imubar) * spinin_up  -                    epsbar * spinin_dn
  							//                                                             - (M_oe)*nrm*(1+imubar)*(M_eo) * spinin_up  -  (M_oe)*nrm*epsbar*(M_eo) * spinin_dn
  
  
  ////////////
  // gamma5 //
  ////////////
  
  // gamma5
  #ifdef RELATIVISTIC_BASIS
    dev_gamma5_rel<<<gridsize3, blocksize3>>>(dev_spin_eo2_up, dev_spin_eo2_up);					// dev_spin_eo2_up = gamma5 * dev_spin_eo2_up 
    dev_gamma5_rel<<<gridsize3, blocksize3>>>(dev_spin_eo2_dn, dev_spin_eo2_dn);					// dev_spin_eo2_dn = gamma5 * dev_spin_eo2_dn  
  #else
    dev_gamma5<<<gridsize3, blocksize3>>>(dev_spin_eo2_up, dev_spin_eo2_up);					// dev_spin_eo2_up = gamma5 * dev_spin_eo2_up 
    dev_gamma5<<<gridsize3, blocksize3>>>(dev_spin_eo2_dn, dev_spin_eo2_dn);					// dev_spin_eo2_dn = gamma5 * dev_spin_eo2_dn
  #endif	
  
  
  
  
  
  ////////////////////
  // (a,b) -> (b,a) //
  ////////////////////
  dev_copy_spinor_field<<<gridsize4, blocksize4>>>(dev_spin_eo2_dn, dev_spin_eo3_up);			// dev_spin_eo3_up = dev_spin_eo2_dn
  dev_copy_spinor_field<<<gridsize4, blocksize4>>>(dev_spin_eo2_up, dev_spin_eo3_dn);			// dev_spin_eo3_dn = dev_spin_eo2_up
  
  
  
  
  
  
  ///////////////////////////////////								///////////////////////
  //        Q_tilde(2x2)           //								// (Q_tilde) * (b,a) //
  ///////////////////////////////////								///////////////////////
  
  
  ////////////////////////////
  // (M_oe) (Mee^-1) (M_eo) //
  ////////////////////////////
  
  		// xchange
  		#ifndef HOPPING_DEBUG
  		  xchange_field_wrapper(dev_spin_eo3_up, 0);
  		#endif
  
  // hopping matrix
  #ifdef USETEXTURE
    bind_texture_spin(dev_spin_eo3_up,1);
  #endif
  
    #ifndef HOPPING_DEBUG
      dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo3_up, dev_spin_eo1_up, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0,0,N_sites);	// dev_spin_eo1_up = (M_eo) * dev_spin_eo3_up
    #else
      Hopping_Matrix_wrapper(0, dev_spin_eo1_up, dev_spin_eo3_up);
    #endif
  
  #ifdef USETEXTURE
    unbind_texture_spin(1);
  #endif
  
  
  		// xchange
  		#ifndef HOPPING_DEBUG
  		  xchange_field_wrapper(dev_spin_eo3_dn, 0);
  		#endif
  
  // hopping matrix
  #ifdef USETEXTURE
    bind_texture_spin(dev_spin_eo3_dn,1);
  #endif
  
    #ifndef HOPPING_DEBUG
      dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo3_dn, dev_spin_eo1_dn, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0,0,N_sites);	// dev_spin_eo1_dn = (M_eo) * dev_spin_eo3_dn
    #else
      Hopping_Matrix_wrapper(0, dev_spin_eo1_dn, dev_spin_eo3_dn);
    #endif
  
  #ifdef USETEXTURE
    unbind_texture_spin(1);
  #endif
  
  
  // imubar, gamma5
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo1_up, dev_spin_eo2_up, -1.0);			// dev_spin_eo2_up  =  (1 - imubar) * dev_spin_eo1_up  =  (1 - imubar)*(M_eo) * dev_spin_eo3_up
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo1_dn, dev_spin_eo2_dn, +1.0);			// dev_spin_eo2_dn  =  (1 + imubar) * dev_spin_eo1_dn  =  (1 + imubar)*(M_eo) * dev_spin_eo3_dn
  
  // linear algebra
  cublasSaxpy_wrapper (N_floats, g_epsbar, (float *) dev_spin_eo1_dn, (float *) dev_spin_eo2_up);										// dev_spin_eo2_up  =  dev_spin_eo2_up  +  epsbar * dev_spin_eo1_dn  =  (1 - imubar)*(M_eo) * dev_spin_eo3_up  +  epsbar * (M_eo) * dev_spin_eo3_dn
  cublasSaxpy_wrapper (N_floats, g_epsbar, (float *) dev_spin_eo1_up, (float *) dev_spin_eo2_dn);
  										// dev_spin_eo2_dn  =  dev_spin_eo2_dn  +  epsbar * dev_spin_eo1_up  =  (1 + imubar)*(M_eo) * dev_spin_eo3_dn  +  epsbar * (M_eo) * dev_spin_eo3_up
  // lineare algebra
  cublasSscal_wrapper (N_floats, nrm, (float *) dev_spin_eo2_up);			// dev_spin_eo2_up  =  nrm * dev_spin_eo2_up  =  nrm*(1-imubar)*(M_eo) * dev_spin_eo3_up  +  nrm*epsbar*(M_eo) * dev_spin_eo3_dn
  cublasSscal_wrapper (N_floats, nrm, (float *) dev_spin_eo2_dn);			// dev_spin_eo2_dn  =  nrm * dev_spin_eo2_dn  =  nrm*(1+imubar)*(M_eo) * dev_spin_eo3_dn  +  nrm*epsbar*(M_eo) * dev_spin_eo3_up
  
  		// xchange
  		#ifndef HOPPING_DEBUG
  		  xchange_field_wrapper(dev_spin_eo2_up, 1);
  		#endif
  
  // hopping matrix
  #ifdef USETEXTURE
    bind_texture_spin(dev_spin_eo2_up,1);									// remember: this is ((M_oe) (Mee^-1) (M_eo)) * (dev_spin_eo3_up, dev_spin_eo3_dn):
  #endif
  
    #ifndef HOPPING_DEBUG
      dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo2_up, dev_spin_eo1_up, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1,0,N_sites);
    #else
      Hopping_Matrix_wrapper(1, dev_spin_eo1_up, dev_spin_eo2_up);
    #endif
  
  #ifdef USETEXTURE
    unbind_texture_spin(1);							// dev_spin_eo1_up  =  (M_oe) * dev_spin_eo2_up  =  (M_oe)*nrm*(1-imubar)*(M_eo) * dev_spin_eo3_up  +  (M_oe)*nrm*epsbar*(M_eo) * dev_spin_eo3_dn
  #endif
  
  		// xchange
  		#ifndef HOPPING_DEBUG
  		  xchange_field_wrapper(dev_spin_eo2_dn, 1);
  		#endif
  
  // hopping matrix
  #ifdef USETEXTURE
    bind_texture_spin(dev_spin_eo2_dn,1);
  #endif
  
    #ifndef HOPPING_DEBUG
      dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo2_dn, dev_spin_eo1_dn, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1,0,N_sites);
    #else
      Hopping_Matrix_wrapper(1, dev_spin_eo1_dn, dev_spin_eo2_dn);
    #endif
  
  #ifdef USETEXTURE
    unbind_texture_spin(1);							// dev_spin_eo1_dn  =  (M_oe) * dev_spin_eo2_dn  =  (M_oe)*nrm*(1+imubar)*(M_eo) * dev_spin_eo3_dn  +  (M_oe)*nrm*epsbar*(M_eo) * dev_spin_eo3_up
  #endif
  
  
  ////////////
  // (M_oo) //
  ////////////
  
  // imubar, gamma5
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo3_up, dev_spin_eo2_up, +1.0);				// dev_spin_eo2_up = (1 + imubar) * dev_spin_eo3_up
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo3_dn, dev_spin_eo2_dn, -1.0);				// dev_spin_eo2_dn = (1 - imubar) * dev_spin_eo3_dn
  
  
  // lineare algebra										// remember: this is (M_oo) * (dev_spin_eo3_up, dev_spin_eo3_dn):
  cublasSaxpy_wrapper (N_floats, -g_epsbar, (float *) dev_spin_eo3_dn, (float *) dev_spin_eo2_up);												// dev_spin_eo2_up  =  dev_spin_eo2_up - epsbar*dev_spin_eo3_dn  =  (1+imubar)*dev_spin_eo3_up - epsbar*dev_spin_eo3_dn
  cublasSaxpy_wrapper (N_floats, -g_epsbar, (float *) dev_spin_eo3_up, (float *) dev_spin_eo2_dn);
  												// dev_spin_eo2_dn  =  dev_spin_eo2_dn - epsbar*dev_spin_eo3_up  =  (1-imubar)*dev_spin_eo3_dn - epsbar*dev_spin_eo3_up
  
  ///////////////////////////////////////
  // (M_oo)  -  (M_oe) (Mee^-1) (M_eo) //
  ///////////////////////////////////////
  
  // lineare algebra										// this is ( (M_oo) - (M_oe) (Mee^-1) (M_eo) ) * (dev_spin_eo3_up, dev_spin_eo3_dn)
  cublasSaxpy_wrapper (N_floats, -1.0, (float *) dev_spin_eo1_up, (float *) dev_spin_eo2_up);
  							// dev_spin_eo2_up  =  dev_spin_eo2_up  -  dev_spin_eo1_up  =                      (1+imubar) * dev_spin_eo3_up  -                    epsbar * dev_spin_eo3_dn						//                                                             - (M_oe)*nrm*(1-imubar)*(M_eo) * dev_spin_eo3_up  -  (M_oe)*nrm*epsbar*(M_eo) * dev_spin_eo3_dn
  cublasSaxpy_wrapper (N_floats, -1.0, (float *) dev_spin_eo1_dn, (float *) dev_spin_eo2_dn);
  							// dev_spin_eo2_dn  =  dev_spin_eo2_dn  -  dev_spin_eo1_dn  =                      (1-imubar) * dev_spin_eo3_dn  -                    epsbar * dev_spin_eo3_up
  							//                                                             - (M_oe)*nrm*(1+imubar)*(M_eo) * dev_spin_eo3_dn  -  (M_oe)*nrm*epsbar*(M_eo) * dev_spin_eo3_up
 
  ////////////
  // gamma5 //
  ////////////
  
  // gamma5
  #ifdef RELATIVISTIC_BASIS
    dev_gamma5_rel<<<gridsize3, blocksize3>>>(dev_spin_eo2_up, dev_spin_eo2_up);					// dev_spin_eo2_up = gamma5 * dev_spin_eo2_up 
    dev_gamma5_rel<<<gridsize3, blocksize3>>>(dev_spin_eo2_dn, dev_spin_eo2_dn);					// dev_spin_eo2_dn = gamma5 * dev_spin_eo2_dn  
  #else
    dev_gamma5<<<gridsize3, blocksize3>>>(dev_spin_eo2_up, dev_spin_eo2_up);					// dev_spin_eo2_up = gamma5 * dev_spin_eo2_up 
    dev_gamma5<<<gridsize3, blocksize3>>>(dev_spin_eo2_dn, dev_spin_eo2_dn);					// dev_spin_eo2_dn = gamma5 * dev_spin_eo2_dn
  #endif

  cublasSscal_wrapper (N_floats, phmc_invmaxev*phmc_invmaxev, (float *) dev_spin_eo2_up);
  cublasSscal_wrapper (N_floats, phmc_invmaxev*phmc_invmaxev, (float *) dev_spin_eo2_dn);   
  
  
  
#elif ASYNC == 1  //ASYNC
  
  HOPPING_ASYNC(dev_gf, spinin_dn, dev_spin_eo1_up, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0, gridsize1, blocksize1);
  
  HOPPING_ASYNC(dev_gf, spinin_up, dev_spin_eo1_dn, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0, gridsize1, blocksize1);

  
  // imubar, gamma5
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo1_up, dev_spin_eo2_up, -1.0);		// dev_spin_eo2_up  =  (1 - imubar) * dev_spin_eo1_up  =  (1 - imubar)*(M_eo) * spinin_dn
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo1_dn, dev_spin_eo2_dn, +1.0);		// dev_spin_eo2_dn  =  (1 + imubar) * dev_spin_eo1_dn  =  (1 + imubar)*(M_eo) * spinin_up
  
  #ifdef USETEXTURE
    bind_texture_spin(dev_spin_eo1_up,1);
    bind_texture_spin_dn(dev_spin_eo1_dn,1);   
  #endif  
  dev_nd_linalg2<<<gridsize2, blocksize2>>>(dev_spin_eo1_up, dev_spin_eo1_dn, dev_spin_eo2_up, dev_spin_eo2_dn, g_epsbar, nrm); 
  #ifdef USETEXTURE
    unbind_texture_spin(1);
    unbind_texture_spin_dn(1);    
  #endif
  
  HOPPING_ASYNC(dev_gf, dev_spin_eo2_up, dev_spin_eo1_up, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1, gridsize1, blocksize1);
  
  HOPPING_ASYNC(dev_gf, dev_spin_eo2_dn, dev_spin_eo1_dn, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1, gridsize1, blocksize1);

  
  ////////////
  // (M_oo) //
  ////////////
  
  
  #ifdef USETEXTURE
    bind_texture_spin(spinin_up,1);
    bind_texture_spin_dn(spinin_dn,1);   
  #endif     
  dev_nd_linalg3_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo1_up, dev_spin_eo1_dn, spinin_up, spinin_dn, dev_spin_eo2_up, dev_spin_eo2_dn, g_epsbar, +1.0 ); 
  #ifdef USETEXTURE
    unbind_texture_spin(1);
    unbind_texture_spin_dn(1);    
  #endif 
  
  
  ////////////////////
  // (a,b) -> (b,a) //
  ////////////////////
  dev_copy_spinor_field<<<gridsize4, blocksize4>>>(dev_spin_eo2_dn, dev_spin_eo3_up);			// dev_spin_eo3_up = dev_spin_eo2_dn
  dev_copy_spinor_field<<<gridsize4, blocksize4>>>(dev_spin_eo2_up, dev_spin_eo3_dn);			// dev_spin_eo3_dn = dev_spin_eo2_up
  
  
  
  
  
  ///////////////////////////////////								///////////////////////
  //        Q_tilde(2x2)           //								// (Q_tilde) * (b,a) //
  ///////////////////////////////////								///////////////////////
  
  
  ////////////////////////////
  // (M_oe) (Mee^-1) (M_eo) //
  ////////////////////////////
  
    	  
  HOPPING_ASYNC(dev_gf, dev_spin_eo3_up, dev_spin_eo1_up, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0, gridsize1, blocksize1);
    	  
  HOPPING_ASYNC(dev_gf, dev_spin_eo3_dn, dev_spin_eo1_dn, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0, gridsize1, blocksize1);

  
  
  // imubar, gamma5
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo1_up, dev_spin_eo2_up, -1.0);			// dev_spin_eo2_up  =  (1 - imubar) * dev_spin_eo1_up  =  (1 - imubar)*(M_eo) * dev_spin_eo3_up
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo1_dn, dev_spin_eo2_dn, +1.0);			// dev_spin_eo2_dn  =  (1 + imubar) * dev_spin_eo1_dn  =  (1 + imubar)*(M_eo) * dev_spin_eo3_dn
  
  
  #ifdef USETEXTURE
    bind_texture_spin(dev_spin_eo1_up,1);
    bind_texture_spin_dn(dev_spin_eo1_dn,1);   
  #endif    
  dev_nd_linalg2<<<gridsize2, blocksize2>>>(dev_spin_eo1_up, dev_spin_eo1_dn, dev_spin_eo2_up, dev_spin_eo2_dn, g_epsbar, nrm); 
  #ifdef USETEXTURE
    unbind_texture_spin(1);
    unbind_texture_spin_dn(1);    
  #endif 
    	  
  HOPPING_ASYNC(dev_gf, dev_spin_eo2_up, dev_spin_eo1_up, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1, gridsize1, blocksize1);

  HOPPING_ASYNC(dev_gf, dev_spin_eo2_dn, dev_spin_eo1_dn, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1, gridsize1, blocksize1);

  
  
  
  ////////////
  // (M_oo) //
  ////////////
  
  #ifdef USETEXTURE
    bind_texture_spin(dev_spin_eo3_up,1);
    bind_texture_spin_dn(dev_spin_eo3_dn,1);   
  #endif    
    dev_nd_linalg4_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo1_up, dev_spin_eo1_dn, dev_spin_eo3_up, dev_spin_eo3_dn, dev_spin_eo2_up, dev_spin_eo2_dn, g_epsbar, +1.0 ); 
  #ifdef USETEXTURE
    unbind_texture_spin(1);
    unbind_texture_spin_dn(1);    
  #endif 

  cublasSscal_wrapper (N_floats, phmc_invmaxev*phmc_invmaxev, (float *) dev_spin_eo2_up);
  cublasSscal_wrapper (N_floats, phmc_invmaxev*phmc_invmaxev, (float *) dev_spin_eo2_dn);
  

#endif  //ASYNC  
  /*
  ////////////
  // output //										// output is already done by setting  dev_spin_eo2_up/dn = spinout_up/dn
  ////////////
  
  dev_copy_spinor_field<<<griddim2, blockdim2 >>>(dev_spin_eo2_up, spinout_up);		// spinout_up = dev_spin_eo2_up
  dev_copy_spinor_field<<<griddim2, blockdim2 >>>(dev_spin_eo2_dn, spinout_dn);		// spinout_dn = dev_spin_eo2_dn
  */

  
  return;
  
}

#endif //MPI





void dev_Qtm_pm_ndpsi (dev_spinor * spinout_up, dev_spinor * spinout_dn,
                              dev_spinor * spinin_up , dev_spinor * spinin_dn , 
                              int gridsize1, int blocksize1, int gridsize2, int blocksize2,
                              int gridsize3, int blocksize3, int gridsize4, int blocksize4){
#ifdef _USE_MPI
  dev_Qtm_pm_ndpsi_mpi (spinout_up, spinout_dn,
                    spinin_up , spinin_dn , 
                    gridsize1, blocksize1, gridsize2, blocksize2,
                    gridsize3, blocksize3, gridsize4, blocksize4);
#else
  dev_Qtm_pm_ndpsi_scalar (spinout_up, spinout_dn,
                    spinin_up , spinin_dn , 
                    gridsize1, blocksize1, gridsize2, blocksize2,
                    gridsize3, blocksize3, gridsize4, blocksize4);
#endif

}







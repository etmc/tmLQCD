

__global__ void dev_tm_dirac_kappa_d(dev_su3_2v_d * gf, dev_spinor_d * sin, dev_spinor_d * sout, int * dev_nn){

    int pos,hoppos;
    double4 shelp1[6], ssum[6];
    dev_su3_d gfsmem;
      

  pos= threadIdx.x + blockDim.x*blockIdx.x;  
  
 #ifdef TEMPORALGAUGE
  int spatialvol = dev_LX*dev_LY*dev_LZ;
 #endif
  

  
  if(pos < dev_VOLUME){
        
          //dev_zero_spinor(&(ssum[0])); // zero sum
          //skalarer Term
          dev_read_spinor_d(&(ssum[0]), &sin[pos]);

          
	  
           //positive direction
            hoppos = dev_nn[8*pos];
             //hoppos = tex1Dfetch(nn_tex,8*pos);
            //color
            
            #ifdef TEMPORALGAUGE
              // gf == ID for t != T-1 => just read the spinor
              #ifdef _USE_MPI
                if ( ((pos) < (dev_T-1)*spatialvol) || (dev_rank < dev_nproc-1) ) {
                //if ((pos) < (dev_T-1)*spatialvol) { // FAKE TEMPORALGAUGE
              #else
                if ((pos/spatialvol) != (dev_T-1) ) {
              #endif
              
		#ifdef RELATIVISTIC_BASIS
		  DIRECTLOAD_REL_UP_D(shelp1, sin, hoppos);
		#else
		  DIRECTLOAD_D(shelp1, sin, hoppos);		
		#endif
              }
              else{
                // gf != ID for t == T-1 => mult spinor with gf
                dev_reconstructgf_2vtexref_d(gf, 4*pos+0, &(gfsmem));
                
                
                #ifdef RELATIVISTIC_BASIS
                    dev_su3MtV_rel_up_d(gfsmem, &(sin[hoppos]), &(shelp1[0]));
                #else
                    dev_su3MtV_d(gfsmem, &(sin[hoppos]), &(shelp1[0]));
                #endif
              }
            #else
              dev_reconstructgf_2vtexref_d(gf, 4*pos+0, &(gfsmem));
              #ifdef RELATIVISTIC_BASIS
                  dev_su3MtV_rel_up_d(gfsmem, &(sin[hoppos]), &(shelp1[0]));
              #else
                  dev_su3MtV_d(gfsmem, &(sin[hoppos]), &(shelp1[0]));
              #endif
            #endif
            
            #ifdef RELATIVISTIC_BASIS
              dev_kappaP0_plus_relativistic_d(&(ssum[0]), &(shelp1[0]), dev_cconj_d(dev_k0_d));
            #else
            //-kappa(r - gamma_mu)
              dev_kappaP0_plus_d(&(ssum[0]), &(shelp1[0]), dev_cconj_d(dev_k0_d));
            #endif
	    
//l==0,t
            //negative direction
            hoppos = dev_nn[8*pos+4]; 
             //hoppos = tex1Dfetch(nn_tex,8*pos+4);
            //color
            #ifdef TEMPORALGAUGE
              // gf == ID for t != T-1 => just read the spinor
              #ifdef _USE_MPI
                if ( ((hoppos) < (dev_T-1)*spatialvol) || (dev_rank > 0) ) {
                //if ((hoppos) < (dev_T-1)*spatialvol) { // FAKE TEMPORALGAUGE
              #else
                if ((hoppos/spatialvol) != (dev_T-1) ) {
              #endif
              
		#ifdef RELATIVISTIC_BASIS
		  DIRECTLOAD_REL_DN_D(shelp1, sin, hoppos);   
		#else 		
		  DIRECTLOAD_D(shelp1, sin, hoppos);
		#endif
              }
              else{
                // gf != ID for t == T-1 => mult spinor with gf
                dev_reconstructgf_2vtexref_dagger_d(gf, 4*hoppos+0,&(gfsmem));
                #ifdef RELATIVISTIC_BASIS
                    dev_su3MtV_rel_dn_d(gfsmem, &(sin[hoppos]), &(shelp1[0]));
                #else
                    dev_su3MtV_d(gfsmem, &(sin[hoppos]), &(shelp1[0]));
                #endif
              }
            #else            
              dev_reconstructgf_2vtexref_dagger_d(gf, 4*hoppos+0, &(gfsmem));
              #ifdef RELATIVISTIC_BASIS
                  dev_su3MtV_rel_dn_d(gfsmem, &(sin[hoppos]), &(shelp1[0]));
              #else
                  dev_su3MtV_d(gfsmem, &(sin[hoppos]), &(shelp1[0]));
              #endif
            #endif
            
            #ifdef RELATIVISTIC_BASIS
              dev_kappaP0_minus_relativistic_d(&(ssum[0]), &(shelp1[0]), dev_k0_d);
            #else
              //-kappa(r + gamma_mu)
              dev_kappaP0_minus_d(&(ssum[0]), &(shelp1[0]), dev_k0_d);
            #endif




//l==3,z 
            //positive direction
            hoppos = dev_nn[8*pos+3];
             //hoppos = tex1Dfetch(nn_tex,8*pos+3);
            //color
            dev_reconstructgf_2vtexref_d(gf, 4*pos+3, &(gfsmem));
            dev_su3MtV_kappaP3_plus_d(gfsmem,&(sin[hoppos]), &(ssum[0]), dev_k3_d.re);
            
            

//l==3,z               
            
            //negative direction
            hoppos = dev_nn[8*pos+7];
             //hoppos = tex1Dfetch(nn_tex,8*pos+7); 
            //color
            dev_reconstructgf_2vtexref_dagger_d(gf, 4*hoppos+3, &(gfsmem));
	    dev_su3MtV_kappaP3_minus_d(gfsmem,&(sin[hoppos]), &(ssum[0]), dev_k3_d.re);
            
            



//l==2,y 
            //positive direction
            hoppos = dev_nn[8*pos+2];
             //hoppos = tex1Dfetch(nn_tex,8*pos+2);
            //color
            dev_reconstructgf_2vtexref_d(gf, 4*pos+2, &(gfsmem));
	    dev_su3MtV_kappaP2_plus_d(gfsmem,&(sin[hoppos]), &(ssum[0]), dev_k2_d.re);
     
            


//l==2,y        

            
            //negative direction
            hoppos = dev_nn[8*pos+6]; 
             //hoppos = tex1Dfetch(nn_tex,8*pos+6);
            //color
            dev_reconstructgf_2vtexref_dagger_d(gf, 4*hoppos+2, &(gfsmem));
	    dev_su3MtV_kappaP2_minus_d(gfsmem,&(sin[hoppos]), &(ssum[0]), dev_k2_d.re);
            
            
            



//l==1,x 
            //positive direction
            hoppos = dev_nn[8*pos+1];
             //hoppos = tex1Dfetch(nn_tex,8*pos+1);
            //color
            dev_reconstructgf_2vtexref_d(gf, 4*pos+1, &(gfsmem));
	    dev_su3MtV_kappaP1_plus_d(gfsmem,&(sin[hoppos]), &(ssum[0]), dev_k1_d.re);
            
            
            


//l==1,x 
            
            //negative direction
            hoppos = dev_nn[8*pos+5]; 
             //hoppos = tex1Dfetch(nn_tex,8*pos+5);
            //color
            dev_reconstructgf_2vtexref_dagger_d(gf, 4*hoppos+1, &(gfsmem));
	    dev_su3MtV_kappaP1_minus_d(gfsmem,&(sin[hoppos]), &(ssum[0]), dev_k1_d.re);
            
     
          
          //gamma5 term
          dev_read_spinor_d(&(shelp1[0]), &sin[pos]);
          
          
          //dev_GammatV(4,&(shelp1[0]));
	  #ifdef RELATIVISTIC_BASIS
            dev_Gamma5_rel_d(&(shelp1[0]));          
          #else
            dev_Gamma5_d(&(shelp1[0]));
	  #endif
          dev_complexmult_add_assign_writetoglobal_spinor_d(&(ssum[0]),dev_initcomplex_d(0.0,2.0*kappa_d*mu_d),&(shelp1[0]), &(sout[pos]));
  }
}














__global__ void dev_swapmu_d(){
  if(blockIdx.x == 0 && threadIdx.x == 0){
    mu_d = - mu_d;
  }
}






void dev_Q_pm_psi_d(dev_spinor_d* spinin, dev_spinor_d* spinout){

  //mu -> -mu
  dev_swapmu_d <<<1,1>>> ();
 
  //D_tm 
  dev_tm_dirac_kappa_d <<<gpu_gd_M_d, gpu_bd_M_d >>> (dev_gf_d, spinin, dev_spin_eo1_d, dev_nn);
  
  //GAMMA5
  #ifdef RELATIVISTIC_BASIS 
     dev_gamma5_rel_d <<<gpu_gd_linalg_d, gpu_bd_linalg_d >>> (dev_spin_eo1_d, dev_spin_eo2_d);
  #else
     dev_gamma5_d <<<gpu_gd_linalg_d, gpu_bd_linalg_d >>> (dev_spin_eo1_d, dev_spin_eo2_d);     
  #endif

  //mu -> -mu
  dev_swapmu_d <<<1,1>>> ();

  //D_tm
  dev_tm_dirac_kappa_d <<<gpu_gd_M_d, gpu_bd_M_d >>> (dev_gf_d, dev_spin_eo2_d, dev_spin_eo1_d, dev_nn);
  
  //GAMMA5
  #ifdef RELATIVISTIC_BASIS    
     dev_gamma5_rel_d <<<gpu_gd_linalg_d, gpu_bd_linalg_d >>>(dev_spin_eo1_d, spinout);
  #else
     dev_gamma5_d <<<gpu_gd_linalg_d, gpu_bd_linalg_d >>> (dev_spin_eo1_d, spinout);     
  #endif
     
}














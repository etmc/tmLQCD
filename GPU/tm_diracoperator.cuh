/***********************************************************************
 *
 * Copyright (C) 2010 Florian Burger
 *
 * This file is part of tmLQCD.
 *
 * tmLQCD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * tmLQCD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with tmLQCD.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  
 * File: tm_diracoperator.cuh
 *
 * CUDA twisted mass dirac operator and its adjoint
 *
 * 
 *
 **************************************************************************/






//applies the full tm Operator
// uses texture cache (spin_tex) for input spinor
// runs through whole lattice for output spinor
// D_psi uses phase_mu and not ka_mu for the boundary conds (vice versa in HoppingMatrix) 
// -> thats why complexmult and complexcgmult are interchanged in dev_HoppingMatrix and in 
// dev_tm_dirac_kappa
__global__ void dev_tm_dirac_kappa(dev_su3_2v * gf, dev_spinor * sin, dev_spinor * sout, int * dev_nn){

    int pos,hoppos;
    dev_spinor shelp1[6], ssum[6];
    dev_su3_pad gfsmem;
      

  pos= threadIdx.x + blockDim.x*blockIdx.x;  
  int gaugevol = dev_VOLUME;
  
 #ifdef TEMPORALGAUGE
  int spatialvol = dev_LX*dev_LY*dev_LZ;
 #endif
  

  
  if(pos < dev_VOLUME){
        
          //dev_zero_spinor(&(ssum[0])); // zero sum
          //skalarer Term
         #ifdef USETEXTURE
          ssum[0] = tex1Dfetch(spin_tex0,pos);
          ssum[1] = tex1Dfetch(spin_tex1,pos);
          ssum[2] = tex1Dfetch(spin_tex2,pos);
          ssum[3] = tex1Dfetch(spin_tex3,pos);
          ssum[4] = tex1Dfetch(spin_tex4,pos);
          ssum[5] = tex1Dfetch(spin_tex5,pos);
	 #else
	  ssum[0] = sin[pos+0*DEVOFF];
          ssum[1] = sin[pos+1*DEVOFF];
          ssum[2] = sin[pos+2*DEVOFF];
          ssum[3] = sin[pos+3*DEVOFF];
          ssum[4] = sin[pos+4*DEVOFF];
          ssum[5] = sin[pos+5*DEVOFF];
	 #endif
          

	  
           //positive direction
            hoppos = dev_nn[8*pos];
             //hoppos = tex1Dfetch(nn_tex,8*pos);
            //color
            
            #ifdef TEMPORALGAUGE
              // gf == ID for t != T-1 => just read the spinor
              #ifdef MPI
                if ( ((pos) < (dev_T-1)*spatialvol) || (dev_rank < dev_nproc-1) ) {
                //if ((pos) < (dev_T-1)*spatialvol) { // FAKE TEMPORALGAUGE
              #else
                if ((pos/spatialvol) != (dev_T-1) ) {
              #endif
              
              #ifdef USETEXTURE 
                shelp1[0] = tex1Dfetch(spin_tex0,hoppos);
                shelp1[1] = tex1Dfetch(spin_tex1,hoppos);
                shelp1[2] = tex1Dfetch(spin_tex2,hoppos);
                #ifdef RELATIVISTIC_BASIS
                  shelp1[3].x = 0.0f; shelp1[3].y = 0.0f; shelp1[3].z = 0.0f; shelp1[3].w = 0.0f;
		  shelp1[4].x = 0.0f; shelp1[4].y = 0.0f; shelp1[4].z = 0.0f; shelp1[4].w = 0.0f;
		  shelp1[5].x = 0.0f; shelp1[5].y = 0.0f; shelp1[5].z = 0.0f; shelp1[5].w = 0.0f;
		#else
		  shelp1[3] = tex1Dfetch(spin_tex3,hoppos);
                  shelp1[4] = tex1Dfetch(spin_tex4,hoppos);
                  shelp1[5] = tex1Dfetch(spin_tex5,hoppos);
		#endif
              #else
                shelp1[0] = sin[hoppos+0*DEVOFF];
                shelp1[1] = sin[hoppos+1*DEVOFF];
                shelp1[2] = sin[hoppos+2*DEVOFF];
                #ifdef RELATIVISTIC_BASIS
                  shelp1[3].x = 0.0f; shelp1[3].y = 0.0f; shelp1[3].z = 0.0f; shelp1[3].w = 0.0f;
		  shelp1[4].x = 0.0f; shelp1[4].y = 0.0f; shelp1[4].z = 0.0f; shelp1[4].w = 0.0f;
		  shelp1[5].x = 0.0f; shelp1[5].y = 0.0f; shelp1[5].z = 0.0f; shelp1[5].w = 0.0f;
                #else
		  shelp1[3] = sin[hoppos+3*DEVOFF];
                  shelp1[4] = sin[hoppos+4*DEVOFF];
                  shelp1[5] = sin[hoppos+5*DEVOFF];
                #endif
	      #endif
              }
              else{
                // gf != ID for t == T-1 => mult spinor with gf
                #ifdef GF_8
                dev_reconstructgf_8texref(gf, pos, 0, gaugevol ,&(gfsmem.m));
                #else
                dev_reconstructgf_2vtexref(gf, pos, 0, gaugevol ,&(gfsmem.m));
                #endif
                
                #ifdef RELATIVISTIC_BASIS
                  #ifdef USETEXTURE
                    dev_su3MtV_spintex_rel_up(gfsmem.m, hoppos, &(shelp1[0]));
                  #else
                    dev_su3MtV_rel_up(gfsmem.m, &(sin[hoppos]), &(shelp1[0]));
                  #endif
                #else
                  #ifdef USETEXTURE
                    dev_su3MtV_spintex(gfsmem.m, hoppos, &(shelp1[0]));
                  #else
                    dev_su3MtV(gfsmem.m, &(sin[hoppos]), &(shelp1[0]));
                  #endif
                #endif
              }
            #else
              #ifdef GF_8
              dev_reconstructgf_8texref(gf, pos, 0, gaugevol ,&(gfsmem.m));
              #else
              dev_reconstructgf_2vtexref(gf, pos, 0, gaugevol ,&(gfsmem.m));
              #endif
              #ifdef RELATIVISTIC_BASIS
                #ifdef USETEXTURE
                  dev_su3MtV_spintex_rel_up(gfsmem.m, hoppos, &(shelp1[0]));
                #else
                  dev_su3MtV_rel_up(gfsmem.m, &(sin[hoppos]), &(shelp1[0]));
                #endif
              #else
                #ifdef USETEXTURE
                  dev_su3MtV_spintex(gfsmem.m, hoppos, &(shelp1[0]));
                #else
                  dev_su3MtV(gfsmem.m, &(sin[hoppos]), &(shelp1[0]));
                #endif
              #endif
            #endif
            
            #ifdef RELATIVISTIC_BASIS
              dev_kappaP0_plus_relativistic(&(ssum[0]), &(shelp1[0]), dev_cconj(dev_k0));
            #else
            //-kappa(r - gamma_mu)
              dev_kappaP0_plus(&(ssum[0]), &(shelp1[0]), dev_cconj(dev_k0));
            #endif
	    
//l==0,t
            //negative direction
            hoppos = dev_nn[8*pos+4]; 
             //hoppos = tex1Dfetch(nn_tex,8*pos+4);
            //color
            #ifdef TEMPORALGAUGE
              // gf == ID for t != T-1 => just read the spinor
              #ifdef MPI
                if ( ((hoppos) < (dev_T-1)*spatialvol) || (dev_rank > 0) ) {
                //if ((hoppos) < (dev_T-1)*spatialvol) { // FAKE TEMPORALGAUGE
              #else
                if ((hoppos/spatialvol) != (dev_T-1) ) {
              #endif
              
               #ifdef USETEXTURE
                #ifdef RELATIVISTIC_BASIS
                  shelp1[0].x = 0.0f; shelp1[0].y = 0.0f; shelp1[0].z = 0.0f; shelp1[0].w = 0.0f;
		  shelp1[1].x = 0.0f; shelp1[1].y = 0.0f; shelp1[1].z = 0.0f; shelp1[1].w = 0.0f;
		  shelp1[2].x = 0.0f; shelp1[2].y = 0.0f; shelp1[2].z = 0.0f; shelp1[2].w = 0.0f;
                #else
                  shelp1[0] = tex1Dfetch(spin_tex0,hoppos);
                  shelp1[1] = tex1Dfetch(spin_tex1,hoppos);
                  shelp1[2] = tex1Dfetch(spin_tex2,hoppos);
                #endif
		shelp1[3] = tex1Dfetch(spin_tex3,hoppos);
                shelp1[4] = tex1Dfetch(spin_tex4,hoppos);
                shelp1[5] = tex1Dfetch(spin_tex5,hoppos);
               #else
                #ifdef RELATIVISTIC_BASIS
                  shelp1[0].x = 0.0f; shelp1[0].y = 0.0f; shelp1[0].z = 0.0f; shelp1[0].w = 0.0f;
		  shelp1[1].x = 0.0f; shelp1[1].y = 0.0f; shelp1[1].z = 0.0f; shelp1[1].w = 0.0f;
		  shelp1[2].x = 0.0f; shelp1[2].y = 0.0f; shelp1[2].z = 0.0f; shelp1[2].w = 0.0f;
                #else
                  shelp1[0] = sin[hoppos+0*DEVOFF];
                  shelp1[1] = sin[hoppos+1*DEVOFF];
                  shelp1[2] = sin[hoppos+2*DEVOFF];
                #endif
		shelp1[3] = sin[hoppos+3*DEVOFF];
                shelp1[4] = sin[hoppos+4*DEVOFF];
                shelp1[5] = sin[hoppos+5*DEVOFF];
               #endif
              }
              else{
                // gf != ID for t == T-1 => mult spinor with gf
                #ifdef GF_8
                dev_reconstructgf_8texref_dagger(gf,hoppos, 0, gaugevol ,&(gfsmem.m));
                #else
                dev_reconstructgf_2vtexref_dagger(gf,hoppos, 0, gaugevol ,&(gfsmem.m));
                #endif
                #ifdef RELATIVISTIC_BASIS
                  #ifdef USETEXTURE
                    dev_su3MtV_spintex_rel_down(gfsmem.m, hoppos, &(shelp1[0]));
                  #else
                    dev_su3MtV_rel_down(gfsmem.m, &(sin[hoppos]), &(shelp1[0]));
                  #endif
                #else
                  #ifdef USETEXTURE
                    dev_su3MtV_spintex(gfsmem.m, hoppos, &(shelp1[0]));
                  #else
                    dev_su3MtV(gfsmem.m, &(sin[hoppos]), &(shelp1[0]));
                  #endif
                #endif
              }
            #else            
              #ifdef GF_8
              dev_reconstructgf_8texref_dagger(gf,hoppos, 0, gaugevol ,&(gfsmem.m));
              #else
              dev_reconstructgf_2vtexref_dagger(gf,hoppos, 0, gaugevol ,&(gfsmem.m));
              #endif
              #ifdef RELATIVISTIC_BASIS
                #ifdef USETEXTURE
                  dev_su3MtV_spintex_rel_down(gfsmem.m, hoppos, &(shelp1[0]));  
                #else
                  dev_su3MtV_rel_down(gfsmem.m, &(sin[hoppos]), &(shelp1[0]));
                #endif 
              #else
                #ifdef USETEXTURE
                  dev_su3MtV_spintex(gfsmem.m, hoppos, &(shelp1[0]));  
                #else
                  dev_su3MtV(gfsmem.m, &(sin[hoppos]), &(shelp1[0]));
                #endif 
              #endif
            #endif
            
            #ifdef RELATIVISTIC_BASIS
              dev_kappaP0_minus_relativistic(&(ssum[0]), &(shelp1[0]), dev_k0);
            #else
              //-kappa(r + gamma_mu)
              dev_kappaP0_minus(&(ssum[0]), &(shelp1[0]), dev_k0);
            #endif




//l==3,z 
            //positive direction
            hoppos = dev_nn[8*pos+3];
             //hoppos = tex1Dfetch(nn_tex,8*pos+3);
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref(gf,pos, 3, gaugevol ,&(gfsmem.m));
            #else
            dev_reconstructgf_2vtexref(gf, pos, 3, gaugevol ,&(gfsmem.m));
            #endif
            #ifdef USETEXTURE
              //dev_su3MtV_spintex(gfsmem.m, hoppos, &(shelp1[0]));
              //-kappa(r - gamma_mu), no spin projection optimization yet  
              //dev_kappaP3_plus(&(ssum[0]), &(shelp1[0]), dev_k3.re);
	      
              dev_su3MtV_kappaP3_plus_spintex(gfsmem.m,hoppos, &(ssum[0]), dev_k3);	      
	    #else
              dev_su3MtV_kappaP3_plus(gfsmem.m,&(sin[hoppos]), &(ssum[0]), dev_k3);
            #endif
            

//l==3,z               
            
            //negative direction
            hoppos = dev_nn[8*pos+7];
             //hoppos = tex1Dfetch(nn_tex,8*pos+7); 
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref_dagger(gf,hoppos, 3, gaugevol ,&(gfsmem.m));
            #else
            dev_reconstructgf_2vtexref_dagger(gf,hoppos, 3, gaugevol ,&(gfsmem.m));
            #endif
            #ifdef USETEXTURE
              //dev_su3MtV_spintex(gfsmem.m, hoppos, &(shelp1[0]));
              //-kappa(r + gamma_mu), no spin projection optimization yet 
              //dev_kappaP3_minus(&(ssum[0]), &(shelp1[0]), dev_k3.re);
	      
	      dev_su3MtV_kappaP3_minus_spintex(gfsmem.m,hoppos, &(ssum[0]), dev_k3);	      
	    #else
	      dev_su3MtV_kappaP3_minus(gfsmem.m,&(sin[hoppos]), &(ssum[0]), dev_k3);
            #endif
            





//l==2,y 
            //positive direction
            hoppos = dev_nn[8*pos+2];
             //hoppos = tex1Dfetch(nn_tex,8*pos+2);
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref(gf,pos, 2, gaugevol ,&(gfsmem.m));
            #else
            dev_reconstructgf_2vtexref(gf,pos, 2, gaugevol ,&(gfsmem.m));
            #endif
            #ifdef USETEXTURE
              //dev_su3MtV_spintex(gfsmem.m, hoppos, &(shelp1[0]));
              //-kappa(r - gamma_mu), no spin projection optimization yet
              //dev_kappaP2_plus(&(ssum[0]), &(shelp1[0]), dev_k2.re);
	      
	      dev_su3MtV_kappaP2_plus_spintex(gfsmem.m,hoppos, &(ssum[0]), dev_k2);	      
	    #else
	      dev_su3MtV_kappaP2_plus(gfsmem.m,&(sin[hoppos]), &(ssum[0]), dev_k2);
            #endif
            


//l==2,y        

            
            //negative direction
            hoppos = dev_nn[8*pos+6]; 
             //hoppos = tex1Dfetch(nn_tex,8*pos+6);
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref_dagger(gf,hoppos, 2, gaugevol ,&(gfsmem.m));
            #else
            dev_reconstructgf_2vtexref_dagger(gf,hoppos, 2, gaugevol ,&(gfsmem.m));
            #endif
            #ifdef USETEXTURE
              //dev_su3MtV_spintex(gfsmem.m, hoppos, &(shelp1[0]));
              //-kappa(r + gamma_mu), no spin projection optimization yet
              //dev_kappaP2_minus(&(ssum[0]), &(shelp1[0]), dev_k2.re);
	      
	      dev_su3MtV_kappaP2_minus_spintex(gfsmem.m,hoppos, &(ssum[0]), dev_k2);	      
	    #else
	      dev_su3MtV_kappaP2_minus(gfsmem.m,&(sin[hoppos]), &(ssum[0]), dev_k2);
            #endif
            
            



//l==1,x 
            //positive direction
            hoppos = dev_nn[8*pos+1];
             //hoppos = tex1Dfetch(nn_tex,8*pos+1);
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref(gf,pos, 1, gaugevol ,&(gfsmem.m));
            #else
            dev_reconstructgf_2vtexref(gf,pos, 1, gaugevol ,&(gfsmem.m));
            #endif
            #ifdef USETEXTURE
              //dev_su3MtV_spintex(gfsmem.m, hoppos, &(shelp1[0]));
              //-kappa(r - gamma_mu), no spin projection optimization yet
              //dev_kappaP1_plus(&(ssum[0]), &(shelp1[0]), dev_k1.re);
	      
	      dev_su3MtV_kappaP1_plus_spintex(gfsmem.m,hoppos, &(ssum[0]), dev_k1);	      
	    #else
	      dev_su3MtV_kappaP1_plus(gfsmem.m,&(sin[hoppos]), &(ssum[0]), dev_k1);
            #endif
            
            


//l==1,x 
            
            //negative direction
            hoppos = dev_nn[8*pos+5]; 
             //hoppos = tex1Dfetch(nn_tex,8*pos+5);
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref_dagger(gf,hoppos, 1, gaugevol ,&(gfsmem.m));
            #else
            dev_reconstructgf_2vtexref_dagger(gf,hoppos, 1, gaugevol ,&(gfsmem.m));
            #endif
            #ifdef USETEXTURE
              //dev_su3MtV_spintex(gfsmem.m, hoppos, &(shelp1[0]));
              //-kappa(r + gamma_mu), no spin projection optimization yet
              //dev_kappaP1_minus(&(ssum[0]), &(shelp1[0]), dev_k1.re);
	      
	      dev_su3MtV_kappaP1_minus_spintex(gfsmem.m,hoppos, &(ssum[0]), dev_k1);	      
	    #else
	      dev_su3MtV_kappaP1_minus(gfsmem.m,&(sin[hoppos]), &(ssum[0]), dev_k1);
            #endif
            	  
	  
	  
	  
          
          
          //gamma5 term
         #ifdef USETEXTURE
          shelp1[0] = tex1Dfetch(spin_tex0,pos);
          shelp1[1] = tex1Dfetch(spin_tex1,pos);
          shelp1[2] = tex1Dfetch(spin_tex2,pos);
          shelp1[3] = tex1Dfetch(spin_tex3,pos);
          shelp1[4] = tex1Dfetch(spin_tex4,pos);
          shelp1[5] = tex1Dfetch(spin_tex5,pos);
         #else
          shelp1[0] = sin[pos+0*DEVOFF];
          shelp1[1] = sin[pos+1*DEVOFF];
          shelp1[2] = sin[pos+2*DEVOFF];
          shelp1[3] = sin[pos+3*DEVOFF];
          shelp1[4] = sin[pos+4*DEVOFF];
          shelp1[5] = sin[pos+5*DEVOFF];
         #endif 
          
          
          //dev_GammatV(4,&(shelp1[0]));
	  #ifdef RELATIVISTIC_BASIS
            dev_Gamma5_rel(&(shelp1[0]));          
          #else
            dev_Gamma5(&(shelp1[0]));
	  #endif
          dev_complexmult_add_assign_writetoglobal_spinor(&(ssum[0]),dev_initcomplex(0.0,2.0*kappa*mu),&(shelp1[0]), &(sout[pos]));
  }
}




__global__ void dev_gamma5(dev_spinor * sin, dev_spinor * sout){
  int pos;
  pos= threadIdx.x + blockDim.x*blockIdx.x;  
  if(pos < dev_VOLUME){
          sout[pos+0*DEVOFF].x = sin[pos+0*DEVOFF].x;
          sout[pos+0*DEVOFF].y = sin[pos+0*DEVOFF].y;
          sout[pos+0*DEVOFF].z = sin[pos+0*DEVOFF].z;
          sout[pos+0*DEVOFF].w = sin[pos+0*DEVOFF].w;
          sout[pos+1*DEVOFF].x = sin[pos+1*DEVOFF].x;
          sout[pos+1*DEVOFF].y = sin[pos+1*DEVOFF].y;
          
          sout[pos+1*DEVOFF].z = sin[pos+1*DEVOFF].z;
          sout[pos+1*DEVOFF].w = sin[pos+1*DEVOFF].w;
          sout[pos+2*DEVOFF].x = sin[pos+2*DEVOFF].x;
          sout[pos+2*DEVOFF].y = sin[pos+2*DEVOFF].y;
          sout[pos+2*DEVOFF].z = sin[pos+2*DEVOFF].z;
          sout[pos+2*DEVOFF].w = sin[pos+2*DEVOFF].w;   
          
          sout[pos+3*DEVOFF].x = -1.0*sin[pos+3*DEVOFF].x;
          sout[pos+3*DEVOFF].y = -1.0*sin[pos+3*DEVOFF].y;
          sout[pos+3*DEVOFF].z = -1.0*sin[pos+3*DEVOFF].z;
          sout[pos+3*DEVOFF].w = -1.0*sin[pos+3*DEVOFF].w;
          sout[pos+4*DEVOFF].x = -1.0*sin[pos+4*DEVOFF].x;
          sout[pos+4*DEVOFF].y = -1.0*sin[pos+4*DEVOFF].y; 

          sout[pos+4*DEVOFF].z = -1.0*sin[pos+4*DEVOFF].z;
          sout[pos+4*DEVOFF].w = -1.0*sin[pos+4*DEVOFF].w;
          sout[pos+5*DEVOFF].x = -1.0*sin[pos+5*DEVOFF].x;
          sout[pos+5*DEVOFF].y = -1.0*sin[pos+5*DEVOFF].y;
          sout[pos+5*DEVOFF].z = -1.0*sin[pos+5*DEVOFF].z;
          sout[pos+5*DEVOFF].w = -1.0*sin[pos+5*DEVOFF].w;                 
  } 
}



  

__global__ void dev_gamma5_rel(dev_spinor * sin, dev_spinor * sout){
  int pos;
  pos= threadIdx.x + blockDim.x*blockIdx.x;  
  if(pos < dev_VOLUME){
  dev_spinor shelp[6];
  
  shelp[0].x = sin[pos+0*DEVOFF].x;
  shelp[0].y = sin[pos+0*DEVOFF].y;
  shelp[0].z = sin[pos+0*DEVOFF].z;
  shelp[0].w = sin[pos+0*DEVOFF].w;
  
  shelp[1].x = sin[pos+1*DEVOFF].x;
  shelp[1].y = sin[pos+1*DEVOFF].y;
  shelp[1].z = sin[pos+1*DEVOFF].z;
  shelp[1].w = sin[pos+1*DEVOFF].w;
  
  shelp[2].x = sin[pos+2*DEVOFF].x;
  shelp[2].y = sin[pos+2*DEVOFF].y;
  shelp[2].z = sin[pos+2*DEVOFF].z;
  shelp[2].w = sin[pos+2*DEVOFF].w; 
  
  shelp[3].x = sin[pos+3*DEVOFF].x;
  shelp[3].y = sin[pos+3*DEVOFF].y;
  shelp[3].z = sin[pos+3*DEVOFF].z;
  shelp[3].w = sin[pos+3*DEVOFF].w;   
  
  shelp[4].x = sin[pos+4*DEVOFF].x;
  shelp[4].y = sin[pos+4*DEVOFF].y;
  shelp[4].z = sin[pos+4*DEVOFF].z;
  shelp[4].w = sin[pos+4*DEVOFF].w;  
  
  shelp[5].x = sin[pos+5*DEVOFF].x;
  shelp[5].y = sin[pos+5*DEVOFF].y;
  shelp[5].z = sin[pos+5*DEVOFF].z;
  shelp[5].w = sin[pos+5*DEVOFF].w;  
  
  
 sout[pos+3*DEVOFF].x  = -shelp[0].x; 
 sout[pos+3*DEVOFF].y  = -shelp[0].y;
 sout[pos+3*DEVOFF].z  = -shelp[0].z;
 sout[pos+3*DEVOFF].w  = -shelp[0].w;
 sout[pos+4*DEVOFF].x  = -shelp[1].x;
 sout[pos+4*DEVOFF].y  = -shelp[1].y;
 

 sout[pos+4*DEVOFF].z  = -shelp[1].z;
 sout[pos+4*DEVOFF].w  = -shelp[1].w;
 sout[pos+5*DEVOFF].x  = -shelp[2].x;
 sout[pos+5*DEVOFF].y  = -shelp[2].y;
 sout[pos+5*DEVOFF].z  = -shelp[2].z;
 sout[pos+5*DEVOFF].w  = -shelp[2].w;

 sout[pos+0*DEVOFF].x = -shelp[3].x;
 sout[pos+0*DEVOFF].y = -shelp[3].y; 
 sout[pos+0*DEVOFF].z = -shelp[3].z;
 sout[pos+0*DEVOFF].w = -shelp[3].w;
 sout[pos+1*DEVOFF].x = -shelp[4].x;
 sout[pos+1*DEVOFF].y = -shelp[4].y;
 
 
 sout[pos+1*DEVOFF].z = -shelp[4].z;
 sout[pos+1*DEVOFF].w = -shelp[4].w; 
 sout[pos+2*DEVOFF].x = -shelp[5].x;
 sout[pos+2*DEVOFF].y = -shelp[5].y; 
 sout[pos+2*DEVOFF].z = -shelp[5].z;
 sout[pos+2*DEVOFF].w = -shelp[5].w;

               
  } 
}




extern "C" void dev_tm_dirac_dagger_kappa(dev_su3_2v * gf,dev_spinor* spinin, dev_spinor* spinout, 
 int *grid, int * nn_grid, float* output, float* erg, int xsize, int ysize){
 int gridsize;
 if( VOLUME >= 128){
   gridsize =VOLUME/128;
 }
 else{
   gridsize=1;
 }
 dim3 griddim2(gridsize,1,1);
 dim3 blockdim2(128,1,1);
 dim3 blockdim(xsize,ysize);
 
  dim3 blockdim3(BLOCK,1,1);
 if( VOLUME >= BLOCK){
   gridsize = (int)(VOLUME/BLOCK) + 1;
 }
 else{
   gridsize=1;
 }
 dim3 griddim3(gridsize,1,1); 
  dev_gamma5 <<<griddim2, blockdim2 >>> (spinin,spinout);
  dev_tm_dirac_kappa <<<griddim3, blockdim3 >>> (gf, spinout, spinin, dev_nn);
  dev_gamma5 <<<griddim2, blockdim2 >>>(spinin,spinout);
}




__global__ void dev_swapmu(){
  if(blockIdx.x == 0 && threadIdx.x == 0){
    mu = - mu;
  }
}






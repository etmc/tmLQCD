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
    __shared__ dev_su3 gfsmem[BLOCK];
    

  pos= threadIdx.x + blockDim.x*blockIdx.x;  
  int ix = threadIdx.x;
  int gaugevol = dev_VOLUME;
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
          
//hopping term                
//l==0,t
            //positive direction
            hoppos = dev_nn[8*pos];
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref(gf,pos, 0, gaugevol ,&(gfsmem[ix]));
            #else
            dev_reconstructgf_2vtexref(gf,pos, 0, gaugevol ,&(gfsmem[ix]));
            #endif
            #ifdef USETEXTURE
              dev_su3MtV_spintex(gfsmem[ix], hoppos, &(shelp1[0]));
            #else
              dev_su3MtV(gfsmem[ix], &(sin[hoppos]), &(shelp1[0]));
            #endif
            //-kappa(r - gamma_mu)
            dev_complexmult_add_assign_spinor(&(ssum[0]),dev_mk0,&(shelp1[0]), &(ssum[0]));
            //dev_GammatV(0,&(shelp1[0]));
            dev_Gamma0(&(shelp1[0]));
            dev_complexmult_add_assign_spinor(&(ssum[0]),dev_k0,&(shelp1[0]), &(ssum[0]));

            //negative direction
            hoppos = dev_nn[8*pos+4];
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref_dagger(gf,hoppos, 0, gaugevol ,&(gfsmem[ix]));
            #else
            dev_reconstructgf_2vtexref_dagger(gf,hoppos, 0, gaugevol ,&(gfsmem[ix]));
            #endif
            #ifdef USETEXTURE
              dev_su3MtV_spintex(gfsmem[ix], hoppos, &(shelp1[0]));  
            #else
              dev_su3MtV(gfsmem[ix], &(sin[hoppos]), &(shelp1[0]));
            #endif    
            //-kappa(r + gamma_mu)
            dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_mk0,&(shelp1[0]), &(ssum[0]));
            //dev_GammatV(0,&(shelp1[0]));
            dev_Gamma0(&(shelp1[0]));
            dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_mk0,&(shelp1[0]), &(ssum[0]));


//l==3,z               
            //positive direction
            hoppos = dev_nn[8*pos+3];
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref(gf,pos, 3, gaugevol ,&(gfsmem[ix]));
            #else
            dev_reconstructgf_2vtexref(gf,pos, 3, gaugevol ,&(gfsmem[ix]));
            #endif
            #ifdef USETEXTURE
              dev_su3MtV_spintex(gfsmem[ix], hoppos, &(shelp1[0]));
            #else
              dev_su3MtV(gfsmem[ix], &(sin[hoppos]), &(shelp1[0]));
            #endif
            //-kappa(r - gamma_mu)
            dev_complexmult_add_assign_spinor(&(ssum[0]),dev_mk3,&(shelp1[0]), &(ssum[0]));
            //dev_GammatV(3,&(shelp1[0]));
            dev_Gamma3(&(shelp1[0]));
            dev_complexmult_add_assign_spinor(&(ssum[0]),dev_k3,&(shelp1[0]), &(ssum[0]));

            //negative direction
            hoppos = dev_nn[8*pos+7];
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref_dagger(gf,hoppos, 3, gaugevol ,&(gfsmem[ix]));
            #else
            dev_reconstructgf_2vtexref_dagger(gf,hoppos, 3, gaugevol ,&(gfsmem[ix]));
            #endif
            #ifdef USETEXTURE
              dev_su3MtV_spintex(gfsmem[ix], hoppos, &(shelp1[0]));
            #else
              dev_su3MtV(gfsmem[ix], &(sin[hoppos]), &(shelp1[0]));
            #endif
            //-kappa(r + gamma_mu)
            dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_mk3,&(shelp1[0]), &(ssum[0]));
            //dev_GammatV(3,&(shelp1[0]));
            dev_Gamma3(&(shelp1[0]));
            dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_mk3,&(shelp1[0]), &(ssum[0]));
         
         
//l==2,y        
            //positive direction
            hoppos = dev_nn[8*pos+2];
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref(gf,pos, 2, gaugevol ,&(gfsmem[ix]));
            #else
            dev_reconstructgf_2vtexref(gf,pos, 2, gaugevol ,&(gfsmem[ix]));
            #endif
            #ifdef USETEXTURE
              dev_su3MtV_spintex(gfsmem[ix], hoppos, &(shelp1[0]));
            #else
              dev_su3MtV(gfsmem[ix], &(sin[hoppos]), &(shelp1[0]));
            #endif
            //-kappa(r - gamma_mu)
            dev_complexmult_add_assign_spinor(&(ssum[0]),dev_mk2,&(shelp1[0]), &(ssum[0]));
            //dev_GammatV(2,&(shelp1[0]));
            dev_Gamma2(&(shelp1[0]));
            dev_complexmult_add_assign_spinor(&(ssum[0]),dev_k2,&(shelp1[0]), &(ssum[0]));
            
            //negative direction
            hoppos = dev_nn[8*pos+6];
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref_dagger(gf,hoppos, 2, gaugevol ,&(gfsmem[ix]));
            #else
            dev_reconstructgf_2vtexref_dagger(gf,hoppos, 2, gaugevol ,&(gfsmem[ix]));
            #endif
            #ifdef USETEXTURE
              dev_su3MtV_spintex(gfsmem[ix], hoppos, &(shelp1[0]));
            #else
              dev_su3MtV(gfsmem[ix], &(sin[hoppos]), &(shelp1[0]));
            #endif
            //-kappa(r + gamma_mu)
            dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_mk2,&(shelp1[0]), &(ssum[0]));
            //dev_GammatV(2,&(shelp1[0]));
            dev_Gamma2(&(shelp1[0]));
            dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_mk2,&(shelp1[0]), &(ssum[0]));


//l==1,x 
            //positive direction
            hoppos = dev_nn[8*pos+1];
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref(gf,pos, 1, gaugevol ,&(gfsmem[ix]));
            #else
            dev_reconstructgf_2vtexref(gf,pos, 1, gaugevol ,&(gfsmem[ix]));
            #endif
            #ifdef USETEXTURE
              dev_su3MtV_spintex(gfsmem[ix], hoppos, &(shelp1[0]));
            #else
              dev_su3MtV(gfsmem[ix], &(sin[hoppos]), &(shelp1[0]));
            #endif
            //-kappa(r - gamma_mu)
            dev_complexmult_add_assign_spinor(&(ssum[0]),dev_mk1,&(shelp1[0]), &(ssum[0]));
            //dev_GammatV(1,&(shelp1[0]));
            dev_Gamma1(&(shelp1[0]));
            dev_complexmult_add_assign_spinor(&(ssum[0]),dev_k1,&(shelp1[0]), &(ssum[0]));

            //negative direction
            hoppos = dev_nn[8*pos+5];
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref_dagger(gf,hoppos, 1, gaugevol ,&(gfsmem[ix]));
            #else
            dev_reconstructgf_2vtexref_dagger(gf,hoppos, 1, gaugevol ,&(gfsmem[ix]));
            #endif
            #ifdef USETEXTURE
              dev_su3MtV_spintex(gfsmem[ix], hoppos, &(shelp1[0]));
            #else
              dev_su3MtV(gfsmem[ix], &(sin[hoppos]), &(shelp1[0]));
            #endif
            //-kappa(r + gamma_mu)
            dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_mk1,&(shelp1[0]), &(ssum[0]));
            //dev_GammatV(1,&(shelp1[0]));
            dev_Gamma1(&(shelp1[0]));
            dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_mk1,&(shelp1[0]), &(ssum[0]));  
          
          
          
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
 
 sout[pos+3*DEVOFF].x  = -sin[pos].x; 
 sout[pos+3*DEVOFF].y  = -sin[pos].y;
 sout[pos+3*DEVOFF].z  = -sin[pos].z;
 sout[pos+3*DEVOFF].w  = -sin[pos].w;
 sout[pos+4*DEVOFF].x  = -sin[pos+1*DEVOFF].x;
 sout[pos+4*DEVOFF].y  = -sin[pos+1*DEVOFF].y;
 

 sout[pos+4*DEVOFF].z  = -sin[pos+1*DEVOFF].z;
 sout[pos+4*DEVOFF].w  = -sin[pos+1*DEVOFF].w;
 sout[pos+5*DEVOFF].x  = -sin[pos+2*DEVOFF].x;
 sout[pos+5*DEVOFF].y  = -sin[pos+2*DEVOFF].y;
 sout[pos+5*DEVOFF].z  = -sin[pos+2*DEVOFF].z;
 sout[pos+5*DEVOFF].w  = -sin[pos+2*DEVOFF].w;


 //this is the lower spinor
 //re and im parts are interchanged w.r.t. above
 //same sign as above 

 sout[pos+0*DEVOFF].x = -sin[pos+3*DEVOFF].x;
 sout[pos+0*DEVOFF].y = -sin[pos+3*DEVOFF].y; 
 sout[pos+0*DEVOFF].z = -sin[pos+3*DEVOFF].z;
 sout[pos+0*DEVOFF].w = -sin[pos+3*DEVOFF].w;
 sout[pos+1*DEVOFF].x = -sin[pos+4*DEVOFF].x;
 sout[pos+1*DEVOFF].y = -sin[pos+4*DEVOFF].y;
 
 
 sout[pos+1*DEVOFF].z = -sin[pos+4*DEVOFF].z;
 sout[pos+1*DEVOFF].w = -sin[pos+4*DEVOFF].w; 
 sout[pos+2*DEVOFF].x = -sin[pos+5*DEVOFF].x;
 sout[pos+2*DEVOFF].y = -sin[pos+5*DEVOFF].y; 
 sout[pos+2*DEVOFF].z = -sin[pos+5*DEVOFF].z;
 sout[pos+2*DEVOFF].w = -sin[pos+5*DEVOFF].w;

               
  } 
}




extern "C" void dev_tm_dirac_dagger_kappa(dev_su3_2v * gf,dev_spinor* spinin, dev_spinor* spinout, 
 int *grid, int * nn_grid, REAL* output,REAL* erg, int xsize, int ysize){
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






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



template<class RealT>
__global__ void dev_tm_dirac_kappa
(
  typename dev_su3_2vT<RealT>::type * gf,
  typename dev_spinorT<RealT>::type * sin,
  typename dev_spinorT<RealT>::type * sout,
  int * dev_nn
){
    int pos,hoppos;
    typename dev_spinorT<RealT>::type shelp1[6], ssum[6];
    __shared__ typename dev_su3T<RealT>::type gfsmem[BLOCK];
    

  pos= threadIdx.x + blockDim.x*blockIdx.x;
  int ix = threadIdx.x;
  if(pos < dev_VOLUME){
        
          //dev_zero_spinor(&(ssum[0])); // zero sum
          //skalarer Term
         #ifdef USETEXTURE
          ssum[0] = tex1Dfetch(spin_tex,6*pos);
          ssum[1] = tex1Dfetch(spin_tex,6*pos+1);
          ssum[2] = tex1Dfetch(spin_tex,6*pos+2);
          ssum[3] = tex1Dfetch(spin_tex,6*pos+3);
          ssum[4] = tex1Dfetch(spin_tex,6*pos+4);
          ssum[5] = tex1Dfetch(spin_tex,6*pos+5);
	 #else
	  ssum[0] = sin[6*pos];
          ssum[1] = sin[6*pos+1];
          ssum[2] = sin[6*pos+2];
          ssum[3] = sin[6*pos+3];
          ssum[4] = sin[6*pos+4];
          ssum[5] = sin[6*pos+5];
	 #endif
          
//hopping term                
//l==0,t
            //positive direction
            hoppos = dev_nn[8*pos];
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref <RealT>(gf,4*pos,&(gfsmem[ix]));
            #else
            dev_reconstructgf_2vtexref<RealT>(gf,4*pos,&(gfsmem[ix]));
            #endif
            #ifdef USETEXTURE
              dev_su3MtV_spintex<RealT>(gfsmem[ix], hoppos, &(shelp1[0]));
            #else
              dev_su3MtV        <RealT>(gfsmem[ix], &(sin[6*hoppos]), &(shelp1[0]));
            #endif
            //-kappa(r - gamma_mu)
            dev_complexmult_add_assign_spinor(&(ssum[0]),dev_complexT<RealT>(dev_mk0),&(shelp1[0]), &(ssum[0]));
            //dev_GammatV(0,&(shelp1[0]));
            dev_Gamma0<RealT>(&(shelp1[0]));
            dev_complexmult_add_assign_spinor(&(ssum[0]),dev_complexT<RealT>(dev_k0),&(shelp1[0]), &(ssum[0]));

            //negative direction
            hoppos = dev_nn[8*pos+4];
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref_dagger <RealT>(gf,4*hoppos,&(gfsmem[ix]));
            #else
            dev_reconstructgf_2vtexref_dagger<RealT>(gf,4*hoppos,&(gfsmem[ix]));
            #endif
            #ifdef USETEXTURE
              dev_su3MtV_spintex<RealT>(gfsmem[ix], hoppos, &(shelp1[0]));  
            #else
              dev_su3MtV        <RealT>(gfsmem[ix], &(sin[6*hoppos]), &(shelp1[0]));
            #endif    
            //-kappa(r + gamma_mu)
            dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_complexT<RealT>(dev_mk0),&(shelp1[0]), &(ssum[0]));
            //dev_GammatV(0,&(shelp1[0]));
            dev_Gamma0<RealT>(&(shelp1[0]));
            dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_complexT<RealT>(dev_mk0),&(shelp1[0]), &(ssum[0]));


//l==3,z               
            //positive direction
            hoppos = dev_nn[8*pos+3];
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref <RealT>(gf,4*pos+(3),&(gfsmem[ix]));
            #else
            dev_reconstructgf_2vtexref<RealT>(gf,4*pos+(3),&(gfsmem[ix]));
            #endif
            #ifdef USETEXTURE
              dev_su3MtV_spintex<RealT>(gfsmem[ix], hoppos, &(shelp1[0]));
            #else
              dev_su3MtV        <RealT>(gfsmem[ix], &(sin[6*hoppos]), &(shelp1[0]));
            #endif
            //-kappa(r - gamma_mu)
            dev_complexmult_add_assign_spinor(&(ssum[0]),dev_complexT<RealT>(dev_mk3),&(shelp1[0]), &(ssum[0]));
            //dev_GammatV(3,&(shelp1[0]));
            dev_Gamma3<RealT>(&(shelp1[0]));
            dev_complexmult_add_assign_spinor(&(ssum[0]),dev_complexT<RealT>(dev_k3),&(shelp1[0]), &(ssum[0]));

            //negative direction
            hoppos = dev_nn[8*pos+7];
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref_dagger <RealT>(gf,4*hoppos+(3),&(gfsmem[ix]));
            #else
            dev_reconstructgf_2vtexref_dagger<RealT>(gf,4*hoppos+(3),&(gfsmem[ix]));
            #endif
            #ifdef USETEXTURE
              dev_su3MtV_spintex<RealT>(gfsmem[ix], hoppos, &(shelp1[0]));
            #else
              dev_su3MtV        <RealT>(gfsmem[ix], &(sin[6*hoppos]), &(shelp1[0]));
            #endif
            //-kappa(r + gamma_mu)
            dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_complexT<RealT>(dev_mk3),&(shelp1[0]), &(ssum[0]));
            //dev_GammatV(3,&(shelp1[0]));
            dev_Gamma3<RealT>(&(shelp1[0]));
            dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_complexT<RealT>(dev_mk3),&(shelp1[0]), &(ssum[0]));
         
         
//l==2,y        
            //positive direction
            hoppos = dev_nn[8*pos+2];
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref <RealT>(gf,4*pos+(2),&(gfsmem[ix]));
            #else
            dev_reconstructgf_2vtexref<RealT>(gf,4*pos+(2),&(gfsmem[ix]));
            #endif
            #ifdef USETEXTURE
              dev_su3MtV_spintex<RealT>(gfsmem[ix], hoppos, &(shelp1[0]));
            #else
              dev_su3MtV        <RealT>(gfsmem[ix], &(sin[6*hoppos]), &(shelp1[0]));
            #endif
            //-kappa(r - gamma_mu)
            dev_complexmult_add_assign_spinor(&(ssum[0]),dev_complexT<RealT>(dev_mk2),&(shelp1[0]), &(ssum[0]));
            //dev_GammatV(2,&(shelp1[0]));
            dev_Gamma2<RealT>(&(shelp1[0]));
            dev_complexmult_add_assign_spinor(&(ssum[0]),dev_complexT<RealT>(dev_k2),&(shelp1[0]), &(ssum[0]));
            
            //negative direction
            hoppos = dev_nn[8*pos+6];
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref_dagger <RealT>(gf,4*hoppos+(2),&(gfsmem[ix]));
            #else
            dev_reconstructgf_2vtexref_dagger<RealT>(gf,4*hoppos+(2),&(gfsmem[ix]));
            #endif
            #ifdef USETEXTURE
              dev_su3MtV_spintex<RealT>(gfsmem[ix], hoppos, &(shelp1[0]));
            #else
              dev_su3MtV        <RealT>(gfsmem[ix], &(sin[6*hoppos]), &(shelp1[0]));
            #endif
            //-kappa(r + gamma_mu)
            dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_complexT<RealT>(dev_mk2),&(shelp1[0]), &(ssum[0]));
            //dev_GammatV(2,&(shelp1[0]));
            dev_Gamma2<RealT>(&(shelp1[0]));
            dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_complexT<RealT>(dev_mk2),&(shelp1[0]), &(ssum[0]));


//l==1,x 
            //positive direction
            hoppos = dev_nn[8*pos+1];
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref <RealT>(gf,4*pos+(1),&(gfsmem[ix]));
            #else
            dev_reconstructgf_2vtexref<RealT>(gf,4*pos+(1),&(gfsmem[ix]));
            #endif
            #ifdef USETEXTURE
              dev_su3MtV_spintex<RealT>(gfsmem[ix], hoppos, &(shelp1[0]));
            #else
              dev_su3MtV        <RealT>(gfsmem[ix], &(sin[6*hoppos]), &(shelp1[0]));
            #endif
            //-kappa(r - gamma_mu)
            dev_complexmult_add_assign_spinor(&(ssum[0]),dev_complexT<RealT>(dev_mk1),&(shelp1[0]), &(ssum[0]));
            //dev_GammatV(1,&(shelp1[0]));
            dev_Gamma1<RealT>(&(shelp1[0]));
            dev_complexmult_add_assign_spinor(&(ssum[0]),dev_complexT<RealT>(dev_k1),&(shelp1[0]), &(ssum[0]));

            //negative direction
            hoppos = dev_nn[8*pos+5];
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref_dagger <RealT>(gf,4*hoppos+(1),&(gfsmem[ix]));
            #else
            dev_reconstructgf_2vtexref_dagger<RealT>(gf,4*hoppos+(1),&(gfsmem[ix]));
            #endif
            #ifdef USETEXTURE
              dev_su3MtV_spintex<RealT>(gfsmem[ix], hoppos, &(shelp1[0]));
            #else
              dev_su3MtV        <RealT>(gfsmem[ix], &(sin[6*hoppos]), &(shelp1[0]));
            #endif
            //-kappa(r + gamma_mu)
            dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_complexT<RealT>(dev_mk1),&(shelp1[0]), &(ssum[0]));
            //dev_GammatV(1,&(shelp1[0]));
            dev_Gamma1<RealT>(&(shelp1[0]));
            dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_complexT<RealT>(dev_mk1),&(shelp1[0]), &(ssum[0]));  
          
          
          
          //gamma5 term
         #ifdef USETEXTURE
          shelp1[0] = tex1Dfetch(spin_tex,6*pos);
          shelp1[1] = tex1Dfetch(spin_tex,6*pos+1);
          shelp1[2] = tex1Dfetch(spin_tex,6*pos+2);
          shelp1[3] = tex1Dfetch(spin_tex,6*pos+3);
          shelp1[4] = tex1Dfetch(spin_tex,6*pos+4);
          shelp1[5] = tex1Dfetch(spin_tex,6*pos+5);
         #else
          shelp1[0] = sin[6*pos];
          shelp1[1] = sin[6*pos+1];
          shelp1[2] = sin[6*pos+2];
          shelp1[3] = sin[6*pos+3];
          shelp1[4] = sin[6*pos+4];
          shelp1[5] = sin[6*pos+5];
         #endif 
          
          
          //dev_GammatV(4,&(shelp1[0]));
          dev_Gamma5<RealT>(&(shelp1[0]));
          dev_complexmult_add_assign_spinor<RealT>(&(ssum[0]),dev_initcomplex<RealT>(0.0,2.0*kappa*mu),&(shelp1[0]), &(sout[6*pos]));
  }
}





template<class RealT>
__global__ void dev_gamma5(typename dev_spinorT<RealT>::type * sin,typename dev_spinorT<RealT>::type * sout){
  int pos;
  pos= threadIdx.x + blockDim.x*blockIdx.x;  
  if(pos < dev_VOLUME){
          sout[6*pos+0].x = sin[6*pos+0].x;
          sout[6*pos+0].y = sin[6*pos+0].y;
          sout[6*pos+0].z = sin[6*pos+0].z;
          sout[6*pos+0].w = sin[6*pos+0].w;
          sout[6*pos+1].x = sin[6*pos+1].x;
          sout[6*pos+1].y = sin[6*pos+1].y;
          
          sout[6*pos+1].z = sin[6*pos+1].z;
          sout[6*pos+1].w = sin[6*pos+1].w;
          sout[6*pos+2].x = sin[6*pos+2].x;
          sout[6*pos+2].y = sin[6*pos+2].y;
          sout[6*pos+2].z = sin[6*pos+2].z;
          sout[6*pos+2].w = sin[6*pos+2].w;   
          
          sout[6*pos+3].x = -1.0*sin[6*pos+3].x;
          sout[6*pos+3].y = -1.0*sin[6*pos+3].y;
          sout[6*pos+3].z = -1.0*sin[6*pos+3].z;
          sout[6*pos+3].w = -1.0*sin[6*pos+3].w;
          sout[6*pos+4].x = -1.0*sin[6*pos+4].x;
          sout[6*pos+4].y = -1.0*sin[6*pos+4].y; 

          sout[6*pos+4].z = -1.0*sin[6*pos+4].z;
          sout[6*pos+4].w = -1.0*sin[6*pos+4].w;
          sout[6*pos+5].x = -1.0*sin[6*pos+5].x;
          sout[6*pos+5].y = -1.0*sin[6*pos+5].y;
          sout[6*pos+5].z = -1.0*sin[6*pos+5].z;
          sout[6*pos+5].w = -1.0*sin[6*pos+5].w;                 
  } 
}





template<class RealT>
void dev_tm_dirac_dagger_kappa
(
  typename dev_su3_2vT<RealT>::type * gf,
  typename dev_spinorT<RealT>::type* spinin,
  typename dev_spinorT<RealT>::type* spinout,
  int *grid, int * nn_grid, RealT* output,RealT* erg, int xsize, int ysize
){
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


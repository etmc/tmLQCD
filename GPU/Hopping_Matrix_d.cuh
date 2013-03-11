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
 * File: Hopping_Matrix_d.cuh
 *
 * CUDA Hopping_Matrix in double precision and associated functions
 *
 * 
 *
 **************************************************************************/



/*********** double constants on GPU *********************/
__device__ double mu_d;
__device__ double kappa_d;
__device__ double twokappamu_d;
__device__ double mubar_d, epsbar_d;

__device__ dev_complex_d dev_k0_d;
__device__ dev_complex_d dev_k1_d;
__device__ dev_complex_d dev_k2_d;
__device__ dev_complex_d dev_k3_d;

__device__ dev_complex_d dev_mk0_d;
__device__ dev_complex_d dev_mk1_d;
__device__ dev_complex_d dev_mk2_d;
__device__ dev_complex_d dev_mk3_d;

__constant__ __device__ dev_complex_d dev_k0c_d;
__constant__ __device__ dev_complex_d dev_k1c_d;
__constant__ __device__ dev_complex_d dev_k2c_d;
__constant__ __device__ dev_complex_d dev_k3c_d;

__constant__ __device__ dev_complex_d dev_mk0c_d;
__constant__ __device__ dev_complex_d dev_mk1c_d;
__constant__ __device__ dev_complex_d dev_mk2c_d;
__constant__ __device__ dev_complex_d dev_mk3c_d;



//-kappa(r - gamma_mu)
__device__ void dev_kappaP0_plus_d(dev_spinor_d * out, dev_spinor_d * in, dev_complex_d kappa){


     (*(out+0)).x -= (*(in+0)).x*kappa.re - (*(in+0)).y*kappa.im;
     (*(out+0)).y -= (*(in+0)).y*kappa.re + (*(in+0)).x*kappa.im;
     (*(out+0)).x -= (*(in+3)).x*kappa.re - (*(in+3)).y*kappa.im;
     (*(out+0)).y -= (*(in+3)).y*kappa.re + (*(in+3)).x*kappa.im;
     
     (*(out+3)).x -= (*(in+3)).x*kappa.re - (*(in+3)).y*kappa.im;
     (*(out+3)).y -= (*(in+3)).y*kappa.re + (*(in+3)).x*kappa.im;   
     (*(out+3)).x -= (*(in+0)).x*kappa.re - (*(in+0)).y*kappa.im;
     (*(out+3)).y -= (*(in+0)).y*kappa.re + (*(in+0)).x*kappa.im; 



     (*(out+0)).z -= (*(in+0)).z*kappa.re - (*(in+0)).w*kappa.im;
     (*(out+0)).w -= (*(in+0)).w*kappa.re + (*(in+0)).z*kappa.im;
     (*(out+0)).z -= (*(in+3)).z*kappa.re - (*(in+3)).w*kappa.im;
     (*(out+0)).w -= (*(in+3)).w*kappa.re + (*(in+3)).z*kappa.im;     
          
     (*(out+3)).z -= (*(in+3)).z*kappa.re - (*(in+3)).w*kappa.im;
     (*(out+3)).w -= (*(in+3)).w*kappa.re + (*(in+3)).z*kappa.im;
     (*(out+3)).z -= (*(in+0)).z*kappa.re - (*(in+0)).w*kappa.im;
     (*(out+3)).w -= (*(in+0)).w*kappa.re + (*(in+0)).z*kappa.im;

 
  
     (*(out+1)).x -= (*(in+1)).x*kappa.re - (*(in+1)).y*kappa.im;
     (*(out+1)).y -= (*(in+1)).y*kappa.re + (*(in+1)).x*kappa.im;
     (*(out+1)).x -= (*(in+4)).x*kappa.re - (*(in+4)).y*kappa.im;
     (*(out+1)).y -= (*(in+4)).y*kappa.re + (*(in+4)).x*kappa.im;     
     
     (*(out+4)).x -= (*(in+4)).x*kappa.re - (*(in+4)).y*kappa.im;
     (*(out+4)).y -= (*(in+4)).y*kappa.re + (*(in+4)).x*kappa.im;    
     (*(out+4)).x -= (*(in+1)).x*kappa.re - (*(in+1)).y*kappa.im;
     (*(out+4)).y -= (*(in+1)).y*kappa.re + (*(in+1)).x*kappa.im; 
 
 
 
     (*(out+1)).z -= (*(in+1)).z*kappa.re - (*(in+1)).w*kappa.im;
     (*(out+1)).w -= (*(in+1)).w*kappa.re + (*(in+1)).z*kappa.im;
     (*(out+1)).z -= (*(in+4)).z*kappa.re - (*(in+4)).w*kappa.im;
     (*(out+1)).w -= (*(in+4)).w*kappa.re + (*(in+4)).z*kappa.im;     
     
     (*(out+4)).z -= (*(in+4)).z*kappa.re - (*(in+4)).w*kappa.im;
     (*(out+4)).w -= (*(in+4)).w*kappa.re + (*(in+4)).z*kappa.im;    
     (*(out+4)).z -= (*(in+1)).z*kappa.re - (*(in+1)).w*kappa.im;
     (*(out+4)).w -= (*(in+1)).w*kappa.re + (*(in+1)).z*kappa.im;      
 

     
     (*(out+2)).x -= (*(in+2)).x*kappa.re - (*(in+2)).y*kappa.im;
     (*(out+2)).y -= (*(in+2)).y*kappa.re + (*(in+2)).x*kappa.im;
     (*(out+2)).x -= (*(in+5)).x*kappa.re - (*(in+5)).y*kappa.im;
     (*(out+2)).y -= (*(in+5)).y*kappa.re + (*(in+5)).x*kappa.im;
                
     (*(out+5)).x -= (*(in+5)).x*kappa.re - (*(in+5)).y*kappa.im;
     (*(out+5)).y -= (*(in+5)).y*kappa.re + (*(in+5)).x*kappa.im;   
     (*(out+5)).x -= (*(in+2)).x*kappa.re - (*(in+2)).y*kappa.im;
     (*(out+5)).y -= (*(in+2)).y*kappa.re + (*(in+2)).x*kappa.im; 
   
   
        
     (*(out+2)).z -= (*(in+2)).z*kappa.re - (*(in+2)).w*kappa.im;
     (*(out+2)).w -= (*(in+2)).w*kappa.re + (*(in+2)).z*kappa.im;
     (*(out+2)).z -= (*(in+5)).z*kappa.re - (*(in+5)).w*kappa.im;
     (*(out+2)).w -= (*(in+5)).w*kappa.re + (*(in+5)).z*kappa.im;
               
     (*(out+5)).z -= (*(in+5)).z*kappa.re - (*(in+5)).w*kappa.im;
     (*(out+5)).w -= (*(in+5)).w*kappa.re + (*(in+5)).z*kappa.im;   
     (*(out+5)).z -= (*(in+2)).z*kappa.re - (*(in+2)).w*kappa.im;
     (*(out+5)).w -= (*(in+2)).w*kappa.re + (*(in+2)).z*kappa.im; 
  
}






//-kappa(r + gamma_mu)
__device__ void dev_kappaP0_minus_d(dev_spinor_d * out, dev_spinor_d * in, dev_complex_d kappa){


     (*(out+0)).x -= (*(in+0)).x*kappa.re - (*(in+0)).y*kappa.im;
     (*(out+0)).y -= (*(in+0)).y*kappa.re + (*(in+0)).x*kappa.im;
     (*(out+0)).x += (*(in+3)).x*kappa.re - (*(in+3)).y*kappa.im;
     (*(out+0)).y += (*(in+3)).y*kappa.re + (*(in+3)).x*kappa.im;
     
     (*(out+3)).x -= (*(in+3)).x*kappa.re - (*(in+3)).y*kappa.im;
     (*(out+3)).y -= (*(in+3)).y*kappa.re + (*(in+3)).x*kappa.im;   
     (*(out+3)).x += (*(in+0)).x*kappa.re - (*(in+0)).y*kappa.im;
     (*(out+3)).y += (*(in+0)).y*kappa.re + (*(in+0)).x*kappa.im; 



     (*(out+0)).z -= (*(in+0)).z*kappa.re - (*(in+0)).w*kappa.im;
     (*(out+0)).w -= (*(in+0)).w*kappa.re + (*(in+0)).z*kappa.im;
     (*(out+0)).z += (*(in+3)).z*kappa.re - (*(in+3)).w*kappa.im;
     (*(out+0)).w += (*(in+3)).w*kappa.re + (*(in+3)).z*kappa.im;     
          
     (*(out+3)).z -= (*(in+3)).z*kappa.re - (*(in+3)).w*kappa.im;
     (*(out+3)).w -= (*(in+3)).w*kappa.re + (*(in+3)).z*kappa.im;
     (*(out+3)).z += (*(in+0)).z*kappa.re - (*(in+0)).w*kappa.im;
     (*(out+3)).w += (*(in+0)).w*kappa.re + (*(in+0)).z*kappa.im;

 
 
     (*(out+1)).x -= (*(in+1)).x*kappa.re - (*(in+1)).y*kappa.im;
     (*(out+1)).y -= (*(in+1)).y*kappa.re + (*(in+1)).x*kappa.im;
     (*(out+1)).x += (*(in+4)).x*kappa.re - (*(in+4)).y*kappa.im;
     (*(out+1)).y += (*(in+4)).y*kappa.re + (*(in+4)).x*kappa.im;     
     
     (*(out+4)).x -= (*(in+4)).x*kappa.re - (*(in+4)).y*kappa.im;
     (*(out+4)).y -= (*(in+4)).y*kappa.re + (*(in+4)).x*kappa.im;    
     (*(out+4)).x += (*(in+1)).x*kappa.re - (*(in+1)).y*kappa.im;
     (*(out+4)).y += (*(in+1)).y*kappa.re + (*(in+1)).x*kappa.im; 
 
 
 
     (*(out+1)).z -= (*(in+1)).z*kappa.re - (*(in+1)).w*kappa.im;
     (*(out+1)).w -= (*(in+1)).w*kappa.re + (*(in+1)).z*kappa.im;
     (*(out+1)).z += (*(in+4)).z*kappa.re - (*(in+4)).w*kappa.im;
     (*(out+1)).w += (*(in+4)).w*kappa.re + (*(in+4)).z*kappa.im;     
     
     (*(out+4)).z -= (*(in+4)).z*kappa.re - (*(in+4)).w*kappa.im;
     (*(out+4)).w -= (*(in+4)).w*kappa.re + (*(in+4)).z*kappa.im;    
     (*(out+4)).z += (*(in+1)).z*kappa.re - (*(in+1)).w*kappa.im;
     (*(out+4)).w += (*(in+1)).w*kappa.re + (*(in+1)).z*kappa.im;      
 

     
     (*(out+2)).x -= (*(in+2)).x*kappa.re - (*(in+2)).y*kappa.im;
     (*(out+2)).y -= (*(in+2)).y*kappa.re + (*(in+2)).x*kappa.im;
     (*(out+2)).x += (*(in+5)).x*kappa.re - (*(in+5)).y*kappa.im;
     (*(out+2)).y += (*(in+5)).y*kappa.re + (*(in+5)).x*kappa.im;
                
     (*(out+5)).x -= (*(in+5)).x*kappa.re - (*(in+5)).y*kappa.im;
     (*(out+5)).y -= (*(in+5)).y*kappa.re + (*(in+5)).x*kappa.im;   
     (*(out+5)).x += (*(in+2)).x*kappa.re - (*(in+2)).y*kappa.im;
     (*(out+5)).y += (*(in+2)).y*kappa.re + (*(in+2)).x*kappa.im; 
   
   
       
     (*(out+2)).z -= (*(in+2)).z*kappa.re - (*(in+2)).w*kappa.im;
     (*(out+2)).w -= (*(in+2)).w*kappa.re + (*(in+2)).z*kappa.im;
     (*(out+2)).z += (*(in+5)).z*kappa.re - (*(in+5)).w*kappa.im;
     (*(out+2)).w += (*(in+5)).w*kappa.re + (*(in+5)).z*kappa.im;
               
     (*(out+5)).z -= (*(in+5)).z*kappa.re - (*(in+5)).w*kappa.im;
     (*(out+5)).w -= (*(in+5)).w*kappa.re + (*(in+5)).z*kappa.im;   
     (*(out+5)).z += (*(in+2)).z*kappa.re - (*(in+2)).w*kappa.im;
     (*(out+5)).w += (*(in+2)).w*kappa.re + (*(in+2)).z*kappa.im; 
  
}







#ifdef RELATIVISTIC_BASIS
//  here comes P0+- for the relativistic basis
//  in this basis we have:
//
//  gamma0 =
//  -1  0  0  0 
//   0 -1  0  0
//   0  0  1  0
//   0  0  0  1
//
//



//-kappa(r - gamma_mu)
__device__ void dev_kappaP0_plus_relativistic_d(dev_spinor_d * out, dev_spinor_d * in, dev_complex_d kappa){


     (*(out+0)).x -= 2.0*( (*(in+0)).x*kappa.re - (*(in+0)).y*kappa.im );
     (*(out+0)).y -= 2.0*( (*(in+0)).y*kappa.re + (*(in+0)).x*kappa.im );
     (*(out+0)).z -= 2.0*( (*(in+0)).z*kappa.re - (*(in+0)).w*kappa.im );
     (*(out+0)).w -= 2.0*( (*(in+0)).w*kappa.re + (*(in+0)).z*kappa.im );
     
     (*(out+1)).x -= 2.0*( (*(in+1)).x*kappa.re - (*(in+1)).y*kappa.im );
     (*(out+1)).y -= 2.0*( (*(in+1)).y*kappa.re + (*(in+1)).x*kappa.im );     
     (*(out+1)).z -= 2.0*( (*(in+1)).z*kappa.re - (*(in+1)).w*kappa.im );
     (*(out+1)).w -= 2.0*( (*(in+1)).w*kappa.re + (*(in+1)).z*kappa.im );     
     
     (*(out+2)).x -= 2.0*( (*(in+2)).x*kappa.re - (*(in+2)).y*kappa.im );
     (*(out+2)).y -= 2.0*( (*(in+2)).y*kappa.re + (*(in+2)).x*kappa.im );
     (*(out+2)).z -= 2.0*( (*(in+2)).z*kappa.re - (*(in+2)).w*kappa.im );
     (*(out+2)).w -= 2.0*( (*(in+2)).w*kappa.re + (*(in+2)).z*kappa.im );

}






//-kappa(r + gamma_mu)
__device__ void dev_kappaP0_minus_relativistic_d(dev_spinor_d * out, dev_spinor_d * in, dev_complex_d kappa){

    
     (*(out+3)).x -= 2.0*( (*(in+3)).x*kappa.re - (*(in+3)).y*kappa.im );
     (*(out+3)).y -= 2.0*( (*(in+3)).y*kappa.re + (*(in+3)).x*kappa.im );     
     (*(out+3)).z -= 2.0*( (*(in+3)).z*kappa.re - (*(in+3)).w*kappa.im );
     (*(out+3)).w -= 2.0*( (*(in+3)).w*kappa.re + (*(in+3)).z*kappa.im );


     (*(out+4)).x -= 2.0*( (*(in+4)).x*kappa.re - (*(in+4)).y*kappa.im );
     (*(out+4)).y -= 2.0*( (*(in+4)).y*kappa.re + (*(in+4)).x*kappa.im );     
     (*(out+4)).z -= 2.0*( (*(in+4)).z*kappa.re - (*(in+4)).w*kappa.im );
     (*(out+4)).w -= 2.0*( (*(in+4)).w*kappa.re + (*(in+4)).z*kappa.im );    

           
     (*(out+5)).x -= 2.0*( (*(in+5)).x*kappa.re - (*(in+5)).y*kappa.im );
     (*(out+5)).y -= 2.0*( (*(in+5)).y*kappa.re + (*(in+5)).x*kappa.im );            
     (*(out+5)).z -= 2.0*( (*(in+5)).z*kappa.re - (*(in+5)).w*kappa.im );
     (*(out+5)).w -= 2.0*( (*(in+5)).w*kappa.re + (*(in+5)).z*kappa.im );   

  
}

//RELATIVISTIC_BASIS
#endif 







//applies the Hopping Part Even-Odd !
//the gauge field is the complete gaugefield!
//the gauge field at the local point is reconstructed by 2*pos+eo where pos is the eo-position
//from 0..VOLUME/2-1, eo = 0 or 1
//the positions in the gauge fields are passed in "gfindex_site" for gf's that are attached at
//the actual positions and in "gfindex_nextsite" for gf's that start at a position of the 
//other eo-sublattice.
//for the hopping positions of the eo-spinor field we use on of the two dedicated eo-nn fields
//the boundary conditions are implemented as in Hopping_Matrix.c
//mult with complex conjugate k0,k1,k2,k3 in positive direction because
// psi(x+mu) != exp(i theta_mu) psi(x)  
// we start from site index start and go to start+size
// by this we can split up bulk and rand in mpi
__global__ void dev_Hopping_Matrix_d(dev_su3_2v_d * gf, const dev_spinor_d * sin, dev_spinor_d * sout, const int * gfindex_site, const int* gfindex_nextsite, const int * nn_evenodd, const int eo, int start, int size){

  int pos,hoppos;
    dev_spinor_d shelp1[6], ssum[6];
    __shared__ dev_su3_d gfsmem[BLOCKD];


  pos = start  +  threadIdx.x + blockDim.x * blockIdx.x;
  int ix = threadIdx.x;
  
  
  if(pos < (start + size)){
  

  dev_zero_spinor_local_d(&(ssum[0])); // zero sum        
 #ifdef TEMPORALGAUGE
  int spatialvol = dev_LX*dev_LY*dev_LZ;
 #endif
  



//hopping term                
//l==0,t
            //positive direction
            hoppos = nn_evenodd[8*pos];
             //hoppos = tex1Dfetch(nn_tex,8*pos);
            //color
            
            #ifdef TEMPORALGAUGE
              // gf == ID for t != T-1 => just read the spinor
              #ifdef MPI
                if ( ((gfindex_site[pos]) < (dev_T-1)*spatialvol) || (dev_rank < dev_nproc-1) ) {
                //if ((gfindex_site[pos]) < (dev_T-1)*spatialvol) { // FAKE TEMPORALGAUGE
              #else
                if ((gfindex_site[pos]/spatialvol) != (dev_T-1) ) {
              #endif
                shelp1[0] = sin[hoppos+0*DEVOFF];
                shelp1[1] = sin[hoppos+1*DEVOFF];
                shelp1[2] = sin[hoppos+2*DEVOFF];
		shelp1[3] = sin[hoppos+3*DEVOFF];
                shelp1[4] = sin[hoppos+4*DEVOFF];
                shelp1[5] = sin[hoppos+5*DEVOFF];
              }
              else{
                // gf != ID for t == T-1 => mult spinor with gf
                dev_reconstructgf_2vtexref_d(gf,4*gfindex_site[pos]+0,&(gfsmem[ix]));
                dev_su3MtV_d(gfsmem[ix], &(sin[hoppos]), &(shelp1[0]));

              }
            #else
              dev_reconstructgf_2vtexref_d(gf, 4*gfindex_site[pos]+ 0 ,&(gfsmem[ix]));
              dev_su3MtV_d(gfsmem[ix], &(sin[hoppos]), &(shelp1[0]));
            #endif
            //-kappa(r - gamma_mu)
              dev_kappaP0_plus_d(&(ssum[0]), &(shelp1[0]), dev_cconj_d(dev_k0_d));
	    
//l==0,t
            //negative direction
            hoppos = nn_evenodd[8*pos+4]; 
             //hoppos = tex1Dfetch(nn_tex,8*pos+4);
            //color
            #ifdef TEMPORALGAUGE
              // gf == ID for t != T-1 => just read the spinor
              #ifdef MPI
                if ( ((gfindex_nextsite[hoppos]) < (dev_T-1)*spatialvol) || (dev_rank > 0) ) {
                //if ((gfindex_nextsite[hoppos]) < (dev_T-1)*spatialvol) { // FAKE TEMPORALGAUGE
              #else
                if ((gfindex_nextsite[hoppos]/spatialvol) != (dev_T-1) ) {
              #endif
              
                shelp1[0] = sin[hoppos+0*DEVOFF];
                shelp1[1] = sin[hoppos+1*DEVOFF];
                shelp1[2] = sin[hoppos+2*DEVOFF];
		shelp1[3] = sin[hoppos+3*DEVOFF];
                shelp1[4] = sin[hoppos+4*DEVOFF];
                shelp1[5] = sin[hoppos+5*DEVOFF];

              }
              else{
                // gf != ID for t == T-1 => mult spinor with gf
                dev_reconstructgf_2vtexref_dagger_d(gf,4*gfindex_nextsite[hoppos]+0, &(gfsmem[ix]));
                dev_su3MtV_d(gfsmem[ix], &(sin[hoppos]), &(shelp1[0]));
              }
            #else            
              dev_reconstructgf_2vtexref_dagger_d(gf,4*gfindex_nextsite[hoppos]+0, &(gfsmem[ix]));
              dev_su3MtV_d(gfsmem[ix], &(sin[hoppos]), &(shelp1[0]));
            #endif
            
            //-kappa(r + gamma_mu)
              dev_kappaP0_minus_d(&(ssum[0]), &(shelp1[0]), dev_k0_d);



//l==3,z 
            //positive direction
            hoppos = nn_evenodd[8*pos+3];
             //hoppos = tex1Dfetch(nn_tex,8*pos+3);
            //color
            dev_reconstructgf_2vtexref_d(gf, 4*gfindex_site[pos]+3, &(gfsmem[ix]));
            dev_su3MtV_kappaP3_plus_d(gfsmem[ix],&(sin[hoppos]), &(ssum[0]), dev_k3_d.re);

//l==3,z               
            //negative direction
            hoppos = nn_evenodd[8*pos+7];
             //hoppos = tex1Dfetch(nn_tex,8*pos+7); 
            //color
            dev_reconstructgf_2vtexref_dagger_d(gf,4*gfindex_nextsite[hoppos]+3, &(gfsmem[ix]));
	    dev_su3MtV_kappaP3_minus_d(gfsmem[ix],&(sin[hoppos]), &(ssum[0]), dev_k3_d.re);


	    
//l==2,y 
            //positive direction
            hoppos = nn_evenodd[8*pos+2];
             //hoppos = tex1Dfetch(nn_tex,8*pos+2);
            //color
            dev_reconstructgf_2vtexref_d(gf,4*gfindex_site[pos]+2, &(gfsmem[ix]));
	    dev_su3MtV_kappaP2_plus_d(gfsmem[ix],&(sin[hoppos]), &(ssum[0]), dev_k2_d.re);

//l==2,y                  
            //negative direction
            hoppos = nn_evenodd[8*pos+6]; 
             //hoppos = tex1Dfetch(nn_tex,8*pos+6);
            //color
            dev_reconstructgf_2vtexref_dagger_d(gf,4*gfindex_nextsite[hoppos]+2, &(gfsmem[ix]));
	    dev_su3MtV_kappaP2_minus_d(gfsmem[ix],&(sin[hoppos]), &(ssum[0]), dev_k2_d.re);
            

	    
//l==1,x 
            //positive direction
            hoppos = nn_evenodd[8*pos+1];
             //hoppos = tex1Dfetch(nn_tex,8*pos+1);
            //color
            dev_reconstructgf_2vtexref_d(gf,4*gfindex_site[pos]+1, &(gfsmem[ix]));
	    dev_su3MtV_kappaP1_plus_d(gfsmem[ix],&(sin[hoppos]), &(ssum[0]), dev_k1_d.re);

//l==1,x 
            //negative direction
            hoppos = nn_evenodd[8*pos+5]; 
             //hoppos = tex1Dfetch(nn_tex,8*pos+5);
            //color
            dev_reconstructgf_2vtexref_dagger_d(gf,4*gfindex_nextsite[hoppos]+1, &(gfsmem[ix]));
	    dev_su3MtV_kappaP1_minus_d(gfsmem[ix],&(sin[hoppos]), &(ssum[0]), dev_k1_d.re);
            
          
        //copy to output spinor
        dev_write_spinor_d(&(ssum[0]),&(sout[pos])); 
  }
}






__global__ void dev_Hopping_Matrix_updn_d(dev_su3_2v_d * gf, const dev_spinor_d * sin_up, const dev_spinor_d * sin_dn, dev_spinor_d * sout_up, dev_spinor_d * sout_dn, const int * gfindex_site, const int* gfindex_nextsite, const int * nn_evenodd, const int eo, int start, int size){

  int pos,hoppos;
    dev_spinor_d shelp1_up[6], ssum_up[6];
    dev_spinor_d shelp1_dn[6], ssum_dn[6];    
    __shared__ dev_su3_d gfsmem[BLOCKD];


  pos = start  +  threadIdx.x + blockDim.x * blockIdx.x;
  int ix = threadIdx.x;
  
  
  if(pos < (start + size)){
  

  dev_zero_spinor_local_d(&(ssum_up[0])); // zero sum 
  dev_zero_spinor_local_d(&(ssum_dn[0])); // zero sum   
 #ifdef TEMPORALGAUGE
  int spatialvol = dev_LX*dev_LY*dev_LZ;
 #endif
  



//hopping term                
//l==0,t
            //positive direction
            hoppos = nn_evenodd[8*pos];
             //hoppos = tex1Dfetch(nn_tex,8*pos);
            //color
            
            #ifdef TEMPORALGAUGE
              // gf == ID for t != T-1 => just read the spinor
              #ifdef MPI
                if ( ((gfindex_site[pos]) < (dev_T-1)*spatialvol) || (dev_rank < dev_nproc-1) ) {
                //if ((gfindex_site[pos]) < (dev_T-1)*spatialvol) { // FAKE TEMPORALGAUGE
              #else
                if ((gfindex_site[pos]/spatialvol) != (dev_T-1) ) {
              #endif
                shelp1_up[0] = sin_up[hoppos+0*DEVOFF];
                shelp1_up[1] = sin_up[hoppos+1*DEVOFF];
                shelp1_up[2] = sin_up[hoppos+2*DEVOFF];
		shelp1_up[3] = sin_up[hoppos+3*DEVOFF];
                shelp1_up[4] = sin_up[hoppos+4*DEVOFF];
                shelp1_up[5] = sin_up[hoppos+5*DEVOFF];
                shelp1_dn[0] = sin_dn[hoppos+0*DEVOFF];
                shelp1_dn[1] = sin_dn[hoppos+1*DEVOFF];
                shelp1_dn[2] = sin_dn[hoppos+2*DEVOFF];
		shelp1_dn[3] = sin_dn[hoppos+3*DEVOFF];
                shelp1_dn[4] = sin_dn[hoppos+4*DEVOFF];
                shelp1_dn[5] = sin_dn[hoppos+5*DEVOFF];		
              }
              else{
                // gf != ID for t == T-1 => mult spinor with gf
                dev_reconstructgf_2vtexref_d(gf,4*gfindex_site[pos]+0,&(gfsmem[ix]));
                dev_su3MtV_d(gfsmem[ix], &(sin_up[hoppos]), &(shelp1_up[0]));
                dev_su3MtV_d(gfsmem[ix], &(sin_dn[hoppos]), &(shelp1_dn[0]));
              }
            #else
              dev_reconstructgf_2vtexref_d(gf, 4*gfindex_site[pos]+ 0 ,&(gfsmem[ix]));
              dev_su3MtV_d(gfsmem[ix], &(sin_up[hoppos]), &(shelp1_up[0]));
              dev_su3MtV_d(gfsmem[ix], &(sin_dn[hoppos]), &(shelp1_dn[0]));	      
            #endif
            //-kappa(r - gamma_mu)
              dev_kappaP0_plus_d(&(ssum_up[0]), &(shelp1_up[0]), dev_cconj_d(dev_k0_d));
              dev_kappaP0_plus_d(&(ssum_dn[0]), &(shelp1_dn[0]), dev_cconj_d(dev_k0_d));	    
//l==0,t
            //negative direction
            hoppos = nn_evenodd[8*pos+4]; 
             //hoppos = tex1Dfetch(nn_tex,8*pos+4);
            //color
            #ifdef TEMPORALGAUGE
              // gf == ID for t != T-1 => just read the spinor
              #ifdef MPI
                if ( ((gfindex_nextsite[hoppos]) < (dev_T-1)*spatialvol) || (dev_rank > 0) ) {
                //if ((gfindex_nextsite[hoppos]) < (dev_T-1)*spatialvol) { // FAKE TEMPORALGAUGE
              #else
                if ((gfindex_nextsite[hoppos]/spatialvol) != (dev_T-1) ) {
              #endif
              
                shelp1_up[0] = sin_up[hoppos+0*DEVOFF];
                shelp1_up[1] = sin_up[hoppos+1*DEVOFF];
                shelp1_up[2] = sin_up[hoppos+2*DEVOFF];
		shelp1_up[3] = sin_up[hoppos+3*DEVOFF];
                shelp1_up[4] = sin_up[hoppos+4*DEVOFF];
                shelp1_up[5] = sin_up[hoppos+5*DEVOFF];
                shelp1_dn[0] = sin_dn[hoppos+0*DEVOFF];
                shelp1_dn[1] = sin_dn[hoppos+1*DEVOFF];
                shelp1_dn[2] = sin_dn[hoppos+2*DEVOFF];
		shelp1_dn[3] = sin_dn[hoppos+3*DEVOFF];
                shelp1_dn[4] = sin_dn[hoppos+4*DEVOFF];
                shelp1_dn[5] = sin_dn[hoppos+5*DEVOFF];
              }
              else{
                // gf != ID for t == T-1 => mult spinor with gf
                dev_reconstructgf_2vtexref_dagger_d(gf,4*gfindex_nextsite[hoppos]+0, &(gfsmem[ix]));
                dev_su3MtV_d(gfsmem[ix], &(sin_up[hoppos]), &(shelp1_up[0]));
                dev_su3MtV_d(gfsmem[ix], &(sin_dn[hoppos]), &(shelp1_dn[0]));		
              }
            #else            
              dev_reconstructgf_2vtexref_dagger_d(gf,4*gfindex_nextsite[hoppos]+0, &(gfsmem[ix]));
              dev_su3MtV_d(gfsmem[ix], &(sin_up[hoppos]), &(shelp1_up[0]));
              dev_su3MtV_d(gfsmem[ix], &(sin_dn[hoppos]), &(shelp1_dn[0]));	      
            #endif
            
            //-kappa(r + gamma_mu)
              dev_kappaP0_minus_d(&(ssum_up[0]), &(shelp1_up[0]), dev_k0_d);
              dev_kappaP0_minus_d(&(ssum_dn[0]), &(shelp1_dn[0]), dev_k0_d);


//l==3,z 
            //positive direction
            hoppos = nn_evenodd[8*pos+3];
             //hoppos = tex1Dfetch(nn_tex,8*pos+3);
            //color
            dev_reconstructgf_2vtexref_d(gf, 4*gfindex_site[pos]+3, &(gfsmem[ix]));
            dev_su3MtV_kappaP3_plus_d(gfsmem[ix],&(sin_up[hoppos]), &(ssum_up[0]), dev_k3_d.re);
            dev_su3MtV_kappaP3_plus_d(gfsmem[ix],&(sin_dn[hoppos]), &(ssum_dn[0]), dev_k3_d.re);
//l==3,z               
            //negative direction
            hoppos = nn_evenodd[8*pos+7];
             //hoppos = tex1Dfetch(nn_tex,8*pos+7); 
            //color
            dev_reconstructgf_2vtexref_dagger_d(gf,4*gfindex_nextsite[hoppos]+3, &(gfsmem[ix]));
	    dev_su3MtV_kappaP3_minus_d(gfsmem[ix],&(sin_up[hoppos]), &(ssum_up[0]), dev_k3_d.re);
	    dev_su3MtV_kappaP3_minus_d(gfsmem[ix],&(sin_dn[hoppos]), &(ssum_dn[0]), dev_k3_d.re);

	    
//l==2,y 
            //positive direction
            hoppos = nn_evenodd[8*pos+2];
             //hoppos = tex1Dfetch(nn_tex,8*pos+2);
            //color
            dev_reconstructgf_2vtexref_d(gf,4*gfindex_site[pos]+2, &(gfsmem[ix]));
	    dev_su3MtV_kappaP2_plus_d(gfsmem[ix],&(sin_up[hoppos]), &(ssum_up[0]), dev_k2_d.re);
	    dev_su3MtV_kappaP2_plus_d(gfsmem[ix],&(sin_dn[hoppos]), &(ssum_dn[0]), dev_k2_d.re);
//l==2,y                  
            //negative direction
            hoppos = nn_evenodd[8*pos+6]; 
             //hoppos = tex1Dfetch(nn_tex,8*pos+6);
            //color
            dev_reconstructgf_2vtexref_dagger_d(gf,4*gfindex_nextsite[hoppos]+2, &(gfsmem[ix]));
	    dev_su3MtV_kappaP2_minus_d(gfsmem[ix],&(sin_up[hoppos]), &(ssum_up[0]), dev_k2_d.re);
	    dev_su3MtV_kappaP2_minus_d(gfsmem[ix],&(sin_dn[hoppos]), &(ssum_dn[0]), dev_k2_d.re);            

	    
//l==1,x 
            //positive direction
            hoppos = nn_evenodd[8*pos+1];
             //hoppos = tex1Dfetch(nn_tex,8*pos+1);
            //color
            dev_reconstructgf_2vtexref_d(gf,4*gfindex_site[pos]+1, &(gfsmem[ix]));
	    dev_su3MtV_kappaP1_plus_d(gfsmem[ix],&(sin_up[hoppos]), &(ssum_up[0]), dev_k1_d.re);
	    dev_su3MtV_kappaP1_plus_d(gfsmem[ix],&(sin_dn[hoppos]), &(ssum_dn[0]), dev_k1_d.re);
//l==1,x 
            //negative direction
            hoppos = nn_evenodd[8*pos+5]; 
             //hoppos = tex1Dfetch(nn_tex,8*pos+5);
            //color
            dev_reconstructgf_2vtexref_dagger_d(gf,4*gfindex_nextsite[hoppos]+1, &(gfsmem[ix]));
	    dev_su3MtV_kappaP1_minus_d(gfsmem[ix],&(sin_up[hoppos]), &(ssum_up[0]), dev_k1_d.re);
	    dev_su3MtV_kappaP1_minus_d(gfsmem[ix],&(sin_dn[hoppos]), &(ssum_dn[0]), dev_k1_d.re);            
          
        //copy to output spinor
        dev_write_spinor_d(&(ssum_up[0]),&(sout_up[pos])); 
        dev_write_spinor_d(&(ssum_dn[0]),&(sout_dn[pos])); 	
  }
}





// computes sout = 1/(1 +- mutilde gamma5) sin = (1 -+ i mutilde gamma5)/(1+mutilde^2) sin
// mutilde = 2 kappa mu
// uses shared local memory for manipulation
__global__ void dev_mul_one_pm_imu_inv_d(dev_spinor_d* sin, dev_spinor_d* sout, const double sign){
   
   dev_spinor_d slocal[6];
   //need the inverse sign in the numerator because of inverse
   dev_complex_d pm_imu = dev_initcomplex_d(0.0,-1.0*sign*twokappamu_d);
   
   double one_plus_musquare_inv = 1.0/(1.0 + twokappamu_d*twokappamu_d);
   int pos;
   pos= threadIdx.x + blockDim.x*blockIdx.x;  

   if(pos < dev_VOLUME){
     dev_skalarmult_gamma5_globalspinor_d(&(slocal[0]), pm_imu, &(sin[pos]) );
     dev_add_globalspinor_assign_d(&(slocal[0]), &(sin[pos])); 
     dev_realmult_spinor_assigntoglobal_d(&(sout[pos]), one_plus_musquare_inv, &(slocal[0]) );
   }
}





// sout = gamma_5*((1\pm i\mutilde \gamma_5)*sin1 - sin2)
// uses shared local memory for manipulation
__global__ void dev_mul_one_pm_imu_sub_mul_gamma5_d(dev_spinor_d* sin1, dev_spinor_d* sin2, dev_spinor_d* sout, const double sign){
   dev_spinor_d slocal[6];
   dev_complex_d pm_imu = dev_initcomplex_d(0.0, sign*twokappamu_d); // i mutilde
   int pos;
   pos= threadIdx.x + blockDim.x*blockIdx.x; 

   if(pos < dev_VOLUME){
     dev_skalarmult_gamma5_globalspinor_d(&(slocal[0]),pm_imu,&(sin1[pos]));
     dev_add_globalspinor_assign_d(&(slocal[0]), &(sin1[pos]));
     dev_sub_globalspinor_assign_d(&(slocal[0]), &(sin2[pos]));
     dev_Gamma5_assigntoglobal_d(&(sout[pos]), &(slocal[0]));
   }   
}



















// aequivalent to Qtm_pm_psi in tm_operators.c, this is NON-MPI version
extern "C" void dev_Qtm_pm_psi_d(dev_spinor_d* spinin, dev_spinor_d* spinout, 
				 dev_spinor_d* spin_eo1_d, dev_spinor_d* spin_eo2_d, 
				 int gridsize, int blocksize, int gridsize2, int blocksize2,
				 int* dev_eoidx_even, int* dev_eoidx_odd, 
				 int* dev_nn_eo, int* dev_nn_oe){
  //spinin == odd
  //spinout == odd
  cudaError_t cudaerr;
  int VolumeEO = VOLUME/2;
  //Q_{-}
  //cudaFuncSetCacheConfig(dev_Hopping_Matrix, cudaFuncCachePreferL1);
  dev_Hopping_Matrix_d<<<gridsize, blocksize>>>
             (dev_gf_d, spinin, spin_eo1_d, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0, 0, VolumeEO); //spin_eo1 == even -> 0           
  dev_mul_one_pm_imu_inv_d<<<gridsize2, blocksize2>>>(spin_eo1_d,spin_eo2_d, -1.);
 

  //cudaFuncSetCacheConfig(dev_Hopping_Matrix, cudaFuncCachePreferL1);
  dev_Hopping_Matrix_d<<<gridsize, blocksize>>>
            (dev_gf_d, spin_eo2_d, spin_eo1_d, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1, 0, VolumeEO); 

  dev_mul_one_pm_imu_sub_mul_gamma5_d<<<gridsize2, blocksize2>>>(spinin, spin_eo1_d,  spin_eo2_d, -1.);

  //Q_{+}
  //cudaFuncSetCacheConfig(dev_Hopping_Matrix, cudaFuncCachePreferL1);
  dev_Hopping_Matrix_d<<<gridsize, blocksize>>>
          (dev_gf_d, spin_eo2_d, spin_eo1_d, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0, 0, VolumeEO); //spin_eo1 == even -> 0
  dev_mul_one_pm_imu_inv_d<<<gridsize2, blocksize2>>>(spin_eo1_d,spinout, +1.);
 

  //cudaFuncSetCacheConfig(dev_Hopping_Matrix, cudaFuncCachePreferL1);
  dev_Hopping_Matrix_d<<<gridsize, blocksize>>>
             (dev_gf_d, spinout, spin_eo1_d, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1, 0, VolumeEO); 
  dev_mul_one_pm_imu_sub_mul_gamma5_d<<<gridsize2, blocksize2>>>(spin_eo2_d, spin_eo1_d,  spinout , +1.); 

  if ((cudaerr=cudaGetLastError())!=cudaSuccess) {
    if (g_proc_id == 0){
      printf("Error in dev_Qtm_pm_psi_d: %s\n", cudaGetErrorString(cudaGetLastError()));
      printf("Error code is: %f\n",cudaerr);
    }
    exit(100);
  }
}



// aequivalent to H_eo_tm_inv_psi in tm_operators.c, this is NON-MPI version
extern "C" void dev_H_eo_tm_inv_psi_d(dev_spinor_d* spinin, dev_spinor_d* spinout,dev_spinor_d* spin_eo1_d, 
				 int gridsize, int blocksize, int gridsize2, int blocksize2,
				 int* dev_eoidx_even, int* dev_eoidx_odd, 
				 int* dev_nn_eo, int* dev_nn_oe, const int ieo, const double sign){
 
  cudaError_t cudaerr;
  int VolumeEO = VOLUME/2;
  
  if(ieo==0){
  //Q_{-}
  //cudaFuncSetCacheConfig(dev_Hopping_Matrix, cudaFuncCachePreferL1);
  dev_Hopping_Matrix_d<<<gridsize, blocksize>>>
             (dev_gf_d, spinin, spin_eo1_d, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0, 0, VolumeEO); //spin_eo1 == even -> 0           
  dev_mul_one_pm_imu_inv_d<<<gridsize2, blocksize2>>>(spin_eo1_d,spinout, sign);
  }
  
  if(ieo==1){
  //cudaFuncSetCacheConfig(dev_Hopping_Matrix, cudaFuncCachePreferL1);
  dev_Hopping_Matrix_d<<<gridsize, blocksize>>>
            (dev_gf_d, spinin, spin_eo1_d, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1, 0, VolumeEO);  //spin_eo1 == odd -> 1 
  dev_mul_one_pm_imu_inv_d<<<gridsize2, blocksize2>>>(spin_eo1_d,spinout, sign);
  }


  if ((cudaerr=cudaGetLastError())!=cudaSuccess) {
    if (g_proc_id == 0){
      printf("Error in dev_H_eo_tm_inv_psi_d: %s\n", cudaGetErrorString(cudaGetLastError()));
      printf("Error code is: %f\n",cudaerr);
    }
    exit(100);
  }
}





// aequivalent to Qtm_minus_psi in tm_operators.c, this is NON-MPI version
extern "C" void dev_Qtm_minus_psi_d(dev_spinor_d* spinin, dev_spinor_d* spinout,dev_spinor_d* spin_eo1_d, 
				 int gridsize, int blocksize, int gridsize2, int blocksize2,
				 int* dev_eoidx_even, int* dev_eoidx_odd, 
				 int* dev_nn_eo, int* dev_nn_oe){
 
  cudaError_t cudaerr;
  int VolumeEO = VOLUME/2;
  
  dev_H_eo_tm_inv_psi_d(spinin, spinout , spin_eo1_d, 
				 gridsize, blocksize, gridsize2, blocksize2,
				 dev_eoidx_even, dev_eoidx_odd, 
				 dev_nn_eo, dev_nn_oe, 0, -1);

  dev_Hopping_Matrix_d<<<gridsize, blocksize>>>
            (dev_gf_d, spinout, spin_eo1_d, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1, 0, VolumeEO);
	    
  dev_mul_one_pm_imu_sub_mul_gamma5_d<<<gridsize2, blocksize2>>>(spinin, spin_eo1_d,  spinout , -1.); 
	    

  if ((cudaerr=cudaGetLastError())!=cudaSuccess) {
    if (g_proc_id == 0){
      printf("Error in dev_Qtm_minus_psi_d: %s\n", cudaGetErrorString(cudaGetLastError()));
      printf("Error code is: %f\n",cudaerr);
    }
    exit(100);
  }
}



// void H_eo_tm_inv_psi(spinor * const l, spinor * const k, 
// 		     const int ieo, const double _sign) {
// 
//   Hopping_Matrix(ieo, l, k);
//   mul_one_pm_imu_inv(l, _sign, VOLUME/2);
// }  
// void Qtm_minus_psi(spinor * const l, spinor * const k) {
//   H_eo_tm_inv_psi(g_spinor_field[DUM_MATRIX+1], k, EO, -1);
//   Hopping_Matrix(OE, g_spinor_field[DUM_MATRIX+2], g_spinor_field[DUM_MATRIX+1]);
//   mul_one_pm_imu_sub_mul_gamma5(l, k, g_spinor_field[DUM_MATRIX+2], -1);
//   //tm_sub_H_eo_gamma5(l, k, g_spinor_field[DUM_MATRIX+1], OE, -1.);
// }



// derived from function  dev_mul_one_pm_imu_inv
// applies (1 +- imubar*gamma5)
// one thread per lattice site
__global__ void dev_mul_one_pm_imubar_gamma5_d (dev_spinor_d * sin,
                                              dev_spinor_d * sout,
                                              double sign         ) {
   
  dev_spinor_d slocal[6];
  dev_complex_d pm_imu = dev_initcomplex_d(0.0, sign * mubar_d);
  int pos = threadIdx.x + blockDim.x*blockIdx.x;
  
  if (pos < dev_VOLUME) {
    dev_skalarmult_gamma5_globalspinor_d(&(slocal[0]), pm_imu, &(sin[pos]) );			// slocal  =  pm_imu * (gamma5) * sin
    dev_add_globalspinor_assign_d(&(slocal[0]), &(sin[pos]));					// slocal  =  slocal + sin  =  pm_imu * (gamma5) * sin + sin
    dev_realmult_spinor_assigntoglobal_d(&(sout[pos]), 1.0, &(slocal[0]) );			// sout    =  slocal
  }
}






__global__ void dev_gamma5_d(dev_spinor_d * sin, dev_spinor_d * sout){
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



// the GPU implementation of  Qtm_pm_ndpsi(...)  from tm_operators_nd.c
// we have the possibility to add a constant shift > 0 
void dev_Qtm_pm_ndpsi_d (dev_spinor_d * spinout_up, dev_spinor_d * spinout_dn,
                              dev_spinor_d * spinin_up , dev_spinor_d * spinin_dn , 
                              int gridsize1, int blocksize1, int gridsize2, int blocksize2,
                              int gridsize3, int blocksize3, int gridsize4, int blocksize4) {
  
  
  // we will use the auxiliary fields  dev_spin_eo{1,2}_up/dn  for working on and buffering
  // and set  dev_spin_eo2_up/dn  equal  spinout_up/dn
  // spinin_up/dn  have to remain unchanged !!
  // spinout_up/dn  can be freely used
  

  int N_sites  =    VOLUME/2;		// #lattice sites


  dev_spin_eo2_up_d = spinout_up;		// need no memory allocated
  dev_spin_eo2_dn_d = spinout_dn;
  
 
  double nrm = 1.0 / (1.0 + g_mubar*g_mubar - g_epsbar*g_epsbar);
  
  dev_Hopping_Matrix_updn_d<<<gridsize1, blocksize1>>>(dev_gf_d, spinin_dn, spinin_up, dev_spin_eo1_up_d, dev_spin_eo1_dn_d,dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0,0,N_sites);	// dev_spin_eo1_up = (M_eo) * spinin_dn

  
  dev_mul_one_pm_imubar_gamma5_d<<<gridsize2, blocksize2>>>(dev_spin_eo1_up_d, dev_spin_eo2_up_d, -1.0);		// dev_spin_eo2_up  =  (1 - imubar) * dev_spin_eo1_up  =  (1 - imubar)*(M_eo) * spinin_dn
  dev_mul_one_pm_imubar_gamma5_d<<<gridsize2, blocksize2>>>(dev_spin_eo1_dn_d, dev_spin_eo2_dn_d, +1.0);		// dev_spin_eo2_dn  =  (1 + imubar) * dev_spin_eo1_dn  =  (1 + imubar)*(M_eo) * spinin_up
  
  dev_axpy_d<<<gridsize2, blocksize2>>>(g_epsbar, dev_spin_eo1_dn_d, dev_spin_eo2_up_d);
  
  dev_axpy_d<<<gridsize2, blocksize2>>>(g_epsbar, dev_spin_eo1_up_d, dev_spin_eo2_dn_d);

  										// dev_spin_eo2_up  =  dev_spin_eo2_up  +  epsbar * dev_spin_eo1_dn  =  (1 - imubar)*(M_eo) * spinin_dn  +  epsbar * (M_eo) * spinin_updev_xpay_d<<<gridsize2, blocksize2>>>(g_epsbar, dev_spin_eo1_up, dev_spin_eo2_dn);
  										// dev_spin_eo2_dn  =  dev_spin_eo2_dn  +  epsbar * dev_spin_eo1_up  =  (1 + imubar)*(M_eo) * spinin_up  +  epsbar * (M_eo) * spinin_dn
  dev_blasscal_d<<<gridsize2, blocksize2>>>(nrm, dev_spin_eo2_up_d);			// dev_spin_eo2_up  =  nrm * dev_spin_eo2_up  =  nrm*(1-imubar)*(M_eo) * spinin_dn  +  nrm*epsbar*(M_eo) * spinin_up
  dev_blasscal_d<<<gridsize2, blocksize2>>>(nrm, dev_spin_eo2_dn_d);			// dev_spin_eo2_dn  =  nrm * dev_spin_eo2_dn  =  nrm*(1+imubar)*(M_eo) * spinin_up  +  nrm*epsbar*(M_eo) * spinin_dn
  
  dev_Hopping_Matrix_updn_d<<<gridsize1, blocksize1>>>(dev_gf_d, dev_spin_eo2_up_d, dev_spin_eo2_dn_d, dev_spin_eo1_up_d, dev_spin_eo1_dn_d, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1,0,N_sites);
								

  dev_mul_one_pm_imubar_gamma5_d<<<gridsize2, blocksize2>>>(spinin_dn, dev_spin_eo2_up_d, +1.0);			// dev_spin_eo2_up = (1 + imubar) * spinin_dn
  dev_mul_one_pm_imubar_gamma5_d<<<gridsize2, blocksize2>>>(spinin_up, dev_spin_eo2_dn_d, -1.0);			// dev_spin_eo2_dn = (1 - imubar) * spinin_up
  

  dev_axpy_d<<<gridsize2, blocksize2>>>(-g_epsbar, spinin_up, dev_spin_eo2_up_d);
  												// dev_spin_eo2_up  =  dev_spin_eo2_up - epsbar*spinin_up  =  (1+imubar)*spinin_dn - epsbar*spinin_up
  dev_axpy_d<<<gridsize2, blocksize2>>>(-g_epsbar, spinin_dn, dev_spin_eo2_dn_d);
  												// dev_spin_eo2_dn  =  dev_spin_eo2_dn - epsbar*spinin_dn  =  (1-imubar)*spinin_up - epsbar*spinin_dn
  
  													// this is ((M_oo) - (M_oe)(Mee^-1)(M_eo)) * (spinin_dn, spinin_up):
  dev_axpy_d<<<gridsize2, blocksize2>>>(-1.0, dev_spin_eo1_up_d, dev_spin_eo2_up_d);
  							// dev_spin_eo2_up  =  dev_spin_eo2_up  -  dev_spin_eo1_up  =                      (1+imubar) * spinin_dn  -                    epsbar * spinin_up
  							//                                                             - (M_oe)*nrm*(1-imubar)*(M_eo) * spinin_dn  -  (M_oe)*nrm*epsbar*(M_eo) * spinin_up
  dev_axpy_d<<<gridsize2, blocksize2>>>(-1.0, dev_spin_eo1_dn_d, dev_spin_eo2_dn_d);
  							// dev_spin_eo2_dn  =  dev_spin_eo2_dn  -  dev_spin_eo1_dn  =                      (1-imubar) * spinin_up  -                    epsbar * spinin_dn
  							//                                                             - (M_oe)*nrm*(1+imubar)*(M_eo) * spinin_up  -  (M_oe)*nrm*epsbar*(M_eo) * spinin_dn
  

  dev_gamma5_d<<<gridsize3, blocksize3>>>(dev_spin_eo2_up_d, dev_spin_eo2_up_d);					// dev_spin_eo2_up = gamma5 * dev_spin_eo2_up 
  dev_gamma5_d<<<gridsize3, blocksize3>>>(dev_spin_eo2_dn_d, dev_spin_eo2_dn_d);					// dev_spin_eo2_dn = gamma5 * dev_spin_eo2_dn

 
  ////////////////////
  // (a,b) -> (b,a) //
  ////////////////////
  
  dev_blascopy_d<<<gridsize4, blocksize4>>>(dev_spin_eo2_dn_d, dev_spin_eo3_up_d);			// dev_spin_eo3_up = dev_spin_eo2_dn
  dev_blascopy_d<<<gridsize4, blocksize4>>>(dev_spin_eo2_up_d, dev_spin_eo3_dn_d);			// dev_spin_eo3_dn = dev_spin_eo2_up
      
  dev_Hopping_Matrix_updn_d<<<gridsize1, blocksize1>>>(dev_gf_d, dev_spin_eo3_up_d, dev_spin_eo3_dn_d, dev_spin_eo1_up_d, dev_spin_eo1_dn_d, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0,0,N_sites);	// dev_spin_eo1_up = (M_eo) * dev_spin_eo3_up

  dev_mul_one_pm_imubar_gamma5_d<<<gridsize2, blocksize2>>>(dev_spin_eo1_up_d, dev_spin_eo2_up_d, -1.0);			// dev_spin_eo2_up  =  (1 - imubar) * dev_spin_eo1_up  =  (1 - imubar)*(M_eo) * dev_spin_eo3_up
  dev_mul_one_pm_imubar_gamma5_d<<<gridsize2, blocksize2>>>(dev_spin_eo1_dn_d, dev_spin_eo2_dn_d, +1.0);			// dev_spin_eo2_dn  =  (1 + imubar) * dev_spin_eo1_dn  =  (1 + imubar)*(M_eo) * dev_spin_eo3_dn
  
  dev_axpy_d<<<gridsize2, blocksize2>>>(g_epsbar, dev_spin_eo1_dn_d, dev_spin_eo2_up_d);
  										// dev_spin_eo2_up  =  dev_spin_eo2_up  +  epsbar * dev_spin_eo1_dn  =  (1 - imubar)*(M_eo) * dev_spin_eo3_up  +  epsbar * (M_eo) * dev_spin_eo3_dn
  dev_axpy_d<<<gridsize2, blocksize2>>>(g_epsbar, dev_spin_eo1_up_d, dev_spin_eo2_dn_d);
  										// dev_spin_eo2_dn  =  dev_spin_eo2_dn  +  epsbar * dev_spin_eo1_up  =  (1 + imubar)*(M_eo) * dev_spin_eo3_dn  +  epsbar * (M_eo) * dev_spin_eo3_up
  
  dev_blasscal_d<<<gridsize2, blocksize2>>>(nrm, dev_spin_eo2_up_d);			// dev_spin_eo2_up  =  nrm * dev_spin_eo2_up  =  nrm*(1-imubar)*(M_eo) * dev_spin_eo3_up  +  nrm*epsbar*(M_eo) * dev_spin_eo3_dn
  dev_blasscal_d<<<gridsize2, blocksize2>>>(nrm, dev_spin_eo2_dn_d);			// dev_spin_eo2_dn  =  nrm * dev_spin_eo2_dn  =  nrm*(1+imubar)*(M_eo) * dev_spin_eo3_dn  +  nrm*epsbar*(M_eo) * dev_spin_eo3_up
  
  

  dev_Hopping_Matrix_updn_d<<<gridsize1, blocksize1>>>(dev_gf_d, dev_spin_eo2_up_d, dev_spin_eo2_dn_d, dev_spin_eo1_up_d,dev_spin_eo1_dn_d, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1,0,N_sites);

  

  dev_mul_one_pm_imubar_gamma5_d<<<gridsize2, blocksize2>>>(dev_spin_eo3_up_d, dev_spin_eo2_up_d, +1.0);				// dev_spin_eo2_up = (1 + imubar) * dev_spin_eo3_up
  dev_mul_one_pm_imubar_gamma5_d<<<gridsize2, blocksize2>>>(dev_spin_eo3_dn_d, dev_spin_eo2_dn_d, -1.0);				// dev_spin_eo2_dn = (1 - imubar) * dev_spin_eo3_dn
  
  											// remember: this is (M_oo) * (dev_spin_eo3_up, dev_spin_eo3_dn):
  dev_axpy_d<<<gridsize2, blocksize2>>>(-g_epsbar, dev_spin_eo3_dn_d, dev_spin_eo2_up_d);
  												// dev_spin_eo2_up  =  dev_spin_eo2_up - epsbar*dev_spin_eo3_dn  =  (1+imubar)*dev_spin_eo3_up - epsbar*dev_spin_eo3_dn
  dev_axpy_d<<<gridsize2, blocksize2>>>(-g_epsbar, dev_spin_eo3_up_d, dev_spin_eo2_dn_d);
  												// dev_spin_eo2_dn  =  dev_spin_eo2_dn - epsbar*dev_spin_eo3_up  =  (1-imubar)*dev_spin_eo3_dn - epsbar*dev_spin_eo3_up
  											// this is ( (M_oo) - (M_oe) (Mee^-1) (M_eo) ) * (dev_spin_eo3_up, dev_spin_eo3_dn)
  dev_axpy_d<<<gridsize2, blocksize2>>>(-1.0, dev_spin_eo1_up_d, dev_spin_eo2_up_d);
  							// dev_spin_eo2_up  =  dev_spin_eo2_up  -  dev_spin_eo1_up  =                      (1+imubar) * dev_spin_eo3_up  -                    epsbar * dev_spin_eo3_dn
  							//                                                             - (M_oe)*nrm*(1-imubar)*(M_eo) * dev_spin_eo3_up  -  (M_oe)*nrm*epsbar*(M_eo) * dev_spin_eo3_dn
  dev_axpy_d<<<gridsize2, blocksize2>>>( -1.0, dev_spin_eo1_dn_d, dev_spin_eo2_dn_d);
  							// dev_spin_eo2_dn  =  dev_spin_eo2_dn  -  dev_spin_eo1_dn  =                      (1-imubar) * dev_spin_eo3_dn  -                    epsbar * dev_spin_eo3_up
  							//                                                             - (M_oe)*nrm*(1+imubar)*(M_eo) * dev_spin_eo3_dn  -  (M_oe)*nrm*epsbar*(M_eo) * dev_spin_eo3_up
  
  dev_gamma5_d<<<gridsize3, blocksize3>>>(dev_spin_eo2_up_d, dev_spin_eo2_up_d);					// dev_spin_eo2_up = gamma5 * dev_spin_eo2_up 
  dev_gamma5_d<<<gridsize3, blocksize3>>>(dev_spin_eo2_dn_d, dev_spin_eo2_dn_d);					// dev_spin_eo2_dn = gamma5 * dev_spin_eo2_dn

  
  return;
  
}//dev_Qtm_pm_ndpsi_d()











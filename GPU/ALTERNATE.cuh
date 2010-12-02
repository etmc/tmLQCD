/**************************************************************************
 *
 * Copyright (C) 2010 Joseph Nagel
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
 **************************************************************************/




/////////////////////
// initializations //
/////////////////////


/*
//#ifdef MPI
  #ifdef ALTERNATE_HOPPING_MATRIX

// allocates memory for the fields for the alternative way of addressing positions in dev_Hopping_Matrix_alternate()

void init_gpu_indexfields() {
  
  size_t size;
  
  		// debug
  		//printf("Test: %p ?= %p\n", *g_iup, g_iup[0]);
  		//printf("Test: %p ?= %p\n", *g_idn, g_idn[0]);
  		//printf("Test: %p ?= %p\n", g_lexic2eo, &g_lexic2eo[0]);
  		//printf("Test: %p ?= %p\n", g_lexic2eosub, &g_lexic2eosub[0]);
  		//printf("Test: %p ?= %p\n", g_eo2lexic, &g_eo2lexic[0]);
  		//printf("Test: %p ?= %p ?= %p\n", ***g_ipt, &g_ipt[0][0][0][0], g_ipt[0][0][0]);
  
  size = 4*(VOLUME+RAND)*sizeof(int);
  cudaMalloc((void **) &dev_g_iup, size);
  cudaMalloc((void **) &dev_g_idn, size);
  cudaMemcpy(dev_g_iup, g_iup[0], size, cudaMemcpyHostToDevice);
  cudaMemcpy(dev_g_idn, g_idn[0], size, cudaMemcpyHostToDevice);
  
  size = (VOLUME+RAND)*sizeof(int);
  cudaMalloc((void **) &dev_g_lexic2eo, size);
  cudaMalloc((void **) &dev_g_lexic2eosub, size);
  cudaMemcpy(dev_g_lexic2eo, g_lexic2eo, size, cudaMemcpyHostToDevice);
  cudaMemcpy(dev_g_lexic2eosub, g_lexic2eosub, size, cudaMemcpyHostToDevice);
  
  size = (VOLUME+RAND)*sizeof(int);
  cudaMalloc((void **) &dev_g_eo2lexic, size);
  cudaMemcpy(dev_g_eo2lexic, g_eo2lexic, size, cudaMemcpyHostToDevice);
  
  size = VOLUME*sizeof(int);
  cudaMalloc((void **) &dev_g_ipt, size);
  cudaMemcpy(dev_g_ipt, g_ipt[0][0][0], size, cudaMemcpyHostToDevice);
  
}




// frees the memory

void free_gpu_indexfields() {

  cudaFree(dev_g_iup);
  cudaFree(dev_g_idn);
  
  cudaFree(dev_g_lexic2eo);
  cudaFree(dev_g_lexic2eosub);
  
  cudaFree(dev_g_eo2lexic);
  
  cudaFree(dev_g_ipt);

}

  #endif
//#endif	// MPI
*/






////////////////////
// hopping matrix //
////////////////////


/*
//#ifdef MPI
  #ifdef ALTERNATE_HOPPING_MATRIX

// applies the Hopping Part Even-Odd !
// the gauge field is the complete gaugefield!
// the gauge field at the local point is reconstructed by 2*pos+eo where pos is the eo-position
// from 0..VOLUME/2-1, eo = 0 or 1
// the positions in the gauge fields are passed in "gfindex_site" for gf's that are attached at
// the actual positions and in "gfindex_nextsite" for gf's that start at a position of the 
// other eo-sublattice.
// for the hopping positions of the eo-spinor field we use on of the two dedicated eo-nn fields
// the boundary conditions are implemented as in Hopping_Matrix.c
// mult with complex conjugate k0,k1,k2,k3 in positive direction because
// psi(x+mu) != exp(i theta_mu) psi(x)

__global__ void dev_Hopping_Matrix_alternate (const dev_su3_2v * gf, const dev_spinor * sin, dev_spinor * sout,
                                              int * dev_iup, int * dev_idn, int * dev_eo2lexic, int * dev_lexic2eosub,
                                              int ieo) {


  // guess: ieo = 0  corresponds to even sites ?!
  
  // USETEXTURE is not likely to work ... not now ...
  // same for TEMPORALGAUGE ...
  

  int pos_eo;
  int pos_global;
  int hoppos_eo;
  int hoppos_global;
  
  dev_spinor shelp1[6], ssum[6];
  __shared__ dev_su3_pad gfsmem[BLOCK];



  pos_eo = threadIdx.x + blockDim.x*blockIdx.x;  
  int ix = threadIdx.x;
  
  
  
  
  //////////
  // main //
  //////////
  
  
  if (pos_eo < dev_VOLUME) {
  
  
    if (ieo == 0)
      pos_global = dev_eo2lexic[pos_eo];
    else
      pos_global = dev_eo2lexic[dev_VOLUMEPLUSRAND/2 + pos_eo];
    
    
    dev_zero_spinor(&(ssum[0])); // zero sum  
    
        
    #ifdef TEMPORALGAUGE
      int spatialvol = dev_LX*dev_LY*dev_LZ;
    #endif
    
    
    
  
    ///////////////
    // l == 0, t //
    ///////////////
  
            // positive direction
            hoppos_global = dev_iup[4*pos_global + 0];
            hoppos_eo     = dev_lexic2eosub[hoppos_global];
            
            #ifdef TEMPORALGAUGE
              // gf == ID for t != T-1 => just read the spinor
              
              if((gfindex_site[pos]/spatialvol) != (dev_T-1) ){
              #ifdef USETEXTURE
                shelp1[0] = tex1Dfetch(spin_tex,6*hoppos);
                shelp1[1] = tex1Dfetch(spin_tex,6*hoppos+1);
                shelp1[2] = tex1Dfetch(spin_tex,6*hoppos+2);
                shelp1[3] = tex1Dfetch(spin_tex,6*hoppos+3);
                shelp1[4] = tex1Dfetch(spin_tex,6*hoppos+4);
                shelp1[5] = tex1Dfetch(spin_tex,6*hoppos+5);
              #else
                shelp1[0] = sin[6*hoppos];
                shelp1[1] = sin[6*hoppos+1];
                shelp1[2] = sin[6*hoppos+2];
                shelp1[3] = sin[6*hoppos+3];
                shelp1[4] = sin[6*hoppos+4];
                shelp1[5] = sin[6*hoppos+5];
              #endif
              }
              else{
                // gf != ID for t == T-1 => mult spinor with gf
                #ifdef GF_8
                  dev_reconstructgf_8texref(gf, 4*(gfindex_site[pos]),&(gfsmem[ix].m));
                #else
                  dev_reconstructgf_2vtexref(gf,4*(gfindex_site[pos]),&(gfsmem[ix].m));
                #endif
                #ifdef USETEXTURE
                  dev_su3MtV_spintex(gfsmem[ix].m, hoppos_eo, &(shelp1[0]));
                #else
                  dev_su3MtV(gfsmem[ix].m, &(sin[6*hoppos_eo]), &(shelp1[0]));
                #endif
              }
            #else
              #ifdef GF_8
                dev_reconstructgf_8texref(gf, 4*hoppos_global, &(gfsmem[ix].m));
              #else
                dev_reconstructgf_2vtexref(gf, 4*hoppos_global, &(gfsmem[ix].m));
              #endif
              #ifdef USETEXTURE
                dev_su3MtV_spintex(gfsmem[ix].m, hoppos_eo, &(shelp1[0]));
              #else
                dev_su3MtV(gfsmem[ix].m, &(sin[6*hoppos_eo]), &(shelp1[0]));
              #endif
            #endif
            
            //-kappa(r - gamma_mu)
            #ifdef GF_8
              dev_kappaP0_plus(&(ssum[0]), &(shelp1[0]), dev_cconj(dev_k0));
            #else
              dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_mk0,&(shelp1[0]), &(ssum[0]));
              dev_Gamma0(&(shelp1[0]));
              dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_k0,&(shelp1[0]), &(ssum[0]));
	    #endif
	    
	    
	    
	    
    ///////////////
    // l == 0, t //
    ///////////////

            // negative direction
            hoppos_global = dev_idn[4*pos_global + 0];
            hoppos_eo     = dev_lexic2eosub[hoppos_global];
            
             //hoppos = tex1Dfetch(nn_tex,8*pos+4);
            //color
            #ifdef TEMPORALGAUGE
              // gf == ID for t != T-1 => just read the spinor
              if((gfindex_nextsite[hoppos]/spatialvol) != (dev_T-1) ){
               #ifdef USETEXTURE
                shelp1[0] = tex1Dfetch(spin_tex,6*hoppos);
                shelp1[1] = tex1Dfetch(spin_tex,6*hoppos+1);
                shelp1[2] = tex1Dfetch(spin_tex,6*hoppos+2);
                shelp1[3] = tex1Dfetch(spin_tex,6*hoppos+3);
                shelp1[4] = tex1Dfetch(spin_tex,6*hoppos+4);
                shelp1[5] = tex1Dfetch(spin_tex,6*hoppos+5);
               #else
                shelp1[0] = sin[6*hoppos];
                shelp1[1] = sin[6*hoppos+1];
                shelp1[2] = sin[6*hoppos+2];
                shelp1[3] = sin[6*hoppos+3];
                shelp1[4] = sin[6*hoppos+4];
                shelp1[5] = sin[6*hoppos+5];
               #endif
              }
              else{
                // gf != ID for t == T-1 => mult spinor with gf
                #ifdef GF_8
                  dev_reconstructgf_8texref_dagger(gf,4*gfindex_nextsite[hoppos],&(gfsmem[ix].m));
                #else
                  dev_reconstructgf_2vtexref_dagger(gf,4*gfindex_nextsite[hoppos],&(gfsmem[ix].m));
                #endif
                #ifdef USETEXTURE
                  dev_su3MtV_spintex(gfsmem[ix].m, hoppos, &(shelp1[0]));
                #else
                  dev_su3MtV(gfsmem[ix].m, &(sin[6*hoppos]), &(shelp1[0]));
                #endif 
              }
            #else            
              #ifdef GF_8
                dev_reconstructgf_8texref_dagger(gf, 4*hoppos_global, &(gfsmem[ix].m));
              #else
                dev_reconstructgf_2vtexref_dagger(gf, 4*hoppos_global, &(gfsmem[ix].m));
              #endif
              #ifdef USETEXTURE
                dev_su3MtV_spintex(gfsmem[ix].m, hoppos_eo, &(shelp1[0]));  
              #else
                dev_su3MtV(gfsmem[ix].m, &(sin[6*hoppos_eo]), &(shelp1[0]));
              #endif 
            #endif
            
            //-kappa(r + gamma_mu)
            #ifdef GF_8
              dev_kappaP0_minus(&(ssum[0]), &(shelp1[0]), dev_k0);
            #else
              dev_complexmult_add_assign_spinor(&(ssum[0]),dev_mk0,&(shelp1[0]), &(ssum[0]));
              dev_Gamma0(&(shelp1[0]));
              dev_complexmult_add_assign_spinor(&(ssum[0]),dev_mk0,&(shelp1[0]), &(ssum[0]));
            #endif




    ///////////////
    // l == 3, z //
    ///////////////

            // positive direction
            hoppos_global = dev_iup[4*pos_global + 3];
            hoppos_eo     = dev_lexic2eosub[hoppos_global];
            
             //hoppos = tex1Dfetch(nn_tex,8*pos+3);
            //color
            #ifdef GF_8
              dev_reconstructgf_8texref(gf, 4*(hoppos_global)+(3), &(gfsmem[ix].m));
            #else
              dev_reconstructgf_2vtexref(gf, 4*(hoppos_global)+(3),&(gfsmem[ix].m));
            #endif
            #ifdef USETEXTURE
              dev_su3MtV_spintex(gfsmem[ix].m, hoppos_eo, &(shelp1[0]));
            #else
              dev_su3MtV(gfsmem[ix].m, &(sin[6*hoppos_eo]), &(shelp1[0]));
            #endif
            //-kappa(r - gamma_mu)    
            #ifdef GF_8
              dev_kappaP3_plus(&(ssum[0]), &(shelp1[0]), dev_k3.re);
            #else
              dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_mk3,&(shelp1[0]), &(ssum[0]));
              dev_Gamma3(&(shelp1[0]));
              dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_k3,&(shelp1[0]), &(ssum[0]));
	    #endif




    ///////////////
    // l == 3, z //
    ///////////////
            
            // negative direction
            hoppos_global = dev_idn[4*pos_global + 3];
            hoppos_eo     = dev_lexic2eosub[hoppos_global];
            
             //hoppos = tex1Dfetch(nn_tex,8*pos+7); 
            //color
            #ifdef GF_8
              dev_reconstructgf_8texref_dagger(gf, 4*hoppos_global+(3), &(gfsmem[ix].m));
            #else
              dev_reconstructgf_2vtexref_dagger(gf, 4*hoppos_global+(3), &(gfsmem[ix].m));
            #endif
            #ifdef USETEXTURE
              dev_su3MtV_spintex(gfsmem[ix].m, hoppos_eo, &(shelp1[0]));
            #else
              dev_su3MtV(gfsmem[ix].m, &(sin[6*hoppos_eo]), &(shelp1[0]));
            #endif
            //-kappa(r + gamma_mu)
            #ifdef GF_8
              dev_kappaP3_minus(&(ssum[0]), &(shelp1[0]), dev_k3.re);
            #else
              dev_complexmult_add_assign_spinor(&(ssum[0]),dev_mk3,&(shelp1[0]), &(ssum[0]));
              dev_Gamma3(&(shelp1[0]));
              dev_complexmult_add_assign_spinor(&(ssum[0]),dev_mk3,&(shelp1[0]), &(ssum[0]));
            #endif




    ///////////////
    // l == 2, y //
    ///////////////

            // positive direction
            hoppos_global = dev_iup[4*pos_global + 2];
            hoppos_eo     = dev_lexic2eosub[hoppos_global];
            
             //hoppos = tex1Dfetch(nn_tex,8*pos+2);
            //color
            #ifdef GF_8
              dev_reconstructgf_8texref(gf, 4*(hoppos_global)+(2), &(gfsmem[ix].m));
            #else
              dev_reconstructgf_2vtexref(gf, 4*(hoppos_global)+(2), &(gfsmem[ix].m));
            #endif
            #ifdef USETEXTURE
              dev_su3MtV_spintex(gfsmem[ix].m, hoppos_eo, &(shelp1[0]));
            #else
              dev_su3MtV(gfsmem[ix].m, &(sin[6*hoppos_eo]), &(shelp1[0]));
            #endif
            //-kappa(r - gamma_mu)
            #ifdef GF_8
              dev_kappaP2_plus(&(ssum[0]), &(shelp1[0]), dev_k2.re);
            #else
              dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_mk2,&(shelp1[0]), &(ssum[0]));
              dev_Gamma2(&(shelp1[0]));
              dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_k2,&(shelp1[0]), &(ssum[0]));
            #endif




    ///////////////
    // l == 2, y //
    ///////////////
            
            // negative direction
            hoppos_global = dev_idn[4*pos_global + 2];
            hoppos_eo     = dev_lexic2eosub[hoppos_global];
            
             //hoppos = tex1Dfetch(nn_tex,8*pos+6);
            //color
            #ifdef GF_8
              dev_reconstructgf_8texref_dagger(gf, 4*(hoppos_global)+(2), &(gfsmem[ix].m));
            #else
              dev_reconstructgf_2vtexref_dagger(gf, 4*(hoppos_global)+(2), &(gfsmem[ix].m));
            #endif
            #ifdef USETEXTURE
              dev_su3MtV_spintex(gfsmem[ix].m, hoppos_eo, &(shelp1[0]));
            #else
              dev_su3MtV(gfsmem[ix].m, &(sin[6*hoppos_eo]), &(shelp1[0]));
            #endif
            //-kappa(r + gamma_mu)
            #ifdef GF_8
              dev_kappaP2_minus(&(ssum[0]), &(shelp1[0]), dev_k2.re);
            #else
              dev_complexmult_add_assign_spinor(&(ssum[0]),dev_mk2,&(shelp1[0]), &(ssum[0]));
              dev_Gamma2(&(shelp1[0]));
              dev_complexmult_add_assign_spinor(&(ssum[0]),dev_mk2,&(shelp1[0]), &(ssum[0]));
	    #endif




    ///////////////
    // l == 1, x //
    ///////////////

            // positive direction
            hoppos_global = dev_iup[4*pos_global + 1];
            hoppos_eo     = dev_lexic2eosub[hoppos_global];
            
             //hoppos = tex1Dfetch(nn_tex,8*pos+1);
            //color
            #ifdef GF_8
              dev_reconstructgf_8texref(gf, 4*(hoppos_global)+(1), &(gfsmem[ix].m));
            #else
              dev_reconstructgf_2vtexref(gf, 4*(hoppos_global)+(1), &(gfsmem[ix].m));
            #endif
            #ifdef USETEXTURE
              dev_su3MtV_spintex(gfsmem[ix].m, hoppos_eo, &(shelp1[0]));
            #else
              dev_su3MtV(gfsmem[ix].m, &(sin[6*hoppos_eo]), &(shelp1[0]));
            #endif
            //-kappa(r - gamma_mu)
            #ifdef GF_8
              dev_kappaP1_plus(&(ssum[0]), &(shelp1[0]), dev_k1.re);
            #else
              dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_mk1,&(shelp1[0]), &(ssum[0]));
              dev_Gamma1(&(shelp1[0]));
              dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_k1,&(shelp1[0]), &(ssum[0]));
	    #endif




    ///////////////
    // l == 1, x //
    ///////////////
            
            // negative direction
            hoppos_global = dev_idn[4*pos_global + 1];
            hoppos_eo     = dev_lexic2eosub[hoppos_global];
            
             //hoppos = tex1Dfetch(nn_tex,8*pos+5);
            //color
            #ifdef GF_8
              dev_reconstructgf_8texref_dagger(gf, 4*(hoppos_global)+(1), &(gfsmem[ix].m));
            #else
              dev_reconstructgf_2vtexref_dagger(gf, 4*(hoppos_global)+(1), &(gfsmem[ix].m));
            #endif
            #ifdef USETEXTURE
              dev_su3MtV_spintex(gfsmem[ix].m, hoppos_eo, &(shelp1[0]));
            #else
              dev_su3MtV(gfsmem[ix].m, &(sin[6*hoppos_eo]), &(shelp1[0]));
            #endif
            //-kappa(r + gamma_mu)
            #ifdef GF_8
              dev_kappaP1_minus(&(ssum[0]), &(shelp1[0]), dev_k1.re);
            #else
              dev_complexmult_add_assign_spinor(&(ssum[0]),dev_mk1,&(shelp1[0]), &(ssum[0]));
              dev_Gamma1(&(shelp1[0]));
              dev_complexmult_add_assign_spinor(&(ssum[0]),dev_mk1,&(shelp1[0]), &(ssum[0]));      
            #endif
 
 
 
 
    /////////////
    // output //
    ////////////
  
        //copy to output spinor
        dev_copy_spinor(&(ssum[0]),&(sout[6*pos_eo])); 
        
  }
  
  
}//dev_Hopping_Matrix_alternate<<<>>>()

  #endif
//#endif	// MPI
*/






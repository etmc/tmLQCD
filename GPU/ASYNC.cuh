

#define EXTERN extern
				// taken from global.h
EXTERN MPI_Status status;
EXTERN MPI_Request req1,req2,req3,req4;
EXTERN MPI_Comm g_cart_grid;
EXTERN MPI_Comm g_mpi_time_slices;
EXTERN MPI_Comm g_mpi_SV_slices;
EXTERN MPI_Comm g_mpi_z_slices;
EXTERN MPI_Comm g_mpi_ST_slices;

/* the next neighbours for MPI */
EXTERN int g_nb_x_up, g_nb_x_dn;
EXTERN int g_nb_y_up, g_nb_y_dn;
EXTERN int g_nb_t_up, g_nb_t_dn;
EXTERN int g_nb_z_up, g_nb_z_dn;




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

__global__ void dev_Hopping_Matrix_ASYNC (const dev_su3_2v * gf, 
                                          const dev_spinor * sin, dev_spinor * sout,
                                          const int * gfindex_site, const int* gfindex_nextsite, const int * nn_evenodd,
                                          const int eo,
                                          int start, int size) {
  
  int pos, hoppos;
  
  
  dev_spinor shelp1[6], ssum[6];
  __shared__ dev_su3_pad gfsmem[BLOCK];
  
  
  
  pos = start  +  threadIdx.x + blockDim.x * blockIdx.x;
  int ix = threadIdx.x;
  
  
  if (pos < start + size) {
  
  
  dev_zero_spinor(&(ssum[0])); // zero sum
  
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
                  dev_su3MtV_spintex(gfsmem[ix].m, hoppos, &(shelp1[0]));
                #else
                  dev_su3MtV(gfsmem[ix].m, &(sin[6*hoppos]), &(shelp1[0]));
                #endif
              }
            #else
              #ifdef GF_8
                dev_reconstructgf_8texref(gf, 4*(gfindex_site[pos]),&(gfsmem[ix].m));
              #else
                dev_reconstructgf_2vtexref(gf, 4*(gfindex_site[pos]),&(gfsmem[ix].m));
              #endif
              #ifdef USETEXTURE
                dev_su3MtV_spintex(gfsmem[ix].m, hoppos, &(shelp1[0]));
              #else
                dev_su3MtV(gfsmem[ix].m, &(sin[6*hoppos]), &(shelp1[0]));
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
                dev_reconstructgf_8texref_dagger(gf,4*gfindex_nextsite[hoppos],&(gfsmem[ix].m));
              #else
                dev_reconstructgf_2vtexref_dagger(gf,4*gfindex_nextsite[hoppos],&(gfsmem[ix].m));
              #endif
              #ifdef USETEXTURE
                dev_su3MtV_spintex(gfsmem[ix].m, hoppos, &(shelp1[0]));  
              #else
                dev_su3MtV(gfsmem[ix].m, &(sin[6*hoppos]), &(shelp1[0]));
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




//l==3,z 
            //positive direction
            hoppos = nn_evenodd[8*pos+3];
             //hoppos = tex1Dfetch(nn_tex,8*pos+3);
            //color
            #ifdef GF_8
              dev_reconstructgf_8texref(gf,4*(gfindex_site[pos])+(3),&(gfsmem[ix].m));
            #else
              dev_reconstructgf_2vtexref(gf, 4*(gfindex_site[pos])+(3),&(gfsmem[ix].m));
            #endif
            #ifdef USETEXTURE
              dev_su3MtV_spintex(gfsmem[ix].m, hoppos, &(shelp1[0]));
            #else
              dev_su3MtV(gfsmem[ix].m, &(sin[6*hoppos]), &(shelp1[0]));
            #endif
            //-kappa(r - gamma_mu)    
            #ifdef GF_8
              dev_kappaP3_plus(&(ssum[0]), &(shelp1[0]), dev_k3.re);
            #else
              dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_mk3,&(shelp1[0]), &(ssum[0]));
              dev_Gamma3(&(shelp1[0]));
              dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_k3,&(shelp1[0]), &(ssum[0]));
	    #endif
//l==3,z               
            
            //negative direction
            hoppos = nn_evenodd[8*pos+7];
             //hoppos = tex1Dfetch(nn_tex,8*pos+7); 
            //color
            #ifdef GF_8
              dev_reconstructgf_8texref_dagger(gf,4*gfindex_nextsite[hoppos]+(3),&(gfsmem[ix].m));
            #else
              dev_reconstructgf_2vtexref_dagger(gf,4*gfindex_nextsite[hoppos]+(3),&(gfsmem[ix].m));
            #endif
            #ifdef USETEXTURE
              dev_su3MtV_spintex(gfsmem[ix].m, hoppos, &(shelp1[0]));
            #else
              dev_su3MtV(gfsmem[ix].m, &(sin[6*hoppos]), &(shelp1[0]));
            #endif
            //-kappa(r + gamma_mu)
            #ifdef GF_8
              dev_kappaP3_minus(&(ssum[0]), &(shelp1[0]), dev_k3.re);
            #else
              dev_complexmult_add_assign_spinor(&(ssum[0]),dev_mk3,&(shelp1[0]), &(ssum[0]));
              dev_Gamma3(&(shelp1[0]));
              dev_complexmult_add_assign_spinor(&(ssum[0]),dev_mk3,&(shelp1[0]), &(ssum[0]));
            #endif




//l==2,y 
            //positive direction
            hoppos = nn_evenodd[8*pos+2];
             //hoppos = tex1Dfetch(nn_tex,8*pos+2);
            //color
            #ifdef GF_8
              dev_reconstructgf_8texref(gf,4*(gfindex_site[pos])+(2),&(gfsmem[ix].m));
            #else
              dev_reconstructgf_2vtexref(gf,4*(gfindex_site[pos])+(2),&(gfsmem[ix].m));
            #endif
            #ifdef USETEXTURE
              dev_su3MtV_spintex(gfsmem[ix].m, hoppos, &(shelp1[0]));
            #else
              dev_su3MtV(gfsmem[ix].m, &(sin[6*hoppos]), &(shelp1[0]));
            #endif
            //-kappa(r - gamma_mu)
            #ifdef GF_8
              dev_kappaP2_plus(&(ssum[0]), &(shelp1[0]), dev_k2.re);
            #else
              dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_mk2,&(shelp1[0]), &(ssum[0]));
              dev_Gamma2(&(shelp1[0]));
              dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_k2,&(shelp1[0]), &(ssum[0]));
            #endif

//l==2,y        

            
            //negative direction
            hoppos = nn_evenodd[8*pos+6]; 
             //hoppos = tex1Dfetch(nn_tex,8*pos+6);
            //color
            #ifdef GF_8
              dev_reconstructgf_8texref_dagger(gf,4*gfindex_nextsite[hoppos]+(2),&(gfsmem[ix].m));
            #else
              dev_reconstructgf_2vtexref_dagger(gf,4*gfindex_nextsite[hoppos]+(2),&(gfsmem[ix].m));
            #endif
            #ifdef USETEXTURE
              dev_su3MtV_spintex(gfsmem[ix].m, hoppos, &(shelp1[0]));
            #else
              dev_su3MtV(gfsmem[ix].m, &(sin[6*hoppos]), &(shelp1[0]));
            #endif
            //-kappa(r + gamma_mu)
            #ifdef GF_8
              dev_kappaP2_minus(&(ssum[0]), &(shelp1[0]), dev_k2.re);
            #else
              dev_complexmult_add_assign_spinor(&(ssum[0]),dev_mk2,&(shelp1[0]), &(ssum[0]));
              dev_Gamma2(&(shelp1[0]));
              dev_complexmult_add_assign_spinor(&(ssum[0]),dev_mk2,&(shelp1[0]), &(ssum[0]));
	    #endif



//l==1,x 
            //positive direction
            hoppos = nn_evenodd[8*pos+1];
             //hoppos = tex1Dfetch(nn_tex,8*pos+1);
            //color
            #ifdef GF_8
              dev_reconstructgf_8texref(gf,4*(gfindex_site[pos])+(1),&(gfsmem[ix].m));
            #else
              dev_reconstructgf_2vtexref(gf,4*(gfindex_site[pos])+(1),&(gfsmem[ix].m));
            #endif
            #ifdef USETEXTURE
              dev_su3MtV_spintex(gfsmem[ix].m, hoppos, &(shelp1[0]));
            #else
              dev_su3MtV(gfsmem[ix].m, &(sin[6*hoppos]), &(shelp1[0]));
            #endif
            //-kappa(r - gamma_mu)
            #ifdef GF_8
              dev_kappaP1_plus(&(ssum[0]), &(shelp1[0]), dev_k1.re);
            #else
              dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_mk1,&(shelp1[0]), &(ssum[0]));
              dev_Gamma1(&(shelp1[0]));
              dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_k1,&(shelp1[0]), &(ssum[0]));
	    #endif


//l==1,x 
            
            //negative direction
            hoppos = nn_evenodd[8*pos+5]; 
             //hoppos = tex1Dfetch(nn_tex,8*pos+5);
            //color
            #ifdef GF_8
              dev_reconstructgf_8texref_dagger(gf,4*gfindex_nextsite[hoppos]+(1),&(gfsmem[ix].m));
            #else
              dev_reconstructgf_2vtexref_dagger(gf,4*gfindex_nextsite[hoppos]+(1),&(gfsmem[ix].m));
            #endif
            #ifdef USETEXTURE
              dev_su3MtV_spintex(gfsmem[ix].m, hoppos, &(shelp1[0]));
            #else
              dev_su3MtV(gfsmem[ix].m, &(sin[6*hoppos]), &(shelp1[0]));
            #endif
            //-kappa(r + gamma_mu)
            #ifdef GF_8
              dev_kappaP1_minus(&(ssum[0]), &(shelp1[0]), dev_k1.re);
            #else
              dev_complexmult_add_assign_spinor(&(ssum[0]),dev_mk1,&(shelp1[0]), &(ssum[0]));
              dev_Gamma1(&(shelp1[0]));
              dev_complexmult_add_assign_spinor(&(ssum[0]),dev_mk1,&(shelp1[0]), &(ssum[0]));      
            #endif
 
        //copy to output spinor
        dev_copy_spinor(&(ssum[0]),&(sout[6*pos])); 
  }
}//dev_Hopping_Matrix_ASYNC()










void HOPPING_ASYNC (dev_su3_2v * gf, 
                    dev_spinor * spinin, dev_spinor * spinout,
                    int * gfindex_site, int * gfindex_nextsite, int * nn_evenodd,
                    int ieo,
                    int gridsize, int blocksize) {
  
  
  // for even/odd
  int tSliceEO = LX*LY*LZ/2;
  int VolumeEO = VOLUME/2;
  
  #ifdef ASYNC_OPTIMIZED
    int offset;
    if (tSliceEO % nStreams == 0) {
      offset = tSliceEO / nStreams;
    }
    else {
      printf("Error in HOPPING_ASYNC(): tSliceEO is not divisible by nStreams!\n");
      exit(-1);
    } 
  #endif
  
  // gridsizes
  int gridsize1;
  int gridsize2;
  
  if ( (VolumeEO-2*tSliceEO) % blocksize == 0 ) {
    gridsize1  = (VolumeEO-2*tSliceEO) / blocksize;
  }
  else {
    gridsize1  = (int) ( ((VolumeEO-2*tSliceEO)/blocksize) + 1);
  }
  
  if ( (tSliceEO) % blocksize == 0 ) {
    gridsize2  = (tSliceEO) / blocksize;
  }
  else {
    gridsize2  = (int) ( ((tSliceEO)/blocksize) + 1);
  }
  
  
  
  
  #ifdef USETEXTURE
    bind_texture_spin(spinin,1);
  #endif
  
  
  
  
  // exchanges the boundaries
  #ifndef ASYNC_OPTIMIZED	// primitive version

	// applies to the parts which don't need communication
  	dev_Hopping_Matrix_ASYNC <<<gridsize1, blocksize>>> ( gf,
                                                              spinin, spinout,
                                                              gfindex_site, gfindex_nextsite, nn_evenodd,
                                                              ieo,
                                                              //2*tSliceEO, VolumeEO-4*tSliceEO );
                                                              tSliceEO, VolumeEO-2*tSliceEO );
  	
  	xchange_field_wrapper(spinin, ieo);			// to be further optimized !!
  	
  	// applies the hopping matrix to remaining parts
  	dev_Hopping_Matrix_ASYNC <<<gridsize2, blocksize>>> ( gf,
  	                                                      spinin, spinout,
  	                                                      gfindex_site, gfindex_nextsite, nn_evenodd,
  	                                                      ieo,
  	                                                      //0, 2*tSliceEO );
  	                                                      0, tSliceEO );
  	
  	dev_Hopping_Matrix_ASYNC <<<gridsize2, blocksize>>> ( gf,
  	                                                      spinin, spinout,
  	                                                      gfindex_site, gfindex_nextsite, nn_evenodd,
  	                                                      ieo,
  	                                                      //VolumeEO-2*tSliceEO, 2*tSliceEO );
  	                                                      VolumeEO-tSliceEO, tSliceEO );	



  
  #else				// optimized version


  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		
		/*
				#ifdef ASYNC_TIMING
  				  cudaEventRecord(start_ALL, 0);
  				#endif
        	
        	
  		// copies first FACE to host
  		cudaMemcpyAsync(RAND1, spinin                      , tSliceEO*6*sizeof(float4), cudaMemcpyDeviceToHost, stream[1]);
  		
  				#ifdef ASYNC_TIMING
  				  cudaEventRecord(stop_D2H_1, stream[1]);
  				#endif
  		
  		
  		// INTERNAL kernel
  		dev_Hopping_Matrix_ASYNC <<<gridsize1, blocksize, 0, stream[0]>>> ( gf,
        	                                                                    spinin, spinout,
        	                                                                    gfindex_site, gfindex_nextsite, nn_evenodd,
        	                                                                    ieo,
        	                                                                    tSliceEO, VolumeEO-2*tSliceEO );
        			#ifdef ASYNC_TIMING
        			  cudaEventRecord(stop_INT_0, stream[0]);
        			#endif
  		
  		
  		// exchanges first FACE
  		cudaStreamSynchronize(stream[1]);												// synchronous
  			//MPI_Irecv(RAND3, 24*tSliceEO, MPI_FLOAT, g_nb_t_up, 0,
  			//          g_cart_grid, &recv_req[0]);
  			//MPI_Isend(RAND1, 24*tSliceEO, MPI_FLOAT, g_nb_t_dn, 0,
  			//          g_cart_grid, &send_req[0]);
  		
  			MPI_Sendrecv(RAND1, 24*tSliceEO, MPI_FLOAT, g_nb_t_dn, 0,
  			             RAND3, 24*tSliceEO, MPI_FLOAT, g_nb_t_up, 0,
  			             g_cart_grid, &stat[0]);
  		
  		
  		// copies first FACE back to device												// order may switched
  		//MPI_Wait(&recv_req[0], &stat[0]);												// synchronous
  		cudaMemcpyAsync(spinin+6*VolumeEO, RAND3, tSliceEO*6*sizeof(float4), cudaMemcpyHostToDevice, stream[1]);
  		
  				#ifdef ASYNC_TIMING
  				  cudaEventRecord(stop_H2D_3, stream[1]);
  				#endif
  		
  			
  		
  		
  		// applies first FACE
  		dev_Hopping_Matrix_ASYNC <<<gridsize2, blocksize, 0, stream[1]>>> ( gf,
  		                                                                      spinin, spinout,
  		                                                                      gfindex_site, gfindex_nextsite, nn_evenodd,
  		                                                                      ieo,
  		                                                                      VolumeEO-tSliceEO, tSliceEO );
  				#ifdef ASYNC_TIMING
  				  cudaEventRecord(stop_EXT_1, stream[1]);
  				#endif
  		
  		
  		// copies second FACE to host
  		cudaMemcpyAsync(RAND2, spinin+6*(VolumeEO-tSliceEO), tSliceEO*6*sizeof(float4), cudaMemcpyDeviceToHost, stream[2]);
  		
  				#ifdef ASYNC_TIMING
  				  cudaEventRecord(stop_D2H_2, stream[2]);
  				#endif
  		
  		
  		// exchanges second FACE
  		cudaStreamSynchronize(stream[2]);												// synchronous
  			//MPI_Irecv(RAND4, 24*tSliceEO, MPI_FLOAT, g_nb_t_dn, 1,
  			//          g_cart_grid, &recv_req[1]);
  			//MPI_Isend(RAND2, 24*tSliceEO, MPI_FLOAT, g_nb_t_up, 1,
  			//          g_cart_grid, &send_req[1]);
  		
  			MPI_Sendrecv(RAND2, 24*tSliceEO, MPI_FLOAT, g_nb_t_up, 1,
  			             RAND4, 24*tSliceEO, MPI_FLOAT, g_nb_t_dn, 1,
  			             g_cart_grid, &stat[1]);
  		
  		
  		
  		
  		
  		// copies second FACE back to device
  		//MPI_Wait(&recv_req[1], &stat[1]);												// synchronous
  		cudaMemcpyAsync(spinin+6*(VolumeEO+tSliceEO), RAND4, tSliceEO*6*sizeof(float4), cudaMemcpyHostToDevice, stream[2]);
  		
  				#ifdef ASYNC_TIMING
  				  cudaEventRecord(stop_H2D_4, stream[2]);
  				#endif
  		
  		
  		// applies second FACE
  		dev_Hopping_Matrix_ASYNC <<<gridsize2, blocksize, 0, stream[2]>>> ( gf,
  		                                                                    spinin, spinout,
  		                                                                    gfindex_site, gfindex_nextsite, nn_evenodd,
  		                                                                    ieo,
  		                                                                    0, tSliceEO );
  				#ifdef ASYNC_TIMING
  				  cudaEventRecord(stop_EXT_2, stream[2]);
  				#endif
  		
  		
  				#ifdef ASYNC_TIMING
  				  cudaEventRecord(stop_ALL, 0);
  				#endif
  		*/
  		
  		
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  		
  		
  		
  				#ifdef ASYNC_TIMING
  				  cudaEventRecord(start_ALL, 0);
  				#endif
  		
  		
  		// copies first FACE to host
  		cudaMemcpyAsync(RAND1, spinin                      , tSliceEO*6*sizeof(float4), cudaMemcpyDeviceToHost, stream[1]);
  		
  				#ifdef ASYNC_TIMING
  				  cudaEventRecord(stop_D2H_1, stream[1]);
  				#endif
  		
  		// copies second FACE to host
  		cudaMemcpyAsync(RAND2, spinin+6*(VolumeEO-tSliceEO), tSliceEO*6*sizeof(float4), cudaMemcpyDeviceToHost, stream[2]);
  		
  				#ifdef ASYNC_TIMING
        			  cudaEventRecord(stop_D2H_2, stream[2]);
        			#endif
  		
  		
  		// INTERNAL kernel
  		dev_Hopping_Matrix_ASYNC <<<gridsize1, blocksize, 0, stream[0]>>> ( gf,
        	                                                                    spinin, spinout,
        	                                                                    gfindex_site, gfindex_nextsite, nn_evenodd,
        	                                                                    ieo,
        	                                                                    tSliceEO, VolumeEO-2*tSliceEO );
        			#ifdef ASYNC_TIMING
        			  cudaEventRecord(stop_INT_0, stream[0]);
        			#endif
        	
  		
  		// first FACE
  		cudaStreamSynchronize(stream[1]);
  		
  		for (int i = 0; i < nStreams; i++) {
  		
  		  MPI_Sendrecv(RAND1+6*i*offset, 24*offset, MPI_FLOAT, g_nb_t_dn, 0,		// NOT asynchronous
  		               RAND3+6*i*offset, 24*offset, MPI_FLOAT, g_nb_t_up, 0,
  		               g_cart_grid, &stat[i]);
  		  
  		  //MPI_Isend(RAND1+6*i*offset, 24*offset, MPI_FLOAT, g_nb_t_dn, i,
  		  //          g_cart_grid, &send_req[i]);
  		  //MPI_Irecv (RAND3+6*i*offset, 24*offset, MPI_FLOAT, g_nb_t_up, i,
  		  //          g_cart_grid, &recv_req[i]);
  		  
  		  //MPI_Wait(&recv_req[i], &stat[i]);
  		           
  		  cudaMemcpyAsync(spinin+6*VolumeEO+6*i*offset, RAND3+6*i*offset, offset*6*sizeof(float4), cudaMemcpyHostToDevice, stream[1+i]);
  		  
  		  		#ifdef ASYNC_TIMING
  		  		  cudaEventRecord(stop_H2D_3, stream[1]);
  		  		#endif
  		  
  		  dev_Hopping_Matrix_ASYNC <<<gridsize2, blocksize, 0, stream[1+i]>>> ( gf,
  		                                                                        spinin, spinout,
  		                                                                        gfindex_site, gfindex_nextsite, nn_evenodd,
  		                                                                        ieo,
  		                                                                        VolumeEO-tSliceEO+i*offset, offset );
  		  		#ifdef ASYNC_TIMING
  		  		  cudaEventRecord(stop_EXT_1, stream[1]);
  		  		#endif
  		  
  		}
  		
  		
  		// second FACE
  		cudaStreamSynchronize(stream[nStreams+1]);
  		
  		for (int i = 0; i < nStreams; i++) {
  		
  		  MPI_Sendrecv(RAND2+6*i*offset, 24*offset, MPI_FLOAT, g_nb_t_up, 1,
  		               RAND4+6*i*offset, 24*offset, MPI_FLOAT, g_nb_t_dn, 1,
  		               g_cart_grid, &stat[nStreams+i]);
		  
  		  //MPI_Isend(RAND2+6*i*offset, 24*offset, MPI_FLOAT, g_nb_t_up, nStreams+i,
  		  //          g_cart_grid, &send_req[nStreams+i]);
  		  //MPI_Irecv (RAND4+6*i*offset, 24*offset, MPI_FLOAT, g_nb_t_dn, nStreams+i,
  		  //          g_cart_grid, &recv_req[nStreams+i]);
  		  
  		  //MPI_Wait(&recv_req[nStreams+i], &stat[nStreams+i]);
  		  
  		  cudaMemcpyAsync(spinin+6*(VolumeEO+tSliceEO)+6*i*offset, RAND4+6*i*offset, offset*6*sizeof(float4), cudaMemcpyHostToDevice, stream[nStreams+1+i]);
  		  
  		  		#ifdef ASYNC_TIMING
  		  		  cudaEventRecord(stop_H2D_4, stream[2]);
  		  		#endif
  		  
  		  dev_Hopping_Matrix_ASYNC <<<gridsize2, blocksize, 0, stream[nStreams+1+i]>>> ( gf,
  		                                                                                 spinin, spinout,
  		                                                                                 gfindex_site, gfindex_nextsite, nn_evenodd,
  		                                                                                 ieo,
  		                                                                                 0+i*offset, offset );
  		  		#ifdef ASYNC_TIMING
  		  		  cudaEventRecord(stop_EXT_2, stream[2]);
  		  		#endif
  		  
		}
  		
  		
  				#ifdef ASYNC_TIMING
  				  cudaEventRecord(stop_ALL, 0);
  				#endif
  		
  		
  		
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  		
  		
  		/*
  				#ifdef ASYNC_TIMING
  				  cudaEventRecord(start_ALL, 0);
  				#endif
  		
  		
  		// copies first FACE to host
  		cudaMemcpyAsync(RAND1, spinin                      , tSliceEO*6*sizeof(float4), cudaMemcpyDeviceToHost, stream[1]);
  		
  				#ifdef ASYNC_TIMING
  				  cudaEventRecord(stop_D2H_1, stream[1]);
  				#endif
  		
  		
  		// copies second FACE to host
  		cudaMemcpyAsync(RAND2, spinin+6*(VolumeEO-tSliceEO), tSliceEO*6*sizeof(float4), cudaMemcpyDeviceToHost, stream[2]);
  		
  				#ifdef ASYNC_TIMING
        			  cudaEventRecord(stop_D2H_2, stream[2]);
        			#endif
  		
  		
  		// INTERNAL kernel
  		dev_Hopping_Matrix_ASYNC <<<gridsize1, blocksize, 0, stream[0]>>> ( gf,
        	                                                                    spinin, spinout,
        	                                                                    gfindex_site, gfindex_nextsite, nn_evenodd,
        	                                                                    ieo,
        	                                                                    tSliceEO, VolumeEO-2*tSliceEO );
  				#ifdef ASYNC_TIMING
        			  cudaEventRecord(stop_INT_0, stream[0]);
        			#endif
  		
  		
  		// first FACE
		cudaStreamSynchronize(stream[1]);
		
  		MPI_Sendrecv(RAND1, 24*tSliceEO, MPI_FLOAT, g_nb_t_dn, 0,
  		             RAND3, 24*tSliceEO, MPI_FLOAT, g_nb_t_up, 0,
  		             g_cart_grid, &stat[0]);
  		
  		//MPI_Isend(RAND1, 24*tSliceEO, MPI_FLOAT, g_nb_t_dn, 0,
  		//          g_cart_grid, &send_req[0]);
  		//MPI_Recv(RAND3, 24*tSliceEO, MPI_FLOAT, g_nb_t_up, 0,
  		//         g_cart_grid, &stat[0]);
  		
  		//MPI_Wait(&recv_request1, &stat[0]);
  		
  		cudaMemcpyAsync(spinin+6*VolumeEO           , RAND3, tSliceEO*6*sizeof(float4), cudaMemcpyHostToDevice, stream[1]);
  		
  				#ifdef ASYNC_TIMING
  		  		  cudaEventRecord(stop_H2D_3, stream[1]);
  		  		#endif
  		
  		
  		dev_Hopping_Matrix_ASYNC <<<gridsize2, blocksize, 0, stream[1]>>> ( gf,
  		                                                                    spinin, spinout,
  		                                                                    gfindex_site, gfindex_nextsite, nn_evenodd,
  		                                                                    ieo,
  		                                                                    VolumeEO-tSliceEO, tSliceEO );
  				#ifdef ASYNC_TIMING
  		  		  cudaEventRecord(stop_EXT_1, stream[1]);
  		  		#endif
  		
  		
  		// second FACE
  		cudaStreamSynchronize(stream[2]);
		
  		MPI_Sendrecv(RAND2, 24*tSliceEO, MPI_FLOAT, g_nb_t_up, 1,
  		             RAND4, 24*tSliceEO, MPI_FLOAT, g_nb_t_dn, 1,
  		             g_cart_grid, &stat[1]);
  		
  		//MPI_Isend(RAND2, 24*tSliceEO, MPI_FLOAT, g_nb_t_up, 1,
  		//          g_cart_grid, &send_req[1]);
  		//MPI_Recv(RAND4, 24*tSliceEO, MPI_FLOAT, g_nb_t_dn, 1,
  		//         g_cart_grid, &stat[1]);
  		
  		//MPI_Wait(&recv_request2, &stat[1]);
  		
  		cudaMemcpyAsync(spinin+6*(VolumeEO+tSliceEO), RAND4, tSliceEO*6*sizeof(float4), cudaMemcpyHostToDevice, stream[2]);
  		
  				#ifdef ASYNC_TIMING
  		  		  cudaEventRecord(stop_H2D_4, stream[2]);
  		  		#endif
		
		
  		dev_Hopping_Matrix_ASYNC <<<gridsize2, blocksize, 0, stream[2]>>> ( gf,
  		                                                                    spinin, spinout,
  		                                                                    gfindex_site, gfindex_nextsite, nn_evenodd,
  		                                                                    ieo,
  		                                                                    0, tSliceEO );
  				#ifdef ASYNC_TIMING
  		  		  cudaEventRecord(stop_EXT_2, stream[2]);
  		  		#endif
  		
  		
  				#ifdef ASYNC_TIMING
  				  cudaEventRecord(stop_ALL, 0);
  				#endif
  		*/
  		
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  		
  		
  		cudaThreadSynchronize();		// test if needed
  
  
  
  #endif
  
  
  #ifdef USETEXTURE
    unbind_texture_spin(1);
  #endif
  
  
}






///////////////////////////
// MATRIX MULTIPLICATION //
///////////////////////////

// the GPU implementation of  Q_Qdagger_ND(...)  from Nondegenerate_Matrix.c
//	Flo's equivalent function for the standard and non-nd case is  dev_Qtm_pm_psi

void matrix_multiplication32_mpi_ASYNC (dev_spinor * spinout_up, dev_spinor * spinout_dn,
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
  
    	  
  HOPPING_ASYNC(dev_gf, spinin_dn, dev_spin_eo1_up, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0, gridsize1, blocksize1);
  
  HOPPING_ASYNC(dev_gf, spinin_up, dev_spin_eo1_dn, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0, gridsize1, blocksize1);

  
  // imubar, gamma5
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo1_up, dev_spin_eo2_up, -1.0);		// dev_spin_eo2_up  =  (1 - imubar) * dev_spin_eo1_up  =  (1 - imubar)*(M_eo) * spinin_dn
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo1_dn, dev_spin_eo2_dn, +1.0);		// dev_spin_eo2_dn  =  (1 + imubar) * dev_spin_eo1_dn  =  (1 + imubar)*(M_eo) * spinin_up
  
  
  // linear algebra
  cublasSaxpy (N_floats, g_epsbar, (float *) dev_spin_eo1_dn, 1, (float *) dev_spin_eo2_up, 1);
  										// dev_spin_eo2_up  =  dev_spin_eo2_up  +  epsbar * dev_spin_eo1_dn  =  (1 - imubar)*(M_eo) * spinin_dn  +  epsbar * (M_eo) * spinin_up
  cublasSaxpy (N_floats, g_epsbar, (float *) dev_spin_eo1_up, 1, (float *) dev_spin_eo2_dn, 1);
  										// dev_spin_eo2_dn  =  dev_spin_eo2_dn  +  epsbar * dev_spin_eo1_up  =  (1 + imubar)*(M_eo) * spinin_up  +  epsbar * (M_eo) * spinin_dn
  
  // linear algebra
  cublasSscal (N_floats, nrm, (float *) dev_spin_eo2_up, 1);			// dev_spin_eo2_up  =  nrm * dev_spin_eo2_up  =  nrm*(1-imubar)*(M_eo) * spinin_dn  +  nrm*epsbar*(M_eo) * spinin_up
  
  cublasSscal (N_floats, nrm, (float *) dev_spin_eo2_dn, 1);			// dev_spin_eo2_dn  =  nrm * dev_spin_eo2_dn  =  nrm*(1+imubar)*(M_eo) * spinin_up  +  nrm*epsbar*(M_eo) * spinin_dn
  
  
  
  HOPPING_ASYNC(dev_gf, dev_spin_eo2_up, dev_spin_eo1_up, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1, gridsize1, blocksize1);
  
  HOPPING_ASYNC(dev_gf, dev_spin_eo2_dn, dev_spin_eo1_dn, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1, gridsize1, blocksize1);

  
  
  
  
  ////////////
  // (M_oo) //
  ////////////
  
  
  // imubar, gamma5
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(spinin_dn, dev_spin_eo2_up, +1.0);			// dev_spin_eo2_up = (1 + imubar) * spinin_dn
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(spinin_up, dev_spin_eo2_dn, -1.0);			// dev_spin_eo2_dn = (1 - imubar) * spinin_up
  
  
  // linear algebra													// remember: this is (M_oo) * (spinin_dn, spinin_up):
  cublasSaxpy (N_floats, -g_epsbar, (float *) spinin_up, 1, (float *) dev_spin_eo2_up, 1);
  												// dev_spin_eo2_up  =  dev_spin_eo2_up - epsbar*spinin_up  =  (1+imubar)*spinin_dn - epsbar*spinin_up
  cublasSaxpy (N_floats, -g_epsbar, (float *) spinin_dn, 1, (float *) dev_spin_eo2_dn, 1);
  												// dev_spin_eo2_dn  =  dev_spin_eo2_dn - epsbar*spinin_dn  =  (1-imubar)*spinin_up - epsbar*spinin_dn
  
  
  
  ///////////////////////////////////////
  // (M_oo)  -  (M_oe) (Mee^-1) (M_eo) //
  ///////////////////////////////////////
  
  // linear algebra													// this is ((M_oo) - (M_oe)(Mee^-1)(M_eo)) * (spinin_dn, spinin_up):
  cublasSaxpy (N_floats, -1.0, (float *) dev_spin_eo1_up, 1, (float *) dev_spin_eo2_up, 1);
  							// dev_spin_eo2_up  =  dev_spin_eo2_up  -  dev_spin_eo1_up  =                      (1+imubar) * spinin_dn  -                    epsbar * spinin_up
  							//                                                             - (M_oe)*nrm*(1-imubar)*(M_eo) * spinin_dn  -  (M_oe)*nrm*epsbar*(M_eo) * spinin_up
  cublasSaxpy (N_floats, -1.0, (float *) dev_spin_eo1_dn, 1, (float *) dev_spin_eo2_dn, 1);
  							// dev_spin_eo2_dn  =  dev_spin_eo2_dn  -  dev_spin_eo1_dn  =                      (1-imubar) * spinin_up  -                    epsbar * spinin_dn
  							//                                                             - (M_oe)*nrm*(1+imubar)*(M_eo) * spinin_up  -  (M_oe)*nrm*epsbar*(M_eo) * spinin_dn
  
  
  ////////////
  // gamma5 //
  ////////////
  
  // gamma5
  dev_gamma5<<<gridsize3, blocksize3>>>(dev_spin_eo2_up, dev_spin_eo2_up);					// dev_spin_eo2_up = gamma5 * dev_spin_eo2_up 
  dev_gamma5<<<gridsize3, blocksize3>>>(dev_spin_eo2_dn, dev_spin_eo2_dn);					// dev_spin_eo2_dn = gamma5 * dev_spin_eo2_dn
  	
  
  
  
  
  
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
  
  
  // linear algebra
  cublasSaxpy (N_floats, g_epsbar, (float *) dev_spin_eo1_dn, 1, (float *) dev_spin_eo2_up, 1);
  										// dev_spin_eo2_up  =  dev_spin_eo2_up  +  epsbar * dev_spin_eo1_dn  =  (1 - imubar)*(M_eo) * dev_spin_eo3_up  +  epsbar * (M_eo) * dev_spin_eo3_dn
  cublasSaxpy (N_floats, g_epsbar, (float *) dev_spin_eo1_up, 1, (float *) dev_spin_eo2_dn, 1);
  										// dev_spin_eo2_dn  =  dev_spin_eo2_dn  +  epsbar * dev_spin_eo1_up  =  (1 + imubar)*(M_eo) * dev_spin_eo3_dn  +  epsbar * (M_eo) * dev_spin_eo3_up
  
  // lineare algebra
  cublasSscal (N_floats, nrm, (float *) dev_spin_eo2_up, 1);			// dev_spin_eo2_up  =  nrm * dev_spin_eo2_up  =  nrm*(1-imubar)*(M_eo) * dev_spin_eo3_up  +  nrm*epsbar*(M_eo) * dev_spin_eo3_dn
  
  cublasSscal (N_floats, nrm, (float *) dev_spin_eo2_dn, 1);			// dev_spin_eo2_dn  =  nrm * dev_spin_eo2_dn  =  nrm*(1+imubar)*(M_eo) * dev_spin_eo3_dn  +  nrm*epsbar*(M_eo) * dev_spin_eo3_up
  
    
    	  
  HOPPING_ASYNC(dev_gf, dev_spin_eo2_up, dev_spin_eo1_up, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1, gridsize1, blocksize1);

  HOPPING_ASYNC(dev_gf, dev_spin_eo2_dn, dev_spin_eo1_dn, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1, gridsize1, blocksize1);

  
  
  
  ////////////
  // (M_oo) //
  ////////////
  
  // imubar, gamma5
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo3_up, dev_spin_eo2_up, +1.0);				// dev_spin_eo2_up = (1 + imubar) * dev_spin_eo3_up
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo3_dn, dev_spin_eo2_dn, -1.0);				// dev_spin_eo2_dn = (1 - imubar) * dev_spin_eo3_dn
  
  
  // lineare algebra										// remember: this is (M_oo) * (dev_spin_eo3_up, dev_spin_eo3_dn):
  cublasSaxpy (N_floats, -g_epsbar, (float *) dev_spin_eo3_dn, 1, (float *) dev_spin_eo2_up, 1);
  												// dev_spin_eo2_up  =  dev_spin_eo2_up - epsbar*dev_spin_eo3_dn  =  (1+imubar)*dev_spin_eo3_up - epsbar*dev_spin_eo3_dn
  cublasSaxpy (N_floats, -g_epsbar, (float *) dev_spin_eo3_up, 1, (float *) dev_spin_eo2_dn, 1);
  												// dev_spin_eo2_dn  =  dev_spin_eo2_dn - epsbar*dev_spin_eo3_up  =  (1-imubar)*dev_spin_eo3_dn - epsbar*dev_spin_eo3_up
  
  
  
  ///////////////////////////////////////
  // (M_oo)  -  (M_oe) (Mee^-1) (M_eo) //
  ///////////////////////////////////////
  
  // lineare algebra										// this is ( (M_oo) - (M_oe) (Mee^-1) (M_eo) ) * (dev_spin_eo3_up, dev_spin_eo3_dn)
  cublasSaxpy (N_floats, -1.0, (float *) dev_spin_eo1_up, 1, (float *) dev_spin_eo2_up, 1);
  							// dev_spin_eo2_up  =  dev_spin_eo2_up  -  dev_spin_eo1_up  =                      (1+imubar) * dev_spin_eo3_up  -                    epsbar * dev_spin_eo3_dn
  							//                                                             - (M_oe)*nrm*(1-imubar)*(M_eo) * dev_spin_eo3_up  -  (M_oe)*nrm*epsbar*(M_eo) * dev_spin_eo3_dn
  cublasSaxpy (N_floats, -1.0, (float *) dev_spin_eo1_dn, 1, (float *) dev_spin_eo2_dn, 1);
  							// dev_spin_eo2_dn  =  dev_spin_eo2_dn  -  dev_spin_eo1_dn  =                      (1-imubar) * dev_spin_eo3_dn  -                    epsbar * dev_spin_eo3_up
  							//                                                             - (M_oe)*nrm*(1+imubar)*(M_eo) * dev_spin_eo3_dn  -  (M_oe)*nrm*epsbar*(M_eo) * dev_spin_eo3_up
  
  
  ////////////
  // gamma5 //
  ////////////
  
  // gamma5
  dev_gamma5<<<gridsize3, blocksize3>>>(dev_spin_eo2_up, dev_spin_eo2_up);					// dev_spin_eo2_up = gamma5 * dev_spin_eo2_up 
  dev_gamma5<<<gridsize3, blocksize3>>>(dev_spin_eo2_dn, dev_spin_eo2_dn);					// dev_spin_eo2_dn = gamma5 * dev_spin_eo2_dn
  
  
  
  
  /*
  ////////////
  // output //										// output is already done by setting  dev_spin_eo2_up/dn = spinout_up/dn
  ////////////
  
  dev_copy_spinor_field<<<griddim2, blockdim2 >>>(dev_spin_eo2_up, spinout_up);		// spinout_up = dev_spin_eo2_up
  dev_copy_spinor_field<<<griddim2, blockdim2 >>>(dev_spin_eo2_dn, spinout_dn);		// spinout_dn = dev_spin_eo2_dn
  */
  
  
  return;
  
}//matrix_multiplication32_mpi_ASYNC()





















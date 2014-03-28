/**************************************************************************
 *
 * Copyright (C) 2010 Joseph Nagel
 *               2010 Florian Burger
 *               
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
 
 
 



#ifndef HALF


//this gathers spinors in the space-time first indexing to a continuous block that can be exchanged via mpi
//starts at start and goes to start+size
//other versions for relativistic basis
__global__ void dev_gather_rand(dev_spinor * sin, dev_spinor * rand, int start, int size){
  int pos,pos2;
  pos2 = threadIdx.x + blockDim.x * blockIdx.x;
  pos = start  + pos2;
  
  if(pos < start + size){
  
  rand[6*pos2].x = sin[pos+0*DEVOFF].x; 
  rand[6*pos2].y = sin[pos+0*DEVOFF].y; 
  rand[6*pos2].z = sin[pos+0*DEVOFF].z;
  rand[6*pos2].w = sin[pos+0*DEVOFF].w;
  
  rand[6*pos2+1].x = sin[pos+1*DEVOFF].x; 
  rand[6*pos2+1].y = sin[pos+1*DEVOFF].y; 
  rand[6*pos2+1].z = sin[pos+1*DEVOFF].z;
  rand[6*pos2+1].w = sin[pos+1*DEVOFF].w; 
  
  rand[6*pos2+2].x = sin[pos+2*DEVOFF].x; 
  rand[6*pos2+2].y = sin[pos+2*DEVOFF].y; 
  rand[6*pos2+2].z = sin[pos+2*DEVOFF].z;
  rand[6*pos2+2].w = sin[pos+2*DEVOFF].w;  
  
  rand[6*pos2+3].x = sin[pos+3*DEVOFF].x; 
  rand[6*pos2+3].y = sin[pos+3*DEVOFF].y; 
  rand[6*pos2+3].z = sin[pos+3*DEVOFF].z;
  rand[6*pos2+3].w = sin[pos+3*DEVOFF].w;   
  
  rand[6*pos2+4].x = sin[pos+4*DEVOFF].x; 
  rand[6*pos2+4].y = sin[pos+4*DEVOFF].y; 
  rand[6*pos2+4].z = sin[pos+4*DEVOFF].z;
  rand[6*pos2+4].w = sin[pos+4*DEVOFF].w; 
  
  rand[6*pos2+5].x = sin[pos+5*DEVOFF].x; 
  rand[6*pos2+5].y = sin[pos+5*DEVOFF].y; 
  rand[6*pos2+5].z = sin[pos+5*DEVOFF].z;
  rand[6*pos2+5].w = sin[pos+5*DEVOFF].w; 
  
  }
}



//this gathers spinors in the space-time first indexing to a continuous block that can be exchanged via mpi
//starts at start and goes to start+size
//RELATIVISTIC_BASIS!!
__global__ void dev_gather_rand_relup(dev_spinor * sin, dev_spinor * rand, int start, int size){
  int pos, pos2;
  pos2 = threadIdx.x + blockDim.x * blockIdx.x;
  pos = start  +  pos2;
  
  if(pos < start + size){
//only fetch upper spinor and store consecutively (NOTE the "3" instead of "6") 
  rand[3*pos2].x = sin[pos+0*DEVOFF].x; 
  rand[3*pos2].y = sin[pos+0*DEVOFF].y; 
  rand[3*pos2].z = sin[pos+0*DEVOFF].z;
  rand[3*pos2].w = sin[pos+0*DEVOFF].w;
  
  rand[3*pos2+1].x = sin[pos+1*DEVOFF].x; 
  rand[3*pos2+1].y = sin[pos+1*DEVOFF].y; 
  rand[3*pos2+1].z = sin[pos+1*DEVOFF].z;
  rand[3*pos2+1].w = sin[pos+1*DEVOFF].w; 
  
  rand[3*pos2+2].x = sin[pos+2*DEVOFF].x; 
  rand[3*pos2+2].y = sin[pos+2*DEVOFF].y; 
  rand[3*pos2+2].z = sin[pos+2*DEVOFF].z;
  rand[3*pos2+2].w = sin[pos+2*DEVOFF].w;  
  
  }
}


//this gathers spinors in the space-time first indexing to a continuous block that can be exchanged via mpi
//starts at start and goes to start+size
//RELATIVISTIC_BASIS!!
__global__ void dev_gather_rand_reldn(dev_spinor * sin, dev_spinor * rand, int start, int size){
  int pos,pos2;
  pos2 = threadIdx.x + blockDim.x * blockIdx.x;
  pos = start  +  pos2;
  
  if(pos < start + size){
//only fetch upper spinor and store consecutively (NOTE the "3" instead of "6") 
  rand[3*pos2].x = sin[pos+3*DEVOFF].x; 
  rand[3*pos2].y = sin[pos+3*DEVOFF].y; 
  rand[3*pos2].z = sin[pos+3*DEVOFF].z;
  rand[3*pos2].w = sin[pos+3*DEVOFF].w;
  
  rand[3*pos2+1].x = sin[pos+4*DEVOFF].x; 
  rand[3*pos2+1].y = sin[pos+4*DEVOFF].y; 
  rand[3*pos2+1].z = sin[pos+4*DEVOFF].z;
  rand[3*pos2+1].w = sin[pos+4*DEVOFF].w; 
  
  rand[3*pos2+2].x = sin[pos+5*DEVOFF].x; 
  rand[3*pos2+2].y = sin[pos+5*DEVOFF].y; 
  rand[3*pos2+2].z = sin[pos+5*DEVOFF].z;
  rand[3*pos2+2].w = sin[pos+5*DEVOFF].w;  
  
  }
}





//this spreads spinors to the space-time first indexing, rand is passed as a continuous block that is exchanged via mpi
//starts at start and goes to start+size
//other versions for relativistic basis
__global__ void dev_spread_rand(dev_spinor * sin, dev_spinor * rand, int start, int size){
  int pos, pos2;
  pos2 = threadIdx.x + blockDim.x * blockIdx.x;
  pos = start  + pos2; 
  
  if(pos < start + size){
  sin[pos+0*DEVOFF].x = rand[6*pos2].x; 
  sin[pos+0*DEVOFF].y = rand[6*pos2].y; 
  sin[pos+0*DEVOFF].z = rand[6*pos2].z;
  sin[pos+0*DEVOFF].w = rand[6*pos2].w;
  
  sin[pos+1*DEVOFF].x = rand[6*pos2+1].x; 
  sin[pos+1*DEVOFF].y = rand[6*pos2+1].y;
  sin[pos+1*DEVOFF].z = rand[6*pos2+1].z;
  sin[pos+1*DEVOFF].w = rand[6*pos2+1].w; 
  
  sin[pos+2*DEVOFF].x = rand[6*pos2+2].x; 
  sin[pos+2*DEVOFF].y = rand[6*pos2+2].y; 
  sin[pos+2*DEVOFF].z = rand[6*pos2+2].z;
  sin[pos+2*DEVOFF].w = rand[6*pos2+2].w;  
  
  sin[pos+3*DEVOFF].x = rand[6*pos2+3].x; 
  sin[pos+3*DEVOFF].y = rand[6*pos2+3].y; 
  sin[pos+3*DEVOFF].z = rand[6*pos2+3].z;
  sin[pos+3*DEVOFF].w = rand[6*pos2+3].w;   
  
  sin[pos+4*DEVOFF].x = rand[6*pos2+4].x; 
  sin[pos+4*DEVOFF].y = rand[6*pos2+4].y; 
  sin[pos+4*DEVOFF].z = rand[6*pos2+4].z;
  sin[pos+4*DEVOFF].w = rand[6*pos2+4].w; 
  
  sin[pos+5*DEVOFF].x = rand[6*pos2+5].x; 
  sin[pos+5*DEVOFF].y = rand[6*pos2+5].y; 
  sin[pos+5*DEVOFF].z = rand[6*pos2+5].z;
  sin[pos+5*DEVOFF].w = rand[6*pos2+5].w; 
  
  }
}



//this spreads spinors to the space-time first indexing, rand is passed as a continuous block that is exchanged via mpi
//starts at start and goes to start+size
//RELATIVISTIC_BASIS !!!
__global__ void dev_spread_rand_relup(dev_spinor * sin, dev_spinor * rand, int start, int size){
  int pos, pos2;
  pos2 = threadIdx.x + blockDim.x * blockIdx.x;
  pos = start  + pos2; 
  
  if(pos < start + size){
//only spread upper spinor and store consecutively (NOTE the "3" instead of "6") 
  sin[pos+0*DEVOFF].x = rand[3*pos2].x; 
  sin[pos+0*DEVOFF].y = rand[3*pos2].y; 
  sin[pos+0*DEVOFF].z = rand[3*pos2].z;
  sin[pos+0*DEVOFF].w = rand[3*pos2].w;
  
  sin[pos+1*DEVOFF].x = rand[3*pos2+1].x; 
  sin[pos+1*DEVOFF].y = rand[3*pos2+1].y;
  sin[pos+1*DEVOFF].z = rand[3*pos2+1].z;
  sin[pos+1*DEVOFF].w = rand[3*pos2+1].w; 
  
  sin[pos+2*DEVOFF].x = rand[3*pos2+2].x; 
  sin[pos+2*DEVOFF].y = rand[3*pos2+2].y; 
  sin[pos+2*DEVOFF].z = rand[3*pos2+2].z;
  sin[pos+2*DEVOFF].w = rand[3*pos2+2].w;  
  
  }
}



//this spreads spinors to the space-time first indexing, rand is passed as a continuous block that is exchanged via mpi
//starts at start and goes to start+size
//RELATIVISTIC_BASIS !!!
__global__ void dev_spread_rand_reldn(dev_spinor * sin, dev_spinor * rand, int start, int size){
  int pos,pos2;
  pos2 = threadIdx.x + blockDim.x * blockIdx.x;
  pos = start  +  pos2;
  
  if(pos < start + size){
//only spread lower spinor and store consecutively (NOTE the "3" instead of "6")  
  
  sin[pos+3*DEVOFF].x = rand[3*pos2].x; 
  sin[pos+3*DEVOFF].y = rand[3*pos2].y; 
  sin[pos+3*DEVOFF].z = rand[3*pos2].z;
  sin[pos+3*DEVOFF].w = rand[3*pos2].w;   
  
  sin[pos+4*DEVOFF].x = rand[3*pos2+1].x; 
  sin[pos+4*DEVOFF].y = rand[3*pos2+1].y; 
  sin[pos+4*DEVOFF].z = rand[3*pos2+1].z;
  sin[pos+4*DEVOFF].w = rand[3*pos2+1].w; 
  
  sin[pos+5*DEVOFF].x = rand[3*pos2+2].x; 
  sin[pos+5*DEVOFF].y = rand[3*pos2+2].y; 
  sin[pos+5*DEVOFF].z = rand[3*pos2+2].z;
  sin[pos+5*DEVOFF].w = rand[3*pos2+2].w; 
  
  }
}









///////////////////////
// boundary exchange //
///////////////////////

#ifdef MPI

// both versions do work:


/*
// copies VOLUME to host, exchanges, copies RAND back to device

void xchange_field_wrapper (dev_spinor * dev_spin, int ieo) {

  size_t size_Volume = VOLUME/2 * 6*sizeof(dev_spinor);
  size_t size_Rand   = RAND/2   * 6*sizeof(dev_spinor);

  to_host_mpi(spinor_xchange, dev_spin, h2d_spin_up, size_Volume, 0, VOLUME/2);
  xchange_field(spinor_xchange, ieo);
  to_device_mpi(dev_spin, spinor_xchange, h2d_spin_dn, size_Rand, VOLUME/2, (VOLUME+RAND)/2);

}
*/




// copies the boundary t-slices t=0 and t=T-1 to host		// will be used in dev_Qtm_pm_ndpsi_mpi(), not ASYNC
//	exchanges						// provides a wrapped version of Carsten's xchange_field()
//		copies RAND back to device			//	and not asynchronous version of ASYNC.cuh

void xchange_field_wrapper (dev_spinor * dev_spin, int ieo) {
  
    
    int tSliceEO = LX*LY*LZ/2;
    int VolumeEO = VOLUME/2;
    
    //this is the same partitioning as for dev_mul_one_pm...
    int gridsize;
    int blockdim = BLOCK2;
    if( tSliceEO % blockdim == 0){
      gridsize = (int) tSliceEO/blockdim;
    }
    else{
      gridsize = (int) tSliceEO/blockdim + 1;
    }
    int griddim = gridsize;
    
    
    #ifdef RELATIVISTIC_BASIS
      //this goes backwards, so for the receiver this is forward
      dev_gather_rand_relup<<<griddim, blockdim >>>(dev_spin,RAND_BW,0,tSliceEO);
      cudaMemcpy(R1, RAND_BW, tSliceEO*3*sizeof(float4), cudaMemcpyDeviceToHost);
    
      //this goes forward, so for the receiver this is backward
      dev_gather_rand_reldn<<<griddim, blockdim >>>(dev_spin,RAND_FW,(VolumeEO-tSliceEO),tSliceEO);
      cudaMemcpy(R2, RAND_FW, tSliceEO*3*sizeof(float4), cudaMemcpyDeviceToHost);    
    #else
      //this goes backwards
      dev_gather_rand<<<griddim, blockdim >>>(dev_spin,RAND_BW,0,tSliceEO);
      cudaMemcpy(R1, RAND_BW, tSliceEO*6*sizeof(float4), cudaMemcpyDeviceToHost);
    
      //this goes forward
      dev_gather_rand<<<griddim, blockdim >>>(dev_spin,RAND_FW,(VolumeEO-tSliceEO),tSliceEO);
      cudaMemcpy(R2, RAND_FW, tSliceEO*6*sizeof(float4), cudaMemcpyDeviceToHost);
    #endif
    
    //we only need to exchange half of the spinors in relativistic basis (upper or lower part)
    int nfloat_per_spin;
    #ifdef RELATIVISTIC_BASIS
      nfloat_per_spin = 12;
    #else
      nfloat_per_spin = 24;
    #endif
    
    MPI_Sendrecv(R1, nfloat_per_spin*tSliceEO, MPI_FLOAT, g_nb_t_dn, 0,
                 R3, nfloat_per_spin*tSliceEO, MPI_FLOAT, g_nb_t_up, 0,
                 g_cart_grid, &stat[0]);
    MPI_Sendrecv(R2, nfloat_per_spin*tSliceEO, MPI_FLOAT, g_nb_t_up, 1,
                 R4, nfloat_per_spin*tSliceEO, MPI_FLOAT, g_nb_t_dn, 1,
                 g_cart_grid, &stat[1]);

    #ifdef RELATIVISTIC_BASIS
      cudaMemcpy(RAND_BW, R3, tSliceEO*3*sizeof(float4), cudaMemcpyHostToDevice);
      dev_spread_rand_relup<<<griddim, blockdim >>>(dev_spin,RAND_BW,VolumeEO,tSliceEO);
    
      cudaMemcpy(RAND_FW, R4, tSliceEO*3*sizeof(float4), cudaMemcpyHostToDevice);
      dev_spread_rand_reldn<<<griddim, blockdim >>>(dev_spin,RAND_FW,(VolumeEO+tSliceEO),tSliceEO);   
    #else 
      cudaMemcpy(RAND_BW, R3, tSliceEO*6*sizeof(float4), cudaMemcpyHostToDevice);
      dev_spread_rand<<<griddim, blockdim >>>(dev_spin,RAND_BW,VolumeEO,tSliceEO);
    
      cudaMemcpy(RAND_FW, R4, tSliceEO*6*sizeof(float4), cudaMemcpyHostToDevice);
      dev_spread_rand<<<griddim, blockdim >>>(dev_spin,RAND_FW,(VolumeEO+tSliceEO),tSliceEO);
    #endif

  
}

#endif // MPI












/*
This is the wrapper function for the device hopping matrix for MPI with support of
CUDA streams in order to parallelize bulk calculation and boundary exchange
*/

void HOPPING_ASYNC (dev_su3_2v * gf, 
                    dev_spinor * spinin, dev_spinor * spinout,
                    int * gfindex_site, int * gfindex_nextsite, int * nn_evenodd,
                    int ieo,
                    int gridsize, dim3 blocksize_in) {
  
  
  // for even/odd
  int tSliceEO = LX*LY*LZ/2;
  int VolumeEO = VOLUME/2;
  
  
  // gridsizes
  int gridsize1;
  int gridsize2;
  
    int blocksize = blocksize_in.x;
  
  #ifndef ASYNC_TSLICES
    int tSlices = 1; 
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
  #else
    int tSlices = ASYNC_TSLICES;
    if ( (VolumeEO-2*tSlices*tSliceEO) % blocksize == 0 ) {
      gridsize1  = (VolumeEO-2*tSlices*tSliceEO) / blocksize;
    }
    else {
      gridsize1  = (int) ( ((VolumeEO-2*tSlices*tSliceEO)/blocksize) + 1);
    }
    
    if ( (tSlices*tSliceEO) % blocksize == 0 ) {
      gridsize2  = (tSlices*tSliceEO) / blocksize;
    }
    else {
      gridsize2  = (int) ( ((tSlices*tSliceEO)/blocksize) + 1);
    }
  #endif
    //this is the same partitioning as for dev_mul_one_pm...
    int gridsize3;
    int blocksize3 = BLOCK2;
    if( tSliceEO % blocksize3 == 0){
      gridsize3 = (int) tSliceEO/blocksize3;
    }
    else{
      gridsize3 = (int) tSliceEO/blocksize3 + 1;
    } 
  
  
  
  #ifdef USETEXTURE
    bind_texture_spin(spinin,1);
  #endif
  
  
  
  #if ASYNC == 0		// primitive version
    		

 cudaError_t cudaerr;

 if((cudaerr=cudaGetLastError()) != cudaSuccess){
    printf("Error in ASYNC: %s\n", cudaGetErrorString(cudaerr));
    printf("gridsize = %d, blocksize = %d\n",gridsize1,blocksize);
    exit(200);
  } 
  
		// applies to the parts which don't need communication
  		dev_Hopping_Matrix <<<gridsize1, blocksize>>> ( gf,
        	                                                      spinin, spinout,
        	                                                      gfindex_site, gfindex_nextsite, nn_evenodd,
        	                                                      ieo,
        	                                                      //2*tSliceEO, VolumeEO-4*tSliceEO );
        	                                                       tSliceEO, VolumeEO-2*tSliceEO );
	                                                       


  		
  		
  		// exchanges the boundaries
  		xchange_field_wrapper(spinin, ieo);			// to be further optimized !!
  		
  		// applies the hopping matrix to remaining parts
  		dev_Hopping_Matrix <<<gridsize2, blocksize>>> ( gf,
  		                                                      spinin, spinout,
  		                                                      gfindex_site, gfindex_nextsite, nn_evenodd,
  		                                                      ieo,
  		                                                      //0, 2*tSliceEO );
  		                                                      0, tSliceEO );
  		
  		dev_Hopping_Matrix <<<gridsize2, blocksize>>> ( gf,
  		                                                      spinin, spinout,
  		                                                      gfindex_site, gfindex_nextsite, nn_evenodd,
  		                                                      ieo,
  		                                                      //VolumeEO-2*tSliceEO, 2*tSliceEO );
  		                                                      VolumeEO-tSliceEO, tSliceEO );	
		
		
		
		
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  
  
  
  #elif ASYNC == 1		// optimized version
                //set the amount of data we need to transfer with mpi
      		#ifdef RELATIVISTIC_BASIS
  		  int flperspin = 12;
  		#else
  		  int flperspin = 24;
  		#endif	
  		
				#ifdef ASYNC_TIMING
  				  cudaEventRecord(start_ALL, 0);
  				  mpi_start_ALL = MPI_Wtime();
  				#endif
        	
        	
        	// copies first FACE to host
        	#ifdef RELATIVISTIC_BASIS
        	  dev_gather_rand_relup<<<gridsize3, blocksize3 ,0,stream[1] >>>(spinin,RAND_BW,0,tSliceEO);
  		  cudaMemcpyAsync(RAND1, RAND_BW, tSliceEO*3*sizeof(float4), cudaMemcpyDeviceToHost, stream[1]);
  		#else
        	  dev_gather_rand<<<gridsize3, blocksize3,0,stream[1] >>>(spinin,RAND_BW,0,tSliceEO);
  		  cudaMemcpyAsync(RAND1, RAND_BW, tSliceEO*6*sizeof(float4), cudaMemcpyDeviceToHost, stream[1]);  		
		  //printf("g_proc_id = %d: R1  %e \n", g_proc_id, RAND1[11].x);
		#endif
  				#ifdef ASYNC_TIMING
  				  cudaEventRecord(stop_D2H_1, stream[1]);
  				#endif

//         	 if((cudaerr=cudaPeekAtLastError()) != cudaSuccess){
//                    printf("Error in ASYNC: %s\n", cudaGetErrorString(cudaerr));
//                    printf("gridsize = %d, blocksize = %d\n",gridsize1,blocksize);
//                    exit(200);
//                  } 


                // copies second FACE to host
  		#ifdef RELATIVISTIC_BASIS
  		  dev_gather_rand_reldn<<<gridsize3, blocksize3, 0, stream[2] >>>(spinin,RAND_FW,(VolumeEO-tSliceEO),tSliceEO);
  		  cudaMemcpyAsync(RAND2, RAND_FW, tSliceEO*3*sizeof(float4), cudaMemcpyDeviceToHost, stream[2]);
  		#else
                  dev_gather_rand<<<gridsize3, blocksize3, 0, stream[2] >>>(spinin,RAND_FW,(VolumeEO-tSliceEO),tSliceEO);
  		  cudaMemcpyAsync(RAND2, RAND_FW, tSliceEO*6*sizeof(float4), cudaMemcpyDeviceToHost, stream[2]);  		
		  //printf("g_proc_id = %d: R2  %e \n", g_proc_id, RAND2[11].x);
		#endif
  				#ifdef ASYNC_TIMING
  				  cudaEventRecord(stop_D2H_2, stream[2]);
  				#endif	
        	

        	
        	
//INTERNAL        	
        	// starts INTERNAL kernel
  		dev_Hopping_Matrix <<<gridsize1, blocksize, 0, stream[0]>>> ( gf,
        	                                                                    spinin, spinout,
        	                                                                    gfindex_site, gfindex_nextsite, nn_evenodd,
        	                                                                    ieo,
        	                                                                    tSlices*tSliceEO, VolumeEO-2*tSlices*tSliceEO );
        			#ifdef ASYNC_TIMING
        			  cudaEventRecord(stop_INT_0, stream[0]);
        			#endif
  		

//FIRST FACE  		
  		// exchanges first FACE
  		cudaStreamSynchronize(stream[1]);				// SYNCPOINT
  		
  				#ifdef ASYNC_TIMING
  				  mpi_start_sendrecv_1 = MPI_Wtime();
  				#endif
  		
  		//MPI_Irecv(RAND3, flperspin*tSliceEO, MPI_FLOAT, g_nb_t_up, 0,
  		//          g_cart_grid, &recv_req[0]);
  		//MPI_Isend(RAND1, flperspin*tSliceEO, MPI_FLOAT, g_nb_t_dn, 0,
  		//          g_cart_grid, &send_req[0]);
  		
  		MPI_Sendrecv(RAND1, flperspin*tSliceEO, MPI_FLOAT, g_nb_t_dn, 0,	// SYNCPOINT
  		             RAND3, flperspin*tSliceEO, MPI_FLOAT, g_nb_t_up, 0,
  		             g_cart_grid, &stat[0]);
  		
  				#ifdef ASYNC_TIMING
  				  mpi_stop_sendrecv_1 = MPI_Wtime();
  				#endif
  		
  		
  		// copies first FACE back to device				
  		
  		
  		// order may switched
  		//MPI_Wait(&recv_req[0], &stat[0]);									
  		// synchronous
                #ifdef RELATIVISTIC_BASIS 
  		  cudaMemcpyAsync(RAND_BW, RAND3, tSliceEO*3*sizeof(float4), cudaMemcpyHostToDevice, stream[1]);
  		  dev_spread_rand_relup<<<gridsize3, blocksize3, 0, stream[1] >>>(spinin,RAND_BW,VolumeEO,tSliceEO);
                #else
                  //printf("g_proc_id = %d:  R3 %e \n", g_proc_id, RAND3[11].x);
                  cudaMemcpyAsync(RAND_BW, RAND3, tSliceEO*6*sizeof(float4), cudaMemcpyHostToDevice, stream[1]);
                  dev_spread_rand<<<gridsize3, blocksize3, 0, stream[1] >>>(spinin,RAND_BW,VolumeEO,tSliceEO);
                #endif

  				#ifdef ASYNC_TIMING
  				  cudaEventRecord(stop_H2D_3, stream[1]);
  				#endif
  		
  		
  		// applies first FACE
  		dev_Hopping_Matrix <<<gridsize2, blocksize, 0, stream[1]>>> ( gf,
  		                                                                      spinin, spinout,
  		                                                                      gfindex_site, gfindex_nextsite, nn_evenodd,
  		                                                                      ieo,
  		                                                                      VolumeEO-tSlices*tSliceEO, tSlices*tSliceEO );
  				#ifdef ASYNC_TIMING
  				  cudaEventRecord(stop_EXT_1, stream[1]);
  				#endif
  		
	
  		
 //SECOND FACE  	 		
  		// exchanges second FACE
  		cudaStreamSynchronize(stream[2]);				// SYNCPOINT
  		
  				#ifdef ASYNC_TIMING
  				  mpi_start_sendrecv_2 = MPI_Wtime();
  				#endif
  			
  		//MPI_Irecv(RAND4, flperspin*tSliceEO, MPI_FLOAT, g_nb_t_dn, 1,
  		//          g_cart_grid, &recv_req[1]);
  		//MPI_Isend(RAND2, flperspin*tSliceEO, MPI_FLOAT, g_nb_t_up, 1,
  		//          g_cart_grid, &send_req[1]);
  		
  		MPI_Sendrecv(RAND2, flperspin*tSliceEO, MPI_FLOAT, g_nb_t_up, 1,	// SYNCPOINT
  		             RAND4, flperspin*tSliceEO, MPI_FLOAT, g_nb_t_dn, 1,
  		             g_cart_grid, &stat[1]);
  		
  				#ifdef ASYNC_TIMING
  				  mpi_stop_sendrecv_2 = MPI_Wtime();
  				#endif
  		
  		// copies second FACE back to device
  		//MPI_Wait(&recv_req[1], &stat[1]); 
  		#ifdef RELATIVISTIC_BASIS 
  		  cudaMemcpyAsync(RAND_FW, RAND4, tSliceEO*3*sizeof(float4), cudaMemcpyHostToDevice, stream[2]);
  		  dev_spread_rand_reldn<<<gridsize3, blocksize3, 0, stream[2] >>>(spinin,RAND_FW,VolumeEO+tSliceEO,tSliceEO);
  		#else
  		  //printf("g_proc_id = %d: R4  %e \n", g_proc_id, RAND4[11].x);
   		  cudaMemcpyAsync(RAND_FW, RAND4, tSliceEO*6*sizeof(float4), cudaMemcpyHostToDevice, stream[2]);
  		  dev_spread_rand<<<gridsize3, blocksize3, 0, stream[2] >>>(spinin,RAND_FW,VolumeEO+tSliceEO,tSliceEO); 		
  		#endif
		
  				#ifdef ASYNC_TIMING
  				  cudaEventRecord(stop_H2D_4, stream[2]);
  				#endif
  		
  		// applies second FACE
  		dev_Hopping_Matrix <<<gridsize2, blocksize, 0, stream[2]>>> ( gf,
  		                                                                    spinin, spinout,
  		                                                                    gfindex_site, gfindex_nextsite, nn_evenodd,
  		                                                                    ieo,
  		                                                                    0, tSlices*tSliceEO );
  				#ifdef ASYNC_TIMING
  				  cudaEventRecord(stop_EXT_2, stream[2]);
  				#endif
  		
		
  		//done
  				#ifdef ASYNC_TIMING
  				  cudaEventRecord(stop_ALL, 0);
  				#endif
  		
  		
  
  
  #endif		// different optimized and non-optimized version
  
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  		
  		
  		cudaThreadSynchronize();		// test if needed	YES IS NEEDED according to Programming Guide
  
  //exit(200);
  #ifdef USETEXTURE
    unbind_texture_spin(1);
  #endif
  
  
}






#else //HALF

void matrix_multiplication32_mpi_ASYNC (dev_spinor * spinout_up, dev_spinor * spinout_dn,
                                        dev_spinor * spinin_up , dev_spinor * spinin_dn ,
                                        int gridsize1, int blocksize1, int gridsize2, int blocksize2,
                                        int gridsize3, int blocksize3, int gridsize4, int blocksize4) {
                                        
  printf("Warning: 'matrix_multiplication32_mpi_ASYNC' has been called from HALF code part. Not impemented yet. Aborting...\n");
  exit(200);                                        
}



void HOPPING_HALF_ASYNC (dev_su3_2v_half * gf, 
                    dev_spinor_half * spinin, float* spinin_norm, dev_spinor_half * spinout,
                    float* spinout_norm, int * gfindex_site, int * gfindex_nextsite, 
                    int * nn_evenodd, int ieo,
                    int gridsize, int blocksize) {
  
  
  // for even/odd
  int tSliceEO = LX*LY*LZ/2;
  int VolumeEO = VOLUME/2;
  
  #if defined ASYNC_OPTIMIZED && ASYNC == 3
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
  
  #ifndef ASYNC_TSLICES
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
  #else
    int tSlices = ASYNC_TSLICES;
    if ( (VolumeEO-2*tSlices*tSliceEO) % blocksize == 0 ) {
      gridsize1  = (VolumeEO-2*tSlices*tSliceEO) / blocksize;
    }
    else {
      gridsize1  = (int) ( ((VolumeEO-2*tSlices*tSliceEO)/blocksize) + 1);
    }
    
    if ( (tSlices*tSliceEO) % blocksize == 0 ) {
      gridsize2  = (tSlices*tSliceEO) / blocksize;
    }
    else {
      gridsize2  = (int) ( ((tSlices*tSliceEO)/blocksize) + 1);
    }
  #endif
  
  
  
  
  #ifdef USETEXTURE
    bind_halfspinor_texture(spinin, spinin_norm);
  #endif
  
  
  
  
  #if ASYNC == 0		// primitive version
  
  
  /*
  
		// applies to the parts which don't need communication
  		dev_Hopping_Matrix_half_ASYNC <<<gridsize1, blocksize>>> ( gf,
        	                                                      spinin, spinout,
        	                                                      gfindex_site, gfindex_nextsite, nn_evenodd,
        	                                                      ieo,
        	                                                      //2*tSliceEO, VolumeEO-4*tSliceEO );
        	                                                      tSliceEO, VolumeEO-2*tSliceEO );
  		
  		// exchanges the boundaries
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
		
    */		
		
		
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  
  
  
  #elif ASYNC == 1		// optimized version
  
  
				#ifdef ASYNC_TIMING
  				  cudaEventRecord(start_ALL, 0);
  				  mpi_start_ALL = MPI_Wtime();
  				#endif
        	
        	
        	// copies first FACE to host
  		cudaMemcpyAsync(RAND1, spinin, tSliceEO*6*sizeof(short4), cudaMemcpyDeviceToHost, stream[1]);
  		cudaMemcpyAsync(RAND1_norm, spinin_norm, tSliceEO*sizeof(float), cudaMemcpyDeviceToHost, stream[1]);
  				#ifdef ASYNC_TIMING
  				  cudaEventRecord(stop_D2H_1, stream[1]);
  				#endif
        	
        	
        	// INTERNAL kernel
  		dev_Hopping_Matrix_half_ASYNC <<<gridsize1, blocksize, 0, stream[0]>>> ( gf,
        	                                                       spinin, spinin_norm, 
        	                                                       spinout, spinout_norm,   
        	                                                       gfindex_site, gfindex_nextsite, 
        	                                                       nn_evenodd, ieo,
        	                                                       tSlices*tSliceEO,
        	                                                       VolumeEO-2*tSlices*tSliceEO );
        			#ifdef ASYNC_TIMING
        			  cudaEventRecord(stop_INT_0, stream[0]);
        			#endif
  		
  		
  		// exchanges first FACE
  		cudaStreamSynchronize(stream[1]);				// SYNCPOINT
  		
  				#ifdef ASYNC_TIMING
  				  mpi_start_sendrecv_1 = MPI_Wtime();
  				#endif
  		
  		
  		// copies second FACE to host
  		cudaMemcpyAsync(RAND2, spinin+6*(VolumeEO-tSliceEO), tSliceEO*6*sizeof(short4), cudaMemcpyDeviceToHost, stream[2]);
  		cudaMemcpyAsync(RAND2_norm, spinin_norm+(VolumeEO-tSliceEO), tSliceEO*sizeof(float), cudaMemcpyDeviceToHost, stream[2]);
  				#ifdef ASYNC_TIMING
  				  cudaEventRecord(stop_D2H_2, stream[2]);
  				#endif
  		
  		
  		//MPI_Irecv(RAND3, 24*tSliceEO, MPI_FLOAT, g_nb_t_up, 0,
  		//          g_cart_grid, &recv_req[0]);
  		//MPI_Isend(RAND1, 24*tSliceEO, MPI_FLOAT, g_nb_t_dn, 0,
  		//          g_cart_grid, &send_req[0]);
  		
  		MPI_Sendrecv(RAND1, 24*tSliceEO, MPI_SHORT, g_nb_t_dn, 0,	// SYNCPOINT
  		             RAND3, 24*tSliceEO, MPI_SHORT, g_nb_t_up, 0,
  		             g_cart_grid, &stat[0]);
  		// send norm             
  		MPI_Sendrecv(RAND1_norm, tSliceEO, MPI_FLOAT, g_nb_t_dn, 0,	// SYNCPOINT
  		             RAND3_norm, tSliceEO, MPI_FLOAT, g_nb_t_up, 0,
  		             g_cart_grid, &stat[0]);
  				#ifdef ASYNC_TIMING
  				  mpi_stop_sendrecv_1 = MPI_Wtime();
  				#endif
  		
  		
  		// copies first FACE back to device												// order may switched
  		//MPI_Wait(&recv_req[0], &stat[0]);												// synchronous
  		cudaMemcpyAsync(spinin+6*VolumeEO, RAND3, tSliceEO*6*sizeof(short4), cudaMemcpyHostToDevice, stream[1]);
  		cudaMemcpyAsync(spinin_norm+VolumeEO, RAND3_norm, tSliceEO*sizeof(float), cudaMemcpyHostToDevice, stream[1]);
  				#ifdef ASYNC_TIMING
  				  cudaEventRecord(stop_H2D_3, stream[1]);
  				#endif
  		
  		
  		// applies first FACE
  		dev_Hopping_Matrix_half_ASYNC <<<gridsize2, blocksize, 0, stream[1]>>> ( gf,
  		                                                                      spinin, spinin_norm,
  		                                                                      spinout, spinout_norm,
  		                                                                      gfindex_site, gfindex_nextsite,
  		                                                                      nn_evenodd, ieo,
  		                                                                      VolumeEO-tSlices*tSliceEO,
  		                                                                      tSlices*tSliceEO );
  				#ifdef ASYNC_TIMING
  				  cudaEventRecord(stop_EXT_1, stream[1]);
  				#endif
  		
  		
  		// exchanges second FACE
  		cudaStreamSynchronize(stream[2]);				// SYNCPOINT
  		
  				#ifdef ASYNC_TIMING
  				  mpi_start_sendrecv_2 = MPI_Wtime();
  				#endif

  		MPI_Sendrecv(RAND2, 24*tSliceEO, MPI_SHORT, g_nb_t_up, 1,	// SYNCPOINT
  		             RAND4, 24*tSliceEO, MPI_SHORT, g_nb_t_dn, 1,
  		             g_cart_grid, &stat[1]);
  		 //send norms            
  		MPI_Sendrecv(RAND2_norm, tSliceEO, MPI_FLOAT, g_nb_t_up, 1,	// SYNCPOINT
  		             RAND4_norm, tSliceEO, MPI_FLOAT, g_nb_t_dn, 1,
  		             g_cart_grid, &stat[1]);
  		             
  				#ifdef ASYNC_TIMING
  				  mpi_stop_sendrecv_2 = MPI_Wtime();
  				#endif
  		
	
  		// copies second FACE back to device
  		//MPI_Wait(&recv_req[1], &stat[1]);
  		cudaMemcpyAsync(spinin+6*(VolumeEO+tSliceEO), RAND4, tSliceEO*6*sizeof(short4), cudaMemcpyHostToDevice, stream[2]);
  		cudaMemcpyAsync(spinin_norm+(VolumeEO+tSliceEO), RAND4_norm, tSliceEO*sizeof(float), cudaMemcpyHostToDevice, stream[2]);
  				#ifdef ASYNC_TIMING
  				  cudaEventRecord(stop_H2D_4, stream[2]);
  				#endif
  		
  		
  		// applies second FACE
  		dev_Hopping_Matrix_half_ASYNC <<<gridsize2, blocksize, 0, stream[2]>>> ( gf,
  		                                                                    spinin, spinin_norm,
  		                                                                    spinout, spinout_norm,
  		                                                                    gfindex_site, gfindex_nextsite,
  		                                                                    nn_evenodd, ieo,
  		                                                                    0, tSlices*tSliceEO );
  				#ifdef ASYNC_TIMING
  				  cudaEventRecord(stop_EXT_2, stream[2]);
  				#endif
  		
  		
  				#ifdef ASYNC_TIMING
  				  cudaEventRecord(stop_ALL, 0);
  				#endif
  		
  		
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  
  
  #endif		// different optimized and non-optimized version
  
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  		
  		
  		cudaThreadSynchronize();		// test if needed	// for timing ...
  
  
  #ifdef USETEXTURE
    unbind_halfspinor_texture();
  #endif
  
}





#endif //HALF















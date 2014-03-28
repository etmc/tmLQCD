/**************************************************************************
 *
 * Copyright (C) 2014 Florian Burger
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
 
 
 

void to_device_mpi_d (dev_spinor_d * device, spinor * host, int size, int start, int end) {

  cudaMemcpy(device+12*start, host+12*start, size, cudaMemcpyHostToDevice);

}


void to_host_mpi_d (spinor * host, dev_spinor_d * device, int size, int start, int end) {

  cudaMemcpy(host+12*start, device+12*start, size, cudaMemcpyDeviceToHost);

}




//this gathers spinors in the space-time first indexing to a continuous block that can be exchanged via mpi
//starts at start and goes to start+size
//other versions for relativistic basis
__global__ void dev_gather_rand_d(dev_spinor_d * sin, dev_spinor_d * rand, int start, int size){
  int pos,pos2;
  pos2 = threadIdx.x + blockDim.x * blockIdx.x;
  pos = start  + pos2;
  
  if(pos < start + size){
  
  rand[12*pos2+0].x = sin[pos+0*DEVOFF].x; 
  rand[12*pos2+0].y = sin[pos+0*DEVOFF].y; 
  rand[12*pos2+1].x = sin[pos+1*DEVOFF].x;
  rand[12*pos2+1].y = sin[pos+1*DEVOFF].y;
  
  rand[12*pos2+2].x = sin[pos+2*DEVOFF].x; 
  rand[12*pos2+2].y = sin[pos+2*DEVOFF].y; 
  rand[12*pos2+3].x = sin[pos+3*DEVOFF].x;
  rand[12*pos2+3].y = sin[pos+3*DEVOFF].y; 
  
  rand[12*pos2+4].x = sin[pos+4*DEVOFF].x; 
  rand[12*pos2+4].y = sin[pos+4*DEVOFF].y; 
  rand[12*pos2+5].x = sin[pos+5*DEVOFF].x;
  rand[12*pos2+5].y = sin[pos+5*DEVOFF].y;  
  
  rand[12*pos2+6].x = sin[pos+6*DEVOFF].x;
  rand[12*pos2+6].y = sin[pos+6*DEVOFF].y;     
  rand[12*pos2+7].x = sin[pos+7*DEVOFF].x; 
  rand[12*pos2+7].y = sin[pos+7*DEVOFF].y;
  
  rand[12*pos2+8].x = sin[pos+8*DEVOFF].x;
  rand[12*pos2+8].y = sin[pos+8*DEVOFF].y;   
  rand[12*pos2+9].x = sin[pos+9*DEVOFF].x; 
  rand[12*pos2+9].y = sin[pos+9*DEVOFF].y;
  
  rand[12*pos2+10].x = sin[pos+10*DEVOFF].x;
  rand[12*pos2+10].y = sin[pos+10*DEVOFF].y; 
  rand[12*pos2+11].x = sin[pos+11*DEVOFF].x;
  rand[12*pos2+11].y = sin[pos+11*DEVOFF].y;  
  }
}


//this gathers spinors in the space-time first indexing to a continuous block that can be exchanged via mpi
//starts at start and goes to start+size
//RELATIVISTIC_BASIS!!
__global__ void dev_gather_rand_relup_d(dev_spinor_d * sin, dev_spinor_d * rand, int start, int size){
  int pos, pos2;
  pos2 = threadIdx.x + blockDim.x * blockIdx.x;
  pos = start  +  pos2;
  
  if(pos < start + size){
//only fetch upper spinor and store consecutively (NOTE the "6" instead of "12") 
  rand[6*pos2].x = sin[pos+0*DEVOFF].x; 
  rand[6*pos2].y = sin[pos+0*DEVOFF].y; 
  rand[6*pos2+1].x = sin[pos+1*DEVOFF].x;
  rand[6*pos2+1].y = sin[pos+1*DEVOFF].y;
  
  rand[6*pos2+2].x = sin[pos+2*DEVOFF].x; 
  rand[6*pos2+2].y = sin[pos+2*DEVOFF].y; 
  rand[6*pos2+3].x = sin[pos+3*DEVOFF].x;
  rand[6*pos2+3].y = sin[pos+3*DEVOFF].y; 
  
  rand[6*pos2+4].x = sin[pos+4*DEVOFF].x; 
  rand[6*pos2+4].y = sin[pos+4*DEVOFF].y; 
  rand[6*pos2+5].x = sin[pos+5*DEVOFF].x;
  rand[6*pos2+5].y = sin[pos+5*DEVOFF].y;  
  
  }
}





//this gathers spinors in the space-time first indexing to a continuous block that can be exchanged via mpi
//starts at start and goes to start+size
//RELATIVISTIC_BASIS!!
__global__ void dev_gather_rand_reldn_d(dev_spinor_d * sin, dev_spinor_d * rand, int start, int size){
  int pos,pos2;
  pos2 = threadIdx.x + blockDim.x * blockIdx.x;
  pos = start  +  pos2;
  
  if(pos < start + size){
//only fetch upper spinor and store consecutively (NOTE the "6" instead of "12") 
  rand[6*pos2+0].x = sin[pos+6*DEVOFF].x; 
  rand[6*pos2+0].y = sin[pos+6*DEVOFF].y; 
  rand[6*pos2+1].x = sin[pos+7*DEVOFF].x;
  rand[6*pos2+1].y = sin[pos+7*DEVOFF].y;
  
  rand[6*pos2+2].x = sin[pos+8*DEVOFF].x; 
  rand[6*pos2+2].y = sin[pos+8*DEVOFF].y; 
  rand[6*pos2+3].x = sin[pos+9*DEVOFF].x;
  rand[6*pos2+3].y = sin[pos+9*DEVOFF].y; 
  
  rand[6*pos2+4].x = sin[pos+10*DEVOFF].x; 
  rand[6*pos2+4].y = sin[pos+10*DEVOFF].y; 
  rand[6*pos2+5].x = sin[pos+11*DEVOFF].x;
  rand[6*pos2+5].y = sin[pos+11*DEVOFF].y;  
  
  }
}





//this spreads spinors to the space-time first indexing, rand is passed as a continuous block that is exchanged via mpi
//starts at start and goes to start+size
//other versions for relativistic basis
__global__ void dev_spread_rand_d(dev_spinor_d * sin, dev_spinor_d * rand, int start, int size){
  int pos, pos2;
  pos2 = threadIdx.x + blockDim.x * blockIdx.x;
  pos = start  + pos2; 
  
  if(pos < start + size){
  sin[pos+0*DEVOFF].x = rand[12*pos2+0].x; 
  sin[pos+0*DEVOFF].y = rand[12*pos2+0].y; 
  sin[pos+1*DEVOFF].x = rand[12*pos2+1].x;
  sin[pos+1*DEVOFF].y = rand[12*pos2+1].y;
  
  sin[pos+2*DEVOFF].x = rand[12*pos2+2].x; 
  sin[pos+2*DEVOFF].y = rand[12*pos2+2].y;
  sin[pos+3*DEVOFF].x = rand[12*pos2+3].x;
  sin[pos+3*DEVOFF].y = rand[12*pos2+3].y; 
  
  sin[pos+4*DEVOFF].x = rand[12*pos2+4].x; 
  sin[pos+4*DEVOFF].y = rand[12*pos2+4].y; 
  sin[pos+5*DEVOFF].x = rand[12*pos2+5].x;
  sin[pos+5*DEVOFF].y = rand[12*pos2+5].y;  
  
  sin[pos+6*DEVOFF].x = rand[12*pos2+6].x; 
  sin[pos+6*DEVOFF].y = rand[12*pos2+6].y; 
  sin[pos+7*DEVOFF].x = rand[12*pos2+7].x;
  sin[pos+7*DEVOFF].y = rand[12*pos2+7].y;   
  
  sin[pos+8*DEVOFF].x = rand[12*pos2+8].x; 
  sin[pos+8*DEVOFF].y = rand[12*pos2+8].y; 
  sin[pos+9*DEVOFF].x = rand[12*pos2+9].x;
  sin[pos+9*DEVOFF].y = rand[12*pos2+9].y; 
  
  sin[pos+10*DEVOFF].x = rand[12*pos2+10].x; 
  sin[pos+10*DEVOFF].y = rand[12*pos2+10].y; 
  sin[pos+11*DEVOFF].x = rand[12*pos2+11].x;
  sin[pos+11*DEVOFF].y = rand[12*pos2+11].y; 
  
  }
}



//this spreads spinors to the space-time first indexing, rand is passed as a continuous block that is exchanged via mpi
//starts at start and goes to start+size
//RELATIVISTIC_BASIS !!!
__global__ void dev_spread_rand_relup_d(dev_spinor_d * sin, dev_spinor_d * rand, int start, int size){
  int pos, pos2;
  pos2 = threadIdx.x + blockDim.x * blockIdx.x;
  pos = start  + pos2; 
  
  if(pos < start + size){
//only spread upper spinor and store consecutively (NOTE the "6" instead of "12") 
  sin[pos+0*DEVOFF].x = rand[6*pos2+0].x; 
  sin[pos+0*DEVOFF].y = rand[6*pos2+0].y; 
  sin[pos+1*DEVOFF].x = rand[6*pos2+1].x;
  sin[pos+1*DEVOFF].y = rand[6*pos2+1].y;
  
  sin[pos+2*DEVOFF].x = rand[6*pos2+2].x; 
  sin[pos+2*DEVOFF].y = rand[6*pos2+2].y;
  sin[pos+3*DEVOFF].x = rand[6*pos2+3].x;
  sin[pos+3*DEVOFF].y = rand[6*pos2+3].y; 
  
  sin[pos+4*DEVOFF].x = rand[6*pos2+4].x; 
  sin[pos+4*DEVOFF].y = rand[6*pos2+4].y; 
  sin[pos+5*DEVOFF].x = rand[6*pos2+5].x;
  sin[pos+5*DEVOFF].y = rand[6*pos2+5].y;  
  
  }
}



//this spreads spinors to the space-time first indexing, rand is passed as a continuous block that is exchanged via mpi
//starts at start and goes to start+size
//RELATIVISTIC_BASIS !!!
__global__ void dev_spread_rand_reldn_d(dev_spinor_d * sin, dev_spinor_d * rand, int start, int size){
  int pos,pos2;
  pos2 = threadIdx.x + blockDim.x * blockIdx.x;
  pos = start  +  pos2;
  
  if(pos < start + size){
//only spread lower spinor and store consecutively (NOTE the "6" instead of "12")  
  
  sin[pos+6*DEVOFF].x = rand[6*pos2+0].x; 
  sin[pos+6*DEVOFF].y = rand[6*pos2+0].y; 
  sin[pos+7*DEVOFF].x = rand[6*pos2+1].x;
  sin[pos+7*DEVOFF].y = rand[6*pos2+1].y;   
  
  sin[pos+8*DEVOFF].x = rand[6*pos2+2].x; 
  sin[pos+8*DEVOFF].y = rand[6*pos2+2].y; 
  sin[pos+9*DEVOFF].x = rand[6*pos2+3].x;
  sin[pos+9*DEVOFF].y = rand[6*pos2+3].y; 
  
  sin[pos+10*DEVOFF].x = rand[6*pos2+4].x; 
  sin[pos+10*DEVOFF].y = rand[6*pos2+4].y; 
  sin[pos+11*DEVOFF].x = rand[6*pos2+5].x;
  sin[pos+11*DEVOFF].y = rand[6*pos2+5].y; 
  
  }
}



// copies the boundary t-slices t=0 and t=T-1 to host		// will be used in dev_Qtm_pm_ndpsi_mpi(), not ASYNC
//	exchanges						// provides a wrapped version of Carsten's xchange_field()
//		copies RAND back to device			//	and not asynchronous version of ASYNC.cuh

void xchange_field_wrapper_d (dev_spinor_d * dev_spin, int ieo) {
  

    
    int tSliceEO = LX*LY*LZ/2;
    int VolumeEO = VOLUME/2;
    
    //this is the same partitioning as for dev_mul_one_pm...
    int gridsize;
    int blockdim = BLOCK2D;
    if( tSliceEO % blockdim == 0){
      gridsize = (int) tSliceEO/blockdim;
    }
    else{
      gridsize = (int) tSliceEO/blockdim + 1;
    }
    int griddim = gridsize;
    
    
    #ifdef RELATIVISTIC_BASIS
      //this goes backwards, so for the receiver this is forward
      dev_gather_rand_relup_d<<<griddim, blockdim >>>(dev_spin,RAND_BW_D,0,tSliceEO);
      cudaMemcpy(R1_D, RAND_BW_D, tSliceEO*6*sizeof(double2), cudaMemcpyDeviceToHost);
    
      //this goes forward, so for the receiver this is backward
      dev_gather_rand_reldn_d<<<griddim, blockdim >>>(dev_spin,RAND_FW_D,(VolumeEO-tSliceEO),tSliceEO);
      cudaMemcpy(R2_D, RAND_FW_D, tSliceEO*6*sizeof(double2), cudaMemcpyDeviceToHost);    
    #else
      //this goes backwards
      dev_gather_rand_d<<<griddim, blockdim >>>(dev_spin,RAND_BW_D,0,tSliceEO);
      cudaMemcpy(R1_D, RAND_BW_D, tSliceEO*12*sizeof(double2), cudaMemcpyDeviceToHost);
    
      //this goes forward
      dev_gather_rand_d<<<griddim, blockdim >>>(dev_spin,RAND_FW_D,(VolumeEO-tSliceEO),tSliceEO);
      cudaMemcpy(R2_D, RAND_FW_D, tSliceEO*12*sizeof(double2), cudaMemcpyDeviceToHost);
    #endif
    
    //we only need to exchange half of the spinors in relativistic basis (upper or lower part)
    int ndouble_per_spin;
    #ifdef RELATIVISTIC_BASIS
      ndouble_per_spin = 12;
    #else
      ndouble_per_spin = 24;
    #endif
    
    MPI_Sendrecv(R1_D, ndouble_per_spin*tSliceEO, MPI_DOUBLE, g_nb_t_dn, 0,
                 R3_D, ndouble_per_spin*tSliceEO, MPI_DOUBLE, g_nb_t_up, 0,
                 g_cart_grid, &stat[0]);
    MPI_Sendrecv(R2_D, ndouble_per_spin*tSliceEO, MPI_DOUBLE, g_nb_t_up, 1,
                 R4_D, ndouble_per_spin*tSliceEO, MPI_DOUBLE, g_nb_t_dn, 1,
                 g_cart_grid, &stat[1]);

    #ifdef RELATIVISTIC_BASIS
      cudaMemcpy(RAND_BW_D, R3_D, tSliceEO*6*sizeof(double2), cudaMemcpyHostToDevice);
      dev_spread_rand_relup_d<<<griddim, blockdim >>>(dev_spin,RAND_BW_D,VolumeEO,tSliceEO);
    
      cudaMemcpy(RAND_FW_D, R4_D, tSliceEO*6*sizeof(double2), cudaMemcpyHostToDevice);
      dev_spread_rand_reldn_d<<<griddim, blockdim >>>(dev_spin,RAND_FW_D,(VolumeEO+tSliceEO),tSliceEO);   
    #else 
      cudaMemcpy(RAND_BW_D, R3_D, tSliceEO*12*sizeof(double2), cudaMemcpyHostToDevice);
      dev_spread_rand_d<<<griddim, blockdim >>>(dev_spin,RAND_BW_D,VolumeEO,tSliceEO);
    
      cudaMemcpy(RAND_FW_D, R4_D, tSliceEO*12*sizeof(double2), cudaMemcpyHostToDevice);
      dev_spread_rand_d<<<griddim, blockdim >>>(dev_spin,RAND_FW_D,(VolumeEO+tSliceEO),tSliceEO);
    #endif

  
}




















/*
This is the wrapper function for the device hopping matrix for MPI with support of
CUDA streams in order to parallelize bulk calculation and boundary exchange
*/
//! relativistic basis needs to be put into dev_Hopping_Matrix_d before we can uncomment it here
void HOPPING_ASYNC_D (dev_su3_2v_d * gf, 
                    dev_spinor_d * spinin, dev_spinor_d * spinout,
                    int * gfindex_site, int * gfindex_nextsite, int * nn_evenodd,
                    int ieo,
                    int gridsize, int blocksize) {
  
  
  // for even/odd
  int tSliceEO = LX*LY*LZ/2;
  int VolumeEO = VOLUME/2;
  
//!FIXME all block and grid inits should be done elsewhere  
  // gridsizes
  int gridsize1;
  int gridsize2;
  

  
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

  
  #if ASYNC == 0		// primitive version
    		

  
		// applies to the parts which don't need communication
  		dev_Hopping_Matrix_d <<<gridsize1, blocksize>>> ( gf,
        	                                                      spinin, spinout,
        	                                                      gfindex_site, gfindex_nextsite, nn_evenodd,
        	                                                      ieo,
        	                                                      //2*tSliceEO, VolumeEO-4*tSliceEO );
        	                                                       tSliceEO, VolumeEO-2*tSliceEO );
	                                                       


  		
  		
  		// exchanges the boundaries
  		xchange_field_wrapper_d(spinin, ieo);			// to be further optimized !!
  		
  		// applies the hopping matrix to remaining parts
  		dev_Hopping_Matrix_d <<<gridsize2, blocksize>>> ( gf,
  		                                                      spinin, spinout,
  		                                                      gfindex_site, gfindex_nextsite, nn_evenodd,
  		                                                      ieo,
  		                                                      //0, 2*tSliceEO );
  		                                                      0, tSliceEO );
  		
  		dev_Hopping_Matrix_d <<<gridsize2, blocksize>>> ( gf,
  		                                                      spinin, spinout,
  		                                                      gfindex_site, gfindex_nextsite, nn_evenodd,
  		                                                      ieo,
  		                                                      //VolumeEO-2*tSliceEO, 2*tSliceEO );
  		                                                      VolumeEO-tSliceEO, tSliceEO );	
		
		
		
		
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  
  
  
  #elif ASYNC == 1		// optimized version
    //this is the same partitioning as for dev_mul_one_pm...
    //we need this here, as tSliceEO is the partition size FIXME
    int gridsize3;
    int blocksize3 = BLOCK2D;
    if( tSliceEO % blocksize3 == 0){
      gridsize3 = (int) tSliceEO/blocksize3;
    }
    else{
      gridsize3 = (int) tSliceEO/blocksize3 + 1;
    } 		
		
                //set the amount of data we need to transfer with mpi
      		#ifdef RELATIVISTIC_BASIS
  		  int dbperspin = 12;
  		#else
  		  int dbperspin = 24;
  		#endif	
  		
				#ifdef ASYNC_TIMING
  				  cudaEventRecord(start_ALL, 0);
  				  mpi_start_ALL = MPI_Wtime();
  				#endif
        	
        	
        	// copies first FACE to host
				  
        	#ifdef RELATIVISTIC_BASIS
        	  dev_gather_rand_relup_d<<<gridsize3, blocksize3 ,0,stream[1] >>>(spinin,RAND_BW_D,0,tSliceEO);
  		  cudaMemcpyAsync(RAND1_D, RAND_BW_D, tSliceEO*6*sizeof(double2), cudaMemcpyDeviceToHost, stream[1]);
  		#else
        	  dev_gather_rand_d<<<gridsize3, blocksize3,0,stream[1] >>>(spinin,RAND_BW_D,0,tSliceEO);
  		  cudaMemcpyAsync(RAND1_D, RAND_BW_D, tSliceEO*12*sizeof(double2), cudaMemcpyDeviceToHost, stream[1]);  		
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
  		  dev_gather_rand_reldn_d<<<gridsize3, blocksize3, 0, stream[2] >>>(spinin,RAND_FW_D,(VolumeEO-tSliceEO),tSliceEO);
  		  cudaMemcpyAsync(RAND2_D, RAND_FW_D, tSliceEO*6*sizeof(double2), cudaMemcpyDeviceToHost, stream[2]);
  		#else
                  dev_gather_rand_d<<<gridsize3, blocksize3, 0, stream[2] >>>(spinin,RAND_FW_D,(VolumeEO-tSliceEO),tSliceEO);
  		  cudaMemcpyAsync(RAND2_D, RAND_FW_D, tSliceEO*12*sizeof(double2), cudaMemcpyDeviceToHost, stream[2]);  		
		  //printf("g_proc_id = %d: R2  %e \n", g_proc_id, RAND2[11].x);
		#endif
  				#ifdef ASYNC_TIMING
  				  cudaEventRecord(stop_D2H_2, stream[2]);
  				#endif	
        	

        	
        	
//INTERNAL        	
        	// starts INTERNAL kernel
  		dev_Hopping_Matrix_d <<<gridsize1, blocksize, 0, stream[0]>>> ( gf,
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
  		
  		
  		MPI_Sendrecv(RAND1_D, dbperspin*tSliceEO, MPI_DOUBLE, g_nb_t_dn, 0,	// SYNCPOINT
  		             RAND3_D, dbperspin*tSliceEO, MPI_DOUBLE, g_nb_t_up, 0,
  		             g_cart_grid, &stat[0]);
  		
  				#ifdef ASYNC_TIMING
  				  mpi_stop_sendrecv_1 = MPI_Wtime();
  				#endif
  		
  		
  		// copies first FACE back to device				
  		
  		
  		// order may switched
  		//MPI_Wait(&recv_req[0], &stat[0]);									
  		// synchronous
                
		#ifdef RELATIVISTIC_BASIS 
  		  cudaMemcpyAsync(RAND_BW_D, RAND3_D, tSliceEO*6*sizeof(double2), cudaMemcpyHostToDevice, stream[1]);
  		  dev_spread_rand_relup_d<<<gridsize3, blocksize3, 0, stream[1] >>>(spinin,RAND_BW_D,VolumeEO,tSliceEO);
                #else
                  //printf("g_proc_id = %d:  R3 %e \n", g_proc_id, RAND3[11].x);
                  cudaMemcpyAsync(RAND_BW_D, RAND3_D, tSliceEO*12*sizeof(double2), cudaMemcpyHostToDevice, stream[1]);
                  dev_spread_rand_d<<<gridsize3, blocksize3, 0, stream[1] >>>(spinin,RAND_BW_D,VolumeEO,tSliceEO);
                #endif

  				#ifdef ASYNC_TIMING
  				  cudaEventRecord(stop_H2D_3, stream[1]);
  				#endif
  		
  		
  		// applies first FACE
  		dev_Hopping_Matrix_d <<<gridsize2, blocksize, 0, stream[1]>>> ( gf,
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

  		
  		MPI_Sendrecv(RAND2_D, dbperspin*tSliceEO, MPI_DOUBLE, g_nb_t_up, 1,	// SYNCPOINT
  		             RAND4_D, dbperspin*tSliceEO, MPI_DOUBLE, g_nb_t_dn, 1,
  		             g_cart_grid, &stat[1]);
  		
  				#ifdef ASYNC_TIMING
  				  mpi_stop_sendrecv_2 = MPI_Wtime();
  				#endif
  		
  		// copies second FACE back to device
  		//MPI_Wait(&recv_req[1], &stat[1]); 
				  
  		#ifdef RELATIVISTIC_BASIS 
  		  cudaMemcpyAsync(RAND_FW_D, RAND4_D, tSliceEO*6*sizeof(double2), cudaMemcpyHostToDevice, stream[2]);
  		  dev_spread_rand_reldn_d<<<gridsize3, blocksize3, 0, stream[2] >>>(spinin,RAND_FW_D,VolumeEO+tSliceEO,tSliceEO);
  		#else
  		  //printf("g_proc_id = %d: R4  %e \n", g_proc_id, RAND4[11].x);
   		  cudaMemcpyAsync(RAND_FW_D, RAND4_D, tSliceEO*12*sizeof(double2), cudaMemcpyHostToDevice, stream[2]);
  		  dev_spread_rand_d<<<gridsize3, blocksize3, 0, stream[2] >>>(spinin,RAND_FW_D,VolumeEO+tSliceEO,tSliceEO); 		
  		#endif
		
  				#ifdef ASYNC_TIMING
  				  cudaEventRecord(stop_H2D_4, stream[2]);
  				#endif
  		
  		// applies second FACE
  		dev_Hopping_Matrix_d <<<gridsize2, blocksize, 0, stream[2]>>> ( gf,
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
  

  
}



//! relativistic basis needs to be put into dev_Hopping_Matrix_d before we can uncomment it here
void HOPPING_ASYNC_UPDN_D (dev_su3_2v_d * gf, 
                    dev_spinor_d * spinin_up, dev_spinor_d * spinin_dn, 
		    dev_spinor_d * spinout_up, dev_spinor_d * spinout_dn,
                    int * gfindex_site, int * gfindex_nextsite, int * nn_evenodd,
                    int ieo,
                    int gridsize, int blocksize) {
  
  
  // for even/odd
  int tSliceEO = LX*LY*LZ/2;
  int VolumeEO = VOLUME/2;
  
//!FIXME all block and grid inits should be done elsewhere  
  // gridsizes
  int gridsize1;
  int gridsize2;
  

  
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

  
  #if ASYNC == 0		// primitive version
    		

  
		// applies to the parts which don't need communication
  		dev_Hopping_Matrix_updn_d <<<gridsize1, blocksize>>> ( gf,
        	                                                      spinin_up, spinin_dn, spinout_up, spinout_dn,
        	                                                      gfindex_site, gfindex_nextsite, nn_evenodd,
        	                                                      ieo,
        	                                                      //2*tSliceEO, VolumeEO-4*tSliceEO );
        	                                                       tSliceEO, VolumeEO-2*tSliceEO );
	                                                       


  		
  		
  		// exchanges the boundaries
  		xchange_field_wrapper_d(spinin_up, ieo);
  		xchange_field_wrapper_d(spinin_dn, ieo);
  		
  		// applies the hopping matrix to remaining parts
  		dev_Hopping_Matrix_updn_d <<<gridsize2, blocksize>>> ( gf,
  		                                                      spinin_up, spinin_dn, spinout_up, spinout_dn,
  		                                                      gfindex_site, gfindex_nextsite, nn_evenodd,
  		                                                      ieo,
  		                                                      //0, 2*tSliceEO );
  		                                                      0, tSliceEO );
  		
  		dev_Hopping_Matrix_updn_d <<<gridsize2, blocksize>>> ( gf,
  		                                                      spinin_up, spinin_dn, spinout_up, spinout_dn,
  		                                                      gfindex_site, gfindex_nextsite, nn_evenodd,
  		                                                      ieo,
  		                                                      //VolumeEO-2*tSliceEO, 2*tSliceEO );
  		                                                      VolumeEO-tSliceEO, tSliceEO );	
		
		
		
		
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  
  
  
  #elif ASYNC == 1		// optimized version
		
    //this is the same partitioning as for dev_mul_one_pm...
    int gridsize3;
    int blocksize3 = BLOCK2;
    if( tSliceEO % blocksize3 == 0){
      gridsize3 = (int) tSliceEO/blocksize3;
    }
    else{
      gridsize3 = (int) tSliceEO/blocksize3 + 1;
    }  
    
                //set the amount of data we need to transfer with mpi
      		#ifdef RELATIVISTIC_BASIS
  		  int dbperspin = 12;
  		#else
  		  int dbperspin = 24;
  		#endif	
  		
				#ifdef ASYNC_TIMING
  				  cudaEventRecord(start_ALL, 0);
  				  mpi_start_ALL = MPI_Wtime();
  				#endif
        	
        	
        	// copies first FACE to host
        	#ifdef RELATIVISTIC_BASIS
        	  dev_gather_rand_relup_d<<<gridsize3, blocksize3 ,0,stream[1] >>>(spinin_up,RAND_BW_UP_D,0,tSliceEO);
  		  cudaMemcpyAsync(RAND1_UP_D, RAND_BW_UP_D, tSliceEO*6*sizeof(double2), cudaMemcpyDeviceToHost, stream[1]);
        	  dev_gather_rand_relup_d<<<gridsize3, blocksize3 ,0,stream[1] >>>(spinin_dn,RAND_BW_DN_D,0,tSliceEO);
  		  cudaMemcpyAsync(RAND1_DN_D, RAND_BW_DN_D, tSliceEO*6*sizeof(double2), cudaMemcpyDeviceToHost, stream[1]);		  
  		#else
        	  dev_gather_rand_d<<<gridsize3, blocksize3,0,stream[1] >>>(spinin_up,RAND_BW_UP_D,0,tSliceEO);
  		  cudaMemcpyAsync(RAND1_UP_D, RAND_BW_UP_D, tSliceEO*12*sizeof(double2), cudaMemcpyDeviceToHost, stream[1]);  		
        	  dev_gather_rand_d<<<gridsize3, blocksize3,0,stream[1] >>>(spinin_dn,RAND_BW_DN_D,0,tSliceEO);
  		  cudaMemcpyAsync(RAND1_DN_D, RAND_BW_DN_D, tSliceEO*12*sizeof(double2), cudaMemcpyDeviceToHost, stream[1]); 
		  //printf("g_proc_id = %d: R1  %e \n", g_proc_id, RAND1[11].x);
		#endif
  				#ifdef ASYNC_TIMING
  				  cudaEventRecord(stop_D2H_1, stream[1]);
  				#endif



                // copies second FACE to host
  		#ifdef RELATIVISTIC_BASIS
  		  dev_gather_rand_reldn_d<<<gridsize3, blocksize3, 0, stream[2] >>>(spinin_up,RAND_FW_UP_D,(VolumeEO-tSliceEO),tSliceEO);
  		  cudaMemcpyAsync(RAND2_UP_D, RAND_FW_UP_D, tSliceEO*6*sizeof(double2), cudaMemcpyDeviceToHost, stream[2]);
  		  dev_gather_rand_reldn_d<<<gridsize3, blocksize3, 0, stream[2] >>>(spinin_dn,RAND_FW_DN_D,(VolumeEO-tSliceEO),tSliceEO);
  		  cudaMemcpyAsync(RAND2_DN_D, RAND_FW_DN_D, tSliceEO*6*sizeof(double2), cudaMemcpyDeviceToHost, stream[2]);		  
  		#else
                  dev_gather_rand_d<<<gridsize3, blocksize3, 0, stream[2] >>>(spinin_up,RAND_FW_UP_D,(VolumeEO-tSliceEO),tSliceEO);
  		  cudaMemcpyAsync(RAND2_UP_D, RAND_FW_UP_D, tSliceEO*12*sizeof(double2), cudaMemcpyDeviceToHost, stream[2]); 
                  dev_gather_rand_d<<<gridsize3, blocksize3, 0, stream[2] >>>(spinin_dn,RAND_FW_DN_D,(VolumeEO-tSliceEO),tSliceEO);
  		  cudaMemcpyAsync(RAND2_DN_D, RAND_FW_DN_D, tSliceEO*12*sizeof(double2), cudaMemcpyDeviceToHost, stream[2]); 		  
		  //printf("g_proc_id = %d: R2  %e \n", g_proc_id, RAND2[11].x);
		#endif
  				#ifdef ASYNC_TIMING
  				  cudaEventRecord(stop_D2H_2, stream[2]);
  				#endif	
        	

        	
        	
//INTERNAL        	
        	// starts INTERNAL kernel
  		dev_Hopping_Matrix_updn_d <<<gridsize1, blocksize, 0, stream[0]>>> ( gf,
        	                                                                    spinin_up, spinin_dn, spinout_up, spinout_dn,
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
  		
  		
  		MPI_Sendrecv(RAND1_UP_D, dbperspin*tSliceEO, MPI_DOUBLE, g_nb_t_dn, 0,	// SYNCPOINT
  		             RAND3_UP_D, dbperspin*tSliceEO, MPI_DOUBLE, g_nb_t_up, 0,
  		             g_cart_grid, &stat[0]);
  		MPI_Sendrecv(RAND1_DN_D, dbperspin*tSliceEO, MPI_DOUBLE, g_nb_t_dn, 0,	// SYNCPOINT
  		             RAND3_DN_D, dbperspin*tSliceEO, MPI_DOUBLE, g_nb_t_up, 0,
  		             g_cart_grid, &stat[0]);
		
  				#ifdef ASYNC_TIMING
  				  mpi_stop_sendrecv_1 = MPI_Wtime();
  				#endif
  		
  		
  		// copies first FACE back to device                
		#ifdef RELATIVISTIC_BASIS 
  		  cudaMemcpyAsync(RAND_BW_UP_D, RAND3_UP_D, tSliceEO*6*sizeof(double2), cudaMemcpyHostToDevice, stream[1]);
  		  dev_spread_rand_relup_d<<<gridsize3, blocksize3, 0, stream[1] >>>(spinin_up,RAND_BW_UP_D,VolumeEO,tSliceEO);
  		  cudaMemcpyAsync(RAND_BW_DN_D, RAND3_DN_D, tSliceEO*6*sizeof(double2), cudaMemcpyHostToDevice, stream[1]);
  		  dev_spread_rand_relup_d<<<gridsize3, blocksize3, 0, stream[1] >>>(spinin_dn,RAND_BW_DN_D,VolumeEO,tSliceEO);
		#else
                  //printf("g_proc_id = %d:  R3 %e \n", g_proc_id, RAND3[11].x);
                  cudaMemcpyAsync(RAND_BW_UP_D, RAND3_UP_D, tSliceEO*12*sizeof(double2), cudaMemcpyHostToDevice, stream[1]);
                  dev_spread_rand_d<<<gridsize3, blocksize3, 0, stream[1] >>>(spinin_up,RAND_BW_UP_D,VolumeEO,tSliceEO);
                  cudaMemcpyAsync(RAND_BW_DN_D, RAND3_DN_D, tSliceEO*12*sizeof(double2), cudaMemcpyHostToDevice, stream[1]);
                  dev_spread_rand_d<<<gridsize3, blocksize3, 0, stream[1] >>>(spinin_dn,RAND_BW_DN_D,VolumeEO,tSliceEO);		  
                #endif

  				#ifdef ASYNC_TIMING
  				  cudaEventRecord(stop_H2D_3, stream[1]);
  				#endif
  		
  		
  		// applies first FACE
  		dev_Hopping_Matrix_updn_d <<<gridsize2, blocksize, 0, stream[1]>>> ( gf,
  		                                                                      spinin_up, spinin_dn, spinout_up, spinout_dn,
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

  		
  		MPI_Sendrecv(RAND2_UP_D, dbperspin*tSliceEO, MPI_DOUBLE, g_nb_t_up, 1,	// SYNCPOINT
  		             RAND4_UP_D, dbperspin*tSliceEO, MPI_DOUBLE, g_nb_t_dn, 1,
  		             g_cart_grid, &stat[1]);
  		MPI_Sendrecv(RAND2_DN_D, dbperspin*tSliceEO, MPI_DOUBLE, g_nb_t_up, 1,	// SYNCPOINT
  		             RAND4_DN_D, dbperspin*tSliceEO, MPI_DOUBLE, g_nb_t_dn, 1,
  		             g_cart_grid, &stat[1]);
		
  				#ifdef ASYNC_TIMING
  				  mpi_stop_sendrecv_2 = MPI_Wtime();
  				#endif
 				  
  		// copies second FACE back to device
  		//MPI_Wait(&recv_req[1], &stat[1]); 
				  
  		#ifdef RELATIVISTIC_BASIS 
  		  cudaMemcpyAsync(RAND_FW_UP_D, RAND4_UP_D, tSliceEO*6*sizeof(double2), cudaMemcpyHostToDevice, stream[2]);
  		  dev_spread_rand_reldn_d<<<gridsize3, blocksize3, 0, stream[2] >>>(spinin_up,RAND_FW_UP_D,VolumeEO+tSliceEO,tSliceEO);
  		  cudaMemcpyAsync(RAND_FW_DN_D, RAND4_DN_D, tSliceEO*6*sizeof(double2), cudaMemcpyHostToDevice, stream[2]);
  		  dev_spread_rand_reldn_d<<<gridsize3, blocksize3, 0, stream[2] >>>(spinin_dn,RAND_FW_DN_D,VolumeEO+tSliceEO,tSliceEO);		  
  		#else
  		  //printf("g_proc_id = %d: R4  %e \n", g_proc_id, RAND4[11].x);
   		  cudaMemcpyAsync(RAND_FW_UP_D, RAND4_UP_D, tSliceEO*12*sizeof(double2), cudaMemcpyHostToDevice, stream[2]);
  		  dev_spread_rand_d<<<gridsize3, blocksize3, 0, stream[2] >>>(spinin_up,RAND_FW_UP_D,VolumeEO+tSliceEO,tSliceEO);
   		  cudaMemcpyAsync(RAND_FW_DN_D, RAND4_DN_D, tSliceEO*12*sizeof(double2), cudaMemcpyHostToDevice, stream[2]);
  		  dev_spread_rand_d<<<gridsize3, blocksize3, 0, stream[2] >>>(spinin_dn,RAND_FW_DN_D,VolumeEO+tSliceEO,tSliceEO);		  
  		#endif
		
  				#ifdef ASYNC_TIMING
  				  cudaEventRecord(stop_H2D_4, stream[2]);
  				#endif
  		
  		// applies second FACE
  		dev_Hopping_Matrix_updn_d <<<gridsize2, blocksize, 0, stream[2]>>> ( gf,
  		                                                                    spinin_up, spinin_dn, spinout_up, spinout_dn,
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
  

  
}











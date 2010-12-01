


// convert spinor to double 
void convert2double_spin (dev_spinor* spin, spinor* h2d) {

  int i, Vol;
  
  //#ifndef MPI
    if (even_odd_flag) {
      Vol = VOLUME/2;
    }
    else {
      Vol = VOLUME;
    }
  //#else
  //  Vol = (VOLUME+RAND)/2;
  //#endif
  
  
  for (i = 0; i < Vol; i++) {
  
        h2d[i].s0.c0.re = (double) spin[6*i+0].x;
        h2d[i].s0.c0.im = (double) spin[6*i+0].y;
        h2d[i].s0.c1.re = (double) spin[6*i+0].z;
        h2d[i].s0.c1.im = (double) spin[6*i+0].w;
        
        h2d[i].s0.c2.re = (double) spin[6*i+1].x;
        h2d[i].s0.c2.im = (double) spin[6*i+1].y;
        h2d[i].s1.c0.re = (double) spin[6*i+1].z;
        h2d[i].s1.c0.im = (double) spin[6*i+1].w;   
        
        h2d[i].s1.c1.re = (double) spin[6*i+2].x;
        h2d[i].s1.c1.im = (double) spin[6*i+2].y;
        h2d[i].s1.c2.re = (double) spin[6*i+2].z;
        h2d[i].s1.c2.im = (double) spin[6*i+2].w;  
        
        h2d[i].s2.c0.re = (double) spin[6*i+3].x;
        h2d[i].s2.c0.im = (double) spin[6*i+3].y;
        h2d[i].s2.c1.re = (double) spin[6*i+3].z;
        h2d[i].s2.c1.im = (double) spin[6*i+3].w;  
        
        h2d[i].s2.c2.re = (double) spin[6*i+4].x;
        h2d[i].s2.c2.im = (double) spin[6*i+4].y;
        h2d[i].s3.c0.re = (double) spin[6*i+4].z;
        h2d[i].s3.c0.im = (double) spin[6*i+4].w; 
        
        h2d[i].s3.c1.re = (double) spin[6*i+5].x;
        h2d[i].s3.c1.im = (double) spin[6*i+5].y;
        h2d[i].s3.c2.re = (double) spin[6*i+5].z;
        h2d[i].s3.c2.im = (double) spin[6*i+5].w; 
        
  }
}





// convert spinor to REAL4 (float4, double4) 
void convert2REAL4_spin(spinor* spin, dev_spinor* h2d){

  int i, Vol;
  
  //#ifndef MPI
    if (even_odd_flag) {
      Vol = VOLUME/2;
    }
    else {
      Vol = VOLUME;
    }
  //#else
  //  Vol = (VOLUME+RAND)/2;
  //#endif
  
  for (i = 0; i < Vol; i++) {
    
        h2d[6*i+0].x = (REAL) spin[i].s0.c0.re;
        h2d[6*i+0].y = (REAL) spin[i].s0.c0.im;
        h2d[6*i+0].z = (REAL) spin[i].s0.c1.re;
        h2d[6*i+0].w = (REAL) spin[i].s0.c1.im;
        
        h2d[6*i+1].x = (REAL) spin[i].s0.c2.re;
        h2d[6*i+1].y = (REAL) spin[i].s0.c2.im;
        h2d[6*i+1].z = (REAL) spin[i].s1.c0.re;
        h2d[6*i+1].w = (REAL) spin[i].s1.c0.im;
        
        h2d[6*i+2].x = (REAL) spin[i].s1.c1.re;
        h2d[6*i+2].y = (REAL) spin[i].s1.c1.im;
        h2d[6*i+2].z = (REAL) spin[i].s1.c2.re;
        h2d[6*i+2].w = (REAL) spin[i].s1.c2.im;
        
        h2d[6*i+3].x = (REAL) spin[i].s2.c0.re;
        h2d[6*i+3].y = (REAL) spin[i].s2.c0.im;
        h2d[6*i+3].z = (REAL) spin[i].s2.c1.re;
        h2d[6*i+3].w = (REAL) spin[i].s2.c1.im;
        
        h2d[6*i+4].x = (REAL) spin[i].s2.c2.re;
        h2d[6*i+4].y = (REAL) spin[i].s2.c2.im;
        h2d[6*i+4].z = (REAL) spin[i].s3.c0.re;
        h2d[6*i+4].w = (REAL) spin[i].s3.c0.im;
        
        h2d[6*i+5].x = (REAL) spin[i].s3.c1.re;
        h2d[6*i+5].y = (REAL) spin[i].s3.c1.im;
        h2d[6*i+5].z = (REAL) spin[i].s3.c2.re;
        h2d[6*i+5].w = (REAL) spin[i].s3.c2.im;
    
  }
}






// host/device interaction

// remark: the host spinors are double precision and therefore need twice the memory !!
//		dev_spinor * device:    dev_spinsize
//		spinor * host:        2*dev_spinsize
//		dev_spinor * auxiliary: dev_spinsize
//         the parameter "size" specifies the memory needed for the spinor n the device !!
//         

void to_device (dev_spinor * device, spinor * host, dev_spinor * auxiliary, int size) {

  convert2REAL4_spin(host, auxiliary);						// auxiliary = (float) host
  cudaMemcpy(device, auxiliary, size, cudaMemcpyHostToDevice);			// device = auxiliary  (on device)

}


void to_host (spinor * host, dev_spinor * device, dev_spinor * auxiliary, int size) {

  cudaMemcpy(auxiliary, device, size, cudaMemcpyDeviceToHost);			// auxiliary = device  (on device)
  convert2double_spin(auxiliary, host);						// host = (double) auxiliary

}




/////////
// MPI //
/////////


#ifdef MPI

// convert spinor to double

void convert2double_spin_mpi (dev_spinor * spin, spinor * h2d, int start, int end) {

  int i;
  
  for (i = start; i < end; i++) {
  
        h2d[i].s0.c0.re = (double) spin[6*i+0].x;
        h2d[i].s0.c0.im = (double) spin[6*i+0].y;
        h2d[i].s0.c1.re = (double) spin[6*i+0].z;
        h2d[i].s0.c1.im = (double) spin[6*i+0].w;
        
        h2d[i].s0.c2.re = (double) spin[6*i+1].x;
        h2d[i].s0.c2.im = (double) spin[6*i+1].y;
        h2d[i].s1.c0.re = (double) spin[6*i+1].z;
        h2d[i].s1.c0.im = (double) spin[6*i+1].w;   
        
        h2d[i].s1.c1.re = (double) spin[6*i+2].x;
        h2d[i].s1.c1.im = (double) spin[6*i+2].y;
        h2d[i].s1.c2.re = (double) spin[6*i+2].z;
        h2d[i].s1.c2.im = (double) spin[6*i+2].w;  
        
        h2d[i].s2.c0.re = (double) spin[6*i+3].x;
        h2d[i].s2.c0.im = (double) spin[6*i+3].y;
        h2d[i].s2.c1.re = (double) spin[6*i+3].z;
        h2d[i].s2.c1.im = (double) spin[6*i+3].w;  
        
        h2d[i].s2.c2.re = (double) spin[6*i+4].x;
        h2d[i].s2.c2.im = (double) spin[6*i+4].y;
        h2d[i].s3.c0.re = (double) spin[6*i+4].z;
        h2d[i].s3.c0.im = (double) spin[6*i+4].w; 
        
        h2d[i].s3.c1.re = (double) spin[6*i+5].x;
        h2d[i].s3.c1.im = (double) spin[6*i+5].y;
        h2d[i].s3.c2.re = (double) spin[6*i+5].z;
        h2d[i].s3.c2.im = (double) spin[6*i+5].w; 
        
  }
}



// convert spinor to REAL4 (float4, double4)

void convert2REAL4_spin_mpi (spinor * spin, dev_spinor * h2d, int start, int end) {

  int i;
  
  for (i = start; i < end; i++) {
    
        h2d[6*i+0].x = (float) spin[i].s0.c0.re;
        h2d[6*i+0].y = (float) spin[i].s0.c0.im;
        h2d[6*i+0].z = (float) spin[i].s0.c1.re;
        h2d[6*i+0].w = (float) spin[i].s0.c1.im;
        
        h2d[6*i+1].x = (float) spin[i].s0.c2.re;
        h2d[6*i+1].y = (float) spin[i].s0.c2.im;
        h2d[6*i+1].z = (float) spin[i].s1.c0.re;
        h2d[6*i+1].w = (float) spin[i].s1.c0.im;
        
        h2d[6*i+2].x = (float) spin[i].s1.c1.re;
        h2d[6*i+2].y = (float) spin[i].s1.c1.im;
        h2d[6*i+2].z = (float) spin[i].s1.c2.re;
        h2d[6*i+2].w = (float) spin[i].s1.c2.im;
        
        h2d[6*i+3].x = (float) spin[i].s2.c0.re;
        h2d[6*i+3].y = (float) spin[i].s2.c0.im;
        h2d[6*i+3].z = (float) spin[i].s2.c1.re;
        h2d[6*i+3].w = (float) spin[i].s2.c1.im;
        
        h2d[6*i+4].x = (float) spin[i].s2.c2.re;
        h2d[6*i+4].y = (float) spin[i].s2.c2.im;
        h2d[6*i+4].z = (float) spin[i].s3.c0.re;
        h2d[6*i+4].w = (float) spin[i].s3.c0.im;
        
        h2d[6*i+5].x = (float) spin[i].s3.c1.re;
        h2d[6*i+5].y = (float) spin[i].s3.c1.im;
        h2d[6*i+5].z = (float) spin[i].s3.c2.re;
        h2d[6*i+5].w = (float) spin[i].s3.c2.im;
    
  }
}



// cudaMemcpy gets  "spinor+6*offset"  because of pointer to float4 and there are 24 floats per site

void to_device_mpi (dev_spinor * device, spinor * host, dev_spinor * auxiliary, int size, int start, int end) {

  convert2REAL4_spin_mpi(host, auxiliary, start, end);					// auxiliary = (float) host
  cudaMemcpy(device+6*start, auxiliary+6*start, size, cudaMemcpyHostToDevice);		// device = auxiliary  (on device)

}


void to_host_mpi (spinor * host, dev_spinor * device, dev_spinor * auxiliary, int size, int start, int end) {

  cudaMemcpy(auxiliary+6*start, device+6*start, size, cudaMemcpyDeviceToHost);		// auxiliary = device  (on device)
  convert2double_spin_mpi(auxiliary, host, start, end);					// host = (double) auxiliary

}






// boundary exchange
//	all three versions do work:


/*
// preliminarily exchanges the full spinor field instead of only the boundaries

void xchange_field_wrapper (dev_spinor * dev_spin, int ieo) {

  size_t size = (VOLUME+RAND)/2 * 6*sizeof(dev_spinor);

  to_host_mpi(spinor_xchange, dev_spin, h2d_spin_up, size, 0, (VOLUME+RAND)/2);
  xchange_field(spinor_xchange, ieo);
  to_device_mpi(dev_spin, spinor_xchange, h2d_spin_dn, size, 0, (VOLUME+RAND)/2);

}
*/




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




// copies the boundary t-slices t=0 and t=T-1 to host
//	exchanges
//		copies RAND back to device

void xchange_field_wrapper (dev_spinor * dev_spin, int ieo) {
  
  #ifndef ALTERNATE_FIELD_XCHANGE
    
    size_t size_tSlice = LX*LY*LZ/2 * 6*sizeof(dev_spinor);
    size_t size_Rand   = RAND/2     * 6*sizeof(dev_spinor);
    
    to_host_mpi(spinor_xchange, dev_spin, h2d_spin_up, size_tSlice, 0 , LX*LY*LZ/2);
    to_host_mpi(spinor_xchange, dev_spin, h2d_spin_dn, size_tSlice, (T-1)*LX*LY*LZ/2, (VOLUME)/2);
    
    xchange_field(spinor_xchange, ieo);
    
    to_device_mpi(dev_spin, spinor_xchange, h2d_spin_up, size_Rand, VOLUME/2, (VOLUME+RAND)/2);
    
  #else
    
    int tSliceEO = LX*LY*LZ/2;
    int VolumeEO = VOLUME/2;
    
    cudaMemcpy(R1, dev_spin                      , tSliceEO*6*sizeof(float4), cudaMemcpyDeviceToHost);
    cudaMemcpy(R2, dev_spin+6*(VolumeEO-tSliceEO), tSliceEO*6*sizeof(float4), cudaMemcpyDeviceToHost);
    
    MPI_Sendrecv(R1, 24*tSliceEO, MPI_FLOAT, g_nb_t_dn, 0,
                 R3, 24*tSliceEO, MPI_FLOAT, g_nb_t_up, 0,
                 g_cart_grid, &stat[0]);
    MPI_Sendrecv(R2, 24*tSliceEO, MPI_FLOAT, g_nb_t_up, 1,
                 R4, 24*tSliceEO, MPI_FLOAT, g_nb_t_dn, 1,
                 g_cart_grid, &stat[1]);
    
    cudaMemcpy(dev_spin+6*VolumeEO           , R3, tSliceEO*6*sizeof(float4), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_spin+6*(VolumeEO+tSliceEO), R4, tSliceEO*6*sizeof(float4), cudaMemcpyHostToDevice);
    
  #endif
  
}

#endif	// MPI






/*

//////////////////////////
// SU(3) reconstruction //
//////////////////////////


// get 2 first rows of gf float4 type
//  
//
void su3to2vf4_mpi (su3** gf, dev_su3_2v* h2d_gf) {

  int i, j;
  
  for (i = 0; i < (VOLUME+RAND); i++){
    for (j = 0; j < 4; j++) {
    
      //first row
      h2d_gf[3*(4*i+j)].x = (float) gf[i][j].c00.re;
      h2d_gf[3*(4*i+j)].y = (float) gf[i][j].c00.im;
      h2d_gf[3*(4*i+j)].z = (float) gf[i][j].c01.re;
      h2d_gf[3*(4*i+j)].w = (float) gf[i][j].c01.im;
      h2d_gf[3*(4*i+j)+1].x = (float) gf[i][j].c02.re;
      h2d_gf[3*(4*i+j)+1].y = (float) gf[i][j].c02.im;      
      //second row
      h2d_gf[3*(4*i+j)+1].z = (float) gf[i][j].c10.re;
      h2d_gf[3*(4*i+j)+1].w = (float) gf[i][j].c10.im;
      h2d_gf[3*(4*i+j)+2].x = (float) gf[i][j].c11.re;
      h2d_gf[3*(4*i+j)+2].y = (float) gf[i][j].c11.im;
      h2d_gf[3*(4*i+j)+2].z = (float) gf[i][j].c12.re;
      h2d_gf[3*(4*i+j)+2].w = (float) gf[i][j].c12.im;      

    } 
  }
}




// bring gf into the form
// a2 a3, theta_a1, theta_c1, b1
// 
void su3to8_mpi (su3** gf, dev_su3_8* h2d_gf) {

  int i, j;
  
  for (i = 0; i < (VOLUME+RAND); i++) {
    for (j = 0; j < 4; j++) {
    
      // a2, a3
      h2d_gf[2*(4*i+j)].x = (float) gf[i][j].c01.re;
      h2d_gf[2*(4*i+j)].y = (float) gf[i][j].c01.im;
      h2d_gf[2*(4*i+j)].z = (float) gf[i][j].c02.re;
      h2d_gf[2*(4*i+j)].w = (float) gf[i][j].c02.im;
      
      // theta_a1, theta_c1
      // use atan2 for this: following the reference, atan2 should give an angle -pi < phi < +pi  
      h2d_gf[2*(4*i+j)+1].x = (float)( atan2((float) gf[i][j].c00.im,(float) gf[i][j].c00.re ));
      h2d_gf[2*(4*i+j)+1].y = (float) ( atan2((float) gf[i][j].c20.im,(float)gf[i][j].c20.re ));
      
      // b1
      h2d_gf[2*(4*i+j)+1].z = (float) gf[i][j].c10.re ;
      h2d_gf[2*(4*i+j)+1].w = (float) gf[i][j].c10.im ;
      
    } 
  }
}

*/





////////////////////
// linear algebra //
////////////////////

#ifdef MPI

// have to rebuilt some linear algebra functions which contain global communication
// can be done as wrappers to appropriate CUBLAS routines



// a wrapper function for cublasSdot() (with the same interface)
// provides the MPI communication via MPI_Allreduce()

float cublasSdot_wrapper(int size, float * A, int incx, float * B, int incy) {

  float result;
  float buffer;
  
  buffer = cublasSdot(size, (float *) A, incx, (float *) B, incy);
  MPI_Allreduce(&buffer, &result, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
  
  return(result);
  
}

#endif



















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
 **************************************************************************
 *
 * 	inspired by: Florian Burger
 * 	             Carsten Urbach
 *
 **************************************************************************/
 
 
 
 
#undef MPI
#undef REAL
  #include <mpi.h>
  //#include "../mpi_init.h"
#define MPI
#define REAL float



// optimization of the communication
#include "ASYNC.cuh"





// PRELIMINARY:
// compile the MPI implementation of the EO, ND mixed solver:
//	#include this file at the end of mixed_solve_eo_nd.cuh or after #including in mixed_solve.cu
//	adjust the makefile appropriately
//	make sure that _GAUGE_COPY and _USE_HALFSPINOR are NOT defined in config.h




	//////////////////////////////////////////////////////////////////
	//								//
	//  this is the MPI implementation of the eo, nd mixed solver	//
	//								//
	//	PARALLELT parallelization				//
	//	no _GAUGE_COPY and no _USE_HALFSPINOR			//
	//								//
	//////////////////////////////////////////////////////////////////






////////////////////////////////////
// host <--> device communication //
////////////////////////////////////


// convert spinor to double

void convert2double_spin_mpi (dev_spinor * spin, spinor * h2d, int start, int end) {

  int i;
  /*
  int Vol;
  
  if(even_odd_flag){
    Vol = (VOLUME+RAND)/2;
  }
  else{
    Vol = (VOLUME+RAND);
  }
  
  for (i = 0; i < Vol; i++) {
  */
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
  /*
  int Vol;
  
  if(even_odd_flag){
    Vol = (VOLUME+RAND)/2;
  }
  else{
    Vol = (VOLUME+RAND);
  }
  
  for (i = 0; i < Vol; i++) {
  */
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




// host/device interaction

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






#ifdef HOPPING_DEBUG

  // applies the hopping matrix on host for debugging purposes
  
  void Hopping_Matrix_wrapper (int ieo, dev_spinor * out, dev_spinor * in) {
  
    #ifdef MPI
      size_t size = (VOLUME+RAND)/2 * 6*sizeof(dev_spinor);
      
      //to_host(g_chi_up_spinor_field[DUM_SOLVER+3], in, h2d_spin_up, size);
      //Hopping_Matrix(ieo, g_chi_dn_spinor_field[DUM_SOLVER+3], g_chi_up_spinor_field[DUM_SOLVER+3]);
      //to_device(out, g_chi_dn_spinor_field[DUM_SOLVER+3], h2d_spin_up, size);
      
      to_host_mpi(spinor_debug_in, in, h2d_spin_up, size, 0, (VOLUME+RAND)/2);
      Hopping_Matrix(ieo, spinor_debug_out, spinor_debug_in);
      to_device_mpi(out, spinor_debug_out, h2d_spin_dn, size, 0, (VOLUME+RAND)/2);
    #else
      size_t size = VOLUME/2 * 6*sizeof(dev_spinor);
      
      to_host(spinor_debug_in, in, h2d_spin_up, size);
      Hopping_Matrix(ieo, spinor_debug_out, spinor_debug_in);
      to_device(out, spinor_debug_out, h2d_spin_dn, size);  
    #endif
    
    
  }

#endif






////////////////////
// linear algebra //
////////////////////

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
		// COMMENT: this is so far the only function with MPI AND(!) CUDA code
		//          either seperate this from the code and compile it appropriately
		//          or use the "--use-bindir"-option to nvcc in order to wrap mpicc with nvcc

		// PRELIMINARILY: was put to communication.c




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






/////////////////////////////////////////////
// geometry- and nearest-neighbour indices //
/////////////////////////////////////////////


// builds an array  iseven[global position]  to check wether is even or odd

void init_iseven() {

  int x0, x1, x2, x3;
  int ix;
  
  for (x0 = -1; x0 < T+1; x0++) {
    for (x1 = 0; x1 < LX; x1++) {
      for (x2 = 0; x2 < LY; x2++) {
        for (x3 = 0; x3 < LZ; x3++) {
          
          ix = Index(x0, x1, x2, x3);
          
	  if ((x0 + x1 + x2 + x3 + g_proc_coords[0]*T  + g_proc_coords[1]*LX + 
		                   g_proc_coords[2]*LY + g_proc_coords[3]*LZ) % 2 == 0) {
	    iseven[ix] = 1;
	  } 
	  else {
	    iseven[ix] = 0; 
	  }
        
        }}}}
     
}


// initialize nearest-neighbour table for gpu with even-odd enabled

void init_nnspinor_eo_mpi() {
									
  int x, y, z, t, ind, nnpos, j;					// mixed_solve_eo(...) allocates 8 integers per even or odd lattice site: size_t nnsize = 8*VOLUME*sizeof(int);
  									
  for (t = 0; t < T; t++) {						// loop goes over all INTERN latice sites !!
    for (x = 0; x < LX; x++) {						// doesn't refer to any EXTERN BOUNDARIES !!  ->  CORRESPONDS TO THE WHOLE LATTICE (I.E. WHEN NO SUBLATTICES ARE ASSIGNED) !!
      for (y = 0; y < LY; y++) {					//						  because of the behaviour of  g_iup[][] in the non-parallel case
        for (z = 0; z < LZ; z++) {
        								// NOTICE: g_ipt, g_iup, g_idn, and g_lexic2eosub  refer to pos. of lin. proj. pos. of the lattice
          ind = g_ipt[t][x][y][z];					// g_ipt[t][x][y][z] 	returns the linearly projected position of (t,x,y,z) of the lattice
          								//	indexes computed in geometry_eo() from geometry_eo.c
          								//	memory for the index array allocated by init_geometry_indices() from init_geometry_indices.c
          if ((t+x+y+z)%2 == 0) { // EVEN
          //if ((t + x + y + z + g_proc_coords[0]*T  + g_proc_coords[1]*LX + 
	  //	               g_proc_coords[2]*LY + g_proc_coords[3]*LZ) % 2 == 0) {
          
            nnpos = g_lexic2eosub[ind];					// g_lexic2eosub[ind] 	returns the position of [ind] in the sub-eo-notation
            										      ////////////////
            for (j = 0; j < 4; j++) { // plus direction			// here are also the  // BOUNDARIES //  included and properly mapped:
            								//		      ////////////////
              nn_eo[8*nnpos+j] = g_lexic2eosub[ g_iup[ind][j] ];	// g_iup[ind][j] 	returns the position of the nearest neighbour of [ind] in direction +[j]
            }								//				-->  for the non-parallized code g_iup[][] maps INTERN !!
            for (j = 0; j < 4; j++) { // minus direction
              nn_eo[8*nnpos+4+j] = g_lexic2eosub[ g_idn[ind][j] ];	// g_idn[ind][j] 	returns the position of the nearest neighbour of [ind] in direction -[j]
            }
          }
          
          else {		  // ODD
          
            nnpos = g_lexic2eosub[ind];
            
            for (j = 0; j < 4; j++) { // plus direction
              nn_oe[8*nnpos+j] = g_lexic2eosub[ g_iup[ind][j] ];	// nn_oe	      will return the nearest neigbours
            }								// nn_eo  and  nn_oe  strictly refer to the 4d-spacetime lattice
            
            for (j = 0; j < 4; j++) { // minus direction
              nn_oe[8*nnpos+4+j] = g_lexic2eosub[ g_idn[ind][j] ];
            }
          }
  }}}} // for loops
}




// the following functions can all be used to properly initialize the fields  eoidx_even[]  and  eoidx_odd[]  for addressing the gauge fields:


void init_idxgauge_mpi() {		// works!

  int t, x, y, z;
  int pos_eo, pos_global;
  
  for (t = -1; t < T+1; t++) {
    for (x = 0; x < LX; x++) {
      for (y = 0; y < LY; y++) {
        for (z = 0; z < LZ; z++) {
        
        //pos_global = g_ipt[t][x][y][z];
        pos_global = Index(t,x,y,z);
        pos_eo     = g_lexic2eosub[pos_global];
        
        //if ((t+x+y+z)%2 == 0) { // EVEN
        if ((t + x + y + z + g_proc_coords[0]*T  + g_proc_coords[1]*LX + 
	  	             g_proc_coords[2]*LY + g_proc_coords[3]*LZ) % 2 == 0) {
          eoidx_even[pos_eo] = g_eo2lexic[pos_eo];
        }
        else  {			// ODD
          eoidx_odd[pos_eo] = g_eo2lexic[(VOLUME+RAND)/2+pos_eo];
        }
  }}}} // for loop over the INTERN lattice
  
  //printf("This was init_idxgauge_mpi().\n");
  
}



/*
void init_idxgauge_mpi() {		// works!

  int t, x, y, z;
  int pos_eo, pos_global;
  
  for (t = -1; t < T+1; t++) {
    for (x = 0; x < LX; x++) {
      for (y = 0; y < LY; y++) {
        for (z = 0; z < LZ; z++) {
        
        //pos_global = g_ipt[t][x][y][z];
        pos_global = Index(t,x,y,z);
        pos_eo     = g_lexic2eosub[pos_global];
        
        //if ((t+x+y+z)%2 == 0) { // EVEN
        if ((t + x + y + z + g_proc_coords[0]*T  + g_proc_coords[1]*LX + 
	  	             g_proc_coords[2]*LY + g_proc_coords[3]*LZ) % 2 == 0) {
          eoidx_even[pos_eo] = pos_global;
        }
        else  {			// ODD
          eoidx_odd[pos_eo] = pos_global;
        }
  }}}} // for loop over the INTERN lattice
  
  //printf("This was init_idxgauge_mpi().\n");
  
}
*/


/*
void init_idxgauge_mpi() {		// works!

  int pos_eo, pos_global_even, pos_global_odd;
  
  for (pos_eo = 0; pos_eo < (VOLUME+RAND)/2; pos_eo++) {
      // even
      pos_global_even = g_eo2lexic[pos_eo];
      eoidx_even[pos_eo] = pos_global_even;
      // odd
      pos_global_odd = g_eo2lexic[(VOLUME+RAND)/2 + pos_eo];
      eoidx_odd[pos_eo] = pos_global_odd;
  }
  
  //printf("This was init_idxgauge_mpi().\n");
  
}
*/


/*
void init_idxgauge_mpi() {		// works!

  int pos_eo, pos_global;
  
  for (pos_global = 0; pos_global < (VOLUME+RAND); pos_global++) {
  
    pos_eo = g_lexic2eosub[pos_global];
    
    if (iseven[pos_global] == 1) {
    //if (pos_global%2 == 0) {
      eoidx_even[pos_eo] = pos_global;
    }
    else {
      eoidx_odd[pos_eo]  = pos_global;
    }
      
  }
  
  //printf("This was init_idxgauge_mpi().\n");
  
}
*/


/*
void init_idxgauge_mpi() {		// works!

  int x, y, z, t;
  int ind;
  int evenpos = 0;
  int oddpos = 0;
  
  for (t = 0; t < T; t++) {
    for (x = 0; x < LX; x++) {
      for (y = 0; y < LY; y++) {
        for (z = 0; z < LZ; z++) {
          ind = g_ipt[t][x][y][z];
          if ((t+x+y+z) % 2 == 0) {
            eoidx_even[evenpos] = ind;
            evenpos++;
          }
          else {
            eoidx_odd[oddpos] = ind;
            oddpos++;
          }
  }}}} // INTERN
  
  
  		t = T;
  		  for (x = 0; x < LX; x++) {
  		    for (y = 0; y < LY; y++) {
  		      for (z = 0; z < LZ; z++) {
  		        ind = VOLUME + z + LZ*y + LZ*LY*x;
  		        //if (iseven[ind] == 1) {
  		        if ((t+x+y+z) % 2 == 0) {
  		          eoidx_even[evenpos] = ind;
  		          evenpos++;
  		        }
  		        else {
  		          eoidx_odd[oddpos] = ind;
  		          oddpos++;
  		        }
  		}}} // EXTERN
  
  
  				t = -1;
  				  for (x = 0; x < LX; x++) {
  				    for (y = 0; y < LY; y++) {
  				      for (z = 0; z < LZ; z++) {
  				        ind = VOLUME + LX*LY*LZ + z + LZ*y + LZ*LY*x;
  				        //if (iseven[ind] == 1) {
  				        if ((t+x+y+z) % 2 == 0) {
  				          eoidx_even[evenpos] = ind;
  				          evenpos++;
  				        }
  				        else {
  				          eoidx_odd[oddpos] = ind;
  				          oddpos++;
  				        }
  				}}} // EXTERN
  				
  //printf("This was init_idxgauge_mpi().\n");
  
}
*/



/////////////////////
// initializations //
/////////////////////


// allocates memory for the fields for the alternative way of addressing positions in dev_Hopping_Matrix_alternate()

void init_gpu_indexfields() {
  
  size_t size;
  
  		// debug
  		/*
  		printf("Test: %p ?= %p\n", *g_iup, g_iup[0]);
  		printf("Test: %p ?= %p\n", *g_idn, g_idn[0]);
  		printf("Test: %p ?= %p\n", g_lexic2eo, &g_lexic2eo[0]);
  		printf("Test: %p ?= %p\n", g_lexic2eosub, &g_lexic2eosub[0]);
  		printf("Test: %p ?= %p\n", g_eo2lexic, &g_eo2lexic[0]);
  		printf("Test: %p ?= %p ?= %p\n", ***g_ipt, &g_ipt[0][0][0][0], g_ipt[0][0][0]);
  		*/
  
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




// puts the additional variables VOLUMEPLUSRAND and RAND on the device
__global__ void he_cg_init_nd_additional_mpi (int param_VOLUMEPLUSRAND, int param_RAND, int rank, int nproc) {

  dev_VOLUMEPLUSRAND  = param_VOLUMEPLUSRAND;
  dev_RAND            = param_RAND;
  
  dev_rank            = rank;
  dev_nproc           = nproc;

}





// code to list available devices, not yet included in main code
// this is copied from the CUDA sdk 
extern "C" int find_devices_mpi() {

int deviceCount, dev;

    cudaGetDeviceCount(&deviceCount);
    if (deviceCount == 0)
        if (g_cart_id == 0) printf("There is no device supporting CUDA\n");
    for (dev = 0; dev < deviceCount; ++dev) {
        cudaDeviceProp deviceProp;
        cudaGetDeviceProperties(&deviceProp, dev);
        if (dev == 0) {
            if (deviceProp.major == 9999 && deviceProp.minor == 9999)
                if (g_cart_id == 0) printf("There is no device supporting CUDA.\n");
            else if (deviceCount == 1)
                if (g_cart_id == 0) printf("There is 1 device supporting CUDA\n");
            else
                if (g_cart_id == 0) printf("There are %d devices supporting CUDA\n", deviceCount);
        }
        if (g_cart_id == 0) {
          printf("\nDevice %d: \"%s\"\n", dev, deviceProp.name);
          printf("  Major revision number:                         %d\n",
               deviceProp.major);
          printf("  Minor revision number:                         %d\n",
               deviceProp.minor);
          printf("  Total amount of global memory:                 %u bytes\n",
               deviceProp.totalGlobalMem);
    #if CUDART_VERSION >= 2000
          printf("  Number of multiprocessors:                     %d\n",
               deviceProp.multiProcessorCount);
          printf("  Number of cores:                               %d\n",
               8 * deviceProp.multiProcessorCount);
    #endif
          printf("  Total amount of constant memory:               %u bytes\n",
               deviceProp.totalConstMem); 
          printf("  Total amount of shared memory per block:       %u bytes\n",
               deviceProp.sharedMemPerBlock);
          printf("  Total number of registers available per block: %d\n",
               deviceProp.regsPerBlock);
          printf("  Warp size:                                     %d\n",
               deviceProp.warpSize);
          printf("  Maximum number of threads per block:           %d\n",
               deviceProp.maxThreadsPerBlock);
          printf("  Maximum sizes of each dimension of a block:    %d x %d x %d\n",
               deviceProp.maxThreadsDim[0],
               deviceProp.maxThreadsDim[1],
               deviceProp.maxThreadsDim[2]);
          printf("  Maximum sizes of each dimension of a grid:     %d x %d x %d\n",
               deviceProp.maxGridSize[0],
               deviceProp.maxGridSize[1],
               deviceProp.maxGridSize[2]);
          printf("  Maximum memory pitch:                          %u bytes\n",
               deviceProp.memPitch);
          printf("  Texture alignment:                             %u bytes\n",
               deviceProp.textureAlignment);
          printf("  Clock rate:                                    %.2f GHz\n",
               deviceProp.clockRate * 1e-6f);
    #if CUDART_VERSION >= 2000
          printf("  Concurrent copy and execution:                 %s\n",
               deviceProp.deviceOverlap ? "Yes" : "No");
    #endif
      }
    }
    return(deviceCount);
}







// initializes and allocates all quantities for the mixed solver
// more precise:
//	puts the gauge field on device as "2 rows" or "8 floats" per SU(3)-matrix
//	allocates memory for all spinor fields
//	puts the nn- and eoidx-fields on device memory
//
// for the MPI purposes we HAVE TO:
//	choose one of several GPUs according to an appropriate mechanism	// still to do !
//	mainly: VOLUME --> (VOLUME+RAND)
//	maybe: grid[5] = (VOLUME+RAND)/2  ??

void init_mixedsolve_eo_nd_mpi(su3** gf) {	// gf is the full gauge field
  
  
  /////////////////////
  // LOCAL VARIABLES //
  /////////////////////
  cudaError_t cudaerr;		// CUDA errors
  int ndev;			// number of devices
  size_t dev_gfsize;		// size of the gauge-field on device memory
  size_t nnsize;		// size of memory for nn-table
  size_t idxsize;		// size of memory for even/odd-positions
  size_t dev_spinsize;		// size of memory for spinors
  int grid[5];			// array for grid specifications
  float * host_output;		// ??
  
  /*
  #ifdef GF_8
    dev_su3_8 * h2d_gf;
  #else
    dev_su3_2v * h2d_gf;
  #endif
  */
  
  /*
  int * nn;
  int * nn_eo;
  int * nn_oe;
  int * eoidx_even;
  int * eoidx_odd;
  */
  
  
  
  
  // get number of devices										// HAVE TO: mechanism to choose the device
  ndev = find_devices_mpi();
  if (ndev == 0) {
    fprintf(stderr, "Process %d of %d: Error: no CUDA devices found. Aborting...\n", g_proc_id, g_nproc);
    exit(300);
  }
  #ifndef DEVICE_EQUAL_RANK
    // try to set active device to device_num given in input file
    // each process gets bounded to the same GPU , preliminary !!
    if (device_num < ndev) {
      printf("Process %d of %d: Setting active device to: %d\n", g_proc_id, g_nproc, device_num);
      cudaSetDevice(device_num);
    }
    else {
      fprintf(stderr, "Process %d of %d: Error: There is no CUDA device with No. %d. Aborting...\n", g_proc_id, g_nproc, device_num);
      exit(301);
    }
  #else
    // device number = mpi rank
    if (g_cart_id < ndev) {
      printf("Process %d of %d: Setting active device to: %d\n", g_proc_id, g_nproc, g_cart_id);
      cudaSetDevice(g_cart_id);
    }
    else {
      fprintf(stderr, "Process %d of %d: Error: There is no CUDA device with No. %d. Aborting...\n", g_proc_id, g_nproc, g_cart_id);
      exit(301);
    }
  #endif
  if ((cudaerr=cudaGetLastError()) != cudaSuccess) {
    printf("Process %d of %d: Error in init_mixedsolve_eo_nd(): Could not set active device. Aborting...\n", g_proc_id, g_nproc);
    exit(302);
  }
  
  
  
  
  #ifdef USETEXTURE
    printf("Process %d of %d: Using texture references.\n", g_proc_id, g_nproc);
  #else
    printf("Process %d of %d: NOT using texture references.\n", g_proc_id, g_nproc);
  #endif
  
  
  
  
  /////////////////
  // GAUGE FIELD //
  /////////////////
  #ifdef GF_8
    printf("Process %d of %d: Using GF 8 reconstruction.\n", g_proc_id, g_nproc);			// dev_su3_8 = float4
    dev_gfsize = 4*(VOLUME+RAND) * 2*sizeof(dev_su3_8);		// allocates for each lattice site and RAND for 4 directions  2*float4 = 8 floats  = 8 real parameters
  #else
    printf("Process %d of %d: Using GF 12 reconstruction.\n", g_proc_id, g_nproc);			// dev_su3_2v = float4
    dev_gfsize = 4*(VOLUME+RAND) * 3*sizeof(dev_su3_2v); 	// allocates for each lattice site and RAND for 4 directions  3*float4 = 12 floats = 2 rows of complex 3-vectors
  #endif
  
  
  if ( (cudaerr = cudaMalloc((void **) &dev_gf, dev_gfsize)) != cudaSuccess ){	// allocates memory for the gauge field dev_gf on device
    printf("Process %d of %d: Error in init_mixedsolve_eo_nd_mpi(): Memory allocation of gauge field failed. Aborting...\n", g_proc_id, g_nproc);
    printf("%s\n", cudaGetErrorString(cudaGetLastError()));
    exit(200);
  }
  else{
    printf("Process %d of %d: Allocated gauge field on device.\n", g_proc_id, g_nproc);
  }  
  
  
  #ifdef GF_8
    h2d_gf = (dev_su3_8 *) malloc(dev_gfsize); 			// allocates on host
    su3to8_mpi(gf, h2d_gf);					// h2d_gf  is the gauge field  gf  with the 8-real-parameter-representation (according to M. Clark, p. 28)
  #else
    h2d_gf = (dev_su3_2v *) malloc(dev_gfsize);			// allocates on host
    su3to2vf4_mpi(gf, h2d_gf);					// h2d_gf  is the gauge field  gf  with the first two rows stored
  #endif
  
  cudaMemcpy(dev_gf, h2d_gf, dev_gfsize, cudaMemcpyHostToDevice);
  								// dev_gf = h2d_gf  on device memory
  
  		// debug	// CUDA
  		#ifdef CUDA_DEBUG
  		  printf("Process %d of %d: ", g_proc_id, g_nproc);
  		  CUDA_CHECK("CUDA error in init_mixedsolve_eo_nd(). Copying dev_gf to device failed.", "Allocated dev_gf on device.");
  		#endif
  
  
  
  
  //////////
  // GRID //
  //////////
  nnsize = 8*(VOLUME)*sizeof(int);				// size of memory for 8*VOLUME integers
  nn = (int *) malloc(nnsize);					//	4*2 global lin. proj. nearest neighbour positions
  nn_eo = (int *) malloc(nnsize/2);
  nn_oe = (int *) malloc(nnsize/2);				// NOTICE: here we don't need VOLUME+RAND !!
  cudaMalloc((void **) &dev_nn, nnsize);
  cudaMalloc((void **) &dev_nn_eo, nnsize/2);
  cudaMalloc((void **) &dev_nn_oe, nnsize/2);
  
  
  idxsize = (VOLUME+RAND)/2*sizeof(int);			// size of memory necessary for VOLUME/2 integers
  eoidx_even = (int *) malloc(idxsize);				//	gobal lin. proj. positions of the even- or odd sublattice
  eoidx_odd = (int *) malloc(idxsize);
  cudaMalloc((void **) &dev_eoidx_even, idxsize);		// NOTICE: here wo do need VOLUME+RAND !!
  cudaMalloc((void **) &dev_eoidx_odd, idxsize);
  
  
  initnn();							// initialize nearest-neighbour table for gpu
  //initnn_eo();
  
  //iseven = (int *) malloc((VOLUME+RAND)*sizeof(int));
  //init_iseven();
  
  init_nnspinor_eo_mpi();					// initialize nearest-neighbour table for gpu with even-odd enabled
  init_idxgauge_mpi();
  
  #ifdef ALTERNATE_HOPPING_MATRIX
    init_gpu_indexfields();
  #endif
  
  
  cudaMemcpy(dev_nn, nn, nnsize, cudaMemcpyHostToDevice);	// copies the previous initialized index-arrays from host to device memory
  cudaMemcpy(dev_nn_eo, nn_eo, nnsize/2, cudaMemcpyHostToDevice);
  cudaMemcpy(dev_nn_oe, nn_oe, nnsize/2, cudaMemcpyHostToDevice);
  cudaMemcpy(dev_eoidx_even, eoidx_even, idxsize, cudaMemcpyHostToDevice);
  cudaMemcpy(dev_eoidx_odd, eoidx_odd, idxsize, cudaMemcpyHostToDevice);
  
  
  free(eoidx_odd);						// deallocates the host memory for the field
  free(eoidx_even);						// they are only on the device
  free(nn_oe);
  free(nn_eo);							// not necessary for locally defined variables ??
  free(nn);
  
  
  		// debug	// CUDA
  		#ifdef CUDA_DEBUG
  		  printf("Process %d of %d: ", g_proc_id, g_nproc);
  		  CUDA_CHECK("CUDA error in init_mixedsolve_eo_nd(). Memory allocation of grid stuff failed.", "Allocated grid stuff on device.");
  		#endif
  
  
  
  
  /////////////
  // SPINORS //							// allocates device memory for the odd part of the spinor fields (24 = 4*6 floats per odd lattice sites)
  /////////////							// now we have to consider 2 flavors: up, dn
  
  dev_spinsize = (VOLUME+RAND)/2 * 6*sizeof(dev_spinor);	// remember: dev_spinor = float4	// NOTICE: this refers to the memory requirements for the device, host needs twice the memory !!
  
  cudaMalloc((void **) &dev_spin1_up, dev_spinsize);   		// allocates device memory for the fields spinor fields used in dev_cg_eo_nd(...)
  cudaMalloc((void **) &dev_spin1_dn, dev_spinsize);		// pointing to device
  cudaMalloc((void **) &dev_spin2_up, dev_spinsize);		// ...
  cudaMalloc((void **) &dev_spin2_dn, dev_spinsize);
  cudaMalloc((void **) &dev_spin3_up, dev_spinsize);
  cudaMalloc((void **) &dev_spin3_dn, dev_spinsize);
  /*
  cudaMalloc((void **) &dev_spin4_up, dev_spinsize);		// not needed
  cudaMalloc((void **) &dev_spin4_dn, dev_spinsize);
  cudaMalloc((void **) &dev_spin5_up, dev_spinsize);
  cudaMalloc((void **) &dev_spin5_dn, dev_spinsize);
  */
  
  cudaMalloc((void **) &dev_spinin_up , dev_spinsize);		// host/device interaction
  cudaMalloc((void **) &dev_spinin_dn , dev_spinsize);		// inner/outer interaction
  cudaMalloc((void **) &dev_spinout_up, dev_spinsize);
  cudaMalloc((void **) &dev_spinout_dn, dev_spinsize);
  
  
  		// debug	// host code
  		if ( (void *) (h2d_spin_up = (dev_spinor *) malloc(dev_spinsize) ) == NULL) {						// MEMORY REQUIREMENTS: these are auxiliary fields for  to_host()  and  to_device()
  		  printf("Process %d of %d: Could not allocate memory for h2d_spin_up. Aborting...\n", g_proc_id, g_nproc);		//                      they have to store floats (not doubles)
  		  exit(200);
  		}
  		if ( (void *) (h2d_spin_dn = (dev_spinor *) malloc(dev_spinsize) ) == NULL) {
  		  printf("Process %d of %d: Could not allocate memory for h2d_spin_dn. Aborting...\n", g_proc_id, g_nproc);
  		  exit(200);
  		}
  // h2d_spin_up = (dev_spinor *) malloc(dev_spinsize);		// for transferring the spin field in double precision on host to single precision on device
  // h2d_spin_dn = (dev_spinor *) malloc(dev_spinsize);		// on host pointing to host
  
  #if defined(ALTERNATE_FIELD_XCHANGE) || defined(ASYNC_OPTIMIZED)
    int tSliceEO = LX*LY*LZ/2;
  #endif
  
  		#ifndef ALTERNATE_FIELD_XCHANGE
  		  // debug	// host code
  		  if ( (void *) (spinor_xchange = (spinor *) malloc(2*dev_spinsize) ) == NULL) {						// MEMORY REQUIREMENTS: auxiliary fields for   xchange_field_wrapper()  and  Hopping_Matrix_wrapper()
  		    printf("Process %d of %d: Could not allocate memory for spinor_xchange. Aborting...\n", g_proc_id, g_nproc);		//                      have to store doubles --> 2*dev_spinsize  !!
  		    exit(200);
  		  }
  		#else
  		  // int tSliceEO = LX*LY*LZ/2;
  		  R1 = (dev_spinor *) malloc(2*tSliceEO*24*sizeof(float));
  		  R2 = R1 + 6*tSliceEO;
  		  R3 = (dev_spinor *) malloc(2*tSliceEO*24*sizeof(float));
  		  R4 = R3 + 6*tSliceEO;
  		#endif
  		
  		#ifdef HOPPING_DEBUG
  		  // debug	// host code
  		  if ( (void *) (spinor_debug_in = (spinor *) malloc(2*dev_spinsize) ) == NULL) {
  		    printf("Process %d of %d: Could not allocate memory for spinor_debug_in. Aborting...\n", g_proc_id, g_nproc);
  		    exit(200);
  		  }
  		  // debug	// host code
  		  if ( (void *) (spinor_debug_out = (spinor *) malloc(2*dev_spinsize) ) == NULL) {
  		    printf("Process %d of %d: Could not allocate memory for spinor_debug_out. Aborting...\n", g_proc_id, g_nproc);
  		    exit(200);
  		  }
  		#endif
  
  
  cudaMalloc((void **) &dev_spin_eo1_up, dev_spinsize);		// used for matrix_multiplication32(...)
  cudaMalloc((void **) &dev_spin_eo1_dn, dev_spinsize);
  cudaMalloc((void **) &dev_spin_eo3_up, dev_spinsize);
  cudaMalloc((void **) &dev_spin_eo3_dn, dev_spinsize);
  /*
  cudaMalloc((void **) &dev_spin_eo2_up, dev_spinsize);		// no memory allocation needed
  cudaMalloc((void **) &dev_spin_eo2_dn, dev_spinsize);		// will point to already allocated memory when used in matrix_multiplication
  */
  
  
  		// debug	// CUDA
  		#ifdef CUDA_DEBUG
  		  printf("Process %d of %d: ", g_proc_id, g_nproc);
  		  CUDA_CHECK("CUDA error in init_mixedsolve_eo_nd(). Memory allocation of spinor fields failed.", "Allocated spinor fields on device.");
  		#endif
  
  
  #ifdef ASYNC_OPTIMIZED					// for exchanging the boundaries in the MPI code
    // int tSliceEO = LX*LY*LZ/2;
    cudaMallocHost(&RAND3, 2*tSliceEO*6*sizeof(float4));
    RAND4 = RAND3 + 6*tSliceEO;
    // RAND1 = (dev_spinor *) malloc(2*tSliceEO*6*sizeof(float4));
    // RAND2 = RAND1 + 6*tSliceEO;
    cudaMallocHost(&RAND1, 2*tSliceEO*6*sizeof(float4));
    RAND2 = RAND1 + 6*tSliceEO;
    
    for (int i = 0; i < 2*nStreams+1; i++) {
      cudaStreamCreate(&stream[i]);
    }
    
    #ifdef ASYNC_TIMING
      cudaEventCreate(&start_ALL);
      cudaEventCreate(&stop_ALL);
      cudaEventCreate(&stop_D2H_1);
      cudaEventCreate(&stop_D2H_2);
      cudaEventCreate(&stop_INT_0);
      cudaEventCreate(&stop_H2D_3);
      cudaEventCreate(&stop_H2D_4);
      cudaEventCreate(&stop_EXT_1);
      cudaEventCreate(&stop_EXT_2);
    #endif
  #endif
  
  
  
  
  ////////////
  // output //						// ??
  ////////////
  output_size = LZ*T*sizeof(float); 			// parallel in t and z direction
  cudaMalloc((void **) &dev_output, output_size);	// output array
  host_output = (float *) malloc(output_size);
  
  
  		// debug	// CUDA
  		#ifdef CUDA_DEBUG
  		  printf("Process %d of %d: ", g_proc_id, g_nproc);
  		  CUDA_CHECK("CUDA error in init_mixedsolve_eo_nd(). Memory allocation output stuff failed.", "Allocated output stuff on device.");
  		#endif
  
  
  
  										// HAVE TO: maybe set grid[5] = (VOLUME+RAND)/2 ??	// no because refers to INTERN lattice sites !!
  ////////////////////////////
  // grid[ ] specifications //							// allocate and initializes the array grid[5] on device
  ////////////////////////////
  grid[0] = LX;									// it contains the dimensions of the lattice and the volume of the eo-sublattice
  grid[1] = LY;
  grid[2] = LZ;
  grid[3] = T;
  grid[4] = VOLUME/2;								// will be used to set dev_VOLUME: dev_VOLUME is half of VOLUME for eo
  
  cudaMalloc((void **) &dev_grid, 5*sizeof(int));				// dev_grid
  cudaMemcpy(dev_grid, &(grid[0]), 5*sizeof(int), cudaMemcpyHostToDevice);
  
  		// debug	// CUDA
  		#ifdef CUDA_DEBUG
  		  printf("Process %d of %d: ", g_proc_id, g_nproc);
  		  CUDA_CHECK("CUDA error in init_mixedsolve_eo_nd(). Memory allocation of grid[] specifications failed.", "Allocated grid[] specifications on device.");
  		#endif
  
  
  // MPI_Barrier(g_cart_grid);
  
  
  
  
  
  
  
}//init_mixedsolve_eo_nd()




// deallocates the previous allocated memory

void finalize_mixedsolve_eo_nd_mpi(void) {
  
  cudaError_t cudaerr;
  
  cudaFree(dev_spin1_up);
  cudaFree(dev_spin1_dn);
  cudaFree(dev_spin2_up);
  cudaFree(dev_spin2_dn);
  cudaFree(dev_spin3_up);
  cudaFree(dev_spin3_dn);
  /*
  cudaFree(dev_spin4_up);
  cudaFree(dev_spin4_dn);
  cudaFree(dev_spin5_up);
  cudaFree(dev_spin5_dn);
  */
  
  cudaFree(dev_spinin_up);
  cudaFree(dev_spinin_dn);
  cudaFree(dev_spinout_up);
  cudaFree(dev_spinout_dn);
  
  free(h2d_spin_up);
  free(h2d_spin_dn);
  #ifndef ALTERNATE_FIELD_XCHANGE
    free(spinor_xchange);
  #else
    free(R1);
    free(R3);
  #endif
  #ifdef HOPPING_DEBUG
    free(spinor_debug_in);
    free(spinor_debug_out);
  #endif
  
  cudaFree(dev_spin_eo1_up);
  cudaFree(dev_spin_eo1_dn);
  cudaFree(dev_spin_eo3_up);
  cudaFree(dev_spin_eo3_dn);
  /*
  cudaFree(dev_spin_eo2_up);
  cudaFree(dev_spin_eo2_dn);
  */
  
  
  cudaFree(dev_nn);
  cudaFree(dev_nn_eo);
  cudaFree(dev_nn_oe);
  cudaFree(dev_eoidx_even);
  cudaFree(dev_eoidx_odd);
  

  cudaFree(dev_gf);
  cudaFree(dev_output);
  cudaFree(dev_grid);
  
  
  free(h2d_gf);
  
  #ifdef ALTERNATE_HOPPING_MATRIX
    free_gpu_indexfields();
  #endif
  
  #ifdef ASYNC_OPTIMIZED
    cudaFreeHost(RAND3);
    // free(RAND1);
    cudaFreeHost(RAND1);

    for (int i = 0; i < 2*nStreams+1; i++) {
      cudaStreamDestroy(stream[i]);
    }
    
    #ifdef ASYNC_TIMING
      cudaEventDestroy(start_ALL);
      cudaEventDestroy(stop_ALL);
      cudaEventDestroy(stop_D2H_1);
      cudaEventDestroy(stop_D2H_2);
      cudaEventDestroy(stop_INT_0);
      cudaEventDestroy(stop_H2D_3);
      cudaEventDestroy(stop_H2D_4);
      cudaEventDestroy(stop_EXT_1);
      cudaEventDestroy(stop_EXT_2);
    #endif
  #endif
  
  
  // Clean up CUDA API for calling thread	// ??
  cudaThreadExit();				// is essential
  
  
  		// debug	// CUDA
  		#ifdef CUDA_DEBUG
  		  CUDA_CHECK("CUDA error in finalize_mixedsolve_eo_nd(). Device memory deallocation failed", "Device memory deallocated.");
  		#endif
  
}//finalize_mixedsolve_eo_nd_mpi()





/*
////////////////////
// hopping matrix //
////////////////////


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
*/








///////////////////////////
// MATRIX MULTIPLICATION //
///////////////////////////

// the GPU implementation of  Q_Qdagger_ND(...)  from Nondegenerate_Matrix.c
//	Flo's equivalent function for the standard and non-nd case is  dev_Qtm_pm_psi

void matrix_multiplication32_mpi (dev_spinor * spinout_up, dev_spinor * spinout_dn,
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
  
  
  ///////////////////////////////////////							/////////////////////////////////
  //        Q_tilde_dagger(2x2)        //							// (a,b) = (Q_tilde) * (dn,up) //
  ///////////////////////////////////////							/////////////////////////////////
  
  
  ////////////////////////////
  // (M_oe) (Mee^-1) (M_eo) //
  ////////////////////////////
  
  
  		// xchange
  		#ifndef HOPPING_DEBUG
  		  xchange_field_wrapper(spinin_dn, 0);
  		#endif
  
  // hopping matrix
  #ifdef USETEXTURE
    bind_texture_spin(spinin_dn,1);
  #endif
  
    #ifndef HOPPING_DEBUG
    	#ifndef ALTERNATE_HOPPING_MATRIX
    	  dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, spinin_dn, dev_spin_eo1_up, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0);	// dev_spin_eo1_up = (M_eo) * spinin_dn
    	#else
    	  dev_Hopping_Matrix_alternate<<<gridsize1, blocksize1>>>(dev_gf, spinin_dn, dev_spin_eo1_up, dev_g_iup, dev_g_idn, dev_g_eo2lexic, dev_g_lexic2eosub, 0);
    	#endif
    #else
      Hopping_Matrix_wrapper(0, dev_spin_eo1_up, spinin_dn);
    #endif
  
  #ifdef USETEXTURE
    unbind_texture_spin(1);
  #endif
  
  
  		// xchange
  		#ifndef HOPPING_DEBUG
  		  xchange_field_wrapper(spinin_up, 0);
  		#endif
  
  // hopping matrix
  #ifdef USETEXTURE
    bind_texture_spin(spinin_up,1);
  #endif
  
    #ifndef HOPPING_DEBUG
    	#ifndef ALTERNATE_HOPPING_MATRIX
    	  dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, spinin_up, dev_spin_eo1_dn, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0);	// dev_spin_eo1_dn = (M_eo) * spinin_up
    	#else
    	  dev_Hopping_Matrix_alternate<<<gridsize1, blocksize1>>>(dev_gf, spinin_up, dev_spin_eo1_dn, dev_g_iup, dev_g_idn, dev_g_eo2lexic, dev_g_lexic2eosub, 0);
    	#endif
    #else
      Hopping_Matrix_wrapper(0, dev_spin_eo1_dn, spinin_up);
    #endif
  
  #ifdef USETEXTURE
    unbind_texture_spin(1);
  #endif
  
  
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
  
  
  		// xchange
  		#ifndef HOPPING_DEBUG
  		  xchange_field_wrapper(dev_spin_eo2_up, 1);
  		#endif
  
  // hopping matrix
  #ifdef USETEXTURE
    bind_texture_spin(dev_spin_eo2_up,1);									// remember: this is ((M_oe)(Mee^-1)(M_eo)) * (spinin_dn, spinin_up):
  #endif
  
    #ifndef HOPPING_DEBUG
    	#ifndef ALTERNATE_HOPPING_MATRIX
    	  dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo2_up, dev_spin_eo1_up, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1);
    	#else
    	  dev_Hopping_Matrix_alternate<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo2_up, dev_spin_eo1_up, dev_g_iup, dev_g_idn, dev_g_eo2lexic, dev_g_lexic2eosub, 1);
    	#endif
    #else
      Hopping_Matrix_wrapper(1, dev_spin_eo1_up, dev_spin_eo2_up);
    #endif
  
  #ifdef USETEXTURE
    unbind_texture_spin(1);							// dev_spin_eo1_up  =  (M_oe) * dev_spin_eo2_up  =  (M_oe)*nrm*(1-imubar)*(M_eo) * spinin_dn  +  (M_oe)*nrm*epsbar*(M_eo) * spinin_up
  #endif									
  
  
  		// xchange
  		#ifndef HOPPING_DEBUG
  		  xchange_field_wrapper(dev_spin_eo2_dn, 1);
  		#endif
  
  // hopping matrix
  #ifdef USETEXTURE
    bind_texture_spin(dev_spin_eo2_dn,1);
  #endif
  
    #ifndef HOPPING_DEBUG
    	#ifndef ALTERNATE_HOPPING_MATRIX
    	  dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo2_dn, dev_spin_eo1_dn, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1);
    	#else
    	  dev_Hopping_Matrix_alternate<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo2_dn, dev_spin_eo1_dn, dev_g_iup, dev_g_idn, dev_g_eo2lexic, dev_g_lexic2eosub, 1);
    	#endif
    #else
      Hopping_Matrix_wrapper(1, dev_spin_eo1_dn, dev_spin_eo2_dn);
    #endif
    
  #ifdef USETEXTURE
    unbind_texture_spin(1);							// dev_spin_eo1_dn  =  (M_oe) * dev_spin_eo2_dn  =  (M_oe)*nrm*(1+imubar)*(M_eo) * spinin_up  +  (M_oe)*nrm*epsbar*(M_eo) * spinin_dn
  #endif
  
  
  
  
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
  
  		// xchange
  		#ifndef HOPPING_DEBUG
  		  xchange_field_wrapper(dev_spin_eo3_up, 0);
  		#endif
  
  // hopping matrix
  #ifdef USETEXTURE
    bind_texture_spin(dev_spin_eo3_up,1);
  #endif
  
    #ifndef HOPPING_DEBUG
    	#ifndef ALTERNATE_HOPPING_MATRIX
    	  dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo3_up, dev_spin_eo1_up, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0);	// dev_spin_eo1_up = (M_eo) * dev_spin_eo3_up
    	#else
    	  dev_Hopping_Matrix_alternate<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo3_up, dev_spin_eo1_up, dev_g_iup, dev_g_idn, dev_g_eo2lexic, dev_g_lexic2eosub, 0);
    	#endif
    #else
      Hopping_Matrix_wrapper(0, dev_spin_eo1_up, dev_spin_eo3_up);
    #endif
  
  #ifdef USETEXTURE
    unbind_texture_spin(1);
  #endif
  
  
  		// xchange
  		#ifndef HOPPING_DEBUG
  		  xchange_field_wrapper(dev_spin_eo3_dn, 0);
  		#endif
  
  // hopping matrix
  #ifdef USETEXTURE
    bind_texture_spin(dev_spin_eo3_dn,1);
  #endif
  
    #ifndef HOPPING_DEBUG
    	#ifndef ALTERNATE_HOPPING_MATRIX
    	  dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo3_dn, dev_spin_eo1_dn, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0);	// dev_spin_eo1_dn = (M_eo) * dev_spin_eo3_dn
    	#else
    	  dev_Hopping_Matrix_alternate<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo3_dn, dev_spin_eo1_dn, dev_g_iup, dev_g_idn, dev_g_eo2lexic, dev_g_lexic2eosub, 0);
    	#endif
    #else
      Hopping_Matrix_wrapper(0, dev_spin_eo1_dn, dev_spin_eo3_dn);
    #endif
  
  #ifdef USETEXTURE
    unbind_texture_spin(1);
  #endif
  
  
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
  
  
  		// xchange
  		#ifndef HOPPING_DEBUG
  		  xchange_field_wrapper(dev_spin_eo2_up, 1);
  		#endif
  
  // hopping matrix
  #ifdef USETEXTURE
    bind_texture_spin(dev_spin_eo2_up,1);									// remember: this is ((M_oe) (Mee^-1) (M_eo)) * (dev_spin_eo3_up, dev_spin_eo3_dn):
  #endif
  
    #ifndef HOPPING_DEBUG
    	#ifndef ALTERNATE_HOPPING_MATRIX
    	  dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo2_up, dev_spin_eo1_up, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1);
    	#else
    	  dev_Hopping_Matrix_alternate<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo2_up, dev_spin_eo1_up, dev_g_iup, dev_g_idn, dev_g_eo2lexic, dev_g_lexic2eosub, 1);
    	#endif
    #else
      Hopping_Matrix_wrapper(1, dev_spin_eo1_up, dev_spin_eo2_up);
    #endif
  
  #ifdef USETEXTURE
    unbind_texture_spin(1);							// dev_spin_eo1_up  =  (M_oe) * dev_spin_eo2_up  =  (M_oe)*nrm*(1-imubar)*(M_eo) * dev_spin_eo3_up  +  (M_oe)*nrm*epsbar*(M_eo) * dev_spin_eo3_dn
  #endif
  
  
  		// xchange
  		#ifndef HOPPING_DEBUG
  		  xchange_field_wrapper(dev_spin_eo2_dn, 1);
  		#endif
  
  // hopping matrix
  #ifdef USETEXTURE
    bind_texture_spin(dev_spin_eo2_dn,1);
  #endif
  
    #ifndef HOPPING_DEBUG
    	#ifndef ALTERNATE_HOPPING_MATRIX
    	  dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo2_dn, dev_spin_eo1_dn, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1);
    	#else
    	  dev_Hopping_Matrix_alternate<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo2_dn, dev_spin_eo1_dn, dev_g_iup, dev_g_idn, dev_g_eo2lexic, dev_g_lexic2eosub, 1);
    	#endif
    #else
      Hopping_Matrix_wrapper(1, dev_spin_eo1_dn, dev_spin_eo2_dn);
    #endif
  
  #ifdef USETEXTURE
    unbind_texture_spin(1);							// dev_spin_eo1_dn  =  (M_oe) * dev_spin_eo2_dn  =  (M_oe)*nrm*(1+imubar)*(M_eo) * dev_spin_eo3_dn  +  (M_oe)*nrm*epsbar*(M_eo) * dev_spin_eo3_up
  #endif
  
  
  
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
  
}//matrix_multiplication32_mpi()







extern "C" void benchmark_eo_nd_mpi (spinor * Q_up, spinor * Q_dn, int N) {

  
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  //													//
  // total FLOPS  =  (#iterations) * (FLOPS/matrix application) * (#lattice sites)			//
  //													//
  //													//
  // FLOPS per lattice site and application of the function,						//
  // count the floating point op's on device:								//
  //													//
  // dev_Hopping_Matrix	          = 4136								//
  // dev_mul_one_pm_imubar_gamma5 = 120									//
  // dev_gamma5                   = 12									//
  //													//
  // cublasSaxpy                  = 24*2 = 48								//
  // cublasSscal                  = 24*1 = 24								//
  //													//
  //													//
  // (FLOPS/matrix application)  =  2 * (4*4136 + 4*120 + 6*48 + 2*24 + 2*12)  =  2 * 17384  =  34768	//
  //													//
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  
  // timing
  double singleTimeElapsed;
  double maxTimeElapsed;
  double startBenchmark;
  double stopBenchmark;
  
  // counter
  int i;
  
  // flop counting
  double realFlopsPerApp = 34768.0;
  double effectiveFlopsPerApp = 23984.0;
  
  double realDeviceFlops;
  double allRealDeviceFlops;
  double realFlops;
  
  double effectiveDeviceFlops;
  double allEffectiveDeviceFlops;
  double effectiveFlops;
  
  // CUDA errors
  cudaError_t cudaerr;
  cublasStatus cublasstatus;
  
  // size of a spinor
  size_t dev_spinsize = 6*(VOLUME+RAND)/2 * sizeof(dev_spinor);
  
  // formal parameters
  int staticsource = 0;		// 1: applies matrix every time on the same source
  				// 0: applies matrix consecutively ...
  
  
  // init_mixedsolve_eo_nd(g_gauge_field);		// only when externally called
  
  
  dev_spinor * A_up;
  dev_spinor * A_dn;
  dev_spinor * B_up;
  dev_spinor * B_dn;
  
  dev_spinor * C_up;
  dev_spinor * C_dn;
  
  cudaMalloc((void **) &A_up, dev_spinsize);
  cudaMalloc((void **) &A_dn, dev_spinsize);
  cudaMalloc((void **) &B_up, dev_spinsize);
  cudaMalloc((void **) &B_dn, dev_spinsize);
  
  		// debug	// CUDA
  		#ifdef CUDA_DEBUG
  		  CUDA_CHECK_NO_SUCCESS_MSG("CUDA error in benchmark_eo_nd(). Memory allocation of spinor fields failed.");
  		#endif
  
  
  /*
  #ifdef USETEXTURE
    bind_texture_gf(dev_gf);
  #endif
  */
  
  
  /*		// only when externally called
  //Initialize some stuff
  dev_complex h0,h1,h2,h3,mh0, mh1, mh2, mh3;
  
  h0.re = (float)ka0.re;    h0.im = -(float)ka0.im;
  h1.re = (float)ka1.re;    h1.im = -(float)ka1.im;
  h2.re = (float)ka2.re;    h2.im = -(float)ka2.im;
  h3.re = (float)ka3.re;    h3.im = -(float)ka3.im;
  
  mh0.re = -(float)ka0.re;    mh0.im = (float)ka0.im;
  mh1.re = -(float)ka1.re;    mh1.im = (float)ka1.im;
  mh2.re = -(float)ka2.re;    mh2.im = (float)ka2.im;
  mh3.re = -(float)ka3.re;    mh3.im = (float)ka3.im;
  
  // try using constant mem for kappas
  cudaMemcpyToSymbol("dev_k0c", &h0, sizeof(dev_complex)); 
  cudaMemcpyToSymbol("dev_k1c", &h1, sizeof(dev_complex)); 
  cudaMemcpyToSymbol("dev_k2c", &h2, sizeof(dev_complex)); 
  cudaMemcpyToSymbol("dev_k3c", &h3, sizeof(dev_complex));
  
  cudaMemcpyToSymbol("dev_mk0c", &mh0, sizeof(dev_complex)); 
  cudaMemcpyToSymbol("dev_mk1c", &mh1, sizeof(dev_complex)); 
  cudaMemcpyToSymbol("dev_mk2c", &mh2, sizeof(dev_complex)); 
  cudaMemcpyToSymbol("dev_mk3c", &mh3, sizeof(dev_complex));
  */
  
  
  
  int blocksize;		// auxiliary
  
  blocksize = BLOCKSIZE1;
  int blockdim1, griddim1;					// here:	dev_zero_spinor_field , dev_copy_spinor_field
  if ( (VOLUME/2) % blocksize == 0 ) {
    blockdim1 = blocksize;
    griddim1  = VOLUME/2/blocksize;
  }
  else {
    blockdim1 = blocksize;
    griddim1  = (int) ((VOLUME/2/blocksize) + 1);
  }
  
  blocksize = BLOCKSIZE2;
  int blockdim2, griddim2;					// passed:	dev_Hopping_Matrix
  if ( (VOLUME/2) % blocksize == 0 ) {
    blockdim2 = blocksize;
    griddim2  = VOLUME/2/blocksize;
  }
  else {
    blockdim2 = blocksize;
    griddim2  = (int) ((VOLUME/2/blocksize) + 1);
  }
  
  blocksize = BLOCKSIZE3;
  int blockdim3, griddim3;					// passed:	dev_mul_one_pm_imubar_gamma5
  if ( (VOLUME/2) % blocksize == 0 ) {
    blockdim3 = blocksize;
    griddim3  = VOLUME/2/blocksize;
  }
  else {
    blockdim3 = blocksize;
    griddim3  = (int) ((VOLUME/2/blocksize) + 1);
  }
  
  blocksize = BLOCKSIZE4;
  int blockdim4, griddim4;					// passed:	dev_gamma5
  if ( (VOLUME/2) % blocksize == 0 ) {
    blockdim4 = blocksize;
    griddim4  = VOLUME/2/blocksize;
  }
  else {
    blockdim4 = blocksize;
    griddim4  = (int) ((VOLUME/2/blocksize) + 1);
  }
  
  blocksize = BLOCKSIZE5;
  int blockdim5, griddim5;					// passed:	dev_copy_spinor_field
  if ( (VOLUME/2) % blocksize == 0 ) {
    blockdim5 = blocksize;
    griddim5  = VOLUME/2/blocksize;
  }
  else {
    blockdim5 = blocksize;
    griddim5  = (int) ((VOLUME/2/blocksize) + 1);
  }
  
  
  		//debug
  		if (g_proc_id == 0) printf("\nStarting a little BENCHMARK. benchmark_eo_nd_mpi().\n");
  
  
  
  
  /*		// only when externally called
  he_cg_init<<< 1, 1 >>> (dev_grid, (float) g_kappa, (float)(g_mu/(2.0*g_kappa)), h0, h1, h2, h3);
  
  		// debug	// kernel
  		#ifdef CUDA_DEBUG
  		  CUDA_KERNEL_CHECK_NO_SUCCESS_MSG("Kernel error in he_cg_init(). Couldn't initialize some stuff.");
  		#endif
  
  
  he_cg_init_nd_additional<<<1,1>>> (g_mubar, g_epsbar);
  
  		// debug	// kernel
  		#ifdef CUDA_DEBUG
  		  CUDA_KERNEL_CHECK_NO_SUCCESS_MSG("Kernel error in he_cg_init_nd_additional(). Couldn't initialize some stuff.");
  		#endif
  */
  
  
  
  		/*
  		// debug	// CUBLAS helper function
  		#ifdef CUDA_DEBUG
  		  CUBLAS_HELPER_CHECK_NO_SUCCESS_MSG(cublasInit(), "CUBLAS error in benchmark_eo_nd(). Couldn't initialize CUBLAS.");
  		#else
  		  cublasInit();
  		#endif
  		*/
  
  
  		// debug
  		if (g_proc_id == 0) printf("Applying the eo-preconditioned matrix %i times.\n", N);
  
  
  to_device_mpi(B_up, Q_up, h2d_spin_up, dev_spinsize, 0, (VOLUME+RAND)/2);
  to_device_mpi(B_dn, Q_dn, h2d_spin_dn, dev_spinsize, 0, (VOLUME+RAND)/2);
  
  
  // timer
  startBenchmark = MPI_Wtime();
  
  
  
  for (i = 0; i < N; i++) {
  
    #ifndef ASYNC
      matrix_multiplication32_mpi(A_up, A_dn,				// A = (matrix)*B
                                  B_up, B_dn,
                                  griddim2, blockdim2,
                                  griddim3, blockdim3,
                                  griddim4, blockdim4,
                                  griddim5, blockdim5);
    #else
      matrix_multiplication32_mpi_ASYNC(A_up, A_dn,				// A = (matrix)*B
                                        B_up, B_dn,
                                        griddim2, blockdim2,
                                        griddim3, blockdim3,
                                        griddim4, blockdim4,
                                        griddim5, blockdim5);
    #endif
    
    if (staticsource = 0) {
      // swaps A and B
      C_up = B_up;
      C_dn = B_dn;
      B_up = A_up;
      B_dn = A_dn;
      A_up = C_up;
      A_dn = C_dn;
    }
    else {
      // do nothing
    }
    
  }
  
  		// debug	// CUDA
  		#ifdef CUDA_DEBUG
  		  // CUDA_CHECK("CUDA error in matrix_muliplication32(). Applying the matrix on GPU failed.", "The matrix was applied on GPU.");
  		  CUDA_CHECK_NO_SUCCESS_MSG("CUDA error in matrix_muliplication32(). Applying the matrix on GPU failed.");
  		#endif
  
  
  
  // timer
  stopBenchmark = MPI_Wtime();
  
  
  singleTimeElapsed = stopBenchmark - startBenchmark;
  MPI_Allreduce(&singleTimeElapsed, &maxTimeElapsed, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  
  realDeviceFlops      = N * VOLUME/2 * realFlopsPerApp;
  MPI_Allreduce(&realDeviceFlops, &allRealDeviceFlops, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  realFlops            = allRealDeviceFlops / maxTimeElapsed / 1.0e9;
  
  effectiveDeviceFlops = N * VOLUME/2 * effectiveFlopsPerApp;
  MPI_Allreduce(&effectiveDeviceFlops, &allEffectiveDeviceFlops, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  effectiveFlops       = allEffectiveDeviceFlops / maxTimeElapsed / 1.0e9;
  
  
  if (g_proc_id == 0) {
    /*
    printf("\tTEST:\n");
    printf("\ttime:        %.2e sec\n", singleTimeElapsed);
    printf("\tflop's:      %.2e flops\n", realDeviceFlops);
    printf("\tperformance: %.2e Gflop/s\n\n", realDeviceFlops / singleTimeElapsed / 1.0e9);
  
    printf("\tREAL:\n");
    printf("\ttime:        %.2e sec\n", maxTimeElapsed);
    printf("\tflop's:      %.2e flops\n", allRealDeviceFlops);
    printf("\tperformance: %.2e Gflop/s\n\n", realFlops);
    */
    printf("\tEFFECTIVE:\n");
    printf("\ttime:        %.2e sec\n", maxTimeElapsed);
    printf("\tflop's:      %.2e flops\n", allEffectiveDeviceFlops);
    printf("\tperformance: %.2e Gflop/s\n\n", effectiveFlops);
    
    #ifdef ASYNC_TIMING
      cudaEventElapsedTime(&time_D2H_1, start_ALL, stop_D2H_1);
      cudaEventElapsedTime(&time_D2H_2, start_ALL, stop_D2H_2);
      cudaEventElapsedTime(&time_INT_0, start_ALL, stop_INT_0);
      cudaEventElapsedTime(&time_H2D_3, start_ALL, stop_H2D_3);
      cudaEventElapsedTime(&time_H2D_4, start_ALL, stop_H2D_4);
      cudaEventElapsedTime(&time_EXT_1, start_ALL, stop_EXT_1);
      cudaEventElapsedTime(&time_EXT_2, start_ALL, stop_EXT_2);
      cudaEventElapsedTime(&time_ALL, start_ALL, stop_ALL);
      mpiTime_sendrecv_1 = mpiTime_stop_sendrecv_1 - mpiTime_start_sendrecv_1;
      mpiTime_sendrecv_2 = mpiTime_stop_sendrecv_2 - mpiTime_start_sendrecv_2;
      
      printf("\tADDITIONAL:\n");
      printf("\ttime_D2H_1         = %.4e sec\n", time_D2H_1/1000);
      printf("\ttime_D2H_2         = %.4e sec\n", time_D2H_2/1000);
      printf("\ttime_INT_0         = %.4e sec\n", time_INT_0/1000);
      printf("\tmpiTime_sendrecv_1 = %.4e sec\n", mpiTime_sendrecv_1);
      printf("\tmpiTime_sendrecv_2 = %.4e sec\n", mpiTime_sendrecv_2);
      printf("\ttime_H2D_3 = %.4e sec\n", time_H2D_3/1000);
      printf("\ttime_H2D_4 = %.4e sec\n", time_H2D_4/1000);
      printf("\ttime_EXT_1 = %.4e sec\n", time_EXT_1/1000);
      printf("\ttime_EXT_2 = %.4e sec\n", time_EXT_2/1000);
      printf("\ttime_ALL   = %.4e sec\n", time_ALL/1000);
    #endif
  
  }
  
  
  cudaFree(A_up);
  cudaFree(A_dn);
  cudaFree(B_up);
  cudaFree(B_dn);
  
  // finalize_mixedsolve_eo_nd();		// only when externally called
  
  /*
  #ifdef USETEXTURE
    unbind_texture_gf();
  #endif
  
  
  		// debug	// CUBLAS helper function
  		#ifdef CUDA_DEBUG
  		  CUBLAS_HELPER_CHECK_NO_SUCCESS_MSG(cublasShutdown(), "CUBLAS error in benchmark_eo_nd(). Couldn't shut down CUBLAS.");
  		#else
  		  cublasShutdown();
  		#endif
  */
  
  
}//benchmark_eo_nd_mpi()






////////////////////////
// CONJUGATE GRADIENT //
////////////////////////

// for the odd field after even/odd-preconditioning
// single precision on GPU

int cg_eo_nd_mpi (dev_su3_2v * gf,
                  dev_spinor * P_up, dev_spinor * P_dn,
                  dev_spinor * Q_up, dev_spinor * Q_dn,
                  int max_iter,
                  int check_abs , int check_rel,
                  double eps_abs, double eps_rel       ) {
  
  // P_up/dn  can be used as auxiliary field to work on, as it is not later used (could be used as initial guess at the very start)
  // Q_up/dn  can be used as feedback, or if not, also as auxiliary field
  
  
  
  
  /////////////////////
  // LOCAL VARIABLES //
  /////////////////////
  
  // CUDA
  cudaError_t cudaerr;
  cublasStatus cublasstatus;
  
  // algorithm
  float rr_up;
  float rr_dn;
  float rr;
  float rr_old;
  float r0r0;
  
  float dAd_up;
  float dAd_dn;
  float dAd;
  
  float alpha;
  float beta;
  
  // (auxiliary) device fields
  dev_spinor *  r_up, *  r_dn,
             * Ad_up, * Ad_dn,
             *  x_up, *  x_dn,
             *  d_up, *  d_dn,
             * Ax_up, * Ax_dn;		// for recalculating the residue
  
  // counting
  int j;				// iteration counter
  
  // formal parameters
  size_t dev_spinsize_int =  6*VOLUME/2*sizeof(dev_spinor);
  int N_sites_int         =    VOLUME/2;
  int N_floats_int        = 24*VOLUME/2;// (single precision) CUBLAS functions get the number of floats as input
  
  size_t dev_spinsize_ext =  6*(VOLUME+RAND)/2*sizeof(dev_spinor);
  int N_sites_ext         =    (VOLUME+RAND)/2;
  int N_floats_ext        = 24*(VOLUME+RAND)/2;
  
  // algorithm control parameters
  int N_recalc_res = 10;		// recalculate residue r(k+1) = b - A*x(k+1) each N_recalc_res iteration
  
  
  
  
  /////////////////////////////////////////////
  // CUDA block- and gridsize specifications //
  /////////////////////////////////////////////
  
  // int gridsize;		// auxiliary
  int blocksize;		// auxiliary
  
  blocksize = BLOCKSIZE1;
  int blockdim1, griddim1;					// here:	dev_zero_spinor_field , dev_copy_spinor_field
  if ( (VOLUME/2) % blocksize == 0 ) {
    blockdim1 = blocksize;
    griddim1  = VOLUME/2/blocksize;
  }
  else {
    blockdim1 = blocksize;
    griddim1  = (int) ((VOLUME/2/blocksize) + 1);
  }
  
  blocksize = BLOCKSIZE2;
  int blockdim2, griddim2;					// passed:	dev_Hopping_Matrix
  if ( (VOLUME/2) % blocksize == 0 ) {
    blockdim2 = blocksize;
    griddim2  = VOLUME/2/blocksize;
  }
  else {
    blockdim2 = blocksize;
    griddim2  = (int) ((VOLUME/2/blocksize) + 1);
  }
  
  blocksize = BLOCKSIZE3;
  int blockdim3, griddim3;					// passed:	dev_mul_one_pm_imubar_gamma5
  if ( (VOLUME/2) % blocksize == 0 ) {
    blockdim3 = blocksize;
    griddim3  = VOLUME/2/blocksize;
  }
  else {
    blockdim3 = blocksize;
    griddim3  = (int) ((VOLUME/2/blocksize) + 1);
  }
  
  blocksize = BLOCKSIZE4;
  int blockdim4, griddim4;					// passed:	dev_gamma5
  if ( (VOLUME/2) % blocksize == 0 ) {
    blockdim4 = blocksize;
    griddim4  = VOLUME/2/blocksize;
  }
  else {
    blockdim4 = blocksize;
    griddim4  = (int) ((VOLUME/2/blocksize) + 1);
  }
  
  blocksize = BLOCKSIZE5;
  int blockdim5, griddim5;					// passed:	dev_copy_spinor_field
  if ( (VOLUME/2) % blocksize == 0 ) {
    blockdim5 = blocksize;
    griddim5  = VOLUME/2/blocksize;
  }
  else {
    blockdim5 = blocksize;
    griddim5  = (int) ((VOLUME/2/blocksize) + 1);
  }
  		/*
  		// debug
  		printf("griddim1 = %i, blockdim1 = %i\n", griddim1, blockdim1);
  		printf("griddim2 = %i, blockdim2 = %i\n", griddim2, blockdim2);
  		printf("griddim3 = %i, blockdim3 = %i\n", griddim3, blockdim3);
  		printf("griddim4 = %i, blockdim4 = %i\n", griddim4, blockdim4);
  		printf("griddim5 = %i, blockdim5 = %i\n", griddim5, blockdim5);
  		*/
  
  
  
  
  /////////////////
  // ASSIGNMENTS //
  /////////////////
  
  x_up  = P_up;							// can use the output spinors also as auxiliary fields
  x_dn  = P_dn;							//	saves copying the output spinor field
  /*
  r_up  = Q_up;							// could possibly be done if  Q_up/dn  is not used as feedback
  r_dn  = Q_dn;							//	would save one field and one copying the field
  */
  r_up  = dev_spin1_up;						// use these pointers to the allocated space on device memory (allocated by init_mixedsolve_eo_nd)
  r_dn  = dev_spin1_dn;
  d_up  = dev_spin2_up;
  d_dn  = dev_spin2_dn;
  Ad_up = dev_spin3_up;
  Ad_dn = dev_spin3_dn;
  Ax_up = Ad_up;						// works as long as no initial guess vector x(0) is passed to cg_eo_nd()
  Ax_dn = Ad_dn;
  
  
  
  
  /////////////////////
  // INITIALIZATIONS //
  /////////////////////
  
  /*		// relocated to mixedsolve_eo_nd(), before here were:
  		// Initialize some stuff ...
  		// try using constant mem for kappas ...
  */
  
  /*
  // bind texture gf
  #ifdef USETEXTURE						// needed for subfunctions of dev_Hopping_Matrix(...)
    bind_texture_gf(gf);					//	e.g. dev_reconstructgf_2vtexref(...), dev_reconstructgf_8texref(...), 
  #endif							//	in general in functions  dev_reconstructgf[...]  with  "tex1Dfetch(gf_tex[...]"
  
  		// debug	// CUDA
    		#ifdef CUDA_DEBUG
    		  CUDA_CHECK("CUDA error in bind_texture_gf(). Binding GF to texture failed.", "GF bound to texture.");
    		#endif
  */
  
  /*		// relocated to mixedsolve_eo_nd(), before here were:
  he_cg_init<<< 1, 1 >>> (dev_grid, (REAL) g_kappa, (REAL)(g_mu/(2.0*g_kappa)), h0, h1, h2, h3);
  
  		// debug	// kernel
  		#ifdef CUDA_DEBUG
  		  CUDA_KERNEL_CHECK("Kernel error in he_cg_init(). Couldn't initialize some stuff.", "he_cg_init() succeeded.");
  		#endif
  
  
  he_cg_init_nd_additional<<<1,1>>> (g_mubar, g_epsbar);
  
  		// debug	// kernel
  		#ifdef CUDA_DEBUG
  		  CUDA_KERNEL_CHECK("Kernel error in he_cg_init_nd_additional(). Couldn't initialize some stuff.", "he_cg_init_nd_additional() succeeded.");
  		#endif
  */
  
  /*
  // cublasInit();			// init CUBLAS
  
  		// debug	// CUBLAS helper function
  		#ifdef CUDA_DEBUG
  		  CUBLAS_HELPER_CHECK(cublasInit(), "CUBLAS error in cublasInit(). Couldn't initialize CUBLAS.", "CUBLAS initialized.");
  		#else
  		  cublasInit();
  		#endif
  */
  
  
  
  
  ///////////////
  // ALGORITHM //
  ///////////////
  
  
  // initialize x(0) = 0	// will be added up
  dev_zero_spinor_field<<<griddim1, blockdim1>>>(x_up);
  dev_zero_spinor_field<<<griddim1, blockdim1>>>(x_dn);
  
  
  // r(0) = b - A*x(0) = b
  dev_copy_spinor_field<<<griddim1, blockdim1>>>(Q_up, r_up);
  dev_copy_spinor_field<<<griddim1, blockdim1>>>(Q_dn, r_dn);
  
  
  // d(0) = r(0)
  dev_copy_spinor_field<<<griddim1, blockdim1>>>(r_up, d_up);
  dev_copy_spinor_field<<<griddim1, blockdim1>>>(r_dn, d_dn);
  
  
  		// debug	// kernel
  		#ifdef CUDA_DEBUG
  		  CUDA_KERNEL_CHECK("Kernel error in cg_eo_nd(). Initializing spinor fields on device failed.", "Spinor fields initialized on device.");
  		#endif
  
  
  
  
  // rr = (r_up)^2 + (r_dn)^2
  rr_up = cublasSdot_wrapper(N_floats_int, (float *) r_up, 1, (float *) r_up, 1);
  rr_dn = cublasSdot_wrapper(N_floats_int, (float *) r_dn, 1, (float *) r_dn, 1);
  rr    = rr_up + rr_dn;
  
  		// benchmark
  		#ifdef GPU_BENCHMARK
  		  flopcount(device_flops, 2*2);
  		  // flopcount(device_flops, 2*2*N_floats);
  		#endif
  
  
  r0r0   = rr;	// for relative precision 
  rr_old = rr;	// for the first iteration
  
  
  
  
  //////////
  // LOOP //
  //////////
  
  		// debug
    		if (g_cart_id == 0) printf("\nEntering inner loop.\n");
  
  		// debug	// CUBLAS core function
		#ifdef CUDA_DEBUG
		  // CUBLAS_CORE_CHECK("CUBLAS error in cg_eo_nd(). Calculating initial residue failed.", "Initial inner residue calculated.");
		  CUBLAS_CORE_CHECK_NO_SUCCESS_MSG("CUBLAS error in cg_eo_nd(). Calculating initial residue failed.");
		#endif
  
  		// debug
  		if (g_cart_id == 0) printf("Initial inner residue: %.6e\n", r0r0);
  
  
  for (j = 0; j < max_iter; j++) {
    
    
    #ifndef MATRIX_DEBUG
    
      // A*d(k)
      #ifndef ASYNC
        matrix_multiplication32_mpi(Ad_up, Ad_dn,										// normally:  matrix_multiplication32_mpi()
                                     d_up,  d_dn,										// debugging: matrix_mpi_debug1/2/3/4()
                                    griddim2, blockdim2,
                                    griddim3, blockdim3,
                                    griddim4, blockdim4,
                                    griddim5, blockdim5);
      #else															// tries to overlap computation and communication
        matrix_multiplication32_mpi_ASYNC(Ad_up, Ad_dn,
                                           d_up,  d_dn,
                                          griddim2, blockdim2,
                                          griddim3, blockdim3,
                                          griddim4, blockdim4,
                                          griddim5, blockdim5);
      #endif
    
  		// debug	// CUDA		// also other stuff ?!
  		#ifdef CUDA_DEBUG
  		  // CUDA_CHECK("CUDA error in matrix_muliplication32(). Applying the matrix on GPU failed.", "The matrix was applied on GPU.");
  		  CUDA_CHECK_NO_SUCCESS_MSG("CUDA error in matrix_muliplication32(). Applying the matrix on GPU failed.");
  		#endif
    		
    		// benchmark
    		#ifdef GPU_BENCHMARK
    		  flopcount(device_flops, 1448);
      		  // flopcount(device_flops, 1448*N_floats);
    		#endif
    
    #else
    
    		// debug	// apply the host matrix on trial
    		
    		// host/device interaction
    		to_host(g_chi_up_spinor_field[DUM_SOLVER+3], d_up, h2d_spin_up, dev_spinsize_ext);
    		to_host(g_chi_dn_spinor_field[DUM_SOLVER+3], d_dn, h2d_spin_dn, dev_spinsize_ext);
    		
    		// matrix multiplication
    		if (g_proc_id == 0) printf("This is Q_Qdagger_ND(). ");
    		Q_Qdagger_ND(g_chi_up_spinor_field[DUM_SOLVER+4], g_chi_dn_spinor_field[DUM_SOLVER+4],			// normally:  Q_Qdagger_ND()
    		             g_chi_up_spinor_field[DUM_SOLVER+3], g_chi_dn_spinor_field[DUM_SOLVER+3] );		// debugging: matrix_mpi_debug10()
    		
    		// host/device interaction
    		to_device(Ad_up, g_chi_up_spinor_field[DUM_SOLVER+4], h2d_spin_up, dev_spinsize_ext);
    		to_device(Ad_dn, g_chi_dn_spinor_field[DUM_SOLVER+4], h2d_spin_dn, dev_spinsize_ext);
    		
    		
    				// debug	// CUDA
  				#ifdef CUDA_DEBUG
  			  	// CUDA_CHECK("CUDA error in cg_eo_nd(). Applying the matrix on CPU failed.", "The matrix was applied on CPU.");
  			  	CUDA_CHECK_NO_SUCCESS_MSG("CUDA error in cg_eo_nd(). Applying the matrix on CPU failed.");
  				#endif
    
    #endif
    
    
    // alpha = r(k)*r(k) / d(k)*A*d(k)
    dAd_up = cublasSdot_wrapper(N_floats_int, (float *) d_up, 1, (float *) Ad_up, 1);
    dAd_dn = cublasSdot_wrapper(N_floats_int, (float *) d_dn, 1, (float *) Ad_dn, 1);
    dAd    = dAd_up + dAd_dn;
    
    		// debug	// is NaN ?
    		if isnan(dAd) {
    		  printf("Error in cg_eo_nd(). dAd is NaN.\n");
    		  exit(-1);
    		}
    
    alpha  = rr_old / dAd;	// rr_old is taken from the last iteration respectively
    
    		// benchmark
    		#ifdef GPU_BENCHMARK
    		  flopcount(device_flops, 2*2);
    		  // flopcount(device_flops, 2*2*N_floats);
    		#endif
    
    
    // x(k+1) = x(k) + alpha*d(k)
    cublasSaxpy(N_floats_int, alpha, (float *) d_up, 1, (float *) x_up, 1);
    cublasSaxpy(N_floats_int, alpha, (float *) d_dn, 1, (float *) x_dn, 1);
    
    		// benchmark
    		#ifdef GPU_BENCHMARK
    		  flopcount(device_flops, 2*2);
    		  // flopcount(device_flops, 2*2*N_floats);
    		#endif
    
    
    // r(k+1)
    if ( (j+1) % N_recalc_res != 0 ) {	// r(k+1) = r(k) - alpha*A*d(k)
      cublasSaxpy(N_floats_int, -1.0*alpha, (float *) Ad_up, 1, (float *) r_up, 1);
      cublasSaxpy(N_floats_int, -1.0*alpha, (float *) Ad_dn, 1, (float *) r_dn, 1);
      
      		// benchmark
      		#ifdef GPU_BENCHMARK
      		  flopcount(device_flops, 2*2);
      		  // flopcount(device_flops, 2*2*N_floats);
      		#endif
    }
    
    else {				// recalculate residue r(k+1) = b - A*x(k+1)
    					//	"feedback"
      		// debug
      		if (g_proc_id == 0) printf("Recalculating the inner residue.\n");
      
      
      #ifndef MATRIX_DEBUG
        // A*x(k+1)
        #ifndef ASYNC
          matrix_multiplication32_mpi(Ax_up, Ax_dn,									// normally:  matrix_multiplication32_mpi()
                                       x_up,  x_dn,									// debugging: matrix_mpi_debug1/2/3/4()
                                      griddim2, blockdim2,
                                      griddim3, blockdim3,
                                      griddim4, blockdim4,
                                      griddim5, blockdim5);
        #else
          matrix_multiplication32_mpi_ASYNC(Ax_up, Ax_dn,									// normally:  matrix_multiplication32_mpi()
                                             x_up,  x_dn,									// debugging: matrix_mpi_debug1/2/3/4()
                                            griddim2, blockdim2,
                                            griddim3, blockdim3,
                                            griddim4, blockdim4,
                                            griddim5, blockdim5);
        #endif
      #else
      		// debug	// apply the host matrix on trial
    		
    		// host/device interaction
    		to_host(g_chi_up_spinor_field[DUM_SOLVER+3], x_up, h2d_spin_up, dev_spinsize_ext);
    		to_host(g_chi_dn_spinor_field[DUM_SOLVER+3], x_dn, h2d_spin_dn, dev_spinsize_ext);
    		
    		// matrix multiplication
    		if (g_proc_id == 0) printf("This is Q_Qdagger_ND(). ");
    		Q_Qdagger_ND(g_chi_up_spinor_field[DUM_SOLVER+4], g_chi_dn_spinor_field[DUM_SOLVER+4],			// normally:  Q_Qdagger_ND()
    		             g_chi_up_spinor_field[DUM_SOLVER+3], g_chi_dn_spinor_field[DUM_SOLVER+3] );		// debugging: matrix_mpi_debug10()
    		
    		// host/device interaction
    		to_device(Ax_up, g_chi_up_spinor_field[DUM_SOLVER+4], h2d_spin_up, dev_spinsize_ext);
    		to_device(Ax_dn, g_chi_dn_spinor_field[DUM_SOLVER+4], h2d_spin_dn, dev_spinsize_ext);
    		
    		
    				// debug	// CUDA
  				#ifdef CUDA_DEBUG
  			  	// CUDA_CHECK("CUDA error in cg_eo_nd(). Applying the matrix on CPU failed.", "The matrix was applied on CPU.");
  			  	CUDA_CHECK_NO_SUCCESS_MSG("CUDA error in cg_eo_nd(). Applying the matrix on CPU failed.");
  				#endif
      #endif
      
      
      		// benchmark
      		#ifdef GPU_BENCHMARK
      		  flopcount(device_flops, 1448);
      		  // flopcount(device_flops, 1448*N_floats);
      		#endif
      
      // r(k+1) = b - A*x(k+1)
      cublasScopy(N_floats_int, (float *) Q_up, 1, (float *) r_up, 1);		// r_up = Q_up
      cublasScopy(N_floats_int, (float *) Q_dn, 1, (float *) r_dn, 1);		// r_dn = Q_dn
      cublasSaxpy(N_floats_int, -1.0, (float *) Ax_up, 1, (float *) r_up, 1);	// r_up = Q_up - Ax_up
      cublasSaxpy(N_floats_int, -1.0, (float *) Ax_dn, 1, (float *) r_dn, 1);	// r_dn = Q_dn - Ax_dn
      
      		// benchmark
      		#ifdef GPU_BENCHMARK
      		  flopcount(device_flops, 2*2);
      		  // flopcount(device_flops, 2*2*N_floats);
      		#endif
    
    }
        
    
    // r(k+1)*r(k+1)
    rr_up  = cublasSdot_wrapper(N_floats_int, (float *) r_up, 1, (float *) r_up, 1);
    rr_dn  = cublasSdot_wrapper(N_floats_int, (float *) r_dn, 1, (float *) r_dn, 1);
    rr     = rr_up + rr_dn;
    
		// debug	// CUBLAS core function
		#ifdef CUDA_DEBUG
		  CUBLAS_CORE_CHECK_NO_SUCCESS_MSG("CUBLAS error in cg_eo_nd(). CUBLAS function failed.");
		#endif
		
		//benchmark
    		#ifdef GPU_BENCHMARK
    		  flopcount(device_flops, 2*2);
    		  // flopcount(device_flops, 2*2*N_floats);
    		#endif
    
    
    		// debug
    		if (g_proc_id == 0) printf("inner iteration j = %i: rr = %.6e\n", j, rr);
    		
    		// debug	// is NaN ?
    		if isnan(rr) {
    		  printf("Error in cg_eo_nd(). Inner residue is NaN.\n");
    		  exit(-1);
    		}
    
    
    // aborting ?? // check wether precision is reached ...
    if ( (check_abs)&&(rr <= eps_abs) || (check_rel)&&(rr <= eps_rel*r0r0) ) {
      
      if (g_cart_id == 0) {
      
      		// debug
      		printf("Finished inner loop because of reached precision.\n");
        
        if ((check_rel)&&(rr <= eps_rel*r0r0)) {
      		// debug
      		printf("Reached relative inner solver precision of eps_rel = %.2e\n", eps_rel);
        }
        if ((check_abs)&&(rr <= eps_abs)) {
      		// debug
      		printf("Reached absolute inner solver precision of eps_abs = %.2e\n", eps_abs);
        }
        
      		//debug
      		printf("Final inner residue: %.6e\n", rr);
      	
      }
      
      /*
      #ifdef USETEXTURE
        unbind_texture_gf();
      #endif
      
      		// debug	// CUDA
  		#ifdef CUDA_DEBUG
  		  CUDA_CHECK("CUDA error in unbind_texture(). Unbindung the GF texture failed.", "GF texture unbound.");
  		#endif
      
      
      // cublasShutdown();			// ends CUBLAS
      
      		// debug	// CUBLAS helper function
  		#ifdef CUDA_DEBUG
  		  CUBLAS_HELPER_CHECK(cublasShutdown(), "CUBLAS error in cublasInit(). Couldn't shut down CUBLAS.", "CUBLAS is shutted down.");
  		#else
  		  cublasShutdown();
  		#endif
      */
      
      return(j+1);
    }
    
    
    // beta = r(k+1)*r(k+1) / r(k)*r(k)
    beta = rr / rr_old;
    
    
    rr_old = rr;  // for next iteration
    
    
    // d(k+1) = r(k+1) + beta*d(k)
    cublasSscal (N_floats_int, beta, (float *) d_up, 1);
    cublasSaxpy (N_floats_int, 1.0 , (float *) r_up, 1, (float *) d_up, 1);
    
    cublasSscal (N_floats_int, beta, (float *) d_dn, 1);
    cublasSaxpy (N_floats_int, 1.0 , (float *) r_dn, 1, (float *) d_dn, 1);
    
    		// debug	// CUBLAS core function
    		#ifdef CUDA_DEBUG
    		  CUBLAS_CORE_CHECK_NO_SUCCESS_MSG("CUBLAS error in cg_eo_nd(). Error in CUBLAS function.");
    		#endif
    		
    		// benchmark
    		#ifdef GPU_BENCHMARK
    		  flopcount(device_flops, 2*3);
    		  // flopcount(device_flops, 2*3*N_floats);
    		#endif
  
  
  }//LOOP
  
  
  		// debug
  		if (g_cart_id == 0) printf("Finished inner loop beacuse of maximal number of inner iterations.\n");
  		if (g_cart_id == 0) printf("Final inner residue: %.6e\n", rr);
  
  /*
  #ifdef USETEXTURE
    unbind_texture_gf();
  #endif
  
  		// debug	// CUDA
  		#ifdef CUDA_DEBUG
  		  CUDA_CHECK("CUDA error in unbind_texture(). Unbindung the GF texture failed.", "GF texture unbound.");
  		#endif
  
  
  // cublasShutdown();
  
  		// debug	// CUBLAS helper function
  		#ifdef CUDA_DEBUG
  		  CUBLAS_HELPER_CHECK(cublasShutdown(), "CUBLAS error in cublasInit(). Couldn't shut down CUBLAS.", "CUBLAS is shutted down.");
  		#else
  		  cublasShutdown();
  		#endif
  */
  
  return(j+1);
  
}//cg_eo_nd()








//////////////////
// MIXED SOLVER //
//////////////////

// iterative refinement, defect correction
// that function is to replace the call of  cg_her_nd()  in  invert_doublet_eo.c
// solves the odd part of the full eo and nd problem
//	more precisely we have to invert  Qhat(2x2)*Qhat(2x2)^dagger
//	multiplying by  Qhat(2x2)^dagger  is done in  invert_doublet_eo.c

extern "C" int mixedsolve_eo_nd_mpi (spinor * P_up, spinor * P_dn,
                                     spinor * Q_up, spinor * Q_dn,
                                     int max_iter, double eps_sq, int rel_prec) {
  
  
  // basically  P_up/dn  and  Q_up/dn  could be used as auxiliary fields
  //	P_up/dn  is the output field (and can be used as initial guess)
  //	Q_up/dn  is not used later in the calling  invert_doublet_eo.c
  //		 but will be used as feedback in r(k+1) = b - A*x(k+1)
  
  
  		// debug
  		if (g_proc_id == 0) {
  		  printf("\n\nmixedsolve_eo_nd_mpi():\n");
  		  
  		  printf("SOLVER PARAMETERS:\n");
  		  
  		  printf("outer:");
  		  printf("\tmaximal iterations: %i\n", max_iter);
  		  printf("\trelative check?:    %i\n", bool(rel_prec));
  		  printf("\tprecision:          %.8e\n", eps_sq);
  		  
  		  printf("inner:");
  		  printf("\tmaximal iterations: %i\n", max_innersolver_it);
  		  printf("\tabsolute check?:    %i\n", bool(innersolver_precision_check_abs));
  		  printf("\trelative check?:    %i\n", bool(innersolver_precision_check_rel));
  		  printf("\tabsolute precision: %.8e\n", innersolver_precision_abs);
  		  printf("\trelative precision: %.8e\n", innersolver_precision_rel);
  		}
  
  
  
  
  /////////////////////
  // LOCAL VARIABLES //
  /////////////////////
  
  // CUDA
  cudaError_t cudaerr;
  cublasStatus cublasstatus;
  
  // algorithm
  double rr_up;
  double rr_dn;
  double rr;
  double rr_old;
  double r0r0;
  double bb;
  
  // counting
  int i = 0;					// iteration counter
  int innercount;				// latest inner solver iterations
  int outercount = 0;				// total inner solver iterations
  double flops;
  #ifdef EFFECTIVE_BENCHMARK
    double effectiveflops;			// will used to count the "effective" flop's (from the algorithmic perspective)
    double hoppingflops = 1488.0;
    double matrixflops  = 2  *  (  2 * ( (2*hoppingflops+12+3) + (2*hoppingflops+3) + (12+2) + 12 )  );
    double allflops;				// flops added for all processes
  #endif
  
  // timing
  clock_t startouter, stopouter;
  clock_t startinner, stopinner;
  // double timeelapsed;
  clock_t innerclocks;
  clock_t totalinnerclocks = 0;
  clock_t totalouterclocks = 0;
  #ifdef EFFECTIVE_BENCHMARK
    double starteffective;
    double stopeffective;
    double singletime;				// time for each process = stopeffective - starteffective
    double maxtime;				// max. parallel process time
  #endif
  
  // (auxiliary) fields
  spinor *  r_up, *  r_dn,
         * Ad_up, * Ad_dn,
         *  x_up, *  x_dn,
         *  d_up, *  d_dn,
         * Ax_up, * Ax_dn;
  
  // formal parameters
  size_t dev_spinsize_int =  6*VOLUME/2*sizeof(dev_spinor);		// 24 floats per spinor per even lattice site
  int N_sites_int         =    VOLUME/2;				// Carsten's functions get the number of lattice points as input
  int N_floats_int        = 24*VOLUME/2;
  
  size_t dev_spinsize_ext =  6*(VOLUME+RAND)/2*sizeof(dev_spinor);
  int N_sites_ext         =    (VOLUME+RAND)/2;
  int N_floats_ext        = 24*(VOLUME+RAND)/2;
  
  // algorithm control parameters
  bool rbAx = true;						// choose how to calculate r(k+1)
  bool initial_guess = false;					// choose if initial guess
  
  
  
  
  //////////////////
  // INITIALIZING //
  //////////////////
  
  
  		//debug
  		if (g_cart_id == 0) printf("init_mixedsolve_eo_nd_mpi():\n");
  
  
  init_mixedsolve_eo_nd_mpi(g_gauge_field);			// initializes and allocates all quantities for the mixed solver
  								// more precise:
  								//	puts the gauge field on device as "2 rows" or "8 floats" per SU(3)-matrix
  								//	allocates memory for all spinor fields
  								//	puts the nn- and eoidx-fields on device memory
    
  
  		//debug
  		if (g_cart_id == 0) printf("mixedsolve_eo_nd_mpi():\n");
  
  
  // the following initializations are moved from cg_eo_nd():
  
  // Initialize some stuff
  dev_complex h0, h1, h2, h3, mh0, mh1, mh2, mh3;
  
  h0.re  =  (float) ka0.re;	h0.im  = -(float) ka0.im;	// ka{0-4} are defined in boundary.c
  h1.re  =  (float) ka1.re;	h1.im  = -(float) ka1.im;	// what is the meaning?
  h2.re  =  (float) ka2.re;	h2.im  = -(float) ka2.im;
  h3.re  =  (float) ka3.re;	h3.im  = -(float) ka3.im;
  
  mh0.re = -(float) ka0.re;	mh0.im =  (float) ka0.im;
  mh1.re = -(float) ka1.re;	mh1.im =  (float) ka1.im;
  mh2.re = -(float) ka2.re;	mh2.im =  (float) ka2.im;
  mh3.re = -(float) ka3.re;	mh3.im =  (float) ka3.im;
  /*
  // try using constant mem for kappas		// constant memory is cached!
  cudaMemcpyToSymbol("dev_k0c", &h0, sizeof(dev_complex));
  cudaMemcpyToSymbol("dev_k1c", &h1, sizeof(dev_complex));
  cudaMemcpyToSymbol("dev_k2c", &h2, sizeof(dev_complex));
  cudaMemcpyToSymbol("dev_k3c", &h3, sizeof(dev_complex));
  
  cudaMemcpyToSymbol("dev_mk0c", &mh0, sizeof(dev_complex));
  cudaMemcpyToSymbol("dev_mk1c", &mh1, sizeof(dev_complex));
  cudaMemcpyToSymbol("dev_mk2c", &mh2, sizeof(dev_complex));
  cudaMemcpyToSymbol("dev_mk3c", &mh3, sizeof(dev_complex));
  */
  
  
  // bind texture gf
  #ifdef USETEXTURE						// needed for subfunctions of dev_Hopping_Matrix(...)
    bind_texture_gf(dev_gf);					//	e.g. dev_reconstructgf_2vtexref(...), dev_reconstructgf_8texref(...), 
  #endif							//	in general in functions  dev_reconstructgf[...]  with  "tex1Dfetch(gf_tex[...]"
  
  		// debug	// CUDA
    		#ifdef CUDA_DEBUG
    		  CUDA_CHECK("CUDA error in bind_texture_gf(). Binding GF to texture failed.", "GF bound to texture.");
    		#endif
  
  
  he_cg_init<<< 1, 1 >>> (dev_grid, (float) g_kappa, (float)(g_mu/(2.0*g_kappa)), h0, h1, h2, h3);
  		// "he" = "host entry"
  		// BEWARE in dev_tm_dirac_kappa we need the true mu (not 2 kappa mu!)	// ??
  		
  		// dev_LX, dev_LY, dev_LZ, dev_T, dev_VOLUME  =  grid[5]  =  dev_grid[5]
  		//	dev_VOLUME  is necessary for many kernel functions as for instance  dev_gamma5()
  		// initializes  mu, kappa and twokappamu  on the device
  		// initializes the strange  dev_k{0-3}, dev_mk{0-3}  as derived from the  ka{0-3}  from boundary.c
  
  		// debug	// kernel
  		#ifdef CUDA_DEBUG
  		  CUDA_KERNEL_CHECK("Kernel error in he_cg_init(). Couldn't initialize some stuff.", "he_cg_init() succeeded.");
  		#endif
  
  		// debug	// check stuff on device
  		#ifdef STUFF_DEBUG
  		if (g_proc_id == 0) {
  		  printf("\tOn host:\n");
  		  printf("\tVOLUME = %i\n", VOLUME);							// checking VOLUME and RAND in the parallel case 
  		  printf("\tRAND   = %i\n", RAND);
  		  printf("\tVOLUME + RAND = %i\n",  VOLUME+RAND);
  		  printf("\tVOLUMEPLUSRAND = %i\n", VOLUMEPLUSRAND);
  		
  		  int host_check_LX, host_check_LY, host_check_LZ, host_check_T, host_check_VOLUME;
  		  cudaMemcpyFromSymbol(&host_check_LX, dev_LX, sizeof(int));
  		  cudaMemcpyFromSymbol(&host_check_LY, dev_LY, sizeof(int));
  		  cudaMemcpyFromSymbol(&host_check_LZ, dev_LZ, sizeof(int));
  		  cudaMemcpyFromSymbol(&host_check_T, dev_T, sizeof(int));
  		  cudaMemcpyFromSymbol(&host_check_VOLUME, dev_VOLUME, sizeof(int));
  		  // printf("\teven_odd_flag = %i\n", even_odd_flag);
  		  printf("\tOn device:\n");
  		  printf("\tdev_LX = %i\n", host_check_LX);
  		  printf("\tdev_LY = %i\n", host_check_LY);
  		  printf("\tdev_LZ = %i\n", host_check_LZ);
  		  printf("\tdev_T = %i\n", host_check_T);
  		  printf("\tdev_VOLUME = %i/2 ?!= %i\n", host_check_LX*host_check_LY*host_check_LZ*host_check_T, host_check_VOLUME);
  		  
  		  float host_check_mu, host_check_kappa, host_check_twokappamu;
  		  cudaMemcpyFromSymbol(&host_check_mu, mu, sizeof(float));
  		  cudaMemcpyFromSymbol(&host_check_kappa, kappa, sizeof(float));
  		  cudaMemcpyFromSymbol(&host_check_twokappamu, twokappamu, sizeof(float));
  		  // printf("\tOn device:\n");
  		  // printf("\tmu = %f\n", host_check_mu);		// not needed for the nd case
  		  printf("\tkappa = %f\n", host_check_kappa);
  		  // printf("\ttwokappamu = %f\n", host_check_twokappamu);
  		  }
  		#endif
  
  
  he_cg_init_nd_additional<<<1,1>>> (g_mubar, g_epsbar);
  
  		// debug	// kernel
  		#ifdef CUDA_DEBUG
  		  CUDA_KERNEL_CHECK("Kernel error in he_cg_init_nd_additional(). Couldn't initialize some stuff.", "he_cg_init_nd_additional() succeeded.");
  		#endif
  
  		// debug	// check mubar and epsbar on host and device
  		#ifdef STUFF_DEBUG
  		if (g_proc_id == 0) {
  		  // printf("\tOn host:\n");
  		  // printf("\tg_mubar = %f\n", g_mubar);
  		  // printf("\tg_epsbar = %f\n", g_epsbar);
  		  
  		  float host_check_mubar, host_check_epsbar;
  		  cudaMemcpyFromSymbol(&host_check_mubar, mubar, sizeof(float));
  		  cudaMemcpyFromSymbol(&host_check_epsbar, epsbar, sizeof(float));
  		  printf("\tOn device:\n");
  		  printf("\tmubar = %f\n", host_check_mubar);
  		  printf("\tepsbar = %f\n", host_check_epsbar);
  		}
  		#endif
  
  
  he_cg_init_nd_additional_mpi<<<1,1>>>(VOLUMEPLUSRAND, RAND, g_cart_id, g_nproc);
  
  		// debug	// kernel
  		#ifdef CUDA_DEBUG
  		  CUDA_KERNEL_CHECK("Kernel error in he_cg_init_nd_additional_mpi(). Couldn't initialize some stuff.", "he_cg_init_nd_additional_mpi() succeeded.");
  		#endif
  
  		// debug	// check dev_VOLUMEPLUSRAND and dev_RAND on device
  		#ifdef STUFF_DEBUG
  		if (g_proc_id == 0) {
  		  int host_check_VOLUMEPLUSRAND, host_check_RAND;
  		  cudaMemcpyFromSymbol(&host_check_VOLUMEPLUSRAND, dev_VOLUMEPLUSRAND, sizeof(int));
  		  cudaMemcpyFromSymbol(&host_check_RAND, dev_RAND, sizeof(int));
  		  printf("\tOn device:\n");
  		  printf("\tdev_VOLUMEPLUSRAND = %i\n", host_check_VOLUMEPLUSRAND);
  		  printf("\tdev_RAND = %i\n", host_check_RAND);
  		}
  		#endif
  
  
  
  
  /*		// necessary ??
  // cublasInit();
  
  		// debug	// CUBLAS helper function
  		#ifdef CUDA_DEBUG
  		  CUBLAS_HELPER_CHECK(cublasInit(), "CUBLAS error in cublasInit(). Couldn't initialize CUBLAS.", "CUBLAS initialized.");
  		#else
  		  cublasInit();
  		#endif
  */
  
  
  
  #ifdef OPERATOR_BENCHMARK
    benchmark_eo_nd_mpi(Q_up, Q_dn, OPERATOR_BENCHMARK);
  #endif
  
  
  
  
  /////////////////
  // ASSIGNMENTS //
  /////////////////
  
  
  x_up = P_up;							// can use the output spinors also as auxiliary fields
  x_dn = P_dn;							//	can use as initial guess at the same time
  
  
  #ifndef CG_DEBUG
  
    r_up  = g_chi_up_spinor_field[DUM_SOLVER];			// use the pre-allocated memory on host memory
    r_dn  = g_chi_dn_spinor_field[DUM_SOLVER];			// allocated by  init_chi_spinor_field.c  and  invert_doublet.c  !?
    d_up  = g_chi_up_spinor_field[DUM_SOLVER+1];		// the fields  g_chi_up/dn_spinor_field[DUM_SOLVER{ , +1, ... , +5}]  are used in  cg_her_nd()
    d_dn  = g_chi_dn_spinor_field[DUM_SOLVER+1];
    Ad_up = g_chi_up_spinor_field[DUM_SOLVER+2];
    Ad_dn = g_chi_dn_spinor_field[DUM_SOLVER+2];
    Ax_up = Ad_up;
    Ax_dn = Ad_dn;
  		// debug
  		if (g_cart_id == 0) printf("Now using the fields g_chi_up/dn_spinor_field[DUM_SOLVER{ , +1, +2}] in the mixedsolve_eo_nd().\n");
  
  #else
  
  		r_up  = (spinor *) malloc(24*N_sites_int*sizeof(double));		// if using cg_her_nd() as the CG, we cannot use the g_chi_up/dn-fields at the same time
  		r_dn  = (spinor *) malloc(24*N_sites_int*sizeof(double));
  		d_up  = (spinor *) malloc(24*N_sites_int*sizeof(double));		// N_sites_int because only fields on which Hopping_Matrix() is directly applied
  		d_dn  = (spinor *) malloc(24*N_sites_int*sizeof(double));		//	have to have the boundaries !!
  		Ad_up = (spinor *) malloc(24*N_sites_int*sizeof(double));
  		Ad_dn = (spinor *) malloc(24*N_sites_int*sizeof(double));
  		Ax_up = Ad_up;
  		Ax_dn = Ad_dn;
  				// debug
  				if (g_cart_id == 0) printf("Now allocating new host space for the fields in mixedsolve_eo_nd().\n");
  
  #endif
  
  
  		// benchmark
  		#ifdef GPU_BENCHMARK
  		  device_flops = 0;
  		#endif
  		
  		#ifdef CPU_BENCHMARK
  		  host_flops = 0;
  		#endif
  
  
  
  
  ///////////////
  // ALGORITHM //
  ///////////////
  
  // timer
  startouter = clock();
  
  #ifdef EFFECTIVE_BENCHMARK
    starteffective = MPI_Wtime();
  #endif
  
  
  // r(0)
  if (!initial_guess) {		// r(0) = b = Q	// for x(0) = 0
    assign(r_up, Q_up, N_sites_int);
    assign(r_dn, Q_dn, N_sites_int);
    if (g_cart_id == 0) printf("x(0) = 0\n");
  }
  else {			// r(0) = b - A*x(0) = Q - A*P
    bb = square_norm(P_up, N_sites_int, 1) + square_norm(P_dn, N_sites_int, 1);			// NOTICE: for the parallel case we have to set "1" in the interface of square_norm()
    		// benchmark
    		#ifdef CPU_BENCHMARK
    		  flopcount(host_flops, 2*2);
    		  // flopcount(host_flops, 2*2*N_floats);
    		#endif
    if (g_cart_id == 0) printf("bb = %.10e\n", bb);
    if (bb == 0) {
      assign(r_up, Q_up, N_sites_int);
      assign(r_dn, Q_dn, N_sites_int);
      if (g_cart_id == 0) printf("x(0) = 0\n");
    }
    else {
      Q_Qdagger_ND(Ax_up, Ax_dn, P_up, P_dn);
      diff(r_up, Q_up, Ax_up, N_sites_int);
      diff(r_dn, Q_dn, Ax_dn, N_sites_int);
      		// benchmark
      		#ifdef CPU_BENCHMARK
      		  flopcount(host_flops, 2*2*(55+2+2+1+55) + 2);
      		  // flopcount(host_flops, 2*2*(55+2+2+1+55)*N_floats + 2*N_floats);
      		#endif
      if (g_cart_id == 0) printf("x(0) != 0\n");
    }
  }
  
  
  // rr = (r_up)^2 + (r_dn)^2
  rr_up = square_norm(r_up, N_sites_int, 1);
  rr_dn = square_norm(r_dn, N_sites_int, 1);
  rr = rr_up + rr_dn;
  
  		// benchmark
  		#ifdef CPU_BENCHMARK
  		  flopcount(host_flops, 2*2);
  		  // flopcount(host_flops, 2*2*N_floats);
  		#endif
  
  
  r0r0   = rr; // for relative precision
  rr_old = rr; // for the first iteration
  
  		// debug
  		if (g_cart_id == 0) printf("Initial outer residue: %.10e\n", rr_old);
  
  
  // set to zero	// x_up, x_dn  will be added up		// as  x_up/dn = P_up/dn  up to here  P_up/dn  was not changed
  zero_spinor_field(x_up, N_sites_int);
  zero_spinor_field(x_dn, N_sites_int);
  
  
  
  
  ////////////////
  // OUTER LOOP //
  ////////////////
  
  		// debug
    		if (g_cart_id == 0) printf("\nEntering outer loop.");
  
  
  do {		// for (i = 0; i < max_iter; i++) {
    
    i++;
  
    		// debug
    		if (g_cart_id == 0) printf("\nouter iteration i = %i\n", i);
    
    
    
    
    #ifndef CG_DEBUG
    
      // host/device interaction
      to_device(dev_spinin_up, r_up, h2d_spin_up, dev_spinsize_ext);
      to_device(dev_spinin_dn, r_dn, h2d_spin_dn, dev_spinsize_ext);
    
    		// debug	// CUDA
    		#ifdef CUDA_DEBUG
    		  // CUDA_CHECK("CUDA error in mixedsolve_eo_nd(). Host to device interaction failed.", "Fields copied to device.");
    		  CUDA_CHECK_NO_SUCCESS_MSG("CUDA error in mixedsolve_eo_nd(). Host to device interaction failed.");
    		#endif
    
    
    
    
      ////////////////////////////////////
      // INNER LOOP, CONJUGATE GRADIENT //
      ////////////////////////////////////
      
      // timer
      startinner = clock();
      
      		// debug
      		if (g_cart_id == 0) printf("cg_eo_nd():\n");
      
      // solves A*p(k+1) = r(k)
      //        A*p(0)   = r(0) = b
      innercount = cg_eo_nd_mpi(dev_gf,
                                dev_spinout_up, dev_spinout_dn,
                                dev_spinin_up , dev_spinin_dn,
                                max_innersolver_it,
                                innersolver_precision_check_abs, innersolver_precision_check_rel,
                                innersolver_precision_abs      , innersolver_precision_rel      );
      
      outercount = outercount + innercount;
      
      // timer
      stopinner = clock();
      innerclocks = stopinner-startinner;
      totalinnerclocks = totalinnerclocks + innerclocks;
      
    		// debug
    		if (g_cart_id == 0) printf("Inner solver done in: %.4e sec\n", double(innerclocks) / double(CLOCKS_PER_SEC));
    
    
      // host/device interaction
      to_host(d_up, dev_spinout_up, h2d_spin_up, dev_spinsize_ext);
      to_host(d_dn, dev_spinout_dn, h2d_spin_dn, dev_spinsize_ext);
    
    		// debug	// CUDA
    		#ifdef CUDA_DEBUG
    		  // CUDA_CHECK("CUDA error in mixedsolve_eo_nd(). Device to host interaction failed.", "Fields copied back to device.");
    		  CUDA_CHECK_NO_SUCCESS_MSG("CUDA error in mixedsolve_eo_nd(). Device to host interaction failed.");
    		#endif
    
    
    #else	// CG_DEBUG
    
    
    				// debug
    				if (g_cart_id == 0) printf("cg_her_nd():\n");
    		
    		innercount = cg_her_nd(d_up, d_dn, r_up, r_dn,
				       1000, eps_sq/2, 0,
				       VOLUME/2, &Q_Qdagger_ND, 0, 1000);
    		
    		outercount = outercount + innercount;
    		
    				// debug
    				if (g_cart_id == 0) printf("cg_her_nd() on host was used for debugging purposes.\n");
    
    
    #endif	// CG_DEBUG
    
    
    
    
    		// debug
    		if (g_cart_id == 0) printf("mixedsolve_eo_nd():\n");
    
    
    // x(k+1) = x(k) + d(k+1)
    add(x_up, x_up, d_up, N_sites_int);
    add(x_dn, x_dn, d_dn, N_sites_int);
    
    		// benchmark
  		#ifdef CPU_BENCHMARK
  		  flopcount(host_flops, 2);
  		  // flopcount(host_flops, 2*N_floats);
  		#endif
    
    
    // r(k+1)
    if (rbAx) {				// r(k+1) = b - A*x(k+1)
      // A*x(k+1)
      Q_Qdagger_ND(Ax_up, Ax_dn, x_up, x_dn);
      		// debug
      		if (g_cart_id == 0) printf("The matrix was applied on CPU in double precision. r = b - Ax\n");
      diff(r_up, Q_up, Ax_up, N_sites_int);
      diff(r_dn, Q_dn, Ax_dn, N_sites_int);
    }
    else {				// r(k+1) = r(k) - A*d(k+1)	// makes actually no sense ;)
      // A*d(k+1)
      Q_Qdagger_ND(Ad_up, Ad_dn, d_up, d_dn);
    		// debug
    		if (g_cart_id == 0) printf("The matrix was applied on CPU in double precision. r = r - Ad\n");
      // r(k+1) = r(k) - A*d(k+1)
      diff(r_up, r_up, Ad_up, N_sites_int);
      diff(r_dn, r_dn, Ad_dn, N_sites_int);
    }
    
    		// benchmark
  		#ifdef CPU_BENCHMARK
  		  flopcount(host_flops, 2*2*(55+2+2+1+55) + 2);
  		  // flopcount(host_flops, 2*2*(55+2+2+1+55)*N_floats + 2*N_floats);
  		#endif
    
    
    // rr = (rr_up)^2 + (r_dn)^2
    rr_up = square_norm(r_up, N_sites_int, 1);
    rr_dn = square_norm(r_dn, N_sites_int, 1);
    rr    = rr_up + rr_dn;
    
    		// debug
    		if (g_cart_id == 0) printf("Outer residue in the outer iteration i = %i after %i total inner iterations : %.10e\n", i, outercount, rr);
    		
    		// debug	// is NaN ?
    		if isnan(rr) {
    		  printf("Error in mixedsolve_eo_nd(). Outer residue is NaN.\n");
    		  exit(-1);
    		}
    		
    		// benchmark
  		#ifdef CPU_BENCHMARK
  		  flopcount(host_flops, 2*2);
  		  // flopcount(host_flops, 2*2*N_floats);
  		#endif
    
    
    // aborting ?? // check wether precision is reached ...
    if ( ((rr <= eps_sq) && (rel_prec == 0))  ||  ((rr <= eps_sq*r0r0) && (rel_prec == 1)) ) {
      
      // timer
      stopouter = clock();
      totalouterclocks = stopouter-startouter - totalinnerclocks;
      
      #ifdef EFFECTIVE_BENCHMARK
        stopeffective = MPI_Wtime();
      #endif
      
      		/*
      		// benchmark
  		#ifdef GPU_BENCHMARK2
  		  device_flops = 0;
  		  int help = ( 4 + outercount*(1448+5*4+6) + outercount/10*1448 ) * N_floats;
  		  flopcount(device_flops, help);			// N_recalcres = 10
  		#endif
  		*/
      
      		// debug
      		if (g_cart_id == 0) {
      		  printf("\nEO inversion done in mixed precision.\n");
      		  if (rel_prec == 0) printf("Finished outer loop because of reached absolute outer solver precision.\n");
      		  if (rel_prec == 1) printf("Finished outer loop because of reached relative outer solver precision.\n");
      		  printf("Total number of inner iterations: %i\n", outercount);
      		  printf("Total number of outer iterations: %i\n", i+1);
      		  printf("Squared residue: %.10e\n", rr); 
      		  printf("Outer solver done in: %.4e sec\n", double(stopouter-startouter) / double(CLOCKS_PER_SEC));
      		}
      		// benchmark
      		#ifdef EFFECTIVE_BENCHMARK
      		  // will now count the number of effective flops
      		  // effectiveflops  =  #(inner iterations)*(matrixflops+linalgflops)*VOLUME/2  +  #(outer iterations)*(matrixflops+linalgflops)*VOLUME/2
      		  // outer loop: linalg  =  flops for calculating  r(k+1) and x(k+1)
      		  // inner loop: linalg  =  flops for calculating  alpha, x(k+1), r(k+1), beta, d(k+1)
      		  singletime = double(stopeffective-starteffective);
      		  effectiveflops = outercount*(matrixflops + 2*2*2*24 + 2*2*24 + 2*2*24 + 2*2*2*24 + 2*2*24)*VOLUME/2   +   i*(matrixflops + 2*24 + 2*24)*VOLUME/2;
      		  MPI_Allreduce(&singletime, &maxtime, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      		  MPI_Allreduce(&effectiveflops, &allflops, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      		  if (g_proc_id == 0) printf("effective BENCHMARK:\n");
      		  if (g_proc_id == 0) printf("\ttotal mixed solver time:   %.2e sec\n", double(maxtime));
      		  if (g_proc_id == 0) printf("\tfloating point operations: %.2e flops\n", double(allflops));
      		  if (g_proc_id == 0) printf("\tinner solver performance:  %.2e Gflop/s\n", double(allflops) / double(maxtime) / 1.0e9);
      		  /*
      		  printf("this is for checking:\n");
      		  printf("\ttotal mixed solver time:   %.2e sec\n", double(stopeffective-starteffective));
      		  printf("\tfloating point operations: %.2e flops\n", effectiveflops);
      		  printf("\tinner solver performance:  %.2e Gflop/s\n", double(effectiveflops) / double(stopeffective-starteffective) / 1.0e9);
      		  */
      		#endif
      		
      		#if defined(GPU_BENCHMARK) || defined(GPU_BENCHMARK2)
      		  // REMARK: device_flops has to be multiplied by N_floats !!
      		  flops = device_flops * N_floats_int / (double(totalinnerclocks)/double(CLOCKS_PER_SEC)) / 1.0e9;
      		  printf("Inner solver BENCHMARK:\n");
      		  printf("\ttotal inner solver time:   %.2e sec\n", double(totalinnerclocks) / double(CLOCKS_PER_SEC));
      		  printf("\tfloating point operations: %.2e flops\n", double(device_flops) * double(N_floats_int));
      		  printf("\tinner solver performance:  %.2e Gflop/s\n", flops);
      		#endif
      		
      		#ifdef CPU_BENCHMARK
      		  // REMARK: host_flops has to be multiplied by N_floats !!
      		  flops = host_flops * N_floats_int / (double(totalouterclocks)/double(CLOCKS_PER_SEC)) / 1.0e9;
      		  printf("Outer solver BENCHMARK:\n");
      		  printf("\ttotal outer solver time:   %.2e sec\n", double(totalouterclocks) / double(CLOCKS_PER_SEC));
      		  printf("\tfloating point operations: %.2e flops\n", double(host_flops) * double(N_floats_int));
      		  printf("\touter solver performance:  %.2e Gflop/s\n", flops);
      		#endif
      
      
      #ifdef USETEXTURE
        unbind_texture_gf();
      #endif
      
      		// debug	// CUDA
  		#ifdef CUDA_DEBUG
  		  CUDA_CHECK("CUDA error in unbind_texture(). Unbindung the GF texture failed.", "GF texture unbound.");
  		#endif
      
      /*
      // cublasShutdown();
      
      		// debug	// CUBLAS helper function
  		#ifdef CUDA_DEBUG
  		  CUBLAS_HELPER_CHECK(cublasShutdown(), "CUBLAS error in cublasShutdown(). Couldn't shut down CUBLAS.", "CUBLAS is shutted down.");
  		#else
  		  cublasShutdown();
  		#endif
      */
      
      		// debug
      		if (g_cart_id == 0) printf("finalize_mixedsolve_eo_nd():\n");
      
      finalize_mixedsolve_eo_nd_mpi();
      
      		// debug
      		if (g_cart_id == 0) printf("\n");
      
      return(outercount);
    }
    
    
  }//OUTER LOOP
  while (outercount <= max_iter);
  
  
  // multiplying with Qhat(2x2)^dagger is done in invert_doublet_eo.c
  
  
  // timer
  stopouter = clock();
  totalouterclocks = stopouter-startouter - totalinnerclocks;
  
  #ifdef EFFECTIVE_BENCHMARK
    stopeffective = MPI_Wtime();
  #endif
  
  		/*
  		// benchmark
  		#ifdef GPU_BENCHMARK2
  		  device_flops = 0;
  		  int help = ( 4 + outercount*(1448+5*4+6) + outercount/10*1448 ) * N_floats;
  		  flopcount(device_flops, help);			// N_recalcres = 10
  		#endif
  		*/
  
  		// debug
  		if (g_cart_id == 0) {
  		  printf("\nEO inversion done in mixed precision.\n");
  		  printf("Finished outer loop, because of maximal number of outer iterations.\n");
      		  printf("Total number of inner iterations: %i\n", outercount);
      		  printf("Total number of outer iterations: %i\n", i+1);
      		  printf("Squared residue: %.10e\n", rr); 
      		  printf("Outer solver done in: %.4e sec\n", double(stopouter-startouter)/CLOCKS_PER_SEC);
      		}
      		// benchmark
      		#ifdef EFFECTIVE_BENCHMARK
      		  // will now count the number of effective flops
      		  // effectiveflops  =  #(inner iterations)*(matrixflops+linalgflops)*VOLUME/2  +  #(outer iterations)*(matrixflops+linalgflops)*VOLUME/2
      		  // outer loop: linalg  =  flops for calculating  r(k+1) and x(k+1)
      		  // inner loop: linalg  =  flops for calculating  alpha, x(k+1), r(k+1), beta, d(k+1)
      		  singletime = double(stopeffective-starteffective);
      		  effectiveflops = outercount*(matrixflops + 2*2*2*24 + 2*2*24 + 2*2*24 + 2*2*2*24 + 2*2*24)*VOLUME/2   +   i*(matrixflops + 2*24 + 2*24)*VOLUME/2;
      		  MPI_Allreduce(&singletime, &maxtime, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      		  MPI_Allreduce(&effectiveflops, &allflops, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      		  if (g_proc_id == 0) printf("effective BENCHMARK:\n");
      		  if (g_proc_id == 0) printf("\ttotal mixed solver time:   %.2e sec\n", double(maxtime));
      		  if (g_proc_id == 0) printf("\tfloating point operations: %.2e flops\n", double(allflops));
      		  if (g_proc_id == 0) printf("\tinner solver performance:  %.2e Gflop/s\n", double(allflops) / double(maxtime) / 1.0e9);
      		  /*
      		  printf("this is for checking:\n");
      		  printf("\ttotal mixed solver time:   %.2e sec\n", double(stopeffective-starteffective));
      		  printf("\tfloating point operations: %.2e flops\n", effectiveflops);
      		  printf("\tinner solver performance:  %.2e Gflop/s\n", double(effectiveflops) / double(stopeffective-starteffective) / 1.0e9);
      		  */
      		#endif
      		
      		#if defined(GPU_BENCHMARK) || defined(GPU_BENCHMARK2)
      		  // REMARK: device_flops has to be multiplied by N_floats !!
      		  flops = device_flops * N_floats_int / (double(totalinnerclocks)/double(CLOCKS_PER_SEC)) / 1.0e9;
      		  printf("Inner solver BENCHMARK:\n");
      		  printf("\ttotal inner solver time:   %.2e sec\n", double(totalinnerclocks) / double(CLOCKS_PER_SEC));
      		  printf("\tfloating point operations: %.2e flops\n", double(device_flops) * double(N_floats_int));
      		  printf("\tinner solver performance:  %.2e Gflop/s\n", flops);
      		#endif
      		
      		#ifdef CPU_BENCHMARK
      		  // REMARK: host_flops has to be multiplied by N_floats !!
      		  flops = host_flops * N_floats_int / (double(totalouterclocks)/double(CLOCKS_PER_SEC)) / 1.0e9;
      		  printf("Outer solver BENCHMARK:\n");
      		  printf("\ttotal outer solver time:   %.2e sec\n", double(totalouterclocks) / double(CLOCKS_PER_SEC));
      		  printf("\tfloating point operations: %.2e flops\n", double(host_flops) * double(N_floats_int));
      		  printf("\touter solver performance:  %.2e Gflop/s\n", flops);
      		#endif
  
  
  #ifdef USETEXTURE
    unbind_texture_gf();
  #endif
      
      		// debug	// CUDA
  		#ifdef CUDA_DEBUG
  		  CUDA_CHECK("CUDA error in unbind_texture(). Unbindung the GF texture failed.", "GF texture unbound.");
  		#endif
      
      /*
      // cublasShutdown();
      
      		// debug	// CUBLAS helper function
  		#ifdef CUDA_DEBUG
  		  CUBLAS_HELPER_CHECK(cublasShutdown(), "CUBLAS error in cublasShutdown(). Couldn't shut down CUBLAS.", "CUBLAS is shutted down.");
  		#else
  		  cublasShutdown();
  		#endif
      */
  
  		// debug
  		if (g_cart_id == 0) printf("finalize_mixedsolve_eo_nd():\n");  
  
  finalize_mixedsolve_eo_nd_mpi();
  
  		// debug
  		if (g_cart_id == 0) printf("\n");
  
  return(outercount);
  
  
}//mixedsolve_eo_nd()















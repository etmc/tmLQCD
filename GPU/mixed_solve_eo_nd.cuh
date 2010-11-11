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




	//////////////////////////////////////////////////////////////////
	//								//
	//    this is the implementation of the EO, ND mixed solver	//
	//								//
	//////////////////////////////////////////////////////////////////




///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




// the debugging functions can be included here via:	#include "./DEBUG/MATRIX_DEBUG.cuh"
//							#include "./DEBUG/MATRIX_MPI_DEBUG.cuh"



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




#ifdef HAVE_CONFIG_H
  #include<config.h>
#endif

extern "C" {
#include "../Nondegenerate_Matrix.h"
#include "../Hopping_Matrix.h"
#include "../solver/cg_her_nd.h"
#ifdef MPI
  #include "xchange.h"
  #include "communication.h"
#endif
}
/*
#ifdef MPI				// have to divide host- and device- and communication- and execution code more strictly for MPI purposes
  #include <mpi.h>			//	code with kernels and kernel launches syntax has to be compiled with nvcc
  //#include "mpi_init.h"
#endif
*/

#include "MACROS.cuh"


// spinor fields (pointing to device)
dev_spinor * dev_spin1_up;		// auxiliary fields for cg_eo_nd()
dev_spinor * dev_spin1_dn;
dev_spinor * dev_spin2_up;
dev_spinor * dev_spin2_dn;
dev_spinor * dev_spin3_up;
dev_spinor * dev_spin3_dn;
/*
dev_spinor * dev_spin4_up;
dev_spinor * dev_spin4_dn;
dev_spinor * dev_spin5_up;
dev_spinor * dev_spin5_dn;
*/

dev_spinor * dev_spinin_up;		// host/device interaction	// mixedsolve_eo_nd()  <-->  cg_eo_nd()
dev_spinor * dev_spinin_dn;		// inner/outer interaction
dev_spinor * dev_spinout_up;
dev_spinor * dev_spinout_dn;

dev_spinor * h2d_spin_up;		// for transferring in double precision on host to single precision on device (pointing to host)
dev_spinor * h2d_spin_dn;

dev_spinor * dev_spin_eo1_up;		// auxiliary for  matrix_multiplication32()  called by  dev_cg_eo_nd()
dev_spinor * dev_spin_eo1_dn;
dev_spinor * dev_spin_eo2_up;
dev_spinor * dev_spin_eo2_dn;
dev_spinor * dev_spin_eo3_up;
dev_spinor * dev_spin_eo3_dn;


// physical parameters (on device)
__device__ float mubar, epsbar;


// benchmark
#if defined(GPU_BENCHMARK) || defined(GPU_BENCHMARK2)
  unsigned long long int device_flops = 0;	// attention: integer overflow ....
#endif

#ifdef CPU_BENCHMARK
  unsigned long long int host_flops = 0;
#endif


#ifdef MPI && PARALLELT				// collecting variables for the MPI implementation
  int * iseven;
  int * dev_g_iup;
  int * dev_g_idn;
  int * dev_g_lexic2eo;
  int * dev_g_lexic2eosub;
  int * dev_g_eo2lexic;
  int * dev_g_ipt;
  spinor * spinor_xchange;			// for xchange_field_wrapper()
  spinor * spinor_debug_in;			// for Hopping_Matrix_wrapper()
  spinor * spinor_debug_out;			// for Hopping_Matrix_wrapper()
  __device__ int dev_RAND;			// not used, maybe later ...
  __device__ int dev_VOLUMEPLUSRAND;		// is now used in dev_Hopping_Matrix_mpi()
  //__device__ int dev_rank;			// was for the moment put to mixed_solve.cu ...
  //__device__ int dev_nproc;
#endif

/*
#if defined(MPI) && defined(PARALLELT)		// put to mixed_solve.cu
  #include "MPI.cuh"
#endif
*/






///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////






// host/device interaction

// remark: the host spinors are double precision and therefore need twice the memory !!
//		dev_spinor * device:    dev_spinsize
//		spinor * host:          2*dev_spinsize
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





// puts the additional nd parameters mubar and epsbar on the device
__global__ void he_cg_init_nd_additional (float param_mubar, float param_epsbar) {

  mubar  = param_mubar;
  epsbar = param_epsbar;

}






// derived from Flo's function  dev_mul_one_pm_imu_inv
//	order of the arguments also like Flo's convention: (spinin, spinout)

// applies (1 +- imubar*gamma5)
// uses shared local memory for manipulation	// really ??	where ??
// one thread per lattice site


__global__ void dev_mul_one_pm_imubar_gamma5 (dev_spinor * sin,
                                              dev_spinor * sout,
                                              float sign         ) {
   
  dev_spinor slocal[6];									// dev_spinor = float4		// 6*float4 = 24 floats		// auxiliary for each thread
  
  dev_complex pm_imu = dev_initcomplex(0.0, sign * mubar);				// dev_complex = struct { REAL re; REAL im; }	// pm_imu.re = 0.0
  																	// pm_imu.im = sign * mubar
  int pos = threadIdx.x + blockDim.x*blockIdx.x;
  
  if (pos < dev_VOLUME) {
    dev_skalarmult_gamma5_spinor(&(slocal[0]), pm_imu, &(sin[6*pos]) );			// slocal  =  pm_imu * (gamma5) * sin
    dev_add_spinor_assign(&(slocal[0]), &(sin[6*pos]));					// slocal  =  slocal + sin  =  pm_imu * (gamma5) * sin + sin
    dev_realmult_spinor_assign(&(sout[6*pos]), 1.0, &(slocal[0]) );			// sout    =  slocal
  }
}





// TEST:
/*
void init_nnspinor_eo_test() {
									
  int x, y, z, t, ind, nnpos, j;					// mixed_solve_eo(...) allocates 8 integers per even or odd lattice site: size_t nnsize = 8*VOLUME*sizeof(int);
  									
  for (t = 0; t < T; t++) {						// loop goes over all INTERN latice sites !!
    for (x = 0; x < LX; x++) {						// doesn't refer to any EXTERN BOUNDARIES !!  ->  CORRESPONDS TO THE WHOLE LATTICE (I.E. WHEN NO SUBLATTICES ARE ASSIGNED) !!
      for (y = 0; y < LY; y++) {					//						  because of the behaviour of  g_iup[][] in the non-parallel case
        for (z = 0; z < LZ; z++) {
        								// NOTICE: g_ipt, g_iup, g_idn, and g_lexic2eosub  refer to pos. of lin. proj. pos. of the lattice
          ind = g_ipt[t][x][y][z];					// g_ipt[t][x][y][z] 	returns the linearly projected position of (t,x,y,z) of the lattice
          								//	indexes computed in geometry_eo() from geometry_eo.c
          								//	memory for the index array allocated by init_geometry_indices() from init_geometry_indices.c
          //if ((t+x+y+z)%2 == 0) { // EVEN
          if ((t + x + y + z + g_proc_coords[0]*T  + g_proc_coords[1]*LX + 
	  	               g_proc_coords[2]*LY + g_proc_coords[3]*LZ) % 2 == 0) {
          
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
*/


/*
void init_idxgauge_test() {		// works

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
  
  printf("This was init_idxgauge_mpi().\n");
  
}
*/


/*
void init_idxgauge_test() {		// works

  int pos_eo, pos_global_even, pos_global_odd;
  
  for (pos_eo = 0; pos_eo < VOLUME/2; pos_eo++) {
      // even
      pos_global_even = g_eo2lexic[pos_eo];
      eoidx_even[pos_eo] = pos_global_even;
      // odd
      pos_global_odd = g_eo2lexic[(VOLUME+RAND)/2 + pos_eo];
      eoidx_odd[pos_eo] = pos_global_odd;
  }
  
  printf("This was init_idxgauge_mpi().\n");
  
}
*/










// initializes and allocates all quantities for the mixed solver
// more precise:
//	puts the gauge field on device as "2 rows" or "8 floats" per SU(3)-matrix
//	allocates memory for all spinor fields
//	puts the nn- and eoidx-fields on device memory

void init_mixedsolve_eo_nd(su3** gf) {	// gf is the full gauge field
  
  
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
  
  
  
  
  // get number of devices
  ndev = find_devices();
  if(ndev == 0){
    fprintf(stderr, "Error: no CUDA devices found. Aborting...\n");
    exit(300);
  }
  // try to set active device to device_num given in input file
  if(device_num < ndev){
    printf("Setting active device to: %d\n", device_num);
    cudaSetDevice(device_num);
  }
  else{
    fprintf(stderr, "Error: There is no CUDA device with No. %d. Aborting...\n",device_num);
    exit(301);
  }
  if((cudaerr=cudaGetLastError())!=cudaSuccess){
    printf("Error in init_mixedsolve_eo_nd(): Could not set active device. Aborting...\n");
    exit(302);
  }
  
  
  
  
  #ifdef USETEXTURE
    printf("Using texture references.\n");
  #else
    printf("NOT using texture references.\n");
  #endif
  
  
  
  
  /////////////////
  // GAUGE FIELD //
  /////////////////
  #ifdef GF_8
    /* allocate 8 floats for gf = 2*4*VOLUME float4's*/
    printf("Using GF 8 reconstruction.\n");			// dev_su3_8 = float4
    dev_gfsize = 4*VOLUME * 2*sizeof(dev_su3_8);		// allocates for each lattice site and for 4 directions  2*float4 = 8 floats  = 8 real parameters
  #else
    /* allocate 2 rows of gf = 3*4*VOLUME float4's*/
    printf("Using GF 12 reconstruction.\n");			// dev_su3_2v = float4
    dev_gfsize = 4*VOLUME * 3*sizeof(dev_su3_2v); 		// allocates for each lattice site and for 4 directions  3*float4 = 12 floats = 2 rows of complex 3-vectors
  #endif
  
  
  if ( (cudaerr = cudaMalloc((void **) &dev_gf, dev_gfsize)) != cudaSuccess ){	// allocates memory for the gauge field dev_gf on device
    printf("Error in init_mixedsolve_eo_nd(): Memory allocation of gauge field failed. Aborting...\n");
    printf("%s\n", cudaGetErrorString(cudaGetLastError()));
    exit(200);
  }
  else{
    printf("Allocated gauge field on device.\n");
  }  
  
  
  #ifdef GF_8
    // Allocate REAL conversion gf on host ??
    h2d_gf = (dev_su3_8 *)malloc(dev_gfsize); 			// allocates on host
    su3to8(gf,h2d_gf);						// h2d_gf  is the gauge field  gf  with the 8-real-parameter-representation (according to M. Clark, p. 28)
  #else
    // Allocate REAL conversion gf on host ??
    h2d_gf = (dev_su3_2v *)malloc(dev_gfsize);			// allocates on host
    su3to2vf4(gf,h2d_gf);					// h2d_gf  is the gauge field  gf  with the first two rows stored
  #endif
  
  cudaMemcpy(dev_gf, h2d_gf, dev_gfsize, cudaMemcpyHostToDevice);
  								// dev_gf = h2d_gf  on device memory
  
  		// debug	// CUDA
  		#ifdef CUDA_DEBUG
  		  CUDA_CHECK("CUDA error in init_mixedsolve_eo_nd(). Copying dev_gf to device failed.", "Allocated dev_gf on device.");
  		#endif
  
  
  
  
  //////////
  // GRID //
  //////////
  nnsize = 8*VOLUME*sizeof(int);				// size of memory for 8*VOLUME integers
  nn = (int *) malloc(nnsize);					// allocate this memory on host
  nn_eo = (int *) malloc(nnsize/2);				// allocate half this memory
  nn_oe = (int *) malloc(nnsize/2);				// allocate half this memory
  cudaMalloc((void **) &dev_nn, nnsize);			// memory on device
  cudaMalloc((void **) &dev_nn_eo, nnsize/2);			// half the memory on device
  cudaMalloc((void **) &dev_nn_oe, nnsize/2);			// half the memory on device
  
  
  idxsize = VOLUME/2*sizeof(int);				// size of memory necessary for VOLUME/2 integers
  eoidx_even = (int *) malloc(idxsize);				// allocate on host
  eoidx_odd = (int *) malloc(idxsize);				// allocate on host
  cudaMalloc((void **) &dev_eoidx_even, idxsize);		// allocate on device
  cudaMalloc((void **) &dev_eoidx_odd, idxsize);		// allocate on device
  
  
  initnn();							// initialize nearest-neighbour table for gpu
  initnn_eo();							// initialize nearest-neighbour table for gpu with even-odd enabled
  // test:
  //init_nnspinor_eo_test();
  //init_idxgauge_test();
  
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
  		  CUDA_CHECK("CUDA error in init_mixedsolve_eo_nd(). Memory allocation of grid stuff failed.", "Allocated grid stuff on device.");
  		#endif
  
  
  
  
  /////////////
  // SPINORS //							// allocates device memory for the odd part of the spinor fields (24 = 4*6 floats per odd lattice sites)
  /////////////							// now we have to consider 2 flavors: up, dn
  
  dev_spinsize = 6*VOLUME/2*sizeof(dev_spinor);			// remember: dev_spinor = float4
  
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
  		if ( (void *) (h2d_spin_up = (dev_spinor *) malloc(dev_spinsize) ) == NULL) {
  		  printf("Could not allocate memory for h2d_spin_up. Aborting...\n");
  		  exit(200);
  		}
  		if ( (void *) (h2d_spin_dn = (dev_spinor *) malloc(dev_spinsize) ) == NULL) {
  		  printf("Could not allocate memory for h2d_spin_dn. Aborting...\n");
  		  exit(200);
  		}
  // h2d_spin_up = (dev_spinor *) malloc(dev_spinsize);		// for transferring the spin field in double precision on host to single precision on device
  // h2d_spin_dn = (dev_spinor *) malloc(dev_spinsize);		// on host pointing to host
  
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
  		  CUDA_CHECK("CUDA error in init_mixedsolve_eo_nd(). Memory allocation of spinor fields failed.", "Allocated spinor fields on device.");
  		#endif
  
  
  
  
  ////////////
  // output //						// ??
  ////////////
  output_size = LZ*T*sizeof(float); 			// parallel in t and z direction
  cudaMalloc((void **) &dev_output, output_size);	// output array
  host_output = (float *) malloc(output_size);
  
  
  		// debug	// CUDA
  		#ifdef CUDA_DEBUG
  		  CUDA_CHECK("CUDA error in init_mixedsolve_eo_nd(). Memory allocation output stuff failed.", "Allocated output stuff on device.");
  		#endif
  
  
  
  
  ////////////////////////////
  // grid[ ] specifications //							// allocate and initializes the array grid[5] on device
  ////////////////////////////
  grid[0]=LX; grid[1]=LY; grid[2]=LZ; grid[3]=T; grid[4]=VOLUME/2;		// it contains the dimensions of the lattice and the volume of the eo-sublattice
  										// dev_VOLUME is half of VOLUME for eo
 
  cudaMalloc((void **) &dev_grid, 5*sizeof(int));				// dev_grid
  cudaMemcpy(dev_grid, &(grid[0]), 5*sizeof(int), cudaMemcpyHostToDevice);
  
  
  		// debug	// CUDA
  		#ifdef CUDA_DEBUG
  		  CUDA_CHECK("CUDA error in init_mixedsolve_eo_nd(). Memory allocation of grid[] specifications failed.", "Allocated grid[] specifications on device.");
  		#endif
  
  
}//init_mixedsolve_eo_nd()





// deallocates the previous allocated memory

void finalize_mixedsolve_eo_nd(void) {
  
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
  
  
  // Clean up CUDA API for calling thread	// ??
  cudaThreadExit();				// is essential
  
  
  		// debug	// CUDA
  		#ifdef CUDA_DEBUG
  		  CUDA_CHECK("CUDA error in finalize_mixedsolve_eo_nd(). Device memory deallocation failed", "Device memory deallocated.");
  		#endif
  
}






///////////////////////////
// MATRIX MULTIPLICATION //
///////////////////////////

// the GPU implementation of  Q_Qdagger_ND(...)  from Nondegenerate_Matrix.c
//	Flo's equivalent function for the standard and non-nd case is  dev_Qtm_pm_psi

void matrix_multiplication32 (dev_spinor * spinout_up, dev_spinor * spinout_dn,
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
  
  // Flo:
  #ifdef USETEXTURE
    bind_texture_spin(spinin_dn,1);
  #endif
    dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, spinin_dn, dev_spin_eo1_up, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0);	// dev_spin_eo1_up = (M_eo) * spinin_dn
  #ifdef USETEXTURE
    unbind_texture_spin(1);
  #endif
  
  #ifdef USETEXTURE
    bind_texture_spin(spinin_up,1);
  #endif
    dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, spinin_up, dev_spin_eo1_dn, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0);	// dev_spin_eo1_dn = (M_eo) * spinin_up
  #ifdef USETEXTURE
    unbind_texture_spin(1);
  #endif
  
  
  // written:
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo1_up, dev_spin_eo2_up, -1.0);		// dev_spin_eo2_up  =  (1 - imubar) * dev_spin_eo1_up  =  (1 - imubar)*(M_eo) * spinin_dn
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo1_dn, dev_spin_eo2_dn, +1.0);		// dev_spin_eo2_dn  =  (1 + imubar) * dev_spin_eo1_dn  =  (1 + imubar)*(M_eo) * spinin_up
  
  
  // CUBLAS:
  cublasSaxpy (N_floats, g_epsbar, (float *) dev_spin_eo1_dn, 1, (float *) dev_spin_eo2_up, 1);
  										// dev_spin_eo2_up  =  dev_spin_eo2_up  +  epsbar * dev_spin_eo1_dn  =  (1 - imubar)*(M_eo) * spinin_dn  +  epsbar * (M_eo) * spinin_up
  cublasSaxpy (N_floats, g_epsbar, (float *) dev_spin_eo1_up, 1, (float *) dev_spin_eo2_dn, 1);
  										// dev_spin_eo2_dn  =  dev_spin_eo2_dn  +  epsbar * dev_spin_eo1_up  =  (1 + imubar)*(M_eo) * spinin_up  +  epsbar * (M_eo) * spinin_dn
  
  // CUBLAS:
  cublasSscal (N_floats, nrm, (float *) dev_spin_eo2_up, 1);			// dev_spin_eo2_up  =  nrm * dev_spin_eo2_up  =  nrm*(1-imubar)*(M_eo) * spinin_dn  +  nrm*epsbar*(M_eo) * spinin_up
  
  cublasSscal (N_floats, nrm, (float *) dev_spin_eo2_dn, 1);			// dev_spin_eo2_dn  =  nrm * dev_spin_eo2_dn  =  nrm*(1+imubar)*(M_eo) * spinin_up  +  nrm*epsbar*(M_eo) * spinin_dn
  
  
  // Flo:
  #ifdef USETEXTURE
    bind_texture_spin(dev_spin_eo2_up,1);									// remember: this is ((M_oe)(Mee^-1)(M_eo)) * (spinin_dn, spinin_up):
  #endif
    dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo2_up, dev_spin_eo1_up, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1);
  #ifdef USETEXTURE
    unbind_texture_spin(1);							// dev_spin_eo1_up  =  (M_oe) * dev_spin_eo2_up  =  (M_oe)*nrm*(1-imubar)*(M_eo) * spinin_dn  +  (M_oe)*nrm*epsbar*(M_eo) * spinin_up
  #endif									
  
  #ifdef USETEXTURE
    bind_texture_spin(dev_spin_eo2_dn,1);
  #endif
    dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo2_dn, dev_spin_eo1_dn, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1);
  #ifdef USETEXTURE
    unbind_texture_spin(1);							// dev_spin_eo1_dn  =  (M_oe) * dev_spin_eo2_dn  =  (M_oe)*nrm*(1+imubar)*(M_eo) * spinin_up  +  (M_oe)*nrm*epsbar*(M_eo) * spinin_dn
  #endif
  
  
  
  
  ////////////
  // (M_oo) //
  ////////////
  
  
  // written:
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(spinin_dn, dev_spin_eo2_up, +1.0);			// dev_spin_eo2_up = (1 + imubar) * spinin_dn
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(spinin_up, dev_spin_eo2_dn, -1.0);			// dev_spin_eo2_dn = (1 - imubar) * spinin_up
  
  
  // CUBLAS:													// remember: this is (M_oo) * (spinin_dn, spinin_up):
  cublasSaxpy (N_floats, -g_epsbar, (float *) spinin_up, 1, (float *) dev_spin_eo2_up, 1);
  												// dev_spin_eo2_up  =  dev_spin_eo2_up - epsbar*spinin_up  =  (1+imubar)*spinin_dn - epsbar*spinin_up
  cublasSaxpy (N_floats, -g_epsbar, (float *) spinin_dn, 1, (float *) dev_spin_eo2_dn, 1);
  												// dev_spin_eo2_dn  =  dev_spin_eo2_dn - epsbar*spinin_dn  =  (1-imubar)*spinin_up - epsbar*spinin_dn
  
  
  
  ///////////////////////////////////////
  // (M_oo)  -  (M_oe) (Mee^-1) (M_eo) //
  ///////////////////////////////////////
  
  // CUBLAS:													// this is ((M_oo) - (M_oe)(Mee^-1)(M_eo)) * (spinin_dn, spinin_up):
  cublasSaxpy (N_floats, -1.0, (float *) dev_spin_eo1_up, 1, (float *) dev_spin_eo2_up, 1);
  							// dev_spin_eo2_up  =  dev_spin_eo2_up  -  dev_spin_eo1_up  =                      (1+imubar) * spinin_dn  -                    epsbar * spinin_up
  							//                                                             - (M_oe)*nrm*(1-imubar)*(M_eo) * spinin_dn  -  (M_oe)*nrm*epsbar*(M_eo) * spinin_up
  cublasSaxpy (N_floats, -1.0, (float *) dev_spin_eo1_dn, 1, (float *) dev_spin_eo2_dn, 1);
  							// dev_spin_eo2_dn  =  dev_spin_eo2_dn  -  dev_spin_eo1_dn  =                      (1-imubar) * spinin_up  -                    epsbar * spinin_dn
  							//                                                             - (M_oe)*nrm*(1+imubar)*(M_eo) * spinin_up  -  (M_oe)*nrm*epsbar*(M_eo) * spinin_dn
  
  
  ////////////
  // gamma5 //
  ////////////
  
  // Flo:
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
  
  // Flo:
  #ifdef USETEXTURE
    bind_texture_spin(dev_spin_eo3_up,1);
  #endif
    dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo3_up, dev_spin_eo1_up, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0);	// dev_spin_eo1_up = (M_eo) * dev_spin_eo3_up
  #ifdef USETEXTURE
    unbind_texture_spin(1);
  #endif
  
  #ifdef USETEXTURE
    bind_texture_spin(dev_spin_eo3_dn,1);
  #endif
    dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo3_dn, dev_spin_eo1_dn, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0);	// dev_spin_eo1_dn = (M_eo) * dev_spin_eo3_dn
  #ifdef USETEXTURE
    unbind_texture_spin(1);
  #endif
  
  
  // written:
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo1_up, dev_spin_eo2_up, -1.0);			// dev_spin_eo2_up  =  (1 - imubar) * dev_spin_eo1_up  =  (1 - imubar)*(M_eo) * dev_spin_eo3_up
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo1_dn, dev_spin_eo2_dn, +1.0);			// dev_spin_eo2_dn  =  (1 + imubar) * dev_spin_eo1_dn  =  (1 + imubar)*(M_eo) * dev_spin_eo3_dn
  
  
  // CUBLAS:
  cublasSaxpy (N_floats, g_epsbar, (float *) dev_spin_eo1_dn, 1, (float *) dev_spin_eo2_up, 1);
  										// dev_spin_eo2_up  =  dev_spin_eo2_up  +  epsbar * dev_spin_eo1_dn  =  (1 - imubar)*(M_eo) * dev_spin_eo3_up  +  epsbar * (M_eo) * dev_spin_eo3_dn
  cublasSaxpy (N_floats, g_epsbar, (float *) dev_spin_eo1_up, 1, (float *) dev_spin_eo2_dn, 1);
  										// dev_spin_eo2_dn  =  dev_spin_eo2_dn  +  epsbar * dev_spin_eo1_up  =  (1 + imubar)*(M_eo) * dev_spin_eo3_dn  +  epsbar * (M_eo) * dev_spin_eo3_up
  
  // CUBLAS:
  cublasSscal (N_floats, nrm, (float *) dev_spin_eo2_up, 1);			// dev_spin_eo2_up  =  nrm * dev_spin_eo2_up  =  nrm*(1-imubar)*(M_eo) * dev_spin_eo3_up  +  nrm*epsbar*(M_eo) * dev_spin_eo3_dn
  
  cublasSscal (N_floats, nrm, (float *) dev_spin_eo2_dn, 1);			// dev_spin_eo2_dn  =  nrm * dev_spin_eo2_dn  =  nrm*(1+imubar)*(M_eo) * dev_spin_eo3_dn  +  nrm*epsbar*(M_eo) * dev_spin_eo3_up
  
  
  // Flo:
  #ifdef USETEXTURE
    bind_texture_spin(dev_spin_eo2_up,1);									// remember: this is ((M_oe) (Mee^-1) (M_eo)) * (dev_spin_eo3_up, dev_spin_eo3_dn):
  #endif
    dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo2_up, dev_spin_eo1_up, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1);
  #ifdef USETEXTURE
    unbind_texture_spin(1);							// dev_spin_eo1_up  =  (M_oe) * dev_spin_eo2_up  =  (M_oe)*nrm*(1-imubar)*(M_eo) * dev_spin_eo3_up  +  (M_oe)*nrm*epsbar*(M_eo) * dev_spin_eo3_dn
  #endif
  
  #ifdef USETEXTURE
    bind_texture_spin(dev_spin_eo2_dn,1);
  #endif
    dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo2_dn, dev_spin_eo1_dn, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1);
  #ifdef USETEXTURE
    unbind_texture_spin(1);							// dev_spin_eo1_dn  =  (M_oe) * dev_spin_eo2_dn  =  (M_oe)*nrm*(1+imubar)*(M_eo) * dev_spin_eo3_dn  +  (M_oe)*nrm*epsbar*(M_eo) * dev_spin_eo3_up
  #endif
  
  
  
  ////////////
  // (M_oo) //
  ////////////
  
  // written:
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo3_up, dev_spin_eo2_up, +1.0);				// dev_spin_eo2_up = (1 + imubar) * dev_spin_eo3_up
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo3_dn, dev_spin_eo2_dn, -1.0);				// dev_spin_eo2_dn = (1 - imubar) * dev_spin_eo3_dn
  
  
  // CUBLAS:											// remember: this is (M_oo) * (dev_spin_eo3_up, dev_spin_eo3_dn):
  cublasSaxpy (N_floats, -g_epsbar, (float *) dev_spin_eo3_dn, 1, (float *) dev_spin_eo2_up, 1);
  												// dev_spin_eo2_up  =  dev_spin_eo2_up - epsbar*dev_spin_eo3_dn  =  (1+imubar)*dev_spin_eo3_up - epsbar*dev_spin_eo3_dn
  cublasSaxpy (N_floats, -g_epsbar, (float *) dev_spin_eo3_up, 1, (float *) dev_spin_eo2_dn, 1);
  												// dev_spin_eo2_dn  =  dev_spin_eo2_dn - epsbar*dev_spin_eo3_up  =  (1-imubar)*dev_spin_eo3_dn - epsbar*dev_spin_eo3_up
  
  
  
  ///////////////////////////////////////
  // (M_oo)  -  (M_oe) (Mee^-1) (M_eo) //
  ///////////////////////////////////////
  
  // CUBLAS:											// this is ( (M_oo) - (M_oe) (Mee^-1) (M_eo) ) * (dev_spin_eo3_up, dev_spin_eo3_dn)
  cublasSaxpy (N_floats, -1.0, (float *) dev_spin_eo1_up, 1, (float *) dev_spin_eo2_up, 1);
  							// dev_spin_eo2_up  =  dev_spin_eo2_up  -  dev_spin_eo1_up  =                      (1+imubar) * dev_spin_eo3_up  -                    epsbar * dev_spin_eo3_dn
  							//                                                             - (M_oe)*nrm*(1-imubar)*(M_eo) * dev_spin_eo3_up  -  (M_oe)*nrm*epsbar*(M_eo) * dev_spin_eo3_dn
  cublasSaxpy (N_floats, -1.0, (float *) dev_spin_eo1_dn, 1, (float *) dev_spin_eo2_dn, 1);
  							// dev_spin_eo2_dn  =  dev_spin_eo2_dn  -  dev_spin_eo1_dn  =                      (1-imubar) * dev_spin_eo3_dn  -                    epsbar * dev_spin_eo3_up
  							//                                                             - (M_oe)*nrm*(1+imubar)*(M_eo) * dev_spin_eo3_dn  -  (M_oe)*nrm*epsbar*(M_eo) * dev_spin_eo3_up
  
  
  ////////////
  // gamma5 //
  ////////////
  
  // Flo:
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
  
}//matrix_multiplication32()






// will be used to count the floating point operations per (24 * (#lattice sites)) = (#floats)
// i.e total has to be multiplied by N_floats

void flopcount(unsigned long long int& total, int add) {

  total = total + add;
  
}







extern "C" void benchmark_eo_nd (spinor * Q_up, spinor * Q_dn, int N) {

  
  
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
  double timeElapsed;
  double startBenchmark;
  double stopBenchmark;
  
  // counter
  int i;
  
  // flop counting
  double realFlopsPerApp = 34768.0;
  double effectiveFlopsPerApp = 23984.0;
  
  double realDeviceFlops;
  double realFlops;
  double effectiveDeviceFlops;
  double effectiveFlops;
  
  // CUDA errors
  cudaError_t cudaerr;
  cublasStatus cublasstatus;
  
  // size of a spinor
  size_t dev_spinsize = 6*VOLUME/2 * sizeof(dev_spinor);
  
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
  		printf("\nStarting a little BENCHMARK. benchmark_eo_nd().\n");
  
  
  
  
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
  		printf("Applying the eo-preconditioned matrix %i times.\n", N);
  
  
  to_device(B_up, Q_up, h2d_spin_up, dev_spinsize);
  to_device(B_dn, Q_dn, h2d_spin_dn, dev_spinsize);
  
  
  // timer
  startBenchmark = double(clock()) / double(CLOCKS_PER_SEC);
  
  
  
  for (i = 0; i < N; i++) {
  
    matrix_multiplication32(A_up, A_dn,				// A = (matrix)*B
                            B_up, B_dn,
                            griddim2, blockdim2,
                            griddim3, blockdim3,
                            griddim4, blockdim4,
                            griddim5, blockdim5);
    
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
  stopBenchmark = double(clock()) / double(CLOCKS_PER_SEC);
  
  
  timeElapsed = stopBenchmark - startBenchmark;
  
  realDeviceFlops      = N * VOLUME/2 * realFlopsPerApp;
  realFlops            = N * VOLUME/2 * realFlopsPerApp / timeElapsed / 1.0e9;
  
  effectiveDeviceFlops = N * VOLUME/2 * effectiveFlopsPerApp;
  effectiveFlops       = N * VOLUME/2 * effectiveFlopsPerApp / timeElapsed / 1.0e9;
  
  
  printf("REAL:\n");
  printf("\ttime:        %.2e sec\n", timeElapsed);
  printf("\tflop's:      %.2e flops\n", realDeviceFlops);
  printf("\tperformance: %.2e Gflop/s\n\n", realFlops);
  
  printf("EFFECTIVE:\n");
  printf("\ttime:        %.2e sec\n", timeElapsed);
  printf("\tflop's:      %.2e flops\n", effectiveDeviceFlops);
  printf("\tperformance: %.2e Gflop/s\n\n", effectiveFlops);
  
  
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
  
  
}//benchmark_eo_nd()










////////////////////////
// CONJUGATE GRADIENT //
////////////////////////

// for the odd field after even/odd-preconditioning
// single precision on GPU

int cg_eo_nd (dev_su3_2v * gf,
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
  size_t dev_spinsize = 6*VOLUME/2*sizeof(dev_spinor);
  int N_sites  =    VOLUME/2;
  int N_floats = 24*VOLUME/2;		// (single precision) CUBLAS functions get the number of floats as input
  
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
  he_cg_init<<< 1, 1 >>> (dev_grid, (float) g_kappa, (float)(g_mu/(2.0*g_kappa)), h0, h1, h2, h3);
  
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
  rr_up = cublasSdot(N_floats, (float *) r_up, 1, (float *) r_up, 1);
  rr_dn = cublasSdot(N_floats, (float *) r_dn, 1, (float *) r_dn, 1);
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
    		printf("\nEntering inner loop.\n");
  
  		// debug	// CUBLAS core function
		#ifdef CUDA_DEBUG
		  // CUBLAS_CORE_CHECK("CUBLAS error in cg_eo_nd(). Calculating initial residue failed.", "Initial inner residue calculated.");
		  CUBLAS_CORE_CHECK_NO_SUCCESS_MSG("CUBLAS error in cg_eo_nd(). Calculating initial residue failed.");
		#endif
  
  		// debug
  		printf("Initial inner residue: %.6e\n", r0r0);
  
  
  for (j = 0; j < max_iter; j++) {
    
    
    #ifndef MATRIX_DEBUG
    
      // A*d(k)
      matrix_multiplication32(Ad_up, Ad_dn,										// normally:  matrix_multiplication32()
                               d_up,  d_dn,										// debugging: matrix_debug1(), matrix_multiplication_test()
                              griddim2, blockdim2,
                              griddim3, blockdim3,
                              griddim4, blockdim4,
                              griddim5, blockdim5);
    
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
    		to_host(g_chi_up_spinor_field[DUM_SOLVER+3], d_up, h2d_spin_up, dev_spinsize);
    		to_host(g_chi_dn_spinor_field[DUM_SOLVER+3], d_dn, h2d_spin_dn, dev_spinsize);
    		
    		// matrix multiplication
    		printf("This is Q_Qdagger_ND(). ");
    		Q_Qdagger_ND(g_chi_up_spinor_field[DUM_SOLVER+4], g_chi_dn_spinor_field[DUM_SOLVER+4],			// normally:  Q_Qdagger_ND()
    		             g_chi_up_spinor_field[DUM_SOLVER+3], g_chi_dn_spinor_field[DUM_SOLVER+3] );		// debugging: matrix_debug2(), Zwitter1(), Zwitter2(), Zwitter3()
    		
    		// host/device interaction
    		to_device(Ad_up, g_chi_up_spinor_field[DUM_SOLVER+4], h2d_spin_up, dev_spinsize);
    		to_device(Ad_dn, g_chi_dn_spinor_field[DUM_SOLVER+4], h2d_spin_dn, dev_spinsize);
    		
    		
    				// debug	// CUDA
  				#ifdef CUDA_DEBUG
  			  	// CUDA_CHECK("CUDA error in cg_eo_nd(). Applying the matrix on CPU failed.", "The matrix was applied on CPU.");
  			  	CUDA_CHECK_NO_SUCCESS_MSG("CUDA error in cg_eo_nd(). Applying the matrix on CPU failed.");
  				#endif
    
    #endif
    
    
    // alpha = r(k)*r(k) / d(k)*A*d(k)
    dAd_up = cublasSdot(N_floats, (float *) d_up, 1, (float *) Ad_up, 1);
    dAd_dn = cublasSdot(N_floats, (float *) d_dn, 1, (float *) Ad_dn, 1);
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
    cublasSaxpy(N_floats, alpha, (float *) d_up, 1, (float *) x_up, 1);
    cublasSaxpy(N_floats, alpha, (float *) d_dn, 1, (float *) x_dn, 1);
    
    		// benchmark
    		#ifdef GPU_BENCHMARK
    		  flopcount(device_flops, 2*2);
    		  // flopcount(device_flops, 2*2*N_floats);
    		#endif
    
    
    // r(k+1)
    if ( (j+1) % N_recalc_res != 0 ) {	// r(k+1) = r(k) - alpha*A*d(k)
      cublasSaxpy(N_floats, -1.0*alpha, (float *) Ad_up, 1, (float *) r_up, 1);
      cublasSaxpy(N_floats, -1.0*alpha, (float *) Ad_dn, 1, (float *) r_dn, 1);
      
      		// benchmark
      		#ifdef GPU_BENCHMARK
      		  flopcount(device_flops, 2*2);
      		  // flopcount(device_flops, 2*2*N_floats);
      		#endif
    }
    
    else {				// recalculate residue r(k+1) = b - A*x(k+1)
    					//	"feedback"
      		// debug
      		printf("Recalculating the inner residue.\n");
      
      
      #ifndef MATRIX_DEBUG
        // A*x(k+1)
        matrix_multiplication32(Ax_up, Ax_dn,
                                 x_up,  x_dn,
                                griddim2, blockdim2,
                                griddim3, blockdim3,
                                griddim4, blockdim4,
                                griddim5, blockdim5);
      #else
      		// debug	// apply the host matrix on trial
    		
    		// host/device interaction
    		to_host(g_chi_up_spinor_field[DUM_SOLVER+3], x_up, h2d_spin_up, dev_spinsize);
    		to_host(g_chi_dn_spinor_field[DUM_SOLVER+3], x_dn, h2d_spin_dn, dev_spinsize);
    		
    		// matrix multiplication
    		printf("This is Q_Qdagger_ND(). ");
    		Q_Qdagger_ND(g_chi_up_spinor_field[DUM_SOLVER+4], g_chi_dn_spinor_field[DUM_SOLVER+4],
    		             g_chi_up_spinor_field[DUM_SOLVER+3], g_chi_dn_spinor_field[DUM_SOLVER+3] );
    		
    		// host/device interaction
    		to_device(Ax_up, g_chi_up_spinor_field[DUM_SOLVER+4], h2d_spin_up, dev_spinsize);
    		to_device(Ax_dn, g_chi_dn_spinor_field[DUM_SOLVER+4], h2d_spin_dn, dev_spinsize);
    		
    		
    				// debug	// CUDA
  				#ifdef CUDA_DEBUG
  			  	// CUDA_CHECK("CUDA error in cg_eo_nd(). Applying the matrix on CPU failed.", "The matrix was applied on CPU.");
  			  	CUDA_CHECK_NO_SUCCESS_MSG("CUDA error in cg_eo_nd(). Applying the matrix on CPU failed.");
  				#endif
      #endif
      
      
      // r(k+1) = b - A*x(k+1)
      cublasScopy(N_floats, (float *) Q_up, 1, (float *) r_up, 1);		// r_up = Q_up
      cublasScopy(N_floats, (float *) Q_dn, 1, (float *) r_dn, 1);		// r_dn = Q_dn
      cublasSaxpy(N_floats, -1.0, (float *) Ax_up, 1, (float *) r_up, 1);	// r_up = Q_up - Ax_up
      cublasSaxpy(N_floats, -1.0, (float *) Ax_dn, 1, (float *) r_dn, 1);	// r_dn = Q_dn - Ax_dn
      
      		// benchmark
      		#ifdef GPU_BENCHMARK
      		  flopcount(device_flops, 2*2);
      		  // flopcount(device_flops, 2*2*N_floats);
      		#endif
    
    }
        
    
    // r(k+1)*r(k+1)
    rr_up  = cublasSdot(N_floats, (float *) r_up, 1, (float *) r_up, 1);
    rr_dn  = cublasSdot(N_floats, (float *) r_dn, 1, (float *) r_dn, 1);
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
    		printf("inner iteration j = %i: rr = %.6e\n", j, rr);
    		
    		// debug	// is NaN ?
    		if isnan(rr) {
    		  printf("Error in cg_eo_nd(). Inner residue is NaN.\n");
    		  exit(-1);
    		}
    
    
    // aborting ?? // check wether precision is reached ...
    if ( (check_abs)&&(rr <= eps_abs) || (check_rel)&&(rr <= eps_rel*r0r0) ) {
      
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
    cublasSscal (N_floats, beta, (float *) d_up, 1);
    cublasSaxpy (N_floats, 1.0 , (float *) r_up, 1, (float *) d_up, 1);
    
    cublasSscal (N_floats, beta, (float *) d_dn, 1);
    cublasSaxpy (N_floats, 1.0 , (float *) r_dn, 1, (float *) d_dn, 1);
    
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
  		printf("Finished inner loop beacuse of maximal number of inner iterations.\n");
  		printf("Final inner residue: %.6e\n", rr);
  
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

extern "C" int mixedsolve_eo_nd (spinor * P_up, spinor * P_dn,
                                 spinor * Q_up, spinor * Q_dn,
                                 int max_iter, double eps_sq, int rel_prec) {
  
  
  // basically  P_up/dn  and  Q_up/dn  could be used as auxiliary fields
  //	P_up/dn  is the output field (and can be used as initial guess)
  //	Q_up/dn  is not used later in the calling  invert_doublet_eo.c
  //		but will be used as feedback in r(k+1) = b - A*x(k+1)
  
  
  		// debug
  		printf("\n\nmixedsolve_eo_nd():\n");
  		
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
  #endif
  
  // timing
  clock_t startouter, stopouter;
  clock_t startinner, stopinner;
  // double timeelapsed;
  clock_t innerclocks;
  clock_t totalinnerclocks = 0;
  clock_t totalouterclocks = 0;
  #ifdef EFFECTIVE_BENCHMARK
    clock_t starteffective;
    clock_t stopeffective;
  #endif
  
  // (auxiliary) fields
  spinor *  r_up, *  r_dn,
         * Ad_up, * Ad_dn,
         *  x_up, *  x_dn,
         *  d_up, *  d_dn,
         * Ax_up, * Ax_dn;
  
  // formal parameters
  size_t dev_spinsize = 6*VOLUME/2*sizeof(dev_spinor);		// 24 floats per spinor per even lattice site
  int N_sites  =    VOLUME/2;					// Carsten's functions get the number of lattice points as input
  int N_floats = 24*VOLUME/2;
  
  // algorithm control parameters
  bool rbAx = true;						// choose how to calculate r(k+1)
  bool initial_guess = false;					// choose if initial guess
  
  
  
  
  //////////////////
  // INITIALIZING //
  //////////////////
  
  
  		//debug
  		printf("init_mixedsolve_eo_nd():\n");
  
  
    init_mixedsolve_eo_nd(g_gauge_field);			// initializes and allocates all quantities for the mixed solver
  								// more precise:
  								//	puts the gauge field on device as "2 rows" or "8 floats" per SU(3)-matrix
  								//	allocates memory for all spinor fields
  								//	puts the nn- and eoidx-fields on device memory

  		//debug
  		printf("mixedsolve_eo_nd():\n");
  
  
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
  		#endif
  
  
  he_cg_init_nd_additional<<<1,1>>> (g_mubar, g_epsbar);
  
  		// debug	// kernel
  		#ifdef CUDA_DEBUG
  		  CUDA_KERNEL_CHECK("Kernel error in he_cg_init_nd_additional(). Couldn't initialize some stuff.", "he_cg_init_nd_additional() succeeded.");
  		#endif
  
  		// debug	// check mubar and epsbar on host and device
  		#ifdef STUFF_DEBUG
  		  // printf("\tOn host:\n");
  		  // printf("\tg_mubar = %f\n", g_mubar);
  		  // printf("\tg_epsbar = %f\n", g_epsbar);
  		  
  		  float host_check_mubar, host_check_epsbar;
  		  cudaMemcpyFromSymbol(&host_check_mubar, mubar, sizeof(float));
  		  cudaMemcpyFromSymbol(&host_check_epsbar, epsbar, sizeof(float));
  		  printf("\tOn device:\n");
  		  printf("\tmubar = %f\n", host_check_mubar);
  		  printf("\tepsbar = %f\n", host_check_epsbar);
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
    benchmark_eo_nd(Q_up, Q_dn, OPERATOR_BENCHMARK);
  #endif
  
  
  
  
  /////////////////
  // ASSIGNMENTS //
  /////////////////
  
  
  x_up = P_up;							// can use the output spinors also as auxiliary fields
  x_dn = P_dn;							//	can use as initial guess at the same time
  
  
  #ifndef CG_DEBUG
  
  r_up  = g_chi_up_spinor_field[DUM_SOLVER];			// use the pre-allocated memory on host memory
  r_dn  = g_chi_dn_spinor_field[DUM_SOLVER];			// allocated by  init_chi_spinor_field.c  and  invert_doublet.c  !?
  d_up  = g_chi_up_spinor_field[DUM_SOLVER+1];			// the fields  g_chi_up/dn_spinor_field[DUM_SOLVER{ , +1, ... , +5}]  are used in  cg_her_nd()
  d_dn  = g_chi_dn_spinor_field[DUM_SOLVER+1];
  Ad_up = g_chi_up_spinor_field[DUM_SOLVER+2];
  Ad_dn = g_chi_dn_spinor_field[DUM_SOLVER+2];
  Ax_up = Ad_up;
  Ax_dn = Ad_dn;
  		// debug
  		printf("Now using the fields g_chi_up/dn_spinor_field[DUM_SOLVER{ , +1, +2}] in the mixedsolve_eo_nd().\n");
  
  #else
  
  		r_up  = (spinor *) malloc(24*N_sites*sizeof(double));		// if using cg_her_nd() as the CG, we cannot use the g_chi_up/dn-fields at the same time
  		r_dn  = (spinor *) malloc(24*N_sites*sizeof(double));
  		d_up  = (spinor *) malloc(24*N_sites*sizeof(double));
  		d_dn  = (spinor *) malloc(24*N_sites*sizeof(double));
  		Ad_up = (spinor *) malloc(24*N_sites*sizeof(double));
  		Ad_dn = (spinor *) malloc(24*N_sites*sizeof(double));
  		Ax_up = Ad_up;
  		Ax_dn = Ad_dn;
  				// debug
  				printf("Now allocating new host space for the fields in mixedsolve_eo_nd().\n");
  
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
    starteffective = ((double)clock()) / ((double)(CLOCKS_PER_SEC));
  #endif
  
  
  // r(0)
  if (!initial_guess) {		// r(0) = b = Q	// for x(0) = 0
    assign(r_up, Q_up, N_sites);
    assign(r_dn, Q_dn, N_sites);
    printf("x(0) = 0\n");
  }
  else {			// r(0) = b - A*x(0) = Q - A*P
    bb = square_norm(P_up, N_sites, 0) + square_norm(P_dn, N_sites, 0);
    		// benchmark
    		#ifdef CPU_BENCHMARK
    		  flopcount(host_flops, 2*2);
    		  // flopcount(host_flops, 2*2*N_floats);
    		#endif
    printf("bb = %.10e\n", bb);
    if (bb == 0) {
      assign(r_up, Q_up, N_sites);
      assign(r_dn, Q_dn, N_sites);
      printf("x(0) = 0\n");
    }
    else {
      Q_Qdagger_ND(Ax_up, Ax_dn, P_up, P_dn);
      diff(r_up, Q_up, Ax_up, N_sites);
      diff(r_dn, Q_dn, Ax_dn, N_sites);
      		// benchmark
      		#ifdef CPU_BENCHMARK
      		  flopcount(host_flops, 2*2*(55+2+2+1+55) + 2);
      		  // flopcount(host_flops, 2*2*(55+2+2+1+55)*N_floats + 2*N_floats);
      		#endif
      printf("x(0) != 0\n");
    }
  }
  
  
  // rr = (r_up)^2 + (r_dn)^2
  rr_up = square_norm(r_up, N_sites, 0);
  rr_dn = square_norm(r_dn, N_sites, 0);
  rr = rr_up + rr_dn;
  
  		// benchmark
  		#ifdef CPU_BENCHMARK
  		  flopcount(host_flops, 2*2);
  		  // flopcount(host_flops, 2*2*N_floats);
  		#endif
  
  
  r0r0   = rr; // for relative precision
  rr_old = rr; // for the first iteration
  
  		// debug
  		printf("Initial outer residue: %.10e\n", rr_old);
  
  
  // set to zero	// x_up, x_dn  will be added up		// as  x_up/dn = P_up/dn  up to here  P_up/dn  was not changed
  zero_spinor_field(x_up, N_sites);
  zero_spinor_field(x_dn, N_sites);
  
  
  
  
  ////////////////
  // OUTER LOOP //
  ////////////////
  
  		// debug
    		printf("\nEntering outer loop.");
  
  
  do {		// for (i = 0; i < max_iter; i++) {
    
    i++;
  
    		// debug
    		printf("\nouter iteration i = %i\n", i);
    
    
    
    
    #ifndef CG_DEBUG
    
    // host/device interaction
    to_device(dev_spinin_up, r_up, h2d_spin_up, dev_spinsize);
    to_device(dev_spinin_dn, r_dn, h2d_spin_dn, dev_spinsize);
    
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
    		printf("cg_eo_nd():\n");
    
    // solves A*p(k+1) = r(k)
    //        A*p(0)   = r(0) = b
    innercount = cg_eo_nd(dev_gf,
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
    		printf("Inner solver done in: %.4e sec\n", double(innerclocks) / double(CLOCKS_PER_SEC));
    
    
    // host/device interaction
    to_host(d_up, dev_spinout_up, h2d_spin_up, dev_spinsize);
    to_host(d_dn, dev_spinout_dn, h2d_spin_dn, dev_spinsize);
    
    		// debug	// CUDA
    		#ifdef CUDA_DEBUG
    		  // CUDA_CHECK("CUDA error in mixedsolve_eo_nd(). Device to host interaction failed.", "Fields copied back to device.");
    		  CUDA_CHECK_NO_SUCCESS_MSG("CUDA error in mixedsolve_eo_nd(). Device to host interaction failed.");
    		#endif
    
    
    #else
    
    
    				// debug
    				printf("cg_her_nd():\n");
    		
    		innercount = cg_her_nd(d_up, d_dn, r_up, r_dn,		// MISTAKE, was: r_up, r_dn, d_up, d_dn,
				           1000, eps_sq/2, 0,
				           VOLUME/2, &Q_Qdagger_ND, 0, 1000);
    		
    		outercount = outercount + innercount;
    		
    				// debug
    				printf("cg_her_nd() on host was used for debugging purposes.\n");
    
    
    #endif
    
    
    		// debug
    		printf("mixedsolve_eo_nd():\n");
    
    
    // x(k+1) = x(k) + d(k+1)
    add(x_up, x_up, d_up, N_sites);
    add(x_dn, x_dn, d_dn, N_sites);
    
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
      		printf("The matrix was applied on CPU in double precision. r = b - Ax\n");
      diff(r_up, Q_up, Ax_up, N_sites);
      diff(r_dn, Q_dn, Ax_dn, N_sites);
    }
    else {				// r(k+1) = r(k) - A*d(k+1)	// makes actually no sense ;)
      // A*d(k+1)
      Q_Qdagger_ND(Ad_up, Ad_dn, d_up, d_dn);
    		// debug
    		printf("The matrix was applied on CPU in double precision. r = r - Ad\n");
      // r(k+1) = r(k) - A*d(k+1)
      diff(r_up, r_up, Ad_up, N_sites);
      diff(r_dn, r_dn, Ad_dn, N_sites);
    }
    
    		// benchmark
  		#ifdef CPU_BENCHMARK
  		  flopcount(host_flops, 2*2*(55+2+2+1+55) + 2);
  		  // flopcount(host_flops, 2*2*(55+2+2+1+55)*N_floats + 2*N_floats);
  		#endif
    
    
    // rr = (rr_up)^2 + (r_dn)^2
    rr_up = square_norm(r_up, N_sites, 0);
    rr_dn = square_norm(r_dn, N_sites, 0);
    rr    = rr_up + rr_dn;
    
    		// debug
    		printf("Outer residue in the outer iteration i = %i after %i total inner iterations : %.10e\n", i, outercount, rr);
    		
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
        stopeffective = ((double)clock()) / ((double)(CLOCKS_PER_SEC));
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
      		printf("\nEO inversion done in mixed precision.\n");
      		if (rel_prec == 0) printf("Finished outer loop because of reached absolute outer solver precision.\n");
      		if (rel_prec == 1) printf("Finished outer loop because of reached relative outer solver precision.\n");
      		printf("Total number of inner iterations: %i\n", outercount);
      		printf("Total number of outer iterations: %i\n", i+1);
      		printf("Squared residue: %.10e\n", rr); 
      		printf("Outer solver done in: %.4e sec\n", double(stopouter-startouter) / double(CLOCKS_PER_SEC));
      		// benchmark
      		#ifdef EFFECTIVE_BENCHMARK
      		  // will now count the number of effective flops
      		  // effectiveflops  =  #(inner iterations)*(matrixflops+linalgflops)*VOLUME/2  +  #(outer iterations)*(matrixflops+linalgflops)*VOLUME/2
      		  // outer loop: linalg  =  flops for calculating  r(k+1) and x(k+1)
      		  // inner loop: linalg  =  flops for calculating  alpha, x(k+1), r(k+1), beta, d(k+1)
      		  effectiveflops = outercount*(matrixflops + 2*2*2*24 + 2*2*24 + 2*2*24 + 2*2*2*24 + 2*2*24)*VOLUME/2   +   i*(matrixflops + 2*24 + 2*24)*VOLUME/2;
      		  printf("effective BENCHMARK:\n");
      		  printf("\ttotal mixed solver time:   %.2e sec\n", double(stopeffective-starteffective));
      		  printf("\tfloating point operations: %.2e flops\n", effectiveflops);
      		  printf("\tinner solver performance:  %.2e Gflop/s\n", double(effectiveflops) / double(stopeffective-starteffective) / 1.0e9);
      		#endif
      		
      		#if defined(GPU_BENCHMARK) || defined(GPU_BENCHMARK2)
      		  // REMARK: device_flops has to be multiplied by N_floats !!
      		  flops = device_flops * N_floats / (double(totalinnerclocks)/double(CLOCKS_PER_SEC)) / 1.0e9;
      		  printf("Inner solver BENCHMARK:\n");
      		  printf("\ttotal inner solver time:   %.2e sec\n", double(totalinnerclocks) / double(CLOCKS_PER_SEC));
      		  printf("\tfloating point operations: %.2e flops\n", double(device_flops) * double(N_floats));
      		  printf("\tinner solver performance:  %.2e Gflop/s\n", flops);
      		#endif
      		
      		#ifdef CPU_BENCHMARK
      		  // REMARK: host_flops has to be multiplied by N_floats !!
      		  flops = host_flops * N_floats / (double(totalouterclocks)/double(CLOCKS_PER_SEC)) / 1.0e9;
      		  printf("Outer solver BENCHMARK:\n");
      		  printf("\ttotal outer solver time:   %.2e sec\n", double(totalouterclocks) / double(CLOCKS_PER_SEC));
      		  printf("\tfloating point operations: %.2e flops\n", double(host_flops) * double(N_floats));
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
      		printf("finalize_mixedsolve_eo_nd():\n");
      
      finalize_mixedsolve_eo_nd();
      
      		// debug
      		printf("\n");
      
      return(outercount);
    }
    
    
  }//OUTER LOOP
  while (outercount <= max_iter);
  
  
  // multiplying with Qhat(2x2)^dagger is done in invert_doublet_eo.c
  
  
  // timer
  stopouter = clock();
  totalouterclocks = stopouter-startouter - totalinnerclocks;
  
  #ifdef EFFECTIVE_BENCHMARK
    stopeffective = ((double)clock()) / ((double)(CLOCKS_PER_SEC));
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
  		printf("\nEO inversion done in mixed precision.\n");
  		printf("Finished outer loop, because of maximal number of outer iterations.\n");
      		printf("Total number of inner iterations: %i\n", outercount);
      		printf("Total number of outer iterations: %i\n", i+1);
      		printf("Squared residue: %.10e\n", rr); 
      		printf("Outer solver done in: %.4e sec\n", double(stopouter-startouter)/CLOCKS_PER_SEC);
      		// benchmark
      		#ifdef EFFECTIVE_BENCHMARK
      		  // will now count the number of effective flops
      		  // effectiveflops  =  #(inner iterations)*(matrixflops+linalgflops)*VOLUME/2  +  #(outer iterations)*(matrixflops+linalgflops)*VOLUME/2
      		  // outer loop: linalg  =  flops for calculating  r(k+1) and x(k+1)
      		  // inner loop: linalg  =  flops for calculating  alpha, x(k+1), r(k+1), beta, d(k+1)
      		  effectiveflops = outercount*(matrixflops + 2*2*2*24 + 2*2*24 + 2*2*24 + 2*2*2*24 + 2*2*24)*VOLUME/2   +   i*(matrixflops + 2*24 + 2*24)*VOLUME/2;
      		  printf("effective BENCHMARK:\n");
      		  printf("\ttotal mixed solver time:   %.2e sec\n", double(stopeffective-starteffective));
      		  printf("\tfloating point operations: %.2e flops\n", effectiveflops);
      		  printf("\tinner solver performance:  %.2e Gflop/s\n", double(effectiveflops) / double(stopeffective-starteffective) / 1.0e9);
      		#endif
      		
      		#if defined(GPU_BENCHMARK) || defined(GPU_BENCHMARK2)
      		  // REMARK: device_flops has to be multiplied by N_floats !!
      		  flops = device_flops * N_floats / (double(totalinnerclocks)/double(CLOCKS_PER_SEC)) / 1.0e9;
      		  printf("Inner solver BENCHMARK:\n");
      		  printf("\ttotal inner solver time:   %.2e sec\n", double(totalinnerclocks) / double(CLOCKS_PER_SEC));
      		  printf("\tfloating point operations: %.2e flops\n", double(device_flops) * double(N_floats));
      		  printf("\tinner solver performance:  %.2e Gflop/s\n", flops);
      		#endif
      		
      		#ifdef CPU_BENCHMARK
      		  // REMARK: host_flops has to be multiplied by N_floats !!
      		  flops = host_flops * N_floats / (double(totalouterclocks)/double(CLOCKS_PER_SEC)) / 1.0e9;
      		  printf("Outer solver BENCHMARK:\n");
      		  printf("\ttotal outer solver time:   %.2e sec\n", double(totalouterclocks) / double(CLOCKS_PER_SEC));
      		  printf("\tfloating point operations: %.2e flops\n", double(host_flops) * double(N_floats));
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
  		printf("finalize_mixedsolve_eo_nd():\n");  
  
  finalize_mixedsolve_eo_nd();
  
  		// debug
  		printf("\n");
  
  return(outercount);
  
  
}//mixedsolve_eo_nd()



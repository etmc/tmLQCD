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
	
	//////////////////////////////////////////////////////////////////
	//								//
	//    and the MPI implementation of the EO, ND mixed solver	//
	//								//
	//	PARALLELT parallelization				//
	//	no _GAUGE_COPY and no _USE_HALFSPINOR			//
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
#include "../solver/solver_field.h"
}
#include "../global.h"

#ifdef MPI
  #undef MPI
  #undef REAL
    #include <mpi.h>
  #define MPI
  #define REAL float
#endif







// global formal parameters
size_t dev_gfsize;
size_t dev_spinsize_int;		// making the structure transparent:							
int N_sites_int;			// _int: internal sites
int N_floats_int;			// _ext: internal sites + additional boundaries
#ifdef MPI
  size_t dev_spinsize_ext;
  int N_sites_ext;
  int N_floats_ext;
#endif


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


#ifdef MPI					// collecting variables for the MPI implementation
  						// put to mixed_solve.cu
  /*
  __device__ int dev_RAND;			// not used, maybe later ...
  __device__ int dev_VOLUMEPLUSRAND;		// is now used in dev_Hopping_Matrix_mpi()
  __device__ int dev_rank;			// was for the moment put to mixed_solve.cu ...
  __device__ int dev_nproc;
  */
  
  int * iseven;
  int * dev_g_iup;
  int * dev_g_idn;
  int * dev_g_lexic2eo;
  int * dev_g_lexic2eosub;
  int * dev_g_eo2lexic;
  int * dev_g_ipt;
  
  #ifdef HOPPING_DEBUG
    spinor * spinor_debug_in;			// for Hopping_Matrix_wrapper()
    spinor * spinor_debug_out;			// for Hopping_Matrix_wrapper()
  #endif
  
  
  #if ASYNC > 0
    #ifdef ASYNC_TIMING
      cudaEvent_t start_ALL;			// CUDA events for timing and profiling
      cudaEvent_t stop_ALL;
      cudaEvent_t stop_D2H_1;
      cudaEvent_t stop_D2H_2;
      cudaEvent_t stop_INT_0;
      cudaEvent_t stop_H2D_3;
      cudaEvent_t stop_H2D_4;
      cudaEvent_t stop_EXT_1;
      cudaEvent_t stop_EXT_2;
      float time_stop_D2H_1;			// CUDA times in milliseconds
      float time_stop_D2H_2;
      float time_stop_INT_0;
      float time_stop_H2D_3;
      float time_stop_H2D_4;
      float time_stop_EXT_1;
      float time_stop_EXT_2;
      float time_stop_ALL;
      double mpi_start_ALL;			// MPI times with arbitrary zero-point for timing and profiling
      double mpi_start_sendrecv_1;
      double mpi_stop_sendrecv_1;
      double mpi_start_sendrecv_2;
      double mpi_stop_sendrecv_2;
      double mpiTime_start_sendrecv_1;		// MPI times in seconds
      double mpiTime_stop_sendrecv_1;
      double mpiTime_start_sendrecv_2;
      double mpiTime_stop_sendrecv_2;
    #endif
  #endif
   
#endif





//#include "communication.cuh"
//#include "index_fields.cuh"






///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////






		///////////////////////////
		//                       //
		//    INITIALIZATIONS    //
		//                       //
		///////////////////////////




////////////////////
// GPU parameters //
////////////////////


// puts the additional nd parameters mubar and epsbar on the device
__global__ void he_cg_init_nd_additional (float param_mubar, float param_epsbar) {

  mubar  = param_mubar;
  epsbar = param_epsbar;

}






#ifdef MPI

// puts the additional variables VOLUMEPLUSRAND and RAND on the device
__global__ void he_cg_init_nd_additional_mpi (int param_VOLUMEPLUSRAND, int param_RAND, int rank, int nproc) {

  dev_VOLUMEPLUSRAND  = param_VOLUMEPLUSRAND;
  dev_RAND            = param_RAND;
  
  dev_rank            = rank;
  dev_nproc           = nproc;

}

#endif






/////////////////////////////////////////////
// geometry- and nearest-neighbour indices //
/////////////////////////////////////////////


#ifdef MPI

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


#endif	// MPI






void set_global_sizes() {
  
  #ifndef MPI
  	#ifdef GF_8
  	  // allocate 8 floats for gf = 2*4*VOLUME float4's			// dev_su3_8 = float4
  	  dev_gfsize = 4*VOLUME * 2*sizeof(dev_su3_8);				// allocates for each lattice site and for 4 directions  2*float4 = 8 floats  = 8 real parameters
  	#else
  	  // allocate 2 rows of gf = 3*4*VOLUME float4's			// dev_su3_2v = float4
  	  dev_gfsize = 4*VOLUME * 3*sizeof(dev_su3_2v); 			// allocates for each lattice site and for 4 directions  3*float4 = 12 floats = 2 rows of complex 3-vectors
  	#endif
  #else
  	#ifdef GF_8								// dev_su3_8 = float4
  	  dev_gfsize = 4*(VOLUME+RAND) * 2*sizeof(dev_su3_8);			// allocates for each lattice site and RAND for 4 directions  2*float4 = 8 floats  = 8 real parameters
  	#else									// dev_su3_2v = float4
  	  dev_gfsize = 4*(VOLUME+RAND) * 3*sizeof(dev_su3_2v); 			// allocates for each lattice site and RAND for 4 directions  3*float4 = 12 floats = 2 rows of complex 3-vectors
  	#endif
  #endif
  
  dev_spinsize_int   =  6*VOLUME/2*sizeof(dev_spinor);				// 24 floats per lattice site
  N_sites_int        =    VOLUME/2;
  N_floats_int       = 24*VOLUME/2;
  #ifdef MPI
    dev_spinsize_ext =  6*(VOLUME+RAND)/2*sizeof(dev_spinor);
    N_sites_ext      =    (VOLUME+RAND)/2;
    N_floats_ext     = 24*(VOLUME+RAND)/2;
  #endif
  
}






////////////////
// ALLOCATING //
////////////////

// initializes and allocates all quantities for the mixed solver
// more precise:
//	puts the gauge field on device as "2 rows" or "8 floats" per SU(3)-matrix
//	allocates memory for all spinor fields
//	puts the nn- and eoidx-fields on device memory

void init_mixedsolve_eo_nd (su3** gf) {	// gf is the full gauge field
  
  
  
  
  typedef REAL RealT;

  //////////////////////
  // GLOBAL VARIABLES //
  //////////////////////
  
  /*
  #ifndef MPI
  	#ifdef GF_8
  	  // allocate 8 floats for gf = 2*4*VOLUME float4's			// dev_su3_8 = float4
  	  dev_gfsize = 4*VOLUME * 2*sizeof(dev_su3_8);				// allocates for each lattice site and for 4 directions  2*float4 = 8 floats  = 8 real parameters
  	#else
  	  // allocate 2 rows of gf = 3*4*VOLUME float4's			// dev_su3_2v = float4
  	  dev_gfsize = 4*VOLUME * 3*sizeof(dev_su3_2v); 			// allocates for each lattice site and for 4 directions  3*float4 = 12 floats = 2 rows of complex 3-vectors
  	#endif
  #else
  	#ifdef GF_8								// dev_su3_8 = float4
  	  dev_gfsize = 4*(VOLUME+RAND) * 2*sizeof(dev_su3_8);			// allocates for each lattice site and RAND for 4 directions  2*float4 = 8 floats  = 8 real parameters
  	#else									// dev_su3_2v = float4
  	  dev_gfsize = 4*(VOLUME+RAND) * 3*sizeof(dev_su3_2v); 			// allocates for each lattice site and RAND for 4 directions  3*float4 = 12 floats = 2 rows of complex 3-vectors
  	#endif
  #endif
  
  dev_spinsize_int   =  6*VOLUME/2*sizeof(dev_spinor);				// 24 floats per lattice site
  N_sites_int        =    VOLUME/2;
  N_floats_int       = 24*VOLUME/2;
  #ifdef MPI
    dev_spinsize_ext =  6*(VOLUME+RAND)/2*sizeof(dev_spinor);
    N_sites_ext      =    (VOLUME+RAND)/2;
    N_floats_ext     = 24*(VOLUME+RAND)/2;
  #endif
  */
  
  set_global_sizes();
  
  
  
  /////////////////////
  // LOCAL VARIABLES //
  /////////////////////
  
  cudaError_t cudaerr;		// CUDA errors
  int ndev;			// number of devices
  //size_t dev_gfsize;		// size of the gauge-field on device memory	// put to global	// non-MPI: VOLUME/2	// MPI: (VOLUME+RAND)/2
  size_t nnsize;		// size of memory for nn-table
  size_t idxsize;		// size of memory for even/odd-positions
  //size_t dev_spinsize;	// size of memory for spinors			// put to global
  int grid[5];			// array for grid specifications
  float * host_output;		// ??
  
  
  
  
  // get number of devices
  
  if (havedevice == 0) {
  
  ndev = find_devices();
  if (ndev == 0) {
    fprintf(stderr, "Error: no CUDA devices found. Aborting...\n");
    exit(300);
  }
  
  #ifndef MPI
      // only if device_num is not the default (-1)
      if(device_num > -1){ 
    	// try to set active device to device_num given in input file
    	if (device_num < ndev) {
    	  printf("Setting active device to: %d\n", device_num);
    	  cudaSetDevice(device_num);
    	}
    	else {
    	  fprintf(stderr, "Error: There is no CUDA device with No. %d. Aborting...\n",device_num);
    	  exit(301);
    	}
    	if ((cudaerr=cudaGetLastError())!=cudaSuccess) {
    	  printf("Error in init_mixedsolve_eo_nd(): Could not set active device. Aborting...\n");
    	  exit(302);
    	}
      }
      else{
        printf("Not setting any active device. Let the driver choose.\n");
      }   
  #else
    	#ifndef DEVICE_EQUAL_RANK
    	  // try to set active device to device_num given in input file
    	  // each process gets bounded to the same GPU
    	  if (device_num < ndev) {
    	    printf("Process %d of %d: Setting active device to: %d\n", g_proc_id, g_nproc, device_num);
    	    //cudaSetDevice(device_num);
    	  }
    	  else {
    	    fprintf(stderr, "Process %d of %d: Error: There is no CUDA device with No. %d. Aborting...\n", g_proc_id, g_nproc, device_num);
    	    exit(301);
    	  }
  	#else
    	  // device number = mpi rank
    	  if (g_cart_id < ndev) {
    	    printf("Process %d of %d: Setting active device to: %d\n", g_proc_id, g_nproc, g_cart_id);
    	    //cudaSetDevice(g_cart_id);
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
  #endif
  
  havedevice = 1;
  }
  
  
  
  
  // output
  #ifdef MPI
    if (g_cart_id == 0) {
  #endif
  
  	#ifdef USETEXTURE
  	  printf("Using texture references.\n");
  	#else
  	  printf("NOT using texture references.\n");
  	#endif
  
  	#ifdef GF_8
 	  printf("Using GF 8 reconstruction.\n");
  	#else
  	  printf("Using GF 12 reconstruction.\n");
  	#endif
  
  #ifdef MPI
    }
  #endif
  
  
  
  
  /////////////////
  // GAUGE FIELD //
  /////////////////
  
  /*									// put to global
  #ifndef MPI
  	#ifdef GF_8
  	  // allocate 8 floats for gf = 2*4*VOLUME float4's		// dev_su3_8 = float4
  	  dev_gfsize = 4*VOLUME * 2*sizeof(dev_su3_8);			// allocates for each lattice site and for 4 directions  2*float4 = 8 floats  = 8 real parameters
  	#else
  	  // allocate 2 rows of gf = 3*4*VOLUME float4's		// dev_su3_2v = float4
  	  dev_gfsize = 4*VOLUME * 3*sizeof(dev_su3_2v); 		// allocates for each lattice site and for 4 directions  3*float4 = 12 floats = 2 rows of complex 3-vectors
  	#endif
  #else
  	#ifdef GF_8							// dev_su3_8 = float4
  	  dev_gfsize = 4*(VOLUME+RAND) * 2*sizeof(dev_su3_8);		// allocates for each lattice site and RAND for 4 directions  2*float4 = 8 floats  = 8 real parameters
  	#else								// dev_su3_2v = float4
  	  dev_gfsize = 4*(VOLUME+RAND) * 3*sizeof(dev_su3_2v); 		// allocates for each lattice site and RAND for 4 directions  3*float4 = 12 floats = 2 rows of complex 3-vectors
  	#endif
  #endif
  */
  
  
  if ( (cudaerr = cudaMalloc((void **) &MixedsolveParameter<RealT>::getGlobalP()->dev_gf, dev_gfsize)) != cudaSuccess ) {	// allocates memory for the gauge field MixedsolveParameter<RealT>::getGlobalP()->dev_gf on device
    printf("Error in init_mixedsolve_eo_nd(): Memory allocation of gauge field failed. Aborting...\n");
    printf("%s\n", cudaGetErrorString(cudaGetLastError()));
    exit(200);
  }
  else {
    #ifndef MPI
      printf("Allocated memory for gauge field on device.\n");
    #else
      if (g_cart_id == 0) printf("Allocated memory for gauge gauge field on devices.\n");
    #endif
  }
  
  
  #ifdef GF_8
    MixedsolveParameter<RealT>::getGlobalP()->h2d_gf = (dev_su3_8 *) malloc(dev_gfsize); 			// allocates on host
    su3to8<RealT>(gf, MixedsolveParameter<RealT>::getGlobalP()->h2d_gf);						// MixedsolveParameter<RealT>::getGlobalP()->h2d_gf  is the gauge field  gf  with the 8-real-parameter-representation (according to M. Clark, p. 28)
  #else
    MixedsolveParameter<RealT>::getGlobalP()->h2d_gf = (dev_su3_2v *) malloc(dev_gfsize);			// allocates on host
    su3to2vf4<RealT>(gf, MixedsolveParameter<RealT>::getGlobalP()->h2d_gf);					// MixedsolveParameter<RealT>::getGlobalP()->h2d_gf  is the gauge field  gf  with the first two rows stored
  #endif
  
  cudaMemcpy(MixedsolveParameter<RealT>::getGlobalP()->dev_gf, MixedsolveParameter<RealT>::getGlobalP()->h2d_gf, dev_gfsize, cudaMemcpyHostToDevice);
  								// MixedsolveParameter<RealT>::getGlobalP()->dev_gf = MixedsolveParameter<RealT>::getGlobalP()->h2d_gf  on device memory
  
  		// debug	// CUDA
  		#ifdef CUDA_DEBUG
  		  #ifndef MPI
  		    CUDA_CHECK("CUDA error in init_mixedsolve_eo_nd(). Copying MixedsolveParameter<RealT>::getGlobalP()->dev_gf to device failed.", "Copied MixedsolveParameter<RealT>::getGlobalP()->dev_gf to device.");
  		  #else
  		    CUDA_CHECK("CUDA error in init_mixedsolve_eo_nd(). Copying MixedsolveParameter<RealT>::getGlobalP()->dev_gf to device failed.", "Copied MixedsolveParameter<RealT>::getGlobalP()->dev_gf to devices.");
  		  #endif
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
  
  
  #ifndef MPI
    idxsize = VOLUME/2*sizeof(int);				// size of memory necessary for VOLUME/2 integers
  #else
    idxsize = (VOLUME+RAND)/2*sizeof(int);
  #endif
  eoidx_even = (int *) malloc(idxsize);				// allocate on host
  eoidx_odd = (int *) malloc(idxsize);				// allocate on host
  cudaMalloc((void **) &dev_eoidx_even, idxsize);		// allocate on device
  cudaMalloc((void **) &dev_eoidx_odd, idxsize);		// allocate on device
  
  
  #ifndef MPI
    initnn();							// initialize nearest-neighbour table for gpu
    initnn_eo();						// initialize nearest-neighbour table for gpu with even-odd enabled
  #else
    init_nnspinor_eo_mpi();					// initialize nearest-neighbour table for gpu with even-odd enabled
    init_idxgauge_mpi();
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
  		  #ifndef MPI
  		    CUDA_CHECK("CUDA error in init_mixedsolve_eo_nd(). Memory allocation of grid stuff failed.", "Allocated grid stuff on device.");
  		  #else
  		    CUDA_CHECK("CUDA error in init_mixedsolve_eo_nd(). Memory allocation of grid stuff failed.", "Allocated grid stuff on devices.");
  		  #endif
  		#endif
  
  
  
  
  /////////////
  // SPINORS //							// allocates device memory for the odd part of the spinor fields (24 = 4*6 floats per odd lattice sites)
  /////////////							// now we have to consider 2 flavors: up, dn
  
  /*
  #ifndef MPI
    dev_spinsize = 6*VOLUME/2*sizeof(dev_spinor);		// remember: dev_spinor = float4
  #else
    dev_spinsize = (VOLUME+RAND)/2 * 6*sizeof(dev_spinor);	// NOTICE: this refers to the memory requirements for the device, host needs twice the memory !!
  #endif
  */
  
  
  #ifndef MPI
  
    cudaMalloc((void **) &dev_spin1_up, dev_spinsize_int);   	// allocates device memory for the fields spinor fields used in dev_cg_eo_nd(...)
    cudaMalloc((void **) &dev_spin1_dn, dev_spinsize_int);	// pointing to device
    cudaMalloc((void **) &dev_spin2_up, dev_spinsize_int);	// ...
    cudaMalloc((void **) &dev_spin2_dn, dev_spinsize_int);
    cudaMalloc((void **) &dev_spin3_up, dev_spinsize_int);
    cudaMalloc((void **) &dev_spin3_dn, dev_spinsize_int);
    /*
    cudaMalloc((void **) &dev_spin4_up, dev_spinsize_int);	// not needed
    cudaMalloc((void **) &dev_spin4_dn, dev_spinsize_int);
    cudaMalloc((void **) &dev_spin5_up, dev_spinsize_int);
    cudaMalloc((void **) &dev_spin5_dn, dev_spinsize_int);
    */
    cudaMalloc((void **) &dev_spinin_up , dev_spinsize_int);	// host/device interaction
    cudaMalloc((void **) &dev_spinin_dn , dev_spinsize_int);	// inner/outer interaction
    cudaMalloc((void **) &dev_spinout_up, dev_spinsize_int);
    cudaMalloc((void **) &dev_spinout_dn, dev_spinsize_int);
  
  #else
  
    cudaMalloc((void **) &dev_spin1_up, dev_spinsize_ext);
    cudaMalloc((void **) &dev_spin1_dn, dev_spinsize_ext);
    cudaMalloc((void **) &dev_spin2_up, dev_spinsize_ext);
    cudaMalloc((void **) &dev_spin2_dn, dev_spinsize_ext);
    cudaMalloc((void **) &dev_spin3_up, dev_spinsize_ext);
    cudaMalloc((void **) &dev_spin3_dn, dev_spinsize_ext);
    /*
    cudaMalloc((void **) &dev_spin4_up, dev_spinsize_ext);
    cudaMalloc((void **) &dev_spin4_dn, dev_spinsize_ext);
    cudaMalloc((void **) &dev_spin5_up, dev_spinsize_ext);
    cudaMalloc((void **) &dev_spin5_dn, dev_spinsize_ext);
    */
    cudaMalloc((void **) &dev_spinin_up , dev_spinsize_ext);
    cudaMalloc((void **) &dev_spinin_dn , dev_spinsize_ext);
    cudaMalloc((void **) &dev_spinout_up, dev_spinsize_ext);
    cudaMalloc((void **) &dev_spinout_dn, dev_spinsize_ext);
  
  #endif
  
  
  #ifndef MPI
  		// debug	// host code
  		if ( (void *) (h2d_spin_up = (dev_spinor *) malloc(dev_spinsize_int) ) == NULL) {
  		  printf("Could not allocate memory for h2d_spin_up. Aborting...\n");
  		  exit(200);
  		}
  		
  		if ( (void *) (h2d_spin_dn = (dev_spinor *) malloc(dev_spinsize_int) ) == NULL) {
  		  printf("Could not allocate memory for h2d_spin_dn. Aborting...\n");
  		  exit(200);
  		}
  #else
  		// debug	// host code
  		if ( (void *) (h2d_spin_up = (dev_spinor *) malloc(dev_spinsize_ext) ) == NULL) {				// MEMORY REQUIREMENTS: these are auxiliary fields for  to_host()  and  to_device()
  		  printf("Process %d of %d: Could not allocate memory for h2d_spin_up. Aborting...\n", g_proc_id, g_nproc);	//                      they have to store floats (not doubles)
  		  exit(200);													//			can use "_int" ...
  		}														//			must use "_ext" when used with to_host_mpi as in xchange_field_wrapper()
  		
  		if ( (void *) (h2d_spin_dn = (dev_spinor *) malloc(dev_spinsize_ext) ) == NULL) {
  		  printf("Process %d of %d: Could not allocate memory for h2d_spin_dn. Aborting...\n", g_proc_id, g_nproc);
  		  exit(200);
  		}
  #endif
  
  
  #ifndef MPI
  
    cudaMalloc((void **) &dev_spin_eo1_up, dev_spinsize_int);		// used for matrix_multiplication32(...)
    cudaMalloc((void **) &dev_spin_eo1_dn, dev_spinsize_int);
    cudaMalloc((void **) &dev_spin_eo3_up, dev_spinsize_int);
    cudaMalloc((void **) &dev_spin_eo3_dn, dev_spinsize_int);
    /*
    cudaMalloc((void **) &dev_spin_eo2_up, dev_spinsize_int);		// no memory allocation needed
    cudaMalloc((void **) &dev_spin_eo2_dn, dev_spinsize_int);		// will point to already allocated memory when used in matrix_multiplication
    */
  
  #else
  
    cudaMalloc((void **) &dev_spin_eo1_up, dev_spinsize_ext);
    cudaMalloc((void **) &dev_spin_eo1_dn, dev_spinsize_ext);
    cudaMalloc((void **) &dev_spin_eo3_up, dev_spinsize_ext);
    cudaMalloc((void **) &dev_spin_eo3_dn, dev_spinsize_ext);
    /*
    cudaMalloc((void **) &dev_spin_eo2_up, dev_spinsize_ext);
    cudaMalloc((void **) &dev_spin_eo2_dn, dev_spinsize_ext);
    */
  
  #endif
  
  		// debug	// CUDA
  		#ifdef CUDA_DEBUG
  		  #ifndef MPI
  		    CUDA_CHECK("CUDA error in init_mixedsolve_eo_nd(). Memory allocation of spinor fields failed.", "Allocated spinor fields on device.");
  		  #else
  		    CUDA_CHECK("CUDA error in init_mixedsolve_eo_nd(). Memory allocation of spinor fields failed.", "Allocated spinor fields on devices.");
  		  #endif
  		#endif
  
  
  
  
  
  #ifdef MPI
  
  	#ifdef HOPPING_DEBUG													// Hopping_Matrix() is applied upon these spinor fields
  		// debug	// host code
  		if ( (void *) (spinor_debug_in = (spinor *) malloc(2*dev_spinsize_ext) ) == NULL) {
  		  printf("Process %d of %d: Could not allocate memory for spinor_debug_in. Aborting...\n", g_proc_id, g_nproc);
  		  exit(200);
  		}
  		// debug	// host code
  		if ( (void *) (spinor_debug_out = (spinor *) malloc(2*dev_spinsize_ext) ) == NULL) {
  		  printf("Process %d of %d: Could not allocate memory for spinor_debug_out. Aborting...\n", g_proc_id, g_nproc);
  		  exit(200);
  		}
  	#endif
  
  
  
  
  	#if defined(ALTERNATE_FIELD_XCHANGE) || ASYNC > 0
  	  int tSliceEO = LX*LY*LZ/2;
  	#endif
        
        
  	#ifndef ALTERNATE_FIELD_XCHANGE												// xchange_field() acts on this spinor field
  		// debug	// host code											// MEMORY REQUIREMENTS:
  		if ( (void *) (spinor_xchange = (spinor *) malloc(2*dev_spinsize_ext) ) == NULL) {				//	auxiliary fields for  xchange_field_wrapper()  and  Hopping_Matrix_wrapper()
  		  printf("Process %d of %d: Could not allocate memory for spinor_xchange. Aborting...\n", g_proc_id, g_nproc);	//	have to store doubles --> 2*dev_spinsize  !!
  		  exit(200);
  		}
  	#else		// xchange procedure comparable to ASYNC
  		R1 = (dev_spinor *) malloc(2*tSliceEO*24*sizeof(float));
  		R2 = R1 + 6*tSliceEO;
  		R3 = (dev_spinor *) malloc(2*tSliceEO*24*sizeof(float));
  		R4 = R3 + 6*tSliceEO;
  	#endif
  
  
  	#if ASYNC > 0	// asynchronous communication and computation
  	
  	  // page-locked memory
  	  cudaMallocHost(&RAND3, 2*tSliceEO*6*sizeof(float4));
  	  RAND4 = RAND3 + 6*tSliceEO;
  	  cudaMallocHost(&RAND1, 2*tSliceEO*6*sizeof(float4));
  	  RAND2 = RAND1 + 6*tSliceEO;
  	  
  	  // CUDA streams and events
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
  
  #endif	// MPI
  
  
  
  
  ////////////
  // output //						// ??
  ////////////
  /*
  output_size = LZ*T*sizeof(float); 			// parallel in t and z direction
  cudaMalloc((void **) &dev_output, output_size);	// output array
  host_output = (float *) malloc(output_size);
  
  
  		// debug	// CUDA
  		#ifdef CUDA_DEBUG
  		  #ifndef MPI
  		    CUDA_CHECK("CUDA error in init_mixedsolve_eo_nd(). Memory allocation output stuff failed.", "Allocated output stuff on device.");
  		  #else
  		    CUDA_CHECK("CUDA error in init_mixedsolve_eo_nd(). Memory allocation output stuff failed.", "Allocated output stuff on devices.");
  		  #endif
  		#endif
  */
  
  
  
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
  		  #ifndef MPI
  		    CUDA_CHECK("CUDA error in init_mixedsolve_eo_nd(). Memory allocation of grid[] specifications failed.", "Allocated grid[] specifications on device.");
  		  #else
  		    CUDA_CHECK("CUDA error in init_mixedsolve_eo_nd(). Memory allocation of grid[] specifications failed.", "Allocated grid[] specifications on devices.");
  		  #endif
  		#endif
  
  
  
  
  // MPI_Barrier(g_cart_grid);
  
  
  
  
}//init_mixedsolve_eo_nd()






		////////////////////////
		//                    //
		//    FINALIZATION    //
		//                    //
		////////////////////////




// deallocates the previous allocated memory

void finalize_mixedsolve_eo_nd(void) {

  typedef REAL RealT; 
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
  
  #ifdef MPI
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
  

  cudaFree(MixedsolveParameter<RealT>::getGlobalP()->dev_gf);
  //cudaFree(dev_output);
  cudaFree(dev_grid);
  
  
  free(MixedsolveParameter<RealT>::getGlobalP()->h2d_gf);
  
  
  #ifdef MPI
  	#ifdef ALTERNATE_HOPPING_MATRIX
  	  free_gpu_indexfields();
  	#endif
  	
  	#if ASYNC > 0
  	  cudaFreeHost(RAND1);
  	  cudaFreeHost(RAND3);
  	  
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
  #endif
  
  
  // Clean up CUDA API for calling thread	// ??
  cudaThreadExit();				// is essential
  
  
  		// debug	// CUDA
  		#ifdef CUDA_DEBUG
  		  CUDA_CHECK("CUDA error in finalize_mixedsolve_eo_nd(). Device memory deallocation failed", "Device memory deallocated.");
  		#endif
  
  
}






		/////////////////////////////////
		//                             //
		//    H <--> D interactions    //
		//                             //
		/////////////////////////////////




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

#endif // MPI






/////////////////////////////
// host/device interaction //
/////////////////////////////

// remark: the host spinors are double precision and therefore need twice the memory !!
//		dev_spinor * device:    dev_spinsize
//		spinor * host:        2*dev_spinsize
//		dev_spinor * auxiliary: dev_spinsize
//         the parameter "size" specifies the memory needed for the spinor n the device !!
//         

void to_device (dev_spinor * device, spinor * host, dev_spinor * auxiliary, int size) {

  typedef REAL RealT;
  convert2REAL4_spin<RealT>(host, auxiliary);						// auxiliary = (float) host
  cudaMemcpy(device, auxiliary, size, cudaMemcpyHostToDevice);			// device = auxiliary  (on device)

}


void to_host (spinor * host, dev_spinor * device, dev_spinor * auxiliary, int size) {

  typedef REAL RealT;
  cudaMemcpy(auxiliary, device, size, cudaMemcpyDeviceToHost);			// auxiliary = device  (on device)
  convert2double_spin<RealT>(auxiliary, host);						// host = (double) auxiliary

}






///////////////////////
// boundary exchange //
///////////////////////

#ifdef MPI

// all three versions do work:

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




// copies the boundary t-slices t=0 and t=T-1 to host		// will be used in matrix_multiplication32_mpi(), not ASYNC
//	exchanges						// provides a wrapped version of Carsten's xchange_field()
//		copies RAND back to device			//	and not asynchronous version of ASYNC.cuh

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

#endif // MPI






////////////////////
// hopping matrix //
////////////////////

#ifdef MPI	// implemented for checking the MPI implementation of the hopping matrix
  #ifdef HOPPING_DEBUG

  // applies the hopping matrix on host for debugging purposes
  
  void Hopping_Matrix_wrapper (int ieo, dev_spinor * out, dev_spinor * in) {
  
      //size_t size = (VOLUME+RAND)/2 * 6*sizeof(dev_spinor);
      //to_host(g_chi_up_spinor_field[DUM_OLVER+3], in, h2d_spin_up, size);
      //Hopping_Matrix(ieo, g_chi_dn_spinor_field[DUM_OLVER+3], g_chi_up_spinor_field[DUM_OLVER+3]);
      //to_device(out, g_chi_dn_spinor_field[DUM_OLVER+3], h2d_spin_up, size);
      
      to_host(spinor_debug_in, in, h2d_spin_up, dev_spinsize_int);
      Hopping_Matrix(ieo, spinor_debug_out, spinor_debug_in);
      to_device(out, spinor_debug_out, h2d_spin_dn, dev_spinsize_int);    
    
  }

  #endif
#endif






////////////////////
// linear algebra //
////////////////////

#ifdef MPI

// have to rebuilt some linear algebra functions which contain global communication
// can be done as wrappers to appropriate CUBLAS routines



// a wrapper function for cublasDot() (with the same interface)
// provides the MPI communication via MPI_Allreduce()

float cublasDot_wrapper(int size, float * A, int incx, float * B, int incy) {

  float result;
  float buffer;
  
  buffer = cublasDot(size, (float *) A, incx, (float *) B, incy);
  MPI_Allreduce(&buffer, &result, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
  
  return(result);
  
}

#endif






		//////////////////////////////////
		//                              //
		//    MATRIX MULTIPLICATIONS    //
		//                              //
		//////////////////////////////////




/////////////
// KERNELS //
/////////////

// derived from Flo's function  dev_mul_one_pm_imu_inv
//	order of the arguments also like Flo's convention: (spinin, spinout)

// applies (1 +- imubar*gamma5)
// uses shared local memory for manipulation	// really ??	where ??
// one thread per lattice site


__global__ void dev_mul_one_pm_imubar_gamma5 (dev_spinor * sin,
                                              dev_spinor * sout,
                                              float sign         ) {
   
  dev_spinor slocal[6];									// dev_spinor = float4		// 6*float4 = 24 floats		// auxiliary for each thread
  
  dev_complex pm_imu = dev_initcomplex<REAL>(0.0, sign * mubar);				// dev_complex = struct { REAL re; REAL im; }	// pm_imu.re = 0.0
  																	// pm_imu.im = sign * mubar
  int pos = threadIdx.x + blockDim.x*blockIdx.x;
  
  if (pos < dev_VOLUME) {
    dev_skalarmult_gamma5_spinor(&(slocal[0]), pm_imu, &(sin[6*pos]) );			// slocal  =  pm_imu * (gamma5) * sin
    dev_add_spinor_assign<REAL>(&(slocal[0]), &(sin[6*pos]));					// slocal  =  slocal + sin  =  pm_imu * (gamma5) * sin + sin
    dev_realmult_spinor_assign<REAL>(&(sout[6*pos]), 1.0, &(slocal[0]) );			// sout    =  slocal
  }
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
  typedef REAL RealT;
  
  
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
    dev_Hopping_Matrix<RealT> <<<gridsize1, blocksize1>>>(MixedsolveParameter<RealT>::getGlobalP()->dev_gf, spinin_dn, dev_spin_eo1_up, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0);	// dev_spin_eo1_up = (M_eo) * spinin_dn
  #ifdef USETEXTURE
    unbind_texture_spin(1);
  #endif
  
  #ifdef USETEXTURE
    bind_texture_spin(spinin_up,1);
  #endif
    dev_Hopping_Matrix<RealT> <<<gridsize1, blocksize1>>>(MixedsolveParameter<RealT>::getGlobalP()->dev_gf, spinin_up, dev_spin_eo1_dn, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0);	// dev_spin_eo1_dn = (M_eo) * spinin_up
  #ifdef USETEXTURE
    unbind_texture_spin(1);
  #endif
  
  
  // written:
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo1_up, dev_spin_eo2_up, -1.0);		// dev_spin_eo2_up  =  (1 - imubar) * dev_spin_eo1_up  =  (1 - imubar)*(M_eo) * spinin_dn
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo1_dn, dev_spin_eo2_dn, +1.0);		// dev_spin_eo2_dn  =  (1 + imubar) * dev_spin_eo1_dn  =  (1 + imubar)*(M_eo) * spinin_up
  
  
  // CUBLAS:
  cublasAxpy (N_floats, g_epsbar, (RealT*)dev_spin_eo1_dn, 1, (RealT*)dev_spin_eo2_up, 1);
  										// dev_spin_eo2_up  =  dev_spin_eo2_up  +  epsbar * dev_spin_eo1_dn  =  (1 - imubar)*(M_eo) * spinin_dn  +  epsbar * (M_eo) * spinin_up
  cublasAxpy (N_floats, g_epsbar, (RealT*)dev_spin_eo1_up, 1, (RealT*)dev_spin_eo2_dn, 1);
  										// dev_spin_eo2_dn  =  dev_spin_eo2_dn  +  epsbar * dev_spin_eo1_up  =  (1 + imubar)*(M_eo) * spinin_up  +  epsbar * (M_eo) * spinin_dn
  
  // CUBLAS:
  cublasScal (N_floats, nrm, (RealT*)dev_spin_eo2_up, 1);			// dev_spin_eo2_up  =  nrm * dev_spin_eo2_up  =  nrm*(1-imubar)*(M_eo) * spinin_dn  +  nrm*epsbar*(M_eo) * spinin_up
  
  cublasScal (N_floats, nrm, (RealT*)dev_spin_eo2_dn, 1);			// dev_spin_eo2_dn  =  nrm * dev_spin_eo2_dn  =  nrm*(1+imubar)*(M_eo) * spinin_up  +  nrm*epsbar*(M_eo) * spinin_dn
  
  
  // Flo:
  #ifdef USETEXTURE
    bind_texture_spin(dev_spin_eo2_up,1);									// remember: this is ((M_oe)(Mee^-1)(M_eo)) * (spinin_dn, spinin_up):
  #endif
    dev_Hopping_Matrix<RealT> <<<gridsize1, blocksize1>>>(MixedsolveParameter<RealT>::getGlobalP()->dev_gf, dev_spin_eo2_up, dev_spin_eo1_up, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1);
  #ifdef USETEXTURE
    unbind_texture_spin(1);							// dev_spin_eo1_up  =  (M_oe) * dev_spin_eo2_up  =  (M_oe)*nrm*(1-imubar)*(M_eo) * spinin_dn  +  (M_oe)*nrm*epsbar*(M_eo) * spinin_up
  #endif									
  
  #ifdef USETEXTURE
    bind_texture_spin(dev_spin_eo2_dn,1);
  #endif
    dev_Hopping_Matrix<RealT> <<<gridsize1, blocksize1>>>(MixedsolveParameter<RealT>::getGlobalP()->dev_gf, dev_spin_eo2_dn, dev_spin_eo1_dn, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1);
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
  cublasAxpy (N_floats, -g_epsbar, (RealT*)spinin_up, 1, (RealT*)dev_spin_eo2_up, 1);
  												// dev_spin_eo2_up  =  dev_spin_eo2_up - epsbar*spinin_up  =  (1+imubar)*spinin_dn - epsbar*spinin_up
  cublasAxpy (N_floats, -g_epsbar, (RealT*)spinin_dn, 1, (RealT*)dev_spin_eo2_dn, 1);
  												// dev_spin_eo2_dn  =  dev_spin_eo2_dn - epsbar*spinin_dn  =  (1-imubar)*spinin_up - epsbar*spinin_dn
  
  
  
  ///////////////////////////////////////
  // (M_oo)  -  (M_oe) (Mee^-1) (M_eo) //
  ///////////////////////////////////////
  
  // CUBLAS:													// this is ((M_oo) - (M_oe)(Mee^-1)(M_eo)) * (spinin_dn, spinin_up):
  cublasAxpy (N_floats, -1.0, (RealT*)dev_spin_eo1_up, 1, (RealT*)dev_spin_eo2_up, 1);
  							// dev_spin_eo2_up  =  dev_spin_eo2_up  -  dev_spin_eo1_up  =                      (1+imubar) * spinin_dn  -                    epsbar * spinin_up
  							//                                                             - (M_oe)*nrm*(1-imubar)*(M_eo) * spinin_dn  -  (M_oe)*nrm*epsbar*(M_eo) * spinin_up
  cublasAxpy (N_floats, -1.0, (RealT*)dev_spin_eo1_dn, 1, (RealT*)dev_spin_eo2_dn, 1);
  							// dev_spin_eo2_dn  =  dev_spin_eo2_dn  -  dev_spin_eo1_dn  =                      (1-imubar) * spinin_up  -                    epsbar * spinin_dn
  							//                                                             - (M_oe)*nrm*(1+imubar)*(M_eo) * spinin_up  -  (M_oe)*nrm*epsbar*(M_eo) * spinin_dn
  
  
  ////////////
  // gamma5 //
  ////////////
  
  // Flo:
  dev_gamma5<RealT> <<<gridsize3, blocksize3>>>(dev_spin_eo2_up, dev_spin_eo2_up);					// dev_spin_eo2_up = gamma5 * dev_spin_eo2_up 
  dev_gamma5<RealT> <<<gridsize3, blocksize3>>>(dev_spin_eo2_dn, dev_spin_eo2_dn);					// dev_spin_eo2_dn = gamma5 * dev_spin_eo2_dn
  	
  
  
  
  
  
  ////////////////////
  // (a,b) -> (b,a) //
  ////////////////////
  dev_copy_spinor_field<RealT,RealT> <<<gridsize4, blocksize4>>>(dev_spin_eo2_dn, dev_spin_eo3_up);			// dev_spin_eo3_up = dev_spin_eo2_dn
  dev_copy_spinor_field<RealT,RealT> <<<gridsize4, blocksize4>>>(dev_spin_eo2_up, dev_spin_eo3_dn);			// dev_spin_eo3_dn = dev_spin_eo2_up
  
  
  
  
  
  
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
    dev_Hopping_Matrix<RealT> <<<gridsize1, blocksize1>>>(MixedsolveParameter<RealT>::getGlobalP()->dev_gf, dev_spin_eo3_up, dev_spin_eo1_up, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0);	// dev_spin_eo1_up = (M_eo) * dev_spin_eo3_up
  #ifdef USETEXTURE
    unbind_texture_spin(1);
  #endif
  
  #ifdef USETEXTURE
    bind_texture_spin(dev_spin_eo3_dn,1);
  #endif
    dev_Hopping_Matrix<RealT> <<<gridsize1, blocksize1>>>(MixedsolveParameter<RealT>::getGlobalP()->dev_gf, dev_spin_eo3_dn, dev_spin_eo1_dn, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0);	// dev_spin_eo1_dn = (M_eo) * dev_spin_eo3_dn
  #ifdef USETEXTURE
    unbind_texture_spin(1);
  #endif
  
  
  // written:
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo1_up, dev_spin_eo2_up, -1.0);			// dev_spin_eo2_up  =  (1 - imubar) * dev_spin_eo1_up  =  (1 - imubar)*(M_eo) * dev_spin_eo3_up
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo1_dn, dev_spin_eo2_dn, +1.0);			// dev_spin_eo2_dn  =  (1 + imubar) * dev_spin_eo1_dn  =  (1 + imubar)*(M_eo) * dev_spin_eo3_dn
  
  
  // CUBLAS:
  cublasAxpy (N_floats, g_epsbar, (RealT*)dev_spin_eo1_dn, 1, (RealT*)dev_spin_eo2_up, 1);
  										// dev_spin_eo2_up  =  dev_spin_eo2_up  +  epsbar * dev_spin_eo1_dn  =  (1 - imubar)*(M_eo) * dev_spin_eo3_up  +  epsbar * (M_eo) * dev_spin_eo3_dn
  cublasAxpy (N_floats, g_epsbar, (RealT*)dev_spin_eo1_up, 1, (RealT*)dev_spin_eo2_dn, 1);
  										// dev_spin_eo2_dn  =  dev_spin_eo2_dn  +  epsbar * dev_spin_eo1_up  =  (1 + imubar)*(M_eo) * dev_spin_eo3_dn  +  epsbar * (M_eo) * dev_spin_eo3_up
  
  // CUBLAS:
  cublasScal (N_floats, nrm, (RealT*)dev_spin_eo2_up, 1);			// dev_spin_eo2_up  =  nrm * dev_spin_eo2_up  =  nrm*(1-imubar)*(M_eo) * dev_spin_eo3_up  +  nrm*epsbar*(M_eo) * dev_spin_eo3_dn
  
  cublasScal (N_floats, nrm, (RealT*)dev_spin_eo2_dn, 1);			// dev_spin_eo2_dn  =  nrm * dev_spin_eo2_dn  =  nrm*(1+imubar)*(M_eo) * dev_spin_eo3_dn  +  nrm*epsbar*(M_eo) * dev_spin_eo3_up
  
  
  // Flo:
  #ifdef USETEXTURE
    bind_texture_spin(dev_spin_eo2_up,1);									// remember: this is ((M_oe) (Mee^-1) (M_eo)) * (dev_spin_eo3_up, dev_spin_eo3_dn):
  #endif
    dev_Hopping_Matrix<RealT> <<<gridsize1, blocksize1>>>(MixedsolveParameter<RealT>::getGlobalP()->dev_gf, dev_spin_eo2_up, dev_spin_eo1_up, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1);
  #ifdef USETEXTURE
    unbind_texture_spin(1);							// dev_spin_eo1_up  =  (M_oe) * dev_spin_eo2_up  =  (M_oe)*nrm*(1-imubar)*(M_eo) * dev_spin_eo3_up  +  (M_oe)*nrm*epsbar*(M_eo) * dev_spin_eo3_dn
  #endif
  
  #ifdef USETEXTURE
    bind_texture_spin(dev_spin_eo2_dn,1);
  #endif
    dev_Hopping_Matrix<RealT> <<<gridsize1, blocksize1>>>(MixedsolveParameter<RealT>::getGlobalP()->dev_gf, dev_spin_eo2_dn, dev_spin_eo1_dn, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1);
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
  cublasAxpy (N_floats, -g_epsbar, (RealT*)dev_spin_eo3_dn, 1, (RealT*)dev_spin_eo2_up, 1);
  												// dev_spin_eo2_up  =  dev_spin_eo2_up - epsbar*dev_spin_eo3_dn  =  (1+imubar)*dev_spin_eo3_up - epsbar*dev_spin_eo3_dn
  cublasAxpy (N_floats, -g_epsbar, (RealT*)dev_spin_eo3_up, 1, (RealT*)dev_spin_eo2_dn, 1);
  												// dev_spin_eo2_dn  =  dev_spin_eo2_dn - epsbar*dev_spin_eo3_up  =  (1-imubar)*dev_spin_eo3_dn - epsbar*dev_spin_eo3_up
  
  
  
  ///////////////////////////////////////
  // (M_oo)  -  (M_oe) (Mee^-1) (M_eo) //
  ///////////////////////////////////////
  
  // CUBLAS:											// this is ( (M_oo) - (M_oe) (Mee^-1) (M_eo) ) * (dev_spin_eo3_up, dev_spin_eo3_dn)
  cublasAxpy (N_floats, -1.0, (RealT*)dev_spin_eo1_up, 1, (RealT*)dev_spin_eo2_up, 1);
  							// dev_spin_eo2_up  =  dev_spin_eo2_up  -  dev_spin_eo1_up  =                      (1+imubar) * dev_spin_eo3_up  -                    epsbar * dev_spin_eo3_dn
  							//                                                             - (M_oe)*nrm*(1-imubar)*(M_eo) * dev_spin_eo3_up  -  (M_oe)*nrm*epsbar*(M_eo) * dev_spin_eo3_dn
  cublasAxpy (N_floats, -1.0, (RealT*)dev_spin_eo1_dn, 1, (RealT*)dev_spin_eo2_dn, 1);
  							// dev_spin_eo2_dn  =  dev_spin_eo2_dn  -  dev_spin_eo1_dn  =                      (1-imubar) * dev_spin_eo3_dn  -                    epsbar * dev_spin_eo3_up
  							//                                                             - (M_oe)*nrm*(1+imubar)*(M_eo) * dev_spin_eo3_dn  -  (M_oe)*nrm*epsbar*(M_eo) * dev_spin_eo3_up
  
  
  ////////////
  // gamma5 //
  ////////////
  
  // Flo:
  dev_gamma5<RealT> <<<gridsize3, blocksize3>>>(dev_spin_eo2_up, dev_spin_eo2_up);					// dev_spin_eo2_up = gamma5 * dev_spin_eo2_up 
  dev_gamma5<RealT> <<<gridsize3, blocksize3>>>(dev_spin_eo2_dn, dev_spin_eo2_dn);					// dev_spin_eo2_dn = gamma5 * dev_spin_eo2_dn
  
  
  
  
  /*
  ////////////
  // output //										// output is already done by setting  dev_spin_eo2_up/dn = spinout_up/dn
  ////////////
  
  dev_copy_spinor_field<<<griddim2, blockdim2 >>>(dev_spin_eo2_up, spinout_up);		// spinout_up = dev_spin_eo2_up
  dev_copy_spinor_field<<<griddim2, blockdim2 >>>(dev_spin_eo2_dn, spinout_dn);		// spinout_dn = dev_spin_eo2_dn
  */
  
  
  return;
  
}//matrix_multiplication32()






#ifdef MPI

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
    	  dev_Hopping_Matrix<RealT> <<<gridsize1, blocksize1>>>(MixedsolveParameter<RealT>::getGlobalP()->dev_gf, spinin_dn, dev_spin_eo1_up, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0);	// dev_spin_eo1_up = (M_eo) * spinin_dn
    	#else
    	  dev_Hopping_Matrix_alternate<<<gridsize1, blocksize1>>>(MixedsolveParameter<RealT>::getGlobalP()->dev_gf, spinin_dn, dev_spin_eo1_up, dev_g_iup, dev_g_idn, dev_g_eo2lexic, dev_g_lexic2eosub, 0);
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
    	  dev_Hopping_Matrix<RealT> <<<gridsize1, blocksize1>>>(MixedsolveParameter<RealT>::getGlobalP()->dev_gf, spinin_up, dev_spin_eo1_dn, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0);	// dev_spin_eo1_dn = (M_eo) * spinin_up
    	#else
    	  dev_Hopping_Matrix_alternate<<<gridsize1, blocksize1>>>(MixedsolveParameter<RealT>::getGlobalP()->dev_gf, spinin_up, dev_spin_eo1_dn, dev_g_iup, dev_g_idn, dev_g_eo2lexic, dev_g_lexic2eosub, 0);
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
  cublasAxpy (N_floats, g_epsbar, (RealT*)dev_spin_eo1_dn, 1, (RealT*)dev_spin_eo2_up, 1);
  										// dev_spin_eo2_up  =  dev_spin_eo2_up  +  epsbar * dev_spin_eo1_dn  =  (1 - imubar)*(M_eo) * spinin_dn  +  epsbar * (M_eo) * spinin_up
  cublasAxpy (N_floats, g_epsbar, (RealT*)dev_spin_eo1_up, 1, (RealT*)dev_spin_eo2_dn, 1);
  										// dev_spin_eo2_dn  =  dev_spin_eo2_dn  +  epsbar * dev_spin_eo1_up  =  (1 + imubar)*(M_eo) * spinin_up  +  epsbar * (M_eo) * spinin_dn
  
  // linear algebra
  cublasScal (N_floats, nrm, (RealT*)dev_spin_eo2_up, 1);			// dev_spin_eo2_up  =  nrm * dev_spin_eo2_up  =  nrm*(1-imubar)*(M_eo) * spinin_dn  +  nrm*epsbar*(M_eo) * spinin_up
  
  cublasScal (N_floats, nrm, (RealT*)dev_spin_eo2_dn, 1);			// dev_spin_eo2_dn  =  nrm * dev_spin_eo2_dn  =  nrm*(1+imubar)*(M_eo) * spinin_up  +  nrm*epsbar*(M_eo) * spinin_dn
  
  
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
    	  dev_Hopping_Matrix<RealT> <<<gridsize1, blocksize1>>>(MixedsolveParameter<RealT>::getGlobalP()->dev_gf, dev_spin_eo2_up, dev_spin_eo1_up, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1);
    	#else
    	  dev_Hopping_Matrix_alternate<<<gridsize1, blocksize1>>>(MixedsolveParameter<RealT>::getGlobalP()->dev_gf, dev_spin_eo2_up, dev_spin_eo1_up, dev_g_iup, dev_g_idn, dev_g_eo2lexic, dev_g_lexic2eosub, 1);
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
    	  dev_Hopping_Matrix<RealT> <<<gridsize1, blocksize1>>>(MixedsolveParameter<RealT>::getGlobalP()->dev_gf, dev_spin_eo2_dn, dev_spin_eo1_dn, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1);
    	#else
    	  dev_Hopping_Matrix_alternate<<<gridsize1, blocksize1>>>(MixedsolveParameter<RealT>::getGlobalP()->dev_gf, dev_spin_eo2_dn, dev_spin_eo1_dn, dev_g_iup, dev_g_idn, dev_g_eo2lexic, dev_g_lexic2eosub, 1);
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
  cublasAxpy (N_floats, -g_epsbar, (RealT*)spinin_up, 1, (RealT*)dev_spin_eo2_up, 1);
  												// dev_spin_eo2_up  =  dev_spin_eo2_up - epsbar*spinin_up  =  (1+imubar)*spinin_dn - epsbar*spinin_up
  cublasAxpy (N_floats, -g_epsbar, (RealT*)spinin_dn, 1, (RealT*)dev_spin_eo2_dn, 1);
  												// dev_spin_eo2_dn  =  dev_spin_eo2_dn - epsbar*spinin_dn  =  (1-imubar)*spinin_up - epsbar*spinin_dn
  
  
  
  ///////////////////////////////////////
  // (M_oo)  -  (M_oe) (Mee^-1) (M_eo) //
  ///////////////////////////////////////
  
  // linear algebra													// this is ((M_oo) - (M_oe)(Mee^-1)(M_eo)) * (spinin_dn, spinin_up):
  cublasAxpy (N_floats, -1.0, (RealT*)dev_spin_eo1_up, 1, (RealT*)dev_spin_eo2_up, 1);
  							// dev_spin_eo2_up  =  dev_spin_eo2_up  -  dev_spin_eo1_up  =                      (1+imubar) * spinin_dn  -                    epsbar * spinin_up
  							//                                                             - (M_oe)*nrm*(1-imubar)*(M_eo) * spinin_dn  -  (M_oe)*nrm*epsbar*(M_eo) * spinin_up
  cublasAxpy (N_floats, -1.0, (RealT*)dev_spin_eo1_dn, 1, (RealT*)dev_spin_eo2_dn, 1);
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
    	  dev_Hopping_Matrix<RealT> <<<gridsize1, blocksize1>>>(MixedsolveParameter<RealT>::getGlobalP()->dev_gf, dev_spin_eo3_up, dev_spin_eo1_up, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0);	// dev_spin_eo1_up = (M_eo) * dev_spin_eo3_up
    	#else
    	  dev_Hopping_Matrix_alternate<<<gridsize1, blocksize1>>>(MixedsolveParameter<RealT>::getGlobalP()->dev_gf, dev_spin_eo3_up, dev_spin_eo1_up, dev_g_iup, dev_g_idn, dev_g_eo2lexic, dev_g_lexic2eosub, 0);
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
    	  dev_Hopping_Matrix<RealT> <<<gridsize1, blocksize1>>>(MixedsolveParameter<RealT>::getGlobalP()->dev_gf, dev_spin_eo3_dn, dev_spin_eo1_dn, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0);	// dev_spin_eo1_dn = (M_eo) * dev_spin_eo3_dn
    	#else
    	  dev_Hopping_Matrix_alternate<<<gridsize1, blocksize1>>>(MixedsolveParameter<RealT>::getGlobalP()->dev_gf, dev_spin_eo3_dn, dev_spin_eo1_dn, dev_g_iup, dev_g_idn, dev_g_eo2lexic, dev_g_lexic2eosub, 0);
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
  cublasAxpy (N_floats, g_epsbar, (RealT*)dev_spin_eo1_dn, 1, (RealT*)dev_spin_eo2_up, 1);
  										// dev_spin_eo2_up  =  dev_spin_eo2_up  +  epsbar * dev_spin_eo1_dn  =  (1 - imubar)*(M_eo) * dev_spin_eo3_up  +  epsbar * (M_eo) * dev_spin_eo3_dn
  cublasAxpy (N_floats, g_epsbar, (RealT*)dev_spin_eo1_up, 1, (RealT*)dev_spin_eo2_dn, 1);
  										// dev_spin_eo2_dn  =  dev_spin_eo2_dn  +  epsbar * dev_spin_eo1_up  =  (1 + imubar)*(M_eo) * dev_spin_eo3_dn  +  epsbar * (M_eo) * dev_spin_eo3_up
  
  // lineare algebra
  cublasScal (N_floats, nrm, (RealT*)dev_spin_eo2_up, 1);			// dev_spin_eo2_up  =  nrm * dev_spin_eo2_up  =  nrm*(1-imubar)*(M_eo) * dev_spin_eo3_up  +  nrm*epsbar*(M_eo) * dev_spin_eo3_dn
  
  cublasScal (N_floats, nrm, (RealT*)dev_spin_eo2_dn, 1);			// dev_spin_eo2_dn  =  nrm * dev_spin_eo2_dn  =  nrm*(1+imubar)*(M_eo) * dev_spin_eo3_dn  +  nrm*epsbar*(M_eo) * dev_spin_eo3_up
  
  
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
    	  dev_Hopping_Matrix<RealT> <<<gridsize1, blocksize1>>>(MixedsolveParameter<RealT>::getGlobalP()->dev_gf, dev_spin_eo2_up, dev_spin_eo1_up, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1);
    	#else
    	  dev_Hopping_Matrix_alternate<<<gridsize1, blocksize1>>>(MixedsolveParameter<RealT>::getGlobalP()->dev_gf, dev_spin_eo2_up, dev_spin_eo1_up, dev_g_iup, dev_g_idn, dev_g_eo2lexic, dev_g_lexic2eosub, 1);
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
    	  dev_Hopping_Matrix<RealT> <<<gridsize1, blocksize1>>>(MixedsolveParameter<RealT>::getGlobalP()->dev_gf, dev_spin_eo2_dn, dev_spin_eo1_dn, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1);
    	#else
    	  dev_Hopping_Matrix_alternate<<<gridsize1, blocksize1>>>(MixedsolveParameter<RealT>::getGlobalP()->dev_gf, dev_spin_eo2_dn, dev_spin_eo1_dn, dev_g_iup, dev_g_idn, dev_g_eo2lexic, dev_g_lexic2eosub, 1);
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
  cublasAxpy (N_floats, -g_epsbar, (RealT*)dev_spin_eo3_dn, 1, (RealT*)dev_spin_eo2_up, 1);
  												// dev_spin_eo2_up  =  dev_spin_eo2_up - epsbar*dev_spin_eo3_dn  =  (1+imubar)*dev_spin_eo3_up - epsbar*dev_spin_eo3_dn
  cublasAxpy (N_floats, -g_epsbar, (RealT*)dev_spin_eo3_up, 1, (RealT*)dev_spin_eo2_dn, 1);
  												// dev_spin_eo2_dn  =  dev_spin_eo2_dn - epsbar*dev_spin_eo3_up  =  (1-imubar)*dev_spin_eo3_dn - epsbar*dev_spin_eo3_up
  
  
  
  ///////////////////////////////////////
  // (M_oo)  -  (M_oe) (Mee^-1) (M_eo) //
  ///////////////////////////////////////
  
  // lineare algebra										// this is ( (M_oo) - (M_oe) (Mee^-1) (M_eo) ) * (dev_spin_eo3_up, dev_spin_eo3_dn)
  cublasAxpy (N_floats, -1.0, (RealT*)dev_spin_eo1_up, 1, (RealT*)dev_spin_eo2_up, 1);
  							// dev_spin_eo2_up  =  dev_spin_eo2_up  -  dev_spin_eo1_up  =                      (1+imubar) * dev_spin_eo3_up  -                    epsbar * dev_spin_eo3_dn
  							//                                                             - (M_oe)*nrm*(1-imubar)*(M_eo) * dev_spin_eo3_up  -  (M_oe)*nrm*epsbar*(M_eo) * dev_spin_eo3_dn
  cublasAxpy (N_floats, -1.0, (RealT*)dev_spin_eo1_dn, 1, (RealT*)dev_spin_eo2_dn, 1);
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


#endif	// MPI






		/////////////////////
		//                 //
		//    BENCHMARK    //
		//                 //
		/////////////////////




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
  #ifndef MPI
    double timeElapsed;
  #else
    double singleTimeElapsed;
    double maxTimeElapsed;
  #endif
  double startBenchmark;
  double stopBenchmark;
  
  // counter
  int i;
  
  // flop counting
  /*
  double realFlopsPerApp = 34768.0;
  */
  // double effectiveFlopsPerApp = 23984.0;	// hopping = 1488
  double effectiveFlopsPerApp = 21296.0;	// per lattice site
  
  #ifndef MPI
    /*
    double realDeviceFlops;
    double realFlops;
    */
    double effectiveDeviceFlops;
    double effectiveFlops;
  #else
    /*
    double realDeviceFlops;
    double allRealDeviceFlops;
    double realFlops;
    */
    double effectiveDeviceFlops;
    double allEffectiveDeviceFlops;
    double effectiveFlops;
  #endif
  
  // CUDA errors
  cudaError_t cudaerr;
  cublasStatus cublasstatus;
  
  // size of a spinor
  /*
  size_t dev_spinsize_int = 6*VOLUME/2 * sizeof(dev_spinor);
  #ifdef MPI
    size_t dev_spinsize_ext = 6*(VOLUME+RAND)/2 * sizeof(dev_spinor);
  #endif
  */
  
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
  
  #ifndef MPI
    cudaMalloc((void **) &A_up, dev_spinsize_int);
    cudaMalloc((void **) &A_dn, dev_spinsize_int);
    cudaMalloc((void **) &B_up, dev_spinsize_int);
    cudaMalloc((void **) &B_dn, dev_spinsize_int);
  #else
    cudaMalloc((void **) &A_up, dev_spinsize_ext);
    cudaMalloc((void **) &A_dn, dev_spinsize_ext);
    cudaMalloc((void **) &B_up, dev_spinsize_ext);
    cudaMalloc((void **) &B_dn, dev_spinsize_ext);
  #endif
  
  		// debug	// CUDA
  		#ifdef CUDA_DEBUG
  		  CUDA_CHECK_NO_SUCCESS_MSG("CUDA error in benchmark_eo_nd(). Memory allocation of spinor fields failed.");
  		#endif
  
  
  /*
  #ifdef USETEXTURE
    bind_texture_gf(MixedsolveParameter<RealT>::getGlobalP()->dev_gf);
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
  		#ifndef MPI
  		  printf("\nStarting a little BENCHMARK. benchmark_eo_nd().\n");
  		#else
  		  if (g_proc_id == 0) printf("\nStarting a little BENCHMARK. benchmark_eo_nd_mpi().\n");
  		#endif
  
  
  
  
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
  		#ifndef MPI
  		  printf("Applying the eo-preconditioned matrix %i times.\n", N);
  		#else
  		  if (g_proc_id == 0) printf("Applying the eo-preconditioned matrix %i times.\n", N);
  		#endif
  
  
  to_device(B_up, Q_up, h2d_spin_up, dev_spinsize_int);
  to_device(B_dn, Q_dn, h2d_spin_dn, dev_spinsize_int);
  
  
  // timer
  #ifndef MPI
    startBenchmark = double(clock()) / double(CLOCKS_PER_SEC);
  #else
    startBenchmark = MPI_Wtime();
  #endif
  
  
  
  
  for (i = 0; i < N; i++) {
  
  
    #ifndef MPI
    	matrix_multiplication32(A_up, A_dn,					// A = (matrix)*B
    	                        B_up, B_dn,
    	                        griddim2, blockdim2,
    	                        griddim3, blockdim3,
    	                        griddim4, blockdim4,
    	                        griddim5, blockdim5);
    #else
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
    #endif
    
    
    if (staticsource == 0) {
      // swaps A and B
      C_up = B_up;
      C_dn = B_dn;
      B_up = A_up;
      B_dn = A_dn;
      A_up = C_up;
      A_dn = C_dn;
    }
    //else {
      // do nothing
    //}
    
  }
  
  		// debug	// CUDA
  		#ifdef CUDA_DEBUG
  		  CUDA_CHECK_NO_SUCCESS_MSG("CUDA error in matrix_muliplication32(). Applying the matrix on GPU failed.");
  		#endif
  
  
  
  // timer
  #ifndef MPI
    stopBenchmark = double(clock()) / double(CLOCKS_PER_SEC);
  #else
    stopBenchmark = MPI_Wtime();
  #endif
  
  
  #ifndef MPI
  
  	timeElapsed = stopBenchmark - startBenchmark;
  	/*
  	realDeviceFlops      = N * VOLUME/2 * realFlopsPerApp;
  	realFlops            = N * VOLUME/2 * realFlopsPerApp / timeElapsed / 1.0e9;
  	*/
  	effectiveDeviceFlops = N * VOLUME/2 * effectiveFlopsPerApp;
  	effectiveFlops       = N * VOLUME/2 * effectiveFlopsPerApp / timeElapsed / 1.0e9;
  	
  	/*
  	printf("REAL:\n");
  	printf("\ttime:        %.2e sec\n", timeElapsed);
  	printf("\tflop's:      %.2e flops\n", realDeviceFlops);
  	printf("\tperformance: %.2e Gflop/s\n\n", realFlops);
  	*/
  	printf("EFFECTIVE:\n");
  	printf("\ttime:        %.4e sec\n", timeElapsed);
  	printf("\tflop's:      %.4e flops\n", effectiveDeviceFlops);
  	printf("\tperformance: %.4e Gflop/s\n\n", effectiveFlops);
  	
  #else
  	
  	singleTimeElapsed = stopBenchmark - startBenchmark;
  	MPI_Allreduce(&singleTimeElapsed, &maxTimeElapsed, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  	/*
  	realDeviceFlops      = N * VOLUME/2 * realFlopsPerApp;
  	MPI_Allreduce(&realDeviceFlops, &allRealDeviceFlops, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  	realFlops            = allRealDeviceFlops / maxTimeElapsed / 1.0e9;
  	*/
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
  	  printf("\ttime:        %.4e sec\n", maxTimeElapsed);
  	  printf("\tflop's:      %.4e flops\n", allEffectiveDeviceFlops);
  	  printf("\tperformance: %.4e Gflop/s\n\n", effectiveFlops);
  	  
  	  #if ASYNC > 0 && defined(ASYNC_TIMING)
  	    // calculate the times from the "beginning"
  	    cudaEventElapsedTime(&time_stop_D2H_1, start_ALL, stop_D2H_1);
  	    cudaEventElapsedTime(&time_stop_D2H_2, start_ALL, stop_D2H_2);
  	    cudaEventElapsedTime(&time_stop_INT_0, start_ALL, stop_INT_0);
  	    cudaEventElapsedTime(&time_stop_H2D_3, start_ALL, stop_H2D_3);
  	    cudaEventElapsedTime(&time_stop_H2D_4, start_ALL, stop_H2D_4);
  	    cudaEventElapsedTime(&time_stop_EXT_1, start_ALL, stop_EXT_1);
  	    cudaEventElapsedTime(&time_stop_EXT_2, start_ALL, stop_EXT_2);
  	    cudaEventElapsedTime(&time_stop_ALL  , start_ALL, stop_ALL);
  	    mpiTime_start_sendrecv_1 = mpi_start_sendrecv_1 - mpi_start_ALL; 
  	    mpiTime_stop_sendrecv_1  = mpi_stop_sendrecv_1  - mpi_start_ALL;
  	    mpiTime_start_sendrecv_2 = mpi_start_sendrecv_2 - mpi_start_ALL;
  	    mpiTime_stop_sendrecv_2  = mpi_stop_sendrecv_2  - mpi_start_ALL;
  	    // outputting the times
  	    #if ASYNC == 1
  	      printf("\tTIMING[sec]:\n");
  	      printf("\tSTART:        %.2e   -       \n", 0.0);
  	      printf("\tD2H_1:                   -   %.2e\n", time_stop_D2H_1/1000);
  	      printf("\tINT_0:                   -   %.2e\n", time_stop_INT_0/1000);
  	      printf("\tSENDRECV_1:   %.2e   -   %.2e\n", mpiTime_start_sendrecv_1, mpiTime_stop_sendrecv_1);
  	      printf("\tH2D_3:        %.2e   -   %.2e\n", mpiTime_stop_sendrecv_1, time_stop_H2D_3/1000);
  	      printf("\tEXT_1:        %.2e   -   %.2e\n", time_stop_H2D_3/1000, time_stop_EXT_1/1000);
  	      printf("\tD2H_2:        %.2e   -   %.2e\n", mpiTime_stop_sendrecv_1, time_stop_D2H_2/1000);
  	      printf("\tSENDRECV_2:   %.2e   -   %.2e\n", mpiTime_start_sendrecv_2, mpiTime_stop_sendrecv_2);
  	      printf("\tH2D_4:        %.2e   -   %.2e\n", mpiTime_stop_sendrecv_2, time_stop_H2D_4/1000);
  	      printf("\tEXT_2:        %.2e   -   %.2e\n", time_stop_H2D_4/1000, time_stop_EXT_2/1000);
  	      printf("\tSTOP:                    -   %.2e\n", time_stop_ALL/1000);
  	    #elif ASYNC == 2
  	      printf("\tTIMING[sec]:\n");
  	      printf("\tSTART:        %.2e   -       \n", 0.0);
  	      printf("\tD2H_1:                   -   %.2e\n", time_stop_D2H_1/1000);
  	      printf("\tD2H_2:                   -   %.2e\n", time_stop_D2H_2/1000);
  	      printf("\tINT_0:                   -   %.2e\n", time_stop_INT_0/1000);
  	      printf("\tSENDRECV_1:   %.2e   -   %.2e\n", mpiTime_start_sendrecv_1, mpiTime_stop_sendrecv_1);
  	      printf("\tH2D_3:        %.2e   -   %.2e\n", mpiTime_stop_sendrecv_1, time_stop_H2D_3/1000);
  	      printf("\tEXT_1:        %.2e   -   %.2e\n", time_stop_H2D_3/1000, time_stop_EXT_1/1000);
  	      printf("\tSENDRECV_2:   %.2e   -   %.2e\n", mpiTime_start_sendrecv_2, mpiTime_stop_sendrecv_2);
  	      printf("\tH2D_4:        %.2e   -   %.2e\n", mpiTime_stop_sendrecv_2, time_stop_H2D_4/1000);
  	      printf("\tEXT_2:        %.2e   -   %.2e\n", time_stop_H2D_4/1000, time_stop_EXT_2/1000);
  	      printf("\tSTOP:                    -   %.2e\n", time_stop_ALL/1000);
  	    #endif
  	  #endif
  	
  	}
  	
  #endif	// MPI
  
  
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
		//                    //
		//    MIXED SOLVER    //
		//                    //
		////////////////////////




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
  
  

  typedef REAL RealT;
  
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
  /*
  size_t dev_spinsize_int   =  6*VOLUME/2*sizeof(dev_spinor);
  int N_sites_int           =    VOLUME/2;
  int N_floats_int          = 24*VOLUME/2;// (single precision) CUBLAS functions get the number of floats as input
  #ifdef MPI
    size_t dev_spinsize_ext =  6*(VOLUME+RAND)/2*sizeof(dev_spinor);
    int N_sites_ext         =    (VOLUME+RAND)/2;
    int N_floats_ext        = 24*(VOLUME+RAND)/2;
  #endif
  */
  
  // algorithm control parameters
  // int N_recalc_res = 10;		// recalculate residue r(k+1) = b - A*x(k+1) each N_recalc_res iteration
  int N_recalc_res = 1000;
  spinor ** up_field = NULL;
  spinor ** dn_field = NULL;
  const int nr_sf = 5;
  
  init_solver_field(&up_field, VOLUMEPLUSRAND/2, nr_sf);
  init_solver_field(&dn_field, VOLUMEPLUSRAND/2, nr_sf);
  
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
  dev_zero_spinor_field<RealT> <<<griddim1, blockdim1>>>(x_up);
  dev_zero_spinor_field<RealT> <<<griddim1, blockdim1>>>(x_dn);
  
  
  // r(0) = b - A*x(0) = b
  dev_copy_spinor_field<RealT,RealT> <<<griddim1, blockdim1>>>(Q_up, r_up);
  dev_copy_spinor_field<RealT,RealT> <<<griddim1, blockdim1>>>(Q_dn, r_dn);
  
  
  // d(0) = r(0)
  dev_copy_spinor_field<RealT,RealT> <<<griddim1, blockdim1>>>(r_up, d_up);
  dev_copy_spinor_field<RealT,RealT> <<<griddim1, blockdim1>>>(r_dn, d_dn);
  
  
  		// debug	// kernel
  		#ifdef CUDA_DEBUG
  		  CUDA_KERNEL_CHECK("Kernel error in cg_eo_nd(). Initializing spinor fields on device failed.", "Spinor fields initialized on device.");
  		#endif
  
  
  
  
  // rr = (r_up)^2 + (r_dn)^2
  #ifndef MPI
    rr_up = cublasDot(N_floats_int, (float *) r_up, 1, (float *) r_up, 1);
    rr_dn = cublasDot(N_floats_int, (float *) r_dn, 1, (float *) r_dn, 1);
  #else
    rr_up = cublasDot_wrapper(N_floats_int, (float *) r_up, 1, (float *) r_up, 1);
    rr_dn = cublasDot_wrapper(N_floats_int, (float *) r_dn, 1, (float *) r_dn, 1);
  #endif
  rr    = rr_up + rr_dn;
  



  r0r0   = rr;	// for relative precision 
  rr_old = rr;	// for the first iteration
  
  
  
  
  //////////
  // LOOP //
  //////////
  
  
  		// debug
  		#ifndef MPI
    		  printf("\nEntering inner loop.\n");
    		#else
    		  if (g_cart_id == 0) printf("\nEntering inner loop.\n");
    		#endif
  
  		// debug	// CUBLAS core function
		#ifdef CUDA_DEBUG
		  // CUBLAS_CORE_CHECK("CUBLAS error in cg_eo_nd(). Calculating initial residue failed.", "Initial inner residue calculated.");
		  CUBLAS_CORE_CHECK_NO_SUCCESS_MSG("CUBLAS error in cg_eo_nd(). Calculating initial residue failed.");
		#endif
  
  		// debug
  		#ifndef MPI
  		  printf("Initial inner residue: %.6e\n", r0r0);
  		#else
  		  if (g_cart_id == 0) printf("Initial inner residue: %.6e\n", r0r0);
  		#endif
  
  
  
  
  for (j = 0; j < max_iter; j++) {
    
    
    #ifndef MATRIX_DEBUG
    
      // A*d(k)
      #ifndef MPI
      		matrix_multiplication32(Ad_up, Ad_dn,										// normally:  matrix_multiplication32()
      		                         d_up,  d_dn,										// debugging: matrix_debug1(), matrix_multiplication_test()
      		                        griddim2, blockdim2,
      		                        griddim3, blockdim3,
      		                        griddim4, blockdim4,
      		                        griddim5, blockdim5);
      #else
      	#ifndef ASYNC
        	matrix_multiplication32_mpi(Ad_up, Ad_dn,									// normally:  matrix_multiplication32_mpi()
        	                             d_up,  d_dn,									// debugging: matrix_mpi_debug1/2/3/4()
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
      #endif	// MPI
      
      
  		// debug	// CUDA		// also other stuff ?!
  		#ifdef CUDA_DEBUG
  		  // CUDA_CHECK("CUDA error in matrix_muliplication32(). Applying the matrix on GPU failed.", "The matrix was applied on GPU.");
  		  CUDA_CHECK_NO_SUCCESS_MSG("CUDA error in matrix_muliplication32(). Applying the matrix on GPU failed.");
  		#endif
    		
    
    #else
    
    		// debug	// apply the host matrix on trial
    		
    		// host/device interaction
    		to_host(up_field[3], d_up, h2d_spin_up, dev_spinsize_int);
    		to_host(dn_field[3], d_dn, h2d_spin_dn, dev_spinsize_int);
    		
    		// matrix multiplication
    		#ifndef MPI
    		  printf("This is Q_Qdagger_ND(). ");
    		#else
    		  if (g_proc_id == 0) printf("This is Q_Qdagger_ND(). ");
    		#endif
    		Q_Qdagger_ND(up_field[4], dn_field[4],			// normally:  Q_Qdagger_ND()
    		             up_field[3], dn_field[3] );		// debugging: matrix_debug2(), Zwitter1(), Zwitter2(), Zwitter3()
    															//       mpi: matrix_mpi_debug10()
    		// host/device interaction
    		to_device(Ad_up, up_field[4], h2d_spin_up, dev_spinsize_int);
    		to_device(Ad_dn, dn_field[4], h2d_spin_dn, dev_spinsize_int);
    		
    		
    				// debug	// CUDA
  				#ifdef CUDA_DEBUG
  			  	// CUDA_CHECK("CUDA error in cg_eo_nd(). Applying the matrix on CPU failed.", "The matrix was applied on CPU.");
  			  	CUDA_CHECK_NO_SUCCESS_MSG("CUDA error in cg_eo_nd(). Applying the matrix on CPU failed.");
  				#endif
    
    #endif	// MATRIX_DEBUG
    
    
    // alpha = r(k)*r(k) / d(k)*A*d(k)
    #ifndef MPI
      dAd_up = cublasDot(N_floats_int, (float *) d_up, 1, (float *) Ad_up, 1);
      dAd_dn = cublasDot(N_floats_int, (float *) d_dn, 1, (float *) Ad_dn, 1);
    #else
      dAd_up = cublasDot_wrapper(N_floats_int, (float *) d_up, 1, (float *) Ad_up, 1);
      dAd_dn = cublasDot_wrapper(N_floats_int, (float *) d_dn, 1, (float *) Ad_dn, 1);
    #endif
    dAd    = dAd_up + dAd_dn;
    
    		// debug	// is NaN ?
    		if isnan(dAd) {
    		  printf("Error in cg_eo_nd(). dAd is NaN.\n");
    		  exit(-1);
    		}
    
    alpha  = rr_old / dAd;	// rr_old is taken from the last iteration respectively
    

    
    
    // x(k+1) = x(k) + alpha*d(k)
    cublasAxpy(N_floats_int, alpha, (RealT*)d_up, 1, (RealT*)x_up, 1);
    cublasAxpy(N_floats_int, alpha, (RealT*)d_dn, 1, (RealT*)x_dn, 1);
    

    
    
    // r(k+1)
    if ( (j+1) % N_recalc_res != 0 ) {	// r(k+1) = r(k) - alpha*A*d(k)
      cublasAxpy(N_floats_int, -1.0*alpha, (RealT*)Ad_up, 1, (RealT*)r_up, 1);
      cublasAxpy(N_floats_int, -1.0*alpha, (RealT*)Ad_dn, 1, (RealT*)r_dn, 1);
    }
    
    else {				// recalculate residue r(k+1) = b - A*x(k+1)
    					//	"feedback"
      		// debug
      		#ifndef MPI
      		  printf("Recalculating the inner residue.\n");
      		#else
      		  if (g_proc_id == 0) printf("Recalculating the inner residue.\n");
      		#endif
      
      
      // A*x(k+1)
      
      #ifndef MATRIX_DEBUG
      
      	#ifndef MPI
        	matrix_multiplication32(Ax_up, Ax_dn,
        	                         x_up,  x_dn,
        	                        griddim2, blockdim2,
        	                        griddim3, blockdim3,
        	                        griddim4, blockdim4,
        	                        griddim5, blockdim5);
        #else
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
        #endif	// MPI
        
      #else
      
      		// debug	// apply the host matrix on trial
    		
    		// host/device interaction
    		to_host(up_field[3], x_up, h2d_spin_up, dev_spinsize_int);
    		to_host(dn_field[3], x_dn, h2d_spin_dn, dev_spinsize_int);
    		
    		// matrix multiplication
    		#ifndef MPI
    		  printf("This is Q_Qdagger_ND(). ");
    		#else
    		  if (g_proc_id == 0) printf("This is Q_Qdagger_ND(). ");
    		#endif
    		Q_Qdagger_ND(up_field[4], dn_field[4],			// normally:       Q_Qdagger_ND()
    		             up_field[3], dn_field[3] );		// debugging, mpi: matrix_mpi_debug10()
    		
    		// host/device interaction
    		to_device(Ax_up, up_field[4], h2d_spin_up, dev_spinsize_int);
    		to_device(Ax_dn, dn_field[4], h2d_spin_dn, dev_spinsize_int);
    		
    		
    				// debug	// CUDA
  				#ifdef CUDA_DEBUG
  			  	// CUDA_CHECK("CUDA error in cg_eo_nd(). Applying the matrix on CPU failed.", "The matrix was applied on CPU.");
  			  	CUDA_CHECK_NO_SUCCESS_MSG("CUDA error in cg_eo_nd(). Applying the matrix on CPU failed.");
  				#endif
      
      #endif	// MATRIX_DEBUG
      
      
      
      
      // r(k+1) = b - A*x(k+1)
      cublasCopy(N_floats_int, (RealT*)Q_up, 1, (RealT*)r_up, 1);		// r_up = Q_up
      cublasCopy(N_floats_int, (RealT*)Q_dn, 1, (RealT*)r_dn, 1);		// r_dn = Q_dn
      cublasAxpy(N_floats_int, -1.0, (RealT*)Ax_up, 1, (RealT*)r_up, 1);	// r_up = Q_up - Ax_up
      cublasAxpy(N_floats_int, -1.0, (RealT*)Ax_dn, 1, (RealT*)r_dn, 1);	// r_dn = Q_dn - Ax_dn
    
    
    } // recalculate residue
    
    
    
    
    // r(k+1)*r(k+1)
    #ifndef MPI
      rr_up  = cublasDot(N_floats_int, (float *) r_up, 1, (float *) r_up, 1);
      rr_dn  = cublasDot(N_floats_int, (float *) r_dn, 1, (float *) r_dn, 1);
    #else
      rr_up  = cublasDot_wrapper(N_floats_int, (float *) r_up, 1, (float *) r_up, 1);
      rr_dn  = cublasDot_wrapper(N_floats_int, (float *) r_dn, 1, (float *) r_dn, 1);
    #endif
    rr     = rr_up + rr_dn;
    
		// debug	// CUBLAS core function
		#ifdef CUDA_DEBUG
		  CUBLAS_CORE_CHECK_NO_SUCCESS_MSG("CUBLAS error in cg_eo_nd(). CUBLAS function failed.");
		#endif
    
    
    		// debug
    		#ifndef MPI
    		  printf("inner iteration j = %i: rr = %.6e\n", j, rr);
    		#else
    		  if (g_proc_id == 0) printf("inner iteration j = %i: rr = %.6e\n", j, rr);
    		#endif
    		
    		// debug	// is NaN ?
    		if isnan(rr) {
    		  printf("Error in cg_eo_nd(). Inner residue is NaN.\n");
    		  exit(-1);
    		}
    
    
    // aborting ?? // check wether precision is reached ...
    if ( (check_abs)&&(rr <= eps_abs) || (check_rel)&&(rr <= eps_rel*r0r0) ) {
    
      #ifdef MPI
        if (g_cart_id == 0) {
      #endif
      
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
      
      #ifdef MPI
        }
      #endif
      
      
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
    cublasScal (N_floats_int, beta, (RealT*)d_up, 1);
    cublasAxpy (N_floats_int, 1.0 , (RealT*)r_up, 1, (RealT*)d_up, 1);
    
    cublasScal (N_floats_int, beta, (RealT*)d_dn, 1);
    cublasAxpy (N_floats_int, 1.0 , (RealT*)r_dn, 1, (RealT*)d_dn, 1);
    
    		// debug	// CUBLAS core function
    		#ifdef CUDA_DEBUG
    		  CUBLAS_CORE_CHECK_NO_SUCCESS_MSG("CUBLAS error in cg_eo_nd(). Error in CUBLAS function.");
    		#endif
    		
  
  }//LOOP
  
  
  		// debug
  		#ifndef MPI
  		  printf("Finished inner loop beacuse of maximal number of inner iterations.\n");
  		  printf("Final inner residue: %.6e\n", rr);
  		#else
  		  if (g_cart_id == 0) printf("Finished inner loop beacuse of maximal number of inner iterations.\n");
  		  if (g_cart_id == 0) printf("Final inner residue: %.6e\n", rr);
  		#endif
  
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
// OUTER SOLVER //
//////////////////

// iterative refinement, defect correction
// that function is to replace the call of  cg_her_nd()  in  invert_doublet_eo.c
// solves the odd part of the full eo and nd problem
//	more precisely we have to invert  Qhat(2x2)*Qhat(2x2)^dagger
//	multiplying by  Qhat(2x2)^dagger  is done in  invert_doublet_eo.c

extern "C" int mixedsolve_eo_nd (spinor * P_up, spinor * P_dn,
                                 spinor * Q_up, spinor * Q_dn,
                                 int max_iter, double eps_sq, int rel_prec) {
  
  typedef REAL RealT;
  
  // basically  P_up/dn  and  Q_up/dn  could be used as auxiliary fields
  //	P_up/dn  is the output field (and can be used as initial guess)
  //	Q_up/dn  is not used later in the calling  invert_doublet_eo.c
  //		but will be used as feedback in r(k+1) = b - A*x(k+1)
  
  
  		// debug
  		#ifdef MPI
  		  if (g_proc_id == 0) {
  		#endif
  		
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
  
  		#ifdef MPI
  		  }
  		#endif
  
  
  
  
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
  #ifdef ALGORITHM_BENCHMARK
    double effectiveflops;			// will used to count the "effective" flop's (from the algorithmic perspective)
    // double hoppingflops = 1488.0;
    double hoppingflops = 1608.0;
    double matrixflops  = 2  *  (  2 * ( (2*hoppingflops+12+3) + (2*hoppingflops+3) + (12+2) + 12 )  );
    #ifdef MPI
      double allflops;				// flops added for all processes
    #endif
  #endif
  
  // timing
  clock_t startouter, stopouter;
  clock_t startinner, stopinner;
  // double timeelapsed;
  clock_t innerclocks;
  clock_t totalinnerclocks = 0;
  clock_t totalouterclocks = 0;
  
  #ifdef ALGORITHM_BENCHMARK
    #ifndef MPI
      clock_t starteffective;
      clock_t stopeffective;
    #else
      double starteffective;
      double stopeffective;
      double singletime;				// time for each process = stopeffective - starteffective
      double maxtime;				// max. parallel process time
    #endif
  #endif
  
  // (auxiliary) fields
  spinor *  r_up, *  r_dn,
         * Ad_up, * Ad_dn,
         *  x_up, *  x_dn,
         *  d_up, *  d_dn,
         * Ax_up, * Ax_dn;
  
  spinor ** up_field = NULL;
  spinor ** dn_field = NULL;
  const int nr_sf = 5;

  init_solver_field(&up_field, VOLUMEPLUSRAND/2, nr_sf);
  init_solver_field(&dn_field, VOLUMEPLUSRAND/2, nr_sf);

  // formal parameters
  /*
  size_t dev_spinsize_int   =  6*VOLUME/2*sizeof(dev_spinor);		// 24 floats per spinor per even lattice site
  int N_sites_int           =    VOLUME/2;				// Carsten's functions get the number of lattice points as input
  int N_floats_int          = 24*VOLUME/2;
  #ifdef MPI
    size_t dev_spinsize_ext =  6*(VOLUME+RAND)/2*sizeof(dev_spinor);
    int N_sites_ext         =    (VOLUME+RAND)/2;
    int N_floats_ext        = 24*(VOLUME+RAND)/2;
  #endif
  */
  
  // algorithm control parameters
  bool rbAx = true;						// choose how to calculate r(k+1)
  bool initial_guess = false;					// choose if initial guess
  
  
  
  
  //////////////////
  // INITIALIZING //
  //////////////////
  
  
  		//debug
  		#ifndef MPI
  		  printf("init_mixedsolve_eo_nd():\n");
  		#else
  		  if (g_cart_id == 0) printf("init_mixedsolve_eo_nd_mpi():\n");
  		#endif
  
  
    init_mixedsolve_eo_nd(g_gauge_field);			// initializes and allocates all quantities for the mixed solver
  								// more precise:
  								//	puts the gauge field on device as "2 rows" or "8 floats" per SU(3)-matrix
  								//	allocates memory for all spinor fields
  								//	puts the nn- and eoidx-fields on device memory

  		//debug
  		#ifndef MPI
  		  printf("mixedsolve_eo_nd():\n");
  		#else
  		  if (g_cart_id == 0) printf("mixedsolve_eo_nd_mpi():\n");
  		#endif
  
  
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
    bind_texture_gf(MixedsolveParameter<RealT>::getGlobalP()->dev_gf);					//	e.g. dev_reconstructgf_2vtexref(...), dev_reconstructgf_8texref(...), 
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
  		
  			#ifdef MPI
  			  if (g_proc_id == 0) {
  			#endif
  			
  			#ifdef MPI
  			  printf("\tOn host:\n");
  			  printf("\tVOLUME = %i\n", VOLUME);							// checking VOLUME and RAND in the parallel case 
  			  printf("\tRAND   = %i\n", RAND);
  			  printf("\tVOLUME + RAND = %i\n",  VOLUME+RAND);
  			#endif
  			
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
  			
  			#ifdef MPI
  			  }
  			#endif
  		
  		#endif
  
  
  he_cg_init_nd_additional<<<1,1>>> (g_mubar, g_epsbar);
  
  		// debug	// kernel
  		#ifdef CUDA_DEBUG
  		  CUDA_KERNEL_CHECK("Kernel error in he_cg_init_nd_additional(). Couldn't initialize some stuff.", "he_cg_init_nd_additional() succeeded.");
  		#endif
  
  		// debug	// check mubar and epsbar on host and device
  		#ifdef STUFF_DEBUG
  		
  			#ifdef MPI
  			  if (g_proc_id == 0) {
  			#endif
  			
  			// printf("\tOn host:\n");
  			// printf("\tg_mubar = %f\n", g_mubar);
  			// printf("\tg_epsbar = %f\n", g_epsbar);
  			
  			float host_check_mubar, host_check_epsbar;
  			cudaMemcpyFromSymbol(&host_check_mubar, mubar, sizeof(float));
  			cudaMemcpyFromSymbol(&host_check_epsbar, epsbar, sizeof(float));
  			printf("\tOn device:\n");
  			printf("\tmubar = %f\n", host_check_mubar);
  			printf("\tepsbar = %f\n", host_check_epsbar);
  			
  			#ifdef MPI
  			  }
  			#endif
  		
  		#endif
  
  
  #ifdef MPI
  
  	he_cg_init_nd_additional_mpi<<<1,1>>>(VOLUMEPLUSRAND, RAND, g_cart_id, g_nproc);
  	
  		// debug	// kernel
  		#ifdef CUDA_DEBUG
  		  CUDA_KERNEL_CHECK("Kernel error in he_cg_init_nd_additional_mpi(). Couldn't initialize some stuff.", "he_cg_init_nd_additional_mpi() succeeded.");
  		#endif
  		
  		// debug
  		#ifdef STUFF_DEBUG
  		
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
  			
  		#endif
  		
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
  
    r_up  = up_field[0];			// use the pre-allocated memory on host memory
    r_dn  = dn_field[0];			// allocated by  init_chi_spinor_field.c  and  invert_doublet.c  !?
    d_up  = up_field[1];		// the fields  g_chi_up/dn_spinor_field[DUM_SOLVER{ , +1, ... , +5}]  are used in  cg_her_nd()
    d_dn  = dn_field[1];
    Ad_up = up_field[2];
    Ad_dn = dn_field[2];
    Ax_up = Ad_up;
    Ax_dn = Ad_dn;
    
  		// debug
  		#ifndef MPI
  		  printf("Now using the fields g_chi_up/dn_spinor_field[DUM_SOLVER{ , +1, +2}] in the mixedsolve_eo_nd().\n");
  		#else
  		  if (g_cart_id == 0) printf("Now using the fields g_chi_up/dn_spinor_field[DUM_SOLVER{ , +1, +2}] in the mixedsolve_eo_nd().\n");
  		#endif
  
  #else
  
  		r_up  = (spinor *) malloc(24*N_sites_int*sizeof(double));		// if using cg_her_nd() as the CG, we cannot use the g_chi_up/dn-fields at the same time
  		r_dn  = (spinor *) malloc(24*N_sites_int*sizeof(double));
  		d_up  = (spinor *) malloc(24*N_sites_int*sizeof(double));
  		d_dn  = (spinor *) malloc(24*N_sites_int*sizeof(double));
  		Ad_up = (spinor *) malloc(24*N_sites_int*sizeof(double));
  		Ad_dn = (spinor *) malloc(24*N_sites_int*sizeof(double));
  		Ax_up = Ad_up;
  		Ax_dn = Ad_dn;
  				// debug
  				#ifndef MPI
  				  printf("Now allocating new host space for the fields in mixedsolve_eo_nd().\n");
  				#else
  				  if (g_cart_id == 0) printf("Now allocating new host space for the fields in mixedsolve_eo_nd().\n");
  				#endif
  
  #endif
  
  
  
  
  ///////////////
  // ALGORITHM //
  ///////////////
  
  // timer
  startouter = clock();
  
  #ifdef ALGORITHM_BENCHMARK
    #ifndef MPI
      starteffective = ((double)clock()) / ((double)(CLOCKS_PER_SEC));
    #else
      starteffective = MPI_Wtime();
    #endif
  #endif
  
  
  // r(0)
  if (!initial_guess) {		// r(0) = b = Q	// for x(0) = 0
    assign(r_up, Q_up, N_sites_int);
    assign(r_dn, Q_dn, N_sites_int);
    #ifndef MPI
      printf("x(0) = 0\n");
    #else
      if (g_cart_id == 0) printf("x(0) = 0\n");
    #endif
  }
  else {			// r(0) = b - A*x(0) = Q - A*P
    bb = square_norm(P_up, N_sites_int, 1) + square_norm(P_dn, N_sites_int, 1);
    #ifndef MPI
      printf("bb = %.10e\n", bb);
    #else
      if (g_cart_id == 0) printf("bb = %.10e\n", bb);
    #endif
    if (bb == 0) {
      assign(r_up, Q_up, N_sites_int);
      assign(r_dn, Q_dn, N_sites_int);
      #ifndef MPI
        printf("x(0) = 0\n");
      #else
        if (g_cart_id == 0) printf("x(0) = 0\n");
      #endif
    }
    else {
      Q_Qdagger_ND(Ax_up, Ax_dn, P_up, P_dn);
      diff(r_up, Q_up, Ax_up, N_sites_int);
      diff(r_dn, Q_dn, Ax_dn, N_sites_int);
      #ifndef MPI
        printf("x(0) != 0\n");
      #else
        if (g_cart_id == 0) printf("x(0) != 0\n");
      #endif
    }
  }
  
  
  // rr = (r_up)^2 + (r_dn)^2
  rr_up = square_norm(r_up, N_sites_int, 1);
  rr_dn = square_norm(r_dn, N_sites_int, 1);
  rr = rr_up + rr_dn;
  
  
  r0r0   = rr; // for relative precision
  rr_old = rr; // for the first iteration
  
  		// debug
  		#ifndef MPI
  		  printf("Initial outer residue: %.10e\n", rr_old);
  		#else
  		  if (g_cart_id == 0) printf("Initial outer residue: %.10e\n", rr_old);
  		#endif
  
  
  // set to zero	// x_up, x_dn  will be added up		// as  x_up/dn = P_up/dn  up to here  P_up/dn  was not changed
  zero_spinor_field(x_up, N_sites_int);
  zero_spinor_field(x_dn, N_sites_int);
  
  
  
  
  ////////////////
  // OUTER LOOP //
  ////////////////
  
  		// debug
  		#ifndef MPI
    		  printf("\nEntering outer loop.");
    		#else
    		  if (g_cart_id == 0) printf("\nEntering outer loop.");
    		#endif
  
  
  do {		// for (i = 0; i < max_iter; i++) {
    
    i++;
  
    		// debug
    		#ifndef MPI
    		  printf("\nouter iteration i = %i\n", i);
    		#else
    		  if (g_cart_id == 0) printf("\nouter iteration i = %i\n", i);
    		#endif
    
    
    
    
    #ifndef CG_DEBUG
    
    // host/device interaction
    to_device(dev_spinin_up, r_up, h2d_spin_up, dev_spinsize_int);		// notice: for MPI communicateion the boundary exchange takes place when the hopping matrix is applied
    to_device(dev_spinin_dn, r_dn, h2d_spin_dn, dev_spinsize_int);
    
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
    		#ifndef MPI
    		  printf("cg_eo_nd():\n");
    		#else
    		  if (g_cart_id == 0) printf("cg_eo_nd():\n");
    		#endif
    
    
    // solves A*p(k+1) = r(k)
    //        A*p(0)   = r(0) = b
    innercount = cg_eo_nd(MixedsolveParameter<RealT>::getGlobalP()->dev_gf,
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
    		#ifndef MPI
    		  printf("Inner solver done in: %.4e sec\n", double(innerclocks) / double(CLOCKS_PER_SEC));
    		#else
    		  if (g_cart_id == 0) printf("Inner solver done in: %.4e sec\n", double(innerclocks) / double(CLOCKS_PER_SEC));
    		#endif
    
    
    // host/device interaction
    to_host(d_up, dev_spinout_up, h2d_spin_up, dev_spinsize_int);
    to_host(d_dn, dev_spinout_dn, h2d_spin_dn, dev_spinsize_int);
    
    		// debug	// CUDA
    		#ifdef CUDA_DEBUG
    		  // CUDA_CHECK("CUDA error in mixedsolve_eo_nd(). Device to host interaction failed.", "Fields copied back to device.");
    		  CUDA_CHECK_NO_SUCCESS_MSG("CUDA error in mixedsolve_eo_nd(). Device to host interaction failed.");
    		#endif
    
    
    #else
    
    
    				// debug
    				#ifndef MPI
    				  printf("cg_her_nd():\n");
    				#else
    				  if (g_cart_id == 0) printf("cg_her_nd():\n");
    				#endif
    		
    		innercount = cg_her_nd(d_up, d_dn, r_up, r_dn,		// MISTAKE, was: r_up, r_dn, d_up, d_dn,
				       1000, eps_sq/2, 0,
				       VOLUME/2, &Q_Qdagger_ND, 0, 1000);
    		
    		outercount = outercount + innercount;
    		
    				// debug
    				#ifndef MPI
    				  printf("cg_her_nd() on host was used for debugging purposes.\n");
    				#else
    				  if (g_cart_id == 0) printf("cg_her_nd() on host was used for debugging purposes.\n");
    				#endif
    
    
    #endif
    
    
    		// debug
    		#ifndef MPI
    		  printf("mixedsolve_eo_nd():\n");
    		#else
    		  if (g_cart_id == 0) printf("mixedsolve_eo_nd():\n");
    		#endif
    
    
    // x(k+1) = x(k) + d(k+1)
    add(x_up, x_up, d_up, N_sites_int);
    add(x_dn, x_dn, d_dn, N_sites_int);
    

    
    
    // r(k+1)
    if (rbAx) {				// r(k+1) = b - A*x(k+1)
      // A*x(k+1)
      Q_Qdagger_ND(Ax_up, Ax_dn, x_up, x_dn);
      		// debug
      		#ifndef MPI
      		  printf("The matrix was applied on CPU in double precision. r = b - Ax\n");
      		#else
      		  if (g_cart_id == 0) printf("The matrix was applied on CPU in double precision. r = b - Ax\n");
      		#endif
      diff(r_up, Q_up, Ax_up, N_sites_int);
      diff(r_dn, Q_dn, Ax_dn, N_sites_int);
    }
    else {				// r(k+1) = r(k) - A*d(k+1)	// makes actually no sense ;)
      // A*d(k+1)
      Q_Qdagger_ND(Ad_up, Ad_dn, d_up, d_dn);
    		// debug
    		#ifndef MPI
    		  printf("The matrix was applied on CPU in double precision. r = r - Ad\n");
    		#else
    		  if (g_cart_id == 0) printf("The matrix was applied on CPU in double precision. r = r - Ad\n");
    		#endif
      // r(k+1) = r(k) - A*d(k+1)
      diff(r_up, r_up, Ad_up, N_sites_int);
      diff(r_dn, r_dn, Ad_dn, N_sites_int);
    }
    
    
    
    
    // rr = (rr_up)^2 + (r_dn)^2
    rr_up = square_norm(r_up, N_sites_int, 1);
    rr_dn = square_norm(r_dn, N_sites_int, 1);
    rr    = rr_up + rr_dn;
    
    		// debug
    		#ifndef MPI
    		  printf("Outer residue in the outer iteration i = %i after %i total inner iterations : %.10e\n", i, outercount, rr);
    		#else
    		  if (g_cart_id == 0) printf("Outer residue in the outer iteration i = %i after %i total inner iterations : %.10e\n", i, outercount, rr);
    		#endif
    		
    		// debug	// is NaN ?
    		if isnan(rr) {
    		  printf("Error in mixedsolve_eo_nd(). Outer residue is NaN.\n");
    		  exit(-1);
    		}
    		

    
    
    // aborting ?? // check wether precision is reached ...
    if ( ((rr <= eps_sq) && (rel_prec == 0))  ||  ((rr <= eps_sq*r0r0) && (rel_prec == 1)) ) {
      
      // timer
      stopouter = clock();
      totalouterclocks = stopouter-startouter - totalinnerclocks;
      
      #ifdef ALGORITHM_BENCHMARK
        #ifndef MPI
          stopeffective = ((double)clock()) / ((double)(CLOCKS_PER_SEC));
        #else
          stopeffective = MPI_Wtime();
        #endif
      #endif
      
      
      		// debug
      		#ifdef MPI
      		  if (g_cart_id == 0) {
      		#endif
      		printf("\nEO inversion done in mixed precision.\n");
      		if (rel_prec == 0) printf("Finished outer loop because of reached absolute outer solver precision.\n");
      		if (rel_prec == 1) printf("Finished outer loop because of reached relative outer solver precision.\n");
      		printf("Total number of inner iterations: %i\n", outercount);
      		printf("Total number of outer iterations: %i\n", i+1);
      		printf("Squared residue: %.10e\n", rr); 
      		printf("Outer solver done in: %.4e sec\n", double(stopouter-startouter) / double(CLOCKS_PER_SEC));
      		#ifdef MPI
      		  }
      		#endif
      		
      		// benchmark
      		#ifdef ALGORITHM_BENCHMARK
      		  // will now count the number of effective flops
      		  // effectiveflops  =  #(inner iterations)*(matrixflops+linalgflops)*VOLUME/2  +  #(outer iterations)*(matrixflops+linalgflops)*VOLUME/2
      		  // outer loop: linalg  =  flops for calculating  r(k+1) and x(k+1)
      		  // inner loop: linalg  =  flops for calculating  alpha, x(k+1), r(k+1), beta, d(k+1)
      		  #ifndef MPI
      		  	effectiveflops = outercount*(matrixflops + 2*2*2*24 + 2*2*24 + 2*2*24 + 2*2*2*24 + 2*2*24)*VOLUME/2   +   i*(matrixflops + 2*24 + 2*24)*VOLUME/2;
      		  	printf("effective BENCHMARK:\n");
      		  	printf("\ttotal mixed solver time:   %.4e sec\n", double(stopeffective-starteffective));
      		  	printf("\tfloating point operations: %.4e flops\n", effectiveflops);
      		  	printf("\tinner solver performance:  %.4e Gflop/s\n", double(effectiveflops) / double(stopeffective-starteffective) / 1.0e9);
      		  #else
      		  	singletime = double(stopeffective-starteffective);
      		  	effectiveflops = outercount*(matrixflops + 2*2*2*24 + 2*2*24 + 2*2*24 + 2*2*2*24 + 2*2*24)*VOLUME/2   +   i*(matrixflops + 2*24 + 2*24)*VOLUME/2;
      		  	MPI_Allreduce(&singletime, &maxtime, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      		  	MPI_Allreduce(&effectiveflops, &allflops, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      		  	if (g_proc_id == 0) printf("effective BENCHMARK:\n");
      		  	if (g_proc_id == 0) printf("\ttotal mixed solver time:   %.4e sec\n", double(maxtime));
      		  	if (g_proc_id == 0) printf("\tfloating point operations: %.4e flops\n", double(allflops));
      		  	if (g_proc_id == 0) printf("\tinner solver performance:  %.4e Gflop/s\n", double(allflops) / double(maxtime) / 1.0e9);
      		  	/*
      		  	printf("this is for checking:\n");
      		  	printf("\ttotal mixed solver time:   %.2e sec\n", double(stopeffective-starteffective));
      		  	printf("\tfloating point operations: %.2e flops\n", effectiveflops);
      		  	printf("\tinner solver performance:  %.2e Gflop/s\n", double(effectiveflops) / double(stopeffective-starteffective) / 1.0e9);
      		  	*/
      		  #endif
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
      		#ifndef MPI
      		  printf("finalize_mixedsolve_eo_nd():\n");
      		#else
      		  if (g_cart_id == 0) printf("finalize_mixedsolve_eo_nd():\n");
      		#endif
      
      finalize_mixedsolve_eo_nd();
      
      		// debug
      		#ifndef MPI
      		  printf("\n");
      		#else
      		  if (g_cart_id == 0) printf("\n");
      		#endif
      finalize_solver(up_field, nr_sf);
      finalize_solver(dn_field, nr_sf); 
      return(outercount);
      
    }
    
    
    
    
  }//OUTER LOOP
  while (outercount <= max_iter);
  
  
  // multiplying with Qhat(2x2)^dagger is done in invert_doublet_eo.c
  
  
  // timer
  stopouter = clock();
  totalouterclocks = stopouter-startouter - totalinnerclocks;
  
  #ifdef ALGORITHM_BENCHMARK
    #ifndef MPI
      stopeffective = ((double)clock()) / ((double)(CLOCKS_PER_SEC));
    #else
      stopeffective = MPI_Wtime();
    #endif
  #endif
  
  
  		// debug
  		#ifdef MPI
  		  if (g_cart_id == 0) {
  		#endif
  		printf("\nEO inversion done in mixed precision.\n");
  		printf("Finished outer loop, because of maximal number of outer iterations.\n");
      		printf("Total number of inner iterations: %i\n", outercount);
      		printf("Total number of outer iterations: %i\n", i+1);
      		printf("Squared residue: %.10e\n", rr); 
      		printf("Outer solver done in: %.4e sec\n", double(stopouter-startouter)/CLOCKS_PER_SEC);
      		#ifdef MPI
      		  }
      		#endif
      		
      		// benchmark
      		#ifdef ALGORITHM_BENCHMARK
      		  // will now count the number of effective flops
      		  // effectiveflops  =  #(inner iterations)*(matrixflops+linalgflops)*VOLUME/2  +  #(outer iterations)*(matrixflops+linalgflops)*VOLUME/2
      		  // outer loop: linalg  =  flops for calculating  r(k+1) and x(k+1)
      		  // inner loop: linalg  =  flops for calculating  alpha, x(k+1), r(k+1), beta, d(k+1)
      		  #ifndef MPI
      		  	effectiveflops = outercount*(matrixflops + 2*2*2*24 + 2*2*24 + 2*2*24 + 2*2*2*24 + 2*2*24)*VOLUME/2   +   i*(matrixflops + 2*24 + 2*24)*VOLUME/2;
      		  	printf("effective BENCHMARK:\n");
      		  	printf("\ttotal mixed solver time:   %.4e sec\n", double(stopeffective-starteffective));
      		  	printf("\tfloating point operations: %.4e flops\n", effectiveflops);
      		  	printf("\tinner solver performance:  %.4e Gflop/s\n", double(effectiveflops) / double(stopeffective-starteffective) / 1.0e9);
      		  #else
      		  	singletime = double(stopeffective-starteffective);
      		  	effectiveflops = outercount*(matrixflops + 2*2*2*24 + 2*2*24 + 2*2*24 + 2*2*2*24 + 2*2*24)*VOLUME/2   +   i*(matrixflops + 2*24 + 2*24)*VOLUME/2;
      		  	MPI_Allreduce(&singletime, &maxtime, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      		  	MPI_Allreduce(&effectiveflops, &allflops, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      		  	if (g_proc_id == 0) printf("effective BENCHMARK:\n");
      		  	if (g_proc_id == 0) printf("\ttotal mixed solver time:   %.4e sec\n", double(maxtime));
      		  	if (g_proc_id == 0) printf("\tfloating point operations: %.4e flops\n", double(allflops));
      		  	if (g_proc_id == 0) printf("\tinner solver performance:  %.4e Gflop/s\n", double(allflops) / double(maxtime) / 1.0e9);
      		  	/*
      		  	printf("this is for checking:\n");
      		  	printf("\ttotal mixed solver time:   %.2e sec\n", double(stopeffective-starteffective));
      		  	printf("\tfloating point operations: %.2e flops\n", effectiveflops);
      		  	printf("\tinner solver performance:  %.2e Gflop/s\n", double(effectiveflops) / double(stopeffective-starteffective) / 1.0e9);
      		  	*/
      		  #endif
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
  		#ifndef MPI
  		  printf("finalize_mixedsolve_eo_nd():\n");  
  		#else
  		  if (g_cart_id == 0) printf("finalize_mixedsolve_eo_nd():\n");
  		#endif
  
  finalize_mixedsolve_eo_nd();
  
  		// debug
  		#ifndef MPI
  		  printf("\n");
  		#else
  		  if (g_cart_id == 0) printf("\n");
  		#endif
  
  finalize_solver(up_field, nr_sf);
  finalize_solver(dn_field, nr_sf);

  return(outercount);
  
  
}//mixedsolve_eo_nd()



















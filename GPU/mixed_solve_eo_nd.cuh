/**************************************************************************
 *
 * Copyright (C) 2010 Joseph Nagel
 *               2012 Florian Burger
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
#include "../operator/tm_operators_nd.h"
#include "../operator/Hopping_Matrix.h"
#include "../solver/cg_her_nd.h"
}
#include "../global.h"
#include "HEADER.h"
#ifdef MPI
    #include <mpi.h>
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
dev_spinor * dev_spin1_up;		// auxiliary fields for dev_cg_eo_nd()
dev_spinor * dev_spin1_dn;
dev_spinor * dev_spin2_up;
dev_spinor * dev_spin2_dn;
dev_spinor * dev_spin3_up;
dev_spinor * dev_spin3_dn;
dev_spinor * dev_spinin_up;		// host/device interaction	// mixedsolve_eo_nd()  <-->  dev_cg_eo_nd()
dev_spinor * dev_spinin_dn;		// inner/outer interaction
dev_spinor * dev_spinout_up;
dev_spinor * dev_spinout_dn;
dev_spinor * h2d_spin_up;		// for transferring in double precision on host to single precision on device (pointing to host)
dev_spinor * h2d_spin_dn;
dev_spinor * dev_spin_eo1_up;		// auxiliary for  dev_Qtm_pm_ndpsi()  called by  dev_cg_eo_nd()
dev_spinor * dev_spin_eo1_dn;
dev_spinor * dev_spin_eo2_up;
dev_spinor * dev_spin_eo2_dn;
dev_spinor * dev_spin_eo3_up;
dev_spinor * dev_spin_eo3_dn;






#ifdef MPI					// collecting variables for the MPI implementation
  						// put to mixed_solve.cu
 
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
__global__ void he_cg_init_nd_additional (double param_mubar, double param_epsbar) {

  mubar  = (float) param_mubar;
  epsbar = (float) param_epsbar;
  mubar_d  = param_mubar;
  epsbar_d = param_epsbar; 
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

// initializes and allocates all quantities for the ND mixed solver
// all degenerate related initializations done in init_mixedsolve_eo
// more precise:
//	puts the gauge field on device as "2 rows" or "8 floats" per SU(3)-matrix
//	allocates memory for all ND spinor fields


void init_mixedsolve_eo_nd (su3** gf) {	// gf is the full gauge field
  
  
  //////////////////////
  // GLOBAL VARIABLES //
  //////////////////////
  
  dev_spinsize_int   =  6*VOLUME/2*sizeof(dev_spinor);				// 24 floats per lattice site
  N_sites_int        =    VOLUME/2;
  N_floats_int       = 24*VOLUME/2;
  #ifdef MPI
    dev_spinsize_ext =  6*(VOLUME+RAND)/2*sizeof(dev_spinor);
    N_sites_ext      =    (VOLUME+RAND)/2;
    N_floats_ext     = 24*(VOLUME+RAND)/2;
  #endif

  
  set_global_sizes();
  
  
  
  /////////////////////
  // LOCAL VARIABLES //
  /////////////////////
  
  cudaError_t cudaerr;		// CUDA errors
  int grid[6];			// array for grid specifications
  
  

  /////////////
  // SPINORS //							// allocates device memory for the odd part of the spinor fields (24 = 4*6 floats per odd lattice sites)
  /////////////							// now we have to consider 2 flavors: up, dn
  

  #ifndef MPI
  
    cudaMalloc((void **) &dev_spin1_up, dev_spinsize_int);   	// allocates device memory for the fields spinor fields used in dev_cg_eo_nd(...)
    cudaMalloc((void **) &dev_spin1_dn, dev_spinsize_int);	// pointing to device
    cudaMalloc((void **) &dev_spin2_up, dev_spinsize_int);	// ...
    cudaMalloc((void **) &dev_spin2_dn, dev_spinsize_int);
    cudaMalloc((void **) &dev_spin3_up, dev_spinsize_int);
    cudaMalloc((void **) &dev_spin3_dn, dev_spinsize_int);
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
    cudaMalloc((void **) &dev_spinin_up , dev_spinsize_ext);
    cudaMalloc((void **) &dev_spinin_dn , dev_spinsize_ext);
    cudaMalloc((void **) &dev_spinout_up, dev_spinsize_ext);
    cudaMalloc((void **) &dev_spinout_dn, dev_spinsize_ext);
  
  #endif
  if((cudaerr=cudaGetLastError())!=cudaSuccess){
    if(g_proc_id==0) printf("Error in init_mixedsolve_eo(): Memory allocation of nd additional double spinor fields failed. Aborting...\n");
    exit(200);
  }
  
  

  #ifdef GPU_DOUBLE
   	  size_t dev_spinsize_d = 12*VOLUME/2 * sizeof(dev_spinor_d); /* double4 */
	  //allocate fields used in dev_Qtm_pm_ndpsi_d
  	  cudaMalloc((void **) &dev_spin_eo1_up_d, dev_spinsize_d);	
  	  cudaMalloc((void **) &dev_spin_eo1_dn_d, dev_spinsize_d);	  
  	  cudaMalloc((void **) &dev_spin_eo3_up_d, dev_spinsize_d);	
  	  cudaMalloc((void **) &dev_spin_eo3_dn_d, dev_spinsize_d);
	  
  	  if((cudaerr=cudaGetLastError())!=cudaSuccess){
              if(g_proc_id==0) printf("Error in init_mixedsolve_eo(): Memory allocation of nd additional double spinor fields failed. Aborting...\n");
              exit(200);
          }
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
  
    cudaMalloc((void **) &dev_spin_eo1_up, dev_spinsize_int);		// used for dev_Qtm_pm_ndpsi(...)
    cudaMalloc((void **) &dev_spin_eo1_dn, dev_spinsize_int);
    cudaMalloc((void **) &dev_spin_eo3_up, dev_spinsize_int);
    cudaMalloc((void **) &dev_spin_eo3_dn, dev_spinsize_int);

  #else
  
    cudaMalloc((void **) &dev_spin_eo1_up, dev_spinsize_ext);
    cudaMalloc((void **) &dev_spin_eo1_dn, dev_spinsize_ext);
    cudaMalloc((void **) &dev_spin_eo3_up, dev_spinsize_ext);
    cudaMalloc((void **) &dev_spin_eo3_dn, dev_spinsize_ext);
 
  #endif
  
  		// debug	// CUDA
  		#ifdef CUDA_DEBUG
  		  #ifndef MPI
  		    CUDA_CHECK("CUDA error in init_mixedsolve_eo_nd(). Memory allocation of spinor fields failed.", "Allocated spinor fields on device.");
  		  #else
  		    CUDA_CHECK("CUDA error in init_mixedsolve_eo_nd(). Memory allocation of spinor fields failed.", "Allocated spinor fields on devices.");
  		  #endif
  		#endif

  for (int i = 0; i < 2; i++) {
        cudaStreamCreate(&stream_nd[i]);
  }   

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
  
  
  ////////////////////////////
  // grid[ ] specifications //							// allocate and initializes the array grid[5] on device
  ////////////////////////////
  
  grid[0] = LX;									// it contains the dimensions of the lattice and the volume of the eo-sublattice
  grid[1] = LY;
  grid[2] = LZ;
  grid[3] = T;
  grid[4] = VOLUME/2;								// will be used to set dev_VOLUME: dev_VOLUME is half of VOLUME for eo
  
  // put dev_Offset accordingly depending on mpi/non-mpi
  #ifdef MPI
   grid[5] = (VOLUME+RAND)/2;
  #else
   grid[5] = VOLUME/2;
  #endif

  //done in init_mixedsolve_eo
  //cudaMalloc((void **) &dev_grid, 6*sizeof(int));				// dev_grid
  cudaMemcpy(dev_grid, &(grid[0]), 6*sizeof(int), cudaMemcpyHostToDevice);
  
  
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
  
  cudaError_t cudaerr;
  
  cudaFree(dev_spin1_up);
  cudaFree(dev_spin1_dn);
  cudaFree(dev_spin2_up);
  cudaFree(dev_spin2_dn);
  cudaFree(dev_spin3_up);
  cudaFree(dev_spin3_dn);
  cudaFree(dev_spinin_up);
  cudaFree(dev_spinin_dn);
  cudaFree(dev_spinout_up);
  cudaFree(dev_spinout_dn);
  #ifdef GPU_DOUBLE  
    cudaFree(dev_spin_eo1_up_d);
    cudaFree(dev_spin_eo1_dn_d);     
    cudaFree(dev_spin_eo3_up_d);
    cudaFree(dev_spin_eo3_dn_d);
  #endif
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
  
  for (int i = 0; i < 2; i++) {
    cudaStreamDestroy(stream_nd[i]);
  }
  
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
  //cudaThreadExit();				// is essential
  
  
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
  
        h2d[i].s0.c0 = spin[6*i+0].x + I* spin[6*i+0].y;
        h2d[i].s0.c1 = spin[6*i+0].z + I* spin[6*i+0].w;
        
        h2d[i].s0.c2 = spin[6*i+1].x + I* spin[6*i+1].y;
        h2d[i].s1.c0 = spin[6*i+1].z + I* spin[6*i+1].w;   
        
        h2d[i].s1.c1 = spin[6*i+2].x + I* spin[6*i+2].y;
        h2d[i].s1.c2 = spin[6*i+2].z + I* spin[6*i+2].w;  
        
        h2d[i].s2.c0 = spin[6*i+3].x + I* spin[6*i+3].y;
        h2d[i].s2.c1 = spin[6*i+3].z + I* spin[6*i+3].w;  
        
        h2d[i].s2.c2 = spin[6*i+4].x + I* spin[6*i+4].y;
        h2d[i].s3.c0 = spin[6*i+4].z + I* spin[6*i+4].w; 
        
        h2d[i].s3.c1 = spin[6*i+5].x + I* spin[6*i+5].y;
        h2d[i].s3.c2 = spin[6*i+5].z + I* spin[6*i+5].w; 
        
  }
}



// convert spinor to REAL4 (float4, double4)

void convert2REAL4_spin_mpi (spinor * spin, dev_spinor * h2d, int start, int end) {

  int i;
  
  for (i = start; i < end; i++) {
    
        h2d[6*i+0].x = creal(spin[i].s0.c0);
        h2d[6*i+0].y = cimag(spin[i].s0.c0);
        h2d[6*i+0].z = creal(spin[i].s0.c1);
        h2d[6*i+0].w = cimag(spin[i].s0.c1);
        
        h2d[6*i+1].x = creal(spin[i].s0.c2);
        h2d[6*i+1].y = cimag(spin[i].s0.c2);
        h2d[6*i+1].z = creal(spin[i].s1.c0);
        h2d[6*i+1].w = cimag(spin[i].s1.c0);
        
        h2d[6*i+2].x = creal(spin[i].s1.c1);
        h2d[6*i+2].y = cimag(spin[i].s1.c1);
        h2d[6*i+2].z = creal(spin[i].s1.c2);
        h2d[6*i+2].w = cimag(spin[i].s1.c2);
        
        h2d[6*i+3].x = creal(spin[i].s2.c0);
        h2d[6*i+3].y = cimag(spin[i].s2.c0);
        h2d[6*i+3].z = creal(spin[i].s2.c1);
        h2d[6*i+3].w = cimag(spin[i].s2.c1);
        
        h2d[6*i+4].x = creal(spin[i].s2.c2);
        h2d[6*i+4].y = cimag(spin[i].s2.c2);
        h2d[6*i+4].z = creal(spin[i].s3.c0);
        h2d[6*i+4].w = cimag(spin[i].s3.c0);
        
        h2d[6*i+5].x = creal(spin[i].s3.c1);
        h2d[6*i+5].y = cimag(spin[i].s3.c1);
        h2d[6*i+5].z = creal(spin[i].s3.c2);
        h2d[6*i+5].w = cimag(spin[i].s3.c2);
    
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

  convert2REAL4_spin(host, auxiliary);						// auxiliary = (float) host
  cudaMemcpy(device, auxiliary, size, cudaMemcpyHostToDevice);			// device = auxiliary  (on device)

}


void to_host (spinor * host, dev_spinor * device, dev_spinor * auxiliary, int size) {

  cudaMemcpy(auxiliary, device, size, cudaMemcpyDeviceToHost);			// auxiliary = device  (on device)
  convert2double_spin(auxiliary, host);						// host = (double) auxiliary

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



// a wrapper function for cublasSdot() (with the same interface)
// provides the MPI communication via MPI_Allreduce()

float cublasSdot_wrapper(int size, float * A, int incx, float * B, int incy) {

  float result;
  float buffer;
  
  buffer = cublasSdot(size, (float *) A, incx, (float *) B, incy);
  //printf("proc no %d: my square is: %.8f\n", g_proc_id, buffer);
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

// derived from function  dev_mul_one_pm_imu_inv
//	order of the arguments (spinin, spinout)
// applies (1 +- imubar*gamma5)
// one thread per lattice site
__global__ void dev_mul_one_pm_imubar_gamma5 (dev_spinor * sin,
                                              dev_spinor * sout,
                                              float sign         ) {
   
  dev_spinor slocal[6];									// dev_spinor = float4		// 6*float4 = 24 floats		// auxiliary for each thread
  dev_complex pm_imu = dev_initcomplex(0.0, sign * mubar);				// dev_complex = struct { REAL re; REAL im; }	// pm_imu.re = 0.0																	// pm_imu.im = sign * mubar
  int pos = threadIdx.x + blockDim.x*blockIdx.x;
  
  if (pos < dev_VOLUME) {
    #ifdef RELATIVISTIC_BASIS
      dev_skalarmult_gamma5_globalspinor_rel(&(slocal[0]), pm_imu, &(sin[pos]) );			// slocal  =  pm_imu * (gamma5) * sin
    #else
      dev_skalarmult_gamma5_globalspinor(&(slocal[0]), pm_imu, &(sin[pos]) );			// slocal  =  pm_imu * (gamma5) * sin
    #endif
    dev_add_globalspinor_assign(&(slocal[0]), &(sin[pos]));					// slocal  =  slocal + sin  =  pm_imu * (gamma5) * sin + sin
    dev_realmult_spinor_assigntoglobal(&(sout[pos]), 1.0, &(slocal[0]) );			// sout    =  slocal
  }
}



//s2_up = s2_up -epsbar s3_dn - s1_up
//s2_dn = s2_dn -epsbar s3_up - s1_dn
__global__ void dev_nd_linalg1 (dev_spinor * s1_up, dev_spinor * s1_dn,
				dev_spinor * s3_up, dev_spinor * s3_dn,  
				dev_spinor * s2_up, dev_spinor * s2_dn,	float epsbar			
				) {
  dev_spinor s1[6], s3[6], sout[6];	
  int pos = threadIdx.x + blockDim.x*blockIdx.x;
  
  if (pos < dev_VOLUME) {
    //upper output spinor
    dev_read_spinor(&(sout[0]), &(s2_up[pos]));
    dev_read_spinor(&(s1[0]), &(s1_up[pos]));    
    dev_read_spinor(&(s3[0]), &(s3_up[pos])); 
    
    sout[0].x -= epsbar*s3[0].x;
    sout[0].x -= s1[0].x;
    sout[0].y -= epsbar*s3[0].y;
    sout[0].y -= s1[0].y;    
    sout[0].z -= epsbar*s3[0].z;
    sout[0].z -= s1[0].z;    
    sout[0].w -= epsbar*s3[0].w;
    sout[0].w -= s1[0].w;       
    
    sout[1].x -= epsbar*s3[1].x;
    sout[1].x -= s1[1].x;
    sout[1].y -= epsbar*s3[1].y;
    sout[1].y -= s1[1].y;    
    sout[1].z -= epsbar*s3[1].z;
    sout[1].z -= s1[1].z;    
    sout[1].w -= epsbar*s3[1].w;
    sout[1].w -= s1[1].w;       
    
    sout[2].x -= epsbar*s3[2].x;
    sout[2].x -= s1[2].x;
    sout[2].y -= epsbar*s3[2].y;
    sout[2].y -= s1[2].y;    
    sout[2].z -= epsbar*s3[2].z;
    sout[2].z -= s1[2].z;    
    sout[2].w -= epsbar*s3[2].w;
    sout[2].w -= s1[2].w;       
    
    sout[3].x -= epsbar*s3[3].x;
    sout[3].x -= s1[3].x;
    sout[3].y -= epsbar*s3[3].y;
    sout[3].y -= s1[3].y;    
    sout[3].z -= epsbar*s3[3].z;
    sout[3].z -= s1[3].z;    
    sout[3].w -= epsbar*s3[3].w;
    sout[3].w -= s1[3].w;   
    
    sout[4].x -= epsbar*s3[4].x;
    sout[4].x -= s1[4].x;
    sout[4].y -= epsbar*s3[4].y;
    sout[4].y -= s1[4].y;    
    sout[4].z -= epsbar*s3[4].z;
    sout[4].z -= s1[4].z;    
    sout[4].w -= epsbar*s3[4].w;
    sout[4].w -= s1[4].w;       
    
    sout[5].x -= epsbar*s3[5].x;
    sout[5].x -= s1[5].x;
    sout[5].y -= epsbar*s3[5].y;
    sout[5].y -= s1[5].y;    
    sout[5].z -= epsbar*s3[5].z;
    sout[5].z -= s1[5].z;    
    sout[5].w -= epsbar*s3[5].w;
    sout[5].w -= s1[5].w;       
    
    dev_write_spinor(&(sout[0]),&(s2_up[pos]));    
    
    
    //lower output spinor
    dev_read_spinor(&(sout[0]), &(s2_dn[pos]));
    dev_read_spinor(&(s1[0]), &(s1_dn[pos]));    
    dev_read_spinor(&(s3[0]), &(s3_dn[pos]));    
    

    sout[0].x -= epsbar*s3[0].x;
    sout[0].x -= s1[0].x;
    sout[0].y -= epsbar*s3[0].y;
    sout[0].y -= s1[0].y;    
    sout[0].z -= epsbar*s3[0].z;
    sout[0].z -= s1[0].z;    
    sout[0].w -= epsbar*s3[0].w;
    sout[0].w -= s1[0].w;       
    
    sout[1].x -= epsbar*s3[1].x;
    sout[1].x -= s1[1].x;
    sout[1].y -= epsbar*s3[1].y;
    sout[1].y -= s1[1].y;    
    sout[1].z -= epsbar*s3[1].z;
    sout[1].z -= s1[1].z;    
    sout[1].w -= epsbar*s3[1].w;
    sout[1].w -= s1[1].w;       
    
    sout[2].x -= epsbar*s3[2].x;
    sout[2].x -= s1[2].x;
    sout[2].y -= epsbar*s3[2].y;
    sout[2].y -= s1[2].y;    
    sout[2].z -= epsbar*s3[2].z;
    sout[2].z -= s1[2].z;    
    sout[2].w -= epsbar*s3[2].w;
    sout[2].w -= s1[2].w;       
    
    sout[3].x -= epsbar*s3[3].x;
    sout[3].x -= s1[3].x;
    sout[3].y -= epsbar*s3[3].y;
    sout[3].y -= s1[3].y;    
    sout[3].z -= epsbar*s3[3].z;
    sout[3].z -= s1[3].z;    
    sout[3].w -= epsbar*s3[3].w;
    sout[3].w -= s1[3].w;   
    
    sout[4].x -= epsbar*s3[4].x;
    sout[4].x -= s1[4].x;
    sout[4].y -= epsbar*s3[4].y;
    sout[4].y -= s1[4].y;    
    sout[4].z -= epsbar*s3[4].z;
    sout[4].z -= s1[4].z;    
    sout[4].w -= epsbar*s3[4].w;
    sout[4].w -= s1[4].w;       
    
    sout[5].x -= epsbar*s3[5].x;
    sout[5].x -= s1[5].x;
    sout[5].y -= epsbar*s3[5].y;
    sout[5].y -= s1[5].y;    
    sout[5].z -= epsbar*s3[5].z;
    sout[5].z -= s1[5].z;    
    sout[5].w -= epsbar*s3[5].w;
    sout[5].w -= s1[5].w;         
    
    
    dev_write_spinor(&(sout[0]),&(s2_dn[pos]));       
  }
}








//s2_up = gamma5*(s2_up -epsbar s3_up - s1_up)
//s2_dn = gamma5*(s2_dn -epsbar s3_dn - s1_dn)
__global__ void dev_nd_linalg1_gamma5 (dev_spinor * s1_up, dev_spinor * s1_dn,
				dev_spinor * s3_up, dev_spinor * s3_dn,  
				dev_spinor * s2_up, dev_spinor * s2_dn,	float epsbar			
				) {
  dev_spinor s1[6], s3[6], sout[6];	
  int pos = threadIdx.x + blockDim.x*blockIdx.x;
  
  if (pos < dev_VOLUME) {
    //upper output spinor
    dev_read_spinor(&(sout[0]), &(s2_up[pos]));
    dev_read_spinor(&(s1[0]), &(s1_up[pos]));    
    dev_read_spinor(&(s3[0]), &(s3_up[pos])); 
    
    sout[0].x -= epsbar*s3[0].x;
    sout[0].x -= s1[0].x;
    sout[0].y -= epsbar*s3[0].y;
    sout[0].y -= s1[0].y;    
    sout[0].z -= epsbar*s3[0].z;
    sout[0].z -= s1[0].z;    
    sout[0].w -= epsbar*s3[0].w;
    sout[0].w -= s1[0].w;       
    
    sout[1].x -= epsbar*s3[1].x;
    sout[1].x -= s1[1].x;
    sout[1].y -= epsbar*s3[1].y;
    sout[1].y -= s1[1].y;    
    sout[1].z -= epsbar*s3[1].z;
    sout[1].z -= s1[1].z;    
    sout[1].w -= epsbar*s3[1].w;
    sout[1].w -= s1[1].w;       
    
    sout[2].x -= epsbar*s3[2].x;
    sout[2].x -= s1[2].x;
    sout[2].y -= epsbar*s3[2].y;
    sout[2].y -= s1[2].y;    
    sout[2].z -= epsbar*s3[2].z;
    sout[2].z -= s1[2].z;    
    sout[2].w -= epsbar*s3[2].w;
    sout[2].w -= s1[2].w;       
    
    sout[3].x -= epsbar*s3[3].x;
    sout[3].x -= s1[3].x;
    sout[3].y -= epsbar*s3[3].y;
    sout[3].y -= s1[3].y;    
    sout[3].z -= epsbar*s3[3].z;
    sout[3].z -= s1[3].z;    
    sout[3].w -= epsbar*s3[3].w;
    sout[3].w -= s1[3].w;   
    
    sout[4].x -= epsbar*s3[4].x;
    sout[4].x -= s1[4].x;
    sout[4].y -= epsbar*s3[4].y;
    sout[4].y -= s1[4].y;    
    sout[4].z -= epsbar*s3[4].z;
    sout[4].z -= s1[4].z;    
    sout[4].w -= epsbar*s3[4].w;
    sout[4].w -= s1[4].w;       
    
    sout[5].x -= epsbar*s3[5].x;
    sout[5].x -= s1[5].x;
    sout[5].y -= epsbar*s3[5].y;
    sout[5].y -= s1[5].y;    
    sout[5].z -= epsbar*s3[5].z;
    sout[5].z -= s1[5].z;    
    sout[5].w -= epsbar*s3[5].w;
    sout[5].w -= s1[5].w;       
    
    #ifdef RELATIVISTIC_BASIS
      dev_Gamma5_assign_rel(&(s1[0]), &(sout[0]));
    #else
      dev_Gamma5_assign(&(s1[0]), &(sout[0]));    
    #endif   
    dev_write_spinor(&(s1[0]),&(s2_up[pos]));    
    
    
    //lower output spinor
    dev_read_spinor(&(sout[0]), &(s2_dn[pos]));
    dev_read_spinor(&(s1[0]), &(s1_dn[pos]));    
    dev_read_spinor(&(s3[0]), &(s3_dn[pos]));    
    

     sout[0].x -= epsbar*s3[0].x;
    sout[0].x -= s1[0].x;
    sout[0].y -= epsbar*s3[0].y;
    sout[0].y -= s1[0].y;    
    sout[0].z -= epsbar*s3[0].z;
    sout[0].z -= s1[0].z;    
    sout[0].w -= epsbar*s3[0].w;
    sout[0].w -= s1[0].w;       
    
    sout[1].x -= epsbar*s3[1].x;
    sout[1].x -= s1[1].x;
    sout[1].y -= epsbar*s3[1].y;
    sout[1].y -= s1[1].y;    
    sout[1].z -= epsbar*s3[1].z;
    sout[1].z -= s1[1].z;    
    sout[1].w -= epsbar*s3[1].w;
    sout[1].w -= s1[1].w;       
    
    sout[2].x -= epsbar*s3[2].x;
    sout[2].x -= s1[2].x;
    sout[2].y -= epsbar*s3[2].y;
    sout[2].y -= s1[2].y;    
    sout[2].z -= epsbar*s3[2].z;
    sout[2].z -= s1[2].z;    
    sout[2].w -= epsbar*s3[2].w;
    sout[2].w -= s1[2].w;       
    
    sout[3].x -= epsbar*s3[3].x;
    sout[3].x -= s1[3].x;
    sout[3].y -= epsbar*s3[3].y;
    sout[3].y -= s1[3].y;    
    sout[3].z -= epsbar*s3[3].z;
    sout[3].z -= s1[3].z;    
    sout[3].w -= epsbar*s3[3].w;
    sout[3].w -= s1[3].w;   
    
    sout[4].x -= epsbar*s3[4].x;
    sout[4].x -= s1[4].x;
    sout[4].y -= epsbar*s3[4].y;
    sout[4].y -= s1[4].y;    
    sout[4].z -= epsbar*s3[4].z;
    sout[4].z -= s1[4].z;    
    sout[4].w -= epsbar*s3[4].w;
    sout[4].w -= s1[4].w;       
    
    sout[5].x -= epsbar*s3[5].x;
    sout[5].x -= s1[5].x;
    sout[5].y -= epsbar*s3[5].y;
    sout[5].y -= s1[5].y;    
    sout[5].z -= epsbar*s3[5].z;
    sout[5].z -= s1[5].z;    
    sout[5].w -= epsbar*s3[5].w;
    sout[5].w -= s1[5].w;         
    
    #ifdef RELATIVISTIC_BASIS
      dev_Gamma5_assign_rel(&(s1[0]), &(sout[0]));
    #else
      dev_Gamma5_assign(&(s1[0]), &(sout[0]));    
    #endif    
    dev_write_spinor(&(s1[0]),&(s2_dn[pos]));       
  }
}



//sout_up = gamma5*((1+imubar gamma5) s2_dn - epsbar s2_up - s1_up)
//sout_dn = gamma5*((1+imubar gamma5) s2_up - epsbar s2_dn - s1_dn)
__global__ void dev_nd_linalg3_gamma5 (dev_spinor * s1_up, dev_spinor * s1_dn,
				dev_spinor * s2_up, dev_spinor * s2_dn,  
				dev_spinor * sout_up, dev_spinor * sout_dn, float epsbar, float sign
				) {
  dev_spinor s1[6], s3[6], sout[6], slocal[6], save_up[6], save_dn[6];	
  int pos = threadIdx.x + blockDim.x*blockIdx.x;
  dev_complex pm_imu = dev_initcomplex(0.0, sign * mubar);
  
  if (pos < dev_VOLUME) {
    //upper output spinor

    dev_read_spinor(&(s1[0]), &(s1_up[pos]));    

    #ifdef USETEXTURE
      dev_read_spinor_tex_up(&(s3[0]), pos);  
      dev_read_spinor_tex_dn(&(slocal[0]), pos);      
    #else
      dev_read_spinor(&(s3[0]), &(s2_up[pos]));  
      dev_read_spinor(&(slocal[0]), &(s2_dn[pos])); 
    #endif
    //store s2_up/dn in local variables
    //dev_copy_spinor_local(&(slocal[0]), &(save_dn[0]));
    //dev_copy_spinor_local(&(s3[0]), &(save_up[0]));    
    
    #ifdef RELATIVISTIC_BASIS
      dev_skalarmult_gamma5_spinor_rel(&(sout[0]), pm_imu, &(slocal[0]) );
    #else
      dev_skalarmult_gamma5_spinor(&(sout[0]), pm_imu, &(slocal[0]) );
    #endif
    dev_add_spinor_assign(&(sout[0]), &(slocal[0]));
    
 
    
    sout[0].x -= epsbar*s3[0].x;
    sout[0].x -= s1[0].x;
    sout[0].y -= epsbar*s3[0].y;
    sout[0].y -= s1[0].y;    
    sout[0].z -= epsbar*s3[0].z;
    sout[0].z -= s1[0].z;    
    sout[0].w -= epsbar*s3[0].w;
    sout[0].w -= s1[0].w;       
    
    sout[1].x -= epsbar*s3[1].x;
    sout[1].x -= s1[1].x;
    sout[1].y -= epsbar*s3[1].y;
    sout[1].y -= s1[1].y;    
    sout[1].z -= epsbar*s3[1].z;
    sout[1].z -= s1[1].z;    
    sout[1].w -= epsbar*s3[1].w;
    sout[1].w -= s1[1].w;       
    
    sout[2].x -= epsbar*s3[2].x;
    sout[2].x -= s1[2].x;
    sout[2].y -= epsbar*s3[2].y;
    sout[2].y -= s1[2].y;    
    sout[2].z -= epsbar*s3[2].z;
    sout[2].z -= s1[2].z;    
    sout[2].w -= epsbar*s3[2].w;
    sout[2].w -= s1[2].w;       
    
    sout[3].x -= epsbar*s3[3].x;
    sout[3].x -= s1[3].x;
    sout[3].y -= epsbar*s3[3].y;
    sout[3].y -= s1[3].y;    
    sout[3].z -= epsbar*s3[3].z;
    sout[3].z -= s1[3].z;    
    sout[3].w -= epsbar*s3[3].w;
    sout[3].w -= s1[3].w;   
    
    sout[4].x -= epsbar*s3[4].x;
    sout[4].x -= s1[4].x;
    sout[4].y -= epsbar*s3[4].y;
    sout[4].y -= s1[4].y;    
    sout[4].z -= epsbar*s3[4].z;
    sout[4].z -= s1[4].z;    
    sout[4].w -= epsbar*s3[4].w;
    sout[4].w -= s1[4].w;       
    
    sout[5].x -= epsbar*s3[5].x;
    sout[5].x -= s1[5].x;
    sout[5].y -= epsbar*s3[5].y;
    sout[5].y -= s1[5].y;    
    sout[5].z -= epsbar*s3[5].z;
    sout[5].z -= s1[5].z;    
    sout[5].w -= epsbar*s3[5].w;
    sout[5].w -= s1[5].w;       
    
    #ifdef RELATIVISTIC_BASIS
      dev_Gamma5_assign_rel(&(s1[0]), &(sout[0]));
    #else
      dev_Gamma5_assign(&(s1[0]), &(sout[0]));    
    #endif   
    dev_write_spinor(&(s1[0]),&(sout_up[pos]));    
    
    
    //lower output spinor
    
    pm_imu = dev_initcomplex(0.0, -sign * mubar);

    dev_read_spinor(&(s1[0]), &(s1_dn[pos]));    
    #ifdef USETEXTURE
      dev_read_spinor_tex_dn(&(s3[0]), pos);  
      dev_read_spinor_tex_up(&(slocal[0]), pos);      
    #else
      dev_read_spinor(&(s3[0]), &(s2_dn[pos]));  
      dev_read_spinor(&(slocal[0]), &(s2_up[pos]));    
    #endif
    //restore s2_up/dn from local variables
    //dev_copy_spinor_local(&(save_dn[0]), &(s3[0]));
    //dev_copy_spinor_local(&(save_up[0]), &(slocal[0])); 
    
    #ifdef RELATIVISTIC_BASIS
      dev_skalarmult_gamma5_spinor_rel(&(sout[0]), pm_imu, &(slocal[0]) );
    #else
      dev_skalarmult_gamma5_spinor(&(sout[0]), pm_imu, &(slocal[0]) );
    #endif
    dev_add_spinor_assign(&(sout[0]), &(slocal[0]));
    

    

    sout[0].x -= epsbar*s3[0].x;
    sout[0].x -= s1[0].x;
    sout[0].y -= epsbar*s3[0].y;
    sout[0].y -= s1[0].y;    
    sout[0].z -= epsbar*s3[0].z;
    sout[0].z -= s1[0].z;    
    sout[0].w -= epsbar*s3[0].w;
    sout[0].w -= s1[0].w;       
    
    sout[1].x -= epsbar*s3[1].x;
    sout[1].x -= s1[1].x;
    sout[1].y -= epsbar*s3[1].y;
    sout[1].y -= s1[1].y;    
    sout[1].z -= epsbar*s3[1].z;
    sout[1].z -= s1[1].z;    
    sout[1].w -= epsbar*s3[1].w;
    sout[1].w -= s1[1].w;       
    
    sout[2].x -= epsbar*s3[2].x;
    sout[2].x -= s1[2].x;
    sout[2].y -= epsbar*s3[2].y;
    sout[2].y -= s1[2].y;    
    sout[2].z -= epsbar*s3[2].z;
    sout[2].z -= s1[2].z;    
    sout[2].w -= epsbar*s3[2].w;
    sout[2].w -= s1[2].w;       
    
    sout[3].x -= epsbar*s3[3].x;
    sout[3].x -= s1[3].x;
    sout[3].y -= epsbar*s3[3].y;
    sout[3].y -= s1[3].y;    
    sout[3].z -= epsbar*s3[3].z;
    sout[3].z -= s1[3].z;    
    sout[3].w -= epsbar*s3[3].w;
    sout[3].w -= s1[3].w;   
    
    sout[4].x -= epsbar*s3[4].x;
    sout[4].x -= s1[4].x;
    sout[4].y -= epsbar*s3[4].y;
    sout[4].y -= s1[4].y;    
    sout[4].z -= epsbar*s3[4].z;
    sout[4].z -= s1[4].z;    
    sout[4].w -= epsbar*s3[4].w;
    sout[4].w -= s1[4].w;       
    
    sout[5].x -= epsbar*s3[5].x;
    sout[5].x -= s1[5].x;
    sout[5].y -= epsbar*s3[5].y;
    sout[5].y -= s1[5].y;    
    sout[5].z -= epsbar*s3[5].z;
    sout[5].z -= s1[5].z;    
    sout[5].w -= epsbar*s3[5].w;
    sout[5].w -= s1[5].w;         
    
    #ifdef RELATIVISTIC_BASIS
      dev_Gamma5_assign_rel(&(s1[0]), &(sout[0]));
    #else
      dev_Gamma5_assign(&(s1[0]), &(sout[0]));    
    #endif    
    dev_write_spinor(&(s1[0]),&(sout_dn[pos]));   
    
  }
}






//sout_up = gamma5*((1+imubar gamma5) s2_up - epsbar s2_dn - s1_up)
//sout_dn = gamma5*((1+imubar gamma5) s2_dn - epsbar s2_up - s1_dn)
__global__ void dev_nd_linalg4_gamma5 (dev_spinor * s1_up, dev_spinor * s1_dn,
				dev_spinor * s2_up, dev_spinor * s2_dn,  
				dev_spinor * sout_up, dev_spinor * sout_dn, float epsbar, float sign
				) {
  dev_spinor s1[6], s3[6], sout[6], slocal[6], save_up[6], save_dn[6];	
  int pos = threadIdx.x + blockDim.x*blockIdx.x;
  dev_complex pm_imu = dev_initcomplex(0.0, sign * mubar);
  
  if (pos < dev_VOLUME) {
    //upper output spinor

    dev_read_spinor(&(s1[0]), &(s1_up[pos]));    
    #ifdef USETEXTURE
      dev_read_spinor_tex_dn(&(s3[0]), pos);  
      dev_read_spinor_tex_up(&(slocal[0]), pos);      
    #else
      dev_read_spinor(&(slocal[0]), &(s2_up[pos]));    
      dev_read_spinor(&(s3[0]), &(s2_dn[pos]));  
    #endif 
    //store s2_up/dn in local variables
    //dev_copy_spinor_local(&(slocal[0]), &(save_up[0]));
    //dev_copy_spinor_local(&(s3[0]), &(save_dn[0]));       
    
    #ifdef RELATIVISTIC_BASIS
      dev_skalarmult_gamma5_spinor_rel(&(sout[0]), pm_imu, &(slocal[0]) );
    #else
      dev_skalarmult_gamma5_spinor(&(sout[0]), pm_imu, &(slocal[0]) );
    #endif
    dev_add_spinor_assign(&(sout[0]), &(slocal[0]));
    
 
    
    sout[0].x -= epsbar*s3[0].x;
    sout[0].x -= s1[0].x;
    sout[0].y -= epsbar*s3[0].y;
    sout[0].y -= s1[0].y;    
    sout[0].z -= epsbar*s3[0].z;
    sout[0].z -= s1[0].z;    
    sout[0].w -= epsbar*s3[0].w;
    sout[0].w -= s1[0].w;       
    
    sout[1].x -= epsbar*s3[1].x;
    sout[1].x -= s1[1].x;
    sout[1].y -= epsbar*s3[1].y;
    sout[1].y -= s1[1].y;    
    sout[1].z -= epsbar*s3[1].z;
    sout[1].z -= s1[1].z;    
    sout[1].w -= epsbar*s3[1].w;
    sout[1].w -= s1[1].w;       
    
    sout[2].x -= epsbar*s3[2].x;
    sout[2].x -= s1[2].x;
    sout[2].y -= epsbar*s3[2].y;
    sout[2].y -= s1[2].y;    
    sout[2].z -= epsbar*s3[2].z;
    sout[2].z -= s1[2].z;    
    sout[2].w -= epsbar*s3[2].w;
    sout[2].w -= s1[2].w;       
    
    sout[3].x -= epsbar*s3[3].x;
    sout[3].x -= s1[3].x;
    sout[3].y -= epsbar*s3[3].y;
    sout[3].y -= s1[3].y;    
    sout[3].z -= epsbar*s3[3].z;
    sout[3].z -= s1[3].z;    
    sout[3].w -= epsbar*s3[3].w;
    sout[3].w -= s1[3].w;   
    
    sout[4].x -= epsbar*s3[4].x;
    sout[4].x -= s1[4].x;
    sout[4].y -= epsbar*s3[4].y;
    sout[4].y -= s1[4].y;    
    sout[4].z -= epsbar*s3[4].z;
    sout[4].z -= s1[4].z;    
    sout[4].w -= epsbar*s3[4].w;
    sout[4].w -= s1[4].w;       
    
    sout[5].x -= epsbar*s3[5].x;
    sout[5].x -= s1[5].x;
    sout[5].y -= epsbar*s3[5].y;
    sout[5].y -= s1[5].y;    
    sout[5].z -= epsbar*s3[5].z;
    sout[5].z -= s1[5].z;    
    sout[5].w -= epsbar*s3[5].w;
    sout[5].w -= s1[5].w;       
    
    #ifdef RELATIVISTIC_BASIS
      dev_Gamma5_assign_rel(&(s1[0]), &(sout[0]));
    #else
      dev_Gamma5_assign(&(s1[0]), &(sout[0]));    
    #endif   
    dev_write_spinor(&(s1[0]),&(sout_up[pos]));    
    
    
    //lower output spinor
    
    pm_imu = dev_initcomplex(0.0, -sign * mubar);

    dev_read_spinor(&(s1[0]), &(s1_dn[pos]));
    #ifdef USETEXTURE
      dev_read_spinor_tex_dn(&(slocal[0]), pos);     
      dev_read_spinor_tex_up(&(s3[0]), pos);  
    #else    
      dev_read_spinor(&(slocal[0]), &(s2_dn[pos]));    
      dev_read_spinor(&(s3[0]), &(s2_up[pos]));  
    #endif
    //restore s2_up/dn from local variables
    //dev_copy_spinor_local(&(save_dn[0]), &(slocal[0]));
    //dev_copy_spinor_local(&(save_up[0]), &(s3[0]));     
    
    
    #ifdef RELATIVISTIC_BASIS
      dev_skalarmult_gamma5_spinor_rel(&(sout[0]), pm_imu, &(slocal[0]) );
    #else
      dev_skalarmult_gamma5_spinor(&(sout[0]), pm_imu, &(slocal[0]) );
    #endif
    dev_add_spinor_assign(&(sout[0]), &(slocal[0]));
    

    

    sout[0].x -= epsbar*s3[0].x;
    sout[0].x -= s1[0].x;
    sout[0].y -= epsbar*s3[0].y;
    sout[0].y -= s1[0].y;    
    sout[0].z -= epsbar*s3[0].z;
    sout[0].z -= s1[0].z;    
    sout[0].w -= epsbar*s3[0].w;
    sout[0].w -= s1[0].w;       
    
    sout[1].x -= epsbar*s3[1].x;
    sout[1].x -= s1[1].x;
    sout[1].y -= epsbar*s3[1].y;
    sout[1].y -= s1[1].y;    
    sout[1].z -= epsbar*s3[1].z;
    sout[1].z -= s1[1].z;    
    sout[1].w -= epsbar*s3[1].w;
    sout[1].w -= s1[1].w;       
    
    sout[2].x -= epsbar*s3[2].x;
    sout[2].x -= s1[2].x;
    sout[2].y -= epsbar*s3[2].y;
    sout[2].y -= s1[2].y;    
    sout[2].z -= epsbar*s3[2].z;
    sout[2].z -= s1[2].z;    
    sout[2].w -= epsbar*s3[2].w;
    sout[2].w -= s1[2].w;       
    
    sout[3].x -= epsbar*s3[3].x;
    sout[3].x -= s1[3].x;
    sout[3].y -= epsbar*s3[3].y;
    sout[3].y -= s1[3].y;    
    sout[3].z -= epsbar*s3[3].z;
    sout[3].z -= s1[3].z;    
    sout[3].w -= epsbar*s3[3].w;
    sout[3].w -= s1[3].w;   
    
    sout[4].x -= epsbar*s3[4].x;
    sout[4].x -= s1[4].x;
    sout[4].y -= epsbar*s3[4].y;
    sout[4].y -= s1[4].y;    
    sout[4].z -= epsbar*s3[4].z;
    sout[4].z -= s1[4].z;    
    sout[4].w -= epsbar*s3[4].w;
    sout[4].w -= s1[4].w;       
    
    sout[5].x -= epsbar*s3[5].x;
    sout[5].x -= s1[5].x;
    sout[5].y -= epsbar*s3[5].y;
    sout[5].y -= s1[5].y;    
    sout[5].z -= epsbar*s3[5].z;
    sout[5].z -= s1[5].z;    
    sout[5].w -= epsbar*s3[5].w;
    sout[5].w -= s1[5].w;         
    
    #ifdef RELATIVISTIC_BASIS
      dev_Gamma5_assign_rel(&(s1[0]), &(sout[0]));
    #else
      dev_Gamma5_assign(&(s1[0]), &(sout[0]));    
    #endif    
    dev_write_spinor(&(s1[0]),&(sout_dn[pos]));   
    
  }
}







//s2_up = nrm*(s2_up +epsbar s1_dn)
//s2_dn = nrm*(s2_dn +epsbar s1_up)
// __global__ void dev_nd_linalg2 (dev_spinor * s1_up, dev_spinor * s1_dn, 
// 				dev_spinor * s2_up, dev_spinor * s2_dn,	float epsbar, float nrm			
// 				) {
//   dev_spinor s1[6], sout[6];	
//   int pos = threadIdx.x + blockDim.x*blockIdx.x;
//   
//   if (pos < dev_VOLUME) {
//     //upper output spinor
//     dev_read_spinor(&(sout[0]), &(s2_up[pos]));
//     #ifdef USETEXTURE
//       dev_read_spinor_tex_dn(&(s1[0]), pos);      
//     #else
//       dev_read_spinor(&(s1[0]), &(s1_dn[pos]));    
//     #endif
//     
//     sout[0].x += epsbar*s1[0].x;
//     sout[0].x *= nrm;
//     sout[0].y += epsbar*s1[0].y;
//     sout[0].y *= nrm;       
//     sout[0].z += epsbar*s1[0].z;
//     sout[0].z *= nrm;      
//     sout[0].w += epsbar*s1[0].w;
//     sout[0].w *= nrm;
//        
//     sout[1].x += epsbar*s1[1].x;
//     sout[1].x *= nrm;
//     sout[1].y += epsbar*s1[1].y;
//     sout[1].y *= nrm;       
//     sout[1].z += epsbar*s1[1].z;
//     sout[1].z *= nrm;      
//     sout[1].w += epsbar*s1[1].w;
//     sout[1].w *= nrm;    
//     
//     sout[2].x += epsbar*s1[2].x;
//     sout[2].x *= nrm;
//     sout[2].y += epsbar*s1[2].y;
//     sout[2].y *= nrm;       
//     sout[2].z += epsbar*s1[2].z;
//     sout[2].z *= nrm;      
//     sout[2].w += epsbar*s1[2].w;
//     sout[2].w *= nrm;    
// 
//     sout[3].x += epsbar*s1[3].x;
//     sout[3].x *= nrm;
//     sout[3].y += epsbar*s1[3].y;
//     sout[3].y *= nrm;       
//     sout[3].z += epsbar*s1[3].z;
//     sout[3].z *= nrm;      
//     sout[3].w += epsbar*s1[3].w;
//     sout[3].w *= nrm;    
// 
//     sout[4].x += epsbar*s1[4].x;
//     sout[4].x *= nrm;
//     sout[4].y += epsbar*s1[4].y;
//     sout[4].y *= nrm;       
//     sout[4].z += epsbar*s1[4].z;
//     sout[4].z *= nrm;      
//     sout[4].w += epsbar*s1[4].w;
//     sout[4].w *= nrm;    
//     
//     sout[5].x += epsbar*s1[5].x;
//     sout[5].x *= nrm;
//     sout[5].y += epsbar*s1[5].y;
//     sout[5].y *= nrm;       
//     sout[5].z += epsbar*s1[5].z;
//     sout[5].z *= nrm;      
//     sout[5].w += epsbar*s1[5].w;
//     sout[5].w *= nrm;    
//     
//     dev_write_spinor(&(sout[0]),&(s2_up[pos]));  
//     
//     
//     //upper output spinor
//     dev_read_spinor(&(sout[0]), &(s2_dn[pos]));
//     #ifdef USETEXTURE
//       dev_read_spinor_tex_up(&(s1[0]), pos);      
//     #else
//      dev_read_spinor(&(s1[0]), &(s1_up[pos]));    
//     #endif
//     
//     sout[0].x += epsbar*s1[0].x;
//     sout[0].x *= nrm;
//     sout[0].y += epsbar*s1[0].y;
//     sout[0].y *= nrm;       
//     sout[0].z += epsbar*s1[0].z;
//     sout[0].z *= nrm;      
//     sout[0].w += epsbar*s1[0].w;
//     sout[0].w *= nrm;
//        
//     sout[1].x += epsbar*s1[1].x;
//     sout[1].x *= nrm;
//     sout[1].y += epsbar*s1[1].y;
//     sout[1].y *= nrm;       
//     sout[1].z += epsbar*s1[1].z;
//     sout[1].z *= nrm;      
//     sout[1].w += epsbar*s1[1].w;
//     sout[1].w *= nrm;    
//     
//     sout[2].x += epsbar*s1[2].x;
//     sout[2].x *= nrm;
//     sout[2].y += epsbar*s1[2].y;
//     sout[2].y *= nrm;       
//     sout[2].z += epsbar*s1[2].z;
//     sout[2].z *= nrm;      
//     sout[2].w += epsbar*s1[2].w;
//     sout[2].w *= nrm;    
// 
//     sout[3].x += epsbar*s1[3].x;
//     sout[3].x *= nrm;
//     sout[3].y += epsbar*s1[3].y;
//     sout[3].y *= nrm;       
//     sout[3].z += epsbar*s1[3].z;
//     sout[3].z *= nrm;      
//     sout[3].w += epsbar*s1[3].w;
//     sout[3].w *= nrm;    
// 
//     sout[4].x += epsbar*s1[4].x;
//     sout[4].x *= nrm;
//     sout[4].y += epsbar*s1[4].y;
//     sout[4].y *= nrm;       
//     sout[4].z += epsbar*s1[4].z;
//     sout[4].z *= nrm;      
//     sout[4].w += epsbar*s1[4].w;
//     sout[4].w *= nrm;    
//     
//     sout[5].x += epsbar*s1[5].x;
//     sout[5].x *= nrm;
//     sout[5].y += epsbar*s1[5].y;
//     sout[5].y *= nrm;       
//     sout[5].z += epsbar*s1[5].z;
//     sout[5].z *= nrm;      
//     sout[5].w += epsbar*s1[5].w;
//     sout[5].w *= nrm;    
//     
// /*    s2_dn[pos+0*DEVOFF].x = sout[0].x;
//     s2_dn[pos+0*DEVOFF].y = sout[0].y;
//     s2_dn[pos+0*DEVOFF].z = sout[0].z;
//     s2_dn[pos+0*DEVOFF].w = sout[0].w;   
//     
//     s2_dn[pos+1*DEVOFF].x = sout[1].x;
//     s2_dn[pos+1*DEVOFF].y = sout[1].y;
//     s2_dn[pos+1*DEVOFF].z = sout[1].z;
//     s2_dn[pos+1*DEVOFF].w = sout[1].w;       
//     
//     s2_dn[pos+2*DEVOFF].x = sout[2].x;
//     s2_dn[pos+2*DEVOFF].y = sout[2].y;
//     s2_dn[pos+2*DEVOFF].z = sout[2].z;
//     s2_dn[pos+2*DEVOFF].w = sout[2].w;      
//     
//     s2_dn[pos+3*DEVOFF].x = sout[3].x;
//     s2_dn[pos+3*DEVOFF].y = sout[3].y;
//     s2_dn[pos+3*DEVOFF].z = sout[3].z;
//     s2_dn[pos+3*DEVOFF].w = sout[3].w;      
//     
//     s2_dn[pos+4*DEVOFF].x = sout[4].x;
//     s2_dn[pos+4*DEVOFF].y = sout[4].y;
//     s2_dn[pos+4*DEVOFF].z = sout[4].z;
//     s2_dn[pos+4*DEVOFF].w = sout[4].w;    
//     
//     s2_dn[pos+5*DEVOFF].x = sout[5].x;
//     s2_dn[pos+5*DEVOFF].y = sout[5].y;
//     s2_dn[pos+5*DEVOFF].z = sout[5].z;
//     s2_dn[pos+5*DEVOFF].w = sout[5].w;    */  
//     
//     dev_write_spinor(&(sout[0]),&(s2_dn[pos]));     
//     
//     
//     
//   }
// }
// 


__global__ void dev_nd_linalg2 (dev_spinor * s1_up, dev_spinor * s1_dn, 
				dev_spinor * s2_up, dev_spinor * s2_dn,	float epsbar, float nrm			
				) {
  dev_spinor s1[6], sout[6];	
  int pos = threadIdx.x + blockDim.x*blockIdx.x;
  
  if (pos < dev_VOLUME) {
    //upper output spinor
    dev_read_spinor(&(sout[0]), &(s2_up[pos]));
    #ifdef USETEXTURE
      dev_read_spinor_tex_dn(&(s1[0]), pos);      
    #else
      dev_read_spinor(&(s1[0]), &(s1_dn[pos]));    
    #endif
    
    sout[0].x += epsbar*s1[0].x;
    s2_up[pos+0*DEVOFF].x = sout[0].x*nrm;    
    sout[0].y += epsbar*s1[0].y;
    s2_up[pos+0*DEVOFF].y = sout[0].y*nrm;    
    sout[0].z += epsbar*s1[0].z;
    s2_up[pos+0*DEVOFF].z = sout[0].z*nrm;      
    sout[0].w += epsbar*s1[0].w;
    s2_up[pos+0*DEVOFF].w = sout[0].w*nrm; 
       
    sout[1].x += epsbar*s1[1].x;
    s2_up[pos+1*DEVOFF].x = sout[1].x*nrm; 
    sout[1].y += epsbar*s1[1].y;
    s2_up[pos+1*DEVOFF].y = sout[1].y*nrm;      
    sout[1].z += epsbar*s1[1].z;
    s2_up[pos+1*DEVOFF].z = sout[1].z*nrm;       
    sout[1].w += epsbar*s1[1].w;
    s2_up[pos+1*DEVOFF].w = sout[1].w*nrm;    
    
    sout[2].x += epsbar*s1[2].x;
    s2_up[pos+2*DEVOFF].x = sout[2].x*nrm; 
    sout[2].y += epsbar*s1[2].y;
    s2_up[pos+2*DEVOFF].y = sout[2].y*nrm;  
    sout[2].z += epsbar*s1[2].z;
    s2_up[pos+2*DEVOFF].z = sout[2].z*nrm;    
    sout[2].w += epsbar*s1[2].w;
    s2_up[pos+2*DEVOFF].w = sout[2].w*nrm;   

    sout[3].x += epsbar*s1[3].x;
    s2_up[pos+3*DEVOFF].x = sout[3].x*nrm;
    sout[3].y += epsbar*s1[3].y;
    s2_up[pos+3*DEVOFF].y = sout[3].y*nrm;   
    sout[3].z += epsbar*s1[3].z;
    s2_up[pos+3*DEVOFF].z = sout[3].z*nrm;     
    sout[3].w += epsbar*s1[3].w;
    s2_up[pos+3*DEVOFF].w = sout[3].w*nrm;   

    sout[4].x += epsbar*s1[4].x;
    s2_up[pos+4*DEVOFF].x = sout[4].x*nrm;
    sout[4].y += epsbar*s1[4].y;
    s2_up[pos+4*DEVOFF].y = sout[4].y*nrm;      
    sout[4].z += epsbar*s1[4].z;
    s2_up[pos+4*DEVOFF].z = sout[4].z*nrm;    
    sout[4].w += epsbar*s1[4].w;
    s2_up[pos+4*DEVOFF].w = sout[4].w*nrm;   
    
    sout[5].x += epsbar*s1[5].x;
    s2_up[pos+5*DEVOFF].x = sout[5].x*nrm;
    sout[5].y += epsbar*s1[5].y;
    s2_up[pos+5*DEVOFF].y = sout[5].y*nrm;      
    sout[5].z += epsbar*s1[5].z;
    s2_up[pos+5*DEVOFF].z = sout[5].z*nrm;   
    sout[5].w += epsbar*s1[5].w;
    s2_up[pos+5*DEVOFF].w = sout[5].w*nrm;   
    

    
    
    //upper output spinor
    dev_read_spinor(&(sout[0]), &(s2_dn[pos]));
    #ifdef USETEXTURE
      dev_read_spinor_tex_up(&(s1[0]), pos);      
    #else
     dev_read_spinor(&(s1[0]), &(s1_up[pos]));    
    #endif
    
    sout[0].x += epsbar*s1[0].x;
    s2_dn[pos+0*DEVOFF].x = sout[0].x*nrm;
    sout[0].y += epsbar*s1[0].y;
    s2_dn[pos+0*DEVOFF].y = sout[0].y*nrm;        
    sout[0].z += epsbar*s1[0].z;
    s2_dn[pos+0*DEVOFF].z = sout[0].z*nrm;             
    sout[0].w += epsbar*s1[0].w;
    s2_dn[pos+0*DEVOFF].w = sout[0].w*nrm;  
       
    sout[1].x += epsbar*s1[1].x;
    s2_dn[pos+1*DEVOFF].x = sout[1].x*nrm;    
    sout[1].y += epsbar*s1[1].y;
    s2_dn[pos+1*DEVOFF].y = sout[1].y*nrm; 
    sout[1].z += epsbar*s1[1].z;
    s2_dn[pos+1*DEVOFF].z = sout[1].z*nrm; 
    sout[1].w += epsbar*s1[1].w;
    s2_dn[pos+1*DEVOFF].w = sout[1].w*nrm; 
    
    sout[2].x += epsbar*s1[2].x;
    s2_dn[pos+2*DEVOFF].x = sout[2].x*nrm; 
    sout[2].y += epsbar*s1[2].y;
    s2_dn[pos+2*DEVOFF].y = sout[2].y*nrm;     
    sout[2].z += epsbar*s1[2].z;
    s2_dn[pos+2*DEVOFF].z = sout[2].z*nrm; 
    sout[2].w += epsbar*s1[2].w;
    s2_dn[pos+2*DEVOFF].w = sout[2].w*nrm;  

    sout[3].x += epsbar*s1[3].x;
    s2_dn[pos+3*DEVOFF].x = sout[3].x*nrm; 
    sout[3].y += epsbar*s1[3].y;
    s2_dn[pos+3*DEVOFF].y = sout[3].y*nrm;  
    sout[3].z += epsbar*s1[3].z;
    s2_dn[pos+3*DEVOFF].z = sout[3].z*nrm;      
    sout[3].w += epsbar*s1[3].w;
    s2_dn[pos+3*DEVOFF].w = sout[3].w*nrm;     

    sout[4].x += epsbar*s1[4].x;
    s2_dn[pos+4*DEVOFF].x = sout[4].x*nrm; 
    sout[4].y += epsbar*s1[4].y;
    s2_dn[pos+4*DEVOFF].y = sout[4].y*nrm;    
    sout[4].z += epsbar*s1[4].z;
    s2_dn[pos+4*DEVOFF].z = sout[4].z*nrm;      
    sout[4].w += epsbar*s1[4].w;
    s2_dn[pos+4*DEVOFF].w = sout[4].w*nrm;    
    
    sout[5].x += epsbar*s1[5].x;
    s2_dn[pos+5*DEVOFF].x = sout[5].x*nrm; 
    sout[5].y += epsbar*s1[5].y;
    s2_dn[pos+5*DEVOFF].y = sout[5].y*nrm; 
    sout[5].z += epsbar*s1[5].z;
    s2_dn[pos+5*DEVOFF].z = sout[5].z*nrm;    
    sout[5].w += epsbar*s1[5].w;
    s2_dn[pos+5*DEVOFF].w = sout[5].w*nrm;    

  }
}










///////////////////////////
// MATRIX MULTIPLICATION //
///////////////////////////

// the GPU implementation of  Qtm_pm_ndpsi(...)  from tm_operators_nd.c
// we have the possibility to add a constant shift > 0 
void dev_Qtm_pm_ndpsi_old (dev_spinor * spinout_up, dev_spinor * spinout_dn,
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
  //   MATCHING with  Qtm_pm_ndpsi  //
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
  
  // Hopping:
  #ifdef USETEXTURE
    bind_texture_spin(spinin_dn,1);
  #endif
    dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, spinin_dn, dev_spin_eo1_up, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0,0,N_sites);	// dev_spin_eo1_up = (M_eo) * spinin_dn
  #ifdef USETEXTURE
    unbind_texture_spin(1);
  #endif
  
  #ifdef USETEXTURE
    bind_texture_spin(spinin_up,1);
  #endif
    dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, spinin_up, dev_spin_eo1_dn, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0,0,N_sites);	// dev_spin_eo1_dn = (M_eo) * spinin_up
  #ifdef USETEXTURE
    unbind_texture_spin(1);
  #endif

  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo1_up, dev_spin_eo2_up, -1.0);		// dev_spin_eo2_up  =  (1 - imubar) * dev_spin_eo1_up  =  (1 - imubar)*(M_eo) * spinin_dn
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo1_dn, dev_spin_eo2_dn, +1.0);		// dev_spin_eo2_dn  =  (1 + imubar) * dev_spin_eo1_dn  =  (1 + imubar)*(M_eo) * spinin_up

  dev_nd_linalg2<<<gridsize2, blocksize2>>>(dev_spin_eo1_up, dev_spin_eo1_dn, dev_spin_eo2_up, dev_spin_eo2_dn, g_epsbar, nrm); 


  // Hopping:
  #ifdef USETEXTURE
    bind_texture_spin(dev_spin_eo2_up,1);									// remember: this is ((M_oe)(Mee^-1)(M_eo)) * (spinin_dn, spinin_up):
  #endif
    dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo2_up, dev_spin_eo1_up, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1,0,N_sites);
  #ifdef USETEXTURE
    unbind_texture_spin(1);							// dev_spin_eo1_up  =  (M_oe) * dev_spin_eo2_up  =  (M_oe)*nrm*(1-imubar)*(M_eo) * spinin_dn  +  (M_oe)*nrm*epsbar*(M_eo) * spinin_up
  #endif									
  
  #ifdef USETEXTURE
    bind_texture_spin(dev_spin_eo2_dn,1);
  #endif
    dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo2_dn, dev_spin_eo1_dn, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1,0,N_sites);
  #ifdef USETEXTURE
    unbind_texture_spin(1);							// dev_spin_eo1_dn  =  (M_oe) * dev_spin_eo2_dn  =  (M_oe)*nrm*(1+imubar)*(M_eo) * spinin_up  +  (M_oe)*nrm*epsbar*(M_eo) * spinin_dn
  #endif
  
 
  ////////////
  // (M_oo) //
  ////////////
  
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(spinin_dn, dev_spin_eo2_up, +1.0);			// dev_spin_eo2_up = (1 + imubar) * spinin_dn
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(spinin_up, dev_spin_eo2_dn, -1.0);			// dev_spin_eo2_dn = (1 - imubar) * spinin_up
  
  ///////////////////////////////////////
  // (M_oo)  -  (M_oe) (Mee^-1) (M_eo) //
  ///////////////////////////////////////
   
  dev_nd_linalg1_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo1_up, dev_spin_eo1_dn, spinin_up, spinin_dn, dev_spin_eo2_up, dev_spin_eo2_dn, g_epsbar ); 
    

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
  
  // Hopping:
  #ifdef USETEXTURE
    bind_texture_spin(dev_spin_eo3_up,1);
  #endif
    dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo3_up, dev_spin_eo1_up, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0,0,N_sites);	// dev_spin_eo1_up = (M_eo) * dev_spin_eo3_up
  #ifdef USETEXTURE
    unbind_texture_spin(1);
  #endif
  
  #ifdef USETEXTURE
    bind_texture_spin(dev_spin_eo3_dn,1);
  #endif
    dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo3_dn, dev_spin_eo1_dn, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0,0,N_sites);	// dev_spin_eo1_dn = (M_eo) * dev_spin_eo3_dn
  #ifdef USETEXTURE
    unbind_texture_spin(1);
  #endif
  
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo1_up, dev_spin_eo2_up, -1.0);			// dev_spin_eo2_up  =  (1 - imubar) * dev_spin_eo1_up  =  (1 - imubar)*(M_eo) * dev_spin_eo3_up
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo1_dn, dev_spin_eo2_dn, +1.0);			// dev_spin_eo2_dn  =  (1 + imubar) * dev_spin_eo1_dn  =  (1 + imubar)*(M_eo) * dev_spin_eo3_dn
  dev_nd_linalg2<<<gridsize2, blocksize2>>>(dev_spin_eo1_up, dev_spin_eo1_dn, dev_spin_eo2_up, dev_spin_eo2_dn, g_epsbar, nrm); 
  
  // Hopping:
  #ifdef USETEXTURE
    bind_texture_spin(dev_spin_eo2_up,1);									// remember: this is ((M_oe) (Mee^-1) (M_eo)) * (dev_spin_eo3_up, dev_spin_eo3_dn):
  #endif
    dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo2_up, dev_spin_eo1_up, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1,0,N_sites);
  #ifdef USETEXTURE
    unbind_texture_spin(1);							// dev_spin_eo1_up  =  (M_oe) * dev_spin_eo2_up  =  (M_oe)*nrm*(1-imubar)*(M_eo) * dev_spin_eo3_up  +  (M_oe)*nrm*epsbar*(M_eo) * dev_spin_eo3_dn
  #endif
  
  #ifdef USETEXTURE
    bind_texture_spin(dev_spin_eo2_dn,1);
  #endif
    dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo2_dn, dev_spin_eo1_dn, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1,0,N_sites);
  #ifdef USETEXTURE
    unbind_texture_spin(1);							// dev_spin_eo1_dn  =  (M_oe) * dev_spin_eo2_dn  =  (M_oe)*nrm*(1+imubar)*(M_eo) * dev_spin_eo3_dn  +  (M_oe)*nrm*epsbar*(M_eo) * dev_spin_eo3_up
  #endif
  
  ////////////
  // (M_oo) //
  ////////////
  
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo3_up, dev_spin_eo2_up, +1.0);				// dev_spin_eo2_up = (1 + imubar) * dev_spin_eo3_up
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo3_dn, dev_spin_eo2_dn, -1.0);				// dev_spin_eo2_dn = (1 - imubar) * dev_spin_eo3_dn
  
  //dev_spin_eo3_up <-> dev_spin_eo3_dn here!! 
  dev_nd_linalg1_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo1_up, dev_spin_eo1_dn, dev_spin_eo3_dn, dev_spin_eo3_up, dev_spin_eo2_up, dev_spin_eo2_dn, g_epsbar ); 
  
  
  ////////////
  // output //	output is already done by setting  dev_spin_eo2_up/dn = spinout_up/dn
  ////////////
 
  return;
  
}//dev_Qtm_pm_ndpsi_old()









void dev_Qtm_pm_ndpsi (dev_spinor * spinout_up, dev_spinor * spinout_dn,
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
  //   MATCHING with  Qtm_pm_ndpsi  //
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
  
  // Hopping:
  #ifdef USETEXTURE
    bind_texture_spin(spinin_dn,1);
  #endif
    dev_Hopping_Matrix_ext<<<gridsize1, blocksize1, 0, stream_nd[0]>>>(dev_gf, spinin_dn, dev_spin_eo1_up, dev_spin_eo2_up, -1.0, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0,0,N_sites);	// dev_spin_eo1_up = (M_eo) * spinin_dn
  #ifdef USETEXTURE
    unbind_texture_spin(1);
  #endif
  
  #ifdef USETEXTURE
    bind_texture_spin(spinin_up,1);
  #endif
    dev_Hopping_Matrix_ext<<<gridsize1, blocksize1, 0, stream_nd[1]>>>(dev_gf, spinin_up, dev_spin_eo1_dn, dev_spin_eo2_dn, 1.0, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0,0,N_sites);	// dev_spin_eo1_dn = (M_eo) * spinin_up
  #ifdef USETEXTURE
    unbind_texture_spin(1);
  #endif

  cudaStreamSynchronize(stream_nd[0]);
  cudaStreamSynchronize(stream_nd[1]);  

  #ifdef USETEXTURE
    bind_texture_spin(dev_spin_eo1_up,1);
    bind_texture_spin_dn(dev_spin_eo1_dn,1);   
  #endif  
  dev_nd_linalg2<<<gridsize2, blocksize2>>>(dev_spin_eo1_up, dev_spin_eo1_dn, dev_spin_eo2_up, dev_spin_eo2_dn, g_epsbar, nrm); 
  #ifdef USETEXTURE
    unbind_texture_spin(1);
    unbind_texture_spin_dn(1);    
  #endif
  

  // Hopping:
  #ifdef USETEXTURE
    bind_texture_spin(dev_spin_eo2_up,1);									// remember: this is ((M_oe)(Mee^-1)(M_eo)) * (spinin_dn, spinin_up):
  #endif
    dev_Hopping_Matrix<<<gridsize1, blocksize1, 0, stream_nd[0]>>>(dev_gf, dev_spin_eo2_up, dev_spin_eo1_up, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1,0,N_sites);
  #ifdef USETEXTURE
    unbind_texture_spin(1);							// dev_spin_eo1_up  =  (M_oe) * dev_spin_eo2_up  =  (M_oe)*nrm*(1-imubar)*(M_eo) * spinin_dn  +  (M_oe)*nrm*epsbar*(M_eo) * spinin_up
  #endif									
  
  #ifdef USETEXTURE
    bind_texture_spin(dev_spin_eo2_dn,1);
  #endif
    dev_Hopping_Matrix<<<gridsize1, blocksize1, 0, stream_nd[1]>>>(dev_gf, dev_spin_eo2_dn, dev_spin_eo1_dn, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1,0,N_sites);
  #ifdef USETEXTURE
    unbind_texture_spin(1);							// dev_spin_eo1_dn  =  (M_oe) * dev_spin_eo2_dn  =  (M_oe)*nrm*(1+imubar)*(M_eo) * spinin_up  +  (M_oe)*nrm*epsbar*(M_eo) * spinin_dn
  #endif
  
 
  ////////////
  // (M_oo) //
  ////////////
  
  ///////////////////////////////////////
  // (M_oo)  -  (M_oe) (Mee^-1) (M_eo) //
  ///////////////////////////////////////  
  
  cudaStreamSynchronize(stream_nd[0]);
  cudaStreamSynchronize(stream_nd[1]);  
  
  #ifdef USETEXTURE
    bind_texture_spin(spinin_up,1);
    bind_texture_spin_dn(spinin_dn,1);   
  #endif     
  dev_nd_linalg3_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo1_up, dev_spin_eo1_dn, spinin_up, spinin_dn, dev_spin_eo2_up, dev_spin_eo2_dn, g_epsbar, +1.0 ); 
  #ifdef USETEXTURE
    unbind_texture_spin(1);
    unbind_texture_spin_dn(1);    
  #endif 
 
  

  ////////////////////
  // (a,b) -> (b,a) //
  ////////////////////
  dev_copy_spinor_field<<<gridsize4, blocksize4, 0, stream_nd[0]>>>(dev_spin_eo2_dn, dev_spin_eo3_up);			// dev_spin_eo3_up = dev_spin_eo2_dn
  dev_copy_spinor_field<<<gridsize4, blocksize4, 0, stream_nd[1]>>>(dev_spin_eo2_up, dev_spin_eo3_dn);			// dev_spin_eo3_dn = dev_spin_eo2_up

  ///////////////////////////////////								///////////////////////
  //        Q_tilde(2x2)           //								// (Q_tilde) * (b,a) //
  ///////////////////////////////////								///////////////////////
  
  ////////////////////////////
  // (M_oe) (Mee^-1) (M_eo) //
  ////////////////////////////
  
  // Hopping:
  #ifdef USETEXTURE
    bind_texture_spin(dev_spin_eo3_up,1);
  #endif
    dev_Hopping_Matrix_ext<<<gridsize1, blocksize1, 0, stream_nd[0]>>>(dev_gf, dev_spin_eo3_up, dev_spin_eo1_up, dev_spin_eo2_up, -1.0, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0,0,N_sites);	// dev_spin_eo1_up = (M_eo) * dev_spin_eo3_up
  #ifdef USETEXTURE
    unbind_texture_spin(1);
  #endif
  
  #ifdef USETEXTURE
    bind_texture_spin(dev_spin_eo3_dn,1);
  #endif
    dev_Hopping_Matrix_ext<<<gridsize1, blocksize1, 0, stream_nd[1]>>>(dev_gf, dev_spin_eo3_dn, dev_spin_eo1_dn, dev_spin_eo2_dn, +1.0, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0,0,N_sites);	// dev_spin_eo1_dn = (M_eo) * dev_spin_eo3_dn
  #ifdef USETEXTURE
    unbind_texture_spin(1);
  #endif
  
  cudaStreamSynchronize(stream_nd[0]);
  cudaStreamSynchronize(stream_nd[1]);  
  
  #ifdef USETEXTURE
    bind_texture_spin(dev_spin_eo1_up,1);
    bind_texture_spin_dn(dev_spin_eo1_dn,1);   
  #endif    
  dev_nd_linalg2<<<gridsize2, blocksize2>>>(dev_spin_eo1_up, dev_spin_eo1_dn, dev_spin_eo2_up, dev_spin_eo2_dn, g_epsbar, nrm); 
  #ifdef USETEXTURE
    unbind_texture_spin(1);
    unbind_texture_spin_dn(1);    
  #endif 
  
  // Hopping:
  #ifdef USETEXTURE
    bind_texture_spin(dev_spin_eo2_up,1);									// remember: this is ((M_oe) (Mee^-1) (M_eo)) * (dev_spin_eo3_up, dev_spin_eo3_dn):
  #endif
    dev_Hopping_Matrix<<<gridsize1, blocksize1, 0, stream_nd[0]>>>(dev_gf, dev_spin_eo2_up, dev_spin_eo1_up, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1,0,N_sites);
  #ifdef USETEXTURE
    unbind_texture_spin(1);							// dev_spin_eo1_up  =  (M_oe) * dev_spin_eo2_up  =  (M_oe)*nrm*(1-imubar)*(M_eo) * dev_spin_eo3_up  +  (M_oe)*nrm*epsbar*(M_eo) * dev_spin_eo3_dn
  #endif
  
  #ifdef USETEXTURE
    bind_texture_spin(dev_spin_eo2_dn,1);
  #endif
    dev_Hopping_Matrix<<<gridsize1, blocksize1, 0, stream_nd[1]>>>(dev_gf, dev_spin_eo2_dn, dev_spin_eo1_dn, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1,0,N_sites);
  #ifdef USETEXTURE
    unbind_texture_spin(1);							// dev_spin_eo1_dn  =  (M_oe) * dev_spin_eo2_dn  =  (M_oe)*nrm*(1+imubar)*(M_eo) * dev_spin_eo3_dn  +  (M_oe)*nrm*epsbar*(M_eo) * dev_spin_eo3_up
  #endif
  
  ////////////
  // (M_oo) //
  ////////////
  cudaStreamSynchronize(stream_nd[0]);
  cudaStreamSynchronize(stream_nd[1]);   

  #ifdef USETEXTURE
    bind_texture_spin(dev_spin_eo3_up,1);
    bind_texture_spin_dn(dev_spin_eo3_dn,1);   
  #endif    
    dev_nd_linalg4_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo1_up, dev_spin_eo1_dn, dev_spin_eo3_up, dev_spin_eo3_dn, dev_spin_eo2_up, dev_spin_eo2_dn, g_epsbar, +1.0 ); 
  #ifdef USETEXTURE
    unbind_texture_spin(1);
    unbind_texture_spin_dn(1);    
  #endif 
  
  cublasSscal (N_floats, phmc_invmaxev*phmc_invmaxev, (float *) dev_spin_eo2_up, 1);
  cublasSscal (N_floats, phmc_invmaxev*phmc_invmaxev, (float *) dev_spin_eo2_dn, 1);
  
  ////////////
  // output //	output is already done by setting  dev_spin_eo2_up/dn = spinout_up/dn
  ////////////
 
  return;
  
}//dev_Qtm_pm_ndpsi()









///////////////////////////
// MATRIX MULTIPLICATION //
///////////////////////////

// the GPU implementation of  Qtm_pm_ndpsi(...)  from tm_operators_nd.c
// we have the possibility to add a constant shift > 0 
void dev_Qtm_pm_ndpsi_updn (dev_spinor * spinout_up, dev_spinor * spinout_dn,
                              dev_spinor * spinin_up , dev_spinor * spinin_dn , 
                              int gridsize1, int blocksize1, int gridsize2, int blocksize2,
                              int gridsize3, int blocksize3, int gridsize4, int blocksize4) {
  
  
  // we will use the auxiliary fields  dev_spin_eo{1,2}_up/dn  for working on and buffering
  // and set  dev_spin_eo2_up/dn  equal  spinout_up/dn
  // spinin_up/dn  have to remain unchanged !!
  // spinout_up/dn  can be freely used
  

  int N_sites  =    VOLUME/2;		// #lattice sites
  int N_floats = 24*VOLUME/2;		// #floats

  dev_spin_eo2_up = spinout_up;		// need no memory allocated
  dev_spin_eo2_dn = spinout_dn;

 
  double nrm = 1.0 / (1.0 + g_mubar*g_mubar - g_epsbar*g_epsbar);
  
  #ifdef USETEXTURE
    bind_texture_spin(spinin_dn,1);
    bind_texture_spin_dn(spinin_up,1);    
  #endif
    dev_Hopping_Matrix_updn<<<gridsize1, blocksize1>>>(dev_gf, spinin_dn, spinin_up, dev_spin_eo1_up, dev_spin_eo1_dn,dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0,0,N_sites);	// dev_spin_eo1_up = (M_eo) * spinin_dn
  #ifdef USETEXTURE
    unbind_texture_spin(1);
    unbind_texture_spin_dn(1);   
  #endif
  
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo1_up, dev_spin_eo2_up, -1.0);		// dev_spin_eo2_up  =  (1 - imubar) * dev_spin_eo1_up  =  (1 - imubar)*(M_eo) * spinin_dn
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo1_dn, dev_spin_eo2_dn, +1.0);		// dev_spin_eo2_dn  =  (1 + imubar) * dev_spin_eo1_dn  =  (1 + imubar)*(M_eo) * spinin_up

  /*
  cublasSaxpy (N_floats, g_epsbar, (float *) dev_spin_eo1_dn, 1, (float *) dev_spin_eo2_up, 1);
  										// dev_spin_eo2_up  =  dev_spin_eo2_up  +  epsbar * dev_spin_eo1_dn  =  (1 - imubar)*(M_eo) * spinin_dn  +  epsbar * (M_eo) * spinin_up
  cublasSaxpy (N_floats, g_epsbar, (float *) dev_spin_eo1_up, 1, (float *) dev_spin_eo2_dn, 1);
  										// dev_spin_eo2_dn  =  dev_spin_eo2_dn  +  epsbar * dev_spin_eo1_up  =  (1 + imubar)*(M_eo) * spinin_up  +  epsbar * (M_eo) * spinin_dn
  

  cublasSscal (N_floats, nrm, (float *) dev_spin_eo2_up, 1);			// dev_spin_eo2_up  =  nrm * dev_spin_eo2_up  =  nrm*(1-imubar)*(M_eo) * spinin_dn  +  nrm*epsbar*(M_eo) * spinin_up
  
  cublasSscal (N_floats, nrm, (float *) dev_spin_eo2_dn, 1);			// dev_spin_eo2_dn  =  nrm * dev_spin_eo2_dn  =  nrm*(1+imubar)*(M_eo) * spinin_up  +  nrm*epsbar*(M_eo) * spinin_dn
  */
  dev_nd_linalg2<<<gridsize2, blocksize2>>>(dev_spin_eo1_up, dev_spin_eo1_dn, dev_spin_eo2_up, dev_spin_eo2_dn, g_epsbar, nrm); 

  #ifdef USETEXTURE
    bind_texture_spin(dev_spin_eo2_up,1);									// remember: this is ((M_oe)(Mee^-1)(M_eo)) * (spinin_dn, spinin_up):
    bind_texture_spin_dn(dev_spin_eo2_dn,1);									// remember: this is ((M_oe)(Mee^-1)(M_eo)) * (spinin_dn, spinin_up):    
  #endif
    dev_Hopping_Matrix_updn<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo2_up, dev_spin_eo2_dn, dev_spin_eo1_up, dev_spin_eo1_dn, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1,0,N_sites);
  #ifdef USETEXTURE
    unbind_texture_spin(1);	
    unbind_texture_spin_dn(1);
  #endif									

  
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(spinin_dn, dev_spin_eo2_up, +1.0);			// dev_spin_eo2_up = (1 + imubar) * spinin_dn
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(spinin_up, dev_spin_eo2_dn, -1.0);			// dev_spin_eo2_dn = (1 - imubar) * spinin_up
  
  /*
  cublasSaxpy (N_floats, -g_epsbar, (float *) spinin_up, 1, (float *) dev_spin_eo2_up, 1);
  												// dev_spin_eo2_up  =  dev_spin_eo2_up - epsbar*spinin_up  =  (1+imubar)*spinin_dn - epsbar*spinin_up
  cublasSaxpy (N_floats, -g_epsbar, (float *) spinin_dn, 1, (float *) dev_spin_eo2_dn, 1);
  												// dev_spin_eo2_dn  =  dev_spin_eo2_dn - epsbar*spinin_dn  =  (1-imubar)*spinin_up - epsbar*spinin_dn
  
  													// this is ((M_oo) - (M_oe)(Mee^-1)(M_eo)) * (spinin_dn, spinin_up):
  cublasSaxpy (N_floats, -1.0, (float *) dev_spin_eo1_up, 1, (float *) dev_spin_eo2_up, 1);
  							// dev_spin_eo2_up  =  dev_spin_eo2_up  -  dev_spin_eo1_up  =                      (1+imubar) * spinin_dn  -                    epsbar * spinin_up
  							//                                                             - (M_oe)*nrm*(1-imubar)*(M_eo) * spinin_dn  -  (M_oe)*nrm*epsbar*(M_eo) * spinin_up
  cublasSaxpy (N_floats, -1.0, (float *) dev_spin_eo1_dn, 1, (float *) dev_spin_eo2_dn, 1);
  							// dev_spin_eo2_dn  =  dev_spin_eo2_dn  -  dev_spin_eo1_dn  =                      (1-imubar) * spinin_up  -                    epsbar * spinin_dn
  							//                                                             - (M_oe)*nrm*(1+imubar)*(M_eo) * spinin_up  -  (M_oe)*nrm*epsbar*(M_eo) * spinin_dn
  */
  
  dev_nd_linalg1<<<gridsize2, blocksize2>>>(dev_spin_eo1_up, dev_spin_eo1_dn, spinin_up, spinin_dn, dev_spin_eo2_up, dev_spin_eo2_dn, g_epsbar ); 
  
  
  #ifdef RELATIVISTIC_BASIS
    dev_gamma5_rel<<<gridsize3, blocksize3>>>(dev_spin_eo2_up, dev_spin_eo2_up);					// dev_spin_eo2_up = gamma5 * dev_spin_eo2_up 
    dev_gamma5_rel<<<gridsize3, blocksize3>>>(dev_spin_eo2_dn, dev_spin_eo2_dn);					// dev_spin_eo2_dn = gamma5 * dev_spin_eo2_dn  
  #else
    dev_gamma5<<<gridsize3, blocksize3>>>(dev_spin_eo2_up, dev_spin_eo2_up);					// dev_spin_eo2_up = gamma5 * dev_spin_eo2_up 
    dev_gamma5<<<gridsize3, blocksize3>>>(dev_spin_eo2_dn, dev_spin_eo2_dn);					// dev_spin_eo2_dn = gamma5 * dev_spin_eo2_dn
  #endif
  
  
  
  ////////////////////
  // (a,b) -> (b,a) //
  ////////////////////
  dev_copy_spinor_field<<<gridsize4, blocksize4>>>(dev_spin_eo2_dn, dev_spin_eo3_up);			// dev_spin_eo3_up = dev_spin_eo2_dn
  dev_copy_spinor_field<<<gridsize4, blocksize4>>>(dev_spin_eo2_up, dev_spin_eo3_dn);			// dev_spin_eo3_dn = dev_spin_eo2_up
  
  
  
  #ifdef USETEXTURE
    bind_texture_spin(dev_spin_eo3_up,1);
    bind_texture_spin_dn(dev_spin_eo3_dn,1);   
  #endif
    dev_Hopping_Matrix_updn<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo3_up, dev_spin_eo3_dn, dev_spin_eo1_up, dev_spin_eo1_dn, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0,0,N_sites);	// dev_spin_eo1_up = (M_eo) * dev_spin_eo3_up
  #ifdef USETEXTURE
    unbind_texture_spin(1);
    unbind_texture_spin_dn(1);    
  #endif
  
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo1_up, dev_spin_eo2_up, -1.0);			// dev_spin_eo2_up  =  (1 - imubar) * dev_spin_eo1_up  =  (1 - imubar)*(M_eo) * dev_spin_eo3_up
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo1_dn, dev_spin_eo2_dn, +1.0);			// dev_spin_eo2_dn  =  (1 + imubar) * dev_spin_eo1_dn  =  (1 + imubar)*(M_eo) * dev_spin_eo3_dn

  /*
  cublasSaxpy (N_floats, g_epsbar, (float *) dev_spin_eo1_dn, 1, (float *) dev_spin_eo2_up, 1);
  										// dev_spin_eo2_up  =  dev_spin_eo2_up  +  epsbar * dev_spin_eo1_dn  =  (1 - imubar)*(M_eo) * dev_spin_eo3_up  +  epsbar * (M_eo) * dev_spin_eo3_dn
  cublasSaxpy (N_floats, g_epsbar, (float *) dev_spin_eo1_up, 1, (float *) dev_spin_eo2_dn, 1);
  										// dev_spin_eo2_dn  =  dev_spin_eo2_dn  +  epsbar * dev_spin_eo1_up  =  (1 + imubar)*(M_eo) * dev_spin_eo3_dn  +  epsbar * (M_eo) * dev_spin_eo3_up
  
  cublasSscal (N_floats, nrm, (float *) dev_spin_eo2_up, 1);			// dev_spin_eo2_up  =  nrm * dev_spin_eo2_up  =  nrm*(1-imubar)*(M_eo) * dev_spin_eo3_up  +  nrm*epsbar*(M_eo) * dev_spin_eo3_dn
  
  cublasSscal (N_floats, nrm, (float *) dev_spin_eo2_dn, 1);			// dev_spin_eo2_dn  =  nrm * dev_spin_eo2_dn  =  nrm*(1+imubar)*(M_eo) * dev_spin_eo3_dn  +  nrm*epsbar*(M_eo) * dev_spin_eo3_up
  */
  
  dev_nd_linalg2<<<gridsize2, blocksize2>>>(dev_spin_eo1_up, dev_spin_eo1_dn, dev_spin_eo2_up, dev_spin_eo2_dn, g_epsbar, nrm); 

  
  #ifdef USETEXTURE
    bind_texture_spin(dev_spin_eo2_up,1);	
    bind_texture_spin_dn(dev_spin_eo2_dn,1);    
  #endif
    dev_Hopping_Matrix_updn<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo2_up, dev_spin_eo2_dn, dev_spin_eo1_up,dev_spin_eo1_dn, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1,0,N_sites);
  #ifdef USETEXTURE
    unbind_texture_spin(1);
    unbind_texture_spin_dn(1);
  #endif
  

  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo3_up, dev_spin_eo2_up, +1.0);				// dev_spin_eo2_up = (1 + imubar) * dev_spin_eo3_up
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo3_dn, dev_spin_eo2_dn, -1.0);				// dev_spin_eo2_dn = (1 - imubar) * dev_spin_eo3_dn
  
  /*											// remember: this is (M_oo) * (dev_spin_eo3_up, dev_spin_eo3_dn):
  cublasSaxpy (N_floats, -g_epsbar, (float *) dev_spin_eo3_dn, 1, (float *) dev_spin_eo2_up, 1);
  												// dev_spin_eo2_up  =  dev_spin_eo2_up - epsbar*dev_spin_eo3_dn  =  (1+imubar)*dev_spin_eo3_up - epsbar*dev_spin_eo3_dn
  cublasSaxpy (N_floats, -g_epsbar, (float *) dev_spin_eo3_up, 1, (float *) dev_spin_eo2_dn, 1);
  												// dev_spin_eo2_dn  =  dev_spin_eo2_dn - epsbar*dev_spin_eo3_up  =  (1-imubar)*dev_spin_eo3_dn - epsbar*dev_spin_eo3_up
  
  

  											// this is ( (M_oo) - (M_oe) (Mee^-1) (M_eo) ) * (dev_spin_eo3_up, dev_spin_eo3_dn)
  cublasSaxpy (N_floats, -1.0, (float *) dev_spin_eo1_up, 1, (float *) dev_spin_eo2_up, 1);
  							// dev_spin_eo2_up  =  dev_spin_eo2_up  -  dev_spin_eo1_up  =                      (1+imubar) * dev_spin_eo3_up  -                    epsbar * dev_spin_eo3_dn
  							//                                                             - (M_oe)*nrm*(1-imubar)*(M_eo) * dev_spin_eo3_up  -  (M_oe)*nrm*epsbar*(M_eo) * dev_spin_eo3_dn
  cublasSaxpy (N_floats, -1.0, (float *) dev_spin_eo1_dn, 1, (float *) dev_spin_eo2_dn, 1);
  							// dev_spin_eo2_dn  =  dev_spin_eo2_dn  -  dev_spin_eo1_dn  =                      (1-imubar) * dev_spin_eo3_dn  -                    epsbar * dev_spin_eo3_up
  							//                                                             - (M_oe)*nrm*(1+imubar)*(M_eo) * dev_spin_eo3_dn  -  (M_oe)*nrm*epsbar*(M_eo) * dev_spin_eo3_up
  */
  dev_nd_linalg1<<<gridsize2, blocksize2>>>(dev_spin_eo1_up, dev_spin_eo1_dn, dev_spin_eo3_dn, dev_spin_eo3_up, dev_spin_eo2_up, dev_spin_eo2_dn, g_epsbar ); 

  
  #ifdef RELATIVISTIC_BASIS
    dev_gamma5_rel<<<gridsize3, blocksize3>>>(dev_spin_eo2_up, dev_spin_eo2_up);					// dev_spin_eo2_up = gamma5 * dev_spin_eo2_up 
    dev_gamma5_rel<<<gridsize3, blocksize3>>>(dev_spin_eo2_dn, dev_spin_eo2_dn);					// dev_spin_eo2_dn = gamma5 * dev_spin_eo2_dn  
  #else
    dev_gamma5<<<gridsize3, blocksize3>>>(dev_spin_eo2_up, dev_spin_eo2_up);					// dev_spin_eo2_up = gamma5 * dev_spin_eo2_up 
    dev_gamma5<<<gridsize3, blocksize3>>>(dev_spin_eo2_dn, dev_spin_eo2_dn);					// dev_spin_eo2_dn = gamma5 * dev_spin_eo2_dn
  #endif
  
  cublasSscal (N_floats, phmc_invmaxev*phmc_invmaxev, (float *) dev_spin_eo2_up, 1);
  cublasSscal (N_floats, phmc_invmaxev*phmc_invmaxev, (float *) dev_spin_eo2_dn, 1);
  
  return;
  
}//dev_Qtm_pm_ndpsi_updn()
















#ifdef MPI

///////////////////////////
// MATRIX MULTIPLICATION //
///////////////////////////

// the GPU implementation of  Qtm_pm_ndpsi(...)  tm_operators_nd.c
//	Flo's equivalent function for the standard and non-nd case is  dev_Qtm_pm_psi

void dev_Qtm_pm_ndpsi_mpi (dev_spinor * spinout_up, dev_spinor * spinout_dn,
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
  //   MATCHING with  Qtm_pm_ndpsi  //
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
    	  dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, spinin_dn, dev_spin_eo1_up, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0,0,N_sites);	// dev_spin_eo1_up = (M_eo) * spinin_dn
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
    	  dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, spinin_up, dev_spin_eo1_dn, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0,0,N_sites);	// dev_spin_eo1_dn = (M_eo) * spinin_up
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
    	  dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo2_up, dev_spin_eo1_up, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1,0,N_sites);
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
    	  dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo2_dn, dev_spin_eo1_dn, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1,0,N_sites);
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
  #ifdef RELATIVISTIC_BASIS
    dev_gamma5_rel<<<gridsize3, blocksize3>>>(dev_spin_eo2_up, dev_spin_eo2_up);					// dev_spin_eo2_up = gamma5 * dev_spin_eo2_up 
    dev_gamma5_rel<<<gridsize3, blocksize3>>>(dev_spin_eo2_dn, dev_spin_eo2_dn);					// dev_spin_eo2_dn = gamma5 * dev_spin_eo2_dn  
  #else
    dev_gamma5<<<gridsize3, blocksize3>>>(dev_spin_eo2_up, dev_spin_eo2_up);					// dev_spin_eo2_up = gamma5 * dev_spin_eo2_up 
    dev_gamma5<<<gridsize3, blocksize3>>>(dev_spin_eo2_dn, dev_spin_eo2_dn);					// dev_spin_eo2_dn = gamma5 * dev_spin_eo2_dn
  #endif	
  
  
  
  
  
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
    	  dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo3_up, dev_spin_eo1_up, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0,0,N_sites);	// dev_spin_eo1_up = (M_eo) * dev_spin_eo3_up
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
    	  dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo3_dn, dev_spin_eo1_dn, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0,0,N_sites);	// dev_spin_eo1_dn = (M_eo) * dev_spin_eo3_dn
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
  cublasSaxpy (N_floats, g_epsbar, (float *) dev_spin_eo1_dn, 1, (float *) dev_spin_eo2_up, 1);										// dev_spin_eo2_up  =  dev_spin_eo2_up  +  epsbar * dev_spin_eo1_dn  =  (1 - imubar)*(M_eo) * dev_spin_eo3_up  +  epsbar * (M_eo) * dev_spin_eo3_dn
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
    	  dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo2_up, dev_spin_eo1_up, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1,0,N_sites);
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
    	  dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo2_dn, dev_spin_eo1_dn, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1,0,N_sites);
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
  cublasSaxpy (N_floats, -g_epsbar, (float *) dev_spin_eo3_dn, 1, (float *) dev_spin_eo2_up, 1);												// dev_spin_eo2_up  =  dev_spin_eo2_up - epsbar*dev_spin_eo3_dn  =  (1+imubar)*dev_spin_eo3_up - epsbar*dev_spin_eo3_dn
  cublasSaxpy (N_floats, -g_epsbar, (float *) dev_spin_eo3_up, 1, (float *) dev_spin_eo2_dn, 1);
  												// dev_spin_eo2_dn  =  dev_spin_eo2_dn - epsbar*dev_spin_eo3_up  =  (1-imubar)*dev_spin_eo3_dn - epsbar*dev_spin_eo3_up
  
  ///////////////////////////////////////
  // (M_oo)  -  (M_oe) (Mee^-1) (M_eo) //
  ///////////////////////////////////////
  
  // lineare algebra										// this is ( (M_oo) - (M_oe) (Mee^-1) (M_eo) ) * (dev_spin_eo3_up, dev_spin_eo3_dn)
  cublasSaxpy (N_floats, -1.0, (float *) dev_spin_eo1_up, 1, (float *) dev_spin_eo2_up, 1);
  							// dev_spin_eo2_up  =  dev_spin_eo2_up  -  dev_spin_eo1_up  =                      (1+imubar) * dev_spin_eo3_up  -                    epsbar * dev_spin_eo3_dn						//                                                             - (M_oe)*nrm*(1-imubar)*(M_eo) * dev_spin_eo3_up  -  (M_oe)*nrm*epsbar*(M_eo) * dev_spin_eo3_dn
  cublasSaxpy (N_floats, -1.0, (float *) dev_spin_eo1_dn, 1, (float *) dev_spin_eo2_dn, 1);
  							// dev_spin_eo2_dn  =  dev_spin_eo2_dn  -  dev_spin_eo1_dn  =                      (1-imubar) * dev_spin_eo3_dn  -                    epsbar * dev_spin_eo3_up
  							//                                                             - (M_oe)*nrm*(1+imubar)*(M_eo) * dev_spin_eo3_dn  -  (M_oe)*nrm*epsbar*(M_eo) * dev_spin_eo3_up
 
  ////////////
  // gamma5 //
  ////////////
  
  // gamma5
  #ifdef RELATIVISTIC_BASIS
    dev_gamma5_rel<<<gridsize3, blocksize3>>>(dev_spin_eo2_up, dev_spin_eo2_up);					// dev_spin_eo2_up = gamma5 * dev_spin_eo2_up 
    dev_gamma5_rel<<<gridsize3, blocksize3>>>(dev_spin_eo2_dn, dev_spin_eo2_dn);					// dev_spin_eo2_dn = gamma5 * dev_spin_eo2_dn  
  #else
    dev_gamma5<<<gridsize3, blocksize3>>>(dev_spin_eo2_up, dev_spin_eo2_up);					// dev_spin_eo2_up = gamma5 * dev_spin_eo2_up 
    dev_gamma5<<<gridsize3, blocksize3>>>(dev_spin_eo2_dn, dev_spin_eo2_dn);					// dev_spin_eo2_dn = gamma5 * dev_spin_eo2_dn
  #endif

  cublasSscal (N_floats, phmc_invmaxev*phmc_invmaxev, (float *) dev_spin_eo2_up, 1);
  cublasSscal (N_floats, phmc_invmaxev*phmc_invmaxev, (float *) dev_spin_eo2_dn, 1);   

  ////////////
  // output //	output is already done by setting  dev_spin_eo2_up/dn = spinout_up/dn
  ////////////

  
  return;
  
}//dev_Qtm_pm_ndpsi_mpi()


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
  // dev_Hopping_Matrix	          = 1608								//
  // dev_mul_one_pm_imubar_gamma5 = 120									//
  // dev_gamma5                   = 12									//
  //													//
  // cublasSaxpy                  = 24*2 = 48								//
  // cublasSscal                  = 24*1 = 24								//
  //													//
  //													//
  // (FLOPS/matrix application)  =  2 * (4*1608 + 4*120 + 6*48 + 2*24 + 2*12)  =  2 * 7224  =  14448	//
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
  //double effectiveFlopsPerApp = 21296.0;	// per lattice site
  double effectiveFlopsPerApp = 14448.0; // hopping = 1608
  
  #ifndef MPI
    double effectiveDeviceFlops;
    double effectiveFlops;
  #else
    double effectiveDeviceFlops;
    double allEffectiveDeviceFlops;
    double effectiveFlops;
  #endif
  
  // CUDA errors
  cudaError_t cudaerr;


  // formal parameters
  int staticsource = 0;		// 1: applies matrix every time on the same source
  				// 0: applies matrix consecutively ...
  

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

  int blocksize;		// auxiliary
  
  blocksize = BLOCK;
  int griddim2;					// passed:	dev_Hopping_Matrix
  dim3 blockdim2 (0,0,0);
  if ( (VOLUME/2) % blocksize == 0 ) {
    blockdim2.x = blocksize;
    griddim2  = VOLUME/2/blocksize;
  }
  else {
    blockdim2.x = blocksize;
    griddim2  = (int) ((VOLUME/2/blocksize) + 1);
  }
  
  blocksize = BLOCK2;
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

  
  		// debug
  		#ifndef MPI
  		  printf("Applying the eo-preconditioned matrix %i times.\n", N);
  		#else
  		  if (g_proc_id == 0) printf("Applying the eo-preconditioned matrix %i times.\n", N);
  		#endif
  
  
  to_device(B_up, Q_up, h2d_spin_up, dev_spinsize_int);
  to_device(B_dn, Q_dn, h2d_spin_dn, dev_spinsize_int);
  
  
  // timer
    startBenchmark = gettime();

  

  for (i = 0; i < N; i++) {
  
  
    #ifndef MPI
    	dev_Qtm_pm_ndpsi(A_up, A_dn,					// A = (matrix)*B
    	                        B_up, B_dn, 
    	                        griddim2, blockdim2.x,
    	                        griddim3, blockdim3,
    	                        griddim4, blockdim4,
    	                        griddim5, blockdim5);
    #else
    	#ifndef ASYNC
    	  dev_Qtm_pm_ndpsi_mpi(A_up, A_dn,				// A = (matrix)*B
    	                              B_up, B_dn, 
    	                              griddim2, blockdim2,
    	                              griddim3, blockdim3,
    	                              griddim4, blockdim4,
    	                              griddim5, blockdim5);
    	#else
    	  dev_Qtm_pm_ndpsi_mpi_ASYNC(A_up, A_dn,				// A = (matrix)*B
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

  }
  
  printf("%s\n", cudaGetErrorString(cudaGetLastError())); 
  printf("Done\n"); 
  cudaThreadSynchronize();
  
  
  		// debug	// CUDA
  		#ifdef CUDA_DEBUG
  		  CUDA_CHECK_NO_SUCCESS_MSG("CUDA error in matrix_muliplication32(). Applying the matrix on GPU failed.");
  		#endif

  // timer
    stopBenchmark = gettime();

  
  
  #ifndef MPI
  
  	timeElapsed = stopBenchmark - startBenchmark;
  	effectiveDeviceFlops = N * VOLUME/2 * effectiveFlopsPerApp;
  	effectiveFlops       = N * VOLUME/2 * effectiveFlopsPerApp / timeElapsed / 1.0e9;
  	printf("EFFECTIVE:\n");
  	printf("\ttime:        %.4e sec\n", timeElapsed);
  	printf("\tflop's:      %.4e flops\n", effectiveDeviceFlops);
  	printf("\tperformance: %.4e Gflop/s\n\n", effectiveFlops);
  	
  #else
  	
  	singleTimeElapsed = stopBenchmark - startBenchmark;
  	MPI_Allreduce(&singleTimeElapsed, &maxTimeElapsed, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  	effectiveDeviceFlops = N * VOLUME/2 * effectiveFlopsPerApp;
  	MPI_Allreduce(&effectiveDeviceFlops, &allEffectiveDeviceFlops, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  	effectiveFlops       = allEffectiveDeviceFlops / maxTimeElapsed / 1.0e9;
  	
  	
  	if (g_proc_id == 0) {

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

int dev_cg_eo_nd (dev_su3_2v * gf,
              dev_spinor * P_up, dev_spinor * P_dn,
              dev_spinor * Q_up, dev_spinor * Q_dn,
	      float shift,
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
  
  // algorithm control parameters
  // int N_recalc_res = 10;		// recalculate residue r(k+1) = b - A*x(k+1) each N_recalc_res iteration
  int N_recalc_res = 1000;

  
  
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
  
  blocksize = BLOCK;
  int blockdim2, griddim2;					// passed:	dev_Hopping_Matrix
  if ( (VOLUME/2) % blocksize == 0 ) {
    blockdim2 = blocksize;
    griddim2  = VOLUME/2/blocksize;
  }
  else {
    blockdim2 = blocksize;
    griddim2  = (int) ((VOLUME/2/blocksize) + 1);
  }
  
  blocksize = BLOCK2;
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

 #ifdef RELATIVISTIC_BASIS 
   //transform to relativistic gamma basis
   to_relativistic_basis<<<griddim4, blockdim4>>> (Q_up);
   to_relativistic_basis<<<griddim4, blockdim4>>> (Q_dn);
   to_relativistic_basis<<<griddim4, blockdim4>>> (P_up);
   to_relativistic_basis<<<griddim4, blockdim4>>> (P_dn);   
   
   if ((cudaerr=cudaGetLastError())!=cudaSuccess) {
     if (g_proc_id == 0) printf("%s\n", cudaGetErrorString(cudaGetLastError()));
   }
   else{
     #ifndef LOWOUTPUT 
     if (g_proc_id == 0) printf("Switched to relativistic basis\n");
     #endif
   }
 #endif

 
  
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
  
  
  		// debug	// CUBLAS helper function
  		#ifdef CUDA_DEBUG
  		  CUBLAS_HELPER_CHECK(cublasInit(), "CUBLAS error in cublasInit(). Couldn't initialize CUBLAS.", "CUBLAS initialized.");
  		#else
		  #ifdef CUDA_45
		    cublasHandle_t handle;
		    cublasCreate(&handle);
		  #else
		    cublasInit();
		  #endif 
  		#endif

  if(g_debug_level > 3) printf("cublasstatus = %f\n", cublasstatus);
  
  
  
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
  		  CUDA_KERNEL_CHECK("Kernel error in dev_cg_eo_nd(). Initializing spinor fields on device failed.", "Spinor fields initialized on device.");
  		#endif
  
  
  
  
  // rr = (r_up)^2 + (r_dn)^2
  #ifndef MPI
    rr_up = cublasSdot(N_floats_int, (float *) r_up, 1, (float *) r_up, 1);
    rr_dn = cublasSdot(N_floats_int, (float *) r_dn, 1, (float *) r_dn, 1);
  #else
    rr_up = cublasSdot_wrapper(N_floats_int, (float *) r_up, 1, (float *) r_up, 1);
    rr_dn = cublasSdot_wrapper(N_floats_int, (float *) r_dn, 1, (float *) r_dn, 1);
  #endif
  rr    = rr_up + rr_dn;
  



  r0r0   = rr;	// for relative precision 
  rr_old = rr;	// for the first iteration
  
  
  
  
  //////////
  // LOOP //
  //////////
  
               
                #ifndef LOWOUTPUT
    		if (g_proc_id == 0) printf("\nEntering inner loop.\n");
	        #endif
  
  		// debug	// CUBLAS core function
		#ifdef CUDA_DEBUG
		  // CUBLAS_CORE_CHECK("CUBLAS error in cg_eo_nd(). Calculating initial residue failed.", "Initial inner residue calculated.");
		  CUBLAS_CORE_CHECK_NO_SUCCESS_MSG("CUBLAS error in dev_cg_eo_nd(). Calculating initial residue failed.");
		#endif
  
  		if (g_proc_id == 0) printf("Initial inner residue: %.6e\n", r0r0);
  
  
  for (j = 0; j < max_iter; j++) {

     
      // A*d(k)
      #ifndef MPI
      		dev_Qtm_pm_ndpsi(Ad_up, Ad_dn,										// normally:  dev_Qtm_pm_ndpsi()
      		                         d_up,  d_dn,										// debugging: matrix_debug1(), matrix_multiplication_test()
      		                        griddim2, blockdim2,
      		                        griddim3, blockdim3,
      		                        griddim4, blockdim4,
      		                        griddim5, blockdim5);

    
      #else
      	#ifndef ASYNC
        	dev_Qtm_pm_ndpsi_mpi(Ad_up, Ad_dn,									// normally:  dev_Qtm_pm_ndpsi_mpi()
        	                             d_up,  d_dn, 									// debugging: matrix_mpi_debug1/2/3/4()
        	                            griddim2, blockdim2,
        	                            griddim3, blockdim3,
        	                            griddim4, blockdim4,
        	                            griddim5, blockdim5);
      	#else															// tries to overlap computation and communication
        	dev_Qtm_pm_ndpsi_mpi_ASYNC(Ad_up, Ad_dn, 
        	                                   d_up,  d_dn, 
        	                                  griddim2, blockdim2,
        	                                  griddim3, blockdim3,
        	                                  griddim4, blockdim4,
        	                                  griddim5, blockdim5);
      	#endif
      #endif	// MPI
      
      
	if(shift != 0.0f){
           //add constant shift if nonzero
           // CUBLAS:
          cublasSaxpy (N_floats_int, shift, (float *) d_up, 1, (float *) Ad_up, 1);
          cublasSaxpy (N_floats_int, shift, (float *) d_dn, 1, (float *) Ad_dn, 1);
        }      
        if((cudaerr=cudaGetLastError()) != cudaSuccess){
           printf("%s\n", cudaGetErrorString(cudaerr));
           exit(200);
        }
      
  		// debug	// CUDA		// also other stuff ?!
  		#ifdef CUDA_DEBUG
  		  CUDA_CHECK("CUDA error in matrix_muliplication32(). Applying the matrix on GPU failed.", "The matrix was applied on GPU.");
  		  //CUDA_CHECK_NO_SUCCESS_MSG("CUDA error in matrix_muliplication32(). Applying the matrix on GPU failed.");
  		#endif

    // alpha = r(k)*r(k) / d(k)*A*d(k)
    #ifndef MPI
      dAd_up = cublasSdot(N_floats_int, (float *) d_up, 1, (float *) Ad_up, 1);
      dAd_dn = cublasSdot(N_floats_int, (float *) d_dn, 1, (float *) Ad_dn, 1);
    #else
      dAd_up = cublasSdot_wrapper(N_floats_int, (float *) d_up, 1, (float *) Ad_up, 1);
      dAd_dn = cublasSdot_wrapper(N_floats_int, (float *) d_dn, 1, (float *) Ad_dn, 1);
    #endif
    dAd    = dAd_up + dAd_dn;
    
    		// debug	// is NaN ?
    		if isnan(dAd) {
    		  printf("Error in dev_cg_eo_nd(). dAd is NaN.\n");
    		  exit(-1);
    		}
    
    alpha  = rr_old / dAd;	// rr_old is taken from the last iteration respectively
    
    // x(k+1) = x(k) + alpha*d(k)
    cublasSaxpy(N_floats_int, alpha, (float *) d_up, 1, (float *) x_up, 1);
    cublasSaxpy(N_floats_int, alpha, (float *) d_dn, 1, (float *) x_dn, 1);
    
    // r(k+1)
    if ( (j+1) % N_recalc_res != 0 ) {	// r(k+1) = r(k) - alpha*A*d(k)
      cublasSaxpy(N_floats_int, -1.0*alpha, (float *) Ad_up, 1, (float *) r_up, 1);
      cublasSaxpy(N_floats_int, -1.0*alpha, (float *) Ad_dn, 1, (float *) r_dn, 1);
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

      	#ifndef MPI
        	dev_Qtm_pm_ndpsi(Ax_up, Ax_dn,
        	                         x_up,  x_dn, 
        	                        griddim2, blockdim2,
        	                        griddim3, blockdim3,
        	                        griddim4, blockdim4,
        	                        griddim5, blockdim5);
        #else
        	#ifndef ASYNC
        	  dev_Qtm_pm_ndpsi_mpi(Ax_up, Ax_dn,									// normally:  dev_Qtm_pm_ndpsi_mpi()
        	                               x_up,  x_dn, 									// debugging: matrix_mpi_debug1/2/3/4()
        	                              griddim2, blockdim2,
        	                              griddim3, blockdim3,
        	                              griddim4, blockdim4,
        	                              griddim5, blockdim5);
        	#else
        	  dev_Qtm_pm_ndpsi_mpi_ASYNC(Ax_up, Ax_dn,									// normally:  dev_Qtm_pm_ndpsi_mpi()
        	                                     x_up,  x_dn, 									// debugging: matrix_mpi_debug1/2/3/4()
        	                                    griddim2, blockdim2,
        	                                    griddim3, blockdim3,
        	                                    griddim4, blockdim4,
        	                                    griddim5, blockdim5);
        	#endif
        #endif	// MPI
	if(shift != 0.0f){
           //add constant shift if nonzero
           // CUBLAS:
          cublasSaxpy (N_floats_int, shift, (float *) x_up, 1, (float *) Ad_up, 1);
          cublasSaxpy (N_floats_int, shift, (float *) x_dn, 1, (float *) Ad_dn, 1);
        }    

      // r(k+1) = b - A*x(k+1)
      cublasScopy(N_floats_int, (float *) Q_up, 1, (float *) r_up, 1);		// r_up = Q_up
      cublasScopy(N_floats_int, (float *) Q_dn, 1, (float *) r_dn, 1);		// r_dn = Q_dn
      cublasSaxpy(N_floats_int, -1.0, (float *) Ax_up, 1, (float *) r_up, 1);	// r_up = Q_up - Ax_up
      cublasSaxpy(N_floats_int, -1.0, (float *) Ax_dn, 1, (float *) r_dn, 1);	// r_dn = Q_dn - Ax_dn
    
    } // recalculate residue
    

    // r(k+1)*r(k+1)
    #ifndef MPI
      rr_up  = cublasSdot(N_floats_int, (float *) r_up, 1, (float *) r_up, 1);
      rr_dn  = cublasSdot(N_floats_int, (float *) r_dn, 1, (float *) r_dn, 1);
    #else
      rr_up  = cublasSdot_wrapper(N_floats_int, (float *) r_up, 1, (float *) r_up, 1);
      rr_dn  = cublasSdot_wrapper(N_floats_int, (float *) r_dn, 1, (float *) r_dn, 1);
    #endif
    rr     = rr_up + rr_dn;
    
		// debug	// CUBLAS core function
		#ifdef CUDA_DEBUG
		  CUBLAS_CORE_CHECK_NO_SUCCESS_MSG("CUBLAS error in dev_cg_eo_nd(). CUBLAS function failed.");
		#endif
    
    
               #ifndef LOWOUTPUT
    		  if (g_proc_id == 0) printf("inner iteration j = %i: rr = %.6e\n", j, rr);
	       #endif
		 
	
    		// debug	// is NaN ?
    		if isnan(rr) {
    		  printf("Error in dev_cg_eo_nd(). Inner residue is NaN.\n");
    		  exit(-1);
    		}
    
    
    // aborting ?? // check wether precision is reached ...
    if ( ((check_abs)&&(rr <= eps_abs)) || ((check_rel)&&(rr <= eps_rel*r0r0)) ) {
    
      #ifdef MPI
        if (g_proc_id == 0) {
      #endif
      
      if ((check_rel)&&(rr <= eps_rel*r0r0)) {
      	printf("Reached relative inner solver precision of eps_rel = %.2e\n", eps_rel);
      }
      if ((check_abs)&&(rr <= eps_abs)) {
      	printf("Reached absolute inner solver precision of eps_abs = %.2e\n", eps_abs);
      }
      
      printf("Final inner residue: %.6e\n", rr);
      
      #ifdef MPI
        }
      #endif
      
      #ifdef RELATIVISTIC_BASIS 
         //transform back to tmlqcd gamma basis
         to_tmlqcd_basis<<<griddim4, blockdim4>>> (x_up);
	 to_tmlqcd_basis<<<griddim4, blockdim4>>> (x_dn);
      #endif
      
      #ifdef CUDA_45  
	cublasDestroy(handle);
      #else
	cublasShutdown();
      #endif 
      
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
    		  CUBLAS_CORE_CHECK_NO_SUCCESS_MSG("CUBLAS error in dev_cg_eo_nd(). Error in CUBLAS function.");
    		#endif
    		
  
  }//LOOP
  
  
  #ifndef LOWOUTPUT
  if (g_proc_id == 0) printf("Finished inner loop because of maximal number of inner iterations.\n");
  #endif
  if (g_proc_id == 0) printf("Final inner residue: %.6e\n", rr);



      #ifdef RELATIVISTIC_BASIS 
         //transform back to tmlqcd gamma basis
         to_tmlqcd_basis<<<griddim4, blockdim4>>> (x_up);
	 to_tmlqcd_basis<<<griddim4, blockdim4>>> (x_dn);
      #endif
    
    #ifdef CUDA_45  
      cublasDestroy(handle);
    #else
      cublasShutdown();
    #endif 
  return(j+1);
  
}//dev_cg_eo_nd()





void test_double_nd_operator(spinor* const Q_up, spinor* const Q_dn, const int N){
   
   size_t dev_spinsize_d = 12*VOLUME/2 * sizeof(dev_spinor_d); // double4 even-odd !   
   int gridsize;
     //this is the partitioning for the HoppingMatrix kernel
     int blockdim3 = BLOCKD;
     if( VOLUME/2 % blockdim3 == 0){
       gridsize = (int) VOLUME/2/blockdim3;
     }
     else{
       gridsize = (int) VOLUME/2/blockdim3 + 1;
     }
     int griddim3 = gridsize;
   
     //this is the partitioning for dev_mul_one_pm...
     int blockdim4 = BLOCK2D;
     if( VOLUME/2 % blockdim4 == 0){
       gridsize = (int) VOLUME/2/blockdim4;
     }
     else{
       gridsize = (int) VOLUME/2/blockdim4 + 1;
     }
     int griddim4 = gridsize;    
     
      dev_spinor_d * x_up_d = dev_spin0_d;
      dev_spinor_d * x_dn_d = dev_spin1_d;
      dev_spinor_d * Ax_up_d = dev_spin2_d;
      dev_spinor_d * Ax_dn_d = dev_spin3_d;
      dev_spinor_d * Q_up_d = dev_spin_eo1_d;
      dev_spinor_d * Q_dn_d = dev_spin_eo2_d;     
      
  spinor ** solver_field_up = NULL;
  spinor ** solver_field_dn = NULL;  
  const int nr_sf = 3;
  init_solver_field(&solver_field_up, VOLUMEPLUSRAND/2, nr_sf);  
  init_solver_field(&solver_field_dn, VOLUMEPLUSRAND/2, nr_sf); 

  //apply cpu matrix

  
   Qtm_pm_ndpsi(solver_field_up[0], solver_field_dn[0], Q_up, Q_dn);
  
  /*
  Hopping_Matrix(EO,solver_field_up[1], Q_up);
  Hopping_Matrix(OE, solver_field_up[0] , solver_field_up[1]); 
  */
  //apply gpu matrix
  order_spin_gpu(Q_up, h2d_spin_d);
  cudaMemcpy(x_up_d, h2d_spin_d, dev_spinsize_d, cudaMemcpyHostToDevice);      
  order_spin_gpu(Q_dn, h2d_spin_d);
  cudaMemcpy(x_dn_d, h2d_spin_d, dev_spinsize_d, cudaMemcpyHostToDevice);	  
  // r_up/dn = Q-A*x_up/dn
  
  dev_Qtm_pm_ndpsi_d(Ax_up_d, Ax_dn_d,  
		      x_up_d, x_dn_d, 
		      griddim3, blockdim3, griddim4, blockdim4,
		      griddim4, blockdim4, griddim4, blockdim4);  
  /*
  dev_Hopp_d(x_dn_d, x_up_d,  
		      griddim3, blockdim3, griddim4, blockdim4,0); 
  dev_Hopp_d(Ax_up_d, x_dn_d,  
		      griddim3, blockdim3, griddim4, blockdim4,1);   
  */
  cudaMemcpy(h2d_spin_d, Ax_up_d, dev_spinsize_d, cudaMemcpyDeviceToHost);
  unorder_spin_gpu(h2d_spin_d, solver_field_up[1]);      

  cudaMemcpy(h2d_spin_d, Ax_dn_d, dev_spinsize_d, cudaMemcpyDeviceToHost);
  unorder_spin_gpu(h2d_spin_d, solver_field_dn[1]);     
  
  diff(solver_field_up[2], solver_field_up[1], solver_field_up[0],N);
  diff(solver_field_dn[2], solver_field_dn[1], solver_field_dn[0],N);  
  
  int at_max = -1;
  int at_min = -1;
  double max_dev = 0.0;
  double min_dev = 1.0;  
  double dev;
  _Complex double cdev;
  spinor * s = solver_field_up[2];
  for(int i=0; i<N; i++){
    
    cdev = conj(s->s0.c0) * s->s0.c0 +
         conj(s->s0.c1) * s->s0.c1 +
         conj(s->s0.c2) * s->s0.c2 +
         conj(s->s1.c0) * s->s1.c0 +
         conj(s->s1.c1) * s->s1.c1 +
         conj(s->s1.c2) * s->s1.c2 +
         conj(s->s2.c0) * s->s2.c0 +
         conj(s->s2.c1) * s->s2.c1 +
         conj(s->s2.c2) * s->s2.c2 +
         conj(s->s3.c0) * s->s3.c0 +
         conj(s->s3.c1) * s->s3.c1 +
         conj(s->s3.c2) * s->s3.c2;   
     dev = creal(cdev);
      
     if(dev > max_dev){
       max_dev = dev;
       at_max=i;
    }
    if(dev < min_dev){
       min_dev = dev;
       at_min=i;
    } 
    s++;
  }
  
  
  double rk_up = square_norm(solver_field_up[2], N, 1);
  double rk_dn = square_norm(solver_field_dn[2], N, 1);  
  double rk = rk_up + rk_dn;
  printf("Testing double matrix:\n");
  printf("cpu: Squared difference is   UP: %.8e\n", rk_up);
  printf("cpu: Squared difference is DOWN: %.8e\n", rk_dn);  
  printf("cpu: Squared difference per spinor component is: %.8e\n", rk/N/24.0);  
  printf("Max. difference at position %i: %.8e\n", at_max, max_dev);
  printf("Min. difference at position %i: %.8e\n", at_min, min_dev);  
  dev_complex_d h0,h1,h2,h3;
  cudaMemcpyFromSymbol(&h0, dev_k0_d, sizeof(dev_complex_d));
  printf("cpu: k0.re: %.16e\t k0.im: %.16e\n", creal(ka0), cimag(ka0));
  printf("gpu: k0.re: %.16e\t k0.im: %.16e\n", h0.re, h0.im);  
  printf("diff: k0.re: %.16e\t k0.im: %.16e\n", (creal(ka0)-h0.re), (cimag(ka0)-h0.im));
  
  printf("\n");
  cudaMemcpyFromSymbol(&h1, dev_k1_d, sizeof(dev_complex_d));  
  printf("cpu: k1.re: %.16e\t k1.im: %.16e\n", creal(ka1), cimag(ka1));
  printf("gpu: k1.re: %.16e\t k1.im: %.16e\n", h1.re, h1.im);  
  printf("diff: k1.re: %.16e\t k1.im: %.16e\n", (creal(ka1)-h1.re), (cimag(ka1)-h1.im));
  
  printf("\n");
  cudaMemcpyFromSymbol(&h2, dev_k2_d, sizeof(dev_complex_d));  
  printf("cpu: k2.re: %.16e\t k2.im: %.16e\n", creal(ka2), cimag(ka2));
  printf("gpu: k2.re: %.16e\t k2.im: %.16e\n", h2.re, h2.im);  
  printf("diff: k2.re: %.16e\t k2.im: %.16e\n", (creal(ka1)-h2.re), (cimag(ka1)-h2.im));
 
  printf("\n");
  cudaMemcpyFromSymbol(&h3, dev_k3_d, sizeof(dev_complex_d));  
  printf("cpu: k3.re: %.16e\t k3.im: %.16e\n", creal(ka3), cimag(ka3));
  printf("gpu: k3.re: %.16e\t k3.im: %.16e\n", h3.re, h3.im);  
  printf("diff: k3.re: %.16e\t k3.im: %.16e\n", (creal(ka1)-h3.re), (cimag(ka1)-h3.im));
   
  finalize_solver(solver_field_up, nr_sf); 
  finalize_solver(solver_field_dn, nr_sf);  
}





//////////////////
// OUTER SOLVER //
//////////////////

// iterative refinement, defect correction
// that function is to replace the call of  cg_her_nd()  in  invert_doublet_eo.c
// solves the odd part of the full eo and nd problem
//	more precisely we have to invert  Qhat(2x2)*Qhat(2x2)^dagger
//	multiplying by  Qhat(2x2)^dagger  is done in  invert_doublet_eo.c

extern "C" int mixedsolve_eo_nd (spinor * P_up, spinor * P_dn,
                                 spinor * Q_up, spinor * Q_dn, double shift,
                                 int max_iter, double eps_sq, int rel_prec) {
  
  if(rel_prec) {
    innersolver_precision_check_rel = 1;
    innersolver_precision_check_abs = 0;
  }
  else {
    innersolver_precision_check_rel = 0;
    innersolver_precision_check_abs = 1;
  }


  // basically  P_up/dn  and  Q_up/dn  could be used as auxiliary fields
  //	P_up/dn  is the output field (and can be used as initial guess)
  //	Q_up/dn  is not used later in the calling  invert_doublet_eo.c
  //		but will be used as feedback in r(k+1) = b - A*x(k+1)
  // with shift a positive shift can be given to the matrix dev_Qtm_pm_ndpsi
  #ifdef STUFF_DEBUG
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
  #ifdef ALGORITHM_BENCHMARK
    double effectiveflops;
    //double hoppingflops = 1608.0;
    double matrixflops  = 14448;
    #ifdef MPI
      double allflops;				// flops added for all processes
    #endif
  #endif
  
  // timing
  double startouter, stopouter;
  double startinner, stopinner;
  // double timeelapsed;
  double innerclocks;
  double totalinnerclocks = 0;
  double totalouterclocks = 0;
  
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
	 

  //////////////////
  // INITIALIZING //
  //////////////////
  

    init_mixedsolve_eo_nd(g_gauge_field);			// initializes and allocates all quantities for the mixed solver
  								// more precise:
  								//	puts the gauge field on device as "2 rows" or "8 floats" per SU(3)-matrix
  								//	allocates memory for all spinor fields
  								//	puts the nn- and eoidx-fields on device memory
  
  
  // the following initializations are moved from dev_cg_eo_nd():
  
  // Initialize some stuff
  dev_complex h0, h1, h2, h3; 

  
  h0.re  =  (float) creal(ka0);	h0.im  = -(float) cimag(ka0);	// ka{0-4} are defined in boundary.c
  h1.re  =  (float) creal(ka1);	h1.im  = -(float) cimag(ka1);	// what is the meaning?
  h2.re  =  (float) creal(ka2);	h2.im  = -(float) cimag(ka2);
  h3.re  =  (float) creal(ka3);	h3.im  = -(float) cimag(ka3);
  
//    dev_complex mh0, mh1, mh2, mh3;
//   mh0.re = -(float) creal(ka0);	mh0.im =  (float) cimag(ka0);
//   mh1.re = -(float) creal(ka1);	mh1.im =  (float) cimag(ka1);
//   mh2.re = -(float) creal(ka2);	mh2.im =  (float) cimag(ka2);
//   mh3.re = -(float) creal(ka3);	mh3.im =  (float) cimag(ka3);

  //update the gpu single gauge_field
  update_gpu_gf(g_gauge_field); 

  #ifdef GPU_DOUBLE
      
      //dev_spin0_d == x_up
      //dev_spin1_d == x_dn
      //dev_spin2_d == Ax_up
      //dev_spin3_d == Ax_dn
      //dev_spin_eo1_d == Q_up
      //dev_spin_eo2_d == Q_dn
      //dev_spin_eo3_up_d  == r_up == diff (Q-Ax) and used in dev_Qtm_pm_ndpsi_d
      //dev_spin_eo3_dn_d  == r_dn == diff (Q-Ax) and used in dev_Qtm_pm_ndpsi_d    
      
      dev_spinor_d * x_up_d = dev_spin0_d;
      dev_spinor_d * x_dn_d = dev_spin1_d;
      dev_spinor_d * Ax_up_d = dev_spin2_d;
      dev_spinor_d * Ax_dn_d = dev_spin3_d;
      dev_spinor_d * Q_up_d = dev_spin_eo1_d;
      dev_spinor_d * Q_dn_d = dev_spin_eo2_d;  
      dev_spinor_d * r_up_d = dev_spin_eo3_up_d;  
      dev_spinor_d * r_dn_d = dev_spin_eo3_dn_d;      
      
       size_t dev_spinsize_d = 12*VOLUME/2 * sizeof(dev_spinor_d); // double4 even-odd !   
   int gridsize;
     //this is the partitioning for the HoppingMatrix kernel
     int blockdim3 = BLOCKD;
     if( VOLUME/2 % blockdim3 == 0){
       gridsize = (int) VOLUME/2/blockdim3;
     }
     else{
       gridsize = (int) VOLUME/2/blockdim3 + 1;
     }
     int griddim3 = gridsize;
   
     //this is the partitioning for dev_mul_one_pm...
     int blockdim4 = BLOCK2D;
     if( VOLUME/2 % blockdim4 == 0){
       gridsize = (int) VOLUME/2/blockdim4;
     }
     else{
       gridsize = (int) VOLUME/2/blockdim4 + 1;
     }
     int griddim4 = gridsize;  
    update_constants_d(dev_grid);
    update_gpu_gf_d(g_gauge_field);
    

  #endif  
  
  
  
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
  		
  			#ifdef MPI
  			  if (g_proc_id == 0) {
  			#endif
  			
  			#ifdef MPI
  			  printf("\tOn host:\n");
  			  printf("\tVOLUME = %i\n", VOLUME);							// checking VOLUME and RAND in the parallel case 
  			  printf("\tRAND   = %i\n", RAND);
  			  printf("\tVOLUME + RAND = %i\n",  VOLUME+RAND);
  			#endif
  			
  			int host_check_LX, host_check_LY, host_check_LZ, host_check_T, host_check_VOLUME, host_check_OFFSET;
  			cudaMemcpyFromSymbol(&host_check_LX, dev_LX, sizeof(int));
  			cudaMemcpyFromSymbol(&host_check_LY, dev_LY, sizeof(int));
  			cudaMemcpyFromSymbol(&host_check_LZ, dev_LZ, sizeof(int));
  			cudaMemcpyFromSymbol(&host_check_T, dev_T, sizeof(int));
  			cudaMemcpyFromSymbol(&host_check_VOLUME, dev_VOLUME, sizeof(int));
			cudaMemcpyFromSymbol(&host_check_OFFSET, dev_Offset, sizeof(int));
  			// printf("\teven_odd_flag = %i\n", even_odd_flag);
  			printf("\tOn device:\n");
  			printf("\tdev_LX = %i\n", host_check_LX);
  			printf("\tdev_LY = %i\n", host_check_LY);
  			printf("\tdev_LZ = %i\n", host_check_LZ);
  			printf("\tdev_T = %i\n", host_check_T);
  			printf("\tdev_VOLUME = %i/2 ?!= %i\n", host_check_LX*host_check_LY*host_check_LZ*host_check_T, host_check_VOLUME);
  			printf("\tdev_Offset = %i\n", host_check_OFFSET);
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
  			
  			printf("\tOn host:\n");
  			printf("\tg_mubar = %f\n", g_mubar);
  			printf("\tg_epsbar = %f\n", g_epsbar);
  			
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
  
  
  
  #ifdef OPERATOR_BENCHMARK
    benchmark_eo_nd(Q_up, Q_dn, OPERATOR_BENCHMARK);
  #endif
  
  #ifdef GPU_DOUBLE
    #ifdef MATRIX_DEBUG
      test_double_nd_operator(Q_up, Q_dn, N_sites_int);
    #endif
  #endif
  /////////////////
  // ASSIGNMENTS //
  /////////////////
  
  
  x_up = P_up;							// can use the output spinors also as auxiliary fields
  x_dn = P_dn;							//	can use as initial guess at the same time
  

    r_up  = up_field[0];			// use the pre-allocated memory on host memory
    r_dn  = dn_field[0];			// allocated by  init_chi_spinor_field.c  and  invert_doublet.c  !?
    d_up  = up_field[1];		// the fields  up_field/dn_field[{0 , 1, ... , 5}]  are used in  cg_her_nd()
    d_dn  = dn_field[1];
    Ad_up = up_field[2];
    Ad_dn = dn_field[2];
    Ax_up = Ad_up;
    Ax_dn = Ad_dn;
    



  ///////////////
  // ALGORITHM //
  ///////////////
  

  printf("phmc_invmaxev = %f\n",phmc_invmaxev);  
  
  // timer
  startouter = gettime();
  
  #ifdef ALGORITHM_BENCHMARK
      starteffective = gettime();
  #endif
  
  #ifdef GPU_DOUBLE
    // r(0)
    // r(0) = b - A*x(0) = Q - A*P
      bb = square_norm(P_up, N_sites_int, 1) + square_norm(P_dn, N_sites_int, 1);
      order_spin_gpu(Q_up, h2d_spin_d);
      cudaMemcpy(Q_up_d, h2d_spin_d, dev_spinsize_d, cudaMemcpyHostToDevice);
      cudaMemcpy(r_up_d, h2d_spin_d, dev_spinsize_d, cudaMemcpyHostToDevice);      
      order_spin_gpu(Q_dn, h2d_spin_d);
      cudaMemcpy(Q_dn_d, h2d_spin_d, dev_spinsize_d, cudaMemcpyHostToDevice);
      cudaMemcpy(r_dn_d, h2d_spin_d, dev_spinsize_d, cudaMemcpyHostToDevice);
      
      //set solution accumulation fields to zero
      dev_zero_spinor_field_d<<<griddim4, blockdim4>>>(x_up_d);
      dev_zero_spinor_field_d<<<griddim4, blockdim4>>>(x_dn_d);      
      if (bb > 0.0) {
	  //set x_up/dn to initial guess in P_up/dn
	  printf("bb = %.16e \n", bb);
          order_spin_gpu(P_up, h2d_spin_d);
          cudaMemcpy(x_up_d, h2d_spin_d, dev_spinsize_d, cudaMemcpyHostToDevice);      
          order_spin_gpu(P_dn, h2d_spin_d);
          cudaMemcpy(x_dn_d, h2d_spin_d, dev_spinsize_d, cudaMemcpyHostToDevice);	  
          // r_up/dn = Q-A*x_up/dn
	  dev_Qtm_pm_ndpsi_d(Ax_up_d, Ax_dn_d,  
		      x_up_d, x_dn_d, 
		      griddim3, blockdim3, griddim4, blockdim4,
		      griddim4, blockdim4, griddim4, blockdim4);
          if(shift != 0.0) {
	    dev_axpy_d<<<griddim4,blockdim4>>>(shift, x_up_d, Ax_up_d);
	    dev_axpy_d<<<griddim4,blockdim4>>>(shift, x_dn_d, Ax_dn_d);
          }        
        
          dev_diff_d<<<griddim4,blockdim4>>>(r_up_d,Q_up_d,Ax_up_d);         
          dev_diff_d<<<griddim4,blockdim4>>>(r_dn_d,Q_dn_d,Ax_dn_d); 
	
     }

    
    // rr = (r_up)^2 + (r_dn)^2
    rr_up = double_dotprod(r_up_d, r_up_d);
    rr_dn = double_dotprod(r_dn_d, r_dn_d);
    rr = rr_up + rr_dn;
    rr_old = rr; // for the first iteration  
    
    r0r0   = double_dotprod(Q_up_d, Q_up_d) 
           + double_dotprod(Q_dn_d, Q_dn_d);    
  #else
      
      // P==0 -> no initial guess
      bb = square_norm(P_up, N_sites_int, 1) + square_norm(P_dn, N_sites_int, 1);
      if (bb > 0.0) {
	Qtm_pm_ndpsi(Ax_up, Ax_dn, P_up, P_dn);
	if(shift != 0.0) {
	  assign_add_mul_r(Ax_up, P_up , shift, N_sites_int);
	  assign_add_mul_r(Ax_dn, P_dn , shift, N_sites_int);
	}
	diff(r_up, Q_up, Ax_up, N_sites_int);
	diff(r_dn, Q_dn, Ax_dn, N_sites_int);
      }
      else {	
	assign(r_up, Q_up, N_sites_int);
	assign(r_dn, Q_dn, N_sites_int);
	//set solution accumulation fields to zero (==P_up/dn in this case)
	assign(x_up, P_up, N_sites_int);
	assign(x_dn, P_dn, N_sites_int);	
      }
    
    
    // rr = (r_up)^2 + (r_dn)^2
    rr_up = square_norm(r_up, N_sites_int, 1);
    rr_dn = square_norm(r_dn, N_sites_int, 1); 
    rr = rr_up + rr_dn;
    rr_old = rr; // for the first iteration
    //relative precision: norm of source
    r0r0   = square_norm(Q_up, N_sites_int, 1) + square_norm(Q_dn, N_sites_int, 1);   
  #endif //GPU_DOUBLE

  if (g_proc_id == 0) printf("Initial outer residue: %.10e\n", rr_old);
  if(rr_old < SP_MIN_EPS){
    if (g_proc_id == 0) printf("Initial residue too small for mixed precision! Stopping Inversion.\n");    
    finalize_mixedsolve_eo_nd();
    finalize_solver(up_field, nr_sf);
    finalize_solver(dn_field, nr_sf); 
    return(-1);
  }

  ////////////////
  // OUTER LOOP //
  ////////////////
  

  do {		// for (i = 0; i < max_iter; i++) {
    
    i++;
  
    #ifndef LOWOUTPUT
    if (g_proc_id == 0) printf("Outer iteration i = %i\n", i);
    #endif

    // host/device interaction    
    #ifdef GPU_DOUBLE
      dev_d2f<<<griddim4,blockdim4>>>(dev_spinin_up, r_up_d);
      dev_d2f<<<griddim4,blockdim4>>>(dev_spinin_dn, r_dn_d);     
    #else
      to_device(dev_spinin_up, r_up, h2d_spin_up, dev_spinsize_int);		// notice: for MPI communicateion the boundary exchange takes place when the hopping matrix is applied
      to_device(dev_spinin_dn, r_dn, h2d_spin_dn, dev_spinsize_int);
    #endif
    		// debug	// CUDA
    		#ifdef CUDA_DEBUG
    		  // CUDA_CHECK("CUDA error in mixedsolve_eo_nd(). Host to device interaction failed.", "Fields copied to device.");
    		  CUDA_CHECK_NO_SUCCESS_MSG("CUDA error in mixedsolve_eo_nd(). Host to device interaction failed.");
    		#endif
    
 
    ////////////////////////////////////
    // INNER LOOP, CONJUGATE GRADIENT //
    ////////////////////////////////////
    
    // SHIFT?
    // solves (A + shift)*p(k+1) = r(k)
    //        (A + shift)*p(0)   = r(0) = b
    float shift_single;
    if(shift==0.0){
      shift_single = 0.0f;
    }
    else{
      shift_single = (float) shift;
    }

    
    // timer
    startinner = gettime(); 
    innercount = dev_cg_eo_nd(dev_gf,
                          dev_spinout_up, dev_spinout_dn,
                          dev_spinin_up , dev_spinin_dn, shift_single,
                          max_innersolver_it,
                          innersolver_precision_check_abs, innersolver_precision_check_rel,
                          innersolver_precision_abs      , innersolver_precision_rel      );
    
    outercount = outercount + innercount;
    
    // timer
    stopinner = gettime();
    innerclocks = stopinner-startinner;
    #ifdef ALGORITHM_BENCHMARK
      effectiveflops = innercount*(matrixflops + 2*2*2*24 + 2*2*24 + 2*2*24 + 2*2*2*24 + 2*2*24)*VOLUME/2;   
      printf("inner solver BENCHMARK:\n");
      printf("\ttotal mixed solver time:   %.4e sec\n", innerclocks);
      printf("\tfloating point operations: %.4e flops\n", effectiveflops);
      printf("\tinner solver performance:  %.4e Gflop/s\n", double(effectiveflops)/innerclocks/ 1.0e9);
    #endif
      
    totalinnerclocks = totalinnerclocks + innerclocks;
    
    #ifndef LOWOUTPUT
    if (g_proc_id == 0) printf("Inner solver done in: %.4e sec\n", double(innerclocks));
    #endif

    #ifdef GPU_DOUBLE
      dev_add_f2d<<<griddim4,blockdim4>>>(x_up_d,x_up_d,dev_spinout_up); 
      dev_add_f2d<<<griddim4,blockdim4>>>(x_dn_d,x_dn_d,dev_spinout_dn);       

      dev_Qtm_pm_ndpsi_d(Ax_up_d, Ax_dn_d,  
		      x_up_d, x_dn_d, 
		      griddim3, blockdim3, griddim4, blockdim4,
		      griddim4, blockdim4, griddim4, blockdim4);
        if(shift != 0.0) {
	  dev_axpy_d<<<griddim4,blockdim4>>>(shift, x_up_d, Ax_up_d);
	  dev_axpy_d<<<griddim4,blockdim4>>>(shift, x_dn_d, Ax_dn_d);
        }        
        
      dev_diff_d<<<griddim4,blockdim4>>>(r_up_d,Q_up_d,Ax_up_d);         
      dev_diff_d<<<griddim4,blockdim4>>>(r_dn_d,Q_dn_d,Ax_dn_d); 
 
      rr_up = double_dotprod(r_up_d, r_up_d);
      rr_dn = double_dotprod(r_dn_d, r_dn_d);
      rr    = rr_up + rr_dn;
      
    #else
      // host/device interaction
      to_host(d_up, dev_spinout_up, h2d_spin_up, dev_spinsize_int);
      to_host(d_dn, dev_spinout_dn, h2d_spin_dn, dev_spinsize_int);
    
    		// debug	// CUDA
    		#ifdef CUDA_DEBUG
    		  // CUDA_CHECK("CUDA error in mixedsolve_eo_nd(). Device to host interaction failed.", "Fields copied back to device.");
    		  CUDA_CHECK_NO_SUCCESS_MSG("CUDA error in mixedsolve_eo_nd(). Device to host interaction failed.");
    		#endif
  
   

    // x(k+1) = x(k) + d(k+1)
      add(x_up, x_up, d_up, N_sites_int);
      add(x_dn, x_dn, d_dn, N_sites_int);

      // r(k+1)
      // r(k+1) = b - (A + shift)*x(k+1)
      // (A + shift)*x(k+1)
        Qtm_pm_ndpsi(Ax_up, Ax_dn, x_up, x_dn);
        if(shift != 0.0) {
	  assign_add_mul_r(Ax_up, x_up , shift, N_sites_int);
	  assign_add_mul_r(Ax_dn, x_dn , shift, N_sites_int);
        }          
        #ifndef LOWOUTPUT
        if (g_proc_id == 0) printf("The matrix was applied on CPU in double precision. r = b - Ax\n");
        #endif
        diff(r_up, Q_up, Ax_up, N_sites_int);
        diff(r_dn, Q_dn, Ax_dn, N_sites_int);

      // rr = (rr_up)^2 + (r_dn)^2
      rr_up = square_norm(r_up, N_sites_int, 1);
      rr_dn = square_norm(r_dn, N_sites_int, 1);
      rr    = rr_up + rr_dn;
    #endif //GPU_DOUBLE 

    if (g_proc_id == 0){ 
      printf("Outer residue at iteration i = %i : %.10e\n", i, rr);
      #ifndef LOWOUTPUT
       printf("No inner solver iterations: %i\n", outercount);
      #endif
    }

    // debug	// is NaN ?
    if isnan(rr) {
    	printf("Error in mixedsolve_eo_nd(). Outer residue is NaN.\n");
    	exit(-1);
    }
    		

    // aborting ?? // check wether precision is reached ...
    if ( ((rr <= eps_sq) && (rel_prec == 0))  ||  ((rr <= eps_sq*r0r0) && (rel_prec == 1)) ) {
 
     #ifdef GPU_DOUBLE
       //for GPU_DOUBLE we have to assign P_up/dn 
       //(for cpu outer solver this is not necessary, as P_up/dn == x_up/dn
        cudaMemcpy(h2d_spin_d, x_up_d, dev_spinsize_d, cudaMemcpyDeviceToHost);
        unorder_spin_gpu(h2d_spin_d, P_up); 
        cudaMemcpy(h2d_spin_d, x_dn_d, dev_spinsize_d, cudaMemcpyDeviceToHost);
        unorder_spin_gpu(h2d_spin_d, P_dn); 	
     #endif       
      
      // timer
      stopouter = gettime();
      totalouterclocks = stopouter-startouter - totalinnerclocks;
      if(g_debug_level > 3) printf("totalouterclocks = %d\n", totalouterclocks);
      #ifdef ALGORITHM_BENCHMARK
          stopeffective = gettime();
      #endif
      
      
      		// debug
      		#ifdef MPI
      		  if (g_proc_id == 0) {
      		#endif
      		printf("\nEO inversion done in mixed precision.\n");
      		if (rel_prec == 0) printf("Finished outer loop because of reached absolute outer solver precision.\n");
      		if (rel_prec == 1) printf("Finished outer loop because of reached relative outer solver precision.\n");
      		printf("Total number of inner iterations: %i\n", outercount);
      		printf("Total number of outer iterations: %i\n", i+1);
      		printf("Squared residue: %.10e\n", rr); 
      		printf("Outer solver done in: %.4e sec\n", double(stopouter-startouter));
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
      		  #endif
      		#endif
      		
      
      #ifdef USETEXTURE
        unbind_texture_gf();
      #endif
      
      		// debug	// CUDA
  		#ifdef CUDA_DEBUG
  		  CUDA_CHECK("CUDA error in unbind_texture(). Unbindung the GF texture failed.", "GF texture unbound.");
  		#endif
 
      finalize_mixedsolve_eo_nd();
      finalize_solver(up_field, nr_sf);
      finalize_solver(dn_field, nr_sf);  
      
      return(outercount);
      
    }
    //check if we can afford another inner solver run
    if(rr < SP_MIN_EPS){
      if (g_proc_id == 0) printf("At iteration %i: residue too small for mixed precision! Stopping Inversion.\n", outercount);    
      finalize_mixedsolve_eo_nd();
      finalize_solver(up_field, nr_sf);
      finalize_solver(dn_field, nr_sf); 
      return(-1);
    }  
  }//OUTER LOOP
  while (outercount <= max_iter);
  
  
  // multiplying with Qhat(2x2)^dagger is done in invert_doublet_eo.c
  
  
  // timer
  stopouter = clock();
  totalouterclocks = stopouter-startouter - totalinnerclocks;
  
  #ifdef ALGORITHM_BENCHMARK
      stopeffective = gettime();
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
      		printf("Outer solver done in: %.4e sec\n", double(stopouter-startouter));
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
      		  #endif
      		#endif
  
  
  #ifdef USETEXTURE
    unbind_texture_gf();
  #endif
      
      		// debug	// CUDA
  		#ifdef CUDA_DEBUG
  		  CUDA_CHECK("CUDA error in unbind_texture(). Unbindung the GF texture failed.", "GF texture unbound.");
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









#ifdef GPU_DOUBLE



void init_doublesolve_eo_nd (su3** gf) {	// gf is the full gauge field
  
  
  //////////////////////
  // GLOBAL VARIABLES //
  //////////////////////
  
  dev_spinsize_int   =  6*VOLUME/2*sizeof(dev_spinor);				// 24 floats per lattice site
  N_sites_int        =    VOLUME/2;
  N_floats_int       = 24*VOLUME/2;

  set_global_sizes();
  
  /////////////////////
  // LOCAL VARIABLES //
  /////////////////////
  
  cudaError_t cudaerr;		// CUDA errors
  int grid[6];			// array for grid specifications
  
  

  /////////////
  // SPINORS //							// allocates device memory for the odd part of the spinor fields (24 = 4*6 floats per odd lattice sites)
  /////////////							// now we have to consider 2 flavors: up, dn
  
  
  size_t dev_spinsize_d = 12*VOLUME/2 * sizeof(dev_spinor_d); /* double4 */
  //allocate fields used in dev_Qtm_pm_ndpsi_d
  cudaMalloc((void **) &dev_spin_eo1_up_d, dev_spinsize_d);	
  cudaMalloc((void **) &dev_spin_eo1_dn_d, dev_spinsize_d);	  
  cudaMalloc((void **) &dev_spin_eo3_up_d, dev_spinsize_d);	
  cudaMalloc((void **) &dev_spin_eo3_dn_d, dev_spinsize_d);
  //fields for cg
  //.. and allocate 6 additional
  cudaMalloc((void **) &dev_spin1_up_d, dev_spinsize_d);	
  cudaMalloc((void **) &dev_spin1_dn_d, dev_spinsize_d);	  
  cudaMalloc((void **) &dev_spin2_up_d, dev_spinsize_d);	
  cudaMalloc((void **) &dev_spin2_dn_d, dev_spinsize_d);	  
  cudaMalloc((void **) &dev_spin3_up_d, dev_spinsize_d);	
  cudaMalloc((void **) &dev_spin3_dn_d, dev_spinsize_d);
  
  cudaMalloc((void **) &dev_spinin_up_d, dev_spinsize_d);	
  cudaMalloc((void **) &dev_spinin_dn_d, dev_spinsize_d);
  cudaMalloc((void **) &dev_spinout_up_d, dev_spinsize_d);	
  cudaMalloc((void **) &dev_spinout_dn_d, dev_spinsize_d);
  
  if((cudaerr=cudaGetLastError())!=cudaSuccess){
      if(g_proc_id==0) printf("Error in init_doublesolve_eo_nd(): Memory allocation of nd additional double spinor fields failed. Aborting...\n");
      exit(200);
  }

  
  // debug	// CUDA
  #ifdef CUDA_DEBUG
    #ifndef MPI
      CUDA_CHECK("CUDA error in init_doublesolve_eo_nd(). Memory allocation of spinor fields failed.", "Allocated spinor fields on device.");
    #else
      CUDA_CHECK("CUDA error in init_doublesolve_eo_nd(). Memory allocation of spinor fields failed.", "Allocated spinor fields on devices.");
    #endif
  #endif

  
  ////////////////////////////
  // grid[ ] specifications //							// allocate and initializes the array grid[5] on device
  ////////////////////////////
  
  grid[0] = LX;									// it contains the dimensions of the lattice and the volume of the eo-sublattice
  grid[1] = LY;
  grid[2] = LZ;
  grid[3] = T;
  grid[4] = VOLUME/2;								// will be used to set dev_VOLUME: dev_VOLUME is half of VOLUME for eo
  
  grid[5] = VOLUME/2;

  //done in init_mixedsolve_eo
  //cudaMalloc((void **) &dev_grid, 6*sizeof(int));				// dev_grid
  cudaMemcpy(dev_grid, &(grid[0]), 6*sizeof(int), cudaMemcpyHostToDevice);
  
  // debug	// CUDA
  #ifdef CUDA_DEBUG
    #ifndef MPI
      CUDA_CHECK("CUDA error in init_mixedsolve_eo_nd(). Memory allocation of grid[] specifications failed.", "Allocated grid[] specifications on device.");
    #else
      CUDA_CHECK("CUDA error in init_mixedsolve_eo_nd(). Memory allocation of grid[] specifications failed.", "Allocated grid[] specifications on devices.");
    #endif
  #endif
  

}//init_doublesolve_eo_nd()



void finalize_doublesolve_eo_nd(void) {
  
  cudaError_t cudaerr;
   
  cudaFree(dev_spin_eo1_up_d);
  cudaFree(dev_spin_eo1_dn_d);     
  cudaFree(dev_spin_eo3_up_d);
  cudaFree(dev_spin_eo3_dn_d);
  cudaFree(dev_spin1_up_d);  
  cudaFree(dev_spin1_dn_d);
  cudaFree(dev_spin2_up_d);  
  cudaFree(dev_spin2_dn_d);
  cudaFree(dev_spin3_up_d);  
  cudaFree(dev_spin3_dn_d); 
  cudaFree(dev_spinin_up_d);  
  cudaFree(dev_spinin_dn_d);
  cudaFree(dev_spinout_up_d);  
  cudaFree(dev_spinout_dn_d);  
  #ifdef CUDA_DEBUG
    CUDA_CHECK("CUDA error in finalize_mixedsolve_eo_nd(). Device memory deallocation failed", "Device memory deallocated.");
  #endif
  
  
}


int dev_cg_eo_nd_d (dev_su3_2v_d * gf,
              dev_spinor_d * P_up, dev_spinor_d * P_dn,
              dev_spinor_d * Q_up, dev_spinor_d * Q_dn,
	      double shift,
              int max_iter,
              int check_abs , int check_rel,
              double eps_abs, double eps_rel, int min_solver_it       ) {
  
  // P_up/dn  can be used as auxiliary field to work on, as it is not later used (could be used as initial guess at the very start)
  // Q_up/dn  can be used as feedback, or if not, also as auxiliary field
  
  
  /////////////////////
  // LOCAL VARIABLES //
  /////////////////////
  
  // CUDA
  cudaError_t cudaerr;
  cublasStatus cublasstatus;
  
  // algorithm
  double rr_up, rr_dn, rr, rr_old, r0r0, dAd_up, dAd_dn, dAd;
  double alpha, beta;
  

  // (auxiliary) device fields
  // for recalculating the residue
  dev_spinor_d *  r_up, *  r_dn, * Ad_up, * Ad_dn, *  x_up, *  x_dn, *  d_up, *  d_dn, * Ax_up, * Ax_dn;		
  
 
  // counting
  int j;				// iteration counter
 
  // algorithm control parameters
  // int N_recalc_res = 10;		
  // recalculate residue r(k+1) = b - A*x(k+1) each N_recalc_res iteration
  int N_recalc_res = 1000;
  

  /////////////////////////////////////////////
  // CUDA block- and gridsize specifications //
  /////////////////////////////////////////////
  
  // int gridsize;		// auxiliary
  int blocksize;		// auxiliary
  
  blocksize = BLOCKD;
  int blockdim2, griddim2;					// passed:	dev_Hopping_Matrix
  if ( (VOLUME/2) % blocksize == 0 ) {
    blockdim2 = blocksize;
    griddim2  = VOLUME/2/blocksize;
  }
  else {
    blockdim2 = blocksize;
    griddim2  = (int) ((VOLUME/2/blocksize) + 1);
  }
  
  blocksize = BLOCK2D;
  int blockdim3, griddim3;					// passed:	dev_mul_one_pm_imubar_gamma5
  if ( (VOLUME/2) % blocksize == 0 ) {
    blockdim3 = blocksize;
    griddim3  = VOLUME/2/blocksize;
  }
  else {
    blockdim3 = blocksize;
    griddim3  = (int) ((VOLUME/2/blocksize) + 1);
  }
 

  /////////////////
  // ASSIGNMENTS //
  /////////////////
  
  x_up  = P_up;							
  // can use the output spinors also as auxiliary fields
  x_dn  = P_dn;							
  //	saves copying the output spinor field
  // use these pointers to the allocated space on device memory 
  // (allocated by init_mixedsolve_eo_nd_d)  
  r_up  = dev_spin1_up_d;						
  r_dn  = dev_spin1_dn_d;
  d_up  = dev_spin2_up_d;
  d_dn  = dev_spin2_dn_d;
  Ad_up = dev_spin3_up_d;
  Ad_dn = dev_spin3_dn_d;
  Ax_up = Ad_up;						
  // works as long as no initial guess vector x(0) is passed to cg_eo_nd_d()
  Ax_dn = Ad_dn;
  
  
 
  /////////////////////
  // INITIALIZATIONS //
  /////////////////////
  
      // debug	// CUBLAS helper function
      #ifdef CUDA_DEBUG
	CUBLAS_HELPER_CHECK(cublasInit(), "CUBLAS error in cublasInit(). Couldn't initialize CUBLAS.", "CUBLAS initialized.");
      #else
	#ifdef CUDA_45
	  cublasHandle_t handle;
	  cublasCreate(&handle);
	#else
	  cublasInit();
	#endif 
      #endif

  if(g_debug_level > 3) printf("cublasstatus = %f\n", cublasstatus);
  
  
  
  ///////////////
  // ALGORITHM //
  ///////////////
  
  
  // initialize x(0) = 0	// will be added up
  dev_zero_spinor_field_d<<<griddim2, blockdim2>>>(x_up);
  dev_zero_spinor_field_d<<<griddim2, blockdim2>>>(x_dn);
  
  // r(0) = b - A*x(0) = b
  dev_copy_spinor_field_d<<<griddim2, blockdim2>>>(Q_up, r_up);
  dev_copy_spinor_field_d<<<griddim2, blockdim2>>>(Q_dn, r_dn);
  
  // d(0) = r(0)
  dev_copy_spinor_field_d<<<griddim2, blockdim2>>>(r_up, d_up);
  dev_copy_spinor_field_d<<<griddim2, blockdim2>>>(r_dn, d_dn);
  
      // debug	// kernel
      #ifdef CUDA_DEBUG
	CUDA_KERNEL_CHECK("Kernel error in dev_cg_eo_nd(). Initializing spinor fields on device failed.", "Spinor fields initialized on device.");
      #endif
  
  
  
  // rr = (r_up)^2 + (r_dn)^2

  #ifdef MPI
    //launch error and exit
    if(g_proc_id == 0){
      printf("Error: pure double GPU ND solver not implemented for MPI. Aborting...\n");
      exit(200);
    }
  #endif
  
    rr_up = cublasDdot(N_floats_int, (double *) r_up, 1, (double *) r_up, 1);
    rr_dn = cublasDdot(N_floats_int, (double *) r_dn, 1, (double *) r_dn, 1);

  rr    = rr_up + rr_dn;

  r0r0   = rr;	// for relative precision 
  rr_old = rr;	// for the first iteration
  
  //////////
  // LOOP //
  //////////
  
      #ifndef LOWOUTPUT
      if (g_proc_id == 0) printf("\nEntering inner loop.\n");
      #endif

      // debug	// CUBLAS core function
      #ifdef CUDA_DEBUG
	// CUBLAS_CORE_CHECK("CUBLAS error in cg_eo_nd(). Calculating initial residue failed.", "Initial inner residue calculated.");
	CUBLAS_CORE_CHECK_NO_SUCCESS_MSG("CUBLAS error in dev_cg_eo_nd_d(). Calculating initial residue failed.");
      #endif

      if (g_proc_id == 0) printf("Initial inner residue: %.6e\n", r0r0);
  
  
  for (j = 0; j < max_iter; j++) {

     
      // A*d(k)

    dev_Qtm_pm_ndpsi_d(Ad_up, Ad_dn, d_up,  d_dn,										// debugging: matrix_debug1(), matrix_multiplication_test()
		      griddim2, blockdim2, griddim3, blockdim3,
		      griddim3, blockdim3, griddim3, blockdim3);

	if(shift != 0.0f){
           //add constant shift if nonzero
           // CUBLAS:
          cublasDaxpy (N_floats_int, shift, (double *) d_up, 1, (double *) Ad_up, 1);
          cublasDaxpy (N_floats_int, shift, (double *) d_dn, 1, (double *) Ad_dn, 1);
        }      
        if((cudaerr=cudaGetLastError()) != cudaSuccess){
           printf("%s\n", cudaGetErrorString(cudaerr));
           exit(200);
        }
      
	#ifdef CUDA_DEBUG
	  CUDA_CHECK("CUDA error in matrix_muliplication32(). Applying the matrix on GPU failed.", "The matrix was applied on GPU.");
	#endif

    // alpha = r(k)*r(k) / d(k)*A*d(k)
      dAd_up = cublasDdot(N_floats_int, (double *) d_up, 1, (double *) Ad_up, 1);
      dAd_dn = cublasDdot(N_floats_int, (double *) d_dn, 1, (double *) Ad_dn, 1);

    dAd    = dAd_up + dAd_dn;
    
    		// debug	// is NaN ?
    		if isnan(dAd) {
    		  printf("Error in dev_cg_eo_nd_d(). dAd is NaN.\n");
    		  exit(-1);
    		}
    
    alpha  = rr_old / dAd;	// rr_old is taken from the last iteration respectively
    
    // x(k+1) = x(k) + alpha*d(k)
    cublasDaxpy(N_floats_int, alpha, (double *) d_up, 1, (double *) x_up, 1);
    cublasDaxpy(N_floats_int, alpha, (double *) d_dn, 1, (double *) x_dn, 1);
    
    // r(k+1)
    if ( (j+1) % N_recalc_res != 0 ) {	// r(k+1) = r(k) - alpha*A*d(k)
      cublasDaxpy(N_floats_int, -1.0*alpha, (double *) Ad_up, 1, (double *) r_up, 1);
      cublasDaxpy(N_floats_int, -1.0*alpha, (double *) Ad_dn, 1, (double *) r_dn, 1);
    }
    else {				// recalculate residue r(k+1) = b - A*x(k+1)
	// debug
	#ifndef MPI
	  printf("Recalculating the inner residue.\n");
	#else
	  if (g_proc_id == 0) printf("Recalculating the inner residue.\n");
	#endif
      
      
      // A*x(k+1)


        dev_Qtm_pm_ndpsi_d(Ax_up, Ax_dn, x_up,  x_dn, 
        	           griddim2, blockdim2, griddim3, blockdim3,
        	           griddim3, blockdim3, griddim3, blockdim3);

	if(shift != 0.0){
           //add constant shift if nonzero
           // CUBLAS:
          cublasDaxpy (N_floats_int, shift, (double *) x_up, 1, (double *) Ad_up, 1);
          cublasDaxpy (N_floats_int, shift, (double *) x_dn, 1, (double *) Ad_dn, 1);
        }    

      // r(k+1) = b - A*x(k+1)
      cublasDcopy(N_floats_int, (double *) Q_up, 1, (double *) r_up, 1);		// r_up = Q_up
      cublasDcopy(N_floats_int, (double *) Q_dn, 1, (double *) r_dn, 1);		// r_dn = Q_dn
      cublasDaxpy(N_floats_int, -1.0, (double *) Ax_up, 1, (double *) r_up, 1);	// r_up = Q_up - Ax_up
      cublasDaxpy(N_floats_int, -1.0, (double *) Ax_dn, 1, (double *) r_dn, 1);	// r_dn = Q_dn - Ax_dn
    
    } // recalculate residue
    

    // r(k+1)*r(k+1)
      rr_up  = cublasDdot(N_floats_int, (double *) r_up, 1, (double *) r_up, 1);
      rr_dn  = cublasDdot(N_floats_int, (double *) r_dn, 1, (double *) r_dn, 1);

    rr     = rr_up + rr_dn;
    
    #ifdef CUDA_DEBUG
      CUBLAS_CORE_CHECK_NO_SUCCESS_MSG("CUBLAS error in dev_cg_eo_nd(). CUBLAS function failed.");
    #endif

    #ifndef LOWOUTPUT
      if (g_proc_id == 0) printf("inner iteration j = %i: rr = %.6e\n", j, rr);
    #endif
		 
	
    // debug	// is NaN ?
    if isnan(rr) {
      printf("Error in dev_cg_eo_nd(). Inner residue is NaN.\n");
      exit(-1);
    }
    
    // aborting ?? // check wether precision is reached ...
    if ( ((check_abs)&&(rr <= eps_abs)&&(j>min_solver_it)) || ((check_rel)&&(rr <= eps_rel*r0r0)&&(j>min_solver_it)) ) {
    

    
    if ((check_rel)&&(rr <= eps_rel*r0r0)) {
      printf("Reached relative inner solver precision of eps_rel = %.2e\n", eps_rel);
    }
    if ((check_abs)&&(rr <= eps_abs)) {
      printf("Reached absolute inner solver precision of eps_abs = %.2e\n", eps_abs);
    }
    
    printf("Final inner residue: %.6e\n", rr);
     
      #ifdef CUDA_45  
	cublasDestroy(handle);
      #else
	cublasShutdown();
      #endif 
      
      return(j+1);
    }
    
    
    // beta = r(k+1)*r(k+1) / r(k)*r(k)
    beta = rr / rr_old;
    
    
    rr_old = rr;  // for next iteration
    
    
    // d(k+1) = r(k+1) + beta*d(k)
    cublasDscal (N_floats_int, beta, (double *) d_up, 1);
    cublasDaxpy (N_floats_int, 1.0 , (double *) r_up, 1, (double *) d_up, 1);
    
    cublasDscal (N_floats_int, beta, (double *) d_dn, 1);
    cublasDaxpy (N_floats_int, 1.0 , (double *) r_dn, 1, (double *) d_dn, 1);
    
    // debug	// CUBLAS core function
    #ifdef CUDA_DEBUG
      CUBLAS_CORE_CHECK_NO_SUCCESS_MSG("CUBLAS error in dev_cg_eo_nd(). Error in CUBLAS function.");
    #endif
  }//LOOP
  
  
  #ifndef LOWOUTPUT
  if (g_proc_id == 0) printf("Finished inner loop because of maximal number of inner iterations.\n");
  #endif
  if (g_proc_id == 0) printf("Final inner residue: %.6e\n", rr);

    
    #ifdef CUDA_45  
      cublasDestroy(handle);
    #else
      cublasShutdown();
    #endif 
  return(j+1);
  
}//dev_cg_eo_nd_d()



//////////////////
// OUTER SOLVER //
//////////////////

// iterative refinement, defect correction
// that function is to replace the call of  cg_her_nd()  in  invert_doublet_eo.c
// solves the odd part of the full eo and nd problem
//	more precisely we have to invert  Qhat(2x2)*Qhat(2x2)^dagger
//	multiplying by  Qhat(2x2)^dagger  is done in  invert_doublet_eo.c

extern "C" int doublesolve_eo_nd (spinor * P_up, spinor * P_dn,
                                 spinor * Q_up, spinor * Q_dn, double shift,
                                 int max_iter, double eps_sq, int rel_prec, int min_solver_it) {
   
  
   if(rel_prec){
    innersolver_precision_check_rel = 1;
    innersolver_precision_check_abs = 0;
   }
   else{
    innersolver_precision_check_rel = 0;
    innersolver_precision_check_abs = 1;     
  }


  /////////////////////
  // LOCAL VARIABLES //
  /////////////////////
  
  // CUDA
  cudaError_t cudaerr;
  cublasStatus cublasstatus;
  
  // algorithm
  double rr_up, rr_dn, rr, rr_old, r0r0, bb;
  // counting
  int i = 0;					// iteration counter
  int innercount;				// latest inner solver iterations
  int outercount = 0;				// total inner solver iterations

  #ifdef ALGORITHM_BENCHMARK
    double effectiveflops;
    //double hoppingflops = 1608.0;
    double matrixflops  = 14448;
  #endif
  
  // timing
  double startouter, stopouter;
  double startinner, stopinner;  
  double innerclocks;
  
  // (auxiliary) fields
  spinor *  r_up, *  r_dn, * Ad_up, * Ad_dn, *  x_up, *  x_dn, *  d_up, *  d_dn, * Ax_up, * Ax_dn;

  spinor ** up_field = NULL;
  spinor ** dn_field = NULL;
  const int nr_sf = 5;

  init_solver_field(&up_field, VOLUMEPLUSRAND/2, nr_sf);
  init_solver_field(&dn_field, VOLUMEPLUSRAND/2, nr_sf);	 
	 

  //////////////////
  // INITIALIZING //
  //////////////////
  
  
      dev_spinor_d * x_up_d = dev_spin0_d;
      dev_spinor_d * x_dn_d = dev_spin1_d;
      dev_spinor_d * Ax_up_d = dev_spin2_d;
      dev_spinor_d * Ax_dn_d = dev_spin3_d;
      dev_spinor_d * Q_up_d = dev_spin_eo1_d;
      dev_spinor_d * Q_dn_d = dev_spin_eo2_d;  
      dev_spinor_d * r_up_d = dev_spin_eo3_up_d;  
      dev_spinor_d * r_dn_d = dev_spin_eo3_dn_d;      
      
       size_t dev_spinsize_d = 12*VOLUME/2 * sizeof(dev_spinor_d); // double4 even-odd !   
   int gridsize;
     //this is the partitioning for the HoppingMatrix kernel
     int blockdim3 = BLOCKD;
     if( VOLUME/2 % blockdim3 == 0){
       gridsize = (int) VOLUME/2/blockdim3;
     }
     else{
       gridsize = (int) VOLUME/2/blockdim3 + 1;
     }
     int griddim3 = gridsize;
   
     //this is the partitioning for dev_mul_one_pm...
     int blockdim4 = BLOCK2D;
     if( VOLUME/2 % blockdim4 == 0){
       gridsize = (int) VOLUME/2/blockdim4;
     }
     else{
       gridsize = (int) VOLUME/2/blockdim4 + 1;
     }
     int griddim4 = gridsize;  
    
     update_constants_d(dev_grid);
     update_gpu_gf_d(g_gauge_field);
    
  //do we need this after all??
  dev_complex h0, h1, h2, h3; 
  h0.re  =  (float) creal(ka0);	h0.im  = -(float) cimag(ka0);	// ka{0-4} are defined in boundary.c
  h1.re  =  (float) creal(ka1);	h1.im  = -(float) cimag(ka1);	// what is the meaning?
  h2.re  =  (float) creal(ka2);	h2.im  = -(float) cimag(ka2);
  h3.re  =  (float) creal(ka3);	h3.im  = -(float) cimag(ka3);
  he_cg_init<<< 1, 1 >>> (dev_grid, (float) g_kappa, (float)(g_mu/(2.0*g_kappa)), h0, h1, h2, h3);

  		
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
	    
	    int host_check_LX, host_check_LY, host_check_LZ, host_check_T, host_check_VOLUME, host_check_OFFSET;
	    cudaMemcpyFromSymbol(&host_check_LX, dev_LX, sizeof(int));
	    cudaMemcpyFromSymbol(&host_check_LY, dev_LY, sizeof(int));
	    cudaMemcpyFromSymbol(&host_check_LZ, dev_LZ, sizeof(int));
	    cudaMemcpyFromSymbol(&host_check_T, dev_T, sizeof(int));
	    cudaMemcpyFromSymbol(&host_check_VOLUME, dev_VOLUME, sizeof(int));
	    cudaMemcpyFromSymbol(&host_check_OFFSET, dev_Offset, sizeof(int));
	    // printf("\teven_odd_flag = %i\n", even_odd_flag);
	    printf("\tOn device:\n");
	    printf("\tdev_LX = %i\n", host_check_LX);
	    printf("\tdev_LY = %i\n", host_check_LY);
	    printf("\tdev_LZ = %i\n", host_check_LZ);
	    printf("\tdev_T = %i\n", host_check_T);
	    printf("\tdev_VOLUME = %i/2 ?!= %i\n", host_check_LX*host_check_LY*host_check_LZ*host_check_T, host_check_VOLUME);
	    printf("\tdev_Offset = %i\n", host_check_OFFSET);
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
  

    // debug	// check mubar and epsbar on host and device
    #ifdef STUFF_DEBUG
    
	    #ifdef MPI
	      if (g_proc_id == 0) {
	    #endif
	    
	    printf("\tOn host:\n");
	    printf("\tg_mubar = %f\n", g_mubar);
	    printf("\tg_epsbar = %f\n", g_epsbar);
	    
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
  

  #ifdef GPU_DOUBLE
   #ifdef MATRIX_DEBUG
    test_double_nd_operator(Q_up, Q_dn, N_sites_int);
   #endif
  #endif
  /////////////////
  // ASSIGNMENTS //
  /////////////////
  
  
  x_up = P_up;							// can use the output spinors also as auxiliary fields
  x_dn = P_dn;							//	can use as initial guess at the same time
//   x_up = up_field[3];							// can use the output spinors also as auxiliary fields
//   x_dn = dn_field[3];

    r_up  = up_field[0];			// use the pre-allocated memory on host memory
    r_dn  = dn_field[0];			// allocated by  init_chi_spinor_field.c  and  invert_doublet.c  !?
    d_up  = up_field[1];		// the fields  up_field/dn_field[{0 , 1, ... , 5}]  are used in  cg_her_nd()
    d_dn  = dn_field[1];
    Ad_up = up_field[2];
    Ad_dn = dn_field[2];
    Ax_up = Ad_up;
    Ax_dn = Ad_dn;
    

  ///////////////
  // ALGORITHM //
  ///////////////
  

  printf("phmc_invmaxev = %f\n",phmc_invmaxev);  
  
  
    // r(0)
    // r(0) = b - A*x(0) = Q - A*P
      bb = square_norm(P_up, VOLUME/2, 1);
      bb += square_norm(P_dn, VOLUME/2, 1);
      printf("bb = %.16e \n", bb);
      order_spin_gpu(Q_up, h2d_spin_d);
      cudaMemcpy(Q_up_d, h2d_spin_d, dev_spinsize_d, cudaMemcpyHostToDevice);      
      order_spin_gpu(Q_dn, h2d_spin_d);
      cudaMemcpy(Q_dn_d, h2d_spin_d, dev_spinsize_d, cudaMemcpyHostToDevice);	
      
//       if (bb > 0.0) {
// 	  printf("Have non-zero initial guess\n");
//       }
      
      //set x_up/dn to initial guess in P_up/dn
      order_spin_gpu(P_up, h2d_spin_d);
      cudaMemcpy(x_up_d, h2d_spin_d, dev_spinsize_d, cudaMemcpyHostToDevice);      
      order_spin_gpu(P_dn, h2d_spin_d);
      cudaMemcpy(x_dn_d, h2d_spin_d, dev_spinsize_d, cudaMemcpyHostToDevice);	  
      // r_up/dn = Q-A*x_up/dn
      dev_Qtm_pm_ndpsi_d(Ax_up_d, Ax_dn_d,  
		  x_up_d, x_dn_d, 
		  griddim3, blockdim3, griddim4, blockdim4,
		  griddim4, blockdim4, griddim4, blockdim4);
      if(shift != 0.0) {
	dev_axpy_d<<<griddim4,blockdim4>>>(shift, x_up_d, Ax_up_d);
	dev_axpy_d<<<griddim4,blockdim4>>>(shift, x_dn_d, Ax_dn_d);
      }        
    
      dev_diff_d<<<griddim4,blockdim4>>>(r_up_d,Q_up_d,Ax_up_d);         
      dev_diff_d<<<griddim4,blockdim4>>>(r_dn_d,Q_dn_d,Ax_dn_d); 
      
      cudaMemcpy(h2d_spin_d, r_up_d , dev_spinsize_d, cudaMemcpyDeviceToHost);      
      unorder_spin_gpu(h2d_spin_d, r_up);
      cudaMemcpy(h2d_spin_d, r_dn_d , dev_spinsize_d, cudaMemcpyDeviceToHost);      
      unorder_spin_gpu(h2d_spin_d, r_dn);
     
    
    // rr = (r_up)^2 + (r_dn)^2
    rr_up = double_dotprod(r_up_d, r_up_d);
    rr_dn = double_dotprod(r_dn_d, r_dn_d);
    rr = rr_up + rr_dn;
    rr_old = rr; // for the first iteration  
    
    r0r0   = double_dotprod(Q_up_d, Q_up_d) 
           + double_dotprod(Q_dn_d, Q_dn_d);    


  if (g_proc_id == 0) printf("Initial outer residue: %.10e\n", rr_old);


  ////////////////
  // OUTER LOOP //
  ////////////////

  do {		// for (i = 0; i < max_iter; i++) {
    
    i++;
  
    #ifndef LOWOUTPUT
    if (g_proc_id == 0) printf("Outer iteration i = %i\n", i);
    #endif

    order_spin_gpu(r_up, h2d_spin_d);
    cudaMemcpy(dev_spinin_up_d, h2d_spin_d, dev_spinsize_d, cudaMemcpyHostToDevice);      
    order_spin_gpu(r_dn, h2d_spin_d);
    cudaMemcpy(dev_spinin_dn_d, h2d_spin_d, dev_spinsize_d, cudaMemcpyHostToDevice);
 
    ////////////////////////////////////
    // INNER LOOP, CONJUGATE GRADIENT //
    ////////////////////////////////////
    

    startinner = gettime();
    innercount = dev_cg_eo_nd_d(dev_gf_d,
                          dev_spinout_up_d, dev_spinout_dn_d,
                          dev_spinin_up_d , dev_spinin_dn_d, shift,
                          max_innersolver_it,
                          innersolver_precision_check_abs, innersolver_precision_check_rel,
                          eps_sq     , eps_sq, min_solver_it    );
    //after first inner solve in double we merely check for absolute tolerance
    innersolver_precision_check_rel = 0;
    innersolver_precision_check_abs = 1;
    
    outercount = outercount + innercount;
    stopinner = gettime();
    innerclocks = stopinner-startinner;
    #ifdef ALGORITHM_BENCHMARK
      effectiveflops = innercount*(matrixflops + 2*2*2*24 + 2*2*24 + 2*2*24 + 2*2*2*24 + 2*2*24)*VOLUME/2;   
      printf("inner solver BENCHMARK:\n");
      printf("\tinner solver performance:  %.4e Gflop/s\n", double(effectiveflops)/innerclocks/ 1.0e9);
    #endif
      
    //copy result back
    cudaMemcpy(h2d_spin_d, dev_spinout_up_d , dev_spinsize_d, cudaMemcpyDeviceToHost);      
    unorder_spin_gpu(h2d_spin_d, d_up);
    cudaMemcpy(h2d_spin_d, dev_spinout_dn_d , dev_spinsize_d, cudaMemcpyDeviceToHost);      
    unorder_spin_gpu(h2d_spin_d, d_dn);
      
    		// debug	// CUDA
    		#ifdef CUDA_DEBUG
    		  // CUDA_CHECK("CUDA error in mixedsolve_eo_nd(). Device to host interaction failed.", "Fields copied back to device.");
    		  CUDA_CHECK_NO_SUCCESS_MSG("CUDA error in mixedsolve_eo_nd(). Device to host interaction failed.");
    		#endif
  
   

    // x(k+1) = x(k) + d(k+1)
      add(x_up, x_up, d_up, N_sites_int);
      add(x_dn, x_dn, d_dn, N_sites_int);

      // r(k+1)
      // r(k+1) = b - (A + shift)*x(k+1)
      // (A + shift)*x(k+1)
        Qtm_pm_ndpsi(Ax_up, Ax_dn, x_up, x_dn);
        if(shift != 0.0) {
	  assign_add_mul_r(Ax_up, x_up , shift, N_sites_int);
	  assign_add_mul_r(Ax_dn, x_dn , shift, N_sites_int);
        }          
        #ifndef LOWOUTPUT
        if (g_proc_id == 0) printf("The matrix was applied on CPU in double precision. r = b - Ax\n");
        #endif
        diff(r_up, Q_up, Ax_up, N_sites_int);
        diff(r_dn, Q_dn, Ax_dn, N_sites_int);

      // rr = (rr_up)^2 + (r_dn)^2
      rr_up = square_norm(r_up, N_sites_int, 1);
      rr_dn = square_norm(r_dn, N_sites_int, 1);
      rr    = rr_up + rr_dn;
      

    if (g_proc_id == 0){ 
      printf("Outer residue at iteration i = %i : %.10e\n", i, rr);
      #ifndef LOWOUTPUT
       printf("No inner solver iterations: %i\n", outercount);
      #endif
    }		
    // debug	// is NaN ?
    if isnan(rr) {
    	printf("Error in mixedsolve_eo_nd(). Outer residue is NaN.\n");
    	exit(-1);
    }
    		
    // aborting ?? // check wether precision is reached ...
    if ( ((rr <= eps_sq) && (rel_prec == 0))  ||  ((rr <= eps_sq*r0r0) && (rel_prec == 1)) ) {
      
      //we need not assign P_up/dn, as P_up/dn == x_up/dn

        
      #ifdef MPI
	if (g_proc_id == 0) {
      #endif
      printf("\nEO inversion done in double precision.\n");
      if (rel_prec == 0) printf("Finished outer loop because of reached absolute outer solver precision.\n");
      if (rel_prec == 1) printf("Finished outer loop because of reached relative outer solver precision.\n");
      printf("Total number of inner iterations: %i\n", outercount);
      printf("Total number of outer iterations: %i\n", i+1);
      printf("Squared residue: %.10e\n", rr); 
      #ifdef MPI
	}
      #endif
		  

      finalize_solver(up_field, nr_sf);
      finalize_solver(dn_field, nr_sf);  
      
      return(outercount);  
    }
    
  
  }//OUTER LOOP
  while (outercount <= max_iter);
  
  
  // multiplying with Qhat(2x2)^dagger is done in invert_doublet_eo.c
    #ifdef MPI
      if (g_cart_id == 0) {
    #endif
    printf("\nEO inversion done in full precision.\n");
    printf("Finished outer loop, because of maximal number of outer iterations.\n");
    printf("Total number of inner iterations: %i\n", outercount);
    printf("Total number of outer iterations: %i\n", i+1);
    printf("Squared residue: %.10e\n", rr); 
    #ifdef MPI
      }
    #endif
    


  finalize_solver(up_field, nr_sf);
  finalize_solver(dn_field, nr_sf);
  return(outercount);
  
  
}//doublesolve_eo_nd()



#endif











/**************************************************************************
 *
 * Copyright (C) 2010 Joseph Nagel
 *               2012, 2014 Florian Burger
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




extern "C" {
#include "../operator/tm_operators_nd.h"
#include "../operator/Hopping_Matrix.h"
#include "../solver/cg_her_nd.h"
}








// global formal parameters
size_t dev_gfsize;
size_t dev_spinsize_int;		// making the structure transparent:							
int N_sites_int;			// _int: internal sites
int N_floats_int;			// _ext: internal sites + additional boundaries
#ifdef _USE_MPI
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






#ifdef _USE_MPI					// collecting variables for the MPI implementation
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






#ifdef _USE_MPI

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


#ifdef _USE_MPI

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
  
  #ifndef _USE_MPI
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
  #ifdef _USE_MPI
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
  #ifdef _USE_MPI
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
  

  #ifndef _USE_MPI
  
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
    if(g_cart_id==0){
      printf("Error in init_mixedsolve_eo_nd(): Memory allocation of nd additional spinor fields failed.\n");
      printf("Error was %d. Aborting...\n", cudaerr);
    }
    exit(200);
  }
  
  
#ifdef GPU_DOUBLE
      #ifdef _USE_MPI
   	size_t dev_spinsize_d = 12*(VOLUME+RAND)/2 * sizeof(dev_spinor_d); /* double2 */
      #else
   	size_t dev_spinsize_d = 12*VOLUME/2 * sizeof(dev_spinor_d); /* double2 */  
      #endif

	  //allocate fields used in dev_Qtm_pm_ndpsi_d
  	  cudaMalloc((void **) &dev_spin_eo1_up_d, dev_spinsize_d);	
  	  cudaMalloc((void **) &dev_spin_eo1_dn_d, dev_spinsize_d);	  
  	  cudaMalloc((void **) &dev_spin_eo3_up_d, dev_spinsize_d);	
  	  cudaMalloc((void **) &dev_spin_eo3_dn_d, dev_spinsize_d);
	  
  	  if((cudaerr=cudaGetLastError())!=cudaSuccess){
              if(g_cart_id==0) printf("Error in init_mixedsolve_eo_nd(): Memory allocation of nd additional double spinor fields failed. Aborting...\n");
              exit(200);
          }
#endif
 
 
 
  #ifndef _USE_MPI
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
  		  printf("Process %d of %d: Could not allocate memory for h2d_spin_up. Aborting...\n", g_cart_id, g_nproc);	//                      they have to store floats (not doubles)
  		  exit(200);													//			can use "_int" ...
  		}														//			must use "_ext" when used with to_host_mpi as in xchange_field_wrapper()
  		
  		if ( (void *) (h2d_spin_dn = (dev_spinor *) malloc(dev_spinsize_ext) ) == NULL) {
  		  printf("Process %d of %d: Could not allocate memory for h2d_spin_dn. Aborting...\n", g_cart_id, g_nproc);
  		  exit(200);
  		}
  #endif
  
  
  #ifndef _USE_MPI
  
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
  		  #ifndef _USE_MPI
  		    CUDA_CHECK("CUDA error in init_mixedsolve_eo_nd(). Memory allocation of spinor fields failed.", "Allocated spinor fields on device.");
  		  #else
  		    CUDA_CHECK("CUDA error in init_mixedsolve_eo_nd(). Memory allocation of spinor fields failed.", "Allocated spinor fields on devices.");
  		  #endif
  		#endif

  for (int i = 0; i < 2; i++) {
        cudaStreamCreate(&stream_nd[i]);
  }   

  #ifdef _USE_MPI
  
  	#ifdef HOPPING_DEBUG													// Hopping_Matrix() is applied upon these spinor fields
  		// debug	// host code
  		if ( (void *) (spinor_debug_in = (spinor *) malloc(2*dev_spinsize_ext) ) == NULL) {
  		  printf("Process %d of %d: Could not allocate memory for spinor_debug_in. Aborting...\n", g_cart_id, g_nproc);
  		  exit(200);
  		}
  		// debug	// host code
  		if ( (void *) (spinor_debug_out = (spinor *) malloc(2*dev_spinsize_ext) ) == NULL) {
  		  printf("Process %d of %d: Could not allocate memory for spinor_debug_out. Aborting...\n", g_cart_id, g_nproc);
  		  exit(200);
  		}
  	#endif

  	int tSliceEO = LX*LY*LZ/2;

  
  	#if ASYNC > 0	// asynchronous communication and computation

  	  
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

	#ifdef GPU_DOUBLE
	  R1_UP_D = (dev_spinor_d *) malloc(2*tSliceEO*24*sizeof(double));
	  R2_UP_D = R1_UP_D + 12*tSliceEO;
	  R3_UP_D = (dev_spinor_d *) malloc(2*tSliceEO*24*sizeof(double));
	  R4_UP_D = R3_UP_D + 12*tSliceEO;


	//for gathering and spreading of indizes of rand in (gather_rand spread_rand called from xchange_field_wrapper)
	    #ifdef RELATIVISTIC_BASIS
	      cudaMalloc((void **) &RAND_FW_UP_D, tSliceEO*6*sizeof(double2));
	      cudaMalloc((void **) &RAND_BW_UP_D, tSliceEO*6*sizeof(double2));
	    #else
	      cudaMalloc((void **) &RAND_FW_UP_D, tSliceEO*12*sizeof(double2));
	      cudaMalloc((void **) &RAND_BW_UP_D, tSliceEO*12*sizeof(double2));      
	    #endif
	    /*  for async communication */
	    // page-locked memory    
	    #ifdef RELATIVISTIC_BASIS
	      int dbperspin = 6;
	    #else
	      int dbperspin = 12;
	    #endif  
	    cudaMallocHost(&RAND3_UP_D, tSliceEO*dbperspin*sizeof(double2));
	    cudaMallocHost(&RAND4_UP_D, tSliceEO*dbperspin*sizeof(double2));
	    cudaMallocHost(&RAND1_UP_D, tSliceEO*dbperspin*sizeof(double2));
	    cudaMallocHost(&RAND2_UP_D, tSliceEO*dbperspin*sizeof(double2));

	  R1_UP_D = (dev_spinor_d *) malloc(2*tSliceEO*24*sizeof(double));
	  R2_UP_D = R1_UP_D + 12*tSliceEO;
	  R3_UP_D = (dev_spinor_d *) malloc(2*tSliceEO*24*sizeof(double));
	  R4_UP_D = R3_UP_D + 12*tSliceEO;


	//for gathering and spreading of indizes of rand in (gather_rand spread_rand called from xchange_field_wrapper)
	    #ifdef RELATIVISTIC_BASIS
	      cudaMalloc((void **) &RAND_FW_DN_D, tSliceEO*6*sizeof(double2));
	      cudaMalloc((void **) &RAND_BW_DN_D, tSliceEO*6*sizeof(double2));
	    #else
	      cudaMalloc((void **) &RAND_FW_DN_D, tSliceEO*12*sizeof(double2));
	      cudaMalloc((void **) &RAND_BW_DN_D, tSliceEO*12*sizeof(double2));      
	    #endif

	    cudaMallocHost(&RAND3_DN_D, tSliceEO*dbperspin*sizeof(double2));
	    cudaMallocHost(&RAND4_DN_D, tSliceEO*dbperspin*sizeof(double2));
	    cudaMallocHost(&RAND1_DN_D, tSliceEO*dbperspin*sizeof(double2));
	    cudaMallocHost(&RAND2_DN_D, tSliceEO*dbperspin*sizeof(double2));
	    
	  R1_DN_D = (dev_spinor_d *) malloc(2*tSliceEO*24*sizeof(double));
	  R2_DN_D = R1_DN_D + 12*tSliceEO;
	  R3_DN_D = (dev_spinor_d *) malloc(2*tSliceEO*24*sizeof(double));
	  R4_DN_D = R3_DN_D + 12*tSliceEO;    
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
  #ifdef _USE_MPI
   grid[5] = (VOLUME+RAND)/2;
  #else
   grid[5] = VOLUME/2;
  #endif

  //done in init_mixedsolve_eo
  //cudaMalloc((void **) &dev_grid, 6*sizeof(int));				// dev_grid
  cudaMemcpy(dev_grid, &(grid[0]), 6*sizeof(int), cudaMemcpyHostToDevice);
  
  
  		// debug	// CUDA
  		#ifdef CUDA_DEBUG
  		  #ifndef _USE_MPI
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
  
  #ifdef _USE_MPI
  	
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
  
  #ifdef _USE_MPI
  	
  	#if ASYNC > 0
  	  
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
	
	#ifdef GPU_DOUBLE
	  cudaFreeHost(RAND1_UP_D);
	  cudaFreeHost(RAND2_UP_D); 
	  cudaFreeHost(RAND3_UP_D);
	  cudaFreeHost(RAND4_UP_D);  
	  cudaFree(RAND_BW_UP_D);
	  cudaFree(RAND_FW_UP_D);
	  free(R1_UP_D);
	  free(R3_UP_D);
	  cudaFreeHost(RAND1_DN_D);
	  cudaFreeHost(RAND2_DN_D); 
	  cudaFreeHost(RAND3_DN_D);
	  cudaFreeHost(RAND4_DN_D);  
	  cudaFree(RAND_BW_DN_D);
	  cudaFree(RAND_FW_DN_D);
	  free(R1_DN_D);
	  free(R3_DN_D); 
	#endif	    
	    
  #endif
  
  
  // Clean up CUDA API for calling thread	// ??
  //cudaThreadExit();				// is essential
  
  
  		// debug	// CUDA
  		#ifdef CUDA_DEBUG
		  cudaError_t cudaerr;	  
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

#ifdef _USE_MPI

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

void to_device (dev_spinor * device, spinor * host, dev_spinor * auxiliary, int evenodd) {
  int size;
  if(evenodd){
    #ifdef _USE_MPI
      size = 6*(VOLUME+RAND)/2*sizeof(dev_spinor);
    #else
      size = 6*VOLUME/2*sizeof(dev_spinor);
    #endif
  }
  else{
    #ifdef _USE_MPI
      size = 6*(VOLUME+RAND)*sizeof(dev_spinor);
    #else
      size = 6*VOLUME*sizeof(dev_spinor);
    #endif    
  }
  convert2REAL4_spin(host, auxiliary);						// auxiliary = (float) host
  cudaMemcpy(device, auxiliary, size, cudaMemcpyHostToDevice);			// device = auxiliary  (on device)

}


void to_host (spinor * host, dev_spinor * device, dev_spinor * auxiliary, int evenodd) {
  int size;
  if(evenodd){
    #ifdef _USE_MPI
      size = 6*(VOLUME+RAND)/2*sizeof(dev_spinor);
    #else
      size = 6*VOLUME/2*sizeof(dev_spinor);
    #endif
  }
  else{
    #ifdef _USE_MPI
      size = 6*(VOLUME+RAND)*sizeof(dev_spinor);
    #else
      size = 6*VOLUME*sizeof(dev_spinor);
    #endif    
  }
  cudaMemcpy(auxiliary, device, size, cudaMemcpyDeviceToHost);			// auxiliary = device  (on device)
  convert2double_spin(auxiliary, host);						// host = (double) auxiliary

}







////////////////////
// hopping matrix //
////////////////////

#ifdef _USE_MPI	// implemented for checking the MPI implementation of the hopping matrix
  #ifdef HOPPING_DEBUG

  // applies the hopping matrix on host for debugging purposes
  
  void Hopping_Matrix_wrapper (int ieo, dev_spinor * out, dev_spinor * in) {
  
      to_host(spinor_debug_in, in, h2d_spin_up, 1);
      Hopping_Matrix(ieo, spinor_debug_out, spinor_debug_in);
      to_device(out, spinor_debug_out, h2d_spin_dn, 1);    
    
  }

  #endif
#endif








void update_bare_constants_nd(){

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
  
  he_cg_init<<< 1, 1 >>> (dev_grid, (float) g_kappa, (float)(g_mu/(2.0*g_kappa)), h0, h1, h2, h3);
  		// "he" = "host entry"
  		// BEWARE in dev_tm_dirac_kappa we need the true mu (not 2 kappa mu!)	// ??
  		
  		// dev_LX, dev_LY, dev_LZ, dev_T, dev_VOLUME  =  grid[5]  =  dev_grid[5]
  		//	dev_VOLUME  is necessary for many kernel functions as for instance  dev_gamma5()
  		// initializes  mu, kappa and twokappamu  on the device
  		// initializes the strange  dev_k{0-3}, dev_mk{0-3}  as derived from the  ka{0-3}  from boundary.c
  #ifdef CUDA_DEBUG
    cudaError_t cudaerr;
    CUDA_KERNEL_CHECK("Kernel error in he_cg_init(). Couldn't initialize some stuff.", "he_cg_init() succeeded.");
  #endif
    
  he_cg_init_nd_additional<<<1,1>>> (g_mubar, g_epsbar);
  #ifdef CUDA_DEBUG
    CUDA_KERNEL_CHECK("Kernel error in he_cg_init_nd_additional(). Couldn't initialize some stuff.", "he_cg_init_nd_additional() succeeded.");
  #endif
  
  
  #ifdef _USE_MPI
    he_cg_init_nd_additional_mpi<<<1,1>>>(VOLUMEPLUSRAND/2, RAND, g_cart_id, g_nproc);	
    #ifdef CUDA_DEBUG
      CUDA_KERNEL_CHECK("Kernel error in he_cg_init_nd_additional_mpi(). Couldn't initialize some stuff.", "he_cg_init_nd_additional_mpi() succeeded.");
    #endif			
  #endif  
      
}	



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
  #ifndef _USE_MPI
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
  
  #ifndef _USE_MPI
    double effectiveDeviceFlops;
    double effectiveFlops;
  #else
    double effectiveDeviceFlops;
    double allEffectiveDeviceFlops;
    double effectiveFlops;
  #endif
  


  // formal parameters
  int staticsource = 0;		// 1: applies matrix every time on the same source
  				// 0: applies matrix consecutively ...
  

  dev_spinor * A_up;
  dev_spinor * A_dn;
  dev_spinor * B_up;
  dev_spinor * B_dn;
  
  dev_spinor * C_up;
  dev_spinor * C_dn;
  
  #ifndef _USE_MPI
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
		  cudaError_t cudaerr;
  		  CUDA_CHECK_NO_SUCCESS_MSG("CUDA error in benchmark_eo_nd(). Memory allocation of spinor fields failed.");
  		#endif

  
  		//debug
  		#ifndef _USE_MPI
  		  printf("\nStarting a little BENCHMARK. benchmark_eo_nd().\n");
  		#else
  		  if (g_cart_id == 0) printf("\nStarting a little BENCHMARK. benchmark_eo_nd_mpi().\n");
  		#endif

  
  		// debug
  		#ifndef _USE_MPI
  		  printf("Applying the eo-preconditioned matrix %i times.\n", N);
  		#else
  		  if (g_cart_id == 0) printf("Applying the eo-preconditioned matrix %i times.\n", N);
  		#endif
  
  
  to_device(B_up, Q_up, h2d_spin_up, 1);
  to_device(B_dn, Q_dn, h2d_spin_dn, 1);
  
  
  // timer
    startBenchmark = gettime();

  

  for (i = 0; i < N; i++) {
  
  
   
    dev_Qtm_pm_ndpsi(A_up, A_dn, B_up, B_dn, 
    	             gpu_gd_M, gpu_bd_M, gpu_gd_linalg, gpu_bd_linalg,
      		     gpu_gd_blas, gpu_bd_blas, gpu_gd_blas, gpu_bd_blas);

    
    
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

  
  
  #ifndef _USE_MPI
  
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
  	
  	
  	if (g_cart_id == 0) {

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



extern "C" void benchmark_eo_nd_d (spinor * Q_up, spinor * Q_dn, int N) {

  

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
  
  size_t dev_spinsize_d;
  

  // timing
  #ifndef _USE_MPI
    double timeElapsed;
    dev_spinsize_d = 12*VOLUME/2 * sizeof(dev_spinor_d); // double4 even-odd ! 
  #else
    double timeElapsed;    
    dev_spinsize_d = 12*(VOLUME+RAND)/2 * sizeof(dev_spinor_d); // double4 even-odd !     
  #endif
  double startBenchmark;
  double stopBenchmark;
  int staticsource = 0;
  // counter
  int i;
  
  // flop counting

  double effectiveFlopsPerApp = 14448.0; // hopping = 1608
  
  #ifndef _USE_MPI
    double effectiveDeviceFlops;
    double effectiveFlops;
  #else
    double effectiveDeviceFlops;
    double effectiveFlops;
  #endif
  


  dev_spinor_d * A_up = dev_spin_eo1_d;
  dev_spinor_d * A_dn = dev_spin_eo2_d;
  dev_spinor_d * B_up = dev_spin_eo1_up_d;
  dev_spinor_d * B_dn = dev_spin_eo1_dn_d;
  
  dev_spinor_d * C_up;
  dev_spinor_d * C_dn;
  
  
  		//debug
  		#ifndef _USE_MPI
  		  printf("\nStarting a little BENCHMARK. benchmark_eo_nd_d().\n");
  		#else
  		  if (g_cart_id == 0) printf("\nStarting a little BENCHMARK. benchmark_eo_nd_mpi().\n");
  		#endif

  
  		// debug
  		#ifndef _USE_MPI
  		  printf("Applying the eo-preconditioned matrix %i times.\n", N);
  		#else
  		  if (g_cart_id == 0) printf("Applying the eo-preconditioned matrix %i times.\n", N);
  		#endif
  
  
  order_spin_gpu(Q_up, h2d_spin_d);
  cudaMemcpy(B_up, h2d_spin_d, dev_spinsize_d, cudaMemcpyHostToDevice);      
  order_spin_gpu(Q_dn, h2d_spin_d);
  cudaMemcpy(B_dn, h2d_spin_d, dev_spinsize_d, cudaMemcpyHostToDevice);	
  order_spin_gpu(Q_up, h2d_spin_d);
  cudaMemcpy(A_up, h2d_spin_d, dev_spinsize_d, cudaMemcpyHostToDevice);      
  order_spin_gpu(Q_dn, h2d_spin_d);
  cudaMemcpy(A_dn, h2d_spin_d, dev_spinsize_d, cudaMemcpyHostToDevice);	  
  
  //startBenchmark = gettime();
  assert((startBenchmark = clock())!=-1);
  cudaThreadSynchronize();
  // timer
    

  for (i = 0; i < N; i++) {
  
  

    dev_Qtm_pm_ndpsi_d(A_up, A_dn, B_up, B_dn, 
    	                        gpu_gd_M_d, gpu_bd_M_d,
      		                gpu_gd_linalg_d, gpu_bd_linalg_d,
      		                gpu_gd_blas_d, gpu_bd_blas_d,
      		                gpu_gd_blas_d, gpu_bd_blas_d);

    
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
  if(g_cart_id==0){
    printf("%s\n", cudaGetErrorString(cudaGetLastError())); 
    printf("Done\n"); 
  }
  cudaThreadSynchronize();
  
  
  		// debug	// CUDA
  		#ifdef CUDA_DEBUG
		  cudaError_t cudaerr;
  		  CUDA_CHECK_NO_SUCCESS_MSG("CUDA error in nd matrix_muliplication(). Applying the matrix on GPU failed.");
  		#endif

  // timer
    //stopBenchmark = gettime();
    assert((stopBenchmark = clock())!=-1);
    
    if(g_cart_id == 0){
  	//timeElapsed = stopBenchmark - startBenchmark;
        timeElapsed = (double) (stopBenchmark-startBenchmark)/CLOCKS_PER_SEC;
  	effectiveDeviceFlops = g_nproc*N * VOLUME/2 * effectiveFlopsPerApp;
  	effectiveFlops       = double(g_nproc*N * VOLUME/2 * effectiveFlopsPerApp) / timeElapsed / 1.0e9;
  	printf("EFFECTIVE:\n");
  	printf("\ttime:        %.4e sec\n", timeElapsed);
  	printf("\tflop's:      %.4e flops\n", effectiveDeviceFlops);
  	printf("\tperformance: %.4e Gflop/s\n\n", effectiveFlops);
    }

  
}//benchmark_eo_nd_d()





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
    cublasStatus cublasstatus;  
    CUBLAS_HELPER_CHECK(cublasInit(), "CUBLAS error in cublasInit(). Couldn't initialize CUBLAS.", "CUBLAS initialized.");
    if(g_debug_level > 3) printf("cublasstatus = %f\n", cublasstatus);
  #else
    #ifdef CUDA_45
      cublasHandle_t handle;
      cublasCreate(&handle);
    #else
      cublasInit();
    #endif 
  #endif
  #ifdef _USE_MPI
      init_blas(VOLUME/2);
  #endif   

  
  
  
  ///////////////
  // ALGORITHM //
  ///////////////
  
  
  // initialize x(0) = 0	// will be added up
  dev_zero_spinor_field<<<gpu_gd_blas, gpu_bd_blas>>>(x_up);
  dev_zero_spinor_field<<<gpu_gd_blas, gpu_bd_blas>>>(x_dn);
  
  
  // r(0) = b - A*x(0) = b
  dev_copy_spinor_field<<<gpu_gd_blas, gpu_bd_blas>>>(Q_up, r_up);
  dev_copy_spinor_field<<<gpu_gd_blas, gpu_bd_blas>>>(Q_dn, r_dn);
  
  
  // d(0) = r(0)
  dev_copy_spinor_field<<<gpu_gd_blas, gpu_bd_blas>>>(r_up, d_up);
  dev_copy_spinor_field<<<gpu_gd_blas, gpu_bd_blas>>>(r_dn, d_dn);
  
  
  		// debug	// kernel
  		#ifdef CUDA_DEBUG
  		  CUDA_KERNEL_CHECK("Kernel error in dev_cg_eo_nd(). Initializing spinor fields on device failed.", "Spinor fields initialized on device.");
  		#endif
  
  
  
  
  // rr = (r_up)^2 + (r_dn)^2
  rr_up = cublasSdot_wrapper(N_floats_int, (float *) r_up, (float *) r_up);
  rr_dn = cublasSdot_wrapper(N_floats_int, (float *) r_dn, (float *) r_dn);

  rr    = rr_up + rr_dn;
  



  r0r0   = rr;	// for relative precision 
  rr_old = rr;	// for the first iteration
  
  
  
  
  //////////
  // LOOP //
  //////////
  
               
                #ifndef LOWOUTPUT
    		if (g_cart_id == 0) printf("\nEntering inner loop.\n");
	        #endif
  
  		// debug	// CUBLAS core function
		#ifdef CUDA_DEBUG
		  // CUBLAS_CORE_CHECK("CUBLAS error in cg_eo_nd(). Calculating initial residue failed.", "Initial inner residue calculated.");
		  CUBLAS_CORE_CHECK_NO_SUCCESS_MSG("CUBLAS error in dev_cg_eo_nd(). Calculating initial residue failed.");
		#endif
  
  		if (g_cart_id == 0) printf("Initial inner residue: %.6e\n", r0r0);
  
  
  for (j = 0; j < max_iter; j++) {

     
      // A*d(k)
      dev_Qtm_pm_ndpsi(Ad_up, Ad_dn, d_up,  d_dn,										// debugging: matrix_debug1(), matrix_multiplication_test()
      		       gpu_gd_M, gpu_bd_M, gpu_gd_linalg, gpu_bd_linalg,
      		       gpu_gd_blas, gpu_bd_blas, gpu_gd_blas, gpu_bd_blas);

      
      
	if(shift != 0.0f){
           //add constant shift if nonzero
           // CUBLAS:
          cublasSaxpy_wrapper (N_floats_int, shift, (float *) d_up, (float *) Ad_up);
          cublasSaxpy_wrapper (N_floats_int, shift, (float *) d_dn, (float *) Ad_dn);
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
    dAd_up = cublasSdot_wrapper(N_floats_int, (float *) d_up, (float *) Ad_up);
    dAd_dn = cublasSdot_wrapper(N_floats_int, (float *) d_dn, (float *) Ad_dn);

    dAd    = dAd_up + dAd_dn;
    
    		// debug	// is NaN ?
    		if isnan(dAd) {
    		  printf("Error in dev_cg_eo_nd(). dAd is NaN.\n");
    		  exit(-1);
    		}
    
    alpha  = rr_old / dAd;	// rr_old is taken from the last iteration respectively
    
    // x(k+1) = x(k) + alpha*d(k)
    cublasSaxpy_wrapper(N_floats_int, alpha, (float *) d_up, (float *) x_up);
    cublasSaxpy_wrapper(N_floats_int, alpha, (float *) d_dn, (float *) x_dn);
    
    // r(k+1)
    if ( (j+1) % N_recalc_res != 0 ) {	// r(k+1) = r(k) - alpha*A*d(k)
      cublasSaxpy_wrapper(N_floats_int, -1.0*alpha, (float *) Ad_up, (float *) r_up);
      cublasSaxpy_wrapper(N_floats_int, -1.0*alpha, (float *) Ad_dn, (float *) r_dn);
    }
    else {				// recalculate residue r(k+1) = b - A*x(k+1)
    					//	"feedback"
      		// debug
      		#ifndef _USE_MPI
      		  printf("Recalculating the inner residue.\n");
      		#else
      		  if (g_cart_id == 0) printf("Recalculating the inner residue.\n");
      		#endif
      
      
      // A*x(k+1)
        dev_Qtm_pm_ndpsi(Ax_up, Ax_dn, x_up,  x_dn, 
        	         gpu_gd_M, gpu_bd_M,gpu_gd_linalg, gpu_bd_linalg,
        	         gpu_gd_blas, gpu_bd_blas, gpu_gd_blas, gpu_bd_blas);
	

	if(shift != 0.0f){
           //add constant shift if nonzero
           // CUBLAS:
          cublasSaxpy_wrapper (N_floats_int, shift, (float *) x_up, (float *) Ad_up);
          cublasSaxpy_wrapper (N_floats_int, shift, (float *) x_dn, (float *) Ad_dn);
        }    

      // r(k+1) = b - A*x(k+1)
      dev_copy_spinor_field<<<gpu_gd_blas, gpu_bd_blas>>>(Q_up, r_up); // r_up = Q_up
      dev_copy_spinor_field<<<gpu_gd_blas, gpu_bd_blas>>>(Q_dn, r_dn); // r_dn = Q_dn 
      cublasSaxpy_wrapper(N_floats_int, -1.0, (float *) Ax_up, (float *) r_up);	// r_up = Q_up - Ax_up
      cublasSaxpy_wrapper(N_floats_int, -1.0, (float *) Ax_dn, (float *) r_dn);	// r_dn = Q_dn - Ax_dn
    
    } // recalculate residue
    

    // r(k+1)*r(k+1)
    rr_up  = cublasSdot_wrapper(N_floats_int, (float *) r_up, (float *) r_up);
    rr_dn  = cublasSdot_wrapper(N_floats_int, (float *) r_dn, (float *) r_dn);

    rr     = rr_up + rr_dn;
    
		// debug	// CUBLAS core function
		#ifdef CUDA_DEBUG
		  CUBLAS_CORE_CHECK_NO_SUCCESS_MSG("CUBLAS error in dev_cg_eo_nd(). CUBLAS function failed.");
		#endif
    
    
               #ifndef LOWOUTPUT
    		  if (g_cart_id == 0) printf("inner iteration j = %i: rr = %.6e\n", j, rr);
	       #endif
		 
	
    		// debug	// is NaN ?
    		if isnan(rr) {
    		  printf("Error in dev_cg_eo_nd(). Inner residue is NaN.\n");
    		  exit(-1);
    		}
    
    
    // aborting ?? // check wether precision is reached ...
    if ( ((check_abs)&&(rr <= eps_abs)) || ((check_rel)&&(rr <= eps_rel*r0r0)) ) {
    
      #ifdef _USE_MPI
        if (g_cart_id == 0) {
      #endif
      
      if ((check_rel)&&(rr <= eps_rel*r0r0)) {
      	printf("Reached relative inner solver precision of eps_rel = %.2e\n", eps_rel);
      }
      if ((check_abs)&&(rr <= eps_abs)) {
      	printf("Reached absolute inner solver precision of eps_abs = %.2e\n", eps_abs);
      }
      
      printf("Final inner residue: %.6e\n", rr);
      
      #ifdef _USE_MPI
        }
      #endif
      
      
      #ifdef CUDA_45  
	cublasDestroy(handle);
      #else
	cublasShutdown();
      #endif 
      #ifdef _USE_MPI
        finalize_blas();
      #endif  
      return(j+1);
    }
    
    
    // beta = r(k+1)*r(k+1) / r(k)*r(k)
    beta = rr / rr_old;
    
    
    rr_old = rr;  // for next iteration
    
    
    // d(k+1) = r(k+1) + beta*d(k)
    cublasSscal_wrapper (N_floats_int, beta, (float *) d_up);
    cublasSaxpy_wrapper (N_floats_int, 1.0 , (float *) r_up, (float *) d_up);
    
    cublasSscal_wrapper (N_floats_int, beta, (float *) d_dn);
    cublasSaxpy_wrapper (N_floats_int, 1.0 , (float *) r_dn, (float *) d_dn);
    
    		// debug	// CUBLAS core function
    		#ifdef CUDA_DEBUG
    		  CUBLAS_CORE_CHECK_NO_SUCCESS_MSG("CUBLAS error in dev_cg_eo_nd(). Error in CUBLAS function.");
    		#endif
    		
  
  }//LOOP
  
  
  #ifndef LOWOUTPUT
  if (g_cart_id == 0) printf("Finished inner loop because of maximal number of inner iterations.\n");
  #endif
  if (g_cart_id == 0) printf("Final inner residue: %.6e\n", rr);

    
  #ifdef CUDA_45  
    cublasDestroy(handle);
  #else
    cublasShutdown();
  #endif 
  #ifdef _USE_MPI
    finalize_blas();
  #endif  
  return(j+1);
  
}//dev_cg_eo_nd()





void test_double_nd_operator(spinor* const Q_up, spinor* const Q_dn, const int N){
   cudaError_t cudaerr;
   size_t dev_spinsize_d;
#ifdef _USE_MPI
  dev_spinsize_d = 12*(VOLUME+RAND)/2 * sizeof(dev_spinor_d); // double2 even-odd !   
#else
  dev_spinsize_d = 12*VOLUME/2 * sizeof(dev_spinor_d);  
#endif 
      dev_spinor_d * x_up_d = dev_spin0_d;
      dev_spinor_d * x_dn_d = dev_spin1_d;
      dev_spinor_d * Ax_up_d = dev_spin2_d;
      dev_spinor_d * Ax_dn_d = dev_spin3_d;    
      
  spinor ** solver_field_up = NULL;
  spinor ** solver_field_dn = NULL;  
  const int nr_sf = 3;
  init_solver_field(&solver_field_up, VOLUMEPLUSRAND/2, nr_sf);  
  init_solver_field(&solver_field_dn, VOLUMEPLUSRAND/2, nr_sf); 

  //apply cpu matrix

  
   Qtm_pm_ndpsi(solver_field_up[0], solver_field_dn[0], Q_up, Q_dn);
  
  /*
  dev_spinor_d * Q_up_d = dev_spin_eo1_d;
  dev_spinor_d * Q_dn_d = dev_spin_eo2_d; 
  Hopping_Matrix(EO,solver_field_up[1], Q_up);
  Hopping_Matrix(OE, solver_field_up[0] , solver_field_up[1]); 
  */
  //apply gpu matrix

  order_spin_gpu(Q_up, h2d_spin_d);
  cudaMemcpy(x_up_d, h2d_spin_d, dev_spinsize_d, cudaMemcpyHostToDevice);      
  order_spin_gpu(Q_dn, h2d_spin_d);
  cudaMemcpy(x_dn_d, h2d_spin_d, dev_spinsize_d, cudaMemcpyHostToDevice);	  
  
  #ifdef RELATIVISTIC_BASIS
      to_relativistic_basis_d<<<gpu_gd_linalg_d, gpu_bd_linalg_d>>> (x_up_d);
      to_relativistic_basis_d<<<gpu_gd_linalg_d, gpu_bd_linalg_d>>> (x_dn_d);
      if ((cudaerr=cudaGetLastError())!=cudaSuccess) {
	if (g_cart_id == 0) printf("%s\n", cudaGetErrorString(cudaerr));
      }
      else{
	#ifndef LOWOUTPUT 
	if (g_cart_id == 0) printf("Switched to relativistic basis\n");
	#endif
      }    
  #endif   
  
  // r_up/dn = Q-A*x_up/dn
  dev_Qtm_pm_ndpsi_d(Ax_up_d, Ax_dn_d,  
		      x_up_d, x_dn_d, 
		      gpu_gd_M_d, gpu_bd_M_d, gpu_gd_linalg_d, gpu_bd_linalg_d,
		      gpu_gd_linalg_d, gpu_bd_linalg_d, gpu_gd_linalg_d, gpu_bd_linalg_d);  
  /*
  dev_Hopp_d(x_dn_d, x_up_d,  
		      gpu_gd_M_d, gpu_bd_M_d, gpu_gd_linalg_d, gpu_bd_linalg_d,0); 
  dev_Hopp_d(Ax_up_d, x_dn_d,  
		      gpu_gd_M_d, gpu_bd_M_d, gpu_gd_linalg_d, gpu_bd_linalg_d,1);   
  */
  #ifdef RELATIVISTIC_BASIS 
    to_tmlqcd_basis_d<<<gpu_gd_linalg_d, gpu_bd_linalg_d>>> (Ax_up_d);
    to_tmlqcd_basis_d<<<gpu_gd_linalg_d, gpu_bd_linalg_d>>> (Ax_dn_d);      
  #endif   
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




void test_single_nd_operator(spinor* const Q_up, spinor* const Q_dn, const int N){
   cudaError_t cudaerr;

      dev_spinor * x_up = dev_spin1_up;
      dev_spinor * x_dn = dev_spin1_dn;
      dev_spinor * Ax_up = dev_spin2_up;
      dev_spinor * Ax_dn = dev_spin2_dn;    
      
  spinor ** solver_field_up = NULL;
  spinor ** solver_field_dn = NULL;  
  const int nr_sf = 3;
  init_solver_field(&solver_field_up, VOLUMEPLUSRAND/2, nr_sf);  
  init_solver_field(&solver_field_dn, VOLUMEPLUSRAND/2, nr_sf); 

  //apply cpu matrix 
   Qtm_pm_ndpsi(solver_field_up[0], solver_field_dn[0], Q_up, Q_dn);
  

  //apply gpu matrix
  to_device(x_up, Q_up, h2d_spin, 1);
  to_device(x_dn, Q_dn, h2d_spin, 1);  
  #ifdef RELATIVISTIC_BASIS
      to_relativistic_basis<<<gpu_gd_linalg, gpu_bd_linalg>>> (x_up);
      to_relativistic_basis<<<gpu_gd_linalg, gpu_bd_linalg>>> (x_dn);
      if ((cudaerr=cudaGetLastError())!=cudaSuccess) {
	if (g_cart_id == 0) printf("%s\n", cudaGetErrorString(cudaerr));
      }
      else{
	#ifndef LOWOUTPUT 
	if (g_cart_id == 0) printf("Switched to relativistic basis\n");
	#endif
      }    
  #endif   
  
  // r_up/dn = Q-A*x_up/dn
  dev_Qtm_pm_ndpsi(Ax_up, Ax_dn,  
		      x_up, x_dn, 
		      gpu_gd_M, gpu_bd_M, gpu_gd_linalg, gpu_bd_linalg,
		      gpu_gd_linalg, gpu_bd_linalg, gpu_gd_linalg, gpu_bd_linalg);  

  #ifdef RELATIVISTIC_BASIS 
    to_tmlqcd_basis<<<gpu_gd_linalg, gpu_bd_linalg>>> (Ax_up);
    to_tmlqcd_basis<<<gpu_gd_linalg, gpu_bd_linalg>>> (Ax_dn);      
  #endif   
  to_host(solver_field_up[1], Ax_up, h2d_spin, 1);
  to_host(solver_field_dn[1], Ax_dn, h2d_spin, 1); 
      
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
  printf("Testing single matrix:\n");
  printf("cpu: Squared difference is   UP: %.8e\n", rk_up);
  printf("cpu: Squared difference is DOWN: %.8e\n", rk_dn);  
  printf("cpu: Squared difference per spinor component is: %.8e\n", rk/N/24.0);  
  printf("Max. difference at position %i: %.8e\n", at_max, max_dev);
  printf("Min. difference at position %i: %.8e\n", at_min, min_dev);  

  finalize_solver(solver_field_up, nr_sf); 
  finalize_solver(solver_field_dn, nr_sf);  
}



void check_nd_mixedsolve_params(){
  
    #ifdef _USE_MPI
      if (g_cart_id == 0) {
    #endif
    
    #ifdef _USE_MPI
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
    
    #ifdef _USE_MPI
      }
    #endif

    // check mubar and epsbar on host and device

    #ifdef _USE_MPI
      if (g_cart_id == 0) {
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
  
    #ifdef _USE_MPI			  
      int host_check_VOLUMEPLUSRAND, host_check_RAND;
      cudaMemcpyFromSymbol(&host_check_VOLUMEPLUSRAND, dev_VOLUMEPLUSRAND, sizeof(int));
      cudaMemcpyFromSymbol(&host_check_RAND, dev_RAND, sizeof(int));
      printf("\tOn device:\n");
      printf("\tdev_VOLUMEPLUSRAND = %i\n", host_check_VOLUMEPLUSRAND);
      printf("\tdev_RAND = %i\n", host_check_RAND);    
    #endif
    #ifdef _USE_MPI
      }
    #endif

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

  
  /////////////////////
  // LOCAL VARIABLES //
  /////////////////////
  
  
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
    #ifdef _USE_MPI
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
    #ifndef _USE_MPI
      clock_t starteffective;
      clock_t stopeffective;
    #else
      double starteffective;
      double stopeffective;
      double singletime;				// time for each process = stopeffective - starteffective
      double maxtime;				// max. parallel process time
    #endif
  #endif
  


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
  

  //update the gpu single gauge_field
  update_gpu_gf(g_gauge_field); 
  set_gpu_work_layout(1); //set block and grid sizes, eo!
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

      size_t dev_spinsize_d;
      #ifdef _USE_MPI
       dev_spinsize_d = 12*(VOLUME+RAND)/2 * sizeof(dev_spinor_d); // double2 even-odd !      
      #else
       dev_spinsize_d = 12*VOLUME/2 * sizeof(dev_spinor_d);   
      #endif
       update_constants_d(dev_grid);
       update_gpu_gf_d(g_gauge_field);
  #endif  
  
  
  
  // bind texture gf
  #ifdef USETEXTURE						// needed for subfunctions of dev_Hopping_Matrix(...)
    bind_texture_gf(dev_gf);					//	e.g. dev_reconstructgf_2vtexref(...), dev_reconstructgf_8texref(...), 
  #endif							//	in general in functions  dev_reconstructgf[...]  with  "tex1Dfetch(gf_tex[...]"
  #ifdef CUDA_DEBUG
    cudaError_t cudaerr;
    CUDA_CHECK("CUDA error in bind_texture_gf(). Binding GF to texture failed.", "GF bound to texture.");
  #endif
    
      
  //update the constants on device
    update_bare_constants_nd();
    
  #ifdef STUFF_DEBUG
    // check params on device
    check_nd_mixedsolve_params();
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
  #ifndef GPU_DOUBLE
  // (auxiliary) fields in case we of no GPU_DOUBLE
    spinor *r_up, *r_dn, *Ad_up, *Ad_dn, *x_up, *x_dn, *d_up, *d_dn, *Ax_up, *Ax_dn;      
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
  #endif  



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
      #ifdef RELATIVISTIC_BASIS
        to_relativistic_basis_d<<<gpu_gd_linalg_d, gpu_bd_linalg_d>>> (Q_up_d);
        to_relativistic_basis_d<<<gpu_gd_linalg_d, gpu_bd_linalg_d>>> (Q_dn_d);   
        to_relativistic_basis_d<<<gpu_gd_linalg_d, gpu_bd_linalg_d>>> (r_up_d);
        to_relativistic_basis_d<<<gpu_gd_linalg_d, gpu_bd_linalg_d>>> (r_dn_d);	
      #endif 
      //set solution accumulation fields to zero
      dev_zero_spinor_field_d<<<gpu_gd_linalg_d, gpu_bd_linalg_d>>>(x_up_d);
      dev_zero_spinor_field_d<<<gpu_gd_linalg_d, gpu_bd_linalg_d>>>(x_dn_d);      
      if (bb > 0.0) {
	  //set x_up/dn to initial guess in P_up/dn
	  printf("bb = %.16e \n", bb);
          order_spin_gpu(P_up, h2d_spin_d);
          cudaMemcpy(x_up_d, h2d_spin_d, dev_spinsize_d, cudaMemcpyHostToDevice);      
          order_spin_gpu(P_dn, h2d_spin_d);
          cudaMemcpy(x_dn_d, h2d_spin_d, dev_spinsize_d, cudaMemcpyHostToDevice);
          #ifdef RELATIVISTIC_BASIS  
            to_relativistic_basis_d<<<gpu_gd_linalg_d, gpu_bd_linalg_d>>> (x_up_d);
            to_relativistic_basis_d<<<gpu_gd_linalg_d, gpu_bd_linalg_d>>> (x_dn_d);	
          #endif 	  
          // r_up/dn = Q-A*x_up/dn
	  dev_Qtm_pm_ndpsi_d(Ax_up_d, Ax_dn_d,  
		      x_up_d, x_dn_d, 
		      gpu_gd_M_d, gpu_bd_M_d, gpu_gd_linalg_d, gpu_bd_linalg_d,
		      gpu_gd_linalg_d, gpu_bd_linalg_d, gpu_gd_linalg_d, gpu_bd_linalg_d);
          if(shift != 0.0) {
	    dev_axpy_d<<<gpu_gd_linalg_d,gpu_bd_linalg_d>>>(shift, x_up_d, Ax_up_d);
	    dev_axpy_d<<<gpu_gd_linalg_d,gpu_bd_linalg_d>>>(shift, x_dn_d, Ax_dn_d);
          }        
        
          dev_diff_d<<<gpu_gd_linalg_d,gpu_bd_linalg_d>>>(r_up_d,Q_up_d,Ax_up_d);         
          dev_diff_d<<<gpu_gd_linalg_d,gpu_bd_linalg_d>>>(r_dn_d,Q_dn_d,Ax_dn_d); 
	
     }

    
    // rr = (r_up)^2 + (r_dn)^2
    rr_up = double_dotprod(r_up_d, r_up_d, N_sites_int);
    rr_dn = double_dotprod(r_dn_d, r_dn_d, N_sites_int);
    rr = rr_up + rr_dn;
    rr_old = rr; // for the first iteration  
    
    r0r0   = double_dotprod(Q_up_d, Q_up_d, N_sites_int) 
           + double_dotprod(Q_dn_d, Q_dn_d, N_sites_int);    
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

  if (g_cart_id == 0) printf("Initial outer residue: %.10e\n", rr_old);
  if(rr_old < SP_MIN_EPS){
    if (g_cart_id == 0) printf("Initial residue too small for mixed precision! Stopping Inversion.\n");    
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
    if (g_cart_id == 0) printf("Outer iteration i = %i\n", i);
    #endif

    // host/device interaction    
    #ifdef GPU_DOUBLE
      dev_d2f<<<gpu_gd_linalg_d,gpu_bd_linalg_d>>>(dev_spinin_up, r_up_d);
      dev_d2f<<<gpu_gd_linalg_d,gpu_bd_linalg_d>>>(dev_spinin_dn, r_dn_d);     
    #else
      to_device(dev_spinin_up, r_up, h2d_spin_up, 1);
      to_device(dev_spinin_dn, r_dn, h2d_spin_dn, 1);
      #ifdef RELATIVISTIC_BASIS  
        to_relativistic_basis<<<gpu_gd_linalg, gpu_bd_linalg>>> (dev_spinin_up);
        to_relativistic_basis<<<gpu_gd_linalg, gpu_bd_linalg>>> (dev_spinin_dn);	
      #endif 
    #endif
    		
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
      if ((g_cart_id == 0) && (g_debug_level > 1)){
      effectiveflops = innercount*(matrixflops + 2*2*2*24 + 2*2*24 + 2*2*24 + 2*2*2*24 + 2*2*24)*VOLUME/2;   
      printf("inner solver BENCHMARK:\n");
      printf("\ttotal mixed solver time:   %.4e sec\n", innerclocks);
      printf("\tfloating point operations: %.4e flops\n", effectiveflops*g_nproc);
      printf("\tinner solver performance:  %.4e Gflop/s\n", double(effectiveflops*g_nproc)/innerclocks/ 1.0e9);
      }
    #endif
      
    totalinnerclocks = totalinnerclocks + innerclocks;
    
    #ifndef LOWOUTPUT
    if (g_cart_id == 0) printf("Inner solver done in: %.4e sec\n", double(innerclocks));
    #endif

    #ifdef GPU_DOUBLE
      dev_add_f2d<<<gpu_gd_linalg_d,gpu_bd_linalg_d>>>(x_up_d,x_up_d,dev_spinout_up); 
      dev_add_f2d<<<gpu_gd_linalg_d,gpu_bd_linalg_d>>>(x_dn_d,x_dn_d,dev_spinout_dn);       

      dev_Qtm_pm_ndpsi_d(Ax_up_d, Ax_dn_d,  
		      x_up_d, x_dn_d, 
		      gpu_gd_M_d, gpu_bd_M_d, gpu_gd_linalg_d, gpu_bd_linalg_d,
		      gpu_gd_linalg_d, gpu_bd_linalg_d, gpu_gd_linalg_d, gpu_bd_linalg_d);
        if(shift != 0.0) {
	  dev_axpy_d<<<gpu_gd_linalg_d,gpu_bd_linalg_d>>>(shift, x_up_d, Ax_up_d);
	  dev_axpy_d<<<gpu_gd_linalg_d,gpu_bd_linalg_d>>>(shift, x_dn_d, Ax_dn_d);
        }        
        
      dev_diff_d<<<gpu_gd_linalg_d,gpu_bd_linalg_d>>>(r_up_d,Q_up_d,Ax_up_d);         
      dev_diff_d<<<gpu_gd_linalg_d,gpu_bd_linalg_d>>>(r_dn_d,Q_dn_d,Ax_dn_d); 
 
      rr_up = double_dotprod(r_up_d, r_up_d, N_sites_int);
      rr_dn = double_dotprod(r_dn_d, r_dn_d, N_sites_int);
      rr    = rr_up + rr_dn;
      
    #else
      // host/device interaction
      #ifdef RELATIVISTIC_BASIS 
        to_tmlqcd_basis<<<gpu_gd_linalg, gpu_bd_linalg>>> (dev_spinout_up);
        to_tmlqcd_basis<<<gpu_gd_linalg, gpu_bd_linalg>>> (dev_spinout_dn);      
      #endif       
      to_host(d_up, dev_spinout_up, h2d_spin_up, 1);
      to_host(d_dn, dev_spinout_dn, h2d_spin_dn, 1);

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
        if (g_cart_id == 0) printf("The matrix was applied on CPU in double precision. r = b - Ax\n");
        #endif
        diff(r_up, Q_up, Ax_up, N_sites_int);
        diff(r_dn, Q_dn, Ax_dn, N_sites_int);

      // rr = (rr_up)^2 + (r_dn)^2
      rr_up = square_norm(r_up, N_sites_int, 1);
      rr_dn = square_norm(r_dn, N_sites_int, 1);
      rr    = rr_up + rr_dn;
    #endif //GPU_DOUBLE 

    if (g_cart_id == 0){ 
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
        #ifdef RELATIVISTIC_BASIS 
          to_tmlqcd_basis_d<<<gpu_gd_linalg_d, gpu_bd_linalg_d>>> (x_up_d);
          to_tmlqcd_basis_d<<<gpu_gd_linalg_d, gpu_bd_linalg_d>>> (x_dn_d);      
        #endif              
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
      		#ifdef _USE_MPI
      		  if (g_cart_id == 0) {
      		#endif
      		printf("EO inversion done in mixed precision.\n");
      		if (rel_prec == 0) printf("Finished outer loop because of reached absolute outer solver precision.\n");
      		if (rel_prec == 1) printf("Finished outer loop because of reached relative outer solver precision.\n");
      		printf("Total number of inner iterations: %i\n", outercount);
      		printf("Total number of outer iterations: %i\n", i+1);
      		printf("Squared residue: %.10e\n", rr); 
      		printf("Outer solver done in: %.4e sec\n\n", double(stopouter-startouter));
      		#ifdef _USE_MPI
      		  }
      		#endif
      		
      		// benchmark
      		#ifdef ALGORITHM_BENCHMARK
      		  // will now count the number of effective flops
      		  // effectiveflops  =  #(inner iterations)*(matrixflops+linalgflops)*VOLUME/2  +  #(outer iterations)*(matrixflops+linalgflops)*VOLUME/2
      		  // outer loop: linalg  =  flops for calculating  r(k+1) and x(k+1)
      		  // inner loop: linalg  =  flops for calculating  alpha, x(k+1), r(k+1), beta, d(k+1)
      		  #ifndef _USE_MPI
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
      		  	if (g_cart_id == 0) printf("effective BENCHMARK:\n");
      		  	if (g_cart_id == 0) printf("\ttotal mixed solver time:   %.4e sec\n", double(maxtime));
      		  	if (g_cart_id == 0) printf("\tfloating point operations: %.4e flops\n", double(allflops));
      		  	if (g_cart_id == 0) printf("\tinner solver performance:  %.4e Gflop/s\n", double(allflops) / double(maxtime) / 1.0e9);
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
      if (g_cart_id == 0) printf("At iteration %i: residue too small for mixed precision! Stopping Inversion.\n", outercount);    
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
  		#ifdef _USE_MPI
  		  if (g_cart_id == 0) {
  		#endif
  		printf("EO inversion done in mixed precision.\n");
  		printf("Finished outer loop, because of maximal number of outer iterations.\n");
      		printf("Total number of inner iterations: %i\n", outercount);
      		printf("Total number of outer iterations: %i\n", i+1);
      		printf("Squared residue: %.10e\n", rr); 
      		printf("Outer solver done in: %.4e sec\n\n", double(stopouter-startouter));
      		#ifdef _USE_MPI
      		  }
      		#endif
      		
      		// benchmark
      		#ifdef ALGORITHM_BENCHMARK
      		  // will now count the number of effective flops
      		  // effectiveflops  =  #(inner iterations)*(matrixflops+linalgflops)*VOLUME/2  +  #(outer iterations)*(matrixflops+linalgflops)*VOLUME/2
      		  // outer loop: linalg  =  flops for calculating  r(k+1) and x(k+1)
      		  // inner loop: linalg  =  flops for calculating  alpha, x(k+1), r(k+1), beta, d(k+1)
      		  #ifndef _USE_MPI
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
      		  	if (g_cart_id == 0) printf("effective BENCHMARK:\n");
      		  	if (g_cart_id == 0) printf("\ttotal mixed solver time:   %.4e sec\n", double(maxtime));
      		  	if (g_cart_id == 0) printf("\tfloating point operations: %.4e flops\n", double(allflops));
      		  	if (g_cart_id == 0) printf("\tinner solver performance:  %.4e Gflop/s\n", double(allflops) / double(maxtime) / 1.0e9);
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
  		#ifndef _USE_MPI
  		  printf("\n");
  		#else
  		  if (g_cart_id == 0) printf("\n");
  		#endif

  finalize_solver(up_field, nr_sf);
  finalize_solver(dn_field, nr_sf);
  return(outercount);
  
  
}//mixedsolve_eo_nd()












void init_doublesolve_eo_nd (su3** gf) {	// gf is the full gauge field
  

  
  /////////////////////
  // LOCAL VARIABLES //
  /////////////////////
  
  int grid[6];			// array for grid specifications
  cudaError_t cudaerr;		// CUDA errors 
  

  /////////////
  // SPINORS //							// allocates device memory for the odd part of the spinor fields (24 = 4*6 floats per odd lattice sites)
  /////////////							// now we have to consider 2 flavors: up, dn
  
  int N;
#ifdef _USE_MPI
  N = (VOLUME+RAND)/2;
#else
  N = VOLUME/2;
#endif  
  size_t dev_spinsize_d = 12*N * sizeof(dev_spinor_d); /* double2 */
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
    if(g_cart_id==0){
      printf("Error in init_doublesolve_eo_nd(): Memory allocation of nd additional double spinor fields failed.\n");
      printf("Error number %d. Aborting...\n", cudaerr); 
    }  
    exit(200);
  }

  
  // debug	// CUDA
  #ifdef CUDA_DEBUG 
    #ifndef _USE_MPI
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
  
  grid[5] = N;

  //done in init_mixedsolve_eo
  //cudaMalloc((void **) &dev_grid, 6*sizeof(int));				// dev_grid
  cudaMemcpy(dev_grid, &(grid[0]), 6*sizeof(int), cudaMemcpyHostToDevice);
  
  // debug	// CUDA
  #ifdef CUDA_DEBUG
    #ifndef _USE_MPI
      CUDA_CHECK("CUDA error in init_doublesolve_eo_nd(). Memory allocation of grid[] specifications failed.", "Allocated grid[] specifications on device.");
    #else
      CUDA_CHECK("CUDA error in init_doublesolve_eo_nd(). Memory allocation of grid[] specifications failed.", "Allocated grid[] specifications on devices.");
    #endif
  #endif
  

}//init_doublesolve_eo_nd()



void finalize_doublesolve_eo_nd(void) {
  
   
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
    cudaError_t cudaerr; 
    CUDA_CHECK("CUDA error in finalize_doublesolve_eo_nd(). Device memory deallocation failed", "Device memory deallocated.");
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
	cublasStatus cublasstatus;  
	CUBLAS_HELPER_CHECK(cublasInit(), "CUBLAS error in cublasInit(). Couldn't initialize CUBLAS.", "CUBLAS initialized.");
	if(g_debug_level > 3) printf("cublasstatus = %f\n", cublasstatus);
      #else
	#ifdef CUDA_45
	  cublasHandle_t handle;
	  cublasCreate(&handle);
	#else
	  cublasInit();
	#endif 
      #endif


  
  
  
  ///////////////
  // ALGORITHM //
  ///////////////
  
  
  // initialize x(0) = 0	// will be added up
  dev_zero_spinor_field_d<<<gpu_gd_linalg_d, gpu_bd_linalg_d>>>(x_up);
  dev_zero_spinor_field_d<<<gpu_gd_linalg_d, gpu_bd_linalg_d>>>(x_dn);
  
  // r(0) = b - A*x(0) = b
  dev_copy_spinor_field_d<<<gpu_gd_linalg_d, gpu_bd_linalg_d>>>(Q_up, r_up);
  dev_copy_spinor_field_d<<<gpu_gd_linalg_d, gpu_bd_linalg_d>>>(Q_dn, r_dn);
  
  // d(0) = r(0)
  dev_copy_spinor_field_d<<<gpu_gd_linalg_d, gpu_bd_linalg_d>>>(r_up, d_up);
  dev_copy_spinor_field_d<<<gpu_gd_linalg_d, gpu_bd_linalg_d>>>(r_dn, d_dn);
  
      // debug	// kernel
      #ifdef CUDA_DEBUG
	CUDA_KERNEL_CHECK("Kernel error in dev_cg_eo_nd(). Initializing spinor fields on device failed.", "Spinor fields initialized on device.");
      #endif
  
  
  
  // rr = (r_up)^2 + (r_dn)^2

  #ifdef _USE_MPI
    //launch error and exit
    if(g_cart_id == 0){
      printf("Error: pure double GPU ND solver not implemented for MPI. Aborting...\n");
      exit(200);
    }
  #endif
  
    rr_up = cublasDdot_wrapper(N_floats_int, (double *) r_up, (double *) r_up);
    rr_dn = cublasDdot_wrapper(N_floats_int, (double *) r_dn, (double *) r_dn);

  rr    = rr_up + rr_dn;

  r0r0   = rr;	// for relative precision 
  rr_old = rr;	// for the first iteration
  
  //////////
  // LOOP //
  //////////
  
      #ifndef LOWOUTPUT
      if (g_cart_id == 0) printf("\nEntering inner loop.\n");
      #endif

      // debug	// CUBLAS core function
      #ifdef CUDA_DEBUG
	// CUBLAS_CORE_CHECK("CUBLAS error in cg_eo_nd(). Calculating initial residue failed.", "Initial inner residue calculated.");
	CUBLAS_CORE_CHECK_NO_SUCCESS_MSG("CUBLAS error in dev_cg_eo_nd_d(). Calculating initial residue failed.");
      #endif

      if (g_cart_id == 0) printf("Initial inner residue: %.6e\n", r0r0);
  
  
  for (j = 0; j < max_iter; j++) {

     
      // A*d(k)

    dev_Qtm_pm_ndpsi_d(Ad_up, Ad_dn, d_up,  d_dn,										// debugging: matrix_debug1(), matrix_multiplication_test()
		      gpu_gd_M_d, gpu_bd_M_d, gpu_gd_linalg_d, gpu_bd_linalg_d,
		      gpu_gd_linalg_d, gpu_bd_linalg_d, gpu_gd_linalg_d, gpu_bd_linalg_d);

	if(shift != 0.0f){
           //add constant shift if nonzero
           // CUBLAS:
          cublasDaxpy_wrapper (N_floats_int, shift, (double *) d_up, (double *) Ad_up);
          cublasDaxpy_wrapper (N_floats_int, shift, (double *) d_dn, (double *) Ad_dn);
        }      
        if((cudaerr=cudaGetLastError()) != cudaSuccess){
           printf("%s\n", cudaGetErrorString(cudaerr));
           exit(200);
        }
      
	#ifdef CUDA_DEBUG
	  CUDA_CHECK("CUDA error in matrix_muliplication32(). Applying the matrix on GPU failed.", "The matrix was applied on GPU.");
	#endif

    // alpha = r(k)*r(k) / d(k)*A*d(k)
      dAd_up = cublasDdot_wrapper(N_floats_int, (double *) d_up, (double *) Ad_up);
      dAd_dn = cublasDdot_wrapper(N_floats_int, (double *) d_dn, (double *) Ad_dn);

    dAd    = dAd_up + dAd_dn;
    
    		// debug	// is NaN ?
    		if isnan(dAd) {
    		  printf("Error in dev_cg_eo_nd_d(). dAd is NaN.\n");
    		  exit(-1);
    		}
    
    alpha  = rr_old / dAd;	// rr_old is taken from the last iteration respectively
    
    // x(k+1) = x(k) + alpha*d(k)
    cublasDaxpy_wrapper(N_floats_int, alpha, (double *) d_up, (double *) x_up);
    cublasDaxpy_wrapper(N_floats_int, alpha, (double *) d_dn, (double *) x_dn);
    
    // r(k+1)
    if ( (j+1) % N_recalc_res != 0 ) {	// r(k+1) = r(k) - alpha*A*d(k)
      cublasDaxpy_wrapper(N_floats_int, -1.0*alpha, (double *) Ad_up, (double *) r_up);
      cublasDaxpy_wrapper(N_floats_int, -1.0*alpha, (double *) Ad_dn, (double *) r_dn);
    }
    else {				// recalculate residue r(k+1) = b - A*x(k+1)
	// debug
	#ifndef _USE_MPI
	  printf("Recalculating the inner residue.\n");
	#else
	  if (g_cart_id == 0) printf("Recalculating the inner residue.\n");
	#endif
      
      
      // A*x(k+1)


        dev_Qtm_pm_ndpsi_d(Ax_up, Ax_dn, x_up,  x_dn, 
        	           gpu_gd_M_d, gpu_bd_M_d, gpu_gd_linalg_d, gpu_bd_linalg_d,
        	           gpu_gd_linalg_d, gpu_bd_linalg_d, gpu_gd_linalg_d, gpu_bd_linalg_d);

	if(shift != 0.0){
           //add constant shift if nonzero
           // CUBLAS:
          cublasDaxpy_wrapper (N_floats_int, shift, (double *) x_up, (double *) Ad_up);
          cublasDaxpy_wrapper (N_floats_int, shift, (double *) x_dn, (double *) Ad_dn);
        }    

      // r(k+1) = b - A*x(k+1)
      //cublasDcopy(N_floats_int, (double *) Q_up, 1, (double *) r_up, 1);		
      //cublasDcopy(N_floats_int, (double *) Q_dn, 1, (double *) r_dn, 1);		
      dev_copy_spinor_field_d<<<gpu_gd_linalg_d, gpu_bd_linalg_d>>>(Q_up, r_up);        // r_up = Q_up
      dev_copy_spinor_field_d<<<gpu_gd_linalg_d, gpu_bd_linalg_d>>>(Q_dn, r_dn);        // r_dn = Q_dn
      cublasDaxpy_wrapper(N_floats_int, -1.0, (double *) Ax_up, (double *) r_up);	// r_up = Q_up - Ax_up
      cublasDaxpy_wrapper(N_floats_int, -1.0, (double *) Ax_dn, (double *) r_dn);	// r_dn = Q_dn - Ax_dn
    
    } // recalculate residue
    

    // r(k+1)*r(k+1)
      rr_up  = cublasDdot_wrapper(N_floats_int, (double *) r_up, (double *) r_up);
      rr_dn  = cublasDdot_wrapper(N_floats_int, (double *) r_dn, (double *) r_dn);

    rr     = rr_up + rr_dn;
    
    #ifdef CUDA_DEBUG
      CUBLAS_CORE_CHECK_NO_SUCCESS_MSG("CUBLAS error in dev_cg_eo_nd(). CUBLAS function failed.");
    #endif

    #ifndef LOWOUTPUT
      if (g_cart_id == 0) printf("inner iteration j = %i: rr = %.6e\n", j, rr);
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
    cublasDscal_wrapper (N_floats_int, beta, (double *) d_up);
    cublasDaxpy_wrapper (N_floats_int, 1.0 , (double *) r_up, (double *) d_up);
    
    cublasDscal_wrapper (N_floats_int, beta, (double *) d_dn);
    cublasDaxpy_wrapper (N_floats_int, 1.0 , (double *) r_dn, (double *) d_dn);
    
    // debug	// CUBLAS core function
    #ifdef CUDA_DEBUG
      CUBLAS_CORE_CHECK_NO_SUCCESS_MSG("CUBLAS error in dev_cg_eo_nd(). Error in CUBLAS function.");
    #endif
  }//LOOP
  
  
  #ifndef LOWOUTPUT
  if (g_cart_id == 0) printf("Finished inner loop because of maximal number of inner iterations.\n");
  #endif
  if (g_cart_id == 0) printf("Final inner residue: %.6e\n", rr);

    
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
      
     size_t dev_spinsize_d;
#ifdef _USE_MPI
   dev_spinsize_d = 12*(VOLUME+RAND)/2 * sizeof(dev_spinor_d); // double2 even-odd ! 
#else
   dev_spinsize_d = 12*VOLUME/2 * sizeof(dev_spinor_d);     
#endif
     set_global_sizes();
     update_constants_d(dev_grid);
     update_gpu_gf_d(g_gauge_field);    
     update_bare_constants_nd();
  
 
    // check mubar and epsbar etc. on host and device
    #ifdef STUFF_DEBUG
      check_nd_mixedsolve_params();
    #endif
  


   #ifdef MATRIX_DEBUG
    test_double_nd_operator(Q_up, Q_dn, N_sites_int);
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
  

  if( (g_cart_id == 0) && (g_debug_level > 2) ) printf("phmc_invmaxev = %f\n",phmc_invmaxev);  
  
  
    // r(0)
    // r(0) = b - A*x(0) = Q - A*P
      bb = square_norm(P_up, VOLUME/2, 1);
      bb += square_norm(P_dn, VOLUME/2, 1);
      if( (g_cart_id == 0) && (g_debug_level > 2) ) printf("Norm of initial guess: %.16e \n", bb);
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
		  gpu_gd_M_d, gpu_bd_M_d, gpu_gd_linalg_d, gpu_bd_linalg_d,
		  gpu_gd_linalg_d, gpu_bd_linalg_d, gpu_gd_linalg_d, gpu_bd_linalg_d);
      if(shift != 0.0) {
	dev_axpy_d<<<gpu_gd_linalg_d,gpu_bd_linalg_d>>>(shift, x_up_d, Ax_up_d);
	dev_axpy_d<<<gpu_gd_linalg_d,gpu_bd_linalg_d>>>(shift, x_dn_d, Ax_dn_d);
      }        
    
      dev_diff_d<<<gpu_gd_linalg_d,gpu_bd_linalg_d>>>(r_up_d,Q_up_d,Ax_up_d);         
      dev_diff_d<<<gpu_gd_linalg_d,gpu_bd_linalg_d>>>(r_dn_d,Q_dn_d,Ax_dn_d); 
      
      cudaMemcpy(h2d_spin_d, r_up_d , dev_spinsize_d, cudaMemcpyDeviceToHost);      
      unorder_spin_gpu(h2d_spin_d, r_up);
      cudaMemcpy(h2d_spin_d, r_dn_d , dev_spinsize_d, cudaMemcpyDeviceToHost);      
      unorder_spin_gpu(h2d_spin_d, r_dn);
     
    
    // rr = (r_up)^2 + (r_dn)^2
    rr_up = double_dotprod(r_up_d, r_up_d, N_sites_int);
    rr_dn = double_dotprod(r_dn_d, r_dn_d, N_sites_int);
    rr = rr_up + rr_dn;
    rr_old = rr; // for the first iteration  
    
    r0r0   = double_dotprod(Q_up_d, Q_up_d, N_sites_int) 
           + double_dotprod(Q_dn_d, Q_dn_d, N_sites_int);    


  if ( (g_cart_id == 0) && (g_debug_level > 1) ) printf("Initial outer residue: %.10e\n", rr_old);


  ////////////////
  // OUTER LOOP //
  ////////////////

  do {		// for (i = 0; i < max_iter; i++) {
    
    i++;
  
    #ifndef LOWOUTPUT
    if ( (g_cart_id == 0) && (g_debug_level > 1) ) printf("Outer iteration i = %i\n", i);
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
      if ( (g_cart_id == 0) && (g_debug_level > 1) ){
        effectiveflops = innercount*(matrixflops + 2*2*2*24 + 2*2*24 + 2*2*24 + 2*2*2*24 + 2*2*24)*VOLUME/2;   
        printf("inner solver BENCHMARK:\n");
        printf("\tinner solver performance:  %.4e Gflop/s\n", double(effectiveflops)/innerclocks/ 1.0e9);
      }
    #endif
      
    //copy result back
    cudaMemcpy(h2d_spin_d, dev_spinout_up_d , dev_spinsize_d, cudaMemcpyDeviceToHost);      
    unorder_spin_gpu(h2d_spin_d, d_up);
    cudaMemcpy(h2d_spin_d, dev_spinout_dn_d , dev_spinsize_d, cudaMemcpyDeviceToHost);      
    unorder_spin_gpu(h2d_spin_d, d_dn);
      
    		// debug	// CUDA
    		#ifdef CUDA_DEBUG
		  cudaError_t cudaerr;    
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
        if (g_cart_id == 0) printf("The matrix was applied on CPU in double precision. r = b - Ax\n");
        #endif
        diff(r_up, Q_up, Ax_up, N_sites_int);
        diff(r_dn, Q_dn, Ax_dn, N_sites_int);

      // rr = (rr_up)^2 + (r_dn)^2
      rr_up = square_norm(r_up, N_sites_int, 1);
      rr_dn = square_norm(r_dn, N_sites_int, 1);
      rr    = rr_up + rr_dn;
      

    if (g_cart_id == 0){ 
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

        
      #ifdef _USE_MPI
	if (g_cart_id == 0) {
      #endif
      printf("EO inversion done in double precision.\n");
      if (rel_prec == 0) printf("Finished outer loop because of reached absolute outer solver precision.\n");
      if (rel_prec == 1) printf("Finished outer loop because of reached relative outer solver precision.\n");
      printf("Total number of inner iterations: %i\n", outercount);
      printf("Total number of outer iterations: %i\n", i+1);
      printf("Squared residue: %.10e\n\n", rr); 
      #ifdef _USE_MPI
	}
      #endif
		  

      finalize_solver(up_field, nr_sf);
      finalize_solver(dn_field, nr_sf);  
      
      return(outercount);  
    }
    
  
  }//OUTER LOOP
  while (outercount <= max_iter);
  
  
  // multiplying with Qhat(2x2)^dagger is done in invert_doublet_eo.c
    #ifdef _USE_MPI
      if (g_cart_id == 0) {
    #endif
    printf("EO inversion done in full precision.\n");
    printf("Finished outer loop, because of maximal number of outer iterations.\n");
    printf("Total number of inner iterations: %i\n", outercount);
    printf("Total number of outer iterations: %i\n", i+1);
    printf("Squared residue: %.10e\n\n", rr); 
    #ifdef _USE_MPI
      }
    #endif
    


  finalize_solver(up_field, nr_sf);
  finalize_solver(dn_field, nr_sf);
  return(outercount);
  
  
}//doublesolve_eo_nd()


dev_spinor_d * _mms_d_up_d_d;
dev_spinor_d * _mms_d_dn_d_d;
dev_spinor_d * _mms_x_up_d_d;
dev_spinor_d * _mms_x_dn_d_d;

dev_spinor_d ** mms_d_up_d;
dev_spinor_d ** mms_d_dn_d;
dev_spinor_d ** mms_x_up_d;
dev_spinor_d ** mms_x_dn_d;



void init_gpu_nd_mms_fields(int Nshift){
  
  cudaError_t cudaerr;

  int N;
#ifdef _USE_MPI
  N = (VOLUME+RAND)/2;
#else
  N = VOLUME/2;
#endif
  
  //here we allocate one spinor pair less than the number of shifts
  //as for the zero'th shift we re-use fields from the usual cg solver
   size_t dev_spinsize_d = (Nshift-1)*12*N * sizeof(dev_spinor_d); /* double2 */  
   cudaMalloc((void **) &_mms_d_up_d_d, dev_spinsize_d);
   cudaMalloc((void **) &_mms_d_dn_d_d, dev_spinsize_d);   
   cudaMalloc((void **) &_mms_x_up_d_d, dev_spinsize_d);
   cudaMalloc((void **) &_mms_x_dn_d_d, dev_spinsize_d);  


  mms_d_up_d = (dev_spinor_d**)malloc(((Nshift-1)*sizeof(dev_spinor_d*)));
  mms_d_dn_d = (dev_spinor_d**)malloc(((Nshift-1)*sizeof(dev_spinor_d*)));
  mms_x_up_d = (dev_spinor_d**)malloc(((Nshift-1)*sizeof(dev_spinor_d*)));
  mms_x_dn_d = (dev_spinor_d**)malloc(((Nshift-1)*sizeof(dev_spinor_d*)));  
  
  mms_d_up_d[0] = _mms_d_up_d_d;
  mms_d_dn_d[0] = _mms_d_dn_d_d;
  mms_x_up_d[0] = _mms_x_up_d_d;
  mms_x_dn_d[0] = _mms_x_dn_d_d;  
  for(int im = 1; im < (Nshift-1); im++) {
    mms_d_up_d[im] = _mms_d_up_d_d + im*12*N;
    mms_d_dn_d[im] = _mms_d_dn_d_d + im*12*N; 
    mms_x_up_d[im] = _mms_x_up_d_d + im*12*N;
    mms_x_dn_d[im] = _mms_x_dn_d_d + im*12*N;     
  }

#ifdef _USE_MPI
  int tSliceEO = LX*LY*LZ/2;
  R1_UP_D = (dev_spinor_d *) malloc(2*tSliceEO*24*sizeof(double));
  R2_UP_D = R1_UP_D + 12*tSliceEO;
  R3_UP_D = (dev_spinor_d *) malloc(2*tSliceEO*24*sizeof(double));
  R4_UP_D = R3_UP_D + 12*tSliceEO;


//for gathering and spreading of indizes of rand in (gather_rand spread_rand called from xchange_field_wrapper)
    #ifdef RELATIVISTIC_BASIS
      cudaMalloc((void **) &RAND_FW_UP_D, tSliceEO*6*sizeof(double2));
      cudaMalloc((void **) &RAND_BW_UP_D, tSliceEO*6*sizeof(double2));
    #else
      cudaMalloc((void **) &RAND_FW_UP_D, tSliceEO*12*sizeof(double2));
      cudaMalloc((void **) &RAND_BW_UP_D, tSliceEO*12*sizeof(double2));      
    #endif
    /*  for async communication */
    // page-locked memory    
    #ifdef RELATIVISTIC_BASIS
      int dbperspin = 6;
    #else
      int dbperspin = 12;
    #endif  
    cudaMallocHost(&RAND3_UP_D, tSliceEO*dbperspin*sizeof(double2));
    cudaMallocHost(&RAND4_UP_D, tSliceEO*dbperspin*sizeof(double2));
    cudaMallocHost(&RAND1_UP_D, tSliceEO*dbperspin*sizeof(double2));
    cudaMallocHost(&RAND2_UP_D, tSliceEO*dbperspin*sizeof(double2));

  R1_UP_D = (dev_spinor_d *) malloc(2*tSliceEO*24*sizeof(double));
  R2_UP_D = R1_UP_D + 12*tSliceEO;
  R3_UP_D = (dev_spinor_d *) malloc(2*tSliceEO*24*sizeof(double));
  R4_UP_D = R3_UP_D + 12*tSliceEO;


//for gathering and spreading of indizes of rand in (gather_rand spread_rand called from xchange_field_wrapper)
    #ifdef RELATIVISTIC_BASIS
      cudaMalloc((void **) &RAND_FW_DN_D, tSliceEO*6*sizeof(double2));
      cudaMalloc((void **) &RAND_BW_DN_D, tSliceEO*6*sizeof(double2));
    #else
      cudaMalloc((void **) &RAND_FW_DN_D, tSliceEO*12*sizeof(double2));
      cudaMalloc((void **) &RAND_BW_DN_D, tSliceEO*12*sizeof(double2));      
    #endif

    cudaMallocHost(&RAND3_DN_D, tSliceEO*dbperspin*sizeof(double2));
    cudaMallocHost(&RAND4_DN_D, tSliceEO*dbperspin*sizeof(double2));
    cudaMallocHost(&RAND1_DN_D, tSliceEO*dbperspin*sizeof(double2));
    cudaMallocHost(&RAND2_DN_D, tSliceEO*dbperspin*sizeof(double2));
    
  R1_DN_D = (dev_spinor_d *) malloc(2*tSliceEO*24*sizeof(double));
  R2_DN_D = R1_DN_D + 12*tSliceEO;
  R3_DN_D = (dev_spinor_d *) malloc(2*tSliceEO*24*sizeof(double));
  R4_DN_D = R3_DN_D + 12*tSliceEO;    
#endif
  
  if((cudaerr=cudaGetLastError())!=cudaSuccess){
    if(g_cart_id==0){
      printf("Error in init_gpu_nd_mms_fields: Memory allocation of spinor fields failed.\n");
      printf("Error number %d. Aborting...", cudaerr);
      exit(200);
    }
  }
  else{
    if(g_cart_id==0) printf("Allocated nd mms double spinor fields on device\n");
  } 
  
}



void finalize_gpu_nd_mms_fields(){
  free(mms_d_up_d);
  free(mms_d_dn_d); 
  free(mms_x_up_d);
  free(mms_x_dn_d);  
  cudaFree(_mms_d_up_d_d);
  cudaFree(_mms_d_dn_d_d);
  cudaFree(_mms_x_up_d_d);
  cudaFree(_mms_x_dn_d_d);  
#ifdef _USE_MPI
  cudaFreeHost(RAND1_UP_D);
  cudaFreeHost(RAND2_UP_D); 
  cudaFreeHost(RAND3_UP_D);
  cudaFreeHost(RAND4_UP_D);  
  cudaFree(RAND_BW_UP_D);
  cudaFree(RAND_FW_UP_D);
  free(R1_UP_D);
  free(R3_UP_D);
  cudaFreeHost(RAND1_DN_D);
  cudaFreeHost(RAND2_DN_D); 
  cudaFreeHost(RAND3_DN_D);
  cudaFreeHost(RAND4_DN_D);  
  cudaFree(RAND_BW_DN_D);
  cudaFree(RAND_FW_DN_D);
  free(R1_DN_D);
  free(R3_DN_D); 
#endif
}


// a pure double mms solver
// arithmetics equal the cpu counterpart in solver/cg_mms_tm_nd.c
// calculation of alphas, betas, zitas ... taken over from there
int dev_cg_mms_eo_nd_d (dev_su3_2v_d * gf,
              dev_spinor_d * P_up, dev_spinor_d * P_dn,
              dev_spinor_d * Q_up, dev_spinor_d * Q_dn,
	      double * shifts, int Nshift,
              int max_iter,
              int check_abs , int check_rel,
              double eps_abs, double eps_rel, int min_solver_it       ) {
  
  // P_up/dn  can be used as auxiliary field to work on, as it is not later used (could be used as initial guess at the very start)
  // Q_up/dn  can be used as feedback, or if not, also as auxiliary field
  
  int noshifts = Nshift;
  
  //allocate cg constants
  double * sigma;
  double * zitam1, * zita;
  double * alphas, * betas;
  double gamma;
  double alpham1;
    sigma = (double*)calloc((noshifts), sizeof(double));
    zitam1 = (double*)calloc((noshifts), sizeof(double));
    zita = (double*)calloc((noshifts), sizeof(double));
    alphas = (double*)calloc((noshifts), sizeof(double));
    betas = (double*)calloc((noshifts), sizeof(double));

  /////////////////////
  // LOCAL VARIABLES //
  /////////////////////
  
  
  // algorithm
  double rr_up, rr_dn, rr, rr_old, r0r0, dAd_up, dAd_dn, dAd;


  // (auxiliary) device fields
  // for recalculating the residue
  dev_spinor_d *  r_up, *  r_dn, * Ad_up, * Ad_dn, *  x_up, *  x_dn, *  d_up, *  d_dn;		
 
  // counting
  int j;				// iteration counter
 

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

  
  
 
  /////////////////////
  // INITIALIZATIONS //
  /////////////////////
  
      // debug	// CUBLAS helper function
      #ifdef CUDA_DEBUG
        cudaError_t cudaerr;
        cublasStatus cublasstatus;
	CUBLAS_HELPER_CHECK(cublasInit(), "CUBLAS error in cublasInit(). Couldn't initialize CUBLAS.", "CUBLAS initialized.");
	if(g_debug_level > 3) printf("cublasstatus = %f\n", cublasstatus);
      #else
	#ifdef CUDA_45
	  cublasHandle_t handle;
	  cublasCreate(&handle);
	#else
	  cublasInit();
	#endif 
      #endif


  
  
  
  ///////////////
  // ALGORITHM //
  ///////////////
  
  
  // initialize x(0) = 0	// will be added up
  dev_zero_spinor_field_d<<<gpu_gd_blas_d, gpu_bd_blas_d>>>(x_up);
  dev_zero_spinor_field_d<<<gpu_gd_blas_d, gpu_bd_blas_d>>>(x_dn);
  
  // r(0) = b - A*x(0) = b
  dev_copy_spinor_field_d<<<gpu_gd_blas_d, gpu_bd_blas_d>>>(Q_up, r_up);
  dev_copy_spinor_field_d<<<gpu_gd_blas_d, gpu_bd_blas_d>>>(Q_dn, r_dn);
  
  // d(0) = r(0)
  dev_copy_spinor_field_d<<<gpu_gd_blas_d, gpu_bd_blas_d>>>(r_up, d_up);
  dev_copy_spinor_field_d<<<gpu_gd_blas_d, gpu_bd_blas_d>>>(r_dn, d_dn);
  
      // debug	// kernel
      #ifdef CUDA_DEBUG
  
	CUDA_KERNEL_CHECK("Kernel error in dev_cg_mms_eo_nd_d(). Initializing spinor fields on device failed.", "Spinor fields initialized on device.");
      #endif
  
  
  
  // rr = (r_up)^2 + (r_dn)^2
  
    rr_up = cublasDdot_wrapper(N_floats_int, (double *) r_up, (double *) r_up);
    rr_dn = cublasDdot_wrapper(N_floats_int, (double *) r_dn, (double *) r_dn);

  rr    = rr_up + rr_dn;

  r0r0   = rr;	// for relative precision 
  rr_old = rr;	// for the first iteration
  
      // debug	// CUBLAS core function
      #ifdef CUDA_DEBUG
	// CUBLAS_CORE_CHECK("CUBLAS error in cg_eo_nd(). Calculating initial residue failed.", "Initial inner residue calculated.");
	CUBLAS_CORE_CHECK_NO_SUCCESS_MSG("CUBLAS error in dev_cg_mms_eo_nd_d(). Calculating initial residue failed.");
      #endif

      if (g_cart_id == 0) printf("Initial mms residue: %.6e\n", r0r0);
  
  alphas[0] = 1.0;
  betas[0] = 0.0;
  sigma[0] = shifts[0]*shifts[0];
  if(g_cart_id == 0 && g_debug_level > 2) printf("# dev_CGMMSND: shift %d is %e\n", 0, sigma[0]);

  /* currently only implemented for P=0 */
  for(int im = 1; im < noshifts; im++) {
    sigma[im] = shifts[im]*shifts[im] - sigma[0];
    if(g_cart_id == 0 && g_debug_level > 2) printf("# dev_CGMMSND: shift %d is %e\n", im, sigma[im]);
    // these will be the result spinor fields
    dev_zero_spinor_field_d<<<gpu_gd_blas_d, gpu_bd_blas_d>>>(mms_x_up_d[im-1]);
    dev_zero_spinor_field_d<<<gpu_gd_blas_d, gpu_bd_blas_d>>>(mms_x_dn_d[im-1]);    

    // these are intermediate fields
    dev_copy_spinor_field_d<<<gpu_gd_blas_d, gpu_bd_blas_d>>>(Q_up, mms_d_up_d[im-1]);
    dev_copy_spinor_field_d<<<gpu_gd_blas_d, gpu_bd_blas_d>>>(Q_dn, mms_d_dn_d[im-1]);
    zitam1[im] = 1.0;
    zita[im] = 1.0;
    alphas[im] = 1.0;
    betas[im] = 0.0;
  }


  //////////
  // LOOP //
  //////////
    
  for (j = 0; j < max_iter; j++) {

      // A*d(k)
    dev_Qtm_pm_ndpsi_d(Ad_up, Ad_dn, d_up,  d_dn,	
		      gpu_gd_M_d, gpu_bd_M_d, gpu_gd_linalg_d, gpu_bd_linalg_d,
		      gpu_gd_linalg_d, gpu_bd_linalg_d, gpu_gd_linalg_d, gpu_bd_linalg_d);

    //add zero'th shift
    cublasDaxpy_wrapper (N_floats_int, sigma[0], (double *) d_up, (double *) Ad_up);
    cublasDaxpy_wrapper (N_floats_int, sigma[0], (double *) d_dn, (double *) Ad_dn);
	

     #ifdef CUDA_DEBUG
	  CUDA_CHECK("CUDA error in matrix_muliplication(). Applying the matrix on GPU failed.", "The matrix was applied on GPU.");
     #endif

    // alpha = r(k)*r(k) / d(k)*A*d(k)
      dAd_up = cublasDdot_wrapper(N_floats_int, (double *) d_up, (double *) Ad_up);
      dAd_dn = cublasDdot_wrapper(N_floats_int, (double *) d_dn, (double *) Ad_dn);

    dAd    = dAd_up + dAd_dn; 
    alpham1 = alphas[0];
    alphas[0]  = rr_old / dAd;	// rr_old is taken from the last iteration respectively
    
    
    // x_j(k+1) = x_j(k) + alpha_j*d_j(k)   
    for(int im = 1; im < noshifts; im++) {
      gamma = zita[im]*alpham1/(alphas[0]*betas[0]*(1.-zita[im]/zitam1[im]) 
				+ alpham1*(1.+sigma[im]*alphas[0]));
      zitam1[im] = zita[im];
      zita[im] = gamma;
      alphas[im] = alphas[0]*zita[im]/zitam1[im];
      cublasDaxpy_wrapper(N_floats_int, alphas[im], (double *) mms_d_up_d[im-1], (double *) mms_x_up_d[im-1]);
      cublasDaxpy_wrapper(N_floats_int, alphas[im], (double *) mms_d_dn_d[im-1], (double *) mms_x_dn_d[im-1]);
      if(j > 0 && (j % 10 == 0) && (im == noshifts-1)) {
	double sn = cublasDdot_wrapper(N_floats_int, (double *) mms_d_up_d[im-1], (double *) mms_d_up_d[im-1]);
	sn += cublasDdot_wrapper(N_floats_int, (double *) mms_d_dn_d[im-1], (double *) mms_d_dn_d[im-1]);
	if(alphas[noshifts-1]*alphas[noshifts-1]*sn <= eps_abs) {
	  noshifts--;
	  if(g_debug_level > 2 && g_cart_id == 0) {
	    printf("# dev_CGMMSND: at iteration %d removed one shift, %d remaining\n", j, noshifts);
	  }
	}
      }
    }   
    // x_0(k+1) = x_0(k) + alpha_0*d_0(k)      
    cublasDaxpy_wrapper(N_floats_int, alphas[0], (double *) d_up, (double *) x_up);
    cublasDaxpy_wrapper(N_floats_int, alphas[0], (double *) d_dn, (double *) x_dn);
    
    // r(k+1)
    cublasDaxpy_wrapper(N_floats_int, -1.0*alphas[0], (double *) Ad_up, (double *) r_up);
    cublasDaxpy_wrapper(N_floats_int, -1.0*alphas[0], (double *) Ad_dn, (double *) r_dn);

    // r(k+1)*r(k+1)
      rr_up  = cublasDdot_wrapper(N_floats_int, (double *) r_up, (double *) r_up);
      rr_dn  = cublasDdot_wrapper(N_floats_int, (double *) r_dn, (double *) r_dn);

    rr     = rr_up + rr_dn;
    
    #ifdef CUDA_DEBUG
      CUBLAS_CORE_CHECK_NO_SUCCESS_MSG("CUBLAS error in dev_cg_eo_nd(). CUBLAS function failed.");
    #endif

    #ifndef LOWOUTPUT
      if (g_cart_id == 0) printf("mms iteration j = %i: rr = %.6e\n", j, rr);
    #endif
		 
	
    // debug	// is NaN ?
    if isnan(rr) {
      printf("Error in dev_cg_mms_eo_nd_d(). Inner residue is NaN.\n");
      exit(-1);
    }
    
    // aborting ?? // check wether precision is reached ...
    if ( ((check_abs)&&(rr <= eps_abs)&&(j>min_solver_it)) || ((check_rel)&&(rr <= eps_rel*r0r0)&&(j>min_solver_it)) ) {
   
      if ((check_rel)&&(rr <= eps_rel*r0r0)) {
	printf("Reached relative solver precision of eps_rel = %.2e\n", eps_rel);
      }
      if ((check_abs)&&(rr <= eps_abs)) {
	printf("Reached absolute solver precision of eps_abs = %.2e\n", eps_abs);
      }
      
      printf("Final mms solver residue: %.6e\n", rr);
      
	#ifdef CUDA_45  
	  cublasDestroy(handle);
	#else
	  cublasShutdown();
	#endif 
	
	//free cg constants
	free(sigma);
	free(zitam1);
	free(zita);
	free(alphas);
	free(betas);    
	return(j);
    }
    

    // beta = r(k+1)*r(k+1) / r(k)*r(k)
    betas[0] = rr / rr_old;
    rr_old = rr;  // for next iteration
    
    
    // d(k+1) = r(k+1) + beta*d(k)
    cublasDscal_wrapper (N_floats_int, betas[0], (double *) d_up);
    cublasDaxpy_wrapper (N_floats_int, 1.0 , (double *) r_up, (double *) d_up);
    
    cublasDscal_wrapper (N_floats_int, betas[0], (double *) d_dn);
    cublasDaxpy_wrapper (N_floats_int, 1.0 , (double *) r_dn, (double *) d_dn);
    
     for(int im = 1; im < noshifts; im++) {
      betas[im] = betas[0]*zita[im]*alphas[im]/(zitam1[im]*alphas[0]);
      cublasDscal_wrapper (N_floats_int, betas[im], (double *) mms_d_up_d[im-1]);
      cublasDaxpy_wrapper (N_floats_int, zita[im] , (double *) r_up, (double *) mms_d_up_d[im-1]);
    
      cublasDscal_wrapper (N_floats_int, betas[im], (double *) mms_d_dn_d[im-1]);
      cublasDaxpy_wrapper (N_floats_int, zita[im] , (double *) r_dn, (double *) mms_d_dn_d[im-1]);
     }   
    
    // debug	// CUBLAS core function
    #ifdef CUDA_DEBUG
      CUBLAS_CORE_CHECK_NO_SUCCESS_MSG("CUBLAS error in dev_cg_eo_nd(). Error in CUBLAS function.");
    #endif
  }//LOOP
  
  
  #ifndef LOWOUTPUT
  if (g_cart_id == 0) printf("Finished because of maximal number of iterations.\n");
  #endif
  if (g_cart_id == 0) printf("Final mms residue: %.6e\n", rr);

    
    #ifdef CUDA_45  
      cublasDestroy(handle);
    #else
      cublasShutdown();
    #endif 
 
      
  //free cg constants
  free(sigma);
  free(zitam1);
  free(zita);
  free(alphas);
  free(betas);    
      
  return(j);
  
}//dev_cg_mms_eo_nd_d()




extern "C" int doublesolve_mms_eo_nd (spinor ** P_up, spinor ** P_dn,
                                 spinor * Q_up, spinor * Q_dn, double * shifts, int Nshift,
                                 int max_iter, double eps_sq, int rel_prec, int min_solver_it) {
   
  
  // CUDA
  cudaError_t cudaerr;
  
  // counting			// latest inner solver iterations
  int outercount = 0;				// total inner solver iterations

  #ifdef ALGORITHM_BENCHMARK
    double effectiveflops;
    //double hoppingflops = 1608.0;
    double matrixflops  = 14448;
  #endif
  
  // timing
  double startinner, stopinner;  
  double innerclocks;
  
  if(rel_prec){
    innersolver_precision_check_rel = 1;
    innersolver_precision_check_abs = 0;
   }
   else{
    innersolver_precision_check_rel = 0;
    innersolver_precision_check_abs = 1;     
  }
  
  
  //////////////////
  // INITIALIZING //
  //////////////////
     
     size_t dev_spinsize_d;
#ifdef _USE_MPI
  dev_spinsize_d  = 12*(VOLUME+RAND)/2 * sizeof(dev_spinor_d); // double2 even-odd !   
#else
  dev_spinsize_d  = 12*VOLUME/2 * sizeof(dev_spinor_d);  
#endif

    
  update_constants_d(dev_grid);
  update_gpu_gf_d(g_gauge_field);
    

  dev_complex h0, h1, h2, h3; 
  h0.re  =  (float) creal(ka0);	h0.im  = -(float) cimag(ka0);	// ka{0-4} are defined in boundary.c
  h1.re  =  (float) creal(ka1);	h1.im  = -(float) cimag(ka1);	// what is the meaning?
  h2.re  =  (float) creal(ka2);	h2.im  = -(float) cimag(ka2);
  h3.re  =  (float) creal(ka3);	h3.im  = -(float) cimag(ka3);
  he_cg_init<<< 1, 1 >>> (dev_grid, (float) g_kappa, (float)(g_mu/(2.0*g_kappa)), h0, h1, h2, h3);

  he_cg_init_nd_additional<<<1,1>>> (g_mubar, g_epsbar);
      
  init_gpu_nd_mms_fields(Nshift);  
  set_global_sizes();
    
   #ifdef MATRIX_DEBUG
    test_double_nd_operator(Q_up, Q_dn, N_sites_int);
   #endif  
   #ifdef OPERATOR_BENCHMARK
    benchmark_eo_nd_d(Q_up, Q_dn, OPERATOR_BENCHMARK);
   #endif
    
    order_spin_gpu(Q_up, h2d_spin_d);
    cudaMemcpy(dev_spinin_up_d, h2d_spin_d, dev_spinsize_d, cudaMemcpyHostToDevice);      
    order_spin_gpu(Q_dn, h2d_spin_d);
    cudaMemcpy(dev_spinin_dn_d, h2d_spin_d, dev_spinsize_d, cudaMemcpyHostToDevice);
    #ifdef RELATIVISTIC_BASIS
      to_relativistic_basis_d<<<gpu_gd_linalg_d, gpu_bd_linalg_d>>> (dev_spinin_up_d);
      to_relativistic_basis_d<<<gpu_gd_linalg_d, gpu_bd_linalg_d>>> (dev_spinin_dn_d);
      if ((cudaerr=cudaGetLastError())!=cudaSuccess) {
	if (g_cart_id == 0) printf("%s\n", cudaGetErrorString(cudaGetLastError()));
      }
      else{
	#ifndef LOWOUTPUT 
	if (g_cart_id == 0) printf("Switched to relativistic basis\n");
	#endif
      }    
    #endif

    startinner = gettime();
    outercount = dev_cg_mms_eo_nd_d(dev_gf_d,
                          dev_spinout_up_d, dev_spinout_dn_d,
                          dev_spinin_up_d , dev_spinin_dn_d, shifts, Nshift,
                          max_iter,
                          innersolver_precision_check_abs, innersolver_precision_check_rel,
                          eps_sq, eps_sq, min_solver_it);
    stopinner = gettime();
    innerclocks = stopinner-startinner;
    #ifdef ALGORITHM_BENCHMARK
     if( (g_cart_id == 0) && (g_debug_level > 1) ){
      effectiveflops = outercount*g_nproc*(matrixflops + 2*2*2*24 + 2*2*24 + 2*2*24 + 2*2*2*24 + 2*2*24)*VOLUME/2;   
      printf("mms double solver BENCHMARK:\n");
      printf("\tsolver performance:  %.4e Gflop/s\n", double(effectiveflops)/innerclocks/ 1.0e9);
     }
    #endif
    
      
    //copy result back
    #ifdef RELATIVISTIC_BASIS 
      to_tmlqcd_basis_d<<<gpu_gd_linalg_d, gpu_bd_linalg_d>>> (dev_spinout_up_d);
      to_tmlqcd_basis_d<<<gpu_gd_linalg_d, gpu_bd_linalg_d>>> (dev_spinout_dn_d);      
    #endif     
    cudaMemcpy(h2d_spin_d, dev_spinout_up_d , dev_spinsize_d, cudaMemcpyDeviceToHost);      
    unorder_spin_gpu(h2d_spin_d, P_up[0]);
    cudaMemcpy(h2d_spin_d, dev_spinout_dn_d , dev_spinsize_d, cudaMemcpyDeviceToHost);      
    unorder_spin_gpu(h2d_spin_d, P_dn[0]);
    for(int im = 1; im < Nshift; im++) {
      #ifdef RELATIVISTIC_BASIS 
	to_tmlqcd_basis_d<<<gpu_gd_linalg_d, gpu_bd_linalg_d>>> (mms_x_up_d[im-1]);
	to_tmlqcd_basis_d<<<gpu_gd_linalg_d, gpu_bd_linalg_d>>> (mms_x_dn_d[im-1]);      
      #endif        
      cudaMemcpy(h2d_spin_d, mms_x_up_d[im-1], dev_spinsize_d, cudaMemcpyDeviceToHost);      
      unorder_spin_gpu(h2d_spin_d, P_up[im]);
      cudaMemcpy(h2d_spin_d, mms_x_dn_d[im-1] , dev_spinsize_d, cudaMemcpyDeviceToHost);      
      unorder_spin_gpu(h2d_spin_d, P_dn[im]);    
    }
  
    if((cudaerr=cudaGetLastError())!=cudaSuccess){
      if(g_cart_id==0) {
	printf("Error in doublesolve_mms_eo_nd().\n");
        printf("Error was %d. Aborting...\n", cudaerr);
      }
      exit(200);
    }    
      #ifdef _USE_MPI
	if (g_cart_id == 0) {
      #endif
      printf("EO MMS inversion done in double precision.\n");
      printf("Total number of iterations: %i\n\n", outercount);
      #ifdef _USE_MPI
	}
      #endif
		  
  finalize_gpu_nd_mms_fields(); 
  

  return(outercount);
  
}//doublesolve_mms_eo_nd()





dev_spinor * _mms_d_up;
dev_spinor * _mms_d_dn;
dev_spinor * _mms_x_up;
dev_spinor * _mms_x_dn;

dev_spinor ** mms_d_up;
dev_spinor ** mms_d_dn;
dev_spinor ** mms_x_up;
dev_spinor ** mms_x_dn;



void init_gpu_single_nd_mms_fields(int Nshift, int N){
  
  cudaError_t cudaerr;
  
  size_t dev_spinsize = (Nshift-1)*6*N*sizeof(dev_spinor); /* float4 */ 
  
  cudaMalloc((void **) &dev_spin1_up, dev_spinsize); 
  cudaMalloc((void **) &dev_spin1_dn, dev_spinsize);
  cudaMalloc((void **) &dev_spin2_up, dev_spinsize);
  cudaMalloc((void **) &dev_spin2_dn, dev_spinsize);
  cudaMalloc((void **) &dev_spin3_up, dev_spinsize);
  cudaMalloc((void **) &dev_spin3_dn, dev_spinsize);
  cudaMalloc((void **) &dev_spinin_up , dev_spinsize);
  cudaMalloc((void **) &dev_spinin_dn , dev_spinsize);
  cudaMalloc((void **) &dev_spinout_up, dev_spinsize);
  cudaMalloc((void **) &dev_spinout_dn, dev_spinsize);
  
  cudaMalloc((void **) &dev_spin_eo1_up, dev_spinsize);	
  cudaMalloc((void **) &dev_spin_eo1_dn, dev_spinsize);
  cudaMalloc((void **) &dev_spin_eo3_up, dev_spinsize);
  cudaMalloc((void **) &dev_spin_eo3_dn, dev_spinsize);  
  
  
  for (int i = 0; i < 2; i++) {
        cudaStreamCreate(&stream_nd[i]);
  }   
  
  
#ifdef MPI
  int tSliceEO = LX*LY*LZ/2;
  R1_UP_D = (dev_spinor_d *) malloc(2*tSliceEO*24*sizeof(double));
  R2_UP_D = R1_UP_D + 12*tSliceEO;
  R3_UP_D = (dev_spinor_d *) malloc(2*tSliceEO*24*sizeof(double));
  R4_UP_D = R3_UP_D + 12*tSliceEO;


//for gathering and spreading of indizes of rand in (gather_rand spread_rand called from xchange_field_wrapper)
    #ifdef RELATIVISTIC_BASIS
      cudaMalloc((void **) &RAND_FW_UP_D, tSliceEO*6*sizeof(double2));
      cudaMalloc((void **) &RAND_BW_UP_D, tSliceEO*6*sizeof(double2));
    #else
      cudaMalloc((void **) &RAND_FW_UP_D, tSliceEO*12*sizeof(double2));
      cudaMalloc((void **) &RAND_BW_UP_D, tSliceEO*12*sizeof(double2));      
    #endif
    /*  for async communication */
    // page-locked memory    
    #ifdef RELATIVISTIC_BASIS
      int dbperspin = 6;
    #else
      int dbperspin = 12;
    #endif  
    cudaMallocHost(&RAND3_UP_D, tSliceEO*dbperspin*sizeof(double2));
    cudaMallocHost(&RAND4_UP_D, tSliceEO*dbperspin*sizeof(double2));
    cudaMallocHost(&RAND1_UP_D, tSliceEO*dbperspin*sizeof(double2));
    cudaMallocHost(&RAND2_UP_D, tSliceEO*dbperspin*sizeof(double2));

  R1_UP_D = (dev_spinor_d *) malloc(2*tSliceEO*24*sizeof(double));
  R2_UP_D = R1_UP_D + 12*tSliceEO;
  R3_UP_D = (dev_spinor_d *) malloc(2*tSliceEO*24*sizeof(double));
  R4_UP_D = R3_UP_D + 12*tSliceEO;


//for gathering and spreading of indizes of rand in (gather_rand spread_rand called from xchange_field_wrapper)
    #ifdef RELATIVISTIC_BASIS
      cudaMalloc((void **) &RAND_FW_DN_D, tSliceEO*6*sizeof(double2));
      cudaMalloc((void **) &RAND_BW_DN_D, tSliceEO*6*sizeof(double2));
    #else
      cudaMalloc((void **) &RAND_FW_DN_D, tSliceEO*12*sizeof(double2));
      cudaMalloc((void **) &RAND_BW_DN_D, tSliceEO*12*sizeof(double2));      
    #endif

    cudaMallocHost(&RAND3_DN_D, tSliceEO*dbperspin*sizeof(double2));
    cudaMallocHost(&RAND4_DN_D, tSliceEO*dbperspin*sizeof(double2));
    cudaMallocHost(&RAND1_DN_D, tSliceEO*dbperspin*sizeof(double2));
    cudaMallocHost(&RAND2_DN_D, tSliceEO*dbperspin*sizeof(double2));
    
  R1_DN_D = (dev_spinor_d *) malloc(2*tSliceEO*24*sizeof(double));
  R2_DN_D = R1_DN_D + 12*tSliceEO;
  R3_DN_D = (dev_spinor_d *) malloc(2*tSliceEO*24*sizeof(double));
  R4_DN_D = R3_DN_D + 12*tSliceEO;    
#endif  
  
  //here we allocate one spinor pair less than the number of shifts
  //as for the zero'th shift we re-use fields from the usual cg solver
   cudaMalloc((void **) &_mms_d_up, dev_spinsize);
   cudaMalloc((void **) &_mms_d_dn, dev_spinsize);   
   cudaMalloc((void **) &_mms_x_up, dev_spinsize);
   cudaMalloc((void **) &_mms_x_dn, dev_spinsize);  


  mms_d_up = (dev_spinor**)malloc(((Nshift-1)*sizeof(dev_spinor*)));
  mms_d_dn = (dev_spinor**)malloc(((Nshift-1)*sizeof(dev_spinor*)));
  mms_x_up = (dev_spinor**)malloc(((Nshift-1)*sizeof(dev_spinor*)));
  mms_x_dn = (dev_spinor**)malloc(((Nshift-1)*sizeof(dev_spinor*)));  
  
  mms_d_up[0] = _mms_d_up;
  mms_d_dn[0] = _mms_d_dn;
  mms_x_up[0] = _mms_x_up;
  mms_x_dn[0] = _mms_x_dn;  
  for(int im = 1; im < (Nshift-1); im++) {
    mms_d_up[im] = _mms_d_up + im*6*N;
    mms_d_dn[im] = _mms_d_dn + im*6*N; 
    mms_x_up[im] = _mms_x_up + im*6*N;
    mms_x_dn[im] = _mms_x_dn + im*6*N;     
  }
  
  if((cudaerr=cudaGetLastError())!=cudaSuccess){
    if(g_cart_id==0){ 
      printf("Error in init_gpu_nd_mms_fields: Memory allocation of spinor fields failed.\n");
      printf("Error was %d.  Aborting...", cudaerr);
    }
    exit(200);
  }
  else{
    if((g_cart_id==0) && (g_debug_level > 2)) printf("Allocated nd mms single spinor fields on device\n");
  } 
  
}



void finalize_gpu_single_nd_mms_fields(){
  free(mms_d_up);
  free(mms_d_dn); 
  free(mms_x_up);
  free(mms_x_dn);  
  cudaFree(_mms_d_up);
  cudaFree(_mms_d_dn);
  cudaFree(_mms_x_up);
  cudaFree(_mms_x_dn);  
  
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

  cudaFree(dev_spin_eo1_up);
  cudaFree(dev_spin_eo1_dn);
  cudaFree(dev_spin_eo3_up);
  cudaFree(dev_spin_eo3_dn);  
  
#ifdef _USE_MPI
  cudaFreeHost(RAND1_UP_D);
  cudaFreeHost(RAND2_UP_D); 
  cudaFreeHost(RAND3_UP_D);
  cudaFreeHost(RAND4_UP_D);  
  cudaFree(RAND_BW_UP_D);
  cudaFree(RAND_FW_UP_D);
  free(R1_UP_D);
  free(R3_UP_D);
  cudaFreeHost(RAND1_DN_D);
  cudaFreeHost(RAND2_DN_D); 
  cudaFreeHost(RAND3_DN_D);
  cudaFreeHost(RAND4_DN_D);  
  cudaFree(RAND_BW_DN_D);
  cudaFree(RAND_FW_DN_D);
  free(R1_DN_D);
  free(R3_DN_D); 
#endif  
  
  for (int i = 0; i < 2; i++) {
    cudaStreamDestroy(stream_nd[i]);
  }
  
}



void checkspin(dev_spinor* s, int N, const char* name){
  printf("spin %s has norm %e\n", name, cublasSdot_wrapper(N, (float *) s, (float *) s));
  return;
}

void checkspin_d(dev_spinor_d* s, int N, const char* name){
  printf("spin %s has norm %e\n", name, cublasDdot_wrapper(N, (double *) s, (double *) s));
  return;
}

// a mixed reliable update mms solver
// arithmetics equal the cpu counterpart in solver/cg_mms_tm_nd.c
// calculation of alphas, betas, zitas ... taken over from there
extern "C" int mixed_cg_mms_eo_nd (spinor ** P_up, spinor ** P_dn,
                                 spinor * Q_up, spinor * Q_dn, double * shifts, int Nshift,
                                 int max_iter, double eps_sq, int rel_prec, int min_solver_it) {
/*
  #ifdef ALGORITHM_BENCHMARK
    double effectiveflops;
    //double hoppingflops = 1608.0;
    double matrixflops  = 14448;
  #endif
  
  // timing
  double startinner, stopinner;  
  double innerclocks;
*/

  cudaError_t cudaerr;

  int check_abs, check_rel;
  if(rel_prec){
    check_rel = 1;
    check_abs = 0;
   }
   else{
    check_rel = 0;
    check_abs = 1;     
  }
  
  /////////////////////
  // INITIALIZATIONS //
  /////////////////////
  
  // FIXME this should be made dependent on if EO is used or not in future
  int N = VOLUME/2;
 
  //blas - only internal volume
  start_blas(N);
  
  int Vol = VOLUMEPLUSRAND/2;
 
  //allocate an auxiliary solver fields 
  spinor ** solver_field = NULL;
  const int nr_sf = 4;
  init_solver_field(&solver_field, VOLUMEPLUSRAND/2, nr_sf);  
 
  
  //set double & single constants
  set_global_sizes();  
  update_constants(dev_grid);  
  update_constants_d(dev_grid);  
  update_bare_constants_nd();
  set_gpu_work_layout(1); //set block and grid sizes, eo!
  
  //spinor fields  
  init_gpu_single_nd_mms_fields(Nshift, Vol);
  #ifdef _USE_MPI
    he_cg_init_nd_additional_mpi<<<1,1>>>(VOLUMEPLUSRAND/2, RAND, g_cart_id, g_nproc);
  #endif
  //->check
  check_mixedsolve_params();  
  
  
  //update double gauge field
  update_gpu_gf_d(g_gauge_field);
  //update single gauge field
  update_gpu_gf(g_gauge_field);  
  //allocate additianal mms fields in single

  #ifdef USETEXTURE	
    bind_texture_gf(dev_gf);
  #endif
  

  
   #ifdef MATRIX_DEBUG
    test_double_nd_operator(Q_up, Q_dn, N_sites_int);
    test_single_nd_operator(Q_up, Q_dn, N_sites_int);
   #endif
     
   size_t dev_spinsize_d; 
   #ifdef _USE_MPI
     dev_spinsize_d = 12*(VOLUME+RAND)/2 * sizeof(dev_spinor_d); // double2 even-odd ! 
   #else
     dev_spinsize_d = 12*VOLUME/2 * sizeof(dev_spinor_d); // double2 even-odd !  
   #endif
   

  
  
  // P_up/dn  can be used as auxiliary field to work on, as it is not later used (could be used as initial guess at the very start)
  // Q_up/dn  can be used as feedback, or if not, also as auxiliary field
  
  int noshifts = Nshift;
  
  //allocate cg constants
  double * sigma;
  double * zitam1, * zita;
  double * alphas, * betas;
  double gamma;
  double alpham1;
    sigma = (double*)calloc((noshifts), sizeof(double));
    zitam1 = (double*)calloc((noshifts), sizeof(double));
    zita = (double*)calloc((noshifts), sizeof(double));
    alphas = (double*)calloc((noshifts), sizeof(double));
    betas = (double*)calloc((noshifts), sizeof(double));


  
  // algorithm
  double rr_up, rr_dn, rr, rr_old, r0r0, dAd_up, dAd_dn, dAd;

  dev_spinor *  r_up, *  r_dn, * Ad_up, * Ad_dn, *  x_up, *  x_dn, *  d_up, *  d_dn;		
  dev_spinor_d * r_up_d, * r_dn_d, * x_up_d, * x_dn_d, * Ax_up_d, * Ax_dn_d, * Q_up_d, * Q_dn_d;
  
 // iteration counter
 int j; 
 
 //reliable update flag
 int rel_update = 0;
 //no of reliable updates done
 int no_rel_update = 0;
 //use reliable update flag
 int use_reliable = 1;
 
 double rel_delta = 1.0e-10;
 int trigger_shift = -1;
 double * res;
 double * res0;
 double * maxres;
 res = (double*)calloc((noshifts), sizeof(double));
 res0 = (double*)calloc((noshifts), sizeof(double));
 maxres = (double*)calloc((noshifts), sizeof(double)); 
    
  /////////////////
  // ASSIGNMENTS //
  /////////////////
  
  x_up  = dev_spinout_up;	
  x_dn  = dev_spinout_dn;	
  r_up  = dev_spin1_up;	
  r_dn  = dev_spin1_dn;
  d_up  = dev_spin2_up;
  d_dn  = dev_spin2_dn;
  Ad_up = dev_spin3_up;
  Ad_dn = dev_spin3_dn;


  Q_up_d = dev_spinin_up_d;
  Q_dn_d = dev_spinin_dn_d;
  x_up_d = dev_spin1_up_d;
  x_dn_d = dev_spin1_dn_d;
  r_up_d = dev_spin2_up_d;
  r_dn_d = dev_spin2_dn_d;
  Ax_up_d = dev_spin3_up_d;
  Ax_dn_d = dev_spin3_dn_d;  
  
  ///////////////
  // ALGORITHM //
  ///////////////
  
  // initialize x(0) = 0	// will be added up
  dev_zero_spinor_field<<<gpu_gd_blas, gpu_bd_blas>>>(x_up);
  dev_zero_spinor_field<<<gpu_gd_blas, gpu_bd_blas>>>(x_dn);
  
  //initialize device double fields
  dev_zero_spinor_field_d<<<gpu_gd_blas_d, gpu_bd_blas_d>>>(x_up_d);
  dev_zero_spinor_field_d<<<gpu_gd_blas_d, gpu_bd_blas_d>>>(x_dn_d); 
  dev_zero_spinor_field_d<<<gpu_gd_blas_d, gpu_bd_blas_d>>>(r_up_d);
  dev_zero_spinor_field_d<<<gpu_gd_blas_d, gpu_bd_blas_d>>>(r_dn_d); 
  dev_zero_spinor_field_d<<<gpu_gd_blas_d, gpu_bd_blas_d>>>(Ax_up_d);
  dev_zero_spinor_field_d<<<gpu_gd_blas_d, gpu_bd_blas_d>>>(Ax_dn_d); 
  order_spin_gpu(Q_up, h2d_spin_d);
  cudaMemcpy(Q_up_d, h2d_spin_d, dev_spinsize_d, cudaMemcpyHostToDevice);
  order_spin_gpu(Q_dn, h2d_spin_d);
  cudaMemcpy(Q_dn_d, h2d_spin_d, dev_spinsize_d, cudaMemcpyHostToDevice);   
  
  #ifdef RELATIVISTIC_BASIS
    to_relativistic_basis_d<<<gpu_gd_linalg_d, gpu_bd_linalg_d>>> (Q_up_d);
    to_relativistic_basis_d<<<gpu_gd_linalg_d, gpu_bd_linalg_d>>> (Q_dn_d);     
  #endif
  
  // r(0) = b
  dev_d2f<<<gpu_gd_blas_d, gpu_bd_blas_d>>>(r_up, Q_up_d);
  dev_d2f<<<gpu_gd_blas_d, gpu_bd_blas_d>>>(r_dn, Q_dn_d); 
  
  // d(0) = b
  dev_d2f<<<gpu_gd_blas, gpu_bd_blas>>>(d_up, Q_up_d);
  dev_d2f<<<gpu_gd_blas, gpu_bd_blas>>>(d_dn, Q_dn_d); 
  
  // norm of source
  rr_up = double_dotprod(Q_up_d, Q_up_d, N_sites_int);
  rr_dn = double_dotprod(Q_dn_d, Q_dn_d, N_sites_int);
  rr    = rr_up + rr_dn;  
  
  r0r0   = rr;	// for relative precision 
  rr_old = rr;	// for the first iteration
  
 if( (g_cart_id == 0 && g_debug_level > 1)) printf("Initial mms residue: %.6e\n", r0r0);
  maxres[0] = rr;
  res[0] = rr;
  res0[0] = rr;
  alphas[0] = 1.0;
  betas[0] = 0.0;
  sigma[0] = shifts[0]*shifts[0];
  if(g_cart_id == 0 && g_debug_level > 2) printf("# dev_CGMMSND_mixed: shift %d is %e\n", 0, sigma[0]);

  // currently only implemented for P=0 
  for(int im = 1; im < noshifts; im++) {
    maxres[im] = rr;
    res[im] = rr;
    res0[im] = rr;    
    sigma[im] = shifts[im]*shifts[im] - sigma[0];
    if(g_cart_id == 0 && g_debug_level > 2) printf("# dev_CGMMSND_mixed: shift %d is %e\n", im, sigma[im]);
    // these will be the result spinor fields
    dev_zero_spinor_field<<<gpu_gd_blas, gpu_bd_blas>>>(mms_x_up[im-1]);
    dev_zero_spinor_field<<<gpu_gd_blas, gpu_bd_blas>>>(mms_x_dn[im-1]);    

    dev_d2f<<<gpu_gd_blas_d, gpu_bd_blas_d>>>(mms_d_up[im-1], Q_up_d);
    dev_d2f<<<gpu_gd_blas_d, gpu_bd_blas_d>>>(mms_d_dn[im-1], Q_dn_d);
    zitam1[im] = 1.0;
    zita[im] = 1.0;
    alphas[im] = 1.0;
    betas[im] = 0.0;
  }

  //zero host fields for solution P_up, P_dn
  for(int im = 0; im < Nshift; im++){
    zero_spinor_field(P_up[im], N);
    zero_spinor_field(P_dn[im], N);    
  }
  
  
  //////////
  // LOOP //
  //////////
    
  for (j = 0; j < max_iter; j++) {   
      // A*d(k)
    dev_Qtm_pm_ndpsi(Ad_up, Ad_dn, d_up,  d_dn,	
		      gpu_gd_M, gpu_bd_M, gpu_gd_linalg, gpu_bd_linalg,
		      gpu_gd_linalg, gpu_bd_linalg, gpu_gd_linalg, gpu_bd_linalg);     
    //add zero'th shift
    cublasSaxpy_wrapper (N_floats_int, sigma[0], (float *) d_up, (float *) Ad_up);
    cublasSaxpy_wrapper (N_floats_int, sigma[0], (float *) d_dn, (float *) Ad_dn);
	
    #ifdef CUDA_DEBUG
      cublasStatus cublasstatus;
      CUDA_CHECK("CUDA error in matrix_muliplication(). Applying the matrix on GPU failed.", "The matrix was applied on GPU.");
    #endif        
    
    // alpha = r(k)*r(k) / d(k)*A*d(k)
    dAd_up = cublasSdot_wrapper(N_floats_int, (float *) d_up, (float *) Ad_up);
    dAd_dn = cublasSdot_wrapper(N_floats_int, (float *) d_dn, (float *) Ad_dn);

    dAd    = dAd_up + dAd_dn; 
    alpham1 = alphas[0];
    alphas[0]  = rr_old / dAd;	// rr_old is taken from the last iteration respectively
    
   
    // r(k+1)
    cublasSaxpy_wrapper(N_floats_int, -1.0*alphas[0], (float *) Ad_up, (float *) r_up);
    cublasSaxpy_wrapper(N_floats_int, -1.0*alphas[0], (float *) Ad_dn, (float *) r_dn);

    // r(k+1)*r(k+1)
    rr_up  = cublasSdot_wrapper(N_floats_int, (float *) r_up, (float *) r_up);
    rr_dn  = cublasSdot_wrapper(N_floats_int, (float *) r_dn, (float *) r_dn);
    rr     = rr_up + rr_dn;
    
    #ifdef CUDA_DEBUG
      CUBLAS_CORE_CHECK_NO_SUCCESS_MSG("CUBLAS error in dev_cg_eo_nd(). CUBLAS function failed.");
    #endif
    // debug	// is NaN ?
    if isnan(rr) {
      if(g_cart_id == 0) printf("Error in mixed_cg_mms_eo_nd(). Inner residue is NaN.\n");
      exit(-1);
    }
      
    #ifndef LOWOUTPUT
      if((g_cart_id == 0) && (g_debug_level > 2)) printf("mms iteration j = %i: rr = %.6e\n", j, rr);
    #endif
		 

    // aborting ?? // check wether precision is reached ...
    if ( ((check_abs)&&(rr <= eps_sq)) || ((check_rel)&&(rr <= eps_sq*r0r0)) ) 
    {
      #ifndef LOWOUTPUT
	if ((check_rel)&&(rr <= eps_sq*r0r0)) {
	  if((g_cart_id == 0) && (g_debug_level > 1)) printf("Reached relative solver precision of eps_rel = %.2e\n", eps_sq);
	}
      #endif 
      break;
   }
    
    // update alphas and zitas  
    // used later
    for(int im = 1; im < noshifts; im++) {
      gamma = zita[im]*alpham1/(alphas[0]*betas[0]*(1.-zita[im]/zitam1[im]) 
				+ alpham1*(1.+sigma[im]*alphas[0]));
      zitam1[im] = zita[im];
      zita[im] = gamma;
      alphas[im] = alphas[0]*zita[im]/zitam1[im];
    }  
    
    //check for reliable update
    res[0] = rr;
    for(int im=1; im<noshifts; im++) res[im] = rr * zita[im]; 
      
    rel_update = 0;
    for(int im = (noshifts-1); im >= 0; im--) {
      if( res[im] > maxres[im] ) maxres[im] = res[im];
      if( (res[im] < rel_delta*res0[im]) && (res0[im]<=maxres[im]) && (use_reliable) ) rel_update=1; 
      if( rel_update && ( trigger_shift == -1) ) trigger_shift = im;
    }     
    
    if(!rel_update)
    {
      // x_j(k+1) = x_j(k) + alpha_j*d_j(k) 
      // alphas are set above
      cublasSaxpy_wrapper(N_floats_int, alphas[0], (float *) d_up, (float *) x_up);   
      cublasSaxpy_wrapper(N_floats_int, alphas[0], (float *) d_dn, (float *) x_dn);
      for(int im = 1; im < noshifts; im++) {
	cublasSaxpy_wrapper(N_floats_int, alphas[im], (float *) mms_d_up[im-1], (float *) mms_x_up[im-1]);   
	cublasSaxpy_wrapper(N_floats_int, alphas[im], (float *) mms_d_dn[im-1], (float *) mms_x_dn[im-1]);  
      }  
   
      // beta = r(k+1)*r(k+1) / r(k)*r(k)
      betas[0] = rr / rr_old;
      rr_old = rr;  // for next iteration
      
      // d_0(k+1) = r(k+1) + beta*d_0(k)
      cublasSscal_wrapper (N_floats_int, betas[0], (float *) d_up);
      cublasSaxpy_wrapper (N_floats_int, 1.0 , (float *) r_up, (float *) d_up);     
      cublasSscal_wrapper (N_floats_int, betas[0], (float *) d_dn);
      cublasSaxpy_wrapper (N_floats_int, 1.0 , (float *) r_dn, (float *) d_dn);     
      // d_j(k+1) = r(k+1) + beta*d_j(k)
      for(int im = 1; im < noshifts; im++) {
	betas[im] = betas[0]*zita[im]*alphas[im]/(zitam1[im]*alphas[0]);
	cublasSscal_wrapper (N_floats_int, betas[im], (float *) mms_d_up[im-1]);
	cublasSaxpy_wrapper (N_floats_int, zita[im] , (float *) r_up, (float *) mms_d_up[im-1]);
      
	cublasSscal_wrapper (N_floats_int, betas[im], (float *) mms_d_dn[im-1]);
	cublasSaxpy_wrapper (N_floats_int, zita[im] , (float *) r_dn, (float *) mms_d_dn[im-1]);
      }   
    }
    else{
      //reliable update
      if( (g_cart_id == 0) && (g_debug_level > 1) ){
	printf("Shift %d with offset squared %e triggered a reliable update\n", trigger_shift, sigma[trigger_shift]);
      }
      //add low prec solutions on host 
      cublasSaxpy_wrapper(N_floats_int, alphas[0], (float *) d_up, (float *) x_up);
      cublasSaxpy_wrapper(N_floats_int, alphas[0], (float *) d_dn, (float *) x_dn);
      add_f2d_host(P_up[0], solver_field[0], x_up, Vol);
      add_f2d_host(P_dn[0], solver_field[0], x_dn, Vol);	    
      for(int im = 1; im < noshifts; im++) {  
	cublasSaxpy_wrapper(N_floats_int, alphas[im], (float *) mms_d_up[im-1], (float *) mms_x_up[im-1]);
	cublasSaxpy_wrapper(N_floats_int, alphas[im], (float *) mms_d_dn[im-1], (float *) mms_x_dn[im-1]);	
	add_f2d_host(P_up[im], solver_field[0], mms_x_up[im-1], Vol);
        add_f2d_host(P_dn[im], solver_field[0], mms_x_dn[im-1], Vol);	
      }
      
      //add low precision on device for shift 0 only
      dev_add_f2d<<<gpu_gd_blas_d, gpu_bd_blas_d >>>(x_up_d, x_up_d, x_up); 
      dev_add_f2d<<<gpu_gd_blas_d, gpu_bd_blas_d >>>(x_dn_d, x_dn_d, x_dn);      
//      checkspin(x_dn, N_floats_int, "x_dn");
//      checkspin_d(x_dn_d, N_floats_int, "x_dn_d");   
      
      dev_Qtm_pm_ndpsi_d(Ax_up_d, Ax_dn_d, x_up_d,  x_dn_d,	
		      gpu_gd_M_d, gpu_bd_M_d, gpu_gd_linalg_d, gpu_bd_linalg_d,
		      gpu_gd_linalg_d, gpu_bd_linalg_d, gpu_gd_linalg_d, gpu_bd_linalg_d);
      //add zero'th shift
      cublasDaxpy_wrapper (N_floats_int, sigma[0], (double *) x_up_d, (double *) Ax_up_d);
      cublasDaxpy_wrapper (N_floats_int, sigma[0], (double *) x_dn_d, (double *) Ax_dn_d);
      #ifdef CUDA_DEBUG
	CUDA_CHECK("CUDA error in matrix_muliplication(). Applying the matrix on GPU failed.", "The matrix was applied on GPU.");
      #endif
      dev_diff_d<<<gpu_gd_linalg_d,gpu_bd_linalg_d>>>(r_up_d,Q_up_d,Ax_up_d);         
      dev_diff_d<<<gpu_gd_linalg_d,gpu_bd_linalg_d>>>(r_dn_d,Q_dn_d,Ax_dn_d); 
 
      rr_up = double_dotprod(r_up_d, r_up_d, N_sites_int);
      rr_dn = double_dotprod(r_dn_d, r_dn_d, N_sites_int);
      rr    = rr_up + rr_dn;
      if ((g_cart_id == 0) && (g_debug_level > 1) ) printf("New residue after reliable update: %.6e\n", rr);
       
      //update res[im]
      res[0] = rr;
      //for(int im=1; im<noshifts; im++) res[im] = rr * zita[im]* zita[im];

       
      if(res[trigger_shift] > res0[trigger_shift]){
	if(g_cart_id == 0) printf("Warning: residue of shift no %d got larger after rel. update\n", trigger_shift);
	//if this is the zero'th shift not getting better -> no further convergence, break
	if(trigger_shift == 0) break;
      }    
      
      //zero float fields
      dev_zero_spinor_field<<<gpu_gd_blas, gpu_bd_blas>>>(x_up);
      dev_zero_spinor_field<<<gpu_gd_blas, gpu_bd_blas>>>(x_dn);        
      for(int im = 1; im < noshifts; im++) {
	dev_zero_spinor_field<<<gpu_gd_blas, gpu_bd_blas>>>(mms_x_up[im-1]);
	dev_zero_spinor_field<<<gpu_gd_blas, gpu_bd_blas>>>(mms_x_dn[im-1]);  
      }
      
      //update the source
      dev_d2f<<<gpu_gd_blas_d, gpu_bd_blas_d>>>(r_up, r_up_d);
      dev_d2f<<<gpu_gd_blas_d, gpu_bd_blas_d>>>(r_dn, r_dn_d); 
      

      
      betas[0] = res[0]/rr_old;
      rr_old = rr;
      // d_0(k+1) = r(k+1) + beta*d_0(k)
      cublasSscal_wrapper (N_floats_int, betas[0], (float *) d_up);
      cublasSaxpy_wrapper (N_floats_int, 1.0 , (float *) r_up, (float *) d_up);     
      cublasSscal_wrapper (N_floats_int, betas[0], (float *) d_dn);
      cublasSaxpy_wrapper (N_floats_int, 1.0 , (float *) r_dn, (float *) d_dn);     
      // d_j(k+1) = r(k+1) + beta*d_j(k)
      for(int im = 1; im < noshifts; im++) {
	betas[im] = betas[0]*zita[im]*alphas[im]/(zitam1[im]*alphas[0]);
	cublasSscal_wrapper (N_floats_int, betas[im], (float *) mms_d_up[im-1]);
	cublasSaxpy_wrapper (N_floats_int, zita[im] , (float *) r_up, (float *) mms_d_up[im-1]);
      
	cublasSscal_wrapper (N_floats_int, betas[im], (float *) mms_d_dn[im-1]);
	cublasSaxpy_wrapper (N_floats_int, zita[im] , (float *) r_dn, (float *) mms_d_dn[im-1]);
      } 
      
      //new maxres for the shift that initiated the reliable update
      res[trigger_shift] = res[0]*zita[trigger_shift]*zita[trigger_shift];
      res0[trigger_shift] = res[trigger_shift];  
      maxres[trigger_shift] = res[trigger_shift];
      trigger_shift = -1;
      no_rel_update ++;
    }	//reliable update	
    
    //check if some shift is converged
    for(int im = 1; im < noshifts; im++) {    
      if(j > 0 && (j % 10 == 0) && (im == noshifts-1)) {
	double sn = cublasSdot_wrapper(N_floats_int, (float *) mms_d_up[im-1], (float *) mms_d_up[im-1]);
	sn += cublasSdot_wrapper(N_floats_int, (float *) mms_d_dn[im-1], (float *) mms_d_dn[im-1]);
	if(alphas[noshifts-1]*alphas[noshifts-1]*sn <= eps_sq) {
	  noshifts--;
	  if( (g_debug_level > 1) && (g_cart_id == 0) ) {
	    printf("# dev_CGMMSND_mixed: at iteration %d removed one shift, %d remaining\n", j, noshifts);
	  }
	  //if removed we add the latest solution vector for this shift on host	  
	  add_f2d_host(P_up[im], solver_field[0], mms_x_up[im-1], Vol);
          add_f2d_host(P_dn[im], solver_field[0], mms_x_dn[im-1], Vol);
	}
      }
    }
    
    // debug	// CUBLAS core function
    #ifdef CUDA_DEBUG
      CUBLAS_CORE_CHECK_NO_SUCCESS_MSG("CUBLAS error in dev_cg_eo_nd(). Error in CUBLAS function.");
    #endif
      
  }//LOOP
  
  if( (g_cart_id == 0) && (g_debug_level > 1) ) printf("Final mms residue: %.6e\n", rr);

  //add the latest solutions on host
  for(int im = 0; im < noshifts; im++) {  
    if(im == 0){   
      add_f2d_host(P_up[0], solver_field[0], x_up, Vol);
      add_f2d_host(P_dn[0], solver_field[0], x_dn, Vol);        
    }
    else{     
      add_f2d_host(P_up[im], solver_field[0], mms_x_up[im-1], Vol);
      add_f2d_host(P_dn[im], solver_field[0], mms_x_dn[im-1], Vol);      
    }
  }  
#ifdef MATRIX_DEBUG
    if(g_cart_id == 0) printf("# dev_CGMMSND_mixed: Checking mms result:\n");
    //loop over all shifts (-> Nshift) 
    for(int im = 0; im < Nshift; im++){
      Qtm_pm_ndpsi(solver_field[0], solver_field[1], P_up[im], P_dn[im]);
      assign_add_mul_r(solver_field[0], P_up[im] , shifts[im]*shifts[im], N);
      assign_add_mul_r(solver_field[1], P_dn[im] , shifts[im]*shifts[im], N);
      diff(solver_field[2], solver_field[0], Q_up, N);
      diff(solver_field[3], solver_field[1], Q_dn, N);
      rr_up = square_norm(solver_field[2], N, 1);
      rr_dn = square_norm(solver_field[3], N, 1);      
      rr = rr_up + rr_dn;
      if(g_cart_id == 0) printf("# dev_CGMMSND_mixed: Shift[%d] squared residue: %e\n", im, rr);
    }
#endif
  
  if((cudaerr=cudaPeekAtLastError())!=cudaSuccess){
    if(g_cart_id==0){
      printf("Error in mixed_cg_mms_eo_nd().\n");
      printf("Error was %d. Aborting...\n", cudaerr);
    }
    exit(200);
  } 
  
  finalize_gpu_single_nd_mms_fields();  
  finalize_solver(solver_field, nr_sf);    

  #ifdef USETEXTURE
    unbind_texture_gf();
  #endif 
  
  stop_blas();
  //free cg constants
  free(sigma); free(zitam1); free(zita); free(alphas); free(betas);    
  
  //free reliable update stuff
  free(res); free(res0); free(maxres);


    
  
  //if not converged -> return(-1)
  if(j<max_iter){
    return(j);
  }
  else{
    return(-1);
  }
}//








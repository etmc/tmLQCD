/***********************************************************************
 *
 * Copyright (C) 2010 Florian Burger
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
 *  
 * File: mixed_solve.cu
 *
 * CUDA GPU mixed_solver for EO and non-EO
 * CUDA kernels for Hopping-Matrix and D_tm
 *
 * The externally accessible functions are
 *
 *
 *   extern "C" int mixed_solve_eo (spinor * const P, spinor * const Q, const int max_iter, 
           double eps, const int rel_prec, const int N)
 *
 *  extern "C" int mixed_solve (spinor * const P, spinor * const Q, const int max_iter, 
           double eps, const int rel_prec,const int N)
 *
 * input:
 *   Q: source
 * inout:
 *   P: initial guess and result
 * 
 *
 **************************************************************************/


#include <cuda.h>
#include <cuda_runtime.h>

#ifdef CUDA_45  
  #include "cublas_v2.h"
#else
  #include "cublas.h"
#endif 



#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

extern "C" {
#include "../global.h"
}
#include "../hamiltonian_field.h"
#include "cudaglobal.h"
#include "../solver/solver.h"
#include "HEADER.h"
#include "cudadefs.h"
#include <math.h>


extern "C" {
#include "../operator/tm_operators.h"
#include "../linalg_eo.h"
#include "../start.h"
#include "../read_input.h"
#include "../geometry_eo.h"
#include "../boundary.h"
#include "../su3.h"
#include "../temporalgauge.h"
#include "../measure_gauge_action.h"
#include "../measure_rectangles.h"
#include "../polyakov_loop.h"
#include "../su3spinor.h"
#include "../solver/solver_field.h"

#include "../gettime.h"
#ifdef _USE_MPI
  #include "../xchange/xchange.h"
#endif 

#include "../phmc.h"
}


#ifdef HAVE_CONFIG_H
  #include<config.h>
#endif


#include "MACROS.cuh"




int g_numofgpu;

#ifdef GF_8
dev_su3_8 * dev_gf;
dev_su3_8 * h2d_gf;
#else
dev_su3_2v * dev_gf;
dev_su3_2v * h2d_gf;
#endif

//double gauge field and deri field
double2 * dev_gf_d;
dev_su3_2v_d * h2d_gf_d;
dev_su3adj* dev_df0_d;
dev_su3adj* h2d_df0_d;

#ifndef HALF
dev_spinor* dev_spin1;
dev_spinor* dev_spin2;
dev_spinor* dev_spin3;
dev_spinor* dev_spin4;
dev_spinor* dev_spin5;
dev_spinor* dev_spinin;
dev_spinor* dev_spinout;
dev_spinor * h2d_spin;

//additional spinors for even-odd
dev_spinor* dev_spin_eo1;
dev_spinor* dev_spin_eo2;


//double spinor fields for outer loop in double on GPU

dev_spinor_d* dev_spin0_d; 
dev_spinor_d* dev_spin1_d;  
dev_spinor_d* dev_spin2_d;  
dev_spinor_d* dev_spin3_d;
dev_spinor_d* dev_spin_eo1_d;	
dev_spinor_d* dev_spin_eo2_d;
dev_spinor_d * h2d_spin_d;

dev_spinor_d * dev_spin_eo1_up_d;
dev_spinor_d * dev_spin_eo1_dn_d;
dev_spinor_d * dev_spin_eo2_up_d;
dev_spinor_d * dev_spin_eo2_dn_d;
dev_spinor_d * dev_spin_eo3_up_d;
dev_spinor_d * dev_spin_eo3_dn_d;


dev_spinor_d * dev_spinin_up_d;
dev_spinor_d * dev_spinin_dn_d;
dev_spinor_d * dev_spinout_up_d;
dev_spinor_d * dev_spinout_dn_d;

dev_spinor_d * dev_spin1_up_d;
dev_spinor_d * dev_spin1_dn_d;
dev_spinor_d * dev_spin2_up_d;
dev_spinor_d * dev_spin2_dn_d;
dev_spinor_d * dev_spin3_up_d;
dev_spinor_d * dev_spin3_dn_d;

//spinor fields for reliable update solver
dev_spinor* dev_spin1_reliable;
dev_spinor* dev_spin2_reliable;
dev_spinor* dev_spin3_reliable;
dev_spinor* h2d_spin_reliable;

#else

dev_spinor_half* dev_spin1;
dev_spinor_half* dev_spin2;
dev_spinor_half* dev_spin3;
dev_spinor_half* dev_spin4;
dev_spinor_half* dev_spin5;
dev_spinor_half* dev_spinin;
dev_spinor_half* dev_spinout;
dev_spinor_half* h2d_spin;
//additional spinors for even-odd



dev_spinor* dev_spin_eo1;
dev_spinor* dev_spin_eo2;

dev_spinor_half* dev_spin_eo1_half;
dev_spinor_half* dev_spin_eo2_half;

//spinor fields for reliable update solver
dev_spinor* dev_spin1_reliable;
dev_spinor* dev_spin2_reliable;
dev_spinor* dev_spin3_reliable;
dev_spinor* h2d_spin_reliable;


float* dev_spin1_norm;
float* dev_spin2_norm;
float* dev_spin3_norm;
float* dev_spin4_norm;
float* dev_spin5_norm;
float* dev_spinin_norm;
float* dev_spinout_norm;
float* h2d_spin_norm;

float* dev_spin_eo1_half_norm;
float* dev_spin_eo2_half_norm;


  // a half precsion gauge field
  #ifdef GF_8
   dev_su3_8_half * dev_gf_half;
  #else
   dev_su3_2v_half * dev_gf_half;
  #endif
#endif 


//nearest neighbour fields for Dirac
int * nn;
int * nn_eo;
int * nn_oe;
int * eoidx_even;
int * eoidx_odd;

int * dev_nn;
int * dev_nn_eo;
int * dev_nn_oe;

int * dev_eoidx_even;
int * dev_eoidx_odd;

//nearest neighbour fields for gauge deri
int * nn2;
int * dev_nn2;


int* dev_grid;
int havedevice = 0;


float hostr;

float hostm;
float hostmu;

/*********** float constants on GPU *********************/
__device__  float m;
__device__  float mu;
__device__  float r=1.0; // this is implicitly assumed to be 1.0 in the host code!!!
__device__  float kappa;
__device__ float twokappamu;

__device__ dev_complex dev_k0;
__device__ dev_complex dev_k1;
__device__ dev_complex dev_k2;
__device__ dev_complex dev_k3;

__device__ dev_complex dev_mk0;
__device__ dev_complex dev_mk1;
__device__ dev_complex dev_mk2;
__device__ dev_complex dev_mk3;


// __constant__ __device__ dev_complex dev_k0c;
// __constant__ __device__ dev_complex dev_k1c;
// __constant__ __device__ dev_complex dev_k2c;
// __constant__ __device__ dev_complex dev_k3c;
// 
// __constant__ __device__ dev_complex dev_mk0c;
// __constant__ __device__ dev_complex dev_mk1c;
// __constant__ __device__ dev_complex dev_mk2c;
// __constant__ __device__ dev_complex dev_mk3c;


// physical parameters (on device)
__device__ float mubar, epsbar;


/*********** double constants on GPU *********************/
__device__ double mu_d;
__device__ double kappa_d;
__device__ double twokappamu_d;
__device__ double mubar_d, epsbar_d;

__device__ dev_complex_d dev_k0_d;
__device__ dev_complex_d dev_k1_d;
__device__ dev_complex_d dev_k2_d;
__device__ dev_complex_d dev_k3_d;

__device__ dev_complex_d dev_mk0_d;
__device__ dev_complex_d dev_mk1_d;
__device__ dev_complex_d dev_mk2_d;
__device__ dev_complex_d dev_mk3_d;




//dev_Offset is the jump in the spinor fields we have to do because of space-time first ordering
__device__  int  dev_LX,dev_LY,dev_LZ,dev_T,dev_VOLUME,dev_Offset,dev_VOLUMEPLUSRAND;

__device__  int  dev_VOL2;


//host_Offset is the jump in the spinor fields we have to do because of space-time first ordering on host
int host_Offset;


int nStreams_nd = 2;
cudaStream_t stream_nd[2];


// here we set the block and grid sizes for kernels
int gpu_gd_M, gpu_bd_M, gpu_gd_linalg, gpu_bd_linalg, gpu_gd_blas, gpu_bd_blas;
int gpu_gd_M_d, gpu_bd_M_d, gpu_gd_linalg_d, gpu_bd_linalg_d, gpu_gd_blas_d, gpu_bd_blas_d;
// for gpu_gauge_derivative
int blockdimgauge, griddimgauge;
// for gpu_deriv_SB
int blockdimdsb,  griddimdsb;

#ifdef _USE_MPI


// from mixed_solve_eo_nd.cuh
__device__ int dev_RAND;                        // not used, maybe later ...
__device__ int dev_rank;
__device__ int dev_nproc;



    dev_spinor * R1;
    dev_spinor * R2;
    dev_spinor * R3;
    dev_spinor * R4;
    dev_spinor_d * R1_D;
    dev_spinor_d * R2_D;
    dev_spinor_d * R3_D;
    dev_spinor_d * R4_D;
    dev_spinor_d * R1_UP_D;
    dev_spinor_d * R2_UP_D;
    dev_spinor_d * R3_UP_D;
    dev_spinor_d * R4_UP_D; 
    dev_spinor_d * R1_DN_D;
    dev_spinor_d * R2_DN_D;
    dev_spinor_d * R3_DN_D;
    dev_spinor_d * R4_DN_D; 


  dev_spinor* RAND_FW;
  dev_spinor* RAND_BW;
  dev_spinor_d* RAND_FW_D;
  dev_spinor_d* RAND_BW_D;
  dev_spinor_d* RAND_FW_UP_D;
  dev_spinor_d* RAND_BW_UP_D;
  dev_spinor_d* RAND_FW_DN_D;
  dev_spinor_d* RAND_BW_DN_D;
  
//also need RAND? for ASYNC==0
//#if ASYNC > 0
    int nStreams = ASYNC_OPTIMIZED;
    cudaStream_t stream[2*ASYNC_OPTIMIZED+1];
    
    dev_spinor_d  * RAND1_D;   // for exchanging the boundaries in ASYNC_D.cuh
    dev_spinor_d  * RAND2_D;
    dev_spinor_d  * RAND3_D; 
    dev_spinor_d  * RAND4_D;
 
    dev_spinor_d  * RAND1_UP_D;   // for exchanging the boundaries in ASYNC_D.cuh
    dev_spinor_d  * RAND2_UP_D;
    dev_spinor_d  * RAND3_UP_D; 
    dev_spinor_d  * RAND4_UP_D;

    dev_spinor_d  * RAND1_DN_D;   // for exchanging the boundaries in ASYNC_D.cuh
    dev_spinor_d  * RAND2_DN_D;
    dev_spinor_d  * RAND3_DN_D; 
    dev_spinor_d  * RAND4_DN_D;
    
   #ifndef HALF
    dev_spinor * RAND1;   // for exchanging the boundaries in ASYNC.cuh
    dev_spinor * RAND2;
    dev_spinor * RAND3; // page-locked memory
    dev_spinor * RAND4;
   #else
     dev_spinor_half * RAND1;   // for exchanging the boundaries in ASYNC.cuh
     dev_spinor_half * RAND2;
     dev_spinor_half * RAND3; // page-locked memory
     dev_spinor_half * RAND4;
     //we also need page-locked norms
      float * RAND1_norm;  
      float * RAND2_norm;
      float * RAND3_norm; 
      float * RAND4_norm;
    #endif
//#endif




  MPI_Status stat[2];
  MPI_Request send_req[2];
  MPI_Request recv_req[2]; 



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

#endif //MPI





// include files with other GPU code as all GPU code has to reside in one file 
// the texture references and functions
#include "textures.cuh"
// if we want to use half precision
#ifdef HALF 
 #include "half.cuh"
#endif
// linear algebra functions and gamma-multiplications
#include "linalg.cuh"
// reconstruction of the gauge field
#include "gauge_reconstruction.cuh"
// the device su3 functions
#include "su3.cuh"
// the plaquette and rectangle routines
#include "observables.cuh"

// linear algebra for double
#include "linalg_d.cuh"

/* include double matrix on device here */
#include "Hopping_Matrix_d.cuh"

//gauge staple calculations
#include "gauge_monomial.cuh"
/* include gauge_update on device here */
#include "hybrid_update.cuh"
/* include deriv_Sb on device here */
#include "deriv_Sb.cuh"
/* include the wilson flow on device*/
//#include "wilson_flow.cuh"



// the device Hopping_Matrix
#include "Hopping_Matrix.cuh"


// eo tm operator in single and half precision
#include "tm_eo.cuh"
// the non-EO twisted mass dirac operator
#include "tm_diracoperator.cuh"
// mixed solver, even/odd, non-degenerate two flavour
#include "mixed_solve_eo_nd.cuh"
// nd eo tm 1+1 operator in single precision
#include "tm_nd_eo.cuh"

#ifdef _USE_MPI
// optimization of the communication
  #include "ASYNC.cuh"
  #include "ASYNC_D.cuh"
#endif 

// nd-mms solver based on single mass solver and polynomial initial guess
#include "nd_mms.cuh"

#include "tmclover.cuh"









// init the gpu inner solver, assigen constants etc.
__global__ void he_cg_init (int* grid, float param_kappa, float param_mu, dev_complex k0, dev_complex k1, dev_complex k2, dev_complex k3){
  dev_LX = grid[0];
  dev_LY = grid[1];
  dev_LZ = grid[2];
  dev_T = grid[3];
  dev_VOLUME = grid[4]; // grid[4] is initialized 1/2 VOLUME for eo
  dev_Offset = grid[5]; //this is the offset for the spinor fields
  dev_VOLUMEPLUSRAND = grid[5]; 
  
  kappa = param_kappa;
  mu = param_mu;
  twokappamu = 2.0*param_kappa*param_mu;
  
  dev_k0.re = k0.re;
  dev_k0.im = k0.im;
  dev_mk0.re = -k0.re;
  dev_mk0.im = -k0.im;
  
  dev_k1.re = k1.re;
  dev_k1.im = k1.im;
  dev_mk1.re = -k1.re;
  dev_mk1.im = -k1.im;
  
  dev_k2.re = k2.re;
  dev_k2.im = k2.im;
  dev_mk2.re = -k2.re;
  dev_mk2.im = -k2.im;
  
  dev_k3.re = k3.re;
  dev_k3.im = k3.im;
  dev_mk3.re = -k3.re;
  dev_mk3.im = -k3.im;
}






// init the gpu, assign dimensions 
__global__ void dev_init_grid (int* grid){
  dev_LX = grid[0];
  dev_LY = grid[1];
  dev_LZ = grid[2];
  dev_T = grid[3];
  dev_VOLUME = grid[4]; // grid[4] is initialized 1/2 VOLUME for eo
  dev_Offset = grid[5];
  dev_VOLUMEPLUSRAND = grid[5];
}




void update_constants(int *grid){
  dev_complex h0,h1,h2,h3;
  
  //Hopping Matrix and tm_Dirac_op are defined with a relative minus in the imaginary parts of kappa
  //both comply with the cpu part
  float sign;
  if(even_odd_flag){
    sign=-1.0;
  }
  else{
    sign=-1.0;
  }
 
  h0.re = (float)creal(ka0);    h0.im = sign*(float)cimag(ka0);
  h1.re = (float)creal(ka1);    h1.im = (float)cimag(ka1);
  h2.re = (float)creal(ka2);    h2.im = (float)cimag(ka2);
  h3.re = (float)creal(ka3);    h3.im = (float)cimag(ka3);
  

  #ifndef LOWOUTPUT
  if(g_cart_id==0){
    printf("ka0.re = %f\n",  h0.re);
    printf("ka0.im = %f\n",  h0.im); 
    printf("ka1.re = %f\n",  h1.re);
    printf("ka1.im = %f\n",  h1.im); 
    printf("ka2.re = %f\n",  h2.re);
    printf("ka2.im = %f\n",  h2.im);    
    printf("ka3.re = %f\n",  h3.re);
    printf("ka3.im = %f\n",  h3.im);      
    
    printf("mu = %f\n", g_mu/(2.0*g_kappa));
    printf("2kappamu = %f\n", g_mu);
  }
  #endif
   
  
  he_cg_init<<< 1, 1 >>> (grid, (float) g_kappa, (float)(g_mu/(2.0*g_kappa)), h0,h1,h2,h3);
  // BEWARE in dev_tm_dirac_kappa we need the true mu (not 2 kappa mu!)
}





// code to list available devices, not yet included in main code
// this is copied from the CUDA sdk 
extern "C" int find_devices() {

  int deviceCount, dev;

  cudaGetDeviceCount(&deviceCount);
    
  #ifdef _USE_MPI
    if (g_cart_id == 0) {
  #endif
    
    if (deviceCount == 0)
        printf("There is no device supporting CUDA\n");
    for (dev = 0; dev < deviceCount; ++dev) {
        cudaDeviceProp deviceProp;
        cudaGetDeviceProperties(&deviceProp, dev);
        if (dev == 0) {
            if (deviceProp.major == 9999 && deviceProp.minor == 9999)
                printf("There is no device supporting CUDA.\n");
            else if (deviceCount == 1)
                printf("There is 1 device supporting CUDA\n");
            else
                printf("There are %d devices supporting CUDA\n", deviceCount);
        }
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
    
    #ifdef _USE_MPI 
      }
    #endif
    
    return(deviceCount);
}




void set_gpu_work_layout(int eoflag){
  
 int N;
 if (eoflag){
   N = VOLUME/2;  
 }
 else{
   N = VOLUME;
 }
 
 // this is the partitioning for the Hopping Matrix Mult.
 // single
 gpu_bd_M = BLOCK;
 if( N % gpu_bd_M == 0){
   gpu_gd_M = (int) N/gpu_bd_M;
 }
 else{
   gpu_gd_M  = (int) N/gpu_bd_M + 1;
 }
 // double
 gpu_bd_M_d = BLOCKD;
 if( N % gpu_bd_M_d == 0){
   gpu_gd_M_d = (int) N/gpu_bd_M_d;
 }
 else{
   gpu_gd_M_d  = (int) N/gpu_bd_M_d + 1;
 } 
 
 // this is the partitioning for the dev_mul_one_pm... 
 // single
 gpu_bd_linalg = BLOCK2;
 if( N % gpu_bd_linalg == 0){
   gpu_gd_linalg = (int) N/gpu_bd_linalg;
 }
 else{
   gpu_gd_linalg  = (int) N/gpu_bd_linalg + 1;
 }
 // double
 gpu_bd_linalg_d = BLOCK2D;
 if( N % gpu_bd_linalg_d == 0){
   gpu_gd_linalg_d = (int) N/gpu_bd_linalg_d;
 }
 else{
   gpu_gd_linalg_d  = (int) N/gpu_bd_linalg_d + 1;
 } 
  
  
 // this is the partitioning for the copying of fields 
 // single and double
 gpu_bd_blas = BLOCK3;
 gpu_bd_blas_d = gpu_bd_blas;
 if( N % gpu_bd_blas == 0){
   gpu_gd_blas = (int) N/gpu_bd_blas;
 }
 else{
   gpu_gd_blas  = (int) N/gpu_bd_blas + 1;
 }
 gpu_gd_blas_d = gpu_gd_blas;
 
}





extern "C" void test_operator(dev_su3_2v * gf,dev_spinor* spinin, dev_spinor* spinout, 
dev_spinor* spin0, dev_spinor* spin1, dev_spinor* spin2, dev_spinor* spin3, dev_spinor* spin4, int *grid, int * nn_grid, float* output,float* erg, int xsize, int ysize){
 

 update_constants(grid);
 
  float scaleparam = sqrt(1.0/(2.0 * (float) g_kappa));
  dev_skalarmult_spinor_field<<<gpu_gd_blas, gpu_bd_blas >>>(spinin,scaleparam*scaleparam, spin4);
 
 #ifdef USETEXTURE
   bind_texture_gf(gf);
   bind_texture_spin(spin4,1);
 #endif 
  // apply D_tm
  dev_tm_dirac_kappa <<<gpu_gd_M, gpu_bd_M >>>(gf, spin4, spinout, nn_grid);

 #ifdef USETEXTURE
  unbind_texture_gf();
  unbind_texture_spin(1);
 #endif
}





// this is the eo version of the device cg inner solver 
// we invert the hermitean D_tm D_tm^{+}
extern "C" int dev_cg(
       dev_su3_2v * gf,
       dev_spinor* spinin, 
       dev_spinor* spinout, 
       dev_spinor* spin0, 
       dev_spinor* spin1, 
       dev_spinor* spin2, 
       dev_spinor* spin3, 
       dev_spinor* spin4, 
       int *grid, int * nn_grid, int rescalekappa){
 
 
 float alpha, beta, dotprod, rk, sourcesquarenorm;
 
 
 cudaError_t cudaerr;
 int i;
 int maxit = max_innersolver_it;
 float eps = (float) innersolver_precision;
 int N_recalcres = 30; // after N_recalcres iterations calculate r = A x_k - b
 
 
 // initialize grid and block, make sure VOLUME is a multiple of blocksize 
 if(VOLUME%DOTPROD_DIM != 0){
   printf("Error: VOLUME is not a multiple of DOTPROD_DIM. Aborting...\n");
   exit(100); 
 }

 
 size_t size2 = sizeof(float4)*6*VOLUME;
 
 #ifdef USETEXTURE
   //Bind texture gf
   bind_texture_gf(gf);
  //Bind texture spinor to spin4 (D_tm is always applied to spin4)
  bind_texture_spin(spin4,1);
 #endif
 
 //Initialize some stuff
 update_constants(grid);
 
 //init blas
 start_blas(VOLUME);


 
 // go over to kappa (if wanted)
 float scaleparam = sqrt(1.0/(2.0 * g_kappa));
 printf("1/2kappa = %.8f\n",scaleparam);
 //dev_skalarmult_spinor_field<<<griddim2, blockdim2 >>>(spinin,scaleparam, spin1);
 //dev_copy_spinor_field<<<griddim2, blockdim2 >>>(spin1, spinin);
 
  #ifdef RELATIVISTIC_BASIS 
   //transform to relativistic gamma basis
   to_relativistic_basis<<<gpu_gd_linalg, gpu_bd_linalg>>> (spinin);

   if ((cudaerr=cudaGetLastError())!=cudaSuccess) {
     if (g_cart_id == 0) printf("%s\n", cudaGetErrorString(cudaGetLastError()));
   }
   else{
     #ifndef LOWOUTPUT 
     if (g_cart_id == 0) printf("Switched to relativistic basis\n");
     #endif
   }
 #endif
 
 dev_copy_spinor_field<<<gpu_gd_blas, gpu_bd_blas >>>(spinin, spin0);
 dev_zero_spinor_field<<<gpu_gd_blas, gpu_bd_blas >>>(spin1); // x_0 = 0
 dev_copy_spinor_field<<<gpu_gd_blas, gpu_bd_blas >>>(spinin, spin2);
 dev_zero_spinor_field<<<gpu_gd_blas, gpu_bd_blas >>>(spin3);
 printf("%s\n", cudaGetErrorString(cudaGetLastError()));
 
 
 
 
 //relative precision -> get initial residue
 sourcesquarenorm = cublasSdot (24*VOLUME, (const float *)spinin, 1, (const float *)spinin, 1);
 rk = sourcesquarenorm; //for use in main loop
 printf("Squarenorm Source:\t%.8e\n", sourcesquarenorm);
 printf("%s\n", cudaGetErrorString(cudaGetLastError()));
 
  printf("Entering cg-loop\n");
 for(i=0;i<maxit;i++){ //MAIN LOOP
  
  // D Ddagger    --   Ddagger = gamma5 D gamma5  for Wilson Dirac Operator
  // mu -> -mu for twisted term
  // DO NOT USE tm_dirac_dagger_kappa here, otherwise spin2 will be overwritten!!!
  #ifdef USETEXTURE
    unbind_texture_spin(1);
  #endif
     // GAMMA5, mu -> -mu
  #ifdef RELATIVISTIC_BASIS 
     dev_gamma5_rel <<<gpu_gd_blas, gpu_bd_blas >>> (spin2,spin4);
  #else
     dev_gamma5 <<<gpu_gd_blas, gpu_bd_blas >>> (spin2,spin4);     
  #endif
     dev_swapmu <<<1,1>>> ();
  #ifdef USETEXTURE
   bind_texture_spin(spin4,1);
  #endif
     //D_tm 
     dev_tm_dirac_kappa <<<gpu_gd_M, gpu_bd_M >>> (gf, spin4, spin3, dev_nn);
  #ifdef USETEXTURE
   unbind_texture_spin(1);
  #endif
     //GAMMA5 mu -> -mu
  #ifdef RELATIVISTIC_BASIS    
     dev_gamma5_rel <<<gpu_gd_blas, gpu_bd_blas >>>(spin3,spin4);
  #else
     dev_gamma5 <<<gpu_gd_blas, gpu_bd_blas >>> (spin3,spin4);     
  #endif
     
     dev_swapmu <<<1,1>>> ();
  #ifdef USETEXTURE
   bind_texture_spin(spin4,1);
  #endif
     //D_tm
     dev_tm_dirac_kappa <<<gpu_gd_M, gpu_bd_M >>> (gf, spin4, spin3, dev_nn);
  
  //Here we have used the output spinor (spinout) to temporarly take the field and to 
  //copy it to the texture field (spin4)!!

  
 //alpha
  dotprod = cublasSdot (24*VOLUME, (const float *) spin2, 1,
            (const float *) spin3, 1);
  alpha = (rk / dotprod); // alpha = r*r/ p M p
   
 //r(k+1)
 cublasSaxpy (24*VOLUME,-1.0*alpha, (const float *) spin3, 1, (float *) spin0, 1);  

 //x(k+1);
 cublasSaxpy (24*VOLUME, alpha, (const float *) spin2,  1, (float *) spin1, 1);

  if((cudaerr=cudaGetLastError()) != cudaSuccess){
    printf("%s\n", cudaGetErrorString(cudaerr));
    exit(200);
  }
  

  //Abbruch?
  dotprod = cublasSdot (24*VOLUME, (const float *) spin0, 1,(const float *) spin0, 1);
  
 if ((dotprod <= eps*sourcesquarenorm)){//error-limit erreicht
   break; 
 }
  printf("iter %d: err = %.8e\n", i, dotprod);
  
 //beta
 beta =dotprod/rk;
 //p(k+1)
 cublasSscal (24*VOLUME, beta, (float *)spin2, 1);
 cublasSaxpy (24*VOLUME, 1.0, (const float *) spin0,  1, (float *) spin2, 1);

 rk = dotprod;
 
 // recalculate residue frome r = b - Ax
 if(((i+1) % N_recalcres) == 0){
    // r_(k+1) = Ax -b 
    printf("Recalculating residue\n");
    
    // D Ddagger   --   Ddagger = gamma5 D gamma5  for Wilson Dirac Operator
    // DO NOT USE tm_dirac_dagger_kappa here, otherwise spin2 will be overwritten!!!
      
      //GAMMA5
    #ifdef USETEXTURE
     unbind_texture_spin(1);
    #endif
    #ifdef RELATIVISTIC_BASIS  
      dev_gamma5_rel <<<gpu_gd_blas, gpu_bd_blas >>> (spin1,spin4);
    #else
      dev_gamma5 <<<gpu_gd_blas, gpu_bd_blas >>> (spin1,spin4);      
    #endif
      dev_swapmu <<<1,1>>> ();
    #ifdef USETEXTURE
     bind_texture_spin(spin4,1);
    #endif
   
      //D_tm GAMMA5, mu -> -mu
      dev_tm_dirac_kappa <<< gpu_gd_M, gpu_bd_M >>> (gf, spin4, spin3, dev_nn);
    #ifdef RELATIVISTIC_BASIS        
      dev_gamma5_rel <<<gpu_gd_blas, gpu_bd_blas >>>(spin3,spinout);
    #else
       dev_gamma5 <<<gpu_gd_blas, gpu_bd_blas >>>(spin3,spinout);     
    #endif
      dev_swapmu <<<1,1>>> ();
  
    //printf("Unbinding texture of spinorfield\n");
    #ifdef USETEXTURE
     unbind_texture_spin(1);
    #endif
    cudaMemcpy(spin4, spinout,size2, cudaMemcpyDeviceToDevice);
    //printf("Rebinding texture to spinorfield\n");
    #ifdef USETEXTURE
     bind_texture_spin(spin4,1);
    #endif
      
      //D_tm
      dev_tm_dirac_kappa<<< gpu_gd_M, gpu_bd_M >>>(gf, spin4, spin3, dev_nn);
    
    // r = b - Ax
    cublasSscal (24*VOLUME, -1.0, (float *)spin3, 1);
    cublasSaxpy (24*VOLUME, 1.0, (const float *) spinin,  1, (float *) spin3, 1);
    cublasScopy (24*VOLUME, (const float *)spin3, 1, (float *)spin0, 1);
    
   }//recalculate residue

 }//MAIN LOOP cg	
  
  
  printf("Final residue: %.6e\n", dotprod);
  // x_result = spin1 !
  
 if(rescalekappa == 1){  //want D^-1 rescaled by 2*kappa
  
//multiply with D^dagger
    #ifdef USETEXTURE
     unbind_texture_spin(1);
    #endif
  #ifdef RELATIVISTIC_BASIS      
      dev_gamma5_rel <<<gpu_gd_blas, gpu_bd_blas >>> (spin1,spin4);
  #else
      dev_gamma5 <<<gpu_gd_blas, gpu_bd_blas >>> (spin1,spin4); 
  #endif
      dev_swapmu <<<1,1>>> ();
    #ifdef USETEXTURE
     bind_texture_spin(spin4,1);
    #endif
      dev_tm_dirac_kappa <<< gpu_gd_M, gpu_bd_M >>> (gf, spin4, spin3, dev_nn);
  #ifdef RELATIVISTIC_BASIS       
      dev_gamma5_rel <<<gpu_gd_blas, gpu_bd_blas >>>(spin3,spin1);
  #else
        dev_gamma5 <<<gpu_gd_blas, gpu_bd_blas >>>(spin3,spin1);    
  #endif
      dev_swapmu <<<1,1>>> ();
    #ifdef USETEXTURE
     unbind_texture_spin(1);
    #endif


 //go over to non-kappa, Ddagger = g5 D g5
 dev_skalarmult_spinor_field<<<gpu_gd_blas, gpu_bd_blas >>>(spin1,1.0/(scaleparam*scaleparam), spinout);  
 // times operator == source ?? 
  #ifdef RELATIVISTIC_BASIS 
   //transform back to tmlqcd gamma basis
   to_tmlqcd_basis<<<gpu_gd_linalg, gpu_bd_linalg>>> (spinout);
  #endif  
  }
  else{
   dev_copy_spinor_field<<<gpu_gd_blas, gpu_bd_blas >>>(spin1,spinout);  
   #ifdef RELATIVISTIC_BASIS 
   //transform back to tmlqcd gamma basis
   to_tmlqcd_basis<<<gpu_gd_linalg, gpu_bd_linalg>>> (spinout);
  #endif   

  }
  
  #ifdef USETEXTURE
   unbind_texture_gf();
  #endif
   
  stop_blas();
  
  return(i);
}



void showspinor(dev_spinor* s){
  int i,j;
  dev_spinor help[6];
  size_t size = 6*sizeof(dev_spinor);
  
  for(i=0; i<VOLUME/2; i++){
    cudaMemcpy(&(help[0]), (s+6*i), size, cudaMemcpyDeviceToHost);
    for(j=0;j<6; j++){
      printf("(%.3f %.3f) (%.3f, %.3f) ", help[j].x, help[j].y, help[j].z, help[j].w);
    }
    printf("\n");
  }
  
}


void check_mixedsolve_params(){
  
    #ifndef LOWOUTPUT
        if (g_cart_id == 0) {
  	  int host_check_VOLUME;
  	  int host_check_Offset;
  	  cudaMemcpyFromSymbol(&host_check_VOLUME, dev_VOLUME, sizeof(int));
  	  cudaMemcpyFromSymbol(&host_check_Offset, dev_Offset, sizeof(int));
  	  printf("\tOn device:\n");
  	  printf("\tdev_VOLUME = %i\n", host_check_VOLUME);
  	  printf("\tdev_Offset = %i\n", host_check_Offset);
	  #ifdef _USE_MPI
	    int host_check_VOLUMEPLUSRAND, host_check_RAND;
	    int host_check_rank, host_check_nproc;	  
	    cudaMemcpyFromSymbol(&host_check_RAND, dev_RAND, sizeof(int));
	    cudaMemcpyFromSymbol(&host_check_VOLUMEPLUSRAND, dev_VOLUMEPLUSRAND, sizeof(int));
	    printf("\tdev_VOLUMEPLUSRAND = %i\n", host_check_VOLUMEPLUSRAND);
	    printf("\tdev_RAND = %i\n", host_check_RAND);
	    cudaMemcpyFromSymbol(&host_check_rank, dev_rank, sizeof(int));
	    cudaMemcpyFromSymbol(&host_check_nproc, dev_nproc, sizeof(int));
	    printf("\tdev_rank = %i\n", host_check_rank);
	    printf("\tdev_nproc = %i\n", host_check_nproc);
	  #endif	  
  	}
  #endif
  
}



#ifndef HALF

// this is the eo version of the device cg inner solver 
// we invert the hermitean Q_{-} Q_{+}
extern "C" int dev_cg_eo(
      dev_su3_2v * gf,
      dev_spinor* spinin, 
      dev_spinor* spinout, 
      dev_spinor* spin0, 
      dev_spinor* spin1, 
      dev_spinor* spin2, 
      dev_spinor* spin3, 
      dev_spinor* spin4, 
      int *grid, int * nn_grid, float epsfinal){
 
 
 float alpha, beta, dotprod, rk, sourcesquarenorm;
 
 int i;
 int maxit = max_innersolver_it;
 float eps = (float) innersolver_precision;
 int N_recalcres = 40; // after N_recalcres iterations calculate r = A x_k - b
 
 cudaError_t cudaerr;
 
  //constants
  update_constants(grid);  
  #ifdef _USE_MPI
    he_cg_init_nd_additional_mpi<<<1,1>>>(VOLUMEPLUSRAND/2, RAND, g_cart_id, g_nproc);
  #endif
  //->check
  check_mixedsolve_params();
  
  #ifdef USETEXTURE
    //Bind texture gf
    bind_texture_gf(gf);
  #endif
 
  start_blas(VOLUME/2);

 //init fields 
 dev_copy_spinor_field<<<gpu_gd_blas, gpu_bd_blas >>>(spinin, spin0);
 dev_zero_spinor_field<<<gpu_gd_blas, gpu_bd_blas >>>(spin1); // x_0 = 0
 dev_copy_spinor_field<<<gpu_gd_blas, gpu_bd_blas >>>(spinin, spin2);
 dev_zero_spinor_field<<<gpu_gd_blas, gpu_bd_blas >>>(spin3);
  

 //relative precision -> get initial residue
 sourcesquarenorm = cublasSdot_wrapper (24*VOLUME/2, (float *)spinin, (float *)spinin);
 rk = sourcesquarenorm; //for use in main loop
 

  if (g_cart_id == 0) {
    #ifndef LOWOUTPUT
      printf("Squarenorm Source:\t%.8e\n", sourcesquarenorm);
      printf("Entering inner solver cg-loop\n");
    #endif
  }
  if ((cudaerr=cudaGetLastError())!=cudaSuccess) {
      printf("%s\n", cudaGetErrorString(cudaGetLastError()));
      exit(200);
  }
    

  #ifdef ALGORITHM_BENCHMARK
    double effectiveflops;	// will used to count the "effective" flop's (from the algorithmic perspective)
    double hoppingflops = 1608.0;
    double matrixflops  = 2*(2*hoppingflops + 120.0 + 132.0); //that is for dev_Qtm_pm_psi
    double starteffective;
    double stopeffective;
   // timer
   starteffective = gettime();
  #endif


 for(i=0;i<maxit;i++){ //MAIN LOOP
  
  // Q_{-}Q{+}
   dev_Qtm_pm_psi(spin2, spin3, gpu_gd_M, gpu_bd_M, gpu_gd_linalg, gpu_bd_linalg);

  // check on Matrix error
  if((cudaerr=cudaGetLastError()) != cudaSuccess){
    if (g_cart_id == 0) printf("%s\n", cudaGetErrorString(cudaerr));
    if (g_cart_id == 0) printf("Error in dev_cg_eo: CUDA error after Matrix application\n", cudaGetErrorString(cudaerr));
    exit(200);
  }
  
  
 //alpha
 dotprod = cublasSdot_wrapper (24*VOLUME/2, (float *) spin2, (float *) spin3);

 alpha = (rk / dotprod); // alpha = r*r/ p M p
   
 //r(k+1)
 cublasSaxpy_wrapper (24*VOLUME/2,-1.0*alpha, (float *) spin3, (float *) spin0);  

 //x(k+1);
 cublasSaxpy_wrapper (24*VOLUME/2, alpha, (float *) spin2, (float *) spin1);

 //break?
 dotprod = cublasSdot_wrapper (24*VOLUME/2, (float *) spin0, (float *) spin0);
 if (((dotprod <= eps*sourcesquarenorm)) || ( dotprod <= epsfinal)){
   break; 
 }
  
 //residue 
 #ifndef LOWOUTPUT 
   if (g_cart_id == 0) printf("iter %d: err = %.8e\n", i, dotprod);
 #endif
  
  
 //beta
 beta = dotprod/rk;
 //p(k+1)
 cublasSscal_wrapper (24*VOLUME/2, beta, (float *)spin2);
 cublasSaxpy_wrapper (24*VOLUME/2, 1.0, (float *) spin0,  (float *) spin2);

 rk = dotprod;
 
 // recalculate residue frome r = b - Ax
 if(((i+1) % N_recalcres) == 0){
    // r_(k+1) = Ax -b 
    #ifndef LOWOUTPUT
    if (g_cart_id == 0) printf("Recalculating residue\n");
    #endif
    // D Ddagger   --   Ddagger = gamma5 D gamma5  for Wilson Dirac Operator
    // DO NOT USE tm_dirac_dagger_kappa here, otherwise spin2 will be overwritten!!!
      
    // Q_{-}Q{+}
     dev_Qtm_pm_psi(spin1, spin3, gpu_gd_M, gpu_bd_M, gpu_gd_linalg, gpu_bd_linalg);

    if((cudaerr=cudaGetLastError()) != cudaSuccess){
      printf("%s\n", cudaGetErrorString(cudaerr));
      exit(200);
    }  
        
    
    // r = b - Ax
      cublasSscal_wrapper (24*VOLUME/2, -1.0, (float *)spin3);
      cublasSaxpy_wrapper (24*VOLUME/2, 1.0, (float *) spinin, (float *) spin3);
      //cublasScopy_wrapper (24*VOLUME/2, (float *)spin3, (float *)spin0);
      dev_copy_spinor_field<<<gpu_gd_blas, gpu_bd_blas >>>(spin3, spin0);

   }//recalculate residue

 }//MAIN LOOP cg	
  
    
    if (g_cart_id == 0) printf("Final residue: %.6e\n", dotprod);
 
  #ifdef ALGORITHM_BENCHMARK
    cudaThreadSynchronize();
    stopeffective = gettime();
      // will now count the number of effective flops
      // effectiveflops  =  #(inner iterations)*(matrixflops+linalgflops)*VOLUME/2  +  #(outer iterations)*(matrixflops+linalgflops)*VOLUME/2
      // outer loop: linalg  =  flops for calculating  r(k+1) and x(k+1)
      // inner loop: linalg  =  flops for calculating  alpha, x(k+1), r(k+1), beta, d(k+1)
     #ifdef _USE_MPI
       int proccount = g_nproc;
     #else
       int proccount = 1;
     #endif 
     if(g_cart_id == 0){
      	effectiveflops = i*proccount*(matrixflops + 2*2*24 + 2*24 + 2*24 + 2*2*24 + 2*24)*VOLUME/2;
      	printf("effective BENCHMARK:\n");
      	printf("\ttotal mixed solver time:   %.4e sec\n", double(stopeffective-starteffective));
      	printf("\tfloating point operations: %.4e flops\n", effectiveflops);
      	printf("\tinner solver performance:  %.4e Gflop/s\n", double(effectiveflops) / double(stopeffective-starteffective) / 1.0e9);
     }
   #endif
  
  // x_result = spin1 !
  
  //no multiplication with D^{dagger} here and no return to non-kappa basis as in dev_cg!
  dev_copy_spinor_field<<<gpu_gd_blas, gpu_bd_blas >>>(spin1,spinout);
  
  
  #ifdef USETEXTURE
   unbind_texture_gf();
  #endif
 
  stop_blas();

  return(i);
}

#endif










//initialize nearest-neighbour table for gpu
void initnn(){
  int t,x,y,z,pos;
  for(t=0;t<T;t++){
   for(x=0; x<LX; x++){
    for(y=0; y<LY; y++){
     for(z=0; z<LZ; z++){   
          pos= z + LZ*(y + LY*(x + LX*t));
          //plus direction
          nn[8*pos+0] = z + LZ*(y + LY*(x + LX*((t+1)%T)));
          nn[8*pos+1] = z + LZ*(y + LY*((x+1)%LX + LX*t));
          nn[8*pos+2] = z + LZ*((y+1)%LY + LY*(x + LX*t));
          nn[8*pos+3] = (z+1)%LZ + LX*(y + LY*(x + LX*t));
          //minus direction
          if(t==0){
            nn[8*pos+4] = z + LZ*(y + LY*(x + LX*((T-1))));
          }
          else{
            nn[8*pos+4] = z + LZ*(y + LY*(x + LX*((t-1))));
          }
          if(x==0){
            nn[8*pos+5] = z + LZ*(y + LY*((LX-1) + LX*t));
          }
          else{
            nn[8*pos+5] = z + LZ*(y + LY*((x-1) + LX*t));
          }
          if(y==0){
            nn[8*pos+6] = z + LZ*((LY-1) + LY*(x + LX*t));
          }
          else{
            nn[8*pos+6] = z + LZ*((y-1) + LY*(x + LX*t));
          }
          if(z==0){
            nn[8*pos+7] = (LZ-1) + LZ*(y + LY*(x + LX*t));
          }
          else{
            nn[8*pos+7] = (z-1) + LZ*(y + LY*(x + LX*t));
          }          
        }
      }
    } 
  }
}





//initialize nearest-neighbour table for gpu with even-odd enabled
//init_nn must have been called before for initialization of nn
void initnn_eo(){
  int x,y,z,t,ind,nnpos,j;
  int evenpos=0;
  int oddpos=0;
  for(t=0;t<T;t++){
    for(x=0;x<LX;x++){
      for(y=0;y<LY;y++){
        for(z=0;z<LZ;z++){
          ind = g_ipt[t][x][y][z];
          
          if(((t+x+y+z)%2 == 0)){
            nnpos = g_lexic2eosub[ind];
            for(j=0;j<4;j++){
              nn_eo[8*nnpos+j] = g_lexic2eosub[ g_iup[ind][j] ];
            }
            for(j=0;j<4;j++){
              nn_eo[8*nnpos+4+j] = g_lexic2eosub[ g_idn[ind][j] ];
            }
            eoidx_even[evenpos] = ind;
            evenpos++;
          }
          else{
            nnpos = g_lexic2eosub[ind];
            for(j=0;j<4;j++){
              nn_oe[8*nnpos+j] = g_lexic2eosub[ g_iup[ind][j] ];
            }
            for(j=0;j<4;j++){
              nn_oe[8*nnpos+4+j] = g_lexic2eosub[ g_idn[ind][j] ];
            }
            eoidx_odd[oddpos] = ind;
            oddpos++;
          }
        }
      }
    }
  }
}




// show the nn table eo
void shownn_eo(){
  int i,pos;
  printf("eo part\n");
  for(pos=0;pos<VOLUME/2;pos++){ 
       printf("p=%d\t", pos);
       for(i=0;i<8;i++){
          printf("%d ",nn_eo[8*pos+i]);
          //lptovec(nn[8*pos+i]);
        }
        printf("\n");
    }
  printf("oe part\n");
  for(pos=0;pos<VOLUME/2;pos++){ 
       printf("p=%d\t", pos);
       for(i=0;i<8;i++){
          printf("%d ",nn_oe[8*pos+i]);
          //lptovec(nn[8*pos+i]);
        }
        printf("\n");
    }
    
  printf("site index even\n");
  for(pos=0;pos<VOLUME/2;pos++){ 
       printf("p=%d\t", pos);
          printf("%d ",eoidx_even[pos]);
          //lptovec(nn[8*pos+i]);
        printf("\n");
  }

  printf("site index odd\n");
  for(pos=0;pos<VOLUME/2;pos++){ 
       printf("p=%d\t", pos);
          printf("%d ",eoidx_odd[pos]);
          //lptovec(nn[8*pos+i]);
        printf("\n");
  }
  printf("checking forward even\n");
  for(pos=0;pos<VOLUME/2;pos++){
    for(i=0;i<4;i++){
      printf("%d = %d\n",pos, nn_oe[8*nn_eo[8*pos+i]+4+i]);
    }
  }

  printf("checking backward even\n");
  for(pos=0;pos<VOLUME/2;pos++){
    for(i=0;i<4;i++){
      printf("%d = %d\n",pos, nn_oe[8*nn_eo[8*pos+4+i]+i]);
    }
  }

  printf("checking forward odd\n");
  for(pos=0;pos<VOLUME/2;pos++){
    for(i=0;i<4;i++){
      printf("%d = %d\n",pos, nn_eo[8*nn_oe[8*pos+i]+4+i]);
    }
  }

  printf("checking backward odd\n");
  for(pos=0;pos<VOLUME/2;pos++){
    for(i=0;i<4;i++){
      printf("%d = %d\n",pos, nn_eo[8*nn_oe[8*pos+4+i]+i]);
    }
  }
}



void lptovec(int k){
  int L3 = L*L*L;
  int L2 = L*L;
  int x0,x1,x2,x3;
  x0 = k/L3;
  k = k-x0*L3; 
  x3 = k/L2;
  k = k-x3*L2;
  x2 = k/L;
  k = k-x2*L;
  x1 = k;
  printf("%d,%d,%d,%d;  ",x0,x3,x2,x1);
}


// show nn table 
void shownn(){
  int t,x,y,z,i,pos;
  int lx,ly,lz,lt;
    lx = LX;
    ly = LY;
    lz = LZ;
    lt =T;  
  for(t=0;t<lt;t++){ 
    for(x=0; x<lx; x++){
      for(y=0; y<ly; y++){
        for(z=0; z<lz; z++){
          pos= z + lz*(y + ly*(x + lx*t));
          printf("p=%d\t", pos);
          for(i=0;i<8;i++){
            printf("%d ",nn[8*pos+i]);
            //lptovec(nn[8*pos+i]);
          }
          printf("\n");
          //compare with geometry fields of hmc
          //might NOT WORK for even-odd? What are geometry indices in case of even-odd?
          printf("%d: %d %d %d %d %d %d %d %d\n",g_ipt[t][x][y][z],g_iup[pos][0],g_iup[pos][1],g_iup[pos][2],g_iup[pos][3],g_idn[pos][0],g_idn[pos][1],g_idn[pos][2],g_idn[pos][3]);
        }
      }
    }
  }
}







extern "C" void init_mixedsolve(su3** gf){
cudaError_t cudaerr;

   // get number of devices
   if(havedevice == 0){
     int ndev = find_devices();
	   if(ndev == 0){
	       fprintf(stderr, "Error: no CUDA devices found. Aborting...\n");
	       exit(300);
	    }
            // only if device_num is not the default (-1)
            if(device_num > -1){ 
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
	      printf("Error in init_mixedsolve(): Could not set active device. Aborting...\n");
	      exit(302);
	    }
           }
           else{
            printf("Not setting any active device. Let the driver choose.\n");
           }        
    havedevice = 1;
    }
  #ifdef GF_8
  /* allocate 8 floats of gf = 2*4*VOLUME float4's*/
  printf("Using GF 8 reconstruction\n");
  size_t dev_gfsize = 2*4*VOLUME * sizeof(dev_su3_8);
  #else
  /* allocate 2 rows of gf = 3*4*VOLUME float4's*/
  printf("Using GF 12 reconstruction\n");
  size_t dev_gfsize = 3*4*VOLUME * sizeof(dev_su3_2v); 
  #endif
  
  #ifdef USETEXTURE
    printf("Using texture references\n");
  #else
    printf("NOT using texture references\n");
  #endif
  if((cudaerr=cudaMalloc((void **) &dev_gf, dev_gfsize)) != cudaSuccess){
    printf("Error in init_mixedsolve(): Memory allocation of gauge field failed. Aborting...\n");
    exit(200);
  }   // Allocate array on device
  else{
    printf("Allocated gauge field on device\n");
  }  
  
  #ifdef GF_8
  h2d_gf = (dev_su3_8 *)malloc(dev_gfsize); // Allocate float conversion gf on host
  su3to8(gf,h2d_gf);  
  #else
  h2d_gf = (dev_su3_2v *)malloc(dev_gfsize); // Allocate float conversion gf on host
  su3to2vf4(gf,h2d_gf);
  #endif
  cudaMemcpy(dev_gf, h2d_gf, dev_gfsize, cudaMemcpyHostToDevice);


//grid 
  size_t nnsize = 8*VOLUME*sizeof(int);
  nn = (int *) malloc(nnsize);
  cudaMalloc((void **) &dev_nn, nnsize);
  
  initnn();
  //shownn();
  //showcompare_gf(T-1, LX-1, LY-1, LZ-1, 3);
  cudaMemcpy(dev_nn, nn, nnsize, cudaMemcpyHostToDevice);
  
  //free again
  free(nn);



// Spinors
  #ifndef HALF
  size_t dev_spinsize = 6*VOLUME * sizeof(dev_spinor); /* float4 */  
  if((void*)(h2d_spin = (dev_spinor *)malloc(dev_spinsize)) == NULL){
    printf("Could not allocate memory for h2d_spin. Aborting...\n");
    exit(200);
  } // Allocate float conversion spinor on host
  #else
  size_t dev_spinsize = 6*VOLUME * sizeof(dev_spinor_half); /*short4*/  
  if((void*)(h2d_spin = (dev_spinor_half *)malloc(dev_spinsize)) == NULL){
    printf("Could not allocate memory for h2d_spin. Aborting...\n");
    exit(200);
  } // Allocate float conversion spinor on host 
  size_t dev_normsize = VOLUME/2 * sizeof(float);
  if((void*)(h2d_spin_norm = (float*)malloc(dev_normsize)) == NULL){
    printf("Could not allocate memory for h2d_spin_norm. Aborting...\n");
    exit(200);
  } // Allocate float conversion norm on host 
  #endif
  
  

  if((cudaerr=cudaGetLastError())!=cudaSuccess){
    printf("Error in init_mixedsolve(): Memory allocation of spinor fields failed. Aborting...\n");
    printf("Error code is: %d\n",cudaerr);    
    exit(200);
  }
  else{
    printf("Allocated spinor fields on device\n");
  }
  

  int grid[6];
  grid[0]=LX; grid[1]=LY; grid[2]=LZ; grid[3]=T; grid[4]=VOLUME; grid[5]=VOLUME;
  // FOR MPI: we have to put grid[5] to VOLUME+RAND !!! 
 
  cudaMalloc((void **) &dev_grid, 6*sizeof(int));
  cudaMemcpy(dev_grid, &(grid[0]), 6*sizeof(int), cudaMemcpyHostToDevice);
  
  if((cudaerr=cudaGetLastError())!=cudaSuccess){
    printf("Error in init_mixedsolve(): grid initialization failed. Aborting...\n");
    exit(200);
  }
  else{
    printf("Allocated grid on device\n");
  }
  //set grid and block sizes, non-eo
  set_gpu_work_layout(0);
}


void update_gpu_gf(su3** gf){
  cudaError_t cudaerr;
  #ifndef LOWOUTPUT
    if(g_cart_id == 0) printf("Updating device single gauge field.\n");
  #endif
  
  #ifndef _USE_MPI
        #ifdef GF_8
          /* allocate 8 floats for gf = 2*4*VOLUME float4's*/
          size_t dev_gsize = 2*4*VOLUME * sizeof(dev_su3_8);
        #else
          /* allocate 2 rows of gf = 3*4*VOLUME float4's*/
          size_t dev_gsize = 3*4*VOLUME * sizeof(dev_su3_2v);
        #endif
  #else
        #ifdef GF_8
          /* allocate 8 floats for gf = 2*4*VOLUME float4's*/
          size_t dev_gsize = 2*4*(VOLUME+RAND) * sizeof(dev_su3_8);
        #else
          /* allocate 2 rows of gf = 3*4*VOLUME float4's*/
          size_t dev_gsize = 3*4*(VOLUME+RAND) * sizeof(dev_su3_2v);
        #endif
  
  #endif  
  #ifdef GF_8
    su3to8(gf,h2d_gf);
  #else
    su3to2vf4(gf,h2d_gf);
  #endif
  //bring to device
  
  if((cudaerr=cudaMemcpy(dev_gf, h2d_gf, dev_gsize, cudaMemcpyHostToDevice))!=cudaSuccess){
          printf("Error in update_gpu_gf(): Could not transfer gf to device. Aborting...\n");
          printf("%s\n", cudaGetErrorString(cudaerr));
	  printf("Error code is: %d\n",cudaerr);
	  exit(200);
   }
}



__global__ void DummyKernel(){}



extern "C" void init_mixedsolve_eo(su3** gf){

  cudaError_t cudaerr;


  if (havedevice == 0) {
  
    // get number of devices
    int ndev = find_devices();
    if(ndev == 0){
      fprintf(stderr, "Error: no CUDA devices found. Aborting...\n");
      exit(300);
    }
    
    // try to set active device to device_num given in input file (or mpi rank)
    #ifndef _USE_MPI
    // only if device_num is not the default (-1)
     if(device_num > -1){ 
    	if(device_num < ndev){
    	  printf("Setting active device to: %d\n", device_num);
    	  //cudaSetDevice(device_num);
    	}
    	else{
   	  fprintf(stderr, "Error: There is no CUDA device with No. %d. Aborting...\n",device_num);
    	  exit(301);
    	}
    	if((cudaerr=cudaGetLastError())!=cudaSuccess){
    	  printf("Error in init_mixedsolve_eo(): Could not set active device. Aborting...\n");
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
    	  if(device_num > -1){ 
            if (device_num < ndev) {
    	      printf("Process %d of %d: Setting active device to: %d\n", g_cart_id, g_nproc, device_num);
    	      cudaSetDevice(device_num);
    	    }
    	    else {
    	      fprintf(stderr, "Process %d of %d: Error: There is no CUDA device with No. %d. Aborting...\n", g_cart_id, g_nproc, device_num);
    	      exit(301);
    	    }
          }
          else{
            printf("Not setting any active device. Let the driver choose.\n");
          } 
  	#else
    	  // device number = mpi rank
          printf("num_gpu_per_node = %d\n", num_gpu_per_node);
	  printf("first device = %d\n", device_num);
    	  if (((g_cart_id + device_num)%num_gpu_per_node ) < ndev) {
    	    printf("Process %d of %d: Setting active device to: %d\n", g_cart_id, g_nproc, (g_cart_id + device_num)%num_gpu_per_node);
    	    cudaSetDevice(((g_cart_id + device_num)%num_gpu_per_node));
	  }
    	  else {
    	    fprintf(stderr, "Process %d of %d: Error: There is no CUDA device with No. %d. Aborting...\n", g_cart_id, g_nproc, g_cart_id);
    	    exit(301);
    	  }
  	#endif
  	if ((cudaerr=cudaGetLastError()) != cudaSuccess) {
  	  printf("Process %d of %d: Error in init_mixedsolve_eo(): Could not set active device. Aborting...\n", g_cart_id, g_nproc);
  	  exit(302);
  	}
    #endif
    
    havedevice=1;
    
  }
//set cache configuration
  //if (g_cart_id == 0) printf("Setting GPU cache configuration to prefer L1 cache\n");
  //cudaFuncSetCacheConfig(DummyKernel, cudaFuncCachePreferL1);
//  
  
  
  
  // output
  #ifdef _USE_MPI
    if (g_cart_id == 0) {
  #endif
    #ifndef LOWOUTPUT
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
    #endif
  #ifdef _USE_MPI
    }
  #endif
  #ifndef LOWOUTPUT
    #ifdef RELATIVISTIC_BASIS
      if (g_cart_id == 0) printf("Using RELATIVISTIC gamma basis.\n");
    #else
      if (g_cart_id == 0) printf("Using TMLQCD gamma basis.\n");    
    #endif
  #endif
  #ifndef _USE_MPI
  	#ifdef GF_8
  	  /* allocate 8 floats for gf = 2*4*VOLUME float4's*/
  	  size_t dev_gfsize = 2*4*VOLUME * sizeof(dev_su3_8);
  	#else
  	  /* allocate 2 rows of gf = 3*4*VOLUME float4's*/
  	  size_t dev_gfsize = 3*4*VOLUME * sizeof(dev_su3_2v);
  	#endif
  #else
  	#ifdef GF_8
  	  /* allocate 8 floats for gf = 2*4*VOLUME float4's*/
  	  size_t dev_gfsize = 2*4*(VOLUME+RAND) * sizeof(dev_su3_8);
  	#else
  	  /* allocate 2 rows of gf = 3*4*VOLUME float4's*/
  	  size_t dev_gfsize = 3*4*(VOLUME+RAND) * sizeof(dev_su3_2v);
  	#endif
  
  #endif
  
  if((cudaerr=cudaMalloc((void **) &dev_gf, dev_gfsize)) != cudaSuccess){
    printf("Error in init_mixedsolve_eo(): Memory allocation of gauge field failed. Aborting...\n");
    exit(200);
  }   // Allocate array on device
  else {
    #ifndef _USE_MPI
     #ifndef LOWOUTPUT
      printf("Allocated memory for gauge field on device.\n");
     #endif
    #else
      if (g_cart_id == 0) printf("Allocated memory for gauge field on devices.\n");
    #endif
  }
  
  #ifdef GF_8
    h2d_gf = (dev_su3_8 *)malloc(dev_gfsize); // Allocate float conversion gf on host
    su3to8(gf,h2d_gf);
  #else
    h2d_gf = (dev_su3_2v *)malloc(dev_gfsize); // Allocate float conversion gf on host
    su3to2vf4(gf,h2d_gf);
  #endif
  //bring to device
  cudaMemcpy(dev_gf, h2d_gf, dev_gfsize, cudaMemcpyHostToDevice);
  
  
  #ifdef HALF
    #ifndef _USE_MPI
      #ifdef GF_8
       /* allocate 8 floats for gf = 2*4*VOLUME float4's*/
        printf("Using half precision GF 8 reconstruction\n");
       dev_gfsize = 2*4*VOLUME * sizeof(dev_su3_8_half); 
      #else
        /* allocate 2 rows of gf = 3*4*VOLUME float4's*/
        printf("Using half precision GF 12 reconstruction\n");
        dev_gfsize = 3*4*VOLUME * sizeof(dev_su3_2v_half); 
      #endif  
    #else // MPI
      #ifdef GF_8
       /* allocate 8 floats for gf = 2*4*VOLUME float4's*/
        printf("Using half precision GF 8 reconstruction\n");
       dev_gfsize = 2*4*(VOLUME+RAND) * sizeof(dev_su3_8_half); 
      #else
        /* allocate 2 rows of gf = 3*4*VOLUME float4's*/
        printf("Using half precision GF 12 reconstruction\n");
        dev_gfsize = 3*4*(VOLUME+RAND) * sizeof(dev_su3_2v_half); 
      #endif      
    #endif //MPI
    if((cudaerr=cudaMalloc((void **) &dev_gf_half, dev_gfsize)) != cudaSuccess){
    printf("Error in init_mixedsolve_eo(): Memory allocation of half precsion gauge field failed. Aborting...\n");
    exit(200);
    }   // Allocate array on device
    else{
      printf("Allocated half precision gauge field on device\n");
    }      
     
  #endif // HALF


//grid 
  size_t nnsize = 8*VOLUME*sizeof(int);
  nn = (int *) malloc(nnsize);
  
  //nn grid for even-odd
  nn_eo = (int *) malloc(nnsize/2);
  nn_oe = (int *) malloc(nnsize/2);
  
  cudaMalloc((void **) &dev_nn, nnsize);
  cudaMalloc((void **) &dev_nn_eo, nnsize/2);
  cudaMalloc((void **) &dev_nn_oe, nnsize/2);
  
  #ifndef _USE_MPI
    size_t idxsize = VOLUME/2*sizeof(int);
  #else
    size_t idxsize = (VOLUME+RAND)/2*sizeof(int);
  #endif
  eoidx_even = (int *) malloc(idxsize);
  eoidx_odd = (int *) malloc(idxsize);
  cudaMalloc((void **) &dev_eoidx_even, idxsize);
  cudaMalloc((void **) &dev_eoidx_odd, idxsize);
  
  #ifndef _USE_MPI
    initnn();
    initnn_eo();
    //shownn_eo();
  #else
    init_nnspinor_eo_mpi();
    init_idxgauge_mpi();

/*
    char filename[50];
    sprintf(filename, "nnfield_proc%d", g_cart_id);
    FILE * outfile = fopen(filename,"w");
    int m1;
     for(m1=0; m1<VOLUME/2; m1++){

        fprintf(outfile,"%d %d %d %d %d %d %d %d\n",nn_eo[8*m1],
                       nn_eo[8*m1+1],nn_eo[8*m1+2],nn_eo[8*m1+3], 
                       nn_eo[8*m1+4],nn_eo[8*m1+5],nn_eo[8*m1+6],
                       nn_eo[8*m1+7] );        
         fprintf(outfile,"%d\n",eoidx_even[m1]);                              
      }
   fclose(outfile);
*/
 
  #endif
  
  //shownn();
  //showcompare_gf(T-1, LX-1, LY-1, LZ-1, 3);
  //check_gauge_reconstruction_8(gf, dev_gf, 0, 0);
  cudaMemcpy(dev_nn, nn, nnsize, cudaMemcpyHostToDevice);
  cudaMemcpy(dev_nn_eo, nn_eo, nnsize/2, cudaMemcpyHostToDevice);
  cudaMemcpy(dev_nn_oe, nn_oe, nnsize/2, cudaMemcpyHostToDevice);
  cudaMemcpy(dev_eoidx_even, eoidx_even, idxsize, cudaMemcpyHostToDevice);
  cudaMemcpy(dev_eoidx_odd, eoidx_odd, idxsize, cudaMemcpyHostToDevice);

  if((cudaerr=cudaGetLastError())!=cudaSuccess){
    if(g_cart_id==0) printf("Error in init_mixedsolve_eo(): NN-fields transfer failed. Aborting...\n");
    exit(200);
  }
  else{
    if(g_cart_id==0) printf("Allocated NN-fields on device\n");
  }  
  //free again
  free(eoidx_odd);
  free(eoidx_even);
  free(nn_oe);
  free(nn_eo);
  free(nn);
  
  
// Spinors
  #ifndef HALF

	#ifdef _USE_MPI
  	  size_t dev_spinsize_ext =  6*(VOLUME+RAND)/2*sizeof(dev_spinor);
	  if((void*)(h2d_spin = (dev_spinor *)malloc(dev_spinsize_ext)) == NULL){
  	    printf("Could not allocate memory for h2d_spin. Aborting...\n");
  	    exit(200);
  	  } // Allocate float conversion spinor on host
	#else  	
  	  size_t dev_spinsize = 6*VOLUME/2 * sizeof(dev_spinor); /* float4 */	
	  if((void*)(h2d_spin = (dev_spinor *)malloc(dev_spinsize)) == NULL){
  	    printf("Could not allocate memory for h2d_spin. Aborting...\n");
  	    exit(200);
  	  } // Allocate float conversion spinor on host	
        #endif
  #else
  	size_t dev_spinsize = 6*VOLUME/2 * sizeof(dev_spinor_half);/*short4*/
  	if((void*)(h2d_spin = (dev_spinor_half *)malloc(dev_spinsize)) == NULL){
  	  printf("Could not allocate memory for h2d_spin. Aborting...\n");
  	  exit(200);
  	} // Allocate float conversion spinor on host 
  	
  	size_t dev_spinsize2 = 6*VOLUME/2 * sizeof(dev_spinor);/*short4*/
  	if((void*)(h2d_spin_reliable = (dev_spinor*)malloc(dev_spinsize2)) == NULL){
  	  printf("Could not allocate memory for h2d_spin_reliable. Aborting...\n");
  	  exit(200);
  	} // Allocate float conversion spinor on host for reliable update
  	
  	size_t dev_normsize = VOLUME/2 * sizeof(float);
  	if((void*)(h2d_spin_norm = (float *)malloc(dev_normsize)) == NULL){
  	  printf("Could not allocate memory for h2d_spin_norm. Aborting...\n");
  	  exit(200);
  	} // Allocate float conversion norm on host 

  #endif
  

#ifdef _USE_MPI
   	size_t dev_spinsize_d = 12*(VOLUME+RAND)/2 * sizeof(dev_spinor_d); /* double2 */
#else
   	size_t dev_spinsize_d = 12*VOLUME/2 * sizeof(dev_spinor_d); /* double2 */  
#endif
	//dev_spin_eo1/2_d we need for gpu_deriv_SB in any case	
  	cudaMalloc((void **) &dev_spin_eo1_d, dev_spinsize_d);	
  	cudaMalloc((void **) &dev_spin_eo2_d, dev_spinsize_d);   
	if((void*)(h2d_spin_d = (dev_spinor_d *)malloc(dev_spinsize_d)) == NULL){
  	    printf("Could not allocate memory for double h2d_spin_d. Aborting...\n");
  	    exit(200);
  	} 
  	  
	#ifdef GPU_DOUBLE   
	  cudaMalloc((void **) &dev_spin0_d, dev_spinsize_d);   
	  cudaMalloc((void **) &dev_spin1_d, dev_spinsize_d);	
  	  cudaMalloc((void **) &dev_spin2_d, dev_spinsize_d);   
  	  cudaMalloc((void **) &dev_spin3_d, dev_spinsize_d);


  	  if((cudaerr=cudaGetLastError())!=cudaSuccess){
              if(g_cart_id==0) printf("Error in init_mixedsolve_eo(): Memory allocation of double spinor fields failed. Aborting...\n");
              exit(200);
          }
        #endif  
        
#ifndef _USE_MPI

	
       #ifdef HALF

  	cudaMalloc((void **) &dev_spin_eo1_half, dev_spinsize);
  	cudaMalloc((void **) &dev_spin_eo2_half, dev_spinsize);
        cudaMalloc((void **) &dev_spin_eo1_half_norm, dev_normsize);
        cudaMalloc((void **) &dev_spin_eo2_half_norm, dev_normsize);
	
	// allocate 3 float fields for reliable update
	// as well as eo fields for the float matrix
	cudaMalloc((void **) &dev_spin1_reliable, dev_spinsize2);  
  	cudaMalloc((void **) &dev_spin2_reliable, dev_spinsize2);
	cudaMalloc((void **) &dev_spin3_reliable, dev_spinsize2);
	cudaMalloc((void **) &dev_spin_eo1, dev_spinsize2);
  	cudaMalloc((void **) &dev_spin_eo2, dev_spinsize2);
      #else
        cudaMalloc((void **) &dev_spin_eo1, dev_spinsize);
  	cudaMalloc((void **) &dev_spin_eo2, dev_spinsize);
      #endif
  
  
#else  //now comes mpi

  	#ifdef HALF
  	  dev_spinsize_ext =  6*(VOLUME+RAND)/2*sizeof(dev_spinor_half);
  	  size_t dev_normsize_ext =  (VOLUME+RAND)/2*sizeof(float);
  	#else
  	  dev_spinsize_ext =  6*(VOLUME+RAND)/2*sizeof(dev_spinor);
	#endif


  	
        #ifdef HALF

         cudaMalloc((void **) &dev_spin_eo1_half, dev_spinsize_ext);
  	 cudaMalloc((void **) &dev_spin_eo2_half, dev_spinsize_ext);
         cudaMalloc((void **) &dev_spin_eo1_half_norm, dev_normsize_ext);
         cudaMalloc((void **) &dev_spin_eo2_half_norm, dev_normsize_ext);
	
	
	// allocate 3 float fields for reliable update
	// as well as eo fields for the float matrix
	dev_spinsize_ext2 =  6*(VOLUME+RAND)/2*sizeof(dev_spinor);
	cudaMalloc((void **) &dev_spin1_reliable, dev_spinsize_ext2);  
  	cudaMalloc((void **) &dev_spin2_reliable, dev_spinsize_ext2);
	cudaMalloc((void **) &dev_spin3_reliable, dev_spinsize_ext2);
	
	cudaMalloc((void **) &dev_spin_eo1, dev_spinsize_ext);
  	cudaMalloc((void **) &dev_spin_eo2, dev_spinsize_ext);
       #else
        cudaMalloc((void **) &dev_spin_eo1, dev_spinsize_ext);
  	cudaMalloc((void **) &dev_spin_eo2, dev_spinsize_ext);
      #endif   	
      
      int tSliceEO = LX*LY*LZ/2;
      #ifndef HALF
  	R1 = (dev_spinor *) malloc(2*tSliceEO*24*sizeof(float));
  	R2 = R1 + 6*tSliceEO;
  	R3 = (dev_spinor *) malloc(2*tSliceEO*24*sizeof(float));
  	R4 = R3 + 6*tSliceEO;
      #else
      
  	// implement this for half?
   
      #endif

//for gathering and spreading of indizes of rand in (gather_rand spread_rand called from xchange_field_wrapper)
    #ifdef RELATIVISTIC_BASIS
      cudaMalloc((void **) &RAND_FW, tSliceEO*3*sizeof(float4));
      cudaMalloc((void **) &RAND_BW, tSliceEO*3*sizeof(float4));
    #else
      cudaMalloc((void **) &RAND_FW, tSliceEO*6*sizeof(float4));
      cudaMalloc((void **) &RAND_BW, tSliceEO*6*sizeof(float4));      
    #endif
    /*  for async communication */
    // page-locked memory    
   #ifndef HALF 
    #ifdef RELATIVISTIC_BASIS
      int flperspin = 3;
    #else
      int flperspin = 6;
    #endif
    
//     cudaMallocHost(&RAND3, 2*tSliceEO*flperspin*sizeof(float4));
//     RAND4 = RAND3 + flperspin*tSliceEO;
//     cudaMallocHost(&RAND1, 2*tSliceEO*flperspin*sizeof(float4));
//     RAND2 = RAND1 + flperspin*tSliceEO;
   
    cudaMallocHost(&RAND3, tSliceEO*flperspin*sizeof(float4));
    cudaMallocHost(&RAND4, tSliceEO*flperspin*sizeof(float4));
    cudaMallocHost(&RAND1, tSliceEO*flperspin*sizeof(float4));
    cudaMallocHost(&RAND2, tSliceEO*flperspin*sizeof(float4));
  #else
    cudaMallocHost(&RAND3, 2*tSliceEO*6*sizeof(short4));
    RAND4 = RAND3 + 6*tSliceEO;
    cudaMallocHost(&RAND1, 2*tSliceEO*6*sizeof(short4));
    RAND2 = RAND1 + 6*tSliceEO;
    //norm page-locked mem
    cudaMallocHost(&RAND3_norm, 2*tSliceEO*sizeof(float));
    RAND4_norm = RAND3_norm + tSliceEO;
    cudaMallocHost(&RAND1_norm, 2*tSliceEO*sizeof(float));
    RAND2_norm = RAND1_norm + tSliceEO;
  #endif  
  //HALF

#ifdef GPU_DOUBLE
  R1_D = (dev_spinor_d *) malloc(2*tSliceEO*24*sizeof(double));
  R2_D = R1_D + 12*tSliceEO;
  R3_D = (dev_spinor_d *) malloc(2*tSliceEO*24*sizeof(double));
  R4_D = R3_D + 12*tSliceEO;


//for gathering and spreading of indizes of rand in (gather_rand spread_rand called from xchange_field_wrapper)
    #ifdef RELATIVISTIC_BASIS
      cudaMalloc((void **) &RAND_FW_D, tSliceEO*6*sizeof(double2));
      cudaMalloc((void **) &RAND_BW_D, tSliceEO*6*sizeof(double2));
    #else
      cudaMalloc((void **) &RAND_FW_D, tSliceEO*12*sizeof(double2));
      cudaMalloc((void **) &RAND_BW_D, tSliceEO*12*sizeof(double2));      
    #endif
    /*  for async communication */
    // page-locked memory    
    #ifdef RELATIVISTIC_BASIS
      int dbperspin = 6;
    #else
      int dbperspin = 12;
    #endif  
    cudaMallocHost(&RAND3_D, tSliceEO*dbperspin*sizeof(double2));
    cudaMallocHost(&RAND4_D, tSliceEO*dbperspin*sizeof(double2));
    cudaMallocHost(&RAND1_D, tSliceEO*dbperspin*sizeof(double2));
    cudaMallocHost(&RAND2_D, tSliceEO*dbperspin*sizeof(double2));
    

#endif    
    
    // CUDA streams and events
    for (int i = 0; i < 2*nStreams+1; i++) {
        cudaStreamCreate(&stream[i]);
    }    
    /* end for async communication */  	

#endif 
//MPI


  if((cudaerr=cudaGetLastError())!=cudaSuccess){
    if(g_cart_id==0) printf("Error in init_mixedsolve_eo(): Memory allocation of spinor fields failed. Aborting...\n");
    exit(200);
  }
  else{
    if(g_cart_id==0) printf("Allocated spinor fields on device\n");
  }
  

 
  if((cudaerr=cudaPeekAtLastError()) != cudaSuccess){
    printf("Error in init_mixedsolve_eo: %s\n", cudaGetErrorString(cudaerr));
    exit(200);
  }  
  

  int grid[6];
  grid[0]=LX; grid[1]=LY; grid[2]=LZ; grid[3]=T; grid[4]=VOLUME/2; 
  // dev_VOLUME is half of VOLUME for eo
  
  // put dev_Offset accordingly depending on mpi/non-mpi
  #ifdef _USE_MPI
   grid[5] = (VOLUME+RAND)/2;
  #else
   grid[5] = VOLUME/2;
  #endif
  cudaMalloc((void **) &dev_grid, 6*sizeof(int));
  cudaMemcpy(dev_grid, &(grid[0]), 6*sizeof(int), cudaMemcpyHostToDevice);
  
  
  /*
  init_dev_observables();
 
  clock_t start, stop; 
  double timeelapsed = 0.0;
  int count;
  
  assert((start = clock())!=-1);
  float devplaq;
  //for(count=0; count<1; count++){
    devplaq = calc_plaquette(dev_gf, dev_nn);
  //}
  assert((stop = clock())!=-1);
  timeelapsed = (double) (stop-start)/CLOCKS_PER_SEC;
  printf("Calculating Plaquette on device: plaq(device) = %.8f\n", devplaq);
  printf("Time spent calculating: %f sec\n", timeelapsed);
  
  assert((start = clock())!=-1);
  float hostplaq;
  int a = 0;
  //for(count=0; count<1; count++){
    hostplaq = (float) measure_gauge_action()/(6.*VOLUME*g_nproc);
  //}
  assert((stop = clock())!=-1);
  timeelapsed = (double) (stop-start)/CLOCKS_PER_SEC;
  printf("Calculating Plaquette on host: plaq(host) = %.8f\n", hostplaq);
  printf("Time spent calculating: %f sec\n", timeelapsed);

  float devrect;
  assert((start = clock())!=-1);
  //for(count=0; count<100; count++){
    devrect = calc_rectangle(dev_gf, dev_nn);
  //}
  assert((stop = clock())!=-1);
  timeelapsed = (double) (stop-start)/CLOCKS_PER_SEC;
  printf("Calculating Rectangles on device: rectangle(device) = %.8f\n", devrect);
  printf("Time spent calculating: %f sec\n", timeelapsed);
  
  float hostrect;
  assert((start = clock())!=-1);
  //for(count=0; count<100; count++){
    hostrect = (float) measure_rectangles()/(12.*VOLUME*g_nproc);
  //}
  assert((stop = clock())!=-1);
  timeelapsed = (double) (stop-start)/CLOCKS_PER_SEC;
  printf("Calculating Rectangles on host: rectangle(host) = %.8f\n", hostrect);
  printf("Time spent calculating: %f sec\n", timeelapsed);
 
 
  float2 ret;

  calc_polyakov_0(&ret, dev_gf, dev_nn);
  printf("Calculating Polyakov loop on device:\n");  
  printf("pl_0 (Re) = %.8e\n",ret.x);
  printf("pl_0 (Im) = %.8e\n",ret.y);
  
  //polyakov_loop_dir(1, 0);
  //printf("Calculating Polyakov loop on host:\n");  
 
  finalize_dev_observables();

  exit(100);
  */
  
  
  //set grid and block sizes, eo!
  set_gpu_work_layout(1);
  
  if((cudaerr=cudaGetLastError())!=cudaSuccess){
    if(g_cart_id==0) printf("Error in init_mixedsolve_eo(): Something went wrong. CUDA error at end of function. Aborting...\n");
    exit(200);
  }
  else{
    if(g_cart_id==0) printf("Finished init_mixedsolve_eo()\n");
  }
  
}



extern "C" void finalize_mixedsolve(){


  cudaFree(dev_gf);
  cudaFree(dev_grid);
  cudaFree(dev_nn);
  
  if(even_odd_flag){
    cudaFree(dev_spin_eo1);
    cudaFree(dev_spin_eo2);
    cudaFree(dev_eoidx_even);
    cudaFree(dev_eoidx_odd);
    cudaFree(dev_nn_eo);
    cudaFree(dev_nn_oe);
  }
  
  #ifdef HALF
    cudaFree(dev_gf_half);
 
    
    //reliable update fields
    cudaFree(dev_spin1_reliable);
    cudaFree(dev_spin2_reliable);
    cudaFree(dev_spin3_reliable);
    free(h2d_spin_reliable);
    
    if(even_odd_flag){
     cudaFree(dev_spin_eo1_half);
     cudaFree(dev_spin_eo2_half); 
     cudaFree(dev_spin_eo1_half_norm);
     cudaFree(dev_spin_eo2_half_norm);
    }
        
  #endif
  
#ifdef _USE_MPI
  cudaFreeHost(RAND1);
  cudaFreeHost(RAND2);  
  cudaFreeHost(RAND3);
  cudaFreeHost(RAND4);  
  cudaFree(RAND_BW);
  cudaFree(RAND_FW);
  free(R1);
  free(R3);
#ifdef GPU_DOUBLE
  cudaFreeHost(RAND1_D);
  cudaFreeHost(RAND2_D); 
  cudaFreeHost(RAND3_D);
  cudaFreeHost(RAND4_D);  
  cudaFree(RAND_BW_D);
  cudaFree(RAND_FW_D);
  free(R1_D);
  free(R3_D);  
#endif
  #ifdef HALF
   cudaFreeHost(RAND1_norm);
   cudaFreeHost(RAND2_norm);   
   cudaFreeHost(RAND3_norm);
   cudaFreeHost(RAND4_norm);   
  #endif
             
  for (int i = 0; i < 2*nStreams+1; i++) {
     cudaStreamDestroy(stream[i]);
  } 
#endif 
  
  
  free(h2d_spin);
  free(h2d_gf);
}


//allocates cg single solver fields
void init_mixedsolve_fields(int eo){
  int N;
  size_t dev_spinsize; 
  cudaError_t cudaerr;
  
#ifdef HALF
  size_t dev_normsize;
#endif
  
  if(eo){
   #ifdef _USE_MPI
     N = (VOLUME+RAND)/2;
   #else
     N = VOLUME/2;
   #endif
  }
  else{
   #ifdef _USE_MPI
     N = (VOLUME+RAND);
   #else
     N = VOLUME;
   #endif  
  }

#ifndef LOWOUTPUT
  if(g_cart_id ==0) printf("Allocating mixed solver single fields on device.\n");
#endif

  #ifdef _USE_MPI
    #ifdef HALF
      dev_spinsize =  6* N * sizeof(dev_spinor_half);
      dev_normsize =  N *sizeof(float);
    #else
      dev_spinsize =  6* N * sizeof(dev_spinor);    
    #endif
  #else
    #ifdef HALF
      dev_spinsize = 6* N * sizeof(dev_spinor_half);  
      dev_normsize =  N *sizeof(float);      
    #else
      dev_spinsize = 6* N * sizeof(dev_spinor);
    #endif
  #endif  
	  
  cudaMalloc((void **) &dev_spin1, dev_spinsize);   // Allocate array spin1 on device
  cudaMalloc((void **) &dev_spin2, dev_spinsize);   // Allocate array spin2 on device
  cudaMalloc((void **) &dev_spin3, dev_spinsize);   // Allocate array spin3 on device
  cudaMalloc((void **) &dev_spin4, dev_spinsize);
  cudaMalloc((void **) &dev_spin5, dev_spinsize);
  cudaMalloc((void **) &dev_spinin, dev_spinsize);
  cudaMalloc((void **) &dev_spinout, dev_spinsize); 
  
  #ifdef HALF
    cudaMalloc((void **) &dev_spin1_norm, dev_normsize);   // Allocate norm spin1 on device
    cudaMalloc((void **) &dev_spin2_norm, dev_normsize);   // Allocate norm spin2 on device
    cudaMalloc((void **) &dev_spin3_norm, dev_normsize);   // Allocate norm spin3 on device
    cudaMalloc((void **) &dev_spin4_norm, dev_normsize);
    cudaMalloc((void **) &dev_spin5_norm, dev_normsize);
    cudaMalloc((void **) &dev_spinin_norm, dev_normsize);
    cudaMalloc((void **) &dev_spinout_norm, dev_normsize);
  #endif
	
  if((cudaerr=cudaGetLastError())!=cudaSuccess){
    if(g_cart_id==0){
      printf("Error in init_mixedsolve_fields(): could not allocate fields\n");
      printf("Error was %d. Aborting...\n", cudaerr);
    }
    exit(200);
  }	
}



//frees cg single solver fields
void finalize_mixedsolve_fields(){

  cudaFree(dev_spin1);
  cudaFree(dev_spin2);
  cudaFree(dev_spin3);
  cudaFree(dev_spin4);
  cudaFree(dev_spin5);
  cudaFree(dev_spinin);
  cudaFree(dev_spinout);
  
#ifdef HALF
    cudaFree(dev_spin1_norm);
    cudaFree(dev_spin2_norm);
    cudaFree(dev_spin3_norm);
    cudaFree(dev_spin4_norm);
    cudaFree(dev_spin5_norm);
    cudaFree(dev_spinin_norm);
    cudaFree(dev_spinout_norm);
#endif
    
}


// include half versions of dev_cg - solvers
#ifdef HALF
  #include "half_solvers.cuh"
#endif





#ifndef HALF




void benchmark(spinor * const Q){
  
  double timeelapsed = 0.0;
  clock_t start, stop;
  int i;
  
  int ibench;
  #ifdef OPERATOR_BENCHMARK
   ibench = OPERATOR_BENCHMARK;
  #else
    ibench = 100;
  #endif
  
  size_t dev_spinsize = 6*VOLUME*sizeof(dev_spinor); // float4 even-odd !
  convert2REAL4_spin(Q,h2d_spin);
  cudaMemcpy(dev_spinin, h2d_spin, dev_spinsize, cudaMemcpyHostToDevice);
  printf("%s\n", cudaGetErrorString(cudaGetLastError()));
  
  #ifndef _USE_MPI
    assert((start = clock())!=-1);
  #else
    start = MPI_Wtime();
  #endif  
  
 

 #ifdef USETEXTURE
  //Bind texture gf
  bind_texture_gf(dev_gf);
 #endif

 //Initialize some stuff
  printf("mu = %f\n", g_mu);
  dev_complex h0,h1,h2,h3;
  h0.re = (float)creal(ka0);    h0.im = -(float)cimag(ka0);
  h1.re = (float)creal(ka1);    h1.im = -(float)cimag(ka1);
  h2.re = (float)creal(ka2);    h2.im = -(float)cimag(ka2);
  h3.re = (float)creal(ka3);    h3.im = -(float)cimag(ka3);
  
  
  he_cg_init<<< 1, 1 >>> (dev_grid, (float) g_kappa, (float)(g_mu/(2.0*g_kappa)), h0,h1,h2,h3); 
  printf("%s\n", cudaGetErrorString(cudaGetLastError()));
  printf("Applying H %d times\n", ibench);
  for(i=0; i<ibench; i++){
  
    #ifdef USETEXTURE
         bind_texture_spin(dev_spinin,1);
    #endif
     dev_tm_dirac_kappa <<<gpu_gd_M, gpu_bd_M >>> (dev_gf, dev_spinout, dev_spinin, dev_nn);
    #ifdef USETEXTURE             
      unbind_texture_spin(1);
    #endif


  }  
  printf("%s\n", cudaGetErrorString(cudaGetLastError())); 
  printf("Done\n"); 
  
  cudaThreadSynchronize();

    assert((stop = clock())!=-1);
    timeelapsed = (double) (stop-start)/CLOCKS_PER_SEC;
    double benchres = 1608.0*(VOLUME)* ibench / timeelapsed / 1.0e9;
    printf("Elapsed time was: %f sec\n", timeelapsed); 
    printf("Benchmark: %f Gflops\n", benchres); 


  #ifdef USETEXTURE
    unbind_texture_gf();
  #endif
}








extern "C" int mixed_solve (spinor * const P, spinor * const Q, const int max_iter, 
	   double eps, const int rel_prec,const int N){
  
  // source in Q, initial solution in P (not yet implemented)
  double rk;
  int outercount=0;
  int totalcount=0;
  clock_t start, stop, startinner, stopinner; 
  double timeelapsed = 0.0;
  double sourcesquarenorm;
  int iter;
  spinor ** solver_field = NULL;
  const int nr_sf = 4;
  init_solver_field(&solver_field, VOLUMEPLUSRAND, nr_sf);  

  size_t dev_spinsize = 6*VOLUME * sizeof(dev_spinor); // float4 
  //update the gpu single gauge_field
  update_gpu_gf(g_gauge_field); 
  
  //allocate solver fields non-eo!
  init_mixedsolve_fields(0);
  
  // Start timer
  assert((start = clock())!=-1);
  
  rk = square_norm(Q, N, 0);
  sourcesquarenorm = rk; // for relative precision
  assign(solver_field[0],Q,N);
  printf("Initial residue: %.16e\n",rk);
  zero_spinor_field(solver_field[1],  N);//spin2 = x_k
  zero_spinor_field(solver_field[2],  N);
  printf("The VOLUME is: %d\n",N);
  
  #ifdef OPERATOR_BENCHMARK
    // small benchmark
      assign(solver_field[0],Q,N);
        benchmark(solver_field[0]);
    // end small benchmark
  #endif
 
  
for(iter=0; iter<max_iter; iter++){

   printf("Applying double precision Dirac-Op...\n");
   
   Q_pm_psi_gpu(solver_field[3], solver_field[2]);
   diff(solver_field[0],solver_field[0],solver_field[3],N);
    // r_k = b - D x_k
   
   rk = square_norm(solver_field[0], N, 0);
  
   #ifdef GF_8
    if(isnan(rk)){
      fprintf(stderr, "Error in mixed_solve: Residue is NaN.\n  May happen with GF 8 reconstruction. Aborting ...\n");
      exit(200);
    }
   #endif
   
   printf("Residue after %d inner solver iterations: %.18e\n",outercount,rk);
   if(((rk <= eps) && (rel_prec == 0)) || ((rk <= eps*sourcesquarenorm) && (rel_prec == 1)))
   {
     printf("Reached solver precision of eps=%.2e\n",eps);
     //multiply with D^dagger
     Q_minus_psi_gpu(solver_field[3], solver_field[1]);
     assign(P, solver_field[3], N);
  

    stop = clock();
    timeelapsed = (double) (stop-start)/CLOCKS_PER_SEC;
    printf("Inversion done in mixed precision.\n Number of iterations in outer solver: %d\n Squared residue: %.8e\n Time elapsed: %.6e sec\n", outercount, rk, timeelapsed);
    finalize_solver(solver_field, nr_sf);
    finalize_mixedsolve_fields();
    return(totalcount);  
   }
   

  //initialize spin fields on device
  convert2REAL4_spin(solver_field[0],h2d_spin);
  
  cudaMemcpy(dev_spinin, h2d_spin, dev_spinsize, cudaMemcpyHostToDevice);
  printf("%s\n", cudaGetErrorString(cudaGetLastError()));

   // solve in single prec on device
   // D p_k = r_k
   printf("Entering inner solver\n");
   assert((startinner = clock())!=-1);
   totalcount += dev_cg(dev_gf, dev_spinin, dev_spinout, dev_spin1, dev_spin2, dev_spin3, dev_spin4, dev_spin5, dev_grid,dev_nn, 0);
   stopinner = clock();
   timeelapsed = (double) (stopinner-startinner)/CLOCKS_PER_SEC;
   printf("Inner solver done\nTime elapsed: %.6e sec\n", timeelapsed);
   
  
   // copy back
   cudaMemcpy(h2d_spin, dev_spinout, dev_spinsize, cudaMemcpyDeviceToHost);
   printf("%s\n", cudaGetErrorString(cudaGetLastError()));
   
   convert2double_spin(h2d_spin, solver_field[2]);
   
   add(solver_field[1],solver_field[1],solver_field[2],N);
   // x_(k+1) = x_k + p_k
   
   outercount ++;
    
}// outer loop 

     printf("Did NOT reach solver precision of eps=%.2e\n",eps);
     //multiply with D^dagger
     Q_minus_psi_gpu(solver_field[3], solver_field[1]);
     assign(P, solver_field[3], N);
  

    stop = clock();
    timeelapsed = (double) (stop-start)/CLOCKS_PER_SEC;
    printf("Inversion done in mixed precision.\n Number of iterations in outer solver: %d\n Squared residue: %.8e\n Time elapsed: %.6e sec\n", outercount, rk, timeelapsed);
    finalize_solver(solver_field, nr_sf);
    finalize_mixedsolve_fields();    
  return(-1);
}






void benchmark_eo(spinor * const Q){
  
  double timeelapsed = 0.0;
  clock_t start, stop;
  int i;
  
  int ibench;
  #ifdef OPERATOR_BENCHMARK
   ibench = OPERATOR_BENCHMARK;
  #else
    ibench = 100;
  #endif
  
  
  size_t dev_spinsize = 6*VOLUME/2 * sizeof(dev_spinor); // float4 even-odd !
  convert2REAL4_spin(Q,h2d_spin);
  cudaMemcpy(dev_spinin, h2d_spin, dev_spinsize, cudaMemcpyHostToDevice);
  printf("%s\n", cudaGetErrorString(cudaGetLastError()));
  
  #ifndef _USE_MPI
    assert((start = clock())!=-1);
  #else
    start = MPI_Wtime();
  #endif  
  
 
  int VolumeEO = VOLUME/2;
  

 #ifdef USETEXTURE
  //Bind texture gf
  bind_texture_gf(dev_gf);
 #endif

 //Initialize some stuff
  printf("mu = %f\n", g_mu);
  dev_complex h0,h1,h2,h3;
  h0.re = (float)creal(ka0);    h0.im = -(float)cimag(ka0);
  h1.re = (float)creal(ka1);    h1.im = -(float)cimag(ka1);
  h2.re = (float)creal(ka2);    h2.im = -(float)cimag(ka2);
  h3.re = (float)creal(ka3);    h3.im = -(float)cimag(ka3);
  

  
  he_cg_init<<< 1, 1 >>> (dev_grid, (float) g_kappa, (float)(g_mu/(2.0*g_kappa)), h0,h1,h2,h3); 
  printf("%s\n", cudaGetErrorString(cudaGetLastError()));
  printf("Applying H %d times\n", ibench);
  for(i=0; i<ibench; i++){
  
      #ifdef _USE_MPI
           xchange_field_wrapper(dev_spinin, 0);
      #endif
      #ifdef USETEXTURE
         bind_texture_spin(dev_spinin,1);
      #endif
       //bind_texture_nn(dev_nn_eo);
      dev_Hopping_Matrix<<<gpu_gd_M, gpu_bd_M >>>
             (dev_gf, dev_spinin, dev_spin_eo1, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0, 0, VolumeEO); //dev_spin_eo1 == even -> 0
       //unbind_texture_nn();
    #ifdef USETEXTURE             
      unbind_texture_spin(1);
    #endif

    #ifdef _USE_MPI
        xchange_field_wrapper(dev_spin_eo1, 0);
    #endif
       bind_texture_spin(dev_spin_eo1,1);
  //bind_texture_nn(dev_nn_oe);
    dev_Hopping_Matrix<<<gpu_gd_M, gpu_bd_M >>>
            (dev_gf, dev_spin_eo1, dev_spinin, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1, 0, VolumeEO); 
  //unbind_texture_nn();
    #ifdef USETEXTURE
      unbind_texture_spin(1);
   #endif
  

  }  
  printf("%s\n", cudaGetErrorString(cudaGetLastError())); 
  printf("Done\n"); 
  
  cudaThreadSynchronize();

  #ifndef _USE_MPI
    assert((stop = clock())!=-1);
    timeelapsed = (double) (stop-start)/CLOCKS_PER_SEC;
    // x2 because 2x Hopping per iteration
    double benchres = 1608.0*2*(VOLUME/2)* ibench / timeelapsed / 1.0e9;
    printf("Elapsed time was: %f sec\n", timeelapsed); 
    printf("Benchmark: %f Gflops\n", benchres); 
  #else
    stop = MPI_Wtime();
    timeelapsed = (double) (stop-start);
    // x2 because 2x Hopping per iteration
    double benchres = 1608.0*2*(g_nproc*VOLUME/2)* ibench / timeelapsed / 1.0e9;
    if (g_cart_id == 0) {
      printf("Benchmark: %f Gflops\n", benchres); 
    }
  #endif  
  
   
   
  #ifdef USETEXTURE
    unbind_texture_gf();
  #endif
}



#ifdef _USE_MPI
void benchmark_eo_mpi(spinor * const Q){
  
  double timeelapsed = 0.0;
  clock_t start, stop;
  int i;
  int ibench;
  #ifdef OPERATOR_BENCHMARK
   ibench = OPERATOR_BENCHMARK;
  #else
    ibench = 100;
  #endif
 
  size_t dev_spinsize = 6*VOLUME/2 * sizeof(dev_spinor); // float4 even-odd !
  convert2REAL4_spin(Q,h2d_spin);
  cudaMemcpy(dev_spinin, h2d_spin, dev_spinsize, cudaMemcpyHostToDevice);
  printf("%s\n", cudaGetErrorString(cudaGetLastError()));
  
  #ifndef _USE_MPI
    assert((start = clock())!=-1);
  #else
    start = MPI_Wtime();
  #endif  
  
  
  

 #ifdef USETEXTURE
  //Bind texture gf
  bind_texture_gf(dev_gf);
 #endif

 //Initialize some stuff
  printf("mu = %f\n", g_mu);
  dev_complex h0,h1,h2,h3;
  h0.re = (float)creal(ka0);    h0.im = -(float)cimag(ka0);
  h1.re = (float)creal(ka1);    h1.im = -(float)cimag(ka1);
  h2.re = (float)creal(ka2);    h2.im = -(float)cimag(ka2);
  h3.re = (float)creal(ka3);    h3.im = -(float)cimag(ka3);
  
  
  he_cg_init<<< 1, 1 >>> (dev_grid, (float) g_kappa, (float)(g_mu/(2.0*g_kappa)), h0,h1,h2,h3); 
  printf("%s\n", cudaGetErrorString(cudaGetLastError()));
  printf("Applying dev_Qtm_pm_psi %d times\n",ibench);
  
  for(i=0; i<ibench; i++){
  
      
   dev_Qtm_pm_psi_mpi(dev_spinin, dev_spin_eo1, gpu_gd_M, gpu_bd_M, gpu_gd_linalg, gpu_bd_linalg);   
   
   dev_Qtm_pm_psi_mpi(dev_spin_eo1, dev_spinin, gpu_gd_M, gpu_bd_M, gpu_gd_linalg, gpu_bd_linalg); 

  }  
  printf("%s\n", cudaGetErrorString(cudaGetLastError())); 
  printf("Done\n"); 
  
  
  
  #ifndef _USE_MPI
    assert((stop = clock())!=-1);
    timeelapsed = (double) (stop-start)/CLOCKS_PER_SEC;
    // x8 because 8x Hopping per iteration
    double benchres = 1608.0*8*(VOLUME/2)* ibench / timeelapsed / 1.0e9;
    printf("Benchmark: %f Gflops\n", benchres); 
  #else
    stop = MPI_Wtime();
    timeelapsed = (double) (stop-start);
    // 8 because 8x Hopping per iteration
    double benchres = 1608.0*8*(g_nproc*VOLUME/2)* ibench / timeelapsed / 1.0e9;
    if (g_cart_id == 0) {
      printf("Benchmark: %f Gflops\n", benchres); 
    }
  #endif  
  
  
  #ifdef USETEXTURE
    unbind_texture_gf();
  #endif
}

#endif









#else // with half we end up here upon calling mixed_solve FIXME
extern "C" int mixed_solve (spinor * const P, spinor * const Q, const int max_iter, 
           double eps, const int rel_prec,const int N){
   printf("WARNING dummy function mixed_solve called\n");
   return(0);           
}

#endif






void test_double_operator(spinor* const Q, const int N){
   
   size_t dev_spinsize_d;
#ifdef _USE_MPI
  dev_spinsize_d  = 12*(VOLUME+RAND)/2 * sizeof(dev_spinor_d); // double2 even-odd !
#else
  dev_spinsize_d  = 12*VOLUME/2 * sizeof(dev_spinor_d); 
#endif   
 
        
  spinor ** solver_field = NULL;
  const int nr_sf = 3;
  init_solver_field(&solver_field, VOLUMEPLUSRAND/2, nr_sf);  

  //apply cpu matrix
  Qtm_pm_psi(solver_field[0], Q);  
  /*
  Hopping_Matrix(EO,solver_field[1], Q);
  Hopping_Matrix(OE, solver_field[0] , solver_field_up[1]); 
  */
  //apply gpu matrix
  order_spin_gpu(Q, h2d_spin_d);
  cudaMemcpy(dev_spin0_d, h2d_spin_d, dev_spinsize_d, cudaMemcpyHostToDevice);
  #ifdef RELATIVISTIC_BASIS
     to_relativistic_basis_d<<<gpu_gd_linalg_d, gpu_bd_linalg_d>>> (dev_spin0_d);
  #endif  
  dev_Qtm_pm_psi_d(dev_spin0_d, dev_spin1_d,  
		      dev_spin_eo1_d, dev_spin_eo2_d, 
		      gpu_gd_M_d, gpu_bd_M_d, gpu_gd_linalg_d, gpu_bd_linalg_d,
		      dev_eoidx_even, dev_eoidx_odd, 
		      dev_nn_eo, dev_nn_oe); 
  
  #ifdef RELATIVISTIC_BASIS
     to_tmlqcd_basis_d<<<gpu_gd_linalg_d, gpu_bd_linalg_d>>> (dev_spin1_d);
  #endif  		
  cudaMemcpy(h2d_spin_d, dev_spin1_d, dev_spinsize_d, cudaMemcpyDeviceToHost);
  unorder_spin_gpu(h2d_spin_d, solver_field[1]);      

  diff(solver_field[2], solver_field[1], solver_field[0],N);
  double rk = square_norm(solver_field[2], N, 1);
    
  printf("Testing double matrix:\n");
  printf("cpu: Squared difference is: %.8e\n", rk);
  printf("cpu: Squared difference per spinor component is: %.8e\n", rk/N/24.0);  
  
  //now test dev_diff...
  order_spin_gpu(solver_field[0], h2d_spin_d);
  cudaMemcpy(dev_spin0_d, h2d_spin_d, dev_spinsize_d, cudaMemcpyHostToDevice);


  dev_diff_d<<<gpu_gd_blas_d, gpu_bd_blas_d>>>(dev_spin2_d,dev_spin0_d,dev_spin1_d);
  double rk_gpu = double_dotprod(dev_spin2_d,dev_spin2_d, N);
  
  printf("Testing double linalg:\n");
  printf("gpu: Squared difference is: %.8e\n", rk_gpu);
  printf("gpu: Squared difference per spinor component is: %.8e\n", rk_gpu/N/24.0); 
  
  //fetch back gpu diff 
  cudaMemcpy(h2d_spin_d, dev_spin2_d, dev_spinsize_d, cudaMemcpyDeviceToHost);
  unorder_spin_gpu(h2d_spin_d, solver_field[0]); 
  
  diff(solver_field[1], solver_field[2], solver_field[0],N);
  rk = square_norm(solver_field[1], N, 1);
  printf("gpu <-> cpu diff operation is: %.8e\n", rk); 
  
  
  finalize_solver(solver_field, nr_sf);  
}










extern "C" int mixed_solve_eo (spinor * const P, spinor * const Q, const int max_iter, 
	   double eps, const int rel_prec, const int N){

  // source in Q, initial solution in P (not yet implemented)
  double rk;
  int outercount=0;
  int totalcount=0;
  clock_t start, stop, startinner, stopinner; 
  double timeelapsed = 0.0;
  double sourcesquarenorm;
  int iter;
  cudaError_t cudaerr;
  #ifndef GPU_DOUBLE
    size_t dev_spinsize;
    
    #ifndef HALF
      #ifndef _USE_MPI
	dev_spinsize = 6*VOLUME/2 * sizeof(dev_spinor); // float4 even-odd !
      #else
	dev_spinsize = 6*VOLUMEPLUSRAND/2 * sizeof(dev_spinor); // float4 even-odd !
      #endif
    #else
    #ifndef _USE_MPI
      dev_spinsize = 6*VOLUME/2 * sizeof(dev_spinor_half); //short4 eo !
      size_t dev_normsize = VOLUME/2 * sizeof(float);
    #else
      dev_spinsize = 6*VOLUMEPLUSRAND/2 * sizeof(dev_spinor_half); //short4 eo !
      size_t dev_normsize = VOLUMEPLUSRAND/2 * sizeof(float);   
    #endif
    #endif  
  #endif
  
  //update the gpu single gauge_field
  update_gpu_gf(g_gauge_field);
  
  //allocate solver fields eo!
  init_mixedsolve_fields(1);

  #ifdef GPU_DOUBLE
    size_t dev_spinsize_d; 
    #ifdef _USE_MPI
      dev_spinsize_d = 12*(VOLUME+RAND)/2 * sizeof(dev_spinor_d); // double2 even-odd ! 
    #else
      dev_spinsize_d = 12*VOLUME/2 * sizeof(dev_spinor_d); // double2 even-odd !  
    #endif
  
    update_constants_d(dev_grid);
    update_gpu_gf_d(g_gauge_field);
  #endif 
  
  //initialize solver fields 
  spinor ** solver_field = NULL;
  const int nr_sf = 4;
  init_solver_field(&solver_field, VOLUMEPLUSRAND/2, nr_sf);  

  
  #ifdef OPERATOR_BENCHMARK
    #ifndef HALF
    // small benchmark
      assign(solver_field[0],Q,N);
      #ifndef _USE_MPI
        benchmark_eo(solver_field[0]);
      #else
        benchmark_eo_mpi(solver_field[0]); 
      #endif
    // end small benchmark
    #endif //not HALF
  #endif
 

  // Start timer
  assert((start = clock())!=-1);
  rk = square_norm(Q, N, 1);
  sourcesquarenorm=rk; // for relative prec
  double finaleps;
  if(rel_prec == 1){
    finaleps = eps * sourcesquarenorm;
  }
  else{
    finaleps = eps;
  }


  
  #ifdef GPU_DOUBLE
    #ifdef MATRIX_DEBUG
      test_double_operator(Q,N);
    #endif

    order_spin_gpu(Q, h2d_spin_d);
    cudaMemcpy(dev_spin0_d, h2d_spin_d, dev_spinsize_d, cudaMemcpyHostToDevice);

    order_spin_gpu(P, h2d_spin_d);
    cudaMemcpy(dev_spin1_d, h2d_spin_d, dev_spinsize_d, cudaMemcpyHostToDevice); 
    cudaMemcpy(dev_spin2_d, h2d_spin_d, dev_spinsize_d, cudaMemcpyHostToDevice);     
    cudaMemcpy(dev_spin3_d, h2d_spin_d, dev_spinsize_d, cudaMemcpyHostToDevice); 
    #ifdef RELATIVISTIC_BASIS
        to_relativistic_basis_d<<<gpu_gd_linalg_d, gpu_bd_linalg_d>>> (dev_spin0_d);
        to_relativistic_basis_d<<<gpu_gd_linalg_d, gpu_bd_linalg_d>>> (dev_spin1_d);   
        to_relativistic_basis_d<<<gpu_gd_linalg_d, gpu_bd_linalg_d>>> (dev_spin2_d);
        to_relativistic_basis_d<<<gpu_gd_linalg_d, gpu_bd_linalg_d>>> (dev_spin3_d);	
    #endif 
  #else  
    assign(solver_field[0],Q,N);
    zero_spinor_field(solver_field[1],  N);//spin2 = x_k
    zero_spinor_field(solver_field[2],  N);
  #endif
    
  #ifndef LOWOUTPUT
    if(g_cart_id==0) printf("Initial residue: %.16e\n",rk);
    if(g_cart_id==0) printf("The VOLUME/2 is: %d\n",N);
  #endif


for(iter=0; iter<max_iter; iter++){
   #ifndef LOWOUTPUT
   if(g_cart_id==0) printf("Applying double precision EO Dirac-Op Q_{-}Q{+}...\n");
   #endif
     // r_k = b - D x_k 
   #ifdef GPU_DOUBLE
     dev_Qtm_pm_psi_d(dev_spin2_d, dev_spin3_d,  
		      dev_spin_eo1_d, dev_spin_eo2_d, 
		      gpu_gd_M_d, gpu_bd_M_d, gpu_gd_linalg_d, gpu_bd_linalg_d,
		      dev_eoidx_even, dev_eoidx_odd, 
		      dev_nn_eo, dev_nn_oe); 
     dev_diff_d<<<gpu_gd_blas_d, gpu_bd_blas_d>>>(dev_spin0_d,dev_spin0_d,dev_spin3_d);
     rk = double_dotprod(dev_spin0_d,dev_spin0_d, N);
   #else   
     Qtm_pm_psi(solver_field[3], solver_field[2]);
     diff(solver_field[0],solver_field[0],solver_field[3],N);

     rk = square_norm(solver_field[0], N, 1);
   #endif
   
   #ifdef GF_8
    if(isnan(rk)){
      fprintf(stderr, "Error in mixed_solve_eo: Residue is NaN.\n  May happen with GF 8 reconstruction. Aborting ...\n");
      exit(200);
    }
   #endif
   
   if(g_cart_id==0) printf("Residue after %d inner solver iterations: %.18e\n",outercount,rk);
   
   if(((rk <= eps) && (rel_prec == 0)) || ((rk <= eps*sourcesquarenorm) && (rel_prec == 1)))
   {
/*     

    #ifdef GPU_DOUBLE
      dev_Qtm_minus_psi_d(dev_spin1_d, dev_spin2_d,  
		      dev_spin_eo1_d, 
		      griddim3, blockdim3, griddim4, blockdim4,
		      dev_eoidx_even, dev_eoidx_odd, 
		      dev_nn_eo, dev_nn_oe);       
      cudaMemcpy(h2d_spin_d, dev_spin2_d, dev_spinsize_d, cudaMemcpyDeviceToHost);
      unorder_spin_gpu(h2d_spin_d, solver_field[3]);  
      assign(P, solver_field[3], N);     
    #else
     //multiply with Qtm_minus_psi (for non gpu done in invert_eo.c)
     Qtm_minus_psi(solver_field[3], solver_field[1]);
     assign(P, solver_field[3], N);
    #endif

*/
    #ifdef GPU_DOUBLE
      #ifdef RELATIVISTIC_BASIS
        to_tmlqcd_basis_d<<<gpu_gd_linalg_d, gpu_bd_linalg_d>>> (dev_spin1_d);
      #endif
      cudaMemcpy(h2d_spin_d, dev_spin1_d, dev_spinsize_d, cudaMemcpyDeviceToHost);
      unorder_spin_gpu(h2d_spin_d, solver_field[1]);
    #endif
     //multiply with Qtm_minus_psi (for non gpu done in invert_eo.c)
     Qtm_minus_psi(solver_field[3], solver_field[1]);
     assign(P, solver_field[3], N);

     stop = clock();
     timeelapsed = (double) (stop-start)/CLOCKS_PER_SEC;
     
     #ifndef LOWOUTPUT
      if(g_cart_id==0) printf("Reached solver precision of eps=%.2e\n",eps);
      if(g_cart_id==0) printf("EO Inversion done in mixed precision.\n Number of iterations in outer solver: %d\n Squared residue: %.8e\n Time elapsed: %.6e sec\n", outercount, rk, timeelapsed);
     #endif
     finalize_solver(solver_field, nr_sf);
     finalize_mixedsolve_fields();
     return(totalcount);  
   }
   
  //initialize spin fields on device
  #ifndef HALF
    #ifdef GPU_DOUBLE
      dev_d2f<<<gpu_gd_blas_d, gpu_bd_blas_d>>>(dev_spinin, dev_spin0_d);
      //cudaThreadSynchronize();
      if ((cudaerr=cudaPeekAtLastError())!=cudaSuccess) {
        if (g_cart_id == 0) printf("Error in mixedsolve_eo_gpu: %s\n", cudaGetErrorString(cudaerr));
        exit(100);
      }      
    #else   
      convert2REAL4_spin(solver_field[0],h2d_spin);
      cudaMemcpy(dev_spinin, h2d_spin, dev_spinsize, cudaMemcpyHostToDevice); 
      #ifdef RELATIVISTIC_BASIS
        to_relativistic_basis<<<gpu_gd_linalg, gpu_bd_linalg>>> (dev_spinin);      
      #endif
    #endif
  #else
    convert2REAL4_spin_half(solver_field[0],h2d_spin, h2d_spin_norm); 
    cudaMemcpy(dev_spinin, h2d_spin, dev_spinsize, cudaMemcpyHostToDevice); 
  // also copy half spinor norm and reliable (float) source
    cudaMemcpy(dev_spinin_norm, h2d_spin_norm, dev_normsize, cudaMemcpyHostToDevice);
    convert2REAL4_spin(solver_field[0],h2d_spin_reliable);
    cudaMemcpy(dev_spin1_reliable, h2d_spin_reliable, dev_spinsize_reliable, cudaMemcpyHostToDevice); 
  #endif

  
  
  if ((cudaerr=cudaGetLastError())!=cudaSuccess) {
    printf("%s\n", cudaGetErrorString(cudaerr));
  }
   // solve in single prec on device
   // D p_k = r_k
   #ifndef LOWOUTPUT
     if(g_cart_id==0) printf("Entering inner solver\n");
   #endif
   assert((startinner = clock())!=-1);
   #ifndef HALF
      totalcount += dev_cg_eo(dev_gf, dev_spinin, dev_spinout, dev_spin1, dev_spin2, dev_spin3, dev_spin4, dev_spin5, dev_grid,dev_nn, (float) finaleps);
   #else
     
     totalcount += dev_cg_eo_half(dev_gf, 
                 dev_spinin, dev_spinin_norm,
                 dev_spinout,dev_spinout_norm,
                 dev_spin1, dev_spin1_norm,
                 dev_spin2, dev_spin2_norm,
                 dev_spin3, dev_spin3_norm,
                 dev_spin4, dev_spin4_norm,
                 dev_spin5, dev_spin5_norm,
                 dev_grid,dev_nn, (float) finaleps); 


		 
   #endif
   stopinner = clock();
   timeelapsed = (double) (stopinner-startinner)/CLOCKS_PER_SEC;
   
   #ifndef LOWOUTPUT
   if(g_cart_id==0) printf("Inner solver done\nTime elapsed: %.6e sec\n", timeelapsed);
   #endif
   // copy back
   #ifdef GPU_DOUBLE
     dev_add_f2d<<<gpu_gd_blas_d, gpu_bd_blas_d >>>(dev_spin1_d,dev_spin1_d,dev_spinout);
     dev_f2d<<<gpu_gd_blas_d, gpu_bd_blas_d >>>(dev_spin2_d,dev_spinout);
   #else 
     #ifdef RELATIVISTIC_BASIS
        to_tmlqcd_basis<<<gpu_gd_linalg, gpu_bd_linalg>>> (dev_spinout);
     #endif
     cudaMemcpy(h2d_spin, dev_spinout, dev_spinsize, cudaMemcpyDeviceToHost);
     #ifdef HALF
      cudaMemcpy(h2d_spin_norm, dev_spinout_norm, dev_normsize, cudaMemcpyDeviceToHost);
     #endif
     if ((cudaerr=cudaGetLastError())!=cudaSuccess) {
       printf("%s\n", cudaGetErrorString(cudaerr));
       printf("Error code is: %d\n",cudaerr);       
     }
     #ifndef HALF
       convert2double_spin(h2d_spin, solver_field[2]);
     #else
       convert2double_spin_half(h2d_spin, h2d_spin_norm, solver_field[2]);
     #endif
   
     // x_(k+1) = x_k + p_k
     add(solver_field[1],solver_field[1],solver_field[2],N);
   #endif
   outercount ++;   
}// outer loop 
    
     if(g_cart_id==0) printf("Did NOT reach solver precision of eps=%.2e\n",eps);
     //multiply with Qtm_minus_psi (for non gpu done in invert_eo.c)
   #ifdef GPU_DOUBLE
      #ifdef RELATIVISTIC_BASIS
        to_tmlqcd_basis_d<<<gpu_gd_linalg_d, gpu_bd_linalg_d>>> (dev_spin1_d);
      #endif    
      cudaMemcpy(h2d_spin_d, dev_spin1_d, dev_spinsize_d, cudaMemcpyDeviceToHost);
      unorder_spin_gpu(h2d_spin_d, solver_field[1]);   
   #endif      
     Qtm_minus_psi(solver_field[3], solver_field[1]);
     assign(P, solver_field[3], N);
    

    assert((stop = clock())!=-1);
    timeelapsed = (double) (stop-start)/CLOCKS_PER_SEC;
    if(g_cart_id==0) printf("Inversion done in mixed precision.\n Number of iterations in outer solver: %d\n Squared residue: %.8e\n Time elapsed: %.6e sec\n", outercount, rk, timeelapsed);

    finalize_solver(solver_field, nr_sf);
    finalize_mixedsolve_fields();
  return(-1);
}









extern "C" int linsolve_eo_gpu (spinor * const P, spinor * const Q, const int max_iter, 
	   double eps, const int rel_prec, const int N){

  // source in Q, initial solution in P
  double rk;
  int outercount=0;
  int totalcount=0;
  clock_t start, stop, startinner, stopinner; 
  double timeelapsed = 0.0;
  double sourcesquarenorm;
  int iter;
  cudaError_t cudaerr;
  
  #ifndef GPU_DOUBLE
    size_t dev_spinsize;
    #ifndef HALF
      #ifndef _USE_MPI
	dev_spinsize = 6*VOLUME/2 * sizeof(dev_spinor); // float4 even-odd !     
      #else
	dev_spinsize = 6*VOLUMEPLUSRAND/2 * sizeof(dev_spinor); // float4 even-odd !
      #endif
    #else
    #ifndef _USE_MPI
      dev_spinsize = 6*VOLUME/2 * sizeof(dev_spinor_half); //short4 eo !
      size_t dev_normsize = VOLUME/2 * sizeof(float);
    #else
      dev_spinsize = 6*VOLUMEPLUSRAND/2 * sizeof(dev_spinor_half); //short4 eo !
      size_t dev_normsize = VOLUMEPLUSRAND/2 * sizeof(float);   
    #endif
    #endif
  #endif  
  //update the gpu single gauge_field
  update_gpu_gf(g_gauge_field);
  //allocate solver fields eo!
  init_mixedsolve_fields(1);
  
  #ifdef GPU_DOUBLE
    size_t dev_spinsize_d;   
    #ifdef _USE_MPI
      dev_spinsize_d = 12*(VOLUME+RAND)/2 * sizeof(dev_spinor_d); // double2 even-odd ! 
    #else
      dev_spinsize_d = 12*VOLUME/2 * sizeof(dev_spinor_d); // double2 even-odd !  
    #endif    
      int grid[6];
      grid[0]=LX; grid[1]=LY; grid[2]=LZ; grid[3]=T; grid[4]=VOLUME/2; 
	// dev_VOLUME is half of VOLUME for eo
	
      // put dev_Offset accordingly depending on mpi/non-mpi
      #ifdef _USE_MPI
      grid[5] = (VOLUME+RAND)/2;
      #else
      grid[5] = VOLUME/2;
      #endif
      cudaMemcpy(dev_grid, &(grid[0]), 6*sizeof(int), cudaMemcpyHostToDevice);    
      
      update_constants_d(dev_grid);
      update_gpu_gf_d(g_gauge_field);
  #endif
 
  // Start timer
  assert((start = clock())!=-1);
  rk = square_norm(Q, N, 1);
  sourcesquarenorm=rk; // for relative prec
  double finaleps;
  if(rel_prec == 1){
    finaleps = eps * sourcesquarenorm;
  }
  else{
    finaleps = eps;
  }

  spinor ** solver_field = NULL;
  const int nr_sf = 4;
  init_solver_field(&solver_field, VOLUMEPLUSRAND/2, nr_sf);  

  #ifdef GPU_DOUBLE
    
    
    #ifdef MATRIX_DEBUG
      test_double_operator(Q,N);
    #endif    
      
    order_spin_gpu(Q, h2d_spin_d);
    cudaMemcpy(dev_spin0_d, h2d_spin_d, dev_spinsize_d, cudaMemcpyHostToDevice);
    //cudaThreadSynchronize();
    

      order_spin_gpu(P, h2d_spin_d);
      cudaMemcpy(dev_spin1_d, h2d_spin_d, dev_spinsize_d, cudaMemcpyHostToDevice); 
      cudaMemcpy(dev_spin2_d, h2d_spin_d, dev_spinsize_d, cudaMemcpyHostToDevice);     
      cudaMemcpy(dev_spin3_d, h2d_spin_d, dev_spinsize_d, cudaMemcpyHostToDevice);
      #ifdef RELATIVISTIC_BASIS
	  to_relativistic_basis_d<<<gpu_gd_linalg_d, gpu_bd_linalg_d>>> (dev_spin0_d);
	  to_relativistic_basis_d<<<gpu_gd_linalg_d, gpu_bd_linalg_d>>> (dev_spin1_d);   
	  to_relativistic_basis_d<<<gpu_gd_linalg_d, gpu_bd_linalg_d>>> (dev_spin2_d);
	  to_relativistic_basis_d<<<gpu_gd_linalg_d, gpu_bd_linalg_d>>> (dev_spin3_d);	
      #endif     
    #else
    assign(solver_field[0],Q,N);
    /* set initial guess*/
    assign(solver_field[1],P,N);    
    assign(solver_field[2],P,N);
  #endif
  
  #ifndef LOWOUTPUT
    if (g_cart_id == 0) printf("Initial residue: %.16e\n",rk);
    if (g_cart_id == 0) printf("The VOLUME/2 is: %d\n",N);
  #endif


for(iter=0; iter<max_iter; iter++){
   #ifndef LOWOUTPUT
   if (g_cart_id == 0) printf("Applying double precision EO Dirac-Op Q_{-}Q{+}...\n");
   #endif
   
   #ifdef GPU_DOUBLE
     dev_Qtm_pm_psi_d(dev_spin2_d, dev_spin3_d,  
		      dev_spin_eo1_d, dev_spin_eo2_d, 
		      gpu_gd_M_d, gpu_bd_M_d, gpu_gd_linalg_d, gpu_bd_linalg_d,
		      dev_eoidx_even, dev_eoidx_odd, 
		      dev_nn_eo, dev_nn_oe); 
     dev_diff_d<<<gpu_gd_blas_d, gpu_bd_blas_d >>>(dev_spin0_d,dev_spin0_d,dev_spin3_d);
     rk = double_dotprod(dev_spin0_d,dev_spin0_d, N);
   #else
     Qtm_pm_psi(solver_field[3], solver_field[2]);
     diff(solver_field[0],solver_field[0],solver_field[3],N);
      // r_k = b - D x_k  
     rk = square_norm(solver_field[0], N, 1);
   #endif
   
   #ifdef GF_8
    if(isnan(rk)){
      if (g_cart_id == 0) fprintf(stderr, "Error in linsolve_eo_gpu: Residue is NaN.\n  May happen with GF 8 reconstruction. Aborting ...\n");
      exit(200);
    }
   #endif
   
    if (g_cart_id == 0) printf("Residue after %d inner solver iterations: %.18e\n",outercount,rk);
   
   if(((rk <= eps) && (rel_prec == 0)) || ((rk <= eps*sourcesquarenorm) && (rel_prec == 1)))
   {
    #ifdef GPU_DOUBLE
      #ifdef RELATIVISTIC_BASIS
        to_tmlqcd_basis_d<<<gpu_gd_linalg_d, gpu_bd_linalg_d>>> (dev_spin1_d);
      #endif
      cudaMemcpy(h2d_spin_d, dev_spin1_d, dev_spinsize_d, cudaMemcpyDeviceToHost);
      unorder_spin_gpu(h2d_spin_d, P);    
    #else  
     assign(P, solver_field[1], N);
    #endif
     stop = clock();
     timeelapsed = (double) (stop-start)/CLOCKS_PER_SEC;
     
     #ifndef LOWOUTPUT
      if (g_cart_id == 0) printf("Reached solver precision of eps=%.2e\n",eps);
      if (g_cart_id == 0) printf("EO Inversion done in mixed precision.\n Number of iterations in outer solver: %d\n Squared residue: %.8e\n Time elapsed: %.6e sec\n", outercount, rk, timeelapsed);
     #endif
     if (g_cart_id == 0) printf("Number of inner (single) solver iterations: %d\n",totalcount);
     finalize_solver(solver_field, nr_sf);
     finalize_mixedsolve_fields();
     return(totalcount);  
   }
   
  //initialize spin fields on device
  #ifndef HALF
    #ifdef GPU_DOUBLE
      dev_d2f<<<gpu_gd_blas_d, gpu_bd_blas_d >>>(dev_spinin, dev_spin0_d);
      //cudaThreadSynchronize();
      if ((cudaerr=cudaPeekAtLastError())!=cudaSuccess) {
        if (g_cart_id == 0) printf("Error in linsolve_eo_gpu: %s\n", cudaGetErrorString(cudaGetLastError()));
        if (g_cart_id == 0) printf("Error code is: %d\n",cudaerr);
	exit(100);
      }      
    #else 
      convert2REAL4_spin(solver_field[0],h2d_spin);
      cudaMemcpy(dev_spinin, h2d_spin, dev_spinsize, cudaMemcpyHostToDevice);    
      #ifdef RELATIVISTIC_BASIS
        to_relativistic_basis<<<gpu_gd_linalg, gpu_bd_linalg>>> (dev_spinin);      
      #endif      
    #endif
  #else
    convert2REAL4_spin_half(solver_field[0],h2d_spin, h2d_spin_norm);
    cudaMemcpy(dev_spinin, h2d_spin, dev_spinsize, cudaMemcpyHostToDevice);    
    // also copy half spinor norm   
      cudaMemcpy(dev_spinin_norm, h2d_spin_norm, dev_normsize, cudaMemcpyHostToDevice);
  #endif

  if ((cudaerr=cudaGetLastError())!=cudaSuccess) {
    if (g_cart_id == 0) printf("%s\n", cudaGetErrorString(cudaGetLastError()));
    if (g_cart_id == 0) printf("Error code is: %d\n",cudaerr);   
  }
   // solve in single prec on device
   // D p_k = r_k
   #ifndef LOWOUTPUT
     if (g_cart_id == 0) printf("Entering inner solver\n");
   #endif
   assert((startinner = clock())!=-1);
   #ifndef HALF
      totalcount += dev_cg_eo(dev_gf, 
			      dev_spinin, dev_spinout, 
			      dev_spin1, dev_spin2, 
			      dev_spin3, dev_spin4, 
			      dev_spin5, dev_grid,
			      dev_nn, (float) finaleps);
   #else
     totalcount += dev_cg_eo_half(dev_gf, 
                 dev_spinin, dev_spinin_norm,
                 dev_spinout,dev_spinout_norm,
                 dev_spin1, dev_spin1_norm,
                 dev_spin2, dev_spin2_norm,
                 dev_spin3, dev_spin3_norm,
                 dev_spin4, dev_spin4_norm,
                 dev_spin5, dev_spin5_norm,
                 dev_grid,dev_nn, (float) finaleps); 
   #endif
   stopinner = clock();
   timeelapsed = (double) (stopinner-startinner)/CLOCKS_PER_SEC;
   
   #ifndef LOWOUTPUT
   if (g_cart_id == 0) printf("Inner solver done\nTime elapsed: %.6e sec\n", timeelapsed);
   #endif

   #ifdef GPU_DOUBLE
     dev_add_f2d<<<gpu_gd_blas_d, gpu_bd_blas_d >>>(dev_spin1_d,dev_spin1_d,dev_spinout);
     dev_f2d<<<gpu_gd_blas_d, gpu_bd_blas_d >>>(dev_spin2_d,dev_spinout);
   #else
     #ifdef RELATIVISTIC_BASIS
        to_tmlqcd_basis<<<gpu_gd_linalg, gpu_bd_linalg>>> (dev_spinout);
     #endif     
     // copy back
     cudaMemcpy(h2d_spin, dev_spinout, dev_spinsize, cudaMemcpyDeviceToHost);
     #ifdef HALF
      cudaMemcpy(h2d_spin_norm, dev_spinout_norm, dev_normsize, cudaMemcpyDeviceToHost);
     #endif

     if ((cudaerr=cudaGetLastError())!=cudaSuccess) {
       if (g_cart_id == 0) printf("%s\n", cudaGetErrorString(cudaerr));
     }
     #ifndef HALF
       convert2double_spin(h2d_spin, solver_field[2]);
     #else
       convert2double_spin_half(h2d_spin, h2d_spin_norm, solver_field[2]);
     #endif
   
     // x_(k+1) = x_k + p_k
     add(solver_field[1],solver_field[1],solver_field[2],N);
   #endif
   outercount ++;   
}// outer loop 
    
     if (g_cart_id == 0) printf("Did NOT reach solver precision of eps=%.2e\n",eps);
   #ifdef GPU_DOUBLE
      #ifdef RELATIVISTIC_BASIS
        to_tmlqcd_basis_d<<<gpu_gd_linalg_d, gpu_bd_linalg_d>>> (dev_spin1_d);
      #endif        
      cudaMemcpy(h2d_spin_d, dev_spin1_d, dev_spinsize_d, cudaMemcpyDeviceToHost);
      unorder_spin_gpu(h2d_spin_d, P);   
   #else 
     assign(P, solver_field[1], N);
   #endif 

    assert((stop = clock())!=-1);
    timeelapsed = (double) (stop-start)/CLOCKS_PER_SEC;
    if (g_cart_id == 0) printf("Inversion done in mixed precision.\n Number of iterations in outer solver: %d\n Squared residue: %.8e\n Time elapsed: %.6e sec\n", outercount, rk, timeelapsed);

    finalize_solver(solver_field, nr_sf);
    finalize_mixedsolve_fields();
  return(-1);
}








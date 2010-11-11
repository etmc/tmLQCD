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
#include "cublas.h"
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "../global.h"
#include "cudaglobal.h"
//#include "mixed_solve.h"
#include "HEADER.h"
#include "cudadefs.h"
#include <math.h>


extern "C" {
#include "../tm_operators.h"
#include "../linalg_eo.h"
#include "../start.h"
#include "../complex.h"
#include "../read_input.h"
#include "../geometry_eo.h"
#include "../boundary.h"
#include "../su3.h"
#include "../temporalgauge.h"
#include "../observables.h"
#include "../measure_rectangles.h"
#include "../polyakov_loop.h"
}



#ifdef HAVE_CONFIG_H
  #include<config.h>
#endif



int g_numofgpu;

#ifdef GF_8
dev_su3_8 * dev_gf;
dev_su3_8 * h2d_gf;
#else
dev_su3_2v * dev_gf;
dev_su3_2v * h2d_gf;
#endif


dev_spinor* dev_spin1;
dev_spinor* dev_spin2;
dev_spinor* dev_spin3;
dev_spinor* dev_spin4;
dev_spinor* dev_spin5;
dev_spinor* dev_spinin;
dev_spinor* dev_spinout;
dev_spinor * h2d_spin;


#ifdef HALF
  // some additional fields for half prec.
  dev_spinor_half* dev_half_aux;
  float* dev_half_norm;
  // a half precsion gauge field
  #ifdef GF_8
   dev_su3_8 * dev_gf_half;
  #else
   dev_su3_2v * dev_gf_half;
  #endif
#endif 

//additional spinors for even-odd
dev_spinor* dev_spin_eo1;
dev_spinor* dev_spin_eo2;

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


size_t output_size;
int* dev_grid;
float * dev_output;
int havedevice = 0;


REAL hostr;
REAL hostkappa;
REAL hostm;
REAL hostmu;


__device__  REAL m;
__device__  REAL mu;
__device__  REAL r=1.0; // this is implicitly assumed to be 1.0 in the host code!!!
__device__  REAL kappa;
__device__ REAL twokappamu;

__device__ dev_complex dev_k0;
__device__ dev_complex dev_k1;
__device__ dev_complex dev_k2;
__device__ dev_complex dev_k3;

__device__ dev_complex dev_mk0;
__device__ dev_complex dev_mk1;
__device__ dev_complex dev_mk2;
__device__ dev_complex dev_mk3;



__constant__ __device__ dev_complex dev_k0c;
__constant__ __device__ dev_complex dev_k1c;
__constant__ __device__ dev_complex dev_k2c;
__constant__ __device__ dev_complex dev_k3c;

__constant__ __device__ dev_complex dev_mk0c;
__constant__ __device__ dev_complex dev_mk1c;
__constant__ __device__ dev_complex dev_mk2c;
__constant__ __device__ dev_complex dev_mk3c;



__device__  int  dev_LX,dev_LY,dev_LZ,dev_T,dev_VOLUME;


#ifdef MPI
__device__ int dev_rank;		// will be put to mixed_solve_eo_nd.cuh ...
__device__ int dev_nproc;
#endif


// include files with other GPU code as all GPU code has to reside in one file 
// the texture references and functions
#include "textures.cuh"
// linear algebra functions and gamma-multiplications
#include "linalg.cuh"
// reconstruction of the gauge field
#include "gauge_reconstruction.cuh"
// the device Hopping_Matrix
#include "Hopping_Matrix.cuh"
// the non-EO twisted mass dirac operator
#include "tm_diracoperator.cuh"
// the device su3 functions
#include "su3.cuh"
// the plaquette and rectangle routines
#include "observables.cuh"

// if we want to use half precision
#ifdef HALF 
 #include "half.cuh"
#endif



// computes sout = 1/(1 +- mutilde gamma5) sin = (1 -+ i mutilde gamma5)/(1+mutilde^2) sin
// mutilde = 2 kappa mu
// uses shared local memory for manipulation
__global__ void dev_mul_one_pm_imu_inv(dev_spinor* sin, dev_spinor* sout, const REAL sign){
   
   dev_spinor slocal[6];
   //need the inverse sign in the numerator because of inverse
   dev_complex pm_imu = dev_initcomplex(0.0,-1.0*sign*twokappamu);
   
   REAL one_plus_musquare_inv = 1.0/(1.0 + twokappamu*twokappamu);
   int pos;
   pos= threadIdx.x + blockDim.x*blockIdx.x;  
   int ix = threadIdx.x;
   if(pos < dev_VOLUME){
     //dev_skalarmult_spinor(&(sin[6*pos]), pm_imu, &(slocal[0]));
     //dev_Gamma5(&(slocal[0]));
     dev_skalarmult_gamma5_spinor(&(slocal[0]), pm_imu, &(sin[6*pos]) );
	 dev_add_spinor_assign(&(slocal[0]), &(sin[6*pos]));
     //dev_realmult_spinor(&(slocal[0]), one_plus_musquare_inv);
     //dev_copy_spinor(&(slocal[0]), &(sout[6*pos])); 
     dev_realmult_spinor_assign(&(sout[6*pos]), one_plus_musquare_inv, &(slocal[0]) );
   }
}





// sout = gamma_5*((1\pm i\mutilde \gamma_5)*sin1 - sin2)
// uses shared local memory for manipulation
__global__ void dev_mul_one_pm_imu_sub_mul_gamma5(dev_spinor* sin1, dev_spinor* sin2, dev_spinor* sout, const REAL sign){
   dev_spinor slocal[6];
   dev_complex pm_imu = dev_initcomplex(0.0, sign*twokappamu); // i mutilde
   int pos;
   pos= threadIdx.x + blockDim.x*blockIdx.x; 
   int ix = threadIdx.x;
   if(pos < dev_VOLUME){
     //dev_skalarmult_spinor(&(sin1[6*pos]), pm_imu, &(slocal[0]));
     //dev_Gamma5(&(slocal[0]));
     dev_skalarmult_gamma5_spinor(&(slocal[0]),pm_imu,&(sin1[6*pos]));
	 dev_add_spinor_assign(&(slocal[0]), &(sin1[6*pos]));
     dev_sub_spinor_assign(&(slocal[0]), &(sin2[6*pos]));
     //dev_Gamma5(&(slocal[0]));
     //dev_copy_spinor(&(slocal[0]), &(sout[6*pos]));
     dev_Gamma5_assign(&(sout[6*pos]), &(slocal[0]));
   }   
}







// aequivalent to Qtm_pm_psi in tm_operators.c
extern "C" void dev_Qtm_pm_psi(dev_spinor* spinin, dev_spinor* spinout, int gridsize, int blocksize, int gridsize2, int blocksize2){
  //spinin == odd
  //spinout == odd
  
  //Q_{-}
  #ifdef USETEXTURE
   #ifndef HALF
    bind_texture_spin(spinin,1);
   #else
    prepare_halfspinor_texture(spinin);
   #endif
  #endif
  //bind_texture_nn(dev_nn_eo);
  //cudaFuncSetCacheConfig(dev_Hopping_Matrix, cudaFuncCachePreferL1);
    dev_Hopping_Matrix<<<gridsize, blocksize>>>
             (dev_gf, spinin, dev_spin_eo1, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0); //dev_spin_eo1 == even -> 0           
  //unbind_texture_nn();           
  #ifdef USETEXTURE
   #ifndef HALF
    unbind_texture_spin(1);
   #else
    release_halfspinor_texture();
   #endif
  #endif
  dev_mul_one_pm_imu_inv<<<gridsize2, blocksize2>>>(dev_spin_eo1,dev_spin_eo2, -1.);
  
  #ifdef USETEXTURE
   #ifndef HALF
    bind_texture_spin(dev_spin_eo2,1);
   #else
    prepare_halfspinor_texture(dev_spin_eo2);
   #endif
  #endif
  //bind_texture_nn(dev_nn_oe);
  //cudaFuncSetCacheConfig(dev_Hopping_Matrix, cudaFuncCachePreferL1);
    dev_Hopping_Matrix<<<gridsize, blocksize>>>
            (dev_gf, dev_spin_eo2, dev_spin_eo1, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1); 
  //unbind_texture_nn();
  #ifdef USETEXTURE
   #ifndef HALF
    unbind_texture_spin(1);
   #else
    release_halfspinor_texture();
   #endif
  #endif
  dev_mul_one_pm_imu_sub_mul_gamma5<<<gridsize2, blocksize2>>>(spinin, dev_spin_eo1,  dev_spin_eo2, -1.);
  
  //Q_{+}
  #ifdef USETEXTURE
   #ifndef HALF
    bind_texture_spin(dev_spin_eo2,1);
   #else
    prepare_halfspinor_texture(dev_spin_eo2);
   #endif
  #endif
  //bind_texture_nn(dev_nn_eo);
  //cudaFuncSetCacheConfig(dev_Hopping_Matrix, cudaFuncCachePreferL1);
    dev_Hopping_Matrix<<<gridsize, blocksize>>>
          (dev_gf, dev_spin_eo2, dev_spin_eo1, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0); //dev_spin_eo1 == even -> 0
  //unbind_texture_nn();      
  #ifdef USETEXTURE  
   #ifndef HALF
    unbind_texture_spin(1);
   #else
    release_halfspinor_texture();
   #endif
  #endif
  dev_mul_one_pm_imu_inv<<<gridsize2, blocksize2>>>(dev_spin_eo1,spinout, +1.);
  
  #ifdef USETEXTURE
   #ifndef HALF
    bind_texture_spin(spinout,1);
   #else
    prepare_halfspinor_texture(spinout);
   #endif
  #endif
  //bind_texture_nn(dev_nn_oe);
  //cudaFuncSetCacheConfig(dev_Hopping_Matrix, cudaFuncCachePreferL1);
    dev_Hopping_Matrix<<<gridsize, blocksize>>>
             (dev_gf, spinout, dev_spin_eo1, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1); 
  //unbind_texture_nn();  
  #ifdef USETEXTURE
   #ifndef HALF
    unbind_texture_spin(1);
   #else
    release_halfspinor_texture();
   #endif
  #endif
  dev_mul_one_pm_imu_sub_mul_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo2, dev_spin_eo1,  spinout , +1.); 
}








__global__ void dev_zero_spinor_field(dev_spinor* s1){
  int pos;
  pos= threadIdx.x + blockDim.x*blockIdx.x;  
  if(pos < dev_VOLUME){
          dev_zero_spinor(&(s1[6*pos]));
  }
}




__global__ void dev_copy_spinor_field(dev_spinor* s1, dev_spinor* s2){
    int pos;
  pos= threadIdx.x + blockDim.x*blockIdx.x;  
  if(pos < dev_VOLUME){
      dev_copy_spinor(&(s1[6*pos]),&(s2[6*pos]));
  } 
}



__global__ void dev_skalarmult_add_assign_spinor_field(dev_spinor* s1, REAL lambda, dev_spinor* s2, dev_spinor* so){
  int pos;
  pos= threadIdx.x + blockDim.x*blockIdx.x;  
  if(pos < dev_VOLUME){
    dev_skalarmult_add_assign_spinor(&(s1[6*pos]), lambda ,&(s2[6*pos]), &(so[6*pos]) );
  }
}



__global__ void dev_skalarmult_spinor_field(dev_spinor* s1, REAL lambda, dev_spinor* so){
  int pos;
  pos= threadIdx.x + blockDim.x*blockIdx.x;  
  if(pos < dev_VOLUME){
    dev_skalarmult_spinor(&(s1[6*pos]), dev_initcomplex(lambda,0.0) , &(so[6*pos]) );
  }
}  



__global__ void dev_complexmult_spinor_field(dev_spinor* s1, dev_complex lambda, dev_spinor* so){
  int pos;
  pos= threadIdx.x + blockDim.x*blockIdx.x;  
  if(pos < dev_VOLUME){
    dev_skalarmult_spinor(&(s1[6*pos]), lambda , &(so[6*pos]) );
  }
}






// init the gpu inner solver, assigen constants etc.
__global__ void he_cg_init (int* grid, REAL param_kappa, REAL param_mu, dev_complex k0, dev_complex k1, dev_complex k2, dev_complex k3){
  dev_LX = grid[0];
  dev_LY = grid[1];
  dev_LZ = grid[2];
  dev_T = grid[3];
  dev_VOLUME = grid[4]; // grid[4] is initialized 1/2 VOLUME for eo
  
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
}





// code to list available devices, not yet included in main code
// this is copied from the CUDA sdk 
extern "C" int find_devices(){
int deviceCount, dev;
    cudaGetDeviceCount(&deviceCount);
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
    return(deviceCount);
}









extern "C" void test_operator(dev_su3_2v * gf,dev_spinor* spinin, dev_spinor* spinout, 
dev_spinor* spin0, dev_spinor* spin1, dev_spinor* spin2, dev_spinor* spin3, dev_spinor* spin4, int *grid, int * nn_grid, REAL* output,REAL* erg, int xsize, int ysize){
 
 int  gridsize;

 dim3 blockdim(1,1);
 dim3 blockdim2(128,1,1);
 if( VOLUME >= 128){
   gridsize =VOLUME/128;
 }
 else{
   gridsize=1;
 }
 dim3 griddim2(gridsize,1,1);
 
 
 dim3 blockdim3(BLOCK,1,1);
 if( VOLUME >= BLOCK){
   gridsize = (int) VOLUME/BLOCK + 1;
 }
 else{
   gridsize=1;
 }
 dim3 griddim3(gridsize,1,1); 
 
 
  dev_complex h0,h1,h2,h3;
  h0.re = (REAL)ka0.re;    h0.im = (REAL)ka0.im;
  h1.re = (REAL)ka1.re;    h1.im = (REAL)ka1.im;
  h2.re = (REAL)ka2.re;    h2.im = (REAL)ka2.im;
  h3.re = (REAL)ka3.re;    h3.im = (REAL)ka3.im;
  he_cg_init<<< 1, 1 >>> (grid, (REAL) g_kappa, (REAL)(g_mu/(2.0*g_kappa)), h0,h1,h2,h3);
 
 
  REAL scaleparam = sqrt(1.0/(2.0 * (REAL) hostkappa));
  dev_skalarmult_spinor_field<<<griddim2, blockdim2 >>>(spinin,scaleparam*scaleparam, spin4);
 
 #ifdef USETEXTURE
   bind_texture_gf(gf);
   bind_texture_spin(spin4,1);
 #endif 
  // apply D_tm
  dev_tm_dirac_kappa <<<griddim3, blockdim3 >>>(gf, spin4, spinout, nn_grid);

 #ifdef USETEXTURE
  unbind_texture_gf();
  unbind_texture_spin(1);
 #endif
}





// this is the eo version of the device cg inner solver 
// we invert the hermitean D_tm D_tm^{+}
extern "C" int dev_cg(dev_su3_2v * gf,dev_spinor* spinin, dev_spinor* spinout, 
dev_spinor* spin0, dev_spinor* spin1, dev_spinor* spin2, dev_spinor* spin3, dev_spinor* spin4, int *grid, int * nn_grid, REAL* output,REAL* erg, int xsize, int ysize, int rescalekappa){
 
 
 REAL host_alpha, host_beta, host_dotprod, host_rk, sourcesquarenorm;
 REAL * dotprod, * dotprod2, * rk, * alpha, *beta;
 
 
 
 int i, gridsize;
 int maxit = max_innersolver_it;
 REAL eps = (REAL) innersolver_precision;
 int N_recalcres = 10; // after N_recalcres iterations calculate r = A x_k - b
 
 
 // initialize grid and block, make sure VOLUME is a multiple of blocksize 
 if(VOLUME%DOTPROD_DIM != 0){
   printf("Error: VOLUME is not a multiple of DOTPROD_DIM. Aborting...\n");
   exit(100); 
 }

 // this is the partitioning for the copying of fields 
 dim3 blockdim(1,1);
 dim3 blockdim2(128,1,1);
 if( VOLUME >= 128){
   gridsize = (int) VOLUME/128 + 1;
 }
 else{
   gridsize=1;
 }
 dim3 griddim2(gridsize,1,1);
 
 // this is the partitioning for the Dirac-Kernel
 dim3 blockdim3(BLOCK,1,1);
 if( VOLUME >= BLOCK){
   gridsize = (int) (VOLUME/BLOCK) +1;
 }
 else{
   gridsize=1;
 }
 dim3 griddim3(gridsize,1,1); 
 
 size_t size2 = sizeof(float4)*6*VOLUME;
 
 #ifdef USETEXTURE
   //Bind texture gf
   bind_texture_gf(gf);
  //Bind texture spinor to spin4 (D_tm is always applied to spin4)
  bind_texture_spin(spin4,1);
 #endif
 
 //Initialize some stuff
  printf("mu = %f\n", g_mu);
  dev_complex h0,h1,h2,h3;
  h0.re = (REAL)ka0.re;    h0.im = (REAL)ka0.im;
  h1.re = (REAL)ka1.re;    h1.im = (REAL)ka1.im;
  h2.re = (REAL)ka2.re;    h2.im = (REAL)ka2.im;
  h3.re = (REAL)ka3.re;    h3.im = (REAL)ka3.im;
  he_cg_init<<< 1, 1 >>> (grid, (REAL) g_kappa, (REAL)(g_mu/(2.0*g_kappa)), h0,h1,h2,h3);
  // BEWARE in dev_tm_dirac_kappa we need the true mu (not 2 kappa mu!)
 
 // Init x,p,r for k=0
 // Allocate some numbers for host <-> device interaction
 cudaMalloc((void **) &dotprod, sizeof(REAL));
 cudaMalloc((void **) &dotprod2, sizeof(REAL));
 cudaMalloc((void **) &rk, sizeof(REAL));
 cudaMalloc((void **) &alpha, sizeof(REAL));
 cudaMalloc((void **) &beta, sizeof(REAL));
 printf("%s\n", cudaGetErrorString(cudaGetLastError())); 
 
 
 //init blas
 cublasInit();
 printf("%s\n", cudaGetErrorString(cudaGetLastError())); 
 printf("have initialized cublas\n");
 
 
 // go over to kappa (if wanted)
 REAL scaleparam = sqrt(1.0/(2.0 * (REAL)hostkappa));
 printf("1/2kappa = %.8f\n",scaleparam);
 //dev_skalarmult_spinor_field<<<griddim2, blockdim2 >>>(spinin,scaleparam, spin1);
 //dev_copy_spinor_field<<<griddim2, blockdim2 >>>(spin1, spinin);
 
 
 dev_copy_spinor_field<<<griddim2, blockdim2 >>>(spinin, spin0);
 dev_zero_spinor_field<<<griddim2, blockdim2 >>>(spin1); // x_0 = 0
 dev_copy_spinor_field<<<griddim2, blockdim2 >>>(spinin, spin2);
 dev_zero_spinor_field<<<griddim2, blockdim2 >>>(spin3);
 printf("%s\n", cudaGetErrorString(cudaGetLastError()));
 
 
 
 
 //relative precision -> get initial residue
 sourcesquarenorm = cublasSdot (24*VOLUME, (const float *)spinin, 1, (const float *)spinin, 1);
 host_rk = sourcesquarenorm; //for use in main loop
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
     dev_gamma5 <<<griddim2, blockdim2 >>> (spin2,spin4);
     dev_swapmu <<<1,1>>> ();
  #ifdef USETEXTURE
   bind_texture_spin(spin4,1);
  #endif
     //D_tm 
     dev_tm_dirac_kappa <<<griddim3, blockdim3 >>> (gf, spin4, spin3, dev_nn);
  #ifdef USETEXTURE
   unbind_texture_spin(1);
  #endif
     //GAMMA5 mu -> -mu
     dev_gamma5 <<<griddim2, blockdim2 >>>(spin3,spin4);
     dev_swapmu <<<1,1>>> ();
  #ifdef USETEXTURE
   bind_texture_spin(spin4,1);
  #endif
     //D_tm
     dev_tm_dirac_kappa <<<griddim3, blockdim3 >>> (gf, spin4, spin3, dev_nn);
  
  //Here we have used the output spinor (spinout) to temporarly take the field and to 
  //copy it to the texture field (spin4)!!

  
 //alpha
  host_dotprod = cublasSdot (24*VOLUME, (const float *) spin2, 1,
            (const float *) spin3, 1);
  host_alpha = (host_rk / host_dotprod); // alpha = r*r/ p M p
   
 //r(k+1)
 cublasSaxpy (24*VOLUME,-1.0*host_alpha, (const float *) spin3, 1, (float *) spin0, 1);  

 //x(k+1);
 cublasSaxpy (24*VOLUME, host_alpha, (const float *) spin2,  1, (float *) spin1, 1);

 printf("%s\n", cudaGetErrorString(cudaGetLastError()));

  //Abbruch?
  host_dotprod = cublasSdot (24*VOLUME, (const float *) spin0, 1,(const float *) spin0, 1);
  
 if ((host_dotprod <= eps*sourcesquarenorm)){//error-limit erreicht
   break; 
 }
  printf("iter %d: err = %.8e\n", i, host_dotprod);
  
 //beta
 host_beta =host_dotprod/host_rk;
 //p(k+1)
 cublasSscal (24*VOLUME, host_beta, (float *)spin2, 1);
 cublasSaxpy (24*VOLUME, 1.0, (const float *) spin0,  1, (float *) spin2, 1);

 host_rk = host_dotprod;
 
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
      dev_gamma5 <<<griddim2, blockdim2 >>> (spin1,spin4);
      dev_swapmu <<<1,1>>> ();
    #ifdef USETEXTURE
     bind_texture_spin(spin4,1);
    #endif
   
      //D_tm GAMMA5, mu -> -mu
      dev_tm_dirac_kappa <<<griddim3, blockdim3 >>> (gf, spin4, spin3, dev_nn);
      dev_gamma5 <<<griddim2, blockdim2 >>>(spin3,spinout);
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
      dev_tm_dirac_kappa<<<griddim3, blockdim3 >>>(gf, spin4, spin3, dev_nn);
    
    // r = b - Ax
    cublasSscal (24*VOLUME, -1.0, (float *)spin3, 1);
    cublasSaxpy (24*VOLUME, 1.0, (const float *) spinin,  1, (float *) spin3, 1);
    cublasScopy (24*VOLUME, (const float *)spin3, 1, (float *)spin0, 1);
    
    //dev_skalarmult_add_assign_spinor_field<<<griddim2, blockdim2 >>>(spinin, -1.0, spin3, spin0);
   }//recalculate residue

 }//MAIN LOOP cg	
  
  
  printf("Final residue: %.6e\n",host_dotprod);
  // x_result = spin1 !
  
 if(rescalekappa == 1){  //want D^-1 rescaled by 2*kappa
  
//multiply with D^dagger
    #ifdef USETEXTURE
     unbind_texture_spin(1);
    #endif
      dev_gamma5 <<<griddim2, blockdim2 >>> (spin1,spin4);
      dev_swapmu <<<1,1>>> ();
    #ifdef USETEXTURE
     bind_texture_spin(spin4,1);
    #endif
      dev_tm_dirac_kappa <<<griddim3, blockdim3 >>> (gf, spin4, spin3, dev_nn);
      dev_gamma5 <<<griddim2, blockdim2 >>>(spin3,spin1);
      dev_swapmu <<<1,1>>> ();
    #ifdef USETEXTURE
     unbind_texture_spin(1);
    #endif


 //go over to non-kappa, Ddagger = g5 D g5
 dev_skalarmult_spinor_field<<<griddim2, blockdim2 >>>(spin1,1.0/(scaleparam*scaleparam), spinout);  
 
  // times operator == source ?? 
  //dev_tm_dirac_kappa<<<griddim3, blockdim3 >>>(gf, spin3, spinout, nn_grid);
  }
  else{
   dev_copy_spinor_field<<<griddim2, blockdim2 >>>(spin1,spinout);
  }
  
  #ifdef USETEXTURE
   unbind_texture_gf();
  #endif
  cudaFree(dotprod);
  cudaFree(dotprod2);
  cudaFree(rk);
  cudaFree(alpha);
  cudaFree(beta);
  cublasShutdown();
  return(i);
}





// this is the eo version of the device cg inner solver 
// we invert the hermitean Q_{-} Q_{+}
extern "C" int dev_cg_eo(dev_su3_2v * gf,dev_spinor* spinin, dev_spinor* spinout, 
dev_spinor* spin0, dev_spinor* spin1, dev_spinor* spin2, dev_spinor* spin3, dev_spinor* spin4, int *grid, int * nn_grid, REAL* output,REAL* erg, int xsize, int ysize, int rescalekappa, REAL epsfinal){
 
 
 REAL host_alpha, host_beta, host_dotprod, host_rk, sourcesquarenorm;
 REAL * dotprod, * dotprod2, * rk, * alpha, *beta;
 
 
 
 int i, gridsize;
 int maxit = max_innersolver_it;
 REAL eps = (REAL) innersolver_precision;
 int N_recalcres = 20; // after N_recalcres iterations calculate r = A x_k - b
 
 cudaError_t cudaerr;
 // this is the partitioning for the copying of fields
 dim3 blockdim(1,1);
 dim3 blockdim2(128,1,1);
 if( VOLUME/2 >= 128){
   gridsize = (int) VOLUME/2/128 + 1;
 }
 else{
   gridsize=1;
 }
 dim3 griddim2(gridsize,1,1);

 
 //this is the partitioning for the HoppingMatrix kernel
 int blockdim3=BLOCK;
 if( VOLUME/2 >= BLOCK){
   gridsize = (int)(VOLUME/2/BLOCK) + 1;
 }
 else{
   gridsize=1;
 }
 printf("gridsize = %d\n", gridsize);
 int griddim3=gridsize; 
 
 //this is the partitioning for dev_mul_one_pm...
 int blockdim4=BLOCK2;
 if( VOLUME/2 >= BLOCK2){
   gridsize = (int)(VOLUME/2/BLOCK2) + 1;
 }
 else{
   gridsize=1;
 }
 int griddim4=gridsize;  
 
 
 size_t size2 = sizeof(float4)*6*VOLUME/2;
 
 #ifdef USETEXTURE
  //Bind texture gf
  bind_texture_gf(gf);
 #endif
 //Initialize some stuff
  printf("mu = %f\n", g_mu);
  dev_complex h0,h1,h2,h3,mh0, mh1, mh2, mh3;
  h0.re = (REAL)ka0.re;    h0.im = -(REAL)ka0.im;
  h1.re = (REAL)ka1.re;    h1.im = -(REAL)ka1.im;
  h2.re = (REAL)ka2.re;    h2.im = -(REAL)ka2.im;
  h3.re = (REAL)ka3.re;    h3.im = -(REAL)ka3.im;
  
  mh0.re = -(REAL)ka0.re;    mh0.im = (REAL)ka0.im;
  mh1.re = -(REAL)ka1.re;    mh1.im = (REAL)ka1.im;
  mh2.re = -(REAL)ka2.re;    mh2.im = (REAL)ka2.im;
  mh3.re = -(REAL)ka3.re;    mh3.im = (REAL)ka3.im;
  
  // try using constant mem for kappas
  cudaMemcpyToSymbol("dev_k0c", &h0, sizeof(dev_complex)) ; 
  cudaMemcpyToSymbol("dev_k1c", &h1, sizeof(dev_complex)) ; 
  cudaMemcpyToSymbol("dev_k2c", &h2, sizeof(dev_complex)) ; 
  cudaMemcpyToSymbol("dev_k3c", &h3, sizeof(dev_complex)) ;
  
  cudaMemcpyToSymbol("dev_mk0c", &mh0, sizeof(dev_complex)) ; 
  cudaMemcpyToSymbol("dev_mk1c", &mh1, sizeof(dev_complex)) ; 
  cudaMemcpyToSymbol("dev_mk2c", &mh2, sizeof(dev_complex)) ; 
  cudaMemcpyToSymbol("dev_mk3c", &mh3, sizeof(dev_complex)) ;  
  
  he_cg_init<<< 1, 1 >>> (grid, (REAL) g_kappa, (REAL)(g_mu/(2.0*g_kappa)), h0,h1,h2,h3);
  // BEWARE in dev_tm_dirac_kappa we need the true mu (not 2 kappa mu!)
 
 
 
 // Init x,p,r for k=0
 // Allocate some numbers for host <-> device interaction
 cudaMalloc((void **) &dotprod, sizeof(REAL));
 cudaMalloc((void **) &dotprod2, sizeof(REAL));
 cudaMalloc((void **) &rk, sizeof(REAL));
 cudaMalloc((void **) &alpha, sizeof(REAL));
 cudaMalloc((void **) &beta, sizeof(REAL));
 printf("%s\n", cudaGetErrorString(cudaGetLastError())); 
 
 
 //init blas
 cublasInit();
 printf("%s\n", cudaGetErrorString(cudaGetLastError())); 
 printf("have initialized cublas\n");
 
 
 // go over to kappa (if wanted)
 REAL scaleparam = sqrt(1.0/(2.0 * (REAL)hostkappa));
 printf("1/2kappa = %.8f\n",scaleparam);
 //dev_skalarmult_spinor_field<<<griddim2, blockdim2 >>>(spinin,scaleparam, spin1);
 //dev_copy_spinor_field<<<griddim2, blockdim2 >>>(spin1, spinin);
 
 
 dev_copy_spinor_field<<<griddim2, blockdim2 >>>(spinin, spin0);
 dev_zero_spinor_field<<<griddim2, blockdim2 >>>(spin1); // x_0 = 0
 dev_copy_spinor_field<<<griddim2, blockdim2 >>>(spinin, spin2);
 dev_zero_spinor_field<<<griddim2, blockdim2 >>>(spin3);
 printf("%s\n", cudaGetErrorString(cudaGetLastError()));
 
 
 

 //relative precision -> get initial residue
 sourcesquarenorm = cublasSdot (24*VOLUME/2, (const float *)spinin, 1, (const float *)spinin, 1);
 host_rk = sourcesquarenorm; //for use in main loop
 printf("Squarenorm Source:\t%.8e\n", sourcesquarenorm);
 printf("%s\n", cudaGetErrorString(cudaGetLastError()));
 
  printf("Entering cg-loop\n");
 for(i=0;i<maxit;i++){ //MAIN LOOP
  
  // Q_{-}Q{+}
  dev_Qtm_pm_psi(spin2, spin3, griddim3, blockdim3, griddim4, blockdim4);
  if((cudaerr=cudaGetLastError()) != cudaSuccess){
    printf("%s\n", cudaGetErrorString(cudaerr));
    exit(200);
  }
  
  
  
 //alpha
  host_dotprod = cublasSdot (24*VOLUME/2, (const float *) spin2, 1,
            (const float *) spin3, 1);
  host_alpha = (host_rk / host_dotprod); // alpha = r*r/ p M p
   
 //r(k+1)
 cublasSaxpy (24*VOLUME/2,-1.0*host_alpha, (const float *) spin3, 1, (float *) spin0, 1);  

 //x(k+1);
 cublasSaxpy (24*VOLUME/2, host_alpha, (const float *) spin2,  1, (float *) spin1, 1);

 printf("%s\n", cudaGetErrorString(cudaGetLastError()));

  //Abbruch?
  host_dotprod = cublasSdot (24*VOLUME/2, (const float *) spin0, 1,(const float *) spin0, 1);
  
 if (((host_dotprod <= eps*sourcesquarenorm) && (i > maxit / 4) ) || ( host_dotprod <= epsfinal/2.)){//error-limit erreicht (epsfinal/2 sollte ausreichen um auch in double precision zu bestehen)
   break; 
 }
  printf("iter %d: err = %.8e\n", i, host_dotprod);
  
 //beta
 host_beta =host_dotprod/host_rk;
 //p(k+1)
 cublasSscal (24*VOLUME/2, host_beta, (float *)spin2, 1);
 cublasSaxpy (24*VOLUME/2, 1.0, (const float *) spin0,  1, (float *) spin2, 1);

 host_rk = host_dotprod;
 
 // recalculate residue frome r = b - Ax
 if(((i+1) % N_recalcres) == 0){
    // r_(k+1) = Ax -b 
    printf("Recalculating residue\n");
    
    // D Ddagger   --   Ddagger = gamma5 D gamma5  for Wilson Dirac Operator
    // DO NOT USE tm_dirac_dagger_kappa here, otherwise spin2 will be overwritten!!!
      
    // Q_{-}Q{+}
    dev_Qtm_pm_psi(spin1, spin3, griddim3, blockdim3, griddim4, blockdim4);
      
        
    
    // r = b - Ax
    cublasSscal (24*VOLUME/2, -1.0, (float *)spin3, 1);
    cublasSaxpy (24*VOLUME/2, 1.0, (const float *) spinin,  1, (float *) spin3, 1);
    cublasScopy (24*VOLUME/2, (const float *)spin3, 1, (float *)spin0, 1);
    //dev_skalarmult_add_assign_spinor_field<<<griddim2, blockdim2 >>>(spinin, -1.0, spin3, spin0);
   }//recalculate residue

 }//MAIN LOOP cg	
  
  
  printf("Final residue: %.6e\n",host_dotprod);
  // x_result = spin1 !
  
  //no multiplication with D^{dagger} here and no return to non-kappa basis as in dev_cg!
  dev_copy_spinor_field<<<griddim2, blockdim2 >>>(spin1,spinout);
  
  #ifdef USETEXTURE
   unbind_texture_gf();
  #endif
  cudaFree(dotprod);
  cudaFree(dotprod2);
  cudaFree(rk);
  cudaFree(alpha);
  cudaFree(beta);
  cublasShutdown();
  return(i);
}












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









// convert spinor to double 
void convert2double_spin(dev_spinor* spin, spinor* h2d){
  int i,Vol;
  if(even_odd_flag){
    Vol = VOLUME/2;
  }
  else{
    Vol = VOLUME;
  }
  for (i=0;i<Vol;i++){
  
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
  int i,Vol;
  if(even_odd_flag){
    Vol = VOLUME/2;
  }
  else{
    Vol = VOLUME;
  }
  for (i=0;i<Vol;i++){
    
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





void init_mixedsolve(su3** gf){
cudaError_t cudaerr;

   // get number of devices
   if(havedevice == 0){
     int ndev = find_devices();
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
	    printf("Error in init_mixedsolve_eo(): Could not set active device. Aborting...\n");
	    exit(302);
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
  h2d_gf = (dev_su3_8 *)malloc(dev_gfsize); // Allocate REAL conversion gf on host
  su3to8(gf,h2d_gf);  
  #else
  h2d_gf = (dev_su3_2v *)malloc(dev_gfsize); // Allocate REAL conversion gf on host
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
  size_t dev_spinsize = 6*VOLUME * sizeof(dev_spinor); /* float4 */

  if((void*)(h2d_spin = (dev_spinor *)malloc(dev_spinsize)) == NULL){
    printf("Could not allocate memory for h2d_spin. Aborting...\n");
    exit(200);
  } // Allocate float conversion spinor on host
  
  cudaMalloc((void **) &dev_spin1, dev_spinsize);   // Allocate array spin1 on device
  cudaMalloc((void **) &dev_spin2, dev_spinsize);   // Allocate array spin2 on device
  cudaMalloc((void **) &dev_spin3, dev_spinsize);   // Allocate array spin3 on device
  cudaMalloc((void **) &dev_spin4, dev_spinsize);
  cudaMalloc((void **) &dev_spin5, dev_spinsize);
  cudaMalloc((void **) &dev_spinin, dev_spinsize);
  cudaMalloc((void **) &dev_spinout, dev_spinsize);
  if((cudaerr=cudaGetLastError())!=cudaSuccess){
    printf("Error in init_mixedsolve(): Memory allocation of spinor fields failed. Aborting...\n");
    exit(200);
  }
  else{
    printf("Allocated spinor fields on device\n");
  }
  
  
  output_size = LZ*T*sizeof(float); // parallel in t and z direction
  cudaMalloc((void **) &dev_output, output_size);   // output array
  float * host_output = (float*) malloc(output_size);

  int grid[5];
  grid[0]=LX; grid[1]=LY; grid[2]=LZ; grid[3]=T; grid[4]=VOLUME;
 
  cudaMalloc((void **) &dev_grid, 5*sizeof(int));
  cudaMemcpy(dev_grid, &(grid[0]), 5*sizeof(int), cudaMemcpyHostToDevice);
  
}






void init_mixedsolve_eo(su3** gf){
cudaError_t cudaerr;
  dev_complex help;

  if(havedevice == 0){
   // get number of devices
     int ndev = find_devices();
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
	    printf("Error in init_mixedsolve_eo(): Could not set active device. Aborting...\n");
	    exit(302);
	    }  
   havedevice=1;
  }
  #ifdef GF_8
  /* allocate 8 floats for gf = 2*4*VOLUME float4's*/
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
  h2d_gf = (dev_su3_8 *)malloc(dev_gfsize); // Allocate REAL conversion gf on host
  su3to8(gf,h2d_gf);
  #else
  h2d_gf = (dev_su3_2v *)malloc(dev_gfsize); // Allocate REAL conversion gf on host
  su3to2vf4(gf,h2d_gf);
  #endif
  //bring to device
  cudaMemcpy(dev_gf, h2d_gf, dev_gfsize, cudaMemcpyHostToDevice);

  #ifdef HALF
    #ifdef GF_8
      /* allocate 8 floats for gf = 2*4*VOLUME float4's*/
      printf("Using half precision GF 8 reconstruction\n");
      dev_gfsize = 2*4*VOLUME * sizeof(dev_su3_8_half); 
    #else
      /* allocate 2 rows of gf = 3*4*VOLUME float4's*/
      printf("Using half precision GF 12 reconstruction\n");
      dev_gfsize = 3*4*VOLUME * sizeof(dev_su3_2v_half); 
    #endif  
    
    if((cudaerr=cudaMalloc((void **) &dev_gf_half, dev_gfsize)) != cudaSuccess){
    printf("Error in init_mixedsolve(): Memory allocation of half precsion gauge field failed. Aborting...\n");
    exit(200);
    }   // Allocate array on device
    else{
      printf("Allocated half precision gauge field on device\n");
    }      
    
    int gridsize;
    // determine gridsize for conversion to half
    if( VOLUME >= BLOCK2){
       gridsize = (int) (VOLUME/BLOCK2) +1;
    }
    else{
      gridsize=1;
    }
   
   printf("Converting gauge to half precision... ");
     float2half_gaugefield <<< gridsize, BLOCK2  >>>(dev_gf, dev_gf_half);
   printf("Done\n");
  #endif


//grid 
  size_t nnsize = 8*VOLUME*sizeof(int);
  nn = (int *) malloc(nnsize);
  
  //nn grid for even-odd
  nn_eo = (int *) malloc(nnsize/2);
  nn_oe = (int *) malloc(nnsize/2);
  
  cudaMalloc((void **) &dev_nn, nnsize);
  cudaMalloc((void **) &dev_nn_eo, nnsize/2);
  cudaMalloc((void **) &dev_nn_oe, nnsize/2);
  
  
  size_t idxsize = VOLUME/2*sizeof(int);
  eoidx_even = (int *) malloc(idxsize);
  eoidx_odd = (int *) malloc(idxsize);
  cudaMalloc((void **) &dev_eoidx_even, idxsize);
  cudaMalloc((void **) &dev_eoidx_odd, idxsize);
  
  initnn();
  initnn_eo();
  //shownn_eo();
  
  //shownn();
  //showcompare_gf(T-1, LX-1, LY-1, LZ-1, 3);
  //check_gauge_reconstruction_8(gf, dev_gf, 0, 0);
  cudaMemcpy(dev_nn, nn, nnsize, cudaMemcpyHostToDevice);
  cudaMemcpy(dev_nn_eo, nn_eo, nnsize/2, cudaMemcpyHostToDevice);
  cudaMemcpy(dev_nn_oe, nn_oe, nnsize/2, cudaMemcpyHostToDevice);
  cudaMemcpy(dev_eoidx_even, eoidx_even, idxsize, cudaMemcpyHostToDevice);
  cudaMemcpy(dev_eoidx_odd, eoidx_odd, idxsize, cudaMemcpyHostToDevice);
  
  //free again
  free(eoidx_odd);
  free(eoidx_even);
  free(nn_oe);
  free(nn_eo);
  free(nn);
  
// Spinors
  size_t dev_spinsize = 6*VOLUME/2 * sizeof(dev_spinor); /* float4 */

  if((void*)(h2d_spin = (dev_spinor *)malloc(dev_spinsize)) == NULL){
    printf("Could not allocate memory for h2d_spin. Aborting...\n");
    exit(200);
  } // Allocate float conversion spinor on host
  
  cudaMalloc((void **) &dev_spin1, dev_spinsize);   // Allocate array spin1 on device
  cudaMalloc((void **) &dev_spin2, dev_spinsize);   // Allocate array spin2 on device
  cudaMalloc((void **) &dev_spin3, dev_spinsize);   // Allocate array spin3 on device
  cudaMalloc((void **) &dev_spin4, dev_spinsize);
  cudaMalloc((void **) &dev_spin5, dev_spinsize);
  cudaMalloc((void **) &dev_spinin, dev_spinsize);
  cudaMalloc((void **) &dev_spinout, dev_spinsize);
  
  cudaMalloc((void **) &dev_spin_eo1, dev_spinsize);
  cudaMalloc((void **) &dev_spin_eo2, dev_spinsize);
  
  if((cudaerr=cudaGetLastError())!=cudaSuccess){
    printf("Error in init_mixedsolve(): Memory allocation of spinor fields failed. Aborting...\n");
    exit(200);
  }
  else{
    printf("Allocated spinor fields on device\n");
  }
  
  #ifdef HALF
    //allocate half spinor
    dev_spinsize = 6*VOLUME/2 * sizeof(dev_spinor_half); /* short4 */
    cudaMalloc((void **) &dev_half_aux, dev_spinsize); 
    //allocate its norm
    dev_spinsize = VOLUME/2 * sizeof(float); /* float */
    cudaMalloc((void **) &dev_half_norm, dev_spinsize); 
  #endif
  if((cudaerr=cudaGetLastError())!=cudaSuccess){
    printf("Error in init_mixedsolve(): Memory allocation of half precision spinor fields failed. Aborting...\n");
    exit(200);
  }
  else{
    printf("Allocated half precision spinor fields on device\n");
  }
  
  
  output_size = LZ*T*sizeof(float); // parallel in t and z direction
  cudaMalloc((void **) &dev_output, output_size);   // output array
  float * host_output = (float*) malloc(output_size);

  int grid[5];
  grid[0]=LX; grid[1]=LY; grid[2]=LZ; grid[3]=T; grid[4]=VOLUME/2; 
  // dev_VOLUME is half of VOLUME for eo
 
  cudaMalloc((void **) &dev_grid, 5*sizeof(int));
  cudaMemcpy(dev_grid, &(grid[0]), 5*sizeof(int), cudaMemcpyHostToDevice);
  
  
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
    g_update_gauge_energy = 1;
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
    g_update_rectangle_energy = 1;
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

}



void finalize_mixedsolve(){

  cudaFree(dev_spin1);
  cudaFree(dev_spin2);
  cudaFree(dev_spin3);
  cudaFree(dev_spin4);
  cudaFree(dev_spin5);
  cudaFree(dev_spinin);
  cudaFree(dev_spinout);
  cudaFree(dev_gf);
  cudaFree(dev_grid);
  cudaFree(dev_output);
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
    cudaFree(dev_half_aux);
    cudaFree(dev_half_norm);
  #endif
  
  free(h2d_spin);
  free(h2d_gf);
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
  
  size_t dev_spinsize = 6*VOLUME * sizeof(dev_spinor); // float4 
  init_mixedsolve(g_gauge_field);
  
  // Start timer
  assert((start = clock())!=-1);
  
  rk = square_norm(Q, N, 1);
  sourcesquarenorm = rk; // for relative precision
  assign(g_spinor_field[DUM_SOLVER],Q,N);
  printf("Initial residue: %.16e\n",rk);
  zero_spinor_field(g_spinor_field[DUM_SOLVER+1],  N);//spin2 = x_k
  zero_spinor_field(g_spinor_field[DUM_SOLVER+2],  N);
  printf("The VOLUME is: %d\n",N);
  
  
  
for(iter=0; iter<max_iter; iter++){

   printf("Applying double precision Dirac-Op...\n");
   
   Q_pm_psi_gpu(g_spinor_field[DUM_SOLVER+3], g_spinor_field[DUM_SOLVER+2]);
   diff(g_spinor_field[DUM_SOLVER],g_spinor_field[DUM_SOLVER],g_spinor_field[DUM_SOLVER+3],N);
    // r_k = b - D x_k
   
   rk = square_norm(g_spinor_field[DUM_SOLVER], N, 0);
  
   #ifdef GF_8
    if(isnan(rk)){
      fprintf(stderr, "Error in mixed_solve_eo: Residue is NaN.\n  May happen with GF 8 reconstruction. Aborting ...\n");
      exit(200);
    }
   #endif
   
   printf("Residue after %d inner solver iterations: %.18e\n",outercount,rk);
   if(((rk <= eps) && (rel_prec == 0)) || ((rk <= eps*sourcesquarenorm) && (rel_prec == 1)))
   {
     printf("Reached solver precision of eps=%.2e\n",eps);
     //multiply with D^dagger
     Q_minus_psi_gpu(g_spinor_field[DUM_SOLVER+3], g_spinor_field[DUM_SOLVER+1]);
     assign(P, g_spinor_field[DUM_SOLVER+3], N);
  

    stop = clock();
    timeelapsed = (double) (stop-start)/CLOCKS_PER_SEC;
    printf("Inversion done in mixed precision.\n Number of iterations in outer solver: %d\n Squared residue: %.8e\n Time elapsed: %.6e sec\n", outercount, rk, timeelapsed);
    finalize_mixedsolve();
    return(totalcount);  
   }
   

  //initialize spin fields on device
  convert2REAL4_spin(g_spinor_field[DUM_SOLVER],h2d_spin);
  
  cudaMemcpy(dev_spinin, h2d_spin, dev_spinsize, cudaMemcpyHostToDevice);
  printf("%s\n", cudaGetErrorString(cudaGetLastError()));

   // solve in single prec on device
   // D p_k = r_k
   printf("Entering inner solver\n");
   assert((startinner = clock())!=-1);
   totalcount += dev_cg(dev_gf, dev_spinin, dev_spinout, dev_spin1, dev_spin2, dev_spin3, dev_spin4, dev_spin5, dev_grid,dev_nn, dev_output,NULL, T, LZ,0);
   stopinner = clock();
   timeelapsed = (double) (stopinner-startinner)/CLOCKS_PER_SEC;
   printf("Inner solver done\nTime elapsed: %.6e sec\n", timeelapsed);
   
  
   // copy back
   cudaMemcpy(h2d_spin, dev_spinout, dev_spinsize, cudaMemcpyDeviceToHost);
   printf("%s\n", cudaGetErrorString(cudaGetLastError()));
   
   convert2double_spin(h2d_spin, g_spinor_field[DUM_SOLVER+2]);
   
   add(g_spinor_field[DUM_SOLVER+1],g_spinor_field[DUM_SOLVER+1],g_spinor_field[DUM_SOLVER+2],N);
   // x_(k+1) = x_k + p_k
   
   outercount ++;
    
}// outer loop 

     printf("Did NOT reach solver precision of eps=%.2e\n",eps);
     //multiply with D^dagger
     Q_minus_psi_gpu(g_spinor_field[DUM_SOLVER+3], g_spinor_field[DUM_SOLVER+1]);
     assign(P, g_spinor_field[DUM_SOLVER+3], N);
  

    stop = clock();
    timeelapsed = (double) (stop-start)/CLOCKS_PER_SEC;
    printf("Inversion done in mixed precision.\n Number of iterations in outer solver: %d\n Squared residue: %.8e\n Time elapsed: %.6e sec\n", outercount, rk, timeelapsed);

  return(-1);
}






void benchmark(spinor * const Q){
  
  double timeelapsed = 0.0;
  clock_t start, stop;
  int i;
  
  size_t dev_spinsize = 6*VOLUME/2 * sizeof(dev_spinor); // float4 even-odd !
  convert2REAL4_spin(Q,h2d_spin);
  cudaMemcpy(dev_spinin, h2d_spin, dev_spinsize, cudaMemcpyHostToDevice);
  printf("%s\n", cudaGetErrorString(cudaGetLastError()));
  assert((start = clock())!=-1);

 #ifdef USETEXTURE
  //Bind texture gf
  bind_texture_gf(dev_gf);
 #endif

 //Initialize some stuff
  printf("mu = %f\n", g_mu);
  dev_complex h0,h1,h2,h3,mh0, mh1, mh2, mh3;
  h0.re = (REAL)ka0.re;    h0.im = -(REAL)ka0.im;
  h1.re = (REAL)ka1.re;    h1.im = -(REAL)ka1.im;
  h2.re = (REAL)ka2.re;    h2.im = -(REAL)ka2.im;
  h3.re = (REAL)ka3.re;    h3.im = -(REAL)ka3.im;
  
  mh0.re = -(REAL)ka0.re;    mh0.im = (REAL)ka0.im;
  mh1.re = -(REAL)ka1.re;    mh1.im = (REAL)ka1.im;
  mh2.re = -(REAL)ka2.re;    mh2.im = (REAL)ka2.im;
  mh3.re = -(REAL)ka3.re;    mh3.im = (REAL)ka3.im;
  
  // try using constant mem for kappas
  cudaMemcpyToSymbol("dev_k0c", &h0, sizeof(dev_complex)) ; 
  cudaMemcpyToSymbol("dev_k1c", &h1, sizeof(dev_complex)) ; 
  cudaMemcpyToSymbol("dev_k2c", &h2, sizeof(dev_complex)) ; 
  cudaMemcpyToSymbol("dev_k3c", &h3, sizeof(dev_complex)) ;
  
  cudaMemcpyToSymbol("dev_mk0c", &mh0, sizeof(dev_complex)) ; 
  cudaMemcpyToSymbol("dev_mk1c", &mh1, sizeof(dev_complex)) ; 
  cudaMemcpyToSymbol("dev_mk2c", &mh2, sizeof(dev_complex)) ; 
  cudaMemcpyToSymbol("dev_mk3c", &mh3, sizeof(dev_complex)) ;  
  
  
  int blockdim3=BLOCK;
  int gridsize;
  if( VOLUME/2 >= BLOCK){
    gridsize = (int)(VOLUME/2/BLOCK) + 1;
  }
  else{
    gridsize=1;
  }
 printf("gridsize = %d\n", gridsize);
 int griddim3=gridsize; 
  
  
  he_cg_init<<< 1, 1 >>> (dev_grid, (REAL) g_kappa, (REAL)(g_mu/(2.0*g_kappa)), h0,h1,h2,h3); 
  printf("%s\n", cudaGetErrorString(cudaGetLastError()));
  printf("Applying H 1000 times\n");
  for(i=0; i<1000; i++){
      #ifdef USETEXTURE
       bind_texture_spin(dev_spinin,1);
      #endif
       //bind_texture_nn(dev_nn_eo);
      //cudaFuncSetCacheConfig(dev_Hopping_Matrix, cudaFuncCachePreferL1);
      dev_Hopping_Matrix<<<griddim3, blockdim3>>>
             (dev_gf, dev_spinin, dev_spin_eo1, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0); //dev_spin_eo1 == even -> 0
       //unbind_texture_nn();
    #ifdef USETEXTURE             
     unbind_texture_spin(1);
    #endif

    #ifdef USETEXTURE
     bind_texture_spin(dev_spin_eo1,1);
    #endif
  //bind_texture_nn(dev_nn_oe);
   // cudaFuncSetCacheConfig(dev_Hopping_Matrix, cudaFuncCachePreferL1);
    dev_Hopping_Matrix<<<griddim3, blockdim3>>>
            (dev_gf, dev_spin_eo1, dev_spinin, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1); 
  //unbind_texture_nn();
    #ifdef USETEXTURE
     unbind_texture_spin(1);
    #endif

  }  
  printf("%s\n", cudaGetErrorString(cudaGetLastError())); 
  printf("Done\n"); 
  
  assert((stop = clock())!=-1);
  timeelapsed = (double) (stop-start)/CLOCKS_PER_SEC;
  // x2 because 2x Hopping per iteration
  double benchres = 1400.0*2*(VOLUME/2)* 1000 / timeelapsed / 1.0e9;
  printf("Benchmark: %f Gflops\n", benchres); 
   
  #ifdef USETEXTURE
    unbind_texture_gf();
  #endif
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
  int iter, retval;
  

  size_t dev_spinsize = 6*VOLUME/2 * sizeof(dev_spinor); // float4 even-odd !
    
  init_mixedsolve_eo(g_gauge_field);
  
  /*
  // small benchmark
    assign(g_spinor_field[DUM_SOLVER],Q,N);
    benchmark(g_spinor_field[DUM_SOLVER]);
  // end small benchmark
  
  exit(100);
  */
 
 


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
  assign(g_spinor_field[DUM_SOLVER],Q,N);
  printf("Initial residue: %.16e\n",rk);
  zero_spinor_field(g_spinor_field[DUM_SOLVER+1],  N);//spin2 = x_k
  zero_spinor_field(g_spinor_field[DUM_SOLVER+2],  N);
  printf("The VOLUME/2 is: %d\n",N);
  
for(iter=0; iter<max_iter; iter++){

   printf("Applying double precision EO Dirac-Op Q_{-}Q{+}...\n");
   
   Qtm_pm_psi(g_spinor_field[DUM_SOLVER+3], g_spinor_field[DUM_SOLVER+2]);
   diff(g_spinor_field[DUM_SOLVER],g_spinor_field[DUM_SOLVER],g_spinor_field[DUM_SOLVER+3],N);
    // r_k = b - D x_k
   
   rk = square_norm(g_spinor_field[DUM_SOLVER], N, 0);
   #ifdef GF_8
    if(isnan(rk)){
      fprintf(stderr, "Error in mixed_solve_eo: Residue is NaN.\n  May happen with GF 8 reconstruction. Aborting ...\n");
      exit(200);
    }
   #endif
   
   printf("Residue after %d inner solver iterations: %.18e\n",outercount,rk);
   
   if(((rk <= eps) && (rel_prec == 0)) || ((rk <= eps*sourcesquarenorm) && (rel_prec == 1)))
   {
     printf("Reached solver precision of eps=%.2e\n",eps);
     //multiply with Qtm_minus_psi (for non gpu done in invert_eo.c)
     Qtm_minus_psi(g_spinor_field[DUM_SOLVER+3], g_spinor_field[DUM_SOLVER+1]);
     assign(P, g_spinor_field[DUM_SOLVER+3], N);

     printf("EO Inversion done in mixed precision.\n Number of iterations in outer solver: %d\n Squared residue: %.8e\n Time elapsed: %.6e sec\n", outercount, rk, timeelapsed);
   
  
     stop = clock();
     timeelapsed = (double) (stop-start)/CLOCKS_PER_SEC;
        
     finalize_mixedsolve();
     return(totalcount);  
   }
   
  //initialize spin fields on device
  convert2REAL4_spin(g_spinor_field[DUM_SOLVER],h2d_spin);
  
  cudaMemcpy(dev_spinin, h2d_spin, dev_spinsize, cudaMemcpyHostToDevice);
  printf("%s\n", cudaGetErrorString(cudaGetLastError()));

   // solve in single prec on device
   // D p_k = r_k
   printf("Entering inner solver\n");
   assert((startinner = clock())!=-1);
   totalcount += dev_cg_eo(dev_gf, dev_spinin, dev_spinout, dev_spin1, dev_spin2, dev_spin3, dev_spin4, dev_spin5, dev_grid,dev_nn, dev_output,NULL, T, LZ,0,(REAL) finaleps);
   stopinner = clock();
   timeelapsed = (double) (stopinner-startinner)/CLOCKS_PER_SEC;
   printf("Inner solver done\nTime elapsed: %.6e sec\n", timeelapsed);
 
   // copy back
   cudaMemcpy(h2d_spin, dev_spinout, dev_spinsize, cudaMemcpyDeviceToHost);
   printf("%s\n", cudaGetErrorString(cudaGetLastError()));
   
   convert2double_spin(h2d_spin, g_spinor_field[DUM_SOLVER+2]);
   // x_(k+1) = x_k + p_k
   add(g_spinor_field[DUM_SOLVER+1],g_spinor_field[DUM_SOLVER+1],g_spinor_field[DUM_SOLVER+2],N);

   outercount ++;   
}// outer loop 
    
     printf("Did NOT reach solver precision of eps=%.2e\n",eps);
     //multiply with Qtm_minus_psi (for non gpu done in invert_eo.c)
     Qtm_minus_psi(g_spinor_field[DUM_SOLVER+3], g_spinor_field[DUM_SOLVER+1]);
     assign(P, g_spinor_field[DUM_SOLVER+3], N);
    

    assert((stop = clock())!=-1);
    timeelapsed = (double) (stop-start)/CLOCKS_PER_SEC;
    printf("Inversion done in mixed precision.\n Number of iterations in outer solver: %d\n Squared residue: %.8e\n Time elapsed: %.6e sec\n", outercount, rk, timeelapsed);

  finalize_mixedsolve();
  return(-1);
}



// mixed solver, even/odd, non-degenerate two flavour
#include "mixed_solve_eo_nd.cuh"


#if defined(MPI) && defined(PARALLELT)
  //#include "./DEBUG/MATRIX_MPI_DEBUG.cuh"		// for debugging
  #include "MPI.cuh"
#endif






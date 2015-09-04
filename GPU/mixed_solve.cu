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
#include "../su3spinor.h"
#include "../solver/solver_field.h"

#ifdef MPI
  #include "../xchange.h"
#endif 

}



#ifdef HAVE_CONFIG_H
  #include<config.h>
#endif


#ifdef MPI
  #undef MPI
  #undef REAL
    #include <mpi.h>
  #define MPI
  #define REAL float
#endif

#include "MACROS.cuh"
#include "cublasWrapper.cuh"



int g_numofgpu;

template<class RealT=REAL>
struct MixedsolveParameter
{//internal variables of mixed solver routine corresponding to fields on the device
  #ifdef GF_8
    dev_su3_8M(RealT)* dev_gf;
    dev_su3_8M(RealT)* h2d_gf;
  #else
    dev_su3_2vM(RealT)* dev_gf;
    dev_su3_2vM(RealT)* h2d_gf;
  #endif

  #ifndef HALF
    dev_spinorM(RealT)* dev_spin1;
    dev_spinorM(RealT)* dev_spin2;
    dev_spinorM(RealT)* dev_spin3;
    dev_spinorM(RealT)* dev_spin4;
    dev_spinorM(RealT)* dev_spin5;
    dev_spinorM(RealT)* dev_spinin;
    dev_spinorM(RealT)* dev_spinout;
    dev_spinorM(RealT)* h2d_spin;

    //additional spinors for even-odd
    dev_spinorM(RealT)* dev_spin_eo1;
    dev_spinorM(RealT)* dev_spin_eo2;
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
    dev_spinor_half* dev_spin_eo1;
    dev_spinor_half* dev_spin_eo2;


    RealT* dev_spin1_norm;
    RealT* dev_spin2_norm;
    RealT* dev_spin3_norm;
    RealT* dev_spin4_norm;
    RealT* dev_spin5_norm;
    RealT* dev_spinin_norm;
    RealT* dev_spinout_norm;
    RealT* h2d_spin_norm;

    RealT* dev_spin_eo1_norm;
    RealT* dev_spin_eo2_norm;


    // a half precsion gauge field
    #ifdef GF_8
       dev_su3_8_half* dev_gf_half;
    #else
       dev_su3_2v_half* dev_gf_half;
    #endif
  #endif 



  // selects global instance of this structure depending on template parameter RealT to determine precision
  static MixedsolveParameter<RealT>* getGlobalP();
};
MixedsolveParameter<REAL > mixedsolveParameter ;
MixedsolveParameter<REALD> mixedsolveParameterD;

template<class RealT> inline MixedsolveParameter<RealT>* MixedsolveParameter<RealT>::getGlobalP() { printf("WARNING: MixedsolveParameter::getGlobal() called with invalid template argument.\n"); return NULL; }
template<           > inline MixedsolveParameter<REAL >* MixedsolveParameter<REAL >::getGlobalP() { return &mixedsolveParameter ; }
template<           > inline MixedsolveParameter<REALD>* MixedsolveParameter<REALD>::getGlobalP() { return &mixedsolveParameterD; }


//{
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
  float* dev_output;


  REALD hostr;
  REALD hostkappa;
  REALD hostm;
  REALD hostmu;
//}


int havedevice = 0;


__device__  REAL m;
__device__  REAL mu;
__device__  REAL r=1.0; // this is implicitly assumed to be 1.0 in the host code!!!
__device__  REAL kappa;
__device__  REAL twokappamu;

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




#ifdef MPI


// from mixed_solve_eo_nd.cuh
__device__ int dev_RAND;                        // not used, maybe later ...
__device__ int dev_VOLUMEPLUSRAND;              // is now used in dev_Hopping_Matrix_mpi()
__device__ int dev_rank;
__device__ int dev_nproc;


  #ifndef ALTERNATE_FIELD_XCHANGE
    spinor * spinor_xchange;                    // for xchange_field_wrapper()
  #else
    dev_spinor * R1;
    dev_spinor * R2;
    dev_spinor * R3;
    dev_spinor * R4;
  #endif


#if ASYNC > 0
    int nStreams = ASYNC_OPTIMIZED;
    cudaStream_t stream[2*ASYNC_OPTIMIZED+1];

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
#endif



#if defined(ALTERNATE_FIELD_XCHANGE) || defined(ASYNC_OPTIMIZED)
  MPI_Status stat[2];
  MPI_Request send_req[2];
  MPI_Request recv_req[2]; 
#endif


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



// the device Hopping_Matrix
#include "Hopping_Matrix.cuh"
// the non-EO twisted mass dirac operator
#include "tm_diracoperator.cuh"
// mixed solver, even/odd, non-degenerate two flavour
#include "mixed_solve_eo_nd.cuh"

#ifdef MPI
// optimization of the communication
  #include "ASYNC.cuh"
#endif 





#ifndef HALF
// computes sout = 1/(1 +- mutilde gamma5) sin = (1 -+ i mutilde gamma5)/(1+mutilde^2) sin
// mutilde = 2 kappa mu
// uses shared local memory for manipulation
template<class RealT>
__global__ void dev_mul_one_pm_imu_inv(dev_spinorM(RealT)* sin, dev_spinorM(RealT)* sout, const RealT sign){
   dev_spinorM(RealT) slocal[6];
   //need the inverse sign in the numerator because of inverse
   dev_complexM(RealT) pm_imu = dev_initcomplex<RealT>(0.0,-1.0*sign*twokappamu);
   
   RealT one_plus_musquare_inv = 1.0/(1.0 + twokappamu*twokappamu);
   int pos;
   pos= threadIdx.x + blockDim.x*blockIdx.x;
   //not referenced: int ix = threadIdx.x;
   if(pos < dev_VOLUME){
     //dev_skalarmult_spinor(&(sin[6*pos]), pm_imu, &(slocal[0]));
     //dev_Gamma5(&(slocal[0]));
     dev_skalarmult_gamma5_spinor(&(slocal[0]), pm_imu, &(sin[6*pos]) );
	 dev_add_spinor_assign<RealT>(&(slocal[0]), &(sin[6*pos]));
     //dev_realmult_spinor(&(slocal[0]), one_plus_musquare_inv);
     //dev_copy_spinor(&(slocal[0]), &(sout[6*pos])); 
     dev_realmult_spinor_assign(&(sout[6*pos]), one_plus_musquare_inv, &(slocal[0]) );
   }
}





// sout = gamma_5*((1\pm i\mutilde \gamma_5)*sin1 - sin2)
// uses shared local memory for manipulation
template<class RealT>
__global__ void dev_mul_one_pm_imu_sub_mul_gamma5(dev_spinorM(RealT)* sin1, dev_spinorM(RealT)* sin2, dev_spinorM(RealT)* sout, const RealT sign){
   dev_spinorM(RealT) slocal[6];
   dev_complexM(RealT) pm_imu = dev_initcomplex<RealT>(0.0, sign*twokappamu); // i mutilde
   int pos;
   pos= threadIdx.x + blockDim.x*blockIdx.x; 
   //not referenced: int ix = threadIdx.x;
   if(pos < dev_VOLUME){
     //dev_skalarmult_spinor(&(sin1[6*pos]), pm_imu, &(slocal[0]));
     //dev_Gamma5(&(slocal[0]));
     dev_skalarmult_gamma5_spinor(&(slocal[0]),pm_imu,&(sin1[6*pos]));
	 dev_add_spinor_assign<RealT>(&(slocal[0]), &(sin1[6*pos]));
     dev_sub_spinor_assign<RealT>(&(slocal[0]), &(sin2[6*pos]));
     //dev_Gamma5(&(slocal[0]));
     //dev_copy_spinor(&(slocal[0]), &(sout[6*pos]));
     dev_Gamma5_assign<RealT>(&(sout[6*pos]), &(slocal[0]));
   }
}













// aequivalent to Qtm_pm_psi in tm_operators.c
template<class RealT>
void dev_Qtm_pm_psi(dev_spinorM(RealT)* spinin, dev_spinorM(RealT)* spinout, int gridsize, int blocksize, int gridsize2, int blocksize2, MixedsolveParameter<RealT>& mixedsolveParameter){
  //spinin == odd
  //spinout == odd
  
  //Q_{-}
  		#ifdef MPI
  		  xchange_field_wrapper(spinin, 0);
  		#endif
  #ifdef USETEXTURE
    bind_texture_spin(spinin,1);
  #endif
  //bind_texture_nn(dev_nn_eo);
  //cudaFuncSetCacheConfig(dev_Hopping_Matrix, cudaFuncCachePreferL1);
    dev_Hopping_Matrix<RealT> <<<gridsize, blocksize>>>
             (mixedsolveParameter.dev_gf, spinin, mixedsolveParameter.dev_spin_eo1, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0); //mixedsolveParameter.dev_spin_eo1 == even -> 0           
  //unbind_texture_nn();           
  #ifdef USETEXTURE
    unbind_texture_spin(1);
  #endif
  dev_mul_one_pm_imu_inv<RealT> <<<gridsize2, blocksize2>>>(mixedsolveParameter.dev_spin_eo1,mixedsolveParameter.dev_spin_eo2, -1.);
  
  		#ifdef MPI
  		  xchange_field_wrapper(mixedsolveParameter.dev_spin_eo2, 1);
  		#endif
  #ifdef USETEXTURE
    bind_texture_spin(mixedsolveParameter.dev_spin_eo2,1);
  #endif
  //bind_texture_nn(dev_nn_oe);
  //cudaFuncSetCacheConfig(dev_Hopping_Matrix, cudaFuncCachePreferL1);
    dev_Hopping_Matrix<RealT> <<<gridsize, blocksize>>>
            (mixedsolveParameter.dev_gf, mixedsolveParameter.dev_spin_eo2, mixedsolveParameter.dev_spin_eo1, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1); 
  //unbind_texture_nn();
  #ifdef USETEXTURE
    unbind_texture_spin(1);
  #endif
  dev_mul_one_pm_imu_sub_mul_gamma5<RealT> <<<gridsize2, blocksize2>>>(spinin, mixedsolveParameter.dev_spin_eo1,  mixedsolveParameter.dev_spin_eo2, -1.);
  
  
  //Q_{+}
  		#ifdef MPI
  		  xchange_field_wrapper(mixedsolveParameter.dev_spin_eo2, 0);
  		#endif
  #ifdef USETEXTURE
    bind_texture_spin(mixedsolveParameter.dev_spin_eo2,1);
  #endif
  //bind_texture_nn(dev_nn_eo);
  //cudaFuncSetCacheConfig(dev_Hopping_Matrix, cudaFuncCachePreferL1);
    dev_Hopping_Matrix<RealT> <<<gridsize, blocksize>>>
          (mixedsolveParameter.dev_gf, mixedsolveParameter.dev_spin_eo2, mixedsolveParameter.dev_spin_eo1, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0); //mixedsolveParameter.dev_spin_eo1 == even -> 0
  //unbind_texture_nn();      
  #ifdef USETEXTURE  
    unbind_texture_spin(1);
  #endif
  dev_mul_one_pm_imu_inv<RealT> <<<gridsize2, blocksize2>>>(mixedsolveParameter.dev_spin_eo1,spinout, +1.);
  
  		#ifdef MPI
  		  xchange_field_wrapper(spinout, 1);
  		#endif
  #ifdef USETEXTURE
    bind_texture_spin(spinout,1);
  #endif
  //bind_texture_nn(dev_nn_oe);
  //cudaFuncSetCacheConfig(dev_Hopping_Matrix, cudaFuncCachePreferL1);
    dev_Hopping_Matrix<RealT> <<<gridsize, blocksize>>>
             (mixedsolveParameter.dev_gf, spinout, mixedsolveParameter.dev_spin_eo1, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1); 
  //unbind_texture_nn();  
  #ifdef USETEXTURE
    unbind_texture_spin(1);
  #endif
  dev_mul_one_pm_imu_sub_mul_gamma5<RealT> <<<gridsize2, blocksize2>>>(mixedsolveParameter.dev_spin_eo2, mixedsolveParameter.dev_spin_eo1,  spinout , +1.); 
}




#ifdef MPI
// aequivalent to Qtm_pm_psi in tm_operators.c
// using HOPPING_ASYNC for mpi
template<class RealT>
void dev_Qtm_pm_psi_mpi(dev_spinorM(RealT)* spinin, dev_spinorM(RealT)* spinout, int gridsize, int blocksize, int gridsize2, int blocksize2,MixedsolveParameter<RealT>& mixedsolveParameter){
  //spinin == odd
  //spinout == odd
  
  //Q_{-}

  //cudaFuncSetCacheConfig(dev_Hopping_Matrix, cudaFuncCachePreferL1);
    HOPPING_ASYNC(mixedsolveParameter.dev_gf, spinin, mixedsolveParameter.dev_spin_eo1, dev_eoidx_even, 
               dev_eoidx_odd, dev_nn_eo, 0,gridsize, blocksize); //mixedsolveParameter.dev_spin_eo1 == even -> 0           
          


  dev_mul_one_pm_imu_inv<RealT> <<<gridsize2, blocksize2>>>(mixedsolveParameter.dev_spin_eo1,mixedsolveParameter.dev_spin_eo2, -1.);
  



  //cudaFuncSetCacheConfig(dev_Hopping_Matrix, cudaFuncCachePreferL1);
    HOPPING_ASYNC(mixedsolveParameter.dev_gf, mixedsolveParameter.dev_spin_eo2, mixedsolveParameter.dev_spin_eo1, 
          dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1,gridsize, 
          blocksize); 

  dev_mul_one_pm_imu_sub_mul_gamma5<RealT> <<<gridsize2, blocksize2>>>(spinin, mixedsolveParameter.dev_spin_eo1,  mixedsolveParameter.dev_spin_eo2, -1.);
  
  
  //Q_{+}

  //cudaFuncSetCacheConfig(dev_Hopping_Matrix, cudaFuncCachePreferL1);
    HOPPING_ASYNC(mixedsolveParameter.dev_gf, mixedsolveParameter.dev_spin_eo2, mixedsolveParameter.dev_spin_eo1, 
         dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0, gridsize, 
         blocksize); //mixedsolveParameter.dev_spin_eo1 == even -> 0

  dev_mul_one_pm_imu_inv<RealT> <<<gridsize2, blocksize2>>>(mixedsolveParameter.dev_spin_eo1,spinout, +1.);
  

  //cudaFuncSetCacheConfig(dev_Hopping_Matrix, cudaFuncCachePreferL1);
    HOPPING_ASYNC(mixedsolveParameter.dev_gf, spinout, mixedsolveParameter.dev_spin_eo1, dev_eoidx_odd, 
           dev_eoidx_even, dev_nn_oe, 1,gridsize, blocksize); 

  dev_mul_one_pm_imu_sub_mul_gamma5<RealT> <<<gridsize2, blocksize2>>>(mixedsolveParameter.dev_spin_eo2, mixedsolveParameter.dev_spin_eo1,  spinout , +1.); 
}
#endif




#else // HALF

// computes sout = 1/(1 +- mutilde gamma5) sin = (1 -+ i mutilde gamma5)/(1+mutilde^2) sin
// mutilde = 2 kappa mu
// uses shared local memory for manipulation
__global__ void dev_mul_one_pm_imu_inv_half(dev_spinor_half* sin, float* sin_norm, dev_spinor_half* sout, float* sout_norm, const REAL sign){
   
   typedef REAL RealT;
   dev_spinor slocal[6];
   dev_spinor s[6];
   float norm;
   
   //need the inverse sign in the numerator because of inverse
   dev_complex pm_imu = dev_initcomplex<RealT>(0.0,-1.0*sign*twokappamu);
   
   REAL one_plus_musquare_inv = 1.0/(1.0 + twokappamu*twokappamu);
   int pos;
   pos= threadIdx.x + blockDim.x*blockIdx.x;  
   int ix = threadIdx.x;
   if(pos < dev_VOLUME){
     norm = sin_norm[pos];
     construct_spinor_fromhalf(s, sin, norm, pos);

     dev_skalarmult_gamma5_spinor(&(slocal[0]), pm_imu, &(s[0]) );
         dev_add_spinor_assign<RealT>(&(slocal[0]), &(s[0]));
     
     dev_realmult_spinor_assign(&(s[0]), one_plus_musquare_inv, &(slocal[0]) );
     
     dev_write_spinor_half(&(s[0]),&(sout[6*pos]), &(sout_norm[pos]));
   }
}





// sout = gamma_5*((1\pm i\mutilde \gamma_5)*sin1 - sin2)
// uses shared local memory for manipulation
__global__ void dev_mul_one_pm_imu_sub_mul_gamma5_half(dev_spinor_half* sin1, float* sin1_norm, dev_spinor_half* sin2, float* sin2_norm, dev_spinor_half* sout, float* sout_norm, const REAL sign){
   typedef REAL RealT;
   dev_spinor slocal[6];
   dev_spinor s1[6];
   dev_spinor s2[6];
   float norm;
   dev_complex pm_imu = dev_initcomplex<RealT>(0.0, sign*twokappamu); // i mutilde
   int pos;
   pos= threadIdx.x + blockDim.x*blockIdx.x; 
   int ix = threadIdx.x;
   if(pos < dev_VOLUME){
     norm = sin1_norm[pos];
     construct_spinor_fromhalf(s1, sin1,norm, pos);
     norm = sin2_norm[pos];
     construct_spinor_fromhalf(s2, sin2, norm, pos);

     dev_skalarmult_gamma5_spinor(&(slocal[0]),pm_imu,&(s1[0]));
         dev_add_spinor_assign<RealT>(&(slocal[0]), &(s1[0]));
     dev_sub_spinor_assign<RealT>(&(slocal[0]), &(s2[0]));
     dev_Gamma5_assign<RealT>(&(s1[0]), &(slocal[0]));
     dev_write_spinor_half(&(s1[0]),&(sout[6*pos]), &(sout_norm[pos]));
   }   
}





// aequivalent to Qtm_pm_psi in tm_operators.c for half precision
extern "C" void dev_Qtm_pm_psi_half(dev_spinor_half* spinin, float* spinin_norm, dev_spinor_half* spinout, float* spinout_norm, int gridsize, int blocksize, int gridsize2, int blocksize2){
  //spinin == odd
  //spinout == odd
  
  //Q_{-}
  #ifdef USETEXTURE
    bind_halfspinor_texture(spinin, spinin_norm);
  #endif
  //cudaFuncSetCacheConfig(dev_Hopping_Matrix_half, cudaFuncCachePreferL1);
    dev_Hopping_Matrix_half<<<gridsize, blocksize>>>
             (mixedsolveParameter.dev_gf_half, spinin, spinin_norm, mixedsolveParameter.dev_spin_eo1, mixedsolveParameter.dev_spin_eo1_norm, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0); //mixedsolveParameter.dev_spin_eo1 == even -> 0  
  #ifdef USETEXTURE
    unbind_halfspinor_texture();
  #endif
  dev_mul_one_pm_imu_inv_half<<<gridsize2, blocksize2>>>(mixedsolveParameter.dev_spin_eo1, mixedsolveParameter.dev_spin_eo1_norm ,mixedsolveParameter.dev_spin_eo2, mixedsolveParameter.dev_spin_eo2_norm, -1.);
  
  #ifdef USETEXTURE
    bind_halfspinor_texture(mixedsolveParameter.dev_spin_eo2, mixedsolveParameter.dev_spin_eo2_norm);
  #endif
  //cudaFuncSetCacheConfig(dev_Hopping_Matrix_half, cudaFuncCachePreferL1);
    dev_Hopping_Matrix_half<<<gridsize, blocksize>>>
            (mixedsolveParameter.dev_gf_half, mixedsolveParameter.dev_spin_eo2, mixedsolveParameter.dev_spin_eo2_norm, mixedsolveParameter.dev_spin_eo1, mixedsolveParameter.dev_spin_eo1_norm, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1); 
  #ifdef USETEXTURE
    unbind_halfspinor_texture();
  #endif
  dev_mul_one_pm_imu_sub_mul_gamma5_half<<<gridsize2, blocksize2>>>(spinin, spinin_norm, mixedsolveParameter.dev_spin_eo1, mixedsolveParameter.dev_spin_eo1_norm,  mixedsolveParameter.dev_spin_eo2, mixedsolveParameter.dev_spin_eo2_norm, -1.);
  
  //Q_{+}
  #ifdef USETEXTURE
    bind_halfspinor_texture(mixedsolveParameter.dev_spin_eo2, mixedsolveParameter.dev_spin_eo2_norm);
  #endif
  //cudaFuncSetCacheConfig(dev_Hopping_Matrix_half, cudaFuncCachePreferL1);
    dev_Hopping_Matrix_half<<<gridsize, blocksize>>>
          (mixedsolveParameter.dev_gf_half, mixedsolveParameter.dev_spin_eo2, mixedsolveParameter.dev_spin_eo2_norm, mixedsolveParameter.dev_spin_eo1, mixedsolveParameter.dev_spin_eo1_norm, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0); //mixedsolveParameter.dev_spin_eo1 == even -> 0    
  #ifdef USETEXTURE  
    unbind_halfspinor_texture();
  #endif
  dev_mul_one_pm_imu_inv_half<<<gridsize2, blocksize2>>>(mixedsolveParameter.dev_spin_eo1, mixedsolveParameter.dev_spin_eo1_norm,spinout, spinout_norm, +1.);
  
  #ifdef USETEXTURE
    bind_halfspinor_texture(spinout, spinout_norm);
  #endif
  //cudaFuncSetCacheConfig(dev_Hopping_Matrix_half, cudaFuncCachePreferL1);
    dev_Hopping_Matrix_half<<<gridsize, blocksize>>>
             (mixedsolveParameter.dev_gf_half, spinout, spinout_norm, mixedsolveParameter.dev_spin_eo1, mixedsolveParameter.dev_spin_eo1_norm, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1);  
  #ifdef USETEXTURE
    unbind_halfspinor_texture();
  #endif
  dev_mul_one_pm_imu_sub_mul_gamma5_half<<<gridsize2, blocksize2>>>(mixedsolveParameter.dev_spin_eo2, mixedsolveParameter.dev_spin_eo2_norm, mixedsolveParameter.dev_spin_eo1, mixedsolveParameter.dev_spin_eo1_norm,  spinout, spinout_norm , +1.); 
}


#ifdef MPI

// aequivalent to Qtm_pm_psi in tm_operators.c for half precision
extern "C" void dev_Qtm_pm_psi_half_mpi(dev_spinor_half* spinin, float* spinin_norm, dev_spinor_half* spinout, float* spinout_norm, int gridsize, int blocksize, int gridsize2, int blocksize2){
  //spinin == odd
  //spinout == odd
  
  //Q_{-}
  HOPPING_HALF_ASYNC(mixedsolveParameter.dev_gf_half, spinin, spinin_norm, mixedsolveParameter.dev_spin_eo1, mixedsolveParameter.dev_spin_eo1_norm, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0,gridsize, blocksize); //mixedsolveParameter.dev_spin_eo1 == even -> 0  

  dev_mul_one_pm_imu_inv_half<<<gridsize2, blocksize2>>>(mixedsolveParameter.dev_spin_eo1, mixedsolveParameter.dev_spin_eo1_norm ,mixedsolveParameter.dev_spin_eo2, mixedsolveParameter.dev_spin_eo2_norm, -1.);
  

    HOPPING_HALF_ASYNC(mixedsolveParameter.dev_gf_half, mixedsolveParameter.dev_spin_eo2, mixedsolveParameter.dev_spin_eo2_norm, mixedsolveParameter.dev_spin_eo1, mixedsolveParameter.dev_spin_eo1_norm, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1,gridsize, blocksize); 

  dev_mul_one_pm_imu_sub_mul_gamma5_half<<<gridsize2, blocksize2>>>(spinin, spinin_norm, mixedsolveParameter.dev_spin_eo1, mixedsolveParameter.dev_spin_eo1_norm,  mixedsolveParameter.dev_spin_eo2, mixedsolveParameter.dev_spin_eo2_norm, -1.);
  
  //Q_{+}
    HOPPING_HALF_ASYNC (mixedsolveParameter.dev_gf_half, mixedsolveParameter.dev_spin_eo2, mixedsolveParameter.dev_spin_eo2_norm, mixedsolveParameter.dev_spin_eo1, mixedsolveParameter.dev_spin_eo1_norm, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0,gridsize, blocksize); //mixedsolveParameter.dev_spin_eo1 == even -> 0    
    
  dev_mul_one_pm_imu_inv_half<<<gridsize2, blocksize2>>>(mixedsolveParameter.dev_spin_eo1, mixedsolveParameter.dev_spin_eo1_norm,spinout, spinout_norm, +1.);
  
    HOPPING_HALF_ASYNC (mixedsolveParameter.dev_gf_half, spinout, spinout_norm, mixedsolveParameter.dev_spin_eo1, mixedsolveParameter.dev_spin_eo1_norm, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1,gridsize, blocksize);  

  dev_mul_one_pm_imu_sub_mul_gamma5_half<<<gridsize2, blocksize2>>>(mixedsolveParameter.dev_spin_eo2, mixedsolveParameter.dev_spin_eo2_norm, mixedsolveParameter.dev_spin_eo1, mixedsolveParameter.dev_spin_eo1_norm,  spinout, spinout_norm , +1.); 
}
#endif // MPI





/*
extern "C" void dev_Qtm_pm_psi(dev_spinor* spinin, dev_spinor* spinout, int gridsize, int blocksize, int gridsize2, int blocksize2){

  printf("WARNING: dummy function 'dev_Qtm_pm_psi' was called\n");
  
}
*/





#endif //HALF



template<class RealT>
__global__ void dev_zero_spinor_field(typename dev_spinorT<RealT>::type* s1){
  int pos;
  pos= threadIdx.x + blockDim.x*blockIdx.x;  
  if(pos < dev_VOLUME){
          dev_zero_spinor<RealT>(&(s1[6*pos]));
  }
}




template<class RealT1,class RealT2>
__global__ void dev_copy_spinor_field(dev_spinorM(RealT1)* s1, dev_spinorM(RealT2)* s2){
  int pos;
  pos= threadIdx.x + blockDim.x*blockIdx.x;  
  if(pos < dev_VOLUME){
      dev_copy_spinor<RealT1,RealT2>(&(s1[6*pos]),&(s2[6*pos]));
  } 
}



template<class RealT>
__global__ void dev_skalarmult_add_assign_spinor_field(dev_spinorM(RealT)* s1, RealT lambda, dev_spinorM(RealT)* s2, dev_spinorM(RealT)* so){
  int pos;
  pos= threadIdx.x + blockDim.x*blockIdx.x;  
  if(pos < dev_VOLUME){
    dev_skalarmult_add_assign_spinor(&(s1[6*pos]), lambda ,&(s2[6*pos]), &(so[6*pos]) );
  }
}



template<class RealT>
__global__ void dev_skalarmult_spinor_field(dev_spinorM(RealT)* s1, RealT lambda, dev_spinorM(RealT)* so){
  int pos;
  pos= threadIdx.x + blockDim.x*blockIdx.x;  
  if(pos < dev_VOLUME){
    dev_skalarmult_spinor(&(s1[6*pos]), dev_initcomplex<RealT>(lambda,0.0) , &(so[6*pos]) );
  }
}  



template<class RealT>
__global__ void dev_complexmult_spinor_field(dev_spinorM(RealT)* s1, dev_complexM(RealT) lambda, dev_spinorM(RealT)* so){
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
extern "C" int find_devices() {

  int deviceCount, dev;

  cudaGetDeviceCount(&deviceCount);
    
  #ifdef MPI
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
    
    #ifdef MPI 
      }
    #endif
    
    return(deviceCount);
}









extern "C" void test_operator(dev_su3_2v * gf,dev_spinor* spinin, dev_spinor* spinout, 
dev_spinor* spin0, dev_spinor* spin1, dev_spinor* spin2, dev_spinor* spin3, dev_spinor* spin4, int *grid, int * nn_grid, REAL* output,REAL* erg, int xsize, int ysize){

 typedef REAL RealT;
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
  dev_tm_dirac_kappa<RealT> <<<griddim3, blockdim3 >>>(gf, spin4, spinout, nn_grid);

 #ifdef USETEXTURE
  unbind_texture_gf();
  unbind_texture_spin(1);
 #endif
}





// this is the eo version of the device cg inner solver 
// we invert the hermitean D_tm D_tm^{+}


/*
///  member definition of CG-interface class ///

template<class RealT>class MixedsolveOperator // interface class
{
public:
  virtual ~MixedsolveOperator() { }

  virtual void gpuInit  (dev_spinorM(RealT)* spinin,dev_spinorM(RealT)* spinTmp,dev_spinorM(RealT)* spinout,dev_su3_2vM(RealT)* gf,int* dev_nn,const dim3& linAlgGriddim,const dim3& linAlgBlockdim) { }
  virtual void gpu      (dev_spinorM(RealT)* spinin,dev_spinorM(RealT)* spinTmp,dev_spinorM(RealT)* spinout,dev_su3_2vM(RealT)* gf,int* dev_nn,const dim3& linAlgGriddim,const dim3& linAlgBlockdim) =0;
  virtual void gpuDeinit(dev_spinorM(RealT)* spininout,dev_spinorM(RealT)* spinTmp,dev_su3_2vM(RealT)* gf,int* dev_nn,const dim3& linAlgGriddim,const dim3& linAlgBlockdim,const RealT scaleparam)   { }

  virtual void checkInit  (spinor* const spinin,spinor* const spinTmp,spinor* const spinout,int Volume)                                 { }
  virtual void check      (spinor* const conjungateBasisPSpinin,spinor* const spinout,const int Volume) =0; 
  virtual void checkDeinit(spinor* const spinin,spinor* const spinTmp,spinor* const spinout,int Volume)                                 { }
};
*/
template<class RealT>MixedsolveOperator<RealT>::~MixedsolveOperator() { }

template<class RealT>void MixedsolveOperator<RealT>::gpuInit  (dev_spinorM(RealT)* spinin,dev_spinorM(RealT)* spinTmp,dev_spinorM(RealT)* spinout,dev_su3_2vM(RealT)* gf,int* dev_nn,const dim3& linAlgGriddim,const dim3& linAlgBlockdim) { }
template<class RealT>void MixedsolveOperator<RealT>::gpuDeinit(dev_spinorM(RealT)* spininout,dev_spinorM(RealT)* spinTmp,dev_su3_2vM(RealT)* gf,int* dev_nn,const dim3& linAlgGriddim,const dim3& linAlgBlockdim,const RealT scaleparam)   { }

template<class RealT>void MixedsolveOperator<RealT>::checkInit  (spinor* const spinin,spinor* const spinTmp,spinor* const spinout,int Volume) { }
template<class RealT>void MixedsolveOperator<RealT>::checkDeinit(spinor* const spinin,spinor* const spinTmp,spinor* const spinout,int Volume) { }



template<class RealT,template<class MixedsolveOperatorSRealT>class MixedsolveOperatorT>
int dev_cg(
       dev_su3_2vM(RealT)* gf,
       dev_spinorM(RealT)* spinin, 
       dev_spinorM(RealT)* spinout, 
       dev_spinorM(RealT)* spin0, 
       dev_spinorM(RealT)* spin1, 
       dev_spinorM(RealT)* spin2, 
       dev_spinorM(RealT)* spin3, 
       dev_spinorM(RealT)* spin4, 
       int* grid, int* nn_grid, MixedsolveOperatorT<RealT>& mixedsolveOperator,
       REALD initial_sourcesquarenorm,bool rel_prec,double finalEps/*,bool& reachedFinalPrecision*/){
 

 RealT host_alpha, host_beta, host_dotprod, host_rk, sourcesquarenorm;
 RealT * dotprod, * dotprod2, * rk, * alpha, *beta;
 
 
 cudaError_t cudaerr;
 int i, gridsize;
 int maxit = max_innersolver_it;
 RealT eps = (RealT) innersolver_precision;
 int N_recalcres = 30; // after N_recalcres iterations calculate r = A x_k - b
 
 
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

 mixedsolveOperator.gpuInit(spin2,spin4,spin3,gf,dev_nn,griddim2,blockdim2);
 
 
 //Initialize some stuff
  printf("mu = %f\n", g_mu);
  dev_complexM(RealT) h0,h1,h2,h3;
  h0.re = (RealT)ka0.re;    h0.im = (RealT)ka0.im;
  h1.re = (RealT)ka1.re;    h1.im = (RealT)ka1.im;
  h2.re = (RealT)ka2.re;    h2.im = (RealT)ka2.im;
  h3.re = (RealT)ka3.re;    h3.im = (RealT)ka3.im;
  he_cg_init<<< 1, 1 >>> (grid, (RealT) g_kappa, (RealT)(g_mu/(2.0*g_kappa)), h0,h1,h2,h3);
  // BEWARE in dev_tm_dirac_kappa we need the true mu (not 2 kappa mu!)
 
 // Init x,p,r for k=0
 // Allocate some numbers for host <-> device interaction
 cudaMalloc((void **) &dotprod, sizeof(RealT));
 cudaMalloc((void **) &dotprod2, sizeof(RealT));
 cudaMalloc((void **) &rk, sizeof(RealT));
 cudaMalloc((void **) &alpha, sizeof(RealT));
 cudaMalloc((void **) &beta, sizeof(RealT));
 printf("%s\n", cudaGetErrorString(cudaGetLastError())); 
 
 
 //init blas
 cublasInit();
 printf("%s\n", cudaGetErrorString(cudaGetLastError())); 
 printf("have initialized cublas\n");
 
 
 // go over to kappa (if wanted)
 RealT scaleparam = sqrt(1.0/(2.0 * (RealT)hostkappa));
 printf("1/2kappa = %.16f\n",scaleparam);
 //dev_skalarmult_spinor_field<<<griddim2, blockdim2 >>>(spinin,scaleparam, spin1);
 //dev_copy_spinor_field<<<griddim2, blockdim2 >>>(spin1, spinin);
 
 
 dev_copy_spinor_field<RealT,RealT> <<<griddim2, blockdim2 >>>(spinin, spin0);
 dev_zero_spinor_field<RealT>       <<<griddim2, blockdim2 >>>(spin1); // x_0 = 0
 dev_copy_spinor_field<RealT,RealT> <<<griddim2, blockdim2 >>>(spinin, spin2);
 dev_zero_spinor_field<RealT>       <<<griddim2, blockdim2 >>>(spin3);
 printf("%s\n", cudaGetErrorString(cudaGetLastError()));
 
 
 
 
 //relative precision -> get initial residue
 sourcesquarenorm = cublasDot (24*VOLUME, (const RealT*)spinin, 1, (const RealT*)spinin, 1);
 host_rk = sourcesquarenorm; //for use in main loop
 printf("Squarenorm Source:\t%.16e\n", sourcesquarenorm);
 printf("%s\n", cudaGetErrorString(cudaGetLastError()));
 
 printf("Entering cg-loop\n");
 for(i=0;i<maxit;i++)
 { //MAIN LOOP
   mixedsolveOperator.gpu(spin2,spin4,spin3,gf,dev_nn,griddim2,blockdim2);// call func
  
  //alpha  /// (spin2=spinin) * (spin3=D_dagger D)
   host_dotprod = cublasDot (24*VOLUME, (const RealT*) spin2, 1, (const RealT*) spin3, 1);
   host_alpha = (host_rk / host_dotprod);
   
  //r(k+1)  /// spin0:=-alpha*(spin3=D_dagger D) + (spin0=spinin(if i==1))
   cublasAxpy (24*VOLUME,-1.0*host_alpha, (const RealT*)spin3, 1, (RealT*)spin0, 1);

  //x(k+1);
   cublasAxpy (24*VOLUME, host_alpha, (const RealT*)spin2,  1, (RealT*)spin1, 1);

   if((cudaerr=cudaGetLastError()) != cudaSuccess)
   {
     printf("%s\n", cudaGetErrorString(cudaerr));
     exit(200);
   }


  //Abbruch?
   host_dotprod = cublasDot (24*VOLUME, (const RealT*)spin0, 1,(const RealT*)spin0, 1);
  
   if ((host_dotprod <= eps*sourcesquarenorm))//error-limit erreicht
     break; 
   printf("iter %d: err = %.16e\n", i, host_dotprod);
  
  //beta
   host_beta =host_dotprod/host_rk;
  //p(k+1)
   cublasScal (24*VOLUME, host_beta, (RealT*)spin2, 1);
   cublasAxpy (24*VOLUME, 1.0, (const RealT*)spin0, 1, (RealT*)spin2, 1);

   host_rk = host_dotprod;
 
  // recalculate residue frome r = b - Ax
   if(((i+0) % N_recalcres) == 0)
   {
    // r_(k+1) = Ax -b 
     printf("Recalculating residue\n");
    
     mixedsolveOperator.gpu(spin1,spin4,spin3,gf,dev_nn,griddim2,blockdim2);// call func
    // r = b - Ax
     cublasScal (24*VOLUME, -1.0, (RealT*)spin3, 1);
     cublasAxpy (24*VOLUME, 1.0, (const RealT*)spinin,  1, (RealT*)spin3, 1);
     cublasCopy (24*VOLUME, (const RealT*)spin3, 1, (RealT*)spin0, 1);

     if( sizeof(RealT)>=sizeof(REALD) && ((host_rk<=eps&&rel_prec==0) || (host_rk<=finalEps*initial_sourcesquarenorm&&rel_prec==1)) )//different from abort criterium some lines above: here we check wether we reached the final desired precision, which only works in double precision 
     {//the final precision is reached
       printf("inner solver: Reached precision of eps=%.2e\n",( rel_prec==0 ? eps : finalEps ));
       break;//escape innner solver if desired prec. is reached: should not happen with singele precision - here only the double prec. outer solver is reliable
     }
   }//recalculate residue
 }//MAIN LOOP cg
  
  
  printf("Final residue: %.6e\n",host_dotprod);
  // x_result = spin1 !
  
  mixedsolveOperator.gpuDeinit(spin1,spin4,gf,dev_nn,griddim2,blockdim2,scaleparam);
  dev_copy_spinor_field<RealT,RealT> <<<griddim2, blockdim2 >>>(spin1,spinout);
  
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



template<class RealT>
void showspinor(dev_spinorM(RealT)* s){
  int i,j;
  dev_spinor help[6];
  size_t size = 6*sizeof(dev_spinorM(RealT));

  for(i=0; i<VOLUME/2; i++){
    cudaMemcpy(&(help[0]), (s+6*i), size, cudaMemcpyDeviceToHost);
    for(j=0;j<6; j++){
      printf("(%.3f %.3f) (%.3f, %.3f) ", help[j].x, help[j].y, help[j].z, help[j].w);
    }
    printf("\n");
  }
}




#ifndef HALF

// this is the eo version of the device cg inner solver 
// we invert the hermitean Q_{-} Q_{+}
template<class RealT>
int dev_cg_eo(
      dev_su3_2vM(RealT)* gf,
      dev_spinorM(RealT)* spinin, 
      dev_spinorM(RealT)* spinout, 
      dev_spinorM(RealT)* spin0, 
      dev_spinorM(RealT)* spin1, 
      dev_spinorM(RealT)* spin2, 
      dev_spinorM(RealT)* spin3, 
      dev_spinorM(RealT)* spin4, 
      int* grid, int* nn_grid, RealT epsfinal, MixedsolveParameter<RealT>& mixedsolveParameter){


 RealT host_alpha, host_beta, host_dotprod, host_rk, sourcesquarenorm;
 RealT * dotprod, * dotprod2, * rk, * alpha, *beta;
 
 
 
 int i, gridsize;
 int maxit = max_innersolver_it;
 RealT eps = (RealT) innersolver_precision;
 int N_recalcres = 40; // after N_recalcres iterations calculate r = A x_k - b
 
 cudaError_t cudaerr;
 
 // this is the partitioning for the copying of fields
 dim3 blockdim(1,1);
 //dim3 blockdim2(128,1,1);
 
 int blockdim2 = BLOCK3;
 if( VOLUME/2 % blockdim2 == 0){
   gridsize = (int) VOLUME/2/blockdim2;
 }
 else{
   gridsize = (int) VOLUME/2/blockdim2 + 1;
 }
 int griddim2 = gridsize;

 
 //this is the partitioning for the HoppingMatrix kernel
 /*
 int blockdim3=BLOCK;
 if( VOLUME/2 >= BLOCK){
   gridsize = (int)(VOLUME/2/BLOCK) + 1;
 }
 else{
   gridsize=1;
 }
 int griddim3=gridsize; 
 */
 int blockdim3 = BLOCK;
 if( VOLUME/2 % blockdim3 == 0){
   gridsize = (int) VOLUME/2/blockdim3;
 }
 else{
   gridsize = (int) VOLUME/2/blockdim3 + 1;
 }
 int griddim3 = gridsize;
   

    if (g_proc_id == 0) { printf("gridsize = %d\nsizeof(Real) = %hi\n", gridsize, sizeof(RealT)); }


 
 //this is the partitioning for dev_mul_one_pm...
 /*
 int blockdim4=BLOCK2;
 if( VOLUME/2 >= BLOCK2){
   gridsize = (int)(VOLUME/2/BLOCK2) + 1;
 }
 else{
   gridsize=1;
 }
 int griddim4=gridsize;  
 */
 int blockdim4 = BLOCK2;
 if( VOLUME/2 % blockdim4 == 0){
   gridsize = (int) VOLUME/2/blockdim4;
 }
 else{
   gridsize = (int) VOLUME/2/blockdim4 + 1;
 }
 int griddim4 = gridsize;
 
 
 //never referenced: size_t size2 = sizeof(dev_spinorM(RealT))*6*VOLUME/2;
 
 
 //Initialize some stuff
 

    if (g_proc_id == 0) printf("mu = %f\n", g_mu);

  
  
  
  dev_complex h0,h1,h2,h3,mh0, mh1, mh2, mh3;
  h0.re = (RealT)ka0.re;    h0.im = -(RealT)ka0.im;
  h1.re = (RealT)ka1.re;    h1.im = -(RealT)ka1.im;
  h2.re = (RealT)ka2.re;    h2.im = -(RealT)ka2.im;
  h3.re = (RealT)ka3.re;    h3.im = -(RealT)ka3.im;
  
  mh0.re = -(RealT)ka0.re;    mh0.im = (RealT)ka0.im;
  mh1.re = -(RealT)ka1.re;    mh1.im = (RealT)ka1.im;
  mh2.re = -(RealT)ka2.re;    mh2.im = (RealT)ka2.im;
  mh3.re = -(RealT)ka3.re;    mh3.im = (RealT)ka3.im;
  
  // try using constant mem for kappas
  cudaMemcpyToSymbol("dev_k0c", &h0, sizeof(h0)) ; 
  cudaMemcpyToSymbol("dev_k1c", &h1, sizeof(h1)) ; 
  cudaMemcpyToSymbol("dev_k2c", &h2, sizeof(h2)) ; 
  cudaMemcpyToSymbol("dev_k3c", &h3, sizeof(h3)) ;
  
  cudaMemcpyToSymbol("dev_mk0c", &mh0, sizeof(mh0)) ; 
  cudaMemcpyToSymbol("dev_mk1c", &mh1, sizeof(mh1)) ; 
  cudaMemcpyToSymbol("dev_mk2c", &mh2, sizeof(mh2)) ; 
  cudaMemcpyToSymbol("dev_mk3c", &mh3, sizeof(mh3)) ;  
  
  he_cg_init<<< 1, 1 >>> (grid, (REAL) g_kappa, (REAL)(g_mu/(2.0*g_kappa)), h0,h1,h2,h3);
  // BEWARE in dev_tm_dirac_kappa we need the true mu (not 2 kappa mu!)
  
  #ifdef MPI
    he_cg_init_nd_additional_mpi<<<1,1>>>(VOLUMEPLUSRAND, RAND, g_cart_id, g_nproc);
    // debug	// check dev_VOLUMEPLUSRAND and dev_RAND on device
  	if (g_proc_id == 0) {
  	  int host_check_VOLUMEPLUSRAND, host_check_RAND;
  	  int host_check_rank, host_check_nproc;
  	  cudaMemcpyFromSymbol(&host_check_VOLUMEPLUSRAND, dev_VOLUMEPLUSRAND, sizeof(int));
  	  cudaMemcpyFromSymbol(&host_check_RAND, dev_RAND, sizeof(int));
  	  printf("\tOn device:\n");
  	  printf("\tdev_VOLUMEPLUSRAND = %i\n", host_check_VOLUMEPLUSRAND);
  	  printf("\tdev_RAND = %i\n", host_check_RAND);
  	  cudaMemcpyFromSymbol(&host_check_rank, dev_rank, sizeof(int));
  	  cudaMemcpyFromSymbol(&host_check_nproc, dev_nproc, sizeof(int));
  	  printf("\tdev_rank = %i\n", host_check_rank);
  	  printf("\tdev_nproc = %i\n", host_check_nproc);
  	}
  #endif
  
  
  #ifdef USETEXTURE
    //Bind texture gf
    bind_texture_gf(gf);
  #endif
 
 
 // Init x,p,r for k=0
 // Allocate some numbers for host <-> device interaction
 cudaMalloc((void **) &dotprod, sizeof(RealT));
 cudaMalloc((void **) &dotprod2, sizeof(RealT));
 cudaMalloc((void **) &rk, sizeof(RealT));
 cudaMalloc((void **) &alpha, sizeof(RealT));
 cudaMalloc((void **) &beta, sizeof(RealT));
 printf("%s\n", cudaGetErrorString(cudaGetLastError())); 
 
 
 //init blas
 cublasInit();
  

    if (g_proc_id == 0) {
      printf("%s\n", cudaGetErrorString(cudaGetLastError())); 
      printf("have initialized cublas\n"); 
    }

 

 
 

 //dev_skalarmult_spinor_field<<<griddim2, blockdim2 >>>(spinin,scaleparam, spin1);
 //dev_copy_spinor_field<<<griddim2, blockdim2 >>>(spin1, spinin);
 
 
 dev_copy_spinor_field<RealT,RealT> <<<griddim2, blockdim2 >>>(spinin, spin0);
 dev_zero_spinor_field<RealT> <<<griddim2, blockdim2 >>>(spin1); // x_0 = 0
 dev_copy_spinor_field<RealT,RealT> <<<griddim2, blockdim2 >>>(spinin, spin2);
 dev_zero_spinor_field<RealT> <<<griddim2, blockdim2 >>>(spin3);
  

    if (g_proc_id == 0) printf("%s\n", cudaGetErrorString(cudaGetLastError()));

 
 
 

 //relative precision -> get initial residue
 #ifndef MPI
   sourcesquarenorm = cublasDot (24*VOLUME/2, (const RealT*)spinin, 1, (const RealT*)spinin, 1);
 #else
   sourcesquarenorm = cublasDot_wrapper (24*VOLUME/2, (RealT*)spinin, 1, (RealT*)spinin, 1);
 #endif
 host_rk = sourcesquarenorm; //for use in main loop
 

    if (g_proc_id == 0) {
      printf("Squarenorm Source:\t%.16e\n", sourcesquarenorm);
      printf("%s\n", cudaGetErrorString(cudaGetLastError()));
      printf("Entering inner solver cg-loop\n");
    }

 
 
 
 
 
 for(i=0;i<maxit;i++){ //MAIN LOOP
  
  // Q_{-}Q{+}
  #ifndef MPI
    dev_Qtm_pm_psi    <RealT>(spin2, spin3, griddim3, blockdim3, griddim4, blockdim4, mixedsolveParameter);
  #else
    dev_Qtm_pm_psi_mpi<RealT>(spin2, spin3, griddim3, blockdim3, griddim4, blockdim4, mixedsolveParameter);
  #endif
  
  if((cudaerr=cudaGetLastError()) != cudaSuccess){
    printf("%s\n", cudaGetErrorString(cudaerr));
    exit(200);
  }
  
  
 //alpha
  #ifndef MPI
    host_dotprod = cublasDot (24*VOLUME/2, (const RealT*) spin2, 1, (const RealT*) spin3, 1);
  #else
    host_dotprod = cublasDot_wrapper (24*VOLUME/2, (RealT*) spin2, 1, (RealT*) spin3, 1);
  #endif
  
  host_alpha = (host_rk / host_dotprod); // alpha = r*r/ p M p
   
 //r(k+1)
 cublasAxpy (24*VOLUME/2,-1.0*host_alpha, (const RealT*)spin3, 1, (RealT*)spin0, 1);  


 //x(k+1);
 cublasAxpy (24*VOLUME/2, host_alpha, (const RealT*)spin2,  1, (RealT*)spin1, 1);
 
  if((cudaerr=cudaGetLastError()) != cudaSuccess){
    printf("%s\n", cudaGetErrorString(cudaerr));
    exit(200);
  }

  //Abbruch?
  #ifndef MPI
    host_dotprod = cublasDot (24*VOLUME/2, (const RealT*) spin0, 1,(const RealT*) spin0, 1);
  #else
    host_dotprod = cublasDot_wrapper (24*VOLUME/2, (RealT*) spin0, 1,(RealT*) spin0, 1);
  #endif
  
 if (((host_dotprod <= eps*sourcesquarenorm) && (i > maxit / 4) ) || ( host_dotprod <= epsfinal/2.)){//error-limit erreicht (epsfinal/2 sollte ausreichen um auch in double precision zu bestehen)
   break; 
 }
  
  
    if (g_proc_id == 0) printf("iter %d: err = %.16e\n", i, host_dotprod);

  
  
 //beta
 host_beta =host_dotprod/host_rk;
 //p(k+1)
 cublasScal (24*VOLUME/2, host_beta, (RealT*)spin2, 1);
 cublasAxpy (24*VOLUME/2, 1.0, (const RealT*)spin0,  1, (RealT*)spin2, 1);

 host_rk = host_dotprod;
 
 // recalculate residue frome r = b - Ax
 if(((i+1) % N_recalcres) == 0){
    // r_(k+1) = Ax -b 

    if (g_proc_id == 0) printf("Recalculating residue\n");

    // D Ddagger   --   Ddagger = gamma5 D gamma5  for Wilson Dirac Operator
    // DO NOT USE tm_dirac_dagger_kappa here, otherwise spin2 will be overwritten!!!

    // Q_{-}Q{+}
    #ifndef MPI
        dev_Qtm_pm_psi    <RealT>(spin1, spin3, griddim3, blockdim3, griddim4, blockdim4, mixedsolveParameter);
    #else
        dev_Qtm_pm_psi_mpi<RealT>(spin1, spin3, griddim3, blockdim3, griddim4, blockdim4, mixedsolveParameter);
    #endif
    if((cudaerr=cudaGetLastError()) != cudaSuccess){
      printf("%s\n", cudaGetErrorString(cudaerr));
      exit(200);
    }  
        
    
    // r = b - Ax
    cublasScal (24*VOLUME/2, -1.0, (RealT*)spin3, 1);
    cublasAxpy (24*VOLUME/2, 1.0, (const RealT*)spinin,  1, (RealT*)spin3, 1);
    cublasCopy (24*VOLUME/2, (const RealT*)spin3, 1, (RealT*)spin0, 1);
    //dev_skalarmult_add_assign_spinor_field<<<griddim2, blockdim2 >>>(spinin, -1.0, spin3, spin0);
   }//recalculate residue

 }//MAIN LOOP cg	


    if (g_proc_id == 0) printf("Final residue: %.16e\n",host_dotprod);


  // x_result = spin1 !

  //no multiplication with D^{dagger} here and no return to non-kappa basis as in dev_cg!
  dev_copy_spinor_field<RealT,RealT> <<<griddim2, blockdim2 >>>(spin1,spinout);

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








// convert spinor to double
template<class RealT>
void convert2double_spin (typename dev_spinorT<RealT>::type* spin, spinor* h2d) {

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
template<class RealT>
void convert2REAL4_spin(spinor* spin, typename dev_spinorT<RealT>::type* h2d){

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

        h2d[6*i+0].x = (RealT) spin[i].s0.c0.re;
        h2d[6*i+0].y = (RealT) spin[i].s0.c0.im;
        h2d[6*i+0].z = (RealT) spin[i].s0.c1.re;
        h2d[6*i+0].w = (RealT) spin[i].s0.c1.im;

        h2d[6*i+1].x = (RealT) spin[i].s0.c2.re;
        h2d[6*i+1].y = (RealT) spin[i].s0.c2.im;
        h2d[6*i+1].z = (RealT) spin[i].s1.c0.re;
        h2d[6*i+1].w = (RealT) spin[i].s1.c0.im;

        h2d[6*i+2].x = (RealT) spin[i].s1.c1.re;
        h2d[6*i+2].y = (RealT) spin[i].s1.c1.im;
        h2d[6*i+2].z = (RealT) spin[i].s1.c2.re;
        h2d[6*i+2].w = (RealT) spin[i].s1.c2.im;

        h2d[6*i+3].x = (RealT) spin[i].s2.c0.re;
        h2d[6*i+3].y = (RealT) spin[i].s2.c0.im;
        h2d[6*i+3].z = (RealT) spin[i].s2.c1.re;
        h2d[6*i+3].w = (RealT) spin[i].s2.c1.im;

        h2d[6*i+4].x = (RealT) spin[i].s2.c2.re;
        h2d[6*i+4].y = (RealT) spin[i].s2.c2.im;
        h2d[6*i+4].z = (RealT) spin[i].s3.c0.re;
        h2d[6*i+4].w = (RealT) spin[i].s3.c0.im;

        h2d[6*i+5].x = (RealT) spin[i].s3.c1.re;
        h2d[6*i+5].y = (RealT) spin[i].s3.c1.im;
        h2d[6*i+5].z = (RealT) spin[i].s3.c2.re;
        h2d[6*i+5].w = (RealT) spin[i].s3.c2.im;

  }
}










template<class RealT>
MixedsolveParameter<RealT>* init_mixedsolve(su3** gf){
  
   cudaError_t cudaerr;
   MixedsolveParameter<RealT>& mixedsolveParameter=*MixedsolveParameter<RealT>::getGlobalP();

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
	      printf("Error in init_mixedsolve_eo(): Could not set active device. Aborting...\n");
	      exit(302);
	    }
           }
           else{
            printf("Not setting any active device. Let the driver choose.\n");
            int device=-1;cudaGetDevice(&device);printf("device=%i",device);
           }        
    havedevice = 1;
    }
  #ifdef GF_8
  /* allocate 8 floats of gf = 2*4*VOLUME float4's*/
  printf("Using GF 8 reconstruction\n");
  size_t dev_gfsize = 2*4*VOLUME * sizeof(dev_su3_8M(RealT));
  #else
  /* allocate 2 rows of gf = 3*4*VOLUME float4's*/
  printf("Using GF 12 reconstruction\n");
  size_t dev_gfsize = 3*4*VOLUME * sizeof(dev_su3_2vM(RealT)); 
  #endif
  
  #ifdef USETEXTURE
    printf("Using texture references\n");
  #else
    printf("NOT using texture references\n");
  #endif
  if((cudaerr=cudaMalloc((void **) &mixedsolveParameter.dev_gf, dev_gfsize)) != cudaSuccess){
    printf("Error in init_mixedsolve(): Memory allocation of gauge field failed. Aborting...\n");
    exit(200);
  }   // Allocate array on device
  else{
    printf("Allocated gauge field on device\n");
  }  
  
  #ifdef GF_8
  mixedsolveParameter.h2d_gf = (dev_su3_8M(RealT)*)malloc(dev_gfsize); // Allocate REAL conversion gf on host
  su3to8<RealT>(gf,mixedsolveParameter.h2d_gf);  
  #else
  mixedsolveParameter.h2d_gf = (dev_su3_2vM(RealT)*)malloc(dev_gfsize); // Allocate REAL conversion gf on host
  su3to2vf4<RealT>(gf,mixedsolveParameter.h2d_gf);
  #endif
  cudaMemcpy(mixedsolveParameter.dev_gf, mixedsolveParameter.h2d_gf, dev_gfsize, cudaMemcpyHostToDevice);


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
  size_t dev_spinsize = 6*VOLUME * sizeof(dev_spinorM(RealT)); /* float4 */  
  if((void*)(mixedsolveParameter.h2d_spin = (dev_spinorM(RealT)*)malloc(dev_spinsize)) == NULL){
    printf("Could not allocate memory for mixedsolveParameter.h2d_spin. Aborting...\n");
    exit(200);
  } // Allocate float conversion spinor on host
  #else
  size_t dev_spinsize = 6*VOLUME * sizeof(dev_spinor_half); /*short4*/  
  if((void*)(mixedsolveParameter.h2d_spin = (dev_spinor_half *)malloc(dev_spinsize)) == NULL){
    printf("Could not allocate memory for mixedsolveParameter.h2d_spin. Aborting...\n");
    exit(200);
  } // Allocate float conversion spinor on host 
  size_t dev_normsize = VOLUME/2 * sizeof(float);
  if((void*)(mixedsolveParameter.h2d_spin_norm = (float*)malloc(dev_normsize)) == NULL){
    printf("Could not allocate memory for mixedsolveParameter.h2d_spin_norm. Aborting...\n");
    exit(200);
  } // Allocate float conversion norm on host 
  #endif
  
  
  cudaMalloc((void **) &mixedsolveParameter.dev_spin1, dev_spinsize);   // Allocate array spin1 on device
  cudaMalloc((void **) &mixedsolveParameter.dev_spin2, dev_spinsize);   // Allocate array spin2 on device
  cudaMalloc((void **) &mixedsolveParameter.dev_spin3, dev_spinsize);   // Allocate array spin3 on device
  cudaMalloc((void **) &mixedsolveParameter.dev_spin4, dev_spinsize);
  cudaMalloc((void **) &mixedsolveParameter.dev_spin5, dev_spinsize);
  cudaMalloc((void **) &mixedsolveParameter.dev_spinin, dev_spinsize);
  cudaMalloc((void **) &mixedsolveParameter.dev_spinout, dev_spinsize);

  #ifdef HALF
   dev_spinsize = VOLUME/2*sizeof(float);
   cudaMalloc((void **) &mixedsolveParameter.dev_spin1_norm, dev_spinsize);   // Allocate norm spin1 on device
   cudaMalloc((void **) &mixedsolveParameter.dev_spin2_norm, dev_spinsize);   // Allocate norm spin2 on device
   cudaMalloc((void **) &mixedsolveParameter.dev_spin3_norm, dev_spinsize);   // Allocate norm spin3 on device
   cudaMalloc((void **) &mixedsolveParameter.dev_spin4_norm, dev_spinsize);
   cudaMalloc((void **) &mixedsolveParameter.dev_spin5_norm, dev_spinsize);
   cudaMalloc((void **) &mixedsolveParameter.dev_spinin_norm, dev_spinsize);
   cudaMalloc((void **) &mixedsolveParameter.dev_spinout_norm, dev_spinsize);
  #endif


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


  return &mixedsolveParameter;  
}





template<class RealT>
MixedsolveParameter<RealT>* init_mixedsolve_eo(su3** gf){

  cudaError_t                 cudaerr;
  MixedsolveParameter<RealT>& mixedsolveParameter=*MixedsolveParameter<RealT>::getGlobalP();

  if (havedevice == 0) {
  
    // get number of devices
    int ndev = find_devices();
    if(ndev == 0){
      fprintf(stderr, "Error: no CUDA devices found. Aborting...\n");
      exit(300);
    }
    
    // try to set active device to device_num given in input file (or mpi rank)
    #ifndef MPI
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
    	      printf("Process %d of %d: Setting active device to: %d\n", g_proc_id, g_nproc, device_num);
    	      cudaSetDevice(device_num);
    	    }
    	    else {
    	      fprintf(stderr, "Process %d of %d: Error: There is no CUDA device with No. %d. Aborting...\n", g_proc_id, g_nproc, device_num);
    	      exit(301);
    	    }
          }
          else{
            printf("Not setting any active device. Let the driver choose.\n");
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
    #endif
    
    havedevice=1;
    
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
  
  #ifndef MPI
  	#ifdef GF_8
  	  /* allocate 8 floats for gf = 2*4*VOLUME float4's*/
  	  size_t dev_gfsize = 2*4*VOLUME * sizeof(dev_su3_8M(RealT));
  	#else
  	  /* allocate 2 rows of gf = 3*4*VOLUME float4's*/
  	  size_t dev_gfsize = 3*4*VOLUME * sizeof(dev_su3_2vM(RealT));
  	#endif
  #else
  	#ifdef GF_8
  	  /* allocate 8 floats for gf = 2*4*VOLUME float4's*/
  	  size_t dev_gfsize = 2*4*(VOLUME+RAND) * sizeof(dev_su3_8M(RealT));
  	#else
  	  /* allocate 2 rows of gf = 3*4*VOLUME float4's*/
  	  size_t dev_gfsize = 3*4*(VOLUME+RAND) * sizeof(dev_su3_2vM(RealT));
  	#endif
  
  #endif
  
  if((cudaerr=cudaMalloc((void **) &mixedsolveParameter.dev_gf, dev_gfsize)) != cudaSuccess){
    printf("Error in init_mixedsolve(): Memory allocation of gauge field failed. Aborting...\n");
    exit(200);
  }   // Allocate array on device
  else {
    #ifndef MPI
      printf("Allocated memory for gauge field on device.\n");
    #else
      if (g_cart_id == 0) printf("Allocated memory for gauge field on devices.\n");
    #endif
  }
  
  #ifdef GF_8
    mixedsolveParameter.h2d_gf = (dev_su3_8M(RealT)*)malloc(dev_gfsize); // Allocate REAL conversion gf on host
    su3to8<RealT>(gf,mixedsolveParameter.h2d_gf);
  #else
    mixedsolveParameter.h2d_gf = (dev_su3_2vM(RealT)*)malloc(dev_gfsize); // Allocate REAL conversion gf on host
    su3to2vf4<RealT>(gf,mixedsolveParameter.h2d_gf);
  #endif
  //bring to device
  cudaMemcpy(mixedsolveParameter.dev_gf, mixedsolveParameter.h2d_gf, dev_gfsize, cudaMemcpyHostToDevice);
  
  
  #ifdef HALF
    #ifndef MPI
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
    if((cudaerr=cudaMalloc((void **) &mixedsolveParameter.dev_gf_half, dev_gfsize)) != cudaSuccess){
    printf("Error in init_mixedsolve(): Memory allocation of half precsion gauge field failed. Aborting...\n");
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
  
  #ifndef MPI
    size_t idxsize = VOLUME/2*sizeof(int);
  #else
    size_t idxsize = (VOLUME+RAND)/2*sizeof(int);
  #endif
  eoidx_even = (int *) malloc(idxsize);
  eoidx_odd = (int *) malloc(idxsize);
  cudaMalloc((void **) &dev_eoidx_even, idxsize);
  cudaMalloc((void **) &dev_eoidx_odd, idxsize);
  
  #ifndef MPI
    initnn();
    initnn_eo();
    //shownn_eo();
  #else
    init_nnspinor_eo_mpi();
    init_idxgauge_mpi();
  #endif
  
  //shownn();
  //showcompare_gf(T-1, LX-1, LY-1, LZ-1, 3);
  //check_gauge_reconstruction_8(gf, mixedsolveParameter.dev_gf, 0, 0);
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
  #ifndef HALF
  	size_t dev_spinsize = 6*VOLUME/2 * sizeof(dev_spinorM(RealT)); /* float4 */
  	if((void*)(mixedsolveParameter.h2d_spin = (dev_spinorM(RealT)*)malloc(dev_spinsize)) == NULL){
  	  printf("Could not allocate memory for mixedsolveParameter.h2d_spin. Aborting...\n");
  	  exit(200);
  	} // Allocate float conversion spinor on host
  	#ifdef MPI
  	  size_t dev_spinsize_ext =  6*(VOLUME+RAND)/2*sizeof(dev_spinorM(RealT));
  	#endif
  #else
  	size_t dev_spinsize = 6*VOLUME/2 * sizeof(dev_spinor_half);/*short4*/
  	if((void*)(mixedsolveParameter.h2d_spin = (dev_spinor_half *)malloc(dev_spinsize)) == NULL){
  	  printf("Could not allocate memory for mixedsolveParameter.h2d_spin. Aborting...\n");
  	  exit(200);
  	} // Allocate float conversion spinor on host 
  	size_t dev_normsize = VOLUME/2 * sizeof(RealT);
  	if((void*)(mixedsolveParameter.h2d_spin_norm = (RealT *)malloc(dev_normsize)) == NULL){
  	  printf("Could not allocate memory for mixedsolveParameter.h2d_spin_norm. Aborting...\n");
  	  exit(200);
  	} // Allocate float conversion norm on host 
  	#ifdef MPI
  	  size_t dev_spinsize_ext =  6*(VOLUME+RAND)/2*sizeof(dev_spinor_half);
  	  size_t dev_normsize_ext =  (VOLUME+RAND)/2*sizeof(float);
  	#endif
  #endif
  
  
  #ifndef MPI
  	cudaMalloc((void **) &mixedsolveParameter.dev_spin1, dev_spinsize);   // Allocate array spin1 on device
  	cudaMalloc((void **) &mixedsolveParameter.dev_spin2, dev_spinsize);   // Allocate array spin2 on device
  	cudaMalloc((void **) &mixedsolveParameter.dev_spin3, dev_spinsize);   // Allocate array spin3 on device
  	cudaMalloc((void **) &mixedsolveParameter.dev_spin4, dev_spinsize);
  	cudaMalloc((void **) &mixedsolveParameter.dev_spin5, dev_spinsize);
  	cudaMalloc((void **) &mixedsolveParameter.dev_spinin, dev_spinsize);
  	cudaMalloc((void **) &mixedsolveParameter.dev_spinout, dev_spinsize);
  	
  	cudaMalloc((void **) &mixedsolveParameter.dev_spin_eo1, dev_spinsize);
  	cudaMalloc((void **) &mixedsolveParameter.dev_spin_eo2, dev_spinsize);
 
 
       #ifdef HALF
         cudaMalloc((void **) &mixedsolveParameter.dev_spin1_norm, dev_spinsize);   // Allocate norm spin1 on device
         cudaMalloc((void **) &mixedsolveParameter.dev_spin2_norm, dev_spinsize);   // Allocate norm spin2 on device
         cudaMalloc((void **) &mixedsolveParameter.dev_spin3_norm, dev_spinsize);   // Allocate norm spin3 on device
         cudaMalloc((void **) &mixedsolveParameter.dev_spin4_norm, dev_spinsize);
         cudaMalloc((void **) &mixedsolveParameter.dev_spin5_norm, dev_spinsize);
         cudaMalloc((void **) &mixedsolveParameter.dev_spinin_norm, dev_spinsize);
         cudaMalloc((void **) &mixedsolveParameter.dev_spinout_norm, dev_spinsize);

        cudaMalloc((void **) &mixedsolveParameter.dev_spin_eo1_norm, dev_spinsize);
        cudaMalloc((void **) &mixedsolveParameter.dev_spin_eo2_norm, dev_spinsize);
      #endif  
  
  
  #else
  	cudaMalloc((void **) &mixedsolveParameter.dev_spin1, dev_spinsize_ext);
  	cudaMalloc((void **) &mixedsolveParameter.dev_spin2, dev_spinsize_ext);
  	cudaMalloc((void **) &mixedsolveParameter.dev_spin3, dev_spinsize_ext);
  	cudaMalloc((void **) &mixedsolveParameter.dev_spin4, dev_spinsize_ext);
  	cudaMalloc((void **) &mixedsolveParameter.dev_spin5, dev_spinsize_ext);
  	cudaMalloc((void **) &mixedsolveParameter.dev_spinin, dev_spinsize_ext);
  	cudaMalloc((void **) &mixedsolveParameter.dev_spinout, dev_spinsize_ext);
  	
  	cudaMalloc((void **) &mixedsolveParameter.dev_spin_eo1, dev_spinsize_ext);
  	cudaMalloc((void **) &mixedsolveParameter.dev_spin_eo2, dev_spinsize_ext);
  	
        #ifdef HALF
         cudaMalloc((void **) &mixedsolveParameter.dev_spin1_norm, dev_normsize_ext);   // Allocate norm spin1 on device
         cudaMalloc((void **) &mixedsolveParameter.dev_spin2_norm, dev_normsize_ext);   // Allocate norm spin2 on device
         cudaMalloc((void **) &mixedsolveParameter.dev_spin3_norm, dev_normsize_ext);   // Allocate norm spin3 on device
         cudaMalloc((void **) &mixedsolveParameter.dev_spin4_norm, dev_normsize_ext);
         cudaMalloc((void **) &mixedsolveParameter.dev_spin5_norm, dev_normsize_ext);
         cudaMalloc((void **) &mixedsolveParameter.dev_spinin_norm, dev_normsize_ext);
         cudaMalloc((void **) &mixedsolveParameter.dev_spinout_norm, dev_normsize_ext);

        cudaMalloc((void **) &mixedsolveParameter.dev_spin_eo1_norm, dev_normsize_ext);
        cudaMalloc((void **) &mixedsolveParameter.dev_spin_eo2_norm, dev_normsize_ext);
      #endif   	
      
      int tSliceEO = LX*LY*LZ/2;
      #ifndef HALF
  	R1 = (dev_spinor *) malloc(2*tSliceEO*24*sizeof(float));
  	R2 = R1 + 6*tSliceEO;
  	R3 = (dev_spinor *) malloc(2*tSliceEO*24*sizeof(float));
  	R4 = R3 + 6*tSliceEO;
      #else
      
  	// implement this for half?
  	// -> ALTERNATE_FIELD_EXCHANGE     
      #endif
  	
  #endif
 
 



  if((cudaerr=cudaGetLastError())!=cudaSuccess){
    printf("Error in init_mixedsolve(): Memory allocation of spinor fields failed. Aborting...\n");
    exit(200);
  }
  else{
    printf("Allocated spinor fields on device\n");
  }
  

  #ifdef MPI
    /*  for async communication */
    // page-locked memory
   #ifndef HALF 
    cudaMallocHost(&RAND3, 2*tSliceEO*6*sizeof(REAL4M(RealT)));
    RAND4 = RAND3 + 6*tSliceEO;
    cudaMallocHost(&RAND1, 2*tSliceEO*6*sizeof(REAL4M(RealT)));
    RAND2 = RAND1 + 6*tSliceEO;
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
          
    // CUDA streams and events
    for (int i = 0; i < 3; i++) {
        cudaStreamCreate(&stream[i]);
    }    
    /* end for async communication */
  #endif
  
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
    devplaq = calc_plaquette(mixedsolveParameter.dev_gf, dev_nn);
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
    devrect = calc_rectangle(mixedsolveParameter.dev_gf, dev_nn);
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

  calc_polyakov_0(&ret, mixedsolveParameter.dev_gf, dev_nn);
  printf("Calculating Polyakov loop on device:\n");  
  printf("pl_0 (Re) = %.8e\n",ret.x);
  printf("pl_0 (Im) = %.8e\n",ret.y);
  
  //polyakov_loop_dir(1, 0);
  //printf("Calculating Polyakov loop on host:\n");  
 
  finalize_dev_observables();

  exit(100);
  */


  return &mixedsolveParameter;  
}



template<class RealT>
void finalize_mixedsolve(MixedsolveParameter<RealT>* mixedsolveParameterP){

  MixedsolveParameter<RealT>& mixedsolveParameter=*mixedsolveParameterP;//use pointer in interface so we can delete mix\.solv\.Param later here

  cudaFree(mixedsolveParameter.dev_spin1);
  cudaFree(mixedsolveParameter.dev_spin2);
  cudaFree(mixedsolveParameter.dev_spin3);
  cudaFree(mixedsolveParameter.dev_spin4);
  cudaFree(mixedsolveParameter.dev_spin5);
  cudaFree(mixedsolveParameter.dev_spinin);
  cudaFree(mixedsolveParameter.dev_spinout);
  cudaFree(mixedsolveParameter.dev_gf);
  cudaFree(dev_grid);
  cudaFree(dev_output);
  cudaFree(dev_nn);
  
  if(even_odd_flag){
    cudaFree(mixedsolveParameter.dev_spin_eo1);
    cudaFree(mixedsolveParameter.dev_spin_eo2);
    cudaFree(dev_eoidx_even);
    cudaFree(dev_eoidx_odd);
    cudaFree(dev_nn_eo);
    cudaFree(dev_nn_oe);
  }
  
  #ifdef HALF
    cudaFree(mixedsolveParameter.dev_gf_half);
 
    cudaFree(mixedsolveParameter.dev_spin1_norm);
    cudaFree(mixedsolveParameter.dev_spin2_norm);
    cudaFree(mixedsolveParameter.dev_spin3_norm);
    cudaFree(mixedsolveParameter.dev_spin4_norm);
    cudaFree(mixedsolveParameter.dev_spin5_norm);
    cudaFree(mixedsolveParameter.dev_spinin_norm);
    cudaFree(mixedsolveParameter.dev_spinout_norm);
    
    if(even_odd_flag){
     cudaFree(mixedsolveParameter.dev_spin_eo1_norm);
     cudaFree(mixedsolveParameter.dev_spin_eo2_norm);
    }
    
    
  #endif
  
#ifdef MPI
  cudaFreeHost(RAND1);
  cudaFreeHost(RAND3);
  
  #ifdef HALF
   cudaFreeHost(RAND1_norm);
   cudaFreeHost(RAND3_norm);
  #endif
             
  for (int i = 0; i < 3; i++) {
     cudaStreamDestroy(stream[i]);
  } 
#endif 
  
  
  
  free(mixedsolveParameter.h2d_spin);
  free(mixedsolveParameter.h2d_gf);
}


// include half versions of dev_cg - solvers
#ifdef HALF
  #include "half_solvers.cuh"
#endif





#ifndef HALF
template<class RealT,template<class MixedsolveOperatorRealT>class MixedsolveOperatorT>
int mixed_solveT(spinor * const P, spinor * const Q, const int max_iter, 
		 double eps, const int rel_prec,const int N, MixedsolveOperatorT<RealT>& mixedsolveOperator){
 
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

  size_t dev_spinsize = 6*VOLUME * sizeof(dev_spinorM(RealT)); // float4 
  MixedsolveParameter<RealT>& mixedsolveParameter=*init_mixedsolve<RealT>(g_gauge_field);
  
  // Start timer
  assert((start = clock())!=-1);

  rk = square_norm(Q, N, 0);
  sourcesquarenorm = rk; // for relative precision
  assign(solver_field[0],Q,N);
  printf("Initial residue: %.16e\n",rk);
  zero_spinor_field(solver_field[1],  N);//spin2 = x_k
  zero_spinor_field(solver_field[2],  N);
  printf("The VOLUME is: %d\n",N);

  mixedsolveOperator.checkInit(solver_field[2],solver_field[3],solver_field[0],N);
  
  
  for(iter=0; iter<max_iter; iter++){

    //"Applying double precision Dirac-Op...\n"
    mixedsolveOperator.check(solver_field[2],solver_field[3],N);//spinTmp ^= solver_field[3], N=volume
    // r_k = b - D x_k
    diff(solver_field[0], solver_field[0], solver_field[3] ,N);//residueRSpininout ^= solver_field[0]

    rk = square_norm(solver_field[0], N, 0);

#ifdef GF_8
    if(isnan(rk)){
      fprintf(stderr, "Error in mixed_solveT: Residue is NaN.\n  May happen with GF 8 reconstruction. Aborting ...\n");
      exit(200);
    }
#endif

    printf("Residue after %d inner solver iterations: %.18e\n",outercount,rk);
    if((rk<=eps && rel_prec==0) || (rk<=eps*sourcesquarenorm && rel_prec==1))
      {
	printf("Reached solver precision of eps=%.2e\n",eps);
	//multiply with D^dagger
	mixedsolveOperator.checkDeinit(solver_field[1],solver_field[3],P,N);


	stop = clock();
	timeelapsed = (double) (stop-start)/CLOCKS_PER_SEC;
	printf("Inversion done in mixed precision.\n Number of iterations in outer solver: %d\n Squared residue: %.16e\n Time elapsed: %.6e sec\n", outercount, rk, timeelapsed);
	finalize_mixedsolve(&mixedsolveParameter);
	finalize_solver(solver_field, nr_sf);
	return(totalcount);  
      }


    //initialize spin fields on device
    convert2REAL4_spin<RealT>(solver_field[0],mixedsolveParameter.h2d_spin);

    cudaMemcpy(mixedsolveParameter.dev_spinin, mixedsolveParameter.h2d_spin, dev_spinsize, cudaMemcpyHostToDevice);
    printf("%s\n", cudaGetErrorString(cudaGetLastError()));

    // solve in single prec on device
    // D p_k = r_k
    printf("Entering inner solver\n");
    assert((startinner = clock())!=-1);
    totalcount += dev_cg<RealT,MixedsolveOperatorT>(mixedsolveParameter.dev_gf, mixedsolveParameter.dev_spinin, mixedsolveParameter.dev_spinout, mixedsolveParameter.dev_spin1, mixedsolveParameter.dev_spin2, mixedsolveParameter.dev_spin3, mixedsolveParameter.dev_spin4, mixedsolveParameter.dev_spin5, dev_grid,dev_nn, mixedsolveOperator, sourcesquarenorm, rel_prec, eps);
    stopinner = clock();
    timeelapsed = (double) (stopinner-startinner)/CLOCKS_PER_SEC;
    printf("Inner solver done\nTime elapsed: %.6e sec\n", timeelapsed);


    // copy back
    cudaMemcpy(mixedsolveParameter.h2d_spin, mixedsolveParameter.dev_spinout, dev_spinsize, cudaMemcpyDeviceToHost);
    printf("%s\n", cudaGetErrorString(cudaGetLastError()));
   
    convert2double_spin<RealT>(mixedsolveParameter.h2d_spin, solver_field[2]);
   
    add(solver_field[1],solver_field[1],solver_field[2],N);
    // x_(k+1) = x_k + p_k
   
    outercount ++;
    
  }// outer loop 

  printf("Did NOT reach solver precision of eps=%.2e\n",eps);
  //multiply with D^dagger
  mixedsolveOperator.checkDeinit(solver_field[1],solver_field[3],P,N);
  finalize_mixedsolve(&mixedsolveParameter);
       

  stop = clock();
  timeelapsed = (double) (stop-start)/CLOCKS_PER_SEC;
  printf("Inversion done in mixed precision.\n Number of iterations in outer solver: %d\n Squared residue: %.16e\n Time elapsed: %.6e sec\n", outercount, rk, timeelapsed);

  finalize_solver(solver_field, nr_sf);
  return(-1);
}


#include "mixedsolveOperator.cuh"


extern "C" int mixed_solve (spinor* const P, spinor* const Q, const int max_iter, double eps, const int rel_prec,const int N)
{ 
  MixedsolveOperatorDirac<REAL> mixedsolveOperator(0);
  return mixed_solveT<REAL ,MixedsolveOperatorDirac>(P,Q,max_iter,eps,rel_prec,N,mixedsolveOperator); 
}
extern "C" int mixed_solveD(spinor* const P, spinor* const Q, const int max_iter, double eps, const int rel_prec,const int N)
{ 
  #ifndef USETEXTURE
    MixedsolveOperatorDirac<REALD> mixedsolveOperator(0);
    return mixed_solveT<REALD,MixedsolveOperatorDirac>(P,Q,max_iter,eps,rel_prec,N,mixedsolveOperator);
  #else
    printf("WARNING: Using GPU/mixed_solve instead of double precision version.");
    return mixed_solve(P,Q,max_iter,eps,rel_prec,N);
  #endif
}

extern "C" int mixed_solve_DiracDaggerDirac (spinor* const P, spinor* const Q, const int max_iter, double eps, const int rel_prec,const int N)
{ 
  MixedsolveOperatorDiracDaggerDirac<REAL> mixedsolveOperator;
  return mixed_solveT<REAL ,MixedsolveOperatorDiracDaggerDirac>(P,Q,max_iter,eps,rel_prec,N,mixedsolveOperator); 
}
extern "C" int mixed_solve_DiracDaggerDiracD(spinor* const P, spinor* const Q, const int max_iter, double eps, const int rel_prec,const int N)
{ 
  #ifndef USETEXTURE
    MixedsolveOperatorDiracDaggerDirac<REALD> mixedsolveOperator;
    return mixed_solveT<REALD,MixedsolveOperatorDiracDaggerDirac>(P,Q,max_iter,eps,rel_prec,N,mixedsolveOperator);
  #else
    printf("WARNING: Using GPU/mixed_solve_DiracDaggerDirac instead of double precision version.");
    return mixed_solve_DiracDaggerDirac(P,Q,max_iter,eps,rel_prec,N);
  #endif
}

extern "C" int mixed_solve_DiracDaggerDiracDiracDaggerDirac (spinor* const P, spinor* const Q, const int max_iter, double eps, const int rel_prec,const int N)
{ 
  MixedsolveOperatorDiracDaggerDiracDiracDaggerDirac<REAL> mixedsolveOperator;
  return mixed_solveT<REAL ,MixedsolveOperatorDiracDaggerDiracDiracDaggerDirac>(P,Q,max_iter,eps,rel_prec,N,mixedsolveOperator); 
}
extern "C" int mixed_solve_DiracDaggerDiracDiracDaggerDiracD(spinor* const P, spinor* const Q, const int max_iter, double eps, const int rel_prec,const int N)
{ 
  #ifndef USETEXTURE
    MixedsolveOperatorDiracDaggerDiracDiracDaggerDirac<REALD> mixedsolveOperator;
    return mixed_solveT<REALD,MixedsolveOperatorDiracDaggerDiracDiracDaggerDirac>(P,Q,max_iter,eps,rel_prec,N,mixedsolveOperator);
  #else
    printf("WARNING: Using GPU/mixed_solve_DiracDaggerDiracDiracDaggerDirac instead of double precision version.");
    return mixed_solve_DiracDaggerDiracDiracDaggerDirac(P,Q,max_iter,eps,rel_prec,N);
  #endif
}









template<class RealT>
void benchmark(spinor * const Q,MixedsolveParameter<RealT>& mixedsolveParameter){
  
  double timeelapsed = 0.0;
  clock_t start, stop;
  int i;
  
  size_t dev_spinsize = 6*VOLUME/2 * sizeof(dev_spinorM(RealT)); // float4 even-odd !
  convert2REAL4_spin<RealT>(Q,mixedsolveParameter.h2d_spin);
  cudaMemcpy(mixedsolveParameter.dev_spinin, mixedsolveParameter.h2d_spin, dev_spinsize, cudaMemcpyHostToDevice);
  printf("%s\n", cudaGetErrorString(cudaGetLastError()));
  
  #ifndef MPI
    assert((start = clock())!=-1);
  #else
    start = MPI_Wtime();
  #endif  
  
  
  

 #ifdef USETEXTURE
  //Bind texture gf
  bind_texture_gf(mixedsolveParameter.dev_gf);
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
  cudaMemcpyToSymbol("dev_k0c", &h0, sizeof(h0)) ; 
  cudaMemcpyToSymbol("dev_k1c", &h1, sizeof(h1)) ; 
  cudaMemcpyToSymbol("dev_k2c", &h2, sizeof(h2)) ; 
  cudaMemcpyToSymbol("dev_k3c", &h3, sizeof(h3)) ;
  
  cudaMemcpyToSymbol("dev_mk0c", &mh0, sizeof(mh0)) ; 
  cudaMemcpyToSymbol("dev_mk1c", &mh1, sizeof(mh1)) ; 
  cudaMemcpyToSymbol("dev_mk2c", &mh2, sizeof(mh2)) ; 
  cudaMemcpyToSymbol("dev_mk3c", &mh3, sizeof(mh3)) ;  
  
  
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
  
      #ifdef MPI
           xchange_field_wrapper(mixedsolveParameter.dev_spinin, 0);
      #endif
      #ifdef USETEXTURE
         bind_texture_spin(mixedsolveParameter.dev_spinin,1);
      #endif
       //bind_texture_nn(dev_nn_eo);
      //cudaFuncSetCacheConfig(dev_Hopping_Matrix, cudaFuncCachePreferL1);
      dev_Hopping_Matrix<RealT> <<<griddim3, blockdim3>>>
             (mixedsolveParameter.dev_gf, mixedsolveParameter.dev_spinin, mixedsolveParameter.dev_spin_eo1, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0); //mixedsolveParameter.dev_spin_eo1 == even -> 0
       //unbind_texture_nn();
    #ifdef USETEXTURE             
      unbind_texture_spin(1);
    #endif

    #ifdef MPI
        xchange_field_wrapper(mixedsolveParameter.dev_spin_eo1, 0);
    #endif
       bind_texture_spin(mixedsolveParameter.dev_spin_eo1,1);
  //bind_texture_nn(dev_nn_oe);
   // cudaFuncSetCacheConfig(dev_Hopping_Matrix, cudaFuncCachePreferL1);
    dev_Hopping_Matrix<RealT> <<<griddim3, blockdim3>>>
            (mixedsolveParameter.dev_gf, mixedsolveParameter.dev_spin_eo1, mixedsolveParameter.dev_spinin, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1); 
  //unbind_texture_nn();
    #ifdef USETEXTURE
      unbind_texture_spin(1);
   #endif

  }  
  printf("%s\n", cudaGetErrorString(cudaGetLastError())); 
  printf("Done\n"); 
  
  
  
  #ifndef MPI
    assert((stop = clock())!=-1);
    timeelapsed = (double) (stop-start)/CLOCKS_PER_SEC;
    // x2 because 2x Hopping per iteration
    double benchres = 1608.0*2*(VOLUME/2)* 1000 / timeelapsed / 1.0e9;
    printf("Benchmark: %f Gflops\n", benchres); 
  #else
    stop = MPI_Wtime();
    timeelapsed = (double) (stop-start);
    // x2 because 2x Hopping per iteration
    double benchres = 1608.0*2*(g_nproc*VOLUME/2)* 1000 / timeelapsed / 1.0e9;
    if (g_proc_id == 0) {
      printf("Benchmark: %f Gflops\n", benchres); 
    }
  #endif  
  
  
  
  
   
  #ifdef USETEXTURE
    unbind_texture_gf();
  #endif
}



#ifdef MPI
template<class RealT>
void benchmark2(spinor * const Q,MixedsolveParameter<RealT>& mixedsolveParameter){
  
  double timeelapsed = 0.0;
  clock_t start, stop;
  int i;
  
  size_t dev_spinsize = 6*VOLUME/2 * sizeof(dev_spinorM(RealT)); // float4 even-odd !
  convert2REAL4_spin(Q,mixedsolveParameter.h2d_spin);
  cudaMemcpy(mixedsolveParameter.dev_spinin, mixedsolveParameter.h2d_spin, dev_spinsize, cudaMemcpyHostToDevice);
  printf("%s\n", cudaGetErrorString(cudaGetLastError()));
  
  #ifndef MPI
    assert((start = clock())!=-1);
  #else
    start = MPI_Wtime();
  #endif  
  
  
  

 #ifdef USETEXTURE
  //Bind texture gf
  bind_texture_gf(mixedsolveParameter.dev_gf);
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
  cudaMemcpyToSymbol("dev_k0c", &h0, sizeof(h0)) ; 
  cudaMemcpyToSymbol("dev_k1c", &h1, sizeof(h1)) ; 
  cudaMemcpyToSymbol("dev_k2c", &h2, sizeof(h2)) ; 
  cudaMemcpyToSymbol("dev_k3c", &h3, sizeof(h3)) ;
  
  cudaMemcpyToSymbol("dev_mk0c", &mh0, sizeof(mh0)) ; 
  cudaMemcpyToSymbol("dev_mk1c", &mh1, sizeof(mh1)) ; 
  cudaMemcpyToSymbol("dev_mk2c", &mh2, sizeof(mh2)) ; 
  cudaMemcpyToSymbol("dev_mk3c", &mh3, sizeof(mh3)) ;  
  
  
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
  
  
 int blockdim4 = BLOCK2;
 if( VOLUME/2 % blockdim4 == 0){
   gridsize = (int) VOLUME/2/blockdim4;
 }
 else{
   gridsize = (int) VOLUME/2/blockdim4 + 1;
 }
 int griddim4 = gridsize;  
  
  
  
  he_cg_init<<< 1, 1 >>> (dev_grid, (REAL) g_kappa, (REAL)(g_mu/(2.0*g_kappa)), h0,h1,h2,h3); 
  printf("%s\n", cudaGetErrorString(cudaGetLastError()));
  printf("Applying dev_Qtm_pm_psi 100 times\n");
  
  for(i=0; i<100; i++){
  
      
   dev_Qtm_pm_psi_mpi(mixedsolveParameter.dev_spinin, mixedsolveParameter.dev_spin_eo1, griddim3,blockdim3, griddim4, blockdim4);   
   
   dev_Qtm_pm_psi_mpi(mixedsolveParameter.dev_spin_eo1, mixedsolveParameter.dev_spinin, griddim3,blockdim3, griddim4, blockdim4); 

  }  
  printf("%s\n", cudaGetErrorString(cudaGetLastError())); 
  printf("Done\n"); 
  
  
  
  #ifndef MPI
    assert((stop = clock())!=-1);
    timeelapsed = (double) (stop-start)/CLOCKS_PER_SEC;
    // x8 because 8x Hopping per iteration
    double benchres = 1608.0*8*(VOLUME/2)* 100 / timeelapsed / 1.0e9;
    printf("Benchmark: %f Gflops\n", benchres); 
  #else
    stop = MPI_Wtime();
    timeelapsed = (double) (stop-start);
    // 8 because 8x Hopping per iteration
    double benchres = 1608.0*8*(g_nproc*VOLUME/2)* 100 / timeelapsed / 1.0e9;
    if (g_proc_id == 0) {
      printf("Benchmark: %f Gflops\n", benchres); 
    }
  #endif  
  
  
  
  
   
  #ifdef USETEXTURE
    unbind_texture_gf();
  #endif
}

#endif









#else
extern "C" int mixed_solve (spinor * const P, spinor * const Q, const int max_iter, 
           double eps, const int rel_prec,const int N){
   printf("WARNING dummy function mixed_solve called\n");
   return(0);           
}

#endif
// WORK TO DO:
// Separate half and non-half inner solvers in a more transparent way!!




template<class RealT>
int mixed_solve_eoT (spinor * const P, spinor * const Q, const int max_iter, 
       	            double eps, const int rel_prec, const int N){

  // source in Q, initial solution in P (not yet implemented)
  double rk;
  int outercount=0;
  int totalcount=0;
  clock_t start, stop, startinner, stopinner; 
  double timeelapsed = 0.0;
  double sourcesquarenorm;
  int iter;//never referenced: , retval;
  spinor ** solver_field = NULL;
  const int nr_sf = 4;
  
  init_solver_field(&solver_field, VOLUMEPLUSRAND/2, nr_sf);
  size_t dev_spinsize;
  #ifndef HALF
    dev_spinsize = 6*VOLUME/2 * sizeof(dev_spinorM(RealT)); // float4 even-odd !
  #else
    dev_spinsize = 6*VOLUME/2 * sizeof(dev_spinor_half); //short4 eo !
    size_t dev_normsize = VOLUME/2 * sizeof(float);
  #endif
  MixedsolveParameter<RealT>& mixedsolveParameter=*init_mixedsolve_eo<RealT>(g_gauge_field);

  
 /* 
  #ifndef HALF
  // small benchmark
    assign(solver_field[0],Q,N);
    #ifndef MPI
      benchmark(solver_field[0]);
    #else
      benchmark2(solver_field[0]); 
    #endif
  // end small benchmark
  
   //exit(100);
  
  #endif //not HALF
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
  assign(solver_field[0],Q,N);
  printf("Initial residue: %.16e\n",rk);
  zero_spinor_field(solver_field[1],  N);//spin2 = x_k
  zero_spinor_field(solver_field[2],  N);
  printf("The VOLUME/2 is: %d\n",N);

  
  double norm = sqrt(_spinor_prod_re(Q[0],Q[0]));
  printf("norm source(0): %f\n", norm);


  //#include "test.sqz"

for(iter=0; iter<max_iter; iter++){

   printf("Applying double precision EO Dirac-Op Q_{-}Q{+}...\n");
   
   Qtm_pm_psi(solver_field[3], solver_field[2]);
   diff(solver_field[0],solver_field[0],solver_field[3],N);
    // r_k = b - D x_k
   
   rk = square_norm(solver_field[0], N, 1);
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
     Qtm_minus_psi(solver_field[3], solver_field[1]);
     assign(P, solver_field[3], N);

     printf("EO Inversion done in mixed precision.\n Number of iterations in outer solver: %d\n Squared residue: %.16e\n Time elapsed: %.6e sec\n", outercount, rk, timeelapsed);
   
  
     stop = clock();
     timeelapsed = (double) (stop-start)/CLOCKS_PER_SEC;
        
     finalize_mixedsolve(&mixedsolveParameter);
     finalize_solver(solver_field, nr_sf);
     return(totalcount);  
   }
   
  //initialize spin fields on device
  #ifndef HALF
    convert2REAL4_spin<RealT>(solver_field[0],mixedsolveParameter.h2d_spin);
  #else
    convert2REAL4_spin_half(solver_field[0],mixedsolveParameter.h2d_spin, mixedsolveParameter.h2d_spin_norm);
  #endif
  cudaMemcpy(mixedsolveParameter.dev_spinin, mixedsolveParameter.h2d_spin, dev_spinsize, cudaMemcpyHostToDevice);

  // also copy half spinor norm
  #ifdef HALF
  cudaMemcpy(mixedsolveParameter.dev_spinin_norm, mixedsolveParameter.h2d_spin_norm, dev_normsize, cudaMemcpyHostToDevice);
  #endif

  printf("%s\n", cudaGetErrorString(cudaGetLastError()));

   // solve in single prec on device
   // D p_k = r_k
   printf("Entering inner solver\n");
   assert((startinner = clock())!=-1);
   #ifndef HALF
      totalcount += dev_cg_eo(mixedsolveParameter.dev_gf, mixedsolveParameter.dev_spinin, mixedsolveParameter.dev_spinout, mixedsolveParameter.dev_spin1, mixedsolveParameter.dev_spin2, mixedsolveParameter.dev_spin3, mixedsolveParameter.dev_spin4, mixedsolveParameter.dev_spin5, dev_grid,dev_nn, (RealT) finaleps, mixedsolveParameter);
   #else
      totalcount += dev_cg_eo_half(mixedsolveParameter.dev_gf, 
                 mixedsolveParameter.dev_spinin, mixedsolveParameter.dev_spinin_norm,
                 mixedsolveParameter.dev_spinout,mixedsolveParameter.dev_spinout_norm,
                 mixedsolveParameter.dev_spin1, mixedsolveParameter.dev_spin1_norm,
                 mixedsolveParameter.dev_spin2, mixedsolveParameter.dev_spin2_norm,
                 mixedsolveParameter.dev_spin3, mixedsolveParameter.dev_spin3_norm,
                 mixedsolveParameter.dev_spin4, mixedsolveParameter.dev_spin4_norm,
                 mixedsolveParameter.dev_spin5, mixedsolveParameter.dev_spin5_norm,
                 dev_grid,dev_nn, (RealT) finaleps,
                 mixedsolveParameter); 
   #endif
   stopinner = clock();
   timeelapsed = (double) (stopinner-startinner)/CLOCKS_PER_SEC;
   printf("Inner solver done\nTime elapsed: %.6e sec\n", timeelapsed);
 
   // copy back
   cudaMemcpy(mixedsolveParameter.h2d_spin, mixedsolveParameter.dev_spinout, dev_spinsize, cudaMemcpyDeviceToHost);
   #ifdef HALF
     cudaMemcpy(mixedsolveParameter.h2d_spin_norm, mixedsolveParameter.dev_spinout_norm, dev_normsize, cudaMemcpyDeviceToHost);
   #endif
   printf("%s\n", cudaGetErrorString(cudaGetLastError()));

   #ifndef HALF
     convert2double_spin<RealT>(mixedsolveParameter.h2d_spin, solver_field[2]);
   #else
     convert2double_spin_half(mixedsolveParameter.h2d_spin, mixedsolveParameter.h2d_spin_norm, solver_field[2]);
   #endif

   // x_(k+1) = x_k + p_k
   add(solver_field[1],solver_field[1],solver_field[2],N);

   outercount ++;
}// outer loop 
    
     printf("Did NOT reach solver precision of eps=%.2e\n",eps);
     //multiply with Qtm_minus_psi (for non gpu done in invert_eo.c)
     Qtm_minus_psi(solver_field[3], solver_field[1]);
     assign(P, solver_field[3], N);
    

    assert((stop = clock())!=-1);
    timeelapsed = (double) (stop-start)/CLOCKS_PER_SEC;
    printf("Inversion done in mixed precision.\n Number of iterations in outer solver: %d\n Squared residue: %.16e\n Time elapsed: %.6e sec\n", outercount, rk, timeelapsed);

  finalize_mixedsolve(&mixedsolveParameter);
  finalize_solver(solver_field, nr_sf);
  return(-1);
}

extern "C" int mixed_solve_eo  (spinor * const P, spinor * const Q, const int max_iter, double eps, const int rel_prec, const int N)
{ return mixed_solve_eoT<REAL >(P,Q,max_iter,eps,rel_prec,N); };
#ifndef HALF
  #ifndef USETEXTURE
    extern "C" int mixed_solve_eoD (spinor * const P, spinor * const Q, const int max_iter, double eps, const int rel_prec, const int N)
    { return mixed_solve_eoT<REALD>(P,Q,max_iter,eps,rel_prec,N); }
  #endif
#endif












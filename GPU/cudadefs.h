
#ifndef _CUDADEF_H
#define _CUDADEF_H

#define ACCUM_N 2048
#define DOTPROD_DIM 128
#define GF_8
#define TEMPORALGAUGE
#define USETEXTURE
//#define HALF
#define GPU_DOUBLE
#define RELATIVISTIC_BASIS 
//#define GPU_3DBLOCK
#define LOWOUTPUT
//#define CUDA_45


#define BLOCK 288 // Block Size	
					// dev_Hopping_Matrix<<<>>>()
#define BLOCKSUB 8


#define BLOCK2 256 // Block Size 2 for dev_mul_one_pm...		// dev_mul_one_pm_imu_inv<<<>>>()
#define BLOCK3 512							// dev_copy_spinor_field<<<>>>(), dev_zero_spinor_field<<<>>>()
#define REDUCTION_N 128 // Block size for reduction operations

#define BLOCKGAUGE 128


#define maxblockdim 512

// offset in the fields (mem-space between two subsequent indices at same space-time point)
#define DEVOFF dev_Offset
#define HOSTOFF host_Offset


//Block sizes for GPU_DOUBLE
#define BLOCKD 288 // Block Size	
#define BLOCK2D 512

//this defines a stopping criterion for the inner mixed solver
#define SP_MIN_EPS 1.0e-20

#endif



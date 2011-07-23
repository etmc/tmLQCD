#ifndef _CUDADEF_H
#define _CUDADEF_H

#define cublasStatus_t cublasStatus//necessary for cublas>=4.0

#define ACCUM_N 2048
#define DOTPROD_DIM 128

//#define GF_8
//#define TEMPORALGAUGE
//#define USETEXTURE
//#define HALF
#define OLD_CUBLAS //cublas older than 4.0 ?

#define REAL   float
#define REALD  double
#define REAL4  float4
#define REAL4D double4


#define BLOCK 64//192 // Block Size						// dev_Hopping_Matrix<<<>>>()
#define BLOCK2 64//320 // Block Size 2 for dev_mul_one_pm...		// dev_mul_one_pm_imu_inv<<<>>>()
#define BLOCK3 64//128							// dev_copy_spinor_field<<<>>>(), dev_zero_spinor_field<<<>>>()
#define REDUCTION_N 64//512 // Block size for reduction operations //old: 512



#define maxblockdim 64//512


#endif



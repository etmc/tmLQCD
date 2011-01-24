
#ifndef _CUDADEF_H
#define _CUDADEF_H

#define ACCUM_N 2048
#define DOTPROD_DIM 128

#define GF_8
//#define TEMPORALGAUGE
//#define USETEXTURE
//#define HALF


#define REAL float
#define BLOCK 192 // Block Size						// dev_Hopping_Matrix<<<>>>()
#define BLOCK2 320 // Block Size 2 for dev_mul_one_pm...		// dev_mul_one_pm_imu_inv<<<>>>()
#define BLOCK3 128							// dev_copy_spinor_field<<<>>>(), dev_zero_spinor_field<<<>>>()
#define REDUCTION_N 512 // Block size for reduction operations



#define maxblockdim 512


#endif



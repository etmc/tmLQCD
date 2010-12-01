
#ifndef _CUDADEF_H
#define _CUDADEF_H

#define ACCUM_N 2048
#define DOTPROD_DIM 128

#define GF_8
#define TEMPORALGAUGE
#define USETEXTURE
//#define HALF


#define REAL float
#define BLOCK 192 // Block Size
#define BLOCK2 320 // Block Size 2 for dev_mul_one_pm...
#define REDUCTION_N 512 // Block size for reduction operations



#define maxblockdim 512


#endif



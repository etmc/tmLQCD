#include "cudadefs.h"

/* GPU Stuff */

typedef struct  {
  REAL re;
  REAL im;
} dev_complex;

/* Device Gauge Fields */
typedef dev_complex dev_su3 [3][3];  /* su(3)-Matrix 3x3 komplexe Einträge DEVICE */

typedef struct{
  dev_su3 m;
  float pad;
} dev_su3_pad;


 
typedef float4 su3_2v ;  /* 2 Zeilen der su(3)-Matrix, 6 komplexe Einträge HOST  
                           3*4*VOLUME in array -> texture */
typedef float4 dev_su3_2v ;  /* 2 Zeilen der su(3)-Matrix 
                             3*2 komplexe Einträge DEVICE 3*4*VOLUME in array -> texture*/
                             

typedef float4 dev_su3_8 ;  /* 8 numbers to reconstruct the gauge field as described in M. Clark */   
                                              
/* Device Spinor Fields */
typedef float4 dev_spinor;
typedef struct dev_spinor_smem{
  dev_spinor spin;
  REAL dummy;
} dev_spinor_smem;
typedef dev_complex dev_propmatrix[12][12];
typedef dev_complex dev_fbyf[4][4];




/* END GPU Stuff */


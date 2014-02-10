#include "cudadefs.h"

/* GPU Stuff */

typedef struct  {
  float re;
  float im;
} dev_complex;



/* Device Gauge Fields */
typedef dev_complex dev_su3 [3][3];  /* su(3)-Matrix 3x3 komplexe Einträge DEVICE */


typedef struct{
  dev_su3 m;
//  char pad; 
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
  float dummy;
} dev_spinor_smem;
typedef dev_complex dev_propmatrix[12][12];
typedef dev_complex dev_fbyf[4][4];

#ifdef HALF
 typedef short4 dev_spinor_half;
 typedef short4 dev_su3_2v_half;
 typedef short4 dev_su3_8_half;
#endif



/* double gauge and momentum fields for gauge_monomial.cuh */
typedef double2 dev_su3_2v_d ;
typedef struct
{
   double d1,d2,d3,d4,d5,d6,d7,d8;
} dev_su3adj;
typedef struct  {
  double re;
  double im;
} dev_complex_d;
typedef dev_complex_d dev_su3_d [3][3]; 


typedef struct {
dev_complex_d c0;
dev_complex_d c1;
dev_complex_d c2;
} dev_su3_vec_d;

typedef dev_complex dev_su3_vec [3];


typedef double2 dev_spinor_d;




/* END GPU Stuff */


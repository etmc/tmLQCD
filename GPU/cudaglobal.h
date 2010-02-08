/***********************************************************************
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
 ***********************************************************************/
/* GPU Stuff */


#define REAL float
#define BLOCK 128 // Block Size 
#define maxblockdim 512



typedef struct __align__(8) {
  REAL re;
  REAL im;
} dev_complex;

/* Device Gauge Fields */
typedef dev_complex dev_su3 [3][3];  /* su(3)-Matrix 3x3 komplexe Einträge DEVICE */
typedef float4 su3_2v ;  /* 2 Zeilen der su(3)-Matrix, 6 komplexe Einträge HOST  
                           3*4*VOLUME in array -> texture */
typedef float4 dev_su3_2v ;  /* 2 Zeilen der su(3)-Matrix 
                             3*2 komplexe Einträge DEVICE 3*4*VOLUME in array -> texture*/
                             

typedef float4 dev_su3_8 ;  /* 8 numbers to reconstruc the gauge field as described in M. Clark */   
                                              
/* Device Spinor Fields */
typedef float4 dev_spinor;
typedef struct dev_spinor_smem{
  dev_spinor spin;
  REAL dummy;
} dev_spinor_smem;
typedef dev_complex dev_propmatrix[12][12];
typedef dev_complex dev_fbyf[4][4];




/* END GPU Stuff */


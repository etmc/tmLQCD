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

/***********************************************************************
 * This provides a NVIDIA CUDA implementation of a mixed-precision solver 
************************************************************************/

#include <cuda.h>
#include <cuda_runtime.h>
#include "cublas.h"
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "global.h"
#include "cudaglobal.h"
#include "mixed_solve.h"
#include <math.h>


extern "C" {
#include "../tm_operators.h"
#include "../linalg_eo.h"
#include "../start.h"
#include "../complex.h"
#include "../read_input.h"
#include "../geometry_eo.h"
#include "../boundary.h"
}

#define ACCUM_N 2048
#define DOTPROD_DIM 128




int g_numofgpu;


dev_su3_2v * dev_gf;
dev_su3_2v * h2d_gf;



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


__device__  int  dev_LX,dev_LY,dev_LZ,dev_T,dev_VOLUME;


 /* texture for gauge field */
 texture<float4,1, cudaReadModeElementType> gf_tex;
 const textureReference* gf_texRefPtr = NULL;
 cudaChannelFormatDesc gf_channelDesc;
 
 /* texture for spinor field */
 texture<float4,1, cudaReadModeElementType> spin_tex;
 const textureReference* spin_texRefPtr = NULL;
 cudaChannelFormatDesc spin_channelDesc;

 /* texture for spinor field 2*/
 texture<float4,1, cudaReadModeElementType> spin_tex2;
 const textureReference* spin_texRefPtr2 = NULL;
 cudaChannelFormatDesc spin_channelDesc2;

__device__ inline dev_complex dev_cconj (dev_complex c){ /*konjugiert komplexe Zahl*/
 dev_complex erg;
 erg.re = c.re;
 erg.im = -1.0*c.im;
return erg;
}

__device__ inline void dev_ccopy(dev_complex* von, dev_complex* nach){/*kopiert complex von nach complex nach*/
  nach->re = von->re;
  nach->im = von->im;
}

__device__ inline REAL dev_cabssquare (dev_complex c){ /*gibt abs^2 einer komplexen Zahl zurück*/
 return c.re*c.re + c.im*c.im;
}

__device__ inline REAL dev_cabsolute (dev_complex c){/*gibt Betrag einer kompl. zahl zurück*/
 return sqrt(c.re*c.re + c.im*c.im);
}


__device__ inline  dev_complex dev_crealmult(dev_complex c1, REAL real){ /*multipliziert c1 mit reeller zahl re*/
  dev_complex erg;
  erg.re = real*c1.re;
  erg.im = real*c1.im;
return erg;
}

__device__ inline dev_complex dev_cmult (dev_complex c1, dev_complex c2){ /*multiplizier zwei komplexe Zahlen*/
  dev_complex erg;
  erg.re = c1.re * c2.re - c1.im * c2.im;
  erg.im = c1.re * c2.im + c1.im * c2.re;
return erg;
}

__device__ inline dev_complex dev_cadd (dev_complex c1, dev_complex c2){ /*addiert zwei komplexe Zahlen */
  dev_complex erg;
  erg.re = c1.re + c2.re;
  erg.im = c1.im + c2.im;
return erg;
}


__device__ inline dev_complex dev_cdiv(dev_complex c1, dev_complex c2) { /* dividiert c1 durch c2 */
  dev_complex erg;
  REAL oneovernenner = 1.0/(c2.re*c2.re + c2.im*c2.im);
  erg.re = oneovernenner*(c1.re*c2.re + c1.im*c2.im);
  erg.im = oneovernenner*(c1.im*c2.re - c1.re*c2.im);
return erg;
}


__device__ inline dev_complex dev_csub(dev_complex c1, dev_complex c2){
   dev_complex erg;
   erg.re = c1.re - c2.re;
   erg.im = c1.im - c2.im;
return erg;
}


__device__ inline dev_complex dev_initcomplex(REAL re, REAL im){/* gibt komplexe Zahl mit Realt re und Imt im zurück*/
    dev_complex erg;
    erg.re = re;
    erg.im = im;
return (erg);
}





__device__ inline void dev_copy_spinor(dev_spinor *i1, dev_spinor *i2){
  int i;
  #pragma unroll 6
  for(i=0;i<6;i++){ //color + spin
    (*(i2+i)).x = (*(i1+i)).x;
    (*(i2+i)).y = (*(i1+i)).y;
    (*(i2+i)).z = (*(i1+i)).z;
    (*(i2+i)).w = (*(i1+i)).w;
  }
}

__device__ inline void dev_zero_spinor(dev_spinor *sin){
  int i;
  #pragma unroll 6
  for(i=0;i<6;i++){ //color + spin
    (*(sin+i)).x = 0.0;
    (*(sin+i)).y = 0.0;
    (*(sin+i)).z = 0.0;
    (*(sin+i)).w = 0.0;
  }
}






//out = in + lambda in2
__device__ inline void dev_skalarmult_add_assign_spinor(dev_spinor *in, REAL lambda,dev_spinor * in2, dev_spinor * out){
  int i; 
  #pragma unroll 6
for(i=0;i<6;i++){ //color + spin
    (*(out+i)).x = (*(in+i)).x + lambda* (*(in2+i)).x;
    (*(out+i)).y = (*(in+i)).y + lambda* (*(in2+i)).y;
    (*(out+i)).z = (*(in+i)).z + lambda* (*(in2+i)).z;
    (*(out+i)).w = (*(in+i)).w + lambda* (*(in2+i)).w;
  }
}




//out = in + lambda in2
__device__ inline void dev_complexmult_add_assign_spinor(dev_spinor * in, dev_complex lambda,dev_spinor * in2, dev_spinor * out){
  int i;
  #pragma unroll 6
  for(i=0;i<6;i++){ //color + spin
    (*(out+i)).x = (*(in+i)).x + ((*(in2+i)).x*lambda.re - (*(in2+i)).y*lambda.im);
    (*(out+i)).y = (*(in+i)).y + ((*(in2+i)).x*lambda.im + (*(in2+i)).y*lambda.re);
    (*(out+i)).z = (*(in+i)).z + ((*(in2+i)).z*lambda.re - (*(in2+i)).w*lambda.im);
    (*(out+i)).w = (*(in+i)).w + ((*(in2+i)).z*lambda.im + (*(in2+i)).w*lambda.re);
  }
}




//out = in + (lambda)* in2
__device__ inline void dev_complexcgmult_add_assign_spinor(dev_spinor * in, dev_complex lambda,dev_spinor * in2, dev_spinor * out){
  int i;
  #pragma unroll 6
  for(i=0;i<6;i++){ //color + spin
    (*(out+i)).x = (*(in+i)).x + ((*(in2+i)).x*lambda.re + (*(in2+i)).y*lambda.im);
    (*(out+i)).y = (*(in+i)).y + (-(*(in2+i)).x*lambda.im + (*(in2+i)).y*lambda.re);
    (*(out+i)).z = (*(in+i)).z + ((*(in2+i)).z*lambda.re + (*(in2+i)).w*lambda.im);
    (*(out+i)).w = (*(in+i)).w + (-(*(in2+i)).z*lambda.im + (*(in2+i)).w*lambda.re);
  }
}



__device__ void inline dev_skalarmult_spinor(dev_spinor * in, dev_complex lambda, dev_spinor * out){
  int i;
  #pragma unroll 6
  for(i=0;i<6;i++){ //color + spin
    //out[i] = dev_cmult(in[i],lambda);
    
    (*(out+i)).x = (*(in+i)).x*lambda.re - (*(in+i)).y*lambda.im;
    (*(out+i)).y = (*(in+i)).y*lambda.re + (*(in+i)).x*lambda.im;
    
    (*(out+i)).z = (*(in+i)).z*lambda.re - (*(in+i)).w*lambda.im;
    (*(out+i)).w = (*(in+i)).w*lambda.re + (*(in+i)).z*lambda.im;
  }
}



__device__ void inline dev_realmult_spinor(dev_spinor * in, REAL lambda){
  int i;
  #pragma unroll 6
  for(i=0;i<6;i++){ //color + spin
    //in[i] = in[i]*lambda;
    (*(in+i)).x = (*(in+i)).x*lambda;
    (*(in+i)).y = (*(in+i)).y*lambda;
    
    (*(in+i)).z = (*(in+i)).z*lambda;
    (*(in+i)).w = (*(in+i)).w*lambda;
  }
}





__device__ inline void dev_add_spinor_assign(dev_spinor * i1, dev_spinor * i2){
  int i;
  #pragma unroll 6
  for(i=0;i<6;i++){ //color + spin
    (*(i1+i)).x = (*(i1+i)).x + (*(i2+i)).x;
    (*(i1+i)).y = (*(i1+i)).y + (*(i2+i)).y;
    (*(i1+i)).z = (*(i1+i)).z + (*(i2+i)).z;
    (*(i1+i)).w = (*(i1+i)).w + (*(i2+i)).w;
  }
}



__device__ inline void dev_sub_spinor_assign(dev_spinor * i1, dev_spinor * i2){
  int i;
  #pragma unroll 6
  for(i=0;i<6;i++){ //color + spin
    (*(i1+i)).x = (*(i1+i)).x - (*(i2+i)).x;
    (*(i1+i)).y = (*(i1+i)).y - (*(i2+i)).y;
    (*(i1+i)).z = (*(i1+i)).z - (*(i2+i)).z;
    (*(i1+i)).w = (*(i1+i)).w - (*(i2+i)).w;
  }
}





//multipliziert su3-Matrix mal Spinor im Dirac-Raum
//code in su3_MtV.txt -- generated with codegen
__device__ void dev_su3MtV_spintex(dev_su3 M, int pos, dev_spinor * out){

dev_spinor s1, s2;

s1 = tex1Dfetch(spin_tex,6*pos);
s2 = tex1Dfetch(spin_tex,6*pos+1);


(*(out+0)).x =  ( M[0][0].re*s1.x - M[0][0].im*s1.y ) + ( M[0][1].re*s1.z - M[0][1].im*s1.w ) + ( M[0][2].re*s2.x - M[0][2].im*s2.y );
(*(out+0)).y = ( M[0][0].re*s1.y + M[0][0].im*s1.x ) + ( M[0][1].re*s1.w + M[0][1].im*s1.z ) + ( M[0][2].re*s2.y + M[0][2].im*s2.x );


/* checked by look of eye */
(*(out+0)).z =  ( M[1][0].re*s1.x - M[1][0].im*s1.y ) + ( M[1][1].re*s1.z - M[1][1].im*s1.w ) + ( M[1][2].re*s2.x - M[1][2].im*s2.y );
/* checked */

(*(out+0)).w =  ( M[1][0].re*s1.y + M[1][0].im*s1.x ) + ( M[1][1].re*s1.w + M[1][1].im*s1.z ) + ( M[1][2].re*s2.y + M[1][2].im*s2.x );


(*(out+1)).x = ( M[2][0].re*s1.x - M[2][0].im*s1.y ) + ( M[2][1].re*s1.z - M[2][1].im*s1.w ) + ( M[2][2].re*s2.x - M[2][2].im*s2.y );
(*(out+1)).y =  ( M[2][0].re*s1.y + M[2][0].im*s1.x ) + ( M[2][1].re*s1.w + M[2][1].im*s1.z ) + ( M[2][2].re*s2.y + M[2][2].im*s2.x );



s1 = tex1Dfetch(spin_tex,6*pos+2);
(*(out+1)).z =  ( M[0][0].re*s2.z - M[0][0].im*s2.w ) + ( M[0][1].re*s1.x - M[0][1].im*s1.y ) + ( M[0][2].re*s1.z - M[0][2].im*s1.w );
(*(out+1)).w =  ( M[0][0].re*s2.w + M[0][0].im*s2.z ) + ( M[0][1].re*s1.y + M[0][1].im*s1.x ) + ( M[0][2].re*s1.w + M[0][2].im*s1.z );


(*(out+2)).x = ( M[1][0].re*s2.z - M[1][0].im*s2.w ) + ( M[1][1].re*s1.x - M[1][1].im*s1.y ) + ( M[1][2].re*s1.z - M[1][2].im*s1.w );
(*(out+2)).y =  ( M[1][0].re*s2.w + M[1][0].im*s2.z ) + ( M[1][1].re*s1.y + M[1][1].im*s1.x ) + ( M[1][2].re*s1.w + M[1][2].im*s1.z );


(*(out+2)).z =  ( M[2][0].re*s2.z - M[2][0].im*s2.w ) + ( M[2][1].re*s1.x - M[2][1].im*s1.y ) + ( M[2][2].re*s1.z - M[2][2].im*s1.w );
(*(out+2)).w =  ( M[2][0].re*s2.w + M[2][0].im*s2.z ) + ( M[2][1].re*s1.y + M[2][1].im*s1.x ) + ( M[2][2].re*s1.w + M[2][2].im*s1.z );



s1 = tex1Dfetch(spin_tex,6*pos+3);
s2 = tex1Dfetch(spin_tex,6*pos+4);
(*(out+3)).x =  ( M[0][0].re*s1.x - M[0][0].im*s1.y ) + ( M[0][1].re*s1.z - M[0][1].im*s1.w ) + ( M[0][2].re*s2.x - M[0][2].im*s2.y );
(*(out+3)).y =   ( M[0][0].re*s1.y + M[0][0].im*s1.x ) + ( M[0][1].re*s1.w + M[0][1].im*s1.z ) + ( M[0][2].re*s2.y + M[0][2].im*s2.x );


(*(out+3)).z =  ( M[1][0].re*s1.x - M[1][0].im*s1.y ) + ( M[1][1].re*s1.z - M[1][1].im*s1.w ) + ( M[1][2].re*s2.x - M[1][2].im*s2.y );
(*(out+3)).w =  ( M[1][0].re*s1.y + M[1][0].im*s1.x ) + ( M[1][1].re*s1.w + M[1][1].im*s1.z ) + ( M[1][2].re*s2.y + M[1][2].im*s2.x );


(*(out+4)).x =  ( M[2][0].re*s1.x - M[2][0].im*s1.y ) + ( M[2][1].re*s1.z - M[2][1].im*s1.w ) + ( M[2][2].re*s2.x - M[2][2].im*s2.y );
(*(out+4)).y =  ( M[2][0].re*s1.y + M[2][0].im*s1.x ) + ( M[2][1].re*s1.w + M[2][1].im*s1.z ) + ( M[2][2].re*s2.y + M[2][2].im*s2.x );



s1 = tex1Dfetch(spin_tex,6*pos+5);
(*(out+4)).z =  ( M[0][0].re*s2.z - M[0][0].im*s2.w ) + ( M[0][1].re*s1.x - M[0][1].im*s1.y ) + ( M[0][2].re*s1.z - M[0][2].im*s1.w );
(*(out+4)).w =   ( M[0][0].re*s2.w + M[0][0].im*s2.z ) + ( M[0][1].re*s1.y + M[0][1].im*s1.x ) + ( M[0][2].re*s1.w + M[0][2].im*s1.z );


(*(out+5)).x = ( M[1][0].re*s2.z - M[1][0].im*s2.w ) + ( M[1][1].re*s1.x - M[1][1].im*s1.y ) + ( M[1][2].re*s1.z - M[1][2].im*s1.w );
(*(out+5)).y =  ( M[1][0].re*s2.w + M[1][0].im*s2.z ) + ( M[1][1].re*s1.y + M[1][1].im*s1.x ) + ( M[1][2].re*s1.w + M[1][2].im*s1.z );


(*(out+5)).z =  ( M[2][0].re*s2.z - M[2][0].im*s2.w ) + ( M[2][1].re*s1.x - M[2][1].im*s1.y ) + ( M[2][2].re*s1.z - M[2][2].im*s1.w );
(*(out+5)).w =  ( M[2][0].re*s2.w + M[2][0].im*s2.z ) + ( M[2][1].re*s1.y + M[2][1].im*s1.x ) + ( M[2][2].re*s1.w + M[2][2].im*s1.z );


}










//multipliziert su3-Matrix mal Spinor im Dirac-Raum
//code in su3_MtV.txt -- generated with codegen
__device__ void dev_su3MtV(dev_su3 M, dev_spinor * s, dev_spinor * out){

(*(out+0)).x =  ( M[0][0].re*(*(s+0)).x - M[0][0].im*(*(s+0)).y ) + ( M[0][1].re*(*(s+0)).z - M[0][1].im*(*(s+0)).w ) + ( M[0][2].re*(*(s+1)).x - M[0][2].im*(*(s+1)).y );
(*(out+0)).y = ( M[0][0].re*(*(s+0)).y + M[0][0].im*(*(s+0)).x ) + ( M[0][1].re*(*(s+0)).w + M[0][1].im*(*(s+0)).z ) + ( M[0][2].re*(*(s+1)).y + M[0][2].im*(*(s+1)).x );


(*(out+0)).z =  ( M[1][0].re*(*(s+0)).x - M[1][0].im*(*(s+0)).y ) + ( M[1][1].re*(*(s+0)).z - M[1][1].im*(*(s+0)).w ) + ( M[1][2].re*(*(s+1)).x - M[1][2].im*(*(s+1)).y );
(*(out+0)).w =  ( M[1][0].re*(*(s+0)).y + M[1][0].im*(*(s+0)).x ) + ( M[1][1].re*(*(s+0)).w + M[1][1].im*(*(s+0)).z ) + ( M[1][2].re*(*(s+1)).y + M[1][2].im*(*(s+1)).x );


(*(out+1)).x = ( M[2][0].re*(*(s+0)).x - M[2][0].im*(*(s+0)).y ) + ( M[2][1].re*(*(s+0)).z - M[2][1].im*(*(s+0)).w ) + ( M[2][2].re*(*(s+1)).x - M[2][2].im*(*(s+1)).y );
(*(out+1)).y =  ( M[2][0].re*(*(s+0)).y + M[2][0].im*(*(s+0)).x ) + ( M[2][1].re*(*(s+0)).w + M[2][1].im*(*(s+0)).z ) + ( M[2][2].re*(*(s+1)).y + M[2][2].im*(*(s+1)).x );


(*(out+1)).z =  ( M[0][0].re*(*(s+1)).z - M[0][0].im*(*(s+1)).w ) + ( M[0][1].re*(*(s+2)).x - M[0][1].im*(*(s+2)).y ) + ( M[0][2].re*(*(s+2)).z - M[0][2].im*(*(s+2)).w );
(*(out+1)).w =  ( M[0][0].re*(*(s+1)).w + M[0][0].im*(*(s+1)).z ) + ( M[0][1].re*(*(s+2)).y + M[0][1].im*(*(s+2)).x ) + ( M[0][2].re*(*(s+2)).w + M[0][2].im*(*(s+2)).z );


(*(out+2)).x = ( M[1][0].re*(*(s+1)).z - M[1][0].im*(*(s+1)).w ) + ( M[1][1].re*(*(s+2)).x - M[1][1].im*(*(s+2)).y ) + ( M[1][2].re*(*(s+2)).z - M[1][2].im*(*(s+2)).w );
(*(out+2)).y =  ( M[1][0].re*(*(s+1)).w + M[1][0].im*(*(s+1)).z ) + ( M[1][1].re*(*(s+2)).y + M[1][1].im*(*(s+2)).x ) + ( M[1][2].re*(*(s+2)).w + M[1][2].im*(*(s+2)).z );


(*(out+2)).z =  ( M[2][0].re*(*(s+1)).z - M[2][0].im*(*(s+1)).w ) + ( M[2][1].re*(*(s+2)).x - M[2][1].im*(*(s+2)).y ) + ( M[2][2].re*(*(s+2)).z - M[2][2].im*(*(s+2)).w );
(*(out+2)).w =  ( M[2][0].re*(*(s+1)).w + M[2][0].im*(*(s+1)).z ) + ( M[2][1].re*(*(s+2)).y + M[2][1].im*(*(s+2)).x ) + ( M[2][2].re*(*(s+2)).w + M[2][2].im*(*(s+2)).z );


(*(out+3)).x =  ( M[0][0].re*(*(s+3)).x - M[0][0].im*(*(s+3)).y ) + ( M[0][1].re*(*(s+3)).z - M[0][1].im*(*(s+3)).w ) + ( M[0][2].re*(*(s+4)).x - M[0][2].im*(*(s+4)).y );
(*(out+3)).y =   ( M[0][0].re*(*(s+3)).y + M[0][0].im*(*(s+3)).x ) + ( M[0][1].re*(*(s+3)).w + M[0][1].im*(*(s+3)).z ) + ( M[0][2].re*(*(s+4)).y + M[0][2].im*(*(s+4)).x );


(*(out+3)).z =  ( M[1][0].re*(*(s+3)).x - M[1][0].im*(*(s+3)).y ) + ( M[1][1].re*(*(s+3)).z - M[1][1].im*(*(s+3)).w ) + ( M[1][2].re*(*(s+4)).x - M[1][2].im*(*(s+4)).y );
(*(out+3)).w =  ( M[1][0].re*(*(s+3)).y + M[1][0].im*(*(s+3)).x ) + ( M[1][1].re*(*(s+3)).w + M[1][1].im*(*(s+3)).z ) + ( M[1][2].re*(*(s+4)).y + M[1][2].im*(*(s+4)).x );


(*(out+4)).x =  ( M[2][0].re*(*(s+3)).x - M[2][0].im*(*(s+3)).y ) + ( M[2][1].re*(*(s+3)).z - M[2][1].im*(*(s+3)).w ) + ( M[2][2].re*(*(s+4)).x - M[2][2].im*(*(s+4)).y );
(*(out+4)).y =  ( M[2][0].re*(*(s+3)).y + M[2][0].im*(*(s+3)).x ) + ( M[2][1].re*(*(s+3)).w + M[2][1].im*(*(s+3)).z ) + ( M[2][2].re*(*(s+4)).y + M[2][2].im*(*(s+4)).x );


(*(out+4)).z =  ( M[0][0].re*(*(s+4)).z - M[0][0].im*(*(s+4)).w ) + ( M[0][1].re*(*(s+5)).x - M[0][1].im*(*(s+5)).y ) + ( M[0][2].re*(*(s+5)).z - M[0][2].im*(*(s+5)).w );
(*(out+4)).w =   ( M[0][0].re*(*(s+4)).w + M[0][0].im*(*(s+4)).z ) + ( M[0][1].re*(*(s+5)).y + M[0][1].im*(*(s+5)).x ) + ( M[0][2].re*(*(s+5)).w + M[0][2].im*(*(s+5)).z );


(*(out+5)).x = ( M[1][0].re*(*(s+4)).z - M[1][0].im*(*(s+4)).w ) + ( M[1][1].re*(*(s+5)).x - M[1][1].im*(*(s+5)).y ) + ( M[1][2].re*(*(s+5)).z - M[1][2].im*(*(s+5)).w );
(*(out+5)).y =  ( M[1][0].re*(*(s+4)).w + M[1][0].im*(*(s+4)).z ) + ( M[1][1].re*(*(s+5)).y + M[1][1].im*(*(s+5)).x ) + ( M[1][2].re*(*(s+5)).w + M[1][2].im*(*(s+5)).z );


(*(out+5)).z =  ( M[2][0].re*(*(s+4)).z - M[2][0].im*(*(s+4)).w ) + ( M[2][1].re*(*(s+5)).x - M[2][1].im*(*(s+5)).y ) + ( M[2][2].re*(*(s+5)).z - M[2][2].im*(*(s+5)).w );
(*(out+5)).w =  ( M[2][0].re*(*(s+4)).w + M[2][0].im*(*(s+4)).z ) + ( M[2][1].re*(*(s+5)).y + M[2][1].im*(*(s+5)).x ) + ( M[2][2].re*(*(s+5)).w + M[2][2].im*(*(s+5)).z );
}





//multipliziert gedaggerte su3-Matrix mal Spinor im Dirac-Raum  -- generated with codegen
__device__ void dev_su3MdaggertV(dev_su3 M, dev_spinor * s, dev_spinor * out){
  dev_complex help1;
help1.re = M[0][0].re*(*(s+0)).x + M[0][0].im*(*(s+0)).y + M[1][0].re*(*(s+0)).z + M[1][0].im*(*(s+0)).w + M[2][0].re*(*(s+1)).x + M[2][0].im*(*(s+1)).y;
(*(out+0)).x = help1.re;
help1.im = M[0][0].re*(*(s+0)).y - M[0][0].im*(*(s+0)).x + M[1][0].re*(*(s+0)).w - M[1][0].im*(*(s+0)).z + M[2][0].re*(*(s+1)).y - M[2][0].im*(*(s+1)).x;
(*(out+0)).y = help1.im;

help1.re = M[0][1].re*(*(s+0)).x + M[0][1].im*(*(s+0)).y + M[1][1].re*(*(s+0)).z + M[1][1].im*(*(s+0)).w + M[2][1].re*(*(s+1)).x + M[2][1].im*(*(s+1)).y;
(*(out+0)).z = help1.re;
help1.im = M[0][1].re*(*(s+0)).y - M[0][1].im*(*(s+0)).x + M[1][1].re*(*(s+0)).w - M[1][1].im*(*(s+0)).z + M[2][1].re*(*(s+1)).y - M[2][1].im*(*(s+1)).x;
(*(out+0)).w = help1.im;

help1.re = M[0][2].re*(*(s+0)).x + M[0][2].im*(*(s+0)).y + M[1][2].re*(*(s+0)).z + M[1][2].im*(*(s+0)).w + M[2][2].re*(*(s+1)).x + M[2][2].im*(*(s+1)).y;
(*(out+1)).x = help1.re;
help1.im = M[0][2].re*(*(s+0)).y - M[0][2].im*(*(s+0)).x + M[1][2].re*(*(s+0)).w - M[1][2].im*(*(s+0)).z + M[2][2].re*(*(s+1)).y - M[2][2].im*(*(s+1)).x;
(*(out+1)).y = help1.im;

help1.re = M[0][0].re*(*(s+1)).z + M[0][0].im*(*(s+1)).w + M[1][0].re*(*(s+2)).x + M[1][0].im*(*(s+2)).y + M[2][0].re*(*(s+2)).z + M[2][0].im*(*(s+2)).w;
(*(out+1)).z = help1.re;
help1.im = M[0][0].re*(*(s+1)).w - M[0][0].im*(*(s+1)).z + M[1][0].re*(*(s+2)).y - M[1][0].im*(*(s+2)).x + M[2][0].re*(*(s+2)).w - M[2][0].im*(*(s+2)).z;
(*(out+1)).w = help1.im;

help1.re = M[0][1].re*(*(s+1)).z + M[0][1].im*(*(s+1)).w + M[1][1].re*(*(s+2)).x + M[1][1].im*(*(s+2)).y + M[2][1].re*(*(s+2)).z + M[2][1].im*(*(s+2)).w;
(*(out+2)).x = help1.re;
help1.im = M[0][1].re*(*(s+1)).w - M[0][1].im*(*(s+1)).z + M[1][1].re*(*(s+2)).y - M[1][1].im*(*(s+2)).x + M[2][1].re*(*(s+2)).w - M[2][1].im*(*(s+2)).z;
(*(out+2)).y = help1.im;

help1.re = M[0][2].re*(*(s+1)).z + M[0][2].im*(*(s+1)).w + M[1][2].re*(*(s+2)).x + M[1][2].im*(*(s+2)).y + M[2][2].re*(*(s+2)).z + M[2][2].im*(*(s+2)).w;
(*(out+2)).z = help1.re;
help1.im = M[0][2].re*(*(s+1)).w - M[0][2].im*(*(s+1)).z + M[1][2].re*(*(s+2)).y - M[1][2].im*(*(s+2)).x + M[2][2].re*(*(s+2)).w - M[2][2].im*(*(s+2)).z;
(*(out+2)).w = help1.im;

help1.re = M[0][0].re*(*(s+3)).x + M[0][0].im*(*(s+3)).y + M[1][0].re*(*(s+3)).z + M[1][0].im*(*(s+3)).w + M[2][0].re*(*(s+4)).x + M[2][0].im*(*(s+4)).y;
(*(out+3)).x = help1.re;
help1.im = M[0][0].re*(*(s+3)).y - M[0][0].im*(*(s+3)).x + M[1][0].re*(*(s+3)).w - M[1][0].im*(*(s+3)).z + M[2][0].re*(*(s+4)).y - M[2][0].im*(*(s+4)).x;
(*(out+3)).y = help1.im;

help1.re = M[0][1].re*(*(s+3)).x + M[0][1].im*(*(s+3)).y + M[1][1].re*(*(s+3)).z + M[1][1].im*(*(s+3)).w + M[2][1].re*(*(s+4)).x + M[2][1].im*(*(s+4)).y;
(*(out+3)).z = help1.re;
help1.im = M[0][1].re*(*(s+3)).y - M[0][1].im*(*(s+3)).x + M[1][1].re*(*(s+3)).w - M[1][1].im*(*(s+3)).z + M[2][1].re*(*(s+4)).y - M[2][1].im*(*(s+4)).x;
(*(out+3)).w = help1.im;

help1.re = M[0][2].re*(*(s+3)).x + M[0][2].im*(*(s+3)).y + M[1][2].re*(*(s+3)).z + M[1][2].im*(*(s+3)).w + M[2][2].re*(*(s+4)).x + M[2][2].im*(*(s+4)).y;
(*(out+4)).x = help1.re;
help1.im = M[0][2].re*(*(s+3)).y - M[0][2].im*(*(s+3)).x + M[1][2].re*(*(s+3)).w - M[1][2].im*(*(s+3)).z + M[2][2].re*(*(s+4)).y - M[2][2].im*(*(s+4)).x;
(*(out+4)).y = help1.im;

help1.re = M[0][0].re*(*(s+4)).z + M[0][0].im*(*(s+4)).w + M[1][0].re*(*(s+5)).x + M[1][0].im*(*(s+5)).y + M[2][0].re*(*(s+5)).z + M[2][0].im*(*(s+5)).w;
(*(out+4)).z = help1.re;
help1.im = M[0][0].re*(*(s+4)).w - M[0][0].im*(*(s+4)).z + M[1][0].re*(*(s+5)).y - M[1][0].im*(*(s+5)).x + M[2][0].re*(*(s+5)).w - M[2][0].im*(*(s+5)).z;
(*(out+4)).w = help1.im;

help1.re = M[0][1].re*(*(s+4)).z + M[0][1].im*(*(s+4)).w + M[1][1].re*(*(s+5)).x + M[1][1].im*(*(s+5)).y + M[2][1].re*(*(s+5)).z + M[2][1].im*(*(s+5)).w;
(*(out+5)).x = help1.re;
help1.im = M[0][1].re*(*(s+4)).w - M[0][1].im*(*(s+4)).z + M[1][1].re*(*(s+5)).y - M[1][1].im*(*(s+5)).x + M[2][1].re*(*(s+5)).w - M[2][1].im*(*(s+5)).z;
(*(out+5)).y = help1.im;

help1.re = M[0][2].re*(*(s+4)).z + M[0][2].im*(*(s+4)).w + M[1][2].re*(*(s+5)).x + M[1][2].im*(*(s+5)).y + M[2][2].re*(*(s+5)).z + M[2][2].im*(*(s+5)).w;
(*(out+5)).z = help1.re;
help1.im = M[0][2].re*(*(s+4)).w - M[0][2].im*(*(s+4)).z + M[1][2].re*(*(s+5)).y - M[1][2].im*(*(s+5)).x + M[2][2].re*(*(s+5)).w - M[2][2].im*(*(s+5)).z;
(*(out+5)).w = help1.im;
}




// Gamma t
__device__ void dev_Gamma0(dev_spinor * in){
  REAL tempre,tempim;
     tempre = (*(in+0)).x;
     tempim = (*(in+0)).y;
     (*(in+0)).x = -1.0*(*(in+3)).x;
     (*(in+0)).y = -1.0*(*(in+3)).y;
     (*(in+3)).x = -1.0*tempre;
     (*(in+3)).y = -1.0*tempim;     
     
     tempre = (*(in+0)).z;
     tempim = (*(in+0)).w;
     (*(in+0)).z = -1.0*(*(in+3)).z;
     (*(in+0)).w = -1.0*(*(in+3)).w;
     (*(in+3)).z = -1.0*tempre;
     (*(in+3)).w = -1.0*tempim; 
 
 
     tempre = (*(in+1)).x;
     tempim = (*(in+1)).y;
     (*(in+1)).x = -1.0*(*(in+4)).x;
     (*(in+1)).y = -1.0*(*(in+4)).y;
     (*(in+4)).x = -1.0*tempre;
     (*(in+4)).y = -1.0*tempim;     
     
     tempre = (*(in+1)).z;
     tempim = (*(in+1)).w;
     (*(in+1)).z = -1.0*(*(in+4)).z;
     (*(in+1)).w = -1.0*(*(in+4)).w;
     (*(in+4)).z = -1.0*tempre;
     (*(in+4)).w = -1.0*tempim;     
     
     
     tempre = (*(in+2)).x;
     tempim = (*(in+2)).y;
     (*(in+2)).x = -1.0*(*(in+5)).x;
     (*(in+2)).y = -1.0*(*(in+5)).y;
     (*(in+5)).x = -1.0*tempre;
     (*(in+5)).y = -1.0*tempim;     
   
   
     tempre = (*(in+2)).z;
     tempim = (*(in+2)).w;
     (*(in+2)).z = -1.0*(*(in+5)).z;
     (*(in+2)).w = -1.0*(*(in+5)).w;
     (*(in+5)).z = -1.0*tempre;
     (*(in+5)).w = -1.0*tempim;
}



//Gamma z
__device__ void dev_Gamma3(dev_spinor * in){
  REAL tempre,tempim;
     tempre = (*(in+0)).x;
     tempim = (*(in+0)).y;
     (*(in+0)).x = (*(in+3)).y;
     (*(in+0)).y = -1.0*(*(in+3)).x;
     (*(in+3)).x = -1.0*tempim;
     (*(in+3)).y = tempre;    
     
     tempre = (*(in+0)).z;
     tempim = (*(in+0)).w;
     (*(in+0)).z = (*(in+3)).w;
     (*(in+0)).w = -1.0*(*(in+3)).z;
     (*(in+3)).z = -1.0*tempim;
     (*(in+3)).w = tempre;    
     
     
     tempre = (*(in+1)).x;
     tempim = (*(in+1)).y;
     (*(in+1)).x = (*(in+4)).y;
     (*(in+1)).y = -1.0*(*(in+4)).x;
     (*(in+4)).x = -1.0*tempim;
     (*(in+4)).y = tempre;     
     
     
     tempre = (*(in+1)).z;
     tempim = (*(in+1)).w;
     (*(in+1)).z = -1.0*(*(in+4)).w;
     (*(in+1)).w = (*(in+4)).z;
     (*(in+4)).z  = tempim;
     (*(in+4)).w  = -1.0*tempre;     
     
     
     
     tempre = (*(in+2)).x;
     tempim = (*(in+2)).y;
     (*(in+2)).x = -1.0*(*(in+5)).y;
     (*(in+2)).y = (*(in+5)).x;
     (*(in+5)).x = tempim;
     (*(in+5)).y = -1.0*tempre;    
     
     
     tempre = (*(in+2)).z;
     tempim = (*(in+2)).w;
     (*(in+2)).z = -1.0*(*(in+5)).w;
     (*(in+2)).w = (*(in+5)).z;
     (*(in+5)).z = tempim;
     (*(in+5)).w = -1.0*tempre;

}



//Gamma y
__device__ void dev_Gamma2(dev_spinor * in){
  REAL tempre,tempim;
     tempre = (*(in+0)).x;
     tempim = (*(in+0)).y;
     (*(in+0)).x = -1.0*(*(in+4)).z;
     (*(in+0)).y = -1.0*(*(in+4)).w;
     (*(in+4)).z = -1.0*tempre;
     (*(in+4)).w = -1.0*tempim;    
     
     tempre = (*(in+0)).z;
     tempim = (*(in+0)).w;
     (*(in+0)).z = -1.0*(*(in+5)).x;
     (*(in+0)).w = -1.0*(*(in+5)).y;
     (*(in+5)).x = -1.0*tempre;
     (*(in+5)).y = -1.0*tempim;     
     
     
     tempre = (*(in+1)).x;
     tempim = (*(in+1)).y;
     (*(in+1)).x = -1.0*(*(in+5)).z;
     (*(in+1)).y = -1.0*(*(in+5)).w;
     (*(in+5)).z = -1.0*tempre;
     (*(in+5)).w = -1.0*tempim;     
     
     tempre = (*(in+1)).z;
     tempim = (*(in+1)).w;
     (*(in+1)).z = (*(in+3)).x;
     (*(in+1)).w = (*(in+3)).y;
     (*(in+3)).x = tempre;
     (*(in+3)).y = tempim;     
     
     tempre = (*(in+2)).x;
     tempim = (*(in+2)).y;
     (*(in+2)).x = (*(in+3)).z;
     (*(in+2)).y = (*(in+3)).w;
     (*(in+3)).z = tempre;
     (*(in+3)).w = tempim;     
     
     
     tempre = (*(in+2)).z;
     tempim = (*(in+2)).w;
     (*(in+2)).z = (*(in+4)).x;
     (*(in+2)).w = (*(in+4)).y;
     (*(in+4)).x = tempre;
     (*(in+4)).y = tempim;
}



//Gamma x
__device__ void dev_Gamma1(dev_spinor * in){
  REAL tempre,tempim;
     tempre = (*(in+0)).x;
     tempim = (*(in+0)).y;
     (*(in+0)).x = (*(in+4)).w;
     (*(in+0)).y = -1.0*(*(in+4)).z;
     (*(in+4)).z  = -1.0*tempim;
     (*(in+4)).w  = tempre;    
     
     tempre = (*(in+0)).z;
     tempim = (*(in+0)).w;
     (*(in+0)).z = (*(in+5)).y;
     (*(in+0)).w = -1.0*(*(in+5)).x;
     (*(in+5)).x = -1.0*tempim;
     (*(in+5)).y = tempre;     
     
     tempre = (*(in+1)).x;
     tempim = (*(in+1)).y;
     (*(in+1)).x = (*(in+5)).w;
     (*(in+1)).y = -1.0*(*(in+5)).z;
     (*(in+5)).z = -1.0*tempim;
     (*(in+5)).w = tempre;     
     
     tempre = (*(in+1)).z;
     tempim = (*(in+1)).w;
     (*(in+1)).z = (*(in+3)).y;
     (*(in+1)).w = -1.0*(*(in+3)).x;
     (*(in+3)).x = -1.0*tempim;
     (*(in+3)).y = tempre;     
     
     tempre = (*(in+2)).x;
     tempim = (*(in+2)).y;
     (*(in+2)).x = (*(in+3)).w;
     (*(in+2)).y = -1.0*(*(in+3)).z;
     (*(in+3)).z = -1.0*tempim;
     (*(in+3)).w = tempre;     
     
     
     tempre = (*(in+2)).z;
     tempim = (*(in+2)).w;
     (*(in+2)).z = (*(in+4)).y;
     (*(in+2)).w = -1.0*(*(in+4)).x;
     (*(in+4)).x = -1.0*tempim;
     (*(in+4)).y = tempre;
  
}



__device__ void dev_Gamma5(dev_spinor * in){
          (*(in+3)).x = -1.0*(*(in+3)).x;
          (*(in+3)).y = -1.0*(*(in+3)).y;
          (*(in+3)).z = -1.0*(*(in+3)).z;
          (*(in+3)).w = -1.0*(*(in+3)).w;
          (*(in+4)).x = -1.0*(*(in+4)).x;
          (*(in+4)).y = -1.0*(*(in+4)).y; 

          (*(in+4)).z = -1.0*(*(in+4)).z;
          (*(in+4)).w = -1.0*(*(in+4)).w;
          (*(in+5)).x = -1.0*(*(in+5)).x;
          (*(in+5)).y = -1.0*(*(in+5)).y;
          (*(in+5)).z = -1.0*(*(in+5)).z;
          (*(in+5)).w = -1.0*(*(in+5)).w;  
}





// older version, all in one function
__device__ void dev_GammatV(int mu, dev_spinor * in){//multipliziert Gamma(mu)*V effizientes ausnutzen der Nullen 
 REAL tempre,tempim;
 /* ORDER: t, z, y, x*/
 switch (mu){
 
 case 0:
     tempre = (*(in+0)).x;
     tempim = (*(in+0)).y;
     (*(in+0)).x = -1.0*(*(in+3)).x;
     (*(in+0)).y = -1.0*(*(in+3)).y;
     (*(in+3)).x = -1.0*tempre;
     (*(in+3)).y = -1.0*tempim;     
     
     tempre = (*(in+0)).z;
     tempim = (*(in+0)).w;
     (*(in+0)).z = -1.0*(*(in+3)).z;
     (*(in+0)).w = -1.0*(*(in+3)).w;
     (*(in+3)).z = -1.0*tempre;
     (*(in+3)).w = -1.0*tempim; 
 
 
     tempre = (*(in+1)).x;
     tempim = (*(in+1)).y;
     (*(in+1)).x = -1.0*(*(in+4)).x;
     (*(in+1)).y = -1.0*(*(in+4)).y;
     (*(in+4)).x = -1.0*tempre;
     (*(in+4)).y = -1.0*tempim;     
     
     tempre = (*(in+1)).z;
     tempim = (*(in+1)).w;
     (*(in+1)).z = -1.0*(*(in+4)).z;
     (*(in+1)).w = -1.0*(*(in+4)).w;
     (*(in+4)).z = -1.0*tempre;
     (*(in+4)).w = -1.0*tempim;     
     
     
     tempre = (*(in+2)).x;
     tempim = (*(in+2)).y;
     (*(in+2)).x = -1.0*(*(in+5)).x;
     (*(in+2)).y = -1.0*(*(in+5)).y;
     (*(in+5)).x = -1.0*tempre;
     (*(in+5)).y = -1.0*tempim;     
   
   
     tempre = (*(in+2)).z;
     tempim = (*(in+2)).w;
     (*(in+2)).z = -1.0*(*(in+5)).z;
     (*(in+2)).w = -1.0*(*(in+5)).w;
     (*(in+5)).z = -1.0*tempre;
     (*(in+5)).w = -1.0*tempim;

 break;
 
 case 1:
     tempre = (*(in+0)).x;
     tempim = (*(in+0)).y;
     (*(in+0)).x = (*(in+3)).y;
     (*(in+0)).y = -1.0*(*(in+3)).x;
     (*(in+3)).x = -1.0*tempim;
     (*(in+3)).y = tempre;    
     
     tempre = (*(in+0)).z;
     tempim = (*(in+0)).w;
     (*(in+0)).z = (*(in+3)).w;
     (*(in+0)).w = -1.0*(*(in+3)).z;
     (*(in+3)).z = -1.0*tempim;
     (*(in+3)).w = tempre;    
     
     
     tempre = (*(in+1)).x;
     tempim = (*(in+1)).y;
     (*(in+1)).x = (*(in+4)).y;
     (*(in+1)).y = -1.0*(*(in+4)).x;
     (*(in+4)).x = -1.0*tempim;
     (*(in+4)).y = tempre;     
     
     
     tempre = (*(in+1)).z;
     tempim = (*(in+1)).w;
     (*(in+1)).z = -1.0*(*(in+4)).w;
     (*(in+1)).w = (*(in+4)).z;
     (*(in+4)).z  = tempim;
     (*(in+4)).w  = -1.0*tempre;     
     
     
     
     tempre = (*(in+2)).x;
     tempim = (*(in+2)).y;
     (*(in+2)).x = -1.0*(*(in+5)).y;
     (*(in+2)).y = (*(in+5)).x;
     (*(in+5)).x = tempim;
     (*(in+5)).y = -1.0*tempre;    
     
     
     tempre = (*(in+2)).z;
     tempim = (*(in+2)).w;
     (*(in+2)).z = -1.0*(*(in+5)).w;
     (*(in+2)).w = (*(in+5)).z;
     (*(in+5)).z = tempim;
     (*(in+5)).w = -1.0*tempre;


 break;
 
 case 2:
     tempre = (*(in+0)).x;
     tempim = (*(in+0)).y;
     (*(in+0)).x = -1.0*(*(in+4)).z;
     (*(in+0)).y = -1.0*(*(in+4)).w;
     (*(in+4)).z = -1.0*tempre;
     (*(in+4)).w = -1.0*tempim;    
     
     tempre = (*(in+0)).z;
     tempim = (*(in+0)).w;
     (*(in+0)).z = -1.0*(*(in+5)).x;
     (*(in+0)).w = -1.0*(*(in+5)).y;
     (*(in+5)).x = -1.0*tempre;
     (*(in+5)).y = -1.0*tempim;     
     
     
     tempre = (*(in+1)).x;
     tempim = (*(in+1)).y;
     (*(in+1)).x = -1.0*(*(in+5)).z;
     (*(in+1)).y = -1.0*(*(in+5)).w;
     (*(in+5)).z = -1.0*tempre;
     (*(in+5)).w = -1.0*tempim;     
     
     tempre = (*(in+1)).z;
     tempim = (*(in+1)).w;
     (*(in+1)).z = (*(in+3)).x;
     (*(in+1)).w = (*(in+3)).y;
     (*(in+3)).x = tempre;
     (*(in+3)).y = tempim;     
     
     tempre = (*(in+2)).x;
     tempim = (*(in+2)).y;
     (*(in+2)).x = (*(in+3)).z;
     (*(in+2)).y = (*(in+3)).w;
     (*(in+3)).z = tempre;
     (*(in+3)).w = tempim;     
     
     
     tempre = (*(in+2)).z;
     tempim = (*(in+2)).w;
     (*(in+2)).z = (*(in+4)).x;
     (*(in+2)).w = (*(in+4)).y;
     (*(in+4)).x = tempre;
     (*(in+4)).y = tempim;

 break; 
 
 case 3:


     tempre = (*(in+0)).x;
     tempim = (*(in+0)).y;
     (*(in+0)).x = (*(in+4)).w;
     (*(in+0)).y = -1.0*(*(in+4)).z;
     (*(in+4)).z  = -1.0*tempim;
     (*(in+4)).w  = tempre;    
     
     tempre = (*(in+0)).z;
     tempim = (*(in+0)).w;
     (*(in+0)).z = (*(in+5)).y;
     (*(in+0)).w = -1.0*(*(in+5)).x;
     (*(in+5)).x = -1.0*tempim;
     (*(in+5)).y = tempre;     
     
     tempre = (*(in+1)).x;
     tempim = (*(in+1)).y;
     (*(in+1)).x = (*(in+5)).w;
     (*(in+1)).y = -1.0*(*(in+5)).z;
     (*(in+5)).z = -1.0*tempim;
     (*(in+5)).w = tempre;     
     
     tempre = (*(in+1)).z;
     tempim = (*(in+1)).w;
     (*(in+1)).z = (*(in+3)).y;
     (*(in+1)).w = -1.0*(*(in+3)).x;
     (*(in+3)).x = -1.0*tempim;
     (*(in+3)).y = tempre;     
     
     tempre = (*(in+2)).x;
     tempim = (*(in+2)).y;
     (*(in+2)).x = (*(in+3)).w;
     (*(in+2)).y = -1.0*(*(in+3)).z;
     (*(in+3)).z = -1.0*tempim;
     (*(in+3)).w = tempre;     
     
     
     tempre = (*(in+2)).z;
     tempim = (*(in+2)).w;
     (*(in+2)).z = (*(in+4)).y;
     (*(in+2)).w = -1.0*(*(in+4)).x;
     (*(in+4)).x = -1.0*tempim;
     (*(in+4)).y = tempre;
     
     
 break;
 
 
 case 4:
  
          (*(in+3)).x = -1.0*(*(in+3)).x;
          (*(in+3)).y = -1.0*(*(in+3)).y;
          (*(in+3)).z = -1.0*(*(in+3)).z;
          (*(in+3)).w = -1.0*(*(in+3)).w;
          (*(in+4)).x = -1.0*(*(in+4)).x;
          (*(in+4)).y = -1.0*(*(in+4)).y; 

          (*(in+4)).z = -1.0*(*(in+4)).z;
          (*(in+4)).w = -1.0*(*(in+4)).w;
          (*(in+5)).x = -1.0*(*(in+5)).x;
          (*(in+5)).y = -1.0*(*(in+5)).y;
          (*(in+5)).z = -1.0*(*(in+5)).z;
          (*(in+5)).w = -1.0*(*(in+5)).w;  
 break;
 }
}




// reconstruction of the link fields from two rows of the su3 matrix
// numbers are fetched from texture cache
__device__ void dev_reconstructgf_2vtexref (int pos, dev_su3* gf){
  dev_complex help1;
  dev_complex help2;
  float4 gfin;
  
  gfin = tex1Dfetch(gf_tex,3*pos);
  //first row
  (*gf)[0][0].re = gfin.x;
  (*gf)[0][0].im = gfin.y;
  (*gf)[0][1].re = gfin.z;
  (*gf)[0][1].im = gfin.w;
  gfin = tex1Dfetch(gf_tex,3*pos+1);
  (*gf)[0][2].re = gfin.x;
  (*gf)[0][2].im = gfin.y;
  //second row
  (*gf)[1][0].re = gfin.z;
  (*gf)[1][0].im = gfin.w;
  gfin = tex1Dfetch(gf_tex,3*pos+2);
  (*gf)[1][1].re = gfin.x;
  (*gf)[1][1].im = gfin.y;
  (*gf)[1][2].re = gfin.z;
  (*gf)[1][2].im = gfin.w;
  
  //third row from cconj(cross product of first and second row)
  help1 = dev_cmult((*gf)[0][1],(*gf)[1][2]);
  help2 = dev_cmult((*gf)[0][2],(*gf)[1][1]);
  help1 = dev_cconj(dev_csub(help1,help2));
  (*gf)[2][0] = help1;

  
  help1 = dev_cmult((*gf)[0][2],(*gf)[1][0]);
  help2 = dev_cmult((*gf)[0][0],(*gf)[1][2]);
  help1 = dev_cconj(dev_csub(help1,help2));
  
  (*gf)[2][1] = help1;
  
  help1 = dev_cmult((*gf)[0][0],(*gf)[1][1]);
  help2 = dev_cmult((*gf)[0][1],(*gf)[1][0]);
  help1 = dev_cconj(dev_csub(help1,help2));
  
  (*gf)[2][2] = help1;
  
  return;
}




// su3 - dagger reconstruction from two rows  
__device__ void dev_reconstructgf_2vtexref_dagger (int pos, dev_su3* gf){
  dev_complex help1;
  dev_complex help2;
  float4 gfin;
  
  
  //first column (minus in im for complex conj.)
  gfin = tex1Dfetch(gf_tex,3*pos);
  (*gf)[0][0].re = gfin.x;
  (*gf)[0][0].im = -gfin.y;
  (*gf)[1][0].re = gfin.z;
  (*gf)[1][0].im = -gfin.w;
  gfin = tex1Dfetch(gf_tex,3*pos+1);
  (*gf)[2][0].re = gfin.x;
  (*gf)[2][0].im = -gfin.y;
  
  //second  column (minus in im for complex conj.)
  (*gf)[0][1].re = gfin.z;
  (*gf)[0][1].im = -gfin.w;
  gfin = tex1Dfetch(gf_tex,3*pos+2);
  (*gf)[1][1].re = gfin.x;
  (*gf)[1][1].im = -gfin.y;
  (*gf)[2][1].re = gfin.z;
  (*gf)[2][1].im = -gfin.w;
  
  //third column from (cross product of cconj(first column) and cconj(second column))
 
  help1 = dev_cconj(dev_cmult((*gf)[1][0],(*gf)[2][1]));
  help2 = dev_cconj(dev_cmult((*gf)[2][0],(*gf)[1][1]));
  help1 = dev_csub(help1,help2);
  (*gf)[0][2] = help1;

  
  help1 = dev_cconj(dev_cmult((*gf)[2][0],(*gf)[0][1]));
  help2 = dev_cconj(dev_cmult((*gf)[0][0],(*gf)[2][1]));
  help1 = dev_csub(help1,help2);
  (*gf)[1][2] = help1;

  
  help1 = dev_cconj(dev_cmult((*gf)[0][0],(*gf)[1][1]));
  help2 = dev_cconj(dev_cmult((*gf)[1][0],(*gf)[0][1]));
  help1 = dev_csub(help1,help2);
  (*gf)[2][2] = help1;
  
  
  /* does this also work?
  help1 = dev_cmult((*gf)[1][0],(*gf)[2][1]);
  help2 = dev_cmult((*gf)[2][0],(*gf)[1][1]);
  help1 = dev_cconj(dev_csub(help1,help2));
  (*gf)[0][2] = help1;

  
  help1 = dev_cmult((*gf)[2][0],(*gf)[0][1]);
  help2 = dev_cmult((*gf)[0][0],(*gf)[2][1]);
  help1 = dev_cconj(dev_csub(help1,help2));
  (*gf)[1][2] = help1;

  
  help1 = dev_cmult((*gf)[0][0],(*gf)[1][1]);
  help2 = dev_cmult((*gf)[1][0],(*gf)[0][1]));
  help1 = dev_cconj(dev_csub(help1,help2));
  (*gf)[2][2] = help1;
  */
}








__device__ void dev_reconstructgf_8texref (int pos, dev_su3* gf){

  float4 gfin;
  REAL N, one_over_N, help;
  dev_complex p1,p2, chelp1, chelp2;
  
  gfin = tex1Dfetch(gf_tex,2*pos);
  // read a2 a3
  (*gf)[0][1].re = gfin.x;
  (*gf)[0][1].im = gfin.y;
  (*gf)[0][2].re = gfin.z;
  (*gf)[0][2].im = gfin.w;  
 
  help = gfin.x*gfin.x + gfin.y*gfin.y + gfin.z*gfin.z + gfin.w*gfin.w; // use later on
  N = sqrtf(help);
  one_over_N = 1.0/N;
  
  // read theta_a1, theta_c1, b1
  gfin = tex1Dfetch(gf_tex,2*pos + 1);
  
  // reconstruct a1
  help = sqrtf(1.0f - help);
  (*gf)[0][0].re = help*cosf(gfin.x);
  (*gf)[0][0].im = help*sinf(gfin.x);
  
  // assign b1
  (*gf)[1][0].re = gfin.z;
  (*gf)[1][0].im = gfin.w;
  
  // p2 = 1/N b1
  p2.re = one_over_N*(*gf)[1][0].re;
  p2.im = one_over_N*(*gf)[1][0].im;  


  // reconstruct c1
  help = sqrtf(1.0 - 
              (*gf)[0][0].re * (*gf)[0][0].re - (*gf)[0][0].im * (*gf)[0][0].im - 
              (*gf)[1][0].re * (*gf)[1][0].re - (*gf)[1][0].im * (*gf)[1][0].im
          );
  (*gf)[2][0].re = help*cosf(gfin.y);
  (*gf)[2][0].im = help*sinf(gfin.y);

  
  // p1 = 1/N*cconj(c1)
  p1.re = one_over_N*(*gf)[2][0].re;
  p1.im = - one_over_N*(*gf)[2][0].im;
  
  
  // calculate b2
  chelp1 = dev_cmult(p1,  dev_cconj( (*gf)[0][2] )   );
  chelp2 = dev_cmult(p2, dev_cmult( dev_cconj((*gf)[0][0]) , (*gf)[0][1] )  );
  chelp1 = dev_cadd(chelp1, chelp2);
  (*gf)[1][1] = dev_crealmult(chelp1, -one_over_N);
  
  // calculate b3
  chelp1 = dev_cmult(p1,  dev_cconj( (*gf)[0][1] )   );
  chelp2 = dev_cmult(p2, dev_cmult( dev_cconj((*gf)[0][0]) , (*gf)[0][2] )  );
  chelp1 = dev_cadd(chelp1, chelp2);
  (*gf)[1][2] = dev_crealmult(chelp1, -one_over_N);
  
  // calculate c2
  chelp1 = dev_cmult(  dev_cconj(p2) ,  dev_cconj( (*gf)[0][2] )   );
  chelp2 = dev_cmult(  dev_cconj(p1) , 
                       dev_cmult(  dev_cconj( (*gf)[0][0] ) , (*gf)[0][1] )
                     );
  chelp1 = dev_csub(chelp1, chelp2);
  (*gf)[2][1] = dev_crealmult(chelp1, one_over_N);
  
  
  // calculate c3
  chelp1 = dev_cmult(  dev_cconj(p2) ,  dev_cconj( (*gf)[0][1] )   );
  chelp2 = dev_cmult(  dev_cconj(p1) , 
                       dev_cmult(  dev_cconj( (*gf)[0][0] ) , (*gf)[0][2] )
                     );
  chelp1 = dev_csub(chelp1, chelp2);
  (*gf)[2][2] = dev_crealmult(chelp1, one_over_N);
                       
}



__device__ void dev_reconstructgf_8texref_dagger (int pos, dev_su3* gf){


  float4 gfin;
  REAL N, one_over_N, help;
  dev_complex p1,p2, chelp1, chelp2;
  
  gfin = tex1Dfetch(gf_tex,2*pos);
  // read a2 a3
  (*gf)[0][1].re = gfin.x;
  (*gf)[0][1].im = gfin.y;
  (*gf)[0][2].re = gfin.z;
  (*gf)[0][2].im = gfin.w;  
 
  help = gfin.x*gfin.x + gfin.y*gfin.y + gfin.z*gfin.z + gfin.w*gfin.w; // use later on
  N = sqrtf(help);
  one_over_N = 1.0/N;
  
  // read theta_a1, theta_c1, b1
  gfin = tex1Dfetch(gf_tex,2*pos + 1);
  
  // reconstruct a1
  help = sqrtf(1.0f - help);
  (*gf)[0][0].re = help*cosf(gfin.x);
  (*gf)[0][0].im = help*sinf(gfin.x);
  
  // assign b1
  (*gf)[1][0].re = gfin.z;
  (*gf)[1][0].im = gfin.w;
  
  // p2 = 1/N b1
  p2.re = one_over_N*(*gf)[1][0].re;
  p2.im = one_over_N*(*gf)[1][0].im;  


  // reconstruct c1
  help = sqrtf(1.0 - 
              (*gf)[0][0].re * (*gf)[0][0].re - (*gf)[0][0].im * (*gf)[0][0].im - 
              (*gf)[1][0].re * (*gf)[1][0].re - (*gf)[1][0].im * (*gf)[1][0].im
          );
  (*gf)[2][0].re = help*cosf(gfin.y);
  (*gf)[2][0].im = help*sinf(gfin.y);

  
  // p1 = 1/N*cconj(c1)
  p1.re = one_over_N*(*gf)[2][0].re;
  p1.im = - one_over_N*(*gf)[2][0].im;
  
  
  // calculate b2
  chelp1 = dev_cmult(p1,  dev_cconj( (*gf)[0][2] )   );
  chelp2 = dev_cmult(p2, dev_cmult( dev_cconj((*gf)[0][0]) , (*gf)[0][1] )  );
  chelp1 = dev_cadd(chelp1, chelp2);
  (*gf)[1][1] = dev_crealmult(chelp1, -one_over_N);
  
  // calculate b3
  chelp1 = dev_cmult(p1,  dev_cconj( (*gf)[0][1] )   );
  chelp2 = dev_cmult(p2, dev_cmult( dev_cconj((*gf)[0][0]) , (*gf)[0][2] )  );
  chelp1 = dev_cadd(chelp1, chelp2);
  (*gf)[1][2] = dev_crealmult(chelp1, -one_over_N);
  
  // calculate c2
  chelp1 = dev_cmult(  dev_cconj(p2) ,  dev_cconj( (*gf)[0][2] )   );
  chelp2 = dev_cmult(  dev_cconj(p1) , 
                       dev_cmult(  dev_cconj( (*gf)[0][0] ) , (*gf)[0][1] )
                     );
  chelp1 = dev_csub(chelp1, chelp2);
  (*gf)[2][1] = dev_crealmult(chelp1, one_over_N);
  
  
  // calculate c3
  chelp1 = dev_cmult(  dev_cconj(p2) ,  dev_cconj( (*gf)[0][1] )   );
  chelp2 = dev_cmult(  dev_cconj(p1) , 
                       dev_cmult(  dev_cconj( (*gf)[0][0] ) , (*gf)[0][2] )
                     );
  chelp1 = dev_csub(chelp1, chelp2);
  (*gf)[2][2] = dev_crealmult(chelp1, one_over_N);




  // till here this was the reconstruction of the link (not daggered)
  // calculate now U^dagger trivially
  // this should be implemented more efficiently !!
  
  
  //diagonal
  (*gf)[0][0] = dev_cconj((*gf)[0][0]);
  (*gf)[1][1] = dev_cconj((*gf)[1][1]);
  (*gf)[2][2] = dev_cconj((*gf)[2][2]);
  
  //off-diagonals
  chelp1 = dev_cconj(  (*gf)[0][1] );
  (*gf)[0][1] = dev_cconj ( (*gf)[1][0]  );
  (*gf)[1][0] = chelp1;
  
  chelp1 = dev_cconj(  (*gf)[0][2] );
  (*gf)[0][2] = dev_cconj ( (*gf)[2][0]  );
  (*gf)[2][0] = chelp1;
  
  chelp1 = dev_cconj(  (*gf)[1][0] );
  (*gf)[1][0] = dev_cconj ( (*gf)[0][1]  );
  (*gf)[0][1] = chelp1;
  
}







__global__ void dev_gamma5(dev_spinor * sin, dev_spinor * sout){
  int pos;
  pos= threadIdx.x + blockDim.x*blockIdx.x;  
  if(pos < dev_VOLUME){
          sout[6*pos+0].x = sin[6*pos+0].x;
          sout[6*pos+0].y = sin[6*pos+0].y;
          sout[6*pos+0].z = sin[6*pos+0].z;
          sout[6*pos+0].w = sin[6*pos+0].w;
          sout[6*pos+1].x = sin[6*pos+1].x;
          sout[6*pos+1].y = sin[6*pos+1].y;
          
          sout[6*pos+1].z = sin[6*pos+1].z;
          sout[6*pos+1].w = sin[6*pos+1].w;
          sout[6*pos+2].x = sin[6*pos+2].x;
          sout[6*pos+2].y = sin[6*pos+2].y;
          sout[6*pos+2].z = sin[6*pos+2].z;
          sout[6*pos+2].w = sin[6*pos+2].w;   
          
          sout[6*pos+3].x = -1.0*sin[6*pos+3].x;
          sout[6*pos+3].y = -1.0*sin[6*pos+3].y;
          sout[6*pos+3].z = -1.0*sin[6*pos+3].z;
          sout[6*pos+3].w = -1.0*sin[6*pos+3].w;
          sout[6*pos+4].x = -1.0*sin[6*pos+4].x;
          sout[6*pos+4].y = -1.0*sin[6*pos+4].y; 

          sout[6*pos+4].z = -1.0*sin[6*pos+4].z;
          sout[6*pos+4].w = -1.0*sin[6*pos+4].w;
          sout[6*pos+5].x = -1.0*sin[6*pos+5].x;
          sout[6*pos+5].y = -1.0*sin[6*pos+5].y;
          sout[6*pos+5].z = -1.0*sin[6*pos+5].z;
          sout[6*pos+5].w = -1.0*sin[6*pos+5].w;                 
  } 
}




__global__ void dev_swapmu(){
  if(blockIdx.x == 0 && threadIdx.x == 0){
    mu = - mu;
  }
}


// computes sout = 1/(1 +- mutilde gamma5) sin = (1 -+ i mutilde gamma5)/(1+mutilde^2) sin
// mutilde = 2 kappa mu
__global__ void dev_mul_one_pm_imu_inv(dev_spinor* sin, dev_spinor* sout, const REAL sign){
  
   //need the inverse sign in the numerator because of inverse
   dev_complex pm_imu = dev_initcomplex(0.0,-1.0*sign*twokappamu);
   
   REAL one_plus_musquare_inv = 1.0/(1.0 + twokappamu*twokappamu);
   int pos;
   pos= threadIdx.x + blockDim.x*blockIdx.x;  
   if(pos < dev_VOLUME){
     dev_skalarmult_spinor(&(sin[6*pos]), pm_imu, &(sout[6*pos]));
     dev_Gamma5(&(sout[6*pos]));
     dev_add_spinor_assign(&(sout[6*pos]), &(sin[6*pos]));
     dev_realmult_spinor(&(sout[6*pos]), one_plus_musquare_inv);
   }
}

// sout = gamma_5*((1\pm i\mutilde \gamma_5)*sin1 - sin2)
__global__ void dev_mul_one_pm_imu_sub_mul_gamma5(dev_spinor* sin1, dev_spinor* sin2, dev_spinor* sout, const REAL sign){

   dev_complex pm_imu = dev_initcomplex(0.0, sign*twokappamu); // i mutilde
   int pos;
   pos= threadIdx.x + blockDim.x*blockIdx.x;  
   if(pos < dev_VOLUME){
     dev_skalarmult_spinor(&(sin1[6*pos]), pm_imu, &(sout[6*pos]));
     dev_Gamma5(&(sout[6*pos]));
     dev_add_spinor_assign(&(sout[6*pos]), &(sin1[6*pos]));
     dev_sub_spinor_assign(&(sout[6*pos]), &(sin2[6*pos]));
     dev_Gamma5(&(sout[6*pos]));
   }   
}




//applies the Hopping Part Even-Odd !
//the gauge field is the complete gaugefield!
//the gauge field at the local point is reconstructed by 2*pos+eo where pos is the eo-position
//from 0..VOLUME/2-1, eo = 0 or 1
//the positions in the gauge fields are passed in "gfindex_site" for gf's that are attached at
//the actual positions and in "gfindex_nextsite" for gf's that start at a position of the 
//other eo-sublattice.
//for the hopping positions of the eo-spinor field we use on of the two dedicated eo-nn fields
//the boundary conditions are implemented as in Hopping_Matrix.c
//mult with complex conjugate k0,k1,k2,k3 in positive direction because
// psi(x+mu) != exp(i theta_mu) psi(x)  
__global__ void dev_Hopping_Matrix(dev_su3_2v * gf, dev_spinor * sin, dev_spinor * sout, int * gfindex_site,int* gfindex_nextsite, int * nn_evenodd, const int eo){

  int pos,hoppos;
    dev_spinor shelp1[6], ssum[6];
    __shared__ dev_su3 gfsmem[BLOCK];
    

  pos= threadIdx.x + blockDim.x*blockIdx.x;  
  int ix = threadIdx.x;
  if(pos < dev_VOLUME){

  dev_zero_spinor(&(ssum[0])); // zero sum        
//hopping term                
//l==0,t
            //positive direction
            hoppos = nn_evenodd[8*pos];
            //color
            dev_reconstructgf_2vtexref(4*(gfindex_site[pos]),&(gfsmem[ix]));
            dev_su3MtV_spintex(gfsmem[ix], hoppos, &(shelp1[0]));
            //-kappa(r - gamma_mu)
            dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_mk0,&(shelp1[0]), &(ssum[0]));
            dev_Gamma0(&(shelp1[0]));
            dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_k0,&(shelp1[0]), &(ssum[0]));

            //negative direction
            hoppos = nn_evenodd[8*pos+4]; 
            //color
            dev_reconstructgf_2vtexref_dagger(4*gfindex_nextsite[hoppos],&(gfsmem[ix]));
            dev_su3MtV_spintex(gfsmem[ix], hoppos, &(shelp1[0]));     
            //-kappa(r + gamma_mu)
            dev_complexmult_add_assign_spinor(&(ssum[0]),dev_mk0,&(shelp1[0]), &(ssum[0]));
            dev_Gamma0(&(shelp1[0]));
            dev_complexmult_add_assign_spinor(&(ssum[0]),dev_mk0,&(shelp1[0]), &(ssum[0]));


//l==3,z               
            //positive direction
            hoppos = nn_evenodd[8*pos+3];
            //color
            dev_reconstructgf_2vtexref(4*(gfindex_site[pos])+(3),&(gfsmem[ix]));
            dev_su3MtV_spintex(gfsmem[ix], hoppos, &(shelp1[0]));
            //-kappa(r - gamma_mu)
            dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_mk3,&(shelp1[0]), &(ssum[0]));
            dev_Gamma3(&(shelp1[0]));
            dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_k3,&(shelp1[0]), &(ssum[0]));

            //negative direction
            hoppos = nn_evenodd[8*pos+7]; 
            //color
            dev_reconstructgf_2vtexref_dagger(4*gfindex_nextsite[hoppos]+(3),&(gfsmem[ix]));
            dev_su3MtV_spintex(gfsmem[ix], hoppos, &(shelp1[0]));
            //-kappa(r + gamma_mu)
            dev_complexmult_add_assign_spinor(&(ssum[0]),dev_mk3,&(shelp1[0]), &(ssum[0]));
            dev_Gamma3(&(shelp1[0]));
            dev_complexmult_add_assign_spinor(&(ssum[0]),dev_mk3,&(shelp1[0]), &(ssum[0]));
         
         
//l==2,y        
            //positive direction
            hoppos = nn_evenodd[8*pos+2];
            //color
            dev_reconstructgf_2vtexref(4*(gfindex_site[pos])+(2),&(gfsmem[ix]));
            dev_su3MtV_spintex(gfsmem[ix], hoppos, &(shelp1[0]));
            //-kappa(r - gamma_mu)
            dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_mk2,&(shelp1[0]), &(ssum[0]));
            dev_Gamma2(&(shelp1[0]));
            dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_k2,&(shelp1[0]), &(ssum[0]));
            
            //negative direction
            hoppos = nn_evenodd[8*pos+6]; 
            //color
            dev_reconstructgf_2vtexref_dagger(4*gfindex_nextsite[hoppos]+(2),&(gfsmem[ix]));
            dev_su3MtV_spintex(gfsmem[ix], hoppos, &(shelp1[0]));
            //-kappa(r + gamma_mu)
            dev_complexmult_add_assign_spinor(&(ssum[0]),dev_mk2,&(shelp1[0]), &(ssum[0]));
            dev_Gamma2(&(shelp1[0]));
            dev_complexmult_add_assign_spinor(&(ssum[0]),dev_mk2,&(shelp1[0]), &(ssum[0]));


//l==1,x 
            //positive direction
            hoppos = nn_evenodd[8*pos+1];
            //color
            dev_reconstructgf_2vtexref(4*(gfindex_site[pos])+(1),&(gfsmem[ix]));
            dev_su3MtV_spintex(gfsmem[ix], hoppos, &(shelp1[0]));
            //-kappa(r - gamma_mu)
            dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_mk1,&(shelp1[0]), &(ssum[0]));
            dev_Gamma1(&(shelp1[0]));
            dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_k1,&(shelp1[0]), &(ssum[0]));

            //negative direction
            hoppos = nn_evenodd[8*pos+5]; 
            //color
            dev_reconstructgf_2vtexref_dagger(4*gfindex_nextsite[hoppos]+(1),&(gfsmem[ix]));
            dev_su3MtV_spintex(gfsmem[ix], hoppos, &(shelp1[0]));
            //-kappa(r + gamma_mu)
            dev_complexmult_add_assign_spinor(&(ssum[0]),dev_mk1,&(shelp1[0]), &(ssum[0]));
            dev_Gamma1(&(shelp1[0]));
            dev_complexmult_add_assign_spinor(&(ssum[0]),dev_mk1,&(shelp1[0]), &(ssum[0]));      
               
 
        //copy to output spinor
        dev_copy_spinor(&(ssum[0]),&(sout[6*pos])); 
  }
}





// aequivalent to Qtm_pm_psi in tm_operators.c
extern "C" void dev_Qtm_pm_psi(dev_spinor* spinin, dev_spinor* spinout, int gridsize, int blocksize){
  //spinin == odd
  //spinout == odd
  
  //Q_{-}
  bind_texture_spin(spinin,1);
    dev_Hopping_Matrix<<<gridsize, blocksize>>>
             (dev_gf, spinin, dev_spin_eo1, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0); //dev_spin_eo1 == even -> 0
  unbind_texture_spin(1);
  dev_mul_one_pm_imu_inv<<<gridsize, blocksize>>>(dev_spin_eo1,dev_spin_eo2, -1.);
  
  bind_texture_spin(dev_spin_eo2,1);
    dev_Hopping_Matrix<<<gridsize, blocksize>>>
            (dev_gf, dev_spin_eo2, dev_spin_eo1, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1); 
  unbind_texture_spin(1);
  dev_mul_one_pm_imu_sub_mul_gamma5<<<gridsize, blocksize>>>(spinin, dev_spin_eo1,  dev_spin_eo2, -1.);
  
  //Q_{+}
  bind_texture_spin(dev_spin_eo2,1);
    dev_Hopping_Matrix<<<gridsize, blocksize>>>
          (dev_gf, dev_spin_eo2, dev_spin_eo1, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0); //dev_spin_eo1 == even -> 0
  unbind_texture_spin(1);
  dev_mul_one_pm_imu_inv<<<gridsize, blocksize>>>(dev_spin_eo1,spinout, +1.);
  
  bind_texture_spin(spinout,1);
    dev_Hopping_Matrix<<<gridsize, blocksize>>>
             (dev_gf, spinout, dev_spin_eo1, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1); 
  unbind_texture_spin(1);
  dev_mul_one_pm_imu_sub_mul_gamma5<<<gridsize, blocksize>>>(dev_spin_eo2, dev_spin_eo1,  spinout , +1.); 
  
}






//applies the full tm Operator
// uses texture cache (spin_tex) for input spinor
// runs through whole lattice for output spinor
// D_psi uses phase_mu and not ka_mu for the boundary conds (vice versa in HoppingMatrix) 
// -> thats why complexmult and complexcgmult are interchanged in dev_HoppingMatrix and in 
// dev_tm_dirac_kappa
__global__ void dev_tm_dirac_kappa(dev_su3_2v * gf, dev_spinor * sin, dev_spinor * sout, int * dev_nn){
    int pos,hoppos;
    dev_spinor shelp1[6], ssum[6];
    __shared__ dev_su3 gfsmem[BLOCK];
    

  pos= threadIdx.x + blockDim.x*blockIdx.x;  
  int ix = threadIdx.x;
  if(pos < dev_VOLUME){
        
          //dev_zero_spinor(&(ssum[0])); // zero sum
          //skalarer Term
          ssum[0] = tex1Dfetch(spin_tex,6*pos);
          ssum[1] = tex1Dfetch(spin_tex,6*pos+1);
          ssum[2] = tex1Dfetch(spin_tex,6*pos+2);
          ssum[3] = tex1Dfetch(spin_tex,6*pos+3);
          ssum[4] = tex1Dfetch(spin_tex,6*pos+4);
          ssum[5] = tex1Dfetch(spin_tex,6*pos+5);

          
//hopping term                
//l==0,t
            //positive direction
            hoppos = dev_nn[8*pos];
            //color
            dev_reconstructgf_2vtexref(4*pos,&(gfsmem[ix]));
            dev_su3MtV_spintex(gfsmem[ix], hoppos, &(shelp1[0]));
            //-kappa(r - gamma_mu)
            dev_complexmult_add_assign_spinor(&(ssum[0]),dev_mk0,&(shelp1[0]), &(ssum[0]));
            //dev_GammatV(0,&(shelp1[0]));
            dev_Gamma0(&(shelp1[0]));
            dev_complexmult_add_assign_spinor(&(ssum[0]),dev_k0,&(shelp1[0]), &(ssum[0]));

            //negative direction
            hoppos = dev_nn[8*pos+4];
            //color
            dev_reconstructgf_2vtexref_dagger(4*hoppos,&(gfsmem[ix]));
            dev_su3MtV_spintex(gfsmem[ix], hoppos, &(shelp1[0]));     
            //-kappa(r + gamma_mu)
            dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_mk0,&(shelp1[0]), &(ssum[0]));
            //dev_GammatV(0,&(shelp1[0]));
            dev_Gamma0(&(shelp1[0]));
            dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_mk0,&(shelp1[0]), &(ssum[0]));


//l==3,z               
            //positive direction
            hoppos = dev_nn[8*pos+3];
            //color
            dev_reconstructgf_2vtexref(4*pos+(3),&(gfsmem[ix]));
            dev_su3MtV_spintex(gfsmem[ix], hoppos, &(shelp1[0]));
            //-kappa(r - gamma_mu)
            dev_complexmult_add_assign_spinor(&(ssum[0]),dev_mk3,&(shelp1[0]), &(ssum[0]));
            //dev_GammatV(3,&(shelp1[0]));
            dev_Gamma3(&(shelp1[0]));
            dev_complexmult_add_assign_spinor(&(ssum[0]),dev_k3,&(shelp1[0]), &(ssum[0]));

            //negative direction
            hoppos = dev_nn[8*pos+7];
            //color
            dev_reconstructgf_2vtexref_dagger(4*hoppos+(3),&(gfsmem[ix]));
            dev_su3MtV_spintex(gfsmem[ix], hoppos, &(shelp1[0]));
            //-kappa(r + gamma_mu)
            dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_mk3,&(shelp1[0]), &(ssum[0]));
            //dev_GammatV(3,&(shelp1[0]));
            dev_Gamma3(&(shelp1[0]));
            dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_mk3,&(shelp1[0]), &(ssum[0]));
         
         
//l==2,y        
            //positive direction
            hoppos = dev_nn[8*pos+2];
            //color
            dev_reconstructgf_2vtexref(4*pos+(2),&(gfsmem[ix]));
            dev_su3MtV_spintex(gfsmem[ix], hoppos, &(shelp1[0]));
            //-kappa(r - gamma_mu)
            dev_complexmult_add_assign_spinor(&(ssum[0]),dev_mk2,&(shelp1[0]), &(ssum[0]));
            //dev_GammatV(2,&(shelp1[0]));
            dev_Gamma2(&(shelp1[0]));
            dev_complexmult_add_assign_spinor(&(ssum[0]),dev_k2,&(shelp1[0]), &(ssum[0]));
            
            //negative direction
            hoppos = dev_nn[8*pos+6];
            //color
            dev_reconstructgf_2vtexref_dagger(4*hoppos+(2),&(gfsmem[ix]));
            dev_su3MtV_spintex(gfsmem[ix], hoppos, &(shelp1[0]));
            //-kappa(r + gamma_mu)
            dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_mk2,&(shelp1[0]), &(ssum[0]));
            //dev_GammatV(2,&(shelp1[0]));
            dev_Gamma2(&(shelp1[0]));
            dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_mk2,&(shelp1[0]), &(ssum[0]));


//l==1,x 
            //positive direction
            hoppos = dev_nn[8*pos+1];
            //color
            dev_reconstructgf_2vtexref(4*pos+(1),&(gfsmem[ix]));
            dev_su3MtV_spintex(gfsmem[ix], hoppos, &(shelp1[0]));
            //-kappa(r - gamma_mu)
            dev_complexmult_add_assign_spinor(&(ssum[0]),dev_mk1,&(shelp1[0]), &(ssum[0]));
            //dev_GammatV(1,&(shelp1[0]));
            dev_Gamma1(&(shelp1[0]));
            dev_complexmult_add_assign_spinor(&(ssum[0]),dev_k1,&(shelp1[0]), &(ssum[0]));

            //negative direction
            hoppos = dev_nn[8*pos+5];
            //color
            dev_reconstructgf_2vtexref_dagger(4*hoppos+(1),&(gfsmem[ix]));
            dev_su3MtV_spintex(gfsmem[ix], hoppos, &(shelp1[0]));
            //-kappa(r + gamma_mu)
            dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_mk1,&(shelp1[0]), &(ssum[0]));
            //dev_GammatV(1,&(shelp1[0]));
            dev_Gamma1(&(shelp1[0]));
            dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_mk1,&(shelp1[0]), &(ssum[0]));  
          
          
          
          //gamma5 term
          shelp1[0] = tex1Dfetch(spin_tex,6*pos);
          shelp1[1] = tex1Dfetch(spin_tex,6*pos+1);
          shelp1[2] = tex1Dfetch(spin_tex,6*pos+2);
          shelp1[3] = tex1Dfetch(spin_tex,6*pos+3);
          shelp1[4] = tex1Dfetch(spin_tex,6*pos+4);
          shelp1[5] = tex1Dfetch(spin_tex,6*pos+5);
          
          
          
          //dev_GammatV(4,&(shelp1[0]));
          dev_Gamma5(&(shelp1[0]));
          dev_complexmult_add_assign_spinor(&(ssum[0]),dev_initcomplex(0.0,2.0*kappa*mu),&(shelp1[0]), &(sout[6*pos]));
  }
}





extern "C" void dev_tm_dirac_dagger_kappa(dev_su3_2v * gf,dev_spinor* spinin, dev_spinor* spinout, 
 int *grid, int * nn_grid, REAL* output,REAL* erg, int xsize, int ysize){
 int gridsize;
 if( VOLUME >= 128){
   gridsize =VOLUME/128;
 }
 else{
   gridsize=1;
 }
 dim3 griddim2(gridsize,1,1);
 dim3 blockdim2(128,1,1);
 dim3 blockdim(xsize,ysize);
 
  dim3 blockdim3(BLOCK,1,1);
 if( VOLUME >= BLOCK){
   gridsize =VOLUME/BLOCK;
 }
 else{
   gridsize=1;
 }
 dim3 griddim3(gridsize,1,1); 
  dev_gamma5 <<<griddim2, blockdim2 >>> (spinin,spinout);
  dev_tm_dirac_kappa <<<griddim3, blockdim3 >>> (gf, spinout, spinin, dev_nn);
  dev_gamma5 <<<griddim2, blockdim2 >>>(spinin,spinout);
}







__device__ inline REAL dev_skalarprod_spinor(dev_spinor * s1, dev_spinor * s2){
  REAL skalprod = 0.0;
  int i;
  #pragma unroll 6
  for(i=0;i<6;i++){ //color + spin
    skalprod += ((*(s1+i)).x*(*(s2+i)).x + (*(s1+i)).y*(*(s2+i)).y + (*(s1+i)).z*(*(s2+i)).z + (*(s1+i)).w*(*(s2+i)).w);
  }
  return skalprod;
}




__device__ inline REAL dev_squarenorm_spinor(dev_spinor * s1){
  REAL skalprod = 0.0;
  int i;
  #pragma unroll 6
  for(i=0;i<6;i++){ //color + spin
    skalprod += ((*(s1+i)).x*(*(s1+i)).x + (*(s1+i)).y*(*(s1+i)).y + (*(s1+i)).z*(*(s1+i)).z + (*(s1+i)).w*(*(s1+i)).w);
  }
  return skalprod;
}



__device__ inline REAL dev_squarenorm_spinor_tex(int pos){
  REAL skalprod = 0.0;
  int i;
  float4 help;
  
  #pragma unroll 6
  for(i=0;i<6;i++){ //color + spin
    help = tex1Dfetch(spin_tex2,6*pos+i);
    skalprod += help.x*help.x + help.y*help.y + help.z*help.z + help.w*help.w;
  }
  return skalprod;
}




//only 1 dim parallel possible, because need __syncthread !
__global__ void dev_skalarprod_spinor_field2(dev_spinor* s1, dev_spinor* s2, REAL* erg){
  __shared__ REAL shrinkarray[ACCUM_N];
  int pos,stepwidth;
  REAL ks,kc,ds,tr,ts,tt;
  
   
   // ADD ERROR HERE if t > maxblockdim
   
   ks=0.0;
   kc=0.0; 
   
   if(blockDim.x > dev_VOLUME){
     stepwidth = 1;  
   }
   else{
     stepwidth = dev_VOLUME/(gridDim.x*blockDim.x);
   }
   
     int start = (blockIdx.x*blockDim.x + threadIdx.x)*stepwidth;
     int end = (blockIdx.x*blockDim.x + threadIdx.x+1)*stepwidth;  
   
   for(pos=start;pos<end; pos++){
     if(pos < dev_VOLUME){
          ds = dev_skalarprod_spinor(&(s1[6*pos]), &(s2[6*pos] ));
      }
      else{
        ds = 0.0;
      } 
          // Kahan summation
          tr=ds+kc;
          ts=tr+ks;
          tt=ts-ks;
          ks=ts;
          kc=tr-tt;
    }
    kc=ks+kc;
    shrinkarray[threadIdx.x] = kc;
    __syncthreads();
    
    
    if(threadIdx.x==0){
      ks=0.0;
      kc=0.0; 
      int k;
      for(k=0; k<blockDim.x; k++){
        ds = shrinkarray[k];
        tr=ds+kc;
        ts=tr+ks;
        tt=ts-ks;
        ks=ts;
        kc=tr-tt;
      }
      kc=ks+kc;
      (erg[blockIdx.x])=kc;
    }//threadIdx==0;
}








//only 1 dim parallel possible, because need __syncthread !
__global__ void dev_squarenorm_spinor_field(dev_spinor* s1, REAL* erg){
  __shared__ REAL shrinkarray[ACCUM_N];
  int pos,stepwidth;
  REAL ks,kc,ds,tr,ts,tt;
  
   
   // ADD ERROR HERE if t > maxblockdim
   
   ks=0.0;
   kc=0.0; 
   
   if(blockDim.x > dev_VOLUME){
     stepwidth = 1;  
   }
   else{
     stepwidth = dev_VOLUME/(gridDim.x*blockDim.x);
   }
   
     int start = (blockIdx.x*blockDim.x + threadIdx.x)*stepwidth;
     int end = (blockIdx.x*blockDim.x + threadIdx.x+1)*stepwidth;  
   
   for(pos=start;pos<end; pos++){
     if(pos < dev_VOLUME){
         //ds = dev_squarenorm_spinor_tex(pos);
         ds = dev_squarenorm_spinor(&(s1[6*pos]));
      }
      else{
        ds = 0.0;
      } 
          // Kahan summation
          tr=ds+kc;
          ts=tr+ks;
          tt=ts-ks;
          ks=ts;
          kc=tr-tt;
    }
    kc=ks+kc;
    shrinkarray[threadIdx.x] = kc;
    __syncthreads();
    
    
    if(threadIdx.x==0){
      ks=0.0;
      kc=0.0; 
      int k;
      for(k=0; k<blockDim.x; k++){
        ds = shrinkarray[k];
        tr=ds+kc;
        ts=tr+ks;
        tt=ts-ks;
        ks=ts;
        kc=tr-tt;
      }
      kc=ks+kc;
      (erg[blockIdx.x])=kc;
    }//threadIdx==0;
}






//only 1 dim parallel, because need __syncthread !
__global__ void dev_skalarprod_spinor_field(dev_spinor* s1, dev_spinor* s2, REAL* erg){
  __shared__ REAL shrinkarray[ACCUM_N];
  int pos,stepwidth, sweepsperthread;
  REAL ks,kc,ds,tr,ts,tt;
   
   // ADD ERROR HERE if t > maxblockdim
   
   ks=0.0;
   kc=0.0; 
   
   if(ACCUM_N > dev_VOLUME){
     stepwidth = 1;
     sweepsperthread = 1;  
   }
   else{
     stepwidth = dev_VOLUME/ACCUM_N;
     sweepsperthread = ACCUM_N/blockDim.x;
   }
    
   
   
 for(int j = 0; j < sweepsperthread; j++){
   
     int start = (threadIdx.x + j*blockDim.x)*stepwidth;
     int end = (threadIdx.x+j*blockDim.x+1)*stepwidth;  
     ks=0.0;
     kc=0.0; 
     
   for(pos=start;pos<end; pos++){
     if(pos < dev_VOLUME){
          ds = dev_skalarprod_spinor(&(s1[6*pos]), &(s2[6*pos] ));
          
      }
      else{
        ds = 0.0;
      } 
          // Kahan summation
          tr=ds+kc;
          ts=tr+ks;
          tt=ts-ks;
          ks=ts;
          kc=tr-tt;
    }
    kc=ks+kc;
    shrinkarray[threadIdx.x+j*blockDim.x] = kc;
  }
  __syncthreads();
   
 
    for(int stride = ACCUM_N / 2; stride > 0; stride >>= 1){
       __syncthreads();
       for(int iAccum = threadIdx.x; iAccum < stride; iAccum += blockDim.x)
           shrinkarray[iAccum] += shrinkarray[stride + iAccum];
    }

    if(threadIdx.x == 0) (*erg) = shrinkarray[0];
    
    
    /*
    if(threadIdx.x==0){
      ks=0.0;
      kc=0.0; 
      int k;
      for(k=0; k<sweepsperthread*blockDim.x; k++){
        ds = shrinkarray[k];
        tr=ds+kc;
        ts=tr+ks;
        tt=ts-ks;
        ks=ts;
        kc=tr-tt;
      }
      kc=ks+kc;
      (*erg)=kc;
    }//threadIdx==0;
    */
    
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






extern "C" int bind_texture_spin(dev_spinor* s, int i){
  
  size_t size;
  if(even_odd_flag){
    size = sizeof(float4)*6*VOLUME/2;
  }
  else{
    size = sizeof(float4)*6*VOLUME;
  }
   
  
  switch(i){
    case 1:
      //printf("Binding texture to spinorfield 1\n");
      spin_texRefPtr = NULL;
      cudaGetTextureReference(&spin_texRefPtr, "spin_tex");
      spin_channelDesc =  cudaCreateChannelDesc<float4>();
      cudaBindTexture(0, spin_texRefPtr, s, &spin_channelDesc, size);
      //printf("%s\n", cudaGetErrorString(cudaGetLastError())); 
      return(0);
   
    case 2:
      //printf("Binding texture to spinorfield 2\n");
      spin_texRefPtr2 = NULL;
      cudaGetTextureReference(&spin_texRefPtr2, "spin_tex2");
      spin_channelDesc2 =  cudaCreateChannelDesc<float4>();
      cudaBindTexture(0, spin_texRefPtr2, s, &spin_channelDesc2, size);
      //printf("%s\n", cudaGetErrorString(cudaGetLastError()));  
      return(0);
  }
return(1);  
}


extern "C" int unbind_texture_spin(int i){
  switch(i){
    case 1:
      //printf("Unbinding texture of spinorfield 1\n");
      cudaUnbindTexture(spin_texRefPtr);
      //printf("%s\n", cudaGetErrorString(cudaGetLastError()));    
      return(0);
    case 2:
      //printf("Unbinding texture of spinorfield 2\n");
      cudaUnbindTexture(spin_texRefPtr2);
      //printf("%s\n", cudaGetErrorString(cudaGetLastError()));    
      return(0);    
  }
  
return(1);
}



extern "C" int bind_texture_gf(dev_su3_2v * gf){
 //printf("Binding texture to gaugefield\n");
 size_t size = sizeof(float4)*3*VOLUME*4;
 cudaGetTextureReference(&gf_texRefPtr, "gf_tex");
 gf_channelDesc =  cudaCreateChannelDesc<float4>();
 cudaBindTexture(0, gf_texRefPtr, gf, &gf_channelDesc, size);
 printf("%s\n", cudaGetErrorString(cudaGetLastError()));    
 return(0);
}


extern "C" int unbind_texture_gf(){
 //printf("Unbinding texture to gaugefield\n");
 cudaUnbindTexture(gf_texRefPtr);
 printf("%s\n", cudaGetErrorString(cudaGetLastError()));    
 return(0);
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
   gridsize =VOLUME/BLOCK;
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
  
 bind_texture_gf(gf);
 bind_texture_spin(spin4,1);
  
  // apply D_tm
  dev_tm_dirac_kappa <<<griddim3, blockdim3 >>>(gf, spin4, spinout, nn_grid);


  unbind_texture_gf();
  unbind_texture_spin(1);
}






extern "C" void dev_cg(dev_su3_2v * gf,dev_spinor* spinin, dev_spinor* spinout, 
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

 dim3 blockdim(1,1);
 dim3 blockdim2(128,1,1);
 if( VOLUME >= 128){
   gridsize =VOLUME/128;
 }
 else{
   gridsize=1;
 }
 dim3 griddim2(gridsize,1,1);
 
 
 if(VOLUME%BLOCK != 0){
   printf("Error: VOLUME is not a multiple of BLOCK. Aborting...\n");
   exit(100);
 }
 dim3 blockdim3(BLOCK,1,1);
 if( VOLUME >= BLOCK){
   gridsize =VOLUME/BLOCK;
 }
 else{
   gridsize=1;
 }
 dim3 griddim3(gridsize,1,1); 
 
 size_t size2 = sizeof(float4)*6*VOLUME;
 
 //Bind texture gf
 bind_texture_gf(gf);
 //Bind texture spinor to spin4 (D_tm is always applied to spin4)
 bind_texture_spin(spin4,1);
 
 
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
  
  unbind_texture_spin(1);
     // GAMMA5, mu -> -mu
     dev_gamma5 <<<griddim2, blockdim2 >>> (spin2,spin4);
     dev_swapmu <<<1,1>>> ();
  bind_texture_spin(spin4,1);
     //D_tm 
     dev_tm_dirac_kappa <<<griddim3, blockdim3 >>> (gf, spin4, spin3, dev_nn);
  unbind_texture_spin(1);
     //GAMMA5 mu -> -mu
     dev_gamma5 <<<griddim2, blockdim2 >>>(spin3,spin4);
     dev_swapmu <<<1,1>>> ();
  bind_texture_spin(spin4,1);
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
    unbind_texture_spin(1);
      dev_gamma5 <<<griddim2, blockdim2 >>> (spin1,spin4);
      dev_swapmu <<<1,1>>> ();
    bind_texture_spin(spin4,1);
   
      //D_tm GAMMA5, mu -> -mu
      dev_tm_dirac_kappa <<<griddim3, blockdim3 >>> (gf, spin4, spin3, dev_nn);
      dev_gamma5 <<<griddim2, blockdim2 >>>(spin3,spinout);
      dev_swapmu <<<1,1>>> ();
  
    //printf("Unbinding texture of spinorfield\n");
    unbind_texture_spin(1);
    cudaMemcpy(spin4, spinout,size2, cudaMemcpyDeviceToDevice);
    //printf("Rebinding texture to spinorfield\n");
    bind_texture_spin(spin4,1);
      
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
    unbind_texture_spin(1);
      dev_gamma5 <<<griddim2, blockdim2 >>> (spin1,spin4);
      dev_swapmu <<<1,1>>> ();
    bind_texture_spin(spin4,1);
      dev_tm_dirac_kappa <<<griddim3, blockdim3 >>> (gf, spin4, spin3, dev_nn);
      dev_gamma5 <<<griddim2, blockdim2 >>>(spin3,spin1);
      dev_swapmu <<<1,1>>> ();
    unbind_texture_spin(1);


 //go over to non-kappa, Ddagger = g5 D g5
 dev_skalarmult_spinor_field<<<griddim2, blockdim2 >>>(spin1,1.0/(scaleparam*scaleparam), spinout);  
 
  // times operator == source ?? 
  //dev_tm_dirac_kappa<<<griddim3, blockdim3 >>>(gf, spin3, spinout, nn_grid);
  }
  else{
   dev_copy_spinor_field<<<griddim2, blockdim2 >>>(spin1,spinout);
  }
  
  
  unbind_texture_gf();
  cudaFree(dotprod);
  cudaFree(dotprod2);
  cudaFree(rk);
  cudaFree(alpha);
  cudaFree(beta);
  cublasShutdown();
}





// this is the eo version of the device cg inner solver 
// we invert the hermitean Q_{-} Q_{+}
extern "C" void dev_cg_eo(dev_su3_2v * gf,dev_spinor* spinin, dev_spinor* spinout, 
dev_spinor* spin0, dev_spinor* spin1, dev_spinor* spin2, dev_spinor* spin3, dev_spinor* spin4, int *grid, int * nn_grid, REAL* output,REAL* erg, int xsize, int ysize, int rescalekappa){
 
 
 REAL host_alpha, host_beta, host_dotprod, host_rk, sourcesquarenorm;
 REAL * dotprod, * dotprod2, * rk, * alpha, *beta;
 
 
 
 int i, gridsize;
 int maxit = max_innersolver_it;
 REAL eps = (REAL) innersolver_precision;
 int N_recalcres = 10; // after N_recalcres iterations calculate r = A x_k - b
 
 

 dim3 blockdim(1,1);
 dim3 blockdim2(128,1,1);
 if( VOLUME/2 >= 128){
   gridsize =VOLUME/2/128;
 }
 else{
   gridsize=1;
 }
 dim3 griddim2(gridsize,1,1);
 
 
 if((VOLUME/2)%BLOCK != 0){
   printf("Error: VOLUME/2 is not a multiple of BLOCK. Aborting...\n");
   exit(100);
 }
 int blockdim3=BLOCK;
 if( VOLUME/2 >= BLOCK){
   gridsize =VOLUME/2/BLOCK;
 }
 else{
   gridsize=1;
 }
 int griddim3=gridsize; 
 
 size_t size2 = sizeof(float4)*6*VOLUME/2;
 
 //Bind texture gf
 bind_texture_gf(gf);
 
 
 //Initialize some stuff
  printf("mu = %f\n", g_mu);
  dev_complex h0,h1,h2,h3;
  h0.re = (REAL)ka0.re;    h0.im = -(REAL)ka0.im;
  h1.re = (REAL)ka1.re;    h1.im = -(REAL)ka1.im;
  h2.re = (REAL)ka2.re;    h2.im = -(REAL)ka2.im;
  h3.re = (REAL)ka3.re;    h3.im = -(REAL)ka3.im;
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
  dev_Qtm_pm_psi(spin2, spin3, griddim3, blockdim3);
  
  
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
  
 if ((host_dotprod <= eps*sourcesquarenorm)){//error-limit erreicht
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
    dev_Qtm_pm_psi(spin1, spin3, griddim3, blockdim3);
      
        
    
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
  
  
  unbind_texture_gf();
  cudaFree(dotprod);
  cudaFree(dotprod2);
  cudaFree(rk);
  cudaFree(alpha);
  cudaFree(beta);
  cublasShutdown();
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




// get 2 first rows of gf float4 type
//  
//
void su3to2vf4(su3** gf, dev_su3_2v* h2d_gf){
  int i,j;
  for (i=0;i<VOLUME;i++){
   for(j=0;j<4;j++){
   //first row
    h2d_gf[3*(4*i+j)].x = (REAL) gf[i][j].c00.re;
    h2d_gf[3*(4*i+j)].y = (REAL) gf[i][j].c00.im;
    h2d_gf[3*(4*i+j)].z = (REAL) gf[i][j].c01.re;
    h2d_gf[3*(4*i+j)].w = (REAL) gf[i][j].c01.im;
    h2d_gf[3*(4*i+j)+1].x = (REAL) gf[i][j].c02.re;
    h2d_gf[3*(4*i+j)+1].y = (REAL) gf[i][j].c02.im;      
   //second row
    h2d_gf[3*(4*i+j)+1].z = (REAL) gf[i][j].c10.re;
    h2d_gf[3*(4*i+j)+1].w = (REAL) gf[i][j].c10.im;
    h2d_gf[3*(4*i+j)+2].x = (REAL) gf[i][j].c11.re;
    h2d_gf[3*(4*i+j)+2].y = (REAL) gf[i][j].c11.im;
    h2d_gf[3*(4*i+j)+2].z = (REAL) gf[i][j].c12.re;
    h2d_gf[3*(4*i+j)+2].w = (REAL) gf[i][j].c12.im;      
  } 
 }
}







// this is to reconstruc the gf on the host from 2 rows of the link
// may be used for tests
void reconstructgf_2v (dev_su3* gf){
  complex help1;
  complex help2;
  //third row from cconj(cross product of first and second row)
  _mult_assign_complex(help1,(*gf)[0][1],(*gf)[1][2]);
  _mult_assign_complex(help2,(*gf)[0][2],(*gf)[1][1]);
  _diff_complex(help1,help2);
  help1.im = -help1.im;
  (*gf)[2][0].re = help1.re;
  (*gf)[2][0].im = help1.im;
  
  _mult_assign_complex(help1,(*gf)[0][2],(*gf)[1][0]);
  _mult_assign_complex(help2,(*gf)[0][0],(*gf)[1][2]);
  _diff_complex(help1,help2);
  help1.im = -help1.im;
  (*gf)[2][1].re = help1.re;
  (*gf)[2][1].im = help1.im;
  
  _mult_assign_complex(help1,(*gf)[0][0],(*gf)[1][1]);
  _mult_assign_complex(help2,(*gf)[0][1],(*gf)[1][0]);
  _diff_complex(help1,help2);
  help1.im = -help1.im;
  (*gf)[2][2].re = help1.re;
  (*gf)[2][2].im = help1.im;
  return;
}



void showcompare_gf(int t, int x, int y, int z, int mu){
   int ind1 = g_ipt[t][x][y][z];
   su3 ** gf1 = g_gauge_field;
   printf("(%f,%f)\t(%f,%f)\t(%f,%f)\n",gf1[ind1][mu].c00.re,
   					gf1[ind1][mu].c00.im,
   					gf1[ind1][mu].c01.re,
   					gf1[ind1][mu].c01.im,
   					gf1[ind1][mu].c02.re,
   					gf1[ind1][mu].c02.im
   );
   printf("(%f,%f)\t(%f,%f)\t(%f,%f)\n",gf1[ind1][mu].c10.re,
   					gf1[ind1][mu].c10.im,
   					gf1[ind1][mu].c11.re,
   					gf1[ind1][mu].c11.im,
   					gf1[ind1][mu].c12.re,
   					gf1[ind1][mu].c12.im
   );
   printf("(%f,%f)\t(%f,%f)\t(%f,%f)\n",gf1[ind1][mu].c20.re,
   					gf1[ind1][mu].c20.im,
   					gf1[ind1][mu].c21.re,
   					gf1[ind1][mu].c21.im,
   					gf1[ind1][mu].c22.re,
   					gf1[ind1][mu].c22.im
   );
   printf("\n\n");

   int ind2 =  z + LZ*(y + LY*(x + LX*t));
   printf("(%f,%f)\t(%f,%f)\t(%f,%f)\n",h2d_gf[3*(4*ind2+mu)].x,
   					h2d_gf[3*(4*ind2+mu)].y,
   					h2d_gf[3*(4*ind2+mu)].z,
   					h2d_gf[3*(4*ind2+mu)].w,
   					h2d_gf[3*(4*ind2+mu)+1].x,
   					h2d_gf[3*(4*ind2+mu)+1].y
   );
   printf("(%f,%f)\t(%f,%f)\t(%f,%f)\n",h2d_gf[3*(4*ind2+mu)+1].z,
   					h2d_gf[3*(4*ind2+mu)+1].w,
   					h2d_gf[3*(4*ind2+mu)+2].x,
   					h2d_gf[3*(4*ind2+mu)+2].y,
   					h2d_gf[3*(4*ind2+mu)+2].z,
   					h2d_gf[3*(4*ind2+mu)+2].w
   );   
   
   dev_su3 help;
   
   help[0][0].re = h2d_gf[3*(4*ind2+mu)].x;
   help[0][0].im = h2d_gf[3*(4*ind2+mu)].y;
   help[0][1].re = h2d_gf[3*(4*ind2+mu)].z;
   help[0][1].im = h2d_gf[3*(4*ind2+mu)].w;

   help[0][2].re = h2d_gf[3*(4*ind2+mu)+1].x;
   help[0][2].im = h2d_gf[3*(4*ind2+mu)+1].y;
   help[1][0].re = h2d_gf[3*(4*ind2+mu)+1].z;
   help[1][0].im = h2d_gf[3*(4*ind2+mu)+1].w;
   
   help[1][1].re = h2d_gf[3*(4*ind2+mu)+2].x;
   help[1][1].im = h2d_gf[3*(4*ind2+mu)+2].y;
   help[1][2].re = h2d_gf[3*(4*ind2+mu)+2].z;
   help[1][2].im = h2d_gf[3*(4*ind2+mu)+2].w;   
   
   reconstructgf_2v (&help); 

   printf("(%f,%f)\t(%f,%f)\t(%f,%f)\n",help[2][0].re,
   					help[2][0].im,
   					help[2][1].re,
   					help[2][1].im,
   					help[2][2].re,
   					help[2][2].im
   );
     
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

  
  /* allocate 2 rows of gf = 3*4*VOLUME float4's*/
  size_t dev_gfsize = 3*4*VOLUME * sizeof(dev_su3_2v); 
  
  if((cudaerr=cudaMalloc((void **) &dev_gf, dev_gfsize)) != cudaSuccess){
    printf("Error in init_mixedsolve(): Memory allocation of gauge field failed. Aborting...\n");
    exit(200);
  }   // Allocate array on device
  else{
    printf("Allocated gauge field on device\n");
  }  
  
  
  h2d_gf = (dev_su3_2v *)malloc(dev_gfsize); // Allocate REAL conversion gf on host
  su3to2vf4(gf,h2d_gf);
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

  /* allocate 2 rows of gf = 3*4*VOLUME float4's*/
  size_t dev_gfsize = 3*4*VOLUME * sizeof(dev_su3_2v); 
  
  
  
  if((cudaerr=cudaMalloc((void **) &dev_gf, dev_gfsize)) != cudaSuccess){
    printf("Error in init_mixedsolve(): Memory allocation of gauge field failed. Aborting...\n");
    exit(200);
  }   // Allocate array on device
  else{
    printf("Allocated gauge field on device\n");
  }  
  
  
  h2d_gf = (dev_su3_2v *)malloc(dev_gfsize); // Allocate REAL conversion gf on host
  su3to2vf4(gf,h2d_gf);
  cudaMemcpy(dev_gf, h2d_gf, dev_gfsize, cudaMemcpyHostToDevice);



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
  
  
  output_size = LZ*T*sizeof(float); // parallel in t and z direction
  cudaMalloc((void **) &dev_output, output_size);   // output array
  float * host_output = (float*) malloc(output_size);

  int grid[5];
  grid[0]=LX; grid[1]=LY; grid[2]=LZ; grid[3]=T; grid[4]=VOLUME/2; 
  // dev_VOLUME is half of VOLUME for eo
 
  cudaMalloc((void **) &dev_grid, 5*sizeof(int));
  cudaMemcpy(dev_grid, &(grid[0]), 5*sizeof(int), cudaMemcpyHostToDevice);
  
  

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
  
  
  
  free(h2d_spin);
  free(h2d_gf);
}







extern "C" int mixed_solve (spinor * const P, spinor * const Q, const int max_iter, 
	   double eps, const int rel_prec,const int N){
  
  // source in Q, initial solution in P (not yet implemented)
  double rk;
  int outercount=0;
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
    return(iter*max_innersolver_it);  // MAYBE ONE SHOULD KEEP TRACK OF REAL INNER SOLVER STEPS
   }
   

  //initialize spin fields on device
  convert2REAL4_spin(g_spinor_field[DUM_SOLVER],h2d_spin);
  
  cudaMemcpy(dev_spinin, h2d_spin, dev_spinsize, cudaMemcpyHostToDevice);
  printf("%s\n", cudaGetErrorString(cudaGetLastError()));

   // solve in single prec on device
   // D p_k = r_k
   printf("Entering inner solver\n");
   assert((startinner = clock())!=-1);
   dev_cg(dev_gf, dev_spinin, dev_spinout, dev_spin1, dev_spin2, dev_spin3, dev_spin4, dev_spin5, dev_grid,dev_nn, dev_output,NULL, T, LZ,0);
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








extern "C" int mixed_solve_eo (spinor * const P, spinor * const Q, const int max_iter, 
	   double eps, const int rel_prec, const int N){
  
  // source in Q, initial solution in P (not yet implemented)
  double rk;
  int outercount=0;
  clock_t start, stop, startinner, stopinner; 
  double timeelapsed = 0.0;
  double sourcesquarenorm;
  int iter;
  
  
  size_t dev_spinsize = 6*VOLUME/2 * sizeof(dev_spinor); // float4 even-odd !
  init_mixedsolve_eo(g_gauge_field);
  
  // Start timer
  assert((start = clock())!=-1);
  
  rk = square_norm(Q, N, 1);
  sourcesquarenorm=rk; // for relative prec
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
   
   printf("Residue after %d inner solver iterations: %.18e\n",outercount,rk);
   
   if(((rk <= eps) && (rel_prec == 0)) || ((rk <= eps*sourcesquarenorm) && (rel_prec == 1)))
   {
     printf("Reached solver precision of eps=%.2e\n",eps);
     //multiply with D^dagger
     Qtm_minus_psi(g_spinor_field[DUM_SOLVER+3], g_spinor_field[DUM_SOLVER+1]);
     assign(P, g_spinor_field[DUM_SOLVER+3], N);
  

     stop = clock();
     timeelapsed = (double) (stop-start)/CLOCKS_PER_SEC;
     printf("EO Inversion done in mixed precision.\n Number of iterations in outer solver: %d\n Squared residue: %.8e\n Time elapsed: %.6e sec\n", outercount, rk, timeelapsed);
   
     finalize_mixedsolve();
     return(iter*max_innersolver_it);  // MAYBE ONE SHOULD KEEP TRACK OF REAL INNER SOLVER STEPS
   }
   
  //initialize spin fields on device
  convert2REAL4_spin(g_spinor_field[DUM_SOLVER],h2d_spin);
  
  cudaMemcpy(dev_spinin, h2d_spin, dev_spinsize, cudaMemcpyHostToDevice);
  printf("%s\n", cudaGetErrorString(cudaGetLastError()));

   // solve in single prec on device
   // D p_k = r_k
   printf("Entering inner solver\n");
   assert((startinner = clock())!=-1);
   dev_cg_eo(dev_gf, dev_spinin, dev_spinout, dev_spin1, dev_spin2, dev_spin3, dev_spin4, dev_spin5, dev_grid,dev_nn, dev_output,NULL, T, LZ,0);
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


  finalize_mixedsolve();
  return(-1);
}








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
 * File: linalg.cuh
 *
 * CUDA linear algebra functions and implementation of gamma-multiplication
 *
 * 
 *
 **************************************************************************/
 





// ******************** complex double operations ***********************

__device__ inline dev_complex_d dev_cconj_d (dev_complex_d c){ /*konjugiert komplexe Zahl*/
 dev_complex_d erg;
 erg.re = c.re;
 erg.im = -1.0*c.im;
return erg;
}

__device__ inline void dev_ccopy_d(dev_complex_d* von, dev_complex_d* nach){/*kopiert complex von nach complex nach*/
  nach->re = von->re;
  nach->im = von->im;
}

__device__ inline double dev_cabssquare_d (dev_complex_d c){ /*gibt abs^2 einer komplexen Zahl zurück*/
 return c.re*c.re + c.im*c.im;
}

__device__ inline double dev_cabsolute_d (dev_complex_d c){/*gibt Betrag einer kompl. zahl zurück*/
 return sqrt(c.re*c.re + c.im*c.im);
}


__device__ inline  dev_complex_d dev_crealmult_d(dev_complex_d c1, double real){ /*multipliziert c1 mit reeller zahl re*/
  dev_complex_d erg;
  erg.re = real*c1.re;
  erg.im = real*c1.im;
return erg;
}

__device__ inline dev_complex_d dev_cmult_d (dev_complex_d c1, dev_complex_d c2){ /*multiplizier zwei komplexe Zahlen*/
  dev_complex_d erg;
  erg.re = c1.re * c2.re;
  erg.re -= c1.im * c2.im;
  erg.im = c1.re * c2.im;
  erg.im += c1.im * c2.re;
return erg;
}

__device__ inline dev_complex_d dev_cadd_d (dev_complex_d c1, dev_complex_d c2){ /*addiert zwei komplexe Zahlen */
  dev_complex_d erg;
  erg.re = c1.re + c2.re;
  erg.im = c1.im + c2.im;
return erg;
}


__device__ inline dev_complex_d dev_cdiv_d(dev_complex_d c1, dev_complex_d c2) { /* dividiert c1 durch c2 */
  dev_complex_d erg;
  double oneovernenner = 1.0/(c2.re*c2.re + c2.im*c2.im);
  erg.re = oneovernenner*(c1.re*c2.re + c1.im*c2.im);
  erg.im = oneovernenner*(c1.im*c2.re - c1.re*c2.im);
return erg;
}


__device__ inline dev_complex_d dev_csub_d(dev_complex_d c1, dev_complex_d c2){
   dev_complex_d erg;
   erg.re = c1.re - c2.re;
   erg.im = c1.im - c2.im;
return erg;
}


__device__ inline dev_complex_d dev_initcomplex_d(double re, double im){/* gibt komplexe Zahl mit Realt re und Imt im zurück*/
    dev_complex_d erg;
    erg.re = re;
    erg.im = im;
return (erg);
}

//********************* complex double operations **********************

 




// this is to copy a spinor in global mem
__device__ inline void dev_copy_spinor_d(dev_spinor_d *i1, dev_spinor_d *i2){
  int i;
  #pragma unroll 6
  for(i=0;i<12;i++){ //color + spin
    (*(i2+i*DEVOFF)).x = (*(i1+i*DEVOFF)).x;
    (*(i2+i*DEVOFF)).y = (*(i1+i*DEVOFF)).y;
  }
}


// this is to copy a spinor in local mem
__device__ inline void dev_copy_spinor_local_d(double4 *i1, double4 *i2){
  int i;
  #pragma unroll 6
  for(i=0;i<6;i++){ //color + spin
    (*(i2+i)).x = (*(i1+i)).x;
    (*(i2+i)).y = (*(i1+i)).y;
    (*(i2+i)).z = (*(i1+i)).z;
    (*(i2+i)).w = (*(i1+i)).w;
  }
}


// this is to store a spinor in local mem (24 numbers adjacent) to global mem
__device__ inline void dev_write_spinor_d(double4 *i1, dev_spinor_d *i2){
  int i;
  #pragma unroll 6
  for(i=0;i<6;i++){ //color + spin
    (*(i2+(2*i)*DEVOFF)).x = (*(i1+i)).x;
    (*(i2+(2*i)*DEVOFF)).y = (*(i1+i)).y;
    (*(i2+(2*i+1)*DEVOFF)).x = (*(i1+i)).z;
    (*(i2+(2*i+1)*DEVOFF)).y = (*(i1+i)).w;
  }
}


// this is to read a spinor from global into local mem (24 numbers adjacent)
__device__ inline void dev_read_spinor_d(double4 *i1, dev_spinor_d *i2){
  int i;
  #pragma unroll 6
  for(i=0;i<6;i++){ //color + spin
    (*(i1+i)).x  = (*(i2+(2*i)*DEVOFF)).x;
    (*(i1+i)).y  = (*(i2+(2*i)*DEVOFF)).y;
    (*(i1+i)).z  = (*(i2+(2*i+1)*DEVOFF)).x;
    (*(i1+i)).w  = (*(i2+(2*i+1)*DEVOFF)).y;
  }
}


// this is to zero a spinor in global mem
__device__ inline void dev_zero_spinor_d(dev_spinor_d *sin){
  int i;
  #pragma unroll 12
  for(i=0;i<12;i++){ //color + spin
    (*(sin+i*DEVOFF)).x = 0.0;
    (*(sin+i*DEVOFF)).y = 0.0;
  }
}


// this is to zero a spinor in local mem
__device__ inline void dev_zero_spinor_local_d(double4 *sin){
  int i;
  #pragma unroll 6
  for(i=0;i<6;i++){ //color + spin
    (*(sin+i)).x = 0.0;
    (*(sin+i)).y = 0.0;
    (*(sin+i)).z = 0.0;
    (*(sin+i)).w = 0.0;
  }
}







__global__ void dev_gamma5_d(dev_spinor_d * sin, dev_spinor_d * sout){
  int pos;
  pos= threadIdx.x + blockDim.x*blockIdx.x;  
  if(pos < dev_VOLUME){
          sout[pos+0*DEVOFF].x = sin[pos+0*DEVOFF].x;
          sout[pos+0*DEVOFF].y = sin[pos+0*DEVOFF].y;
          sout[pos+1*DEVOFF].x = sin[pos+1*DEVOFF].x;
          sout[pos+1*DEVOFF].y = sin[pos+1*DEVOFF].y;
          sout[pos+2*DEVOFF].x = sin[pos+2*DEVOFF].x;
          sout[pos+2*DEVOFF].y = sin[pos+2*DEVOFF].y;
          
          sout[pos+3*DEVOFF].x = sin[pos+3*DEVOFF].x;
          sout[pos+3*DEVOFF].y = sin[pos+3*DEVOFF].y;
          sout[pos+4*DEVOFF].x = sin[pos+4*DEVOFF].x;
          sout[pos+4*DEVOFF].y = sin[pos+4*DEVOFF].y;
          sout[pos+5*DEVOFF].x = sin[pos+5*DEVOFF].x;
          sout[pos+5*DEVOFF].y = sin[pos+5*DEVOFF].y;   
          
          sout[pos+6*DEVOFF].x = -1.0*sin[pos+6*DEVOFF].x;
          sout[pos+6*DEVOFF].y = -1.0*sin[pos+6*DEVOFF].y;
          sout[pos+7*DEVOFF].x = -1.0*sin[pos+7*DEVOFF].x;
          sout[pos+7*DEVOFF].y = -1.0*sin[pos+7*DEVOFF].y;
          sout[pos+8*DEVOFF].x = -1.0*sin[pos+8*DEVOFF].x;
          sout[pos+8*DEVOFF].y = -1.0*sin[pos+8*DEVOFF].y; 

          sout[pos+9*DEVOFF].x = -1.0*sin[pos+9*DEVOFF].x;
          sout[pos+9*DEVOFF].y = -1.0*sin[pos+9*DEVOFF].y;
          sout[pos+10*DEVOFF].x = -1.0*sin[pos+10*DEVOFF].x;
          sout[pos+10*DEVOFF].y = -1.0*sin[pos+10*DEVOFF].y;
          sout[pos+11*DEVOFF].x = -1.0*sin[pos+11*DEVOFF].x;
          sout[pos+11*DEVOFF].y = -1.0*sin[pos+11*DEVOFF].y;                 
  } 
}



/* OBSOLETE??
__device__ void inline dev_skalarmult_gamma5_spinor_d(dev_spinor_d * out, dev_complex_d lambda, dev_spinor_d * in){

 dev_spinor_d shelp, tempout;

shelp = *(in);
 tempout.x = shelp.x*lambda.re;
 tempout.x -= shelp.y*lambda.im;
 
 tempout.y = shelp.y*lambda.re;
 tempout.y += shelp.x*lambda.im;
 
 tempout.z = shelp.z*lambda.re;
 tempout.z -= shelp.w*lambda.im;

 tempout.w = shelp.w*lambda.re;
 tempout.w += shelp.z*lambda.im;
(*(out)) = tempout;
 
 
shelp = *(in+1);
 tempout.x = shelp.x*lambda.re;
 tempout.x -= shelp.y*lambda.im;
 
 tempout.y = shelp.y*lambda.re;
 tempout.y += shelp.x*lambda.im;
 
 tempout.z = shelp.z*lambda.re;
 tempout.z -= shelp.w*lambda.im;

 tempout.w = shelp.w*lambda.re;
 tempout.w += shelp.z*lambda.im;
(*(out+1)) = tempout; 
 
 
shelp = *(in+2);
 tempout.x = shelp.x*lambda.re;
 tempout.x -= shelp.y*lambda.im;

 tempout.y = shelp.y*lambda.re;
 tempout.y += shelp.x*lambda.im;
 
 tempout.z = shelp.z*lambda.re;
 tempout.z -= shelp.w*lambda.im;

 tempout.w = shelp.w*lambda.re;
 tempout.w += shelp.z*lambda.im;
(*(out+2)) = tempout; 

 
shelp = *(in+3);
 tempout.x = shelp.y*lambda.im;
 tempout.x -= shelp.x*lambda.re;

 tempout.y = - shelp.x*lambda.im;
 tempout.y -= shelp.y*lambda.re;
 
 tempout.z = shelp.w*lambda.im;
 tempout.z -= shelp.z*lambda.re;

 tempout.w = -shelp.z*lambda.im;
 tempout.w -= shelp.w*lambda.re;
(*(out+3)) = tempout;

shelp = *(in+4);
 tempout.x = shelp.y*lambda.im;
 tempout.x -= shelp.x*lambda.re;

 tempout.y = - shelp.x*lambda.im;
 tempout.y -= shelp.y*lambda.re;
 
 tempout.z = shelp.w*lambda.im;
 tempout.z -= shelp.z*lambda.re;

 tempout.w = -shelp.z*lambda.im;
 tempout.w -= shelp.w*lambda.re;
(*(out+4)) = tempout;

shelp = *(in+5);
 tempout.x = shelp.y*lambda.im;
 tempout.x -= shelp.x*lambda.re;

 tempout.y = - shelp.x*lambda.im;
 tempout.y -= shelp.y*lambda.re;
 
 tempout.z = shelp.w*lambda.im;
 tempout.z -= shelp.z*lambda.re;

 tempout.w = -shelp.z*lambda.im;
 tempout.w -= shelp.w*lambda.re;
(*(out+5)) = tempout;
}

*/




__device__ void inline dev_skalarmult_gamma5_globalspinor_d(double4 * out, dev_complex_d lambda, dev_spinor_d * in){

 dev_spinor_d shelp;
 double4 tempout;

shelp = *(in);
 tempout.x = shelp.x*lambda.re;
 tempout.x -= shelp.y*lambda.im;
 
 tempout.y = shelp.y*lambda.re;
 tempout.y += shelp.x*lambda.im;

shelp = *(in+1*DEVOFF); 
 tempout.z = shelp.x*lambda.re;
 tempout.z -= shelp.y*lambda.im;

 tempout.w = shelp.y*lambda.re;
 tempout.w += shelp.x*lambda.im;
 (*(out)) = tempout;
 
 
shelp = *(in+2*DEVOFF);
 tempout.x = shelp.x*lambda.re;
 tempout.x -= shelp.y*lambda.im;
 
 tempout.y = shelp.y*lambda.re;
 tempout.y += shelp.x*lambda.im;

shelp = *(in+3*DEVOFF); 
 tempout.z = shelp.x*lambda.re;
 tempout.z -= shelp.y*lambda.im;

 tempout.w = shelp.y*lambda.re;
 tempout.w += shelp.x*lambda.im;
(*(out+1)) = tempout; 
 
 
shelp = *(in+4*DEVOFF);
 tempout.x = shelp.x*lambda.re;
 tempout.x -= shelp.y*lambda.im;

 tempout.y = shelp.y*lambda.re;
 tempout.y += shelp.x*lambda.im;

shelp = *(in+5*DEVOFF); 
 tempout.z = shelp.x*lambda.re;
 tempout.z -= shelp.y*lambda.im;

 tempout.w = shelp.y*lambda.re;
 tempout.w += shelp.x*lambda.im;
(*(out+2)) = tempout; 

 
shelp = *(in+6*DEVOFF);
 tempout.x = shelp.y*lambda.im;
 tempout.x -= shelp.x*lambda.re;

 tempout.y = - shelp.x*lambda.im;
 tempout.y -= shelp.y*lambda.re;

shelp = *(in+7*DEVOFF); 
 tempout.z = shelp.y*lambda.im;
 tempout.z -= shelp.x*lambda.re;

 tempout.w = -shelp.x*lambda.im;
 tempout.w -= shelp.y*lambda.re;
(*(out+3)) = tempout;


shelp = *(in+8*DEVOFF);
 tempout.x = shelp.y*lambda.im;
 tempout.x -= shelp.x*lambda.re;

 tempout.y = - shelp.x*lambda.im;
 tempout.y -= shelp.y*lambda.re;

shelp = *(in+9*DEVOFF); 
 tempout.z = shelp.y*lambda.im;
 tempout.z -= shelp.x*lambda.re;

 tempout.w = -shelp.x*lambda.im;
 tempout.w -= shelp.y*lambda.re;
(*(out+4)) = tempout;


shelp = *(in+10*DEVOFF);
 tempout.x = shelp.y*lambda.im;
 tempout.x -= shelp.x*lambda.re;

 tempout.y = - shelp.x*lambda.im;
 tempout.y -= shelp.y*lambda.re;

shelp = *(in+11*DEVOFF); 
 tempout.z = shelp.y*lambda.im;
 tempout.z -= shelp.x*lambda.re;

 tempout.w = -shelp.x*lambda.im;
 tempout.w -= shelp.y*lambda.re;
(*(out+5)) = tempout;
}








/* OBSOLETE??
__device__ void inline dev_realmult_spinor_d(dev_spinor_d * in, double lambda){
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


__device__ void inline dev_realmult_spinor_assign_d(dev_spinor_d* out, double lambda, dev_spinor_d* in){
int i;
#pragma unroll 6
  for(i=0;i<6;i++){ //color + spin
  //out[i] = in[i]*lambda;
      (*(out+i)).x = (*(in+i)).x*lambda;
	  (*(out+i)).y = (*(in+i)).y*lambda;
      
	  (*(out+i)).z = (*(in+i)).z*lambda;
      (*(out+i)).w = (*(in+i)).w*lambda;
  }
}
*/



__device__ void inline dev_realmult_spinor_assigntoglobal_d(dev_spinor_d* out, double lambda, double4* in){
int i;
#pragma unroll 6
  for(i=0;i<6;i++){ //color + spin
  //out[i] = in[i]*lambda;
      (*(out+(2*i)*DEVOFF)).x = (*(in+i)).x*lambda;
	  (*(out+(2*i)*DEVOFF)).y = (*(in+i)).y*lambda;
      
	  (*(out+(2*i+1)*DEVOFF)).x = (*(in+i)).z*lambda;
      (*(out+(2*i+1)*DEVOFF)).y = (*(in+i)).w*lambda;
  }
}


/* OBSOLETE??
__device__ void dev_assign_realmult_add_spinor_d(dev_spinor_d* out, double lambda, dev_spinor_d* in1,  dev_spinor_d* in2){
int i;
double help;
//out = lambda*(in1 + in2)
#pragma unroll 6
  for(i=0;i<6;i++){ //color + spin

      help = (*(in1+i)).x*lambda;
      help += (*(in2+i)).x*lambda;
      (*(out+i)).x = help;
      
      help = (*(in1+i)).y*lambda;
      help += (*(in2+i)).y*lambda;
      (*(out+i)).y = help;      

      help = (*(in1+i)).z*lambda;
      help += (*(in2+i)).z*lambda;
      (*(out+i)).z = help;

      help = (*(in1+i)).w*lambda;
      help += (*(in2+i)).w*lambda;
      (*(out+i)).w = help;
  } 
}

*/

__device__ inline void dev_add_spinor_assign_d(double4 * i1, double4 * i2){
  int i;
  #pragma unroll 6
  for(i=0;i<6;i++){ //color + spin
    (*(i1+i)).x = (*(i1+i)).x + (*(i2+i)).x;
    (*(i1+i)).y = (*(i1+i)).y + (*(i2+i)).y;
    (*(i1+i)).z = (*(i1+i)).z + (*(i2+i)).z;
    (*(i1+i)).w = (*(i1+i)).w + (*(i2+i)).w;
  }
}



__device__ inline void dev_add_globalspinor_assign_d(double4 * i1, dev_spinor_d * i2){
  int i;
  #pragma unroll 6
  for(i=0;i<6;i++){ //color + spin
    (*(i1+i)).x = (*(i1+i)).x + (*(i2+(2*i)*DEVOFF)).x;
    (*(i1+i)).y = (*(i1+i)).y + (*(i2+(2*i)*DEVOFF)).y;
    (*(i1+i)).z = (*(i1+i)).z + (*(i2+(2*i+1)*DEVOFF)).x;
    (*(i1+i)).w = (*(i1+i)).w + (*(i2+(2*i+1)*DEVOFF)).y;
  }
}



__device__ inline void dev_sub_spinor_assign_d(double4 * i1, double4 * i2){
  int i;
  #pragma unroll 6
  for(i=0;i<6;i++){ //color + spin
    (*(i1+i)).x = (*(i1+i)).x - (*(i2+i)).x;
    (*(i1+i)).y = (*(i1+i)).y - (*(i2+i)).y;
    (*(i1+i)).z = (*(i1+i)).z - (*(i2+i)).z;
    (*(i1+i)).w = (*(i1+i)).w - (*(i2+i)).w;
  }
}

__device__ inline void dev_sub_globalspinor_assign_d(double4 * i1, dev_spinor_d * i2){
  int i;
  #pragma unroll 6
  for(i=0;i<6;i++){ //color + spin
    (*(i1+i)).x = (*(i1+i)).x - (*(i2+(2*i)*DEVOFF)).x;
    (*(i1+i)).y = (*(i1+i)).y - (*(i2+(2*i)*DEVOFF)).y;
    (*(i1+i)).z = (*(i1+i)).z - (*(i2+(2*i+1)*DEVOFF)).x;
    (*(i1+i)).w = (*(i1+i)).w - (*(i2+(2*i+1)*DEVOFF)).y;
  }
}





//multipliziert su3-Matrix mal Spinor im Dirac-Raum
//code in su3_MtV.txt -- generated with codegen
__device__ void dev_su3MtV_d(dev_su3_d M, const dev_spinor_d * s, double4 * out){

(*(out+0)).x =  M[0][0].re*(*(s+0*DEVOFF)).x;
(*(out+0)).x -= M[0][0].im*(*(s+0*DEVOFF)).y;
(*(out+0)).x += M[0][1].re*(*(s+1*DEVOFF)).x;
(*(out+0)).x -= M[0][1].im*(*(s+1*DEVOFF)).y;
(*(out+0)).x += M[0][2].re*(*(s+2*DEVOFF)).x;
(*(out+0)).x -= M[0][2].im*(*(s+2*DEVOFF)).y;

(*(out+0)).y =  M[0][0].re*(*(s+0*DEVOFF)).y;
(*(out+0)).y += M[0][0].im*(*(s+0*DEVOFF)).x;
(*(out+0)).y += M[0][1].re*(*(s+1*DEVOFF)).y;
(*(out+0)).y += M[0][1].im*(*(s+1*DEVOFF)).x;
(*(out+0)).y += M[0][2].re*(*(s+2*DEVOFF)).y;
(*(out+0)).y += M[0][2].im*(*(s+2*DEVOFF)).x;


(*(out+0)).z =  M[1][0].re*(*(s+0*DEVOFF)).x;
(*(out+0)).z -= M[1][0].im*(*(s+0*DEVOFF)).y;
(*(out+0)).z += M[1][1].re*(*(s+1*DEVOFF)).x;
(*(out+0)).z -= M[1][1].im*(*(s+1*DEVOFF)).y;
(*(out+0)).z += M[1][2].re*(*(s+2*DEVOFF)).x;
(*(out+0)).z -= M[1][2].im*(*(s+2*DEVOFF)).y;

(*(out+0)).w =  M[1][0].re*(*(s+0*DEVOFF)).y;
(*(out+0)).w += M[1][0].im*(*(s+0*DEVOFF)).x;
(*(out+0)).w += M[1][1].re*(*(s+1*DEVOFF)).y;
(*(out+0)).w += M[1][1].im*(*(s+1*DEVOFF)).x;
(*(out+0)).w += M[1][2].re*(*(s+2*DEVOFF)).y;
(*(out+0)).w += M[1][2].im*(*(s+2*DEVOFF)).x;


(*(out+1)).x =  M[2][0].re*(*(s+0*DEVOFF)).x;
(*(out+1)).x -= M[2][0].im*(*(s+0*DEVOFF)).y;
(*(out+1)).x += M[2][1].re*(*(s+1*DEVOFF)).x;
(*(out+1)).x -= M[2][1].im*(*(s+1*DEVOFF)).y;
(*(out+1)).x += M[2][2].re*(*(s+2*DEVOFF)).x;
(*(out+1)).x -= M[2][2].im*(*(s+2*DEVOFF)).y;

(*(out+1)).y =  M[2][0].re*(*(s+0*DEVOFF)).y;
(*(out+1)).y += M[2][0].im*(*(s+0*DEVOFF)).x;
(*(out+1)).y += M[2][1].re*(*(s+1*DEVOFF)).y;
(*(out+1)).y += M[2][1].im*(*(s+1*DEVOFF)).x;
(*(out+1)).y += M[2][2].re*(*(s+2*DEVOFF)).y;
(*(out+1)).y += M[2][2].im*(*(s+2*DEVOFF)).x;


(*(out+1)).z =  M[0][0].re*(*(s+3*DEVOFF)).x;
(*(out+1)).z -= M[0][0].im*(*(s+3*DEVOFF)).y;
(*(out+1)).z += M[0][1].re*(*(s+4*DEVOFF)).x;
(*(out+1)).z -= M[0][1].im*(*(s+4*DEVOFF)).y;
(*(out+1)).z += M[0][2].re*(*(s+5*DEVOFF)).x;
(*(out+1)).z -= M[0][2].im*(*(s+5*DEVOFF)).y;

(*(out+1)).w =  M[0][0].re*(*(s+3*DEVOFF)).y;
(*(out+1)).w += M[0][0].im*(*(s+3*DEVOFF)).x;
(*(out+1)).w += M[0][1].re*(*(s+4*DEVOFF)).y;
(*(out+1)).w += M[0][1].im*(*(s+4*DEVOFF)).x;
(*(out+1)).w += M[0][2].re*(*(s+5*DEVOFF)).y;
(*(out+1)).w += M[0][2].im*(*(s+5*DEVOFF)).x;


(*(out+2)).x =  M[1][0].re*(*(s+3*DEVOFF)).x;
(*(out+2)).x -= M[1][0].im*(*(s+3*DEVOFF)).y;
(*(out+2)).x += M[1][1].re*(*(s+4*DEVOFF)).x;
(*(out+2)).x -= M[1][1].im*(*(s+4*DEVOFF)).y;
(*(out+2)).x += M[1][2].re*(*(s+5*DEVOFF)).x;
(*(out+2)).x -= M[1][2].im*(*(s+5*DEVOFF)).y;

(*(out+2)).y =  M[1][0].re*(*(s+3*DEVOFF)).y;
(*(out+2)).y += M[1][0].im*(*(s+3*DEVOFF)).x;
(*(out+2)).y += M[1][1].re*(*(s+4*DEVOFF)).y;
(*(out+2)).y += M[1][1].im*(*(s+4*DEVOFF)).x;
(*(out+2)).y += M[1][2].re*(*(s+5*DEVOFF)).y;
(*(out+2)).y += M[1][2].im*(*(s+5*DEVOFF)).x;


(*(out+2)).z =  M[2][0].re*(*(s+3*DEVOFF)).x;
(*(out+2)).z -= M[2][0].im*(*(s+3*DEVOFF)).y;
(*(out+2)).z += M[2][1].re*(*(s+4*DEVOFF)).x;
(*(out+2)).z -= M[2][1].im*(*(s+4*DEVOFF)).y;
(*(out+2)).z += M[2][2].re*(*(s+5*DEVOFF)).x;
(*(out+2)).z -= M[2][2].im*(*(s+5*DEVOFF)).y;

(*(out+2)).w =  M[2][0].re*(*(s+3*DEVOFF)).y;
(*(out+2)).w += M[2][0].im*(*(s+3*DEVOFF)).x;
(*(out+2)).w += M[2][1].re*(*(s+4*DEVOFF)).y;
(*(out+2)).w += M[2][1].im*(*(s+4*DEVOFF)).x;
(*(out+2)).w += M[2][2].re*(*(s+5*DEVOFF)).y;
(*(out+2)).w += M[2][2].im*(*(s+5*DEVOFF)).x;


(*(out+3)).x =  M[0][0].re*(*(s+6*DEVOFF)).x;
(*(out+3)).x -= M[0][0].im*(*(s+6*DEVOFF)).y;
(*(out+3)).x += M[0][1].re*(*(s+7*DEVOFF)).x;
(*(out+3)).x -= M[0][1].im*(*(s+7*DEVOFF)).y;
(*(out+3)).x += M[0][2].re*(*(s+8*DEVOFF)).x;
(*(out+3)).x -= M[0][2].im*(*(s+8*DEVOFF)).y;

(*(out+3)).y =  M[0][0].re*(*(s+6*DEVOFF)).y;
(*(out+3)).y += M[0][0].im*(*(s+6*DEVOFF)).x;
(*(out+3)).y += M[0][1].re*(*(s+7*DEVOFF)).y;
(*(out+3)).y += M[0][1].im*(*(s+7*DEVOFF)).x;
(*(out+3)).y += M[0][2].re*(*(s+8*DEVOFF)).y;
(*(out+3)).y += M[0][2].im*(*(s+8*DEVOFF)).x;


(*(out+3)).z =  M[1][0].re*(*(s+6*DEVOFF)).x;
(*(out+3)).z -= M[1][0].im*(*(s+6*DEVOFF)).y;
(*(out+3)).z += M[1][1].re*(*(s+7*DEVOFF)).x;
(*(out+3)).z -= M[1][1].im*(*(s+7*DEVOFF)).y;
(*(out+3)).z += M[1][2].re*(*(s+8*DEVOFF)).x;
(*(out+3)).z -= M[1][2].im*(*(s+8*DEVOFF)).y;

(*(out+3)).w =  M[1][0].re*(*(s+6*DEVOFF)).y;
(*(out+3)).w += M[1][0].im*(*(s+6*DEVOFF)).x;
(*(out+3)).w += M[1][1].re*(*(s+7*DEVOFF)).y;
(*(out+3)).w += M[1][1].im*(*(s+7*DEVOFF)).x;
(*(out+3)).w += M[1][2].re*(*(s+8*DEVOFF)).y;
(*(out+3)).w += M[1][2].im*(*(s+8*DEVOFF)).x;


(*(out+4)).x =  M[2][0].re*(*(s+6*DEVOFF)).x;
(*(out+4)).x -= M[2][0].im*(*(s+6*DEVOFF)).y;
(*(out+4)).x += M[2][1].re*(*(s+7*DEVOFF)).x;
(*(out+4)).x -= M[2][1].im*(*(s+7*DEVOFF)).y;
(*(out+4)).x += M[2][2].re*(*(s+8*DEVOFF)).x;
(*(out+4)).x -= M[2][2].im*(*(s+8*DEVOFF)).y;

(*(out+4)).y =  M[2][0].re*(*(s+6*DEVOFF)).y;
(*(out+4)).y += M[2][0].im*(*(s+6*DEVOFF)).x;
(*(out+4)).y += M[2][1].re*(*(s+7*DEVOFF)).y;
(*(out+4)).y += M[2][1].im*(*(s+7*DEVOFF)).x;
(*(out+4)).y += M[2][2].re*(*(s+8*DEVOFF)).y;
(*(out+4)).y += M[2][2].im*(*(s+8*DEVOFF)).x;


(*(out+4)).z =  M[0][0].re*(*(s+9*DEVOFF)).x;
(*(out+4)).z -= M[0][0].im*(*(s+9*DEVOFF)).y;
(*(out+4)).z += M[0][1].re*(*(s+10*DEVOFF)).x;
(*(out+4)).z -= M[0][1].im*(*(s+10*DEVOFF)).y;
(*(out+4)).z += M[0][2].re*(*(s+11*DEVOFF)).x;
(*(out+4)).z -= M[0][2].im*(*(s+11*DEVOFF)).y;

(*(out+4)).w =  M[0][0].re*(*(s+9*DEVOFF)).y;
(*(out+4)).w += M[0][0].im*(*(s+9*DEVOFF)).x;
(*(out+4)).w += M[0][1].re*(*(s+10*DEVOFF)).y;
(*(out+4)).w += M[0][1].im*(*(s+10*DEVOFF)).x;
(*(out+4)).w += M[0][2].re*(*(s+11*DEVOFF)).y;
(*(out+4)).w += M[0][2].im*(*(s+11*DEVOFF)).x;


(*(out+5)).x =  M[1][0].re*(*(s+9*DEVOFF)).x;
(*(out+5)).x -= M[1][0].im*(*(s+9*DEVOFF)).y;
(*(out+5)).x += M[1][1].re*(*(s+10*DEVOFF)).x;
(*(out+5)).x -= M[1][1].im*(*(s+10*DEVOFF)).y;
(*(out+5)).x += M[1][2].re*(*(s+11*DEVOFF)).x;
(*(out+5)).x -= M[1][2].im*(*(s+11*DEVOFF)).y;

(*(out+5)).y =  M[1][0].re*(*(s+9*DEVOFF)).y;
(*(out+5)).y += M[1][0].im*(*(s+9*DEVOFF)).x;
(*(out+5)).y += M[1][1].re*(*(s+10*DEVOFF)).y;
(*(out+5)).y += M[1][1].im*(*(s+10*DEVOFF)).x;
(*(out+5)).y += M[1][2].re*(*(s+11*DEVOFF)).y;
(*(out+5)).y += M[1][2].im*(*(s+11*DEVOFF)).x;


(*(out+5)).z =  M[2][0].re*(*(s+9*DEVOFF)).x;
(*(out+5)).z -= M[2][0].im*(*(s+9*DEVOFF)).y;
(*(out+5)).z += M[2][1].re*(*(s+10*DEVOFF)).x;
(*(out+5)).z -= M[2][1].im*(*(s+10*DEVOFF)).y;
(*(out+5)).z += M[2][2].re*(*(s+11*DEVOFF)).x;
(*(out+5)).z -= M[2][2].im*(*(s+11*DEVOFF)).y;

(*(out+5)).w =  M[2][0].re*(*(s+9*DEVOFF)).y;
(*(out+5)).w += M[2][0].im*(*(s+9*DEVOFF)).x;
(*(out+5)).w += M[2][1].re*(*(s+10*DEVOFF)).y;
(*(out+5)).w += M[2][1].im*(*(s+10*DEVOFF)).x;
(*(out+5)).w += M[2][2].re*(*(s+11*DEVOFF)).y;
(*(out+5)).w += M[2][2].im*(*(s+11*DEVOFF)).x;
}






//-kappa(r - gamma_mu)
__device__ void dev_kappaP0_plus_d(double4 * out, double4 * in, dev_complex_d kappa){


     (*(out+0)).x -= (*(in+0)).x*kappa.re;
     (*(out+0)).x += (*(in+0)).y*kappa.im;     
     (*(out+0)).y -= (*(in+0)).y*kappa.re;
     (*(out+0)).y -= (*(in+0)).x*kappa.im;     
     (*(out+0)).x -= (*(in+3)).x*kappa.re;
     (*(out+0)).x += (*(in+3)).y*kappa.im;     
     (*(out+0)).y -= (*(in+3)).y*kappa.re;
     (*(out+0)).y -= (*(in+3)).x*kappa.im;

     
     (*(out+3)).x -= (*(in+3)).x*kappa.re;
     (*(out+3)).x += (*(in+3)).y*kappa.im;     
     (*(out+3)).y -= (*(in+3)).y*kappa.re;
     (*(out+3)).y -= (*(in+3)).x*kappa.im;          
     (*(out+3)).x -= (*(in+0)).x*kappa.re;
     (*(out+3)).x += (*(in+0)).y*kappa.im;     
     (*(out+3)).y -= (*(in+0)).y*kappa.re; 
     (*(out+3)).y -= (*(in+0)).x*kappa.im;


     (*(out+0)).z -= (*(in+0)).z*kappa.re;
     (*(out+0)).z += (*(in+0)).w*kappa.im;     
     (*(out+0)).w -= (*(in+0)).w*kappa.re;
     (*(out+0)).w -= (*(in+0)).z*kappa.im;     
     (*(out+0)).z -= (*(in+3)).z*kappa.re;
     (*(out+0)).z += (*(in+3)).w*kappa.im;          
     (*(out+0)).w -= (*(in+3)).w*kappa.re; 
     (*(out+0)).w -= (*(in+3)).z*kappa.im;     

     
     (*(out+3)).z -= (*(in+3)).z*kappa.re;
     (*(out+3)).z += (*(in+3)).w*kappa.im;     
     (*(out+3)).w -= (*(in+3)).w*kappa.re;
     (*(out+3)).w -= (*(in+3)).z*kappa.im;     
     (*(out+3)).z -= (*(in+0)).z*kappa.re;
     (*(out+3)).z += (*(in+0)).w*kappa.im;     
     (*(out+3)).w -= (*(in+0)).w*kappa.re;
     (*(out+3)).w -= (*(in+0)).z*kappa.im;
 
  
     (*(out+1)).x -= (*(in+1)).x*kappa.re;
     (*(out+1)).x += (*(in+1)).y*kappa.im;     
     (*(out+1)).y -= (*(in+1)).y*kappa.re;
     (*(out+1)).y -= (*(in+1)).x*kappa.im;     
     (*(out+1)).x -= (*(in+4)).x*kappa.re;
     (*(out+1)).x += (*(in+4)).y*kappa.im;     
     (*(out+1)).y -= (*(in+4)).y*kappa.re;     
     (*(out+1)).y -= (*(in+4)).x*kappa.im; 
     
     
     (*(out+4)).x -= (*(in+4)).x*kappa.re;
     (*(out+4)).x += (*(in+4)).y*kappa.im;     
     (*(out+4)).y -= (*(in+4)).y*kappa.re; 
     (*(out+4)).y -= (*(in+4)).x*kappa.im;     
     (*(out+4)).x -= (*(in+1)).x*kappa.re;
     (*(out+4)).x += (*(in+1)).y*kappa.im;     
     (*(out+4)).y -= (*(in+1)).y*kappa.re; 
     (*(out+4)).y -= (*(in+1)).x*kappa.im; 
 
 
     (*(out+1)).z -= (*(in+1)).z*kappa.re;
     (*(out+1)).z += (*(in+1)).w*kappa.im;     
     (*(out+1)).w -= (*(in+1)).w*kappa.re;
     (*(out+1)).w -= (*(in+1)).z*kappa.im;     
     (*(out+1)).z -= (*(in+4)).z*kappa.re;
     (*(out+1)).z += (*(in+4)).w*kappa.im;     
     (*(out+1)).w -= (*(in+4)).w*kappa.re; 
     (*(out+1)).w -= (*(in+4)).z*kappa.im;    
     
     
     (*(out+4)).z -= (*(in+4)).z*kappa.re;
     (*(out+4)).z += (*(in+4)).w*kappa.im;     
     (*(out+4)).w -= (*(in+4)).w*kappa.re;
     (*(out+4)).w -= (*(in+4)).z*kappa.im;     
     (*(out+4)).z -= (*(in+1)).z*kappa.re;
     (*(out+4)).z += (*(in+1)).w*kappa.im;     
     (*(out+4)).w -= (*(in+1)).w*kappa.re;      
     (*(out+4)).w -= (*(in+1)).z*kappa.im; 

     
     (*(out+2)).x -= (*(in+2)).x*kappa.re;
     (*(out+2)).x += (*(in+2)).y*kappa.im;     
     (*(out+2)).y -= (*(in+2)).y*kappa.re;
     (*(out+2)).y -= (*(in+2)).x*kappa.im;     
     (*(out+2)).x -= (*(in+5)).x*kappa.re;
     (*(out+2)).x += (*(in+5)).y*kappa.im;     
     (*(out+2)).y -= (*(in+5)).y*kappa.re;
     (*(out+2)).y -= (*(in+5)).x*kappa.im;
     
     
     (*(out+5)).x -= (*(in+5)).x*kappa.re;
     (*(out+5)).x += (*(in+5)).y*kappa.im;     
     (*(out+5)).y -= (*(in+5)).y*kappa.re;
     (*(out+5)).y -= (*(in+5)).x*kappa.im;     
     (*(out+5)).x -= (*(in+2)).x*kappa.re;
     (*(out+5)).x += (*(in+2)).y*kappa.im;     
     (*(out+5)).y -= (*(in+2)).y*kappa.re; 
     (*(out+5)).y -= (*(in+2)).x*kappa.im; 
     
           
     (*(out+2)).z -= (*(in+2)).z*kappa.re;
     (*(out+2)).z += (*(in+2)).w*kappa.im;     
     (*(out+2)).w -= (*(in+2)).w*kappa.re;
     (*(out+2)).w -= (*(in+2)).z*kappa.im;     
     (*(out+2)).z -= (*(in+5)).z*kappa.re;
     (*(out+2)).z += (*(in+5)).w*kappa.im;     
     (*(out+2)).w -= (*(in+5)).w*kappa.re;
     (*(out+2)).w -= (*(in+5)).z*kappa.im;
     
               
     (*(out+5)).z -= (*(in+5)).z*kappa.re;
     (*(out+5)).z += (*(in+5)).w*kappa.im;     
     (*(out+5)).w -= (*(in+5)).w*kappa.re;
     (*(out+5)).w -= (*(in+5)).z*kappa.im;     
     (*(out+5)).z -= (*(in+2)).z*kappa.re;
     (*(out+5)).z += (*(in+2)).w*kappa.im;     
     (*(out+5)).w -= (*(in+2)).w*kappa.re;
     (*(out+5)).w -= (*(in+2)).z*kappa.im;     
  
}






//-kappa(r + gamma_mu)
__device__ void dev_kappaP0_minus_d(double4 * out, double4 * in, dev_complex_d kappa){

/*
 
     (*(out+0)).x -= (*(in+0)).x*kappa.re;
     (*(out+0)).x += (*(in+0)).y*kappa.im;     
     (*(out+0)).y -= (*(in+0)).y*kappa.re;
     (*(out+0)).y -= (*(in+0)).x*kappa.im;     
     (*(out+0)).x += (*(in+3)).x*kappa.re;
     (*(out+0)).x -= (*(in+3)).y*kappa.im;     
     (*(out+0)).y += (*(in+3)).y*kappa.re;
     (*(out+0)).y += (*(in+3)).x*kappa.im;     

     
     (*(out+3)).x -= (*(in+3)).x*kappa.re;
     (*(out+3)).x += (*(in+3)).y*kappa.im;     
     (*(out+3)).y -= (*(in+3)).y*kappa.re;
     (*(out+3)).y -= (*(in+3)).x*kappa.im;     
     (*(out+3)).x += (*(in+0)).x*kappa.re;
     (*(out+3)).x -= (*(in+0)).y*kappa.im;     
     (*(out+3)).y += (*(in+0)).y*kappa.re; 
     (*(out+3)).y += (*(in+0)).x*kappa.im; 


     (*(out+0)).z -= (*(in+0)).z*kappa.re;
     (*(out+0)).z += (*(in+0)).w*kappa.im;     
     (*(out+0)).w -= (*(in+0)).w*kappa.re;
     (*(out+0)).w -= (*(in+0)).z*kappa.im;     
     (*(out+0)).z += (*(in+3)).z*kappa.re;
     (*(out+0)).z -= (*(in+3)).w*kappa.im;     
     (*(out+0)).w += (*(in+3)).w*kappa.re;
     (*(out+0)).w += (*(in+3)).z*kappa.im;
     
          
     (*(out+3)).z -= (*(in+3)).z*kappa.re;
     (*(out+3)).z += (*(in+3)).w*kappa.im;     
     (*(out+3)).w -= (*(in+3)).w*kappa.re;
     (*(out+3)).w -= (*(in+3)).z*kappa.im;     
     (*(out+3)).z += (*(in+0)).z*kappa.re;
     (*(out+3)).z -= (*(in+0)).w*kappa.im;     
     (*(out+3)).w += (*(in+0)).w*kappa.re;
     (*(out+3)).w += (*(in+0)).z*kappa.im;
 
 
     (*(out+1)).x -= (*(in+1)).x*kappa.re;
     (*(out+1)).x += (*(in+1)).y*kappa.im;
     (*(out+1)).y -= (*(in+1)).y*kappa.re;
     (*(out+1)).y -= (*(in+1)).x*kappa.im;     
     (*(out+1)).x += (*(in+4)).x*kappa.re;
     (*(out+1)).x -= (*(in+4)).y*kappa.im;     
     (*(out+1)).y += (*(in+4)).y*kappa.re;
     (*(out+1)).y += (*(in+4)).x*kappa.im;     

     
     (*(out+4)).x -= (*(in+4)).x*kappa.re;
     (*(out+4)).x += (*(in+4)).y*kappa.im;     
     (*(out+4)).y -= (*(in+4)).y*kappa.re; 
     (*(out+4)).y -= (*(in+4)).x*kappa.im;     
     (*(out+4)).x += (*(in+1)).x*kappa.re;
     (*(out+4)).x -= (*(in+1)).y*kappa.im;     
     (*(out+4)).y += (*(in+1)).y*kappa.re; 
     (*(out+4)).y += (*(in+1)).x*kappa.im; 
 
 
     (*(out+1)).z -= (*(in+1)).z*kappa.re;
     (*(out+1)).z += (*(in+1)).w*kappa.im;     
     (*(out+1)).w -= (*(in+1)).w*kappa.re;
     (*(out+1)).w -= (*(in+1)).z*kappa.im;     
     (*(out+1)).z += (*(in+4)).z*kappa.re;
     (*(out+1)).z -= (*(in+4)).w*kappa.im;     
     (*(out+1)).w += (*(in+4)).w*kappa.re;
     (*(out+1)).w += (*(in+4)).z*kappa.im;     

     
     (*(out+4)).z -= (*(in+4)).z*kappa.re;
     (*(out+4)).z += (*(in+4)).w*kappa.im;     
     (*(out+4)).w -= (*(in+4)).w*kappa.re;
     (*(out+4)).w -= (*(in+4)).z*kappa.im;     
     (*(out+4)).z += (*(in+1)).z*kappa.re;
     (*(out+4)).z -= (*(in+1)).w*kappa.im;     
     (*(out+4)).w += (*(in+1)).w*kappa.re;      
     (*(out+4)).w += (*(in+1)).z*kappa.im; 

     
     (*(out+2)).x -= (*(in+2)).x*kappa.re;
     (*(out+2)).x += (*(in+2)).y*kappa.im;     
     (*(out+2)).y -= (*(in+2)).y*kappa.re;
     (*(out+2)).y -= (*(in+2)).x*kappa.im;     
     (*(out+2)).x += (*(in+5)).x*kappa.re;
     (*(out+2)).x -= (*(in+5)).y*kappa.im;     
     (*(out+2)).y += (*(in+5)).y*kappa.re;
     (*(out+2)).y += (*(in+5)).x*kappa.im;

     
     (*(out+5)).x -= (*(in+5)).x*kappa.re;
     (*(out+5)).x += (*(in+5)).y*kappa.im;     
     (*(out+5)).y -= (*(in+5)).y*kappa.re; 
     (*(out+5)).y -= (*(in+5)).x*kappa.im;     
     (*(out+5)).x += (*(in+2)).x*kappa.re;
     (*(out+5)).x -= (*(in+2)).y*kappa.im;     
     (*(out+5)).y += (*(in+2)).y*kappa.re; 
     (*(out+5)).y += (*(in+2)).x*kappa.im;   
   
       
     (*(out+2)).z -= (*(in+2)).z*kappa.re;
     (*(out+2)).z += (*(in+2)).w*kappa.im;     
     (*(out+2)).w -= (*(in+2)).w*kappa.re;
     (*(out+2)).w -= (*(in+2)).z*kappa.im;     
     (*(out+2)).z += (*(in+5)).z*kappa.re;
     (*(out+2)).z -= (*(in+5)).w*kappa.im;     
     (*(out+2)).w += (*(in+5)).w*kappa.re;
     (*(out+2)).w += (*(in+5)).z*kappa.im;
     
     
     (*(out+5)).z -= (*(in+5)).z*kappa.re;
     (*(out+5)).z += (*(in+5)).w*kappa.im;     
     (*(out+5)).w -= (*(in+5)).w*kappa.re;
     (*(out+5)).w -= (*(in+5)).z*kappa.im;     
     (*(out+5)).z += (*(in+2)).z*kappa.re;
     (*(out+5)).z -= (*(in+2)).w*kappa.im;      
     (*(out+5)).w += (*(in+2)).w*kappa.re; 
     (*(out+5)).w += (*(in+2)).z*kappa.im;
      
*/
     (*(out+0)).x -= (*(in+0)).x*kappa.re - (*(in+0)).y*kappa.im;
     (*(out+0)).y -= (*(in+0)).y*kappa.re + (*(in+0)).x*kappa.im;
     (*(out+0)).x += (*(in+3)).x*kappa.re - (*(in+3)).y*kappa.im;
     (*(out+0)).y += (*(in+3)).y*kappa.re + (*(in+3)).x*kappa.im;
     
     (*(out+3)).x -= (*(in+3)).x*kappa.re - (*(in+3)).y*kappa.im;
     (*(out+3)).y -= (*(in+3)).y*kappa.re + (*(in+3)).x*kappa.im;   
     (*(out+3)).x += (*(in+0)).x*kappa.re - (*(in+0)).y*kappa.im;
     (*(out+3)).y += (*(in+0)).y*kappa.re + (*(in+0)).x*kappa.im; 



     (*(out+0)).z -= (*(in+0)).z*kappa.re - (*(in+0)).w*kappa.im;
     (*(out+0)).w -= (*(in+0)).w*kappa.re + (*(in+0)).z*kappa.im;
     (*(out+0)).z += (*(in+3)).z*kappa.re - (*(in+3)).w*kappa.im;
     (*(out+0)).w += (*(in+3)).w*kappa.re + (*(in+3)).z*kappa.im;     
          
     (*(out+3)).z -= (*(in+3)).z*kappa.re - (*(in+3)).w*kappa.im;
     (*(out+3)).w -= (*(in+3)).w*kappa.re + (*(in+3)).z*kappa.im;
     (*(out+3)).z += (*(in+0)).z*kappa.re - (*(in+0)).w*kappa.im;
     (*(out+3)).w += (*(in+0)).w*kappa.re + (*(in+0)).z*kappa.im;

 
 
     (*(out+1)).x -= (*(in+1)).x*kappa.re - (*(in+1)).y*kappa.im;
     (*(out+1)).y -= (*(in+1)).y*kappa.re + (*(in+1)).x*kappa.im;
     (*(out+1)).x += (*(in+4)).x*kappa.re - (*(in+4)).y*kappa.im;
     (*(out+1)).y += (*(in+4)).y*kappa.re + (*(in+4)).x*kappa.im;     
     
     (*(out+4)).x -= (*(in+4)).x*kappa.re - (*(in+4)).y*kappa.im;
     (*(out+4)).y -= (*(in+4)).y*kappa.re + (*(in+4)).x*kappa.im;    
     (*(out+4)).x += (*(in+1)).x*kappa.re - (*(in+1)).y*kappa.im;
     (*(out+4)).y += (*(in+1)).y*kappa.re + (*(in+1)).x*kappa.im; 
 
 
 
     (*(out+1)).z -= (*(in+1)).z*kappa.re - (*(in+1)).w*kappa.im;
     (*(out+1)).w -= (*(in+1)).w*kappa.re + (*(in+1)).z*kappa.im;
     (*(out+1)).z += (*(in+4)).z*kappa.re - (*(in+4)).w*kappa.im;
     (*(out+1)).w += (*(in+4)).w*kappa.re + (*(in+4)).z*kappa.im;     
     
     (*(out+4)).z -= (*(in+4)).z*kappa.re - (*(in+4)).w*kappa.im;
     (*(out+4)).w -= (*(in+4)).w*kappa.re + (*(in+4)).z*kappa.im;    
     (*(out+4)).z += (*(in+1)).z*kappa.re - (*(in+1)).w*kappa.im;
     (*(out+4)).w += (*(in+1)).w*kappa.re + (*(in+1)).z*kappa.im;      
 

     
     (*(out+2)).x -= (*(in+2)).x*kappa.re - (*(in+2)).y*kappa.im;
     (*(out+2)).y -= (*(in+2)).y*kappa.re + (*(in+2)).x*kappa.im;
     (*(out+2)).x += (*(in+5)).x*kappa.re - (*(in+5)).y*kappa.im;
     (*(out+2)).y += (*(in+5)).y*kappa.re + (*(in+5)).x*kappa.im;
                
     (*(out+5)).x -= (*(in+5)).x*kappa.re - (*(in+5)).y*kappa.im;
     (*(out+5)).y -= (*(in+5)).y*kappa.re + (*(in+5)).x*kappa.im;   
     (*(out+5)).x += (*(in+2)).x*kappa.re - (*(in+2)).y*kappa.im;
     (*(out+5)).y += (*(in+2)).y*kappa.re + (*(in+2)).x*kappa.im; 
   
   
       
     (*(out+2)).z -= (*(in+2)).z*kappa.re - (*(in+2)).w*kappa.im;
     (*(out+2)).w -= (*(in+2)).w*kappa.re + (*(in+2)).z*kappa.im;
     (*(out+2)).z += (*(in+5)).z*kappa.re - (*(in+5)).w*kappa.im;
     (*(out+2)).w += (*(in+5)).w*kappa.re + (*(in+5)).z*kappa.im;
               
     (*(out+5)).z -= (*(in+5)).z*kappa.re - (*(in+5)).w*kappa.im;
     (*(out+5)).w -= (*(in+5)).w*kappa.re + (*(in+5)).z*kappa.im;   
     (*(out+5)).z += (*(in+2)).z*kappa.re - (*(in+2)).w*kappa.im;
     (*(out+5)).w += (*(in+2)).w*kappa.re + (*(in+2)).z*kappa.im; 
  
}







#ifdef RELATIVISTIC_BASIS
//  here comes P0+- for the relativistic basis
//  in this basis we have:
//
//  gamma0 =
//  -1  0  0  0 
//   0 -1  0  0
//   0  0  1  0
//   0  0  0  1
//
//

__device__ void dev_su3MtV_rel_up_d(dev_su3_d M, const dev_spinor_d * s, double4 * out){

(*(out+0)).x =  M[0][0].re*(*(s+0*DEVOFF)).x;
(*(out+0)).x -= M[0][0].im*(*(s+0*DEVOFF)).y;
(*(out+0)).x += M[0][1].re*(*(s+1*DEVOFF)).x;
(*(out+0)).x -= M[0][1].im*(*(s+1*DEVOFF)).y;
(*(out+0)).x += M[0][2].re*(*(s+2*DEVOFF)).x;
(*(out+0)).x -= M[0][2].im*(*(s+2*DEVOFF)).y;

(*(out+0)).y =  M[0][0].re*(*(s+0*DEVOFF)).y;
(*(out+0)).y += M[0][0].im*(*(s+0*DEVOFF)).x;
(*(out+0)).y += M[0][1].re*(*(s+1*DEVOFF)).y;
(*(out+0)).y += M[0][1].im*(*(s+1*DEVOFF)).x;
(*(out+0)).y += M[0][2].re*(*(s+2*DEVOFF)).y;
(*(out+0)).y += M[0][2].im*(*(s+2*DEVOFF)).x;


(*(out+0)).z =  M[1][0].re*(*(s+0*DEVOFF)).x;
(*(out+0)).z -= M[1][0].im*(*(s+0*DEVOFF)).y;
(*(out+0)).z += M[1][1].re*(*(s+1*DEVOFF)).x;
(*(out+0)).z -= M[1][1].im*(*(s+1*DEVOFF)).y;
(*(out+0)).z += M[1][2].re*(*(s+2*DEVOFF)).x;
(*(out+0)).z -= M[1][2].im*(*(s+2*DEVOFF)).y;

(*(out+0)).w =  M[1][0].re*(*(s+0*DEVOFF)).y;
(*(out+0)).w += M[1][0].im*(*(s+0*DEVOFF)).x;
(*(out+0)).w += M[1][1].re*(*(s+1*DEVOFF)).y;
(*(out+0)).w += M[1][1].im*(*(s+1*DEVOFF)).x;
(*(out+0)).w += M[1][2].re*(*(s+2*DEVOFF)).y;
(*(out+0)).w += M[1][2].im*(*(s+2*DEVOFF)).x;


(*(out+1)).x =  M[2][0].re*(*(s+0*DEVOFF)).x;
(*(out+1)).x -= M[2][0].im*(*(s+0*DEVOFF)).y;
(*(out+1)).x += M[2][1].re*(*(s+1*DEVOFF)).x;
(*(out+1)).x -= M[2][1].im*(*(s+1*DEVOFF)).y;
(*(out+1)).x += M[2][2].re*(*(s+2*DEVOFF)).x;
(*(out+1)).x -= M[2][2].im*(*(s+2*DEVOFF)).y;

(*(out+1)).y =  M[2][0].re*(*(s+0*DEVOFF)).y;
(*(out+1)).y += M[2][0].im*(*(s+0*DEVOFF)).x;
(*(out+1)).y += M[2][1].re*(*(s+1*DEVOFF)).y;
(*(out+1)).y += M[2][1].im*(*(s+1*DEVOFF)).x;
(*(out+1)).y += M[2][2].re*(*(s+2*DEVOFF)).y;
(*(out+1)).y += M[2][2].im*(*(s+2*DEVOFF)).x;


(*(out+1)).z =  M[0][0].re*(*(s+3*DEVOFF)).x;
(*(out+1)).z -= M[0][0].im*(*(s+3*DEVOFF)).y;
(*(out+1)).z += M[0][1].re*(*(s+4*DEVOFF)).x;
(*(out+1)).z -= M[0][1].im*(*(s+4*DEVOFF)).y;
(*(out+1)).z += M[0][2].re*(*(s+5*DEVOFF)).x;
(*(out+1)).z -= M[0][2].im*(*(s+5*DEVOFF)).y;

(*(out+1)).w =  M[0][0].re*(*(s+3*DEVOFF)).y;
(*(out+1)).w += M[0][0].im*(*(s+3*DEVOFF)).x;
(*(out+1)).w += M[0][1].re*(*(s+4*DEVOFF)).y;
(*(out+1)).w += M[0][1].im*(*(s+4*DEVOFF)).x;
(*(out+1)).w += M[0][2].re*(*(s+5*DEVOFF)).y;
(*(out+1)).w += M[0][2].im*(*(s+5*DEVOFF)).x;


(*(out+2)).x =  M[1][0].re*(*(s+3*DEVOFF)).x;
(*(out+2)).x -= M[1][0].im*(*(s+3*DEVOFF)).y;
(*(out+2)).x += M[1][1].re*(*(s+4*DEVOFF)).x;
(*(out+2)).x -= M[1][1].im*(*(s+4*DEVOFF)).y;
(*(out+2)).x += M[1][2].re*(*(s+5*DEVOFF)).x;
(*(out+2)).x -= M[1][2].im*(*(s+5*DEVOFF)).y;

(*(out+2)).y =  M[1][0].re*(*(s+3*DEVOFF)).y;
(*(out+2)).y += M[1][0].im*(*(s+3*DEVOFF)).x;
(*(out+2)).y += M[1][1].re*(*(s+4*DEVOFF)).y;
(*(out+2)).y += M[1][1].im*(*(s+4*DEVOFF)).x;
(*(out+2)).y += M[1][2].re*(*(s+5*DEVOFF)).y;
(*(out+2)).y += M[1][2].im*(*(s+5*DEVOFF)).x;


(*(out+2)).z =  M[2][0].re*(*(s+3*DEVOFF)).x;
(*(out+2)).z -= M[2][0].im*(*(s+3*DEVOFF)).y;
(*(out+2)).z += M[2][1].re*(*(s+4*DEVOFF)).x;
(*(out+2)).z -= M[2][1].im*(*(s+4*DEVOFF)).y;
(*(out+2)).z += M[2][2].re*(*(s+5*DEVOFF)).x;
(*(out+2)).z -= M[2][2].im*(*(s+5*DEVOFF)).y;

(*(out+2)).w =  M[2][0].re*(*(s+3*DEVOFF)).y;
(*(out+2)).w += M[2][0].im*(*(s+3*DEVOFF)).x;
(*(out+2)).w += M[2][1].re*(*(s+4*DEVOFF)).y;
(*(out+2)).w += M[2][1].im*(*(s+4*DEVOFF)).x;
(*(out+2)).w += M[2][2].re*(*(s+5*DEVOFF)).y;
(*(out+2)).w += M[2][2].im*(*(s+5*DEVOFF)).x;


(*(out+3)).x = 0.0;
(*(out+3)).y = 0.0;
(*(out+3)).z = 0.0;
(*(out+3)).w = 0.0;


(*(out+4)).x = 0.0;
(*(out+4)).y = 0.0;
(*(out+4)).z = 0.0;
(*(out+4)).w = 0.0;


(*(out+5)).x = 0.0;
(*(out+5)).y = 0.0;
(*(out+5)).z = 0.0;
(*(out+5)).w = 0.0;
}



__device__ void dev_su3MtV_rel_dn_d(dev_su3_d M, const dev_spinor_d * s, double4 * out){
  
(*(out+0)).x = 0.0;
(*(out+0)).y = 0.0;
(*(out+0)).z = 0.0;
(*(out+0)).w = 0.0; 


(*(out+1)).x = 0.0;
(*(out+1)).y = 0.0;
(*(out+1)).z = 0.0;
(*(out+1)).w = 0.0;


(*(out+2)).x = 0.0;
(*(out+2)).y = 0.0;
(*(out+2)).z = 0.0;
(*(out+2)).w = 0.0;

(*(out+3)).x =  M[0][0].re*(*(s+6*DEVOFF)).x;
(*(out+3)).x -= M[0][0].im*(*(s+6*DEVOFF)).y;
(*(out+3)).x += M[0][1].re*(*(s+7*DEVOFF)).x;
(*(out+3)).x -= M[0][1].im*(*(s+7*DEVOFF)).y;
(*(out+3)).x += M[0][2].re*(*(s+8*DEVOFF)).x;
(*(out+3)).x -= M[0][2].im*(*(s+8*DEVOFF)).y;

(*(out+3)).y =  M[0][0].re*(*(s+6*DEVOFF)).y;
(*(out+3)).y += M[0][0].im*(*(s+6*DEVOFF)).x;
(*(out+3)).y += M[0][1].re*(*(s+7*DEVOFF)).y;
(*(out+3)).y += M[0][1].im*(*(s+7*DEVOFF)).x;
(*(out+3)).y += M[0][2].re*(*(s+8*DEVOFF)).y;
(*(out+3)).y += M[0][2].im*(*(s+8*DEVOFF)).x;


(*(out+3)).z =  M[1][0].re*(*(s+6*DEVOFF)).x;
(*(out+3)).z -= M[1][0].im*(*(s+6*DEVOFF)).y;
(*(out+3)).z += M[1][1].re*(*(s+7*DEVOFF)).x;
(*(out+3)).z -= M[1][1].im*(*(s+7*DEVOFF)).y;
(*(out+3)).z += M[1][2].re*(*(s+8*DEVOFF)).x;
(*(out+3)).z -= M[1][2].im*(*(s+8*DEVOFF)).y;

(*(out+3)).w =  M[1][0].re*(*(s+6*DEVOFF)).y;
(*(out+3)).w += M[1][0].im*(*(s+6*DEVOFF)).x;
(*(out+3)).w += M[1][1].re*(*(s+7*DEVOFF)).y;
(*(out+3)).w += M[1][1].im*(*(s+7*DEVOFF)).x;
(*(out+3)).w += M[1][2].re*(*(s+8*DEVOFF)).y;
(*(out+3)).w += M[1][2].im*(*(s+8*DEVOFF)).x;


(*(out+4)).x =  M[2][0].re*(*(s+6*DEVOFF)).x;
(*(out+4)).x -= M[2][0].im*(*(s+6*DEVOFF)).y;
(*(out+4)).x += M[2][1].re*(*(s+7*DEVOFF)).x;
(*(out+4)).x -= M[2][1].im*(*(s+7*DEVOFF)).y;
(*(out+4)).x += M[2][2].re*(*(s+8*DEVOFF)).x;
(*(out+4)).x -= M[2][2].im*(*(s+8*DEVOFF)).y;

(*(out+4)).y =  M[2][0].re*(*(s+6*DEVOFF)).y;
(*(out+4)).y += M[2][0].im*(*(s+6*DEVOFF)).x;
(*(out+4)).y += M[2][1].re*(*(s+7*DEVOFF)).y;
(*(out+4)).y += M[2][1].im*(*(s+7*DEVOFF)).x;
(*(out+4)).y += M[2][2].re*(*(s+8*DEVOFF)).y;
(*(out+4)).y += M[2][2].im*(*(s+8*DEVOFF)).x;


(*(out+4)).z =  M[0][0].re*(*(s+9*DEVOFF)).x;
(*(out+4)).z -= M[0][0].im*(*(s+9*DEVOFF)).y;
(*(out+4)).z += M[0][1].re*(*(s+10*DEVOFF)).x;
(*(out+4)).z -= M[0][1].im*(*(s+10*DEVOFF)).y;
(*(out+4)).z += M[0][2].re*(*(s+11*DEVOFF)).x;
(*(out+4)).z -= M[0][2].im*(*(s+11*DEVOFF)).y;

(*(out+4)).w =  M[0][0].re*(*(s+9*DEVOFF)).y;
(*(out+4)).w += M[0][0].im*(*(s+9*DEVOFF)).x;
(*(out+4)).w += M[0][1].re*(*(s+10*DEVOFF)).y;
(*(out+4)).w += M[0][1].im*(*(s+10*DEVOFF)).x;
(*(out+4)).w += M[0][2].re*(*(s+11*DEVOFF)).y;
(*(out+4)).w += M[0][2].im*(*(s+11*DEVOFF)).x;


(*(out+5)).x =  M[1][0].re*(*(s+9*DEVOFF)).x;
(*(out+5)).x -= M[1][0].im*(*(s+9*DEVOFF)).y;
(*(out+5)).x += M[1][1].re*(*(s+10*DEVOFF)).x;
(*(out+5)).x -= M[1][1].im*(*(s+10*DEVOFF)).y;
(*(out+5)).x += M[1][2].re*(*(s+11*DEVOFF)).x;
(*(out+5)).x -= M[1][2].im*(*(s+11*DEVOFF)).y;

(*(out+5)).y =  M[1][0].re*(*(s+9*DEVOFF)).y;
(*(out+5)).y += M[1][0].im*(*(s+9*DEVOFF)).x;
(*(out+5)).y += M[1][1].re*(*(s+10*DEVOFF)).y;
(*(out+5)).y += M[1][1].im*(*(s+10*DEVOFF)).x;
(*(out+5)).y += M[1][2].re*(*(s+11*DEVOFF)).y;
(*(out+5)).y += M[1][2].im*(*(s+11*DEVOFF)).x;


(*(out+5)).z =  M[2][0].re*(*(s+9*DEVOFF)).x;
(*(out+5)).z -= M[2][0].im*(*(s+9*DEVOFF)).y;
(*(out+5)).z += M[2][1].re*(*(s+10*DEVOFF)).x;
(*(out+5)).z -= M[2][1].im*(*(s+10*DEVOFF)).y;
(*(out+5)).z += M[2][2].re*(*(s+11*DEVOFF)).x;
(*(out+5)).z -= M[2][2].im*(*(s+11*DEVOFF)).y;

(*(out+5)).w =  M[2][0].re*(*(s+9*DEVOFF)).y;
(*(out+5)).w += M[2][0].im*(*(s+9*DEVOFF)).x;
(*(out+5)).w += M[2][1].re*(*(s+10*DEVOFF)).y;
(*(out+5)).w += M[2][1].im*(*(s+10*DEVOFF)).x;
(*(out+5)).w += M[2][2].re*(*(s+11*DEVOFF)).y;
(*(out+5)).w += M[2][2].im*(*(s+11*DEVOFF)).x;
}



//-kappa(r - gamma_mu)
__device__ void dev_kappaP0_plus_relativistic_d(double4 * out, double4 * in, dev_complex_d kappa){


     (*(out+0)).x -= 2.0*( (*(in+0)).x*kappa.re - (*(in+0)).y*kappa.im );
     (*(out+0)).y -= 2.0*( (*(in+0)).y*kappa.re + (*(in+0)).x*kappa.im );
     (*(out+0)).z -= 2.0*( (*(in+0)).z*kappa.re - (*(in+0)).w*kappa.im );
     (*(out+0)).w -= 2.0*( (*(in+0)).w*kappa.re + (*(in+0)).z*kappa.im );
     
     (*(out+1)).x -= 2.0*( (*(in+1)).x*kappa.re - (*(in+1)).y*kappa.im );
     (*(out+1)).y -= 2.0*( (*(in+1)).y*kappa.re + (*(in+1)).x*kappa.im );     
     (*(out+1)).z -= 2.0*( (*(in+1)).z*kappa.re - (*(in+1)).w*kappa.im );
     (*(out+1)).w -= 2.0*( (*(in+1)).w*kappa.re + (*(in+1)).z*kappa.im );     
     
     (*(out+2)).x -= 2.0*( (*(in+2)).x*kappa.re - (*(in+2)).y*kappa.im );
     (*(out+2)).y -= 2.0*( (*(in+2)).y*kappa.re + (*(in+2)).x*kappa.im );
     (*(out+2)).z -= 2.0*( (*(in+2)).z*kappa.re - (*(in+2)).w*kappa.im );
     (*(out+2)).w -= 2.0*( (*(in+2)).w*kappa.re + (*(in+2)).z*kappa.im );

}






//-kappa(r + gamma_mu)
__device__ void dev_kappaP0_minus_relativistic_d(double4 * out, double4 * in, dev_complex_d kappa){

    
     (*(out+3)).x -= 2.0*( (*(in+3)).x*kappa.re - (*(in+3)).y*kappa.im );
     (*(out+3)).y -= 2.0*( (*(in+3)).y*kappa.re + (*(in+3)).x*kappa.im );     
     (*(out+3)).z -= 2.0*( (*(in+3)).z*kappa.re - (*(in+3)).w*kappa.im );
     (*(out+3)).w -= 2.0*( (*(in+3)).w*kappa.re + (*(in+3)).z*kappa.im );


     (*(out+4)).x -= 2.0*( (*(in+4)).x*kappa.re - (*(in+4)).y*kappa.im );
     (*(out+4)).y -= 2.0*( (*(in+4)).y*kappa.re + (*(in+4)).x*kappa.im );     
     (*(out+4)).z -= 2.0*( (*(in+4)).z*kappa.re - (*(in+4)).w*kappa.im );
     (*(out+4)).w -= 2.0*( (*(in+4)).w*kappa.re + (*(in+4)).z*kappa.im );    

           
     (*(out+5)).x -= 2.0*( (*(in+5)).x*kappa.re - (*(in+5)).y*kappa.im );
     (*(out+5)).y -= 2.0*( (*(in+5)).y*kappa.re + (*(in+5)).x*kappa.im );            
     (*(out+5)).z -= 2.0*( (*(in+5)).z*kappa.re - (*(in+5)).w*kappa.im );
     (*(out+5)).w -= 2.0*( (*(in+5)).w*kappa.re + (*(in+5)).z*kappa.im );   

  
}

//RELATIVISTIC_BASIS
#endif 








// multiplies su3-Matrix, kappa, spin-projector onto spinor
// this is +z 
// uses spin projection reduction 
// with 1-gamma_3
//
// | 1  0  +i  0 |     s0      s0 + i s2
// | 0  1   0 -i |     s1      s1 - i s3
// |-i  0   1  0 |  X  s2   =  -i(s0 + i s2)
// | 0 +i   0  1 |     s3       i(s1 - i s3)
//
// second half spinor is proportional to first one!!!
//
//-kappa(r - gamma_mu) kappa reell !!!!
__device__ void dev_su3MtV_kappaP3_plus_d(dev_su3_d M, const dev_spinor_d * s, double4 * out, double kappa){

  
  double4 sh0,sh1,sh2,zh0,zh1,zh2;
  
  
  //first apply Projector on upper halfspinor 
     sh0.x = kappa*( (*(s+0*DEVOFF)).x - (*(s+6*DEVOFF)).y);
     sh0.y = kappa*( (*(s+0*DEVOFF)).y + (*(s+6*DEVOFF)).x);
   
     sh0.z = kappa*( (*(s+1*DEVOFF)).x - (*(s+7*DEVOFF)).y);
     sh0.w = kappa*( (*(s+1*DEVOFF)).y + (*(s+7*DEVOFF)).x);
   
     sh1.x = kappa*( (*(s+2*DEVOFF)).x - (*(s+8*DEVOFF)).y);
     sh1.y = kappa*( (*(s+2*DEVOFF)).y + (*(s+8*DEVOFF)).x); 
     
  

//Multiply by gauge field  
zh0.x =   M[0][0].re*sh0.x;
zh0.x -=  M[0][0].im*sh0.y; 
zh0.x +=  M[0][1].re*sh0.z;
zh0.x -=  M[0][1].im*sh0.w;
zh0.x +=  M[0][2].re*sh1.x;
zh0.x -=  M[0][2].im*sh1.y;
zh0.y =   M[0][0].re*sh0.y;
zh0.y +=  M[0][0].im*sh0.x;
zh0.y +=  M[0][1].re*sh0.w;
zh0.y +=  M[0][1].im*sh0.z;
zh0.y +=  M[0][2].re*sh1.y;
zh0.y +=  M[0][2].im*sh1.x;
(*(out+0)).x -= zh0.x;
(*(out+0)).y -= zh0.y;

zh0.z =  M[1][0].re*sh0.x;
zh0.z -= M[1][0].im*sh0.y;
zh0.z += M[1][1].re*sh0.z;
zh0.z -= M[1][1].im*sh0.w;
zh0.z += M[1][2].re*sh1.x;
zh0.z -= M[1][2].im*sh1.y;
zh0.w =  M[1][0].re*sh0.y;
zh0.w += M[1][0].im*sh0.x;
zh0.w += M[1][1].re*sh0.w;
zh0.w += M[1][1].im*sh0.z;
zh0.w += M[1][2].re*sh1.y;
zh0.w += M[1][2].im*sh1.x;
(*(out+0)).z -= zh0.z;
(*(out+0)).w -= zh0.w;

zh1.x =  M[2][0].re*sh0.x;
zh1.x -= M[2][0].im*sh0.y;
zh1.x += M[2][1].re*sh0.z;
zh1.x -= M[2][1].im*sh0.w;
zh1.x += M[2][2].re*sh1.x;
zh1.x -= M[2][2].im*sh1.y;
zh1.y =  M[2][0].re*sh0.y;
zh1.y += M[2][0].im*sh0.x;
zh1.y += M[2][1].re*sh0.w;
zh1.y += M[2][1].im*sh0.z;
zh1.y += M[2][2].re*sh1.y;
zh1.y += M[2][2].im*sh1.x;
(*(out+1)).x -= zh1.x;
(*(out+1)).y -= zh1.y;


     sh1.z = kappa*( (*(s+3*DEVOFF)).x + (*(s+9*DEVOFF)).y);
     sh1.w = kappa*( (*(s+3*DEVOFF)).y - (*(s+9*DEVOFF)).x);
         
     sh2.x = kappa*( (*(s+4*DEVOFF)).x + (*(s+10*DEVOFF)).y);
     sh2.y = kappa*( (*(s+4*DEVOFF)).y - (*(s+10*DEVOFF)).x);
       
     sh2.z = kappa*( (*(s+5*DEVOFF)).x + (*(s+11*DEVOFF)).y);
     sh2.w = kappa*( (*(s+5*DEVOFF)).y - (*(s+11*DEVOFF)).x);




zh1.z =  M[0][0].re*sh1.z;
zh1.z -= M[0][0].im*sh1.w;
zh1.z += M[0][1].re*sh2.x;
zh1.z -= M[0][1].im*sh2.y;
zh1.z += M[0][2].re*sh2.z;
zh1.z -= M[0][2].im*sh2.w;
zh1.w =  M[0][0].re*sh1.w;
zh1.w += M[0][0].im*sh1.z;
zh1.w += M[0][1].re*sh2.y;
zh1.w += M[0][1].im*sh2.x;
zh1.w += M[0][2].re*sh2.w;
zh1.w += M[0][2].im*sh2.z;
(*(out+1)).z -= zh1.z;
(*(out+1)).w -= zh1.w;

zh2.x =  M[1][0].re*sh1.z;
zh2.x -= M[1][0].im*sh1.w;
zh2.x += M[1][1].re*sh2.x;
zh2.x -= M[1][1].im*sh2.y;
zh2.x += M[1][2].re*sh2.z;
zh2.x -= M[1][2].im*sh2.w;
zh2.y =  M[1][0].re*sh1.w;
zh2.y += M[1][0].im*sh1.z;
zh2.y += M[1][1].re*sh2.y;
zh2.y += M[1][1].im*sh2.x;
zh2.y += M[1][2].re*sh2.w;
zh2.y += M[1][2].im*sh2.z;
(*(out+2)).x -= zh2.x;
(*(out+2)).y -= zh2.y;

zh2.z =  M[2][0].re*sh1.z;
zh2.z -= M[2][0].im*sh1.w;
zh2.z += M[2][1].re*sh2.x;
zh2.z -= M[2][1].im*sh2.y;
zh2.z += M[2][2].re*sh2.z;
zh2.z -= M[2][2].im*sh2.w;
zh2.w =  M[2][0].re*sh1.w;
zh2.w += M[2][0].im*sh1.z;
zh2.w += M[2][1].re*sh2.y;
zh2.w += M[2][1].im*sh2.x;
zh2.w += M[2][2].re*sh2.w;
zh2.w += M[2][2].im*sh2.z;
(*(out+2)).z -= zh2.z;
(*(out+2)).w -= zh2.w;


//Reconstruct lower half spinor
//this is just a multiplication by "-i"
     (*(out+3)).x -= zh0.y;
     (*(out+3)).y += zh0.x; 

     (*(out+3)).z -= zh0.w;
     (*(out+3)).w += zh0.z; 
     
     (*(out+4)).x -= zh1.y;
     (*(out+4)).y += zh1.x;       


//this is just a multiplication by "i"     
     (*(out+4)).z += zh1.w;
     (*(out+4)).w -= zh1.z;       

     (*(out+5)).x += zh2.y;
     (*(out+5)).y -= zh2.x;  
     
     (*(out+5)).z += zh2.w;
     (*(out+5)).w -= zh2.z;
     
}


__device__ void dev_su3MtV_kappaP3_plus_spintex_d(dev_su3_d M, int pos, double4 * out, double kappa, int updn){

  
  double4 sh0,sh1,sh2,zh0,zh1,zh2;
  double2 s0,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11;
  
  if(updn==1){
    s0 = fetch1D_spin_d(spin_d_tex,pos+0*DEVOFF);
    s1 = fetch1D_spin_d(spin_d_tex,pos+1*DEVOFF);
    s2 = fetch1D_spin_d(spin_d_tex,pos+2*DEVOFF);
    s6 = fetch1D_spin_d(spin_d_tex,pos+6*DEVOFF);
    s7 = fetch1D_spin_d(spin_d_tex,pos+7*DEVOFF);
    s8 = fetch1D_spin_d(spin_d_tex,pos+8*DEVOFF);  
  }
  else{
    s0 = fetch1D_spin_d(spin_d_tex_dn,pos+0*DEVOFF);
    s1 = fetch1D_spin_d(spin_d_tex_dn,pos+1*DEVOFF);
    s2 = fetch1D_spin_d(spin_d_tex_dn,pos+2*DEVOFF);
    s6 = fetch1D_spin_d(spin_d_tex_dn,pos+6*DEVOFF);
    s7 = fetch1D_spin_d(spin_d_tex_dn,pos+7*DEVOFF); 
    s8 = fetch1D_spin_d(spin_d_tex_dn,pos+8*DEVOFF);     
  }
  //first apply Projector on upper halfspinor 
     sh0.x = kappa*( s0.x - s6.y);
     sh0.y = kappa*( s0.y + s6.x);
   
     sh0.z = kappa*( s1.x - s7.y);
     sh0.w = kappa*( s1.y + s7.x);
   
     sh1.x = kappa*( s2.x - s8.y);
     sh1.y = kappa*( s2.y + s8.x); 
     
  

//Multiply by gauge field  
zh0.x =   M[0][0].re*sh0.x;
zh0.x -=  M[0][0].im*sh0.y; 
zh0.x +=  M[0][1].re*sh0.z;
zh0.x -=  M[0][1].im*sh0.w;
zh0.x +=  M[0][2].re*sh1.x;
zh0.x -=  M[0][2].im*sh1.y;
zh0.y =   M[0][0].re*sh0.y;
zh0.y +=  M[0][0].im*sh0.x;
zh0.y +=  M[0][1].re*sh0.w;
zh0.y +=  M[0][1].im*sh0.z;
zh0.y +=  M[0][2].re*sh1.y;
zh0.y +=  M[0][2].im*sh1.x;
(*(out+0)).x -= zh0.x;
(*(out+0)).y -= zh0.y;

zh0.z =  M[1][0].re*sh0.x;
zh0.z -= M[1][0].im*sh0.y;
zh0.z += M[1][1].re*sh0.z;
zh0.z -= M[1][1].im*sh0.w;
zh0.z += M[1][2].re*sh1.x;
zh0.z -= M[1][2].im*sh1.y;
zh0.w =  M[1][0].re*sh0.y;
zh0.w += M[1][0].im*sh0.x;
zh0.w += M[1][1].re*sh0.w;
zh0.w += M[1][1].im*sh0.z;
zh0.w += M[1][2].re*sh1.y;
zh0.w += M[1][2].im*sh1.x;
(*(out+0)).z -= zh0.z;
(*(out+0)).w -= zh0.w;

zh1.x =  M[2][0].re*sh0.x;
zh1.x -= M[2][0].im*sh0.y;
zh1.x += M[2][1].re*sh0.z;
zh1.x -= M[2][1].im*sh0.w;
zh1.x += M[2][2].re*sh1.x;
zh1.x -= M[2][2].im*sh1.y;
zh1.y =  M[2][0].re*sh0.y;
zh1.y += M[2][0].im*sh0.x;
zh1.y += M[2][1].re*sh0.w;
zh1.y += M[2][1].im*sh0.z;
zh1.y += M[2][2].re*sh1.y;
zh1.y += M[2][2].im*sh1.x;
(*(out+1)).x -= zh1.x;
(*(out+1)).y -= zh1.y;

  if(updn==1){
    s3 = fetch1D_spin_d(spin_d_tex,pos+3*DEVOFF);
    s4 = fetch1D_spin_d(spin_d_tex,pos+4*DEVOFF);
    s5 = fetch1D_spin_d(spin_d_tex,pos+5*DEVOFF);
    s9 = fetch1D_spin_d(spin_d_tex,pos+9*DEVOFF);
    s10 = fetch1D_spin_d(spin_d_tex,pos+10*DEVOFF);
    s11 = fetch1D_spin_d(spin_d_tex,pos+11*DEVOFF);  
  }
  else{
    s3 = fetch1D_spin_d(spin_d_tex_dn,pos+3*DEVOFF);
    s4 = fetch1D_spin_d(spin_d_tex_dn,pos+4*DEVOFF);
    s5 = fetch1D_spin_d(spin_d_tex_dn,pos+5*DEVOFF);
    s9 = fetch1D_spin_d(spin_d_tex_dn,pos+9*DEVOFF);
    s10 = fetch1D_spin_d(spin_d_tex_dn,pos+10*DEVOFF); 
    s11 = fetch1D_spin_d(spin_d_tex_dn,pos+11*DEVOFF);     
  }
  
     sh1.z = kappa*( s3.x + s9.y);
     sh1.w = kappa*( s3.y - s9.x);
         
     sh2.x = kappa*( s4.x + s10.y);
     sh2.y = kappa*( s4.y - s10.x);
       
     sh2.z = kappa*( s5.x + s11.y);
     sh2.w = kappa*( s5.y - s11.x);




zh1.z =  M[0][0].re*sh1.z;
zh1.z -= M[0][0].im*sh1.w;
zh1.z += M[0][1].re*sh2.x;
zh1.z -= M[0][1].im*sh2.y;
zh1.z += M[0][2].re*sh2.z;
zh1.z -= M[0][2].im*sh2.w;
zh1.w =  M[0][0].re*sh1.w;
zh1.w += M[0][0].im*sh1.z;
zh1.w += M[0][1].re*sh2.y;
zh1.w += M[0][1].im*sh2.x;
zh1.w += M[0][2].re*sh2.w;
zh1.w += M[0][2].im*sh2.z;
(*(out+1)).z -= zh1.z;
(*(out+1)).w -= zh1.w;

zh2.x =  M[1][0].re*sh1.z;
zh2.x -= M[1][0].im*sh1.w;
zh2.x += M[1][1].re*sh2.x;
zh2.x -= M[1][1].im*sh2.y;
zh2.x += M[1][2].re*sh2.z;
zh2.x -= M[1][2].im*sh2.w;
zh2.y =  M[1][0].re*sh1.w;
zh2.y += M[1][0].im*sh1.z;
zh2.y += M[1][1].re*sh2.y;
zh2.y += M[1][1].im*sh2.x;
zh2.y += M[1][2].re*sh2.w;
zh2.y += M[1][2].im*sh2.z;
(*(out+2)).x -= zh2.x;
(*(out+2)).y -= zh2.y;

zh2.z =  M[2][0].re*sh1.z;
zh2.z -= M[2][0].im*sh1.w;
zh2.z += M[2][1].re*sh2.x;
zh2.z -= M[2][1].im*sh2.y;
zh2.z += M[2][2].re*sh2.z;
zh2.z -= M[2][2].im*sh2.w;
zh2.w =  M[2][0].re*sh1.w;
zh2.w += M[2][0].im*sh1.z;
zh2.w += M[2][1].re*sh2.y;
zh2.w += M[2][1].im*sh2.x;
zh2.w += M[2][2].re*sh2.w;
zh2.w += M[2][2].im*sh2.z;
(*(out+2)).z -= zh2.z;
(*(out+2)).w -= zh2.w;


//Reconstruct lower half spinor
//this is just a multiplication by "-i"
     (*(out+3)).x -= zh0.y;
     (*(out+3)).y += zh0.x; 

     (*(out+3)).z -= zh0.w;
     (*(out+3)).w += zh0.z; 
     
     (*(out+4)).x -= zh1.y;
     (*(out+4)).y += zh1.x;       


//this is just a multiplication by "i"     
     (*(out+4)).z += zh1.w;
     (*(out+4)).w -= zh1.z;       

     (*(out+5)).x += zh2.y;
     (*(out+5)).y -= zh2.x;  
     
     (*(out+5)).z += zh2.w;
     (*(out+5)).w -= zh2.z;
     
}



// multiplies su3-Matrix, kappa, spin-projector onto spinor
// this is -z 
// uses spin projection reduction 
// with 1+gamma_3
//
// | 1  0  -i  0 |     s0      s0 - i s2
// | 0  1   0 +i |     s1      s1 + i s3
// |+i  0   1  0 |  X  s2   =   i(s0 - i s2)
// | 0 -i   0  1 |     s3      -i(s1 + i s3)
//
// second half spinor is proportional to first one!!!
//
//-kappa(r + gamma_mu) kappa reell !!!
__device__ void dev_su3MtV_kappaP3_minus_d(dev_su3_d M, const dev_spinor_d * s, double4 * out, double kappa){

  double4 sh0,sh1,sh2,zh0,zh1,zh2;
  
  
  //first apply Projector on upper halfspinor 
     sh0.x = kappa*( (*(s+0*DEVOFF)).x + (*(s+6*DEVOFF)).y);
     sh0.y = kappa*( (*(s+0*DEVOFF)).y - (*(s+6*DEVOFF)).x);
   
     sh0.z = kappa*( (*(s+1*DEVOFF)).x + (*(s+7*DEVOFF)).y);
     sh0.w = kappa*( (*(s+1*DEVOFF)).y - (*(s+7*DEVOFF)).x);
   
     sh1.x = kappa*( (*(s+2*DEVOFF)).x + (*(s+8*DEVOFF)).y);
     sh1.y = kappa*( (*(s+2*DEVOFF)).y - (*(s+8*DEVOFF)).x);
        


  

//Multiply by gauge field  
zh0.x  = M[0][0].re*sh0.x;
zh0.x -= M[0][0].im*sh0.y;
zh0.x += M[0][1].re*sh0.z;
zh0.x -= M[0][1].im*sh0.w;
zh0.x += M[0][2].re*sh1.x;
zh0.x -= M[0][2].im*sh1.y;
zh0.y  = M[0][0].re*sh0.y;
zh0.y += M[0][0].im*sh0.x;
zh0.y += M[0][1].re*sh0.w;
zh0.y += M[0][1].im*sh0.z;
zh0.y += M[0][2].re*sh1.y;
zh0.y += M[0][2].im*sh1.x;
(*(out+0)).x -= zh0.x;
(*(out+0)).y -= zh0.y;

zh0.z  = M[1][0].re*sh0.x;
zh0.z -= M[1][0].im*sh0.y;
zh0.z += M[1][1].re*sh0.z;
zh0.z -= M[1][1].im*sh0.w;
zh0.z += M[1][2].re*sh1.x;
zh0.z -= M[1][2].im*sh1.y;
zh0.w  = M[1][0].re*sh0.y;
zh0.w += M[1][0].im*sh0.x;
zh0.w += M[1][1].re*sh0.w;
zh0.w += M[1][1].im*sh0.z;
zh0.w += M[1][2].re*sh1.y;
zh0.w += M[1][2].im*sh1.x;
(*(out+0)).z -= zh0.z;
(*(out+0)).w -= zh0.w;

zh1.x  = M[2][0].re*sh0.x;
zh1.x -= M[2][0].im*sh0.y;
zh1.x += M[2][1].re*sh0.z;
zh1.x -= M[2][1].im*sh0.w;
zh1.x += M[2][2].re*sh1.x;
zh1.x -= M[2][2].im*sh1.y;
zh1.y  = M[2][0].re*sh0.y;
zh1.y += M[2][0].im*sh0.x;
zh1.y += M[2][1].re*sh0.w;
zh1.y += M[2][1].im*sh0.z;
zh1.y += M[2][2].re*sh1.y;
zh1.y += M[2][2].im*sh1.x;
(*(out+1)).x -= zh1.x;
(*(out+1)).y -= zh1.y;


     sh1.z = kappa*( (*(s+3*DEVOFF)).x - (*(s+9*DEVOFF)).y);
     sh1.w = kappa*( (*(s+3*DEVOFF)).y + (*(s+9*DEVOFF)).x);
         
     sh2.x = kappa*( (*(s+4*DEVOFF)).x - (*(s+10*DEVOFF)).y);
     sh2.y = kappa*( (*(s+4*DEVOFF)).y + (*(s+10*DEVOFF)).x);
       
     sh2.z = kappa*( (*(s+5*DEVOFF)).x - (*(s+11*DEVOFF)).y);
     sh2.w = kappa*( (*(s+5*DEVOFF)).y + (*(s+11*DEVOFF)).x);

     
zh1.z  = M[0][0].re*sh1.z;
zh1.z -= M[0][0].im*sh1.w;
zh1.z += M[0][1].re*sh2.x;
zh1.z -= M[0][1].im*sh2.y;
zh1.z += M[0][2].re*sh2.z;
zh1.z -= M[0][2].im*sh2.w;
zh1.w  = M[0][0].re*sh1.w;
zh1.w += M[0][0].im*sh1.z;
zh1.w += M[0][1].re*sh2.y;
zh1.w += M[0][1].im*sh2.x;
zh1.w += M[0][2].re*sh2.w;
zh1.w += M[0][2].im*sh2.z;
(*(out+1)).z -= zh1.z;
(*(out+1)).w -= zh1.w;

zh2.x  = M[1][0].re*sh1.z;
zh2.x -= M[1][0].im*sh1.w;
zh2.x += M[1][1].re*sh2.x;
zh2.x -= M[1][1].im*sh2.y;
zh2.x += M[1][2].re*sh2.z;
zh2.x -= M[1][2].im*sh2.w;
zh2.y  = M[1][0].re*sh1.w;
zh2.y += M[1][0].im*sh1.z;
zh2.y += M[1][1].re*sh2.y;
zh2.y += M[1][1].im*sh2.x;
zh2.y += M[1][2].re*sh2.w;
zh2.y += M[1][2].im*sh2.z;
(*(out+2)).x -= zh2.x;
(*(out+2)).y -= zh2.y;

zh2.z  = M[2][0].re*sh1.z;
zh2.z -= M[2][0].im*sh1.w;
zh2.z += M[2][1].re*sh2.x;
zh2.z -= M[2][1].im*sh2.y;
zh2.z += M[2][2].re*sh2.z;
zh2.z -= M[2][2].im*sh2.w;
zh2.w  = M[2][0].re*sh1.w;
zh2.w += M[2][0].im*sh1.z;
zh2.w += M[2][1].re*sh2.y;
zh2.w += M[2][1].im*sh2.x;
zh2.w += M[2][2].re*sh2.w;
zh2.w += M[2][2].im*sh2.z;
(*(out+2)).z -= zh2.z;
(*(out+2)).w -= zh2.w;


//Reconstruct lower half spinor
//this is just a multiplication by "i"
     (*(out+3)).x += zh0.y;
     (*(out+3)).y -= zh0.x; 

     (*(out+3)).z += zh0.w;
     (*(out+3)).w -= zh0.z; 
     
     (*(out+4)).x += zh1.y;
     (*(out+4)).y -= zh1.x;       


//this is just a multiplication by "-i"     
     (*(out+4)).z -= zh1.w;
     (*(out+4)).w += zh1.z;       

     (*(out+5)).x -= zh2.y;
     (*(out+5)).y += zh2.x;  
     
     (*(out+5)).z -= zh2.w;
     (*(out+5)).w += zh2.z;

}



__device__ void dev_su3MtV_kappaP3_minus_spintex_d(dev_su3_d M, int pos, double4 * out, double kappa, int updn){

  double4 sh0,sh1,sh2,zh0,zh1,zh2;
  double2 s0,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11;
  
  if(updn==1){
    s0 = fetch1D_spin_d(spin_d_tex,pos+0*DEVOFF);
    s1 = fetch1D_spin_d(spin_d_tex,pos+1*DEVOFF);
    s2 = fetch1D_spin_d(spin_d_tex,pos+2*DEVOFF);
    s6 = fetch1D_spin_d(spin_d_tex,pos+6*DEVOFF);
    s7 = fetch1D_spin_d(spin_d_tex,pos+7*DEVOFF);
    s8 = fetch1D_spin_d(spin_d_tex,pos+8*DEVOFF);  
  }
  else{
    s0 = fetch1D_spin_d(spin_d_tex_dn,pos+0*DEVOFF);
    s1 = fetch1D_spin_d(spin_d_tex_dn,pos+1*DEVOFF);
    s2 = fetch1D_spin_d(spin_d_tex_dn,pos+2*DEVOFF);
    s6 = fetch1D_spin_d(spin_d_tex_dn,pos+6*DEVOFF);
    s7 = fetch1D_spin_d(spin_d_tex_dn,pos+7*DEVOFF); 
    s8 = fetch1D_spin_d(spin_d_tex_dn,pos+8*DEVOFF);     
  }  
  
  //first apply Projector on upper halfspinor 
     sh0.x = kappa*( s0.x + s6.y);
     sh0.y = kappa*( s0.y - s6.x);
   
     sh0.z = kappa*( s1.x + s7.y);
     sh0.w = kappa*( s1.y - s7.x);
   
     sh1.x = kappa*( s2.x + s8.y);
     sh1.y = kappa*( s2.y - s8.x);
        


  

//Multiply by gauge field  
zh0.x  = M[0][0].re*sh0.x;
zh0.x -= M[0][0].im*sh0.y;
zh0.x += M[0][1].re*sh0.z;
zh0.x -= M[0][1].im*sh0.w;
zh0.x += M[0][2].re*sh1.x;
zh0.x -= M[0][2].im*sh1.y;
zh0.y  = M[0][0].re*sh0.y;
zh0.y += M[0][0].im*sh0.x;
zh0.y += M[0][1].re*sh0.w;
zh0.y += M[0][1].im*sh0.z;
zh0.y += M[0][2].re*sh1.y;
zh0.y += M[0][2].im*sh1.x;
(*(out+0)).x -= zh0.x;
(*(out+0)).y -= zh0.y;

zh0.z  = M[1][0].re*sh0.x;
zh0.z -= M[1][0].im*sh0.y;
zh0.z += M[1][1].re*sh0.z;
zh0.z -= M[1][1].im*sh0.w;
zh0.z += M[1][2].re*sh1.x;
zh0.z -= M[1][2].im*sh1.y;
zh0.w  = M[1][0].re*sh0.y;
zh0.w += M[1][0].im*sh0.x;
zh0.w += M[1][1].re*sh0.w;
zh0.w += M[1][1].im*sh0.z;
zh0.w += M[1][2].re*sh1.y;
zh0.w += M[1][2].im*sh1.x;
(*(out+0)).z -= zh0.z;
(*(out+0)).w -= zh0.w;

zh1.x  = M[2][0].re*sh0.x;
zh1.x -= M[2][0].im*sh0.y;
zh1.x += M[2][1].re*sh0.z;
zh1.x -= M[2][1].im*sh0.w;
zh1.x += M[2][2].re*sh1.x;
zh1.x -= M[2][2].im*sh1.y;
zh1.y  = M[2][0].re*sh0.y;
zh1.y += M[2][0].im*sh0.x;
zh1.y += M[2][1].re*sh0.w;
zh1.y += M[2][1].im*sh0.z;
zh1.y += M[2][2].re*sh1.y;
zh1.y += M[2][2].im*sh1.x;
(*(out+1)).x -= zh1.x;
(*(out+1)).y -= zh1.y;

  if(updn==1){
    s3 = fetch1D_spin_d(spin_d_tex,pos+3*DEVOFF);
    s4 = fetch1D_spin_d(spin_d_tex,pos+4*DEVOFF);
    s5 = fetch1D_spin_d(spin_d_tex,pos+5*DEVOFF);
    s9 = fetch1D_spin_d(spin_d_tex,pos+9*DEVOFF);
    s10 = fetch1D_spin_d(spin_d_tex,pos+10*DEVOFF);
    s11 = fetch1D_spin_d(spin_d_tex,pos+11*DEVOFF);  
  }
  else{
    s3 = fetch1D_spin_d(spin_d_tex_dn,pos+3*DEVOFF);
    s4 = fetch1D_spin_d(spin_d_tex_dn,pos+4*DEVOFF);
    s5 = fetch1D_spin_d(spin_d_tex_dn,pos+5*DEVOFF);
    s9 = fetch1D_spin_d(spin_d_tex_dn,pos+9*DEVOFF);
    s10 = fetch1D_spin_d(spin_d_tex_dn,pos+10*DEVOFF); 
    s11 = fetch1D_spin_d(spin_d_tex_dn,pos+11*DEVOFF);     
  }
  
     sh1.z = kappa*( s3.x - s9.y);
     sh1.w = kappa*( s3.y + s9.x);
         
     sh2.x = kappa*( s4.x - s10.y);
     sh2.y = kappa*( s4.y + s10.x);
       
     sh2.z = kappa*( s5.x - s11.y);
     sh2.w = kappa*( s5.y + s11.x);

     
zh1.z  = M[0][0].re*sh1.z;
zh1.z -= M[0][0].im*sh1.w;
zh1.z += M[0][1].re*sh2.x;
zh1.z -= M[0][1].im*sh2.y;
zh1.z += M[0][2].re*sh2.z;
zh1.z -= M[0][2].im*sh2.w;
zh1.w  = M[0][0].re*sh1.w;
zh1.w += M[0][0].im*sh1.z;
zh1.w += M[0][1].re*sh2.y;
zh1.w += M[0][1].im*sh2.x;
zh1.w += M[0][2].re*sh2.w;
zh1.w += M[0][2].im*sh2.z;
(*(out+1)).z -= zh1.z;
(*(out+1)).w -= zh1.w;

zh2.x  = M[1][0].re*sh1.z;
zh2.x -= M[1][0].im*sh1.w;
zh2.x += M[1][1].re*sh2.x;
zh2.x -= M[1][1].im*sh2.y;
zh2.x += M[1][2].re*sh2.z;
zh2.x -= M[1][2].im*sh2.w;
zh2.y  = M[1][0].re*sh1.w;
zh2.y += M[1][0].im*sh1.z;
zh2.y += M[1][1].re*sh2.y;
zh2.y += M[1][1].im*sh2.x;
zh2.y += M[1][2].re*sh2.w;
zh2.y += M[1][2].im*sh2.z;
(*(out+2)).x -= zh2.x;
(*(out+2)).y -= zh2.y;

zh2.z  = M[2][0].re*sh1.z;
zh2.z -= M[2][0].im*sh1.w;
zh2.z += M[2][1].re*sh2.x;
zh2.z -= M[2][1].im*sh2.y;
zh2.z += M[2][2].re*sh2.z;
zh2.z -= M[2][2].im*sh2.w;
zh2.w  = M[2][0].re*sh1.w;
zh2.w += M[2][0].im*sh1.z;
zh2.w += M[2][1].re*sh2.y;
zh2.w += M[2][1].im*sh2.x;
zh2.w += M[2][2].re*sh2.w;
zh2.w += M[2][2].im*sh2.z;
(*(out+2)).z -= zh2.z;
(*(out+2)).w -= zh2.w;


//Reconstruct lower half spinor
//this is just a multiplication by "i"
     (*(out+3)).x += zh0.y;
     (*(out+3)).y -= zh0.x; 

     (*(out+3)).z += zh0.w;
     (*(out+3)).w -= zh0.z; 
     
     (*(out+4)).x += zh1.y;
     (*(out+4)).y -= zh1.x;       


//this is just a multiplication by "-i"     
     (*(out+4)).z -= zh1.w;
     (*(out+4)).w += zh1.z;       

     (*(out+5)).x -= zh2.y;
     (*(out+5)).y += zh2.x;  
     
     (*(out+5)).z -= zh2.w;
     (*(out+5)).w += zh2.z;

}







// multiplies su3-Matrix, kappa, spin-projector onto spinor
// this is +y 
// uses spin projection reduction 
// with 1-gamma_3
//
// | 1  0   0  1 |     s0      s0 + s3
// | 0  1  -1  0 |     s1      s1 - s2
// | 0 -1   1  0 |  X  s2   =  -(s1 - s2)
// | 1  0   0  1 |     s3       (s1 + s3)
//
// second half spinor is proportional to first one!!!
//
//-kappa(r - gamma_mu) kappa reell !!!!
__device__ void dev_su3MtV_kappaP2_plus_d(dev_su3_d M, const dev_spinor_d * s, double4 * out, double kappa){

  
  double4 sh0,sh1,sh2,zh0,zh1,zh2;


     sh0.x = kappa*( (*(s+0*DEVOFF)).x + (*(s+9*DEVOFF)).x);
     sh0.y = kappa*( (*(s+0*DEVOFF)).y + (*(s+9*DEVOFF)).y);
   
     sh0.z = kappa*( (*(s+1*DEVOFF)).x + (*(s+10*DEVOFF)).x);
     sh0.w = kappa*( (*(s+1*DEVOFF)).y + (*(s+10*DEVOFF)).y);
       
     sh1.x = kappa*( (*(s+2*DEVOFF)).x + (*(s+11*DEVOFF)).x);
     sh1.y = kappa*( (*(s+2*DEVOFF)).y + (*(s+11*DEVOFF)).y);
     


     
     
//Multiply by gauge field  
zh0.x =  M[0][0].re*sh0.x;
zh0.x -= M[0][0].im*sh0.y;
zh0.x += M[0][1].re*sh0.z;
zh0.x -= M[0][1].im*sh0.w;
zh0.x += M[0][2].re*sh1.x;
zh0.x -= M[0][2].im*sh1.y;
zh0.y =  M[0][0].re*sh0.y;
zh0.y += M[0][0].im*sh0.x;
zh0.y += M[0][1].re*sh0.w;
zh0.y += M[0][1].im*sh0.z;
zh0.y += M[0][2].re*sh1.y;
zh0.y += M[0][2].im*sh1.x;
(*(out+0)).x -= zh0.x;
(*(out+0)).y -= zh0.y;

zh0.z =  M[1][0].re*sh0.x;
zh0.z -= M[1][0].im*sh0.y;
zh0.z += M[1][1].re*sh0.z;
zh0.z -= M[1][1].im*sh0.w;
zh0.z += M[1][2].re*sh1.x;
zh0.z -= M[1][2].im*sh1.y;
zh0.w =  M[1][0].re*sh0.y;
zh0.w += M[1][0].im*sh0.x;
zh0.w += M[1][1].re*sh0.w;
zh0.w += M[1][1].im*sh0.z;
zh0.w += M[1][2].re*sh1.y;
zh0.w += M[1][2].im*sh1.x;
(*(out+0)).z -= zh0.z;
(*(out+0)).w -= zh0.w;

zh1.x =  M[2][0].re*sh0.x;
zh1.x -= M[2][0].im*sh0.y;
zh1.x += M[2][1].re*sh0.z;
zh1.x -= M[2][1].im*sh0.w;
zh1.x += M[2][2].re*sh1.x;
zh1.x -= M[2][2].im*sh1.y;
zh1.y =  M[2][0].re*sh0.y;
zh1.y += M[2][0].im*sh0.x;
zh1.y += M[2][1].re*sh0.w;
zh1.y += M[2][1].im*sh0.z;
zh1.y += M[2][2].re*sh1.y;
zh1.y += M[2][2].im*sh1.x;
(*(out+1)).x -= zh1.x;
(*(out+1)).y -= zh1.y;


     sh1.z = kappa*( (*(s+3*DEVOFF)).x - (*(s+6*DEVOFF)).x);
     sh1.w = kappa*( (*(s+3*DEVOFF)).y - (*(s+6*DEVOFF)).y);
     
     sh2.x = kappa*( (*(s+4*DEVOFF)).x - (*(s+7*DEVOFF)).x);
     sh2.y = kappa*( (*(s+4*DEVOFF)).y - (*(s+7*DEVOFF)).y);
          
     sh2.z = kappa*( (*(s+5*DEVOFF)).x - (*(s+8*DEVOFF)).x);
     sh2.w = kappa*( (*(s+5*DEVOFF)).y - (*(s+8*DEVOFF)).y);

     
zh1.z =  M[0][0].re*sh1.z;
zh1.z -= M[0][0].im*sh1.w;
zh1.z += M[0][1].re*sh2.x;
zh1.z -= M[0][1].im*sh2.y;
zh1.z += M[0][2].re*sh2.z;
zh1.z -= M[0][2].im*sh2.w;
zh1.w =  M[0][0].re*sh1.w;
zh1.w += M[0][0].im*sh1.z;
zh1.w += M[0][1].re*sh2.y;
zh1.w += M[0][1].im*sh2.x;
zh1.w += M[0][2].re*sh2.w;
zh1.w += M[0][2].im*sh2.z;
(*(out+1)).z -= zh1.z;
(*(out+1)).w -= zh1.w;

zh2.x =  M[1][0].re*sh1.z;
zh2.x -= M[1][0].im*sh1.w;
zh2.x += M[1][1].re*sh2.x;
zh2.x -= M[1][1].im*sh2.y;
zh2.x += M[1][2].re*sh2.z;
zh2.x -= M[1][2].im*sh2.w;
zh2.y =  M[1][0].re*sh1.w;
zh2.y += M[1][0].im*sh1.z;
zh2.y += M[1][1].re*sh2.y;
zh2.y += M[1][1].im*sh2.x;
zh2.y += M[1][2].re*sh2.w;
zh2.y += M[1][2].im*sh2.z;
(*(out+2)).x -= zh2.x;
(*(out+2)).y -= zh2.y;

zh2.z =  M[2][0].re*sh1.z;
zh2.z -= M[2][0].im*sh1.w;
zh2.z += M[2][1].re*sh2.x;
zh2.z -= M[2][1].im*sh2.y;
zh2.z += M[2][2].re*sh2.z;
zh2.z -= M[2][2].im*sh2.w;
zh2.w =  M[2][0].re*sh1.w;
zh2.w += M[2][0].im*sh1.z;
zh2.w += M[2][1].re*sh2.y;
zh2.w += M[2][1].im*sh2.x;
zh2.w += M[2][2].re*sh2.w;
zh2.w += M[2][2].im*sh2.z;
(*(out+2)).z -= zh2.z;
(*(out+2)).w -= zh2.w;


//here we use "minus component 1"
     
     (*(out+3)).x += zh1.z;
     (*(out+3)).y += zh1.w; 
     
     (*(out+3)).z += zh2.x;
     (*(out+3)).w += zh2.y;
   
     (*(out+4)).x += zh2.z;
     (*(out+4)).y += zh2.w;

     
//here we use "component 0"     
     (*(out+4)).z -= zh0.x;
     (*(out+4)).w -= zh0.y;
     
     (*(out+5)).x -= zh0.z;
     (*(out+5)).y -= zh0.w;

     (*(out+5)).z -= zh1.x;
     (*(out+5)).w -= zh1.y;     
     
}



__device__ void dev_su3MtV_kappaP2_plus_spintex_d(dev_su3_d M, int pos, double4 * out, double kappa, int updn){

  
  double4 sh0,sh1,sh2,zh0,zh1,zh2;
  double2 s0,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11;
  
  if(updn==1){
    s0 = fetch1D_spin_d(spin_d_tex,pos+0*DEVOFF);
    s1 = fetch1D_spin_d(spin_d_tex,pos+1*DEVOFF);
    s2 = fetch1D_spin_d(spin_d_tex,pos+2*DEVOFF);
    s9 = fetch1D_spin_d(spin_d_tex,pos+9*DEVOFF);
    s10 = fetch1D_spin_d(spin_d_tex,pos+10*DEVOFF);
    s11 = fetch1D_spin_d(spin_d_tex,pos+11*DEVOFF);  
  }
  else{
    s0 = fetch1D_spin_d(spin_d_tex_dn,pos+0*DEVOFF);
    s1 = fetch1D_spin_d(spin_d_tex_dn,pos+1*DEVOFF);
    s2 = fetch1D_spin_d(spin_d_tex_dn,pos+2*DEVOFF);
    s9 = fetch1D_spin_d(spin_d_tex_dn,pos+9*DEVOFF);
    s10 = fetch1D_spin_d(spin_d_tex_dn,pos+10*DEVOFF); 
    s11 = fetch1D_spin_d(spin_d_tex_dn,pos+11*DEVOFF);     
  } 

     sh0.x = kappa*( s0.x + s9.x);
     sh0.y = kappa*( s0.y + s9.y);
   
     sh0.z = kappa*( s1.x + s10.x);
     sh0.w = kappa*( s1.y + s10.y);
       
     sh1.x = kappa*( s2.x + s11.x);
     sh1.y = kappa*( s2.y + s11.y);
     


     
     
//Multiply by gauge field  
zh0.x =  M[0][0].re*sh0.x;
zh0.x -= M[0][0].im*sh0.y;
zh0.x += M[0][1].re*sh0.z;
zh0.x -= M[0][1].im*sh0.w;
zh0.x += M[0][2].re*sh1.x;
zh0.x -= M[0][2].im*sh1.y;
zh0.y =  M[0][0].re*sh0.y;
zh0.y += M[0][0].im*sh0.x;
zh0.y += M[0][1].re*sh0.w;
zh0.y += M[0][1].im*sh0.z;
zh0.y += M[0][2].re*sh1.y;
zh0.y += M[0][2].im*sh1.x;
(*(out+0)).x -= zh0.x;
(*(out+0)).y -= zh0.y;

zh0.z =  M[1][0].re*sh0.x;
zh0.z -= M[1][0].im*sh0.y;
zh0.z += M[1][1].re*sh0.z;
zh0.z -= M[1][1].im*sh0.w;
zh0.z += M[1][2].re*sh1.x;
zh0.z -= M[1][2].im*sh1.y;
zh0.w =  M[1][0].re*sh0.y;
zh0.w += M[1][0].im*sh0.x;
zh0.w += M[1][1].re*sh0.w;
zh0.w += M[1][1].im*sh0.z;
zh0.w += M[1][2].re*sh1.y;
zh0.w += M[1][2].im*sh1.x;
(*(out+0)).z -= zh0.z;
(*(out+0)).w -= zh0.w;

zh1.x =  M[2][0].re*sh0.x;
zh1.x -= M[2][0].im*sh0.y;
zh1.x += M[2][1].re*sh0.z;
zh1.x -= M[2][1].im*sh0.w;
zh1.x += M[2][2].re*sh1.x;
zh1.x -= M[2][2].im*sh1.y;
zh1.y =  M[2][0].re*sh0.y;
zh1.y += M[2][0].im*sh0.x;
zh1.y += M[2][1].re*sh0.w;
zh1.y += M[2][1].im*sh0.z;
zh1.y += M[2][2].re*sh1.y;
zh1.y += M[2][2].im*sh1.x;
(*(out+1)).x -= zh1.x;
(*(out+1)).y -= zh1.y;

  if(updn==1){
    s3 = fetch1D_spin_d(spin_d_tex,pos+3*DEVOFF);
    s4 = fetch1D_spin_d(spin_d_tex,pos+4*DEVOFF);
    s5 = fetch1D_spin_d(spin_d_tex,pos+5*DEVOFF);
    s6 = fetch1D_spin_d(spin_d_tex,pos+6*DEVOFF);
    s7 = fetch1D_spin_d(spin_d_tex,pos+7*DEVOFF);
    s8 = fetch1D_spin_d(spin_d_tex,pos+8*DEVOFF);  
  }
  else{
    s3 = fetch1D_spin_d(spin_d_tex_dn,pos+3*DEVOFF);
    s4 = fetch1D_spin_d(spin_d_tex_dn,pos+4*DEVOFF);
    s5 = fetch1D_spin_d(spin_d_tex_dn,pos+5*DEVOFF);
    s6 = fetch1D_spin_d(spin_d_tex_dn,pos+6*DEVOFF);
    s7 = fetch1D_spin_d(spin_d_tex_dn,pos+7*DEVOFF); 
    s8 = fetch1D_spin_d(spin_d_tex_dn,pos+8*DEVOFF);     
  }
  
     sh1.z = kappa*( s3.x - s6.x);
     sh1.w = kappa*( s3.y - s6.y);
     
     sh2.x = kappa*( s4.x - s7.x);
     sh2.y = kappa*( s4.y - s7.y);
          
     sh2.z = kappa*( s5.x - s8.x);
     sh2.w = kappa*( s5.y - s8.y);

     
zh1.z =  M[0][0].re*sh1.z;
zh1.z -= M[0][0].im*sh1.w;
zh1.z += M[0][1].re*sh2.x;
zh1.z -= M[0][1].im*sh2.y;
zh1.z += M[0][2].re*sh2.z;
zh1.z -= M[0][2].im*sh2.w;
zh1.w =  M[0][0].re*sh1.w;
zh1.w += M[0][0].im*sh1.z;
zh1.w += M[0][1].re*sh2.y;
zh1.w += M[0][1].im*sh2.x;
zh1.w += M[0][2].re*sh2.w;
zh1.w += M[0][2].im*sh2.z;
(*(out+1)).z -= zh1.z;
(*(out+1)).w -= zh1.w;

zh2.x =  M[1][0].re*sh1.z;
zh2.x -= M[1][0].im*sh1.w;
zh2.x += M[1][1].re*sh2.x;
zh2.x -= M[1][1].im*sh2.y;
zh2.x += M[1][2].re*sh2.z;
zh2.x -= M[1][2].im*sh2.w;
zh2.y =  M[1][0].re*sh1.w;
zh2.y += M[1][0].im*sh1.z;
zh2.y += M[1][1].re*sh2.y;
zh2.y += M[1][1].im*sh2.x;
zh2.y += M[1][2].re*sh2.w;
zh2.y += M[1][2].im*sh2.z;
(*(out+2)).x -= zh2.x;
(*(out+2)).y -= zh2.y;

zh2.z =  M[2][0].re*sh1.z;
zh2.z -= M[2][0].im*sh1.w;
zh2.z += M[2][1].re*sh2.x;
zh2.z -= M[2][1].im*sh2.y;
zh2.z += M[2][2].re*sh2.z;
zh2.z -= M[2][2].im*sh2.w;
zh2.w =  M[2][0].re*sh1.w;
zh2.w += M[2][0].im*sh1.z;
zh2.w += M[2][1].re*sh2.y;
zh2.w += M[2][1].im*sh2.x;
zh2.w += M[2][2].re*sh2.w;
zh2.w += M[2][2].im*sh2.z;
(*(out+2)).z -= zh2.z;
(*(out+2)).w -= zh2.w;


//here we use "minus component 1"
     
     (*(out+3)).x += zh1.z;
     (*(out+3)).y += zh1.w; 
     
     (*(out+3)).z += zh2.x;
     (*(out+3)).w += zh2.y;
   
     (*(out+4)).x += zh2.z;
     (*(out+4)).y += zh2.w;

     
//here we use "component 0"     
     (*(out+4)).z -= zh0.x;
     (*(out+4)).w -= zh0.y;
     
     (*(out+5)).x -= zh0.z;
     (*(out+5)).y -= zh0.w;

     (*(out+5)).z -= zh1.x;
     (*(out+5)).w -= zh1.y;     
     
}



// multiplies su3-Matrix, kappa, spin-projector onto spinor
// this is +y 
// uses spin projection reduction 
// with 1-gamma_3
//
// | 1  0   0 -1 |     s0      s0 - s3
// | 0  1   1  0 |     s1      s1 + s2
// | 0  1   1  0 |  X  s2   =  (s1 - s2)
// |-1  0   0  1 |     s3      -(s1 - s3)
//
// second half spinor is proportional to first one!!!
//
//-kappa(r + gamma_mu) kappa reell !!!!
__device__ void dev_su3MtV_kappaP2_minus_d(dev_su3_d M, const dev_spinor_d * s, double4 * out, double kappa){

  
  double4 sh0,sh1,sh2,zh0,zh1,zh2;


     sh0.x = kappa*( (*(s+0*DEVOFF)).x - (*(s+9*DEVOFF)).x);
     sh0.y = kappa*( (*(s+0*DEVOFF)).y - (*(s+9*DEVOFF)).y);
   
     sh0.z = kappa*( (*(s+1*DEVOFF)).x - (*(s+10*DEVOFF)).x);
     sh0.w = kappa*( (*(s+1*DEVOFF)).y - (*(s+10*DEVOFF)).y);
       
     sh1.x = kappa*( (*(s+2*DEVOFF)).x - (*(s+11*DEVOFF)).x);
     sh1.y = kappa*( (*(s+2*DEVOFF)).y - (*(s+11*DEVOFF)).y);
     


     
     
//Multiply by gauge field  
zh0.x =  M[0][0].re*sh0.x;
zh0.x -= M[0][0].im*sh0.y;
zh0.x += M[0][1].re*sh0.z;
zh0.x -= M[0][1].im*sh0.w;
zh0.x += M[0][2].re*sh1.x;
zh0.x -= M[0][2].im*sh1.y;
zh0.y =  M[0][0].re*sh0.y;
zh0.y += M[0][0].im*sh0.x;
zh0.y += M[0][1].re*sh0.w;
zh0.y += M[0][1].im*sh0.z;
zh0.y += M[0][2].re*sh1.y;
zh0.y += M[0][2].im*sh1.x;
(*(out+0)).x -= zh0.x;
(*(out+0)).y -= zh0.y;

zh0.z =  M[1][0].re*sh0.x;
zh0.z -= M[1][0].im*sh0.y;
zh0.z += M[1][1].re*sh0.z;
zh0.z -= M[1][1].im*sh0.w;
zh0.z += M[1][2].re*sh1.x;
zh0.z -= M[1][2].im*sh1.y;
zh0.w =  M[1][0].re*sh0.y;
zh0.w += M[1][0].im*sh0.x;
zh0.w += M[1][1].re*sh0.w;
zh0.w += M[1][1].im*sh0.z;
zh0.w += M[1][2].re*sh1.y;
zh0.w += M[1][2].im*sh1.x;
(*(out+0)).z -= zh0.z;
(*(out+0)).w -= zh0.w;

zh1.x =  M[2][0].re*sh0.x;
zh1.x -= M[2][0].im*sh0.y;
zh1.x += M[2][1].re*sh0.z;
zh1.x -= M[2][1].im*sh0.w;
zh1.x += M[2][2].re*sh1.x;
zh1.x -= M[2][2].im*sh1.y;
zh1.y =  M[2][0].re*sh0.y;
zh1.y += M[2][0].im*sh0.x;
zh1.y += M[2][1].re*sh0.w;
zh1.y += M[2][1].im*sh0.z;
zh1.y += M[2][2].re*sh1.y;
zh1.y += M[2][2].im*sh1.x;
(*(out+1)).x -= zh1.x;
(*(out+1)).y -= zh1.y;


     sh1.z = kappa*( (*(s+3*DEVOFF)).x + (*(s+6*DEVOFF)).x);
     sh1.w = kappa*( (*(s+3*DEVOFF)).y + (*(s+6*DEVOFF)).y);
     
     sh2.x = kappa*( (*(s+4*DEVOFF)).x + (*(s+7*DEVOFF)).x);
     sh2.y = kappa*( (*(s+4*DEVOFF)).y + (*(s+7*DEVOFF)).y);
          
     sh2.z = kappa*( (*(s+5*DEVOFF)).x + (*(s+8*DEVOFF)).x);
     sh2.w = kappa*( (*(s+5*DEVOFF)).y + (*(s+8*DEVOFF)).y);

     
zh1.z =  M[0][0].re*sh1.z;
zh1.z -= M[0][0].im*sh1.w;
zh1.z += M[0][1].re*sh2.x;
zh1.z -= M[0][1].im*sh2.y;
zh1.z += M[0][2].re*sh2.z;
zh1.z -= M[0][2].im*sh2.w;
zh1.w =  M[0][0].re*sh1.w;
zh1.w += M[0][0].im*sh1.z;
zh1.w += M[0][1].re*sh2.y;
zh1.w += M[0][1].im*sh2.x;
zh1.w += M[0][2].re*sh2.w;
zh1.w += M[0][2].im*sh2.z;  
(*(out+1)).z -= zh1.z;
(*(out+1)).w -= zh1.w;

zh2.x =  M[1][0].re*sh1.z;
zh2.x -= M[1][0].im*sh1.w;
zh2.x += M[1][1].re*sh2.x;
zh2.x -= M[1][1].im*sh2.y;
zh2.x += M[1][2].re*sh2.z;
zh2.x -= M[1][2].im*sh2.w;
zh2.y =  M[1][0].re*sh1.w;
zh2.y += M[1][0].im*sh1.z;
zh2.y += M[1][1].re*sh2.y;
zh2.y += M[1][1].im*sh2.x;
zh2.y += M[1][2].re*sh2.w;
zh2.y += M[1][2].im*sh2.z;
(*(out+2)).x -= zh2.x;
(*(out+2)).y -= zh2.y;

zh2.z =  M[2][0].re*sh1.z;
zh2.z -= M[2][0].im*sh1.w;
zh2.z += M[2][1].re*sh2.x;
zh2.z -= M[2][1].im*sh2.y;
zh2.z += M[2][2].re*sh2.z;
zh2.z -= M[2][2].im*sh2.w;
zh2.w =  M[2][0].re*sh1.w;
zh2.w += M[2][0].im*sh1.z;
zh2.w += M[2][1].re*sh2.y;
zh2.w += M[2][1].im*sh2.x;
zh2.w += M[2][2].re*sh2.w;
zh2.w += M[2][2].im*sh2.z;
(*(out+2)).z -= zh2.z;
(*(out+2)).w -= zh2.w;


//here we use "component 1"
     
     (*(out+3)).x -= zh1.z;
     (*(out+3)).y -= zh1.w; 
     
     (*(out+3)).z -= zh2.x;
     (*(out+3)).w -= zh2.y;
   
     (*(out+4)).x -= zh2.z;
     (*(out+4)).y -= zh2.w;

     
//here we use "minus component 0"     
     (*(out+4)).z += zh0.x;
     (*(out+4)).w += zh0.y;
     
     (*(out+5)).x += zh0.z;
     (*(out+5)).y += zh0.w;

     (*(out+5)).z += zh1.x;
     (*(out+5)).w += zh1.y;     
     
}







// multiplies su3-Matrix, kappa, spin-projector onto spinor
// this is +x 
// uses spin projection reduction 
// with 1-gamma_1
//
// | 1  0   0  i |     s0      s0 + i s3
// | 0  1   i  0 |     s1      s1 + i s2
// | 0 -i   1  0 |  X  s2   =  -i(s1 + i s2)
// |-i  0   0  1 |     s3      -i(s0 + i s3)
//
// second half spinor is proportional to first one!!!
//
//-kappa(r - gamma_mu) kappa reell !!!!
__device__ void dev_su3MtV_kappaP1_plus_d(dev_su3_d M, const dev_spinor_d * s, double4 * out, double kappa){

  
  double4 sh0,sh1,sh2,zh0,zh1,zh2;


     sh0.x = kappa*( (*(s+0*DEVOFF)).x - (*(s+9*DEVOFF)).y);
     sh0.y = kappa*( (*(s+0*DEVOFF)).y + (*(s+9*DEVOFF)).x);
     
     sh0.z = kappa*( (*(s+1*DEVOFF)).x - (*(s+10*DEVOFF)).y);
     sh0.w = kappa*( (*(s+1*DEVOFF)).y + (*(s+10*DEVOFF)).x);    

     sh1.x = kappa*((*(s+2*DEVOFF)).x - (*(s+11*DEVOFF)).y);
     sh1.y = kappa*((*(s+2*DEVOFF)).y + (*(s+11*DEVOFF)).x); 
     
     
     

     
     
//Multiply by gauge field  
zh0.x =  M[0][0].re*sh0.x;
zh0.x -= M[0][0].im*sh0.y;
zh0.x += M[0][1].re*sh0.z;
zh0.x -= M[0][1].im*sh0.w;
zh0.x += M[0][2].re*sh1.x;
zh0.x -= M[0][2].im*sh1.y;
zh0.y =  M[0][0].re*sh0.y;
zh0.y += M[0][0].im*sh0.x;
zh0.y += M[0][1].re*sh0.w;
zh0.y += M[0][1].im*sh0.z;
zh0.y += M[0][2].re*sh1.y;
zh0.y += M[0][2].im*sh1.x;
(*(out+0)).x -= zh0.x;
(*(out+0)).y -= zh0.y;

zh0.z =  M[1][0].re*sh0.x;
zh0.z -= M[1][0].im*sh0.y;
zh0.z += M[1][1].re*sh0.z;
zh0.z -= M[1][1].im*sh0.w;
zh0.z += M[1][2].re*sh1.x;
zh0.z -= M[1][2].im*sh1.y;
zh0.w =  M[1][0].re*sh0.y;
zh0.w += M[1][0].im*sh0.x;
zh0.w += M[1][1].re*sh0.w;
zh0.w += M[1][1].im*sh0.z;
zh0.w += M[1][2].re*sh1.y;
zh0.w += M[1][2].im*sh1.x;
(*(out+0)).z -= zh0.z;
(*(out+0)).w -= zh0.w;

zh1.x =  M[2][0].re*sh0.x;
zh1.x -= M[2][0].im*sh0.y;
zh1.x += M[2][1].re*sh0.z;
zh1.x -= M[2][1].im*sh0.w;
zh1.x += M[2][2].re*sh1.x;
zh1.x -= M[2][2].im*sh1.y;
zh1.y =  M[2][0].re*sh0.y;
zh1.y += M[2][0].im*sh0.x;
zh1.y += M[2][1].re*sh0.w;
zh1.y += M[2][1].im*sh0.z;
zh1.y += M[2][2].re*sh1.y;
zh1.y += M[2][2].im*sh1.x;
(*(out+1)).x -= zh1.x;
(*(out+1)).y -= zh1.y;


     sh1.z = kappa*((*(s+3*DEVOFF)).x - (*(s+6*DEVOFF)).y);
     sh1.w = kappa*((*(s+3*DEVOFF)).y + (*(s+6*DEVOFF)).x); 
     
     sh2.x = kappa*((*(s+4*DEVOFF)).x - (*(s+7*DEVOFF)).y);
     sh2.y = kappa*((*(s+4*DEVOFF)).y + (*(s+7*DEVOFF)).x);
     
     sh2.z = kappa*((*(s+5*DEVOFF)).x - (*(s+8*DEVOFF)).y);
     sh2.w = kappa*((*(s+5*DEVOFF)).y + (*(s+8*DEVOFF)).x);

     
zh1.z =  M[0][0].re*sh1.z;
zh1.z -= M[0][0].im*sh1.w;
zh1.z += M[0][1].re*sh2.x;
zh1.z -= M[0][1].im*sh2.y;
zh1.z += M[0][2].re*sh2.z;
zh1.z -= M[0][2].im*sh2.w;
zh1.w =  M[0][0].re*sh1.w;
zh1.w += M[0][0].im*sh1.z;
zh1.w += M[0][1].re*sh2.y;
zh1.w += M[0][1].im*sh2.x;
zh1.w += M[0][2].re*sh2.w;
zh1.w += M[0][2].im*sh2.z;
(*(out+1)).z -= zh1.z;
(*(out+1)).w -= zh1.w;

zh2.x =  M[1][0].re*sh1.z;
zh2.x -= M[1][0].im*sh1.w;
zh2.x += M[1][1].re*sh2.x;
zh2.x -= M[1][1].im*sh2.y;
zh2.x += M[1][2].re*sh2.z;
zh2.x -= M[1][2].im*sh2.w;
zh2.y =  M[1][0].re*sh1.w;
zh2.y += M[1][0].im*sh1.z;
zh2.y += M[1][1].re*sh2.y;
zh2.y += M[1][1].im*sh2.x;
zh2.y += M[1][2].re*sh2.w;
zh2.y += M[1][2].im*sh2.z;
(*(out+2)).x -= zh2.x;
(*(out+2)).y -= zh2.y;

zh2.z =  M[2][0].re*sh1.z;
zh2.z -= M[2][0].im*sh1.w;
zh2.z += M[2][1].re*sh2.x;
zh2.z -= M[2][1].im*sh2.y;
zh2.z += M[2][2].re*sh2.z;
zh2.z -= M[2][2].im*sh2.w;
zh2.w =  M[2][0].re*sh1.w;
zh2.w += M[2][0].im*sh1.z;
zh2.w += M[2][1].re*sh2.y;
zh2.w += M[2][1].im*sh2.x;
zh2.w += M[2][2].re*sh2.w;
zh2.w += M[2][2].im*sh2.z;
(*(out+2)).z -= zh2.z;
(*(out+2)).w -= zh2.w;


//here we use "-i times component 1"
     
     (*(out+3)).x -= zh1.w;
     (*(out+3)).y += zh1.z; 
     
     (*(out+3)).z -= zh2.y;
     (*(out+3)).w += zh2.x;
   
     (*(out+4)).x -= zh2.w;
     (*(out+4)).y += zh2.z;

     
//here we use "-i times component 0"     
     (*(out+4)).z -= zh0.y;
     (*(out+4)).w += zh0.x;
     
     (*(out+5)).x -= zh0.w;
     (*(out+5)).y += zh0.z;

     (*(out+5)).z -= zh1.y;
     (*(out+5)).w += zh1.x;     
     
}




// multiplies su3-Matrix, kappa, spin-projector onto spinor
// this is +x 
// uses spin projection reduction 
// with 1+gamma_1
//
// | 1  0   0 -i |     s0      s0 - i s3
// | 0  1  -i  0 |     s1      s1 - i s2
// | 0  i   1  0 |  X  s2   =  i(s1 - i s2)
// | i  0   0  1 |     s3      i(s0 - i s3)
//
// second half spinor is proportional to first one!!!
//
//-kappa(r + gamma_mu) kappa reell !!!!
__device__ void dev_su3MtV_kappaP1_minus_d(dev_su3_d M, const dev_spinor_d * s, double4 * out, double kappa){

  
  double4 sh0,sh1,sh2,zh0,zh1,zh2;


     sh0.x = kappa*( (*(s+0*DEVOFF)).x + (*(s+9*DEVOFF)).y);
     sh0.y = kappa*( (*(s+0*DEVOFF)).y - (*(s+9*DEVOFF)).x);
     
     sh0.z = kappa*( (*(s+1*DEVOFF)).x + (*(s+10*DEVOFF)).y);
     sh0.w = kappa*( (*(s+1*DEVOFF)).y - (*(s+10*DEVOFF)).x);    

     sh1.x = kappa*((*(s+2*DEVOFF)).x + (*(s+11*DEVOFF)).y);
     sh1.y = kappa*((*(s+2*DEVOFF)).y - (*(s+11*DEVOFF)).x); 
     
     

     
     
//Multiply by gauge field  
zh0.x =  M[0][0].re*sh0.x;
zh0.x -= M[0][0].im*sh0.y;
zh0.x += M[0][1].re*sh0.z;
zh0.x -= M[0][1].im*sh0.w;
zh0.x += M[0][2].re*sh1.x;
zh0.x -= M[0][2].im*sh1.y;
zh0.y =  M[0][0].re*sh0.y;
zh0.y += M[0][0].im*sh0.x;
zh0.y += M[0][1].re*sh0.w;
zh0.y += M[0][1].im*sh0.z;
zh0.y += M[0][2].re*sh1.y;
zh0.y += M[0][2].im*sh1.x;
(*(out+0)).x -= zh0.x;
(*(out+0)).y -= zh0.y;

zh0.z =  M[1][0].re*sh0.x;
zh0.z -= M[1][0].im*sh0.y;
zh0.z += M[1][1].re*sh0.z;
zh0.z -= M[1][1].im*sh0.w;
zh0.z += M[1][2].re*sh1.x;
zh0.z -= M[1][2].im*sh1.y;
zh0.w =  M[1][0].re*sh0.y;
zh0.w += M[1][0].im*sh0.x;
zh0.w += M[1][1].re*sh0.w;
zh0.w += M[1][1].im*sh0.z;
zh0.w += M[1][2].re*sh1.y;
zh0.w += M[1][2].im*sh1.x;
(*(out+0)).z -= zh0.z;
(*(out+0)).w -= zh0.w;

zh1.x =  M[2][0].re*sh0.x;
zh1.x -= M[2][0].im*sh0.y;
zh1.x += M[2][1].re*sh0.z;
zh1.x -= M[2][1].im*sh0.w;
zh1.x += M[2][2].re*sh1.x;
zh1.x -= M[2][2].im*sh1.y;
zh1.y =  M[2][0].re*sh0.y;
zh1.y += M[2][0].im*sh0.x;
zh1.y += M[2][1].re*sh0.w;
zh1.y += M[2][1].im*sh0.z;
zh1.y += M[2][2].re*sh1.y;
zh1.y += M[2][2].im*sh1.x;
(*(out+1)).x -= zh1.x;
(*(out+1)).y -= zh1.y;


     sh1.z = kappa*((*(s+3*DEVOFF)).x + (*(s+6*DEVOFF)).y);
     sh1.w = kappa*((*(s+3*DEVOFF)).y - (*(s+6*DEVOFF)).x); 
     
     sh2.x = kappa*((*(s+4*DEVOFF)).x + (*(s+7*DEVOFF)).y);
     sh2.y = kappa*((*(s+4*DEVOFF)).y - (*(s+7*DEVOFF)).x);
     
     sh2.z = kappa*((*(s+5*DEVOFF)).x + (*(s+8*DEVOFF)).y);
     sh2.w = kappa*((*(s+5*DEVOFF)).y - (*(s+8*DEVOFF)).x);


zh1.z =  M[0][0].re*sh1.z;
zh1.z -= M[0][0].im*sh1.w;
zh1.z += M[0][1].re*sh2.x;
zh1.z -= M[0][1].im*sh2.y;
zh1.z += M[0][2].re*sh2.z;
zh1.z -= M[0][2].im*sh2.w;
zh1.w =  M[0][0].re*sh1.w;
zh1.w += M[0][0].im*sh1.z;
zh1.w += M[0][1].re*sh2.y;
zh1.w += M[0][1].im*sh2.x;
zh1.w += M[0][2].re*sh2.w;
zh1.w += M[0][2].im*sh2.z; 
(*(out+1)).z -= zh1.z;
(*(out+1)).w -= zh1.w;

zh2.x =  M[1][0].re*sh1.z;
zh2.x -= M[1][0].im*sh1.w;
zh2.x += M[1][1].re*sh2.x;
zh2.x -= M[1][1].im*sh2.y;
zh2.x += M[1][2].re*sh2.z;
zh2.x -= M[1][2].im*sh2.w;
zh2.y =  M[1][0].re*sh1.w;
zh2.y += M[1][0].im*sh1.z;
zh2.y += M[1][1].re*sh2.y;
zh2.y += M[1][1].im*sh2.x;
zh2.y += M[1][2].re*sh2.w;
zh2.y += M[1][2].im*sh2.z;
(*(out+2)).x -= zh2.x;
(*(out+2)).y -= zh2.y;

zh2.z =  M[2][0].re*sh1.z;
zh2.z -= M[2][0].im*sh1.w;
zh2.z += M[2][1].re*sh2.x;
zh2.z -= M[2][1].im*sh2.y;
zh2.z += M[2][2].re*sh2.z;
zh2.z -= M[2][2].im*sh2.w;
zh2.w =  M[2][0].re*sh1.w;
zh2.w += M[2][0].im*sh1.z;
zh2.w += M[2][1].re*sh2.y;
zh2.w += M[2][1].im*sh2.x;
zh2.w += M[2][2].re*sh2.w;
zh2.w += M[2][2].im*sh2.z;
(*(out+2)).z -= zh2.z;
(*(out+2)).w -= zh2.w;


//here we use "i times component 1"
     
     (*(out+3)).x += zh1.w;
     (*(out+3)).y -= zh1.z; 
     
     (*(out+3)).z += zh2.y;
     (*(out+3)).w -= zh2.x;
   
     (*(out+4)).x -= -zh2.w;
     (*(out+4)).y -= zh2.z;

     
//here we use "i times component 0"     
     (*(out+4)).z += zh0.y;
     (*(out+4)).w -= zh0.x;
     
     (*(out+5)).x += zh0.w;
     (*(out+5)).y -= zh0.z;

     (*(out+5)).z += zh1.y;
     (*(out+5)).w -= zh1.x;     
     
}









__device__ void dev_Gamma5_d(double4 * in){
          (*(in+3)).x = -(*(in+3)).x;
          (*(in+3)).y = -(*(in+3)).y;
          (*(in+3)).z = -(*(in+3)).z;
          (*(in+3)).w = -(*(in+3)).w;
          (*(in+4)).x = -(*(in+4)).x;
          (*(in+4)).y = -(*(in+4)).y; 

          (*(in+4)).z = -(*(in+4)).z;
          (*(in+4)).w = -(*(in+4)).w;
          (*(in+5)).x = -(*(in+5)).x;
          (*(in+5)).y = -(*(in+5)).y;
          (*(in+5)).z = -(*(in+5)).z;
          (*(in+5)).w = -(*(in+5)).w;  
}


__device__ void dev_Gamma5_assign_d(double4* out, double4* in){
  (*(out)).x = (*(in)).x;
  (*(out)).y = (*(in)).y;
  (*(out)).z = (*(in)).z;
  (*(out)).w = (*(in)).w;
  (*(out+1)).x = (*(in+1)).x;
  (*(out+1)).y = (*(in+1)).y;

  (*(out+1)).z = (*(in+1)).z;
  (*(out+1)).w = (*(in+1)).w;
  (*(out+2)).x = (*(in+2)).x;
  (*(out+2)).y = (*(in+2)).y;
  (*(out+2)).z = (*(in+2)).z;
  (*(out+2)).w = (*(in+2)).w;

  (*(out+3)).x = -1.0*(*(in+3)).x;
  (*(out+3)).y = -1.0*(*(in+3)).y;
  (*(out+3)).z = -1.0*(*(in+3)).z;
  (*(out+3)).w = -1.0*(*(in+3)).w;
  (*(out+4)).x = -1.0*(*(in+4)).x;
  (*(out+4)).y = -1.0*(*(in+4)).y;

  (*(out+4)).z = -1.0*(*(in+4)).z;
  (*(out+4)).w = -1.0*(*(in+4)).w;
  (*(out+5)).x = -1.0*(*(in+5)).x;
  (*(out+5)).y = -1.0*(*(in+5)).y;
  (*(out+5)).z = -1.0*(*(in+5)).z;
  (*(out+5)).w = -1.0*(*(in+5)).w;
}



__device__ void dev_Gamma5_assigntoglobal_d(dev_spinor_d* out, double4* in){
  (*(out)).x = (*(in)).x;
  (*(out)).y = (*(in)).y;
  (*(out+1*DEVOFF)).x = (*(in)).z;
  (*(out+1*DEVOFF)).y = (*(in)).w;
  (*(out+2*DEVOFF)).x = (*(in+1)).x;
  (*(out+2*DEVOFF)).y = (*(in+1)).y;

  (*(out+3*DEVOFF)).x = (*(in+1)).z;
  (*(out+3*DEVOFF)).y = (*(in+1)).w;
  (*(out+4*DEVOFF)).x = (*(in+2)).x;
  (*(out+4*DEVOFF)).y = (*(in+2)).y;
  (*(out+5*DEVOFF)).x = (*(in+2)).z;
  (*(out+5*DEVOFF)).y = (*(in+2)).w;

  (*(out+6*DEVOFF)).x = -1.0*(*(in+3)).x;
  (*(out+6*DEVOFF)).y = -1.0*(*(in+3)).y;
  (*(out+7*DEVOFF)).x = -1.0*(*(in+3)).z;
  (*(out+7*DEVOFF)).y = -1.0*(*(in+3)).w;
  (*(out+8*DEVOFF)).x = -1.0*(*(in+4)).x;
  (*(out+8*DEVOFF)).y = -1.0*(*(in+4)).y;

  (*(out+9*DEVOFF)).x = -1.0*(*(in+4)).z;
  (*(out+9*DEVOFF)).y = -1.0*(*(in+4)).w;
  (*(out+10*DEVOFF)).x = -1.0*(*(in+5)).x;
  (*(out+10*DEVOFF)).y = -1.0*(*(in+5)).y;
  (*(out+11*DEVOFF)).x = -1.0*(*(in+5)).z;
  (*(out+11*DEVOFF)).y = -1.0*(*(in+5)).w;
}





// **************** double su3 operations **************************

__device__ void dev_unitsu3_d(dev_su3_d * g){
  (*g)[0][0].re = 1.0;
  (*g)[0][0].im = 0.0;
  (*g)[0][1].re = 0.0;
  (*g)[0][1].im = 0.0;
  (*g)[0][2].re = 0.0;
  (*g)[0][2].im = 0.0;
  
  (*g)[1][0].re = 0.0;
  (*g)[1][0].im = 0.0;
  (*g)[1][1].re = 1.0;
  (*g)[1][1].im = 0.0;
  (*g)[1][2].re = 0.0;
  (*g)[1][2].im = 0.0; 
   
  (*g)[2][0].re = 0.0;
  (*g)[2][0].im = 0.0;
  (*g)[2][1].re = 0.0;
  (*g)[2][1].im = 0.0;
  (*g)[2][2].re = 1.0;
  (*g)[2][2].im = 0.0;
   
}


__device__ void dev_zerosu3_d(dev_su3_d * g){
  (*g)[0][0].re = 0.0;
  (*g)[0][0].im = 0.0;
  (*g)[0][1].re = 0.0;
  (*g)[0][1].im = 0.0;
  (*g)[0][2].re = 0.0;
  (*g)[0][2].im = 0.0;
  
  (*g)[1][0].re = 0.0;
  (*g)[1][0].im = 0.0;
  (*g)[1][1].re = 0.0;
  (*g)[1][1].im = 0.0;
  (*g)[1][2].re = 0.0;
  (*g)[1][2].im = 0.0; 
   
  (*g)[2][0].re = 0.0;
  (*g)[2][0].im = 0.0;
  (*g)[2][1].re = 0.0;
  (*g)[2][1].im = 0.0;
  (*g)[2][2].re = 0.0;
  (*g)[2][2].im = 0.0;
   
}


// u = v * w
__device__ void dev_su3_ti_su3_d(dev_su3_d* u, dev_su3_d * v, dev_su3_d * w){
  dev_complex_d help1, help2;
  dev_complex_d zero = dev_initcomplex_d(0.0,0.0);
  int i,j,k;
  #pragma unroll 3
  for(i=0; i<3;i++){
    #pragma unroll 3
    for(j=0; j<3; j++){
    
      help2 = zero;
      #pragma unroll 3
      for(k=0; k<3; k++){
          help1 = dev_cmult_d((*v)[i][k],(*w)[k][j]);
          help2 = dev_cadd_d(help1, help2);
        }
        (*u)[i][j] = help2;    
    }
  }
}






// DEBUG!!!! THIS IS NOT WORKING!!!
// u = v * w
__device__ void _dev_su3_ti_su3_d(dev_su3_d* u, dev_su3_d * v, dev_su3_d * w){
  
//first row
  (*u)[0][0].re = (*v)[0][0].re* (*w)[0][0].re 
                - (*v)[0][0].im* (*w)[0][0].im
                + (*v)[0][1].re* (*w)[1][0].re 
                - (*v)[0][1].im* (*w)[1][0].im
                + (*v)[0][2].re* (*w)[2][0].re 
                - (*v)[0][2].im* (*w)[2][0].im;
  (*u)[0][0].im = (*v)[0][0].re* (*w)[0][0].im           
                + (*v)[0][0].im* (*w)[0][0].re
                + (*v)[0][1].re* (*w)[1][0].im           
                + (*v)[0][1].im* (*w)[1][0].re
                + (*v)[0][2].re* (*w)[2][0].im           
                + (*v)[0][2].im* (*w)[2][0].re;
                
  (*u)[0][1].re = (*v)[0][0].re* (*w)[0][1].re 
                - (*v)[0][0].im* (*w)[0][1].im
                + (*v)[0][1].re* (*w)[1][1].re 
                - (*v)[0][1].im* (*w)[1][1].im
                + (*v)[0][2].re* (*w)[2][1].re 
                - (*v)[0][2].im* (*w)[2][1].im;
  (*u)[0][1].im = (*v)[0][0].re* (*w)[0][1].im           
                + (*v)[0][0].im* (*w)[0][1].re
                + (*v)[0][1].re* (*w)[1][1].im           
                + (*v)[0][1].im* (*w)[1][1].re
                + (*v)[0][2].re* (*w)[2][1].im           
                + (*v)[0][2].im* (*w)[2][1].re; 
                            
  (*u)[0][2].re = (*v)[0][0].re* (*w)[0][2].re 
                - (*v)[0][0].im* (*w)[0][2].im
                + (*v)[0][1].re* (*w)[1][2].re 
                - (*v)[0][1].im* (*w)[1][2].im
                + (*v)[0][2].re* (*w)[2][2].re 
                - (*v)[0][2].im* (*w)[2][2].im;
  (*u)[0][2].im = (*v)[0][0].re* (*w)[0][2].im           
                + (*v)[0][0].im* (*w)[0][2].re
                + (*v)[0][1].re* (*w)[1][2].im           
                + (*v)[0][1].im* (*w)[1][2].re
                + (*v)[0][2].re* (*w)[2][2].im           
                + (*v)[0][2].im* (*w)[2][2].re;                            
                
                
                
                
//second row                
  (*u)[1][0].re = (*v)[1][0].re* (*w)[0][0].re 
                - (*v)[1][0].im* (*w)[0][0].im
                + (*v)[1][1].re* (*w)[1][0].re 
                - (*v)[1][1].im* (*w)[1][0].im
                + (*v)[1][2].re* (*w)[2][0].re 
                - (*v)[1][2].im* (*w)[2][0].im;
  (*u)[1][0].im = (*v)[1][0].re* (*w)[0][0].im           
                + (*v)[1][0].im* (*w)[0][0].re
                + (*v)[1][1].re* (*w)[1][0].im           
                + (*v)[1][1].im* (*w)[1][0].re
                + (*v)[1][2].re* (*w)[2][0].im           
                + (*v)[1][2].im* (*w)[2][0].re;
                
  (*u)[1][1].re = (*v)[1][0].re* (*w)[0][1].re 
                - (*v)[1][0].im* (*w)[0][1].im
                + (*v)[1][1].re* (*w)[1][1].re 
                - (*v)[1][1].im* (*w)[1][1].im
                + (*v)[1][2].re* (*w)[2][1].re 
                - (*v)[1][2].im* (*w)[2][1].im;
  (*u)[1][1].im = (*v)[1][0].re* (*w)[0][1].im           
                + (*v)[1][0].im* (*w)[0][1].re
                + (*v)[1][1].re* (*w)[1][1].im           
                + (*v)[1][1].im* (*w)[1][1].re
                + (*v)[1][2].re* (*w)[2][1].im           
                + (*v)[1][2].im* (*w)[2][1].re; 
                            
  (*u)[1][2].re = (*v)[1][0].re* (*w)[0][2].re 
                - (*v)[1][0].im* (*w)[0][2].im
                + (*v)[1][1].re* (*w)[1][2].re 
                - (*v)[1][1].im* (*w)[1][2].im
                + (*v)[1][2].re* (*w)[2][2].re 
                - (*v)[1][2].im* (*w)[2][2].im;
  (*u)[1][2].im = (*v)[1][0].re* (*w)[0][2].im           
                + (*v)[1][0].im* (*w)[0][2].re
                + (*v)[1][1].re* (*w)[1][2].im           
                + (*v)[1][1].im* (*w)[1][2].re
                + (*v)[1][2].re* (*w)[2][2].im           
                + (*v)[1][2].im* (*w)[2][2].re;                
                
                
//third row                
  (*u)[2][0].re = (*v)[2][0].re* (*w)[0][0].re 
                - (*v)[2][0].im* (*w)[0][0].im
                + (*v)[2][1].re* (*w)[1][0].re 
                - (*v)[2][1].im* (*w)[1][0].im
                + (*v)[2][2].re* (*w)[2][0].re 
                - (*v)[2][2].im* (*w)[2][0].im;
  (*u)[2][0].im = (*v)[2][0].re* (*w)[0][0].im           
                + (*v)[2][0].im* (*w)[0][0].re
                + (*v)[2][1].re* (*w)[1][0].im           
                + (*v)[2][1].im* (*w)[1][0].re
                + (*v)[2][2].re* (*w)[2][0].im           
                + (*v)[2][2].im* (*w)[2][0].re;
                
  (*u)[2][1].re = (*v)[2][0].re* (*w)[0][1].re 
                - (*v)[2][0].im* (*w)[0][1].im
                + (*v)[2][1].re* (*w)[1][1].re 
                - (*v)[2][1].im* (*w)[1][1].im
                + (*v)[2][2].re* (*w)[2][1].re 
                - (*v)[2][2].im* (*w)[2][1].im;
  (*u)[2][1].im = (*v)[2][0].re* (*w)[0][1].im           
                + (*v)[2][0].im* (*w)[0][1].re
                + (*v)[2][1].re* (*w)[1][1].im           
                + (*v)[2][1].im* (*w)[1][1].re
                + (*v)[2][2].re* (*w)[2][1].im           
                + (*v)[2][2].im* (*w)[2][1].re; 
                            
  (*u)[2][2].re = (*v)[2][0].re* (*w)[0][2].re 
                - (*v)[2][0].im* (*w)[0][2].im
                + (*v)[2][1].re* (*w)[1][2].re 
                - (*v)[2][1].im* (*w)[1][2].im
                + (*v)[2][2].re* (*w)[2][2].re 
                - (*v)[2][2].im* (*w)[2][2].im;
  (*u)[2][2].im = (*v)[2][0].re* (*w)[0][2].im           
                + (*v)[2][0].im* (*w)[0][2].re
                + (*v)[2][1].re* (*w)[1][2].im           
                + (*v)[2][1].im* (*w)[1][2].re
                + (*v)[2][2].re* (*w)[2][2].im           
                + (*v)[2][2].im* (*w)[2][2].re;                
}



// gfield[6*pos] = v * w
//only calculate the first 2 out of 3 rows
__device__ void dev_store_su3_ti_su3_d(int pos, dev_su3_2v_d* gfield , dev_su3_d * v, dev_su3_d * w){
  dev_complex_d help1, help2;
  dev_complex_d zero = dev_initcomplex_d(0.0,0.0);
  dev_su3_d u;
  int i,j,k;
  //only calculate the first 2 out of 3 rows
  #pragma unroll 2
  for(i=0; i<2;i++){
    #pragma unroll 3
    for(j=0; j<3; j++){
    
      help2 = zero;
      #pragma unroll 3
      for(k=0; k<3; k++){
          help1 = dev_cmult_d((*v)[i][k],(*w)[k][j]);
          help2 = dev_cadd_d(help1, help2);
        }
        (u)[i][j] = help2;    
    }
  }
  
  //store directly
   gfield[6*pos].x = (u)[0][0].re;
   gfield[6*pos].y = (u)[0][0].im;
   gfield[6*pos+1].x = (u)[0][1].re;
   gfield[6*pos+1].y = (u)[0][1].im;
   
   gfield[6*pos+2].x = (u)[0][2].re;
   gfield[6*pos+2].y = (u)[0][2].im;
   gfield[6*pos+3].x = (u)[1][0].re;
   gfield[6*pos+3].y = (u)[1][0].im;
   
   gfield[6*pos+4].x = (u)[1][1].re;
   gfield[6*pos+4].y = (u)[1][1].im;
   gfield[6*pos+5].x = (u)[1][2].re;
   gfield[6*pos+5].y = (u)[1][2].im;
  
}




// gfield[6*pos] = v * w
//only calculate the first 2 out of 3 rows
__device__ void _dev_store_su3_ti_su3_d(int pos, dev_su3_2v_d* gfield , dev_su3_d * v, dev_su3_d * w){
  dev_su3_d u;
  //only calculate the first 2 out of 3 rows

//first row
  (u)[0][0].re = (*v)[0][0].re* (*w)[0][0].re 
                - (*v)[0][0].im* (*w)[0][0].im
                + (*v)[0][1].re* (*w)[1][0].re 
                - (*v)[0][1].im* (*w)[1][0].im
                + (*v)[0][2].re* (*w)[2][0].re 
                - (*v)[0][2].im* (*w)[2][0].im;
  (u)[0][0].im = (*v)[0][0].re* (*w)[0][0].im           
                + (*v)[0][0].im* (*w)[0][0].re
                + (*v)[0][1].re* (*w)[1][0].im           
                + (*v)[0][1].im* (*w)[1][0].re
                + (*v)[0][2].re* (*w)[2][0].im           
                + (*v)[0][2].im* (*w)[2][0].re;
                
  (u)[0][1].re = (*v)[0][0].re* (*w)[0][1].re 
                - (*v)[0][0].im* (*w)[0][1].im
                + (*v)[0][1].re* (*w)[1][1].re 
                - (*v)[0][1].im* (*w)[1][1].im
                + (*v)[0][2].re* (*w)[2][1].re 
                - (*v)[0][2].im* (*w)[2][1].im;
  (u)[0][1].im = (*v)[0][0].re* (*w)[0][1].im           
                + (*v)[0][0].im* (*w)[0][1].re
                + (*v)[0][1].re* (*w)[1][1].im           
                + (*v)[0][1].im* (*w)[1][1].re
                + (*v)[0][2].re* (*w)[2][1].im           
                + (*v)[0][2].im* (*w)[2][1].re; 
                            
  (u)[0][2].re = (*v)[0][0].re* (*w)[0][2].re 
                - (*v)[0][0].im* (*w)[0][2].im
                + (*v)[0][1].re* (*w)[1][2].re 
                - (*v)[0][1].im* (*w)[1][2].im
                + (*v)[0][2].re* (*w)[2][2].re 
                - (*v)[0][2].im* (*w)[2][2].im;
  (u)[0][2].im = (*v)[0][0].re* (*w)[0][2].im           
                + (*v)[0][0].im* (*w)[0][2].re
                + (*v)[0][1].re* (*w)[1][2].im           
                + (*v)[0][1].im* (*w)[1][2].re
                + (*v)[0][2].re* (*w)[2][2].im           
                + (*v)[0][2].im* (*w)[2][2].re;                            
                
  //store directly
   gfield[6*pos].x = (u)[0][0].re;
   gfield[6*pos].y = (u)[0][0].im;
   gfield[6*pos+1].x = (u)[0][1].re;
   gfield[6*pos+1].y = (u)[0][1].im;
   gfield[6*pos+2].x = (u)[0][2].re;
   gfield[6*pos+2].y = (u)[0][2].im;                
                
                
//second row                
  (u)[1][0].re = (*v)[1][0].re* (*w)[0][0].re 
                - (*v)[1][0].im* (*w)[0][0].im
                + (*v)[1][1].re* (*w)[1][0].re 
                - (*v)[1][1].im* (*w)[1][0].im
                + (*v)[1][2].re* (*w)[2][0].re 
                - (*v)[1][2].im* (*w)[2][0].im;
  (u)[1][0].im = (*v)[1][0].re* (*w)[0][0].im           
                + (*v)[1][0].im* (*w)[0][0].re
                + (*v)[1][1].re* (*w)[1][0].im           
                + (*v)[1][1].im* (*w)[1][0].re
                + (*v)[1][2].re* (*w)[2][0].im           
                + (*v)[1][2].im* (*w)[2][0].re;
                
  (u)[1][1].re = (*v)[1][0].re* (*w)[0][1].re 
                - (*v)[1][0].im* (*w)[0][1].im
                + (*v)[1][1].re* (*w)[1][1].re 
                - (*v)[1][1].im* (*w)[1][1].im
                + (*v)[1][2].re* (*w)[2][1].re 
                - (*v)[1][2].im* (*w)[2][1].im;
  (u)[1][1].im = (*v)[1][0].re* (*w)[0][1].im           
                + (*v)[1][0].im* (*w)[0][1].re
                + (*v)[1][1].re* (*w)[1][1].im           
                + (*v)[1][1].im* (*w)[1][1].re
                + (*v)[1][2].re* (*w)[2][1].im           
                + (*v)[1][2].im* (*w)[2][1].re; 
                            
  (u)[1][2].re = (*v)[1][0].re* (*w)[0][2].re 
                - (*v)[1][0].im* (*w)[0][2].im
                + (*v)[1][1].re* (*w)[1][2].re 
                - (*v)[1][1].im* (*w)[1][2].im
                + (*v)[1][2].re* (*w)[2][2].re 
                - (*v)[1][2].im* (*w)[2][2].im;
  (u)[1][2].im = (*v)[1][0].re* (*w)[0][2].im           
                + (*v)[1][0].im* (*w)[0][2].re
                + (*v)[1][1].re* (*w)[1][2].im           
                + (*v)[1][1].im* (*w)[1][2].re
                + (*v)[1][2].re* (*w)[2][2].im           
                + (*v)[1][2].im* (*w)[2][2].re; 

   //store directly
   gfield[6*pos+3].x = (u)[1][0].re;
   gfield[6*pos+3].y = (u)[1][0].im;
   gfield[6*pos+4].x = (u)[1][1].re;
   gfield[6*pos+4].y = (u)[1][1].im;
   gfield[6*pos+5].x = (u)[1][2].re;
   gfield[6*pos+5].y = (u)[1][2].im;
  
}










// u = v* w^+
__device__ void dev_su3_ti_su3d_d(dev_su3_d* u, dev_su3_d * v, dev_su3_d * w){
  dev_complex_d help1, help2;
  dev_complex_d zero = dev_initcomplex_d(0.0,0.0);
  int i,j,k;
  #pragma unroll 3
  for(i=0; i<3;i++){
    #pragma unroll 3
    for(j=0; j<3; j++){ 
      help2 = zero;
      #pragma unroll 3
      for(k=0; k<3; k++){
          help1 = dev_cmult_d((*v)[i][k],dev_cconj_d( (*w)[j][k] ) );
          help2 = dev_cadd_d(help1, help2);
        }
        (*u)[i][j] = help2;   
    }
  }
}


// a = a+b
__device__ void dev_su3_add_d(dev_su3_d* a, dev_su3_d* b){
  int i,j;
  #pragma unroll 3
  for(i=0; i<3;i++){
    #pragma unroll 3
    for(j=0; j<3; j++){ 
      (*a)[i][j] = dev_cadd_d((*a)[i][j], (*b)[i][j]);
    }
  }
}

// **************** double su3 operations **************************







__device__ void dev_add_tr_lambda(dev_su3adj* r, dev_su3_d* a){
 (*r).d1+=-(*a)[1][0].im-(*a)[0][1].im; 
 (*r).d2+=+(*a)[1][0].re-(*a)[0][1].re; 
 (*r).d3+=-(*a)[0][0].im+(*a)[1][1].im; 
 (*r).d4+=-(*a)[2][0].im-(*a)[0][2].im; 
 (*r).d5+=+(*a)[2][0].re-(*a)[0][2].re; 
 (*r).d6+=-(*a)[2][1].im-(*a)[1][2].im; 
 (*r).d7+=+(*a)[2][1].re-(*a)[1][2].re; 
 (*r).d8+=(-(*a)[0][0].im-(*a)[1][1].im 
           + 2.0*(*a)[2][2].im)*0.577350269189625;
}


__device__ void dev_tr_lambda(dev_su3adj* r, dev_su3_d* a){
 (*r).d1=-(*a)[1][0].im-(*a)[0][1].im; 
 (*r).d2=+(*a)[1][0].re-(*a)[0][1].re; 
 (*r).d3=-(*a)[0][0].im+(*a)[1][1].im; 
 (*r).d4=-(*a)[2][0].im-(*a)[0][2].im; 
 (*r).d5=+(*a)[2][0].re-(*a)[0][2].re; 
 (*r).d6=-(*a)[2][1].im-(*a)[1][2].im; 
 (*r).d7=+(*a)[2][1].re-(*a)[1][2].re; 
 (*r).d8=(-(*a)[0][0].im-(*a)[1][1].im 
           + 2.0*(*a)[2][2].im)*0.577350269189625;
}



__device__ void dev_madd_tr_lambda(double d, dev_su3adj* r ,dev_su3_d *  a){
 r->d1+=(d)* (-(*a)[1][0].im-(*a)[0][1].im);
 r->d2+=(d)*(+(*a)[1][0].re-(*a)[0][1].re); 
 r->d3+=(d)*(-(*a)[0][0].im+(*a)[1][1].im); 
 r->d4+=(d)*(-(*a)[2][0].im-(*a)[0][2].im); 
 r->d5+=(d)*(+(*a)[2][0].re-(*a)[0][2].re); 
 r->d6+=(d)*(-(*a)[2][1].im-(*a)[1][2].im); 
 r->d7+=(d)*(+(*a)[2][1].re-(*a)[1][2].re); 
 r->d8+=(d)*((-(*a)[0][0].im-(*a)[1][1].im 
           + 2.0*(*a)[2][2].im)*0.577350269189625); 
}







// reconstruction of the link fields from two rows of the su3 matrix
// numbers are fetched from texture cache
__device__ void dev_reconstructgf_2vtexref_d (dev_su3_2v_d * field, int pos, dev_su3_d* gf){

  double2 gfin;
  
  gfin = field[6*pos];  
  //first row
  (*gf)[0][0].re = gfin.x;
  (*gf)[0][0].im = gfin.y;

  gfin = field[6*pos+1];  
  (*gf)[0][1].re = gfin.x;
  (*gf)[0][1].im = gfin.y;
  
  gfin = field[6*pos+2];  
  (*gf)[0][2].re = gfin.x;
  (*gf)[0][2].im = gfin.y;

  gfin = field[6*pos+3];    
  //second row
  (*gf)[1][0].re = gfin.x;
  (*gf)[1][0].im = gfin.y;
    
  gfin = field[6*pos+4];  
  (*gf)[1][1].re = gfin.x;
  (*gf)[1][1].im = gfin.y;
  
  gfin = field[6*pos+5];    
  (*gf)[1][2].re = gfin.x;
  (*gf)[1][2].im = gfin.y;

  
  
  //third row from cconj(cross product of first and second row)

  (*gf)[2][0].re = (*gf)[0][1].re * (*gf)[1][2].re;
  (*gf)[2][0].re -= (*gf)[0][1].im * (*gf)[1][2].im;
  (*gf)[2][0].re -= (*gf)[0][2].re * (*gf)[1][1].re;
  (*gf)[2][0].re += (*gf)[0][2].im * (*gf)[1][1].im;
  
  (*gf)[2][0].im = -(*gf)[0][1].re * (*gf)[1][2].im;
  (*gf)[2][0].im -= (*gf)[0][1].im * (*gf)[1][2].re;
  (*gf)[2][0].im += (*gf)[0][2].re * (*gf)[1][1].im;
  (*gf)[2][0].im += (*gf)[0][2].im * (*gf)[1][1].re;
  


  (*gf)[2][1].re = (*gf)[0][2].re * (*gf)[1][0].re;
  (*gf)[2][1].re -= (*gf)[0][2].im * (*gf)[1][0].im;
  (*gf)[2][1].re -= (*gf)[0][0].re * (*gf)[1][2].re;
  (*gf)[2][1].re += (*gf)[0][0].im * (*gf)[1][2].im;
  
  (*gf)[2][1].im = -(*gf)[0][2].re * (*gf)[1][0].im;
  (*gf)[2][1].im -= (*gf)[0][2].im * (*gf)[1][0].re;
  (*gf)[2][1].im += (*gf)[0][0].re * (*gf)[1][2].im;
  (*gf)[2][1].im += (*gf)[0][0].im * (*gf)[1][2].re;


  
  (*gf)[2][2].re = (*gf)[0][0].re * (*gf)[1][1].re;
  (*gf)[2][2].re -= (*gf)[0][0].im * (*gf)[1][1].im;
  (*gf)[2][2].re -= (*gf)[0][1].re * (*gf)[1][0].re;
  (*gf)[2][2].re += (*gf)[0][1].im * (*gf)[1][0].im;
  
  (*gf)[2][2].im = -(*gf)[0][0].re * (*gf)[1][1].im;
  (*gf)[2][2].im -= (*gf)[0][0].im * (*gf)[1][1].re;
  (*gf)[2][2].im += (*gf)[0][1].re * (*gf)[1][0].im;
  (*gf)[2][2].im += (*gf)[0][1].im * (*gf)[1][0].re;  
  
  return;
}




// su3 - dagger reconstruction from two rows  
__device__ void dev_reconstructgf_2vtexref_dagger_d (dev_su3_2v_d * field, int pos, dev_su3_d* gf){

  double2 gfin; 
  
  //first column (minus in im for complex conj.)

  gfin = field[6*pos];  
  (*gf)[0][0].re = gfin.x;
  (*gf)[0][0].im = -gfin.y;
  
  gfin = field[6*pos+1];    
  (*gf)[1][0].re = gfin.x;
  (*gf)[1][0].im = -gfin.y;

  gfin = field[6*pos+2];
  (*gf)[2][0].re = gfin.x;
  (*gf)[2][0].im = -gfin.y;


  gfin = field[6*pos+3];  
  //second  column (minus in im for complex conj.)
  (*gf)[0][1].re = gfin.x;
  (*gf)[0][1].im = -gfin.y;
  

  gfin = field[6*pos+4]; 
  (*gf)[1][1].re = gfin.x;
  (*gf)[1][1].im = -gfin.y;
  
  gfin = field[6*pos+5];   
  (*gf)[2][1].re = gfin.x;
  (*gf)[2][1].im = -gfin.y;
  
  //third column from (cross product of cconj(first column) and cconj(second column))
 
  (*gf)[0][2].re = (*gf)[1][0].re * (*gf)[2][1].re;
  (*gf)[0][2].re -= (*gf)[1][0].im * (*gf)[2][1].im;
  (*gf)[0][2].re -= (*gf)[2][0].re * (*gf)[1][1].re;
  (*gf)[0][2].re += (*gf)[2][0].im * (*gf)[1][1].im;
  
  (*gf)[0][2].im = -(*gf)[1][0].re* (*gf)[2][1].im;
  (*gf)[0][2].im -= (*gf)[1][0].im* (*gf)[2][1].re;
  (*gf)[0][2].im += (*gf)[2][0].re*(*gf)[1][1].im;
  (*gf)[0][2].im += (*gf)[2][0].im*(*gf)[1][1].re;
  

  (*gf)[1][2].re = (*gf)[2][0].re*(*gf)[0][1].re;
  (*gf)[1][2].re -= (*gf)[2][0].im*(*gf)[0][1].im;
  (*gf)[1][2].re -= (*gf)[0][0].re*(*gf)[2][1].re;
  (*gf)[1][2].re += (*gf)[0][0].im*(*gf)[2][1].im;
  
  (*gf)[1][2].im = -(*gf)[2][0].re * (*gf)[0][1].im;
  (*gf)[1][2].im -= (*gf)[2][0].im * (*gf)[0][1].re;
  (*gf)[1][2].im += (*gf)[0][0].re * (*gf)[2][1].im;
  (*gf)[1][2].im += (*gf)[0][0].im * (*gf)[2][1].re;
  
  (*gf)[2][2].re = (*gf)[0][0].re * (*gf)[1][1].re;
  (*gf)[2][2].re -= (*gf)[0][0].im * (*gf)[1][1].im;
  (*gf)[2][2].re -= (*gf)[1][0].re * (*gf)[0][1].re;
  (*gf)[2][2].re += (*gf)[1][0].im * (*gf)[0][1].im;
  
  (*gf)[2][2].im = -(*gf)[0][0].re * (*gf)[1][1].im;
  (*gf)[2][2].im -= (*gf)[0][0].im * (*gf)[1][1].re;
  (*gf)[2][2].im += (*gf)[1][0].re * (*gf)[0][1].im;
  (*gf)[2][2].im += (*gf)[1][0].im * (*gf)[0][1].re;  
  
}





/* this reconstructs the 3rd row of 3x3 gauge field from first 2 rows */
void reconstructgf_2v_host (su3* gf){
  
 
  double complex help1;
  //third row from cconj(cross product of first and second row)
  help1 = (*gf).c01 * (*gf).c12;
  help1 -= (*gf).c02 * (*gf).c11;
  (*gf).c20 = conj(help1);

  
  help1 = (*gf).c02 * (*gf).c10;
  help1 -= (*gf).c00 * (*gf).c12;
  (*gf).c21 = conj(help1);

  
  help1 = (*gf).c00 * (*gf).c11;
  help1 -= (*gf).c01 * (*gf).c10;
  (*gf).c22 = conj(help1);
  //normalize third row, too
  
  double norm = creal( 
                (*gf).c20 *conj((*gf).c20)
              + (*gf).c21 *conj((*gf).c21)
              + (*gf).c22 *conj((*gf).c22)
              );
  norm=1.0/sqrt(norm);
  
  (*gf).c20*=norm; 
  (*gf).c21*=norm;  
  (*gf).c22*=norm; 
  
  return;
}





__device__ __inline__ void dev_get_matrix_dagger(dev_su3_2v_d* gf, int pos, int dir, dev_su3_d* matrix){
   dev_reconstructgf_2vtexref_dagger_d(gf, (4*pos+dir),matrix);
}

__device__ __inline__ void dev_get_matrix(dev_su3_2v_d* gf, int pos, int dir, dev_su3_d* matrix){
   dev_reconstructgf_2vtexref_d(gf, (4*pos+dir),matrix);
}








////////////////////////  BLAS KERNELS DOUBLE ////////////////////////////////////////////////////




  int blas_gridsize_d;
  int blas_blocksize_d; // kernel parameters for the half_dot and axpy kernels
  double * dev_blas_redfield_d;  //this is the reduction field for the
                               //blas reduction kernels
  double * dev_blas_sredfield_d; // this is the small reduction field after one sweep of reduction
  double * blas_sredfield_d;                             
  int blas_redblocks_d;  // the number of blocks of the reduction kernel
                     // VOLUME/REDUCTION_N
                     // also the size of the final sum (of reduction)
                     // performed on host



void init_blas_d(int vol){
  cudaError_t cudaerr;
  
  blas_blocksize_d=BLOCK2D;
  if( vol >= BLOCK2D){
   blas_gridsize_d = (int)(vol/BLOCK2D) + 1;
  }
  else{
    blas_gridsize_d=1;
  }
  
  size_t size = vol * sizeof(double);
  
  if((cudaerr=cudaMalloc((void **) &dev_blas_redfield_d, size)) != cudaSuccess){
    if(g_proc_id==0) printf("Error in init_blas_d(): Memory allocation of double reduction field failed. Aborting...\n");
    exit(200);
  }   // Allocate array on device
  else{
    #ifndef LOWOUTPUT 
    if(g_proc_id==0) printf("Allocated blas double reduction field on device\n");
    #endif
  }  


  // IMPLEMENT THIS FOR ALL LATTICE SIZES !!!!!!!!!!!!!!!!!!!!
  if((vol%REDUCTION_N) == 0){
    blas_redblocks_d = vol/REDUCTION_N;
  }
  else{
    if(g_proc_id==0) {
      fprintf(stderr,"Error in init_blas_d(): Volume is not a multiple of REDUCTION_N (%d). Aborting...\n", REDUCTION_N);
      fprintf(stderr,"Error code is: %f\n",cudaerr);
    }
      exit(100);
  }
  
  // initialize small redfields
  size = blas_redblocks_d * sizeof(double);
  if((cudaerr=cudaMalloc((void **) &dev_blas_sredfield_d, size)) != cudaSuccess){
    if(g_proc_id==0) {
      printf("Error in init_blas_d(): Memory allocation of small double reduction field failed. Aborting...\n");
      printf("Error code is: %f\n",cudaerr);
    }
    exit(200);
  }   // Allocate array on device
  else{
    #ifndef LOWOUTPUT 
    if(g_proc_id==0) printf("Allocated blas small double reduction field on device\n");
    #endif
  }  
  
  if((void*)(blas_sredfield_d = (double *)malloc(size)) == NULL){
    if(g_proc_id==0){
      printf("Error in init_blas_d(): Could not allocate memory for double blas small redfield on host. Aborting...\n");
      printf("Error code is: %f\n",cudaerr);
    }
    exit(200);
  } 
  
  
  
}


void finalize_blas_d(){
  cudaFree(dev_blas_redfield_d);
  cudaFree(dev_blas_sredfield_d);
  free(blas_sredfield_d);
}




// this is a reduction algorithm for double based on the CUDA SDK 
__global__ void reduce_double(double *g_idata, double *g_odata, unsigned int n)
{
    extern __shared__ double sdata_d[];

    // load shared mem
    unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x*blockDim.x + threadIdx.x;
    
    sdata_d[tid] = (i < n) ? g_idata[i] : 0;
    
    __syncthreads();

    // do reduction in shared mem
    for(unsigned int s=blockDim.x/2; s>0; s>>=1) 
    {
        if (tid < s) 
        {
            sdata_d[tid] += sdata_d[tid + s];
        }
        __syncthreads();
    }

    // write result for this block to global mem
    if (tid == 0) g_odata[blockIdx.x] = sdata_d[0];
}




// this is the version for double2 
__global__ void reduce_double2(double2 *g_idata, double2 *g_odata, unsigned int n)
{
    extern __shared__ double2 sdata2_d[];

    // load shared mem
    unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x*blockDim.x + threadIdx.x;
    
    sdata2_d[tid].x = (i < n) ? g_idata[i].x : 0;
    sdata2_d[tid].y = (i < n) ? g_idata[i].y : 0;
    
    __syncthreads();

    // do reduction in shared mem
    for(unsigned int s=blockDim.x/2; s>0; s>>=1) 
    {
        if (tid < s) 
        {
            sdata2_d[tid].x += sdata2_d[tid + s].x;
            sdata2_d[tid].y += sdata2_d[tid + s].y;
        }
        __syncthreads();
    }

    // write result for this block to global mem
    if (tid == 0) {
      g_odata[blockIdx.x].x = sdata2_d[0].x;
      g_odata[blockIdx.x].y = sdata2_d[0].y;
    }
}







// y = alpha*x + y 
// x is not read from texture
// y is not read from texture
__global__ void dev_axpy_d (double alpha, dev_spinor_d* x, dev_spinor_d* y){
   int pos= threadIdx.x + blockDim.x*blockIdx.x;
   double4 xhelp[6]; 
   double4 erghelp[6];
   int i;

   
   if(pos < dev_VOLUME){
   
   //load y
   dev_read_spinor_d(&(erghelp[0]), &(y[pos]));
   //load x
   dev_read_spinor_d(&(xhelp[0]), &(x[pos]));

    
    #pragma unroll 6
    for(i=0; i<6; i++){
       erghelp[i].x += alpha*xhelp[i].x;
       erghelp[i].y += alpha*xhelp[i].y;
       erghelp[i].z += alpha*xhelp[i].z;
       erghelp[i].w += alpha*xhelp[i].w;
    }
        
     //write out spinors
     dev_write_spinor_d(&(erghelp[0]),&(y[pos])); 
   }//dev_VOLUME
}


// y = x + alpha*y 
// x is not read from texture
// y is not read from texture
__global__ void dev_xpay_d (double alpha, dev_spinor_d* x, dev_spinor_d* y){
   int pos= threadIdx.x + blockDim.x*blockIdx.x;
   double4 yhelp[6]; 
   double4 erghelp[6];
   int i;

   
   if(pos < dev_VOLUME){
   
   //load x
   dev_read_spinor_d(&(erghelp[0]), &(x[pos]));
   //load y
   dev_read_spinor_d(&(yhelp[0]), &(y[pos]));

    
    #pragma unroll 6
    for(i=0; i<6; i++){
       erghelp[i].x += alpha*yhelp[i].x;
       erghelp[i].y += alpha*yhelp[i].y;
       erghelp[i].z += alpha*yhelp[i].z;
       erghelp[i].w += alpha*yhelp[i].w;
    }
        
     //write out spinors
     dev_write_spinor_d(&(erghelp[0]),&(y[pos])); 
   }//dev_VOLUME
}








// y = alpha*y 
// y is not read from texture
__global__ void dev_blasscal_d (double alpha, dev_spinor_d* y){
   int pos= threadIdx.x + blockDim.x*blockIdx.x;
   double4 yhelp[6]; 
   double4 erghelp[6];
   int i;

   
   if(pos < dev_VOLUME){
   
   //load y
   dev_read_spinor_d(&(yhelp[0]), &(y[pos]));

    
    #pragma unroll 6
    for(i=0; i<6; i++){
       erghelp[i].x = alpha*yhelp[i].x;
       erghelp[i].y = alpha*yhelp[i].y;
       erghelp[i].z = alpha*yhelp[i].z;
       erghelp[i].w = alpha*yhelp[i].w;
    }
        
     //write out spinors
     dev_write_spinor_d(&(erghelp[0]),&(y[pos])); 
   }//dev_VOLUME
}



// y = x 
// x is not read from texture
__global__ void dev_blascopy_d (dev_spinor_d* x, dev_spinor_d* y){
   int pos= threadIdx.x + blockDim.x*blockIdx.x;
   double4 xhelp[6]; 

   
   if(pos < dev_VOLUME){
   
   //load y
     dev_read_spinor_d(&(xhelp[0]), &(x[pos]));
   //write out spinor
     dev_write_spinor_d(&(xhelp[0]),&(y[pos])); 
   }//dev_VOLUME
}




//this is a dotprod implementation for double
//x*y at spacetime point x is put into redfield at pos
__global__ void dev_dot_d( double* redfield, dev_spinor_d* x,dev_spinor_d* y){
   int pos= threadIdx.x + blockDim.x*blockIdx.x;
   double4 xhelp[6],yhelp[6];
   int i;
   double dotp = 0.0;
   
   if(pos < dev_VOLUME){
    // this is the loop over the 6 double4 forming one spinor
    
    //load y
      dev_read_spinor_d(&(yhelp[0]), &(y[pos]));
    //load x
      dev_read_spinor_d(&(xhelp[0]), &(x[pos]));

    
    #pragma unroll 6
    for(i=0; i<6; i++){
       //xhelp = tex1Dfetch(spinhalf_tex, 6*pos+i);
       //xnhelp = tex1Dfetch(spinnormhalf_tex, pos);
       
       dotp += xhelp[i].x * yhelp[i].x;
       dotp += xhelp[i].y * yhelp[i].y;
       dotp += xhelp[i].z * yhelp[i].z;
       dotp += xhelp[i].w * yhelp[i].w;
    }
    // write sum_i (x_i y_i) to reduction field 
    redfield[pos] = dotp;
   }//dev_VOLUME
}





// calculates the dot product of x and y
extern "C" double double_dotprod(dev_spinor_d* x, dev_spinor_d* y){
   int i;
   #ifdef MPI  
   double result;
   #endif
   cudaError_t cudaerr;
   
   dev_dot_d<<< blas_gridsize_d, blas_blocksize_d >>> 
                      (dev_blas_redfield_d, x, y);
   if((cudaerr=cudaGetLastError()) != cudaSuccess){
      if(g_proc_id==0) printf("%s\n", cudaGetErrorString(cudaerr));
      if(g_proc_id==0) printf("Error in double_dotprod. dev_dot_d kernel erroneous. Aborting...\n");
      exit(200);
   } 
   //reduce reductionfield on device 
   reduce_double<<< blas_redblocks_d, REDUCTION_N, 
                REDUCTION_N*sizeof(double) >>> 
                ( dev_blas_redfield_d, dev_blas_sredfield_d,  VOLUME);
   //this reduction always takes the VOLUME (also for mpi)     

   if((cudaerr=cudaGetLastError()) != cudaSuccess){
      if(g_proc_id==0) printf("%s\n", cudaGetErrorString(cudaerr));
      if(g_proc_id==0) printf("Error in double_dotprod. reduction kernel erroneous. Aborting...\n");
      exit(200);
   } 

   //copy back
   cudaMemcpy(blas_sredfield_d, dev_blas_sredfield_d, (size_t)(blas_redblocks_d*sizeof(double)), cudaMemcpyDeviceToHost);
           
   //do final reduction on host
   double finalsum=0.0;
   for(i=0; i<blas_redblocks_d; i++){
     finalsum += blas_sredfield_d[i];
   }
   #ifdef MPI
     //printf("proc %d : %f\n",g_proc_id,finalsum);
     MPI_Allreduce(&finalsum, &result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
     finalsum=result;
   #endif
   return(finalsum);
}





__global__ void dev_zero_spinor_field_d(dev_spinor_d* s1){
  int pos;
  pos= threadIdx.x + blockDim.x*blockIdx.x;  
  if(pos < dev_VOLUME){
          dev_zero_spinor_d(&(s1[pos]));
  }
}



__global__ void dev_copy_spinor_field_d(dev_spinor_d* s1, dev_spinor_d* s2){
  int pos;
  pos= threadIdx.x + blockDim.x*blockIdx.x;  
  if(pos < dev_VOLUME){
      dev_copy_spinor_d(&(s1[pos]),&(s2[pos]));
  } 
}



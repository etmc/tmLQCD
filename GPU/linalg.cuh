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
 
 
 
 
__device__ inline dev_complex dev_cconj (dev_complex c){ /*konjugiert komplexe Zahl*/
 dev_complex erg;
 erg.re = c.re;
 erg.im = -1.0f*c.im;
return erg;
}

__device__ inline void dev_ccopy(dev_complex* von, dev_complex* nach){/*kopiert complex von nach complex nach*/
  nach->re = von->re;
  nach->im = von->im;
}

__device__ inline float dev_cabssquare (dev_complex c){ /*gibt abs^2 einer komplexen Zahl zurück*/
 return c.re*c.re + c.im*c.im;
}

__device__ inline float dev_cabsolute (dev_complex c){/*gibt Betrag einer kompl. zahl zurück*/
 return sqrt(c.re*c.re + c.im*c.im);
}


__device__ inline  dev_complex dev_crealmult(dev_complex c1, float real){ /*multipliziert c1 mit reeller zahl re*/
  dev_complex erg;
  erg.re = real*c1.re;
  erg.im = real*c1.im;
return erg;
}

__device__ inline dev_complex dev_cmult (dev_complex c1, dev_complex c2){ /*multiplizier zwei komplexe Zahlen*/
  dev_complex erg;
  erg.re = c1.re * c2.re;
  erg.re -= c1.im * c2.im;
  erg.im = c1.re * c2.im;
  erg.im += c1.im * c2.re;
return erg;
}

__device__ inline dev_complex dev_mcmult (dev_complex c1, dev_complex c2){ /*minus mult zwei komplexe Zahlen*/
  dev_complex erg;
  erg.re = -c1.re * c2.re;
  erg.re += c1.im * c2.im;
  erg.im = -c1.re * c2.im;
  erg.im -= c1.im * c2.re;
return erg;
}

__device__ inline dev_complex dev_cmult_conj (dev_complex c1, dev_complex cc2){ /*out = c1*cc2^+ */
  dev_complex erg;
  erg.re = c1.re * cc2.re;
  erg.re += c1.im * cc2.im;
  erg.im = -c1.re * cc2.im;
  erg.im += c1.im * cc2.re;
return erg;
}

__device__ inline dev_complex dev_mcmult_conj (dev_complex c1, dev_complex cc2){ /* out = -c1* cc2^+ */
  dev_complex erg;
  erg.re = -c1.re * cc2.re;
  erg.re -= c1.im * cc2.im;
  erg.im = c1.re * cc2.im;
  erg.im -= c1.im * cc2.re;
return erg;
}

// a+= b*c
#define dev_cmult_addto(a, b, c)		\
  a.re += b.re*c.re;				\
  a.re -= b.im*c.im;				\
  a.im += b.re*c.im;				\
  a.im += b.im*c.re;
// a+= b*cc^+
#define dev_cmult_conj_addto(a, b, cc)		\
  a.re += b.re*cc.re;				\
  a.re += b.im*cc.im;				\
  a.im -= b.re*cc.im;				\
  a.im += b.im*cc.re;
// a-= b*c
#define dev_cmult_subfrom(a, b, c)		\
  a.re -= b.re*c.re;				\
  a.re += b.im*c.im;				\
  a.im -= b.re*c.im;				\
  a.im -= b.im*c.re;
// a-= b*cc^+
#define dev_cmult_conj_subfrom(a, b, cc)	\
  a.re -= b.re*cc.re;				\
  a.re -= b.im*cc.im;				\
  a.im += b.re*cc.im;				\
  a.im -= b.im*cc.re;



__device__ inline dev_complex dev_cadd (dev_complex c1, dev_complex c2){ /*addiert zwei komplexe Zahlen */
  dev_complex erg;
  erg.re = c1.re;
  erg.re += c2.re;
  erg.im = c1.im;
  erg.im += c2.im;
return erg;
}


__device__ inline dev_complex dev_cdiv(dev_complex c1, dev_complex c2) { /* dividiert c1 durch c2 */
  dev_complex erg;
  float oneovernenner = 1.0f/(c2.re*c2.re + c2.im*c2.im);
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


__device__ inline dev_complex dev_initcomplex(float re, float im){/* gibt komplexe Zahl mit Realt re und Imt im zurück*/
    dev_complex erg;
    erg.re = re;
    erg.im = im;
return (erg);
}












// this is to copy a spinor in global mem
__device__ inline void dev_copy_spinor(dev_spinor *i1, dev_spinor *i2){
  int i;
  #pragma unroll 6
  for(i=0;i<6;i++){ //color + spin
    (*(i2+i*DEVOFF)).x = (*(i1+i*DEVOFF)).x;
    (*(i2+i*DEVOFF)).y = (*(i1+i*DEVOFF)).y;
    (*(i2+i*DEVOFF)).z = (*(i1+i*DEVOFF)).z;
    (*(i2+i*DEVOFF)).w = (*(i1+i*DEVOFF)).w;
  }
}


// this is to copy a spinor in local mem
__device__ inline void dev_copy_spinor_local(dev_spinor *i1, dev_spinor *i2){
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
__device__ inline void dev_write_spinor(dev_spinor *i1, dev_spinor *i2){


    (*(i2+0*DEVOFF)).x = (*(i1+0)).x;
    (*(i2+0*DEVOFF)).y = (*(i1+0)).y;
    (*(i2+0*DEVOFF)).z = (*(i1+0)).z;
    (*(i2+0*DEVOFF)).w = (*(i1+0)).w;

    (*(i2+1*DEVOFF)).x = (*(i1+1)).x;
    (*(i2+1*DEVOFF)).y = (*(i1+1)).y;
    (*(i2+1*DEVOFF)).z = (*(i1+1)).z;
    (*(i2+1*DEVOFF)).w = (*(i1+1)).w;
    
    (*(i2+2*DEVOFF)).x = (*(i1+2)).x;
    (*(i2+2*DEVOFF)).y = (*(i1+2)).y;
    (*(i2+2*DEVOFF)).z = (*(i1+2)).z;
    (*(i2+2*DEVOFF)).w = (*(i1+2)).w;
    
    (*(i2+3*DEVOFF)).x = (*(i1+3)).x;
    (*(i2+3*DEVOFF)).y = (*(i1+3)).y;
    (*(i2+3*DEVOFF)).z = (*(i1+3)).z;
    (*(i2+3*DEVOFF)).w = (*(i1+3)).w;
    
    (*(i2+4*DEVOFF)).x = (*(i1+4)).x;
    (*(i2+4*DEVOFF)).y = (*(i1+4)).y;
    (*(i2+4*DEVOFF)).z = (*(i1+4)).z;
    (*(i2+4*DEVOFF)).w = (*(i1+4)).w;
    
    (*(i2+5*DEVOFF)).x = (*(i1+5)).x;
    (*(i2+5*DEVOFF)).y = (*(i1+5)).y;
    (*(i2+5*DEVOFF)).z = (*(i1+5)).z;
    (*(i2+5*DEVOFF)).w = (*(i1+5)).w;    
    
    
}


// this is to read a spinor from global into local mem (24 numbers adjacent)
__device__ inline void dev_read_spinor(dev_spinor *i1, dev_spinor *i2){



    (*(i1+0)).x  = (*(i2+0*DEVOFF)).x;
    (*(i1+0)).y  = (*(i2+0*DEVOFF)).y;
    (*(i1+0)).z  = (*(i2+0*DEVOFF)).z;
    (*(i1+0)).w  = (*(i2+0*DEVOFF)).w;

    (*(i1+1)).x  = (*(i2+1*DEVOFF)).x;
    (*(i1+1)).y  = (*(i2+1*DEVOFF)).y;
    (*(i1+1)).z  = (*(i2+1*DEVOFF)).z;
    (*(i1+1)).w  = (*(i2+1*DEVOFF)).w;
    
    (*(i1+2)).x  = (*(i2+2*DEVOFF)).x;
    (*(i1+2)).y  = (*(i2+2*DEVOFF)).y;
    (*(i1+2)).z  = (*(i2+2*DEVOFF)).z;
    (*(i1+2)).w  = (*(i2+2*DEVOFF)).w; 
    
    (*(i1+3)).x  = (*(i2+3*DEVOFF)).x;
    (*(i1+3)).y  = (*(i2+3*DEVOFF)).y;
    (*(i1+3)).z  = (*(i2+3*DEVOFF)).z;
    (*(i1+3)).w  = (*(i2+3*DEVOFF)).w;

    (*(i1+4)).x  = (*(i2+4*DEVOFF)).x;
    (*(i1+4)).y  = (*(i2+4*DEVOFF)).y;
    (*(i1+4)).z  = (*(i2+4*DEVOFF)).z;
    (*(i1+4)).w  = (*(i2+4*DEVOFF)).w;
    
    (*(i1+5)).x  = (*(i2+5*DEVOFF)).x;
    (*(i1+5)).y  = (*(i2+5*DEVOFF)).y;
    (*(i1+5)).z  = (*(i2+5*DEVOFF)).z;
    (*(i1+5)).w  = (*(i2+5*DEVOFF)).w;     
    
}




// this is to read a spinor from texture into local mem (24 numbers adjacent)
__device__ inline void dev_read_spinor_tex_up(dev_spinor *i1, int pos){

dev_spinor s0,s1,s2,s3,s4,s5;

  s0 = tex1Dfetch(spin_tex0,pos);
  s1 = tex1Dfetch(spin_tex1,pos);
  s2 = tex1Dfetch(spin_tex2,pos);
  s3 = tex1Dfetch(spin_tex3,pos);
  s4 = tex1Dfetch(spin_tex4,pos); 
  s5 = tex1Dfetch(spin_tex5,pos); 

    (*(i1+0)).x  = s0.x;
    (*(i1+0)).y  = s0.y;
    (*(i1+0)).z  = s0.z;
    (*(i1+0)).w  = s0.w;

    (*(i1+1)).x  = s1.x;
    (*(i1+1)).y  = s1.y;
    (*(i1+1)).z  = s1.z;
    (*(i1+1)).w  = s1.w;
    
    (*(i1+2)).x  = s2.x;
    (*(i1+2)).y  = s2.y;
    (*(i1+2)).z  = s2.z;
    (*(i1+2)).w  = s2.w; 
    
    (*(i1+3)).x  = s3.x;
    (*(i1+3)).y  = s3.y;
    (*(i1+3)).z  = s3.z;
    (*(i1+3)).w  = s3.w;

    (*(i1+4)).x  = s4.x;
    (*(i1+4)).y  = s4.y;
    (*(i1+4)).z  = s4.z;
    (*(i1+4)).w  = s4.w;
    
    (*(i1+5)).x  = s5.x;
    (*(i1+5)).y  = s5.y;
    (*(i1+5)).z  = s5.z;
    (*(i1+5)).w  = s5.w;     
    
}



// this is to read a spinor from texture into local mem (24 numbers adjacent)
__device__ inline void dev_read_spinor_tex_dn(dev_spinor *i1, int pos){

dev_spinor s0,s1,s2,s3,s4,s5;

  s0 = tex1Dfetch(spin_tex_dn0,pos);
  s1 = tex1Dfetch(spin_tex_dn1,pos);
  s2 = tex1Dfetch(spin_tex_dn2,pos);
  s3 = tex1Dfetch(spin_tex_dn3,pos);
  s4 = tex1Dfetch(spin_tex_dn4,pos); 
  s5 = tex1Dfetch(spin_tex_dn5,pos); 

    (*(i1+0)).x  = s0.x;
    (*(i1+0)).y  = s0.y;
    (*(i1+0)).z  = s0.z;
    (*(i1+0)).w  = s0.w;

    (*(i1+1)).x  = s1.x;
    (*(i1+1)).y  = s1.y;
    (*(i1+1)).z  = s1.z;
    (*(i1+1)).w  = s1.w;
    
    (*(i1+2)).x  = s2.x;
    (*(i1+2)).y  = s2.y;
    (*(i1+2)).z  = s2.z;
    (*(i1+2)).w  = s2.w; 
    
    (*(i1+3)).x  = s3.x;
    (*(i1+3)).y  = s3.y;
    (*(i1+3)).z  = s3.z;
    (*(i1+3)).w  = s3.w;

    (*(i1+4)).x  = s4.x;
    (*(i1+4)).y  = s4.y;
    (*(i1+4)).z  = s4.z;
    (*(i1+4)).w  = s4.w;
    
    (*(i1+5)).x  = s5.x;
    (*(i1+5)).y  = s5.y;
    (*(i1+5)).z  = s5.z;
    (*(i1+5)).w  = s5.w;     
    
}








// this is to zero a spinor in global mem
__device__ inline void dev_zero_spinor(dev_spinor *sin){
  int i;
  #pragma unroll 6
  for(i=0;i<6;i++){ //color + spin
    (*(sin+i*DEVOFF)).x = 0.0f;
    (*(sin+i*DEVOFF)).y = 0.0f;
    (*(sin+i*DEVOFF)).z = 0.0f;
    (*(sin+i*DEVOFF)).w = 0.0f;
  }
}


// this is to zero a spinor in local mem
__device__ inline void dev_zero_spinor_local(dev_spinor *sin){
  int i;
  #pragma unroll 6
  for(i=0;i<6;i++){ //color + spin
    (*(sin+i)).x = 0.0f;
    (*(sin+i)).y = 0.0f;
    (*(sin+i)).z = 0.0f;
    (*(sin+i)).w = 0.0f;
  }
}


//out = in + lambda in2
__device__ inline void dev_skalarmult_add_assign_spinor(dev_spinor *in, float lambda,dev_spinor * in2, dev_spinor * out){
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


//out(global mem) = in + lambda in2(local mem)
//have to care about different indices!!
__device__ inline void dev_complexmult_add_assign_writetoglobal_spinor(dev_spinor * in, dev_complex lambda,dev_spinor * in2, dev_spinor * out){
  int i;
  #pragma unroll 6
  for(i=0;i<6;i++){ //color + spin
    (*(out+i*DEVOFF)).x = (*(in+i)).x + ((*(in2+i)).x*lambda.re - (*(in2+i)).y*lambda.im);
    (*(out+i*DEVOFF)).y = (*(in+i)).y + ((*(in2+i)).x*lambda.im + (*(in2+i)).y*lambda.re);
    (*(out+i*DEVOFF)).z = (*(in+i)).z + ((*(in2+i)).z*lambda.re - (*(in2+i)).w*lambda.im);
    (*(out+i*DEVOFF)).w = (*(in+i)).w + ((*(in2+i)).z*lambda.im + (*(in2+i)).w*lambda.re);
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




//out(global mem) = in + (lambda)* in2 (local mem)
//have to care about different indices!!
__device__ inline void dev_complexcgmult_add_assign_writetoglobal_spinor(dev_spinor * in, dev_complex lambda,dev_spinor * in2, dev_spinor * out){
  int i;
  #pragma unroll 6
  for(i=0;i<6;i++){ //color + spin
    (*(out+i*DEVOFF)).x = (*(in+i)).x + ((*(in2+i)).x*lambda.re + (*(in2+i)).y*lambda.im);
    (*(out+i*DEVOFF)).y = (*(in+i)).y + (-(*(in2+i)).x*lambda.im + (*(in2+i)).y*lambda.re);
    (*(out+i*DEVOFF)).z = (*(in+i)).z + ((*(in2+i)).z*lambda.re + (*(in2+i)).w*lambda.im);
    (*(out+i*DEVOFF)).w = (*(in+i)).w + (-(*(in2+i)).z*lambda.im + (*(in2+i)).w*lambda.re);
  }
}



__device__ void inline dev_skalarmult_spinor(dev_spinor * in, dev_complex lambda, dev_spinor * out){
  int i;
  #pragma unroll 6
  for(i=0;i<6;i++){ //color + spin
    //out[i] = dev_cmult(in[i],lambda);
    
    (*(out+i*DEVOFF)).x = (*(in+i*DEVOFF)).x*lambda.re - (*(in+i*DEVOFF)).y*lambda.im;
    (*(out+i*DEVOFF)).y = (*(in+i*DEVOFF)).y*lambda.re + (*(in+i*DEVOFF)).x*lambda.im;
    
    (*(out+i*DEVOFF)).z = (*(in+i*DEVOFF)).z*lambda.re - (*(in+i*DEVOFF)).w*lambda.im;
    (*(out+i*DEVOFF)).w = (*(in+i*DEVOFF)).w*lambda.re + (*(in+i*DEVOFF)).z*lambda.im;
  }
}




__device__ void inline dev_skalarmult_gamma5_spinor(dev_spinor * out, dev_complex lambda, dev_spinor * in){

 dev_spinor shelp, tempout;

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






__device__ void inline dev_skalarmult_gamma5_globalspinor(dev_spinor * out, dev_complex lambda, dev_spinor * in){

 dev_spinor shelp, tempout;

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
 
 
shelp = *(in+1*DEVOFF);
 tempout.x = shelp.x*lambda.re;
 tempout.x -= shelp.y*lambda.im;
 
 tempout.y = shelp.y*lambda.re;
 tempout.y += shelp.x*lambda.im;
 
 tempout.z = shelp.z*lambda.re;
 tempout.z -= shelp.w*lambda.im;

 tempout.w = shelp.w*lambda.re;
 tempout.w += shelp.z*lambda.im;
(*(out+1)) = tempout; 
 
 
shelp = *(in+2*DEVOFF);
 tempout.x = shelp.x*lambda.re;
 tempout.x -= shelp.y*lambda.im;

 tempout.y = shelp.y*lambda.re;
 tempout.y += shelp.x*lambda.im;
 
 tempout.z = shelp.z*lambda.re;
 tempout.z -= shelp.w*lambda.im;

 tempout.w = shelp.w*lambda.re;
 tempout.w += shelp.z*lambda.im;
(*(out+2)) = tempout; 

 
shelp = *(in+3*DEVOFF);
 tempout.x = shelp.y*lambda.im;
 tempout.x -= shelp.x*lambda.re;

 tempout.y = - shelp.x*lambda.im;
 tempout.y -= shelp.y*lambda.re;
 
 tempout.z = shelp.w*lambda.im;
 tempout.z -= shelp.z*lambda.re;

 tempout.w = -shelp.z*lambda.im;
 tempout.w -= shelp.w*lambda.re;
(*(out+3)) = tempout;

shelp = *(in+4*DEVOFF);
 tempout.x = shelp.y*lambda.im;
 tempout.x -= shelp.x*lambda.re;

 tempout.y = - shelp.x*lambda.im;
 tempout.y -= shelp.y*lambda.re;
 
 tempout.z = shelp.w*lambda.im;
 tempout.z -= shelp.z*lambda.re;

 tempout.w = -shelp.z*lambda.im;
 tempout.w -= shelp.w*lambda.re;
(*(out+4)) = tempout;

shelp = *(in+5*DEVOFF);
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













//this is in the relativistic basis
//gamma5 NOT diagonal
//we read a spinor from global mem and store in lacal mem!
__device__ void inline dev_skalarmult_gamma5_globalspinor_rel(dev_spinor * out, dev_complex lambda, dev_spinor * in){


 (*(out+3)).x  = -(*(in)).x*lambda.re;
 (*(out+3)).x += (*(in)).y*lambda.im; 
 (*(out+3)).y  = -(*(in)).y*lambda.re;
 (*(out+3)).y -= (*(in)).x*lambda.im; 
 (*(out+3)).z  = -(*(in)).z*lambda.re;
 (*(out+3)).z += (*(in)).w*lambda.im;
 (*(out+3)).w  = -(*(in)).w*lambda.re;
 (*(out+3)).w -= (*(in)).z*lambda.im;
 (*(out+4)).x  = -(*(in+1*DEVOFF)).x*lambda.re;
 (*(out+4)).x += (*(in+1*DEVOFF)).y*lambda.im; 
 (*(out+4)).y  = -(*(in+1*DEVOFF)).y*lambda.re;
 (*(out+4)).y -= (*(in+1*DEVOFF)).x*lambda.im;
 

 (*(out+4)).z  = -(*(in+1*DEVOFF)).z*lambda.re;
 (*(out+4)).z += (*(in+1*DEVOFF)).w*lambda.im;
 (*(out+4)).w  = -(*(in+1*DEVOFF)).w*lambda.re;
 (*(out+4)).w -= (*(in+1*DEVOFF)).z*lambda.im;
 (*(out+5)).x  = -(*(in+2*DEVOFF)).x*lambda.re;
 (*(out+5)).x += (*(in+2*DEVOFF)).y*lambda.im; 
 (*(out+5)).y  = -(*(in+2*DEVOFF)).y*lambda.re;
 (*(out+5)).y -= (*(in+2*DEVOFF)).x*lambda.im;
 (*(out+5)).z  = -(*(in+2*DEVOFF)).z*lambda.re;
 (*(out+5)).z += (*(in+2*DEVOFF)).w*lambda.im;
 (*(out+5)).w  = -(*(in+2*DEVOFF)).w*lambda.re;
 (*(out+5)).w -= (*(in+2*DEVOFF)).z*lambda.im;

 //this is the lower spinor
 //re and im parts are interchanged w.r.t. above
 //same sign as above 
 (*(out+0)).x  = (*(in+3*DEVOFF)).y*lambda.im;
 (*(out+0)).x -= (*(in+3*DEVOFF)).x*lambda.re;
 (*(out+0)).y  = -(*(in+3*DEVOFF)).x*lambda.im;
 (*(out+0)).y -= (*(in+3*DEVOFF)).y*lambda.re; 
 (*(out+0)).z  = (*(in+3*DEVOFF)).w*lambda.im;
 (*(out+0)).z -= (*(in+3*DEVOFF)).z*lambda.re;
 (*(out+0)).w  = -(*(in+3*DEVOFF)).z*lambda.im;
 (*(out+0)).w -= (*(in+3*DEVOFF)).w*lambda.re;
 (*(out+1)).x  = (*(in+4*DEVOFF)).y*lambda.im;
 (*(out+1)).x -= (*(in+4*DEVOFF)).x*lambda.re;
 (*(out+1)).y  = -(*(in+4*DEVOFF)).x*lambda.im;
 (*(out+1)).y -= (*(in+4*DEVOFF)).y*lambda.re;
 
 
 
 (*(out+1)).z  = (*(in+4*DEVOFF)).w*lambda.im;
 (*(out+1)).z -= (*(in+4*DEVOFF)).z*lambda.re;
 (*(out+1)).w  = -(*(in+4*DEVOFF)).z*lambda.im;
 (*(out+1)).w -= (*(in+4*DEVOFF)).w*lambda.re; 
 (*(out+2)).x  = (*(in+5*DEVOFF)).y*lambda.im;
 (*(out+2)).x -= (*(in+5*DEVOFF)).x*lambda.re;
 (*(out+2)).y  = -(*(in+5*DEVOFF)).x*lambda.im;
 (*(out+2)).y -= (*(in+5*DEVOFF)).y*lambda.re; 
 (*(out+2)).z  = (*(in+5*DEVOFF)).w*lambda.im;
 (*(out+2)).z -= (*(in+5*DEVOFF)).z*lambda.re;
 (*(out+2)).w  = -(*(in+5*DEVOFF)).z*lambda.im;
 (*(out+2)).w -= (*(in+5*DEVOFF)).w*lambda.re;

}


//this is in the relativistic basis
//gamma5 NOT diagonal
__device__ void inline dev_skalarmult_gamma5_spinor_rel(dev_spinor * out, dev_complex lambda, dev_spinor * in){


 (*(out+3)).x  = -(*(in)).x*lambda.re;
 (*(out+3)).x += (*(in)).y*lambda.im; 
 (*(out+3)).y  = -(*(in)).y*lambda.re;
 (*(out+3)).y -= (*(in)).x*lambda.im; 
 (*(out+3)).z  = -(*(in)).z*lambda.re;
 (*(out+3)).z += (*(in)).w*lambda.im;
 (*(out+3)).w  = -(*(in)).w*lambda.re;
 (*(out+3)).w -= (*(in)).z*lambda.im;
 (*(out+4)).x  = -(*(in+1)).x*lambda.re;
 (*(out+4)).x += (*(in+1)).y*lambda.im; 
 (*(out+4)).y  = -(*(in+1)).y*lambda.re;
 (*(out+4)).y -= (*(in+1)).x*lambda.im;
 

 (*(out+4)).z  = -(*(in+1)).z*lambda.re;
 (*(out+4)).z += (*(in+1)).w*lambda.im;
 (*(out+4)).w  = -(*(in+1)).w*lambda.re;
 (*(out+4)).w -= (*(in+1)).z*lambda.im;
 (*(out+5)).x  = -(*(in+2)).x*lambda.re;
 (*(out+5)).x += (*(in+2)).y*lambda.im; 
 (*(out+5)).y  = -(*(in+2)).y*lambda.re;
 (*(out+5)).y -= (*(in+2)).x*lambda.im;
 (*(out+5)).z  = -(*(in+2)).z*lambda.re;
 (*(out+5)).z += (*(in+2)).w*lambda.im;
 (*(out+5)).w  = -(*(in+2)).w*lambda.re;
 (*(out+5)).w -= (*(in+2)).z*lambda.im;

 //this is the lower spinor
 //re and im parts are interchanged w.r.t. above
 //same sign as above 
 (*(out+0)).x  = (*(in+3)).y*lambda.im;
 (*(out+0)).x -= (*(in+3)).x*lambda.re;
 (*(out+0)).y  = -(*(in+3)).x*lambda.im;
 (*(out+0)).y -= (*(in+3)).y*lambda.re; 
 (*(out+0)).z  = (*(in+3)).w*lambda.im;
 (*(out+0)).z -= (*(in+3)).z*lambda.re;
 (*(out+0)).w  = -(*(in+3)).z*lambda.im;
 (*(out+0)).w -= (*(in+3)).w*lambda.re;
 (*(out+1)).x  = (*(in+4)).y*lambda.im;
 (*(out+1)).x -= (*(in+4)).x*lambda.re;
 (*(out+1)).y  = -(*(in+4)).x*lambda.im;
 (*(out+1)).y -= (*(in+4)).y*lambda.re;
 
 
 
 (*(out+1)).z  = (*(in+4)).w*lambda.im;
 (*(out+1)).z -= (*(in+4)).z*lambda.re;
 (*(out+1)).w  = -(*(in+4)).z*lambda.im;
 (*(out+1)).w -= (*(in+4)).w*lambda.re; 
 (*(out+2)).x  = (*(in+5)).y*lambda.im;
 (*(out+2)).x -= (*(in+5)).x*lambda.re;
 (*(out+2)).y  = -(*(in+5)).x*lambda.im;
 (*(out+2)).y -= (*(in+5)).y*lambda.re; 
 (*(out+2)).z  = (*(in+5)).w*lambda.im;
 (*(out+2)).z -= (*(in+5)).z*lambda.re;
 (*(out+2)).w  = -(*(in+5)).z*lambda.im;
 (*(out+2)).w -= (*(in+5)).w*lambda.re;

}





__device__ void inline dev_realmult_spinor(dev_spinor * in, float lambda){
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


__device__ void inline dev_realmult_spinor_assign(dev_spinor* out, float lambda, dev_spinor* in){
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


__device__ void inline dev_realmult_spinor_assigntoglobal(dev_spinor* out, float lambda, dev_spinor* in){
int i;

  //out[i] = in[i]*lambda;
      (*(out+0*DEVOFF)).x = (*(in+0)).x*lambda;
      (*(out+0*DEVOFF)).y = (*(in+0)).y*lambda;     
      (*(out+0*DEVOFF)).z = (*(in+0)).z*lambda;
      (*(out+0*DEVOFF)).w = (*(in+0)).w*lambda;
  
      (*(out+1*DEVOFF)).x = (*(in+1)).x*lambda;
      (*(out+1*DEVOFF)).y = (*(in+1)).y*lambda;     
      (*(out+1*DEVOFF)).z = (*(in+1)).z*lambda;
      (*(out+1*DEVOFF)).w = (*(in+1)).w*lambda;
      
      (*(out+2*DEVOFF)).x = (*(in+2)).x*lambda;
      (*(out+2*DEVOFF)).y = (*(in+2)).y*lambda;     
      (*(out+2*DEVOFF)).z = (*(in+2)).z*lambda;
      (*(out+2*DEVOFF)).w = (*(in+2)).w*lambda;
      
      (*(out+3*DEVOFF)).x = (*(in+3)).x*lambda;
      (*(out+3*DEVOFF)).y = (*(in+3)).y*lambda;     
      (*(out+3*DEVOFF)).z = (*(in+3)).z*lambda;
      (*(out+3*DEVOFF)).w = (*(in+3)).w*lambda;
      
      (*(out+4*DEVOFF)).x = (*(in+4)).x*lambda;
      (*(out+4*DEVOFF)).y = (*(in+4)).y*lambda;     
      (*(out+4*DEVOFF)).z = (*(in+4)).z*lambda;
      (*(out+4*DEVOFF)).w = (*(in+4)).w*lambda;
      
      (*(out+5*DEVOFF)).x = (*(in+5)).x*lambda;
      (*(out+5*DEVOFF)).y = (*(in+5)).y*lambda;     
      (*(out+5*DEVOFF)).z = (*(in+5)).z*lambda;
      (*(out+5*DEVOFF)).w = (*(in+5)).w*lambda;     
  
}

__device__ void dev_assign_realmult_add_spinor(dev_spinor* out, float lambda, dev_spinor* in1,  dev_spinor* in2){
int i;
float help;
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


__device__ inline void dev_add_spinor_assign(dev_spinor * i1, dev_spinor * i2){
  int i;

    (*(i1+0)).x = (*(i1+0)).x + (*(i2+0)).x;
    (*(i1+0)).y = (*(i1+0)).y + (*(i2+0)).y;
    (*(i1+0)).z = (*(i1+0)).z + (*(i2+0)).z;
    (*(i1+0)).w = (*(i1+0)).w + (*(i2+0)).w;
    
    (*(i1+1)).x = (*(i1+1)).x + (*(i2+1)).x;
    (*(i1+1)).y = (*(i1+1)).y + (*(i2+1)).y;
    (*(i1+1)).z = (*(i1+1)).z + (*(i2+1)).z;
    (*(i1+1)).w = (*(i1+1)).w + (*(i2+1)).w;
    
    (*(i1+2)).x = (*(i1+2)).x + (*(i2+2)).x;
    (*(i1+2)).y = (*(i1+2)).y + (*(i2+2)).y;
    (*(i1+2)).z = (*(i1+2)).z + (*(i2+2)).z;
    (*(i1+2)).w = (*(i1+2)).w + (*(i2+2)).w;

    (*(i1+3)).x = (*(i1+3)).x + (*(i2+3)).x;
    (*(i1+3)).y = (*(i1+3)).y + (*(i2+3)).y;
    (*(i1+3)).z = (*(i1+3)).z + (*(i2+3)).z;
    (*(i1+3)).w = (*(i1+3)).w + (*(i2+3)).w;
    
    (*(i1+4)).x = (*(i1+4)).x + (*(i2+4)).x;
    (*(i1+4)).y = (*(i1+4)).y + (*(i2+4)).y;
    (*(i1+4)).z = (*(i1+4)).z + (*(i2+4)).z;
    (*(i1+4)).w = (*(i1+4)).w + (*(i2+4)).w;    
    
    (*(i1+5)).x = (*(i1+5)).x + (*(i2+5)).x;
    (*(i1+5)).y = (*(i1+5)).y + (*(i2+5)).y;
    (*(i1+5)).z = (*(i1+5)).z + (*(i2+5)).z;
    (*(i1+5)).w = (*(i1+5)).w + (*(i2+5)).w;  
}


__device__ inline void dev_add_globalspinor_assign(dev_spinor * i1, dev_spinor * i2){
  int i;

    (*(i1+0)).x = (*(i1+0)).x + (*(i2+0*DEVOFF)).x;
    (*(i1+0)).y = (*(i1+0)).y + (*(i2+0*DEVOFF)).y;
    (*(i1+0)).z = (*(i1+0)).z + (*(i2+0*DEVOFF)).z;
    (*(i1+0)).w = (*(i1+0)).w + (*(i2+0*DEVOFF)).w;

    (*(i1+1)).x = (*(i1+1)).x + (*(i2+1*DEVOFF)).x;
    (*(i1+1)).y = (*(i1+1)).y + (*(i2+1*DEVOFF)).y;
    (*(i1+1)).z = (*(i1+1)).z + (*(i2+1*DEVOFF)).z;
    (*(i1+1)).w = (*(i1+1)).w + (*(i2+1*DEVOFF)).w;    

    (*(i1+2)).x = (*(i1+2)).x + (*(i2+2*DEVOFF)).x;
    (*(i1+2)).y = (*(i1+2)).y + (*(i2+2*DEVOFF)).y;
    (*(i1+2)).z = (*(i1+2)).z + (*(i2+2*DEVOFF)).z;
    (*(i1+2)).w = (*(i1+2)).w + (*(i2+2*DEVOFF)).w;
    
    (*(i1+3)).x = (*(i1+3)).x + (*(i2+3*DEVOFF)).x;
    (*(i1+3)).y = (*(i1+3)).y + (*(i2+3*DEVOFF)).y;
    (*(i1+3)).z = (*(i1+3)).z + (*(i2+3*DEVOFF)).z;
    (*(i1+3)).w = (*(i1+3)).w + (*(i2+3*DEVOFF)).w;
    
    (*(i1+4)).x = (*(i1+4)).x + (*(i2+4*DEVOFF)).x;
    (*(i1+4)).y = (*(i1+4)).y + (*(i2+4*DEVOFF)).y;
    (*(i1+4)).z = (*(i1+4)).z + (*(i2+4*DEVOFF)).z;
    (*(i1+4)).w = (*(i1+4)).w + (*(i2+4*DEVOFF)).w;
    
    (*(i1+5)).x = (*(i1+5)).x + (*(i2+5*DEVOFF)).x;
    (*(i1+5)).y = (*(i1+5)).y + (*(i2+5*DEVOFF)).y;
    (*(i1+5)).z = (*(i1+5)).z + (*(i2+5*DEVOFF)).z;
    (*(i1+5)).w = (*(i1+5)).w + (*(i2+5*DEVOFF)).w;
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

__device__ inline void dev_sub_globalspinor_assign(dev_spinor * i1, dev_spinor * i2){
  int i;
  #pragma unroll 6
  for(i=0;i<6;i++){ //color + spin
    (*(i1+i)).x = (*(i1+i)).x - (*(i2+i*DEVOFF)).x;
    (*(i1+i)).y = (*(i1+i)).y - (*(i2+i*DEVOFF)).y;
    (*(i1+i)).z = (*(i1+i)).z - (*(i2+i*DEVOFF)).z;
    (*(i1+i)).w = (*(i1+i)).w - (*(i2+i*DEVOFF)).w;
  }
}




// this is to store a spinor in local mem (24 numbers adjacent) to global mem
__device__ inline void dev_write_spinor_d2(dev_spinor_d *i1, dev_spinor_d *i2){
  int i;
  #pragma unroll 6
  for(i=0;i<6;i++){ //color + spin
    (*(i2+i*DEVOFF)).x = (*(i1+i)).x;
    (*(i2+i*DEVOFF)).y = (*(i1+i)).y;
    (*(i2+i*DEVOFF)).z = (*(i1+i)).z;
    (*(i2+i*DEVOFF)).w = (*(i1+i)).w;
  }
}


// this is to store a spinor in local mem (24 numbers adjacent) to global mem
__device__ inline void dev_write_spinor_d2f(dev_spinor_d *i1, dev_spinor *i2){
  int i;
  #pragma unroll 6
  for(i=0;i<6;i++){ //color + spin
    (*(i2+i*DEVOFF)).x = __double2float_rn((*(i1+i)).x);
    (*(i2+i*DEVOFF)).y = __double2float_rn((*(i1+i)).y);
    (*(i2+i*DEVOFF)).z = __double2float_rn((*(i1+i)).z);
    (*(i2+i*DEVOFF)).w = __double2float_rn((*(i1+i)).w);
  }
}



// this is to store a spinor in local mem (24 numbers adjacent) to global mem
__device__ inline void dev_write_spinor_f2d(dev_spinor *i1, dev_spinor_d *i2){
  int i;
  #pragma unroll 6
  for(i=0;i<6;i++){ //color + spin
    (*(i2+i*DEVOFF)).x = (double) (*(i1+i)).x;
    (*(i2+i*DEVOFF)).y = (double) (*(i1+i)).y;
    (*(i2+i*DEVOFF)).z = (double) (*(i1+i)).z;
    (*(i2+i*DEVOFF)).w = (double) (*(i1+i)).w;
  }
}



// this is to read a spinor from global into local mem (24 numbers adjacent)
__device__ inline void dev_read_spinor_d2(dev_spinor_d *i1, dev_spinor_d *i2){
  int i;
  #pragma unroll 6
  for(i=0;i<6;i++){ //color + spin
    (*(i1+i)).x  = (*(i2+i*DEVOFF)).x;
    (*(i1+i)).y  = (*(i2+i*DEVOFF)).y;
    (*(i1+i)).z  = (*(i2+i*DEVOFF)).z;
    (*(i1+i)).w  = (*(i2+i*DEVOFF)).w;
  }
}




__global__ void dev_d2f(dev_spinor* spinfloat, dev_spinor_d* spindouble){

   int pos;
   pos= threadIdx.x + blockDim.x*blockIdx.x;  
   
   double4 help[6]; 
   
   if(pos < dev_VOLUME){
   
   dev_read_spinor_d2(&(help[0]), &(spindouble[pos]));
   dev_write_spinor_d2f(&(help[0]),&(spinfloat[pos])); 


//        spinfloat[6*pos+0].x = 1.0f;
//        spinfloat[6*pos+0].y = 1.0f;
//        spinfloat[6*pos+0].z = 1.0f;
//        spinfloat[6*pos+0].w = 1.0f;      
// 
//        spinfloat[6*pos+1].x = 1.0f;
//        spinfloat[6*pos+1].y = 1.0f;
//        spinfloat[6*pos+1].z = 1.0f;
//        spinfloat[6*pos+1].w = 1.0f;
//        
//        spinfloat[6*pos+2].x = 1.0f;
//        spinfloat[6*pos+2].y = 1.0f;
//        spinfloat[6*pos+2].z = 1.0f;
//        spinfloat[6*pos+2].w = 1.0f;  
//        
//        spinfloat[6*pos+3].x = 1.0f;
//        spinfloat[6*pos+3].y = 1.0f;
//        spinfloat[6*pos+3].z = 1.0f;
//        spinfloat[6*pos+3].w = 1.0f; 
//        
//        spinfloat[6*pos+4].x = 1.0f;
//        spinfloat[6*pos+4].y = 1.0f;
//        spinfloat[6*pos+4].z = 1.0f;
//        spinfloat[6*pos+4].w = 1.0f; 
//        
//        spinfloat[6*pos+5].x = 1.0f;
//        spinfloat[6*pos+5].y = 1.0f;
//        spinfloat[6*pos+5].z = 1.0f;
//        spinfloat[6*pos+5].w = 1.0f; 
       
/*       spinfloat[6*pos+0].x = __double2float_rn(spindouble[6*pos+0].x);
       spinfloat[6*pos+0].y = __double2float_rn(spindouble[6*pos+0].y);
       spinfloat[6*pos+0].z = __double2float_rn(spindouble[6*pos+0].z);
       spinfloat[6*pos+0].w = __double2float_rn(spindouble[6*pos+0].w);      

       spinfloat[6*pos+1].x = __double2float_rn(spindouble[6*pos+1].x);
       spinfloat[6*pos+1].y = __double2float_rn(spindouble[6*pos+1].y);
       spinfloat[6*pos+1].z = __double2float_rn(spindouble[6*pos+1].z);
       spinfloat[6*pos+1].w = __double2float_rn(spindouble[6*pos+1].w);
       
       spinfloat[6*pos+2].x = __double2float_rn(spindouble[6*pos+2].x);
       spinfloat[6*pos+2].y = __double2float_rn(spindouble[6*pos+2].y);
       spinfloat[6*pos+2].z = __double2float_rn(spindouble[6*pos+2].z);
       spinfloat[6*pos+2].w = __double2float_rn(spindouble[6*pos+2].w);  
       
       spinfloat[6*pos+3].x = __double2float_rn(spindouble[6*pos+3].x);
       spinfloat[6*pos+3].y = __double2float_rn(spindouble[6*pos+3].y);
       spinfloat[6*pos+3].z = __double2float_rn(spindouble[6*pos+3].z);
       spinfloat[6*pos+3].w = __double2float_rn(spindouble[6*pos+3].w); 
       
       spinfloat[6*pos+4].x = __double2float_rn(spindouble[6*pos+4].x);
       spinfloat[6*pos+4].y = __double2float_rn(spindouble[6*pos+4].y);
       spinfloat[6*pos+4].z = __double2float_rn(spindouble[6*pos+4].z);
       spinfloat[6*pos+4].w = __double2float_rn(spindouble[6*pos+4].w); 
       
       spinfloat[6*pos+5].x = __double2float_rn(spindouble[6*pos+5].x);
       spinfloat[6*pos+5].y = __double2float_rn(spindouble[6*pos+5].y);
       spinfloat[6*pos+5].z = __double2float_rn(spindouble[6*pos+5].z);
       spinfloat[6*pos+5].w = __double2float_rn(spindouble[6*pos+5].w);      */  
       
       
   }
   
}


__global__ void dev_f2d(dev_spinor_d* spindouble, dev_spinor* spinfloat){

   int pos;
   pos= threadIdx.x + blockDim.x*blockIdx.x;  
   
   float4 help[6]; 
   
   if(pos < dev_VOLUME){
     dev_read_spinor(&(help[0]), &(spinfloat[pos]));
     dev_write_spinor_f2d(&(help[0]),&(spindouble[pos])); 
   }
   
}






// out = x + (float) y 
// x is not read from texture
// y is not read from texture
__global__ void dev_add_f2d (dev_spinor_d* out, dev_spinor_d* x, dev_spinor* y){
   int pos= threadIdx.x + blockDim.x*blockIdx.x;
   double4 xhelp[6]; 
   float4 yhelp[6]; 
   double4 erghelp[6];
   int i;

   if(pos < dev_VOLUME){
   
   //load x
   dev_read_spinor_d2(&(xhelp[0]), &(x[pos]));
   //load y
   dev_read_spinor(&(yhelp[0]), &(y[pos]));

    #pragma unroll 6
    for(i=0; i<6; i++){
       erghelp[i].x = xhelp[i].x + (double) yhelp[i].x;
       erghelp[i].y = xhelp[i].y + (double) yhelp[i].y;
       erghelp[i].z = xhelp[i].z + (double) yhelp[i].z;
       erghelp[i].w = xhelp[i].w + (double) yhelp[i].w;
    }
        
     //write out spinors
     dev_write_spinor_d2(&(erghelp[0]),&(out[pos])); 
   }//dev_VOLUME
}


// out = x - y 
// x is not read from texture
// y is not read from texture
__global__ void dev_diff_d (dev_spinor_d* out, dev_spinor_d* x, dev_spinor_d* y){
   int pos= threadIdx.x + blockDim.x*blockIdx.x;
   double4 xhelp[6];
   double4 yhelp[6];   
   double4 erghelp[6];
   int i;

   if(pos < dev_VOLUME){
   
   //load x
   dev_read_spinor_d2(&(xhelp[0]), &(x[pos]));
   //load y
   dev_read_spinor_d2(&(yhelp[0]), &(y[pos]));

    #pragma unroll 6
    for(i=0; i<6; i++){
       erghelp[i].x = xhelp[i].x - yhelp[i].x;
       erghelp[i].y = xhelp[i].y - yhelp[i].y;
       erghelp[i].z = xhelp[i].z - yhelp[i].z;
       erghelp[i].w = xhelp[i].w - yhelp[i].w;
    }
        
     //write out spinors
     dev_write_spinor_d2(&(erghelp[0]),&(out[pos])); 
   }//dev_VOLUME
}


/*
//multipliziert su3-Matrix mal Spinor im Dirac-Raum
//code in su3_MtV.txt -- generated with codegen
__device__ void dev_su3MtV_spintex(dev_su3 M, int pos, dev_spinor * out){

dev_spinor s1, s2;

s1 = tex1Dfetch(spin_tex,6*pos);
s2 = tex1Dfetch(spin_tex,6*pos+1);


//(*(out+0)).x =  ( M[0][0].re*s1.x - M[0][0].im*s1.y ) + ( M[0][1].re*s1.z - M[0][1].im*s1.w ) + ( M[0][2].re*s2.x - M[0][2].im*s2.y );
//(*(out+0)).y = ( M[0][0].re*s1.y + M[0][0].im*s1.x ) + ( M[0][1].re*s1.w + M[0][1].im*s1.z ) + ( M[0][2].re*s2.y + M[0][2].im*s2.x );

(*(out+0)).x =  M[0][0].re*s1.x;
  (*(out+0)).y =  M[0][0].re*s1.y;
(*(out+0)).x -= M[0][0].im*s1.y;
  (*(out+0)).y += M[0][0].im*s1.x;
(*(out+0)).x += M[0][1].re*s1.z;
  (*(out+0)).y += M[0][1].re*s1.w;
(*(out+0)).x -= M[0][1].im*s1.w;
  (*(out+0)).y += M[0][1].im*s1.z;
(*(out+0)).x += M[0][2].re*s2.x;
  (*(out+0)).y += M[0][2].re*s2.y;
(*(out+0)).x -= M[0][2].im*s2.y;
  (*(out+0)).y += M[0][2].im*s2.x;



//(*(out+0)).z =  ( M[1][0].re*s1.x - M[1][0].im*s1.y ) + ( M[1][1].re*s1.z - M[1][1].im*s1.w ) + ( M[1][2].re*s2.x - M[1][2].im*s2.y );
//(*(out+0)).w =  ( M[1][0].re*s1.y + M[1][0].im*s1.x ) + ( M[1][1].re*s1.w + M[1][1].im*s1.z ) + ( M[1][2].re*s2.y + M[1][2].im*s2.x );


(*(out+0)).z = M[1][0].re*s1.x;
  (*(out+0)).w = M[1][0].re*s1.y;
(*(out+0)).z -= M[1][0].im*s1.y;
  (*(out+0)).w += M[1][0].im*s1.x;
(*(out+0)).z += M[1][1].re*s1.z;
  (*(out+0)).w += M[1][1].re*s1.w;
(*(out+0)).z -= M[1][1].im*s1.w;
  (*(out+0)).w += M[1][1].im*s1.z;
(*(out+0)).z += M[1][2].re*s2.x;
  (*(out+0)).w += M[1][2].re*s2.y;
(*(out+0)).z -= M[1][2].im*s2.y;
  (*(out+0)).w += M[1][2].im*s2.x;



//(*(out+1)).x = ( M[2][0].re*s1.x - M[2][0].im*s1.y ) + ( M[2][1].re*s1.z - M[2][1].im*s1.w ) + ( M[2][2].re*s2.x - M[2][2].im*s2.y );
//(*(out+1)).y =  ( M[2][0].re*s1.y + M[2][0].im*s1.x ) + ( M[2][1].re*s1.w + M[2][1].im*s1.z ) + ( M[2][2].re*s2.y + M[2][2].im*s2.x );


(*(out+1)).x = M[2][0].re*s1.x;
  (*(out+1)).y = M[2][0].re*s1.y;
(*(out+1)).x -= M[2][0].im*s1.y;
  (*(out+1)).y += M[2][0].im*s1.x;
(*(out+1)).x += M[2][1].re*s1.z;
  (*(out+1)).y += M[2][1].re*s1.w;
(*(out+1)).x -= M[2][1].im*s1.w;
  (*(out+1)).y += M[2][1].im*s1.z;
(*(out+1)).x += M[2][2].re*s2.x;
  (*(out+1)).y += M[2][2].re*s2.y;
(*(out+1)).x -= M[2][2].im*s2.y;
  (*(out+1)).y += M[2][2].im*s2.x;





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
*/





//multipliziert su3-Matrix mal Spinor im Dirac-Raum
//code in su3_MtV.txt -- generated with codegen
__device__ void dev_su3MtV_spintex(dev_su3 M, int pos, dev_spinor * out){

dev_spinor s1, s2;

#ifndef HALF
 s1 = tex1Dfetch(spin_tex0,pos);
#else
 s1 = tex1Dfetch(spinhalf_tex0,pos);
 float norm = tex1Dfetch(spinnormhalf_tex,pos);
 s1.x *= norm; 
 s1.y *= norm; 
 s1.z *= norm;
 s1.w *= norm;
#endif

#ifndef HALF
 s2 = tex1Dfetch(spin_tex1,pos);
#else
 s2 = tex1Dfetch(spinhalf_tex1,pos);
 s2.x *= norm; 
 s2.y *= norm; 
 s2.z *= norm; 
 s2.w *= norm;
#endif

(*(out+0)).x =  ( M[0][0].re*s1.x - M[0][0].im*s1.y ) + ( M[0][1].re*s1.z - M[0][1].im*s1.w ) + ( M[0][2].re*s2.x - M[0][2].im*s2.y );
(*(out+0)).y = ( M[0][0].re*s1.y + M[0][0].im*s1.x ) + ( M[0][1].re*s1.w + M[0][1].im*s1.z ) + ( M[0][2].re*s2.y + M[0][2].im*s2.x );



(*(out+0)).z =  ( M[1][0].re*s1.x - M[1][0].im*s1.y ) + ( M[1][1].re*s1.z - M[1][1].im*s1.w ) + ( M[1][2].re*s2.x - M[1][2].im*s2.y );
(*(out+0)).w =  ( M[1][0].re*s1.y + M[1][0].im*s1.x ) + ( M[1][1].re*s1.w + M[1][1].im*s1.z ) + ( M[1][2].re*s2.y + M[1][2].im*s2.x );


(*(out+1)).x = ( M[2][0].re*s1.x - M[2][0].im*s1.y ) + ( M[2][1].re*s1.z - M[2][1].im*s1.w ) + ( M[2][2].re*s2.x - M[2][2].im*s2.y );
(*(out+1)).y =  ( M[2][0].re*s1.y + M[2][0].im*s1.x ) + ( M[2][1].re*s1.w + M[2][1].im*s1.z ) + ( M[2][2].re*s2.y + M[2][2].im*s2.x );


#ifndef HALF
 s1 = tex1Dfetch(spin_tex2,pos);
#else
 s1 = tex1Dfetch(spinhalf_tex2,pos);
 s1.x *= norm; 
 s1.y *= norm; 
 s1.z *= norm;
 s1.w *= norm;
#endif

(*(out+1)).z =  ( M[0][0].re*s2.z - M[0][0].im*s2.w ) + ( M[0][1].re*s1.x - M[0][1].im*s1.y ) + ( M[0][2].re*s1.z - M[0][2].im*s1.w );
(*(out+1)).w =  ( M[0][0].re*s2.w + M[0][0].im*s2.z ) + ( M[0][1].re*s1.y + M[0][1].im*s1.x ) + ( M[0][2].re*s1.w + M[0][2].im*s1.z );


(*(out+2)).x = ( M[1][0].re*s2.z - M[1][0].im*s2.w ) + ( M[1][1].re*s1.x - M[1][1].im*s1.y ) + ( M[1][2].re*s1.z - M[1][2].im*s1.w );
(*(out+2)).y =  ( M[1][0].re*s2.w + M[1][0].im*s2.z ) + ( M[1][1].re*s1.y + M[1][1].im*s1.x ) + ( M[1][2].re*s1.w + M[1][2].im*s1.z );


(*(out+2)).z =  ( M[2][0].re*s2.z - M[2][0].im*s2.w ) + ( M[2][1].re*s1.x - M[2][1].im*s1.y ) + ( M[2][2].re*s1.z - M[2][2].im*s1.w );
(*(out+2)).w =  ( M[2][0].re*s2.w + M[2][0].im*s2.z ) + ( M[2][1].re*s1.y + M[2][1].im*s1.x ) + ( M[2][2].re*s1.w + M[2][2].im*s1.z );


#ifndef HALF
 s1 = tex1Dfetch(spin_tex3,pos);
#else
 s1 = tex1Dfetch(spinhalf_tex3,pos);
 s1.x *= norm; 
 s1.y *= norm; 
 s1.z *= norm;
 s1.w *= norm;
#endif

#ifndef HALF
 s2 = tex1Dfetch(spin_tex4,pos);
#else
 s2 = tex1Dfetch(spinhalf_tex4,pos);
 s2.x *= norm; 
 s2.y *= norm; 
 s2.z *= norm; 
 s2.w *= norm;
#endif
(*(out+3)).x =  ( M[0][0].re*s1.x - M[0][0].im*s1.y ) + ( M[0][1].re*s1.z - M[0][1].im*s1.w ) + ( M[0][2].re*s2.x - M[0][2].im*s2.y );
(*(out+3)).y =   ( M[0][0].re*s1.y + M[0][0].im*s1.x ) + ( M[0][1].re*s1.w + M[0][1].im*s1.z ) + ( M[0][2].re*s2.y + M[0][2].im*s2.x );


(*(out+3)).z =  ( M[1][0].re*s1.x - M[1][0].im*s1.y ) + ( M[1][1].re*s1.z - M[1][1].im*s1.w ) + ( M[1][2].re*s2.x - M[1][2].im*s2.y );
(*(out+3)).w =  ( M[1][0].re*s1.y + M[1][0].im*s1.x ) + ( M[1][1].re*s1.w + M[1][1].im*s1.z ) + ( M[1][2].re*s2.y + M[1][2].im*s2.x );


(*(out+4)).x =  ( M[2][0].re*s1.x - M[2][0].im*s1.y ) + ( M[2][1].re*s1.z - M[2][1].im*s1.w ) + ( M[2][2].re*s2.x - M[2][2].im*s2.y );
(*(out+4)).y =  ( M[2][0].re*s1.y + M[2][0].im*s1.x ) + ( M[2][1].re*s1.w + M[2][1].im*s1.z ) + ( M[2][2].re*s2.y + M[2][2].im*s2.x );


#ifndef HALF
 s1 = tex1Dfetch(spin_tex5,pos);
#else
 s1 = tex1Dfetch(spinhalf_tex5,pos);
 s1.x *= norm; 
 s1.y *= norm; 
 s1.z *= norm;
 s1.w *= norm;
#endif
(*(out+4)).z =  ( M[0][0].re*s2.z - M[0][0].im*s2.w ) + ( M[0][1].re*s1.x - M[0][1].im*s1.y ) + ( M[0][2].re*s1.z - M[0][2].im*s1.w );
(*(out+4)).w =   ( M[0][0].re*s2.w + M[0][0].im*s2.z ) + ( M[0][1].re*s1.y + M[0][1].im*s1.x ) + ( M[0][2].re*s1.w + M[0][2].im*s1.z );


(*(out+5)).x = ( M[1][0].re*s2.z - M[1][0].im*s2.w ) + ( M[1][1].re*s1.x - M[1][1].im*s1.y ) + ( M[1][2].re*s1.z - M[1][2].im*s1.w );
(*(out+5)).y =  ( M[1][0].re*s2.w + M[1][0].im*s2.z ) + ( M[1][1].re*s1.y + M[1][1].im*s1.x ) + ( M[1][2].re*s1.w + M[1][2].im*s1.z );


(*(out+5)).z =  ( M[2][0].re*s2.z - M[2][0].im*s2.w ) + ( M[2][1].re*s1.x - M[2][1].im*s1.y ) + ( M[2][2].re*s1.z - M[2][2].im*s1.w );
(*(out+5)).w =  ( M[2][0].re*s2.w + M[2][0].im*s2.z ) + ( M[2][1].re*s1.y + M[2][1].im*s1.x ) + ( M[2][2].re*s1.w + M[2][2].im*s1.z );


}


__device__ void dev_su3MtV_spintex2(dev_su3 M, int pos, dev_spinor * out){

dev_spinor s1, s2;


 s1 = tex1Dfetch(spin_tex_dn0,pos);
 s2 = tex1Dfetch(spin_tex_dn1,pos);


(*(out+0)).x =  ( M[0][0].re*s1.x - M[0][0].im*s1.y ) + ( M[0][1].re*s1.z - M[0][1].im*s1.w ) + ( M[0][2].re*s2.x - M[0][2].im*s2.y );
(*(out+0)).y = ( M[0][0].re*s1.y + M[0][0].im*s1.x ) + ( M[0][1].re*s1.w + M[0][1].im*s1.z ) + ( M[0][2].re*s2.y + M[0][2].im*s2.x );



(*(out+0)).z =  ( M[1][0].re*s1.x - M[1][0].im*s1.y ) + ( M[1][1].re*s1.z - M[1][1].im*s1.w ) + ( M[1][2].re*s2.x - M[1][2].im*s2.y );
(*(out+0)).w =  ( M[1][0].re*s1.y + M[1][0].im*s1.x ) + ( M[1][1].re*s1.w + M[1][1].im*s1.z ) + ( M[1][2].re*s2.y + M[1][2].im*s2.x );


(*(out+1)).x = ( M[2][0].re*s1.x - M[2][0].im*s1.y ) + ( M[2][1].re*s1.z - M[2][1].im*s1.w ) + ( M[2][2].re*s2.x - M[2][2].im*s2.y );
(*(out+1)).y =  ( M[2][0].re*s1.y + M[2][0].im*s1.x ) + ( M[2][1].re*s1.w + M[2][1].im*s1.z ) + ( M[2][2].re*s2.y + M[2][2].im*s2.x );



 s1 = tex1Dfetch(spin_tex_dn2,pos);


(*(out+1)).z =  ( M[0][0].re*s2.z - M[0][0].im*s2.w ) + ( M[0][1].re*s1.x - M[0][1].im*s1.y ) + ( M[0][2].re*s1.z - M[0][2].im*s1.w );
(*(out+1)).w =  ( M[0][0].re*s2.w + M[0][0].im*s2.z ) + ( M[0][1].re*s1.y + M[0][1].im*s1.x ) + ( M[0][2].re*s1.w + M[0][2].im*s1.z );


(*(out+2)).x = ( M[1][0].re*s2.z - M[1][0].im*s2.w ) + ( M[1][1].re*s1.x - M[1][1].im*s1.y ) + ( M[1][2].re*s1.z - M[1][2].im*s1.w );
(*(out+2)).y =  ( M[1][0].re*s2.w + M[1][0].im*s2.z ) + ( M[1][1].re*s1.y + M[1][1].im*s1.x ) + ( M[1][2].re*s1.w + M[1][2].im*s1.z );


(*(out+2)).z =  ( M[2][0].re*s2.z - M[2][0].im*s2.w ) + ( M[2][1].re*s1.x - M[2][1].im*s1.y ) + ( M[2][2].re*s1.z - M[2][2].im*s1.w );
(*(out+2)).w =  ( M[2][0].re*s2.w + M[2][0].im*s2.z ) + ( M[2][1].re*s1.y + M[2][1].im*s1.x ) + ( M[2][2].re*s1.w + M[2][2].im*s1.z );



 s1 = tex1Dfetch(spin_tex_dn3,pos);
 s2 = tex1Dfetch(spin_tex_dn4,pos);

(*(out+3)).x =  ( M[0][0].re*s1.x - M[0][0].im*s1.y ) + ( M[0][1].re*s1.z - M[0][1].im*s1.w ) + ( M[0][2].re*s2.x - M[0][2].im*s2.y );
(*(out+3)).y =   ( M[0][0].re*s1.y + M[0][0].im*s1.x ) + ( M[0][1].re*s1.w + M[0][1].im*s1.z ) + ( M[0][2].re*s2.y + M[0][2].im*s2.x );


(*(out+3)).z =  ( M[1][0].re*s1.x - M[1][0].im*s1.y ) + ( M[1][1].re*s1.z - M[1][1].im*s1.w ) + ( M[1][2].re*s2.x - M[1][2].im*s2.y );
(*(out+3)).w =  ( M[1][0].re*s1.y + M[1][0].im*s1.x ) + ( M[1][1].re*s1.w + M[1][1].im*s1.z ) + ( M[1][2].re*s2.y + M[1][2].im*s2.x );


(*(out+4)).x =  ( M[2][0].re*s1.x - M[2][0].im*s1.y ) + ( M[2][1].re*s1.z - M[2][1].im*s1.w ) + ( M[2][2].re*s2.x - M[2][2].im*s2.y );
(*(out+4)).y =  ( M[2][0].re*s1.y + M[2][0].im*s1.x ) + ( M[2][1].re*s1.w + M[2][1].im*s1.z ) + ( M[2][2].re*s2.y + M[2][2].im*s2.x );



 s1 = tex1Dfetch(spin_tex_dn5,pos);

(*(out+4)).z =  ( M[0][0].re*s2.z - M[0][0].im*s2.w ) + ( M[0][1].re*s1.x - M[0][1].im*s1.y ) + ( M[0][2].re*s1.z - M[0][2].im*s1.w );
(*(out+4)).w =   ( M[0][0].re*s2.w + M[0][0].im*s2.z ) + ( M[0][1].re*s1.y + M[0][1].im*s1.x ) + ( M[0][2].re*s1.w + M[0][2].im*s1.z );


(*(out+5)).x = ( M[1][0].re*s2.z - M[1][0].im*s2.w ) + ( M[1][1].re*s1.x - M[1][1].im*s1.y ) + ( M[1][2].re*s1.z - M[1][2].im*s1.w );
(*(out+5)).y =  ( M[1][0].re*s2.w + M[1][0].im*s2.z ) + ( M[1][1].re*s1.y + M[1][1].im*s1.x ) + ( M[1][2].re*s1.w + M[1][2].im*s1.z );


(*(out+5)).z =  ( M[2][0].re*s2.z - M[2][0].im*s2.w ) + ( M[2][1].re*s1.x - M[2][1].im*s1.y ) + ( M[2][2].re*s1.z - M[2][2].im*s1.w );
(*(out+5)).w =  ( M[2][0].re*s2.w + M[2][0].im*s2.z ) + ( M[2][1].re*s1.y + M[2][1].im*s1.x ) + ( M[2][2].re*s1.w + M[2][2].im*s1.z );


}







//multipliziert su3-Matrix mal Spinor im Dirac-Raum
//we do not load the lower spinor components as we are working in the relativistic basis
//code in su3_MtV.txt -- generated with codegen
__device__ void dev_su3MtV_spintex_rel_up(dev_su3 M, int pos, dev_spinor * out){

dev_spinor s1, s2;

#ifndef HALF
 s1 = tex1Dfetch(spin_tex0,pos);
#else
 s1 = tex1Dfetch(spinhalf_tex0,pos);
 float norm = tex1Dfetch(spinnormhalf_tex,pos);
 s1.x *= norm; 
 s1.y *= norm; 
 s1.z *= norm;
 s1.w *= norm;
#endif

#ifndef HALF
 s2 = tex1Dfetch(spin_tex1,pos);
#else
 s2 = tex1Dfetch(spinhalf_tex1,pos);
 s2.x *= norm; 
 s2.y *= norm; 
 s2.z *= norm; 
 s2.w *= norm;
#endif


volatile float help;

help  = M[0][0].re*s1.x;
help -= M[0][0].im*s1.y;
help += M[0][1].re*s1.z;
help -= M[0][1].im*s1.w;
help += M[0][2].re*s2.x;
help -= M[0][2].im*s2.y;
(*(out+0)).x = help; 

help  = M[0][0].re*s1.y;
help += M[0][0].im*s1.x;
help += M[0][1].re*s1.w;
help += M[0][1].im*s1.z;
help += M[0][2].re*s2.y;
help += M[0][2].im*s2.x;
(*(out+0)).y = help;


help  = M[1][0].re*s1.x;
help -= M[1][0].im*s1.y;
help += M[1][1].re*s1.z;
help -= M[1][1].im*s1.w;
help += M[1][2].re*s2.x;
help -= M[1][2].im*s2.y;
(*(out+0)).z = help;

help  = M[1][0].re*s1.y;
help += M[1][0].im*s1.x;
help += M[1][1].re*s1.w;
help += M[1][1].im*s1.z;
help += M[1][2].re*s2.y;
help += M[1][2].im*s2.x;
(*(out+0)).w = help; 


help  = M[2][0].re*s1.x;
help -= M[2][0].im*s1.y;
help += M[2][1].re*s1.z;
help -= M[2][1].im*s1.w;
help += M[2][2].re*s2.x;
help -= M[2][2].im*s2.y;
(*(out+1)).x = help;

help  = M[2][0].re*s1.y;
help += M[2][0].im*s1.x;
help += M[2][1].re*s1.w;
help += M[2][1].im*s1.z;
help += M[2][2].re*s2.y;
help += M[2][2].im*s2.x;
(*(out+1)).y = help;


#ifndef HALF
 s1 = tex1Dfetch(spin_tex2,pos);
#else
 s1 = tex1Dfetch(spinhalf_tex2,pos);
 s1.x *= norm; 
 s1.y *= norm; 
 s1.z *= norm;
 s1.w *= norm;
#endif

help  = M[0][0].re*s2.z;
help -= M[0][0].im*s2.w;
help += M[0][1].re*s1.x;
help -= M[0][1].im*s1.y;
help += M[0][2].re*s1.z;
help -= M[0][2].im*s1.w;
(*(out+1)).z = help;

help  = M[0][0].re*s2.w;
help += M[0][0].im*s2.z;
help += M[0][1].re*s1.y;
help += M[0][1].im*s1.x;
help += M[0][2].re*s1.w;
help += M[0][2].im*s1.z;
(*(out+1)).w = help;


help  = M[1][0].re*s2.z;
help -= M[1][0].im*s2.w;
help += M[1][1].re*s1.x;
help -= M[1][1].im*s1.y;
help += M[1][2].re*s1.z;
help -= M[1][2].im*s1.w;
(*(out+2)).x = help; 

help  = M[1][0].re*s2.w;
help += M[1][0].im*s2.z;
help += M[1][1].re*s1.y;
help += M[1][1].im*s1.x;
help += M[1][2].re*s1.w;
help += M[1][2].im*s1.z;
(*(out+2)).y = help; 


help  = M[2][0].re*s2.z;
help -= M[2][0].im*s2.w;
help += M[2][1].re*s1.x;
help -= M[2][1].im*s1.y;
help += M[2][2].re*s1.z;
help -= M[2][2].im*s1.w;
(*(out+2)).z = help; 

help  = M[2][0].re*s2.w;
help += M[2][0].im*s2.z;
help += M[2][1].re*s1.y;
help += M[2][1].im*s1.x;
help += M[2][2].re*s1.w;
help += M[2][2].im*s1.z;
(*(out+2)).w = help; 




(*(out+3)).x = 0.0f;
(*(out+3)).y = 0.0f;  
(*(out+3)).z = 0.0f; 
(*(out+3)).w = 0.0f; 


(*(out+4)).x = 0.0f; 
(*(out+4)).y = 0.0f; 
(*(out+4)).z = 0.0f; 
(*(out+4)).w = 0.0f; 


(*(out+5)).x = 0.0f;
(*(out+5)).y = 0.0f;
(*(out+5)).z = 0.0f; 
(*(out+5)).w = 0.0f; 


}



__device__ void dev_su3MtV_spintex2_rel_up(dev_su3 M, int pos, dev_spinor * out){

dev_spinor s1, s2;


 s1 = tex1Dfetch(spin_tex_dn0,pos);
 s2 = tex1Dfetch(spin_tex_dn1,pos);


volatile float help;

help  = M[0][0].re*s1.x;
help -= M[0][0].im*s1.y;
help += M[0][1].re*s1.z;
help -= M[0][1].im*s1.w;
help += M[0][2].re*s2.x;
help -= M[0][2].im*s2.y;
(*(out+0)).x = help; 

help  = M[0][0].re*s1.y;
help += M[0][0].im*s1.x;
help += M[0][1].re*s1.w;
help += M[0][1].im*s1.z;
help += M[0][2].re*s2.y;
help += M[0][2].im*s2.x;
(*(out+0)).y = help;


help  = M[1][0].re*s1.x;
help -= M[1][0].im*s1.y;
help += M[1][1].re*s1.z;
help -= M[1][1].im*s1.w;
help += M[1][2].re*s2.x;
help -= M[1][2].im*s2.y;
(*(out+0)).z = help;

help  = M[1][0].re*s1.y;
help += M[1][0].im*s1.x;
help += M[1][1].re*s1.w;
help += M[1][1].im*s1.z;
help += M[1][2].re*s2.y;
help += M[1][2].im*s2.x;
(*(out+0)).w = help; 


help  = M[2][0].re*s1.x;
help -= M[2][0].im*s1.y;
help += M[2][1].re*s1.z;
help -= M[2][1].im*s1.w;
help += M[2][2].re*s2.x;
help -= M[2][2].im*s2.y;
(*(out+1)).x = help;

help  = M[2][0].re*s1.y;
help += M[2][0].im*s1.x;
help += M[2][1].re*s1.w;
help += M[2][1].im*s1.z;
help += M[2][2].re*s2.y;
help += M[2][2].im*s2.x;
(*(out+1)).y = help;



 s1 = tex1Dfetch(spin_tex_dn2,pos);


help  = M[0][0].re*s2.z;
help -= M[0][0].im*s2.w;
help += M[0][1].re*s1.x;
help -= M[0][1].im*s1.y;
help += M[0][2].re*s1.z;
help -= M[0][2].im*s1.w;
(*(out+1)).z = help;

help  = M[0][0].re*s2.w;
help += M[0][0].im*s2.z;
help += M[0][1].re*s1.y;
help += M[0][1].im*s1.x;
help += M[0][2].re*s1.w;
help += M[0][2].im*s1.z;
(*(out+1)).w = help;


help  = M[1][0].re*s2.z;
help -= M[1][0].im*s2.w;
help += M[1][1].re*s1.x;
help -= M[1][1].im*s1.y;
help += M[1][2].re*s1.z;
help -= M[1][2].im*s1.w;
(*(out+2)).x = help; 

help  = M[1][0].re*s2.w;
help += M[1][0].im*s2.z;
help += M[1][1].re*s1.y;
help += M[1][1].im*s1.x;
help += M[1][2].re*s1.w;
help += M[1][2].im*s1.z;
(*(out+2)).y = help; 


help  = M[2][0].re*s2.z;
help -= M[2][0].im*s2.w;
help += M[2][1].re*s1.x;
help -= M[2][1].im*s1.y;
help += M[2][2].re*s1.z;
help -= M[2][2].im*s1.w;
(*(out+2)).z = help; 

help  = M[2][0].re*s2.w;
help += M[2][0].im*s2.z;
help += M[2][1].re*s1.y;
help += M[2][1].im*s1.x;
help += M[2][2].re*s1.w;
help += M[2][2].im*s1.z;
(*(out+2)).w = help; 




(*(out+3)).x = 0.0f;
(*(out+3)).y = 0.0f;  
(*(out+3)).z = 0.0f; 
(*(out+3)).w = 0.0f; 


(*(out+4)).x = 0.0f; 
(*(out+4)).y = 0.0f; 
(*(out+4)).z = 0.0f; 
(*(out+4)).w = 0.0f; 


(*(out+5)).x = 0.0f;
(*(out+5)).y = 0.0f;
(*(out+5)).z = 0.0f; 
(*(out+5)).w = 0.0f; 


}




//multipliziert su3-Matrix mal Spinor im Dirac-Raum
//we do not load the upper spinor components as we are working in the relativistic basis
//code in su3_MtV.txt -- generated with codegen
__device__ void dev_su3MtV_spintex_rel_down(dev_su3 M, int pos, dev_spinor * out){

dev_spinor s1, s2;


(*(out+0)).x = 0.0f;
(*(out+0)).y = 0.0f;
(*(out+0)).z = 0.0f;
(*(out+0)).w = 0.0f;


(*(out+1)).x = 0.0f;
(*(out+1)).y = 0.0f;
(*(out+1)).z = 0.0f;
(*(out+1)).w = 0.0f;


(*(out+2)).x = 0.0f;
(*(out+2)).y = 0.0f;
(*(out+2)).z = 0.0f;
(*(out+2)).w = 0.0f;


#ifndef HALF
 s1 = tex1Dfetch(spin_tex3,pos);
#else
 float norm = tex1Dfetch(spinnormhalf_tex,pos);
 s1 = tex1Dfetch(spinhalf_tex3,pos);
 s1.x *= norm; 
 s1.y *= norm; 
 s1.z *= norm;
 s1.w *= norm;
#endif

#ifndef HALF
 s2 = tex1Dfetch(spin_tex4,pos);
#else
 s2 = tex1Dfetch(spinhalf_tex4,pos);
 s2.x *= norm; 
 s2.y *= norm; 
 s2.z *= norm; 
 s2.w *= norm;
#endif

volatile float help;

help  = M[0][0].re*s1.x;
help -= M[0][0].im*s1.y;
help += M[0][1].re*s1.z;
help -= M[0][1].im*s1.w;
help += M[0][2].re*s2.x;
help -= M[0][2].im*s2.y;
(*(out+3)).x =  help; 

help  = M[0][0].re*s1.y;
help += M[0][0].im*s1.x;
help += M[0][1].re*s1.w;
help += M[0][1].im*s1.z;
help += M[0][2].re*s2.y;
help += M[0][2].im*s2.x;
(*(out+3)).y = help; 


help  = M[1][0].re*s1.x;
help -= M[1][0].im*s1.y;
help += M[1][1].re*s1.z;
help -= M[1][1].im*s1.w;
help += M[1][2].re*s2.x;
help -= M[1][2].im*s2.y;
(*(out+3)).z = help;

help  = M[1][0].re*s1.y;
help += M[1][0].im*s1.x;
help += M[1][1].re*s1.w;
help += M[1][1].im*s1.z;
help += M[1][2].re*s2.y;
help += M[1][2].im*s2.x;
(*(out+3)).w = help;


help  = M[2][0].re*s1.x;
help -= M[2][0].im*s1.y;
help += M[2][1].re*s1.z;
help -= M[2][1].im*s1.w;
help += M[2][2].re*s2.x;
help -= M[2][2].im*s2.y;
(*(out+4)).x = help; 

help  = M[2][0].re*s1.y;
help += M[2][0].im*s1.x;
help += M[2][1].re*s1.w;
help += M[2][1].im*s1.z;
help += M[2][2].re*s2.y;
help += M[2][2].im*s2.x;
(*(out+4)).y = help; 


#ifndef HALF
 s1 = tex1Dfetch(spin_tex5,pos);
#else
 s1 = tex1Dfetch(spinhalf_tex5,pos);
 s1.x *= norm; 
 s1.y *= norm; 
 s1.z *= norm;
 s1.w *= norm;
#endif

help  = M[0][0].re*s2.z;
help -= M[0][0].im*s2.w;
help += M[0][1].re*s1.x;
help -= M[0][1].im*s1.y;
help += M[0][2].re*s1.z;
help -= M[0][2].im*s1.w;
(*(out+4)).z = help; 

help  = M[0][0].re*s2.w;
help += M[0][0].im*s2.z;
help += M[0][1].re*s1.y;
help += M[0][1].im*s1.x;
help += M[0][2].re*s1.w;
help += M[0][2].im*s1.z;
(*(out+4)).w = help;  


help  = M[1][0].re*s2.z;
help -= M[1][0].im*s2.w;
help += M[1][1].re*s1.x;
help -= M[1][1].im*s1.y;
help += M[1][2].re*s1.z;
help -= M[1][2].im*s1.w;
(*(out+5)).x = help;

help  = M[1][0].re*s2.w;
help += M[1][0].im*s2.z;
help += M[1][1].re*s1.y;
help += M[1][1].im*s1.x;
help += M[1][2].re*s1.w;
help += M[1][2].im*s1.z;
(*(out+5)).y = help; 


help  = M[2][0].re*s2.z;
help -= M[2][0].im*s2.w;
help += M[2][1].re*s1.x;
help -= M[2][1].im*s1.y;
help += M[2][2].re*s1.z;
help -= M[2][2].im*s1.w;
(*(out+5)).z = help; 

help  = M[2][0].re*s2.w;
help += M[2][0].im*s2.z;
help += M[2][1].re*s1.y;
help += M[2][1].im*s1.x;
help += M[2][2].re*s1.w;
help += M[2][2].im*s1.z;
(*(out+5)).w = help; 


}




__device__ void dev_su3MtV_spintex2_rel_down(dev_su3 M, int pos, dev_spinor * out){

dev_spinor s1, s2;


(*(out+0)).x = 0.0f;
(*(out+0)).y = 0.0f;
(*(out+0)).z = 0.0f;
(*(out+0)).w = 0.0f;


(*(out+1)).x = 0.0f;
(*(out+1)).y = 0.0f;
(*(out+1)).z = 0.0f;
(*(out+1)).w = 0.0f;


(*(out+2)).x = 0.0f;
(*(out+2)).y = 0.0f;
(*(out+2)).z = 0.0f;
(*(out+2)).w = 0.0f;



 s1 = tex1Dfetch(spin_tex_dn3,pos);
 s2 = tex1Dfetch(spin_tex_dn4,pos);


volatile float help;

help  = M[0][0].re*s1.x;
help -= M[0][0].im*s1.y;
help += M[0][1].re*s1.z;
help -= M[0][1].im*s1.w;
help += M[0][2].re*s2.x;
help -= M[0][2].im*s2.y;
(*(out+3)).x =  help; 

help  = M[0][0].re*s1.y;
help += M[0][0].im*s1.x;
help += M[0][1].re*s1.w;
help += M[0][1].im*s1.z;
help += M[0][2].re*s2.y;
help += M[0][2].im*s2.x;
(*(out+3)).y = help; 


help  = M[1][0].re*s1.x;
help -= M[1][0].im*s1.y;
help += M[1][1].re*s1.z;
help -= M[1][1].im*s1.w;
help += M[1][2].re*s2.x;
help -= M[1][2].im*s2.y;
(*(out+3)).z = help;

help  = M[1][0].re*s1.y;
help += M[1][0].im*s1.x;
help += M[1][1].re*s1.w;
help += M[1][1].im*s1.z;
help += M[1][2].re*s2.y;
help += M[1][2].im*s2.x;
(*(out+3)).w = help;


help  = M[2][0].re*s1.x;
help -= M[2][0].im*s1.y;
help += M[2][1].re*s1.z;
help -= M[2][1].im*s1.w;
help += M[2][2].re*s2.x;
help -= M[2][2].im*s2.y;
(*(out+4)).x = help; 

help  = M[2][0].re*s1.y;
help += M[2][0].im*s1.x;
help += M[2][1].re*s1.w;
help += M[2][1].im*s1.z;
help += M[2][2].re*s2.y;
help += M[2][2].im*s2.x;
(*(out+4)).y = help; 



 s1 = tex1Dfetch(spin_tex_dn5,pos);


help  = M[0][0].re*s2.z;
help -= M[0][0].im*s2.w;
help += M[0][1].re*s1.x;
help -= M[0][1].im*s1.y;
help += M[0][2].re*s1.z;
help -= M[0][2].im*s1.w;
(*(out+4)).z = help; 

help  = M[0][0].re*s2.w;
help += M[0][0].im*s2.z;
help += M[0][1].re*s1.y;
help += M[0][1].im*s1.x;
help += M[0][2].re*s1.w;
help += M[0][2].im*s1.z;
(*(out+4)).w = help;  


help  = M[1][0].re*s2.z;
help -= M[1][0].im*s2.w;
help += M[1][1].re*s1.x;
help -= M[1][1].im*s1.y;
help += M[1][2].re*s1.z;
help -= M[1][2].im*s1.w;
(*(out+5)).x = help;

help  = M[1][0].re*s2.w;
help += M[1][0].im*s2.z;
help += M[1][1].re*s1.y;
help += M[1][1].im*s1.x;
help += M[1][2].re*s1.w;
help += M[1][2].im*s1.z;
(*(out+5)).y = help; 


help  = M[2][0].re*s2.z;
help -= M[2][0].im*s2.w;
help += M[2][1].re*s1.x;
help -= M[2][1].im*s1.y;
help += M[2][2].re*s1.z;
help -= M[2][2].im*s1.w;
(*(out+5)).z = help; 

help  = M[2][0].re*s2.w;
help += M[2][0].im*s2.z;
help += M[2][1].re*s1.y;
help += M[2][1].im*s1.x;
help += M[2][2].re*s1.w;
help += M[2][2].im*s1.z;
(*(out+5)).w = help; 


}





//multipliziert su3-Matrix mal Spinor im Dirac-Raum
//code in su3_MtV.txt -- generated with codegen
__device__ void dev_su3MtV(dev_su3 M, const dev_spinor * s, dev_spinor * out){

(*(out+0)).x =  ( M[0][0].re*(*(s+0*DEVOFF)).x - M[0][0].im*(*(s+0*DEVOFF)).y ) + ( M[0][1].re*(*(s+0*DEVOFF)).z - M[0][1].im*(*(s+0*DEVOFF)).w ) + ( M[0][2].re*(*(s+1*DEVOFF)).x - M[0][2].im*(*(s+1*DEVOFF)).y );
(*(out+0)).y = ( M[0][0].re*(*(s+0*DEVOFF)).y + M[0][0].im*(*(s+0*DEVOFF)).x ) + ( M[0][1].re*(*(s+0*DEVOFF)).w + M[0][1].im*(*(s+0*DEVOFF)).z ) + ( M[0][2].re*(*(s+1*DEVOFF)).y + M[0][2].im*(*(s+1*DEVOFF)).x );


(*(out+0)).z =  ( M[1][0].re*(*(s+0*DEVOFF)).x - M[1][0].im*(*(s+0*DEVOFF)).y ) + ( M[1][1].re*(*(s+0*DEVOFF)).z - M[1][1].im*(*(s+0*DEVOFF)).w ) + ( M[1][2].re*(*(s+1*DEVOFF)).x - M[1][2].im*(*(s+1*DEVOFF)).y );
(*(out+0)).w =  ( M[1][0].re*(*(s+0*DEVOFF)).y + M[1][0].im*(*(s+0*DEVOFF)).x ) + ( M[1][1].re*(*(s+0*DEVOFF)).w + M[1][1].im*(*(s+0*DEVOFF)).z ) + ( M[1][2].re*(*(s+1*DEVOFF)).y + M[1][2].im*(*(s+1*DEVOFF)).x );


(*(out+1)).x = ( M[2][0].re*(*(s+0*DEVOFF)).x - M[2][0].im*(*(s+0*DEVOFF)).y ) + ( M[2][1].re*(*(s+0*DEVOFF)).z - M[2][1].im*(*(s+0*DEVOFF)).w ) + ( M[2][2].re*(*(s+1*DEVOFF)).x - M[2][2].im*(*(s+1*DEVOFF)).y );
(*(out+1)).y =  ( M[2][0].re*(*(s+0*DEVOFF)).y + M[2][0].im*(*(s+0*DEVOFF)).x ) + ( M[2][1].re*(*(s+0*DEVOFF)).w + M[2][1].im*(*(s+0*DEVOFF)).z ) + ( M[2][2].re*(*(s+1*DEVOFF)).y + M[2][2].im*(*(s+1*DEVOFF)).x );


(*(out+1)).z =  ( M[0][0].re*(*(s+1*DEVOFF)).z - M[0][0].im*(*(s+1*DEVOFF)).w ) + ( M[0][1].re*(*(s+2*DEVOFF)).x - M[0][1].im*(*(s+2*DEVOFF)).y ) + ( M[0][2].re*(*(s+2*DEVOFF)).z - M[0][2].im*(*(s+2*DEVOFF)).w );
(*(out+1)).w =  ( M[0][0].re*(*(s+1*DEVOFF)).w + M[0][0].im*(*(s+1*DEVOFF)).z ) + ( M[0][1].re*(*(s+2*DEVOFF)).y + M[0][1].im*(*(s+2*DEVOFF)).x ) + ( M[0][2].re*(*(s+2*DEVOFF)).w + M[0][2].im*(*(s+2*DEVOFF)).z );


(*(out+2)).x = ( M[1][0].re*(*(s+1*DEVOFF)).z - M[1][0].im*(*(s+1*DEVOFF)).w ) + ( M[1][1].re*(*(s+2*DEVOFF)).x - M[1][1].im*(*(s+2*DEVOFF)).y ) + ( M[1][2].re*(*(s+2*DEVOFF)).z - M[1][2].im*(*(s+2*DEVOFF)).w );
(*(out+2)).y =  ( M[1][0].re*(*(s+1*DEVOFF)).w + M[1][0].im*(*(s+1*DEVOFF)).z ) + ( M[1][1].re*(*(s+2*DEVOFF)).y + M[1][1].im*(*(s+2*DEVOFF)).x ) + ( M[1][2].re*(*(s+2*DEVOFF)).w + M[1][2].im*(*(s+2*DEVOFF)).z );


(*(out+2)).z =  ( M[2][0].re*(*(s+1*DEVOFF)).z - M[2][0].im*(*(s+1*DEVOFF)).w ) + ( M[2][1].re*(*(s+2*DEVOFF)).x - M[2][1].im*(*(s+2*DEVOFF)).y ) + ( M[2][2].re*(*(s+2*DEVOFF)).z - M[2][2].im*(*(s+2*DEVOFF)).w );
(*(out+2)).w =  ( M[2][0].re*(*(s+1*DEVOFF)).w + M[2][0].im*(*(s+1*DEVOFF)).z ) + ( M[2][1].re*(*(s+2*DEVOFF)).y + M[2][1].im*(*(s+2*DEVOFF)).x ) + ( M[2][2].re*(*(s+2*DEVOFF)).w + M[2][2].im*(*(s+2*DEVOFF)).z );


(*(out+3)).x =  ( M[0][0].re*(*(s+3*DEVOFF)).x - M[0][0].im*(*(s+3*DEVOFF)).y ) + ( M[0][1].re*(*(s+3*DEVOFF)).z - M[0][1].im*(*(s+3*DEVOFF)).w ) + ( M[0][2].re*(*(s+4*DEVOFF)).x - M[0][2].im*(*(s+4*DEVOFF)).y );
(*(out+3)).y =   ( M[0][0].re*(*(s+3*DEVOFF)).y + M[0][0].im*(*(s+3*DEVOFF)).x ) + ( M[0][1].re*(*(s+3*DEVOFF)).w + M[0][1].im*(*(s+3*DEVOFF)).z ) + ( M[0][2].re*(*(s+4*DEVOFF)).y + M[0][2].im*(*(s+4*DEVOFF)).x );


(*(out+3)).z =  ( M[1][0].re*(*(s+3*DEVOFF)).x - M[1][0].im*(*(s+3*DEVOFF)).y ) + ( M[1][1].re*(*(s+3*DEVOFF)).z - M[1][1].im*(*(s+3*DEVOFF)).w ) + ( M[1][2].re*(*(s+4*DEVOFF)).x - M[1][2].im*(*(s+4*DEVOFF)).y );
(*(out+3)).w =  ( M[1][0].re*(*(s+3*DEVOFF)).y + M[1][0].im*(*(s+3*DEVOFF)).x ) + ( M[1][1].re*(*(s+3*DEVOFF)).w + M[1][1].im*(*(s+3*DEVOFF)).z ) + ( M[1][2].re*(*(s+4*DEVOFF)).y + M[1][2].im*(*(s+4*DEVOFF)).x );


(*(out+4)).x =  ( M[2][0].re*(*(s+3*DEVOFF)).x - M[2][0].im*(*(s+3*DEVOFF)).y ) + ( M[2][1].re*(*(s+3*DEVOFF)).z - M[2][1].im*(*(s+3*DEVOFF)).w ) + ( M[2][2].re*(*(s+4*DEVOFF)).x - M[2][2].im*(*(s+4*DEVOFF)).y );
(*(out+4)).y =  ( M[2][0].re*(*(s+3*DEVOFF)).y + M[2][0].im*(*(s+3*DEVOFF)).x ) + ( M[2][1].re*(*(s+3*DEVOFF)).w + M[2][1].im*(*(s+3*DEVOFF)).z ) + ( M[2][2].re*(*(s+4*DEVOFF)).y + M[2][2].im*(*(s+4*DEVOFF)).x );


(*(out+4)).z =  ( M[0][0].re*(*(s+4*DEVOFF)).z - M[0][0].im*(*(s+4*DEVOFF)).w ) + ( M[0][1].re*(*(s+5*DEVOFF)).x - M[0][1].im*(*(s+5*DEVOFF)).y ) + ( M[0][2].re*(*(s+5*DEVOFF)).z - M[0][2].im*(*(s+5*DEVOFF)).w );
(*(out+4)).w =   ( M[0][0].re*(*(s+4*DEVOFF)).w + M[0][0].im*(*(s+4*DEVOFF)).z ) + ( M[0][1].re*(*(s+5*DEVOFF)).y + M[0][1].im*(*(s+5*DEVOFF)).x ) + ( M[0][2].re*(*(s+5*DEVOFF)).w + M[0][2].im*(*(s+5*DEVOFF)).z );


(*(out+5)).x = ( M[1][0].re*(*(s+4*DEVOFF)).z - M[1][0].im*(*(s+4*DEVOFF)).w ) + ( M[1][1].re*(*(s+5*DEVOFF)).x - M[1][1].im*(*(s+5*DEVOFF)).y ) + ( M[1][2].re*(*(s+5*DEVOFF)).z - M[1][2].im*(*(s+5*DEVOFF)).w );
(*(out+5)).y =  ( M[1][0].re*(*(s+4*DEVOFF)).w + M[1][0].im*(*(s+4*DEVOFF)).z ) + ( M[1][1].re*(*(s+5*DEVOFF)).y + M[1][1].im*(*(s+5*DEVOFF)).x ) + ( M[1][2].re*(*(s+5*DEVOFF)).w + M[1][2].im*(*(s+5*DEVOFF)).z );


(*(out+5)).z =  ( M[2][0].re*(*(s+4*DEVOFF)).z - M[2][0].im*(*(s+4*DEVOFF)).w ) + ( M[2][1].re*(*(s+5*DEVOFF)).x - M[2][1].im*(*(s+5*DEVOFF)).y ) + ( M[2][2].re*(*(s+5*DEVOFF)).z - M[2][2].im*(*(s+5*DEVOFF)).w );
(*(out+5)).w =  ( M[2][0].re*(*(s+4*DEVOFF)).w + M[2][0].im*(*(s+4*DEVOFF)).z ) + ( M[2][1].re*(*(s+5*DEVOFF)).y + M[2][1].im*(*(s+5*DEVOFF)).x ) + ( M[2][2].re*(*(s+5*DEVOFF)).w + M[2][2].im*(*(s+5*DEVOFF)).z );
}







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
__device__ void dev_su3MtV_kappaP3_plus(dev_su3 M, const dev_spinor * s, dev_spinor * out, dev_complex kappa){

  
  dev_spinor sh0,sh1,sh2,zh0,zh1,zh2;
  
  
  //first apply Projector on upper halfspinor 
     sh0.x = kappa.re*( (*(s+0*DEVOFF)).x - (*(s+3*DEVOFF)).y);
     sh0.x -= kappa.im*( (*(s+0*DEVOFF)).y + (*(s+3*DEVOFF)).x);
     sh0.y = kappa.re*( (*(s+0*DEVOFF)).y + (*(s+3*DEVOFF)).x);
     sh0.y += kappa.im*( (*(s+0*DEVOFF)).x - (*(s+3*DEVOFF)).y);
   
     sh0.z = kappa.re*( (*(s+0*DEVOFF)).z - (*(s+3*DEVOFF)).w);
     sh0.z -= kappa.im*( (*(s+0*DEVOFF)).w + (*(s+3*DEVOFF)).z);
     sh0.w = kappa.re*( (*(s+0*DEVOFF)).w + (*(s+3*DEVOFF)).z);
     sh0.w += kappa.im*( (*(s+0*DEVOFF)).z - (*(s+3*DEVOFF)).w);
     
     sh1.x = kappa.re*( (*(s+1*DEVOFF)).x - (*(s+4*DEVOFF)).y);
     sh1.x -= kappa.im*( (*(s+1*DEVOFF)).y + (*(s+4*DEVOFF)).x);
     sh1.y = kappa.re*( (*(s+1*DEVOFF)).y + (*(s+4*DEVOFF)).x);
     sh1.y += kappa.im*( (*(s+1*DEVOFF)).x - (*(s+4*DEVOFF)).y);  


  

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


     sh1.z = kappa.re*( (*(s+1*DEVOFF)).z + (*(s+4*DEVOFF)).w);
     sh1.z -= kappa.im*( (*(s+1*DEVOFF)).w - (*(s+4*DEVOFF)).z);
     sh1.w = kappa.re*( (*(s+1*DEVOFF)).w - (*(s+4*DEVOFF)).z);
     sh1.w += kappa.im*( (*(s+1*DEVOFF)).z + (*(s+4*DEVOFF)).w);
         
     sh2.x = kappa.re*( (*(s+2*DEVOFF)).x + (*(s+5*DEVOFF)).y);
     sh2.x -= kappa.im*( (*(s+2*DEVOFF)).y - (*(s+5*DEVOFF)).x);
     sh2.y = kappa.re*( (*(s+2*DEVOFF)).y - (*(s+5*DEVOFF)).x);
     sh2.y += kappa.im*( (*(s+2*DEVOFF)).x + (*(s+5*DEVOFF)).y);
       
     sh2.z = kappa.re*( (*(s+2*DEVOFF)).z + (*(s+5*DEVOFF)).w);
     sh2.z -= kappa.im*( (*(s+2*DEVOFF)).w - (*(s+5*DEVOFF)).z);
     sh2.w = kappa.re*( (*(s+2*DEVOFF)).w - (*(s+5*DEVOFF)).z);
     sh2.w += kappa.im*( (*(s+2*DEVOFF)).z + (*(s+5*DEVOFF)).w);

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

__device__ void dev_su3MtV_kappaP3_plus_spintex(dev_su3 M, int pos, dev_spinor * out, dev_complex kappa){

  
  dev_spinor sh0,sh1,sh2,zh0,zh1,zh2;
  
  dev_spinor s0, s1, s2, s3, s4, s5;
  s0 = tex1Dfetch(spin_tex0,pos);
  s1 = tex1Dfetch(spin_tex1,pos);
  s2 = tex1Dfetch(spin_tex2,pos);
  s3 = tex1Dfetch(spin_tex3,pos);
  s4 = tex1Dfetch(spin_tex4,pos); 
  s5 = tex1Dfetch(spin_tex5,pos);  
  //first apply Projector on upper halfspinor 
     sh0.x = kappa.re*(s0.x - s3.y);
     sh0.x -= kappa.im*(s0.y + s3.x);
     sh0.y = kappa.re*(s0.y + s3.x);
     sh0.y += kappa.im*(s0.x - s3.y);
     
     sh0.z = kappa.re*(s0.z - s3.w);
     sh0.z -= kappa.im*(s0.w + s3.z);
     sh0.w = kappa.re*(s0.w + s3.z);
     sh0.w += kappa.im*(s0.z - s3.w);
   
     sh1.x = kappa.re*(s1.x - s4.y);
     sh1.x -= kappa.im*(s1.y + s4.x);
     sh1.y = kappa.re*(s1.y + s4.x);
     sh1.y += kappa.im*(s1.x - s4.y);   


  

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


     sh1.z = kappa.re*(s1.z + s4.w);
     sh1.z -= kappa.im*(s1.w - s4.z);
     sh1.w = kappa.re*(s1.w - s4.z);
     sh1.w += kappa.im*(s1.z + s4.w);
         
     sh2.x = kappa.re*(s2.x + s5.y);
     sh2.x -= kappa.im*(s2.y - s5.x);
     sh2.y = kappa.re*(s2.y - s5.x);
     sh2.y += kappa.im*(s2.x + s5.y);
       
     sh2.z = kappa.re*(s2.z + s5.w);
     sh2.z -= kappa.im*(s2.w - s5.z);
     sh2.w = kappa.re*(s2.w - s5.z);
     sh2.w += kappa.im*(s2.z + s5.w);

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



__device__ void dev_su3MtV_kappaP3_plus_spintex2(dev_su3 M, int pos, dev_spinor * out, dev_complex kappa){

  
  dev_spinor sh0,sh1,sh2,zh0,zh1,zh2;
  
  dev_spinor s0, s1, s2, s3, s4, s5;
  s0 = tex1Dfetch(spin_tex_dn0,pos);
  s1 = tex1Dfetch(spin_tex_dn1,pos);
  s2 = tex1Dfetch(spin_tex_dn2,pos);
  s3 = tex1Dfetch(spin_tex_dn3,pos);
  s4 = tex1Dfetch(spin_tex_dn4,pos); 
  s5 = tex1Dfetch(spin_tex_dn5,pos);  
  //first apply Projector on upper halfspinor 
     sh0.x = kappa.re*(s0.x - s3.y);
     sh0.x -= kappa.im*(s0.y + s3.x);
     sh0.y = kappa.re*(s0.y + s3.x);
     sh0.y += kappa.im*(s0.x - s3.y);
     
     sh0.z = kappa.re*(s0.z - s3.w);
     sh0.z -= kappa.im*(s0.w + s3.z);
     sh0.w = kappa.re*(s0.w + s3.z);
     sh0.w += kappa.im*(s0.z - s3.w);
   
     sh1.x = kappa.re*(s1.x - s4.y);
     sh1.x -= kappa.im*(s1.y + s4.x);
     sh1.y = kappa.re*(s1.y + s4.x);
     sh1.y += kappa.im*(s1.x - s4.y);   


  

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


     sh1.z = kappa.re*(s1.z + s4.w);
     sh1.z -= kappa.im*(s1.w - s4.z);
     sh1.w = kappa.re*(s1.w - s4.z);
     sh1.w += kappa.im*(s1.z + s4.w);
         
     sh2.x = kappa.re*(s2.x + s5.y);
     sh2.x -= kappa.im*(s2.y - s5.x);
     sh2.y = kappa.re*(s2.y - s5.x);
     sh2.y += kappa.im*(s2.x + s5.y);
       
     sh2.z = kappa.re*(s2.z + s5.w);
     sh2.z -= kappa.im*(s2.w - s5.z);
     sh2.w = kappa.re*(s2.w - s5.z);
     sh2.w += kappa.im*(s2.z + s5.w);

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





__device__ void dev_su3MtV_kappaP3_plus_spintex_ud(dev_su3 M, int pos, dev_spinor * out_u, dev_spinor * out_d , dev_complex kappa){

  
  dev_spinor sh0_u,sh1_u,sh2_u,zh0_u,zh1_u,zh2_u;
  dev_spinor sh0_d,sh1_d,sh2_d,zh0_d,zh1_d,zh2_d;
  
  dev_spinor s0_u, s1_u, s2_u, s3_u, s4_u, s5_u;
  dev_spinor s0_d, s1_d, s2_d, s3_d, s4_d, s5_d;  

  s0_u = tex1Dfetch(spin_tex0,pos);
  s1_u = tex1Dfetch(spin_tex1,pos);
  s2_u = tex1Dfetch(spin_tex2,pos);
  s3_u = tex1Dfetch(spin_tex3,pos);
  s4_u = tex1Dfetch(spin_tex4,pos); 
  s5_u = tex1Dfetch(spin_tex5,pos);    
  
  s0_d = tex1Dfetch(spin_tex_dn0,pos);
  s1_d = tex1Dfetch(spin_tex_dn1,pos);
  s2_d = tex1Dfetch(spin_tex_dn2,pos);
  s3_d = tex1Dfetch(spin_tex_dn3,pos);
  s4_d = tex1Dfetch(spin_tex_dn4,pos); 
  s5_d = tex1Dfetch(spin_tex_dn5,pos);  
  
  //first apply Projector on upper halfspinor 
     sh0_u.x = kappa.re*(s0_u.x - s3_u.y);
     sh0_u.x -= kappa.im*(s0_u.y + s3_u.x);
     sh0_u.y = kappa.re*(s0_u.y + s3_u.x);
     sh0_u.y += kappa.im*(s0_u.x - s3_u.y);
     
     sh0_u.z = kappa.re*(s0_u.z - s3_u.w);
     sh0_u.z -= kappa.im*(s0_u.w + s3_u.z);
     sh0_u.w = kappa.re*(s0_u.w + s3_u.z);
     sh0_u.w += kappa.im*(s0_u.z - s3_u.w);
   
     sh1_u.x = kappa.re*(s1_u.x - s4_u.y);
     sh1_u.x -= kappa.im*(s1_u.y + s4_u.x);
     sh1_u.y = kappa.re*(s1_u.y + s4_u.x);
     sh1_u.y += kappa.im*(s1_u.x - s4_u.y);   

     sh0_d.x = kappa.re*(s0_d.x - s3_d.y);
     sh0_d.x -= kappa.im*(s0_d.y + s3_d.x);
     sh0_d.y = kappa.re*(s0_d.y + s3_d.x);
     sh0_d.y += kappa.im*(s0_d.x - s3_d.y);
     
     sh0_d.z = kappa.re*(s0_d.z - s3_d.w);
     sh0_d.z -= kappa.im*(s0_d.w + s3_d.z);
     sh0_d.w = kappa.re*(s0_d.w + s3_d.z);
     sh0_d.w += kappa.im*(s0_d.z - s3_d.w);
   
     sh1_d.x = kappa.re*(s1_d.x - s4_d.y);
     sh1_d.x -= kappa.im*(s1_d.y + s4_d.x);
     sh1_d.y = kappa.re*(s1_d.y + s4_d.x);
     sh1_d.y += kappa.im*(s1_d.x - s4_d.y);   

  

//Multiply by gauge field 
zh0_u.x  = M[0][0].re*sh0_u.x;
zh0_u.x -= M[0][0].im*sh0_u.y;
zh0_u.x += M[0][1].re*sh0_u.z;
zh0_u.x -= M[0][1].im*sh0_u.w;
zh0_u.x += M[0][2].re*sh1_u.x;
zh0_u.x -= M[0][2].im*sh1_u.y;


zh0_u.y  = M[0][0].re*sh0_u.y;
zh0_u.y += M[0][0].im*sh0_u.x;
zh0_u.y += M[0][1].re*sh0_u.w;
zh0_u.y += M[0][1].im*sh0_u.z;
zh0_u.y += M[0][2].re*sh1_u.y;
zh0_u.y += M[0][2].im*sh1_u.x;
(*(out_u+0)).x -= zh0_u.x;
(*(out_u+0)).y -= zh0_u.y;


zh0_d.x  = M[0][0].re*sh0_d.x;
zh0_d.x -= M[0][0].im*sh0_d.y;
zh0_d.x += M[0][1].re*sh0_d.z;
zh0_d.x -= M[0][1].im*sh0_d.w;
zh0_d.x += M[0][2].re*sh1_d.x;
zh0_d.x -= M[0][2].im*sh1_d.y;


zh0_d.y  = M[0][0].re*sh0_d.y;
zh0_d.y += M[0][0].im*sh0_d.x;
zh0_d.y += M[0][1].re*sh0_d.w;
zh0_d.y += M[0][1].im*sh0_d.z;
zh0_d.y += M[0][2].re*sh1_d.y;
zh0_d.y += M[0][2].im*sh1_d.x;
(*(out_d+0)).x -= zh0_d.x;
(*(out_d+0)).y -= zh0_d.y;


zh0_u.z  = M[1][0].re*sh0_u.x;
zh0_u.z -= M[1][0].im*sh0_u.y;
zh0_u.z += M[1][1].re*sh0_u.z;
zh0_u.z -= M[1][1].im*sh0_u.w;
zh0_u.z += M[1][2].re*sh1_u.x;
zh0_u.z -= M[1][2].im*sh1_u.y;

zh0_u.w  = M[1][0].re*sh0_u.y;
zh0_u.w += M[1][0].im*sh0_u.x;
zh0_u.w += M[1][1].re*sh0_u.w;
zh0_u.w += M[1][1].im*sh0_u.z;
zh0_u.w += M[1][2].re*sh1_u.y;
zh0_u.w += M[1][2].im*sh1_u.x;
(*(out_u+0)).z -= zh0_u.z;
(*(out_u+0)).w -= zh0_u.w;


zh0_d.z  = M[1][0].re*sh0_d.x;
zh0_d.z -= M[1][0].im*sh0_d.y;
zh0_d.z += M[1][1].re*sh0_d.z;
zh0_d.z -= M[1][1].im*sh0_d.w;
zh0_d.z += M[1][2].re*sh1_d.x;
zh0_d.z -= M[1][2].im*sh1_d.y;

zh0_d.w  = M[1][0].re*sh0_d.y;
zh0_d.w += M[1][0].im*sh0_d.x;
zh0_d.w += M[1][1].re*sh0_d.w;
zh0_d.w += M[1][1].im*sh0_d.z;
zh0_d.w += M[1][2].re*sh1_d.y;
zh0_d.w += M[1][2].im*sh1_d.x;
(*(out_d+0)).z -= zh0_d.z;
(*(out_d+0)).w -= zh0_d.w;


zh1_u.x  = M[2][0].re*sh0_u.x;
zh1_u.x -= M[2][0].im*sh0_u.y;
zh1_u.x += M[2][1].re*sh0_u.z;
zh1_u.x -= M[2][1].im*sh0_u.w;
zh1_u.x += M[2][2].re*sh1_u.x;
zh1_u.x -= M[2][2].im*sh1_u.y;

zh1_u.y  = M[2][0].re*sh0_u.y;
zh1_u.y += M[2][0].im*sh0_u.x;
zh1_u.y += M[2][1].re*sh0_u.w;
zh1_u.y += M[2][1].im*sh0_u.z;
zh1_u.y += M[2][2].re*sh1_u.y;
zh1_u.y += M[2][2].im*sh1_u.x;
(*(out_u+1)).x -= zh1_u.x;
(*(out_u+1)).y -= zh1_u.y;

zh1_d.x  = M[2][0].re*sh0_d.x;
zh1_d.x -= M[2][0].im*sh0_d.y;
zh1_d.x += M[2][1].re*sh0_d.z;
zh1_d.x -= M[2][1].im*sh0_d.w;
zh1_d.x += M[2][2].re*sh1_d.x;
zh1_d.x -= M[2][2].im*sh1_d.y;

zh1_d.y  = M[2][0].re*sh0_d.y;
zh1_d.y += M[2][0].im*sh0_d.x;
zh1_d.y += M[2][1].re*sh0_d.w;
zh1_d.y += M[2][1].im*sh0_d.z;
zh1_d.y += M[2][2].re*sh1_d.y;
zh1_d.y += M[2][2].im*sh1_d.x;
(*(out_d+1)).x -= zh1_d.x;
(*(out_d+1)).y -= zh1_d.y;


     sh1_u.z = kappa.re*(s1_u.z + s4_u.w);
     sh1_u.z -= kappa.im*(s1_u.w - s4_u.z);
     sh1_u.w = kappa.re*(s1_u.w - s4_u.z);
     sh1_u.w += kappa.im*(s1_u.z + s4_u.w);
         
     sh2_u.x = kappa.re*(s2_u.x + s5_u.y);
     sh2_u.x -= kappa.im*(s2_u.y - s5_u.x);
     sh2_u.y = kappa.re*(s2_u.y - s5_u.x);
     sh2_u.y += kappa.im*(s2_u.x + s5_u.y);
       
     sh2_u.z = kappa.re*(s2_u.z + s5_u.w);
     sh2_u.z -= kappa.im*(s2_u.w - s5_u.z);
     sh2_u.w = kappa.re*(s2_u.w - s5_u.z);
     sh2_u.w += kappa.im*(s2_u.z + s5_u.w);
     
     
     sh1_d.z = kappa.re*(s1_d.z + s4_d.w);
     sh1_d.z -= kappa.im*(s1_d.w - s4_d.z);
     sh1_d.w = kappa.re*(s1_d.w - s4_d.z);
     sh1_d.w += kappa.im*(s1_d.z + s4_d.w);
         
     sh2_d.x = kappa.re*(s2_d.x + s5_d.y);
     sh2_d.x -= kappa.im*(s2_d.y - s5_d.x);
     sh2_d.y = kappa.re*(s2_d.y - s5_d.x);
     sh2_d.y += kappa.im*(s2_d.x + s5_d.y);
       
     sh2_d.z = kappa.re*(s2_d.z + s5_d.w);
     sh2_d.z -= kappa.im*(s2_d.w - s5_d.z);
     sh2_d.w = kappa.re*(s2_d.w - s5_d.z);
     sh2_d.w += kappa.im*(s2_d.z + s5_d.w);     
     

zh1_u.z  = M[0][0].re*sh1_u.z;
zh1_u.z -= M[0][0].im*sh1_u.w;
zh1_u.z += M[0][1].re*sh2_u.x;
zh1_u.z -= M[0][1].im*sh2_u.y;
zh1_u.z += M[0][2].re*sh2_u.z;
zh1_u.z -= M[0][2].im*sh2_u.w;

zh1_u.w  = M[0][0].re*sh1_u.w;
zh1_u.w += M[0][0].im*sh1_u.z;
zh1_u.w += M[0][1].re*sh2_u.y;
zh1_u.w += M[0][1].im*sh2_u.x;
zh1_u.w += M[0][2].re*sh2_u.w;
zh1_u.w += M[0][2].im*sh2_u.z; 
(*(out_u+1)).z -= zh1_u.z;
(*(out_u+1)).w -= zh1_u.w;


zh1_d.z  = M[0][0].re*sh1_d.z;
zh1_d.z -= M[0][0].im*sh1_d.w;
zh1_d.z += M[0][1].re*sh2_d.x;
zh1_d.z -= M[0][1].im*sh2_d.y;
zh1_d.z += M[0][2].re*sh2_d.z;
zh1_d.z -= M[0][2].im*sh2_d.w;

zh1_d.w  = M[0][0].re*sh1_d.w;
zh1_d.w += M[0][0].im*sh1_d.z;
zh1_d.w += M[0][1].re*sh2_d.y;
zh1_d.w += M[0][1].im*sh2_d.x;
zh1_d.w += M[0][2].re*sh2_d.w;
zh1_d.w += M[0][2].im*sh2_d.z; 
(*(out_d+1)).z -= zh1_d.z;
(*(out_d+1)).w -= zh1_d.w;


zh2_u.x  = M[1][0].re*sh1_u.z;
zh2_u.x -= M[1][0].im*sh1_u.w;
zh2_u.x += M[1][1].re*sh2_u.x;
zh2_u.x -= M[1][1].im*sh2_u.y;
zh2_u.x += M[1][2].re*sh2_u.z;
zh2_u.x -= M[1][2].im*sh2_u.w;

zh2_u.y  = M[1][0].re*sh1_u.w;
zh2_u.y += M[1][0].im*sh1_u.z;
zh2_u.y += M[1][1].re*sh2_u.y;
zh2_u.y += M[1][1].im*sh2_u.x;
zh2_u.y += M[1][2].re*sh2_u.w;
zh2_u.y += M[1][2].im*sh2_u.z;
(*(out_u+2)).x -= zh2_u.x;
(*(out_u+2)).y -= zh2_u.y;


zh2_d.x  = M[1][0].re*sh1_d.z;
zh2_d.x -= M[1][0].im*sh1_d.w;
zh2_d.x += M[1][1].re*sh2_d.x;
zh2_d.x -= M[1][1].im*sh2_d.y;
zh2_d.x += M[1][2].re*sh2_d.z;
zh2_d.x -= M[1][2].im*sh2_d.w;

zh2_d.y  = M[1][0].re*sh1_d.w;
zh2_d.y += M[1][0].im*sh1_d.z;
zh2_d.y += M[1][1].re*sh2_d.y;
zh2_d.y += M[1][1].im*sh2_d.x;
zh2_d.y += M[1][2].re*sh2_d.w;
zh2_d.y += M[1][2].im*sh2_d.z;
(*(out_d+2)).x -= zh2_d.x;
(*(out_d+2)).y -= zh2_d.y;


zh2_u.z  = M[2][0].re*sh1_u.z;
zh2_u.z -= M[2][0].im*sh1_u.w;
zh2_u.z += M[2][1].re*sh2_u.x;
zh2_u.z -= M[2][1].im*sh2_u.y;
zh2_u.z += M[2][2].re*sh2_u.z;
zh2_u.z -= M[2][2].im*sh2_u.w; 

zh2_u.w  = M[2][0].re*sh1_u.w;
zh2_u.w += M[2][0].im*sh1_u.z;
zh2_u.w += M[2][1].re*sh2_u.y;
zh2_u.w += M[2][1].im*sh2_u.x;
zh2_u.w += M[2][2].re*sh2_u.w;
zh2_u.w += M[2][2].im*sh2_u.z; 
(*(out_u+2)).z -= zh2_u.z;
(*(out_u+2)).w -= zh2_u.w;


zh2_d.z  = M[2][0].re*sh1_d.z;
zh2_d.z -= M[2][0].im*sh1_d.w;
zh2_d.z += M[2][1].re*sh2_d.x;
zh2_d.z -= M[2][1].im*sh2_d.y;
zh2_d.z += M[2][2].re*sh2_d.z;
zh2_d.z -= M[2][2].im*sh2_d.w; 

zh2_d.w  = M[2][0].re*sh1_d.w;
zh2_d.w += M[2][0].im*sh1_d.z;
zh2_d.w += M[2][1].re*sh2_d.y;
zh2_d.w += M[2][1].im*sh2_d.x;
zh2_d.w += M[2][2].re*sh2_d.w;
zh2_d.w += M[2][2].im*sh2_d.z; 
(*(out_d+2)).z -= zh2_d.z;
(*(out_d+2)).w -= zh2_d.w;


//Reconstruct lower half spinor
//this is just a multiplication by "-i"
     (*(out_u+3)).x -= zh0_u.y;
     (*(out_u+3)).y += zh0_u.x; 

     (*(out_u+3)).z -= zh0_u.w;
     (*(out_u+3)).w += zh0_u.z; 
     
     (*(out_u+4)).x -= zh1_u.y;
     (*(out_u+4)).y += zh1_u.x;       


//this is just a multiplication by "i"     
     (*(out_u+4)).z += zh1_u.w;
     (*(out_u+4)).w -= zh1_u.z;       

     (*(out_u+5)).x += zh2_u.y;
     (*(out_u+5)).y -= zh2_u.x;  
     
     (*(out_u+5)).z += zh2_u.w;
     (*(out_u+5)).w -= zh2_u.z;
  
//Reconstruct lower half spinor
//this is just a multiplication by "-i"
     (*(out_d+3)).x -= zh0_d.y;
     (*(out_d+3)).y += zh0_d.x; 

     (*(out_d+3)).z -= zh0_d.w;
     (*(out_d+3)).w += zh0_d.z; 
     
     (*(out_d+4)).x -= zh1_d.y;
     (*(out_d+4)).y += zh1_d.x;       


//this is just a multiplication by "i"     
     (*(out_d+4)).z += zh1_d.w;
     (*(out_d+4)).w -= zh1_d.z;       

     (*(out_d+5)).x += zh2_d.y;
     (*(out_d+5)).y -= zh2_d.x;  
     
     (*(out_d+5)).z += zh2_d.w;
     (*(out_d+5)).w -= zh2_d.z;     
     
     
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
//-conj(kappa)(r + gamma_mu) kappa reell !!!
__device__ void dev_su3MtV_kappaP3_minus(dev_su3 M, const dev_spinor * s, dev_spinor * out, dev_complex kappa){

  dev_spinor sh0,sh1,sh2,zh0,zh1,zh2;
  
  
  //first apply Projector on upper halfspinor 
     sh0.x = kappa.re*( (*(s+0*DEVOFF)).x + (*(s+3*DEVOFF)).y);
     sh0.x += kappa.im*( (*(s+0*DEVOFF)).y - (*(s+3*DEVOFF)).x);
     sh0.y = kappa.re*( (*(s+0*DEVOFF)).y - (*(s+3*DEVOFF)).x);
     sh0.y -= kappa.im*( (*(s+0*DEVOFF)).x + (*(s+3*DEVOFF)).y);
   
     sh0.z = kappa.re*( (*(s+0*DEVOFF)).z + (*(s+3*DEVOFF)).w);
     sh0.z += kappa.im*( (*(s+0*DEVOFF)).w - (*(s+3*DEVOFF)).z);
     sh0.w = kappa.re*( (*(s+0*DEVOFF)).w - (*(s+3*DEVOFF)).z);
     sh0.w -= kappa.im*( (*(s+0*DEVOFF)).z + (*(s+3*DEVOFF)).w);
   
     sh1.x = kappa.re*( (*(s+1*DEVOFF)).x + (*(s+4*DEVOFF)).y);
     sh1.x += kappa.im*( (*(s+1*DEVOFF)).y - (*(s+4*DEVOFF)).x);
     sh1.y = kappa.re*( (*(s+1*DEVOFF)).y - (*(s+4*DEVOFF)).x);
     sh1.y -= kappa.im*( (*(s+1*DEVOFF)).x + (*(s+4*DEVOFF)).y);  


  

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


     sh1.z = kappa.re*( (*(s+1*DEVOFF)).z - (*(s+4*DEVOFF)).w);
     sh1.z += kappa.im*( (*(s+1*DEVOFF)).w + (*(s+4*DEVOFF)).z);
     sh1.w = kappa.re*( (*(s+1*DEVOFF)).w + (*(s+4*DEVOFF)).z);
     sh1.w -= kappa.im*( (*(s+1*DEVOFF)).z - (*(s+4*DEVOFF)).w);
     
     sh2.x = kappa.re*( (*(s+2*DEVOFF)).x - (*(s+5*DEVOFF)).y);
     sh2.x += kappa.im*( (*(s+2*DEVOFF)).y + (*(s+5*DEVOFF)).x);
     sh2.y = kappa.re*( (*(s+2*DEVOFF)).y + (*(s+5*DEVOFF)).x);
     sh2.y -= kappa.im*( (*(s+2*DEVOFF)).x - (*(s+5*DEVOFF)).y);  
     
     sh2.z = kappa.re*( (*(s+2*DEVOFF)).z - (*(s+5*DEVOFF)).w);
     sh2.z += kappa.im*( (*(s+2*DEVOFF)).w + (*(s+5*DEVOFF)).z);
     sh2.w = kappa.re*( (*(s+2*DEVOFF)).w + (*(s+5*DEVOFF)).z);
     sh2.w -= kappa.im*( (*(s+2*DEVOFF)).z - (*(s+5*DEVOFF)).w);
     
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


__device__ void dev_su3MtV_kappaP3_minus_spintex(dev_su3 M, int pos, dev_spinor * out, dev_complex kappa){

  dev_spinor sh0,sh1,sh2,zh0,zh1,zh2;
  
  dev_spinor s0, s1, s2, s3, s4, s5;
  s0 = tex1Dfetch(spin_tex0,pos);
  s1 = tex1Dfetch(spin_tex1,pos);
  s2 = tex1Dfetch(spin_tex2,pos);
  s3 = tex1Dfetch(spin_tex3,pos);
  s4 = tex1Dfetch(spin_tex4,pos); 
  s5 = tex1Dfetch(spin_tex5,pos); 
  
  //first apply Projector on upper halfspinor 
     sh0.x = kappa.re*(s0.x + s3.y);
     sh0.x += kappa.im*(s0.y - s3.x);
     sh0.y = kappa.re*(s0.y - s3.x);
     sh0.y -= kappa.im*(s0.x + s3.y);
   
     sh0.z = kappa.re*(s0.z + s3.w);
     sh0.z += kappa.im*(s0.w - s3.z);
     sh0.w = kappa.re*(s0.w - s3.z);
     sh0.w -= kappa.im*(s0.z + s3.w);
   
     sh1.x = kappa.re*(s1.x + s4.y);
     sh1.x += kappa.im*(s1.y - s4.x);
     sh1.y = kappa.re*(s1.y - s4.x);
     sh1.y -= kappa.im*(s1.x + s4.y);   


  

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


     sh1.z = kappa.re*(s1.z - s4.w);
     sh1.z += kappa.im*(s1.w + s4.z);
     sh1.w = kappa.re*(s1.w + s4.z);
     sh1.w -= kappa.im*(s1.z - s4.w);
         
     sh2.x = kappa.re*(s2.x - s5.y);
     sh2.x += kappa.im*(s2.y + s5.x);
     sh2.y = kappa.re*(s2.y + s5.x);
     sh2.y -= kappa.im*(s2.x - s5.y);
       
     sh2.z = kappa.re*(s2.z - s5.w);
     sh2.z += kappa.im*(s2.w + s5.z);
     sh2.w = kappa.re*(s2.w + s5.z);
     sh2.w -= kappa.im*(s2.z - s5.w);

     
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



__device__ void dev_su3MtV_kappaP3_minus_spintex2(dev_su3 M, int pos, dev_spinor * out, dev_complex kappa){

  dev_spinor sh0,sh1,sh2,zh0,zh1,zh2;
  
  dev_spinor s0, s1, s2, s3, s4, s5;
  s0 = tex1Dfetch(spin_tex_dn0,pos);
  s1 = tex1Dfetch(spin_tex_dn1,pos);
  s2 = tex1Dfetch(spin_tex_dn2,pos);
  s3 = tex1Dfetch(spin_tex_dn3,pos);
  s4 = tex1Dfetch(spin_tex_dn4,pos); 
  s5 = tex1Dfetch(spin_tex_dn5,pos); 
  
  //first apply Projector on upper halfspinor 
     sh0.x = kappa.re*(s0.x + s3.y);
     sh0.x += kappa.im*(s0.y - s3.x);
     sh0.y = kappa.re*(s0.y - s3.x);
     sh0.y -= kappa.im*(s0.x + s3.y);
   
     sh0.z = kappa.re*(s0.z + s3.w);
     sh0.z += kappa.im*(s0.w - s3.z);
     sh0.w = kappa.re*(s0.w - s3.z);
     sh0.w -= kappa.im*(s0.z + s3.w);
   
     sh1.x = kappa.re*(s1.x + s4.y);
     sh1.x += kappa.im*(s1.y - s4.x);
     sh1.y = kappa.re*(s1.y - s4.x);
     sh1.y -= kappa.im*(s1.x + s4.y);   


  

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


     sh1.z = kappa.re*(s1.z - s4.w);
     sh1.z += kappa.im*(s1.w + s4.z);
     sh1.w = kappa.re*(s1.w + s4.z);
     sh1.w -= kappa.im*(s1.z - s4.w);
         
     sh2.x = kappa.re*(s2.x - s5.y);
     sh2.x += kappa.im*(s2.y + s5.x);
     sh2.y = kappa.re*(s2.y + s5.x);
     sh2.y -= kappa.im*(s2.x - s5.y);
       
     sh2.z = kappa.re*(s2.z - s5.w);
     sh2.z += kappa.im*(s2.w + s5.z);
     sh2.w = kappa.re*(s2.w + s5.z);
     sh2.w -= kappa.im*(s2.z - s5.w);

     
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





__device__ void dev_su3MtV_kappaP3_minus_spintex_ud(dev_su3 M, int pos, dev_spinor * out_u, dev_spinor * out_d, dev_complex kappa){

  dev_spinor sh0_u,sh1_u,sh2_u,zh0_u,zh1_u,zh2_u;
  dev_spinor s0_u, s1_u, s2_u, s3_u, s4_u, s5_u;
  
  dev_spinor sh0_d,sh1_d,sh2_d,zh0_d,zh1_d,zh2_d;
  dev_spinor s0_d, s1_d, s2_d, s3_d, s4_d, s5_d; 

  
  s0_u = tex1Dfetch(spin_tex0,pos);
  s1_u = tex1Dfetch(spin_tex1,pos);
  s2_u = tex1Dfetch(spin_tex2,pos);
  s3_u = tex1Dfetch(spin_tex3,pos);
  s4_u = tex1Dfetch(spin_tex4,pos); 
  s5_u = tex1Dfetch(spin_tex5,pos); 
  
  s0_d = tex1Dfetch(spin_tex_dn0,pos);
  s1_d = tex1Dfetch(spin_tex_dn1,pos);
  s2_d = tex1Dfetch(spin_tex_dn2,pos);
  s3_d = tex1Dfetch(spin_tex_dn3,pos);
  s4_d = tex1Dfetch(spin_tex_dn4,pos); 
  s5_d = tex1Dfetch(spin_tex_dn5,pos);  
  
  
  //first apply Projector on upper halfspinor 
     sh0_u.x = kappa.re*(s0_u.x + s3_u.y);
     sh0_u.x += kappa.im*(s0_u.y - s3_u.x);
     sh0_u.y = kappa.re*(s0_u.y - s3_u.x);
     sh0_u.y -= kappa.im*(s0_u.x + s3_u.y);
   
     sh0_u.z = kappa.re*(s0_u.z + s3_u.w);
     sh0_u.z += kappa.im*(s0_u.w - s3_u.z);
     sh0_u.w = kappa.re*(s0_u.w - s3_u.z);
     sh0_u.w -= kappa.im*(s0_u.z + s3_u.w);
   
     sh1_u.x = kappa.re*(s1_u.x + s4_u.y);
     sh1_u.x += kappa.im*(s1_u.y - s4_u.x);
     sh1_u.y = kappa.re*(s1_u.y - s4_u.x);
     sh1_u.y -= kappa.im*(s1_u.x + s4_u.y);   

     sh0_d.x = kappa.re*(s0_d.x + s3_d.y);
     sh0_d.x += kappa.im*(s0_d.y - s3_d.x);
     sh0_d.y = kappa.re*(s0_d.y - s3_d.x);
     sh0_d.y -= kappa.im*(s0_d.x + s3_d.y);
   
     sh0_d.z = kappa.re*(s0_d.z + s3_d.w);
     sh0_d.z += kappa.im*(s0_d.w - s3_d.z);
     sh0_d.w = kappa.re*(s0_d.w - s3_d.z);
     sh0_d.w -= kappa.im*(s0_d.z + s3_d.w);
   
     sh1_d.x = kappa.re*(s1_d.x + s4_d.y);
     sh1_d.x += kappa.im*(s1_d.y - s4_d.x);
     sh1_d.y = kappa.re*(s1_d.y - s4_d.x);
     sh1_d.y -= kappa.im*(s1_d.x + s4_d.y);   
  

//Multiply by gauge field  
zh0_u.x  = M[0][0].re*sh0_u.x;
zh0_u.x -= M[0][0].im*sh0_u.y;
zh0_u.x += M[0][1].re*sh0_u.z;
zh0_u.x -= M[0][1].im*sh0_u.w;
zh0_u.x += M[0][2].re*sh1_u.x;
zh0_u.x -= M[0][2].im*sh1_u.y;

zh0_u.y  = M[0][0].re*sh0_u.y;
zh0_u.y += M[0][0].im*sh0_u.x;
zh0_u.y += M[0][1].re*sh0_u.w;
zh0_u.y += M[0][1].im*sh0_u.z;
zh0_u.y += M[0][2].re*sh1_u.y;
zh0_u.y += M[0][2].im*sh1_u.x;
(*(out_u+0)).x -= zh0_u.x;
(*(out_u+0)).y -= zh0_u.y;

zh0_d.x  = M[0][0].re*sh0_d.x;
zh0_d.x -= M[0][0].im*sh0_d.y;
zh0_d.x += M[0][1].re*sh0_d.z;
zh0_d.x -= M[0][1].im*sh0_d.w;
zh0_d.x += M[0][2].re*sh1_d.x;
zh0_d.x -= M[0][2].im*sh1_d.y;

zh0_d.y  = M[0][0].re*sh0_d.y;
zh0_d.y += M[0][0].im*sh0_d.x;
zh0_d.y += M[0][1].re*sh0_d.w;
zh0_d.y += M[0][1].im*sh0_d.z;
zh0_d.y += M[0][2].re*sh1_d.y;
zh0_d.y += M[0][2].im*sh1_d.x;
(*(out_d+0)).x -= zh0_d.x;
(*(out_d+0)).y -= zh0_d.y;


zh0_u.z  = M[1][0].re*sh0_u.x;
zh0_u.z -= M[1][0].im*sh0_u.y;
zh0_u.z += M[1][1].re*sh0_u.z;
zh0_u.z -= M[1][1].im*sh0_u.w;
zh0_u.z += M[1][2].re*sh1_u.x;
zh0_u.z -= M[1][2].im*sh1_u.y;

zh0_u.w  = M[1][0].re*sh0_u.y;
zh0_u.w += M[1][0].im*sh0_u.x;
zh0_u.w += M[1][1].re*sh0_u.w;
zh0_u.w += M[1][1].im*sh0_u.z;
zh0_u.w += M[1][2].re*sh1_u.y;
zh0_u.w += M[1][2].im*sh1_u.x;
(*(out_u+0)).z -= zh0_u.z;
(*(out_u+0)).w -= zh0_u.w;


zh0_d.z  = M[1][0].re*sh0_d.x;
zh0_d.z -= M[1][0].im*sh0_d.y;
zh0_d.z += M[1][1].re*sh0_d.z;
zh0_d.z -= M[1][1].im*sh0_d.w;
zh0_d.z += M[1][2].re*sh1_d.x;
zh0_d.z -= M[1][2].im*sh1_d.y;

zh0_d.w  = M[1][0].re*sh0_d.y;
zh0_d.w += M[1][0].im*sh0_d.x;
zh0_d.w += M[1][1].re*sh0_d.w;
zh0_d.w += M[1][1].im*sh0_d.z;
zh0_d.w += M[1][2].re*sh1_d.y;
zh0_d.w += M[1][2].im*sh1_d.x;
(*(out_d+0)).z -= zh0_d.z;
(*(out_d+0)).w -= zh0_d.w;


zh1_u.x  = M[2][0].re*sh0_u.x;
zh1_u.x -= M[2][0].im*sh0_u.y;
zh1_u.x += M[2][1].re*sh0_u.z;
zh1_u.x -= M[2][1].im*sh0_u.w;
zh1_u.x += M[2][2].re*sh1_u.x;
zh1_u.x -= M[2][2].im*sh1_u.y;

zh1_u.y  = M[2][0].re*sh0_u.y;
zh1_u.y += M[2][0].im*sh0_u.x;
zh1_u.y += M[2][1].re*sh0_u.w;
zh1_u.y += M[2][1].im*sh0_u.z;
zh1_u.y += M[2][2].re*sh1_u.y;
zh1_u.y += M[2][2].im*sh1_u.x;
(*(out_u+1)).x -= zh1_u.x;
(*(out_u+1)).y -= zh1_u.y;


zh1_d.x  = M[2][0].re*sh0_d.x;
zh1_d.x -= M[2][0].im*sh0_d.y;
zh1_d.x += M[2][1].re*sh0_d.z;
zh1_d.x -= M[2][1].im*sh0_d.w;
zh1_d.x += M[2][2].re*sh1_d.x;
zh1_d.x -= M[2][2].im*sh1_d.y;

zh1_d.y  = M[2][0].re*sh0_d.y;
zh1_d.y += M[2][0].im*sh0_d.x;
zh1_d.y += M[2][1].re*sh0_d.w;
zh1_d.y += M[2][1].im*sh0_d.z;
zh1_d.y += M[2][2].re*sh1_d.y;
zh1_d.y += M[2][2].im*sh1_d.x;
(*(out_d+1)).x -= zh1_d.x;
(*(out_d+1)).y -= zh1_d.y;

     sh1_u.z = kappa.re*(s1_u.z - s4_u.w);
     sh1_u.z += kappa.im*(s1_u.w + s4_u.z);
     sh1_u.w = kappa.re*(s1_u.w + s4_u.z);
     sh1_u.w -= kappa.im*(s1_u.z - s4_u.w);
         
     sh2_u.x = kappa.re*(s2_u.x - s5_u.y);
     sh2_u.x += kappa.im*(s2_u.y + s5_u.x);
     sh2_u.y = kappa.re*(s2_u.y + s5_u.x);
     sh2_u.y -= kappa.im*(s2_u.x - s5_u.y);
       
     sh2_u.z = kappa.re*(s2_u.z - s5_u.w);
     sh2_u.z += kappa.im*(s2_u.w + s5_u.z);
     sh2_u.w = kappa.re*(s2_u.w + s5_u.z);
     sh2_u.w -= kappa.im*(s2_u.z - s5_u.w);

     sh1_d.z = kappa.re*(s1_d.z - s4_d.w);
     sh1_d.z += kappa.im*(s1_d.w + s4_d.z);
     sh1_d.w = kappa.re*(s1_d.w + s4_d.z);
     sh1_d.w -= kappa.im*(s1_d.z - s4_d.w);
         
     sh2_d.x = kappa.re*(s2_d.x - s5_d.y);
     sh2_d.x += kappa.im*(s2_d.y + s5_d.x);
     sh2_d.y = kappa.re*(s2_d.y + s5_d.x);
     sh2_d.y -= kappa.im*(s2_d.x - s5_d.y);
       
     sh2_d.z = kappa.re*(s2_d.z - s5_d.w);
     sh2_d.z += kappa.im*(s2_d.w + s5_d.z);
     sh2_d.w = kappa.re*(s2_d.w + s5_d.z);
     sh2_d.w -= kappa.im*(s2_d.z - s5_d.w);
     
     
zh1_u.z  = M[0][0].re*sh1_u.z;
zh1_u.z -= M[0][0].im*sh1_u.w;
zh1_u.z += M[0][1].re*sh2_u.x;
zh1_u.z -= M[0][1].im*sh2_u.y;
zh1_u.z += M[0][2].re*sh2_u.z;
zh1_u.z -= M[0][2].im*sh2_u.w;

zh1_u.w  = M[0][0].re*sh1_u.w;
zh1_u.w += M[0][0].im*sh1_u.z;
zh1_u.w += M[0][1].re*sh2_u.y;
zh1_u.w += M[0][1].im*sh2_u.x;
zh1_u.w += M[0][2].re*sh2_u.w;
zh1_u.w += M[0][2].im*sh2_u.z;
(*(out_u+1)).z -= zh1_u.z;
(*(out_u+1)).w -= zh1_u.w;

zh1_d.z  = M[0][0].re*sh1_d.z;
zh1_d.z -= M[0][0].im*sh1_d.w;
zh1_d.z += M[0][1].re*sh2_d.x;
zh1_d.z -= M[0][1].im*sh2_d.y;
zh1_d.z += M[0][2].re*sh2_d.z;
zh1_d.z -= M[0][2].im*sh2_d.w;

zh1_d.w  = M[0][0].re*sh1_d.w;
zh1_d.w += M[0][0].im*sh1_d.z;
zh1_d.w += M[0][1].re*sh2_d.y;
zh1_d.w += M[0][1].im*sh2_d.x;
zh1_d.w += M[0][2].re*sh2_d.w;
zh1_d.w += M[0][2].im*sh2_d.z;
(*(out_d+1)).z -= zh1_d.z;
(*(out_d+1)).w -= zh1_d.w;

zh2_u.x  = M[1][0].re*sh1_u.z;
zh2_u.x -= M[1][0].im*sh1_u.w;
zh2_u.x += M[1][1].re*sh2_u.x;
zh2_u.x -= M[1][1].im*sh2_u.y;
zh2_u.x += M[1][2].re*sh2_u.z;
zh2_u.x -= M[1][2].im*sh2_u.w;

zh2_u.y  = M[1][0].re*sh1_u.w;
zh2_u.y += M[1][0].im*sh1_u.z;
zh2_u.y += M[1][1].re*sh2_u.y;
zh2_u.y += M[1][1].im*sh2_u.x;
zh2_u.y += M[1][2].re*sh2_u.w;
zh2_u.y += M[1][2].im*sh2_u.z;
(*(out_u+2)).x -= zh2_u.x;
(*(out_u+2)).y -= zh2_u.y;

zh2_d.x  = M[1][0].re*sh1_d.z;
zh2_d.x -= M[1][0].im*sh1_d.w;
zh2_d.x += M[1][1].re*sh2_d.x;
zh2_d.x -= M[1][1].im*sh2_d.y;
zh2_d.x += M[1][2].re*sh2_d.z;
zh2_d.x -= M[1][2].im*sh2_d.w;

zh2_d.y  = M[1][0].re*sh1_d.w;
zh2_d.y += M[1][0].im*sh1_d.z;
zh2_d.y += M[1][1].re*sh2_d.y;
zh2_d.y += M[1][1].im*sh2_d.x;
zh2_d.y += M[1][2].re*sh2_d.w;
zh2_d.y += M[1][2].im*sh2_d.z;
(*(out_d+2)).x -= zh2_d.x;
(*(out_d+2)).y -= zh2_d.y;

zh2_u.z  = M[2][0].re*sh1_u.z;
zh2_u.z -= M[2][0].im*sh1_u.w;
zh2_u.z += M[2][1].re*sh2_u.x;
zh2_u.z -= M[2][1].im*sh2_u.y;
zh2_u.z += M[2][2].re*sh2_u.z;
zh2_u.z -= M[2][2].im*sh2_u.w;

zh2_u.w  = M[2][0].re*sh1_u.w;
zh2_u.w += M[2][0].im*sh1_u.z;
zh2_u.w += M[2][1].re*sh2_u.y;
zh2_u.w += M[2][1].im*sh2_u.x;
zh2_u.w += M[2][2].re*sh2_u.w;
zh2_u.w += M[2][2].im*sh2_u.z;
(*(out_u+2)).z -= zh2_u.z;
(*(out_u+2)).w -= zh2_u.w;

zh2_d.z  = M[2][0].re*sh1_d.z;
zh2_d.z -= M[2][0].im*sh1_d.w;
zh2_d.z += M[2][1].re*sh2_d.x;
zh2_d.z -= M[2][1].im*sh2_d.y;
zh2_d.z += M[2][2].re*sh2_d.z;
zh2_d.z -= M[2][2].im*sh2_d.w;

zh2_d.w  = M[2][0].re*sh1_d.w;
zh2_d.w += M[2][0].im*sh1_d.z;
zh2_d.w += M[2][1].re*sh2_d.y;
zh2_d.w += M[2][1].im*sh2_d.x;
zh2_d.w += M[2][2].re*sh2_d.w;
zh2_d.w += M[2][2].im*sh2_d.z;
(*(out_d+2)).z -= zh2_d.z;
(*(out_d+2)).w -= zh2_d.w;

//Reconstruct lower half spinor
//this is just a multiplication by "i"
     (*(out_u+3)).x += zh0_u.y;
     (*(out_u+3)).y -= zh0_u.x; 

     (*(out_u+3)).z += zh0_u.w;
     (*(out_u+3)).w -= zh0_u.z; 
     
     (*(out_u+4)).x += zh1_u.y;
     (*(out_u+4)).y -= zh1_u.x;       


//this is just a multiplication by "-i"     
     (*(out_u+4)).z -= zh1_u.w;
     (*(out_u+4)).w += zh1_u.z;       

     (*(out_u+5)).x -= zh2_u.y;
     (*(out_u+5)).y += zh2_u.x;  
     
     (*(out_u+5)).z -= zh2_u.w;
     (*(out_u+5)).w += zh2_u.z;

     
//Reconstruct lower half spinor
//this is just a multiplication by "i"
     (*(out_d+3)).x += zh0_d.y;
     (*(out_d+3)).y -= zh0_d.x; 

     (*(out_d+3)).z += zh0_d.w;
     (*(out_d+3)).w -= zh0_d.z; 
     
     (*(out_d+4)).x += zh1_d.y;
     (*(out_d+4)).y -= zh1_d.x;       


//this is just a multiplication by "-i"     
     (*(out_d+4)).z -= zh1_d.w;
     (*(out_d+4)).w += zh1_d.z;       

     (*(out_d+5)).x -= zh2_d.y;
     (*(out_d+5)).y += zh2_d.x;  
     
     (*(out_d+5)).z -= zh2_d.w;
     (*(out_d+5)).w += zh2_d.z;     
     
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
__device__ void dev_su3MtV_kappaP2_plus(dev_su3 M, const dev_spinor * s, dev_spinor * out, dev_complex kappa){

  
  dev_spinor sh0,sh1,sh2,zh0,zh1,zh2;


     sh0.x = kappa.re*( (*(s+0*DEVOFF)).x + (*(s+4*DEVOFF)).z);
     sh0.x -= kappa.im*( (*(s+0*DEVOFF)).y + (*(s+4*DEVOFF)).w);
     sh0.y = kappa.re*( (*(s+0*DEVOFF)).y + (*(s+4*DEVOFF)).w);
     sh0.y += kappa.im*( (*(s+0*DEVOFF)).x + (*(s+4*DEVOFF)).z);
   
     
     sh0.z = kappa.re*( (*(s+0*DEVOFF)).z + (*(s+5*DEVOFF)).x);
     sh0.z -= kappa.im*( (*(s+0*DEVOFF)).w + (*(s+5*DEVOFF)).y);
     sh0.w = kappa.re*( (*(s+0*DEVOFF)).w + (*(s+5*DEVOFF)).y);
     sh0.w += kappa.im*( (*(s+0*DEVOFF)).z + (*(s+5*DEVOFF)).x);
       
     
     sh1.x = kappa.re*( (*(s+1*DEVOFF)).x + (*(s+5*DEVOFF)).z);
     sh1.x -= kappa.im*( (*(s+1*DEVOFF)).y + (*(s+5*DEVOFF)).w);
     sh1.y = kappa.re*( (*(s+1*DEVOFF)).y + (*(s+5*DEVOFF)).w);
     sh1.y += kappa.im*( (*(s+1*DEVOFF)).x + (*(s+5*DEVOFF)).z);


     
     
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


     sh1.z = kappa.re*( (*(s+1*DEVOFF)).z - (*(s+3*DEVOFF)).x);
     sh1.z -= kappa.im*( (*(s+1*DEVOFF)).w - (*(s+3*DEVOFF)).y);
     sh1.w = kappa.re*( (*(s+1*DEVOFF)).w - (*(s+3*DEVOFF)).y);
     sh1.w += kappa.im*( (*(s+1*DEVOFF)).z - (*(s+3*DEVOFF)).x);
     
     
     sh2.x = kappa.re*( (*(s+2*DEVOFF)).x - (*(s+3*DEVOFF)).z);
     sh2.x -= kappa.im*( (*(s+2*DEVOFF)).y - (*(s+3*DEVOFF)).w);
     sh2.y = kappa.re*( (*(s+2*DEVOFF)).y - (*(s+3*DEVOFF)).w);
     sh2.y += kappa.im*( (*(s+2*DEVOFF)).x - (*(s+3*DEVOFF)).z);
          
     
     sh2.z = kappa.re*( (*(s+2*DEVOFF)).z - (*(s+4*DEVOFF)).x);
     sh2.z -= kappa.im*( (*(s+2*DEVOFF)).w - (*(s+4*DEVOFF)).y);
     sh2.w = kappa.re*( (*(s+2*DEVOFF)).w - (*(s+4*DEVOFF)).y);
     sh2.w += kappa.im*( (*(s+2*DEVOFF)).z - (*(s+4*DEVOFF)).x);
     
     
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


__device__ void dev_su3MtV_kappaP2_plus_spintex(dev_su3 M, int pos, dev_spinor * out, dev_complex kappa){

  
  dev_spinor sh0,sh1,sh2,zh0,zh1,zh2;
  
  dev_spinor s0, s1, s2, s3, s4, s5;
  s0 = tex1Dfetch(spin_tex0,pos);
  s1 = tex1Dfetch(spin_tex1,pos);
  s2 = tex1Dfetch(spin_tex2,pos);
  s3 = tex1Dfetch(spin_tex3,pos);
  s4 = tex1Dfetch(spin_tex4,pos); 
  s5 = tex1Dfetch(spin_tex5,pos); 

     sh0.x = kappa.re*(s0.x + s4.z);
     sh0.x -= kappa.im*(s0.y + s4.w);
     sh0.y = kappa.re*(s0.y + s4.w);
     sh0.y += kappa.im*(s0.x + s4.z);
   
     
     sh0.z = kappa.re*(s0.z + s5.x);
     sh0.z -= kappa.im*(s0.w + s5.y);
     sh0.w = kappa.re*(s0.w + s5.y);
     sh0.w += kappa.im*(s0.z + s5.x);
       
     
     sh1.x = kappa.re*(s1.x + s5.z);
     sh1.x -= kappa.im*(s1.y + s5.w);
     sh1.y = kappa.re*(s1.y + s5.w);
     sh1.y += kappa.im*(s1.x + s5.z);


     
     
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


     sh1.z = kappa.re*(s1.z - s3.x);
     sh1.z -= kappa.im*(s1.w - s3.y);
     sh1.w = kappa.re*(s1.w - s3.y);
     sh1.w += kappa.im*(s1.z - s3.x);
     
     
     sh2.x = kappa.re*(s2.x - s3.z);
     sh2.x -= kappa.im*(s2.y - s3.w);
     sh2.y = kappa.re*(s2.y - s3.w);
     sh2.y += kappa.im*(s2.x - s3.z);
     
          
     sh2.z = kappa.re*(s2.z - s4.x);
     sh2.z -= kappa.im*(s2.w - s4.y);
     sh2.w = kappa.re*(s2.w - s4.y);
     sh2.w += kappa.im*(s2.z - s4.x); 

     
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


__device__ void dev_su3MtV_kappaP2_plus_spintex2(dev_su3 M, int pos, dev_spinor * out, dev_complex kappa){

  
  dev_spinor sh0,sh1,sh2,zh0,zh1,zh2;
  
  dev_spinor s0, s1, s2, s3, s4, s5;
  s0 = tex1Dfetch(spin_tex_dn0,pos);
  s1 = tex1Dfetch(spin_tex_dn1,pos);
  s2 = tex1Dfetch(spin_tex_dn2,pos);
  s3 = tex1Dfetch(spin_tex_dn3,pos);
  s4 = tex1Dfetch(spin_tex_dn4,pos); 
  s5 = tex1Dfetch(spin_tex_dn5,pos); 

     sh0.x = kappa.re*(s0.x + s4.z);
     sh0.x -= kappa.im*(s0.y + s4.w);
     sh0.y = kappa.re*(s0.y + s4.w);
     sh0.y += kappa.im*(s0.x + s4.z);
   
     
     sh0.z = kappa.re*(s0.z + s5.x);
     sh0.z -= kappa.im*(s0.w + s5.y);
     sh0.w = kappa.re*(s0.w + s5.y);
     sh0.w += kappa.im*(s0.z + s5.x);
       
     
     sh1.x = kappa.re*(s1.x + s5.z);
     sh1.x -= kappa.im*(s1.y + s5.w);
     sh1.y = kappa.re*(s1.y + s5.w);
     sh1.y += kappa.im*(s1.x + s5.z);


     
     
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


     sh1.z = kappa.re*(s1.z - s3.x);
     sh1.z -= kappa.im*(s1.w - s3.y);
     sh1.w = kappa.re*(s1.w - s3.y);
     sh1.w += kappa.im*(s1.z - s3.x);
     
     
     sh2.x = kappa.re*(s2.x - s3.z);
     sh2.x -= kappa.im*(s2.y - s3.w);
     sh2.y = kappa.re*(s2.y - s3.w);
     sh2.y += kappa.im*(s2.x - s3.z);
     
          
     sh2.z = kappa.re*(s2.z - s4.x);
     sh2.z -= kappa.im*(s2.w - s4.y);
     sh2.w = kappa.re*(s2.w - s4.y);
     sh2.w += kappa.im*(s2.z - s4.x); 

     
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
//-cconj(kappa)(r + gamma_mu) kappa reell !!!!
__device__ void dev_su3MtV_kappaP2_minus(dev_su3 M, const dev_spinor * s, dev_spinor * out, dev_complex kappa){

  
  dev_spinor sh0,sh1,sh2,zh0,zh1,zh2;


     sh0.x = kappa.re*( (*(s+0*DEVOFF)).x - (*(s+4*DEVOFF)).z);
     sh0.x += kappa.im*( (*(s+0*DEVOFF)).y - (*(s+4*DEVOFF)).w);
     sh0.y = kappa.re*( (*(s+0*DEVOFF)).y - (*(s+4*DEVOFF)).w);
     sh0.y -= kappa.im*( (*(s+0*DEVOFF)).x - (*(s+4*DEVOFF)).z);
   
     
     sh0.z = kappa.re*( (*(s+0*DEVOFF)).z - (*(s+5*DEVOFF)).x);
     sh0.z += kappa.im*( (*(s+0*DEVOFF)).w - (*(s+5*DEVOFF)).y);
     sh0.w = kappa.re*( (*(s+0*DEVOFF)).w - (*(s+5*DEVOFF)).y);
     sh0.w -= kappa.im*( (*(s+0*DEVOFF)).z - (*(s+5*DEVOFF)).x);
     
       
     sh1.x = kappa.re*( (*(s+1*DEVOFF)).x - (*(s+5*DEVOFF)).z);
     sh1.x += kappa.im*( (*(s+1*DEVOFF)).y - (*(s+5*DEVOFF)).w);
     sh1.y = kappa.re*( (*(s+1*DEVOFF)).y - (*(s+5*DEVOFF)).w);
     sh1.y -= kappa.im*( (*(s+1*DEVOFF)).x - (*(s+5*DEVOFF)).z);


     
     
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


     sh1.z = kappa.re*( (*(s+1*DEVOFF)).z + (*(s+3*DEVOFF)).x);
     sh1.z += kappa.im*( (*(s+1*DEVOFF)).w + (*(s+3*DEVOFF)).y);
     sh1.w = kappa.re*( (*(s+1*DEVOFF)).w + (*(s+3*DEVOFF)).y);
     sh1.w -= kappa.im*( (*(s+1*DEVOFF)).z + (*(s+3*DEVOFF)).x);
     
     sh2.x = kappa.re*( (*(s+2*DEVOFF)).x + (*(s+3*DEVOFF)).z);
     sh2.x += kappa.im*( (*(s+2*DEVOFF)).y + (*(s+3*DEVOFF)).w);
     sh2.y = kappa.re*( (*(s+2*DEVOFF)).y + (*(s+3*DEVOFF)).w);
     sh2.y -= kappa.im*( (*(s+2*DEVOFF)).x + (*(s+3*DEVOFF)).z);
          
     
     sh2.z = kappa.re*( (*(s+2*DEVOFF)).z + (*(s+4*DEVOFF)).x);
     sh2.z += kappa.im*( (*(s+2*DEVOFF)).w + (*(s+4*DEVOFF)).y);
     sh2.w = kappa.re*( (*(s+2*DEVOFF)).w + (*(s+4*DEVOFF)).y);
     sh2.w -= kappa.im*( (*(s+2*DEVOFF)).z + (*(s+4*DEVOFF)).x);

     

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


__device__ void dev_su3MtV_kappaP2_minus_spintex(dev_su3 M, int pos, dev_spinor * out, dev_complex kappa){

  
  dev_spinor sh0,sh1,sh2,zh0,zh1,zh2;

  dev_spinor s0, s1, s2, s3, s4, s5;
  s0 = tex1Dfetch(spin_tex0,pos);
  s1 = tex1Dfetch(spin_tex1,pos);
  s2 = tex1Dfetch(spin_tex2,pos);
  s3 = tex1Dfetch(spin_tex3,pos);
  s4 = tex1Dfetch(spin_tex4,pos); 
  s5 = tex1Dfetch(spin_tex5,pos); 
  
     sh0.x = kappa.re*(s0.x - s4.z);
     sh0.x += kappa.im*(s0.y - s4.w);
     sh0.y = kappa.re*(s0.y - s4.w);
     sh0.y -= kappa.im*(s0.x - s4.z);
     
   
     sh0.z = kappa.re*(s0.z - s5.x);
     sh0.z += kappa.im*(s0.w - s5.y); 
     sh0.w = kappa.re*(s0.w - s5.y);
     sh0.w -= kappa.im*(s0.z - s5.x);
     
       
     sh1.x = kappa.re*(s1.x - s5.z);
     sh1.x += kappa.im*(s1.y - s5.w);
     sh1.y = kappa.re*(s1.y - s5.w);
     sh1.y -= kappa.im*(s1.x - s5.z);
     


     
     
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


     sh1.z = kappa.re*(s1.z + s3.x);
     sh1.z += kappa.im*(s1.w + s3.y);
     sh1.w = kappa.re*(s1.w + s3.y);
     sh1.w -= kappa.im*(s1.z + s3.x);
     
     
     sh2.x = kappa.re*(s2.x + s3.z);
     sh2.x += kappa.im*(s2.y + s3.w);
     sh2.y = kappa.re*(s2.y + s3.w);
     sh2.y -= kappa.im*(s2.x + s3.z);
     
          
     sh2.z = kappa.re*(s2.z + s4.x);
     sh2.z += kappa.im*(s2.w + s4.y);
     sh2.w = kappa.re*(s2.w + s4.y);
     sh2.w -= kappa.im*(s2.z + s4.x);

     
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


__device__ void dev_su3MtV_kappaP2_minus_spintex2(dev_su3 M, int pos, dev_spinor * out, dev_complex kappa){

  
  dev_spinor sh0,sh1,sh2,zh0,zh1,zh2;

  dev_spinor s0, s1, s2, s3, s4, s5;
  s0 = tex1Dfetch(spin_tex_dn0,pos);
  s1 = tex1Dfetch(spin_tex_dn1,pos);
  s2 = tex1Dfetch(spin_tex_dn2,pos);
  s3 = tex1Dfetch(spin_tex_dn3,pos);
  s4 = tex1Dfetch(spin_tex_dn4,pos); 
  s5 = tex1Dfetch(spin_tex_dn5,pos); 
  
     sh0.x = kappa.re*(s0.x - s4.z);
     sh0.x += kappa.im*(s0.y - s4.w);
     sh0.y = kappa.re*(s0.y - s4.w);
     sh0.y -= kappa.im*(s0.x - s4.z);
     
   
     sh0.z = kappa.re*(s0.z - s5.x);
     sh0.z += kappa.im*(s0.w - s5.y); 
     sh0.w = kappa.re*(s0.w - s5.y);
     sh0.w -= kappa.im*(s0.z - s5.x);
     
       
     sh1.x = kappa.re*(s1.x - s5.z);
     sh1.x += kappa.im*(s1.y - s5.w);
     sh1.y = kappa.re*(s1.y - s5.w);
     sh1.y -= kappa.im*(s1.x - s5.z);
     


     
     
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


     sh1.z = kappa.re*(s1.z + s3.x);
     sh1.z += kappa.im*(s1.w + s3.y);
     sh1.w = kappa.re*(s1.w + s3.y);
     sh1.w -= kappa.im*(s1.z + s3.x);
     
     
     sh2.x = kappa.re*(s2.x + s3.z);
     sh2.x += kappa.im*(s2.y + s3.w);
     sh2.y = kappa.re*(s2.y + s3.w);
     sh2.y -= kappa.im*(s2.x + s3.z);
     
          
     sh2.z = kappa.re*(s2.z + s4.x);
     sh2.z += kappa.im*(s2.w + s4.y);
     sh2.w = kappa.re*(s2.w + s4.y);
     sh2.w -= kappa.im*(s2.z + s4.x);

     
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
__device__ void dev_su3MtV_kappaP1_plus(dev_su3 M, const dev_spinor * s, dev_spinor * out, dev_complex kappa){

  
  dev_spinor sh0,sh1,sh2,zh0,zh1,zh2;


     sh0.x = kappa.re*( (*(s+0*DEVOFF)).x - (*(s+4*DEVOFF)).w);
     sh0.x -= kappa.im*( (*(s+0*DEVOFF)).y + (*(s+4*DEVOFF)).z); 
     sh0.y = kappa.re*( (*(s+0*DEVOFF)).y + (*(s+4*DEVOFF)).z);
     sh0.y += kappa.im*( (*(s+0*DEVOFF)).x - (*(s+4*DEVOFF)).w);
     
     
     sh0.z = kappa.re*( (*(s+0*DEVOFF)).z - (*(s+5*DEVOFF)).y);
     sh0.z -= kappa.im*( (*(s+0*DEVOFF)).w + (*(s+5*DEVOFF)).x); 
     sh0.w = kappa.re*( (*(s+0*DEVOFF)).w + (*(s+5*DEVOFF)).x); 
     sh0.w += kappa.im*( (*(s+0*DEVOFF)).z - (*(s+5*DEVOFF)).y);
     

     sh1.x = kappa.re*((*(s+1*DEVOFF)).x - (*(s+5*DEVOFF)).w);
     sh1.x -= kappa.im*((*(s+1*DEVOFF)).y + (*(s+5*DEVOFF)).z);
     sh1.y = kappa.re*((*(s+1*DEVOFF)).y + (*(s+5*DEVOFF)).z); 
     sh1.y += kappa.im*((*(s+1*DEVOFF)).x - (*(s+5*DEVOFF)).w);
     
     
     

     
     
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


     sh1.z = kappa.re*((*(s+1*DEVOFF)).z - (*(s+3*DEVOFF)).y);
     sh1.z -= kappa.im*((*(s+1*DEVOFF)).w + (*(s+3*DEVOFF)).x);
     sh1.w = kappa.re*((*(s+1*DEVOFF)).w + (*(s+3*DEVOFF)).x); 
     sh1.w += kappa.im*((*(s+1*DEVOFF)).z - (*(s+3*DEVOFF)).y);
     
     
     sh2.x = kappa.re*((*(s+2*DEVOFF)).x - (*(s+3*DEVOFF)).w);
     sh2.x -= kappa.im*((*(s+2*DEVOFF)).y + (*(s+3*DEVOFF)).z);
     sh2.y = kappa.re*((*(s+2*DEVOFF)).y + (*(s+3*DEVOFF)).z);
     sh2.y += kappa.im*((*(s+2*DEVOFF)).x - (*(s+3*DEVOFF)).w);
     
     
     sh2.z = kappa.re*((*(s+2*DEVOFF)).z - (*(s+4*DEVOFF)).y);
     sh2.z -= kappa.im*((*(s+2*DEVOFF)).w + (*(s+4*DEVOFF)).x);
     sh2.w = kappa.re*((*(s+2*DEVOFF)).w + (*(s+4*DEVOFF)).x);
     sh2.w += kappa.im*((*(s+2*DEVOFF)).z - (*(s+4*DEVOFF)).y);

     
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


__device__ void dev_su3MtV_kappaP1_plus_spintex(dev_su3 M, int pos, dev_spinor * out, dev_complex kappa){

  
  dev_spinor sh0,sh1,sh2,zh0,zh1,zh2;

  dev_spinor s0, s1, s2, s3, s4, s5;
  s0 = tex1Dfetch(spin_tex0,pos);
  s1 = tex1Dfetch(spin_tex1,pos);
  s2 = tex1Dfetch(spin_tex2,pos);
  s3 = tex1Dfetch(spin_tex3,pos);
  s4 = tex1Dfetch(spin_tex4,pos); 
  s5 = tex1Dfetch(spin_tex5,pos); 
  
     sh0.x = kappa.re*(s0.x - s4.w);
     sh0.x -= kappa.im*(s0.y + s4.z);
     sh0.y = kappa.re*(s0.y + s4.z);
     sh0.y += kappa.im*(s0.x - s4.w);
     
     
     sh0.z = kappa.re*(s0.z - s5.y);
     sh0.z -= kappa.im*(s0.w + s5.x);
     sh0.w = kappa.re*(s0.w + s5.x);    
     sh0.w += kappa.im*(s0.z - s5.y);
     

     sh1.x = kappa.re*(s1.x - s5.w);
     sh1.x -= kappa.im*(s1.y + s5.z); 
     sh1.y = kappa.re*(s1.y + s5.z); 
     sh1.y += kappa.im*(s1.x - s5.w);
     
     
     
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


     sh1.z = kappa.re*(s1.z - s3.y);
     sh1.z -= kappa.im*(s1.w + s3.x); 
     sh1.w = kappa.re*(s1.w + s3.x); 
     sh1.w += kappa.im*(s1.z - s3.y);
     
     
     sh2.x = kappa.re*(s2.x - s3.w);
     sh2.x -= kappa.im*(s2.y + s3.z);
     sh2.y = kappa.re*(s2.y + s3.z);
     sh2.y += kappa.im*(s2.x - s3.w);
     
     
     sh2.z = kappa.re*(s2.z - s4.y);
     sh2.z -= kappa.im*(s2.w + s4.x);
     sh2.w = kappa.re*(s2.w + s4.x);
     sh2.w += kappa.im*(s2.z - s4.y);

     
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



__device__ void dev_su3MtV_kappaP1_plus_spintex2(dev_su3 M, int pos, dev_spinor * out, dev_complex kappa){

  
  dev_spinor sh0,sh1,sh2,zh0,zh1,zh2;

  dev_spinor s0, s1, s2, s3, s4, s5;
  s0 = tex1Dfetch(spin_tex_dn0,pos);
  s1 = tex1Dfetch(spin_tex_dn1,pos);
  s2 = tex1Dfetch(spin_tex_dn2,pos);
  s3 = tex1Dfetch(spin_tex_dn3,pos);
  s4 = tex1Dfetch(spin_tex_dn4,pos); 
  s5 = tex1Dfetch(spin_tex_dn5,pos); 
  
     sh0.x = kappa.re*(s0.x - s4.w);
     sh0.x -= kappa.im*(s0.y + s4.z);
     sh0.y = kappa.re*(s0.y + s4.z);
     sh0.y += kappa.im*(s0.x - s4.w);
     
     
     sh0.z = kappa.re*(s0.z - s5.y);
     sh0.z -= kappa.im*(s0.w + s5.x);
     sh0.w = kappa.re*(s0.w + s5.x);    
     sh0.w += kappa.im*(s0.z - s5.y);
     

     sh1.x = kappa.re*(s1.x - s5.w);
     sh1.x -= kappa.im*(s1.y + s5.z); 
     sh1.y = kappa.re*(s1.y + s5.z); 
     sh1.y += kappa.im*(s1.x - s5.w);
     
     
     
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


     sh1.z = kappa.re*(s1.z - s3.y);
     sh1.z -= kappa.im*(s1.w + s3.x); 
     sh1.w = kappa.re*(s1.w + s3.x); 
     sh1.w += kappa.im*(s1.z - s3.y);
     
     
     sh2.x = kappa.re*(s2.x - s3.w);
     sh2.x -= kappa.im*(s2.y + s3.z);
     sh2.y = kappa.re*(s2.y + s3.z);
     sh2.y += kappa.im*(s2.x - s3.w);
     
     
     sh2.z = kappa.re*(s2.z - s4.y);
     sh2.z -= kappa.im*(s2.w + s4.x);
     sh2.w = kappa.re*(s2.w + s4.x);
     sh2.w += kappa.im*(s2.z - s4.y);

     
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
//-cconj(kappa)(r + gamma_mu) kappa reell !!!!
__device__ void dev_su3MtV_kappaP1_minus(dev_su3 M, const dev_spinor * s, dev_spinor * out, dev_complex kappa){

  
  dev_spinor sh0,sh1,sh2,zh0,zh1,zh2;


     sh0.x = kappa.re*( (*(s+0*DEVOFF)).x + (*(s+4*DEVOFF)).w);
     sh0.x += kappa.im*( (*(s+0*DEVOFF)).y - (*(s+4*DEVOFF)).z);
     sh0.y = kappa.re*( (*(s+0*DEVOFF)).y - (*(s+4*DEVOFF)).z);
     sh0.y -= kappa.im*( (*(s+0*DEVOFF)).x + (*(s+4*DEVOFF)).w);
     
     
     sh0.z = kappa.re*( (*(s+0*DEVOFF)).z + (*(s+5*DEVOFF)).y);
     sh0.z += kappa.im*( (*(s+0*DEVOFF)).w - (*(s+5*DEVOFF)).x);  
     sh0.w = kappa.re*( (*(s+0*DEVOFF)).w - (*(s+5*DEVOFF)).x);    
     sh0.w -= kappa.im*( (*(s+0*DEVOFF)).z + (*(s+5*DEVOFF)).y);
     

     sh1.x = kappa.re*((*(s+1*DEVOFF)).x + (*(s+5*DEVOFF)).w);
     sh1.x += kappa.im*((*(s+1*DEVOFF)).y - (*(s+5*DEVOFF)).z);
     sh1.y = kappa.re*((*(s+1*DEVOFF)).y - (*(s+5*DEVOFF)).z); 
     sh1.y -= kappa.im*((*(s+1*DEVOFF)).x + (*(s+5*DEVOFF)).w);
     

     
     
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


     sh1.z = kappa.re*((*(s+1*DEVOFF)).z + (*(s+3*DEVOFF)).y);
     sh1.z += kappa.im*((*(s+1*DEVOFF)).w - (*(s+3*DEVOFF)).x);
     sh1.w = kappa.re*((*(s+1*DEVOFF)).w - (*(s+3*DEVOFF)).x);
     sh1.w -= kappa.im*((*(s+1*DEVOFF)).z + (*(s+3*DEVOFF)).y);
     
     
     sh2.x = kappa.re*((*(s+2*DEVOFF)).x + (*(s+3*DEVOFF)).w);
     sh2.x += kappa.im*((*(s+2*DEVOFF)).y - (*(s+3*DEVOFF)).z);
     sh2.y = kappa.re*((*(s+2*DEVOFF)).y - (*(s+3*DEVOFF)).z);
     sh2.y -= kappa.im*((*(s+2*DEVOFF)).x + (*(s+3*DEVOFF)).w);
     
     
     sh2.z = kappa.re*((*(s+2*DEVOFF)).z + (*(s+4*DEVOFF)).y);
     sh2.z += kappa.im*((*(s+2*DEVOFF)).w - (*(s+4*DEVOFF)).x);
     sh2.w = kappa.re*((*(s+2*DEVOFF)).w - (*(s+4*DEVOFF)).x);
     sh2.w -= kappa.im*((*(s+2*DEVOFF)).z + (*(s+4*DEVOFF)).y);

     

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


__device__ void dev_su3MtV_kappaP1_minus_spintex(dev_su3 M, int pos, dev_spinor * out, dev_complex kappa){

  
  dev_spinor sh0,sh1,sh2,zh0,zh1,zh2;


  dev_spinor s0, s1, s2, s3, s4, s5;
  s0 = tex1Dfetch(spin_tex0,pos);
  s1 = tex1Dfetch(spin_tex1,pos);
  s2 = tex1Dfetch(spin_tex2,pos);
  s3 = tex1Dfetch(spin_tex3,pos);
  s4 = tex1Dfetch(spin_tex4,pos); 
  s5 = tex1Dfetch(spin_tex5,pos); 
  
     sh0.x = kappa.re*(s0.x + s4.w);
     sh0.x += kappa.im*(s0.y - s4.z);
     sh0.y = kappa.re*(s0.y - s4.z);
     sh0.y -= kappa.im*(s0.x + s4.w);
     
     
     sh0.z = kappa.re*(s0.z + s5.y);
     sh0.z += kappa.im*(s0.w - s5.x); 
     sh0.w = kappa.re*(s0.w - s5.x);
     sh0.w -= kappa.im*(s0.z + s5.y);
     

     sh1.x = kappa.re*(s1.x + s5.w);
     sh1.x += kappa.im*(s1.y - s5.z);
     sh1.y = kappa.re*(s1.y - s5.z); 
     sh1.y -= kappa.im*(s1.x + s5.w);
     

     
     
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


     sh1.z = kappa.re*(s1.z + s3.y);
     sh1.z += kappa.im*(s1.w - s3.x);
     sh1.w = kappa.re*(s1.w - s3.x); 
     sh1.w -= kappa.im*(s1.z + s3.y);
     
     
     sh2.x = kappa.re*(s2.x + s3.w);
     sh2.x += kappa.im*(s2.y - s3.z);
     sh2.y = kappa.re*(s2.y - s3.z);
     sh2.y -= kappa.im*(s2.x + s3.w);
     
     
     sh2.z = kappa.re*(s2.z + s4.y);
     sh2.z += kappa.im*(s2.w - s4.x);
     sh2.w = kappa.re*(s2.w - s4.x);
     sh2.w -= kappa.im*(s2.z + s4.y);


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


__device__ void dev_su3MtV_kappaP1_minus_spintex2(dev_su3 M, int pos, dev_spinor * out, dev_complex kappa){

  
  dev_spinor sh0,sh1,sh2,zh0,zh1,zh2;


  dev_spinor s0, s1, s2, s3, s4, s5;
  s0 = tex1Dfetch(spin_tex_dn0,pos);
  s1 = tex1Dfetch(spin_tex_dn1,pos);
  s2 = tex1Dfetch(spin_tex_dn2,pos);
  s3 = tex1Dfetch(spin_tex_dn3,pos);
  s4 = tex1Dfetch(spin_tex_dn4,pos); 
  s5 = tex1Dfetch(spin_tex_dn5,pos); 
  
     sh0.x = kappa.re*(s0.x + s4.w);
     sh0.x += kappa.im*(s0.y - s4.z);
     sh0.y = kappa.re*(s0.y - s4.z);
     sh0.y -= kappa.im*(s0.x + s4.w);
     
     
     sh0.z = kappa.re*(s0.z + s5.y);
     sh0.z += kappa.im*(s0.w - s5.x); 
     sh0.w = kappa.re*(s0.w - s5.x);
     sh0.w -= kappa.im*(s0.z + s5.y);
     

     sh1.x = kappa.re*(s1.x + s5.w);
     sh1.x += kappa.im*(s1.y - s5.z);
     sh1.y = kappa.re*(s1.y - s5.z); 
     sh1.y -= kappa.im*(s1.x + s5.w);
     

     
     
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


     sh1.z = kappa.re*(s1.z + s3.y);
     sh1.z += kappa.im*(s1.w - s3.x);
     sh1.w = kappa.re*(s1.w - s3.x); 
     sh1.w -= kappa.im*(s1.z + s3.y);
     
     
     sh2.x = kappa.re*(s2.x + s3.w);
     sh2.x += kappa.im*(s2.y - s3.z);
     sh2.y = kappa.re*(s2.y - s3.z);
     sh2.y -= kappa.im*(s2.x + s3.w);
     
     
     sh2.z = kappa.re*(s2.z + s4.y);
     sh2.z += kappa.im*(s2.w - s4.x);
     sh2.w = kappa.re*(s2.w - s4.x);
     sh2.w -= kappa.im*(s2.z + s4.y);


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



//multipliziert su3-Matrix mal Spinor im Dirac-Raum
//we do not load the lower spinor components as we are working in the relativistic basis
__device__ void dev_su3MtV_rel_up(dev_su3 M, const dev_spinor * s, dev_spinor * out){

float help;

help =  M[0][0].re*(*(s+0*DEVOFF)).x;
help -= M[0][0].im*(*(s+0*DEVOFF)).y;
help += M[0][1].re*(*(s+0*DEVOFF)).z;
help -= M[0][1].im*(*(s+0*DEVOFF)).w;
help += M[0][2].re*(*(s+1*DEVOFF)).x;
help -= M[0][2].im*(*(s+1*DEVOFF)).y;
(*(out+0)).x = help;		 		 
help =  M[0][0].re*(*(s+0*DEVOFF)).y;
help += M[0][0].im*(*(s+0*DEVOFF)).x; 
help += M[0][1].re*(*(s+0*DEVOFF)).w;
help += M[0][1].im*(*(s+0*DEVOFF)).z;
help += M[0][2].re*(*(s+1*DEVOFF)).y;
help += M[0][2].im*(*(s+1*DEVOFF)).x;
(*(out+0)).y = help;

help =  M[1][0].re*(*(s+0*DEVOFF)).x;
help -= M[1][0].im*(*(s+0*DEVOFF)).y;
help += M[1][1].re*(*(s+0*DEVOFF)).z;
help -= M[1][1].im*(*(s+0*DEVOFF)).w;
help += M[1][2].re*(*(s+1*DEVOFF)).x;
help -= M[1][2].im*(*(s+1*DEVOFF)).y;
(*(out+0)).z = help; 		  		  
help =  M[1][0].re*(*(s+0*DEVOFF)).y;
help += M[1][0].im*(*(s+0*DEVOFF)).x;
help += M[1][1].re*(*(s+0*DEVOFF)).w;
help += M[1][1].im*(*(s+0*DEVOFF)).z;
help += M[1][2].re*(*(s+1*DEVOFF)).y;
help += M[1][2].im*(*(s+1*DEVOFF)).x;
(*(out+0)).w = help;

help =  M[2][0].re*(*(s+0*DEVOFF)).x;
help -= M[2][0].im*(*(s+0*DEVOFF)).y;
help += M[2][1].re*(*(s+0*DEVOFF)).z;
help -= M[2][1].im*(*(s+0*DEVOFF)).w;
help += M[2][2].re*(*(s+1*DEVOFF)).x;
help -= M[2][2].im*(*(s+1*DEVOFF)).y;
(*(out+1)).x = help;
help =  M[2][0].re*(*(s+0*DEVOFF)).y;
help += M[2][0].im*(*(s+0*DEVOFF)).x;
help += M[2][1].re*(*(s+0*DEVOFF)).w;
help += M[2][1].im*(*(s+0*DEVOFF)).z;
help += M[2][2].re*(*(s+1*DEVOFF)).y;
help += M[2][2].im*(*(s+1*DEVOFF)).x;
(*(out+1)).y = help;


help =  M[0][0].re*(*(s+1*DEVOFF)).z;
help -= M[0][0].im*(*(s+1*DEVOFF)).w;
help += M[0][1].re*(*(s+2*DEVOFF)).x;
help -= M[0][1].im*(*(s+2*DEVOFF)).y;
help += M[0][2].re*(*(s+2*DEVOFF)).z;
help -= M[0][2].im*(*(s+2*DEVOFF)).w;
(*(out+1)).z = help;		 
help =  M[0][0].re*(*(s+1*DEVOFF)).w;
help += M[0][0].im*(*(s+1*DEVOFF)).z;
help += M[0][1].re*(*(s+2*DEVOFF)).y;
help += M[0][1].im*(*(s+2*DEVOFF)).x;
help += M[0][2].re*(*(s+2*DEVOFF)).w;
help += M[0][2].im*(*(s+2*DEVOFF)).z;
(*(out+1)).w = help;


help =  M[1][0].re*(*(s+1*DEVOFF)).z;
help -= M[1][0].im*(*(s+1*DEVOFF)).w;
help += M[1][1].re*(*(s+2*DEVOFF)).x;
help -= M[1][1].im*(*(s+2*DEVOFF)).y;
help += M[1][2].re*(*(s+2*DEVOFF)).z;
help -= M[1][2].im*(*(s+2*DEVOFF)).w;
(*(out+2)).x = help;		 		 
help =  M[1][0].re*(*(s+1*DEVOFF)).w;
help += M[1][0].im*(*(s+1*DEVOFF)).z;
help += M[1][1].re*(*(s+2*DEVOFF)).y;
help += M[1][1].im*(*(s+2*DEVOFF)).x;
help += M[1][2].re*(*(s+2*DEVOFF)).w;
help += M[1][2].im*(*(s+2*DEVOFF)).z;
(*(out+2)).y = help;


help =  M[2][0].re*(*(s+1*DEVOFF)).z;
help -= M[2][0].im*(*(s+1*DEVOFF)).w;
help += M[2][1].re*(*(s+2*DEVOFF)).x;
help -= M[2][1].im*(*(s+2*DEVOFF)).y;
help += M[2][2].re*(*(s+2*DEVOFF)).z;
help -= M[2][2].im*(*(s+2*DEVOFF)).w;
(*(out+2)).z = help;		 		 
help =  M[2][0].re*(*(s+1*DEVOFF)).w;
help += M[2][0].im*(*(s+1*DEVOFF)).z;
help += M[2][1].re*(*(s+2*DEVOFF)).y;
help += M[2][1].im*(*(s+2*DEVOFF)).x;
help += M[2][2].re*(*(s+2*DEVOFF)).w;
help += M[2][2].im*(*(s+2*DEVOFF)).z;
(*(out+2)).w = help;

(*(out+3)).x = 0.0f;
(*(out+3)).y = 0.0f;
(*(out+3)).z = 0.0f;
(*(out+3)).w = 0.0f;


(*(out+4)).x = 0.0f;
(*(out+4)).y = 0.0f;
(*(out+4)).z = 0.0f;
(*(out+4)).w = 0.0f;


(*(out+5)).x = 0.0f;
(*(out+5)).y = 0.0f;
(*(out+5)).z = 0.0f;
(*(out+5)).w = 0.0f;
}



//multipliziert su3-Matrix mal Spinor im Dirac-Raum
//we do not load the upper spinor components as we are working in the relativistic basis
//code in su3_MtV.txt -- generated with codegen
__device__ void dev_su3MtV_rel_down(dev_su3 M, const dev_spinor * s, dev_spinor * out){

float help;
(*(out+0)).x = 0.0f;
(*(out+0)).y = 0.0f;
(*(out+0)).z = 0.0f;
(*(out+0)).w = 0.0f; 


(*(out+1)).x = 0.0f;
(*(out+1)).y = 0.0f;
(*(out+1)).z = 0.0f;
(*(out+1)).w = 0.0f;


(*(out+2)).x = 0.0f;
(*(out+2)).y = 0.0f;
(*(out+2)).z = 0.0f;
(*(out+2)).w = 0.0f;


help =  M[0][0].re*(*(s+3*DEVOFF)).x;
help -= M[0][0].im*(*(s+3*DEVOFF)).y;
help += M[0][1].re*(*(s+3*DEVOFF)).z;
help -= M[0][1].im*(*(s+3*DEVOFF)).w;
help += M[0][2].re*(*(s+4*DEVOFF)).x;
help -= M[0][2].im*(*(s+4*DEVOFF)).y;
(*(out+3)).x = help;
help =  M[0][0].re*(*(s+3*DEVOFF)).y;
help += M[0][0].im*(*(s+3*DEVOFF)).x;
help += M[0][1].re*(*(s+3*DEVOFF)).w;
help += M[0][1].im*(*(s+3*DEVOFF)).z;
help += M[0][2].re*(*(s+4*DEVOFF)).y;
help += M[0][2].im*(*(s+4*DEVOFF)).x;
(*(out+3)).y = help;


help =  M[1][0].re*(*(s+3*DEVOFF)).x;
help -= M[1][0].im*(*(s+3*DEVOFF)).y;
help += M[1][1].re*(*(s+3*DEVOFF)).z;
help -= M[1][1].im*(*(s+3*DEVOFF)).w;
help += M[1][2].re*(*(s+4*DEVOFF)).x;
help -= M[1][2].im*(*(s+4*DEVOFF)).y;
(*(out+3)).z = help;		 
help =  M[1][0].re*(*(s+3*DEVOFF)).y;
help += M[1][0].im*(*(s+3*DEVOFF)).x;
help += M[1][1].re*(*(s+3*DEVOFF)).w;
help += M[1][1].im*(*(s+3*DEVOFF)).z;
help += M[1][2].re*(*(s+4*DEVOFF)).y;
help += M[1][2].im*(*(s+4*DEVOFF)).x;
(*(out+3)).w = help;


help =  M[2][0].re*(*(s+3*DEVOFF)).x;
help -= M[2][0].im*(*(s+3*DEVOFF)).y;
help += M[2][1].re*(*(s+3*DEVOFF)).z;
help -= M[2][1].im*(*(s+3*DEVOFF)).w;
help += M[2][2].re*(*(s+4*DEVOFF)).x;
help -= M[2][2].im*(*(s+4*DEVOFF)).y;
(*(out+4)).x = help;		  
help =  M[2][0].re*(*(s+3*DEVOFF)).y;
help += M[2][0].im*(*(s+3*DEVOFF)).x;
help += M[2][1].re*(*(s+3*DEVOFF)).w;
help += M[2][1].im*(*(s+3*DEVOFF)).z;
help += M[2][2].re*(*(s+4*DEVOFF)).y;
help += M[2][2].im*(*(s+4*DEVOFF)).x;
(*(out+4)).y = help;


help =  M[0][0].re*(*(s+4*DEVOFF)).z;
help -= M[0][0].im*(*(s+4*DEVOFF)).w;
help += M[0][1].re*(*(s+5*DEVOFF)).x;
help -= M[0][1].im*(*(s+5*DEVOFF)).y;
help += M[0][2].re*(*(s+5*DEVOFF)).z;
help -= M[0][2].im*(*(s+5*DEVOFF)).w;
(*(out+4)).z = help;		 		 
help =  M[0][0].re*(*(s+4*DEVOFF)).w;
help += M[0][0].im*(*(s+4*DEVOFF)).z;
help += M[0][1].re*(*(s+5*DEVOFF)).y;
help += M[0][1].im*(*(s+5*DEVOFF)).x;
help += M[0][2].re*(*(s+5*DEVOFF)).w;
help += M[0][2].im*(*(s+5*DEVOFF)).z;
(*(out+4)).w = help;


help =  M[1][0].re*(*(s+4*DEVOFF)).z;
help -= M[1][0].im*(*(s+4*DEVOFF)).w;
help += M[1][1].re*(*(s+5*DEVOFF)).x;
help -= M[1][1].im*(*(s+5*DEVOFF)).y;
help += M[1][2].re*(*(s+5*DEVOFF)).z;
help -= M[1][2].im*(*(s+5*DEVOFF)).w;
(*(out+5)).x = help;
help =  M[1][0].re*(*(s+4*DEVOFF)).w;
help += M[1][0].im*(*(s+4*DEVOFF)).z;
help += M[1][1].re*(*(s+5*DEVOFF)).y;
help += M[1][1].im*(*(s+5*DEVOFF)).x;
help += M[1][2].re*(*(s+5*DEVOFF)).w;
help += M[1][2].im*(*(s+5*DEVOFF)).z;
(*(out+5)).y = help;


help =  M[2][0].re*(*(s+4*DEVOFF)).z;
help -= M[2][0].im*(*(s+4*DEVOFF)).w;
help += M[2][1].re*(*(s+5*DEVOFF)).x;
help -= M[2][1].im*(*(s+5*DEVOFF)).y;
help += M[2][2].re*(*(s+5*DEVOFF)).z;
help -= M[2][2].im*(*(s+5*DEVOFF)).w;
(*(out+5)).z = help;
help =  M[2][0].re*(*(s+4*DEVOFF)).w;
help += M[2][0].im*(*(s+4*DEVOFF)).z;
help += M[2][1].re*(*(s+5*DEVOFF)).y;
help += M[2][1].im*(*(s+5*DEVOFF)).x;
help += M[2][2].re*(*(s+5*DEVOFF)).w;
help += M[2][2].im*(*(s+5*DEVOFF)).z;
(*(out+5)).w = help;	

}















#ifdef HALF

//multipliziert su3-Matrix mal Spinor im Dirac-Raum
//code in su3_MtV.txt -- generated with codegen
__device__ void dev_su3MtV_half(dev_su3 M, const dev_spinor_half * s, const float * s_norm, dev_spinor * out){
float norm = * s_norm;

(*(out+0)).x = ( M[0][0].re*half2fl((*(s+0*DEVOFF)).x,norm) -
                 M[0][0].im*half2fl((*(s+0*DEVOFF)).y,norm) ) 
             + ( M[0][1].re*half2fl((*(s+0*DEVOFF)).z,norm) - 
                 M[0][1].im*half2fl((*(s+0*DEVOFF)).w,norm) ) 
             + ( M[0][2].re*half2fl((*(s+1*DEVOFF)).x,norm) -
                 M[0][2].im*half2fl((*(s+1*DEVOFF)).y,norm) );
(*(out+0)).y = ( M[0][0].re*half2fl((*(s+0*DEVOFF)).y,norm) +
                 M[0][0].im*half2fl((*(s+0*DEVOFF)).x,norm) ) 
             + ( M[0][1].re*half2fl((*(s+0*DEVOFF)).w,norm) +
                 M[0][1].im*half2fl((*(s+0*DEVOFF)).z,norm) ) 
             + ( M[0][2].re*half2fl((*(s+1*DEVOFF)).y,norm) +
                 M[0][2].im*half2fl((*(s+1*DEVOFF)).x,norm) );


(*(out+0)).z =  ( M[1][0].re*half2fl((*(s+0*DEVOFF)).x,norm) -
                  M[1][0].im*half2fl((*(s+0*DEVOFF)).y,norm) ) 
             +  ( M[1][1].re*half2fl((*(s+0*DEVOFF)).z,norm) - 
                  M[1][1].im*half2fl((*(s+0*DEVOFF)).w,norm) ) 
             +  ( M[1][2].re*half2fl((*(s+1*DEVOFF)).x,norm) -
                  M[1][2].im*half2fl((*(s+1*DEVOFF)).y,norm) );
(*(out+0)).w =  ( M[1][0].re*half2fl((*(s+0*DEVOFF)).y,norm) +
                  M[1][0].im*half2fl((*(s+0*DEVOFF)).x,norm) ) 
             +  ( M[1][1].re*half2fl((*(s+0*DEVOFF)).w,norm) +
                  M[1][1].im*half2fl((*(s+0*DEVOFF)).z,norm) ) 
             +  ( M[1][2].re*half2fl((*(s+1*DEVOFF)).y,norm) + 
                  M[1][2].im*half2fl((*(s+1*DEVOFF)).x,norm) );


(*(out+1)).x = ( M[2][0].re*half2fl((*(s+0*DEVOFF)).x, norm) -
                 M[2][0].im*half2fl((*(s+0*DEVOFF)).y, norm) ) 
             + ( M[2][1].re*half2fl((*(s+0*DEVOFF)).z, norm) -
                 M[2][1].im*half2fl((*(s+0*DEVOFF)).w, norm) ) 
             + ( M[2][2].re*half2fl((*(s+1*DEVOFF)).x, norm) -
                 M[2][2].im*half2fl((*(s+1*DEVOFF)).y, norm) );
(*(out+1)).y = ( M[2][0].re*half2fl((*(s+0*DEVOFF)).y, norm) +
                 M[2][0].im*half2fl((*(s+0*DEVOFF)).x, norm) ) 
             + ( M[2][1].re*half2fl((*(s+0*DEVOFF)).w, norm) +
                 M[2][1].im*half2fl((*(s+0*DEVOFF)).z, norm) ) 
             + ( M[2][2].re*half2fl((*(s+1*DEVOFF)).y, norm) +
                 M[2][2].im*half2fl((*(s+1*DEVOFF)).x, norm) );


(*(out+1)).z = ( M[0][0].re*half2fl((*(s+1*DEVOFF)).z, norm) -
                 M[0][0].im*half2fl((*(s+1*DEVOFF)).w, norm) ) 
             + ( M[0][1].re*half2fl((*(s+2*DEVOFF)).x, norm) -
                 M[0][1].im*half2fl((*(s+2*DEVOFF)).y, norm) ) 
             + ( M[0][2].re*half2fl((*(s+2*DEVOFF)).z, norm) -
                 M[0][2].im*half2fl((*(s+2*DEVOFF)).w, norm) );
(*(out+1)).w = ( M[0][0].re*half2fl((*(s+1*DEVOFF)).w, norm) +
                 M[0][0].im*half2fl((*(s+1*DEVOFF)).z, norm) ) 
             + ( M[0][1].re*half2fl((*(s+2*DEVOFF)).y, norm) +
                 M[0][1].im*half2fl((*(s+2*DEVOFF)).x, norm) ) 
             + ( M[0][2].re*half2fl((*(s+2*DEVOFF)).w, norm) +
                 M[0][2].im*half2fl((*(s+2*DEVOFF)).z, norm) );


(*(out+2)).x = ( M[1][0].re*half2fl((*(s+1*DEVOFF)).z, norm) -
                 M[1][0].im*half2fl((*(s+1*DEVOFF)).w, norm) ) 
             + ( M[1][1].re*half2fl((*(s+2*DEVOFF)).x, norm) -
                 M[1][1].im*half2fl((*(s+2*DEVOFF)).y, norm) ) 
             + ( M[1][2].re*half2fl((*(s+2*DEVOFF)).z, norm) -
                 M[1][2].im*half2fl((*(s+2*DEVOFF)).w, norm) );
(*(out+2)).y = ( M[1][0].re*half2fl((*(s+1*DEVOFF)).w, norm) +
                 M[1][0].im*half2fl((*(s+1*DEVOFF)).z, norm) ) 
             + ( M[1][1].re*half2fl((*(s+2*DEVOFF)).y, norm) +
                 M[1][1].im*half2fl((*(s+2*DEVOFF)).x, norm) ) 
             + ( M[1][2].re*half2fl((*(s+2*DEVOFF)).w, norm) +
                 M[1][2].im*half2fl((*(s+2*DEVOFF)).z, norm) );


(*(out+2)).z = ( M[2][0].re*half2fl((*(s+1*DEVOFF)).z, norm) -
                 M[2][0].im*half2fl((*(s+1*DEVOFF)).w, norm) ) 
             + ( M[2][1].re*half2fl((*(s+2*DEVOFF)).x, norm) -
                 M[2][1].im*half2fl((*(s+2*DEVOFF)).y, norm) ) 
             + ( M[2][2].re*half2fl((*(s+2*DEVOFF)).z, norm) -
                 M[2][2].im*half2fl((*(s+2*DEVOFF)).w, norm) );
(*(out+2)).w = ( M[2][0].re*half2fl((*(s+1*DEVOFF)).w, norm) +
                 M[2][0].im*half2fl((*(s+1*DEVOFF)).z, norm) ) 
             + ( M[2][1].re*half2fl((*(s+2*DEVOFF)).y, norm) +
                 M[2][1].im*half2fl((*(s+2*DEVOFF)).x, norm) ) 
             + ( M[2][2].re*half2fl((*(s+2*DEVOFF)).w, norm) +
                 M[2][2].im*half2fl((*(s+2*DEVOFF)).z, norm) );


(*(out+3)).x = ( M[0][0].re*half2fl((*(s+3*DEVOFF)).x, norm) -
                 M[0][0].im*half2fl((*(s+3*DEVOFF)).y, norm) ) 
             + ( M[0][1].re*half2fl((*(s+3*DEVOFF)).z, norm) -
                 M[0][1].im*half2fl((*(s+3*DEVOFF)).w, norm) ) 
             + ( M[0][2].re*half2fl((*(s+4*DEVOFF)).x, norm) -
                 M[0][2].im*half2fl((*(s+4*DEVOFF)).y, norm) );
(*(out+3)).y = ( M[0][0].re*half2fl((*(s+3*DEVOFF)).y, norm) +
                 M[0][0].im*half2fl((*(s+3*DEVOFF)).x, norm) ) 
             + ( M[0][1].re*half2fl((*(s+3*DEVOFF)).w, norm) +
                 M[0][1].im*half2fl((*(s+3*DEVOFF)).z, norm) ) 
             + ( M[0][2].re*half2fl((*(s+4*DEVOFF)).y, norm) +
                 M[0][2].im*half2fl((*(s+4*DEVOFF)).x, norm) );


(*(out+3)).z = ( M[1][0].re*half2fl((*(s+3*DEVOFF)).x, norm) -
                 M[1][0].im*half2fl((*(s+3*DEVOFF)).y, norm) ) 
             + ( M[1][1].re*half2fl((*(s+3*DEVOFF)).z, norm) -
                 M[1][1].im*half2fl((*(s+3*DEVOFF)).w, norm) ) 
             + ( M[1][2].re*half2fl((*(s+4*DEVOFF)).x, norm) -
                 M[1][2].im*half2fl((*(s+4*DEVOFF)).y, norm) );
(*(out+3)).w = ( M[1][0].re*half2fl((*(s+3*DEVOFF)).y, norm) +
                 M[1][0].im*half2fl((*(s+3*DEVOFF)).x, norm) ) 
             + ( M[1][1].re*half2fl((*(s+3*DEVOFF)).w, norm) +
                 M[1][1].im*half2fl((*(s+3*DEVOFF)).z, norm) ) 
             + ( M[1][2].re*half2fl((*(s+4*DEVOFF)).y, norm) +
                 M[1][2].im*half2fl((*(s+4*DEVOFF)).x, norm) );


(*(out+4)).x = ( M[2][0].re*half2fl((*(s+3*DEVOFF)).x, norm) -
                 M[2][0].im*half2fl((*(s+3*DEVOFF)).y, norm) ) 
             + ( M[2][1].re*half2fl((*(s+3*DEVOFF)).z, norm) -
                 M[2][1].im*half2fl((*(s+3*DEVOFF)).w, norm) ) 
             + ( M[2][2].re*half2fl((*(s+4*DEVOFF)).x, norm) -
                 M[2][2].im*half2fl((*(s+4*DEVOFF)).y, norm) );
(*(out+4)).y = ( M[2][0].re*half2fl((*(s+3*DEVOFF)).y, norm) +
                 M[2][0].im*half2fl((*(s+3*DEVOFF)).x, norm) ) 
             + ( M[2][1].re*half2fl((*(s+3*DEVOFF)).w, norm) +
                 M[2][1].im*half2fl((*(s+3*DEVOFF)).z, norm) ) 
             + ( M[2][2].re*half2fl((*(s+4*DEVOFF)).y, norm) +
                 M[2][2].im*half2fl((*(s+4*DEVOFF)).x, norm) );


(*(out+4)).z = ( M[0][0].re*half2fl((*(s+4*DEVOFF)).z, norm) -
                 M[0][0].im*half2fl((*(s+4*DEVOFF)).w, norm) ) 
             + ( M[0][1].re*half2fl((*(s+5*DEVOFF)).x, norm) - 
                 M[0][1].im*half2fl((*(s+5*DEVOFF)).y, norm) ) 
             + ( M[0][2].re*half2fl((*(s+5*DEVOFF)).z, norm) -
                 M[0][2].im*half2fl((*(s+5*DEVOFF)).w, norm) );
(*(out+4)).w = ( M[0][0].re*half2fl((*(s+4*DEVOFF)).w, norm) +
                 M[0][0].im*half2fl((*(s+4*DEVOFF)).z, norm) ) 
             + ( M[0][1].re*half2fl((*(s+5*DEVOFF)).y, norm) +
                 M[0][1].im*half2fl((*(s+5*DEVOFF)).x, norm) ) 
             + ( M[0][2].re*half2fl((*(s+5*DEVOFF)).w, norm) +
                 M[0][2].im*half2fl((*(s+5*DEVOFF)).z, norm) );


(*(out+5)).x = ( M[1][0].re*half2fl((*(s+4*DEVOFF)).z, norm) -
                 M[1][0].im*half2fl((*(s+4*DEVOFF)).w, norm) ) 
             + ( M[1][1].re*half2fl((*(s+5*DEVOFF)).x, norm) -
                 M[1][1].im*half2fl((*(s+5*DEVOFF)).y, norm) ) 
             + ( M[1][2].re*half2fl((*(s+5*DEVOFF)).z, norm) -
                 M[1][2].im*half2fl((*(s+5*DEVOFF)).w, norm) );
(*(out+5)).y = ( M[1][0].re*half2fl((*(s+4*DEVOFF)).w, norm) +
                 M[1][0].im*half2fl((*(s+4*DEVOFF)).z, norm) ) 
             + ( M[1][1].re*half2fl((*(s+5*DEVOFF)).y, norm) +
                 M[1][1].im*half2fl((*(s+5*DEVOFF)).x, norm) ) 
             + ( M[1][2].re*half2fl((*(s+5*DEVOFF)).w, norm) +
                 M[1][2].im*half2fl((*(s+5*DEVOFF)).z, norm) );


(*(out+5)).z = ( M[2][0].re*half2fl((*(s+4*DEVOFF)).z, norm) - 
                 M[2][0].im*half2fl((*(s+4*DEVOFF)).w, norm) ) 
             + ( M[2][1].re*half2fl((*(s+5*DEVOFF)).x, norm) -
                 M[2][1].im*half2fl((*(s+5*DEVOFF)).y, norm) ) 
             + ( M[2][2].re*half2fl((*(s+5*DEVOFF)).z, norm) -
                 M[2][2].im*half2fl((*(s+5*DEVOFF)).w, norm) );
(*(out+5)).w = ( M[2][0].re*half2fl((*(s+4*DEVOFF)).w, norm) +
                 M[2][0].im*half2fl((*(s+4*DEVOFF)).z, norm) ) 
             + ( M[2][1].re*half2fl((*(s+5*DEVOFF)).y, norm) +
                 M[2][1].im*half2fl((*(s+5*DEVOFF)).x, norm) ) 
             + ( M[2][2].re*half2fl((*(s+5*DEVOFF)).w, norm) +
                 M[2][2].im*half2fl((*(s+5*DEVOFF)).z, norm) );
}






//multipliziert su3-Matrix mal Spinor im Dirac-Raum
//we do not load the lower spinor components as we are working in the relativistic basis
//code in su3_MtV.txt -- generated with codegen
__device__ void dev_su3MtV_half_rel_up(dev_su3 M, const dev_spinor_half * s, const float * s_norm, dev_spinor * out){
float norm = * s_norm;

(*(out+0)).x = ( M[0][0].re*half2fl((*(s+0*DEVOFF)).x,norm) -
                 M[0][0].im*half2fl((*(s+0*DEVOFF)).y,norm) ) 
             + ( M[0][1].re*half2fl((*(s+0*DEVOFF)).z,norm) - 
                 M[0][1].im*half2fl((*(s+0*DEVOFF)).w,norm) ) 
             + ( M[0][2].re*half2fl((*(s+1*DEVOFF)).x,norm) -
                 M[0][2].im*half2fl((*(s+1*DEVOFF)).y,norm) );
(*(out+0)).y = ( M[0][0].re*half2fl((*(s+0*DEVOFF)).y,norm) +
                 M[0][0].im*half2fl((*(s+0*DEVOFF)).x,norm) ) 
             + ( M[0][1].re*half2fl((*(s+0*DEVOFF)).w,norm) +
                 M[0][1].im*half2fl((*(s+0*DEVOFF)).z,norm) ) 
             + ( M[0][2].re*half2fl((*(s+1*DEVOFF)).y,norm) +
                 M[0][2].im*half2fl((*(s+1*DEVOFF)).x,norm) );


(*(out+0)).z =  ( M[1][0].re*half2fl((*(s+0*DEVOFF)).x,norm) -
                  M[1][0].im*half2fl((*(s+0*DEVOFF)).y,norm) ) 
             +  ( M[1][1].re*half2fl((*(s+0*DEVOFF)).z,norm) - 
                  M[1][1].im*half2fl((*(s+0*DEVOFF)).w,norm) ) 
             +  ( M[1][2].re*half2fl((*(s+1*DEVOFF)).x,norm) -
                  M[1][2].im*half2fl((*(s+1*DEVOFF)).y,norm) );
(*(out+0)).w =  ( M[1][0].re*half2fl((*(s+0*DEVOFF)).y,norm) +
                  M[1][0].im*half2fl((*(s+0*DEVOFF)).x,norm) ) 
             +  ( M[1][1].re*half2fl((*(s+0*DEVOFF)).w,norm) +
                  M[1][1].im*half2fl((*(s+0*DEVOFF)).z,norm) ) 
             +  ( M[1][2].re*half2fl((*(s+1*DEVOFF)).y,norm) + 
                  M[1][2].im*half2fl((*(s+1*DEVOFF)).x,norm) );


(*(out+1)).x = ( M[2][0].re*half2fl((*(s+0*DEVOFF)).x, norm) -
                 M[2][0].im*half2fl((*(s+0*DEVOFF)).y, norm) ) 
             + ( M[2][1].re*half2fl((*(s+0*DEVOFF)).z, norm) -
                 M[2][1].im*half2fl((*(s+0*DEVOFF)).w, norm) ) 
             + ( M[2][2].re*half2fl((*(s+1*DEVOFF)).x, norm) -
                 M[2][2].im*half2fl((*(s+1*DEVOFF)).y, norm) );
(*(out+1)).y = ( M[2][0].re*half2fl((*(s+0*DEVOFF)).y, norm) +
                 M[2][0].im*half2fl((*(s+0*DEVOFF)).x, norm) ) 
             + ( M[2][1].re*half2fl((*(s+0*DEVOFF)).w, norm) +
                 M[2][1].im*half2fl((*(s+0*DEVOFF)).z, norm) ) 
             + ( M[2][2].re*half2fl((*(s+1*DEVOFF)).y, norm) +
                 M[2][2].im*half2fl((*(s+1*DEVOFF)).x, norm) );


(*(out+1)).z = ( M[0][0].re*half2fl((*(s+1*DEVOFF)).z, norm) -
                 M[0][0].im*half2fl((*(s+1*DEVOFF)).w, norm) ) 
             + ( M[0][1].re*half2fl((*(s+2*DEVOFF)).x, norm) -
                 M[0][1].im*half2fl((*(s+2*DEVOFF)).y, norm) ) 
             + ( M[0][2].re*half2fl((*(s+2*DEVOFF)).z, norm) -
                 M[0][2].im*half2fl((*(s+2*DEVOFF)).w, norm) );
(*(out+1)).w = ( M[0][0].re*half2fl((*(s+1*DEVOFF)).w, norm) +
                 M[0][0].im*half2fl((*(s+1*DEVOFF)).z, norm) ) 
             + ( M[0][1].re*half2fl((*(s+2*DEVOFF)).y, norm) +
                 M[0][1].im*half2fl((*(s+2*DEVOFF)).x, norm) ) 
             + ( M[0][2].re*half2fl((*(s+2*DEVOFF)).w, norm) +
                 M[0][2].im*half2fl((*(s+2*DEVOFF)).z, norm) );


(*(out+2)).x = ( M[1][0].re*half2fl((*(s+1*DEVOFF)).z, norm) -
                 M[1][0].im*half2fl((*(s+1*DEVOFF)).w, norm) ) 
             + ( M[1][1].re*half2fl((*(s+2*DEVOFF)).x, norm) -
                 M[1][1].im*half2fl((*(s+2*DEVOFF)).y, norm) ) 
             + ( M[1][2].re*half2fl((*(s+2*DEVOFF)).z, norm) -
                 M[1][2].im*half2fl((*(s+2*DEVOFF)).w, norm) );
(*(out+2)).y = ( M[1][0].re*half2fl((*(s+1*DEVOFF)).w, norm) +
                 M[1][0].im*half2fl((*(s+1*DEVOFF)).z, norm) ) 
             + ( M[1][1].re*half2fl((*(s+2*DEVOFF)).y, norm) +
                 M[1][1].im*half2fl((*(s+2*DEVOFF)).x, norm) ) 
             + ( M[1][2].re*half2fl((*(s+2*DEVOFF)).w, norm) +
                 M[1][2].im*half2fl((*(s+2*DEVOFF)).z, norm) );


(*(out+2)).z = ( M[2][0].re*half2fl((*(s+1*DEVOFF)).z, norm) -
                 M[2][0].im*half2fl((*(s+1*DEVOFF)).w, norm) ) 
             + ( M[2][1].re*half2fl((*(s+2*DEVOFF)).x, norm) -
                 M[2][1].im*half2fl((*(s+2*DEVOFF)).y, norm) ) 
             + ( M[2][2].re*half2fl((*(s+2*DEVOFF)).z, norm) -
                 M[2][2].im*half2fl((*(s+2*DEVOFF)).w, norm) );
(*(out+2)).w = ( M[2][0].re*half2fl((*(s+1*DEVOFF)).w, norm) +
                 M[2][0].im*half2fl((*(s+1*DEVOFF)).z, norm) ) 
             + ( M[2][1].re*half2fl((*(s+2*DEVOFF)).y, norm) +
                 M[2][1].im*half2fl((*(s+2*DEVOFF)).x, norm) ) 
             + ( M[2][2].re*half2fl((*(s+2*DEVOFF)).w, norm) +
                 M[2][2].im*half2fl((*(s+2*DEVOFF)).z, norm) );


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



//multipliziert su3-Matrix mal Spinor im Dirac-Raum
//we do not load the upper spinor components as we are working in the relativistic basis
//code in su3_MtV.txt -- generated with codegen
__device__ void dev_su3MtV_half_rel_down(dev_su3 M, const dev_spinor_half * s, const float * s_norm, dev_spinor * out){
float norm = * s_norm;

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


(*(out+3)).x = ( M[0][0].re*half2fl((*(s+3*DEVOFF)).x, norm) -
                 M[0][0].im*half2fl((*(s+3*DEVOFF)).y, norm) ) 
             + ( M[0][1].re*half2fl((*(s+3*DEVOFF)).z, norm) -
                 M[0][1].im*half2fl((*(s+3*DEVOFF)).w, norm) ) 
             + ( M[0][2].re*half2fl((*(s+4*DEVOFF)).x, norm) -
                 M[0][2].im*half2fl((*(s+4*DEVOFF)).y, norm) );
(*(out+3)).y = ( M[0][0].re*half2fl((*(s+3*DEVOFF)).y, norm) +
                 M[0][0].im*half2fl((*(s+3*DEVOFF)).x, norm) ) 
             + ( M[0][1].re*half2fl((*(s+3*DEVOFF)).w, norm) +
                 M[0][1].im*half2fl((*(s+3*DEVOFF)).z, norm) ) 
             + ( M[0][2].re*half2fl((*(s+4*DEVOFF)).y, norm) +
                 M[0][2].im*half2fl((*(s+4*DEVOFF)).x, norm) );


(*(out+3)).z = ( M[1][0].re*half2fl((*(s+3*DEVOFF)).x, norm) -
                 M[1][0].im*half2fl((*(s+3*DEVOFF)).y, norm) ) 
             + ( M[1][1].re*half2fl((*(s+3*DEVOFF)).z, norm) -
                 M[1][1].im*half2fl((*(s+3*DEVOFF)).w, norm) ) 
             + ( M[1][2].re*half2fl((*(s+4*DEVOFF)).x, norm) -
                 M[1][2].im*half2fl((*(s+4*DEVOFF)).y, norm) );
(*(out+3)).w = ( M[1][0].re*half2fl((*(s+3*DEVOFF)).y, norm) +
                 M[1][0].im*half2fl((*(s+3*DEVOFF)).x, norm) ) 
             + ( M[1][1].re*half2fl((*(s+3*DEVOFF)).w, norm) +
                 M[1][1].im*half2fl((*(s+3*DEVOFF)).z, norm) ) 
             + ( M[1][2].re*half2fl((*(s+4*DEVOFF)).y, norm) +
                 M[1][2].im*half2fl((*(s+4*DEVOFF)).x, norm) );


(*(out+4)).x = ( M[2][0].re*half2fl((*(s+3*DEVOFF)).x, norm) -
                 M[2][0].im*half2fl((*(s+3*DEVOFF)).y, norm) ) 
             + ( M[2][1].re*half2fl((*(s+3*DEVOFF)).z, norm) -
                 M[2][1].im*half2fl((*(s+3*DEVOFF)).w, norm) ) 
             + ( M[2][2].re*half2fl((*(s+4*DEVOFF)).x, norm) -
                 M[2][2].im*half2fl((*(s+4*DEVOFF)).y, norm) );
(*(out+4)).y = ( M[2][0].re*half2fl((*(s+3*DEVOFF)).y, norm) +
                 M[2][0].im*half2fl((*(s+3*DEVOFF)).x, norm) ) 
             + ( M[2][1].re*half2fl((*(s+3*DEVOFF)).w, norm) +
                 M[2][1].im*half2fl((*(s+3*DEVOFF)).z, norm) ) 
             + ( M[2][2].re*half2fl((*(s+4*DEVOFF)).y, norm) +
                 M[2][2].im*half2fl((*(s+4*DEVOFF)).x, norm) );


(*(out+4)).z = ( M[0][0].re*half2fl((*(s+4*DEVOFF)).z, norm) -
                 M[0][0].im*half2fl((*(s+4*DEVOFF)).w, norm) ) 
             + ( M[0][1].re*half2fl((*(s+5*DEVOFF)).x, norm) - 
                 M[0][1].im*half2fl((*(s+5*DEVOFF)).y, norm) ) 
             + ( M[0][2].re*half2fl((*(s+5*DEVOFF)).z, norm) -
                 M[0][2].im*half2fl((*(s+5*DEVOFF)).w, norm) );
(*(out+4)).w = ( M[0][0].re*half2fl((*(s+4*DEVOFF)).w, norm) +
                 M[0][0].im*half2fl((*(s+4*DEVOFF)).z, norm) ) 
             + ( M[0][1].re*half2fl((*(s+5*DEVOFF)).y, norm) +
                 M[0][1].im*half2fl((*(s+5*DEVOFF)).x, norm) ) 
             + ( M[0][2].re*half2fl((*(s+5*DEVOFF)).w, norm) +
                 M[0][2].im*half2fl((*(s+5*DEVOFF)).z, norm) );


(*(out+5)).x = ( M[1][0].re*half2fl((*(s+4*DEVOFF)).z, norm) -
                 M[1][0].im*half2fl((*(s+4*DEVOFF)).w, norm) ) 
             + ( M[1][1].re*half2fl((*(s+5*DEVOFF)).x, norm) -
                 M[1][1].im*half2fl((*(s+5*DEVOFF)).y, norm) ) 
             + ( M[1][2].re*half2fl((*(s+5*DEVOFF)).z, norm) -
                 M[1][2].im*half2fl((*(s+5*DEVOFF)).w, norm) );
(*(out+5)).y = ( M[1][0].re*half2fl((*(s+4*DEVOFF)).w, norm) +
                 M[1][0].im*half2fl((*(s+4*DEVOFF)).z, norm) ) 
             + ( M[1][1].re*half2fl((*(s+5*DEVOFF)).y, norm) +
                 M[1][1].im*half2fl((*(s+5*DEVOFF)).x, norm) ) 
             + ( M[1][2].re*half2fl((*(s+5*DEVOFF)).w, norm) +
                 M[1][2].im*half2fl((*(s+5*DEVOFF)).z, norm) );


(*(out+5)).z = ( M[2][0].re*half2fl((*(s+4*DEVOFF)).z, norm) - 
                 M[2][0].im*half2fl((*(s+4*DEVOFF)).w, norm) ) 
             + ( M[2][1].re*half2fl((*(s+5*DEVOFF)).x, norm) -
                 M[2][1].im*half2fl((*(s+5*DEVOFF)).y, norm) ) 
             + ( M[2][2].re*half2fl((*(s+5*DEVOFF)).z, norm) -
                 M[2][2].im*half2fl((*(s+5*DEVOFF)).w, norm) );
(*(out+5)).w = ( M[2][0].re*half2fl((*(s+4*DEVOFF)).w, norm) +
                 M[2][0].im*half2fl((*(s+4*DEVOFF)).z, norm) ) 
             + ( M[2][1].re*half2fl((*(s+5*DEVOFF)).y, norm) +
                 M[2][1].im*half2fl((*(s+5*DEVOFF)).x, norm) ) 
             + ( M[2][2].re*half2fl((*(s+5*DEVOFF)).w, norm) +
                 M[2][2].im*half2fl((*(s+5*DEVOFF)).z, norm) );
}




#endif







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
  float tempre,tempim;
     tempre = (*(in+0)).x;
     tempim = (*(in+0)).y;
     (*(in+0)).x = -(*(in+3)).x;
     (*(in+0)).y = -(*(in+3)).y;
     (*(in+3)).x = -tempre;
     (*(in+3)).y = -tempim;     
     
     tempre = (*(in+0)).z;
     tempim = (*(in+0)).w;
     (*(in+0)).z = -(*(in+3)).z;
     (*(in+0)).w = -(*(in+3)).w;
     (*(in+3)).z = -tempre;
     (*(in+3)).w = -tempim; 
 
 
     tempre = (*(in+1)).x;
     tempim = (*(in+1)).y;
     (*(in+1)).x = -(*(in+4)).x;
     (*(in+1)).y = -(*(in+4)).y;
     (*(in+4)).x = -tempre;
     (*(in+4)).y = -tempim;     
     
     tempre = (*(in+1)).z;
     tempim = (*(in+1)).w;
     (*(in+1)).z = -(*(in+4)).z;
     (*(in+1)).w = -(*(in+4)).w;
     (*(in+4)).z = -tempre;
     (*(in+4)).w = -tempim;     
     
     
     tempre = (*(in+2)).x;
     tempim = (*(in+2)).y;
     (*(in+2)).x = -(*(in+5)).x;
     (*(in+2)).y = -(*(in+5)).y;
     (*(in+5)).x = -tempre;
     (*(in+5)).y = -tempim;     
   
   
     tempre = (*(in+2)).z;
     tempim = (*(in+2)).w;
     (*(in+2)).z = -(*(in+5)).z;
     (*(in+2)).w = -(*(in+5)).w;
     (*(in+5)).z = -tempre;
     (*(in+5)).w = -tempim;
}



//Gamma z
__device__ void dev_Gamma3(dev_spinor * in){
  float tempre,tempim;
     tempre = (*(in+0)).x;
     tempim = (*(in+0)).y;
     (*(in+0)).x = (*(in+3)).y;
     (*(in+0)).y = -(*(in+3)).x;
     (*(in+3)).x = -tempim;
     (*(in+3)).y = tempre;    
     
     tempre = (*(in+0)).z;
     tempim = (*(in+0)).w;
     (*(in+0)).z = (*(in+3)).w;
     (*(in+0)).w = -(*(in+3)).z;
     (*(in+3)).z = -tempim;
     (*(in+3)).w = tempre;    
     
     
     tempre = (*(in+1)).x;
     tempim = (*(in+1)).y;
     (*(in+1)).x = (*(in+4)).y;
     (*(in+1)).y = -(*(in+4)).x;
     (*(in+4)).x = -tempim;
     (*(in+4)).y = tempre;     
     
     
     tempre = (*(in+1)).z;
     tempim = (*(in+1)).w;
     (*(in+1)).z = -(*(in+4)).w;
     (*(in+1)).w = (*(in+4)).z;
     (*(in+4)).z  = tempim;
     (*(in+4)).w  = -tempre;     
     
     
     
     tempre = (*(in+2)).x;
     tempim = (*(in+2)).y;
     (*(in+2)).x = -(*(in+5)).y;
     (*(in+2)).y = (*(in+5)).x;
     (*(in+5)).x = tempim;
     (*(in+5)).y = -tempre;    
     
     
     tempre = (*(in+2)).z;
     tempim = (*(in+2)).w;
     (*(in+2)).z = -(*(in+5)).w;
     (*(in+2)).w = (*(in+5)).z;
     (*(in+5)).z = tempim;
     (*(in+5)).w = -tempre;

}



//Gamma y
__device__ void dev_Gamma2(dev_spinor * in){
  float tempre,tempim;
     tempre = (*(in+0)).x;
     tempim = (*(in+0)).y;
     (*(in+0)).x = -(*(in+4)).z;
     (*(in+0)).y = -(*(in+4)).w;
     (*(in+4)).z = -tempre;
     (*(in+4)).w = -tempim;    
     
     tempre = (*(in+0)).z;
     tempim = (*(in+0)).w;
     (*(in+0)).z = -(*(in+5)).x;
     (*(in+0)).w = -(*(in+5)).y;
     (*(in+5)).x = -tempre;
     (*(in+5)).y = -tempim;     
     
     
     tempre = (*(in+1)).x;
     tempim = (*(in+1)).y;
     (*(in+1)).x = -(*(in+5)).z;
     (*(in+1)).y = -(*(in+5)).w;
     (*(in+5)).z = -tempre;
     (*(in+5)).w = -tempim;     
     
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
  float tempre,tempim;
     tempre = (*(in+0)).x;
     tempim = (*(in+0)).y;
     (*(in+0)).x = (*(in+4)).w;
     (*(in+0)).y = -(*(in+4)).z;
     (*(in+4)).z  = -tempim;
     (*(in+4)).w  = tempre;    
     
     tempre = (*(in+0)).z;
     tempim = (*(in+0)).w;
     (*(in+0)).z = (*(in+5)).y;
     (*(in+0)).w = -(*(in+5)).x;
     (*(in+5)).x = -tempim;
     (*(in+5)).y = tempre;     
     
     tempre = (*(in+1)).x;
     tempim = (*(in+1)).y;
     (*(in+1)).x = (*(in+5)).w;
     (*(in+1)).y = -(*(in+5)).z;
     (*(in+5)).z = -tempim;
     (*(in+5)).w = tempre;     
     
     tempre = (*(in+1)).z;
     tempim = (*(in+1)).w;
     (*(in+1)).z = (*(in+3)).y;
     (*(in+1)).w = -(*(in+3)).x;
     (*(in+3)).x = -tempim;
     (*(in+3)).y = tempre;     
     
     tempre = (*(in+2)).x;
     tempim = (*(in+2)).y;
     (*(in+2)).x = (*(in+3)).w;
     (*(in+2)).y = -(*(in+3)).z;
     (*(in+3)).z = -tempim;
     (*(in+3)).w = tempre;     
     
     
     tempre = (*(in+2)).z;
     tempim = (*(in+2)).w;
     (*(in+2)).z = (*(in+4)).y;
     (*(in+2)).w = -(*(in+4)).x;
     (*(in+4)).x = -tempim;
     (*(in+4)).y = tempre;
  
}



__device__ void dev_Gamma5(dev_spinor * in){
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



//this works in the relativistic basis
//gamma5 NOT diagonal!!
__device__ void dev_Gamma5_rel(dev_spinor* in){

dev_spinor help[6]; 
  help[0].x = -1.0*(*(in+3)).x;
  help[0].y = -1.0*(*(in+3)).y;
  help[0].z = -1.0*(*(in+3)).z;
  help[0].w = -1.0*(*(in+3)).w;
  help[1].x = -1.0*(*(in+4)).x;
  help[1].y = -1.0*(*(in+4)).y;

  help[1].z = -1.0*(*(in+4)).z;
  help[1].w = -1.0*(*(in+4)).w;
  help[2].x = -1.0*(*(in+5)).x;
  help[2].y = -1.0*(*(in+5)).y;
  help[2].z = -1.0*(*(in+5)).z;
  help[2].w = -1.0*(*(in+5)).w;

  help[3].x = -1.0*(*(in+0)).x;
  help[3].y = -1.0*(*(in+0)).y;
  help[3].z = -1.0*(*(in+0)).z;
  help[3].w = -1.0*(*(in+0)).w;
  help[4].x = -1.0*(*(in+1)).x;
  help[4].y = -1.0*(*(in+1)).y;

  help[4].z = -1.0*(*(in+1)).z;
  help[4].w = -1.0*(*(in+1)).w;
  help[5].x = -1.0*(*(in+2)).x;
  help[5].y = -1.0*(*(in+2)).y;
  help[5].z = -1.0*(*(in+2)).z;
  help[5].w = -1.0*(*(in+2)).w;
  
  (*(in+0)).x = help[0].x;
  (*(in+0)).y = help[0].y; 
  (*(in+0)).z = help[0].z;  
  (*(in+0)).w = help[0].w;  
  
  (*(in+1)).x = help[1].x;
  (*(in+1)).y = help[1].y; 
  (*(in+1)).z = help[1].z;  
  (*(in+1)).w = help[1].w;   
  
  (*(in+2)).x = help[2].x;
  (*(in+2)).y = help[2].y; 
  (*(in+2)).z = help[2].z;  
  (*(in+2)).w = help[2].w;   
  
  (*(in+3)).x = help[3].x;
  (*(in+3)).y = help[3].y; 
  (*(in+3)).z = help[3].z;  
  (*(in+3)).w = help[3].w;    

  (*(in+4)).x = help[4].x;
  (*(in+4)).y = help[4].y; 
  (*(in+4)).z = help[4].z;  
  (*(in+4)).w = help[4].w;    
  
  (*(in+5)).x = help[5].x;
  (*(in+5)).y = help[5].y; 
  (*(in+5)).z = help[5].z;  
  (*(in+5)).w = help[5].w;    
  
}






__device__ void dev_Gamma5_assign(dev_spinor* out, dev_spinor* in){
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



__device__ void dev_Gamma5_assigntoglobal(dev_spinor* out, dev_spinor* in){
  (*(out)).x = (*(in)).x;
  (*(out)).y = (*(in)).y;
  (*(out)).z = (*(in)).z;
  (*(out)).w = (*(in)).w;
  (*(out+1*DEVOFF)).x = (*(in+1)).x;
  (*(out+1*DEVOFF)).y = (*(in+1)).y;

  (*(out+1*DEVOFF)).z = (*(in+1)).z;
  (*(out+1*DEVOFF)).w = (*(in+1)).w;
  (*(out+2*DEVOFF)).x = (*(in+2)).x;
  (*(out+2*DEVOFF)).y = (*(in+2)).y;
  (*(out+2*DEVOFF)).z = (*(in+2)).z;
  (*(out+2*DEVOFF)).w = (*(in+2)).w;

  (*(out+3*DEVOFF)).x = -1.0*(*(in+3)).x;
  (*(out+3*DEVOFF)).y = -1.0*(*(in+3)).y;
  (*(out+3*DEVOFF)).z = -1.0*(*(in+3)).z;
  (*(out+3*DEVOFF)).w = -1.0*(*(in+3)).w;
  (*(out+4*DEVOFF)).x = -1.0*(*(in+4)).x;
  (*(out+4*DEVOFF)).y = -1.0*(*(in+4)).y;

  (*(out+4*DEVOFF)).z = -1.0*(*(in+4)).z;
  (*(out+4*DEVOFF)).w = -1.0*(*(in+4)).w;
  (*(out+5*DEVOFF)).x = -1.0*(*(in+5)).x;
  (*(out+5*DEVOFF)).y = -1.0*(*(in+5)).y;
  (*(out+5*DEVOFF)).z = -1.0*(*(in+5)).z;
  (*(out+5*DEVOFF)).w = -1.0*(*(in+5)).w;
}



//this works in the relativistic basis
//gamma5 NOT diagonal!!
__device__ void dev_Gamma5_assign_rel(dev_spinor* out, dev_spinor* in){
  (*(out+0)).x = -1.0*(*(in+3)).x;
  (*(out+0)).y = -1.0*(*(in+3)).y;
  (*(out+0)).z = -1.0*(*(in+3)).z;
  (*(out+0)).w = -1.0*(*(in+3)).w;
  (*(out+1)).x = -1.0*(*(in+4)).x;
  (*(out+1)).y = -1.0*(*(in+4)).y;

  (*(out+1)).z = -1.0*(*(in+4)).z;
  (*(out+1)).w = -1.0*(*(in+4)).w;
  (*(out+2)).x = -1.0*(*(in+5)).x;
  (*(out+2)).y = -1.0*(*(in+5)).y;
  (*(out+2)).z = -1.0*(*(in+5)).z;
  (*(out+2)).w = -1.0*(*(in+5)).w;

  (*(out+3)).x = -1.0*(*(in+0)).x;
  (*(out+3)).y = -1.0*(*(in+0)).y;
  (*(out+3)).z = -1.0*(*(in+0)).z;
  (*(out+3)).w = -1.0*(*(in+0)).w;
  (*(out+4)).x = -1.0*(*(in+1)).x;
  (*(out+4)).y = -1.0*(*(in+1)).y;

  (*(out+4)).z = -1.0*(*(in+1)).z;
  (*(out+4)).w = -1.0*(*(in+1)).w;
  (*(out+5)).x = -1.0*(*(in+2)).x;
  (*(out+5)).y = -1.0*(*(in+2)).y;
  (*(out+5)).z = -1.0*(*(in+2)).z;
  (*(out+5)).w = -1.0*(*(in+2)).w;

}



//this works in the relativistic basis
//gamma5 NOT diagonal!!
__device__ void dev_Gamma5_assigntoglobal_rel(dev_spinor* out, dev_spinor* in){
  (*(out+0*DEVOFF)).x = -1.0*(*(in+3)).x;
  (*(out+0*DEVOFF)).y = -1.0*(*(in+3)).y;
  (*(out+0*DEVOFF)).z = -1.0*(*(in+3)).z;
  (*(out+0*DEVOFF)).w = -1.0*(*(in+3)).w;
  (*(out+1*DEVOFF)).x = -1.0*(*(in+4)).x;
  (*(out+1*DEVOFF)).y = -1.0*(*(in+4)).y;

  (*(out+1*DEVOFF)).z = -1.0*(*(in+4)).z;
  (*(out+1*DEVOFF)).w = -1.0*(*(in+4)).w;
  (*(out+2*DEVOFF)).x = -1.0*(*(in+5)).x;
  (*(out+2*DEVOFF)).y = -1.0*(*(in+5)).y;
  (*(out+2*DEVOFF)).z = -1.0*(*(in+5)).z;
  (*(out+2*DEVOFF)).w = -1.0*(*(in+5)).w;

  (*(out+3*DEVOFF)).x = -1.0*(*(in+0)).x;
  (*(out+3*DEVOFF)).y = -1.0*(*(in+0)).y;
  (*(out+3*DEVOFF)).z = -1.0*(*(in+0)).z;
  (*(out+3*DEVOFF)).w = -1.0*(*(in+0)).w;
  (*(out+4*DEVOFF)).x = -1.0*(*(in+1)).x;
  (*(out+4*DEVOFF)).y = -1.0*(*(in+1)).y;

  (*(out+4*DEVOFF)).z = -1.0*(*(in+1)).z;
  (*(out+4*DEVOFF)).w = -1.0*(*(in+1)).w;
  (*(out+5*DEVOFF)).x = -1.0*(*(in+2)).x;
  (*(out+5*DEVOFF)).y = -1.0*(*(in+2)).y;
  (*(out+5*DEVOFF)).z = -1.0*(*(in+2)).z;
  (*(out+5*DEVOFF)).w = -1.0*(*(in+2)).w;

}





// older version, all in one function
__device__ void dev_GammatV(int mu, dev_spinor * in){//multipliziert Gamma(mu)*V effizientes ausnutzen der Nullen 
 float tempre,tempim;
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





/////////////////////  gamma basis transformations //////////////////////////////////
//
// to go to relativistic basis we have:
// g0_rel = t1 . g0_tmlqcd  . t1^+                                    
// g5_rel = t1 . g5_tmlqcd  . t1^+
//
// t1 =
//              |  1  0  1  0 |
//  1/sqrt(2) * |  0  1  0  1 |
//              | -1  0  1  0 |
//              |  0 -1  0  1 |
//



__global__ void to_relativistic_basis(dev_spinor* spinin){
   int pos= threadIdx.x + blockDim.x*blockIdx.x;
   float4 inhelp[6];
   float4 outhelp[6];

   const float sq2 = rsqrtf(2.0f);
   
   if(pos < dev_VOLUME){
    dev_read_spinor(&(inhelp[0]), &(spinin[pos]));
    
    outhelp[0].x = sq2*(inhelp[0].x + inhelp[3].x);
    outhelp[0].y = sq2*(inhelp[0].y + inhelp[3].y);
    outhelp[0].z = sq2*(inhelp[0].z + inhelp[3].z);
    outhelp[0].w = sq2*(inhelp[0].w + inhelp[3].w);
    outhelp[1].x = sq2*(inhelp[1].x + inhelp[4].x);
    outhelp[1].y = sq2*(inhelp[1].y + inhelp[4].y);
    
    
    outhelp[1].z = sq2*(inhelp[1].z + inhelp[4].z);
    outhelp[1].w = sq2*(inhelp[1].w + inhelp[4].w);
    outhelp[2].x = sq2*(inhelp[2].x + inhelp[5].x);
    outhelp[2].y = sq2*(inhelp[2].y + inhelp[5].y);
    outhelp[2].z = sq2*(inhelp[2].z + inhelp[5].z);
    outhelp[2].w = sq2*(inhelp[2].w + inhelp[5].w); 
    
    
    outhelp[3].x = sq2*(inhelp[3].x - inhelp[0].x);
    outhelp[3].y = sq2*(inhelp[3].y - inhelp[0].y);
    outhelp[3].z = sq2*(inhelp[3].z - inhelp[0].z);
    outhelp[3].w = sq2*(inhelp[3].w - inhelp[0].w);
    outhelp[4].x = sq2*(inhelp[4].x - inhelp[1].x);
    outhelp[4].y = sq2*(inhelp[4].y - inhelp[1].y);   
    
    
    outhelp[4].z = sq2*(inhelp[4].z - inhelp[1].z);
    outhelp[4].w = sq2*(inhelp[4].w - inhelp[1].w);
    outhelp[5].x = sq2*(inhelp[5].x - inhelp[2].x);
    outhelp[5].y = sq2*(inhelp[5].y - inhelp[2].y);
    outhelp[5].z = sq2*(inhelp[5].z - inhelp[2].z);
    outhelp[5].w = sq2*(inhelp[5].w - inhelp[2].w);    

    //copy to output spinor
      dev_write_spinor(&(outhelp[0]),&(spinin[pos])); 
   }//dev_VOLUME

}



__device__ void to_relativistic_basis_spinor(dev_spinor* spinin){
   float4 outhelp[6];

   const float sq2 = rsqrtf(2.0f);
   
   
    outhelp[0].x = sq2*(spinin[0].x + spinin[3].x);
    outhelp[0].y = sq2*(spinin[0].y + spinin[3].y);
    outhelp[0].z = sq2*(spinin[0].z + spinin[3].z);
    outhelp[0].w = sq2*(spinin[0].w + spinin[3].w);
    outhelp[1].x = sq2*(spinin[1].x + spinin[4].x);
    outhelp[1].y = sq2*(spinin[1].y + spinin[4].y);
    
    
    outhelp[1].z = sq2*(spinin[1].z + spinin[4].z);
    outhelp[1].w = sq2*(spinin[1].w + spinin[4].w);
    outhelp[2].x = sq2*(spinin[2].x + spinin[5].x);
    outhelp[2].y = sq2*(spinin[2].y + spinin[5].y);
    outhelp[2].z = sq2*(spinin[2].z + spinin[5].z);
    outhelp[2].w = sq2*(spinin[2].w + spinin[5].w); 
    
    
    outhelp[3].x = sq2*(spinin[3].x - spinin[0].x);
    outhelp[3].y = sq2*(spinin[3].y - spinin[0].y);
    outhelp[3].z = sq2*(spinin[3].z - spinin[0].z);
    outhelp[3].w = sq2*(spinin[3].w - spinin[0].w);
    outhelp[4].x = sq2*(spinin[4].x - spinin[1].x);
    outhelp[4].y = sq2*(spinin[4].y - spinin[1].y);   
    
    
    outhelp[4].z = sq2*(spinin[4].z - spinin[1].z);
    outhelp[4].w = sq2*(spinin[4].w - spinin[1].w);
    outhelp[5].x = sq2*(spinin[5].x - spinin[2].x);
    outhelp[5].y = sq2*(spinin[5].y - spinin[2].y);
    outhelp[5].z = sq2*(spinin[5].z - spinin[2].z);
    outhelp[5].w = sq2*(spinin[5].w - spinin[2].w);    

   dev_copy_spinor_local(&(outhelp[0]), spinin);

}


__global__ void to_tmlqcd_basis(dev_spinor* spinin){
  
   int pos= threadIdx.x + blockDim.x*blockIdx.x;
   float4 inhelp[6];
   float4 outhelp[6];

   const float sq2 = rsqrtf(2.0f);
   
   if(pos < dev_VOLUME){
    dev_read_spinor(&(inhelp[0]), &(spinin[pos]));
    
    outhelp[0].x = sq2*(inhelp[0].x - inhelp[3].x);
    outhelp[0].y = sq2*(inhelp[0].y - inhelp[3].y);
    outhelp[0].z = sq2*(inhelp[0].z - inhelp[3].z);
    outhelp[0].w = sq2*(inhelp[0].w - inhelp[3].w);
    outhelp[1].x = sq2*(inhelp[1].x - inhelp[4].x);
    outhelp[1].y = sq2*(inhelp[1].y - inhelp[4].y);
    
    
    outhelp[1].z = sq2*(inhelp[1].z - inhelp[4].z);
    outhelp[1].w = sq2*(inhelp[1].w - inhelp[4].w);
    outhelp[2].x = sq2*(inhelp[2].x - inhelp[5].x);
    outhelp[2].y = sq2*(inhelp[2].y - inhelp[5].y);
    outhelp[2].z = sq2*(inhelp[2].z - inhelp[5].z);
    outhelp[2].w = sq2*(inhelp[2].w - inhelp[5].w); 
    
    
    outhelp[3].x = sq2*(inhelp[3].x + inhelp[0].x);
    outhelp[3].y = sq2*(inhelp[3].y + inhelp[0].y);
    outhelp[3].z = sq2*(inhelp[3].z + inhelp[0].z);
    outhelp[3].w = sq2*(inhelp[3].w + inhelp[0].w);
    outhelp[4].x = sq2*(inhelp[4].x + inhelp[1].x);
    outhelp[4].y = sq2*(inhelp[4].y + inhelp[1].y);   
    
    
    outhelp[4].z = sq2*(inhelp[4].z + inhelp[1].z);
    outhelp[4].w = sq2*(inhelp[4].w + inhelp[1].w);
    outhelp[5].x = sq2*(inhelp[5].x + inhelp[2].x);
    outhelp[5].y = sq2*(inhelp[5].y + inhelp[2].y);
    outhelp[5].z = sq2*(inhelp[5].z + inhelp[2].z);
    outhelp[5].w = sq2*(inhelp[5].w + inhelp[2].w);    

    //copy to output spinor
      dev_write_spinor(&(outhelp[0]),&(spinin[pos]));     
   }//dev_VOLUME


}


__device__ void to_tmlqcd_basis_spinor(dev_spinor* spinin){
   float4 outhelp[6];
   
   const float sq2 = rsqrtf(2.0f);
   
    outhelp[0].x = sq2*(spinin[0].x - spinin[3].x);
    outhelp[0].y = sq2*(spinin[0].y - spinin[3].y);
    outhelp[0].z = sq2*(spinin[0].z - spinin[3].z);
    outhelp[0].w = sq2*(spinin[0].w - spinin[3].w);
    outhelp[1].x = sq2*(spinin[1].x - spinin[4].x);
    outhelp[1].y = sq2*(spinin[1].y - spinin[4].y);
    
    
    outhelp[1].z = sq2*(spinin[1].z - spinin[4].z);
    outhelp[1].w = sq2*(spinin[1].w - spinin[4].w);
    outhelp[2].x = sq2*(spinin[2].x - spinin[5].x);
    outhelp[2].y = sq2*(spinin[2].y - spinin[5].y);
    outhelp[2].z = sq2*(spinin[2].z - spinin[5].z);
    outhelp[2].w = sq2*(spinin[2].w - spinin[5].w); 
    
    
    outhelp[3].x = sq2*(spinin[3].x + spinin[0].x);
    outhelp[3].y = sq2*(spinin[3].y + spinin[0].y);
    outhelp[3].z = sq2*(spinin[3].z + spinin[0].z);
    outhelp[3].w = sq2*(spinin[3].w + spinin[0].w);
    outhelp[4].x = sq2*(spinin[4].x + spinin[1].x);
    outhelp[4].y = sq2*(spinin[4].y + spinin[1].y);   
    
    
    outhelp[4].z = sq2*(spinin[4].z + spinin[1].z);
    outhelp[4].w = sq2*(spinin[4].w + spinin[1].w);
    outhelp[5].x = sq2*(spinin[5].x + spinin[2].x);
    outhelp[5].y = sq2*(spinin[5].y + spinin[2].y);
    outhelp[5].z = sq2*(spinin[5].z + spinin[2].z);
    outhelp[5].w = sq2*(spinin[5].w + spinin[2].w);    

    //copy to output spinor
    dev_copy_spinor_local(&(outhelp[0]), spinin);
}



////////////////////////  BLAS KERNELS FLOAT ////////////////////////////////////////////////////




  int blas_gridsize;
  int blas_blocksize; // kernel parameters for the half_dot and axpy kernels
  float * dev_blas_redfield;  //this is the reduction field for the
                               //blas reduction kernels
  float * dev_blas_sredfield; // this is the small reduction field after one sweep of reduction
  float * blas_sredfield;                             
  int blas_redblocks;  // the number of blocks of the reduction kernel
                     // VOLUME/REDUCTION_N
                     // also the size of the final sum (of reduction)
                     // performed on host



void init_blas(int vol){
  cudaError_t cudaerr;
  
  blas_blocksize=BLOCK2;
  if( vol >= BLOCK2){
   blas_gridsize = (int)(vol/BLOCK2) + 1;
  }
  else{
    blas_gridsize=1;
  }
  
  size_t size = vol * sizeof(float);
  
  if((cudaerr=cudaMalloc((void **) &dev_blas_redfield, size)) != cudaSuccess){
    if(g_proc_id==0){
      printf("Error in init_blas(): Memory allocation of reduction field failed. Aborting...\n");
      printf("Error code is: %f\n",cudaerr);
    }
    exit(200);
  }   // Allocate array on device
  else{
    #ifndef LOWOUTPUT 
    if(g_proc_id==0) printf("Allocated blas reduction field on device\n");
    #endif
  }  


  // IMPLEMENT THIS FOR ALL LATTICE SIZES !!!!!!!!!!!!!!!!!!!!
  if((vol%REDUCTION_N) == 0){
    blas_redblocks = vol/REDUCTION_N;
  }
  else{
    if(g_proc_id==0) fprintf(stderr,"Error: Volume is not a multiple of REDUCTION_N (%d). Aborting...\n", REDUCTION_N);
    exit(100);
  }
  
  // initialize small redfields
  size = blas_redblocks * sizeof(float);
  if((cudaerr=cudaMalloc((void **) &dev_blas_sredfield, size)) != cudaSuccess){
    if(g_proc_id==0) {
      printf("Error in init_blas(): Memory allocation of small reduction field failed. Aborting...\n");
      printf("Error code is: %f\n",cudaerr);
    }
    exit(200);
  }   // Allocate array on device
  else{
    #ifndef LOWOUTPUT 
    if(g_proc_id==0) printf("Allocated blas small reduction field on device\n");
    #endif
  }  
  
  if((void*)(blas_sredfield = (float *)malloc(size)) == NULL){
    if(g_proc_id==0) {
      printf("Could not allocate memory for blas small redfield on host. Aborting...\n");
      printf("Error code is: %f\n",cudaerr);
    }
    exit(200);
  } 
  
  
  
}


void finalize_blas(){
  cudaFree(dev_blas_redfield);
  cudaFree(dev_blas_sredfield);
  free(blas_sredfield);
}




// this is a reduction algorithm for float based on the CUDA SDK 
__global__ void reduce_float(float *g_idata, float *g_odata, unsigned int n)
{
    extern __shared__ float sdata[];

    // load shared mem
    unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x*blockDim.x + threadIdx.x;
    
    sdata[tid] = (i < n) ? g_idata[i] : 0;
    
    __syncthreads();

    // do reduction in shared mem
    for(unsigned int s=blockDim.x/2; s>0; s>>=1) 
    {
        if (tid < s) 
        {
            sdata[tid] += sdata[tid + s];
        }
        __syncthreads();
    }

    // write result for this block to global mem
    if (tid == 0) g_odata[blockIdx.x] = sdata[0];
}




// this is the version for float2 
__global__ void reduce_float2(float2 *g_idata, float2 *g_odata, unsigned int n)
{
    extern __shared__ float2 sdata2[];

    // load shared mem
    unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x*blockDim.x + threadIdx.x;
    
    sdata2[tid].x = (i < n) ? g_idata[i].x : 0;
    sdata2[tid].y = (i < n) ? g_idata[i].y : 0;
    
    __syncthreads();

    // do reduction in shared mem
    for(unsigned int s=blockDim.x/2; s>0; s>>=1) 
    {
        if (tid < s) 
        {
            sdata2[tid].x += sdata2[tid + s].x;
            sdata2[tid].y += sdata2[tid + s].y;
        }
        __syncthreads();
    }

    // write result for this block to global mem
    if (tid == 0) {
      g_odata[blockIdx.x].x = sdata2[0].x;
      g_odata[blockIdx.x].y = sdata2[0].y;
    }
}







// y = alpha*x + y 
// x is not read from texture
// y is not read from texture
__global__ void dev_axpy (float alpha, dev_spinor* x, dev_spinor* y){
   int pos= threadIdx.x + blockDim.x*blockIdx.x;
   float4 xhelp[6]; 
   float4 erghelp[6];
   int i;

   
   if(pos < dev_VOLUME){
   
   //load y
   dev_read_spinor(&(erghelp[0]), &(y[pos]));
   //load x
   dev_read_spinor(&(xhelp[0]), &(x[pos]));

    
    #pragma unroll 6
    for(i=0; i<6; i++){
       erghelp[i].x += alpha*xhelp[i].x;
       erghelp[i].y += alpha*xhelp[i].y;
       erghelp[i].z += alpha*xhelp[i].z;
       erghelp[i].w += alpha*xhelp[i].w;
    }
        
     //write out spinors
     dev_write_spinor(&(erghelp[0]),&(y[pos])); 
   }//dev_VOLUME
}


// y = x + alpha*y 
// x is not read from texture
// y is not read from texture
__global__ void dev_xpay (float alpha, dev_spinor* x, dev_spinor* y){
   int pos= threadIdx.x + blockDim.x*blockIdx.x;
   float4 yhelp[6]; 
   float4 erghelp[6];
   int i;

   
   if(pos < dev_VOLUME){
   
   //load x
   dev_read_spinor(&(erghelp[0]), &(x[pos]));
   //load y
   dev_read_spinor(&(yhelp[0]), &(y[pos]));

    
    #pragma unroll 6
    for(i=0; i<6; i++){
       erghelp[i].x += alpha*yhelp[i].x;
       erghelp[i].y += alpha*yhelp[i].y;
       erghelp[i].z += alpha*yhelp[i].z;
       erghelp[i].w += alpha*yhelp[i].w;
    }
        
     //write out spinors
     dev_write_spinor(&(erghelp[0]),&(y[pos])); 
   }//dev_VOLUME
}



// y = alpha*y 
// y is not read from texture
__global__ void dev_blasscal (float alpha, dev_spinor* y){
   int pos= threadIdx.x + blockDim.x*blockIdx.x;
   float4 yhelp[6]; 
   float4 erghelp[6];
   int i;

   
   if(pos < dev_VOLUME){
   
   //load y
   dev_read_spinor(&(yhelp[0]), &(y[pos]));

    
    #pragma unroll 6
    for(i=0; i<6; i++){
       erghelp[i].x = alpha*yhelp[i].x;
       erghelp[i].y = alpha*yhelp[i].y;
       erghelp[i].z = alpha*yhelp[i].z;
       erghelp[i].w = alpha*yhelp[i].w;
    }
        
     //write out spinors
     dev_write_spinor(&(erghelp[0]),&(y[pos])); 
   }//dev_VOLUME
}



// y = x 
// x is not read from texture
__global__ void dev_blascopy (dev_spinor* x, dev_spinor* y){
   int pos= threadIdx.x + blockDim.x*blockIdx.x;
   float4 xhelp[6]; 

   
   if(pos < dev_VOLUME){
   
   //load y
     dev_read_spinor(&(xhelp[0]), &(x[pos]));
   //write out spinor
     dev_write_spinor(&(xhelp[0]),&(y[pos])); 
   }//dev_VOLUME
}


//this is a dotprod implementation for float
//x*y at spacetime point x is put into redfield at pos
__global__ void dev_dot( float* redfield, dev_spinor* x,dev_spinor* y){
   int pos= threadIdx.x + blockDim.x*blockIdx.x;
   float4 xhelp[6],yhelp[6];
   int i;
   float dotp = 0.0f;
   
   if(pos < dev_VOLUME){
    // this is the loop over the 6 float4 forming one spinor
    
    //load y
      dev_read_spinor(&(yhelp[0]), &(y[pos]));
    //load x
      dev_read_spinor(&(xhelp[0]), &(x[pos]));

    
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
float float_dotprod(dev_spinor* x, dev_spinor* y){
   int i;

   cudaError_t cudaerr;
   
   dev_dot<<< blas_gridsize, blas_blocksize >>> 
                      (dev_blas_redfield, x, y);
   if((cudaerr=cudaGetLastError()) != cudaSuccess){
      if(g_proc_id==0) printf("%s\n", cudaGetErrorString(cudaerr));
      if(g_proc_id==0) printf("Error in float_dotprod. dev_dot kernel erroneous. Aborting...\n");
      exit(200);
   } 
   //reduce reductionfield on device 
   reduce_float <<< blas_redblocks, REDUCTION_N, 
                REDUCTION_N*sizeof(float) >>> 
                ( dev_blas_redfield, dev_blas_sredfield,  VOLUME);
   //this reduction always takes the VOLUME (also for mpi)     
   
   //copy back
   cudaMemcpy(blas_sredfield, dev_blas_sredfield, (size_t)(blas_redblocks*sizeof(float)), cudaMemcpyDeviceToHost);
           
   //do final reduction on host
   float finalsum=0.0f;
   for(i=0; i<blas_redblocks; i++){
     finalsum += blas_sredfield[i];
   }
   #ifdef MPI
     float result;
     //printf("proc %d : %f\n",g_proc_id,finalsum);
     MPI_Allreduce(&finalsum, &result, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
     finalsum=result;
   #endif
   return(finalsum);
}





// convert spinor to double 
void convert2double_spin (dev_spinor* spin, spinor* h2d) {

  int i, Vol, offset;
  
  #ifndef MPI
    if (even_odd_flag) {
      Vol = VOLUME/2;
    }
    else {
      Vol = VOLUME;
    }
  #else
   if (even_odd_flag) {
     Vol = (VOLUME+RAND)/2;
   }
   else{
     Vol = (VOLUME+RAND);
   }
  #endif
  offset = Vol;
  
  for (i = 0; i < Vol; i++) {
  
        h2d[i].s0.c0 = (double) spin[i+0*offset].x + I*(double) spin[i+0*offset].y;
        h2d[i].s0.c1 = (double) spin[i+0*offset].z + I*(double) spin[i+0*offset].w;
        
        h2d[i].s0.c2 = (double) spin[i+1*offset].x + I*(double) spin[i+1*offset].y;
        h2d[i].s1.c0 = (double) spin[i+1*offset].z + I*(double) spin[i+1*offset].w;   
        
        h2d[i].s1.c1 = (double) spin[i+2*offset].x + I*(double) spin[i+2*offset].y;
        h2d[i].s1.c2 = (double) spin[i+2*offset].z + I*(double) spin[i+2*offset].w;  
        
        h2d[i].s2.c0 = (double) spin[i+3*offset].x + I*(double) spin[i+3*offset].y;
        h2d[i].s2.c1 = (double) spin[i+3*offset].z + I*(double) spin[i+3*offset].w;  
        
        h2d[i].s2.c2 = (double) spin[i+4*offset].x + I*(double) spin[i+4*offset].y;
        h2d[i].s3.c0 = (double) spin[i+4*offset].z + I*(double) spin[i+4*offset].w; 
        
        h2d[i].s3.c1 = (double) spin[i+5*offset].x + I*(double) spin[i+5*offset].y;
        h2d[i].s3.c2 = (double) spin[i+5*offset].z + I*(double) spin[i+5*offset].w; 
        
  }
}





// convert spinor to REAL4 (float4, double4) 
void convert2REAL4_spin(spinor* spin, dev_spinor* h2d){

  int i, Vol, offset;
  
  #ifndef MPI
    if (even_odd_flag) {
      Vol = VOLUME/2;
    }
    else {
      Vol = VOLUME;
    }
  #else
   if (even_odd_flag) {
     Vol = (VOLUME+RAND)/2;
   }
   else{
     Vol = (VOLUME+RAND);
   }
  #endif
  //set the offset in fields
  offset = Vol;
  
  for (i = 0; i < Vol; i++) {
    
        h2d[i+0*offset].x = (float) creal(spin[i].s0.c0);
        h2d[i+0*offset].y = (float) cimag(spin[i].s0.c0);
        h2d[i+0*offset].z = (float) creal(spin[i].s0.c1);
        h2d[i+0*offset].w = (float) cimag(spin[i].s0.c1);
        
        h2d[i+1*offset].x = (float) creal(spin[i].s0.c2);
        h2d[i+1*offset].y = (float) cimag(spin[i].s0.c2);
        h2d[i+1*offset].z = (float) creal(spin[i].s1.c0);
        h2d[i+1*offset].w = (float) cimag(spin[i].s1.c0);
        
        h2d[i+2*offset].x = (float) creal(spin[i].s1.c1);
        h2d[i+2*offset].y = (float) cimag(spin[i].s1.c1);
        h2d[i+2*offset].z = (float) creal(spin[i].s1.c2);
        h2d[i+2*offset].w = (float) cimag(spin[i].s1.c2);
        
        h2d[i+3*offset].x = (float) creal(spin[i].s2.c0);
        h2d[i+3*offset].y = (float) cimag(spin[i].s2.c0);
        h2d[i+3*offset].z = (float) creal(spin[i].s2.c1);
        h2d[i+3*offset].w = (float) cimag(spin[i].s2.c1);
        
        h2d[i+4*offset].x = (float) creal(spin[i].s2.c2);
        h2d[i+4*offset].y = (float) cimag(spin[i].s2.c2);
        h2d[i+4*offset].z = (float) creal(spin[i].s3.c0);
        h2d[i+4*offset].w = (float) cimag(spin[i].s3.c0);
        
        h2d[i+5*offset].x = (float) creal(spin[i].s3.c1);
        h2d[i+5*offset].y = (float) cimag(spin[i].s3.c1);
        h2d[i+5*offset].z = (float) creal(spin[i].s3.c2);
        h2d[i+5*offset].w = (float) cimag(spin[i].s3.c2);
    
  }
}




// orders a double spinor on host according to the ordering used on device 
void order_spin_gpu(spinor* spin, dev_spinor_d* h2d){

  int i, Vol, offset;
  
  #ifndef MPI
    if (even_odd_flag) {
      Vol = VOLUME/2;
    }
    else {
      Vol = VOLUME;
    }
  #else
   if (even_odd_flag) {
     Vol = (VOLUME+RAND)/2;
   }
   else{
     Vol = (VOLUME+RAND);
   }
  #endif
  //set the offset in fields
  offset = Vol;
  
  for (i = 0; i < Vol; i++) {
    
        h2d[i+0*offset].x = creal(spin[i].s0.c0);
        h2d[i+0*offset].y = cimag(spin[i].s0.c0);
        h2d[i+0*offset].z = creal(spin[i].s0.c1);
        h2d[i+0*offset].w = cimag(spin[i].s0.c1);
        
        h2d[i+1*offset].x = creal(spin[i].s0.c2);
        h2d[i+1*offset].y = cimag(spin[i].s0.c2);
        h2d[i+1*offset].z = creal(spin[i].s1.c0);
        h2d[i+1*offset].w = cimag(spin[i].s1.c0);
        
        h2d[i+2*offset].x = creal(spin[i].s1.c1);
        h2d[i+2*offset].y = cimag(spin[i].s1.c1);
        h2d[i+2*offset].z = creal(spin[i].s1.c2);
        h2d[i+2*offset].w = cimag(spin[i].s1.c2);
        
        h2d[i+3*offset].x = creal(spin[i].s2.c0);
        h2d[i+3*offset].y = cimag(spin[i].s2.c0);
        h2d[i+3*offset].z = creal(spin[i].s2.c1);
        h2d[i+3*offset].w = cimag(spin[i].s2.c1);
        
        h2d[i+4*offset].x = creal(spin[i].s2.c2);
        h2d[i+4*offset].y = cimag(spin[i].s2.c2);
        h2d[i+4*offset].z = creal(spin[i].s3.c0);
        h2d[i+4*offset].w = cimag(spin[i].s3.c0);
        
        h2d[i+5*offset].x = creal(spin[i].s3.c1);
        h2d[i+5*offset].y = cimag(spin[i].s3.c1);
        h2d[i+5*offset].z = creal(spin[i].s3.c2);
        h2d[i+5*offset].w = cimag(spin[i].s3.c2);
    
  }
}


// orders a double spinor on host according to the ordering used on host 
void unorder_spin_gpu (dev_spinor_d* spin, spinor* h2d) {

  int i, Vol, offset;
  
  #ifndef MPI
    if (even_odd_flag) {
      Vol = VOLUME/2;
    }
    else {
      Vol = VOLUME;
    }
  #else
   if (even_odd_flag) {
     Vol = (VOLUME+RAND)/2;
   }
   else{
     Vol = (VOLUME+RAND);
   }
  #endif
  offset = Vol;
  
  for (i = 0; i < Vol; i++) {
  
        h2d[i].s0.c0 = spin[i+0*offset].x + I* spin[i+0*offset].y;
        h2d[i].s0.c1 = spin[i+0*offset].z + I* spin[i+0*offset].w;
        
        h2d[i].s0.c2 = spin[i+1*offset].x + I* spin[i+1*offset].y;
        h2d[i].s1.c0 = spin[i+1*offset].z + I* spin[i+1*offset].w;   
        
        h2d[i].s1.c1 = spin[i+2*offset].x + I* spin[i+2*offset].y;
        h2d[i].s1.c2 = spin[i+2*offset].z + I* spin[i+2*offset].w;  
        
        h2d[i].s2.c0 = spin[i+3*offset].x + I* spin[i+3*offset].y;
        h2d[i].s2.c1 = spin[i+3*offset].z + I* spin[i+3*offset].w;  
        
        h2d[i].s2.c2 = spin[i+4*offset].x + I* spin[i+4*offset].y;
        h2d[i].s3.c0 = spin[i+4*offset].z + I* spin[i+4*offset].w; 
        
        h2d[i].s3.c1 = spin[i+5*offset].x + I* spin[i+5*offset].y;
        h2d[i].s3.c2 = spin[i+5*offset].z + I* spin[i+5*offset].w; 
        
  }
}



__global__ void dev_zero_spinor_field(dev_spinor* s1){
  int pos;
  pos= threadIdx.x + blockDim.x*blockIdx.x;  
  if(pos < dev_VOLUME){
          dev_zero_spinor(&(s1[pos]));
  }
}




__global__ void dev_copy_spinor_field(dev_spinor* s1, dev_spinor* s2){
    int pos;
  pos= threadIdx.x + blockDim.x*blockIdx.x;  
  if(pos < dev_VOLUME){
      dev_copy_spinor(&(s1[pos]),&(s2[pos]));
  } 
}



__global__ void dev_skalarmult_add_assign_spinor_field(dev_spinor* s1, float lambda, dev_spinor* s2, dev_spinor* so){
  int pos;
  pos= threadIdx.x + blockDim.x*blockIdx.x;  
  if(pos < dev_VOLUME){
    dev_skalarmult_add_assign_spinor(&(s1[pos]), lambda ,&(s2[pos]), &(so[pos]) );
  }
}



__global__ void dev_skalarmult_spinor_field(dev_spinor* s1, float lambda, dev_spinor* so){
  int pos;
  pos= threadIdx.x + blockDim.x*blockIdx.x;  
  if(pos < dev_VOLUME){
    dev_skalarmult_spinor(&(s1[pos]), dev_initcomplex(lambda,0.0) , &(so[pos]) );
  }
}  



__global__ void dev_complexmult_spinor_field(dev_spinor* s1, dev_complex lambda, dev_spinor* so){
  int pos;
  pos= threadIdx.x + blockDim.x*blockIdx.x;  
  if(pos < dev_VOLUME){
    dev_skalarmult_spinor(&(s1[pos]), lambda , &(so[pos]) );
  }
}





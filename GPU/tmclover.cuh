/***********************************************************************
 *
 * Copyright (C) 2013 Florian Burger
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
 * File: tmclover.cuh
 *
 * CUDA GPU kernels for tm+clover
 *
 **************************************************************************/


extern "C" {
#include "../operator/clovertm_operators.h"
} 

//the clover fields
float2 * dev_sw;
float2 * h2d_sw;
float2 * dev_sw_inv;
float2 * h2d_sw_inv;







// load sw_inv field
// pos is site index
// a is first index
// b is second index
// numbers can be fetched from texture cache float2 field (one float2 == complex number)
// adjacent threads read linearly, as pos is fastest index
// difference btw even and odd (+ and - mu) is managed by pos
// pos needs to be set pos+8*2*VOLUME*9
// we need 9 float2 to store a 3x3 complex matrix
__device__ void dev_load_sw_inv (dev_su3* gf, int pos, int a, int b, int VOL, const float2* field){
  float2 gfin;
  
  int ASIZE = 4;
  int SWSIZE = 9;


  //first row
  #ifdef USETEXTURE
    gfin = tex1Dfetch(sw_inv_tex,pos+VOL*(0+SWSIZE*(a + ASIZE*b)));
  #else
    gfin = field[pos+VOL*(0+SWSIZE*(a + ASIZE*b))];
  #endif
  (*gf)[0][0].re = gfin.x;
  (*gf)[0][0].im = gfin.y;

  #ifdef USETEXTURE
    gfin = tex1Dfetch(sw_inv_tex,pos+VOL*(1+SWSIZE*(a + ASIZE*b)));
  #else
    gfin = field[pos+VOL*(1+SWSIZE*(a + ASIZE*b))];
  #endif
  (*gf)[0][1].re = gfin.x;
  (*gf)[0][1].im = gfin.y;

  #ifdef USETEXTURE
    gfin = tex1Dfetch(sw_inv_tex,pos+VOL*(2+SWSIZE*(a + ASIZE*b)));
  #else
    gfin = field[pos+VOL*(2+SWSIZE*(a + ASIZE*b))];
  #endif
  (*gf)[0][2].re = gfin.x;
  (*gf)[0][2].im = gfin.y;
  //*******




  //second row
  #ifdef USETEXTURE
    gfin = tex1Dfetch(sw_inv_tex,pos+VOL*(3+SWSIZE*(a + ASIZE*b)));
  #else
    gfin = field[pos+VOL*(3+SWSIZE*(a + ASIZE*b))];
  #endif
  (*gf)[1][0].re = gfin.x;
  (*gf)[1][0].im = gfin.y;

  #ifdef USETEXTURE
    gfin = tex1Dfetch(sw_inv_tex,pos+VOL*(4+SWSIZE*(a + ASIZE*b)));
  #else
    gfin = field[pos+VOL*(4+SWSIZE*(a + ASIZE*b))];
  #endif
  (*gf)[1][1].re = gfin.x;
  (*gf)[1][1].im = gfin.y;

  #ifdef USETEXTURE
    gfin = tex1Dfetch(sw_inv_tex,pos+VOL*(5+SWSIZE*(a + ASIZE*b)));
  #else
    gfin = field[pos+VOL*(5+SWSIZE*(a + ASIZE*b))];
  #endif
  (*gf)[1][2].re = gfin.x;
  (*gf)[1][2].im = gfin.y;
  //*******  




  //third row
  #ifdef USETEXTURE
    gfin = tex1Dfetch(sw_inv_tex,pos+VOL*(6+SWSIZE*(a + ASIZE*b)));
  #else
    gfin = field[pos+VOL*(6+SWSIZE*(a + ASIZE*b))];
  #endif
  (*gf)[2][0].re = gfin.x;
  (*gf)[2][0].im = gfin.y;

  #ifdef USETEXTURE
    gfin = tex1Dfetch(sw_inv_tex,pos+VOL*(7+SWSIZE*(a + ASIZE*b)));
  #else
    gfin = field[pos+VOL*(7+SWSIZE*(a + ASIZE*b))];
  #endif
  (*gf)[2][1].re = gfin.x;
  (*gf)[2][1].im = gfin.y;

  #ifdef USETEXTURE
    gfin = tex1Dfetch(sw_inv_tex,pos+VOL*(8+SWSIZE*(a + ASIZE*b)));
  #else
    gfin = field[pos+VOL*(8+SWSIZE*(a + ASIZE*b))];
  #endif
  (*gf)[2][2].re = gfin.x;
  (*gf)[2][2].im = gfin.y;
  //*******  

  return;
}


//double version
__device__ void dev_load_sw_inv_d (dev_su3_d* gf, int pos, int a, int b, int VOL, const double2* field){
  double2 gfin;
  
  int ASIZE = 4;
  int SWSIZE = 9;


  //first row
  gfin = field[pos+VOL*(0+SWSIZE*(a + ASIZE*b))];
  (*gf)[0][0].re = gfin.x;
  (*gf)[0][0].im = gfin.y;

  gfin = field[pos+VOL*(1+SWSIZE*(a + ASIZE*b))];
  (*gf)[0][1].re = gfin.x;
  (*gf)[0][1].im = gfin.y;


  gfin = field[pos+VOL*(2+SWSIZE*(a + ASIZE*b))];
  (*gf)[0][2].re = gfin.x;
  (*gf)[0][2].im = gfin.y;
  //*******


  //second row
  gfin = field[pos+VOL*(3+SWSIZE*(a + ASIZE*b))];
  (*gf)[1][0].re = gfin.x;
  (*gf)[1][0].im = gfin.y;

  gfin = field[pos+VOL*(4+SWSIZE*(a + ASIZE*b))];
  (*gf)[1][1].re = gfin.x;
  (*gf)[1][1].im = gfin.y;

  gfin = field[pos+VOL*(5+SWSIZE*(a + ASIZE*b))];
  (*gf)[1][2].re = gfin.x;
  (*gf)[1][2].im = gfin.y;
  //*******  


  //third row
  gfin = field[pos+VOL*(6+SWSIZE*(a + ASIZE*b))];
  (*gf)[2][0].re = gfin.x;
  (*gf)[2][0].im = gfin.y;

  gfin = field[pos+VOL*(7+SWSIZE*(a + ASIZE*b))];
  (*gf)[2][1].re = gfin.x;
  (*gf)[2][1].im = gfin.y;

  gfin = field[pos+VOL*(8+SWSIZE*(a + ASIZE*b))];
  (*gf)[2][2].re = gfin.x;
  (*gf)[2][2].im = gfin.y;
  //*******  

  return;
}






// load sw field
// pos is site index
// a is first index
// b is second index
// numbers can be fetched from texture cache float2 field (one float2 == complex number)
// adjacent threads read linearly, as pos is fastest index
// difference btw even and odd (+ and - mu) is managed by pos
// pos needs to be set pos+6*2*VOLUME*9
// we need 9 float2 to store a 3x3 complex matrix
__device__ void dev_load_sw(dev_su3* gf, int pos, int a, int b, int VOL, const float2* field ){
  float2 gfin;
  
  int ASIZE = 3;
  int SWSIZE = 9;



  //first row
  #ifdef USETEXTURE
    gfin = tex1Dfetch(sw_tex,pos+VOL*(0+SWSIZE*(a + ASIZE*b)));
  #else
    gfin = field[pos+VOL*(0+SWSIZE*(a + ASIZE*b))];
  #endif
  (*gf)[0][0].re = gfin.x;
  (*gf)[0][0].im = gfin.y;

  #ifdef USETEXTURE
    gfin = tex1Dfetch(sw_tex,pos+VOL*(1+SWSIZE*(a + ASIZE*b)));
  #else
    gfin = field[pos+VOL*(1+SWSIZE*(a + ASIZE*b))];
  #endif
  (*gf)[0][1].re = gfin.x;
  (*gf)[0][1].im = gfin.y;

  #ifdef USETEXTURE
    gfin = tex1Dfetch(sw_tex,pos+VOL*(2+SWSIZE*(a + ASIZE*b)));
  #else
    gfin = field[pos+VOL*(2+SWSIZE*(a + ASIZE*b))];
  #endif
  (*gf)[0][2].re = gfin.x;
  (*gf)[0][2].im = gfin.y;
  //*******




  //second row
  #ifdef USETEXTURE
    gfin = tex1Dfetch(sw_tex,pos+VOL*(3+SWSIZE*(a + ASIZE*b)));
  #else
    gfin = field[pos+VOL*(3+SWSIZE*(a + ASIZE*b))];
  #endif
  (*gf)[1][0].re = gfin.x;
  (*gf)[1][0].im = gfin.y;

  #ifdef USETEXTURE
    gfin = tex1Dfetch(sw_tex,pos+VOL*(4+SWSIZE*(a + ASIZE*b)));
  #else
    gfin = field[pos+VOL*(4+SWSIZE*(a + ASIZE*b))];
  #endif
  (*gf)[1][1].re = gfin.x;
  (*gf)[1][1].im = gfin.y;

  #ifdef USETEXTURE
    gfin = tex1Dfetch(sw_tex,pos+VOL*(5+SWSIZE*(a + ASIZE*b)));
  #else
    gfin = field[pos+VOL*(5+SWSIZE*(a + ASIZE*b))];
  #endif
  (*gf)[1][2].re = gfin.x;
  (*gf)[1][2].im = gfin.y;
  //*******  




  //third row
  #ifdef USETEXTURE
    gfin = tex1Dfetch(sw_tex,pos+VOL*(6+SWSIZE*(a + ASIZE*b)));
  #else
    gfin = field[pos+VOL*(6+SWSIZE*(a + ASIZE*b))];
  #endif
  (*gf)[2][0].re = gfin.x;
  (*gf)[2][0].im = gfin.y;

  #ifdef USETEXTURE
    gfin = tex1Dfetch(sw_tex,pos+VOL*(7+SWSIZE*(a + ASIZE*b)));
  #else
    gfin = field[pos+VOL*(7+SWSIZE*(a + ASIZE*b))];
  #endif
  (*gf)[2][1].re = gfin.x;
  (*gf)[2][1].im = gfin.y;

  #ifdef USETEXTURE
    gfin = tex1Dfetch(sw_tex,pos+VOL*(8+SWSIZE*(a + ASIZE*b)));
  #else
    gfin = field[pos+VOL*(8+SWSIZE*(a + ASIZE*b))];
  #endif
  (*gf)[2][2].re = gfin.x;
  (*gf)[2][2].im = gfin.y;
  //*******  

  return;
}





__device__ void dev_load_sw_d(dev_su3_d* gf, int pos, int a, int b, int VOL, const double2* field ){
  double2 gfin;
  
  int ASIZE = 3;
  int SWSIZE = 9;



  //first row
  gfin = field[pos+VOL*(0+SWSIZE*(a + ASIZE*b))];
  (*gf)[0][0].re = gfin.x;
  (*gf)[0][0].im = gfin.y;

  gfin = field[pos+VOL*(1+SWSIZE*(a + ASIZE*b))];
  (*gf)[0][1].re = gfin.x;
  (*gf)[0][1].im = gfin.y;

  gfin = field[pos+VOL*(2+SWSIZE*(a + ASIZE*b))];
  (*gf)[0][2].re = gfin.x;
  (*gf)[0][2].im = gfin.y;
  //*******


  //second row
  gfin = field[pos+VOL*(3+SWSIZE*(a + ASIZE*b))];
  (*gf)[1][0].re = gfin.x;
  (*gf)[1][0].im = gfin.y;

  gfin = field[pos+VOL*(4+SWSIZE*(a + ASIZE*b))];
  (*gf)[1][1].re = gfin.x;
  (*gf)[1][1].im = gfin.y;

  gfin = field[pos+VOL*(5+SWSIZE*(a + ASIZE*b))];
  (*gf)[1][2].re = gfin.x;
  (*gf)[1][2].im = gfin.y;
  //*******  




  //third row
  gfin = field[pos+VOL*(6+SWSIZE*(a + ASIZE*b))];
  (*gf)[2][0].re = gfin.x;
  (*gf)[2][0].im = gfin.y;

  gfin = field[pos+VOL*(7+SWSIZE*(a + ASIZE*b))];
  (*gf)[2][1].re = gfin.x;
  (*gf)[2][1].im = gfin.y;

  gfin = field[pos+VOL*(8+SWSIZE*(a + ASIZE*b))];
  (*gf)[2][2].re = gfin.x;
  (*gf)[2][2].im = gfin.y;
  //*******  

  return;
}







// *************************  -> linalg.cuh *****************************

__device__ inline void dev_get_su3_vec0(dev_vector out, float4* in){

   out[0].re = in[0].x;  out[0].im = in[0].y;
   out[1].re = in[0].z;  out[1].im = in[0].w;
   out[2].re = in[1].x;  out[2].im = in[1].y;
}



__device__ inline void dev_get_su3_vec1(dev_vector out, float4* in){

   out[0].re = in[1].z;  out[0].im = in[1].w;
   out[1].re = in[2].x;  out[1].im = in[2].y;
   out[2].re = in[2].z;  out[2].im = in[2].w;
}


__device__ inline void dev_get_su3_vec2(dev_vector out, float4* in){

   out[0].re = in[3].x;  out[0].im = in[3].y;
   out[1].re = in[3].z;  out[1].im = in[3].w;
   out[2].re = in[4].x;  out[2].im = in[4].y;
}


__device__ inline void dev_get_su3_vec3(dev_vector out, float4* in){

   out[0].re = in[4].z;  out[0].im = in[4].w;
   out[1].re = in[5].x;  out[1].im = in[5].y;
   out[2].re = in[5].z;  out[2].im = in[5].w;
}



// *************************  -> linalg_d.cuh *****************************

__device__ inline void dev_get_su3_vec0_d(dev_vector_d out, double4* in){

   out[0].re = in[0].x;  out[0].im = in[0].y;
   out[1].re = in[0].z;  out[1].im = in[0].w;
   out[2].re = in[1].x;  out[2].im = in[1].y;
}



__device__ inline void dev_get_su3_vec1_d(dev_vector_d out, double4* in){

   out[0].re = in[1].z;  out[0].im = in[1].w;
   out[1].re = in[2].x;  out[1].im = in[2].y;
   out[2].re = in[2].z;  out[2].im = in[2].w;
}


__device__ inline void dev_get_su3_vec2_d(dev_vector_d out, double4* in){

   out[0].re = in[3].x;  out[0].im = in[3].y;
   out[1].re = in[3].z;  out[1].im = in[3].w;
   out[2].re = in[4].x;  out[2].im = in[4].y;
}


__device__ inline void dev_get_su3_vec3_d(dev_vector_d out, double4* in){

   out[0].re = in[4].z;  out[0].im = in[4].w;
   out[1].re = in[5].x;  out[1].im = in[5].y;
   out[2].re = in[5].z;  out[2].im = in[5].w;
}






__device__ inline void dev_store_spinor_from_vec(dev_spinor* out, 
                                                 dev_vector in0, dev_vector in1,
                                                 dev_vector in2, dev_vector in3){
   out[0].x = in0[0].re; 
   out[0].y = in0[0].im;
   out[0].z = in0[1].re; 
   out[0].w = in0[1].im;
   out[1].x = in0[2].re; 
   out[1].y = in0[2].im;

   out[1].z = in1[0].re; 
   out[1].w = in1[0].im;
   out[2].x = in1[1].re; 
   out[2].y = in1[1].im;
   out[2].z = in1[2].re; 
   out[2].w = in1[2].im;

   out[3].x = in2[0].re; 
   out[3].y = in2[0].im;
   out[3].z = in2[1].re; 
   out[3].w = in2[1].im;
   out[4].x = in2[2].re; 
   out[4].y = in2[2].im;

   out[4].z = in3[0].re; 
   out[4].w = in3[0].im;
   out[5].x = in3[1].re; 
   out[5].y = in3[1].im;
   out[5].z = in3[2].re; 
   out[5].w = in3[2].im;

}



__device__ inline void dev_su3vec_add(dev_vector out, dev_vector in1, dev_vector in2){

   out[0].re  = in1[0].re;
   out[0].re += in2[0].re;
   out[0].im  = in1[0].im;
   out[0].im += in2[0].im;
   out[1].re  = in1[1].re;
   out[1].re += in2[1].re;
   out[1].im  = in1[1].im;
   out[1].im += in2[1].im;
   out[2].re  = in1[2].re;
   out[2].re += in2[2].re;
   out[2].im  = in1[2].im;
   out[2].im += in2[2].im;

}



/* out += I * c * in (c real) */
__device__ inline void dev_su3vec_add_i_mul_old(dev_vector out, dev_complex c, dev_vector in1){

   out[0].re += c.re*in1[0].re;
   out[0].re -= c.im*in1[0].im;
   out[0].im += c.re*in1[0].im;
   out[0].im += c.im*in1[0].re;

   out[1].re += c.re*in1[1].re;
   out[1].re -= c.im*in1[1].im;
   out[1].im += c.re*in1[1].im;
   out[1].im += c.im*in1[1].re;

   out[2].re += c.re*in1[2].re;
   out[2].re -= c.im*in1[2].im;
   out[2].im += c.re*in1[2].im;
   out[2].im += c.im*in1[2].re;
}


__device__ inline void dev_su3vec_add_i_mul(dev_vector out, float c, dev_vector in1){
  
   out[0].re -= c*in1[0].im; 
   out[0].im += c*in1[0].re;

   out[1].re -= c*in1[1].im;
   out[1].im += c*in1[1].re;

   out[2].re -= c*in1[2].im;
   out[2].im += c*in1[2].re;

}


__device__ inline void dev_su3vec_add_assign(dev_vector out, dev_vector in1){

   out[0].re += in1[0].re;
   out[0].im += in1[0].im;
   out[1].re += in1[1].re;
   out[1].im += in1[1].im;
   out[2].re += in1[2].re;
   out[2].im += in1[2].im;

}



__device__ inline void dev_su3vec_sub(dev_vector out, dev_vector in1, dev_vector in2){

   out[0].re  = in1[0].re;
   out[0].re -= in2[0].re;
   out[0].im  = in1[0].im;
   out[0].im -= in2[0].im;
   out[1].re  = in1[1].re;
   out[1].re -= in2[1].re;
   out[1].im  = in1[1].im;
   out[1].im -= in2[1].im;
   out[2].re  = in1[2].re;
   out[2].re -= in2[2].re;
   out[2].im  = in1[2].im;
   out[2].im -= in2[2].im;

}



__device__ inline void dev_su3_multiply(dev_vector r, dev_su3 u, dev_vector s){					        
   r[0].re=   u[0][0].re*s[0].re-u[0][0].im*s[0].im  
             +u[0][1].re*s[1].re-u[0][1].im*s[1].im  
             +u[0][2].re*s[2].re-u[0][2].im*s[2].im; 
   r[0].im=   u[0][0].re*s[0].im+u[0][0].im*s[0].re  
             +u[0][1].re*s[1].im+u[0][1].im*s[1].re  
             +u[0][2].re*s[2].im+u[0][2].im*s[2].re; 
   r[1].re=   u[1][0].re*s[0].re-u[1][0].im*s[0].im  
             +u[1][1].re*s[1].re-u[1][1].im*s[1].im  
             +u[1][2].re*s[2].re-u[1][2].im*s[2].im; 
   r[1].im=   u[1][0].re*s[0].im+u[1][0].im*s[0].re  
             +u[1][1].re*s[1].im+u[1][1].im*s[1].re  
             +u[1][2].re*s[2].im+u[1][2].im*s[2].re; 
   r[2].re=   u[2][0].re*s[0].re-u[2][0].im*s[0].im  
             +u[2][1].re*s[1].re-u[2][1].im*s[1].im  
             +u[2][2].re*s[2].re-u[2][2].im*s[2].im; 
   r[2].im=   u[2][0].re*s[0].im+u[2][0].im*s[0].re  
             +u[2][1].re*s[1].im+u[2][1].im*s[1].re  
             +u[2][2].re*s[2].im+u[2][2].im*s[2].re;
}


__device__ inline void dev_su3_inverse_multiply(dev_vector r, dev_su3 u, dev_vector s){	

   r[0].re= u[0][0].re*s[0].re+u[0][0].im*s[0].im  
           +u[1][0].re*s[1].re+u[1][0].im*s[1].im  
           +u[2][0].re*s[2].re+u[2][0].im*s[2].im; 
   r[0].im= u[0][0].re*s[0].im-u[0][0].im*s[0].re  
           +u[1][0].re*s[1].im-u[1][0].im*s[1].re  
           +u[2][0].re*s[2].im-u[2][0].im*s[2].re; 
   r[1].re= u[0][1].re*s[0].re+u[0][1].im*s[0].im  
           +u[1][1].re*s[1].re+u[1][1].im*s[1].im  
           +u[2][1].re*s[2].re+u[2][1].im*s[2].im; 
   r[1].im= u[0][1].re*s[0].im-u[0][1].im*s[0].re  
           +u[1][1].re*s[1].im-u[1][1].im*s[1].re  
           +u[2][1].re*s[2].im-u[2][1].im*s[2].re; 
   r[2].re= u[0][2].re*s[0].re+u[0][2].im*s[0].im  
           +u[1][2].re*s[1].re+u[1][2].im*s[1].im  
           +u[2][2].re*s[2].re+u[2][2].im*s[2].im; 
   r[2].im= u[0][2].re*s[0].im-u[0][2].im*s[0].re  
           +u[1][2].re*s[1].im-u[1][2].im*s[1].re  
           +u[2][2].re*s[2].im-u[2][2].im*s[2].re;
}





// !!!!*************************  -> linalg.cuh *****************************





void order_sw_gpu(float2* h2d){
  int i,a,b;
  int VOL = VOLUME;
  int SWSIZE = 9;  
  int ASIZE = 3;
  int pos=0;
  for(i=0; i<VOLUME; i++){
    for(a=0;a<3;a++){
      for(b=0;b<2;b++){
        h2d[pos+VOL*(0+SWSIZE*(a + ASIZE*b))].x = (float)creal(sw[pos][a][b].c00);
        h2d[pos+VOL*(0+SWSIZE*(a + ASIZE*b))].y = (float)cimag(sw[pos][a][b].c00);
        h2d[pos+VOL*(1+SWSIZE*(a + ASIZE*b))].x = (float)creal(sw[pos][a][b].c01);
        h2d[pos+VOL*(1+SWSIZE*(a + ASIZE*b))].y = (float)cimag(sw[pos][a][b].c01);
        h2d[pos+VOL*(2+SWSIZE*(a + ASIZE*b))].x = (float)creal(sw[pos][a][b].c02);
        h2d[pos+VOL*(2+SWSIZE*(a + ASIZE*b))].y = (float)cimag(sw[pos][a][b].c02);
	
        h2d[pos+VOL*(3+SWSIZE*(a + ASIZE*b))].x = (float)creal(sw[pos][a][b].c10);
        h2d[pos+VOL*(3+SWSIZE*(a + ASIZE*b))].y = (float)cimag(sw[pos][a][b].c10);
        h2d[pos+VOL*(4+SWSIZE*(a + ASIZE*b))].x = (float)creal(sw[pos][a][b].c11);
        h2d[pos+VOL*(4+SWSIZE*(a + ASIZE*b))].y = (float)cimag(sw[pos][a][b].c11);
        h2d[pos+VOL*(5+SWSIZE*(a + ASIZE*b))].x = (float)creal(sw[pos][a][b].c12);
        h2d[pos+VOL*(5+SWSIZE*(a + ASIZE*b))].y = (float)cimag(sw[pos][a][b].c12);	
	
        h2d[pos+VOL*(6+SWSIZE*(a + ASIZE*b))].x = (float)creal(sw[pos][a][b].c20);
        h2d[pos+VOL*(6+SWSIZE*(a + ASIZE*b))].y = (float)cimag(sw[pos][a][b].c20);
        h2d[pos+VOL*(7+SWSIZE*(a + ASIZE*b))].x = (float)creal(sw[pos][a][b].c21);
        h2d[pos+VOL*(7+SWSIZE*(a + ASIZE*b))].y = (float)cimag(sw[pos][a][b].c21);
        h2d[pos+VOL*(8+SWSIZE*(a + ASIZE*b))].x = (float)creal(sw[pos][a][b].c22);
        h2d[pos+VOL*(8+SWSIZE*(a + ASIZE*b))].y = (float)cimag(sw[pos][a][b].c22);	
      }
    }
    pos++;
  }
}


void order_sw_inv_gpu(float2* h2d){
  int i,a,b;
  int VOL = VOLUME;
  int SWSIZE = 9;  
  int ASIZE = 4;
  int pos=0;
  for(i=0; i<VOLUME; i++){
    for(a=0;a<4;a++){
      for(b=0;b<2;b++){
        h2d[pos+VOL*(0+SWSIZE*(a + ASIZE*b))].x = (float)creal(sw_inv[pos][a][b].c00);
        h2d[pos+VOL*(0+SWSIZE*(a + ASIZE*b))].y = (float)cimag(sw_inv[pos][a][b].c00);
        h2d[pos+VOL*(1+SWSIZE*(a + ASIZE*b))].x = (float)creal(sw_inv[pos][a][b].c01);
        h2d[pos+VOL*(1+SWSIZE*(a + ASIZE*b))].y = (float)cimag(sw_inv[pos][a][b].c01);
        h2d[pos+VOL*(2+SWSIZE*(a + ASIZE*b))].x = (float)creal(sw_inv[pos][a][b].c02);
        h2d[pos+VOL*(2+SWSIZE*(a + ASIZE*b))].y = (float)cimag(sw_inv[pos][a][b].c02);
	
        h2d[pos+VOL*(3+SWSIZE*(a + ASIZE*b))].x = (float)creal(sw_inv[pos][a][b].c10);
        h2d[pos+VOL*(3+SWSIZE*(a + ASIZE*b))].y = (float)cimag(sw_inv[pos][a][b].c10);
        h2d[pos+VOL*(4+SWSIZE*(a + ASIZE*b))].x = (float)creal(sw_inv[pos][a][b].c11);
        h2d[pos+VOL*(4+SWSIZE*(a + ASIZE*b))].y = (float)cimag(sw_inv[pos][a][b].c11);
        h2d[pos+VOL*(5+SWSIZE*(a + ASIZE*b))].x = (float)creal(sw_inv[pos][a][b].c12);
        h2d[pos+VOL*(5+SWSIZE*(a + ASIZE*b))].y = (float)cimag(sw_inv[pos][a][b].c12);	
	
        h2d[pos+VOL*(6+SWSIZE*(a + ASIZE*b))].x = (float)creal(sw_inv[pos][a][b].c20);
        h2d[pos+VOL*(6+SWSIZE*(a + ASIZE*b))].y = (float)cimag(sw_inv[pos][a][b].c20);
        h2d[pos+VOL*(7+SWSIZE*(a + ASIZE*b))].x = (float)creal(sw_inv[pos][a][b].c21);
        h2d[pos+VOL*(7+SWSIZE*(a + ASIZE*b))].y = (float)cimag(sw_inv[pos][a][b].c21);
        h2d[pos+VOL*(8+SWSIZE*(a + ASIZE*b))].x = (float)creal(sw_inv[pos][a][b].c22);
        h2d[pos+VOL*(8+SWSIZE*(a + ASIZE*b))].y = (float)cimag(sw_inv[pos][a][b].c22);	
      }
    }
    pos++;
  }
}


//assumes that sw and sw_inv are initialized on host!!
//we are performing a reordering of the fields in which
//spacetime is the fastest index
void init_gpu_clover_fields(){

  cudaError_t cudaerr;
  if(g_cart_id == 0) printf("Initializing clover fields on gpu...\n");
  //sw - field
  size_t sw_size = 6*VOLUME*9*sizeof(float2); /* float2 */
  if((void*)(h2d_sw = (float2*)malloc(sw_size)) == NULL){
    if(g_cart_id == 0) printf("Could not allocate memory for h2d_sw. Aborting...\n");
    exit(201);
  } // Allocate conversion sw field on host 
  
  if((cudaerr=cudaMalloc((void **) &dev_sw, sw_size)) != cudaSuccess){
    if(g_cart_id == 0) printf("Error in init_gpu_clover_fields(): Memory allocation of sw field failed with error %d. Aborting...\n", cudaerr);
    exit(201);
  }   // Allocate array on device 
  
  if(g_cart_id == 0) printf("Rearranging sw field on gpu...\n");
  order_sw_gpu(h2d_sw);
  cudaMemcpy(dev_sw, h2d_sw, sw_size, cudaMemcpyHostToDevice);
  
  //sw_inv - field
  size_t sw_inv_size = 8*VOLUME*9*sizeof(float2); /* float2 */
  if((void*)(h2d_sw_inv = (float2*)malloc(sw_inv_size)) == NULL){
    if(g_cart_id == 0) printf("Could not allocate memory for h2d_sw_inv. Aborting...\n");
    exit(201);
  } // Allocate conversion sw_inv field on host 
  
  if((cudaerr=cudaMalloc((void **) &dev_sw_inv, sw_inv_size)) != cudaSuccess){
    if(g_cart_id == 0) printf("Error in init_gpu_clover_fields(): Memory allocation of sw field failed with error %d. Aborting...\n", cudaerr);
    exit(201);
  }   // Allocate array on device 
  
  if(g_cart_id == 0) printf("Rearranging sw_inv field on gpu...\n");  
  order_sw_inv_gpu(h2d_sw_inv);
  cudaMemcpy(dev_sw_inv, h2d_sw_inv, sw_inv_size, cudaMemcpyHostToDevice);   

  
  if(g_cart_id == 0) printf("Finished clover gpu initialisation\n");   
}



void finalize_gpu_clover_fields(){
  free(h2d_sw);
  free(h2d_sw_inv);
  cudaFree(dev_sw);
  cudaFree(dev_sw_inv); 
}











// computes same as clover_inv
// sin is read and stored after mult with sw_inv
__global__ void dev_clover_inv(dev_spinor* sin,  float2* sw_inv, const float mu_in){
   
   dev_spinor slocal[6];
   dev_vector psi, chi, phi0, phi1, phi2, phi3, r0, r1, r2, r3;
    __shared__ dev_su3 mat[BLOCK];
   int offset = 0;
   //dev_VOLUME is VOLUME/2 for EO, therefor no /2 here!!
   if(mu_in < 0.0) offset = dev_VOLUME;  


   int pos, ix;
   pos= threadIdx.x + blockDim.x*blockIdx.x;  
   ix= threadIdx.x;
   if(pos < dev_VOLUME){

     dev_read_spinor(&(slocal[0]), &(sin[pos]));
#ifdef RELATIVISTIC_BASIS 
     //rotate to tmlqcd basis, as clover term is defined in that basis
     to_tmlqcd_basis_spinor(&(slocal[0])); 
#endif     
     dev_get_su3_vec0(phi0, &(slocal[0]));
     dev_get_su3_vec1(phi1, &(slocal[0]));
     dev_get_su3_vec2(phi2, &(slocal[0]));
     dev_get_su3_vec3(phi3, &(slocal[0]));

/*
    _vector_assign(phi1,(*rn).s0);
    _vector_assign(phi3,(*rn).s2);

    w1=&sw_inv[icy][0][0];
    w2=w1+2;  // &sw_inv[icy][1][0]; 
    w3=w1+4;  // &sw_inv[icy][2][0]; 
    w4=w1+6;  // &sw_inv[icy][3][0]; 
    _su3_multiply(psi,*w1,phi1); 
    _su3_multiply(chi,*w2,(*rn).s1);
    _vector_add((*rn).s0,psi,chi);
    _su3_multiply(psi,*w4,phi1); 
    _su3_multiply(chi,*w3,(*rn).s1);
    _vector_add((*rn).s1,psi,chi);
*/

     //sw_inv[?][0]
     dev_load_sw_inv (&(mat[ix]), (offset+pos), 0, 0, 2*dev_VOLUME, sw_inv);
     dev_su3_multiply(psi,mat[ix],phi0);

     dev_load_sw_inv (&(mat[ix]), (offset+pos), 1, 0, 2*dev_VOLUME, sw_inv);
     dev_su3_multiply(chi,mat[ix],phi1);
     dev_su3vec_add(r0,psi,chi);

     dev_load_sw_inv (&(mat[ix]), (offset+pos), 3, 0, 2*dev_VOLUME, sw_inv);
     dev_su3_multiply(psi,mat[ix],phi0);

     dev_load_sw_inv (&(mat[ix]), (offset+pos), 2, 0, 2*dev_VOLUME, sw_inv);
     dev_su3_multiply(chi,mat[ix],phi1);
     dev_su3vec_add(r1,psi,chi);

/*
    w1++; // &sw_inv[icy][0][1]; 
    w2++; // &sw_inv[icy][1][1]; 
    w3++; // &sw_inv[icy][2][1]; 
    w4++; // &sw_inv[icy][3][1]; 
    _su3_multiply(psi,*w1,phi3); 
    _su3_multiply(chi,*w2,(*rn).s3);
    _vector_add((*rn).s2,psi,chi);
    _su3_multiply(psi,*w4,phi3); 
    _su3_multiply(chi,*w3,(*rn).s3);
    _vector_add((*rn).s3,psi,chi);
*/

     //sw_inv[?][1]
     dev_load_sw_inv (&(mat[ix]), (offset+pos), 0, 1, 2*dev_VOLUME, sw_inv);
     dev_su3_multiply(psi,mat[ix],phi2);

     dev_load_sw_inv (&(mat[ix]), (offset+pos), 1, 1, 2*dev_VOLUME, sw_inv);
     dev_su3_multiply(chi,mat[ix],phi3);
     dev_su3vec_add(r2,psi,chi);

     dev_load_sw_inv (&(mat[ix]), (offset+pos), 3, 1, 2*dev_VOLUME, sw_inv);
     dev_su3_multiply(psi,mat[ix],phi2);

     dev_load_sw_inv (&(mat[ix]), (offset+pos), 2, 1, 2*dev_VOLUME, sw_inv);
     dev_su3_multiply(chi,mat[ix],phi3);
     dev_su3vec_add(r3,psi,chi);

     dev_store_spinor_from_vec(&(slocal[0]), r0, r1, r2, r3);
     
#ifdef RELATIVISTIC_BASIS 
     //rotate to tmlqcd basis, as clover term is defined in that basis
     to_relativistic_basis_spinor(&(slocal[0])); 
#endif      
     
     dev_write_spinor(&(slocal[0]),&(sin[pos]));    

   }
}





// computes same as clover_gamma5
// j and k are read and l is stored
__global__ void dev_clover_gamma5(const int ieo, const int * index_site, dev_spinor* l, dev_spinor* k ,dev_spinor* j, float2* sw, const float mu_in){
   
   dev_spinor slocal[6];
   //rho? is output
   dev_vector psi1, psi2, chi, phi0, phi1, phi2, phi3,  
                         tau0, tau1, tau2, tau3,
                         rho0, rho1, rho2, rho3;                         
   __shared__ dev_su3 mat[BLOCK];

   int pos, ix;
   
   pos= threadIdx.x + blockDim.x*blockIdx.x;  
   ix= threadIdx.x; 
   
   if(pos < dev_VOLUME){

/*  
    r = l + icx-ioff;
    s = k + icx-ioff;
    t = j + icx-ioff;
    
    w1=&sw[ix][0][0];
    w2=w1+2; //&sw[ix][1][0];
    w3=w1+4; //&sw[ix][2][0];
    _su3_multiply(psi1,*w1,(*s).s0); 
    _su3_multiply(chi,*w2,(*s).s1);
    _vector_add_assign(psi1,chi);
    _su3_inverse_multiply(psi2,*w2,(*s).s0); 
    _su3_multiply(chi,*w3,(*s).s1);
    _vector_add_assign(psi2,chi); 
    // add in the twisted mass term (plus in the upper components)
    _vector_add_i_mul(psi1, mu, (*s).s0);
    _vector_add_i_mul(psi2, mu, (*s).s1);

    _vector_sub((*r).s0,psi1,(*t).s0);
    _vector_sub((*r).s1,psi2,(*t).s1);
*/
     //load k
     dev_read_spinor(&(slocal[0]), &(k[pos]));
     
#ifdef RELATIVISTIC_BASIS 
     //rotate to tmlqcd basis, as clover term is defined in that basis
     to_tmlqcd_basis_spinor(&(slocal[0])); 
#endif        
     
     dev_get_su3_vec0(phi0, &(slocal[0]));
     dev_get_su3_vec1(phi1, &(slocal[0]));
     dev_get_su3_vec2(phi2, &(slocal[0]));
     dev_get_su3_vec3(phi3, &(slocal[0]));

     dev_load_sw (&(mat[ix]), index_site[pos], 0, 0, 2*dev_VOLUME, sw);
     dev_su3_multiply(psi1,mat[ix],phi0);
     
     dev_load_sw (&(mat[ix]), index_site[pos], 1, 0, 2*dev_VOLUME, sw);
     dev_su3_multiply(chi,mat[ix],phi1);
     dev_su3vec_add_assign(psi1,chi);
     dev_su3_inverse_multiply(psi2,mat[ix],phi0); 

     dev_load_sw (&(mat[ix]), index_site[pos], 2, 0, 2*dev_VOLUME, sw);     
     dev_su3_multiply(chi,mat[ix],phi1);
     dev_su3vec_add_assign(psi2,chi);

     dev_su3vec_add_i_mul(psi1, mu_in, phi0);
     dev_su3vec_add_i_mul(psi2, mu_in, phi1);
    
     //load j
     dev_read_spinor(&(slocal[0]), &(j[pos]));
     
#ifdef RELATIVISTIC_BASIS 
     //rotate to tmlqcd basis, as clover term is defined in that basis
     to_tmlqcd_basis_spinor(&(slocal[0])); 
#endif         
     
     dev_get_su3_vec0(tau0, &(slocal[0]));
     dev_get_su3_vec1(tau1, &(slocal[0]));
     dev_get_su3_vec2(tau2, &(slocal[0]));
     dev_get_su3_vec3(tau3, &(slocal[0]));

     dev_su3vec_sub(rho0, psi1, tau0);
     dev_su3vec_sub(rho1, psi2, tau1);

/*
    w1++; //=&sw[ix][0][1];
    w2++; //=&sw[ix][1][1];
    w3++; //=&sw[ix][2][1];
    _su3_multiply(psi1,*w1,(*s).s2); 
    _su3_multiply(chi,*w2,(*s).s3);
    _vector_add_assign(psi1,chi); 
    _su3_inverse_multiply(psi2,*w2,(*s).s2); 
    _su3_multiply(chi,*w3,(*s).s3);
    _vector_add_assign(psi2,chi); 
    // add in the twisted mass term (minus from g5 in the lower components)
    _vector_add_i_mul(psi1, -mu, (*s).s2);
    _vector_add_i_mul(psi2, -mu, (*s).s3);

    // *************** multiply with  gamma5 included *****************************
    _vector_sub((*r).s2,(*t).s2,psi1);
    _vector_sub((*r).s3,(*t).s3,psi2);
    // ******************************* end of loop ********************************
*/

     dev_load_sw (&(mat[ix]), index_site[pos], 0, 1, 2*dev_VOLUME, sw);
     dev_su3_multiply(psi1,mat[ix],phi2);
     dev_load_sw (&(mat[ix]), index_site[pos], 1, 1, 2*dev_VOLUME, sw);
     dev_su3_multiply(chi,mat[ix],phi3);
     dev_su3vec_add_assign(psi1,chi);
     dev_su3_inverse_multiply(psi2,mat[ix],phi2); 

     dev_load_sw (&(mat[ix]), index_site[pos], 2, 1, 2*dev_VOLUME, sw);
     dev_su3_multiply(chi,mat[ix],phi3);
     dev_su3vec_add_assign(psi2,chi);

     dev_su3vec_add_i_mul(psi1, -mu_in, phi2);
     dev_su3vec_add_i_mul(psi2, -mu_in, phi3);

     dev_su3vec_sub(rho2, tau2, psi1);
     dev_su3vec_sub(rho3, tau3, psi2);

     //write to l
     dev_store_spinor_from_vec(&(slocal[0]), rho0, rho1, rho2, rho3);
     
 #ifdef RELATIVISTIC_BASIS 
     //rotate to tmlqcd basis, as clover term is defined in that basis
     to_relativistic_basis_spinor(&(slocal[0])); 
#endif       
     
     dev_write_spinor(&(slocal[0]),&(l[pos]));  

   }
}







// aequivalent to Qsw_pm_psi in clovertm_operators.c, this is NON-MPI version
extern "C" void dev_Qsw_pm_psi(dev_spinor* spinin, dev_spinor* spinout){
  //spinin == odd
  //spinout == odd
  #ifndef _USE_MPI
    int VolumeEO = VOLUME/2;
  #endif
    
  #ifdef USETEXTURE
    bind_texture_sw(dev_sw);
    bind_texture_sw_inv(dev_sw_inv);  
  #endif
    
  //Q_{-}
  #ifdef USETEXTURE
    bind_texture_spin(spinin,1);
  #endif
    
    
  #ifdef _USE_MPI
    HOPPING_ASYNC(dev_gf, spinin, dev_spin_eo1, 
		  dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 
		  0, gpu_gd_M, gpu_bd_M); //dev_spin_eo1 == even -> 0     
  #else
    dev_Hopping_Matrix<<<gpu_gd_M, gpu_bd_M>>>
             (dev_gf, spinin, dev_spin_eo1, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0, 0, VolumeEO); //dev_spin_eo1 == even -> 0           
  #endif 
	     
  #ifdef USETEXTURE
    unbind_texture_spin(1);
  #endif
  dev_clover_inv<<<gpu_gd_M, gpu_bd_M>>>(dev_spin_eo1, dev_sw_inv, -1.0);
  

  #ifdef USETEXTURE
    bind_texture_spin(dev_spin_eo1,1);
  #endif

   #ifdef _USE_MPI
    HOPPING_ASYNC(dev_gf, dev_spin_eo1, dev_spin_eo2, 
		  dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 
		  1, gpu_gd_M, gpu_bd_M);     
   #else
    dev_Hopping_Matrix<<<gpu_gd_M, gpu_bd_M>>>
            (dev_gf, dev_spin_eo1, dev_spin_eo2, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1, 0, VolumeEO); 
   #endif
  #ifdef USETEXTURE
    unbind_texture_spin(1);
  #endif
  dev_clover_gamma5<<<gpu_gd_M, gpu_bd_M>>>(1,dev_eoidx_odd, dev_spin_eo2, spinin, dev_spin_eo2,  dev_sw, (float)(-(g_mu+g_mu3))); 
  
  //Q_{+}

  #ifdef USETEXTURE
    bind_texture_spin(dev_spin_eo2,1);
  #endif

  #ifdef _USE_MPI
    HOPPING_ASYNC(dev_gf, dev_spin_eo2, spinout, 
		  dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 
		  0, gpu_gd_M, gpu_bd_M); //dev_spin_eo1 == even -> 0    
  #else
    dev_Hopping_Matrix<<<gpu_gd_M, gpu_bd_M>>>
          (dev_gf, dev_spin_eo2, spinout, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0, 0, VolumeEO); //dev_spin_eo1 == even -> 0
  #endif     
  #ifdef USETEXTURE  
    unbind_texture_spin(1);
  #endif
  dev_clover_inv<<<gpu_gd_M, gpu_bd_M>>>(spinout, dev_sw_inv, +1.0);
  

  #ifdef USETEXTURE
    bind_texture_spin(spinout,1);
  #endif

  #ifdef _USE_MPI
    HOPPING_ASYNC(dev_gf, spinout, dev_spin_eo1, 
		  dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 
		  1, gpu_gd_M, gpu_bd_M);     
  #else
    dev_Hopping_Matrix<<<gpu_gd_M, gpu_bd_M>>>
             (dev_gf, spinout, dev_spin_eo1, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1, 0, VolumeEO); 
  #endif  
  #ifdef USETEXTURE
    unbind_texture_spin(1);
  #endif
  dev_clover_gamma5<<<gpu_gd_M, gpu_bd_M>>>(1,dev_eoidx_odd, spinout, dev_spin_eo2, dev_spin_eo1,  dev_sw, (float)(+(g_mu+g_mu3))); 
 
  #ifdef USETEXTURE
    unbind_texture_sw();
    unbind_texture_sw_inv();  
  #endif
    
}

















extern "C" void clover_gamma5(const int ieo, 
		   spinor * const l, const spinor * const k, const spinor * const j,
		   const double mu);




__global__ void check_sw_reconst(float2* sw, dev_su3* to, int vol){
   
   int pos;
   pos= threadIdx.x + blockDim.x*blockIdx.x;  

   if(pos < 1){
     dev_load_sw(to, 1, 1, 0, 2*dev_VOLUME, sw);
   }
}











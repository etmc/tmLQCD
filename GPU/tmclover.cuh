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
extern "C" void dev_Qsw_pm_psi(dev_spinor* spinin, dev_spinor* spinout, int gridsize, dim3 blocksize, int gridsize2, int blocksize2){
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
		  0, gridsize, blocksize); //dev_spin_eo1 == even -> 0     
  #else
    dev_Hopping_Matrix<<<gridsize, blocksize>>>
             (dev_gf, spinin, dev_spin_eo1, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0, 0, VolumeEO); //dev_spin_eo1 == even -> 0           
  #endif 
	     
  #ifdef USETEXTURE
    unbind_texture_spin(1);
  #endif
  dev_clover_inv<<<gridsize, blocksize>>>(dev_spin_eo1, dev_sw_inv, -1.0);
  

  #ifdef USETEXTURE
    bind_texture_spin(dev_spin_eo1,1);
  #endif

   #ifdef _USE_MPI
    HOPPING_ASYNC(dev_gf, dev_spin_eo1, dev_spin_eo2, 
		  dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 
		  1, gridsize, blocksize);     
   #else
    dev_Hopping_Matrix<<<gridsize, blocksize>>>
            (dev_gf, dev_spin_eo1, dev_spin_eo2, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1, 0, VolumeEO); 
   #endif
  #ifdef USETEXTURE
    unbind_texture_spin(1);
  #endif
  dev_clover_gamma5<<<gridsize, blocksize>>>(1,dev_eoidx_odd, dev_spin_eo2, spinin, dev_spin_eo2,  dev_sw, (float)(-(g_mu+g_mu3))); 
  
  //Q_{+}

  #ifdef USETEXTURE
    bind_texture_spin(dev_spin_eo2,1);
  #endif

  #ifdef _USE_MPI
    HOPPING_ASYNC(dev_gf, dev_spin_eo2, spinout, 
		  dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 
		  0, gridsize, blocksize); //dev_spin_eo1 == even -> 0    
  #else
    dev_Hopping_Matrix<<<gridsize, blocksize>>>
          (dev_gf, dev_spin_eo2, spinout, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0, 0, VolumeEO); //dev_spin_eo1 == even -> 0
  #endif     
  #ifdef USETEXTURE  
    unbind_texture_spin(1);
  #endif
  dev_clover_inv<<<gridsize, blocksize>>>(spinout, dev_sw_inv, +1.0);
  

  #ifdef USETEXTURE
    bind_texture_spin(spinout,1);
  #endif

  #ifdef _USE_MPI
    HOPPING_ASYNC(dev_gf, spinout, dev_spin_eo1, 
		  dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 
		  1, gridsize, blocksize);     
  #else
    dev_Hopping_Matrix<<<gridsize, blocksize>>>
             (dev_gf, spinout, dev_spin_eo1, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1, 0, VolumeEO); 
  #endif  
  #ifdef USETEXTURE
    unbind_texture_spin(1);
  #endif
  dev_clover_gamma5<<<gridsize, blocksize>>>(1,dev_eoidx_odd, spinout, dev_spin_eo2, dev_spin_eo1,  dev_sw, (float)(+(g_mu+g_mu3))); 
 
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






/*FIXME have to get rid of all that comes below */
void update_constants(int *grid);

	   



void test_clover_operator(spinor* const Q, const int N){
   
   size_t dev_spinsize = 6*VOLUME/2 * sizeof(dev_spinor); // float4 even-odd !   
   int gridsize;
     //this is the partitioning for the HoppingMatrix kernel
     int blockdim3 = BLOCK;
     if( VOLUME/2 % blockdim3 == 0){
       gridsize = (int) VOLUME/2/blockdim3;
     }
     else{
       gridsize = (int) VOLUME/2/blockdim3 + 1;
     }
     int griddim3 = gridsize;
   
     //this is the partitioning for dev_mul_one_pm...
     int blockdim4 = BLOCK2;
     if( VOLUME/2 % blockdim4 == 0){
       gridsize = (int) VOLUME/2/blockdim4;
     }
     else{
       gridsize = (int) VOLUME/2/blockdim4 + 1;
     }   
     
  int grid[6];
  grid[0]=LX; grid[1]=LY; grid[2]=LZ; grid[3]=T; grid[4]=VOLUME/2; 
  // dev_VOLUME is half of VOLUME for eo
  
  // put dev_Offset accordingly depending on mpi/non-mpi
  #ifdef _USE_MPI
   grid[5] = (VOLUME+RAND)/2;
  #else
   grid[5] = VOLUME/2;
  #endif
  cudaMalloc((void **) &dev_grid, 6*sizeof(int));
  cudaMemcpy(dev_grid, &(grid[0]), 6*sizeof(int), cudaMemcpyHostToDevice);
  update_constants(dev_grid);    
  
  
  #ifdef USETEXTURE
    bind_texture_sw(dev_sw);
    bind_texture_sw_inv(dev_sw_inv);  
  #endif
  
//   size_t idxsize = VOLUME/2*sizeof(int);
//   size_t nnsize = 8*VOLUME*sizeof(int);
//   initnn();
//   initnn_eo();   
//   cudaMemcpy(dev_nn, nn, nnsize, cudaMemcpyHostToDevice);
//   cudaMemcpy(dev_nn_eo, nn_eo, nnsize/2, cudaMemcpyHostToDevice);
//   cudaMemcpy(dev_nn_oe, nn_oe, nnsize/2, cudaMemcpyHostToDevice);
//   cudaMemcpy(dev_eoidx_even, eoidx_even, idxsize, cudaMemcpyHostToDevice);
//   cudaMemcpy(dev_eoidx_odd, eoidx_odd, idxsize, cudaMemcpyHostToDevice);
//    
        if (g_cart_id == 0) {
  	  int host_check_VOLUMEPLUSRAND;
  	  int host_check_VOLUME;
  	  int host_check_Offset;
  	  cudaMemcpyFromSymbol(&host_check_VOLUMEPLUSRAND, dev_VOLUMEPLUSRAND, sizeof(int));
  	  cudaMemcpyFromSymbol(&host_check_VOLUME, dev_VOLUME, sizeof(int));
  	  cudaMemcpyFromSymbol(&host_check_Offset, dev_Offset, sizeof(int));
  	  printf("\tOn device:\n");
  	  printf("\tdev_VOLUMEPLUSRAND = %i\n", host_check_VOLUMEPLUSRAND);
  	  printf("\tdev_VOLUME = %i\n", host_check_VOLUME);
  	  printf("\tdev_Offset = %i\n", host_check_Offset);
  	} 
  
  int test[20],i ;
  cudaMemcpy(&test[0], dev_eoidx_even, 20*sizeof(int), cudaMemcpyDeviceToHost);     
  for(i=0; i<20; i++) printf("%d ", test[i]);
  printf("\n");
 
     
  spinor ** solver_field = NULL;
  const int nr_sf = 3;
  init_solver_field(&solver_field, VOLUMEPLUSRAND/2, nr_sf);  

  //apply cpu matrix
  ////complete
  //Qsw_pm_psi(solver_field[0], Q);  

  ////clover_inv
  //assign(solver_field[0], Q, N);
  //clover_inv(EE,solver_field[0], -g_mu); 

  ////clover_gamma5 
  clover_gamma5(OO, solver_field[0], Q, Q, +(g_mu + g_mu3));

  
  //apply gpu matrix
  convert2REAL4_spin(Q,h2d_spin);
  cudaMemcpy(dev_spin2, h2d_spin, dev_spinsize, cudaMemcpyHostToDevice);

  ////complete
  //dev_Qsw_pm_psi(dev_spin2, dev_spin1,griddim3, blockdim3, griddim4, blockdim4); 
  //cudaMemcpy(h2d_spin, dev_spin1, dev_spinsize, cudaMemcpyDeviceToHost);

  ////clover_inv
  //dev_clover_inv<<<griddim3, blockdim3>>>(dev_spin2, dev_sw_inv, -1.0);
  //cudaMemcpy(h2d_spin, dev_spin2, dev_spinsize, cudaMemcpyDeviceToHost);

  ////clover_gamma5
  dev_clover_gamma5<<<griddim3, blockdim3>>>(1,dev_eoidx_odd, dev_spin1, dev_spin2, dev_spin2,  dev_sw, (float)(+(g_mu+g_mu3))); 
  cudaMemcpy(h2d_spin, dev_spin1, dev_spinsize, cudaMemcpyDeviceToHost);  
  
  printf("%s\n", cudaGetErrorString(cudaGetLastError()));
  convert2double_spin(h2d_spin, solver_field[1]);

  diff(solver_field[2], solver_field[1], solver_field[0],N);
  double rk = square_norm(solver_field[2], N, 1);
    
  printf("%.6e %.6e\n", creal(solver_field[2][0].s0.c0), cimag(solver_field[2][0].s0.c0));
  printf("%.6e %.6e\n", creal(solver_field[2][0].s0.c1), cimag(solver_field[2][0].s0.c1));  
  printf("%.6e %.6e\n", creal(solver_field[2][0].s0.c2), cimag(solver_field[2][0].s0.c2)); 
  
  printf("%.6e %.6e\n", creal(solver_field[2][0].s1.c0), cimag(solver_field[2][0].s1.c0));
  printf("%.6e %.6e\n", creal(solver_field[2][0].s1.c1), cimag(solver_field[2][0].s1.c1));  
  printf("%.6e %.6e\n", creal(solver_field[2][0].s1.c2), cimag(solver_field[2][0].s1.c2));  
  
  printf("%.6e %.6e\n", creal(solver_field[2][0].s2.c0), cimag(solver_field[2][0].s2.c0));
  printf("%.6e %.6e\n", creal(solver_field[2][0].s2.c1), cimag(solver_field[2][0].s2.c1));  
  printf("%.6e %.6e\n", creal(solver_field[2][0].s2.c2), cimag(solver_field[2][0].s2.c2));  
  
  printf("%.6e %.6e\n", creal(solver_field[2][0].s3.c0), cimag(solver_field[2][0].s3.c0));
  printf("%.6e %.6e\n", creal(solver_field[2][0].s3.c1), cimag(solver_field[2][0].s3.c1));  
  printf("%.6e %.6e\n", creal(solver_field[2][0].s3.c2), cimag(solver_field[2][0].s3.c2));   
    
  
  printf("\n\n");
  
  
  printf("Testing clover matrix:\n");
  printf("cpu: Squared difference is: %.8e\n", rk);
  printf("cpu: Squared difference per spinor component is: %.8e\n\n\n", rk/N/24.0);  
  
  
  printf("Checking clover first entry...\n");
  float first_real;
  cudaMemcpy(&first_real, dev_sw, sizeof(float), cudaMemcpyDeviceToHost);     
  printf("%.8f\n\n", first_real);  

  
    
  
  printf("Checking clover reconstruction...\n");
  float hostfield[18];
  check_sw_reconst<<<1,1>>>(dev_sw, (dev_su3*) dev_spin2, VOLUME); 
  cudaMemcpy(&hostfield, dev_spin2, 18*sizeof(float), cudaMemcpyDeviceToHost);    


  
  printf("On device:\n");  
  for(int i=0; i<18; i++) printf("%.6e ", hostfield[i]);

  printf("\n");   
  printf("On host:\n");  
  printf("%.6e %.6e\n", creal(sw[1][1][0].c00), cimag(sw[1][1][0].c00));
  printf("%.6e %.6e\n", creal(sw[1][1][0].c01), cimag(sw[1][1][0].c01));  
  printf("%.6e %.6e\n", creal(sw[1][1][0].c02), cimag(sw[1][1][0].c02));
  
  printf("%.6e %.6e\n", creal(sw[1][1][0].c10), cimag(sw[1][1][0].c10));
  printf("%.6e %.6e\n", creal(sw[1][1][0].c11), cimag(sw[1][1][0].c11));  
  printf("%.6e %.6e\n", creal(sw[1][1][0].c12), cimag(sw[1][1][0].c12)); 
  
  printf("%.6e %.6e\n", creal(sw[1][1][0].c20), cimag(sw[1][1][0].c20));
  printf("%.6e %.6e\n", creal(sw[1][1][0].c21), cimag(sw[1][1][0].c21));  
  printf("%.6e %.6e\n", creal(sw[1][1][0].c22), cimag(sw[1][1][0].c22));   
  
  
  
  
    
  #ifdef USETEXTURE
    unbind_texture_sw();
    unbind_texture_sw_inv();  
  #endif
  
  
  finalize_solver(solver_field, nr_sf);  
}
















/*FIXME this must be merged with tm eo innver solver and made transparent */
// this is the eo version of the device cg inner solver 
// we invert the hermitean Qsw_pm_psi
extern "C" int dev_cg_eo_clover(
      dev_su3_2v * gf,
      dev_spinor* spinin, 
      dev_spinor* spinout, 
      dev_spinor* spin0, 
      dev_spinor* spin1, 
      dev_spinor* spin2, 
      dev_spinor* spin3, 
      dev_spinor* spin4, 
      int *grid, int * nn_grid, float epsfinal){
 
 
 float host_alpha, host_beta, host_dotprod, host_rk, sourcesquarenorm;
 float * dotprod, * dotprod2, * rk, * alpha, *beta;
 
 
 
 int i, gridsize;
 int maxit = max_innersolver_it;
 float eps = (float) innersolver_precision;
 int N_recalcres = 40; // after N_recalcres iterations calculate r = A x_k - b
 
 cudaError_t cudaerr;
 
 // this is the partitioning for the copying of fields
 dim3 blockdim(1,1);
 //dim3 blockdim2(128,1,1);
 
 int blockdim2 = BLOCK3;
 if( VOLUME/2 % blockdim2 == 0){
   gridsize = (int) VOLUME/2/blockdim2;
 }
 else{
   gridsize = (int) VOLUME/2/blockdim2 + 1;
 }
 int griddim2 = gridsize;

 
 //this is the partitioning for the HoppingMatrix kernel
  dim3 blockdim3(BLOCK);
  int blocksize = BLOCK;   
  if( VOLUME/2 % blocksize == 0){
    gridsize = (int) VOLUME/2/blocksize;
  }
  else{
    gridsize = (int) VOLUME/2/blocksize + 1;
  }
  int griddim3 = gridsize;


 
 //this is the partitioning for dev_mul_one_pm...
 int blockdim4 = BLOCK2;
 if( VOLUME/2 % blockdim4 == 0){
   gridsize = (int) VOLUME/2/blockdim4;
 }
 else{
   gridsize = (int) VOLUME/2/blockdim4 + 1;
 }
 int griddim4 = gridsize;
 
 #ifndef LOWOUTPUT
    if (g_proc_id == 0) {
      printf("griddim3 = %d\n", griddim3);
      printf("blockdim3 = %d\n", blockdim3.x);
      printf("griddim4 = %d\n", griddim4); 
      printf("blockdim4 = %d\n", blockdim4);          
    }
 #endif
 
 
 
 //Initialize some stuff
    //if (g_proc_id == 0) printf("mu = %f\n", g_mu);

  update_constants(grid);
  
  #ifdef _USE_MPI
    he_cg_init_nd_additional_mpi<<<1,1>>>(VOLUMEPLUSRAND/2, RAND, g_cart_id, g_nproc);
    // debug	// check dev_VOLUMEPLUSRAND and dev_RAND on device
    #ifndef LOWOUTPUT
        if (g_proc_id == 0) {
  	  int host_check_VOLUMEPLUSRAND, host_check_RAND;
  	  int host_check_rank, host_check_nproc;
  	  int host_check_VOLUME;
  	  int host_check_Offset;
  	  cudaMemcpyFromSymbol(&host_check_VOLUMEPLUSRAND, dev_VOLUMEPLUSRAND, sizeof(int));
  	  cudaMemcpyFromSymbol(&host_check_RAND, dev_RAND, sizeof(int));
  	  cudaMemcpyFromSymbol(&host_check_VOLUME, dev_VOLUME, sizeof(int));
  	  cudaMemcpyFromSymbol(&host_check_Offset, dev_Offset, sizeof(int));
  	  printf("\tOn device:\n");
  	  printf("\tdev_VOLUMEPLUSRAND = %i\n", host_check_VOLUMEPLUSRAND);
  	  printf("\tdev_VOLUME = %i\n", host_check_VOLUME);
  	  printf("\tdev_Offset = %i\n", host_check_Offset);
  	  printf("\tdev_RAND = %i\n", host_check_RAND);
  	  cudaMemcpyFromSymbol(&host_check_rank, dev_rank, sizeof(int));
  	  cudaMemcpyFromSymbol(&host_check_nproc, dev_nproc, sizeof(int));
  	  printf("\tdev_rank = %i\n", host_check_rank);
  	  printf("\tdev_nproc = %i\n", host_check_nproc);
  	}
    #endif
  #endif
  
  
  #ifdef USETEXTURE
    //Bind texture gf
    bind_texture_gf(gf);
  #endif
 
 
 // Init x,p,r for k=0
 // Allocate some numbers for host <-> device interaction
 cudaMalloc((void **) &dotprod, sizeof(float));
 cudaMalloc((void **) &dotprod2, sizeof(float));
 cudaMalloc((void **) &rk, sizeof(float));
 cudaMalloc((void **) &alpha, sizeof(float));
 cudaMalloc((void **) &beta, sizeof(float));
 #ifndef LOWOUTPUT 
 if ((cudaerr=cudaGetLastError())!=cudaSuccess) {
   if (g_proc_id == 0) printf("%s\n", cudaGetErrorString(cudaGetLastError())); 
 }
 #endif
 //init blas
 start_blas(VOLUME/2);

    if (g_proc_id == 0) {
      if ((cudaerr=cudaGetLastError())!=cudaSuccess) {
        printf("%s\n", cudaGetErrorString(cudaGetLastError())); 
      }
      #ifndef LOWOUTPUT 
        printf("have initialized cublas\n"); 
      #endif
    }

 
 #ifdef RELATIVISTIC_BASIS 
   //transform to relativistic gamma basis
   to_relativistic_basis<<<griddim4, blockdim4>>> (spinin);

   if ((cudaerr=cudaGetLastError())!=cudaSuccess) {
     if (g_proc_id == 0) printf("%s\n", cudaGetErrorString(cudaGetLastError()));
   }
   else{
     #ifndef LOWOUTPUT 
     if (g_proc_id == 0) printf("Switched to relativistic basis\n");
     #endif
   }
 #endif
 

 dev_copy_spinor_field<<<griddim2, blockdim2 >>>(spinin, spin0);
 dev_zero_spinor_field<<<griddim2, blockdim2 >>>(spin1); // x_0 = 0
 dev_copy_spinor_field<<<griddim2, blockdim2 >>>(spinin, spin2);
 dev_zero_spinor_field<<<griddim2, blockdim2 >>>(spin3);
  
   if ((cudaerr=cudaGetLastError())!=cudaSuccess) {
     if (g_proc_id == 0) printf("%s\n", cudaGetErrorString(cudaGetLastError()));
   }
 

 //relative precision -> get initial residue
 #ifndef _USE_MPI
   sourcesquarenorm = cublasSdot (24*VOLUME/2, (const float *)spinin, 1, (const float *)spinin, 1);
 #else
   sourcesquarenorm = float_dotprod(spinin, spinin,24*VOLUME/2);
 #endif
 host_rk = sourcesquarenorm; //for use in main loop
 


    if (g_proc_id == 0) {
      #ifndef LOWOUTPUT
        printf("Squarenorm Source:\t%.8e\n", sourcesquarenorm);
        printf("Entering inner solver cg-loop\n");
      #endif
      if ((cudaerr=cudaGetLastError())!=cudaSuccess) {
        printf("%s\n", cudaGetErrorString(cudaGetLastError()));
      }
    }


  #ifdef ALGORITHM_BENCHMARK
    double effectiveflops;	// will used to count the "effective" flop's (from the algorithmic perspective)
    double hoppingflops = 1608.0;
    double matrixflops  = 2*(2*(hoppingflops + 984)); //that is for dev_Qsw_pm_psi
    double starteffective;
    double stopeffective;
   // timer
   starteffective = gettime();
  #endif


 
 for(i=0;i<maxit;i++){ //MAIN LOOP
  
  // Q_{-}Q{+}
    dev_Qsw_pm_psi(spin2, spin3, griddim3, blockdim3, griddim4, blockdim4);

  
  if((cudaerr=cudaGetLastError()) != cudaSuccess){
    if (g_proc_id == 0) printf("%s\n", cudaGetErrorString(cudaerr));
    if (g_proc_id == 0) printf("Error in dev_cg_eo_clover: CUDA error after Matrix application\n", cudaGetErrorString(cudaerr));
    exit(200);
  }
  
  
 //alpha
  #ifndef _USE_MPI
    host_dotprod = cublasSdot (24*VOLUME/2, (const float *) spin2, 1, (const float *) spin3, 1);
  #else
    host_dotprod =  float_dotprod(spin2, spin3, 24*VOLUME/2);
  #endif
  
  host_alpha = (host_rk / host_dotprod); // alpha = r*r/ p M p
   
 //r(k+1)
 #ifndef _USE_MPI
   cublasSaxpy (24*VOLUME/2,-1.0*host_alpha, (const float *) spin3, 1, (float *) spin0, 1);  
 #else
   dev_axpy<<<griddim4, blockdim4>>> (-1.0*host_alpha, spin3, spin0);
 #endif 

 //x(k+1);
 #ifndef _USE_MPI 
   cublasSaxpy (24*VOLUME/2, host_alpha, (const float *) spin2,  1, (float *) spin1, 1);
 #else
   dev_axpy<<<griddim4, blockdim4>>> (host_alpha, spin2, spin1);
 #endif
 
 
  if((cudaerr=cudaGetLastError()) != cudaSuccess){
    printf("%s\n", cudaGetErrorString(cudaerr));
    exit(200);
  }

  //Abbruch?
  #ifndef _USE_MPI
    host_dotprod = cublasSdot (24*VOLUME/2, (const float *) spin0, 1,(const float *) spin0, 1);
  #else
    host_dotprod = float_dotprod(spin0, spin0, 24*VOLUME/2);
  #endif
  
 if (((host_dotprod <= eps*sourcesquarenorm)) || ( host_dotprod <= epsfinal/2.)){//error-limit erreicht (epsfinal/2 sollte ausreichen um auch in double precision zu bestehen)
   break; 
 }
  
  
    #ifndef LOWOUTPUT 
    if (g_proc_id == 0) printf("iter %d: err = %.8e\n", i, host_dotprod);
    #endif
  
  
 //beta
 host_beta =host_dotprod/host_rk;
 //p(k+1)
 #ifndef _USE_MPI
   cublasSscal (24*VOLUME/2, host_beta, (float *)spin2, 1);
   cublasSaxpy (24*VOLUME/2, 1.0, (const float *) spin0,  1, (float *) spin2, 1);
 #else
   dev_blasscal<<<griddim4, blockdim4>>> (host_beta, spin2);
   dev_axpy<<<griddim4, blockdim4>>> (1.0, spin0, spin2);
 #endif
 host_rk = host_dotprod;
 
 // recalculate residue frome r = b - Ax
 if(((i+1) % N_recalcres) == 0){
    // r_(k+1) = Ax -b 
    #ifndef LOWOUTPUT
    if (g_proc_id == 0) printf("Recalculating residue\n");
    #endif
    // D Ddagger   --   Ddagger = gamma5 D gamma5  for Wilson Dirac Operator
    // DO NOT USE tm_dirac_dagger_kappa here, otherwise spin2 will be overwritten!!!
      
    // Q_{-}Q{+}
        dev_Qsw_pm_psi(spin1, spin3, griddim3, blockdim3, griddim4, blockdim4);

    if((cudaerr=cudaGetLastError()) != cudaSuccess){
      printf("%s\n", cudaGetErrorString(cudaerr));
      exit(200);
    }  
        
    
    // r = b - Ax
    #ifndef _USE_MPI
      cublasSscal (24*VOLUME/2, -1.0, (float *)spin3, 1);
      cublasSaxpy (24*VOLUME/2, 1.0, (const float *) spinin,  1, (float *) spin3, 1);
      cublasScopy (24*VOLUME/2, (const float *)spin3, 1, (float *)spin0, 1);
    #else
      dev_blasscal<<<griddim4, blockdim4>>> (-1.0, spin3);
      dev_axpy<<<griddim4, blockdim4>>> (1.0, spinin, spin3);
      dev_blascopy<<<griddim4, blockdim4>>> (spin3, spin0);    
    #endif

   }//recalculate residue

 }//MAIN LOOP cg	
  
    
    if (g_proc_id == 0) printf("Final residue: %.6e\n",host_dotprod);
 
  #ifdef ALGORITHM_BENCHMARK
    cudaThreadSynchronize();
    stopeffective = gettime();
      // will now count the number of effective flops
      // effectiveflops  =  #(inner iterations)*(matrixflops+linalgflops)*VOLUME/2  +  #(outer iterations)*(matrixflops+linalgflops)*VOLUME/2
      // outer loop: linalg  =  flops for calculating  r(k+1) and x(k+1)
      // inner loop: linalg  =  flops for calculating  alpha, x(k+1), r(k+1), beta, d(k+1)
     #ifdef _USE_MPI
       int proccount = g_nproc;
     #else
       int proccount = 1;
     #endif 
     if(g_proc_id == 0){
      	effectiveflops = i*proccount*(matrixflops + 2*2*24 + 2*24 + 2*24 + 2*2*24 + 2*24)*VOLUME/2;
      	printf("effective BENCHMARK:\n");
      	printf("\ttotal mixed solver time:   %.4e sec\n", double(stopeffective-starteffective));
      	printf("\tfloating point operations: %.4e flops\n", effectiveflops);
      	printf("\tinner solver performance:  %.4e Gflop/s\n", double(effectiveflops) / double(stopeffective-starteffective) / 1.0e9);
     }
#endif
  
  // x_result = spin1 !
  
  //no multiplication with D^{dagger} here and no return to non-kappa basis as in dev_cg!
  dev_copy_spinor_field<<<griddim2, blockdim2 >>>(spin1,spinout);
  
  #ifdef RELATIVISTIC_BASIS 
   //transform back to tmlqcd gamma basis
   to_tmlqcd_basis<<<griddim4, blockdim4>>> (spinout);
  #endif
  
  #ifdef USETEXTURE
   unbind_texture_gf();
  #endif
  cudaFree(dotprod);
  cudaFree(dotprod2);
  cudaFree(rk);
  cudaFree(alpha);
  cudaFree(beta);
  stop_blas();

  return(i);
}













/*FIXME this must be merged with tm eo outer solver and made transparent */
extern "C" int mixed_solve_eo_clover (spinor * const P, spinor * const Q, const int max_iter, 
	   double eps, const int rel_prec, const int N){

  // source in Q, initial solution in P (not yet implemented)
  double rk;
  int outercount=0;
  int totalcount=0;
  clock_t start, stop, startinner, stopinner; 
  double timeelapsed = 0.0;
  double sourcesquarenorm;
  int iter;
  cudaError_t cudaerr;
  
  size_t dev_spinsize;
  
  #ifndef HALF
    #ifndef _USE_MPI
      dev_spinsize = 6*VOLUME/2 * sizeof(dev_spinor); // float4 even-odd !
    #else
      dev_spinsize = 6*VOLUMEPLUSRAND/2 * sizeof(dev_spinor); // float4 even-odd !
    #endif
  #else
   #ifndef _USE_MPI
    dev_spinsize = 6*VOLUME/2 * sizeof(dev_spinor_half); //short4 eo !
    size_t dev_normsize = VOLUME/2 * sizeof(float);
   #else
    dev_spinsize = 6*VOLUMEPLUSRAND/2 * sizeof(dev_spinor_half); //short4 eo !
    size_t dev_normsize = VOLUMEPLUSRAND/2 * sizeof(float);   
   #endif
  #endif  
  
  
  //update the gpu single gauge_field
  update_gpu_gf(g_gauge_field);

  
  //initialize solver fields 
  spinor ** solver_field = NULL;
  const int nr_sf = 4;
  init_solver_field(&solver_field, VOLUMEPLUSRAND/2, nr_sf);  
  //allocate solver fields non-eo!
  init_mixedsolve_fields(0);
  
  //FIXME this has to go somewhere else
  init_gpu_clover_fields();

  // FIXME implement temporal gauge
  // temporarily we fix it with a breakout test
  #ifdef TEMPORALGAUGE
    printf("Error: GPU clover and TEMPORALGAUGE (defined in GPU/cudadefs.h) not implemented.\n");
    printf("Please recompile without TEMPORALGAUGE. Aborting...\n");
    exit(-1);
  #endif
  
  //test_clover_operator(Q, N);
  
  #ifdef OPERATOR_BENCHMARK
    #ifndef HALF
    // small benchmark
      assign(solver_field[0],Q,N);
      #ifndef _USE_MPI
        benchmark_eo(solver_field[0]);
      #else
        benchmark_eo_mpi(solver_field[0]); 
      #endif
    // end small benchmark
    #endif //not HALF
  #endif
 

  // Start timer
  assert((start = clock())!=-1);
  rk = square_norm(Q, N, 1);
  sourcesquarenorm=rk; // for relative prec
  double finaleps;
  if(rel_prec == 1){
    finaleps = eps * sourcesquarenorm;
  }
  else{
    finaleps = eps;
  }


  
    assign(solver_field[0],Q,N);
    zero_spinor_field(solver_field[1],  N);//spin2 = x_k
    zero_spinor_field(solver_field[2],  N);
    
  #ifndef LOWOUTPUT
    if(g_proc_id==0) printf("Initial residue: %.16e\n",rk);
    if(g_proc_id==0) printf("The VOLUME/2 is: %d\n",N);
  #endif


for(iter=0; iter<max_iter; iter++){
   #ifndef LOWOUTPUT
   if(g_proc_id==0) printf("Applying double precision EO Dirac-Op Q_{-}Q{+}...\n");
   #endif
   
   // r_k = b - D x_k   
   Qsw_pm_psi(solver_field[3], solver_field[2]);
   diff(solver_field[0],solver_field[0],solver_field[3],N);
   rk = square_norm(solver_field[0], N, 1);
   
   #ifdef GF_8
    if(isnan(rk)){
      fprintf(stderr, "Error in mixed_solve_eo: Residue is NaN.\n  May happen with GF 8 reconstruction. Aborting ...\n");
      exit(200);
    }
   #endif
   
   if(g_proc_id==0) printf("Residue after %d inner solver iterations: %.18e\n",outercount,rk);
   
   if(((rk <= eps) && (rel_prec == 0)) || ((rk <= eps*sourcesquarenorm) && (rel_prec == 1)))
   {

     //multiply with Qtm_minus_psi (for non gpu done in invert_eo.c)
     Qsw_minus_psi(solver_field[3], solver_field[1]);
     assign(P, solver_field[3], N);

     stop = clock();
     timeelapsed = (double) (stop-start)/CLOCKS_PER_SEC;
     
     #ifndef LOWOUTPUT
      if(g_proc_id==0) printf("Reached solver precision of eps=%.2e\n",eps);
      if(g_proc_id==0) printf("EO Inversion done in mixed precision.\n Number of iterations in outer solver: %d\n Squared residue: %.8e\n Time elapsed: %.6e sec\n", outercount, rk, timeelapsed);
     #endif
     finalize_solver(solver_field, nr_sf);
     finalize_gpu_clover_fields();
     finalize_mixedsolve_fields();     
     return(totalcount);  
   }
   
  convert2REAL4_spin(solver_field[0],h2d_spin);
  cudaMemcpy(dev_spinin, h2d_spin, dev_spinsize, cudaMemcpyHostToDevice);     


  
  
  if ((cudaerr=cudaGetLastError())!=cudaSuccess) {
    printf("%s\n", cudaGetErrorString(cudaGetLastError()));
  }
   // solve in single prec on device
   // D p_k = r_k
   #ifndef LOWOUTPUT
     if(g_proc_id==0) printf("Entering inner solver\n");
   #endif
   assert((startinner = clock())!=-1);
   #ifndef HALF
      totalcount += dev_cg_eo_clover(dev_gf, dev_spinin, dev_spinout, dev_spin1, dev_spin2, dev_spin3, dev_spin4, dev_spin5, dev_grid,dev_nn, (float) finaleps);
   #else
     
     totalcount += dev_cg_eo_half(dev_gf, 
                 dev_spinin, dev_spinin_norm,
                 dev_spinout,dev_spinout_norm,
                 dev_spin1, dev_spin1_norm,
                 dev_spin2, dev_spin2_norm,
                 dev_spin3, dev_spin3_norm,
                 dev_spin4, dev_spin4_norm,
                 dev_spin5, dev_spin5_norm,
                 dev_grid,dev_nn, (float) finaleps); 


		 
   #endif
   stopinner = clock();
   timeelapsed = (double) (stopinner-startinner)/CLOCKS_PER_SEC;
   
   #ifndef LOWOUTPUT
   if(g_proc_id==0) printf("Inner solver done\nTime elapsed: %.6e sec\n", timeelapsed);

   printf("g_mu3 = %f \n", g_mu3);
   #endif
   // copy back
    
  cudaMemcpy(h2d_spin, dev_spinout, dev_spinsize, cudaMemcpyDeviceToHost);
  if ((cudaerr=cudaGetLastError())!=cudaSuccess) {
    printf("%s\n", cudaGetErrorString(cudaGetLastError()));
    printf("Error code is: %f\n",cudaerr);       
  }
  convert2double_spin(h2d_spin, solver_field[2]);


    // x_(k+1) = x_k + p_k
    add(solver_field[1],solver_field[1],solver_field[2],N);

   outercount ++;   
}// outer loop 
    
    if(g_proc_id==0) printf("Did NOT reach solver precision of eps=%.2e\n",eps);
    //multiply with Qtm_minus_psi (for non gpu done in invert_eo.c)      
    Qsw_minus_psi(solver_field[3], solver_field[1]);
    assign(P, solver_field[3], N);
    

    assert((stop = clock())!=-1);
    timeelapsed = (double) (stop-start)/CLOCKS_PER_SEC;
    if(g_proc_id==0) printf("Inversion done in mixed precision.\n Number of iterations in outer solver: %d\n Squared residue: %.8e\n Time elapsed: %.6e sec\n", outercount, rk, timeelapsed);

    finalize_solver(solver_field, nr_sf);
    finalize_gpu_clover_fields();
    finalize_mixedsolve_fields();     
  return(-1);
}







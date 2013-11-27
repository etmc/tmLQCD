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
 * File: half.cuh
 *
 * CUDA half precision conversions and BLAS kernels
 *
 * 
 *
 **************************************************************************/
 
 
#define pi_float 3.141592654f




//////// short <-> float conversion /////////////
#define SHORT_LEN 65536
//#define SCALE ((SHORT_LEN-1) * 0.5)
//#define SHIFT (-1.0f/(SHORT_LEN-1))
const float SHIFT = (-1.0f/(SHORT_LEN-1));
const float SCALE = (float)((SHORT_LEN-1) * 0.5);
const float INVSCALE = 1.0f/SCALE;


/*
#define fl2sh(f) ((short)(((f)+SHIFT)*SCALE))
#define sh2fl(s) ((float)((s)/SCALE) - SHIFT)
#define half2fl(s,norm) (norm*((float)((s)/SCALE) - SHIFT))
*/

#define fl2sh(f) ((short)(((f)+SHIFT)*SCALE))
#define sh2fl(s) ((float)((s)*INVSCALE) - SHIFT)
#define half2fl(s,norm) (norm*((float)((s)*INVSCALE) - SHIFT))

short fl2sh_host(float f) {
  short ret = (short)((f+SHIFT)*SCALE);
  return ret;
}

float sh2fl_host(short s) {
  return ((float)(s/SCALE) - SHIFT);
}


float half2float_host(short in, float innorm){
  return(sh2fl_host(in)*innorm);
}



#define construct_spinor_fromhalf(sf, sh, shn, pos){  \
   (sf)[0].x = shn*sh2fl(sh[(pos)].x);                        \
   (sf)[0].y = shn*sh2fl(sh[(pos)].y);                        \
   (sf)[0].z = shn*sh2fl(sh[(pos)].z);                        \
   (sf)[0].w = shn*sh2fl(sh[(pos)].w);                        \
   (sf)[1].x = shn*sh2fl(sh[(pos)+1*DEVOFF].x);                        \
   (sf)[1].y = shn*sh2fl(sh[(pos)+1*DEVOFF].y);                        \
   (sf)[1].z = shn*sh2fl(sh[(pos)+1*DEVOFF].z);                        \
   (sf)[1].w = shn*sh2fl(sh[(pos)+1*DEVOFF].w);                        \
   (sf)[2].x = shn*sh2fl(sh[(pos)+2*DEVOFF].x);                        \
   (sf)[2].y = shn*sh2fl(sh[(pos)+2*DEVOFF].y);                        \
   (sf)[2].z = shn*sh2fl(sh[(pos)+2*DEVOFF].z);                        \
   (sf)[2].w = shn*sh2fl(sh[(pos)+2*DEVOFF].w);                        \
   (sf)[3].x = shn*sh2fl(sh[(pos)+3*DEVOFF].x);                        \
   (sf)[3].y = shn*sh2fl(sh[(pos)+3*DEVOFF].y);                        \
   (sf)[3].z = shn*sh2fl(sh[(pos)+3*DEVOFF].z);                        \
   (sf)[3].w = shn*sh2fl(sh[(pos)+3*DEVOFF].w);                        \
   (sf)[4].x = shn*sh2fl(sh[(pos)+4*DEVOFF].x);                        \
   (sf)[4].y = shn*sh2fl(sh[(pos)+4*DEVOFF].y);                        \
   (sf)[4].z = shn*sh2fl(sh[(pos)+4*DEVOFF].z);                        \
   (sf)[4].w = shn*sh2fl(sh[(pos)+4*DEVOFF].w);                        \
   (sf)[5].x = shn*sh2fl(sh[(pos)+5*DEVOFF].x);                        \
   (sf)[5].y = shn*sh2fl(sh[(pos)+5*DEVOFF].y);                        \
   (sf)[5].z = shn*sh2fl(sh[(pos)+5*DEVOFF].z);                        \
   (sf)[5].w = shn*sh2fl(sh[(pos)+5*DEVOFF].w);                       }\
   
   
#define get_half_norm(n,s){ \
  float c0 = fmaxf(fabsf((s[0]).x), fabsf((s[0]).y));                   \
  float c1 = fmaxf(fabsf((s[0]).z), fabsf((s[0]).w));                   \
  float c2 = fmaxf(fabsf((s[1]).x), fabsf((s[1]).y));                   \
  float c3 = fmaxf(fabsf((s[1]).z), fabsf((s[1]).w));                   \
  float c4 = fmaxf(fabsf((s[2]).x), fabsf((s[2]).y));                   \
  float c5 = fmaxf(fabsf((s[2]).z), fabsf((s[2]).w));                   \
  float c6 = fmaxf(fabsf((s[3]).x), fabsf((s[3]).y));                   \
  float c7 = fmaxf(fabsf((s[3]).z), fabsf((s[3]).w));                   \
  float c8 = fmaxf(fabsf((s[4]).x), fabsf((s[4]).y));                   \
  float c9 = fmaxf(fabsf((s[4]).z), fabsf((s[4]).w));                   \
  float c10 = fmaxf(fabsf((s[5]).x), fabsf((s[5]).y));                  \
  float c11 = fmaxf(fabsf((s[5]).z), fabsf((s[5]).w));                  \
  c0 = fmaxf(c0, c1); c1 = fmaxf(c2, c3); c2 = fmaxf(c4, c5);           \
  c3 = fmaxf(c6, c7); c4 = fmaxf(c8, c9); c5 = fmaxf(c10, c11);         \
  c0 = fmaxf(c0, c1); c1 = fmaxf(c2, c3); c2 = fmaxf(c4, c5);           \
  c0 = fmaxf(c0, c1); c0 = fmaxf(c0, c2);                               \
  n = c0;                        \
}\
   
   
#define get_half_norm_from_pos(n,s,pos){                                             \
  float c0 = fmaxf(fabsf((s[pos+0*DEVOFF]).x), fabsf((s[pos+0*DEVOFF]).y));                   \
  float c1 = fmaxf(fabsf((s[pos+0*DEVOFF]).z), fabsf((s[pos+0*DEVOFF]).w));                   \
  float c2 = fmaxf(fabsf((s[pos+1*DEVOFF]).x), fabsf((s[pos+1*DEVOFF]).y));                   \
  float c3 = fmaxf(fabsf((s[pos+1*DEVOFF]).z), fabsf((s[pos+1*DEVOFF]).w));                   \
  float c4 = fmaxf(fabsf((s[pos+2*DEVOFF]).x), fabsf((s[pos+2*DEVOFF]).y));                   \
  float c5 = fmaxf(fabsf((s[pos+2*DEVOFF]).z), fabsf((s[pos+2*DEVOFF]).w));                   \
  float c6 = fmaxf(fabsf((s[pos+3*DEVOFF]).x), fabsf((s[pos+3*DEVOFF]).y));                   \
  float c7 = fmaxf(fabsf((s[pos+3*DEVOFF]).z), fabsf((s[pos+3*DEVOFF]).w));                   \
  float c8 = fmaxf(fabsf((s[pos+4*DEVOFF]).x), fabsf((s[pos+4*DEVOFF]).y));                   \
  float c9 = fmaxf(fabsf((s[pos+4*DEVOFF]).z), fabsf((s[pos+4*DEVOFF]).w));                   \
  float c10 = fmaxf(fabsf((s[pos+5*DEVOFF]).x), fabsf((s[pos+5*DEVOFF]).y));                  \
  float c11 = fmaxf(fabsf((s[pos+5*DEVOFF]).z), fabsf((s[pos+5*DEVOFF]).w));                  \
  c0 = fmaxf(c0, c1); c1 = fmaxf(c2, c3); c2 = fmaxf(c4, c5);           \
  c3 = fmaxf(c6, c7); c4 = fmaxf(c8, c9); c5 = fmaxf(c10, c11);         \
  c0 = fmaxf(c0, c1); c1 = fmaxf(c2, c3); c2 = fmaxf(c4, c5);           \
  c0 = fmaxf(c0, c1); c0 = fmaxf(c0, c2);                               \
  n = c0;                        \
}\
   
   
   
#define get_half_norm_from_pos_host(n,s,pos){                                             \
  float c0 = fmaxf(fabsf((s[pos].s0.c0.re)), fabsf((s[pos].s0.c0.im))); \
  float c1 = fmaxf(fabsf((s[pos].s0.c1.re)), fabsf((s[pos].s0.c1.im))); \
  float c2 = fmaxf(fabsf((s[pos].s0.c2.re)), fabsf((s[pos].s0.c2.im))); \
  float c3 = fmaxf(fabsf((s[pos].s1.c0.re)), fabsf((s[pos].s1.c0.im))); \
  float c4 = fmaxf(fabsf((s[pos].s1.c1.re)), fabsf((s[pos].s1.c1.im))); \
  float c5 = fmaxf(fabsf((s[pos].s1.c2.re)), fabsf((s[pos].s1.c2.im))); \
  float c6 = fmaxf(fabsf((s[pos].s2.c0.re)), fabsf((s[pos].s2.c0.im))); \
  float c7 = fmaxf(fabsf((s[pos].s2.c1.re)), fabsf((s[pos].s2.c1.im))); \
  float c8 = fmaxf(fabsf((s[pos].s2.c2.re)), fabsf((s[pos].s2.c2.im))); \
  float c9 = fmaxf(fabsf((s[pos].s3.c0.re)), fabsf((s[pos].s3.c0.im))); \
  float c10 = fmaxf(fabsf((s[pos].s3.c1.re)), fabsf((s[pos].s3.c1.im))); \
  float c11 = fmaxf(fabsf((s[pos].s3.c2.re)), fabsf((s[pos].s3.c2.im))); \
  c0 = fmaxf(c0, c1); c1 = fmaxf(c2, c3); c2 = fmaxf(c4, c5);           \
  c3 = fmaxf(c6, c7); c4 = fmaxf(c8, c9); c5 = fmaxf(c10, c11);         \
  c0 = fmaxf(c0, c1); c1 = fmaxf(c2, c3); c2 = fmaxf(c4, c5);           \
  c0 = fmaxf(c0, c1); c0 = fmaxf(c0, c2);                               \
  n = c0;                        \
}\
   
   
//////////////////////////////




// these are textures for the half spinor fields - maybe move to textures.h if possible

 /* texture for spinor field */
 texture<short4,1, cudaReadModeNormalizedFloat> spinhalf_tex0;
 texture<short4,1, cudaReadModeNormalizedFloat> spinhalf_tex1;
 texture<short4,1, cudaReadModeNormalizedFloat> spinhalf_tex2;
 texture<short4,1, cudaReadModeNormalizedFloat> spinhalf_tex3;
 texture<short4,1, cudaReadModeNormalizedFloat> spinhalf_tex4;
 texture<short4,1, cudaReadModeNormalizedFloat> spinhalf_tex5;
 
 /* texture for norm of spinor field 1*/
 texture<float,1, cudaReadModeElementType> spinnormhalf_tex;
 


  int blas_half_gridsize;
  int blas_half_blocksize; // kernel parameters for the half_dot and axpy kernels
  float * dev_blas_half_redfield;  //this is the reduction field for the
                               //blas reduction kernels
  float * dev_blas_half_sredfield; // this is the small reduction field after one sweep of reduction
  float * blas_half_sredfield;                             
  int blas_half_redblocks;  // the number of blocks of the reduction kernel
                     // VOLUME/REDUCTION_N
                     // also the size of the final sum (of reduction)
                     // performed on host
              
              
              
              
// write float spinor in to half spinor out and out_norm          
// the input spinor is a local spinor with 24 contiguous numbers
// the output spinor has to be reordered 
__device__ void dev_write_spinor_half(dev_spinor* in, dev_spinor_half* out, float* out_norm){
   float norm = 0.0f;
   int i;
   
    get_half_norm(norm, in);
   
    //store unit direction vector
    *out_norm = norm;
    float invnorm = 1.0f/norm;
    if (norm != 0.0f){
      //store norm
      #pragma unroll 6
      for(i=0; i<6; i++){
         (*(out+i*DEVOFF)).x = fl2sh(in[i].x*invnorm);
         (*(out+i*DEVOFF)).y = fl2sh(in[i].y*invnorm);
         (*(out+i*DEVOFF)).z = fl2sh(in[i].z*invnorm);
         (*(out+i*DEVOFF)).w = fl2sh(in[i].w*invnorm);
         } 
       }
      else{
       //store norm
       #pragma unroll 6
       for(i=0; i<6; i++){
         (*(out+i*DEVOFF)).x = fl2sh(0.0f);
         (*(out+i*DEVOFF)).y = fl2sh(0.0f);
         (*(out+i*DEVOFF)).z = fl2sh(0.0f);
         (*(out+i*DEVOFF)).w = fl2sh(0.0f);
       }
     }
}              
              
              
              
              
              
              

// stores the float spinor field in s into the half spinor field sh and the norm into shnorm
__global__ void float2half_spinorfield(dev_spinor_half* sh, float* shnorm, dev_spinor* s ){
  int pos=threadIdx.x + blockDim.x*blockIdx.x;
  //__shared__ float4 slocal[6];
  int i;
  float norm = 0.0;
  if(pos < dev_VOLUME){
   
     
     get_half_norm_from_pos(norm, s, pos);
     
    shnorm[pos] = norm;
    
    float invnorm = 1.0f/norm;
    //store unit direction vector
    if (norm != 0.0f){
      //store norm
      #pragma unroll 6
      for(i=0; i<6; i++){
         sh[pos+i*DEVOFF].x = fl2sh(s[6*pos+i].x*invnorm);
         sh[pos+i*DEVOFF].y = fl2sh(s[6*pos+i].y*invnorm);
         sh[pos+i*DEVOFF].z = fl2sh(s[6*pos+i].z*invnorm);
         sh[pos+i*DEVOFF].w = fl2sh(s[6*pos+i].w*invnorm);
          } 
        }
      else{
       //store norm
       #pragma unroll 6
       for(i=0; i<6; i++){
         sh[pos+i*DEVOFF].x = fl2sh(0.0f);
         sh[pos+i*DEVOFF].y = fl2sh(0.0f);
         sh[pos+i*DEVOFF].z = fl2sh(0.0f);
         sh[pos+i*DEVOFF].w = fl2sh(0.0f);
       }
     }
     
     
  }
}



// stores the half spinor field in sh, shnorm into the float spinor field s 
__global__ void half2float_spinorfield(dev_spinor* s, dev_spinor_half* sh, float* shnorm ){
  int pos=threadIdx.x + blockDim.x*blockIdx.x;
  //__shared__ float4 slocal[6];
  int i;
  float norm = 0.0;
  if(pos < dev_VOLUME){
     
    norm = shnorm[pos];
    //store unit direction vector
    if (norm != 0.0f){
      //store norm
      #pragma unroll 6
      for(i=0; i<6; i++){
         s[6*pos+i].x = sh2fl(sh[pos+i*DEVOFF].x)*norm;
         s[6*pos+i].y = sh2fl(sh[pos+i*DEVOFF].y)*norm;
         s[6*pos+i].z = sh2fl(sh[pos+i*DEVOFF].z)*norm;
         s[6*pos+i].w = sh2fl(sh[pos+i*DEVOFF].w)*norm;
          } 
        }
      else{
       //store norm
       #pragma unroll 6
       for(i=0; i<6; i++){
         s[6*pos+i].x = sh2fl(0.0f);
         s[6*pos+i].y = sh2fl(0.0f);
         s[6*pos+i].z = sh2fl(0.0f);
         s[6*pos+i].w = sh2fl(0.0f);
       }
     }
     
     
  }
}




// adds the half spinor field in sh, shnorm to the float spinor field s 
__global__ void addhalf2float_spinorfield(dev_spinor* s, dev_spinor_half* sh, float* shnorm ){
  int pos=threadIdx.x + blockDim.x*blockIdx.x;
  //__shared__ float4 slocal[6];
  int i;
  float norm = 0.0;
  if(pos < dev_VOLUME){
     
    norm = shnorm[pos];
    //store unit direction vector

      #pragma unroll 6
      for(i=0; i<6; i++){
         s[6*pos+i].x += sh2fl(sh[pos+i*DEVOFF].x)*norm;
         s[6*pos+i].y += sh2fl(sh[pos+i*DEVOFF].y)*norm;
         s[6*pos+i].z += sh2fl(sh[pos+i*DEVOFF].z)*norm;
         s[6*pos+i].w += sh2fl(sh[pos+i*DEVOFF].w)*norm;
      } 
  }     
  
}





/*
// reads half spinor from texture "spinhalf_tex" and the norm from "spinnorm_tex" and stores it into float spinor to 
__global__ void half2float_spinorfield_tex(dev_spinor* to){

int pos=threadIdx.x + blockDim.x*blockIdx.x;
  int i;
  float norm = 0.0;
  float4 help;
  if(pos < dev_VOLUME){
    norm = tex1Dfetch(spinnormhalf_tex,pos);
    for(i=0; i<6; i++){
      help = tex1Dfetch(spinhalf_tex,6*pos+i);
      to[6*pos+i].x = help.x*norm;
      to[6*pos+i].y = help.y*norm;
      to[6*pos+i].z = help.z*norm;
      to[6*pos+i].w = help.w*norm;
    }

  }
}
*/

// stores the float4 gauge field gf into the half gauge field gfh
// for GF_8 we have to be careful, as we have two angles in -Pi .. Pi
// so we have to divide them by (Pi) This is taken care of in the gauge
// reconstruction routines
// the volume is given explicitly (vol) here, to make sure alway the complete 
// VOLUME is processed and not VOLUME/2 (eo)
__global__ void float2half_gaugefield(dev_su3_2v* gf, dev_su3_2v_half* gfh, int vol){

  int pos=threadIdx.x + blockDim.x*blockIdx.x;
  int nf4,mu;
  if(pos < vol){
  for(mu=0; mu<4; mu++){
  #ifdef GF_8
    nf4 = 2;
    gfh[pos+vol*(0 + nf4*mu)].x = fl2sh(gf[pos+vol*(0 + nf4*mu)].x);
    gfh[pos+vol*(0 + nf4*mu)].y = fl2sh(gf[pos+vol*(0 + nf4*mu)].y);
    gfh[pos+vol*(0 + nf4*mu)].z = fl2sh(gf[pos+vol*(0 + nf4*mu)].z);
    gfh[pos+vol*(0 + nf4*mu)].w = fl2sh(gf[pos+vol*(0 + nf4*mu)].w);
    
    gfh[pos+vol*(1 + nf4*mu)].x = fl2sh(gf[pos+vol*(1 + nf4*mu)].x/pi_float);
    gfh[pos+vol*(1 + nf4*mu)].y = fl2sh(gf[pos+vol*(1 + nf4*mu)].y/pi_float);
    gfh[pos+vol*(1 + nf4*mu)].z = fl2sh(gf[pos+vol*(1 + nf4*mu)].z);
    gfh[pos+vol*(1 + nf4*mu)].w = fl2sh(gf[pos+vol*(1 + nf4*mu)].w);
  #else
    int i;
    nf4 = 3;
    
       gfh[pos+vol*(0 + nf4*mu)].x = fl2sh(gf[pos+vol*(0 + nf4*mu)].x);
       gfh[pos+vol*(0 + nf4*mu)].y = fl2sh(gf[pos+vol*(0 + nf4*mu)].y);
       gfh[pos+vol*(0 + nf4*mu)].z = fl2sh(gf[pos+vol*(0 + nf4*mu)].z);
       gfh[pos+vol*(0 + nf4*mu)].w = fl2sh(gf[pos+vol*(0 + nf4*mu)].w);
     
       gfh[pos+vol*(1 + nf4*mu)].x = fl2sh(gf[pos+vol*(1 + nf4*mu)].x);
       gfh[pos+vol*(1 + nf4*mu)].y = fl2sh(gf[pos+vol*(1 + nf4*mu)].y);
       gfh[pos+vol*(1 + nf4*mu)].z = fl2sh(gf[pos+vol*(1 + nf4*mu)].z);
       gfh[pos+vol*(1 + nf4*mu)].w = fl2sh(gf[pos+vol*(1 + nf4*mu)].w); 
       
       gfh[pos+vol*(2 + nf4*mu)].x = fl2sh(gf[pos+vol*(2 + nf4*mu)].x);
       gfh[pos+vol*(2 + nf4*mu)].y = fl2sh(gf[pos+vol*(2 + nf4*mu)].y);
       gfh[pos+vol*(2 + nf4*mu)].z = fl2sh(gf[pos+vol*(2 + nf4*mu)].z);
       gfh[pos+vol*(2 + nf4*mu)].w = fl2sh(gf[pos+vol*(2 + nf4*mu)].w);       
       
  #endif
  }//mu
  }
}





__device__ inline void dev_copy_spinor_half(dev_spinor_half *i1, float* i1_norm, dev_spinor_half *i2, float* i2_norm){
  int i;
  (*i2_norm) = (*i1_norm);
  #pragma unroll 6
  for(i=0;i<6;i++){ //color + spin
    (*(i2+i*DEVOFF)).x = (*(i1+i*DEVOFF)).x;
    (*(i2+i*DEVOFF)).y = (*(i1+i*DEVOFF)).y;
    (*(i2+i*DEVOFF)).z = (*(i1+i*DEVOFF)).z;
    (*(i2+i*DEVOFF)).w = (*(i1+i*DEVOFF)).w;
  }
}

__device__ inline void dev_zero_spinor_half(dev_spinor_half *sin, float* sin_norm){
  int i;
  *sin_norm = 0.0f;
  #pragma unroll 6
  for(i=0;i<6;i++){ //color + spin
    (*(sin+i*DEVOFF)).x = fl2sh(0.0f);
    (*(sin+i*DEVOFF)).y = fl2sh(0.0f);
    (*(sin+i*DEVOFF)).z = fl2sh(0.0f);
    (*(sin+i*DEVOFF)).w = fl2sh(0.0f);
  }
}





__global__ void dev_zero_spinor_field_half(
        dev_spinor_half* s1, float* s1_norm){
  int pos;
  pos= threadIdx.x + blockDim.x*blockIdx.x;  
  if(pos < dev_VOLUME){
          dev_zero_spinor_half(&(s1[pos]), &(s1_norm[pos]));
  }
}




__global__ void dev_copy_spinor_field_half(
    dev_spinor_half* s1, float* s1_norm, 
    dev_spinor_half* s2, float* s2_norm
){
    int pos;
  pos= threadIdx.x + blockDim.x*blockIdx.x;  
  if(pos < dev_VOLUME){
      dev_copy_spinor_half( &(s1[pos]), &(s1_norm[pos]),
                       &(s2[pos]), &(s2_norm[pos]));
  } 
}











extern "C" int bind_texture_spin(dev_spinor* s, int i){
 printf("Warning: DUMMY ROUTINE 'bind_texture_spin' called\n");
 return(-1);
}

extern "C" int unbind_texture_spin(int i){
  printf("Warning: DUMMY ROUTINE 'unbind_texture_spin' called\n");
  return(-1);
}


//extern "C" int bind_halfspinor_texture(dev_spinor* s, dev_spinor_half* sh, float* shnorm){
  
extern "C" int bind_halfspinor_texture(dev_spinor_half* sh, float* shnorm){
  size_t size, sizenorm;
  int gridsize, offset;
  
  #ifdef MPI
    if(even_odd_flag){
      size = sizeof(short4)*(VOLUME+RAND)/2;
      sizenorm = sizeof(float)*(VOLUME+RAND)/2;
      offset = (VOLUME+RAND)/2;
    }
    else{
      size = sizeof(short4)*(VOLUME+RAND);
      sizenorm = sizeof(float)*(VOLUME+RAND);
      offset = (VOLUME+RAND);
    }
  #else
    if(even_odd_flag){
      size = sizeof(short4)*VOLUME/2;
      sizenorm = sizeof(float)*VOLUME/2;
      offset = VOLUME/2;
    }
    else{
      size = sizeof(short4)*VOLUME;
      sizenorm = sizeof(float)*VOLUME;
      offset = VOLUME;
    }
  #endif
   /*
   // determine gridsize for conversion to half
   if( VOLUME >= BLOCK2){
      gridsize = (int) (VOLUME/BLOCK2) +1;
    }
    else{
      gridsize=1;
    }
  
     
   //DEBUG_FLO
   int i;
    
   float4 blub[6];
   for(i=0; i<6; i++){
     blub[i].x = (float) 0.1;
     blub[i].y = (float) 0.2;
     blub[i].z = (float) 0.3;
     blub[i].w = (float) 0.4;
   }
   for(i=0; i<6; i++){
     printf("%d:x of float test vector: %f\n",i,blub[i].x);
     printf("%d:y of float test vector: %f\n",i,blub[i].y);
     printf("%d:z of float test vector: %f\n",i,blub[i].z);
     printf("%d:w of float test vector: %f\n",i,blub[i].w);   
   }


   cudaMemcpy(s, &(blub[0]) , 6*sizeof(float4), cudaMemcpyHostToDevice);
   
   //END DEBUG_FLO


   //printf("Converting spinor to half precision... ");
     float2half_spinorfield <<< gridsize, BLOCK2  >>>(s, sh, shnorm);
   //printf("Done\n");
    cudaError_t cudaerr;
    if((cudaerr=cudaGetLastError()) != cudaSuccess){
       printf("%s\n", cudaGetErrorString(cudaerr));
       exit(200);
     }


   //DEBUG_FLO
   short4 testnorm_half[6];
   float thenorm;
   float4 testnorm[6];
   
   cudaMemcpy(&(testnorm_half), sh, 6*sizeof(short4), cudaMemcpyDeviceToHost);
   cudaMemcpy(&(thenorm), shnorm, sizeof(float), cudaMemcpyDeviceToHost);
   
   printf("norm of float test vector: %f\n",thenorm);
   printf("%f\n", sh2fl_host((short)(-32767)));
   printf("%f\n", sh2fl_host((short)(32767)));
   
   printf("%d\n", fl2sh_host(-1.0));
   printf("%d\n", fl2sh_host(1.0));
   for(i=0; i<6; i++){
     testnorm[i].x = half2float_host(testnorm_half[i].x, thenorm);
     testnorm[i].y = half2float_host(testnorm_half[i].y, thenorm);
     testnorm[i].z = half2float_host(testnorm_half[i].z, thenorm);
     testnorm[i].w = half2float_host(testnorm_half[i].w, thenorm);
     printf("%d:x of float test vector: %f\n",i,testnorm[i].x);
     printf("%d:y of float test vector: %f\n",i,testnorm[i].y);
     printf("%d:z of float test vector: %f\n",i,testnorm[i].z);
     printf("%d:w of float test vector: %f\n",i,testnorm[i].w);   
   }
    cudaBindTexture(0, spinhalf_tex,sh, size);
    // bind texture for norm
    cudaBindTexture(0, spinnormhalf_tex, shnorm,  sizenorm);
    
    if((cudaerr=cudaGetLastError()) != cudaSuccess){
       printf("%s\n", cudaGetErrorString(cudaerr));
       exit(200);
     }
   
   half2float_spinorfield_tex <<< gridsize, BLOCK2 >>>(dev_spin4);
    
    if((cudaerr=cudaGetLastError()) != cudaSuccess){
       printf("%s\n", cudaGetErrorString(cudaerr));
       exit(200);
     }
   
   cudaUnbindTexture(spinhalf_tex);
   cudaUnbindTexture(spinnormhalf_tex);
   
   cudaMemcpy(&(testnorm), dev_spin4, 6*sizeof(float4), cudaMemcpyDeviceToHost);
   for(i=0; i<6; i++){
     printf("%d:x of float test vector: %f\n",i,testnorm[i].x);
     printf("%d:y of float test vector: %f\n",i,testnorm[i].y);
     printf("%d:z of float test vector: %f\n",i,testnorm[i].z);
     printf("%d:w of float test vector: %f\n",i,testnorm[i].w);   
   }

   //exit(100);
   //END DEBUG
   */

   //printf("Binding textures to half spinorfield\n");
    // bind texture for vector
    cudaBindTexture(0, spinhalf_tex0, sh, size);
    cudaBindTexture(0, spinhalf_tex1, sh+offset, size);
    cudaBindTexture(0, spinhalf_tex2, sh+2*offset, size);
    cudaBindTexture(0, spinhalf_tex3, sh+3*offset, size);
    cudaBindTexture(0, spinhalf_tex4, sh+4*offset, size);
    cudaBindTexture(0, spinhalf_tex5, sh+5*offset, size);
    
    
    // bind texture for norm
    cudaBindTexture(0, spinnormhalf_tex, shnorm, sizenorm);
  //printf("%s\n", cudaGetErrorString(cudaGetLastError())); 
 

return(0);  
}


extern "C" int unbind_halfspinor_texture(){
   //printf("Unbinding textures of half spinorfield\n");
   cudaUnbindTexture(spinhalf_tex0);
   cudaUnbindTexture(spinhalf_tex1);
   cudaUnbindTexture(spinhalf_tex2);
   cudaUnbindTexture(spinhalf_tex3);
   cudaUnbindTexture(spinhalf_tex4);
   cudaUnbindTexture(spinhalf_tex5);
   cudaUnbindTexture(spinnormhalf_tex);
   //printf("%s\n", cudaGetErrorString(cudaGetLastError()));    
return(0);    
}




extern "C" int bind_texture_gf_half(dev_su3_2v_half * gf){
 //printf("Binding texture to gaugefield\n");
 
  #ifdef MPI
    #ifdef GF_8
     size_t size = sizeof(short4)*2*(VOLUME+RAND)*4;
    #else
     size_t size = sizeof(short4)*3*(VOLUME+RAND)*4;
    #endif
  #else
    #ifdef GF_8
     size_t size = sizeof(short4)*2*VOLUME*4;
    #else
     size_t size = sizeof(short4)*3*VOLUME*4;
    #endif
  #endif
 
 cudaGetTextureReference(&gf_texRefPtr, "gf_tex");
 gf_channelDesc =  cudaCreateChannelDesc<short4>();
 cudaBindTexture(0, gf_texRefPtr, gf, &gf_channelDesc, size);
 printf("%s\n", cudaGetErrorString(cudaGetLastError()));    
 return(0);
}



extern "C" int unbind_texture_gf_half(){
 //printf("Unbinding texture to gaugefield\n");
 cudaUnbindTexture(gf_texRefPtr);
 printf("%s\n", cudaGetErrorString(cudaGetLastError()));    
 return(0);
}



// convert spinor to REAL4 (float4, double4) 
void convert2REAL4_spin_half(spinor* spin, dev_spinor_half* h2d, float* h2d_norm){
  int i,Vol,offset;
  float norm;
  if(even_odd_flag){
    Vol = VOLUME/2;
    offset = VOLUME/2;
  }
  else{
    Vol = VOLUME;
    offset = VOLUME;
  }
  for (i=0;i<Vol;i++){
    get_half_norm_from_pos_host(norm, spin, i);
    h2d_norm[i] = norm;
        h2d[i+0*offset].x = fl2sh_host( (float)spin[i].s0.c0.re/norm);
        h2d[i+0*offset].y = fl2sh_host( (float)spin[i].s0.c0.im/norm);
        h2d[i+0*offset].z = fl2sh_host( (float)spin[i].s0.c1.re/norm);
        h2d[i+0*offset].w = fl2sh_host( (float)spin[i].s0.c1.im/norm);
        
        h2d[i+1*offset].x = fl2sh_host( (float)spin[i].s0.c2.re/norm);
        h2d[i+1*offset].y = fl2sh_host( (float)spin[i].s0.c2.im/norm);
        h2d[i+1*offset].z = fl2sh_host( (float)spin[i].s1.c0.re/norm);
        h2d[i+1*offset].w = fl2sh_host( (float)spin[i].s1.c0.im/norm);
        
        h2d[i+2*offset].x = fl2sh_host( (float)spin[i].s1.c1.re/norm);
        h2d[i+2*offset].y = fl2sh_host( (float)spin[i].s1.c1.im/norm);
        h2d[i+2*offset].z = fl2sh_host( (float)spin[i].s1.c2.re/norm);
        h2d[i+2*offset].w = fl2sh_host( (float)spin[i].s1.c2.im/norm);
        
        h2d[i+3*offset].x = fl2sh_host( (float)spin[i].s2.c0.re/norm);
        h2d[i+3*offset].y = fl2sh_host( (float)spin[i].s2.c0.im/norm);
        h2d[i+3*offset].z = fl2sh_host( (float)spin[i].s2.c1.re/norm);
        h2d[i+3*offset].w = fl2sh_host( (float)spin[i].s2.c1.im/norm);
        
        h2d[i+4*offset].x = fl2sh_host( (float)spin[i].s2.c2.re/norm);
        h2d[i+4*offset].y = fl2sh_host( (float)spin[i].s2.c2.im/norm);
        h2d[i+4*offset].z = fl2sh_host( (float)spin[i].s3.c0.re/norm);
        h2d[i+4*offset].w = fl2sh_host( (float)spin[i].s3.c0.im/norm);
        
        h2d[i+5*offset].x = fl2sh_host( (float)spin[i].s3.c1.re/norm);
        h2d[i+5*offset].y = fl2sh_host( (float)spin[i].s3.c1.im/norm);
        h2d[i+5*offset].z = fl2sh_host( (float)spin[i].s3.c2.re/norm);
        h2d[i+5*offset].w = fl2sh_host( (float)spin[i].s3.c2.im/norm);
    
  }
}




// convert spinor to double 
void convert2double_spin_half(dev_spinor_half* spin, float* spin_norm, spinor* h2d){
  int i,Vol, offset;
  double norm;
  if(even_odd_flag){
    Vol = VOLUME/2;
    offset = VOLUME/2;
  }
  else{
    Vol = VOLUME;
    offset = VOLUME;
  }
  for (i=0;i<Vol;i++){
     norm=(double) spin_norm[i];
        h2d[i].s0.c0.re = (double) (sh2fl_host(spin[i+0*offset].x))*norm;
        h2d[i].s0.c0.im = (double) (sh2fl_host(spin[i+0*offset].y))*norm;
        h2d[i].s0.c1.re = (double) (sh2fl_host(spin[i+0*offset].z))*norm;
        h2d[i].s0.c1.im = (double) (sh2fl_host(spin[i+0*offset].w))*norm;
        
        h2d[i].s0.c2.re = (double) (sh2fl_host(spin[i+1*offset].x))*norm;
        h2d[i].s0.c2.im = (double) (sh2fl_host(spin[i+1*offset].y))*norm;
        h2d[i].s1.c0.re = (double) (sh2fl_host(spin[i+1*offset].z))*norm;
        h2d[i].s1.c0.im = (double) (sh2fl_host(spin[i+1*offset].w))*norm;   
        
        h2d[i].s1.c1.re = (double) (sh2fl_host(spin[i+2*offset].x))*norm;
        h2d[i].s1.c1.im = (double) (sh2fl_host(spin[i+2*offset].y))*norm;
        h2d[i].s1.c2.re = (double) (sh2fl_host(spin[i+2*offset].z))*norm;
        h2d[i].s1.c2.im = (double) (sh2fl_host(spin[i+2*offset].w))*norm;  
        
        h2d[i].s2.c0.re = (double) (sh2fl_host(spin[i+3*offset].x))*norm;
        h2d[i].s2.c0.im = (double) (sh2fl_host(spin[i+3*offset].y))*norm;
        h2d[i].s2.c1.re = (double) (sh2fl_host(spin[i+3*offset].z))*norm;
        h2d[i].s2.c1.im = (double) (sh2fl_host(spin[i+3*offset].w))*norm;  
        
        h2d[i].s2.c2.re = (double) (sh2fl_host(spin[i+4*offset].x))*norm;
        h2d[i].s2.c2.im = (double) (sh2fl_host(spin[i+4*offset].y))*norm;
        h2d[i].s3.c0.re = (double) (sh2fl_host(spin[i+4*offset].z))*norm;
        h2d[i].s3.c0.im = (double) (sh2fl_host(spin[i+4*offset].w))*norm; 
        
        h2d[i].s3.c1.re = (double) (sh2fl_host(spin[i+5*offset].x))*norm;
        h2d[i].s3.c1.im = (double) (sh2fl_host(spin[i+5*offset].y))*norm;
        h2d[i].s3.c2.re = (double) (sh2fl_host(spin[i+5*offset].z))*norm;
        h2d[i].s3.c2.im = (double) (sh2fl_host(spin[i+5*offset].w))*norm; 
        
  }
}










/////    SOME BLAS KERNELs - work in progress //////////



// y(half) = alpha*x(half) + y(half) 
// x is not read from texture
// y is not read from texture
__global__ void axpy_half (float alpha, dev_spinor_half* x, float* x_norm, dev_spinor_half* y, float* y_norm){
   int pos= threadIdx.x + blockDim.x*blockIdx.x;
   float4 xhelp[6]; 
   float4 erghelp[6];
   int i;
   float xnhelp, ynhelp;
   
   
   if(pos < dev_VOLUME){
    xnhelp = x_norm[pos];
    construct_spinor_fromhalf(xhelp, x,  xnhelp, pos);
    ynhelp = y_norm[pos];    
    construct_spinor_fromhalf(erghelp, y,  ynhelp, pos);
    
    for(i=0; i<6; i++){
       erghelp[i].x += alpha*xhelp[i].x;
       erghelp[i].y += alpha*xhelp[i].y;
       erghelp[i].z += alpha*xhelp[i].z;
       erghelp[i].w += alpha*xhelp[i].w;
    }
        
	//calculate norm of resulting spinor
	float ergnorm = 0.0f;


//CHECK THIS
//is this needed or just double work?
        get_half_norm(ergnorm, erghelp);
	y_norm[pos] = ergnorm;
//CHECK THIS
	
	
	//write out normalized spinors in half
    if (ergnorm != 0.0f){
      dev_write_spinor_half(&(erghelp[0]),&(y[pos]), &(y_norm[pos]));
      }
      else{
       dev_write_spinor_half(&(erghelp[0]),&(y[pos]), &(y_norm[pos]));
     }
   }//dev_VOLUME
}




// x = alpha*x(half) 
// x is not read from texture
__global__ void scal_half (float alpha, dev_spinor_half* x, float* x_norm){
   int pos= threadIdx.x + blockDim.x*blockIdx.x;
   float4 xhelp;
   float4 erghelp[6];
   int i;
   float xnhelp;
   
   
   if(pos < dev_VOLUME){
    // this is the loop over the 6 float4 forming one spinor
    #pragma unroll 6
    for(i=0; i<6; i++){
       //xhelp = tex1Dfetch(spinhalf_tex, 6*pos+i);
       //xnhelp = tex1Dfetch(spinnormhalf_tex, pos);
       
       xnhelp = x_norm[pos];
       xhelp.x = sh2fl(x[pos+i*DEVOFF].x)*xnhelp;
       xhelp.y = sh2fl(x[pos+i*DEVOFF].y)*xnhelp;
       xhelp.z = sh2fl(x[pos+i*DEVOFF].z)*xnhelp;
       xhelp.w = sh2fl(x[pos+i*DEVOFF].w)*xnhelp;
        
       erghelp[i].x = alpha*xhelp.x;
       erghelp[i].y = alpha*xhelp.y;
       erghelp[i].z = alpha*xhelp.z;
       erghelp[i].w = alpha*xhelp.w;
    }
    
        //calculate norm of resulting spinor
        float ergnorm = 0.0f;
        
        get_half_norm(ergnorm, erghelp);
        x_norm[pos] = ergnorm;
	float invnorm = 1.0f/ergnorm;

        //write out normalized spinors in half
    if (ergnorm != 0.0f){
      #pragma unroll 6
      for(i=0; i<6; i++){
         x[pos+i*DEVOFF].x = fl2sh(erghelp[i].x*invnorm);
         x[pos+i*DEVOFF].y = fl2sh(erghelp[i].y*invnorm);
         x[pos+i*DEVOFF].z = fl2sh(erghelp[i].z*invnorm);
         x[pos+i*DEVOFF].w = fl2sh(erghelp[i].w*invnorm);
          } 
        }
      else{
       #pragma unroll 6
       for(i=0; i<6; i++){
         x[pos+i*DEVOFF].x = fl2sh(0.0f);
         x[pos+i*DEVOFF].y = fl2sh(0.0f);
         x[pos+i*DEVOFF].z = fl2sh(0.0f);
         x[pos+i*DEVOFF].w = fl2sh(0.0f);
       }
     }


   }//dev_VOLUME
}












void init_blas_half(int vol){
  cudaError_t cudaerr;
  
  blas_half_blocksize=BLOCK2;
  if( vol >= BLOCK2){
   blas_half_gridsize = (int)(vol/BLOCK2) + 1;
  }
  else{
    blas_half_gridsize=1;
  }
  
  size_t size = vol * sizeof(float);
  
  if((cudaerr=cudaMalloc((void **) &dev_blas_half_redfield, size)) != cudaSuccess){
    printf("Error in init_blas_half(): Memory allocation of reduction field failed. Aborting...\n");
    exit(200);
  }   // Allocate array on device
  else{
    printf("Allocated blas reduction field on device\n");
  }  


  // IMPLEMENT THIS FOR ALL LATTICE SIZES !!!!!!!!!!!!!!!!!!!!
  if((vol%REDUCTION_N) == 0){
    blas_half_redblocks = vol/REDUCTION_N;
  }
  else{
    fprintf(stderr,"Error: Volume is not a multiple of REDUCTION_N (%d). Aborting...\n", REDUCTION_N);
    exit(100);
  }
  
  // initialize small redfields
  size = blas_half_redblocks * sizeof(float);
  if((cudaerr=cudaMalloc((void **) &dev_blas_half_sredfield, size)) != cudaSuccess){
    printf("Error in init_blas_half(): Memory allocation of small reduction field failed. Aborting...\n");
    exit(200);
  }   // Allocate array on device
  else{
    printf("Allocated blas small reduction field on device\n");
  }  
  
  if((void*)(blas_half_sredfield = (float *)malloc(size)) == NULL){
    printf("Could not allocate memory for blas small redfield on host. Aborting...\n");
    exit(200);
  } 
  
  
  
}


void finalize_blas_half(){
  cudaFree(dev_blas_half_redfield);
  cudaFree(dev_blas_half_sredfield);
  free(blas_half_sredfield);
}










__global__ void dot_half ( float* redfield, dev_spinor_half* x, float* x_norm, dev_spinor_half* y, float* y_norm){
   int pos= threadIdx.x + blockDim.x*blockIdx.x;
   float4 xhelp,yhelp;
   int i;
   float xnhelp, ynhelp;
   float dotp = 0.0f;
   
   if(pos < dev_VOLUME){
    // this is the loop over the 6 float4 forming one spinor
    #pragma unroll 6
    for(i=0; i<6; i++){
       //xhelp = tex1Dfetch(spinhalf_tex, 6*pos+i);
       //xnhelp = tex1Dfetch(spinnormhalf_tex, pos);
       
       xnhelp = x_norm[pos];
       xhelp.x = sh2fl(x[pos+i*DEVOFF].x)*xnhelp;
       xhelp.y = sh2fl(x[pos+i*DEVOFF].y)*xnhelp;
       xhelp.z = sh2fl(x[pos+i*DEVOFF].z)*xnhelp;
       xhelp.w = sh2fl(x[pos+i*DEVOFF].w)*xnhelp;
      
       ynhelp = y_norm[pos];
       yhelp.x = sh2fl(y[pos+i*DEVOFF].x)*ynhelp;
       yhelp.y = sh2fl(y[pos+i*DEVOFF].y)*ynhelp;
       yhelp.z = sh2fl(y[pos+i*DEVOFF].z)*ynhelp;
       yhelp.w = sh2fl(y[pos+i*DEVOFF].w)*ynhelp;
        
       dotp += xhelp.x * yhelp.x;
       dotp += xhelp.y * yhelp.y;
       dotp += xhelp.z * yhelp.z;
       dotp += xhelp.w * yhelp.w;
    }
    // write sum_i (x_i y_i) to reduction field 
    redfield[pos] = dotp;
   }//dev_VOLUME
}



// kernel for the square of the norm of a half spinor
// local squared norms are written to reduction field redfield
__global__ void sqnorm_half (float* redfield, dev_spinor_half* x, float* x_norm){
   int pos= threadIdx.x + blockDim.x*blockIdx.x;
   float xnhelp;
   float dotp = 0.0;
   int i;
   float4 xhelp;
   if(pos < dev_VOLUME){
    #pragma unroll 6
    for(i=0; i<6; i++){
       xnhelp = x_norm[pos];
       xhelp.x = sh2fl(x[pos+i*DEVOFF].x)*xnhelp;
       xhelp.y = sh2fl(x[pos+i*DEVOFF].y)*xnhelp;
       xhelp.z = sh2fl(x[pos+i*DEVOFF].z)*xnhelp;
       xhelp.w = sh2fl(x[pos+i*DEVOFF].w)*xnhelp;
       
       dotp += xhelp.x * xhelp.x;
       dotp += xhelp.y * xhelp.y;
       dotp += xhelp.z * xhelp.z;
       dotp += xhelp.w * xhelp.w;
    }
    redfield[pos] = dotp;
   }//dev_VOLUME
}







// calculates the dot product of x and y
float dotprod_half(dev_spinor_half* x, float* x_norm, dev_spinor_half* y, float* y_norm){
   int i;
   float result;
   cudaError_t cudaerr;
   
   dot_half <<< blas_half_gridsize, blas_half_blocksize >>> 
                      (dev_blas_half_redfield, x, x_norm, y, y_norm);
   if((cudaerr=cudaGetLastError()) != cudaSuccess){
      printf("%s\n", cudaGetErrorString(cudaerr));
      exit(200);
   } 
   //reduce reductionfield on device 
   reduce_float <<< blas_half_redblocks, REDUCTION_N, 
                REDUCTION_N*sizeof(float) >>> 
                ( dev_blas_half_redfield, dev_blas_half_sredfield,  VOLUME);
   //this reduction always takes the VOLUME (also for mpi)     
   
   //copy back
   cudaMemcpy(blas_half_sredfield, dev_blas_half_sredfield, (size_t)(blas_half_redblocks*sizeof(float)), cudaMemcpyDeviceToHost);
           
   //do final reduction on host
   float finalsum=0.0f;
   for(i=0; i<blas_half_redblocks; i++){
     finalsum += blas_half_sredfield[i];
   }
   #ifdef MPI
     MPI_Allreduce(&finalsum, &result, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
     finalsum=result;
   #endif
   return(finalsum);
}



// calculates the norm^2 of a half spinor with norm x 
float squarenorm_half(dev_spinor_half* x, float * xnorm){
   int i;
   float result;
   cudaError_t cudaerr;
   
   sqnorm_half <<< blas_half_gridsize, blas_half_blocksize >>> 
                      (dev_blas_half_redfield, x, xnorm);
   if((cudaerr=cudaGetLastError()) != cudaSuccess){
      printf("%s\n", cudaGetErrorString(cudaerr));
      exit(200);
   } 
   //reduce reductionfield on device 
   reduce_float <<< blas_half_redblocks, REDUCTION_N, 
                REDUCTION_N*sizeof(float) >>> 
                ( dev_blas_half_redfield, dev_blas_half_sredfield,  VOLUME);
   //this reduction always takes the VOLUME (also for mpi)     
   
   //copy back
   cudaMemcpy(blas_half_sredfield, dev_blas_half_sredfield, (size_t)(blas_half_redblocks*sizeof(float)), cudaMemcpyDeviceToHost);
           
   //do final reduction on host
   float finalsum=0.0f;
   for(i=0; i<blas_half_redblocks; i++){
     finalsum += blas_half_sredfield[i];
   }
   #ifdef MPI
     MPI_Allreduce(&finalsum, &result, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
     finalsum=result;
   #endif
   return(finalsum);
}






/////////////////////  BLAS KERNELs  ////////////////////////////








//////////////////   TEST ROUTINES /////////////////////////////


__global__ void testhalf (dev_su3_2v* gf, dev_su3 * to){

  int pos=threadIdx.x + blockDim.x*blockIdx.x;
  if (pos==0){
   #ifdef GF_8
    dev_reconstructgf_8texref(gf, 2 , to);
   #else
    dev_reconstructgf_2vtexref(gf, 2 , to);
   #endif
  }
}


void testhalf_gf(dev_su3_2v_half * gf){
  dev_su3 * testfield;
  dev_su3 hosttestfield;
  size_t size =  sizeof(dev_su3);
  dev_su3 hostmatrix;

  
  cudaMalloc((void **) &testfield, size);
  #ifdef USETEXTURE
    //Bind texture gf
    bind_texture_gf_half(dev_gf_half);
  #endif  
  testhalf <<< 1, 1 >>> (dev_gf, testfield);
  cudaMemcpy(&(hosttestfield), testfield, size, cudaMemcpyDeviceToHost);  
  show_dev_su3(hosttestfield);
  
  
  #ifdef USETEXTURE
   unbind_texture_gf_half();
  #endif
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



__global__ void to_relativistic_basis_half(
     dev_spinor_half* spinin, float* spinin_norm){
   int pos= threadIdx.x + blockDim.x*blockIdx.x;
   float4 inhelp[6];
   float4 outhelp[6];
   int i;
   float norminhelp, normouthelp;
   const float sq2 = rsqrtf(2.0f);
   
   if(pos < dev_VOLUME){
    norminhelp = spinin_norm[pos];
    construct_spinor_fromhalf(inhelp, spinin,  norminhelp, pos);
    

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

    



     //write out normalized spinors in half
     dev_write_spinor_half(&(outhelp[0]),&(spinin[pos]), &(spinin_norm[pos]));

   }//dev_VOLUME

}




__global__ void to_tmlqcd_basis_half(
     dev_spinor_half* spinin, float* spinin_norm){
  
   int pos= threadIdx.x + blockDim.x*blockIdx.x;
   float4 inhelp[6];
   float4 outhelp[6];
   int i;
   float norminhelp, normouthelp;
   const float sq2 = rsqrtf(2.0f);
   
   if(pos < dev_VOLUME){
    norminhelp = spinin_norm[pos];
    construct_spinor_fromhalf(inhelp, spinin,  norminhelp, pos);
    

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

    



     //write out normalized spinors in half
     dev_write_spinor_half(&(outhelp[0]),&(spinin[pos]), &(spinin_norm[pos]));

   }//dev_VOLUME


}



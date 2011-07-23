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
#define SCALE ((SHORT_LEN-1) * 0.5)
#define SHIFT (-1.0f/(SHORT_LEN-1))

/*
__device__ short fl2sh(float f) {
  short ret = (short)((f+SHIFT)*SCALE);
  return ret;
}

__device__ float sh2fl(short s) {
  return ((float)(s/SCALE) - SHIFT);
}


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

*/


#define fl2sh(f) ((short)(((f)+SHIFT)*SCALE))
#define sh2fl(s) ((float)((s)/SCALE) - SHIFT)
#define half2fl(s,norm) (norm*((float)((s)/SCALE) - SHIFT))

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
   (sf)[0].x = shn*sh2fl(sh[6*(pos)].x);                        \
   (sf)[0].y = shn*sh2fl(sh[6*(pos)].y);                        \
   (sf)[0].z = shn*sh2fl(sh[6*(pos)].z);                        \
   (sf)[0].w = shn*sh2fl(sh[6*(pos)].w);                        \
   (sf)[1].x = shn*sh2fl(sh[6*(pos)+1].x);                        \
   (sf)[1].y = shn*sh2fl(sh[6*(pos)+1].y);                        \
   (sf)[1].z = shn*sh2fl(sh[6*(pos)+1].z);                        \
   (sf)[1].w = shn*sh2fl(sh[6*(pos)+1].w);                        \
   (sf)[2].x = shn*sh2fl(sh[6*(pos)+2].x);                        \
   (sf)[2].y = shn*sh2fl(sh[6*(pos)+2].y);                        \
   (sf)[2].z = shn*sh2fl(sh[6*(pos)+2].z);                        \
   (sf)[2].w = shn*sh2fl(sh[6*(pos)+2].w);                        \
   (sf)[3].x = shn*sh2fl(sh[6*(pos)+3].x);                        \
   (sf)[3].y = shn*sh2fl(sh[6*(pos)+3].y);                        \
   (sf)[3].z = shn*sh2fl(sh[6*(pos)+3].z);                        \
   (sf)[3].w = shn*sh2fl(sh[6*(pos)+3].w);                        \
   (sf)[4].x = shn*sh2fl(sh[6*(pos)+4].x);                        \
   (sf)[4].y = shn*sh2fl(sh[6*(pos)+4].y);                        \
   (sf)[4].z = shn*sh2fl(sh[6*(pos)+4].z);                        \
   (sf)[4].w = shn*sh2fl(sh[6*(pos)+4].w);                        \
   (sf)[5].x = shn*sh2fl(sh[6*(pos)+5].x);                        \
   (sf)[5].y = shn*sh2fl(sh[6*(pos)+5].y);                        \
   (sf)[5].z = shn*sh2fl(sh[6*(pos)+5].z);                        \
   (sf)[5].w = shn*sh2fl(sh[6*(pos)+5].w);                       }\
   
   
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
  float c0 = fmaxf(fabsf((s[6*pos+0]).x), fabsf((s[6*pos+0]).y));                   \
  float c1 = fmaxf(fabsf((s[6*pos+0]).z), fabsf((s[6*pos+0]).w));                   \
  float c2 = fmaxf(fabsf((s[6*pos+1]).x), fabsf((s[6*pos+1]).y));                   \
  float c3 = fmaxf(fabsf((s[6*pos+1]).z), fabsf((s[6*pos+1]).w));                   \
  float c4 = fmaxf(fabsf((s[6*pos+2]).x), fabsf((s[6*pos+2]).y));                   \
  float c5 = fmaxf(fabsf((s[6*pos+2]).z), fabsf((s[6*pos+2]).w));                   \
  float c6 = fmaxf(fabsf((s[6*pos+3]).x), fabsf((s[6*pos+3]).y));                   \
  float c7 = fmaxf(fabsf((s[6*pos+3]).z), fabsf((s[6*pos+3]).w));                   \
  float c8 = fmaxf(fabsf((s[6*pos+4]).x), fabsf((s[6*pos+4]).y));                   \
  float c9 = fmaxf(fabsf((s[6*pos+4]).z), fabsf((s[6*pos+4]).w));                   \
  float c10 = fmaxf(fabsf((s[6*pos+5]).x), fabsf((s[6*pos+5]).y));                  \
  float c11 = fmaxf(fabsf((s[6*pos+5]).z), fabsf((s[6*pos+5]).w));                  \
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
 texture<short4,1, cudaReadModeNormalizedFloat> spinhalf_tex;
 
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
__device__ void dev_write_spinor_half(dev_spinor* in, dev_spinor_half* out, float* out_norm){
   float norm = 0.0f;
   int i;
   
    get_half_norm(norm, in);
   
    //store unit direction vector
    *out_norm = norm;
    if (norm != 0.0f){
      //store norm
      #pragma unroll 6
      for(i=0; i<6; i++){
         out[i].x = fl2sh(in[i].x/norm);
         out[i].y = fl2sh(in[i].y/norm);
         out[i].z = fl2sh(in[i].z/norm);
         out[i].w = fl2sh(in[i].w/norm);
         } 
       }
      else{
       //store norm
       #pragma unroll 6
       for(i=0; i<6; i++){
         out[i].x = fl2sh(0.0f);
         out[i].y = fl2sh(0.0f);
         out[i].z = fl2sh(0.0f);
         out[i].w = fl2sh(0.0f);
       }
     }
}              
              
              
              
              
              
              

// stores the float spinor field in s into the half spinor field sh and the norm into shnorm
__global__ void float2half_spinorfield(dev_spinor* s, dev_spinor_half* sh, float* shnorm){
  int pos=threadIdx.x + blockDim.x*blockIdx.x;
  //__shared__ float4 slocal[6];
  int i;
  float norm = 0.0;
  if(pos < dev_VOLUME){
   
   /*  BEWARE THIS IS NOT WORKING FOR SOME REASON dev_copy_spinor fails because slocal is shared???
   dev_copy_spinor(&(s[6*pos]), &(slocal[0]));
    // calculate norm
   
   
    for(i=0; i<6; i++){
       norm += slocal[i].x*slocal[i].x + slocal[i].y*slocal[i].y 
             + slocal[i].z*slocal[i].z + slocal[i].w*slocal[i].w; 
     }
   
    */ 
     
     get_half_norm_from_pos(norm, s, pos);
     
    shnorm[pos] = norm;
    //store unit direction vector
    if (norm != 0.0f){
      //store norm
      #pragma unroll 6
      for(i=0; i<6; i++){
         sh[6*pos+i].x = fl2sh(s[6*pos+i].x/norm);
         sh[6*pos+i].y = fl2sh(s[6*pos+i].y/norm);
         sh[6*pos+i].z = fl2sh(s[6*pos+i].z/norm);
         sh[6*pos+i].w = fl2sh(s[6*pos+i].w/norm);
          } 
        }
      else{
       //store norm
       #pragma unroll 6
       for(i=0; i<6; i++){
         sh[6*pos+i].x = fl2sh(0.0f);
         sh[6*pos+i].y = fl2sh(0.0f);
         sh[6*pos+i].z = fl2sh(0.0f);
         sh[6*pos+i].w = fl2sh(0.0f);
       }
     }
     
     
  }
}


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
    gfh[nf4*(4*pos+mu)].x = fl2sh(gf[nf4*(4*pos+mu)].x);
    gfh[nf4*(4*pos+mu)].y = fl2sh(gf[nf4*(4*pos+mu)].y);
    gfh[nf4*(4*pos+mu)].z = fl2sh(gf[nf4*(4*pos+mu)].z);
    gfh[nf4*(4*pos+mu)].w = fl2sh(gf[nf4*(4*pos+mu)].w);
    
    gfh[nf4*(4*pos+mu)+1].x = fl2sh(gf[nf4*(4*pos+mu)+1].x/pi_float);
    gfh[nf4*(4*pos+mu)+1].y = fl2sh(gf[nf4*(4*pos+mu)+1].y/pi_float);
    gfh[nf4*(4*pos+mu)+1].z = fl2sh(gf[nf4*(4*pos+mu)+1].z);
    gfh[nf4*(4*pos+mu)+1].w = fl2sh(gf[nf4*(4*pos+mu)+1].w);
  #else
    int i;
    nf4 = 3;
    for(i=0; i<nf4; i++){
       gfh[nf4*(4*pos+mu)+i].x = fl2sh(gf[nf4*(4*pos+mu)+i].x);
       gfh[nf4*(4*pos+mu)+i].y = fl2sh(gf[nf4*(4*pos+mu)+i].y);
       gfh[nf4*(4*pos+mu)+i].z = fl2sh(gf[nf4*(4*pos+mu)+i].z);
       gfh[nf4*(4*pos+mu)+i].w = fl2sh(gf[nf4*(4*pos+mu)+i].w);
    } 
  #endif
  }//mu
  }
}





__device__ inline void dev_copy_spinor_half(dev_spinor_half *i1, float* i1_norm, dev_spinor_half *i2, float* i2_norm){
  int i;
  (*i2_norm) = (*i1_norm);
  #pragma unroll 6
  for(i=0;i<6;i++){ //color + spin
    (*(i2+i)).x = (*(i1+i)).x;
    (*(i2+i)).y = (*(i1+i)).y;
    (*(i2+i)).z = (*(i1+i)).z;
    (*(i2+i)).w = (*(i1+i)).w;
  }
}

__device__ inline void dev_zero_spinor_half(dev_spinor_half *sin, float* sin_norm){
  int i;
  *sin_norm = 0.0f;
  #pragma unroll 6
  for(i=0;i<6;i++){ //color + spin
    (*(sin+i)).x = fl2sh(0.0f);
    (*(sin+i)).y = fl2sh(0.0f);
    (*(sin+i)).z = fl2sh(0.0f);
    (*(sin+i)).w = fl2sh(0.0f);
  }
}





__global__ void dev_zero_spinor_field_half(
        dev_spinor_half* s1, float* s1_norm){
  int pos;
  pos= threadIdx.x + blockDim.x*blockIdx.x;  
  if(pos < dev_VOLUME){
          dev_zero_spinor_half(&(s1[6*pos]), &(s1_norm[pos]));
  }
}




__global__ void dev_copy_spinor_field_half(
    dev_spinor_half* s1, float* s1_norm, 
    dev_spinor_half* s2, float* s2_norm
){
    int pos;
  pos= threadIdx.x + blockDim.x*blockIdx.x;  
  if(pos < dev_VOLUME){
      dev_copy_spinor_half( &(s1[6*pos]), &(s1_norm[pos]),
                       &(s2[6*pos]), &(s2_norm[pos]));
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
  int gridsize;
  
  #ifdef MPI
    if(even_odd_flag){
      size = sizeof(short4)*6*(VOLUME+RAND)/2;
      sizenorm = sizeof(float)*(VOLUME+RAND)/2;
    }
    else{
      size = sizeof(short4)*6*(VOLUME+RAND);
      sizenorm = sizeof(float)*(VOLUME+RAND);
    }
  #else
    if(even_odd_flag){
      size = sizeof(short4)*6*VOLUME/2;
      sizenorm = sizeof(float)*VOLUME/2;
    }
    else{
      size = sizeof(short4)*6*VOLUME;
      sizenorm = sizeof(float)*VOLUME;
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
    cudaBindTexture(0, spinhalf_tex, sh, size);
  
    // bind texture for norm
    cudaBindTexture(0, spinnormhalf_tex, shnorm, sizenorm);
  //printf("%s\n", cudaGetErrorString(cudaGetLastError())); 
 

return(0);  
}


extern "C" int unbind_halfspinor_texture(){
   //printf("Unbinding textures of half spinorfield\n");
   cudaUnbindTexture(spinhalf_tex);
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
  int i,Vol;
  float norm;
  if(even_odd_flag){
    Vol = VOLUME/2;
  }
  else{
    Vol = VOLUME;
  }
  for (i=0;i<Vol;i++){
    get_half_norm_from_pos_host(norm, spin, i);
    h2d_norm[i] = norm;
        h2d[6*i+0].x = fl2sh_host( (float)spin[i].s0.c0.re/norm);
        h2d[6*i+0].y = fl2sh_host( (float)spin[i].s0.c0.im/norm);
        h2d[6*i+0].z = fl2sh_host( (float)spin[i].s0.c1.re/norm);
        h2d[6*i+0].w = fl2sh_host( (float)spin[i].s0.c1.im/norm);
        
        h2d[6*i+1].x = fl2sh_host( (float)spin[i].s0.c2.re/norm);
        h2d[6*i+1].y = fl2sh_host( (float)spin[i].s0.c2.im/norm);
        h2d[6*i+1].z = fl2sh_host( (float)spin[i].s1.c0.re/norm);
        h2d[6*i+1].w = fl2sh_host( (float)spin[i].s1.c0.im/norm);
        
        h2d[6*i+2].x = fl2sh_host( (float)spin[i].s1.c1.re/norm);
        h2d[6*i+2].y = fl2sh_host( (float)spin[i].s1.c1.im/norm);
        h2d[6*i+2].z = fl2sh_host( (float)spin[i].s1.c2.re/norm);
        h2d[6*i+2].w = fl2sh_host( (float)spin[i].s1.c2.im/norm);
        
        h2d[6*i+3].x = fl2sh_host( (float)spin[i].s2.c0.re/norm);
        h2d[6*i+3].y = fl2sh_host( (float)spin[i].s2.c0.im/norm);
        h2d[6*i+3].z = fl2sh_host( (float)spin[i].s2.c1.re/norm);
        h2d[6*i+3].w = fl2sh_host( (float)spin[i].s2.c1.im/norm);
        
        h2d[6*i+4].x = fl2sh_host( (float)spin[i].s2.c2.re/norm);
        h2d[6*i+4].y = fl2sh_host( (float)spin[i].s2.c2.im/norm);
        h2d[6*i+4].z = fl2sh_host( (float)spin[i].s3.c0.re/norm);
        h2d[6*i+4].w = fl2sh_host( (float)spin[i].s3.c0.im/norm);
        
        h2d[6*i+5].x = fl2sh_host( (float)spin[i].s3.c1.re/norm);
        h2d[6*i+5].y = fl2sh_host( (float)spin[i].s3.c1.im/norm);
        h2d[6*i+5].z = fl2sh_host( (float)spin[i].s3.c2.re/norm);
        h2d[6*i+5].w = fl2sh_host( (float)spin[i].s3.c2.im/norm);
    
  }
}




// convert spinor to double 
void convert2double_spin_half(dev_spinor_half* spin, float* spin_norm, spinor* h2d){
  int i,Vol;
  double norm;
  if(even_odd_flag){
    Vol = VOLUME/2;
  }
  else{
    Vol = VOLUME;
  }
  for (i=0;i<Vol;i++){
     norm=(double) spin_norm[i];
        h2d[i].s0.c0.re = (double) (sh2fl_host(spin[6*i+0].x))*norm;
        h2d[i].s0.c0.im = (double) (sh2fl_host(spin[6*i+0].y))*norm;
        h2d[i].s0.c1.re = (double) (sh2fl_host(spin[6*i+0].z))*norm;
        h2d[i].s0.c1.im = (double) (sh2fl_host(spin[6*i+0].w))*norm;
        
        h2d[i].s0.c2.re = (double) (sh2fl_host(spin[6*i+1].x))*norm;
        h2d[i].s0.c2.im = (double) (sh2fl_host(spin[6*i+1].y))*norm;
        h2d[i].s1.c0.re = (double) (sh2fl_host(spin[6*i+1].z))*norm;
        h2d[i].s1.c0.im = (double) (sh2fl_host(spin[6*i+1].w))*norm;   
        
        h2d[i].s1.c1.re = (double) (sh2fl_host(spin[6*i+2].x))*norm;
        h2d[i].s1.c1.im = (double) (sh2fl_host(spin[6*i+2].y))*norm;
        h2d[i].s1.c2.re = (double) (sh2fl_host(spin[6*i+2].z))*norm;
        h2d[i].s1.c2.im = (double) (sh2fl_host(spin[6*i+2].w))*norm;  
        
        h2d[i].s2.c0.re = (double) (sh2fl_host(spin[6*i+3].x))*norm;
        h2d[i].s2.c0.im = (double) (sh2fl_host(spin[6*i+3].y))*norm;
        h2d[i].s2.c1.re = (double) (sh2fl_host(spin[6*i+3].z))*norm;
        h2d[i].s2.c1.im = (double) (sh2fl_host(spin[6*i+3].w))*norm;  
        
        h2d[i].s2.c2.re = (double) (sh2fl_host(spin[6*i+4].x))*norm;
        h2d[i].s2.c2.im = (double) (sh2fl_host(spin[6*i+4].y))*norm;
        h2d[i].s3.c0.re = (double) (sh2fl_host(spin[6*i+4].z))*norm;
        h2d[i].s3.c0.im = (double) (sh2fl_host(spin[6*i+4].w))*norm; 
        
        h2d[i].s3.c1.re = (double) (sh2fl_host(spin[6*i+5].x))*norm;
        h2d[i].s3.c1.im = (double) (sh2fl_host(spin[6*i+5].y))*norm;
        h2d[i].s3.c2.re = (double) (sh2fl_host(spin[6*i+5].z))*norm;
        h2d[i].s3.c2.im = (double) (sh2fl_host(spin[6*i+5].w))*norm; 
        
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

        get_half_norm(ergnorm, erghelp);
	y_norm[pos] = ergnorm;

	//write out normalized spinors in half
    if (ergnorm != 0.0f){
      dev_write_spinor_half(&(erghelp[0]),&(y[6*pos]), &(y_norm[pos]));
      }
      else{
       dev_write_spinor_half(&(erghelp[0]),&(y[6*pos]), &(y_norm[pos]));
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
       xhelp.x = sh2fl(x[6*pos+i].x)*xnhelp;
       xhelp.y = sh2fl(x[6*pos+i].y)*xnhelp;
       xhelp.z = sh2fl(x[6*pos+i].z)*xnhelp;
       xhelp.w = sh2fl(x[6*pos+i].w)*xnhelp;
        
       erghelp[i].x = alpha*xhelp.x;
       erghelp[i].y = alpha*xhelp.y;
       erghelp[i].z = alpha*xhelp.z;
       erghelp[i].w = alpha*xhelp.w;
    }
    
        //calculate norm of resulting spinor
        float ergnorm = 0.0f;
        
        get_half_norm(ergnorm, erghelp);
        x_norm[pos] = ergnorm;

        //write out normalized spinors in half
    if (ergnorm != 0.0f){
      #pragma unroll 6
      for(i=0; i<6; i++){
         x[6*pos+i].x = fl2sh(erghelp[i].x/ergnorm);
         x[6*pos+i].y = fl2sh(erghelp[i].y/ergnorm);
         x[6*pos+i].z = fl2sh(erghelp[i].z/ergnorm);
         x[6*pos+i].w = fl2sh(erghelp[i].w/ergnorm);
          } 
        }
      else{
       #pragma unroll 6
       for(i=0; i<6; i++){
         x[6*pos+i].x = fl2sh(0.0f);
         x[6*pos+i].y = fl2sh(0.0f);
         x[6*pos+i].z = fl2sh(0.0f);
         x[6*pos+i].w = fl2sh(0.0f);
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
       xhelp.x = sh2fl(x[6*pos+i].x)*xnhelp;
       xhelp.y = sh2fl(x[6*pos+i].y)*xnhelp;
       xhelp.z = sh2fl(x[6*pos+i].z)*xnhelp;
       xhelp.w = sh2fl(x[6*pos+i].w)*xnhelp;
      
       ynhelp = y_norm[pos];
       yhelp.x = sh2fl(y[6*pos+i].x)*ynhelp;
       yhelp.y = sh2fl(y[6*pos+i].y)*ynhelp;
       yhelp.z = sh2fl(y[6*pos+i].z)*ynhelp;
       yhelp.w = sh2fl(y[6*pos+i].w)*ynhelp;
        
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
       xhelp.x = sh2fl(x[6*pos+i].x)*xnhelp;
       xhelp.y = sh2fl(x[6*pos+i].y)*xnhelp;
       xhelp.z = sh2fl(x[6*pos+i].z)*xnhelp;
       xhelp.w = sh2fl(x[6*pos+i].w)*xnhelp;
       
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

  typedef REAL RealT;
  int pos=threadIdx.x + blockDim.x*blockIdx.x;
  if (pos==0){
   #ifdef GF_8
    dev_reconstructgf_8texref <RealT>(gf, 2 , to);
   #else
    dev_reconstructgf_2vtexref<RealT>(gf, 2 , to);
   #endif
  }
}


void testhalf_gf(dev_su3_2v_half * gf,MixedsolveParameter<REAL>& mixedsolveParameter){
  dev_su3 * testfield;
  dev_su3 hosttestfield;
  size_t size =  sizeof(dev_su3);
  dev_su3 hostmatrix;

  
  cudaMalloc((void **) &testfield, size);
  #ifdef USETEXTURE
    //Bind texture gf
    bind_texture_gf_half(mixedsolveParameter.dev_gf_half);
  #endif  
  testhalf <<< 1, 1 >>> (mixedsolveParameter.dev_gf, testfield);
  cudaMemcpy(&(hosttestfield), testfield, size, cudaMemcpyDeviceToHost);  
  show_dev_su3(hosttestfield);
  
  
  #ifdef USETEXTURE
   unbind_texture_gf_half();
  #endif
}








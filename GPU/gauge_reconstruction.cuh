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
 * File: gauge_reconstruction.cuh
 *
 * CUDA gauge reconstruction functions
 *
 * 
 *
 **************************************************************************/





////////////////////// DEVICE FUNCTIONS FOR GAUGE RECONSTRUCTION ////////////////////


#ifdef HALF
  #define pi_float 3.141592654f
  #define sh4tofl4(fl) make_float4(sh2fl(fl.x), sh2fl(fl.y), sh2fl(fl.z), sh2fl(fl.w))
#else
  #define sh4tofl4(fl) (fl)
#endif 









// reconstruction of the link fields from two rows of the su3 matrix
// numbers are fetched from texture cache
__device__ void dev_reconstructgf_2vtexref (const dev_su3_2v* field, int pos, int mu, int VOL, dev_su3* gf){
  float4 gfin;
  int COL = 3;
  
  #ifdef USETEXTURE
    gfin = tex1Dfetch(gf_tex,pos+VOL*(0 + COL*mu));
  #else
    gfin = field[pos+VOL*(0 + COL*mu)];
  #endif
  //first row
  (*gf)[0][0].re = gfin.x;
  (*gf)[0][0].im = gfin.y;
  (*gf)[0][1].re = gfin.z;
  (*gf)[0][1].im = gfin.w;
  #ifdef USETEXTURE
    gfin = tex1Dfetch(gf_tex,pos+VOL*(1 + COL*mu));
  #else
    gfin = field[pos+VOL*(1 + COL*mu)];
  #endif
  (*gf)[0][2].re = gfin.x;
  (*gf)[0][2].im = gfin.y;
  //second row
  (*gf)[1][0].re = gfin.z;
  (*gf)[1][0].im = gfin.w;
  #ifdef USETEXTURE
    gfin = tex1Dfetch(gf_tex,pos+VOL*(2 + COL*mu));
  #else
    gfin = field[pos+VOL*(2 + COL*mu)];
  #endif
  (*gf)[1][1].re = gfin.x;
  (*gf)[1][1].im = gfin.y;
  (*gf)[1][2].re = gfin.z;
  (*gf)[1][2].im = gfin.w;
  
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
__device__ void dev_reconstructgf_2vtexref_dagger (const dev_su3_2v* field, int pos, int mu, int VOL, dev_su3* gf){
  //dev_complex help1;
  //dev_complex help2;
  float4 gfin;
  int COL = 3;
  
  //first column (minus in im for complex conj.)
  #ifdef USETEXTURE
    gfin = tex1Dfetch(gf_tex,pos+VOL*(0 + COL*mu));
  #else
    gfin = field[pos+VOL*(0 + COL*mu)];
  #endif
  (*gf)[0][0].re = gfin.x;
  (*gf)[0][0].im = -gfin.y;
  (*gf)[1][0].re = gfin.z;
  (*gf)[1][0].im = -gfin.w;
  #ifdef USETEXTURE
    gfin = tex1Dfetch(gf_tex,pos+VOL*(1 + COL*mu));
  #else
    gfin = field[pos+VOL*(1 + COL*mu)];
  #endif
  (*gf)[2][0].re = gfin.x;
  (*gf)[2][0].im = -gfin.y;
  
  //second  column (minus in im for complex conj.)
  (*gf)[0][1].re = gfin.z;
  (*gf)[0][1].im = -gfin.w;
  #ifdef USETEXTURE
    gfin = tex1Dfetch(gf_tex,pos+VOL*(2 + COL*mu));
  #else
    gfin = field[pos+VOL*(2 + COL*mu)];
  #endif
  (*gf)[1][1].re = gfin.x;
  (*gf)[1][1].im = -gfin.y;
  (*gf)[2][1].re = gfin.z;
  (*gf)[2][1].im = -gfin.w;
  


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


#define NEW
#ifdef NEW



/*
// reconstruction of the gf using 8 real parameters as 
// described in the appendix of hep-lat 0911.3191 (M.Clark et al.)
// optimized twice
__device__ void dev_reconstructgf_8texref (const dev_su3_2v * field, int pos, int mu, int VOL, dev_su3* gf){

  float4 gfin;
  float help;
  dev_complex p1,p2;
  int COL = 2;
  #ifdef USETEXTURE
    gfin = tex1Dfetch(gf_tex,pos+VOL*(0 + COL*mu));
  #else
    gfin = field[pos+VOL*(0 + COL*mu)];
  #endif
  // read a2 a3
  (*gf)[0][1].re = gfin.x;
  (*gf)[0][1].im = gfin.y;
  (*gf)[0][2].re = gfin.z;
  (*gf)[0][2].im = gfin.w;  
 
  float a2a3sum = gfin.x*gfin.x + gfin.y*gfin.y + gfin.z*gfin.z + gfin.w*gfin.w; // use later on
  

  // read theta_a1, theta_c1, b1
  #ifdef USETEXTURE
    gfin = tex1Dfetch(gf_tex,pos+VOL*(1 + COL*mu));
  #else
    gfin = field[pos+VOL*(1 + COL*mu)];
  #endif
  
  // reconstruct a1 use x*rsqrt(x) for sqrt(x)
  help = 1.0f - a2a3sum;
  help *= rsqrtf(help>0?help:1.0e14);
  #ifdef HALF
    // we have to multiply by two pi because normalization to -1..1
    gfin.x *= pi_float;
    gfin.y *= pi_float;
  #endif
  __sincosf(gfin.x, &(*gf)[0][0].im, &(*gf)[0][0].re);
  (*gf)[0][0].re *= help;
  (*gf)[0][0].im *= help;
  
   // add a1^2 to sum of first row and assign 1/N^2
   float one_over_N2 = 1.0f/a2a3sum; //__fdividef(1.0f, a2a3sum);
  
  // assign b1
  (*gf)[1][0].re = gfin.z;
  (*gf)[1][0].im = gfin.w;
  

  // reconstruct c1 use x*rsqrt(x) for sqrt(x)
  help =1.0f - 
              (*gf)[0][0].re * (*gf)[0][0].re - (*gf)[0][0].im * (*gf)[0][0].im - 
              (*gf)[1][0].re * (*gf)[1][0].re - (*gf)[1][0].im * (*gf)[1][0].im;

   help *= rsqrtf(help>0?help:1.0e14);
  __sincosf(gfin.y, &(*gf)[2][0].im, &(*gf)[2][0].re);
  (*gf)[2][0].re *= help; 
  (*gf)[2][0].im *= help;
  
  
   
   // p2 = 1/N^2 b1
  p2.re = one_over_N2*(*gf)[1][0].re;
  p2.im = one_over_N2*(*gf)[1][0].im;  

 
  // p1 = 1/N^2*cconj(c1)
  p1.re = one_over_N2*(*gf)[2][0].re;
  p1.im = - one_over_N2*(*gf)[2][0].im;
  
  
  // calculate b2
  (*gf)[1][1] = dev_mcmult_conj(p1,(*gf)[0][2]); 
  dev_cmult_conj_subfrom((*gf)[1][1], p2, (*gf)[0][1]);
 
  // calculate b3
  (*gf)[1][2] = dev_cmult_conj(p1,(*gf)[0][1]);
  dev_complex chelp=dev_cmult_conj((*gf)[0][2], (*gf)[0][0]);
  dev_cmult_subfrom((*gf)[1][2], p2, chelp);  
  
  // calculate c2
  (*gf)[2][1] = dev_cmult(p2,(*gf)[0][2]); 
  chelp=dev_cmult_conj( (*gf)[0][0], (*gf)[0][1]);
  dev_cmult_subfrom((*gf)[2][1], p1, chelp);  
  (*gf)[2][1].im = -(*gf)[2][1].im;
  
  // calculate c3
  (*gf)[2][2] = dev_mcmult(p2,(*gf)[0][1]);
  chelp=dev_cmult_conj( (*gf)[0][0], (*gf)[0][2]);     
  dev_cmult_subfrom((*gf)[2][2], p1, chelp);
  (*gf)[2][2].im = -(*gf)[2][2].im;   
                     
}
*/



#ifdef USETEXTURE
#define dev_reconstructgf_8texref(field, pos, mu, VOL, gf)\
{  float4 gfin;\
  REAL one_over_N, help;\
  dev_complex p1,p2;\
  int COL = 2;\
  gfin = tex1Dfetch(gf_tex,pos+VOL*(0 + COL*mu));\
  (*gf)[0][1].re = gfin.x;\
  (*gf)[0][1].im = gfin.y;\
  (*gf)[0][2].re = gfin.z;\
  (*gf)[0][2].im = gfin.w;  \
  p1.re = gfin.x*gfin.x + gfin.y*gfin.y + gfin.z*gfin.z + gfin.w*gfin.w; \
  one_over_N = rsqrtf(p1.re);\
  gfin = tex1Dfetch(gf_tex,pos+VOL*(1 + COL*mu));\
  help = 1.0f - p1.re;\
  p1.re = help*rsqrtf(help>0?help:1.0e14);  \
  __sincosf(gfin.x, &(*gf)[0][0].im, &(*gf)[0][0].re);\
  (*gf)[0][0].re *=  p1.re;\
  (*gf)[0][0].im *=  p1.re;\
  (*gf)[1][0].re = gfin.z;\
  (*gf)[1][0].im = gfin.w;\
  p2.re = one_over_N*(*gf)[1][0].re;\
  p2.im = one_over_N*(*gf)[1][0].im;  \
  help =1.0f - \
              (*gf)[0][0].re * (*gf)[0][0].re - (*gf)[0][0].im * (*gf)[0][0].im - \
              (*gf)[1][0].re * (*gf)[1][0].re - (*gf)[1][0].im * (*gf)[1][0].im;\
   p1.re = help*rsqrtf(help>0?help:1.0e14);\
  __sincosf(gfin.y, &(*gf)[2][0].im, &(*gf)[2][0].re);\
  (*gf)[2][0].re = (*gf)[2][0].re * p1.re; \
  (*gf)[2][0].im = (*gf)[2][0].im * p1.re;\
  p1.re = one_over_N*(*gf)[2][0].re;\
  p1.im = - one_over_N*(*gf)[2][0].im;\
  (*gf)[1][1].re = p1.re*(*gf)[0][2].re + p1.im*(*gf)[0][2].im;\
  (*gf)[1][1].im = p1.im*(*gf)[0][2].re - p1.re*(*gf)[0][2].im;\
  (*gf)[2][2].re = (*gf)[0][0].re * (*gf)[0][1].re + (*gf)[0][0].im * (*gf)[0][1].im;\
  (*gf)[2][2].im = (*gf)[0][0].re * (*gf)[0][1].im - (*gf)[0][0].im * (*gf)[0][1].re;\
  (*gf)[2][2] = dev_cmult(p2, (*gf)[2][2]);  \
  (*gf)[1][1].re = -one_over_N*( (*gf)[1][1].re + (*gf)[2][2].re);\
  (*gf)[1][1].im = -one_over_N*((*gf)[1][1].im + (*gf)[2][2].im);\
  (*gf)[1][2].re = p1.re*(*gf)[0][1].re + p1.im*(*gf)[0][1].im;\
  (*gf)[1][2].im = p1.im*(*gf)[0][1].re - p1.re*(*gf)[0][1].im;\
  (*gf)[2][2].re = (*gf)[0][0].re*(*gf)[0][2].re + (*gf)[0][0].im*(*gf)[0][2].im;\
  (*gf)[2][2].im = (*gf)[0][0].re*(*gf)[0][2].im - (*gf)[0][0].im*(*gf)[0][2].re;\
  (*gf)[2][2] = dev_cmult(p2,(*gf)[2][2]);\
  (*gf)[1][2].re = one_over_N*( (*gf)[1][2].re - (*gf)[2][2].re);\
  (*gf)[1][2].im = one_over_N*( (*gf)[1][2].im - (*gf)[2][2].im);\
  (*gf)[2][1].re = p2.re*(*gf)[0][2].re - p2.im*(*gf)[0][2].im;\
  (*gf)[2][1].im = -p2.re*(*gf)[0][2].im - p2.im*(*gf)[0][2].re;\
  (*gf)[2][2].re = (*gf)[0][0].re*(*gf)[0][1].re + (*gf)[0][0].im*(*gf)[0][1].im;\
  (*gf)[2][2].im = (*gf)[0][0].re* (*gf)[0][1].im - (*gf)[0][0].im* (*gf)[0][1].re;\
  help = (*gf)[2][2].re;\
  (*gf)[2][2].re = p1.re*(*gf)[2][2].re + p1.im*(*gf)[2][2].im;\
  (*gf)[2][2].im = p1.re*(*gf)[2][2].im - p1.im*help;\
  (*gf)[2][1].re = one_over_N*((*gf)[2][1].re - (*gf)[2][2].re);\
  (*gf)[2][1].im = one_over_N*((*gf)[2][1].im - (*gf)[2][2].im);\
  (*gf)[2][2].re = p2.re * (*gf)[0][1].re - p2.im * (*gf)[0][1].im;\
  (*gf)[2][2].im = - p2.im*(*gf)[0][1].re - p2.re*(*gf)[0][1].im;\
  p2.re = (*gf)[0][0].re * (*gf)[0][2].re + (*gf)[0][0].im * (*gf)[0][2].im;\
  p2.im = (*gf)[0][0].re * (*gf)[0][2].im - (*gf)[0][0].im * (*gf)[0][2].re;\
  p2 = dev_cmult(  dev_cconj(p1) , p2);\
  (*gf)[2][2] = dev_cadd((*gf)[2][2], p2);\
  (*gf)[2][2] = dev_crealmult((*gf)[2][2], -one_over_N);}



#define dev_reconstructgf_8texref_dagger(field, pos, mu, VOL, gf)\
{  float4 gfin;\
  REAL one_over_N, help;\
  dev_complex p1,p2;\
  int COL = 2;\
  gfin = tex1Dfetch(gf_tex,pos+VOL*(0 + COL*mu));\
  (*gf)[1][0].re = gfin.x;\
  (*gf)[1][0].im = -gfin.y;\
  (*gf)[2][0].re = gfin.z;\
  (*gf)[2][0].im = -gfin.w;  \
  p1.re = gfin.x*gfin.x + gfin.y*gfin.y + gfin.z*gfin.z + gfin.w*gfin.w; \
  one_over_N = rsqrtf(p1.re); \
  gfin = tex1Dfetch(gf_tex,pos+VOL*(1 + COL*mu));\
  help = 1.0f - p1.re;\
  p1.re = help*rsqrtf(help>0 ?help:1.0e14 );  \
  __sincosf(gfin.x, &(*gf)[0][0].im, &(*gf)[0][0].re);\
  (*gf)[0][0].re *= p1.re;\
  (*gf)[0][0].im *= -p1.re;\
  (*gf)[0][1].re = gfin.z;\
  (*gf)[0][1].im = -gfin.w;\
  p2.re = one_over_N*(*gf)[0][1].re;\
  p2.im = -one_over_N*(*gf)[0][1].im;  \
  help = 1.0f - \
              (*gf)[0][0].re * (*gf)[0][0].re - (*gf)[0][0].im * (*gf)[0][0].im - \
              (*gf)[0][1].re * (*gf)[0][1].re - (*gf)[0][1].im * (*gf)[0][1].im;\
  p1.re = help*rsqrtf(help>0?help:1.0e14);\
  __sincosf(gfin.y, &(*gf)[0][2].im, &(*gf)[0][2].re);\
  (*gf)[0][2].re = (*gf)[0][2].re * p1.re;\
  (*gf)[0][2].im = -(*gf)[0][2].im * p1.re;\
  p1.re = one_over_N*(*gf)[0][2].re;\
  p1.im = one_over_N*(*gf)[0][2].im;\
  (*gf)[1][1].re  = p1.re*(*gf)[2][0].re - p1.im*(*gf)[2][0].im;\
  (*gf)[1][1].im  = p1.re*(*gf)[2][0].im + p1.im*(*gf)[2][0].re;\
  dev_complex ch;\
  ch.re  = (*gf)[0][0].re*(*gf)[1][0].re + (*gf)[0][0].im*(*gf)[1][0].im;\
  ch.im  = -(*gf)[0][0].re*(*gf)[1][0].im + (*gf)[0][0].im*(*gf)[1][0].re;\
  (*gf)[2][2].re  = p2.re*ch.re - p2.im*ch.im;\
  (*gf)[2][2].im  = p2.re*ch.im + p2.im*ch.re;\
  (*gf)[1][1] = dev_cadd((*gf)[1][1], (*gf)[2][2]);\
  (*gf)[1][1].re = -one_over_N * (*gf)[1][1].re;\
  (*gf)[1][1].im = one_over_N * (*gf)[1][1].im;\
  (*gf)[2][1] = dev_cmult(p1,   (*gf)[1][0]    );\
  ch.re = (*gf)[0][0].re*(*gf)[2][0].re + (*gf)[0][0].im*(*gf)[2][0].im;\
  ch.im = -(*gf)[0][0].re*(*gf)[2][0].im + (*gf)[0][0].im*(*gf)[2][0].re;\
  (*gf)[2][2].re = p2.re*ch.re - p2.im*ch.im;\
  (*gf)[2][2].im = p2.re*ch.im + p2.im*ch.re;\
  (*gf)[2][1] = dev_csub((*gf)[2][1], (*gf)[2][2]);\
  (*gf)[2][1] = dev_cconj(dev_crealmult((*gf)[2][1], one_over_N));\
  (*gf)[1][2] = dev_cmult(  dev_cconj(p2) ,  (*gf)[2][0]    );\
  ch.re = (*gf)[0][0].re*(*gf)[1][0].re + (*gf)[0][0].im*(*gf)[1][0].im;\
  ch.im = -(*gf)[0][0].re*(*gf)[1][0].im + (*gf)[0][0].im*(*gf)[1][0].re;\
  (*gf)[2][2].re = p1.re*ch.re + p1.im*ch.im;\
  (*gf)[2][2].im = p1.re*ch.im - p1.im*ch.re;\
  (*gf)[1][2] = dev_csub((*gf)[1][2], (*gf)[2][2]);\
  (*gf)[1][2] = dev_cconj(dev_crealmult((*gf)[1][2], one_over_N));\
  (*gf)[2][2] = dev_cmult(  dev_cconj(p2) ,   (*gf)[1][0]    );\
  ch.re = (*gf)[0][0].re*(*gf)[2][0].re + (*gf)[0][0].im*(*gf)[2][0].im;\
  ch.im = (*gf)[0][0].im*(*gf)[2][0].re - (*gf)[0][0].re*(*gf)[2][0].im;\
  p2.re = p1.re*ch.re + p1.im*ch.im;\
  p2.im = p1.re*ch.im - p1.im*ch.re;\
  (*gf)[2][2] = dev_cadd((*gf)[2][2], p2);\
  (*gf)[2][2] = dev_cconj(dev_crealmult((*gf)[2][2], -one_over_N));}



#define dev_reconstructgf_2vtexref(field, pos, mu, VOL, gf){ \
  float4 gfin; \
  int COL = 3; \
  gfin = tex1Dfetch(gf_tex,pos+VOL*(0 + COL*mu)); \
  (*gf)[0][0].re = gfin.x; \
  (*gf)[0][0].im = gfin.y; \
  (*gf)[0][1].re = gfin.z; \
  (*gf)[0][1].im = gfin.w; \
  gfin = tex1Dfetch(gf_tex,pos+VOL*(1 + COL*mu)); \
  (*gf)[0][2].re = gfin.x; \
  (*gf)[0][2].im = gfin.y; \
  (*gf)[1][0].re = gfin.z; \
  (*gf)[1][0].im = gfin.w; \
  gfin = tex1Dfetch(gf_tex,pos+VOL*(2 + COL*mu)); \
  (*gf)[1][1].re = gfin.x; \
  (*gf)[1][1].im = gfin.y; \
  (*gf)[1][2].re = gfin.z; \
  (*gf)[1][2].im = gfin.w; \
  (*gf)[2][0].re = (*gf)[0][1].re * (*gf)[1][2].re; \
  (*gf)[2][0].re -= (*gf)[0][1].im * (*gf)[1][2].im; \
  (*gf)[2][0].re -= (*gf)[0][2].re * (*gf)[1][1].re; \
  (*gf)[2][0].re += (*gf)[0][2].im * (*gf)[1][1].im; \
  (*gf)[2][0].im = -(*gf)[0][1].re * (*gf)[1][2].im; \
  (*gf)[2][0].im -= (*gf)[0][1].im * (*gf)[1][2].re; \
  (*gf)[2][0].im += (*gf)[0][2].re * (*gf)[1][1].im; \
  (*gf)[2][0].im += (*gf)[0][2].im * (*gf)[1][1].re; \
  (*gf)[2][1].re = (*gf)[0][2].re * (*gf)[1][0].re; \
  (*gf)[2][1].re -= (*gf)[0][2].im * (*gf)[1][0].im; \
  (*gf)[2][1].re -= (*gf)[0][0].re * (*gf)[1][2].re; \
  (*gf)[2][1].re += (*gf)[0][0].im * (*gf)[1][2].im; \
  (*gf)[2][1].im = -(*gf)[0][2].re * (*gf)[1][0].im; \
  (*gf)[2][1].im -= (*gf)[0][2].im * (*gf)[1][0].re; \
  (*gf)[2][1].im += (*gf)[0][0].re * (*gf)[1][2].im; \
  (*gf)[2][1].im += (*gf)[0][0].im * (*gf)[1][2].re; \
  (*gf)[2][2].re = (*gf)[0][0].re * (*gf)[1][1].re; \
  (*gf)[2][2].re -= (*gf)[0][0].im * (*gf)[1][1].im; \
  (*gf)[2][2].re -= (*gf)[0][1].re * (*gf)[1][0].re; \
  (*gf)[2][2].re += (*gf)[0][1].im * (*gf)[1][0].im; \
  (*gf)[2][2].im = -(*gf)[0][0].re * (*gf)[1][1].im; \
  (*gf)[2][2].im -= (*gf)[0][0].im * (*gf)[1][1].re; \
  (*gf)[2][2].im += (*gf)[0][1].re * (*gf)[1][0].im; \
  (*gf)[2][2].im += (*gf)[0][1].im * (*gf)[1][0].re;}




#define dev_reconstructgf_2vtexref_dagger(field, pos, mu, VOL, gf){ \
  float4 gfin; \
  int COL = 3; \
  gfin = tex1Dfetch(gf_tex,pos+VOL*(0 + COL*mu)); \
  (*gf)[0][0].re = gfin.x; \
  (*gf)[0][0].im = -gfin.y; \
  (*gf)[1][0].re = gfin.z; \
  (*gf)[1][0].im = -gfin.w; \
  gfin = tex1Dfetch(gf_tex,pos+VOL*(1 + COL*mu)); \
  (*gf)[2][0].re = gfin.x; \
  (*gf)[2][0].im = -gfin.y; \
   \
  (*gf)[0][1].re = gfin.z; \
  (*gf)[0][1].im = -gfin.w; \
  gfin = tex1Dfetch(gf_tex,pos+VOL*(2 + COL*mu)); \
  (*gf)[1][1].re = gfin.x; \
  (*gf)[1][1].im = -gfin.y; \
  (*gf)[2][1].re = gfin.z; \
  (*gf)[2][1].im = -gfin.w; \
  (*gf)[0][2].re = (*gf)[1][0].re * (*gf)[2][1].re; \
  (*gf)[0][2].re -= (*gf)[1][0].im * (*gf)[2][1].im; \
  (*gf)[0][2].re -= (*gf)[2][0].re * (*gf)[1][1].re; \
  (*gf)[0][2].re += (*gf)[2][0].im * (*gf)[1][1].im; \
  (*gf)[0][2].im = -(*gf)[1][0].re* (*gf)[2][1].im; \
  (*gf)[0][2].im -= (*gf)[1][0].im* (*gf)[2][1].re; \
  (*gf)[0][2].im += (*gf)[2][0].re*(*gf)[1][1].im; \
  (*gf)[0][2].im += (*gf)[2][0].im*(*gf)[1][1].re; \
  (*gf)[1][2].re = (*gf)[2][0].re*(*gf)[0][1].re; \
  (*gf)[1][2].re -= (*gf)[2][0].im*(*gf)[0][1].im; \
  (*gf)[1][2].re -= (*gf)[0][0].re*(*gf)[2][1].re; \
  (*gf)[1][2].re += (*gf)[0][0].im*(*gf)[2][1].im; \
  (*gf)[1][2].im = -(*gf)[2][0].re * (*gf)[0][1].im; \
  (*gf)[1][2].im -= (*gf)[2][0].im * (*gf)[0][1].re; \
  (*gf)[1][2].im += (*gf)[0][0].re * (*gf)[2][1].im; \
  (*gf)[1][2].im += (*gf)[0][0].im * (*gf)[2][1].re; \
  (*gf)[2][2].re = (*gf)[0][0].re * (*gf)[1][1].re; \
  (*gf)[2][2].re -= (*gf)[0][0].im * (*gf)[1][1].im; \
  (*gf)[2][2].re -= (*gf)[1][0].re * (*gf)[0][1].re; \
  (*gf)[2][2].re += (*gf)[1][0].im * (*gf)[0][1].im; \
  (*gf)[2][2].im = -(*gf)[0][0].re * (*gf)[1][1].im; \
  (*gf)[2][2].im -= (*gf)[0][0].im * (*gf)[1][1].re; \
  (*gf)[2][2].im += (*gf)[1][0].re * (*gf)[0][1].im; \
  (*gf)[2][2].im += (*gf)[1][0].im * (*gf)[0][1].re;} 


#else
#define dev_reconstructgf_8texref(field, pos, mu, VOL, gf)\
{  float4 gfin;\
  REAL one_over_N, help;\
  dev_complex p1,p2;\
  int COL = 2;\
  gfin = field[pos+VOL*(0 + COL*mu)];\
  (*gf)[0][1].re = gfin.x;\
  (*gf)[0][1].im = gfin.y;\
  (*gf)[0][2].re = gfin.z;\
  (*gf)[0][2].im = gfin.w;  \
  p1.re = gfin.x*gfin.x + gfin.y*gfin.y + gfin.z*gfin.z + gfin.w*gfin.w; \
  one_over_N = rsqrtf(p1.re);\
  gfin = field[pos+VOL*(1 + COL*mu)];\
  help = 1.0f - p1.re;\
  p1.re = help*rsqrtf(help>0?help:1.0e14);  \
  __sincosf(gfin.x, &(*gf)[0][0].im, &(*gf)[0][0].re);\
  (*gf)[0][0].re *=  p1.re;\
  (*gf)[0][0].im *=  p1.re;\
  (*gf)[1][0].re = gfin.z;\
  (*gf)[1][0].im = gfin.w;\
  p2.re = one_over_N*(*gf)[1][0].re;\
  p2.im = one_over_N*(*gf)[1][0].im;  \
  help =1.0f - \
              (*gf)[0][0].re * (*gf)[0][0].re - (*gf)[0][0].im * (*gf)[0][0].im - \
              (*gf)[1][0].re * (*gf)[1][0].re - (*gf)[1][0].im * (*gf)[1][0].im;\
   p1.re = help*rsqrtf(help>0?help:1.0e14);\
  __sincosf(gfin.y, &(*gf)[2][0].im, &(*gf)[2][0].re);\
  (*gf)[2][0].re = (*gf)[2][0].re * p1.re; \
  (*gf)[2][0].im = (*gf)[2][0].im * p1.re;\
  p1.re = one_over_N*(*gf)[2][0].re;\
  p1.im = - one_over_N*(*gf)[2][0].im;\
  (*gf)[1][1].re = p1.re*(*gf)[0][2].re + p1.im*(*gf)[0][2].im;\
  (*gf)[1][1].im = p1.im*(*gf)[0][2].re - p1.re*(*gf)[0][2].im;\
  (*gf)[2][2].re = (*gf)[0][0].re * (*gf)[0][1].re + (*gf)[0][0].im * (*gf)[0][1].im;\
  (*gf)[2][2].im = (*gf)[0][0].re * (*gf)[0][1].im - (*gf)[0][0].im * (*gf)[0][1].re;\
  (*gf)[2][2] = dev_cmult(p2, (*gf)[2][2]);  \
  (*gf)[1][1].re = -one_over_N*( (*gf)[1][1].re + (*gf)[2][2].re);\
  (*gf)[1][1].im = -one_over_N*((*gf)[1][1].im + (*gf)[2][2].im);\
  (*gf)[1][2].re = p1.re*(*gf)[0][1].re + p1.im*(*gf)[0][1].im;\
  (*gf)[1][2].im = p1.im*(*gf)[0][1].re - p1.re*(*gf)[0][1].im;\
  (*gf)[2][2].re = (*gf)[0][0].re*(*gf)[0][2].re + (*gf)[0][0].im*(*gf)[0][2].im;\
  (*gf)[2][2].im = (*gf)[0][0].re*(*gf)[0][2].im - (*gf)[0][0].im*(*gf)[0][2].re;\
  (*gf)[2][2] = dev_cmult(p2,(*gf)[2][2]);\
  (*gf)[1][2].re = one_over_N*( (*gf)[1][2].re - (*gf)[2][2].re);\
  (*gf)[1][2].im = one_over_N*( (*gf)[1][2].im - (*gf)[2][2].im);\
  (*gf)[2][1].re = p2.re*(*gf)[0][2].re - p2.im*(*gf)[0][2].im;\
  (*gf)[2][1].im = -p2.re*(*gf)[0][2].im - p2.im*(*gf)[0][2].re;\
  (*gf)[2][2].re = (*gf)[0][0].re*(*gf)[0][1].re + (*gf)[0][0].im*(*gf)[0][1].im;\
  (*gf)[2][2].im = (*gf)[0][0].re* (*gf)[0][1].im - (*gf)[0][0].im* (*gf)[0][1].re;\
  help = (*gf)[2][2].re;\
  (*gf)[2][2].re = p1.re*(*gf)[2][2].re + p1.im*(*gf)[2][2].im;\
  (*gf)[2][2].im = p1.re*(*gf)[2][2].im - p1.im*help;\
  (*gf)[2][1].re = one_over_N*((*gf)[2][1].re - (*gf)[2][2].re);\
  (*gf)[2][1].im = one_over_N*((*gf)[2][1].im - (*gf)[2][2].im);\
  (*gf)[2][2].re = p2.re * (*gf)[0][1].re - p2.im * (*gf)[0][1].im;\
  (*gf)[2][2].im = - p2.im*(*gf)[0][1].re - p2.re*(*gf)[0][1].im;\
  p2.re = (*gf)[0][0].re * (*gf)[0][2].re + (*gf)[0][0].im * (*gf)[0][2].im;\
  p2.im = (*gf)[0][0].re * (*gf)[0][2].im - (*gf)[0][0].im * (*gf)[0][2].re;\
  p2 = dev_cmult(  dev_cconj(p1) , p2);\
  (*gf)[2][2] = dev_cadd((*gf)[2][2], p2);\
  (*gf)[2][2] = dev_crealmult((*gf)[2][2], -one_over_N);}



#define dev_reconstructgf_8texref_dagger(field, pos, mu, VOL, gf)\
{  float4 gfin;\
  REAL one_over_N, help;\
  dev_complex p1,p2;\
  int COL = 2;\
  gfin = field[pos+VOL*(0 + COL*mu)];\
  (*gf)[1][0].re = gfin.x;\
  (*gf)[1][0].im = -gfin.y;\
  (*gf)[2][0].re = gfin.z;\
  (*gf)[2][0].im = -gfin.w;  \
  p1.re = gfin.x*gfin.x + gfin.y*gfin.y + gfin.z*gfin.z + gfin.w*gfin.w; \
  one_over_N = rsqrtf(p1.re); \
  gfin = field[pos+VOL*(1 + COL*mu)];\
  help = 1.0f - p1.re;\
  p1.re = help*rsqrtf(help>0 ?help:1.0e14 );  \
  __sincosf(gfin.x, &(*gf)[0][0].im, &(*gf)[0][0].re);\
  (*gf)[0][0].re *= p1.re;\
  (*gf)[0][0].im *= -p1.re;\
  (*gf)[0][1].re = gfin.z;\
  (*gf)[0][1].im = -gfin.w;\
  p2.re = one_over_N*(*gf)[0][1].re;\
  p2.im = -one_over_N*(*gf)[0][1].im;  \
  help = 1.0f - \
              (*gf)[0][0].re * (*gf)[0][0].re - (*gf)[0][0].im * (*gf)[0][0].im - \
              (*gf)[0][1].re * (*gf)[0][1].re - (*gf)[0][1].im * (*gf)[0][1].im;\
  p1.re = help*rsqrtf(help>0?help:1.0e14);\
  __sincosf(gfin.y, &(*gf)[0][2].im, &(*gf)[0][2].re);\
  (*gf)[0][2].re = (*gf)[0][2].re * p1.re;\
  (*gf)[0][2].im = -(*gf)[0][2].im * p1.re;\
  p1.re = one_over_N*(*gf)[0][2].re;\
  p1.im = one_over_N*(*gf)[0][2].im;\
  (*gf)[1][1].re  = p1.re*(*gf)[2][0].re - p1.im*(*gf)[2][0].im;\
  (*gf)[1][1].im  = p1.re*(*gf)[2][0].im + p1.im*(*gf)[2][0].re;\
  dev_complex ch;\
  ch.re  = (*gf)[0][0].re*(*gf)[1][0].re + (*gf)[0][0].im*(*gf)[1][0].im;\
  ch.im  = -(*gf)[0][0].re*(*gf)[1][0].im + (*gf)[0][0].im*(*gf)[1][0].re;\
  (*gf)[2][2].re  = p2.re*ch.re - p2.im*ch.im;\
  (*gf)[2][2].im  = p2.re*ch.im + p2.im*ch.re;\
  (*gf)[1][1] = dev_cadd((*gf)[1][1], (*gf)[2][2]);\
  (*gf)[1][1].re = -one_over_N * (*gf)[1][1].re;\
  (*gf)[1][1].im = one_over_N * (*gf)[1][1].im;\
  (*gf)[2][1] = dev_cmult(p1,   (*gf)[1][0]    );\
  ch.re = (*gf)[0][0].re*(*gf)[2][0].re + (*gf)[0][0].im*(*gf)[2][0].im;\
  ch.im = -(*gf)[0][0].re*(*gf)[2][0].im + (*gf)[0][0].im*(*gf)[2][0].re;\
  (*gf)[2][2].re = p2.re*ch.re - p2.im*ch.im;\
  (*gf)[2][2].im = p2.re*ch.im + p2.im*ch.re;\
  (*gf)[2][1] = dev_csub((*gf)[2][1], (*gf)[2][2]);\
  (*gf)[2][1] = dev_cconj(dev_crealmult((*gf)[2][1], one_over_N));\
  (*gf)[1][2] = dev_cmult(  dev_cconj(p2) ,  (*gf)[2][0]    );\
  ch.re = (*gf)[0][0].re*(*gf)[1][0].re + (*gf)[0][0].im*(*gf)[1][0].im;\
  ch.im = -(*gf)[0][0].re*(*gf)[1][0].im + (*gf)[0][0].im*(*gf)[1][0].re;\
  (*gf)[2][2].re = p1.re*ch.re + p1.im*ch.im;\
  (*gf)[2][2].im = p1.re*ch.im - p1.im*ch.re;\
  (*gf)[1][2] = dev_csub((*gf)[1][2], (*gf)[2][2]);\
  (*gf)[1][2] = dev_cconj(dev_crealmult((*gf)[1][2], one_over_N));\
  (*gf)[2][2] = dev_cmult(  dev_cconj(p2) ,   (*gf)[1][0]    );\
  ch.re = (*gf)[0][0].re*(*gf)[2][0].re + (*gf)[0][0].im*(*gf)[2][0].im;\
  ch.im = (*gf)[0][0].im*(*gf)[2][0].re - (*gf)[0][0].re*(*gf)[2][0].im;\
  p2.re = p1.re*ch.re + p1.im*ch.im;\
  p2.im = p1.re*ch.im - p1.im*ch.re;\
  (*gf)[2][2] = dev_cadd((*gf)[2][2], p2);\
  (*gf)[2][2] = dev_cconj(dev_crealmult((*gf)[2][2], -one_over_N));}



#define dev_reconstructgf_2vtexref(field, pos, mu, VOL, gf){ \
  float4 gfin; \
  int COL = 3; \
  gfin = field[pos+VOL*(0 + COL*mu)]; \
  (*gf)[0][0].re = gfin.x; \
  (*gf)[0][0].im = gfin.y; \
  (*gf)[0][1].re = gfin.z; \
  (*gf)[0][1].im = gfin.w; \
  gfin = field[pos+VOL*(1 + COL*mu)]; \
  (*gf)[0][2].re = gfin.x; \
  (*gf)[0][2].im = gfin.y; \
  (*gf)[1][0].re = gfin.z; \
  (*gf)[1][0].im = gfin.w; \
  gfin = field[pos+VOL*(2 + COL*mu)]; \
  (*gf)[1][1].re = gfin.x; \
  (*gf)[1][1].im = gfin.y; \
  (*gf)[1][2].re = gfin.z; \
  (*gf)[1][2].im = gfin.w; \
  (*gf)[2][0].re = (*gf)[0][1].re * (*gf)[1][2].re; \
  (*gf)[2][0].re -= (*gf)[0][1].im * (*gf)[1][2].im; \
  (*gf)[2][0].re -= (*gf)[0][2].re * (*gf)[1][1].re; \
  (*gf)[2][0].re += (*gf)[0][2].im * (*gf)[1][1].im; \
  (*gf)[2][0].im = -(*gf)[0][1].re * (*gf)[1][2].im; \
  (*gf)[2][0].im -= (*gf)[0][1].im * (*gf)[1][2].re; \
  (*gf)[2][0].im += (*gf)[0][2].re * (*gf)[1][1].im; \
  (*gf)[2][0].im += (*gf)[0][2].im * (*gf)[1][1].re; \
  (*gf)[2][1].re = (*gf)[0][2].re * (*gf)[1][0].re; \
  (*gf)[2][1].re -= (*gf)[0][2].im * (*gf)[1][0].im; \
  (*gf)[2][1].re -= (*gf)[0][0].re * (*gf)[1][2].re; \
  (*gf)[2][1].re += (*gf)[0][0].im * (*gf)[1][2].im; \
  (*gf)[2][1].im = -(*gf)[0][2].re * (*gf)[1][0].im; \
  (*gf)[2][1].im -= (*gf)[0][2].im * (*gf)[1][0].re; \
  (*gf)[2][1].im += (*gf)[0][0].re * (*gf)[1][2].im; \
  (*gf)[2][1].im += (*gf)[0][0].im * (*gf)[1][2].re; \
  (*gf)[2][2].re = (*gf)[0][0].re * (*gf)[1][1].re; \
  (*gf)[2][2].re -= (*gf)[0][0].im * (*gf)[1][1].im; \
  (*gf)[2][2].re -= (*gf)[0][1].re * (*gf)[1][0].re; \
  (*gf)[2][2].re += (*gf)[0][1].im * (*gf)[1][0].im; \
  (*gf)[2][2].im = -(*gf)[0][0].re * (*gf)[1][1].im; \
  (*gf)[2][2].im -= (*gf)[0][0].im * (*gf)[1][1].re; \
  (*gf)[2][2].im += (*gf)[0][1].re * (*gf)[1][0].im; \
  (*gf)[2][2].im += (*gf)[0][1].im * (*gf)[1][0].re;}




#define dev_reconstructgf_2vtexref_dagger(field, pos, mu, VOL, gf){ \
  float4 gfin; \
  int COL = 3; \
  gfin = field[pos+VOL*(0 + COL*mu)]; \
  (*gf)[0][0].re = gfin.x; \
  (*gf)[0][0].im = -gfin.y; \
  (*gf)[1][0].re = gfin.z; \
  (*gf)[1][0].im = -gfin.w; \
  gfin = field[pos+VOL*(1 + COL*mu)]; \
  (*gf)[2][0].re = gfin.x; \
  (*gf)[2][0].im = -gfin.y; \
  (*gf)[0][1].re = gfin.z; \
  (*gf)[0][1].im = -gfin.w; \
  gfin = field[pos+VOL*(2 + COL*mu)]; \
  (*gf)[1][1].re = gfin.x; \
  (*gf)[1][1].im = -gfin.y; \
  (*gf)[2][1].re = gfin.z; \
  (*gf)[2][1].im = -gfin.w; \
  (*gf)[0][2].re = (*gf)[1][0].re * (*gf)[2][1].re; \
  (*gf)[0][2].re -= (*gf)[1][0].im * (*gf)[2][1].im; \
  (*gf)[0][2].re -= (*gf)[2][0].re * (*gf)[1][1].re; \
  (*gf)[0][2].re += (*gf)[2][0].im * (*gf)[1][1].im; \
  (*gf)[0][2].im = -(*gf)[1][0].re* (*gf)[2][1].im; \
  (*gf)[0][2].im -= (*gf)[1][0].im* (*gf)[2][1].re; \
  (*gf)[0][2].im += (*gf)[2][0].re*(*gf)[1][1].im; \
  (*gf)[0][2].im += (*gf)[2][0].im*(*gf)[1][1].re; \
  (*gf)[1][2].re = (*gf)[2][0].re*(*gf)[0][1].re; \
  (*gf)[1][2].re -= (*gf)[2][0].im*(*gf)[0][1].im; \
  (*gf)[1][2].re -= (*gf)[0][0].re*(*gf)[2][1].re; \
  (*gf)[1][2].re += (*gf)[0][0].im*(*gf)[2][1].im; \
  (*gf)[1][2].im = -(*gf)[2][0].re * (*gf)[0][1].im; \
  (*gf)[1][2].im -= (*gf)[2][0].im * (*gf)[0][1].re; \
  (*gf)[1][2].im += (*gf)[0][0].re * (*gf)[2][1].im; \
  (*gf)[1][2].im += (*gf)[0][0].im * (*gf)[2][1].re; \
  (*gf)[2][2].re = (*gf)[0][0].re * (*gf)[1][1].re; \
  (*gf)[2][2].re -= (*gf)[0][0].im * (*gf)[1][1].im; \
  (*gf)[2][2].re -= (*gf)[1][0].re * (*gf)[0][1].re; \
  (*gf)[2][2].re += (*gf)[1][0].im * (*gf)[0][1].im; \
  (*gf)[2][2].im = -(*gf)[0][0].re * (*gf)[1][1].im; \
  (*gf)[2][2].im -= (*gf)[0][0].im * (*gf)[1][1].re; \
  (*gf)[2][2].im += (*gf)[1][0].re * (*gf)[0][1].im; \
  (*gf)[2][2].im += (*gf)[1][0].im * (*gf)[0][1].re;}


#endif //USETEXTURE




#else



// reconstruction of the gf using 8 real parameters as 
// described in the appendix of hep-lat 0911.3191 (M.Clark et al.)
// optimized once
__device__ void dev_reconstructgf_8texref (const dev_su3_2v * field, int pos, int mu, int VOL, dev_su3* gf){

  float4 gfin;
  REAL one_over_N, help;
  dev_complex p1,p2;
  int COL = 2;
  #ifdef USETEXTURE
    gfin = tex1Dfetch(gf_tex,pos+VOL*(0 + COL*mu));
  #else
    gfin = field[pos+VOL*(0 + COL*mu)];
  #endif
  // read a2 a3
  (*gf)[0][1].re = gfin.x;
  (*gf)[0][1].im = gfin.y;
  (*gf)[0][2].re = gfin.z;
  (*gf)[0][2].im = gfin.w;  
 
  p1.re = gfin.x*gfin.x + gfin.y*gfin.y + gfin.z*gfin.z + gfin.w*gfin.w; // use later on
  one_over_N = rsqrtf(p1.re); //reciprocal sqrt

  // read theta_a1, theta_c1, b1
  #ifdef USETEXTURE
    gfin = tex1Dfetch(gf_tex,pos+VOL*(1 + COL*mu));
  #else
    gfin = field[pos+VOL*(1 + COL*mu)];
  #endif
  
  // reconstruct a1 use sqrt instead of sin
  help = 1.0f - p1.re;
  p1.re = sqrtf(help>0?help:0);

  #ifdef HALF
    // we have to multiply by two pi because normalization to -1..1
    gfin.x = gfin.x*pi_float;
    gfin.y = gfin.y*pi_float;
  #endif
  __sincosf(gfin.x, &(*gf)[0][0].im, &(*gf)[0][0].re);
  (*gf)[0][0].re = (*gf)[0][0].re * p1.re;
  (*gf)[0][0].im = (*gf)[0][0].im * p1.re;
  
  
  
  // assign b1
  (*gf)[1][0].re = gfin.z;
  (*gf)[1][0].im = gfin.w;
  
  // p2 = 1/N b1
  p2.re = one_over_N*(*gf)[1][0].re;
  p2.im = one_over_N*(*gf)[1][0].im;  


  // reconstruct c1 use sqrt instead of sin
  help =1.0f - 
              (*gf)[0][0].re * (*gf)[0][0].re - (*gf)[0][0].im * (*gf)[0][0].im - 
              (*gf)[1][0].re * (*gf)[1][0].re - (*gf)[1][0].im * (*gf)[1][0].im;

   p1.re = sqrtf(help>0?help:0);

  __sincosf(gfin.y, &(*gf)[2][0].im, &(*gf)[2][0].re);
  (*gf)[2][0].re = (*gf)[2][0].re * p1.re; 
  (*gf)[2][0].im = (*gf)[2][0].im * p1.re;
   
  
  
  // p1 = 1/N*cconj(c1)
  p1.re = one_over_N*(*gf)[2][0].re;
  p1.im = - one_over_N*(*gf)[2][0].im;
  
  
  
  //use the last reconstructed gf component gf[2][2] (c3) as a help variable for b2,b3 and c2
  //this is in order to save registers and to prevent extra loading and storing from global mem
  // calculate b2
  
  (*gf)[1][1].re = p1.re*(*gf)[0][2].re + p1.im*(*gf)[0][2].im;
  (*gf)[1][1].im = p1.im*(*gf)[0][2].re - p1.re*(*gf)[0][2].im;
  (*gf)[2][2].re = (*gf)[0][0].re * (*gf)[0][1].re + (*gf)[0][0].im * (*gf)[0][1].im;
  (*gf)[2][2].im = (*gf)[0][0].re * (*gf)[0][1].im - (*gf)[0][0].im * (*gf)[0][1].re;
  (*gf)[2][2] = dev_cmult(p2, (*gf)[2][2]);  
  (*gf)[1][1].re = -one_over_N*( (*gf)[1][1].re + (*gf)[2][2].re);
  (*gf)[1][1].im = -one_over_N*((*gf)[1][1].im + (*gf)[2][2].im);
  
  
  // calculate b3
  (*gf)[1][2].re = p1.re*(*gf)[0][1].re + p1.im*(*gf)[0][1].im;
  (*gf)[1][2].im = p1.im*(*gf)[0][1].re - p1.re*(*gf)[0][1].im;
  (*gf)[2][2].re = (*gf)[0][0].re*(*gf)[0][2].re + (*gf)[0][0].im*(*gf)[0][2].im;
  (*gf)[2][2].im = (*gf)[0][0].re*(*gf)[0][2].im - (*gf)[0][0].im*(*gf)[0][2].re;
  (*gf)[2][2] = dev_cmult(p2,(*gf)[2][2]);
  (*gf)[1][2].re = one_over_N*( (*gf)[1][2].re - (*gf)[2][2].re);
  (*gf)[1][2].im = one_over_N*( (*gf)[1][2].im - (*gf)[2][2].im);
  
  // calculate c2
  (*gf)[2][1].re = p2.re*(*gf)[0][2].re - p2.im*(*gf)[0][2].im;
  (*gf)[2][1].im = -p2.re*(*gf)[0][2].im - p2.im*(*gf)[0][2].re;
  (*gf)[2][2].re = (*gf)[0][0].re*(*gf)[0][1].re + (*gf)[0][0].im*(*gf)[0][1].im;
  (*gf)[2][2].im = (*gf)[0][0].re* (*gf)[0][1].im - (*gf)[0][0].im* (*gf)[0][1].re;
  help = (*gf)[2][2].re;
  (*gf)[2][2].re = p1.re*(*gf)[2][2].re + p1.im*(*gf)[2][2].im;
  (*gf)[2][2].im = p1.re*(*gf)[2][2].im - p1.im*help;
  (*gf)[2][1].re = one_over_N*((*gf)[2][1].re - (*gf)[2][2].re);
  (*gf)[2][1].im = one_over_N*((*gf)[2][1].im - (*gf)[2][2].im);
  
  // now we have to use p2 and p1 as a help variable, as this is not 
  // needed any more after the first
  // step
  // calculate c3
  (*gf)[2][2].re = p2.re * (*gf)[0][1].re - p2.im * (*gf)[0][1].im;
  (*gf)[2][2].im = - p2.im*(*gf)[0][1].re - p2.re*(*gf)[0][1].im;

  
  p2.re = (*gf)[0][0].re * (*gf)[0][2].re + (*gf)[0][0].im * (*gf)[0][2].im;
  p2.im = (*gf)[0][0].re * (*gf)[0][2].im - (*gf)[0][0].im * (*gf)[0][2].re;
  p2 = dev_cmult(  dev_cconj(p1) , p2);
  
  (*gf)[2][2] = dev_cadd((*gf)[2][2], p2);
  (*gf)[2][2] = dev_crealmult((*gf)[2][2], -one_over_N);
                      
}








__device__ void dev_reconstructgf_8texref_dagger (const dev_su3_2v* field,int pos, int mu, int VOL, dev_su3* gf){


  float4 gfin;
  REAL one_over_N, help;
  dev_complex p1,p2;
  int COL = 2;
  
  #ifdef USETEXTURE
    gfin = tex1Dfetch(gf_tex,pos+VOL*(0 + COL*mu));
  #else
    gfin = field[pos+VOL*(0 + COL*mu)];
  #endif
  // read a2 a3
  (*gf)[1][0].re = gfin.x;
  (*gf)[1][0].im = -gfin.y;
  (*gf)[2][0].re = gfin.z;
  (*gf)[2][0].im = -gfin.w;  
 
  p1.re = gfin.x*gfin.x + gfin.y*gfin.y + gfin.z*gfin.z + gfin.w*gfin.w; // use later on
  one_over_N = rsqrtf(p1.re);  // reciprocal sqrt

  
  // read theta_a1, theta_c1, b1
  #ifdef USETEXTURE
    gfin = tex1Dfetch(gf_tex,pos+VOL*(1 + COL*mu));
  #else
    gfin = field[pos+VOL*(1 + COL*mu)];
  #endif
  
  // reconstruct a1
  help = 1.0f - p1.re;
  p1.re = sqrtf(help>0 ?help:0 );   

  //(*gf)[0][0].re = p1.re*cosf(gfin.x);
  //(*gf)[0][0].im = -p1.re*sinf(gfin.x);
  
  #ifdef HALF
    // we have to multiply by two pi because normalization to -1..1
    gfin.x = gfin.x*pi_float;
    gfin.y = gfin.y*pi_float;
  #endif  
  
  __sincosf(gfin.x, &(*gf)[0][0].im, &(*gf)[0][0].re);
  (*gf)[0][0].re = (*gf)[0][0].re * p1.re;
  (*gf)[0][0].im = -(*gf)[0][0].im * p1.re;
    

  // assign b1
  (*gf)[0][1].re = gfin.z;
  (*gf)[0][1].im = -gfin.w;
  
  // p2 = 1/N b1
  p2.re = one_over_N*(*gf)[0][1].re;
  p2.im = -one_over_N*(*gf)[0][1].im;  


  // reconstruct c1
  help = 1.0f - 
              (*gf)[0][0].re * (*gf)[0][0].re - (*gf)[0][0].im * (*gf)[0][0].im - 
              (*gf)[0][1].re * (*gf)[0][1].re - (*gf)[0][1].im * (*gf)[0][1].im;
  p1.re = sqrtf(help>0?help:0);

  //(*gf)[0][2].re = p1.re*cosf(gfin.y);
  //(*gf)[0][2].im = -p1.re*sinf(gfin.y);
  
  __sincosf(gfin.y, &(*gf)[0][2].im, &(*gf)[0][2].re);
  (*gf)[0][2].re = (*gf)[0][2].re * p1.re;
  (*gf)[0][2].im = -(*gf)[0][2].im * p1.re;
     
  
  // p1 = 1/N*cconj(c1)
  p1.re = one_over_N*(*gf)[0][2].re;
  p1.im = one_over_N*(*gf)[0][2].im;
  
  //use the last reconstructed gf component gf[2][2] (c3) as a help variable for b2,b3 and c2
  //this is in order to save registers and to prevent extra loading and storing from global mem
  // calculate b2
  (*gf)[1][1].re  = p1.re*(*gf)[2][0].re - p1.im*(*gf)[2][0].im;
  (*gf)[1][1].im  = p1.re*(*gf)[2][0].im + p1.im*(*gf)[2][0].re;
  dev_complex ch;
  ch.re  = (*gf)[0][0].re*(*gf)[1][0].re + (*gf)[0][0].im*(*gf)[1][0].im;
  ch.im  = -(*gf)[0][0].re*(*gf)[1][0].im + (*gf)[0][0].im*(*gf)[1][0].re;
  (*gf)[2][2].re  = p2.re*ch.re - p2.im*ch.im;
  (*gf)[2][2].im  = p2.re*ch.im + p2.im*ch.re;
  (*gf)[1][1] = dev_cadd((*gf)[1][1], (*gf)[2][2]);
  (*gf)[1][1].re = -one_over_N * (*gf)[1][1].re;
  (*gf)[1][1].im = one_over_N * (*gf)[1][1].im;
  
  // calculate b3
  (*gf)[2][1] = dev_cmult(p1,   (*gf)[1][0]    );
  ch.re = (*gf)[0][0].re*(*gf)[2][0].re + (*gf)[0][0].im*(*gf)[2][0].im;
  ch.im = -(*gf)[0][0].re*(*gf)[2][0].im + (*gf)[0][0].im*(*gf)[2][0].re;
  (*gf)[2][2].re = p2.re*ch.re - p2.im*ch.im;
  (*gf)[2][2].im = p2.re*ch.im + p2.im*ch.re;
  (*gf)[2][1] = dev_csub((*gf)[2][1], (*gf)[2][2]);
  (*gf)[2][1] = dev_cconj(dev_crealmult((*gf)[2][1], one_over_N));
  
  // calculate c2
  (*gf)[1][2] = dev_cmult(  dev_cconj(p2) ,  (*gf)[2][0]    );
  ch.re = (*gf)[0][0].re*(*gf)[1][0].re + (*gf)[0][0].im*(*gf)[1][0].im;
  ch.im = -(*gf)[0][0].re*(*gf)[1][0].im + (*gf)[0][0].im*(*gf)[1][0].re;
  (*gf)[2][2].re = p1.re*ch.re + p1.im*ch.im;
  (*gf)[2][2].im = p1.re*ch.im - p1.im*ch.re;
  
  
  (*gf)[1][2] = dev_csub((*gf)[1][2], (*gf)[2][2]);
  (*gf)[1][2] = dev_cconj(dev_crealmult((*gf)[1][2], one_over_N));
  
  // use p2 as help variable after the first step
  // calculate c3
  (*gf)[2][2] = dev_cmult(  dev_cconj(p2) ,   (*gf)[1][0]    );
  ch.re = (*gf)[0][0].re*(*gf)[2][0].re + (*gf)[0][0].im*(*gf)[2][0].im;
  ch.im = (*gf)[0][0].im*(*gf)[2][0].re - (*gf)[0][0].re*(*gf)[2][0].im;
  p2.re = p1.re*ch.re + p1.im*ch.im;
  p2.im = p1.re*ch.im - p1.im*ch.re;
		     
  (*gf)[2][2] = dev_cadd((*gf)[2][2], p2);
  (*gf)[2][2] = dev_cconj(dev_crealmult((*gf)[2][2], -one_over_N));

}



#endif //NEW



#ifdef HALF // for half precision 

// reconstruction of the link fields from two rows of the su3 matrix
// numbers are fetched from texture cache
__device__ void dev_reconstructgf_2vtexref_half (const dev_su3_2v_half* field, int pos, int mu, int VOL, dev_su3* gf){
  float4 gfin;
  int COL = 3;
  
  #ifdef USETEXTURE
    gfin = tex1Dfetch(gf_tex,pos+VOL*(0 + COL*mu));
  #else
    gfin = sh4tofl4(field[pos+VOL*(0 + COL*mu)]);
  #endif
  //first row
  (*gf)[0][0].re = gfin.x;
  (*gf)[0][0].im = gfin.y;
  (*gf)[0][1].re = gfin.z;
  (*gf)[0][1].im = gfin.w;
  #ifdef USETEXTURE
    gfin = tex1Dfetch(gf_tex,pos+VOL*(1 + COL*mu));
  #else
    gfin = sh4tofl4(field[pos+VOL*(1 + COL*mu)]);
  #endif
  (*gf)[0][2].re = gfin.x;
  (*gf)[0][2].im = gfin.y;
  //second row
  (*gf)[1][0].re = gfin.z;
  (*gf)[1][0].im = gfin.w;
  #ifdef USETEXTURE
    gfin = tex1Dfetch(gf_tex,pos+VOL*(2 + COL*mu));
  #else
    gfin = sh4tofl4(field[pos+VOL*(2 + COL*mu)]);
  #endif
  (*gf)[1][1].re = gfin.x;
  (*gf)[1][1].im = gfin.y;
  (*gf)[1][2].re = gfin.z;
  (*gf)[1][2].im = gfin.w;
  
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
__device__ void dev_reconstructgf_2vtexref_dagger_half (const dev_su3_2v_half* field, int pos, int mu, int VOL, dev_su3* gf){
  //dev_complex help1;
  //dev_complex help2;
  float4 gfin;
  int COL = 3;
  
  //first column (minus in im for complex conj.)
  #ifdef USETEXTURE
    gfin = tex1Dfetch(gf_tex,pos+VOL*(0 + COL*mu));
  #else
    gfin = sh4tofl4(field[pos+VOL*(0 + COL*mu)]);
  #endif
  (*gf)[0][0].re = gfin.x;
  (*gf)[0][0].im = -gfin.y;
  (*gf)[1][0].re = gfin.z;
  (*gf)[1][0].im = -gfin.w;
  #ifdef USETEXTURE
    gfin = tex1Dfetch(gf_tex,pos+VOL*(1 + COL*mu));
  #else
    gfin = sh4tofl4(field[pos+VOL*(1 + COL*mu)]);
  #endif
  (*gf)[2][0].re = gfin.x;
  (*gf)[2][0].im = -gfin.y;
  
  //second  column (minus in im for complex conj.)
  (*gf)[0][1].re = gfin.z;
  (*gf)[0][1].im = -gfin.w;
  #ifdef USETEXTURE
    gfin = tex1Dfetch(gf_tex,pos+VOL*(2 + COL*mu));
  #else
    gfin = sh4tofl4(field[pos+VOL*(2 + COL*mu)]);
  #endif
  (*gf)[1][1].re = gfin.x;
  (*gf)[1][1].im = -gfin.y;
  (*gf)[2][1].re = gfin.z;
  (*gf)[2][1].im = -gfin.w;
  


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






// reconstruction of the gf using 8 real parameters as 
// described in the appendix of hep-lat 0911.3191 (M.Clark et al.)
// optimized once
__device__ void dev_reconstructgf_8texref_half (const dev_su3_2v_half * field, int pos, int mu, int VOL, dev_su3* gf){

  float4 gfin;
  REAL one_over_N, help;
  dev_complex p1,p2;
  int COL = 2;
  #ifdef USETEXTURE
    gfin = tex1Dfetch(gf_tex,pos+VOL*(0 + COL*mu));
  #else
    gfin = sh4tofl4(field[pos+VOL*(0 + COL*mu)]);
  #endif
  // read a2 a3
  (*gf)[0][1].re = gfin.x;
  (*gf)[0][1].im = gfin.y;
  (*gf)[0][2].re = gfin.z;
  (*gf)[0][2].im = gfin.w;  
 
  p1.re = gfin.x*gfin.x + gfin.y*gfin.y + gfin.z*gfin.z + gfin.w*gfin.w; // use later on
  one_over_N = rsqrtf(p1.re); //reciprocal sqrt

  // read theta_a1, theta_c1, b1
  #ifdef USETEXTURE
    gfin = tex1Dfetch(gf_tex,pos+VOL*(1 + COL*mu));
  #else
    gfin = sh4tofl4(field[pos+VOL*(1 + COL*mu)]);
  #endif
  
  // reconstruct a1 use sqrt instead of sin
  help = 1.0f - p1.re;
  p1.re = sqrtf(help>0?help:0);

  #ifdef HALF
    // we have to multiply by two pi because normalization to -1..1
    gfin.x = gfin.x*pi_float;
    gfin.y = gfin.y*pi_float;
  #endif
  sincos(gfin.x, &(*gf)[0][0].im, &(*gf)[0][0].re);
  (*gf)[0][0].re = (*gf)[0][0].re * p1.re;
  (*gf)[0][0].im = (*gf)[0][0].im * p1.re;
  
  
  
  // assign b1
  (*gf)[1][0].re = gfin.z;
  (*gf)[1][0].im = gfin.w;
  
  // p2 = 1/N b1
  p2.re = one_over_N*(*gf)[1][0].re;
  p2.im = one_over_N*(*gf)[1][0].im;  


  // reconstruct c1 use sqrt instead of sin
  help =1.0f - 
              (*gf)[0][0].re * (*gf)[0][0].re - (*gf)[0][0].im * (*gf)[0][0].im - 
              (*gf)[1][0].re * (*gf)[1][0].re - (*gf)[1][0].im * (*gf)[1][0].im;
  p1.re = sqrtf(help>0?help:0);

  sincos(gfin.y, &(*gf)[2][0].im, &(*gf)[2][0].re);
  (*gf)[2][0].re = (*gf)[2][0].re * p1.re; 
  (*gf)[2][0].im = (*gf)[2][0].im * p1.re;
   
  
  
  // p1 = 1/N*cconj(c1)
  p1.re = one_over_N*(*gf)[2][0].re;
  p1.im = - one_over_N*(*gf)[2][0].im;
  
  
  
  //use the last reconstructed gf component gf[2][2] (c3) as a help variable for b2,b3 and c2
  //this is in order to save registers and to prevent extra loading and storing from global mem
  // calculate b2
  
  (*gf)[1][1].re = p1.re*(*gf)[0][2].re;
  (*gf)[1][1].re += p1.im*(*gf)[0][2].im;
  (*gf)[1][1].im = p1.im*(*gf)[0][2].re;
  (*gf)[1][1].im -= p1.re*(*gf)[0][2].im;
  
  (*gf)[2][2].re = (*gf)[0][0].re * (*gf)[0][1].re;
  (*gf)[2][2].re += (*gf)[0][0].im * (*gf)[0][1].im;
  
  (*gf)[2][2].im = (*gf)[0][0].re * (*gf)[0][1].im;
  (*gf)[2][2].im -= (*gf)[0][0].im * (*gf)[0][1].re;
  (*gf)[2][2] = dev_cmult(p2, (*gf)[2][2]);
  
  (*gf)[1][1].re = -one_over_N*( (*gf)[1][1].re + (*gf)[2][2].re);
  (*gf)[1][1].im = -one_over_N*((*gf)[1][1].im + (*gf)[2][2].im);
  
  
  
  
  
  // calculate b3
  (*gf)[1][2].re = p1.re*(*gf)[0][1].re;
  (*gf)[1][2].re += p1.im*(*gf)[0][1].im;
  (*gf)[1][2].im = p1.im*(*gf)[0][1].re;
  (*gf)[1][2].im -= p1.re*(*gf)[0][1].im;
  
  (*gf)[2][2].re = (*gf)[0][0].re*(*gf)[0][2].re;
  (*gf)[2][2].re += (*gf)[0][0].im*(*gf)[0][2].im;
  (*gf)[2][2].im = (*gf)[0][0].re*(*gf)[0][2].im;
  (*gf)[2][2].im -= (*gf)[0][0].im*(*gf)[0][2].re;
  (*gf)[2][2] = dev_cmult(p2,(*gf)[2][2]);
  
  (*gf)[1][2].re = one_over_N*( (*gf)[1][2].re - (*gf)[2][2].re);
  (*gf)[1][2].im = one_over_N*( (*gf)[1][2].im - (*gf)[2][2].im);
  
  
  // calculate c2
  (*gf)[2][1].re = p2.re*(*gf)[0][2].re;
  (*gf)[2][1].re -= p2.im*(*gf)[0][2].im;
  (*gf)[2][1].im = -p2.re*(*gf)[0][2].im;
  (*gf)[2][1].im -= p2.im*(*gf)[0][2].re;
  
  

  (*gf)[2][2].re = (*gf)[0][0].re*(*gf)[0][1].re;
  (*gf)[2][2].re += (*gf)[0][0].im*(*gf)[0][1].im;
  (*gf)[2][2].im = (*gf)[0][0].re* (*gf)[0][1].im;
  (*gf)[2][2].im -= (*gf)[0][0].im* (*gf)[0][1].re;
  help = (*gf)[2][2].re;
  (*gf)[2][2].re = p1.re*(*gf)[2][2].re;
  (*gf)[2][2].re += p1.im*(*gf)[2][2].im;
  (*gf)[2][2].im = p1.re*(*gf)[2][2].im - p1.im*help;
  
  
  (*gf)[2][1].re = one_over_N*((*gf)[2][1].re - (*gf)[2][2].re);
  (*gf)[2][1].im = one_over_N*((*gf)[2][1].im - (*gf)[2][2].im);
  
  // now we have to use p2 and p1 as a help variable, as this is not 
  // needed any more after the first
  // step
  // calculate c3
  (*gf)[2][2].re = p2.re * (*gf)[0][1].re;
  (*gf)[2][2].re -= p2.im * (*gf)[0][1].im;
  (*gf)[2][2].im = - p2.im*(*gf)[0][1].re;
  (*gf)[2][2].im -= p2.re*(*gf)[0][1].im;
  
  p2.re = (*gf)[0][0].re * (*gf)[0][2].re;
  p2.re += (*gf)[0][0].im * (*gf)[0][2].im;
  p2.im = (*gf)[0][0].re * (*gf)[0][2].im;
  p2.im -= (*gf)[0][0].im * (*gf)[0][2].re;
  p2 = dev_cmult(  dev_cconj(p1) , p2);
  
  (*gf)[2][2] = dev_cadd((*gf)[2][2], p2);
  (*gf)[2][2] = dev_crealmult((*gf)[2][2], -one_over_N);
                      
}








__device__ void dev_reconstructgf_8texref_dagger_half (const dev_su3_2v_half* field,int pos, int mu, int VOL, dev_su3* gf){


  float4 gfin;
  REAL one_over_N, help;
  dev_complex p1,p2;
  int COL = 2;
  #ifdef USETEXTURE
    gfin = tex1Dfetch(gf_tex,pos+VOL*(0 + COL*mu));
  #else
    gfin = sh4tofl4(field[pos+VOL*(0 + COL*mu)]);
  #endif
  // read a2 a3
  (*gf)[1][0].re = gfin.x;
  (*gf)[1][0].im = -gfin.y;
  (*gf)[2][0].re = gfin.z;
  (*gf)[2][0].im = -gfin.w;  
 
  p1.re = gfin.x*gfin.x + gfin.y*gfin.y + gfin.z*gfin.z + gfin.w*gfin.w; // use later on
  one_over_N = rsqrtf(p1.re);  // reciprocal sqrt

  
  // read theta_a1, theta_c1, b1
  #ifdef USETEXTURE
    gfin = tex1Dfetch(gf_tex,pos+VOL*(1 + COL*mu));
  #else
    gfin = sh4tofl4(field[pos+VOL*(1 + COL*mu)]);
  #endif
  
  // reconstruct a1
  help = 1.0f - p1.re;
  p1.re = sqrtf(help>0?help:0);   

  //(*gf)[0][0].re = p1.re*cosf(gfin.x);
  //(*gf)[0][0].im = -p1.re*sinf(gfin.x);
  
  #ifdef HALF
    // we have to multiply by two pi because normalization to -1..1
    gfin.x = gfin.x*pi_float;
    gfin.y = gfin.y*pi_float;
  #endif  
  
  sincos(gfin.x, &(*gf)[0][0].im, &(*gf)[0][0].re);
  (*gf)[0][0].re = (*gf)[0][0].re * p1.re;
  (*gf)[0][0].im = -(*gf)[0][0].im * p1.re;
    

  // assign b1
  (*gf)[0][1].re = gfin.z;
  (*gf)[0][1].im = -gfin.w;
  
  // p2 = 1/N b1
  p2.re = one_over_N*(*gf)[0][1].re;
  p2.im = -one_over_N*(*gf)[0][1].im;  


  // reconstruct c1
  help = 1.0f - 
              (*gf)[0][0].re * (*gf)[0][0].re - (*gf)[0][0].im * (*gf)[0][0].im - 
              (*gf)[0][1].re * (*gf)[0][1].re - (*gf)[0][1].im * (*gf)[0][1].im;
  p1.re = sqrtf(help>0?help:0);

  //(*gf)[0][2].re = p1.re*cosf(gfin.y);
  //(*gf)[0][2].im = -p1.re*sinf(gfin.y);
  
  sincos(gfin.y, &(*gf)[0][2].im, &(*gf)[0][2].re);
  (*gf)[0][2].re = (*gf)[0][2].re * p1.re;
  (*gf)[0][2].im = -(*gf)[0][2].im * p1.re;
     
  
  // p1 = 1/N*cconj(c1)
  p1.re = one_over_N*(*gf)[0][2].re;
  p1.im = one_over_N*(*gf)[0][2].im;
  
  //use the last reconstructed gf component gf[2][2] (c3) as a help variable for b2,b3 and c2
  //this is in order to save registers and to prevent extra loading and storing from global mem
  // calculate b2
  (*gf)[1][1] = dev_cmult(p1,   (*gf)[2][0]    );
  (*gf)[2][2] = dev_cmult(p2, dev_cmult( (*gf)[0][0] , dev_cconj((*gf)[1][0] ))  );
  (*gf)[1][1] = dev_cadd((*gf)[1][1], (*gf)[2][2]);
  (*gf)[1][1] = dev_cconj(dev_crealmult((*gf)[1][1], -one_over_N));
  
  // calculate b3
  (*gf)[2][1] = dev_cmult(p1,   (*gf)[1][0]    );
  (*gf)[2][2] = dev_cmult(p2, dev_cmult( (*gf)[0][0] , dev_cconj((*gf)[2][0] ))  );
  (*gf)[2][1] = dev_csub((*gf)[2][1], (*gf)[2][2]);
  (*gf)[2][1] = dev_cconj(dev_crealmult((*gf)[2][1], one_over_N));
  
  // calculate c2
  (*gf)[1][2] = dev_cmult(  dev_cconj(p2) ,  (*gf)[2][0]    );
  (*gf)[2][2] = dev_cmult(  dev_cconj(p1) , 
                       dev_cmult(   (*gf)[0][0]  , dev_cconj( (*gf)[1][0]) )
                     );
  (*gf)[1][2] = dev_csub((*gf)[1][2], (*gf)[2][2]);
  (*gf)[1][2] = dev_cconj(dev_crealmult((*gf)[1][2], one_over_N));
  
  // use p2 as help variable after the first step
  // calculate c3
  (*gf)[2][2] = dev_cmult(  dev_cconj(p2) ,   (*gf)[1][0]    );
  p2 = dev_cmult(  dev_cconj(p1) , 
                       dev_cmult(   (*gf)[0][0]  , dev_cconj((*gf)[2][0] ) )
                     );
  (*gf)[2][2] = dev_cadd((*gf)[2][2], p2);
  (*gf)[2][2] = dev_cconj(dev_crealmult((*gf)[2][2], -one_over_N));

}


#endif // HALF









// !!!!!!!!!!!!!!!!!!!  BEWARE  !!!!!!!!!!!!!!!!!!!!!!!!!!
// NOT WORKING ANY MORE -> have to change indexing !!
__global__ void dev_check_gauge_reconstruction_8(dev_su3_2v* gf, int pos, dev_su3 * outgf1, dev_su3* outgf2){
  //dev_reconstructgf_8texref (gf,pos, outgf1);
  //dev_reconstructgf_8texref_dagger (gf,pos, outgf2);
}






////////////////////// HOST FUNCTIONS FOR GAUGE RECONSTRUCTION ////////////////////


// get 2 first rows of gf float4 type
//  
//
void su3to2vf4(su3** gf, dev_su3_2v* h2d_gf){
  int i,j;
  int pos=0;
  int VOL;
  int COL=3; // 3 float4 per gauge link
  #ifndef MPI
    VOL = VOLUME;
    for (i = 0; i < VOLUME; i++) {  
  #else
    VOL = (VOLUME+RAND);
    for (i = 0; i < (VOLUME+RAND); i++) {  
  #endif
   for(j=0;j<4;j++){
   //first row
    h2d_gf[pos+VOL*(0 + COL*j)].x = (REAL) creal(gf[i][j].c00);
    h2d_gf[pos+VOL*(0 + COL*j)].y = (REAL) cimag(gf[i][j].c00);
    h2d_gf[pos+VOL*(0 + COL*j)].z = (REAL) creal(gf[i][j].c01);
    h2d_gf[pos+VOL*(0 + COL*j)].w = (REAL) cimag(gf[i][j].c01);
    h2d_gf[pos+VOL*(1 + COL*j)].x = (REAL) creal(gf[i][j].c02);
    h2d_gf[pos+VOL*(1 + COL*j)].y = (REAL) cimag(gf[i][j].c02);      
   //second row
    h2d_gf[pos+VOL*(1 + COL*j)].z = (REAL) creal(gf[i][j].c10);
    h2d_gf[pos+VOL*(1 + COL*j)].w = (REAL) cimag(gf[i][j].c10);
    h2d_gf[pos+VOL*(2 + COL*j)].x = (REAL) creal(gf[i][j].c11);
    h2d_gf[pos+VOL*(2 + COL*j)].y = (REAL) cimag(gf[i][j].c11);
    h2d_gf[pos+VOL*(2 + COL*j)].z = (REAL) creal(gf[i][j].c12);
    h2d_gf[pos+VOL*(2 + COL*j)].w = (REAL) cimag(gf[i][j].c12);      
  }
  pos++;
 }
}




// bring gf into the form
// a2 a3, theta_a1, theta_c1, b1
// 
void su3to8(su3** gf, dev_su3_8* h2d_gf){
  int i,j;
  int pos=0;
  int VOL;
  int COL=2; // 2 float4 per gauge link
  #ifndef MPI
    VOL = VOLUME;
    for (i = 0; i < VOLUME; i++) {
  #else
    VOL = (VOLUME+RAND);
    for (i = 0; i < (VOLUME+RAND); i++) {  
  #endif
      for(j=0;j<4;j++){
      // a2, a3 
      h2d_gf[pos+VOL*(0 + COL*j)].x = (REAL) creal(gf[i][j].c01);
      h2d_gf[pos+VOL*(0 + COL*j)].y = (REAL) cimag(gf[i][j].c01);
      h2d_gf[pos+VOL*(0 + COL*j)].z = (REAL) creal(gf[i][j].c02);
      h2d_gf[pos+VOL*(0 + COL*j)].w = (REAL) cimag(gf[i][j].c02);
      
      // theta_a1, theta_c1
      // use atan2 for this: following the reference, atan2 should give an angle -pi < phi < +pi  
      h2d_gf[pos+VOL*(1 + COL*j)].x = (REAL)( atan2((REAL) cimag(gf[i][j].c00),(REAL) creal(gf[i][j].c00) ));
      h2d_gf[pos+VOL*(1 + COL*j)].y = (REAL) ( atan2((REAL) cimag(gf[i][j].c20),(REAL) creal(gf[i][j].c20) ));
      
      // b1
      h2d_gf[pos+VOL*(1 + COL*j)].z = (REAL) creal(gf[i][j].c10) ;
      h2d_gf[pos+VOL*(1 + COL*j)].w = (REAL) cimag(gf[i][j].c10) ;
     
     } 
   pos++;
 }
}





/*
//FIXME!!
//below code cannot be compiled any more due to c99 complex used now

// this is to reconstruct the gf on the host from 2 rows of the link
// may be used for tests
void reconstructgf_2v (dev_su3* gf){
  dev_complex help1;
  dev_complex help2;
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




// this is to reconstruct the gf on the host from 8 numbers
// !!!!!!!!!!!!!!!!!!!  BEWARE  !!!!!!!!!!!!!!!!!!!!!!!!!!
// NOT WORKING ANY MORE -> have to change indexing !!
void reconstructgf_8 (dev_su3_8 * h2d_gf, dev_su3* gf){

  float4 gfin;
  REAL N, one_over_N, help;
  complex p1,p2, chelp1, chelp2, chelp3, chelpconj, chelpconj2;
  
  gfin = h2d_gf[0];
  // read a2 a3
  (*gf)[0][1].re = gfin.x;
  (*gf)[0][1].im = gfin.y;
  (*gf)[0][2].re = gfin.z;
  (*gf)[0][2].im = gfin.w;  
 
  help = gfin.x*gfin.x + gfin.y*gfin.y + gfin.z*gfin.z + gfin.w*gfin.w; // use later on
  N = sqrt(help);
  one_over_N = 1.0f/N;
  
  // read theta_a1, theta_c1, b1
  gfin = h2d_gf[1];
  
  // reconstruct a1
  help = sqrt(1.0f - help);
  (*gf)[0][0].re = help*cos(gfin.x);
  (*gf)[0][0].im = help*sin(gfin.x);
  
  // assign b1
  (*gf)[1][0].re = gfin.z;
  (*gf)[1][0].im = gfin.w;
  
  // p2 = 1/N b1
  p2.re = one_over_N*(*gf)[1][0].re;
  p2.im = one_over_N*(*gf)[1][0].im;  


  // reconstruct c1
  help = sqrt(1.0f - 
              (*gf)[0][0].re * (*gf)[0][0].re - (*gf)[0][0].im * (*gf)[0][0].im - 
              (*gf)[1][0].re * (*gf)[1][0].re - (*gf)[1][0].im * (*gf)[1][0].im
          );
  (*gf)[2][0].re = help*cos(gfin.y);
  (*gf)[2][0].im = help*sin(gfin.y);

  
  // p1 = 1/N*cconj(c1)
  p1.re = one_over_N*(*gf)[2][0].re;
  p1.im = - one_over_N*(*gf)[2][0].im;
  
  
  float temp = p1.re*p1.re + p1.im*p1.im + p2.re*p2.re + p2.im*p2.im;
  printf("p1**2 + p2**2 = %f\n", temp);
  
  
  // calculate b2
  _complex_conj(chelpconj, (*gf)[0][2] );
  _mult_assign_complex(chelp1, p1, chelpconj   );
  _complex_conj(chelpconj, (*gf)[0][0]);
  _mult_assign_complex(chelp3, chelpconj , (*gf)[0][1] ); 
  _mult_assign_complex(chelp2, p2, chelp3);
  _add_complex(chelp1, chelp2);
  _mult_real((*gf)[1][1], chelp1, -one_over_N);

  
  // calculate b3
  _complex_conj(chelpconj, (*gf)[0][1] );
  _mult_assign_complex(chelp1, p1,  chelpconj  );
  _complex_conj(chelpconj, (*gf)[0][0]);
  _mult_assign_complex(chelp3, chelpconj  , (*gf)[0][2] );  
  _mult_assign_complex(chelp2, p2, chelp3 );
  _diff_complex(chelp1, chelp2);
  _mult_real((*gf)[1][2],chelp1, one_over_N);

  
  // calculate c2
  _complex_conj(chelpconj, p2);
  _complex_conj(chelpconj2, (*gf)[0][2]);
  _mult_assign_complex(chelp1, chelpconj , chelpconj2 );
  _complex_conj(chelpconj,(*gf)[0][0]);
  _mult_assign_complex(chelp3, chelpconj  , (*gf)[0][1] );
  _complex_conj(chelpconj2,p1);
  _mult_assign_complex(chelp2, chelpconj2  , chelp3);
  _diff_complex(chelp1, chelp2);
  _mult_real((*gf)[2][1],chelp1, one_over_N);
  
  
  // calculate c3
  _complex_conj(chelpconj, p2);
  _complex_conj(chelpconj2, (*gf)[0][1] );
  _mult_assign_complex(chelp1, chelpconj  , chelpconj2   );
  _complex_conj(chelpconj,(*gf)[0][0]);
  _mult_assign_complex(chelp3, chelpconj ,(*gf)[0][2]);
  _complex_conj(chelpconj,p1);
  _mult_assign_complex( chelp2, chelpconj  , chelp3 );
  _add_complex(chelp1, chelp2);
  _mult_real((*gf)[2][2], chelp1, -one_over_N);
                  
}








void show_su3(su3 gf1){
   printf("(%f,%f)\t(%f,%f)\t(%f,%f)\n",gf1.c00.re,
   					gf1.c00.im,
   					gf1.c01.re,
   					gf1.c01.im,
   					gf1.c02.re,
   					gf1.c02.im
   );
   printf("(%f,%f)\t(%f,%f)\t(%f,%f)\n",gf1.c10.re,
   					gf1.c10.im,
   					gf1.c11.re,
   					gf1.c11.im,
   					gf1.c12.re,
   					gf1.c12.im
   );
   printf("(%f,%f)\t(%f,%f)\t(%f,%f)\n",gf1.c20.re,
   					gf1.c20.im,
   					gf1.c21.re,
   					gf1.c21.im,
   					gf1.c22.re,
   					gf1.c22.im
   ); 
}


void show_dev_su3(dev_su3 gf1){
   printf("(%f,%f)\t(%f,%f)\t(%f,%f)\n",gf1[0][0].re,
   					gf1[0][0].im,
   					gf1[0][1].re,
   					gf1[0][1].im,
   					gf1[0][2].re,
   					gf1[0][2].im
   );   
   printf("(%f,%f)\t(%f,%f)\t(%f,%f)\n",gf1[1][0].re,
   					gf1[1][0].im,
   					gf1[1][1].re,
   					gf1[1][1].im,
   					gf1[1][2].re,
   					gf1[1][2].im
   );
   printf("(%f,%f)\t(%f,%f)\t(%f,%f)\n",gf1[2][0].re,
   					gf1[2][0].im,
   					gf1[2][1].re,
   					gf1[2][1].im,
   					gf1[2][2].re,
   					gf1[2][2].im
   ); 

}



// !!!!!!!!!!!!!!!!!!!  BEWARE  !!!!!!!!!!!!!!!!!!!!!!!!!!
// NOT WORKING ANY MORE -> have to change indexing !!
void check_gauge_reconstruction_8(su3 ** gf1, dev_su3_2v * gf2, int ind1, int mu){
  dev_su3 * reconst_g , * reconst_g_dagger;
  dev_su3  result, result_dagger;
   printf("Checking 8 paramater reconstruction of gauge field:\n");  
  su3 gfdagger;
  #ifdef USETEXTURE
    bind_texture_gf(gf2);
  #endif
  printf("\n");
  size_t cpsize = sizeof(dev_su3); // parallel in t and z direction
  cudaMalloc((void **) &reconst_g, cpsize); 
  cudaMalloc((void **) &reconst_g_dagger, cpsize); 
  
  show_su3(gf1[ind1][mu]);
  printf("\n");
  
  dev_check_gauge_reconstruction_8  <<< 1 , 1 >>> (dev_gf,4*ind1 + mu, reconst_g, reconst_g_dagger);
  cudaMemcpy(&result, reconst_g, cpsize, cudaMemcpyDeviceToHost);
  cudaMemcpy(&result_dagger, reconst_g_dagger, cpsize, cudaMemcpyDeviceToHost);

  show_dev_su3(result);
  printf("\n");
  
  _su3_dagger(gfdagger,gf1[ind1][mu]);
  show_su3(gfdagger);
  printf("\n");
  show_dev_su3(result_dagger);



  #ifdef USETEXTURE
   unbind_texture_gf();
  #endif
  cudaFree(reconst_g);
}










// !!!!!!!!!!!!!!!!!!!  BEWARE  !!!!!!!!!!!!!!!!!!!!!!!!!!
// NOT WORKING ANY MORE -> have to change indexing !!
// compare host gauge-field with gauge-field that is reconstructed (on host)
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
#ifdef GF_8
   printf("8-field:\t(%f,%f,%f,%f) (%f,%f,%f,%f)\n",
     h2d_gf[2*(4*ind2+mu)].x,
     h2d_gf[2*(4*ind2+mu)].y,
     h2d_gf[2*(4*ind2+mu)].z,
     h2d_gf[2*(4*ind2+mu)].w,
     h2d_gf[2*(4*ind2+mu)+1].x,
     h2d_gf[2*(4*ind2+mu)+1].y,
     h2d_gf[2*(4*ind2+mu)+1].z,
     h2d_gf[2*(4*ind2+mu)+1].w
   );
   dev_su3 help; 
   reconstructgf_8( &(h2d_gf[2*(4*ind2+mu)]) , &help );
   printf("(%f,%f)\t(%f,%f)\t(%f,%f)\n",help[0][0].re,
   					help[0][0].im,
   					help[0][1].re,
   					help[0][1].im,
   					help[0][2].re,
   					help[0][2].im
   );   
   printf("(%f,%f)\t(%f,%f)\t(%f,%f)\n",help[1][0].re,
   					help[1][0].im,
   					help[1][1].re,
   					help[1][1].im,
   					help[1][2].re,
   					help[1][2].im
   );
   printf("(%f,%f)\t(%f,%f)\t(%f,%f)\n",help[2][0].re,
   					help[2][0].im,
   					help[2][1].re,
   					help[2][1].im,
   					help[2][2].re,
   					help[2][2].im
   );   
   
#else
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
#endif 
}




*/








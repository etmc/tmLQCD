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
template<class RealT>
__device__ void dev_reconstructgf_2vtexref (const typename dev_su3_2vT<RealT>::type* field, int pos, typename dev_su3T<RealT>::type* gf){
  typename REAL4T<RealT>::type gfin;
  
  #ifdef USETEXTURE
    gfin = tex1Dfetch(gf_tex,3*pos);
  #else
    gfin = field[3*pos];
  #endif
  //first row
  (*gf)[0][0].re = gfin.x;
  (*gf)[0][0].im = gfin.y;
  (*gf)[0][1].re = gfin.z;
  (*gf)[0][1].im = gfin.w;
  #ifdef USETEXTURE
    gfin = tex1Dfetch(gf_tex,3*pos+1);
  #else
    gfin = field[3*pos + 1];
  #endif
  (*gf)[0][2].re = gfin.x;
  (*gf)[0][2].im = gfin.y;
  //second row
  (*gf)[1][0].re = gfin.z;
  (*gf)[1][0].im = gfin.w;
  #ifdef USETEXTURE
    gfin = tex1Dfetch(gf_tex,3*pos+2);
  #else
    gfin = field[3*pos + 2];
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
template<class RealT>
__device__ void dev_reconstructgf_2vtexref_dagger (const typename dev_su3_2vT<RealT>::type* field, int pos, typename dev_su3T<RealT>::type* gf){
  //dev_complex help1;
  //dev_complex help2;
  typename REAL4T<RealT>::type gfin;
  
  
  //first column (minus in im for complex conj.)
  #ifdef USETEXTURE
    gfin = tex1Dfetch(gf_tex,3*pos);
  #else
    gfin = field[3*pos];
  #endif
  (*gf)[0][0].re = gfin.x;
  (*gf)[0][0].im = -gfin.y;
  (*gf)[1][0].re = gfin.z;
  (*gf)[1][0].im = -gfin.w;
  #ifdef USETEXTURE
    gfin = tex1Dfetch(gf_tex,3*pos+1);
  #else
    gfin = field[3*pos +1];
  #endif
  (*gf)[2][0].re = gfin.x;
  (*gf)[2][0].im = -gfin.y;
  
  //second  column (minus in im for complex conj.)
  (*gf)[0][1].re = gfin.z;
  (*gf)[0][1].im = -gfin.w;
  #ifdef USETEXTURE
    gfin = tex1Dfetch(gf_tex,3*pos+2);
  #else
    gfin = field[3*pos +2];
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
template<class RealT>
__device__ void dev_reconstructgf_8texref (const typename dev_su3_2vT<RealT>::type* field, int pos, typename dev_su3T<RealT>::type* gf){

  typename REAL4T<RealT>::type gfin;
  RealT one_over_N, help;
  dev_complexT<RealT> p1,p2;
  
  #ifdef USETEXTURE
    gfin = tex1Dfetch(gf_tex,2*pos);
  #else
    gfin = field[2*pos];
  #endif
  // read a2 a3
  (*gf)[0][1].re = gfin.x;
  (*gf)[0][1].im = gfin.y;
  (*gf)[0][2].re = gfin.z;
  (*gf)[0][2].im = gfin.w;  
 
  p1.re = gfin.x*gfin.x + gfin.y*gfin.y + gfin.z*gfin.z + gfin.w*gfin.w; // use later on
  one_over_N = rsqrt(p1.re); //reciprocal sqrt

  // read theta_a1, theta_c1, b1
  #ifdef USETEXTURE
    gfin = tex1Dfetch(gf_tex,2*pos + 1);
  #else
    gfin = field[2*pos + 1];
  #endif
  
  // reconstruct a1 use sqrt instead of sin
  help = 1.0f - p1.re;
  if(help > 0.0f){
     p1.re = sqrt(help);
  }
  else{
    p1.re = 0.0f;
  }
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
  if(help > 0.0f){
    p1.re = sqrt(help);
  }   
  else{      
    p1.re = 0.0f;  
  }
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







template<class RealT>
__device__ void dev_reconstructgf_8texref_dagger (const typename dev_su3_2vT<RealT>::type* field,int pos, typename dev_su3T<RealT>::type* gf){


  typename REAL4T<RealT>::type gfin;
  RealT one_over_N, help;
  dev_complexT<RealT> p1,p2;
  
  #ifdef USETEXTURE
    gfin = tex1Dfetch(gf_tex,2*pos);
  #else
    gfin = field[2*pos];
  #endif
  // read a2 a3
  (*gf)[1][0].re = gfin.x;
  (*gf)[1][0].im = -gfin.y;
  (*gf)[2][0].re = gfin.z;
  (*gf)[2][0].im = -gfin.w;  
 
  p1.re = gfin.x*gfin.x + gfin.y*gfin.y + gfin.z*gfin.z + gfin.w*gfin.w; // use later on
  one_over_N = rsqrt(p1.re);  // reciprocal sqrt

  
  // read theta_a1, theta_c1, b1
  #ifdef USETEXTURE
    gfin = tex1Dfetch(gf_tex,2*pos + 1);
  #else
    gfin = field[2*pos + 1];
  #endif
  
  // reconstruct a1
  help = 1.0f - p1.re;
  if(help > 0.0f){
     p1.re = sqrt(help);   
  }
  else{
    p1.re = 0.0f;
  }
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
  if(help > 0.0f){
    p1.re = sqrt(help);
  }
  else{
    p1.re = 0.0f;
  }
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




#ifdef HALF // for half precision 

// reconstruction of the link fields from two rows of the su3 matrix
// numbers are fetched from texture cache
__device__ void dev_reconstructgf_2vtexref_half (const dev_su3_2v_half* field, int pos, dev_su3* gf){
  float4 gfin;
  
  #ifdef USETEXTURE
    gfin = tex1Dfetch(gf_tex,3*pos);
  #else
    gfin = sh4tofl4(field[3*pos]);
  #endif
  //first row
  (*gf)[0][0].re = gfin.x;
  (*gf)[0][0].im = gfin.y;
  (*gf)[0][1].re = gfin.z;
  (*gf)[0][1].im = gfin.w;
  #ifdef USETEXTURE
    gfin = tex1Dfetch(gf_tex,3*pos+1);
  #else
    gfin = sh4tofl4(field[3*pos + 1]);
  #endif
  (*gf)[0][2].re = gfin.x;
  (*gf)[0][2].im = gfin.y;
  //second row
  (*gf)[1][0].re = gfin.z;
  (*gf)[1][0].im = gfin.w;
  #ifdef USETEXTURE
    gfin = tex1Dfetch(gf_tex,3*pos+2);
  #else
    gfin = sh4tofl4(field[3*pos + 2]);
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
__device__ void dev_reconstructgf_2vtexref_dagger_half (const dev_su3_2v_half* field, int pos, dev_su3* gf){
  //dev_complex help1;
  //dev_complex help2;
  float4 gfin;
  
  
  //first column (minus in im for complex conj.)
  #ifdef USETEXTURE
    gfin = tex1Dfetch(gf_tex,3*pos);
  #else
    gfin = sh4tofl4(field[3*pos]);
  #endif
  (*gf)[0][0].re = gfin.x;
  (*gf)[0][0].im = -gfin.y;
  (*gf)[1][0].re = gfin.z;
  (*gf)[1][0].im = -gfin.w;
  #ifdef USETEXTURE
    gfin = tex1Dfetch(gf_tex,3*pos+1);
  #else
    gfin = sh4tofl4(field[3*pos +1]);
  #endif
  (*gf)[2][0].re = gfin.x;
  (*gf)[2][0].im = -gfin.y;
  
  //second  column (minus in im for complex conj.)
  (*gf)[0][1].re = gfin.z;
  (*gf)[0][1].im = -gfin.w;
  #ifdef USETEXTURE
    gfin = tex1Dfetch(gf_tex,3*pos+2);
  #else
    gfin = sh4tofl4(field[3*pos +2]);
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
__device__ void dev_reconstructgf_8texref_half (const dev_su3_2v_half * field, int pos, dev_su3* gf){

  float4 gfin;
  REAL one_over_N, help;
  dev_complex p1,p2;
  
  #ifdef USETEXTURE
    gfin = tex1Dfetch(gf_tex,2*pos);
  #else
    gfin = sh4tofl4(field[2*pos]);
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
    gfin = tex1Dfetch(gf_tex,2*pos + 1);
  #else
    gfin = sh4tofl4(field[2*pos + 1]);
  #endif
  
  // reconstruct a1 use sqrt instead of sin
  help = 1.0f - p1.re;
  if(help > 0.0f){
     p1.re = sqrtf(help);
  }
  else{
    p1.re = 0.0f;
  }
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
  if(help > 0.0f){
    p1.re = sqrtf(help);
  }   
  else{      
    p1.re = 0.0f;  
  }
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








__device__ void dev_reconstructgf_8texref_dagger_half (const dev_su3_2v_half* field,int pos, dev_su3* gf){


  float4 gfin;
  REAL one_over_N, help;
  dev_complex p1,p2;
  
  #ifdef USETEXTURE
    gfin = tex1Dfetch(gf_tex,2*pos);
  #else
    gfin = sh4tofl4(field[2*pos]);
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
    gfin = tex1Dfetch(gf_tex,2*pos + 1);
  #else
    gfin = sh4tofl4(field[2*pos + 1]);
  #endif
  
  // reconstruct a1
  help = 1.0f - p1.re;
  if(help > 0.0f){
     p1.re = sqrtf(help);   
  }
  else{
    p1.re = 0.0f;
  }
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
  if(help > 0.0f){
    p1.re = sqrtf(help);
  }
  else{
    p1.re = 0.0f;
  }
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









template<class RealT>
__global__ void dev_check_gauge_reconstruction_8(typename dev_su3_2vT<RealT>::type* gf, int pos, typename dev_su3T<RealT>::type * outgf1, typename dev_su3T<RealT>::type* outgf2){
  dev_reconstructgf_8texref        <RealT>(gf,pos, outgf1);
  dev_reconstructgf_8texref_dagger <RealT>(gf,pos, outgf2);
}






////////////////////// HOST FUNCTIONS FOR GAUGE RECONSTRUCTION ////////////////////


// get 2 first rows of gf float4 type
//  
//
template<class RealT>
void su3to2vf4(su3** gf, typename dev_su3_2vT<RealT>::type* h2d_gf){
  int i,j;
  #ifndef MPI
    for (i = 0; i < VOLUME; i++) {
  #else
    for (i = 0; i < (VOLUME+RAND); i++) {
  #endif
   for(j=0;j<4;j++){
   //first row
    h2d_gf[3*(4*i+j)].x = (RealT) gf[i][j].c00.re;
    h2d_gf[3*(4*i+j)].y = (RealT) gf[i][j].c00.im;
    h2d_gf[3*(4*i+j)].z = (RealT) gf[i][j].c01.re;
    h2d_gf[3*(4*i+j)].w = (RealT) gf[i][j].c01.im;
    h2d_gf[3*(4*i+j)+1].x = (RealT) gf[i][j].c02.re;
    h2d_gf[3*(4*i+j)+1].y = (RealT) gf[i][j].c02.im;      
   //second row
    h2d_gf[3*(4*i+j)+1].z = (RealT) gf[i][j].c10.re;
    h2d_gf[3*(4*i+j)+1].w = (RealT) gf[i][j].c10.im;
    h2d_gf[3*(4*i+j)+2].x = (RealT) gf[i][j].c11.re;
    h2d_gf[3*(4*i+j)+2].y = (RealT) gf[i][j].c11.im;
    h2d_gf[3*(4*i+j)+2].z = (RealT) gf[i][j].c12.re;
    h2d_gf[3*(4*i+j)+2].w = (RealT) gf[i][j].c12.im;      
  } 
 }
}




// bring gf into the form
// a2 a3, theta_a1, theta_c1, b1
// 
template<class RealT>
void su3to8(su3** gf, typename dev_su3_8T<RealT>::type* h2d_gf){
  int i,j;
  #ifndef MPI
    for (i = 0; i < VOLUME; i++) {
  #else
    for (i = 0; i < (VOLUME+RAND); i++) {
  #endif
      for(j=0;j<4;j++){
      // a2, a3
      h2d_gf[2*(4*i+j)].x = (RealT) gf[i][j].c01.re;
      h2d_gf[2*(4*i+j)].y = (RealT) gf[i][j].c01.im;
      h2d_gf[2*(4*i+j)].z = (RealT) gf[i][j].c02.re;
      h2d_gf[2*(4*i+j)].w = (RealT) gf[i][j].c02.im;
      
      // theta_a1, theta_c1
      // use atan2 for this: following the reference, atan2 should give an angle -pi < phi < +pi  
      h2d_gf[2*(4*i+j)+1].x = (RealT)( atan2((RealT) gf[i][j].c00.im,(RealT) gf[i][j].c00.re ));
      h2d_gf[2*(4*i+j)+1].y = (RealT) ( atan2((RealT) gf[i][j].c20.im,(RealT)gf[i][j].c20.re ));
      
      // b1
      h2d_gf[2*(4*i+j)+1].z = (RealT) gf[i][j].c10.re ;
      h2d_gf[2*(4*i+j)+1].w = (RealT) gf[i][j].c10.im ;
  } 
 }
}






// this is to reconstruct the gf on the host from 2 rows of the link
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




// this is to reconstruct the gf on the host from 2 rows of the link
// may be used for tests
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



template<class RealT>
void check_gauge_reconstruction_8(su3 ** gf1, dev_su3_2vM(RealT) * gf2, int ind1, int mu, MixedsolveParameter<RealT>& mixedsolveParameter){
  dev_su3M(RealT) * reconst_g , * reconst_g_dagger;
  dev_su3M(RealT)  result, result_dagger;
   printf("Checking 8 paramater reconstruction of gauge field:\n");  
  su3 gfdagger;
  #ifdef USETEXTURE
    bind_texture_gf(gf2);
  #endif
  printf("\n");
  size_t cpsize = sizeof(dev_su3M(RealT)); // parallel in t and z direction
  cudaMalloc((void **) &reconst_g, cpsize); 
  cudaMalloc((void **) &reconst_g_dagger, cpsize); 
  
  show_su3(gf1[ind1][mu]);
  printf("\n");
  
  dev_check_gauge_reconstruction_8<REAL>  <<< 1 , 1 >>> (mixedsolveParameter.dev_gf,4*ind1 + mu, reconst_g, reconst_g_dagger);
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











// compare host gauge-field with gauge-field that is reconstructed (on host)
template<class RealT>
void showcompare_gf(int t, int x, int y, int z, int mu, MixedsolveParameter<RealT>& mixedsolveParameter){
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
     mixedsolveParameter.h2d_gf[2*(4*ind2+mu)].x,
     mixedsolveParameter.h2d_gf[2*(4*ind2+mu)].y,
     mixedsolveParameter.h2d_gf[2*(4*ind2+mu)].z,
     mixedsolveParameter.h2d_gf[2*(4*ind2+mu)].w,
     mixedsolveParameter.h2d_gf[2*(4*ind2+mu)+1].x,
     mixedsolveParameter.h2d_gf[2*(4*ind2+mu)+1].y,
     mixedsolveParameter.h2d_gf[2*(4*ind2+mu)+1].z,
     mixedsolveParameter.h2d_gf[2*(4*ind2+mu)+1].w
   );
   dev_su3M(RealT) help; 
   reconstructgf_8( &(mixedsolveParameter.h2d_gf[2*(4*ind2+mu)]) , &help );
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
   printf("(%f,%f)\t(%f,%f)\t(%f,%f)\n",mixedsolveParameter.h2d_gf[3*(4*ind2+mu)].x,
   					mixedsolveParameter.h2d_gf[3*(4*ind2+mu)].y,
   					mixedsolveParameter.h2d_gf[3*(4*ind2+mu)].z,
   					mixedsolveParameter.h2d_gf[3*(4*ind2+mu)].w,
   					mixedsolveParameter.h2d_gf[3*(4*ind2+mu)+1].x,
   					mixedsolveParameter.h2d_gf[3*(4*ind2+mu)+1].y
   );
   printf("(%f,%f)\t(%f,%f)\t(%f,%f)\n",mixedsolveParameter.h2d_gf[3*(4*ind2+mu)+1].z,
   					mixedsolveParameter.h2d_gf[3*(4*ind2+mu)+1].w,
   					mixedsolveParameter.h2d_gf[3*(4*ind2+mu)+2].x,
   					mixedsolveParameter.h2d_gf[3*(4*ind2+mu)+2].y,
   					mixedsolveParameter.h2d_gf[3*(4*ind2+mu)+2].z,
   					mixedsolveParameter.h2d_gf[3*(4*ind2+mu)+2].w
   );   
   
   dev_su3M(RealT) help;
   
   help[0][0].re = mixedsolveParameter.h2d_gf[3*(4*ind2+mu)].x;
   help[0][0].im = mixedsolveParameter.h2d_gf[3*(4*ind2+mu)].y;
   help[0][1].re = mixedsolveParameter.h2d_gf[3*(4*ind2+mu)].z;
   help[0][1].im = mixedsolveParameter.h2d_gf[3*(4*ind2+mu)].w;

   help[0][2].re = mixedsolveParameter.h2d_gf[3*(4*ind2+mu)+1].x;
   help[0][2].im = mixedsolveParameter.h2d_gf[3*(4*ind2+mu)+1].y;
   help[1][0].re = mixedsolveParameter.h2d_gf[3*(4*ind2+mu)+1].z;
   help[1][0].im = mixedsolveParameter.h2d_gf[3*(4*ind2+mu)+1].w;
   
   help[1][1].re = mixedsolveParameter.h2d_gf[3*(4*ind2+mu)+2].x;
   help[1][1].im = mixedsolveParameter.h2d_gf[3*(4*ind2+mu)+2].y;
   help[1][2].re = mixedsolveParameter.h2d_gf[3*(4*ind2+mu)+2].z;
   help[1][2].im = mixedsolveParameter.h2d_gf[3*(4*ind2+mu)+2].w;   
   
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













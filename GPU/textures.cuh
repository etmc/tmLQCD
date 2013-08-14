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
 * File: textures.cuh
 *
 * CUDA texture functions and references
 *
 * 
 *
 **************************************************************************/


#ifdef HAVE_CONFIG_H
  #include<config.h>
#endif
 
 /* texture for nearest neighbours*/
 texture<int,1, cudaReadModeElementType> nn_tex;
 const textureReference* nn_texRefPtr = NULL;
 cudaChannelFormatDesc nn_channelDesc;

 /* textures for spinor field */
 texture<float4,1, cudaReadModeElementType> spin_tex0;
 texture<float4,1, cudaReadModeElementType> spin_tex1;
 texture<float4,1, cudaReadModeElementType> spin_tex2;
 texture<float4,1, cudaReadModeElementType> spin_tex3; 
 texture<float4,1, cudaReadModeElementType> spin_tex4;
 texture<float4,1, cudaReadModeElementType> spin_tex5;
 
  /* textures for down type spinor field */
 texture<float4,1, cudaReadModeElementType> spin_tex_dn0;
 texture<float4,1, cudaReadModeElementType> spin_tex_dn1;
 texture<float4,1, cudaReadModeElementType> spin_tex_dn2;
 texture<float4,1, cudaReadModeElementType> spin_tex_dn3; 
 texture<float4,1, cudaReadModeElementType> spin_tex_dn4;
 texture<float4,1, cudaReadModeElementType> spin_tex_dn5;
 
 texture<float2,1,cudaReadModeElementType> sw_tex;
 texture<float2,1,cudaReadModeElementType> sw_inv_tex;
 
#ifndef HALF
 /* texture for gauge field */
 texture<float4,1, cudaReadModeElementType> gf_tex;
 const textureReference* gf_texRefPtr = NULL;
 cudaChannelFormatDesc gf_channelDesc;
 

extern "C" int bind_texture_spin(dev_spinor* s, int i){
  
  size_t size;
  int offset;
  #ifdef MPI
    if(even_odd_flag){
      size = sizeof(float4)*(VOLUME+RAND)/2;
      offset = (VOLUME+RAND)/2;
    }
    else{
      size = sizeof(float4)*(VOLUME+RAND);
      offset = (VOLUME+RAND);
    }
  #else
    if(even_odd_flag){
      size = sizeof(float4)*VOLUME/2;
      offset = VOLUME/2;
    }
    else{
      size = sizeof(float4)*VOLUME;
      offset = VOLUME;
    }
  #endif
   

  //printf("Binding texture to spinorfield 1\n");
    cudaBindTexture(0, spin_tex0, s, size);
    cudaBindTexture(0, spin_tex1, s+offset, size);
    cudaBindTexture(0, spin_tex2, s+2*offset, size);
    cudaBindTexture(0, spin_tex3, s+3*offset, size);
    cudaBindTexture(0, spin_tex4, s+4*offset, size);
    cudaBindTexture(0, spin_tex5, s+5*offset, size);  
  //printf("%s\n", cudaGetErrorString(cudaGetLastError())); 
  
  return(0);
}


extern "C" int bind_texture_sw(float2* sw){
  
  size_t size;
  int offset;
  #ifdef MPI
      size = sizeof(float2)*6*9*(VOLUME+RAND);
  #else
      size = sizeof(float2)*6*9*VOLUME;
  #endif
   

  //printf("Binding texture to spinorfield 1\n");
    cudaBindTexture(0, sw_tex, sw, size);
  //printf("%s\n", cudaGetErrorString(cudaGetLastError())); 
  
  return(0);
}


extern "C" int unbind_texture_sw(){
  cudaUnbindTexture(sw_tex);
}




extern "C" int bind_texture_sw_inv(float2* sw_inv){
  
  size_t size;
  int offset;
  #ifdef MPI
      size = sizeof(float2)*8*9*(VOLUME+RAND);
  #else
      size = sizeof(float2)*8*9*VOLUME;
  #endif
   

  //printf("Binding texture to spinorfield 1\n");
    cudaBindTexture(0, sw_inv_tex, sw_inv, size);
  //printf("%s\n", cudaGetErrorString(cudaGetLastError())); 
  
  return(0);
}

extern "C" int unbind_texture_sw_inv(){
  cudaUnbindTexture(sw_inv_tex);
}





extern "C" int unbind_texture_spin(int i){
  cudaUnbindTexture(spin_tex0);
  cudaUnbindTexture(spin_tex1);
  cudaUnbindTexture(spin_tex2);
  cudaUnbindTexture(spin_tex3);
  cudaUnbindTexture(spin_tex4);
  cudaUnbindTexture(spin_tex5);
return(1);
}



extern "C" int bind_texture_spin_dn(dev_spinor* s, int i){
  
  size_t size;
  int offset;
  #ifdef MPI
    if(even_odd_flag){
      size = sizeof(float4)*(VOLUME+RAND)/2;
      offset = (VOLUME+RAND)/2;
    }
    else{
      size = sizeof(float4)*(VOLUME+RAND);
      offset = (VOLUME+RAND);
    }
  #else
    if(even_odd_flag){
      size = sizeof(float4)*VOLUME/2;
      offset = VOLUME/2;
    }
    else{
      size = sizeof(float4)*VOLUME;
      offset = VOLUME;
    }
  #endif
   

  //printf("Binding texture to spinorfield 1\n");
    cudaBindTexture(0, spin_tex_dn0, s, size);
    cudaBindTexture(0, spin_tex_dn1, s+offset, size);
    cudaBindTexture(0, spin_tex_dn2, s+2*offset, size);
    cudaBindTexture(0, spin_tex_dn3, s+3*offset, size);
    cudaBindTexture(0, spin_tex_dn4, s+4*offset, size);
    cudaBindTexture(0, spin_tex_dn5, s+5*offset, size);  
  //printf("%s\n", cudaGetErrorString(cudaGetLastError())); 
  
  return(0);
}


extern "C" int unbind_texture_spin_dn(int i){
  cudaUnbindTexture(spin_tex_dn0);
  cudaUnbindTexture(spin_tex_dn1);
  cudaUnbindTexture(spin_tex_dn2);
  cudaUnbindTexture(spin_tex_dn3);
  cudaUnbindTexture(spin_tex_dn4);
  cudaUnbindTexture(spin_tex_dn5);
return(1);
}







#else 

 /* texture for gauge field */
 texture<short4,1, cudaReadModeNormalizedFloat> gf_tex;
 const textureReference* gf_texRefPtr = NULL;
 cudaChannelFormatDesc gf_channelDesc;

  // the textures for the half spinors are defined in half.cuh

#endif // NOT HALF







extern "C" int bind_texture_gf(dev_su3_2v * gf){
 //printf("Binding texture to gaugefield\n");
 
  #ifdef MPI
    #ifdef GF_8
     size_t size = sizeof(float4)*2*(VOLUME+RAND)*4;
    #else
     size_t size = sizeof(float4)*3*(VOLUME+RAND)*4;
    #endif
  #else
    #ifdef GF_8
     size_t size = sizeof(float4)*2*VOLUME*4;
    #else
     size_t size = sizeof(float4)*3*VOLUME*4;
    #endif
  #endif
 /*
 cudaGetTextureReference(&gf_texRefPtr, "gf_tex");
 gf_channelDesc =  cudaCreateChannelDesc<float4>();
 cudaBindTexture(0, gf_texRefPtr, gf, &gf_channelDesc, size);
 */
 cudaBindTexture(0, gf_tex, gf, size);
 //printf("%s\n", cudaGetErrorString(cudaGetLastError()));    
 return(0);
}


extern "C" int unbind_texture_gf(){
 //printf("Unbinding texture to gaugefield\n");
 /*
 cudaUnbindTexture(gf_texRefPtr);
 */
 cudaUnbindTexture(gf_tex);
 //printf("%s\n", cudaGetErrorString(cudaGetLastError()));    
 return(0);
}







extern "C" int bind_texture_nn(int* nn){
 //printf("Binding texture to nn field\n");
  size_t size;
  
  #ifdef MPI
    if(even_odd_flag){
      size = sizeof(int)*8*(VOLUME+RAND)/2;
    }
    else{
      size = sizeof(int)*8*(VOLUME+RAND);
    }
  #else
    if(even_odd_flag){
      size = sizeof(int)*8*VOLUME/2;
    }
    else{
      size = sizeof(int)*8*VOLUME;
    }
  #endif
 

 cudaGetTextureReference(&nn_texRefPtr, "nn_tex");
 nn_channelDesc =  cudaCreateChannelDesc<int>();
 cudaBindTexture(0, nn_texRefPtr, nn, &nn_channelDesc, size);
 //printf("%s\n", cudaGetErrorString(cudaGetLastError()));    
 return(0);
}


extern "C" int unbind_texture_nn(){
 //printf("Unbinding texture to nn field\n");
 cudaUnbindTexture(nn_texRefPtr);
 //printf("%s\n", cudaGetErrorString(cudaGetLastError()));    
 return(0);
}














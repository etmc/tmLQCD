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

 /* texture for spinor field */
 texture<float4,1, cudaReadModeElementType> spin_tex;
 const textureReference* spin_texRefPtr = NULL;
 cudaChannelFormatDesc spin_channelDesc;

 /* texture for spinor field 2*/
 texture<float4,1, cudaReadModeElementType> spin_tex2;
 const textureReference* spin_texRefPtr2 = NULL;
 cudaChannelFormatDesc spin_channelDesc2;


#ifndef HALF
 /* texture for gauge field */
 texture<float4,1, cudaReadModeElementType> gf_tex;
 const textureReference* gf_texRefPtr = NULL;
 cudaChannelFormatDesc gf_channelDesc;
 

extern "C" int bind_texture_spin(dev_spinor* s, int i){
  
  size_t size;
  
  #ifdef MPI
    if(even_odd_flag){
      size = sizeof(float4)*6*(VOLUME+RAND)/2;
    }
    else{
      size = sizeof(float4)*6*(VOLUME+RAND);
    }
  #else
    if(even_odd_flag){
      size = sizeof(float4)*6*VOLUME/2;
    }
    else{
      size = sizeof(float4)*6*VOLUME;
    }
  #endif
   
  
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
 
 cudaGetTextureReference(&gf_texRefPtr, "gf_tex");
 gf_channelDesc =  cudaCreateChannelDesc<float4>();
 cudaBindTexture(0, gf_texRefPtr, gf, &gf_channelDesc, size);
 //printf("%s\n", cudaGetErrorString(cudaGetLastError()));    
 return(0);
}


extern "C" int unbind_texture_gf(){
 //printf("Unbinding texture to gaugefield\n");
 cudaUnbindTexture(gf_texRefPtr);
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














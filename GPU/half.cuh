
#define twopi_float 6.283185307f

// stores the float spinor field in s into the half spinor field sh and the norm into shnorm
__global__ void float2half_spinorfield(dev_spinor* s, dev_spinor_half* sh, float* shnorm){
  int pos=threadIdx.x + blockDim.x*blockIdx.x;
  __shared__ dev_spinor slocal[6];
  int i;
  float norm = 0.0;
  if(pos < dev_VOLUME){
    dev_copy_spinor(&(s[6*pos]), &(slocal[0]));
    // calculate norm
    #pragma unroll
    for(i=0; i<6; i++){
       norm += slocal[i].x*slocal[i].x + slocal[i].y*slocal[i].y 
             + slocal[i].z*slocal[i].z + slocal[i].w*slocal[i].w; 
     }
    //store norm
    shnorm[pos] = sqrt(norm);
    //store unit direction vector
    #pragma unroll
    for(i=0; i<6; i++){
       sh[6*pos+i].x = __float2half(slocal[i].x/norm);
       sh[6*pos+i].y = __float2half(slocal[i].y/norm);
       sh[6*pos+i].z = __float2half(slocal[i].z/norm);
       sh[6*pos+i].w = __float2half(slocal[i].w/norm);
    } 
  }
}



// stores the float4 gauge field gf into the half gauge field gfh
// for GF_8 we have to be careful, as we have two angles in -Pi .. Pi
// so we have to divide them by (2 Pi) This is taken care of in the gauge
// reconstruction routines
__global__ float2half_gaugefield(dev_su3_2v* gf, dev_su3_2v_half* gfh){

  int pos=threadIdx.x + blockDim.x*blockIdx.x;
  int nf4;
  if(pos < dev_VOLUME){
  #ifdef GF_8
    gfh[nf4*pos].x = __float2half(gf[nf4*pos].x);
    gfh[nf4*pos].y = __float2half(gf[nf4*pos].y);
    gfh[nf4*pos].z = __float2half(gf[nf4*pos].z);
    gfh[nf4*pos].w = __float2half(gf[nf4*pos].w);
    
    gfh[nf4*pos+1].x = __float2half(gf[nf4*pos+1].x/twopi_float);
    gfh[nf4*pos+1].y = __float2half(gf[nf4*pos+1].y/twopi_float);
    gfh[nf4*pos+1].z = __float2half(gf[nf4*pos+1].z);
    gfh[nf4*pos+1].w = __float2half(gf[nf4*pos+1].w);
  #else
    int nf4 = 3;
    #pragma unroll
    for(i=0; i<nf4; i++){
       gfh[nf4*pos+i].x = __float2half(gf[nf4*pos+i].x);
       gfh[nf4*pos+i].y = __float2half(gf[nf4*pos+i].y);
       gfh[nf4*pos+i].z = __float2half(gf[nf4*pos+i].z);
       gfh[nf4*pos+i].w = __float2half(gf[nf4*pos+i].w);
    } 
  #endif
  }
}




extern "C" int prepare_halfspinor_texture(dev_spinor* s, dev_spinor_half* sh, float* shnorm){
  
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
   
   // determine gridsize for conversion to half
   if( VOLUME >= BLOCK2){
      gridsize = (int) (VOLUME/BLOCK2) +1;
    }
    else{
      gridsize=1;
    }
   
   printf("Converting spinor to half precision... ");
     float2half_spinorfield <<< gridsize, BLOCK2  >>>(s, sh, shnorm);
   printf("Done\n");
   
   printf("Binding textures to half spinorfield\n");
    // bind texture for vector
    spin_texRefPtr = NULL;
    cudaGetTextureReference(&spin_texRefPtr, "spin_tex");
    spin_channelDesc =  cudaCreateChannelDesc<short4>();
    cudaBindTexture(0, spin_texRefPtr, s, &spin_channelDesc, size);
  
    // bind texture for norm
    spinnorm_texRefPtr = NULL;
    cudaGetTextureReference(&spinnorm_texRefPtr, "spinnorm_tex");
    spinnorm_channelDesc =  cudaCreateChannelDesc<float>();
    cudaBindTexture(0, spinnorm_texRefPtr, s, &spinnorm_channelDesc, sizenorm);
  //printf("%s\n", cudaGetErrorString(cudaGetLastError())); 
 

return(0);  
}


extern "C" int release_halfspinor_texture(){
   printf("Unbinding textures of half spinorfield\n");
   cudaUnbindTexture(spin_texRefPtr);
   cudaUnbindTexture(spinnorm_texRefPtr);
   //printf("%s\n", cudaGetErrorString(cudaGetLastError()));    
return(0);    
}







/////    SOME BLAS KERNELs - work in progress //////////

__global__ float2half ( dev_halfspinor* erg, float* erg_norm, dev_spinor* in){
  int pos=threadIdx.x + blockDim.x*blockIdx.x;
  if(pos < dev_VOLUME){
     float norm=0.0;
     float4 help;
     //calculate norm over float4 spinor
     for(i=0; i<6; i++){
       help = in[6*pos+i];
       norm += help.x*help.x +  help.y*help.y 
             + help.z*help.z + help.w*help.w; 
     }
     norm = sqrt(norm);
     erg_norm[pos] = norm;
     // normalize and save to erg
     for(i=0; i<6; i++){
       help = in[6*pos+i];
       erg.x[6*pos+i].x = help.x/norm;
       erg.y[6*pos+i].y = help.y/norm;
       erg.z[6*pos+i].z = help.z/norm;
       erg.w[6*pos+i].w = help.w/norm;
     }
  }
}







// erg(float) = alpha*x(half) + y(half) 
__global__ void saxpy_half_kernel (float alpha, dev_spinor* erg, dev_halfspinor* x, float* x_norm, dev_halfspinor* y, float* y_norm){
   int pos= threadIdx.x + blockDim.x*blockIdx.x;
   float4 xhelp,yhelp;
   int i;
   float xnhelp, ynhelp;
   if(pos < dev_VOLUME){
    // this is the loop over the 6 float4 forming one spinor
    #pragma unroll
    for(i=0; i<6; i++){
       xhelp = tex1Dfetch(x, 6*pos+i);
       xnhelp = x_norm[pos];
       xhelp.x = xhelp.x*xnhelp;
       xhelp.y = xhelp.y*xnhelp;
       xhelp.z = xhelp.z*xnhelp;
       xhelp.w = xhelp.w*xnhelp;

       yhelp = tex1Dfetch(y, 6*pos+i);
       ynhelp = y_norm[pos];
       yhelp.x = yhelp.x*ynhelp;
       yhelp.y = yhelp.y*ynhelp;
       yhelp.z = yhelp.z*ynhelp;
       yhelp.w = yhelp.w*ynhelp;   
   
       erg[6*pos+i].x = alpha*xhelp.x + yhelp.x;
       erg[6*pos+i].y = alpha*xhelp.y + yhelp.y;
       erg[6*pos+i].z = alpha*xhelp.z + yhelp.z;
       erg[6*pos+i].w = alpha*xhelp.w + yhelp.w;
    }
   }
}



// helperg = alpha*x + y 
// y = normalize_half (helperg)
void saxpy_half(int N, float alpha, dev_spinor* helperg, dev_halfspinor* x, float* x_norm, dev_halfspinor* y, float* ynorm,  ){
  // N is total number of floats Vol*6*4
  // -> number of ushort4's = N/4
  int Neff = N/4/6;

}


/////////////////////  BLAS KERNELs  ////////////////////////////











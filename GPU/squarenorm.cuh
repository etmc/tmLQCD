


__device__ inline REAL dev_skalarprod_spinor(dev_spinor * s1, dev_spinor * s2){
  REAL skalprod = 0.0;
  int i;
  #pragma unroll 6
  for(i=0;i<6;i++){ //color + spin
    skalprod += ((*(s1+i)).x*(*(s2+i)).x + (*(s1+i)).y*(*(s2+i)).y + (*(s1+i)).z*(*(s2+i)).z + (*(s1+i)).w*(*(s2+i)).w);
  }
  return skalprod;
}




__device__ inline REAL dev_squarenorm_spinor(dev_spinor * s1){
  REAL skalprod = 0.0;
  int i;
  #pragma unroll 6
  for(i=0;i<6;i++){ //color + spin
    skalprod += ((*(s1+i)).x*(*(s1+i)).x + (*(s1+i)).y*(*(s1+i)).y + (*(s1+i)).z*(*(s1+i)).z + (*(s1+i)).w*(*(s1+i)).w);
  }
  return skalprod;
}



__device__ inline REAL dev_squarenorm_spinor_tex(int pos){
  REAL skalprod = 0.0;
  int i;
  float4 help;
  
  #pragma unroll 6
  for(i=0;i<6;i++){ //color + spin
    help = tex1Dfetch(spin_tex2,6*pos+i);
    skalprod += help.x*help.x + help.y*help.y + help.z*help.z + help.w*help.w;
  }
  return skalprod;
}




//only 1 dim parallel possible, because need __syncthread !
__global__ void dev_skalarprod_spinor_field2(dev_spinor* s1, dev_spinor* s2, REAL* erg){
  __shared__ REAL shrinkarray[ACCUM_N];
  int pos,stepwidth;
  REAL ks,kc,ds,tr,ts,tt;
  
   
   // ADD ERROR HERE if t > maxblockdim
   
   ks=0.0;
   kc=0.0; 
   
   if(blockDim.x > dev_VOLUME){
     stepwidth = 1;  
   }
   else{
     stepwidth = dev_VOLUME/(gridDim.x*blockDim.x);
   }
   
     int start = (blockIdx.x*blockDim.x + threadIdx.x)*stepwidth;
     int end = (blockIdx.x*blockDim.x + threadIdx.x+1)*stepwidth;  
   
   for(pos=start;pos<end; pos++){
     if(pos < dev_VOLUME){
          ds = dev_skalarprod_spinor(&(s1[6*pos]), &(s2[6*pos] ));
      }
      else{
        ds = 0.0;
      } 
          // Kahan summation
          tr=ds+kc;
          ts=tr+ks;
          tt=ts-ks;
          ks=ts;
          kc=tr-tt;
    }
    kc=ks+kc;
    shrinkarray[threadIdx.x] = kc;
    __syncthreads();
    
    
    if(threadIdx.x==0){
      ks=0.0;
      kc=0.0; 
      int k;
      for(k=0; k<blockDim.x; k++){
        ds = shrinkarray[k];
        tr=ds+kc;
        ts=tr+ks;
        tt=ts-ks;
        ks=ts;
        kc=tr-tt;
      }
      kc=ks+kc;
      (erg[blockIdx.x])=kc;
    }//threadIdx==0;
}








//only 1 dim parallel possible, because need __syncthread !
__global__ void dev_squarenorm_spinor_field(dev_spinor* s1, REAL* erg){
  __shared__ REAL shrinkarray[ACCUM_N];
  int pos,stepwidth;
  REAL ks,kc,ds,tr,ts,tt;
  
   
   // ADD ERROR HERE if t > maxblockdim
   
   ks=0.0;
   kc=0.0; 
   
   if(blockDim.x > dev_VOLUME){
     stepwidth = 1;  
   }
   else{
     stepwidth = dev_VOLUME/(gridDim.x*blockDim.x);
   }
   
     int start = (blockIdx.x*blockDim.x + threadIdx.x)*stepwidth;
     int end = (blockIdx.x*blockDim.x + threadIdx.x+1)*stepwidth;  
   
   for(pos=start;pos<end; pos++){
     if(pos < dev_VOLUME){
         //ds = dev_squarenorm_spinor_tex(pos);
         ds = dev_squarenorm_spinor(&(s1[6*pos]));
      }
      else{
        ds = 0.0;
      } 
          // Kahan summation
          tr=ds+kc;
          ts=tr+ks;
          tt=ts-ks;
          ks=ts;
          kc=tr-tt;
    }
    kc=ks+kc;
    shrinkarray[threadIdx.x] = kc;
    __syncthreads();
    
    
    if(threadIdx.x==0){
      ks=0.0;
      kc=0.0; 
      int k;
      for(k=0; k<blockDim.x; k++){
        ds = shrinkarray[k];
        tr=ds+kc;
        ts=tr+ks;
        tt=ts-ks;
        ks=ts;
        kc=tr-tt;
      }
      kc=ks+kc;
      (erg[blockIdx.x])=kc;
    }//threadIdx==0;
}






//only 1 dim parallel, because need __syncthread !
__global__ void dev_skalarprod_spinor_field(dev_spinor* s1, dev_spinor* s2, REAL* erg){
  __shared__ REAL shrinkarray[ACCUM_N];
  int pos,stepwidth, sweepsperthread;
  REAL ks,kc,ds,tr,ts,tt;
   
   // ADD ERROR HERE if t > maxblockdim
   
   ks=0.0;
   kc=0.0; 
   
   if(ACCUM_N > dev_VOLUME){
     stepwidth = 1;
     sweepsperthread = 1;  
   }
   else{
     stepwidth = dev_VOLUME/ACCUM_N;
     sweepsperthread = ACCUM_N/blockDim.x;
   }
    
   
   
 for(int j = 0; j < sweepsperthread; j++){
   
     int start = (threadIdx.x + j*blockDim.x)*stepwidth;
     int end = (threadIdx.x+j*blockDim.x+1)*stepwidth;  
     ks=0.0;
     kc=0.0; 
     
   for(pos=start;pos<end; pos++){
     if(pos < dev_VOLUME){
          ds = dev_skalarprod_spinor(&(s1[6*pos]), &(s2[6*pos] ));
          
      }
      else{
        ds = 0.0;
      } 
          // Kahan summation
          tr=ds+kc;
          ts=tr+ks;
          tt=ts-ks;
          ks=ts;
          kc=tr-tt;
    }
    kc=ks+kc;
    shrinkarray[threadIdx.x+j*blockDim.x] = kc;
  }
  __syncthreads();
   
 
    for(int stride = ACCUM_N / 2; stride > 0; stride >>= 1){
       __syncthreads();
       for(int iAccum = threadIdx.x; iAccum < stride; iAccum += blockDim.x)
           shrinkarray[iAccum] += shrinkarray[stride + iAccum];
    }

    if(threadIdx.x == 0) (*erg) = shrinkarray[0];
    
    
    /*
    if(threadIdx.x==0){
      ks=0.0;
      kc=0.0; 
      int k;
      for(k=0; k<sweepsperthread*blockDim.x; k++){
        ds = shrinkarray[k];
        tr=ds+kc;
        ts=tr+ks;
        tt=ts-ks;
        ks=ts;
        kc=tr-tt;
      }
      kc=ks+kc;
      (*erg)=kc;
    }//threadIdx==0;
    */
    
}




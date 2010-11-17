
#define twopi_float 6.283185307f




//////// short <-> float conversion /////////////
#define SHORT_LEN 65536
#define SCALE ((SHORT_LEN-1) * 0.5)
#define SHIFT (-1.0f/(SHORT_LEN-1))

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
              

// stores the float spinor field in s into the half spinor field sh and the norm into shnorm
__global__ void float2half_spinorfield(dev_spinor* s, dev_spinor_half* sh, float* shnorm){
  int pos=threadIdx.x + blockDim.x*blockIdx.x;
  __shared__ float4 slocal[6];
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
     
     for(i=0; i<6; i++){
       norm += s[6*pos+i].x*s[6*pos+i].x 
             + s[6*pos+i].y*s[6*pos+i].y 
             + s[6*pos+i].z*s[6*pos+i].z 
             + s[6*pos+i].w*s[6*pos+i].w; 
     }
   
   
    //store unit direction vector
    if (norm != 0.0f){
      norm = sqrtf(norm);
      shnorm[pos] = norm;
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
       shnorm[pos] = norm;
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
// so we have to divide them by (2 Pi) This is taken care of in the gauge
// reconstruction routines
// the volume is given explicitly (vol) here, to make sure alway the complete 
// VOLUME is processed and not VOLUME/2 (eo)
__global__ void float2half_gaugefield(dev_su3_2v* gf, dev_su3_2v_half* gfh, int vol){

  int pos=threadIdx.x + blockDim.x*blockIdx.x;
  int nf4,i;
  if(pos < vol){
  #ifdef GF_8
    nf4 = 2;
    gfh[nf4*pos].x = fl2sh(gf[nf4*pos].x);
    gfh[nf4*pos].y = fl2sh(gf[nf4*pos].y);
    gfh[nf4*pos].z = fl2sh(gf[nf4*pos].z);
    gfh[nf4*pos].w = fl2sh(gf[nf4*pos].w);
    
    gfh[nf4*pos+1].x = fl2sh(gf[nf4*pos+1].x/twopi_float);
    gfh[nf4*pos+1].y = fl2sh(gf[nf4*pos+1].y/twopi_float);
    gfh[nf4*pos+1].z = fl2sh(gf[nf4*pos+1].z);
    gfh[nf4*pos+1].w = fl2sh(gf[nf4*pos+1].w);
  #else
    nf4 = 3;
    for(i=0; i<nf4; i++){
       gfh[nf4*pos+i].x = fl2sh(gf[nf4*pos+i].x);
       gfh[nf4*pos+i].y = fl2sh(gf[nf4*pos+i].y);
       gfh[nf4*pos+i].z = fl2sh(gf[nf4*pos+i].z);
       gfh[nf4*pos+i].w = fl2sh(gf[nf4*pos+i].w);
    } 
  #endif
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


extern "C" int bind_halfspinor_texture(dev_spinor* s, dev_spinor_half* sh, float* shnorm){
  
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
  
   /*  
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
   */

   //printf("Converting spinor to half precision... ");
     float2half_spinorfield <<< gridsize, BLOCK2  >>>(s, sh, shnorm);
   //printf("Done\n");
    cudaError_t cudaerr;
    if((cudaerr=cudaGetLastError()) != cudaSuccess){
       printf("%s\n", cudaGetErrorString(cudaerr));
       exit(200);
     }

   /*
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







/////    SOME BLAS KERNELs - work in progress //////////





// y(half) = alpha*x(half) + y(half) 
// x is not read from texture
// y is not read from texture
__global__ void axpy_half (float alpha, dev_spinor* erg, dev_spinor_half* x, float* x_norm, dev_spinor_half* y, float* y_norm){
   int pos= threadIdx.x + blockDim.x*blockIdx.x;
   float4 xhelp,yhelp;
   __shared__ float4 erghelp[6];
   int i;
   float xnhelp, ynhelp;
   
   
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
        
       erghelp[6*pos+i].x = alpha*xhelp.x + yhelp.x;
       erghelp[6*pos+i].y = alpha*xhelp.y + yhelp.y;
       erghelp[6*pos+i].z = alpha*xhelp.z + yhelp.z;
       erghelp[6*pos+i].w = alpha*xhelp.w + yhelp.w;
    }
    
	//calculate norm of resulting spinor
	float ergnorm = 0.0f;

	for(i=0; i<6; i++){
                  ergnorm+= erghelp[6*pos+i].x*erghelp[6*pos+i].x
	                  + erghelp[6*pos+i].y*erghelp[6*pos+i].y
			  + erghelp[6*pos+i].z*erghelp[6*pos+i].z
			  + erghelp[6*pos+i].w*erghelp[6*pos+i].w;
	}
	ergnorm = sqrt(ergnorm);
	y_norm[pos] = ergnorm;

	//write out normalized spinors in half
    if (ergnorm != 0.0f){
      #pragma unroll 6
      for(i=0; i<6; i++){
         y[6*pos+i].x = fl2sh(erghelp[6*pos+i].x/ergnorm);
         y[6*pos+i].y = fl2sh(erghelp[6*pos+i].y/ergnorm);
         y[6*pos+i].z = fl2sh(erghelp[6*pos+i].z/ergnorm);
         y[6*pos+i].w = fl2sh(erghelp[6*pos+i].w/ergnorm);
          } 
        }
      else{
       #pragma unroll 6
       for(i=0; i<6; i++){
         y[6*pos+i].x = fl2sh(0.0f);
         y[6*pos+i].y = fl2sh(0.0f);
         y[6*pos+i].z = fl2sh(0.0f);
         y[6*pos+i].w = fl2sh(0.0f);
       }
     }


   }//dev_VOLUME
}




void init_blas_half(int vol){
  cudaError_t cudaerr;
  
  blas_half_blocksize=BLOCK2;
  if( VOLUME >= BLOCK2){
   blas_half_gridsize = (int)(VOLUME/BLOCK2) + 1;
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
__global__ void squarenorm_half (float* redfield, float* x_norm){
   int pos= threadIdx.x + blockDim.x*blockIdx.x;
   int i;
   float xnhelp;
   if(pos < dev_VOLUME){
    xnhelp = x_norm[pos];
    redfield[pos] = xnhelp*xnhelp;
   }//dev_VOLUME
}







// calculates the dot product of x and y
float dotprod_half(dev_spinor_half* x, float* x_norm, dev_spinor_half* y, float* y_norm){
   int i;
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
   return(finalsum);
}



// calculates the norm^2 of a half spinor with norm x 
float squarenorm_half(float * xnorm){
   int i;
   cudaError_t cudaerr;
   
   squarenorm_half <<< blas_half_gridsize, blas_half_blocksize >>> 
                      (dev_blas_half_redfield, xnorm);
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
   
   return(finalsum);
}






/////////////////////  BLAS KERNELs  ////////////////////////////








//////////////////   TEST ROUTINES /////////////////////////////


__global__ void testhalf (dev_su3_2v* gf, dev_su3 * to){

  int pos=threadIdx.x + blockDim.x*blockIdx.x;
  if (pos==0){
    dev_reconstructgf_2vtexref(gf, 2 , to);
  }
}


void testhalf_gf(dev_su3_2v * gf){
  dev_su3 * testfield;
  dev_su3 hosttestfield;
  size_t size =  sizeof(dev_su3);
  dev_su3 hostmatrix;

  
  cudaMalloc((void **) &testfield, size);
  #ifdef USETEXTURE
    //Bind texture gf
    bind_texture_gf(gf);
  #endif  
  testhalf <<< 1, 1 >>> (dev_gf, testfield);
  cudaMemcpy(&(hosttestfield), testfield, size, cudaMemcpyDeviceToHost);  
  show_dev_su3(hosttestfield);
  
  
  #ifdef USETEXTURE
   unbind_texture_gf();
  #endif
}








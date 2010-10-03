





/// This uses all the single precision device gauge field
/// Think of porting this to double !


float*  redfield_plaq;
float * dev_redfield_plaq;

float*  redfield_rect;
float * dev_redfield_rect;




int init_dev_observables(){
  cudaError_t cudaerr;
  
  cudaMalloc((void **) &dev_redfield_plaq, T*sizeof(float));
  if((redfield_plaq = (float*)malloc(T*sizeof(float)))==(void*)NULL){
    fprintf(stderr,"Error in init_dev_observables: malloc error(plaq)\n");
    return(1);
  }

  cudaMalloc((void **) &dev_redfield_rect, T*sizeof(float));
  if((redfield_rect = (float*)malloc(T*sizeof(float)))==(void*)NULL){
    fprintf(stderr,"Error in init_dev_observables: malloc error(rect)\n");
    return(1);
  }

  if((cudaerr=cudaGetLastError())!=cudaSuccess){
    fprintf(stderr, "Error in init_dev_observables(): GPU memory allocation of reduction fields failed. Aborting...\n");
    return(2);
  }
  
 return(0);
}





void finalize_dev_observables(){

  free(redfield_plaq);
  free(redfield_rect);
  cudaFree(dev_redfield_plaq);
  cudaFree(dev_redfield_rect);
  
}




// update the global device gaugefield with host gaugefield given by gf
// this is all single precision!
void update_dev_gaugefield(su3** gf){

  size_t dev_gfsize;
  #ifdef GF_8
    dev_gfsize = 2*4*VOLUME * sizeof(dev_su3_8);
    su3to8(gf,h2d_gf);  
  #else
    dev_gfsize = 3*4*VOLUME * sizeof(dev_su3_2v);
    su3to2vf4(gf,h2d_gf);
  #endif
  cudaMemcpy(dev_gf, h2d_gf, dev_gfsize, cudaMemcpyHostToDevice);
}












// calculates the mean plaquette of the gauge field
// uses 2d parallelization:
// N_grid = T, N_block = LZ
// CPU has to do last summation over T
// BLOCKPLAQ MUST be larger than LZ !!! -> define correctly in cudaglobal.h
__global__ void dev_mean_plaq(float* reductionfield_plaq, int * dev_nn, dev_su3_2v * gf){
  float mplaq = 0.0;
  int x0pos, x1pos, x2pos ; /* x0pos = basepoint of plaquette, x1pos = x0pos + e_mu, x2pos = x0pos + e_nu */
  int iz,x,y,z,t,mu,nu;
  dev_su3 su3matrix,su3matrix2, M1,M2,M3,M4;
  dev_complex chelp;
  

  __shared__ float output[BLOCKPLAQ];
  t = blockIdx.x;
  z = threadIdx.x;
  

      for(y=0; y<dev_LY; y++){
        for(x=0; x<dev_LX; x++){
          for(nu=0;nu <3; nu++){
            for(mu =nu+1; mu < 4; mu++){
              x0pos = x + dev_LX*(y + dev_LY*(z + dev_LZ*t));
              x1pos = dev_nn[8*x0pos + mu];
              x2pos = dev_nn[8*x0pos + nu];          

/* U_nu(x) */
            #ifdef GF_8
              dev_reconstructgf_8texref(gf, (4*x0pos+mu),&M1);
            #else
              dev_reconstructgf_2vtexref(gf, (4*x0pos+mu),&M1);
            #endif
/* U_mu(x+e_mu) */
            #ifdef GF_8
              dev_reconstructgf_8texref(gf, (4*x1pos+nu),&M2);
            #else
              dev_reconstructgf_2vtexref(gf, (4*x1pos+nu),&M2);
            #endif

/* Udagger_nu(x+e_nu) */
            #ifdef GF_8
              dev_reconstructgf_8texref_dagger(gf, (4*x2pos+mu),&M3);
            #else
              dev_reconstructgf_2vtexref_dagger(gf, (4*x2pos+mu),&M3);
            #endif
/* Udagger_mu(x)*/
            #ifdef GF_8
              dev_reconstructgf_8texref_dagger(gf, (4*x0pos+nu),&M4);
            #else
              dev_reconstructgf_2vtexref_dagger(gf, (4*x0pos+nu),&M4);
            #endif
  
              /* multiply these and store in su3matrix*/
              dev_su3_ti_su3(&su3matrix, &M3,&M4);
              dev_su3_ti_su3(&su3matrix2, &M2,&su3matrix);
              dev_su3_ti_su3(&su3matrix, &M1,&su3matrix2);

               mplaq += dev_su3Retrace(&su3matrix)/3.0;
            }
          }
           
        }
      } 
    output[z] = mplaq;

  __syncthreads();
  
  if(threadIdx.x == 0){
    
    /* normieren */
    float accum = 0.0;
    for(iz=0; iz < dev_LZ; iz++){
      accum += output[iz];  
    }
    accum = accum*(1.0/(6.0*dev_VOLUME));
    reductionfield_plaq[t] = accum;
  }
  __syncthreads();
  
}













float calc_plaquette(dev_su3_2v * U){
   float erg=0.0;
   int j;
   dev_mean_plaq <<< T , LZ >>> (dev_redfield_plaq, dev_nn, U) ;
   cudaMemcpy(redfield_plaq, dev_redfield_plaq, (size_t)(T*sizeof(float)), cudaMemcpyDeviceToHost);
   // we have to add up the final sum on host
   for(j=0; j<T; j++){
     erg+=redfield_plaq[j];
     //printf("%e\n", redfield_plaq[j]);
   }
   printf("PLAQ = %.8f\n",erg);
   printf("%s\n", cudaGetErrorString(cudaGetLastError()));
   return(erg);
}





//calculates 1 rectangle situated at x,y,z,t in direction mu, nu
//
//          -> nu
//        ^
//        |
//        
//        ^ mu
//        |
//        
__device__ float dev_onerect(int * dev_nn, dev_su3_2v * gf,
 				int t, int x, int y, int z, int mu, int nu ){
              int x0pos, x1pos, x2pos;
              dev_su3 su3matrix,su3matrix2, M1,M2;
             
             
              x0pos = x + dev_LX*(y + dev_LY*(z + dev_LZ*t));
              x1pos = dev_nn[8*x0pos + mu];
              x2pos = dev_nn[8*x1pos + nu];          

/* U_muu(x) */
            #ifdef GF_8
              dev_reconstructgf_8texref(gf, (4*x0pos+mu),&M1);
            #else
              dev_reconstructgf_2vtexref(gf, (4*x0pos+mu),&M1);
            #endif
/* U_mu(x+e_mu) */
            #ifdef GF_8
              dev_reconstructgf_8texref(gf, (4*x1pos+mu),&M2);
            #else
              dev_reconstructgf_2vtexref(gf, (4*x1pos+mu),&M2);
            #endif
            dev_su3_ti_su3(&su3matrix2, &M1,&M2);
            
/* U_nu(x+ 2 e_mu) */
            #ifdef GF_8
              dev_reconstructgf_8texref(gf, (4*x2pos+nu),&M1);
            #else
              dev_reconstructgf_2vtexref(gf, (4*x2pos+nu),&M1);
            #endif
            dev_su3_ti_su3(&su3matrix, &M1, &su3matrix2);



       //x0pos = x + dev_LX*(y + dev_LY*(z + dev_LZ*t));
         x1pos = dev_nn[8*x0pos + nu];
         x2pos = dev_nn[8*x1pos + mu];




/* Udagger_mu(x+e_nu+e_mu)*/
            #ifdef GF_8
              dev_reconstructgf_8texref_dagger(gf, (4*x2pos+mu),&M1);
            #else
              dev_reconstructgf_2vtexref_dagger(gf, (4*x2pos+mu),&M1);
            #endif
            dev_su3_ti_su3(&su3matrix2, &su3matrix , &M1);
 /* Udagger_mu(x+e_nu) */
            #ifdef GF_8
              dev_reconstructgf_8texref_dagger(gf, (4*x1pos+mu),&M2);
            #else
              dev_reconstructgf_2vtexref_dagger(gf, (4*x1pos+mu),&M2);
            #endif 
            dev_su3_ti_su3(&su3matrix, &su3matrix2 , &M2);
/* Udagger_nu(x)*/
            #ifdef GF_8
              dev_reconstructgf_8texref_dagger(gf, (4*x0pos+nu),&M1);
            #else
              dev_reconstructgf_2vtexref_dagger(gf, (4*x0pos+nu),&M1);
            #endif  
            dev_su3_ti_su3(&su3matrix2, &su3matrix , &M1);
  

              
            float chelp = dev_su3Retrace(&su3matrix2);
return(chelp);  
}




// calculates the rectangles of the gauge field
// uses 2d parallelization:
// N_grid = T, N_block = LZ
// CPU has to do last summation over T
// BLOCKPLAQ MUST be larger than LZ !!! -> define correctly in cudaglobal.h
__global__ void dev_rectangle(float* reductionfield_rect, 
            int * dev_nn, dev_su3_2v * gf){
  float mrect = 0.0;
  int iz,x,y,z,t,mu,nu;
  dev_complex chelp;
  
  //space for the sum of blocks
  __shared__ float output[BLOCKPLAQ];
  
  t = blockIdx.x;
  z = threadIdx.x;
  

      for(y=0; y<dev_LY; y++){
        for(x=0; x<dev_LX; x++){
          
          mrect += dev_onerect(dev_nn, gf, t, x, y, z, 0, 1);
          mrect += dev_onerect(dev_nn, gf, t, x, y, z, 0, 2);
          mrect += dev_onerect(dev_nn, gf, t, x, y, z, 0, 3);
          
          mrect += dev_onerect(dev_nn, gf, t, x, y, z, 1, 0);
          mrect += dev_onerect(dev_nn, gf, t, x, y, z, 1, 2);
          mrect += dev_onerect(dev_nn, gf, t, x, y, z, 1, 3);
          
          mrect += dev_onerect(dev_nn, gf, t, x, y, z, 2, 0);
          mrect += dev_onerect(dev_nn, gf, t, x, y, z, 2, 1);
          mrect += dev_onerect(dev_nn, gf, t, x, y, z, 2, 3);
          
          mrect += dev_onerect(dev_nn, gf, t, x, y, z, 3, 0);
          mrect += dev_onerect(dev_nn, gf, t, x, y, z, 3, 1);
          mrect += dev_onerect(dev_nn, gf, t, x, y, z, 3, 2);

           
        }
      } 
    output[z] = mrect;

  __syncthreads();
  
  if(threadIdx.x == 0){
    
    /* normieren */
    float accum = 0.0;
    for(iz=0; iz < dev_LZ; iz++){
      accum += output[iz];  
    }
    accum = accum*(1.0/(6.0*dev_VOLUME));
    reductionfield_rect[t] = accum;
  }
  __syncthreads();
  
}







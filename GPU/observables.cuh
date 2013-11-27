





/// This uses all the single precision device gauge field
/// Think of porting this to double !


// reduction field on host and device
float*  redfield;
float * dev_redfield;
float * dev_sredfield;
float * dev_ssredfield;

// same for polyakov loop
float2 * dev_polyredfield;
float2 * dev_polysredfield;
float2 * dev_polyssredfield;
float2 * polyredfield;

// size of the first (small) and second (smallsmall) reduction fields 
int sredsize, ssredsize;

// same for polyakov loop
int polysredsize, polyssredsize;


int init_dev_observables(){
  cudaError_t cudaerr;
  
  
  // IMPLEMENT THIS FOR ALL LATTICE SIZES !!!!!!!!!!!!!!!!!!!!
  if((VOLUME%REDUCTION_N) == 0){
    sredsize = VOLUME/REDUCTION_N;
  }
  else{
    fprintf(stderr,"Error: Volume is not a multiple of REDUCTION_N (%d). Aborting...\n", REDUCTION_N);
    exit(100);
  }
  
  if(sredsize < REDUCTION_N){
    ssredsize = 1;
  }
  else{
    if(sredsize%REDUCTION_N == 0){
      ssredsize = sredsize/REDUCTION_N;
    }
    else{
      ssredsize = sredsize/REDUCTION_N + 1;
    }
  }
  
  
  //VOLUME * float on device
  cudaMalloc((void **) &dev_redfield, VOLUME*sizeof(float));

  if((redfield = (float*)malloc(sredsize*sizeof(float)))==(void*)NULL){
    fprintf(stderr,"Error in init_dev_observables: malloc error(plaq)\n");
    return(1);
  }
  cudaMalloc((void **) &dev_sredfield, sredsize*sizeof(float));
  cudaMalloc((void **) &dev_ssredfield, ssredsize*sizeof(float));

  if((cudaerr=cudaGetLastError())!=cudaSuccess){
    if(g_proc_id==0) {
      fprintf(stderr, "Error in init_dev_observables(): GPU memory allocation of reduction fields failed. Aborting...\n");
      fprintf(stderr, "Error code is: %f\n",cudaerr);
    }
    return(2);
  }
  


  int spatialvol = LX*LY*LZ;
  if((spatialvol%REDUCTION_N) == 0){
    polysredsize = spatialvol/REDUCTION_N;
  }
  else{
    fprintf(stderr,"Error: spatial Volume is not a multiple of REDUCTION_N (%d). Aborting...\n", REDUCTION_N);
    exit(100);
  }
  
  if(polysredsize < REDUCTION_N){
    polyssredsize = 1;
  }
  else{
    if(polysredsize%REDUCTION_N == 0){
      polyssredsize = polysredsize/REDUCTION_N;
    }
    else{
      polyssredsize = polysredsize/REDUCTION_N + 1;
    }
  }
  

  // spatial volume*2 (->complex) field for Polyakov loop data
  cudaMalloc((void **) &dev_polyredfield, spatialvol*sizeof(float2));
  if((polyredfield = (float2*)malloc(spatialvol*sizeof(float2)))==(void*)NULL){
    fprintf(stderr,"Error in init_dev_observables: malloc error(poly)\n");
    return(1);
  }
  cudaMalloc((void **) &dev_polysredfield, polysredsize*sizeof(float2));//complex !!
  cudaMalloc((void **) &dev_polyssredfield, polyssredsize*sizeof(float2));//complex !!

  if((cudaerr=cudaGetLastError())!=cudaSuccess){
    fprintf(stderr, "Error in init_dev_observables(): GPU memory allocation of poly reduction fields failed. Aborting...\n");
    return(2);
  }

  
  cudaMemcpyToSymbol("dev_VOLUME", &VOLUME, sizeof(int)) ; 
  cudaMemcpyToSymbol("dev_LX", &LX, sizeof(int)) ;
  cudaMemcpyToSymbol("dev_LY", &LY, sizeof(int)) ;
  cudaMemcpyToSymbol("dev_LZ", &LZ, sizeof(int)) ;
  cudaMemcpyToSymbol("dev_T", &T, sizeof(int)) ;


 return(0);
}





void finalize_dev_observables(){
  free(redfield);
  cudaFree(dev_redfield);
  cudaFree(dev_sredfield);
  cudaFree(dev_ssredfield);
  cudaFree(dev_polyredfield);
  cudaFree(dev_polysredfield);
  cudaFree(dev_polyssredfield); 
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
__global__ void dev_mean_plaq(float* reductionfield, int * dev_nn, dev_su3_2v * gf){
  float mplaq = 0.0;
  int x0pos, x1pos, x2pos ; /* x0pos = basepoint of plaquette, x1pos = x0pos + e_mu, x2pos = x0pos + e_nu */
  int mu,nu;
  dev_su3 su3matrix,su3matrix2, M1,M2,M3,M4;
  
 x0pos = threadIdx.x + blockDim.x*blockIdx.x;  
 #ifdef TEMPORALGAUGE
  int spatialvol = dev_LX*dev_LY*dev_LZ;
  int t;
 #endif

 int gaugevol = dev_VOLUME;
 
 
 if(x0pos < dev_VOLUME){
    
    
    //nu == T-direction -> beware temporal gauge and GF8
      nu = 0;
      for(mu =nu+1; mu < 4; mu++){
          x1pos = dev_nn[8*x0pos + mu];
          x2pos = dev_nn[8*x0pos + nu];          

/* U_mu(x) */
            #ifdef GF_8
              dev_reconstructgf_8texref(gf, x0pos, mu, gaugevol ,&M1);
            #else
              dev_reconstructgf_2vtexref(gf, x0pos, mu, gaugevol ,&M1);
            #endif
/* U_nu(x+e_mu) */
     #ifdef TEMPORALGAUGE
       t = x0pos/spatialvol; // are we on timeslice T-1? no -> U==ID
       if(t != (dev_T-1) ){ 
         dev_unit_su3(&M2);
       }
       else{
         #ifdef GF_8
           dev_reconstructgf_8texref(gf, x1pos, nu, gaugevol ,&M2);
         #else
           dev_reconstructgf_2vtexref(gf, x1pos, nu, gaugevol ,&M2);
         #endif
       }
     #else
            #ifdef GF_8
              dev_reconstructgf_8texref(gf, x1pos, nu, gaugevol ,&M2);
            #else
              dev_reconstructgf_2vtexref(gf, x1pos, nu, gaugevol ,&M2);
            #endif
     #endif /*TEMPORALGAUGE*/
     dev_su3_ti_su3(&su3matrix, &M1,&M2);
     
/* Udagger_mu(x+e_nu) */
            #ifdef GF_8
              dev_reconstructgf_8texref_dagger(gf, x2pos, mu, gaugevol ,&M3);
            #else
              dev_reconstructgf_2vtexref_dagger(gf, x2pos, mu, gaugevol ,&M3);
            #endif
    dev_su3_ti_su3(&su3matrix2, &su3matrix,&M3);
/* Udagger_nu(x)*/
     #ifdef TEMPORALGAUGE
       if(t != (dev_T-1) ){ 
          dev_unit_su3(&M4);
       }
       else{
         #ifdef GF_8
           dev_reconstructgf_8texref_dagger(gf, x0pos, nu, gaugevol ,&M4);
         #else
           dev_reconstructgf_2vtexref_dagger(gf, x0pos, nu, gaugevol,&M4);
         #endif
       }
     #else
       #ifdef GF_8
         dev_reconstructgf_8texref_dagger(gf, x0pos, nu, gaugevol ,&M4);
       #else
         dev_reconstructgf_2vtexref_dagger(gf, x0pos, nu, gaugevol,&M4);
       #endif
     #endif  /*TEMPORALGAUGE*/
     dev_su3_ti_su3(&su3matrix, &su3matrix2,&M4);
     mplaq += dev_su3Retrace(&su3matrix)/3.0;
     
     /* multiply these and store in su3matrix*/
     /*
     dev_su3_ti_su3(&su3matrix, &M3,&M4);
              dev_su3_ti_su3(&su3matrix2, &M2,&su3matrix);
              dev_su3_ti_su3(&su3matrix, &M1,&su3matrix2);

               mplaq += dev_su3Retrace(&su3matrix)/3.0;
      */
      }
                   
          
  // nu != T-direction -> no problem with temporal gauge and GF8
          for(nu=1;nu <3; nu++){
            for(mu =nu+1; mu < 4; mu++){
              x1pos = dev_nn[8*x0pos + mu];
              x2pos = dev_nn[8*x0pos + nu];          

/* U_nu(x) */
            #ifdef GF_8
              dev_reconstructgf_8texref(gf, x0pos, mu, gaugevol ,&M1);
            #else
              dev_reconstructgf_2vtexref(gf, x0pos, mu, gaugevol ,&M1);
            #endif
/* U_mu(x+e_mu) */
            #ifdef GF_8
              dev_reconstructgf_8texref(gf, x1pos, nu, gaugevol ,&M2);
            #else
              dev_reconstructgf_2vtexref(gf, x1pos, nu, gaugevol ,&M2);
            #endif
            dev_su3_ti_su3(&su3matrix, &M1,&M2);
/* Udagger_nu(x+e_nu) */
            #ifdef GF_8
              dev_reconstructgf_8texref_dagger(gf, x2pos, mu, gaugevol ,&M3);
            #else
              dev_reconstructgf_2vtexref_dagger(gf, x2pos, mu, gaugevol ,&M3);
            #endif
            dev_su3_ti_su3(&su3matrix2, &su3matrix,&M3);
/* Udagger_mu(x)*/
            #ifdef GF_8
              dev_reconstructgf_8texref_dagger(gf, x0pos, nu, gaugevol ,&M4);
            #else
              dev_reconstructgf_2vtexref_dagger(gf, x0pos, nu, gaugevol ,&M4);
            #endif
            dev_su3_ti_su3(&su3matrix, &su3matrix2,&M4);
            mplaq += dev_su3Retrace(&su3matrix)/3.0;
              /* multiply these and store in su3matrix*/
              /*
              dev_su3_ti_su3(&su3matrix, &M3,&M4);
              dev_su3_ti_su3(&su3matrix2, &M2,&su3matrix);
              dev_su3_ti_su3(&su3matrix, &M1,&su3matrix2);

              mplaq += dev_su3Retrace(&su3matrix)/3.0;
              */
            }
          }

    reductionfield[x0pos] = mplaq;
  }
}







float calc_plaquette(dev_su3_2v * U, int* nn){
  float erg=0.0;
  int j;
  #ifdef USETEXTURE
   //Bind texture gf
   bind_texture_gf(U);
  #endif 
  
  int gridsize;
  int blocksize=BLOCK2;
  if( VOLUME >= BLOCK2){
   gridsize = (int)(VOLUME/BLOCK2) + 1;
  }
  else{
    gridsize=1;
  }
    
  dev_mean_plaq <<< gridsize , blocksize >>> (dev_redfield, nn, U) ;
  printf("Plaquette calculation on device: %s\n", cudaGetErrorString(cudaGetLastError()));
  
  #ifdef USETEXTURE
    unbind_texture_gf();
  #endif
  
  int redblocks;
  if(sredsize > 1){
    redblocks = sredsize; // VOLUME/REDUCTION_N
    reduce_float <<< redblocks, REDUCTION_N, REDUCTION_N*sizeof(float) >>> 
        ( dev_redfield, dev_sredfield,  VOLUME);
   printf("Reduction 1 of data: %s\n", cudaGetErrorString(cudaGetLastError()));
  }
  if(ssredsize > 1){
    redblocks = ssredsize; 
    reduce_float <<< redblocks, REDUCTION_N, REDUCTION_N*sizeof(float) >>> 
        ( dev_sredfield, dev_ssredfield,  sredsize );
    printf("Reduction 2 of data: %s\n", cudaGetErrorString(cudaGetLastError()));
  
    cudaMemcpy(redfield, dev_ssredfield, (size_t)(redblocks*sizeof(float)), cudaMemcpyDeviceToHost);
  }
  else{
    cudaMemcpy(redfield, dev_sredfield, (size_t)(redblocks*sizeof(float)), cudaMemcpyDeviceToHost);
  }
  
  // we have to add up the final sum on host
  for(j=0; j<redblocks; j++){
    erg+=redfield[j];
    //printf("%e\n", redfield[j]);
  }

  
/*
  cudaMemcpy(redfield, dev_sredfield, (size_t)(redblocks*sizeof(float)), cudaMemcpyDeviceToHost);
  
  // we have to add up the final sum on host
  for(j=0; j<redblocks; j++){
    erg+=redfield[j];
    //printf("%e\n", redfield[j]);
  }
*/
  
  erg=erg/6.0/VOLUME;
  
  
  //printf("PLAQ = %.8f\n",erg);
  //printf("%s\n", cudaGetErrorString(cudaGetLastError()));
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
//        x
__device__ float dev_onerect(int * dev_nn, dev_su3_2v * gf, int tid,
 				int x0pos, int mu, int nu ){
            int x1pos, x2pos;
            dev_su3  M1,M2, su3matrix;
            int gaugevol = dev_VOLUME; 

            __shared__ dev_su3 su3matrix2[BLOCK]; 
            
            

              x1pos = dev_nn[8*x0pos + mu];
              x2pos = dev_nn[8*x1pos + mu];          

/* U_mu(x) */
            #ifdef GF_8
              dev_reconstructgf_8texref(gf, x0pos, mu, gaugevol ,&M1);
            #else
              dev_reconstructgf_2vtexref(gf, x0pos, mu, gaugevol ,&M1);
            #endif
/* U_mu(x+e_mu) */
            #ifdef GF_8
              dev_reconstructgf_8texref(gf, x1pos, mu, gaugevol ,&M2);
            #else
              dev_reconstructgf_2vtexref(gf, x1pos, mu, gaugevol,&M2);
            #endif
            dev_su3_ti_su3(&su3matrix2[tid], &M1,&M2);
            
/* U_nu(x+ 2 e_mu) */
            #ifdef GF_8
              dev_reconstructgf_8texref(gf, x2pos, nu, gaugevol ,&M1);
            #else
              dev_reconstructgf_2vtexref(gf, x2pos, nu, gaugevol ,&M1);
            #endif
            dev_su3_ti_su3(&su3matrix, &su3matrix2[tid], &M1);



       //x0pos = x + dev_LX*(y + dev_LY*(z + dev_LZ*t));
         x1pos = dev_nn[8*x0pos + mu];
         x2pos = dev_nn[8*x1pos + nu];
         
         x1pos = dev_nn[8*x0pos + nu];

/* Udagger_mu(x+e_nu+e_mu)*/
            #ifdef GF_8
              dev_reconstructgf_8texref_dagger(gf, x2pos, mu, gaugevol ,&M1);
            #else
              dev_reconstructgf_2vtexref_dagger(gf, x2pos, mu, gaugevol ,&M1);
            #endif
            dev_su3_ti_su3(&su3matrix2[tid], &su3matrix , &M1);
 /* Udagger_mu(x+e_nu) */
            #ifdef GF_8
              dev_reconstructgf_8texref_dagger(gf, x1pos, mu, gaugevol ,&M2);
            #else
              dev_reconstructgf_2vtexref_dagger(gf, x1pos, mu, gaugevol ,&M2);
            #endif 
            dev_su3_ti_su3(&su3matrix, &su3matrix2[tid] , &M2);
/* Udagger_nu(x)*/
            #ifdef GF_8
              dev_reconstructgf_8texref_dagger(gf, x0pos, nu, gaugevol ,&M1);
            #else
              dev_reconstructgf_2vtexref_dagger(gf, x0pos, nu, gaugevol ,&M1);
            #endif  
            dev_su3_ti_su3(&su3matrix2[tid], &su3matrix , &M1);
  

              
            float help = dev_su3Retrace(&su3matrix2[tid])/3.0;
return(help);  
}




// calculates the rectangles of the gauge field
// uses 2d parallelization:
__global__ void dev_rectangle(float* reductionfield, 
            int * dev_nn, dev_su3_2v * gf){
  float mrect = 0.0;
  int x0pos,ix;

  ix = threadIdx.x;
  x0pos = threadIdx.x + blockDim.x*blockIdx.x;  

 if(x0pos < dev_VOLUME){
  
          
          mrect += dev_onerect(dev_nn, gf, ix, x0pos, 0, 1);
          mrect += dev_onerect(dev_nn, gf, ix, x0pos, 0, 2);
          mrect += dev_onerect(dev_nn, gf, ix, x0pos, 0, 3);
          
          mrect += dev_onerect(dev_nn, gf, ix, x0pos, 1, 0);
          mrect += dev_onerect(dev_nn, gf, ix, x0pos, 1, 2);
          mrect += dev_onerect(dev_nn, gf, ix, x0pos, 1, 3);
          
          mrect += dev_onerect(dev_nn, gf, ix, x0pos, 2, 0);
          mrect += dev_onerect(dev_nn, gf, ix, x0pos, 2, 1);
          mrect += dev_onerect(dev_nn, gf, ix, x0pos, 2, 3);
          
          mrect += dev_onerect(dev_nn, gf, ix, x0pos, 3, 0);
          mrect += dev_onerect(dev_nn, gf, ix, x0pos, 3, 1);
          mrect += dev_onerect(dev_nn, gf, ix, x0pos, 3, 2);


    
    reductionfield[x0pos] = mrect;
  
  }
  
}





float calc_rectangle(dev_su3_2v * U, int* nn){
  float erg=0.0;
  int j;
  #ifdef USETEXTURE
   //Bind texture gf
   bind_texture_gf(U);
  #endif 
  
  int gridsize;
  int blocksize=BLOCK;
  if( VOLUME >= BLOCK){
   gridsize = (int)(VOLUME/BLOCK) + 1;
  }
  else{
    gridsize=1;
  }
    
  dev_rectangle <<< gridsize , blocksize >>> (dev_redfield, nn, U) ;
  printf("Rectangle calculation on device: %s\n", cudaGetErrorString(cudaGetLastError()));
  
  #ifdef USETEXTURE
    unbind_texture_gf();
  #endif
  
  int redblocks;
  if(sredsize > 0){
    redblocks = sredsize; // VOLUME/REDUCTION_N
  
    reduce_float <<< redblocks, REDUCTION_N, REDUCTION_N*sizeof(float) >>> 
        ( dev_redfield, dev_sredfield,  VOLUME);
    printf("Reduction 1 of data: %s\n", cudaGetErrorString(cudaGetLastError()));
  }
  if(ssredsize > 0){
    redblocks = ssredsize; 
    reduce_float <<< redblocks, REDUCTION_N, REDUCTION_N*sizeof(float) >>> 
        ( dev_sredfield, dev_ssredfield,  sredsize );
    printf("Reduction 2 of data: %s\n", cudaGetErrorString(cudaGetLastError()));
  
    cudaMemcpy(redfield, dev_ssredfield, (size_t)(redblocks*sizeof(float)), cudaMemcpyDeviceToHost); 
  }
  else{
    cudaMemcpy(redfield, dev_sredfield, (size_t)(redblocks*sizeof(float)), cudaMemcpyDeviceToHost);
  }
  // we have to add up the final sum on host
  for(j=0; j<redblocks; j++){
    erg+=redfield[j];
    //printf("%e\n", redfield[j]);
  }

  erg=erg/12.0/VOLUME;
  
  //printf("RECT = %.8f\n",erg);
  //printf("%s\n", cudaGetErrorString(cudaGetLastError()));
  return(erg);
}










// calculates the polyakov loop in time direction
__global__ void dev_polyakov_0(float2* reductionfield, int * dev_nn, dev_su3_2v * gf){
  float2 poly;
  poly.x = 0.0; poly.y = 0.0;
  int gaugevol = dev_VOLUME;

  int t, spatialpos;
  dev_su3 M1,  gather, tmp;
  
 spatialpos = threadIdx.x + blockDim.x*blockIdx.x;  
 int spatialvol = dev_LX*dev_LY*dev_LZ;

 if(spatialpos < spatialvol){
   //initialize the polyakov loop su3 matrix 
   dev_unit_su3(&gather);
   
   // t is slowest index, such that the starting gauge field at x,y,z,t=0
   // is actually located at spatialpos
   // -> set actualpos to spatialpos
   int actualpos = spatialpos; 
   
    for(t=0; t < dev_T; t++){
                 
/* U_0(x) */
            #ifdef GF_8
              dev_reconstructgf_8texref(gf, actualpos, 0, gaugevol ,&M1);
            #else
              dev_reconstructgf_2vtexref(gf, actualpos, 0, gaugevol ,&M1);
            #endif
          
        //multiply 
        dev_su3_ti_su3(&tmp, &gather, &M1);
        //store again gather
        dev_su3_assign(&gather, &tmp);

        //go one step in 0-direction
        actualpos = dev_nn[8*actualpos];  
    }
    dev_su3trace(&poly, &gather);
    reductionfield[spatialpos].x = poly.x/3.0;
    reductionfield[spatialpos].y = poly.y/3.0;
  }
}









void calc_polyakov_0(float2* ret, dev_su3_2v * U, int* nn){
  int j;
  float2 erg;
  #ifdef USETEXTURE
   //Bind texture gf
   bind_texture_gf(U);
  #endif 
  
  int gridsize;
  int blocksize=BLOCK2;
  int spatialvol = LX*LY*LZ;
  if( spatialvol >= BLOCK2){
   gridsize = (int)(spatialvol/BLOCK2) + 1;
  }
  else{
    gridsize=1;
  }
    
  dev_polyakov_0 <<< gridsize , blocksize >>> (dev_polyredfield, nn, U) ;
  printf("Polyakov loop calculation on device: %s\n", cudaGetErrorString(cudaGetLastError()));
  
  #ifdef USETEXTURE
    unbind_texture_gf();
  #endif
  
  
  
  
  int redblocks = polysredsize; // VOLUME/REDUCTION_N
  cudaMemcpy(polyredfield, dev_polyredfield, (size_t)(redblocks*sizeof(float2)), cudaMemcpyDeviceToHost);
  /* write to file */
  
  if(polysredsize > 1){
    reduce_float2 <<< redblocks, REDUCTION_N, REDUCTION_N*sizeof(float2) >>> 
        ( dev_polyredfield, dev_polysredfield,  spatialvol);
    printf("Reduction 1 of data: %s\n", cudaGetErrorString(cudaGetLastError()));
  }
  if(polyssredsize > 1){
    redblocks = polyssredsize; 
    reduce_float2 <<< redblocks, REDUCTION_N, REDUCTION_N*sizeof(float2) >>> 
        ( dev_polysredfield, dev_polyssredfield,  polysredsize );
    printf("Reduction 2 of data: %s\n", cudaGetErrorString(cudaGetLastError()));
  
    cudaMemcpy(polyredfield, dev_polyssredfield, (size_t)(redblocks*sizeof(float2)), cudaMemcpyDeviceToHost);
  }
  else{
     cudaMemcpy(polyredfield, dev_polysredfield, (size_t)(redblocks*sizeof(float2)), cudaMemcpyDeviceToHost);
  }  
    
    
  // we have to add up the final sum on host
  for(j=0; j<redblocks; j++){
    erg.x +=polyredfield[j].x;
    erg.y +=polyredfield[j].y;
    //printf("%e\n", polyredfield[j].x);
  }

 
  
  erg.x=erg.x/spatialvol;
  erg.y=erg.y/spatialvol;
  
  
  //printf("pl_0 (Re) = %.8f\n",erg.x);
  //printf("pl_0 (Im) = %.8f\n",erg.y);
  
  (*ret).x = erg.x;
  (*ret).y = erg.y;
  printf("%s\n", cudaGetErrorString(cudaGetLastError()));
 
}










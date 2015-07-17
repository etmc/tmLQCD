

void test_spinor_normalization(dev_spinor_half* s, float* sn){
  dev_spinor_half s_host[6];
  int i;
  int offset = VOLUME/2;
  size_t size = sizeof(dev_spinor_half);
  for(i=0; i<6; i++){
    cudaMemcpy( &(s_host[0]),s+i*offset , sizeof(short4), cudaMemcpyDeviceToHost);
  }

  for(i=0; i<6; i++){
    float helpx = sh2fl_host(s_host[i].x);
    float helpy = sh2fl_host(s_host[i].y);
    float helpz = sh2fl_host(s_host[i].z);
    float helpw = sh2fl_host(s_host[i].w);

	printf("%f, %f, %f, %f\n", helpx, helpy, helpz, helpw);
  }

}



void showspinor_half(dev_spinor_half* s, float* snorm){
  int i,j;
  
  dev_spinor_half help[6];
  dev_spinor help2[6];
  float norm;
  
  size_t size = sizeof(dev_spinor_half);
  int offset = VOLUME/2;  

  for(i=0; i<VOLUME/2; i++){
    for(j=0; j<6; j++){
      cudaMemcpy(&(help[j]), (s+i+j*offset), size, cudaMemcpyDeviceToHost);
    }
    cudaMemcpy(&(norm), (snorm+i), sizeof(float), cudaMemcpyDeviceToHost);
    
    for(j=0; j<6; j++){
      help2[j].x = half2float_host(help[j].x, norm);
      help2[j].y = half2float_host(help[j].y, norm);
      help2[j].z = half2float_host(help[j].z, norm);
      help2[j].w = half2float_host(help[j].w, norm);
    }
    
    for(j=0;j<6; j++){
      printf("(%.3f %.3f) (%.3f, %.3f) ", help2[j].x, help2[j].y, help2[j].z, help2[j].w);
    }
    printf("\n");
  }
  
}










void benchmark_half(dev_spinor_half* spin1, float* spin1_norm, dev_spinor_half* spin2, float* spin2_norm, int gridsize, int blocksize){
  
  double timeelapsed = 0.0;
  clock_t start, stop;
  int i;
  
  assert((start = clock())!=-1);
  printf("%s\n", cudaGetErrorString(cudaGetLastError()));
  printf("Applying H 1000 times\n");
  for(i=0; i<1000; i++){
  //Q_{-}
  #ifdef USETEXTURE
    bind_halfspinor_texture(spin1, spin1_norm);
  #endif
  //cudaFuncSetCacheConfig(dev_Hopping_Matrix_half, cudaFuncCachePreferL1);
    dev_Hopping_Matrix_half<<<gridsize, blocksize>>>
             (dev_gf_half, spin1, spin1_norm, spin2, spin2_norm, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0); //dev_spin_eo1 == even -> 0  
  #ifdef USETEXTURE
    unbind_halfspinor_texture();
  #endif

  #ifdef USETEXTURE
    bind_halfspinor_texture(spin2, spin2_norm);
  #endif
  //cudaFuncSetCacheConfig(dev_Hopping_Matrix_half, cudaFuncCachePreferL1);
    dev_Hopping_Matrix_half<<<gridsize, blocksize>>>
            (dev_gf_half, spin2, spin2_norm, spin1, spin1_norm, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1); 
  #ifdef USETEXTURE
    unbind_halfspinor_texture();
  #endif
 
  }  
  printf("%s\n", cudaGetErrorString(cudaGetLastError())); 
  printf("Done\n"); 
  
  assert((stop = clock())!=-1);
  timeelapsed = (double) (stop-start)/CLOCKS_PER_SEC;
  // x2 because 2x Hopping per iteration
  double benchres = 1400.0*2*(VOLUME/2)* 1000 / timeelapsed / 1.0e9;
  printf("Benchmark: %f Gflops\n", benchres); 
}











// this is the HALF eo version of the device cg inner solver 
// we invert the hermitean Q_{-} Q_{+}
extern "C" int dev_cg_eo_half(
     dev_su3_2v * gf,
     dev_spinor_half* spinin, float* spinin_norm,
     dev_spinor_half* spinout, float* spinout_norm, 
     dev_spinor_half* spin0, float* spin0_norm,
     dev_spinor_half* spin1, float* spin1_norm,
     dev_spinor_half* spin2, float* spin2_norm,
     dev_spinor_half* spin3, float* spin3_norm,
     dev_spinor_half* spin4, float* spin4_norm,
     int *grid, int * nn_grid,
     float epsfinal){
 
 
 float host_alpha, host_beta, host_dotprod, host_rk, sourcesquarenorm;
 float * dotprod, * dotprod2, * rk, * alpha, *beta;
 
 
 
 int i, gridsize;
 int maxit = max_innersolver_it;
 float eps = (float) innersolver_precision;
 int N_recalcres = 20; // after N_recalcres iterations calculate r = A x_k - b
 
 cudaError_t cudaerr;
 // this is the partitioning for the copying of fields
 dim3 blockdim(1,1);
 dim3 blockdim2(128,1,1);
 if( VOLUME/2 >= 128){
   gridsize = (int) VOLUME/2/128 + 1;
 }
 else{
   gridsize=1;
 }
 dim3 griddim2(gridsize,1,1);

 
 //this is the partitioning for the HoppingMatrix kernel
 int blockdim3=BLOCK;
 if( VOLUME/2 >= BLOCK){
   gridsize = (int)(VOLUME/2/BLOCK) + 1;
 }
 else{
   gridsize=1;
 }
 printf("gridsize = %d\n", gridsize);
 int griddim3=gridsize; 
 
 //this is the partitioning for dev_mul_one_pm...
 int blockdim4=BLOCK2;
 if( VOLUME/2 >= BLOCK2){
   gridsize = (int)(VOLUME/2/BLOCK2) + 1;
 }
 else{
   gridsize=1;
 }
 int griddim4=gridsize;  
 

 
 //Initialize some stuff
  printf("mu = %f\n", g_mu);
  dev_complex h0,h1,h2,h3,mh0, mh1, mh2, mh3;
  h0.re = (float)ka0.re;    h0.im = -(float)ka0.im;
  h1.re = (float)ka1.re;    h1.im = -(float)ka1.im;
  h2.re = (float)ka2.re;    h2.im = -(float)ka2.im;
  h3.re = (float)ka3.re;    h3.im = -(float)ka3.im;
  
  mh0.re = -(float)ka0.re;    mh0.im = (float)ka0.im;
  mh1.re = -(float)ka1.re;    mh1.im = (float)ka1.im;
  mh2.re = -(float)ka2.re;    mh2.im = (float)ka2.im;
  mh3.re = -(float)ka3.re;    mh3.im = (float)ka3.im;
  
  // try using constant mem for kappas
  cudaMemcpyToSymbol("dev_k0c", &h0, sizeof(dev_complex)) ; 
  cudaMemcpyToSymbol("dev_k1c", &h1, sizeof(dev_complex)) ; 
  cudaMemcpyToSymbol("dev_k2c", &h2, sizeof(dev_complex)) ; 
  cudaMemcpyToSymbol("dev_k3c", &h3, sizeof(dev_complex)) ;
  
  cudaMemcpyToSymbol("dev_mk0c", &mh0, sizeof(dev_complex)) ; 
  cudaMemcpyToSymbol("dev_mk1c", &mh1, sizeof(dev_complex)) ; 
  cudaMemcpyToSymbol("dev_mk2c", &mh2, sizeof(dev_complex)) ; 
  cudaMemcpyToSymbol("dev_mk3c", &mh3, sizeof(dev_complex)) ;  
  
  he_cg_init<<< 1, 1 >>> (grid, (float) g_kappa, (float)(g_mu/(2.0*g_kappa)), h0,h1,h2,h3);
  // BEWARE in dev_tm_dirac_kappa we need the true mu (not 2 kappa mu!)


  //use full volume here as we need the complete gauge field!!!
  int Vol;
   #ifndef _USE_MPI
     Vol = VOLUME;
   #else
     Vol = VOLUME+RAND;
   #endif
    
    if( Vol >= BLOCK2){
       gridsize = (int)(Vol/BLOCK2) + 1;
    }
    else{
       gridsize=1;
    }

   printf("Converting gauge to half precision... ");
   

     float2half_gaugefield <<< gridsize, BLOCK2  >>>(dev_gf, dev_gf_half, Vol);
   printf("Done\n"); 
   
   //testhalf_gf(dev_gf_half);
  
 
 
  #ifdef USETEXTURE
    //Bind texture gf
    bind_texture_gf_half(dev_gf_half);
  #endif
 
 
 // Init x,p,r for k=0
 // Allocate some numbers for host <-> device interaction
 cudaMalloc((void **) &dotprod, sizeof(float));
 cudaMalloc((void **) &dotprod2, sizeof(float));
 cudaMalloc((void **) &rk, sizeof(float));
 cudaMalloc((void **) &alpha, sizeof(float));
 cudaMalloc((void **) &beta, sizeof(float));
 printf("%s\n", cudaGetErrorString(cudaGetLastError())); 
 
 
 //init blas
 init_blas_half(VOLUME/2);
 printf("Have initialized blas for half precision\n");
 

 #ifdef RELATIVISTIC_BASIS 
   //transform to relativistic gamma basis
   to_relativistic_basis_half<<<griddim4, blockdim4>>> (spinin, spinin_norm);
 #endif

 dev_copy_spinor_field_half
    <<<griddim2, blockdim2 >>>(spinin, spinin_norm, spin0, spin0_norm);
 dev_zero_spinor_field_half
     <<<griddim2, blockdim2 >>>(spin1,spin1_norm); // x_0 = 0
 dev_copy_spinor_field_half
     <<<griddim2, blockdim2 >>>(spinin, spinin_norm, spin2, spin2_norm);
 dev_zero_spinor_field_half
     <<<griddim2, blockdim2 >>>(spin3, spin3_norm);
 printf("%s\n", cudaGetErrorString(cudaGetLastError()));
 
 
 //test_spinor_normalization(spin2, spin2_norm);
 //showspinor_half(spinin, spinin_norm);
 
 
 //relative precision -> get initial residue
 sourcesquarenorm = squarenorm_half(spinin, spinin_norm);
 printf("with squarenorm: %f\n", sourcesquarenorm);
 sourcesquarenorm = dotprod_half(spinin, spinin_norm,spinin, spinin_norm);
 printf("with dotprod: %f\n", sourcesquarenorm);
 
 host_rk = sourcesquarenorm; //for use in main loop
 printf("Squarenorm Source:\t%.8e\n", sourcesquarenorm);
 printf("%s\n", cudaGetErrorString(cudaGetLastError()));
 
 
 /*
 // small benchmark for half /////////////
 benchmark_half(spin2, spin2_norm, spin3, spin3_norm, griddim3,blockdim3);
 //exit(0);
 /////////////////////////////////////////
 */
 
 
 printf("Entering inner solver cg-loop\n");
 for(i=0;i<maxit;i++){ //MAIN LOOP
  
  // Q_{-}Q{+}
  #ifndef _USE_MPI
    dev_Qtm_pm_psi_half(spin2, spin2_norm, spin3, spin3_norm, griddim3, blockdim3, griddim4, blockdim4);
  #else
    dev_Qtm_pm_psi_half_mpi(spin2, spin2_norm, spin3, spin3_norm, griddim3, blockdim3, griddim4, blockdim4);
  #endif
  if((cudaerr=cudaGetLastError()) != cudaSuccess){
    printf("%s\n", cudaGetErrorString(cudaerr));
    exit(200);
  }
  //printf("Q\n");
  //test_spinor_normalization(spin2, spin2_norm);
  //test_spinor_normalization(spin3, spin3_norm);
  
  //showspinor_half(spin3, spin3_norm);
  //exit(200);
  
 //alpha
  //host_dotprod = cublasSdot (24*VOLUME/2, (const float *) spin2, 1, (const float *) spin3, 1);
  host_dotprod = dotprod_half( spin2, spin2_norm, spin3, spin3_norm);
  host_alpha = (host_rk / host_dotprod); // alpha = r*r/ p M p
  //printf("alpha = %f\n",host_alpha);

 //r(k+1)
 //cublasSaxpy (24*VOLUME/2,-1.0*host_alpha, (const float *) spin3, 1, (float *) spin0, 1);  
  axpy_half<<<griddim4, blockdim4>>> 
       (-1.0*host_alpha, spin3, spin3_norm, spin0, spin0_norm);
  
  //showspinor_half(spin0, spin0_norm);
  //exit(200);
 
  //printf("r(k+1)\n");
  //test_spinor_normalization(spin3, spin3_norm);
  //test_spinor_normalization(spin0, spin0_norm);
 //x(k+1);
 //cublasSaxpy (24*VOLUME/2, host_alpha, (const float *) spin2,  1, (float *) spin1, 1);
 axpy_half<<<griddim4, blockdim4>>> 
       (host_alpha, spin2, spin2_norm, spin1, spin1_norm);
  
  //printf("x(k+1)\n");
  //test_spinor_normalization(spin1, spin1_norm); 
  //test_spinor_normalization(spin2, spin2_norm);
 
 if((cudaerr=cudaGetLastError()) != cudaSuccess){
    printf("%s\n", cudaGetErrorString(cudaerr));
    exit(200);
  }

  //Abbruch?
  host_dotprod = squarenorm_half(spin0, spin0_norm);
  
 if (((host_dotprod <= eps*sourcesquarenorm) && (i > maxit / 4) ) || ( host_dotprod <= epsfinal/2.)){//error-limit erreicht (epsfinal/2 sollte ausreichen um auch in double precision zu bestehen)
   break; 
 }
  printf("iter %d: err = %.8e\n", i, host_dotprod);
  
 //beta
 host_beta =host_dotprod/host_rk;
 //printf("beta = %f\n",host_beta);
 //p(k+1)
 //cublasSscal (24*VOLUME/2, host_beta, (float *)spin2, 1);
 scal_half<<<griddim4, blockdim4>>>
        (host_beta, spin2, spin2_norm);
 //printf("scal p\n");
 //test_spinor_normalization(spin2, spin2_norm);

 //cublasSaxpy (24*VOLUME/2, 1.0, (const float *) spin0,  1, (float *) spin2, 1);
 axpy_half<<<griddim4, blockdim4>>>
        (1.0, spin0, spin0_norm, spin2, spin2_norm);
 //printf("axpy p\n");
 //test_spinor_normalization(spin2, spin2_norm);
 
 host_rk = host_dotprod;
 
 // recalculate residue frome r = b - Ax
 if(((i+1) % N_recalcres) == 0){
    // r_(k+1) = Ax -b 
    printf("Recalculating residue\n");
    
    // D Ddagger   --   Ddagger = gamma5 D gamma5  for Wilson Dirac Operator
    // DO NOT USE tm_dirac_dagger_kappa here, otherwise spin2 will be overwritten!!!
      
    // Q_{-}Q{+}
    #ifndef _USE_MPI
      dev_Qtm_pm_psi_half(spin1, spin1_norm, spin3, spin3_norm, griddim3, blockdim3, griddim4, blockdim4);
    #else
      dev_Qtm_pm_psi_half_mpi(spin1, spin1_norm, spin3, spin3_norm, griddim3, blockdim3, griddim4, blockdim4);
    #endif
    if((cudaerr=cudaGetLastError()) != cudaSuccess){
      printf("%s\n", cudaGetErrorString(cudaerr));
      exit(200);
    }  
        
    
    // r = b - Ax
    //cublasSscal (24*VOLUME/2, -1.0, (float *)spin3, 1);
    scal_half<<<griddim4, blockdim4>>>
               (-1.0, spin3, spin3_norm);
    
    //cublasSaxpy (24*VOLUME/2, 1.0, (const float *) spinin,  1, (float *) spin3, 1);
    axpy_half<<<griddim4, blockdim4>>>
            (1.0, spinin, spinin_norm, spin3, spin3_norm);
    
    //cublasScopy (24*VOLUME/2, (const float *)spin3, 1, (float *)spin0, 1);
    dev_copy_spinor_field_half
     <<<griddim2, blockdim2 >>>(spin3, spin3_norm, spin0, spin0_norm);
    
   }//recalculate residue

 }//MAIN LOOP cg        
  
  
  printf("Final residue: %.6e\n",host_dotprod);
  // x_result = spin1 !
  
  //no multiplication with D^{dagger} here and no return to non-kappa basis as in dev_cg!
  dev_copy_spinor_field_half<<<griddim2, blockdim2 >>>(spin1, spin1_norm,spinout, spinout_norm);
  
  #ifdef RELATIVISTIC_BASIS 
   //transform back to tmlqcd gamma basis
   to_tmlqcd_basis_half<<<griddim4, blockdim4>>> (spinout, spinout_norm);
  #endif


  #ifdef USETEXTURE
   unbind_texture_gf_half();
  #endif
  cudaFree(dotprod);
  cudaFree(dotprod2);
  cudaFree(rk);
  cudaFree(alpha);
  cudaFree(beta);
  finalize_blas_half();
  
  return(i);
}




// this is the HALF eo version of the device cg inner solver 
// we invert the hermitean Q_{-} Q_{+}
// we use the reliable update technique and update the residue periodically in single 
// precision
extern "C" int dev_cg_eo_half_reliable(
     dev_su3_2v * gf,
     dev_spinor_half* spinin, float* spinin_norm,
     dev_spinor_half* spinout, float* spinout_norm, 
     dev_spinor_half* spin0, float* spin0_norm,
     dev_spinor_half* spin1, float* spin1_norm,
     dev_spinor_half* spin2, float* spin2_norm,
     dev_spinor_half* spin3, float* spin3_norm,
     dev_spinor_half* spin4, float* spin4_norm,
     dev_spinor* spinin_f,
     dev_spinor* spin0_f,
     dev_spinor* spin1_f,
     int *grid, int * nn_grid,
     float epsfinal, float delta){
 
 
 float host_alpha, host_beta, host_dotprod, host_rk, sourcesquarenorm;
 float * dotprod, * dotprod2, * rk, * alpha, *beta;
 
 
 
 int i, gridsize;
 int maxit = max_innersolver_it;
 float eps = (float) innersolver_precision;
 int N_recalcres = 20; // after N_recalcres iterations calculate r = A x_k - b
 
 float max = 0.0f;


 cudaError_t cudaerr;
 // this is the partitioning for the copying of fields
 dim3 blockdim(1,1);
 dim3 blockdim2(128,1,1);
 if( VOLUME/2 >= 128){
   gridsize = (int) VOLUME/2/128 + 1;
 }
 else{
   gridsize=1;
 }
 dim3 griddim2(gridsize,1,1);

 
 //this is the partitioning for the HoppingMatrix kernel
 int blockdim3=BLOCK;
 if( VOLUME/2 >= BLOCK){
   gridsize = (int)(VOLUME/2/BLOCK) + 1;
 }
 else{
   gridsize=1;
 }
 printf("gridsize = %d\n", gridsize);
 int griddim3=gridsize; 
 
 //this is the partitioning for dev_mul_one_pm...
 int blockdim4=BLOCK2;
 if( VOLUME/2 >= BLOCK2){
   gridsize = (int)(VOLUME/2/BLOCK2) + 1;
 }
 else{
   gridsize=1;
 }
 int griddim4=gridsize;  
 

 
 //Initialize some stuff
  printf("mu = %f\n", g_mu);
  dev_complex h0,h1,h2,h3,mh0, mh1, mh2, mh3;
  h0.re = (float)ka0.re;    h0.im = -(float)ka0.im;
  h1.re = (float)ka1.re;    h1.im = -(float)ka1.im;
  h2.re = (float)ka2.re;    h2.im = -(float)ka2.im;
  h3.re = (float)ka3.re;    h3.im = -(float)ka3.im;
  
  mh0.re = -(float)ka0.re;    mh0.im = (float)ka0.im;
  mh1.re = -(float)ka1.re;    mh1.im = (float)ka1.im;
  mh2.re = -(float)ka2.re;    mh2.im = (float)ka2.im;
  mh3.re = -(float)ka3.re;    mh3.im = (float)ka3.im;
  
  // try using constant mem for kappas
  cudaMemcpyToSymbol("dev_k0c", &h0, sizeof(dev_complex)) ; 
  cudaMemcpyToSymbol("dev_k1c", &h1, sizeof(dev_complex)) ; 
  cudaMemcpyToSymbol("dev_k2c", &h2, sizeof(dev_complex)) ; 
  cudaMemcpyToSymbol("dev_k3c", &h3, sizeof(dev_complex)) ;
  
  cudaMemcpyToSymbol("dev_mk0c", &mh0, sizeof(dev_complex)) ; 
  cudaMemcpyToSymbol("dev_mk1c", &mh1, sizeof(dev_complex)) ; 
  cudaMemcpyToSymbol("dev_mk2c", &mh2, sizeof(dev_complex)) ; 
  cudaMemcpyToSymbol("dev_mk3c", &mh3, sizeof(dev_complex)) ;  
  
  he_cg_init<<< 1, 1 >>> (grid, (float) g_kappa, (float)(g_mu/(2.0*g_kappa)), h0,h1,h2,h3);
  // BEWARE in dev_tm_dirac_kappa we need the true mu (not 2 kappa mu!)


  //use full volume here as we need the complete gauge field!!!
  int Vol;
   #ifndef _USE_MPI
     Vol = VOLUME;
   #else
     Vol = VOLUME+RAND;
   #endif
    
    if( Vol >= BLOCK2){
       gridsize = (int)(Vol/BLOCK2) + 1;
    }
    else{
       gridsize=1;
    }


   #ifndef LOWOUTPUT 
     printf("Converting gauge to half precision... ");
   #endif 
     float2half_gaugefield <<< gridsize, BLOCK2  >>>(dev_gf, dev_gf_half, Vol);
   #ifndef LOWOUTPUT
    printf("Done\n"); 
   #endif

  
  #ifdef USETEXTURE
    //Bind texture gf
    bind_texture_gf_half(dev_gf_half);
  #endif
 
 
 // Init x,p,r for k=0
 // Allocate some numbers for host <-> device interaction
 cudaMalloc((void **) &dotprod, sizeof(float));
 cudaMalloc((void **) &dotprod2, sizeof(float));
 cudaMalloc((void **) &rk, sizeof(float));
 cudaMalloc((void **) &alpha, sizeof(float));
 cudaMalloc((void **) &beta, sizeof(float));
#ifndef LOWOUTPUT
 printf("%s\n", cudaGetErrorString(cudaGetLastError())); 
#endif 
 
 //init blas half
 init_blas_half(VOLUME/2);
#ifndef LOWOUTPUT
 printf("Have initialized blas for half precision\n");
#endif
 //init blas single
 cublasInit();
 if((cudaerr=cudaGetLastError()) != cudaSuccess){
    printf("%s\n", cudaGetErrorString(cudaerr));
    exit(200);
  }

 #ifdef RELATIVISTIC_BASIS 
   //transform to relativistic gamma basis
   to_relativistic_basis_half<<<griddim4, blockdim4>>> (spinin, spinin_norm);
 #endif

  dev_copy_spinor_field_half
    <<<griddim2, blockdim2 >>>(spinin, spinin_norm, spin0, spin0_norm);
  dev_copy_spinor_field_half
     <<<griddim2, blockdim2 >>>(spinin, spinin_norm, spin2, spin2_norm);
  dev_zero_spinor_field_half
     <<<griddim2, blockdim2 >>>(spin1,spin1_norm); // x_0 = 0 
  dev_zero_spinor_field_half
     <<<griddim2, blockdim2 >>>(spin3, spin3_norm);
  dev_zero_spinor_field
     <<<griddim2, blockdim2 >>>(spin0_f);
  dev_zero_spinor_field
     <<<griddim2, blockdim2 >>>(spin1_f);
#ifndef LOWOUTPUT 
  printf("%s\n", cudaGetErrorString(cudaGetLastError()));
#endif 
 
 
 //relative precision -> get initial residue
 sourcesquarenorm = dotprod_half(spinin, spinin_norm,spinin, spinin_norm);
 host_rk = sourcesquarenorm; //for use in main loop
 max = sourcesquarenorm;
#ifndef LOWOUTPUT 
 printf("Squarenorm Source:\t%.8e\n", sourcesquarenorm);
 printf("%s\n", cudaGetErrorString(cudaGetLastError()));
#endif 
 

 
 printf("Entering inner solver cg-loop\n");
 for(i=0;i<maxit;i++){ //MAIN LOOP
  
  // Q_{-}Q{+}
  #ifndef _USE_MPI
    dev_Qtm_pm_psi_half(spin2, spin2_norm, spin3, spin3_norm, griddim3, blockdim3, griddim4, blockdim4);
  #else
    dev_Qtm_pm_psi_half_mpi(spin2, spin2_norm, spin3, spin3_norm, griddim3, blockdim3, griddim4, blockdim4);
  #endif
  if((cudaerr=cudaGetLastError()) != cudaSuccess){
    printf("%s\n", cudaGetErrorString(cudaerr));
    exit(200);
  }

 //alpha
  host_dotprod = dotprod_half( spin2, spin2_norm, spin3, spin3_norm);
  host_alpha = (host_rk / host_dotprod); // alpha = r*r/ p M p
 //r(k+1)
  axpy_half<<<griddim4, blockdim4>>> 
       (-1.0*host_alpha, spin3, spin3_norm, spin0, spin0_norm);
 

 //x(k+1);
 axpy_half<<<griddim4, blockdim4>>> 
       (host_alpha, spin2, spin2_norm, spin1, spin1_norm);
  
 if((cudaerr=cudaGetLastError()) != cudaSuccess){
    printf("%s\n", cudaGetErrorString(cudaerr));
    exit(200);
  }

  //Abbruch?
  host_dotprod = squarenorm_half(spin0, spin0_norm);
  
  //update max since last residual update
  if(host_dotprod > max){
   max=host_dotprod;
  }  

 if (((host_dotprod <= eps*sourcesquarenorm) && (i > maxit / 8) ) || ( host_dotprod <= epsfinal/2.)){
 //error-limit erreicht (epsfinal/2 sollte ausreichen um auch in double precision zu bestehen)
   break; 
 }
#ifndef LOWOUTPUT  
 printf("iter %d: err = %.8e\n", i, host_dotprod);
#endif  
 

  //beta
 host_beta =host_dotprod/host_rk;
 scal_half<<<griddim4, blockdim4>>>
        (host_beta, spin2, spin2_norm);
 //p(k+1)
 axpy_half<<<griddim4, blockdim4>>>
        (1.0, spin0, spin0_norm, spin2, spin2_norm);

 
 
 // recalculate residue from r = b - Ax in single precision
 if( (host_dotprod < delta*max) && (i > maxit/8) ){
    // r_(k+1) = Ax -b 
    printf("Recalculating residue\n");
   
   #ifdef RELATIVISTIC_BASIS 
   //transform to relativistic gamma basis
   to_tmlqcd_basis_half<<<griddim4, blockdim4>>> (spin1, spin1_norm);
   #endif
    
   addhalf2float_spinorfield<<<griddim4, blockdim4>>>(spin0_f, spin1, spin1_norm);
  
    // Q_{-}Q{+}
    #ifndef _USE_MPI
      dev_Qtm_pm_psi(spin0_f, spin1_f, griddim3, blockdim3, griddim4, blockdim4);
    #else
      dev_Qtm_pm_psi_mpi(spin0_f, spin1_f, griddim3, blockdim3, griddim4, blockdim4);
    #endif
    if((cudaerr=cudaGetLastError()) != cudaSuccess){
      printf("%s\n", cudaGetErrorString(cudaerr));
      exit(200);
    }  
        
    
    // r = b - Ax
    cublasSscal (24*VOLUME/2, -1.0, (float *)spin1_f, 1);
    cublasSaxpy (24*VOLUME/2, 1.0, (const float *) spinin_f,  1, (float *) spin1_f, 1);
    float2half_spinorfield<<<griddim4, blockdim4>>>(spin0, spin0_norm,spin1_f);
    #ifdef RELATIVISTIC_BASIS 
      //transform to relativistic gamma basis
      to_relativistic_basis_half<<<griddim4, blockdim4>>> (spin0, spin0_norm);
    #endif  

    
    // x = 0
    dev_zero_spinor_field_half<<<griddim2, blockdim2 >>>(spin1,spin1_norm); 


    // update max since this residual update
    host_dotprod = squarenorm_half(spin0, spin0_norm);
    max = host_dotprod;

dev_zero_spinor_field_half
     <<<griddim2, blockdim2 >>>(spin3, spin3_norm);


/*
    float host_gamma = dotprod_half( spin0, spin0_norm, spin2, spin2_norm);
    //p = p - (p,r) r
    axpy_half<<<griddim4, blockdim4>>>
        (-host_gamma/host_dotprod, spin0, spin0_norm, spin2, spin2_norm);

    host_alpha = 0.0;
*/
/*
    // r = b
    dev_copy_spinor_field_half
    <<<griddim2, blockdim2 >>>(spinin, spinin_norm, spin0, spin0_norm);
*/  
   }//recalculate residue
  

  host_rk = host_dotprod;
 
 }//MAIN LOOP cg        
  
  
  printf("Final residue: %.6e\n",host_dotprod);
  // x_result = spin1 !
  
  //no multiplication with D^{dagger} here and no return to non-kappa basis as in dev_cg!
  dev_copy_spinor_field_half<<<griddim2, blockdim2 >>>(spin1, spin1_norm,spinout, spinout_norm);
  #ifdef RELATIVISTIC_BASIS 
   //transform back to tmlqcd gamma basis
   to_tmlqcd_basis_half<<<griddim4, blockdim4>>> (spinout, spinout_norm);
  #endif

  //add latest x_k (low precision) to high precision accumulate
  #ifdef RELATIVISTIC_BASIS 
   //transform back to tmlqcd gamma basis
   to_tmlqcd_basis_half<<<griddim4, blockdim4>>> (spin1, spin1_norm);
  #endif 
  addhalf2float_spinorfield<<<griddim4, blockdim4>>>(spin0_f, spin1, spin1_norm);

  #ifdef USETEXTURE
   unbind_texture_gf_half();
  #endif
  cudaFree(dotprod);
  cudaFree(dotprod2);
  cudaFree(rk);
  cudaFree(alpha);
  cudaFree(beta);
  finalize_blas_half();
  cublasShutdown();
  return(i);
}








// this is the HALF eo version of the device cg inner solver 
// we invert the hermitean Q_{-} Q_{+}
extern "C" int dev_cg_half_reliable_update(
     dev_su3_2v * gf,
     dev_spinor_half* spinin, float* spinin_norm,
     dev_spinor_half* spinout, float* spinout_norm, 
     dev_spinor_half* spin0, float* spin0_norm,
     dev_spinor_half* spin1, float* spin1_norm,
     dev_spinor_half* spin2, float* spin2_norm,
     dev_spinor_half* spin3, float* spin3_norm,
     dev_spinor_half* spin4, float* spin4_norm,
     int *grid, int * nn_grid,
     float epsfinal){
 
 
 float host_alpha, host_beta, host_dotprod, host_rk, sourcesquarenorm;
 float * dotprod, * dotprod2, * rk, * alpha, *beta;
 
 
 
 int i, gridsize;
 int maxit = max_innersolver_it;
 float eps = (float) innersolver_precision;
 int N_recalcres = 1000; // after N_recalcres iterations calculate r = A x_k - b
 
 cudaError_t cudaerr;
 // this is the partitioning for the copying of fields
 dim3 blockdim(1,1);
 dim3 blockdim2(128,1,1);
 if( VOLUME/2 >= 128){
   gridsize = (int) VOLUME/2/128 + 1;
 }
 else{
   gridsize=1;
 }
 dim3 griddim2(gridsize,1,1);

 
 //this is the partitioning for the HoppingMatrix kernel
 int blockdim3=BLOCK;
 if( VOLUME/2 >= BLOCK){
   gridsize = (int)(VOLUME/2/BLOCK) + 1;
 }
 else{
   gridsize=1;
 }
 printf("gridsize = %d\n", gridsize);
 int griddim3=gridsize; 
 
 //this is the partitioning for dev_mul_one_pm...
 int blockdim4=BLOCK2;
 if( VOLUME/2 >= BLOCK2){
   gridsize = (int)(VOLUME/2/BLOCK2) + 1;
 }
 else{
   gridsize=1;
 }
 int griddim4=gridsize;  
 

 
 //Initialize some stuff
  printf("mu = %f\n", g_mu);
  dev_complex h0,h1,h2,h3,mh0, mh1, mh2, mh3;
  h0.re = (float)ka0.re;    h0.im = -(float)ka0.im;
  h1.re = (float)ka1.re;    h1.im = -(float)ka1.im;
  h2.re = (float)ka2.re;    h2.im = -(float)ka2.im;
  h3.re = (float)ka3.re;    h3.im = -(float)ka3.im;
  
  mh0.re = -(float)ka0.re;    mh0.im = (float)ka0.im;
  mh1.re = -(float)ka1.re;    mh1.im = (float)ka1.im;
  mh2.re = -(float)ka2.re;    mh2.im = (float)ka2.im;
  mh3.re = -(float)ka3.re;    mh3.im = (float)ka3.im;
  
  // try using constant mem for kappas
  cudaMemcpyToSymbol("dev_k0c", &h0, sizeof(dev_complex)) ; 
  cudaMemcpyToSymbol("dev_k1c", &h1, sizeof(dev_complex)) ; 
  cudaMemcpyToSymbol("dev_k2c", &h2, sizeof(dev_complex)) ; 
  cudaMemcpyToSymbol("dev_k3c", &h3, sizeof(dev_complex)) ;
  
  cudaMemcpyToSymbol("dev_mk0c", &mh0, sizeof(dev_complex)) ; 
  cudaMemcpyToSymbol("dev_mk1c", &mh1, sizeof(dev_complex)) ; 
  cudaMemcpyToSymbol("dev_mk2c", &mh2, sizeof(dev_complex)) ; 
  cudaMemcpyToSymbol("dev_mk3c", &mh3, sizeof(dev_complex)) ;  
  
  he_cg_init<<< 1, 1 >>> (grid, (float) g_kappa, (float)(g_mu/(2.0*g_kappa)), h0,h1,h2,h3);
  // BEWARE in dev_tm_dirac_kappa we need the true mu (not 2 kappa mu!)


  //use full volume here as we need the complete gauge field!!!
  if( VOLUME >= BLOCK2){
     gridsize = (int)(VOLUME/BLOCK2) + 1;
  }
  else{
         gridsize=1;
  }

   printf("Converting gauge to half precision... ");
     float2half_gaugefield <<< gridsize, BLOCK2  >>>(dev_gf, dev_gf_half, VOLUME);
   printf("Done\n"); 
   
   //testhalf_gf(dev_gf_half);
  
 
 
  #ifdef USETEXTURE
    //Bind texture gf
    bind_texture_gf_half(dev_gf_half);
  #endif
 
 
 // Init x,p,r for k=0
 // Allocate some numbers for host <-> device interaction
 cudaMalloc((void **) &dotprod, sizeof(float));
 cudaMalloc((void **) &dotprod2, sizeof(float));
 cudaMalloc((void **) &rk, sizeof(float));
 cudaMalloc((void **) &alpha, sizeof(float));
 cudaMalloc((void **) &beta, sizeof(float));
 printf("%s\n", cudaGetErrorString(cudaGetLastError())); 
 
 
 //init blas
 init_blas_half(VOLUME/2);
 printf("Have initialized blas for half precision\n");
 
 

 dev_copy_spinor_field_half
    <<<griddim2, blockdim2 >>>(spinin, spinin_norm, spin0, spin0_norm);
 dev_zero_spinor_field_half
     <<<griddim2, blockdim2 >>>(spin1,spin1_norm); // x_0 = 0
 dev_copy_spinor_field_half
     <<<griddim2, blockdim2 >>>(spinin, spinin_norm, spin2, spin2_norm);
 dev_zero_spinor_field_half
     <<<griddim2, blockdim2 >>>(spin3, spin3_norm);
 printf("%s\n", cudaGetErrorString(cudaGetLastError()));
 
 
 //test_spinor_normalization(spin2, spin2_norm);
 //showspinor_half(spinin, spinin_norm);
 
 
 //relative precision -> get initial residue
 sourcesquarenorm = squarenorm_half(spinin, spinin_norm);
 printf("with squarenorm: %f\n", sourcesquarenorm);
 sourcesquarenorm = dotprod_half(spinin, spinin_norm,spinin, spinin_norm);
 printf("with dotprod: %f\n", sourcesquarenorm);
 
 host_rk = sourcesquarenorm; //for use in main loop
 printf("Squarenorm Source:\t%.8e\n", sourcesquarenorm);
 printf("%s\n", cudaGetErrorString(cudaGetLastError()));
 
 
 
 printf("Entering inner solver cg-loop\n");
 for(i=0;i<maxit;i++){ //MAIN LOOP
  
  // Q_{-}Q{+}
  dev_Qtm_pm_psi_half(spin2, spin2_norm, spin3, spin3_norm, griddim3, blockdim3, griddim4, blockdim4);
  if((cudaerr=cudaGetLastError()) != cudaSuccess){
    printf("%s\n", cudaGetErrorString(cudaerr));
    exit(200);
  }
  //printf("Q\n");
  //test_spinor_normalization(spin2, spin2_norm);
  //test_spinor_normalization(spin3, spin3_norm);
  
  
 //alpha
  //host_dotprod = cublasSdot (24*VOLUME/2, (const float *) spin2, 1, (const float *) spin3, 1);
  host_dotprod = dotprod_half( spin2, spin2_norm, spin3, spin3_norm);
  host_alpha = (host_rk / host_dotprod); // alpha = r*r/ p M p
  //printf("alpha = %f\n",host_alpha);

 //r(k+1)
 //cublasSaxpy (24*VOLUME/2,-1.0*host_alpha, (const float *) spin3, 1, (float *) spin0, 1);  
  axpy_half<<<griddim4, blockdim4>>> 
       (-1.0*host_alpha, spin3, spin3_norm, spin0, spin0_norm);
  
  //showspinor_half(spin0, spin0_norm);
  //exit(200);
 
  //printf("r(k+1)\n");
  //test_spinor_normalization(spin3, spin3_norm);
  //test_spinor_normalization(spin0, spin0_norm);
 //x(k+1);
 //cublasSaxpy (24*VOLUME/2, host_alpha, (const float *) spin2,  1, (float *) spin1, 1);
 axpy_half<<<griddim4, blockdim4>>> 
       (host_alpha, spin2, spin2_norm, spin1, spin1_norm);
  
  //printf("x(k+1)\n");
  //test_spinor_normalization(spin1, spin1_norm); 
  //test_spinor_normalization(spin2, spin2_norm);
 
 if((cudaerr=cudaGetLastError()) != cudaSuccess){
    printf("%s\n", cudaGetErrorString(cudaerr));
    exit(200);
  }

  //Abbruch?
  host_dotprod = squarenorm_half(spin0, spin0_norm);
  
 if (((host_dotprod <= eps*sourcesquarenorm) && (i > maxit / 4) ) || ( host_dotprod <= epsfinal/2.)){//error-limit erreicht (epsfinal/2 sollte ausreichen um auch in double precision zu bestehen)
   break; 
 }
  printf("iter %d: err = %.8e\n", i, host_dotprod);
  
 //beta
 host_beta =host_dotprod/host_rk;
 //printf("beta = %f\n",host_beta);
 //p(k+1)
 //cublasSscal (24*VOLUME/2, host_beta, (float *)spin2, 1);
 scal_half<<<griddim4, blockdim4>>>
        (host_beta, spin2, spin2_norm);
 //printf("scal p\n");
 //test_spinor_normalization(spin2, spin2_norm);

 //cublasSaxpy (24*VOLUME/2, 1.0, (const float *) spin0,  1, (float *) spin2, 1);
 axpy_half<<<griddim4, blockdim4>>>
        (1.0, spin0, spin0_norm, spin2, spin2_norm);
 //printf("axpy p\n");
 //test_spinor_normalization(spin2, spin2_norm);
 
 host_rk = host_dotprod;
 
 // recalculate residue frome r = b - Ax
 if(((i+1) % N_recalcres) == 0){
    // r_(k+1) = Ax -b 
    printf("Recalculating residue\n");
    
    // D Ddagger   --   Ddagger = gamma5 D gamma5  for Wilson Dirac Operator
    // DO NOT USE tm_dirac_dagger_kappa here, otherwise spin2 will be overwritten!!!
      
    // Q_{-}Q{+}
    dev_Qtm_pm_psi_half(spin1, spin1_norm, spin3, spin3_norm, griddim3, blockdim3, griddim4, blockdim4);
    if((cudaerr=cudaGetLastError()) != cudaSuccess){
      printf("%s\n", cudaGetErrorString(cudaerr));
      exit(200);
    }  
        
    
    // r = b - Ax
    //cublasSscal (24*VOLUME/2, -1.0, (float *)spin3, 1);
    scal_half<<<griddim4, blockdim4>>>
               (-1.0, spin3, spin3_norm);
    
    //cublasSaxpy (24*VOLUME/2, 1.0, (const float *) spinin,  1, (float *) spin3, 1);
    axpy_half<<<griddim4, blockdim4>>>
            (1.0, spinin, spinin_norm, spin3, spin3_norm);
    
    //cublasScopy (24*VOLUME/2, (const float *)spin3, 1, (float *)spin0, 1);
    dev_copy_spinor_field_half
     <<<griddim2, blockdim2 >>>(spin3, spin3_norm, spin0, spin0_norm);
    
   }//recalculate residue

 }//MAIN LOOP cg        
  
  
  printf("Final residue: %.6e\n",host_dotprod);
  // x_result = spin1 !
  
  //no multiplication with D^{dagger} here and no return to non-kappa basis as in dev_cg!
  dev_copy_spinor_field_half<<<griddim2, blockdim2 >>>(spin1, spin1_norm,spinout, spinout_norm);
  
  #ifdef USETEXTURE
   unbind_texture_gf_half();
  #endif
  cudaFree(dotprod);
  cudaFree(dotprod2);
  cudaFree(rk);
  cudaFree(alpha);
  cudaFree(beta);
  finalize_blas_half();
  
  return(i);
}








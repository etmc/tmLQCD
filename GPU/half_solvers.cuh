








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
     REAL epsfinal){
 
 
 REAL host_alpha, host_beta, host_dotprod, host_rk, sourcesquarenorm;
 REAL * dotprod, * dotprod2, * rk, * alpha, *beta;
 
 
 
 int i, gridsize;
 int maxit = max_innersolver_it;
 REAL eps = (REAL) innersolver_precision;
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
 
 
 size_t size2 = sizeof(float4)*6*VOLUME/2;
 
 
 //Initialize some stuff
  printf("mu = %f\n", g_mu);
  dev_complex h0,h1,h2,h3,mh0, mh1, mh2, mh3;
  h0.re = (REAL)ka0.re;    h0.im = -(REAL)ka0.im;
  h1.re = (REAL)ka1.re;    h1.im = -(REAL)ka1.im;
  h2.re = (REAL)ka2.re;    h2.im = -(REAL)ka2.im;
  h3.re = (REAL)ka3.re;    h3.im = -(REAL)ka3.im;
  
  mh0.re = -(REAL)ka0.re;    mh0.im = (REAL)ka0.im;
  mh1.re = -(REAL)ka1.re;    mh1.im = (REAL)ka1.im;
  mh2.re = -(REAL)ka2.re;    mh2.im = (REAL)ka2.im;
  mh3.re = -(REAL)ka3.re;    mh3.im = (REAL)ka3.im;
  
  // try using constant mem for kappas
  cudaMemcpyToSymbol("dev_k0c", &h0, sizeof(dev_complex)) ; 
  cudaMemcpyToSymbol("dev_k1c", &h1, sizeof(dev_complex)) ; 
  cudaMemcpyToSymbol("dev_k2c", &h2, sizeof(dev_complex)) ; 
  cudaMemcpyToSymbol("dev_k3c", &h3, sizeof(dev_complex)) ;
  
  cudaMemcpyToSymbol("dev_mk0c", &mh0, sizeof(dev_complex)) ; 
  cudaMemcpyToSymbol("dev_mk1c", &mh1, sizeof(dev_complex)) ; 
  cudaMemcpyToSymbol("dev_mk2c", &mh2, sizeof(dev_complex)) ; 
  cudaMemcpyToSymbol("dev_mk3c", &mh3, sizeof(dev_complex)) ;  
  
  he_cg_init<<< 1, 1 >>> (grid, (REAL) g_kappa, (REAL)(g_mu/(2.0*g_kappa)), h0,h1,h2,h3);
  // BEWARE in dev_tm_dirac_kappa we need the true mu (not 2 kappa mu!)
 
   printf("Converting gauge to half precision... ");
     float2half_gaugefield <<< gridsize, BLOCK2  >>>(dev_gf, dev_gf_half, VOLUME);
   printf("Done\n"); 
   
   testhalf_gf(gf);
  
 
 
  #ifdef USETEXTURE
    //Bind texture gf
    bind_texture_gf(gf);
  #endif
 
 
 // Init x,p,r for k=0
 // Allocate some numbers for host <-> device interaction
 cudaMalloc((void **) &dotprod, sizeof(REAL));
 cudaMalloc((void **) &dotprod2, sizeof(REAL));
 cudaMalloc((void **) &rk, sizeof(REAL));
 cudaMalloc((void **) &alpha, sizeof(REAL));
 cudaMalloc((void **) &beta, sizeof(REAL));
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
 
 
 

 //relative precision -> get initial residue
 sourcesquarenorm = squarenorm_half(spinin_norm);
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
  
  
  
 //alpha
  //host_dotprod = cublasSdot (24*VOLUME/2, (const float *) spin2, 1, (const float *) spin3, 1);
  host_dotprod = dotprod_half( spin2, spin2_norm, spin3, spin3_norm);
  host_alpha = (host_rk / host_dotprod); // alpha = r*r/ p M p
   
 //r(k+1)
 //cublasSaxpy (24*VOLUME/2,-1.0*host_alpha, (const float *) spin3, 1, (float *) spin0, 1);  
  axpy_half<<<griddim4, blockdim4>>> 
       (-1.0*host_alpha, spin3, spin3_norm, spin0, spin0_norm);

 //x(k+1);
 //cublasSaxpy (24*VOLUME/2, host_alpha, (const float *) spin2,  1, (float *) spin1, 1);
 axpy_half<<<griddim4, blockdim4>>> 
       (-1.0*host_alpha, spin2, spin2_norm, spin1, spin1_norm);
 
 
 if((cudaerr=cudaGetLastError()) != cudaSuccess){
    printf("%s\n", cudaGetErrorString(cudaerr));
    exit(200);
  }

  //Abbruch?
  host_dotprod = squarenorm_half(spin0_norm);
  
 if (((host_dotprod <= eps*sourcesquarenorm) && (i > maxit / 4) ) || ( host_dotprod <= epsfinal/2.)){//error-limit erreicht (epsfinal/2 sollte ausreichen um auch in double precision zu bestehen)
   break; 
 }
  printf("iter %d: err = %.8e\n", i, host_dotprod);
  
 //beta
 host_beta =host_dotprod/host_rk;
 //p(k+1)
 //cublasSscal (24*VOLUME/2, host_beta, (float *)spin2, 1);
 scal_half<<<griddim4, blockdim4>>>
        (host_beta, spin2, spin2_norm);
 //cublasSaxpy (24*VOLUME/2, 1.0, (const float *) spin0,  1, (float *) spin2, 1);
 axpy_half<<<griddim4, blockdim4>>>
        (1.0, spin0, spin0_norm, spin2, spin2_norm);

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
   unbind_texture_gf();
  #endif
  cudaFree(dotprod);
  cudaFree(dotprod2);
  cudaFree(rk);
  cudaFree(alpha);
  cudaFree(beta);
  finalize_blas_half();
  
  return(i);
}



















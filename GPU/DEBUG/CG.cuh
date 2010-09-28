
// this is an exact implementation of the CG according to cg_her_nd()


////////////////////////
// CONJUGATE GRADIENT //
////////////////////////

// for the odd field after even/odd-preconditioning
// single precision on GPU

int alt_cg_eo_nd (dev_su3_2v * gf,
                  dev_spinor * P_up, dev_spinor * P_dn,
                  dev_spinor * Q_up, dev_spinor * Q_dn,
                  double eps_sq) {
  
  // P_up/dn  can be used as auxiliary field to work on, as it is not later used
  // Q_up/dn  can be used as feedback or initial guess, or if not, also as auxiliary field
  
  printf("This is the ALTERNATIVE CG on the device!\n");
  
  
  /////////////////////
  // LOCAL VARIABLES //			// single precision
  /////////////////////
  
  // some variables
  int iteration;			// iteration counter
  cudaError_t cudaerr;			// CUDA errors
  cublasStatus cublasstatus;		// CUBLAS status
  
  
  
  
  // some variables
  int maxit = 1000;
  float eps_rel = 0.0001;
  // int maxit = max_innersolver_it;		// maximal number of inner iterations per one outer iteration
  // float eps_rel = (float) innersolver_precision;// precision for the inner solver
  float eps_abs = (float) eps_sq/2.0;
  
  int N_sites  =    VOLUME/2;
  int N_floats = 24*VOLUME/2;			// (single precision) CUBLAS functions get the number of floats as input
  
  int N_recalcres = 40;
  
  size_t dev_spinsize = 6*VOLUME/2*sizeof(dev_spinor);
  
  
  
  /////////////////////////////////////////////
  // CUDA block- and gridsize specifications //			// why here initialized and passed ?? 	// could also be done locally in  matrix_multiplication32
  /////////////////////////////////////////////
  
  int gridsize;			// auxiliary
  int blocksize;		// auxiliary
  
  blocksize = 128;
  int blockdim1 = blocksize;					// here:	dev_zero_spinor_field , dev_copy_spinor_field
  int griddim1  = (int) (VOLUME/2/blocksize) + 1;
  
  blocksize = 128;
  int blockdim2 = blocksize;					// passed:	dev_Hopping_Matrix
  int griddim2  = (int) (VOLUME/2/blocksize) + 1;
  
  blocksize = 128;
  int blockdim3 = blocksize;					// passed:	dev_mul_one_pm_imubar_gamma5
  int griddim3  = (int) (VOLUME/2/blocksize) + 1;
  
  blocksize = 128;
  int blockdim4 = blocksize;					// passed:	dev_gamma5
  int griddim4  = (int) (VOLUME/2/blocksize) + 1;
  
  blocksize = 128;
  int blockdim5 = blocksize;					// passed:	dev_copy_spinor_field
  int griddim5  = (int) (VOLUME/2/blocksize) + 1;
  
  
  
  
  /////////////////
  // ASSIGNMENTS //
  /////////////////


  
  
  /////////////////////
  // INITIALIZATIONS //
  /////////////////////
  
  
  // Initialize some stuff
  dev_complex h0, h1, h2, h3, mh0, mh1, mh2, mh3;
  
  h0.re  =  (REAL) ka0.re;	h0.im  = -(REAL) ka0.im;	// ka{0-4} are defined in boundary.c
  h1.re  =  (REAL) ka1.re;	h1.im  = -(REAL) ka1.im;	// what is the meaning?
  h2.re  =  (REAL) ka2.re;	h2.im  = -(REAL) ka2.im;
  h3.re  =  (REAL) ka3.re;	h3.im  = -(REAL) ka3.im;
  
  mh0.re = -(REAL) ka0.re;	mh0.im =  (REAL) ka0.im;
  mh1.re = -(REAL) ka1.re;	mh1.im =  (REAL) ka1.im;
  mh2.re = -(REAL) ka2.re;	mh2.im =  (REAL) ka2.im;
  mh3.re = -(REAL) ka3.re;	mh3.im =  (REAL) ka3.im;
  
  // try using constant mem for kappas		// constant memory is cached!
  cudaMemcpyToSymbol("dev_k0c", &h0, sizeof(dev_complex));
  cudaMemcpyToSymbol("dev_k1c", &h1, sizeof(dev_complex));
  cudaMemcpyToSymbol("dev_k2c", &h2, sizeof(dev_complex));
  cudaMemcpyToSymbol("dev_k3c", &h3, sizeof(dev_complex));
  
  cudaMemcpyToSymbol("dev_mk0c", &mh0, sizeof(dev_complex));
  cudaMemcpyToSymbol("dev_mk1c", &mh1, sizeof(dev_complex));
  cudaMemcpyToSymbol("dev_mk2c", &mh2, sizeof(dev_complex));
  cudaMemcpyToSymbol("dev_mk3c", &mh3, sizeof(dev_complex));
  
  		// debug	// CUDA
    		#ifdef CUDA_DEBUG
    		  CUDA_CHECK("CUDA error in cg_eo_nd(). Trying to use constant memory for strange stuff failed.", "Using constant memory for strange stuff.");
    		#endif
  		
  
  // bind texture gf
  bind_texture_gf(gf);						// needed for subfunctions of dev_Hopping_Matrix(...)
  								//	e.g. dev_reconstructgf_2vtexref(...), dev_reconstructgf_8texref(...), 
  								//	in general in functions  dev_reconstructgf[...]  with  "tex1Dfetch(gf_tex[...]"
  		// debug	// CUDA
    		#ifdef CUDA_DEBUG
    		  CUDA_CHECK("CUDA error in bind_texture_gf(). Binding GF to texture failed.", "GF bound to texture.");
    		#endif
  
  
  // "he" = "host entry"
  he_cg_init<<< 1, 1 >>> (dev_grid, (REAL) g_kappa, (REAL)(g_mu/(2.0*g_kappa)), h0, h1, h2, h3);
  // BEWARE in dev_tm_dirac_kappa we need the true mu (not 2 kappa mu!)	// ??
  
  	// dev_LX, dev_LY, dev_LZ, dev_T, dev_VOLUME  =  grid[5]  =  dev_grid[5]
  	//	dev_VOLUME  is necessary for many kernel functions as for instance  dev_gamma5()
  	// initializes  mu, kappa and twokappamu  on the device
  	// initializes the strange  dev_k{0-3}, dev_mk{0-3}  as  derived from the  ka{0-3}  from boundary.c
  
  		// debug	// kernel
  		#ifdef CUDA_DEBUG
  		  CUDA_KERNEL_CHECK("Kernel error in he_cg_init(). Couldn't initialize some stuff.", "he_cg_init() succeeded.");
  		#endif
  
  		// debug	// check stuff on device
  		#ifdef STUFF_DEBUG
  		  int host_check_LX, host_check_LY, host_check_LZ, host_check_T, host_check_VOLUME;
  		  cudaMemcpyFromSymbol(&host_check_LX, dev_LX, sizeof(int));
  		  cudaMemcpyFromSymbol(&host_check_LY, dev_LY, sizeof(int));
  		  cudaMemcpyFromSymbol(&host_check_LZ, dev_LZ, sizeof(int));
  		  cudaMemcpyFromSymbol(&host_check_T, dev_T, sizeof(int));
  		  cudaMemcpyFromSymbol(&host_check_VOLUME, dev_VOLUME, sizeof(int));
  		  printf("\teven_odd_flag = %i\n", even_odd_flag);
  		  printf("\tOn device:\n");
  		  printf("\tdev_LX = %i\n", host_check_LX);
  		  printf("\tdev_LY = %i\n", host_check_LY);
  		  printf("\tdev_LZ = %i\n", host_check_LZ);
  		  printf("\tdev_T = %i\n", host_check_T);
  		  printf("\tdev_VOLUME = %i/2 ?!= %i\n", host_check_LX*host_check_LY*host_check_LZ*host_check_T, host_check_VOLUME);
  		  
  		  float host_check_mu, host_check_kappa, host_check_twokappamu;
  		  cudaMemcpyFromSymbol(&host_check_mu, mu, sizeof(float));
  		  cudaMemcpyFromSymbol(&host_check_kappa, kappa, sizeof(float));
  		  cudaMemcpyFromSymbol(&host_check_twokappamu, twokappamu, sizeof(float));
  		  // printf("\tOn device:\n");
  		  // printf("\tmu = %f\n", host_check_mu);		// not needed for the nd case
  		  printf("\tkappa = %f\n", host_check_kappa);
  		  // printf("\ttwokappamu = %f\n", host_twokappamu);
  		#endif
  
  // additional
  dev_spinor * dev_spin_up;
  dev_spinor * dev_spin_dn;
  dev_spinor * dev_spin4_up;
  dev_spinor * dev_spin4_dn;
  dev_spinor * dev_spin5_up;
  dev_spinor * dev_spin5_dn;
  cudaMalloc((void **) &dev_spin_up, dev_spinsize);
  cudaMalloc((void **) &dev_spin_dn, dev_spinsize);
  cudaMalloc((void **) &dev_spin4_up, dev_spinsize);
  cudaMalloc((void **) &dev_spin4_dn, dev_spinsize);
  cudaMalloc((void **) &dev_spin5_up, dev_spinsize);
  cudaMalloc((void **) &dev_spin5_dn, dev_spinsize);
  
  
  he_cg_init_nd_additional<<<1,1>>> (g_mubar, g_epsbar);
  
  		// debug	// kernel
  		#ifdef CUDA_DEBUG
  		  CUDA_KERNEL_CHECK("Kernel error in he_cg_init_nd_additional(). Couldn't initialize some stuff.", "he_cg_init_nd_additional() succeeded.");
  		#endif
  
  		// debug	// check mubar and epsbar on host and device
  		#ifdef STUFF_DEBUG
  		  // printf("\tOn host:\n");
  		  // printf("\tg_mubar = %f\n", g_mubar);
  		  // printf("\tg_epsbar = %f\n", g_epsbar);
  		  
  		  float host_check_mubar, host_check_epsbar;
  		  cudaMemcpyFromSymbol(&host_check_mubar, mubar, sizeof(float));
  		  cudaMemcpyFromSymbol(&host_check_epsbar, epsbar, sizeof(float));
  		  printf("\tOn device:\n");
  		  printf("\tmubar = %f\n", host_check_mubar);
  		  printf("\tepsbar = %f\n", host_check_epsbar);
  		#endif
  
  
  
  
  
  // init CUBLAS
  // cublasInit();
  
  		// debug	// CUBLAS helper function
  		#ifdef CUDA_DEBUG
  		  CUBLAS_HELPER_CHECK(cublasInit(), "CUBLAS error in cublasInit(). Couldn't initialize CUBLAS.", "CUBLAS initialized.");
  		#endif
  
  
  
  float squarenorm_up, squarenorm_dn, squarenorm;
  float normsp_up, normsp_dn, normsp;				// just to check if zero
  float normsq_up, normsq_dn, normsq;				// for the algorithm
  float pro_up, pro_dn, pro;					// 
  float alpha_cg, beta_cg;					// alpha, beta
  float err_up, err_dn, err;					// r(k+1)*r(k+1)
  
  
  
  								//////////////////////////
  ///////////////						//  dev_spin  = x(k)	//
  // ALGORITHM //						//  dev_spin1 = r(k)	//
  ///////////////						//  dev_spin2 = d(k)	//
  								//  dev_spin4 = A*d(k)	//
  								//////////////////////////
  
  
  // P = 0
  dev_zero_spinor_field<<<griddim1, blockdim1>>>(P_up);				// P_up = 0
  dev_zero_spinor_field<<<griddim1, blockdim1>>>(P_dn);				// P_dn = 0
  
  
  // squarenorm = (Q_up)^2 + (Q_dn)^2
  squarenorm_up = cublasSdot(N_floats, (float *) Q_up, 1, (float *) Q_up, 1);
  squarenorm_dn = cublasSdot(N_floats, (float *) Q_dn, 1, (float *) Q_dn, 1);
  squarenorm = squarenorm_up + squarenorm_dn;
  
  // x(0) = P = 0
  dev_copy_spinor_field<<<griddim1, blockdim1>>>(P_up, dev_spin_up);		// dev_spin_up = P_up
  dev_copy_spinor_field<<<griddim1, blockdim1>>>(P_dn, dev_spin_dn);		// dev_spin_dn = P_dn
  
  // normsp = (P_up)^2 + (P_dn)^2
  normsp_up = cublasSdot(N_floats, (float *) P_up, 1, (float *) P_up, 1);
  normsp_dn = cublasSdot(N_floats, (float *) P_dn, 1, (float *) P_dn, 1);
  normsp = normsp_up + normsp_dn;
  
  // why and what ??
  dev_copy_spinor_field<<<griddim1, blockdim1>>>(Q_up, dev_spin5_up);		// dev_spin5_up = Q_up
  dev_copy_spinor_field<<<griddim1, blockdim1>>>(Q_dn, dev_spin5_dn);		// dev_spin5_dn = Q_dn
  
  
  if (normsp == 0) {
  
    printf("Yes, normsp = 0!\n");
    
    // r(0) = p(0) = Q
    dev_copy_spinor_field<<<griddim1, blockdim1>>>(dev_spin5_up, dev_spin1_up);	// dev_spin1_up = dev_spin5_up = Q_up
    dev_copy_spinor_field<<<griddim1, blockdim1>>>(dev_spin5_dn, dev_spin1_dn);	// dev_spin1_dn = dev_spin5_dn = Q_dn
    dev_copy_spinor_field<<<griddim1, blockdim1>>>(dev_spin5_up, dev_spin2_up);	// dev_spin2_up = dev_spin5_up = Q_up
    dev_copy_spinor_field<<<griddim1, blockdim1>>>(dev_spin5_dn, dev_spin2_dn);	// dev_spin2_dn = dev_spin5_dn = Q_dn
    
    // normsq = (Q_up)^2 + (Q_dn)^2	// first residual for (initial guess = 0)
    normsq_up = cublasSdot(N_floats, (float *) Q_up, 1, (float *) Q_up, 1);
    normsq_dn = cublasSdot(N_floats, (float *) Q_dn, 1, (float *) Q_dn, 1);
    normsq = normsq_up + normsq_dn;
    
  }
  
  		// debug	// CUBLAS core function
		#ifdef CUDA_DEBUG
		  CUBLAS_CORE_CHECK("CUBLAS error in alt_cg_eo_nd(). Calculating initial residue failed.", "Initial residue calculated.");
		#endif
  
  		// debug
  		printf("Initial inner residue: %.8e\n", squarenorm);
  
  
  
  //////////
  // LOOP //
  //////////
  
  
  		// debug
  		printf("Entering inner CG loop.\n");
  
  
  for (iteration = 0; iteration < maxit; iteration++) {
    
    
    // A*d(k)
    #ifndef MATRIX_DEBUG
    
    matrix_multiplication32(dev_spin4_up, dev_spin4_dn,
                            dev_spin2_up, dev_spin2_dn,
                            griddim2, blockdim2,	// for calling some kernels as subfunctions
                            griddim3, blockdim3,	// ...
                            griddim4, blockdim4,
                            griddim5, blockdim5);
    
    
  		// debug	// CUDA		// also other stuff ?!
  		#ifdef CUDA_DEBUG
  		  CUDA_CHECK("CUDA error in matrix_muliplication32(). Applying the matrix on GPU failed.", "The matrix was applied on GPU.");
  		#endif
  		
    #else
    
    // debug	// apply the host matrix on trial
    
    // host/device interaction
    cudaMemcpy(h2d_spin_up, dev_spin2_up, dev_spinsize, cudaMemcpyDeviceToHost);
    cudaMemcpy(h2d_spin_dn, dev_spin2_dn, dev_spinsize, cudaMemcpyDeviceToHost);
    convert2double_spin(h2d_spin_up, g_chi_up_spinor_field[DUM_SOLVER+3]);
    convert2double_spin(h2d_spin_up, g_chi_dn_spinor_field[DUM_SOLVER+3]);
    
    // matrix multiplication
    Q_Qdagger_ND(g_chi_up_spinor_field[DUM_SOLVER+4], g_chi_dn_spinor_field[DUM_SOLVER+4],
                 g_chi_up_spinor_field[DUM_SOLVER+3], g_chi_dn_spinor_field[DUM_SOLVER+3] );
    
    // host/device interaction
    convert2REAL4_spin(g_chi_up_spinor_field[DUM_SOLVER+4], h2d_spin_up);
    convert2REAL4_spin(g_chi_dn_spinor_field[DUM_SOLVER+4], h2d_spin_dn);
    cudaMemcpy(dev_spin4_up, h2d_spin_up, dev_spinsize, cudaMemcpyHostToDevice);
    cudaMemcpy(dev_spin4_dn, h2d_spin_dn, dev_spinsize, cudaMemcpyHostToDevice);
    
    		// debug	// CUDA
  		#ifdef CUDA_DEBUG
  		  CUDA_CHECK("CUDA error in cg_eo_nd(). Applying the matrix on CPU failed.", "The matrix was applied on CPU.");
  		#endif
    
    #endif
    
    
    // pro = d*A*d
    pro_up = cublasSdot(N_floats, (float *) dev_spin2_up, 1, (float *) dev_spin4_up, 1);			// what about the imaginary and real parts !?
    pro_dn = cublasSdot(N_floats, (float *) dev_spin2_dn, 1, (float *) dev_spin4_dn, 1);			// does that really work ?
    pro = pro_up + pro_dn;
    
    // alpha = r(k-1)*r(k-1) / d*A*d
    alpha_cg = normsq / pro;
    
    // x(k+1) = x(k) + alpha*d(k)
    cublasSaxpy(N_floats, alpha_cg, (float *) dev_spin2_up, 1, (float *) dev_spin_up, 1);			// dev_spin_up = dev_spin_up + alpha * dev_spin2_up
    cublasSaxpy(N_floats, alpha_cg, (float *) dev_spin2_dn, 1, (float *) dev_spin_dn, 1);			// dev_spin_dn = dev_spin_dn + alpha * dev_spin2_dn
    
    // r(k+1) = r(k) - alpha*A*d(k)
    cublasSaxpy(N_floats, -1.0*alpha_cg, (float *) dev_spin4_up, 1, (float *) dev_spin1_up, 1);			// dev_spin1_up = dev_spin1_up - alpha*dev_spin4_up
    cublasSaxpy(N_floats, -1.0*alpha_cg, (float *) dev_spin4_dn, 1, (float *) dev_spin1_dn, 1);			// dev_spin1_dn = dev_spin1_dn - alpha*dev_spin4_dn
    
    // err = r(k) * r(k)
    err_up = cublasSdot(N_floats, (float *) dev_spin1_up, 1, (float *) dev_spin1_up, 1);			// err_up = (dev_spin1_up)^2
    err_dn = cublasSdot(N_floats, (float *) dev_spin1_dn, 1, (float *) dev_spin1_dn, 1);			// err_dn = (dev_spin1_dn)^2
    err = err_up + err_dn;
    
    		// debug	// CUBLAS core function
		#ifdef CUDA_DEBUG
		  CUBLAS_CORE_CHECK_NO_SUCCESS_MSG("CUBLAS error in cg_eo_nd(). CUBLAS function failed.");
		#endif
    
    		// debug
    		printf("inner iteration j = %i: err = %.8e\n", iteration, err);
    		
    		// debug	// is NaN ?
    		if isnan(err) {
    		  printf("Error in cg_eo_nd(). Inner residue is NaN.\n");
    		  exit(-1);
    		}
    
    // check wether precision is reached
    if ( (err <= eps_rel*squarenorm) || (err <= eps_abs) ) {
    
      // output
      dev_copy_spinor_field<<<griddim1, blockdim1>>>(dev_spin_up, P_up);					// P_up = dev_spin_up
      dev_copy_spinor_field<<<griddim1, blockdim1>>>(dev_spin_dn, P_dn);					// P_dn = dev_spin_dn
      
      g_sloppy_precision = 0;
      return(iteration+1);
    }
    
    // beta = r(k+1)*r(k+1) / r(k)*r(k)
    beta_cg = err / normsq;
    
    // d(k+1) = r(k+1) + beta*d(k)
    cublasSscal (N_floats, beta_cg, (float *) dev_spin2_up, 1);							// dev_spin2_up  = beta * dev_spin2_up
    cublasSaxpy (N_floats, 1.0 , (float *) dev_spin1_up, 1, (float *) dev_spin2_up, 1);				// dev_spin2_up += dev_spin1_up
    cublasSscal (N_floats, beta_cg, (float *) dev_spin2_dn, 1);							// dev_spin2_dn  = beta * dev_spin2_dn
    cublasSaxpy (N_floats, 1.0 , (float *) dev_spin1_dn, 1, (float *) dev_spin2_dn, 1);				// dev_spin2_dn += dev_spin1_dn
    
    normsq = err;	// for the next iteration
    
    
  }//LOOP
  
  
  // output
  dev_copy_spinor_field<<<griddim1, blockdim1>>>(dev_spin_up, P_up);						// P_up = dev_spin_up
  dev_copy_spinor_field<<<griddim1, blockdim1>>>(dev_spin_dn, P_dn);						// P_dn = dev_spin_dn
  
  g_sloppy_precision = 0;  
  
  // additional
  cudaFree(dev_spin_up);
  cudaFree(dev_spin_dn);
  cudaFree(dev_spin4_up);
  cudaFree(dev_spin4_dn);
  cudaFree(dev_spin5_up);
  cudaFree(dev_spin5_dn);
  
  return(-1);
  
   
}//alt_cg_eo_nd()
  
  
  
  
  
  

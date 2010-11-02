
// matrix_mpi_debug1() and matrix_mpi_debug1() replace matrix_multiplication32_mpi()
// matrix_mpi_debug3() replaces Q_Qdagger_ND()



extern "C" {
  //#ifdef HAVE_CONFIG_H
  //# include<config.h>
  //#endif
  //#include <stdlib.h>
  //#include <stdio.h>
  //#include "../global.h"
  //#include "../su3.h"
  #include "../../Hopping_Matrix.h"
  #include "../../phmc.h"
  #include "../../gamma.h"
  //#include "../linsolve.h"
  //#include "../linalg_eo.h"
  //#include "../Nondegenerate_Matrix.h"
}






//#define CHECK_HOPPING_MATRIX
//#define CHECK_IMUGAMMA5
//#define CHECK_GAMMA5
//#define CHECK_CUBLAS1
//#define CHECK_CUBLAS2
//#define CHECK_CUBLAS3
//#define CHECK_COPY
//#define CHECK_MAXEV









// replaces matrix_multiplication32_mpi() for debugging !!

void matrix_mpi_debug1 (dev_spinor * spinout_up, dev_spinor * spinout_dn,
                        dev_spinor * spinin_up , dev_spinor * spinin_dn ,
                        int gridsize1, int blocksize1, int gridsize2, int blocksize2,
                        int gridsize3, int blocksize3, int gridsize4, int blocksize4) {
  
  
  // we will use the auxiliary fields  dev_spin_eo{1,2}_up/dn  for working on and buffering
  // and set  dev_spin_eo2_up/dn  equal  spinout_up/dn
  // spinin_up/dn  have to remain unchanged !!
  // spinout_up/dn  can be freely used
  
  
  
  
  /////////////////////
  // LOCAL VARIABLES //
  /////////////////////
  
  int N_sites  =    VOLUME/2;		// #lattice sites
  int N_floats = 24*VOLUME/2;		// #floats
  
  // additional for the debugging purposes
  size_t dev_spinsize = 6*(VOLUME+RAND)/2*sizeof(dev_spinor);
  
  
  
  
  ////////////////////////////////////
  //   MATCHING with Q_Qdagger_ND   //
  ////////////////////////////////////
  //                                //
  // _strange = _up                 //
  // _charm   = _dn                 //
  //                                //
  // DUM_MATRIX   = dev_spin_eo1_up //
  // DUM_MATRIX+1 = dev_spin_eo1_dn //
  //                                //
  // DUM_MATRIX+2 = dev_spin_eo2_up //
  // DUM_MATRIX+3 = dev_spin_eo2_dn //
  //                                //
  ////////////////////////////////////
  
  // can savely use the following spinors on host:	g_spinor_field[DUM_MATRIX{  , +1, ... , +7}]
  
  
  
  
  ///////////////////////////////////
  // INITIALIZATIONS & ASSIGNMENTS //	// have to use (one) other auxiliary field(s) than the calling function dev_cg_eo_nd
  ///////////////////////////////////
  
  dev_spin_eo2_up = spinout_up;		// need no memory allocated
  dev_spin_eo2_dn = spinout_dn;
  												///////////// THEORY ////////////////////////////////////////////////////////////////
  												//                                                                                 //
  												//  (Q_tilde) = gamma5 * ((M_oo) - (M_oe)(Mee^-1)(M_eo))                           //
  												//  (Q_tilde)(Q_tilde_dagger) * (up,dn) = (Q_tilde) * (b,a)                        //
  ///////////////										//                                                    (a,b) = (Q_tilde) * (dn,up)  //
  // MAIN BODY //										//                                                                                 //
  ///////////////										/////////////////////////////////////////////////////////////////////////////////////
  
  
  double nrm = 1.0 / (1.0 + g_mubar*g_mubar - g_epsbar*g_epsbar);
  
  
  if (g_proc_id == 0)
    printf("This is matrix_mpi_debug1(). ");
  
  
  spinor * l_strange = (spinor *) malloc(2*dev_spinsize);
  spinor * l_charm   = (spinor *) malloc(2*dev_spinsize);
  spinor * k_strange = (spinor *) malloc(2*dev_spinsize);
  spinor * k_charm   = (spinor *) malloc(2*dev_spinsize);
  
  
  ///////////////////////////////////////							/////////////////////////////////
  //        Q_tilde_dagger(2x2)        //							// (a,b) = (Q_tilde) * (dn,up) //
  ///////////////////////////////////////							/////////////////////////////////
  
  
  ////////////////////////////
  // (M_oe) (Mee^-1) (M_eo) //
  ////////////////////////////
  
  #ifndef CHECK_HOPPING_MATRIX
  			// xchange
  			#ifndef HOPPING_DEBUG
  			  xchange_field_wrapper(spinin_dn, 0);
  			#endif
  	  #ifndef HOPPING_DEBUG
  	    dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, spinin_dn, dev_spin_eo1_up, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0);	// dev_spin_eo1_up = (M_eo) * spinin_dn
  	  #else
  	    Hopping_Matrix_wrapper(0, dev_spin_eo1_up, spinin_dn);
  	  #endif
  			// xchange
  			#ifndef HOPPING_DEBUG
  			  xchange_field_wrapper(spinin_up, 0);
  			#endif
  	  #ifndef HOPPING_DEBUG
  	    dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, spinin_up, dev_spin_eo1_dn, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0);	// dev_spin_eo1_dn = (M_eo) * spinin_up
  	  #else
  	    Hopping_Matrix_wrapper(0, dev_spin_eo1_dn, spinin_up);
  	  #endif
  #else
    	to_host(k_charm, spinin_dn, h2d_spin_up, dev_spinsize);
    	Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX], k_charm);				// g_spinor_field[DUM_MATRIX]   = (M_eo) * k_charm
    	to_device(dev_spin_eo1_up, g_spinor_field[DUM_MATRIX], h2d_spin_up, dev_spinsize);
    	
    	to_host(k_strange, spinin_up, h2d_spin_dn, dev_spinsize);
    	Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX+1], k_strange);				// g_spinor_field[DUM_MATRIX+1] = (M_eo) * k_strange
    	to_device(dev_spin_eo1_dn, g_spinor_field[DUM_MATRIX+1], h2d_spin_dn, dev_spinsize);
  #endif
  
  
  
  
  #ifndef CHECK_IMUGAMMA5
  	dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo1_up, dev_spin_eo2_up, -1.0);		// dev_spin_eo2_up  =  (1 - imubar) * dev_spin_eo1_up  =  (1 - imubar)*(M_eo) * spinin_dn
  	dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo1_dn, dev_spin_eo2_dn, +1.0);		// dev_spin_eo2_dn  =  (1 + imubar) * dev_spin_eo1_dn  =  (1 + imubar)*(M_eo) * spinin_up
  #else
    	to_host(g_spinor_field[DUM_MATRIX], dev_spin_eo1_up, h2d_spin_up, dev_spinsize);
    	mul_one_minus_imubar(g_spinor_field[DUM_MATRIX+2], g_spinor_field[DUM_MATRIX]);		// g_spinor_field[DUM_MATRIX+2] = (1 - imubar)*(M_eo) * g_spinor_field[DUM_MATRIX]
    	to_device(dev_spin_eo2_up, g_spinor_field[DUM_MATRIX+2], h2d_spin_up, dev_spinsize);
    	
    	to_host(g_spinor_field[DUM_MATRIX+1], dev_spin_eo1_dn, h2d_spin_dn, dev_spinsize);
    	mul_one_plus_imubar (g_spinor_field[DUM_MATRIX+3], g_spinor_field[DUM_MATRIX+1]);	// g_spinor_field[DUM_MATRIX+3] = (1 + imubar)*(M_eo) * g_spinor_field[DUM_MATRIX+1]
    	to_device(dev_spin_eo2_dn, g_spinor_field[DUM_MATRIX+3], h2d_spin_dn, dev_spinsize);
  #endif
  
  
  
  
  #ifndef CHECK_CUBLAS1
  	cublasSaxpy (N_floats, g_epsbar, (float *) dev_spin_eo1_dn, 1, (float *) dev_spin_eo2_up, 1);
  										// dev_spin_eo2_up  =  dev_spin_eo2_up  +  epsbar * dev_spin_eo1_dn  =  (1 - imubar)*(M_eo) * spinin_dn  +  epsbar * (M_eo) * spinin_up
  	cublasSaxpy (N_floats, g_epsbar, (float *) dev_spin_eo1_up, 1, (float *) dev_spin_eo2_dn, 1);
  										// dev_spin_eo2_dn  =  dev_spin_eo2_dn  +  epsbar * dev_spin_eo1_up  =  (1 + imubar)*(M_eo) * spinin_up  +  epsbar * (M_eo) * spinin_dn
  
  	cublasSscal (N_floats, nrm, (float *) dev_spin_eo2_up, 1);			// dev_spin_eo2_up  =  nrm * dev_spin_eo2_up  =  nrm*(1-imubar)*(M_eo) * spinin_dn  +  nrm*epsbar*(M_eo) * spinin_up
  
  	cublasSscal (N_floats, nrm, (float *) dev_spin_eo2_dn, 1);			// dev_spin_eo2_dn  =  nrm * dev_spin_eo2_dn  =  nrm*(1+imubar)*(M_eo) * spinin_up  +  nrm*epsbar*(M_eo) * spinin_dn
  #else
    	to_host(g_spinor_field[DUM_MATRIX+1], dev_spin_eo1_dn, h2d_spin_up, dev_spinsize);
    	to_host(g_spinor_field[DUM_MATRIX+2], dev_spin_eo2_up, h2d_spin_up, dev_spinsize);
    	assign_add_mul_r(g_spinor_field[DUM_MATRIX+2], g_spinor_field[DUM_MATRIX+1], g_epsbar, VOLUME/2);
  												// g_spinor_field[DUM_MATRIX+2] = g_spinor_field[DUM_MATRIX+2]     +  epsbar*g_spinor_field[DUM_MATRIX+1]
    	mul_r(g_spinor_field[DUM_MATRIX+2], nrm, g_spinor_field[DUM_MATRIX+2], VOLUME/2);	// g_spinor_field[DUM_MATRIX+2] = nrm * g_spinor_field[DUM_MATRIX+2]
    	to_device(dev_spin_eo2_up, g_spinor_field[DUM_MATRIX+2], h2d_spin_up, dev_spinsize);
    	
    	
    	to_host(g_spinor_field[DUM_MATRIX], dev_spin_eo1_up, h2d_spin_dn, dev_spinsize);
    	to_host(g_spinor_field[DUM_MATRIX+3], dev_spin_eo2_dn, h2d_spin_dn, dev_spinsize);
    	assign_add_mul_r(g_spinor_field[DUM_MATRIX+3], g_spinor_field[DUM_MATRIX], g_epsbar, VOLUME/2);
  												// g_spinor_field[DUM_MATRIX+3] = g_spinor_field[DUM_MATRIX+3]     +  epsbar*g_spinor_field[DUM_MATRIX]
  	 mul_r(g_spinor_field[DUM_MATRIX+3], nrm, g_spinor_field[DUM_MATRIX+3], VOLUME/2);	// g_spinor_field[DUM_MATRIX+3] = nrm * g_spinor_field[DUM_MATRIX+3]
  	 to_device(dev_spin_eo2_dn, g_spinor_field[DUM_MATRIX+3] , h2d_spin_dn, dev_spinsize);
  #endif
  
  
  
  
  #ifndef CHECK_HOPPING_MATRIX
  			// xchange
  			#ifndef HOPPING_DEBUG
  			  xchange_field_wrapper(dev_spin_eo2_up, 1);
  			#endif
  	  #ifndef HOPPING_DEBUG
  	    dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo2_up, dev_spin_eo1_up, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1);
  	  #else
  	    Hopping_Matrix_wrapper(1, dev_spin_eo1_up, dev_spin_eo2_up);
  	  #endif							
  			// xchange
  			#ifndef HOPPING_DEBUG
  			  xchange_field_wrapper(dev_spin_eo2_dn, 1);
  			#endif
  	  #ifndef HOPPING_DEBUG
  	    dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo2_dn, dev_spin_eo1_dn, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1);
  	  #else
  	    Hopping_Matrix_wrapper(1, dev_spin_eo1_dn, dev_spin_eo2_dn);
  	  #endif
  #else
    	to_host(g_spinor_field[DUM_MATRIX+2], dev_spin_eo2_up, h2d_spin_up, dev_spinsize);
    	Hopping_Matrix(OE, g_spinor_field[DUM_MATRIX]  , g_spinor_field[DUM_MATRIX+2]);		// g_spinor_field[DUM_MATRIX]   = (M_oe) * g_spinor_field[DUM_MATRIX+2]
    	to_device(dev_spin_eo1_up, g_spinor_field[DUM_MATRIX], h2d_spin_up, dev_spinsize);
    	
    	to_host(g_spinor_field[DUM_MATRIX+3], dev_spin_eo2_dn, h2d_spin_dn, dev_spinsize);
    	Hopping_Matrix(OE, g_spinor_field[DUM_MATRIX+1], g_spinor_field[DUM_MATRIX+3]);		// g_spinor_field[DUM_MATRIX+1] = (M_oe) * g_spinor_field[DUM_MATRIX+3]
    	to_device(dev_spin_eo1_dn, g_spinor_field[DUM_MATRIX+1], h2d_spin_dn, dev_spinsize);
  #endif
  
  
  
  
  ////////////
  // (M_oo) //
  ////////////
  
  
  #ifndef CHECK_IMUGAMMA5
  	dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(spinin_dn, dev_spin_eo2_up, +1.0);			// dev_spin_eo2_up = (1 + imubar) * spinin_dn
  	dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(spinin_up, dev_spin_eo2_dn, -1.0);			// dev_spin_eo2_dn = (1 - imubar) * spinin_up
  #else
    	to_host(k_charm, spinin_dn, h2d_spin_up, dev_spinsize);
    	mul_one_plus_imubar (g_spinor_field[DUM_MATRIX+2], k_charm);				// g_spinor_field[DUM_MATRIX+2] = (1 + imubar) * k_charm
    	to_device(dev_spin_eo2_up, g_spinor_field[DUM_MATRIX+2], h2d_spin_up, dev_spinsize);
    	
    	to_host(k_strange, spinin_up, h2d_spin_dn, dev_spinsize);
    	mul_one_minus_imubar(g_spinor_field[DUM_MATRIX+3], k_strange);				// g_spinor_field[DUM_MATRIX+3] = (1 - imubar) * k_strange
    	to_device(dev_spin_eo2_dn, g_spinor_field[DUM_MATRIX+3], h2d_spin_dn, dev_spinsize);
  #endif
  
  
  
  
  #ifndef CHECK_CUBLAS2													// remember: this is (M_oo) * (spinin_dn, spinin_up):
  	cublasSaxpy (N_floats, -g_epsbar, (float *) spinin_up, 1, (float *) dev_spin_eo2_up, 1);
  												// dev_spin_eo2_up  =  dev_spin_eo2_up - epsbar*spinin_up  =  (1+imubar)*spinin_dn - epsbar*spinin_up
  	cublasSaxpy (N_floats, -g_epsbar, (float *) spinin_dn, 1, (float *) dev_spin_eo2_dn, 1);
  												// dev_spin_eo2_dn  =  dev_spin_eo2_dn - epsbar*spinin_dn  =  (1-imubar)*spinin_up - epsbar*spinin_dn
  #else
    	to_host(k_strange, spinin_up, h2d_spin_up, dev_spinsize);
    	to_host(g_spinor_field[DUM_MATRIX+2], dev_spin_eo2_up, h2d_spin_up, dev_spinsize);
    	assign_add_mul_r(g_spinor_field[DUM_MATRIX+2], k_strange, -g_epsbar, VOLUME/2);	// g_spinor_field[DUM_MATRIX+2] = g_spinor_field[DUM_MATRIX+2] - epsbar*k_strange
    	to_device(dev_spin_eo2_up, g_spinor_field[DUM_MATRIX+2], h2d_spin_up, dev_spinsize);
    	
    	to_host(k_charm, spinin_dn, h2d_spin_dn, dev_spinsize);
    	to_host(g_spinor_field[DUM_MATRIX+3], dev_spin_eo2_dn, h2d_spin_dn, dev_spinsize);
    	assign_add_mul_r(g_spinor_field[DUM_MATRIX+3], k_charm  , -g_epsbar, VOLUME/2);	// g_spinor_field[DUM_MATRIX+3] = g_spinor_field[DUM_MATRIX+3] - epsbar*k_charm
    	to_device(dev_spin_eo2_dn, g_spinor_field[DUM_MATRIX+3], h2d_spin_dn, dev_spinsize);
  #endif
  
  
  
  
  ///////////////////////////////////////
  // (M_oo)  -  (M_oe) (Mee^-1) (M_eo) //
  ///////////////////////////////////////
  
  #ifndef CHECK_CUBLAS3													// this is ((M_oo) - (M_oe)(Mee^-1)(M_eo)) * (spinin_dn, spinin_up):
  	cublasSaxpy (N_floats, -1.0, (float *) dev_spin_eo1_up, 1, (float *) dev_spin_eo2_up, 1);
  							// dev_spin_eo2_up  =  dev_spin_eo2_up  -  dev_spin_eo1_up  =                      (1+imubar) * spinin_dn  -                    epsbar * spinin_up
  							//                                                             - (M_oe)*nrm*(1-imubar)*(M_eo) * spinin_dn  -  (M_oe)*nrm*epsbar*(M_eo) * spinin_up
  	cublasSaxpy (N_floats, -1.0, (float *) dev_spin_eo1_dn, 1, (float *) dev_spin_eo2_dn, 1);
  							// dev_spin_eo2_dn  =  dev_spin_eo2_dn  -  dev_spin_eo1_dn  =                      (1-imubar) * spinin_up  -                    epsbar * spinin_dn
  							//                                                             - (M_oe)*nrm*(1+imubar)*(M_eo) * spinin_up  -  (M_oe)*nrm*epsbar*(M_eo) * spinin_dn
  #else
    	to_host(g_spinor_field[DUM_MATRIX], dev_spin_eo1_up, h2d_spin_up, dev_spinsize);
    	to_host(g_spinor_field[DUM_MATRIX+2], dev_spin_eo2_up, h2d_spin_up, dev_spinsize);
    	diff(g_spinor_field[DUM_MATRIX+4], g_spinor_field[DUM_MATRIX+2], g_spinor_field[DUM_MATRIX]  , VOLUME/2);
    												// g_spinor_field[DUM_MATRIX+4] = g_spinor_field[DUM_MATRIX+2] - g_spinor_field[DUM_MATRIX]
    	to_device(dev_spin_eo2_up, g_spinor_field[DUM_MATRIX+4], h2d_spin_up, dev_spinsize);
    	
    	
    	to_host(g_spinor_field[DUM_MATRIX+1], dev_spin_eo1_dn, h2d_spin_dn, dev_spinsize);
    	to_host(g_spinor_field[DUM_MATRIX+3], dev_spin_eo2_dn, h2d_spin_dn, dev_spinsize);
    	diff(g_spinor_field[DUM_MATRIX+5], g_spinor_field[DUM_MATRIX+3], g_spinor_field[DUM_MATRIX+1], VOLUME/2);
    												// g_spinor_field[DUM_MATRIX+5] = g_spinor_field[DUM_MATRIX+3] - g_spinor_field[DUM_MATRIX+1]
    	to_device(dev_spin_eo2_dn, g_spinor_field[DUM_MATRIX+5], h2d_spin_dn, dev_spinsize);
  #endif
  
  
  
  
  ////////////
  // gamma5 //
  ////////////
  
  #ifndef CHECK_GAMMA5
  	dev_gamma5<<<gridsize3, blocksize3>>>(dev_spin_eo2_up, dev_spin_eo2_up);					// dev_spin_eo2_up = gamma5 * dev_spin_eo2_up 
  	dev_gamma5<<<gridsize3, blocksize3>>>(dev_spin_eo2_dn, dev_spin_eo2_dn);					// dev_spin_eo2_dn = gamma5 * dev_spin_eo2_dn
  #else
    	to_host(g_spinor_field[DUM_MATRIX+4], dev_spin_eo2_up, h2d_spin_up, dev_spinsize);
    	gamma5(g_spinor_field[DUM_MATRIX+2], g_spinor_field[DUM_MATRIX+4], VOLUME/2);		// g_spinor_field[DUM_MATRIX+2] = gamma5 * g_spinor_field[DUM_MATRIX+4]
    	to_device(dev_spin_eo2_up, g_spinor_field[DUM_MATRIX+2], h2d_spin_up, dev_spinsize);
    	
    	to_host(g_spinor_field[DUM_MATRIX+5], dev_spin_eo2_dn, h2d_spin_dn, dev_spinsize);
    	gamma5(g_spinor_field[DUM_MATRIX+3], g_spinor_field[DUM_MATRIX+5], VOLUME/2);		// g_spinor_field[DUM_MATRIX+3] = gamma5 * g_spinor_field[DUM_MATRIX+5]
    	to_device(dev_spin_eo2_dn, g_spinor_field[DUM_MATRIX+3], h2d_spin_dn, dev_spinsize);
  #endif	
  
  
  
  
  
  
  ////////////////////
  // (a,b) -> (b,a) //
  ////////////////////
  
  #ifndef CHECK_COPY
  	dev_copy_spinor_field<<<gridsize4, blocksize4>>>(dev_spin_eo2_dn, dev_spin_eo3_up);			// dev_spin_eo3_up = dev_spin_eo2_dn
  	dev_copy_spinor_field<<<gridsize4, blocksize4>>>(dev_spin_eo2_up, dev_spin_eo3_dn);			// dev_spin_eo3_dn = dev_spin_eo2_up
  #else
    	to_host(g_spinor_field[DUM_MATRIX+2], dev_spin_eo2_dn, h2d_spin_up, dev_spinsize);
    	assign(g_spinor_field[DUM_MATRIX+6], g_spinor_field[DUM_MATRIX+2], VOLUME/2);		// g_spinor_field[DUM_MATRIX+6] = g_spinor_field[DUM_MATRIX+2]
    	to_device(dev_spin_eo3_up, g_spinor_field[DUM_MATRIX+6], h2d_spin_up, dev_spinsize);
    	
    	to_host(g_spinor_field[DUM_MATRIX+3], dev_spin_eo2_up, h2d_spin_dn, dev_spinsize);
    	assign(g_spinor_field[DUM_MATRIX+7], g_spinor_field[DUM_MATRIX+3], VOLUME/2);		// g_spinor_field[DUM_MATRIX+7] = g_spinor_field[DUM_MATRIX+3]
    	to_device(dev_spin_eo3_dn, g_spinor_field[DUM_MATRIX+7], h2d_spin_dn, dev_spinsize);
  #endif
  
  
  
  
  
  
  ///////////////////////////////////								///////////////////////
  //        Q_tilde(2x2)           //								// (Q_tilde) * (b,a) //
  ///////////////////////////////////								///////////////////////
  
  
  ////////////////////////////
  // (M_oe) (Mee^-1) (M_eo) //
  ////////////////////////////
  
  #ifndef CHECK_HOPPING_MATRIX
  			// xchange
  			#ifndef HOPPING_DEBUG
  			  xchange_field_wrapper(dev_spin_eo3_up, 0);
  			#endif
  	  #ifndef HOPPING_DEBUG
  	    dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo3_up, dev_spin_eo1_up, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0);	// dev_spin_eo1_up = (M_eo) * dev_spin_eo3_up
  	  #else
  	    Hopping_Matrix_wrapper(0, dev_spin_eo1_up, dev_spin_eo3_up);
  	  #endif
  			// xchange
  			#ifndef HOPPING_DEBUG
  			  xchange_field_wrapper(dev_spin_eo3_dn, 0);
  			#endif
  	  #ifndef HOPPING_DEBUG
  	    dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo3_dn, dev_spin_eo1_dn, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0);	// dev_spin_eo1_dn = (M_eo) * dev_spin_eo3_dn
  	  #else
  	    Hopping_Matrix_wrapper(0, dev_spin_eo1_dn, dev_spin_eo3_dn);
  	  #endif
  #else
    	to_host(g_spinor_field[DUM_MATRIX+7], dev_spin_eo3_up, h2d_spin_up, dev_spinsize);
    	Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX]  , g_spinor_field[DUM_MATRIX+7]);	// g_spinor_field[DUM_MATRIX]   = (M_eo) * g_spinor_field[DUM_MATRIX+7]
    	to_device(dev_spin_eo1_up, g_spinor_field[DUM_MATRIX], h2d_spin_up, dev_spinsize);
    	
    	to_host(g_spinor_field[DUM_MATRIX+6], dev_spin_eo3_dn, h2d_spin_dn, dev_spinsize);
    	Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX+1], g_spinor_field[DUM_MATRIX+6]);	// g_spinor_field[DUM_MATRIX+1] = (M_eo) * g_spinor_field[DUM_MATRIX+6]
    	to_device(dev_spin_eo1_dn, g_spinor_field[DUM_MATRIX+1], h2d_spin_dn, dev_spinsize);
  #endif
  
  
  
  
  #ifndef CHECK_IMUGAMMA5
  	dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo1_up, dev_spin_eo2_up, -1.0);			// dev_spin_eo2_up  =  (1 - imubar) * dev_spin_eo1_up  =  (1 - imubar)*(M_eo) * dev_spin_eo3_up
  	dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo1_dn, dev_spin_eo2_dn, +1.0);			// dev_spin_eo2_dn  =  (1 + imubar) * dev_spin_eo1_dn  =  (1 + imubar)*(M_eo) * dev_spin_eo3_dn
  #else
    	to_host(g_spinor_field[DUM_MATRIX], dev_spin_eo1_up, h2d_spin_up, dev_spinsize);
    	mul_one_minus_imubar(g_spinor_field[DUM_MATRIX+2], g_spinor_field[DUM_MATRIX]);			// g_spinor_field[DUM_MATRIX+2] = (1 - imubar) * g_spinor_field[DUM_MATRIX]
    	to_device(dev_spin_eo2_up, g_spinor_field[DUM_MATRIX+2], h2d_spin_up, dev_spinsize);
    	
    	to_host(g_spinor_field[DUM_MATRIX+1], dev_spin_eo1_dn, h2d_spin_dn, dev_spinsize);
    	mul_one_plus_imubar (g_spinor_field[DUM_MATRIX+3], g_spinor_field[DUM_MATRIX+1]);		// g_spinor_field[DUM_MATRIX+3] = (1 + imubar) * g_spinor_field[DUM_MATRIX+1]
    	to_device(dev_spin_eo2_dn, g_spinor_field[DUM_MATRIX+3], h2d_spin_dn, dev_spinsize);
  #endif
  
  
  
  
  #ifndef CHECK_CUBLAS1
  	cublasSaxpy (N_floats, g_epsbar, (float *) dev_spin_eo1_dn, 1, (float *) dev_spin_eo2_up, 1);
  										// dev_spin_eo2_up  =  dev_spin_eo2_up  +  epsbar * dev_spin_eo1_dn  =  (1 - imubar)*(M_eo) * dev_spin_eo3_up  +  epsbar * (M_eo) * dev_spin_eo3_dn
  	cublasSaxpy (N_floats, g_epsbar, (float *) dev_spin_eo1_up, 1, (float *) dev_spin_eo2_dn, 1);
  										// dev_spin_eo2_dn  =  dev_spin_eo2_dn  +  epsbar * dev_spin_eo1_up  =  (1 + imubar)*(M_eo) * dev_spin_eo3_dn  +  epsbar * (M_eo) * dev_spin_eo3_up
  
  	cublasSscal (N_floats, nrm, (float *) dev_spin_eo2_up, 1);			// dev_spin_eo2_up  =  nrm * dev_spin_eo2_up  =  nrm*(1-imubar)*(M_eo) * dev_spin_eo3_up  +  nrm*epsbar*(M_eo) * dev_spin_eo3_dn
  
  	cublasSscal (N_floats, nrm, (float *) dev_spin_eo2_dn, 1);			// dev_spin_eo2_dn  =  nrm * dev_spin_eo2_dn  =  nrm*(1+imubar)*(M_eo) * dev_spin_eo3_dn  +  nrm*epsbar*(M_eo) * dev_spin_eo3_up
  #else
    	to_host(g_spinor_field[DUM_MATRIX+1], dev_spin_eo1_dn, h2d_spin_up, dev_spinsize);
    	to_host(g_spinor_field[DUM_MATRIX+2], dev_spin_eo2_up, h2d_spin_up, dev_spinsize);
    	assign_add_mul_r(g_spinor_field[DUM_MATRIX+2], g_spinor_field[DUM_MATRIX+1], g_epsbar, VOLUME/2);
  												// g_spinor_field[DUM_MATRIX+2] = g_spinor_field[DUM_MATRIX+2]   +  epsbar * g_spinor_field[DUM_MATRIX+1]
    	mul_r(g_spinor_field[DUM_MATRIX+2], nrm, g_spinor_field[DUM_MATRIX+2], VOLUME/2);	// g_spinor_field[DUM_MATRIX+2] = nrm * g_spinor_field[DUM_MATRIX+2]
    	to_device(dev_spin_eo2_up, g_spinor_field[DUM_MATRIX+2], h2d_spin_up, dev_spinsize);
    	
    	
    	to_host(g_spinor_field[DUM_MATRIX], dev_spin_eo1_up, h2d_spin_dn, dev_spinsize);
    	to_host(g_spinor_field[DUM_MATRIX+3], dev_spin_eo2_dn, h2d_spin_dn, dev_spinsize);
    	assign_add_mul_r(g_spinor_field[DUM_MATRIX+3], g_spinor_field[DUM_MATRIX]  , g_epsbar, VOLUME/2);
  												// g_spinor_field[DUM_MATRIX+3] = g_spinor_field[DUM_MATRIX+3]   +  epsbar * g_spinor_field[DUM_MATRIX]
  	mul_r(g_spinor_field[DUM_MATRIX+3], nrm, g_spinor_field[DUM_MATRIX+3], VOLUME/2);	// g_spinor_field[DUM_MATRIX+3] = nrm * g_spinor_field[DUM_MATRIX+3]
    	to_device(dev_spin_eo2_dn, g_spinor_field[DUM_MATRIX+3] , h2d_spin_dn, dev_spinsize);
  #endif
  
  
  
  
  #ifndef CHECK_HOPPING_MATRIX
  			// xchange
  			#ifndef HOPPING_DEBUG
  			  xchange_field_wrapper(dev_spin_eo2_up, 1);
  			#endif
  	  #ifndef HOPPING_DEBUG
  	    dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo2_up, dev_spin_eo1_up, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1);
  	  #else
  	    Hopping_Matrix_wrapper(1, dev_spin_eo1_up, dev_spin_eo2_up);
  	  #endif
  			// xchange
  			#ifndef HOPPING_DEBUG
  			  xchange_field_wrapper(dev_spin_eo2_dn, 1);
  			#endif
  	  #ifndef HOPPING_DEBUG
  	    dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo2_dn, dev_spin_eo1_dn, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1);
  	  #else
  	    Hopping_Matrix_wrapper(1, dev_spin_eo1_dn, dev_spin_eo2_dn);
  	  #endif
  #else
    	to_host(g_spinor_field[DUM_MATRIX+2], dev_spin_eo2_up, h2d_spin_up, dev_spinsize);
    	Hopping_Matrix(OE, l_strange, g_spinor_field[DUM_MATRIX+2]);				// l_strange = (M_oe) * g_spinor_field[DUM_MATRIX+2]
    	to_device(dev_spin_eo1_up, l_strange, h2d_spin_up, dev_spinsize);
    	
    	to_host(g_spinor_field[DUM_MATRIX+3], dev_spin_eo2_dn, h2d_spin_dn, dev_spinsize);
    	Hopping_Matrix(OE, l_charm  , g_spinor_field[DUM_MATRIX+3]);				// l_charm   = (M_oe) * g_spinor_field[DUM_MATRIX+3]
    	to_device(dev_spin_eo1_dn, l_charm, h2d_spin_dn, dev_spinsize);
  #endif
  
  
  
  
  ////////////
  // (M_oo) //
  ////////////
  
  #ifndef CHECK_IMUGAMMA5
  	dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo3_up, dev_spin_eo2_up, +1.0);				// dev_spin_eo2_up = (1 + imubar) * dev_spin_eo3_up
  	dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo3_dn, dev_spin_eo2_dn, -1.0);				// dev_spin_eo2_dn = (1 - imubar) * dev_spin_eo3_dn
  #else
    	to_host(g_spinor_field[DUM_MATRIX+7], dev_spin_eo3_up, h2d_spin_up, dev_spinsize);
    	mul_one_plus_imubar (g_spinor_field[DUM_MATRIX]  , g_spinor_field[DUM_MATRIX+7]);		// g_spinor_field[DUM_MATRIX]   = (1 + imubar) * g_spinor_field[DUM_MATRIX+7]
    	to_device(dev_spin_eo2_up, g_spinor_field[DUM_MATRIX], h2d_spin_up, dev_spinsize);
    	
    	to_host(g_spinor_field[DUM_MATRIX+6], dev_spin_eo3_dn, h2d_spin_dn, dev_spinsize);
    	mul_one_minus_imubar(g_spinor_field[DUM_MATRIX+1], g_spinor_field[DUM_MATRIX+6]);		// g_spinor_field[DUM_MATRIX+1] = (1 - imubar) * g_spinor_field[DUM_MATRIX+6]
    	to_device(dev_spin_eo2_dn, g_spinor_field[DUM_MATRIX+1], h2d_spin_dn, dev_spinsize);  
  #endif
  
  
  
  
  #ifndef CHECK_CUBLAS2											// remember: this is (M_oo) * (dev_spin_eo3_up, dev_spin_eo3_dn):
  	cublasSaxpy (N_floats, -g_epsbar, (float *) dev_spin_eo3_dn, 1, (float *) dev_spin_eo2_up, 1);
  												// dev_spin_eo2_up  =  dev_spin_eo2_up - epsbar*dev_spin_eo3_dn  =  (1+imubar)*dev_spin_eo3_up - epsbar*dev_spin_eo3_dn
  	cublasSaxpy (N_floats, -g_epsbar, (float *) dev_spin_eo3_up, 1, (float *) dev_spin_eo2_dn, 1);
  												// dev_spin_eo2_dn  =  dev_spin_eo2_dn - epsbar*dev_spin_eo3_up  =  (1-imubar)*dev_spin_eo3_dn - epsbar*dev_spin_eo3_up
  #else
    	to_host(g_spinor_field[DUM_MATRIX+6], dev_spin_eo3_dn, h2d_spin_up, dev_spinsize);
    	to_host(g_spinor_field[DUM_MATRIX], dev_spin_eo2_up, h2d_spin_up, dev_spinsize);
    	assign_add_mul_r(g_spinor_field[DUM_MATRIX]  , g_spinor_field[DUM_MATRIX+6], -g_epsbar, VOLUME/2);
  												// g_spinor_field[DUM_MATRIX]   = g_spinor_field[DUM_MATRIX]    -  epsbar * g_spinor_field[DUM_MATRIX+6]
    	to_device(dev_spin_eo2_up, g_spinor_field[DUM_MATRIX], h2d_spin_up, dev_spinsize);
    	
    	
    	to_host(g_spinor_field[DUM_MATRIX+7], dev_spin_eo3_up, h2d_spin_dn, dev_spinsize);
    	to_host(g_spinor_field[DUM_MATRIX+1], dev_spin_eo2_dn, h2d_spin_dn, dev_spinsize);
    	assign_add_mul_r(g_spinor_field[DUM_MATRIX+1], g_spinor_field[DUM_MATRIX+7], -g_epsbar, VOLUME/2);
  												// g_spinor_field[DUM_MATRIX+1] = g_spinor_field[DUM_MATRIX+1]  -  epsbar * g_spinor_field[DUM_MATRIX+7]
    	to_device(dev_spin_eo2_dn, g_spinor_field[DUM_MATRIX+1], h2d_spin_dn, dev_spinsize);
  #endif
  
  
  
  
  ///////////////////////////////////////
  // (M_oo)  -  (M_oe) (Mee^-1) (M_eo) //
  ///////////////////////////////////////
  
  #ifndef CHECK_CUBLAS3											// this is ( (M_oo) - (M_oe) (Mee^-1) (M_eo) ) * (dev_spin_eo3_up, dev_spin_eo3_dn)
  	cublasSaxpy (N_floats, -1.0, (float *) dev_spin_eo1_up, 1, (float *) dev_spin_eo2_up, 1);
  							// dev_spin_eo2_up  =  dev_spin_eo2_up  -  dev_spin_eo1_up  =                      (1+imubar) * dev_spin_eo3_up  -                    epsbar * dev_spin_eo3_dn
  							//                                                             - (M_oe)*nrm*(1-imubar)*(M_eo) * dev_spin_eo3_up  -  (M_oe)*nrm*epsbar*(M_eo) * dev_spin_eo3_dn
  	cublasSaxpy (N_floats, -1.0, (float *) dev_spin_eo1_dn, 1, (float *) dev_spin_eo2_dn, 1);
  							// dev_spin_eo2_dn  =  dev_spin_eo2_dn  -  dev_spin_eo1_dn  =                      (1-imubar) * dev_spin_eo3_dn  -                    epsbar * dev_spin_eo3_up
  							//                                                             - (M_oe)*nrm*(1+imubar)*(M_eo) * dev_spin_eo3_dn  -  (M_oe)*nrm*epsbar*(M_eo) * dev_spin_eo3_up
  #else
    	to_host(l_strange, dev_spin_eo1_up, h2d_spin_up, dev_spinsize);
    	to_host(g_spinor_field[DUM_MATRIX], dev_spin_eo2_up, h2d_spin_up, dev_spinsize);
    	diff(l_strange, g_spinor_field[DUM_MATRIX], l_strange, VOLUME/2);			// l_strange = g_spinor_field[DUM_MATRIX]   - l_strange
    	to_device(dev_spin_eo2_up, l_strange, h2d_spin_up, dev_spinsize);
    	
    	to_host(l_charm, dev_spin_eo1_dn, h2d_spin_dn, dev_spinsize);
    	to_host(g_spinor_field[DUM_MATRIX+1], dev_spin_eo2_dn, h2d_spin_dn, dev_spinsize);
    	diff(l_charm, g_spinor_field[DUM_MATRIX+1], l_charm  , VOLUME/2);			// l_charm   = g_spinor_field[DUM_MATRIX+1] - l_charm
    	to_device(dev_spin_eo2_dn, l_charm, h2d_spin_dn, dev_spinsize);
  #endif
  
  
  
  ////////////
  // gamma5 //
  ////////////
  
  #ifndef CHECK_GAMMA5
  	dev_gamma5<<<gridsize3, blocksize3>>>(dev_spin_eo2_up, dev_spin_eo2_up);					// dev_spin_eo2_up = gamma5 * dev_spin_eo2_up 
  	dev_gamma5<<<gridsize3, blocksize3>>>(dev_spin_eo2_dn, dev_spin_eo2_dn);					// dev_spin_eo2_dn = gamma5 * dev_spin_eo2_dn
  #else
    	to_host(l_strange, dev_spin_eo2_up, h2d_spin_up, dev_spinsize);
    	gamma5(l_strange, l_strange, VOLUME/2);						// l_strange = gamma5 * l_strange
    	to_device(dev_spin_eo2_up, l_strange, h2d_spin_up, dev_spinsize);
    	
    	to_host(l_charm, dev_spin_eo2_dn, h2d_spin_dn, dev_spinsize);
    	gamma5(l_charm  , l_charm  , VOLUME/2);						// l_charm   = gamma5 * l_charm
    	to_device(dev_spin_eo2_dn, l_charm, h2d_spin_dn, dev_spinsize);
  #endif
  
  
  
  /*
  ////////////
  // output //										// output is already done by setting  dev_spin_eo2_up/dn = spinout_up/dn
  ////////////
  
  dev_copy_spinor_field<<<griddim2, blockdim2 >>>(dev_spin_eo2_up, spinout_up);		// spinout_up = dev_spin_eo2_up
  dev_copy_spinor_field<<<griddim2, blockdim2 >>>(dev_spin_eo2_dn, spinout_dn);		// spinout_dn = dev_spin_eo2_dn
  */
  
  
  return;
  
}//matrix_mpi_debug1()










// replaces matrix_multiplication32_mpi() for debugging !!

// RESULT: twice the memory requirements for some of the auxiliary host/device-interaction fields !!

void matrix_mpi_debug2 (dev_spinor * spinout_up, dev_spinor * spinout_dn,
                        dev_spinor * spinin_up , dev_spinor * spinin_dn ,
                        int gridsize1, int blocksize1, int gridsize2, int blocksize2,
                        int gridsize3, int blocksize3, int gridsize4, int blocksize4) {
  
  
  // we will use the auxiliary fields  dev_spin_eo{1,2}_up/dn  for working on and buffering
  // and set  dev_spin_eo2_up/dn  equal  spinout_up/dn
  // spinin_up/dn  have to remain unchanged !!
  // spinout_up/dn  can be freely used
  
  
  
  
  /////////////////////
  // LOCAL VARIABLES //
  /////////////////////
  
  int N_sites  =    VOLUME/2;		// #lattice sites
  int N_floats = 24*VOLUME/2;		// #floats
  
  // additional for the debugging purposes
  size_t dev_spinsize_int = VOLUME/2 * 6*sizeof(dev_spinor);
  size_t dev_spinsize_ext = (VOLUME+RAND)/2 * 6*sizeof(dev_spinor);
  size_t dev_spinsize     = dev_spinsize_int;
  
  
  
  ////////////////////////////////////
  //   MATCHING with Q_Qdagger_ND   //
  ////////////////////////////////////
  //                                //
  // _strange = _up                 //
  // _charm   = _dn                 //
  //                                //
  // DUM_MATRIX   = dev_spin_eo1_up //
  // DUM_MATRIX+1 = dev_spin_eo1_dn //
  //                                //
  // DUM_MATRIX+2 = dev_spin_eo2_up //
  // DUM_MATRIX+3 = dev_spin_eo2_dn //
  //                                //
  ////////////////////////////////////
  
  
  
  
  ///////////////////////////////////
  // INITIALIZATIONS & ASSIGNMENTS //	// have to use (one) other auxiliary field(s) than the calling function dev_cg_eo_nd
  ///////////////////////////////////
  
  dev_spin_eo2_up = spinout_up;		// need no memory allocated
  dev_spin_eo2_dn = spinout_dn;
  												///////////// THEORY ////////////////////////////////////////////////////////////////
  												//                                                                                 //
  												//  (Q_tilde) = gamma5 * ((M_oo) - (M_oe)(Mee^-1)(M_eo))                           //
  												//  (Q_tilde)(Q_tilde_dagger) * (up,dn) = (Q_tilde) * (b,a)                        //
  ///////////////										//                                                    (a,b) = (Q_tilde) * (dn,up)  //
  // MAIN BODY //										//                                                                                 //
  ///////////////										/////////////////////////////////////////////////////////////////////////////////////
  
  
  double nrm = 1.0 / (1.0 + g_mubar*g_mubar - g_epsbar*g_epsbar);
  
  
  if (g_proc_id == 0)
    printf("This is matrix_mpi_debug2(). ");
  
  
  spinor * help1  =  (spinor *) malloc((VOLUME+RAND) * 6*sizeof(dev_spinor));			// "/2" is missing because of "double"
  spinor * help2  =  (spinor *) malloc((VOLUME+RAND) * 6*sizeof(dev_spinor));
  spinor * help3  =  (spinor *) malloc((VOLUME+RAND) * 6*sizeof(dev_spinor));
  spinor * help4  =  (spinor *) malloc((VOLUME+RAND) * 6*sizeof(dev_spinor));

  
  
  ///////////////////////////////////////							/////////////////////////////////
  //        Q_tilde_dagger(2x2)        //							// (a,b) = (Q_tilde) * (dn,up) //
  ///////////////////////////////////////							/////////////////////////////////
  
  
  ////////////////////////////
  // (M_oe) (Mee^-1) (M_eo) //
  ////////////////////////////
  
  
  #ifndef CHECK_HOPPING_MATRIX
  			// xchange
  			#ifndef HOPPING_DEBUG
  			  xchange_field_wrapper(spinin_dn, 0);
  			#endif
  	// Flo:
  	#ifdef USETEXTURE
  	  bind_texture_spin(spinin_dn,1);
  	#endif
  	  #ifndef HOPPING_DEBUG
  	    dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, spinin_dn, dev_spin_eo1_up, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0);	// dev_spin_eo1_up = (M_eo) * spinin_dn
  	  #else
  	    Hopping_Matrix_wrapper(0, dev_spin_eo1_up, spinin_dn);
  	  #endif
  	#ifdef USETEXTURE
  	  unbind_texture_spin(1);
  	#endif
  			// xchange
  			#ifndef HOPPING_DEBUG
  			  xchange_field_wrapper(spinin_up, 0);
  			#endif
  	#ifdef USETEXTURE
  	  bind_texture_spin(spinin_up,1);
  	#endif
  	  #ifndef HOPPING_DEBUG
  	    dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, spinin_up, dev_spin_eo1_dn, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0);	// dev_spin_eo1_dn = (M_eo) * spinin_up
  	  #else
  	    Hopping_Matrix_wrapper(0, dev_spin_eo1_dn, spinin_up);
  	  #endif
  	#ifdef USETEXTURE
  	  unbind_texture_spin(1);
  	#endif
  #else
    	to_host(help1, spinin_dn, h2d_spin_up, dev_spinsize);
    	Hopping_Matrix(EO, help2, help1);				// g_spinor_field[DUM_MATRIX]   = (M_eo) * k_charm
    	to_device(dev_spin_eo1_up, help2, h2d_spin_up, dev_spinsize);
    	
    	to_host(help3, spinin_up, h2d_spin_dn, dev_spinsize);
    	Hopping_Matrix(EO, help4, help3);				// g_spinor_field[DUM_MATRIX+1] = (M_eo) * k_strange
    	to_device(dev_spin_eo1_dn, help4, h2d_spin_dn, dev_spinsize);
  #endif
  
  
  // written:
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo1_up, dev_spin_eo2_up, -1.0);		// dev_spin_eo2_up  =  (1 - imubar) * dev_spin_eo1_up  =  (1 - imubar)*(M_eo) * spinin_dn
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo1_dn, dev_spin_eo2_dn, +1.0);		// dev_spin_eo2_dn  =  (1 + imubar) * dev_spin_eo1_dn  =  (1 + imubar)*(M_eo) * spinin_up
  
  
  // CUBLAS:
  cublasSaxpy (N_floats, g_epsbar, (float *) dev_spin_eo1_dn, 1, (float *) dev_spin_eo2_up, 1);
  										// dev_spin_eo2_up  =  dev_spin_eo2_up  +  epsbar * dev_spin_eo1_dn  =  (1 - imubar)*(M_eo) * spinin_dn  +  epsbar * (M_eo) * spinin_up
  cublasSaxpy (N_floats, g_epsbar, (float *) dev_spin_eo1_up, 1, (float *) dev_spin_eo2_dn, 1);
  										// dev_spin_eo2_dn  =  dev_spin_eo2_dn  +  epsbar * dev_spin_eo1_up  =  (1 + imubar)*(M_eo) * spinin_up  +  epsbar * (M_eo) * spinin_dn
  
  // CUBLAS:
  cublasSscal (N_floats, nrm, (float *) dev_spin_eo2_up, 1);			// dev_spin_eo2_up  =  nrm * dev_spin_eo2_up  =  nrm*(1-imubar)*(M_eo) * spinin_dn  +  nrm*epsbar*(M_eo) * spinin_up
  
  cublasSscal (N_floats, nrm, (float *) dev_spin_eo2_dn, 1);			// dev_spin_eo2_dn  =  nrm * dev_spin_eo2_dn  =  nrm*(1+imubar)*(M_eo) * spinin_up  +  nrm*epsbar*(M_eo) * spinin_dn
  
  
  #ifndef CHECK_HOPPING_MATRIX
  			// xchange
  			#ifndef HOPPING_DEBUG
  			  xchange_field_wrapper(dev_spin_eo2_up, 1);
  			#endif
  	// Flo:
  	#ifdef USETEXTURE
  	  bind_texture_spin(dev_spin_eo2_up,1);									// remember: this is ((M_oe)(Mee^-1)(M_eo)) * (spinin_dn, spinin_up):
  	#endif
  	  #ifndef HOPPING_DEBUG
  	    dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo2_up, dev_spin_eo1_up, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1);
  	  #else
  	    Hopping_Matrix_wrapper(1, dev_spin_eo1_up, dev_spin_eo2_up);
  	  #endif
  	#ifdef USETEXTURE
  	  unbind_texture_spin(1);							// dev_spin_eo1_up  =  (M_oe) * dev_spin_eo2_up  =  (M_oe)*nrm*(1-imubar)*(M_eo) * spinin_dn  +  (M_oe)*nrm*epsbar*(M_eo) * spinin_up
  	#endif									
  			// xchange
  			#ifndef HOPPING_DEBUG
  			  xchange_field_wrapper(dev_spin_eo2_dn, 1);
  			#endif
  	#ifdef USETEXTURE
  	  bind_texture_spin(dev_spin_eo2_dn,1);
  	#endif
  	  #ifndef HOPPING_DEBUG
  	    dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo2_dn, dev_spin_eo1_dn, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1);
  	  #else
  	    Hopping_Matrix_wrapper(1, dev_spin_eo1_dn, dev_spin_eo2_dn);
  	  #endif
  	#ifdef USETEXTURE
  	  unbind_texture_spin(1);							// dev_spin_eo1_dn  =  (M_oe) * dev_spin_eo2_dn  =  (M_oe)*nrm*(1+imubar)*(M_eo) * spinin_up  +  (M_oe)*nrm*epsbar*(M_eo) * spinin_dn
  	#endif
  #else
    	to_host(help1, dev_spin_eo2_up, h2d_spin_up, dev_spinsize);
    	Hopping_Matrix(OE, help2, help1);		// g_spinor_field[DUM_MATRIX]   = (M_oe) * g_spinor_field[DUM_MATRIX+2]
    	to_device(dev_spin_eo1_up, help2, h2d_spin_up, dev_spinsize);
    	
    	to_host(help3, dev_spin_eo2_dn, h2d_spin_dn, dev_spinsize);
    	Hopping_Matrix(OE, help4, help3);		// g_spinor_field[DUM_MATRIX+1] = (M_oe) * g_spinor_field[DUM_MATRIX+3]
    	to_device(dev_spin_eo1_dn, help4, h2d_spin_dn, dev_spinsize);
  #endif
  
  
  
  ////////////
  // (M_oo) //
  ////////////
  
  
  // written:
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(spinin_dn, dev_spin_eo2_up, +1.0);			// dev_spin_eo2_up = (1 + imubar) * spinin_dn
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(spinin_up, dev_spin_eo2_dn, -1.0);			// dev_spin_eo2_dn = (1 - imubar) * spinin_up
  
  
  // CUBLAS:													// remember: this is (M_oo) * (spinin_dn, spinin_up):
  cublasSaxpy (N_floats, -g_epsbar, (float *) spinin_up, 1, (float *) dev_spin_eo2_up, 1);
  												// dev_spin_eo2_up  =  dev_spin_eo2_up - epsbar*spinin_up  =  (1+imubar)*spinin_dn - epsbar*spinin_up
  cublasSaxpy (N_floats, -g_epsbar, (float *) spinin_dn, 1, (float *) dev_spin_eo2_dn, 1);
  												// dev_spin_eo2_dn  =  dev_spin_eo2_dn - epsbar*spinin_dn  =  (1-imubar)*spinin_up - epsbar*spinin_dn
  
  
  
  ///////////////////////////////////////
  // (M_oo)  -  (M_oe) (Mee^-1) (M_eo) //
  ///////////////////////////////////////
  
  // CUBLAS:													// this is ((M_oo) - (M_oe)(Mee^-1)(M_eo)) * (spinin_dn, spinin_up):
  cublasSaxpy (N_floats, -1.0, (float *) dev_spin_eo1_up, 1, (float *) dev_spin_eo2_up, 1);
  							// dev_spin_eo2_up  =  dev_spin_eo2_up  -  dev_spin_eo1_up  =                      (1+imubar) * spinin_dn  -                    epsbar * spinin_up
  							//                                                             - (M_oe)*nrm*(1-imubar)*(M_eo) * spinin_dn  -  (M_oe)*nrm*epsbar*(M_eo) * spinin_up
  cublasSaxpy (N_floats, -1.0, (float *) dev_spin_eo1_dn, 1, (float *) dev_spin_eo2_dn, 1);
  							// dev_spin_eo2_dn  =  dev_spin_eo2_dn  -  dev_spin_eo1_dn  =                      (1-imubar) * spinin_up  -                    epsbar * spinin_dn
  							//                                                             - (M_oe)*nrm*(1+imubar)*(M_eo) * spinin_up  -  (M_oe)*nrm*epsbar*(M_eo) * spinin_dn
  
  
  ////////////
  // gamma5 //
  ////////////
  
  // Flo:
  dev_gamma5<<<gridsize3, blocksize3>>>(dev_spin_eo2_up, dev_spin_eo2_up);					// dev_spin_eo2_up = gamma5 * dev_spin_eo2_up 
  dev_gamma5<<<gridsize3, blocksize3>>>(dev_spin_eo2_dn, dev_spin_eo2_dn);					// dev_spin_eo2_dn = gamma5 * dev_spin_eo2_dn
  	
  
  
  
  
  
  ////////////////////
  // (a,b) -> (b,a) //
  ////////////////////
  dev_copy_spinor_field<<<gridsize4, blocksize4>>>(dev_spin_eo2_dn, dev_spin_eo3_up);			// dev_spin_eo3_up = dev_spin_eo2_dn
  dev_copy_spinor_field<<<gridsize4, blocksize4>>>(dev_spin_eo2_up, dev_spin_eo3_dn);			// dev_spin_eo3_dn = dev_spin_eo2_up
  
  
  
  
  
  
  ///////////////////////////////////								///////////////////////
  //        Q_tilde(2x2)           //								// (Q_tilde) * (b,a) //
  ///////////////////////////////////								///////////////////////
  
  
  ////////////////////////////
  // (M_oe) (Mee^-1) (M_eo) //
  ////////////////////////////
  
  #ifndef CHECK_HOPPING_MATRIX
  			// xchange
  			#ifndef HOPPING_DEBUG
  			  xchange_field_wrapper(dev_spin_eo3_up, 0);
  			#endif
  	// Flo:
  	#ifdef USETEXTURE
  	  bind_texture_spin(dev_spin_eo3_up,1);
  	#endif
  	  #ifndef HOPPING_DEBUG
  	    dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo3_up, dev_spin_eo1_up, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0);	// dev_spin_eo1_up = (M_eo) * dev_spin_eo3_up
  	  #else
  	    Hopping_Matrix_wrapper(0, dev_spin_eo1_up, dev_spin_eo3_up);
  	  #endif
  	#ifdef USETEXTURE
  	  unbind_texture_spin(1);
  	#endif
  			// xchange
  			#ifndef HOPPING_DEBUG
  			  xchange_field_wrapper(dev_spin_eo3_dn, 0);
  			#endif
  	#ifdef USETEXTURE
  	  bind_texture_spin(dev_spin_eo3_dn,1);
  	#endif
  	  #ifndef HOPPING_DEBUG
  	    dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo3_dn, dev_spin_eo1_dn, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0);	// dev_spin_eo1_dn = (M_eo) * dev_spin_eo3_dn
  	  #else
  	    Hopping_Matrix_wrapper(0, dev_spin_eo1_dn, dev_spin_eo3_dn);
  	  #endif
  	#ifdef USETEXTURE
  	  unbind_texture_spin(1);
  	#endif
  #else
    	to_host(help1, dev_spin_eo3_up, h2d_spin_up, dev_spinsize);
    	Hopping_Matrix(EO, help2, help1);	// g_spinor_field[DUM_MATRIX]   = (M_eo) * g_spinor_field[DUM_MATRIX+7]
    	to_device(dev_spin_eo1_up, help2, h2d_spin_up, dev_spinsize);
    	
    	to_host(help3, dev_spin_eo3_dn, h2d_spin_dn, dev_spinsize);
    	Hopping_Matrix(EO, help4, help3);	// g_spinor_field[DUM_MATRIX+1] = (M_eo) * g_spinor_field[DUM_MATRIX+6]
    	to_device(dev_spin_eo1_dn, help4, h2d_spin_dn, dev_spinsize);
  #endif
  
  
  // written:
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo1_up, dev_spin_eo2_up, -1.0);			// dev_spin_eo2_up  =  (1 - imubar) * dev_spin_eo1_up  =  (1 - imubar)*(M_eo) * dev_spin_eo3_up
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo1_dn, dev_spin_eo2_dn, +1.0);			// dev_spin_eo2_dn  =  (1 + imubar) * dev_spin_eo1_dn  =  (1 + imubar)*(M_eo) * dev_spin_eo3_dn
  
  
  // CUBLAS:
  cublasSaxpy (N_floats, g_epsbar, (float *) dev_spin_eo1_dn, 1, (float *) dev_spin_eo2_up, 1);
  										// dev_spin_eo2_up  =  dev_spin_eo2_up  +  epsbar * dev_spin_eo1_dn  =  (1 - imubar)*(M_eo) * dev_spin_eo3_up  +  epsbar * (M_eo) * dev_spin_eo3_dn
  cublasSaxpy (N_floats, g_epsbar, (float *) dev_spin_eo1_up, 1, (float *) dev_spin_eo2_dn, 1);
  										// dev_spin_eo2_dn  =  dev_spin_eo2_dn  +  epsbar * dev_spin_eo1_up  =  (1 + imubar)*(M_eo) * dev_spin_eo3_dn  +  epsbar * (M_eo) * dev_spin_eo3_up
  
  // CUBLAS:
  cublasSscal (N_floats, nrm, (float *) dev_spin_eo2_up, 1);			// dev_spin_eo2_up  =  nrm * dev_spin_eo2_up  =  nrm*(1-imubar)*(M_eo) * dev_spin_eo3_up  +  nrm*epsbar*(M_eo) * dev_spin_eo3_dn
  
  cublasSscal (N_floats, nrm, (float *) dev_spin_eo2_dn, 1);			// dev_spin_eo2_dn  =  nrm * dev_spin_eo2_dn  =  nrm*(1+imubar)*(M_eo) * dev_spin_eo3_dn  +  nrm*epsbar*(M_eo) * dev_spin_eo3_up
  
  
  #ifndef CHECK_HOPPING_MATRIX
  			// xchange
  			#ifndef HOPPING_DEBUG
  			  xchange_field_wrapper(dev_spin_eo2_up, 1);
  			#endif
  	// Flo:
  	#ifdef USETEXTURE
  	  bind_texture_spin(dev_spin_eo2_up,1);									// remember: this is ((M_oe) (Mee^-1) (M_eo)) * (dev_spin_eo3_up, dev_spin_eo3_dn):
  	#endif
  	  #ifndef HOPPING_DEBUG
  	    dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo2_up, dev_spin_eo1_up, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1);
  	  #else
  	    Hopping_Matrix_wrapper(1, dev_spin_eo1_up, dev_spin_eo2_up);
  	  #endif
  	#ifdef USETEXTURE
  	  unbind_texture_spin(1);							// dev_spin_eo1_up  =  (M_oe) * dev_spin_eo2_up  =  (M_oe)*nrm*(1-imubar)*(M_eo) * dev_spin_eo3_up  +  (M_oe)*nrm*epsbar*(M_eo) * dev_spin_eo3_dn
  	#endif
  			// xchange
  			#ifndef HOPPING_DEBUG
  			  xchange_field_wrapper(dev_spin_eo2_dn, 1);
  			#endif
  	#ifdef USETEXTURE
  	  bind_texture_spin(dev_spin_eo2_dn,1);
  	#endif
  	  #ifndef HOPPING_DEBUG
  	    dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo2_dn, dev_spin_eo1_dn, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1);
  	  #else
  	    Hopping_Matrix_wrapper(1, dev_spin_eo1_dn, dev_spin_eo2_dn);
  	  #endif
  	#ifdef USETEXTURE
  	  unbind_texture_spin(1);							// dev_spin_eo1_dn  =  (M_oe) * dev_spin_eo2_dn  =  (M_oe)*nrm*(1+imubar)*(M_eo) * dev_spin_eo3_dn  +  (M_oe)*nrm*epsbar*(M_eo) * dev_spin_eo3_up
  	#endif
  #else
    	to_host(help1, dev_spin_eo2_up, h2d_spin_up, dev_spinsize);
    	Hopping_Matrix(OE, help2, help1);				// l_strange = (M_oe) * g_spinor_field[DUM_MATRIX+2]
    	to_device(dev_spin_eo1_up, help2, h2d_spin_up, dev_spinsize);
    	
    	to_host(help3, dev_spin_eo2_dn, h2d_spin_dn, dev_spinsize);
    	Hopping_Matrix(OE, help4, help3);				// l_charm   = (M_oe) * g_spinor_field[DUM_MATRIX+3]
    	to_device(dev_spin_eo1_dn, help4, h2d_spin_dn, dev_spinsize);
  #endif
  
  
  ////////////
  // (M_oo) //
  ////////////
  
  // written:
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo3_up, dev_spin_eo2_up, +1.0);				// dev_spin_eo2_up = (1 + imubar) * dev_spin_eo3_up
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo3_dn, dev_spin_eo2_dn, -1.0);				// dev_spin_eo2_dn = (1 - imubar) * dev_spin_eo3_dn
  
  
  // CUBLAS:											// remember: this is (M_oo) * (dev_spin_eo3_up, dev_spin_eo3_dn):
  cublasSaxpy (N_floats, -g_epsbar, (float *) dev_spin_eo3_dn, 1, (float *) dev_spin_eo2_up, 1);
  												// dev_spin_eo2_up  =  dev_spin_eo2_up - epsbar*dev_spin_eo3_dn  =  (1+imubar)*dev_spin_eo3_up - epsbar*dev_spin_eo3_dn
  cublasSaxpy (N_floats, -g_epsbar, (float *) dev_spin_eo3_up, 1, (float *) dev_spin_eo2_dn, 1);
  												// dev_spin_eo2_dn  =  dev_spin_eo2_dn - epsbar*dev_spin_eo3_up  =  (1-imubar)*dev_spin_eo3_dn - epsbar*dev_spin_eo3_up
  
  
  
  ///////////////////////////////////////
  // (M_oo)  -  (M_oe) (Mee^-1) (M_eo) //
  ///////////////////////////////////////
  
  // CUBLAS:											// this is ( (M_oo) - (M_oe) (Mee^-1) (M_eo) ) * (dev_spin_eo3_up, dev_spin_eo3_dn)
  cublasSaxpy (N_floats, -1.0, (float *) dev_spin_eo1_up, 1, (float *) dev_spin_eo2_up, 1);
  							// dev_spin_eo2_up  =  dev_spin_eo2_up  -  dev_spin_eo1_up  =                      (1+imubar) * dev_spin_eo3_up  -                    epsbar * dev_spin_eo3_dn
  							//                                                             - (M_oe)*nrm*(1-imubar)*(M_eo) * dev_spin_eo3_up  -  (M_oe)*nrm*epsbar*(M_eo) * dev_spin_eo3_dn
  cublasSaxpy (N_floats, -1.0, (float *) dev_spin_eo1_dn, 1, (float *) dev_spin_eo2_dn, 1);
  							// dev_spin_eo2_dn  =  dev_spin_eo2_dn  -  dev_spin_eo1_dn  =                      (1-imubar) * dev_spin_eo3_dn  -                    epsbar * dev_spin_eo3_up
  							//                                                             - (M_oe)*nrm*(1+imubar)*(M_eo) * dev_spin_eo3_dn  -  (M_oe)*nrm*epsbar*(M_eo) * dev_spin_eo3_up
  
  
  ////////////
  // gamma5 //
  ////////////
  
  // Flo:
  dev_gamma5<<<gridsize3, blocksize3>>>(dev_spin_eo2_up, dev_spin_eo2_up);					// dev_spin_eo2_up = gamma5 * dev_spin_eo2_up 
  dev_gamma5<<<gridsize3, blocksize3>>>(dev_spin_eo2_dn, dev_spin_eo2_dn);					// dev_spin_eo2_dn = gamma5 * dev_spin_eo2_dn
  
  
  
  
  /*
  ////////////
  // output //										// output is already done by setting  dev_spin_eo2_up/dn = spinout_up/dn
  ////////////
  
  dev_copy_spinor_field<<<griddim2, blockdim2 >>>(dev_spin_eo2_up, spinout_up);		// spinout_up = dev_spin_eo2_up
  dev_copy_spinor_field<<<griddim2, blockdim2 >>>(dev_spin_eo2_dn, spinout_dn);		// spinout_dn = dev_spin_eo2_dn
  */
  
  
  return;
  
}//matrix_mpi_debug2()








// replaces matrix_multiplication32_mpi() for debugging !!

// is for the moment an identical copy of matrix_mpi_debug2()
//	will be changed for testing the alternative version of the hopping matrix ...
//	will for now test the xchange_fields_mpi() procedure
//
//	one therefore has to turn of the xchange_fields() in the CPU's Hopping_Matrix() routine !!
//		we are most likely use version 5 of Hopping_Matrix() in Hopping_Matrix.c

#define check_xchange

void matrix_mpi_debug3 (dev_spinor * spinout_up, dev_spinor * spinout_dn,
                        dev_spinor * spinin_up , dev_spinor * spinin_dn ,
                        int gridsize1, int blocksize1, int gridsize2, int blocksize2,
                        int gridsize3, int blocksize3, int gridsize4, int blocksize4) {
  
  
  // we will use the auxiliary fields  dev_spin_eo{1,2}_up/dn  for working on and buffering
  // and set  dev_spin_eo2_up/dn  equal  spinout_up/dn
  // spinin_up/dn  have to remain unchanged !!
  // spinout_up/dn  can be freely used
  
  
  
  
  /////////////////////
  // LOCAL VARIABLES //
  /////////////////////
  
  int N_sites  =    VOLUME/2;		// #lattice sites
  int N_floats = 24*VOLUME/2;		// #floats
  
  // additional for the debugging purposes
  size_t dev_spinsize_int = VOLUME/2 * 6*sizeof(dev_spinor);
  size_t dev_spinsize_ext = (VOLUME+RAND)/2 * 6*sizeof(dev_spinor);
  size_t dev_spinsize     = dev_spinsize_int;
  
  
  
  ////////////////////////////////////
  //   MATCHING with Q_Qdagger_ND   //
  ////////////////////////////////////
  //                                //
  // _strange = _up                 //
  // _charm   = _dn                 //
  //                                //
  // DUM_MATRIX   = dev_spin_eo1_up //
  // DUM_MATRIX+1 = dev_spin_eo1_dn //
  //                                //
  // DUM_MATRIX+2 = dev_spin_eo2_up //
  // DUM_MATRIX+3 = dev_spin_eo2_dn //
  //                                //
  ////////////////////////////////////
  
  
  
  
  ///////////////////////////////////
  // INITIALIZATIONS & ASSIGNMENTS //	// have to use (one) other auxiliary field(s) than the calling function dev_cg_eo_nd
  ///////////////////////////////////
  
  dev_spin_eo2_up = spinout_up;		// need no memory allocated
  dev_spin_eo2_dn = spinout_dn;
  												///////////// THEORY ////////////////////////////////////////////////////////////////
  												//                                                                                 //
  												//  (Q_tilde) = gamma5 * ((M_oo) - (M_oe)(Mee^-1)(M_eo))                           //
  												//  (Q_tilde)(Q_tilde_dagger) * (up,dn) = (Q_tilde) * (b,a)                        //
  ///////////////										//                                                    (a,b) = (Q_tilde) * (dn,up)  //
  // MAIN BODY //										//                                                                                 //
  ///////////////										/////////////////////////////////////////////////////////////////////////////////////
  
  
  double nrm = 1.0 / (1.0 + g_mubar*g_mubar - g_epsbar*g_epsbar);
  
  
  if (g_proc_id == 0)
    printf("This is matrix_mpi_debug3(). ");
  
  
  spinor * help1  =  (spinor *) malloc((VOLUME+RAND) * 6*sizeof(dev_spinor));			// "/2" is missing because of "double"
  spinor * help2  =  (spinor *) malloc((VOLUME+RAND) * 6*sizeof(dev_spinor));
  spinor * help3  =  (spinor *) malloc((VOLUME+RAND) * 6*sizeof(dev_spinor));
  spinor * help4  =  (spinor *) malloc((VOLUME+RAND) * 6*sizeof(dev_spinor));

  
  
  ///////////////////////////////////////							/////////////////////////////////
  //        Q_tilde_dagger(2x2)        //							// (a,b) = (Q_tilde) * (dn,up) //
  ///////////////////////////////////////							/////////////////////////////////
  
  
  ////////////////////////////
  // (M_oe) (Mee^-1) (M_eo) //
  ////////////////////////////
  
  
  #ifndef CHECK_HOPPING_MATRIX

  	  #ifndef HOPPING_DEBUG
  	    xchange_field_wrapper(spinin_dn, 0);
  	    dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, spinin_dn, dev_spin_eo1_up, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0);	// dev_spin_eo1_up = (M_eo) * spinin_dn
  	  #else
  	    printf("Does it work?\n");
  	    xchange_field_wrapper(spinin_dn, 0);
  	    Hopping_Matrix_wrapper(0, dev_spin_eo1_up, spinin_dn);
  	  #endif

  	  #ifndef HOPPING_DEBUG
  	    xchange_field_wrapper(spinin_up, 0);
  	    dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, spinin_up, dev_spin_eo1_dn, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0);	// dev_spin_eo1_dn = (M_eo) * spinin_up
  	  #else
  	    xchange_field_wrapper(spinin_up, 0);
  	    Hopping_Matrix_wrapper(0, dev_spin_eo1_dn, spinin_up);
  	  #endif
  #else
    	to_host(help1, spinin_dn, h2d_spin_up, dev_spinsize);
    	Hopping_Matrix(EO, help2, help1);				// g_spinor_field[DUM_MATRIX]   = (M_eo) * k_charm
    	to_device(dev_spin_eo1_up, help2, h2d_spin_up, dev_spinsize);
    	
    	to_host(help3, spinin_up, h2d_spin_dn, dev_spinsize);
    	Hopping_Matrix(EO, help4, help3);				// g_spinor_field[DUM_MATRIX+1] = (M_eo) * k_strange
    	to_device(dev_spin_eo1_dn, help4, h2d_spin_dn, dev_spinsize);
  #endif
  
  
  // written:
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo1_up, dev_spin_eo2_up, -1.0);		// dev_spin_eo2_up  =  (1 - imubar) * dev_spin_eo1_up  =  (1 - imubar)*(M_eo) * spinin_dn
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo1_dn, dev_spin_eo2_dn, +1.0);		// dev_spin_eo2_dn  =  (1 + imubar) * dev_spin_eo1_dn  =  (1 + imubar)*(M_eo) * spinin_up
  
  
  // CUBLAS:
  cublasSaxpy (N_floats, g_epsbar, (float *) dev_spin_eo1_dn, 1, (float *) dev_spin_eo2_up, 1);
  										// dev_spin_eo2_up  =  dev_spin_eo2_up  +  epsbar * dev_spin_eo1_dn  =  (1 - imubar)*(M_eo) * spinin_dn  +  epsbar * (M_eo) * spinin_up
  cublasSaxpy (N_floats, g_epsbar, (float *) dev_spin_eo1_up, 1, (float *) dev_spin_eo2_dn, 1);
  										// dev_spin_eo2_dn  =  dev_spin_eo2_dn  +  epsbar * dev_spin_eo1_up  =  (1 + imubar)*(M_eo) * spinin_up  +  epsbar * (M_eo) * spinin_dn
  
  // CUBLAS:
  cublasSscal (N_floats, nrm, (float *) dev_spin_eo2_up, 1);			// dev_spin_eo2_up  =  nrm * dev_spin_eo2_up  =  nrm*(1-imubar)*(M_eo) * spinin_dn  +  nrm*epsbar*(M_eo) * spinin_up
  
  cublasSscal (N_floats, nrm, (float *) dev_spin_eo2_dn, 1);			// dev_spin_eo2_dn  =  nrm * dev_spin_eo2_dn  =  nrm*(1+imubar)*(M_eo) * spinin_up  +  nrm*epsbar*(M_eo) * spinin_dn
  
  
  #ifndef CHECK_HOPPING_MATRIX

  	  #ifndef HOPPING_DEBUG
  	    xchange_field_wrapper(dev_spin_eo2_up, 1);
  	    dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo2_up, dev_spin_eo1_up, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1);
  	  #else
  	    xchange_field_wrapper(dev_spin_eo2_up, 1);
  	    Hopping_Matrix_wrapper(1, dev_spin_eo1_up, dev_spin_eo2_up);
  	  #endif							

  	  #ifndef HOPPING_DEBUG
  	    xchange_field_wrapper(dev_spin_eo2_dn, 1);
  	    dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo2_dn, dev_spin_eo1_dn, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1);
  	  #else
  	    xchange_field_wrapper(dev_spin_eo2_dn, 1);
  	    Hopping_Matrix_wrapper(1, dev_spin_eo1_dn, dev_spin_eo2_dn);
  	  #endif
  #else
    	to_host(help1, dev_spin_eo2_up, h2d_spin_up, dev_spinsize);
    	Hopping_Matrix(OE, help2, help1);		// g_spinor_field[DUM_MATRIX]   = (M_oe) * g_spinor_field[DUM_MATRIX+2]
    	to_device(dev_spin_eo1_up, help2, h2d_spin_up, dev_spinsize);
    	
    	to_host(help3, dev_spin_eo2_dn, h2d_spin_dn, dev_spinsize);
    	Hopping_Matrix(OE, help4, help3);		// g_spinor_field[DUM_MATRIX+1] = (M_oe) * g_spinor_field[DUM_MATRIX+3]
    	to_device(dev_spin_eo1_dn, help4, h2d_spin_dn, dev_spinsize);
  #endif
  
  
  
  ////////////
  // (M_oo) //
  ////////////
  
  
  // written:
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(spinin_dn, dev_spin_eo2_up, +1.0);			// dev_spin_eo2_up = (1 + imubar) * spinin_dn
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(spinin_up, dev_spin_eo2_dn, -1.0);			// dev_spin_eo2_dn = (1 - imubar) * spinin_up
  
  
  // CUBLAS:													// remember: this is (M_oo) * (spinin_dn, spinin_up):
  cublasSaxpy (N_floats, -g_epsbar, (float *) spinin_up, 1, (float *) dev_spin_eo2_up, 1);
  												// dev_spin_eo2_up  =  dev_spin_eo2_up - epsbar*spinin_up  =  (1+imubar)*spinin_dn - epsbar*spinin_up
  cublasSaxpy (N_floats, -g_epsbar, (float *) spinin_dn, 1, (float *) dev_spin_eo2_dn, 1);
  												// dev_spin_eo2_dn  =  dev_spin_eo2_dn - epsbar*spinin_dn  =  (1-imubar)*spinin_up - epsbar*spinin_dn
  
  
  
  ///////////////////////////////////////
  // (M_oo)  -  (M_oe) (Mee^-1) (M_eo) //
  ///////////////////////////////////////
  
  // CUBLAS:													// this is ((M_oo) - (M_oe)(Mee^-1)(M_eo)) * (spinin_dn, spinin_up):
  cublasSaxpy (N_floats, -1.0, (float *) dev_spin_eo1_up, 1, (float *) dev_spin_eo2_up, 1);
  							// dev_spin_eo2_up  =  dev_spin_eo2_up  -  dev_spin_eo1_up  =                      (1+imubar) * spinin_dn  -                    epsbar * spinin_up
  							//                                                             - (M_oe)*nrm*(1-imubar)*(M_eo) * spinin_dn  -  (M_oe)*nrm*epsbar*(M_eo) * spinin_up
  cublasSaxpy (N_floats, -1.0, (float *) dev_spin_eo1_dn, 1, (float *) dev_spin_eo2_dn, 1);
  							// dev_spin_eo2_dn  =  dev_spin_eo2_dn  -  dev_spin_eo1_dn  =                      (1-imubar) * spinin_up  -                    epsbar * spinin_dn
  							//                                                             - (M_oe)*nrm*(1+imubar)*(M_eo) * spinin_up  -  (M_oe)*nrm*epsbar*(M_eo) * spinin_dn
  
  
  ////////////
  // gamma5 //
  ////////////
  
  // Flo:
  dev_gamma5<<<gridsize3, blocksize3>>>(dev_spin_eo2_up, dev_spin_eo2_up);					// dev_spin_eo2_up = gamma5 * dev_spin_eo2_up 
  dev_gamma5<<<gridsize3, blocksize3>>>(dev_spin_eo2_dn, dev_spin_eo2_dn);					// dev_spin_eo2_dn = gamma5 * dev_spin_eo2_dn
  	
  
  
  
  
  
  ////////////////////
  // (a,b) -> (b,a) //
  ////////////////////
  dev_copy_spinor_field<<<gridsize4, blocksize4>>>(dev_spin_eo2_dn, dev_spin_eo3_up);			// dev_spin_eo3_up = dev_spin_eo2_dn
  dev_copy_spinor_field<<<gridsize4, blocksize4>>>(dev_spin_eo2_up, dev_spin_eo3_dn);			// dev_spin_eo3_dn = dev_spin_eo2_up
  
  
  
  
  
  
  ///////////////////////////////////								///////////////////////
  //        Q_tilde(2x2)           //								// (Q_tilde) * (b,a) //
  ///////////////////////////////////								///////////////////////
  
  
  ////////////////////////////
  // (M_oe) (Mee^-1) (M_eo) //
  ////////////////////////////
  
  #ifndef CHECK_HOPPING_MATRIX

  	  #ifndef HOPPING_DEBUG
  	    xchange_field_wrapper(dev_spin_eo3_up, 0);
  	    dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo3_up, dev_spin_eo1_up, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0);	// dev_spin_eo1_up = (M_eo) * dev_spin_eo3_up
  	  #else
  	    xchange_field_wrapper(dev_spin_eo3_up, 0);
  	    Hopping_Matrix_wrapper(0, dev_spin_eo1_up, dev_spin_eo3_up);
  	  #endif

  	  #ifndef HOPPING_DEBUG
  	    xchange_field_wrapper(dev_spin_eo3_dn, 0);
  	    dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo3_dn, dev_spin_eo1_dn, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0);	// dev_spin_eo1_dn = (M_eo) * dev_spin_eo3_dn
  	  #else
  	    xchange_field_wrapper(dev_spin_eo3_dn, 0);
  	    Hopping_Matrix_wrapper(0, dev_spin_eo1_dn, dev_spin_eo3_dn);
  	  #endif
  #else
    	to_host(help1, dev_spin_eo3_up, h2d_spin_up, dev_spinsize);
    	Hopping_Matrix(EO, help2, help1);	// g_spinor_field[DUM_MATRIX]   = (M_eo) * g_spinor_field[DUM_MATRIX+7]
    	to_device(dev_spin_eo1_up, help2, h2d_spin_up, dev_spinsize);
    	
    	to_host(help3, dev_spin_eo3_dn, h2d_spin_dn, dev_spinsize);
    	Hopping_Matrix(EO, help4, help3);	// g_spinor_field[DUM_MATRIX+1] = (M_eo) * g_spinor_field[DUM_MATRIX+6]
    	to_device(dev_spin_eo1_dn, help4, h2d_spin_dn, dev_spinsize);
  #endif
  
  
  // written:
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo1_up, dev_spin_eo2_up, -1.0);			// dev_spin_eo2_up  =  (1 - imubar) * dev_spin_eo1_up  =  (1 - imubar)*(M_eo) * dev_spin_eo3_up
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo1_dn, dev_spin_eo2_dn, +1.0);			// dev_spin_eo2_dn  =  (1 + imubar) * dev_spin_eo1_dn  =  (1 + imubar)*(M_eo) * dev_spin_eo3_dn
  
  
  // CUBLAS:
  cublasSaxpy (N_floats, g_epsbar, (float *) dev_spin_eo1_dn, 1, (float *) dev_spin_eo2_up, 1);
  										// dev_spin_eo2_up  =  dev_spin_eo2_up  +  epsbar * dev_spin_eo1_dn  =  (1 - imubar)*(M_eo) * dev_spin_eo3_up  +  epsbar * (M_eo) * dev_spin_eo3_dn
  cublasSaxpy (N_floats, g_epsbar, (float *) dev_spin_eo1_up, 1, (float *) dev_spin_eo2_dn, 1);
  										// dev_spin_eo2_dn  =  dev_spin_eo2_dn  +  epsbar * dev_spin_eo1_up  =  (1 + imubar)*(M_eo) * dev_spin_eo3_dn  +  epsbar * (M_eo) * dev_spin_eo3_up
  
  // CUBLAS:
  cublasSscal (N_floats, nrm, (float *) dev_spin_eo2_up, 1);			// dev_spin_eo2_up  =  nrm * dev_spin_eo2_up  =  nrm*(1-imubar)*(M_eo) * dev_spin_eo3_up  +  nrm*epsbar*(M_eo) * dev_spin_eo3_dn
  
  cublasSscal (N_floats, nrm, (float *) dev_spin_eo2_dn, 1);			// dev_spin_eo2_dn  =  nrm * dev_spin_eo2_dn  =  nrm*(1+imubar)*(M_eo) * dev_spin_eo3_dn  +  nrm*epsbar*(M_eo) * dev_spin_eo3_up
  
  
  #ifndef CHECK_HOPPING_MATRIX

  	  #ifndef HOPPING_DEBUG
  	    xchange_field_wrapper(dev_spin_eo2_up, 1);
  	    dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo2_up, dev_spin_eo1_up, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1);
  	  #else
  	    xchange_field_wrapper(dev_spin_eo2_up, 1);
  	    Hopping_Matrix_wrapper(1, dev_spin_eo1_up, dev_spin_eo2_up);
  	  #endif

  	  #ifndef HOPPING_DEBUG
  	    xchange_field_wrapper(dev_spin_eo2_dn, 1);
  	    dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo2_dn, dev_spin_eo1_dn, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1);
  	  #else
  	    xchange_field_wrapper(dev_spin_eo2_dn, 1);
  	    Hopping_Matrix_wrapper(1, dev_spin_eo1_dn, dev_spin_eo2_dn);
  	  #endif
  #else
    	to_host(help1, dev_spin_eo2_up, h2d_spin_up, dev_spinsize);
    	Hopping_Matrix(OE, help2, help1);				// l_strange = (M_oe) * g_spinor_field[DUM_MATRIX+2]
    	to_device(dev_spin_eo1_up, help2, h2d_spin_up, dev_spinsize);
    	
    	to_host(help3, dev_spin_eo2_dn, h2d_spin_dn, dev_spinsize);
    	Hopping_Matrix(OE, help4, help3);				// l_charm   = (M_oe) * g_spinor_field[DUM_MATRIX+3]
    	to_device(dev_spin_eo1_dn, help4, h2d_spin_dn, dev_spinsize);
  #endif
  
  
  ////////////
  // (M_oo) //
  ////////////
  
  // written:
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo3_up, dev_spin_eo2_up, +1.0);				// dev_spin_eo2_up = (1 + imubar) * dev_spin_eo3_up
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo3_dn, dev_spin_eo2_dn, -1.0);				// dev_spin_eo2_dn = (1 - imubar) * dev_spin_eo3_dn
  
  
  // CUBLAS:											// remember: this is (M_oo) * (dev_spin_eo3_up, dev_spin_eo3_dn):
  cublasSaxpy (N_floats, -g_epsbar, (float *) dev_spin_eo3_dn, 1, (float *) dev_spin_eo2_up, 1);
  												// dev_spin_eo2_up  =  dev_spin_eo2_up - epsbar*dev_spin_eo3_dn  =  (1+imubar)*dev_spin_eo3_up - epsbar*dev_spin_eo3_dn
  cublasSaxpy (N_floats, -g_epsbar, (float *) dev_spin_eo3_up, 1, (float *) dev_spin_eo2_dn, 1);
  												// dev_spin_eo2_dn  =  dev_spin_eo2_dn - epsbar*dev_spin_eo3_up  =  (1-imubar)*dev_spin_eo3_dn - epsbar*dev_spin_eo3_up
  
  
  
  ///////////////////////////////////////
  // (M_oo)  -  (M_oe) (Mee^-1) (M_eo) //
  ///////////////////////////////////////
  
  // CUBLAS:											// this is ( (M_oo) - (M_oe) (Mee^-1) (M_eo) ) * (dev_spin_eo3_up, dev_spin_eo3_dn)
  cublasSaxpy (N_floats, -1.0, (float *) dev_spin_eo1_up, 1, (float *) dev_spin_eo2_up, 1);
  							// dev_spin_eo2_up  =  dev_spin_eo2_up  -  dev_spin_eo1_up  =                      (1+imubar) * dev_spin_eo3_up  -                    epsbar * dev_spin_eo3_dn
  							//                                                             - (M_oe)*nrm*(1-imubar)*(M_eo) * dev_spin_eo3_up  -  (M_oe)*nrm*epsbar*(M_eo) * dev_spin_eo3_dn
  cublasSaxpy (N_floats, -1.0, (float *) dev_spin_eo1_dn, 1, (float *) dev_spin_eo2_dn, 1);
  							// dev_spin_eo2_dn  =  dev_spin_eo2_dn  -  dev_spin_eo1_dn  =                      (1-imubar) * dev_spin_eo3_dn  -                    epsbar * dev_spin_eo3_up
  							//                                                             - (M_oe)*nrm*(1+imubar)*(M_eo) * dev_spin_eo3_dn  -  (M_oe)*nrm*epsbar*(M_eo) * dev_spin_eo3_up
  
  
  ////////////
  // gamma5 //
  ////////////
  
  // Flo:
  dev_gamma5<<<gridsize3, blocksize3>>>(dev_spin_eo2_up, dev_spin_eo2_up);					// dev_spin_eo2_up = gamma5 * dev_spin_eo2_up 
  dev_gamma5<<<gridsize3, blocksize3>>>(dev_spin_eo2_dn, dev_spin_eo2_dn);					// dev_spin_eo2_dn = gamma5 * dev_spin_eo2_dn
  
  
  
  
  /*
  ////////////
  // output //										// output is already done by setting  dev_spin_eo2_up/dn = spinout_up/dn
  ////////////
  
  dev_copy_spinor_field<<<griddim2, blockdim2 >>>(dev_spin_eo2_up, spinout_up);		// spinout_up = dev_spin_eo2_up
  dev_copy_spinor_field<<<griddim2, blockdim2 >>>(dev_spin_eo2_dn, spinout_dn);		// spinout_dn = dev_spin_eo2_dn
  */
  
  
  return;
  
}//matrix_mpi_debug3()







// replaces matrix_multiplication32_mpi() for debugging !!

// is for the moment an identical copy of matrix_mpi_debug2()
//	will be changed for testing the alternative version of the hopping matrix ...
//	will for now test the xchange_fields_mpi() procedure
//
//	one therefore has to turn of the xchange_fields() in the CPU's Hopping_Matrix() routine !!
//		we are most likely use version 5 of Hopping_Matrix() in Hopping_Matrix.c

void matrix_mpi_debug4 (dev_spinor * spinout_up, dev_spinor * spinout_dn,
                        dev_spinor * spinin_up , dev_spinor * spinin_dn ,
                        int gridsize1, int blocksize1, int gridsize2, int blocksize2,
                        int gridsize3, int blocksize3, int gridsize4, int blocksize4) {
  
  
  // we will use the auxiliary fields  dev_spin_eo{1,2}_up/dn  for working on and buffering
  // and set  dev_spin_eo2_up/dn  equal  spinout_up/dn
  // spinin_up/dn  have to remain unchanged !!
  // spinout_up/dn  can be freely used
  
  
  
  
  /////////////////////
  // LOCAL VARIABLES //
  /////////////////////
  
  int N_sites  =    VOLUME/2;		// #lattice sites
  int N_floats = 24*VOLUME/2;		// #floats
  
  // additional for the debugging purposes
  size_t dev_spinsize_int = VOLUME/2 * 6*sizeof(dev_spinor);
  size_t dev_spinsize_ext = (VOLUME+RAND)/2 * 6*sizeof(dev_spinor);
  size_t dev_spinsize     = dev_spinsize_int;
  
  
  
  ////////////////////////////////////
  //   MATCHING with Q_Qdagger_ND   //
  ////////////////////////////////////
  //                                //
  // _strange = _up                 //
  // _charm   = _dn                 //
  //                                //
  // DUM_MATRIX   = dev_spin_eo1_up //
  // DUM_MATRIX+1 = dev_spin_eo1_dn //
  //                                //
  // DUM_MATRIX+2 = dev_spin_eo2_up //
  // DUM_MATRIX+3 = dev_spin_eo2_dn //
  //                                //
  ////////////////////////////////////
  
  
  
  
  ///////////////////////////////////
  // INITIALIZATIONS & ASSIGNMENTS //	// have to use (one) other auxiliary field(s) than the calling function dev_cg_eo_nd
  ///////////////////////////////////
  
  dev_spin_eo2_up = spinout_up;		// need no memory allocated
  dev_spin_eo2_dn = spinout_dn;
  												///////////// THEORY ////////////////////////////////////////////////////////////////
  												//                                                                                 //
  												//  (Q_tilde) = gamma5 * ((M_oo) - (M_oe)(Mee^-1)(M_eo))                           //
  												//  (Q_tilde)(Q_tilde_dagger) * (up,dn) = (Q_tilde) * (b,a)                        //
  ///////////////										//                                                    (a,b) = (Q_tilde) * (dn,up)  //
  // MAIN BODY //										//                                                                                 //
  ///////////////										/////////////////////////////////////////////////////////////////////////////////////
  
  
  double nrm = 1.0 / (1.0 + g_mubar*g_mubar - g_epsbar*g_epsbar);
  
  
  if (g_proc_id == 0)
    printf("This is matrix_mpi_debug4(). ");
  
  
  spinor * help1  =  (spinor *) malloc((VOLUME+RAND) * 6*sizeof(dev_spinor));			// "/2" is missing because of "double"
  spinor * help2  =  (spinor *) malloc((VOLUME+RAND) * 6*sizeof(dev_spinor));
  spinor * help3  =  (spinor *) malloc((VOLUME+RAND) * 6*sizeof(dev_spinor));
  spinor * help4  =  (spinor *) malloc((VOLUME+RAND) * 6*sizeof(dev_spinor));

  
  
  ///////////////////////////////////////							/////////////////////////////////
  //        Q_tilde_dagger(2x2)        //							// (a,b) = (Q_tilde) * (dn,up) //
  ///////////////////////////////////////							/////////////////////////////////
  
  
  ////////////////////////////
  // (M_oe) (Mee^-1) (M_eo) //
  ////////////////////////////
  
  
  xchange_test = 1;		// is set to "1", the xchange_field() in the host Hopping_Matrix() will be turned off
  				// procedure:	define global variable in Hopping_Matrix.c which can disable the xchange of the fields
  				//		set this variable appropriately ... 
  
  #ifndef CHECK_HOPPING_MATRIX

  	  #ifndef HOPPING_DEBUG
  	    xchange_field_wrapper(spinin_dn, 0);
  	    dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, spinin_dn, dev_spin_eo1_up, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0);	// dev_spin_eo1_up = (M_eo) * spinin_dn
  	  #else
  	    //printf("Does it work?\n");
  	    xchange_field_wrapper(spinin_dn, 0);
  	    Hopping_Matrix_wrapper(0, dev_spin_eo1_up, spinin_dn);
  	  #endif

  	  #ifndef HOPPING_DEBUG
  	    xchange_field_wrapper(spinin_up, 0);
  	    dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, spinin_up, dev_spin_eo1_dn, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0);	// dev_spin_eo1_dn = (M_eo) * spinin_up
  	  #else
  	    xchange_field_wrapper(spinin_up, 0);
  	    Hopping_Matrix_wrapper(0, dev_spin_eo1_dn, spinin_up);
  	  #endif
  #else
    	to_host(help1, spinin_dn, h2d_spin_up, dev_spinsize);
    	Hopping_Matrix(EO, help2, help1);				// g_spinor_field[DUM_MATRIX]   = (M_eo) * k_charm
    	to_device(dev_spin_eo1_up, help2, h2d_spin_up, dev_spinsize);
    	
    	to_host(help3, spinin_up, h2d_spin_dn, dev_spinsize);
    	Hopping_Matrix(EO, help4, help3);				// g_spinor_field[DUM_MATRIX+1] = (M_eo) * k_strange
    	to_device(dev_spin_eo1_dn, help4, h2d_spin_dn, dev_spinsize);
  #endif
  
  
  // written:
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo1_up, dev_spin_eo2_up, -1.0);		// dev_spin_eo2_up  =  (1 - imubar) * dev_spin_eo1_up  =  (1 - imubar)*(M_eo) * spinin_dn
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo1_dn, dev_spin_eo2_dn, +1.0);		// dev_spin_eo2_dn  =  (1 + imubar) * dev_spin_eo1_dn  =  (1 + imubar)*(M_eo) * spinin_up
  
  
  // CUBLAS:
  cublasSaxpy (N_floats, g_epsbar, (float *) dev_spin_eo1_dn, 1, (float *) dev_spin_eo2_up, 1);
  										// dev_spin_eo2_up  =  dev_spin_eo2_up  +  epsbar * dev_spin_eo1_dn  =  (1 - imubar)*(M_eo) * spinin_dn  +  epsbar * (M_eo) * spinin_up
  cublasSaxpy (N_floats, g_epsbar, (float *) dev_spin_eo1_up, 1, (float *) dev_spin_eo2_dn, 1);
  										// dev_spin_eo2_dn  =  dev_spin_eo2_dn  +  epsbar * dev_spin_eo1_up  =  (1 + imubar)*(M_eo) * spinin_up  +  epsbar * (M_eo) * spinin_dn
  
  // CUBLAS:
  cublasSscal (N_floats, nrm, (float *) dev_spin_eo2_up, 1);			// dev_spin_eo2_up  =  nrm * dev_spin_eo2_up  =  nrm*(1-imubar)*(M_eo) * spinin_dn  +  nrm*epsbar*(M_eo) * spinin_up
  
  cublasSscal (N_floats, nrm, (float *) dev_spin_eo2_dn, 1);			// dev_spin_eo2_dn  =  nrm * dev_spin_eo2_dn  =  nrm*(1+imubar)*(M_eo) * spinin_up  +  nrm*epsbar*(M_eo) * spinin_dn
  
  
  #ifndef CHECK_HOPPING_MATRIX

  	  #ifndef HOPPING_DEBUG
  	    xchange_field_wrapper(dev_spin_eo2_up, 1);
  	    dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo2_up, dev_spin_eo1_up, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1);
  	  #else
  	    xchange_field_wrapper(dev_spin_eo2_up, 1);
  	    Hopping_Matrix_wrapper(1, dev_spin_eo1_up, dev_spin_eo2_up);
  	  #endif							

  	  #ifndef HOPPING_DEBUG
  	    xchange_field_wrapper(dev_spin_eo2_dn, 1);
  	    dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo2_dn, dev_spin_eo1_dn, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1);
  	  #else
  	    xchange_field_wrapper(dev_spin_eo2_dn, 1);
  	    Hopping_Matrix_wrapper(1, dev_spin_eo1_dn, dev_spin_eo2_dn);
  	  #endif
  #else
    	to_host(help1, dev_spin_eo2_up, h2d_spin_up, dev_spinsize);
    	Hopping_Matrix(OE, help2, help1);		// g_spinor_field[DUM_MATRIX]   = (M_oe) * g_spinor_field[DUM_MATRIX+2]
    	to_device(dev_spin_eo1_up, help2, h2d_spin_up, dev_spinsize);
    	
    	to_host(help3, dev_spin_eo2_dn, h2d_spin_dn, dev_spinsize);
    	Hopping_Matrix(OE, help4, help3);		// g_spinor_field[DUM_MATRIX+1] = (M_oe) * g_spinor_field[DUM_MATRIX+3]
    	to_device(dev_spin_eo1_dn, help4, h2d_spin_dn, dev_spinsize);
  #endif
  
  
  
  ////////////
  // (M_oo) //
  ////////////
  
  
  // written:
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(spinin_dn, dev_spin_eo2_up, +1.0);			// dev_spin_eo2_up = (1 + imubar) * spinin_dn
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(spinin_up, dev_spin_eo2_dn, -1.0);			// dev_spin_eo2_dn = (1 - imubar) * spinin_up
  
  
  // CUBLAS:													// remember: this is (M_oo) * (spinin_dn, spinin_up):
  cublasSaxpy (N_floats, -g_epsbar, (float *) spinin_up, 1, (float *) dev_spin_eo2_up, 1);
  												// dev_spin_eo2_up  =  dev_spin_eo2_up - epsbar*spinin_up  =  (1+imubar)*spinin_dn - epsbar*spinin_up
  cublasSaxpy (N_floats, -g_epsbar, (float *) spinin_dn, 1, (float *) dev_spin_eo2_dn, 1);
  												// dev_spin_eo2_dn  =  dev_spin_eo2_dn - epsbar*spinin_dn  =  (1-imubar)*spinin_up - epsbar*spinin_dn
  
  
  
  ///////////////////////////////////////
  // (M_oo)  -  (M_oe) (Mee^-1) (M_eo) //
  ///////////////////////////////////////
  
  // CUBLAS:													// this is ((M_oo) - (M_oe)(Mee^-1)(M_eo)) * (spinin_dn, spinin_up):
  cublasSaxpy (N_floats, -1.0, (float *) dev_spin_eo1_up, 1, (float *) dev_spin_eo2_up, 1);
  							// dev_spin_eo2_up  =  dev_spin_eo2_up  -  dev_spin_eo1_up  =                      (1+imubar) * spinin_dn  -                    epsbar * spinin_up
  							//                                                             - (M_oe)*nrm*(1-imubar)*(M_eo) * spinin_dn  -  (M_oe)*nrm*epsbar*(M_eo) * spinin_up
  cublasSaxpy (N_floats, -1.0, (float *) dev_spin_eo1_dn, 1, (float *) dev_spin_eo2_dn, 1);
  							// dev_spin_eo2_dn  =  dev_spin_eo2_dn  -  dev_spin_eo1_dn  =                      (1-imubar) * spinin_up  -                    epsbar * spinin_dn
  							//                                                             - (M_oe)*nrm*(1+imubar)*(M_eo) * spinin_up  -  (M_oe)*nrm*epsbar*(M_eo) * spinin_dn
  
  
  ////////////
  // gamma5 //
  ////////////
  
  // Flo:
  dev_gamma5<<<gridsize3, blocksize3>>>(dev_spin_eo2_up, dev_spin_eo2_up);					// dev_spin_eo2_up = gamma5 * dev_spin_eo2_up 
  dev_gamma5<<<gridsize3, blocksize3>>>(dev_spin_eo2_dn, dev_spin_eo2_dn);					// dev_spin_eo2_dn = gamma5 * dev_spin_eo2_dn
  	
  
  
  
  
  
  ////////////////////
  // (a,b) -> (b,a) //
  ////////////////////
  dev_copy_spinor_field<<<gridsize4, blocksize4>>>(dev_spin_eo2_dn, dev_spin_eo3_up);			// dev_spin_eo3_up = dev_spin_eo2_dn
  dev_copy_spinor_field<<<gridsize4, blocksize4>>>(dev_spin_eo2_up, dev_spin_eo3_dn);			// dev_spin_eo3_dn = dev_spin_eo2_up
  
  
  
  
  
  
  ///////////////////////////////////								///////////////////////
  //        Q_tilde(2x2)           //								// (Q_tilde) * (b,a) //
  ///////////////////////////////////								///////////////////////
  
  
  ////////////////////////////
  // (M_oe) (Mee^-1) (M_eo) //
  ////////////////////////////
  
  #ifndef CHECK_HOPPING_MATRIX

  	  #ifndef HOPPING_DEBUG
  	    xchange_field_wrapper(dev_spin_eo3_up, 0);
  	    dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo3_up, dev_spin_eo1_up, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0);	// dev_spin_eo1_up = (M_eo) * dev_spin_eo3_up
  	  #else
  	    xchange_field_wrapper(dev_spin_eo3_up, 0);
  	    Hopping_Matrix_wrapper(0, dev_spin_eo1_up, dev_spin_eo3_up);
  	  #endif

  	  #ifndef HOPPING_DEBUG
  	    xchange_field_wrapper(dev_spin_eo3_dn, 0);
  	    dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo3_dn, dev_spin_eo1_dn, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0);	// dev_spin_eo1_dn = (M_eo) * dev_spin_eo3_dn
  	  #else
  	    xchange_field_wrapper(dev_spin_eo3_dn, 0);
  	    Hopping_Matrix_wrapper(0, dev_spin_eo1_dn, dev_spin_eo3_dn);
  	  #endif
  #else
    	to_host(help1, dev_spin_eo3_up, h2d_spin_up, dev_spinsize);
    	Hopping_Matrix(EO, help2, help1);	// g_spinor_field[DUM_MATRIX]   = (M_eo) * g_spinor_field[DUM_MATRIX+7]
    	to_device(dev_spin_eo1_up, help2, h2d_spin_up, dev_spinsize);
    	
    	to_host(help3, dev_spin_eo3_dn, h2d_spin_dn, dev_spinsize);
    	Hopping_Matrix(EO, help4, help3);	// g_spinor_field[DUM_MATRIX+1] = (M_eo) * g_spinor_field[DUM_MATRIX+6]
    	to_device(dev_spin_eo1_dn, help4, h2d_spin_dn, dev_spinsize);
  #endif
  
  
  // written:
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo1_up, dev_spin_eo2_up, -1.0);			// dev_spin_eo2_up  =  (1 - imubar) * dev_spin_eo1_up  =  (1 - imubar)*(M_eo) * dev_spin_eo3_up
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo1_dn, dev_spin_eo2_dn, +1.0);			// dev_spin_eo2_dn  =  (1 + imubar) * dev_spin_eo1_dn  =  (1 + imubar)*(M_eo) * dev_spin_eo3_dn
  
  
  // CUBLAS:
  cublasSaxpy (N_floats, g_epsbar, (float *) dev_spin_eo1_dn, 1, (float *) dev_spin_eo2_up, 1);
  										// dev_spin_eo2_up  =  dev_spin_eo2_up  +  epsbar * dev_spin_eo1_dn  =  (1 - imubar)*(M_eo) * dev_spin_eo3_up  +  epsbar * (M_eo) * dev_spin_eo3_dn
  cublasSaxpy (N_floats, g_epsbar, (float *) dev_spin_eo1_up, 1, (float *) dev_spin_eo2_dn, 1);
  										// dev_spin_eo2_dn  =  dev_spin_eo2_dn  +  epsbar * dev_spin_eo1_up  =  (1 + imubar)*(M_eo) * dev_spin_eo3_dn  +  epsbar * (M_eo) * dev_spin_eo3_up
  
  // CUBLAS:
  cublasSscal (N_floats, nrm, (float *) dev_spin_eo2_up, 1);			// dev_spin_eo2_up  =  nrm * dev_spin_eo2_up  =  nrm*(1-imubar)*(M_eo) * dev_spin_eo3_up  +  nrm*epsbar*(M_eo) * dev_spin_eo3_dn
  
  cublasSscal (N_floats, nrm, (float *) dev_spin_eo2_dn, 1);			// dev_spin_eo2_dn  =  nrm * dev_spin_eo2_dn  =  nrm*(1+imubar)*(M_eo) * dev_spin_eo3_dn  +  nrm*epsbar*(M_eo) * dev_spin_eo3_up
  
  
  #ifndef CHECK_HOPPING_MATRIX

  	  #ifndef HOPPING_DEBUG
  	    xchange_field_wrapper(dev_spin_eo2_up, 1);
  	    dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo2_up, dev_spin_eo1_up, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1);
  	  #else
  	    xchange_field_wrapper(dev_spin_eo2_up, 1);
  	    Hopping_Matrix_wrapper(1, dev_spin_eo1_up, dev_spin_eo2_up);
  	  #endif

  	  #ifndef HOPPING_DEBUG
  	    xchange_field_wrapper(dev_spin_eo2_dn, 1);
  	    dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo2_dn, dev_spin_eo1_dn, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1);
  	  #else
  	    xchange_field_wrapper(dev_spin_eo2_dn, 1);
  	    Hopping_Matrix_wrapper(1, dev_spin_eo1_dn, dev_spin_eo2_dn);
  	  #endif
  #else
    	to_host(help1, dev_spin_eo2_up, h2d_spin_up, dev_spinsize);
    	Hopping_Matrix(OE, help2, help1);				// l_strange = (M_oe) * g_spinor_field[DUM_MATRIX+2]
    	to_device(dev_spin_eo1_up, help2, h2d_spin_up, dev_spinsize);
    	
    	to_host(help3, dev_spin_eo2_dn, h2d_spin_dn, dev_spinsize);
    	Hopping_Matrix(OE, help4, help3);				// l_charm   = (M_oe) * g_spinor_field[DUM_MATRIX+3]
    	to_device(dev_spin_eo1_dn, help4, h2d_spin_dn, dev_spinsize);
  #endif
  
  
  ////////////
  // (M_oo) //
  ////////////
  
  // written:
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo3_up, dev_spin_eo2_up, +1.0);				// dev_spin_eo2_up = (1 + imubar) * dev_spin_eo3_up
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo3_dn, dev_spin_eo2_dn, -1.0);				// dev_spin_eo2_dn = (1 - imubar) * dev_spin_eo3_dn
  
  
  // CUBLAS:											// remember: this is (M_oo) * (dev_spin_eo3_up, dev_spin_eo3_dn):
  cublasSaxpy (N_floats, -g_epsbar, (float *) dev_spin_eo3_dn, 1, (float *) dev_spin_eo2_up, 1);
  												// dev_spin_eo2_up  =  dev_spin_eo2_up - epsbar*dev_spin_eo3_dn  =  (1+imubar)*dev_spin_eo3_up - epsbar*dev_spin_eo3_dn
  cublasSaxpy (N_floats, -g_epsbar, (float *) dev_spin_eo3_up, 1, (float *) dev_spin_eo2_dn, 1);
  												// dev_spin_eo2_dn  =  dev_spin_eo2_dn - epsbar*dev_spin_eo3_up  =  (1-imubar)*dev_spin_eo3_dn - epsbar*dev_spin_eo3_up
  
  
  
  ///////////////////////////////////////
  // (M_oo)  -  (M_oe) (Mee^-1) (M_eo) //
  ///////////////////////////////////////
  
  // CUBLAS:											// this is ( (M_oo) - (M_oe) (Mee^-1) (M_eo) ) * (dev_spin_eo3_up, dev_spin_eo3_dn)
  cublasSaxpy (N_floats, -1.0, (float *) dev_spin_eo1_up, 1, (float *) dev_spin_eo2_up, 1);
  							// dev_spin_eo2_up  =  dev_spin_eo2_up  -  dev_spin_eo1_up  =                      (1+imubar) * dev_spin_eo3_up  -                    epsbar * dev_spin_eo3_dn
  							//                                                             - (M_oe)*nrm*(1-imubar)*(M_eo) * dev_spin_eo3_up  -  (M_oe)*nrm*epsbar*(M_eo) * dev_spin_eo3_dn
  cublasSaxpy (N_floats, -1.0, (float *) dev_spin_eo1_dn, 1, (float *) dev_spin_eo2_dn, 1);
  							// dev_spin_eo2_dn  =  dev_spin_eo2_dn  -  dev_spin_eo1_dn  =                      (1-imubar) * dev_spin_eo3_dn  -                    epsbar * dev_spin_eo3_up
  							//                                                             - (M_oe)*nrm*(1+imubar)*(M_eo) * dev_spin_eo3_dn  -  (M_oe)*nrm*epsbar*(M_eo) * dev_spin_eo3_up
  
  
  ////////////
  // gamma5 //
  ////////////
  
  // Flo:
  dev_gamma5<<<gridsize3, blocksize3>>>(dev_spin_eo2_up, dev_spin_eo2_up);					// dev_spin_eo2_up = gamma5 * dev_spin_eo2_up 
  dev_gamma5<<<gridsize3, blocksize3>>>(dev_spin_eo2_dn, dev_spin_eo2_dn);					// dev_spin_eo2_dn = gamma5 * dev_spin_eo2_dn
  
  
  
  
  /*
  ////////////
  // output //										// output is already done by setting  dev_spin_eo2_up/dn = spinout_up/dn
  ////////////
  
  dev_copy_spinor_field<<<griddim2, blockdim2 >>>(dev_spin_eo2_up, spinout_up);		// spinout_up = dev_spin_eo2_up
  dev_copy_spinor_field<<<griddim2, blockdim2 >>>(dev_spin_eo2_dn, spinout_dn);		// spinout_dn = dev_spin_eo2_dn
  */
  
  
  xchange_test = 0;
  
  
  return;
  
}//matrix_mpi_debug4()









// this replaces Q_Qdagger_ND() for debugging purposes in the parallel case !!

// first debugging RESULT:
//	the host Hopping_Matrix() works for the halfspinor/gaugecopy-option enabled or disabled !!
//	contrary to the use of Hopping_Matrix() in the original code
//	fatal error or what?! :D
//	the error has to concern the fields (lengths, usage, order etc.)

// the Hopping_Matrix on device is not implemented here

void matrix_mpi_debug10 (spinor * const l_strange, spinor * const l_charm,	// output
                        spinor * const k_strange, spinor * const k_charm) {	// input
  
  
  int N_sites  =    VOLUME/2;		// #lattice sites
  int N_floats = 24*VOLUME/2;		// #floats
  
  // additional for the debugging purposes
  size_t dev_spinsize = 6*(VOLUME+RAND)/2*sizeof(dev_spinor);
  
  
  int gridsize;			// auxiliary
  int blocksize;		// auxiliary
  
  blocksize = 128;
  int blocksize1 = blocksize;					// here:	dev_zero_spinor_field , dev_copy_spinor_field
  int gridsize1  = (int) (VOLUME/2/blocksize) + 1;
  
  blocksize = 128;
  int blocksize2 = blocksize;					// passed:	dev_Hopping_Matrix
  int gridsize2  = (int) (VOLUME/2/blocksize) + 1;
  
  blocksize = 128;
  int blocksize3 = blocksize;					// passed:	dev_mul_one_pm_imubar_gamma5
  int gridsize3  = (int) (VOLUME/2/blocksize) + 1;
  
  blocksize = 128;
  int blocksize4 = blocksize;					// passed:	dev_gamma5
  int gridsize4  = (int) (VOLUME/2/blocksize) + 1;
  
  blocksize = 128;
  int blocksize5 = blocksize;					// passed:	dev_copy_spinor_field
  int gridsize5  = (int) (VOLUME/2/blocksize) + 1;
  
  dev_spinor * spinin_up;
  dev_spinor * spinin_dn;
  dev_spinor * spinout_up;
  dev_spinor * spinout_dn;
  
  
  cudaMalloc((void **) &dev_spin_eo2_up, dev_spinsize);
  cudaMalloc((void **) &dev_spin_eo2_dn, dev_spinsize);
  cudaMalloc((void **) &spinin_up, dev_spinsize);
  cudaMalloc((void **) &spinin_dn, dev_spinsize);
  cudaMalloc((void **) &spinout_up, dev_spinsize);
  cudaMalloc((void **) &spinout_dn, dev_spinsize);
  
  
  
  double nrm = 1./(1. + g_mubar*g_mubar - g_epsbar*g_epsbar);				// nrm = (1 + mubar^2 - epsbar^2)^-1
  
  
  
  
  if (g_proc_id == 0)
    printf("This is matrix_mpi_debug10(). ");
  
  


  /* FIRST THE  Qhat(2x2)^dagger  PART */						// we will apply Qhat(2x2) with charme and strange interchanged
  											// which is equivalent to apply Qhat(2x2)^dagger
  
  /* Here the  M_oe Mee^-1 M_eo  implementation  */
  
  #ifndef CHECK_HOPPING_MATRIX
    Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX]  , k_charm);				// g_spinor_field[DUM_MATRIX]   = (M_eo) * k_charm				// notice the order
    Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX+1], k_strange);			// g_spinor_field[DUM_MATRIX+1] = (M_eo) * k_strange				//   of k_charm and k_strange
  #else
    				// xchange
  				#ifndef HOPPING_DEBUG
  				  xchange_field_wrapper(spinin_dn, 0);
  				#endif
  		  #ifndef HOPPING_DEBUG
  		    dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, spinin_dn, dev_spin_eo1_up, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0);	// dev_spin_eo1_up = (M_eo) * spinin_dn
  		  #else
  		    Hopping_Matrix_wrapper(0, dev_spin_eo1_up, spinin_dn);
  		  #endif
  				// xchange
  				#ifndef HOPPING_DEBUG
  				  xchange_field_wrapper(spinin_up, 0);
  				#endif
  		  #ifndef HOPPING_DEBUG
  		    dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, spinin_up, dev_spin_eo1_dn, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0);	// dev_spin_eo1_dn = (M_eo) * spinin_up
  		  #else
  		    Hopping_Matrix_wrapper(0, dev_spin_eo1_dn, spinin_up);
  		  #endif
  #endif




  #ifndef CHECK_IMUGAMMA5
											//				remark: here the factor  GAMMA5  is not written:
    mul_one_minus_imubar(g_spinor_field[DUM_MATRIX+2], g_spinor_field[DUM_MATRIX]);	// g_spinor_field[DUM_MATRIX+2] = (1 - imubar)*(M_eo) * g_spinor_field[DUM_MATRIX]
    											//                              = (1 - imubar)*(M_eo) * k_charm
    mul_one_plus_imubar (g_spinor_field[DUM_MATRIX+3], g_spinor_field[DUM_MATRIX+1]);	// g_spinor_field[DUM_MATRIX+3] = (1 + imubar)*(M_eo) * g_spinor_field[DUM_MATRIX+1]
  											//                              = (1 + imubar)*(M_eo) * k_strange
  #else
    to_device(dev_spin_eo1_up, g_spinor_field[DUM_MATRIX], h2d_spin_up, dev_spinsize);
      dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo1_up, dev_spin_eo2_up, -1.0);		// dev_spin_eo2_up  =  (1 - imubar) * dev_spin_eo1_up
    to_host(g_spinor_field[DUM_MATRIX+2], dev_spin_eo2_up, h2d_spin_up, dev_spinsize);
    
    to_device(dev_spin_eo1_dn, g_spinor_field[DUM_MATRIX+1], h2d_spin_dn, dev_spinsize);
      dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo1_dn, dev_spin_eo2_dn, +1.0);		// dev_spin_eo2_dn  =  (1 + imubar) * dev_spin_eo1_dn
    to_host(g_spinor_field[DUM_MATRIX+3], dev_spin_eo2_dn, h2d_spin_dn, dev_spinsize);
  #endif
  
  
  
  
  #ifndef CHECK_CUBLAS1
    assign_add_mul_r(g_spinor_field[DUM_MATRIX+2], g_spinor_field[DUM_MATRIX+1], g_epsbar, VOLUME/2);
    											// g_spinor_field[DUM_MATRIX+2] = g_spinor_field[DUM_MATRIX+2]     +  epsbar*g_spinor_field[DUM_MATRIX+1]
    											//                              = (1 - imubar)*(M_eo) * k_charm    +  epsbar*(M_eo) * k_strange
    assign_add_mul_r(g_spinor_field[DUM_MATRIX+3], g_spinor_field[DUM_MATRIX]  , g_epsbar, VOLUME/2);
    											// g_spinor_field[DUM_MATRIX+3] = g_spinor_field[DUM_MATRIX+3]     +  epsbar*g_spinor_field[DUM_MATRIX]
    											//                              = (1 + imubar)*(M_eo) * k_strange  +  epsbar*(M_eo) * k_charm
    
    mul_r(g_spinor_field[DUM_MATRIX+2], nrm, g_spinor_field[DUM_MATRIX+2], VOLUME/2);	// g_spinor_field[DUM_MATRIX+2] = nrm * g_spinor_field[DUM_MATRIX+2]
    											//                              = nrm * ( (1 - imubar)*(M_eo) * k_charm    +  epsbar*(M_eo) * k_strange )
    											//                              = nrm*(1 - imubar)*(M_eo)*k_charm    +  nrm*epsbar*(M_eo)*k_strange
    											
    mul_r(g_spinor_field[DUM_MATRIX+3], nrm, g_spinor_field[DUM_MATRIX+3], VOLUME/2);	// g_spinor_field[DUM_MATRIX+3] = nrm * g_spinor_field[DUM_MATRIX+3]
    											//                              = nrm * ( (1 + imubar)*(M_eo) * k_strange  +  epsbar*(M_eo) * k_charm )
    											//                              = nrm*(1 + imubar)*(M_eo)*k_strange  +  nrm*epsbar*(M_eo)*k_charm
  #else
    to_device(dev_spin_eo1_dn, g_spinor_field[DUM_MATRIX+1], h2d_spin_up, dev_spinsize);
    to_device(dev_spin_eo2_up, g_spinor_field[DUM_MATRIX+2], h2d_spin_up, dev_spinsize);
      cublasSaxpy (N_floats, g_epsbar, (float *) dev_spin_eo1_dn, 1, (float *) dev_spin_eo2_up, 1);
    										// dev_spin_eo2_up  =  dev_spin_eo2_up  +  epsbar * dev_spin_eo1_dn
      cublasSscal (N_floats, nrm, (float *) dev_spin_eo2_up, 1);		// dev_spin_eo2_up  =  nrm * dev_spin_eo2_up
    to_host(g_spinor_field[DUM_MATRIX+2], dev_spin_eo2_up, h2d_spin_up, dev_spinsize);
    
    
    to_device(dev_spin_eo1_up, g_spinor_field[DUM_MATRIX], h2d_spin_dn, dev_spinsize);
    to_device(dev_spin_eo2_dn, g_spinor_field[DUM_MATRIX+3], h2d_spin_dn, dev_spinsize);
      cublasSaxpy (N_floats, g_epsbar, (float *) dev_spin_eo1_up, 1, (float *) dev_spin_eo2_dn, 1);
    										// dev_spin_eo2_dn  =  dev_spin_eo2_dn  +  epsbar * dev_spin_eo1_up
      cublasSscal (N_floats, nrm, (float *) dev_spin_eo2_dn, 1);		// dev_spin_eo2_dn  =  nrm * dev_spin_eo2_dn
    to_host(g_spinor_field[DUM_MATRIX+3], dev_spin_eo2_dn, h2d_spin_dn, dev_spinsize);
  #endif
  
  
  
  
  #ifndef CHECK_HOPPING_MATRIX
    Hopping_Matrix(OE, g_spinor_field[DUM_MATRIX]  , g_spinor_field[DUM_MATRIX+2]);	// g_spinor_field[DUM_MATRIX]   = (M_oe) * g_spinor_field[DUM_MATRIX+2]
    											//                              = (M_oe)*nrm*(1 - imubar)*(M_eo)*k_charm  +  (M_oe)*nrm*epsbar*(M_eo)*k_strange
    Hopping_Matrix(OE, g_spinor_field[DUM_MATRIX+1], g_spinor_field[DUM_MATRIX+3]);	// g_spinor_field[DUM_MATRIX+1] = (M_oe) * g_spinor_field[DUM_MATRIX+3]
    											//                              = (M_oe)*nrm*(1 + imubar)*(M_eo)*k_strange  +  (M_oe)*nrm*epsbar*(M_eo) * k_charm
  #else
  				// xchange
  				#ifndef HOPPING_DEBUG
  				  xchange_field_wrapper(dev_spin_eo2_up, 1);
  				#endif
  		  #ifndef HOPPING_DEBUG
  		    dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo2_up, dev_spin_eo1_up, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1);
  		  #else
  		    Hopping_Matrix_wrapper(1, dev_spin_eo1_up, dev_spin_eo2_up);
  		  #endif								
  				// xchange
  				#ifndef HOPPING_DEBUG
  				  xchange_field_wrapper(dev_spin_eo2_dn, 1);
  				#endif
  		  #ifndef HOPPING_DEBUG
  		    dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo2_dn, dev_spin_eo1_dn, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1);
  		  #else
  		    Hopping_Matrix_wrapper(1, dev_spin_eo1_dn, dev_spin_eo2_dn);
  		  #endif
  #endif
  
  
  
  
  /* Here the M_oo  implementation  */
  
  #ifndef CHECK_IMUGAMMA5
    mul_one_plus_imubar (g_spinor_field[DUM_MATRIX+2], k_charm);			// g_spinor_field[DUM_MATRIX+2] = (1 + imubar) * k_charm
    mul_one_minus_imubar(g_spinor_field[DUM_MATRIX+3], k_strange);			// g_spinor_field[DUM_MATRIX+3] = (1 - imubar) * k_strange
  #else
    to_device(spinin_dn, k_charm, h2d_spin_up, dev_spinsize);
      dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(spinin_dn, dev_spin_eo2_up, +1.0);			// dev_spin_eo2_up = (1 + imubar) * spinin_dn
    to_host(g_spinor_field[DUM_MATRIX+2], dev_spin_eo2_up, h2d_spin_up, dev_spinsize);
    
    to_device(spinin_up, k_strange, h2d_spin_dn, dev_spinsize);
      dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(spinin_up, dev_spin_eo2_dn, -1.0);			// dev_spin_eo2_dn = (1 - imubar) * spinin_up
    to_host(g_spinor_field[DUM_MATRIX+3], dev_spin_eo2_dn, h2d_spin_dn, dev_spinsize);
  #endif




  #ifndef CHECK_CUBLAS2
    assign_add_mul_r(g_spinor_field[DUM_MATRIX+2], k_strange, -g_epsbar, VOLUME/2);	// g_spinor_field[DUM_MATRIX+2] = g_spinor_field[DUM_MATRIX+2] - epsbar*k_strange
    											//                              = (1 + imubar) * k_charm       - epsbar*k_strange
    assign_add_mul_r(g_spinor_field[DUM_MATRIX+3], k_charm  , -g_epsbar, VOLUME/2);	// g_spinor_field[DUM_MATRIX+3] = g_spinor_field[DUM_MATRIX+3] - epsbar*k_charm
    											//                              = (1 - imubar) * k_strange     - epsbar*k_charm
  #else
    to_device(spinin_up, k_strange, h2d_spin_up, dev_spinsize);
    to_device(dev_spin_eo2_up, g_spinor_field[DUM_MATRIX+2], h2d_spin_up, dev_spinsize);
      cublasSaxpy (N_floats, -g_epsbar, (float *) spinin_up, 1, (float *) dev_spin_eo2_up, 1);	// dev_spin_eo2_up  =  dev_spin_eo2_up - epsbar*spinin_up
    to_host(g_spinor_field[DUM_MATRIX+2], dev_spin_eo2_up, h2d_spin_up, dev_spinsize);
    
    to_device(spinin_dn, k_charm, h2d_spin_dn, dev_spinsize);
    to_device(dev_spin_eo2_dn, g_spinor_field[DUM_MATRIX+3], h2d_spin_dn, dev_spinsize);
      cublasSaxpy (N_floats, -g_epsbar, (float *) spinin_dn, 1, (float *) dev_spin_eo2_dn, 1);	// dev_spin_eo2_dn  =  dev_spin_eo2_dn - epsbar*spinin_dn
    to_host(g_spinor_field[DUM_MATRIX+3], dev_spin_eo2_dn, h2d_spin_dn, dev_spinsize);
  #endif
  
  
  
  
  // here the (M_oo - M_oe Mee^-1 M_eo) implementation
  
  #ifndef CHECK_CUBLAS3
    diff(g_spinor_field[DUM_MATRIX+4], g_spinor_field[DUM_MATRIX+2], g_spinor_field[DUM_MATRIX]  , VOLUME/2);
    											// g_spinor_field[DUM_MATRIX+4] = g_spinor_field[DUM_MATRIX+2] - g_spinor_field[DUM_MATRIX]
    											//                              =                          (1 + imubar) * k_charm  -                    epsbar * k_strange
    											//                                     - (M_oe)*nrm*(1 - imubar)*(M_eo) * k_charm  -  (M_oe)*nrm*epsbar*(M_eo) * k_strange
    diff(g_spinor_field[DUM_MATRIX+5], g_spinor_field[DUM_MATRIX+3], g_spinor_field[DUM_MATRIX+1], VOLUME/2);
    											// g_spinor_field[DUM_MATRIX+5] = g_spinor_field[DUM_MATRIX+3] - g_spinor_field[DUM_MATRIX+1]
    											//                              =                          (1 - imubar) * k_strange  -                    epsbar * k_charm
    											//                                     - (M_oe)*nrm*(1 + imubar)*(M_eo) * k_strange  -  (M_oe)*nrm*epsbar*(M_eo) * k_charm
  #else
    to_device(dev_spin_eo1_up, g_spinor_field[DUM_MATRIX], h2d_spin_up, dev_spinsize);
    to_device(dev_spin_eo2_up, g_spinor_field[DUM_MATRIX+2], h2d_spin_up, dev_spinsize);
      cublasSaxpy (N_floats, -1.0, (float *) dev_spin_eo1_up, 1, (float *) dev_spin_eo2_up, 1);	// dev_spin_eo2_up  =  dev_spin_eo2_up  -  dev_spin_eo1_up
    to_host(g_spinor_field[DUM_MATRIX+4], dev_spin_eo2_up, h2d_spin_up, dev_spinsize);
    
    
    to_device(dev_spin_eo1_dn, g_spinor_field[DUM_MATRIX+1], h2d_spin_dn, dev_spinsize);
    to_device(dev_spin_eo2_dn, g_spinor_field[DUM_MATRIX+3], h2d_spin_dn, dev_spinsize);
      cublasSaxpy (N_floats, -1.0, (float *) dev_spin_eo1_dn, 1, (float *) dev_spin_eo2_dn, 1);	// dev_spin_eo2_dn  =  dev_spin_eo2_dn  -  dev_spin_eo1_dn
    to_host(g_spinor_field[DUM_MATRIX+5], dev_spin_eo2_dn, h2d_spin_dn, dev_spinsize);
  #endif
  
  
  

  /* and finally the  GAMMA5  multiplication  */
  
  #ifndef CHECK_GAMMA5
    gamma5(g_spinor_field[DUM_MATRIX+2], g_spinor_field[DUM_MATRIX+4], VOLUME/2);		// g_spinor_field[DUM_MATRIX+2] = gamma5 * g_spinor_field[DUM_MATRIX+4] ?!= l_charm'
    gamma5(g_spinor_field[DUM_MATRIX+3], g_spinor_field[DUM_MATRIX+5], VOLUME/2);		// g_spinor_field[DUM_MATRIX+3] = gamma5 * g_spinor_field[DUM_MATRIX+5] ?!= l_strange'
  #else
    to_device(dev_spin_eo2_up, g_spinor_field[DUM_MATRIX+4], h2d_spin_up, dev_spinsize);
      dev_gamma5<<<gridsize3, blocksize3>>>(dev_spin_eo2_up, dev_spin_eo2_up);					// dev_spin_eo2_up = gamma5 * dev_spin_eo2_up 
    to_host(g_spinor_field[DUM_MATRIX+2], dev_spin_eo2_up, h2d_spin_up, dev_spinsize);
    
    to_device(dev_spin_eo2_dn, g_spinor_field[DUM_MATRIX+5], h2d_spin_dn, dev_spinsize);
      dev_gamma5<<<gridsize3, blocksize3>>>(dev_spin_eo2_dn, dev_spin_eo2_dn);					// dev_spin_eo2_dn = gamma5 * dev_spin_eo2_dn
    to_host(g_spinor_field[DUM_MATRIX+3], dev_spin_eo2_dn, h2d_spin_dn, dev_spinsize);
  #endif
  
  
  

  /* The normalisation by the max. eigenvalue is done twice at the end */		// what ??


  /* We have to reassigin as follows to avoid overwriting */
  /* Recall in fact that   Q^hat = tau_1 Q tau_1  , hence  */

  /*  ABOVE: dum_matrix+2  is  l_charm   goes to  dum_matrix+6 :BELOW */
  /*  ABOVE: dum_matrix+3  is  l_strange   goes to  dum_matrix+7 :BELOW */
  
  
  
  
  #ifndef CHECK_COPY
    assign(g_spinor_field[DUM_MATRIX+6], g_spinor_field[DUM_MATRIX+2], VOLUME/2);		// g_spinor_field[DUM_MATRIX+6] = g_spinor_field[DUM_MATRIX+2] ?!= l_charm'
    assign(g_spinor_field[DUM_MATRIX+7], g_spinor_field[DUM_MATRIX+3], VOLUME/2);		// g_spinor_field[DUM_MATRIX+7] = g_spinor_field[DUM_MATRIX+3] ?!= l_strange'
  #else
    to_device(dev_spin_eo2_dn, g_spinor_field[DUM_MATRIX+2], h2d_spin_up, dev_spinsize);
      dev_copy_spinor_field<<<gridsize4, blocksize4>>>(dev_spin_eo2_dn, spinin_up);		// spinin_up = dev_spin_eo2_dn
    to_host(g_spinor_field[DUM_MATRIX+6], spinin_up, h2d_spin_up, dev_spinsize);
    
    to_device(dev_spin_eo2_up, g_spinor_field[DUM_MATRIX+3], h2d_spin_dn, dev_spinsize);
      dev_copy_spinor_field<<<gridsize4, blocksize4>>>(dev_spin_eo2_up, spinin_dn);		// spinin_dn = dev_spin_eo2_up
    to_host(g_spinor_field[DUM_MATRIX+7], spinin_dn, h2d_spin_dn, dev_spinsize);
  #endif




  /* AND THEN THE  Qhat(2x2)  PART */							// notice the swapping !


  /* Here the  M_oe Mee^-1 M_eo  implementation  */		// SWAP:
  
  #ifndef CHECK_HOPPING_MATRIX
    Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX]  , g_spinor_field[DUM_MATRIX+7]);	// g_spinor_field[DUM_MATRIX]   = (M_eo) * g_spinor_field[DUM_MATRIX+7] = (M_eo) * l_strange'		// notice the order
    Hopping_Matrix(EO, g_spinor_field[DUM_MATRIX+1], g_spinor_field[DUM_MATRIX+6]);	// g_spinor_field[DUM_MATRIX+1] = (M_eo) * g_spinor_field[DUM_MATRIX+6] = (M_eo) * l_charm'		//   of l_strange and l_charm
  #else
  				// xchange
  				#ifndef HOPPING_DEBUG
  				  xchange_field_wrapper(dev_spin_eo3_up, 0);
  				#endif
  		  #ifndef HOPPING_DEBUG
  		    dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo3_up, dev_spin_eo1_up, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0);	// dev_spin_eo1_up = (M_eo) * dev_spin_eo3_up
  		  #else
  		    Hopping_Matrix_wrapper(0, dev_spin_eo1_up, dev_spin_eo3_up);
  		  #endif
  				// xchange
  				#ifndef HOPPING_DEBUG
  				  xchange_field_wrapper(dev_spin_eo3_dn, 0);
  				#endif
  		  #ifndef HOPPING_DEBUG
  		    dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo3_dn, dev_spin_eo1_dn, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0);	// dev_spin_eo1_dn = (M_eo) * dev_spin_eo3_dn
  		  #else
  		    Hopping_Matrix_wrapper(0, dev_spin_eo1_dn, dev_spin_eo3_dn);
  		  #endif
  #endif
  
  
  

  #ifndef CHECK_IMUGAMMA5
    // remark: here we don't need  g_mu = -g_mu						//				remark: here the factor  GAMMA5  is not written:
    mul_one_minus_imubar(g_spinor_field[DUM_MATRIX+2], g_spinor_field[DUM_MATRIX]);	// g_spinor_field[DUM_MATRIX+2] = (1 - imubar) * g_spinor_field[DUM_MATRIX]   = (1 - imubar)*(M_eo) * l_strange
    mul_one_plus_imubar (g_spinor_field[DUM_MATRIX+3], g_spinor_field[DUM_MATRIX+1]);	// g_spinor_field[DUM_MATRIX+3] = (1 + imubar) * g_spinor_field[DUM_MATRIX+1] = (1 + imubar)*(M_eo) * l_charm
  #else
    to_device(dev_spin_eo1_up, g_spinor_field[DUM_MATRIX], h2d_spin_up, dev_spinsize);
      dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo1_up, dev_spin_eo2_up, -1.0);			// dev_spin_eo2_up  =  (1 - imubar) * dev_spin_eo1_up
    to_host(g_spinor_field[DUM_MATRIX+2], dev_spin_eo2_up, h2d_spin_up, dev_spinsize);
    
    to_device(dev_spin_eo1_dn, g_spinor_field[DUM_MATRIX+1], h2d_spin_dn, dev_spinsize);
      dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo1_dn, dev_spin_eo2_dn, +1.0);			// dev_spin_eo2_dn  =  (1 + imubar) * dev_spin_eo1_dn
    to_host(g_spinor_field[DUM_MATRIX+3], dev_spin_eo2_dn, h2d_spin_dn, dev_spinsize);
  #endif
  
  

  #ifndef CHECK_CUBLAS1
    assign_add_mul_r(g_spinor_field[DUM_MATRIX+2], g_spinor_field[DUM_MATRIX+1], g_epsbar, VOLUME/2);
    											// g_spinor_field[DUM_MATRIX+2] = g_spinor_field[DUM_MATRIX+2]   +  epsbar * g_spinor_field[DUM_MATRIX+1]
    											//                              = (1 - imubar)*(M_eo)*l_strange  +  epsbar * (M_eo) * l_charm
    assign_add_mul_r(g_spinor_field[DUM_MATRIX+3], g_spinor_field[DUM_MATRIX]  , g_epsbar, VOLUME/2);
    											// g_spinor_field[DUM_MATRIX+3] = g_spinor_field[DUM_MATRIX+3]   +  epsbar * g_spinor_field[DUM_MATRIX]
    											//                              = (1 + imubar)*(M_eo)*l_charm    +  epsbar * (M_eo) * l_strange
    
    mul_r(g_spinor_field[DUM_MATRIX+2], nrm, g_spinor_field[DUM_MATRIX+2], VOLUME/2);	// g_spinor_field[DUM_MATRIX+2] = nrm * g_spinor_field[DUM_MATRIX+2] = nrm*(1 - imubar)*(M_eo)*l_strange + nrm*epsbar*(M_eo)*l_charm
    mul_r(g_spinor_field[DUM_MATRIX+3], nrm, g_spinor_field[DUM_MATRIX+3], VOLUME/2);	// g_spinor_field[DUM_MATRIX+3] = nrm * g_spinor_field[DUM_MATRIX+3] = nrm*(1 + imubar)*(M_eo)*l_charm   + nrm*epsbar*(M_eo)*l_strange
  #else
    to_device(dev_spin_eo1_dn, g_spinor_field[DUM_MATRIX+1], h2d_spin_up, dev_spinsize);
    to_device(dev_spin_eo2_up, g_spinor_field[DUM_MATRIX+2], h2d_spin_up, dev_spinsize);
      cublasSaxpy (N_floats, g_epsbar, (float *) dev_spin_eo1_dn, 1, (float *) dev_spin_eo2_up, 1);
    										// dev_spin_eo2_up  =  dev_spin_eo2_up  +  epsbar * dev_spin_eo1_dn
      cublasSscal (N_floats, nrm, (float *) dev_spin_eo2_up, 1);		// dev_spin_eo2_up  =  nrm * dev_spin_eo2_up
    to_host(g_spinor_field[DUM_MATRIX+2], dev_spin_eo2_up, h2d_spin_up, dev_spinsize);
    
    
    to_device(dev_spin_eo1_up, g_spinor_field[DUM_MATRIX], h2d_spin_dn, dev_spinsize);
    to_device(dev_spin_eo2_dn, g_spinor_field[DUM_MATRIX+3], h2d_spin_dn, dev_spinsize);
      cublasSaxpy (N_floats, g_epsbar, (float *) dev_spin_eo1_up, 1, (float *) dev_spin_eo2_dn, 1);
    										// dev_spin_eo2_dn  =  dev_spin_eo2_dn  +  epsbar * dev_spin_eo1_up
      cublasSscal (N_floats, nrm, (float *) dev_spin_eo2_dn, 1);		// dev_spin_eo2_dn  =  nrm * dev_spin_eo2_dn
    to_host(g_spinor_field[DUM_MATRIX+3], dev_spin_eo2_dn, h2d_spin_dn, dev_spinsize);
  #endif
  
  
  
   
  #ifndef CHECK_HOPPING_MATRIX
    Hopping_Matrix(OE, l_strange, g_spinor_field[DUM_MATRIX+2]);				// l_strange = (M_oe) * g_spinor_field[DUM_MATRIX+2] = (M_oe)*nrm*(1 - imubar)*(M_eo)*l_strange + (M_oe)*nrm*epsbar*(M_eo)*l_charm
    Hopping_Matrix(OE, l_charm  , g_spinor_field[DUM_MATRIX+3]);				// l_charm   = (M_oe) * g_spinor_field[DUM_MATRIX+3] = (M_oe)*nrm*(1 + imubar)*(M_eo)*l_charm   + (M_oe)*nrm*epsbar*(M_eo)*l_strange
  #else
  				// xchange
  				#ifndef HOPPING_DEBUG
  				  xchange_field_wrapper(dev_spin_eo2_up, 1);
  				#endif
  		  #ifndef HOPPING_DEBUG
  		    dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo2_up, dev_spin_eo1_up, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1);
  		  #else
  		    Hopping_Matrix_wrapper(1, dev_spin_eo1_up, dev_spin_eo2_up);
  		  #endif
  				// xchange
  				#ifndef HOPPING_DEBUG
  				  xchange_field_wrapper(dev_spin_eo2_dn, 1);
  				#endif
  		  #ifndef HOPPING_DEBUG
  		    dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo2_dn, dev_spin_eo1_dn, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1);
  		  #else
  		    Hopping_Matrix_wrapper(1, dev_spin_eo1_dn, dev_spin_eo2_dn);
  		  #endif
  #endif




  /* Here the M_oo  implementation  */
  
  #ifndef CHECK_IMUGAMMA5  
    mul_one_plus_imubar (g_spinor_field[DUM_MATRIX]  , g_spinor_field[DUM_MATRIX+7]);	// g_spinor_field[DUM_MATRIX]   = (1 + imubar) * g_spinor_field[DUM_MATRIX+7] = (1 + imubar) * l_strange
    mul_one_minus_imubar(g_spinor_field[DUM_MATRIX+1], g_spinor_field[DUM_MATRIX+6]);	// g_spinor_field[DUM_MATRIX+1] = (1 - imubar) * g_spinor_field[DUM_MATRIX+6] = (1 - imubar) * l_charm
  #else
    to_device(spinin_up, g_spinor_field[DUM_MATRIX+7], h2d_spin_up, dev_spinsize);
      dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(spinin_up, dev_spin_eo2_up, +1.0);	// dev_spin_eo2_up = (1 + imubar) * spinin_up
    to_host(g_spinor_field[DUM_MATRIX], dev_spin_eo2_up, h2d_spin_up, dev_spinsize);
    
    to_device(spinin_dn, g_spinor_field[DUM_MATRIX+6], h2d_spin_dn, dev_spinsize);
      dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(spinin_dn, dev_spin_eo2_dn, -1.0);	// dev_spin_eo2_dn = (1 - imubar) * spinin_dn
    to_host(g_spinor_field[DUM_MATRIX+1], dev_spin_eo2_dn, h2d_spin_dn, dev_spinsize);
  #endif
  
  
  

  #ifndef CHECK_CUBLAS2
    assign_add_mul_r(g_spinor_field[DUM_MATRIX]  , g_spinor_field[DUM_MATRIX+6], -g_epsbar, VOLUME/2);
    											// g_spinor_field[DUM_MATRIX]   = g_spinor_field[DUM_MATRIX]    -  epsbar * g_spinor_field[DUM_MATRIX+6]
    											//                              = (1 + imubar) * l_strange      -  epsbar * l_charm
    assign_add_mul_r(g_spinor_field[DUM_MATRIX+1], g_spinor_field[DUM_MATRIX+7], -g_epsbar, VOLUME/2);
    											// g_spinor_field[DUM_MATRIX+1] = g_spinor_field[DUM_MATRIX+1]  -  epsbar * g_spinor_field[DUM_MATRIX+7]
    											//                              = (1 - imubar) * l_charm        -  epsbar * l_strange
  #else
    to_device(spinin_dn, g_spinor_field[DUM_MATRIX+6], h2d_spin_up, dev_spinsize);
    to_device(dev_spin_eo2_up, g_spinor_field[DUM_MATRIX], h2d_spin_up, dev_spinsize);
      cublasSaxpy (N_floats, -g_epsbar, (float *) spinin_dn, 1, (float *) dev_spin_eo2_up, 1);	// dev_spin_eo2_up  =  dev_spin_eo2_up - epsbar*spinin_dn
    to_host(g_spinor_field[DUM_MATRIX], dev_spin_eo2_up, h2d_spin_up, dev_spinsize);
    
    to_device(spinin_up, g_spinor_field[DUM_MATRIX+7], h2d_spin_dn, dev_spinsize);
    to_device(dev_spin_eo2_dn, g_spinor_field[DUM_MATRIX+1], h2d_spin_dn, dev_spinsize);
      cublasSaxpy (N_floats, -g_epsbar, (float *) spinin_up, 1, (float *) dev_spin_eo2_dn, 1);	// dev_spin_eo2_dn  =  dev_spin_eo2_dn - epsbar*spinin_up
    to_host(g_spinor_field[DUM_MATRIX+1], dev_spin_eo2_dn, h2d_spin_dn, dev_spinsize);
  #endif
  
  
  
  
  // here the (M_oo - M_oe Mee^-1 M_eo) implementation
  
  #ifndef CHECK_CUBLAS3
    diff(l_strange, g_spinor_field[DUM_MATRIX]  , l_strange, VOLUME/2);			// l_strange = g_spinor_field[DUM_MATRIX]   - l_strange
    											//           =                         (1 + imubar) * l_strange -                   epsbar * l_charm
    											//                -  (M_oe)*nrm*(1 - imubar)*(M_eo) * l_strange + (M_oe)*nrm*epsbar*(M_eo) * l_charm
    											
    diff(l_charm  , g_spinor_field[DUM_MATRIX+1], l_charm  , VOLUME/2);			// l_charm   = g_spinor_field[DUM_MATRIX+1] - l_charm
    											//           =                         (1 - imubar) * l_charm   -                   epsbar * l_strange
    											//                -  (M_oe)*nrm*(1 + imubar)*(M_eo) * l_charm   + (M_oe)*nrm*epsbar*(M_eo) * l_strange
  #else
    to_device(dev_spin_eo1_up, l_strange, h2d_spin_up, dev_spinsize);
    to_device(dev_spin_eo2_up, g_spinor_field[DUM_MATRIX], h2d_spin_up, dev_spinsize);
      cublasSaxpy (N_floats, -1.0, (float *) dev_spin_eo1_up, 1, (float *) dev_spin_eo2_up, 1);	// dev_spin_eo2_up  =  dev_spin_eo2_up  -  dev_spin_eo1_up
    to_host(l_strange, dev_spin_eo2_up, h2d_spin_up, dev_spinsize);
    
    to_device(dev_spin_eo1_dn, l_charm, h2d_spin_dn, dev_spinsize);
    to_device(dev_spin_eo2_dn, g_spinor_field[DUM_MATRIX+1], h2d_spin_dn, dev_spinsize);
      cublasSaxpy (N_floats, -1.0, (float *) dev_spin_eo1_dn, 1, (float *) dev_spin_eo2_dn, 1);	// dev_spin_eo2_dn  =  dev_spin_eo2_dn  -  dev_spin_eo1_dn
    to_host(l_charm, dev_spin_eo2_dn, h2d_spin_dn, dev_spinsize);
  #endif
  
  
  
  
  /* and finally the  GAMMA5  multiplication  */
  
  #ifndef CHECK_GAMMA5
    gamma5(l_strange, l_strange, VOLUME/2);						// l_strange = gamma5 * l_strange
    gamma5(l_charm  , l_charm  , VOLUME/2);						// l_charm   = gamma5 * l_charm
  #else
    to_device(dev_spin_eo2_up, l_strange, h2d_spin_up, dev_spinsize);
      dev_gamma5<<<gridsize3, blocksize3>>>(dev_spin_eo2_up, dev_spin_eo2_up);		// dev_spin_eo2_up = gamma5 * dev_spin_eo2_up
    to_host(l_strange, dev_spin_eo2_up, h2d_spin_up, dev_spinsize);
    
    to_device(dev_spin_eo2_dn, l_charm, h2d_spin_dn, dev_spinsize);
      dev_gamma5<<<gridsize3, blocksize3>>>(dev_spin_eo2_dn, dev_spin_eo2_dn);		// dev_spin_eo2_dn = gamma5 * dev_spin_eo2_dn
    to_host(l_charm, dev_spin_eo2_dn, h2d_spin_dn, dev_spinsize);
  #endif
  
  


  /* At the end, the normalisation by the max. eigenvalue  */ 
  /* Twice  phmc_invmaxev  since we consider here  D Ddag  !!! */
  mul_r(l_charm, phmc_invmaxev*phmc_invmaxev, l_charm, VOLUME/2);
  mul_r(l_strange, phmc_invmaxev*phmc_invmaxev, l_strange, VOLUME/2);
  return;
  
}//matrix__mpi_debug10()











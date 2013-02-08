

//the additional momentum fields we
dev_su3adj* dev_df1_d;
dev_su3adj* dev_df2_d;
dev_su3adj* dev_mom_update;


extern "C" {
#include "../fieldstrength.h"
}



//this allocates two more momentum fields on device which we need for the Wilson Flow
extern "C" void wf_init_additional_gpu_fields(){
  cudaError_t cudaerr;

  size_t dev_momsize = 4*VOLUME * sizeof(dev_su3adj);
  if((cudaerr=cudaMalloc((void **) &dev_df1_d, dev_momsize)) != cudaSuccess){
    printf("Error in init_gpu_fields(): Memory allocation of double momentum field failed. Aborting...\n");
    exit(200);
  }   // Allocate array on device
  else{
    printf("Allocated double momentum field dev_df1_d on device\n");
  } 
  

  if((cudaerr=cudaMalloc((void **) &dev_df2_d, dev_momsize)) != cudaSuccess){
    printf("Error in init_gpu_fields(): Memory allocation of double momentum field failed. Aborting...\n");
    exit(200);
  }   // Allocate array on device
  else{
    printf("Allocated double momentum field dev_df2_d on device\n");
  } 


  if((cudaerr=cudaMalloc((void **) &dev_mom_update, dev_momsize)) != cudaSuccess){
    printf("Error in init_gpu_fields(): Memory allocation of double momentum field failed. Aborting...\n");
    exit(200);
  }   // Allocate array on device
  else{
    printf("Allocated double momentum field dev_mom_update on device\n");
  } 


  if ((cudaerr=cudaGetLastError())!=cudaSuccess) {
   printf("%s\n", cudaGetErrorString(cudaGetLastError()));
  }  
}

extern "C" void wf_finalize_additional_gpu_fields(){
  cudaFree(dev_df1_d);
  cudaFree(dev_df2_d);
}



// this is supposed to calculate the derivative of the gauge field and
// to store it in the momentum field mf
// withrectangles==1 : rectangles are also calculated
// withrectangles==0 : rectangles are NOT  calculated
void gpu_deriv_gauge(int withrectangles, double c, dev_su3adj* mf){
cudaError_t cudaerr;
//printf("gpu_dev_gauge:\n");
//plaquette
cudaFuncSetCacheConfig(dev_gauge_derivative, cudaFuncCachePreferL1);
dev_gauge_derivative <<<griddimgauge, blockdimgauge >>> 
            (dev_nn2, dev_gf_d, mf);          

//rectangles            
if(withrectangles){         
  int mu;
  //we had to split this up into three kernels to keep a single
  //kernel small (nvcc seems to have some problems with large
  //kernels)
  
  
  for(mu=0; mu<4; mu++)
  {
  //here we do
  //    __ __                    __ __
  //  ||__ __|         and      |__ __||
  //                             
    dev_rect_deri_one <<<griddimgauge, blockdimgauge >>> 
            (dev_nn2, dev_gf_d, mf, c, mu);  
  //here we do
  //    __                         __
  //   |  |                      ||  | 
  //  ||__|            and        |__|
  //                          
    dev_rect_deri_two <<<griddimgauge, blockdimgauge >>> 
              (dev_nn2, dev_gf_d, mf, c, mu);  
  //here we do
  //   __                         __    
  //  |  |            and        |  ||    
  //  |__||                      |__|   
  //     
    dev_rect_deri_three <<<griddimgauge, blockdimgauge >>> 
              (dev_nn2, dev_gf_d, mf, c, mu);  
  }
}                

}




void gpu_upd_gauge(double step, dev_su3adj* mf){
  dev_gauge_update <<<griddimgauge, blockdimgauge >>> 
            (dev_gf_d, mf, step);   
}




void zero_moment_field(su3adj** mf){
int ix,mu;
  /* set ddummy to zero */
  for(ix = 0; ix < VOLUME+RAND; ix++){
    for(mu=0; mu<4; mu++){
      mf[ix][mu].d1=0.;
      mf[ix][mu].d2=0.;
      mf[ix][mu].d3=0.;
      mf[ix][mu].d4=0.;
      mf[ix][mu].d5=0.;
      mf[ix][mu].d6=0.;
      mf[ix][mu].d7=0.;
      mf[ix][mu].d8=0.;
    }
  }

}


extern "C" void wf_init(){
  zero_moment_field(df0);
  //X0 = U_mu
  update_gpu_fields(g_gauge_field, df0,1);
  printf("Have updated gpu fields\n");
}

extern "C" void wf_finalize(){
  to_host_gf(g_gauge_field, h2d_gf_d);
}




void show_su3_2(su3 gf1){
   printf("(%f,%f)\t(%f,%f)\t(%f,%f)\n",gf1.c00.re,
   					gf1.c00.im,
   					gf1.c01.re,
   					gf1.c01.im,
   					gf1.c02.re,
   					gf1.c02.im
   );
   printf("(%f,%f)\t(%f,%f)\t(%f,%f)\n",gf1.c10.re,
   					gf1.c10.im,
   					gf1.c11.re,
   					gf1.c11.im,
   					gf1.c12.re,
   					gf1.c12.im
   );
   printf("(%f,%f)\t(%f,%f)\t(%f,%f)\n",gf1.c20.re,
   					gf1.c20.im,
   					gf1.c21.re,
   					gf1.c21.im,
   					gf1.c22.re,
   					gf1.c22.im
   ); 
}


void check_unitarity(su3 gf1){
   double row1squared = gf1.c00.re*gf1.c00.re + gf1.c00.im*gf1.c00.im + 
                        gf1.c01.re*gf1.c01.re + gf1.c01.im*gf1.c01.im +
                        gf1.c02.re*gf1.c02.re + gf1.c02.im*gf1.c02.im;

   double row2squared = gf1.c10.re*gf1.c10.re + gf1.c10.im*gf1.c10.im + 
                        gf1.c11.re*gf1.c11.re + gf1.c11.im*gf1.c11.im +
                        gf1.c12.re*gf1.c12.re + gf1.c12.im*gf1.c12.im;

   double row3squared = gf1.c20.re*gf1.c20.re + gf1.c20.im*gf1.c20.im + 
                        gf1.c21.re*gf1.c21.re + gf1.c21.im*gf1.c21.im +
                        gf1.c22.re*gf1.c22.re + gf1.c22.im*gf1.c22.im;

   printf("1st row squared: %e\t2nd row squared: %e\t3rd row squared: %e\n",
                                          row1squared,row2squared,row3squared);

}

extern "C" void wf_rk_step(int withrectangles,double eps){
cudaError_t cudaerr;
int i;
       

//set everything to 0 because in gpu_deriv_gauge the derivative is added to the given momentum field
  dev_zero_mom_field<<<griddimgauge, blockdimgauge >>>(dev_mom_update);  
  dev_zero_mom_field<<<griddimgauge, blockdimgauge >>>(dev_df0_d);
  dev_zero_mom_field<<<griddimgauge, blockdimgauge >>>(dev_df1_d);
  dev_zero_mom_field<<<griddimgauge, blockdimgauge >>>(dev_df2_d);


  //printf("1st update of flow with stepsize %e\n",eps);
  //X1 = exp(1/4 Z0) X0 
  gpu_deriv_gauge(withrectangles, eps, dev_df0_d);
  dev_assign_const_times_mom_field<<<griddimgauge, blockdimgauge >>>(dev_mom_update, 1./4, dev_df0_d);
  gpu_upd_gauge(1.0, dev_mom_update);

  //printf("2nd update of flow with stepsize %e\n",eps);
  //X2 = exp(8/9 Z1 - 17/36 Z0) X1 
  gpu_deriv_gauge(withrectangles, eps, dev_df1_d);
  dev_assign_const_times_mom_field<<<griddimgauge, blockdimgauge >>>(dev_mom_update, -17./36., dev_df0_d);  
  dev_add_const_times_mom_field<<<griddimgauge, blockdimgauge >>>(dev_mom_update, 8./9., dev_df1_d);
  gpu_upd_gauge(1.0, dev_mom_update); 


//   dev_su3adj help[6];
//   size_t helpsize = 6* sizeof(dev_su3adj);
//   cudaMemcpy(&help, dev_mom_update, helpsize, cudaMemcpyDeviceToHost);
//   for(i=0; i<6; i++){
//       printf("%e, ",help[i].d1);
//       printf("%e, ",help[i].d2);
//       printf("%e, ",help[i].d3);
//       printf("%e, ",help[i].d4);
//       printf("%e, ",help[i].d5);
//       printf("%e, ",help[i].d6);
//       printf("%e, ",help[i].d7);
//       printf("%e, ",help[i].d8);
//       printf("\n");
//   }
//  
//   to_host_gf(g_gauge_field, h2d_gf_d);
//   show_su3_2(g_gauge_field[0][0]);
//   check_unitarity(g_gauge_field[0][0]);


  //printf("3rd update of flow with stepsize %e\n",eps);
  //U_mu = exp(3/4 Z2 - 8/9 Z1 + 17/36 Z0) X2 
  gpu_deriv_gauge(withrectangles, eps, dev_df2_d);
  dev_assign_const_times_mom_field<<<griddimgauge, blockdimgauge >>>(dev_mom_update, +17./36., dev_df0_d);  
  dev_add_const_times_mom_field<<<griddimgauge, blockdimgauge >>>(dev_mom_update, -8./9., dev_df1_d);
  dev_add_const_times_mom_field<<<griddimgauge, blockdimgauge >>>(dev_mom_update, +3./4., dev_df2_d);
  gpu_upd_gauge(1.0, dev_mom_update);  


  if ((cudaerr=cudaPeekAtLastError())!=cudaSuccess) {
    printf("%s\n", cudaGetErrorString(cudaPeekAtLastError()));
  }           
       
  to_host_gf(g_gauge_field, h2d_gf_d);
  check_unitarity(g_gauge_field[0][0]);
}











//here is adopted the gauge_update as performed in hybrid_update.c in double precision







#define _make_su3_dev(v,p) \
(v)[0][0].im= (p).d3+0.5773502691896258*(p).d8; \
(v)[0][1].im= (p).d1; \
(v)[0][2].im= (p).d4; \
(v)[1][0].im= (p).d1; \
(v)[1][1].im=-(p).d3+0.5773502691896258*(p).d8; \
(v)[1][2].im= (p).d6; \
(v)[2][0].im= (p).d4; \
(v)[2][1].im= (p).d6; \
(v)[2][2].im=-1.154700538379252*(p).d8; \
(v)[0][0].re= 0.0; \
(v)[0][1].re= (p).d2; \
(v)[0][2].re= (p).d5; \
(v)[1][0].re=-(p).d2; \
(v)[1][1].re= 0.0; \
(v)[1][2].re= (p).d7; \
(v)[2][0].re=-(p).d5; \
(v)[2][1].re=-(p).d7; \
(v)[2][2].re= 0.0;


#define _assign_const_times_mom_dev(res,c,in) \
(res).d1=(c)*(in).d1; \
(res).d2=(c)*(in).d2; \
(res).d3=(c)*(in).d3; \
(res).d4=(c)*(in).d4; \
(res).d5=(c)*(in).d5; \
(res).d6=(c)*(in).d6; \
(res).d7=(c)*(in).d7; \
(res).d8=(c)*(in).d8; 


__device__  void dev_zero_mom(dev_su3adj *res){
  (*res).d1=0.; 
  (*res).d2=0.; 
  (*res).d3=0.; 
  (*res).d4=0.;
  (*res).d5=0.; 
  (*res).d6=0.; 
  (*res).d7=0.; 
  (*res).d8=0.; 
}


__device__  void dev_assign_const_times_mom(dev_su3adj *res, double c, dev_su3adj *in){
  (*res).d1=(c)*(*in).d1; 
  (*res).d2=(c)*(*in).d2; 
  (*res).d3=(c)*(*in).d3; 
  (*res).d4=(c)*(*in).d4; 
  (*res).d5=(c)*(*in).d5; 
  (*res).d6=(c)*(*in).d6; 
  (*res).d7=(c)*(*in).d7; 
  (*res).d8=(c)*(*in).d8; 
}


__device__  void dev_add_const_times_mom(dev_su3adj *res, double c, dev_su3adj *in){
  (*res).d1+=(c)*(*in).d1; 
  (*res).d2+=(c)*(*in).d2; 
  (*res).d3+=(c)*(*in).d3; 
  (*res).d4+=(c)*(*in).d4; 
  (*res).d5+=(c)*(*in).d5; 
  (*res).d6+=(c)*(*in).d6; 
  (*res).d7+=(c)*(*in).d7; 
  (*res).d8+=(c)*(*in).d8; 
}


__device__ void dev_storegf_2v_d(int pos, dev_su3_2v_d* gfield , dev_su3_d* U){

   gfield[6*pos].x = (*U)[0][0].re;
   gfield[6*pos].y = (*U)[0][0].im;
   gfield[6*pos+1].x = (*U)[0][1].re;
   gfield[6*pos+1].y = (*U)[0][1].im;
   
   gfield[6*pos+2].x = (*U)[0][2].re;
   gfield[6*pos+2].y = (*U)[0][2].im;
   gfield[6*pos+3].x = (*U)[1][0].re;
   gfield[6*pos+3].y = (*U)[1][0].im;
   
   gfield[6*pos+4].x = (*U)[1][1].re;
   gfield[6*pos+4].y = (*U)[1][1].im;
   gfield[6*pos+5].x = (*U)[1][2].re;
   gfield[6*pos+5].y = (*U)[1][2].im;
   
}







__device__ void dev_exposu3(dev_su3_d * vr, dev_su3adj *p) {

  int i;
  dev_su3_d v,v2;
  double fac,r;
  double a,b;
  dev_complex_d a0,a1,a2,a1p;

  /* it writes 'p=vec(h_{j,mu})' in matrix form 'v' */  
  _make_su3_dev(v,*p);
  /* calculates v^2 */
  dev_su3_ti_su3_d(&(v2), &(v), &(v));
  /* */
  a=0.5*(v2[0][0].re+v2[1][1].re+v2[2][2].re);
  /* 1/3 imaginary part of tr v*v2 */
  b = 0.33333333333333333*
    (v[0][0].re*v2[0][0].im+v[0][0].im*v2[0][0].re
     +v[0][1].re*v2[1][0].im+v[0][1].im*v2[1][0].re
     +v[0][2].re*v2[2][0].im+v[0][2].im*v2[2][0].re
     +v[1][0].re*v2[0][1].im+v[1][0].im*v2[0][1].re
     +v[1][1].re*v2[1][1].im+v[1][1].im*v2[1][1].re
     +v[1][2].re*v2[2][1].im+v[1][2].im*v2[2][1].re
     +v[2][0].re*v2[0][2].im+v[2][0].im*v2[0][2].re
     +v[2][1].re*v2[1][2].im+v[2][1].im*v2[1][2].re
     +v[2][2].re*v2[2][2].im+v[2][2].im*v2[2][2].re  );
  a0.re=0.16059043836821615e-9;    /*  1/13! */
  a0.im=0.0;
  a1.re=0.11470745597729725e-10;   /*  1/14! */
  a1.im=0.0;
  a2.re=0.76471637318198165e-12;   /*  1/15! */
  a2.im=0.0;
  fac=0.20876756987868099e-8;      /*  1/12! */
  r=12.0;
  for(i = 3; i <= 15; i++) {
    a1p.re = a0.re + a * a2.re;
    a1p.im = a0.im + a * a2.im;
    a0.re = fac - b * a2.im;
    a0.im =     + b * a2.re;
    a2.re = a1.re; 
    a2.im = a1.im;
    a1.re = a1p.re; 
    a1.im = a1p.im;
    fac *= r;  
    r -= 1.0;
  }
  /* vr = a0 + a1*v + a2*v2 */
  (*vr)[0][0].re = a0.re + a1.re*v[0][0].re - a1.im*v[0][0].im + a2.re*v2[0][0].re - a2.im*v2[0][0].im;
  (*vr)[0][0].im = a0.im + a1.re*v[0][0].im + a1.im*v[0][0].re + a2.re*v2[0][0].im + a2.im*v2[0][0].re;
  (*vr)[0][1].re =         a1.re*v[0][1].re - a1.im*v[0][1].im + a2.re*v2[0][1].re - a2.im*v2[0][1].im;
  (*vr)[0][1].im =         a1.re*v[0][1].im + a1.im*v[0][1].re + a2.re*v2[0][1].im + a2.im*v2[0][1].re;
  (*vr)[0][2].re =         a1.re*v[0][2].re - a1.im*v[0][2].im + a2.re*v2[0][2].re - a2.im*v2[0][2].im;
  (*vr)[0][2].im =         a1.re*v[0][2].im + a1.im*v[0][2].re + a2.re*v2[0][2].im + a2.im*v2[0][2].re;
  (*vr)[1][0].re =         a1.re*v[1][0].re - a1.im*v[1][0].im + a2.re*v2[1][0].re - a2.im*v2[1][0].im;
  (*vr)[1][0].im =         a1.re*v[1][0].im + a1.im*v[1][0].re + a2.re*v2[1][0].im + a2.im*v2[1][0].re;
  (*vr)[1][1].re = a0.re + a1.re*v[1][1].re - a1.im*v[1][1].im + a2.re*v2[1][1].re - a2.im*v2[1][1].im;
  (*vr)[1][1].im = a0.im + a1.re*v[1][1].im + a1.im*v[1][1].re + a2.re*v2[1][1].im + a2.im*v2[1][1].re;
  (*vr)[1][2].re =         a1.re*v[1][2].re - a1.im*v[1][2].im + a2.re*v2[1][2].re - a2.im*v2[1][2].im;
  (*vr)[1][2].im =         a1.re*v[1][2].im + a1.im*v[1][2].re + a2.re*v2[1][2].im + a2.im*v2[1][2].re;
  (*vr)[2][0].re =         a1.re*v[2][0].re - a1.im*v[2][0].im + a2.re*v2[2][0].re - a2.im*v2[2][0].im;
  (*vr)[2][0].im =         a1.re*v[2][0].im + a1.im*v[2][0].re + a2.re*v2[2][0].im + a2.im*v2[2][0].re;
  (*vr)[2][1].re =         a1.re*v[2][1].re - a1.im*v[2][1].im + a2.re*v2[2][1].re - a2.im*v2[2][1].im;
  (*vr)[2][1].im =         a1.re*v[2][1].im + a1.im*v[2][1].re + a2.re*v2[2][1].im + a2.im*v2[2][1].re;
  (*vr)[2][2].re = a0.re + a1.re*v[2][2].re - a1.im*v[2][2].im + a2.re*v2[2][2].re - a2.im*v2[2][2].im;
  (*vr)[2][2].im = a0.im + a1.re*v[2][2].im + a1.im*v[2][2].re + a2.re*v2[2][2].im + a2.im*v2[2][2].re;
}







__device__ void dev_restoresu3(dev_su3_d* vr, const dev_su3_d * u) {
  double n1,n2;
  
  /* normalize rows 1 and 2 */
  n1= (*u)[0][0].re * (*u)[0][0].re + (*u)[0][0].im * (*u)[0][0].im
    + (*u)[0][1].re * (*u)[0][1].re + (*u)[0][1].im * (*u)[0][1].im
    + (*u)[0][2].re * (*u)[0][2].re + (*u)[0][2].im * (*u)[0][2].im;
  n1 = 1.0/sqrt(n1);
  n2= (*u)[1][0].re * (*u)[1][0].re + (*u)[1][0].im * (*u)[1][0].im
    + (*u)[1][1].re * (*u)[1][1].re + (*u)[1][1].im * (*u)[1][1].im
    + (*u)[1][2].re * (*u)[1][2].re + (*u)[1][2].im * (*u)[1][2].im;
  n2= 1.0/sqrt(n2);
  
  (*vr)[0][0].re=n1*(*u)[0][0].re;  (*vr)[0][0].im=n1*(*u)[0][0].im;
  (*vr)[0][1].re=n1*(*u)[0][1].re;  (*vr)[0][1].im=n1*(*u)[0][1].im;
  (*vr)[0][2].re=n1*(*u)[0][2].re;  (*vr)[0][2].im=n1*(*u)[0][2].im;
  
  (*vr)[1][0].re=n2*(*u)[1][0].re;  (*vr)[1][0].im=n2*(*u)[1][0].im;
  (*vr)[1][1].re=n2*(*u)[1][1].re;  (*vr)[1][1].im=n2*(*u)[1][1].im;
  (*vr)[1][2].re=n2*(*u)[1][2].re;  (*vr)[1][2].im=n2*(*u)[1][2].im;
  

  // 1 = 2 3  - 3 2
  (*vr)[2][0].re= (*vr)[0][1].re*(*vr)[1][2].re-(*vr)[0][2].re*(*vr)[1][1].re
    -(*vr)[0][1].im*(*vr)[1][2].im+(*vr)[0][2].im*(*vr)[1][1].im;
  (*vr)[2][0].im=-(*vr)[0][1].re*(*vr)[1][2].im+(*vr)[0][2].re*(*vr)[1][1].im
    -(*vr)[0][1].im*(*vr)[1][2].re+(*vr)[0][2].im*(*vr)[1][1].re;
  // 2 = 3 1  - 1 3
  (*vr)[2][1].re= (*vr)[0][2].re*(*vr)[1][0].re-(*vr)[0][0].re*(*vr)[1][2].re
    -(*vr)[0][2].im*(*vr)[1][0].im+(*vr)[0][0].im*(*vr)[1][2].im;
  (*vr)[2][1].im=-(*vr)[0][2].re*(*vr)[1][0].im+(*vr)[0][0].re*(*vr)[1][2].im
    -(*vr)[0][2].im*(*vr)[1][0].re+(*vr)[0][0].im*(*vr)[1][2].re;
  // 3 = 1 2  - 2 1
  (*vr)[2][2].re= (*vr)[0][0].re*(*vr)[1][1].re-(*vr)[0][1].re*(*vr)[1][0].re
    -(*vr)[0][0].im*(*vr)[1][1].im+(*vr)[0][1].im*(*vr)[1][0].im;
  (*vr)[2][2].im=-(*vr)[0][0].re*(*vr)[1][1].im+(*vr)[0][1].re*(*vr)[1][0].im
    -(*vr)[0][0].im*(*vr)[1][1].re+(*vr)[0][1].im*(*vr)[1][0].re;

  
}








__global__ void dev_gauge_update(dev_su3_2v_d * gf, dev_su3adj* m, double step){

  int pos, mu;
  dev_su3adj deriv;
  dev_su3_d expo;
  dev_su3_d v,w,x;
  pos= threadIdx.x + blockDim.x*blockIdx.x;  

  if(pos < dev_VOL2){
    mu=0;
    for(mu=0; mu<4; mu++){
     
      dev_assign_const_times_mom(&deriv,step,&(m[4*pos+mu]));
      dev_exposu3(&(expo), &(deriv)); 
      dev_restoresu3(&(x), &(expo)); 
      dev_get_matrix(gf, pos, mu, &(v));
      
      
      //unoptimized but working
      //dev_su3_ti_su3_d(&(w), &(x), &(v));
      //dev_storegf_2v_d((4*pos+mu), gf , &(w));
      
      dev_su3_ti_su3_d(&(w), &(x), &(v));
      dev_restoresu3(&(x), &(w));
      dev_storegf_2v_d((4*pos+mu), gf , &(x));
      
      /*
      //optimized and working
      dev_store_su3_ti_su3_d((4*pos+mu), gf , &(x), &(v));
      */
    }
  }
}



__global__ void dev_zero_mom_field(dev_su3adj* out){

  int pos, mu;
  pos= threadIdx.x + blockDim.x*blockIdx.x;  

  if(pos < dev_VOL2){
    mu=0;
    for(mu=0; mu<4; mu++){
     
      dev_zero_mom(&(out[4*pos+mu]));

    }
  }
}


__global__ void dev_add_const_times_mom_field(dev_su3adj* out, double c, dev_su3adj* in){

  int pos, mu;
  pos= threadIdx.x + blockDim.x*blockIdx.x;  

  if(pos < dev_VOL2){
    mu=0;
    for(mu=0; mu<4; mu++){
     
      dev_add_const_times_mom(&(out[4*pos+mu]),c,&(in[4*pos+mu]));

    }
  }
}



__global__ void dev_assign_const_times_mom_field(dev_su3adj* out, double c, dev_su3adj* in){

  int pos, mu;
  pos= threadIdx.x + blockDim.x*blockIdx.x;  

  if(pos < dev_VOL2){
    mu=0;
    for(mu=0; mu<4; mu++){
     
      dev_assign_const_times_mom(&(out[4*pos+mu]),c,&(in[4*pos+mu]));

    }
  }
}



/* copy 2v double gauge field to host and reconstruct third row on host*/
void to_host_gf(su3** gf, dev_su3_2v_d* h2d){
  int i,j;
  cudaError_t cudaerr;
  
  //MPI: INTERNAL volume correct here  
  size_t dev_gfsize = 6*4*VOLUME * sizeof(dev_su3_2v_d);
  cudaMemcpy(h2d, dev_gf_d, dev_gfsize, cudaMemcpyDeviceToHost);
  if ((cudaerr=cudaGetLastError())!=cudaSuccess) {
     printf("%s\n", cudaGetErrorString(cudaGetLastError()));
  }

  for (i=0;i<VOLUME;i++){
    for(j=0; j<4; j++){
       
     //first row
      gf[i][j].c00.re = h2d[6*(4*i+j)].x;
      gf[i][j].c00.im = h2d[6*(4*i+j)].y;
      gf[i][j].c01.re = h2d[6*(4*i+j)+1].x;
      gf[i][j].c01.im = h2d[6*(4*i+j)+1].y;
      gf[i][j].c02.re = h2d[6*(4*i+j)+2].x;
      gf[i][j].c02.im = h2d[6*(4*i+j)+2].y;      
   //second row
      gf[i][j].c10.re = h2d[6*(4*i+j)+3].x;
      gf[i][j].c10.im = h2d[6*(4*i+j)+3].y;
      gf[i][j].c11.re = h2d[6*(4*i+j)+4].x;
      gf[i][j].c11.im = h2d[6*(4*i+j)+4].y;
      gf[i][j].c12.re = h2d[6*(4*i+j)+5].x;
      gf[i][j].c12.im = h2d[6*(4*i+j)+5].y;      
   //restore the third row -> reconstructgf_2v defined in gauge_reconstruction.cuh
      reconstructgf_2v_host(&gf[i][j]);
      
   }
   
 }
}






extern "C" void gpu_gauge_update(double step){
cudaError_t cudaerr;

update_gpu_fields(g_gauge_field, moment,1);
dev_gauge_update <<<griddimgauge, blockdimgauge >>> 
            (dev_gf_d, dev_df0_d, step);   
if ((cudaerr=cudaGetLastError())!=cudaSuccess) {
   printf("%s\n", cudaGetErrorString(cudaGetLastError()));
}           
       
to_host_gf(g_gauge_field, h2d_gf_d);

}









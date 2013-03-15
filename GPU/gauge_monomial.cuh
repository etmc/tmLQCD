
#ifdef MPI
extern "C"{
 #include "../xchange/xchange_gauge.h"
 #include "../geometry_eo.h"
 }
#endif //MPI

extern "C"{
 #include "../boundary.h"
 }


int blockdimgauge;
int griddimgauge;

int * nn2;
int * dev_nn2;


double2 * dev_gf_d;
dev_su3_2v_d * h2d_gf_d;
dev_su3adj* dev_df0_d;
dev_su3adj* h2d_df0_d;

__device__  int  dev_VOL2;




#include "linalg_d.cuh"
/* include double matrix on device here */
#include "Hopping_Matrix_d.cuh"



//initialize nearest-neighbour table for gpu
void initnn2(){
  int t,x,y,z,pos;
  for(t=0;t<T;t++){
   for(x=0; x<LX; x++){
    for(y=0; y<LY; y++){
     for(z=0; z<LZ; z++){   
          pos= z + LZ*(y + LY*(x + LX*t));
          //plus direction
          nn2[8*pos+0] = z + LZ*(y + LY*(x + LX*((t+1)%T)));
          nn2[8*pos+1] = z + LZ*(y + LY*((x+1)%LX + LX*t));
          nn2[8*pos+2] = z + LZ*((y+1)%LY + LY*(x + LX*t));
          nn2[8*pos+3] = (z+1)%LZ + LX*(y + LY*(x + LX*t));
          //minus direction
          if(t==0){
            nn2[8*pos+4] = z + LZ*(y + LY*(x + LX*((T-1))));
          }
          else{
            nn2[8*pos+4] = z + LZ*(y + LY*(x + LX*((t-1))));
          }
          if(x==0){
            nn2[8*pos+5] = z + LZ*(y + LY*((LX-1) + LX*t));
          }
          else{
            nn2[8*pos+5] = z + LZ*(y + LY*((x-1) + LX*t));
          }
          if(y==0){
            nn2[8*pos+6] = z + LZ*((LY-1) + LY*(x + LX*t));
          }
          else{
            nn2[8*pos+6] = z + LZ*((y-1) + LY*(x + LX*t));
          }
          if(z==0){
            nn2[8*pos+7] = (LZ-1) + LZ*(y + LY*(x + LX*t));
          }
          else{
            nn2[8*pos+7] = (z-1) + LZ*(y + LY*(x + LX*t));
          }          
        }
      }
    } 
  }
}


void initnn2_mpi(){
  int x, y, z, t, ind, j;
  for (t = 0; t < T; t++) {
    for (x = 0; x < LX; x++) {
      for (y = 0; y < LY; y++) {	
        for (z = 0; z < LZ; z++) {
            ind = g_ipt[t][x][y][z];
	    
	    //ind= Index(t,x,y,z);
	    
	    for (j = 0; j < 4; j++) { // plus direction	
                nn2[8*ind+j] = g_iup[ind][j] ;	// g_iup[ind][j] 	returns the position of the nearest neighbour of [ind] in direction +[j]
            }								//				-->  for the non-parallized code g_iup[][] maps INTERN !!
            for (j = 0; j < 4; j++) { // minus direction
              nn2[8*ind+4+j] =  g_idn[ind][j] ;	// g_idn[ind][j] 	returns the position of the nearest neighbour of [ind] in direction -[j]
            }
  }}}} // for loops

}




void su3to2vf4_d(su3** gf, dev_su3_2v_d* h2d){
  int i,j;
  int Vol;
  #ifdef MPI
   Vol = VOLUME + RAND;
  #else
   Vol = VOLUME;
  #endif
  for (i = 0; i < Vol; i++) {
   for(j=0;j<4;j++){
   //first row
    h2d[6*(4*i+j)].x = creal(gf[i][j].c00);
    h2d[6*(4*i+j)].y = cimag(gf[i][j].c00);
    h2d[6*(4*i+j)+1].x = creal(gf[i][j].c01);
    h2d[6*(4*i+j)+1].y = cimag(gf[i][j].c01);
    h2d[6*(4*i+j)+2].x = creal(gf[i][j].c02);
    h2d[6*(4*i+j)+2].y = cimag(gf[i][j].c02);      
   //second row
    h2d[6*(4*i+j)+3].x = creal(gf[i][j].c10);
    h2d[6*(4*i+j)+3].y = cimag(gf[i][j].c10);
    h2d[6*(4*i+j)+4].x = creal(gf[i][j].c11);
    h2d[6*(4*i+j)+4].y = cimag(gf[i][j].c11);
    h2d[6*(4*i+j)+5].x = creal(gf[i][j].c12);
    h2d[6*(4*i+j)+5].y = cimag(gf[i][j].c12);      
  } 
 }
}



// show nn table 
void shownn2(){
  int t,x,y,z,i,pos;
  int lx,ly,lz,lt;  
  int lx1,ly1,lz1,lt1;
//     lx = LX;
//     ly = LY;
//     lz = LZ;
//     lt =T;  
    lx1 = LX;
    ly1 = LY;
    lz1 = LZ;
    lt1 = T;
    
    lx = 0;
    ly = 0;
    lz = 0;
    lt = 0;   
    
  for(t=lt;t<lt1;t++){ 
    for(x=lx; x<lx1; x++){
      for(y=ly; y<ly1; y++){
        for(z=lz; z<lz1; z++){
          pos= g_ipt[t][x][y][z];
          printf("proc %d: pos=%d\t",g_proc_id, pos);
          for(i=0;i<8;i++){
            printf("%d ",nn2[8*pos+i]);
            //lptovec(nn[8*pos+i]);
          }
          printf("\n");
          //compare with geometry fields of hmc
          //might NOT WORK for even-odd? What are geometry indices in case of even-odd?
          //printf("%d: %d %d %d %d %d %d %d %d\n",g_ipt[t][x][y][z],g_iup[pos][0],g_iup[pos][1],g_iup[pos][2],g_iup[pos][3],g_idn[pos][0],g_idn[pos][1],g_idn[pos][2],g_idn[pos][3]);
        }
      }
    }
  }
}









void update_gpu_fields(su3** gf, su3adj** mom, int need_momenta){
  cudaError_t cudaerr;
  int i, mu,Vol;
  printf("updating double gpu fields...");
  #ifdef MPI
    Vol = VOLUME + RAND;
  #else
    Vol = VOLUME;
  #endif
  su3to2vf4_d(gf, h2d_gf_d);
  size_t dev_gfsize = 6*4*Vol * sizeof(dev_su3_2v_d);
  cudaMemcpy(dev_gf_d, h2d_gf_d, dev_gfsize, cudaMemcpyHostToDevice);
  if ((cudaerr=cudaGetLastError())!=cudaSuccess) {
     printf("%s\n", cudaGetErrorString(cudaGetLastError()));
     printf("Error code is: %f\n",cudaerr);
  }
  
  if(need_momenta){
    //note: we only need the INTERNAL lattice size here as the momentum is only updated internally
    size_t dev_momsize = 4*VOLUME * sizeof(dev_su3adj);
    for(i=0; i<VOLUME; i++){
      for(mu=0; mu<4; mu++){
        h2d_df0_d[4*i+mu].d1 = mom[i][mu].d1;
        h2d_df0_d[4*i+mu].d2 = mom[i][mu].d2;
        h2d_df0_d[4*i+mu].d3 = mom[i][mu].d3;
        h2d_df0_d[4*i+mu].d4 = mom[i][mu].d4;
        h2d_df0_d[4*i+mu].d5 = mom[i][mu].d5;
        h2d_df0_d[4*i+mu].d6 = mom[i][mu].d6;
        h2d_df0_d[4*i+mu].d7 = mom[i][mu].d7;
        h2d_df0_d[4*i+mu].d8 = mom[i][mu].d8;
      }
    }
    cudaMemcpy(dev_df0_d, h2d_df0_d, dev_momsize, cudaMemcpyHostToDevice);
  }
  
  if ((cudaerr=cudaGetLastError())!=cudaSuccess) {
     printf("%s\n", cudaGetErrorString(cudaGetLastError()));
     printf("Error code is: %f\n",cudaerr);
  }
  printf("done.\n");
}



__global__ void he_cg_init_d (int* grid, double param_kappa, double param_mu, dev_complex_d k0, dev_complex_d k1, dev_complex_d k2, dev_complex_d k3){

  dev_LX = grid[0];
  dev_LY = grid[1];
  dev_LZ = grid[2];
  dev_T = grid[3];
  dev_VOLUME = grid[4]; // grid[4] is initialized 1/2 VOLUME for eo
  dev_Offset = grid[5]; //this is the offset for the spinor fields
  dev_VOLUMEPLUSRAND = grid[5]; 
  
  kappa_d = param_kappa;
  mu_d = param_mu;
  twokappamu_d = 2.0*param_kappa*param_mu;
  
  dev_k0_d.re = k0.re;
  dev_k0_d.im = k0.im;
  dev_mk0_d.re = -k0.re;
  dev_mk0_d.im = -k0.im;
  
  dev_k1_d.re = k1.re;
  dev_k1_d.im = k1.im;
  dev_mk1_d.re = -k1.re;
  dev_mk1_d.im = -k1.im;
  
  dev_k2_d.re = k2.re;
  dev_k2_d.im = k2.im;
  dev_mk2_d.re = -k2.re;
  dev_mk2_d.im = -k2.im;
  
  dev_k3_d.re = k3.re;
  dev_k3_d.im = k3.im;
  dev_mk3_d.re = -k3.re;
  dev_mk3_d.im = -k3.im;
}


extern "C" void update_constants_d(int *grid){
  dev_complex_d h0,h1,h2,h3,mh0, mh1, mh2, mh3;
  
  h0.re = (double)creal(ka0);    h0.im = -(double)cimag(ka0);
  h1.re = (double)creal(ka1);    h1.im = -(double)cimag(ka1);
  h2.re = (double)creal(ka2);    h2.im = -(double)cimag(ka2);
  h3.re = (double)creal(ka3);    h3.im = -(double)cimag(ka3);
  
  mh0.re = -(double)creal(ka0);    mh0.im = (double)cimag(ka0);
  mh1.re = -(double)creal(ka1);    mh1.im = (double)cimag(ka1);
  mh2.re = -(double)creal(ka2);    mh2.im = (double)cimag(ka2);
  mh3.re = -(double)creal(ka3);    mh3.im = (double)cimag(ka3);
  
  // try using constant mem for kappas
  cudaMemcpyToSymbol("dev_k0c_d", &h0, sizeof(dev_complex_d)) ; 
  cudaMemcpyToSymbol("dev_k1c_d", &h1, sizeof(dev_complex_d)) ; 
  cudaMemcpyToSymbol("dev_k2c_d", &h2, sizeof(dev_complex_d)) ; 
  cudaMemcpyToSymbol("dev_k3c_d", &h3, sizeof(dev_complex_d)) ;
  
  cudaMemcpyToSymbol("dev_mk0c_d", &mh0, sizeof(dev_complex_d)) ; 
  cudaMemcpyToSymbol("dev_mk1c_d", &mh1, sizeof(dev_complex_d)) ; 
  cudaMemcpyToSymbol("dev_mk2c_d", &mh2, sizeof(dev_complex_d)) ; 
  cudaMemcpyToSymbol("dev_mk3c_d", &mh3, sizeof(dev_complex_d)) ;  
  
  he_cg_init_d<<< 1, 1 >>> (grid, (double) g_kappa, (double)(g_mu/(2.0*g_kappa)), h0,h1,h2,h3);
  // BEWARE in dev_tm_dirac_kappa we need the true mu (not 2 kappa mu!)
  
  
      #ifndef LOWOUTPUT
        if (g_proc_id == 0) {
  	  int host_check_VOLUME;
  	  int host_check_Offset;
  	  cudaMemcpyFromSymbol(&host_check_VOLUME, dev_VOLUME, sizeof(int));
  	  cudaMemcpyFromSymbol(&host_check_Offset, dev_Offset, sizeof(int));
  	  printf("\tOn device:\n");
  	  printf("\tdev_VOLUME = %i\n", host_check_VOLUME);
  	  printf("\tdev_Offset = %i\n", host_check_Offset);
  	}
    #endif
  
}




//only update the gaugefield on device (not the momentum field)
//accessible from outside
extern "C" void update_gpu_gf_d(su3** gf){
  cudaError_t cudaerr;
  int Vol;
  #ifdef MPI
    Vol = VOLUME + RAND;
  #else
    Vol = VOLUME;
  #endif 
  
  su3to2vf4_d(gf, h2d_gf_d);
  size_t dev_gfsize = 6*4*Vol * sizeof(dev_su3_2v_d);
  cudaMemcpy(dev_gf_d, h2d_gf_d, dev_gfsize, cudaMemcpyHostToDevice);
  if ((cudaerr=cudaGetLastError())!=cudaSuccess) {
     printf("%s\n", cudaGetErrorString(cudaGetLastError()));
     printf("Error code is: %f\n",cudaerr);
  }
  
 
}



extern "C" void init_gpu_fields(int need_momenta){
  cudaError_t cudaerr;
  int Vol;
  #ifdef MPI
    Vol = VOLUME + RAND;
  #else
    Vol = VOLUME;
  #endif  
  
  size_t dev_gfsize = 6*4*Vol * sizeof(dev_su3_2v_d);
  if((cudaerr=cudaMalloc((void **) &dev_gf_d, dev_gfsize)) != cudaSuccess){
    printf("Error in init_gpu_fields(): Memory allocation of double gauge field failed. Aborting...\n");
    exit(200);
  }   // Allocate array on device
  else{
    printf("Allocated double gauge field on device\n");
  }  

  if((h2d_gf_d = (dev_su3_2v_d *)malloc(dev_gfsize))==(dev_su3_2v_d*) NULL){
    printf("Error in init_gpu_fields(): Memory allocation of double h2d field failed. Aborting...\n");
    exit(200);
  }

  if(need_momenta){
    //only need INTERNAL here
    size_t dev_momsize = 4*VOLUME * sizeof(dev_su3adj);
    if((cudaerr=cudaMalloc((void **) &dev_df0_d, dev_momsize)) != cudaSuccess){
      printf("Error in init_gpu_fields(): Memory allocation of double momentum field failed. Aborting...\n");
      exit(200);
    }   // Allocate array on device
    else{
      printf("Allocated double momentum field on device\n");
    } 

     if((h2d_df0_d = (dev_su3adj *)malloc(dev_momsize))==(dev_su3adj *) NULL){
      printf("Error in init_gpu_fields(): Memory allocation of double momentum h2d field failed. Aborting...\n");
      exit(200);
    }  
  }  
 
  // initialize partition
  int gridsize;
  if(VOLUME%BLOCKGAUGE != 0){
   printf("Error: VOLUME is not a multiple of BLOCKGAUGE. Aborting...\n");
   exit(200);
  }
  blockdimgauge = BLOCKGAUGE;
  if( VOLUME >= BLOCKGAUGE){
   gridsize = (int) VOLUME/BLOCKGAUGE;
  }
  else{
   gridsize=1;
  }
  griddimgauge = gridsize; 
  printf("Have initialized double fields on device.\n");
  
  
  int gvol = VOLUME;
  if((cudaerr=cudaMemcpyToSymbol(dev_VOL2, &(gvol), sizeof(int)))!=cudaSuccess){
    printf("Error in init_gpu_fields(): Could not copy dev_VOL2 to device. Aborting...\n");
    exit(200);
  } 
  
  int host_check_VOL2;
  cudaMemcpyFromSymbol(&host_check_VOL2, dev_VOL2, sizeof(int)); 
  printf("\tOn device:\n");
  printf("\tdev_VOL2 = %i\n", host_check_VOL2);  
  
  //grid 
  size_t nnsize = 8*VOLUME*sizeof(int);
  nn2 = (int *) malloc(nnsize);
  cudaMalloc((void **) &dev_nn2, nnsize);
  #ifdef MPI
    initnn2_mpi();
    //shownn2(); //for tests
  #else
    initnn2();
  #endif
  if((cudaerr=cudaMemcpy(dev_nn2, nn2, nnsize, cudaMemcpyHostToDevice)) !=cudaSuccess){
    printf("Error in init_gpu_fields(): Could not transfer nn2 fields to device. Aborting...\n");
    exit(200);
  } 
  
  //store them immediately
  update_gpu_fields(g_gauge_field, df0,need_momenta);
  
  if ((cudaerr=cudaPeekAtLastError())!=cudaSuccess) {
    printf("%s\n", cudaGetErrorString(cudaPeekAtLastError()));
    printf("Error code is: %f\n",cudaerr);
  } 
  #ifdef MPI
   MPI_Barrier(g_cart_grid);
  #endif
  
  
  //init blas kernels for double
  #ifdef GPU_DOUBLE
    //initialize blas with 
    init_blas_d(Vol/2);
  #endif
}



extern "C" void finalize_gpu_fields(){
  cudaFree(dev_gf_d);
  cudaFree(dev_df0_d);
  cudaFree(dev_nn2);
  free(h2d_gf_d);
  free(h2d_df0_d);
  free(nn2);
  #ifdef GPU_DOUBLE
    finalize_blas_d();
  #endif
}






void to_host_mom(hamiltonian_field_t * const hf){
  cudaError_t cudaerr;
  int i, mu;
  size_t dev_momsize = 4*VOLUME * sizeof(dev_su3adj);
  cudaMemcpy(h2d_df0_d, dev_df0_d, dev_momsize, cudaMemcpyDeviceToHost);
  if ((cudaerr=cudaGetLastError())!=cudaSuccess) {
     printf("%s\n", cudaGetErrorString(cudaGetLastError()));
     printf("Error code is: %f\n",cudaerr);
  }
  
  for(i=0; i<VOLUME; i++){
    for(mu=0; mu<4; mu++){
      hf->derivative[i][mu].d1 = h2d_df0_d[4*i+mu].d1;
      hf->derivative[i][mu].d2 = h2d_df0_d[4*i+mu].d2;
      hf->derivative[i][mu].d3 = h2d_df0_d[4*i+mu].d3;
      hf->derivative[i][mu].d4 = h2d_df0_d[4*i+mu].d4;
      hf->derivative[i][mu].d5 = h2d_df0_d[4*i+mu].d5;
      hf->derivative[i][mu].d6 = h2d_df0_d[4*i+mu].d6;
      hf->derivative[i][mu].d7 = h2d_df0_d[4*i+mu].d7;
      hf->derivative[i][mu].d8 = h2d_df0_d[4*i+mu].d8;
    }
  }
}



void to_host_mom_field(su3adj** mom, dev_su3adj* dev_mom){
  cudaError_t cudaerr;
  int i, mu;
  size_t dev_momsize = 4*VOLUME * sizeof(dev_su3adj);
  cudaMemcpy(h2d_df0_d, dev_mom, dev_momsize, cudaMemcpyDeviceToHost);
  if ((cudaerr=cudaGetLastError())!=cudaSuccess) {
     printf("%s\n", cudaGetErrorString(cudaGetLastError()));
     printf("Error code is: %f\n",cudaerr);
  }
  
  for(i=0; i<VOLUME; i++){
    for(mu=0; mu<4; mu++){
      mom[i][mu].d1 = h2d_df0_d[4*i+mu].d1;
      mom[i][mu].d2 = h2d_df0_d[4*i+mu].d2;
      mom[i][mu].d3 = h2d_df0_d[4*i+mu].d3;
      mom[i][mu].d4 = h2d_df0_d[4*i+mu].d4;
      mom[i][mu].d5 = h2d_df0_d[4*i+mu].d5;
      mom[i][mu].d6 = h2d_df0_d[4*i+mu].d6;
      mom[i][mu].d7 = h2d_df0_d[4*i+mu].d7;
      mom[i][mu].d8 = h2d_df0_d[4*i+mu].d8;
    }
  }
  if ((cudaerr=cudaPeekAtLastError())!=cudaSuccess) {   
    printf("%s\n", cudaGetErrorString(cudaPeekAtLastError()));
  }
}


extern "C" void compare_momenta(){
int i;
for(i=VOLUME-1; i<VOLUME; i++){
      printf("*******\n");
      printf("%.16e\t%.16e\n", df0[i][0].d1, h2d_df0_d[4*i+0].d1);
      printf("%.16e\t%.16e\n", df0[i][1].d1, h2d_df0_d[4*i+1].d1);
      printf("%.16e\t%.16e\n", df0[i][2].d1, h2d_df0_d[4*i+2].d1);
      printf("%.16e\t%.16e\n", df0[i][3].d1, h2d_df0_d[4*i+3].d1);
      printf("*******\n");
      printf("%.16e\t%.16e\n", df0[i][0].d2, h2d_df0_d[4*i+0].d2);
      printf("%.16e\t%.16e\n", df0[i][1].d2, h2d_df0_d[4*i+1].d2);
      printf("%.16e\t%.16e\n", df0[i][2].d2, h2d_df0_d[4*i+2].d2);
      printf("%.16e\t%.16e\n", df0[i][3].d2, h2d_df0_d[4*i+3].d2);
      printf("*******\n");
      printf("%.16e\t%.16e\n", df0[i][0].d3, h2d_df0_d[4*i+0].d3);
      printf("%.16e\t%.16e\n", df0[i][1].d3, h2d_df0_d[4*i+1].d3);
      printf("%.16e\t%.16e\n", df0[i][2].d3, h2d_df0_d[4*i+2].d3);
      printf("%.16e\t%.16e\n", df0[i][3].d3, h2d_df0_d[4*i+3].d3);
      printf("*******\n");
      printf("%.16e\t%.16e\n", df0[i][0].d4, h2d_df0_d[4*i+0].d4);
      printf("%.16e\t%.16e\n", df0[i][1].d4, h2d_df0_d[4*i+1].d4);
      printf("%.16e\t%.16e\n", df0[i][2].d4, h2d_df0_d[4*i+2].d4);
      printf("%.16e\t%.16e\n", df0[i][3].d4, h2d_df0_d[4*i+3].d4);
      printf("*******\n");
      printf("%.16e\t%.16e\n", df0[i][0].d5, h2d_df0_d[4*i+0].d5);
      printf("%.16e\t%.16e\n", df0[i][1].d5, h2d_df0_d[4*i+1].d5);
      printf("%.16e\t%.16e\n", df0[i][2].d5, h2d_df0_d[4*i+2].d5);
      printf("%.16e\t%.16e\n", df0[i][3].d5, h2d_df0_d[4*i+3].d5);
      printf("*******\n");
      printf("%.16e\t%.16e\n", df0[i][0].d6, h2d_df0_d[4*i+0].d6);
      printf("%.16e\t%.16e\n", df0[i][1].d6, h2d_df0_d[4*i+1].d6);
      printf("%.16e\t%.16e\n", df0[i][2].d6, h2d_df0_d[4*i+2].d6);
      printf("%.16e\t%.16e\n", df0[i][3].d6, h2d_df0_d[4*i+3].d6);
      printf("*******\n");
      printf("%.16e\t%.16e\n", df0[i][0].d7, h2d_df0_d[4*i+0].d7);
      printf("%.16e\t%.16e\n", df0[i][1].d7, h2d_df0_d[4*i+1].d7);
      printf("%.16e\t%.16e\n", df0[i][2].d7, h2d_df0_d[4*i+2].d7);
      printf("%.16e\t%.16e\n", df0[i][3].d7, h2d_df0_d[4*i+3].d7);
      printf("*******\n");
      printf("%.16e\t%.16e\n", df0[i][0].d8, h2d_df0_d[4*i+0].d8);
      printf("%.16e\t%.16e\n", df0[i][1].d8, h2d_df0_d[4*i+1].d8);
      printf("%.16e\t%.16e\n", df0[i][2].d8, h2d_df0_d[4*i+2].d8);
      printf("%.16e\t%.16e\n", df0[i][3].d8, h2d_df0_d[4*i+3].d8);
      printf("*******\n");
    
    }
}









__global__ void dev_gauge_derivative(int * dev_nn, dev_su3_2v_d * gf, dev_su3adj* mom, double c){
  
  int pos, ix, mu, nu; /* pos = basepoint of loop */
  int help, gather, old, thispos, nextpos;
  dev_su3_d M[2];
  dev_su3_d newM;
  dev_su3_d w;
  __shared__ dev_su3adj newmom[BLOCKGAUGE];

  pos= threadIdx.x + blockDim.x*blockIdx.x;  
  ix= threadIdx.x;
  if(pos < dev_VOL2){
  
  
  // GAUGE STAPLES
  for(mu=0; mu<4; mu++){
   dev_zerosu3_d(&(w));
   
   for(nu=0; nu<4; nu++){
     if(mu!=nu){
     
     /*forward staples*/
     thispos = pos;
     nextpos = thispos;
     //init old matrix
     gather=1;
     old=0;
     dev_unitsu3_d(&(M[old]));
   

     //load next U
     dev_get_matrix(gf, nextpos, nu, &(newM));
     //multiply
     dev_su3_ti_su3_d(&(M[gather]), &(M[old]), &(newM));
     //get new nextpos
     thispos=nextpos;
     nextpos=dev_nn[8*thispos+nu];
     //swap old and gather
     help=gather;
     gather=old;
     old=help;
   
 
     //load next U
     dev_get_matrix(gf, nextpos, mu, &(newM));
     //multiply
     dev_su3_ti_su3_d(&(M[gather]), &(M[old]), &(newM));
     //get new nextpos
     thispos=nextpos;
     nextpos=dev_nn[8*thispos+mu];
     //swap old and gather
     help=gather;
     gather=old;
     old=help;
   
  
     //get new nextpos
     //we have to do this BEFORE loading the new U matrix to keep 
     //gauge invariance
     thispos=nextpos;
     nextpos=dev_nn[8*thispos+nu+4]; // NEGATIVE DIR !!! 
     //load next U
     dev_get_matrix_dagger(gf, nextpos, nu, &(newM));
     //multiply
     dev_su3_ti_su3_d(&(M[gather]), &(M[old]), &(newM));
     //swap old and gather
     help=gather;
     gather=old;
     old=help;
     
     //accumulate
     dev_su3_add_d(&w, &(M[old]));
     
     
     
  /*backward staples*/
     thispos = pos;
     nextpos = thispos;
     //init old matrix
     gather=1;
     old=0;
     dev_unitsu3_d(&(M[old]));
   
   
     //get new nextpos
     //we have to do this BEFORE loading the new U matrix to keep 
     //gauge invariance
     thispos=nextpos;
     nextpos=dev_nn[8*thispos+nu+4]; // NEGATIVE DIR !!! 
     //load next U
     dev_get_matrix_dagger(gf, nextpos, nu, &(newM));
     //multiply
     dev_su3_ti_su3_d(&(M[gather]), &(M[old]), &(newM));
     //swap old and gather
     help=gather;
     gather=old;
     old=help;


     //load next U
     dev_get_matrix(gf, nextpos, mu, &(newM));
     //multiply
     dev_su3_ti_su3_d(&(M[gather]), &(M[old]), &(newM));
     //get new nextpos
     thispos=nextpos;
     nextpos=dev_nn[8*thispos+mu];
     //swap old and gather
     help=gather;
     gather=old;
     old=help;
   
   
     //load next U
     dev_get_matrix(gf, nextpos, nu, &(newM));
     //multiply
     dev_su3_ti_su3_d(&(M[gather]), &(M[old]), &(newM));
     //get new nextpos
     thispos=nextpos;
     nextpos=dev_nn[8*thispos+nu];
     //swap old and gather
     help=gather;
     gather=old;
     old=help; 

     //accumulate
     dev_su3_add_d(&w, &(M[old]));     
     
     }// mu!=nu
   }//nu
   
   
   dev_get_matrix(gf, pos, mu, &(newM));
   dev_su3_ti_su3d_d(&(M[gather]), &(newM), &(w) );
   dev_tr_lambda(&(newmom[ix]), &(M[gather]));
   
   mom[4*pos+mu].d1 += c*newmom[ix].d1;
   mom[4*pos+mu].d2 += c*newmom[ix].d2;
   mom[4*pos+mu].d3 += c*newmom[ix].d3;
   mom[4*pos+mu].d4 += c*newmom[ix].d4;
   mom[4*pos+mu].d5 += c*newmom[ix].d5;
   mom[4*pos+mu].d6 += c*newmom[ix].d6;
   mom[4*pos+mu].d7 += c*newmom[ix].d7;
   mom[4*pos+mu].d8 += c*newmom[ix].d8;
   

 
  }//mu
  }//pos < dev_VOL2
  
}







 
__global__ void dev_rect_deri_one(int * dev_nn, dev_su3_2v_d * gf, dev_su3adj* mom,  double c, int mu){
  
  int pos, nu; // pos = basepoint of loop 
  int help, gather, old, thispos, nextpos;
  dev_su3_d M[2];
  dev_su3_d newM;
  dev_su3_d w;
  int ix;

  pos= threadIdx.x + blockDim.x*blockIdx.x;  
  ix= threadIdx.x;

  if(pos < dev_VOL2){
  
  
  
  //for(mu=0; mu<4; mu++){
   dev_zerosu3_d(&(w));
   // RECTANGLE STAPLES
    dev_zerosu3_d(&(w));
    for(nu=0; nu<4; nu++){
     if(mu!=nu){
     
  //    __ __
  //  ||__ __|
  // 
     
     
     // forward staples
     thispos = pos;
     nextpos=pos;
     //init old matrix
     gather=1;
     old=0;
     dev_unitsu3_d(&(M[old]));
   
   //two steps in "nu"
     //load next U
     dev_get_matrix(gf, nextpos, nu, &(newM));
     //multiply
     dev_su3_ti_su3_d(&(M[gather]), &(M[old]), &(newM));
     //get new nextpos
     thispos=nextpos;
     nextpos=dev_nn[8*thispos+nu];
     //swap old and gather
     help=gather;
     gather=old;
     old=help;
     
     //load next U
     dev_get_matrix(gf, nextpos, nu, &(newM));
     //multiply
     dev_su3_ti_su3_d(&(M[gather]), &(M[old]), &(newM));
     //get new nextpos
     thispos=nextpos;
     nextpos=dev_nn[8*thispos+nu];
     //swap old and gather
     help=gather;
     gather=old;
     old=help;     
     
   
   
   //one step in "mu"
     //load next U
     dev_get_matrix(gf, nextpos, mu, &(newM));
     //multiply
     dev_su3_ti_su3_d(&(M[gather]), &(M[old]), &(newM));
     //get new nextpos
     thispos=nextpos;
     nextpos=dev_nn[8*thispos+mu];
     //swap old and gather
     help=gather;
     gather=old;
     old=help;     
   
     
   
   //two steps BACK in "nu"
   
     //get new nextpos
     //we have to do this BEFORE loading the new U matrix to keep 
     //gauge invariance
     thispos=nextpos;
     nextpos=dev_nn[8*thispos+nu+4]; // NEGATIVE DIR !!! 
     //load next U
     dev_get_matrix_dagger(gf, nextpos, nu, &(newM));
     //multiply
     dev_su3_ti_su3_d(&(M[gather]), &(M[old]), &(newM));
     //swap old and gather
     help=gather;
     gather=old;
     old=help;
     
     //get new nextpos
     //we have to do this BEFORE loading the new U matrix to keep 
     //gauge invariance
     thispos=nextpos;
     nextpos=dev_nn[8*thispos+nu+4]; // NEGATIVE DIR !!! 
     //load next U
     dev_get_matrix_dagger(gf, nextpos, nu, &(newM));
     //multiply
     dev_su3_ti_su3_d(&(M[gather]), &(M[old]), &(newM));
     //swap old and gather
     help=gather;
     gather=old;
     old=help; 
     //accumulate
     dev_su3_add_d(&w, &(M[old]));
     
   
  //   __ __
  //  |__ __||
  // 
     
   
   
   // backward staples
     thispos = pos;
     nextpos=pos;
     //init old matrix
     gather=1;
     old=0;
     dev_unitsu3_d(&(M[old]));
   
   //two steps BACK in "nu"
     //get new nextpos
     //we have to do this BEFORE loading the new U matrix to keep 
     //gauge invariance
     thispos=nextpos;
     nextpos=dev_nn[8*thispos+nu+4]; // NEGATIVE DIR !!! 
     //load next U
     dev_get_matrix_dagger(gf, nextpos, nu, &(newM));
     //multiply
     dev_su3_ti_su3_d(&(M[gather]), &(M[old]), &(newM));
     //swap old and gather
     help=gather;
     gather=old;
     old=help;
     
     //get new nextpos
     //we have to do this BEFORE loading the new U matrix to keep 
     //gauge invariance
     thispos=nextpos;
     nextpos=dev_nn[8*thispos+nu+4]; // NEGATIVE DIR !!! 
     //load next U
     dev_get_matrix_dagger(gf, nextpos, nu, &(newM));
     //multiply
     dev_su3_ti_su3_d(&(M[gather]), &(M[old]), &(newM));
     //swap old and gather
     help=gather;
     gather=old;
     old=help; 
     
     
    //one step in "mu"
     //load next U
     dev_get_matrix(gf, nextpos, mu, &(newM));
     //multiply
     dev_su3_ti_su3_d(&(M[gather]), &(M[old]), &(newM));
     //get new nextpos
     thispos=nextpos;
     nextpos=dev_nn[8*thispos+mu];
     //swap old and gather
     help=gather;
     gather=old;
     old=help;
         
         
             
   //two steps in "nu"
     //load next U
     dev_get_matrix(gf, nextpos, nu, &(newM));
     //multiply
     dev_su3_ti_su3_d(&(M[gather]), &(M[old]), &(newM));
     //get new nextpos
     thispos=nextpos;
     nextpos=dev_nn[8*thispos+nu];
     //swap old and gather
     help=gather;
     gather=old;
     old=help;
     
     //load next U
     dev_get_matrix(gf, nextpos, nu, &(newM));
     //multiply
     dev_su3_ti_su3_d(&(M[gather]), &(M[old]), &(newM));
     //get new nextpos
     thispos=nextpos;
     nextpos=dev_nn[8*thispos+nu];
     //swap old and gather
     help=gather;
     gather=old;
     old=help;     
   
     //accumulate
     dev_su3_add_d(&w, &(M[old]));
     
     }// mu!=nu
   }//nu
   
   

   __shared__ dev_su3adj newmom[BLOCKGAUGE];
   dev_get_matrix(gf, pos, mu, &(newM));
   dev_su3_ti_su3d_d(&(M[old]), &(newM), &(w) );
   dev_tr_lambda(&(newmom[ix]), &(M[old]));
   
   mom[4*pos+mu].d1 += c*newmom[ix].d1;
   mom[4*pos+mu].d2 += c*newmom[ix].d2;
   mom[4*pos+mu].d3 += c*newmom[ix].d3;
   mom[4*pos+mu].d4 += c*newmom[ix].d4;
   mom[4*pos+mu].d5 += c*newmom[ix].d5;
   mom[4*pos+mu].d6 += c*newmom[ix].d6;
   mom[4*pos+mu].d7 += c*newmom[ix].d7;
   mom[4*pos+mu].d8 += c*newmom[ix].d8;


  //}//mu
  }//pos < dev_VOL2 
  
}






__global__ void dev_rect_deri_two(int * dev_nn, dev_su3_2v_d * gf, dev_su3adj* mom, double c, int mu){
  
  int pos, nu; // pos = basepoint of loop 
  int help, gather, old, thispos, nextpos;
  dev_su3_d M[2];
  dev_su3_d newM;
  dev_su3_d w;
  int ix;
  pos= threadIdx.x + blockDim.x*blockIdx.x;  
  ix = threadIdx.x;
  
  if(pos < dev_VOL2){
   
  //for(mu=0; mu<4; mu++){
   dev_zerosu3_d(&(w));
   
  //    __
  //   |  |
  //  ||__|
  //   
   
   
   for(nu=0; nu<4; nu++){
     if(mu!=nu){     
     // forward staples
     thispos = pos;
     nextpos=pos;
     //init old matrix
     gather=1;
     old=0;
     dev_unitsu3_d(&(M[old]));
   
   //one step in "nu"
     //load next U
     dev_get_matrix(gf, nextpos, nu, &(newM));
     //multiply
     dev_su3_ti_su3_d(&(M[gather]), &(M[old]), &(newM));
     //get new nextpos
     thispos=nextpos;
     nextpos=dev_nn[8*thispos+nu];
     //swap old and gather
     help=gather;
     gather=old;
     old=help;
     
   //two steps in "mu"  
     //load next U
     dev_get_matrix(gf, nextpos, mu, &(newM));
     //multiply
     dev_su3_ti_su3_d(&(M[gather]), &(M[old]), &(newM));
     //get new nextpos
     thispos=nextpos;
     nextpos=dev_nn[8*thispos+mu];
     //swap old and gather
     help=gather;
     gather=old;
     old=help;     
     
     //load next U
     dev_get_matrix(gf, nextpos, mu, &(newM));
     //multiply
     dev_su3_ti_su3_d(&(M[gather]), &(M[old]), &(newM));
     //get new nextpos
     thispos=nextpos;
     nextpos=dev_nn[8*thispos+mu];
     //swap old and gather
     help=gather;
     gather=old;
     old=help;     
     
   
   //one step in NEGATIVE "nu"
     //get new nextpos
     //we have to do this BEFORE loading the new U matrix to keep 
     //gauge invariance
     thispos=nextpos;
     nextpos=dev_nn[8*thispos+nu+4]; // NEGATIVE DIR !!! 
     //load next U
     dev_get_matrix_dagger(gf, nextpos, nu, &(newM));
     //multiply
     dev_su3_ti_su3_d(&(M[gather]), &(M[old]), &(newM));
     //swap old and gather
     help=gather;
     gather=old;
     old=help;
     
     

    //one step in NEGATIVE "mu"
     //get new nextpos
     //we have to do this BEFORE loading the new U matrix to keep 
     //gauge invariance
     thispos=nextpos;
     nextpos=dev_nn[8*thispos+mu+4]; // NEGATIVE DIR !!! 
     //load next U
     dev_get_matrix_dagger(gf, nextpos, mu, &(newM));
     //multiply
     dev_su3_ti_su3_d(&(M[gather]), &(M[old]), &(newM));
     //swap old and gather
     help=gather;
     gather=old;
     old=help; 
     
     //accumulate
     dev_su3_add_d(&w, &(M[old]));

     }// mu!=nu
   }//nu

   __shared__ dev_su3adj newmom[BLOCKGAUGE];
   dev_get_matrix(gf, pos, mu, &(newM));
   dev_su3_ti_su3d_d(&(M[old]), &(newM), &(w) );
   dev_tr_lambda(&(newmom[ix]), &(M[old]));
   
   mom[4*pos+mu].d1 += c*newmom[ix].d1;
   mom[4*pos+mu].d2 += c*newmom[ix].d2;
   mom[4*pos+mu].d3 += c*newmom[ix].d3;
   mom[4*pos+mu].d4 += c*newmom[ix].d4;
   mom[4*pos+mu].d5 += c*newmom[ix].d5;
   mom[4*pos+mu].d6 += c*newmom[ix].d6;
   mom[4*pos+mu].d7 += c*newmom[ix].d7;
   mom[4*pos+mu].d8 += c*newmom[ix].d8;


  //    
  //    __
  //  ||  |
  //   |__|

  dev_zerosu3_d(&(w));
  for(nu=0; nu<4; nu++){
    if(mu!=nu){
   //backward staples
    thispos = pos;
     nextpos=pos;
     //init old matrix
     gather=1;
     old=0;
     dev_unitsu3_d(&(M[old]));
   
   //one step in NEGATIVE "mu"
     //get new nextpos
     //we have to do this BEFORE loading the new U matrix to keep 
     //gauge invariance
     thispos=nextpos;
     nextpos=dev_nn[8*thispos+mu+4]; // NEGATIVE DIR !!! 
     //load next U
     dev_get_matrix_dagger(gf, nextpos, mu, &(newM));
     //multiply
     dev_su3_ti_su3_d(&(M[gather]), &(M[old]), &(newM));
     //swap old and gather
     help=gather;
     gather=old;
     old=help;
          

    //one step in "nu"
     //load next U
     dev_get_matrix(gf, nextpos, nu, &(newM));
     //multiply
     dev_su3_ti_su3_d(&(M[gather]), &(M[old]), &(newM));
     //get new nextpos
     thispos=nextpos;
     nextpos=dev_nn[8*thispos+nu];
     //swap old and gather
     help=gather;
     gather=old;
     old=help;  
    
    
    
    //two steps in "mu"
     //load next U
     dev_get_matrix(gf, nextpos, mu, &(newM));
     //multiply
     dev_su3_ti_su3_d(&(M[gather]), &(M[old]), &(newM));
     //get new nextpos
     thispos=nextpos;
     nextpos=dev_nn[8*thispos+mu];
     //swap old and gather
     help=gather;
     gather=old;
     old=help;     
    
     //load next U
     dev_get_matrix(gf, nextpos, mu, &(newM));
     //multiply
     dev_su3_ti_su3_d(&(M[gather]), &(M[old]), &(newM));
     //get new nextpos
     thispos=nextpos;
     nextpos=dev_nn[8*thispos+mu];
     //swap old and gather
     help=gather;
     gather=old;
     old=help;  
         
         
   //one step in NEGATIVE "nu"
     //get new nextpos
     //we have to do this BEFORE loading the new U matrix to keep 
     //gauge invariance
     thispos=nextpos;
     nextpos=dev_nn[8*thispos+nu+4]; // NEGATIVE DIR !!! 
     //load next U
     dev_get_matrix_dagger(gf, nextpos, nu, &(newM));
     //multiply
     dev_su3_ti_su3_d(&(M[gather]), &(M[old]), &(newM));
     //swap old and gather
     help=gather;
     gather=old;
     old=help;    
    
     //accumulate
     dev_su3_add_d(&w, &(M[old]));
    
     }// mu!=nu
   }//nu
   
   

   dev_get_matrix(gf, pos, mu, &(newM));
   dev_su3_ti_su3d_d(&(M[old]), &(newM), &(w) );
   dev_tr_lambda(&(newmom[ix]), &(M[old]));
   
   mom[4*pos+mu].d1 += c*newmom[ix].d1;
   mom[4*pos+mu].d2 += c*newmom[ix].d2;
   mom[4*pos+mu].d3 += c*newmom[ix].d3;
   mom[4*pos+mu].d4 += c*newmom[ix].d4;
   mom[4*pos+mu].d5 += c*newmom[ix].d5;
   mom[4*pos+mu].d6 += c*newmom[ix].d6;
   mom[4*pos+mu].d7 += c*newmom[ix].d7;
   mom[4*pos+mu].d8 += c*newmom[ix].d8;

  //}//mu
  }//pos < dev_VOL2
  
}







__global__ void dev_rect_deri_three(int * dev_nn, dev_su3_2v_d * gf, dev_su3adj* mom, double c, int mu){
  
  int pos, nu; // pos = basepoint of loop 
  int help, gather, old, thispos, nextpos;
  dev_su3_d M[2];
  dev_su3_d newM;
  dev_su3_d w;
  int ix;
  pos= threadIdx.x + blockDim.x*blockIdx.x;  
  ix = threadIdx.x;
  
  if(pos < dev_VOL2){
   
  //for(mu=0; mu<4; mu++){
   dev_zerosu3_d(&(w));
   
  //   __
  //  |  |
  //  |__||
  //   
   
   
   for(nu=0; nu<4; nu++){
     if(mu!=nu){     
     // forward staples
     thispos = pos;
     nextpos=pos;
     //init old matrix
     gather=1;
     old=0;
     dev_unitsu3_d(&(M[old]));
   
     
    //one step in NEGATIVE "nu"
     //get new nextpos
     //we have to do this BEFORE loading the new U matrix to keep 
     //gauge invariance
     thispos=nextpos;
     nextpos=dev_nn[8*thispos+nu+4]; // NEGATIVE DIR !!! 
     //load next U
     dev_get_matrix_dagger(gf, nextpos, nu, &(newM));
     //multiply
     dev_su3_ti_su3_d(&(M[gather]), &(M[old]), &(newM));
     //swap old and gather
     help=gather;
     gather=old;
     old=help;
     
   //two steps in "mu"
     //load next U
     dev_get_matrix(gf, nextpos, mu, &(newM));
     //multiply
     dev_su3_ti_su3_d(&(M[gather]), &(M[old]), &(newM));
     //get new nextpos
     thispos=nextpos;
     nextpos=dev_nn[8*thispos+mu];
     //swap old and gather
     help=gather;
     gather=old;
     old=help;   
   
     //load next U
     dev_get_matrix(gf, nextpos, mu, &(newM));
     //multiply
     dev_su3_ti_su3_d(&(M[gather]), &(M[old]), &(newM));
     //get new nextpos
     thispos=nextpos;
     nextpos=dev_nn[8*thispos+mu];
     //swap old and gather
     help=gather;
     gather=old;
     old=help;   
       
     
    //one step in "nu"
     //load next U
     dev_get_matrix(gf, nextpos, nu, &(newM));
     //multiply
     dev_su3_ti_su3_d(&(M[gather]), &(M[old]), &(newM));
     //get new nextpos
     thispos=nextpos;
     nextpos=dev_nn[8*thispos+nu];
     //swap old and gather
     help=gather;
     gather=old;
     old=help;
     
    //one step in NEGATIVE "mu"
     //get new nextpos
     //we have to do this BEFORE loading the new U matrix to keep 
     //gauge invariance
     thispos=nextpos;
     nextpos=dev_nn[8*thispos+mu+4]; // NEGATIVE DIR !!! 
     //load next U
     dev_get_matrix_dagger(gf, nextpos, mu, &(newM));
     //multiply
     dev_su3_ti_su3_d(&(M[gather]), &(M[old]), &(newM));
     //swap old and gather
     help=gather;
     gather=old;
     old=help;     
     
     
     //accumulate
     dev_su3_add_d(&w, &(M[old]));

     }// mu!=nu
   }//nu

   __shared__ dev_su3adj newmom[BLOCKGAUGE];
   dev_get_matrix(gf, pos, mu, &(newM));
   dev_su3_ti_su3d_d(&(M[old]), &(newM), &(w) );
   dev_tr_lambda(&(newmom[ix]), &(M[old]));
   
   mom[4*pos+mu].d1 += c*newmom[ix].d1;
   mom[4*pos+mu].d2 += c*newmom[ix].d2;
   mom[4*pos+mu].d3 += c*newmom[ix].d3;
   mom[4*pos+mu].d4 += c*newmom[ix].d4;
   mom[4*pos+mu].d5 += c*newmom[ix].d5;
   mom[4*pos+mu].d6 += c*newmom[ix].d6;
   mom[4*pos+mu].d7 += c*newmom[ix].d7;
   mom[4*pos+mu].d8 += c*newmom[ix].d8;


  //  __
  // |  ||
  // |__|  

  dev_zerosu3_d(&(w));
  for(nu=0; nu<4; nu++){
    if(mu!=nu){
   //backward staples
     thispos = pos;
     nextpos=pos;
     //init old matrix
     gather=1;
     old=0;
     dev_unitsu3_d(&(M[old]));
   
   
   //one step in NEGATIVE "mu"
     //get new nextpos
     //we have to do this BEFORE loading the new U matrix to keep 
     //gauge invariance
     thispos=nextpos;
     nextpos=dev_nn[8*thispos+mu+4]; // NEGATIVE DIR !!! 
     //load next U
     dev_get_matrix_dagger(gf, nextpos, mu, &(newM));
     //multiply
     dev_su3_ti_su3_d(&(M[gather]), &(M[old]), &(newM));
     //swap old and gather
     help=gather;
     gather=old;
     old=help;
          
     
   //one step in NEGATIVE "nu"
     //get new nextpos
     //we have to do this BEFORE loading the new U matrix to keep 
     //gauge invariance
     thispos=nextpos;
     nextpos=dev_nn[8*thispos+nu+4]; // NEGATIVE DIR !!! 
     //load next U
     dev_get_matrix_dagger(gf, nextpos, nu, &(newM));
     //multiply
     dev_su3_ti_su3_d(&(M[gather]), &(M[old]), &(newM));
     //swap old and gather
     help=gather;
     gather=old;
     old=help;

    //two steps in "mu"
    //load next U
     dev_get_matrix(gf, nextpos, mu, &(newM));
     //multiply
     dev_su3_ti_su3_d(&(M[gather]), &(M[old]), &(newM));
     //get new nextpos
     thispos=nextpos;
     nextpos=dev_nn[8*thispos+mu];
     //swap old and gather
     help=gather;
     gather=old;
     old=help; 
    
    //load next U
     dev_get_matrix(gf, nextpos, mu, &(newM));
     //multiply
     dev_su3_ti_su3_d(&(M[gather]), &(M[old]), &(newM));
     //get new nextpos
     thispos=nextpos;
     nextpos=dev_nn[8*thispos+mu];
     //swap old and gather
     help=gather;
     gather=old;
     old=help; 
        
        
    //one step in "nu"
     //load next U
     dev_get_matrix(gf, nextpos, nu, &(newM));
     //multiply
     dev_su3_ti_su3_d(&(M[gather]), &(M[old]), &(newM));
     //get new nextpos
     thispos=nextpos;
     nextpos=dev_nn[8*thispos+nu];
     //swap old and gather
     help=gather;
     gather=old;
     old=help; 

    
     //accumulate
     dev_su3_add_d(&w, &(M[old]));
    
     }// mu!=nu
   }//nu
   
   

   dev_get_matrix(gf, pos, mu, &(newM));
   dev_su3_ti_su3d_d(&(M[old]), &(newM), &(w) );
   dev_tr_lambda(&(newmom[ix]), &(M[old]));
   
   mom[4*pos+mu].d1 += c*newmom[ix].d1;
   mom[4*pos+mu].d2 += c*newmom[ix].d2;
   mom[4*pos+mu].d3 += c*newmom[ix].d3;
   mom[4*pos+mu].d4 += c*newmom[ix].d4;
   mom[4*pos+mu].d5 += c*newmom[ix].d5;
   mom[4*pos+mu].d6 += c*newmom[ix].d6;
   mom[4*pos+mu].d7 += c*newmom[ix].d7;
   mom[4*pos+mu].d8 += c*newmom[ix].d8;

  //}//mu
  }//pos < dev_VOL2
  
}









extern "C" void gpu_gauge_derivative(int withrectangles, hamiltonian_field_t * const hf, double c_gauge, double c_rect){
cudaError_t cudaerr;
printf("GPU gauge derivative..\n");
int host_check_VOL2;

#ifdef MPI
 if(g_proc_id==0) printf("Exchanging gauge...");
 if(g_proc_id==0) printf("done\n");
#endif
  
  if((cudaerr=cudaMemcpyToSymbol(dev_VOL2, &(VOLUME), sizeof(int)))!=cudaSuccess){
    printf("gpu_gauge_derivative(): Could not copy dev_VOL2 to device. Aborting...\n");
    exit(200);
  } 
  cudaMemcpyFromSymbol(&host_check_VOL2, dev_VOL2, sizeof(int)); 
  printf("\tOn device:\n");
  printf("\tdev_VOL2 = %i\n", host_check_VOL2);  

  update_gpu_fields(hf->gaugefield, hf->derivative,1);
  if ((cudaerr=cudaPeekAtLastError())!=cudaSuccess) {
        printf("%s\n", cudaGetErrorString(cudaPeekAtLastError()));
  }
//printf("griddim (gauge_monomial): %d, blockdim: %d\n", griddimgauge, blockdimgauge);

 

//plaquette
cudaFuncSetCacheConfig(dev_gauge_derivative, cudaFuncCachePreferL1);
dev_gauge_derivative <<<griddimgauge, blockdimgauge >>> 
            (dev_nn2, dev_gf_d, dev_df0_d, c_gauge);   
if ((cudaerr=cudaGetLastError())!=cudaSuccess) {
   printf("%s\n", cudaGetErrorString(cudaGetLastError()));
   printf("Error code is: %f\n",cudaerr);
}           


  
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
            (dev_nn2, dev_gf_d, dev_df0_d, c_rect, mu);  
    if ((cudaerr=cudaGetLastError())!=cudaSuccess) {
       printf("%s\n", cudaGetErrorString(cudaGetLastError()));
    }
  //here we do
  //    __                         __
  //   |  |                      ||  | 
  //  ||__|            and        |__|
  //                          
    dev_rect_deri_two <<<griddimgauge, blockdimgauge >>> 
              (dev_nn2, dev_gf_d, dev_df0_d, c_rect, mu);  
    if ((cudaerr=cudaGetLastError())!=cudaSuccess) {
       printf("%s\n", cudaGetErrorString(cudaGetLastError()));
    }
  //here we do
  //   __                         __    
  //  |  |            and        |  ||    
  //  |__||                      |__|   
  //     
    dev_rect_deri_three <<<griddimgauge, blockdimgauge >>> 
              (dev_nn2, dev_gf_d, dev_df0_d, c_rect, mu);  
    if ((cudaerr=cudaGetLastError())!=cudaSuccess) {
       printf("%s\n", cudaGetErrorString(cudaGetLastError()));
    }
  }
}   

  cudaMemcpyFromSymbol(&host_check_VOL2, dev_VOL2, sizeof(int)); 
  printf("\tOn device:\n");
  printf("\tdev_VOL2 = %i\n", host_check_VOL2); 
  to_host_mom(hf);

  if ((cudaerr=cudaPeekAtLastError())!=cudaSuccess) {
        printf("%s\n", cudaGetErrorString(cudaPeekAtLastError()));
        printf("Error code is: %f\n",cudaerr);
  }
printf("finished: GPU gauge derivative..\n");

}





extern "C" void benchmark_gpu_gauge_derivative(int withrectangles, hamiltonian_field_t * const hf,double c){
cudaError_t cudaerr;

  clock_t start, stop; 
  double timeelapsed = 0.0;


update_gpu_fields(hf->gaugefield, hf->derivative,1);
//printf("griddim (gauge_monomial): %d, blockdim: %d\n", griddimgauge, blockdimgauge);

int i;
printf("Applying gpu_gauge_derivative 100 times...");
assert((start = clock())!=-1);
for(i=0;i<100;i++){
//plaquette
cudaFuncSetCacheConfig(dev_gauge_derivative, cudaFuncCachePreferL1);
dev_gauge_derivative <<<griddimgauge, blockdimgauge >>> 
            (dev_nn2, dev_gf_d, dev_df0_d, 1.0);   
if ((cudaerr=cudaGetLastError())!=cudaSuccess) {
   printf("%s\n", cudaGetErrorString(cudaGetLastError()));
}           

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
            (dev_nn2, dev_gf_d, dev_df0_d, c, mu);  
    if ((cudaerr=cudaGetLastError())!=cudaSuccess) {
       printf("%s\n", cudaGetErrorString(cudaGetLastError()));
    }
  //here we do
  //    __                         __
  //   |  |                      ||  | 
  //  ||__|            and        |__|
  //        
    dev_rect_deri_two <<<griddimgauge, blockdimgauge >>> 
              (dev_nn2, dev_gf_d, dev_df0_d, c, mu);  
    if ((cudaerr=cudaGetLastError())!=cudaSuccess) {
       printf("%s\n", cudaGetErrorString(cudaGetLastError()));
    }
  //here we do
  //   __                         __    
  //  |  |            and        |  ||    
  //  |__||                      |__|   
  //     
    dev_rect_deri_three <<<griddimgauge, blockdimgauge >>> 
              (dev_nn2, dev_gf_d, dev_df0_d, c, mu);  
    if ((cudaerr=cudaGetLastError())!=cudaSuccess) {
       printf("%s\n", cudaGetErrorString(cudaGetLastError()));
    }
  }
}  

}//100 x

printf("Error code of last error check is: %f\n",cudaerr);
assert((stop = clock())!=-1);
printf(" done\n");
timeelapsed = (double) (stop-start)/CLOCKS_PER_SEC;
printf("Time elapsed: %e sec\n",timeelapsed);

to_host_mom(hf);
exit(0);
}





/* include gauge_update on device here */
#include "hybrid_update.cuh"
/* include deriv_Sb on device here */
#include "deriv_Sb.cuh"
/* include the wilson flow on device*/
//#include "wilson_flow.cuh"
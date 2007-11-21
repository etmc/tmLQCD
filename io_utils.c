/* $Id$ */

#include"lime.h" 
#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<time.h>
#include<sys/time.h> 
#include<sys/types.h>
#ifdef MPI
# include<unistd.h> 
#endif
#include<math.h>
#include"global.h"
#include"su3.h"
#include"lime.h" 
#include"io_utils.h"

int isnan_f  (float       x) { return x != x; }
int isnan_d  (double      x) { return x != x; }
int isnan_ld (long double x) { return x != x; }


int write_xlf_info(const double plaq, const int counter, char * filename) {
  FILE * ofs;
  LimeWriter * limewriter = NULL;
  LimeRecordHeader * limeheader = NULL;
  /* Message end and Message begin flag */
  int ME_flag=0, MB_flag=0, status=0;
#ifdef MPI
  MPI_Status mpi_status;
#endif
  char message[500];
  n_uint64_t bytes;
  struct timeval t1;
  

  gettimeofday(&t1,NULL);
  if(g_kappa > 0. || g_kappa < 0.) {
    sprintf(message,"\n plaquette = %e\n trajectory nr = %d\n beta = %f, kappa = %f, mu = %f, c2_rec = %f\n time = %ld\n hmcversion = %s\n mubar = %f\n epsilonbar = %f\n ", 
	    plaq, counter, g_beta, g_kappa, g_mu/2./g_kappa, g_rgi_C1,t1.tv_sec, PACKAGE_VERSION, 
	    g_mubar/2./g_kappa, g_epsbar/2./g_kappa);
  }
  else {
    sprintf(message,"\n plaquette = %e\n trajectory nr = %d\n beta = %f, kappa = %f, 2*kappa*mu = %f, c2_rec = %f", 
	    plaq, counter, g_beta, g_kappa, g_mu, g_rgi_C1);
  }
  bytes = strlen( message );
  if(g_cart_id == 0) {
    ofs = fopen(filename, "a");
    if(ofs == (FILE*)NULL) {
      fprintf(stderr, "Could not open file %s for writing!\n Aborting...\n", filename);
#ifdef MPI
      MPI_Abort(MPI_COMM_WORLD, 1);
      MPI_Finalize();
#endif
      exit(500);
    }
    limewriter = limeCreateWriter( ofs );
    if(limewriter == (LimeWriter*)NULL) {
      fprintf(stderr, "LIME error in file %s for writing!\n Aborting...\n", filename);
#ifdef MPI
      MPI_Abort(MPI_COMM_WORLD, 1);
      MPI_Finalize();
#endif
      exit(500);
    }
    
    limeheader = limeCreateHeader(MB_flag, ME_flag, "xlf-info", bytes);
    status = limeWriteRecordHeader( limeheader, limewriter);
    if(status < 0 ) {
      fprintf(stderr, "LIME write header error %d\n", status);
#ifdef MPI
      MPI_Abort(MPI_COMM_WORLD, 1);
      MPI_Finalize();
#endif
      exit(500);
    }
    limeDestroyHeader( limeheader );
    limeWriteRecordData(message, &bytes, limewriter);
    
    limeDestroyWriter( limewriter );
    fclose(ofs);
    fflush(ofs);
    
  }
  return(0);
}



int big_endian(){
  union{
    int l;
    char c[sizeof(int)];
  } u;

  u.l=1;
  return(u.c[sizeof(int) - 1] == 1);
}

void write_su3(su3 * up, FILE * f) {
  fprintf(f,"%f %f %f %f %f %f \n%f %f %f %f %f %f \n%f %f %f %f %f %f %d\n\n",
	     (*up).c00.re, (*up).c00.im, (*up).c01.re, (*up).c01.im,
	     (*up).c02.re, (*up).c02.im, (*up).c10.re, (*up).c10.im,
	     (*up).c11.re, (*up).c11.im, (*up).c12.re, (*up).c12.im,
	     (*up).c20.re, (*up).c20.im, (*up).c21.re, (*up).c21.im,
	     (*up).c22.re, (*up).c22.im, g_cart_id);
}

void byte_swap(void * ptr, int nmemb){
  int j;
  char char_in[4];
  char * in_ptr;
  int * int_ptr;

  for(j = 0, int_ptr = (int *) ptr; j < nmemb; j++, int_ptr++){
    in_ptr = (char *) int_ptr;
    
    char_in[0] = in_ptr[0];
    char_in[1] = in_ptr[1];
    char_in[2] = in_ptr[2];
    char_in[3] = in_ptr[3];

    in_ptr[0] = char_in[3];
    in_ptr[1] = char_in[2];
    in_ptr[2] = char_in[1];
    in_ptr[3] = char_in[0];
  }
}

void byte_swap_assign(void * out_ptr, void * in_ptr, int nmemb){
  int j;
  char * char_in_ptr, * char_out_ptr;
  double * double_in_ptr, * double_out_ptr;

  double_in_ptr = (double *) in_ptr;
  double_out_ptr = (double *) out_ptr;
  for(j = 0; j < nmemb; j++){
    char_in_ptr = (char *) double_in_ptr;
    char_out_ptr = (char *) double_out_ptr;
    
    char_out_ptr[7] = char_in_ptr[0];
    char_out_ptr[6] = char_in_ptr[1];
    char_out_ptr[5] = char_in_ptr[2];
    char_out_ptr[4] = char_in_ptr[3];
    char_out_ptr[3] = char_in_ptr[4];
    char_out_ptr[2] = char_in_ptr[5];
    char_out_ptr[1] = char_in_ptr[6];
    char_out_ptr[0] = char_in_ptr[7];
    double_in_ptr++;
    double_out_ptr++;
  }
}

void single2double(void * out_ptr, void * in_ptr, int nmemb) {
  int i;
  float * float_ptr = (float*) in_ptr;
  double * double_ptr = (double*) out_ptr;

  for(i = 0; i < nmemb; i++) {
    (*double_ptr) = (double) (*float_ptr);

    float_ptr++;
    double_ptr++;
  }

}

void double2single(void * out_ptr, void * in_ptr, int nmemb) {
  int i;
  float * float_ptr = (float*) out_ptr;
  double * double_ptr = (double*) in_ptr;

  for(i = 0; i < nmemb; i++) {
    (*float_ptr) = (float) (*double_ptr);

    float_ptr++;
    double_ptr++;
  }

}

void byte_swap_assign_single2double(void * out_ptr, void * in_ptr, int nmemb){
  int j;
  char * char_in_ptr, * char_out_ptr;
  double * double_out_ptr;
  float * float_in_ptr;
  float tmp;

  float_in_ptr = (float *) in_ptr;
  double_out_ptr = (double *) out_ptr;
  char_out_ptr = (char *) &tmp;
  for(j = 0; j < nmemb; j++){
    char_in_ptr = (char *) float_in_ptr;
    
    char_out_ptr[3] = char_in_ptr[0];
    char_out_ptr[2] = char_in_ptr[1];
    char_out_ptr[1] = char_in_ptr[2];
    char_out_ptr[0] = char_in_ptr[3];
    (*double_out_ptr) = (double) tmp;
    float_in_ptr++;
    double_out_ptr++;
  }
}

void byte_swap_assign_double2single(void * out_ptr, void * in_ptr, int nmemb){
  int j;
  char * char_in_ptr, * char_out_ptr;
  double * double_in_ptr;
  float * float_out_ptr;
  float tmp;

  float_out_ptr = (float *) out_ptr;
  double_in_ptr = (double *) in_ptr;
  char_in_ptr = (char *) &tmp;
  for(j = 0; j < nmemb; j++){
    tmp = (float) (*double_in_ptr);
    char_out_ptr = (char*) float_out_ptr;

    char_out_ptr[3] = char_in_ptr[0];
    char_out_ptr[2] = char_in_ptr[1];
    char_out_ptr[1] = char_in_ptr[2];
    char_out_ptr[0] = char_in_ptr[3];

    float_out_ptr++;
    double_in_ptr++;
  }
}


void single2double_cm(spinor * const R, float * const S) {
  (*R).s0.c0.re = (double) S[0];
  (*R).s0.c0.im = (double) S[1];
  (*R).s0.c1.re = (double) S[2];
  (*R).s0.c1.im = (double) S[3];
  (*R).s0.c2.re = (double) S[4];
  (*R).s0.c2.im = (double) S[5];
  (*R).s1.c0.re = (double) S[6];
  (*R).s1.c0.im = (double) S[7];
  (*R).s1.c1.re = (double) S[8];
  (*R).s1.c1.im = (double) S[9];
  (*R).s1.c2.re = (double) S[10];
  (*R).s1.c2.im = (double) S[11];
  (*R).s2.c0.re = (double) S[12];
  (*R).s2.c0.im = (double) S[13];
  (*R).s2.c1.re = (double) S[14];
  (*R).s2.c1.im = (double) S[15];
  (*R).s2.c2.re = (double) S[16];
  (*R).s2.c2.im = (double) S[17];
  (*R).s3.c0.re = (double) S[18];
  (*R).s3.c0.im = (double) S[19];
  (*R).s3.c1.re = (double) S[20];
  (*R).s3.c1.im = (double) S[21];
  (*R).s3.c2.re = (double) S[22];
  (*R).s3.c2.im = (double) S[23];
}

void double2single_cm(float * const S, spinor * const R) {
  S[0]  = (float) (*R).s0.c0.re ;
  S[1]  = (float) (*R).s0.c0.im ;
  S[2]  = (float) (*R).s0.c1.re ;
  S[3]  = (float) (*R).s0.c1.im ;
  S[4]  = (float) (*R).s0.c2.re ;
  S[5]  = (float) (*R).s0.c2.im ;
  S[6]  = (float) (*R).s1.c0.re ;
  S[7]  = (float) (*R).s1.c0.im ;
  S[8]  = (float) (*R).s1.c1.re ;
  S[9]  = (float) (*R).s1.c1.im ;
  S[10] = (float) (*R).s1.c2.re ;
  S[11] = (float) (*R).s1.c2.im ;
  S[12] = (float) (*R).s2.c0.re ;
  S[13] = (float) (*R).s2.c0.im ;
  S[14] = (float) (*R).s2.c1.re ;
  S[15] = (float) (*R).s2.c1.im ;
  S[16] = (float) (*R).s2.c2.re ;
  S[17] = (float) (*R).s2.c2.im ;
  S[18] = (float) (*R).s3.c0.re ;
  S[19] = (float) (*R).s3.c0.im ;
  S[20] = (float) (*R).s3.c1.re ;
  S[21] = (float) (*R).s3.c1.im ;
  S[22] = (float) (*R).s3.c2.re ;
  S[23] = (float) (*R).s3.c2.im ;
}

void zero_spinor(spinor * const R) {
  (*R).s0.c0.re = 0.;
  (*R).s0.c0.im = 0.;
  (*R).s0.c1.re = 0.;
  (*R).s0.c1.im = 0.;
  (*R).s0.c2.re = 0.;
  (*R).s0.c2.im = 0.;
  (*R).s1.c0.re = 0.;
  (*R).s1.c0.im = 0.;
  (*R).s1.c1.re = 0.;
  (*R).s1.c1.im = 0.;
  (*R).s1.c2.re = 0.;
  (*R).s1.c2.im = 0.;
  (*R).s2.c0.re = 0.;
  (*R).s2.c0.im = 0.;
  (*R).s2.c1.re = 0.;
  (*R).s2.c1.im = 0.;
  (*R).s2.c2.re = 0.;
  (*R).s2.c2.im = 0.;
  (*R).s3.c0.re = 0.;
  (*R).s3.c0.im = 0.;
  (*R).s3.c1.re = 0.;
  (*R).s3.c1.im = 0.;
  (*R).s3.c2.re = 0.;
  (*R).s3.c2.im = 0.;
}

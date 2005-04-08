/* $Id$ */

/****************************************************
 * IO routines:
 *
 * write_gauge_field(char * filename)
 *   writes gauge field configuration to file
 *   with name filename.
 *
 * read_gauge_field(char * filename)
 *   reads gauge field configuration from file
 *   with name filename.
 *
 * int big_endian()
 *   returns 1 if data is in big endian order, 0 otherwise
 *
 * void byte_swap(void *ptr, int nmemb)
 *   Swap bytes in an array of nmemb integers ptr
 *
 * void byte_swap_assign(void * out_ptr, void * in_ptr, int nmemb)
 *   Swap bytes in an array of nmemb doubles in_ptr and assigns it to out_ptr
 *
 * Autor: Ines Wetzorke <wetzorke@ifh.de>
 *        Carsten Urbach <urbach@ifh.de>
 *
 ****************************************************/

#include<stdlib.h>
#include<stdio.h>
#ifdef MPI
#include<unistd.h>
#endif
#include"global.h"
#include"su3.h"
#include"io.h"

void byte_swap(void *ptr, int nmemb);
void byte_swap_assign(void * out_ptr, void * in_ptr, int nmemb);
int big_endian();


#ifdef MPI
int write_gauge_field_time_p(char * filename){
  FILE * ofs = NULL;
  int tag=0, t, x, y, z, id, X, tt;
  MPI_Status status;
  su3 tmp[4];
  int coords[2];
  if(g_cart_id == 0){
    ofs = fopen(filename, "w");
    if(ofs != NULL ){
      fprintf(ofs,"%f %d %d\n", g_beta, g_nproc_x*LX, g_nproc_t*T);
    }
    else{
      fprintf(stderr, "Warning! Could not open file %s in routine write_gauge_field\n", filename);
/*       errorhandler(100, filename); */
      return(1);
    }
  }
  
  for(x = 0; x < LX*g_nproc_x; x++){
    X = x - g_proc_coords[1]*LX; 
    coords[1] = x / LX;
    for(y = 0; y < LY; y++){
      for(z = 0; z < LZ; z++){
	tag = 0;
	for(t = 0; t < T*g_nproc_t; t++){
	  tt = t - g_proc_coords[0]*T; 
	  coords[0] = t / T;
	  MPI_Cart_rank(g_cart_grid, coords, &id);
	  if(g_cart_id == 0){
	    if(g_cart_id == id){
#ifdef LITTLE_ENDIAN
	      byte_swap_assign(tmp, g_gauge_field[ g_ipt[tt][X][y][z] ], 4*sizeof(su3)/8);
	      fwrite(tmp, sizeof(su3), 4, ofs);
#else
	      fwrite(g_gauge_field[ g_ipt[tt][X][y][z] ], 4*sizeof(su3), 1, ofs);
#endif
	    }
	    else{
	      MPI_Recv(tmp, 4*sizeof(su3)/8, MPI_DOUBLE, id, tag, g_cart_grid, &status);
	      fwrite(tmp, sizeof(su3), 4, ofs);
	    }
	  }
	  else{
	    if(g_cart_id == id){
#ifdef LITTLE_ENDIAN
	      byte_swap_assign(tmp, g_gauge_field[ g_ipt[tt][X][y][z] ], 4*sizeof(su3)/8);
	      MPI_Send((void*) tmp, 4*sizeof(su3)/8, MPI_DOUBLE, 0, tag, g_cart_grid);
#else
	      MPI_Send((void*) g_gauge_field[ g_ipt[tt][X][y][z] ], 4*sizeof(su3)/8, MPI_DOUBLE, 0, tag, g_cart_grid);
#endif
	    }
	  }
	  tag++;
	}
	MPI_Barrier(g_cart_grid);
      }
    }
  }
  if(g_cart_id == 0){
    if(ferror(ofs)){
      printf("Error! Failed to write all data to file %s in routine read_gauge_field_time_p\n or could not open file\n",filename);
      return(1);
/*       errorhandler(101, filename); */
    }
    fclose(ofs);
  }
  return(0);
}
#else
int write_gauge_field_time_p(char * filename){
  FILE * ofs;
  int t, x, y, z;
#ifdef LITTLE_ENDIAN
  su3 tmp[4];
#endif
  ofs = fopen(filename, "w");
  if(ofs != NULL ){
    fprintf(ofs,"%f %d %d\n", g_beta, LX, T);
    for(x = 0; x < LX; x++){
      for(y = 0; y < LY; y++){
	for(z = 0; z < LZ; z++){
	  for(t = 0; t < T; t++){
#ifdef LITTLE_ENDIAN
 	    byte_swap_assign(tmp, g_gauge_field[ g_ipt[t][x][y][z] ], 4*sizeof(su3)/8); 
	    fwrite(tmp, sizeof(su3), 4, ofs);
#else
	    fwrite(g_gauge_field[ g_ipt[t][x][y][z] ], sizeof(su3), 4, ofs);
#endif
	  }
	}
      }
    }
    if(ferror(ofs)){
/*       errorhandler(100, filename); */
    }
  }
  else{
/*     errorhandler(100, filename); */
  }
  fclose(ofs);
  return(0);
}
#endif


int read_gauge_field_time_p(char * filename){
  FILE * ifs;
  double beta;
  int l, t, x, y, z;
#ifdef MPI
  int position;
#endif
#ifdef LITTLE_ENDIAN
  su3 tmp[4];
#endif

  ifs = fopen(filename, "r");
  if(ifs != NULL ){
    fscanf(ifs,"%lf %d %d\n",&beta,&l,&t);
#ifdef MPI
    position = ftell(ifs);
#endif
    if(beta!=g_beta){
      fprintf(stderr, "Warning! Configuration %s was produced with a different beta!\n", filename);
/*       errorhandler(112,filename); */
    }
    if((l!=g_nproc_x*LX)||(t!=g_nproc_t*T)){
      printf("Error! Configuration %s was produced with a different lattice size\n Aborting...\n", filename);
      exit(1);
/*       errorhandler(114,filename); */
    }
    /* We do not need to seek here      */
    /* because we never have PARALLELXT */
    /* without PARALLELT                */
    for(x = 0; x < LX; x++){
      for(y = 0; y < LY; y++){
	for(z = 0; z < LZ; z++){
#if (defined MPI && (defined PARALLELT || defined PARALLELXT))
	  fseek(ifs, position +
		(g_proc_coords[0]*T+
		 (((g_proc_coords[1]*LX+x)*LY+y)*LZ+z)*T*g_nproc_t)*4*sizeof(su3),
		SEEK_SET);
#endif
	  for(t = 0; t < T; t++){
#ifdef LITTLE_ENDIAN
	    fread(tmp, sizeof(su3), 4, ifs);
	    byte_swap_assign(g_gauge_field[ g_ipt[t][x][y][z] ], tmp, 4*sizeof(su3)/8);
#else
	    fread(g_gauge_field[ g_ipt[t][x][y][z] ], 4*sizeof(su3), 1, ifs);
#endif
	  }
	}
      }
    }
    if((feof(ifs)) || (ferror(ifs))){
      printf("Error! Could not read all data from %s or could not open file!\n Aborting!\n", filename);
      exit(1);
/*       errorhandler(101, filename);  */
    }
  }
  else{
    printf("Error! Could not read all data from %s or could not open file!\n Aborting!\n", filename);
    exit(1);
/*     errorhandler(101, filename);  */
  }
  fclose(ifs);
  return(0);
}

int write_spinorfield_eo_time_p(spinor * const s, spinor * const r, char * filename, const int append){
  FILE * ofs = NULL;
  int x, X, y, z, t, t0, tag=0, id=0, i=0;
  spinor tmp[1];
  int coords[2];
  spinor * p = NULL;
#ifdef MPI
  MPI_Status status;
#endif

  if(g_cart_id == 0){
    if(append == 0) {
      ofs = fopen(filename, "w");
      if(ofs != NULL ){
	fprintf(ofs,"%f %f %f %d %d\n",g_beta, g_kappa, g_mu, LX*g_nproc_x, T*g_nproc_t);
      }
      else{
	/*       errorhandler(106, filename); */
      }
    }
    else {
      ofs = fopen(filename, "a");
      if(ofs == NULL ) {
	printf("Could not open file %s!\n", filename);
      }
    }
  }
  for(x = 0; x < LX*g_nproc_x; x++){
    X = x - g_proc_coords[1]*LX;
    coords[1] = x / LX;
    for(y = 0; y < LY; y++){
      for(z = 0; z < LZ; z++){
	for(t0 = 0; t0 < T*g_nproc_t; t0++){
	  t = t0 - T*g_proc_coords[0];
	  coords[0] = t0 / T;
#ifdef MPI
	  MPI_Cart_rank(g_cart_grid, coords, &id);
#endif
	  if(g_cart_id == id) {
	    i = g_lexic2eosub[ g_ipt[t][X][y][z] ];
	    if((t+X+y+z+g_proc_coords[0]*T+g_proc_coords[1]*LX)%2 == 0) {
	      p = s;
	    }
	    else {
	      p = r;
	    }
	  }
	  if(g_cart_id == 0){
	    if(g_cart_id == id){
#ifdef LITTLE_ENDIAN
	      byte_swap_assign(tmp, p + i , sizeof(spinor)/8);
	      fwrite(tmp, sizeof(spinor), 1, ofs);
#else
	      fwrite(p + i, sizeof(spinor), 1, ofs);
#endif
	    }
#ifdef MPI
	    else{
	      MPI_Recv(tmp, sizeof(spinor)/8, MPI_DOUBLE, id, tag, g_cart_grid, &status);
	      fwrite(tmp, sizeof(spinor), 1, ofs);
	    }
#endif
	  }
#ifdef MPI
	  else{
	    if(g_cart_id == id){
#ifdef LITTLE_ENDIAN
	      byte_swap_assign(tmp, p + i, sizeof(spinor)/8);
	      MPI_Send((void*) tmp, sizeof(spinor)/8, MPI_DOUBLE, 0, tag, g_cart_grid);
#else
	      MPI_Send((void*) (p + i), sizeof(spinor)/8, MPI_DOUBLE, 0, tag, g_cart_grid);
#endif		  
	    }
	  }
#endif
	  tag++;
	}
#ifdef MPI
 	MPI_Barrier(g_cart_grid); 
#endif
	tag=0;
      }
    }
  }
  if(g_cart_id == 0){
    if(ferror(ofs)){
/*       errorhandler(106, filename); */
    }
    fclose(ofs);
  }
  return(0);
}


int read_spinorfield_eo_time(spinor * const s, spinor * const r, char * filename){
  FILE * ifs;
  double beta,kappa,mu;
  int l, t, x, y , z, i = 0;
  spinor * p = NULL;
#ifdef MPI
  int position;
#endif
#ifdef LITTLE_ENDIAN
  spinor tmp[1];
#endif

  ifs = fopen(filename, "r");
  if(ifs != NULL ){
    fscanf(ifs,"%lf %lf %lf %d %d\n", &beta, &kappa, &mu, &l, &t);
#ifdef MPI
    position = ftell(ifs);
#endif
    if((beta!=g_beta)||(g_kappa!=kappa)||(g_mu!=mu)){
      fprintf(stderr, "Warning! Parameters beta, kappa or mu are inconsistent with file %s!\n", filename);
/*       errorhandler(113,filename); */
    }
    if((l!=LX*g_nproc_x)||(t!=T*g_nproc_t)){
      printf("Error! spinorfield %s was produced for a different lattice size!\nAborting!\n", filename);
      exit(1);
/*       errorhandler(115,filename); */
    }
    /* We do not need a seek here, see write_gauge_field... */
    for(x = 0; x < LX; x++){
#if (defined MPI && defined PARALLEL2)
      fseek(ifs, position +
	    ((g_proc_coords[0]*LX+x)*LY*LZ*T 
	     + g_proc_coords[1]*LY*LZ*T)*sizeof(spinor),
	    SEEK_SET);
#endif
      for(y = 0; y < LY; y++){
#if (defined MPI && defined PARALLEL3)
	fseek(ifs, position +
	      ((g_proc_coords[0]*LX+x)*LY*LZ*T 
	       + (g_proc_coords[1]*LY+y)*LZ*T 
	       + g_proc_coords[2]*LZ*T)*sizeof(spinor),
	      SEEK_SET);
#endif
	for(z = 0; z < LZ; z++){
#if (defined MPI && (defined PARALLELT || defined PARALLELXT))
	  fseek(ifs, position +
		(g_proc_coords[0]*T+
		 (((g_proc_coords[1]*LX+x)*LY+y)*LZ+z)*T*g_nproc_t)*sizeof(spinor),
		SEEK_SET);
#endif
	  for(t = 0; t < T; t++){
	    i = g_lexic2eosub[ g_ipt[t][x][y][z] ];
	    if((t+x+y+z+g_proc_coords[0]*T+g_proc_coords[1]*LX)%2==0) {
	      p = s;
	    }
	    else {
	      p = r;
	    }
#ifdef LITTLE_ENDIAN
	    fread(tmp, sizeof(spinor), 1, ifs);
	    byte_swap_assign(p + i, tmp, sizeof(spinor)/8);
#else
	    fread(p + i, sizeof(spinor), 1, ifs);
#endif
	  }
	}
      }
    }
    if((feof(ifs)) || (ferror(ifs))){
/*       errorhandler(107, filename); */
      printf("Error while reading from filed %s!\nAborting!\n", filename);
      exit(1);
    }
  }
  else{
/*     errorhandler(107, filename); */
    printf("Error while reading from filed %s!\nAborting!\n", filename);
    exit(1);
  }
  fclose(ifs);
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

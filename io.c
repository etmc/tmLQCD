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
  int tag=0, t, x, y, z, id;
  int coords[3];
  MPI_Status status;
  su3 tmp[4];
  if(g_proc_id == 0){
    ofs = fopen(filename, "w");
    if(ofs != NULL ){
      fprintf(ofs,"%f %d %d\n", g_beta, L, g_nproc*T);
    }
    else{
      printf("Warning! Could not open file %s in routine write_gauge_field\n", filename);
/*       errorhandler(100, filename); */
      return(1);
    }
  }
  
  for(x = 0; x < LX; x++){
    for(y = 0; y < LY; y++){
      for(z = 0; z < LZ; z++){
	tag = 0;
	for(id = 0; id < g_nproc; id++){
	  for(t = 0; t < T; t++){
	    if(g_proc_id == 0){
	      if(g_proc_id == id){
#ifdef LITTLE_ENDIAN
		byte_swap_assign(tmp, g_gauge_field[ g_ipt[t][x][y][z] ], 4*sizeof(su3)/8);
		fwrite(tmp, sizeof(su3), 4, ofs);
#else
		fwrite(g_gauge_field[ g_ipt[t][x][y][z] ], 4*sizeof(su3), 1, ofs);
#endif
	      }
	      else{
		MPI_Recv(tmp, 4*sizeof(su3)/8, MPI_DOUBLE, id, tag, g_cart_grid, &status);
		fwrite(tmp, sizeof(su3), 4, ofs);
	      }
	    }
	    else{
	      if(g_proc_id == id){
#ifdef LITTLE_ENDIAN
		byte_swap_assign(tmp, g_gauge_field[ g_ipt[t][x][y][z] ], 4*sizeof(su3)/8);
		MPI_Ssend((void*) tmp, 4*sizeof(su3)/8, MPI_DOUBLE, 0, tag, g_cart_grid);
#else
		MPI_Ssend((void*) g_gauge_field[ g_ipt[t][x][y][z] ], 4*sizeof(su3)/8, MPI_DOUBLE, 0, tag, g_cart_grid);
#endif
	      }
	    }
	    tag++;
	  }
	}
	MPI_Barrier(g_cart_grid);
      }
    }
  }
  if(g_proc_id == 0){
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
    fprintf(ofs,"%f %d %d\n", g_beta, L, T);
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
      errorhandler(100, filename);
    }
  }
  else{
    errorhandler(100, filename);
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
      printf("Warning! Configuration %s was produced with a different beta!\n", filename);
/*       errorhandler(112,filename); */
    }
    if((l!=L)||(t!=T)){
      printf("Error! Configuration %s was produced with a differen lattice size\n Aborting...\n", filename);
      exit(1);
/*       errorhandler(114,filename); */
    }
    for(x = 0; x < LX; x++){
      for(y = 0; y < LY; y++){
	for(z = 0; z < LZ; z++){
#if (defined MPI && defined PARALLELT)
	    fseek(ifs, position +
		  (g_proc_id*T+x*y*z*T*g_nproc)*4*sizeof(su3),
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


int big_endian(){
  union{
    int l;
    char c[sizeof(int)];
  } u;

  u.l=1;
  return(u.c[sizeof(int) - 1] == 1);
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

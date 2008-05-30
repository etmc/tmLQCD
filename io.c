/* $Id$ */

/****************************************************
 * IO routines:
 *
 * write_lime_gauge_field_old(char * filename, const double plaq, const int counter)
 *   write gauge field in ILDG format
 *
 * read_lime_gauge_field_old(char * filename)
 *   read gauge field in ILDG format
 *
 * write_spinorfield_eo_time_p(spinor * const s, spinor * const r, char * filename, const int append)
 *   write a full spinor from one even and one odd spinorfield
 *
 * read_spinorfield_eo_time(spinor * const s, spinor * const r, char * filename)
 *   read a full spinor into one even and one odd spinorfields
 *
 * write_gauge_field(char * filename)
 *   writes gauge field configuration to file
 *   with name filename.
 *
 * read_gauge_field(char * filename)
 *   reads gauge field configuration from file
 *   with name filename.
 *
 * int read_eospinor(spinor * const s, char * filename)
 *   Read an even or odd spinor from file (eigenvector)
 *
 * int write_eospinor(spinor * const s, char * filename, const double evalue, const double prec, const int nstore)
 *   Read an even or odd spinor from file (eigenvector)
 *
 * Autor: Ines Wetzorke <wetzorke@ifh.de>
 *        Carsten Urbach <urbach@ifh.de>
 *
 * Addition of isnan in function read_spinorfield_cm_single to detect a
 * wrong read order and creation of read_spinorfield_cm_swap_single to
 * swap while reading if necessary.
 *
 * Author: Remi Baron <remi.baron@cea.fr> Jul 2007
 *
 ****************************************************/

/*
 * Note:
 * Required version of lime: >= 1.2.3
 * n_uint64_t is a lime defined type!!
 *
 */

#define _FILE_OFFSET_BITS 64

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
#include"read_input.h"
#include"io_utils.h"
#include"io.h"

#define MAXBUF 1048576

n_uint64_t file_size(FILE *fp);

#ifdef MPI
int write_lime_gauge_field_old(char * filename, const double plaq, const int counter){
  FILE * ofs;
  LimeWriter * limewriter = NULL;
  LimeRecordHeader * limeheader = NULL;
  /* Message end and Message begin flag */
  int ME_flag=0, MB_flag=0, status=0;
  int tag=0, t=0, x=0, y=0, z=0, id=0, X=0, tt=0, Y=0, Z=0;
  MPI_Status mpi_status;
  char message[500];
  su3 tmp[4];
  int coords[4];
  n_uint64_t bytes;

  write_xlf_info(plaq, counter, filename, 0);

  if(g_cart_id == 0) {
    ofs = fopen(filename, "a");
    if(ofs == (FILE*)NULL) {
      fprintf(stderr, "Could not open file %s for writing!\n Aboring...\n", filename);
      MPI_Abort(MPI_COMM_WORLD, 1);
      MPI_Finalize();
      exit(500);
    }
    limewriter = limeCreateWriter( ofs );
    if(limewriter == (LimeWriter*)NULL) {
      fprintf(stderr, "LIME error in file %s for writing!\n Aboring...\n", filename);
      MPI_Abort(MPI_COMM_WORLD, 1);
      MPI_Finalize();
      exit(500);
    }
    write_ildg_format_xml("temp.xml", limewriter, 0);
    
    bytes = LX*g_nproc_x*LY*g_nproc_y*LZ*g_nproc_z*T*g_nproc_t*4*sizeof(su3);
    MB_flag=1; ME_flag=1;
    limeheader = limeCreateHeader(MB_flag, ME_flag, "ildg-binary-data", bytes);
    status = limeWriteRecordHeader( limeheader, limewriter);
    if(status < 0 ) {
      fprintf(stderr, "LIME write header error %d\n", status);
      MPI_Abort(MPI_COMM_WORLD, 1);
      MPI_Finalize();
      exit(500);
    }
    limeDestroyHeader( limeheader );
  }

  bytes = sizeof(su3);
  for(t = 0; t < T*g_nproc_t; t++) {
    tt = t - g_proc_coords[0]*T;
    coords[0] = t / T;
    for(z = 0; z < LZ*g_nproc_z; z++) {
      Z = z - g_proc_coords[3]*LZ;
      coords[3] = z / LZ;
      for(y = 0; y < LY*g_nproc_y; y++) {
	tag = 0;
	Y = y - g_proc_coords[2]*LY;
	coords[2] = y / LY;
	for(x = 0; x < LX*g_nproc_x; x++) {
	  X = x - g_proc_coords[1]*LX; 
	  coords[1] = x / LX;

	  MPI_Cart_rank(g_cart_grid, coords, &id);
	  if(g_cart_id == 0) {
	    if(g_cart_id == id) {
#ifndef WORDS_BIGENDIAN
	      byte_swap_assign(tmp, g_gauge_field[ g_ipt[tt][X][Y][Z] ], 4*sizeof(su3)/8); 
	      status = limeWriteRecordData((void*)&tmp[1], &bytes, limewriter);
	      status = limeWriteRecordData((void*)&tmp[2], &bytes, limewriter);
	      status = limeWriteRecordData((void*)&tmp[3], &bytes, limewriter);
	      status = limeWriteRecordData((void*)&tmp[0], &bytes, limewriter);
#else
	      status = limeWriteRecordData((void*)&g_gauge_field[ g_ipt[tt][X][Y][Z] ][1], &bytes, limewriter);
	      status = limeWriteRecordData((void*)&g_gauge_field[ g_ipt[tt][X][Y][Z] ][2], &bytes, limewriter);
	      status = limeWriteRecordData((void*)&g_gauge_field[ g_ipt[tt][X][Y][Z] ][3], &bytes, limewriter);
	      status = limeWriteRecordData((void*)&g_gauge_field[ g_ipt[tt][X][Y][Z] ][0], &bytes, limewriter);
#endif
	    }
	    else {
	      MPI_Recv(tmp, 4*sizeof(su3)/8, MPI_DOUBLE, id, tag, g_cart_grid, &mpi_status);
	      status = limeWriteRecordData((void*)&tmp[1], &bytes, limewriter);
	      status = limeWriteRecordData((void*)&tmp[2], &bytes, limewriter);
	      status = limeWriteRecordData((void*)&tmp[3], &bytes, limewriter);
	      status = limeWriteRecordData((void*)&tmp[0], &bytes, limewriter);
	    }
	    if(status < 0 ) {
	      fprintf(stderr, "LIME write error %d\n", status);
	      fprintf(stderr, "x %d, y %d, z %d, t %d (%d,%d,%d,%d)\n",x,y,z,t,X,Y,Z,tt);
	      MPI_Abort(MPI_COMM_WORLD, 1);
	      MPI_Finalize();
	      exit(500);
	    }
	  }
	  else {
	    if(g_cart_id == id){
#ifndef WORDS_BIGENDIAN
	      byte_swap_assign(tmp, g_gauge_field[ g_ipt[tt][X][Y][Z] ], 4*sizeof(su3)/8);
	      MPI_Send((void*) tmp, 4*sizeof(su3)/8, MPI_DOUBLE, 0, tag, g_cart_grid);
#else
	      MPI_Send((void*) g_gauge_field[ g_ipt[tt][X][Y][Z] ], 4*sizeof(su3)/8, MPI_DOUBLE, 0, tag, g_cart_grid);
#endif
	    }
	  }
	  tag++;
	}
	MPI_Barrier(g_cart_grid);
      }
    }
  }
  if(g_cart_id == 0) {
    limeDestroyWriter( limewriter );
    fflush(ofs);
    fclose(ofs);
  }
  return(0);
}
#else
int write_lime_gauge_field_old(char * filename, const double plaq, const int counter){
  FILE * ofs;
  LimeWriter * limewriter;
  LimeRecordHeader * limeheader;
  /* Message end and Message begin flag */
  int ME_flag=0, MB_flag=0, status=0;
  int t=0, x=0, y=0, z=0;
  char message[100];
  n_uint64_t bytes;
#ifndef WORDS_BIGENDIAN
  su3 tmp[4];
#endif

  write_xlf_info(plaq, counter, filename, 0);

  ofs = fopen(filename, "a");
  if(ofs == (FILE*)NULL) {
    fprintf(stderr, "Could not open file %s for writing!\n Aboring...\n", filename);
    exit(500);
  }
  limewriter = limeCreateWriter( ofs );
  if(limewriter == (LimeWriter*)NULL) {
    fprintf(stderr, "LIME error in file %s for writing!\n Aboring...\n", filename);
    exit(500);
  }
  write_ildg_format_xml("temp.xml", limewriter, 0);

  bytes = (n_uint64_t)(LX*LY*LZ*T*4*sizeof(su3));
  MB_flag=0; ME_flag=1;
  limeheader = limeCreateHeader(MB_flag, ME_flag, "ildg-binary-data", bytes);
  status = limeWriteRecordHeader( limeheader, limewriter);
  if(status < 0 ) {
    fprintf(stderr, "LIME write header error %d\n", status);
    exit(500);
  }
  limeDestroyHeader( limeheader );

  bytes = sizeof(su3);
  for(t = 0; t < T; t++){
    for(z = 0; z < LZ; z++){
      for(y = 0; y < LY; y++){
	for(x = 0; x < LX; x++){
#ifndef WORDS_BIGENDIAN
	  byte_swap_assign(tmp, g_gauge_field[ g_ipt[t][x][y][z] ], 4*sizeof(su3)/8); 
	  status = limeWriteRecordData((void*)&tmp[1], &bytes, limewriter);
	  status = limeWriteRecordData((void*)&tmp[2], &bytes, limewriter);
	  status = limeWriteRecordData((void*)&tmp[3], &bytes, limewriter);
	  status = limeWriteRecordData((void*)&tmp[0], &bytes, limewriter);
#else
	  status = limeWriteRecordData((void*)&g_gauge_field[ g_ipt[t][x][y][z] ][1], &bytes, limewriter);
	  status = limeWriteRecordData((void*)&g_gauge_field[ g_ipt[t][x][y][z] ][2], &bytes, limewriter);
	  status = limeWriteRecordData((void*)&g_gauge_field[ g_ipt[t][x][y][z] ][3], &bytes, limewriter);
	  status = limeWriteRecordData((void*)&g_gauge_field[ g_ipt[t][x][y][z] ][0], &bytes, limewriter);
#endif
	  if(status < 0 ) {
	    fprintf(stderr, "LIME write error %d\n", status);
	    exit(500);
	  }
	}
      }
    }
  }
  limeDestroyWriter( limewriter );
  fflush(ofs);
  fclose(ofs);
  return(0);
}
#endif

#ifdef MPI
int write_lime_gauge_field_singleprec(char * filename, const double plaq, 
				      const int counter){
  FILE * ofs;
  LimeWriter * limewriter = NULL;
  LimeRecordHeader * limeheader = NULL;
  /* Message end and Message begin flag */
  int ME_flag=0, MB_flag=0, status=0;
  int tag=0, t=0, x=0, y=0, z=0, id=0, X=0, tt=0, Y=0, Z=0;
  MPI_Status mpi_status;
  char message[500];
  float tmp[72];
  int coords[4];
  n_uint64_t bytes;

  write_xlf_info(plaq, counter, filename, 0);

  if(g_cart_id == 0) {
    ofs = fopen(filename, "a");
    if(ofs == (FILE*)NULL) {
      fprintf(stderr, "Could not open file %s for writing!\n Aboring...\n", filename);
      MPI_Abort(MPI_COMM_WORLD, 1);
      MPI_Finalize();
      exit(500);
    }
    limewriter = limeCreateWriter( ofs );
    if(limewriter == (LimeWriter*)NULL) {
      fprintf(stderr, "LIME error in file %s for writing!\n Aboring...\n", filename);
      MPI_Abort(MPI_COMM_WORLD, 1);
      MPI_Finalize();
      exit(500);
    }
    write_ildg_format_xml("temp.xml", limewriter, 1);
    
    bytes = LX*g_nproc_x*LY*g_nproc_y*LZ*g_nproc_z*T*g_nproc_t*4*sizeof(su3)/2;
    MB_flag=0; ME_flag=1;
    limeheader = limeCreateHeader(MB_flag, ME_flag, "ildg-binary-data", bytes);
    status = limeWriteRecordHeader( limeheader, limewriter);
    if(status < 0 ) {
      fprintf(stderr, "LIME write header error %d\n", status);
      MPI_Abort(MPI_COMM_WORLD, 1);
      MPI_Finalize();
      exit(500);
    }
    limeDestroyHeader( limeheader );
  }

  bytes = sizeof(su3)/2;
  for(t = 0; t < T*g_nproc_t; t++) {
    tt = t - g_proc_coords[0]*T;
    coords[0] = t / T;
    for(z = 0; z < LZ*g_nproc_z; z++) {
      Z = z - g_proc_coords[3]*LZ;
      coords[3] = z / LZ;
      for(y = 0; y < LY*g_nproc_y; y++) {
	tag = 0;
	Y = y - g_proc_coords[2]*LY;
	coords[2] = y / LY;
	for(x = 0; x < LX*g_nproc_x; x++) {
	  X = x - g_proc_coords[1]*LX; 
	  coords[1] = x / LX;

	  MPI_Cart_rank(g_cart_grid, coords, &id);
	  if(g_cart_id == 0) {
	    if(g_cart_id == id) {
#ifndef WORDS_BIGENDIAN
	      byte_swap_assign_double2single(tmp, g_gauge_field[ g_ipt[tt][X][Y][Z] ], 
					     4*sizeof(su3)/8); 
#else
	      double2single(tmp, g_gauge_field[ g_ipt[tt][X][Y][Z] ], 
			    4*sizeof(su3)/8); 
#endif
	      status = limeWriteRecordData((void*)&tmp[1*18], &bytes, limewriter);
	      status = limeWriteRecordData((void*)&tmp[2*18], &bytes, limewriter);
	      status = limeWriteRecordData((void*)&tmp[3*18], &bytes, limewriter);
	      status = limeWriteRecordData((void*)&tmp[0*18], &bytes, limewriter);
	    }
	    else {
	      MPI_Recv(tmp, 4*sizeof(su3)/8, MPI_FLOAT, id, tag, g_cart_grid, &mpi_status);
	      status = limeWriteRecordData((void*)&tmp[1*18], &bytes, limewriter);
	      status = limeWriteRecordData((void*)&tmp[2*18], &bytes, limewriter);
	      status = limeWriteRecordData((void*)&tmp[3*18], &bytes, limewriter);
	      status = limeWriteRecordData((void*)&tmp[0*18], &bytes, limewriter);
	    }
	    if(status < 0 ) {
	      fprintf(stderr, "LIME write error %d\n", status);
	      fprintf(stderr, "x %d, y %d, z %d, t %d (%d,%d,%d,%d)\n",x,y,z,t,X,Y,Z,tt);
	      MPI_Abort(MPI_COMM_WORLD, 1);
	      MPI_Finalize();
	      exit(500);
	    }
	  }
	  else {
	    if(g_cart_id == id){
#ifndef WORDS_BIGENDIAN
	      byte_swap_assign_double2single(tmp, g_gauge_field[ g_ipt[tt][X][Y][Z] ],
					     4*sizeof(su3)/8);
#else
	      double2single(tmp, g_gauge_field[ g_ipt[tt][X][Y][Z] ],
			    4*sizeof(su3)/8);
#endif
	      MPI_Send((void*) tmp, 4*sizeof(su3)/8, MPI_FLOAT, 0, tag, g_cart_grid);
	    }
	  }
	  tag++;
	}
	MPI_Barrier(g_cart_grid);
      }
    }
  }
  if(g_cart_id == 0) {
    limeDestroyWriter( limewriter );
    fflush(ofs);
    fclose(ofs);
  }
  return(0);
}
#else
int write_lime_gauge_field_singleprec(char * filename, 
				      const double plaq, const int counter){
  FILE * ofs;
  LimeWriter * limewriter;
  LimeRecordHeader * limeheader;
  /* Message end and Message begin flag */
  int ME_flag=0, MB_flag=0, status=0;
  int t=0, x=0, y=0, z=0;
  char message[500];
  n_uint64_t bytes;
  float tmp[72];

  write_xlf_info(plaq, counter, filename, 0);

  ofs = fopen(filename, "a");
  if(ofs == (FILE*)NULL) {
    fprintf(stderr, "Could not open file %s for writing!\n Aboring...\n", filename);
    exit(500);
  }
  limewriter = limeCreateWriter( ofs );
  if(limewriter == (LimeWriter*)NULL) {
    fprintf(stderr, "LIME error in file %s for writing!\n Aboring...\n", filename);
    exit(500);
  }
  /* Prepare it for 32 Bit write */
  write_ildg_format_xml("temp.xml", limewriter, 1);

  bytes = (n_uint64_t)(LX*LY*LZ*T*4*sizeof(su3))/2;
  MB_flag=0; ME_flag=1;
  limeheader = limeCreateHeader(MB_flag, ME_flag, "ildg-binary-data", bytes);
  status = limeWriteRecordHeader( limeheader, limewriter);
  if(status < 0 ) {
    fprintf(stderr, "LIME write header error %d\n", status);
    exit(500);
  }
  limeDestroyHeader( limeheader );

  bytes = sizeof(su3)/2;
  for(t = 0; t < T; t++) {
    for(z = 0; z < LZ; z++) {
      for(y = 0; y < LY; y++) {
	for(x = 0; x < LX; x++) {
#ifndef WORDS_BIGENDIAN
	  byte_swap_assign_double2single(tmp, g_gauge_field[ g_ipt[t][x][y][z] ], 4*sizeof(su3)/8); 
#else
	  double2single(tmp, g_gauge_field[ g_ipt[t][x][y][z] ], 4*sizeof(su3)/8); 
#endif
	  status = limeWriteRecordData((void*)&tmp[1*18], &bytes, limewriter);
	  status = limeWriteRecordData((void*)&tmp[2*18], &bytes, limewriter);
	  status = limeWriteRecordData((void*)&tmp[3*18], &bytes, limewriter);
	  status = limeWriteRecordData((void*)&tmp[0*18], &bytes, limewriter);
	  if(status < 0 ) {
	    fprintf(stderr, "LIME write error %d\n", status);
	    exit(500);
	  }
	}
      }
    }
  }
  limeDestroyWriter( limewriter );
  fflush(ofs);
  fclose(ofs);
  return(0);
}
#endif

#ifdef MPI
int write_gauge_field_time_p(char * filename){
  FILE * ofs = NULL;
  int tag=0, t, x, y, z, id, X=0, tt=0, Y=0, Z=0;
  MPI_Status status;
  su3 tmp[4];
  int coords[4];
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
    for(y = 0; y < g_nproc_y*LY; y++){
      Y = y - g_proc_coords[2]*LY;
      coords[2] = y / LY;
      for(z = 0; z < LZ*g_nproc_z; z++){
	Z = z - g_proc_coords[3]*LZ;
	coords[3] = z / LZ;
	tag = 0;
	for(t = 0; t < T*g_nproc_t; t++){
	  tt = t - g_proc_coords[0]*T; 
	  coords[0] = t / T;
	  MPI_Cart_rank(g_cart_grid, coords, &id);
	  if(g_cart_id == 0){
	    if(g_cart_id == id){
#ifndef WORDS_BIGENDIAN
	      byte_swap_assign(tmp, g_gauge_field[ g_ipt[tt][X][Y][Z] ], 4*sizeof(su3)/8);
	      fwrite(tmp, sizeof(su3), 4, ofs);
#else
	      fwrite(g_gauge_field[ g_ipt[tt][X][Y][Z] ], 4*sizeof(su3), 1, ofs);
#endif
	    }
	    else{
	      MPI_Recv(tmp, 4*sizeof(su3)/8, MPI_DOUBLE, id, tag, g_cart_grid, &status);
	      fwrite(tmp, sizeof(su3), 4, ofs);
	    }
	  }
	  else{
	    if(g_cart_id == id){
#ifndef WORDS_BIGENDIAN
	      byte_swap_assign(tmp, g_gauge_field[ g_ipt[tt][X][Y][Z] ], 4*sizeof(su3)/8);
	      MPI_Send((void*) tmp, 4*sizeof(su3)/8, MPI_DOUBLE, 0, tag, g_cart_grid);
#else
	      MPI_Send((void*) g_gauge_field[ g_ipt[tt][X][Y][Z] ], 4*sizeof(su3)/8, MPI_DOUBLE, 0, tag, g_cart_grid);
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
#ifndef WORDS_BIGENDIAN
  su3 tmp[4];
#endif
  ofs = fopen(filename, "w");
  if(ofs != NULL ){
    fprintf(ofs,"%f %d %d\n", g_beta, LX, T);
    for(x = 0; x < LX; x++){
      for(y = 0; y < LY; y++){
	for(z = 0; z < LZ; z++){
	  for(t = 0; t < T; t++){
#ifndef WORDS_BIGENDIAN
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
#ifndef WORDS_BIGENDIAN
  su3 tmp[4];
#endif

  ifs = fopen(filename, "r");
  if(ifs != NULL ){
    fscanf(ifs,"%lf %d %d\n",&beta,&l,&t);
#ifdef MPI
    position = ftell(ifs);
#endif
    if(beta!=g_beta){
      if(g_proc_id == 0) {
	fprintf(stderr, "Warning! Configuration %s was produced with a different beta!\n", filename);
      }
    }
    if((l!=g_nproc_x*LX)||(t!=g_nproc_t*T)){
      if(g_proc_id == 0) {
	printf("Error! Configuration %s was produced with a different lattice size\n Aborting...\n", filename);
      }
      exit(1);
/*       errorhandler(114,filename); */
    }
    /* We do not need to seek here      */
    /* because we never have PARALLELXT */
    /* without PARALLELT                */
    for(x = 0; x < LX; x++){
      for(y = 0; y < LY; y++){
	for(z = 0; z < LZ; z++){
#if (defined MPI)
	  fseek(ifs, position +
		(g_proc_coords[0]*T+
		 (((g_proc_coords[1]*LX+x)*g_nproc_y*LY + g_proc_coords[2]*LY + y)*g_nproc_z*LZ
		  + g_proc_coords[3]*LZ+z)*T*g_nproc_t)*4*sizeof(su3),
		SEEK_SET);
#endif
	  for(t = 0; t < T; t++){
#ifndef WORDS_BIGENDIAN
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

int read_lime_gauge_field_old(char * filename){
  FILE * ifs;
  int t, x, y, z, status;
  n_uint64_t bytes;
  char * header_type;
  LimeReader * limereader;
#ifndef WORDS_BIGENDIAN
  su3 tmp[4];
#endif

  ifs = fopen(filename, "r");
  if(ifs == (FILE *)NULL) {
    fprintf(stderr, "Could not open file %s\n Aborting...\n", filename);
#ifdef MPI
    MPI_Abort(MPI_COMM_WORLD, 1);
    MPI_Finalize();
#endif
    exit(500);
  }
  limereader = limeCreateReader( ifs );
  if( limereader == (LimeReader *)NULL ) {
    fprintf(stderr, "Unable to open LimeReader\n");
#ifdef MPI
    MPI_Abort(MPI_COMM_WORLD, 1);
    MPI_Finalize();
#endif
    exit(500);
  }
  while( (status = limeReaderNextRecord(limereader)) != LIME_EOF ) {
    if(status != LIME_SUCCESS ) {
      fprintf(stderr, "limeReaderNextRecord returned error with status = %d!\n", status);
      status = LIME_EOF;
      break;
    }
    header_type = limeReaderType(limereader);
    if(!strcmp("ildg-binary-data",header_type)) break;
  }
  if(status == LIME_EOF) {
    if(g_proc_id == 0) {
      fprintf(stderr, "no ildg-binary-data record found in file %s\n",filename);
      fprintf(stderr, "trying old deprecated file format!\n");
    }
    limeDestroyReader(limereader);
    fclose(ifs);
    read_gauge_field_time_p(filename);
    return(0);
  }
  bytes = limeReaderBytes(limereader);
  if((int)bytes != LX*g_nproc_x*LY*g_nproc_y*LZ*g_nproc_z*T*g_nproc_t*4*sizeof(su3)) {
    if((int)bytes == (LX*g_nproc_x*LY*g_nproc_y*LZ*g_nproc_z*T*g_nproc_t*4*sizeof(su3))/2) {
      return(read_lime_gauge_field_singleprec(filename));
    }
    fprintf(stderr, "Probably wrong lattice size or precision (bytes=%d) in file %s\n", (int)bytes, filename);
    fprintf(stderr, "Aborting...!\n");
    fflush( stdout );
#ifdef MPI
    MPI_Abort(MPI_COMM_WORLD, 1);
    MPI_Finalize();
#endif
    exit(501);
  }

  bytes = (n_uint64_t)4*sizeof(su3);
#ifdef WORDS_BIGENDIAN
  bytes = sizeof(su3);
#endif
  for(t = 0; t < T; t++){
    for(z = 0; z < LZ; z++){
      for(y = 0; y < LY; y++){
#if (defined MPI)
	limeReaderSeek(limereader,(n_uint64_t) 
		       (g_proc_coords[1]*LX + 
			(((g_proc_coords[0]*T+t)*g_nproc_z*LZ+g_proc_coords[3]*LZ+z)*g_nproc_y*LY 
			 + g_proc_coords[2]*LY+y)*LX*g_nproc_x)*4*sizeof(su3),
		       SEEK_SET);
#endif
	for(x = 0; x < LX; x++){
#ifndef WORDS_BIGENDIAN
	  status = limeReaderReadData(tmp, &bytes, limereader);
	  byte_swap_assign(&g_gauge_field[ g_ipt[t][x][y][z] ][0], &tmp[3], sizeof(su3)/8);
	  byte_swap_assign(&g_gauge_field[ g_ipt[t][x][y][z] ][1], &tmp[0], sizeof(su3)/8);
	  byte_swap_assign(&g_gauge_field[ g_ipt[t][x][y][z] ][2], &tmp[1], sizeof(su3)/8);
	  byte_swap_assign(&g_gauge_field[ g_ipt[t][x][y][z] ][3], &tmp[2], sizeof(su3)/8);
#else
	  status = limeReaderReadData(&g_gauge_field[ g_ipt[t][x][y][z] ][1], &bytes, limereader);
	  status = limeReaderReadData(&g_gauge_field[ g_ipt[t][x][y][z] ][2], &bytes, limereader);
	  status = limeReaderReadData(&g_gauge_field[ g_ipt[t][x][y][z] ][3], &bytes, limereader);
	  status = limeReaderReadData(&g_gauge_field[ g_ipt[t][x][y][z] ][0], &bytes, limereader);
#endif
	  if(status < 0 && status != LIME_EOR) {
	    fprintf(stderr, "LIME read error occured with status = %d while reading file %s!\n Aborting...\n", status, filename);
#ifdef MPI
	    MPI_Abort(MPI_COMM_WORLD, 1);
	    MPI_Finalize();
#endif
	    exit(500);
	  }
	}
      }
    }
  }
  limeDestroyReader(limereader);
  fclose(ifs);
  return(0);
}


int read_lime_gauge_field_singleprec(char * filename){
  FILE * ifs;
  int t, x, y, z, status;
  n_uint64_t bytes;
  char * header_type;
  LimeReader * limereader;
  float tmp[72];

  ifs = fopen(filename, "r");
  if(ifs == (FILE *)NULL) {
    fprintf(stderr, "Could not open file %s\n Aborting...\n", filename);
#ifdef MPI
    MPI_Abort(MPI_COMM_WORLD, 1);
    MPI_Finalize();
#endif
    exit(500);
  }
  limereader = limeCreateReader( ifs );
  if( limereader == (LimeReader *)NULL ) {
    fprintf(stderr, "Unable to open LimeReader\n");
#ifdef MPI
    MPI_Abort(MPI_COMM_WORLD, 1);
    MPI_Finalize();
#endif
    exit(500);
  }
  while( (status = limeReaderNextRecord(limereader)) != LIME_EOF ) {
    if(status != LIME_SUCCESS ) {
      fprintf(stderr, "limeReaderNextRecord returned error with status = %d!\n", status);
      status = LIME_EOF;
      break;
    }
    header_type = limeReaderType(limereader);
    if(!strcmp("ildg-binary-data",header_type)) break;
  }
  if(status == LIME_EOF) {
    if(g_proc_id == 0) {
      fprintf(stderr, "no ildg-binary-data record found in file %s\n",filename);
      fprintf(stderr, "trying old deprecated file format!\n");
    }
    limeDestroyReader(limereader);
    fclose(ifs);
    read_gauge_field_time_p(filename);
    return(0);
  }
  bytes = limeReaderBytes(limereader);
  if((int)bytes != LX*g_nproc_x*LY*g_nproc_y*LZ*g_nproc_z*T*g_nproc_t*4*sizeof(su3)/2) {
    if((int)bytes == (LX*g_nproc_x*LY*g_nproc_y*LZ*g_nproc_z*T*g_nproc_t*4*sizeof(su3))) {
      return(read_lime_gauge_field_old(filename));
    }
    fprintf(stderr, "Probably wrong lattice size or precision (bytes=%d) in file %s\n", (int)bytes, filename);
    fprintf(stderr, "Aborting...!\n");
    fflush( stdout );
#ifdef MPI
    MPI_Abort(MPI_COMM_WORLD, 1);
    MPI_Finalize();
#endif
    exit(501);
  }

  bytes = (n_uint64_t)4*sizeof(su3)/2;
  for(t = 0; t < T; t++){
    for(z = 0; z < LZ; z++){
      for(y = 0; y < LY; y++){
#if (defined MPI)
	limeReaderSeek(limereader,(n_uint64_t) 
		       (g_proc_coords[1]*LX + 
			(((g_proc_coords[0]*T+t)*g_nproc_z*LZ+g_proc_coords[3]*LZ+z)*g_nproc_y*LY 
			 + g_proc_coords[2]*LY+y)*LX*g_nproc_x)*4*sizeof(su3)/2,
		       SEEK_SET);
#endif
	for(x = 0; x < LX; x++) {
	  status = limeReaderReadData(tmp, &bytes, limereader);
#ifndef WORDS_BIGENDIAN
	  byte_swap_assign_single2double(&g_gauge_field[ g_ipt[t][x][y][z] ][0], 
					 &tmp[3*18], sizeof(su3)/8);
	  byte_swap_assign_single2double(&g_gauge_field[ g_ipt[t][x][y][z] ][1], 
					 &tmp[0*18], sizeof(su3)/8);
	  byte_swap_assign_single2double(&g_gauge_field[ g_ipt[t][x][y][z] ][2], 
					 &tmp[1*18], sizeof(su3)/8);
	  byte_swap_assign_single2double(&g_gauge_field[ g_ipt[t][x][y][z] ][3], 
					 &tmp[2*18], sizeof(su3)/8);
#else
	  single2double(&g_gauge_field[ g_ipt[t][x][y][z] ][0], 
			&tmp[3*18], sizeof(su3)/8);
	  single2double(&g_gauge_field[ g_ipt[t][x][y][z] ][1], 
			&tmp[0*18], sizeof(su3)/8);
	  single2double(&g_gauge_field[ g_ipt[t][x][y][z] ][2], 
			&tmp[1*18], sizeof(su3)/8);
	  single2double(&g_gauge_field[ g_ipt[t][x][y][z] ][3], 
			&tmp[2*18], sizeof(su3)/8);
#endif
	  if(status < 0 && status != LIME_EOR) {
	    fprintf(stderr, "LIME read error occured with status = %d while reading file %s!\n Aborting...\n", status, filename);
#ifdef MPI
	    MPI_Abort(MPI_COMM_WORLD, 1);
	    MPI_Finalize();
#endif
	    exit(500);
	  }
	}
      }
    }
  }
  limeDestroyReader(limereader);
  fclose(ifs);
  return(0);
}


int write_spinorfield_eo_time_p(spinor * const s, spinor * const r, char * filename, const int append){
  FILE * ofs = NULL;
  int x, X, y, Y, z, Z, t, t0, tag=0, id=0, i=0;
  spinor tmp[1];
  int coords[4];
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
	fprintf(stderr, "Could not open file %s!\n", filename);
	return(-1);
      }
    }
  }
  for(x = 0; x < LX*g_nproc_x; x++) {
    X = x - g_proc_coords[1]*LX;
    coords[1] = x / LX;
    for(y = 0; y < LY*g_nproc_y; y++) {
      Y = y - g_proc_coords[2]*LY;
      coords[2] = y / LY;
      for(z = 0; z < LZ*g_nproc_z; z++) {
      Z = z - g_proc_coords[3]*LZ;
      coords[3] = z / LZ;
	for(t0 = 0; t0 < T*g_nproc_t; t0++) {
	  t = t0 - T*g_proc_coords[0];
	  coords[0] = t0 / T;
#ifdef MPI
	  MPI_Cart_rank(g_cart_grid, coords, &id);
#endif
	  if(g_cart_id == id) {
	    i = g_lexic2eosub[ g_ipt[t][X][Y][Z] ];
	    if((t+X+Y+Z+g_proc_coords[3]*LZ+g_proc_coords[2]*LY 
		+ g_proc_coords[0]*T+g_proc_coords[1]*LX)%2 == 0) {
	      p = s;
	    }
	    else {
	      p = r;
	    }
	  }
	  if(g_cart_id == 0) {
	    if(g_cart_id == id) {
#ifndef WORDS_BIGENDIAN
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
#  ifndef WORDS_BIGENDIAN
	      byte_swap_assign(tmp, p + i, sizeof(spinor)/8);
	      MPI_Send((void*) tmp, sizeof(spinor)/8, MPI_DOUBLE, 0, tag, g_cart_grid);
#  else
	      MPI_Send((void*) (p + i), sizeof(spinor)/8, MPI_DOUBLE, 0, tag, g_cart_grid);
#  endif		  
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
#ifndef WORDS_BIGENDIAN
  spinor tmp[1];
#endif

  ifs = fopen(filename, "r");
  if(ifs != NULL ){
    fscanf(ifs,"%lf %lf %lf %d %d\n", &beta, &kappa, &mu, &l, &t);
#ifdef MPI
    position = ftell(ifs);
#endif
    if(((beta!=g_beta)||(g_kappa!=kappa)||(fabs(g_mu-mu)>=1e-6)) && (g_proc_id == 0)){
      fprintf(stderr, "Warning! Parameters beta, kappa or mu are inconsistent with file %s!\n", filename);
    }
    if((l!=LX*g_nproc_x)||(t!=T*g_nproc_t)){
      fprintf(stderr, "Error! spinorfield %s was produced for a different lattice size!\nAborting!\n", filename);
      return(-3);
    }
    /* We do not need a seek here, see write_gauge_field... */
    for(x = 0; x < LX; x++) {
      for(y = 0; y < LY; y++) {
	for(z = 0; z < LZ; z++) {
#if (defined MPI)
	  fseek(ifs, position +
		(g_proc_coords[0]*T+
		 (((g_proc_coords[1]*LX+x)*g_nproc_y*LY+g_proc_coords[2]*LY+y)*g_nproc_z*LZ
		  + g_proc_coords[3]*LZ+z)*T*g_nproc_t)*sizeof(spinor),
		SEEK_SET);
#endif
	  for(t = 0; t < T; t++){
	    i = g_lexic2eosub[ g_ipt[t][x][y][z] ];
	    if((t+x+y+z+
		g_proc_coords[3]*LZ+g_proc_coords[2]*LY
		+g_proc_coords[0]*T+g_proc_coords[1]*LX)%2==0) {
	      p = s;
	    }
	    else {
	      p = r;
	    }
#ifndef WORDS_BIGENDIAN
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
      fclose(ifs);
      return(-1);
    }
  }
  else{
    fclose(ifs);
    return(-2);
  }
  fclose(ifs);
  return(0);
}

/*************************************************
 *
 * This routine writes an even or odd spinor-field
 * so really of size VOLUME/2
 *
 * used for instance for storing eigenvectors
 * of the precoditioned matrix
 *
 *************************************************/

int write_eospinor(spinor * const s, char * filename, 
		   const double evalue, const double prec, const int nstore) {
  FILE * ofs = NULL;
  LimeWriter * limewriter = NULL;
  LimeRecordHeader * limeheader = NULL;
  int x, X, y, Y, z, Z, t, t0, tag=0, id=0, i=0;
  int ME_flag=0, MB_flag=0, status=0;
  spinor tmp[1];
  int coords[4];
  char message[500];
  n_uint64_t bytes;
#ifdef MPI
  MPI_Status mpistatus;
#endif

  if(g_cart_id == 0){  
    if(g_kappa > 0. || g_kappa < 0.) {
      sprintf(message,"\n eigenvalue = %e\n prec = %e\n conf nr = %d\n beta = %f, kappa = %f, mu = %f, c2_rec = %f\n hmcversion = %s", 
	      evalue, prec, nstore, g_beta, g_kappa, g_mu/2./g_kappa, g_rgi_C1, PACKAGE_VERSION);
    }
    else {
      sprintf(message,"\n eigenvalue = %e\n prec = %e\n conf nr = %d\n beta = %f, kappa = %f, 2*kappa*mu = %f, c2_rec = %f\n hmcversion = %s", 
	      evalue, prec, nstore, g_beta, g_kappa, g_mu, g_rgi_C1, PACKAGE_VERSION);
    }
    bytes = strlen( message );
    
    if((ofs = fopen(filename, "w")) == (FILE*)NULL) {
      fprintf(stderr, "Error writing eigenvector to file %s!\n", filename);
      return(-1);
    }
    limewriter = limeCreateWriter( ofs );
    if(limewriter == (LimeWriter*)NULL) {
      fprintf(stderr, "LIME error in file %s for writing!\n Aboring...\n", filename);
#ifdef MPI
      MPI_Abort(MPI_COMM_WORLD, 1);
      MPI_Finalize();
#endif
      exit(500);
    }

    limeheader = limeCreateHeader(MB_flag, ME_flag, "xlf-info", bytes);
    status = limeWriteRecordHeader( limeheader, limewriter);
    if(status < 0 ) {
      fprintf(stderr, "LIME write header (xlf-info) error %d\n", status);
#ifdef MPI
      MPI_Abort(MPI_COMM_WORLD, 1);
      MPI_Finalize();
#endif
      exit(500);
    }
    limeDestroyHeader( limeheader );
    limeWriteRecordData(message, &bytes, limewriter);

    bytes = LX*g_nproc_x*LY*g_nproc_y*LZ*g_nproc_z*T*g_nproc_t*sizeof(spinor)/2;
    MB_flag=0; ME_flag=1;
    limeheader = limeCreateHeader(MB_flag, ME_flag, "eospinor-binary-data", bytes);
    status = limeWriteRecordHeader( limeheader, limewriter);
    if(status < 0 ) {
      fprintf(stderr, "LIME write header (eospinor-binary-data) error %d\n", status);
#ifdef MPI
      MPI_Abort(MPI_COMM_WORLD, 1);
      MPI_Finalize();
#endif
      exit(500);
    }
    limeDestroyHeader( limeheader );
  }

  bytes = sizeof(spinor);
  for(x = 0; x < LX*g_nproc_x; x++){
    X = x - g_proc_coords[1]*LX;
    coords[1] = x / LX;
    for(y = 0; y < LY*g_nproc_y; y++){
      Y = y - g_proc_coords[2]*LY;
      coords[2] = y / LY;
      for(z = 0; z < LZ*g_nproc_z; z++){
	Z = z - g_proc_coords[3]*LZ;
	coords[3] = z / LZ;
	for(t0 = 0; t0 < T*g_nproc_t; t0++){
	  t = t0 - T*g_proc_coords[0];
	  coords[0] = t0 / T;
#ifdef MPI
	  MPI_Cart_rank(g_cart_grid, coords, &id);
#endif
	  i = g_lexic2eosub[ g_ipt[t][X][Y][Z] ];
	  if((t+X+Y+Z+g_proc_coords[3]*LZ+g_proc_coords[2]*LY 
	      + g_proc_coords[0]*T+g_proc_coords[1]*LX)%2 == 0) {
	    if(g_cart_id == 0) {
	      if(g_cart_id == id) {
#ifndef WORDS_BIGENDIAN
		byte_swap_assign(tmp, s + i , sizeof(spinor)/8);
		status = limeWriteRecordData((void*)tmp, &bytes, limewriter);
#else
		status = limeWriteRecordData((void*)(s+i), &bytes, limewriter);
#endif
	      }
#ifdef MPI
	      else {
		MPI_Recv(tmp, sizeof(spinor)/8, MPI_DOUBLE, id, tag, g_cart_grid, &mpistatus);
		status = limeWriteRecordData((void*)tmp, &bytes, limewriter);
	      }
#endif
	      if(status < 0 ) {
		fprintf(stderr, "LIME write error %d\n", status);
#ifdef MPI
		MPI_Abort(MPI_COMM_WORLD, 1);
		MPI_Finalize();
#endif
		exit(500);
	      }
	    }
#ifdef MPI
	    else {
	      if(g_cart_id == id) {
#  ifndef WORDS_BIGENDIAN
		byte_swap_assign(tmp, s + i, sizeof(spinor)/8);
		MPI_Send((void*) tmp, sizeof(spinor)/8, MPI_DOUBLE, 0, tag, g_cart_grid);
#  else
		MPI_Send((void*) (s + i), sizeof(spinor)/8, MPI_DOUBLE, 0, tag, g_cart_grid);
#  endif		  
	      }
	    }
#endif
	    tag++;
	  }
	}
#ifdef MPI
 	MPI_Barrier(g_cart_grid); 
#endif
	tag=0;
      }
    }
  }
  if(g_cart_id == 0) {
    if(ferror(ofs)) {
      fprintf(stderr, "Warning! Error while writing to file %s \n", filename);
    }
    limeDestroyWriter( limewriter );
    fflush(ofs);
    fclose(ofs);
  }
  return(0);
}

int read_eospinor(spinor * const s, char * filename) {
  FILE * ifs;
  int t, x, y , z, i = 0, status=0;
  n_uint64_t bytes;
  char * header_type;
  LimeReader * limereader;
#ifdef MPI
  int position;
#endif
#ifndef WORDS_BIGENDIAN
  spinor tmp[1];
#endif
  
  if((ifs = fopen(filename, "r")) == (FILE*)NULL) {
    if(g_proc_id == 0) {
      fprintf(stderr, "Error opening file %s\n", filename);
    }
    return(-1);
  }

  limereader = limeCreateReader( ifs );
  if( limereader == (LimeReader *)NULL ) {
    if(g_proc_id == 0) {
      fprintf(stderr, "Unable to open LimeReader\n");
    }
    return(-1);
  }
  while( (status = limeReaderNextRecord(limereader)) != LIME_EOF ) {
    if(status != LIME_SUCCESS ) {
      fprintf(stderr, "limeReaderNextRecord returned error with status = %d!\n", status);
      status = LIME_EOF;
      break;
    }
    header_type = limeReaderType(limereader);
    if(!strcmp("eospinor-binary-data",header_type)) break;
  }
  if(status == LIME_EOF) {
    if(g_proc_id == 0) {
      fprintf(stderr, "no eospinor-binary-data record found in file %s\n",filename);
    }
    limeDestroyReader(limereader);
    fclose(ifs);
    return(-1);
  }
  bytes = limeReaderBytes(limereader);
  if((int)bytes != LX*g_nproc_x*LY*g_nproc_y*LZ*g_nproc_z*T*g_nproc_t*sizeof(spinor)/2) {
    if(g_proc_id == 0) {
      fprintf(stderr, "wrong length in eospinor: %d. Aborting read!\n", (int)bytes);
    }
    return(-1);
  }

  bytes = sizeof(spinor);
  for(x = 0; x < LX; x++) {
    for(y = 0; y < LY; y++) {
      for(z = 0; z < LZ; z++) {
#if (defined MPI)
	limeReaderSeek(limereader, (n_uint64_t)
		       (g_proc_coords[0]*T+
			(((g_proc_coords[1]*LX+x)*g_nproc_y*LY+g_proc_coords[2]*LY+y)*g_nproc_z*LZ
			 + g_proc_coords[3]*LZ+z)*T*g_nproc_t)*sizeof(spinor)/2,
		       SEEK_SET);
#endif
	for(t = 0; t < T; t++){
	  i = g_lexic2eosub[ g_ipt[t][x][y][z] ];
	  if((t+x+y+z+
	      g_proc_coords[3]*LZ+g_proc_coords[2]*LY
	      +g_proc_coords[0]*T+g_proc_coords[1]*LX)%2==0) {
	    
#ifndef WORDS_BIGENDIAN
	    status = limeReaderReadData(tmp, &bytes, limereader);
	    byte_swap_assign(s + i, tmp, sizeof(spinor)/8);
#else
	    status = limeReaderReadData((s+i), &bytes, limereader);
#endif
	    if(status < 0 && status != LIME_EOR) {
	      fprintf(stderr, "LIME read error occured with status = %d while reading file %s!\n Aborting...\n", status, filename);
#ifdef MPI
	      MPI_Abort(MPI_COMM_WORLD, 1);
	      MPI_Finalize();
#endif
	      exit(500);
	    }
	  }
	}
      }
    }
  }
  limeDestroyReader(limereader);
  fclose(ifs);
  return(0);
}



int write_rlxd_state(char * filename, int * const _state, const int rlxdsize) {
  n_uint64_t len;
  LimeWriter * limewriter = NULL;
  LimeRecordHeader * limeheader = NULL;
  int * state;
  int i, status;
  FILE * ofs;

  state = (int*)calloc(1, rlxdsize*sizeof(int));
  for(i = 0; i < rlxdsize; i++) {
    state[i] = _state[i];
  }
  ofs = fopen(filename, "a");
  if(ofs == (FILE*)NULL) {
    fprintf(stderr, "Could not open file %s for writing!\n Aboring...\n", filename);
#ifdef MPI
    MPI_Abort(MPI_COMM_WORLD, 1);
    MPI_Finalize();
#endif
    exit(501);
  }
  limewriter = limeCreateWriter( ofs );
  if(limewriter == (LimeWriter*)NULL) {
    fprintf(stderr, "LIME error in file %s for writing!\n Aboring...\n", filename);
#ifdef MPI
    MPI_Abort(MPI_COMM_WORLD, 1);
    MPI_Finalize();
#endif
    exit(501);
  }
  
  len = rlxdsize*sizeof(int);
  limeheader = limeCreateHeader(1, 1, "ranluxd-state", len);
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
#ifndef WORDS_BIGENDIAN
  byte_swap(state, rlxdsize);
#endif
  status = limeWriteRecordData(state, &len, limewriter);
  if(status < 0 ) {
    fprintf(stderr, "LIME write error %d\n", status);
#ifdef MPI
    MPI_Abort(MPI_COMM_WORLD, 1);
    MPI_Finalize();
#endif
    exit(500);
  }
  fflush(ofs);
  fclose(ofs);
  free(state);
  return(0);
}

int read_rlxd_state(char * filename, int * state, const int rlxdsize) {
  int status;
  FILE * ifs;
  n_uint64_t len;
  char * header_type;
  LimeReader * limereader;

  ifs = fopen(filename, "r");
  if(ifs == (FILE *)NULL) {
    fprintf(stderr, "Could not open file %s\n Aborting...\n", filename);
#ifdef MPI
    MPI_Abort(MPI_COMM_WORLD, 1);
    MPI_Finalize();
#endif
    exit(500);
  }
  limereader = limeCreateReader( ifs );
  if( limereader == (LimeReader *)NULL ) {
    fprintf(stderr, "Unable to open LimeReader\n");
#ifdef MPI
    MPI_Abort(MPI_COMM_WORLD, 1);
    MPI_Finalize();
#endif
    exit(500);
  }
  while( (status = limeReaderNextRecord(limereader)) != LIME_EOF ) {
    if(status != LIME_SUCCESS ) {
      fprintf(stderr, "limeReaderNextRecord returned error with status = %d!\n", status);
      status = LIME_EOF;
      break;
    }
    header_type = limeReaderType(limereader);
    if(!strcmp("ranluxd-state",header_type)) break;
  }
  if(status == LIME_EOF) {
    fprintf(stderr, "no ildg-binary-data record found in file %s\n",filename);
    fprintf(stderr, "trying old deprecated file format!\n");
    limeDestroyReader(limereader);
    fclose(ifs);
    return(-1);
  }
  len = limeReaderBytes(limereader);
  if((int)len != rlxdsize*sizeof(int)) {
    fprintf(stderr, "Wrong size for ranluxd-state (len=%d!=%d)  %s\n", 
	    (int)len, (int)rlxdsize*sizeof(int), filename);
    fprintf(stderr, "Aborting...!\n");
    fflush( stdout );
#ifdef MPI
    MPI_Abort(MPI_COMM_WORLD, 1);
    MPI_Finalize();
#endif
    exit(501);
  }
  status = limeReaderReadData(state, &len, limereader);
#ifndef WORDS_BIGENDIAN
  byte_swap(state, rlxdsize);
#endif
  if(status < 0 && status != LIME_EOR) {
    fprintf(stderr, "LIME read error occured with status = %d while reading file %s!\n Aborting...\n", 
	    status, filename);
#ifdef MPI
    MPI_Abort(MPI_COMM_WORLD, 1);
    MPI_Finalize();
#endif
    exit(500);
  }
  limeDestroyReader(limereader);
  fclose(ifs);
  return(0);
}


int read_spinorfield_cm_single(spinor * const s, spinor * const r, char * filename, 
			       const int ts, const int vol) {
  /*
   * ts is the number of the timeslice to be used
   *    if ts < 0 read a volume source
   *
   * if ts >= 0 and vol > 0 the file is a volume file
   * but only one timeslice should be read
   */

  FILE * ifs;
  int t, x, y , z, i = 0;
  spinor * p = NULL;
  float tmp[24];

  ifs = fopen(filename, "r");
  if(ifs == (FILE *)NULL) {
    return(-1);
  }

  for(x = 0; x < LX; x++) {
    for(y = 0; y < LY; y++) {
      for(z = 0; z < LZ; z++) {
#if (defined MPI)
	fseek(ifs,
	      (g_proc_coords[0]*T+
	       (((g_proc_coords[1]*LX+x)*g_nproc_y*LY+g_proc_coords[2]*LY+y)*g_nproc_z*LZ
		+ g_proc_coords[3]*LZ+z)*T*g_nproc_t)*sizeof(spinor)/2,
	      SEEK_SET);
#endif
	for(t = 0; t < T; t++) {

	  i = g_lexic2eosub[ g_ipt[t][x][y][z] ];
	  if((t+x+y+z+
	      g_proc_coords[0]*T+g_proc_coords[1]*LX+
	      g_proc_coords[2]*LY+g_proc_coords[3]*LZ)%2==0) {
	    p = s;
	  }
	  else {
	    p = r;
	  }
	  
	  if(ts == t || ts < 0 || ts >= T){
	    /* Read the data */
	    fread(tmp, sizeof(spinor)/2, 1, ifs);

            /* Test if we read the data with the correct endian order */
            if(isnan(tmp[0]) || isnan(tmp[1]) || isnan(tmp[2]) || isnan(tmp[3]) || isnan(tmp[4]) || isnan(tmp[5]) ||
            isnan(tmp[6]) || isnan(tmp[7]) || isnan(tmp[8]) || isnan(tmp[9]) || isnan(tmp[10]) || isnan(tmp[11]) ||
            isnan(tmp[12]) || isnan(tmp[13]) || isnan(tmp[14]) || isnan(tmp[15]) || isnan(tmp[16]) || isnan(tmp[17]) ||
            isnan(tmp[18]) || isnan(tmp[19]) || isnan(tmp[20]) || isnan(tmp[21]) || isnan(tmp[22]) || isnan(tmp[23]))
            {
              if(g_proc_id == 0)
              {
                if(big_endian())
                  printf("\nBig endian order gives some NaN. Trying little endian order instead...\n\n");
                else
                  printf("\nLittle endian order gives some NaN. Trying big endian order instead...\n\n");
              }

              fclose(ifs);
              return read_spinorfield_cm_swap_single(s,r,filename,ts,vol);
            }
	    single2double_cm(p+i, tmp);
	  }
	  else {
	    if(vol > 0) {
	      fread(tmp, sizeof(spinor)/2, 1, ifs);
	    }
	    /* Padding with zeros */
	    zero_spinor(p+i);
	  }
	}
      }
    }
  }
  fclose(ifs);
  return(0);
}

int read_spinorfield_cm_swap_single(spinor * const s, spinor * const r, char * filename,
                               const int ts, const int vol) {
  /*
   * ts is the number of the timeslice to be used
   *    if ts < 0 read a volume source
   *
   * if ts >= 0 and vol > 0 the file is a volume file
   * but only one timeslice should be read
   */

  FILE * ifs;
  int t, x, y , z, i = 0;
  spinor * p = NULL;
  float tmp[24];

  ifs = fopen(filename, "r");
  if(ifs == (FILE *)NULL) {
    fprintf(stderr, "Could not open file %s\n Aborting...\n", filename);
#ifdef MPI
    MPI_Abort(MPI_COMM_WORLD, 1);
    MPI_Finalize();
#endif
    exit(500);
  }

  for(x = 0; x < LX; x++) {
    for(y = 0; y < LY; y++) {
      for(z = 0; z < LZ; z++) {
#if (defined MPI)
        fseek(ifs,
              (g_proc_coords[0]*T+
               (((g_proc_coords[1]*LX+x)*g_nproc_y*LY+g_proc_coords[2]*LY+y)*g_nproc_z*LZ
                + g_proc_coords[3]*LZ+z)*T*g_nproc_t)*sizeof(spinor)/2,
              SEEK_SET);
#endif
        for(t = 0; t < T; t++) {

          i = g_lexic2eosub[ g_ipt[t][x][y][z] ];
          if((t+x+y+z+
              g_proc_coords[0]*T+g_proc_coords[1]*LX+
              g_proc_coords[2]*LY+g_proc_coords[3]*LZ)%2==0) {
            p = s;
          }
          else {
            p = r;
          }

          if(ts == t || ts < 0 || ts >= T){
            /* Read the data */
            fread(tmp, sizeof(spinor)/2, 1, ifs);

            /* Swap and convert from single to double precision */
            byte_swap_assign_single2double(p+i, tmp, sizeof(spinor)/8);
          }
          else {
            if(vol > 0) {
              fread(tmp, sizeof(spinor)/2, 1, ifs);
            }
            /* Padding with zeros */
            zero_spinor(p+i);
          }
        }
      }
    }
  }
  fclose(ifs);
  return(0);
}


int write_spinorfield_cm_single(spinor * const s, spinor * const r, char * filename) {

  FILE * ofs;
  int t, x, y , z, i = 0;
  int t0, X, Y, Z, id = 0;
  spinor * p = NULL;
  float tmp[24];
  int coords[4];
#ifdef MPI
  int  tag = 0;
  MPI_Status status;
#endif
  
  if(g_cart_id == 0) {
    ofs = fopen(filename, "w");
  }

  for(x = 0; x < LX*g_nproc_x; x++) {
    X = x - LX*g_proc_coords[1];
    coords[1] = x / LX;
    for(y = 0; y < LY*g_nproc_y; y++) {
      Y = y - LY*g_proc_coords[2];
      coords[2] = y / LY;
      for(z = 0; z < LZ*g_nproc_z; z++) {
	Z = z - LZ*g_proc_coords[3];
	coords[3] = z / LZ;
	for(t0 = 0; t0 < T*g_nproc_t; t0++) {
	  t = t0 - T*g_proc_coords[0];
	  coords[0] = t0 / T;
#ifdef MPI
	  MPI_Cart_rank(g_cart_grid, coords, &id);
#endif
	  if(g_cart_id == id) {
	    i = g_lexic2eosub[ g_ipt[t][X][Y][Z] ];
	    if((t+X+Y+Z+g_proc_coords[3]*LZ+g_proc_coords[2]*LY 
		+ g_proc_coords[0]*T+g_proc_coords[1]*LX)%2 == 0) {
	      p = s;
	    }
	    else {
	      p = r;
	    }
	  }
	  if(g_cart_id == 0){
	    if(g_cart_id == id) {
	      double2single_cm(tmp, p + i);
	    }
#ifdef MPI
	    else {
	      MPI_Recv(tmp, sizeof(spinor)/8, MPI_FLOAT, id, tag, g_cart_grid, &status);
	    }
#endif
#ifndef WORDS_BIGENDIAN
            byte_swap(tmp,sizeof(spinor)/8);
#endif
	    fwrite(tmp, sizeof(float), 24, ofs);
	  }
#ifdef MPI
	  else {
	    if(g_cart_id == id) {
	      double2single_cm(tmp, p + i);
	      MPI_Send((void*) tmp, sizeof(spinor)/8, MPI_FLOAT, 0, tag, g_cart_grid);
	    }
	  }
	  tag++;
#endif
	}
#ifdef MPI
	MPI_Barrier(g_cart_grid); 
	tag=0;
#endif
      }
    }
  }
  if(g_cart_id == 0) {
    fclose(ofs);
  }
  return(0);
}




int write_first_messages(FILE * parameterfile, const int integtyp, const int inv) {

  if(inv != 1) {
    printf("# This is the hmc code for twisted Mass Wilson QCD\n\nVersion %s\n", PACKAGE_VERSION);
    fprintf(parameterfile, 
	    "# This is the hmc code for twisted Mass Wilson QCD\n\nVersion %s\n", PACKAGE_VERSION);
  }
  else {
    printf("# This is the invert code for twisted Mass Wilson QCD\n\nVersion %s\n", PACKAGE_VERSION);
    fprintf(parameterfile, 
	    "# This is the invert code for twisted Mass Wilson QCD\n\nVersion %s\n", PACKAGE_VERSION);
  }
#ifdef SSE
  printf("# The code was compiled with SSE instructions\n");
  fprintf(parameterfile, 
	  "# The code was compiled with SSE instructions\n");
#endif
#ifdef SSE2
  printf("# The code was compiled with SSE2 instructions\n");
  fprintf(parameterfile, 
	  "# The code was compiled with SSE2 instructions\n");
#endif
#ifdef SSE3
  printf("# The code was compiled with SSE3 instructions\n");
  fprintf(parameterfile, 
	  "# The code was compiled with SSE3 instructions\n");
#endif
#ifdef P4
  printf("# The code was compiled for Pentium4\n");
  fprintf(parameterfile, 
	  "# The code was compiled for Pentium4\n");
#endif
#if (defined BGL && !defined BGP)
  printf("# The code was compiled for Blue Gene/L\n");
  fprintf(parameterfile, 
	  "# The code was compiled for Blue Gene/L\n");
#  if (defined _USE_BGLDRAM)
  printf("# The code was compiled for Blue Gene/L dram window\n");
  fprintf(parameterfile, 
	  "# The code was compiled for Blue Gene/L dram window\n");
#  endif
#endif
#ifdef BGP
  printf("# The code was compiled for Blue Gene/P\n");
  fprintf(parameterfile,
          "# The code was compiled for Blue Gene/P\n");
#endif
#ifdef OPTERON
  printf("# The code was compiled for AMD Opteron\n");
  fprintf(parameterfile,
	  "# The code was compiled for AMD Opteron\n");
#endif
#ifdef _GAUGE_COPY
  printf("# The code was compiled with -D_GAUGE_COPY\n");
  fprintf(parameterfile,
	  "# The code was compiled with -D_GAUGE_COPY\n");
#endif
#ifdef _USE_HALFSPINOR
  printf("# The code was compiled with -D_USE_HALFSPINOR\n");
  fprintf(parameterfile,
	  "# The code was compiled with -D_USE_HALFSPINOR\n");
#endif
#ifdef _USE_SHMEM
  printf("# the code was compiled with -D_USE_SHMEM\n");
  fprintf(parameterfile,
         "# the code was compiled with -D_USE_SHMEM\n");
#  ifdef _PERSISTENT
  printf("# the code was compiled for persistent MPI calls (halfspinor only)\n");
  fprintf(parameterfile,
         "# the code was compiled for persistent MPI calls (halfspinor only)\n");
#  endif
#endif
#ifdef MPI
#  ifdef _NON_BLOCKING
  printf("# the code was compiled for non-blocking MPI calls (spinor and gauge)\n");
  fprintf(parameterfile,
         "# the code was compiled for non-blocking MPI calls (spinor and gauge)\n");
#  endif

#endif
  printf("# The lattice size is %d x %d x %d x %d\n",
	 (int)(T*g_nproc_t), (int)(LX*g_nproc_x), (int)(LY*g_nproc_y), (int)(LZ*g_nproc_z));
  printf("# The local lattice size is %d x %d x %d x %d\n", 
      (int)(T), (int)(LX), (int)(LY),(int) LZ);
  if(even_odd_flag)
  {
    printf("# Even/odd preconditioning was used\n");
    fprintf(parameterfile, "# Even/odd preconditioning was used\n");
  }
  else
  {
    printf("# Even/odd preconditioning was not used\n");
    fprintf(parameterfile, "# Even/odd preconditioning was not used\n");
  }
  printf("# beta = %f , kappa= %f\n", g_beta, g_kappa);
  printf("# boundary conditions for fermion fields (t,x,y,z) * pi: %f %f %f %f \n",X0,X1,X2,X3);
  if(inv != 1) {
    printf("# mus = %f, %f, %f\n", g_mu1, g_mu2, g_mu3);
    printf("# int_n_gauge = %d, int_n_ferm1 = %d, int_n_ferm2 = %d, int_n_ferm3 = %d\n", 
	   int_n[0], int_n[1], int_n[2], int_n[3]);
    printf("# g_rgi_C0 = %f, g_rgi_C1 = %f\n", g_rgi_C0, g_rgi_C1);
    printf("# Number of pseudo-fermion fields: %d\n", g_nr_of_psf);
    printf("# g_eps_sq_force = %e, g_eps_sq_acc = %e\n", g_eps_sq_force, g_eps_sq_acc);
    printf("# Integration scheme: ");
    if(integtyp == 1) printf("leap-frog (single time scale)\n");
    if(integtyp == 2) printf("Sexton-Weingarten (single time scale)\n");
    if(integtyp == 3) printf("leap-frog (multiple time scales)\n");
    if(integtyp == 4) printf("Sexton-Weingarten (multiple time scales)\n");
    if(integtyp == 5) printf("higher order and leap-frog (multiple time scales)\n");
    if(integtyp == 6) printf("second order minimal norm (velocity version, multiple time scales)\n");
    if(integtyp == 7) printf("second order minimal norm (position version, multiple time scales)\n");
    printf("# Using %s precision for the inversions!\n", 
	   g_relative_precision_flag ? "relative" : "absolute");
    printf("# Using in chronological inverter for spinor_field 1,2,3 a history of %d, %d, %d, respectively\n", 
	   g_csg_N[0], g_csg_N[2], g_csg_N[4]);
  }

  fprintf(parameterfile, "# The lattice size is %d x %d x %d x %d\n", (int)(g_nproc_t*T), (int)(g_nproc_x*LX), 
	  (int)(g_nproc_y*LY), (int)(g_nproc_z*LZ));
  fprintf(parameterfile, "# The local lattice size is %d x %d x %d x %d\n", (int)(T), (int)(LX), (int)(LY), (int)(LZ));
  fprintf(parameterfile, "# g_beta = %f , g_kappa= %f, g_kappa*csw/8= %f \n",g_beta,g_kappa,g_ka_csw_8);
  fprintf(parameterfile, "# boundary conditions for fermion fields (t,x,y,z) * pi: %f %f %f %f \n",X0,X1,X2,X3);
  if(inv != 1) {
    fprintf(parameterfile, "# ITER_MAX_BCG=%d, EPS_SQ0=%e, EPS_SQ1=%e EPS_SQ2=%e, EPS_SQ3=%e \n"
	    ,ITER_MAX_BCG,EPS_SQ0,EPS_SQ1,EPS_SQ2,EPS_SQ3);
    fprintf(parameterfile, "# g_eps_sq_force = %e, g_eps_sq_acc = %e\n", g_eps_sq_force, g_eps_sq_acc);
    fprintf(parameterfile, "# dtau=%f, Nsteps=%d, Nmeas=%d, Nskip=%d, integtyp=%d, nsmall=%d \n",
	    dtau,Nsteps,Nmeas,Nskip,integtyp,nsmall);
    fprintf(parameterfile, "# mu = %f, mu2=%f, mu3=%f\n ", g_mu, g_mu2, g_mu3);
    fprintf(parameterfile, "# int_n_gauge = %d, int_n_ferm1 = %d, int_n_ferm2 = %d, int_n_ferm3 = %d\n ", 
	    int_n[0], int_n[1], int_n[2], int_n[3]);
    fprintf(parameterfile, "g_rgi_C0 = %f, g_rgi_C1 = %f\n", g_rgi_C0, g_rgi_C1);
    fprintf(parameterfile, "# Number of pseudo-fermion fields: %d\n", g_nr_of_psf);
    fprintf(parameterfile, "# Integration scheme: ");
    if(integtyp == 1) fprintf(parameterfile, "leap-frog (single time scale)\n");
    if(integtyp == 2) fprintf(parameterfile, "Sexton-Weingarten (single time scale)\n");
    if(integtyp == 3) fprintf(parameterfile, "leap-frog (multiple time scales)\n");
    if(integtyp == 4) fprintf(parameterfile, "Sexton-Weingarten (multiple time scales)\n");
    if(integtyp == 5) fprintf(parameterfile, "higher order and leap-frog (multiple time scales)\n");
    if(integtyp == 6) fprintf(parameterfile, "second order minimal norm (velocity version, multiple time scales)\n");
    if(integtyp == 7) fprintf(parameterfile, "second order minimal norm (position version, multiple time scales)\n");
    fprintf(parameterfile, "# Using %s precision for the inversions!\n", 
	    g_relative_precision_flag ? "relative" : "absolute");
    fprintf(parameterfile, "# Using in chronological inverter for spinor_field 1,2,3 a history of %d, %d, %d, respectively\n", 
	    g_csg_N[0], g_csg_N[2], g_csg_N[4]);
  }
  if(inv == 1) {
    printf("# beta = %f, mu = %f, kappa = %f\n", g_beta, g_mu/2./g_kappa, g_kappa);
    fprintf(parameterfile,
	    "# beta = %f, mu = %f, kappa = %f\n", g_beta, g_mu/2./g_kappa, g_kappa);
  }
  fflush(stdout); fflush(parameterfile);
  return(0);
}


int write_ildg_format_xml(char *filename, LimeWriter * limewriter, const int prec){
  FILE * ofs;
  n_uint64_t bytes, bytes_left, bytes_to_copy;
  int MB_flag=1, ME_flag=0, status=0;
  LimeRecordHeader * limeheader;
  char buf[MAXBUF];

  ofs = fopen(filename, "w");

  fprintf(ofs, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
  fprintf(ofs, "<ildgFormat xmlns=\"http://www.lqcd.org/ildg\"\n");
  fprintf(ofs, "            xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n");
  fprintf(ofs, "            xsi:schemaLocation=\"http://www.lqcd.org/ildg filefmt.xsd\">\n");
  fprintf(ofs, "  <version> 1.0 </version>\n");
  fprintf(ofs, "  <field> su3gauge </field>\n");
  if(prec == 0) {
    fprintf(ofs, "  <precision> 64 </precision>\n");
  }
  else {
    fprintf(ofs, "  <precision> 32 </precision>\n");
  }
  fprintf(ofs, "  <lx> %d </lx>\n", LX*g_nproc_x);
  fprintf(ofs, "  <ly> %d </ly>\n", LY*g_nproc_y);
  fprintf(ofs, "  <lz> %d </lz>\n", LZ*g_nproc_z);
  fprintf(ofs, "  <lt> %d </lt>\n", T*g_nproc_t);
  fprintf(ofs, "</ildgFormat>");
  fclose( ofs );
  ofs = fopen(filename, "r");
  bytes = file_size(ofs);
  
  limeheader = limeCreateHeader(MB_flag, ME_flag, "ildg-format", bytes);
  if(limeheader == (LimeRecordHeader*)NULL) {
    fprintf(stderr, "LIME create header ildg-format error\n Aborting...\n");
    exit(500);
  }
  status = limeWriteRecordHeader( limeheader, limewriter);
  if(status < 0 ) {
    fprintf(stderr, "LIME write header ildg-format error %d\n Aborting...\n", status);
    exit(500);
  }
  limeDestroyHeader( limeheader );

  /* Buffered copy */
  bytes_left = bytes;
  while(bytes_left > 0){
    if(MAXBUF < bytes_left) {
      bytes_to_copy = MAXBUF;
    }
    else bytes_to_copy = bytes_left;
    if( bytes_to_copy != fread(buf,1,bytes_to_copy,ofs))
      {
	fprintf(stderr, "Error reading %s\n", filename);
	return EXIT_FAILURE;
      }
    
    status = limeWriteRecordData(buf, &bytes_to_copy, limewriter);
    if (status != 0) {
      fprintf(stderr, "LIME error writing ildg-format status = %d\n Aborting...\n", status);
      exit(500);
    }
    bytes_left -= bytes_to_copy;
  }
  

  fclose( ofs );
  return(0);
}

n_uint64_t file_size(FILE *fp)
{
  n_uint64_t oldpos = ftello(fp);
  n_uint64_t length;
  
  if (fseeko(fp, 0L,SEEK_END) == -1)
    return -1;
  
  length = ftello(fp);
  
  return ( fseeko(fp,oldpos,SEEK_SET) == -1 ) ? -1 : length;
  
}


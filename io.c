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
#ifdef MPI
# include<unistd.h> 
#endif
#include"global.h"
#include"su3.h"
#include"lime.h" 
#include"io.h"

#define MAXBUF 1048576

#define off_t n_uint64_t 

void byte_swap(void *ptr, int nmemb);
void byte_swap_assign(void * out_ptr, void * in_ptr, int nmemb);
int big_endian();
int write_ildg_format_xml(char *filename, LimeWriter * limewriter);
off_t file_size(FILE *fp);
void single2double_cm(spinor * const R, float * const S);
void double2single_cm(float * const S, spinor * const R);
void zero_spinor(spinor * const R);

#ifdef MPI
int write_lime_gauge_field(char * filename, const double plaq, const int counter){
  FILE * ofs;
  LimeWriter * limewriter = NULL;
  LimeRecordHeader * limeheader = NULL;
  /* Message end and Message begin flag */
  int ME_flag=0, MB_flag=0, status=0;
  int tag=0, t=0, x=0, y=0, z=0, id=0, X=0, tt=0, Y=0;
  MPI_Status mpi_status;
  char message[500];
  su3 tmp[4];
  int coords[3];
  off_t bytes;
  struct timeval t1;

  gettimeofday(&t1,NULL);
  if(g_kappa > 0. || g_kappa < 0.) {
    sprintf(message,"\n plaquette = %e\n trajectory nr = %d\n beta = %f, kappa = %f, mu = %f, c2_rec = %f\n time = %ld\n hmcversion = %s", 
	    plaq, counter, g_beta, g_kappa, g_mu/2./g_kappa, g_rgi_C1,t1.tv_sec, PACKAGE_VERSION);
  }
  else {
    sprintf(message,"\n plaquette = %e\n trajectory nr = %d\n beta = %f, kappa = %f, 2*kappa*mu = %f, c2_rec = %f", 
	    plaq, counter, g_beta, g_kappa, g_mu, g_rgi_C1);
  }
  bytes = strlen( message );
  if(g_cart_id == 0) {
    ofs = fopen(filename, "w");
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
    write_ildg_format_xml("temp.xml", limewriter);
    
    limeheader = limeCreateHeader(MB_flag, ME_flag, "xlf-info", bytes);
    status = limeWriteRecordHeader( limeheader, limewriter);
    if(status < 0 ) {
      fprintf(stderr, "LIME write header error %d\n", status);
      MPI_Abort(MPI_COMM_WORLD, 1);
      MPI_Finalize();
      exit(500);
    }
    limeDestroyHeader( limeheader );
    limeWriteRecordData(message, &bytes, limewriter);

    bytes = LX*g_nproc_x*LY*g_nproc_y*LZ*T*g_nproc_t*4*sizeof(su3);
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

  bytes = sizeof(su3);
  for(t = 0; t < T*g_nproc_t; t++) {
    tt = t - g_proc_coords[0]*T;
    coords[0] = t / T;
    for(z = 0; z < LZ; z++) {
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
	      byte_swap_assign(tmp, g_gauge_field[ g_ipt[tt][X][Y][z] ], 4*sizeof(su3)/8); 
	      status = limeWriteRecordData((void*)&tmp[1], &bytes, limewriter);
	      status = limeWriteRecordData((void*)&tmp[2], &bytes, limewriter);
	      status = limeWriteRecordData((void*)&tmp[3], &bytes, limewriter);
	      status = limeWriteRecordData((void*)&tmp[0], &bytes, limewriter);
#else
	      status = limeWriteRecordData((void*)&g_gauge_field[ g_ipt[tt][X][Y][z] ][1], &bytes, limewriter);
	      status = limeWriteRecordData((void*)&g_gauge_field[ g_ipt[tt][X][Y][z] ][2], &bytes, limewriter);
	      status = limeWriteRecordData((void*)&g_gauge_field[ g_ipt[tt][X][Y][z] ][3], &bytes, limewriter);
	      status = limeWriteRecordData((void*)&g_gauge_field[ g_ipt[tt][X][Y][z] ][0], &bytes, limewriter);
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
	      fprintf(stderr, "x %d, y %d, z %d, t %d (%d,%d,%d,%d)\n",x,y,z,t,X,Y,z,tt);
	      MPI_Abort(MPI_COMM_WORLD, 1);
	      MPI_Finalize();
	      exit(500);
	    }
	  }
	  else {
	    if(g_cart_id == id){
#ifndef WORDS_BIGENDIAN
	      byte_swap_assign(tmp, g_gauge_field[ g_ipt[tt][X][Y][z] ], 4*sizeof(su3)/8);
	      MPI_Send((void*) tmp, 4*sizeof(su3)/8, MPI_DOUBLE, 0, tag, g_cart_grid);
#else
	      MPI_Send((void*) g_gauge_field[ g_ipt[tt][X][Y][z] ], 4*sizeof(su3)/8, MPI_DOUBLE, 0, tag, g_cart_grid);
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
    fclose(ofs);
  }
  return(0);
}
#else
int write_lime_gauge_field(char * filename, const double plaq, const int counter){
  FILE * ofs;
  LimeWriter * limewriter;
  LimeRecordHeader * limeheader;
  /* Message end and Message begin flag */
  int ME_flag=0, MB_flag=0, status=0;
  int t=0, x=0, y=0, z=0;
  char message[100];
  off_t bytes;
#ifndef WORDS_BIGENDIAN
  su3 tmp[4];
#endif

  if(g_kappa > 0. || g_kappa < 0.) {
    sprintf(message,"\n plaquette = %e\n trajectory nr = %d\n beta = %f, kappa = %f, mu = %f, c2_rec = %f", 
	    plaq, counter, g_beta, g_kappa, g_mu/2./g_kappa, g_rgi_C1);
  }
  else {
    sprintf(message,"\n plaquette = %e\n trajectory nr = %d\n beta = %f, kappa = %f, 2*kappa*mu = %f, c2_rec = %f", 
	    plaq, counter, g_beta, g_kappa, g_mu, g_rgi_C1);
  }
  bytes = (off_t)strlen( message );
  ofs = fopen(filename, "w");
  if(ofs == (FILE*)NULL) {
    fprintf(stderr, "Could not open file %s for writing!\n Aboring...\n", filename);
    exit(500);
  }
  limewriter = limeCreateWriter( ofs );
  if(limewriter == (LimeWriter*)NULL) {
    fprintf(stderr, "LIME error in file %s for writing!\n Aboring...\n", filename);
    exit(500);
  }
  write_ildg_format_xml("temp.xml", limewriter);

  limeheader = limeCreateHeader(MB_flag, ME_flag, "xlf-info", bytes);
  status = limeWriteRecordHeader( limeheader, limewriter);
  if(status < 0 ) {
    fprintf(stderr, "LIME write header error %d\n", status);
    exit(500);
  }
  limeDestroyHeader( limeheader );
  limeWriteRecordData(message, &bytes, limewriter);

  bytes = (off_t)(LX*LY*LZ*T*4*sizeof(su3));
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
  fclose(ofs);
  return(0);
}
#endif

#ifdef MPI
int write_gauge_field_time_p(char * filename){
  FILE * ofs = NULL;
  int tag=0, t, x, y, z, id, X=0, tt=0, Y=0;
  MPI_Status status;
  su3 tmp[4];
  int coords[3];
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
      for(z = 0; z < LZ; z++){
	tag = 0;
	for(t = 0; t < T*g_nproc_t; t++){
	  tt = t - g_proc_coords[0]*T; 
	  coords[0] = t / T;
	  MPI_Cart_rank(g_cart_grid, coords, &id);
	  if(g_cart_id == 0){
	    if(g_cart_id == id){
#ifndef WORDS_BIGENDIAN
	      byte_swap_assign(tmp, g_gauge_field[ g_ipt[tt][X][Y][z] ], 4*sizeof(su3)/8);
	      fwrite(tmp, sizeof(su3), 4, ofs);
#else
	      fwrite(g_gauge_field[ g_ipt[tt][X][Y][z] ], 4*sizeof(su3), 1, ofs);
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
	      byte_swap_assign(tmp, g_gauge_field[ g_ipt[tt][X][Y][z] ], 4*sizeof(su3)/8);
	      MPI_Send((void*) tmp, 4*sizeof(su3)/8, MPI_DOUBLE, 0, tag, g_cart_grid);
#else
	      MPI_Send((void*) g_gauge_field[ g_ipt[tt][X][Y][z] ], 4*sizeof(su3)/8, MPI_DOUBLE, 0, tag, g_cart_grid);
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
#if (defined MPI)
	  fseek(ifs, position +
		(g_proc_coords[0]*T+
		 (((g_proc_coords[1]*LX+x)*g_nproc_y*LY + g_proc_coords[2]*LY + y)*LZ+z)*T*g_nproc_t)*4*sizeof(su3),
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

int read_lime_gauge_field(char * filename){
  FILE * ifs;
  int t, x, y, z, status;
  off_t bytes;
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
    fprintf(stderr, "no ildg-binary-data record found in file %s\n",filename);
    fprintf(stderr, "trying old deprecated file format!\n");
    limeDestroyReader(limereader);
    fclose(ifs);
    read_gauge_field_time_p(filename);
    return(0);
  }
  bytes = limeReaderBytes(limereader);
  if((int)bytes != LX*g_nproc_x*LY*g_nproc_y*LZ*T*g_nproc_t*4*sizeof(su3)) {
    fprintf(stderr, "Probably wrong lattice size or precision (bytes=%d) in file %s\n", (int)bytes, filename);
    fprintf(stderr, "Aborting...!\n");
    fflush( stdout );
#ifdef MPI
    MPI_Abort(MPI_COMM_WORLD, 1);
    MPI_Finalize();
#endif
    exit(501);
  }

  bytes = (off_t)4*sizeof(su3);
#ifdef WORDS_BIGENDIAN
  bytes = sizeof(su3);
#endif
  for(t = 0; t < T; t++){
    for(z = 0; z < LZ; z++){
      for(y = 0; y < LY; y++){
#if (defined MPI)
	limeReaderSeek(limereader,(off_t) 
		       (g_proc_coords[1]*LX + 
			(((g_proc_coords[0]*T+t)*LZ+z)*g_nproc_y*LY+g_proc_coords[2]*LY+y)*LX*g_nproc_x)*4*sizeof(su3),
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


int write_spinorfield_eo_time_p(spinor * const s, spinor * const r, char * filename, const int append){
  FILE * ofs = NULL;
  int x, X, y, Y, z, t, t0, tag=0, id=0, i=0;
  spinor tmp[1];
  int coords[3];
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
    for(y = 0; y < LY*g_nproc_y; y++){
      Y = y - g_proc_coords[2]*LY;
      coords[2] = y / LY;
      for(z = 0; z < LZ; z++){
	for(t0 = 0; t0 < T*g_nproc_t; t0++){
	  t = t0 - T*g_proc_coords[0];
	  coords[0] = t0 / T;
#ifdef MPI
	  MPI_Cart_rank(g_cart_grid, coords, &id);
#endif
	  if(g_cart_id == id) {
	    i = g_lexic2eosub[ g_ipt[t][X][Y][z] ];
	    if((t+X+Y+z+g_proc_coords[2]*LY+g_proc_coords[0]*T+g_proc_coords[1]*LX)%2 == 0) {
	      p = s;
	    }
	    else {
	      p = r;
	    }
	  }
	  if(g_cart_id == 0){
	    if(g_cart_id == id){
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
#ifndef WORDS_BIGENDIAN
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
#ifndef WORDS_BIGENDIAN
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
    for(x = 0; x < LX; x++) {
      for(y = 0; y < LY; y++) {
	for(z = 0; z < LZ; z++) {
#if (defined MPI)
	  fseek(ifs, position +
		(g_proc_coords[0]*T+
		 (((g_proc_coords[1]*LX+x)*g_nproc_y*LY+g_proc_coords[2]*LY+y)*LZ+z)*T*g_nproc_t)*sizeof(spinor),
		SEEK_SET);
#endif
	  for(t = 0; t < T; t++){
	    i = g_lexic2eosub[ g_ipt[t][x][y][z] ];
	    if((t+x+y+z+g_proc_coords[2]*LY+g_proc_coords[0]*T+g_proc_coords[1]*LX)%2==0) {
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
/*       errorhandler(107, filename); */
      fprintf(stderr, "Error while reading from file %s!\nAborting!\n", filename);
      fprintf(stderr, "%d %d\n", feof(ifs), ferror(ifs));
#ifdef MPI
      MPI_Abort(MPI_COMM_WORLD, 1);
      MPI_Finalize();
#endif
      exit(1);
    }
  }
  else{
/*     errorhandler(107, filename); */
    printf("Error while reading from file %s!\nAborting!\n", filename);
    exit(1);
  }
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

  for(x = 0; x < LX; x++) {
    for(y = 0; y < LY; y++) {
      for(z = 0; z < LZ; z++) {
	for(t = 0; t < T; t++) {

	  i = g_lexic2eosub[ g_ipt[t][x][y][z] ];
	  if((t+x+y+z+g_proc_coords[0]*T+g_proc_coords[1]*LX)%2==0) {
	    p = s;
	  }
	  else {
	    p = r;
	  }
	  
	  if(ts == t || ts < 0 || ts >= T){
	    /* Read the data */
	    fread(tmp, sizeof(spinor)/2, 1, ifs);
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

int write_spinorfield_cm_single(spinor * const s, spinor * const r, char * filename) {

  FILE * ofs;
  int t, x, y , z, i = 0;
  spinor * p = NULL;
  float tmp[24];

  ofs = fopen(filename, "w");

  for(x = 0; x < LX; x++) {
    for(y = 0; y < LY; y++) {
      for(z = 0; z < LZ; z++) {
	for(t = 0; t < T; t++) {
  	  i = g_lexic2eosub[ g_ipt[t][x][y][z] ];
	  if((t+x+y+z+g_proc_coords[0]*T+g_proc_coords[1]*LX)%2==0) {
	    p = s;
	  }
	  else {
	    p = r;
	  }
	  double2single_cm(tmp, p + i);
	  fwrite(tmp, sizeof(float), 24, ofs);
	}
      }
    }
  }
  fclose(ofs);
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

int write_ildg_format_xml(char *filename, LimeWriter * limewriter){
  FILE * ofs;
  off_t bytes, bytes_left, bytes_to_copy;
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
  fprintf(ofs, "  <precision> 64 </precision>\n");
  fprintf(ofs, "  <lx> %d </lx>\n", LX*g_nproc_x);
  fprintf(ofs, "  <ly> %d </ly>\n", LY);
  fprintf(ofs, "  <lz> %d </lz>\n", LZ);
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

off_t file_size(FILE *fp)
{
  off_t oldpos = ftello(fp);
  off_t length;
  
  if (fseeko(fp, 0L,SEEK_END) == -1)
    return -1;
  
  length = ftello(fp);
  
  return ( fseeko(fp,oldpos,SEEK_SET) == -1 ) ? -1 : length;
  
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

/* $Id$ */

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


int write_binary_spinor_data(spinor * const s, spinor * const r, LimeWriter * limewriter) {

  int x, X, y, Y, z, Z, t, t0, tag=0, id=0, i=0, status=0;
  spinor * p = NULL;
  spinor tmp[1];
  int coords[4];
  n_uint64_t bytes;

  bytes = sizeof(spinor);
  for(t0 = 0; t0 < T*g_nproc_t; t0++) {
    t = t0 - T*g_proc_coords[0];
    coords[0] = t0 / T;  
    for(z = 0; z < LZ*g_nproc_z; z++) {
      Z = z - g_proc_coords[3]*LZ;
      coords[3] = z / LZ;
      for(y = 0; y < LY*g_nproc_y; y++) {
	Y = y - g_proc_coords[2]*LY;
	coords[2] = y / LY;
	for(x = 0; x < LX*g_nproc_x; x++) {
	  X = x - g_proc_coords[1]*LX;
	  coords[1] = x / LX;
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
	      status = limeWriteRecordData((void*)tmp, &bytes, limewriter);
#else
	      status = limeWriteRecordData((void*)(p + i), &bytes, limewriter);
#endif
	    }
#ifdef MPI
	    else{
	      MPI_Recv(tmp, sizeof(spinor)/8, MPI_DOUBLE, id, tag, g_cart_grid, &status);
	      status = limeWriteRecordData((void*)tmp, &bytes, limewriter);
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
  return(0);
}

int read_binary_spinor_data(spinor * const s, spinor * const r, LimeReader * limereader) {
  int t, x, y , z, i = 0, status=0;
  n_uint64_t bytes;
  spinor * p = NULL;
#ifndef WORDS_BIGENDIAN
  spinor tmp[1];
#endif

  bytes = sizeof(spinor);
  for(t = 0; t < T; t++){
    for(z = 0; z < LZ; z++){
      for(y = 0; y < LY; y++){
#if (defined MPI)
	limeReaderSeek(limereader,(n_uint64_t) 
		       (g_proc_coords[1]*LX + 
			(((g_proc_coords[0]*T+t)*g_nproc_z*LZ+g_proc_coords[3]*LZ+z)*g_nproc_y*LY 
			 + g_proc_coords[2]*LY+y)*LX*g_nproc_x)*sizeof(spinor),
		       SEEK_SET);
#endif
	for(x = 0; x < LX; x++){
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
	  status = limeReaderReadData(tmp, &bytes, limereader);
	  byte_swap_assign(p + i, tmp, sizeof(spinor)/8);
#else
	  status = limeReaderReadData((p+i), &bytes, limereader);
#endif
	  if(status < 0 && status != LIME_EOR) {
	    return(-1);
	  }
	}
      }
    }
  }
  return(0);
}

int write_propagator_type(const int type, char * filename) {

  FILE * ofs = NULL;
  LimeWriter * limewriter = NULL;
  LimeRecordHeader * limeheader = NULL;
  int status = 0;
  int ME_flag=1, MB_flag=1;
  char message[500];
  n_uint64_t bytes;

  if(g_cart_id == 0) {
    ofs = fopen(filename, "w");

    if(ofs == (FILE*)NULL) {
      fprintf(stderr, "Could not open file %s for writing!\n Aboring...\n", filename);
#ifdef MPI
      MPI_Abort(MPI_COMM_WORLD, 1);
      MPI_Finalize();
#endif
      exit(500);
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


    if(type == 0) {
      sprintf(message,"DiracFermion_Sinks");
      bytes = strlen( message );
    }
    else if (type == 1) {
      sprintf(message,"DiracFermion_Source_Sink_Pairs");
      bytes = strlen( message );
    }
    else if (type == 2) {
      sprintf(message,"DiracFermion_ScalarSource_TwelveSink");
      bytes = strlen( message );
    }
    else if (type == 3) {
      sprintf(message,"DiracFermion_ScalarSource_FourSink");
      bytes = strlen( message );
    }

    limeheader = limeCreateHeader(MB_flag, ME_flag, "etmc-propagator-type", bytes);
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

int write_propagator_format(char * filename, const int prec) {
  FILE * ofs = NULL;
  LimeWriter * limewriter = NULL;
  LimeRecordHeader * limeheader = NULL;
  int status = 0;
  int ME_flag=1, MB_flag=1;
  char message[500];
  n_uint64_t bytes;
/*   char * message; */


  if(g_cart_id == 0) {
    ofs = fopen(filename, "a");

    if(ofs == (FILE*)NULL) {
      fprintf(stderr, "Could not open file %s for writing!\n Aboring...\n", filename);
#ifdef MPI
      MPI_Abort(MPI_COMM_WORLD, 1);
      MPI_Finalize();
#endif
      exit(500);
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

    sprintf(message, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<etmcFormat>\n<field>diracFermion</field>\n<precision>%d</precision>\n<lx>%d</lx>\n<ly>%d</ly>\n<lz>%d</lz>\n<lt>%d</lt>\n</etmcFormat>", prec, LX*g_nproc_x, LY*g_nproc_y, LZ*g_nproc_z, T*g_nproc_t);
    bytes = strlen( message );
    limeheader = limeCreateHeader(MB_flag, ME_flag, "etmc-propagator-type", bytes);
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

int write_lime_spinor(spinor * const s, spinor * const r, char * filename, 
		      const int append, const int prec) {

  FILE * ofs = NULL;
  LimeWriter * limewriter = NULL;
  LimeRecordHeader * limeheader = NULL;
  int status = 0;
  int ME_flag=0, MB_flag=0;
  n_uint64_t bytes;
#ifdef MPI
  MPI_Status mpistatus;
#endif

  write_propagator_format(filename, 64);

  if(g_cart_id == 0) {
    if(append) {
      ofs = fopen(filename, "a");
    }
    else {
      ofs = fopen(filename, "w");
    }
    if(ofs == (FILE*)NULL) {
      fprintf(stderr, "Could not open file %s for writing!\n Aboring...\n", filename);
#ifdef MPI
      MPI_Abort(MPI_COMM_WORLD, 1);
      MPI_Finalize();
#endif
      exit(500);
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
  
    bytes = LX*g_nproc_x*LY*g_nproc_y*LZ*g_nproc_z*T*g_nproc_t*sizeof(spinor);
    MB_flag=1; ME_flag=1;
    limeheader = limeCreateHeader(MB_flag, ME_flag, "scidac-binary-data", bytes);
    status = limeWriteRecordHeader( limeheader, limewriter);
    if(status < 0 ) {
      fprintf(stderr, "LIME write header (scidac-binary-data) error %d\n", status);
#ifdef MPI
      MPI_Abort(MPI_COMM_WORLD, 1);
      MPI_Finalize();
#endif
      exit(500);
    }
    limeDestroyHeader( limeheader );
  }

  write_binary_spinor_data(s, r, limewriter);

  if(g_cart_id == 0) {
    if(ferror(ofs)) {
      fprintf(stderr, "Warning! Error while writing to file %s \n", filename);
    }
    limeDestroyWriter( limewriter );
    fclose(ofs);
    fflush(ofs);
  }
  return(0);
}

int read_lime_spinor(spinor * const s, spinor * const r, char * filename) {
  FILE * ifs;
  int status=0;
  n_uint64_t bytes;
  char * header_type;
  LimeReader * limereader;
  
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
    if(!strcmp("scidac-binary-data",header_type)) break;
  }
  if(status == LIME_EOF) {
    if(g_proc_id == 0) {
      fprintf(stderr, "no scidac-binary-data record found in file %s\n",filename);
    }
    limeDestroyReader(limereader);
    fclose(ifs);
    return(-1);
  }
  bytes = limeReaderBytes(limereader);
  if((int)bytes != LX*g_nproc_x*LY*g_nproc_y*LZ*g_nproc_z*T*g_nproc_t*sizeof(spinor)) {
    if(g_proc_id == 0) {
      fprintf(stderr, "wrong length in eospinor: %d. Aborting read!\n", (int)bytes);
    }
    return(-1);
  }

  status = read_binary_spinor_data(s, r, limereader);

  if(status < 0) {
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

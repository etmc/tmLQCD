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
#include"propagator_io.h"
#include"dml.h"
#include"io.h"

/* write a one flavour propagator to file */

int write_propagator(spinor * const s, spinor * const r, char * filename, 
		      const int append, const int prec) {
  int err = 0;

  write_propagator_format(filename, prec, 1);
  err = write_lime_spinor(s, r, filename, append, prec);
  return(err);
}

/* write a one flavour source to file */

int write_source(spinor * const s, spinor * const r, char * filename, 
		 const int append, const int prec) {
  int err = 0;

  write_source_format(filename, prec, 1, T, LX, LY, LZ, 4, 3);
  err = write_lime_spinor(s, r, filename, append, prec);
  return(err);
}

/* write two flavour operator to file */

int write_double_propagator(spinor * const s, spinor * const r, 
			    spinor * const p, spinor * const q,
			    char * filename, const int append, const int prec) {
  int err = 0;

  /* we store strange component first, then charm */
  /* strange -> (mu_sigma - mu_delta)             */
  /* charm   -> (mu_sigma + mu_delta)             */
  /* they are in interchaged order in our code    */

  write_propagator_format(filename, prec, 2);
  err = write_lime_spinor(p, q, filename, append, prec);
  err += write_lime_spinor(s, r, filename, append, prec);
  return(err);
}

int write_binary_spinor_data(spinor * const s, spinor * const r, LimeWriter * limewriter,
				      const int prec, DML_Checksum * ans) {
  
  int x, X, y, Y, z, Z, t, t0, tag=0, id=0, i=0, status=0;
  spinor * p = NULL;
  spinor tmp[1];
  float tmp2[24];
  int coords[4];
  n_uint64_t bytes;
  DML_SiteRank rank;
#ifdef MPI
  MPI_Status mstatus;
#endif
  DML_checksum_init(ans);

  if(prec == 32) bytes = (n_uint64_t)sizeof(spinor)/2;
  else bytes = (n_uint64_t)sizeof(spinor);
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
	  if(g_cart_id == id) {
	    rank = (DML_SiteRank) (((t0*LZ*g_nproc_z + z)*LY*g_nproc_y + y)*LX*g_nproc_x + x);
	  }
	  if(g_cart_id == 0) {
	    if(g_cart_id == id) {
#ifndef WORDS_BIGENDIAN
	      if(prec == 32) {
		byte_swap_assign_double2single((float*)tmp2, p + i, sizeof(spinor)/8); 
		DML_checksum_accum(ans,rank,(char *) tmp2,sizeof(spinor)/2);
		status = limeWriteRecordData((void*)tmp2, &bytes, limewriter);
	      }
	      else {
		byte_swap_assign(tmp, p + i , sizeof(spinor)/8);
		DML_checksum_accum(ans,rank,(char *) tmp,sizeof(spinor));
		status = limeWriteRecordData((void*)tmp, &bytes, limewriter);
	      }
#else
	      if(prec == 32) {
		double2single((float*)tmp2, (p + i), sizeof(spinor)/8); 
		DML_checksum_accum(ans,rank,(char *) tmp2,sizeof(spinor)/2);
		status = limeWriteRecordData((void*)tmp2, &bytes, limewriter);
	      }
	      else {
		status = limeWriteRecordData((void*)(p + i), &bytes, limewriter);
		DML_checksum_accum(ans,rank,(char *) (p + i), sizeof(spinor));
	      }
#endif
	    }
#ifdef MPI
	    else{
	      if(prec == 32) {
		MPI_Recv((void*)tmp2, sizeof(spinor)/8, MPI_FLOAT, id, tag, g_cart_grid, &mstatus);
		DML_checksum_accum(ans,rank,(char *) tmp2, sizeof(spinor)/2);
		status = limeWriteRecordData((void*)tmp2, &bytes, limewriter);
	      }
	      else {
		MPI_Recv((void*)tmp, sizeof(spinor)/8, MPI_DOUBLE, id, tag, g_cart_grid, &mstatus);
		DML_checksum_accum(ans,rank,(char *) tmp, sizeof(spinor));
		status = limeWriteRecordData((void*)tmp, &bytes, limewriter);
	      }
	    }
#endif
	  }
#ifdef MPI
	  else{
	    if(g_cart_id == id){
#  ifndef WORDS_BIGENDIAN
	      if(prec == 32) {
		byte_swap_assign_double2single((float*)tmp2, p + i, sizeof(spinor)/8);
		MPI_Send((void*) tmp2, sizeof(spinor)/8, MPI_FLOAT, 0, tag, g_cart_grid);
	      }
	      else {
		byte_swap_assign(tmp, p + i, sizeof(spinor)/8);
		MPI_Send((void*) tmp, sizeof(spinor)/8, MPI_DOUBLE, 0, tag, g_cart_grid);
	      }
#  else
	      if(prec == 32) {
		double2single((float*)tmp2, (p + i), sizeof(spinor)/8); 
		MPI_Send((void*) tmp2, sizeof(spinor)/8, MPI_FLOAT, 0, tag, g_cart_grid);
	      }
	      else {
		MPI_Send((void*) (p + i), sizeof(spinor)/8, MPI_DOUBLE, 0, tag, g_cart_grid);
	      }
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

int read_binary_spinor_data(spinor * const s, spinor * const r, LimeReader * limereader, 
			    const double prec, DML_Checksum * ans) {
  int t, x, y , z, i = 0, status=0;
  n_uint64_t bytes;
  spinor * p = NULL;
  spinor tmp[1];
  float tmp2[24];
  DML_SiteRank rank;

  DML_checksum_init(ans);
  
  if(prec == 32) bytes = sizeof(spinor)/2;
  else bytes = sizeof(spinor);
  for(t = 0; t < T; t++){
    for(z = 0; z < LZ; z++){
      for(y = 0; y < LY; y++){
#if (defined MPI)
	limeReaderSeek(limereader,(n_uint64_t) 
		       (g_proc_coords[1]*LX + 
			(((g_proc_coords[0]*T+t)*g_nproc_z*LZ+g_proc_coords[3]*LZ+z)*g_nproc_y*LY 
			 + g_proc_coords[2]*LY+y)*LX*g_nproc_x)*bytes,
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
	  rank = (DML_SiteRank) (g_proc_coords[1]*LX + 
				 (((g_proc_coords[0]*T+t)*g_nproc_z*LZ+g_proc_coords[3]*LZ+z)*g_nproc_y*LY 
				  + g_proc_coords[2]*LY+y)*LX*g_nproc_x + x);
	  if(prec == 32) {
	    status = limeReaderReadData(tmp2, &bytes, limereader);
	    DML_checksum_accum(ans,rank,(char *) tmp2, sizeof(spinor)/2);
	  }
	  else {
	    status = limeReaderReadData(tmp, &bytes, limereader);
	    DML_checksum_accum(ans,rank,(char *) tmp, sizeof(spinor));
	  }
#ifndef WORDS_BIGENDIAN
	  if(prec == 32) {
	    byte_swap_assign_single2double(p+i, (float*)tmp2, sizeof(spinor)/8);
	  }
	  else {
	    byte_swap_assign(p + i, tmp, sizeof(spinor)/8);
	  }
#else
	  if(prec == 32) {
	    single2double(p + i, (float*)tmp2, sizeof(spinor)/8);
	  }
	  else memcpy(p+i, tmp, sizeof(spinor));
#endif
	  if(status < 0 && status != LIME_EOR) {
	    return(-1);
	  }
	}
      }
    }
  }
#ifdef MPI
  DML_checksum_combine(ans);
#endif
  return(0);
}

int write_source_type(const int type, char * filename) {
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
      fprintf(stderr, "LIME error in file %s for writing!\n Aborting...\n", filename);
#ifdef MPI
      MPI_Abort(MPI_COMM_WORLD, 1);
      MPI_Finalize();
#endif
      exit(500);
    }


    if(type == 0) {
      sprintf(message,"DiracFermion_Source");
      bytes = strlen( message );
    }
    else if (type == 1) {
      sprintf(message,"DiracFermion_ScalarSource");
      bytes = strlen( message );
    }
    else if (type == 2) {
      sprintf(message,"DiracFermion_FourScalarSource");
      bytes = strlen( message );
    }
    else if (type == 3) {
      sprintf(message,"DiracFermion_TwelveScalarSource");
      bytes = strlen( message );
    }

    limeheader = limeCreateHeader(MB_flag, ME_flag, "source-type", bytes);
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
      fprintf(stderr, "LIME error in file %s for writing!\n Aborting...\n", filename);
#ifdef MPI
      MPI_Abort(MPI_COMM_WORLD, 1);
      MPI_Finalize();
#endif
      exit(500);
    }


    if(type == 0) {
      sprintf(message,"DiracFermion_Sink");
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

    limeheader = limeCreateHeader(MB_flag, ME_flag, "propagator-type", bytes);
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

int write_propagator_format(char * filename, const int prec, const int no_flavours) {
  FILE * ofs = NULL;
  LimeWriter * limewriter = NULL;
  LimeRecordHeader * limeheader = NULL;
  int status = 0;
  int ME_flag=0, MB_flag=1;
  char message[500];
  n_uint64_t bytes;
  /*   char * message; */


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

    sprintf(message, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<etmcFormat>\n<field>diracFermion</field>\n<precision>%d</precision>\n<flavours>%d</flavours>\n<lx>%d</lx>\n<ly>%d</ly>\n<lz>%d</lz>\n<lt>%d</lt>\n</etmcFormat>", prec, no_flavours, LX*g_nproc_x, LY*g_nproc_y, LZ*g_nproc_z, T*g_nproc_t);
    bytes = strlen( message );
    limeheader = limeCreateHeader(MB_flag, ME_flag, "etmc-propagator-format", bytes);
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

int write_source_format(char * filename, const int prec, const int no_flavours,
			const int t, const int lx, const int ly, const int lz,
			const int is, const int ic) {
  FILE * ofs = NULL;
  LimeWriter * limewriter = NULL;
  LimeRecordHeader * limeheader = NULL;
  int status = 0;
  int ME_flag=0, MB_flag=1;
  char message[500];
  n_uint64_t bytes;
  /*   char * message; */


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

    sprintf(message, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<etmcFormat>\n<field>diracFermion</field>\n<precision>%d</precision>\n<flavours>%d</flavours>\n<lx>%d</lx>\n<ly>%d</ly>\n<lz>%d</lz>\n<lt>%d</lt>\n<spin>%d</spin>\n<colour>%d</colour>\n</etmcFormat>", prec, no_flavours, lx, ly, lz, t, is, ic);
    bytes = strlen( message );
    limeheader = limeCreateHeader(MB_flag, ME_flag, "etmc-source-format", bytes);
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
  DML_Checksum checksum;
#ifdef MPI
  MPI_Status mpistatus;
#endif


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
  
    bytes = (n_uint64_t)LX*g_nproc_x*LY*g_nproc_y*LZ*g_nproc_z*T*g_nproc_t*(n_uint64_t)(sizeof(spinor)*prec/64);
    MB_flag=0; ME_flag=1;
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

  status = write_binary_spinor_data(s, r, limewriter, prec, &checksum);
  if(g_proc_id == 0) {
    printf("# checksum for DiracFermion field written to file %s is %#x %#x\n", 
	   filename, checksum.suma, checksum.sumb);
  }

  if(g_cart_id == 0) {
    if(ferror(ofs)) {
      fprintf(stderr, "Warning! Error while writing to file %s \n", filename);
    }
    limeDestroyWriter( limewriter );
    fclose(ofs);
    fflush(ofs);
  }
  write_checksum(filename, &checksum);
  return(0);
}

int get_propagator_type(char * filename) {
  FILE * ifs;
  int status=0;
  n_uint64_t bytes;
  char * header_type;
  char tmp[500];
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
    if(strcmp("propagator-type",header_type) == 0) break;
  }
  if(status == LIME_EOF) {
    if(g_proc_id == 0) {
      fprintf(stderr, "no propagator-type record found in file %s\n",filename);
    }
    limeDestroyReader(limereader);
    fclose(ifs);
    return(-1);
  }
  bytes = limeReaderBytes(limereader);
  status = limeReaderReadData(tmp, &bytes, limereader);
  limeDestroyReader(limereader);
  fclose(ifs);
  if(strcmp("DiracFermion_Sink", tmp)) return(0);
  else if(strcmp("DiracFermion_Source_Sink_Pairs", tmp)) return(1);
  else if(strcmp("DiracFermion_ScalarSource_TwelveSink", tmp)) return(2);
  else if(strcmp("DiracFermion_ScalarSource_FourSink", tmp)) return(3);
  else return(-1);
}

int get_source_type(char * filename) {
  FILE * ifs;
  int status=0;
  n_uint64_t bytes;
  char * header_type;
  char tmp[500];
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
    if(strcmp("source-type",header_type) == 0) break;
  }
  if(status == LIME_EOF) {
    if(g_proc_id == 0) {
      fprintf(stderr, "no source-type record found in file %s\n",filename);
    }
    limeDestroyReader(limereader);
    fclose(ifs);
    return(-1);
  }
  bytes = limeReaderBytes(limereader);
  status = limeReaderReadData(tmp, &bytes, limereader);
  limeDestroyReader(limereader);
  fclose(ifs);
  if(strcmp("DiracFermion_Source", tmp)) return(0);
  else if(strcmp("DiracFermion_ScalarSource", tmp)) return(1);
  else if(strcmp("DiracFermion_FourScalarSource", tmp)) return(2);
  else if(strcmp("DiracFermion_TwelveScalarSource", tmp)) return(3);
  else return(-1);
}

 
int read_lime_spinor(spinor * const s, spinor * const r, char * filename, const int position) {
  FILE * ifs;
  int status=0, getpos=-1;
  n_uint64_t bytes;
  char * header_type;
  LimeReader * limereader;
  int prec = 32;
  DML_Checksum checksum;
  
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
    if(!strcmp("scidac-binary-data",header_type)) getpos++;
    if(getpos == position) break;
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
  if((int)bytes == LX*g_nproc_x*LY*g_nproc_y*LZ*g_nproc_z*T*g_nproc_t*sizeof(spinor)) prec = 64;
  else if((int)bytes == LX*g_nproc_x*LY*g_nproc_y*LZ*g_nproc_z*T*g_nproc_t*sizeof(spinor)/2) prec = 32;
  else {
    if(g_proc_id == 0) {
      fprintf(stderr, "wrong length in eospinor: bytes = %d, not %d. Aborting read!\n", (int)bytes, LX*g_nproc_x*LY*g_nproc_y*LZ*g_nproc_z*T*g_nproc_t*sizeof(spinor)/2);
    }
    return(-1);
  }
  if(g_proc_id == 0 && g_debug_level > 2) {
    printf("# %d Bit precision read\n", prec);
  }

  status = read_binary_spinor_data(s, r, limereader, prec, &checksum);

  if(g_proc_id == 0) {
    printf("# checksum for DiracFermion field in file %s position %d is %#x %#x\n", 
	   filename, position, checksum.suma, checksum.sumb);
  }

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

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
# include<mpi.h>
# include<unistd.h> 
#endif
#include<math.h>
#include"global.h"
#include"su3.h"
#include"lime.h" 
#include"read_input.h"
#include"io_utils.h"
#include"dml.h"
#include"io.h"
#include"gauge_io.h"

/* #define MAXBUF 1048576 */


int write_binary_gauge_data(LimeWriter * limewriter,
			    const int prec, DML_Checksum * ans) {
  
  int x, X, y, Y, z, Z, tt, t0, tag=0, id=0, status=0;
  su3 tmp[4];
  su3 tmp3[4];
  float tmp2[72];
  int coords[4];
/*   n_uint64_t bytes; */
  n_uint64_t bytes; 
  DML_SiteRank rank;
#ifdef MPI
  MPI_Status mpi_status;
#endif

  DML_checksum_init(ans);

  if(prec == 32) bytes = (n_uint64_t)2*sizeof(su3);
  else bytes = (n_uint64_t)4*sizeof(su3);
  for(t0 = 0; t0 < T*g_nproc_t; t0++) {
    tt = t0 - g_proc_coords[0]*T;
    coords[0] = t0 / T;
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
#ifdef MPI
	  MPI_Cart_rank(g_cart_grid, coords, &id);
#endif
	  if(g_cart_id == id) {
	    rank = (DML_SiteRank) (((t0*LZ*g_nproc_z + z)*LY*g_nproc_y + y)*LX*g_nproc_x + x);
	  }
	  if(g_cart_id == 0) {
	    if(g_cart_id == id) {
	      memcpy(&tmp3[0], &g_gauge_field[ g_ipt[tt][X][Y][Z] ][1], sizeof(su3));
	      memcpy(&tmp3[1], &g_gauge_field[ g_ipt[tt][X][Y][Z] ][2], sizeof(su3));
	      memcpy(&tmp3[2], &g_gauge_field[ g_ipt[tt][X][Y][Z] ][3], sizeof(su3));
	      memcpy(&tmp3[3], &g_gauge_field[ g_ipt[tt][X][Y][Z] ][0], sizeof(su3));

#ifndef WORDS_BIGENDIAN
	      if(prec == 32) {
		byte_swap_assign_double2single(tmp2, tmp3, 4*sizeof(su3)/8); 
		DML_checksum_accum(ans, rank, (char*) tmp2, 4*sizeof(su3)/2);
		status = limeWriteRecordData((void*)&tmp2, &bytes, limewriter);
	      }
	      else {
		byte_swap_assign(tmp, tmp3, 4*sizeof(su3)/8);
		DML_checksum_accum(ans, rank, (char*) tmp, 4*sizeof(su3));
		status = limeWriteRecordData((void*)&tmp, &bytes, limewriter);
	      }
#else
	      if(prec == 32) {
		double2single(tmp2, tmp3, 4*sizeof(su3)/8); 
		DML_checksum_accum(ans, rank, (char*) tmp2, 4*sizeof(su3)/2);
		status = limeWriteRecordData((void*)&tmp2, &bytes, limewriter);
	      }
	      else {
		DML_checksum_accum(ans, rank, (char*) tmp3, 4*sizeof(su3));
		status = limeWriteRecordData((void*)&tmp3, &bytes, limewriter);
	      }
#endif
	    }
#ifdef MPI
	    else {
	      if(prec == 32) {
		MPI_Recv(tmp2, 4*sizeof(su3)/8, MPI_FLOAT, id, tag, g_cart_grid, &mpi_status);
		DML_checksum_accum(ans, rank, (char*) tmp2, 4*sizeof(su3)/2);
		status = limeWriteRecordData((void*)&tmp2, &bytes, limewriter);
	      }
	      else {
		MPI_Recv(tmp, 4*sizeof(su3)/8, MPI_DOUBLE, id, tag, g_cart_grid, &mpi_status);
		DML_checksum_accum(ans, rank, (char*) tmp, 4*sizeof(su3));
		status = limeWriteRecordData((void*)&tmp, &bytes, limewriter);
	      }
	    }
#endif
	    if(status < 0 ) {
	      fprintf(stderr, "LIME write error %d\n", status);
	      fprintf(stderr, "x %d, y %d, z %d, t %d (%d,%d,%d,%d)\n",x,y,z,tt,X,Y,Z,tt);
	      fprintf(stderr, "id = %d, bytes = %lu, size = %d\n", g_proc_id, bytes,  4*sizeof(su3)/8); 
#ifdef MPI
	      MPI_Abort(MPI_COMM_WORLD, 1);
	      MPI_Finalize();
#endif
	      exit(500);
	    }
	  }
#ifdef MPI
	  else {
	    if(g_cart_id == id){
	      memcpy(&tmp3[0], &g_gauge_field[ g_ipt[tt][X][Y][Z] ][1], sizeof(su3));
	      memcpy(&tmp3[1], &g_gauge_field[ g_ipt[tt][X][Y][Z] ][2], sizeof(su3));
	      memcpy(&tmp3[2], &g_gauge_field[ g_ipt[tt][X][Y][Z] ][3], sizeof(su3));
	      memcpy(&tmp3[3], &g_gauge_field[ g_ipt[tt][X][Y][Z] ][0], sizeof(su3));
#  ifndef WORDS_BIGENDIAN
	      if(prec == 32) {
		byte_swap_assign_double2single(tmp2, tmp3, 4*sizeof(su3)/8);
		MPI_Send((void*) tmp2, 4*sizeof(su3)/8, MPI_FLOAT, 0, tag, g_cart_grid);
	      }
	      else {
		byte_swap_assign(tmp, tmp3, 4*sizeof(su3)/8);
		MPI_Send((void*) tmp, 4*sizeof(su3)/8, MPI_DOUBLE, 0, tag, g_cart_grid);
	      }
#  else
	      if(prec == 32) {
		double2single(tmp2, tmp3, 4*sizeof(su3)/8);
		MPI_Send((void*) tmp2, 4*sizeof(su3)/8, MPI_FLOAT, 0, tag, g_cart_grid);
	      }
	      else {
		MPI_Send((void*) tmp3, 4*sizeof(su3)/8, MPI_DOUBLE, 0, tag, g_cart_grid);
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
      }
    }
  }
  return(0);
}


int read_binary_gauge_data(LimeReader * limereader, 
			   const double prec, DML_Checksum * ans) {
  int t, x, y , z, status=0;
/*   n_uint64_t bytes; */
  n_uint64_t bytes;
  su3 tmp[4];
  float tmp2[72];
  DML_SiteRank rank;

  DML_checksum_init(ans);

  if(prec == 32) bytes = (n_uint64_t) 4*sizeof(su3)/2;
  else bytes = (n_uint64_t) 4*sizeof(su3);
  for(t = 0; t < T; t++){
    for(z = 0; z < LZ; z++){
      for(y = 0; y < LY; y++){
#if (defined MPI)
	limeReaderSeek(limereader,(n_uint64_t) 
		       (((n_uint64_t) g_proc_coords[1]*LX) + 
			((n_uint64_t) (((g_proc_coords[0]*T+t)*g_nproc_z*LZ+g_proc_coords[3]*LZ+z)*g_nproc_y*LY 
			 + g_proc_coords[2]*LY+y)*LX*g_nproc_x))*bytes,
		       SEEK_SET);
#endif
	for(x = 0; x < LX; x++) {
	  rank = (DML_SiteRank) (g_proc_coords[1]*LX + 
				 (((g_proc_coords[0]*T+t)*g_nproc_z*LZ+g_proc_coords[3]*LZ+z)*g_nproc_y*LY 
				  + g_proc_coords[2]*LY+y)*((DML_SiteRank)LX*g_nproc_x) + x);
	  if(prec == 32) {
	    status = limeReaderReadData(tmp2, &bytes, limereader);
 	    DML_checksum_accum(ans, rank, (char *) tmp2, bytes);
	  }
	  else {
	    status = limeReaderReadData(tmp, &bytes, limereader);
	    DML_checksum_accum(ans, rank, (char *) tmp, bytes);
	  }
	  if(status < 0 && status != LIME_EOR) {
	    fprintf(stderr, "LIME read error occured with status = %d while reading!\n Aborting...\n", status);
#ifdef MPI
	    MPI_Abort(MPI_COMM_WORLD, 1);
	    MPI_Finalize();
#endif
	    exit(500);
	  }
#ifndef WORDS_BIGENDIAN
	  if(prec == 32) {
	    byte_swap_assign_single2double(&g_gauge_field[ g_ipt[t][x][y][z] ][0], &tmp2[3*18], sizeof(su3)/8);
	    byte_swap_assign_single2double(&g_gauge_field[ g_ipt[t][x][y][z] ][1], &tmp2[0*18], sizeof(su3)/8);
	    byte_swap_assign_single2double(&g_gauge_field[ g_ipt[t][x][y][z] ][2], &tmp2[1*18], sizeof(su3)/8);
	    byte_swap_assign_single2double(&g_gauge_field[ g_ipt[t][x][y][z] ][3], &tmp2[2*18], sizeof(su3)/8);
	  }
	  else {
	    byte_swap_assign(&g_gauge_field[ g_ipt[t][x][y][z] ][0], &tmp[3], sizeof(su3)/8);
	    byte_swap_assign(&g_gauge_field[ g_ipt[t][x][y][z] ][1], &tmp[0], sizeof(su3)/8);
	    byte_swap_assign(&g_gauge_field[ g_ipt[t][x][y][z] ][2], &tmp[1], sizeof(su3)/8);
	    byte_swap_assign(&g_gauge_field[ g_ipt[t][x][y][z] ][3], &tmp[2], sizeof(su3)/8);
	  }
#else
	  if(prec == 32) {
	    single2double(&g_gauge_field[ g_ipt[t][x][y][z] ][0], &tmp2[3*18], sizeof(su3)/8);
	    single2double(&g_gauge_field[ g_ipt[t][x][y][z] ][1], &tmp2[0], sizeof(su3)/8);
	    single2double(&g_gauge_field[ g_ipt[t][x][y][z] ][2], &tmp2[18], sizeof(su3)/8);
	    single2double(&g_gauge_field[ g_ipt[t][x][y][z] ][3], &tmp2[2*18], sizeof(su3)/8);
	  }
	  else {
	    memcpy(&g_gauge_field[ g_ipt[t][x][y][z] ][0], &tmp[3], sizeof(su3));
	    memcpy(&g_gauge_field[ g_ipt[t][x][y][z] ][1], &tmp[0], sizeof(su3));
	    memcpy(&g_gauge_field[ g_ipt[t][x][y][z] ][2], &tmp[1], sizeof(su3));
	    memcpy(&g_gauge_field[ g_ipt[t][x][y][z] ][3], &tmp[2], sizeof(su3));
	  }
#endif
	}
      }
    }
  }
#ifdef MPI
  DML_checksum_combine(ans);
#endif
  return(0);
}


int write_lime_gauge_field(char * filename, const double plaq, const int counter, const int prec) {
  FILE * ofs = NULL;
  LimeWriter * limewriter = NULL;
  LimeRecordHeader * limeheader = NULL;
  /* Message end and Message begin flag */
  int ME_flag=0, MB_flag=0, status=0;
  n_uint64_t bytes;
  DML_Checksum checksum;
  
  write_xlf_info(plaq, counter, filename, 0);

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
    write_ildg_format_xml("temp.xml", limewriter, 0);
    
    bytes = ((n_uint64_t)LX*g_nproc_x)*((n_uint64_t)LY*g_nproc_y)*((n_uint64_t)LZ*g_nproc_z)*((n_uint64_t)T*g_nproc_t)*((n_uint64_t)4*sizeof(su3));
    if(prec == 32) bytes = bytes/((n_uint64_t)2);
    MB_flag=0; ME_flag=0;
    limeheader = limeCreateHeader(MB_flag, ME_flag, "ildg-binary-data", bytes);
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
  }

  write_binary_gauge_data(limewriter, prec, &checksum);
  if(g_proc_id == 0) {
    printf("# checksum for Gauge field written to file %s is %#x %#x\n", 
	   filename, checksum.suma, checksum.sumb);
  }

  if(g_cart_id == 0) {
    limeDestroyWriter( limewriter );
    fflush(ofs);
    fclose(ofs);
  }
  write_checksum(filename, &checksum);

  return(0);
}

int read_lime_gauge_field(char * filename) {
  FILE * ifs;
  int status, prec = 64;
  n_uint64_t bytes;
  char * header_type;
  LimeReader * limereader;
  DML_Checksum checksum;

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
    if(strcmp("ildg-binary-data",header_type) == 0) break;
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

  if(bytes == ((n_uint64_t)LX*g_nproc_x)*((n_uint64_t)LY*g_nproc_y)*((n_uint64_t)LZ*g_nproc_z)*((n_uint64_t)T*g_nproc_t)*((n_uint64_t)4*sizeof(su3))) prec = 64;
  else if(bytes == ((n_uint64_t)LX*g_nproc_x)*((n_uint64_t)LY*g_nproc_y)*((n_uint64_t)LZ*g_nproc_z)*((n_uint64_t)T*g_nproc_t)*((n_uint64_t)4*sizeof(su3)/2)) prec = 32;
  else {
    fprintf(stderr, "Probably wrong lattice size or precision (bytes=%lu) in file %s\n", bytes, filename);
    fprintf(stderr, "Aborting...!\n");
    fflush( stdout );
#ifdef MPI
    MPI_Abort(MPI_COMM_WORLD, 1);
    MPI_Finalize();
#endif
    exit(501);
  }
  if(g_proc_id == 0 && g_debug_level > 2) {
    printf("# %d Bit precision read\n", prec);
  }
  read_binary_gauge_data(limereader, prec, &checksum);

  limeDestroyReader(limereader);
  fclose(ifs);
  if(g_proc_id == 0) {
    printf("# checksum for gaugefield %s is %#x %#x\n", filename, checksum.suma, checksum.sumb);
  }
  return(0);
}

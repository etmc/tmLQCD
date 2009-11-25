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
#include"propagator_io.h"
#include<io/dml.h>
#include"io.h"
#include "io/spinor.h"

#include "propagator_io.h"


int read_source(spinor * const s, spinor * const r, char *filename,
                const int format, const int position) {
  int err = 0;

  if(format == 11) {
    /* cmi format */
    err = read_spinorfield_cm_single(s, r, filename,  -1, 1);
  }
  else if(format == 10) {
    /* GWC format */
    err = read_spinorfield_eo_time(s, r, filename);
  }
  else
  {
    /* ETMC standard format */
    read_spinor(g_spinor_field[0], g_spinor_field[1], filename, position);
  }

  if(err != 0) {
    if(g_proc_id == 0) {
      printf("Error reading source! Aborting...\n");
    }
#ifdef MPI
    MPI_Abort(MPI_COMM_WORLD, 1);
    MPI_Finalize();
#endif
    exit(-1);
  }
  return(0);
}



int get_propagator_type(char * filename) {
  FILE * ifs;
  int ret=-1, status=0;
  n_uint64_t bytes;
  char * tmp;
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
    if(strcmp("propagator-type", limeReaderType(limereader)) == 0) break;
  }

  tmp = (char*)calloc(500, sizeof(char));
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

  if(strcmp("DiracFermion_Sink", tmp) == 0) ret = 0;
  else if(strcmp("DiracFermion_Source_Sink_Pairs", tmp) == 0) ret =1;
  else if(strcmp("DiracFermion_ScalarSource_TwelveSink", tmp) == 0) ret = 2;
  else if(strcmp("DiracFermion_ScalarSource_FourSink", tmp) == 0) ret = 3;
  free(tmp);
  return(ret);
}

int get_source_type(char * filename) {
  FILE * ifs;
  int ret=-1, status=0;
  n_uint64_t bytes;
  char * tmp;
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
    if(strcmp("source-type", limeReaderType(limereader)) == 0) break;
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
  tmp = (char*)calloc(500, sizeof(char));
  status = limeReaderReadData(tmp, &bytes, limereader);
  limeDestroyReader(limereader);
  fclose(ifs);
  if(strcmp("DiracFermion_Source", tmp) == 0) ret =0;
  else if(strcmp("DiracFermion_ScalarSource", tmp) == 0) ret =1;
  else if(strcmp("DiracFermion_FourScalarSource", tmp) == 0) ret =2;
  else if(strcmp("DiracFermion_TwelveScalarSource", tmp) == 0) ret =3;
  free(tmp);
  return(ret);
}

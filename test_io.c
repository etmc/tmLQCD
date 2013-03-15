/***********************************************************************
 *
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008 Carsten Urbach
 *               2013                               Bartosz Kostrzewa
 *
 * This file is part of tmLQCD.
 *
 * tmLQCD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * tmLQCD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with tmLQCD.  If not, see <http://www.gnu.org/licenses/>.
 *
 *******************************************************************************/

#define MAIN_PROGRAM
#include "lime.h"
#if HAVE_CONFIG_H
#include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include <signal.h>
#ifdef MPI
# include <mpi.h>
#endif
#ifdef OMP
# include <omp.h>
#endif
#include "global.h"
#include "git_hash.h"
#include <io/params.h>
#include <io/gauge.h>
#include <io/dml.h>
#include "getopt.h"
#include "ranlxd.h"
#include "geometry_eo.h"
#include "start.h"
#include "measure_gauge_action.h"
#include "measure_rectangles.h"
#ifdef MPI
# include "xchange/xchange.h"
#endif
#include "read_input.h"
#include "mpi_init.h"
#include "sighandler.h"
#include "init/init.h"
#include "test/check_geometry.h"

#include "buffers/gauge.h"
#include "buffers/utils.h"

#include "dirty_shameful_business.h"


int const rlxdsize = 105;

static void usage();
static void process_args(int argc, char *argv[], char ** input_filename, char ** filename);
static void set_default_filenames(char ** input_filename, char ** filename);

enum enum_failure_enum
{
  FAIL_READ,
  FAIL_READ_CHKSUM,
  FAIL_REREAD,
  FAIL_REREAD_CHKSUM,
  FAIL_COMPARE_CHKSUM,
  FAIL_WRITE
};
   
char *failure_names[6] = { 
  "read\0",
  "read checksum\0", 
  "reread \0", 
  "reread checksum\0",
  "compare checksum\0", 
  "write\0" };

typedef enum enum_failure_enum enum_failure_t;

typedef struct 
{
  char filename_orig[200];
  char filename_copy[200];
  DML_Checksum checksum_orig;
  DML_Checksum checksum_copy;
} test_conf_t;

typedef struct
{
  int fail_iteration;
  int fail_sub_iteration;
  enum_failure_t fail_type;
} failure_t;

typedef struct
{
  int length;
  failure_t* ptr;
} failure_flex_array_t;

static void add_failure(failure_flex_array_t*, const enum_failure_t, const int iteration, const int sub_iteration);
static void output_failures(const failure_flex_array_t* const);
static void shuffle(int* array,size_t n);


#define ITERATIONS 1
#define NUM_TESTCONFS 50
#define NUM_READS 20
#define NUM_REREADS 20
#define MIN_DELAY 0
#define MAX_DELAY 5

int main(int argc,char *argv[]) {

  char *filename = NULL;
  char gauge_filename[50];
  char tmp_filename[50];
  char *input_filename = NULL;
  int status = 0;
  
  failure_flex_array_t failures;
  failures.ptr = NULL;
  failures.length = 0;
  
  test_conf_t test_confs[NUM_TESTCONFS];
  gauge_field_t gauge_fields[NUM_TESTCONFS];
  srand(time(NULL));
  int conf_indices[NUM_TESTCONFS];

  char* testconf_filename_base = "conf";

  /* Energy corresponding to the Gauge part */
  double plaquette_energy = 0.;

  paramsXlfInfo *xlfInfo;

#if (defined SSE || defined SSE2 || SSE3)
  signal(SIGILL,&catch_ill_inst);
#endif

  strcpy(gauge_filename,"conf.save");
  strcpy(tmp_filename, ".conf.tmp");

  verbose = 1;
  g_use_clover_flag = 0;

#ifdef MPI

#  ifdef OMP
  int mpi_thread_provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_SERIALIZED, &mpi_thread_provided);
#  else
  MPI_Init(&argc, &argv);
#  endif

  MPI_Comm_rank(MPI_COMM_WORLD, &g_proc_id);
#else
  g_proc_id = 0;
#endif

  process_args(argc,argv,&input_filename,&filename);
  set_default_filenames(&input_filename,&filename);

  /* Read the input file */
  if( (status = read_input(input_filename)) != 0) {
    fprintf(stderr, "Could not find input file: %s\nAborting...\n", input_filename);
    exit(-1);
  }

#ifdef OMP
  init_openmp();
#endif

  tmlqcd_mpi_init(argc, argv);
 
  initialize_gauge_buffers(NUM_TESTCONFS+1); 
  
  for( int i = 0; i < NUM_TESTCONFS; ++i) {
    snprintf(test_confs[i].filename_orig,200,"%s.%04d",testconf_filename_base,i+10);
    snprintf(test_confs[i].filename_copy,200,"%s.%04d.copy",testconf_filename_base,i+10);
    DML_checksum_init(&test_confs[i].checksum_orig);
    DML_checksum_init(&test_confs[i].checksum_copy);
    gauge_fields[i] = get_gauge_field();
    conf_indices[i] = i;
  }
  
#ifndef MPI
  g_dbw2rand = 0;
#endif
  
  g_mu = g_mu1;
  
#ifdef _GAUGE_COPY
  status = init_gauge_field(VOLUMEPLUSRAND + g_dbw2rand, 1);
#else
  status = init_gauge_field(VOLUMEPLUSRAND + g_dbw2rand, 0);
#endif
  if (status != 0) {
    fprintf(stderr, "Not enough memory for gauge_fields! Aborting...\n");
    exit(0);
  }
  status  = init_geometry_indices(VOLUMEPLUSRAND + g_dbw2rand);
  if (status != 0) {
    fprintf(stderr, "Not enough memory for geometry_indices! Aborting...\n");
    exit(0);
  }

  /* define the geometry */
  geometry();

  status = check_geometry();

  if (status != 0) {
    fprintf(stderr, "Checking of geometry failed. Unable to proceed.\nAborting....\n");
    exit(1);
  }

  /* Initialise random number generator */
  start_ranlux(rlxd_level, random_seed );

  /*For parallelization: exchange the gaugefield (not really necessary in this test...) */
#ifdef MPI
  exchange_gauge_field(&g_gf);
  //xchange_gauge(&g_gf);
#endif

  /* Loop for tests */
  for(int j = 0; j < ITERATIONS; j++) {
    if(g_proc_id == 0) {
      printf("#\n# Starting test iteration %d\n", j);
    }

    for(int confnum = 0; confnum < NUM_TESTCONFS; ++confnum) {
      ohnohack_remap_g_gauge_field(gauge_fields[confnum]);
      if( g_proc_id == 0 )
        printf("\nReading gauge field %s. Iteration %d\n",test_confs[confnum].filename_orig,j);
      if( (status = read_gauge_field_checksum( test_confs[confnum].filename_orig, 
                                               &test_confs[confnum].checksum_orig)) != 0) {
        if( g_proc_id == 0 )      
          fprintf(stdout, "Error %d while reading gauge field from %s\n", status, test_confs[confnum].filename_orig);
        add_failure(&failures,FAIL_READ,j,-1);
      }
    }

    int delay = MIN_DELAY; 
    for(int num_rereads = 0; num_rereads < NUM_REREADS; ++num_rereads) {

      /* write copies in random order */
      shuffle(conf_indices,NUM_TESTCONFS);
      for(int i = 0; i < NUM_TESTCONFS; ++i){
        int confnum = conf_indices[i];
        ohnohack_remap_g_gauge_field(gauge_fields[confnum]);
        plaquette_energy = measure_gauge_action(gauge_fields[confnum]);
        xlfInfo = construct_paramsXlfInfo(plaquette_energy/(6.*VOLUME*g_nproc), num_rereads);
        if (g_proc_id == 0) {
          fprintf(stdout, "\n# Writing gauge field to %s. Iteration %d, reread %d\n", test_confs[confnum].filename_copy,j,num_rereads);
        }
        if((status = write_gauge_field(test_confs[confnum].filename_copy, gauge_precision_write_flag, xlfInfo) != 0 )) {
          /* Writing the gauge field failed directly */
          if(g_proc_id==0)
            fprintf(stdout, "Error %d while writing gauge field to %s\n", status, test_confs[confnum].filename_copy);
          add_failure(&failures,FAIL_WRITE,j,num_rereads);
        } else {
          if (g_proc_id == 0) {
            fprintf(stdout, "# Write completed.\n");
          }
        }
        free(xlfInfo);
      }

      if( delay > 0) {
#ifdef MPI
        MPI_Barrier(MPI_COMM_WORLD);
#endif
        sleep(delay);
#ifdef MPI
        MPI_Barrier(MPI_COMM_WORLD);
#endif
      }
                  
      if(delay < MAX_DELAY)
        ++delay;
      else
        delay = 0;

      /* now reread in random order */
      shuffle(conf_indices,NUM_TESTCONFS);
      for(int i = 0; i < NUM_TESTCONFS; ++i) {
        int confnum = conf_indices[i];
        if(g_proc_id == 0)
          printf("#\n\n RANDOM reread test %d, iteration %d\n",num_rereads,j);
        ohnohack_remap_g_gauge_field(gauge_fields[confnum]);
        ohnohack_remap_g_gauge_field(gauge_fields[confnum]);
        plaquette_energy = measure_gauge_action(gauge_fields[confnum]);
        xlfInfo = construct_paramsXlfInfo(plaquette_energy/(6.*VOLUME*g_nproc), num_rereads);

        if( (status = read_gauge_field_checksum(test_confs[confnum].filename_copy,&test_confs[confnum].checksum_copy)) != 0) {
              if( g_proc_id == 0 )
                fprintf(stdout, "WARNING, verification of %s discovered errors.\n", test_confs[confnum].filename_copy);
            add_failure(&failures,FAIL_REREAD_CHKSUM,j,num_rereads);
        } else {
          if (g_proc_id == 0)
            fprintf(stdout, "# Write successfully verified.\n");
        } 
        if( test_confs[confnum].checksum_orig.suma != test_confs[confnum].checksum_copy.suma ||
            test_confs[confnum].checksum_orig.sumb != test_confs[confnum].checksum_copy.sumb ) {
            if( g_proc_id == 0 ) 
              printf("# Write verification successful but new checksum does not match the checksum originally computed!\n");
            add_failure(&failures,FAIL_COMPARE_CHKSUM,j,num_rereads);
        }
        free(xlfInfo);
      }
    }
  } /* end of loop over test iterations */
  
  /* TEST add_failure 
  add_failure(&failures,FAIL_READ_CHKSUM,0,22);
  add_failure(&failures,FAIL_READ_PLAQ,20,11);
  add_failure(&failures,FAIL_REREAD_CHKSUM,55,33);
  add_failure(&failures,FAIL_REREAD_PLAQ,121,44);
  add_failure(&failures,FAIL_WRITE,77,55); // */
  
  output_failures(&failures);

  for( int i = 0; i < NUM_TESTCONFS; ++i){
    return_gauge_field( &gauge_fields[i] );
  }

#ifdef MPI
  MPI_Finalize();
#endif
#ifdef OMP
  free_omp_accumulators();
#endif
  free_gauge_tmp();
  free_gauge_field();
  free_geometry_indices();
  return_gauge_field(&g_gf);
  finalize_gauge_buffers();
  free(input_filename);
  free(filename);
  free(failures.ptr);
  return(0);
#ifdef _KOJAK_INST
#pragma pomp inst end(main)
#endif
}

static void usage(){
  fprintf(stdout, "IO test for LIME and LEMON configuration reading, writing and rereading\n");
  fprintf(stdout, "Version %s %s \n\n", PACKAGE_VERSION, git_hash);
  fprintf(stdout, "Please send bug reports to %s\n", PACKAGE_BUGREPORT);
  fprintf(stdout, "Usage:   test_io [options]\n");
  fprintf(stdout, "Options: [-f input-filename]  default: hmc.input\n");
  fprintf(stdout, "         [-v] more verbosity\n");
  fprintf(stdout, "         [-V] print version information and exit\n");
  fprintf(stdout, "         [-h|-? this help]\n");
  exit(0);
}

static void process_args(int argc, char *argv[], char ** input_filename, char ** filename) {
  int c;
  while ((c = getopt(argc, argv, "h?vVf:o:")) != -1) {
    switch (c) {
      case 'f':
        *input_filename = calloc(200, sizeof(char));
        strncpy(*input_filename, optarg, 200);
        break;
      case 'o':
        *filename = calloc(200, sizeof(char));
        strncpy(*filename, optarg, 200);
        break;
      case 'v':
        verbose = 1;
        break;
      case 'V':
        if(g_proc_id == 0) {
          fprintf(stdout,"%s %s\n",PACKAGE_STRING,git_hash);
        }
        exit(0);
        break;
      case 'h':
      case '?':
      default:
        if( g_proc_id == 0 ) {
          usage();
        }
        break;
    }
  }
}

static void set_default_filenames(char ** input_filename, char ** filename) {
  if( *input_filename == NULL ) {
    *input_filename = calloc(14, sizeof(char));
    strcpy(*input_filename,"test_io.input");
  }
  
  if( *filename == NULL ) {
    *filename = calloc(7, sizeof(char));
    strcpy(*filename,"output");
  } 
}

static void add_failure(failure_flex_array_t* failures, const enum_failure_t fail_type, const int iteration, const int sub_iteration) {
  // expand the flexible array
  ++failures->length;
  failures->ptr = realloc(failures->ptr, failures->length * sizeof(failure_t));
  
  failure_t* element = &failures->ptr[failures->length-1];
  
  element->fail_iteration = iteration;
  element->fail_sub_iteration = sub_iteration;
  element->fail_type = fail_type;
}

static void output_failures(const failure_flex_array_t* const failures) {
  if( g_proc_id == 0 ) {
    if( failures->length > 0 ) {

      printf("Failures:\n");
      for(int i = 0; i < failures->length; ++i) {
        printf("%s at iteration %d, sub iteration %d\n", failure_names[ failures->ptr[i].fail_type ], 
                                                         failures->ptr[i].fail_iteration,
                                                         failures->ptr[i].fail_sub_iteration );
      }
    } else {
      printf("No failures!\n");
    }
  }
}

static void shuffle(int *array, size_t n)
{
  if (n > 1) 
  {
    size_t i;
    for (i = 0; i < n - 1; i++) 
    {
      size_t j = i + rand() / (RAND_MAX / (n - i) + 1);
      int t = array[j];
      array[j] = array[i];
      array[i] = t;
    }
  }
}


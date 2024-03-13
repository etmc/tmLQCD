/***********************************************************************
 * Copyright (C) 2018 Bartosz Kostrzewa
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
 ***********************************************************************/

#ifndef OMP_ACCUMULATOR_H
#define OMP_ACCUMULATOR_H

#ifdef HAVE_CONFIG_H
# include<tmlqcd_config.h>
#endif
#ifdef TM_USE_OMP
#include <omp.h>
#endif
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <assert.h>
#include "global.h"

#define TM_CHECK_INIT_OMP_ACC(x) assert( (x)->init == 1 );

#define TM_CHECK_BOUNDS_OMP_ACC(x, count) assert( (x)->num_values == count );

#ifdef TM_USE_OMP
#define TM_CHECK_NUM_THREADS_OMP_ACC(x) assert( (x)->num_threads == omp_get_num_threads() );
#else
#define TM_CHECK_NUM_THREADS_OMP_ACC(x) assert( (x)->num_threads == 1 ); 
#endif

#ifdef TM_USE_OMP
#define TM_CHECK_SINGLE_THREAD assert( omp_in_parallel() != 1 );
#else
#define TM_CHECK_SINGLE_THREAD
#endif

typedef struct omp_re_acc_t {
  double * mem;
  double ** acc;

  int num_threads;
  int num_values;

  int init;

} omp_re_acc_t;

typedef struct omp_cplx_acc_t {
  complex double * mem;
  complex double ** acc;

  int num_threads;
  int num_values;

  int init;

} omp_cplx_acc_t;

static inline omp_re_acc_t new_omp_re_acc(void) {
  omp_re_acc_t ret;
  ret.mem = (double*) NULL;
  ret.acc = (double**) NULL;
  ret.init = 0;
  ret.num_threads = 0;
  ret.num_values = 0;
  return(ret);
}

static inline omp_cplx_acc_t new_omp_cplx_acc(void) {
  omp_cplx_acc_t ret;
  ret.mem = (complex double*) NULL;
  ret.acc = (complex double**) NULL;
  ret.init = 0;
  ret.num_threads = 0;
  ret.num_values = 0;
  return(ret);
}

static inline void omp_re_acc_free(omp_re_acc_t * const acc ){
  free( acc->mem );
  free( acc->acc );
  acc->init = 0;
}

static inline void omp_cplx_acc_free(omp_cplx_acc_t * const acc ){
  free( acc->mem );
  free( acc->acc );
  acc->init = 0;
}

static inline void omp_re_acc_reset(omp_re_acc_t * const acc){
  TM_CHECK_SINGLE_THREAD;
  TM_CHECK_INIT_OMP_ACC(acc);
  memset(acc->mem, 0, acc->num_threads*acc->num_values*sizeof(double)); 
}

static inline void omp_cplx_acc_reset(omp_cplx_acc_t * const acc){
  TM_CHECK_SINGLE_THREAD;
  TM_CHECK_INIT_OMP_ACC(acc);
  memset(acc->mem, 0, acc->num_threads*acc->num_values*sizeof(complex double)); 
}

static inline void omp_re_acc_init(omp_re_acc_t * const acc, const int num_values) {
  TM_CHECK_SINGLE_THREAD;
  if( acc->init != 0 ){
    omp_re_acc_free(acc);
  }
#ifdef TM_USE_OMP
  acc->num_threads = omp_num_threads;
#else
  acc->num_threads = 1;
#endif
  acc->num_values = num_values;
  acc->acc = (double**) calloc( acc->num_threads, sizeof( double* ) );
  acc->mem = (double*) calloc( acc->num_threads*acc->num_values, sizeof(double) );
  for( int thread = 0; thread < acc->num_threads; ++thread ){
    acc->acc[thread] = acc->mem + thread*acc->num_values;
  }
  acc->init = 1;
  omp_re_acc_reset(acc);
}

static inline void omp_cplx_acc_init(omp_cplx_acc_t * const acc, const int num_values) {
  TM_CHECK_SINGLE_THREAD;
  if( acc->init != 0 ){
    omp_cplx_acc_free(acc);
  }
#ifdef TM_USE_OMP
  acc->num_threads = omp_num_threads;
#else
  acc->num_threads = 1;
#endif
  acc->num_values = num_values;
  acc->acc = (complex double**) calloc( acc->num_threads, sizeof( complex double* ) );
  acc->mem = (complex double*) calloc( acc->num_threads*acc->num_values, sizeof(complex double) );
  for( int thread = 0; thread < acc->num_threads; ++thread ){
    acc->acc[thread] = acc->mem + thread*acc->num_values;
  }
  acc->init = 1;
  omp_cplx_acc_reset(acc);
}

static inline void omp_re_acc_add(omp_re_acc_t * const acc, const double * const values, const int num_values){
  TM_CHECK_INIT_OMP_ACC(acc);
  TM_CHECK_BOUNDS_OMP_ACC(acc, num_values);
  TM_CHECK_NUM_THREADS_OMP_ACC(acc);
  for(int i = 0; i < num_values; ++i){ 
#ifdef TM_USE_OMP
    acc->acc[omp_get_thread_num()][i] += values[i];
#else
    acc->acc[0][i] += values[i];
#endif
  }
}

static inline void omp_cplx_acc_add(omp_cplx_acc_t * const acc, const complex double * const values, const int num_values){
  TM_CHECK_INIT_OMP_ACC(acc);
  TM_CHECK_BOUNDS_OMP_ACC(acc, num_values);
  TM_CHECK_NUM_THREADS_OMP_ACC(acc);
  for(int i = 0; i < num_values; ++i){ 
#ifdef TM_USE_OMP
    acc->acc[omp_get_thread_num()][i] += values[i];
#else
    acc->acc[0][i] += values[i];
#endif
  }
}

static inline void omp_re_acc_reduce(double * const ret, omp_re_acc_t * const acc){
  TM_CHECK_SINGLE_THREAD;
  TM_CHECK_INIT_OMP_ACC(acc);
  for(int thread = 0; thread < acc->num_threads; thread++){
    for(int value = 0; value < acc->num_values; value++){
      ret[value] += acc->acc[thread][value];
    }
  }
}

static inline void omp_cplx_acc_reduce(complex double * const ret, omp_cplx_acc_t * const acc){
  TM_CHECK_SINGLE_THREAD;
  TM_CHECK_INIT_OMP_ACC(acc);
  for(int thread = 0; thread < acc->num_threads; thread++){
    for(int value = 0; value < acc->num_values; value++){
      ret[value] += acc->acc[thread][value];
    }
  }
}


#endif  // header guard

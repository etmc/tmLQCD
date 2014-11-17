/***********************************************************************
 * 
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008 Carsten Urbach
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
#ifndef _START_H
#define _START_H


/* functions requesting random numbers can request different distributions by calling _rn_switch
   with the first argument set to a random number type as defined below and a function pointer
   (see start.c for examples) */

#define _rn_switch(type,rn_fn_ptr) \
  switch( type ) { \
    case RN_Z2: \
      rn_fn_ptr = z2_vector; \
      break; \
    case RN_UNIF: \
      rn_fn_ptr = ranlxd; \
      break; \
    case RN_GAUSS: \
    default: \
      rn_fn_ptr = gauss_vector; \
      break; \
  } \

/* RN_GAUSS: gaussian ditributed random numbers
   RN_UNIF:  random numbers drawn from a uniform distribution (this is a simple call to ranlxd!)
   RN_Z2:    z2 noise */

enum RN_TYPE { RN_GAUSS, RN_UNIF, RN_Z2 };

void unit_spinor_field(const int k);
void zero_spinor_field(spinor * const k, const int N);
void constant_spinor_field(spinor * const k, const int p, const int N);

void random_spinor_field_lexic(spinor * const k, const int repro, const enum RN_TYPE rn_type);
void random_spinor_field_eo(spinor * const k, const int repro, const enum RN_TYPE rn_type);

void unit_g_gauge_field(void);

void random_gauge_field(const int repro, su3 ** const gf);

double random_su3adj_field(const int repro, su3adj ** const momenta);

void set_spinor_field(int k, const double c);
void set_gauge_field(const double c);
void set_spinor_point(spinor * s, const double c);
su3 set_su3(const double c);

void source_spinor_field(spinor * const P, spinor * const Q, int is, int ic);
void source_spinor_field_point_from_file(spinor * const P, spinor * const Q, int is, int ic, int source_indx);

void start_ranlux(int level,int seed);

/* This function allows initializing RANLUX from a saved state. IMPORTANT NOTE BELOW:
 * Because of the way that random numbers are used in tmLQCD, this function should only ever be
 * used in "reproducible random numbers" mode with the full understanding that all processes
 * will have exactly the same random number generators!!
 * The main routines in start.c all accomodate this by making every process generate random
 * numbers for the whole volume and only using those relevant for the local volume while throwing
 * all others away.
 * If some function requests random numbers via any of the utility functions or ranlx[d,s] directly
 * without taking this fact into account, all processes will generate the same ones! You have been warned. */
void start_ranlux_from_file(char * const rlxd_state_filename, char * const rlxs_state_filename);
void store_ranlux_state(char * const rlxd_state_filename, char * const rlxs_state_filename);

void gen_test_spinor_field(spinor * const k , const int eoflag);
void write_test_spinor_field(spinor * const k , const int eoflag, char * postfix);
#endif

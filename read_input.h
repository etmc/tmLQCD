/***********************************************************************
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

/* 
 * This is the function to parse the input file.
 * No default values for any paramter will be set
 *
 * read_inputg expects the filename of the input file
 * as an input parameter.
 *
 * read_input returns 2 if the input file did not exist 
 */

#ifndef _PARSER_H
#define _PARSER_H

#define COLD 0
#define HOT 1
#define RESTART 2
#define CONTINUE 3

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */
  
  /* input parameters defined in */
  /* read_input.h */
  extern int verbose;
  extern int startoption;
  extern int Ntherm;
  extern int Nmeas;
  extern int Nsave;
  extern int gmres_m_parameter, gmresdr_nr_ev;
  extern int write_cp_flag;
  extern int cp_interval;
  extern int nstore;
  extern int crylov_space_dim;
  extern char rlxd_input_filename[500];
  extern char gauge_input_filename[500];
  extern int subforwilson_flag;
  extern int eigenvalue_method_flag;
  extern int eigenvalue_max_iterations;
  extern double eigenvalue_precision;
  extern int index_start;
  extern int index_end;
  extern int random_seed;
  extern int rlxd_level;
  extern double X0, X1, X2, X3;
  extern int read_source_flag;
  extern int return_check_flag;
  extern int return_check_interval;
  extern int gauge_precision_read_flag;
  extern int gauge_precision_write_flag;
  extern int reproduce_randomnumber_flag;
  extern double stout_rho;
  extern int stout_no_iter;
  extern int use_stout_flag;
  extern int phmc_no_flavours;
  extern int phmc_heavy_timescale;
  extern int phmc_compute_evs;
  extern int phmc_exact_poly;
  extern int compute_evs;
  extern int no_eigenvalues;
  extern double eigenvalue_precision;
  extern double stilde_max;
  extern double stilde_min;
  extern int degree_of_p;
  extern int source_location;
  extern int sub_evs_cg_flag;
  extern int even_odd_flag;
  extern int bc_flag;
  extern int online_measurement_flag;
  extern int online_measurement_freq;
  extern int reweighting_flag;
  extern int reweighting_samples; 
  extern int no_samples;
  extern int compute_modenumber;
  extern int compute_topsus;
  extern double mstarsq;
  extern int no_sources_z2;
  extern int usegpu_flag;
  extern int max_innersolver_it;
  extern double innersolver_precision;
  extern int device_num;
  
  extern double innersolver_precision_abs;
  extern double innersolver_precision_rel;
  extern int innersolver_precision_check_abs;
  extern int innersolver_precision_check_rel;

  extern int omp_num_threads;

  extern int use_preconditioning;

  extern int nblocks_t, nblocks_x, nblocks_y, nblocks_z;

  int read_input(char *);
  int reread_input(char *);
  
#ifdef __cplusplus
}
#endif

#endif

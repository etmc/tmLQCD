/* $Id$ */

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
  extern int Nskip;
  extern int solver_flag;
  extern int operator_flag;
  extern int matrix_element_flag;
  extern int save_config_flag;
  extern int save_prop_flag;
  extern int save_prop_g2_flag;
  extern int write_cp_flag;
  extern int cp_interval;
  extern int nstore;
  extern int crylov_space_dim;
  extern char rlxd_input_filename[100];
  extern char gauge_input_filename[100];
  extern int subforwilson_flag;
  extern int eigenvalue_method_flag;
  extern int eigenvalue_max_iterations;
  extern double eigenvalue_precision;
  extern int index_start;
  extern int index_end;
  extern int first_prop_flag;
  extern double dtau;
  extern int Nsteps;
  extern double q_off;
  extern double q_off2;
  extern int random_seed;
  extern int integtyp,nsmall;
  extern int ITER_MAX_BCG;
  extern int ITER_MAX_CG;
  extern double X0;
  
  int read_input(char *);
  
#ifdef __cplusplus
}
#endif

#endif

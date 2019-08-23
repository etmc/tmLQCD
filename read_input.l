/**********************************************************************
 *
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008 Carsten Urbach
 *
 * Modified by Jenifer Gonzalez Lopez 2009/03/27
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
 * This is the parser. (Dec 2002)
 * The .c-file is generated from .l using flex.
 * Please edit read_input.l instead of read_input.c!
 * flex should be said to be case insensitive!
 *
 * After modifiing read_input.l please update read_input.c with
 * flex -Ptmlqcd -i -t read_input.l > read_input.c
 *
 * Autor: Carsten Urbach
 *        urbach@physik.fu-berlin.de
 ***********************************************************************/

SPC [[:blank:]]+
CMD [:][[:space:]]+
RLN [1-9(10)(11)(12)(13)(14)(15)(16)][:]
DIGIT [[:digit:]]
ZT [0-9(10)(11)]
IDXEX ("-"{DIGIT}+)
SIGN ("+"|"-")
FLT {SIGN}?{DIGIT}+(".")?{DIGIT}*(e("-"|"+")?{DIGIT}+)?
FILENAME [a-zA-Z0-9_".""-""/"][a-zA-z0-9"."_"-""/"]+
NAME [a-zA-Z0-9_]+
CSTR \"[a-zA-Z0-9\-\._]+\"
TYPE [0-9A-Z]+
STRLIST [0-9"."","" ""_""-"a-zA-Z]*
EQL {SPC}*={SPC}*

%{
#ifdef HAVE_CONFIG_H
#  include<tmlqcd_config.h>
#endif
#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<ctype.h>
#define INIT_GLOBALS
#include"global.h"
#include"read_input.h"
#include"default_input_values.h"
#include"monomial/monomial.h"
#include"solver/solver_types.h"
#include"meas/measurements.h"
#include"integrator.h"
#include"operator.h"
#include"phmc.h"
#include<io/params.h>
#include "qphix_types.h"
#include "quda_types.h"
#include "misc_types.h"

// some general parsing helpers

static inline void rmQuotes(char *str){
  char* strsave=str;

  while(*str== ' ' || *str == '\t') str++;

  if(*str=='\"') *(strsave++)=*( ++str );
  else fprintf(stderr,"Error in removing quotes from string:\n %s  \n" ,  str);
  ++str;

  while( !( *str=='\0' || *str == '\"' || *str == '\n') ) *(strsave++)=*(str++);
  *strsave='\0';
}

/* tokenize the comma-delimited list 'str' and return the first token,
   which is the name of the list 'list = token1, token2, token3, ..',
   via the 'paramname' parameter 
   the length of the array holding 'paramname' parameter must be passed as well
   the memory pointed to by 'input_copy' has to be freed afterwards!*/
static inline void strlist_tokenize(char const * const input,
                                    char * input_copy,
                                    char * const paramname,
                                    const int paramname_length)
{
  char error_message[200];
  char * token = (char*)NULL;
  
  const int input_length = strlen(input) + 1;
  if( input_length > 20000 ){
    yy_fatal_error("Likely overflow bug in 'strlist_tokenize' of 'read_input.l'.\n"
                   "Encountered string of length exceeding 20k characters!\n");
  }
  input_copy = malloc(input_length);
  if( input_copy == NULL ){
    yy_fatal_error("Could not allocate memory for 'input_copy' in 'strlist_tokenize' of 'read_input.l'!\n");
  }
  strncpy(input_copy, input, input_length);

  // the first token is the parameter name, we will return it, 
  // but first check if tokenization was successful
  token = strtok(input_copy, "\n\t =,\\");
  if( (void*)token == NULL ){
    snprintf(error_message, 200, "Unable to parse '%s' in strlist_tokenize!\n", paramname);
    yy_fatal_error(error_message);
  } else {
    // return what we parsed
    snprintf(paramname, paramname_length, "%s", token);
  }
}

/* obtain the next token in the list of comma-separated strings
   'list_end' is a parameter which is set to 1 if the end of the list is reached
   'token' will hold the value of the token
   'token_length' is the length of the array to which the token should be written */
static inline void strlist_next_token(int * const list_end, char * const token_out, const int token_length){
  char * token = (char*)NULL;
  token = strtok(NULL," =,\t");
  if( (void*)token == NULL ){
    *list_end = 1;
  } else {
    *list_end = 0;
    snprintf(token_out, token_length, "%s", token);
  }
}

/* when the list of tokens is a list of doubles or integers, this function
   can be used to obtain the next token in the list directly as a double
   'list_end' is a parameter which is set to 1 if the end of the list 
   has been reached */
static inline double fltlist_next_token(int * const list_end){
  double retval;
  char * token = (char*)NULL;
  token = strtok(NULL," ,=\t");
  if( (void*)token == NULL ){
    *list_end = 1;
    return(0.0);
  } else {
    *list_end = 0;
    sscanf(token,"%lf",&retval);
    return(retval);
  }
}

  /* Name of the parsing routine */
#define YY_DECL         int parse_config()
#define YY_NO_UNPUT
  /* helper vars */
  char *cstring_to_parse=NULL;
  int cstring_caller;
  int solver_caller;

  /* declaration of input parameters */
  int i=0;
  int line_of_file = 1;
  int current_monomial = -1;
  int current_measurement=-1;
  int current_operator = -1;
  extern int no_monomials;
  extern int no_measurements;
  monomial * mnl = NULL;
  measurement * meas = NULL;
  operator * optr = NULL;
  int comment_caller;
  int name_caller;
  int a,b;
  double c;
  float cs;
  int reread = 0;
  char name[100];
  char * type;

  int verbose = 0;
  int myverbose = 0;
  int startoption;
  int Ntherm;
  int Nmeas;
  int Nsave;
  int gmres_m_parameter, gmresdr_nr_ev;
  int write_cp_flag;
  int cp_interval;
  int nstore;
  int index_start, index_end;
  int random_seed;
  int rlxd_level;
  char rlxd_input_filename[500];
  char gauge_input_filename[500];
  int read_source_flag;
  int return_check_flag, return_check_interval;
  int gauge_precision_read_flag;
  int gauge_precision_write_flag;
  int gmres_m_parameter, gmresdr_nr_ev;
  int reproduce_randomnumber_flag;
  double stout_rho;
  int stout_no_iter;
  int use_stout_flag;
  int phmc_no_flavours;
  int phmc_heavy_timescale;
  int phmc_exact_poly;
  int compute_evs;
  int phmc_compute_evs;
  double stilde_max;
  double stilde_min;
  int degree_of_p;
  int source_location;
  int no_eigenvalues;
  double eigenvalue_precision;
  int sub_evs_cg_flag;
  int even_odd_flag;
  int bc_flag;
  int online_measurement_flag;
  int online_measurement_freq;
  int reweighting_flag;
  int reweighting_samples;
  int no_samples;
  int compute_modenumber;
  int compute_topsus;
  double mstarsq;
  int no_sources_z2;
  int usegpu_flag;
  int max_innersolver_it;
  int device_num;
  double innersolver_precision;

  double innersolver_precision_abs;
  double innersolver_precision_rel;
  int innersolver_precision_check_abs;
  int innersolver_precision_check_rel;

  int max_mms_shifts;
  int use_mixed_mms;
  int min_innersolver_it;

  double mixcg_innereps;
  int mixcg_maxinnersolverit;

  int propagator_comparison;
  int nb_cores;

  int omp_num_threads;

  int nblocks_t;
  int nblocks_x;
  int nblocks_y;
  int nblocks_z;

  int dfl_field_iter;
  int dfl_poly_iter;
  double kappa_dflgen;
  double mu_dflgen;
  double kappa_dfl;
  double mu_dfl;
  double kappa_Msap;
  double mu_Msap;

  int use_preconditioning;

  /* macro control over read_input.l is poor. In order to support compilation without
   * a given external library, we need to have these declared here as a dummy.
   * Note that no parameters will actually be passed, the input file reader will
   * terminate the program if one attempts to set the parameters but tmLQCD has
   * been compiled without the interface in question. */
#ifdef TM_USE_QUDA
  extern tm_QudaParams_t quda_input;
#else
  tm_QudaParams_t quda_input;
#endif

#ifdef TM_USE_QPHIX
  extern tm_QPhiXParams_t qphix_input;
#else
  tm_QPhiXParams_t qphix_input;
#endif

  // specific parsing helpers
  static inline void parse_quda_mg_blocksize(char const * const input, const int dim){
    char paramname[100];
    char error_message[200];
    int list_end = 0;
    int level = 0;
    int blocksize = 0;
    
    char * input_copy = (char*)NULL;
    strlist_tokenize(input, input_copy, paramname, 100);
    blocksize = (int)fltlist_next_token(&list_end);
    while( list_end != 1 ){
      if( level >= QUDA_MAX_MG_LEVEL ){
        snprintf(error_message, 200, "Exceeded maximum number of levels (%d-1) parsing %s!\n", QUDA_MAX_MG_LEVEL, paramname);
        yy_fatal_error(error_message);
      }
      
      quda_input.mg_blocksize[level][dim] = blocksize;
      if(myverbose){ 
        printf("  %s, level %d set to %d line %d\n", paramname,
                level, quda_input.mg_blocksize[level][dim], line_of_file);
      }
      
      level++;
      blocksize = (int)fltlist_next_token(&list_end);
    }
    free(input_copy);
  }
  
  static inline void parse_quda_mg_int_par_array(char const * const input, int * par_array){
    char paramname[100];
    char error_message[200];
    int list_end = 0;
    int level = 0;
    char * input_copy = (char*)NULL;
    strlist_tokenize(input, input_copy, paramname, 100);
    int parval = (int)fltlist_next_token(&list_end);
    while( list_end != 1 ){
      if( level >= QUDA_MAX_MG_LEVEL ){
        snprintf(error_message, 200, "Exceeded maximum number of levels (%d) parsing %s!\n", QUDA_MAX_MG_LEVEL, paramname);
        yy_fatal_error(error_message);
      }
      
      par_array[level] = parval;
      if(myverbose){ 
        printf("  %s, level %d set to %d line %d\n", paramname,
                level, par_array[level], line_of_file);
      }
      
      level++;
      parval = (int)fltlist_next_token(&list_end);
    }
    free(input_copy);
  }
  
  static inline void parse_quda_mg_dbl_par_array(char const * const input, double * par_array){
    char paramname[100];
    char error_message[200];
    int list_end = 0;
    int level = 0;
    
    char * input_copy = (char*)NULL;
    strlist_tokenize(input, input_copy, paramname, 100);
    double parval = fltlist_next_token(&list_end);
    while( list_end != 1 ){
      if( level >= QUDA_MAX_MG_LEVEL ){
        snprintf(error_message, 200, "Exceeded maximum number of levels (%d) parsing %s!\n", QUDA_MAX_MG_LEVEL, paramname);
        yy_fatal_error(error_message);
      }
      
      par_array[level] = parval;
      if(myverbose){ 
        printf("  %s, level %d set to %f line %d\n", paramname,
                level, par_array[level], line_of_file);
      }
      
      level++;
      parval = fltlist_next_token(&list_end);
    }
    free(input_copy);
  }

  static inline void parse_quda_mg_bool_par_array(char const * const input, int * par_array){
    char paramname[100];
    char token[100];
    char error_message[200];
    int list_end = 0;
    int level = 0;
    char * input_copy = (char*)NULL;
    strlist_tokenize(input, input_copy, paramname, 100);
    strlist_next_token(&list_end, token, 100);
    while( list_end != 1 ){
      if( level >= QUDA_MAX_MG_LEVEL ){
        snprintf(error_message, 200, "Exceeded maximum number of levels (%d) parsing %s!\n", QUDA_MAX_MG_LEVEL, paramname);
        yy_fatal_error(error_message);
      }
      for(int i = 0; token[i]; i++){
        token[i] = tolower(token[i]);
      }
      if( strcmp(token, "yes") == 0 ){
        par_array[level] = QUDA_BOOLEAN_YES;
      }
      else if( strcmp(token, "no") == 0 ){
        par_array[level] = QUDA_BOOLEAN_NO;
      }
      else {
        snprintf(error_message, 200, "ERROR: Only 'yes' and 'no' are possible for %s!\n", paramname);
        yy_fatal_error(error_message);
      }
      if(myverbose) printf("  %s on level %d set to %s line %d\n", paramname, level, token, line_of_file);
      level++;
      strlist_next_token(&list_end, token, 100);
    }
    free(input_copy);
  }
  
  static inline void parse_quda_mg_verbosity_par_array(char const * const input, QudaVerbosity * par_array){
    char paramname[100];
    char token[100];
    char error_message[200];
    int list_end = 0;
    int level = 0;
    char * input_copy = (char*)NULL;
    strlist_tokenize(input, input_copy, paramname, 100);
    strlist_next_token(&list_end, token, 100);
    while( list_end != 1 ){
      if( level >= QUDA_MAX_MG_LEVEL ){
        snprintf(error_message, 200, "Exceeded maximum number of levels (%d) parsing %s!\n", QUDA_MAX_MG_LEVEL, paramname);
        yy_fatal_error(error_message);
      }
      for(int i = 0; token[i]; i++){
        token[i] = tolower(token[i]);
      }
      if( strcmp(token, "silent") == 0 ){
        par_array[level] = QUDA_SILENT;
      }
      else if( strcmp(token, "summarize") == 0 ){
        par_array[level] = QUDA_SUMMARIZE;
      }
      else if( strcmp(token, "verbose") == 0 ){
        par_array[level] = QUDA_VERBOSE;
      }
      else if( strcmp(token, "debug_verbose") == 0 ){
        par_array[level] = QUDA_DEBUG_VERBOSE;
      }
      else {
        snprintf(error_message, 200, 
                 "ERROR: Unknown verbosity level '%s' for '%s'!"
                   " Only 'silent', 'summarize', 'verbose' and 'debug_verbose' supported!\n",
                 token, paramname);
        yy_fatal_error(error_message);
      }
      if(myverbose) printf("  %s on level %d set to %s line %d\n", paramname, level, token, line_of_file);
      level++;
      strlist_next_token(&list_end, token, 100);
    }
    free(input_copy);
  }
  
%}

%option never-interactive

%x STARTCOND
%x THERMSWEEPS
%x NMEAS
%x KAPPA
%x MUBAR
%x EPSBAR
%x MU
%x CSW
%x SEED
%x RLXDLEVEL
%x NSAVE
%x RLXDINPUTFILE
%x GAUGEINPUTFILE
%x GAUGERPREC
%x GAUGEWPREC
%x DSBLIOCHECK
%x DSBLSRCIOCHECK
%x PRECON
%x WRITECP
%x CPINT
%x NSTORE
%x TT
%x LL
%x LLX
%x LLY
%x LLZ
%x NPROCX
%x NPROCY
%x NPROCZ
%x IOPROC
%x IDX
%x BOUNDT
%x BOUNDX
%x BOUNDY
%x BOUNDZ
%X READSOURCE
%x SOURCEFORMAT
%x SOURCEFILE
%x PROPFILE
%x SOURCETS
%x SOURCETYPE
%x PROPSPLIT
%x NOSAMPLES
%x RELPREC
%x REVCHECK
%x REVINT
%x DEBUG
%x GMRESM
%x GMRESDRNEV
%x REPRORND
%x SLOPPYPREC
%x USESTOUT
%x STOUTRHO
%x STOUTITER
%x COMPUTEEVS
%x SRCLOC
%x SUBEVCG
%x NOEV
%x PRECEV
%x EO
%x BC
%x WRPROPFLAG
%x PROPTYPE
%x COMPUTEMN
%x COMPUTETS
%x MSTARSQ
%x NOSOURCESZ2
%x INITMEASUREMENT
%x ONLINEMEAS
%x PIONNORMMEAS
%x PLOOP
%x ORIENTEDPLAQUETTESMEAS
%x GRADIENTFLOWMEAS

%x REWEIGH
%x REWSAMPLES

%x INITINTEGRATOR
%x INTEGRATOR

%x GPU
%x INITGPU

%x DEFLATION
%x INITDEFLATION

%x MULTIGRID
%x INITMULTIGRID

%x INITOPERATOR
%x TMOP
%x DBTMOP
%x WILSONOP
%x OVERLAPOP
%x CLOVEROP
%x DBCLOVEROP
%x TMSOLVER
%x DBTMSOLVER
%x CSWSOLVER
%x OVSOLVER

%x INITMONOMIAL
%x DETMONOMIAL
%x CLDETMONOMIAL
%x CLDETRATMONOMIAL
%x CLDETRATRWMONOMIAL
%x GAUGEMONOMIAL
%x NDPOLYMONOMIAL
%x NDRATMONOMIAL
%x RATMONOMIAL
%x CLRATMONOMIAL
%x RATCORMONOMIAL
%x CLRATCORMONOMIAL
%x NDCLRATMONOMIAL
%x NDRATCORMONOMIAL
%x NDDETRATMONOMIAL
%x NDCLDETRATMONOMIAL
%x NDCLRATCORMONOMIAL
%x POLYMONOMIAL
%x CLPOLYMONOMIAL
%x MNAME
%x MCSTR
%x MSOLVER
%x RATMSOLVER
%x NDMSOLVER
%x GTYPE

%x COMMENT
%x ERROR

%x PCOMP
%x NBCORES

%x OMPNUMTHREADS

%x DFLNBLOCKT
%x DFLNBLOCKX
%x DFLNBLOCKY
%x DFLNBLOCKZ

%x DFLFIELDITER
%x DFLPOLYITER

%x PRECONDITIONING

%x QUDAINVERTER
%x QPHIXINVERTER
%x COMPRESSION

%x MIXCGEPS
%x MIXCGIT

%x LOWMEM
%x SUBPROCESS 

%x INITEXTERNALINVERTER

%%
^SourceFilename{EQL}               BEGIN(SOURCEFILE);
^PropagatorFilename{EQL}           BEGIN(PROPFILE);
^T{EQL}                            BEGIN(TT);
^L{EQL}                            BEGIN(LL);
^LX{EQL}                           BEGIN(LLX);
^LY{EQL}                           BEGIN(LLY);
^LZ{EQL}                           BEGIN(LLZ);
^NRXProcs{EQL}                     BEGIN(NPROCX);
^NRYProcs{EQL}                     BEGIN(NPROCY);
^NRZProcs{EQL}                     BEGIN(NPROCZ);
^kappa{EQL}                        BEGIN(KAPPA);
^csw{EQL}                          BEGIN(CSW);
^2KappaMu{EQL}                     BEGIN(MU);
^2KappaMubar{EQL}                  BEGIN(MUBAR);
^2KappaEpsBar{EQL}                 BEGIN(EPSBAR);
^NoEigenvalues{EQL}                BEGIN(NOEV);
^EigenvaluePrecision{EQL}          BEGIN(PRECEV);
^seed{EQL}                         BEGIN(SEED);
^StartCondition{EQL}               BEGIN(STARTCOND);
^ThermalisationSweeps{EQL}         BEGIN(THERMSWEEPS);
^Measurements{EQL}                 BEGIN(NMEAS);
^NSave{EQL}                        BEGIN(NSAVE);
^GaugeFieldInFile{EQL}             BEGIN(GAUGEINPUTFILE);
^RlxdStateInFile{EQL}              BEGIN(RLXDINPUTFILE);
^SubtractEVForCG{EQL}              BEGIN(SUBEVCG);
^WriteCheckpoints{EQL}             BEGIN(WRITECP);
^CheckpointInterval{EQL}           BEGIN(CPINT);
^GaugeConfigInputFile{EQL}         BEGIN(GAUGEINPUTFILE);
^RlxdInputFile{EQL}                BEGIN(RLXDINPUTFILE);
^InitialStoreCounter{EQL}          BEGIN(NSTORE);
^StdIOProcessor{EQL}               BEGIN(IOPROC);
^Indices{EQL}                      BEGIN(IDX);
^BCAngleT{EQL}                     BEGIN(BOUNDT);
^ThetaT{EQL}                       BEGIN(BOUNDT);
^ThetaX{EQL}                       BEGIN(BOUNDX);
^ThetaY{EQL}                       BEGIN(BOUNDY);
^ThetaZ{EQL}                       BEGIN(BOUNDZ);
^ReadSource{EQL}                   BEGIN(READSOURCE);
^UseRelativePrecision{EQL}         BEGIN(RELPREC);
^ReversibilityCheck{EQL}           BEGIN(REVCHECK);
^ReversibilityCheckIntervall{EQL}  BEGIN(REVINT);
^DebugLevel{EQL}                   BEGIN(DEBUG);
^GMRESMParameter{EQL}              BEGIN(GMRESM);
^GMRESDRNrEv{EQL}                  BEGIN(GMRESDRNEV);
^GaugeConfigReadPrecision{EQL}     BEGIN(GAUGERPREC);
^GaugeConfigWritePrecision{EQL}    BEGIN(GAUGEWPREC);
^DisableIOChecks{EQL}              BEGIN(DSBLIOCHECK);
^DisableGaugeIOChecks{EQL}         BEGIN(DSBLIOCHECK);
^DisableSourceIOChecks{EQL}        BEGIN(DSBLSRCIOCHECK);
^ReproduceRandomNumbers{EQL}       BEGIN(REPRORND);
^UseSloppyPrecision{EQL}           BEGIN(SLOPPYPREC);
^UseStoutSmearing{EQL}             BEGIN(USESTOUT);
^StoutRho{EQL}                     BEGIN(STOUTRHO);
^StoutNoIterations{EQL}            BEGIN(STOUTITER);
^ComputeEVs{EQL}                   BEGIN(COMPUTEEVS);
^SourceLocation{EQL}               BEGIN(SRCLOC);
^UseEvenOdd{EQL}                   BEGIN(EO);
^Bc{EQL}                           BEGIN(BC);
^WritePropagatorFormat{EQL}        BEGIN(WRPROPFLAG);
^PropagatorType{EQL}               BEGIN(WRPROPFLAG);
^RanluxdLevel{EQL}                 BEGIN(RLXDLEVEL);
^GCRPreconditioner{EQL}            BEGIN(PRECON);
^ComputeReweightingFactor{EQL}     BEGIN(REWEIGH);
^NoReweightingSamples{EQL}         BEGIN(REWSAMPLES);
^SourceTimeSlice{EQL}              BEGIN(SOURCETS);
^SourceType{EQL}                   BEGIN(SOURCETYPE);
^NoSamples{EQL}                    BEGIN(NOSAMPLES);
^SplittedPropagator{EQL}           BEGIN(PROPSPLIT);
^UsePreconditioning{EQL}           BEGIN(PRECONDITIONING);
^UseCompression{EQL}               BEGIN(COMPRESSION);

^BeginMeasurement{SPC}+            BEGIN(INITMEASUREMENT);
^ComputeModeNumber{EQL}            BEGIN(COMPUTEMN);
^ComputeTopSus{EQL}                BEGIN(COMPUTETS);
^MStarSq{EQL}                      BEGIN(MSTARSQ);
^NoSourcesZ2{EQL}                  BEGIN(NOSOURCESZ2);
^BeginMonomial{SPC}+               BEGIN(INITMONOMIAL);
^BeginInt                          BEGIN(INITINTEGRATOR);
^BeginOperator{SPC}+               BEGIN(INITOPERATOR);

^BeginExternalInverter{SPC}+       BEGIN(INITEXTERNALINVERTER);

^PropagatorComparison{EQL}         BEGIN(PCOMP);
^NbCoresPerNode{EQL}               BEGIN(NBCORES);

^OMPNumThreads{EQL}                BEGIN(OMPNUMTHREADS);

^NoBlocksT{EQL}                    BEGIN(DFLNBLOCKT);
^NoBlocksX{EQL}                    BEGIN(DFLNBLOCKX);
^NoBlocksY{EQL}                    BEGIN(DFLNBLOCKY);
^NoBlocksZ{EQL}                    BEGIN(DFLNBLOCKZ);

^DflFieldIter{EQL}                 BEGIN(DFLFIELDITER);
^DflPolyIter{EQL}                  BEGIN(DFLPOLYITER);

^BeginGPU                          BEGIN(INITGPU);
^BeginDeflation		                 BEGIN(INITDEFLATION);
^BeginDDalpha       	             BEGIN(INITMULTIGRID);
^MixCGInnerEps{EQL}                BEGIN(MIXCGEPS);
^MixCGMaxIter{EQL}                 BEGIN(MIXCGIT);

^EnableLowmem{EQL}                 BEGIN(LOWMEM);
^EnableSubprocess{EQL}             BEGIN(SUBPROCESS);


<INITGPU>Init{SPC}* {
 if(myverbose) printf("Initialising GPU line %d\n", line_of_file); 
 usegpu_flag = 1;
  if(myverbose!=0) printf("  Using help of GPU for inversions\n"); 
 BEGIN(GPU);
 }
<GPU>{
  {SPC}*InnersolverPrecision{EQL}{FLT} {
    sscanf(yytext, " %[2a-zA-Z] = %lf", name, &c);
    innersolver_precision = c;
    if(myverbose) printf("  Inner solver precision set to %lf line %d\n", c, line_of_file);
  }
  {SPC}*MaxInnersolverIteration{EQL}{DIGIT}+ {
    sscanf(yytext, " %[2a-zA-Z] = %d", name, &a);
    max_innersolver_it = a;
    if(myverbose) printf("  Inner solver iterations set to %d line %d\n", a, line_of_file);
  }
  {SPC}*DeviceNum{EQL}{DIGIT}+ {
      sscanf(yytext, " %[2a-zA-Z] = %d", name, &a);
      device_num = a;
      if(myverbose) printf("  Device Number set to %d line %d\n", a, line_of_file);
  }
  {SPC}*InnersolverPrecisionAbs{EQL}{FLT} {
    sscanf(yytext, " %[2a-zA-Z] = %lf", name, &c);
    innersolver_precision_abs = c;
    if(myverbose) printf("  innersolver_precision_abs set to %lf line %d\n", c, line_of_file);
  }
  {SPC}*InnersolverPrecisionRel{EQL}{FLT} {
    sscanf(yytext, " %[2a-zA-Z] = %lf", name, &c);
    innersolver_precision_rel = c;
    if(myverbose) printf("  innersolver_precision_rel set to %lf line %d\n", c, line_of_file);
  }
  {SPC}*InnersolverPrecisionCheckAbs{EQL}{DIGIT}+ {
      sscanf(yytext, " %[2a-zA-Z] = %d", name, &a);
      innersolver_precision_check_abs = a;
      if(myverbose) printf("  innersolver_precision_check_abs set to %d line %d\n", a, line_of_file);
  }
  {SPC}*InnersolverPrecisionCheckRel{EQL}{DIGIT}+ {
      sscanf(yytext, " %[2a-zA-Z] = %d", name, &a);
      innersolver_precision_check_rel = a;
      if(myverbose) printf("  innersolver_precision_check_rel set to %d line %d\n", a, line_of_file);
  }
  {SPC}*MinInnerSolverIterations{EQL}{DIGIT}+ {
      sscanf(yytext, " %[2a-zA-Z] = %d", name, &a);
      min_innersolver_it = a;
      if(myverbose) printf("  min_innersolver_it set to %d line %d\n", a, line_of_file);
  }
  {SPC}*MaxMmsShifts{EQL}{DIGIT}+ {
      sscanf(yytext, " %[2a-zA-Z] = %d", name, &a);
      max_mms_shifts = a;
      if(myverbose) printf("  max_mms_shifts set to %d line %d\n", a, line_of_file);
  }
  {SPC}*UseMixedMms{EQL}yes {
    sscanf(yytext, " %[2a-zA-Z] = %lf", name, &c);
    use_mixed_mms = 1;
    if(myverbose) printf("  Using mixed solver for smallest shifts\n");  }
  {SPC}*UseMixedMms{EQL}no {
    sscanf(yytext, " %[2a-zA-Z] = %lf", name, &c);
    use_mixed_mms = 0;
    if(myverbose) printf("  Not using mixed solver for smallest shifts\n");
  }
  EndGPUInit{SPC}* {
  if(myverbose) printf("GPU parsed in line %d\n\n", line_of_file);
  BEGIN(0);
  }
}

<INITDEFLATION>Init{SPC}* {
 if(myverbose) printf("Initialising DEFLATION line %d\n", line_of_file); 
 BEGIN(DEFLATION);
 }
<DEFLATION>{
  {SPC}*DeflationSubspaceDimension{EQL}{DIGIT}+ {
    sscanf(yytext, " %[2a-zA-Z] = %d", name, &a);
    g_N_s = a;
    if(myverbose) printf("Number of global approximate eigenvectors (deflation subspace dimension) set to %d line %d\n", a, line_of_file);
  }
  {SPC}*NiterMsap{EQL}{DIGIT}+ {
    sscanf(yytext, " %[2a-zA-Z] = %d", name, &a);
    NiterMsap = a;
    if(myverbose) printf("NiterMsap for solver preconditioner set to %d line %d\n", a, line_of_file);
  }
  {SPC}*NcycleMsap{EQL}{DIGIT}+ {
    sscanf(yytext, " %[2a-zA-Z] = %d", name, &a);
    NcycleMsap = a;
    if(myverbose) printf("NcycleMsap for solver preconditioner set to %d line %d\n", a, line_of_file);
  }
  {SPC}*KappaMsap{EQL}{FLT}+ {
    sscanf(yytext, " %[2a-zA-Z] = %lf", name, &c);
    kappa_Msap = c;
    if(myverbose) printf("kappa for Msap preconditioner set to %e line %d\n", c, line_of_file);
  }
  {SPC}*2KappaMuMsap{EQL}{FLT}+ {
    sscanf(yytext, " %[2a-zA-Z] = %lf", name, &c);
    mu_Msap = c;
    if(myverbose) printf("2KappaMu for Msap preconditioner set to %e line %d\n", c, line_of_file);
  }
  {SPC}*NiterMsapSubspace{EQL}{DIGIT}+ {
    sscanf(yytext, " %[2a-zA-Z] = %d", name, &a);
    NiterMsap_dflgen = a;
    if(myverbose) printf("NiterMsapDfl for subspace generation set to %d line %d\n", a, line_of_file);
  }
  {SPC}*NcycleMsapSubspace{EQL}{DIGIT}+ {
    sscanf(yytext, " %[2a-zA-Z] = %d", name, &a);
    NcycleMsap_dflgen = a;
    if(myverbose) printf("NcycleMsapDfl for subspace generation set to %d line %d\n", a, line_of_file);
  }
  {SPC}*NsmoothSubspace{EQL}{DIGIT}+ {
    sscanf(yytext, " %[2a-zA-Z] = %d", name, &a);
    NsmoothMsap_dflgen = a;
    if(myverbose) printf("NsmoothMsapDfl for subspace generation set to %d line %d\n", a, line_of_file);
  }
  {SPC}*LittleGMRESMParameter{EQL}{DIGIT}+ {
    sscanf(yytext, " %[2a-zA-Z] = %d", name, &a);
    little_gmres_m_parameter = a;
    if(myverbose) printf("LittleGMRESMParameter set to %d line %d\n", a, line_of_file);
  }
  {SPC}*LittleSolverLowPrecision{EQL}{FLT} {
    sscanf(yytext, " %[2a-zA-Z] = %lf", name, &c);
    little_solver_low_prec = c;
    if(myverbose) printf("Low precision for little Dirac solver set to %e line %d\n", c, line_of_file);
  }
  {SPC}*LittleSolverHighPrecision{EQL}{FLT} {
    sscanf(yytext, " %[2a-zA-Z] = %lf", name, &c);
    little_solver_high_prec = c;
    if(myverbose) printf("High precision for little Dirac solver set to %e line %d\n", c, line_of_file);
  }
  {SPC}*LittleSolverMaxIter{EQL}{DIGIT}+ {
    sscanf(yytext, " %[2a-zA-Z] = %d", name, &a);
    little_solver_max_iter = a;
    if(myverbose) printf("LittleSolverMaxIter set to %d line %d\n", a, line_of_file);
  }
  {SPC}*LittleSolver{EQL}(gcr|mcr|mr|fgmres) {
    little_solver = 0;
    if(myverbose) printf("Littlesolver set to gcr in line %d (currently no other solver available)\n", line_of_file);
  }
  {SPC}*useLittleLittleD{EQL}yes {
    usePL = 1;
    if(myverbose) printf("UsePL set to true line %d\n", line_of_file);
  }
  {SPC}*useLittleLittleD{EQL}no {
    usePL = 0;
    if(myverbose) printf("UsePL set to true line %d\n", line_of_file);
  }
  {SPC}*LittleEvenOdd{EQL}yes {
    little_evenodd = 1;
    if(myverbose) printf("LittleEvenOdd set to true line %d\n", line_of_file);
  }
  {SPC}*LittleEvenOdd{EQL}no {
    little_evenodd = 0;
    if(myverbose) printf("LittleEvenOdd set to true line %d\n", line_of_file);
  }
  {SPC}*KappaSubspace{EQL}{FLT} {
    sscanf(yytext, " %[2a-zA-Z] = %lf", name, &c);
    kappa_dflgen = c;
    if(myverbose) printf("kappa for subspace generation set to %e line %d\n", c, line_of_file);
  }
  {SPC}*2KappaMuSubspace{EQL}{FLT} {
    sscanf(yytext, " %[2a-zA-Z] = %lf", name, &c);
    mu_dflgen = c;
    if(myverbose) printf("2KappaMu for subspace generation set to %e line %d\n", c, line_of_file);
  }
  {SPC}*kappa{EQL}{FLT} {
    sscanf(yytext, " %[2a-zA-Z] = %lf", name, &c);
    kappa_dfl = c;
    if(myverbose) printf("kappa for preconditioner set to %e line %d\n", c, line_of_file);
  }
  {SPC}*2KappaMu{EQL}{FLT} {
    sscanf(yytext, " %[2a-zA-Z] = %lf", name, &c);
    mu_dfl = c;
    if(myverbose) printf("2KappaMu for preconditioner set to %e line %d\n", c, line_of_file);
  }
  EndDeflationInit{SPC}* {
  if(myverbose) printf("DEFLATION parsed in line %d\n\n", line_of_file);
  BEGIN(0);
  }
}

<INITMULTIGRID>AMG{SPC}* {
#ifdef DDalphaAMG
 if(myverbose) printf("Initialising DDalphaAMG line %d\n", line_of_file); 
 BEGIN(MULTIGRID);
#else
 printf("ERROR line %d: DDalphaAMG library not included\n", line_of_file);
 exit(1);
#endif
 }
<MULTIGRID>{
  {SPC}*MGSetupIter{EQL}{DIGIT}+ {
    sscanf(yytext, " %[a-zA-Z] = %d", name, &a);
    mg_setup_iter=a;
    if(myverbose) printf("  MG_Setup_Iter set to %d line %d operator %d\n", mg_setup_iter, line_of_file, current_operator);
  }
  {SPC}*MGCoarseSetupIter{EQL}{DIGIT}+ {
    sscanf(yytext, " %[a-zA-Z] = %d", name, &a);
    mg_coarse_setup_iter=a;
    if(myverbose) printf("  MG_Coarse_Setup_Iter set to %d line %d operator %d\n", mg_coarse_setup_iter, line_of_file, current_operator);
  }
  {SPC}*MGNumberOfVectors{EQL}{DIGIT}+ {
    sscanf(yytext, " %[a-zA-Z] = %d", name, &a);
    mg_Nvec=a;
    if(myverbose) printf("  MGNumberOfVectors set to %d line %d operator %d\n", mg_Nvec, line_of_file, current_operator);
  }
  {SPC}*MGNumberOfLevels{EQL}{DIGIT}+ {
    sscanf(yytext, " %[a-zA-Z] = %d", name, &a);
    mg_lvl=a;
    if(myverbose) printf("  MGNumberOfLevels set to %d line %d operator %d\n", mg_lvl, line_of_file, current_operator);
  }
  {SPC}*MGBlockX{EQL}{DIGIT}+ {
    sscanf(yytext, " %[a-zA-Z] = %d", name, &a);
    mg_blk[3]=a;
    if(myverbose) printf("  MGBlockX set to %d line %d operator %d\n", mg_blk[3], line_of_file, current_operator);
  }
  {SPC}*MGBlockY{EQL}{DIGIT}+ {
    sscanf(yytext, " %[a-zA-Z] = %d", name, &a);
    mg_blk[2]=a;
    if(myverbose) printf("  MGBlockY set to %d line %d operator %d\n", mg_blk[2], line_of_file, current_operator);
  }
  {SPC}*MGBlockZ{EQL}{DIGIT}+ {
    sscanf(yytext, " %[a-zA-Z] = %d", name, &a);
    mg_blk[1]=a;
    if(myverbose) printf("  MGBlockZ set to %d line %d operator %d\n", mg_blk[1], line_of_file, current_operator);
  }
  {SPC}*MGBlockT{EQL}{DIGIT}+ {
    sscanf(yytext, " %[a-zA-Z] = %d", name, &a);
    mg_blk[0]=a;
    if(myverbose) printf("  MGBlockT set to %d line %d operator %d\n", mg_blk[0], line_of_file, current_operator);
  }
  {SPC}*MGSetup2KappaMu{EQL}{FLT}+ {
    sscanf(yytext, " %[2a-zA-Z] = %lf", name, &c);
    mg_setup_mu=c;
    mg_setup_mu_set=1;
    if(myverbose) printf("  MGSetup2KappaMu set to %f line %d operator %d\n", mg_setup_mu, line_of_file, current_operator);
  }
  {SPC}*MGCoarseMuFactor{EQL}{FLT}+ {
    sscanf(yytext, " %[a-zA-Z] = %lf", name, &c);
    mg_cmu_factor=c;
    if(myverbose) printf("  MGCoarseMuFactor set to %f line %d operator %d\n", mg_cmu_factor, line_of_file, current_operator);
  }
  {SPC}*MGMixedPrecision{EQL}yes {
    sscanf(yytext, " %[a-zA-Z] = %lf", name, &c);
    mg_mixed_prec = 1;
    if(myverbose) printf("  MGMixedPrecision set to YES line %d operator %d\n",  line_of_file, current_operator);
  }
  {SPC}*MGMixedPrecision{EQL}no {
    sscanf(yytext, " %[a-zA-Z] = %lf", name, &c);
    mg_mixed_prec = 0;
    if(myverbose) printf("  MGMixedPrecision set to NO line %d operator %d\n",  line_of_file, current_operator);
  }
  {SPC}*MGdtauUpdate{EQL}{FLT}+ {
    sscanf(yytext, " %[a-zA-Z] = %lf", name, &c);
    mg_dtau_update=c;
    if(myverbose) printf("  MGdtauUpdate set to %f line %d operator %d\n", mg_dtau_update, line_of_file, current_operator);
  }
  {SPC}*MGrhoUpdate{EQL}{FLT}+ {
    sscanf(yytext, " %[a-zA-Z] = %lf", name, &c);
    mg_rho_update=c;
    if(myverbose) printf("  MG_rho_Update set to %f line %d operator %d\n", mg_rho_update, line_of_file, current_operator);
  }
  {SPC}*MGUpdateSetupIter{EQL}{DIGIT}+ {
    sscanf(yytext, " %[a-zA-Z] = %d", name, &a);
    mg_update_setup_iter=a;
    if(myverbose) printf("  MG_Update_Setup_Iter set to %d line %d operator %d\n", mg_update_setup_iter, line_of_file, current_operator);
  }
  {SPC}*MGOMPNumThreads{EQL}{DIGIT}+ {
    sscanf(yytext, " %[a-zA-Z] = %d", name, &a);
    mg_omp_num_threads=a;
    if(myverbose) printf("  MG_omp_num_threads set to %d line %d operator %d\n", mg_omp_num_threads, line_of_file, current_operator);
  }
  {SPC}*MGNumberOfShifts{EQL}{DIGIT}+ {
    sscanf(yytext, " %[a-zA-Z] = %d", name, &a);
    mg_no_shifts=a;
    // when the number of shifts is specified, mg_mss_mass must be set to zero!
    mg_mms_mass=0;
    if(myverbose) printf("  MG_Num_of_shifts set to %d line %d operator %d\n", mg_no_shifts, line_of_file, current_operator);
  }
  {SPC}*MGMMSMass{EQL}{FLT}+ {
    sscanf(yytext, " %[a-zA-Z] = %lf", name, &c);
    mg_mms_mass=c;
    // when mg_mms_mass is specified, mg_no_shifts should be set to zero!
    mg_no_shifts=0;
    if(myverbose) printf("  MG_MMS_Mass set to %f line %d operator %d\n", mg_mms_mass, line_of_file, current_operator);
  }
  EndDDalphaAMG{SPC}* {
  if(myverbose) printf("DDalphaAMG parsed in line %d\n\n", line_of_file);
  BEGIN(0);
  }
}


<INITOPERATOR>{TYPE} {
  current_operator++;
  optr = &operator_list[current_operator];
  optr->id = current_operator;
  optr->initialised = 0;
  if(strcmp(yytext, "WILSON")==0) {
    optr->type = WILSON;
  }
  else if(strcmp(yytext, "TMWILSON")==0) {
    optr->type = TMWILSON;
  }
  else if(strcmp(yytext, "CLOVER")==0) {
    optr->type = CLOVER;
  }
  else if(strcmp(yytext, "DBCLOVER")==0) {
    optr->type = DBCLOVER;
  }
  else if(strcmp(yytext, "DBTMWILSON")==0) {
    optr->type = DBTMWILSON;
  }
  else if(strcmp(yytext, "OVERLAP")==0) {
    optr->type = OVERLAP;
  }
  else {
    fprintf(stderr, "Unknown operator type %s in line %d\n", yytext, line_of_file);
    exit(1);
  }
  if(!reread) {
    if(add_operator(optr->type) < 0) {
      fprintf(stderr, "Something went wrong in adding operators\nAborting...!\n");
      exit(1);
    }
  }
  if(myverbose) printf("initialising operator with type %s (%d) line %d\n", yytext, optr->type, line_of_file);
  if(myverbose) printf("operator has id %d\n", current_operator);

  if(optr->type == WILSON) BEGIN(WILSONOP);
  else if(optr->type == CLOVER) BEGIN(CLOVEROP);
  else if(optr->type == TMWILSON) BEGIN(TMOP);
  else if(optr->type == DBTMWILSON) BEGIN(DBTMOP);
  else if(optr->type == DBCLOVER) BEGIN(DBCLOVEROP);
  else BEGIN(OVERLAPOP);
}

<WILSONOP,TMOP,OVERLAPOP,DBTMOP,CLOVEROP,DBCLOVEROP>{
  {SPC}*kappa{EQL}{FLT} {
    sscanf(yytext, " %[2a-zA-Z] = %lf", name, &c);
    optr->kappa = c;
    if(myverbose) printf("  kappa set to %f line %d operator %d\n", c, line_of_file, current_operator);
  }
  {SPC}*MaxSolverIterations{EQL}{DIGIT}+ {
    sscanf(yytext, " %[a-zA-Z] = %d", name, &a);
    optr->maxiter = a;
    if(myverbose) printf("  MaxSolverIterations set to %d line %d operator %d\n", a, line_of_file, current_operator);
  }
  {SPC}*PropagatorPrecision{EQL}32 {
    optr->prop_precision = 32;
    PropInfo.precision = 32;
    if(myverbose) printf("  PropagatorPrecision set to 32 line %d operator %d\n", line_of_file, current_operator);
  }
  {SPC}*PropagatorPrecision{EQL}64 {
    optr->prop_precision = 64;
    PropInfo.precision = 64;
    if(myverbose) printf("  PropagatorPrecision set to 64 line %d operator %d\n", line_of_file, current_operator);
  }
  {SPC}*WritePropagator{EQL}yes {
    sscanf(yytext, " %[2a-zA-Z] = %lf", name, &c);
    optr->write_prop_flag = 1;
    if(myverbose) printf("  WritePropagator set to YES line %d operator %d\n",  line_of_file, current_operator);
  }
  {SPC}*WritePropagator{EQL}no {
    sscanf(yytext, " %[2a-zA-Z] = %lf", name, &c);
    optr->write_prop_flag = 0;
    if(myverbose) printf("  WritePropagator set to NO line %d operator %d\n",  line_of_file, current_operator);
  }
  {SPC}*SolverPrecision{EQL}{FLT} {
    sscanf(yytext, " %[2a-zA-Z] = %lf", name, &c);
    optr->eps_sq = c;
    (optr->solver_params).eigcg_tolsq = c;
    if(myverbose) printf("  SolverPrecision set to %lf line %d operator %d\n", c, line_of_file, current_operator);
    if(myverbose) printf("  EigCGtolsq set to %lf line %d operator %d\n", c, line_of_file, current_operator);
  }
  {SPC}*SolverRelativePrecision{EQL}yes {
    sscanf(yytext, " %[2a-zA-Z] = %lf", name, &c);
    optr->rel_prec = 1;
    if(myverbose) printf("  SolverRelativePrecision set to YES line %d operator %d\n",  line_of_file, current_operator);
  }
  {SPC}*SolverRelativePrecision{EQL}no {
    sscanf(yytext, " %[2a-zA-Z] = %lf", name, &c);
    optr->rel_prec = 0;
    if(myverbose) printf("  SolverRelativePrecision set to NO line %d operator %d\n",  line_of_file, current_operator);
  }
  {SPC}*mcgdelta{EQL}{FLT} {
    sscanf(yytext, " %[a-zA-Z1] = %lf", name, &c);
    (optr->solver_params).mcg_delta = c;
    if(myverbose) printf("  mcg_delta set to %lf line %d operator %d\n", c, line_of_file, current_operator);
  }
  {SPC}*EigCGnrhs{EQL}{DIGIT}+ {
    sscanf(yytext, " %[a-zA-Z] = %d", name, &a);
    (optr->solver_params).eigcg_nrhs = a;
    if(myverbose) printf("  EigCGnrhs set to %d line %d operator %d\n", a, line_of_file, current_operator);
  }
  {SPC}*EigCGnrhs1{EQL}{DIGIT}+ {
    sscanf(yytext, " %[a-zA-Z1] = %d", name, &a);
    (optr->solver_params).eigcg_nrhs1 = a;
    if(myverbose) printf("  EigCGnrhs1 set to %d line %d operator %d\n", a, line_of_file, current_operator);
  }
  {SPC}*EigCGnev{EQL}{DIGIT}+ {
    sscanf(yytext, " %[a-zA-Z] = %d", name, &a);
    (optr->solver_params).eigcg_nev = a;
    if(myverbose) printf("  EigCGnev set to %d line %d operator %d\n", a, line_of_file, current_operator);
  }
  {SPC}*EigCGvmax{EQL}{DIGIT}+ {
    sscanf(yytext, " %[a-zA-Z] = %d", name, &a);
    (optr->solver_params).eigcg_vmax = a;
    if(myverbose) printf("  EigCGvmax set to %d line %d operator %d\n", a, line_of_file, current_operator);
  }
  {SPC}*EigCGldh{EQL}{DIGIT}+ {
    sscanf(yytext, " %[a-zA-Z] = %d", name, &a);
    (optr->solver_params).eigcg_ldh = a;
    if(myverbose) printf("  EigCGldh set to %d line %d operator %d\n", a, line_of_file, current_operator);
  }
  {SPC}*EigCGtolsq1{EQL}{FLT} {
    sscanf(yytext, " %[a-zA-Z1] = %lf", name, &c);
    (optr->solver_params).eigcg_tolsq1 = c;
    if(myverbose) printf("  EigCGtolsq1 set to %lf line %d operator %d\n", c, line_of_file, current_operator);
  }
  {SPC}*EigCGrestolsq{EQL}{FLT} {
    sscanf(yytext, " %[a-zA-Z] = %lf", name, &c);
    (optr->solver_params).eigcg_restolsq = c;
    if(myverbose) printf("  EigCGrestolsq set to %lf line %d operator %d\n", c, line_of_file, current_operator);
  }
  {SPC}*EigCGRandGuessOpt{EQL}{DIGIT}+ {
    sscanf(yytext, " %[a-zA-Z] = %d", name, &a);
    (optr->solver_params).eigcg_rand_guess_opt = a;
    if(myverbose) printf("  EigCGrand_guess_opt set to %d line %d operator %d\n", a, line_of_file, current_operator);
  }  
  {SPC}*ExtraMasses{EQL}{STRLIST} {
    char * token = NULL;
    optr->no_extra_masses = 0;
    double mass;
    
    const int input_length = strlen(yytext) + 1;
    if( input_length > 20000 ){
      yy_fatal_error("Likely overflow bug in parsing of 'ExtraMasses{STRLIST}' of 'read_input.l'.\n"
                     "Encountered string of length exceeding 20k characters!\n");
    }
    char * input_copy = malloc(input_length);
    if( input_copy == NULL ){
      yy_fatal_error("Could not allocate memory for 'input_copy' in parsing of 'ExtraMasses{STRLIST}' of 'read_input.l'!\n");
    }
    strncpy(input_copy, yytext, input_length);
    if( strtok(input_copy,"\n\t =,\\") != NULL ) {
       /* Implicitly drop the first token, it is "ExtraMasses" */
      token = strtok(NULL," =,\t");

      while( token != NULL ) {
        if( optr->no_extra_masses >= MAX_EXTRA_MASSES ) {
          yy_fatal_error("  CGMMS maximum number of extra masses reached. Increase MAX_EXTRA_MASSES in global.h!");
        }
        sscanf(token,"%lf",&mass);
        optr->extra_masses[optr->no_extra_masses] = mass;
        if(myverbose) printf("  CGMMS extra mass %d = %lf (2*kappa*mu) added line %d operator %d\n",optr->no_extra_masses,mass,line_of_file,current_operator);
        ++(optr->no_extra_masses);
        token = strtok(NULL," =,\t");
      }
    } else {
      yy_fatal_error("  CGMMS Failed to read ExtraMasses input line in configuration file!");
    }

    if( optr->no_extra_masses == 0 ) {
      yy_fatal_error("  CGMMS Failed to add any extra masses. Input line must be malformed.");
    } else {    
      if(myverbose) printf("  CGMMS %d extra masses added line %d operator %d\n",optr->no_extra_masses,line_of_file,current_operator);
    }
    free(input_copy);
  }
  {SPC}*ExtraMasses{EQL}{FILENAME} {
    FILE * ifs;
    double tempmass = 0.0;
    optr->no_extra_masses = 0;
    char * token = NULL;
    const int input_length = strlen(yytext) + 1;
    if( input_length > 20000 ){
      yy_fatal_error("Likely overflow bug in parsing of 'ExtraMasses{FILENAME}' of 'read_input.l'.\n"
                     "Encountered string of length exceeding 20k characters!\n");
    }
    char * input_copy = malloc(input_length);
    if( input_copy == NULL ){
      yy_fatal_error("Could not allocate memory for 'input_copy' in parsing of 'ExtraMasses{FILENAME}' of 'read_input.l'!\n");
    }
    strncpy(input_copy, yytext, input_length);
    if( strtok(input_copy,"\n\t =,\\") != NULL ) {
      /* drop the first token, it is ExtraMasses */
      token = strtok(NULL," =\t");
      if( token != NULL ) {
        printf("  CGMMS Reading extra masses input file %s\n",token);
        if ((ifs = fopen(token, "r")) != NULL) {
          while ( fscanf( ifs, "%lf", &tempmass ) != EOF ) {
            if( optr->no_extra_masses >= MAX_EXTRA_MASSES ) {
              yy_fatal_error("  CGMMS maximum number of extra masses reached. Increase MAX_EXTRA_MASSES in global.h!");
            } 
            optr->extra_masses[optr->no_extra_masses] = tempmass;
            if (myverbose) {
              printf("  CGMMS Extra mass %d = %lf (2*kappa*mu) added line %d operator %d\n",
                optr->no_extra_masses, optr->extra_masses[optr->no_extra_masses],line_of_file,current_operator);
            }
            ++(optr->no_extra_masses);
          }
          if( optr->no_extra_masses == 0 ) {
            yy_fatal_error("  CGMMS Failed to add any extra masses. Extra masses input file must be malformed!");
          }
          else {
            if(myverbose) 
              printf("  CGMMS %d extra masses added line %d operator %d\n",optr->no_extra_masses,line_of_file,current_operator);
          }
          fclose(ifs);
        }
        else {
          fprintf(stderr, "Could not open extra masses input file %s\n",token);
          optr->no_extra_masses=0;
        }
      }
    }
    free(input_copy);
  }
  ^EndOperator{SPC}* {
    if(myverbose) printf("operator %d parsed line %d\n\n", current_operator, line_of_file);
    BEGIN(0);
  }
}

<WILSONOP,TMOP>{
  {SPC}*Solver{EQL} {
    name_caller = YY_START; 
    BEGIN(TMSOLVER);
  }
}

<CLOVEROP>{
  {SPC}*Solver{EQL} {
    name_caller = YY_START; 
    BEGIN(CSWSOLVER);
  }
}

<WILSONOP,TMOP,CLOVEROP>{
  {SPC}*UseEvenOdd{EQL}yes {
    if(myverbose) printf("  Use even/odd preconditioning line %d operator %d\n", line_of_file, current_operator);
    optr->even_odd_flag = 1;
  }
  {SPC}*UseEvenOdd{EQL}no {
    if(myverbose) printf("  Do not use even/odd preconditioning line %d operator %d\n", line_of_file, current_operator);
    optr->even_odd_flag = 0;
  }
}

<DBTMOP,DBCLOVEROP>{
  {SPC}*2KappaMubar{EQL}{FLT} {
    sscanf(yytext, " %[2a-zA-Z] = %lf", name, &c);
    optr->mubar = c;
    if(myverbose) printf("  2KappaMubar set to %f line %d operator %d\n", c, line_of_file, current_operator);
  }
  {SPC}*2KappaEpsbar{EQL}{FLT} {
    sscanf(yytext, " %[2a-zA-Z] = %lf", name, &c);
    optr->epsbar = c;
    if(myverbose) printf("  2KappaEpsbar set to %f line %d operator %d\n", c, line_of_file, current_operator);
  }
  {SPC}*Solver{EQL} {
    name_caller = YY_START; 
    BEGIN(DBTMSOLVER);
  }
}

<TMOP,CLOVEROP>{
  {SPC}*2KappaMu{EQL}{FLT} {
    sscanf(yytext, " %[2a-zA-Z] = %lf", name, &c);
    optr->mu = c;
    if(myverbose) printf("  2KappaMu set to %f line %d operator %d\n", c, line_of_file, current_operator);
  }
  {SPC}*AddDownPropagator{EQL}yes {
    optr->DownProp = 1;
    if(myverbose) printf("  Invert for + and - mu set in line %d operator %d\n", line_of_file, current_operator);
  }
  {SPC}*AddDownPropagator{EQL}no {
    optr->DownProp = 0;
    if(myverbose) printf("  Don't invert for + and - mu set in line %d operator %d\n", line_of_file, current_operator);
  }
}

<CLOVEROP,DBCLOVEROP>{
  {SPC}*csw{EQL}{FLT} {
    sscanf(yytext, " %[2a-zA-Z] = %lf", name, &c);
    optr->c_sw = c;
    if(myverbose) printf("  Set c_sw set to %lf in line %d for operator %d\n", optr->c_sw, line_of_file, current_operator);
  }
}

<WILSONOP,TMOP,DBTMOP,CLOVEROP,DBCLOVEROP>{  
  {SPC}*UseExternalInverter{EQL}quda {
    if(myverbose) printf("  Use Quda inverter line %d operator %d\n", line_of_file, current_operator);
    optr->external_inverter = QUDA_INVERTER;
  }
  {SPC}*UseExternalInverter{EQL}qphix {
    if(myverbose) printf("  Use QPhiX inverter line %d operator %d\n", line_of_file, current_operator);
    optr->external_inverter = QPHIX_INVERTER;
  }
  {SPC}*UseExternalInverter{EQL}no {
    if(myverbose) printf("  Do not use external inverter line %d operator %d\n", line_of_file, current_operator);
    optr->external_inverter = NO_EXT_INV;
  }
  {SPC}*UseSloppyPrecision{EQL}yes {
    if(myverbose) printf("  Use use sloppy precision (single) in the inverter (if supported by the inverter) line %d operator %d\n", line_of_file, current_operator);
    optr->sloppy_precision = SLOPPY_SINGLE;
  }
  {SPC}*UseSloppyPrecision{EQL}float {
    if(myverbose) printf("  Use use sloppy precision (single) in the inverter (if supported by the inverter) line %d operator %d\n", line_of_file, current_operator);
    optr->sloppy_precision = SLOPPY_SINGLE;
  }
  {SPC}*UseSloppyPrecision{EQL}single {
    if(myverbose) printf("  Use use sloppy precision (single) in the inverter (if supported by the inverter) line %d operator %d\n", line_of_file, current_operator);
    optr->sloppy_precision = SLOPPY_SINGLE;
  }
  {SPC}*UseSloppyPrecision{EQL}no {
    if(myverbose) printf("  Use use sloppy precision (single) in the inverter line %d operator %d\n", line_of_file, current_operator);
    optr->sloppy_precision = SLOPPY_DOUBLE;
  }
  {SPC}*UseSloppyPrecision{EQL}double {
    if(myverbose) printf("  Use use sloppy precision (single) in the inverter line %d operator %d\n", line_of_file, current_operator);
    optr->sloppy_precision = SLOPPY_DOUBLE;
  }
  {SPC}*UseSloppyPrecision{EQL}half {
    if(myverbose) printf("  Use use sloppy precision (half) in the inverter (if supported by the inverter) line %d operator %d\n", line_of_file, current_operator);
    optr->sloppy_precision = SLOPPY_HALF;
  }
  {SPC}*UseCompression{EQL}12 {
    if(myverbose) printf("  Use 12 compression in the inverter (if supported) line %d operator %d\n", line_of_file, current_operator);
    optr->compression_type = COMPRESSION_12;
  }
  {SPC}*UseCompression{EQL}8 {
    if(myverbose) printf("  Use 8 compression in the inverter (if supported) line %d operator %d\n", line_of_file, current_operator);
    optr->compression_type = COMPRESSION_8;
  }
  {SPC}*UseCompression{EQL}18 {
    if(myverbose) printf("  Not using compression in the inverter line %d operator %d\n", line_of_file, current_operator);
    optr->compression_type = NO_COMPRESSION;
  }
}

<OVERLAPOP>{
  {SPC}*Solver{EQL} { 
   BEGIN(OVSOLVER); 
  }
  {SPC}*m{EQL}{FLT} {
    sscanf(yytext, " %[2a-zA-Z] = %lf", name, &c);
    optr->m = c;
    if(myverbose) printf("  m set to %f line %d operator %d\n", c, line_of_file, current_operator);
  }
  {SPC}*s{EQL}{FLT} {
    sscanf(yytext, " %[2a-zA-Z] = %lf", name, &c);
    optr->s = c;
    if(myverbose) printf("  s set to %f line %d operator %d\n", c, line_of_file, current_operator);
  }
  {SPC}*DegreeOfPolynomial{EQL}{DIGIT}+ {
    sscanf(yytext, " %[a-zA-Z] = %d", name, &a);
    optr->deg_poly = a;
    if(myverbose) printf("  DegreeOfPolynomial set to %d line %d operator %d\n", a, line_of_file, current_operator);
  }
  {SPC}*NoKernelEigenvalues{EQL}{DIGIT}+ {
    sscanf(yytext, " %[a-zA-Z] = %d", name, &a);
    optr->no_ev = a;
    if(myverbose) printf("  NoKernelEigenvalues set to %d line %d operator %d\n", a, line_of_file, current_operator);
  }
  {SPC}*KernelEigenvaluePrecision{EQL}{FLT} {
    sscanf(yytext, " %[2a-zA-Z] = %lf", name, &c);
    optr->ev_prec = c;
    if(myverbose) printf("  KernelEigenvaluePrecision set to %f line %d operator %d\n", c, line_of_file, current_operator);
  }
  {SPC}*KernelEigenvectorsReadWrite{EQL}yes {
    optr->ev_readwrite = 1;
    if(myverbose) printf("  KernelEigenvectorsReadWrite set to 1 line %d operator %d\n", line_of_file, current_operator);
  }
  {SPC}*KernelEigenvectorsReadWrite{EQL}no {
    optr->ev_readwrite = 0;
    if(myverbose) printf("  KernelEigenvectorsReadWrite set to 0 line %d operator %d\n", line_of_file, current_operator);
  }
}

<TMSOLVER,CSWSOLVER>{
  mg {
    optr->solver=MG;
    if(myverbose) printf("  Solver set to MG line %d operator %d\n", line_of_file, current_operator);
    BEGIN(name_caller);
  }
}

<DBTMSOLVER,TMSOLVER,CSWSOLVER>{
  cg {
    optr->solver=CG;
    if(myverbose) printf("  Solver set to CG line %d operator %d\n", line_of_file, current_operator);
    BEGIN(name_caller);
  }
  bicgstab {
    optr->solver=BICGSTAB;
    if(myverbose) printf("  Solver set to BiCGstab line %d operator %d\n", line_of_file, current_operator);
    BEGIN(name_caller);
  }
  mixedbicgstab {
    optr->solver=MIXEDBICGSTAB;
    if(myverbose) printf("  Solver set to MixedBiCGstab line %d operator %d\n", line_of_file, current_operator);
    BEGIN(name_caller);
  }
  mixedcg {
    optr->solver=MIXEDCG;
    if(myverbose) printf("  Solver set to MixedCG line %d operator %d\n", line_of_file, current_operator);
    BEGIN(name_caller);
  }
  rgmixedcg {
    optr->solver=RGMIXEDCG;
    if(myverbose) printf("  Solver set to RGMixedCG line %d operator %d\n", line_of_file, current_operator);
    BEGIN(name_caller);
  }
  DDalphaAMG {
#ifdef DDalphaAMG
    optr->solver = MG;
    if(myverbose) printf("  Solver set to DDalphaAMG line %d operator %d\n", line_of_file, current_operator);
    BEGIN(name_caller);
#else
    printf("ERROR line %d operator %d: DDalphaAMG library not included\n", line_of_file, current_operator);
    exit(1);
#endif
  }
  dummyhermtest {
    optr->solver=DUMMYHERMTEST;
    if(myverbose) printf("  Solver set to DummyHermTest line %d operator %d\n", line_of_file, current_operator);
    BEGIN(name_caller);
  }
}

<TMSOLVER>{
  mixedcg {
    optr->solver=MIXEDCG;
    if(myverbose) printf("  Solver set to MixedCG line %d operator %d\n", line_of_file, current_operator);
    BEGIN(name_caller);
  }
  bicgstab {
    optr->solver=BICGSTAB;
    if(myverbose) printf("  Solver set to BiCGstab line %d operator %d\n", line_of_file, current_operator);
    BEGIN(name_caller);
  }
  mixedbicgstab {
    optr->solver=MIXEDBICGSTAB;
    if(myverbose) printf("  Solver set to MixedBiCGstab line %d operator %d\n", line_of_file, current_operator);
    BEGIN(name_caller);
  }
  bicg {
    optr->solver=BICG;
    if(myverbose) printf("  Solver set to BiCG line %d operator %d\n", line_of_file, current_operator);
    BEGIN(name_caller);
  }
  pcg {
    optr->solver=PCG;
    if(myverbose) printf("  Solver set to PCG line %d operator %d\n", line_of_file, current_operator);
    BEGIN(name_caller);
  }
  gmres {
    optr->solver=GMRES;
    if(myverbose) printf("  Solver set to GMRES line %d operator %d\n", line_of_file, current_operator);
    BEGIN(name_caller);
  }
  gcr {
    optr->solver=GCR;
    if(myverbose) printf("  Solver set to GCR line %d operator %d\n", line_of_file, current_operator);
    BEGIN(name_caller);
  }
  mcr {
    optr->solver=MCR;
    if(myverbose) printf("  Solver set to mCR line %d operator %d\n", line_of_file, current_operator);
    BEGIN(name_caller);
  }
  cr {
    optr->solver=CR;
    if(myverbose) printf("  Solver set to CR line %d operator %d\n", line_of_file, current_operator);
    BEGIN(name_caller);
  }
  gmresdr {
    optr->solver=GMRESDR;
    if(myverbose) printf("  Solver set to GMRES-DR line %d operator %d\n", line_of_file, current_operator);
    BEGIN(name_caller);
  }
  cgs {
    optr->solver=CGS;
    if(myverbose) printf("  Solver set to CGS line %d operator %d\n", line_of_file, current_operator);
    BEGIN(name_caller);
  }
  mr {
    optr->solver=MR;
    if(myverbose) printf("  Solver set to MR line %d operator %d\n", line_of_file, current_operator);
    BEGIN(name_caller);
  }
  fgmres {
    optr->solver=FGMRES;
    if(myverbose) printf("  Solver set to FGMRES line %d operator %d\n", line_of_file, current_operator);
    BEGIN(name_caller);
  }
  dflgcr {
    optr->solver=DFLGCR;
    g_dflgcr_flag = 1;
    if(myverbose) printf("  Solver set to DFL-GCR line %d operator %d\n", line_of_file, current_operator);
    BEGIN(name_caller);
  }
  dflfgmres {
    optr->solver=DFLFGMRES;
    g_dflgcr_flag = 1;
    if(myverbose) printf("  Solver set to DFL-FGMRES line %d operator %d\n", line_of_file, current_operator);
    BEGIN(name_caller);
  }
  cgmms {
    optr->solver = CGMMS;
    if(myverbose) printf("  Solver set to CGMMS line %d operator %d\n", line_of_file, current_operator);
    BEGIN(name_caller);
  }
  increigcg {
    optr->solver = INCREIGCG;
    if(myverbose) printf("  Solver set to INCR-EIG-CG line %d operator %d\n", line_of_file, current_operator);
    BEGIN(name_caller);
  }
  DDalphaAMG {
#ifdef DDalphaAMG
    optr->solver = MG;
    if(myverbose) printf("  Solver set to DDalphaAMG line %d operator %d\n", line_of_file, current_operator);
    BEGIN(name_caller);
#else
    printf("ERROR line %d operator %d: DDalphaAMG library not included\n", line_of_file, current_operator);
    exit(1);
#endif
  }
}

<CSWSOLVER>{
  bicgstab {
    optr->solver = BICGSTAB;
    if(myverbose) printf("  Solver set to BiCGstab line %d operator %d\n", line_of_file, current_operator);
    BEGIN(name_caller);
  }
  mixedbicgstab {
    optr->solver=MIXEDBICGSTAB;
    if(myverbose) printf("  Solver set to MixedBiCGstab line %d operator %d\n", line_of_file, current_operator);
    BEGIN(name_caller);
  }
  increigcg {
    optr->solver = INCREIGCG;
    if(myverbose) printf("  Solver set to INCR-EIG-CG line %d operator %d\n", line_of_file, current_operator);
    BEGIN(name_caller);
  }
  dflgcr {
    optr->solver=DFLGCR;
    g_dflgcr_flag = 1;
    if(myverbose) printf("  Solver set to DFL-GCR line %d operator %d\n", line_of_file, current_operator);
    BEGIN(name_caller);
  }
  dflfgmres {
    optr->solver=DFLFGMRES;
    g_dflgcr_flag = 1;
    if(myverbose) printf("  Solver set to DFL-FGMRES line %d operator %d\n", line_of_file, current_operator);
    BEGIN(name_caller);
  }
  mixedcg {
    optr->solver=MIXEDCG;
    if(myverbose) printf("  Solver set to MixedCG line %d operator %d\n", line_of_file, current_operator);
    BEGIN(name_caller);
  }
  DDalphaAMG {
#ifdef DDalphaAMG
    optr->solver = MG;
    if(myverbose) printf("  Solver set to DDalphaAMG line %d operator %d\n", line_of_file, current_operator);
    BEGIN(name_caller);
#else
    printf("ERROR line %d operator %d: DDalphaAMG library not included\n", line_of_file, current_operator);
    exit(1);
#endif
  }
}

<OVSOLVER>{
  sumr {
    optr->solver = SUMR;
    if(myverbose) printf("  Solver set to SUMR in line %d operator %d\n", line_of_file, current_operator);
    BEGIN(OVERLAPOP);
  }
  cg {
    optr->solver = CG;
    if(myverbose) printf("  Solver set to CG in line %d operator %d\n", line_of_file, current_operator);
    BEGIN(OVERLAPOP);
  }
}

<INITEXTERNALINVERTER>{TYPE} {
  if(strcmp(yytext, "QPHIX")==0) {
#ifdef TM_USE_QPHIX
    if(myverbose) printf("Setting QPHIX external inverter parameters line %d\n", line_of_file);
    BEGIN(QPHIXINVERTER);
#else
    printf("tmLQCD built withtout QPHIX support. Exiting!\n");
    exit(1);
#endif
  } else if ( strcmp(yytext, "QUDA")==0) {
#ifdef TM_USE_QUDA
    if(myverbose) printf("Setting QUDA external inverter parameters line %d\n", line_of_file);
    BEGIN(QUDAINVERTER);
#else
    printf("tmLQCD built withtout QUDA support. Exiting!\n");
    exit(1);
#endif
  } else {
    printf("Error on line %d, inverter '%s' is unknown. Exiting!\n", line_of_file, yytext);
    exit(1);
  }
}
<QUDAINVERTER>{
  {SPC}*FermionBC{EQL}theta {
    quda_input.fermionbc = TM_QUDA_THETABC;
    if(myverbose) printf("  Force theta fermionic temporal boundary conditions for QUDA inversions, line %d\n", line_of_file);
  }
  {SPC}*FermionBC{EQL}apbc {
    quda_input.fermionbc = TM_QUDA_APBC;
    if(myverbose) printf("  Force anti-periodic temporal fermionic boundary conditions for QUDA inversions, line %d\n", line_of_file);
  }
  {SPC}*FermionBC{EQL}pbc {
    quda_input.fermionbc = TM_QUDA_PBC;
    if(myverbose) printf("  Force periodic temporal fermionic boundary conditions for QUDA inversions, line %d\n", line_of_file);
  }
  {SPC}*Pipeline{EQL}{DIGIT}+ {
    sscanf(yytext, " %[a-zA-Z] = %d", name, &a);
    quda_input.pipeline = a;
    if(myverbose) printf("  Pipeline length set to %d line %d\n", quda_input.pipeline, line_of_file);
  }
  {SPC}*MGNumberOfLevels{EQL}{DIGIT}+ {
    sscanf(yytext, " %[a-zA-Z] = %d", name, &a);
    quda_input.mg_n_level = a;
    if(myverbose) printf("  MGNumberOfLevels set to %d line %d\n", quda_input.mg_n_level, line_of_file);
  }
  {SPC}*MGVerbosity{EQL}{STRLIST} {
    parse_quda_mg_verbosity_par_array(yytext, &(quda_input.mg_verbosity[0]));
  }
  {SPC}*MGNumberOfVectors{EQL}{STRLIST} {
    parse_quda_mg_int_par_array(yytext, &(quda_input.mg_n_vec[0]) );
  }
  {SPC}*MGBlockSizesX{EQL}{STRLIST} {
    // note the dimensions here are in QUDA order (XYZT!)
    parse_quda_mg_blocksize(yytext, 0);
  }
  {SPC}*MGBlockSizesY{EQL}{STRLIST} {
    parse_quda_mg_blocksize(yytext, 1);
  }
  {SPC}*MGBlockSizesZ{EQL}{STRLIST} {
    parse_quda_mg_blocksize(yytext, 2);
  }
  {SPC}*MGBlockSizesT{EQL}{STRLIST} {
    parse_quda_mg_blocksize(yytext, 3);
  }
  {SPC}*MGCoarseMuFactor{EQL}{STRLIST} {
    parse_quda_mg_dbl_par_array(yytext, &(quda_input.mg_mu_factor[0]));
  }
  {SPC}*MGCoarseSolverType{EQL}{STRLIST} {
    char paramname[100];
    char token[100];
    char error_message[200];
    int list_end = 0;
    int level = 0;
    char * input_copy = (char*)NULL;
    strlist_tokenize(yytext, input_copy, paramname, 100);
    strlist_next_token(&list_end, token, 100);
    while( list_end != 1 ){
      if( level >= QUDA_MAX_MG_LEVEL ){
        snprintf(error_message, 200, "Exceeded maximum number of levels (%d) parsing %s!\n", QUDA_MAX_MG_LEVEL, paramname);
        yy_fatal_error(error_message);
      }
      /* Make token to lower case */
      for(int i = 0; token[i]; i++){
        token[i] = tolower(token[i]);
      }
      if( strcmp(token, "cagcr") == 0 ){
        quda_input.mg_coarse_solver_type[level] = QUDA_CA_GCR_INVERTER;
      }
      else if( strcmp(token, "gcr") == 0 ){
        quda_input.mg_coarse_solver_type[level] = QUDA_GCR_INVERTER;
      }
      else if( strcmp(token, "mr") == 0){
        quda_input.mg_coarse_solver_type[level] = QUDA_MR_INVERTER;
      }
      else {
        snprintf(error_message, 200, "Unknown inverter %s passed in the parsing of %s!\n", token, paramname);
        yy_fatal_error(error_message);
      }
      if(myverbose) printf("  %s on level %d set to %s line %d\n", paramname, level, token, line_of_file);
      level++;
      strlist_next_token(&list_end, token, 100);
    }
    free(input_copy);
  }
  {SPC}*MGSmootherType{EQL}{STRLIST} {
    char paramname[100];
    char token[100];
    char error_message[200];
    int list_end = 0;
    int level = 0;
    char * input_copy = (char*)NULL;
    strlist_tokenize(yytext, input_copy, paramname, 100);
    strlist_next_token(&list_end, token, 100);
    while( list_end != 1 ){
      if( level >= QUDA_MAX_MG_LEVEL ){
        snprintf(error_message, 200, "Exceeded maximum number of levels (%d) parsing %s!\n", QUDA_MAX_MG_LEVEL, paramname);
        yy_fatal_error(error_message);
      }
      /* Make token to lower case */
      for(int i = 0; token[i]; i++){
        token[i] = tolower(token[i]);
      }
      if( strcmp(token, "cagcr") == 0 ){
        quda_input.mg_smoother_type[level] = QUDA_CA_GCR_INVERTER;
      }
      else if( strcmp(token, "gcr") == 0 ){
        quda_input.mg_smoother_type[level] = QUDA_GCR_INVERTER;
      }
      else if( strcmp(token, "mr") == 0){
        quda_input.mg_smoother_type[level] = QUDA_MR_INVERTER;
      }
      else {
        snprintf(error_message, 200, "Unknown inverter %s passed in the parsing of %s!\n", token, paramname);
        yy_fatal_error(error_message);
      }
      if(myverbose) printf("  %s on level %d set to %s line %d\n", paramname, level, token, line_of_file);
      level++;
      strlist_next_token(&list_end, token, 100);
    }
    free(input_copy);
  }
  {SPC}*MGSetupSolver{EQL}cg {
    quda_input.mg_setup_inv_type = QUDA_CG_INVERTER;
    if(myverbose) printf("  MGSetupSolver set to CG line %d\n", line_of_file);
  }
  {SPC}*MGSetupSolver{EQL}bicgstab {
    quda_input.mg_setup_inv_type = QUDA_BICGSTAB_INVERTER;
    if(myverbose) printf("  MGSetupSolver set to BiCGstab line %d\n", line_of_file);
  }
  {SPC}*MGSetupSolverTolerance{EQL}{STRLIST} {
    parse_quda_mg_dbl_par_array(yytext, &(quda_input.mg_setup_tol[0]));
  }
  {SPC}*MGSetupMaxSolverIterations{EQL}{STRLIST} {
    parse_quda_mg_int_par_array(yytext, &(quda_input.mg_setup_maxiter[0]));
  }
  {SPC}*MGCoarseSolverTolerance{EQL}{STRLIST} {
    parse_quda_mg_dbl_par_array(yytext, &(quda_input.mg_coarse_solver_tol[0]));
  }
  {SPC}*MGCoarseMaxSolverIterations{EQL}{STRLIST} {
    parse_quda_mg_int_par_array(yytext, &(quda_input.mg_coarse_solver_maxiter[0]));
  }
  {SPC}*MGSmootherTolerance{EQL}{STRLIST} {
    parse_quda_mg_dbl_par_array(yytext, &(quda_input.mg_smoother_tol[0]));
  }
  {SPC}*MGSmootherPreIterations{EQL}{STRLIST} {
    parse_quda_mg_int_par_array(yytext, &(quda_input.mg_nu_pre[0]));
  }
  {SPC}*MGSmootherPostIterations{EQL}{STRLIST} {
    parse_quda_mg_int_par_array(yytext, &(quda_input.mg_nu_post[0]));
  }
  {SPC}*MGOverUnderRelaxationFactor{EQL}{STRLIST} {
    parse_quda_mg_dbl_par_array(yytext, &(quda_input.mg_omega[0]));
  }
  {SPC}*MGRunVerify{EQL}yes {
    quda_input.mg_run_verify = QUDA_BOOLEAN_YES;
    if(myverbose) printf("  MGRunVerify set to YES in line %d\n", line_of_file);
  }
  {SPC}*MGRunVerify{EQL}no {
    quda_input.mg_run_verify = QUDA_BOOLEAN_NO;
    if(myverbose) printf("  MGRunVerify set to NO in line %d\n", line_of_file);
  }
  {SPC}*MGRunLowModeCheck{EQL}yes {
    quda_input.mg_run_low_mode_check = QUDA_BOOLEAN_YES;
    if(myverbose) printf("  MGRunLowModeCheck set to YES in line %d\n", line_of_file);
  }
  {SPC}*MGRunLowModeCheck{EQL}no {
    quda_input.mg_run_low_mode_check = QUDA_BOOLEAN_NO;
    if(myverbose) printf("  MGRunLowModeCheck set to NO in line %d\n", line_of_file);
  }
  {SPC}*MGRunObliqueProjectionCheck{EQL}yes {
    quda_input.mg_run_oblique_proj_check = QUDA_BOOLEAN_YES;
    if(myverbose) printf("  MGRunObliqueProjectionCheck set to YES in line %d\n", line_of_file);
  }
  {SPC}*MGRunObliqueProjectionCheck{EQL}no {
    quda_input.mg_run_oblique_proj_check = QUDA_BOOLEAN_NO;
    if(myverbose) printf("  MGRunObliqueProjectionCheck set to NO in line %d\n", line_of_file);
  }
  {SPC}*MGEnableSizeThreeBlocks{EQL}yes {
    quda_input.mg_enable_size_three_blocks = 1;
    if(myverbose) printf("  MGEnableSizeThreeBlocks set to YES in line %d\n", line_of_file);
  }
  {SPC}*MGEnableSizeThreeBlocks{EQL}no {
    quda_input.mg_enable_size_three_blocks = 0;
    if(myverbose) printf("  MGEnableSizeThreeBlocks set to NO in line %d\n", line_of_file);
  }
  {SPC}*MGResetSetupThreshold{EQL}{FLT}+ {
    sscanf(yytext, " %[a-zA-Z] = %lf", name, &c);
    quda_input.mg_reset_setup_threshold = c;
    if(myverbose) printf("  MGResetSetupThreshold set to %f line %d\n", quda_input.mg_reset_setup_threshold, line_of_file);
  }
  {SPC}*MGUseEigSolver{EQL}{STRLIST} {
    parse_quda_mg_bool_par_array(yytext, &(quda_input.mg_use_eig_solver[0]));
  }
  {SPC}*MGEigSolverType{EQL}{STRLIST} {
    char paramname[100];
    char token[100];
    char error_message[200];
    int list_end = 0;
    int level = 0;
    char * input_copy = (char*)NULL;
    strlist_tokenize(yytext, input_copy, paramname, 100);
    strlist_next_token(&list_end, token, 100);
    while( list_end != 1 ){
      if( level >= QUDA_MAX_MG_LEVEL ){
        snprintf(error_message, 200, "Exceeded maximum number of levels (%d) parsing '%s'!\n", QUDA_MAX_MG_LEVEL, paramname);
        yy_fatal_error(error_message);
      }
      /* Make token to lower case */
      for(int i = 0; token[i]; i++){
        token[i] = tolower(token[i]);
      }
      if( strcmp(token, "tr_lanczos") == 0 ){
        quda_input.mg_eig_type[level] = QUDA_EIG_TR_LANCZOS;
      }
      else if( strcmp(token, "ir_lanczos") == 0 ){
        quda_input.mg_eig_type[level] = QUDA_EIG_IR_LANCZOS;
      }
      else if( strcmp(token, "ir_arnoldi") == 0){
        quda_input.mg_eig_type[level] = QUDA_EIG_IR_ARNOLDI;
      }
      else {
        snprintf(error_message, 200, "Unknown type '%s' passed in the parsing of '%s'!\n", token, paramname);
        yy_fatal_error(error_message);
      }
      if(myverbose) printf("  MGEigSolverType on level %d set to '%s' line %d\n", level, token, line_of_file);
      level++;
      strlist_next_token(&list_end, token, 100);
    }
    free(input_copy);
  }
  {SPC}*MGEigSolverSpectrum{EQL}{STRLIST} {
    char paramname[100];
    char token[100];
    char error_message[200];
    int list_end = 0;
    int level = 0;
    char * input_copy = (char*)NULL;
    strlist_tokenize(yytext, input_copy, paramname, 100);
    strlist_next_token(&list_end, token, 100);
    while( list_end != 1 ){
      if( level >= QUDA_MAX_MG_LEVEL ){
        snprintf(error_message, 200, "Exceeded maximum number of levels (%d) parsing '%s'!\n", QUDA_MAX_MG_LEVEL, paramname);
        yy_fatal_error(error_message);
      }
      /* Make token to lower case */
      for(int i = 0; token[i]; i++){
        token[i] = tolower(token[i]);
      }
      if( strcmp(token, "smallest_real") == 0 ){
        quda_input.mg_eig_spectrum[level] = QUDA_SPECTRUM_SR_EIG;
      }
      else if( strcmp(token, "largest_real") == 0 ){
        quda_input.mg_eig_spectrum[level] = QUDA_SPECTRUM_LR_EIG;
      }
      else if( strcmp(token, "smallest_modulus") == 0){
        quda_input.mg_eig_spectrum[level] = QUDA_SPECTRUM_SM_EIG;
      }
      else if( strcmp(token, "largest_modulus") == 0){
        quda_input.mg_eig_spectrum[level] = QUDA_SPECTRUM_LM_EIG;
      }
      else if( strcmp(token, "smallest_imaginary") == 0){
        quda_input.mg_eig_spectrum[level] = QUDA_SPECTRUM_SI_EIG;
      }
      else if( strcmp(token, "largest_imaginary") == 0){
        quda_input.mg_eig_spectrum[level] = QUDA_SPECTRUM_LI_EIG;
      }
      else {
        snprintf(error_message, 200, "Unknown spectrum type '%s' passed in the parsing of '%s'!\n", token, paramname);
        yy_fatal_error(error_message);
      }
      if(myverbose) printf("  %s on level %d set to '%s' line %d\n", paramname, level, token, line_of_file);
      level++;
      strlist_next_token(&list_end, token, 100);
    }
    free(input_copy);
  }
  {SPC}*MGEigSolverCoarseGuess{EQL}yes {
    quda_input.mg_coarse_guess = QUDA_BOOLEAN_YES;
    if(myverbose) printf("  MGEigSolverCoarseGuess set to YES in line %d\n", line_of_file);
  }
  {SPC}*MGEigSolverCoarseGuess{EQL}no {
    quda_input.mg_coarse_guess = QUDA_BOOLEAN_NO;
    if(myverbose) printf("  MGEigSolverCoarseGuess set to NO in line %d\n", line_of_file);
  }
  {SPC}*MGEigSolverNumberOfVectors{EQL}{STRLIST} {
    parse_quda_mg_int_par_array(yytext, &(quda_input.mg_eig_nEv[0]));
  }
  {SPC}*MGEigSolverKrylovSubspaceSize{EQL}{STRLIST} {
    parse_quda_mg_int_par_array(yytext, &(quda_input.mg_eig_nKr[0]));
  }
  {SPC}*MGEigSolverRequireConvergence{EQL}{STRLIST} {
    parse_quda_mg_bool_par_array(yytext, &(quda_input.mg_eig_require_convergence[0]));
  }
  {SPC}*MGEigSolverCheckConvergenceInterval{EQL}{STRLIST} {
    parse_quda_mg_int_par_array(yytext, &(quda_input.mg_eig_check_interval[0]));
  }
  {SPC}*MGEigSolverMaxRestarts{EQL}{STRLIST} {
    parse_quda_mg_int_par_array(yytext, &(quda_input.mg_eig_max_restarts[0]));
  }
  {SPC}*MGEigSolverTolerance{EQL}{STRLIST} {
    parse_quda_mg_dbl_par_array(yytext, &(quda_input.mg_eig_tol[0]));
  }
  {SPC}*MGEigSolverUseNormOp{EQL}{STRLIST} {
    parse_quda_mg_bool_par_array(yytext, &(quda_input.mg_eig_use_normop[0]));
  }
  {SPC}*MGEigSolverUseDagger{EQL}{STRLIST} {
    parse_quda_mg_bool_par_array(yytext, &(quda_input.mg_eig_use_dagger[0]));
  }
  {SPC}*MGEigSolverUsePolynomialAcceleration{EQL}{STRLIST} {
    parse_quda_mg_bool_par_array(yytext, &(quda_input.mg_eig_use_poly_acc[0]));
  }
  {SPC}*MGEigSolverPolynomialDegree{EQL}{STRLIST} {
    parse_quda_mg_int_par_array(yytext, &(quda_input.mg_eig_poly_deg[0]));
  }
  {SPC}*MGEigSolverPolyMin{EQL}{STRLIST} {
    parse_quda_mg_dbl_par_array(yytext, &(quda_input.mg_eig_amin[0]));
  }
  {SPC}*MGEigSolverPolyMax{EQL}{STRLIST} {
    parse_quda_mg_dbl_par_array(yytext, &(quda_input.mg_eig_amax[0]));
  }
  ^EndExternalInverter{SPC}*  { 
    if(myverbose) printf("QUDA external inverter parameters parsed on line %d\n\n", line_of_file);
    BEGIN(0);
  }
}
<QPHIXINVERTER>{
  {SPC}*NCores{EQL}{DIGIT}+ {
    sscanf(yytext, " %[a-zA-Z] = %d", name, &a);
    qphix_input.NCores = a;
    if(myverbose) printf("  Ncores = %d, line %d\n", qphix_input.NCores, line_of_file );
  }
  {SPC}*By{EQL}{DIGIT}+ {
    sscanf(yytext, " %[a-zA-Z] = %d", name, &a);
    qphix_input.By = a;
    if(myverbose) printf("  By = %d, line %d\n", qphix_input.By, line_of_file );
  }
  {SPC}*Bz{EQL}{DIGIT}+ {
    sscanf(yytext, " %[a-zA-Z] = %d", name, &a);
    qphix_input.Bz = a;
    if(myverbose) printf("  Bz = %d, line %d\n", qphix_input.Bz, line_of_file );
  }
  {SPC}*Sy{EQL}{DIGIT}+ {
    sscanf(yytext, " %[a-zA-Z] = %d", name, &a);
    qphix_input.Sy = a;
    if(myverbose) printf("  Sy = %d, line %d\n", qphix_input.Sy, line_of_file );
  }
  {SPC}*Sz{EQL}{DIGIT}+ {
    sscanf(yytext, " %[a-zA-Z] = %d", name, &a);
    qphix_input.Sz = a;
    if(myverbose) printf("  Sz = %d, line %d\n", qphix_input.Sz, line_of_file );
  }
  {SPC}*PadXY{EQL}{DIGIT}+ {
    sscanf(yytext, " %[a-zA-Z] = %d", name, &a);
    qphix_input.PadXY = a;
    if(myverbose) printf("  PadXY = %d, line %d\n", qphix_input.PadXY, line_of_file );
  }
  {SPC}*PadXYZ{EQL}{DIGIT}+ {
    sscanf(yytext, " %[a-zA-Z] = %d", name, &a);
    qphix_input.PadXYZ = a;
    if(myverbose) printf("  PadXYZ = %d, line %d\n", qphix_input.PadXYZ, line_of_file );
  }
  {SPC}*MinCt{EQL}{DIGIT}+ {
    sscanf(yytext, " %[a-zA-Z] = %d", name, &a);
    qphix_input.MinCt = a;
    if(myverbose) printf("  MinCt = %d, line %d\n", qphix_input.MinCt, line_of_file );
  }
  ^EndExternalInverter{SPC}*  {
    if(myverbose) printf("QPHIX external inverter parameters parsed on line %d\n\n", line_of_file);
    BEGIN(0);
  }
}

<INITMONOMIAL>{TYPE} {
  current_monomial++;
  mnl = &monomial_list[current_monomial];
  mnl->id = current_monomial;
  if(strcmp(yytext, "DET")==0) {
    mnl->type = DET;
    strcpy((*mnl).name, "DET");
  }
  else if(strcmp(yytext, "CLOVERDET")==0) {
    mnl->type = CLOVERDET;
    strcpy((*mnl).name, "CLOVERDET");
  }
  else if(strcmp(yytext, "CLOVERDETRATIO")==0) {
    mnl->type = CLOVERDETRATIO;
    strcpy((*mnl).name, "CLOVERDETRATIO");
  }
  else if(strcmp(yytext, "CLOVERDETRATIORW")==0) {
    mnl->type = CLOVERDETRATIORW;
    strcpy((*mnl).name, "CLOVERDETRATIORW");
  }
  else if(strcmp(yytext, "DETRATIO")==0) {
    mnl->type = DETRATIO;
    strcpy((*mnl).name, "DETRATIO");
  }
  else if(strcmp(yytext, "NDDETRATIO")==0) {
    mnl->type = NDDETRATIO;
    strcpy((*mnl).name, "NDDETRATIO");
    g_running_phmc = 1;
  }
  else if(strcmp(yytext, "NDCLOVERDETRATIO")==0) {
    mnl->type = NDCLOVERDETRATIO;
    strcpy((*mnl).name, "NDDCLOVERETRATIO");
    g_running_phmc = 1;
  }
  else if(strcmp(yytext, "NDPOLY")==0) {
    mnl->type = NDPOLY;
    strcpy((*mnl).name, "NDPOLY");
    g_running_phmc = 1;
  }
  else if(strcmp(yytext, "NDRAT")==0) {
    mnl->type = NDRAT;
    strcpy((*mnl).name, "NDRAT");
    g_running_phmc = 1;
  }
  else if(strcmp(yytext, "RAT")==0) {
    mnl->type = RAT;
    strcpy((*mnl).name, "RAT");
    g_running_phmc = 1;
  }
  else if(strcmp(yytext, "CLOVERRAT")==0) {
    mnl->type = CLOVERRAT;
    strcpy((*mnl).name, "CLOVERRAT");
    g_running_phmc = 1;
  }
  else if(strcmp(yytext, "RATCOR")==0) {
    mnl->type = RATCOR;
    strcpy((*mnl).name, "RATCOR");
    g_running_phmc = 1;
  }
  else if(strcmp(yytext, "CLOVERRATCOR")==0) {
    mnl->type = CLOVERRATCOR;
    strcpy((*mnl).name, "CLOVERRATCOR");
    g_running_phmc = 1;
  }
  else if(strcmp(yytext, "NDCLOVERRAT")==0) {
    mnl->type = NDCLOVERRAT;
    strcpy((*mnl).name, "NDCLOVERRAT");
    g_running_phmc = 1;
  }
  else if(strcmp(yytext, "NDRATCOR")==0) {
    mnl->type = NDRATCOR;
    strcpy((*mnl).name, "NDRATCOR");
    g_running_phmc = 1;
  }
  else if(strcmp(yytext, "NDCLOVERRATCOR")==0) {
    mnl->type = NDCLOVERRATCOR;
    strcpy((*mnl).name, "NDCLOVERRATCOR");
    g_running_phmc = 1;
  }
  else if(strcmp(yytext, "NDCLOVER")==0) {
    mnl->type = NDCLOVER;
    strcpy((*mnl).name, "NDCLOVER");
    g_running_phmc = 1;
  }
  else if(strcmp(yytext, "POLY")==0) {
    mnl->type = POLY;
    strcpy((*mnl).name, "POLY");
  }
  else if(strcmp(yytext, "POLYDETRATIO")==0) {
    mnl->type = POLYDETRATIO;
    strcpy((*mnl).name, "POLYDETRATIO");
  }
  else if(strcmp(yytext, "GAUGE")==0) {
    mnl->type = GAUGE;
    mnl->gtype = 3;
    strcpy((*mnl).name, "GAUGE");
  }
  else {
    fprintf(stderr, "Unknown monomial type %s in line %d\n", yytext, line_of_file);
    exit(1);
  }
  if(!reread) {
    if(add_monomial(mnl->type) < 0) {
      fprintf(stderr, "Something went wrong in adding monomials\nAborting...!\n");
      exit(1);
    }
  }
  if(myverbose) printf("initialising monomial with type %s %d line %d\n", yytext, mnl->type, line_of_file);
  if(myverbose) printf("monomial has id %d\n", current_monomial);

  if(mnl->type == GAUGE) BEGIN(GAUGEMONOMIAL);
  else if(mnl->type == NDPOLY) BEGIN(NDPOLYMONOMIAL);
  else if(mnl->type == NDCLOVER) BEGIN(CLPOLYMONOMIAL);
  else if(mnl->type == NDRAT) BEGIN(NDRATMONOMIAL);
  else if(mnl->type == RAT) BEGIN(RATMONOMIAL);
  else if(mnl->type == NDCLOVERRAT) BEGIN(NDCLRATMONOMIAL);
  else if(mnl->type == NDDETRATIO) BEGIN(NDDETRATMONOMIAL);
  else if(mnl->type == CLOVERRAT) BEGIN(CLRATMONOMIAL);
  else if(mnl->type == NDRATCOR) BEGIN(NDRATCORMONOMIAL);
  else if(mnl->type == RATCOR) BEGIN(RATCORMONOMIAL);
  else if(mnl->type == NDCLOVERRATCOR) BEGIN(NDCLRATCORMONOMIAL);
  else if(mnl->type == NDDETRATIO) BEGIN(NDDETRATMONOMIAL);
  else if(mnl->type == NDCLOVERDETRATIO) BEGIN(NDCLDETRATMONOMIAL);
  else if(mnl->type == CLOVERRATCOR) BEGIN(CLRATCORMONOMIAL);
  else if(mnl->type == POLY || mnl->type == POLYDETRATIO)  {
          fprintf(stderr,"starting to parse poly(detratio) monomial\n");
          BEGIN(POLYMONOMIAL); 
  }
  else if(mnl->type == CLOVERDET) BEGIN(CLDETMONOMIAL);
  else if(mnl->type == CLOVERDETRATIO) BEGIN(CLDETRATMONOMIAL);
  else if(mnl->type == CLOVERDETRATIORW) BEGIN(CLDETRATRWMONOMIAL);
  else BEGIN(DETMONOMIAL);
}



<DETMONOMIAL,GAUGEMONOMIAL,NDPOLYMONOMIAL,POLYMONOMIAL,CLDETMONOMIAL,CLDETRATMONOMIAL,CLDETRATRWMONOMIAL,CLPOLYMONOMIAL,NDRATMONOMIAL,NDRATCORMONOMIAL,NDCLRATMONOMIAL,NDCLRATCORMONOMIAL,RATMONOMIAL,CLRATMONOMIAL,RATCORMONOMIAL,CLRATCORMONOMIAL,NDDETRATMONOMIAL,NDCLDETRATMONOMIAL>{
  {SPC}*Timescale{EQL}{DIGIT}+ {
    if(mnl->type == NDDETRATIO || mnl->type == NDCLOVERDETRATIO) {
      mnl->timescale = -5;
      if(myverbose) printf("  timescales set to %d line %d monomial %d since NDDETRATIO is not for MD evolution\n", a, line_of_file, current_monomial);
    }
    else {
      sscanf(yytext, " %[a-zA-Z] = %d", name, &a);
      mnl->timescale = a;
      if(myverbose) printf("  timescales set to %d line %d monomial %d\n", a, line_of_file, current_monomial);
    }
  }
  {SPC}*Name{EQL} {
    name_caller = YY_START;
    BEGIN(MNAME);
  }
  ^EndMonomial{SPC}*  {
    if(myverbose) printf("monomial %d parsed line %d\n\n", current_monomial, line_of_file);
    BEGIN(0);
  }
}

<CLDETMONOMIAL,CLDETRATMONOMIAL,CLDETRATRWMONOMIAL,CLPOLYMONOMIAL,NDCLRATMONOMIAL,NDCLRATCORMONOMIAL,CLRATMONOMIAL,CLRATCORMONOMIAL,NDCLDETRATMONOMIAL>{
  {SPC}*CSW{EQL}{FLT} {
    sscanf(yytext, " %[a-zA-Z] = %lf", name, &c);
    mnl->c_sw = c;
    if(myverbose) printf("  CSW set to %f line %d monomial %d\n", c, line_of_file, current_monomial);
  }
}

<CLDETMONOMIAL,CLDETRATMONOMIAL,CLDETRATRWMONOMIAL>{
  {SPC}*rho{EQL}{FLT} {
    sscanf(yytext, " %[a-zA-Z] = %lf", name, &c);
    mnl->rho = c;
    if(myverbose) printf("  mass shift rho set to %f line %d monomial %d\n", c, line_of_file, current_monomial);
  }
}

<CLDETRATMONOMIAL,CLDETRATRWMONOMIAL>{
  {SPC}*rho2{EQL}{FLT} {
    sscanf(yytext, " %[a-zA-Z]2 = %lf", name, &c);
    mnl->rho2 = c;
    if(myverbose) printf("  mass shift rho2 set to %f line %d monomial %d\n", c, line_of_file, current_monomial);
  }
}

<DETMONOMIAL,POLYMONOMIAL,NDPOLYMONOMIAL,CLPOLYMONOMIAL,NDRATMONOMIAL,NDRATCORMONOMIAL,NDCLRATMONOMIAL,NDCLRATCORMONOMIAL,CLDETMONOMIAL,CLDETRATMONOMIAL,CLDETRATRWMONOMIAL,RATMONOMIAL,CLRATMONOMIAL,RATCORMONOMIAL,CLRATCORMONOMIAL,NDDETRATMONOMIAL,NDCLDETRATMONOMIAL>{
  {SPC}*Kappa{EQL}{FLT} {
    sscanf(yytext, " %[a-zA-Z] = %lf", name, &c);
    mnl->kappa = c;
    if(myverbose) printf("  Kappa set to %f line %d monomial %d\n", c, line_of_file, current_monomial);
  }
}

<DETMONOMIAL,POLYMONOMIAL,CLDETRATRWMONOMIAL>{
  {SPC}*2KappaMu2{EQL}{FLT} {
    sscanf(yytext, " %[2a-zA-Z] = %lf", name, &c);
    mnl->mu2 = c;
    if(myverbose) printf("  2KappaMu2 set to %f line %d monomial %d\n", c, line_of_file, current_monomial);
  }
  {SPC}*Kappa2{EQL}{FLT} {
    sscanf(yytext, " %[2a-zA-Z] = %lf", name, &c);
    mnl->kappa2 = c;
    if(myverbose) printf("  kappa2 set to %f line %d monomial %d\n", c, line_of_file, current_monomial);
  }
}

<NDDETRATMONOMIAL,NDCLDETRATMONOMIAL>{
  {SPC}*Kappa2{EQL}{FLT} {
    sscanf(yytext, " %[2a-zA-Z] = %lf", name, &c);
    mnl->kappa2 = c;
    if(myverbose) printf("  kappa2 set to %f line %d monomial %d\n", c, line_of_file, current_monomial);
  }
}

<NDCLRATMONOMIAL,CLRATMONOMIAL>{
  {SPC}*AddTrLog{EQL}yes {
    mnl->trlog = 1;
    if(myverbose) printf("  Added trlog monomial in line %d to monomial %d\n", line_of_file, current_monomial);
  }
  {SPC}*AddTrLog{EQL}no {
    mnl->trlog = 0;
    if(myverbose) printf("  No trlog monomial (default) in line %d for monomial %d\n", line_of_file, current_monomial);
  }
}

<NDPOLYMONOMIAL,CLPOLYMONOMIAL,NDRATMONOMIAL,NDRATCORMONOMIAL,NDCLRATMONOMIAL,NDCLRATCORMONOMIAL,NDDETRATMONOMIAL,NDCLDETRATMONOMIAL>{
  {SPC}*2KappaMubar{EQL}{FLT} {
    sscanf(yytext, " %[2a-zA-Z] = %lf", name, &c);
    mnl->mubar = c;
    if(myverbose) printf("  2KappaMubar set to %f line %d monomial %d\n", c, line_of_file, current_monomial);
  }
  {SPC}*2KappaEpsbar{EQL}{FLT} {
    sscanf(yytext, " %[2a-zA-Z] = %lf", name, &c);
    mnl->epsbar = c;
    if(myverbose) printf("  2KappaEpsbar set to %f line %d monomial %d\n", c, line_of_file, current_monomial);
  }
}

<NDDETRATMONOMIAL,NDCLDETRATMONOMIAL>{
  {SPC}*2KappaMubar2{EQL}{FLT} {
    sscanf(yytext, " %[2a-zA-Z] = %lf", name, &c);
    mnl->mubar2 = c;
    if(myverbose) printf("  2KappaMubar2 set to %f line %d monomial %d\n", c, line_of_file, current_monomial);
  }
  {SPC}*2KappaEpsbar2{EQL}{FLT} {
    sscanf(yytext, " %[2a-zA-Z] = %lf", name, &c);
    mnl->epsbar2 = c;
    if(myverbose) printf("  2KappaEpsbar2 set to %f line %d monomial %d\n", c, line_of_file, current_monomial);
  }
}
<DETMONOMIAL,CLDETMONOMIAL,CLDETRATMONOMIAL,CLDETRATRWMONOMIAL,NDRATMONOMIAL,NDRATCORMONOMIAL,NDCLRATMONOMIAL,NDCLRATCORMONOMIAL,RATMONOMIAL,RATCORMONOMIAL,CLRATMONOMIAL,CLRATCORMONOMIAL,NDDETRATMONOMIAL,NDCLDETRATMONOMIAL>{
  {SPC}*UseExternalInverter{EQL}quda {
    if(myverbose) printf("  Use Quda inverter line %d monomial %d\n", line_of_file, current_monomial);
    mnl->solver_params.external_inverter = QUDA_INVERTER;
  }
  {SPC}*UseExternalInverter{EQL}qphix {
    if(myverbose) printf("  Use QPhiX inverter line %d monomial %d\n", line_of_file, current_monomial);
    mnl->solver_params.external_inverter = QPHIX_INVERTER;
  }
  {SPC}*UseExternalInverter{EQL}no {
    if(myverbose) printf("  Use QPhiX inverter line %d monomial %d\n", line_of_file, current_monomial);
    mnl->solver_params.external_inverter = NO_EXT_INV;
  }
  {SPC}*UseSloppyPrecision{EQL}yes {
    if(myverbose) printf("  Use use sloppy precision (single) in the inverter (if supported by the inverter) line %d monomial %d\n", line_of_file, current_monomial);
    mnl->solver_params.sloppy_precision = SLOPPY_SINGLE;
  }
  {SPC}*UseSloppyPrecision{EQL}float {
    if(myverbose) printf("  Use use sloppy precision (single) in the inverter (if supported by the inverter) line %d monomial %d\n", line_of_file, current_monomial);
    mnl->solver_params.sloppy_precision = SLOPPY_SINGLE;
  }
  {SPC}*UseSloppyPrecision{EQL}single {
    if(myverbose) printf("  Use use sloppy precision (single) in the inverter (if supported by the inverter) line %d monomial %d\n", line_of_file, current_monomial);
    mnl->solver_params.sloppy_precision = SLOPPY_SINGLE;
  }
  {SPC}*UseSloppyPrecision{EQL}no {
    if(myverbose) printf("  Use use sloppy precision (single) in the inverter line %d monomial %d\n", line_of_file, current_monomial);
    mnl->solver_params.sloppy_precision = SLOPPY_DOUBLE;
  }
  {SPC}*UseSloppyPrecision{EQL}double {
    if(myverbose) printf("  Use use sloppy precision (single) in the inverter line %d monomial %d\n", line_of_file, current_monomial);
    mnl->solver_params.sloppy_precision = SLOPPY_DOUBLE;
  }
  {SPC}*UseSloppyPrecision{EQL}half {
    if(myverbose) printf("  Use use sloppy precision (half) in the inverter (if supported by the inverter) line %d monomial %d\n", line_of_file, current_monomial);
    mnl->solver_params.sloppy_precision = SLOPPY_HALF;
  }
  {SPC}*UseCompression{EQL}12 {
    if(myverbose) printf("  Use 12 compression in the inverter (if supported) line %d monomial %d\n", line_of_file, current_monomial);
    mnl->solver_params.compression_type = COMPRESSION_12;
  }
  {SPC}*UseCompression{EQL}8 {
    if(myverbose) printf("  Use 8 compression in the inverter (if supported) line %d monomial %d\n", line_of_file, current_monomial);
    mnl->solver_params.compression_type = COMPRESSION_8;
  }
  {SPC}*UseCompression{EQL}18 {
    if(myverbose) printf("  Not using compression in the inverter line %d monomial %d\n", line_of_file, current_monomial);
    mnl->solver_params.compression_type = NO_COMPRESSION;
  }
}

<DETMONOMIAL,POLYMONOMIAL,CLDETMONOMIAL,CLDETRATMONOMIAL,CLDETRATRWMONOMIAL,NDRATMONOMIAL,NDRATCORMONOMIAL,NDCLRATMONOMIAL,NDCLRATCORMONOMIAL,RATMONOMIAL,RATCORMONOMIAL,CLRATMONOMIAL,CLRATCORMONOMIAL,NDDETRATMONOMIAL,NDCLDETRATMONOMIAL>{
  {SPC}*ForcePrecision{EQL}{FLT} {
    sscanf(yytext, " %[a-zA-Z] = %lf",name , &c);
    mnl->forceprec = c;
    if(myverbose) printf("  ForcePrecision set to %e line %d monomial %d\n", c, line_of_file, current_monomial);
  }
  {SPC}*AcceptancePrecision{EQL}{FLT} {
    sscanf(yytext, " %[a-zA-Z] = %lf",name , &c);
    mnl->accprec = c;
    if(myverbose) printf("  AcceptancePrecision set to %e line %d monomial %d\n", c, line_of_file, current_monomial);
  }
  {SPC}*MaxSolverIterations{EQL}{DIGIT}+ {
    sscanf(yytext, " %[a-zA-Z] = %d", name, &a);
    mnl->maxiter = a;
    if(myverbose) printf("  MaxSolverIterations set to %d line %d monomial %d\n", a, line_of_file, current_monomial);
  }
  {SPC}*mcgdelta{EQL}{FLT} {
    sscanf(yytext, " %[a-zA-Z1] = %lf", name, &c);
    (mnl->solver_params).mcg_delta = c;
    if(myverbose) printf("  mcg_delta set to %lf line %d monomial %d\n", c, line_of_file, current_monomial);
  }
}

<NDRATMONOMIAL,NDRATCORMONOMIAL,NDCLRATMONOMIAL,NDCLRATCORMONOMIAL,NDDETRATMONOMIAL,NDCLDETRATMONOMIAL>{
  {SPC}*Solver{EQL} {
   solver_caller=YY_START;
   BEGIN(NDMSOLVER);
  }
}

<RATMONOMIAL,RATCORMONOMIAL,CLRATMONOMIAL,CLRATCORMONOMIAL>{
  {SPC}*Solver{EQL} {
   solver_caller=YY_START;
   BEGIN(RATMSOLVER);
  }
}

<DETMONOMIAL,POLYMONOMIAL,CLDETMONOMIAL,CLDETRATMONOMIAL,CLDETRATRWMONOMIAL,RATMONOMIAL,RATCORMONOMIAL,CLRATMONOMIAL,CLRATCORMONOMIAL>{
  {SPC}*2KappaMu{EQL}{FLT} {
    sscanf(yytext, " %[2a-zA-Z] = %lf", name, &c);
    mnl->mu = c;
    if(myverbose) printf("  2KappaMu set to %f line %d monomial %d\n", c, line_of_file, current_monomial);
  }
  {SPC}*CSGHistory{EQL}{DIGIT}+ {
    sscanf(yytext, " %[a-zA-Z] = %d", name, &a);
    mnl->csg_N = a;
    if(myverbose) printf("  csg history length set to %d line %d monomial %d\n", a, line_of_file, current_monomial);
  }
  {SPC}*CSGHistory2{EQL}{DIGIT}+ {
    sscanf(yytext, " %[a-zA-Z2] = %d", name , &a);
    mnl->csg_N2 = a;
    if(myverbose) printf("  csg history2 length (for bicgstab) set to %d line %d monomial %d\n", 
                       a, line_of_file, current_monomial);
  }
  {SPC}*Solver{EQL} {
   solver_caller=YY_START;
   BEGIN(MSOLVER);
  }
}

<GAUGEMONOMIAL>{
  {SPC}*Type{EQL} BEGIN(GTYPE);
  {SPC}*UseRectangleStaples{EQL}yes {
    mnl->use_rectangles = 1;
    g_dbw2rand = 1;
    if(myverbose) printf("  UseRectangleStaples set to true line %d monomial %d\n", line_of_file, current_monomial);
  }
  {SPC}*UseRectangleStaples{EQL}no {
    mnl->use_rectangles = 0;
    /* g_dbw2rand = 0; */
    mnl->c1 = 0.;
    g_rgi_C1 = 0.;
    if(myverbose) printf("  UseRectangleStaples set to false line %d monomial %d\n", line_of_file, current_monomial);
  }
  {SPC}*Beta{EQL}{FLT} {
    sscanf(yytext, " %[a-zA-Z] = %lf",name , &c);
    mnl->beta = c;
    g_beta = c;
    if(myverbose) printf("  beta set to %e line %d monomial %d\n", c, line_of_file, current_monomial);
  }
  {SPC}*RectangleCoefficient{EQL}{FLT} {
    sscanf(yytext, " %[a-zA-Z] = %lf",name , &c);
    mnl->c1 = c;
    g_rgi_C1 = c;
    if(myverbose) printf("  RectangleCoefficient c1  set to %e line %d monomial %d\n", c, line_of_file, current_monomial);
  }
  {SPC}*Lambda{EQL}{FLT} {
    sscanf(yytext, " %[a-zA-Z] = %lf",name , &c);
    mnl->glambda = c;
    if(myverbose) printf("  Gauge lambda parameter (Wilson plaquette only) set to %e line %d monomial %d\n", c, line_of_file, current_monomial);
  }
}

<NDPOLYMONOMIAL>{
  {SPC}*ExactPolynomial{EQL}yes {
    phmc_exact_poly = 1;
    if(myverbose!=0) printf("  phmc_exact_poly set to true line %d monomial %d\n", line_of_file, current_monomial);
  }
  {SPC}*ExactPolynomial{EQL}no {
    phmc_exact_poly = 0;
    if(myverbose!=0) printf("  phmc_exact_poly set to false line %d monomial %d\n", line_of_file, current_monomial);
  }
}
<NDRATMONOMIAL,NDCLRATMONOMIAL,RATMONOMIAL,CLRATMONOMIAL>{
  {SPC}*Cmin{EQL}{DIGIT}+ {
    sscanf(yytext, " %[a-zA-Z] = %d", name, &a);
    mnl->rat.crange[0] = a;
    if(myverbose!=0) printf("  Coefficient range of rational starts at coefficient %d line %d monomial %d\n", a, line_of_file, current_monomial);
  }
  {SPC}*Cmax{EQL}{DIGIT}+ {
    sscanf(yytext, " %[a-zA-Z] = %d", name, &a);
    mnl->rat.crange[1] = a;
    if(myverbose!=0) printf("  Coefficient range of rational ends at coefficient %d line %d monomial %d\n", a, line_of_file, current_monomial);
  }
}

<NDRATMONOMIAL,NDRATCORMONOMIAL,NDCLRATMONOMIAL,NDCLRATCORMONOMIAL,RATMONOMIAL,RATCORMONOMIAL,CLRATMONOMIAL,CLRATCORMONOMIAL>{
  {SPC}*DegreeOfRational{EQL}{DIGIT}+ {
    sscanf(yytext, " %[a-zA-Z] = %d", name, &a);
    mnl->rat.order = a;
    if(myverbose!=0) printf("  Degree of rational approximation set to %d line %d monomial %d\n", a, line_of_file, current_monomial);
  }
}
<NDPOLYMONOMIAL,CLPOLYMONOMIAL,NDRATMONOMIAL,NDRATCORMONOMIAL,NDCLRATMONOMIAL,NDCLRATCORMONOMIAL,RATMONOMIAL,RATCORMONOMIAL,CLRATMONOMIAL,CLRATCORMONOMIAL>{
  {SPC}*StildeMax{EQL}{FLT} {
    sscanf(yytext, " %[a-zA-Z] = %lf",name , &c);
    stilde_max = c;
    mnl->StildeMax = c;
    mnl->rat.range[1] = c;
    if(myverbose!=0) printf("  Stilde max set to %e line %d monomial %d\n", c, line_of_file, current_monomial);
  }
  {SPC}*StildeMin{EQL}{FLT} {
    sscanf(yytext, " %[a-zA-Z] = %lf",name , &c);
    mnl->StildeMin = c;
    stilde_min = c;
    mnl->rat.range[0] = c;
    if(myverbose!=0) printf("  Stilde min set to %e line %d monomial %d\n", c, line_of_file, current_monomial);
  }
  {SPC}*ComputeEVFreq{EQL}{DIGIT}+ {
    sscanf(yytext, " %[a-zA-Z] = %d", name, &a);
    mnl->rec_ev = a;
    if(myverbose!=0) printf("  Frequency for computing EV's set to %d in line %d monomial %d\n", mnl->rec_ev, line_of_file, current_monomial);
  }
}
<NDPOLYMONOMIAL,CLPOLYMONOMIAL>{
  {SPC}*MaxPtildeDegree{EQL}{DIGIT}+ {
    sscanf(yytext, " %[a-zA-Z] = %d", name, &a);
    mnl->MaxPtildeDegree = a;
    phmc_max_ptilde_degree = a;
    if(myverbose!=0) printf("  Maximal Degree of Ptilde set to %d line %d monomial %d\n", mnl->MaxPtildeDegree, line_of_file, current_monomial);
  }
  {SPC}*DegreeOfMDPolynomial{EQL}{DIGIT}+ {
    sscanf(yytext, " %[a-zA-Z] = %d", name, &a);
    mnl->MDPolyDegree = a;
    degree_of_p = a;
    if(myverbose!=0) printf("  Degree of MD polynomial set to %d line %d monomial %d\n", a, line_of_file, current_monomial);
  }
  {SPC}*PrecisionPtilde{EQL}{FLT} {
    sscanf(yytext, " %[a-zA-Z] = %lf",name , &c);
    mnl->PrecisionPtilde = c;
    if(myverbose!=0) printf("  Precision for Ptilde set to %e line %d monomial %d\n", c, line_of_file, current_monomial);
  }
  {SPC}*PrecisionHfinal{EQL}{FLT} {
    sscanf(yytext, " %[a-zA-Z] = %lf",name , &c);
    mnl->PrecisionHfinal = c;
    if(myverbose!=0) printf("  Precision for final H set to %e line %d monomial %d\n", c, line_of_file, current_monomial);
  }
  {SPC}*ComputeOnlyEVs{EQL}yes {
    phmc_compute_evs=1;
    if(myverbose!=0) printf("  Compute only heavy EVs set to true line %d monomial %d\n", line_of_file, current_monomial);
  }
  {SPC}*ComputeOnlyEVs{EQL}no {
    phmc_compute_evs=0;
    if(myverbose!=0) printf("  Compute only heavy EVs set to false line %d monomial %d\n", line_of_file, current_monomial);
  }
}

<POLYMONOMIAL>{
  {SPC}*Degree{EQL}{DIGIT}+ {
    sscanf(yytext, " %[a-zA-Z] = %d", name, &a);
    mnl->MDPolyDegree = a;
    if(myverbose!=0) printf("  Degree of degenerate MD polynomial set to %d line %d monomial %d\n", mnl->MDPolyDegree, line_of_file, current_monomial);
  }
  {SPC}*Lmin{EQL}{FLT} {
    sscanf(yytext, " %[a-zA-Z] = %lf", name, &c);
    mnl->MDPolyLmin = c;
    if(myverbose!=0)
      printf("  lower bound of degenerate MD polynomial set to %f line %d monomial %d\n",
        mnl->MDPolyLmin, line_of_file, current_monomial);
  }
  {SPC}*Lmax{EQL}{FLT} {
    sscanf(yytext, " %[a-zA-Z] = %lf", name, &c);
    mnl->MDPolyLmax = c;
    if(myverbose!=0)
      printf("  upper bound of degenerate MD polynomial set to %f line %d monomial %d\n",
       mnl->MDPolyLmax, line_of_file, current_monomial);
  }
}


<POLYMONOMIAL,NDPOLYMONOMIAL,CLPOLYMONOMIAL>{
  {SPC}*LocNormConst{EQL}{FLT} {
    sscanf(yytext, " %[a-zA-Z] = %lf", name, &c);
    mnl->MDPolyLocNormConst = c;
    if(myverbose!=0)
      printf("  local normalisation constant MD polynomial set to %f line %d monomial %d\n",
       mnl->MDPolyLocNormConst, line_of_file, current_monomial);
  }
  {SPC}*RootsFile{EQL} {
    cstring_to_parse=mnl->MDPolyRootsFile;
    cstring_caller = YY_START;
    BEGIN(MCSTR);
  }
}


<MNAME>{NAME} {
  if(myverbose) printf("  monomial named \"%s\" line %d monomial %d\n", yytext, line_of_file, current_monomial);
  strcpy((*mnl).name, yytext);
  BEGIN(name_caller);
}

<MCSTR>{CSTR} {
  if(myverbose) printf("  monomial named \"%s\" line %d monomial %d\n", yytext, line_of_file, current_monomial);
  strcpy(cstring_to_parse, yytext);
  rmQuotes(cstring_to_parse);
  /* reset variable */
  cstring_to_parse=NULL;
  BEGIN(cstring_caller);
}


<MSOLVER>{
  CG {
    if(myverbose) printf("  Solver set to \"%s\" line %d monomial %d\n", yytext, line_of_file, current_monomial);
    mnl->solver = CG;
    BEGIN(solver_caller);
  }
  mixedCG {
    if(myverbose) printf("  Solver set to \"%s\" line %d monomial %d\n", yytext, line_of_file, current_monomial);
    mnl->solver = MIXEDCG;
    BEGIN(solver_caller);
  }
  rgmixedCG {
    if(myverbose) printf("  Solver set to \"%s\" line %d monomial %d\n", yytext, line_of_file, current_monomial);
    mnl->solver = RGMIXEDCG;
    BEGIN(solver_caller);
  }
  BICGSTAB {
    if(myverbose) printf("  Solver set to \"%s\" line %d monomial %d\n", yytext, line_of_file, current_monomial);
    mnl->solver = BICGSTAB;
    BEGIN(solver_caller);
  }
  mixedbicgstab {
    if(myverbose) printf("  Solver set to \"%s\" line %d monomial %d\n", yytext, line_of_file, current_monomial);
    mnl->solver=MIXEDBICGSTAB;
    BEGIN(name_caller);
  }
  DDalphaAMG {
#ifdef DDalphaAMG
    if(myverbose) printf("  Solver set to DDalphaAMG line %d operator %d\n", line_of_file, current_operator);
    mnl->solver = MG;
    BEGIN(solver_caller);
#else
    printf("ERROR line %d operator %d: DDalphaAMG library not included\n", line_of_file, current_operator);
    exit(1);
#endif
  }
}

<RATMSOLVER>{
  rgmixedCG {
    if(myverbose) printf("  Solver set to \"%s\" line %d monomial %d\n", yytext, line_of_file, current_monomial);
    mnl->solver = RGMIXEDCG;
    BEGIN(solver_caller);
  }
  cgmms {
    if(myverbose) printf("  Solver set to \"%s\" line %d monomial %d\n", yytext, line_of_file, current_monomial);
    mnl->solver = CGMMS;
    BEGIN(solver_caller);
  }
  DDalphaAMG {
#ifdef DDalphaAMG
    if(myverbose) printf("  Solver set to DDalphaAMG line %d operator %d\n", line_of_file, current_operator);
    mnl->solver = MG;
    BEGIN(solver_caller);
#else
    printf("ERROR line %d operator %d: DDalphaAMG library not included\n", line_of_file, current_operator);
    exit(1);
#endif
  }
}

<NDMSOLVER>{
  cgmmsnd {
    if(myverbose) printf("  Solver set to \"%s\" line %d monomial %d\n", yytext, line_of_file, current_monomial);
    mnl->solver = CGMMSND;
    BEGIN(solver_caller);
  }
  rgmixedcg {
    if(myverbose) printf("  Solver set to \"%s\" line %d monomial %d\n", yytext, line_of_file, current_monomial);
    mnl->solver = RGMIXEDCG;
    BEGIN(solver_caller);
  }
  mixedCGmmsnd {
    if(myverbose) printf("  Solver set to \"%s\" line %d monomial %d\n", yytext, line_of_file, current_monomial);
    mnl->solver = MIXEDCGMMSND;
    BEGIN(solver_caller);
  }
  MCR {
    if(myverbose) printf("  Solver set to \"%s\" line %d monomial %d\n", yytext, line_of_file, current_monomial);
    mnl->solver = 14;
    BEGIN(solver_caller);
  }
  DDalphaAMG {
#ifdef DDalphaAMG
    if(myverbose) printf("  Solver set to DDalphaAMG line %d operator %d\n", line_of_file, current_operator);
    mnl->solver = MG;
    BEGIN(solver_caller);
#else
    printf("ERROR line %d operator %d: DDalphaAMG library not included\n", line_of_file, current_operator);
    exit(1);
#endif
  }
}

<GTYPE>{
  Wilson {
    mnl->gtype = 0;
    mnl->c1 = 0.;
    mnl->use_rectangles = 0;
    g_rgi_C1 = 0.;
    /* g_dbw2rand = 0; */
    BEGIN(GAUGEMONOMIAL);
  }
  tlsym {
    mnl->gtype = 1;
    mnl->c1 = -0.083333333;
    g_rgi_C1 = -0.083333333;
    mnl->use_rectangles = 1;
    g_dbw2rand = 1;
    BEGIN(GAUGEMONOMIAL);
  }
  Iwasaki {
    mnl->gtype = 2;
    mnl->c1 = -0.331;
    g_rgi_C1 = -0.331;
    mnl->use_rectangles = 1;
    g_dbw2rand = 1;
    BEGIN(GAUGEMONOMIAL);
  }
  user {
    mnl->gtype = 3;
    BEGIN(GAUGEMONOMIAL);
  }
  DBW2 {
    mnl->gtype = 4;
    mnl->c1 = -1.4088;
    g_rgi_C1 = -1.4088;
    g_dbw2rand = 1;
    mnl->use_rectangles = 1;
    BEGIN(GAUGEMONOMIAL);
  }
}

<INITINTEGRATOR>egrator{SPC}* {
  Integrator.no_timescales = -1;
  Integrator.tau = 1.;
  Integrator.monitor_forces = 0;
  for(i = 0; i < 10; i++) {
    Integrator.lambda[i] = _default_2mn_lambda;
    Integrator.type[i] = MN2;
  }
  if(myverbose) printf("initialising integrator line %d\n", line_of_file);
  BEGIN(INTEGRATOR);
}
<INTEGRATOR>{
  {SPC}*Type{DIGIT}{EQL}{TYPE} {
    type = (char*)malloc(100*sizeof(char));
    sscanf(yytext, " %[a-zA-Z]%d = %s", name, &a, type);
    if(strcmp(type, "LEAPFROG")==0) {
      Integrator.type[a] = LEAPFROG;
    }
    else if(strcmp(type, "2MN")==0) {
      Integrator.type[a] = MN2;
    }
    else if(strcmp(type, "2MNPOSITION")==0) {
      Integrator.type[a] = MN2p;
    }
    else if(strcmp(type, "OMF4")==0) {
      Integrator.type[a] = OMF4;
    }
    else if(strcmp(type, "2MNFG")==0) {
      Integrator.type[a] = MN2FG;
    }
    else {
      fprintf(stderr, "Unknown integrator type %s in line %d\n", yytext, line_of_file);
      exit(1);
    }

    if(myverbose) printf("  timescale %d type = %s line %d\n", a, type, line_of_file);
    free(type);
  }
  {SPC}*MonitorForces{EQL}yes {
    Integrator.monitor_forces = 1;
    if(myverbose) printf("  Force monitoring switched on in line %d\n", line_of_file);
    BEGIN(INTEGRATOR);
  }
  {SPC}*MonitorForces{EQL}no {
    Integrator.monitor_forces = 0;
    if(myverbose) printf("  Force monitoring switched off in line %d\n", line_of_file);
    BEGIN(INTEGRATOR);
  }
  {SPC}*IntegrationSteps{DIGIT}{EQL}{DIGIT}+ {
    sscanf(yytext, " %[a-zA-Z]%d = %d", name, &a, &b);
    if(myverbose) printf("  timescale %d steps=%d line %d\n", a, b, line_of_file);
    Integrator.n_int[a] = b;
    BEGIN(INTEGRATOR);
  }
  {SPC}*Lambda{DIGIT}{EQL}{FLT} {
    sscanf(yytext, " %[a-zA-Z]%d = %lf", name, &a, &c);
    Integrator.lambda[a] = c;
    if(myverbose) printf("  timescale %d Lambda=%f line %d\n", a, c, line_of_file);
    BEGIN(INTEGRATOR);
  }
  {SPC}*NumberOfTimescales{EQL}{DIGIT}+ {
    sscanf(yytext, " %[a-zA-Z] = %d", name, &a);
    if(myverbose) printf("  Number of timescales set to %d line %d\n", a, line_of_file);
    if(a > 10) {
      if(g_proc_id == 0) fprintf(stderr, "maximal number of timescales is 10! Aborting...!\n");
      exit(-1);
    }
    Integrator.no_timescales = a;
    BEGIN(INTEGRATOR);
  }
  {SPC}*Tau{EQL}{FLT} {
    sscanf(yytext, " %[a-zA-Z] = %lf", name, &c);
    if(myverbose) printf("  tau set to %e line %d\n", c, line_of_file);
    Integrator.tau = c;
    BEGIN(INTEGRATOR);
  }
  EndIntegrator{SPC}* {
    if(Integrator.no_timescales == -1) {
      fprintf(stderr, "NumberOfTimescales must be specified!\n");
      exit(1);
    }
    if(myverbose) printf("Integrators parsed line %d\n\n", line_of_file);
    BEGIN(0);
  }
}

<SOURCETYPE>{
  Point {
    SourceInfo.type = SRC_TYPE_POINT;
    if(myverbose) printf("Using Point Sources\n");
  }
  Volume {
    SourceInfo.type = SRC_TYPE_VOL;
    if(myverbose) printf("Using Volume Sources\n");
  }
  TimeSlice {
    SourceInfo.type = SRC_TYPE_TS;
    if(myverbose) printf("Using TimeSlice Sources\n");
  }
  PionTimeSlice {
    SourceInfo.type = SRC_TYPE_PION_TS;
    if(myverbose) printf("Using PionTimeSlice Sources\n");
  }
  GenPionTimeSlice {
    SourceInfo.type = SRC_TYPE_GEN_PION_TS;
    if(myverbose) printf("Using GenPionTimeSlice Sources\n");
  }
}
<PROPSPLIT>{
  yes {
    PropInfo.splitted = 1;
    if(myverbose) printf("Writing Propagators in seperate files\n");
  }
  no {
    PropInfo.splitted = 0;
    if(myverbose) printf("Writing all Propagators per gauge into one single files\n");
  }
}

<NOSAMPLES>{DIGIT}+ {
  no_samples = atoi(yytext);
  if(myverbose) printf("Inverting for %d samples\n", no_samples);
}

<INITMEASUREMENT>{TYPE} {
 if(myverbose) printf("initialising measurements line %d\n", line_of_file);
  current_measurement++;
  meas = &measurement_list[current_measurement];
  meas->id = current_measurement;
  meas->direction = 0;
  meas->max_iter = 15000;
  if(strcmp(yytext, "CORRELATORS")==0) {
    meas->type = ONLINE;
    strcpy((*meas).name, "CORRELATORS");
  }
  else if(strcmp(yytext, "PIONNORM")==0) {
    meas->type = PIONNORM;
    strcpy((*meas).name, "PIONNORM");
  }
  else if(strcmp(yytext, "POLYAKOVLOOP")==0) {
    meas->type = POLYAKOV;
    strcpy((*meas).name, "POLYAKOV");
  } 
  else if(strcmp(yytext, "ORIENTEDPLAQUETTES")==0) {
    meas->type = ORIENTED_PLAQUETTES;
    strcpy(meas->name, "ORIENTEDPLAQUETTES");
  }
  else if(strcmp(yytext, "GRADIENTFLOW")==0) {
    meas->type = GRADIENT_FLOW;
    strcpy(meas->name, "GRADIENTFLOW");
  }
  else {
    fprintf(stderr, "Unknown measurement type %s in line %d\n", yytext, line_of_file);
    exit(1);
  }
  // set default measurement parameters. This has to be done here, since init_measurements would
  // override the inputs
  meas->freq = _default_measurement_freq;

  meas->no_samples = 1;
  meas->all_time_slices = 0;

  meas->gf_eps = _default_gf_eps;
  meas->gf_tmax = _default_gf_tmax;

  if(!reread) {
    if(add_measurement(meas->type) < 0) {
      fprintf(stderr, "Something went wrong in adding measurements\nAborting...!\n");
      exit(1);
    }
  }
  if(myverbose) printf("initialising measurement with type %s %d line %d\n", yytext, meas->type, line_of_file);
  if(myverbose) printf("measurement has id %d\n", current_measurement);

  if(meas->type == ONLINE) BEGIN(ONLINEMEAS);
  else if(meas->type == PIONNORM) BEGIN(PIONNORMMEAS);
  else if(meas->type == POLYAKOV) BEGIN(PLOOP);
  else if(meas->type == ORIENTED_PLAQUETTES) BEGIN(ORIENTEDPLAQUETTESMEAS);
  else if(meas->type == GRADIENT_FLOW) BEGIN(GRADIENTFLOWMEAS);
}

<ONLINEMEAS,PIONNORMMEAS,PLOOP,ORIENTEDPLAQUETTESMEAS,GRADIENTFLOWMEAS>{
  ^EndMeasurement{SPC}* {
    if(myverbose) printf("Measurement with id %d parsed in line %d\n\n", meas->id, line_of_file);
    BEGIN(0);
  }
  {SPC}*Frequency{EQL}{DIGIT}+ {
    sscanf(yytext, " %[a-zA-Z] = %d", name, &a);
    meas->freq = a;
    if(myverbose!=0) printf("  Frequency for measurement with id %d set to %d\n", meas->id, meas->freq);
  }
}

<ONLINEMEAS,PIONNORMMEAS>{
  {SPC}*MaxSolverIterations{EQL}{DIGIT}+ {
    sscanf(yytext, " %[a-zA-Z] = %d", name, &a);
    meas->max_iter = a;
    if(myverbose) printf("  MaxSolverIterations set to %d line %d measurement id=%d\n", a, line_of_file, meas->id);
  }
  {SPC}*AllTimeSlices{EQL}yes {
    meas->all_time_slices = 1;
    if(myverbose) printf("  Inversions for correlators will be done on all time-slices, line %d, measurement id=%d\n", line_of_file, meas->id);
  }
  {SPC}*AllTimeSlices{EQL}no {
    meas->all_time_slices = 0;
    if(myverbose) printf("  Inversions for correlators will NOT be done on all time-slices, line %d, measurement id=%d\n", line_of_file, meas->id);
  }
  {SPC}*NoSamples{EQL}{DIGIT}+ {
    sscanf(yytext, " %[a-zA-Z] = %d", name, &a);
    meas->no_samples = a;
    if(myverbose) printf("  Number of samples for correlators measurement set to %d, line %d measurement id=%d\n", a, line_of_file, meas->id);
  }
}

<GRADIENTFLOWMEAS>{
  {SPC}*StepSize{EQL}{FLT} {
    sscanf(yytext, " %[2a-zA-Z] = %lf", name, &c);
    meas->gf_eps = c;
    if(myverbose) printf("  Gradient flow step size set to %lf line %d, measurement id=%d\n", meas->gf_eps, line_of_file, meas->id);
  }
  {SPC}*MaxFlowTime{EQL}{FLT} {
    sscanf(yytext, " %[2a-zA-Z] = %lf", name, &c);
    meas->gf_tmax = c;
    if(myverbose) printf("  Maximum gradient flow time size set to %lf line %d, measurement id=%d\n", meas->gf_tmax, line_of_file, meas->id);
  }
} 

<PLOOP>{
  {SPC}*Direction{EQL}[03] {
    sscanf(yytext, " %[a-zA-Z] = %d", name, &a);
    meas->direction = a;
    if(myverbose!=0) fprintf(stderr, "  Direction for polyakov loop set to %d\n", meas->direction);
  }
}

<TT>{DIGIT}+                  {
#ifndef FIXEDVOLUME
  T_global = atoi(yytext);
  if(myverbose!=0) printf("T =%s\n", yytext);
#endif
}
<LL>{DIGIT}+                  {
#ifndef FIXEDVOLUME
  L = atoi(yytext);
  if(myverbose!=0) printf("L =%s\n", yytext);
#endif
}
<LLX>{DIGIT}+                  {
#ifndef FIXEDVOLUME
  LX = atoi(yytext);
  if(myverbose!=0) printf("LX =%s\n", yytext);
#endif
}
<LLY>{DIGIT}+                  {
#ifndef FIXEDVOLUME
  LY = atoi(yytext);
  if(myverbose!=0) printf("LY =%s\n", yytext);
#endif
}
<LLZ>{DIGIT}+                  {
#ifndef FIXEDVOLUME
  LZ = atoi(yytext);
  if(myverbose!=0) printf("LZ =%s\n", yytext);
#endif
}
<NPROCX>{DIGIT}+              {
#ifndef FIXEDVOLUME
  N_PROC_X = atoi(yytext);
  if(myverbose!=0) printf("Nr of processors in x direction = %s\n", yytext);
#endif
}
<NPROCY>{DIGIT}+              {
#ifndef FIXEDVOLUME
  N_PROC_Y = atoi(yytext);
  if(myverbose!=0) printf("Nr of processors in y direction = %s\n", yytext);
#endif
}
<NPROCZ>{DIGIT}+              {
#ifndef FIXEDVOLUME
  N_PROC_Z = atoi(yytext);
  if(myverbose!=0) printf("Nr of processors in z direction = %s\n", yytext);
#endif
}
<PCOMP>{DIGIT}+               {
  propagator_comparison=atoi(yytext);
  if(myverbose!=0) printf("propagator_comparison = %s \n", yytext);
}
<NBCORES>{DIGIT}+               {
  nb_cores=atoi(yytext);
  if(myverbose!=0) printf("nb_cores = %s \n", yytext);
}
<OMPNUMTHREADS>{DIGIT}+ {
  omp_num_threads=atoi(yytext);
  if(myverbose!=0) printf("omp_num_threads = %s \n", yytext);
}
<DFLNBLOCKT>{DIGIT}+               {
  nblocks_t=atoi(yytext);
  if(myverbose!=0) printf("nblocks_t = %s \n", yytext);
}
<DFLNBLOCKX>{DIGIT}+               {
  nblocks_x=atoi(yytext);
  if(myverbose!=0) printf("nblocks_x = %s \n", yytext);
}
<DFLNBLOCKY>{DIGIT}+               {
  nblocks_y=atoi(yytext);
  if(myverbose!=0) printf("nblocks_y = %s \n", yytext);
}
<DFLNBLOCKZ>{DIGIT}+               {
  nblocks_z=atoi(yytext);
  if(myverbose!=0) printf("nblocks_z = %s \n", yytext);
}
<DFLFIELDITER>{DIGIT}+               {
  dfl_field_iter=atoi(yytext);
  if(myverbose!=0) printf("dfl_fields_iter = %s \n", yytext);
}
<DFLPOLYITER>{DIGIT}+               {
  dfl_poly_iter=atoi(yytext);
  if(myverbose!=0) printf("dfl_poly_iter = %s \n", yytext);
}
<SEED>{DIGIT}+               {
  random_seed=atoi(yytext);
  if(myverbose!=0) printf("seed=%s \n", yytext);
}
<RLXDLEVEL>[12] {
  rlxd_level = atoi(yytext);
  if(myverbose!=0) printf("RanluxdLevel set to %d \n", rlxd_level);
}
<KAPPA>{FLT}  {
  g_kappa=atof(yytext);
  if(myverbose!=0) printf("kappa=%s \n", yytext);
}
<MUBAR>{FLT}  {
  g_mubar=atof(yytext);
  if(myverbose!=0) printf("2 kappa mubar=%s \n", yytext);
}
<EPSBAR>{FLT}  {
  g_epsbar=atof(yytext);
  if(myverbose!=0) printf("2 kappa epsbar=%s \n", yytext);
}
<MU>{FLT}  {
  g_mu1=atof(yytext);
  if(myverbose!=0) printf("2 kappa mu=%s \n", yytext);
}
<CSW>{FLT}  {
  g_c_sw=atof(yytext);
  if(myverbose!=0) printf("c_sw=%s \n", yytext);
}
<STARTCOND>{
  cold {
    startoption=0;
    if(myverbose!=0) printf("Start Condition is %s \n",yytext);
  }
  hot {
    startoption=1;
    if(myverbose!=0) printf("Start Condition is %s \n",yytext);
  }
  restart {
    startoption=2;
    if(myverbose!=0) printf("Start Condition is %s \n",yytext);
  }
  continue {
    startoption=3;
    if(myverbose!=0) printf("Start Condition is %s \n",yytext);
  }
}
<THERMSWEEPS>{DIGIT}+ {
  Ntherm=atoi(yytext);
  if(myverbose!=0) printf("Nterm= %s \n",yytext);
}
<NMEAS>{DIGIT}+ {
  Nmeas=atoi(yytext); 
  if(myverbose!=0) printf("Nmeas= %s \n",yytext);
}
<NSAVE>{DIGIT}+ {
  Nsave=atoi(yytext);
  if(myverbose!=0) printf("Nsave= %s \n",yytext);
}
<GMRESM>{DIGIT}+ {
  gmres_m_parameter = atoi(yytext);
  if(myverbose!=0) printf("Use Krylov Space of size %d in GMRES \n", gmres_m_parameter);
}
<GMRESDRNEV>{DIGIT}+ {
  gmresdr_nr_ev = atoi(yytext);
  if(myverbose!=0) printf("Deflate %d eigenvectors in GMRES-DR \n", gmresdr_nr_ev);
}
<PRECON>{
  none {
    if(myverbose!=0) printf("Using no right preconditioner \n");
  }
  polynomial {
    if(myverbose!=0) printf("Using polynomial as right preconditioner \n");
  }
  cg {
    if(myverbose!=0) printf("Using cg as right preconditioner \n");
  }
}
<WRITECP>yes     {
  write_cp_flag=1;
  if(myverbose!=0) printf("Write Checkpoints\n");
}
<WRITECP>no     {
  write_cp_flag=0;
  if(myverbose!=0) printf("Don't write Checkpoints\n");
}
<SUBPROCESS>yes     {
  subprocess_flag=1;
  if(myverbose!=0) printf("Running in \'subprocess\' mode\n");
}
<SUBPROCESS>no     {
  subprocess_flag=0;
  if(myverbose!=0) printf("Not running in \'subprocess\' mode\n");
}
<LOWMEM>yes     {
  lowmem_flag=1;
  if(g_proc_id==0) printf("!!! WARNING: running in \'lowmem\' mode, most functionality will NOT work without explicit memory allocations!\n");
}
<LOWMEM>no     {
  lowmem_flag=0;
  if(myverbose!=0) printf("Not running in \'lowmem\' mode\n");
}
<DSBLIOCHECK>yes {
  g_disable_IO_checks = 1;
  if(myverbose!=0) printf("Disable IO checks (and readback in case of Lemon IO)\n");
}
<DSBLIOCHECK>no {
  g_disable_IO_checks = 0;
  if(myverbose!=0) printf("Enable IO checks (and readback in case of Lemon IO)\n");
}
<DSBLSRCIOCHECK>yes {
  g_disable_src_IO_checks = 1;
  if(myverbose!=0) printf("Disable IO checks for sources\n");
}
<DSBLSRCIOCHECK>no {
  g_disable_src_IO_checks = 0;
  if(myverbose!=0) printf("Enable IO checks for sources\n");
}
<CPINT>{DIGIT}+   {
  cp_interval=atoi(yytext);
  if(myverbose!=0) printf("Write Checkpoint all %s measurements\n",yytext);
}
<GAUGEINPUTFILE>{FILENAME} {
  int length = strlen(yytext)+1;
  if(length >= sizeof(gauge_input_filename)/sizeof(char) )
    yy_fatal_error("Filename GaugeConfigInputFile too long! (see read_input.l)\n");

  strcpy(gauge_input_filename,yytext);
  if(myverbose!=0) printf("Gauge Configuration input filename set to %s\n",yytext);
}
<NSTORE>{DIGIT}+   {
  nstore=atoi(yytext);
  if(myverbose!=0) printf("Initial store counter set to %s\n",yytext);
}
<NSTORE>readin {
  nstore=-1;
  if(myverbose!=0) printf("Trying to read InitialStoreCounter from file .nstore_counter\n");
}
<IOPROC>all         {
  g_stdio_proc = -1;
  if(myverbose!=0) printf("All processors will give output to stdout\n");
}
<IOPROC>no          {
  g_stdio_proc = -2;
  if(myverbose!=0) printf("No processor will give output to stdout\n");
}
<IOPROC>{DIGIT}+    {
  g_stdio_proc = atoi(yytext);
  if(myverbose!=0) printf("processor %s will give output to stdout\n", yytext);
}
<IDX>{DIGIT}+ {
  index_start = atoi(yytext);
  index_end = index_start+1;
  if((index_start < 0)||(index_start > 11)){
    printf("Error in line %d! index_start must be in [0,11]! Exiting...!\n", line_of_file);
    exit(1);
  }
  if(myverbose!=0) printf("inverting for index %s\n", yytext);
}
<IDX>{IDXEX}  {
  sscanf(yytext, "-%d", &index_end);
  if((index_end < 0)||(index_end > 11)){
    printf("Error in line %d! index_end must be in [0,11]! Exiting...!\n", line_of_file);
    exit(1);
  }
  if(myverbose!=0) printf("inverting up to spin-color index %d\n", index_end);
  index_end+=1;
}
<BOUNDT>{FLT} {
  X0 = atof(yytext);
  if(myverbose != 0) printf("X0 for boundary cond. in t-direction set to %e times pi\n", X0);
}
<BOUNDX>{FLT} {
  X1 = atof(yytext);
  if(myverbose != 0) printf("X1 for boundary cond. in x-direction set to %e times pi\n", X1);
}
<BOUNDY>{FLT} {
  X2 = atof(yytext);
  if(myverbose != 0) printf("X2 for boundary cond. in y-direction set to %e times pi\n", X2);
}
<BOUNDZ>{FLT} {
  X3 = atof(yytext);
  if(myverbose != 0) printf("X3 for boundary cond. in z-direction set to %e times pi\n", X3);
}
<READSOURCE>yes     {
  read_source_flag=1;
  if(myverbose!=0) printf("Read inversion source from file\n");
}
<READSOURCE>no     {
  read_source_flag=0;
  if(myverbose!=0) printf("Don't read inversion source from file\n");
}
<READSOURCE>nobutsave {
  read_source_flag=2;
  if(myverbose!=0) printf("Don't read inversion source from file, but save the one generated\n");
}
<SOURCEFILE>{FILENAME} {
  if(SourceInfo.basename != NULL) free(SourceInfo.basename);
  SourceInfo.basename = (char*)malloc((strlen(yytext)+1)*sizeof(char));
  strcpy(SourceInfo.basename, yytext);
  if(myverbose!=0) printf("source input filename set to %s of length %lu\n",yytext,strlen(yytext));
}
<PROPFILE>{FILENAME} {
  if(PropInfo.basename != NULL) free(PropInfo.basename);
  PropInfo.basename = (char*)malloc((strlen(yytext)+1)*sizeof(char));
  strcpy(PropInfo.basename, yytext);
  if(myverbose!=0) printf("propagator output filename set to %s\n",yytext);
}
<SOURCEFORMAT>etmc      {
  SourceInfo.format = 0;
  if(myverbose!=0) printf("Using standard ETMC binary format for source input file\n");
}
<SOURCEFORMAT>cmi      {
  SourceInfo.format = 11;
  if(myverbose!=0) printf("Using CM format for source input file\n");
}
<SOURCEFORMAT>gwc      {
  SourceInfo.format = 10;
  if(myverbose!=0) printf("Using GWC format for source input file\n");
}
<SOURCETS>{DIGIT}+   {
  SourceInfo.t = atoi(yytext);
  SourceInfo.automaticTS = 0;
  if(myverbose!=0) printf("Using only timeslice %s of the source, padding the rest with zeros\n", yytext);
}
<SOURCETS>detect {
  SourceInfo.automaticTS = 1;
  if(myverbose!=0) printf("Try to detect source timeslice automatically\n");
}
<RELPREC>yes  {
  g_relative_precision_flag = 1;
  if(myverbose!=0) printf("Using relative precision\n");
}
<RELPREC>no  {
  g_relative_precision_flag = 0;
  if(myverbose!=0) printf("Using absolute precision\n");
}
<REVCHECK>yes {
  return_check_flag = 1;
  if(myverbose!=0) printf("Perform checks of Reversibility\n");
}
<REVCHECK>no {
  return_check_flag = 0;
  if(myverbose!=0) printf("Don't perform checks of Reversibility\n");
}
<REVINT>{DIGIT}+ {
  return_check_interval = atoi(yytext);
  if(myverbose!=0) printf("Check reversibility all %d trajectories\n", return_check_interval);
}
<DEBUG>{DIGIT}+ {
  g_debug_level = atoi(yytext);
  if(myverbose!=0) printf("Debug level = %d\n", g_debug_level);
}
<GAUGERPREC>32 {
  gauge_precision_read_flag = 32;
  if(myverbose!=0) printf("Read gauges in 32 Bit precision!\n");
}
<GAUGERPREC>64 {
  gauge_precision_read_flag = 64;
  if(myverbose!=0) printf("Read gauges in 64 Bit precision!\n");
}
<GAUGEWPREC>32 {
  gauge_precision_write_flag = 32;
  if(myverbose!=0) printf("Save gauges in 32 Bit precision!\n");
}
<GAUGEWPREC>64 {
  gauge_precision_write_flag = 64;
  if(myverbose!=0) printf("Save gauges in 64 Bit precision!\n");
}
<REPRORND>yes {
  reproduce_randomnumber_flag = 1;
  if(myverbose!=0) printf("Use reproducable randomnumbers!\n");
}
<REPRORND>no {
  reproduce_randomnumber_flag = 0;
  if(myverbose!=0) printf("Use a different seed for each process in ranlxd!\n");
}
<SLOPPYPREC>yes {
  g_sloppy_precision_flag = 1;
  if(myverbose!=0) printf("Use sloppy precision if available!\n");
}
<SLOPPYPREC>no {
  g_sloppy_precision_flag = 0;
  if(myverbose!=0) printf("Don't use sloppy precision!\n");
}
<USESTOUT>yes {
  use_stout_flag = 1;
  if(myverbose!=0) printf("Use stout smearing for invert!\n");
}
<USESTOUT>no {
  use_stout_flag = 0;
  if(myverbose!=0) printf("Don't use stout smearing for invert!\n");
}
<PRECONDITIONING>yes {
  use_preconditioning = 1;
  if(myverbose) printf("Using project \"QCD-Preconditioning\"\n");
}
<PRECONDITIONING>no {
  use_preconditioning = 0;
  if(myverbose) printf("not using project \"QCD-Preconditioning\"\n");
}
<STOUTRHO>{FLT} {
  stout_rho=atof(yytext);
  if(myverbose!=0) printf("use stout rho=%e!\n", stout_rho);
}
<STOUTITER>{DIGIT}+ {
  stout_no_iter=atoi(yytext);
  if(myverbose!=0) printf("make %d stout iterations!\n", stout_no_iter);
}
<COMPUTEEVS>yes {
  compute_evs=1;
  if(myverbose!=0) printf("Compute Eigenvalues in invert.");
}
<COMPUTEEVS>no {
  compute_evs=0;
  if(myverbose!=0) printf("Do not compute Eigenvalues in invert.");
}
<COMPUTEEVS>readin {
  compute_evs=2;
  if(myverbose!=0) printf("Try to only read in eigenvalues and vectors in invert.");
}
<COMPUTEMN>yes {
  compute_modenumber=1;
  if(myverbose!=0) printf("Compute Mode Number using Spectral Projectors in invert.\n");
}
<COMPUTEMN>no {
  compute_modenumber=0;
  if(myverbose!=0) printf("Do not compute Mode Number using Spectral Projectors in invert.\n");
}
<COMPUTETS>yes {
  compute_topsus=1;
  if(myverbose!=0) printf("Compute Topological Susceptibility using Spectral Projectors in invert.\n");
}
<COMPUTETS>no {
  compute_topsus=0;
  if(myverbose!=0) printf("Do not compute Topological Susceptibility using Spectral Projectors in invert.\n");
}
<MSTARSQ>{FLT} { 
  mstarsq = atof(yytext);
  if(myverbose!=0) printf("Mstar^2 = %f \n", mstarsq);
}
<NOSOURCESZ2>{DIGIT}+ {
  no_sources_z2 = atoi(yytext);
  if(myverbose!=0) printf("no of Z2 random sources used for the spectral projectors method = %d\n", no_sources_z2);
}
<SRCLOC>{DIGIT}+ {
  source_location=atoi(yytext);
  if(myverbose!=0) printf("source_location = %s\n",yytext);
}
<PRECEV>{FLT} {
  eigenvalue_precision = atof(yytext);
  if(myverbose!=0) printf("precision for eigenvalues = %e\n", eigenvalue_precision);
}
<NOEV>{DIGIT}+ {
  no_eigenvalues = atoi(yytext);
  if(myverbose!=0) printf("no of eigenvalues = %d\n", no_eigenvalues);
}
<SUBEVCG>yes {
  sub_evs_cg_flag = 1;
  if(myverbose!=0) printf("project out eigenvector subspace\n");
}
<SUBEVCG>no {
  sub_evs_cg_flag = 0;
  if(myverbose!=0) printf("Do no project out eigenvector subspace\n");
}
<EO>yes {
  even_odd_flag = 1;
  if(myverbose) printf("Use even/odd preconditioning\n");
}
<EO>no {
  even_odd_flag = 0;
  if(myverbose) printf("Do not use even/odd preconditioning\n");
}
<BC>yes {
  bc_flag = 1;
  if(verbose) printf("Schroedinger Functional bc\n");
}
<BC>no {
  bc_flag = 0;
  if(verbose) printf("Periodic bc\n");
}
<WRPROPFLAG>gwc      {
  PropInfo.format = 10;
  if(myverbose!=0) fprintf(stderr, "GWC format no longer supported for writing propagators\n");
}
<WRPROPFLAG>cmi      {
  PropInfo.format = 11;
  if(myverbose!=0) fprintf(stderr, "CM format no longer supported for writing propagators\n");
}
<WRPROPFLAG>DiracFermion_Sink {
  PropInfo.format = 0;
  if(myverbose!=0) printf("Propagator type: DiracFermion_Sinks\n");
}
<WRPROPFLAG>DiracFermion_Source_Sink_Pairs {
  PropInfo.format = 1;
  if(myverbose!=0) printf("Propagator type: DiracFermion_Source_Sink_Pairs\n");
}
<WRPROPFLAG>DiracFermion_ScalarSource_TwelveSink {
  PropInfo.format = 1;
  fprintf(stderr, "Propagator type: DiracFermion_ScalarSource_TwelveSink, not yet supported\n");
}
<WRPROPFLAG>DiracFermion_ScalarSource_FourSink {
  PropInfo.format = 1;
  fprintf(stderr, "Propagator type: DiracFermion_ScalarSource_FourSink, not yet supported\n");
}
<REWEIGH>yes {
  reweighting_flag = 1;
  if(myverbose!=0) printf("Compute reweighting factor\n");
}
<REWEIGH>no {
  reweighting_flag = 0;
  if(myverbose!=0) printf("Do not compute reweighting factor\n");
}
<REWSAMPLES>{DIGIT}+ {
  reweighting_samples = atoi(yytext);
  if(myverbose!=0) printf("Number of reweighting samples set to %d\n", reweighting_samples);
}

<MIXCGIT>{DIGIT}+ {
  mixcg_maxinnersolverit = atoi(yytext);
  if(myverbose) printf("MixedCG: setting maximal inner solver iterations to %d\n", mixcg_maxinnersolverit);
}

<MIXCGEPS>{FLT}  {
  mixcg_innereps=atof(yytext);
  if(myverbose!=0) printf("MixedCG: setting inner solver eps to %s \n", yytext);
}

<*>^#   {
   comment_caller = YY_START;
   BEGIN(COMMENT);
}
<*>{SPC}*#    {
   comment_caller = YY_START;
   BEGIN(COMMENT);
}
<COMMENT>[^\n]*             {
  BEGIN(comment_caller);
}

<INITMONOMIAL,DETMONOMIAL,CLDETMONOMIAL,CLDETRATMONOMIAL,CLDETRATRWMONOMIAL,NDPOLYMONOMIAL,NDRATMONOMIAL,NDRATCORMONOMIAL,NDCLRATMONOMIAL,NDCLRATCORMONOMIAL,CLPOLYMONOMIAL,GAUGEMONOMIAL,INTEGRATOR,INITINTEGRATOR,INITMEASUREMENT,PIONNORMMEAS,ONLINEMEAS,ORIENTEDPLAQUETTESMEAS,GRADIENTFLOWMEAS,INITOPERATOR,TMOP,DBTMOP,OVERLAPOP,WILSONOP,CLOVEROP,DBCLOVEROP,POLYMONOMIAL,PLOOP,INITGPU,GPU,RATMONOMIAL,RATCORMONOMIAL,CLRATMONOMIAL,CLRATCORMONOMIAL,INITDEFLATION,DEFLATION,INITMULTIGRID,MULTIGRID,INITEXTERNALINVERTER,QUDAINVERTER,QPHIXINVERTER,NDDETRATMONOMIAL,NDCLDETRATMONOMIAL>{SPC}*\n   {
  line_of_file++;
}
<*>{SPC}*\n                       {
  line_of_file++;
  BEGIN(0);
}

<*>. {
  BEGIN(ERROR);
}
<ERROR>[^\t\n]*             {
  if(g_proc_id == 0) fprintf(stderr, "Parsing error in line %d\nAborting...!\n", line_of_file);
  if(g_proc_id == 0) fprintf(stderr, "Could not make sense out off: %s\n", yytext);
  exit(1);
}


%%

/*
 *  Dummy (but not dumb) routine - well, function
 */

int yywrap()
{
  return(1);
}

/*
 * This is the function to parse the input file.
 * default values for all paramters will be set
 * correspondig to settings in
 * default_input_values.h
 *
 * read_input expects the filename of the input file
 * as an input parameter.
 *
 * read_input returns 2 if the input file did not exist 
 */

int read_input(char * conf_file){

  /********************************************
   * Setting default values!
   ********************************************/
  reread = 0;
#ifndef FIXEDVOLUME
  T_global = _default_T_global;
  L = _default_L;
  LX = _default_LX;
  LY = _default_LY;
  LZ = _default_LZ;
  N_PROC_X = _default_N_PROC_X;
  N_PROC_Y = _default_N_PROC_Y;
  N_PROC_Z = _default_N_PROC_Z;
#endif
  propagator_comparison = 0;
  nb_cores = 1;

  omp_num_threads=_default_omp_num_threads;

  nblocks_t = 1;
  nblocks_x = 1;
  nblocks_y = 1;
  nblocks_z = 1;

  dfl_field_iter = 80;
  dfl_poly_iter = 20;

  g_kappa = _default_g_kappa;
  g_mubar = _default_g_mubar;
  g_epsbar = _default_g_epsbar;
  g_mu = _default_g_mu;
  g_c_sw = _default_c_sw;
  g_mu1 = _default_g_mu1;
  g_mu2 = _default_g_mu2;
  g_mu3 = _default_g_mu3;
  g_shift = _default_g_shift;
  g_dbw2rand = 0;
  g_running_phmc = 0;
  g_beta = _default_g_beta;
  g_N_s = _default_g_N_s;
  g_dflgcr_flag = _default_g_dflgcr_flag;

  Msap_precon = _default_Msap_precon;
  NiterMsap = _default_NiterMsap;
  NcycleMsap = _default_NcycleMsap;
  kappa_Msap = _default_kappa_Msap;
  mu_Msap = _default_mu_Msap;

  NsmoothMsap_dflgen = _default_NsmoothMsap_dflgen;
  NiterMsap_dflgen = _default_NiterMsap_dflgen;
  NcycleMsap_dflgen = _default_NcycleMsap_dflgen;
  usePL = _default_usePL;
  mu_dflgen = _default_mu_dflgen;
  kappa_dflgen = _default_kappa_dflgen;
  mu_dfl = _default_mu_dfl;
  kappa_dfl = _default_kappa_dfl;
  little_solver_low_prec = _default_little_solver_low_prec;
  little_solver_high_prec = _default_little_solver_high_prec;
  little_solver_max_iter = _default_little_solver_max_iter;
  little_evenodd = _default_little_evenodd;
  little_solver = _default_little_solver;
  little_gmres_m_parameter = _default_little_gmres_m_parameter;

  random_seed = _default_random_seed;
  rlxd_level = _default_rlxd_level;
  startoption = _default_startoption;
  Ntherm = _default_Ntherm;
  Nmeas = _default_Nmeas;
  Nsave = _default_Nsave;
  write_cp_flag = _default_write_cp_flag;
  cp_interval = _default_cp_interval;
  nstore = _default_nstore;
  strcpy(rlxd_input_filename, _default_rlxd_input_filename);
  strcpy(gauge_input_filename, _default_gauge_input_filename);
  g_stdio_proc = _default_g_stdio_proc;
  index_start = _default_index_start;
  index_end = _default_index_end;
  X0 = _default_X0;
  X1 = _default_X1;
  X2 = _default_X2;
  X3 = _default_X3;
  g_rgi_C1 = _default_g_rgi_C1;
  read_source_flag= _default_read_source_flag;
  SourceInfo.basename = (char*)malloc(strlen(_default_source_filename)+10);
  strcpy(SourceInfo.basename, _default_source_filename);
  PropInfo.basename = (char*)malloc(strlen(_default_source_filename)+10);
  strcpy(PropInfo.basename, _default_source_filename);
  PropInfo.splitted = _default_propagator_splitted;
  SourceInfo.splitted = _default_source_splitted;
  g_relative_precision_flag = _default_g_relative_precision_flag;
  return_check_flag = _default_return_check_flag;
  return_check_interval = _default_return_check_interval;
  g_debug_level = _default_g_debug_level;
  SourceInfo.t = _default_source_time_slice;
  SourceInfo.automaticTS = _default_automaticTS;
  gmres_m_parameter = _default_gmres_m_parameter;
  gmresdr_nr_ev = _default_gmresdr_nr_ev;
  gauge_precision_read_flag = _default_gauge_precision_read_flag;
  gauge_precision_write_flag = _default_gauge_precision_write_flag;
  g_disable_IO_checks = _default_g_disable_IO_checks;
  g_disable_src_IO_checks = _default_g_disable_IO_checks;
  reproduce_randomnumber_flag = _default_reproduce_randomnumber_flag;
  g_sloppy_precision_flag = _default_g_sloppy_precision_flag;
  use_stout_flag = _default_use_stout_flag;
  use_preconditioning = _default_use_preconditioning;
  stout_rho = _default_stout_rho;
  stout_no_iter = _default_stout_no_iter;

  /* check for reread ! */ 
  phmc_compute_evs = _default_phmc_compute_evs;
  compute_evs = _default_compute_evs;
  stilde_min = _default_stilde_min;
  stilde_max = _default_stilde_max;
  degree_of_p = _default_degree_of_p;
  source_location = _default_source_location;
  eigenvalue_precision = _default_eigenvalue_precision;
  no_eigenvalues = _default_no_eigenvalues;
  sub_evs_cg_flag = _default_sub_evs_cg_flag;
  phmc_exact_poly = _default_phmc_exact_poly;
  even_odd_flag = _default_even_odd_flag;
  bc_flag = _default_bc_flag;
  SourceInfo.type = _default_source_type_flag;
  no_samples = _default_no_samples;
  compute_modenumber = _default_compute_modenumber;
  compute_topsus = _default_compute_topsus;
  mstarsq = _default_mstarsq;
  no_sources_z2 = _default_no_sources_z2;
  device_num = _default_device_num;
  min_innersolver_it = _default_min_innersolver_it;
  max_mms_shifts = _default_max_mms_shifts;
  use_mixed_mms = 0;
  innersolver_precision_rel = 1.e-4;
  innersolver_precision_abs = 1.e-4;

  mixcg_innereps = _default_mixcg_innereps;
  mixcg_maxinnersolverit = _default_mixcg_maxinnersolverit;

  lowmem_flag = _default_lowmem_flag;
  subprocess_flag = _default_subprocess_flag;
 
  /* default parameters for QUDA */
  quda_input.fermionbc = TM_QUDA_THETABC;
  quda_input.pipeline = 0;
  quda_input.mg_n_level = _default_quda_mg_n_level;
  for( int level = 0; level < QUDA_MAX_MG_LEVEL; ++level){
    quda_input.mg_verbosity[level] = QUDA_SILENT;
    quda_input.mg_omega[level] = _default_quda_mg_omega;
    quda_input.mg_setup_maxiter[level] = _default_quda_mg_setup_maxiter;
    quda_input.mg_setup_tol[level] = _default_quda_mg_setup_tol;
    quda_input.mg_coarse_solver_tol[level] = _default_quda_mg_coarse_solver_tol;
    quda_input.mg_coarse_solver_maxiter[level] = _default_quda_mg_coarse_solver_maxiter;
    quda_input.mg_nu_pre[level] = _default_quda_mg_nu_pre;
    quda_input.mg_nu_post[level] = _default_quda_mg_nu_post;
    quda_input.mg_smoother_tol[level] = _default_quda_mg_smoother_tol;

    quda_input.mg_n_vec[level] = _default_quda_mg_n_vec;
    quda_input.mg_mu_factor[level] = 1.0;
    quda_input.mg_coarse_solver_type[level] = QUDA_GCR_INVERTER;
    quda_input.mg_smoother_type[level] = QUDA_MR_INVERTER ;

    quda_input.mg_use_eig_solver[level] = QUDA_BOOLEAN_NO;
    quda_input.mg_eig_tol[level] = 1.0e-3;
    quda_input.mg_eig_require_convergence[level] = QUDA_BOOLEAN_YES;
    quda_input.mg_eig_type[level] = QUDA_EIG_TR_LANCZOS;
    quda_input.mg_eig_spectrum[level] = QUDA_SPECTRUM_SR_EIG;
    quda_input.mg_eig_check_interval[level] = 5;
    quda_input.mg_eig_max_restarts[level] = 100;
    quda_input.mg_eig_use_normop[level] = QUDA_BOOLEAN_NO;
    quda_input.mg_eig_use_dagger[level] = QUDA_BOOLEAN_NO;
    quda_input.mg_eig_use_poly_acc[level] = QUDA_BOOLEAN_YES;
    quda_input.mg_eig_poly_deg[level] = 100;
    quda_input.mg_eig_amin[level] = 1.0;
    quda_input.mg_eig_amax[level] = 5.0;

    quda_input.mg_setup_ca_basis[level] = QUDA_POWER_BASIS;
    quda_input.mg_setup_ca_basis_size[level] = 4;
    quda_input.mg_setup_ca_lambda_min[level] = 0.0;
    quda_input.mg_setup_ca_lambda_max[level] = -1.0;
  
    quda_input.mg_coarse_solver_ca_basis[level] = QUDA_POWER_BASIS;
    quda_input.mg_coarse_solver_ca_basis_size[level] = 4;
    quda_input.mg_coarse_solver_ca_lambda_min[level] = 0.0;
    quda_input.mg_coarse_solver_ca_lambda_max[level] = -1.0;

    /* note: when the user does not specify any blocking parameters,
     * a reasonable set will be computed automatically in the MG setup
     * of the QUDA interface */
    for( int dim = 0; dim < 4; ++dim ){
      quda_input.mg_blocksize[level][dim] = 0;
    }
  }
  quda_input.mg_setup_inv_type = QUDA_CG_INVERTER;
  quda_input.mg_run_verify = QUDA_BOOLEAN_YES;
  quda_input.mg_enable_size_three_blocks = _default_quda_mg_enable_size_three_blocks;
  quda_input.mg_reset_setup_threshold = _default_quda_mg_reset_setup_threshold;
  quda_input.mg_run_low_mode_check = QUDA_BOOLEAN_NO;
  quda_input.mg_run_oblique_proj_check = QUDA_BOOLEAN_NO; 

  /* Put -1 in PropInfo.format to see if parse_config() will
  change the value. If not then set it to source_format_flag */
  PropInfo.format = -1;
  /********************************************/

  if(verbose && g_proc_id == 0) {
    myverbose = 1;
  }
  if ((yyin = fopen(conf_file, "rt")) == NULL){
    return(2);
  }
  yyout = fopen("/dev/null", "w");

  parse_config();  
#ifndef FIXEDVOLUME
  if(LX == 0) {
    LX = L;
  }
  if(LY == 0) {
    LY = L;
  }
  if(LZ == 0) {
    LZ = L;
  }
#endif

  if(PropInfo.format == -1) PropInfo.format = SourceInfo.format;
  g_rgi_C0 = 1. - 8.*g_rgi_C1;


  fclose(yyout);
  fclose(yyin);
  return(0);
}


/*
 * This is the function to parse the input file 
 * again. Only parameters are changed, that
 * are specified in the input file.
 * default values for paramters will not be set.
 *
 * reread_input expects the filename of the input file
 * as an input parameter.
 *
 * reread_input returns 2 if the input file did not exist 
 */

int reread_input(char * conf_file){
#ifndef FIXEDVOLUME
  int tt=T, ll=L, lx = LX, ly = LY, lz = LZ, 
      np=N_PROC_X, npy = N_PROC_Y;
#endif
  int nst=nstore;

  if(verbose && g_proc_id == 0) {
    myverbose = 1;
  }
  current_monomial = -1;
  reread = 1;

  /********************************************
   * Setting default values!
   ********************************************/

  /********************************************/

  if ((yyin = fopen(conf_file, "rt")) == NULL){
    return(2);
  }
  yyout = fopen("/dev/null", "w");

  parse_config();  

#ifndef FIXEDVOLUME
  T = tt;
  L = ll;
  LX = lx;
  LY = ly;
  LZ = lz;
  N_PROC_X = np;
  N_PROC_Y = npy;
#endif


  if(g_dbw2rand == 0) {
    g_rgi_C1 = 0.;
  }
  nstore = nst;

  g_rgi_C0 = 1. - 8.*g_rgi_C1;
  g_mu = g_mu1;

  fclose(yyout);
  fclose(yyin);
  return(0);
}

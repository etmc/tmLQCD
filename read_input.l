/***********************************************************************
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
FLTLIST [0-9"."","" "]*
EQL {SPC}*={SPC}*

%{
#ifdef HAVE_CONFIG_H
#  include<config.h>
#endif
#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include"global.h"
#include"read_input.h"
#include"default_input_values.h"
#include"monomial/monomial.h"
#include"measurements.h"
#include"integrator.h"
#include"operator.h"
#include"phmc.h"
#include<io/params.h>

inline void rmQuotes(char *str){
  char* strsave=str;

  while(*str== ' ' || *str == '\t') str++;

  if(*str=='\"') *(strsave++)=*( ++str );
  else fprintf(stderr,"Error in removing quotes from string:\n %s  \n" ,  str);
  ++str;

  while( !( *str=='\0' || *str == '\"' || *str == '\n') ) *(strsave++)=*(str++);
  *strsave='\0';
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
  int g_disable_IO_checks;
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

  int propagator_comparison;
  int nb_cores;

  int omp_num_threads;

  int nblocks_t;
  int nblocks_x;
  int nblocks_y;
  int nblocks_z;

  int dfl_field_iter;
  int dfl_poly_iter;

  int use_preconditioning;
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
%x DFLSP
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

%x REWEIGH
%x REWSAMPLES

%x INITINTEGRATOR
%x INTEGRATOR

%x GPU
%x INITGPU

%x INITOPERATOR
%x TMOP
%x DBTMOP
%x WILSONOP
%x OVERLAPOP
%x CLOVEROP
%x DBCLOVEROP
%x TMSOLVER
%x DBTMSOLVER
%x OVSOLVER

%x INITMONOMIAL
%x DETMONOMIAL
%x CLDETMONOMIAL
%x CLDETRATMONOMIAL
%x GAUGEMONOMIAL
%x NDPOLYMONOMIAL
%x NDRATMONOMIAL
%x RATMONOMIAL
%x CLRATMONOMIAL
%x RATCORMONOMIAL
%x CLRATCORMONOMIAL
%x NDCLRATMONOMIAL
%x NDRATCORMONOMIAL
%x NDCLRATCORMONOMIAL
%x POLYMONOMIAL
%x CLPOLYMONOMIAL
%x MNAME
%x MCSTR
%x MSOLVER
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


%%
^SourceFilename{EQL}               BEGIN(SOURCEFILE);
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
^DeflationSubspaceDimension{EQL}   BEGIN(DFLSP);
^GCRPreconditioner{EQL}            BEGIN(PRECON);
^ComputeReweightingFactor{EQL}     BEGIN(REWEIGH);
^NoReweightingSamples{EQL}         BEGIN(REWSAMPLES);
^SourceTimeSlice{EQL}              BEGIN(SOURCETS);
^SourceType{EQL}                   BEGIN(SOURCETYPE);
^NoSamples{EQL}                    BEGIN(NOSAMPLES);
^SplittedPropagator{EQL}           BEGIN(PROPSPLIT);
^UsePreconditioning{EQL}           BEGIN(PRECONDITIONING);

^BeginMeasurement{SPC}+            BEGIN(INITMEASUREMENT);
^ComputeModeNumber{EQL}            BEGIN(COMPUTEMN);
^ComputeTopSus{EQL}                BEGIN(COMPUTETS);
^MStarSq{EQL}                      BEGIN(MSTARSQ);
^NoSourcesZ2{EQL}                  BEGIN(NOSOURCESZ2);
^BeginMonomial{SPC}+               BEGIN(INITMONOMIAL);
^BeginInt                          BEGIN(INITINTEGRATOR);
^BeginOperator{SPC}+               BEGIN(INITOPERATOR);

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
  EndGPUInit{SPC}* {
  if(myverbose) printf("GPU parsed in line %d\n\n", line_of_file);
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
  {SPC}*SolverPrecision{EQL}{FLT} {
    sscanf(yytext, " %[2a-zA-Z] = %lf", name, &c);
    optr->eps_sq = c;
    if(myverbose) printf("  SolverPrecision set to %f line %d operator %d\n", c, line_of_file, current_operator);
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
  {SPC}*ExtraMasses{EQL}{FLTLIST} {
    char * token = NULL;
    optr->no_extra_masses = 0;
    double mass;
    if( strtok(yytext,"\n\t =,\\") != NULL ) {
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
  }
  {SPC}*ExtraMasses{EQL}{FILENAME} {
    FILE * ifs;
    double tempmass = 0.0;
    optr->no_extra_masses = 0;
    char * token = NULL;
    if( strtok(yytext,"\n\t =,\\") != NULL ) {
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

<WILSONOP,TMOP,DBTMOP>{
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

<DBTMSOLVER,TMSOLVER>{
  cg {
    optr->solver=1;
    if(myverbose) printf("  Solver set to CG line %d operator %d\n", line_of_file, current_operator);
    BEGIN(name_caller);
  }
}

<TMSOLVER>{
  mixedcg {
    optr->solver=13;
    if(myverbose) printf("  Solver set to MixedCG line %d operator %d\n", line_of_file, current_operator);
    BEGIN(name_caller);
  }
  bicgstab {
    optr->solver=0;
    if(myverbose) printf("  Solver set to BiCGstab line %d operator %d\n", line_of_file, current_operator);
    BEGIN(name_caller);
  }
  pcg {
    optr->solver=9;
    if(myverbose) printf("  Solver set to PCG line %d operator %d\n", line_of_file, current_operator);
    BEGIN(name_caller);
  }
  gmres {
    optr->solver=2;
    if(myverbose) printf("  Solver set to GMRES line %d operator %d\n", line_of_file, current_operator);
    BEGIN(name_caller);
  }
  gcr {
    optr->solver=7;
    if(myverbose) printf("  Solver set to GCR line %d operator %d\n", line_of_file, current_operator);
    BEGIN(name_caller);
  }
  gmresdr {
    optr->solver=8;
    if(myverbose) printf("  Solver set to GMRES-DR line %d operator %d\n", line_of_file, current_operator);
    BEGIN(name_caller);
  }
  cgs {
    optr->solver=3;
    if(myverbose) printf("  Solver set to CGS line %d operator %d\n", line_of_file, current_operator);
    BEGIN(name_caller);
  }
  mr {
    optr->solver=4;
    if(myverbose) printf("  Solver set to MR line %d operator %d\n", line_of_file, current_operator);
    BEGIN(name_caller);
  }
  fgmres {
    optr->solver=6;
    if(myverbose) printf("  Solver set to FGMRES line %d operator %d\n", line_of_file, current_operator);
    BEGIN(name_caller);
  }
  dflgcr {
    optr->solver=10;
    g_dflgcr_flag = 1;
    if(myverbose) printf("  Solver set to DFL-GCR line %d operator %d\n", line_of_file, current_operator);
    BEGIN(name_caller);
  }
  dflfgmres {
    optr->solver=11;
    g_dflgcr_flag = 1;
    if(myverbose) printf("  Solver set to DFL-FGMRES line %d operator %d\n", line_of_file, current_operator);
    BEGIN(name_caller);
  }
  cgmms {
    optr->solver = 12;
    if(myverbose) printf("  Solver set to CGMMS line %d operator %d\n", line_of_file, current_operator);
    BEGIN(name_caller);
  }
}

<OVSOLVER>{
  sumr {
    optr->solver = 13;
    if(myverbose) printf("  Solver set to SUMR line %d operator %d\n", line_of_file, current_operator);
    BEGIN(OVERLAPOP);
  }
  cg {
    optr->solver = 1;
    if(myverbose) printf("  Solver set to SUMR line %d operator %d\n", line_of_file, current_operator);
    BEGIN(OVERLAPOP);
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
  else if(strcmp(yytext, "DETRATIO")==0) {
    mnl->type = DETRATIO;
    strcpy((*mnl).name, "DETRATIO");
  }
  else if(strcmp(yytext, "NDDETRATIO")==0) {
    mnl->type = NDDETRATIO;
    strcpy((*mnl).name, "NDDETRATIO");
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
  else if(mnl->type == CLOVERRAT) BEGIN(CLRATMONOMIAL);
  else if(mnl->type == NDRATCOR) BEGIN(NDRATCORMONOMIAL);
  else if(mnl->type == RATCOR) BEGIN(RATCORMONOMIAL);
  else if(mnl->type == NDCLOVERRATCOR) BEGIN(NDCLRATCORMONOMIAL);
  else if(mnl->type == CLOVERRATCOR) BEGIN(CLRATCORMONOMIAL);
  else if(mnl->type == POLY || mnl->type == POLYDETRATIO)  {
          fprintf(stderr,"starting to parse poly(detratio) monomial\n");
          BEGIN(POLYMONOMIAL); 
  }
  else if(mnl->type == CLOVERDET) BEGIN(CLDETMONOMIAL);
  else if(mnl->type == CLOVERDETRATIO) BEGIN(CLDETRATMONOMIAL);
  else BEGIN(DETMONOMIAL);
}



<DETMONOMIAL,GAUGEMONOMIAL,NDPOLYMONOMIAL,POLYMONOMIAL,CLDETMONOMIAL,CLDETRATMONOMIAL,CLPOLYMONOMIAL,NDRATMONOMIAL,NDRATCORMONOMIAL,NDCLRATMONOMIAL,NDCLRATCORMONOMIAL,RATMONOMIAL,CLRATMONOMIAL,RATCORMONOMIAL,CLRATCORMONOMIAL>{
  {SPC}*Timescale{EQL}{DIGIT}+ {
    if(mnl->type == NDDETRATIO) {
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

<CLDETMONOMIAL,CLDETRATMONOMIAL,CLPOLYMONOMIAL,NDCLRATMONOMIAL,NDCLRATCORMONOMIAL,CLRATMONOMIAL,CLRATCORMONOMIAL>{
  {SPC}*CSW{EQL}{FLT} {
    sscanf(yytext, " %[a-zA-Z] = %lf", name, &c);
    mnl->c_sw = c;
    if(myverbose) printf("  CSW set to %f line %d monomial %d\n", c, line_of_file, current_monomial);
  }
}

<CLDETMONOMIAL,CLDETRATMONOMIAL>{
  {SPC}*rho{EQL}{FLT} {
    sscanf(yytext, " %[a-zA-Z] = %lf", name, &c);
    mnl->rho = c;
    if(myverbose) printf("  mass shift rho set to %f line %d monomial %d\n", c, line_of_file, current_monomial);
  }
}

<CLDETRATMONOMIAL>{
  {SPC}*rho2{EQL}{FLT} {
    sscanf(yytext, " %[a-zA-Z]2 = %lf", name, &c);
    mnl->rho2 = c;
    if(myverbose) printf("  mass shift rho2 set to %f line %d monomial %d\n", c, line_of_file, current_monomial);
  }
}

<DETMONOMIAL,POLYMONOMIAL,NDPOLYMONOMIAL,CLPOLYMONOMIAL,NDRATMONOMIAL,NDRATCORMONOMIAL,NDCLRATMONOMIAL,NDCLRATCORMONOMIAL,CLDETMONOMIAL,CLDETRATMONOMIAL,RATMONOMIAL,CLRATMONOMIAL,RATCORMONOMIAL,CLRATCORMONOMIAL>{
  {SPC}*Kappa{EQL}{FLT} {
    sscanf(yytext, " %[a-zA-Z] = %lf", name, &c);
    mnl->kappa = c;
    if(myverbose) printf("  Kappa set to %f line %d monomial %d\n", c, line_of_file, current_monomial);
  }
}

<DETMONOMIAL,POLYMONOMIAL>{
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

<NDPOLYMONOMIAL,CLPOLYMONOMIAL,NDRATMONOMIAL,NDRATCORMONOMIAL,NDCLRATMONOMIAL,NDCLRATCORMONOMIAL>{
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

<DETMONOMIAL,POLYMONOMIAL,CLDETMONOMIAL,CLDETRATMONOMIAL,NDRATMONOMIAL,NDRATCORMONOMIAL,NDCLRATMONOMIAL,NDCLRATCORMONOMIAL,RATMONOMIAL,RATCORMONOMIAL,CLRATMONOMIAL,CLRATCORMONOMIAL>{
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
}

<DETMONOMIAL,POLYMONOMIAL,CLDETMONOMIAL,CLDETRATMONOMIAL>{
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
    mnl->solver = 1;
    BEGIN(solver_caller);
  }
  bicgstab {
    if(myverbose) printf("  Solver set to \"%s\" line %d monomial %d\n", yytext, line_of_file, current_monomial);
    mnl->solver = 0;
    BEGIN(solver_caller);
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
    SourceInfo.type = 0;
    if(myverbose) printf("Using Point Sources\n");
  }
  Volume {
    SourceInfo.type = 1;
    if(myverbose) printf("Using Volume Sources\n");
  }
  TimeSlice {
    SourceInfo.type = 2;
    if(myverbose) printf("Using TimeSlice Sources\n");
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
  else {
    fprintf(stderr, "Unknown measurement type %s in line %d\n", yytext, line_of_file);
    exit(1);
  }
  /*set default frequency here, in case it is not specified
    in the input file */  
  meas->freq = _default_measurement_freq;
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
}

<ONLINEMEAS,PIONNORMMEAS,PLOOP,ORIENTEDPLAQUETTESMEAS>{
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
  if(myverbose!=0) printf("c_sw=%lf \n", g_c_sw);
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
<DFLSP>{DIGIT}+ {
  g_N_s = atoi(yytext);
  if(myverbose!=0) printf("Deflation subspace dimension set to %d \n", g_N_s);
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
<DSBLIOCHECK>yes {
  g_disable_IO_checks = 1;
  if(myverbose!=0) printf("Disable IO checks (and readback in case of Lemon IO)\n");
}
<DSBLIOCHECK>no {
  g_disable_IO_checks = 0;
  if(myverbose!=0) printf("Enable IO checks (and readback in case of Lemon IO)\n");
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
  if(myverbose != 0) printf("X0 for boundary cond. in time set to %e\n", X0);
}
<BOUNDX>{FLT} {
  X1 = atof(yytext);
  if(myverbose != 0) printf("X1 for boundary cond. in time set to %e\n", X0);
}
<BOUNDY>{FLT} {
  X2 = atof(yytext);
  if(myverbose != 0) printf("X2 for boundary cond. in time set to %e\n", X0);
}
<BOUNDZ>{FLT} {
  X3 = atof(yytext);
  if(myverbose != 0) printf("X3 for boundary cond. in time set to %e\n", X0);
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
  if(SourceInfo.basename == NULL) free(SourceInfo.basename);
  SourceInfo.basename = (char*)malloc((strlen(yytext)+1)*sizeof(char));
  strcpy(SourceInfo.basename, yytext);
  if(PropInfo.basename == NULL) free(PropInfo.basename);
  PropInfo.basename = (char*)malloc((strlen(yytext)+1)*sizeof(char));
  strcpy(PropInfo.basename, yytext);
  if(myverbose!=0) printf("source input filename set to %s\n",yytext);
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
  if(myverbose!=0) fprintf(stderr, "Compute reweighting factor\n");
}
<REWEIGH>no {
  reweighting_flag = 0;
  if(myverbose!=0) fprintf(stderr, "Do not compute reweighting factor\n");
}
<REWSAMPLES>{DIGIT}+ {
  reweighting_samples = atoi(yytext);
  if(myverbose!=0) fprintf(stderr, "Number of reweighting samples set to %d\n", reweighting_samples);
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


<INITMONOMIAL,DETMONOMIAL,CLDETMONOMIAL,CLDETRATMONOMIAL,NDPOLYMONOMIAL,NDRATMONOMIAL,NDRATCORMONOMIAL,NDCLRATMONOMIAL,NDCLRATCORMONOMIAL,CLPOLYMONOMIAL,GAUGEMONOMIAL,INTEGRATOR,INITINTEGRATOR,INITMEASUREMENT,PIONNORMMEAS,ONLINEMEAS,ORIENTEDPLAQUETTESMEAS,INITOPERATOR,TMOP,DBTMOP,OVERLAPOP,WILSONOP,CLOVEROP,DBCLOVEROP,POLYMONOMIAL,PLOOP,INITGPU,GPU,RATMONOMIAL,RATCORMONOMIAL,CLRATMONOMIAL,CLRATCORMONOMIAL>{SPC}*\n   {
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
  fprintf(stderr, "Parsing error in line %d\nAborting...!\n", line_of_file);
  fprintf(stderr, "Could not make sense out off: %s\n", yytext);
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
  g_dbw2rand = 0;
  g_running_phmc = 0;
  g_beta = _default_g_beta;
  g_N_s = _default_g_N_s;
  g_dflgcr_flag = _default_g_dflgcr_flag;
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
  if(SourceInfo.basename == NULL) SourceInfo.basename = (char*)malloc(100*sizeof(char));
  strcpy(SourceInfo.basename, _default_source_filename);
  if(PropInfo.basename == NULL) PropInfo.basename = (char*)malloc(100*sizeof(char));
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

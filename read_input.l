/*
 * $Id$
 *
 * This is the parser. (Dec 2002)
 * The .c-file is generated from .l using flex.
 * Please edit read_input.l instead of read_input.c!
 * flex should be said to be case insensitive!
 *
 * After modifiing read_input.l please call once
 * make flex_read_input
 * to update read_input.c
 *
 * Autor: Carsten Urbach
 *        urbach@physik.fu-berlin.de
 */

SPC [[:space:]]+
CMD [:][[:space:]]+
RLN [1-9(10)(11)(12)(13)(14)(15)(16)][:]
DIGIT [[:digit:]]
ZT [0-9(10)(11)]
IDXEX ("-"{DIGIT}+)
SIGN ("+"|"-")
FLT {SIGN}?{DIGIT}*+"."{DIGIT}*(e("-"|"+")?{DIGIT}+)?
FILENAME [a-zA-z0-9_][a-zA-z0-9"."_"-"]+

%{
#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include"global.h"
#include"read_input.h"
#include"default_input_values.h"

  /* Name of the parsing routine */
#define YY_DECL         int parse_config()
#define YY_NO_UNPUT

  /* declaration of input parameters */
  int line_of_file=1;
  int verbose=0;
  int startoption;
  int Ntherm;
  int Nmeas;
  int Nskip;
  int integtyp;
  int nsmall;
  int solver_flag;
  int operator_flag;
  int matrix_element_flag;
  int save_config_flag;
  int save_prop_flag;
  int save_prop_g2_flag;
  int write_cp_flag;
  int cp_interval;
  int nstore;
  int index_start, index_end;
  int random_seed;
  double dtau;
  int Nsteps;
  double q_off, q_off2;
  char rlxd_input_filename[100];
  char gauge_input_filename[100];
  int first_prop_flag;
  int max_solver_iterations;
  double solver_precision;
  int mass_number;
%}

%option never-interactive

%x BETA
%x STARTCOND
%x THERMSWEEPS
%x NMEAS
%x KAPPA
%x MU
%x MU2
%x MU3
%x SEED
%x Q1
%x Q2
%x DTAU
%x NSTEPS
%x CSW
%x INTTYP
%x NSMALL
%x NSKIP
%x RLXDINPUTFILE
%x GAUGEINPUTFILE
%x SOLVFLAG
%x OPFLAG
%x MEFLAG
%x SAVECONF
%x SAVEPROP
%x SAVEPRG2
%x WRITECP
%x CPINT
%x NSTORE
%x TT
%x LL
%x IOPROC
%x IDX
%x FPROP
%x CGMAX
%x BCGMAX
%x BOUND
%x SITER
%x SPREC
%x MNR
%x RGIC

%x COMMENT
%x ERROR

%%
^T{SPC}*={SPC}*                          BEGIN(TT);
^L{SPC}*={SPC}*                          BEGIN(LL);
^kappa{SPC}*={SPC}*                      BEGIN(KAPPA);
^mu{SPC}*={SPC}*                         BEGIN(MU);
^mu2{SPC}*={SPC}*                        BEGIN(MU2);
^mu3{SPC}*={SPC}*                        BEGIN(MU3);
^beta{SPC}*={SPC}*                       BEGIN(BETA);
^seed{SPC}*={SPC}*                       BEGIN(SEED);
^StartCondition{SPC}*={SPC}*             BEGIN(STARTCOND);
^ThermalisationSweeps{SPC}*={SPC}*       BEGIN(THERMSWEEPS);
^Measurements{SPC}*={SPC}*               BEGIN(NMEAS);
^Nskip{SPC}*={SPC}*                      BEGIN(NSKIP);
^GaugeFieldInFile{SPC}*={SPC}*           BEGIN(GAUGEINPUTFILE);
^RlxdStateInFile{SPC}*={SPC}*            BEGIN(RLXDINPUTFILE);
^SolverFlag{SPC}*={SPC}*                 BEGIN(SOLVFLAG);
^OperatorFlag{SPC}*={SPC}*               BEGIN(OPFLAG);
^ComputeMatrixElements{SPC}*={SPC}*      BEGIN(MEFLAG);
^SaveConfigurations{SPC}*={SPC}*         BEGIN(SAVECONF);
^SavePropagators{SPC}*={SPC}*            BEGIN(SAVEPROP);
^SaveGeneralizedPropagators{SPC}*={SPC}* BEGIN(SAVEPRG2);
^WriteCheckpoints{SPC}*={SPC}*           BEGIN(WRITECP);
^CheckpointInterval{SPC}*={SPC}*         BEGIN(CPINT);
^GaugeConfigInputFile{SPC}*={SPC}*       BEGIN(GAUGEINPUTFILE);
^RlxdInputFile{SPC}*={SPC}*              BEGIN(RLXDINPUTFILE);
^InitialStoreCounter{SPC}*={SPC}*        BEGIN(NSTORE);
^StdIOProcessor{SPC}*={SPC}*             BEGIN(IOPROC);
^Indices{SPC}*={SPC}*                    BEGIN(IDX);
^FirstPropagator{SPC}*={SPC}*            BEGIN(FPROP);
^CSW{SPC}*={SPC}*                        BEGIN(CSW);
^Q1{SPC}*={SPC}*                         BEGIN(Q1);
^Q2{SPC}*={SPC}*                         BEGIN(Q2);
^Integrator{SPC}*={SPC}*                 BEGIN(INTTYP);
^NSmall{SPC}*={SPC}*                     BEGIN(NSMALL);
^DTau{SPC}*={SPC}*                       BEGIN(DTAU);
^NSteps{SPC}*={SPC}*                     BEGIN(NSTEPS);
^BCGstabMaxIter{SPC}*={SPC}*             BEGIN(BCGMAX);
^CGMaxIter{SPC}*={SPC}*                  BEGIN(CGMAX);
^BoundaryCond{SPC}*={SPC}*               BEGIN(BOUND);
^BCAngleT{SPC}*={SPC}*                   BEGIN(BOUND);
^MaxSolverIterations{SPC}*={SPC}*        BEGIN(SITER);
^SolverPrecision{SPC}*={SPC}*            BEGIN(SPREC);
^MassNumber{SPC}*={SPC}*                 BEGIN(MNR);
^RGIC1{SPC}*={SPC}*                      BEGIN(RGIC);

<TT>{DIGIT}+                  {
  if(verbose!=0) printf("T =%s\n", yytext);
}
<LL>{DIGIT}+                  {
  if(verbose!=0) printf("L =%s\n", yytext);
}
<SEED>{DIGIT}+               {
  random_seed=atoi(yytext);
  if(verbose!=0) printf("seed=%s \n", yytext);
}
<KAPPA>{FLT}  {
  g_kappa=atof(yytext);
  if(verbose!=0) printf("kappa=%s \n", yytext);
}
<MU>{FLT}  {
  g_mu1=atof(yytext);
  if(verbose!=0) printf("mu=%s \n", yytext);
}
<MU2>{FLT}  {
  g_mu2=atof(yytext);
  if(verbose!=0) printf("mu2=%s \n", yytext);
}
<MU3>{FLT}  {
  g_mu3=atof(yytext);
  if(verbose!=0) printf("mu3=%s \n", yytext);
}
<BETA>{FLT}   {
  g_beta=atof(yytext);
  if(verbose!=0) printf("beta=%s \n",yytext);
}
<STARTCOND>cold {
  startoption=0; 
  if(verbose!=0) printf("Start Condition is %s \n",yytext);
}
<STARTCOND>hot {
  startoption=1;
  if(verbose!=0) printf("Start Condition is %s \n",yytext);
}
<STARTCOND>restart {
  startoption=2;
  if(verbose!=0) printf("Start Condition is %s \n",yytext);
}
<STARTCOND>continue {
  startoption=3;
  if(verbose!=0) printf("Start Condition is %s \n",yytext);
}
<THERMSWEEPS>{DIGIT}+ {
  Ntherm=atoi(yytext);
  if(verbose!=0) printf("Nterm= %s \n",yytext);
}
<NMEAS>{DIGIT}+ {
  Nmeas=atoi(yytext); 
  if(verbose!=0) printf("Nmeas= %s \n",yytext);
}
<NSKIP>{DIGIT}+ {
  Nskip=atoi(yytext);
  if(verbose!=0) printf("Nskip= %s \n",yytext);
}
<SOLVFLAG>bicgstab {
  solver_flag=0;
  if(verbose!=0) printf("Use BiCGStab Solver");
}
<SOLVFLAG>cg {
  solver_flag=1;
  if(verbose!=0) printf("Use CG Solver\n");
}
<SOLVFLAG>gmres {
  solver_flag=2;
  if(verbose!=0) printf("Use GMRES Solver\n");
}
<SOLVFLAG>cgs {
  solver_flag=3;
  if(verbose!=0) printf("Use CGS Solver\n");
}
<SOLVFLAG>mr {
  solver_flag=4;
  if(verbose!=0) printf("Use MR Solver \n");
}
<SITER>{DIGIT}+ {
  max_solver_iterations = atoi(yytext);
  if(verbose!=0) printf("Use %d iterations in the solvers!\n", max_solver_iterations);
}
<SPREC>{FLT} {
  solver_precision = atof(yytext);
  if(verbose!=0) printf("Use %e as convergence precision for the solvers!\n", solver_precision);
}
<OPFLAG>Dtm {
  operator_flag=2;
  if(verbose!=0) printf("Operator Flag is set to %s\n",yytext);
}
<OPFLAG>DWilson {
  operator_flag=1;
  if(verbose!=0) printf("Operator Flag is set to %s\n",yytext);
}
<OPFLAG>DOverlap {
  operator_flag=0;
  if(verbose!=0) printf("Operator Flag is set to %s\n",yytext);
}
<MEFLAG>yes {
  matrix_element_flag=1;
  if(verbose!=0) printf("Compute Matrix Elements: %s\n", yytext);
}
<MEFLAG>no {
  matrix_element_flag=0;
  if(verbose!=0) printf("Compute Matrix Elements: %s\n", yytext);
}
<SAVECONF>yes    {
  save_config_flag=1;
  if(verbose!=0) printf("Save configurations\n");
}
<SAVECONF>no     {
  save_config_flag=0;
  if(verbose!=0) printf("Don't save configurations\n");
}
<SAVEPROP>yes    {
  save_prop_flag=1;
  if(verbose!=0) printf("Save propagators\n");
}
<SAVEPROP>no     {
  save_prop_flag=0;
  if(verbose!=0) printf("Don't save propagators\n");
}
<SAVEPRG2>yes    {
  save_prop_g2_flag=1;
  if(verbose!=0) printf("Save generalized propagators\n");
}
<SAVEPRG2>no     {
  save_prop_g2_flag=0;
  if(verbose!=0) printf("Don't save generalized propagators\n");
}
<WRITECP>yes     {
  write_cp_flag=1;
  if(verbose!=0) printf("Write Checkpoints\n");
}
<WRITECP>no     {
  write_cp_flag=0;
  if(verbose!=0) printf("Don't write Checkpoints\n");
}
<CPINT>{DIGIT}+   {
  cp_interval=atoi(yytext);
  if(verbose!=0) printf("Write Checkpoint all %s measurements\n",yytext);
}
<RLXDINPUTFILE>{FILENAME} {
  strcpy(rlxd_input_filename,yytext);
  if(verbose!=0) printf("Ranluxd input filename set to %s\n",yytext);
}
<GAUGEINPUTFILE>{FILENAME} {
  strcpy(gauge_input_filename,yytext);
  if(verbose!=0) printf("Gauge Configuration input filename set to %s\n",yytext);
}
<NSTORE>{DIGIT}+   {
  nstore=atoi(yytext);
  if(verbose!=0) printf("Initial store counter set to %s\n",yytext);
}
<NSTORE>readin {
  nstore=-1;
  if(verbose!=0) printf("Trying to read InitialStoreCounter from file .nstore_counter\n");
}
<IOPROC>all         {
  g_stdio_proc = -1;
  if(verbose!=0) printf("All processors will give output to stdout\n");
}
<IOPROC>no          {
  g_stdio_proc = -2;
  if(verbose!=0) printf("No processor will give output to stdout\n");
}
<IOPROC>{DIGIT}+    {
  g_stdio_proc = atoi(yytext);
  if(verbose!=0) printf("processor %s will give output to stdout\n", yytext);
}
<IDX>{DIGIT}+ {
  index_start = atoi(yytext);
  index_end = index_start+1;
  if((index_start < 0)||(index_start >11)){
    printf("Error in line %d! index_start must be in [0,11]! Exiting...!\n", line_of_file);
    exit(1);
  }
  if(verbose!=0) printf("inverting for index %s\n", yytext);
}
<IDX>{IDXEX}  {
  sscanf(yytext, "-%d", &index_end);
  if((index_end < 0)||(index_end >11)){
    printf("Error in line %d! index_end must be in [0,11]! Exiting...!\n", line_of_file);
    exit(1);
  }
  if((index_end < 0)||(index_end >11)){
    printf("Warnig! index_end bigger than index_start. Will compute no propagator!\n");
  }
  if(verbose!=0) printf("inverting up to color index %d\n", index_end);
  index_end+=1;
}
<FPROP>no {
  first_prop_flag = -1;
  if(verbose!=0) printf("Do not compute the first propagator (default)\n");
}
<FPROP>compute {
  first_prop_flag = 0;
  if(verbose!=0) printf("Computing the first propagator (default)\n");
}
<FPROP>readin {
  first_prop_flag = 1;
  if(verbose!=0) printf("Reading in the first propagator\n");
}
<Q1>{FLT}   {
  q_off = atof(yytext);
  if(verbose!=0) printf("q_off set to %e\n", q_off);
}
<Q2>{FLT}   {
  q_off2 = atof(yytext);
  if(verbose!=0) printf("q_off2 set to %e\n", q_off2);
}
<INTTYP>leapfrog {
  integtyp = 1;
  if(verbose!=0) printf("Using Leap Frog integrator!\n");
}
<INTTYP>SextonWeingarten {
  integtyp = 2;
  if(verbose!=0) printf("Using SW integrator!\n");
}
<NSMALL>{DIGIT}+ {
  nsmall = atoi(yytext);
  if(verbose!=0) printf("nsmall set to %d\n", nsmall);
}
<CSW>{FLT} {
  g_c_sw = atof(yytext);
  if(verbose!=0) printf("c_sw set to %e\n", g_c_sw);
}
<DTAU>{FLT} {
  dtau = atof(yytext);
  if(verbose!=0) printf("dtau set to %e\n", dtau);
}
<NSTEPS>{DIGIT}+ {
  Nsteps = atoi(yytext);
  if(verbose!=0) printf("NSteps set to %d\n", Nsteps);
}
<BCGMAX>{DIGIT}+ {
  ITER_MAX_BCG = atoi(yytext);
  if(verbose != 0) printf("Maximal number of iterations for BCGstab set ro %d\n", ITER_MAX_BCG);
}
<CGMAX>{DIGIT}+  {
  ITER_MAX_CG = atoi(yytext);
  if(verbose != 0) printf("Maximal number of iterations for CG set ro %d\n", ITER_MAX_CG);
}
<BOUND>{FLT} {
  X0 = atof(yytext);
  if(verbose != 0) printf("X0 for boundary cond. in time set to %e\n", X0);
}
<MNR>{DIGIT}+ {
  mass_number = atoi(yytext);
  if(verbose != 0) printf("Setting mass number to %s\n", yytext);
}
<RGIC>{FLT}  {
  g_rgi_C1=atof(yytext);
  if(verbose!=0) printf("g_rgi_C1=%s \n", yytext);
}

<*>^#                       BEGIN(COMMENT);
<SEED,KAPPA,MU,MU2,MU3,BETA,STARTCOND,NMEAS,THERMSWEEPS,NSKIP,OPFLAG,SOLVFLAG,MEFLAG,SAVECONF,SAVEPROP,WRITECP,CPINT,RLXDINPUTFILE,GAUGEINPUTFILE,SAVEPRG2,NSTORE,TT,LL,IOPROC,IDX,FPROP,DTAU,Q1,Q2,NSTEPS,INTTYP,NSMALL,CSW,BCGMAX,CGMAX,BOUND,SITER,SPREC,MNR,RGIC>#      BEGIN(COMMENT);
<COMMENT>[^\n]*             {
  ;
}

<*>\n                       {
  line_of_file++;
  BEGIN(0);
}

<SEED>. {
  printf("Unknown seed in line %d.\n Must be an integer. Exiting...!\n",line_of_file);
  exit(1);
}
<KAPPA>. {
  printf("Unknown kappa in line %d.\n Must be a floating point number. Exiting...!\n",line_of_file);
  exit(1);
}
<MU>. {
  printf("Unknown mu in line %d.\n Must be a floating point number. Exiting...!\n",line_of_file);
  exit(1);
}
<MU2>. {
  printf("Unknown mu in line %d.\n Must be a floating point number. Exiting...!\n",line_of_file);
  exit(1);
}
<MU3>. {
  printf("Unknown mu in line %d.\n Must be a floating point number. Exiting...!\n",line_of_file);
  exit(1);
}
<BETA>. {
  printf("Unknown beta in line %d.\n Must be a floating point number. Exiting...!\n",line_of_file);
  exit(1);
}
<STARTCOND>.          {
  printf("Unknown Startcondition in line %d! \n Must be hot, cold, continue or restart. Exiting...!\n", line_of_file);
  exit(1);
}
<THERMSWEEPS>. {
  printf("Unknown number of TermSteps in line %d! \n Must be an integer. Exiting...!\n", line_of_file);
  exit(1);
}
<NMEAS>. {
  printf("Unknown number of MeasSteps in line %d! \n Must be an integer. Exiting...!\n", line_of_file);
  exit(1);
}
<NSKIP>. {
  printf("Unknown number of Sweeps to skip in line %d! \n Must be an integer. Exiting...!\n", line_of_file);
  exit(1);
}
<SOLVFLAG>. {
  printf("Unknown value for solver_flag in line %d! \n Must be bicgstab, cg, cgs, mr or gmres. Exiting...!\n", line_of_file);
  exit(1);
}
<OPFLAG>. {
  printf("Unknown value for operator_flag in line %d! \n Must be an integer. Exiting...!\n", line_of_file);
  exit(1);
}
<MEFLAG>. {
  printf("Unknown value for matrix_element_flag in line %d! \n Must be yes or no. Exiting...!\n", line_of_file);
  exit(1);
}
<SAVECONF>. {
  printf("Unknown value for save_config_flag in line %d! \n Must be yes or no! Exiting...!\n", line_of_file);
  exit(1);
}
<SAVEPROP>. {
  printf("Unknown value for save_prop_flag in line %d! \n Must be yes or no! Exiting...!\n", line_of_file);
  exit(1);
}
<SAVEPRG2>. {
  printf("Unknown value for save_prop_g2_flag in line %d! \n Must be yes or no! Exiting...!\n", line_of_file);
  exit(1);
}
<WRITECP>. {
  printf("Unknown value for write_checkpoint_flag in line %d! \n Must be yes or no! Exiting...!\n", line_of_file);
  exit(1);
}
<CPINT>.   {
  printf("Unknown value for checkpoint interval in line %d! \n Must be an integer! Exiting...!\n", line_of_file);
  exit(1);
}
<NSTORE>.  {
  printf("Unknown value for Initial store counter in line %d! \n Must be an integer! Exiting...!\n",line_of_file);
  exit(1);
}
<TT>.         {
  printf("Unknown value for T in line %d!\n Must be an integer value! Exiting...!\n", line_of_file);
  exit(1);
}
<LL>.         {
  printf("Unknown value for L in line %d!\n Must be an integer value! Exiting...!\n", line_of_file);
  exit(1);
}
<IOPROC>.    {
  printf("Unknown value for StdIOProcessor in line %d!\n Must be all, no or an integer value! Exiting...!\n", line_of_file);
  exit(1);
}
<IDX>. {
  printf("Error in line %d! Exiting...!\n",line_of_file);
  exit(1);
}
<FPROP>. {
  printf("Error in line %d! Must be compute or readin! Exiting...!\n", line_of_file);
  exit(1);
}
<Q1,Q2,DTAU,NSTEPS,INTTYP,NSMALL,CSW,CGMAX,BCGMAX,BOUND>. {
  printf("Error in line %d!\n", line_of_file);
  exit(1);
}
<SITER>. {
  printf("Unknown value for MaxSolverIterations in line %d! Must be an integer. Exiting...!\n", line_of_file);
  exit(1);
}
<SPREC>. {
  printf("Unknown value for SolverPrecision in line %d! Must be a floating point number. Exiting...!\n", line_of_file);
  exit(1);
}
<MNR>. {
  printf("Unknown value for MassNumber in line %d! Must be an integer. Exiting...!\n", line_of_file);
  exit(1);
}
<RGIC>. {
  printf("Unknown value for RGIC1 in line %d! Must be a floating point number. Exiting...!\n", line_of_file);
  exit(1);
}

<INITIAL>.                        BEGIN(ERROR);
<ERROR>[^\t\n]*             {
  printf("Error in line %d: %s \n",line_of_file,yytext);
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

  g_kappa = _default_g_kappa;
  g_mu = _default_g_mu;
  g_mu1 = _default_g_mu1;
  g_mu2 = _default_g_mu2;
  g_mu3 = _default_g_mu3;
  g_beta = _default_g_beta;
  g_c_sw = _default_g_c_sw;
  q_off = _default_q_off;
  q_off2 = _default_q_off2;
  dtau = _default_dtau;
  Nsteps = _default_Nsteps;
  nsmall = _default_nsmall;
  integtyp = _default_integtyp;
  random_seed = _default_random_seed;
  matrix_element_flag = _default_matrix_element_flag;
  solver_flag = _default_solver_flag;
  operator_flag = _default_operator_flag;
  startoption = _default_startoption;
  Ntherm = _default_Ntherm;
  Nmeas = _default_Nmeas;
  Nskip = _default_Nskip;
  save_config_flag = _default_save_config_flag;
  save_prop_flag = _default_save_prop_flag;
  save_prop_g2_flag = _default_save_prop_g2_flag;
  write_cp_flag = _default_write_cp_flag;
  cp_interval = _default_cp_interval;
  nstore = _default_nstore;
  strcpy(rlxd_input_filename,_default_rlxd_input_filename);
  strcpy(gauge_input_filename,_default_gauge_input_filename);
  g_stdio_proc = _default_g_stdio_proc;
  index_start = _default_index_start;
  index_end = _default_index_end;
  first_prop_flag = _default_first_prop_flag;
  ITER_MAX_CG = _default_ITER_MAX_CG;
  ITER_MAX_BCG = _default_ITER_MAX_BCG;
  X0 = _default_X0;
  max_solver_iterations = _default_max_solver_iterations;
  solver_precision = _default_solver_precision;
  mass_number = _default_mass_number;
  g_rgi_C1 = _default_g_rgi_C1;

  /********************************************/

  if ((yyin = fopen(conf_file, "rt")) == NULL){
    return(2);
  }
  yyout = fopen("/dev/null", "w");

  parse_config();  

  g_rgi_C0 = 1. - 8.*g_rgi_C1;
  g_ka_csw_8 = g_kappa*g_c_sw/8.;
  
  fclose(yyout);
  fclose(yyin);
  return(0);
}

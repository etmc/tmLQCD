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
#ifdef HAVE_CONFIG_H
#  include<config.h>
#endif
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
  int int_n[4];
  double lambda[4];
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
  double dtau, tau;
  int Nsteps;
  char rlxd_input_filename[100];
  char gauge_input_filename[100];
  int first_prop_flag;
  int max_solver_iterations;
  double solver_precision;
  int mass_number;
  int read_source_flag;
  char source_input_filename[100];
  int return_check_flag, return_check_interval;
  int source_format_flag;
  int source_time_slice;
%}

%option never-interactive

%x BETA
%x STARTCOND
%x THERMSWEEPS
%x NMEAS
%x KAPPA
%x MUBAR
%x EPSBAR
%x MU
%x MU2
%x MU3
%x SEED
%x Q1
%x Q2
%x DTAU
%x TAU
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
%x LLX
%x LLY
%x LLZ
%x NPROCX
%x NPROCY
%x NPROCZ
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
%X READSOURCE
%x SOURCEFORMAT
%x SOURCEFILE
%x SOURCETS
%x INT0
%x INT1
%x INT2
%x INT3
%x INT4
%x LAMBDA0
%x LAMBDA1
%x LAMBDA2
%x LAMBDA3
%x LAMBDA4
%x RELPREC
%x FORCEPREC
%x FORCEPREC1
%x FORCEPREC2
%x FORCEPREC3
%x ACCPREC
%x ACCPREC1
%x ACCPREC2
%x ACCPREC3
%x REVCHECK
%x REVINT
%x DEBUG
%x CSGN1
%x CSGN2
%x CSGN3


%x COMMENT
%x ERROR

%%
^T{SPC}*={SPC}*                          BEGIN(TT);
^L{SPC}*={SPC}*                          BEGIN(LL);
^LX{SPC}*={SPC}*                         BEGIN(LLX);
^LY{SPC}*={SPC}*                         BEGIN(LLY);
^LZ{SPC}*={SPC}*                         BEGIN(LLZ);
^NRXProcs{SPC}*={SPC}*                   BEGIN(NPROCX);
^NRYProcs{SPC}*={SPC}*                   BEGIN(NPROCY);
^NRZProcs{SPC}*={SPC}*                   BEGIN(NPROCZ);
^kappa{SPC}*={SPC}*                      BEGIN(KAPPA);
^mubar{SPC}*={SPC}*                      BEGIN(MUBAR);
^epsbar{SPC}*={SPC}*                     BEGIN(EPSBAR);
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
^Tau{SPC}*={SPC}*                        BEGIN(TAU);
^NSteps{SPC}*={SPC}*                     BEGIN(NSTEPS);
^BCGstabMaxIter{SPC}*={SPC}*             BEGIN(BCGMAX);
^CGMaxIter{SPC}*={SPC}*                  BEGIN(CGMAX);
^BoundaryCond{SPC}*={SPC}*               BEGIN(BOUND);
^BCAngleT{SPC}*={SPC}*                   BEGIN(BOUND);
^MaxSolverIterations{SPC}*={SPC}*        BEGIN(SITER);
^SolverPrecision{SPC}*={SPC}*            BEGIN(SPREC);
^MassNumber{SPC}*={SPC}*                 BEGIN(MNR);
^RGIC1{SPC}*={SPC}*                      BEGIN(RGIC);
^ReadSource{SPC}*={SPC}*                 BEGIN(READSOURCE);
^SourceInputFilename{SPC}*={SPC}*        BEGIN(SOURCEFILE);
^SourceFormat{SPC}*={SPC}*               BEGIN(SOURCEFORMAT);
^SourceTimeSlice{SPC}*={SPC}*            BEGIN(SOURCETS);
^ExtIntStepsGauge{SPC}*={SPC}*           BEGIN(INT0);
^N_G{SPC}*={SPC}*                        BEGIN(INT0);
^ExtIntStepsMu0{SPC}*={SPC}*             BEGIN(INT1);
^IntegrationStepsMu{SPC}*={SPC}*         BEGIN(INT1);
^ExtIntStepsMu1{SPC}*={SPC}*             BEGIN(INT2);
^IntegrationStepsMu2{SPC}*={SPC}*        BEGIN(INT2);
^ExtIntStepsMu2{SPC}*={SPC}*             BEGIN(INT3);
^IntegrationStepsMu3{SPC}*={SPC}*        BEGIN(INT3);
^ExtIntStepsMu3{SPC}*={SPC}*             BEGIN(INT4);
^IntegrationStepsMu4{SPC}*={SPC}*        BEGIN(INT4);
^2MNLambdaGauge{SPC}*={SPC}*             BEGIN(LAMBDA0);
^2MNLambdaMu{SPC}*={SPC}*                BEGIN(LAMBDA1);
^2MNLambdaMu2{SPC}*={SPC}*               BEGIN(LAMBDA2);
^2MNLambdaMu3{SPC}*={SPC}*               BEGIN(LAMBDA3);
^2MNLambdaMu4{SPC}*={SPC}*               BEGIN(LAMBDA4);
^UseRelativePrecision{SPC}*={SPC}*       BEGIN(RELPREC);
^ForcePrecision{SPC}*={SPC}*             BEGIN(FORCEPREC);
^ForcePrecisionMu{SPC}*={SPC}*           BEGIN(FORCEPREC1);
^ForcePrecisionMu2{SPC}*={SPC}*          BEGIN(FORCEPREC2);
^ForcePrecisionMu3{SPC}*={SPC}*          BEGIN(FORCEPREC3);
^AcceptancePrecision{SPC}*={SPC}*        BEGIN(ACCPREC);
^AcceptancePrecisionMu{SPC}*={SPC}*      BEGIN(ACCPREC1);
^AcceptancePrecisionMu2{SPC}*={SPC}*     BEGIN(ACCPREC2);
^AcceptancePrecisionMu3{SPC}*={SPC}*     BEGIN(ACCPREC3);
^ReversibilityCheck{SPC}*={SPC}*         BEGIN(REVCHECK);
^ReversibilityCheckIntervall{SPC}*={SPC}* BEGIN(REVINT);
^DebugLevel{SPC}*={SPC}*                 BEGIN(DEBUG);
^CSGHistMu{SPC}*={SPC}*                  BEGIN(CSGN1);
^CSGHistoryMu{SPC}*={SPC}*               BEGIN(CSGN1);
^CSGHistMu2{SPC}*={SPC}*                 BEGIN(CSGN2);
^CSGHistoryMu2{SPC}*={SPC}*              BEGIN(CSGN2);
^CSGHistMu3{SPC}*={SPC}*                 BEGIN(CSGN3);
^CSGHistoryMu3{SPC}*={SPC}*              BEGIN(CSGN3);

<TT>{DIGIT}+                  {
#ifndef FIXEDVOLUME
  T_global = atoi(yytext);
#endif
  if(verbose!=0) printf("T =%s\n", yytext);
}
<LL>{DIGIT}+                  {
#ifndef FIXEDVOLUME
  L = atoi(yytext);
#endif
  if(verbose!=0) printf("L =%s\n", yytext);
}
<LLX>{DIGIT}+                  {
#ifndef FIXEDVOLUME
  LX = atoi(yytext);
#endif
  if(verbose!=0) printf("LX =%s\n", yytext);
}
<LLY>{DIGIT}+                  {
#ifndef FIXEDVOLUME
  LY = atoi(yytext);
#endif
  if(verbose!=0) printf("LY =%s\n", yytext);
}
<LLZ>{DIGIT}+                  {
#ifndef FIXEDVOLUME
  LZ = atoi(yytext);
#endif
  if(verbose!=0) printf("LZ =%s\n", yytext);
}
<NPROCX>{DIGIT}+              {
#ifndef FIXEDVOLUME
  N_PROC_X = atoi(yytext);
#endif
  if(verbose!=0) printf("Nr of processors in x direction = %s\n", yytext);
}
<NPROCY>{DIGIT}+              {
#ifndef FIXEDVOLUME
  N_PROC_Y = atoi(yytext);
#endif
  if(verbose!=0) printf("Nr of processors in y direction = %s\n", yytext);
}
<NPROCZ>{DIGIT}+              {
#ifndef FIXEDVOLUME
  N_PROC_Z = atoi(yytext);
#endif
  if(verbose!=0) printf("Nr of processors in z direction = %s\n", yytext);
}
<SEED>{DIGIT}+               {
  random_seed=atoi(yytext);
  if(verbose!=0) printf("seed=%s \n", yytext);
}
<KAPPA>{FLT}  {
  g_kappa=atof(yytext);
  if(verbose!=0) printf("kappa=%s \n", yytext);
}
<MUBAR>{FLT}  {
  g_mubar=atof(yytext);
  if(verbose!=0) printf("mubar=%s \n", yytext);
}
<EPSBAR>{FLT}  {
  g_epsbar=atof(yytext);
  if(verbose!=0) printf("epsbar=%s \n", yytext);
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
<INTTYP>leapfrog {
  integtyp = 1;
  if(verbose!=0) printf("Using Leap Frog integrator!\n");
}
<INTTYP>SextonWeingarten {
  integtyp = 2;
  if(verbose!=0) printf("Using SW integrator!\n");
}
<INTTYP>ExtLeapFrog {
  integtyp = 3;
  if(verbose!=0) printf("Using multiple time scale Leapfrog integrator!\n");
}
<INTTYP>ExtSextonWeingarten {
  integtyp = 4;
  if(verbose!=0) printf("Using multiple time scale Sexton-Weingarten integrator!\n");
}
<INTTYP>ImprLeapFrog {
  integtyp = 5;
  if(verbose!=0) printf("Using higher order Leapfrog integrator!\n");
}
<INTTYP>2MN {
  integtyp = 6;
  if(verbose!=0) printf("Using Second order Minimal norm integrator!\n");
}
<INTTYP>2MNposition {
  integtyp = 7;
  if(verbose!=0) printf("Using Second order Minimal norm integrator (position version)!\n");
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
<TAU>{FLT} {
  tau = atof(yytext);
  if(verbose!=0) printf("tau set to %e\n", tau);
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
<READSOURCE>yes     {
  read_source_flag=1;
  if(verbose!=0) printf("Read inversion source from file\n");
}
<READSOURCE>no     {
  read_source_flag=0;
  if(verbose!=0) printf("Don't read inversion source from file\n");
}
<SOURCEFILE>{FILENAME} {
  strcpy(source_input_filename,yytext);
  if(verbose!=0) printf("source input filename set to %s\n",yytext);
}
<SOURCEFORMAT>cmi      {
  source_format_flag = 1;
  if(verbose!=0) printf("Using CM format for source input file\n");
}
<SOURCETS>{DIGIT}+   {
  source_time_slice = atoi(yytext);
  if(verbose!=0) printf("Using only timeslice %s of the source, padding the rest with zeros\n", yytext);
}
<INT0>{DIGIT}+ {
  int_n[0] = atoi(yytext);
  if(verbose!=0) printf("Number of steps in ExtLeapFrog integrator for gauge set to %d!\n", int_n[0]);
}
<INT1>{DIGIT}+ {
  int_n[1] = atoi(yytext);
  if(verbose!=0) printf("Number of steps in ExtLeapFrog integrator for psf 1 (mu) set to %d!\n", int_n[1]);
}
<INT2>{DIGIT}+ {
  int_n[2] = atoi(yytext);
  if(verbose!=0) printf("Number of steps in ExtLeapFrog integrator for psf 2 (mu2) set to %d!\n", int_n[2]);
}
<INT3>{DIGIT}+ {
  int_n[3] = atoi(yytext);
  if(verbose!=0) printf("Number of steps in ExtLeapFrog integrator for psf 3 (mu3) set to %d!\n", int_n[3]);
}
<INT4>{DIGIT}+ {
  if(verbose!=0) printf("Number of steps in ExtLeapFrog integrator for psf 4 (mu4) set to %d!\n", int_n[1]);
}
<LAMBDA0>{FLT} {
  lambda[0] = atof(yytext);
  if(verbose!=0) printf("Set lambda parameter for gauge fields (in the 2MN integrator) to %f!\n", lambda[0]);
}
<LAMBDA1>{FLT} {
  lambda[1] = atof(yytext);
  if(verbose!=0) printf("Set lambda parameter for psf 1 (in the 2MN integrator) to %f!\n", lambda[0]);
}
<LAMBDA2>{FLT} {
  lambda[2] = atof(yytext);
  if(verbose!=0) printf("Set lambda parameter for psf 2 (in the 2MN integrator) to %f!\n", lambda[0]);
}
<LAMBDA3>{FLT} {
  lambda[3] = atof(yytext);
  if(verbose!=0) printf("Set lambda parameter for psf 3 (in the 2MN integrator) to %f!\n", lambda[0]);
}
<LAMBDA4>{FLT} {
  if(verbose!=0) printf("Set lambda parameter for psf 4 (in the 2MN integrator) to %f! (not yet implemented)\n", lambda[0]);
}
<FORCEPREC>{FLT}  {
  g_eps_sq_force=atof(yytext);
  if(verbose!=0) printf("g_eps_sq_force=%s Residual for inversions in the force computation\n", yytext);
}
<FORCEPREC1>{FLT}  {
  g_eps_sq_force1=atof(yytext);
  if(verbose!=0) printf("g_eps_sq_force(mu)=%s Residual for inversions in the force computation\n", yytext);
}
<FORCEPREC2>{FLT}  {
  g_eps_sq_force2=atof(yytext);
  if(verbose!=0) printf("g_eps_sq_force(mu2)=%s Residual for inversions in the force computation\n", yytext);
}
<FORCEPREC3>{FLT}  {
  g_eps_sq_force3=atof(yytext);
  if(verbose!=0) printf("g_eps_sq_force(mu3)=%s Residual for inversions in the force computation\n", yytext);
}
<ACCPREC>{FLT}  {
  g_eps_sq_acc=atof(yytext);
  if(verbose!=0) printf("g_eps_sq_acc=%s Residual for inversions in the acceptance step\n", yytext);
}
<ACCPREC1>{FLT}  {
  g_eps_sq_acc1=atof(yytext);
  if(verbose!=0) printf("g_eps_sq_acc(mu)=%s Residual for inversions in the acceptance step\n", yytext);
}
<ACCPREC2>{FLT}  {
  g_eps_sq_acc2=atof(yytext);
  if(verbose!=0) printf("g_eps_sq_acc(mu2)=%s Residual for inversions in the acceptance step\n", yytext);
}
<ACCPREC3>{FLT}  {
  g_eps_sq_acc3=atof(yytext);
  if(verbose!=0) printf("g_eps_sq_acc(mu3)=%s Residual for inversions in the acceptance step\n", yytext);
}
<RELPREC>yes  {
  g_relative_precision_flag = 1;
  if(verbose!=0) printf("Using relative precision\n");
}
<RELPREC>no  {
  g_relative_precision_flag = 0;
  if(verbose!=0) printf("Using absolute precision\n");
}
<REVCHECK>yes {
  return_check_flag = 1;
  if(verbose!=0) printf("Perform checks of Reversibility\n");
}
<REVCHECK>no {
  return_check_flag = 0;
  if(verbose!=0) printf("Don't perform checks of Reversibility\n");
}
<REVINT>{DIGIT}+ {
  return_check_interval = atoi(yytext);
  if(verbose!=0) printf("Check reversibility all %d trajectories\n", return_check_interval);
}
<DEBUG>{DIGIT}+ {
  g_debug_level = atoi(yytext);
  if(verbose!=0) printf("Debug level = %d\n", g_debug_level);
}
<CSGN1>{DIGIT}+ {
  g_csg_N[0] = atoi(yytext);
  if(verbose!=0) printf("Chronological Invertier history length for mu set to %d\n", g_csg_N[0]);
}
<CSGN2>{DIGIT}+ {
  g_csg_N[2] = atoi(yytext);
  if(verbose!=0) printf("Chronological Invertier history length for mu set to %d\n", g_csg_N[2]);
}
<CSGN3>{DIGIT}+ {
  g_csg_N[4] = atoi(yytext);
  if(verbose!=0) printf("Chronological Invertier history length for mu set to %d\n", g_csg_N[4]);
}

<*>^#                       BEGIN(COMMENT);
<SEED,KAPPA,MUBAR,EPSBAR,MU,MU2,MU3,BETA,STARTCOND,NMEAS,THERMSWEEPS,NSKIP,OPFLAG,SOLVFLAG>#   BEGIN(COMMENT);
<MEFLAG,SAVECONF,SAVEPROP,WRITECP,CPINT,RLXDINPUTFILE,GAUGEINPUTFILE,SAVEPRG2>#   BEGIN(COMMENT);
<NSTORE,TT,LL,LLX,LLY,LLZ,NPROCX,NPROCY,NPROCZ,IOPROC,IDX,FPROP,DTAU,TAU,Q1,Q2,NSTEPS,INTTYP>#   BEGIN(COMMENT);
<NSMALL,CSW,BCGMAX,CGMAX,BOUND,SITER,SPREC,MNR,RGIC,READSOURCE,SOURCEFILE,INT0>#   BEGIN(COMMENT);
<INT1,INT2,INT3,INT4,RELPREC,ACCPREC,FORCEPREC,FORCEPREC1,FORCEPREC2,FORCEPREC3>#   BEGIN(COMMENT);
<ACCPREC1,ACCPREC2,ACCPREC3,REVCHECK,REVINT,DEBUG,CSGN1,CSGN2,CSGN3,SOURCEFORMAT>#      BEGIN(COMMENT);
<SOURCETS>#   BEGIN(COMMENT);
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
<MUBAR>. {
  printf("Unknown mubar in line %d.\n Must be a floating point number. Exiting...!\n",line_of_file);
  exit(1);
}
<EPSBAR>. {
  printf("Unknown epsbar in line %d.\n Must be a floating point number. Exiting...!\n",line_of_file);
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
<LLX>.         {
  printf("Unknown value for LX in line %d!\n Must be an integer value! Exiting...!\n", line_of_file);
  exit(1);
}
<LLY>.         {
  printf("Unknown value for LY in line %d!\n Must be an integer value! Exiting...!\n", line_of_file);
  exit(1);
}
<LLZ>.         {
  printf("Unknown value for LZ in line %d!\n Must be an integer value! Exiting...!\n", line_of_file);
  exit(1);
}
<NPROCX>.    {
  printf("Unknown value for NRXProcs in line %d!\n Must be an integer value! Exiting...!\n", line_of_file);
  exit(1);
}
<NPROCY>.    {
  printf("Unknown value for NRYProcs in line %d!\n Must be an integer value! Exiting...!\n", line_of_file);
  exit(1);
}
<NPROCZ>.    {
  printf("Unknown value for NRYProcs in line %d!\n Must be an integer value! Exiting...!\n", line_of_file);
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
<Q1,Q2,DTAU,TAU,NSTEPS,INTTYP,NSMALL,CSW,CGMAX,BCGMAX,BOUND,SOURCEFILE,READSOURCE,INT0,INT1,INT2,INT3,INT4,REVCHECK,REVINT,DEBUG,SOURCEFORMAT,SOURCETS>. {
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
<RELPREC>. {
  printf("Should be yes or no for relative precision in line %d! Exiting...!\n", line_of_file);
  exit(1);
}
<FORCEPREC,FORCEPREC1,FORCEPREC2,FORCEPREC3>. {
  printf("Unknown value for ForcePrecision in line %d! Must be a floating point number. Exiting...!\n", line_of_file);
  exit(1);
}
<ACCPREC,ACCPREC1,ACCPREC2,ACCPREC3>. {
  printf("Unknown value for AcceptancePrecision in line %d! Must be a floating point number. Exiting...!\n", line_of_file);
  exit(1);
}
<CSGN1,CSGN2,CSGN3>. {
  printf("Unknown value for CSGHistMu in line %d! Must be an integer number. Exiting...!\n", line_of_file);
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
  g_kappa = _default_g_kappa;
  g_mubar = _default_g_mubar;
  g_epsbar = _default_g_epsbar;
  g_mu = _default_g_mu;
  g_mu1 = _default_g_mu1;
  g_mu2 = _default_g_mu2;
  g_mu3 = _default_g_mu3;
  g_beta = _default_g_beta;
  g_c_sw = _default_g_c_sw;
  dtau = _default_dtau;
  tau = _default_tau;
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
  strcpy(rlxd_input_filename, _default_rlxd_input_filename);
  strcpy(gauge_input_filename, _default_gauge_input_filename);
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
  read_source_flag= _default_read_source_flag;
  strcpy(source_input_filename, _default_source_filename);
  g_eps_sq_force = _default_g_eps_sq_force;
  g_eps_sq_acc = _default_g_eps_sq_acc;
  g_eps_sq_force1 = _default_g_eps_sq_force1;
  g_eps_sq_acc1 = _default_g_eps_sq_acc1;
  g_eps_sq_force2 = _default_g_eps_sq_force2;
  g_eps_sq_acc2 = _default_g_eps_sq_acc2;
  g_eps_sq_force3 = _default_g_eps_sq_force3;
  g_eps_sq_acc3 = _default_g_eps_sq_acc3;
  g_relative_precision_flag = _default_g_relative_precision_flag;
  return_check_flag = _default_return_check_flag;
  return_check_interval = _default_return_check_interval;
  g_debug_level = _default_g_debug_level;
  g_csg_N[0] = _default_g_csg_N;
  g_csg_N[2] = _default_g_csg_N;
  g_csg_N[4] = _default_g_csg_N;
  g_csg_N[6] = _default_g_csg_N;
  lambda[0] = _default_2mn_lambda;
  lambda[1] = _default_2mn_lambda;
  lambda[2] = _default_2mn_lambda;
  lambda[3] = _default_2mn_lambda;
  source_format_flag = _default_source_format_flag;
  source_time_slice = _default_source_time_slice;

  /********************************************/

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
  
  if(g_eps_sq_force1 < 0) g_eps_sq_force1 = g_eps_sq_force;
  if(g_eps_sq_force2 < 0) g_eps_sq_force2 = g_eps_sq_force;
  if(g_eps_sq_force3 < 0) g_eps_sq_force3 = g_eps_sq_force;
  if(g_eps_sq_acc1 < 0) g_eps_sq_acc1 = g_eps_sq_acc;
  if(g_eps_sq_acc2 < 0) g_eps_sq_acc2 = g_eps_sq_acc;
  if(g_eps_sq_acc3 < 0) g_eps_sq_acc3 = g_eps_sq_acc;


  g_rgi_C0 = 1. - 8.*g_rgi_C1;
  g_ka_csw_8 = g_kappa*g_c_sw/8.;
  
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
  int nst=nstore, j=0;
  double m2 = g_mu2, m3 = g_mu3;
  int n1 = g_csg_N[0], n2 = g_csg_N[2], n3 = g_csg_N[4], n4 = g_csg_N[6];
  double x;

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
  g_csg_N[0] = n1;
  g_csg_N[2] = n2;
  g_csg_N[4] = n3;
  g_csg_N[6] = n4;


  if(g_dbw2rand == 0) {
    g_rgi_C1 = 0.;
  }
  nstore = nst;

  g_rgi_C0 = 1. - 8.*g_rgi_C1;
  g_ka_csw_8 = g_kappa*g_c_sw/8.;
  
  if(g_mu3 > 0. && g_mu3 != m3) {
    g_mu = g_mu1;
    g_mu1 = g_mu3;
    g_mu3 = g_mu;

    j = int_n[1];
    int_n[1] = int_n[3];
    int_n[3] = j;

    x = lambda[1];
    lambda[1] = lambda[3];
    lambda[3] = x;

    g_nr_of_psf = 3;
  }
  else if(g_mu2 > 0. && g_mu2 != m2) {
    g_mu = g_mu1;
    g_mu1 = g_mu2;
    g_mu2 = g_mu;

    int_n[3] = int_n[1];
    int_n[1] = int_n[2];
    int_n[2] = int_n[3];

    lambda[3] = lambda[1];
    lambda[1] = lambda[2];
    lambda[2] = lambda[3];

    g_nr_of_psf = 2;
  }
  for(j = 0; j < g_nr_of_psf+1; j++) {
    if(int_n[j] == 0) int_n[j] = 1;
  }
  if(g_nr_of_psf == 3) {
    g_eps_sq_force = g_eps_sq_force1;
    g_eps_sq_force1 = g_eps_sq_force3;
    g_eps_sq_force3 = g_eps_sq_force;
    g_eps_sq_acc = g_eps_sq_acc1;
    g_eps_sq_acc1 = g_eps_sq_acc3;
    g_eps_sq_acc3 = g_eps_sq_acc;
  }
  if(g_nr_of_psf == 2) {
    g_eps_sq_force = g_eps_sq_force1;
    g_eps_sq_force1 = g_eps_sq_force2;
    g_eps_sq_force2 = g_eps_sq_force;
    g_eps_sq_acc = g_eps_sq_acc1;
    g_eps_sq_acc1 = g_eps_sq_acc2;
    g_eps_sq_acc2 = g_eps_sq_acc;
  }
  g_mu = g_mu1;
  g_eps_sq_acc = g_eps_sq_acc1;
  g_eps_sq_force = g_eps_sq_force1;

  fclose(yyout);
  fclose(yyin);
  return(0);
}

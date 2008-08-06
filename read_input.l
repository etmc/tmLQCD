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

SPC [[:blank:]]+
CMD [:][[:space:]]+
RLN [1-9(10)(11)(12)(13)(14)(15)(16)][:]
DIGIT [[:digit:]]
ZT [0-9(10)(11)]
IDXEX ("-"{DIGIT}+)
SIGN ("+"|"-")
FLT {SIGN}?{DIGIT}*+"."{DIGIT}*(e("-"|"+")?{DIGIT}+)?
FILENAME [a-zA-Z0-9_".""-""/"][a-zA-z0-9"."_"-""/"]+
NAME [a-zA-Z0-9_]+
TYPE [0-9A-Z]+

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
#include"monomial.h"
#include"integrator.h"

  /* Name of the parsing routine */
#define YY_DECL         int parse_config()
#define YY_NO_UNPUT

  /* declaration of input parameters */
  int i=0;
  int line_of_file=1;
  int current_monomial=-1;
  extern int no_monomials;
  monomial * mnl;
  int comment_caller;
  int a,b;
  float c;
  int reread = 0;
  char name[100];
  char * type;

  int verbose = 0;
  int myverbose = 0;
  int startoption;
  int Ntherm;
  int Nmeas;
  int Nskip;
  int integtyp;
  int int_n[4];
  double lambda[4];
  int nsmall;
  int solver_flag;
  int gmres_m_parameter, gmresdr_nr_ev;
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
  int gauge_precision_read_flag;
  int gauge_precision_write_flag;
  int prop_precision_flag;
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
  int propagator_splitted;
  int source_splitted;
  int source_location;
  int no_eigenvalues;
  double eigenvalue_precision;
  int sub_evs_cg_flag;
  int even_odd_flag;
  int write_prop_format_flag;
  int online_measurement_flag;
  int online_measurement_freq;
%}

%option never-interactive

%x BETA
%x STARTCOND
%x THERMSWEEPS
%x NMEAS
%x KAPPA
%x ACCPTILDE
%x ACCHFIN
%x RECEV
%x MUBAR
%x EPSBAR
%x MU
%x SEED
%x Q1
%x Q2
%x NSKIP
%x RLXDINPUTFILE
%x GAUGEINPUTFILE
%x GAUGERPREC
%x GAUGEWPREC
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
%x PHMCFLAV
%x COMPUTEEVS
%x PCOMPUTEEVS
%x PPP
%x SMAX
%x SMIN
%x DEGP
%x SPLITPROP
%x SPLITSOURCE
%x SRCLOC
%x SUBEVCG
%x NOEV
%x PRECEV
%x HEAVYTS
%x EO
%x WRPROPFLAG
%x PROPPREC
%x PROPTYPE
%x ONMEAS
%x ONFREQ

%x INITINTEGRATOR
%x INTEGRATOR

%x INITMONOMIAL
%x MONOMIAL
%x MNAME
%x MSOLVER


%x COMMENT
%x ERROR

%%
^T{SPC}*={SPC}*                            BEGIN(TT);
^L{SPC}*={SPC}*                            BEGIN(LL);
^LX{SPC}*={SPC}*                           BEGIN(LLX);
^LY{SPC}*={SPC}*                           BEGIN(LLY);
^LZ{SPC}*={SPC}*                           BEGIN(LLZ);
^NRXProcs{SPC}*={SPC}*                     BEGIN(NPROCX);
^NRYProcs{SPC}*={SPC}*                     BEGIN(NPROCY);
^NRZProcs{SPC}*={SPC}*                     BEGIN(NPROCZ);
^kappa{SPC}*={SPC}*                        BEGIN(KAPPA);
^2KappaMu{SPC}*={SPC}*                     BEGIN(MU);
^mubar{SPC}*={SPC}*                        BEGIN(MUBAR);
^Phmc2KappaMubar{SPC}*={SPC}*              BEGIN(MUBAR);
^PhmcExactPolynomial{SPC}*={SPC}*          BEGIN(PPP);
^epsbar{SPC}*={SPC}*                       BEGIN(EPSBAR);
^Phmc2KappaEpsBar{SPC}*={SPC}*             BEGIN(EPSBAR);
^beta{SPC}*={SPC}*                         BEGIN(BETA);
^Acc_Ptilde{SPC}*={SPC}*                   BEGIN(ACCPTILDE);
^PhmcPrecisionPtilde{SPC}*={SPC}*          BEGIN(ACCPTILDE);
^Acc_Hfin{SPC}*={SPC}*                     BEGIN(ACCHFIN);
^PhmcPrecisionHfin{SPC}*={SPC}*            BEGIN(ACCHFIN);
^recev{SPC}*={SPC}*                        BEGIN(RECEV);
^PhmcRecEVInterval{SPC}*={SPC}*            BEGIN(RECEV);
^NoEigenvalues{SPC}*={SPC}*                BEGIN(NOEV);
^EigenvaluePrecision{SPC}*={SPC}*          BEGIN(PRECEV);
^seed{SPC}*={SPC}*                         BEGIN(SEED);
^StartCondition{SPC}*={SPC}*               BEGIN(STARTCOND);
^ThermalisationSweeps{SPC}*={SPC}*         BEGIN(THERMSWEEPS);
^Measurements{SPC}*={SPC}*                 BEGIN(NMEAS);
^Nskip{SPC}*={SPC}*                        BEGIN(NSKIP);
^GaugeFieldInFile{SPC}*={SPC}*             BEGIN(GAUGEINPUTFILE);
^RlxdStateInFile{SPC}*={SPC}*              BEGIN(RLXDINPUTFILE);
^SolverFlag{SPC}*={SPC}*                   BEGIN(SOLVFLAG);
^SubtractEVForCG{SPC}*={SPC}*              BEGIN(SUBEVCG);
^OperatorFlag{SPC}*={SPC}*                 BEGIN(OPFLAG);
^ComputeMatrixElements{SPC}*={SPC}*        BEGIN(MEFLAG);
^SaveConfigurations{SPC}*={SPC}*           BEGIN(SAVECONF);
^SavePropagators{SPC}*={SPC}*              BEGIN(SAVEPROP);
^SaveGeneralizedPropagators{SPC}*={SPC}*   BEGIN(SAVEPRG2);
^WriteCheckpoints{SPC}*={SPC}*             BEGIN(WRITECP);
^CheckpointInterval{SPC}*={SPC}*           BEGIN(CPINT);
^GaugeConfigInputFile{SPC}*={SPC}*         BEGIN(GAUGEINPUTFILE);
^RlxdInputFile{SPC}*={SPC}*                BEGIN(RLXDINPUTFILE);
^InitialStoreCounter{SPC}*={SPC}*          BEGIN(NSTORE);
^StdIOProcessor{SPC}*={SPC}*               BEGIN(IOPROC);
^Indices{SPC}*={SPC}*                      BEGIN(IDX);
^FirstPropagator{SPC}*={SPC}*              BEGIN(FPROP);
^Q1{SPC}*={SPC}*                           BEGIN(Q1);
^Q2{SPC}*={SPC}*                           BEGIN(Q2);
^BCGstabMaxIter{SPC}*={SPC}*               BEGIN(BCGMAX);
^CGMaxIter{SPC}*={SPC}*                    BEGIN(CGMAX);
^BoundaryCond{SPC}*={SPC}*                 BEGIN(BOUND);
^BCAngleT{SPC}*={SPC}*                     BEGIN(BOUND);
^MaxSolverIterations{SPC}*={SPC}*          BEGIN(SITER);
^SolverPrecision{SPC}*={SPC}*              BEGIN(SPREC);
^MassNumber{SPC}*={SPC}*                   BEGIN(MNR);
^ReadSource{SPC}*={SPC}*                   BEGIN(READSOURCE);
^SourceInputFilename{SPC}*={SPC}*          BEGIN(SOURCEFILE);
^SourceFormat{SPC}*={SPC}*                 BEGIN(SOURCEFORMAT);
^SourceTimeSlice{SPC}*={SPC}*              BEGIN(SOURCETS);
^UseRelativePrecision{SPC}*={SPC}*         BEGIN(RELPREC);
^ReversibilityCheck{SPC}*={SPC}*           BEGIN(REVCHECK);
^ReversibilityCheckIntervall{SPC}*={SPC}*  BEGIN(REVINT);
^DebugLevel{SPC}*={SPC}*                   BEGIN(DEBUG);
^GMRESMParameter{SPC}*={SPC}*              BEGIN(GMRESM);
^GMRESDRNrEv{SPC}*={SPC}*                  BEGIN(GMRESDRNEV);
^GaugeConfigReadPrecision{SPC}*={SPC}*     BEGIN(GAUGERPREC);
^GaugeConfigWritePrecision{SPC}*={SPC}*    BEGIN(GAUGEWPREC);
^PropagatorPrecision{SPC}*={SPC}*          BEGIN(PROPPREC);
^ReproduceRandomNumbers{SPC}*={SPC}*       BEGIN(REPRORND);
^UseSloppyPrecision{SPC}*={SPC}*           BEGIN(SLOPPYPREC);
^UseStoutSmearing{SPC}*={SPC}*             BEGIN(USESTOUT);
^StoutRho{SPC}*={SPC}*                     BEGIN(STOUTRHO);
^StoutNoIterations{SPC}*={SPC}*            BEGIN(STOUTITER);
^PhmcNoFlavours{SPC}*={SPC}*               BEGIN(PHMCFLAV);
^PhmcComputeOnlyEVs{SPC}*={SPC}*           BEGIN(PCOMPUTEEVS);
^ComputeEVs{SPC}*={SPC}*                   BEGIN(COMPUTEEVS);
^PhmcStildeMax{SPC}*={SPC}*                BEGIN(SMAX);
^PhmcStildeMin{SPC}*={SPC}*                BEGIN(SMIN);
^PhmcDegreeOfP{SPC}*={SPC}*                BEGIN(DEGP);
^SplittedPropagator{SPC}*={SPC}*           BEGIN(SPLITPROP);
^SplittedSource{SPC}*={SPC}*               BEGIN(SPLITSOURCE);
^SourceLocation{SPC}*={SPC}*               BEGIN(SRCLOC);
^UseEvenOdd{SPC}*={SPC}*                   BEGIN(EO);
^WritePropagatorFormat{SPC}*={SPC}*        BEGIN(WRPROPFLAG);
^PropagatorType{SPC}*={SPC}*               BEGIN(WRPROPFLAG);
^PerformOnlineMeasurements{SPC}*={SPC}*    BEGIN(ONMEAS);
^OnlineMeasurementsFreq{SPC}*={SPC}*       BEGIN(ONFREQ);

^BeginMonomial{SPC}+                       BEGIN(INITMONOMIAL);
^BeginInt                                  BEGIN(INITINTEGRATOR);

<INITMONOMIAL>{TYPE} {
  current_monomial++;
  mnl = &monomial_list[current_monomial];
  mnl->id = current_monomial;
  if(strcmp(yytext, "DET")==0) {
    mnl->type = DET;
    strcpy((*mnl).name, "DET");
  }
  else if(strcmp(yytext, "DETRATIO")==0) {
    mnl->type = DETRATIO;
    strcpy((*mnl).name, "DETRATIO");
  }
  else if(strcmp(yytext, "NDPOLY")==0) {
    mnl->type = NDPOLY;
    strcpy((*mnl).name, "NDPOLY");
    g_running_phmc = 1;
  }
  else if(strcmp(yytext, "GAUGE")==0) {
    mnl->type = GAUGE;
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

  /* here we need to set default values */
  BEGIN(MONOMIAL);
}
<MONOMIAL>{
  {SPC}*Timescale{SPC}*={SPC}*{DIGIT}+ {
    sscanf(yytext, " %[a-zA-Z] = %d", name, &a);
    mnl->timescale = a;
    if(myverbose) printf("  timescales set to %d line %d monomial %d\n", a, line_of_file, current_monomial);
  }
  {SPC}*2KappaMu{SPC}*={SPC}*{FLT} {
    sscanf(yytext, " %[2a-zA-Z] = %f", name, &c);
    mnl->mu = c;
    if(myverbose) printf("  2kappamu set to %f line %d monomial %d\n", c, line_of_file, current_monomial);
  }
  {SPC}*2KappaMu2{SPC}*={SPC}*{FLT} {
    sscanf(yytext, " %[2a-zA-Z] = %f", name, &c);
    mnl->mu2 = c;
    if(myverbose) printf("  2kappamu set to %f line %d monomial %d\n", c, line_of_file, current_monomial);
  }
  {SPC}*Kappa{SPC}*={SPC}*{FLT} {
    sscanf(yytext, " %[a-zA-Z] = %f", name, &c);
    mnl->kappa = c;
    if(myverbose) printf("  Kappa set to %f line %d monomial %d\n", c, line_of_file, current_monomial);
  }
  {SPC}*Kappa2{SPC}*={SPC}*{FLT} {
    sscanf(yytext, " %[2a-zA-Z] = %f", name, &c);
    mnl->kappa2 = c;
    if(myverbose) printf("  Kappa set to %f line %d monomial %d\n", c, line_of_file, current_monomial);
  }
  {SPC}*CSGHistory{SPC}*={SPC}*{DIGIT}+ {
    sscanf(yytext, " %[a-zA-Z] = %d", name, &a);
    mnl->csg_N = a;
    if(myverbose) printf("  csg history length set to %d line %d monomial %d\n", a, line_of_file, current_monomial);
  }
  {SPC}*CSGHistory2{SPC}*={SPC}*{DIGIT}+ {
    sscanf(yytext, " %[a-zA-Z2] = %d", name , &a);
    mnl->csg_N2 = a;
    if(myverbose) printf("  csg history2 length (for bicgstab) set to %d line %d monomial %d\n", 
                       a, line_of_file, current_monomial);
  }
  {SPC}*ForcePrecision{SPC}*={SPC}*{FLT} {
    sscanf(yytext, " %[a-zA-Z] = %f",name , &c);
    mnl->forceprec = c;
    if(myverbose) printf("  ForcePrecision set to %e line %d monomial %d\n", c, line_of_file, current_monomial);
  }
  {SPC}*AcceptancePrecision{SPC}*={SPC}*{FLT} {
    sscanf(yytext, " %[a-zA-Z] = %f",name , &c);
    mnl->accprec = c;
    if(myverbose) printf("  AcceptancePrecision set to %e line %d monomial %d\n", c, line_of_file, current_monomial);
  }
  {SPC}*MaxSolverIterations{SPC}*={SPC}*{DIGIT} {
    sscanf(yytext, " %[a-zA-Z] = %d", name, &a);
    mnl->maxiter = a;
    if(myverbose) printf("  MaxSolverIterations set to %d line %d monomial %d\n", a, line_of_file, current_monomial);
  }
  {SPC}*UseRectangleStaples{SPC}*={SPC}*yes {
    mnl->use_rectangles = 1;
    g_dbw2rand = 1;
    if(myverbose) printf("  UseRectangleStaples set to true line %d monomial %d\n", line_of_file, current_monomial);
  }
  {SPC}*UseRectangleStaples{SPC}*={SPC}*no {
    mnl->use_rectangles = 0;
    if(myverbose) printf("  UseRectangleStaples set to false line %d monomial %d\n", line_of_file, current_monomial);
  }
  {SPC}*Beta{SPC}*={SPC}*{FLT} {
    sscanf(yytext, " %[a-zA-Z] = %f",name , &c);
    mnl->beta = c;
    g_beta = c;
    if(myverbose) printf("  beta set to %e line %d monomial %d\n", c, line_of_file, current_monomial);
  }
  {SPC}*RectangleCoefficient{SPC}*={SPC}*{FLT} {
    sscanf(yytext, " %[a-zA-Z] = %f",name , &c);
    mnl->c1 = c;
    g_rgi_C1 = c;
    if(myverbose) printf("  RectangleCoefficient c1  set to %e line %d monomial %d\n", c, line_of_file, current_monomial);
  }
  ^EndMonomial{SPC}*  {
    if(myverbose) printf("monomial %d parsed line %d\n\n", current_monomial, line_of_file);
    BEGIN(0);
  }
  {SPC}*Name{SPC}*={SPC}* BEGIN(MNAME);
  {SPC}*Solver{SPC}*={SPC}* BEGIN(MSOLVER);
}
<MNAME>{NAME} {
  if(myverbose) printf("  monomial named \"%s\" line %d monomial %d\n", yytext, line_of_file, current_monomial);
  strcpy((*mnl).name, yytext);
  BEGIN(MONOMIAL);
}
<MSOLVER>{
  CG {
    if(myverbose) printf("  solver set to \"%s\" line %d monomial %d\n", yytext, line_of_file, current_monomial);
    mnl->solver = 1;
    BEGIN(MONOMIAL);
  }
  bicgstab {
    if(myverbose) printf("  solver set to \"%s\" line %d monomial %d\n", yytext, line_of_file, current_monomial);
    mnl->solver = 0;
    BEGIN(MONOMIAL);
  }
}

<INITINTEGRATOR>egrator{SPC}* {
  Integrator.no_timescales = -1;
  Integrator.tau = 1.;
  for(i = 0; i < 10; i++) {
    Integrator.lambda[i] = _default_2mn_lambda;
    Integrator.type[i] = MN2;
  }
  if(myverbose) printf("initialising integrator of type %s %d line %d\n", yytext, Integrator.type, line_of_file);
  BEGIN(INTEGRATOR);
}
<INTEGRATOR>{
  {SPC}*Type{DIGIT}{SPC}*={SPC}*{TYPE} {
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
    else {
      fprintf(stderr, "Unknown integrator type %s in line %d\n", yytext, line_of_file);
      exit(1);
    }

    if(myverbose) printf("  timescale %d type = %s line %d\n", a, type, line_of_file);
    free(type);
  }
  {SPC}*IntegrationSteps{DIGIT}{SPC}*={SPC}*{DIGIT}+ {
    sscanf(yytext, " %[a-zA-Z]%d = %d", name, &a, &b);
    if(myverbose) printf("  timescale %d steps=%d line %d\n", a, b, line_of_file);
    Integrator.n_int[a] = b;
    BEGIN(INTEGRATOR);
  }
  {SPC}*Lambda{DIGIT}{SPC}*={SPC}*{FLT} {
    sscanf(yytext, " %[a-zA-Z]%d = %f", name, &a, &c);
    Integrator.lambda[a] = c;
    if(myverbose) printf("  timescale %d Lambda=%f line %d\n", a, c, line_of_file);
    BEGIN(INTEGRATOR);
  }
  {SPC}*NumberOfTimescales{SPC}*={SPC}*{DIGIT}+ {
    sscanf(yytext, " %[a-zA-Z] = %d", name, &a);
    if(myverbose) printf("  Number of timescales set to %d line %d\n", a, line_of_file);
    if(a > 10) {
      if(g_proc_id == 0) fprintf(stderr, "maximal number of timescales is 10! Aborting...!\n");
      exit(-1);
    }
    Integrator.no_timescales = a;
    BEGIN(INTEGRATOR);
  }
  {SPC}*Tau{SPC}*={SPC}*{FLT} {
    sscanf(yytext, " %[a-zA-Z] = %f", name, &c);
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
<SEED>{DIGIT}+               {
  random_seed=atoi(yytext);
  if(myverbose!=0) printf("seed=%s \n", yytext);
}
<KAPPA>{FLT}  {
  g_kappa=atof(yytext);
  if(myverbose!=0) printf("kappa=%s \n", yytext);
}
<ACCPTILDE>{FLT}  {
  g_acc_Ptilde=atof(yytext);
  if(myverbose!=0) printf("Acc_Ptilde=%s \n", yytext);
}
<ACCHFIN>{FLT}  {
  g_acc_Hfin=atof(yytext);
  if(myverbose!=0) printf("Acc_Hfin=%s \n", yytext);
}
<RECEV>{DIGIT}+               {
  g_rec_ev = atoi(yytext);
  if(myverbose!=0) printf("Rec_EV=%s \n", yytext);
}
<MUBAR>{FLT}  {
  g_mubar=atof(yytext);
  if(myverbose!=0) printf("mubar=%s \n", yytext);
}
<EPSBAR>{FLT}  {
  g_epsbar=atof(yytext);
  if(myverbose!=0) printf("epsbar=%s \n", yytext);
}
<MU>{FLT}  {
  g_mu1=atof(yytext);
  if(myverbose!=0) printf("mu=%s \n", yytext);
}
<BETA>{FLT}   {
  g_beta=atof(yytext);
  if(myverbose!=0) printf("beta=%s \n",yytext);
}
<STARTCOND>cold {
  startoption=0; 
  if(myverbose!=0) printf("Start Condition is %s \n",yytext);
}
<STARTCOND>{
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
<NSKIP>{DIGIT}+ {
  Nskip=atoi(yytext);
  if(myverbose!=0) printf("Nskip= %s \n",yytext);
}
<SOLVFLAG>{
  bicgstab {
    solver_flag=0;
    if(myverbose!=0) printf("Use BiCGStab Solver");
  }
  cg {
    solver_flag=1;
    if(myverbose!=0) printf("Use CG Solver\n");
  }
  pcg {
    solver_flag=9;
    if(myverbose!=0) printf("Use PCG Solver (eigenvectors needed) \n");
  }
  gmres {
    solver_flag=2;
    if(myverbose!=0) printf("Use GMRES Solver\n");
  }
  gcr {
    solver_flag=7;
    if(myverbose!=0) printf("Use GCR Solver\n");
  }
  gmresdr {
    solver_flag=8;
    if(myverbose!=0) printf("Use GMRES-DR Solver\n");
  }
  cgs {
    solver_flag=3;
    if(myverbose!=0) printf("Use CGS Solver\n");
  }
  mr {
    solver_flag=4;
    if(myverbose!=0) printf("Use MR Solver \n");
  }
  bicgstabell {
    solver_flag=5;
    if(myverbose!=0) printf("Use BiCGstab(2) Solver \n");
  }
  fgmres {
    solver_flag=6;
    if(myverbose!=0) printf("Use FGMRES solver (eigenvectors needed) \n");
  }
}
<GMRESM>{DIGIT}+ {
  gmres_m_parameter = atoi(yytext);
  if(myverbose!=0) printf("Use Krylov Space of size %d in GMRES \n", gmres_m_parameter);
}
<GMRESDRNEV>{DIGIT}+ {
  gmresdr_nr_ev = atoi(yytext);
  if(myverbose!=0) printf("Deflate %d eigenvectors in GMRES-DR \n", gmresdr_nr_ev);
}
<SITER>{DIGIT}+ {
  max_solver_iterations = atoi(yytext);
  if(myverbose!=0) printf("Use %d iterations in the solvers!\n", max_solver_iterations);
}
<SPREC>{FLT} {
  solver_precision = atof(yytext);
  if(myverbose!=0) printf("Use %e as convergence precision for the solvers!\n", solver_precision);
}
<OPFLAG>Dtm {
  operator_flag=2;
  if(myverbose!=0) printf("Operator Flag is set to %s\n",yytext);
}
<OPFLAG>DWilson {
  operator_flag=1;
  if(myverbose!=0) printf("Operator Flag is set to %s\n",yytext);
}
<OPFLAG>DOverlap {
  operator_flag=0;
  if(myverbose!=0) printf("Operator Flag is set to %s\n",yytext);
}
<MEFLAG>yes {
  matrix_element_flag=1;
  if(myverbose!=0) printf("Compute Matrix Elements: %s\n", yytext);
}
<MEFLAG>no {
  matrix_element_flag=0;
  if(myverbose!=0) printf("Compute Matrix Elements: %s\n", yytext);
}
<SAVECONF>yes    {
  save_config_flag=1;
  if(myverbose!=0) printf("Save configurations\n");
}
<SAVECONF>no     {
  save_config_flag=0;
  if(myverbose!=0) printf("Don't save configurations\n");
}
<SAVEPROP>yes    {
  save_prop_flag=1;
  if(myverbose!=0) printf("Save propagators\n");
}
<SAVEPROP>no     {
  save_prop_flag=0;
  if(myverbose!=0) printf("Don't save propagators\n");
}
<SAVEPRG2>yes    {
  save_prop_g2_flag=1;
  if(myverbose!=0) printf("Save generalized propagators\n");
}
<SAVEPRG2>no     {
  save_prop_g2_flag=0;
  if(myverbose!=0) printf("Don't save generalized propagators\n");
}
<WRITECP>yes     {
  write_cp_flag=1;
  if(myverbose!=0) printf("Write Checkpoints\n");
}
<WRITECP>no     {
  write_cp_flag=0;
  if(myverbose!=0) printf("Don't write Checkpoints\n");
}
<CPINT>{DIGIT}+   {
  cp_interval=atoi(yytext);
  if(myverbose!=0) printf("Write Checkpoint all %s measurements\n",yytext);
}
<GAUGEINPUTFILE>{FILENAME} {
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
  if((index_start < 0)||(index_start >11)){
    printf("Error in line %d! index_start must be in [0,11]! Exiting...!\n", line_of_file);
    exit(1);
  }
  if(myverbose!=0) printf("inverting for index %s\n", yytext);
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
  if(myverbose!=0) printf("inverting up to color index %d\n", index_end);
  index_end+=1;
}
<FPROP>no {
  first_prop_flag = -1;
  if(myverbose!=0) printf("Do not compute the first propagator (default)\n");
}
<FPROP>compute {
  first_prop_flag = 0;
  if(myverbose!=0) printf("Computing the first propagator (default)\n");
}
<FPROP>readin {
  first_prop_flag = 1;
  if(myverbose!=0) printf("Reading in the first propagator\n");
}
<BCGMAX>{DIGIT}+ {
  ITER_MAX_BCG = atoi(yytext);
  if(myverbose != 0) printf("Maximal number of iterations for BCGstab set ro %d\n", ITER_MAX_BCG);
}
<CGMAX>{DIGIT}+  {
  ITER_MAX_CG = atoi(yytext);
  if(myverbose != 0) printf("Maximal number of iterations for CG set ro %d\n", ITER_MAX_CG);
}
<BOUND>{FLT} {
  X0 = atof(yytext);
  if(myverbose != 0) printf("X0 for boundary cond. in time set to %e\n", X0);
}
<READSOURCE>yes     {
  read_source_flag=1;
  if(myverbose!=0) printf("Read inversion source from file\n");
}
<READSOURCE>no     {
  read_source_flag=0;
  if(myverbose!=0) printf("Don't read inversion source from file\n");
}
<SOURCEFILE>{FILENAME} {
  strcpy(source_input_filename,yytext);
  if(myverbose!=0) printf("source input filename set to %s\n",yytext);
}
<SOURCEFORMAT>etmc      {
  source_format_flag = 0;
  if(myverbose!=0) printf("Using standard ETMC binary format for source input file\n");
}
<SOURCEFORMAT>cmi      {
  source_format_flag = 1;
  if(myverbose!=0) printf("Using CM format for source input file\n");
}
<SOURCEFORMAT>gwc      {
  source_format_flag = 2;
  if(myverbose!=0) printf("Using GWC format for source input file\n");
}
<SOURCETS>{DIGIT}+   {
  source_time_slice = atoi(yytext);
  if(myverbose!=0) printf("Using only timeslice %s of the source, padding the rest with zeros\n", yytext);
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
<PROPPREC>32 {
  prop_precision_flag = 32;
  if(myverbose!=0) printf("Save propagators in 32 Bit precision!\n");
}
<PROPPREC>64 {
  prop_precision_flag = 64;
  if(myverbose!=0) printf("Save propagators in 64 Bit precision!\n");
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
<STOUTRHO>{FLT} {
  stout_rho=atof(yytext);
  if(myverbose!=0) printf("use stout rho=%e!\n", stout_rho);
}
<STOUTITER>{DIGIT}+ {
  stout_no_iter=atoi(yytext);
  if(myverbose!=0) printf("make %d stout iterations!\n", stout_no_iter);
}
<PHMCFLAV>2"+"1"+"1 {
  phmc_no_flavours=4;
  if(myverbose!=0) printf("Simulate 2+1+1 flavours (1+1 PHMC).\n");
}
<PHMCFLAV>1"+"1 {
  phmc_no_flavours=2;
  if(myverbose!=0) printf("Simulate 1+1 flavours only (1+1 PHMC).\n");
}
<PCOMPUTEEVS>yes {
  phmc_compute_evs=1;
  if(myverbose!=0) printf("Compute Eigenvalues and exit.");
}
<PCOMPUTEEVS>no {
  phmc_compute_evs=0;
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
<PPP>no {
 phmc_exact_poly = 0;
 if(myverbose!=0) printf("Run the PHMC as usual.");
}
<PPP>yes {
 phmc_exact_poly = 1;
 if(myverbose!=0) printf("Run the PHMC only with usage of the less accurate polynomial.");
}

<SMAX>{FLT} {
  stilde_max = atof(yytext);
  if(myverbose!=0) printf("Stilde max for PHMC set to %e.\n", stilde_max);
}
<SMIN>{FLT} {
  stilde_min = atof(yytext);
  if(myverbose!=0) printf("Stilde min for PHMC set to %e.\n", stilde_min);
}
<DEGP>{DIGIT}+ {
  degree_of_p = atoi(yytext);
  if(myverbose!=0) printf("Degree for less precise polynomial P set to %d \n", degree_of_p);
}
<SPLITPROP>yes {
  propagator_splitted=1;
  if(myverbose!=0) printf("Split the propagator in several files! (invert)\n");
}
<SPLITPROP>no {
  propagator_splitted=0;
  if(myverbose!=0) printf("Do not split the propagator in several files (default) (invert)!\n");
}
<SPLITSOURCE>yes {
  source_splitted=1;
  if(myverbose!=0) printf("Expect source to be split in several files (invert)!\n");
}
<SPLITSOURCE>no {
  source_splitted=0;
  if(myverbose!=0) printf("Do not expect source to be split in several files (default) (invert)!\n");
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
<WRPROPFLAG>gwc      {
  write_prop_format_flag = 10;
  if(myverbose!=0) fprintf(stderr, "GWC format no longer supported for writing propagators\n");
}
<WRPROPFLAG>cmi      {
  write_prop_format_flag = 11;
  if(myverbose!=0) fprintf(stderr, "CM format no longer supported for writing propagators\n");
}
<WRPROPFLAG>DiracFermion_Sink {
  write_prop_format_flag = 0;
  if(myverbose!=0) printf("Propagator type: DiracFermion_Sinks\n");
}
<WRPROPFLAG>DiracFermion_Source_Sink_Pairs {
  write_prop_format_flag = 1;
  if(myverbose!=0) printf("Propagator type: DiracFermion_Source_Sink_Pairs\n");
}
<WRPROPFLAG>DiracFermion_ScalarSource_TwelveSink {
  write_prop_format_flag = 1;
  fprintf(stderr, "Propagator type: DiracFermion_ScalarSource_TwelveSink, not yet supported\n");
}
<WRPROPFLAG>DiracFermion_ScalarSource_FourSink {
  write_prop_format_flag = 1;
  fprintf(stderr, "Propagator type: DiracFermion_ScalarSource_FourSink, not yet supported\n");
}
<ONMEAS>yes {
  online_measurement_flag = 1;
  if(myverbose!=0) fprintf(stderr, "Switched on online measurements\n");
}
<ONMEAS>no {
  online_measurement_flag = 0;
  if(myverbose!=0) fprintf(stderr, "Online measurements not switched on\n");
}
<ONFREQ>{DIGIT}+ {
  online_measurement_freq = atoi(yytext);
  if(myverbose!=0) fprintf(stderr, "Frequency for online measurements set to %s\n", yytext);
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


<INITMONOMIAL,MONOMIAL,INTEGRATOR,INITINTEGRATOR>\n   {
  line_of_file++;
}
<*>\n                       {
  line_of_file++;
  BEGIN(0);
}

<*>. {
  BEGIN(ERROR);
}
<ERROR>[^\t\n]*             {
  printf("Parsing error in line %d\nAborting...!\n", line_of_file);
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
  g_kappa = _default_g_kappa;
  g_acc_Ptilde = _default_g_acc_Ptilde;
  g_acc_Hfin = _default_g_acc_Hfin;
  g_rec_ev = _default_g_rec_ev;
  g_mubar = _default_g_mubar;
  g_epsbar = _default_g_epsbar;
  g_mu = _default_g_mu;
  g_mu1 = _default_g_mu1;
  g_mu2 = _default_g_mu2;
  g_mu3 = _default_g_mu3;
  g_dbw2rand = 0;
  g_running_phmc = 0;
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
  g_relative_precision_flag = _default_g_relative_precision_flag;
  return_check_flag = _default_return_check_flag;
  return_check_interval = _default_return_check_interval;
  g_debug_level = _default_g_debug_level;
  source_format_flag = _default_source_format_flag;
  source_time_slice = _default_source_time_slice;
  gmres_m_parameter = _default_gmres_m_parameter;
  gmresdr_nr_ev = _default_gmresdr_nr_ev;
  gauge_precision_read_flag = _default_gauge_precision_read_flag;
  gauge_precision_write_flag = _default_gauge_precision_write_flag;
  prop_precision_flag = _default_prop_precision_flag;
  reproduce_randomnumber_flag = _default_reproduce_randomnumber_flag;
  g_sloppy_precision_flag = _default_g_sloppy_precision_flag;
  use_stout_flag = _default_use_stout_flag;
  stout_rho = _default_stout_rho;
  stout_no_iter = _default_stout_no_iter;
  /* check for reread ! */ 
  phmc_no_flavours = _default_phmc_no_flavours;
  phmc_compute_evs = _default_phmc_compute_evs;
  compute_evs = _default_compute_evs;
  stilde_min = _default_stilde_min;
  stilde_max = _default_stilde_max;
  degree_of_p = _default_degree_of_p;
  propagator_splitted = _default_propagator_splitted;
  source_splitted = _default_source_splitted;
  source_location = _default_source_location;
  eigenvalue_precision = _default_eigenvalue_precision;
  no_eigenvalues = _default_no_eigenvalues;
  sub_evs_cg_flag = _default_sub_evs_cg_flag;
  phmc_heavy_timescale = _default_phmc_heavy_timescale;
  phmc_exact_poly = _default_phmc_exact_poly;
  even_odd_flag = _default_even_odd_flag;
  online_measurement_flag = _default_online_measurement_flag;
  online_measurement_freq = _default_online_measurement_freq;

  /* Put -1 in write_prop_format_flag to see if parse_config() will
  change the value. If not then set it to source_format_flag */
  write_prop_format_flag = -1;
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
  
  if(write_prop_format_flag == -1) write_prop_format_flag = source_format_flag;
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
  double x;

  if(verbose && g_proc_id == 0) {
    myverbose = 1;
  }
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
  g_ka_csw_8 = g_kappa*g_c_sw/8.;
  
  g_mu = g_mu1;

  fclose(yyout);
  fclose(yyin);
  return(0);
}

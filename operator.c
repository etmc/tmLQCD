/***********************************************************************
 *
 * Copyright (C) 2009 Carsten Urbach
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

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <time.h>
#ifdef MPI
# include <mpi.h>
#endif
#include "global.h"
#include "default_input_values.h"
#include "read_input.h"
#include "su3.h"
#include "tm_operators.h"
#include "linalg_eo.h"
#include "D_psi.h"
#include "Dov_psi.h"
#include "Nondegenerate_Matrix.h"
#include "invert_eo.h"
#include "invert_doublet_eo.h"
#include "invert_overlap.h"
#include "observables.h"
#include "boundary.h"
#include "init_chi_spinor_field.h"
#include "start.h"
#include "solver/eigenvalues.h"
#include "solver/solver.h"
#include <io/params.h>
#include <io/gauge.h>
#include <io/spinor.h>
#include <io/utils.h>
#include "test/overlaptests.h"
#include "solver/index_jd.h"
#include "operator.h"

void dummy_D(spinor * const, spinor * const);
void dummy_DbD(spinor * const s, spinor * const r, spinor * const p, spinor * const q);
void op_invert(const int op_id, const int index_start);
void op_write_prop(const int op_id, const int index_start, const int append_);

operator operator_list[max_no_operators];

int no_operators = 0;

int add_operator(const int type) {

  operator * optr = &operator_list[no_operators];
  if(no_operators == max_no_operators) {
    fprintf(stderr, "maximal number of operators %d exceeded!\n", max_no_operators);
    exit(-1);
  }
  optr->type = type;
  optr->kappa = _default_g_kappa;
  optr->mu = _default_g_mu;
  optr->c_sw = _default_c_sw;
  optr->sloppy_precision = _default_g_sloppy_precision_flag;
  optr->coefs = NULL;
  optr->rel_prec = _default_g_relative_precision_flag;
  optr->eps_sq = _default_solver_precision;
  optr->maxiter = _default_max_solver_iterations;
  optr->even_odd_flag = _default_even_odd_flag;
  optr->solver = _default_solver_flag;
  optr->mubar = _default_g_mubar;
  optr->epsbar = _default_g_epsbar;
  optr->sr0 = NULL;
  optr->sr1 = NULL;
  optr->sr2 = NULL;
  optr->sr3 = NULL;
  optr->prop0 = NULL;
  optr->prop1 = NULL;
  optr->prop2 = NULL;
  optr->prop3 = NULL;
  optr->error_code = 0;
  optr->prop_precision = _default_prop_precision_flag;
  optr->no_flavours = 1;
  optr->DownProp = 0;
  optr->conf_input = _default_gauge_input_filename;

  optr->applyM = &dummy_D;
  optr->applyQ = &dummy_D;
  optr->applyQp = &dummy_D;
  optr->applyQm = &dummy_D;
  optr->applyMp = &dummy_D;
  optr->applyMm = &dummy_D;
  optr->applyQsq = &dummy_D;
  optr->applyDbQsq = &dummy_DbD;

  optr->inverter = &op_invert;
  optr->write_prop = &op_write_prop;

  /* Overlap needs special treatment */
  if(optr->type == OVERLAP) {
    optr->even_odd_flag = 0;
    optr->solver = 13;
    optr->no_ev = 10;
    optr->no_ev_index = 8;
    optr->ev_prec = 1.e-15;
    optr->ev_readwrite = 0;
    optr->deg_poly = 50;
    optr->s = 0.6;
    optr->m = 0.;
    optr->inverter = &op_invert;
  }
  if(optr->type == DBTMWILSON) {
    optr->no_flavours = 2;
    g_running_phmc = 1;
  }
  
  optr->precWS=NULL;

  optr->initialised = 1;

  no_operators++;
  return(no_operators);
}

int init_operators() {
  FILE * ifs;
  int i;
  operator * optr;
  for(i = 0; i < no_operators; i++) {
    optr = operator_list + i;
    /* This is a hack, it should be set on an operator basis. */
    optr->rel_prec = g_relative_precision_flag;
    if(optr->type == TMWILSON || optr->type == WILSON) {
      if(optr->even_odd_flag) {
        optr->applyQp = &Qtm_plus_psi;
        optr->applyQm = &Qtm_minus_psi;
        optr->applyQsq = &Qtm_pm_psi;
        optr->applyMp = &Mtm_plus_psi;
        optr->applyMm = &Mtm_minus_psi;
      }
      else {
        optr->applyQp = &Q_plus_psi;
        optr->applyQm = &Q_minus_psi;
        optr->applyQsq = &Q_pm_psi;
        optr->applyMp = &D_psi;
        optr->applyMm = &D_psi;
      }
      if(optr->solver == 12) {
        if (g_cart_id == 0 && optr->even_odd_flag == 1)
          fprintf(stderr, "CG Multiple mass solver works only without even/odd! Forcing!\n");
        optr->even_odd_flag = 0;
        if (g_cart_id == 0 && optr->DownProp)
          fprintf(stderr, "CGMMS doesn't need AddDownPropagator! Switching Off!\n");
        optr->DownProp = 0;

        /* this is for the extra masses of the CGMMS */
        if (g_no_extra_masses > 0) {
          if ((ifs = fopen("extra_masses.input", "r")) != NULL) {
            for (i = 0; i < g_no_extra_masses; i++) {
              /* Code added below mainly to stop the compiler from whining! */
              if (fscanf(ifs, "%lf", &g_extra_masses[i]) == EOF) {
                if (g_cart_id == 0) {
                  fprintf(stderr, "Expected %d extra masses for multiple mass solver, found only %d in extra_masses.input.\n", g_no_extra_masses, i);
                  fprintf(stderr, "Reducing the number of extra masses to %d due to lack of input values.\n", i);
                }
                g_no_extra_masses = i;
                break;
              }
              if (g_cart_id == 0 && g_debug_level > 0) {
                printf("# Extra mass %d = %lf\n", i, g_extra_masses[i]/(2 * optr->kappa));
              }
            }
            fclose(ifs);
          }
          else {
            fprintf(stderr, "Could not open file extra_masses.input!\n");
            g_no_extra_masses = 0;
          }
        }
      }
    }
    else if(optr->type == OVERLAP) {
      optr->even_odd_flag = 0;
      optr->applyM = &Dov_psi;
      optr->applyQ = &Qov_psi;
    }
    else if(optr->type == DBTMWILSON) {
      optr->even_odd_flag = 1;
      optr->applyDbQsq = &Q_Qdagger_ND;
      /* TODO: this should be here!       */
      /* Chi`s-spinors  memory allocation */
      /*       if(init_chi_spinor_field(VOLUMEPLUSRAND/2, 20) != 0) { */
      /* 	fprintf(stderr, "Not enough memory for 20 NDPHMC Chi fields! Aborting...\n"); */
      /* 	exit(0); */
      /*       } */
    }
  }
  return(0);
}

void dummy_D(spinor * const s, spinor * const r) {
  if(g_proc_id == 0) {
    fprintf(stderr, "dummy_D was called. Was that really intended?\n");
  } 
  return;
}

void dummy_DbD(spinor * const s, spinor * const r, spinor * const p, spinor * const q) {
  if(g_proc_id == 0) {
    fprintf(stderr, "dummy_DbD was called. Was that really intended?\n");
  }
  return;
}

void op_invert(const int op_id, const int index_start) {
  operator * optr = &operator_list[op_id];
  double atime = 0., etime = 0., nrm1 = 0., nrm2 = 0.;
  int i;
  optr->iterations = 0;
  optr->reached_prec = -1.;
  g_kappa = optr->kappa;
  boundary(g_kappa);

#ifdef MPI
  atime = MPI_Wtime();
#else
  atime = (double)clock() / (double)(CLOCKS_PER_SEC);
#endif
  if(optr->type == TMWILSON || optr->type == WILSON) {
    for(i = 0; i < 2; i++) {
      g_mu = optr->mu;
      if (g_cart_id == 0) {
        printf("#\n# mu = %e\n", g_mu);
      }

      if(use_preconditioning==1){
        g_precWS=(void*)optr->precWS;
      }
      else {
        g_precWS=NULL;
      }
      // This is for the CGMMS: it header here because binary data will be wrote inside cg_mms
      //if(optr->solver == 12) {
      //  write_cgmms_headers(op_id,index_start,i);
      //}

      optr->iterations = invert_eo( optr->prop0, optr->prop1, optr->sr0, optr->sr1,
                                    optr->eps_sq, optr->maxiter,
                                    optr->solver, optr->rel_prec,
                                    0, optr->even_odd_flag);

      /* check result */
      M_full(g_spinor_field[4], g_spinor_field[5], optr->prop0, optr->prop1);

      diff(g_spinor_field[4], g_spinor_field[4], optr->sr0, VOLUME / 2);
      diff(g_spinor_field[5], g_spinor_field[5], optr->sr1, VOLUME / 2);

      nrm1 = square_norm(g_spinor_field[4], VOLUME / 2, 1);
      nrm2 = square_norm(g_spinor_field[5], VOLUME / 2, 1);
      optr->reached_prec = nrm1 + nrm2;

      /* convert to standard normalisation  */
      /* we have to mult. by 2*kappa        */
      if (optr->kappa != 0.) {
        mul_r(optr->prop0, (2*optr->kappa), optr->prop0, VOLUME / 2);
        mul_r(optr->prop1, (2*optr->kappa), optr->prop1, VOLUME / 2);
      }
      if (optr->solver != CGMMS) /* CGMMS handles its own I/O */
        optr->write_prop(op_id, index_start, i);
      if(optr->DownProp) {
        optr->mu = -optr->mu;
      } else 
        break;
    }
  }
  else if(optr->type == DBTMWILSON) {
    g_mubar = optr->mubar;
    g_epsbar = optr->epsbar;
    for(i = 0; i < SourceInfo.no_flavours; i++) {
      optr->iterations = invert_doublet_eo( optr->prop0, optr->prop1, optr->prop2, optr->prop3, 
                                            optr->sr0, optr->sr1, optr->sr2, optr->sr3,
                                            optr->eps_sq, optr->maxiter,
                                            optr->solver, optr->rel_prec);

      g_mu = optr->mubar;
      M_full(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI+2], optr->prop0, optr->prop1); 
      assign_add_mul_r(g_spinor_field[DUM_DERI+1], optr->prop2, -optr->epsbar, VOLUME/2);
      assign_add_mul_r(g_spinor_field[DUM_DERI+2], optr->prop3, -optr->epsbar, VOLUME/2);

      g_mu = -g_mu;
      M_full(g_spinor_field[DUM_DERI+3], g_spinor_field[DUM_DERI+4], optr->prop2, optr->prop3); 
      assign_add_mul_r(g_spinor_field[DUM_DERI+3], optr->prop0, -optr->epsbar, VOLUME/2);
      assign_add_mul_r(g_spinor_field[DUM_DERI+4], optr->prop1, -optr->epsbar, VOLUME/2);

      diff(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI+1], optr->sr0, VOLUME/2); 
      diff(g_spinor_field[DUM_DERI+2], g_spinor_field[DUM_DERI+2], optr->sr1, VOLUME/2); 
      diff(g_spinor_field[DUM_DERI+3], g_spinor_field[DUM_DERI+3], optr->sr2, VOLUME/2); 
      diff(g_spinor_field[DUM_DERI+4], g_spinor_field[DUM_DERI+4], optr->sr3, VOLUME/2); 

      nrm1  = square_norm(g_spinor_field[DUM_DERI+1], VOLUME/2, 1); 
      nrm1 += square_norm(g_spinor_field[DUM_DERI+2], VOLUME/2, 1); 
      nrm1 += square_norm(g_spinor_field[DUM_DERI+3], VOLUME/2, 1); 
      nrm1 += square_norm(g_spinor_field[DUM_DERI+4], VOLUME/2, 1); 
      optr->reached_prec = nrm1;
      g_mu = g_mu1;
      /* For standard normalisation */
      /* we have to mult. by 2*kappa */
      mul_r(g_spinor_field[DUM_DERI], (2*optr->kappa), optr->prop0, VOLUME/2);
      mul_r(g_spinor_field[DUM_DERI+1], (2*optr->kappa), optr->prop1, VOLUME/2);
      mul_r(g_spinor_field[DUM_DERI+2], (2*optr->kappa), optr->prop2, VOLUME/2);
      mul_r(g_spinor_field[DUM_DERI+3], (2*optr->kappa), optr->prop3, VOLUME/2);
      /* the final result should be stored in the convention used in */
      /* hep-lat/0606011                                             */
      /* this requires multiplication of source with                 */
      /* (1+itau_2)/sqrt(2) and the result with (1-itau_2)/sqrt(2)   */

      mul_one_pm_itau2(optr->prop0, optr->prop2, g_spinor_field[DUM_DERI], 
                       g_spinor_field[DUM_DERI+2], -1., VOLUME/2);
      mul_one_pm_itau2(optr->prop1, optr->prop3, g_spinor_field[DUM_DERI+1], 
                       g_spinor_field[DUM_DERI+3], -1., VOLUME/2);
      /* write propagator */
      optr->write_prop(op_id, index_start, i);

      mul_r(optr->prop0, 1./(2*optr->kappa), g_spinor_field[DUM_DERI], VOLUME/2);
      mul_r(optr->prop1, 1./(2*optr->kappa), g_spinor_field[DUM_DERI+1], VOLUME/2);
      mul_r(optr->prop2, 1./(2*optr->kappa), g_spinor_field[DUM_DERI+2], VOLUME/2);
      mul_r(optr->prop3, 1./(2*optr->kappa), g_spinor_field[DUM_DERI+3], VOLUME/2);

      /* mirror source, but not for volume sources */
      if(i == 0 && SourceInfo.no_flavours == 2 && SourceInfo.type != 1) {
        if (g_cart_id == 0) {
          fprintf(stdout, "# Inversion done in %d iterations, squared residue = %e!\n",
                  optr->iterations, optr->reached_prec);
        }
        mul_one_pm_itau2(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+2], optr->sr0, optr->sr2, -1., VOLUME/2);
        mul_one_pm_itau2(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI+3], optr->sr1, optr->sr3, -1., VOLUME/2);

        mul_one_pm_itau2(optr->sr0, optr->sr2, g_spinor_field[DUM_DERI+2], g_spinor_field[DUM_DERI], +1., VOLUME/2);
        mul_one_pm_itau2(optr->sr1, optr->sr3, g_spinor_field[DUM_DERI+3], g_spinor_field[DUM_DERI+1], +1., VOLUME/2);

      }
      /* volume sources need only one inversion */
      else if(SourceInfo.type == 1) i++;
    }
  }
  else if(optr->type == OVERLAP) {
    g_mu = 0.;
    m_ov=optr->m;
    eigenvalues(&optr->no_ev, 5000, optr->ev_prec, 0, optr->ev_readwrite, nstore, optr->even_odd_flag);
/*     ov_check_locality(); */
/*      index_jd(&optr->no_ev_index, 5000, 1.e-12, optr->conf_input, nstore, 4); */
    ov_n_cheby=optr->deg_poly;

    if(use_preconditioning==1)
      g_precWS=(void*)optr->precWS;
    else
      g_precWS=NULL;


    if(g_debug_level > 3) ov_check_ginsparg_wilson_relation_strong(); 

    invert_overlap(op_id, index_start); 

    optr->write_prop(op_id, index_start, 0);
  }
#ifdef MPI
  etime = MPI_Wtime();
#else
  etime = (double)clock() / (double)(CLOCKS_PER_SEC);
#endif
  if (g_cart_id == 0 && g_debug_level > 0) {
    fprintf(stdout, "# Inversion done in %d iterations, squared residue = %e!\n",
            optr->iterations, optr->reached_prec);
    fprintf(stdout, "# Inversion done in %1.2e sec. \n", etime - atime);
  }
  return;
}


void op_write_prop(const int op_id, const int index_start, const int append_) {
  operator * optr = &operator_list[op_id];
  char filename[100];
  char ending[15];
  WRITER *writer = NULL;
  int append = 0;
  int status = 0;

  paramsSourceFormat *sourceFormat = NULL;
  paramsPropagatorFormat *propagatorFormat = NULL;
  paramsInverterInfo *inverterInfo = NULL;
  if(optr->type == DBTMWILSON) {
    strcpy(ending, "hinverted");
  }
  else if(optr->type == OVERLAP) {
    strcpy(ending, "ovinverted");
  }
  else {
    strcpy(ending, "inverted");
  }

  if(SourceInfo.type != 1) {
    if (PropInfo.splitted) {
      sprintf(filename, "%s.%.4d.%.2d.%.2d.%s", SourceInfo.basename, SourceInfo.nstore, SourceInfo.t, SourceInfo.ix, ending);
    }
    else {
      sprintf(filename, "%s.%.4d.%.2d.%s", SourceInfo.basename, SourceInfo.nstore, SourceInfo.t, ending);
    }
  }
  else {
    sprintf(filename, "%s.%.4d.%.5d.%s", SourceInfo.basename, SourceInfo.nstore, SourceInfo.sample, ending);
  }

  if(!PropInfo.splitted || append_)
    append = 1;
  /* the 1 is for appending */
  construct_writer(&writer, filename, append);
  if (PropInfo.splitted || SourceInfo.ix == index_start) {
    inverterInfo = construct_paramsInverterInfo(optr->reached_prec, optr->iterations, 
						optr->solver, optr->no_flavours);
    write_spinor_info(writer, PropInfo.format, inverterInfo, append);
    free(inverterInfo);
  }
  /* write the source depending on format */
  /* to be fixed for 2 fl tmwilson        */
  if (PropInfo.format == 1) {
    sourceFormat = construct_paramsSourceFormat(SourceInfo.precision, optr->no_flavours, 4, 3);
    write_source_format(writer, sourceFormat);
    status = write_spinor(writer, &operator_list[op_id].sr0, &operator_list[op_id].sr1, 
			  1, SourceInfo.precision);
    if(optr->no_flavours == 2) {
      status = write_spinor(writer, &operator_list[op_id].sr2, &operator_list[op_id].sr3, 
			    1, SourceInfo.precision);
    }
    free(sourceFormat);
  }
  propagatorFormat = construct_paramsPropagatorFormat(optr->prop_precision, optr->no_flavours);
  write_propagator_format(writer, propagatorFormat);
  free(propagatorFormat);

  if(optr->no_flavours == 2) {
    status = write_spinor(writer, &operator_list[op_id].prop2, &operator_list[op_id].prop3, 1, optr->prop_precision);
  }
  status = write_spinor(writer, &operator_list[op_id].prop0, &operator_list[op_id].prop1, 1, optr->prop_precision);
  destruct_writer(writer);
  return;
}

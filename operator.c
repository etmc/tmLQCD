/***********************************************************************
 * $Id$ 
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
#include <math.h>
#include <errno.h>
#include <time.h>
#ifdef MPI
# include <mpi.h>
#endif
#include "global.h"
#include "default_input_values.h"
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
#include <io/params.h>
#include <io/gauge.h>
#include <io/spinor.h>
#include <io/utils.h>
#include "operator.h"

void dummy_D(spinor * const, spinor * const);
void dummy_DbD(spinor * const s, spinor * const r, spinor * const p, spinor * const q);
void op_invert(const int op_id);
void op_write_prop(const int nstore, const int isample, const int ix, const int op_id, 
		   const int source_time_slice, const int propagator_splitted, 
		   const int index_start, const int write_prop_format_flag,
		   char * source_input_filename, char *gaugelfn, DML_Checksum *gaugecksum);

operator operator_list[max_no_operators];

int no_operators = 0;

int add_operator(const int type) {

  if(no_operators == max_no_operators) {
    fprintf(stderr, "maximal number of operators %d exceeded!\n", max_no_operators);
    exit(-1);
  }
  operator_list[no_operators].type = type;
  operator_list[no_operators].kappa = _default_g_kappa;
  operator_list[no_operators].mu = _default_g_mu;
  operator_list[no_operators].sloppy_precision = _default_g_sloppy_precision_flag;
  operator_list[no_operators].coefs = NULL;
  operator_list[no_operators].rel_prec = _default_g_relative_precision_flag;
  operator_list[no_operators].eps_sq = _default_solver_precision;
  operator_list[no_operators].maxiter = _default_max_solver_iterations;
  operator_list[no_operators].even_odd_flag = _default_even_odd_flag;
  operator_list[no_operators].solver = _default_solver_flag;
  operator_list[no_operators].mubar = _default_g_mubar;
  operator_list[no_operators].epsbar = _default_g_epsbar;
  operator_list[no_operators].sr0 = NULL;
  operator_list[no_operators].sr1 = NULL;
  operator_list[no_operators].sr2 = NULL;
  operator_list[no_operators].sr3 = NULL;
  operator_list[no_operators].prop0 = NULL;
  operator_list[no_operators].prop1 = NULL;
  operator_list[no_operators].prop2 = NULL;
  operator_list[no_operators].prop3 = NULL;
  operator_list[no_operators].error_code = 0;
  operator_list[no_operators].prop_precision = _default_prop_precision_flag;

  operator_list[no_operators].applyM = &dummy_D;
  operator_list[no_operators].applyQ = &dummy_D;
  operator_list[no_operators].applyQp = &dummy_D;
  operator_list[no_operators].applyQm = &dummy_D;
  operator_list[no_operators].applyMp = &dummy_D;
  operator_list[no_operators].applyMm = &dummy_D;
  operator_list[no_operators].applyQsq = &dummy_D;
  operator_list[no_operators].applyDbQsq = &dummy_DbD;

  operator_list[no_operators].inverter = &op_invert;
  operator_list[no_operators].write_prop = &op_write_prop;

  /* Overlap needs special treatment */
  if(operator_list[no_operators].type == OVERLAP) {
    operator_list[no_operators].even_odd_flag = 0;
    operator_list[no_operators].solver = 13;
    operator_list[no_operators].no_ev = 10;
    operator_list[no_operators].ev_prec = 1.e-15;
    operator_list[no_operators].deg_poly = 50;
    operator_list[no_operators].s = 0.6;
    operator_list[no_operators].m = 0.;
    operator_list[no_operators].inverter = &invert_overlap;
  }


  operator_list[no_operators].initialised = 1;

  no_operators++;
  return(no_operators);
}

int init_operators() {
  int i;
  operator * optr;
  for(i = 0; i < no_operators; i++) {
    optr = operator_list + i;
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
    }
    else if(optr->type == OVERLAP) {
      optr->even_odd_flag = 0;
      optr->applyM = &Dov_psi;
      optr->applyQ = &Qov_psi;
    }
    else {
      optr->even_odd_flag = 1;
      optr->applyDbQsq = &Q_Qdagger_ND;
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

void op_invert(const int op_id) {
  operator * optr = &operator_list[op_id];
  double atime = 0., etime = 0., nrm1 = 0., nrm2 = 0.;
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
    g_mu = optr->mu;
    if (g_cart_id == 0) {
      printf("mu = %e\n", g_mu);
    }
    
    optr->iterations = invert_eo(optr->prop0, optr->prop1, optr->sr0, optr->sr1,
				 optr->eps_sq, optr->maxiter,
				 optr->solver, optr->rel_prec,
				 0, optr->even_odd_flag);
    
    /* check result */
    M_full(g_spinor_field[4], g_spinor_field[5], optr->prop0, optr->prop1);
    
    if (optr->kappa != 0.) {
      mul_r(g_spinor_field[4], 1. / (2*optr->kappa), g_spinor_field[4], VOLUME / 2);
      mul_r(g_spinor_field[5], 1. / (2*optr->kappa), g_spinor_field[5], VOLUME / 2);
    }
    
    diff(g_spinor_field[4], g_spinor_field[4], optr->sr0, VOLUME / 2);
    diff(g_spinor_field[5], g_spinor_field[5], optr->sr1, VOLUME / 2);
    
    nrm1 = square_norm(g_spinor_field[4], VOLUME / 2, 1);
    nrm2 = square_norm(g_spinor_field[5], VOLUME / 2, 1);
    optr->reached_prec = nrm1 + nrm2;

    if (g_cart_id == 0) {
      fprintf(stdout, "Inversion done in %d iterations, squared residue = %e!\n",
	      optr->iterations, optr->reached_prec);
    }

  }
  else if(optr->type == DBTMWILSON) {
    g_mubar = optr->mubar;
    g_epsbar = optr->epsbar;
    optr->iterations = invert_doublet_eo(optr->prop0, optr->prop1, optr->prop2, optr->prop3, 
					 optr->sr0, optr->sr1, optr->sr2, optr->sr3,
					 optr->eps_sq, optr->maxiter,
					 optr->solver, optr->rel_prec);
  }
  else if(optr->type == OVERLAP) {
    invert_overlap(op_id);
  }
#ifdef MPI
  etime = MPI_Wtime();
#else
  etime = (double)clock() / (double)(CLOCKS_PER_SEC);
#endif
  if (g_cart_id == 0) {
    printf("Inversion done in %1.2e sec. \n", etime - atime);
  }
  return;
}


void op_write_prop(const int nstore, const int isample, const int ix, const int op_id, 
		   const int source_time_slice, const int propagator_splitted, 
		   const int index_start, const int write_prop_format_flag,
		   char * source_input_filename, char *gaugelfn, DML_Checksum *gaugecksum) {
  operator * optr = &operator_list[op_id];
  double ratime = 0., retime = 0., plaquette_energy;
  char conf_filename[100];

  WRITER *writer = NULL;

  paramsXlfInfo *xlfInfo = NULL;
  paramsSourceFormat *sourceFormat = NULL;
  paramsPropagatorFormat *propagatorFormat = NULL;
  paramsInverterInfo *inverterInfo = NULL;
  
  
  if (propagator_splitted) {
    sprintf(conf_filename, "%s.%.4d.%.2d.%.2d.inverted", source_input_filename, nstore, source_time_slice, ix);
  }
  else {
    sprintf(conf_filename, "%s.%.4d.%.2d.inverted", source_input_filename, nstore, source_time_slice);
  }
  
  construct_writer(&writer, conf_filename);
  
  if(optr->type == TMWILSON || optr->type == WILSON) {    
    if (propagator_splitted || ix == index_start) {
      plaquette_energy = measure_gauge_action();
      xlfInfo = construct_paramsXlfInfo(plaquette_energy / (6.*VOLUME*g_nproc), nstore);
      inverterInfo = construct_paramsInverterInfo(optr->reached_prec, optr->iterations, optr->solver, 1);
      
      write_spinor_info(writer, xlfInfo, write_prop_format_flag, inverterInfo, gaugelfn, gaugecksum);
      
      free(xlfInfo);
      free(inverterInfo);
    }
    
    /* write the source depending on format */
    if (write_prop_format_flag == 1) {
      sourceFormat = construct_paramsSourceFormat(32, 1, 4, 3);
      
      write_source_format(writer, sourceFormat);
      write_spinor(writer, &operator_list[op_id].sr0, &operator_list[op_id].sr1, 1, 32);
      
      free(sourceFormat);
    }
    
#ifdef MPI
    ratime = MPI_Wtime();
#else
    ratime = (double)clock() / (double)(CLOCKS_PER_SEC);
#endif
    propagatorFormat = construct_paramsPropagatorFormat(optr->prop_precision, 1);
    
    write_propagator_format(writer, propagatorFormat);
    write_spinor(writer, &operator_list[op_id].prop0, &operator_list[op_id].prop1, 1, optr->prop_precision);
    
    free(propagatorFormat);
    
#ifdef MPI
    retime = MPI_Wtime();
#else
    retime = (double)clock() / (double)(CLOCKS_PER_SEC);
#endif
    if (g_cart_id == 0) {
      printf("time for writing prop was %e seconds\n", retime - ratime);
    }
  }
  else {

  }
  destruct_writer(writer);
  
  return;
}

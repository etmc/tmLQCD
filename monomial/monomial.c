/***********************************************************************
 *
 * Copyright (C) 2008,2011,2012 Carsten Urbach
 *               2009 Jenifer Gonzalez Lopez
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
#include <string.h>
#include "global.h"
#include "boundary.h"
#include "su3.h"
#include "su3adj.h"
#include "su3spinor.h"
#include "operator/tm_operators.h"
#include "operator/tm_operators_32.h"
#include "operator/clovertm_operators.h"
#include "operator/clovertm_operators_32.h"
#include "operator/clover_leaf.h"
#include "ranlxd.h"
#include "sse.h"
#include "linalg_eo.h"
#include "default_input_values.h"
#include "read_input.h"
#include "monomial/monomial.h"

monomial monomial_list[max_no_monomials];
int no_monomials = 0;
int no_gauge_monomials = 0;
int clover_monomials[max_no_monomials];
int clovernd_monomials[max_no_monomials];
int no_clover_monomials = 0;
int no_clovernd_monomials = 0;
static spinor * _pf;
spinor ** w_fields;
const int no_wfields = 6;

int add_monomial(const int type) {
  
  if(no_monomials == max_no_monomials) {
    fprintf(stderr, "maximal number of monomials %d exceeded!\n", max_no_monomials);
    exit(-1);
  }
  monomial_list[no_monomials].hbfunction = &dummy_heatbath;
  monomial_list[no_monomials].accfunction = &dummy_acc;
  monomial_list[no_monomials].derivativefunction = &dummy_derivative;

  monomial_list[no_monomials].pf = NULL;
  monomial_list[no_monomials].pf2 = NULL;
  monomial_list[no_monomials].w_fields = NULL;
  monomial_list[no_monomials].csg_field = NULL;
  monomial_list[no_monomials].csg_field2 = NULL;
  monomial_list[no_monomials].csg_index_array = NULL;
  monomial_list[no_monomials].csg_index_array2 = NULL;
  monomial_list[no_monomials].no_wfields = no_wfields;
  monomial_list[no_monomials].csg_N = 0;
  monomial_list[no_monomials].csg_N2 = 0;
  monomial_list[no_monomials].csg_n = 1;
  monomial_list[no_monomials].csg_n2 = 1;
  monomial_list[no_monomials].kappa = _default_g_kappa;
  monomial_list[no_monomials].kappa2 = _default_g_kappa;
  monomial_list[no_monomials].mu = _default_g_mu;
  monomial_list[no_monomials].mu2 = _default_g_mu;
  monomial_list[no_monomials].c_sw = _default_c_sw;
  monomial_list[no_monomials].rho = _default_rho;
  monomial_list[no_monomials].rho2 = _default_rho2;
  monomial_list[no_monomials].mubar = _default_g_mubar;
  monomial_list[no_monomials].mubar2 = _default_g_mubar;
  monomial_list[no_monomials].epsbar = _default_g_epsbar;
  monomial_list[no_monomials].epsbar2 = _default_g_epsbar;
  monomial_list[no_monomials].epsilon = _default_g_epsbar;
  monomial_list[no_monomials].timescale = _default_timescale;
  monomial_list[no_monomials].accprec = _default_g_eps_sq_acc;
  monomial_list[no_monomials].forceprec = _default_g_eps_sq_force;
  monomial_list[no_monomials].maxiter = _default_max_solver_iterations;
  if((monomial_list[no_monomials].type == NDRAT) ||
     (monomial_list[no_monomials].type == NDRATCOR) ||
     (monomial_list[no_monomials].type == NDCLOVERRAT) ||
     (monomial_list[no_monomials].type == NDCLOVERRATCOR)
  ) {
    monomial_list[no_monomials].solver = _default_nd_solver_flag;    
  }
  else{
    monomial_list[no_monomials].solver = _default_solver_flag;
  }
  monomial_list[no_monomials].solver_params.mcg_delta = _default_mixcg_innereps;
  monomial_list[no_monomials].solver_params.solution_type = TM_SOLUTION_M_MDAG;
  // the defaut is 1 because the QPhiX interface is generalised in such a way
  // that normal solves correspond to solves with one shift, this does not 
  // affect the used parameters in any way!
  monomial_list[no_monomials].solver_params.no_shifts = 1;
  monomial_list[no_monomials].solver_params.compression_type = _default_compression_type;
  monomial_list[no_monomials].solver_params.external_inverter = _default_external_inverter;
  monomial_list[no_monomials].solver_params.sloppy_precision = _default_operator_sloppy_precision_flag;
  monomial_list[no_monomials].even_odd_flag = _default_even_odd_flag;
  monomial_list[no_monomials].forcefactor = 1.;
  monomial_list[no_monomials].use_rectangles = 0;
  monomial_list[no_monomials].c1 = _default_g_rgi_C1;
  monomial_list[no_monomials].c0 = 1.;
  monomial_list[no_monomials].beta = _default_g_beta;
  monomial_list[no_monomials].glambda = 0.;
  monomial_list[no_monomials].rngrepro = _default_reproduce_randomnumber_flag;
  monomial_list[no_monomials].trlog = 0;
  /* poly monomial */
  monomial_list[no_monomials].rec_ev = _default_g_rec_ev;
  monomial_list[no_monomials].MDPolyDegree = _default_MDPolyDegree;
  monomial_list[no_monomials].MDPolyLmin = _default_MDPolyLmin;
  monomial_list[no_monomials].MDPolyLmax = _default_MDPolyLmax;
  strcpy(monomial_list[no_monomials].MDPolyRootsFile,_default_MDPolyRootsFile);
  monomial_list[no_monomials].MDPolyRoots = NULL;
  monomial_list[no_monomials].MDPoly_chi_spinor_fields = (spinor**)NULL;
  monomial_list[no_monomials].MDPolyLocNormConst = _default_MDPolyLocNormConst;
  monomial_list[no_monomials].MDPolyDetRatio = _default_MDPolyDetRatio;
  monomial_list[no_monomials].MaxPtildeDegree = NTILDE_CHEBYMAX;
  monomial_list[no_monomials].StildeMin = _default_stilde_min;
  monomial_list[no_monomials].StildeMax = _default_stilde_max;
  monomial_list[no_monomials].PrecisionHfinal = _default_g_acc_Hfin;
  monomial_list[no_monomials].PrecisionPtilde = _default_g_acc_Ptilde;

  monomial_list[no_monomials].rat.order = 12;
  monomial_list[no_monomials].rat.range[0] = _default_stilde_min;
  monomial_list[no_monomials].rat.range[1] = _default_stilde_max;
  monomial_list[no_monomials].rat.crange[0] = 0;
  monomial_list[no_monomials].rat.crange[1] = 11;

  monomial_list[no_monomials].initialised = 1;
  if(monomial_list[no_monomials].type == NDDETRATIO || monomial_list[no_monomials].type == NDCLOVERDETRATIO || monomial_list[no_monomials].type == CLOVERDETRATIORW) {
    monomial_list[no_monomials].timescale = -5;
  }

  no_monomials++;
  return(no_monomials);
}


int init_monomials(const int V, const int even_odd_flag) {
  int no=0;
  int retval;
  spinor * __pf = NULL;
  double sw_mu=0., sw_k=0., sw_c=0.;
  double swn_mubar=0., swn_epsbar = 0., swn_k=0., swn_c=0.;

  if (g_exposu3_no_c == 0) init_exposu3();
  
  for(int i = 0; i < no_monomials; i++) {
    if((monomial_list[i].type != GAUGE) && (monomial_list[i].type != SFGAUGE)) no++;
    /* non-degenerate monomials need two pseudo fermion fields */
    if((monomial_list[i].type == NDPOLY) || (monomial_list[i].type == NDDETRATIO) || (monomial_list[i].type == NDCLOVERDETRATIO) || 
       (monomial_list[i].type == NDCLOVER) || (monomial_list[i].type == NDRAT)||
       (monomial_list[i].type == NDRATCOR) || (monomial_list[i].type == NDCLOVERRATCOR) ||
       (monomial_list[i].type == NDCLOVERRAT)) no++;
  }
  if(no_monomials > 0) {
    if((void*)(_pf = (spinor*)calloc((no+no_wfields)*V+1, sizeof(spinor))) == NULL) {
      printf ("malloc errno in monomial pf fields: %d\n",errno); 
      errno = 0;
      return(1);
    }
    else {
      __pf = (spinor*)(((unsigned long int)(_pf)+ALIGN_BASE)&~ALIGN_BASE);
    }
    if((void*)(w_fields = (spinor**)calloc(no_wfields, sizeof(spinor*))) == NULL) {
      printf ("malloc errno in monomial  w_fields: %d\n",errno); 
      errno = 0;
      return(1);
    }
    for(int i = 0; i < no_wfields; i++) {
      w_fields[i] = __pf+(no+i)*V;
    }
  }
  
  no = 0;
  for(int i = 0; i < no_monomials; i++) {
    monomial_list[i].rngrepro = reproduce_randomnumber_flag;
    if((monomial_list[i].type != GAUGE) && (monomial_list[i].type != SFGAUGE)) {
      monomial_list[i].w_fields = w_fields;
      monomial_list[i].pf = __pf+no*V;
      no++;
      
      if(monomial_list[i].type == DET) {
	monomial_list[i].hbfunction = &det_heatbath;
	monomial_list[i].accfunction = &det_acc;
	monomial_list[i].derivativefunction = &det_derivative;
	if(even_odd_flag) {
	  monomial_list[i].Qsq = &Qtm_pm_psi;
	  monomial_list[i].Qp = &Qtm_plus_psi;
	  monomial_list[i].Qm = &Qtm_minus_psi;
	}
	else {
	  monomial_list[i].Qsq = &Q_pm_psi;
	  monomial_list[i].Qp = &Q_plus_psi;
	  monomial_list[i].Qm = &Q_minus_psi;
	}
	if(g_proc_id == 0 && g_debug_level > 1) {
	  printf("# Initialised monomial of type DET, no_monomials= %d\n", no_monomials);
	}
      }
      else if(monomial_list[i].type == CLOVERDET) {
	monomial_list[i].hbfunction = &cloverdet_heatbath;
	monomial_list[i].accfunction = &cloverdet_acc;
	monomial_list[i].derivativefunction = &cloverdet_derivative;
	//monomial_list[i].derivativefunction = &det_derivative;
	if(even_odd_flag) {
	  monomial_list[i].Qsq = &Qsw_pm_psi;
	  monomial_list[i].Qp = &Qsw_plus_psi;
	  monomial_list[i].Qm = &Qsw_minus_psi;
	}
	else {
	  monomial_list[i].Qsq = &Qsw_full_pm_psi;
	  monomial_list[i].Qp = &Qsw_full_plus_psi;
	  monomial_list[i].Qm = &Qsw_full_minus_psi;
	}
	init_swpm(VOLUME);
	clover_monomials[no_clover_monomials] = i;
	no_clover_monomials++;
	if(g_proc_id == 0 && g_debug_level > 1) {
	  printf("# Initialised monomial of type CLOVERDET, no_monomials= %d\n", no_monomials);
	}
      }
      else if(monomial_list[i].type == CLOVERDETRATIO) {
	monomial_list[i].hbfunction = &cloverdetratio_heatbath;
	monomial_list[i].accfunction = &cloverdetratio_acc;
	monomial_list[i].derivativefunction = &cloverdetratio_derivative;
	monomial_list[i].even_odd_flag = 1;
	monomial_list[i].Qsq = &Qsw_pm_psi;
	monomial_list[i].Qp = &Qsw_plus_psi;
	monomial_list[i].Qm = &Qsw_minus_psi;
	init_swpm(VOLUME);
	if(g_proc_id == 0 && g_debug_level > 1) {
	  printf("# Initialised monomial of type CLOVERDETRATIO, no_monomials= %d\n", no_monomials);
	}
      }
      else if(monomial_list[i].type == CLOVERDETRATIORW) {
	monomial_list[i].accfunction = &cloverdetratio_rwacc;
	monomial_list[i].even_odd_flag = 1;
	monomial_list[i].Qsq = &Qsw_pm_psi;
	monomial_list[i].Qp = &Qsw_plus_psi;
	monomial_list[i].Qm = &Qsw_minus_psi;
	init_swpm(VOLUME);
	if(g_proc_id == 0 && g_debug_level > 1) {
	  printf("# Initialised monomial of type CLOVERDETRATIORW, no_monomials= %d, currently only available for reweighting!\n", no_monomials);
	}
      }
      else if(monomial_list[i].type == DETRATIO) {
	monomial_list[i].hbfunction = &detratio_heatbath;
	monomial_list[i].accfunction = &detratio_acc;
	monomial_list[i].derivativefunction = &detratio_derivative;
	monomial_list[i].Qsq = &Qtm_pm_psi;
	monomial_list[i].Qsq32 = &Qtm_pm_psi_32;	
	monomial_list[i].Qp = &Qtm_plus_psi;
	monomial_list[i].Qm = &Qtm_minus_psi;
	if(g_proc_id == 0 && g_debug_level > 1) {
	  printf("# Initialised monomial of type DETRATIO, no_monomials= %d\n", no_monomials);
	}
      }
      else if(monomial_list[i].type == POLY) {
	monomial_list[i].hbfunction = &poly_heatbath;
	monomial_list[i].accfunction = &poly_acc;
	monomial_list[i].derivativefunction = &poly_derivative;
	retval=init_poly_monomial(V,i);
	if(retval != 0) {
	  return retval;
	}
      	if(g_proc_id == 0 && g_debug_level > 1) {
	  printf("# Initialised monomial of type POLY, no_monomials= %d\n", no_monomials);
	}
      }
      else if(monomial_list[i].type == POLYDETRATIO) {
	monomial_list[i].hbfunction = &poly_heatbath;
	monomial_list[i].accfunction = &poly_acc;
	monomial_list[i].derivativefunction = &poly_derivative;
	monomial_list[i].MDPolyDetRatio = 1;
	retval=init_poly_monomial(V,i);
	if(retval!=0) return retval;
      	if(g_proc_id == 0 && g_debug_level > 1) {
	  printf("# Initialised monomial of type POLYDETRATIO, no_monomials= %d\n", no_monomials);
	}
      }
      else if(monomial_list[i].type == NDPOLY) {
	monomial_list[i].hbfunction = &ndpoly_heatbath;
	monomial_list[i].accfunction = &ndpoly_acc;
	monomial_list[i].derivativefunction = &ndpoly_derivative;
	monomial_list[i].even_odd_flag = 1;
	monomial_list[i].pf2 = __pf+no*V;
	no++;
	retval = init_ndpoly_monomial(i);
	if(g_proc_id == 0 && g_debug_level > 1) {
	  printf("# Initialised monomial of type NDPOLY, no_monomials= %d\n", no_monomials);
	}
      }
      else if(monomial_list[i].type == NDCLOVER) {
	init_swpm(VOLUME);
	monomial_list[i].hbfunction = &cloverndpoly_heatbath;
	monomial_list[i].accfunction = &cloverndpoly_acc;
	monomial_list[i].derivativefunction = &cloverndpoly_derivative;
	monomial_list[i].pf2 = __pf+no*V;
	monomial_list[i].even_odd_flag = 1;
	clovernd_monomials[no_clovernd_monomials] = i;
	no_clovernd_monomials++;
	//monomial_list[i].Qsq = &Qsw_pm_ndpsi;
	//monomial_list[i].Qp = &Qsw_ndpsi;
	//monomial_list[i].Qm = &Qsw_dagger_ndpsi;
	no++;
	retval = init_ndpoly_monomial(i);
	if(g_proc_id == 0 && g_debug_level > 1) {
	  printf("# Initialised monomial of type NDCLOVER, no_monomials= %d\n", no_monomials);
	}
      }
      else if(monomial_list[i].type == NDRAT) {
	monomial_list[i].hbfunction = &ndrat_heatbath;
	monomial_list[i].accfunction = &ndrat_acc;
	monomial_list[i].derivativefunction = &ndrat_derivative;
	monomial_list[i].even_odd_flag = 1;
	monomial_list[i].pf2 = __pf+no*V;
	no++;
	retval = init_ndrat_monomial(i);
	if(g_proc_id == 0 && g_debug_level > 1) {
	  printf("# Initialised monomial of type NDRAT, no_monomials= %d\n", no_monomials);
	}
      }
      else if(monomial_list[i].type == RAT) {
	monomial_list[i].hbfunction = &rat_heatbath;
	monomial_list[i].accfunction = &rat_acc;
	monomial_list[i].derivativefunction = &rat_derivative;
	monomial_list[i].even_odd_flag = 1;
	monomial_list[i].mu = 0.;
	monomial_list[i].Qsq = &Qtm_pm_psi;
	monomial_list[i].Qp = &Qtm_plus_psi;
	monomial_list[i].Qm = &Qtm_minus_psi;
	retval = init_ndrat_monomial(i);
	if(g_proc_id == 0 && g_debug_level > 1) {
	  printf("# Initialised monomial of type RAT, no_monomials= %d\n", no_monomials);
	}
      }
      else if(monomial_list[i].type == CLOVERRAT) {
	monomial_list[i].hbfunction = &rat_heatbath;
	monomial_list[i].accfunction = &rat_acc;
	monomial_list[i].derivativefunction = &rat_derivative;
	monomial_list[i].even_odd_flag = 1;
	monomial_list[i].mu = 0.;
	monomial_list[i].Qsq = &Qsw_pm_psi;
	monomial_list[i].Qp = &Qsw_plus_psi;
	monomial_list[i].Qm = &Qsw_minus_psi;
	init_swpm(VOLUME);
	if(monomial_list[i].trlog) {
	  clover_monomials[no_clover_monomials] = i;
	  no_clover_monomials++;
	}
	retval = init_ndrat_monomial(i);
	if(g_proc_id == 0 && g_debug_level > 1) {
	  printf("# Initialised monomial of type CLOVERRAT, no_monomials= %d\n", no_monomials);
	}
      }
      else if(monomial_list[i].type == NDCLOVERRAT) {
	init_swpm(VOLUME);
	monomial_list[i].hbfunction = &ndrat_heatbath;
	monomial_list[i].accfunction = &ndrat_acc;
	monomial_list[i].derivativefunction = &ndrat_derivative;
	monomial_list[i].even_odd_flag = 1;
	monomial_list[i].pf2 = __pf+no*V;
	no++;
	if(monomial_list[i].trlog) {
	  clovernd_monomials[no_clovernd_monomials] = i;
	  no_clovernd_monomials++;
	}
	retval = init_ndrat_monomial(i);
	if(g_proc_id == 0 && g_debug_level > 1) {
	  printf("# Initialised monomial of type NDCLOVERRAT, no_monomials= %d\n", no_monomials);
	}
      }
      else if(monomial_list[i].type == NDRATCOR) {
	monomial_list[i].hbfunction = &ndratcor_heatbath;
	monomial_list[i].accfunction = &ndratcor_acc;
	monomial_list[i].derivativefunction = NULL;
	monomial_list[i].even_odd_flag = 1;
	monomial_list[i].pf2 = __pf+no*V;
	monomial_list[i].rat.crange[0] = 0;
        monomial_list[i].rat.crange[1] = monomial_list[i].rat.order-1;
	
	no++;
	retval = init_ndrat_monomial(i);
	if(g_proc_id == 0 && g_debug_level > 1) {
	  printf("# Initialised monomial of type NDRATCOR, no_monomials= %d\n", no_monomials);
	}
      }
      else if(monomial_list[i].type == NDCLOVERRATCOR) {
	init_swpm(VOLUME);
	monomial_list[i].hbfunction = &ndratcor_heatbath;
	monomial_list[i].accfunction = &ndratcor_acc;
	monomial_list[i].derivativefunction = NULL;
	monomial_list[i].even_odd_flag = 1;
	monomial_list[i].pf2 = __pf+no*V;
	monomial_list[i].rat.crange[0] = 0;
        monomial_list[i].rat.crange[1] = monomial_list[i].rat.order-1;
	
	no++;
	retval = init_ndrat_monomial(i);
	if(g_proc_id == 0 && g_debug_level > 1) {
	  printf("# Initialised monomial of type NDCLOVERRATCOR, no_monomials= %d\n", no_monomials);
	}
      }
      else if(monomial_list[i].type == RATCOR) {
	monomial_list[i].hbfunction = &ratcor_heatbath;
	monomial_list[i].accfunction = &ratcor_acc;
	monomial_list[i].derivativefunction = NULL;
	monomial_list[i].even_odd_flag = 1;
	monomial_list[i].Qsq = &Qtm_pm_psi;
	monomial_list[i].Qp = &Qtm_plus_psi;
	monomial_list[i].Qm = &Qtm_minus_psi;	
	monomial_list[i].rat.crange[0] = 0;
        monomial_list[i].rat.crange[1] = monomial_list[i].rat.order-1;
	monomial_list[i].mu = 0.;
	retval = init_ndrat_monomial(i);
	if(g_proc_id == 0 && g_debug_level > 1) {
	  printf("# Initialised monomial of type RATCOR, no_monomials= %d\n", no_monomials);
	}
      }
      else if(monomial_list[i].type == CLOVERRATCOR) {
	init_swpm(VOLUME);
	monomial_list[i].hbfunction = &ratcor_heatbath;
	monomial_list[i].accfunction = &ratcor_acc;
	monomial_list[i].derivativefunction = NULL;
	monomial_list[i].even_odd_flag = 1;
	monomial_list[i].Qsq = &Qsw_pm_psi;
	monomial_list[i].Qp = &Qsw_plus_psi;
	monomial_list[i].Qm = &Qsw_minus_psi;
	monomial_list[i].mu = 0.;
	monomial_list[i].rat.crange[0] = 0;
        monomial_list[i].rat.crange[1] = monomial_list[i].rat.order-1;
	retval = init_ndrat_monomial(i);
	if(g_proc_id == 0 && g_debug_level > 1) {
	  printf("# Initialised monomial of type CLOVERRATCOR, no_monomials= %d\n", no_monomials);
	}
      }
      else if(monomial_list[i].type == NDDETRATIO) {
	monomial_list[i].hbfunction = &dummy_heatbath;
	monomial_list[i].accfunction = &nddetratio_acc;
	monomial_list[i].derivativefunction = NULL;
	monomial_list[i].pf2 = __pf+no*V;
	monomial_list[i].timescale = -5;
	no++;
	if(g_proc_id == 0 && g_debug_level > 1) {
	  printf("# Initialised monomial of type NDDETRATIO, no_monomials= %d, currently only available for reweighting!\n", no_monomials);
	}
      }
      else if(monomial_list[i].type == NDCLOVERDETRATIO) {
	monomial_list[i].hbfunction = &dummy_heatbath;
	monomial_list[i].accfunction = &nddetratio_acc;
	monomial_list[i].derivativefunction = NULL;
	monomial_list[i].pf2 = __pf+no*V;
	monomial_list[i].timescale = -5;
	no++;
	if(g_proc_id == 0 && g_debug_level > 1) {
	  printf("# Initialised monomial of type NDCLOVERDETRATIO, no_monomials= %d, currently only available for reweighting!\n", no_monomials);
	}
      }
    }
    else {
      monomial_list[i].pf = NULL;
      if(no_gauge_monomials > 0) {
	fprintf(stderr, "maximal number of gauge monomials exceeded!\n");
	exit(-1);
      }      
      else if(monomial_list[i].type == GAUGE) {
	monomial_list[i].hbfunction = &gauge_heatbath;
	monomial_list[i].accfunction = &gauge_acc;
	monomial_list[i].derivativefunction = &gauge_derivative;
	no_gauge_monomials++;
	if(fabs( monomial_list[i].glambda) > 0) {
	  monomial_list[i].derivativefunction = &gauge_EMderivative;
	}
	if(!monomial_list[i].use_rectangles) {
	  monomial_list[i].c1 = 0.;
	  monomial_list[i].c0 = 1.;
	}
	g_rgi_C1 = monomial_list[i].c1;
	monomial_list[i].c0 = 1. - 8.*monomial_list[i].c1;
	g_rgi_C0 = monomial_list[i].c0;
	if(g_proc_id == 0 && g_debug_level > 1) {
	  printf("# Initialised monomial of type GAUGE, no_monomials= %d\n", no_monomials);
	}
      }
    }
    monomial_list[i].id = i;
    monomial_list[i].even_odd_flag = even_odd_flag;
  }
  /* initialize clovertrlog and cloverndtrlog monomials for all clover and clovernd monomials*/
  if( even_odd_flag ) {
    for( int j = 0; j < no_clover_monomials; j++ ) {
      monomial_list[no_monomials].type = CLOVERTRLOG;
      strcpy( monomial_list[no_monomials].name, "CLOVERTRLOG");
      add_monomial(CLOVERTRLOG);
      monomial_list[no_monomials-1].pf = NULL;
      monomial_list[no_monomials-1].id = no_monomials-1;
      monomial_list[no_monomials-1].rngrepro = reproduce_randomnumber_flag;
      // set the parameters according to cloverdet monomial
      // this need alltogether a more general approach
      monomial_list[no_monomials-1].c_sw = monomial_list[clover_monomials[j]].c_sw;
      monomial_list[no_monomials-1].mu = monomial_list[clover_monomials[j]].mu;
      monomial_list[no_monomials-1].kappa = monomial_list[clover_monomials[j]].kappa;
      monomial_list[no_monomials-1].hbfunction = &clover_trlog_heatbath;
      monomial_list[no_monomials-1].accfunction = &clover_trlog_acc;
      monomial_list[no_monomials-1].derivativefunction = NULL;
      monomial_list[no_monomials-1].timescale = 0;
      monomial_list[no_monomials-1].even_odd_flag = even_odd_flag;
      if(g_proc_id == 0 && g_debug_level > 1) {
        printf("# Initialised clover_trlog_monomial, no_monomials= %d\n", no_monomials);
      }
    }
    for( int j = 0; j < no_clovernd_monomials; j++ ) { 
      monomial_list[no_monomials].type = CLOVERNDTRLOG;
      strcpy( monomial_list[no_monomials].name, "CLOVERNDTRLOG");
      add_monomial(CLOVERNDTRLOG);
      monomial_list[no_monomials-1].pf = NULL;
      monomial_list[no_monomials-1].id = no_monomials-1;
      monomial_list[no_monomials-1].rngrepro = reproduce_randomnumber_flag;
      // set the parameters according to cloverdet monomial
      // this need alltogether a more general approach
      monomial_list[no_monomials-1].c_sw = monomial_list[clovernd_monomials[j]].c_sw;
      monomial_list[no_monomials-1].mubar = monomial_list[clovernd_monomials[j]].mubar;
      monomial_list[no_monomials-1].epsbar = monomial_list[clovernd_monomials[j]].epsbar;
      monomial_list[no_monomials-1].kappa = monomial_list[clovernd_monomials[j]].kappa;
      monomial_list[no_monomials-1].hbfunction = &clovernd_trlog_heatbath;
      monomial_list[no_monomials-1].accfunction = &clovernd_trlog_acc;
      monomial_list[no_monomials-1].derivativefunction = NULL;
      monomial_list[no_monomials-1].timescale = 0;
      monomial_list[no_monomials-1].even_odd_flag = 1;
      if(g_proc_id == 0 && g_debug_level > 1) {
        printf("# Initialised clovernd_trlog_monomial, no_monomials= %d\n", no_monomials);
      }
    }
  }
  return(0);
}

void free_monomials() {
  
  free(_pf);
  return;
}


int init_poly_monomial(const int V, const int id){
  
  monomial * mnl = &monomial_list[id];
  int i,j,k;
  FILE* rootsFile=NULL;
  char title[101];
  char filename[257];
  FILE* constFile;
  int errcode;
  double eps;

  spinor *_pf=(spinor*)NULL;
  
  if((void*)(_pf = (spinor*)calloc((mnl->MDPolyDegree/2+2)*V+1, sizeof(spinor))) == NULL) {
    printf ("malloc errno in init_poly_monomial pf fields: %d\n",errno); 
    errno = 0;
    return(1);
  }
  
  if((void*)(mnl->MDPoly_chi_spinor_fields=(spinor**)calloc(mnl->MDPolyDegree/2+2,sizeof(spinor*))) ==NULL ){
    printf ("malloc errno in init_poly_monomial pf fields: %d\n",errno); 
    errno = 0;
    return(2);
  }
  
  (mnl->MDPoly_chi_spinor_fields)[0] = (spinor*)(((unsigned long int)(_pf)+ALIGN_BASE)&~ALIGN_BASE);
  
  for(i = 1; i < (mnl->MDPolyDegree/2+2); i++){
    mnl->MDPoly_chi_spinor_fields[i] = mnl->MDPoly_chi_spinor_fields[i-1]+V;
  }
  
  if(strlen(monomial_list[id].MDPolyRootsFile)==0){
    sprintf(monomial_list[id].MDPolyRootsFile,
	    "%s_deg_%d_eps_%1.16e.roots",
	    "1overX_poly",
	    monomial_list[id].MDPolyDegree,
	    monomial_list[id].MDPolyLmin/monomial_list[id].MDPolyLmax
	    );
    fprintf(stderr,"Warning you didnt specify a rootsfilename -> guessing:\n%s\n",filename);
  }
  if(monomial_list[id].MDPolyLocNormConst==-1.0){
    eps=monomial_list[id].MDPolyLmin/monomial_list[id].MDPolyLmax;
    sprintf(filename,
	    "%s_deg_%d_eps_%1.16e.const",
	    "1overX_poly",
	    monomial_list[id].MDPolyDegree,
	    eps
	    );
    fprintf(stderr,"Warning you didnt specify a local normalization: trying to read it from\n%s\n",filename);
    if((constFile=fopen(filename,"r"))!=NULL) {
      errcode = fscanf(constFile,"%lf\n",&(mnl->MDPolyLocNormConst));
      fclose(constFile);
      fprintf(stderr, "normierung local succesfully read -> lnc =  %e \n", mnl->MDPolyLocNormConst);
    } 
    else {
      fprintf(stderr,"Reading local normalization from file FAILED\n Borting Ab\n");
#ifdef TM_USE_MPI
      MPI_Finalize();
#endif
      exit(6);
    }
  }
  
  /* read in the roots from the given file */
  
  if((void*)(mnl->MDPolyRoots=(_Complex double*)calloc(mnl->MDPolyDegree,sizeof(_Complex double))) ==NULL ){
    printf ("malloc errno in init_poly_monomial roots array: %d\n",errno); 
    errno = 0;
    return(3);
  }
  
  printf("reading roots...!\n");
  if((rootsFile=fopen(mnl->MDPolyRootsFile,"r")) != (FILE*)NULL) {
    if (fgets(title, 100, rootsFile) == NULL) {
      fprintf(stderr, "Cant read Roots file: %s Aborting...\n", mnl->MDPolyRootsFile);
#ifdef TM_USE_MPI
      MPI_Finalize();
#endif
      exit(6);
    }
    
    /* Here we read in the 2n roots needed for the polinomial in sqrt(s) */
    for(j = 0; j < (mnl->MDPolyDegree); j++) {
      errcode = fscanf(rootsFile," %d %lf %lf \n", &k, (double*)&(mnl->MDPolyRoots[j]), (double*)&(mnl->MDPolyRoots[j]) + 1);
    }
    fclose(rootsFile);
  }
  else {
    fprintf(stderr, "Roots File %s is missing! Aborting...\n", mnl->MDPolyRootsFile );
#ifdef TM_USE_MPI
    MPI_Finalize();
#endif
    exit(6);
  }
  
  if(g_proc_id == 0 && g_debug_level > 2) {
    printf("# the root are:\n");
    for(j=0; j<(mnl->MDPolyDegree); j++){
      printf("# %lf %lf\n",  creal(mnl->MDPolyRoots[j]), cimag(mnl->MDPolyRoots[j]));
    }
  }
  
  return 0;
  
}

void dummy_derivative(const int id, hamiltonian_field_t * const hf) {
  if(g_proc_id == 0) {
    fprintf(stderr, "dummy_derivative was called. Was that really intended?\n");
    fprintf(stderr, "callers monomial ID was %d\n", id);
  }
  return;
}

void dummy_heatbath(const int id, hamiltonian_field_t * const hf) {
  if(g_proc_id == 0) {
    fprintf(stderr, "dummy_heatbath was called. Was that really intended?\n");
    fprintf(stderr, "callers monomial ID was %d\n", id);
  }
  return;
}

double dummy_acc(const int id, hamiltonian_field_t * const hf) {
  if(g_proc_id == 0) {
    fprintf(stderr, "dummy_acc was called. Was that really intended?\n");
    fprintf(stderr, "callers monomial ID was %d\n", id);
  }
  return(0.);
}

void mnl_backup_restore_globals(const backup_restore_t mode){
  static double backup_kappa;
  static double backup_mu;
  static double backup_mu1;
  static double backup_mu2;
  static double backup_mu3;
  static double backup_c_sw;
  static double backup_mubar;
  static double backup_epsbar;
  if( mode == TM_BACKUP_GLOBALS ){
    backup_kappa  = g_kappa;
    backup_c_sw   = g_c_sw;
    backup_mu     = g_mu;
    backup_mu1    = g_mu1;
    backup_mu2    = g_mu2;
    backup_mu3    = g_mu3;
    backup_mubar  = g_mubar;
    backup_epsbar = g_epsbar;
  } else {
    g_kappa  = backup_kappa;
    g_c_sw   = backup_c_sw;
    g_mu     = backup_mu;
    g_mu1    = backup_mu1;
    g_mu2    = backup_mu2;
    g_mu3    = backup_mu3;
    g_mubar  = backup_mubar;
    g_epsbar = backup_epsbar;
    boundary(g_kappa);
  }
}

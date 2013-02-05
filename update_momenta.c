/***********************************************************************
 *
 * Copyright (C) 2001 Martin Hasebusch
 *               2002,2003,2004,2005,2006,2007,2008,2012 Carsten Urbach
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
#include <time.h>
#include "global.h"
#include "su3.h"
#include "su3adj.h"
#include "su3spinor.h"
#include "monomial/monomial.h"
#include "xchange/xchange.h"
#include "operator/clover_leaf.h"
#include "read_input.h"
#include "hamiltonian_field.h"
#include "update_momenta.h"
#include "gettime.h"

#include <buffers/utils.h>

#include <dirty_shameful_business.h>

/* Updates the momenta: equation 16 of Gottlieb */
void update_momenta(int * mnllist, double step, const int no, hamiltonian_field_t * const hf)
{
  double atime=0., etime=0.;
  int *relevant_smearings = malloc(no_smearings_monomial * sizeof(int));
  int no_relevant_smearings = 0;
  
  for(int k = 0; k < no; k++)
  {
    int current_idx = 0;
    while ((current_idx < no_relevant_smearings) && (monomial_list[ mnllist[k] ].smearing != relevant_smearings[current_idx]))
      ++current_idx;
    if (current_idx == no_relevant_smearings)
    {
      relevant_smearings[current_idx] = monomial_list[ mnllist[k] ].smearing;
      ++no_relevant_smearings;
    }
  }
  
  adjoint_field_t tmp_derivative = get_adjoint_field();
  
  zero_adjoint_field(&df);

  /* FIXME Further debugging -- replace the field by something that is manageable analytically */
  for (int loc = 0; loc < (VOLUMEPLUSRAND + g_dbw2rand); ++loc)
    for (int nu = 0; nu < 4; ++nu)
    {
      _su3_one(g_gf[loc][nu]);
    }

//  /* FIXME Enter some arbitrary random SU(3). Random to avoid accidental agreements. */
//   g_gf[0][0].c00 = -0.281674601208749 - 0.196755446881004 * I;
//   g_gf[0][0].c01 =  0.097906948456113 + 0.909672287338323 * I;
//   g_gf[0][0].c02 =  0.085554250669676 + 0.193746593367662 * I;
//   g_gf[0][0].c10 =  0.255396409335753 - 0.282878184149928 * I;
//   g_gf[0][0].c11 = -0.057616138763708 - 0.068813232780069 * I;
//   g_gf[0][0].c12 =  0.919668327134945 + 0.030131290107975 * I;
//   g_gf[0][0].c20 =  0.054229296057200 - 0.856597048580097 * I;
//   g_gf[0][0].c21 =  0.329775646781283 - 0.214717725004200 * I;
//   g_gf[0][0].c22 = -0.279641903356002 + 0.173912067032302 * I;

  /* FIXME Enter the simplest deviation from identity that will work */
   g_gf[0][0].c00 = I;
   g_gf[0][0].c11 = -I;
  
  //rnd_gauge_trafo(&g_gf, g_gf);

  ohnohack_remap_df0(tmp_derivative); /* FIXME Such that we can aggregate results per smearing type. */
  for (int s_ctr = 0; s_ctr < no_relevant_smearings; ++s_ctr)
  {
    int s_type = relevant_smearings[s_ctr];

    smear(smearing_control_monomial[s_type], g_gf);
    zero_adjoint_field(&tmp_derivative);
    ohnohack_remap_g_gauge_field(smearing_control_monomial[s_type]->result);

    for(int k = 0; k < no; k++)
    {
      if (monomial_list[ mnllist[k] ].smearing == s_type)
      {
        if(monomial_list[ mnllist[k] ].derivativefunction != NULL)
        {
#ifdef MPI
          atime = MPI_Wtime();
#else
          atime = (double)clock()/(double)(CLOCKS_PER_SEC);
#endif

          monomial_list[ mnllist[k] ].derivativefunction(mnllist[k], hf);

#ifdef MPI
          etime = MPI_Wtime();
#else
          etime = (double)clock()/(double)(CLOCKS_PER_SEC);
#endif
        }
      }
    }

    smear_forces(smearing_control_monomial[s_type], tmp_derivative);

    for(int i = 0; i < (VOLUMEPLUSRAND + g_dbw2rand); ++i)
    { 
      for(int mu = 0; mu < 4; ++mu)
      {
        _add_su3adj(df[i][mu], smearing_control_monomial[s_type]->force_result[i][mu]);
      }
    }
  }
  ohnohack_remap_df0(df);
  ohnohack_remap_g_gauge_field(g_gf);

#ifdef MPI
  xchange_deri(hf->derivative);
#endif

  /* Debugging code inserted here */
  su3adj num_eval;
  int xref = 0;
  int muref = 0;
  calculate_forces_numerically(&num_eval, xref, muref, mnllist, no);
  double *ar_df = (double*)&df[xref][muref];
  double *ar_ne = (double*)&num_eval;
  fprintf(stderr, "[DEBUG] Comparison of force calculation at (%d, %d)!\n", xref, muref);
  fprintf(stderr, "        Involves monomial types: ");
  for (int monnum = 0; monnum < no; ++monnum)
    fprintf(stderr, "%d ", monomial_list[ mnllist[monnum] ].type);
  fprintf(stderr, "\n");
  fprintf(stderr, "           Analytical  vs.  numerical\n");
  for (int component = 0; component < 8; ++component)
    fprintf(stderr, "    [%d]  %+12.7f <-> %+12.7f\n", component, ar_df[component], ar_ne[component]);

#ifdef OMP
#pragma omp parallel for
#endif
  for(int i = 0; i < VOLUME; i++) {
    for(int mu = 0; mu < 4; mu++) {
      /* the minus comes from an extra minus in trace_lambda */
      _su3adj_minus_const_times_su3adj(hf->momenta[i][mu], step, hf->derivative[i][mu]); 
    }
  }
  return_adjoint_field(&tmp_derivative);
  free(relevant_smearings);
  return;
}


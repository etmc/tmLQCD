/***********************************************************************
 *
 * Copyright (C) 2008 Albert Deuzeman, Siebren Reker, Carsten Urbach
 *                    Claude Tadonki
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
#include <config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "global.h"
#include "gettime.h"
#include "su3.h"
#include <complex.h>
#include "read_input.h"
#include "start.h"
#include "gamma.h"
#include "ranlxs.h"
#include "init/init.h"
#include "operator/D_psi.h"
#include "operator/tm_operators.h"
#include "operator/clover_leaf.h"
#include "operator/clovertm_operators.h"
#include "poly_precon.h"
#include "Msap.h"
#include "gmres_precon.h"
#include "linalg_eo.h"
#include "gram-schmidt.h"
#include "lu_solve.h"
#include "block.h"
#include "little_D.h"
#include "gcr4complex.h"
#include "cgne4complex.h"
#include "boundary.h"
#include <io/params.h>
#include <io/gauge.h>
#include <io/spinor.h>
#include <io/utils.h>
#include "solver/solver.h"
#include "solver_field.h"
#include "dfl_projector.h"
#include "generate_dfl_subspace.h"

int init_little_dfl_subspace(const int N_s);
void compute_little_little_D(const int N_s);

spinor ** dfl_fields = NULL;
static spinor * _dfl_fields = NULL;
_Complex double ** little_dfl_fields = NULL;
static _Complex double *_little_dfl_fields = NULL;
_Complex double ** little_dfl_fields_eo = NULL;
static _Complex double *_little_dfl_fields_eo = NULL;
static int init_subspace = 0;
static int init_little_subspace = 0;

static void random_fields(const int Ns) {

  for (int i = 0; i < Ns; i++) {
    random_spinor_field_lexic(dfl_fields[i], 0, RN_PM1UNIF);
  }
  return;
}


// this routine updates the deflation subspace

int update_dfl_subspace(const int Ns, const int N, const int Nsmooth) {
  int vpr = VOLUMEPLUSRAND*sizeof(spinor)/sizeof(_Complex double),
    vol = VOLUME*sizeof(spinor)/sizeof(_Complex double);
  double musave = g_mu;
  double kappasave = g_kappa;
  
  if(kappa_dflgen > 0) {
    g_kappa = kappa_dflgen;
  }
  if(mu_dflgen > -10) {
    g_mu = mu_dflgen;
    if(g_mu*musave < 0) g_mu *= -1.;
  }
  boundary(g_kappa);
  
  double nrm; 
  for(int j = 0; j < Nsmooth; j++) {
    //if(j - loop_SAP > 4) {
    //  usePL = usePLsave;
    //}
    // build little D automatically when little_D is called and dfl_subspace_updated == 1
    for (int i = 0; i < Ns; i++) {
      /* add it to the basis */
      split_global_field_GEN_ID(block_list, i, dfl_fields[i], nb_blocks);
    }
    /* perform local orthonormalization */
    for(int i = 0; i < nb_blocks; i++) {
      block_orthonormalize(block_list+i);
    }
    dfl_subspace_updated = 1;
    compute_little_little_D(Ns);
    
    // why can't I use g_spinor_field[1] in mg_precon??
    for(int i = 0; i < Ns; i++) {
      g_sloppy_precision = 1;

      // v <- C v
      mg_precon(g_spinor_field[0], dfl_fields[i]);
      assign(dfl_fields[i], g_spinor_field[0], VOLUME);
      
      g_sloppy_precision = 0;
      ModifiedGS((_Complex double*)g_spinor_field[0], vol, i, (_Complex double*)dfl_fields[0], vpr);
      nrm = sqrt(square_norm(g_spinor_field[0], N, 1));
      mul_r(dfl_fields[i], 1./nrm, g_spinor_field[0], N);
      
    }
  }
  for(int i = 0; i < Ns; i++) {
    /* add it to the basis */
    split_global_field_GEN_ID(block_list, i, dfl_fields[i], nb_blocks);
  }
  /* perform local orthonormalization */
  for(int i = 0; i < nb_blocks; i++) {
    block_orthonormalize(block_list+i);
  }
  dfl_subspace_updated = 1;
  g_mu = musave;
  g_kappa = kappasave;
  boundary(g_kappa);

  return(0);
}

// this routine generates the deflation subspace

int generate_dfl_subspace(const int Ns, const int N, const int repro) {
  int p=0, vpr = VOLUMEPLUSRAND*sizeof(spinor)/sizeof(_Complex double),
    vol = VOLUME*sizeof(spinor)/sizeof(_Complex double);
  const int loop_SAP = 3;
  spinor **psi;
  double nrm, atime, etime;
  _Complex double s;
  _Complex double * work;
  double musave = g_mu;
  double kappasave = g_kappa;
  int usePLsave = usePL;
  spinor ** work_fields = NULL;
  const int nr_wf = 2;
  g_mu2 = 0.;
  if(kappa_dflgen > 0) {
    g_kappa = kappa_dflgen;
  }
  if(mu_dflgen > -10) {
    g_mu = mu_dflgen;
    if(g_mu*musave < 0) g_mu *= -1.;
  }
  if(g_c_sw > 0) {
    if (g_cart_id == 0 && g_debug_level > 1) {
      printf("#\n# csw = %e, computing clover leafs\n", g_c_sw);
    }
    init_sw_fields(VOLUME);
    sw_term( (const su3**) g_gauge_field, g_kappa, g_c_sw);
    // this must be EE = 0, this is needed for Msap_eo!?
    sw_invert(0, g_mu);
    /* now copy double sw and sw_inv fields to 32bit versions */
    copy_32_sw_fields();
  }

  // currently set to 0 during subspace generation
  usePL = 0;
  atime = gettime();

  init_solver_field(&work_fields, VOLUMEPLUSRAND, nr_wf);
  work = (_Complex double*)malloc(nb_blocks*9*Ns*sizeof(_Complex double));
  psi = (spinor **)calloc(nb_blocks, sizeof(spinor *));
  psi[0] = calloc(VOLUME + nb_blocks, sizeof(spinor));
  for(int i = 1; i < nb_blocks; i++) psi[i] = psi[i-1] + (VOLUME / nb_blocks) + 1;

  if((g_proc_id == 0) && (g_debug_level > 0)) {
    printf("# Initialising subspaces\n");
    fflush(stdout);
  }

  if(init_subspace == 0) p = init_dfl_subspace(Ns);

  if(init_little_subspace == 0) p = init_little_dfl_subspace(Ns);

  if((g_proc_id == 0) && (g_debug_level > 0)) {
    printf("# Generating random fields...");
  }
  double ta = gettime();
  random_fields(Ns);
  double tb = gettime();
  if((g_proc_id == 0) && (g_debug_level > 0)) {
    printf(" done in %e seconds\n", tb-ta);
    fflush(stdout);
  }

  boundary(g_kappa);

  if((g_proc_id == 0) && (p < Ns) && (g_debug_level > 0)) {
    printf("# Compute approximate eigenvectors from scratch\n");
    printf("# Using kappa= %e and mu = %e for the subspace generation\n", g_kappa, g_mu/g_kappa/2.);
  }

  for(int j = 0; j < loop_SAP; j++) {
    for(int i = 0; i < Ns; i++) {
      zero_spinor_field(g_spinor_field[0], VOLUME);  
      g_sloppy_precision = 1;
      Msap_eo(g_spinor_field[0], dfl_fields[i], NcycleMsap_dflgen, NiterMsap_dflgen); 
      
      assign(dfl_fields[i], g_spinor_field[0], VOLUME);
      
      g_sloppy_precision = 0;
      ModifiedGS((_Complex double*)g_spinor_field[0], vol, i, (_Complex double*)dfl_fields[0], vpr);
      nrm = sqrt(square_norm(g_spinor_field[0], N, 1));
      mul_r(dfl_fields[i], 1./nrm, g_spinor_field[0], N);
    }
  }

  update_dfl_subspace(Ns, N, NsmoothMsap_dflgen-loop_SAP);
  
  compute_little_D(0);
  compute_little_little_D(Ns);
  dfl_subspace_updated = 0;
  g_mu = musave;
  g_kappa = kappasave;
  boundary(g_kappa);
  usePL = usePLsave;
  if(g_debug_level > 0 && g_proc_id == 0) {
    printf("# Switched to target parameters kappa= %e, mu=%e\n", g_kappa, g_mu/2/g_kappa);
  }
  if(g_debug_level > 1) {
    for(int i = 0; i < Ns; i++) {
      /* test quality */
      D_psi(work_fields[0], dfl_fields[i]);
      nrm = sqrt(square_norm(work_fields[0], N, 1));
      if(g_proc_id == 0) {
	printf(" ||D psi_%d||/||psi_%d|| = %1.15e\n", i, i, nrm);
      }
    }
  }
  
  if((g_proc_id == 0) && (g_debug_level > 0)) {
    printf("# Approximate eigenvectors generated\n");
  }
  
  g_mu1 = 0.;
  g_mu2 = 0.;

  g_sloppy_precision = 0;

  if(g_debug_level > 4) {
    printf("Checking orthonormality of dfl_fields\n");
    for(int i = 0; i < Ns; i++) {
      for(int j = 0; j < Ns; j++) {
        s = scalar_prod(dfl_fields[i], dfl_fields[j], N, 1);
        if(g_proc_id == 0) {
          printf("<%d, %d> = %1.3e +i %1.3e\n", i, j, creal(s), cimag(s));
        }
      }
    }
  }

  dfl_subspace_updated = 1;

  compute_little_little_D(Ns);
  
  etime = gettime();

  if(g_proc_id == 0) {
    printf("# time for subspace generation %1.3e s\n", etime-atime);
    fflush(stdout);
  }

  finalize_solver(work_fields, nr_wf);
  free_dfl_subspace();
  free(work);
  free(psi[0]);
  free(psi);

  /* Cross-checks */
  if (g_debug_level > 3) {
    check_projectors(reproduce_randomnumber_flag);
    check_local_D(reproduce_randomnumber_flag);
    check_little_D_inversion(reproduce_randomnumber_flag);
  }

  return(0);
}

int generate_dfl_subspace_free(const int Ns, const int N) {
  int i,j, vpr = VOLUMEPLUSRAND*sizeof(spinor)/sizeof(_Complex double), 
    vol = VOLUME*sizeof(spinor)/sizeof(_Complex double);
  double nrm;
  _Complex double s;
  spinor ** work_fields = NULL;
  const int nr_wf = 1;
  init_solver_field(&work_fields, VOLUMEPLUSRAND, nr_wf);

  if(init_subspace == 0) init_dfl_subspace(Ns);

  for(i = 0; i < 12; i++) {
    constant_spinor_field(dfl_fields[i], i, N);
    ModifiedGS((_Complex double*)dfl_fields[i], vol, i, (_Complex double*)dfl_fields[0], vpr);
    nrm = sqrt(square_norm(dfl_fields[i], N, 1));
    mul_r(dfl_fields[i], 1./nrm, dfl_fields[i], N);

    /* test quality */
    if(g_debug_level > -1) {
      D_psi(work_fields[0], dfl_fields[i]);
      nrm = sqrt(square_norm(work_fields[0], N, 1));
      if(g_proc_id == 0) {
        printf(" ||D psi_%d||/||psi_%d|| = %1.5e\n", i, i, nrm); 
      }
    }
  }

  if(g_debug_level > 4) {
    for(i = 0; i < 12; i++) {
      for(j = 0; j < 12; j++) {
        s = scalar_prod(dfl_fields[i], dfl_fields[j], N, 1);
        if(g_proc_id == 0) {
          printf("<%d, %d> = %1.3e +i %1.3e\n", i, j, creal(s), cimag(s));
        }
      }
    }
  }
  finalize_solver(work_fields, nr_wf);
  return(0);
}


void compute_little_little_D(const int Ns) {
  int i_o;
  _Complex double s;
  if(usePL) {
    double atime = gettime();
    _Complex double * work = (_Complex double*)malloc(nb_blocks*9*Ns*sizeof(_Complex double));
    spinor ** psi = (spinor **)calloc(nb_blocks, sizeof(spinor *));
    psi[0] = calloc(VOLUME + nb_blocks, sizeof(spinor));
    for(int i = 1; i < nb_blocks; i++) psi[i] = psi[i-1] + (VOLUME / nb_blocks) + 1;
    
    if(!little_evenodd) {
      // compute the little little basis
      for(int j = 0; j < Ns; j++) {
	for(int i = 0; i < nb_blocks*9*Ns; i++) {
	  (little_dfl_fields[j][i]) = 0.0;
	  (work[i]) = 0.0;
	}
      }
      
      for(int i = 0; i < Ns; i++) {
	split_global_field_GEN(psi, dfl_fields[i], nb_blocks);
	// now take the local scalar products
	for(int j = 0; j < Ns; j++) {
	  for(int blk = 0; blk < nb_blocks; blk++) {
	    little_dfl_fields[i][j + blk*Ns] = scalar_prod(block_list[blk].basis[j], psi[blk], block_list[0].volume, 0);
	  }
	}
      }
      
      // orthonormalise 
      for(int i = 0; i < Ns; i++) {
	for(int j = 0; j < i; j++) {
	  s = lscalar_prod(little_dfl_fields[j], little_dfl_fields[i], nb_blocks*Ns, 1);
	  lassign_diff_mul(little_dfl_fields[i], little_dfl_fields[j], s, nb_blocks*Ns);
	}
	s = lsquare_norm(little_dfl_fields[i], nb_blocks*Ns, 1);
	lmul_r(little_dfl_fields[i], 1./sqrt(creal(s)), little_dfl_fields[i], nb_blocks*Ns);
      }
      if(g_debug_level > 4) {
	printf("Checking orthonormality of little dfl fields\n");
	for(int i = 0; i < Ns; i++) {
	  for(int j = 0; j < Ns; j++) {
	    s = lscalar_prod(little_dfl_fields[i], little_dfl_fields[j], nb_blocks*Ns, 1);
	    if(g_proc_id == 0) {
	      printf("<%d, %d> = %1.3e +i %1.3e\n", i, j, creal(s), cimag(s));
	    }
	  }
	}
      }
      
      for(int i = 0; i < Ns; i++) {
	little_D(work, little_dfl_fields[i]);
	for(int j = 0; j < Ns; j++) {
	  little_A[i * Ns + j]  = lscalar_prod(little_dfl_fields[j], work, nb_blocks*Ns, 1);
	  if(g_proc_id == 0 && g_debug_level > 4) {
	    printf("%1.3e %1.3ei, ", creal(little_A[i * Ns + j]), cimag(little_A[i * Ns + j]));
	  }
	}
	if(g_proc_id == 0 && g_debug_level > 4) printf("\n");
      }
      if(g_proc_id == 0 && g_debug_level > 4) printf("\n");
      LUInvert(Ns, little_A, Ns);
      // inverse of little little D now in little_A 
    }
    else { // if little_evenodd = true
      // compute the little little basis for the even/odd case 
      for(int j = 0; j < Ns; j++) {
	for(int i = 0; i < nb_blocks*9*Ns; i++) {
	  (little_dfl_fields_eo[j][i]) = 0.0;
	  (work[i]) = 0.0;
	}
      }
      
      for(int i = 0; i < Ns; i++) {
	split_global_field_GEN(psi, dfl_fields[i], nb_blocks);
	// now take the local scalar products
	for(int j = 0; j < Ns; j++) {
	  i_o=0;
	  for(int blk = 0; blk < nb_blocks; blk++) {
	    if (block_list[blk].evenodd==1) {
	      little_dfl_fields_eo[i][j + (nb_blocks/2+i_o)*Ns] = scalar_prod(block_list[blk].basis[j], psi[blk], block_list[0].volume, 0);
	      i_o++;
	    }       
	  }
	}
      }  
      
      // orthonormalise
      for(int i = 0; i < Ns; i++) {
	for (int j = 0; j < i; j++) {
	  s = lscalar_prod(little_dfl_fields_eo[j], little_dfl_fields_eo[i], nb_blocks*Ns, 1);
	  lassign_diff_mul(little_dfl_fields_eo[i], little_dfl_fields_eo[j], s, nb_blocks*Ns);
	}
	s = lsquare_norm(little_dfl_fields_eo[i], nb_blocks*Ns, 1);
	lmul_r(little_dfl_fields_eo[i], 1./sqrt(creal(s)), little_dfl_fields_eo[i], nb_blocks*Ns);
      }
      if(g_debug_level > 4) {
	printf("Checking orthonormality of littel_dfl_fields_eo\n");
	for(int i = 0; i < Ns; i++) {
	  for(int j = 0; j < Ns; j++) {
	    s = lscalar_prod(little_dfl_fields_eo[i], little_dfl_fields_eo[j], nb_blocks*Ns, 1);
	    if(g_proc_id == 0) {
	      printf("<%d, %d> = %1.3e +i %1.3e\n", i, j, creal(s), cimag(s));
	    }
	  }
	}
      }
      
      for(int i = 0; i < Ns; i++) {  
	little_D_sym(work, little_dfl_fields_eo[i]);
	for(int j = 0; j < Ns; j++) {
	  little_A_eo[i * Ns + j]  = lscalar_prod(little_dfl_fields_eo[j], work, nb_blocks*Ns, 1);
	  if(g_proc_id == 0 && g_debug_level > 4) {
	    printf("%1.3e %1.3ei, ", creal(little_A_eo[i * Ns + j]), cimag(little_A_eo[i * Ns + j])); 
	  }
	}
	if(g_proc_id == 0 && g_debug_level > 4) printf("\n");
      }
      if(g_proc_id == 0 && g_debug_level > 4) printf("\n");
      LUInvert(Ns, little_A_eo, Ns);
      // inverse of eo little little D now in little_A_eo
    }
    
    double etime = gettime();
    if(g_proc_id == 0 && g_debug_level > 0) {
      printf("# time for little little D computation %1.3e s\n", etime-atime);
      fflush(stdout);
    }
    free(work);
    free(psi[0]);
    free(psi);
  }

  return;
}



int init_little_dfl_subspace(const int N_s) {
  int i;
  if(init_little_subspace == 0) {
    if((void*)(_little_dfl_fields = (_Complex double*)calloc((N_s)*nb_blocks*9*N_s+4, sizeof(_Complex double))) == NULL) {
      return(1);
    }
    if((void*)(little_dfl_fields = (_Complex double**)calloc(N_s, sizeof(_Complex double*))) == NULL) {
      return(1);
    }
    if((void*)(_little_dfl_fields_eo = (_Complex double*)calloc((N_s)*nb_blocks*9*N_s+4, sizeof(_Complex double))) == NULL) {
      return(1);
    }
    if((void*)(little_dfl_fields_eo = (_Complex double**)calloc(N_s, sizeof(_Complex double*))) == NULL) {
      return(1);
    }
    little_dfl_fields[0] = (_Complex double*)(((unsigned long int)(_little_dfl_fields)+ALIGN_BASE)&~ALIGN_BASE);
    little_dfl_fields_eo[0] = (_Complex double*)(((unsigned long int)(_little_dfl_fields_eo)+ALIGN_BASE)&~ALIGN_BASE);
    for (i = 1; i < N_s; i++) {
      little_dfl_fields[i] = little_dfl_fields[i-1] + nb_blocks*9*N_s;
      little_dfl_fields_eo[i] = little_dfl_fields_eo[i-1] + nb_blocks*9*N_s;
    }
    if((void*)(little_A = (_Complex double*)calloc(N_s*N_s, sizeof(_Complex double))) == NULL) {
      return(1);
    }
    if((void*)(little_A_eo = (_Complex double*)calloc(N_s*N_s, sizeof(_Complex double))) == NULL) {
      return(1);
    }   
    init_little_subspace = 1;
  }
  return(0);
}

int init_dfl_subspace(const int N_s) {
  int i;
  init_subspace = 1;
  if((void*)(_dfl_fields = calloc((N_s)*VOLUMEPLUSRAND+1, sizeof(spinor))) == NULL) {
    return(1);
  }
  if ((void*)(dfl_fields = calloc((N_s), sizeof(spinor *))) == NULL) {
    return(1);
  }
  dfl_fields[0] = (spinor*)(((unsigned long int)(_dfl_fields)+ALIGN_BASE)&~ALIGN_BASE);
  for (i = 1; i < N_s; ++i) {
    dfl_fields[i] = dfl_fields[i-1] + VOLUMEPLUSRAND;
  }
  return 0;
}

int free_dfl_subspace() {
  if(init_subspace == 1) {
    free(dfl_fields);
    free(_dfl_fields);
    init_subspace = 0;
  }
  return 0;
}

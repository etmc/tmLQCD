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

 This file was modified according to a flexible number of blocks
 by Claude Tadonki - PetaQCD - April 2010 ( claude.tadonki@lal.in2p3.fr )
 ***********************************************************************/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "global.h"
#include "su3.h"
#include <complex.h>
#include "start.h"
#include "ranlxs.h"
#include "operator/D_psi.h"
#include "poly_precon.h"
#include "Msap.h"
#include "gmres_precon.h"
#include "linalg_eo.h"
#include "gram-schmidt.h"
#include "lu_solve.h"
#include "block.h"
#include "little_D.h"
#include "gcr4complex.h"
#include "boundary.h"
#include <io/params.h>
#include <io/gauge.h>
#include <io/spinor.h>
#include <io/utils.h>
#include "solver_field.h"
#include "generate_dfl_subspace.h"

int init_little_dfl_subspace(const int N_s);

spinor ** dfl_fields = NULL;
static spinor * _dfl_fields = NULL;
_Complex double ** little_dfl_fields = NULL;
static _Complex double *_little_dfl_fields = NULL;
_Complex double ** little_dfl_fields_eo = NULL;
static _Complex double *_little_dfl_fields_eo = NULL;
static int init_subspace = 0;
static int init_little_subspace = 0;

static void random_fields(const int Ns) {
  
  int i, j, ix;
  float r,s[24];
  double *t;
  
  r=(float)(1.0/sqrt(24.0*(double)(VOLUME)));
  
  for (i = 0; i < Ns; i++) {
    t=(double*)(dfl_fields[i]);
    for (ix = 0; ix < VOLUME; ix++){
      ranlxs(s,24);
      for (j = 0; j < 24; j++) {
 	(*t) = (double)(r*(s[j]-0.5f)); 
	(*t) = 1.; 
 	t += 1; 
      }
    }
  }
  return;
}

int generate_dfl_subspace(const int Ns, const int N, const int repro) {
  int ix, i_o,i, j, k, p, blk, vpr = VOLUMEPLUSRAND*sizeof(spinor)/sizeof(_Complex double),
    vol = VOLUME*sizeof(spinor)/sizeof(_Complex double);
  spinor **psi;
  double nrm, e = 0.3, d = 1.1, atime, etime;
  _Complex double s;
  _Complex double * work;
  WRITER *writer = NULL;  
  FILE *fp_dfl_fields; 
  char file_name[500]; // CT
  double musave = g_mu;
  spinor ** work_fields = NULL;
  const int nr_wf = 2;

#ifdef MPI
  atime = MPI_Wtime();
#else
  atime = (double)clock()/(double)(CLOCKS_PER_SEC);
#endif
  init_solver_field(&work_fields, VOLUMEPLUSRAND, nr_wf);
  work = (_Complex double*)malloc(nb_blocks*9*Ns*sizeof(_Complex double));
  psi = (spinor **)calloc(nb_blocks, sizeof(spinor *));
  psi[0] = calloc(VOLUME + nb_blocks, sizeof(spinor));
  for(i = 1; i < nb_blocks; i++) psi[i] = psi[i-1] + (VOLUME / nb_blocks) + 1;
  
  if(init_subspace == 0) i = init_dfl_subspace(Ns);
  
  if(init_little_subspace == 0) i = init_little_dfl_subspace(Ns);
  
  random_fields(Ns);
  if(g_debug_level > 4) {
    for(e = 0.; e < 1.; e=e+0.05) {
      random_spinor_field_lexic(dfl_fields[0], repro, RN_GAUSS);
      nrm = sqrt(square_norm(dfl_fields[0], N, 1));
      mul_r(dfl_fields[0], 1./nrm, dfl_fields[0], N);
      d = 1.1;
      /*       gmres_precon(work_fields[0], dfl_fields[0], 20, 1, 1.e-20, 0, N, &D_psi); */
      poly_nonherm_precon(work_fields[0], dfl_fields[0], e, d, 30, N);
      D_psi(work_fields[1], work_fields[0]);
      diff(work_fields[0], work_fields[1], dfl_fields[0], N);
      nrm = square_norm(work_fields[0], N, 1);
      if(g_proc_id == 0) {
	printf(" e= %f d= %f nrm = %1.5e\n", e, d, nrm);
      }
    }
    d = 1.1;
    e=0.3;
  }
  
  boundary(g_kappa);
  g_mu = 0.;
  /*
    CT: We try to read dfl_fields[i] from file if it exists, 
    otherwise we recalculate it                               
  */
  /* CU: reading and writing should be done with lemon! */
  for(p = 0; p < Ns; p++) {
    sprintf(file_name,"dfl_fields.%.2d", p);
    if((fp_dfl_fields = fopen(file_name,  "r")) == NULL) {
      break;
    }
    else {
      fclose(fp_dfl_fields);
      if((i = read_spinor(dfl_fields[p], NULL, file_name, 0)) != 0) {
	if(g_proc_id == 0) {
	  fprintf(stderr, "Could not read from file %s err = %d\n", file_name, i);
	}
	break;
      }
    }
  }
  
  if((g_proc_id == 0) && (p < Ns) && (g_debug_level > 0))  printf("Compute remaining fields from scratch\n");
  /*CT: We do Ns x 80 x 20 evaluation of Dpsi */
  /*      ModifiedGS((_Complex double*)dfl_fields[i], vol, i, (_Complex double*)dfl_fields[0], vpr); */
  /*      nrm = sqrt(square_norm(dfl_fields[i], N, 1)); */
  /*      mul_r(dfl_fields[i], 1./nrm, dfl_fields[i], N); */
  if(p < Ns) {
    if(1) {
      for(i = 0; i < Ns; i++) {
	/*    ModifiedGS((_Complex double*)dfl_fields[i], vol, i, (_Complex double*)dfl_fields[0], vpr);
	      nrm = sqrt(square_norm(dfl_fields[i], N, 1));
	      mul_r(dfl_fields[i], 1./nrm, dfl_fields[i], N);
	*/
	for(j = 0; j < 20; j++) {
	  zero_spinor_field(g_spinor_field[0],VOLUME);  
	  g_sloppy_precision = 1;
	  Msap_eo(g_spinor_field[0], dfl_fields[i], j+1); 
	  /*      poly_nonherm_precon(g_spinor_field[0], dfl_fields[i], e, d, 2, N);*/
	  /*       gmres_precon(work_fields[0], dfl_fields[i], 20, 1, 1.e-20, 0, N, &D_psi); */
	  
	  for (ix=0;ix<VOLUME;ix++) {
	    _spinor_assign((*(dfl_fields[i] + ix)),(*(g_spinor_field[0]+ix)));
	  }
	  
	  g_sloppy_precision = 0;
	  /*       for (i=0;i<Ns; i++) { */
	  ModifiedGS((_Complex double*)g_spinor_field[0], vol, i, (_Complex double*)dfl_fields[0], vpr);
	  nrm = sqrt(square_norm(g_spinor_field[0], N, 1));
	  mul_r(dfl_fields[i], 1./nrm, g_spinor_field[0], N);
	}
      }
  
      for (i=0; i<Ns; i++) {
	/* test quality */
	if(g_debug_level > -1) {
	  D_psi(work_fields[0], dfl_fields[i]);
	  nrm = sqrt(square_norm(work_fields[0], N, 1));
	  if(g_proc_id == 0) {
	    printf(" ||D psi_%d||/||psi_%d|| = %1.5e\n", i, i, nrm*nrm);
	  }
	}
      }
    }

    if(0) {
      for(j = 0; j < 4; j++) {/*dfl_field_iter = 80  by default */
	for(i = p; i < Ns; i++) {
	  ModifiedGS((_Complex double*)dfl_fields[i], vol, i, (_Complex double*)dfl_fields[0], vpr);
	  nrm = sqrt(square_norm(dfl_fields[i], N, 1));
	  mul_r(dfl_fields[i], 1./nrm, dfl_fields[i], N);
	  for(k = 0; k < 3; k++) {
	    g_sloppy_precision = 1;
	    /* dfl_poly_iter = 20 by default */
	    zero_spinor_field(g_spinor_field[0],VOLUME);
	    Msap_eo(g_spinor_field[0], dfl_fields[i], 4);
	    /* poly_nonherm_precon(g_spinor_field[0], dfl_fields[i], e, d, 4, N);  */
	    g_sloppy_precision = 0;
	    ModifiedGS((_Complex double*)g_spinor_field[0], vol, i, (_Complex double*)dfl_fields[0], vpr);
	    nrm = sqrt(square_norm(g_spinor_field[0], N, 1));
	    mul_r(dfl_fields[i], 1./nrm, g_spinor_field[0], N);
	  }
	  
	  /* test quality */
	  if(g_debug_level > -1) {
	    D_psi(work_fields[0], dfl_fields[i]);
	    nrm = sqrt(square_norm(work_fields[0], N, 1));
	    if(g_proc_id == 0) {
	      printf(" ||D psi_%d||/||psi_%d|| = %1.5e\n", i, i, nrm);
	    }
	  }
	}
      }
    }
    for(i = 0; i < Ns; i++) {
      /*
	CT: We save dfl_fields[i] in a binary file, 
	using a generic nomenclature proc_i__dfl_fields for later reads                               
      */
      sprintf(file_name,"dfl_fields.%.2d", i);
      construct_writer(&writer, file_name, 0);
      write_propagator_type(writer, 4);
      write_spinor(writer, &dfl_fields[i], NULL, 1, 64);
      destruct_writer(writer);
    }
  }
  g_mu = musave;
  g_sloppy_precision = 0;
  boundary(g_kappa);
  if(g_debug_level > 2) {
    for(i = 0; i < Ns; i++) {
      for(j = 0; j < Ns; j++) {
	s = scalar_prod(dfl_fields[i], dfl_fields[j], N, 1);
	if(g_proc_id == 0) {
	  printf("<%d, %d> = %1.3e +i %1.3e\n", i, j, creal(s), cimag(s));
	}
      }
    }
  }
  for (i = 0; i < Ns; i++) {
    /* add it to the basis */
    /* split_global_field(block_list[0].basis[i], block_list[1].basis[i], dfl_fields[i]); */
    split_global_field_GEN_ID(block_list, i, dfl_fields[i], nb_blocks);
  }
  
  /* perform local orthonormalization */
  for(i = 0; i < nb_blocks; i++) block_orthonormalize(block_list+i);
  /* block_orthonormalize(block_list+1); */
  
  dfl_subspace_updated = 1;
  
  for(j = 0; j < Ns; j++) {
    for(i = 0; i < nb_blocks*9*Ns; i++) {
      (little_dfl_fields[j][i]) = 0.0;
      (work[i]) = 0.0;
    }
  }
  
  /* compute the little little basis */
  /* r = work_fields[0]; */
  /* q = g_spinor_field[DUM__SOLVER+1]; */
  
  for(i = 0; i < Ns; i++) {
    /* split_global_field(r, q,  dfl_fields[i]); */
    split_global_field_GEN(psi, dfl_fields[i], nb_blocks);
    /* now take the local scalar products */
    for(j = 0; j < Ns; j++) {
      //p = r;
      for(blk = 0; blk < nb_blocks; blk++) {
	//if(blk == 0) p = r; else p = q;
	little_dfl_fields[i][j + blk*Ns] = scalar_prod(block_list[blk].basis[j], psi[blk], block_list[0].volume, 0);
      }
    }
  }
  
  /* orthonormalise */
  for(i = 0; i < Ns; i++) {
    for (j = 0; j < i; j++) {
      s = lscalar_prod(little_dfl_fields[j], little_dfl_fields[i], nb_blocks*Ns, 1);
      lassign_diff_mul(little_dfl_fields[i], little_dfl_fields[j], s, nb_blocks*Ns);
    }
    s = lsquare_norm(little_dfl_fields[i], nb_blocks*Ns, 1);
    lmul_r(little_dfl_fields[i], 1./sqrt(creal(s)), little_dfl_fields[i], nb_blocks*Ns);
  }
  if(g_debug_level > 0) {
    for(i = 0; i < Ns; i++) {
      for(j = 0; j < Ns; j++) {
	s = lscalar_prod(little_dfl_fields[i], little_dfl_fields[j], nb_blocks*Ns, 1);
	if(g_proc_id == 0) {
	  printf("<%d, %d> = %1.3e +i %1.3e\n", i, j, creal(s), cimag(s));
	}
      }
    }
  }
  
  for(i = 0; i < Ns; i++) {
    little_D(work, little_dfl_fields[i]);
    for(j = 0; j < Ns; j++) {
      little_A[i * Ns + j]  = lscalar_prod(little_dfl_fields[j], work, nb_blocks*Ns, 1);
      if(g_proc_id == 0 && g_debug_level > 4) {
	printf("%1.3e %1.3ei, ", creal(little_A[i * Ns + j]), cimag(little_A[i * Ns + j]));
      }
    }
    if(g_proc_id == 0 && g_debug_level > 4) printf("\n");
  }
  if(g_proc_id == 0 && g_debug_level > 4) printf("\n");
  /* the precision in the inversion is not yet satisfactory! */
  LUInvert(Ns, little_A, Ns);
  /* inverse of little little D now in little_A */


  for(j = 0; j < Ns; j++) {
    for(i = 0; i < nb_blocks*9*Ns; i++) {
      (little_dfl_fields_eo[j][i]) = 0.0;
      (work[i]) = 0.0;
    }
  }

  /* compute the eo little little basis */
  /* r = work_fields[0]; */
  /* q = g_spinor_field[DUM__SOLVER+1]; */
      
  for(i = 0; i < Ns; i++) {
    /* split_global_field(r, q,  dfl_fields[i]); */
    split_global_field_GEN(psi, dfl_fields[i], nb_blocks);
    /* now take the local scalar products */
    for(j = 0; j < Ns; j++) {
      i_o=0;
      for(blk = 0; blk < nb_blocks; blk++) {
         if (block_list[blk].evenodd==1) {
	 little_dfl_fields_eo[i][j + (nb_blocks/2+i_o)*Ns] = scalar_prod(block_list[blk].basis[j], psi[blk], block_list[0].volume, 0);
         i_o++;
	 }	
      }
    }
  }  
     
  /* orthonormalise */
  for(i = 0; i < Ns; i++) {
    for (j = 0; j < i; j++) {
      s = lscalar_prod(little_dfl_fields_eo[j], little_dfl_fields_eo[i], nb_blocks*Ns, 1);
      lassign_diff_mul(little_dfl_fields_eo[i], little_dfl_fields_eo[j], s, nb_blocks*Ns);
    }
    s = lsquare_norm(little_dfl_fields_eo[i], nb_blocks*Ns, 1);
    lmul_r(little_dfl_fields_eo[i], 1./sqrt(creal(s)), little_dfl_fields_eo[i], nb_blocks*Ns);
  }
  if(g_debug_level > 0) {
    for(i = 0; i < Ns; i++) {
      for(j = 0; j < Ns; j++) {
        s = lscalar_prod(little_dfl_fields_eo[i], little_dfl_fields_eo[j], nb_blocks*Ns, 1);
        if(g_proc_id == 0) {
          printf("<%d, %d> = %1.3e +i %1.3e\n", i, j, creal(s), cimag(s));
        }
      }
    }
  }
  
  for(i = 0; i < Ns; i++) {  
    little_D_sym(work, little_dfl_fields_eo[i]);
    for(j = 0; j < Ns; j++) {
      little_A_eo[i * Ns + j]  = lscalar_prod(little_dfl_fields_eo[j], work, nb_blocks*Ns, 1);
      if(g_proc_id == 0 && g_debug_level > 4) {
        printf("%1.3e %1.3ei, ", creal(little_A_eo[i * Ns + j]), cimag(little_A_eo[i * Ns + j])); 
      }
    }
    if(g_proc_id == 0 && g_debug_level > 4) printf("\n");
  }
  if(g_proc_id == 0 && g_debug_level > 4) printf("\n");
  /* the precision in the inversion is not yet satisfactory! */
  LUInvert(Ns, little_A_eo, Ns);
  /* inverse of eo little little D now in little_A_eo */


  
#ifdef MPI
  etime = MPI_Wtime();
#else
  etime = (double)clock()/(double)(CLOCKS_PER_SEC);
#endif
  if(g_proc_id == 0) {
    printf("time for subspace generation %1.3e s\n", etime-atime);
    fflush(stdout);
  }

  finalize_solver(work_fields, nr_wf);
  free_dfl_subspace();
  free(work);
  free(psi[0]);
  free(psi);
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
#if ( defined SSE || defined SSE2 || defined SSE3)
    little_dfl_fields[0] = (_Complex double*)(((unsigned long int)(_little_dfl_fields)+ALIGN_BASE)&~ALIGN_BASE);
    little_dfl_fields_eo[0] = (_Complex double*)(((unsigned long int)(_little_dfl_fields_eo)+ALIGN_BASE)&~ALIGN_BASE);
#else
    little_dfl_fields[0] = _little_dfl_fields;
    little_dfl_fields_eo[0] = _little_dfl_fields_eo;
#endif
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
#if ( defined SSE || defined SSE2 || defined SSE3)
  dfl_fields[0] = (spinor*)(((unsigned long int)(_dfl_fields)+ALIGN_BASE)&~ALIGN_BASE);
#else
  dfl_fields[0] = _dfl_fields;
#endif
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

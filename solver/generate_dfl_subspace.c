/* $Id$ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "global.h"
#include "su3.h"
#include "complex.h"
#include "start.h"
#include "D_psi.h"
#include "poly_precon.h"
#include "gmres_precon.h"
#include "linalg_eo.h"
#include "gram-schmidt.h"
#include "generate_dfl_subspace.h"

spinor ** dfl_fields;
static spinor * _dfl_fields;
static int init_subspace = 0;

int generate_dfl_subspace(const int Ns, const int N) {
  int i,j, vpr = VOLUMEPLUSRAND*sizeof(spinor)/sizeof(complex), 
    vol = VOLUME*sizeof(spinor)/sizeof(complex);
  double nrm, e = 0.5, d = 1., atime, etime;
  complex s;

#ifdef MPI
  atime = MPI_Wtime();
#else
  atime = (double)clock()/(double)(CLOCKS_PER_SEC);
#endif

  if(init_subspace == 0) init_dfl_subspace(Ns);

  for(i = 0; i < Ns; i++) {
    random_spinor_field(dfl_fields[i], N, 0);
/*     random_spinor_field_lexic(dfl_fields[i]); */
    ModifiedGS((complex*)dfl_fields[i], vol, i, (complex*)dfl_fields[0], vpr);
    nrm = sqrt(square_norm(dfl_fields[i], N));
    mul_r(dfl_fields[i], 1./nrm, dfl_fields[i], N);
    for(j = 0; j < 40; j++) {
      poly_nonherm_precon(g_spinor_field[DUM_SOLVER], dfl_fields[i], e, d, 20, N);
/*       gmres_precon(g_spinor_field[DUM_SOLVER], dfl_fields[i], 20, 1, 1.e-20, 0, N, &D_psi); */
      ModifiedGS((complex*)g_spinor_field[DUM_SOLVER], vol, i, (complex*)dfl_fields[0], vpr);
      nrm = sqrt(square_norm(g_spinor_field[DUM_SOLVER], N));
      mul_r(dfl_fields[i], 1./nrm, g_spinor_field[DUM_SOLVER], N);
    }
    /* test quality */
    D_psi(g_spinor_field[DUM_SOLVER], dfl_fields[i]);
    nrm = sqrt(square_norm(g_spinor_field[DUM_SOLVER], N));
    if(g_proc_id == 0 && g_debug_level > -1) {
      printf(" ||D psi_%d||/||psi_%d|| = %1.5e\n", i, i, nrm); 
    }
  }

  if(g_debug_level > 4) {
    for(i = 0; i < Ns; i++) {
      for(j = 0; j < Ns; j++) {
	s = scalar_prod(dfl_fields[i], dfl_fields[j], N);
	if(g_proc_id == 0) {
	  printf("<%d, %d> = %1.3e +i %1.3e\n", i, j, s.re, s.im);
	}
      }
    }
  }
#ifdef MPI
  etime = MPI_Wtime();
#else
  etime = (double)clock()/(double)(CLOCKS_PER_SEC);
#endif
  if(g_proc_id == 0) {
    printf("time for subspace generation %1.3e s\n", etime-atime);
    fflush(stdout);
  }
  return(0);
}

int generate_dfl_subspace_free(const int Ns, const int N) {
  int i,j, vpr = VOLUMEPLUSRAND*sizeof(spinor)/sizeof(complex), 
    vol = VOLUME*sizeof(spinor)/sizeof(complex);
  double nrm, e = 0.5, d = 1.;
  complex s;

  if(init_subspace == 0) init_dfl_subspace(Ns);

  for(i = 0; i < 12; i++) {
    constant_spinor_field(dfl_fields[i], i, N);
    ModifiedGS((complex*)dfl_fields[i], vol, i, (complex*)dfl_fields[0], vpr);
    nrm = sqrt(square_norm(dfl_fields[i], N));
    mul_r(dfl_fields[i], 1./nrm, dfl_fields[i], N);

    /* test quality */
    if(g_debug_level > -1) {
      D_psi(g_spinor_field[DUM_SOLVER], dfl_fields[i]);
      nrm = sqrt(square_norm(g_spinor_field[DUM_SOLVER], N));
      if(g_proc_id == 0) {
	printf(" ||D psi_%d||/||psi_%d|| = %1.5e\n", i, i, nrm); 
      }
    }
  }

  if(g_debug_level > 4) {
    for(i = 0; i < 12; i++) {
      for(j = 0; j < 12; j++) {
	s = scalar_prod(dfl_fields[i], dfl_fields[j], N);
	if(g_proc_id == 0) {
	  printf("<%d, %d> = %1.3e +i %1.3e\n", i, j, s.re, s.im);
	}
      }
    }
  }

  return(0);
}

int init_dfl_subspace(const int N_s) {
  int i;
  init_subspace = 1;
  if((void*)(_dfl_fields = calloc(N_s*VOLUMEPLUSRAND+1, sizeof(spinor))) == NULL) {
    return(1);
  }
  if ((void*)(dfl_fields = calloc(N_s, sizeof(spinor *))) == NULL) {
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

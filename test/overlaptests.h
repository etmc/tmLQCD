#ifndef OVERLAPTESTS_INCLUDE_GUARD
#define OVERLAPTESTS_INCLUDE_GUARD

#include"lime.h"
#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <signal.h>
#ifdef MPI
#include <mpi.h>
#endif
#include "global.h"
#include "getopt.h"
#include "linalg_eo.h"
#include "geometry_eo.h"
#include "start.h"
#include "observables.h"
#ifdef MPI
#include "xchange.h"
#endif
#include "io.h"
#include "io_utils.h"
#include "propagator_io.h"
#include "gauge_io.h"
#include "read_input.h"
#include "mpi_init.h"
#include "sighandler.h"
#include "boundary.h"
#include "solver/solver.h"
#include "init_gauge_field.h"
#include "init_geometry_indices.h"
#include "init_spinor_field.h"
#include "init_moment_field.h"
#include "init_dirac_halfspinor.h"
#include "init_bispinor_field.h"
#include "init_chi_spinor_field.h"
#include "xchange_halffield.h"
#include "stout_smear.h"
#include "invert_eo.h"
#include "monomial.h"
#include "ranlxd.h"
#include "phmc.h"
#include "D_psi.h"
#include "little_D.h"
#include "reweighting_factor.h"
#include "linalg/convert_eo_to_lexic.h"
#include "block.h"
#include "sighandler.h"
#include "Dov_psi.h"

void ov_check_operator(int t, int x, int y, int z);
void ov_check_locality();
void ov_check_ginsparg_wilson_relation(void);
void ov_check_ginsparg_wilson_relation_strong(void);
void ov_compare_4x4(const char * pFileName);
void ov_compare_12x12(const char * pFileName);
void ov_save_12x12(const char * pFileName);

typedef complex matrix4x4[4][4];
typedef complex matrix12x12[12][12];

#define _spinor_norm_l1(d,s) \
   d = 0.; \
   d = _complex_norm((s).s0.c0) + _complex_norm((s).s0.c1) + \
       _complex_norm((s).s0.c2) + _complex_norm((s).s1.c0) + \
       _complex_norm((s).s1.c1) + _complex_norm((s).s1.c2) + \
       _complex_norm((s).s2.c0) + _complex_norm((s).s2.c1) + \
       _complex_norm((s).s2.c2) + _complex_norm((s).s3.c0) + \
       _complex_norm((s).s3.c1) + _complex_norm((s).s3.c2)


#endif

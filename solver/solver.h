/* $Id$ */

#ifndef _SOLVER_H
#define _SOLVER_H

#define BICGSTAB 0
#define CG 1
#define GMRES 2
#define CGS 3
#define MR 4
#define BICGSTABELL 5
#define FGMRES 6
#define GCR 7
#define GMRESDR 8
#define PCG 9
#define DFLGCR 10
#define DFLFGMRES 11

#include"solver/matrix_mult_typedef.h"

#include"solver/gmres.h"
#include"solver/gmres_dr.h"
#include"solver/fgmres.h"
#include"solver/bicgstab_complex.h"
#include"solver/cgs_real.h"
#include"solver/bicgstabell.h"
#include"solver/bicgstab2.h"
#include"solver/cg_her.h"
#include"solver/pcg_her.h"
#include"solver/mr.h"
#include"solver/gcr.h"
#include"solver/eigenvalues.h"

#include"solver/sub_low_ev.h"
#include"solver/gmres_precon.h"
#include"solver/poly_precon.h"

#include "solver/matrix_mult_typedef_bi.h"
#include "solver/bicgstab_complex_bi.h"
#include "solver/cg_her_bi.h"

#include "solver/matrix_mult_typedef_nd.h"
#include "solver/cg_her_nd.h"

#include "solver/generate_dfl_subspace.h"
#endif

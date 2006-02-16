/* $Id$ */

#ifndef _SOLVER_H
#define _SOLVER_H

#define BICGSTAB 0
#define CG 1
#define GMRES 2
#define CGS 3
#define MR 4
#define BICGSTABELL 5

#include"solver/matrix_mult_typedef.h"
#include"solver/gmres.h"
#include"solver/bicgstab_complex.h"
#include"solver/cgs_real.h"
#include"solver/bicgstabell.h"
#include"solver/bicgstab2.h"
#include"solver/cg_her.h"
#include"solver/mr.h"

#include "solver/matrix_mult_typedef_bi.h"
#include "solver/bicgstab_complex_bi.h"
#include "solver/cg_her_bi.h"

#endif

/* $Id$ */

#ifndef _SOLVER_H
#define _SOLVER_H

#define BICGSTAB 0
#define CG 1
#define GMRES 2
#define CGS 3
#define MR 4

#include"solver/matrix_mult_typedef.h"
#include"solver/gmres.h"
#include"solver/bicgstab_complex.h"
#include"solver/cgs_real.h"
#include"solver/bicgstabell.h"
#include"solver/cg_her.h"

#endif

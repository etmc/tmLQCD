#ifndef _INVERT_CLOVER_EO_H
#define _INVERT_CLOVER_EO_H

#include "global.h"
#include "su3.h"
#include "misc_types.h"
#include "solver/matrix_mult_typedef.h"
#include "solver/solver_params.h"

int invert_clover_eo(spinor * const Even_new, spinor * const Odd_new, 
                     spinor * const Even, spinor * const Odd,
                     const double precision, const int max_iter,
                     const int solver_flag, const int rel_prec,
		     const int even_odd_flag, solver_params_t solver_params,
                     su3 *** gf, matrix_mult Qsq, matrix_mult Qm,
                     const ExternalInverter external_inverter, const SloppyPrecision sloppy, const CompressionType compression);

#endif

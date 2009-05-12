#ifndef _CG_MMS_TM_H
#define _CG_MMS_TM_H

#include"matrix_mult_typedef.h"
#include"su3.h"

int cg_mms_tm(spinor * const P,spinor * const Q, const int max_iter, 
	      double eps_sq, const int rel_prec, const int N, matrix_mult f);

#endif

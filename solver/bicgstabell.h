/* $Id$ */

#ifndef _BICGSTABELL_H
#define _BICGSTABELL_H

int bicgstabell(spinor * const x0, spinor * const b, const int max_iter, 
		double eps_sq, const int rel_prec, const int _l, const int N, matrix_mult f);

#endif

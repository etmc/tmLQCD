/* $Id$ */

#ifndef _BICGSTABELL_H
#define _BICGSTABELL_H

int bicgstabell(spinor * const l, spinor * const k, const int max_iter, 
		double eps_sq, const int _l, matrix_mult f);

#endif

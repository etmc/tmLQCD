/* $Id$ */

#ifndef _BICGSTAB2_H
#define _BICGSTAB2_H

int bicgstab2(spinor * const x0, spinor * const b, const int max_iter, 
		double eps_sq, const int rel_prec, const int N, matrix_mult f);

#endif

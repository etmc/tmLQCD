/***********************************************************************
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008 Carsten Urbach
 *
 * This file is part of tmLQCD.
 *
 * tmLQCD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * tmLQCD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with tmLQCD.  If not, see <http://www.gnu.org/licenses/>.
 ***********************************************************************/

void _PSWITCH(ldiff)(_C_TYPE * Q, _C_TYPE * const R, _C_TYPE * const S, const int N);
void _PSWITCH(lassign)(_C_TYPE * Q, _C_TYPE * const R, const int N);
void _PSWITCH(ladd)(_C_TYPE * Q, _C_TYPE * const R, _C_TYPE * const S, const int N);
_F_TYPE _PSWITCH(lsquare_norm)(_C_TYPE * const Q, const int N, const int parallel);
_C_TYPE _PSWITCH(lscalar_prod)(_C_TYPE * const R, _C_TYPE * const S, const int N, const int parallel);
_F_TYPE _PSWITCH(lscalar_prod_r)(_C_TYPE * const R, _C_TYPE * const S, const int N, const int parallel);
void _PSWITCH(lmul_r)(_C_TYPE * const R, const _F_TYPE c, _C_TYPE * const S, const int N);
void _PSWITCH(lmul)(_C_TYPE * const R, const _C_TYPE c, _C_TYPE * const S, const int N);
void _PSWITCH(lassign_diff_mul)(_C_TYPE * const R, _C_TYPE * const S, const _C_TYPE c, const int N);
void _PSWITCH(lassign_add_mul)(_C_TYPE * const R, _C_TYPE * const S, const _C_TYPE c, const int N);
void _PSWITCH(lassign_add_mul_r)(_C_TYPE * const R, _C_TYPE * const S, const _F_TYPE c, const int N);
void _PSWITCH(lassign_mul_add_r)(_C_TYPE * const R, const _F_TYPE c, _C_TYPE * const S, const int N);
void _PSWITCH(ldiff_assign)(_C_TYPE * const Q, _C_TYPE * const S, 
		  const int N);
void _PSWITCH(ladd_assign)(_C_TYPE * const Q, _C_TYPE * const S, 
		  const int N);


int _PSWITCH(gcr4complex)(_C_TYPE * const P, _C_TYPE * const Q, 
			  const int m, const int max_restarts,
			  const double eps_sq, const int rel_prec,
			  const int N, const int parallel,
			  const int lda, const int precon, _PSWITCH(c_matrix_mult) f);



/***********************************************************************
 *
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
 *
 * The externally accessible functions are
 *
 *   int bicg_complex(spinor * const, spinor * const, const int, double, matrix_mult, matrix_mult_dagg)
 *     BiCG solver
 * 
 **************************************************************************/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "global.h"
#include "su3.h"
#include "linalg_eo.h"
#include "start.h"
#include "solver_field.h"
#include "bicg_complex.h"

/* P inout (guess for the solving spinor)
   Q input
 */
int bicg_complex(spinor * const P,spinor * const Q, const int max_iter, 
		double eps_sq, const int rel_prec, 
		const int N, matrix_mult f, matrix_mult fdagg){
	double err, squarenorm;
	_Complex double rho0, rho1, alpha, beta, alphastar, betastar, denom;
	int i;
	//spinor * r, * p, * v, *hatr, * s, * t;
	spinor * p, * phat, * r, * rhat, *tmp, *tmp2;
	spinor ** solver_field = NULL;
	const int nr_sf = 6;

	if(N == VOLUME) {
		init_solver_field(&solver_field, VOLUMEPLUSRAND, nr_sf);
	}
	else {
		init_solver_field(&solver_field, VOLUMEPLUSRAND/2, nr_sf);
	}
	r = solver_field[0];
	rhat = solver_field[1];
	p = solver_field[2];
	phat = solver_field[3];
	tmp = solver_field[4];
	tmp2 = solver_field[5];

	f(tmp, P);
	diff(r, Q, tmp, N); // r = Q - AP
	assign(p, r, N);

	//fdagg(tmp2, P);
	//diff(rhat, Q, tmp2, N); //rhat = Q - Adagg P
	//assign(phat, rhat, N);

	// make rhat different from r, otherwise it won't work
	//random_spinor_field(tmp2, N, 1);
        random_spinor_field_eo(tmp2, 0, RN_GAUSS);
	assign(rhat, tmp2, N);  
	assign(phat, tmp2, N);  

	rho0 = scalar_prod(rhat, r, N, 1);
	squarenorm = square_norm(Q, N, 1);

	printf("rho0 = %f + %fI, squarenorm = %f\n", creal(rho0), cimag(rho0), squarenorm);

	for(i = 0; i < max_iter; i++){
		err = square_norm(r, N, 1);
		if(g_proc_id == g_stdio_proc && g_debug_level > 1) {
			printf("%d %e\n", i, err);
			fflush(stdout);
		}

		if((((err <= eps_sq) && (rel_prec == 0)) || ((err <= eps_sq*squarenorm) && (rel_prec == 1))) && i>0) {
			finalize_solver(solver_field, nr_sf);
			return(i);
		}


		f(tmp, p);
		fdagg(tmp2, phat);
		denom = scalar_prod(phat, tmp, N, 1);
		alpha = rho0/denom;
		alphastar = conj(alpha);

		assign_add_mul(P, p, alpha, N);
		assign_diff_mul(r, tmp, alpha, N);
		assign_diff_mul(rhat, tmp2, alphastar, N);

		rho1 = scalar_prod(rhat, r, N, 1);
		if(fabs(creal(rho1)) < 1.e-25 && fabs(cimag(rho1)) < 1.e-25)
		{
			finalize_solver(solver_field, nr_sf);
			return(-1);
		}
		beta = rho1/rho0;
		betastar = conj(beta);
		mul(tmp, beta, p, N);
		add(p, r, tmp, N);
		mul(tmp2, betastar, phat, N);
		add(phat, rhat, tmp2, N);

		rho0 = rho1;
	}
	finalize_solver(solver_field, nr_sf);
	return -1;
}



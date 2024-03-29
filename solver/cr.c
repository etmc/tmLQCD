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

#ifdef HAVE_CONFIG_H
# include<tmlqcd_config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "global.h"
#include "su3.h"
#include "linalg_eo.h"
#include "solver/gmres_precon.h"
#include "start.h"
#include "operator/tm_operators.h"
#include "solver/poly_precon.h"
#include "solver/cg_her.h"
#include "operator/D_psi.h"
#include "Msap.h"
#include "dfl_projector.h"
#include "solver_field.h"
#include "cr.h"
#include "time.h"
#include "gettime.h"

int cr(spinor * const P, spinor * const Q, 
        const int m, const int max_restarts,
        const double eps_sq, const int rel_prec,
        const int N, const int precon, matrix_mult f) {

    int k, l, restart, i, iter = 0;
    double norm_sq, err;
    spinor * xi, * Axi, * chi, * Achi, *tmp;
    _Complex double alpha, beta;
    static _Complex double one = 1.0;
    double norm, rAr, newrAr;
    double atime, etime;
    spinor ** solver_field = NULL;
    const int nr_sf = 5;
    int save_sloppy = g_sloppy_precision;

    if(N == VOLUME) {
        init_solver_field(&solver_field, VOLUMEPLUSRAND, nr_sf);
    }
    else {
        init_solver_field(&solver_field, VOLUMEPLUSRAND/2, nr_sf);
    }

    atime = gettime();

    xi = solver_field[0];
    Axi = solver_field[1];
    chi = solver_field[2];
    Achi = solver_field[3];
    tmp = solver_field[4];

    norm_sq = square_norm(Q, N, 1);
    if(norm_sq < 1.e-32) {
        norm_sq = 1.;
    }

    dfl_sloppy_prec = 0;
    f(tmp, P);
    diff(chi, Q, tmp, N);
    assign(xi, chi, N);
    f(Axi, xi);
    f(Achi, chi);
    rAr = scalar_prod(chi, Achi, N, 1);
    err = square_norm(chi, N, 1);
    if(((err <= eps_sq) && (rel_prec == 0)) || ((err <= eps_sq*norm_sq) && (rel_prec == 1))) {
        finalize_solver(solver_field, nr_sf);
        return(iter);
    }
    

    for(k = 0; k < m; k++) {

        dfl_sloppy_prec = 1;

        norm = square_norm(Axi, N, 1);
        alpha = rAr/norm;
        assign_add_mul(P, xi, alpha, N);
        /* get the new residual */
        assign_diff_mul(chi, Axi, alpha, N);

        err = square_norm(chi, N, 1);
        iter ++;
        etime = gettime();
        if(g_proc_id == g_stdio_proc && g_debug_level > 3){
            printf("# CR: %d\t%g iterated residue, time spent %f s\n", iter, err, (etime - atime)); 
            fflush(stdout);
        }
        /* Precision reached? */
        if((k == m-1) || ((err <= eps_sq) && (rel_prec == 0)) || ((err <= eps_sq*norm_sq) && (rel_prec == 1))) {
            break;
        }

#ifdef _USE_HALFSPINOR
        if(((err*err <= eps_sq) && (rel_prec == 0)) || ((err*err <= eps_sq*norm_sq) && (rel_prec == 1))) {
            if (g_sloppy_precision_flag == 1) {
                g_sloppy_precision = 1;
                if(g_debug_level > 2 && g_proc_id == g_stdio_proc) {
                    printf("sloppy precision on\n"); fflush( stdout);
                }
            }
        }
#endif

        f(Achi, chi); 

        newrAr = scalar_prod(chi, Achi, N, 1); 
        beta = newrAr/rAr;
        assign_mul_add_mul(xi, beta, chi, one, N);
        assign_mul_add_mul(Axi,beta, Achi, one, N);
        rAr = newrAr;

    }

    g_sloppy_precision = save_sloppy;
    finalize_solver(solver_field, nr_sf);
    return(-1);
}



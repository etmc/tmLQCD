/*****************************************************************************
 * Copyright (C) 2008,2009,2010,2011,2012  
 * Andreas Stathopoulos, Kostas Orginos, Abdou M. Abdel-Rehim
 *
 * This program is based on interfacing the eigCG solver to the tmLQCD code.
 * It was written by Abdou M. Abdel-Rehim based on the original code written
 * by Andreas Stathopoulos and Kostas Orginos and uses functions written in 
 * tmLQCD by Carsten Urbach 
 *
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
 *
 * Gram-Shmidt orthogonalization
 ****************************************************************************/
#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#ifdef MPI
# include <mpi.h>
#endif
#include "global.h"
#include "linalg_eo.h"
#include "start.h"
#include "linalg/blas.h"
#include "linalg/lapack.h"

#include "ortho.h"

/****************************************************************************************************/
int ortho_new_vectors(spinor **Vecs, int N, int nv_old, int nv_new, double orthtol)
{

    //modified Gram-Schmidt orthogonalization 
    int i,j,k;
    int parallel;
    _Complex double alpha;
    int nadded=0;
    double tmpd;
    
    #ifdef MPI
    parallel=1;
    #else
    parallel=0;
    #endif
   
    for(i=nv_old; i< (nv_old+nv_new); i++)
    {
      for(j=0; j<i; j++)
      {
	       alpha=scalar_prod(Vecs[j], Vecs[i],N,parallel);
	       assign_diff_mul(Vecs[i],Vecs[j],alpha,N);
      }
         
      alpha = square_norm(Vecs[i],N,parallel);
      if(creal(alpha) > orthtol*orthtol)
      {
            /* normalize Vecs[i]*/
            tmpd=1.0e+00/sqrt(creal(alpha));
            mul_r(Vecs[i],tmpd,Vecs[i],N);
            nadded= nadded+1;
      }

    }

    return nadded;

} 
/****************************************************************************************************/


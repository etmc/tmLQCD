/***********************************************************************
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
/* 
 *  routine for the initialization of the jocobi field (for use in LapH_ev)
 *  Authors Luigi Scorzato, Marco Cristoforetti
 *
 *
 *******************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include "global.h"
#include "su3.h"
#include "start.h"
#include "xchange/xchange.h"
#include "init_jacobi_field.h"

#ifdef WITHLAPH

su3_vector  *jacobi_field = NULL;

int init_jacobi_field(const int V, const int nr)
{
int i=0;

	if((void*)(jacobi_field = (su3_vector*)calloc(nr*V+1, sizeof(su3_vector))) == NULL)
	{
		printf("malloc errno : %d\n",errno);
		errno = 0;
		return(1);
	}
	if((void*)(g_jacobi_field = (su3_vector**)malloc(nr*sizeof(su3_vector*))) == NULL)
	{
		printf("malloc errno : %d\n",errno);
		errno = 0;
		return(2);
	}

	g_jacobi_field[0] = jacobi_field;
	for(i=1; i<nr; i++)
	{
		g_jacobi_field[i] = g_jacobi_field[i-1]+V;
	}

	return(0);
}


void free_jacobi_field(){
	
	free(jacobi_field);
}


void random_gauss_jacobi_field(su3_vector * const k, const int V)
{
int ix;
su3_vector *s;
double v[6];

 for (ix=0; ix<V ;ix++) {
     s=k+ix;
     gauss_vector(v,6);
     s->c0 = v[0] + v[1] * I;
     s->c1 = v[2] + v[3] * I;
     s->c2 = v[4] + v[5] * I;
 }
#ifdef MPI
 xchange_jacobi(k);
#endif
}

void random_jacobi_field(su3_vector * const k, const int V)
{
int ix;
su3_vector *s;
double v[6];

 for (ix=0; ix<V ;ix++)
   {
     s=k+ix;
     *s=unif_su3_vector();
   }
#ifdef MPI
 xchange_jacobi(k);
#endif
}
#endif // WITHLAPH

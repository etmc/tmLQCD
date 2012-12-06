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
#ifndef _INIT_JACOBI_FIELD_H
#define _INIT_JACOBI_FIELD_H

# ifdef WITHLAPH
int init_jacobi_field(const int V, const int nr);
void free_jacobi_field();
void random_gauss_jacobi_field(su3_vector * const k, const int V);
void random_jacobi_field(su3_vector * const k, const int V);
# endif
#endif

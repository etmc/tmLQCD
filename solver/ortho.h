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


#ifndef _ORTHO_NEW_H
#define _ORTHO_NEW_H

#include "su3.h"

/* Given a set of orthonormal vectors Vecs such that the first nv_old vectors
 * are othonormal, it applies gram-schmit orthogonalization for the new vectors
 * nv_new such that the whole set of nv_old+nv_new are orthonormal. It returns how
 * many vectors were actually added. That could be less than nv_new because of possible
 * linear dependence. If the new orthognalized vector has norm less than orthtol, it is not 
 * added.*/

int ortho_new_vectors(
   spinor **Vecs,  /*the set of vectors*/
   int N,         /* Length of the vectors */
   int nv_old,    /* number of orthonormal vectors */
   int nv_new,    /* number of new vectors*/
   double orthtol /* smallest value of norm of a vector that could be added to Vecs */
);


#endif

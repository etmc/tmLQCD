/***********************************************************************
 * Copyright (C) 2012 Bartosz Kostrzewa
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

/* The two arrays 
 *
 *   g_omp_acc_re
 *   g_omp_acc_cp
 *
 * have as many elements as there are threads (set by ompnumthreads input parameter,
 * stored in omp_num_threads configuration variable). They are initialiazed
 * upon program launch and serve to hold thread-local values over the boundaries
 * of parallel sections, such as for Kahan summations. _re is of type
 * "double" while _cp is of type "_Complex double". They are declared in global.h */

#ifndef _INIT_OMP_ACCUMULATORS_H
#define _INIT_OMP_ACCUMULATORS_H

int init_omp_accumulators(int num);
void free_omp_accumulators();

#endif

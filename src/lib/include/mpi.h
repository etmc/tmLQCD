/***********************************************************************
 *
 * Copyright (C) 2026 Roman Gruber
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
 * Simple MPI header wrapper
 *
 * Author: Roman Gruber
 *         roman.gruber@unibe.ch
 *
 *******************************************************************************/

#ifndef MY_MPI_WRAPPER_H
#define MY_MPI_WRAPPER_H


// include *real* MPI header
#include_next <mpi.h>


/**
 * @brief      MPI context
 *
 * @var        comm MPI communicator
 */
typedef struct {
    MPI_Comm comm;
} MPIContext;


/**
 * @brief      The global application context struct
 *
 * @var        mpi MPI context
 */
typedef struct {
    MPIContext mpi;
} AppContext;


/**
 * @brief      Return the global application context struct
 *
 * @return     Global application context struct
 */
const AppContext* app(void);


/**
 * @brief      Initialize application context
 *
 * @param[in]  comm  The MPI communicator to use throughout the application
 */
void app_context_init(const MPI_Comm comm);


/**
 * @brief      Finalize application context
 */
void app_context_finalize(void);


#endif
